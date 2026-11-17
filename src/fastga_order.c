// fastga_order.c — give each same-strand self pair the value FastGA's
// alignment loop would have had on entry to it.
//
// Background. For a self hit that straddles the main diagonal, FastGA's loop
// (align_contigs) skips the aligner and then reads path->bepos, still set by
// the previous real Local_Alignment call. FastGA aligns pairs in a fixed
// order, so that value is fixed there. In DistributedPanGA it would depend on
// which pairs a rank owns. This file reconstructs FastGA's value:
//
//   order     contigs by the permutation stored in <stem>.gix (longest first);
//             anchor A in that order, then partner B in that order
//   ranges    FastGA cuts the contigs into parts (by length) and each part
//             into thread ranges (by seed count); the value does not carry
//             across a range start (FastGA.c main, RSDsort.c rmsd_sort)
//   v(c)      bepos of the last real aligner call before pair (c,c) in its
//             range: the nearest earlier partner of c that made a call, else
//             the last call of the nearest earlier contig's row
//
// Pairs of two different contigs never read the leftover value, so their
// final bepos is the same on any rank. Each rank therefore sends the result
// of every pair (x,b), b before x in FastGA's order, to the rank that owns
// the self pair (x,x): one small message per pair. The owner knows, from a
// map of which pairs have seeds, the earlier partners of x in order of
// nearness, and can start (x,x) once it has heard from the nearest partner
// that made a call and from every partner nearer than that. Ranks align
// their earlier-partner pairs nearest first, so those answers come early; a
// self pair whose value is not in yet is set aside and aligned after the
// rank's other pairs. Messages are received in poll(), which the alignment
// loop calls often.
//
// The rare value that comes from an earlier contig's row is resolved at the
// end with two reductions, indexed by position x in FastGA's order:
//
//   L[x]  (pos(b)+1)<<32 | bepos for the latest partner of x, other than x,
//         whose pair made a call (MPI_MAX picks the latest)
//   S[x]  self pair (x,x): 0 none, 1 pending, 2 done without a call,
//         3+bepos done with a call
//
// Where FastGA has no earlier call in the range it reads a field it never
// set. That cannot be derived. A diagnostic build of FastGA 178afcb on
// Perlmutter showed 32764..32767 there for the first range of the run (main
// thread; the high half of a stack address) and 0 for other ranges, on four
// datasets. We use PGA_FASTGA_UNSET_BEPOS (default 32767) and 0.

#include "fastga_order.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MSG_TAG 7701

static MPI_Comm Comm;
static int      My_Rank;
static int      N;               // contigs
static int     *Pos_Of;          // contig -> position in FastGA's order
static int     *Contig_At;       // position -> contig (the .gix permutation)
static int     *Range_Start;     // position -> first position of its thread range
static uint8_t *Has;             // Has[xa*N+xb], xb < xa: the pair has same-strand seeds
static const int8_t *Owner;      // rank_for_pair[a*N+b], by contig
static int64_t *Clen;            // contig length by position
static int    **Answer;          // rows whose self pair this rank owns: per earlier
static int     *NNear;           //   partner, nearest first: -1 unknown, -2 no call, else bepos
static int64_t  Expect, Got;     // messages this rank must receive in all
static int64_t *L_loc, *L_glob;
static int64_t *S_loc, *S_glob;
static int      Unset_Main = 32767;
static int      Log_Values = 0;
static int      Round      = 0;
static int64_t  Last_Remaining = -1;

static struct {
    int64_t direct, late, fallback, unset, sent;
    double  wait_s, rounds_s;
    int     rounds;
} St;

static SelfOrderHooks Hooks;

static int64_t make_key(int pos_b, int bepos)
{ return (((int64_t) (pos_b+1)) << 32) | (uint32_t) bepos; }

static int key_pos(int64_t key)   { return (int) (key >> 32) - 1; }
static int key_bepos(int64_t key) { return (int) (key & 0xffffffffll); }

// Nearness of earlier partner xb to xa: how many earlier partners with seeds
// lie between them. 0 = the nearest.
static int nearness(int xa, int xb)
{ const uint8_t *row = Has + (int64_t) xa*N;
  int y, k = 0;

  for (y = xb+1; y < xa; y++)
    k += row[y];
  return k;
}

static void store_answer(int xa, int k, int called, int bepos)
{ Got += 1;
  if (Answer[xa] != NULL)
    Answer[xa][k] = called ? bepos : -2;
}

static void receive_one(MPI_Status *st)
{ int msg[4];

  MPI_Recv(msg,4,MPI_INT,st->MPI_SOURCE,MSG_TAG,Comm,MPI_STATUS_IGNORE);
  store_answer(msg[0],msg[1],msg[2],msg[3]);
}

// The entry value of self pair c from the answers so far.
//   1  known (*bepos set)      0  an answer it needs has not arrived
//  -1  no earlier partner made a call and c is not first in its FastGA range:
//      the value comes from an earlier contig's row (final rounds)
static int from_answers(int c, int *bepos, const char *how)
{ int x = Pos_Of[c];
  int k;

  for (k = 0; k < NNear[x]; k++)
    { int a = Answer[x][k];
      if (a == -1)
        return 0;
      if (a >= 0)
        { *bepos = a;
          if (Log_Values)
            { int y, j = k;                 // the k-th nearest earlier partner
              for (y = x-1; y >= 0; y--)
                if (Has[(int64_t) x*N+y] && j-- == 0)
                  break;
              fprintf(stderr,"[fga-self] contig=%d pos=%d bepos=%d from=partner:%d nearness=%d %s\n",
                      c,x,*bepos,Contig_At[y],k,how);
            }
          return 1;
        }
    }
  if (Range_Start[x] != x)
    return -1;
  *bepos = (x == 0) ? Unset_Main : 0;
  St.unset += 1;
  fprintf(stderr,"[fga-self] contig=%d pos=%d bepos=%d from=unset (first call of the FastGA "
                 "range starting at pos %d) %s\n",c,x,*bepos,x,how);
  return 1;
}

// --- hooks -------------------------------------------------------------------

static int h_pair_class(int a, int b)
{ if (a == b) return 0;
  return (Pos_Of[b] < Pos_Of[a]) ? -1 : 1;
}

// Earlier-partner pairs: nearest first, and among equals the longer anchor
// first, since its self pair is the more expensive one to keep waiting.
static int64_t h_order_key(int a, int b)
{ int xa = Pos_Of[a];
  return (((int64_t) nearness(xa,Pos_Of[b])) << 36) | ((((int64_t) 1) << 36) - 1 - Clen[xa]);
}

static void drain(void)
{ MPI_Status st;
  int flag;

  while (1)
    { MPI_Iprobe(MPI_ANY_SOURCE,MSG_TAG,Comm,&flag,&st);
      if (!flag) break;
      receive_one(&st);
    }
}

// Called from the alignment loop for every aligner call: look for messages
// on every 16th call only, to keep the cost out of the loop.
static void h_poll(void)
{ static unsigned calls = 0;

  if ((++calls & 15) == 0)
    drain();
}

static void h_record(int a, int b, int calls, int bepos)
{ int xa = Pos_Of[a];

  if (a == b)
    { S_loc[xa] = (calls > 0) ? 3 + (int64_t) bepos : 2;
      return;
    }
  if (calls > 0)
    { int64_t key = make_key(Pos_Of[b],bepos);
      if (key > L_loc[xa]) L_loc[xa] = key;
    }
  if (Pos_Of[b] < xa)
    { int owner = Owner[(int64_t) a*N + a];
      int k     = nearness(xa,Pos_Of[b]);

      if (owner == My_Rank)
        store_answer(xa,k,calls > 0,bepos);
      else
        { int msg[4] = { xa, k, calls > 0, bepos };
          MPI_Send(msg,4,MPI_INT,owner,MSG_TAG,Comm);
          St.sent += 1;
        }
    }
}

static int h_self_value(int c, int *bepos)
{ S_loc[Pos_Of[c]] = 1;
  drain();
  if (from_answers(c,bepos,"(in order)") == 1)
    { St.direct += 1;
      return 1;
    }
  return 0;
}

static int h_self_ready(int c, int *bepos)
{ int r;

  drain();
  r = from_answers(c,bepos,"(set aside)");
  if (r == 1) St.late += 1;
  return r;
}

static void h_wait_message(void)
{ MPI_Status st;
  double t0 = MPI_Wtime();

  MPI_Probe(MPI_ANY_SOURCE,MSG_TAG,Comm,&st);
  receive_one(&st);
  St.wait_s += MPI_Wtime()-t0;
}

static void h_round_begin(void)
{ double t0 = MPI_Wtime();

  if (Round == 0)
    { MPI_Status st;                      // every rank has sent all its answers
      while (Got < Expect)
        { MPI_Probe(MPI_ANY_SOURCE,MSG_TAG,Comm,&st);
          receive_one(&st);
        }
      MPI_Allreduce(L_loc,L_glob,N,MPI_INT64_T,MPI_MAX,Comm);
    }
  MPI_Allreduce(S_loc,S_glob,N,MPI_INT64_T,MPI_MAX,Comm);
  St.rounds_s += MPI_Wtime()-t0;
  Round += 1;
}

// A self pair with no earlier partner that made a call, not first in its
// range: walk back through earlier contigs of the range to the nearest row
// with a real call. In row y the order is: partners before y, then (y,y),
// then partners after y; so a partner after y is the row's last call if
// there is one, else the self pair's, else a partner before y.
static int h_resolve(int c, int *bepos)
{ int x = Pos_Of[c];
  int y, r;

  r = from_answers(c,bepos,"(final)");
  if (r == 1) { St.late += 1; return 1; }
  if (r == 0) return 0;

  for (y = x-1; y >= Range_Start[x]; y--)
    { int64_t k = L_glob[y];
      int64_t s = S_glob[y];
      int     from = -1;

      if (k != 0 && key_pos(k) > y)
        { *bepos = key_bepos(k); from = Contig_At[key_pos(k)]; }
      else if (s == 1)
        return 0;                       // that self pair has not run yet
      else if (s >= 3)
        { *bepos = (int) (s-3); from = Contig_At[y]; }
      else if (k != 0)
        { *bepos = key_bepos(k); from = Contig_At[key_pos(k)]; }
      else
        continue;

      St.fallback += 1;
      fprintf(stderr,"[fga-self] contig=%d pos=%d bepos=%d from=row:%d,%d\n",
              c,x,*bepos,Contig_At[y],from);
      return 1;
    }

  *bepos = (Range_Start[x] == 0) ? Unset_Main : 0;
  St.unset += 1;
  fprintf(stderr,"[fga-self] contig=%d pos=%d bepos=%d from=unset (first call of the FastGA "
                 "range starting at pos %d)\n",c,x,*bepos,Range_Start[x]);
  return 1;
}

static int h_round_end(int64_t remaining)
{ int64_t all = 0;
  double  t0 = MPI_Wtime();

  MPI_Allreduce(&remaining,&all,1,MPI_INT64_T,MPI_SUM,Comm);
  St.rounds_s += MPI_Wtime()-t0;
  St.rounds = Round;
  if (all > 0 && all == Last_Remaining)
    { fprintf(stderr,"rank %d: fastga_order: %lld self pairs cannot be resolved\n",
              My_Rank,(long long) all);
      MPI_Abort(Comm,5);
    }
  Last_Remaining = all;
  return (all > 0);
}

// --- setup -------------------------------------------------------------------

// The contig permutation in the .gix stub (layout: FastGA.c Open_Post_List).
static int read_gix_perm(const char *stem, int n, int *perm)
{ char    path[4096];
  int     f, pbyte, cbyte, nfile, freq, nctg, i, ok;
  int64_t maxp;
  char   *seen;

  snprintf(path,sizeof(path),"%s.gix",stem);
  f = open(path,O_RDONLY);
  if (f < 0)
    { fprintf(stderr,"[fga-order] cannot open %s\n",path); return 0; }
  ok = (lseek(f,4*sizeof(int)+0x1000000*sizeof(int64_t),SEEK_SET) >= 0);
  ok = ok && read(f,&pbyte,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&cbyte,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&nfile,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&maxp,sizeof(int64_t)) == sizeof(int64_t);
  ok = ok && read(f,&freq,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&nctg,sizeof(int)) == sizeof(int);
  if (!ok || nctg != n)
    { fprintf(stderr,"[fga-order] %s: unexpected header (contigs %d, GDB has %d)\n",
              path,ok ? nctg : -1,n);
      close(f);
      return 0;
    }
  ok = (read(f,perm,sizeof(int)*(size_t) n) == (ssize_t) (sizeof(int)*(size_t) n));
  close(f);
  if (!ok)
    { fprintf(stderr,"[fga-order] %s: short read of the contig order\n",path); return 0; }
  seen = calloc(n,1);
  for (i = 0; i < n; i++)
    { if (perm[i] < 0 || perm[i] >= n || seen[perm[i]])
        { fprintf(stderr,"[fga-order] %s: contig order is not a permutation\n",path);
          free(seen);
          return 0;
        }
      seen[perm[i]] = 1;
    }
  free(seen);
  return 1;
}

// FastGA's parts (FastGA.c, "Compute GDB split into NTHREADS parts") and, in
// each part, rmsd_sort's thread ranges over the per-anchor seed bytes.
static int compute_ranges(GDB *gdb, const int64_t *count, int swide, int T, int *nrange_out)
{ int     *split = malloc((T+1)*sizeof(int));
  int64_t  npost, cum, t;
  int      p, r, x, i, nparts, nrange;

  npost = gdb->seqtot;
  split[0] = 0;
  p = 0;
  r = T;
  t = npost/T;
  cum = gdb->contigs[Contig_At[0]].clen;
  for (x = 1; x < N; x++)
    { if (cum >= t && x >= r)
        { p += 1;
          split[p] = x;
          t = (npost*(p+1))/T;
          r += T;
        }
      cum += gdb->contigs[Contig_At[x]].clen;
    }
  nparts = p+1;
  split[nparts] = N;

  for (x = 0; x < N; x++)
    Range_Start[x] = x;

  nrange = 0;
  for (i = 0; i < nparts; i++)
    { int64_t asize, thr, sum;
      int     n, beg, y;

      asize = 0;
      for (x = split[i]; x < split[i+1]; x++)
        asize += count[Contig_At[x]]*swide;
      if (asize == 0) continue;

      n   = 0;
      thr = asize/T;
      sum = 0;
      for (x = split[i]; x < split[i+1]; x++)
        if (count[Contig_At[x]] > 0)
          break;
      beg = x;
      for ( ; x < split[i+1]; x++)
        if (count[Contig_At[x]] > 0)
          { sum += count[Contig_At[x]]*swide;
            if (sum >= thr)
              { for (y = beg; y <= x; y++)
                  Range_Start[y] = beg;
                if (Log_Values && My_Rank == 0)
                  fprintf(stderr,"[fga-range] beg=%d end=%d\n",beg,x+1);
                n  += 1;
                thr = (asize*(n+1))/T;
                beg = x+1;
              }
          }
      nrange += n;
    }
  free(split);
  *nrange_out = nrange;
  return nparts;
}

static int64_t read_field(const uint8_t *p, int nbytes)
{ int64_t v = 0;
  memcpy(&v,p,nbytes);
  return v;
}

int fastga_order_setup(const char *gdb_stem, GDB *gdb, const int64_t *n_seed_count,
                       int swide, const FastgaOrderSeeds *seeds,
                       const int8_t *rank_for_pair, MPI_Comm comm)
{ int      ok, T, nparts, nrange, i;
  int64_t *gcount;
  const char *env;

  Comm = comm;
  MPI_Comm_rank(comm,&My_Rank);
  N = gdb->ncontig;

  Contig_At = malloc(N*sizeof(int));
  ok = 0;
  if (My_Rank == 0)
    ok = read_gix_perm(gdb_stem,N,Contig_At);
  MPI_Bcast(&ok,1,MPI_INT,0,comm);
  if (!ok)
    { if (My_Rank == 0)
        fprintf(stderr,"[fga-order] disabled: pairs are aligned in sorted order as before\n");
      free(Contig_At);
      return 0;
    }
  MPI_Bcast(Contig_At,N,MPI_INT,0,comm);

  Pos_Of      = malloc(N*sizeof(int));
  Range_Start = malloc(N*sizeof(int));
  for (i = 0; i < N; i++)
    Pos_Of[Contig_At[i]] = i;

  T = 32;
  env = getenv("PGA_FASTGA_THREADS");
  if (env != NULL && atoi(env) > 0) T = atoi(env);
  env = getenv("PGA_FASTGA_UNSET_BEPOS");
  if (env != NULL) Unset_Main = atoi(env);
  env = getenv("PGA_FASTGA_ORDER_LOG");
  Log_Values = (env != NULL && strcmp(env,"0") != 0);

  gcount = malloc(N*sizeof(int64_t));
  MPI_Allreduce(n_seed_count,gcount,N,MPI_INT64_T,MPI_SUM,comm);
  nparts = compute_ranges(gdb,gcount,swide,T,&nrange);
  free(gcount);

  // Which pairs (x,b), b before x, have same-strand seeds: one scan of this
  // rank's seeds, then a reduction. Every rank then knows each anchor's
  // earlier partners in order of nearness.
  if ((int64_t) N*N > 0x7fffffffll)
    { if (My_Rank == 0)
        fprintf(stderr,"[fga-order] disabled: %d contigs is too many for the pair map\n",N);
      return 0;
    }
  { uint8_t *has_loc = calloc((size_t) N*N,1);
    int      unit  = 1 + seeds->apos_bytes + seeds->actg_bytes
                       + seeds->bpos_bytes + seeds->bctg_bytes;
    int      aoff  = 1 + seeds->apos_bytes;
    int      boff  = aoff + seeds->actg_bytes + seeds->bpos_bytes;
    int64_t  amask = (((int64_t) 1) << (8*seeds->actg_bytes-1)) - 1;
    int64_t  bmask = (((int64_t) 1) << (8*seeds->bctg_bytes-1)) - 1;
    const uint8_t *p   = seeds->data;
    const uint8_t *end = p + seeds->len;

    Has = calloc((size_t) N*N,1);
    for ( ; p + unit <= end; p += unit)
      { int a  = (int) (read_field(p+aoff,seeds->actg_bytes) & amask);
        int b  = (int) (read_field(p+boff,seeds->bctg_bytes) & bmask);
        int xa = Pos_Of[a];
        int xb = Pos_Of[b];
        if (xb < xa)
          has_loc[(int64_t) xa*N+xb] = 1;
      }
    MPI_Allreduce(has_loc,Has,N*N,MPI_BYTE,MPI_BOR,comm);
    free(has_loc);
  }

  Owner  = rank_for_pair;
  Clen   = malloc(N*sizeof(int64_t));
  Answer = calloc(N,sizeof(int *));
  NNear  = calloc(N,sizeof(int));
  Expect = Got = 0;
  for (i = 0; i < N; i++)
    { int c = Contig_At[i], y, k;

      Clen[i] = gdb->contigs[c].clen;
      for (y = 0; y < i; y++)
        NNear[i] += Has[(int64_t) i*N+y];
      if (rank_for_pair[(int64_t) c*N + c] != My_Rank)
        continue;
      Expect += NNear[i];
      Answer[i] = malloc((NNear[i]+1)*sizeof(int));
      for (k = 0; k < NNear[i]; k++)
        Answer[i][k] = -1;
    }

  L_loc   = calloc(N,sizeof(int64_t));
  L_glob  = calloc(N,sizeof(int64_t));
  S_loc   = calloc(N,sizeof(int64_t));
  S_glob  = calloc(N,sizeof(int64_t));
  memset(&St,0,sizeof(St));

  Hooks.pair_class   = h_pair_class;
  Hooks.order_key    = h_order_key;
  Hooks.record       = h_record;
  Hooks.poll         = h_poll;
  Hooks.self_value   = h_self_value;
  Hooks.self_ready   = h_self_ready;
  Hooks.wait_message = h_wait_message;
  Hooks.round_begin  = h_round_begin;
  Hooks.resolve      = h_resolve;
  Hooks.round_end    = h_round_end;
  Self_Order = &Hooks;

  if (My_Rank == 0)
    fprintf(stderr,"[fga-order] on: FastGA order from %s.gix, emulating %d threads: "
                   "%d parts, %d thread ranges; unset value %d (first range) or 0\n",
            gdb_stem,T,nparts,nrange,Unset_Main);
  return 1;
}

void fastga_order_report(void)
{ if (Self_Order == NULL) return;
  fprintf(stderr,"[fga-order] rank %d: self pairs in_order=%lld set_aside=%lld "
                 "from_earlier_row=%lld unset=%lld answers_sent=%lld answers_received=%lld "
                 "wait_s=%.6f rounds=%d rounds_s=%.6f\n",
          My_Rank,(long long) St.direct,(long long) St.late,
          (long long) St.fallback,(long long) St.unset,(long long) St.sent,(long long) Got,
          St.wait_s,St.rounds,St.rounds_s);
  Self_Order = NULL;
}
