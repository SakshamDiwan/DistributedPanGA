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
// final bepos is the same on any rank. Ranks publish them with reductions:
//
//   E[x]  latest partner before x (in order) whose pair (x,b) made a call
//   L[x]  latest partner of x, other than x, whose pair made a call
//   S[x]  state of self pair (x,x): 0 none, 1 pending, 2 done without a
//         call, 3+bepos done with a call
//
// all indexed by position x in FastGA's order, with keys (pos(b)+1)<<32 | bepos
// so that MPI_MAX picks the latest partner.
//
// Where FastGA has no earlier call in the range it reads a field it never
// set. That cannot be derived; PGA_FASTGA_UNSET_BEPOS (default 32767, the
// value observed in FastGA 178afcb on Perlmutter) is used.

#include "fastga_order.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static MPI_Comm Comm;
static int      My_Rank;
static int      N;               // contigs
static int     *Pos_Of;          // contig -> position in FastGA's order
static int     *Contig_At;       // position -> contig (the .gix permutation)
static int     *Range_Start;     // position -> first position of its thread range
static int64_t *E_loc, *E_send, *E_glob;
static int64_t *L_loc, *L_glob;
static int64_t *S_loc, *S_glob;
static MPI_Request E_req = MPI_REQUEST_NULL;
static int      Unset_Bepos = 32767;   // what FastGA 178afcb read there on Perlmutter (yeast7, 8 and 32 threads)
static int      Log_Values  = 0;
static int      Round       = 0;
static int64_t  Last_Remaining = -1;

static struct {
    int64_t simple, fallback, unset;
    double  wait_pass1_s, rounds_s;
    int     rounds;
} St;

static SelfOrderHooks Hooks;

static int64_t make_key(int pos_b, int bepos)
{ return (((int64_t) (pos_b+1)) << 32) | (uint32_t) bepos; }

static int key_pos(int64_t key)   { return (int) (key >> 32) - 1; }
static int key_bepos(int64_t key) { return (int) (key & 0xffffffffll); }

// --- hooks -------------------------------------------------------------------

static int h_pair_class(int a, int b)
{ if (a == b) return 0;
  return (Pos_Of[b] < Pos_Of[a]) ? -1 : 1;
}

static void h_record(int a, int b, int calls, int bepos)
{ int xa = Pos_Of[a];

  if (a == b)
    { S_loc[xa] = (calls > 0) ? 3 + (int64_t) bepos : 2;
      return;
    }
  if (calls > 0)
    { int     xb  = Pos_Of[b];
      int64_t key = make_key(xb,bepos);
      if (key > L_loc[xa]) L_loc[xa] = key;
      if (xb < xa && key > E_loc[xa]) E_loc[xa] = key;
    }
}

static void h_pass1_done(void)
{ memcpy(E_send,E_loc,N*sizeof(int64_t));
  MPI_Iallreduce(E_send,E_glob,N,MPI_INT64_T,MPI_MAX,Comm,&E_req);
}

static void h_pass2_done(void)
{ double t0 = MPI_Wtime();
  MPI_Wait(&E_req,MPI_STATUS_IGNORE);
  St.wait_pass1_s += MPI_Wtime()-t0;
}

static int h_self_value(int c, int *bepos)
{ int x = Pos_Of[c];

  S_loc[x] = 1;
  if (E_glob[x] != 0)
    { *bepos = key_bepos(E_glob[x]);
      St.simple += 1;
      if (Log_Values)
        fprintf(stderr,"[fga-self] contig=%d pos=%d bepos=%d from=partner:%d\n",
                c,x,*bepos,Contig_At[key_pos(E_glob[x])]);
      return 1;
    }
  return 0;
}

static void h_round_begin(void)
{ double t0 = MPI_Wtime();
  if (Round == 0)
    MPI_Allreduce(L_loc,L_glob,N,MPI_INT64_T,MPI_MAX,Comm);
  MPI_Allreduce(S_loc,S_glob,N,MPI_INT64_T,MPI_MAX,Comm);
  St.rounds_s += MPI_Wtime()-t0;
  Round += 1;
}

// Walk back from (c,c) through earlier contigs of the same thread range to
// the nearest row with a real call. In row y the order is: partners before
// y, then (y,y), then partners after y; so a partner after y is the row's
// last call if there is one, else the self pair's, else a partner before y.
static int h_resolve(int c, int *bepos)
{ int x = Pos_Of[c];
  int y;

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

  *bepos = Unset_Bepos;
  St.unset += 1;
  fprintf(stderr,"[fga-self] contig=%d pos=%d bepos=%d from=unset (first call of FastGA "
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
    { fprintf(stderr,"rank %d: fastga_order: %lld deferred self pairs cannot be resolved\n",
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

int fastga_order_setup(const char *gdb_stem, GDB *gdb, const int64_t *n_seed_count,
                       int swide, MPI_Comm comm)
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
  if (env != NULL) Unset_Bepos = atoi(env);
  env = getenv("PGA_FASTGA_ORDER_LOG");
  Log_Values = (env != NULL && strcmp(env,"0") != 0);

  gcount = malloc(N*sizeof(int64_t));
  MPI_Allreduce(n_seed_count,gcount,N,MPI_INT64_T,MPI_SUM,comm);
  nparts = compute_ranges(gdb,gcount,swide,T,&nrange);
  free(gcount);

  E_loc  = calloc(N,sizeof(int64_t));
  E_send = calloc(N,sizeof(int64_t));
  E_glob = calloc(N,sizeof(int64_t));
  L_loc  = calloc(N,sizeof(int64_t));
  L_glob = calloc(N,sizeof(int64_t));
  S_loc  = calloc(N,sizeof(int64_t));
  S_glob = calloc(N,sizeof(int64_t));
  memset(&St,0,sizeof(St));

  Hooks.pair_class  = h_pair_class;
  Hooks.record      = h_record;
  Hooks.pass1_done  = h_pass1_done;
  Hooks.pass2_done  = h_pass2_done;
  Hooks.self_value  = h_self_value;
  Hooks.round_begin = h_round_begin;
  Hooks.resolve     = h_resolve;
  Hooks.round_end   = h_round_end;
  Self_Order = &Hooks;

  if (My_Rank == 0)
    { fprintf(stderr,"[fga-order] on: FastGA order from %s.gix, emulating %d threads: "
                     "%d parts, %d thread ranges; unset value %d\n",
              gdb_stem,T,nparts,nrange,Unset_Bepos);
    }
  return 1;
}

void fastga_order_report(void)
{ if (Self_Order == NULL) return;
  fprintf(stderr,"[fga-order] rank %d: self pairs from_partner=%lld from_earlier_row=%lld "
                 "unset=%lld wait_pass1_s=%.6f rounds=%d rounds_s=%.6f\n",
          My_Rank,(long long) St.simple,(long long) St.fallback,(long long) St.unset,
          St.wait_pass1_s,St.rounds,St.rounds_s);
  Self_Order = NULL;
}
