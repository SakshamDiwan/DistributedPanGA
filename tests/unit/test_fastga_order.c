// ============================================================================
// F-01..F-06 -- Reproducing FastGA's order of self pairs (src/fastga_order_core.c).
// See tests/CATALOG.md.
// ============================================================================
//
// WHAT IS BEING PROTECTED.  FastGA's alignment loop, after skipping a
// same-strand self hit on the main diagonal, continues from path->bepos,
// which the previous real aligner call left.  FastGA aligns pairs in a fixed
// order, so that value is fixed there.  DistributedPanGA gives each self pair
// the value FastGA would have had: it recomputes FastGA's contig order and
// thread ranges, collects the final bepos of the right earlier pair from
// whichever rank aligned it, and falls back to earlier rows of the range in
// rare cases.  Get any rule wrong and the output silently differs from
// FastGA on some datasets or rank counts.
//
// HOW.  The rules are MPI-free (fastga_order_core.c), so this binary checks
// them in one process:
//
//   F-01  the .gix reader: header layout, wrong contig count, a broken
//         permutation, a missing file
//   F-02  parts and thread ranges against a verbatim copy of FastGA's two
//         splitting loops (FastGA.c main, RSDsort.c rmsd_sort), on random
//         inputs including contigs without seeds
//   F-03  nearness helpers
//   F-04  answers: waiting, known, unset, and "needs earlier rows"
//   F-05  the whole rule set against a sequential model of FastGA: random
//         worlds of contigs, pairs and outcomes, walked in FastGA's order,
//         compared with what the distributed rules derive -- including
//         self pairs that depend on other self pairs over several rounds,
//         which no dataset run has exercised
//   F-06  the rounds always finish, and a chain of k dependent self pairs
//         takes k rounds

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>

#include "check.h"
#include "fastga_order_core.h"

static uint64_t rng = 0x9e3779b97f4a7c15ull;
static uint64_t rnd(void)
{ rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static int rndi(int n) { return (int) (rnd() % (uint64_t) n); }

static const char *tmpdir(void)
{ const char *t = getenv("PGA_TEST_TMP");
  return (t && *t) ? t : "tests/bin";
}

// ---------------------------------------------------------------------------
// F-01  .gix reader
// ---------------------------------------------------------------------------

static void write_gix(const char *stem, int nctg, const int *perm)
{ char path[512];
  int  f, zero = 0, nfile = 1, freq = 10;
  int64_t maxp = 100;

  snprintf(path,sizeof(path),"%s.gix",stem);
  f = open(path,O_WRONLY|O_CREAT|O_TRUNC,0644);
  CHECK(f >= 0);
  // 4 ints and the 2^24-entry index: leave them as a hole (sparse file)
  CHECK(lseek(f,4*sizeof(int)+0x1000000*sizeof(int64_t),SEEK_SET) >= 0);
  CHECK(write(f,&zero,sizeof(int)) == sizeof(int));     // pbyte
  CHECK(write(f,&zero,sizeof(int)) == sizeof(int));     // cbyte
  CHECK(write(f,&nfile,sizeof(int)) == sizeof(int));
  CHECK(write(f,&maxp,sizeof(int64_t)) == sizeof(int64_t));
  CHECK(write(f,&freq,sizeof(int)) == sizeof(int));
  CHECK(write(f,&nctg,sizeof(int)) == sizeof(int));
  CHECK(write(f,perm,sizeof(int)*nctg) == (ssize_t) (sizeof(int)*nctg));
  close(f);
}

static void f01_gix_reader(void)
{ char stem[400], path[512], err[600];
  int  perm[5] = { 3, 0, 4, 1, 2 }, bad[5] = { 3, 0, 3, 1, 2 }, got[5];

  snprintf(stem,sizeof(stem),"%s/f01_%d",tmpdir(),(int) getpid());
  snprintf(path,sizeof(path),"%s.gix",stem);

  write_gix(stem,5,perm);
  CHECK_EQ_I(fo_read_gix_perm(stem,5,got,err,sizeof(err)),1);
  CHECK(memcmp(got,perm,sizeof(perm)) == 0);

  CHECK_MSG(fo_read_gix_perm(stem,6,got,err,sizeof(err)) == 0,
            "a GDB with 6 contigs must not accept an index of 5");
  CHECK(strstr(err,"unexpected header") != NULL);

  write_gix(stem,5,bad);
  CHECK_MSG(fo_read_gix_perm(stem,5,got,err,sizeof(err)) == 0,
            "a repeated contig is not a permutation");
  CHECK(strstr(err,"not a permutation") != NULL);

  unlink(path);
  CHECK(fo_read_gix_perm(stem,5,got,err,sizeof(err)) == 0);
  CHECK(strstr(err,"cannot open") != NULL);
}

// ---------------------------------------------------------------------------
// F-02  parts and thread ranges against FastGA's own loops
// ---------------------------------------------------------------------------

// FastGA.c main, "Compute GDB split into NTHREADS parts" (178afcb), with
// gdb1->contigs[Perm1[x]].clen as clen[x].
static int ref_parts(int n, const int64_t *clen, int64_t npost, int T, int *split)
{ int64_t cum, t;
  int     p, r, x;

  split[0] = 0;
  p = 0;
  r = T;
  t = npost/T;
  cum = clen[0];
  for (x = 1; x < n; x++)
    { if (cum >= t && x >= r)
        { p += 1;
          split[p] = x;
          t = (npost*(p+1))/T;
          r += T;
        }
      cum += clen[x];
    }
  split[p+1] = n;
  return p+1;
}

// RSDsort.c rmsd_sort (178afcb), the thread split only: part[] is the panel
// of seed bytes per contig (zero outside the current part), asize their sum.
typedef struct { int beg, end; } Rng;
static int ref_rmsd_split(int nparts, const int64_t *part, int64_t asize, int nthreads, Rng *parms)
{ int64_t sum, thr;
  int     x, n, beg;

  n   = 0;
  thr = asize / nthreads;
  sum = 0;
  for (x = 0; x < nparts; x++)
    if (part[x] > 0)
      break;
  beg = x;
  for (; x < nparts; x++)
    if (part[x] > 0)
      { sum += part[x];
        if (sum >= thr)
          { parms[n].end = x+1;
            parms[n].beg = beg;
            n  += 1;
            thr = (asize * (n+1))/nthreads;
            beg = x+1;
          }
      }
  return n;
}

static void f02_ranges_match_fastga(void)
{ int trial, mism = 0, cases = 0;

  for (trial = 0; trial < 3000; trial++)
    { int      n = 1 + rndi(trial < 1500 ? 40 : 400);
      int      T = (int []) { 1, 2, 3, 8, 16, 32 }[rndi(6)];
      int64_t  swide = 8 + rndi(20);
      int64_t *clen  = malloc(n*sizeof(int64_t));
      int64_t *seeds = malloc(n*sizeof(int64_t));
      int64_t *panel = calloc(n,sizeof(int64_t));
      int     *split = malloc((T+2)*sizeof(int));
      int     *rs    = malloc(n*sizeof(int));
      int     *rb    = malloc(n*sizeof(int));
      int     *re    = malloc(n*sizeof(int));
      Rng     *ref   = malloc((n+T+1)*sizeof(Rng));
      int64_t  npost = 0;
      int      x, i, np, nref, nours, nparts_ours;

      for (x = 0; x < n; x++)            // longest first, as in FastGA's order
        { clen[x] = (x == 0) ? 1 + rndi(1000000) : 1 + rndi((int) clen[x-1]);
          npost  += clen[x];
          seeds[x] = (rndi(4) == 0) ? 0 : rndi(5000);
        }

      nref = 0;
      np = ref_parts(n,clen,npost,T,split);
      for (i = 0; i < np; i++)
        { int64_t asize = 0;
          memset(panel,0,n*sizeof(int64_t));
          for (x = split[i]; x < split[i+1]; x++)
            { panel[x] = seeds[x]*swide;
              asize   += panel[x];
            }
          if (asize > 0)
            nref += ref_rmsd_split(n,panel,asize,T,ref+nref);
        }

      nours = fo_compute_ranges(n,clen,npost,seeds,swide,T,rs,rb,re,&nparts_ours);
      cases += nref;
      if (nparts_ours != np || nours != nref)
        mism += 1;
      else
        for (i = 0; i < nref; i++)
          if (rb[i] != ref[i].beg || re[i] != ref[i].end)
            { mism += 1; break; }
      for (i = 0; i < nours; i++)        // range_start agrees with the ranges
        for (x = rb[i]; x < re[i]; x++)
          if (rs[x] != rb[i]) { mism += 1; i = nours; break; }

      free(clen); free(seeds); free(panel); free(split); free(rs); free(rb); free(re); free(ref);
    }
  CHECK_MSG(mism == 0,"%d of 3000 random inputs split differently from FastGA",mism);
  CHECK_MSG(cases > 10000,"only %d ranges compared",cases);
}

// ---------------------------------------------------------------------------
// F-03  nearness helpers
// ---------------------------------------------------------------------------

static void f03_nearness(void)
{ // row 5 has earlier partners at positions 0, 2 and 4
  int     n = 6;
  uint8_t has[36] = { 0 };

  has[5*6+0] = has[5*6+2] = has[5*6+4] = 1;
  CHECK_EQ_I(fo_count_earlier(has,n,5),3);
  CHECK_EQ_I(fo_nearness(has,n,5,4),0);
  CHECK_EQ_I(fo_nearness(has,n,5,2),1);
  CHECK_EQ_I(fo_nearness(has,n,5,0),2);
  CHECK_EQ_I(fo_kth_earlier(has,n,5,0),4);
  CHECK_EQ_I(fo_kth_earlier(has,n,5,2),0);
  CHECK_EQ_I(fo_kth_earlier(has,n,5,3),-1);
  CHECK_EQ_I(fo_count_earlier(has,n,3),0);
}

// ---------------------------------------------------------------------------
// F-04  answers
// ---------------------------------------------------------------------------

static void f04_answers(void)
{ int bepos, k;

  { int a[3] = { FO_ANSWER_NOCALL, 777, FO_ANSWER_UNKNOWN };     // decisive before unknown
    CHECK_EQ_I(fo_from_answers(a,3,9,0,32767,&bepos,&k),1);
    CHECK_EQ_I(bepos,777); CHECK_EQ_I(k,1);
  }
  { int a[2] = { FO_ANSWER_UNKNOWN, 5 };                         // nearest not in yet
    CHECK_EQ_I(fo_from_answers(a,2,9,0,32767,&bepos,&k),0);
  }
  { int a[2] = { FO_ANSWER_NOCALL, FO_ANSWER_NOCALL };            // no call, not first in range
    CHECK_EQ_I(fo_from_answers(a,2,9,4,32767,&bepos,&k),-1);
  }
  { int a[1] = { FO_ANSWER_NOCALL };                             // no call, first of a later range
    CHECK_EQ_I(fo_from_answers(a,1,9,9,32767,&bepos,&k),1);
    CHECK_EQ_I(bepos,0); CHECK_EQ_I(k,-1);
  }
  CHECK_EQ_I(fo_from_answers(NULL,0,0,0,32767,&bepos,&k),1);    // first of the run's first range
  CHECK_EQ_I(bepos,32767); CHECK_EQ_I(k,-1);
  { int a[1] = { 0 };                                            // bepos 0 is a real value
    CHECK_EQ_I(fo_from_answers(a,1,9,0,32767,&bepos,&k),1);
    CHECK_EQ_I(bepos,0); CHECK_EQ_I(k,0);
  }
}

// ---------------------------------------------------------------------------
// F-05 / F-06  the rules against a sequential model of FastGA
// ---------------------------------------------------------------------------
//
// A world: n positions in FastGA's order, thread ranges, and for every
// ordered pair (x,b), b != x, whether it has same-strand seeds and, if so,
// whether it reaches the aligner and with which final bepos.  Self pair
// (x,x) exists or not; if it reaches the aligner its final bepos is its own
// (selfend[x]), otherwise the leftover passes through unchanged.
//
// FastGA walks each range in order -- rows x ascending, and in row x the
// partners b ascending, with (x,x) at b = x -- carrying one leftover value,
// which starts unset at each range.  The entry value of (x,x) is the
// leftover when the walk reaches it.

typedef struct {
  int      n, unset_main;
  int     *range_start;
  uint8_t *has;        // has[x*n+b]: pair (x,b) has same-strand seeds (any b != x)
  uint8_t *called;     // called[x*n+b]: it reached the aligner
  int     *bepos;      // bepos[x*n+b]: final bepos if called
  uint8_t *self_has, *self_called;
  int     *self_end;
} World;

static void world_free(World *w)
{ free(w->range_start); free(w->has); free(w->called); free(w->bepos);
  free(w->self_has); free(w->self_called); free(w->self_end);
}

static void world_random(World *w, int n, double p_has, double p_call, double p_self_call)
{ int x, b, nb;

  w->n = n;
  w->unset_main  = 32767;
  w->range_start = malloc(n*sizeof(int));
  w->has         = calloc((size_t) n*n,1);
  w->called      = calloc((size_t) n*n,1);
  w->bepos       = calloc((size_t) n*n,sizeof(int));
  w->self_has    = calloc(n,1);
  w->self_called = calloc(n,1);
  w->self_end    = calloc(n,sizeof(int));
  nb = 0;
  for (x = 0; x < n; x++)
    { if (x == 0 || rndi(4) == 0) nb = x;      // a new thread range starts here
      w->range_start[x] = nb;
      for (b = 0; b < n; b++)
        if (b != x && (rnd() % 1000) < p_has*1000)
          { w->has[x*n+b] = 1;
            if ((rnd() % 1000) < p_call*1000)
              { w->called[x*n+b] = 1;
                w->bepos[x*n+b]  = rndi(1 << 30);
              }
          }
      if (rndi(3) != 0)
        { w->self_has[x] = 1;
          if ((rnd() % 1000) < p_self_call*1000)
            { w->self_called[x] = 1;
              w->self_end[x] = rndi(1 << 30);
            }
        }
    }
}

// FastGA's sequential walk: the entry value of every self pair.
static void world_reference(const World *w, int *want)
{ int x, b, cur = 0, n = w->n;

  for (x = 0; x < n; x++)
    { if (w->range_start[x] == x)
        cur = fo_unset_value(x,w->unset_main);
      for (b = 0; b < n; b++)
        if (b == x)
          { if (w->self_has[x])
              { want[x] = cur;
                if (w->self_called[x]) cur = w->self_end[x];
              }
          }
        else if (w->has[x*n+b] && w->called[x*n+b])
          cur = w->bepos[x*n+b];
    }
}

typedef struct { int from_partner, unset_direct, from_rows, rows_self, rows_unset, max_rounds; } Cover;

// What the distributed rules derive once every answer is in: answers by
// nearness, then rounds over the L and S snapshots until no self pair is left.
// Returns the number of rounds.
static int world_rules(const World *w, int *got, Cover *cv)
{ int      n = w->n, x, b, k, rounds = 0, left = 0;
  uint8_t *earlier = calloc((size_t) n*n,1);
  int64_t *L  = calloc(n,sizeof(int64_t));
  int64_t *S  = calloc(n,sizeof(int64_t));
  int64_t *Ss = calloc(n,sizeof(int64_t));
  int     *ans = malloc((n+1)*sizeof(int));

  for (x = 0; x < n; x++)
    for (b = 0; b < x; b++)
      earlier[x*n+b] = w->has[x*n+b];
  for (x = 0; x < n; x++)
    for (b = 0; b < n; b++)
      if (b != x && w->has[x*n+b] && w->called[x*n+b])
        { int64_t key = fo_make_key(b,w->bepos[x*n+b]);
          if (key > L[x]) L[x] = key;
        }

  for (x = 0; x < n; x++)
    { int nn, r, bep, ku;

      if (!w->self_has[x]) continue;
      nn = fo_count_earlier(earlier,n,x);
      for (k = 0; k < nn; k++)
        { int y = fo_kth_earlier(earlier,n,x,k);
          CHECK(fo_nearness(earlier,n,x,y) == k);
          ans[k] = w->called[x*n+y] ? w->bepos[x*n+y] : FO_ANSWER_NOCALL;
        }
      r = fo_from_answers(ans,nn,x,w->range_start[x],w->unset_main,&bep,&ku);
      CHECK_MSG(r != 0,"all answers are in, so nothing may wait");
      if (r == 1)
        { got[x] = bep;
          if (ku >= 0) cv->from_partner += 1; else cv->unset_direct += 1;
          S[x] = w->self_called[x] ? FO_S_CALLED + (int64_t) w->self_end[x] : FO_S_NOCALL;
        }
      else
        { S[x] = FO_S_PENDING;
          left += 1;
        }
    }

  while (left > 0)
    { int progress = 0;

      memcpy(Ss,S,n*sizeof(int64_t));         // what a round's reduction sees
      rounds += 1;
      for (x = 0; x < n; x++)
        if (S[x] == FO_S_PENDING)
          { int bep, row, from;
            if (fo_from_rows(L,Ss,x,w->range_start[x],w->unset_main,&bep,&row,&from))
              { got[x] = bep;
                if (row < 0) cv->rows_unset += 1;
                else { cv->from_rows += 1; if (from == row) cv->rows_self += 1; }
                S[x] = w->self_called[x] ? FO_S_CALLED + (int64_t) w->self_end[x] : FO_S_NOCALL;
                left -= 1;
                progress = 1;
              }
          }
      CHECK_MSG(progress,"a round must resolve at least one self pair");
      if (!progress) break;
    }
  if (rounds > cv->max_rounds) cv->max_rounds = rounds;

  free(earlier); free(L); free(S); free(Ss); free(ans);
  return rounds;
}

static void f05_rules_match_sequential_fastga(void)
{ Cover cv = { 0 };
  int   trial, wrong = 0, compared = 0;

  for (trial = 0; trial < 4000; trial++)
    { World w;
      int   n = 1 + rndi(30), x;
      // sparse worlds make the rare fallbacks common
      double p_has  = (trial % 3 == 0) ? 0.05 : ((trial % 3 == 1) ? 0.25 : 0.7);
      double p_call = (trial % 2 == 0) ? 0.2 : 0.8;
      int   *want, *got;

      world_random(&w,n,p_has,p_call,(trial % 4 == 0) ? 0.1 : 0.6);
      want = malloc(n*sizeof(int));
      got  = malloc(n*sizeof(int));
      world_reference(&w,want);
      world_rules(&w,got,&cv);
      for (x = 0; x < n; x++)
        if (w.self_has[x])
          { compared += 1;
            if (got[x] != want[x]) wrong += 1;
          }
      free(want); free(got);
      world_free(&w);
    }
  CHECK_MSG(wrong == 0,"%d of %d self pairs got a different value from FastGA's walk",
            wrong,compared);
  // the test must actually reach every rule, or it proves less than it claims
  CHECK_MSG(cv.from_partner > 1000,"from an earlier partner: %d",cv.from_partner);
  CHECK_MSG(cv.unset_direct > 100,"unset, first in range: %d",cv.unset_direct);
  CHECK_MSG(cv.from_rows > 100,"from an earlier row: %d",cv.from_rows);
  CHECK_MSG(cv.rows_self > 20,"from an earlier row's self pair: %d",cv.rows_self);
  CHECK_MSG(cv.rows_unset > 20,"unset after walking the rows: %d",cv.rows_unset);
  CHECK_MSG(cv.max_rounds >= 3,"longest chain of rounds: %d",cv.max_rounds);
  printf("    F-05: %d self pairs; partner %d, unset %d, rows %d (self %d, unset %d), "
         "max rounds %d\n",compared,cv.from_partner,cv.unset_direct,cv.from_rows,
         cv.rows_self,cv.rows_unset,cv.max_rounds);
}

// A chain: one range, rows 0..k, each row only its self pair, row 0's
// self pair reaches the aligner and the others do not.  Each self pair's
// value comes from the one before it, so k rounds are needed.
static void f06_chain_rounds(void)
{ int k;

  for (k = 1; k <= 12; k++)
    { World w;
      Cover cv = { 0 };
      int   n = k+1, x, rounds;
      int  *want = malloc(n*sizeof(int)), *got = malloc(n*sizeof(int));

      world_random(&w,n,0,0,0);
      for (x = 0; x < n; x++)
        { w.range_start[x] = 0;
          w.self_has[x]    = 1;
          w.self_called[x] = (x == 0);
          w.self_end[x]    = 4242;
        }
      world_reference(&w,want);
      rounds = world_rules(&w,got,&cv);
      CHECK_EQ_I(rounds,k);
      for (x = 1; x < n; x++)
        CHECK_EQ_I(got[x],4242);
      CHECK_EQ_I(got[0],32767);
      free(want); free(got);
      world_free(&w);
    }
}

int main(void)
{ RUN(f01_gix_reader);
  RUN(f02_ranges_match_fastga);
  RUN(f03_nearness);
  RUN(f04_answers);
  RUN(f05_rules_match_sequential_fastga);
  RUN(f06_chain_rounds);
  return t_report("test_fastga_order");
}
