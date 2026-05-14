#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <sys/resource.h>

#include "libfastk.h"
#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "fastga_pipeline.h"
#include "kmer_adapter.h"

// constants.h (pulled in via kmer_adapter.h -> record.h) defines KMER as the
// integer literal 40. FastGA's pipeline uses KMER as a runtime int variable.
// Drop the macro here so `extern int KMER` and `vlcp[KMER+1]` parse correctly.
#undef KMER

// rmsd_sort lives in lib/RSDsort.c. Its 8th arg is `Range *` but Range is
// a private typedef there; we pass a `Range` of identical layout (defined
// below in this file) and cast through void *.
extern int rmsd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int nparts, int64 *part, int nthreads, void *parms);

#ifndef PTR_SIZE
#define PTR_SIZE  sizeof(void *)
#endif
#ifndef OVL_SIZE
#define OVL_SIZE  sizeof(Overlap)
#endif
#ifndef EXO_SIZE
#define EXO_SIZE  (OVL_SIZE - PTR_SIZE)
#endif

#undef  CALL_ALIGNER
#define CALL_ALIGNER

#undef  BOX_ELIM
#define BOX_ELIM

#define MAX_INT64    0x7fffffffffffffffll
#define TSPACE       100
#define BUCK_SHIFT   6
#define BUCK_WIDTH   64
#define BUCK_ANTI    128
#define BOX_FUZZ     10

#define POST_BUF_LEN  0x1000
#define POST_BUF_MASK 0x0fff

#define MEMORY 4000

#define ELIMINATED   0x4
#define OWNS_MEMORY  0x8
#define RESET_FLAGS  0x3

// Redirect Kmer_Stream operations to our adapter
#define First_Kmer_Entry(T)    adapter_First_Kmer_Entry((KmerStreamAdapter *)(T))
#define Next_Kmer_Entry(T)     adapter_Next_Kmer_Entry((KmerStreamAdapter *)(T))
#define GoTo_Kmer_Index(T, i)  adapter_GoTo_Kmer_Index((KmerStreamAdapter *)(T), (i))
#define Clone_Kmer_Stream(T)   ((Kmer_Stream *)adapter_Clone_Kmer_Stream((KmerStreamAdapter *)(T)))
#define Free_Kmer_Stream(T)    adapter_Free_Kmer_Stream((KmerStreamAdapter *)(T))
#define Current_Entry(T, ent)  adapter_Current_Entry((KmerStreamAdapter *)(T), (ent))

// ============================================================================
// Externs (set by pga-merge.c)
// ============================================================================

extern int    NTHREADS;
extern int    NPARTS;
extern int    SELF;

// Phase 4b: per-destination-rank seed routing.
// If RankForPair is non-NULL, new_self_merge_thread routes seeds by
// RankForPair[icont * NCONTS + jcont] instead of Select[icont]. NUNITS is
// the number of slots in N_Units / C_Units per local thread (= world_size
// in distributed mode, = NPARTS in single-rank mode).
extern int8_t *RankForPair;
extern int     NUNITS;
extern int    SOFT_MASK;
extern int    VERBOSE;
extern int    KMER;

extern int    FREQ;
extern int    CHAIN_BREAK;
extern int    CHAIN_MIN;
extern int    ALIGN_MIN;
extern double ALIGN_RATE;
extern char  *Prog_Name;

extern int    IBYTE, JBYTE;
extern int    ICONT, JCONT;
extern int    IPOST, JPOST;
extern int    ISIGN, JSIGN;
extern int    KBYTE, CBYTE, LBYTE, PAYOFF;
extern int    ESHIFT;

extern int    NCONTS;
extern int64  AMXPOS, BMXPOS, MAXDAG;
extern int    DBYTE;

extern int   *Select;
extern int   *IDBsplit;
extern int   *Perm1, *Perm2;
extern char  *SORT_PATH;
extern char  *PAIR_NAME;
extern char  *ALGN_UNIQ;
extern char  *ALGN_PAIR;

// ============================================================================
// Types
// ============================================================================

// IOBuffer typedef lives in fastga_pipeline.h (shared with the driver).

extern IOBuffer *N_Units;
extern IOBuffer *C_Units;

typedef struct {
    int   beg;
    int   end;
    int64 off;
} Range;

// Post_List typedef lives in fastga_pipeline.h.

typedef struct {
    Kmer_Stream *T1;
    Kmer_Stream *T2;
    Post_List   *P1;
    Post_List   *P2;
    int          flip;
    int          tid;
    int          pbeg, pend;
    uint8       *cache;
    IOBuffer    *nunit;
    IOBuffer    *cunit;
    int64        nhits;
    int64        tseed;
} SP;

typedef struct {
    int       in;
    int       swide;
    int       comp;
    int       inum;
    GDB      *gdb1;
    GDB      *gdb2;
    int64    *buck;
    uint8    *buffer;
    uint8    *sarr;
    Range    *range;
} RP;

typedef struct {
    int       tid;
    int       swide;
    int       comp;
    int64    *panel;
    uint8    *sarr;
    Range    *range;
    GDB       gdb1;
    GDB       gdb2;
    FILE     *ofile;
    FILE     *tfile;
    int64     nhits;
    int64     nlass;
    int64     nlive;
    int64     nlcov;
    int64     nmemo;
    int       tmaxl;
} TP;

typedef struct {
    int         tid;
    GDB        *gdb1, *gdb2;
    FILE       *ofile;
    FILE       *tfile;
    int64       nhits;
    int64       nlass;
    int64       nlive;
    int64       nlcov;
    int64       nmemo;
    Work_Data  *work;
    Align_Spec *spec;
    Alignment   align;
    Overlap     ovl;
} Contig_Bundle;

typedef struct {
    FILE   *stream;
    void   *block;
    void   *ptr;
    void   *top;
    int64   count;
} IO_block;

// ============================================================================
// LCP tables for prefix matching
// ============================================================================

static int cbyte[41] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
                         3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
                         6, 6, 6, 6, 7 };

static int mbyte[41] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03,
                         0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03,
                         0xc0, 0x30, 0x0c, 0x03, 0xc0 };

// ============================================================================
// Forward declarations
// ============================================================================

static int  entwine(Path *jpath, uint8 *jtrace, Path *kpath, uint8 *ktrace, int *where, int show);
static void align_contigs(uint8 *beg, uint8 *end, int swide, int ctg1, int ctg2, Contig_Bundle *pair);
static int  ALIGN_SORT(const void *l, const void *r);
static int  SORT_MAP(const void *x, const void *y);
static void ovl_reload(IO_block *in, int64 bsize);
static void maheap(int s, Overlap **heap, int hsize);

// ============================================================================
// new_self_merge_thread
// ============================================================================

void *new_self_merge_thread(void *args)
{ SP *parm = (SP *) args;
  int tid          = parm->tid;
  uint8 *cache     = parm->cache;
  IOBuffer  *nunit = parm->nunit;
  IOBuffer  *cunit = parm->cunit;
  Kmer_Stream *T1  = parm->T1;

  int64   tbeg, tend;
  int     cpre;
  uint8  *ctop, *suf1;
  int     kbyte, kfreq;
  int     eorun, plen, mlen;
  uint8  *rcur, *rend;
  uint8  *vlcp[KMER+1];
  uint8  *low, *hgh, *top;
  int64   icont;
  uint8  *iptr = (uint8 *) (&icont);
  int     qcnt, pcnt;
  int64   nhits, tseed;

  { int j;
    for (j = 0; j < NUNITS; j++) {
      nunit[j].bend = nunit[j].bufr + (1000000-(IBYTE+JBYTE+1));
      nunit[j].btop = nunit[j].bufr;
      cunit[j].bend = cunit[j].bufr + (1000000-(IBYTE+JBYTE+1));
      cunit[j].btop = cunit[j].bufr;
    }
  }

  ctop  = cache;
  nhits = 0;
  tseed = 0;
  kbyte = T1->pbyte;
  kfreq = FREQ * kbyte;
  icont = 0;
  First_Kmer_Entry(T1);

  if (tid != 0)
    GoTo_Kmer_Index(T1,T1->index[(parm->pbeg<<8) | 0xff]);
  tend = T1->index[(parm->pend<<8) | 0xff];
  tbeg = T1->cidx;

  mlen = KMER+1;
  plen = 12;
  vlcp[12] = rcur = rend = cache;
  eorun = 0;
  qcnt = -1;

  for (suf1 = ctop; 1; suf1 += kbyte)
    { if (suf1 >= ctop)
        { uint8 *cp;
          int    i;

          if (VERBOSE && tid == 0)
            { if (tbeg == tend)
                pcnt = 100;
              else
                pcnt = ((T1->cidx - tbeg) * 100) / (tend-tbeg);
              if (pcnt > qcnt)
                { fprintf(stderr,"\r    Completed %3d%%",pcnt);
                  fflush(stderr);
                }
              qcnt = pcnt;
            }

          if (T1->cidx >= tend)
            break;

          cpre = T1->cpre;
          for (cp = cache; T1->cpre == cpre; cp += kbyte)
            { memcpy(cp,T1->csuf,kbyte);
              Next_Kmer_Entry(T1);
            }
          ctop = cp;
          ctop[LBYTE] = 11;

          plen = 12;
          vlcp[plen] = rcur = rend = suf1 = cache;
          rend += kbyte;
          plen = rend[LBYTE];
          for (i = rcur[LBYTE]; i <= plen; i++)
            vlcp[i] = rcur;
          eorun = (plen <= 11);
        }
      else
        { int i;
          if (eorun)
            plen = suf1[LBYTE];
          rend += kbyte;
          if (rend[LBYTE] < plen)
            eorun = 1;
          else if (rend[LBYTE] == plen)
            eorun = 0;
          else
            { rcur = rend-kbyte;
              for (i = plen+1; i <= rend[LBYTE]; i++)
                vlcp[i] = rcur;
              eorun = 0;
              plen = rend[LBYTE];
            }
        }

      low = vlcp[plen];
      hgh = rend;
      top = low+kfreq;
      if (!eorun)
        { do
            { hgh += kbyte;
              if (hgh > top)
                break;
            }
          while (hgh[LBYTE] >= plen);
        }

      if (hgh >= top) continue;
      if (SOFT_MASK) mlen = plen;
      if (suf1[CBYTE] >= mlen) continue;

      { int       idest, isign, jsign;
        IOBuffer *ou;
        uint8    *p, *pay1, *btop;
        int64     jcont;
        uint8    *jptr = (uint8 *) (&jcont);
        int64     jcont_mask = ((int64) 1 << (8 * JCONT - 1)) - 1;

        pay1 = suf1+PAYOFF;
        isign = (pay1[ISIGN] & 0x80);
        if (isign) pay1[ISIGN] &= 0x7f;
        memcpy(iptr,pay1+IPOST,ICONT);
        idest = Select[icont];   // fallback for non-distributed mode

        for (p = low+PAYOFF; p < hgh; p += kbyte)
          { if (p == pay1) continue;
            if (p[-2] >= mlen) continue;
            jsign = (p[ISIGN] & 0x80);

            // Phase 4b: route by RankForPair[icont, jcont] if available.
            // Otherwise keep original Select[icont] routing.
            if (RankForPair != NULL) {
              jcont = 0;
              memcpy(jptr, p+JPOST, JCONT);
              jcont &= jcont_mask;   // strip jsign bit
              idest = RankForPair[icont * (int64) NCONTS + jcont];
            }

            if (isign == jsign)
              ou = nunit + idest;
            else
              ou = cunit + idest;
            btop = ou->btop;
            *btop++ = plen;
            memcpy(btop,pay1,IBYTE);
            btop += IBYTE;
            memcpy(btop,p,JBYTE);
            btop += JBYTE;

            nhits += 1;
            tseed += plen;
            ou->buck[icont] += 1;

            if (btop >= ou->bend)
              { if (write(ou->file,ou->bufr,btop-ou->bufr) < 0)
                  { fprintf(stderr,"%s: IO write to file %s/%s.%d.%c failed\n",
                                   Prog_Name,SORT_PATH,PAIR_NAME,
                                   ou->inum,(isign == jsign) ? 'N' : 'C');
                    exit (1);
                  }
                ou->btop = ou->bufr;
              }
            else
              ou->btop = btop;
          }
        if (isign) pay1[ISIGN] |= 0x80;
      }
    }

  { int j;
    for (j = 0; j < NUNITS; j++)
      { if (nunit[j].btop > nunit[j].bufr)
          if (write(nunit[j].file,nunit[j].bufr,nunit[j].btop-nunit[j].bufr) < 0)
            { fprintf(stderr,"%s: IO write to file %s/%s.%d.N failed\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,nunit[j].inum);
              exit (1);
            }
        if (cunit[j].btop > cunit[j].bufr)
          if (write(cunit[j].file,cunit[j].bufr,cunit[j].btop-cunit[j].bufr) < 0)
            { fprintf(stderr,"%s: IO write to file %s/%s.%d.C failed\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,cunit[j].inum);
              exit (1);
            }
      }
  }

  parm->nhits = nhits/2;
  parm->tseed = tseed/2;
  return (NULL);
}

// ============================================================================
// self_adaptamer_merge
// ============================================================================

void self_adaptamer_merge(Kmer_Stream *T1, Post_List *P1, int64 g1len)
{ SP         parm[NTHREADS];
  pthread_t  threads[NTHREADS];
  uint8     *cache;
  int64      nhits, tseed;
  int        i;

  if (VERBOSE)
    { fprintf(stderr,"  Starting adaptive seed merge\n");
      fflush(stderr);
    }

  { uint8 *ent;
    int    t;
    int64  p;

    ent = Current_Entry(T1,NULL);
    parm[0].pbeg = 0;
    for (t = 1; t < NTHREADS; t++)
      { p = (T1->nels * t) / NTHREADS;
        GoTo_Kmer_Index(T1,p);
        if (p >= T1->nels)
          parm[t].pbeg = 0xffff;
        else
          { ent = Current_Entry(T1,ent);
            parm[t].pbeg = (T1->cpre >> 8);
          }
      }
    for (t = 0; t < NTHREADS-1; t++)
      parm[t].pend = parm[t+1].pbeg;
    parm[NTHREADS-1].pend = 0xffff;
  }

  parm[0].T1 = T1;
  parm[0].P1 = P1;
  for (i = 1; i < NTHREADS; i++)
    { parm[i].T1 = Clone_Kmer_Stream(T1);
      parm[i].P1 = P1;
    }

  cache = Malloc(NTHREADS*(P1->maxp+1)*KBYTE,"Allocating cache");
  if (cache == NULL) exit (1);

  for (i = 0; i < NTHREADS; i++)
    { IOBuffer *nu, *cu;
      parm[i].tid   = i;
      parm[i].cache = cache + i * (P1->maxp+1) * KBYTE;
      parm[i].nunit = nu = N_Units + i * NUNITS;
      parm[i].cunit = cu = C_Units + i * NUNITS;
      bzero(nu[0].buck,sizeof(int64)*NCONTS);
      bzero(cu[0].buck,sizeof(int64)*NCONTS);
    }

  for (i = 1; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,new_self_merge_thread,parm+i);
  new_self_merge_thread(parm);
  for (i = 1; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  if (VERBOSE)
    { fprintf(stderr,"\r    Completed 100%%\n");
      fflush(stderr);
    }

  free(cache);
  for (i = NTHREADS-1; i >= 1; i--)
    Free_Kmer_Stream(parm[i].T1);
  Free_Kmer_Stream(T1);

  nhits = tseed = 0;
  for (i = 0; i < NTHREADS; i++)
    { nhits += parm[i].nhits;
      tseed += parm[i].tseed;
    }

  if (VERBOSE)
    { fprintf(stderr,"\n  Total seeds = %lld, ave. len = %.1f, seeds per G1 position = %.1f\n",
                     nhits,(1.*tseed)/nhits,(1.*nhits)/g1len);
      fflush(stderr);
    }
}

// ============================================================================
// entwine
// ============================================================================

static int entwine(Path *jpath, uint8 *jtrace, Path *kpath, uint8 *ktrace, int *where, int show)
{ int ac, b2, y2, yp, ae;
  int i, j, k;
  int num, den, min;

  (void) show;
  *where = -1;

  y2 = jpath->bbpos;
  b2 = kpath->bbpos;
  j  = jpath->abpos/TSPACE;
  k  = kpath->abpos/TSPACE;
  ac = k*TSPACE;
  j = 1 + 2*(k-j);
  k = 1;

  for (i = 1; i < j; i += 2)
    y2 += jtrace[i];

  if (j == 1)
    yp = y2 + (jtrace[j] * (kpath->abpos - jpath->abpos)) / (ac+TSPACE - jpath->abpos);
  else
    yp = y2 + (jtrace[j] * (kpath->abpos - ac)) / TSPACE;

  num = b2-yp;
  den = 1;
  min = num;

  ae = jpath->aepos;
  if (ae > kpath->aepos) ae = kpath->aepos;

  for (ac += TSPACE; ac < ae; ac += TSPACE)
    { y2 += jtrace[j];
      b2 += ktrace[k];
      j += 2;
      k += 2;
      i = b2-y2;
      num += i;
      den += 1;
      if (min < 0 && min < i)
        { if (i >= 0) min = 0; else min = i; }
      else if (min > 0 && min > i)
        { if (i <= 0) min = 0; else min = i; }
      if (i == 0) *where = ac;
    }

  ac -= TSPACE;
  if (ae == jpath->aepos)
    { y2 = jpath->bepos;
      if (kpath->aepos >= ac)
        b2 += (ktrace[k] * (ae - ac)) / TSPACE;
      else
        b2 += (ktrace[k] * (ae - ac)) / (kpath->aepos - ac);
    }
  else
    { b2 = kpath->bepos;
      if (jpath->aepos >= ac)
        y2 += (jtrace[j] * (ae - ac)) / TSPACE;
      else
        y2 += (jtrace[j] * (ae - ac)) / (jpath->aepos - ac);
    }

  i = b2-y2;
  num += i;
  den += 1;
  if (min < 0 && min < i)
    { if (i >= 0) min = 0; else min = i; }
  else if (min > 0 && min > i)
    { if (i <= 0) min = 0; else min = i; }

  (void) den;
  return (min);
}

// ============================================================================
// align_contigs
// ============================================================================

static void align_contigs(uint8 *beg, uint8 *end, int swide, int ctg1, int ctg2,
                          Contig_Bundle *pair)
{ Overlap    *ovl   = &(pair->ovl);
  int         comp  = (ovl->flags != 0);
  Work_Data  *work  = pair->work;
  Align_Spec *spec  = pair->spec;
  Alignment  *align = &(pair->align);
  Path       *path  = align->path;
  FILE       *ofile = pair->ofile;
  FILE       *tfile = pair->tfile;

  uint8 *b, *m, *e;
  int    alnMin;
  double alnRate;
  int64  nhit, nlas, nmem, nliv, ncov;
  int64  alen, blen, mlen;
  int64  aoffset, doffset;
  int    new, aux;
  int64  ndiag, cdiag;
  uint8 *_ndiag = (uint8 *) (&ndiag);
  int    self;
  int64  ipost, apost;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_apost = (uint8 *) (&apost);

  ctg1 = Perm1[ctg1];
  ctg2 = Perm2[ctg2];

  if (pair->gdb1->contigs[ctg1].boff < 0 || pair->gdb2->contigs[ctg2].boff < 0)
    return;

  alnMin  = ALIGN_MIN - 50;
  alnRate = ALIGN_RATE + .05;

  // Zero int64 slots before partial-byte memcpy fills only their LSBs.
  // (FastGA.c:3016-3018)
  ndiag = 0;
  ipost = 0;
  apost = 0;

  blen   = pair->gdb2->contigs[ctg2].clen;
  alen   = pair->gdb1->contigs[ctg1].clen;
  mlen   = alen+blen;

  nhit = nlas = nmem = nliv = ncov = 0;

  if (SELF && ctg1 == ctg2 && !comp)
    self = 1;
  else
    self = 0;

  doffset = alen - MAXDAG;
  aoffset = alen - AMXPOS;

  b = e = beg + (DBYTE+2);
  memcpy(_ndiag,e,DBYTE);
  cdiag = ndiag;
  while (ndiag == cdiag && e < end)
    { e += swide;
      memcpy(_ndiag,e,DBYTE);
    }
  new = 1;

  while (1)
    { m = e;
      aux = 0;
      while (ndiag == cdiag+1 && e < end)
        { e += swide;
          memcpy(_ndiag,e,DBYTE);
          aux = 1;
        }

      if (new || aux)
        { int    go, lcp, wch, mix, cov;
          int64  ahgh, alow, amid, alast;
          int64  anti, eant;
          int    dgmin, dgmax, dg;
          uint8 *s, *t;

          alast = -1;
          e -= DBYTE;
          m -= DBYTE;

          s = b-DBYTE;
          memcpy(_ipost,s,DBYTE);
          t = m;
          if (aux) memcpy(_apost,t,DBYTE);
          else apost = MAX_INT64;

          dgmin = 2*BUCK_WIDTH;
          dgmax = 0;
          ahgh = -CHAIN_BREAK;
          if (apost < ipost) alow = apost;
          else alow = ipost;
          cov = 0;
          go = 1;
          mix = 0;

          while (go)
            { if (apost < ipost)
                { lcp = t[-2];
                  dg = t[-1] + BUCK_WIDTH;
                  anti = apost;
                  t += swide;
                  if (t >= e) apost = MAX_INT64;
                  else memcpy(_apost,t,DBYTE);
                  wch = 0x2;
                }
              else
                { lcp = s[-2];
                  dg = s[-1];
                  anti = ipost;
                  s += swide;
                  if (s >= m)
                    { if (s > m) go = 0;
                      else ipost = MAX_INT64;
                    }
                  else memcpy(_ipost,s,DBYTE);
                  wch = 0x1;
                }
              lcp <<= 1;

              if (anti < ahgh + CHAIN_BREAK)
                { int64 cps;
                  cps = anti + lcp;
                  if (cps > ahgh)
                    { if (anti >= ahgh) cov += lcp;
                      else cov += cps-ahgh;
                      ahgh = cps;
                    }
                  mix |= wch;
                  if (dg < dgmin) dgmin = dg;
                  else if (dg > dgmax) dgmax = dg;
                }
              else
                { if (cov >= CHAIN_MIN && (mix != 1 || new))
                    { nhit += 1;

                      if (ctg1 != ovl->aread)
                        { if (Get_Contig(pair->gdb1,ctg1,NUMERIC,align->aseq) == NULL) exit (1);
                          align->alen = alen;
                          ovl->aread  = ctg1;
                          if (comp) Complement_Seq(align->aseq,align->alen);
                        }
                      if (ctg2 != ovl->bread)
                        { if (Get_Contig(pair->gdb2,ctg2,NUMERIC,align->bseq) == NULL) exit (1);
                          align->blen = blen;
                          ovl->bread  = ctg2;
                        }

                      dgmin += (cdiag<<BUCK_SHIFT);
                      dgmax += (cdiag<<BUCK_SHIFT);
                      if (comp)
                        { dgmin += doffset;
                          dgmax += doffset;
                          alow  += aoffset;
                          ahgh  += aoffset;
                        }
                      else
                        { dgmin -= BMXPOS;
                          dgmax -= BMXPOS;
                        }

                      if (ahgh > alast)
                        { int rlen;
                          if (alow < alast) alow = alast;
                          ahgh -= BUCK_ANTI;
                          do {
                            amid = alow + BUCK_ANTI;
                            if (amid > ahgh)
                              { amid = ahgh;
                                if (amid + dgmin < 0)
                                  { dgmin = -amid;
                                    if (dgmin > dgmax) break;
                                  }
                              }
                            if (self)
                              { if (dgmin > 0)
                                  { if (Local_Alignment(align,work,spec,
                                                        dgmin,dgmax,amid,dgmin-1,-1)) exit (1); }
                                else if (dgmax < 0)
                                  { if (Local_Alignment(align,work,spec,
                                                        dgmin,dgmax,amid,-1,-(dgmax+1))) exit (1); }
                                else
                                  path->abpos = path->aepos = 0;
                              }
                            else
                              { if (Local_Alignment(align,work,spec,dgmin,dgmax,amid,-1,-1)) exit(1); }

                            rlen = path->aepos - path->abpos;
                            if (rlen >= alnMin && alnRate*rlen >= path->diffs)
                              { Compress_TraceTo8(ovl,0);
                                if (fwrite(ovl,OVL_SIZE,1,tfile) != 1)
                                  { fprintf(stderr,"%s: Cannot write overlap gather file %s/%s.%d.las\n",
                                           Prog_Name,SORT_PATH,ALGN_PAIR,pair->tid); exit (1); }
                                if (fwrite(ovl->path.trace,ovl->path.tlen,1,tfile) != 1)
                                  { fprintf(stderr,"%s: Cannot write overlap gather file %s/%s.%d.las\n",
                                           Prog_Name,SORT_PATH,ALGN_PAIR,pair->tid); exit (1); }
                                nlas += 1;
                                nmem += path->tlen + OVL_SIZE;
                              }
                            if (comp)
                              eant = mlen-(path->abpos+path->bbpos);
                            else
                              eant = path->aepos+path->bepos;
                            if (eant <= alow) alow = amid;
                            else alow = eant;
                          } while (alow < ahgh);
                          alast = alow;
                        }
                    }
                  if (go)
                    { cov = lcp;
                      ahgh = anti + lcp;
                      mix = wch;
                      alow = anti;
                      dgmin = dgmax = dg;
                    }
                }
            }
          e += DBYTE;
          m += DBYTE;

          // Reset partial-byte int64s before the next chain reads into them.
          // (FastGA.c:3384)
          ipost = apost = 0;
        }

      if (e >= end) break;

      if (aux)
        { b = m;
          cdiag += 1;
          new = 0;
        }
      else
        { b = e;
          cdiag = ndiag;
          while (ndiag == cdiag && e < end)
            { e += swide;
              memcpy(_ndiag,e,DBYTE);
            }
          new = 1;
        }
    }

  if (nlas > 0)
    { void    *oblock;
      Overlap **perm;
      int      j, k, where, dist;
      Path     tpath;

      oblock = Malloc(nmem,"Allocating overlap block");
      perm   = Malloc(nlas*sizeof(Overlap *),"Allocating permutation array");
      if (oblock == NULL || perm == NULL) exit (1);

      rewind(tfile);
      if (fread(oblock,nmem,1,tfile) != 1)
        { fprintf(stderr,"\n%s: Cannot read overlap gather file %s/%s.%d.las\n",
                         Prog_Name,SORT_PATH,ALGN_PAIR,pair->tid); exit (1); }

      { void *off;
        off = oblock;
        for (j = 0; j < nlas; j++)
          { perm[j] = (Overlap *) off;
            off += OVL_SIZE + ((Overlap *) off)->path.tlen;
          }
      }

      qsort(perm,nlas,sizeof(Overlap *),ALIGN_SORT);

      for (j = nlas-1; j >= 0; j--)
        { Overlap *o = perm[j];
          Path *op = &(o->path);
          for (k = j+1; k < nlas; k++)
            { Overlap *w = perm[k];
              Path *wp = &(w->path);
              if (op->aepos <= wp->abpos) break;
              if (w->flags & ELIMINATED) continue;
              if (op->abpos == wp->abpos && op->bbpos == wp->bbpos)
                if (op->aepos == wp->aepos && op->bepos == wp->bepos)
                  { if (op->diffs < wp->aepos)
                      { w->flags |= ELIMINATED; continue; }
                    else
                      { o->flags |= ELIMINATED; break; }
                  }
                else
                  { if (op->aepos > wp->aepos)
                      { w->flags |= ELIMINATED; continue; }
                    else
                      { o->flags |= ELIMINATED; break; }
                  }
              else if (op->aepos == wp->aepos && op->bepos == wp->bepos)
                { if (op->abpos < wp->abpos)
                    { w->flags |= ELIMINATED; continue; }
                  else
                    { o->flags |= ELIMINATED; break; }
                }
            }
        }

      for (j = nlas-1; j >= 0; j--)
        { Overlap *o = perm[j];
          Path *op = &(o->path);
          if (o->flags & ELIMINATED) continue;
          for (k = j+1; k < nlas; k++)
            { Overlap *w = perm[k];
              Path *wp = &(w->path);
              uint8 *otrace, *wtrace;
              if (op->aepos <= wp->abpos) break;
              if (w->flags & ELIMINATED) continue;
              if (op->bepos <= wp->bbpos || op->bbpos >= wp->bepos) continue;

              if (o->flags & OWNS_MEMORY) otrace = (uint8 *) op->trace;
              else otrace = (uint8 *) (o+1);
              if (w->flags & OWNS_MEMORY) wtrace = (uint8 *) wp->trace;
              else wtrace = (uint8 *) (w+1);

              dist = entwine(op,otrace,wp,wtrace,&where,0);
              if (where != -1)
                { uint8 *ntrace;
                  int ocut, wcut, d, h, g;
                  ocut = 2 * (((where-op->abpos)-1)/TSPACE+1);
                  wcut = 2 * (((where-wp->abpos)-1)/TSPACE+1);
                  op->tlen = ocut + (wp->tlen-wcut);
                  ntrace = (uint8 *) Malloc(op->tlen,"Allocating new trace");
                  if (ntrace == NULL) exit (1);
                  d = 0; h = 0;
                  for (g = 0; g < ocut; g += 2)
                    { d += (ntrace[h] = otrace[g]);
                      ntrace[h+1] = otrace[g+1];
                      h += 2;
                    }
                  for (g = wcut; g < wp->tlen; g += 2)
                    { d += (ntrace[h] = wtrace[g]);
                      ntrace[h+1] = wtrace[g+1];
                      h += 2;
                    }
                  if (o->flags & OWNS_MEMORY) free(otrace);
                  if (w->flags & OWNS_MEMORY) free(wtrace);
                  op->diffs = d;
                  op->aepos = wp->aepos;
                  op->bepos = wp->bepos;
                  w->flags |= ELIMINATED;
                  o->flags |= OWNS_MEMORY;
                  op->trace = ntrace;
                  continue;
                }
#ifdef BOX_ELIM
              if (dist != 0)
                { if ((op->aepos - op->abpos) + BOX_FUZZ >= wp->aepos - wp->abpos)
                    { if (wp->aepos <= op->aepos+BOX_FUZZ && wp->bbpos >= op->bbpos-BOX_FUZZ &&
                          wp->bepos <= op->bepos+BOX_FUZZ)
                        { w->flags |= ELIMINATED; continue; }
                    }
                  else
                    { if (op->aepos <= wp->aepos+BOX_FUZZ && op->bbpos >= wp->bbpos-BOX_FUZZ &&
                          op->bepos <= wp->bepos+BOX_FUZZ && op->abpos >= wp->abpos-BOX_FUZZ)
                        { o->flags |= ELIMINATED; continue; }
                    }
                }
#endif
            }
        }

      nmem = 0;
      for (j = 0; j < nlas; j++)
        { Overlap *o = perm[j];
          int hasmem;
          if (o->flags & ELIMINATED) continue;
          hasmem = (o->flags & OWNS_MEMORY);
          o->flags &= RESET_FLAGS;
          if (fwrite( ((char *) o)+PTR_SIZE, EXO_SIZE, 1, ofile) != 1)
            { fprintf(stderr,"%s: Could not write to overlap block file %s/%s.%d.las\n",
                             Prog_Name,SORT_PATH,ALGN_UNIQ,pair->tid); exit (1); }
          if (hasmem)
            { if (fwrite(o->path.trace, o->path.tlen, 1, ofile) != 1)
                { fprintf(stderr,"%s: Could not write to overlap block file %s/%s.%d.las\n",
                                 Prog_Name,SORT_PATH,ALGN_UNIQ,pair->tid); exit (1); }
              free(o->path.trace);
            }
          else
            { if (fwrite( (char *) (o+1), o->path.tlen, 1, ofile) != 1)
                { fprintf(stderr,"%s: Could not write to overlap block file %s/%s.%d.las\n",
                                 Prog_Name,SORT_PATH,ALGN_UNIQ,pair->tid); exit (1); }
            }
          nliv += 1;
          ncov += o->path.aepos - o->path.abpos;
          nmem += EXO_SIZE + o->path.tlen;
        }
      rewind(tfile);
      free(perm);
      free(oblock);
    }
  else
    nmem = 0;

  pair->nhits += nhit;
  pair->nlass += nlas;
  pair->nlive += nliv;
  pair->nlcov += ncov;
  pair->nmemo += nmem;
}

static int ALIGN_SORT(const void *l, const void *r)
{ Overlap *ol = *((Overlap **) l);
  Overlap *or = *((Overlap **) r);
  return (ol->path.abpos - or->path.abpos);
}

// ============================================================================
// search_seeds
// ============================================================================

void *search_seeds(void *args)
{ TP *parm = (TP *) args;
  int      swide  = parm->swide;
  int      comp   = parm->comp;
  int64   *panel  = parm->panel;
  uint8   *sarray = parm->sarr;
  Range   *range  = parm->range;
  int      beg    = range->beg;
  int      end    = range->end;
  GDB     *gdb1   = &(parm->gdb1);
  GDB     *gdb2   = &(parm->gdb2);
  int      foffs  = swide-JCONT;
  FILE    *ofile  = parm->ofile;
  FILE    *tfile  = parm->tfile;

  int    icrnt;
  int64  jcrnt;
  uint8 *_jcrnt = (uint8 *) (&jcrnt);
  Contig_Bundle _pair, *pair = &_pair;
  uint8 *x, *e, *b;

  jcrnt = 0;
  pair->tid  = parm->tid;
  pair->gdb1 = gdb1;
  pair->gdb2 = gdb2;
  pair->align.aseq = New_Contig_Buffer(gdb1);
  pair->align.bseq = New_Contig_Buffer(gdb2);
  if (pair->align.bseq == NULL || pair->align.bseq == NULL) exit (1);
  pair->align.path = &(pair->ovl.path);
  if (comp)
    { pair->ovl.flags = COMP_FLAG;
      pair->align.flags = ACOMP_FLAG;
    }
  else
    { pair->ovl.flags = 0;
      pair->align.flags = 0;
    }
  pair->ovl.aread = -1;
  pair->ovl.bread = -1;
  pair->work = New_Work_Data();
  pair->spec = New_Align_Spec(1.-ALIGN_RATE,100,gdb1->freq,0);
  if (pair->work == NULL || pair->spec == NULL) exit (1);
  pair->ofile = ofile;
  pair->tfile = tfile;
  pair->nhits = 0;
  pair->nlass = 0;
  pair->nlive = 0;
  pair->nlcov = 0;
  pair->nmemo = 0;

  x = sarray + range->off;
  for (icrnt = beg; icrnt < end; icrnt++)
    { e = x + panel[icrnt];
      if (e > x)
        { memcpy(_jcrnt,x+foffs,JCONT);
          b = x;
          for (x += swide; x < e; x += swide)
            if (memcmp(_jcrnt,x+foffs,JCONT))
              { align_contigs(b,x,swide,icrnt,(int) jcrnt,pair);
                memcpy(_jcrnt,x+foffs,JCONT);
                b = x;
              }
          align_contigs(b,x,swide,icrnt,jcrnt,pair);
        }
    }

  Free_Align_Spec(pair->spec);
  Free_Work_Data(pair->work);
  free(pair->align.aseq-1);
  free(pair->align.bseq-1);

  parm->nhits += pair->nhits;
  parm->nlass += pair->nlass;
  parm->nlive += pair->nlive;
  parm->nlcov += pair->nlcov;
  parm->nmemo += pair->nmemo;
  return (NULL);
}

// ============================================================================
// reimport_thread
// ============================================================================

void *reimport_thread(void *args)
{ RP *parm = (RP *) args;
  int    swide  = parm->swide;
  int    in     = parm->in;
  int    comp   = parm->comp;
  uint8 *sarr   = parm->sarr;
  uint8 *bufr   = parm->buffer;
  int64 *buck   = parm->buck;

  int64  ipost, jpost, icont, jcont, band, anti;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_jpost = (uint8 *) (&jpost);
  uint8 *_icont = (uint8 *) (&icont);
  uint8 *_jcont = (uint8 *) (&jcont);
  uint8 *_band  = (uint8 *) (&band);
  uint8 *_anti  = (uint8 *) (&anti);

  int    iamt;
  uint8 *x;
  int    iolen, iunit, lcp;
  int64  diag, flag, mask;
  uint8 *bend, *btop, *b;

  iolen = 2*NPARTS*1000000;
  iunit = IBYTE + JBYTE + 1;

  iamt = read(in,bufr,iolen);
  if (iamt < 0)
    { fprintf(stderr,"%s: IO read error for file %s/%s.%d.%c\n",
                     Prog_Name,SORT_PATH,PAIR_NAME,parm->inum,comp?'C':'N');
      exit (1);
    }
  bend = bufr + iamt;
  if (bend-bufr < iolen) btop = bend;
  else btop = bend-iunit;
  b = bufr;

  // Zero the int64 slots before partial-byte memcpy fills only their LSBs.
  // (FastGA.c:2681-2684)
  ipost = 0;
  jpost = 0;
  icont = 0;
  jcont = 0;

  flag = (0x1ll << (8*JCONT-1));
  mask = flag-1;

  if (bend > bufr)
  while (1)
    { lcp = *b++;
      memcpy(_ipost,b,IPOST); b += IPOST;
      memcpy(_icont,b,ICONT); b += ICONT;
      memcpy(_jpost,b,JPOST); b += JPOST;
      memcpy(_jcont,b,JCONT); b += JCONT;
      jcont &= mask;

      x = sarr + swide * buck[icont]++;
      *x++ = lcp;
      if (comp)
        { diag = MAXDAG - (ipost + jpost);
          anti = AMXPOS - (ipost - jpost);
        }
      else
        { diag = BMXPOS + (ipost - jpost);
          anti = ipost + jpost;
        }
      band = (diag >> BUCK_SHIFT);
      *x++ = diag-(band<<BUCK_SHIFT);
      memcpy(x,_anti,DBYTE); x += DBYTE;
      memcpy(x,_band,DBYTE); x += DBYTE;
      memcpy(x,_jcont,JCONT); x += JCONT;

      if (b >= btop)
        { int ex = bend-b;
          memcpy(bufr,b,ex);
          bend = bufr+ex;
          iamt = read(in,bend,iolen-ex);
          if (iamt < 0)
            { fprintf(stderr,"%s: IO read error for file %s/%s.%d.%c\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,parm->inum,comp?'C':'N');
              exit (1);
            }
          bend += iamt;
          if (bend == bufr) break;
          if (bend-bufr < iolen) btop = bend;
          else btop = bend-iunit;
          b = bufr;
        }
    }
  close(in);
  return (NULL);
}

// ============================================================================
// Heap and sort helpers
// ============================================================================

#define MAPARE(lp,rp)				\
  if (lp->aread > rp->aread)			\
    bigger = 1;					\
  else if (lp->aread < rp->aread)		\
    bigger = 0;					\
  else if (lp->path.abpos > rp->path.abpos)	\
    bigger = 1;					\
  else if (lp->path.abpos < rp->path.abpos)	\
    bigger = 0;					\
  else if (lp > rp)				\
    bigger = 1;					\
  else						\
    bigger = 0;

static void maheap(int s, Overlap **heap, int hsize)
{ int      c, l, r;
  int      bigger;
  Overlap *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        bigger = 1;
      else
        { hr = heap[r];
          MAPARE(hr,hl)
        }
      if (bigger)
        { MAPARE(hs,hl)
          if (bigger) { heap[c] = hl; c = l; }
          else break;
        }
      else
        { MAPARE(hs,hr)
          if (bigger) { heap[c] = hr; c = r; }
          else break;
        }
    }
  if (c != s) heap[c] = hs;
}

static int SORT_MAP(const void *x, const void *y)
{ Overlap *ol = *((Overlap **) x);
  Overlap *or = *((Overlap **) y);
  int al = ol->aread, ar = or->aread;
  if (al != ar) return (al-ar);
  int pl = ol->path.abpos, pr = or->path.abpos;
  if (pl != pr) return (pl-pr);
  int bl = ol->bread, br = or->bread;
  if (bl != br) return (bl-br);
  int cl = COMP(ol->flags), cr = COMP(or->flags);
  if (cl != cr) return (cl-cr);
  if (ol < or) return (-1);
  else if (ol > or) return (1);
  else return (0);
}

// ============================================================================
// la_sort
// ============================================================================

void *la_sort(void *args)
{ TP *parm = (TP *) args;
  FILE *fid  = parm->ofile;
  int64 novl = parm->nlive;
  int64 size = parm->nmemo;
  void     *iblock, *off;
  Overlap **perm;
  int       j, tmaxl;

  if (novl == 0) return (NULL);

  iblock = Malloc(size+PTR_SIZE,"Allocating overlap block");
  perm   = Malloc(sizeof(Overlap *)*novl,"Allocating permutation array");
  if (iblock == NULL || perm == NULL) exit (1);
  iblock += PTR_SIZE;

  rewind(fid);
  if (fread(iblock,size,1,fid) != 1)
    { fprintf(stderr,"\n%s: Cannot not read overlap block file %s/%s.%d.las\n",
                     Prog_Name,SORT_PATH,ALGN_UNIQ,parm->tid); exit (1); }
  rewind(fid);

  off = iblock-PTR_SIZE;
  for (j = 0; j < novl; j++)
    { perm[j] = (Overlap *) off;
      off += EXO_SIZE + ((Overlap *) off)->path.tlen;
    }

  qsort(perm,novl,sizeof(Overlap *),SORT_MAP);

  tmaxl = 0;
  for (j = 0; j < novl; j++)
    { Overlap *o = perm[j];
      if (fwrite( ((void *) o)+PTR_SIZE, EXO_SIZE, 1, fid) != 1)
        { fprintf(stderr,"\n%s: Cannot not write sorted overlap block file %s/%s.%d.las\n",
                         Prog_Name,SORT_PATH,ALGN_UNIQ,parm->tid); exit (1); }
      if (fwrite( (void *) (o+1), o->path.tlen, 1, fid) != 1)
        { fprintf(stderr,"\n%s: Cannot not write sorted overlap block file %s/%s.%d.las\n",
                         Prog_Name,SORT_PATH,ALGN_UNIQ,parm->tid); exit (1); }
      if (o->path.tlen > tmaxl) tmaxl = o->path.tlen;
    }

  rewind(fid);
  free(perm);
  free(iblock-PTR_SIZE);
  parm->tmaxl = tmaxl;
  return (NULL);
}

// ============================================================================
// ovl_reload and la_merge
// ============================================================================

static void ovl_reload(IO_block *in, int64 bsize)
{ int64 remains = in->top - in->ptr;
  if (remains > 0) memmove(in->block, in->ptr, remains);
  in->ptr = in->block;
  in->top = in->block + remains;
  in->top += fread(in->top,1,bsize-remains,in->stream);
}

static int la_merge(TP *parm)
{ IO_block *in;
  int64     bsize;
  char     *block;
  int       i, c;
  Overlap **heap;
  int       hsize, tmaxl;
  Overlap  *ovls;
  int64     totl;
  OneFile  *of;
  int64    *trace64;

  bsize  = (MEMORY*1000000ll)/NTHREADS;
  block  = (char *) Malloc(bsize*NTHREADS+PTR_SIZE,"Allocating LAmerge blocks");
  in     = (IO_block *) Malloc(sizeof(IO_block)*NTHREADS,"Allocating LAmerge IO-reacords");
  if (block == NULL || in == NULL) return (1);
  block += PTR_SIZE;

  tmaxl = 0;
  totl = 0;
  for (c = 0; c < NTHREADS; c++)
    { void *iblock;
      in[c].stream = parm[c].ofile;
      in[c].block  = iblock = block+c*bsize;
      in[c].ptr    = iblock;
      in[c].top    = iblock + fread(iblock,1,bsize,parm[c].ofile);
      in[c].count  = 0;
      totl += parm[c].nlive;
      if (parm[c].tmaxl > tmaxl) tmaxl = parm[c].tmaxl;
    }

  heap = (Overlap **) Malloc(sizeof(Overlap *)*(NTHREADS+1),"Allocating heap");
  ovls = (Overlap *) Malloc(sizeof(Overlap)*NTHREADS,"Allocating heap");
  trace64 = (int64 *) Malloc(sizeof(int64)*(tmaxl/2),"Allocating int64 trace vector");
  if (heap == NULL || ovls == NULL) return (1);

  hsize = 0;
  for (i = 0; i < NTHREADS; i++)
    { if (in[i].ptr < in[i].top)
        { ovls[i]     = *((Overlap *) (in[i].ptr - PTR_SIZE));
          in[i].ptr  += EXO_SIZE;
          hsize      += 1;
          heap[hsize] = ovls + i;
        }
    }

  if (hsize > 3)
    for (i = hsize/2; i > 1; i--)
      maheap(i,heap,hsize);

  { char *db_name;
    char *cpath;
    db_name = Strdup(Catenate(SORT_PATH,"/",ALGN_UNIQ,".1aln"),"db_name");
    cpath = getcwd(NULL,0);
    of = open_Aln_Write(Catenate(SORT_PATH,"/",ALGN_UNIQ,".1aln"), 1,
                        "pga-merge", "0.1", "pga-merge", TSPACE,
                        db_name, NULL, cpath);
    Write_Skeleton(of,&parm->gdb1);
    free(cpath);
    free(db_name);
  }

  while (hsize > 0)
    { Overlap  *ov;
      IO_block *src;
      int64     tsize, span;

      maheap(1,heap,hsize);
      ov  = heap[1];
      src = in + (ov - ovls);
      src->count += 1;
      tsize = ov->path.tlen;
      span  = EXO_SIZE + tsize;
      if (src->ptr + span > src->top) ovl_reload(src,bsize);
      Write_Aln_Overlap(of, ov);
      Write_Aln_Trace(of, src->ptr, tsize, trace64, 0);
      src->ptr += tsize;
      if (src->ptr >= src->top)
        { heap[1] = heap[hsize];
          hsize  -= 1;
          continue;
        }
      *ov       = *((Overlap *) (src->ptr - PTR_SIZE));
      src->ptr += EXO_SIZE;
    }

  oneFileClose(of);
  for (i = 0; i < NTHREADS; i++) fclose(parm[i].ofile);
  for (i = 0; i < NTHREADS; i++) totl -= in[i].count;
  if (totl != 0)
    { fprintf(stderr,"%s: Did not write all records to %s/%s.1aln (%lld)\n",
                     Prog_Name,SORT_PATH,ALGN_UNIQ,totl); return (1); }

  free(trace64);
  free(ovls);
  free(heap);
  free(in);
  free(block-PTR_SIZE);
  return (0);
}

// ============================================================================
// pair_sort_search
// ============================================================================

void pair_sort_search(GDB *gdb1, GDB *gdb2)
{ uint8 *sarray;
  int    swide;
  int64  nels;
  RP     rarm[NTHREADS];
  TP     tarm[NTHREADS];
  pthread_t threads[NTHREADS];
  int64    *panel;
  Range     range[NTHREADS];
  IOBuffer *unit[2], *nu;
  int       nused;
  int       i, p, j, u;

  if (VERBOSE)
    { fprintf(stderr,"\n  Starting seed sort and alignment search, %d parts\n",2*NPARTS);
      fflush(stderr);
    }

  unit[0] = N_Units;
  unit[1] = C_Units;

  { int64 cum, nelmax;
    nelmax = 0;
    for (u = 0; u < 2; u++)
      { cum = 0;
        nu = unit[u];
        for (j = 0; j < NCONTS; j++)
          { for (i = 0; i < NTHREADS; i++)
              { cum += nu[i].buck[j];
                nu[i].buck[j] = cum;
              }
            if (j+1 == NCONTS || Select[j] != Select[j+1])
              { if (cum > nelmax) nelmax = cum;
                cum = 0;
              }
          }
        for (j = NCONTS-1; j >= 0; j--)
          { for (i = NTHREADS-1; i >= 1; i--)
              nu[i].buck[j] = nu[i-1].buck[j];
            if (j == 0 || Select[j] != Select[j-1])
              nu[0].buck[j] = 0;
            else
              nu[0].buck[j] = nu[NTHREADS-1].buck[j-1];
          }
      }

    swide  = 2*DBYTE + JCONT + 2;
    sarray = Malloc((nelmax+1)*swide,"Sort Array");
    panel  = Malloc(NCONTS*sizeof(int64),"Bucket Array");
    if (sarray == NULL || panel == NULL) exit (1);
  }

  for (p = 0; p < NTHREADS; p++)
    { rarm[p].swide  = swide;
      rarm[p].sarr   = sarray;
      rarm[p].buffer = N_Units[p].bufr;
      rarm[p].range  = range+p;
      rarm[p].gdb1   = gdb1;
      rarm[p].gdb2   = gdb2;

      tarm[p].tid    = p;
      tarm[p].swide  = swide;
      tarm[p].sarr   = sarray;
      tarm[p].panel  = panel;
      tarm[p].range  = range+p;
      tarm[p].gdb1   = *gdb1;
      tarm[p].gdb2   = *gdb2;
      if (p > 0)
        { if (gdb1->seqstate == EXTERNAL)
            { tarm[p].gdb1.seqs = fopen(gdb1->seqpath,"r");
              if (tarm[p].gdb1.seqs == NULL)
                { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name); exit (1); }
            }
          if (gdb2->seqstate == EXTERNAL)
            { tarm[p].gdb2.seqs = fopen(gdb2->seqpath,"r");
              if (tarm[p].gdb2.seqs == NULL)
                { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name); exit (1); }
            }
        }
      tarm[p].nhits = 0;
      tarm[p].nlass = 0;
      tarm[p].nlive = 0;
      tarm[p].nlcov = 0;
      tarm[p].nmemo = 0;
      tarm[p].ofile = fopen(Catenate(SORT_PATH,"/",ALGN_UNIQ,Numbered_Suffix(".",p,".las")),"w+");
      if (tarm[p].ofile == NULL)
        { fprintf(stderr,"%s: Cannot open %s/%s.%d.las for writing\n",
                         Prog_Name,SORT_PATH,ALGN_UNIQ,p); exit (1); }
      unlink(Catenate(SORT_PATH,"/",ALGN_UNIQ,Numbered_Suffix(".",p,".las")));
      tarm[p].tfile = fopen(Catenate(SORT_PATH,"/",ALGN_PAIR,Numbered_Suffix(".",p,".las")),"w+");
      if (tarm[p].tfile == NULL)
        { fprintf(stderr,"%s: Cannot open %s/%s.%d.las for reading & writing\n",
                         Prog_Name,SORT_PATH,ALGN_PAIR,p); exit (1); }
      unlink(Catenate(SORT_PATH,"/",ALGN_PAIR,Numbered_Suffix(".",p,".las")));
    }

  for (u = 0; u < 2; u++)
   for (i = 0; i < NPARTS; i++)
    { nu = unit[u] + i*NTHREADS;
      if (VERBOSE)
        { fprintf(stderr,"\r    Loading seeds for part %d  ",u*NPARTS+i+1);
          fflush(stderr);
        }
      for (p = 0; p < NTHREADS; p++)
        { rarm[p].in = nu[p].file;
          lseek(nu[p].file,0,SEEK_SET);
          rarm[p].buck = nu[p].buck;
          rarm[p].comp = u;
          rarm[p].inum = nu[p].inum;
        }
      for (p = 1; p < NTHREADS; p++)
        pthread_create(threads+p,NULL,reimport_thread,rarm+p);
      reimport_thread(rarm);
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);

      { int64 prev, next;
        bzero(panel,sizeof(int64)*NCONTS);
        prev = 0; next = 0;
        for (j = IDBsplit[i]; j < IDBsplit[i+1]; j++)
          { next = nu[NTHREADS-1].buck[j];
            panel[j] = (next - prev)*swide;
            prev = next;
          }
        nels = next;
        if (VERBOSE)
          { fprintf(stderr,"\r    Sorting seeds for part %d  ",u*NPARTS+i+1);
            fflush(stderr);
          }
        nused = rmsd_sort(sarray,nels,swide,swide,NCONTS,panel,NTHREADS,range);
      }
      if (VERBOSE)
        { fprintf(stderr,"\r    Searching seeds for part %d",u*NPARTS+i+1);
          fflush(stderr);
        }
      for (p = 0; p < nused; p++) tarm[p].comp = u;
      for (p = 1; p < nused; p++)
        pthread_create(threads+p,NULL,search_seeds,tarm+p);
      search_seeds(tarm);
      for (p = 1; p < nused; p++)
        pthread_join(threads[p],NULL);
    }

  free(panel);
  free(sarray);
  for (p = 0; p < NTHREADS; p++) fclose(tarm[p].tfile);
  for (p = 1; p < NTHREADS; p++)
    { if (gdb2->seqstate == EXTERNAL) fclose(tarm[p].gdb2.seqs);
      if (gdb1->seqstate == EXTERNAL) fclose(tarm[p].gdb1.seqs);
    }

  if (VERBOSE)
    { int64 nhit, nlas, nliv, ncov;
      fprintf(stderr,"\r    Done                        \n");
      nhit = nlas = nliv = ncov = 0;
      for (p = 0; p < NTHREADS; p++)
        { nhit += tarm[p].nhits;
          nlas += tarm[p].nlass;
          nliv += tarm[p].nlive;
          ncov += tarm[p].nlcov;
        }
      if (nliv == 0)
        fprintf(stderr,"\n  Total hits over %dbp = %lld, %lld aln's, 0 %s\n",
                       CHAIN_MIN/2,nhit,nlas,"non-redundant aln's of ave len 0");
      else
        fprintf(stderr,"\n  Total hits over %dbp = %lld, %lld aln's, %lld %s %lld\n",
                       CHAIN_MIN/2,nhit,nlas,nliv,"non-redundant aln's of ave len",ncov/nliv);
      fflush(stderr);
    }

  if (VERBOSE)
    { fprintf(stderr,"\n  Sorting and merging alignments\n");
      fflush(stderr);
    }

  for (p = 1; p < NTHREADS; p++)
    pthread_create(threads+p,NULL,la_sort,tarm+p);
  la_sort(tarm);
  for (p = 1; p < NTHREADS; p++)
    pthread_join(threads[p],NULL);

  if (la_merge(tarm)) exit (1);
}