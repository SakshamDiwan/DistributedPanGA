// fastga_order_core.c — MPI-free rules behind fastga_order.c; see the header.

#include "fastga_order_core.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int fo_read_gix_perm(const char *stem, int n, int *perm, char *err, int errlen)
{ char    path[4096];
  int     f, pbyte, cbyte, nfile, freq, nctg, i, ok;
  int64_t maxp;
  char   *seen;

  snprintf(path,sizeof(path),"%s.gix",stem);
  f = open(path,O_RDONLY);
  if (f < 0)
    { snprintf(err,errlen,"cannot open %s",path); return 0; }
  ok = (lseek(f,4*sizeof(int)+0x1000000*sizeof(int64_t),SEEK_SET) >= 0);
  ok = ok && read(f,&pbyte,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&cbyte,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&nfile,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&maxp,sizeof(int64_t)) == sizeof(int64_t);
  ok = ok && read(f,&freq,sizeof(int)) == sizeof(int);
  ok = ok && read(f,&nctg,sizeof(int)) == sizeof(int);
  if (!ok || nctg != n)
    { snprintf(err,errlen,"%s: unexpected header (contigs %d, GDB has %d)",path,ok ? nctg : -1,n);
      close(f);
      return 0;
    }
  ok = (read(f,perm,sizeof(int)*(size_t) n) == (ssize_t) (sizeof(int)*(size_t) n));
  close(f);
  if (!ok)
    { snprintf(err,errlen,"%s: short read of the contig order",path); return 0; }
  seen = calloc(n > 0 ? n : 1,1);
  for (i = 0; i < n; i++)
    { if (perm[i] < 0 || perm[i] >= n || seen[perm[i]])
        { snprintf(err,errlen,"%s: contig order is not a permutation",path);
          free(seen);
          return 0;
        }
      seen[perm[i]] = 1;
    }
  free(seen);
  return 1;
}

int fo_compute_ranges(int n, const int64_t *clen, int64_t npost, const int64_t *seeds,
                      int64_t swide, int T, int *range_start, int *beg_out, int *end_out,
                      int *nparts_out)
{ int     *split = malloc((T+1)*sizeof(int));
  int64_t  cum, t;
  int      p, r, x, i, nparts, nrange;

  split[0] = 0;
  p = 0;
  r = T;
  t = npost/T;
  cum = (n > 0) ? clen[0] : 0;
  for (x = 1; x < n; x++)
    { if (cum >= t && x >= r)
        { p += 1;
          split[p] = x;
          t = (npost*(p+1))/T;
          r += T;
        }
      cum += clen[x];
    }
  nparts = p+1;
  split[nparts] = n;

  for (x = 0; x < n; x++)
    range_start[x] = x;

  nrange = 0;
  for (i = 0; i < nparts; i++)
    { int64_t asize, thr, sum;
      int     k, beg, y;

      asize = 0;
      for (x = split[i]; x < split[i+1]; x++)
        asize += seeds[x]*swide;
      if (asize == 0) continue;

      k   = 0;
      thr = asize/T;
      sum = 0;
      for (x = split[i]; x < split[i+1]; x++)
        if (seeds[x] > 0)
          break;
      beg = x;
      for ( ; x < split[i+1]; x++)
        if (seeds[x] > 0)
          { sum += seeds[x]*swide;
            if (sum >= thr)
              { for (y = beg; y <= x; y++)
                  range_start[y] = beg;
                if (beg_out != NULL)
                  { beg_out[nrange] = beg;
                    end_out[nrange] = x+1;
                  }
                nrange += 1;
                k  += 1;
                thr = (asize*(k+1))/T;
                beg = x+1;
              }
          }
    }
  free(split);
  *nparts_out = nparts;
  return nrange;
}

int fo_nearness(const uint8_t *has, int n, int xa, int xb)
{ const uint8_t *row = has + (int64_t) xa*n;
  int y, k = 0;

  for (y = xb+1; y < xa; y++)
    k += (row[y] != 0);
  return k;
}

int fo_count_earlier(const uint8_t *has, int n, int xa)
{ const uint8_t *row = has + (int64_t) xa*n;
  int y, k = 0;

  for (y = 0; y < xa; y++)
    k += (row[y] != 0);
  return k;
}

int fo_kth_earlier(const uint8_t *has, int n, int xa, int k)
{ const uint8_t *row = has + (int64_t) xa*n;
  int y;

  for (y = xa-1; y >= 0; y--)
    if (row[y] && k-- == 0)
      return y;
  return -1;
}

int fo_unset_value(int range_start_x, int unset_main)
{ return (range_start_x == 0) ? unset_main : 0; }

int fo_from_answers(const int *answer, int nnear, int x, int range_start_x, int unset_main,
                    int *bepos, int *k_used)
{ int k;

  for (k = 0; k < nnear; k++)
    { if (answer[k] == FO_ANSWER_UNKNOWN)
        return 0;
      if (answer[k] >= 0)
        { *bepos  = answer[k];
          *k_used = k;
          return 1;
        }
    }
  if (range_start_x != x)
    return -1;
  *bepos  = fo_unset_value(range_start_x,unset_main);
  *k_used = -1;
  return 1;
}

int fo_from_rows(const int64_t *L, const int64_t *S, int x, int range_start_x, int unset_main,
                 int *bepos, int *row, int *from_pos)
{ int y;

  for (y = x-1; y >= range_start_x; y--)
    { int64_t k = L[y];
      int64_t s = S[y];

      if (k != 0 && fo_key_pos(k) > y)
        { *bepos = fo_key_bepos(k); *from_pos = fo_key_pos(k); }
      else if (s == FO_S_PENDING)
        return 0;
      else if (s >= FO_S_CALLED)
        { *bepos = (int) (s-FO_S_CALLED); *from_pos = y; }
      else if (k != 0)
        { *bepos = fo_key_bepos(k); *from_pos = fo_key_pos(k); }
      else
        continue;
      *row = y;
      return 1;
    }
  *bepos    = fo_unset_value(range_start_x,unset_main);
  *row      = -1;
  *from_pos = -1;
  return 1;
}
