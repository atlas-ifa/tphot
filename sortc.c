/* A few sorting routines */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Timings on the Dell D830 laptop: */
/*    qsort8:  t ~  75[usec] * (N log_2(N) / 10000) */
/*    qsort:   t ~ 530[usec] * (N log_2(N) / 10000) */
/*    median:  t ~  16[usec] * (N/1000) + 510[usec] * (N/1000)^2 */
/*    hmedian:  10% slower than qsort8 */
/*    median_sift:  x2 slower than qsort8 for N<50;  
 *                  x2 faster than qsort8 for N>500 (2 pass at N=10) */
/*    quickmedian ~ 1.2x faster than median_sift for N > 200 */

typedef double DATA;


//#define TEST		/* A main to do timing testing */
#ifdef TEST
#include <sys/time.h>
#endif

#define MAXSTACK 256
//#define NSTOP 15
#define NSTOP 8

/* Quicksort program for floats: 1980 (John Tonry) */
int qsort4(int n, float *x)
{
  float key, kl, kr, km, temp;
  int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], sp;
  register int i, j, k;
  int mgtl, lgtr, rgtm;

  sp = 0;
  lstack[sp] = 0;
  rstack[sp] = n-1;

  while(n > NSTOP && sp >= 0) {
/* Sort a subrecord off the stack */
    l = lstack[sp];
    r = rstack[sp];
    sp--;
    m = (l + r) / 2;
/* Set KEY = median of X(L), X(M), X(R) */
    kl = x[l];
    km = x[m];
    kr = x[r];
    mgtl = km > kl;
    rgtm = kr > km;
    lgtr = kl > kr;
#ifndef STD_KEY
    if(mgtl ^ rgtm) {
      if(mgtl ^ lgtr) key = kr;
      else            key = kl;
    } else {
      key = km;
    }
#else
/* Curiously enough, this non-median seems to work as well or better */
    if(mgtl ^ rgtm) {
      key = km;
    } else {
      if(mgtl ^ lgtr) key = kl;
      else            key = kr;
    }
#endif
    i = l;
    j = r;
    while(1) {
/* Find a big record on the left */
      while(x[i] < key) i++;

/* Find a small record on the right */
      while(x[j] > key) j--;

      if(i >= j) break;
/* Exchange records */
      temp = x[i];
      x[i] = x[j];
      x[j] = temp;
      i++;
      j--;
    }

/* Subfile is partitioned into two halves, left .le. right */
/* Push the two halves on the stack */
    if(j-l+1 > NSTOP) {
      lstack[++sp] = l;
      rstack[sp] = j;
    }
    if(sp >= MAXSTACK) break;
    if(r-j > NSTOP) {
      lstack[++sp] = j+1;
      rstack[sp] = r;
    }
    if(sp >= MAXSTACK) break;
  }

/* Insertion sorting routine to finish up */
  for(j=n-2; j>=0; j--) {
    for(i=j+1, k=j; i<n; i++) {
      if(x[j] <= x[i]) break;
      k = i;
    }
    if(k != j) {
      temp = x[j];
      for(i=j+1; i<=k; i++) x[i-1] = x[i];
      x[k] = temp;
    }
  }
  return(0);
}

/* Quicksort program for doubles: 1980 (John Tonry) */
int qsort8(int n, DATA *x)
{
  DATA key, kl, kr, km, temp;
  int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], sp;
  register int i, j, k;
  int mgtl, lgtr, rgtm;

  sp = 0;
  lstack[sp] = 0;
  rstack[sp] = n-1;

  while(n > NSTOP && sp >= 0) {
/* Sort a subrecord off the stack */
    l = lstack[sp];
    r = rstack[sp--];
    m = (l + r) / 2;
/* Set KEY = median of X(L), X(M), X(R) */
    kl = x[l];
    km = x[m];
    kr = x[r];
    mgtl = km > kl;
    rgtm = kr > km;
    lgtr = kl > kr;
#ifndef STD_KEY
    if(mgtl ^ rgtm) {
      if(mgtl ^ lgtr) key = kr;
      else            key = kl;
    } else {
      key = km;
    }
#else
/* Curiously enough, this non-median seems to work as well or better */
    if(mgtl ^ rgtm) {
      key = km;
    } else {
      if(mgtl ^ lgtr) key = kl;
      else            key = kr;
    }
#endif
    i = l;
    j = r;
    while(1) {
/* Find a big record on the left */
      while(x[i] < key) i++;

/* Find a small record on the right */
      while(x[j] > key) j--;

      if(i >= j) break;
/* Exchange records */
      temp = x[i];
      x[i++] = x[j];
      x[j--] = temp;
    }

/* Subfile is partitioned into two halves, left .le. right */
/* Push the two halves on the stack */
    if(j-l+1 > NSTOP) {
      lstack[++sp] = l;
      rstack[sp] = j;
    }
    if(sp >= MAXSTACK) break;
    if(r-j > NSTOP) {
      lstack[++sp] = j+1;
      rstack[sp] = r;
    }
    if(sp >= MAXSTACK) break;
  }

/* Insertion sorting routine to finish up */
  for(j=n-2; j>=0; j--) {
    for(i=j+1, k=j; i<n; i++) {
      if(x[j] <= x[i]) break;
      k = i;
    }
    if(k != j) {
      temp = x[j];
      for(i=j+1; i<=k; i++) x[i-1] = x[i];
      x[k] = temp;
    }
  }
  return(0);
}

/* Quicksort program for doubles, sort index: 1980 (John Tonry) */
int qsort2(int n, DATA *x, int *idx)
{
  DATA key, kl, kr, km, temp;
  int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], sp, itemp;
  register int i, j, k;
  int mgtl, lgtr, rgtm;

  sp = 0;
  lstack[sp] = 0;
  rstack[sp] = n-1;

  while(n > NSTOP && sp >= 0) {
/* Sort a subrecord off the stack */
    l = lstack[sp];
    r = rstack[sp--];
    m = (l + r) / 2;
/* Set KEY = median of X(L), X(M), X(R) */
    kl = x[l];
    km = x[m];
    kr = x[r];
    mgtl = km > kl;
    rgtm = kr > km;
    lgtr = kl > kr;
#ifndef STD_KEY
    if(mgtl ^ rgtm) {
      if(mgtl ^ lgtr) key = kr;
      else            key = kl;
    } else {
      key = km;
    }
#else
/* Curiously enough, this non-median seems to work as well or better */
    if(mgtl ^ rgtm) {
      key = km;
    } else {
      if(mgtl ^ lgtr) key = kl;
      else            key = kr;
    }
#endif
    i = l;
    j = r;
    while(1) {
/* Find a big record on the left */
      while(x[i] < key) i++;

/* Find a small record on the right */
      while(x[j] > key) j--;

      if(i >= j) break;
/* Exchange records */
      temp = x[i];
      x[i] = x[j];
      x[j] = temp;
      itemp = idx[i];
      idx[i] = idx[j];
      idx[j] = itemp;
      i++;
      j--;
    }

/* Subfile is partitioned into two halves, left .le. right */
/* Push the two halves on the stack */
    if(j-l+1 > NSTOP) {
      lstack[++sp] = l;
      rstack[sp] = j;
    }
    if(sp >= MAXSTACK) break;
    if(r-j > NSTOP) {
      lstack[++sp] = j+1;
      rstack[sp] = r;
    }
    if(sp >= MAXSTACK) break;
  }

/* Insertion sorting routine to finish up */
  for(j=n-2; j>=0; j--) {
    for(i=j+1, k=j; i<n; i++) {
      if(x[j] <= x[i]) break;
      k = i;
    }
    if(k != j) {
      temp = x[j];
      itemp = idx[j];
      for(i=j+1; i<=k; i++) {
	 x[i-1] = x[i];
	 idx[i-1] = idx[i];
      }
      x[k] = temp;
      idx[k] = itemp;
    }
  }
  return(0);
}

/* median.c : simple bubble sort, carry another array, return median
 *
 * 020130 John Tonry
 */
DATA median(int n, DATA *key, int *idx)
{
   register int i, j, k;
   int itmp=0;
   DATA tmp;

   if(n == 0) return(0.0);
   if(n == 1) return(key[0]);
   if(idx == NULL) {
      for(j=n-2; j>=0; j--) {
	 for(i=j+1, k=j; i<n; i++) {
	    if(key[j] <= key[i]) break;
	    k = i;
	 }
	 if(k != j) {
	    tmp = key[j];
	    for(i=j+1; i<=k; i++) key[i-1] = key[i];
	    key[k] = tmp;
	 }
      }
   } else {
      for(j=n-2; j>=0; j--) {
	 for(i=j+1, k=j; i<n && key[j] > key[i]; i++) k = i;
	 if(k != j) {
	    tmp = key[j];
	    itmp = idx[j];
	    for(i=j+1; i<=k; i++) { key[i-1] = key[i]; idx[i-1] = idx[i]; }
	    key[k] = tmp;
	    idx[k] = itmp;
	 }
      }
   }
   tmp = 0.5 * (key[n/2]+key[(n-1)/2]);		/* Median */
   return(tmp);
}


/* median_sift() moves central points which straddle goal to start of array */
/* Typical usage:
 *
 *    median_sift(n, buf, n/2, 10, &nnew, &nbelow);
 *    median(nnew, buf, NULL);
 *    med = 0.5*(buf[(n-1)/2-nbelow] + buf[n/2-nbelow]);
 *
 * Note that median_sift() can be called multiple times, e.g.
 *
 *    median_sift(n, buf, n/2, 10, &nnew, &nbelow);
 *    median_sift(nnew, buf,  n/2-nbelow, 10, &nnew, &nbelow2);
 *    median(nnew, buf, NULL);
 * 
 * The median will be found at goal-nbelow (and prev for even), e.g.
 *
 *    med = 0.5*(buf[(n-1)/2-nbelow-nbelow2] + buf[n/2-nbelow-nbelow2]);
 */

/* 090315 John Tonry */
static DATA *sift_posts=NULL;
static int *sift_ngelt=NULL;
static int nsift=0;
int median_sift(int n, DATA *key, int goal, int nbin,
		int *nnew, int *nbelow)
{
   int i, j, k;
   DATA tmp;

   if(n < nbin+1) {
      *nnew = n;
      *nbelow = 0;
      return(0);
   }

/* Allocate some buffer space */
   if(nbin > nsift) {
      if(sift_posts != NULL) free(sift_posts);
      if(sift_ngelt != NULL) free(sift_ngelt);
      sift_posts = (DATA *)calloc(nbin, sizeof(DATA));
      sift_ngelt = (int *)calloc(nbin, sizeof(int));
      nsift = nbin;
   }

/* Sample some points */
   for(i=1; i<nbin; i++) sift_posts[i] = key[(i*(n-1))/nbin];
/* and sort them */
   median(nbin-1, sift_posts+1, NULL);

/* Here's what we're calculating
 *    sift_posts[0]           = minimum value
 *    sift_posts[1:nbin-1] = suggested division points amongst the values
 *    sift_posts[nbin]     = maximum value
 *    
 *    sift_ngelt[0]           = number below sift_posts[1]
 *    sift_ngelt[1:nbin-1] = number >= sift_posts[j] and < sift_posts[j+1]
 *      (i.e. equal to sift_posts[j] goes into sift_ngelt[j])
 */

   sift_posts[0] = sift_posts[1];
   sift_posts[nbin] = sift_posts[nbin-1];
   for(i=0; i<=nbin; i++) sift_ngelt[i] = 0;
/* Count up all the points and find the min and max */
   for(i=0; i<n; i++) {
      if(key[i] < sift_posts[1]) {
	 if(key[i] < sift_posts[0]) sift_posts[0] = key[i];
	 sift_ngelt[0] += 1;
      } else if(key[i] >= sift_posts[nbin-1]) {
	 if(key[i] >= sift_posts[nbin]) sift_posts[nbin] = key[i];
	 sift_ngelt[nbin-1] += 1;
      } else {
	 for(j=1; j<nbin-1; j++) {
	    if(key[i] >= sift_posts[j] && key[i] < sift_posts[j+1]) {
	       sift_ngelt[j] += 1;
	       break;
	    }
	 }
      }
   }

/* Figure out the points which can be ignored for the median */
/* The points between sift_posts[j] and sift_posts[k] contain the goal point */
/* There are nbelow and nabove=i points which don't need to be examined */
   for(j=0, *nbelow=0; j<nbin; j++) {
      if(*nbelow + sift_ngelt[j] >= goal) break;
      *nbelow += sift_ngelt[j];
   }
   for(k=nbin, i=0; k>=1; k--) {
      if(i + sift_ngelt[k-1] >= n-goal-1) break;
      i += sift_ngelt[k-1];
   }
   

/* Swap the points between sift_posts[j] and sift_posts[k] to the beginning */
   for(i=0, *nnew=0; i<n; i++) {
      if(key[i] >= sift_posts[j] && key[i] < sift_posts[k]) {
	 tmp = key[i];
	 key[i] = key[*nnew];
	 key[*nnew] = tmp;
	 *nnew += 1;
      }
   }
#if 0
   if(test) {
      printf("\n%5d %5d %5d %5d %8.3f %8.3f %10d %10d %10d\n", 
	     n, goal, j, k, sift_posts[j], sift_posts[k], *nnew, *nbelow, 
	     n-*nnew-*nbelow);
      for(i=0; i<=nbin; i++) 
	 printf("%3d %10.6f %5d\n", i, sift_posts[i], sift_ngelt[i]);

   }
#endif

   return(0);
}

#ifdef TEST

DATA *qbuf1=NULL, *qbuf2=NULL;
int qnbuf=0;

/* Fast? median for small numbers.  Does not alter input array. */
/* Actually only appears to beat qsort8 for N > 50 */
/* Not thoroughly checked, not optimized at all! */
DATA quickmedian(int n, DATA *key)
{
   int i, j, m;
   int nlt, ngt, neq, ngoal1, ngoal2;
   DATA test;
   DATA *from, *to;

/* Quick return for trivial cases */
   if(n < 3) {
      if(n == 0) return(0.0);
      if(n == 1) return(key[0]);
      if(n == 2) return(0.5*(key[0]+key[1]));
   }

/* Allocate space for median accumulation */
   if(n > qnbuf) {
      if(qbuf1 != NULL) free(qbuf1);
      if(qbuf2 != NULL) free(qbuf2);
      qbuf1 = (DATA *)calloc(n, sizeof(DATA));
      qbuf2 = (DATA *)calloc(n, sizeof(DATA));
      qnbuf = n;
   }

   ngoal1 = (n-1)/2;
   ngoal2 = n/2;
   nlt = ngt = neq = 0;
   from = key;
   to = qbuf1;
/* Cast out values not straddling median */
   while(n > 2) {
      test = from[n/2];
      for(i=m=0; i<n; i++) if(from[i] < test) m++;
      if(nlt+m-1 >= ngoal1) {
	 for(i=m=0; i<n; i++) if(from[i] < test) to[m++] = from[i];
	 if(nlt+m-1 == ngoal1) to[m++] = test;
	 ngt += n-m;
      } else {
	 neq = m;
	 for(i=m=0; i<n; i++) {
	    if(from[i] > test) to[m++] = from[i];
	    if(from[i] == test && neq) {
	       to[m++] = from[i];
	       neq++;
	    }
	 }
	 nlt += n-m;
      }
      n = m;
      from = to;
      to = (from == qbuf1) ? qbuf2 : qbuf1;
   }
   if(n == 2 && from[0] > from[1]) {
      test = from[0];
      from[0] = from[1];
      from[1] = test;
   }
   if(ngoal1 == ngoal2) return(from[ngoal1-nlt]);
   else	                return(0.5*(from[ngoal1-nlt]+from[ngoal2-nlt]));
}

/* Median using Heapsort algorithm (based on Num.Rec, ex sextractor). */
DATA hmedian(DATA *ra, int n)
{
   int		l, j, ir, i;
   DATA	rra;
   if (n<2) return *ra;
   ra--;	/* Bleah, fake up fortran indexing */
   for (l = ((ir=n)>>1)+1;;) {
      if (l>1) {
	 rra = ra[--l];
      } else {
	 rra = ra[ir];
	 ra[ir] = ra[1];
	 if (--ir == 1) {
	    ra[1] = rra;
	    return n&1? ra[n/2+1] : 0.5*(ra[n/2]+ra[n/2+1]);
	 }
      }
      for (j = (i=l)<<1; j <= ir;) {
	 if (j < ir && ra[j] < ra[j+1]) ++j;
	 if (rra < ra[j]) {
	    ra[i] = ra[j];
	    j += (i=j);
	 } else {
	    j = ir + 1;
	 }
      }
      ra[i] = rra;
   }
}

//#define ONEPASS

int dblcmp(DATA *v1, DATA *v2)
{
   if(*v1 < *v2) return(-1);
   if(*v1 > *v2) return(1);
   return(0);
}

main(int argc, char **argv)
{
   int i, k, n, niter, nnew, nbelow, nnew2, nbelow2;
   double time, time0, dummy, med, medok, medok1, medok2;
   DATA *std, *buf;
   struct timeval tv0, tv1;
   DATA median(int n, DATA *key, int *idx);

   n = 1000;
   niter = 1000;
   if(argc > 1) sscanf(argv[1], "%d", &n);
   if(argc > 2) sscanf(argv[2], "%d", &niter);

   std = (DATA *)calloc(n*niter, sizeof(DATA));
   for(i=0; i<n*niter; i++) std[i] = rand() / ((DATA)RAND_MAX);
   buf = (DATA *)calloc(n, sizeof(DATA));

/* How much time just to loop and copy the buffer? */
   gettimeofday(&tv0, NULL);
   for(k=0, dummy=0.0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      dummy += k * buf[(k*n)/niter];
   }
   gettimeofday(&tv1, NULL);
   time0 = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

   printf("n= %5d", n);
/* Test qsort8 */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      qsort8(n, buf);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   medok = 0.5*(buf[(n-1)/2]+buf[n/2]);
   printf("  qsort8= %8.3f %8.5f", 1e6*(time-time0)/niter, medok);
   medok1 = buf[(n-1)/2];
   medok2 = buf[n/2];

/* Test qsort */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      qsort(buf, n, sizeof(DATA), dblcmp);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf[(n-1)/2]+buf[n/2]);
   printf("  qsort= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      median(n, buf, NULL);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf[(n-1)/2]+buf[n/2]);
   printf("  median= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test heap sort median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      hmedian(buf, n);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   med = 0.5*(buf[(n-1)/2]+buf[n/2]);
   printf("  hmedian= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test quick median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
      med = quickmedian(n, buf);
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
   printf("  quickmedian= %8.3f %8.5f", 1e6*(time-time0)/niter, med);

/* Test decimated median */
   gettimeofday(&tv0, NULL);
   for(k=0; k<niter; k++) {
      memcpy(buf, std+k*n, n*sizeof(DATA));
/* Note that median_sift can be chained... */
      median_sift(n, buf, n/2, 10, &nnew, &nbelow);
#ifdef ONEPASS
      median(nnew, buf, NULL);
#else
      median_sift(nnew, buf, n/2-nbelow, 10, &nnew2, &nbelow2);
//      if(k==0) for(i=0; i<nnew; i++) printf("%4d %8.5f\n", i, buf[i]);
      median(nnew2, buf, NULL);
#endif
   }
   gettimeofday(&tv1, NULL);
   time = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);
#ifdef ONEPASS
   med = 0.5*(buf[(n-1)/2-nbelow]+buf[n/2-nbelow]);
#else
   med = 0.5*(buf[(n-1)/2-nbelow-nbelow2]+buf[n/2-nbelow-nbelow2]);
#endif
   printf("  medeci= %8.3f %8.5f [usec] %5d   %.1e\n", 
	  1e6*(time-time0)/niter, med, nnew, dummy);


   if(med != medok) {
      fprintf(stderr, "Difference in median at n=%d %8.5f %8.5f %8.5f %8.5f %5d %5d %5d %5d %5d %5d\n", 
	      n, med, medok, medok1, medok2, nnew, nbelow, nnew2, nbelow2,
	      (n-1)/2-nbelow-nbelow2, n/2-nbelow-nbelow2);
      for(i=(n-1)/2-nbelow-nbelow2-2; i<n/2-nbelow-nbelow2+2; i++) {
	 fprintf(stderr, "%d %8.5f\n", i, buf[i]);
      }
   }

}

#if 0
// Some test code...
( for i in {0,1}{0,1,2,3,4,5,6,7,8,9} ; do sortc $i 20000 ; done ; for i in {1,2,3,5}{0,1,00,01} ; do sortc $i 5000 ; done ; for i in {1,2,3,5}000 ; do sortc $i 300 ; done ) > /tmp/foo
#endif

#endif
