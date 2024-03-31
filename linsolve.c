/* linsolve.c - solve a set of linear equations */
/*
 * Solves the matrix equation   Y = A X,   returns X.
 *        where y[j=0,n-1], x[i=0,n-1], a[i+j*n].
 * 
 * linsolven will solve NY versions of the equations at once; same A
 * but different Y.  In this case x,y = x,y[ny][0:n-1]
 *
 */
/* 021007 - Rev 1.1 add NY option */
/* 020211 - Rev 1.0 John Tonry */

#include <stdio.h>
#include <math.h>

/* Invert a matrix, return the determinant */
int invert(int n, double *a, double *det);

/* Solves the matrix equation   Y = A X,   returns X */
int linsolve(int n, double *y, double *a, double *x);

/* Solves the matrix equation   Y = A X,   returns X for ny versions */
int linsolven(int n, int ny, double *y, double *a, double *x);

#define ABS(x) ((x)<0?(-(x)):(x))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MAXEQ 2000

int linsolven(int n, int ny, double *y, double *a, double *x)
{
   int i, j, k, iy;
   int rowstatus[MAXEQ], row[MAXEQ];
   double rat;

   if(n > MAXEQ) {
      fprintf(stderr, "WHOA: too many equations %d req %d max\n", n, MAXEQ);
      return(-1);
   }

   for(i=0; i<n; i++) rowstatus[i] = 0;

/* Solve matrix by Gaussian elimination */
   for(j=0; j<n; j++) {		/* j = column to be zero'ed out */

/* Find a good equation to work on (pivot row) */
      for(i=k=0, rat=0.0; i<n; i++) {
	 if(rowstatus[i]) continue;
	 if(ABS(a[j+i*n]) > rat) {
	    k = i;
	    rat = ABS(a[j+k*n]);
	 }
      }
      rowstatus[k] = 1;
      row[j] = k;
      if(rat == 0.0) {
//	 fprintf(stderr, "WHOA: singular matrix, iter %d, row %d\n", j, k);
	 return(1);
      }
/* Subtract away matrix below diagonal */
      for(i=0; i<n; i++) {	/* i = subsequent equations */
	 if(rowstatus[i]) continue;
	 rat = a[j+i*n] / a[j+row[j]*n];
	 for(iy=0; iy<ny; iy++) y[i+iy*n] -= rat * y[row[j]+iy*n];
 	 for(k=j+1; k<n; k++) a[k+i*n] -= rat * a[k+row[j]*n];
      }
   }

/* Back substitute into upper diagonal matrix to solve for parameters */
   for(j=n-1; j>=0; j--) {
      for(iy=0; iy<ny; iy++) {
	 x[j+iy*n] = y[row[j]+iy*n];
	 for(k=j+1; k<n; k++) x[j+iy*n] -= a[k+row[j]*n] * x[k+iy*n];
	 x[j+iy*n] /= a[j+row[j]*n];
      }
   }
   return(0);
}

int linsolve(int n, double *y, double *a, double *x)
{
   return(linsolven(n, 1, y, a, x));
}

#define MAXDEC 1e10	/* Maximum decrement in pivot value for singular */

/* Invert a matrix, return the determinant */
/* 770902 John Tonry */
int invert(int n, double *a, double *det)
{
   double save, pivot, onrow, cprev, cnow, decr;
   short int rst[MAXEQ][2];
   int i, j, k, l, mrank, isign, nrow=0, ncol=0;

   if(n > MAXEQ) {
      fprintf(stderr, "Matrix size too big: %d > %d\n", n, MAXEQ);
      return(-2);
   }
   mrank = 0;
   isign = 1;
   *det = 0.0;
   for(j=0; j<n; j++) rst[j][0] = rst[j][1] = -1;
/* Loop over columns, reducing each */
   for(i=0; i<n; i++) {

/* Find the pivot element */
      pivot = -1.0;
      for(j=0; j<n; j++) {
	 if(rst[j][0] != -1) continue;
	 for(k=0; k<n; k++) {
	    if(rst[k][0] != -1) continue;
	    if(pivot >= ABS(a[j+k*n])) continue;
	    pivot = ABS(a[j+k*n]);
	    nrow = j;
	    ncol = k;
	 }
      }
      pivot = a[nrow+ncol*n];
      if(pivot == 0.0) {
	 *det = 0;
	 return(-1);
      }
      rst[ncol][0] = nrow;
      rst[ncol][1] = i;
/* Swap pivot element onto the diagonal */
      for(k=0; k<n; k++) {
	 save = a[nrow+k*n];
	 a[nrow+k*n] = a[ncol+k*n];
	 a[ncol+k*n] = save;
      }
/*   Reduce pivot column */
      for(j=0; j<n; j++) a[j+ncol*n] = -a[j+ncol*n]/pivot;
      a[ncol+ncol*n] = 1/pivot;

/*   Reduce other columns */
      for(k=0; k<n; k++) {
	 if(k == ncol) continue;
/*     Find maximum of column to check for singularity */
	 cprev = 0;
	 for(j=0; j<n; j++) cprev = MAX(cprev,ABS(a[j+k*n]));
/*     Reduce the column */
	 onrow = a[ncol+k*n];
	 a[ncol+k*n] = 0;
	 for(j=0; j<n; j++) a[j+k*n] = a[j+k*n] + onrow*a[j+ncol*n];
/*     Find the new maximum of the column */
	 cnow = 0;
	 for(j=0; j<n; j++) cnow = MAX(cnow, ABS(a[j+k*n]));

/*     Quit if too many figures accuracy were lost (singular) */
	 decr = cprev / cnow;
	 if(cnow == 0.0 || decr > MAXDEC) {
	    *det = 0;
	    return(-1);
	 }
      }
      *det = *det + log(ABS(pivot));
      if(pivot < 0) isign *= -1;
      mrank++;
   }

/*     Now untangle the mess */
   for(j=0; j<n; j++) {
      for(k=0; k<n; k++) {
	 if(rst[k][1] != (n-1-j)) continue;
	 ncol = rst[k][0];
	 if(ncol == k) break;
	 for(l=0; l<n; l++) {
	    save = a[l+ncol*n];
	    a[l+ncol*n] = a[l+k*n];
	    a[l+k*n] = save;
	 }
	 break;
      }
   }
	 
   if(ABS(*det) < 88) *det = isign * exp(*det);
   return(0);
}
