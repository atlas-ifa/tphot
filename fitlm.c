/* Fit a model to data, minimizing chi^2 using Levenberg-Marquardt */
/* 121001: Reworked to fit a function to data instead of minimization
 * 020904: Transcribed to C
 * 11/17/82: Revision 2.0, John Tonry
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Define for a test main() */
// #define MAINTEST

#define MAXPAR  (128)		/* Max number of parameters (HUGE) */
#define LAMBASE (10.0)		/* Base for lambda */

//#define SHOWVAL		/* Print model and diff values for each iter */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

/* fitlm control parameters */
#define FITLM_VERBOSE 0x01	/* Verbose? */
#define FITLM_DERIV   0x02	/* Analytic first derivatives? */
#define FITLM_LINEAR  0x04	/* Linear function? */
#define FITLM_WGT     0x08	/* Data point weights provided? */

static int MAXITER=20;		/* Max iterations */
static double QUITFRAC=1e-8;	/* Fractional change for quitting */
static double PARVARY=1e-5;	/* Fractional parameter change for derivs */
static int LAMINIT=-2;		/* Initial value for lambda */


/* Definition of function that fitlm() fits to data */
/* Arguments:
 *
 * -Inputs-
 *   kpix    = a pixel counter, to be used as index of data values.  Nothing
 *             requires the data or model to be 1D.  kpix runs over the full
 *             pixel count and it's up to this function to interpret its
 *             meaning.
 *   dataptr = pointer into data to be fitted
 *
 *   usr1    = pointer to 1st auxiliary information (for user convenience)
 *             This can be used as a pointer to an independent variable,
 *             for example, e.g. if y_i ~ ymodel(x_i; par).  y_i comes from
 *             dataptr and x_i can come from usr1.
 *   usr2    = pointer to 2nd auxiliary information (for user convenience)
 *
 *   npar    = number of parameters in model (not nec being varied by fitlm)
 *   par     = model parameters
 *
 * -Outputs-
 *   dataval = data value evaluated from dataptr at kpix
 *
 *   model   = model value evaluated at kpix
 *
 *   wgt     = data weight at kpix for chi^2: wgt[kpix] = 1/sigma_i^2;
 *             wgt[kpix] <0 indicates "do not use this point;
 *             wgt = NULL implies use uniform, unity weight 
 *
 *   deriv   = if non-NULL, return vector of model partial derivs wrt
 *             parameters at kpix.  This vector has the full npar
 *             entries, although only parameters that are currently
 *             being varied are used.  It's OK not to fill in entries
 *             that are not being varied.  If deriv is NULL don't
 *             touch!, fitlm will get derivs from model differences.
 *          
 *   return  = 0/error
 */
typedef int LMFUNC(int kpix, void *dataptr, void *usr1, void *usr2, 
		   double *dataval, double *wgt, 
		   int npar, double *par, double *model, double *deriv);

/* Fit a model to data, minimizing chi^2 using Levenberg-Marquardt */
int fitlm(int ndata, void *dataptr, void *usr1, void *usr2, LMFUNC func,
	  int npar, double *par, int *usepar, 
	  double *chi, double *cov, int *ctrl);

/* Give user access to QUITFRAC LAMINIT, etc.  Not for the novice! */
int fitlm_param(int maxiter, int laminit, double quit, double parvary);

/* Solves the matrix equation   Y = A X,   returns X. */
int linsolve(int n, double *y, double *a, double *x);

/* Invert a matrix */
int invert(int n, double *a, double *det);

/* Debug print of status and variables */
int vprint(int n, double *a, double f, int niter, int lambda);

/* Sum chi^2/N, data gradient, and data curvature for this set of parameters */
double chilm(int ndata, void *dataptr, void *usr1, void *usr2, LMFUNC func,
	     int npar, double *par, int *usepar, 
	     int *ctrl, double *alpha, double *beta);

/* Fit a model to data, minimizing chi^2 using Levenberg-Marquardt */
int fitlm(int ndata, void *dataptr, void *usr1, void *usr2, LMFUNC func,
	  int npar, double *par, int *usepar, 
	  double *chi, double *cov, int *ctrl)
{
/* fitlm arguments:
 *
 *    ndata	= number of data points
 *    dataptr   = pointer into data to be fitted
 *
 *    usr1	= pointer to 1st auxiliary information (for user convenience)
 *    usr2	= pointer to 2nd auxiliary information (for user convenience)
 *
 *    func	= LMFUNC function to be evaluated at each pixel, providing
 *                data, uncertainty, and model values, and (possibly)
 *                derivatives of model with respect to parameters
 *
 *    npar	= number of parameters in model
 *    par	= model parameters, initial values input, fitted values output
 *
 *    usepar	= 0/1 indicating whether that parameter is to be varied
 *
 *    chi	= chi^2 for this fit
 *
 *    cov	= error and covariance matrix for fitted parameters.  Note
 *                that this is close-packed and non-varied parameters are
 *                not represented.
 *    ctrl	= control word:
 * 
 *     Input:
 *       ctrl & 0x01: requests reports for each iteration
 *       ctrl & 0x02: LMFUNC returns a vector of analytic first derivatives
 *       ctrl & 0x04: LMFUNC is linear
 *
 *     Output:
 *       ctrl & 0x00ff: final lambda value, offset by 128
 *       ctrl & 0xff00: number of iterations required
 *
 *    return    = 0 successful completion
 *               -1 error from invert: singular matrix
 *               -2 error from invert: matrix too large
 *               -3 too many parameters
 *                N maximum iteration reached (N)
 */
   int j, k, lam, n, nused, on[MAXPAR], err, verbose, doderiv, linfunc;
   double beta[MAXPAR], alpha[MAXPAR*MAXPAR], dpar[MAXPAR];
   double betasave[MAXPAR], alphasave[MAXPAR*MAXPAR];
   double chisq, chinew, det;

   *chi = 0.0;

   verbose = (*ctrl) & FITLM_VERBOSE;
   doderiv = (*ctrl) & FITLM_DERIV;
   linfunc = (*ctrl) & FITLM_LINEAR;

   if(npar > MAXPAR) {
      fprintf(stderr, "fitlm: Parameter count %d exceeds max %d\n", 
	      npar, MAXPAR);
      return(-3);
   }

/* Index to which parameters to vary */
   for(j=nused=0; j<npar; j++) {
      on[nused] = j;
      nused += usepar[j];
   }

   lam = LAMINIT;

/* Initial computation of chi^2 and derivatives */
   chisq = chilm(ndata, dataptr, usr1, usr2, func,
		 npar, par, usepar, ctrl, alpha, beta);
   if(verbose) vprint(npar, par, chisq, 0, lam);

/* Iterate */
   for(n=0; n<MAXITER; n++) {

/* Save a copy of chi^2 derivatives */
      memcpy(alphasave, alpha, nused*nused*sizeof(double));
      memcpy(betasave, beta, nused*sizeof(double));

/* L-M increment to parameters */
      memcpy(cov, alpha, nused*nused*sizeof(double));
      for(j=0; j<nused; j++) {
	 cov[j+j*nused] = ABS((1+pow(LAMBASE, (double)lam)) * cov[j+j*nused]);
      }
      if( (err = linsolve(nused, beta, cov, dpar)) ) {
//	 printf("fitlm: linsolve error %d\n", err);
	 return(err);
      }

/* Evaluate new Chi^2 */
      for(j=0; j<nused; j++) par[on[j]] += dpar[j];
      chinew = chilm(ndata, dataptr, usr1, usr2, func,
		     npar, par, usepar, ctrl, alpha, beta);
      if(verbose) vprint(npar, par, chinew, n, lam);

/* Linear case, we should be done */
      if(linfunc) {
	 *chi = chinew;
	 break;
      }

/* Is it better? */
//      printf(" %d %.7f %.7f %.7f %.7f\n", 
//      lam, chisq, chinew, chisq-chinew, QUITFRAC*chisq);
      if(chinew < chisq) {	/* Yes, keep new params and diminish lam */
/* Are we done? */
	 *chi = chinew;
	 if(chisq - chinew < QUITFRAC*chisq) break;
	 chisq = chinew;
	 lam--;

      } else {			/* No, put parameters back and increase lam */
	 *chi = chisq;
	 for(j=0; j<nused; j++) par[on[j]] -= dpar[j];
	 memcpy(alpha, alphasave, nused*nused*sizeof(double));
	 memcpy(beta, betasave, nused*sizeof(double));
	 lam++;
      }
   }

/* Load up control word with iteration count and final lambda */
   *ctrl = n*256 + (lam+128);

/* Compute the parameter errors and covariance matrix */
   memcpy(cov, alpha, nused*nused*sizeof(double));

   if( (err=invert(nused, cov, &det)) != 0) {
//      fprintf(stderr, "Invert error, det= %.4e\n", det);
      return(err);
   }
   for(j=0; j<nused; j++) {
      cov[j+j*nused] = (cov[j+j*nused] > 0) ? 
	 sqrt(cov[j+j*nused]) : -sqrt(-cov[j+j*nused]);
   }

   for(j=1; j<nused; j++) {
      for(k=0; k<j; k++) cov[j+k*nused] = cov[k+j*nused] = 
			    cov[k+j*nused] / (cov[k+k*nused]*cov[j+j*nused]);
   }

   return( (n==MAXITER) ? n : 0);
}

/* Sum chi^2, gradient, and curvature for this set of parameters */
/*
 * chilm    = Sum wgt_i (y_i-model(x_i))^2
 * beta_j   = Sum wgt_i (y_i-model(x_i)) dchi/da_j packed for parameters in use
 * alpha_jk = Sum wgt_i dchi/da_j dchi/da_k        packed for parameters in use
 */

double chilm(int ndata, void *dataptr, void *usr1, void *usr2, LMFUNC func,
	     int npar, double *par, int *usepar, 
	     int *ctrl, double *alpha, double *beta)
{
   int i, j, k, nused, on[MAXPAR], err, ndof, linfunc, doderiv, verbose, dowgt;
   double chi, dataval, model, dmodel, dfda[MAXPAR], dfwt, df, dp;
   double *wptr=NULL, wgt=1.0;

   verbose = (*ctrl) & FITLM_VERBOSE;
   doderiv = (*ctrl) & FITLM_DERIV;
   linfunc = (*ctrl) & FITLM_LINEAR;
   dowgt   = (*ctrl) & FITLM_WGT;
   if(dowgt) wptr = &wgt;

/* Count the parameters in use, make an index of them, zero alpha and beta */
   for(j=nused=0; j<npar; j++) {
      on[nused] = j;
      nused += usepar[j];
      beta[j] = 0.0;
      for(k=0; k<npar; k++) alpha[j+k*npar] = 0.0;
   }

/* Sum all the data points */
   chi = 0.0;
   ndof = -nused;
   for(i=0; i<ndata; i++) {

/* This is all we need if analytic derivatives are provided */
      if(doderiv) {
	 err = func(i, dataptr, usr1, usr2, 
		    &dataval, wptr, npar, par, &model, dfda);

/* Numerical first derivatives of function */
      } else {
	 err = func(i, dataptr, usr1, usr2, 
		    &dataval, wptr, npar, par, &model, NULL);
	 for(j=0; j<nused; j++) {
	    dp = PARVARY * ((par[on[j]] == 0.0) ? 1.0 : par[on[j]]);
	    par[on[j]] += dp;
	    err = func(i, dataptr, usr1, usr2, 
		       &dataval, wptr, npar, par, &dmodel, NULL);
	    par[on[j]] -= dp;
	    dfda[on[j]] = (dmodel-model) / dp;
	 }
      }

/* Sum derivatives of chi^2 wrt parameters (drop factor of 2) */
      if(err == 0 && wgt > 0) {
	 ndof++;
	 df = dataval - model;
	 for(j=0; j<nused; j++) {
	    dfwt = wgt * dfda[on[j]];
	    for(k=0; k<=j; k++) alpha[j+k*nused] += dfwt * dfda[on[k]];
	    beta[j] += dfwt * df;
	 }
	 chi += wgt * df*df;
      }
   }
/* Symmetrize second derivative matrix */
   for(k=1; k<nused; k++) {
      for(j=0; j<k; j++) alpha[j+k*nused] = alpha[k+j*nused];
   }

#ifdef SHOWVAL
   int *aux=(int *)usr1;
   for(j=aux[0]-1; j>=0; j--) {
      printf("%5d", j);
      for(i=0; i<aux[1]; i++) {
	 k = i + aux[1]*j;
	 err = func(k, dataptr, usr1, usr2, 
		    &dataval, wptr, npar, par, &model, dfda);
	 printf(" %5.0f", model);
      }
      for(i=0; i<aux[1]; i++) {
	 k = i + aux[1]*j;
	 err = func(k, dataptr, usr1, usr2, 
		    &dataval, wptr, npar, par, &model, dfda);
	 printf(" %5.0f", dataval-model);
//	 printf(" %5.0f", wgt*(dataval-model)*(dataval-model));
//	 printf(" %5.0f", dataval);
      }
      printf("\n");
   }
#endif

/* Return chi^2 */
   return(chi);
}

/* Give user access to QUITFRAC LAMINIT, etc via a call? */
int fitlm_param(int maxiter, int laminit, double quit, double parvary)
{
   MAXITER = maxiter;
   LAMINIT = laminit;
   QUITFRAC = quit;
   PARVARY = parvary;
   return(0);
}


int vprint(int n, double *a, double f, int niter, int lambda)
{
   int i;
   for(i=0; i<n; i++) {
      if((i%7) == 0) printf(" A(I) =");
      printf("%11.4g", a[i]);
      if((i%7) == 6) printf("\n");
   }
   if(((i-1)%7) != 6) printf("\n");
   printf(" F = %14.7g     ITER = %3d    LAMBDA = %3d\n", f, niter, lambda);
   return(0);
}

/* Sample main() to illustrate useage */
#ifdef MAINTEST
#define MAXPT 200

/* Function that is fitted to the data */
int gauss(int kpix, void *dataptr, void *usr1, void *usr2, 
	     double *dataval, double *wgt, 
	     int npar, double *par, double *model, double *deriv);

/* Show us the covariance matrix */
int printmat(int n, double *a, double *y)
{
   int i, j;
   for(i=0; i<n; i++) {
      for(j=0; j<n; j++) printf(" %12.4g", a[j+i*n]);
      printf("      %12.4g\n", y[i]);
   }
   return(0);
}

/* Create some data and fit a Gaussian to it */
int main(int argc, char **argv)
{
   int i, k, n, err, npt=MAXPT, npar, usepar[4], ctrl;
   float x[MAXPT], y[MAXPT];
   double sig, ctr, ampl, off, par[4], cov[4*4], *dummy, grand, init[4], chi;

////////////////////////////////////////////////////////////////
/* Create data from four params:  off + ampl * exp(-0.5 * ((x-ctr)/sig)^2) */
   off = 3.0;
   ampl = 9.0;
   ctr = 47.0;
   sig = 18.0;
   for(i=0; i<npt; i++) {
/* Reasonable spread of independent variable */
      x[i] = 10*sig*(random()/(double)RAND_MAX-0.5) + ctr;
/* Gaussian noise of RMS=1 */
      for(k=0, grand=-6; k<12; k++) grand += random()/(double)RAND_MAX;
/* Data approximating our Gaussian model */
      y[i] = off + ampl * exp(-0.5*pow((x[i]-ctr)/sig, 2.0)) + 2 * grand;
/* Tell us about the points? */
//      printf("%8.3f %8.3f %5d\n", x[i], y[i], i);
   }
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/* How many parameters does the function have? */
   npar = 4;
   for(i=0; i<npar; i++) usepar[i] = 1; /* Fit all parameters */
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/* Initialize parameters -- this is worth doing reasonably well... */
   init[0] = init[1] = y[0];
   init[2] = x[0];
   for(i=1; i<npt; i++) {
      if(init[1] < y[i]) {
	 init[1] = y[i];		/* init[1] = ampl = max */
	 init[2] = x[i];		/* init[2] = ctr = location of max */
      }
      if(init[0] > y[i]) init[0] = y[i];	/* init[0] = offset = min */
   }

   init[1] -= init[0];	/* init[1] = ampl = max - min */
   init[3] = 0.0;	/* init[3] = sig = <RMS> of all points above half max */
   for(i=n=0; i<npt; i++) {
      if(y[i] > 0.5*init[1] + init[0]) {
	 init[3] += (x[i]-init[2]) * (x[i]-init[2]);
	 n++;
      }
   }
   init[3] = (init[3] > 0) ? sqrt(init[3]/n) : 1.0;

   printf("True parameters:    %8.3f %8.3f %8.3f %8.3f\n", off, ampl, ctr, sig);
   printf("Initial parameters: %8.3f %8.3f %8.3f %8.3f\n", 
	  init[0], init[1], init[2], init[3]);
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/* Fit the parameters to the data, no weights */
   printf("Numerical derivatives, no weights:\n\n");
   for(i=0; i<npar; i++) par[i] = init[i];
   
   ctrl = FITLM_VERBOSE;
   err = fitlm(npt, (void *)y, (void *)x, (void *)dummy, gauss,
	       npar, par, usepar, &chi, cov, &ctrl);
   printf("  Error/covariance matrix and params:\n");
   printmat(npar, cov, par);
   printf("  chi= %.2f, lam= %d, niter= %d\n", 
	  chi, (ctrl&0xff)-128, ctrl/256);
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/* Fit the parameters to the data, weights */
   printf("\nAnalytic derivatives, weights:\n\n");
   for(i=0; i<npar; i++) par[i] = init[i];
   ctrl = FITLM_VERBOSE + FITLM_DERIV + FITLM_WGT;
   err = fitlm(npt, (void *)y, (void *)x, (void *)dummy, gauss,
	       npar, par, usepar, &chi, cov, &ctrl);
   printf("  Error/covariance matrix and params:\n");
   printmat(npar, cov, par);
   printf("  chi= %.2f, lam= %d, niter= %d\n", 
	  chi, (ctrl&0xff)-128, ctrl/256);
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/* Fit the parameters to the data, no weights */
   printf("\nKeep one parameter fixed, no weights:\n\n");
   for(i=0; i<npar; i++) par[i] = init[i];
   par[2] = ctr + 0.5;
   usepar[2] = 0;
   ctrl = FITLM_VERBOSE + FITLM_DERIV;
   err = fitlm(npt, (void *)y, (void *)x, (void *)dummy, gauss,
	       npar, par, usepar, &chi, cov, &ctrl);
   printf("  chi= %.2f, lam= %d, niter= %d\n", 
	  chi, (ctrl&0xff)-128, ctrl/256);
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/* Fit the parameters to the data, weights */
   printf("\nLinear case (hold non-linear parameters fixed):\n\n");
   for(i=0; i<npar; i++) par[i] = init[i];
   par[2] = ctr;
   par[3] = sig;
   usepar[2] = usepar[3] = 0;
   ctrl = FITLM_VERBOSE + FITLM_DERIV + FITLM_LINEAR + FITLM_WGT;
   err = fitlm(npt, (void *)y, (void *)x, (void *)dummy, gauss,
	       npar, par, usepar, &chi, cov, &ctrl);
   printf("  chi= %.2f, lam= %d, niter= %d\n", 
	  chi, (ctrl&0xff)-128, ctrl/256);
////////////////////////////////////////////////////////////////

   exit(0);
}


int gauss(int kpix, void *dataptr, void *usr1, void *usr2, 
	     double *dataval, double *wgt, 
	     int npar, double *par, double *model, double *deriv)
{
   float *data=(float *)dataptr, *x=(float *)usr1;
   double gexp, dx, sig=par[3];
   *dataval = data[kpix];
   if(wgt != NULL) *wgt = 1/(2.0*2.0); 	/* Secret info, could use usr2 */
   dx = x[kpix] - par[2];
   gexp = exp(-0.5*dx*dx/sig/sig);
   *model =  par[0] + par[1] * gexp;

   if(deriv != NULL) {
      deriv[0] = 1.0;
      deriv[1] = gexp;
      deriv[2] = par[1] * gexp * dx/sig/sig;
      deriv[3] = par[1] * gexp * dx*dx/sig/sig/sig;
   }

   return(0);
}

#endif
