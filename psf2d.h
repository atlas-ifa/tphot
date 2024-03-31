/* psf2d.h -- header file for psf2d.c: fit a 2D Waussian profile */
/* 140621 add Bayer pixel terms */
/* 120929 v1.0 John Tonry */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))
#define NINT(x) (x<0?(int)((x)-0.5):(int)((x)+0.5))

typedef struct {
   double x0;		/* x0 of fit (TSK convention!) */
   double y0;		/* y0 of fit */
   double peak;		/* peak of Waussian */
   double sky;		/* sky of Waussian fit */
   double sx2;		/* sx2 */
   double sxy;		/* sxy */
   double sy2;		/* sy2 */
   double beta4;	/* beta4 */
   double beta6;	/* beta6 */
   double major;	/* major axis fwhm of Waussian fit */
   double minor;	/* minor axis fwhm of Waussian fit */
   double phi;		/* angle of major axis (CCW from x axis) [rad] */
   double chin;		/* chi^2/Ndof of fit */
   double rms;		/* rms residual within a radius of 2 sigma */
   double absresid;	/* average abs residual within a radius of 2 sigma */
   double maxresid;	/* max residual within a radius of 2 sigma */
   double snr;		/* detection SNR */
   double snrms;	/* detection waffle-stomp noise */
   double dx0;		/* uncertainty in x0 of fit (TSK convention!) */
   double dy0;		/* uncertainty in y0 of fit */
   double dpeak;	/* uncertainty in peak of Waussian */
   double dsky;		/* uncertainty in sky of Waussian fit */
   double dsx2;		/* uncertainty in sx2 */
   double dsxy;		/* uncertainty in sxy */
   double dsy2;		/* uncertainty in sy2 */
   double dbeta4;	/* uncertainty in beta4 */
   double dbeta6;	/* uncertainty in beta6 */
   int ninit;		/* How many params well initialized? 0/2/5 */
   int niter;		/* Number of iterations in fit */
   int nignore;		/* Number of nearby pixels with ignore value */
   int bayer;		/* Bayer pixel? 0 => G in 00,11; 1 => G in 01,10 */
   int fitclr;		/* Bayer pixel? 0 => leave fixed; 1 => fit */
} PSF2D_PARAM;

typedef struct {
   double flux;		/* flux */
   double dflux;	/* flux uncertainty */
   double sky;		/* sky level */
   double dsky;		/* sky uncertainty */
   double skyrms;	/* rms in sky */
} PSF2D_FLUX;

typedef struct {
   double cx;		/* [pix] mean center x0 */
   double cy;		/* [pix] mean center y0 */
   int xpk;		/* [pix] highest pixel x */
   int ypk;		/* [pix] highest pixel y */
   double rkron;	/* [pix] Kron radius (linear) */
   double rvar;		/* [pix] Variance radius */
   double flux;		/* [ADU] Flux within 2.5 R_Kron */
   double q1cos;	/* [pix^2] linear offset wrt peak */
   double q1sin;	/* [pix^2] linear offset wrt peak */
   double q2;		/* [pix^2] quadrupole = qxx+qyy */
   double q2cos;	/* [pix^2] quadrupole at phi=0 = qxx-qyy */
   double q2sin;	/* [pix^2] quadrupole at phi=45 = 2*qxy */
   double q3cos;	/* [pix^2] trefoil at phi=0 */
   double q3sin;	/* [pix^2] trefoil at phi=30 */
   double crad;		/* [pix^2] center wrt radial from img center */
   double qrad;		/* [pix^2] quadrupole wrt radial from img center */
   double trad;		/* [pix^2] trefoil wrt radial from img center */
} PSF2D_MOM;

/* Return error code 0: no error, 
 *                  -1: off image, 
 *                   1: fit didn't converge,
 *                N*10: N zero'ed pixels within FWHM
 */
int wpsf(int nx, int ny, float *a, float *wgt,
	 double eadu, double sat, float badata,
          int aprad, int skyrad, int nwpar, 
          PSF2D_PARAM *wpar, PSF2D_FLUX *flux);

/* Fit a trailed PSF */
int trailpsf(int nx, int ny, float *a, float *wgt,
	     double eadu, double sat, float badata,
	     int aprad, int skyrad, int nwpar, 
	     PSF2D_PARAM *wpar, PSF2D_FLUX *flux);
