/* Use cfitsio to read an image into a floating array */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "fitsio.h"

#define NINT(x) (x<0?(int)((x)-0.5):(int)((x)+0.5))
#define ABS(x)  ((x)<0?(-(x)):(x))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

// Macro from robospect to do all the boring fits error handling
#define FRE { fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); } ; status = 0; }

#define FFATAL { fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); exit(1); }  }


/* Print out cfitsio error messages and exit program */
void printerror(int status);

/* Read a 16-bit FITS file */
int cf_readreal(char **head, int *nx, int *ny, int *nz, float **data,
		int rd3d, char *file);

/* Read a FITS file */
void cf_read(fitsfile **fptr, char *fname, int type, 
	     char **head, int *nx, int *ny, int *nz, int rd3d, void **data);

/* Read a FITS file into a floating array */
int cf_readreal(char **head, int *nx, int *ny, int *nz, float **data,
		int rd3d, char *file)
{
   fitsfile *ff;   /* pointer to the FITS file */
   
/* Read the FITS file */
   cf_read(&ff, file, TFLOAT, head, nx, ny, nz, rd3d, (void **)data);

   return(0);
}

// JTFIXME
static int VERBOSE=0;


/* Read a FITS file */
void cf_read(fitsfile **fptr, char *fname, int type, 
	     char **head, int *nx, int *ny, int *nz, int rd3d, void **data)
{
   long int naxes[10], fpixel[10], npix;
   LONGLONG nelements;
   int bitpix, status, naxis, anynull;
   int nkeys, keypos, i;
   void *dataptr;

/* Read the original file: use open_image vs open_file in case compressed */
   status = 0;
   if ( fits_open_image(fptr, fname, READONLY, &status) ) FFATAL;

/* What about it? */
   if( fits_get_img_param(*fptr, 3, &bitpix, &naxis, naxes, &status) ) FFATAL;

   if(VERBOSE > 1) {
      if(naxis <= 2) {
	 printf("Input file %s has %d axes %ld x %ld and bitpix %d\n", 
		fname, naxis, naxes[0], naxes[1], bitpix);
      } else {
	 printf("Input file %s has %d axes %ld x %ld x %ld and bitpix %d\n", 
		fname, naxis, naxes[0], naxes[1], naxes[2], bitpix);
      }
   }

   *nx = naxes[0];
   *ny = naxes[1];
   *nz = 0;

   nelements = naxes[0] * naxes[1];          /* number of pixels to read */
/* Do a monolithic read of both data and variance */
   npix = nelements;
   if(rd3d && naxis == 3) {
      *nz = naxes[2];
      npix = nelements * naxes[2];
   }

/* Create input array and read it */
   if(type == TUSHORT) {
      *data = (ushort *)calloc(npix, sizeof(ushort));
   } else if(type == TSHORT) {
      *data = (short int *)calloc(npix, sizeof(short int));
   } else if(type == TFLOAT) {
      *data = (float *)calloc(npix, sizeof(float));
   } else {
      fprintf(stderr, "Unknown type %d\n", type);
      exit(1);
   }

   fpixel[0] = fpixel[1] = fpixel[2] = 1;    	     /* first pixel to read */
   if(! (rd3d && naxis==3) ) {
      if( fits_read_pix(*fptr, type, fpixel, nelements, NULL,
			*data, &anynull, &status) ) {
	 if(status != NUM_OVERFLOW) FFATAL;
	 if(VERBOSE > 1) FRE;
	 status = 0;
      }
   } else {
      for(i=0; i<naxes[2]; i++) {
	 fpixel[2] = i+1;
	 dataptr = (float *)(*data)+i*nelements;
	 if(type == TUSHORT || type == TSHORT) {
	    dataptr = (ushort *)(*data)+i*nelements;
	 }
	 if( fits_read_pix(*fptr, type, fpixel, nelements, NULL,
			   dataptr, &anynull, &status) ) {
	    if(status != NUM_OVERFLOW) FFATAL;
	    if(VERBOSE > 1) FRE;
	    status = 0;
	 }
      }

/* Alternate read scheme */
#if 0
/* cfitsio does not object to reading beyond nx*ny, so long as data exist */
      if( fits_read_pix(*fptr, type, fpixel, naxes[2]*nelements, NULL,
			*data, &anynull, &status) ) {
	 if(status != NUM_OVERFLOW) FFATAL;
	 if(VERBOSE > 1) FRE;
	 status = 0;
      }
#endif
   }

/* Get number of keywords */
   if (fits_get_hdrpos(*fptr, &nkeys, &keypos, &status) ) {
      if(VERBOSE > 1) FRE;
      status = 0;
   }

/* Read the header */
   *head = malloc((nkeys+2)*80);
   for(i=1; i<=nkeys; i++)  {
      if ( fits_read_record(*fptr, i, *head+(i-1)*80, &status) ) {
	 if(VERBOSE > 1) FRE;
	 status = 0;
      }
   }

   if(status) exit(1);
}

/* Read a 16-bit FITS file */
/* Print out cfitsio error messages and exit program */
void printerror(int status)
{
   if (status) {
      fits_report_error(stderr, status); /* print error report */
      exit( status );    /* terminate the program, returning error status */
   }
   return;
}
