/*
 * rwfits 2.2 - 031218 JT Fix a little bug in addfitsheader(); add rfitshead()
 * rwfits 2.1 - 11/02/99 JT Merge in bug fixes from Saurabh Jha
 * rwfits 2.0 - 9/21/99 JT Add header and convenience routines like rfitsreal()
 * rwfits 1.1 - 9/07/99 Add code for overlapping swab()
 * rwfits 1.0 - 7/08/92  initial rev
 *
 * Routines to read and write disk FITS format data.
 * These files have a header which consists of 36*n 80 byte header records
 * which are a duplicate of the FITS tape format, with a mandatory END record.
 * This is followed by N*2880 bytes of data.
 *
 * The conversion routines look to see whether LOWENDIAN and VAXFP are defined
 * 051113: I've stuck in Sidik's clever endian test instead; LOWENDIAN not nec.
 * Set them if necessary!
 * E.g.    VAX:        -DLOWENDIAN -DVAXFP
 *         DEC MIPS:   -DLOWENDIAN
 *         Intel:      -DLOWENDIAN
 *         Sun:
 *
 * NOTE: problems have arisen with DEC Alpha's running OSF1 because OSF1
 * has a swab() which doesn't work in place!  You can define OSF_ALPHA to
 * cause a swab() to be compiled here which works right:
 *         OSF1 ALPHA:  -DLOWENDIAN -DOSF_ALPHA
 *
 *
 * rwfits.c: John Tonry 7/8/92
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

double atof();

/* Simple run-time test for the host byte-order, and swapping routines. */
static int _ENDIAN_TEST = 1;
#define LOWENDIAN() (*(char*)&_ENDIAN_TEST)

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define NFITS (2880)
#define MAXBUF (65536)

static char buf[MAXBUF];

#define ERR_CANT_OPEN_FILE		1
#define ERR_CANT_READ_HEADER		2
#define ERR_CANT_READ_DATA		3
#define ERR_CANT_WRITE_HEADER		4
#define ERR_CANT_WRITE_DATA		5
#define ERR_INSUFFICIENT_HEADER		6
#define ERR_NO_NAXIS			7
#define ERR_NO_END			8
#define ERR_INCONSISTENT_NAXIS		9
#define ERR_BUFFER_OVERFLOW	       10
#define ERR_NOT_A_FITS_FILE	       11

int testfitslen_(char *file, int *mefits, int *bitpix, 
		 int *naxis, int *dims, int *nhead, int filelen);
int openc_(int *fd, char *file, int *mode, int filelen);
int closec_(int *fd);
int writebytes_(int *fd, char *data, int *nbyte);
int readbytes_(int *fd, char *data, int *nbyte);
int testfits_(char *file, int *mefits, int *bitpix, int filelen);
int rhead_(int *fd, int *nhead, char *head, int headlen);
int whead_(int *fd, int *nhead, char *head, int headlen);
int parsehead_(int *nhead, char *head, int *bitpix, int *nx, int *ny, 
	       double *scale, double *zero, int headlen);
// int ieeefp_(int *npix, char *data, char *fp);
// int fpieee_(int *npix, char *fp, char *data);
// int FITSorder(int bitpix, int npix, short *data);
int counthead(char *header);


/*
 * rfitsreal() reads FITS image from disk and returns a header and real array
 */
int rfitsreal(head, nx,ny, nz, data, rd3d, file)
char **head;		/* header */
char *file;		/* file name */
float **data;		/* image data */
int *nx, *ny, *nz;	/* NAXIS1 and NAXIS2, possibly NAXIS3 from header */
int rd3d;		/* Read 3D NAXIS3 if it exists? */
{
  int fd, mode=0, err, nbyte, nhead, npix;
  double scale, zero;	/* BSCALE and BZERO parameters from header */
  int bitpix;		/* BITPIX for disk and resultant array */
  int mefits, naxis, dims[10], i;

  testfitslen_(file, &mefits, &bitpix, &naxis, dims, &nhead, strlen(file));

  if(mefits != 1) {
    fprintf(stderr, "File '%s' does not appear to be a FITS file!\n", file);
    return(ERR_NOT_A_FITS_FILE);
  }

/* Create an array for the FITS header (allocate 87 extra as well) */
  *head = (char *) malloc(80*(nhead+87));
  nhead = nhead + 87;
  bzero(*head, 80*nhead);

  if( (err=openc_(&fd, file, &mode, strlen(file)+1)) != 0) return(err);

  if( (err=rhead_(&fd, &nhead, *head, strlen(*head)+1)) != 0) {close(fd); return(err);}
  if( (err=parsehead_(&nhead, *head, &bitpix, nx,ny, &scale, &zero,strlen(*head)+1)) != 0) {
    close(fd);    return(err);  }

  if(ABS(bitpix) > 32) {
    close(fd); 
    fprintf(stderr, "Error incompatible BITPIX = %d\n", bitpix);
    return(ERR_CANT_READ_DATA);
  }

  *nz = 0;
  if(rd3d) {
     for(i=0; i<80*nhead; i+=80) {
	if( strncmp(*head+i, "NAXIS3", 6) == 0) *nz = atoi(*head+i+10);
     }
  }

  npix = (*nx) * (*ny);
  if(rd3d && *nz>0) npix *= (*nz);
  nbyte = npix * (ABS(bitpix)/8);
  if(bitpix == 1) nbyte = 2*((npix+15)/16);

  *data = (float *) calloc(npix+17, sizeof(float));

  if( (err=read(fd, *data, nbyte)) != nbyte) {
    close(fd); 
    fprintf(stderr, "Error: rwfits read %d but got %d\n", nbyte, err);
    return(ERR_CANT_READ_DATA);
  }
  close(fd);

  if(bitpix == -32) {
    ieeefp_(&npix, *data, *data);

  } else if(bitpix == 16) {
    shortfp_(&npix, &scale, &zero, *data, *data);

  } else if(bitpix == -16) {
    unsignedfp_(&npix, &scale, &zero, *data, *data);

  } else {
    fprintf(stderr,"I cannot deal with BITPIX = %d\n", bitpix);
    return(1);
  }

/* Update the BITPIX to reflect the actual image */
  err = chfitshead(&i, *head, "BITPIX  ", "INTEGER", -32, 0.0);

  return(0);
}

/*
 * wfitsreal() will write a real FITS image onto disk.
 */
int wfitsreal(head, data, file)
char *head;		/* header */
char *file;		/* file name */
char *data;		/* image data */
{
  int fd, mode=1, err, nx, ny, bitpix;
  int nwrite, nbyte, nwrit;
  long long int imagebyte;
  int nhead;
  double scale, zero;

  if( (err=openc_(&fd, file, &mode, strlen(file)+1)) != 0) return(err);

  nhead = counthead(head);

  if( (err=parsehead_(&nhead,head,&bitpix,&nx,&ny,&scale,&zero,strlen(head)+1)) != 0) {
    close(fd);    return(err);  }

  if( (err=whead_(&fd, &nhead, head, strlen(head)+1)) != 0) {
    close(fd);  return(err); }

  imagebyte = nx*ny * (ABS(bitpix)/8);
  if(bitpix == 1) imagebyte = 2*((nx*ny+15)/16);
  nbyte = imagebyte;		       	/* Bytes to write */

  nwrit = 0;
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF,nbyte);
/* Convert from machine specific floating point to FITS IEEE format */
    nx = nwrite/4;
    fpieee_(&nx, &data[nwrit], buf);

    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
    nbyte -= nwrite;
    nwrit += nwrite;
  }

  nwrite = NFITS - (imagebyte%NFITS);
  if(nwrite < NFITS) {
    bzero(buf,nwrite);
    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
  }

  close(fd);
  return(0);
}

/* wfitsshort() writes a real image to a 16 bit FITS image on disk */
/* BITPIX = 16 or -16 for signed/unsigned 16 bit integers */
int wfitsshort(char *head, char *data, char *file)
{
  int fd, mode=1, err, nx, ny, bitpix;
  int nwrite, nbyte, imagebyte, nwrit;
  int nhead;
  double scale, zero;

  if( (err=openc_(&fd, file, &mode, strlen(file)+1)) != 0) return(err);

  nhead = counthead(head);

  if( (err=parsehead_(&nhead,head,&bitpix,&nx,&ny,&scale,&zero,strlen(head)+1)) != 0) {
    close(fd);    return(err);  }

  if( (err=whead_(&fd, &nhead, head, strlen(head)+1)) != 0) {
    close(fd);  return(err); }

  imagebyte = nx*ny * (ABS(bitpix)/8);
  if(bitpix == 1) imagebyte = 2*((nx*ny+15)/16);
  nbyte = imagebyte;		       	/* Bytes to write */

  nwrit = 0;
/* Invert scale and zero for conversion routines */
  zero = -zero;
  scale = 1.0 / scale;
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF, nbyte);

/* Convert from machine specific floating point to 16 bit format */
    nx = nwrite/2;
    if(bitpix == 16) {
      fpshort_(&nx, &scale, &zero, &data[2*nwrit], buf);
    } else if(bitpix == -16) {
      fpunsigned_(&nx, &scale, &zero, &data[2*nwrit], buf);
    }

    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
    nbyte -= nwrite;
    nwrit += nwrite;
  }

  nwrite = NFITS - (imagebyte%NFITS);
  if(nwrite < NFITS) {
    bzero(buf,nwrite);
    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
  }

  close(fd);
  return(0);
}

/*
 * wfitsappend() will write or append to a real FITS image onto disk.
 */
int wfitsappend(head, data, file, fd)
char *head;		/* header */
char *file;		/* file name */
char *data;		/* image data */
int *fd;		/* -1 to write a new file; fd to append */
{
  int mode=1, err, nx, ny, bitpix;
  int nwrite, nbyte, imagebyte, nwrit;
  int nhead;
  double scale, zero;

  if( *fd < 0) {
     if( (err=openc_(fd, file, &mode, strlen(file)+1)) != 0) return(err);
  }

  nhead = counthead(head);

  if( (err=parsehead_(&nhead,head,&bitpix,&nx,&ny,&scale,&zero,strlen(head)+1)) != 0) {
     if( err != ERR_NO_NAXIS) {
	close(*fd);
	return(err);
     }
  }

  if( (err=whead_(&*fd, &nhead, head, strlen(head)+1)) != 0) {
    close(*fd);  return(err); }

  imagebyte = nx*ny * (ABS(bitpix)/8);
  if(bitpix == 1) imagebyte = 2*((nx*ny+15)/16);
  nbyte = imagebyte;		       	/* Bytes to write */

  nwrit = 0;
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF,nbyte);
/* Convert from machine specific floating point to FITS IEEE format */
    nx = nwrite/4;
    fpieee_(&nx, &data[nwrit], buf);

    if( (err=write(*fd, buf, nwrite)) != nwrite) {
      close(*fd); return(ERR_CANT_WRITE_DATA); }
    nbyte -= nwrite;
    nwrit += nwrite;
  }

  nwrite = NFITS - (imagebyte%NFITS);
  if(nwrite < NFITS) {
    bzero(buf,nwrite);
    if( (err=write(*fd, buf, nwrite)) != nwrite) {
      close(*fd); return(ERR_CANT_WRITE_DATA); }
  }
  return(0);
}

/**************************************/
/* Below are header altering routines */
/**************************************/

/*
 * newfitshead() will create a minimal FITS header
 */
int newfitshead(header,bitpix,nx,ny,object)
char **header;		/* header */
int bitpix, nx, ny;	/* BITPIX, NAXIS1, NAXIS2 header entries */
char *object;		/* object name */
{
  int i;
  *header = (char *) malloc(NFITS+1);
  for(i=0; i<NFITS; i++) (*header)[i] = ' ';
  (*header)[NFITS] = '\0';	/* Just to keep a strlen() happy */
  wfitem(0, *header, "SIMPLE  ", "T", 0, 0.0);
  wfitem(1, *header, "BITPIX  ", NULL, bitpix, 0.0);
  wfitem(2, *header, "NAXIS   ", NULL, 2, 0.0);
  wfitem(3, *header, "NAXIS1  ", NULL, nx, 0.0);
  wfitem(4, *header, "NAXIS2  ", NULL, ny, 0.0);
  wfitem(5, *header, "OBJECT  ", object, 0, 0.0);
  wfitem(6, *header, "END     ", NULL, 0, 0.0);
  return(0);
}

/*
 * addfitshead() adds a new line to a FITS header
 */
int addfitshead(n, header, keyword, cvalue, ivalue, rvalue)
int *n;		/* Number of header entries */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *cvalue;	/* (possible) character value to set (if non NULL) */
int ivalue;	/* (possible) integer value to set (if cvalue = "INTEGER") */
double rvalue;	/* (possible) real value to set (if cvalue = "FLOAT") */
{
  if(strncmp(&header[80*(*n-1)], "END     ", 8) != 0) {
    fprintf(stderr,"addfitshead: header count is wrong; adjusting count\n");
    for(*n=0; *n<strlen(header)/80; *n++) {
      if(strncmp(&header[80*(*n)], "END     ", 8) == 0) break;
    }
    *n += 1;
  }
  wfitem(*n-1, header, keyword, cvalue, ivalue, rvalue);
  wfitem(*n,   header, "END     ", NULL, 0, 0.0);
  *n += 1;
  return(0);
}

/*
 * chfitshead() changes a line in a FITS header
 */
int chfitshead(n, header, keyword, cvalue, ivalue, rvalue)
int *n;		/* Line where keyword was found in header (-1 if not found) */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *cvalue;	/* (possible) character value to set (if non NULL) */
int ivalue;	/* (possible) integer value to set (if cvalue = "INTEGER") */
double rvalue;	/* (possible) real value to set (if cvalue = "FLOAT") */
{
  int i;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) {
      *n = -1;
      return(1);
    }
  }
  wfitem(i, header, keyword, cvalue, ivalue, rvalue);
  *n = i;
  return(0);
}

/*
 * ifitshead() reads an integer component of a FITS header
 */
int ifitshead(header, keyword, ivalue)
char *header;	/* Header */
char *keyword;  /* Keyword to write */
int *ivalue;	/* integer value retrieved */
{
  int i;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) return(1);
  }
  *ivalue = atoi(&header[80*i+9]);
  return(0);
}

/*
 * rfitshead() reads a double component of a FITS header
 */
int rfitshead(header, keyword, rvalue)
char *header;	/* Header */
char *keyword;  /* Keyword to write */
double *rvalue;	/* double value retrieved */
{
  int i;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) return(1);
  }
  *rvalue = atof(&header[80*i+9]);
  return(0);
}

/*
 * cfitshead() reads a character component of a FITS header
 */
int cfitshead(header, keyword, value)
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *value;	/* character string retrieved */
{
  int i;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) return(1);
  }
  i = sscanf(header+80*i+9, "%s", value);
  return(0);
}

/*
 * chfitscomment() changes the comment of a line in a FITS header
 */
int chfitscomment(n, header, keyword, comment)
int *n;		/* Line where keyword was found in header (-1 if not found) */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *comment;	/* comment to append (30 data + ' / ' + <=47 comment */
{
  int i;
  char *c;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) {
      *n = -1;
      return(1);
    }
  }
  *n = i;
  c = &header[80*i] + 30;
  *c++ = ' ';    *c++ = '/';    *c++ = ' ';
  for(i=0; i<47; i++) *c++ = i < strlen(comment) ? comment[i] : ' ';
  return(0);
}

/*
 * wfitem() writes a specific line to a FITS header
 */
int wfitem(n, header, keyword, cvalue, ivalue, rvalue)
int n;		/* Index of header entry to overwrite */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *cvalue;	/* (possible) character value to set (if non NULL) */
int ivalue;	/* (possible) integer value to set (if cvalue = "INTEGER") */
double rvalue;	/* (possible) real value to set (if cvalue = "FLOAT") */
{
  int i;
  if(strncmp(keyword, "END     ", 8) == 0) {
    i = sprintf(&header[80*n], "%-8s", keyword);
  } else if(strncmp(keyword, "SIMPLE  ", 8) == 0) {
    i = sprintf(&header[80*n], "%-8s= %20s", keyword,cvalue);
  } else {
    if(cvalue == NULL || (strncmp(cvalue, "INTEGER", 7)==0)) {
      i = sprintf(&header[80*n], "%-8s= %20d", keyword, ivalue);
    } else if(strcmp(cvalue, "FLOAT") == 0 || strcmp(cvalue, "REAL") == 0) {
      i = sprintf(&header[80*n], "%-8s= %20.8g", keyword, rvalue);
    } else {
      i = sprintf(&header[80*n], "%-8s= %20s", keyword, cvalue);
    }
  }
  if (i < 0) {
    fprintf(stderr,"wfitem: Error writing header string '%s'!\n", keyword);
    return(1);
  }
  while (i < 80) {    /* overwrite any lurking \0's */
    if (header[80*n+i] == '\0') header[80*n+i] = ' ';
    i++;
  }
  return(0);
}

/*
 * counthead() returns the number of lines in a FITS header
 */
int counthead(header)
char *header;
{
  int i;
  for(i=0; i<strlen(header); i+=80) {
    if(strncmp(&header[i], "END     ", 8) == 0) break;
  }
  return(i/80+1);
}




/*******************************************/
/* Below are the basic read/write routines */
/*******************************************/

/*
 * rfits_() will read a FITS file into head and data.  The data are as
 * stored on the disk, and need to be converted to whatever format is
 * desired by FITSorder() or (e.g.) unsignedfp().
 */

int rfits_(file,nhead,head,bitpix,nx,ny,scale,zero,data,filelen,headlen)
char *file;		/* file name */
int *nhead;		/* number of header entries */
char *head;		/* header */
int *bitpix;		/* BITPIX for disk and resultant array */
double *scale, *zero;	/* BSCALE and BZERO parameters from header */
int *nx, *ny;		/* NAXIS1 and NAXIS2 from header */
char *data;		/* image data */
int headlen, filelen;	/* (Fortran) length of character strings */
{
  int fd, mode=0, err, nbyte;
  if( (err=openc_(&fd, file, &mode, filelen)) != 0) return(err);
  if( (err=rhead_(&fd, nhead, head, headlen)) != 0) {close(fd); return(err);}
  if( (err=parsehead_(nhead,head,bitpix,nx,ny,scale,zero,headlen)) != 0) {
    close(fd);    return(err);  }
  nbyte = (*nx)*(*ny) * (ABS(*bitpix)/8);
  if(*bitpix == 1) nbyte = 2*(((*nx)*(*ny)+15)/16);

  if( (err=read(fd, data, nbyte)) != nbyte) {close(fd); return(ERR_CANT_READ_DATA);}
  close(fd);
  return(0);
}

/*
 * wfits_() will write an integer FITS file onto disk.  It is assumed that
 * the data have been converted to FITS byte order via (say) fpshort or
 * fitsorder(). Legal values for bitpix are 1, 16, and 32.
 */

int wfits_(file,nhead,head,bitpix,nx,ny,data,filelen,headlen)
char *file;		/* file name */
int *bitpix;		/* BITPIX for disk and resultant array */
int *nhead;		/* number of header entries */
char *head;		/* header */
int *nx, *ny;		/* NAXIS1 and NAXIS2 from header */
char *data;		/* image data */
int headlen, filelen;	/* (Fortran) length of character strings */
{
  int fd, mode=1, err, nxh, nyh, bph;
  int nwrite, nbyte, imagebyte, nwrit;
  double scale, zero;

  if( (err=openc_(&fd, file, &mode, filelen)) != 0) return(err);

  if( (err=parsehead_(nhead,head,&bph,&nxh,&nyh,&scale,&zero,headlen)) != 0) {
    close(fd);    return(err);  }
  if( nxh != (*nx) || nyh != (*ny) || bph != (*bitpix) ) {
    close(fd);    return(ERR_INCONSISTENT_NAXIS);  }
  if( (err=whead_(&fd, nhead, head, headlen)) != 0) {
    close(fd);  return(err); }

  imagebyte = (*nx)*(*ny) * (ABS(*bitpix)/8);
  if(*bitpix == 1) imagebyte = 2*(((*nx)*(*ny)+15)/16);
  nbyte = imagebyte;		       	/* Bytes to write */

  nwrit = 0;
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF,nbyte);
    memmove(buf, &data[nwrit], nwrite);
    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
    nbyte -= nwrite;
    nwrit += nwrite;
  }

  nwrite = NFITS - (imagebyte%NFITS);
  if(nwrite < NFITS) {
    bzero(buf,nwrite);
    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
  }

  close(fd);
  return(0);
}

int openc_(fd, file, mode, filelen)
int  *fd;     /* File descriptor */ 
int  *mode;   /* Mode: 0/1 for open/create */
char *file;   /* File name */
int  filelen; /* (Fortran) length of character string file */
{
  int i=filelen-1;
  while(i > 0 && (file[i]==' ' || file[i]==0)) i--;
  memmove(buf, file, i+1);
  buf[i+1] = '\0';
  if(*mode == 1) *fd = creat(buf, 0644);
  else *fd = open(buf, 0);
  if(*fd == -1) {
/*    fprintf(stderr,"Can't open ****%s****\n",buf); */
    return(ERR_CANT_OPEN_FILE);
  }
  return (0);
}

int closec_(fd)
int *fd;
{
  close(*fd);
  return(0);
}

int writebytes_(fd,data,nbyte)
int  *fd;		/* Input file descriptor */
char *data;		/* Output buffer */
int  *nbyte;		/* Length of buffer (bytes) */
{
  if (write (*fd, data, *nbyte) != *nbyte)  {
    fprintf(stderr,"Error writing data...\n"); 
    return (ERR_CANT_WRITE_DATA);
  } 
  else return (0);
}

int readbytes_(fd,data,nbyte)
int  *fd;		/* Input file descriptor */
char *data;		/* Input buffer	*/
int  *nbyte;		/* Length of buffer (bytes) */
{
  int i;
  if ((i=read (*fd, data, *nbyte)) != *nbyte) {
    fprintf(stderr,"%d requested, %d read...\n",*nbyte,i); 
    return (ERR_CANT_READ_DATA);
  } 
  else return (0);
}

/* Read just a little bit to see if the first bytes are "SIMPLE" */
int testfits_(file, mefits, bitpix, filelen)
char *file;
int *mefits, *bitpix, filelen;
{
  int fd, i, nread, mode=0;
  if(openc_(&fd, file, &mode, filelen) != 0) return(ERR_CANT_OPEN_FILE);
  if((nread=read(fd,buf,160)) != 160) {
/*    fprintf(stderr,"testfits: error sniffing header\n"); */
    return(ERR_CANT_READ_HEADER);
  }
  i = strncmp(buf,"SIMPLE",6);
  if(i == 0) {
    *mefits = 1;
    *bitpix = atoi(buf+89);
  } else {
    *mefits = 0;
  }
  close(fd);
  return(0);
}

/* Read just a little bit to see if the first bytes are "SIMPLE" */
int testfitslen_(file, mefits, bitpix, naxis, dims, nhead, filelen)
char *file;
int *mefits, *bitpix, filelen;
int *naxis, *dims, *nhead;
{
  int fd, i, nread, mode=0, n;
  *mefits = 0;
  if(openc_(&fd, file, &mode, filelen) != 0) return(ERR_CANT_OPEN_FILE);
  if((nread=read(fd,buf,160)) != 160) {
/*    fprintf(stderr,"testfits: error sniffing header\n"); */
    close(fd);
    return(ERR_CANT_READ_HEADER);
  }
  i = strncmp(buf,"SIMPLE",6);
  if(i == 0) {
    *mefits = 1;
    *bitpix = atoi(buf+89);
  } else {
    return(0);
  }
/* Now look for NAXIS, NAXISn, and END */
  for(i=0; i<10000; i++) {
    if((nread=read(fd, buf, 80)) != 80) {
      close(fd);
      return(ERR_CANT_READ_HEADER);
    }
/*    printf("%3d %.20s\n", i, buf); */
    if(strncmp(buf, "END     ", 8) == 0) {
      *nhead = i + 3;
      close(fd);
      return(0);
    } else if(strncmp(buf, "NAXIS   ", 8) == 0) {
      *naxis = atoi(buf+10);
    } else if( (strncmp(buf, "NAXIS", 5) == 0) &&
	     (strncmp(buf+6, "  ", 2) == 0) &&
	     (sscanf(buf+5, "%1d", &n) == 1) && n <= *naxis && n > 0) {
      dims[n-1] = atoi(buf+10);
    }
  }
  close(fd);
  return(ERR_NO_END);
}

/* Read just the file's header */
int rhead_(fd,nhead,head,headlen)
int *fd, *nhead, headlen;
char *head;
{
  int i, nread, extra, end;

  end = 0;
  for(i=0;i<(*nhead)*80;i+=80) {
    if((nread=read(*fd,head+i,80)) != 80) {
      fprintf(stderr,"rfhead: error reading header\n"); 
      return(ERR_CANT_READ_HEADER); }
    if(strncmp(&head[i],"END     ",8) == 0) {
      end = 1; break;
    }
  }
  if(end == 0) {
    fprintf(stderr,"rfhead: insufficient buffer for header\n"); 
    return(ERR_INSUFFICIENT_HEADER); }

  i+= 80;
  extra = NFITS - (i%NFITS);
  if((i%NFITS) !=0) {
    if((nread = read(*fd,buf,extra)) != extra) {
      fprintf(stderr,"rfhead: error reading header\n"); 
      return(ERR_CANT_READ_HEADER); }
  }
  *nhead = i / 80;
  return(0);
}

/* Write just the file's header */
int whead_(fd,nhead,head,headlen)
int *fd, *nhead, headlen;
char *head;
{
  int i, nwrite, extra, end;

  end = 0;
  for(i=0; i<(*nhead)*80; i+=80) {

    if((nwrite=write(*fd,head+i,80)) != 80) {
      fprintf(stderr,"wfhead: error writing header\n"); 
      return(ERR_CANT_WRITE_HEADER); }
    if(strncmp(&head[i],"END     ",8) == 0) {
      end = 1; break;
    }
  }

  if(end == 0) {
    fprintf(stderr,"wfhead: END missing from header\n"); 
    return(ERR_NO_END); }
  i+= 80;

  extra = NFITS - (i%NFITS);
  if((i%NFITS) !=0) {
    for(i=0;i<extra;i++) buf[i] = ' ';
    if((nwrite = write(*fd,buf,extra)) != extra) { 
      fprintf(stderr,"wfhead: error writing header\n"); 
      return(ERR_CANT_WRITE_HEADER); }
  }
  return(0);
}

int parsehead_(nhead,head,bitpix,nx,ny,scale,zero,headlen)
int *nhead, *bitpix, *nx, *ny, headlen;
char *head;
double *scale, *zero;
{
  int i, setscale=0, setzero=0;
  double atof();
  *bitpix = *nx = *ny = 0;
  for(i=0;i<80*(*nhead);i+=80) {
    if( strncmp(head+i, "BITPIX", 6) == 0) *bitpix = atoi(head+i+10);
    if( strncmp(head+i, "NAXIS1", 6) == 0) *nx = atoi(head+i+10);
    if( strncmp(head+i, "NAXIS2", 6) == 0) *ny = atoi(head+i+10);
    if( strncmp(head+i, "BSCALE", 6) == 0) {
      *scale = atof(head+i+10);
      setscale += 1; }
    if( strncmp(head+i, "BZERO", 5) == 0)  {
      *zero = atof(head+i+10);
      setzero += 1; }
  }
  if((*nx) == 0 || (*ny) == 0 || (*bitpix) == 0) return(ERR_NO_NAXIS);
  if( setscale != 1 ) *scale = 1.0;
  if( setzero != 1) *zero = 0.0;
  return(0);
}

/**********************************/
/* Basic type conversion routines */
/**********************************/

/*
 * Various format conversion routines
 * Routines converting to FP (or bit to short) work from the end of 
 * the array so the input and output arrays may be identical.
 */

/* Convert bitmap to short */
int bitshort_(npix,data,sh)
int *npix;
unsigned char *data;
short *sh;
{
  int i;
  unsigned char b=0;
/* It seems silly always to SWAB data, even on BIGENDIAN machines, but my
 * current definition of a bitmap is in terms of 16-bit words, with LOW bit
 * first.  It's wrong, but it's the law.
 */
/* #ifdef LOWENDIAN */
  FITSorder(1,*npix,data);
/* #endif */
  if( ((*npix-1)%8) != 7 ) b = data[(*npix-1)/8];
  for(i=(*npix-1);i>=0;i--) {
    if( (i%8) == 7 ) b = data[i/8];
    sh[i] = (b & 0x80) ? 1 : 0;
    b = b << 1;
  }
  return(0);
}

/* Convert bitmap to floating point */
int bitfp_(npix,data,fp)
int *npix;
unsigned char *data;
float *fp;
{
  int i;
  unsigned char b=0;
/* It seems silly always to SWAB data, even on BIGENDIAN machines, but my
 * current definition of a bitmap is in terms of 16-bit words, with LOW bit
 * first.  It's wrong, but it's the law.
 */
/* #ifdef LOWENDIAN */
  FITSorder(1,*npix,data);
/* #endif */
  if( ((*npix-1)%8) != 7 ) b = data[(*npix-1)/8];
  for(i=(*npix-1);i>=0;i--) {
    if( (i%8) == 7 ) b = data[i/8];
    fp[i] = (b & 0x80) ? 1.0 : 0.0;
    b = b << 1;
  }
  return(0);
}

/* Convert short to bitmap */
int shortbit_(npix,sh,data)
int *npix;
unsigned char *data;
short *sh;
{
  unsigned char *dp=data;
  int i;
  *dp = 0;
  for(i=0;i<*npix;i++) {
    *dp = (*dp >> 1) | ((*sh++ == 0) ? 0x00 : 0x80);
    if( (i%8) == 7 ) dp++;
  }
/* It seems silly always to SWAB data, even on BIGENDIAN machines, but my
 * current definition of a bitmap is in terms of 16-bit words, with LOW bit
 * first.  It's wrong, but it's the law.
 */
/* #ifdef LOWENDIAN */
  FITSorder(1,*npix,data);
/* #endif */
  return(0);
}

/* Convert floating point to bitmap */
int fpbit_(npix,fp,data)
int *npix;
unsigned char *data;
float *fp;
{
  unsigned char *dp=data;
  int i;
  *dp = 0;
  for(i=0;i<*npix;i++) {
    *dp = (*dp >> 1) | ((*fp++ == 0.0) ? 0x00 : 0x80);
    if( (i%8) == 7 ) dp++;
  }
/* It seems silly always to SWAB data, even on BIGENDIAN machines, but my
 * current definition of a bitmap is in terms of 16-bit words, with LOW bit
 * first.  It's wrong, but it's the law.
 */
/* #ifdef LOWENDIAN */
  FITSorder(1,*npix,data);
/* #endif */
  return(0);
}

/* Convert short integers to floating point */
int shortfp_(npix,scale,zero,data,fp)
int *npix;
short int *data;
float *fp;
double *scale, *zero;
{
  int i;
  float tol=0.00001*(*scale);
  if(LOWENDIAN()) FITSorder(16,*npix,data);
  if(*scale == 1.0 && *zero == 0.0) {
    for(i=(*npix-1);i>=0;i--) fp[i] = (float)data[i];
  } else {
    for(i=(*npix-1);i>=0;i--) {
      fp[i] = (float)((*zero) + (*scale)*(double)data[i]);
/* If the value is REALLY close to zero, it probably is supposed to be zero */
      if(fp[i] > -tol && fp[i] < tol) fp[i] = 0.0;
    }
  }
  return(0);
}

#define NINT(x) ( ((x) < 0.0) ? (int)((x)-0.5) : (int)((x)+0.5) )
#define TRUNCATE(low,high,x) ( ((x) > (low)) ? (((x) < (high)) ? (x) : (high)) : (low) )

/* Convert floating point to short integers */
int fpshort_(npix,scale,zero,fp,data)
int *npix;
short int *data;
float *fp;
double *scale, *zero;
{
  int i;
  double temp;
  short int *dp=data;
  if(*scale == 1.0 && *zero == 0.0) {
    for(i=0;i<*npix;i++) {
      temp = *fp++;
      temp = TRUNCATE(-32768.5,32767.5,temp);
      *dp++ = NINT(temp);
    }
  } else {
    for(i=0;i<*npix;i++) {
      temp = (*scale) * ((*zero) + (double)(*fp++));
      temp = TRUNCATE(-32768.5,32767.5,temp);
      *dp++ = NINT(temp);
    }
  }
  if(LOWENDIAN()) FITSorder(16,*npix,data);
  return(0);
}

/* Convert unsigned short integers to floating point */
int unsignedfp_(npix,scale,zero,data,fp)
int *npix;
unsigned short int *data;
float *fp;
double *scale, *zero;
{
  int i;
  float tol=0.00001*(*scale);
  if(LOWENDIAN()) FITSorder(16,*npix,data);
  if(*scale == 1.0 && *zero == 0.0) {
    for(i=(*npix-1);i>=0;i--) fp[i] = (float)data[i];
  } else {
    for(i=(*npix-1);i>=0;i--) {
      fp[i] = (float)((*zero) + (*scale)*(double)data[i]);
/* If the value is REALLY close to zero, it probably is supposed to be zero */
      if(fp[i] > -tol && fp[i] < tol) fp[i] = 0.0;
    }
  }
  return(0);
}

/* Convert floating point to unsigned short integers */
int fpunsigned_(npix,scale,zero,fp,data)
int *npix;
unsigned short int *data;
float *fp;
double *scale, *zero;
{
  int i;
  double temp;
  unsigned short int *dp=data;
  if(*scale == 1.0 && *zero == 0.0) {
    for(i=0;i<*npix;i++) {
      temp = *fp++;
      temp = TRUNCATE(-0.5,65535.5,temp);
      *dp++ = NINT(temp);
    }
  } else {
    for(i=0;i<*npix;i++) {
      temp = (*scale) * ((*zero) + (*fp++));
      temp = TRUNCATE(-0.5,65535.5,temp);
      *dp++ = NINT(temp);
    }
  }
  if(LOWENDIAN()) FITSorder(16,*npix,data);
  return(0);
}

/* Convert long integers to floating point */
int longfp_(npix,scale,zero,data,fp)
int *npix;
long int *data;
float *fp;
double *scale, *zero;
{
  float tol=0.00001*(*scale);
  int i;
  if(LOWENDIAN()) FITSorder(32,*npix,data);
  if(*scale == 1.0 && *zero == 0.0) {
    for(i=(*npix-1);i>=0;i--) fp[i] = (float)data[i];
  } else {
    for(i=(*npix-1);i>=0;i--) {
      fp[i] = (*zero) + (*scale)*data[i];
/* If the value is REALLY close to zero, it probably is supposed to be zero */
      if(fp[i] > -tol && fp[i] < tol) fp[i] = 0.0;
    }
  }
  return(0);
}

/* Convert floating point to long integers */
int fplong_(npix,scale,zero,fp,data)
int *npix;
long int *data;
float *fp;
double *scale, *zero;
{
  int i;
  double temp;
  long int *dp=data;
  if(*scale == 1.0 && *zero == 0.0) {
    for(i=0;i<*npix;i++) {
      temp = *fp++;
      temp = TRUNCATE(-2147483648.5,2147483647.5,temp);
      *dp++ = NINT(temp);
    }
  } else {
    for(i=0;i<*npix;i++) {
      temp = (*scale) * ((*zero) + (double)(*fp++));
      temp = TRUNCATE(-2147483648.5,2147483647.5,temp);
      *dp++ = NINT(temp);
    }
  }
  if(LOWENDIAN()) FITSorder(32,*npix,data);
  return(0);
}

/* Convert IEEE format floating point to floating point */
int ieeefp_(npix,data,fp)
int *npix;
char *data;
char *fp;		/* Really FP, but may need to access bytes */
{
#ifdef VAXFP
  register char *dp=data, temp;
  for(i=0;i<*npix;i++) {
    temp  = (*dp == 0) ? 0 : (*dp+1) ;
    *fp++ = *(dp+1);
    *fp++ = temp;
    temp  = *(dp+2);
    *fp++ = *(dp+3);
    *fp++ = temp;
    dp += 4;
  }
#else
  if(LOWENDIAN()) FITSorder(32,*npix,data);
  memmove(fp, data, 4*(*npix));
#endif
  return(0);
}

/* Convert floating point to IEEE format floating point */
int fpieee_(npix,fp,data)
int *npix;
char *data;
char *fp;		/* Really FP, but may need to access bytes */
{
#ifdef VAXFP
  register char *dp=data, temp;
  for(i=0;i<*npix;i++) {
    temp  = *(fp);
    *dp++ = (*(fp+1) == 0) ? 0 : (*(fp+1)-1);
    *dp++ = temp;
    temp  = *(fp+2);
    *dp++ = *(fp+3);
    *dp++ = temp;
    fp += 4;
  }
#else
  memmove(data, fp, 4*(*npix));
  if(LOWENDIAN()) FITSorder(32,*npix,data);
#endif
  return(0);
}

/* Swap bytes if bitpix <= 16; reflect 4 bytes if bitpix = +/-32 */
int fitsorder_(bitpix,npix,data)
int *bitpix, *npix;
short *data;
{
  if(LOWENDIAN()) FITSorder(*bitpix,*npix,data);
  return(0);
}

#ifdef OSF_ALPHA
int swab(src, dest, n)
char *src, *dest;
int n;
{
  register char temp, *s=src, *d=dest;
  while(n>1) {
    temp = *s++;
    *d++ = *s++;
    *d++ = temp;
    n -= 2;
  }
  return(0);
}
#endif

int FITSorder(bitpix,npix,data)
int bitpix, npix;
short *data;
{
  int nbyte;
  register short temp, *dp=data;
  register int i=npix;

/* Swap bytes first; work on an even number of bytes */

  nbyte = npix * (ABS(bitpix)/8);
  if(bitpix == 1) nbyte = 2*((npix+15)/16);

  swab(data,data, nbyte);
  if(ABS(bitpix) <= 16) return(0);

/* Swap words to complete the 32 bit byte swap */
  while(i--) {
    temp = *dp;
    *dp = *(dp+1);
    dp++;
    *dp++ = temp;
  }
  return(0);
}
