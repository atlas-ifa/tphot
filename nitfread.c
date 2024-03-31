/* Read a NITF file, format metadata into FITS header */
/* Derived from definition of NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.1 in 
 * MIL-STD-2500B and updates on extensions */
/* v1.1 120707 bust out as a library */
/* v1.0 120704 John Tonry */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/file.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include "nitf.h"

int nitfread(char *nitfile, int bitpix, 
	     char **fitshead, int *nx, int *ny, float **data, int verbose)
{
   FILE *fp;
   int i, j, k;
   int nbx, nby, nxblk, nyblk, nbpp;
   void *nitfdata;

/* Open the NITF file */
   if( (fp = fopen(nitfile, "r")) == NULL) {
      fprintf(stderr, "Cannot open %s for input\n", nitfile);
      return(-1);
   }

/* Read all the primary header quantitities */
   if(verbose > 0) printf("=File Header=\n");
   for(i=0; i<nitf_file_len; i++) {
      readfield(nitf_file_head+i, fp, verbose);
   }

/* Read all the data segment sub-header quantitities */
   if(verbose > 0) printf("\n=Image Header=\n");
   for(i=0; i<nitf_image_len; i++) {
      readfield(nitf_image_head+i, fp, verbose);
   }

/* Get the data size parameters */
   image_size(nx, ny, &nbx, &nby, &nxblk, &nyblk, &nbpp);

   if(verbose > 0) {
      printf("%d %d %d %d %d %d %d\n", *nx, *ny, nbx, nby, nxblk, nyblk, nbpp);
   }

   mkfitshead(fitshead, bitpix, *nx, *ny);

/* Read the data or just stop at the header? */
   if(data == NULL) return(0);

/* Read and reformat into a contiguous array */
   nitfdata = (void *)calloc(nxblk*nyblk, nbpp/8);
   *data = (void *)calloc((*nx)*(*ny), ABS(bitpix)/8);
   i = idxfield(nitf_image_len, nitf_image_head, "PVTYPE");
   if(nitf_image_head[i].val[0] == 'R') nbpp = -nbpp;

   for(j=0; j<nby; j++) {
      for(i=0; i<nbx; i++) {
	 if(verbose > 1) printf("%d %d %d %d\n", i, j, i*nxblk, (*ny)-1-j*nyblk);
	 k = fread(nitfdata, ABS(nbpp)/8, nxblk*nyblk, fp);
	 if( (k=reform(nxblk, nyblk, nbpp, nitfdata, i*nxblk, (*ny)-1-j*nyblk,
		       *nx, *ny, bitpix, *data, 1))) return(k);
      }
   }
   free(nitfdata);

   return(0);
}

/* Parse the nitf_image structure for data size parameters */
void image_size(int *nx, int *ny, int *nbx, int *nby, 
		int *nxblk, int *nyblk, int *nbpp)
{
   int i;
   i = idxfield(nitf_image_len, nitf_image_head, "NCOLS");
   *nx = atoi(nitf_image_head[i].val);	// Number of columns
   i = idxfield(nitf_image_len, nitf_image_head, "NROWS");
   *ny = atoi(nitf_image_head[i].val);	// Number of rows

   i = idxfield(nitf_image_len, nitf_image_head, "NBPR");
   *nbx = atoi(nitf_image_head[i].val);	// Number of blocks per row
   i = idxfield(nitf_image_len, nitf_image_head, "NBPC");
   *nby = atoi(nitf_image_head[i].val);	// Number of blocks per column

   i = idxfield(nitf_image_len, nitf_image_head, "NPPBH");
   *nxblk = atoi(nitf_image_head[i].val);// Pixels per block horizontally
   i = idxfield(nitf_image_len, nitf_image_head, "NPPBV");
   *nyblk = atoi(nitf_image_head[i].val);// Pixels per block vertically

   i = idxfield(nitf_image_len, nitf_image_head, "NBPP");
   *nbpp = atoi(nitf_image_head[i].val);
}

/* Return the index of a field in a keyval array */
int idxfield(int n, NITF_KEYVAL *kv, char *field)
{
   int i;
   for(i=0; i<n; i++) {
      if(strcmp(kv[i].field, field) == 0) return(i);
   }
   return(-1);
}

/* Read TRE header data, saving the ones that we know about, skipping others */
/* fp should be pointed at the start of the TRE(s) */
int readtre(int nxhd, FILE *fp, int verbose)
{
   int i, j, k, ntre;
   char line[256];
   for(j=0; j<nxhd-11; ) {
      k = fread(line, sizeof(char), 6+5, fp);
      if(verbose > 0) printf("\n=TRE= %6.6s\n", line);
      fseek(fp, (long)(-6-5), SEEK_CUR);

      ntre = atoi(line+6) + 6+5;

      if(strncmp(line, "STDIDC", 6) == 0) {
	 for(i=0; i<stdidc_len; i++) {
	    readfield(stdidc_tre+i, fp, verbose);
	 }

      } else if(strncmp(line, "USE00A", 6) == 0) {
	 for(i=0; i<use00a_len; i++) {
	    readfield(use00a_tre+i, fp, verbose);
	 }

      } else if(strncmp(line, "RPC00", 5) == 0) {
	 for(i=0; i<rpc00b_len; i++) {
	    readfield(rpc00b_tre+i, fp, verbose);
	 }

      } else {
	 fseek(fp, (long)(ntre), SEEK_CUR);

      }
      j += ntre;
   }
   return(0);
}

/* Read a field's value from fp */
int readfield(NITF_KEYVAL *kv, FILE *fp, int verb)
{
   int i, l=kv->flen, n, l1, l2, ntre, novfl;
   char line[1024], fmt[1024], type;
   kv->val = malloc(l+1);
   if(fread(kv->val, sizeof(char), l, fp) != l) return(-1);
   kv->val[l] = '\0';

/* Evaluate the value? */
   type = kv->fmt[strlen(kv->fmt)-1];
   if(type != 's') {
      if(type == 'd') {
	 strcpy(fmt, "%d");
      } else {
	 strcpy(fmt, "%lf");
      }
      kv->data = (void *)calloc(1, sizeof(double));
      sscanf(kv->val, fmt, kv->data);
   }

/* Are header/data length subfields called for? */
   if(kv->subflen > 0) {
      n = *(int *)(kv->data);
      l1 = kv->subflen % FLBASE;
      l2 = (kv->subflen) / FLBASE;
      kv->sbhlen = (int *)calloc(n, sizeof(int));
      kv->seglen = (int *)calloc(n, sizeof(int));
      for(i=0; i<n; i++) {
	 if(fread(line, sizeof(char), l1, fp) != l1) return(-1);
	 line[l1] = '\0';
	 kv->sbhlen[i] = atoi(line);
	 if(fread(line, sizeof(char), l2, fp) != l2) return(-1);
	 line[l2] = '\0';
	 kv->seglen[i] = atoi(line);
      }

/* Is this an overflow set of TREs? */
   } else if(kv->subflen < 0 && (ntre = atoi(kv->val)) > 0) {

/* (Ignore) novfl and the DES possibility of an overflow, ugh */
      if(fread(line, sizeof(char),-kv->subflen,fp) != -kv->subflen) return(-1);
      line[-kv->subflen] = '\0';
      novfl = atoi(line);

      readtre(ntre-3, fp, verb);
   }

/* Tell us about it? */
   if(verb > 0) {
      if(strcmp(kv->units, "text") == 0 || kv->val[0] == ' ') {
	 printf("%-18s= '%s' (%d)", kv->field, kv->val, kv->flen);
      } else {
	 printf("%-18s= %s (%d)", kv->field, kv->val, kv->flen);
      }
      if(kv->data == NULL || kv->subflen <= 0 || *(int *)kv->data == 0) {
	 printf("\n");
      } else {
	 printf(" %3d", *(int *)kv->data);
	 for(i=0; i<*(int *)kv->data; i++) {
	    printf(" %6d %10d", kv->sbhlen[i], kv->seglen[i]);
	 }
	 printf("\n");
      }
   }
   return(0);
}

/* Reformat one block of data into another */
int reform(int mx, int my, int mbpp, void *mdata, int ix, int iy, 
	   int nx, int ny, int bitpix, void *data, int yflip)
{
   int j, l, k1, k2;

#if 0
   v = *((unsigned char *)(mdata+512+512*mx));
   printf("%4d %4d %3d   %5d %5d   %5d %5d %3d  %4d\n", mx, my, mbpp, ix, iy, 
	  nx, ny, bitpix, v);
#endif

   for(j=0; j<my; j++) {
      l = yflip ? iy-j :  iy+j;
      if(l > ny-1 || l < 0) continue;
      k1 = MAX(0, ix);
      k2 = MIN(nx-1, ix+mx-1);

      if(mbpp == 8 && bitpix == 16) {
	 charshort(k2-k1+1, ((unsigned char *)mdata)+k1-ix+j*mx,
		    ((unsigned short *)data)+k1+l*nx);
      } else if(mbpp == 16 && bitpix == 16) {
	 memcpy(((unsigned short *)data)+k1+l*nx,
		((unsigned short *)mdata)+k1-ix+j*mx, sizeof(short int)*(k2-k1+1));
	 FITSorder(bitpix, k2-k1+1, ((unsigned short int *)data)+k1+l*nx);
      } else if(mbpp == 32 && bitpix == 16) {
	 longshort(k2-k1+1, ((int *)mdata)+k1-ix+j*mx,
		    ((unsigned short *)data)+k1+l*nx);
      } else if(mbpp == -32 && bitpix == -32) {
	 memcpy(((float *)data)+k1+l*nx,
		((float *)mdata)+k1-ix+j*mx, sizeof(float)*(k2-k1+1));
	 FITSorder(bitpix, k2-k1+1, ((float *)data)+k1+l*nx);
      } else if(mbpp == 8 && bitpix == -32) {
	 charfloat(k2-k1+1, ((unsigned char *)mdata)+k1-ix+j*mx,
		    ((float *)data)+k1+l*nx);
      } else {
	 fprintf(stderr, "Error: reform cannot reformat bpp %d to %d\n", 
		 mbpp, bitpix);
	 return(-2);
      }
   }
   return(0);
}

/* addfitshead() adds a new line to a FITS header */
int addfitshead(int n, char *header, char *keyword, 
	    char *cvalue, int ivalue, double rvalue);

int stuffield(int *nhead, char *head, int n, NITF_KEYVAL *kv, 
	      char *nitfkey, char *fitskey)
{
   int i;
   i = idxfield(n, kv, nitfkey);
   if(i < 0) return(-1);
//   printf("%d %s\n", i, nitf_file_head[i].val);
   wfchar(*nhead, head, fitskey, kv[i].val);
   (*nhead) += 1;
   return(0);
}

/* wfchar writes a specific line to a FITS header */
int wfchar(int n, char *header, char *keyword, char *cvalue)
{
   int i, l;
   char fmt[80];
   for(l=MIN(79, strlen(cvalue)); l>0; l--) if(cvalue[l-1] != ' ') break;
   sprintf(fmt, "%%-8s= '%%-%d.%ds'", l, l);
//   printf("%d %s %s %s\n", l, fmt, keyword, cvalue);
   i = sprintf(&header[80*n], fmt, keyword, cvalue);
   if (i < 0) {
      fprintf(stderr,"wfitem: Error writing header string '%s'!\n", keyword);
      return(1);
   }
   for( ; i<80; i++) if (header[80*n+i] == '\0') header[80*n+i] = ' ';
   return(0);
}

typedef struct {
   int lngd;
   int lngm;
   int lngs;
   char ew;
   int latd;
   int latm;
   int lats;
   char ns;
   double lng;
   double lat;
} LNGLAT;

int stuffgeo(int *nhead, char *head, int nx, int ny)
{
   int i, n;
   LNGLAT UL, UR, LR, LL;
   i = idxfield(nitf_image_len, nitf_image_head, "IGEOLO");
   if(i > 0) {
      n = sscanf(nitf_image_head[i].val, "%2d%2d%2d%c%3d%2d%2d%c%2d%2d%2d%c%3d%2d%2d%c%2d%2d%2d%c%3d%2d%2d%c%2d%2d%2d%c%3d%2d%2d%c",
	     &UL.latd, &UL.latm, &UL.lats, &UL.ns, &UL.lngd, &UL.lngm, &UL.lngs, &UL.ew,
	     &UR.latd, &UR.latm, &UR.lats, &UR.ns, &UR.lngd, &UR.lngm, &UR.lngs, &UR.ew,
	     &LR.latd, &LR.latm, &LR.lats, &LR.ns, &LR.lngd, &LR.lngm, &LR.lngs, &LR.ew,
	     &LL.latd, &LL.latm, &LL.lats, &LL.ns, &LL.lngd, &LL.lngm, &LL.lngs, &LL.ew);
      if(n != 32) {
	 printf("Cannot scan IGEOLO from %s\n", nitf_image_head[i].val);
	 return(-1);
      } else {
	 UL.lat = UL.latd + UL.latm/60.0 + UL.latm/3600.0;
	 UL.lat *= UL.ns == 'N' ? +1 : -1;
	 UL.lng = UL.lngd + UL.lngm/60.0 + UL.lngm/3600.0;
	 UL.lng *= UL.ew == 'E' ? +1 : -1;
	 UR.lat = UR.latd + UR.latm/60.0 + UR.latm/3600.0;
	 UR.lat *= UR.ns == 'N' ? +1 : -1;
	 UR.lng = UR.lngd + UR.lngm/60.0 + UR.lngm/3600.0;
	 UR.lng *= UR.ew == 'E' ? +1 : -1;
	 LR.lat = LR.latd + LR.latm/60.0 + LR.latm/3600.0;
	 LR.lat *= LR.ns == 'N' ? +1 : -1;
	 LR.lng = LR.lngd + LR.lngm/60.0 + LR.lngm/3600.0;
	 LR.lng *= LR.ew == 'E' ? +1 : -1;
	 LL.lat = LL.latd + LL.latm/60.0 + LL.latm/3600.0;
	 LL.lat *= LL.ns == 'N' ? +1 : -1;
	 LL.lng = LL.lngd + LL.lngm/60.0 + LL.lngm/3600.0;
	 LL.lng *= LL.ew == 'E' ? +1 : -1;
	 wfitem(*nhead, head, "CRPIX1  ", NULL, 0, 0.0);
	 (*nhead) += 1;
	 wfitem(*nhead, head, "CRPIX2  ", NULL, 0, 0.0);
	 (*nhead) += 1;
	 wfitem(*nhead, head, "CRVAL1  ", "FLOAT", 0, LL.lng);
	 (*nhead) += 1;
	 wfitem(*nhead, head, "CRVAL2  ", "FLOAT", 0, LL.lat);
	 (*nhead) += 1;
	 wfitem(*nhead, head, "CDELT1  ", "FLOAT", 0, (LR.lng-LL.lng)/nx);
	 (*nhead) += 1;
	 wfitem(*nhead, head, "CDELT2  ", "FLOAT", 0, (UL.lat-LL.lat)/ny);
	 (*nhead) += 1;
      }
   }
   return(0);
}

char *rpct[20]={
   "1", "a", "d", "h", "a*d", "a*h", "d*h", "a*a", "d*d", "h*h",
   "a*d*h", "a*a*a", "a*d*d", "a*h*h", "a*a*d", "d*d*d", "d*h*h", "a*a*h",
   "d*d*h", "h*h*h"};

/* Fill RPC coefficients from NITF TRE into FITS header */
int stuffrpc(int *nhead, char *head)
{
   int rpca, i, k;
   char name[256];
   double bias[2], off[5], scl[5], rn[20], rd[20], cn[20], cd[20];
   if(rpc00b_tre[0].val == NULL) return(-1);
   rpca = strncmp(rpc00b_tre[0].val, "RPC00A", 6) == 0;

/* Get the coefficients */
   for(i=0; i<20; i++) {
      sprintf(name, "LINE_NUM_COEFF_%d", i+1);
      if( (k=idxfield(rpc00b_len, rpc00b_tre, name)) < 0) return(-1);
      rn[i] = *(double *)(rpc00b_tre[k].data);

      sprintf(name, "LINE_DEN_COEFF_%d", i+1);
      if( (k=idxfield(rpc00b_len, rpc00b_tre, name)) < 0) return(-1);
      rd[i] = *(double *)(rpc00b_tre[k].data);

      sprintf(name, "SAMP_NUM_COEFF_%d", i+1);
      if( (k=idxfield(rpc00b_len, rpc00b_tre, name)) < 0) return(-1);
      cn[i] = *(double *)(rpc00b_tre[k].data);

      sprintf(name, "SAMP_DEN_COEFF_%d", i+1);
      if( (k=idxfield(rpc00b_len, rpc00b_tre, name)) < 0) return(-1);
      cd[i] = *(double *)(rpc00b_tre[k].data);
   }

/* Systematic and random error */
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "ERR_BIAS")) < 0) return(-1);
   bias[0] = *(double *)(rpc00b_tre[k].data);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "ERR_RAND")) < 0) return(-1);
   bias[1] = *(double *)(rpc00b_tre[k].data);

/* Offsets */
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "LINE_OFF")) < 0) return(-1);
   off[0] = atoi(rpc00b_tre[k].val);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "SAMP_OFF")) < 0) return(-1);
   off[1] = atoi(rpc00b_tre[k].val);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "LONG_OFF")) < 0) return(-1);
   off[2] = *(double *)(rpc00b_tre[k].data);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "LAT_OFF")) < 0) return(-1);
   off[3] = *(double *)(rpc00b_tre[k].data);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "HEIGHT_OFF")) < 0) return(-1);
   off[4] = atoi(rpc00b_tre[k].val);

/* Scales */
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "LINE_SCALE")) < 0) return(-1);
   scl[0] = atoi(rpc00b_tre[k].val);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "SAMP_SCALE")) < 0) return(-1);
   scl[1] = atoi(rpc00b_tre[k].val);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "LONG_SCALE")) < 0) return(-1);
   scl[2] = *(double *)(rpc00b_tre[k].data);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "LAT_SCALE")) < 0) return(-1);
   scl[3] = *(double *)(rpc00b_tre[k].data);
   if( (k=idxfield(rpc00b_len, rpc00b_tre, "HEIGHT_SCALE")) < 0) return(-1);
   scl[4] = atoi(rpc00b_tre[k].val);

/* Convert to RPC00B if necessary */
   if(rpca) {
      permute_rpcA2B(rn);
      permute_rpcA2B(rd);
      permute_rpcA2B(cn);
      permute_rpcA2B(cd);
   }

/* Make the FITS cards */
   sprintf(head+(*nhead+0)*80, 
	   "RPC_C0  = %20.1f / RPC sample offset x=(col-C0)/C1", off[1]);
   sprintf(head+(*nhead+1)*80, 
	   "RPC_C1  = %20.1f / RPC sample scale", scl[1]);
   sprintf(head+(*nhead+2)*80, 
	   "RPC_R0  = %20.1f / RPC line offset y=(row-R0)/R1", off[0]);
   sprintf(head+(*nhead+3)*80, 
	   "RPC_R1  = %20.1f / RPC line scale", scl[0]);
   sprintf(head+(*nhead+4)*80, 
	   "RPC_A0  = %20.6f / RPC longitude offset a=(lng-A0)/A1", off[2]);
   sprintf(head+(*nhead+5)*80, 
	   "RPC_A1  = %20.6f / RPC longitude scale", scl[2]);
   sprintf(head+(*nhead+6)*80, 
	   "RPC_D0  = %20.6f / RPC latitude offset d=(lat-D0)/D1", off[3]);
   sprintf(head+(*nhead+7)*80, 
	   "RPC_D1  = %20.6f / RPC latitude scale", scl[3]);
   sprintf(head+(*nhead+8)*80, 
	   "RPC_H0  = %20.1f / RPC height offset h=(lat-H0)/H1", off[4]);
   sprintf(head+(*nhead+9)*80, 
	   "RPC_H1  = %20.1f / RPC height scale", scl[4]);
   (*nhead) += 10;

   sprintf(head+(*nhead+0)*80, 
	   "RPC_BIAS= %20.6f / RPC bias error", bias[0]);
   sprintf(head+(*nhead+1)*80, 
	   "RPC_RAND= %20.6f / RPC random error", bias[1]);
   (*nhead) += 2;

   for(i=0; i<20; i++) {
      sprintf(head+(*nhead+i)*80, "RPC_CN%02d= %20.6e / RPC col num %02d %s", 
	      i, cn[i], i, rpct[i]);
   }
   (*nhead) += 20;
   
   for(i=0; i<20; i++) {
      sprintf(head+(*nhead+i)*80, "RPC_CD%02d= %20.6e / RPC col den %02d %s", 
	      i, cd[i], i, rpct[i]);
   }
   (*nhead) += 20;
   
   for(i=0; i<20; i++) {
      sprintf(head+(*nhead+i)*80, "RPC_RN%02d= %20.6e / RPC row num %02d %s", 
	      i, rn[i], i, rpct[i]);
   }
   (*nhead) += 20;
   
   for(i=0; i<20; i++) {
      sprintf(head+(*nhead+i)*80, "RPC_RD%02d= %20.6e / RPC row den %02d %s", 
	      i, rd[i], i, rpct[i]);
   }
   (*nhead) += 20;
   
   return(0);
}

#define FITSHEADBUF 5

/* Build a FITS header from NITF quantities */
void mkfitshead(char **head, int bitpix, int nx, int ny)
{
   int i, nhead=0, yr=0,mo=0,dy=0,h=0,m=0,s=0;
   double mjd=0.0;
   *head = (char *) malloc(FITSHEADBUF*NFITS+1);
   for(i=0; i<FITSHEADBUF*NFITS; i++) (*head)[i] = ' ';
   (*head)[FITSHEADBUF*NFITS] = '\0';	/* Just to keep a strlen() happy */
   wfitem(nhead++, *head, "SIMPLE  ", "T", 0, 0.0);
   wfitem(nhead++, *head, "BITPIX  ", NULL, bitpix, 0.0);
   wfitem(nhead++, *head, "NAXIS   ", NULL, 2, 0.0);
   wfitem(nhead++, *head, "NAXIS1  ", NULL, nx, 0.0);
   wfitem(nhead++, *head, "NAXIS2  ", NULL, ny, 0.0);
   wfitem(nhead++, *head, "CNPIX1  ", NULL, 0, 0.0);
   wfitem(nhead++, *head, "CNPIX2  ", NULL, 0, 0.0);
   stuffield(&nhead, *head, nitf_file_len,nitf_file_head, "FTITLE", "OBJECT  ");

   stuffield(&nhead, *head, nitf_image_len,nitf_image_head,"IDATIM","IDATIM  ");
   stuffield(&nhead, *head, nitf_image_len,nitf_image_head,"ICAT","IMGCAT  ");
   stuffield(&nhead, *head, nitf_image_len,nitf_image_head,"IID1","IMGID   ");

   stuffgeo(&nhead, *head, nx, ny);

   i = idxfield(nitf_image_len, nitf_image_head, "IDATIM");
   if(sscanf(nitf_image_head[i].val, "%4d%2d%2d%2d%2d%2d", 
	     &yr,&mo,&dy, &h,&m,&s) == 6) {
      mjd = jdate(yr, mo, dy, (double)(h+m/60.0+s/3600.0)) - 2400000.5;
   } else if(sscanf(nitf_image_head[i].val, "%4d%2d%2d", 
	     &yr,&mo,&dy) == 3) {
      mjd = jdate(yr, mo, dy, 0.0) - 2400000.5;
   }
   sprintf(*head+(nhead++)*80, 
	   "DATE-OBS= '%04d-%02d-%02d'         / Observation UT date",yr,mo,dy);
   sprintf(*head+(nhead++)*80, 
	   "UT      = '%02d:%02d:%02d'           / Observation UT", h,m,s);
   sprintf(*head+(nhead++)*80, 
	   "MJD     = %20.5f / Observation modified Julian date", mjd);

   stuffrpc(&nhead, *head);

   wfitem(nhead++, *head, "END     ", NULL, 0, 0.0);
   for(i=0; i<FITSHEADBUF*NFITS; i++) if((*head)[i] == '\0') (*head)[i] = ' ';
}

/* Convert a string of bytes to short ints */
void charshort(int n, unsigned char *in, unsigned short *out)
{
   register unsigned char *i=in+n-1;
   register unsigned short *o=out+n-1;
   while(n--) *o-- = *i--;
}

/* Convert a string of bytes to floats */
void charfloat(int n, unsigned char *in, float *out)
{
   register unsigned char *i=in+n-1;
   register float *o=out+n-1;
   while(n--) *o-- = *i--;
}

/* Convert a string of ints to short ints */
void longshort(int n, int *in, unsigned short *out)
{
   register int *i=in+n-1, npix=n;
   register unsigned short *o=out+n-1;
   while(npix--) *o-- = *i--;
   FITSorder(16, n, out);
}

/* Convert a string of floats to short ints */
void floatshort(int n, float *in, unsigned short *out)
{
   register float *i=in+n-1;
   register int npix=n;
   register unsigned short *o=out+n-1;
   while(npix--) *o-- = *i--;
   FITSorder(16, n, out);
}

/* Convert an RPC00A into RPC00B */
int permute_rpcA2B(double *rpc)
{
   int i, A2B[]={0,1,2,3,4,5,6,10,7,8,9,11,14,17,12,15,18,13,16,19};
   double tmp[20];
   for(i=0; i<20; i++) tmp[A2B[i]] = rpc[i];
   for(i=0; i<20; i++) rpc[i] = tmp[i];
   return(0);
}

double jdate(int year, int month, int day, double ut)
{
   double tmp;
   tmp = 367.0*year + 0.5 + ut/24.0;
   tmp = tmp - ((7*(year+((month+9)/12)))/4) + ((275*month)/9) + day + 1721013;
   return(tmp);
}
