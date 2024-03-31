/* Header for NITF read/write routines */
/* Derived from definition of NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.1 in 
 * MIL-STD-2500B and updates on extensions */
/* v1.0 120704 John Tonry */

#define ABS(a) (((a) > 0) ? (a) : -(a))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define NINT(x) ((x)<0?(int)((x)-0.5):(int)((x)+0.5))
#define FRAC(x) ((x) - (int)(x))

// #define TEST			/* Lots of test output? */

#define NFITS (2880)

#define FLBASE 10

typedef struct NITF_KEYVEL_ELEM {
   char *field;			/* Name of field */
   int flen;			/* Size of field */
   int subflen;			/* Size of subfields: hdrlen+FLBASE*seglen */
   char *fmt;			/* Read/write format */
   char *units;			/* Units */
   char *descrip;		/* Description of field */
   char *val;			/* Character value */
   void *data;			/* Values */
   int *sbhlen;			/* segment subheader lengths */
   int *seglen;			/* segment data lengths */
} NITF_KEYVAL;

/* One NITF 2.1 file header */
NITF_KEYVAL nitf_file_head[]={
   {"FHDR",    4, 0, "%-4.4s",   "text", "File profile name", NULL, NULL, NULL, NULL},
   {"FVER",    5, 0, "%-5.5s",   "text", "File profile name", NULL, NULL, NULL, NULL},
   {"CLEVEL",  2, 0, "%02d",     "int",  "Complexity level", NULL, NULL, NULL, NULL},
   {"STYPE",   4, 0, "%-4.4s",   "text", "Standard type", NULL, NULL, NULL, NULL},
   {"OSTAID", 10, 0, "%-10.10s", "text", "Originating station ID", NULL, NULL, NULL, NULL},
   {"FDT",    14, 0, "%-14.14s", "text", "File date and time CCYYMMDDhhmmss", NULL, NULL, NULL, NULL},
   {"FTITLE", 80, 0, "%-80.80s", "text", "File title", NULL, NULL, NULL, NULL},
   {"FSCLAS",  1, 0, "%-1.1s",   "text", "Security classification", NULL, NULL, NULL, NULL},
   {"FSCLSY",  2, 0, "%-2.2s",   "text", "Security classification scheme", NULL, NULL, NULL, NULL},
   {"FSCODE", 11, 0, "%-11.11s", "text", "File codewords", NULL, NULL, NULL, NULL},
   {"FSCTLH",  2, 0, "%-2.2s",   "text", "File control and handling", NULL, NULL, NULL, NULL},
   {"FSREL",  20, 0, "%-20.20s", "text", "File releasing instructions", NULL, NULL, NULL, NULL},
   {"FSDCTP",  2, 0, "%-2.2s",   "text", "File declass type", NULL, NULL, NULL, NULL},
   {"FSDCDT",  8, 0, "%-8.8s",   "text", "File declass date", NULL, NULL, NULL, NULL},
   {"FSDCXM",  4, 0, "%-4.4s",   "text", "File declass exemption", NULL, NULL, NULL, NULL},
   {"FSDG",    1, 0, "%-1.1s",   "text", "File downgrade", NULL, NULL, NULL, NULL},
   {"FSDGDT",  8, 0, "%-8.8s",   "text", "File downgrade date", NULL, NULL, NULL, NULL},
   {"FSCLTX", 43, 0, "%-43.43s", "text", "Classification text", NULL, NULL, NULL, NULL},
   {"FSCATP",  1, 0, "%-1.1s",   "text", "File class authority type", NULL, NULL, NULL, NULL},
   {"FSCAUT", 40, 0, "%-40.40s", "text", "Classification authority", NULL, NULL, NULL, NULL},
   {"FSCRSN",  1, 0, "%-1.1s",   "text", "File classification reason", NULL, NULL, NULL, NULL},
   {"FSSRDT",  8, 0, "%-8.8s",   "text", "File security source date", NULL, NULL, NULL, NULL},
   {"FSCTLN", 15, 0, "%-15.15s", "text", "File security control number", NULL, NULL, NULL, NULL},
   {"FSCOP",   5, 0, "%05d",     "int",  "File copy number", NULL, NULL, NULL, NULL},
   {"FSCPYS",  5, 0, "%05d",     "text", "File number of copies", NULL, NULL, NULL, NULL},
   {"ENCRYP",  1, 0, "%01d",     "text", "Encryption?", NULL, NULL, NULL, NULL},
   {"FBKGC",   3, 0, "%-3.3s",   "text", "File background color", NULL, NULL, NULL, NULL},
   {"ONAME",  24, 0, "%-24.24s", "text", "Originator's name", NULL, NULL, NULL, NULL},
   {"OPHONE", 18, 0, "%-18.18s", "text", "Originator's phone", NULL, NULL, NULL, NULL},
   {"FL",     12, 0, "%012d",    "int",  "File length", NULL, NULL, NULL, NULL},
   {"HL",      6, 0, "%06d",     "int",  "NITF header length", NULL, NULL, NULL, NULL},
   {"NUMI",    3,106, "%03d",     "int",  "Number of images NUMIxLI...", NULL, NULL, NULL, NULL},
   {"NUMS",    3, 64, "%03d",     "int",  "Number of graphics NUMSxLS...", NULL, NULL, NULL, NULL},
   {"NUMX",    3, 0, "%03d",     "int",  "Reserved for future", NULL, NULL, NULL, NULL},
   {"NUMT",    3, 54, "%03d",     "int",  "Number of text files NUMTxLT...", NULL, NULL, NULL, NULL},
   {"NUMDES",  3, 94, "%03d",     "int",  "Number of data extension NUMD*LD...", NULL, NULL, NULL, NULL},
   {"NUMRES",  3, 74, "%03d",     "int",  "Number of reserved extension NUMR*LR...", NULL, NULL, NULL, NULL},
   {"UDHDL",   5, -3, "%05d",     "int",  "User defined header length", NULL, NULL, NULL, NULL},
// If UDHDL > 0
//   {"UDHD",   N, 0, "%03d",     "int",  "User defined header data", NULL, NULL, NULL, NULL};
   {"XHDL",    5, 0, "%05d",     "int",  "Extended header data length", NULL, NULL, NULL, NULL}
};
static int nitf_file_len=(sizeof(nitf_file_head)/sizeof(NITF_KEYVAL));


NITF_KEYVAL nitf_image_head[]={
   {"IM",          2, 0, "%-2.2s",   "text", "File part type", NULL, NULL, NULL, NULL},
   {"IID1",       10, 0, "%-10.10s", "text", "Image ID1", NULL, NULL, NULL, NULL},
   {"IDATIM",     14, 0, "%-14.14s", "text", "Image date/time CCYYMMDDhhmmss", NULL, NULL, NULL, NULL},
   {"TGTID",      17, 0, "%-17.17s", "text", "Target ID", NULL, NULL, NULL, NULL},
   {"IID2",       80, 0, "%-80.80s", "text", "Image ID2", NULL, NULL, NULL, NULL},
   {"ISCLAS",      1, 0, "%-1.1s",   "text", "Image security classification", NULL, NULL, NULL, NULL},
   {"ISCLSY",      2, 0, "%-2.2s",   "text", "Classification system", NULL, NULL, NULL, NULL},
   {"ISCODE",     11, 0, "%-11.11s", "text", "Image codewords", NULL, NULL, NULL, NULL},
   {"ISCTLH",      2, 0, "%-2.2s",   "text", "Image control", NULL, NULL, NULL, NULL},
   {"ISREL",      20, 0, "%-20.20s", "text", "Image release instructions", NULL, NULL, NULL, NULL},
   {"ISDCTP",      2, 0, "%-2.2s",   "text", "Declass type", NULL, NULL, NULL, NULL},
   {"ISDCDT",      8, 0, "%-8.8s",   "text", "Declass date", NULL, NULL, NULL, NULL},
   {"ISDCXM",      4, 0, "%-4.4s",   "text", "Declass exemption", NULL, NULL, NULL, NULL},
   {"ISDG",        1, 0, "%-1.1s",   "text", "Image downgrade", NULL, NULL, NULL, NULL},
   {"ISDGDT",      8, 0, "%-8.8s",   "text", "Downgrade date", NULL, NULL, NULL, NULL},
   {"ISCLTX",     43, 0, "%-43.43s", "text", "class txt", NULL, NULL, NULL, NULL},
   {"ISCATP",      1, 0, "%-1.1s",   "text", "Image class authority type", NULL, NULL, NULL, NULL},
   {"ISCAUT",     40, 0, "%-40.40s", "text", "Authority", NULL, NULL, NULL, NULL},
   {"ISCRSN",      1, 0, "%-1.1s",   "text", "Classification reason", NULL, NULL, NULL, NULL},
   {"ISSRDT",      8, 0, "%-8.8s",   "text", "Security source date", NULL, NULL, NULL, NULL},
   {"ISCTLN",     15, 0, "%-15.15s", "text", "Security control number", NULL, NULL, NULL, NULL},
   {"ENCRYP",      1, 0, "%01d",     "int",  "Encryption", NULL, NULL, NULL, NULL},
   {"ISORCE",     42, 0, "%-42.42s", "text", "Image source", NULL, NULL, NULL, NULL},
   {"NROWS",       8, 0, "%08d",     "int",  "Number of rows", NULL, NULL, NULL, NULL},
   {"NCOLS",       8, 0, "%08d",     "int",  "Number of columns", NULL, NULL, NULL, NULL},
   {"PVTYPE",      3, 0, "%-3.3s",   "text", "Pixel value type = INT,B,SI,R,C", NULL, NULL, NULL, NULL},
   {"IREP",        8, 0, "%-8.8s",   "text", "Pixel value type", NULL, NULL, NULL, NULL},
   {"ICAT",        8, 0, "%-8.8s",   "text", "Image category", NULL, NULL, NULL, NULL},
   {"ABPP",        2, 0, "%02d",     "int",  "Actual bits per pixel", NULL, NULL, NULL, NULL},
   {"PJUST",       1, 0, "%-1.1s",   "text", "Pixel justification if ABPP<NBPP", NULL, NULL, NULL, NULL},
   {"ICORDS",      1, 0, "%-1.1s",   "text", "Image coordinate system", NULL, NULL, NULL, NULL},
   {"IGEOLO",     60, 0, "%-60.60s", "text", "Image geographic location", NULL, NULL, NULL, NULL},
   {"NICOM",       1, 0, "%01d",     "int",  "Number of 80-byte image comments", NULL, NULL, NULL, NULL},
   {"IC",          2, 0, "%-2.2s",   "text", "Image compression", NULL, NULL, NULL, NULL},
// If IC != "NM" or "NC"
//   {"COMRAT",      4, 0, "%-4.4s",   "text", "Compression rate code", NULL, NULL, NULL, NULL},
   {"NBANDS",      1, 0, "%01d",     "int",  "Number of bands", NULL, NULL, NULL, NULL},
// If NBANDS > 9
//   {"XBANDS",      5, 0, "%05d",     "int",  "Number of multi-spectral bands", NULL, NULL, NULL, NULL},
// Repeat NBANDS (or XBANDS) times:
   {"IREPBAND01",  2, 0, "%-2.2s",   "text", "Band 1 representation", NULL, NULL, NULL, NULL},
   {"ISUBCAT01",   6, 0, "%-6.6s",   "text", "Band 1 subcategory", NULL, NULL, NULL, NULL},
   {"IFC01",       1, 0, "%-1.1s",   "text", "Band 1 filter condition", NULL, NULL, NULL, NULL},
   {"IMFLT01",     3, 0, "%-3.3s",   "text", "Band 1 filter code", NULL, NULL, NULL, NULL},
   {"NLUTS01",     1, 0, "%01d",     "int",  "Band 1 number of LUTS", NULL, NULL, NULL, NULL},
// Lots of complication in here for multiband...
// end of NBAND/XBAND repeat
   {"ISYNC",       1, 0, "%01d",     "int",  "Image sync code", NULL, NULL, NULL, NULL},
   {"IMODE",       1, 0, "%-1.1s",   "text", "Image mode", NULL, NULL, NULL, NULL},
   {"NBPR",        4, 0, "%04d",     "int",  "Number of blocks per row", NULL, NULL, NULL, NULL},
   {"NBPC",        4, 0, "%04d",     "int",  "Number of blocks per column", NULL, NULL, NULL, NULL},
   {"NPPBH",       4, 0, "%04d",     "int",  "Pix/block horizontal", NULL, NULL, NULL, NULL},
   {"NPPBV",       4, 0, "%04d",     "int",  "Pix/block vertical", NULL, NULL, NULL, NULL},
   {"NBPP",        2, 0, "%02d",     "int",  "Number of bits per pixel", NULL, NULL, NULL, NULL},
   {"IDLVL",       3, 0, "%03d",     "int",  "Display level", NULL, NULL, NULL, NULL},
   {"IALVL",       3, 0, "%03d",     "int",  "Attachment level", NULL, NULL, NULL, NULL},
   {"ILOC",       10, 0, "%010d",    "int",  "Image location", NULL, NULL, NULL, NULL},
   {"IMAG",        4, 0, "%-4.4s",   "text", "Image magnification", NULL, NULL, NULL, NULL},
   {"UDIDL",       5,-3, "%05d",     "int",  "User-defined data length", NULL, NULL, NULL, NULL},
// If UDIDL > 0
//   {"UDOFL",       3, 0, "%03d",     "int",  "User-defined overflow", NULL, NULL, NULL, NULL},
//   {"UDID",        N, 0, "%s",       "int",  "User-defined Image Data", NULL, NULL, NULL, NULL},
   {"IXSHDL",      5,-3, "%05d",     "int",  "Extended subheader data length", NULL, NULL, NULL, NULL}
// If IXSHDL > 0
//   {"IXSOFL",      3, 0, "%03d",     "int",  "Extended subheader overflow", NULL, NULL, NULL, NULL}
//   {"IXSHD",       x, 0, "%03d",     "int",  "Image Extended Subheader Data", NULL, NULL, NULL, NULL}
};

static int nitf_image_len=(sizeof(nitf_image_head)/sizeof(NITF_KEYVAL));

/* NITF RPC TRE */
NITF_KEYVAL rpc00b_tre[]={
      {"CETAG",              6, 0, "%-6.6s",  "text", "TRE identifier", NULL, NULL, NULL, NULL},
      {"CEL",                5, 0, "%05d",    "int",  "Length of tagged record", NULL, NULL, NULL, NULL},
      {"SUCCESS",            1, 0, "%1d",     "int",  "success?", NULL, NULL, NULL, NULL},
      {"ERR_BIAS",           7, 0, "%7.2f",   "m",    "RMS systematic error", NULL, NULL, NULL, NULL},
      {"ERR_RAND",           7, 0, "%7.2f",   "m",    "RMS random error", NULL, NULL, NULL, NULL},
      {"LINE_OFF",           6, 0, "%06d",    "pix",  "Line offset", NULL, NULL, NULL, NULL},
      {"SAMP_OFF",           5, 0, "%05d",    "pix",  "Sample offset", NULL, NULL, NULL, NULL},
      {"LAT_OFF",            8, 0, "%8.4f",   "deg",  "Latitude offset", NULL, NULL, NULL, NULL},
      {"LONG_OFF",           9, 0, "%9.4f",   "deg",  "Longitude offset", NULL, NULL, NULL, NULL},
      {"HEIGHT_OFF",         5, 0, "%05d",    "m"  ,  "Height offset", NULL, NULL, NULL, NULL},
      {"LINE_SCALE",         6, 0, "%06d",    "pix",  "Line scale", NULL, NULL, NULL, NULL},
      {"SAMP_SCALE",         5, 0, "%05d",    "pix",  "Sample scale", NULL, NULL, NULL, NULL},
      {"LAT_SCALE",          8, 0, "%8.4f",   "deg",  "Latitude scale", NULL, NULL, NULL, NULL},
      {"LONG_SCALE",         9, 0, "%9.4f",   "deg",  "Longitude scale", NULL, NULL, NULL, NULL},
      {"HEIGHT_SCALE",       5, 0, "%5d",     "m",    "Height scale", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_1",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_2",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_3",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_4",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_5",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_6",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_7",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_8",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_9",  12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_10", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_11", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_12", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_13", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_14", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_15", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_16", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_17", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_18", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_19", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_NUM_COEFF_20", 12, 0, "%12.6e",   "RPC", "Line numerator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_1",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_2",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_3",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_4",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_5",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_6",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_7",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_8",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_9",  12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_10", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_11", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_12", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_13", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_14", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_15", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_16", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_17", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_18", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_19", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"LINE_DEN_COEFF_20", 12, 0, "%12.6e",   "RPC", "Line denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_1",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_2",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_3",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_4",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_5",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_6",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_7",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_8",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_9",  12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_10", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_11", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_12", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_13", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_14", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_15", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_16", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_17", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_18", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_19", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_NUM_COEFF_20", 12, 0, "%12.6e",   "RPC", "Sample numerator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_1",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_2",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_3",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_4",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_5",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_6",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_7",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_8",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_9",  12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_10", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_11", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_12", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_13", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_14", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_15", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_16", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_17", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_18", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_19", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL},
      {"SAMP_DEN_COEFF_20", 12, 0, "%12.6e",   "RPC", "Sample denominator coeff", NULL, NULL, NULL, NULL}
};
static int rpc00b_len=(sizeof(rpc00b_tre)/sizeof(NITF_KEYVAL));

/* NITF STDIDC TRE */
NITF_KEYVAL stdidc_tre[]={
      {"CETAG",              6, 0, "%-6.6s",   "text", "TRE identifier", NULL, NULL, NULL, NULL},
      {"CEL",                5, 0, "%05d",     "int",  "Length of tagged record", NULL, NULL, NULL, NULL},
      {"ACQUISITION_DATE",  14, 0, "%-14.14s", "text", "Acquisition Date", NULL, NULL, NULL, NULL},
      {"MISSION",           14, 0, "%-14.14s", "text", "Mission identification", NULL, NULL, NULL, NULL},
      {"PASS",               2, 0, "%-4.4s",   "text", "Pass number", NULL, NULL, NULL, NULL},
      {"OP_NUM",             3, 0, "%03d",     "int",  "Image operation number", NULL, NULL, NULL, NULL},
      {"START_SEGMENT",      2, 0, "%-2.2s",   "text", "Start Segment ID", NULL, NULL, NULL, NULL},
      {"REPRO_NUM",          2, 0, "%02d",     "int",  "Reprocess Number", NULL, NULL, NULL, NULL},
      {"REPLAY_REGEN",       3, 0, "%-3.3s",   "text", "Replay", NULL, NULL, NULL, NULL},
      {"BLANK_FILL",         1, 0, "%1.1s",    "text", "Blank Fill", NULL, NULL, NULL, NULL},
      {"START_COLUMN",       3, 0, "%03d",     "int",  "Starting Column Block", NULL, NULL, NULL, NULL},
      {"START_ROW",          5, 0, "%05d",     "int",  "Starting Row Block", NULL, NULL, NULL, NULL},
      {"END_SEGMENT",        2, 0, "%-2.2s",   "text", "Ending Segment ID of this file", NULL, NULL, NULL, NULL},
      {"END_COLUMN",         3, 0, "%03d",     "int",  "Ending Column Block", NULL, NULL, NULL, NULL},
      {"END_ROW",            5, 0, "%05d",     "int",  "Ending Row Block", NULL, NULL, NULL, NULL},
      {"COUNTRY",            2, 0, "%-2.2s",   "text", "Country Code", NULL, NULL, NULL, NULL},
      {"WAC",                4, 0, "%04d",     "int",  "World Aeronautical Chart", NULL, NULL, NULL, NULL},
      {"LOCATION",          11, 0, "%-11.11s", "text", "Location", NULL, NULL, NULL, NULL},
      {"reserved",           5, 0, "%-5.5s",   "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           8, 0, "%-8.8s",   "text", "reserved", NULL, NULL, NULL, NULL},
};
static int stdidc_len=(sizeof(stdidc_tre)/sizeof(NITF_KEYVAL));

/* NITF USE00A TRE */
NITF_KEYVAL use00a_tre[]={
      {"CETAG",              6, 0, "%-6.6s", "text", "TRE identifier", NULL, NULL, NULL, NULL},
      {"CEL",                5, 0, "%05d",   "int",  "Length of tagged record", NULL, NULL, NULL, NULL},
      {"ANGLE_TO_NORTH",     3, 0, "%03d",   "deg",  "Angle to North.", NULL, NULL, NULL, NULL},
      {"MEAN_GSD",           5, 0, "%05d",   "inch", "Mean Ground Sample Distance.", NULL, NULL, NULL, NULL},
      {"reserved",           1, 0, "%-1.1s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"DYNAMIC_RANGE",      5, 0, "%05d",   "ADU",  "Dynamic Range", NULL, NULL, NULL, NULL},
      {"reserved",           3, 0, "%-3.3s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           1, 0, "%-1.1s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           3, 0, "%-3.3s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"OBL_ANG",            5, 0, "%5.2f",  "deg",  "Obliquity Angle", NULL, NULL, NULL, NULL},
      {"ROLL_ANG",           6, 0, "%6.2f",  "deg",  "Roll Angle", NULL, NULL, NULL, NULL},
      {"reserved",          12, 0, "%-12.12s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",          15, 0, "%-15.15s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           4, 0, "%-4.4s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           1, 0, "%-1.1s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           3, 0, "%-3.3s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           1, 0, "%-1.1s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           1, 0, "%-1.1s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"N_REF",              2, 0, "%02d",   "int",  "Number of Reference Lines", NULL, NULL, NULL, NULL},
      {"REV_NUM",            5, 0, "%05d",   "int",  "Revolution Number", NULL, NULL, NULL, NULL},
      {"N_SEG",              3, 0, "%03d",   "int",  "Number of Segments", NULL, NULL, NULL, NULL},
      {"MAX_LP_SEG",         6, 0, "%06d",   "int",  "Maximum Lines Per Segment", NULL, NULL, NULL, NULL},
      {"reserved",           6, 0, "%-6.6s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"reserved",           6, 0, "%-6.6s", "text", "reserved", NULL, NULL, NULL, NULL},
      {"SUN_EL",             5, 0, "%5.1f",  "deg",  "Sun Elevation.", NULL, NULL, NULL, NULL},
      {"SUN_AZ",             5, 0, "%5.1f",  "deg",  "Sun Azimuth.", NULL, NULL, NULL, NULL}
};
static int use00a_len=(sizeof(use00a_tre)/sizeof(NITF_KEYVAL));

/* Prototypes */
/* Return the index of a field in a keyval array */
int idxfield(int n, NITF_KEYVAL *kv, char *field);
/* Convert a string of bytes to short ints */
void charshort(int n, unsigned char *in, unsigned short *out);
/* Convert a string of bytes to floats */
void charfloat(int n, unsigned char *in, float *out);
/* Convert a string of ints to short ints */
void longshort(int n, int *in, unsigned short *out);
/* Convert a string of floats to short ints */
void floatshort(int n, float *in, unsigned short *out);
/* Reformat one block of data into another */
int reform(int mx, int my, int mbpp, void *mdata, int ix, int iy, 
	   int nx, int ny, int bitpix, void *data, int yflip);
/* Read a field's value from fp */
int readfield(NITF_KEYVAL *kv, FILE *fp, int verb);
/* Parse the nitf_image structure for data size parameters */
void image_size(int *nx, int *ny, int *nbx, int *nby, 
		int *nxblk, int *nyblk, int *nbpp);

/* Convert an RPC00A into RPC00B */
int permute_rpcA2B(double *rpc);

/* Build a FITS header from NITF quantities */
void mkfitshead(char **fitshead, int bitpix, int nx, int ny);

double jdate(int year, int month, int day, double ut);
