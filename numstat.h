/*
* numstat.h
*
* Copyright 2018 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
*
* The programs and source code of the vertools collection are free software.
* They are distributed in the hope that they will be useful,
* WITHOUT ANY WARRANTY OF FITNESS FOR ANY PARTICULAR PURPOSE.  
* 
* Permission is granted for research, educational, and commercial use and 
* modification as long as 1) Code and any derived works are not redistributed
* for any fee, and 2) Proper credit is given to the authors. If you wish to 
* include this software in a product, please contact the authors.
*
* See https://www.verdascend.com/ for more
*
*/


#define VERSION_S "NumStat Version 2.31"


#define MAX_WIDTH       10000   /* Size of lines to eat */
#define HPLOT_SLEN      40      /* Width of histogram plotting space */
#define MIN_HIS_BIN     0.00001 /* Minimum histogram bin size */
#define DEF_HMB         20      /* Default max number of histogram bins */
#define TRUNC_FOLD      2.0     /* truncate single bins if it's X-fold > than next */

#define DEF_COL         2       /* default column */

/* Percentile reporting format string */
#define PERCENTILE_FORM_S   "perc_%02d \t%8.4f\t%6.2f\n" 
#define PERC_SL_FORM_S      "\t%8.4f" 
#define DEF_H_PTMF_S        "\t%8.4f" 
#define DEF_H_SEP_S         "\t" 


#define NUMSTAT_ID      5101
typedef struct NUMSTAT
{
    int ID;
    char inname[NSIZE];     /* Input file name */
    FILE *in;               /* Input file */
    char outname[NSIZE];    /* Output file name */
    FILE *out;              /* Outpuf file */
    int col,scol;           /* Column and second column */
    int col_set;            /* Flag that column is set */
    struct NUMLIST *vals;   /* Data values */
    struct NUMLIST *svals;  /* Second set of data values */
    int num;                /* Number of values */
    int lines;              /* Lines in input file */
    char prls[BBUFF_SIZE];  /* Percentile series string */
    DOUB prsd;              /* Percentile steps */
    DOUB min,max,av,sd,sum; /* Run-time stat values */
    int h_maxbin;           /* List max bins */
    int h_mbin_set;         /* Flag that hist max bins is set */
    int h_pwide;            /* Historam plot width */
    DOUB h_bin;             /* Histogram specified bin */         
    DOUB h_lo, h_hi;        /* Histogram specified low and hi */
    DOUB h_abin;            /* Histogram auto bin */         
    int h_limbin;           /* Histogram auto bin */         
    DOUB h_alo, h_ahi;      /* Histogram auto low and hi */
    char h_pfmt[DEF_BS];    /* Histogram print format string */
    char hvsep[DEF_BS];     /* Histogram value sep string */
    int do_perc, do_hist;   /* Run-time output modifiers */
    int do_hntb;            /* Flag for Histogram Not Trim Bin */
    DOUB htb_xfold;         /* Fist trim bin if highest > xfold */
    int do_sk;              /* Flag to skip missing data lines */
    int do_sl;              /* Flag for single line */
    int do_efi;             /* Flag for echo */
    int do_ic;              /* Flag to inore chars up to this position */
    int do_hplot;           /* Flag to plot hist */
    int do_splot;           /* Flag to plot scatter */
    int do_hends;           /* Flag for no "ends" in hist plotting */
    int do_ploti;           /* Flag to plot hist integral */
    int do_dif;             /* Report diff between C1 and C2 */
    char echo[NSIZE];       /* String for output echo */
}NUMSTAT;

#define CHECK_NUMSTAT(fp)   if(fp){DestroyNumstatI(fp);fp=NULL;}


/*********************** ppp ********************
* C function listing generated by gen_prot
* Wed Mar 16 11:10:29 2016
*/
/****************************************************************
* numstat.c
*/
int main(int argc, char **argv);
void NumstatUse(void);
int NumstatI(int argc, char **argv);
int NumstatHandleHistI(NUMSTAT *nsPO);
void OneHistLineNumsOut(NUMSTAT *nsPO, int n, DOUB stD, char *bmS, int ncum );
void OneHistLinePlotOut(NUMSTAT *nsPO, int n, int ncum, int nmax);
void NumstatHistHeader(NUMSTAT *nsPO, HISTOGRAM *hisPO, FILE *outPF);
void NumstatHistStatHeader(NUMSTAT *nsPO, HISTOGRAM *hisPO, char *shamS, int two, 
    FILE *outPF);
void NumstatHistLabHeader(NUMSTAT *nsPO, FILE *outPF);
int NumstatHandleSplotI(NUMSTAT *nsPO);
NUMSTAT *CreateNumstatPO(void);
int DestroyNumstatI(NUMSTAT *nsPO);
int AddNumstatSecValsI(NUMSTAT *nsPO) ;
void InitNumstat(NUMSTAT *nsPO);
int LoadColValsFromFileI(NUMSTAT *nsPO);
int OpenNumstatFilesI(NUMSTAT *nsPO) ;
int CheckNumstatOptionsI(NUMSTAT *nsPO);
int HandleHisMaxBinsI(NUMSTAT *nsPO);
int GetLineDataValI(NUMSTAT *nsPO, char *bufS, DOUB *rPD, DOUB *sPD);
int GetWordDataValI(char *linePC, int col, int skip, DOUB *vPD);
void NumstatOneLineOut(NUMSTAT *nsPO);
void NumstatHandleHeader(NUMSTAT *nsPO);
void NumstatReportStats(NUMSTAT *nsPO, int scol, char *exS, FILE *outPF);
void HandleListPercentiles(NUMSTAT *nsPO, DOUB *rvalsPD, int num, FILE *outPF);
void HandleStepPercentiles(NUMSTAT *nsPO, DOUB *rvalsPD, int num, FILE *outPF);
void ReportSingPercentile(NUMSTAT *nsPO, int p, DOUB vD, int num, FILE *outPF);
int CalcColDifsI(NUMSTAT *nsPO);

