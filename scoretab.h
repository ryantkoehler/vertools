/*
* scoretab.h
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
*
* The programs and source code of the vertools collection are free software.
* They are distributed in the hope that they will be useful,
* WITHOUT ANY WARRANTY OF FITNESS FOR ANY PARTICULAR PURPOSE.  
* 
* Permission is granted for research, educational, and possibly commercial use 
*   and modification as long as 1) Code and any derived works are not 
*   redistributed for any fee, and 2) Proper credit is given to the authors. 
*   If you wish to include this software in a product, or use it commercially,
*   please contact the authors.
*
* See https://www.verdascend.com/ for more
*
*/

#define VERSION_S   "ScoreTab Version 1.13"

#define SCORETAB_ID     5131
typedef struct SCORETAB
{
    int ID;
    char inname[NSIZE];     /* nput filename */
    char sdfname[NSIZE];    /* core definition filename */
    FILE *out;              /* utput file */
    char outname[NSIZE];    /* utput filename */
    int owhat;              /* utput what? */
    char mergname[NSIZE];   /* erge filename */
    char macpre[NSIZE];     /* erge append col prefix */
    char macsuf[NSIZE];     /* erge append col suffix */
    int rlab,clab;          /* lag to expect table with row/colum lables */
    int corn;               /* lag that table has no row-col label "corner" */
    int do_skp;             /* lag to skip problem input lines */
    int do_mlis;            /* lag to treat merge as list of tables */
    int do_msub;            /* lag to subtract merged tables */
    int do_mmul;            /* lag to multiply merged tables */
    int do_mdiv;            /* lag to div merged tables */
    int do_mmin;            /* lag to take min of merged tables */
    int do_mmax;            /* lag to take max of merged tables */
    int do_mapc;            /* erge by appanding columns */
    int do_mrow;            /* eport row merges */
    int do_abs;             /* lag to replace with absolute value */
    int do_symu,do_symd;    /* lag to symmetrize up (max) or down (min) */
    int do_tran;            /* transpose tables */
    int do_stru,do_strd;    /* flags to sort rows up / down */
    int ccor;               /* correlate this column to others */
    int cinfo;              /* calculate info in this column relative to others */
    int inf_maxhb;          /* information calc max histogram bins */
    int minrow,maxrow;      /* bounds on rows to consider */
    int mincol,maxcol;      /* bounds on cols to consider */
    char rlisname[NSIZE];   /* file with row subset list */
    char clisname[NSIZE];   /* file with col subset list */
    int do_kc;              /* flag to keep case in name comparisions */
    int do_wst;             /* flag for word start only for name compare */
    int do_wsub;            /* flag for word substr only for name compare */
    int do_igd;             /* flag to ignore diagonal */
    int do_wfm;             /* flag to do weak-first matching */
    int do_dwfm;            /* flag to dump weak-first matching matrix */
    int do_ncv;             /* flag to normalize column values */
    int do_nrd;             /* flag to normalize rows to diagonal values */
    int do_qcv;             /* flag / number to quantize column values */
    int do_srow;            /* flag to smooth rows */
    int do_scol;            /* flag to smooth rows */
    int partalg;            /* partition algorithm */
    int psize;              /* flag to partition into pools of this size */
    DOUB pmin;              /* partition minimum compatibility value allowed */
    char dumpbase[NSIZE];   /* partition matrix dumping file basename */
    char pvsep[NSIZE];      /* printed column separator string */
    int do_ocsv, do_ossv;   /* lag for output comma-sep-value or space-sep-value */
    int usp_form;           /* flag for user specific formatting */
    DOUB mval,sval;         /* value multipy and shift coefficients */
    DOUB blval,bhval;       /* bound low and high values */
    DOUB exp;               /* exponent */
    struct TABLE *tab;      /* table with numbers */
    struct TABLE *stab;     /* second table structure */
    struct NUMLIST *tvals1; /* array to hold temp row or col values */ 
    struct NUMLIST *tvals2; /* array to hold temp row or col values */ 
    int nscores;            /* number of score fields */
    struct SCFIELD **scores;/* array to pointers to score fields */
    int do_scg;             /* flag to apply score transforms globally */
    int quiet;
    DOUB flglo,flghi;       /* flagging low and high values */
    int do_not;             /* invert flagging criteria */
    int seed;               /* random number seed */
    int do_gax;             /* flag for some GA crossover operation */
    DOUB gaxr;              /* GA row crossover fraction */
    DOUB gaxc;              /* GA col crossover fraction */
    int do_gam;             /* flag for some GA mutation operation */
    DOUB gamf;              /* GA mutation fraction */
    DOUB gamg;              /* GA mutation gaussian */
}SCORETAB;

#define CHECK_SCORETAB(tu)    if(tu){DestroyScoretabI(tu); tu=NULL;}

/***
*   Default output formatting
*/
#define DEF_PRLFORM_S   "%s"
#define DEF_PVFORM_S    "%3.2f"
#define DEF_PVSEP_S     "\t"

/***
*    Output codes
*/
#define SCTO_FULL    61    /* ull dump */
#define SCTO_RPRO    65    /* ow product */
#define SCTO_RSTAT   68    /* ow stats */
#define SCTO_CSTAT   79    /* ol stats */
#define SCTO_FSTAT   87    /* ull stats */
#define SCTO_CCOR    90    /* ol correlation */
#define SCTO_FCCM    91    /* ul Col correlation matrix */
#define SCTO_CINFO   92    /* ol information measures */
#define SCTO_RCL     93    /* ow-Col list */
#define SCTO_CRL     94    /* ol-Row list */
#define SCTO_FLAG    95    /* lagged subset of values */
#define SCTO_NONE    99    /* othing at all */


#define INFO_DEF_MBIN   5   /* efault (max) bins for mutual information calc */


#define PALG_BTS    222    /* artition algorithm block-then-swap */

#define DEF_PALG    PALG_BTS    /* efault partition algoritm */



/*********************** ppp ********************
* C function listing generated by gen_prot
* Sat Jul 12 13:13:05 2014
*/
/****************************************************************
* sctab_ga.c
*/
int HandleSctGACrossOverI(SCORETAB *stPO);
int HandleSctGAMutationI(SCORETAB *stPO);

/****************************************************************
* sctab_mix.c
*/
int HandleSctTableMergingI(SCORETAB *stPO);
int HandleTableMergeI(SCORETAB *stPO, TABLE *tabPO, char *tabS);
int HandleTableClabPreSufI(TABLE *tabPO, char *preS, char *sufS);
int FigureMergeMixI(SCORETAB *stPO, char *whatS);
int HandleTransposeI(TABLE *tabPO, TABLE **tranPPO);
int HandleRowSortI(TABLE *tabPO, int sdir, int mask);

/****************************************************************
* sctab_part.c
*/
int HandleSetPartitioningI(SCORETAB *stPO, TABLE *tabPO);
int SetStartingPoolsI(SCORETAB *stPO,TABLE *tabPO, TABLE *poolsPO, int npools,
    int psize);
int GetPoolMemTabIndexI(TABLE *poolsPO, int p1, int m1, int *indPI);
int GetPoolPairValueI(TABLE *poolsPO, int p1, int m1, int p2, int m2,
    TABLE *tabPO, DOUB *vPD);
int ScreenPoolI(TABLE *poolsPO,int pool,TABLE *tabPO,DOUB minD, FILE *outPF);
int CanSwapOutI(TABLE *poolsPO,int pool,int i,int j,TABLE *tabPO, DOUB minD);
int HandlePoolSwapI(TABLE *poolsPO,int pool1,int m1,int pool2,int m2);
int OkInPoolI(TABLE *poolsPO, int old, int mem, int new, TABLE *tabPO, DOUB minD);
int KillThisGuyI(int pi,TABLE *poolsPO,int pool,TABLE *tabPO, FILE *outPF);
int DumpPoolMembersI(TABLE *poolsPO, int pool, TABLE *tabPO, FILE *outPF);
void ReportConflict(TABLE *tabPO,int pool, int m1,int m2,DOUB vD,FILE *outPF);
int MaskPoolMemValsI(TABLE *poolsPO,TABLE *tabPO,TABLE *dupPO, int pool);
int ScanShamsI(TABLE *poolsPO,int pool, TABLE *tabPO, DOUB minD);
void WritePartitionHeader(SCORETAB *stPO, TABLE *tabPO, int npools, FILE *outPF);
int FillPartitionAlgoStringI(int alg, char *nameS);

/****************************************************************
* sctab_sc.c
*/
int HandleScoreTransformI(SCORETAB *stPO, TABLE *tabPO);

/****************************************************************
* sctab_str.c
*/
SCORETAB *CreateScoretabPO(void);
int DestroyScoretabI(SCORETAB *stPO);
void InitScoretab(SCORETAB *stPO);
void SetOutFormatting(SCORETAB *stPO,int oftw,int oftp);
int CheckSctOutOptionsI(SCORETAB *stPO);
int SetUpSctAuxDataI(SCORETAB *stPO);
int SetUpSctTableSpaceI(SCORETAB *stPO);
int RestrictTableRangesI(SCORETAB *stPO, TABLE *tabPO);
int HandleMaskTableRowsI(SCORETAB *stPO, TABLE *tabPO);
int HandleMaskTableColsI(SCORETAB *stPO, TABLE *tabPO);
int HandleTabRowSubsetI(SCORETAB *stPO, TABLE *tabPO);
int HandleTabColSubsetI(SCORETAB *stPO, TABLE *tabPO);
int AdjustTableValuesI(SCORETAB *stPO,TABLE *tabPO);
void HandleSymmetrization(SCORETAB *stPO, TABLE *tabPO);

/****************************************************************
* scoretab.c
*/
int main(int argc, char **argv);
void ScoreTabUse(void);
int ScoreTabI(int argc, char **argv);
int HandleScoreRowMergeI(SCORETAB *stPO,TABLE *tabPO);
int HandleNormalizationI(SCORETAB *stPO,TABLE *tabPO);
int HandleSmoothingI(SCORETAB *stPO, TABLE *tabPO);
void HandleSctOutput(SCORETAB *stPO,FILE *outPF);
void FillColOutputHeadLine(TABLE *tabPO, int col, char *bufS);
void HandleSctOutList(SCORETAB *stPO,TABLE *tabPO,FILE *outPF);
void PrintTabLabsValLineI(TABLE *tabPO, char *fS, char *sS, DOUB vD, FILE *outPF);
void HandleSctOutRows(SCORETAB *stPO,TABLE *tabPO,FILE *outPF);
void HandleSctOutCols(SCORETAB *stPO,TABLE *tabPO,FILE *outPF);
void HandleSctOutStats(SCORETAB *stPO,TABLE *tabPO,int what,int which,
    FILE *outPF);
void HandleSctSingColInfo(SCORETAB *stPO,TABLE *tabPO,int c1, int c2,
    FILE *outPF);
void HandleSctSingColCors(SCORETAB *stPO,TABLE *tabPO,int c1, int c2,
    FILE *outPF);
void HandleSctColCors(SCORETAB *stPO,TABLE *tabPO,FILE *outPF);
void HandleSctFlaggedOut(SCORETAB *stPO,TABLE *tabPO,FILE *outPF);
void HandleSctWfmOut(SCORETAB *stPO,TABLE *vtabPO,TABLE *tabPO,FILE *outPF);

