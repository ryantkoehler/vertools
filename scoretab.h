/*
* scoretab.h
*
* Copyright 2017 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
*/

#define VERSION_S   "ScoreTab Version 1.1"

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
    int do_tran;            /* ranspose tables */
    int do_stru,do_strd;    /* lags to sort rows up / down */
    int ccor;               /* orrelate this column to others */
    int cinfo;              /* alculate info in this column relative to others */
    int inf_maxhb;          /* nformation calc max histogram bins */
    int minrow,maxrow;      /* ounds on rows to consider */
    int mincol,maxcol;      /* ounds on cols to consider */
    char rlisname[NSIZE];   /* ile with row subset list */
    char clisname[NSIZE];   /* ile with col subset list */
    int do_kc;              /* lag to keep case in name comparisions */
    int do_wst;             /* lag for word start only for name compare */
    int do_wsub;            /* lag for word substr only for name compare */
    int do_igd;             /* lag to ignore diagonal */
    int do_wfm;             /* lag to do weak-first matching */
    int do_dwfm;            /* lag to dump weak-first matching matrix */
    int do_ncv;             /* lag to normalize column values */
    int do_qcv;             /* lag / number to quantize column values */
    int do_srow;            /* lag to smooth rows */
    int do_scol;            /* lag to smooth rows */
    int partalg;            /* artition algorithm */
    int psize;              /* lag to partition into pools of this size */
    DOUB pmin;              /* artition minimum compatibility value allowed */
    char dumpbase[NSIZE];   /* artition matrix dumping file basename */
    char pvsep[NSIZE];      /* rinted column separator string */
    int do_ocsv, do_ossv;   /* lag for output comma-sep-value or space-sep-value */
    int usp_form;           /* lag for user specific formatting */
    DOUB mval,sval;         /* alue multipy and shift coefficients */
    DOUB blval,bhval;       /* ound low and high values */
    DOUB exp;               /* xponent */
    struct TABLE *tab;      /* able with numbers */
    struct TABLE *stab;     /* econd table structure */
    struct NUMLIST *tvals1; /* rray to hold temp row or col values */ 
    struct NUMLIST *tvals2; /* rray to hold temp row or col values */ 
    int nscores;            /* umber of score fields */
    struct SCFIELD **scores;/* rray to pointers to score fields */
    int do_scg;             /* lag to apply score transforms globally */
    int quiet;
    DOUB flglo,flghi;       /* lagging low and high values */
    int do_not;             /* nvert flagging criteria */
    int seed;               /* andom number seed */
    int do_gax;             /* lag for some GA crossover operation */
    DOUB gaxr;              /* A row crossover fraction */
    DOUB gaxc;              /* A col crossover fraction */
    int do_gam;             /* lag for some GA mutation operation */
    DOUB gamf;              /* A mutation fraction */
    DOUB gamg;              /* A mutation gaussian */
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

