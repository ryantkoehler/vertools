/*
* venpipe.h
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


#include "dna.h"

#define VEN_VERSION_S   "VenPipe Version 1.11"
#define VLIB_VERSION_S  "Vienna Library: ViennaRNA-2.3.1"

#define MAX_VSLEN   10000   /* Sequence buffer size; Max seq len */


#define VENPIPE_ID      5161
typedef struct VENPIPE
{
    int ID;
    char inname[NSIZE]; /* Input filename */
    FILE *in;               /* Input file */
    int iform;              /* Input format */
    int cleanseq;           /* Input sequence cleaning level */
    char outname[NSIZE];    /* Output log name */
    FILE *out;              /* Log file */
    int ofas;               /* Flag to output sequence */
    int oraw;               /* Flag to output sequence */
    int do_pfe;             /* Flag to calculate partition free energy; not min */
    int do_ddb;             /* Flag to dump dot-bracket seq structrure */
    int do_dmb;             /* Flag to dump number of matching bases */
    int do_ksapar;          /* Flag to keep salt adjusted parameter file */
    int firstb,lastb;       /* Base restrictions? */
    int rre;                /* Range relative to end flag */
    int do_ds;              /* flag to dump seq */
    int do_mask;            /* flag to mask range with N */
    int do_not;             /* flag to invert mask range */
    int dmb_f,dmb_l,mrre;   /* Range restrictions for dmb reporting */
    int do_mbtab;           /* flag to dump match base table */
    DOUB mst,men,mj;        /* Melt temperature start, end, and jump */
    /***
    *   Tm_pars / Vienna stuff
    */
    DOUB temp;              /* Temperature */
    DOUB salt;              /* Salt */
    DOUB tcon,pcon;         /* Target and probe concentrations */
    struct VENERGY *ven;    /* Vienna-specific settings structure */
    char vparfile[NSIZE];   /* Vienna parameter file name */
    int do_saltcorrect;     /* Flag to salt correct vienna parameters */
    /***
    *   per seq stuff
    */
    struct SEQ *seq;        /* SEQUENCE Object */
    char tname[NSIZE];      /* target name string */
    char tseq[MAX_VSLEN];   /* sequence string buffer (after cleaned up) */
    char tseq2[MAX_VSLEN];  /* sequence string buffer; no masking */
    int tlen;               /* sequence length */
    char tss[MAX_VSLEN];    /* sec structure string */
    char tss2[MAX_VSLEN];   /* sec structure string 2 */
    DOUB ten, ten2;         /* structure energies */
}VENPIPE;

/*
*   Structure to hold various vienna settings and parameters
*/
#define VENERGY_ID      4082    
typedef struct VENERGY
{
    int ID;
    char vparfile[NSIZE];   /* Vienna parameter file */
    int do_saltcorrect;     /* Flag to salt correct vienna parameters */
    int do_saltkeep;        /* Flag to keep salt-corrected par file */
    DOUB temp;              /* Temperature */
    DOUB salt;              /* Salt, Mg */

/*
    int do_backtrack;       
    int nogu;
    int noclosinggu;
    int nolonelypairs;
    int tetra_loop;
    int energy_set;
    int fold_constrained;
    int dangles;
    int logml;
*/
}VENERGY;


#define CHECK_VENPIPE(t)    if(t){DestroyVenpipeI(t); t=NULL;}
#define CHECK_VENERGY(t)    if(t){DestroyVenergyI(t); t=NULL;}


/*********************************************************************
*   Global flag vars 
*/
#ifdef __VENERGY__
int vienna_okGI = 0;
#else
extern int vienna_okGI;
#endif


/***
*   Defaults 
*/
#define DEF_TEMP            37      /* Default temperature (C) */
#define DEF_SALT            1.00    /* Default salt conc (M) */

#define DEF_SALTCORRECT     TRUE    /* Flag to salt-correct parameter file */
#define SALTCOR_MAXFILE     20      /* Max number for salt-correct temp file */

#define DEF_MELT_JUMP_D     2.5     /* elta for "melting" */
#define MIN_MELT_JUMP_D     1.0     /* in delta for "melting" */



/*********************** ppp ********************
* C function listing generated by gen_prot
* Thu Jan 21 17:12:46 2016
*/
/****************************************************************
* ven_io.c
*/
void WriteVenpipeHeader(VENPIPE *vPO, char *nameS, FILE *oPF);
void DumpVenergy(VENERGY *vePO, FILE *oPF);
void HandleVenpipeOut(VENPIPE *venPO, FILE *outPF);
int DumpMatchBaseTablesI(VENPIPE *venPO, char *nameS, int len, int *sscPI, 
    int *psscPI, FILE *outPF);
int DumpOneMatchBaseTabI(VENPIPE *venPO, char *nameS, int len, int *sscPI, 
    FILE *outPF);
int GetMatchBaseSeqCoordsI(VENPIPE *vpPO, int len, int *stPI, int *enPI, int *rPI);
int FillSeqStructMatchOutI(VENPIPE *vpPO, int *sscPI, int *psscPI, int len, char *matS);
int StrucStringToArrayI(char *seqS, int len, int *aPI);

/****************************************************************
* ven_engy.c
*/
void DumpViennaGlobals(FILE *fPF);
int SetUpViennaEnvI(VENERGY *vePO);
int SetViennaTemperatureI(VENERGY *vePO, DOUB tempD);
int SetViennaSaltI(VENERGY *vePO, DOUB saltD, int use, int keep);
int SetViennaParametersI(VENERGY *vePO, char *parS);
int CleanViennaEnvI(void);
int GetViennaStructI(char *rseqS,int slen, char *ssS, DOUB *enPD,
                        int do_pfold, char *pssS, DOUB *penPD);
int SetSeqViennaVarsI(char *seqS, DOUB enD) ;
int GetViennaStructEnergyI(char *rseqS,char *ssS,DOUB *enPD);

/****************************************************************
* ven_str.c
*/
VENPIPE *CreateVenpipePO(void);
int DestroyVenpipeI(VENPIPE *vpPO);
int InitVenpipeI(VENPIPE *vpPO, int full);
VENERGY *CreateVenergyPO(void);
int DestroyVenergyI(VENERGY *vePO);
int InitVenergyI(VENERGY *vePO);
int RealizeVenpipeI(VENPIPE *vpPO);
int CopyWorkingSeqI(VENPIPE *venPO, SEQ *seqPO);
void ChangeUtoT(char *seqS,int len);
int PrepViennaInputSeqI(char *seqS, int len,char *newS);
int GetVenpipeStructureI(VENPIPE *vpPO);

/****************************************************************
* ven_salt.c
*/
int SaltCorrectViennaParsI(char *inS, DOUB saltD, char *newS);
int OpenSaltCorOutputFileI(char *inS, char *newS, FILE **outPPF);
int ParseViennaBlockTokenI(char *wordS);
DOUB SaltCorrectVenNumberD(DOUB oD,int type,DOUB saltD);

/****************************************************************
* venpipe.c
*/
int main(int argc, char **argv);
void VenPipeUse(void);
int VenPipeI(int argc, char **argv);
int OkVenpipeOptionsI(VENPIPE *vpPO);
int HandleVenpipeMeltingI(VENPIPE *venPO,FILE *outPF);

