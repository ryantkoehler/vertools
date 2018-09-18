/*
* blastout.h
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
*/


#define VERSION_S   "BlastOut version 0.7"

#define BLASTOUT_ID 5031
#define BLASTANS_ID 6005

#define BLBSIZE     1000     /* blastout buff size */
#define HISTDIM     300      /* max identities to tabulate in hist dim */


typedef struct BLASTOUT
{
    int ID;
    char outname[NSIZE];    /* Output file name */
    FILE *out;              /* Output file */
    int owhat;              /* Flag for output processing */
    int dhit;               /* Flag to dump out all hits (summary) */
    int dseq;               /* Flag to dump seqs for all hits */
    char nsqh[NSIZE];       /* ame output seqs query + XXX + hit */
    int do_smu, do_sml;     /* Flags to Sequence Mark Case on sequence dump */
    int do_align;           /* Flag to dump alignment with sequences */
    int dhc;                /* Flag to dump query hit coordinates */
    int dsco;               /* Flag to dump hit score */
    int dmbc;               /* Flag to dump match base counts */
    int dsum;               /* Flag to dump summary / query */
    int do_dfml;            /* Flag to dump fraction map list */
    int do_con;             /* Flag to dump contiguous match */
    int do_co3;             /* Flag contiguous match from 3' end */
    int do_dci;             /* Flag to dump condensed record */
    int dhis;               /* Flag/num to dump histogram of hits */
    int phis;               /* Padding flag for histogram output */
    int chis;               /* Flag for cumulative histogram */
    int firstb,lastb;       /* First and last query bases to consider */
    int rre;                /* Flag for base range relative to end */
    int do_rkc;             /* Flag for range keep case; def = UP/low */
    int mid;                /* Minimum ident to qualify */
    REAL mif;               /* Minimum ident fraction qualify */
    int do_mmm;             /* Filter with more than this many mismatches */
    int do_mgap;            /* Filter with more than this many gaps */
    int do_frq;             /* Flag to filter relative to query (not align / bran) */
    int do_fnot;            /* Flag for NOT logic on minimum match qualifications */
    int firstq,lastq;       /* First and last query to report */
    int firsth,lasth;       /* First and last hit to report */
    char opq[NSIZE];        /* Output per query file extension string */
    struct BLASTANS *ans;   /* Anwsers structure */
}BLASTOUT;

/***
*   Run time stuff; blast answers
*/
typedef struct BLASTANS
{
    int ID;
    char input[NSIZE];  /* Input file (i.e. log) name */
    int itype;          /* Flag for input type (i.e. blast output) */
    FILE *in;           /* Input file */
    char query[NSIZE];  /* Name of (current) query */
    int qlen;           /* Query length */
    int nhits;          /* Number of hits for (current) query */
    int nok;            /* Number of OK hits for (current) query */
    int maxhit;         /* Maximum (OK) hit for (current) query */
    int asize;          /* Array size for hits */
    char *mask;         /* Mask of OK (qualifying) hits */
    char *hits;         /* Hit names[nhits] */
    char *scos;         /* Hit score records[nhits] */ 
    char *idens;        /* Hit identity records[nhits] */ 
    char *qseqs;        /* Query hit sequences[nhits] */
    int *qhsc;          /* Query hit starting coords */
    int *qhec;          /* Query hit ending coords */
    char *sseqs;        /* Subject hit sequences[nhits] */
    int *shsc;          /* Subject hit starting coords */
    int *shec;          /* Subject hit ending coords */
    int ihist[HISTDIM]; /* Identity histogram array */
    int *alens;         /* Alignment hit lengths */
    int *amats;         /* Alignment match counts */
    int *amms;          /* Alignment mismatch counts */
    int *agaps;         /* Alignment gap counts */
    int *hlens;         /* Hit match lengths */
    int *hnums;         /* Hit matching base numbers */
    REAL *hfmb;         /* Hit matching base fraction */
}BLASTANS;

#define CHECK_BLASTOUT(wf) if(wf){ DestroyBlastoutI(wf); wf=NULL; }
#define CHECK_BLASTANS(wf) if(wf){ DestroyBlastansI(wf); wf=NULL; }


/***
*   IO format codes 
*/
#define BOT_NCBI    100     /* NCBI Blast+ output format; New */
#define BOT_NCBOLD  101     /* NCBI Blast output format; Old */
#define BOT_WU      102     /* WashU Blast output format */
#define BOT_SMWU    103     /* Smart WashU Blast output format */

#define HEAD_CHECK  10000   /* Max lines to guess type */

#define BOT_PERQ    150     /* Per-query output */
#define BOT_PERHIT  155     /* Per-hit output */

#define DEF_MAXHITS 2500    /* ef max hits / record */

#define NAME_MAX        50
#define QNAME_FORM_S    "%s%-30s\t%s\n"
#define SNAME_FORM_S    "%s%-30s\t%s\n"

/*********************** ppp ********************
* C function listing generated by gen_prot
* Wed May 30 09:33:00 2018
*/
/****************************************************************
* blast_io.c
*/
int LoadBlastoutRecordI(BLASTANS *aPO,int warn);
int LoadBlastRecordI(BLASTANS *aPO,FILE *inPF,int warn);
int LoadBlastQueryLineData(BLASTANS *aPO, char *inbufS, FILE *inPF);
int BlastRecordEndLine(BLASTANS *aPO, char *bufS);
int AddBlastHitSetI(char *headS, FILE *inPF, BLASTANS *aPO, int hit, int warn);
int SaveOneHitRecI(BLASTANS *aPO, int hit, char *hitS, char *scoreS, char *identS,
    char *qseqS, char *sseqS, int *coordPI);
void DumpBlastout(BLASTOUT *bPO, BLASTANS *aPO, int qn, FILE *outPF);
void DumpBlastoutHitSeqs(BLASTOUT *bPO, BLASTANS *aPO, int hit, int mod_case, FILE *outPF);
void ModSeqAlignCase(char *qseqS, char *sseqS, int mod_case);
void ModSeqBrangeCase(BLASTOUT *bPO, BLASTANS *aPO, int hit, char *qseqS, char *sseqS);
void FillQuerySubjAlignString(char *qseqS, char *sseqS, char *alignS);
void DumpNameSeq(char *nameS, char *fmtS, char *sepS, char *seqS, int newline, FILE *outPF);
int FillBlastRecHitNameI(BLASTANS *aPO, int r, char *bufS);
int FillBlastRecHitCoordsLineI(BLASTANS *aPO, int r, char *bufS);
int FillBlastRecHitStatLineI(BLASTANS *aPO, int r, char *bufS);
void DumpBlastFracMatchList(BLASTOUT *bPO,BLASTANS *aPO,FILE *outPF);
void DumpBlastMatchBaseCounts(BLASTOUT *bPO,BLASTANS *aPO,FILE *outPF);
int GuessBlastInputTypeI(FILE *inPF);
void FillBlastFormatString(int type,char *nameS);
int HandleBlastoutOutfileI(BLASTOUT *bPO,BLASTANS *aPO,int q,int head);
int SetBlastansFileI(char *inS,BLASTANS *aPO);

/****************************************************************
* blasthit.c
*/
void ProcHitList(BLASTOUT *bPO, BLASTANS *aPO);
int FillHitHistI(BLASTOUT *bPO, BLASTANS *aPO);
void IntegrateHist(BLASTANS *aPO);
int HitMatchCountI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int *denPI);
int Cont3pEndMatchI(BLASTANS *aPO, int len, char *qPC, char *sPC);
int FillHitAlignQcoordMapI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int *mapPI);
int FillHitAlignQseqMaskI(BLASTOUT *bPO, BLASTANS *aPO, int hit, char *maskPC);
int SetAlignCountsI(BLASTOUT *bPO, BLASTANS *aPO, int hit);
int AdjustHitLenI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int len);

/****************************************************************
* blastout.c
*/
int main(int argc, char **argv);
void BlastOutUse(void);
int BlastOutI(int argc, char **argv);
int OkBlastoutOptsI(BLASTOUT *bPO,int max);
void WriteHeader(BLASTOUT *bPO, BLASTANS *aPO, FILE *outPF);

/****************************************************************
* blaststr.c
*/
BLASTOUT *CreateBlastoutPO(void);
int DestroyBlastoutI(BLASTOUT *bPO);
int AddHitSpaceI(int num,BLASTANS **aPPO);
BLASTANS *CreateBlastansPO(int max);
int DestroyBlastansI(BLASTANS *aPO);
void InitBlastout(BLASTOUT *bPO);
void InitBlastans(BLASTANS *aPO, int full);
int AdjustHitLenI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int len);

