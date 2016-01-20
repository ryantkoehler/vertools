/*
* dna.h
*
* Copyright 2016 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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


#ifndef __DNAH__
#define __DNAH__


#define MAX_PSLEN       500 /* Max probset amplicon seq length */
#define MAX_IPS         5   /* Max probset internal priming sites */

#define SEQ_ID          4052
#define SEQSET_ID       4054
#define SEQCOMP_ID      4056

#define SSIZE_MIN       100     /* String size (allocated) min size */
#define SSN_BLOCK       100     /* SEQSET nsize (allocated) block increment */
#define SEQ_BLOCK       1000    /* SEQ allocate block increment */

typedef struct SEQ      /* Sequence datatype */
{
    int ID;
    int ind;            /* Index if in set */
    struct SEQSET *par; /* Possible parent sequence set */
    char name[NSIZE];   /* Name string */
    char *seq;          /* Sequence string */
    int ssize;          /* Allocated sequence string length */
    int len;            /* Sequence length */
    int nsnp;           /* Number of SNP records */
    int flag;           /* Sequence bitfield */
}SEQ;

typedef struct SEQSET   /* Sequence set datatype */
{
    int ID;
    int type;           /* Type flag */
    int ambigs;         /* Flag for inclusion of any ambiguous bases */
    char name[NSIZE];   /* Name of whole set */
    char source[NSIZE]; /* Filename */
    int n,nsize;        /* Number of seqs, Allocated array length */
    struct SEQ **seqs;  /* Array of SEQ structures */
    char *mask;         /* Sequence masking field [n] */
}SEQSET;

typedef struct SEQCOMP
{
    int ID;
    int slen,nbase;         /* Seq length, number of non-ambig bases */
    int ra,rc,rg,rt;        /* Maximum base rows */
    int rs,rw,rr,ry,rk,rm;  /* Max degenerate base rows */
    int na,nc,ng,nt;        /* Number of bases */
    int dinuc[16];          /* inucleotide numbers */
    int n_dinuc;            /* umber of dinucleotides */
    REAL fa,fc,fg,ft;       /* raction of bases */
}SEQCOMP;

#define CHECK_SEQ(ob)       if(ob){DestroySeqI(ob);ob=NULL;}
#define CHECK_SEQSET(ob)    if(ob){DestroySeqsetI(ob);ob=NULL;}
#define CHECK_SEQCOMP(ps) if(ps){DestroySeqcompI(ps);ps=NULL;} 

/***
*   Bit field flags for seqs
*/
#define SEQ_SNP (1<<1)      /* Seq has explicit SNP(s) */
#define SEQ_AMB (1<<2)      /* Seq has ambiguous bases; N */
#define SEQ_DEG (1<<3)      /* Seq has degenerate bases; SW,RY... */
#define SEQ_LC  (1<<4)      /* Seq has lower case */
#define SEQ_NS  (1<<5)      /* Seq has non-standard bases; X,L... */
#define SEQ_INS (1<<6)      /* Seq has insertion (-) indicated */
#define SEQ_DEL (1<<7)      /* Seq has deletion (*) indicated */
#define SEQ_COAX (1<<8)     /* Seq has coaxial break (/) indicated */
#define SEQ_HASU (1<<9)     /* Seq has U (i.e. not T) */

#define IS_SEQ_SNP(v)   ((v) & SEQ_SNP)
#define IS_SEQ_AMB(v)   ((v) & SEQ_AMB)
#define IS_SEQ_DEG(v)   ((v) & SEQ_DEG)
#define IS_SEQ_LC(v)    ((v) & SEQ_LC)
#define IS_SEQ_NS(v)    ((v) & SEQ_NS)
#define IS_SEQ_INS(v)   ((v) & SEQ_INS)
#define IS_SEQ_DEL(v)   ((v) & SEQ_DEL)
#define IS_SEQ_COAX(v)  ((v) & SEQ_COAX)
#define IS_SEQ_HASU(v)  ((v) & SEQ_HASU)

#define SET_SEQ_SNP(v)  ((v) |= SEQ_SNP)
#define SET_SEQ_AMB(v)  ((v) |= SEQ_AMB)
#define SET_SEQ_DEG(v)  ((v) |= SEQ_DEG)
#define SET_SEQ_LC(v)   ((v) |= SEQ_LC)
#define SET_SEQ_NS(v)   ((v) |= SEQ_NS)
#define SET_SEQ_INS(v)  ((v) |= SEQ_INS)
#define SET_SEQ_DEL(v)  ((v) |= SEQ_DEL)
#define SET_SEQ_COAX(v) ((v) |= SEQ_COAX)
#define SET_SEQ_HASU(v) ((v) |= SEQ_HASU)

#define OFF_SEQ_SNP(v)  ((v) &= ~SEQ_SNP)
#define OFF_SEQ_AMB(v)  ((v) &= ~SEQ_AMB)
#define OFF_SEQ_DEG(v)  ((v) &= ~SEQ_DEG)
#define OFF_SEQ_LC(v)   ((v) &= ~SEQ_LC)
#define OFF_SEQ_NS(v)   ((v) &= ~SEQ_NS)
#define OFF_SEQ_INS(v)  ((v) &= ~SEQ_INS)
#define OFF_SEQ_DEL(v)  ((v) &= ~SEQ_DEL)
#define OFF_SEQ_COAX(v) ((v) &= ~SEQ_COAX)
#define OFF_SEQ_HASU(v) ((v) &= ~SEQ_HASU)

/***
*   Direction / join codes
*/
#define FORWARD 1
#define REVERSE -1
#define APPEND  1
#define PREPEND -1

/***
*   Alphabet to number mapping
*/
#define DNA_A   0
#define DNA_C   1
#define DNA_G   2
#define DNA_T   3
#define DNA_N   4   /* Ambigous */
#define DNA_O   5   /* Other */
#define DNA_NUM 6

/***
*   File format codes
*/
#define SEQFM_RAW       81  /* Raw format; <name> <seq> on one line */
#define SEQFM_REX       82  /* Raw with extras; <name> <seq> <???>... */
#define SEQFM_FASTA     85  /* Fasta format; > on one line, followed by seq */
#define SEQFM_NFAS      86  /* "nice" Fasta format (seq wrapped) */
#define SEQFM_SEQ       88  /* Sort of like gen-bank? */

/***
*   Formating stuff
*/
#define COMP_S          "_CP"   /* ppended name for compliment seqs */
#define NICE_BLEN       10      /* number of chars / block in "nice" output */
#define NICE_LLEN       50      /* number of chars / line in "nice" output */
#define RAW_PFORM_S     "%-15s\t%s\n"   /* raw file output format string */

#define DEF_SEQ_NAME_S  "No-Name"

/***
*   Sequence "cleaning" levels; should be low < mid < hi
*/
#define SCLEAN_LOW      90  /* Sequence cleaning to any letter */
#define SCLEAN_MID      91  /* Sequence cleaning to ACGTN + IUPAC degens */
#define SCLEAN_HI       92  /* Sequence cleaning to ACGTN */

/***
*   Sequence error codes
*/
#define EC_SEQ_LEN      230     /* eq length */
#define EC_SEQ_BCHAR    231     /* eq bad character */
#define EC_SEQ_WS       232     /* eq with space */
#define EC_SEQ_BSNP     235     /* NP problem */



/*********************** ppp ********************
* C function listing generated by gen_prot
* Fri Jun 27 13:08:46 2014
*/
/****************************************************************
* dna.c
*/
SEQ *CreateSeqPO(int len, char *seqS, char *nameS);
int DestroySeqI(SEQ *seqPO);
int AdjustSeqSizeI(SEQ *seqPO,int len,int error);
void DumpSeq(SEQ *seqPO,FILE *outPF);
void InitSeq(SEQ *seqPO,int par,int mem);
int CopySeqI(SEQ *fseqPO, SEQ *sseqPO, int st, int len);
int NarrowSeqI(SEQ *seqPO,int start,int len,int dir,int fits);
SEQSET *CreateSeqsetPO(int n);
int DestroySeqsetI(SEQSET *ssPO);
int AdjustSeqsetSizeI(SEQSET *ssPO,int size);
void DumpSeqset(SEQSET *ssPO,int all,FILE *outPF);
int AddNamedSequenceToSeqsetI(char *seqS, char *nameS, SEQSET *ssPO,
    int *indPI);
int AddSeqToSeqsetI(SEQ *seqPO,SEQSET *ssPO);
void InitSeqset(SEQSET *ssPO);
void SeqsetUnmaskDims(SEQSET *ssPO,int *nPI,int *sPI,int *mPI, int *snPI);
int FinishSeqSettingsI(SEQ *seqPO,int clean,int error);
SEQCOMP *CreateSeqcompPO(void);
int DestroySeqcompI(SEQCOMP *scPO);
void InitSeqcomp(SEQCOMP *scPO);
void DumpSeqcomp(SEQCOMP *scPO,FILE *outPF);

/****************************************************************
* dna_char.c
*/
int GoodDNABaseI(char c);
int GoodDNACompBasesI(char c,char s);
int DNABaseIndexI(char seqC);
char DNAIndexBaseC(int ind);
int SeqPairIndexI(char *seqS);
int SeqBasePairIndexI(char c1, char c2);
char CompDNABaseC(char c);
int CompDNASeqI(char *seqS,int len,char *compS);
int InverseDNASeqI(char *seqS,int len,char *compS);
int ReverseDNASeqI(char *seqS,int len,char *compS);
int RandomDNASeqI(char *seqS, int len, int *fracsPI);
int RandTweakDNASeqI(char *seqS,int len,int tweak);
int CountNonStandBasesI(char *seqS,int len);
int ExtractForOrRevSubSeqI(char *seqS,int slen,int st,int en,int dir,char *subS);
int GetReducedSeqI(char *seqS,int len,int fir,int las,int rre,char *newS);
int CleanUpSeqI(char *inS,int slen,char *outS,int ols,int mlc);

/****************************************************************
* dna_comp.c
*/
int GetSeqCompositionI(char *seqS, int slen, int clean, SEQCOMP *scPO);

/****************************************************************
* dna_file.c
*/
int FigureSeqFileTypeI(int iraw, int ifas, char *fnameS, int error);
int GuessSeqFileTypeI(char *fnameS, int error);
int ParseSeqTypeI(char *extS,int exact);
void FillSeqFtypeExtString(int type,char *typeS);
void FillSeqFtypeDescString(int type,char *typeS);
int OkSeqInFormatI(int type,char *nameS,int error);

/****************************************************************
* dna_in.c
*/
int GuessAndGetSeqsetI(char *fnameS, SEQSET **ssPPO,int clean,int error);
int ReadInSeqsetI(char *fnameS, int type, int clean, SEQSET **ssPPO, int error);
SEQSET *GetSeqsetPO(FILE *fPF,int type,int clean,int error);
int ParseSeqI(FILE *inPF,int iform,int clean,int error,SEQ *seqPO);
int ParseRawSeqI(FILE *fPF,int error,SEQ *seqPO);
int ParseFastaSeqI(FILE *fPF,int error,SEQ *seqPO);
int TrimSeqTrailingCharsI(SEQ *seqPO, int what);

/****************************************************************
* dna_out.c
*/
void WriteSeqset(SEQSET *ssPO,int oform,int head,FILE *outPF);
void WriteSeqsetHeader(SEQSET *ssPO,char *snameS,char *formS,FILE *outPF);
int WriteSeqsetSeqI(SEQSET *ssPO,int ind,int oform,FILE *outPF);
void WriteSeq(SEQ *seqPO,int oform,FILE *outPF);
void WriteMultiLineSeq(char *seqS,int slen, int nice, int line, int block, FILE *outPF);
void WriteNiceSeqLine(char *preS,char *seqS,int len,FILE *outPF);
void WriteSeqsetNameList(SEQSET *ssPO,FILE *oPF);

/****************************************************************
* dna_seqs.c
*/
int GetSeqsetSeqI(SEQSET *ssPO,int ind, SEQ **seqPPO);
int GetSeqsetSeqStringI(SEQSET *ssPO,int ind, char **seqPPC);
int FillSeqsetSeqStringI(SEQSET *ssPO,int ind,char *seqS,int max);
int SetSeqsetSeqNameI(SEQSET *ssPO,int ind, char *nS, int nlen);
int GetSeqsetSeqNameI(SEQSET *ssPO,int ind, char **namePPC);
int FillSeqsetSeqNameI(SEQSET *ssPO,int ind,char *nS,int max);
int GetSeqsetSeqLenI(SEQSET *ssPO,int ind);
int GetSeqsetSeqSnpCountI(SEQSET *ssPO,int ind);
int GetSeqsetNumI(SEQSET *ssPO);
int SeqsetMinLenI(SEQSET *ssPO);
int SeqsetMaxLenI(SEQSET *ssPO);
int SeqsetMinMaxLenI(SEQSET *ssPO,int *minPI,int *maxPI);
void SetSeqsetName(SEQSET *ssPO,char *nameS);
int FillSeqsetNameStringI(SEQSET *ssPO,char *nameS,int max);
void SetSeqsetSource(SEQSET *ssPO,char *nameS);
int FillSeqsetSourceStringI(SEQSET *ssPO,char *nameS,int max);
int FindNamedSeqInSeqsetI(SEQSET *ssPO, char *nameS, int kc, char *tPC,
    SEQ **seqPPO);
int SetSeqSequenceI(SEQ *seqPO,char *seqS,int len);
int AppendSeqSequenceI(SEQ *seqPO,char *seqS,int len);
int AppendSeqCharI(SEQ *seqPO, char c, int error);
void SetSeqName(SEQ *seqPO,char *nameS);
int GetSeqLenI(SEQ *seqPO);
int GetSeqSnpCountI(SEQ *seqPO);
int AnySeqAmbigsI(SEQ *seqPO);
int GetSeqSeqI(SEQ *seqPO, char **seqPPC);
int FillSeqSeqStringI(SEQ *seqPO,char *seqS,int max);
int FillSeqNameStringI(SEQ *seqPO,char *nameS,int max);
int UppercaseSeqSeqI(SEQ *seqPO);
int UppercaseSeqsetSeqsI(SEQSET *ssPO);
int ReverseCompSeqSeqI(SEQ *seqPO);
int ReverseCompSeqsetSeqsI(SEQSET *ssPO);

/****************************************************************
* dna_snp.c
*/
int ExpandSeqSnpI(SEQ *seqPO, int snp, int al, char *seqS, int up, int max,
    int *startPI, int *targPI, char *alS, int *asizePI);
int GetSeqsetSnpWindowI(SEQSET *seqsPO, int ind, int snp, int up, int down, 
    int max, char *seqS, int *sposPI);
int GetSeqSnpWindowI(SEQ *seqPO, int snp, int up, int down, int max, 
    char *seqS, int *sposPI);
int GetSeqSnpBlockI(SEQ *seqPO, int snp, char *seqS, int max);
int GetSeqSnpAlleleI(SEQ *seqPO, int snp, int al, char *seqS, int max);
int CleanUpSnpseqLineI(SEQ *seqPO);
int SnpHasDelimI(char *seqS,int max);
int ExpandSeqSingBaseSNPsI(SEQ *seqPO);

/****************************************************************
* dna_splc.c
*/
int SpliceTwoSeqsI(char *firS, char *secS, char *upS, int dir, 
    int over, int ext, char *inS, char *newS);
int NeedExtensionI(char *seqS,int ov, char *upS, int dir, int ext);
int GetExtensionSeqI(char *upS,int dir,int ext,char *inS);
int ProbeTargOverlapI(char *probS,int plen,char *targS,int tlen,int dir);
int MergeSeqsI(SEQ *seqPO,SEQ *addPO,int how);

/****************************************************************
* snp_char.c
*/
int CollapseSnpStringI(char *snpS,char *scS);
int ExpandDegBaseI(char bC,char *snpS);
int BaseDegeneracyI(char bC);
int IsSnpSeqExpandedI(char *snpS);
int IsDegenBaseI(char bC, char dC);
int CountSeqAmbigsI(char *sS,int s,int e);
int CountSeqAmbigDegensI(char *sS,int s,int e,int type);
int CountSeqSnpSitesI(char *sS,int s,int e);


#endif
