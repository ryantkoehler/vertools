/*
* prim.h
*
* Copyright 2015 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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


#define RTK_S   "Ryan Koehler, ryan@verdascend.com"
#define BD_S    "Build date"


/**************** TYPES ******************/
#define REAL double
#define DOUB double

/* eys to keep functions out of headers; i.e. keep them private */
#define PRIV_I  int
#define PRIV_V  void    

typedef unsigned long UNLG;
typedef unsigned int *OBJPTR;


/* ast macros */
#define CHAR(x)     (char)(x)
#define INT(x)      (int)(x)
#define SHORT(x)    (short int)(x)
#define LONG(x)     (long int) (x)
#define RNUM(x)     (REAL)(x)
#define DNUM(x)     (double)(x)

/* ype codes */
#define IS_CHAR     10      /* A character */
#define IS_INT      20      /* An integer */
#define IS_SHORT    21      /* A short integer */
#define IS_REAL     30      /* A "real" number */
#define IS_DOUB     32      /* A double precision number */

#define UPPER(x)    toupper((int)(x))

#define TRUE        1
#define FALSE       0
#define BOGUS       -12345666

#define OK          0   /* xit code; e.g. program is OK */
#define NOT_OK      1   /* xit code; e.g. program not OK */

#define IS_BOG(b)       ((b)==BOGUS)
#define ERR(fn,ms)      ErrorMsg(fn, __LINE__, __FILE__, ms);
#define BOG_CHECK(co)   if(co) ERR("BOG_CHECK","BOGOCITY DISCOVERED!")
#define VALIDATE(ob,ty) if((ob)==NULL) \
        { ValidObj(NULL,0,ty,__LINE__,__FILE__); } else \
        { ValidObj((OBJPTR)ob, ob->ID, ty, __LINE__, __FILE__); }
#define ALLOC(nm,sz)    AllocPO(nm,sz, __LINE__, __FILE__)
#define REALLOC(ob,nm,sz) ReAllocPO((UNLG *)(ob),nm,sz, __LINE__, __FILE__)
#define FREE(ob)        {FreeI((UNLG *)(ob), __LINE__, __FILE__); ob = NULL;}
#define CHECK_FREE(ob)  if((ob) != NULL)FREE(ob);
#define FILECLOSE(fl)   FileClose(fl);fl=NULL
#define CHECK_FILE(fl)  if(((fl)!=NULL)&&((fl)!=stdout)){ FILECLOSE(fl); }
#define NFILECLOSE(fl,nm)   NewFileClose(fl,nm);fl=NULL
#define CHECK_NFILE(fl,nm)  if(((fl)!=NULL)&&((fl)!=stdout))\
                                { NFILECLOSE(fl,nm);} 
#define HAND_NFILE(file) if(!(file)){file=stdout;}
#define NOT_YET         printf("\n---- NOT YET IMPLIMENTED ----\n%s (%d)\n\n",\
                        __FILE__,__LINE__)
#define SORT_ASCEND 1
#define SORT_DECEND -1
#define ALLOC_BLOCK     5000    /* llocate in blocks this size */


/**************** STRINGS **************** sss **/
#define LOWER_CASE      2
#define UPPER_CASE      3
#define CAP_CASE        4
#define INIT_S(st)      st[0] = '\0'
#define NO_S(st)        (st[0] == '\0')
#define DEF_BS          255         /* Default buffer size */
#define LINEGRAB        (DEF_BS-1)  
#define BBUFF_SIZE      2000        /* Big buffer size */
#define BLINEGRAB       (BBUFF_SIZE-1)  
#define NSIZE           1000        /* Name buffer size */
#define COM_LINE(bf)    ((bf[0] == '!')||(bf[0] == '#'))
#define EQSTRING(sa,sb,sz) (strncmp((sa),(sb),(sz))==0)
#define STDIN_STR(sa)   (strcmp((sa),"-")==0)
#define PASS_BLANK(st)  while( (*st != '\0') && (!isgraph(INT(*st))) )  \
                        {   st++;  }
#define PASS_WORD(st)   while( (*st != '\0') && (isgraph(INT(*st))) )       \
                        {   st++;  }
#define NEXT_WORD(st)   PASS_WORD(st) PASS_BLANK(st)

/** Replace isprint which has problems with tab */ 
#define ISLINE(c)       ( isgraph(INT(c)) || ((c)==' ') || ((c)=='\t') )

/**************** Standardized warnings / errors **/
#define WARN_LEV    1
#define PROB_LEV    2
#define ABORT_LEV   3

#define LINEBAR_S \
    "# ---------------------------------------------------------------------\n"
#define LINESEP_S \
    "# *********************************************************************\n"
#define NEWLINE         PrintI("\n")
#define BARLINE         PrintI(LINEBAR_S)
#define SEPLINE         PrintI(LINESEP_S)
#define WARNLINE        ErrorBanner(WARN_LEV,NULL)
#define PROBLINE        ErrorBanner(PROB_LEV,NULL)
#define ABORTLINE       ErrorBanner(ABORT_LEV,NULL)
#define ERRLINE         DB_PrI("\nERROR AT: %s, line %d\n",__FILE__,__LINE__)


/**************** MATH ******************/
#ifndef PI
#define PI              3.14159265358979323846
#endif
#define TOO_BIG         1234566600 
#define TOO_BIG_R       1234566600.1
#define TOO_BIG_D       1000000000000.1
#define TINY_R          0.0000000001
#define MAX_PRT_PREC    10                  /* ax precision for print */
#define DEGTORAD        (PI/180.000000)
#define RADTODEG        (180.000000/PI)
#define LOG_E_TEN_R     2.302585092994045901093613792909
#define LOG_E_TWO_R     0.693147180559945286226763982995
#define COS_R(an)       RNUM(cos((double)((an) * DEGTORAD)))
#define SIN_R(an)       RNUM(sin((double)((an) * DEGTORAD)))
#define TAN_R(an)       RNUM(tan((double)((an) * DEGTORAD)))
#define LOG(v)          RNUM(log(DNUM(v)))
#define LOG_10(v)       RNUM(log(DNUM(v))/LOG_E_TEN_R)
#define LOG_2(v)        RNUM(log(DNUM(v))/LOG_E_TWO_R)
#define CEIL(v)         RNUM(ceil(DNUM(v)))
#define FLOOR(v)        RNUM(floor(DNUM(v)))
#define SQRT_R(x)       RNUM(sqrt((double)(x)))

#define NUM_SIGN(r)     ( (RNUM(r) < 0.0) ? (-1) : (1) )
#define EQUAL(r,s,d)    (ABS_VAL(RNUM(r)-RNUM(s))<=RNUM(d))
#define ABS_VAL(r)      ( (RNUM(r) < 0.0) ? (-(r)) : (r) )
#define ROUND(x)        INT( ((x) > 0.0) ? (x) + 0.5 : (x) - 0.5)
#define PERCENT_R(f,w)  (RNUM(f)/RNUM(w)*100.0)
#define MAX_NUM(r,y)    (((r)>(y))?(r):(y))
#define MIN_NUM(r,y)    (((r)<(y))?(r):(y))
#define ODD_NUM(i)      ((int)(i)&1)
#define EVEN_NUM(i)     (!ODD_NUM(i))
#define LIMIT_NUM(vi,lo,hi) {vi=MAX_NUM(vi,lo);vi=MIN_NUM(vi,hi);}
#define BAD_I               INT(-TOO_BIG - 1)
#define BAD_R               RNUM(TOO_BIG_R * -1.000001)
#define BAD_D               RNUM(TOO_BIG_D * -1.000001)
#define BAD_INT(r)          (INT(r) < -TOO_BIG)
#define BAD_REAL(r)         (RNUM(r) < -TOO_BIG_R)
#define BAD_DOUB(r)         (DNUM(r) < -TOO_BIG_D)


/*********************** Globals ********************************* 
*   dbflags is a CHARACTER ARRAY, so the number flags must be < 256 
*/
#define NUM_DB_FLAGS 257
#define DB dbflagsGC


#ifdef __MAIN__

char dbflagsGC[NUM_DB_FLAGS];
int filecountGI = 0;
int numbytesGI = 0;
int numblockGI = 0;
int errlevelGI = 0;
FILE *outfileGPF = NULL;

#else

extern char dbflagsGC[NUM_DB_FLAGS];
extern int filecountGI;
extern int numbytesGI;
extern int numblockGI;
extern int errlevelGI;
extern FILE *outfileGPF;

#endif

#define __PRIMH__

/****************** Object IDs ************************************ iii ****/
/***
*   Geometry
*/
#define COORD_ID        50
#define MATRIX_ID       51
#define COORSYS_ID      52
/***
*   Sets, Lists, Histograms, grids,...
*/
#define SET_ID          60
#define LIST_ID         62
/***
*   Graphics things 
*/
#define COLOR_ID        75
#define PS_PARAM_ID     76
#define GIF_PARAM_ID    77
#define MOL_PLOT_PAR_ID 80


#define MATH_ADD    1
#define MATH_SUB    2
#define MATH_MUL    3
#define MATH_DIV    4
#define MATH_MIN    5
#define MATH_MAX    6


/**
*   Which random number functions to use?
*/
#define RAND        random
#define SRAND       srandom
#define RAND_DEN    2147483648.0
/* if rand() and srand()
*   define RAND_DEN (RAND_MAX + 1)
*/


/*********************** ppp ********************
* C function listing generated by gen_prot
* Thu May  8 15:58:23 2014
*/
/****************************************************************
* autil.c
*/
int InitArrayI(void *aPO, int vt, int start, int end, DOUB valD);
int ScaleArrayI(void *aPO,int vt, int start, int end, DOUB valD);
int ShiftArrayI(void *aPO,int vt, int start, int end, DOUB valD);
int NormArrayI(void *aPO,int vt, int start, int end);
int BoundArrayI(void *aPO, int vt, int start, int end, DOUB lowD, DOUB hiD);
int NumArrayValsI(void *aPO, int vt, int start, int end, DOUB lowD, DOUB hiD);
void CopyArray(void *fPO, void *sPO, int vt, int start, int end);
void SmoothArray(void *fPO, void *sPO, int vt, int n, int win);
void MixArrays(void *fPO, void*sPO, void *tPO, int vt, int n, int mix);
void DumpArray(void *aPO, int vt, int start, int end, char *formS, FILE *oPF);
void ArraySum(void *aPO, int vt, int start, int end, DOUB *sumPD);
int ArrayStatsI(void *aPO, int vt, int start, int end, DOUB *loPD, DOUB *hiPD, 
    DOUB *avPD, DOUB *sdPD);
int ArrayHistI(void *aPO, int vt, int start, int end, DOUB bsD, 
    DOUB lvD, DOUB hvD, int *ncPI, DOUB **hPPD);
int ArrayRandSequenceI(int *seq,int len,int max,int unique);
void InvertMask(char *maskPC,int len);
int MaskCountI(char *maskPC,int len);
int MaskRandSubsetI(char *maskPC, int len, DOUB fracR);
void SortArray(void *aPO, int vt, int n, int dir);

/****************************************************************
* dfutil.c
*/
int DataTableDimsI(FILE *inPF,int hr,int *ncolPI,int *nrowPI);
int FileStatsI(FILE *fPF, int *linePI, int *minPI, int *maxPI, int *comPI, 
    int *blankPI);
int FileLinesI(FILE *fPF, int ig_com, int ig_blank);
int EatOneLineI(FILE *fPF);

/****************************************************************
* mutil.c
*/
double PowD(double baseR, double expR);
void Srand(int seed);
int RandI(int range);
REAL RandR(REAL rangeR);
DOUB RandD(DOUB rangeD);
void FillRandSeedString(int seed, char *sS);
int SetRandSequenceI(int *seq,int len,int max,int unique);
DOUB LogFactorialD(int nI);
DOUB FactorialD(int nI);
DOUB SterlingLogFactorialD(int nI);
DOUB RandGaussD(REAL mD,REAL sD);

/****************************************************************
* putil.c
*/
int ParseArgsI(int argc, char **argv, char *formPC, ...);

/****************************************************************
* sutil.c
*/
void FillDateString(char *dateS);
void FillNumDateString(char *dateS);
int PrintTimeStoryI(char *stimeS,char *etimeS,int ve, char *preS);
int GetDateI(char *bufS);
int DaysI(int sec);
int HoursI(int sec);
int MinutesI(int sec);
int SecondsI(int sec);
void CleanString(char *stringS,int num);
int BlankStringI(char *lS);
void ReplaceChars(char sC, char *sS, char rC, char *rS);
void ReplaceSomeChars(char sC, char *sS, char rC, char *rS, int lim);
void RemoveChars(char *listS, char *sS, char *rS);
int RemoveThisCharI(char sC, char *sS, char *rS);
int RemoveSomeCharsI(char sC, char *sS, char *rS, int lim);
void PadString(char *sS,char pC,int len);
int TruncateStringI(char *bufS, char tC, char *newS);
void PrintString(char *sourceS, int len, FILE *outPF);
void Upperize(char *stringS);
void UpperizeToLen(char *stringS, int n);
void Lowerize(char *stringS);
void LowerizeToLen(char *stringS, int n) ;
int CountStringCaseI(char *sS, int len, int *lowPI, int *upPI);
void KillTrailStringBlanks(char *sS);
void Chomp(char *sS);
void StripFilePath(char *fnameS);
int GetFilePartsI(char *fnameS, char *pathS, char *baseS, char *extS);
int FnameFromPartsI(char *pathS, char *baseS, char *extS, char *fnameS);
int DelimSubStringI(char *inS, char sC, char eC, char *subS);
int GetNthWordI(char *bufS,int n,char *wordS);
int CopyNthWordI(char *bufS,int n,char *wordS,int max);
void ReverseString(char *inS,int len,char *outS);
int ParseTrueFalseI(char *wordS,int warn);
void FillTrueFalseString(int v,char *sS);
int FillOptionalNameStringI(char *nameS,char *deS);
void ReportParseErrorLine(char *bufS, char *funcS, char *storyS);
int SameStringsI(char *fS, char *sS, int kc) ;
void ReportWordMatchMods(char *prefS, int kc, int st, int sub, FILE *outPF);
int WordStringMatchI(char *fS, char *sS, int kc, int st, int sub) ;
void FloatFormatString(int w,int p,char *fmtS);
int NumPrecisionI(DOUB vD, int *nidPI, int *nfdPI, int *pwPI, char *pformS);
int DoubArrayPrecisionI(DOUB *valsPD, int n, int *pPI, int *wPI, char *formS);

/****************************************************************
* sysinfo.c
*/
void VersionSplash(FILE *outPF, char *verS, char *preS, int bars);
void GetSystemInfo(char *userS, char *hostS, char *osS, char *verS, 
    char *archS);
void PrintSysInfo(FILE *outPF,char *preS);
int DeRefSymLinkI(char *lnS, char *drS, int max);

/****************************************************************
* util.c
*/
void Init(int argc, char **argv);
void InitGvars(void);
int AllDoneI(int status,char *progS);
int CheckGvarsI(int err);
void DBInit(void);
void DBComlineCheck(int argc, char **argv);
void DBClear(void);
void DBReport(void);
void ErrorMsg(char *funS, int line, char *fileS, char *mesS);
void ErrorBanner(int level,FILE *oPF);
void ValidObj(OBJPTR obPO, int is, int type, int line, char *fileS);
int PrintI(char *formPS, ...);
int DB_PrI(char *formPS, ...);
void TimeStamp(char *preS, FILE *fPF);
FILE *FileOpenPF(char *nameS, char *modeS, int error);
void FileClose(FILE *fPF);
void NewFileClose(FILE *fPF,char *nameS);
FILE *OpenUFilePF(char *nameS, char *mS, char *fnameS);
int ExpandEnvVarI(char *inS, char *outS);
void *AllocPO(int num, int size,int line, char *fileS);
void *ReAllocPO(void *obPO, int num, int size, int line, char *fileS);
int FreeI(void *obPO, int line, char *fileS);
void ShowMemoryStatus(char *whatS);

