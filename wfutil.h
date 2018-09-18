/*
* wfutil.h
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


#include "dna.h"

#define DEF_WSIZE   3   /* Default word size */
#define ALPHDIM     4   /* Default alphabet dimension; ACGT */

#define WORDFREQ_ID     4092    
typedef struct WORDFREQ
{
    int ID;
    char name[NSIZE];   /* Name field */
    int size;           /* Word size */
    int ald;            /* Alphabet dimension */
    int pmax;           /* Indexing number = n/ald */
    int n;              /* Number of different words; freqs */
    struct FREC *freqs; /* Frequency count record */
}WORDFREQ;

typedef struct FREC
{
    int id; DOUB n;
}FREC;

#define CHECK_WORDFREQ(wf) if(wf){ DestroyWordfreqI(wf); wf=NULL; }

#define ADD_WF      71  /* Add for merge values */
#define SUB_WF      72  /* Subtract for merge */
#define DSR_WF      73  /* Diff-Sum-Ratio for merge */
#define FDSR_WF     74  /* Frequency-weighted-Diff-Sum-Ratio for merge */

#define SORT_HILO   81  /* Code to sort hi to low */
#define SORT_IND    83  /* Code to sort on indices */

/*********************** ppp ********************
* C function listing generated by gen_prot
* Fri Jun 10 07:01:29 2011
*/
/****************************************************************
* wfutil.c
*/
WORDFREQ *CreateWordfreqPO(int size,int ald);
int DestroyWordfreqI(WORDFREQ *wfPO);
void InitFrecIDs(FREC *frecsPO,int nw);
int CalcInDimI(int size,int aldim);
int TallyWordsI(SEQ *seqPO,WORDFREQ *wfPO,int step);
int DumpWordsI(WORDFREQ *wfPO, int dgen, DOUB loD, DOUB hiD, FILE *outPF);
void WordfreqSummary(WORDFREQ *wfPO, DOUB loD, DOUB hiD, FILE *outPF);
DOUB TotalWordsD(WORDFREQ *wfPO,DOUB loD,DOUB hiD);
int TotalDifWordsI(WORDFREQ *wfPO,DOUB loD,DOUB hiD);
int SeqFromIndexI(int ind,int pmax,int aldim,char *wordS);
int IndexFromSeqI(char *seqS,int len,int pmax,int aldim);
int FillWordfreqSeqStringI(WORDFREQ *wfPO, int ind, char *seqS);
int CompIndexI(int ind,int pmax,int aldim);
void CollapseDegenRecs(WORDFREQ *wfPO);
void SortWordFreqs(WORDFREQ *wfPO,int how);
int qSortFrecsI(const void *e1, const void *e2);
int qSortIndsI(const void *e1, const void *e2);
int GetWordFreqsI(char *nameS,WORDFREQ **wfPPO);
int LoadWordFreqDataI(FILE *inPF,WORDFREQ **wfPPO);
void NormalizeFrecs(WORDFREQ *wfPO);
void  LogFrecs(WORDFREQ *wfPO);
int CompatWordfreqsI(WORDFREQ *fPO, WORDFREQ *sPO,int verb);
int MergeWordfreqsI(WORDFREQ *fPO, WORDFREQ *sPO,int how, WORDFREQ *aPO);

