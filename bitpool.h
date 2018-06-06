/*
* bitpool.h
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

#ifndef __BUTUITLH__
#define __BUTUITLH__
#include "wordlist.h"


/***
*   Type for bitwise operations; int32_t is 32 bits
*/
typedef int32_t BITPTR;
#define BU_WORDLEN      32      /* Word length; 32 bits */
#define BU_BLK_SIZE     8       /* block size; Process 8 bits at a time */


#define BITPOOL_ID  4002

typedef struct BITPOOL
{
    int ID;
    char name[NSIZE];       /* Name */
    char fname[NSIZE];      /* Input file name */
    int num;                /* Number of records */
    BITPTR *blocks;         /* Bit array */
    int n_blocks;           /* Dimension of blocks */
    int bsize;              /* Bitstring size for records */
    int bkprec;             /* Blocks per record */
    struct WORDLIST *mnames;/* Member names */
    char prlform[NSIZE];    /* Print formatting string for labels */
    DOUB coef;              /* Coefficient for weighting */
}BITPOOL;

#define CHECK_BITPOOL(ob) if(ob){DestroyBitpoolI(ob); ob=NULL;}


#define BP_NSIZE            100     /* ize for bit pool member name space */

#define DEF_BP_PRLFORM_S    "%s"    /* efault row label format string */


/***
*   Input format codes
*/
#define BPF_BITS    200     /* BitString input format as 0100110.. etc */
#define BPF_SBITS   201     /* BitString input format as 0 1 0 0 1 .. etc */
#define BPF_HEX     202     /* BitString input format as hex numbers */
#define BPF_TAB     203     /* BitString input format as table */

/***
*   Pairwise operator codes
*/
#define BIT_NOT     10     
#define BIT_AND     11
#define BIT_OR      12     
#define BIT_XOR     13     
#define BIT_ON      14
#define BIT_OFF     15



/*********************** ppp ********************
* C function listing generated by gen_prot
* Sat Apr 26 19:34:17 2014
*/
/****************************************************************
* bitpool.c
*/
BITPOOL *CreateBitpoolPO(int bsize, int num);
int DestroyBitpoolI(BITPOOL *bpPO);
int InitBitpoolI(BITPOOL *bpPO);
int SetBitpoolDimsI(BITPOOL *bpPO, int bsize, int num);
int GetBitpoolDimsI(BITPOOL *bpPO, int *bsizePI, int *numPI);
int HandleBitpoolSpaceI(BITPOOL *bpPO, int num);
void SetBitpoolCoef(BITPOOL *bpPO,DOUB coD);
int SetBitpoolNamesI(BITPOOL *bpPO, char *nameS, char *fnameS, int max);
int GetBitpoolNamesI(BITPOOL *bpPO, char *nameS, char *fnameS, int max);
int AutoBitpoolOutFormattingI(BITPOOL *bpPO);
int SetBitpoolPrintFormI(BITPOOL *bpPO, char *prlformS);
int GetBitpoolPrintFormI(BITPOOL *bpPO, char *prlformS);
int CopyThisBitpoolMemberI(BITPOOL *fbPO, int fm, BITPOOL *sbPO, int sm);
int AddThisBitpoolMemberI(BITPOOL *bpPO, int m, char *bitS, char *nameS);
int SetThisBitpoolMemberI(BITPOOL *bpPO, int m, char *bitS, char *nameS);
int SetThisBitpoolNameI(BITPOOL *bpPO, int m, char *nameS);
int GetThisBitpoolNameI(BITPOOL *bpPO, int m, char *nameS);
int SetThisBitpoolBitstringI(BITPOOL *bpPO, int m, char *bitS);
int SetThisBitpoolBitI(BITPOOL *bpPO, int m, int b, int on);
int GetThisBitpoolPtrBitsI(BITPOOL *bpPO, int m, BITPTR **bPPI, int *nPI);
int GetThisBitpoolBitI(BITPOOL *bpPO, int m, int b, int *onPI);
int GetThisBitpoolOnCountI(BITPOOL *bpPO, int m, int *onPI);
int BitpoolBitwiseOpCountI(BITPOOL *fbPO, int fm, BITPOOL *sbPO, int sm, int op, 
    int *onPI);
int ModThisBitpoolI(BITPOOL *bpPO, int st, int en, int op) ;
int ModThisBitpoolBitsI(BITPOOL *bpPO, int m, int st, int en, int op) ;
int SameBitpoolDimsI(BITPOOL *fbPO, BITPOOL *sbPO, int size, int num);
void DumpBitpool(BITPOOL *bpPO, int st, int en, char *preS, FILE *outPF);
void DumpBitpoolDescription(BITPOOL *bpPO, char *pS, FILE *outPF);
void DumpBitpoolMembers(BITPOOL *bpPO, int st, int en, char *pS, FILE *outPF);
int DumpThisBitpoolMemberI(BITPOOL *bpPO, int m, int oform, char *pS, FILE *outPF);
int FigureBitFileTypeI(int ibit, int isbit, int ihex, int itab, char *fnameS, int error);
int GuessBitFileTypeI(char *nameS, int error);
int OkBitFileFormatI(int form, char *nameS, int error);
int OkBitwiseOpI(int op, char *opS);
int GetBitpoolI(char *fileS, int iform, int error, BITPOOL **bpPPO);
int ParseBitpoolI(FILE *inPF, int iform, int error, BITPOOL **bpPPO);
int ParseBitpoolOneZeroI(FILE *inPF, int error, BITPOOL **bpPPO);
int BitStringOneZeroLineDimI(char *lineS, int *bsizePI);
int SetThisBitI(BITPTR *bitsPI,int bit,int on);
int GetThisBitI(BITPTR *bitsPI,int bit);
int ArrayToBitsI(int *aPI, int len, BITPTR *bitsPI);
int BitsToArrayI(BITPTR *bitsPI, int len, int clean, int *aPI);
int OnBitCountI(BITPTR *bitsPI, int len);
int DifBitCountI(BITPTR *bitsPI, BITPTR *sbitsPI, int len);
int BitwiseOpCountI(BITPTR *bitsPI, BITPTR *sbitsPI, int len, int op);
int DifByMaxBitsI(BITPTR *bitsPI, BITPTR *sbitsPI, int len, int max);
REAL BitTanCoefR(BITPTR *bitsPI, BITPTR *sbitsPI, int len);


#endif
