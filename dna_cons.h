/*
* dna_cons.h
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

 

#define IN_CONS_MAX     500 /* Max length to be evaluated with IN_CONS */

/***
*   Data structure for Intrinsic Constraints
*/
#define IN_CONS_ID  6015

typedef struct IN_CONS
{
    int ID;
    /***
    *   Sequence intrinsic constraints
    */
    char parname[100];  /* Parameter file name */
    int maxlen;         /* Sequence length */
    int cont_con;       /* Any content constraints? */
    REAL minAc,maxAc;   /* Min / max A content */
    REAL minCc,maxCc;   /* Min / max C content */
    REAL minGc,maxGc;   /* Min / max G content */
    REAL minTc,maxTc;   /* Min / max T content */
    REAL minSc,maxSc;   /* Min / max S (G,C) content */
    REAL minWc,maxWc;   /* Min / max W (A,T) content */
    REAL minRc,maxRc;   /* Min / max R (A,G) content */
    REAL minYc,maxYc;   /* Min / max Y (C,T) content */
    REAL minKc,maxKc;   /* Min / max K (G,T) content */
    REAL minMc,maxMc;   /* Min / max M (A,C) content */
    int row_con;        /* Any row constraints? */
    int max_rowA;       /* Max A in a row */
    int max_rowC;       /* Max C in a row */
    int max_rowG;       /* Max G in a row */
    int max_rowT;       /* Max T in a row */
    int max_rowS;       /* Max S (G,C) in a row */
    int max_rowW;       /* Max "W"eak (A,T) in a row */
    int max_rowR;       /* Max R (purine A,G) in a row */
    int max_rowY;       /* Max Y (pyridine C,T) in a row */
    int max_rowM;       /* Max M (amine A,C) in a row */
    int max_rowK;       /* Max K (ketone G,T) in a row */
    int rep_con;        /* Any repeat constraints? */
    int max_repRY;      /* Max RY repeats in a row */
    int max_repYR;      /* Max YR repeats in a row */
    int max_repMK;      /* Max MK repeates in a row */
    int max_repKM;      /* Max KM repeats in a row */
    /***
    *   Run time constraint vars 
    */
    int max_A;              /* Max A count */
    int max_C;              /* Max C count */
    int max_G;              /* Max G count */
    int max_T;              /* Max T count */
    int max_S,max_W;        /* Max S (G,C) W (A,T) count */
    int max_R,max_Y;        /* Max R (A,G) Y (C,T) count */
    int max_K,max_M;        /* Max K (G,T) M (A,C) count */
    int min_A[IN_CONS_MAX]; /* Min A count as a function of length */
    int min_C[IN_CONS_MAX]; /* Min C count as a function of length */
    int min_G[IN_CONS_MAX]; /* Min G count as a function of length */
    int min_T[IN_CONS_MAX]; /* Min T count as a function of length */
    int min_S[IN_CONS_MAX]; /* Min S count as a function of length */
    int min_W[IN_CONS_MAX]; /* Min W count as a function of length */
    int min_R[IN_CONS_MAX]; /* Min R count as a function of length */
    int min_Y[IN_CONS_MAX]; /* Min Y count as a function of length */
    int min_K[IN_CONS_MAX]; /* Min K count as a function of length */
    int min_M[IN_CONS_MAX]; /* Min M count as a function of length */
}IN_CONS;

#define CHECK_IN_CONS(ic) if(ic){DestroyInConsI((IN_CONS *)ic); ic=NULL;    }



/*********************** ppp ********************
* C function listing generated by gen_prot
* Fri Feb 23 11:12:26 2001
*/
/****************************************************************
* dna_cons.c
*/
IN_CONS *CreateInConsPO(void);
int DestroyInConsI(IN_CONS *iconPO);
void InitInCons(IN_CONS *iconPO);
int LoadInConsParsI(char *fnameS, IN_CONS *pxPO);
int ParseInConsParsI(FILE *fPF,IN_CONS *iconPO);
int PrepareInConsI(IN_CONS *iconPO);
int ConsistInConsI(IN_CONS *iconPO);
int ConsistPercentParsI(int len,REAL minR,REAL maxR,char *whatS);
void DumpInCons(IN_CONS *iconPO,FILE *oPF);
int ReportMinMaxCompI(char *nameS,REAL minR,REAL maxR,FILE *oPF);
int ReportMaxRowI(char *nameS,int max,int mlen,FILE *oPF);
int SeqInConsOkI(char *seqS, int len, int full, IN_CONS *iconPO);
int FullSeqInConsOkI(char *seqS,int len,int full,IN_CONS *iconPO);
int SeqContConsOkI(char *seqS, int len, int full, IN_CONS *iconPO);
int SeqRowConsOkI(char *seqS, int len, IN_CONS *iconPO);
int SetInConsMaxlenI(IN_CONS *iconPO, int len);
int GetInConsMaxlenI(IN_CONS *iconPO);

