/*
* filepick.h
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

#define VERSION_S       "FilePick version 0.53"

#define DEF_STOK        "UID"

#define FILEPICK_ID     5061
typedef struct FILEPICK
{
    int ID;
    char inname[NSIZE];     /* Input file name */
    FILE *in;               /* Input file */
    char outname[NSIZE];    /* Output file name */
    FILE *out;              /* Outpuf file */
    char stok[NSIZE];       /* Start token for block delineation */
    char etok[NSIZE];       /* End token for block delineation */
    int first, last;        /* irst / last records to qualify */
    char wlisname[NSIZE];   /* ame for subset file */
    struct WORDLIST *wlis;  /* ubset names list structure */
    int do_ner;             /* lag to not dump end record line */
    int do_kc;              /* lag to keep case for name comparison */
    int do_wst;             /* lag to check start for word comparison */
    int do_wsub;            /* lag to check substring for word comparison */
    int do_not;             /* lag to invert record qualifcations */
    int do_dh;              /* lag to dump header (before first record) */
    int do_oh;              /* lag to dump only header */
    int do_snl;             /* lag to separate records with new lines */
    int do_srn;             /* lag to separate records with record number line */
    int do_rsany;
    int do_esf;             /* Flag to extract separate files */
    char esbase[NSIZE];     /* Extracted file base name */
    char esext[NSIZE];      /* Extracted file extension */
}FILEPICK;

#define CHECK_FILEPICK(fp)  if(fp){DestroyFilepickI(fp);fp=NULL;}

#define REC_SEP_S   "----------------------------"


/*********************** ppp ********************
* C function listing generated by gen_prot
* Sat Feb 22 10:10:21 2014
*/
/****************************************************************
* filepick.c
*/
int main(int argc, char **argv);
void FilePickUse(void);
int FilePickI(int argc, char **argv);
FILEPICK *CreateFilepickPO(void);
int DestroyFilepickI(FILEPICK *fpPO);
void InitFilepick(FILEPICK *fpPO);
int CheckFilepickOptions(FILEPICK *fpPO) ;
int LineMatchRecTokenI(FILEPICK *fpPO, char *bufS, char *tokS, char *nameS);
int MatchingTokenI(char *fS,char *sS,int kc);
int IsRecordOkI(FILEPICK *fpPO,int record,char *curecS);
int FpOutputThisLineI(FILEPICK *fpPO,int record,int in_rec,int ok_rec);
int HandleFilepickOutfileI(FILEPICK *fpPO,int record,char *recS);
int FpOutOneLineI(FILEPICK *fpPO, char *lineS, FILE *outPF);

