/*
* table.h
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


#ifndef __TABLEH__
#define __TABLEH__

#include "numlist.h"
#include "wordlist.h"

#define TABLE_ID    4062
typedef struct TABLE
{
    int ID;
    char name[NSIZE+1];     /* able name */
    char fname[NSIZE+1];    /* ile name */
    char prefix[NSIZE];     /* rint line prefix string */
    char pvform[NSIZE];     /* rint formatting string for values */
    char pvform2[NSIZE];    /* rint formatting string for values, more precision */
    char prlform[NSIZE];    /* rint formatting string for labels */
    char pvsep[NSIZE];      /* rint value separator */
    int nrow, ncol;         /* umber of rows, cols */
    struct NUMLIST *vals;   /* alue field (the DATA!)  */
    struct WORDLIST *clabs; /* olumn lables */
    struct WORDLIST *rlabs; /* ow lables */
    int rlab,clab;          /* lags for row / col lables */
    char corner[NSIZE];     /* ame for 'corner' between row and col labs */
    char *cmask;            /* ol mask */
    char *rmask;            /* ow mask */
    int prenum;             /* recision number of digits; e.g. 3 = 1.23400000 */
    DOUB prefac;            /* recision working factor */
}TABLE;

#define CHECK_TABLE(wf) if(wf){ DestroyTableI(wf); wf=NULL; }


/*
define TABBUFF_SIZE 10000   
*/
#define TABBUFF_SIZE        200000  /* arsing buffer size */
#define DEF_TAB_PVFORM_S    "%6.2f" /* efault value-dump format string */
#define DEF_TAB_PVFORM2_S   "%6.4f" /* efault value-dump format string */
#define DEF_TAB_PRLFORM_S   "%s"    /* efault row label format string */
#define DEF_TAB_PVSEP_S     " "     /* efault value-dump delimiter string */

#define DEF_TAB_PRENUM      4       /* efault value precision (digits) */
#define MAX_TAB_PRENUM      12      /* ax value precision */
#define MIN_TAB_PRENUM      0       /* in value precision */

#define TABLE_ROW           333     /* ode to indicate row(s) */
#define TABLE_COL           334     /* ode to indicate col(s) */
#define TABLE_FULL          335     /* ode to indicate full table */



/*********************** ppp ********************
* C function listing generated by gen_prot
* Wed Apr  4 20:38:01 2018
*/
/****************************************************************
* table.c
*/
int SetTablePrecisionI(TABLE *tabPO, int prenum);
int SetTableValPrecisionI(TABLE *tabPO, DOUB valD, DOUB *newPD);
int SameDimTablesI(TABLE *tabPO, TABLE *stabPO, int mask,int error);
int SameDimTableRowsI(TABLE *tabPO, TABLE *stabPO, int mask,int error);
int SameDimTableColsI(TABLE *tabPO, TABLE *stabPO, int mask,int error);
int SetTableRowMaskI(TABLE *tabPO,int row,int what);
int SetTableColMaskI(TABLE *tabPO,int col,int what);
int SetTableRowRangeMaskI(TABLE *tabPO,int srow, int erow, int what);
int SetTableColRangeMaskI(TABLE *tabPO,int scol, int ecol, int what);
int GetTableRowMaskI(TABLE *tabPO,int row);
int GetTableColMaskI(TABLE *tabPO,int col);
int GetTableDimsI(TABLE *tabPO, int *rowPI, int *colPI, int mask);
int GetTableRowsI(TABLE *tabPO,int mask);
int GetTableColsI(TABLE *tabPO,int mask);
int SetTableNamesI(TABLE *tabPO, char *nameS, char *fnameS, int max);
int GetTableNamesI(TABLE *tabPO, char *nameS, char *fnameS, int max);
int FindTableNamedRowI(TABLE *tabPO, int off, char *nameS, int kc, int nc);
int FindTableNamedColI(TABLE *tabPO, int off, char *nameS, int kc, int nc);
int AutoTableOutFormattingI(TABLE *tabPO, int val, int row);
int SetTablePrintformI(TABLE *tabPO, char *valS, char *val2S, char *sepS, 
    char *rowS, char *preS);
int GetTablePrintformI(TABLE *tabPO, char *valS, char *val2S, char *sepS, 
    char *rowS, char *preS);
int IsTableCSVFormatI(TABLE *tabPO);
int SetTableRowLabI(TABLE *tabPO,int row,char *nameS);
int SetTableColLabI(TABLE *tabPO,int col,char *nameS);
int GetTableRowLabI(TABLE *tabPO,int row,char *nameS,int max);
int GetTableColLabI(TABLE *tabPO,int col,char *nameS,int max);
int HideTableRowLabsI(TABLE *tabPO);
int HideTableColLabsI(TABLE *tabPO);
int MatchingTableRowsI(TABLE *t1PO, TABLE *t2PO, int kc);

/****************************************************************
* table_io.c
*/
int GetTableI(char *nameS, int rlab, int clab, int corn, int skp, TABLE **tPPO);
int ParseTableI(FILE *inPF, int rlab, int clab, int corn, int skp, TABLE **tPPO);
int ParseLabelRowFromBuffI(char *bufS, int line, int corn, char *cornS, WORDLIST *labsPO);
int ParseDataRowFromBuffI(char *bufS, int line, int row, WORDLIST *labsPO, NUMLIST *valsPO);
int WriteTableFileI(TABLE *tabPO, int head, int mask, char *nameS);
void DumpTable(TABLE *tabPO, int settings, int mask, FILE *outPF);
void DumpTableColHeadings(TABLE *tabPO,int mask,FILE *outPF);
void DumpTableSettings(TABLE *tabPO,int mask,FILE *outPF);
void DumpTableDescription(TABLE *tabPO, int mask, char *preS, FILE *outPF);

/****************************************************************
* table_ops.c
*/
int MixTablesI(TABLE *ftabPO, TABLE *stabPO, TABLE *ttabPO, int mask, int mix,
    int error);
int NormTableI(TABLE *tabPO, int mask, int what);
int SmoothTableI(TABLE *tabPO, int what, int win, int mask);
int QuantileTableColsI(TABLE *tabPO, int parts, int mask, int verbose);

/****************************************************************
* table_stat.c
*/
int TableColStatsI(TABLE *tabPO, int col, int mask, DOUB *loPD, DOUB *hiPD, 
    DOUB *avPD, DOUB *sdPD);
int TableRowStatsI(TABLE *tabPO, int row, int mask, DOUB *loPD, DOUB *hiPD, 
    DOUB *avPD, DOUB *sdPD);
int TableStatsI(TABLE *tabPO, int rs, int re, int cs, int ce, int mask,
    DOUB *loPD, DOUB *hiPD, DOUB *avPD, DOUB *sdPD);

/****************************************************************
* table_str.c
*/
TABLE *CreateTablePO(int row,int col);
int TableAllocAuxI(TABLE *tabPO);
int DestroyTableI(TABLE *tabPO);
int InitTableI(TABLE *tabPO, int vals);
void SetTableMasks(TABLE *tabPO, int what);
void SetTableRowMask(TABLE *tabPO, int what);
void SetTableColMask(TABLE *tabPO, int what);
int InitTableValsI(TABLE *tabPO, DOUB vD, int mask);
int InitTableMatchingValsI(TABLE *tabPO, DOUB targD, DOUB newD, int mask);
int CopyTableI(TABLE *tabPO, int mask, int tran, TABLE **newPPO);
int CopyTableAttribsI(TABLE *tabPO, TABLE *cpPO);
int CopyTableValsI(TABLE *tabPO, TABLE *cpPO, int mask, int tran);
TABLE *NewPastedTablePO(TABLE *ftabPO, TABLE *stabPO, int mask, int rowname, int kc);

/****************************************************************
* table_val.c
*/
int RowColTableIndexI(TABLE *tabPO,int row,int col, int *inPI);
int AddTableValI(TABLE *tabPO,int row,int col,DOUB valD);
int SetTableValI(TABLE *tabPO,int row,int col,DOUB valD);
int GetTableValI(TABLE *tabPO,int row,int col, DOUB *valPD);
int ModTableValI(TABLE *tabPO,int row,int col,DOUB valD,int mix);
int SetTableDiagValI(TABLE *tabPO,DOUB valD);
int TableNumlistSameDimsI(TABLE *tabPO, NUMLIST *nlPO, int row, int error);
int GetTableRowValsI(TABLE *tabPO, int row, NUMLIST *valsPO, int mask);
int GetTableColValsI(TABLE *tabPO, int col, NUMLIST *valsPO, int mask);
int SetTableRowValsI(TABLE *tabPO, int row, NUMLIST *valsPO, int n, int mask);
int SetTableColValsI(TABLE *tabPO, int col, NUMLIST *valsPO, int n, int mask);
int GetTableNumlistI(TABLE *tabPO, NUMLIST **nlPPO, int *nPI);

/****************************************************************
* table_wf.c
*/
int WeakFirstMappingI(TABLE *tabPO, TABLE *mapPO, int mask, int psize, int verb);


#endif
