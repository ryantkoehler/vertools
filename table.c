/*
* table.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "prim.h"
#include "table.h"

#define DB_TAB      if(DB[70])
#define DB_TAB_LO   if(DB[71])

/**************************************************************************
*   Set value precision limit for table
*/
int SetTablePrecisionI(TABLE *tabPO, int prenum)
{
    int i,r,c;
    DOUB mD,vD;

    DB_TAB DB_PrI(">> SetTablePrecisionI %p prenum=%d\n",tabPO,prenum);
    VALIDATE(tabPO,TABLE_ID);
    if( (prenum>MAX_TAB_PRENUM) || (prenum<MIN_TAB_PRENUM) ) {
        printf("Table precision limited %d to %d digits (after decimal)\n",
            MIN_TAB_PRENUM, MAX_TAB_PRENUM);
        printf("Value %d won't work; ignored\n",prenum);
        return(FALSE);
    }
    tabPO->prenum = prenum;
    tabPO->prefac = 1.0;
    /***
    *   Now set working values 
    */
    if(prenum<0) {
        mD = 0.1;   
    }
    else {
        mD = 10.0;  
    }
    /***
    *   Now set working values by incrementing n steps
    */
    i = 0;
    while(i<ABS_VAL(tabPO->prenum)) {
        tabPO->prefac *= mD;
        i++;
    }
    DB_TAB DB_PrI("+ prefac=%1.12f\n",tabPO->prefac);
    /***
    *   Update all values with new settings
    */
    for(r=0;r<tabPO->nrow;r++)
    {
        for(c=0;c<tabPO->ncol;c++)
        {
            GetTableValI(tabPO,r,c,&vD);
            SetTableValPrecisionI(tabPO, vD, &vD);
            SetTableValI(tabPO,r,c,vD);
        }
    }
    DB_TAB DB_PrI("<< SetTablePrecisionI %p TRUE\n",tabPO);
    return(TRUE);   
}
/***************************************************************************
*   Limits precision of val based on table settings and sets via pointer 
*   new value
*/
int SetTableValPrecisionI(TABLE *tabPO, DOUB valD, DOUB *newPD)
{
    DOUB intD,fracD;

    DB_TAB_LO DB_PrI(">> SetTableValPrecisionI val=%1.12f\n",valD);
    VALIDATE(tabPO,TABLE_ID);
    DB_TAB_LO DB_PrI("+ prenum=%d prefac=%1.12f\n",tabPO->prenum,tabPO->prefac);
    /***
    *   Only mess with the fractional part
    */
    fracD = modf(valD,&intD);
    fracD = DNUM( ROUND(fracD * tabPO->prefac) ) / tabPO->prefac;
    *newPD = intD + fracD;
    DB_TAB_LO DB_PrI("<< SetTableValPrecisionI new=%1.12f TRUE\n",*newPD);
    return(TRUE);
}
/**************************************************************************
*   Same dimension tables?
*/
int SameDimTablesI(TABLE *tabPO, TABLE *stabPO, int mask,int error)
{
    int ok;

    ok = ( SameDimTableRowsI(tabPO, stabPO, mask, error) &&
           SameDimTableColsI(tabPO, stabPO, mask, error) );
    return(ok);
}
/**************************************************************************/
int SameDimTableRowsI(TABLE *tabPO, TABLE *stabPO, int mask,int error)
{
    int fr,sr;

    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(stabPO,TABLE_ID);
    fr = GetTableRowsI(tabPO,mask);
    sr = GetTableRowsI(stabPO,mask);
    if( fr != sr ) {
        if(error) {
            printf("Table ROW COUNTS DIFFER!!!\n");
            printf("   First  = %s = %d\n",tabPO->name,fr);
            printf("   Second = %s = %d\n",stabPO->name,sr);
        }
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************/
int SameDimTableColsI(TABLE *tabPO, TABLE *stabPO, int mask,int error)
{
    int fc,sc;
    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(stabPO,TABLE_ID);
    fc = GetTableColsI(tabPO,mask);
    sc = GetTableColsI(stabPO,mask);
    if( fc != sc ) {
        if(error) {
            printf("Table COLUMN COUNTS DIFFER!!!\n");
            printf("   First  = %s = %d\n",tabPO->name,fc);
            printf("   Second = %s = %d\n",stabPO->name,sc);
        }
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************
*   Set masking for row in table 
*/
int SetTableRowMaskI(TABLE *tabPO,int row,int what)
{
    VALIDATE(tabPO,TABLE_ID);
    if( (row>=0) && (row<tabPO->nrow) ) {
        if(what)
        {   tabPO->rmask[row] = TRUE; }
        else
        {   tabPO->rmask[row] = FALSE; }
        return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************
*   Set masking for col in table 
*/
int SetTableColMaskI(TABLE *tabPO,int col,int what)
{
    VALIDATE(tabPO,TABLE_ID);
    if( (col>=0) && (col<tabPO->ncol) )
    {
        if(what)
        {   tabPO->cmask[col] = TRUE; }
        else
        {   tabPO->cmask[col] = FALSE; }
        return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************
*   Set masking for range of rows, start to end (exclusive)
*/
int SetTableRowRangeMaskI(TABLE *tabPO,int srow, int erow, int what)
{
    int i,any;
    
    any = 0;
    for(i=srow; i<erow; i++)
    {
        if(SetTableRowMaskI(tabPO,i,what)) {
            any++;
        }
    }
    return(any);
}
/**************************************************************************
*   Set masking for range of cols, start to end (exclusive)
*/
int SetTableColRangeMaskI(TABLE *tabPO,int scol, int ecol, int what)
{
    int i,any;
    
    any = 0;
    for(i=scol; i<ecol; i++)
    {
        if(SetTableColMaskI(tabPO,i,what)) {
            any++;
        }
    }
    return(any);
}
/**************************************************************************
*   Get masking status for row in table 
*/
int GetTableRowMaskI(TABLE *tabPO,int row)
{
    VALIDATE(tabPO,TABLE_ID);
    if( (row>=0) && (row<tabPO->nrow) ) {
        return(tabPO->rmask[row]);
    }
    return(BOGUS);
}
/**************************************************************************
*   Get masking status for col in table 
*/
int GetTableColMaskI(TABLE *tabPO,int col)
{
    VALIDATE(tabPO,TABLE_ID);
    if( (col>=0) && (col<tabPO->ncol) ) {
        return(tabPO->cmask[col]);
    }
    return(BOGUS);
}
/*************************************************************************/
int GetTableDimsI(TABLE *tabPO, int *rowPI, int *colPI, int mask)
{
    VALIDATE(tabPO,TABLE_ID);
    if(rowPI) {
        *rowPI = GetTableRowsI(tabPO,mask);
    }
    if(colPI) {
        *colPI = GetTableColsI(tabPO,mask);
    }
    return(TRUE);
}
/**************************************************************************
*   How many rows in table? If mask is true, only masked rows are counted
*/
int GetTableRowsI(TABLE *tabPO,int mask)
{
    int r;

    VALIDATE(tabPO,TABLE_ID);
    if(mask) {
        r = MaskCountI(tabPO->rmask,tabPO->nrow);
    }
    else {
        r = tabPO->nrow;
    }
    return(r);
}
/**************************************************************************
*   How many cols in table? If mask is true, only masked cols are counted
*/
int GetTableColsI(TABLE *tabPO,int mask)
{
    int c;

    VALIDATE(tabPO,TABLE_ID);
    if(mask) {
        c = MaskCountI(tabPO->cmask,tabPO->ncol);
    }
    else {
        c = tabPO->ncol;
    }
    return(c);
}
/**************************************************************************
*   Set the name of the table 
*/
int SetTableNamesI(TABLE *tabPO, char *nameS, char *fnameS, int max)
{
    int n;

    VALIDATE(tabPO,TABLE_ID);
    if(nameS){
        n = (max < 0) ? strlen(nameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(tabPO->name,nameS,n);
        tabPO->name[n] = '\0';
    }
    if(fnameS) {
        n = (max < 0) ? strlen(fnameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(tabPO->fname,fnameS,n);
        tabPO->fname[n] = '\0';
    }
    return(TRUE);
}
/**************************************************************************
*   Copy name of table into nameS up to max
*/
int GetTableNamesI(TABLE *tabPO, char *nameS, char *fnameS, int max)
{
    int n;

    VALIDATE(tabPO,TABLE_ID);
    if(nameS) {
        n = (max < 0) ? strlen(nameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(nameS,tabPO->name,n);
        nameS[n] = '\0';
    }
    if(fnameS) {
        n = (max < 0) ? strlen(fnameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(fnameS,tabPO->fname,n);
        fnameS[n] = '\0';
    }
    return(TRUE);
}
/**************************************************************************
*   Find named row
*/
int FindTableNamedRowI(TABLE *tabPO, int off, char *nameS, int kc, int nc)
{
    int i,r;
    char rowS[NSIZE];

    r = BOGUS;
    /***
    *   No row labels then no chance
    */
    if(!tabPO->rlab) {
        return(r);
    }
    /***
    *   Each row, after offset
    */
    LIMIT_NUM(off,0,tabPO->nrow);
    for(i=off; i<tabPO->nrow; i++) {
        GetTableRowLabI(tabPO,i,rowS,-1);
        if(kc) {
            if( ((nc>0) && (!strncmp(rowS,nameS,nc))) || (!strcmp(rowS,nameS)) ) {
                r = i;
                break;
            }
        }
        else {
            if( ((nc>0) && (!strcasecmp(rowS,nameS))) || (!strcasecmp(rowS,nameS)) ) {
                r = i;
                break;
            }
        }
    }
    return(r);
}
/**************************************************************************
*   Find named col
*/
int FindTableNamedColI(TABLE *tabPO, int off, char *nameS, int kc, int nc)
{
    int i,c;
    char colS[NSIZE];

    c = BOGUS;
    /***
    *   No col labels then no chance
    */
    if(!tabPO->clab) {
        return(c);
    }
    /***
    *   Each row, after offset
    */
    LIMIT_NUM(off,0,tabPO->ncol);
    for(i=off; i<tabPO->ncol; i++) {
        GetTableColLabI(tabPO,i,colS,NSIZE-1);
        if(kc) {
            if( ((nc>0) && (!strncmp(colS,nameS,nc))) || (!strcmp(colS,nameS)) ) {
                c = i;
                break;
            }
        }
        else {
            if( ((nc>0) && (!strcasecmp(colS,nameS))) || (!strcasecmp(colS,nameS)) ) {
                c = i;
                break;
            }
        }
    }
    return(c);
}
/*************************************************************************/
int AutoTableOutFormattingI(TABLE *tabPO, int val, int row)
{
    char valS[DEF_BS], val2S[DEF_BS], rowS[DEF_BS];

    VALIDATE(tabPO,TABLE_ID);
    INIT_S(valS);
    INIT_S(val2S);
    INIT_S(rowS);
    if(val) {
/*
xxx
*/
        NumlistAutoFormatStringI(tabPO->vals, NULL, NULL, valS);
    }
    if(row) {
        WordlistAutoFormatStringI(tabPO->rlabs, NULL, rowS);
    }
    SetTablePrintformI(tabPO, valS, val2S, NULL, rowS, NULL);
    return(TRUE); 
}
/**************************************************************************
*   Set the print format strings for when table is dumped
*       valS = Format string to dump (number) values        (e.g. "%5.2f")
*       sepS = Separator between items                      (e.g. "\t")
*       rowS = Format string to dump (string) row lables    (e.g. "%-12s")
*       preS = Prefix string for lines                      (e.g. "# ")
*/
int SetTablePrintformI(TABLE *tabPO, char *valS, char *val2S, char *sepS, 
    char *rowS, char *preS)
{
    VALIDATE(tabPO,TABLE_ID);
    if(valS && (!NO_S(valS))) {
        strcpy(tabPO->pvform,valS);
    }
    if(val2S && (!NO_S(val2S))) {
        strcpy(tabPO->pvform2,val2S);
    }
    if(sepS && (!NO_S(sepS))) {
        strcpy(tabPO->pvsep,sepS);
    }
    if(rowS && (!NO_S(rowS))) {
        strcpy(tabPO->prlform,rowS);
    }
    if(preS && (!NO_S(preS))) {
        strcpy(tabPO->prefix,preS);
    }
    return(TRUE);
}
/*************************************************************************/
int GetTablePrintformI(TABLE *tabPO, char *valS, char *val2S, char *sepS, 
    char *rowS, char *preS)
{
    VALIDATE(tabPO,TABLE_ID);
    if(valS) {
        strcpy(valS,tabPO->pvform);
    }
    if(val2S) {
        strcpy(val2S,tabPO->pvform2);
    }
    if(sepS) {
        strcpy(sepS,tabPO->pvsep);
    }
    if(rowS) {
        strcpy(rowS,tabPO->prlform);
    }
    if(preS) {
        strcpy(preS,tabPO->prefix);
    }
    return(TRUE);
}
/*************************************************************************/
int IsTableCSVFormatI(TABLE *tabPO)
{
    VALIDATE(tabPO,TABLE_ID);
    return(tabPO->pvsep[0] == ',');
}
/**************************************************************************
*   Set row label to nameS, Also sets that table has row lables 
*/
int SetTableRowLabI(TABLE *tabPO,int row,char *nameS)
{
    int ok;

    VALIDATE(tabPO,TABLE_ID);
    ok = FALSE;
    if( (row>=0) && (row<tabPO->nrow) )
    {
        if(AddWordlistWordI(tabPO->rlabs, row, nameS, TRUE)) {
            ok++;
            tabPO->rlab = TRUE;
        }
    }
    return(ok);
}
/**************************************************************************
*   Set col label to nameS, Also sets that table has col labels
*/
int SetTableColLabI(TABLE *tabPO,int col,char *nameS)
{
    int ok;

    VALIDATE(tabPO,TABLE_ID);
    ok = FALSE;
    if( (col>=0) && (col<tabPO->ncol) )
    {
        if(AddWordlistWordI(tabPO->clabs, col, nameS, TRUE)) {
            ok++;
            tabPO->clab = TRUE;
        }
    }
    return(ok);
}
/**************************************************************************
*   Get row label from table and copy it into nameS up to max
*/
int GetTableRowLabI(TABLE *tabPO,int row,char *nameS,int max)
{
    VALIDATE(tabPO,TABLE_ID);
    INIT_S(nameS);
    return(GetWordlistWordI(tabPO->rlabs, row, nameS, max));
}
/**************************************************************************
*   Get col label from table and copy it into nameS up to max
*/
int GetTableColLabI(TABLE *tabPO,int col,char *nameS,int max)
{
    VALIDATE(tabPO,TABLE_ID);
    INIT_S(nameS);
    return(GetWordlistWordI(tabPO->clabs, col, nameS, max));
}
/**************************************************************************
*   Turn off flag for row lables
*/
int HideTableRowLabsI(TABLE *tabPO)
{
    VALIDATE(tabPO,TABLE_ID);
    tabPO->rlab = FALSE;
    return(TRUE);
}
/**************************************************************************
*   Turn off flag for col lables
*/
int HideTableColLabsI(TABLE *tabPO)
{
    VALIDATE(tabPO,TABLE_ID);
    tabPO->clab = FALSE;
    return(TRUE);
}
/****************************************************************************
*   Compare row number and names in two tables
*/
int MatchingTableRowsI(TABLE *t1PO, TABLE *t2PO, int kc)
{
    int i;
    char n1S[NSIZE],n2S[NSIZE];

    VALIDATE(t1PO,TABLE_ID);
    VALIDATE(t2PO,TABLE_ID);
    if(t1PO->nrow != t2PO->nrow) {
        return(FALSE);
    }
    /***
    *   Either one has no row labels to compare
    */
    if( (!t1PO->rlab) || (!t2PO->rlab) ) {
        return(TRUE);
    }
    /***
    *   Compare row lables
    */
    for(i=0;i<t1PO->nrow;i++)
    {
        GetTableRowLabI(t1PO,i,n1S,-1);
        GetTableRowLabI(t2PO,i,n2S,-1);
        if(!SameStringsI(n1S,n2S,kc)) {
            return(FALSE);
        }
    }
    return(TRUE);
}
