/*
* table_str.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "prim.h"
#include "table.h"

#define DB_TAB      if(DB[70])
#define DB_TAB_LO   if(DB[71])

/****************************************************************************
*   Alloc table structure and call to initialize
*/
TABLE *CreateTablePO(int row,int col)
{
    TABLE *tabPO;

    DB_TAB DB_PrI(">> CreateTablePO row=%d col=%d\n",row,col);
    tabPO = (TABLE*)ALLOC(1,sizeof(TABLE));
    if(!tabPO) {
        ERR("CreateTablePO","Failed to allocate TABLE");
        return(NULL);
    }
    tabPO->ID = TABLE_ID;
    /***
    *   If row / col non-negative, then allocate space and init
    */
    if( (row>0) && (col>0) ) { 
        DB_TAB DB_PrI("+ allocating and initilizing\n");
        tabPO->nrow = row;
        tabPO->ncol = col;
        if(!TableAllocAuxI(tabPO)) {
            CHECK_TABLE(tabPO);
            ERR("CreateTablePO","Failed to allocate Aux TABLE");
            return(NULL);
        }
        InitTableI(tabPO, TRUE);
    }
    DB_TAB DB_PrI("<< CreateTablePO %p\n",tabPO);
    return(tabPO);
}
/*************************************************************************
*   Allocate auxillary space for table
*/
int TableAllocAuxI(TABLE *tabPO)
{
    int row,col;

    DB_TAB DB_PrI(">> TableAllocAuxI %p\n",tabPO);
    VALIDATE(tabPO,TABLE_ID);
    row = tabPO->nrow;
    col = tabPO->ncol;
    DB_TAB DB_PrI("+ row=%d col=%d\n",row,col);
    if( (row<1) || (col<1) ) {
        return(FALSE);
    }
    /***
    *   Only alloc if don't already here
    */
    if( !tabPO->vals ) {
        DB_TAB DB_PrI("+ Creating vals\n");
        tabPO->vals = CreateNumlistPO(IS_DOUB,NULL,row*col);
    }
    if( !tabPO->clabs ) {
        DB_TAB DB_PrI("+ Creating clabs\n");
        tabPO->clabs = CreateWordlistPO(NULL,NSIZE);
    }
    if( !tabPO->rlabs ) {
        DB_TAB DB_PrI("+ Creating rlabs\n");
        tabPO->rlabs = CreateWordlistPO(NULL,NSIZE);
    }
    if( !tabPO->cmask ) {
        DB_TAB DB_PrI("+ Creating cmask\n");
        tabPO->cmask = (char*)ALLOC(col,sizeof(char));
    }
    if( !tabPO->rmask ) {
        DB_TAB DB_PrI("+ Creating rmask\n");
        tabPO->rmask = (char*)ALLOC(row,sizeof(char));
    }
    if( (!tabPO->clabs) || (!tabPO->rlabs) || (!tabPO->cmask) || 
        (!tabPO->rmask) || (!tabPO->vals) )
    {
        printf("TableAllocAuxI attempted %d row, %d col = %d\n",row,col,row*col);
        printf(" vals=%p, clabs=%p, $rlabs=%p, cmask=%p, rmask=%p\n",
                tabPO->vals, tabPO->clabs, tabPO->rlabs, tabPO->cmask, tabPO->rmask);
        DB_TAB DB_PrI("<< TableAllocAuxI FALSE\n");
        return(FALSE);
    }
    DB_TAB DB_PrI("<< TableAllocAuxI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Free table structure
*/
int DestroyTableI(TABLE *tabPO)
{
    DB_TAB DB_PrI(">> DestroyTableI %p\n",tabPO);
    VALIDATE(tabPO,TABLE_ID);
    CHECK_WORDLIST(tabPO->clabs);
    CHECK_WORDLIST(tabPO->rlabs);
    CHECK_NUMLIST(tabPO->vals);
    CHECK_FREE(tabPO->cmask);
    CHECK_FREE(tabPO->rmask);
    CHECK_FREE(tabPO);
    DB_TAB DB_PrI("<< DestroyTableI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Initialize settings in table
*/
int InitTableI(TABLE *tabPO, int vals)
{
    DB_TAB DB_PrI(">> InitTableI %p vals=%d\n",tabPO,vals);
    VALIDATE(tabPO,TABLE_ID);
    INIT_S(tabPO->prefix);
    strcpy(tabPO->pvform2,DEF_TAB_PVFORM_S);
    strcpy(tabPO->pvform,DEF_TAB_PVFORM_S);
    strcpy(tabPO->prlform,DEF_TAB_PRLFORM_S);
    strcpy(tabPO->pvsep,DEF_TAB_PVSEP_S);
    SetTableMasks(tabPO, TRUE);
    tabPO->prefac = 1.0;
    if(vals) {
        InitTableValsI(tabPO, 0.0, FALSE);
    }
    SetTablePrecisionI(tabPO,DEF_TAB_PRENUM);
    SetNumlistNamesI(tabPO->vals, "Table_values", NULL, NSIZE);
    DB_TAB DB_PrI("<< InitTableI\n");
    return(TRUE);
}
/*************************************************************************
*   Set row / col masks to 'what' 
*/
void SetTableMasks(TABLE *tabPO, int what)
{
    InitArrayI(tabPO->rmask,IS_CHAR,0,tabPO->nrow,what);
    InitArrayI(tabPO->cmask,IS_CHAR,0,tabPO->ncol,what);
}
void SetTableRowMask(TABLE *tabPO, int what)
{
    InitArrayI(tabPO->rmask,IS_CHAR,0,tabPO->nrow,what);
}
void SetTableColMask(TABLE *tabPO, int what)
{
    InitArrayI(tabPO->cmask,IS_CHAR,0,tabPO->ncol,what);
}
/*************************************************************************
*   Initialize table values
*   If mask is TRUE, only masked ROWS/COLS are set
*/
int InitTableValsI(TABLE *tabPO, DOUB vD, int mask)
{
    int r,c;

    DB_TAB DB_PrI(">> InitTableValsI %p vD=%f mask=%d\n",tabPO,vD,mask);
    VALIDATE(tabPO,TABLE_ID);
    for(r=0;r<tabPO->nrow;r++) {
        if( (mask) && (!tabPO->rmask[r]) ) {
            continue;
        }
        for(c=0;c<tabPO->ncol;c++) {
            if( (mask) && (!tabPO->cmask[c]) ) {
                continue;
            }
            SetTableValI(tabPO,r,c,vD);
        }
    }
    DB_TAB DB_PrI("<< InitTableValsI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Reset table values ONLY if they match target value
*   If mask is TRUE, only masked ROWS/COLS are set
*/
int InitTableMatchingValsI(TABLE *tabPO, DOUB targD, DOUB newD, int mask)
{
    int r,c;
    DOUB vD;

    DB_TAB DB_PrI(">> InitTableValsI %p targD=%f new=%f mask=%d\n",tabPO,targD,newD,mask);
    VALIDATE(tabPO,TABLE_ID);
    for(r=0;r<tabPO->nrow;r++) {
        if( (mask) && (!tabPO->rmask[r]) ) {
            continue;
        }
        for(c=0;c<tabPO->ncol;c++) {
            if( (mask) && (!tabPO->cmask[c]) ) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            if(vD == targD) {
                SetTableValI(tabPO,r,c,newD);
            }
        }
    }
    DB_TAB DB_PrI("<< InitTableValsI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Returns a copied table of tabPO;
*   If mask is TRUE, only masked ROWS/COLS are copied
*   If transpose is TRUE, rows and cols are transposed
*/
int CopyTableI(TABLE *tabPO, int mask, int tran, TABLE **newPPO)
{
    int nrow,ncol;
    TABLE *cpPO;

    DB_TAB DB_PrI(">> CopyTableI %p mask=%d tran=%d\n",tabPO,mask,tran);
    VALIDATE(tabPO,TABLE_ID);
    if(*newPPO) {
        *newPPO = NULL;
    }
    /***
    *   How big will the new guy be?
    */
    nrow = GetTableRowsI(tabPO,mask);
    ncol = GetTableColsI(tabPO,mask);
    DB_TAB DB_PrI("+ table dims: row=%d col=%d\n",nrow,ncol);
    if( (nrow*ncol) < 1) {
        DB_TAB DB_PrI("<< CopyTableI no cells FALSE\n");
        return(FALSE);
    }
    /***
    *   Create new copy, maybe transposed dimenisons 
    */
    if(tran) {
        cpPO = CreateTablePO(ncol,nrow);
    }
    else {
        cpPO = CreateTablePO(nrow,ncol);
    }
    if(!cpPO) {
        printf("FAILED TO ALLOCATE NEW TABLE: %d X %d\n",nrow,ncol);
        return(FALSE);
    }
    /***
    *   Copy attributes then cell values
    */
    CopyTableAttribsI(tabPO, cpPO);
    if(!CopyTableValsI(tabPO,cpPO,mask,tran)) {
        DB_TAB DB_PrI("<< CopyTableI failed to copy vals FALSE\n");
        CHECK_TABLE(cpPO);
        return(FALSE);
    }
    *newPPO = cpPO;
    DB_TAB DB_PrI("<< CopyTableI TRUE %p new=%p\n",tabPO,cpPO);
    return(TRUE);
}
/**************************************************************************
*   Copy table attribute settings; no row-col values
*/
int CopyTableAttribsI(TABLE *tabPO, TABLE *cpPO)
{
    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(cpPO,TABLE_ID);
    DB_TAB DB_PrI(">> CopyTableAttribsI %p new=%p\n",tabPO,cpPO);
    /***
    *   Copy over
    */
    cpPO->rlab = tabPO->rlab;
    cpPO->clab = tabPO->clab;
    strcpy(cpPO->pvform,tabPO->pvform);
    strcpy(cpPO->prlform,tabPO->prlform);
    strcpy(cpPO->pvsep,tabPO->pvsep);
    strcpy(cpPO->name,tabPO->name);
    strcpy(cpPO->fname,tabPO->fname);
    SetTablePrecisionI(cpPO,tabPO->prenum);
    DB_TAB DB_PrI("<< CopyTableAttribsI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Copy row-col lables & data from table tabPO to copy table cpPO
*   If mask is TRUE, only masked ROW/COLs are copied
*   If tran is TRUE, copy is Transpose of starting
*/
int CopyTableValsI(TABLE *tabPO, TABLE *cpPO, int mask, int tran)
{
    int ok,r,c,r2,c2;
    char nameS[NSIZE];
    DOUB vD;

    DB_TAB DB_PrI(">> CopyTableValsI %p %p mask=%d\n",tabPO,cpPO,mask);
    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(cpPO,TABLE_ID);
    /***
    *   Check dims match
    */
    r = GetTableRowsI(tabPO,mask);
    c = GetTableColsI(tabPO,mask);
    r2 = GetTableRowsI(cpPO,mask);
    c2 = GetTableColsI(cpPO,mask);
    ok = TRUE;
    if(tran) {
        if( (r != c2) || (c != r2) ) {
            ok = FALSE;
        }
    }
    else {
        if( (r != r2) || (c != c2) ) {
            ok = FALSE;
        }
    }
    if(!ok) {
        printf("Problem copying values between different tables\n");
        printf("First:     %d rows, %d cols\n",r,c);
        printf("Second:    %d rows, %d cols\n",r2,c2);
        return(FALSE);
    }
    DB_TAB DB_PrI("+ r=%d c=%d r2=%d c2=%d\n",r,c,r2,c2);
    /***
    *   For each row
    */
    r2 = c2 = 0;
    for(r=0;r<tabPO->nrow;r++)
    {
        /***
        *   If masking and not set, ignore row
        */
        if( (mask) && (!tabPO->rmask[r]) ) {
            continue;
        }
        /***
        *   Row lables 
        */
        if(tabPO->rlab) {
            GetTableRowLabI(tabPO,r,nameS,NSIZE-1);
            if(tran) {
                SetTableColLabI(cpPO,r2,nameS);
            }
            else {
                SetTableRowLabI(cpPO,r2,nameS);
            }
        }
        /***
        *   For each column
        */
        c2 = 0;
        for(c=0;c<tabPO->ncol;c++)
        {
            /***
            *   If masking and not set ignore col
            */
            if( (mask) && (!tabPO->cmask[c]) ) {
                continue;
            }
            /***
            *   Col lables, but only needed first pass on row loop
            */
            if( (tabPO->clab) && (r2==0) ) {
                GetTableColLabI(tabPO,c,nameS,NSIZE-1);
                if(tran) {
                    SetTableRowLabI(cpPO,c2,nameS);
                }
                else {
                    SetTableColLabI(cpPO,c2,nameS);
                }
            }
            /***
            *   Actual values to copy, direct or transpose
            */
            GetTableValI(tabPO,r,c,&vD);
            if(tran) {
                SetTableValI(cpPO,c2,r2,vD);
            }
            else {
                SetTableValI(cpPO,r2,c2,vD);
            }
            c2++;
        }
        r2++;
    }
    DB_TAB DB_PrI("<< CopyTableValsI %p %p TRUE\n",tabPO,cpPO);
    return(TRUE);
}
/**************************************************************************
*   Paste two tables side-by-side and return new one based on first
*   If mask is TRUE, only masked ROW/COLs are copied
*   If rowname > 0, check that row names match; If kc, keep case for comparison
*/
TABLE *NewPastedTablePO(TABLE *ftabPO, TABLE *stabPO, int mask, int rowname, int kc)
{
    int c,r,nr,nc1,nc2;
    char nameS[NSIZE], snameS[NSIZE];
    TABLE *newPO;
    NUMLIST *tcolPO;
    
    VALIDATE(ftabPO,TABLE_ID);
    VALIDATE(stabPO,TABLE_ID);
    if(!SameDimTableRowsI(ftabPO,stabPO,mask,TRUE)) {
        return(NULL);
    }
    /***
    *   New table and temp column array
    */
    nr = GetTableRowsI(ftabPO,mask);
    nc1 = GetTableColsI(ftabPO,mask);
    nc2 = GetTableColsI(stabPO,mask);
    newPO = CreateTablePO(nr, nc1 + nc2);
    tcolPO = CreateNumlistPO(IS_DOUB, NULL, nr);
    if( (!newPO) || (!tcolPO) ) {
        CHECK_TABLE(newPO);     
        CHECK_NUMLIST(tcolPO);
        return(NULL);
    }
    CopyTableAttribsI(ftabPO, newPO);
    /***
    *   Copy columns into new; First table 1 then table 2 
    */
    for(c=0; c<nc1; c++) {
        GetTableColValsI(ftabPO,c,tcolPO,mask);
        SetTableColValsI(newPO,c,tcolPO,nr,mask);
        GetTableColLabI(ftabPO,c,nameS,NSIZE-1);
        SetTableColLabI(newPO,c,nameS);
    }
    for(c=0; c<nc2; c++) {
        GetTableColValsI(stabPO,c,tcolPO,mask);
        SetTableColValsI(newPO,c+nc1,tcolPO,nr,mask);
        GetTableColLabI(stabPO,c,nameS,NSIZE-1);
        SetTableColLabI(newPO,c+nc1,nameS);
    }
    CHECK_NUMLIST(tcolPO);
    /* Row labels */
    for(r=0; r<nr; r++) {
        GetTableRowLabI(ftabPO,r,nameS,NSIZE-1);
        if(rowname > 0) {
            GetTableRowLabI(stabPO,r,snameS,NSIZE-1);
            if(!SameStringsI(nameS,snameS,kc)) {
                printf("# WARNING: table paste row[%d] name mismatch: |%s| v |%s|\n",
                    r,nameS,snameS);
            }
        }
        SetTableRowLabI(newPO,r,nameS);
    }
    return(newPO);
}
