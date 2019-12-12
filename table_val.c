/*
* table_val.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
*
* The programs and source code of the vertools collection are free software.
* They are distributed in the hope that they will be useful,
* WITHOUT ANY WARRANTY OF FITNESS FOR ANY PARTICULAR PURPOSE.  
* 
* Permission is granted for research, educational, and possibly commercial use 
*   and modification as long as 1) Code and any derived works are not 
*   redistributed for any fee, and 2) Proper credit is given to the authors. 
*   If you wish to include this software in a product, or use it commercially,
*   please contact the authors.
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
*   Index for cell given table + row + col; 
*   Table has to have the number of columns set
*/
int RowColTableIndexI(TABLE *tabPO,int row,int col, int *inPI)
{
    int i;

    VALIDATE(tabPO,TABLE_ID);
    if( (row>=0) && (col>=0) && (col<tabPO->ncol) ) {
        i = row * tabPO->ncol + col;
        if(inPI) {
            *inPI = i;
        }
        return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************
*   Set table value, adding space if needed
*/
int AddTableValI(TABLE *tabPO,int row,int col,DOUB valD)
{
    int i;

    VALIDATE(tabPO,TABLE_ID);
    if(RowColTableIndexI(tabPO,row,col,&i)) {
        SetTableValPrecisionI(tabPO,valD,&valD);
        if(AddNumlistDoubI(tabPO->vals,i,valD)) {
            /***
            *   Update number of rows and cols that we now have
            */
            tabPO->nrow = MAX_NUM(tabPO->nrow, row);
            tabPO->ncol = MAX_NUM(tabPO->ncol, col);
            return(TRUE);
        }
    }
    return(FALSE);
}
/**************************************************************************
*   Set value for cell [row][col] in table 
*   First limits to specificed precision
*/
int SetTableValI(TABLE *tabPO,int row,int col,DOUB valD)
{
    int i;

    VALIDATE(tabPO,TABLE_ID);
    if( (row>=0) && (row<tabPO->nrow) && (col>=0) && (col<tabPO->ncol) )
    {
        if(RowColTableIndexI(tabPO,row,col,&i)) {
            SetTableValPrecisionI(tabPO,valD,&valD);
            if(SetNumlistDoubI(tabPO->vals,i,valD)) {
                return(TRUE);
            }
        }
    }
    return(FALSE);
}
/**************************************************************************
*   Get value for cell [row][col] in table 
*/
int GetTableValI(TABLE *tabPO,int row,int col, DOUB *valPD)
{
    int i;

    VALIDATE(tabPO,TABLE_ID);
    if( (row>=0) && (row<tabPO->nrow) && (col>=0) && (col<tabPO->ncol) )
    {
        if(RowColTableIndexI(tabPO,row,col,&i)) {
            if(GetNumlistDoubI(tabPO->vals,i,valPD)) {
                return(TRUE);
            }
        }
    }
    return(FALSE);
}
/**************************************************************************
*   Modify table value with passed arg and operation code
*/
int ModTableValI(TABLE *tabPO,int row,int col,DOUB valD,int mix)
{
    DOUB fD;

    VALIDATE(tabPO,TABLE_ID);
    if( (row>=0) && (row<tabPO->nrow) && (col>=0) && (col<tabPO->ncol) )
    {
        GetTableValI(tabPO,row,col,&fD);
        switch(mix)
        {
            case MATH_ADD:      fD += valD; break;
            case MATH_SUB:      fD -= valD; break;
            case MATH_MUL:      fD *= valD; break;
            case MATH_DIV:      
                if(valD == 0.0)
                {   fD = 0.0;   }
                else
                {   fD = fD / valD; }
                break;
            case MATH_MIN:      fD = MIN_NUM(fD,valD);  break;
            case MATH_MAX:      fD = MAX_NUM(fD,valD);  break;
            default:
                printf("Bogus mix code: %d\n",mix);
                ERR("ModTableValI","bad mixing operation code");
                return(FALSE);
        }
        SetTableValI(tabPO,row,col,fD);
        return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************
*   Set value for diagonal in table 
*/
int SetTableDiagValI(TABLE *tabPO,DOUB valD)
{
    int r;

    VALIDATE(tabPO,TABLE_ID);
    for(r=0;r<tabPO->nrow;r++) {
        if(r<tabPO->ncol) {
            SetTableValI(tabPO,r,r,valD);
        }
    }
    return(r);
}
/*************************************************************************
*   Check that table (row or col) dimensions are same as numlist
*   If error and dimensions don't agree, report values
*/
int TableNumlistSameDimsI(TABLE *tabPO, NUMLIST *nlPO, int row, int error)
{
    int t,n,ok;

    n = GetNumlistLengthI(nlPO);
    if( row ) {
        t = GetTableRowsI(tabPO,FALSE);
    }
    else {
        t = GetTableColsI(tabPO,FALSE);
    }
    ok = (t == n);
    if( (!ok) && error ) {
        if(row) {
            printf("Table row count: %d\n",t);
        }
        else {
            printf("Table col count: %d\n",t);
        }
        printf("Number count:    %d\n",n);
    }
    return(ok);
}
/**************************************************************************
*   Copy vals from table row into value array. 
*   If mask is TRUE, only masked cols are copied
*   Return number of values copied
*/
int GetTableRowValsI(TABLE *tabPO, int row, NUMLIST *valsPO, int mask)
{
    int c,n;
    DOUB vD;

    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(valsPO,NUMLIST_ID);
    if( (row<0) || (row>=tabPO->nrow) ) {
        return(FALSE);
    }
    n = 0;
    for(c=0;c<tabPO->ncol;c++) {
        if( mask && (!tabPO->cmask[c]) ) {
            continue;
        }
        GetTableValI(tabPO,row,c,&vD);
        /***
        *   Add rather than Set; Add increases space as needed
        */
        AddNumlistDoubI(valsPO,n++,vD);
    }
    SetNumlistLengthI(valsPO,n);
    return(n);
}
/**************************************************************************
*   Copy vals from table col into value array. 
*   If mask is TRUE, only masked rows are copied
*   Return number of values copied
*/
int GetTableColValsI(TABLE *tabPO, int col, NUMLIST *valsPO, int mask)
{
    int r,n;
    DOUB vD;

    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(valsPO,NUMLIST_ID);
    if( (col<0) || (col>=tabPO->ncol) ) {
        return(FALSE);
    }
    n = 0;
    for(r=0; r<tabPO->nrow; r++) {
        if( mask && (!tabPO->rmask[r]) ) {
            continue;
        }
        GetTableValI(tabPO,r,col,&vD);
        /***
        *   Add rather than Set; Add increases space as needed
        */
        AddNumlistDoubI(valsPO,n++,vD);
    }
    SetNumlistLengthI(valsPO,n);
    return(n);
}
/**************************************************************************
*   Copy vals from passed arry to table row. 
*   Return number of values copied
*/
int SetTableRowValsI(TABLE *tabPO, int row, NUMLIST *valsPO, int n, int mask)
{
    int c,s;
    DOUB vD;

    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(valsPO,NUMLIST_ID);
    if( (row<0) || (row>=tabPO->nrow) ) {
        return(FALSE);
    }
    s = 0;
    n = MIN_NUM(n,tabPO->ncol);
    for(c=0;c<n;c++) {
        if( mask && (!tabPO->cmask[c]) ) {
            continue;
        }
        GetNumlistDoubI(valsPO,s++,&vD);
        SetTableValI(tabPO,row,c,vD);
    }
    return(s);
}
/**************************************************************************
*   Copy vals from passed array to table col 
*   Return number of values copied
*/
int SetTableColValsI(TABLE *tabPO, int col, NUMLIST *valsPO, int n, int mask)
{
    int r,s;
    DOUB vD;

    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(valsPO,NUMLIST_ID);
    if( (col<0) || (col>=tabPO->ncol) ) {
        return(FALSE);
    }
    s = 0;
    n = MIN_NUM(n,tabPO->nrow);
    for(r=0;r<n;r++) {
        if( mask && (!tabPO->rmask[r]) ) {
            continue;
        }
        GetNumlistDoubI(valsPO,s++,&vD);
        SetTableValI(tabPO,r,col,vD);
    }
    return(s);
}
/*************************************************************************/
int GetTableNumlistI(TABLE *tabPO, NUMLIST **nlPPO, int *nPI)
{
    VALIDATE(tabPO,TABLE_ID);
    *nlPPO = tabPO->vals;
    if(nPI) {
        *nPI = GetNumlistLengthI(tabPO->vals);
    }
    return(TRUE);
}
