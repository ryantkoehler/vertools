/*
* table_ops.c
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prim.h"
#include "table.h"
#include "score.h"

#define DB_TAB      if(DB[70])
#define DB_TAB_LO   if(DB[71])

/**************************************************************************
*   Mix values in first two tables into third
*   All dimensions must be compatible
*/
int MixTablesI(TABLE *ftabPO, TABLE *stabPO, TABLE *ttabPO, int mask, int mix,
    int error)
{
    int r,c;
    DOUB fD,sD,vD;

    VALIDATE(ftabPO,TABLE_ID);
    VALIDATE(stabPO,TABLE_ID);
    VALIDATE(ttabPO,TABLE_ID);
    if(!SameDimTablesI(ftabPO, stabPO, mask, error)) {
        return(FALSE);
    }
    if(!SameDimTablesI(ftabPO, ttabPO, mask, error)) {
        return(FALSE);
    }
    /***
    *   Get and set new value for each cell
    */
    for(c=0;c<ftabPO->ncol;c++)
    {
        if( (mask) && (!ftabPO->cmask[c]) ) {
            continue;
        }
        for(r=0;r<ftabPO->nrow;r++)
        {
            if( (mask) && (!ftabPO->rmask[r]) ) {
                continue;
            }
            GetTableValI(ftabPO,r,c,&fD);
            GetTableValI(stabPO,r,c,&sD);
            switch(mix)
            {
                case MATH_ADD:      vD = fD + sD;   break;
                case MATH_SUB:      vD = fD - sD;   break;
                case MATH_MUL:      vD = fD * sD;   break;
                case MATH_DIV:      
                    if(sD == 0.0)
                    {   vD = 0.0;   }
                    else
                    {   vD = fD / sD;   }
                    break;
                case MATH_MIN:      vD = MIN_NUM(fD,sD);    break;
                case MATH_MAX:      vD = MAX_NUM(fD,sD);    break;
                default:
                    printf("Bogus mix code: %d\n",mix);
                    ERR("MixTablesI","bad mixing operation code");
                    return(FALSE);
            }
            SetTableValI(ttabPO,r,c,vD);
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Handle normalization of table
*/
int NormTableI(TABLE *tabPO, int mask, int what)
{
    int nr,nc,r,c;
    DOUB vD,avD,sdD;

    VALIDATE(tabPO,TABLE_ID);
    nr = GetTableRowsI(tabPO, mask);
    nc = GetTableColsI(tabPO, mask);
    switch(what)
    {
        case TABLE_COL:
            for(c=0;c<nc;c++)
            {
                if( (mask) && (!tabPO->cmask[c]) ) {
                    continue;
                }
                if(!TableColStatsI(tabPO,c,mask,NULL,NULL,&avD,&sdD)) {
                    avD = 0.0; 
                }
                if(sdD == 0.0) {
                    sdD = 1.0;
                }
                for(r=0;r<nr;r++)
                {
                    if( (mask) && (!tabPO->rmask[r]) ) {
                        continue;
                    }
                    GetTableValI(tabPO,r,c,&vD);
                    vD = (vD - avD) / sdD;
                    SetTableValI(tabPO,r,c,vD);
                }
            }
            break;
        case TABLE_ROW:
            for(r=0;r<nr;r++)
            {
                if( (mask) && (!tabPO->rmask[r]) ) {
                     continue;
                }
                if(!TableRowStatsI(tabPO,r,mask,NULL,NULL,&avD,&sdD)) {
                    avD = 0.0; 
                }
                if(sdD == 0.0) {
                    sdD = 1.0;
                }
                for(c=0;c<nc;c++) {
                    if( (mask) && (!tabPO->cmask[c]) ) {
                        continue;
                    }
                    GetTableValI(tabPO,r,c,&vD);
                    vD = (vD - avD) / sdD;
                    SetTableValI(tabPO,r,c,vD);
                }
            }
            break;
        case TABLE_FULL:
            if(!TableStatsI(tabPO,-1,-1,-1,-1,mask,NULL,NULL,&avD,&sdD)) {
                avD = 0.0; 
            }
            if(sdD == 0.0) {
                sdD = 1.0;
            }
            for(r=0;r<nr;r++) {
                if( (mask) && (!tabPO->rmask[r]) ) {
                    continue;
                }
                for(c=0;c<nc;c++) {
                    if( (mask) && (!tabPO->cmask[c]) ) {
                        continue;
                    }
                    GetTableValI(tabPO,r,c,&vD);
                    vD = (vD - avD) / sdD;
                    SetTableValI(tabPO,r,c,vD);
                }
            }
            break;
        case TABLE_DIAG:
            /* Norm rows by diagonal value */
            for(r=0;r<nr;r++) {
                if( (mask) && (!tabPO->rmask[r]) ) {
                    continue;
                }
                if(r >= nc) {
                    break;
                }
                /* Diag; If zero, ignore this row (sham? do something else?) */
                GetTableValI(tabPO,r,r,&avD);
                if(avD == 0.0) {
                    continue;
                }
                for(c=0;c<nc;c++) {
                    if( (mask) && (!tabPO->cmask[c]) ) {
                        continue;
                    }
                    GetTableValI(tabPO,r,c,&vD);
                    vD = vD / avD;
                    SetTableValI(tabPO,r,c,vD);
                }
            }
            break;
        default:
            printf("Bogus table code: %d\n",what);
            ERR("NormTableI","Bad table code");
    }
    return(TRUE);
}
/**************************************************************************
*   Handle normalization of table
*/
int SmoothTableI(TABLE *tabPO, int what, int win, int mask)
{
    int i, n, nr, nc;
    NUMLIST *n1PO, *n2PO;

    VALIDATE(tabPO,TABLE_ID);
    nr = GetTableRowsI(tabPO,FALSE);
    nc = GetTableColsI(tabPO,FALSE);
    /***
    *   Get num lists to grab rows / cols out of table for smoothing
    */
    n1PO = CreateNumlistPO(IS_DOUB, NULL, 0);
    n2PO = CreateNumlistPO(IS_DOUB, NULL, 0);
    if( (!n1PO) || (!n2PO) ) {
        printf("Failed to allocate temp arrays\n");
        ERR("SmoothTableI","Alloc failed");
        CHECK_NUMLIST(n1PO);
        CHECK_NUMLIST(n2PO);
        return(FALSE);
    }
    /***
    *   Each row or col   
    */
    switch(what)
    {
        case TABLE_COL:
            for(i=0;i<nc;i++) {
                n = GetTableColValsI(tabPO,i,n1PO,mask);
                SmoothNumlistI(n1PO, n2PO, win);
                SetTableColValsI(tabPO,i,n2PO,n,mask);
            }
            break;
        case TABLE_ROW:
            for(i=0;i<nr;i++) {
                n = GetTableRowValsI(tabPO,i,n1PO,mask);
                SmoothNumlistI(n1PO, n2PO, win);
                SetTableRowValsI(tabPO,i,n2PO,n,mask);
            }
            break;
    }
    CHECK_NUMLIST(n1PO);
    CHECK_NUMLIST(n2PO);
    return(TRUE);
}
/**************************************************************************
*   Quantil-ize table columns
*/
int QuantileTableColsI(TABLE *tabPO, int parts, int mask, int verbose)
{
    int c,r,nr,nc,quant,qi,lastb;
    DOUB vD,qvalD,bvalD,jumpD;
    SCOREC *cvalsPO;
    char cnameS[NSIZE],vformS[DEF_BS];

    VALIDATE(tabPO,TABLE_ID);
    nr = GetTableRowsI(tabPO,mask);
    nc = GetTableColsI(tabPO,mask);
    GetTablePrintformI(tabPO,vformS,NULL,NULL,NULL,NULL);
    /***
    *   Cannot split into more parts than rows
    *   Alloc sorting structure for n-rows
    */
    LIMIT_NUM(parts,0,nr);
    if(parts < 2) {
        return(TRUE);
    }
    if(!(cvalsPO = CreateScorecsPO(nr))) {
        PROBLINE;
        printf("Failed to allocate temp struct for %d rows (QuantileTableColsI)\n",nr);        
        return(FALSE);
    }
    /***
    *   Jump for each quantile; Tell story for whole set; 
    */
    jumpD = DNUM(nr) / DNUM(parts);
    if(verbose) {
        printf("# QuanTab splitting %d rows into %d parts (%0.3f)\n",nr,parts,jumpD);
    }
    lastb = INT(DNUM(nr) * DNUM(parts-1) / DNUM(parts)) + 1;
    /***
    *   Each column
    */
    for(c=0;c<nc;c++)
    {
        if( (mask) && (!tabPO->cmask[c]) ) {
            continue;
        }
        GetTableColLabI(tabPO,c,cnameS,-1);
        /*** 
        *   Copy values into structure then sort
        */
        for(r=0;r<nr;r++)
        {
            GetTableValI(tabPO,r,c,&vD);
            cvalsPO[r].sc = vD;
        }
        SortScorecVals(cvalsPO,nr,1);
        /***
        *   Initialize; "bins" are at least jump wide, but can be wider if ties
        */
        quant = 1;
        bvalD = jumpD;
        qi = INT(bvalD) - 1;
        BOG_CHECK(qi >= nr);
        qvalD = cvalsPO[qi].sc;
/*
printf("first bval=%f qi=%d qval=%f\n",bvalD,qi,qvalD);
*/
        if(verbose) {
            printf("# QuanTab col_%02d %s quant %d ind %d val <= ",c,cnameS,quant,qi);
            printf(vformS,qvalD);
            printf("\n");
        }
        /***
        *   Now replace records with quantiles 
        *   ... Not clear correct partitioning (sham?); No "little" last bin
        *   ... Calc (logic?) for moving to next bin also a mess (sham?)
        */
        for(r=0;r<nr;r++)
        {
            /*if( (cvalsPO[r].sc > qvalD) ){ */
            if( (cvalsPO[r].sc > qvalD) && (r<=lastb) ){ 
                quant++;
                bvalD = (DNUM(r) > bvalD) ? DNUM(r) : bvalD;
                bvalD += jumpD;
                qi = INT(bvalD) - 1;
                LIMIT_NUM(qi,0,(nr-1));
                qvalD = cvalsPO[qi].sc;
/*
printf("quant=%d r=%d bval=%f qi=%d qval=%f\n",quant,r,bvalD,qi,qvalD);
*/
                if(verbose) {
                    printf("# QuanTab col_%02d %s quant %d ind %d val <= ",c,cnameS,quant,qi);
                    printf(vformS,qvalD);
                    printf("\n");
                }
            }
            GetTableValI(tabPO,cvalsPO[r].id,c,&vD);
            vD = DNUM(quant);
/*
            SetTableValI(tabPO,cvalsPO[r].id,c,qvalD);
*/
            SetTableValI(tabPO,cvalsPO[r].id,c,vD);
        }
    }
    CHECK_SCOREC(cvalsPO);
    return(TRUE);
}
/**************************************************************************
*
*/

