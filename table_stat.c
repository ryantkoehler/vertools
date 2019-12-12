/*
* table_stat.c
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prim.h"
#include "table.h"

#define DB_TAB      if(DB[70])
#define DB_TAB_LO   if(DB[71])

/****************************************************************************/
int TableColStatsI(TABLE *tabPO, int col, int mask, DOUB *loPD, DOUB *hiPD, 
    DOUB *avPD, DOUB *sdPD)
{
    return(TableStatsI(tabPO,0,tabPO->nrow,col,col+1,mask,loPD,hiPD,avPD,sdPD));
}
/****************************************************************************/
int TableRowStatsI(TABLE *tabPO, int row, int mask, DOUB *loPD, DOUB *hiPD, 
    DOUB *avPD, DOUB *sdPD)
{
    return(TableStatsI(tabPO,row,row+1,0,tabPO->ncol,mask,loPD,hiPD,avPD,sdPD));
}
/****************************************************************************
*   Get stats for table 
*   May be subset of rows rs to re and cols cs to ce unless these are < 0, 
*       then full table is scanned
*/
int TableStatsI(TABLE *tabPO, int rs, int re, int cs, int ce, int mask,
    DOUB *loPD, DOUB *hiPD, DOUB *avPD, DOUB *sdPD)
{
    int r,c,n;
    DOUB vD,avD,loD,hiD,sdD;

    VALIDATE(tabPO,TABLE_ID);
    /***
    *   If bogus bounds, get full table
    */
    if(rs<0) {
        rs = 0;
    }
    if(re<0) {
        re = tabPO->nrow;
    }
    if(cs<0) {
        cs = 0;
    }
    if(ce<0) {
        ce = tabPO->ncol;
    }
    /***
    *   Scan table
    */
    loD = TOO_BIG_R;
    hiD = -TOO_BIG_R;
    avD = 0;
    n = 0;
    for(r=rs;r<re;r++) {
        if( mask && (!tabPO->rmask[r]) ) {
            continue;
        }
        for(c=cs;c<ce;c++) {
            if( mask && (!tabPO->cmask[c]) ) {
                continue;
            }
            if(!GetTableValI(tabPO,r,c,&vD)) {
                return(FALSE);
            }
            loD = MIN_NUM(vD,loD);
            hiD = MAX_NUM(vD,hiD);
            avD += vD;
            n++;
        }
    }
    if(n<1) {
        return(n);
    }
    avD /= DNUM(n);
    /***
    *   Set all but sd here
    */
    if(loPD) {
        *loPD = loD;
    }
    if(hiPD) {
        *hiPD = hiD;
    }
    if(avPD) {
        *avPD = avD;
    }
    if(!sdPD) {
        return(n);
    }
    /***
    *   Easy case with no variance, SD = 0
    */
    if(loD == hiD) {
        *sdPD = 0.0;
        return(n);
    }
    /***
    *   Pass two for sd
    */
    sdD = 0.0;
    for(r=rs;r<re;r++) {
        if( mask && (!tabPO->rmask[r]) ) {
            continue;
        }
        for(c=cs;c<ce;c++) {
            if( mask && (!tabPO->cmask[c]) ) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            sdD += ( (vD-avD) * (vD-avD) );
        }
    }
    sdD /= DNUM(n);
    *sdPD = DNUM(sqrt(sdD));    
    return(n);
}
