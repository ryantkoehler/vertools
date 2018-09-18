/*
* sctab_ga.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prim.h"
#include "numlist.h"
#include "wordlist.h"
#include "table.h"
#include "scoretab.h"


/**************************************************************************
*   Handle cross over  
*/
int HandleSctGACrossOverI(SCORETAB *stPO)
{
    int r,c;
    char seedS[DEF_BS];
    DOUB vD;
    TABLE *p1PO, *p2PO;

    /***
    *   Get other parent to merge
    */
    p1PO = stPO->tab;
    p2PO = NULL;
    if(!GetTableI(stPO->mergname,stPO->rlab,stPO->clab,stPO->corn,stPO->do_skp,&p2PO)) {
        PROBLINE;
        printf("Failed to get table for cross over\n");
        return(FALSE);
    }
    if(!SameDimTablesI(p1PO, p2PO, FALSE, TRUE)) {
        PROBLINE;
        printf("Cross over tables don't match\n");
        CHECK_TABLE(p2PO);
        return(FALSE);
    }
    /***
    *   Get random masks ...
    */
    r = c = 0;
    if(stPO->gaxr > 0.0) {
        r = MaskRandSubsetI(p2PO->rmask, p2PO->nrow, stPO->gaxr);
    }
    if(stPO->gaxc > 0.0) {
        c = MaskRandSubsetI(p2PO->cmask, p2PO->ncol, stPO->gaxc);
    }
    printf("# GS Cross over wtih %s\n", stPO->mergname);
    FillRandSeedString(stPO->seed,seedS);
    printf("#    Random seed: %s\n",seedS);
    printf("#    Rows: %f = %d\n",stPO->gaxr,r);
    printf("#      RowMask: "); 
    DumpArray(p2PO->rmask, IS_CHAR, 0,p2PO->nrow, "%d ", NULL);
    printf("\n");
    printf("#    Cols: %f = %d\n",stPO->gaxc,c);
    printf("#      ColMask: "); 
    DumpArray(p2PO->cmask, IS_CHAR, 0,p2PO->ncol, "%d ", NULL);
    printf("\n");
    /***
    *   Copy ...
    */
    if(stPO->gaxr > 0.0) {
        for(r=0; r<p2PO->nrow; r++) 
        {
            if(!p2PO->rmask[r]) {
                continue;
            }
            for(c=0; c<p2PO->ncol; c++) 
            {
                GetTableValI(p2PO,r,c,&vD);
                SetTableValI(p1PO,r,c,vD);
            }
        }
    }
    if(stPO->gaxc > 0.0) {
        for(c=0; c<p2PO->ncol; c++) 
        {
            if(!p2PO->cmask[c]) {
                continue;
            }
            for(r=0; r<p2PO->nrow; r++) 
            {
                GetTableValI(p2PO,r,c,&vD);
                SetTableValI(p1PO,r,c,vD);
            }
        }
    }
    CHECK_TABLE(p2PO);
    return(TRUE);
}
/**************************************************************************
*   Handle cross over  
*/
int HandleSctGAMutationI(SCORETAB *stPO)
{
    int r,c,n;
    char seedS[DEF_BS];
    DOUB vD,gD;
    TABLE *tabPO;
    
    tabPO = stPO->tab;
    printf("# GS mutation frequency %f +/- gaussian s.d. %f\n",stPO->gamf,stPO->gamg);
    FillRandSeedString(stPO->seed,seedS);
    printf("#    Random seed: %s\n",seedS);
    /***
    *
    */
    n = 0;
    for(r=0; r<tabPO->nrow; r++) 
    {
        for(c=0; c<tabPO->ncol; c++) 
        {
            if(RandD(1.0) >= stPO->gamf) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            gD = RandGaussD(1.0, stPO->gamg);
            vD *= gD;
            SetTableValI(tabPO,r,c,vD);
            n++;
        }
    }
    printf("#    Total elements tweaked: %d\n",n);
    return(TRUE);
}
