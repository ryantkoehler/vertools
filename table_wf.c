/*
* table_wf.c
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
#include "table.h"

#define DB_WFM if(DB[73])

/****************************************************************************
*   Perform weak-first mapping algorithm with values in tabPO; Rows are 
*       mapped to columns, with high scores winning (maximize not minimize)
*   Mapping results are given in mapPO, as valid array indices; For 
*       non-overloaded case (ncol >= nrow) the values are 0; For overloaded
*       cases, the values are the paritioned set numbers (order dependent for
*       tie scores), or "pool numbers".
*   If mask is TRUE, non-masked rows / columns will be ignored
*   Poolsize should be <= number of columns; if negative it is set to this
*   To have periodic update of cycle number, use verb > 0
*/
int WeakFirstMappingI(TABLE *tabPO, TABLE *mapPO, int mask, int psize, int verb)
{
    int r,c,v,cmax,nrow,ncol,rreal,creal,pnum,cycle,weak,best,vnum;
    NUMLIST *rmaskPO, *cmaskPO, *pmaxPO;
    DOUB vD,bestD,worstD;
    char rnameS[DEF_BS],cnameS[DEF_BS];
    FILE *oPF;

    DB_WFM DB_PrI(">> WeakFirstMappingI mask=%d psize=%d verb=%d\n",mask,
        psize, verb);
    /***
    *   SHAM xxx; someday this should be a real file?
    */
    oPF = NULL;
    HAND_NFILE(oPF);
    /***
    *   Check that tables are cool and same size
    */
    VALIDATE(tabPO,TABLE_ID);
    VALIDATE(mapPO,TABLE_ID);
    if(!SameDimTablesI(tabPO,mapPO,mask,TRUE)) {
        PROBLINE;
        printf("Differing diminsions for value and mapping tables:\n");
        return(FALSE);
    }
    DB_WFM {
        DumpTable(tabPO,TRUE,FALSE,stdout);
        fflush(stdout);
    }
    /***
    *   Count rows and cols (masked for real counts here)
    *   Check that there's a possible outcome
    */
    rreal = GetTableRowsI(tabPO,mask);
    creal = GetTableColsI(tabPO,mask);
    if( (rreal<1) || (creal<1) ) {
        PROBLINE;
        printf("Too few rows / cols (mask=%d)\n",mask);
        printf(" N-rows = %d   N-cols = %d\n",rreal,creal);
        return(FALSE);
    }
    if( creal<psize ) {
        PROBLINE;
        printf("Too few cols for pool size (mask=%d)\n",mask);
        printf(" N-cols = %d   Psize = %d\n",creal,psize);
        return(FALSE);
    }
    DB_WFM DB_PrI("+ rreal=%d creal=%d\n",rreal,creal);
    /***
    *   Count again without masking for accounting arrays 
    */
    nrow = GetTableRowsI(tabPO,FALSE);
    ncol = GetTableColsI(tabPO,FALSE);
    DB_WFM DB_PrI("+ nrow=%d ncol=%d\n",nrow,ncol);
    /***
    *   Check psize, set column max 
    */
    if(psize < 0) {
        psize = creal;
    }
    cmax = INT(ceil(DNUM(rreal)/DNUM(psize)));
    DB_WFM DB_PrI("+ rreal=%d psize=%d cmax=%d\n",rreal,psize,cmax);
    /***
    *   Allocate working arrays
    */
    rmaskPO = CreateNumlistPO(IS_INT, NULL, nrow);
    cmaskPO = CreateNumlistPO(IS_INT, NULL, ncol);
    pmaxPO = CreateNumlistPO(IS_INT, NULL, cmax);
    if( (!rmaskPO) || (!cmaskPO) || (!pmaxPO) ) {
        PROBLINE;
        printf("Failed to allocate working masks for weak-first mapping\n");
        CHECK_NUMLIST(rmaskPO); 
        CHECK_NUMLIST(cmaskPO); 
        CHECK_NUMLIST(pmaxPO);
        return(FALSE);
    }
    /***
    *   Mapping table gets bogus to indicate not mapped
    *   Rows and cols all get zero to indicate they can be used
    */
    InitTableValsI(mapPO,-1.0,mask);
    SetNumlistRangeIntsI(rmaskPO,-1,-1,0);
    SetNumlistRangeIntsI(cmaskPO,-1,-1,0);
    SetNumlistRangeIntsI(pmaxPO,-1,-1,0);
    DB_WFM DB_PrI("+ initialized nrow=%d ncol=%d cmax=%d\n",nrow,ncol,cmax);
    /***
    *   Iterate over all rows
    */
    if(verb>0) {
        fprintf(oPF,"#\n");
        fprintf(oPF,"# Starting weak-first mapping on %d rows to %d cols\n",
            rreal,creal);
        fprintf(oPF,"#  Poolsize = %d\n",psize);
        fflush(oPF);
    }
    /***
    *   Cycle to find weakest row at each iteration
    */
    vnum = rreal;
    cycle = 0;
    while(rreal > 0)
    {
        /***
        *   Feedback story
        */
        cycle++;
        if(verb > 0)
        {
            if( (cycle%verb) == 0 )
            {
                fprintf(oPF,"#  W-F cycle %4d  %6.2f%% done\n",
                    cycle, PERCENT_R(vnum-rreal,vnum));
                fflush(oPF);
            }
        }
        DB_WFM DB_PrI("+ %d rows remaining\n",rreal);
        /***
        *   Find current, non-mapped, non-masked weakest row; 
        */
        weak = -1;
        worstD = TOO_BIG_R;
        for(r=0;r<tabPO->nrow;r++)
        {
            /***
            *   If using masking and this row not set = ignore
            */
            if( (mask) && (!tabPO->rmask[r]) )
            {
                continue;
            }
            /***
            *   Already mapped, so continue
            */
            GetNumlistIntI(rmaskPO,r,&v);
            if(v) {
                continue;
            }
            /***
            *   Look at all non-masked, non-used column values for this row 
            *   to find best, which is the row's score
            */
            bestD = -TOO_BIG_R;
            for(c=0;c<tabPO->ncol;c++)
            {
                /***
                *   If masking and this col is unset = ignore
                */
                if( (mask) && (!tabPO->cmask[c]) ) {
                    continue;
                }
                /***
                *   Column maxxed out already, so continue
                */
                GetNumlistIntI(cmaskPO, c, &v);
                if( v >= cmax ) {
                    continue;
                }
                /***
                *   Remember best value in this row
                */
                GetTableValI(tabPO,r,c,&vD);
                if(vD > bestD) {
                    bestD = vD;
                }
            }
            /***
            *   Is this the row the weakest (row best less than current worst)?
            */
            if(bestD < worstD) {
                weak = r;
                worstD = bestD;
            }
        }
        BOG_CHECK(weak<0);
        DB_WFM {
            GetTableRowLabI(tabPO,weak,rnameS,DEF_BS-1);
            DB_PrI("+   Weakest=%d\t%s\twith %f\n",weak,rnameS,worstD);
        }
        /***
        *   Now go back and map the best column for weakest row
        */
        best = -1;
        bestD = -TOO_BIG_R;
        for(c=0;c<tabPO->ncol;c++)
        {
            if( (mask) && (!tabPO->cmask[c]) ) {
                continue;
            }
            GetNumlistIntI(cmaskPO, c, &v);
            if( v >= cmax) {
                continue;
            }
            GetTableValI(tabPO,weak,c,&vD);
            if(vD > bestD) {
                best = c;
                bestD = vD;
            }
        }
        BOG_CHECK(best<0);
        /***
        *   Assign mapping of r-to-c
        *   The actual number assigned is the (first allowed) pool number.
        */
        DB_WFM {
            GetTableColLabI(tabPO,best,cnameS,DEF_BS-1);
            DB_PrI("+   Best=%d\t%s\twith %f\n",best,cnameS,bestD);
        }
        GetNumlistIntI(cmaskPO, best, &pnum);
        DB_WFM DB_PrI("+   Pool=%d",pnum);
        GetNumlistIntI(pmaxPO, pnum, &v);
        while( (v == psize) && (pnum < cmax) )
        {
            pnum++;
            GetNumlistIntI(pmaxPO, pnum, &v);
        }
        SetNumlistIntI(cmaskPO, best, pnum);
        DB_WFM DB_PrI(" bumped=%d\n",pnum);
        vD = DNUM(pnum);
        SetTableValI(mapPO,weak,best,vD);
        DB_WFM DB_PrI("+   Mapped R=%d C=%d (%s to %s) with %1.0f\n",
            weak,best,rnameS,cnameS,vD);
        /***
        *   Update accounting arrays
        *   Row has been mapped, so ignore later
        *   Column count updated 
        *   Pool count updated
        */
        SetNumlistIntI(rmaskPO, weak, TRUE);
        GetNumlistIntI(cmaskPO, best, &v);
        SetNumlistIntI(cmaskPO, best, (v+1) );
        GetNumlistIntI(pmaxPO, pnum, &v);
        SetNumlistIntI(pmaxPO, pnum, (v+1) );
        rreal--;
    }   
    /***
    *   Final feedback / debug?
    */
    if(verb > 0)
    {   
        fprintf(oPF,"#\n");
        fflush(oPF);
    }
    DB_WFM 
    {
        DumpTable(mapPO,TRUE,FALSE,stdout);
        fflush(stdout);
    }
    /***
    *   All done
    */
    CHECK_NUMLIST(rmaskPO); 
    CHECK_NUMLIST(cmaskPO); 
    CHECK_NUMLIST(pmaxPO);
    DB_WFM DB_PrI("<< WeakFirstMappingI TRUE\n");
    return(TRUE);
}
