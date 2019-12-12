/*
* sctab_sc.c
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
#include <math.h>
#include "prim.h"
#include "score.h"
#include "table.h"
#include "scoretab.h"


/**************************************************************************
*   Handle score transforms 
*/
int HandleScoreTransformI(SCORETAB *stPO, TABLE *tabPO)
{
    int r,c;
    DOUB vD,sD;

    c = MaskCountI(tabPO->cmask,tabPO->ncol);
    /***
    *   If not globally applying, must have scores for each column
    */
    if( (!stPO->do_scg) && (stPO->nscores<c) ) {
        PROBLINE;
        printf(
          "Table has %d (unmasked) columns but only %d scores were loaded\n",
                c,stPO->nscores);
        printf("NOT TRANSFORMING VALUES\n");
        return(FALSE);
    }
    /***
    *   Go through table transforming values
    */
    for(r=0;r<tabPO->nrow;r++) {
        if(!tabPO->rmask[r]) {
            continue;
        }
        for(c=0;c<tabPO->ncol;c++) {
            if(!tabPO->cmask[c]) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            if(stPO->do_scg) {
                sD = DNUM(EvalScfieldScoreR(stPO->scores[0],vD));
            }
            else {
                sD = DNUM(EvalScfieldScoreR(stPO->scores[c],vD));
            }
            SetTableValI(tabPO,r,c,sD);
        }
    }
    return(TRUE);
}
