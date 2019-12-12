/*
* blaststr.c
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
#include "prim.h"
#include "blastout.h"


/****************************************************************************
*   Create blastout structure
*/
BLASTOUT *CreateBlastoutPO()
{
    BLASTOUT *bPO;

    if(! (bPO = (BLASTOUT *)ALLOC(1,sizeof(BLASTOUT)) ) ) {
        printf("# Failed to allocate working BLASTOUT object\n"); 
        return(NULL);
    }
    bPO->ID = BLASTOUT_ID;
    InitBlastout(bPO);
    return(bPO);
}
/****************************************************************************
*   Free blastout structure
*/
int DestroyBlastoutI(BLASTOUT *bPO)
{
    VALIDATE(bPO,BLASTOUT_ID);
    CHECK_BLASTANS(bPO->ans);
    FREE(bPO);
    return(TRUE);
}
/****************************************************************************
*   Free blastout structure
*/
int AddHitSpaceI(int num,BLASTANS **aPPO)
{
    BLASTANS *aPO;

    if(num<1) {
        num = DEF_MAXHITS;
    }
    if( ! (aPO = CreateBlastansPO(num)) ) {
        printf("# Failed to allocate BLASTANS object (%d hits)\n",num); 
        return(FALSE);
    }
    *aPPO = aPO;
    return(TRUE);
}
/****************************************************************************
*   Create blastans structure
*/
BLASTANS *CreateBlastansPO(int max)
{
    BLASTANS *aPO;

    if(! (aPO = (BLASTANS *)ALLOC(1,sizeof(BLASTANS)) ) ) {
        printf("# failed to allocate working BLASTANS object\n"); 
        return(NULL);
    }
    aPO->ID = BLASTANS_ID;
    if(max <= 0) {
        max = DEF_MAXHITS;
    }
    aPO->asize = max;
    /* Filter qualifying mask */
    aPO->mask = (char *)ALLOC(max, sizeof(char));
    /***
    *   Hit info arrays
    */
    aPO->hits = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    aPO->scos = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    aPO->idens = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    if( (!aPO->mask) || (!aPO->hits) || (!aPO->scos) || (!aPO->idens) ) {
        printf("Failed to allocate hit space for %d\n",max);
        CHECK_BLASTANS(aPO);
        return(NULL);
    }
    aPO->qseqs = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    aPO->qhsc = (int *)ALLOC(max, sizeof(int));
    aPO->qhec = (int *)ALLOC(max, sizeof(int));
    if( (!aPO->qseqs) || (!aPO->qhsc) || (!aPO->qhec) ) {
        printf("Failed to allocate Q space for %d\n",max);
        CHECK_BLASTANS(aPO);
        return(NULL);
    }
    aPO->sseqs = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    aPO->shsc = (int *)ALLOC(max, sizeof(int));
    aPO->shec = (int *)ALLOC(max, sizeof(int));
    if( (!aPO->sseqs) || (!aPO->shsc) || (!aPO->shec) ) {
        printf("Failed to allocate S space for %d\n",max);
        CHECK_BLASTANS(aPO);
        return(NULL);
    }
    /* Alignment related numbers */
    aPO->alens = (int *)ALLOC(max, sizeof(int));
    aPO->amats = (int *)ALLOC(max, sizeof(int));
    aPO->amms = (int *)ALLOC(max, sizeof(int));
    aPO->agaps = (int *)ALLOC(max, sizeof(int));
    if( (!aPO->alens) || (!aPO->amats) || (!aPO->amms) || (!aPO->agaps) ) { 
        printf("Failed to allocate Align space for %d\n",max);
        CHECK_BLASTANS(aPO);
        return(NULL);
    }
    /* Processed hits */
    aPO->hlens = (int *)ALLOC(max, sizeof(int));
    aPO->hnums = (int *)ALLOC(max, sizeof(int));
    aPO->hfmb = (REAL *)ALLOC(max, sizeof(REAL));
    if( (!aPO->hlens) || (!aPO->hnums) || (!aPO->hfmb) ) {
        printf("Failed to allocate Hit space for %d\n",max);
        CHECK_BLASTANS(aPO);
        return(NULL);
    }
    InitBlastans(aPO,TRUE);
    return(aPO);
}
/****************************************************************************
*   Free blastans structure
*/
int DestroyBlastansI(BLASTANS *aPO)
{
    VALIDATE(aPO,BLASTANS_ID);
    CHECK_FILE(aPO->in);
    CHECK_FREE(aPO->mask);
    CHECK_FREE(aPO->hits);
    CHECK_FREE(aPO->scos);
    CHECK_FREE(aPO->idens);
    CHECK_FREE(aPO->qseqs);
    CHECK_FREE(aPO->qhsc);
    CHECK_FREE(aPO->qhec);
    CHECK_FREE(aPO->sseqs);
    CHECK_FREE(aPO->shsc);
    CHECK_FREE(aPO->shec);
    CHECK_FREE(aPO->hlens);
    CHECK_FREE(aPO->hnums);
    CHECK_FREE(aPO->hfmb);
    CHECK_FREE(aPO->alens);
    CHECK_FREE(aPO->amats);
    CHECK_FREE(aPO->amms);
    CHECK_FREE(aPO->agaps);
    FREE(aPO);
    return(TRUE);
}
/****************************************************************************
*   Init data structure; 
*/
void InitBlastout(BLASTOUT *bPO)
{
    bPO->owhat = 0;
    bPO->dhit = FALSE;
    bPO->dsum = FALSE;
    bPO->do_dfml = 0;
    bPO->dseq = FALSE;
    bPO->do_smu = FALSE;
    bPO->do_sml = FALSE;
    bPO->do_align = FALSE;
    bPO->dsco = FALSE;
    bPO->dhc = FALSE;
    bPO->dmbc = FALSE;
    bPO->do_con = FALSE;
    bPO->do_dci = FALSE;
    bPO->dhis = BOGUS;
    bPO->phis = BOGUS;
    bPO->chis = FALSE;
    bPO->mid = 0;
    bPO->mif = 0.0;
    bPO->do_mmm = -1;
    bPO->do_mgap = -1;
    bPO->do_fnot = FALSE;
    bPO->do_frq = FALSE;
    bPO->firstb = 0;
    bPO->lastb = TOO_BIG;
    bPO->rre = FALSE;
    bPO->do_rkc = FALSE;
    bPO->firstq = 0;
    bPO->lastq = TOO_BIG;
    bPO->firsth = 0;
    bPO->lasth = TOO_BIG;
    INIT_S(bPO->nsqh);
    INIT_S(bPO->opq);
    INIT_S(bPO->outname);
    bPO->out = NULL;
}
/****************************************************************************
*   Initialize blast answers struct
*/
void InitBlastans(BLASTANS *aPO, int full)
{
    if(full) {
        INIT_S(aPO->input);
        aPO->itype = 0;
    }
    INIT_S(aPO->query);
    aPO->qlen = 0;
    aPO->nhits = 0;
    aPO->nok = 0;
    aPO->maxhit = 0;
}
