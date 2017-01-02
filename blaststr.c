/*
* blaststr.c
*
* Copyright 2017 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
    CHECK_BLASTANS(bPO->sans);
    CHECK_BLASTANS(bPO->mans);
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
    /***
    *   Hit info arrays
    */
    aPO->hits = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    aPO->scos = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    aPO->idens = (char *)ALLOC(max * BLBSIZE, sizeof(char));
    if( (!aPO->hits) || (!aPO->scos) || (!aPO->idens) ) {
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
    aPO->hlens = (int *)ALLOC(max, sizeof(int));
    aPO->hnums = (int *)ALLOC(max, sizeof(int));
    aPO->hfmb = (REAL *)ALLOC(max, sizeof(REAL));
    if( (!aPO->hlens) || (!aPO->hnums) || (!aPO->hfmb) ) {
        printf("Failed to allocate x space for %d\n",max);
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
    bPO->do_smc = FALSE;
    bPO->do_sml = FALSE;
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
    bPO->do_mfq = FALSE;
    bPO->do_mnot = FALSE;
    bPO->firstb = 0;
    bPO->lastb = TOO_BIG;
    bPO->rre = FALSE;
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
void InitBlastans(BLASTANS *aPO,int full)
{
    if(full) {
        INIT_S(aPO->input);
        aPO->itype = 0;
    }
    INIT_S(aPO->query);
    aPO->qlen = 0;
    aPO->nhits = 0;
    aPO->maxhit = 0;
}
/****************************************************************************
*   Set alignment length for normalization 
*/
int AdjustHitLenI(BLASTOUT *bPO,BLASTANS *aPO,int hit, int len)
{
    int qs,qe;

    if(bPO->do_mfq) {
        len = aPO->qlen;
    }
    else if(bPO->firstb > 0) {
        qs = MIN_NUM(aPO->qhsc[hit], aPO->qhec[hit]);
        qe = MAX_NUM(aPO->qhsc[hit], aPO->qhec[hit]);
        if(bPO->rre) {
            if( qe < (aPO->qlen - bPO->firstb + 1) ) {
                len += (aPO->qlen - bPO->firstb + 1 - qe);
            }
            if( qs > (aPO->qlen - bPO->lastb + 1) ) {
                len += (aPO->qlen - bPO->lastb + 1 - qs);
            }
        }
        else {
            if( qs > bPO->firstb ) {
                len += (qs - bPO->firstb);
            }
            if( qe < bPO->lastb ) {
                len += (bPO->lastb - qe);
            }
        }
    }
    return(len);
}
/****************************************************************************
*   Get query portion of sequence for hit 
*/
int GetBlastansQseqI(BLASTANS *aPO,int hit,char *seqS,int max)
{
    VALIDATE(aPO,BLASTANS_ID);
    if( (hit<0) || (hit>=aPO->nhits) )
    {
        return(FALSE);
    }
    LIMIT_NUM(max,1,BLBSIZE-1);
    strncpy(seqS,&aPO->qseqs[hit*BLBSIZE],max);
    seqS[max] =  '\0';
    return(TRUE);
}
/****************************************************************************
*   Get subject portion of sequence for hit 
*/
int GetBlastansSseqI(BLASTANS *aPO,int hit,char *seqS,int max)
{
    VALIDATE(aPO,BLASTANS_ID);
    if( (hit<0) || (hit>=aPO->nhits) )
    {
        return(FALSE);
    }
    LIMIT_NUM(max,1,BLBSIZE-1);
    strncpy(seqS,&aPO->sseqs[hit*BLBSIZE],max);
    seqS[max] =  '\0';
    return(TRUE);
}
