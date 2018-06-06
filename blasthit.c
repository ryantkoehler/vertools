/*
* blasthit.c
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

#define DB_BIO  if(DB[133])

/****************************************************************************
*   Process raw hit list, setting mask with OK passing ones
*/
void ProcHitList(BLASTOUT *bPO, BLASTANS *aPO)
{
    int i,nok,ok,num,max,mc,den,alen;
    REAL fracR;

    DB_BIO DB_PrI(">> ProcHitList\n");
    /* Init counts and mask array */
    nok = max = 0;
    InitArrayI(aPO->mask, IS_CHAR, 0, aPO->asize, 0);
    for(i=0; i<aPO->nhits; i++) 
    {
        if( ((i+1) < bPO->firsth) || ((i+1) > bPO->lasth) ) {
            continue;
        }
        SetAlignCountsI(bPO, aPO, i);
        /* Get what qualifies as match count */
        num = HitMatchCountI(bPO, aPO, i, &den);
        max = MAX_NUM(num,max);
        /* Qualify flag */
        ok = TRUE;
        /* Min match */
        if( ok && (num < bPO->mid) ) {
            ok = FALSE;
        }
        /* Min fraction match */
        alen = AdjustHitLenI(bPO,aPO,i, den);
        fracR = RNUM(num)/RNUM(alen);
        if( ok && (fracR < bPO->mif) ) {
            ok = FALSE;
        }
        /* Max mismatch */
        if( ok && (bPO->do_mmm >= 0) ) {
            /* Mismatch count is total or just in alignment */
            if(bPO->do_frq) {
                /* Max of explicit mismatch or len - match */
                mc = MAX_NUM(aPO->qlen - aPO->amats[i], aPO->amms[i]);
            }
            else {
                mc = aPO->amms[i];
            }
            if(mc > bPO->do_mmm) {
                ok = FALSE;
            }
        }
        /* Max gaps */
        if( ok && (bPO->do_mgap >= 0) ) {
            if(aPO->agaps[i] > bPO->do_mgap) {
                ok = FALSE;
            }
        }
        /* Logic inversion; Not */
        ok = bPO->do_fnot ? (!ok) : ok;
        /* Count and mark */
        if (ok) {
            nok++;
            aPO->mask[i] = TRUE;
        }
    }
    aPO->nok = nok;
    aPO->maxhit = max;
    DB_BIO DB_PrI("<< ProcHitList\n");
}
/****************************************************************************
*   Fill histogram with matching base counts
*/
int FillHitHistI(BLASTOUT *bPO, BLASTANS *aPO)
{
    int i,num,max;

    InitArrayI(aPO->ihist, IS_INT, 0, HISTDIM, 0);
    max = 0;
/**
printf("nhits=%d\n",aPO->nhits);
*/
    for(i=0;i<aPO->nhits;i++)
    {
        num = HitMatchCountI(bPO, aPO, i, NULL);
        max = MAX_NUM(num,max);
        if(num >= HISTDIM) {
            continue;
        }
        aPO->ihist[num] += 1;
/**
printf("i=%d num=%d ihis=%d max=%d\n",i,num,aPO->ihist[num],max);
**/
    }
    return(max);
}
/****************************************************************************
*   Fill histogram with matching base counts
*/
void IntegrateHist(BLASTANS *aPO)
{
    int i,num;

    num = 0;
    for(i=HISTDIM-1;i>0;i--)
    {
        num += aPO->ihist[i];
        aPO->ihist[i] = num;
    }
}
/****************************************************************************
*   Returns the count of (qualified) matching bases; Can be contig, 3' or total
*   If last arg, set this to (length) denominator 
*/
int HitMatchCountI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int *denPI)
{
    int i, n, c, mat;
    char qseqS[BLBSIZE+1], sseqS[BLBSIZE+1], maskS[BLBSIZE];

    sprintf(qseqS,"%s",&aPO->qseqs[hit*BLBSIZE]);
    sprintf(sseqS,"%s",&aPO->sseqs[hit*BLBSIZE]);
    n = FillHitAlignQseqMaskI(bPO, aPO, hit, maskS);
    mat = 0;
    /* Contig 3' end */
    if( bPO->do_co3 ) {
        /***
        *   Has to be full length; Short circuit loop if not
        */
        n = (aPO->qhec[hit] == aPO->qlen) ? n : 0;
        /* Loop backwards 3' end first */
        for(i=n-1;i>=0;i--)
        {
            if(qseqS[i] != sseqS[i]) {
                break;
            }
            mat++;
        }
    }
    /* General max contig */
    else if( bPO->do_con ) {
        c = 0;
        for(i=0;i<n;i++)
        {
            if(qseqS[i] != sseqS[i]) {
                mat = MAX_NUM(c,mat);
                c = 0;
            }
            else {
                c++;
            }
        }
        mat = MAX_NUM(c,mat);
    }
    /* total matches */
    else {
        mat = aPO->amats[hit];
    }
    /* Denominator? set to len */
    if(denPI) {
        if(bPO->do_frq) {
            *denPI = aPO->qlen;
        }
        else {
            *denPI = n;
        }
    }
    return(mat);
}
/**************************************************************************
*   Count contiguous matches from 3' end
*/
int Cont3pEndMatchI(BLASTANS *aPO, int len, char *qPC, char *sPC)
{
    int n,i;

    n = 0;
    i = len;
/*
printf("xxx i=%d len=%d %p %p\n",i,len,qPC,sPC);
*/
    while(i>=0) {
/*
        printf("\t%c %c %d\n", qPC[i], sPC[i], i);
*/
        if(qPC[i] != sPC[i]) {
            break;
        }
        n++;
        i--;
    }
    return(n);
}
/**************************************************************************
*   Fill passed map array with query coords for hit alignment;
*   Basically, count along seq but skip gaps (i.e. '-' in aligment)
*   Return length of array actually filled
*/
int FillHitAlignQcoordMapI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int *mapPI)
{
    int i,b,n;
    char qseqS[BLBSIZE+1];

    sprintf(qseqS,"%s",&aPO->qseqs[hit*BLBSIZE]);
    if(bPO->rre) {
        b = aPO->qhec[hit];
        i = strlen(qseqS) - 1;
    }
    else {
        b = aPO->qhsc[hit];
        i = 0;
    }
    n = 0;
    while(isgraph(qseqS[i]))
    {
        mapPI[i] = b;
        if(bPO->rre) {
            if(qseqS[i] != '-') {
                b--;
            }
            i--;
        }
        else {
            if(qseqS[i] != '-') {
                b++;
            }
            i++;
        }
        n++;
    }
    return(n);
}
/**************************************************************************
*   Fill passed char mask to mark in-range alignment query string 
*/
int FillHitAlignQseqMaskI(BLASTOUT *bPO, BLASTANS *aPO, int hit, char *maskPC)
{
    int i, n, b, qclis[BLBSIZE];
    char qseqS[BLBSIZE+1];

    /* Init full mask */
    sprintf(qseqS,"%s",&aPO->qseqs[hit*BLBSIZE]);
    n = strlen(qseqS);
    InitArrayI(maskPC, IS_CHAR, 0, n, TRUE);
    /* If base range, unmask out-of-range positions */
    if(bPO->firstb > 0) {
        FillHitAlignQcoordMapI(bPO, aPO, hit, qclis);
        for(i=0; i<n; i++)
        {
            /* Current base from end or start? */
            if(bPO->rre) {
                b = aPO->qlen - qclis[i] + 1;
            }
            else {
                b = qclis[i];
            }
            if( (b < bPO->firstb) || (b > bPO->lastb) ) {
                maskPC[i] = FALSE;
            }
        }
    }
    return(n);
}
/**************************************************************************
*   Set various counts for hit alignment 
*/
int SetAlignCountsI(BLASTOUT *bPO, BLASTANS *aPO, int hit)
{
    int i, n, nmat, nmm, ngap;
    char qseqS[BLBSIZE+1], sseqS[BLBSIZE+1], maskS[BLBSIZE];

    if(hit >= aPO->asize) {
        return(FALSE);
    }
    sprintf(qseqS,"%s",&aPO->qseqs[hit*BLBSIZE]);
    sprintf(sseqS,"%s",&aPO->sseqs[hit*BLBSIZE]);
    /***
    *   Get mask of query positions to look at, then count
    */
    FillHitAlignQseqMaskI(bPO, aPO, hit, maskS);
    n = nmat = nmm = ngap = 0;
    for(i=0;i<strlen(qseqS);i++)
    {
        if(! maskS[i]) {
            continue;
        }
        n++;
        if(sseqS[i] == qseqS[i]) {
            nmat++;
        }
        else {
            nmm++;
            if( (sseqS[i] == '-') || (qseqS[i] == '-') ) {
                ngap++;
            }
        }
    }
    aPO->alens[hit] = n;
    aPO->amats[hit] = nmat;
    aPO->amms[hit] = nmm;
    aPO->agaps[hit] = ngap;
    return(TRUE);
}
/****************************************************************************
*   Set alignment length for normalization
*/
int AdjustHitLenI(BLASTOUT *bPO, BLASTANS *aPO, int hit, int len)
{
    int qs,qe;

    if(bPO->do_frq) {
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

