/*
* stat.c
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
#include <math.h>
#include "prim.h"
#include "score.h"
#include "numlist.h"
#include "stat.h"

#define DB_STAT     if(DB[17])

/**************************************************************************
*   Pearson's correlation coefficient
*/
DOUB PearsonsCorD(void *aPO, void *bPO, int n, int ct)
{
    int i, *iPI,*jPI;
    REAL *iPR,*jPR;
    DOUB rD;
    DOUB *iPD,*jPD,syyD,sxyD,sxxD,axD,ayD,xtD,ytD;
    
    DB_STAT DB_PrI(">> PearsonsCorD, n=%d, ct=%d\n",n,ct);
    /***
    *   Get array averages
    */
    ArrayStatsI(aPO, ct, 0, n, NULL, NULL, &axD, NULL);
    ArrayStatsI(bPO, ct, 0, n, NULL, NULL, &ayD, NULL);
    DB_STAT DB_PrI("+ x %f, y %f\n",axD,ayD);
    /***
    *   Tally difference sums
    */
    syyD = sxyD = sxxD = 0.0;
    switch(ct)
    {
        case IS_INT:
            iPI = (int *)aPO;
            jPI = (int *)bPO;
            for(i=0; i<n; i++)
            {
                xtD = DNUM(iPI[i]) - axD;
                ytD = DNUM(jPI[i]) - ayD;
                sxxD += (xtD * xtD);
                syyD += (ytD * ytD);
                sxyD += (xtD * ytD);
            }
            break;
        case IS_REAL:
            iPR = (REAL *)aPO;
            jPR = (REAL *)bPO;
            for(i=0; i<n; i++)
            {
                xtD = DNUM(iPR[i]) - axD;
                ytD = DNUM(jPR[i]) - ayD;
                sxxD += (xtD * xtD);
                syyD += (ytD * ytD);
                sxyD += (xtD * ytD);
            }
            break;
        case IS_DOUB:
            iPD = (DOUB *)aPO;
            jPD = (DOUB *)bPO;
            for(i=0; i<n; i++)
            {
                xtD = iPD[i] - axD;
                ytD = jPD[i] - ayD;
                sxxD += (xtD * xtD);
                syyD += (ytD * ytD);
                sxyD += (xtD * ytD);
            }
            break;
        default:
            printf("Bogus ct code=%d\n",ct);    
            ERR("PearsonsCorD","Bogus ct code");
            return(-TOO_BIG_R);
    }
    DB_STAT DB_PrI("+ sx %f, sy %f, sxy\n",sxxD,syyD,sxyD);
    rD = DNUM(sxyD / (sqrt(sxxD*syyD) + PEARSON_TINY) );
    DB_STAT DB_PrI("<< PearsonsCorD %f\n",rD);
    return(rD);
}
/**************************************************************************
*   Rank correlation coefficient
*   ALLOCATES AND FREES TEMP ARRAYS HERE
*/
DOUB RankCorD(void *aPO, void *bPO, int n, int ct)
{
    int i, *iPI,*jPI;
    REAL *iPR,*jPR;
    DOUB rD,*iPD,*jPD;
    DOUB *rank1PD,*rank2PD;
    SCOREC *sc1PO, *sc2PO;
    
    DB_STAT DB_PrI(">> RankCorD, n=%d, ct=%d\n",n,ct);
    /***
    *   Attempt to allocate arrays
    */
    sc1PO = CreateScorecsPO(n);
    sc2PO = CreateScorecsPO(n);
    rank1PD = (DOUB *)ALLOC(n,sizeof(DOUB));
    rank2PD = (DOUB *)ALLOC(n,sizeof(DOUB));
    if( (!sc1PO) || (!sc1PO) || (!rank1PD) || (!rank2PD) )
    {
        printf("Problem allocating temp space for rank correlations!\n");
        CHECK_FREE(sc1PO); CHECK_FREE(sc2PO);
        CHECK_FREE(rank1PD); CHECK_FREE(rank2PD);
        return(-TOO_BIG_R);
    }
    /***
    *   Copy values from passed arrays SCORECs for sorting
    */
    switch(ct)
    {
        case IS_INT:
            iPI = (int *)aPO;
            jPI = (int *)bPO;
            for(i=0;i<n;i++)
            {
                sc1PO[i].sc = DNUM(iPI[i]);
                sc2PO[i].sc = DNUM(jPI[i]);
            }
            break;
        case IS_REAL:
            iPR = (REAL *)aPO;
            jPR = (REAL *)bPO;
            for(i=0;i<n;i++)
            {
                sc1PO[i].sc = DNUM(iPR[i]);
                sc2PO[i].sc = DNUM(jPR[i]);
            }
            break;
        case IS_DOUB:
            iPD = (DOUB *)aPO;
            jPD = (DOUB *)bPO;
            for(i=0;i<n;i++)
            {
                sc1PO[i].sc = iPD[i];
                sc2PO[i].sc = jPD[i];
            }
            break;
        default:
            CHECK_FREE(sc1PO); CHECK_FREE(sc2PO);
            CHECK_FREE(rank1PD); CHECK_FREE(rank2PD);
            printf("Bogus ct code=%d\n",ct);    
            ERR("RankCorD","Bogus ct code");
            return(-TOO_BIG_R);
    }
    /***
    *   Sort arrays and get rank arrays
    */
    SortScorecVals(sc1PO,n,SORT_ASCEND);
    SortScorecVals(sc2PO,n,SORT_ASCEND);
    FillScorecRankArrayI(sc1PO,n,rank1PD);
    FillScorecRankArrayI(sc2PO,n,rank2PD);
    /***
    *   Now get normal coorelation on rank arrays
    */
    rD = PearsonsCorD(rank1PD,rank2PD,n,IS_DOUB);
    /***
    *   Clean up and bail
    */
    CHECK_FREE(sc1PO); CHECK_FREE(sc2PO);
    CHECK_FREE(rank1PD); CHECK_FREE(rank2PD);
    return(rD);
}
/**************************************************************************
*   Fill valsPD with rank values from sorted SCOREC
*/
int FillScorecRankArrayI(SCOREC *scPO,int n,DOUB *valsPD)
{
    int i,j,tie;
    DOUB pvD,vD;

    i = 0;
    while(i<n)
    {
        tie = 1;
        pvD = scPO[i].sc;
        while( (scPO[i+tie].sc == pvD) && ((i+tie)<n) )
        {
            tie++;
        }
        vD = DNUM(i) + DNUM(tie)/2.0 + 0.5;
        for(j=0;j<tie;j++)
        {
            valsPD[scPO[i+j].id] = vD;
        }
        i += tie;
    }
    return(TRUE);
}
/***************************************************************************
*   Fill line with symbols to plot histogram
*   bval = bin count
*   ncum = cumulative bin count
*   max = max bin (i.e. mode) 
*   win = win = windows size (i.e. char width) 
*/
void FillHisPlotLine(int bval, int ncum, int max, int win, char *bufS)
{
    int i,j,n,m,clip;

    if(ncum > 0) {
        n = ROUND( DNUM(ncum)/DNUM(max) * DNUM(win) );
        m = ROUND( DNUM(bval)/DNUM(max) * DNUM(win) );
    }
    else {
        n = m = ROUND( DNUM(bval)/DNUM(max) * DNUM(win) );
    }
    /***
    *   If bar would be bigger than window size, clip
    */
    clip = (n>win) ? TRUE : FALSE;
    LIMIT_NUM(n,0,win);
    j = 0;
    if(bval>0) {
        bufS[j] = '|';
    }
    else {
        bufS[j] = ' ';
    }
    j++;
    for(i=0;i<n;i++)
    {
        if( (clip) && ( (i == (n-3)) || (i== (n-5)) ) ) {
            bufS[j] = '/';
        }
        else if( (clip) && (i == (n-4)) ) {
            bufS[j] = ' ';
        }
        else if(i>=m) {
            bufS[j] = '.';
        }
        else {
            bufS[j] = 'X';
        }
        j++;
    }
    bufS[j] = '\0';
}
/**************************************************************************
*   Given z-score, returns p-value in passed pointer
*   Modified from code from Hadar Avi-Itzhak
*/
int ProbabilityIntegralI(DOUB zD, DOUB *pPD) 
{
    int i;
    DOUB y,x,delta;

    y = x = 0.0;
    if(zD < 100.0) {
        delta = zD / (double)PROBINTSTEPS;
        for(i=0;i<PROBINTSTEPS;i++) 
        {
            y += ( exp(-x*x/2.0) +  exp(-(x+delta)*(x+delta)/2.0) );
            x += delta;
        }
        if(y<0.0) {
            y = 0.0;
        }
        y *=  delta / (2.0*sqrt(2.0*PI));
        y = 0.5 - y;
    }
    *pPD = y;
    return(TRUE);
}
/***************************************************************************
*   Make discrete value from any real
*/
int MakeDiscreteValI(DOUB buckD, DOUB *disPD)
{
    int n;
    DOUB disD,mulD;

    DB_STAT DB_PrI(">> MakeDiscreteValI buck=%f\n",buckD);
    /***
    *   Get multiplier to grow, shrink, or keep in single-digit range
    */
    disD = buckD;
    if(disD < 1.0) {
        mulD = 10.0;
    }
    else if(disD > 10.0) {
        mulD = 0.1;
    }
    else {
        mulD = 1.0;
    }
    /***
    *   Loop until in single digit range
    */
    n = 0;
    while( (disD<1.0) || (disD>10.0) ) {
        disD = disD * mulD;
        n++;
    }
    /***
    *   Set discrete value
    */
    if(disD < 2.0) {
        disD = 1.0;
    }
    else if(disD < 5.0) {
        disD = 2.0;
    }
    else {
        disD = 5.0;
    }
    DB_STAT DB_PrI("+ mul=%f dis=%f n=%d\n",mulD,disD,n);
    /***
    *   Back to original scale
    */
    while(n>0) {
        disD = disD / mulD;
        n--;
    }
    BOG_CHECK(!disPD);
    *disPD = disD;
    DB_STAT DB_PrI("<< MakeDiscreteValI %f\n",disD);
    return(TRUE);
}
/*************************************************************************
*   Figure out discrete bin and start values
*/
int GetDiscreteBinStartI(DOUB binD, DOUB *binPD, DOUB loD, DOUB *loPD)
{
    DOUB mD;

    MakeDiscreteValI(binD,&binD);
    mD = DNUM( INT(loD / binD) ) * binD;
    if(mD > loD) {
        mD -= binD;
    }
    loD = mD;
    if(binPD) {
        *binPD = binD;
    }
    if(loPD) {
        *loPD = loD;
    }
    return(TRUE);
}
/**************************************************************************
*   Pick "natural" histogram bin
*/
int NumlistNaturalHistBinI(NUMLIST *nlPO, int max_bin, int dhis, 
    DOUB *binPD, DOUB *lowPD, DOUB *hiPD)
{
    int n;
    DOUB binD,minD,maxD,min_difD;

    DB_STAT DB_PrI(">> NumlistNaturalHistBinI max=%d, dhis=%d\n",max_bin,dhis);
    VALIDATE(nlPO,NUMLIST_ID);
/*
    n = GetNumlistLengthI(nlPO);
    if(n<2) {
        return(FALSE);
    }
*/
    if(!NumlistStatsI(nlPO,-1,-1, &minD, &maxD, NULL, NULL)) {
        DB_STAT DB_PrI("<< NumlistNaturalHistBinI Stats failed FALSE\n");
        return(FALSE);
    }
    if(!NumlistDifStatsI(nlPO, &min_difD, NULL, NULL)) {
        DB_STAT DB_PrI("<< NumlistNaturalHistBinI DifStats failed FALSE\n");
        return(FALSE);
    }
    /***
    *   Flat = one bin, else min diff
    */
    if( minD == maxD ) {
        binD = 1.0;
    }
    else {
        binD = min_difD;
    }
    DB_STAT DB_PrI("+ min=%f max=%f min_dif=%f bin=%f\n",minD,maxD,min_difD,binD);
    /***
    *   Shrink if too many bins
    */
    n = INT( (maxD-minD)/binD);
    DB_STAT DB_PrI("+ n=%d ... ",n);
    if(max_bin>0) {
        while(n > max_bin)
        {
            binD = binD * 2.0;
            n = n/2;
        }
    }
    DB_STAT DB_PrI("n=%d bin=%f\n",n,binD);
    /***
    *   Discrete values?
    */
    if(dhis) {
        DB_STAT DB_PrI("+ calling to make bin and min discrete\n");
        GetDiscreteBinStartI(binD, &binD, minD, &minD);
    }
    n = INT( (maxD-minD)/binD);
    maxD = minD + DNUM(n) * binD;
    DB_STAT DB_PrI("+ n=%d bin=%f min=%f max=%f\n",n,binD,minD,maxD);
    if(binPD) {
        *binPD = binD;
    }
    if(lowPD) {
        *lowPD = minD;
    }
    if(hiPD) {
        *hiPD = maxD;
    }
    DB_STAT DB_PrI("<< NumlistNaturalHistBinI %d\n",n);
    return(n);
}
/***************************************************************************
*   Pick "natural" histogram bin based on minium spacing between values,
*       truncated to max wanted bins
*   The maximum number of bins is specified by max
*   If dhis is true, discrete-bounded histogram bins are chosen
*   Sets starting bin value via lowPD
*   Sets bin size value via binPD
*   Returns number of different values
*/
int NaturalHistBinI(DOUB *valsPD, int num, int max, int dhis,
    DOUB *lowPD, DOUB *binPD)
{
    int i,n;
    DOUB pD,minD,maxD,binD,nminD;
    SCOREC *scPO;

    DB_STAT DB_PrI(">> NaturalHistBinI, num=%d, max=%d, dhis=%d\n",num,max,dhis);
    /***
    *   Init
    */
    if(lowPD) {
        *lowPD = BAD_D;
    }
    if(binPD) {
        *binPD = BAD_D;
    }
    if( (num<1) || (!(scPO=CreateScorecsPO(num))) ) {
        DB_STAT DB_PrI("<< NaturalHistBinI FALSE\n");
        return(FALSE);
    }
    /***
    *   Fill values then sort if needed
    */
    minD = TOO_BIG_R;
    maxD = -TOO_BIG_R;
    for(i=0;i<num;i++)
    {
        minD = MIN_NUM(valsPD[i], minD);
        maxD = MAX_NUM(valsPD[i], maxD);
        scPO[i].sc = valsPD[i];
    }
    /***
    *   Find min between-value spacing 
    */
    if( minD == maxD ) {
        binD = 1.0;
    }
    else {
        SortScorecVals(scPO,num,SORT_ASCEND);
        binD = maxD - minD;
        pD = scPO[0].sc;    
        for(i=1; i<num; i++) 
        {
            /***
            *   If current != previous, bin = save min diff
            */
            if(scPO[i].sc != pD) {
                binD = MIN_NUM(binD, (scPO[i].sc - pD) );
                pD = scPO[i].sc;
            }
        }
    }
    CHECK_SCOREC(scPO);
    DB_STAT DB_PrI("+ min=%f max=%f bin=%f\n",minD,maxD,binD);
    /***
    *   Shrink if too many bins
    */
    n = INT( (maxD-minD)/binD);
    if(max>0) {
        while(n>max)
        {
            binD = binD * 2.0;
            n = n/2;
        }
    }
    /***
    *   Discrete values
    */
    if(dhis) {
        MakeDiscreteValI(binD,&binD);
        nminD = DNUM( INT(minD/binD) ) * binD;
        if(nminD > minD) {
            nminD -= binD;
        }
        minD = nminD;
    }

    DB_STAT DB_PrI("+ min=%f bin=%f\n",minD,binD);
    if(lowPD) {
        *lowPD = minD;
    }
    if(binPD) {
        *binPD = binD;
    }
    DB_STAT DB_PrI("<< NaturalHistBinI %d\n",n);
    return(n);
}
/*************************************************************************
*   Find stats for between-value differences 
*/
int NumlistDifStatsI(NUMLIST *nlPO, DOUB *minPD, DOUB *maxPD, DOUB *avPD)
{
    int len,i,nd;
    DOUB v1D,v2D, delD, minD, maxD, avD;
    NUMLIST *newPO;

    VALIDATE(nlPO,NUMLIST_ID);
    len = GetNumlistLengthI(nlPO);
    if(len<2) {
        PROBLINE;
        printf("No dif stats for array with %d elements!\n",len);
        return(FALSE);
    }
    /***
    *   Duplicate so we can sort without messing up given order
    */
    if(!(newPO = DuplicateNumlistPO(nlPO,-1))) {
        PROBLINE;
        printf("Failed to duplicate numbers for dif stats\n");
        return(FALSE);
    }
    SortNumlistI(newPO, SORT_ASCEND);
    /***
    *   Init, then party through rest of list looking at difs
    */
    minD = TOO_BIG_D;
    maxD = avD = 0.0;
    GetNumlistIntDoubI(newPO,0,NULL,&v1D);
    nd = 0;
    for(i=1;i<len;i++) 
    {
        GetNumlistIntDoubI(newPO,i,NULL,&v2D);
        delD = v2D - v1D;
        if(v2D != v1D) {
            nd++;
            if(delD < minD) {
                minD = delD;
            }
            if(delD > maxD) {
                maxD = delD;
            }
        }
        avD += delD;
        v1D = v2D;
    }
    /***
    *   Special case where all the same, so no difs ever evaluated
    *   Set min and max to 0; Av = any value = last one (v2)
    */
    if(nd < 1) {
        minD = maxD = 0.0;
        avD = v2D;
    }
    if(minPD) {
        *minPD = minD;
    }
    if(maxPD) {
        *maxPD = maxD;
    }
    if(avPD) {
        if(nd > 0) {
            avD = avD / DNUM(len);
        }
        *avPD = avD;
    }
    CHECK_NUMLIST(newPO);
    return(len);
}
/***************************************************************************
*   Numlist interface to raw coorelation functions   
*/
int NumlistPairCorI(NUMLIST *flPO, NUMLIST *slPO, int rank, DOUB *rPD)
{
    int ok, n, *fPI, *sPI;
    DOUB rD, *fPD, *sPD;

    DB_STAT DB_PrI(">> NumlistPairCorI %p v %p rank=%d\n",flPO,slPO,rank);
    VALIDATE(flPO,NUMLIST_ID);
    VALIDATE(slPO,NUMLIST_ID);
    if( !NumlistAreSameI(flPO, slPO) ) {
        DB_STAT DB_PrI("<< NumlistPairCorI Numlists differ FALSE\n");
        return(FALSE);
    }
    /***
    *   Get pointers to raw arrays and coorelate these
    */
    ok = FALSE;
    if(flPO->type == IS_INT) {
        GetNumlistPtrIntsI(flPO, &fPI, &n);
        GetNumlistPtrIntsI(slPO, &sPI, NULL);
        if(rank) {
            rD = RankCorD(fPI, sPI, n, IS_INT);
        }
        else {
            rD = PearsonsCorD(fPI, sPI, n, IS_INT);
        }
        ok++;
    }
    else if(flPO->type == IS_DOUB) {
        GetNumlistPtrDoubsI(flPO, &fPD, &n);
        GetNumlistPtrDoubsI(slPO, &sPD, NULL);
        if(rank) {
            rD = RankCorD(fPD, sPD, n, IS_DOUB);
        }
        else {
            rD = PearsonsCorD(fPD, sPD, n, IS_DOUB);
        }
        ok++;
    } 
    if( ok && rPD) {
        *rPD = rD;
    }
    DB_STAT DB_PrI("<< NumlistPairCorI %d\n",ok);
    return(ok);
}
