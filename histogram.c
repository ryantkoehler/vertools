/*
* histogram.c
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
* See https://www.verdascend.com/ for more
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prim.h"
#include "numlist.h"
#include "stat.h"
#include "histogram.h"

#define DB_HHI      if(DB[35])
#define DB_HLOW     if(DB[36])

/**************************************************************************/
HISTOGRAM *CreateHistogramPO(NUMLIST *valsPO)
{
    int ok;
    HISTOGRAM *hisPO;

    DB_HHI DB_PrI(">> CreateHistogramPO vals %p\n",valsPO);
    if(!(hisPO=(HISTOGRAM *)ALLOC(1,sizeof(HISTOGRAM)))) {
        return(NULL);
    }
    hisPO->ID = HISTOGRAM_ID;
    InitHistogram(hisPO);
    ok = TRUE;
    if(valsPO) {
        ok = SetHistogramValueSourceI(hisPO, valsPO);
        if(ok) {
            ok = SetUpHistogramBinningI(hisPO, -1, -1);
        }
    }
    if(!ok) {
        CHECK_HISTOGRAM(hisPO);
        hisPO = NULL;
    }
    DB_HHI DB_PrI("<< CreateHistogramPO %p\n",hisPO);
    return(hisPO);
}
/*************************************************************************/
int DestroyHistogramI(HISTOGRAM *hisPO)
{
    DB_HHI DB_PrI(">> DestroyHistogramI %p\n",hisPO);
    VALIDATE(hisPO,HISTOGRAM_ID);
    CHECK_NUMLIST(hisPO->bins);
    FREE(hisPO);
    DB_HHI DB_PrI("<< DestroyHistogramI TRUE\n");
    return(TRUE);
}
/*************************************************************************/
void InitHistogram(HISTOGRAM *hisPO)
{
    DB_HHI DB_PrI(">> InitHistogram %p\n",hisPO);
    VALIDATE(hisPO,HISTOGRAM_ID);
    hisPO->n = 0;
    hisPO->nov = hisPO->nun = 0;
    hisPO->bsize = BAD_D;
    hisPO->lo = hisPO->hi = BAD_D;
    hisPO->val_n = 0;
    hisPO->val_lo = hisPO->val_hi = BAD_D;
    DB_HHI DB_PrI("<< InitHistogram\n");
    return;
}
/*************************************************************************/
void DumpHistogram(HISTOGRAM *hisPO, int bins, FILE *outPF)
{
    VALIDATE(hisPO,HISTOGRAM_ID);
    HAND_NFILE(outPF);
    fprintf(outPF,"Histogram at %p\n",hisPO);
    fprintf(outPF,"Source values at %p\n",hisPO->vals);
    fprintf(outPF,"Values from %5.4f to %5.4f\n",hisPO->val_lo,hisPO->val_hi);
    fprintf(outPF,"Histogram from %5.4f to %5.4f, %5.4f\n",hisPO->lo,hisPO->hi,hisPO->ran);
    fprintf(outPF,"Histogram bin size %5.4f\n",hisPO->bsize);
    fprintf(outPF,"Histogram bin number %d\n",hisPO->n);
    fprintf(outPF,"Histogram bins at %p\n",hisPO->bins);
    if(hisPO->bins) {
        if(bins) {
            DumpNumlist(hisPO->bins,-1,-1,outPF);
        }
        else {
            DumpNumlist(hisPO->bins,0,0,outPF);
        }
    }
    return;
}
/*************************************************************************/
int SetHistogramValueSourceI(HISTOGRAM *hisPO, NUMLIST *valsPO)
{
    int ok;

    DB_HHI DB_PrI(">> SetHistogramValueSourceI his %p numlist %p\n",hisPO,valsPO);
    VALIDATE(hisPO,HISTOGRAM_ID);
    VALIDATE(valsPO,NUMLIST_ID);
    ok = NumlistStatsI(valsPO, -1, -1, &hisPO->val_lo, &hisPO->val_hi, NULL,NULL);
    if(ok) {
        hisPO->vals = valsPO;
        hisPO->val_n = GetNumlistLengthI(valsPO);
        hisPO->val_ran = hisPO->val_hi - hisPO->val_lo;
        DB_HHI DB_PrI("+ vals %p n=%d, %f to %f r=%f\n",hisPO->val_n,hisPO->val_lo,hisPO->val_hi,hisPO->val_ran);
    }
    DB_HHI DB_PrI("<< SetHistogramValueSourceI %d\n",ok);
    return(ok);
}
/**************************************************************************/
int SetHistogramParamsI(HISTOGRAM *hisPO, DOUB bsizeD, DOUB loD, DOUB hiD)
{
    DB_HHI DB_PrI(">> SetHistogramParamsI b=%f l=%f h=%f\n",bsizeD,loD,hiD);
    VALIDATE(hisPO,HISTOGRAM_ID);
    if(bsizeD > 0.0) {
        hisPO->bsize = bsizeD;
    }
    if(!BAD_DOUB(loD)) {
        hisPO->lo = loD;
    }
    if(!BAD_DOUB(hiD)) {
        hisPO->hi = hiD;
    }
    hisPO->ran = hisPO->hi - hisPO->lo;
    LIMIT_NUM(hisPO->ran, 0.0, TOO_BIG_D);
    DB_HHI DB_PrI("+ bsize=%f lo=%f hi=%f ran=%f\n",hisPO->bsize,hisPO->lo,hisPO->hi,hisPO->ran);
    DB_HHI DB_PrI("<< SetHistogramParamsI TRUE\n");
    return(TRUE);
}
/*************************************************************************/
int GetHistogramParamsI(HISTOGRAM *hisPO, DOUB *binPD, DOUB *loPD, DOUB *hiPD)
{
    VALIDATE(hisPO,HISTOGRAM_ID);
    if(binPD) {
        *binPD = hisPO->bsize;
    }
    if(loPD) {
        *loPD = hisPO->lo;
    }
    if(hiPD) {
        *hiPD = hisPO->hi;
    }
    return(TRUE);
}
/***************************************************************************/
int GetHistogramNumBinsI(HISTOGRAM *hisPO) 
{
    VALIDATE(hisPO,HISTOGRAM_ID);
    return(hisPO->n);
}
/***************************************************************************/
int SetUpHistogramBinningI(HISTOGRAM *hisPO, int max, int dis)
{
    int n;
    DOUB loD,hiD;

    DB_HHI DB_PrI(">> SetUpHistogramBinningI max=%d dis=%d\n",max,dis);
    VALIDATE(hisPO,HISTOGRAM_ID);
    /***
    *   Make sure hist has values and these are set
    */
    if(!hisPO->vals) {
        DB_HHI DB_PrI("<< SetUpHistogramBinningI No vals FALSE\n");
        return(FALSE);
    }
    if(!SetHistogramValueSourceI(hisPO, hisPO->vals)) {
        DB_HHI DB_PrI("<< SetUpHistogramBinningI setting source failed FALSE\n");
        return(FALSE);
    }
    /***
    *   Auto figure bin size; Also maybe adjust starting point loD?
    */
    loD = hisPO->val_lo;
    hiD = hisPO->val_hi;
    if(hisPO->bsize < 0.0) {
        max = (max<0) ? HIS_NUM_BINS : max;
        dis = (dis<0) ? TRUE : dis;
        DB_HHI DB_PrI("+ figuring bin size; max=%d dis=%d\n",max,dis);
        NumlistNaturalHistBinI(hisPO->vals, max, dis, &hisPO->bsize,&loD,&hiD);
        DB_HHI DB_PrI("+ natural sizes bin=%f lo=%f hi=%f\n",hisPO->bsize,loD,hiD);
    }
    if(BAD_DOUB(hisPO->lo)) {
        if(BAD_DOUB(loD)) {
            loD = hisPO->val_lo;
        }
        hisPO->lo = loD;
        DB_HHI DB_PrI("+ set lo=%f\n",hisPO->lo);
    }
    if(BAD_DOUB(hisPO->hi)) {
        hisPO->hi = hiD;
        DB_HHI DB_PrI("+ set hi=%f\n",hisPO->hi);
    }
/*
    hisPO->ran = hisPO->hi - hisPO->lo;
*/
    hisPO->ran = hisPO->hi - hisPO->lo + TINY_R;
    LIMIT_NUM(hisPO->ran, 0.0, TOO_BIG_D);
    DB_HHI DB_PrI("+ lo=%f hi=%f ran=%f\n",hisPO->lo,hisPO->hi,hisPO->ran);
    /***
    *   Figure the number of bins and make/set space for these 
    *   If all bin size is less or equal to range, bump up number of bins
    */
    n = INT( hisPO->ran / hisPO->bsize );
/*
    if((DNUM(n) * hisPO->bsize) <= hisPO->ran) { 
*/
    if( ( DNUM(n) * hisPO->bsize) <= hisPO->ran ) {
        n++;
    }
    DB_HHI DB_PrI("+ n bins = %d = (%f / %f)\n",n,hisPO->ran,hisPO->bsize);
    if(n > hisPO->n) {
        DB_HHI DB_PrI("+ need space, creating new numlist\n");
        CHECK_NUMLIST(hisPO->bins);
        if(!(hisPO->bins = CreateNumlistPO(IS_INT, NULL, 0))) {
            return(FALSE);
        }
    }
    DB_HHI DB_PrI("+ setting bin len=%d\n",n);
    SetNumlistLengthI(hisPO->bins, n);
    SetNumlistRangeIntsI(hisPO->bins, 0, n, 0);
    hisPO->n = n;
    hisPO->hi = hisPO->lo + DNUM(n) * hisPO->bsize;
    DB_HHI DB_PrI("<< SetUpHistogramBinningI n=%d TRUE\n",hisPO->n);
    return(TRUE);
}
/***************************************************************************
*   Tally bin counts for histogram; Must have values and bins set already
*/
int TallyHistogramI(HISTOGRAM *hisPO)
{
    int i,b,v;
    DOUB vD,stD,enD;
    
    DB_HLOW DB_PrI("\n>> TallyHistogramI\n");
    VALIDATE(hisPO,HISTOGRAM_ID);
    /***
    *   Make sure have values and bins
    */
    if( (!hisPO->vals) || (!hisPO->bins) ) {
        DB_HLOW DB_PrI("<< TallyHistogramI no vals or bins FALSE\n");
        return(FALSE);
    }
    SetNumlistRangeIntsI(hisPO->bins, -1, -1, 0);
    hisPO->nov = hisPO->nun = 0;
    /***
    *   Go through values calculating bins and updating counts
    */
    for(i=0;i<hisPO->val_n;i++) 
    {
        GetNumlistIntDoubI(hisPO->vals,i,NULL,&vD);
        HistogramBinForValueI(hisPO, vD, &b, &stD, &enD);
        DB_HLOW DB_PrI("+ [%d] v=%f => %f %f",i,vD,stD,enD);
        if(b<0) {
            hisPO->nun += 1;
            DB_HLOW DB_PrI(" un  %d\n",hisPO->nun);
        }
        else if(b >= hisPO->n) {
            hisPO->nov += 1;
            DB_HLOW DB_PrI(" ov  %d\n",hisPO->nov);
        }
        else {
            GetNumlistIntI(hisPO->bins,b,&v);
            v += 1;
            DB_HLOW DB_PrI(" [%d] = %d\n",b,v);
            SetNumlistIntI(hisPO->bins,b,v);
        }
    }
    DB_HLOW DB_PrI("<< TallyHistogramI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Counts of values under / over range limits
*/
int GetHistogramUnOvCountsI(HISTOGRAM *hisPO, int *unPI, int *ovPI)
{
    VALIDATE(hisPO,HISTOGRAM_ID);
    if(!hisPO->bins) {
        return(FALSE);
    }
    if(unPI) {
        *unPI = hisPO->nun;
    }
    if(ovPI) {
        *ovPI = hisPO->nov;
    }
    return(TRUE);
}
/***************************************************************************
*   Given histogram bin index, get values: Bin count, start, end
*/
int HistogramValuesForBinI(HISTOGRAM *hisPO, int bin, int *numPI, DOUB *stPD, 
    DOUB *enPD)
{
    DOUB stD, enD;
    int n;

    VALIDATE(hisPO,HISTOGRAM_ID);
    if(!hisPO->bins) {
        return(FALSE);
    }
    /***
    *   Handle special under/over end cases or actually in hist 
    */
    if(bin < 0) {
        stD = hisPO->val_lo;
        enD = hisPO->lo;
        n = hisPO->nun;
    }
    else if(bin >= hisPO->n) {
        stD = hisPO->hi;
        enD = hisPO->val_hi;
        n = hisPO->nov;
    }
    else {
        stD = hisPO->lo + DNUM(bin) * hisPO->bsize;
        enD = stD + hisPO->bsize;
        GetNumlistIntI(hisPO->bins, bin, &n);
    }
    if(numPI) {
        *numPI = n;
    }
    if(stPD) {
        *stPD = stD;
    }
    if(enPD) {
        *enPD = enD;
    }
    return(TRUE); 
}
/***************************************************************************
*   Given value, figure histogram bin index, start, end
*/
int HistogramBinForValueI(HISTOGRAM *hisPO, DOUB vD, int *bPI, DOUB *stPD, DOUB *enPD)
{
    int b,ok;

    DB_HLOW DB_PrI(">> HistogramBinForValueI %f\n",vD);
    VALIDATE(hisPO,HISTOGRAM_ID);
    if(hisPO->bsize <= 0.0) {
        return(FALSE);
    }
    ok = TRUE;
    if(vD < hisPO->lo) {
        b = -1;
    }
    else {
        b = INT( (vD - hisPO->lo) / hisPO->bsize );
    }
    if(bPI) {
        *bPI = b;
    }
    if(stPD || enPD) {
        ok = HistogramValuesForBinI(hisPO, b, NULL, stPD, enPD);
    }
    DB_HLOW DB_PrI("<< HistogramBinForValueI (b=%d) %d\n",b,ok);
    return(ok);
}
/***************************************************************************
*   Get format string for histogram based on bin and end values
*/
int HistogramAutoFormatStringI(HISTOGRAM *hisPO, char *formS)
{
    DOUB valsD[3];

    VALIDATE(hisPO,HISTOGRAM_ID);
    if(!hisPO->bins) {
        return(FALSE);
    }
    /***
    *   Want all three, as bins can be small (precision) but ends big (width)
    */
    valsD[0] = hisPO->bsize;
    valsD[1] = hisPO->lo;
    valsD[2] = hisPO->hi;
    DoubArrayPrecisionI(valsD, 3, NULL, NULL, formS);
    return(TRUE);
}
/***************************************************************************
*   Find maximum and next-highest count bin values for histogram
*/
int GetHistogramMaxTwoBinsI(HISTOGRAM *hisPO, int *maxPI, int *smaxPI) 
{
    int i,max,smax,v1,v2;

    VALIDATE(hisPO,HISTOGRAM_ID);
    max = smax = 0;
    /***    
    *   Max
    */
    GetHistogramUnOvCountsI(hisPO, &v1, &v2);
    max = MAX_NUM(v1,v2);
    for(i=0;i<hisPO->n;i++)
    {
        HistogramValuesForBinI(hisPO, i, &v1, NULL, NULL); 
        max = MAX_NUM(v1,max);
    }
    /***
    *   Second pass for second max
    */
    GetHistogramUnOvCountsI(hisPO, &v1, &v2);
    if(v1 < max) {
        smax = MAX_NUM(smax,v1);
    }
    if(v2 < max) {
        smax = MAX_NUM(smax,v2);
    }
    for(i=0;i<hisPO->n;i++)
    {
        HistogramValuesForBinI(hisPO, i, &v1, NULL, NULL); 
        if(v1 < max) {
            smax = MAX_NUM(v1,smax);
        }
    }
    if(maxPI) {
        *maxPI = max;
    }
    if(smaxPI) {
        *smaxPI = smax;
    }
    return(TRUE);
}
/***************************************************************************
*   Get histogram structure from number list
*/
int NumlistToHistogramI(NUMLIST *nlPO, DOUB binD, DOUB loD, DOUB hiD, 
    HISTOGRAM **hisPPO)
{
    HISTOGRAM *hisPO;

    DB_HHI DB_PrI(">> NumlistToHistogramI bin=%f lo=%f hi=%f\n",binD,loD,hiD);
    VALIDATE(nlPO,NUMLIST_ID);
/*
    if(!(hisPO = CreateHistogramPO(nlPO))) {
*/
    if(!(hisPO = CreateHistogramPO(NULL))) {
        DB_HHI DB_PrI("<< NumlistToHistogramI no hist = FALSE\n");
        return(FALSE);
    }
    SetHistogramParamsI(hisPO,binD,loD,hiD);
    SetHistogramValueSourceI(hisPO, nlPO);
    SetUpHistogramBinningI(hisPO,-1,-1);
    TallyHistogramI(hisPO); 
    *hisPPO = hisPO;
    DB_HHI DB_PrI("<< NumlistToHistogramI his=%p TRUE\n",hisPO);
    return(TRUE);
}

