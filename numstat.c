/*
* numstat.c
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
#include <string.h>
#include <ctype.h>
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "numlist.h"
#include "histogram.h"
#include "table.h"
#include "score.h"
#include "stat.h"
#include "numstat.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(NumstatI(argc,argv),NULL) ); }
/**************************************************************************/
void NumstatUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("    <infile>   File containing data\n");
    printf("    -out XXX   Output file XXX\n");
    printf("    -col #     Take data from column #; Default is column %d\n",DEF_COL);
    printf("    -scol #    Second data col #\n");
    printf("    -dif       Report difference [scol] - [col]\n");
    printf("    -ic #      Ignore characters up to position # on data lines\n");
    printf("    -sk        Skip bogus / missing data; Default is abort\n");
    printf("    -hb #      Histogram data in bins of #\n");
    printf("    -his # #   Histogram data in bins of # starting at #\n");
    printf("    -hlr # #   Histogram limit reported range # to #\n");
    printf("    -hne       Histogram no end reporting (i.e. Under / Over range)\n");
    printf("    -hmb #     Histogram max bins # (for auto); Default = %d\n",DEF_HMB);
    printf("    -hp -pi    Histogram plot / plot with integral\n");
    printf("    -hwd #     Histogram plot width # (in chars); Default = %d\n",HPLOT_SLEN);
    printf("    -hntb      Histogram No Truncating Bins; Default is clip single big one\n");
    printf("    -htxf #    Histogram Truncate (single dominant) if > #-fold next bin\n");
    printf("    -sp        Scatter plot with second col data\n");
    printf("    -prs #     Percentile series in steps of #\n");
    printf("    -prl XXX   Percentile list; comma delimited as \"10,50,90\"\n");
    printf("    -sl        Single line output\n");
    printf("    -efi       Echo filename in output\n");
    printf("    -echo XXX  Echo string XXX in output\n");
    printf("\n");
    return;
}
/**************************************************************************/
int NumstatI(int argc, char **argv)
{
    NUMSTAT *nsPO;
    int num;
    DOUB *rvalsPD;

    nsPO = CreateNumstatPO();
    if(!ParseArgsI(argc,argv,
        "S -sl B -his D2 -echo S -efi B -col I -out S -hlr D2\
        -sk B -ic I -hp B -pi B -hb D -prs D\
        -prl S -hmb I -hwd I -hne B -scol I -sp B -hntb B -htxf D -dif B",
        nsPO->inname, &nsPO->do_sl, &nsPO->h_bin,&nsPO->h_lo, nsPO->echo, 
        &nsPO->do_efi, &nsPO->col,
        nsPO->outname, &nsPO->h_lo,&nsPO->h_hi, &nsPO->do_sk, &nsPO->do_ic, 
        &nsPO->do_hplot, &nsPO->do_ploti, &nsPO->h_bin, &nsPO->prsd, nsPO->prls, 
        &nsPO->h_maxbin, &nsPO->h_pwide, &nsPO->do_hends, 
        &nsPO->scol, &nsPO->do_splot, &nsPO->do_hntb, &nsPO->htb_xfold, &nsPO->do_dif,
        (int *)NULL))
    {
        NumstatUse();
        CHECK_NUMSTAT(nsPO); 
        return(TRUE);
    }
    /***
    *   Arg consistency check 
    */
    if(!CheckNumstatOptionsI(nsPO)) {
        CHECK_NUMSTAT(nsPO); 
        return(FALSE);
    }
    /***
    *   Files
    */
    if(!OpenNumstatFilesI(nsPO)) {
        CHECK_NUMSTAT(nsPO); 
        return(FALSE);
    }
    /*****
    *   Load data lines
    */
    if(!LoadColValsFromFileI(nsPO)) {
        CHECK_NUMSTAT(nsPO); 
        return(FALSE);
    }
    if(!nsPO->do_sl) {
        fprintf(nsPO->out,"%s# %d lines read from %s, %d values\n", nsPO->echo,
                nsPO->lines,nsPO->inname,nsPO->num);
        fprintf(nsPO->out,"%s# (column %d)\n",nsPO->echo, nsPO->col);
    }
    /***
    *   Doing difference? Set these values
    */
    if(nsPO->do_dif) {
        CalcColDifsI(nsPO);
    }
    /***
    *   Get pointer to raw double array and length. 
    *   This is easier to use with array functions; Also sort for percentile
    */
    GetNumlistPtrDoubsI(nsPO->vals, &rvalsPD, &num);
    /***
    *   Call to get stats
    *   If only single-line output and not percentiles, report and bail here
    */
    NumlistStatsI(nsPO->vals, -1, -1, &nsPO->min,&nsPO->max,&nsPO->av,&nsPO->sd);
    NumlistSumI(nsPO->vals, -1, -1, &nsPO->sum);
    if( (nsPO->do_sl) && (!nsPO->do_perc) ) {
        NumstatOneLineOut(nsPO);
        CHECK_NUMSTAT(nsPO); 
        return(TRUE);
    }
    /***
    *   Automatic histogram best binning (report even if no hist)
    */
    HandleHisMaxBinsI(nsPO);
    NumlistNaturalHistBinI(nsPO->vals, nsPO->h_maxbin, TRUE, &nsPO->h_abin, 
        &nsPO->h_alo, &nsPO->h_ahi);
    /***
    *   Header story
    */
    if( !nsPO->do_sl ) {
        NumstatHandleHeader(nsPO);
    }
    /***
    *   Percentile shams? 
    */
    if(nsPO->do_perc) {
        SortArray(rvalsPD,IS_DOUB,num,SORT_ASCEND);
        if(!NO_S(nsPO->prls)) {
            HandleListPercentiles(nsPO, rvalsPD, num, nsPO->out);
        }
        if(!BAD_DOUB(nsPO->prsd)) {
            HandleStepPercentiles(nsPO, rvalsPD, num, nsPO->out);
        }
    }
    /***
    *   Doing histogram?
    */
    if(nsPO->do_hist) {
        NumstatHandleHistI(nsPO);
    }
    if(nsPO->do_splot) {
        NumstatHandleSplotI(nsPO);
    }
    CHECK_NUMSTAT(nsPO);
    return(TRUE);
}
/*************************************************************************/
int NumstatHandleHistI(NUMSTAT *nsPO)
{
    int i,n,ncum,nun,nov,nmax,m1,m2;
    DOUB stD,enD,modeD;
    HISTOGRAM *hisPO;

    /***
    *   If not specified, set binning and bounds to auto values
    */
    if(BAD_DOUB(nsPO->h_bin)) {
        nsPO->h_bin = nsPO->h_abin;
    }
    if(BAD_DOUB(nsPO->h_lo)) {
        nsPO->h_lo = nsPO->h_alo;
    }
    if(BAD_DOUB(nsPO->h_hi)) {
        nsPO->h_hi = nsPO->h_ahi;
    }
    NumlistToHistogramI(nsPO->vals, nsPO->h_bin, nsPO->h_lo, nsPO->h_hi, &hisPO);
    /***
    *   Get max total value or max his bin to normalize 
    *   Also under and over counts 
    */
    if(nsPO->do_ploti) {
        NumlistSumI(hisPO->bins, -1, -1, &modeD);
        nmax = INT(modeD);
    }
    else {
        NumlistStatsI(hisPO->bins, -1, -1, NULL,&modeD,NULL,NULL);
        nmax = INT(modeD);
    }
    GetHistogramUnOvCountsI(hisPO, &nun, &nov);
    if(nsPO->do_hends) {
        nmax = MAX_NUM(nmax ,nun);
        nmax = MAX_NUM(nmax ,nov);
    }
    /***
    *   Check if one bin dominates and truncate back if allowed
    */
    if( (! nsPO->do_hntb) && (hisPO->n > 1) ) {
        GetHistogramMaxTwoBinsI(hisPO, &m1, &m2);
        if ( DNUM(m1) > (DNUM(m2) * nsPO->htb_xfold) )  {
            nmax = ROUND( DNUM(m2) * nsPO->htb_xfold );
        }
    }
    /***
    *   Dump histogram values; Under, body, Over
    */
    NumstatHistHeader(nsPO, hisPO, nsPO->out);
    ncum = nun;
    if(nsPO->do_hends && (nun > 0)) {
        HistogramValuesForBinI(hisPO,-1,NULL,NULL,&enD);
        NumPrecisionI(enD, NULL,NULL,NULL, nsPO->h_pfmt);
        OneHistLineNumsOut(nsPO, nun, enD, "<", ncum);
        OneHistLinePlotOut(nsPO, nun, ncum, nmax);
        fprintf(nsPO->out,"\n");
    } 
    HistogramAutoFormatStringI(hisPO, nsPO->h_pfmt);
    for(i=0; i < hisPO->n; i++)
    {
        HistogramValuesForBinI(hisPO,i,&n,&stD,&enD);
        ncum += n;
        OneHistLineNumsOut(nsPO, n, stD, NULL, ncum);
        OneHistLinePlotOut(nsPO, n, ncum, nmax);
        fprintf(nsPO->out,"\n");
    }
    if(nsPO->do_hends && (nov > 0)) {
        HistogramValuesForBinI(hisPO,hisPO->n,NULL,&stD,NULL);
        NumPrecisionI(stD, NULL,NULL,NULL, nsPO->h_pfmt);
        ncum += nov;
        OneHistLineNumsOut(nsPO, nov, stD, ">", ncum);
        OneHistLinePlotOut(nsPO, nov, ncum, nmax);
        fprintf(nsPO->out,"\n");
    } 
    CHECK_HISTOGRAM(hisPO);
    return(TRUE);
}
/**************************************************************************/
void OneHistLineNumsOut(NUMSTAT *nsPO, int n, DOUB stD, char *bmS, int ncum )
{
    fprintf(nsPO->out,"%s",nsPO->echo);
    if(bmS && (!NO_S(bmS))) {
        fprintf(nsPO->out,"%s",bmS);
    }
    fprintf(nsPO->out,nsPO->h_pfmt,stD);
    fprintf(nsPO->out,"%s",nsPO->hvsep);
    fprintf(nsPO->out,"%-6d",n);
    fprintf(nsPO->out,"%s",nsPO->hvsep);
    fprintf(nsPO->out,"%6.3f", PERCENT_R(n,nsPO->num));
    fprintf(nsPO->out,"%s",nsPO->hvsep);
    fprintf(nsPO->out,"%7.3f", PERCENT_R(ncum,nsPO->num));
}
/**************************************************************************/
void OneHistLinePlotOut(NUMSTAT *nsPO, int n, int ncum, int nmax)
{
    char bufS[BBUFF_SIZE];

    if( nsPO->do_hplot | nsPO->do_ploti ) {
        if(nsPO->do_ploti) {
            FillHisPlotLine(n, ncum, nmax, nsPO->h_pwide, bufS);
        }
        else {
            FillHisPlotLine(n, BOGUS, nmax, nsPO->h_pwide, bufS);
        }
        fprintf(nsPO->out," %s",bufS);
    }
    return;
}
/**************************************************************************/
void NumstatHistHeader(NUMSTAT *nsPO, HISTOGRAM *hisPO, FILE *outPF)
{
    NumstatHistStatHeader(nsPO, hisPO, "Histogram in bins of ", FALSE, outPF);
    NumstatHistLabHeader(nsPO, outPF);
}
/**************************************************************************/
void NumstatHistStatHeader(NUMSTAT *nsPO, HISTOGRAM *hisPO, char *shamS, int two, 
    FILE *outPF)
{
    DOUB binD,loD,hiD;
    char formS[DEF_BS];

    HAND_NFILE(outPF);
    GetHistogramParamsI(hisPO, &binD, &loD, &hiD);
    HistogramAutoFormatStringI(hisPO,formS);
    fprintf(outPF,"%s#\n",nsPO->echo);
    if(two) {
        fprintf(outPF,"%s# %s", nsPO->echo,shamS);
        fprintf(outPF,formS,binD); 
        fprintf(outPF,"\n");
        fprintf(outPF,"%s# Reported range ",nsPO->echo);
    }
    else {
        fprintf(outPF,"%s# %s",nsPO->echo,shamS);
        fprintf(outPF,formS,binD); 
        fprintf(outPF," from ");
    }
    fprintf(outPF,formS,loD); 
    fprintf(outPF," to "); 
    fprintf(outPF,formS,hiD); 
    fprintf(outPF,"\n");
}
/**************************************************************************/
void NumstatHistLabHeader(NUMSTAT *nsPO, FILE *outPF)
{
    fprintf(outPF,"%s#\n",nsPO->echo);
    fprintf(outPF,"%s# Bin",nsPO->echo);
    fprintf(outPF,"%s",nsPO->hvsep);
    fprintf(outPF,"Number");
    fprintf(outPF,"%s",nsPO->hvsep);
    fprintf(outPF,"Perc");
    fprintf(outPF,"%s",nsPO->hvsep);
    fprintf(outPF,"cPerc");
    if( nsPO->do_hplot || nsPO->do_ploti ) {
        fprintf(outPF,"%s",nsPO->hvsep);
        fprintf(nsPO->out,"Plot");
    }
    fprintf(outPF,"\n");
    return;
}
/*************************************************************************
*   Scatter plot
*/
int NumstatHandleSplotI(NUMSTAT *nsPO)
{
    int i,r,c,nrow,ncol;
    DOUB fD,sD;
    char labS[DEF_BS],numS[DEF_BS],formS[DEF_BS];
    HISTOGRAM *fhisPO, *shisPO;
    TABLE *tabPO;

    NumlistToHistogramI(nsPO->vals, -1.0, BAD_D, BAD_D, &fhisPO);
    NumlistToHistogramI(nsPO->svals, -1.0, BAD_D, BAD_D, &shisPO);
    ncol = GetHistogramNumBinsI(fhisPO);
    nrow = GetHistogramNumBinsI(shisPO);
    tabPO = CreateTablePO(nrow, ncol);
    if(!tabPO) {
        CHECK_HISTOGRAM(fhisPO); CHECK_HISTOGRAM(shisPO);
        return(FALSE);
    }
    /***
    *   Set lables
    */
    HistogramAutoFormatStringI(fhisPO,formS);
    for(c=0;c<ncol;c++) 
    {
        HistogramValuesForBinI(fhisPO, c, NULL, &fD, NULL);
        sprintf(numS,formS,fD);
        sprintf(labS,"C_%s",numS);
        ReplaceChars(' ',labS,'_',labS);
        SetTableColLabI(tabPO,c,labS);
    }
    HistogramAutoFormatStringI(shisPO,formS);
    for(r=0; r<nrow; r++)
    {
        HistogramValuesForBinI(shisPO, r, NULL, &fD, NULL);
        sprintf(numS,formS,fD);
        sprintf(labS,"Row_%s",numS);
        ReplaceChars(' ',labS,'_',labS);
        /* HAM? Why -2  *
        SetTableRowLabI(tabPO,nrow - r -2,labS);
        */
        SetTableRowLabI(tabPO,nrow - r -1,labS);
    }
    /***
    *   Fill table values
    */
    for(i=0; i<nsPO->num; i++)
    {
        GetNumlistDoubI(nsPO->vals,i,&fD);
        GetNumlistDoubI(nsPO->svals,i,&sD);
        HistogramBinForValueI(fhisPO, fD, &c, NULL,NULL);
        HistogramBinForValueI(shisPO, sD, &r, NULL,NULL);
        /* HAM? Why -2  *
        r = nrow - r - 2;
        */
        r = nrow - r - 1;
        GetTableValI(tabPO,r,c,&fD);
        fD += 1.0;
        SetTableValI(tabPO,r,c,fD);
    }
    fprintf(nsPO->out,"%s#\n",nsPO->echo);
    AutoTableOutFormattingI(tabPO, TRUE, FALSE);
    DumpTable(tabPO,FALSE,FALSE,nsPO->out);
    /***
    *   Clean up 
    */
    CHECK_TABLE(tabPO);
    CHECK_HISTOGRAM(fhisPO); CHECK_HISTOGRAM(shisPO);
    return(TRUE);
}
/*************************************************************************/
NUMSTAT *CreateNumstatPO()
{
    NUMSTAT *nsPO;

    if(!(nsPO=(NUMSTAT *)ALLOC(1,sizeof(NUMSTAT)))) {
        return(NULL);
    }
    nsPO->ID = NUMSTAT_ID;
    InitNumstat(nsPO);
    /***    
    *   Create empty numlist; will populate later
    */
    if(!(nsPO->vals = CreateNumlistPO(IS_DOUB, NULL, 0)) ) {
        CHECK_NUMSTAT(nsPO);
        return(NULL);
    }
    return(nsPO);
}
/*************************************************************************/
int DestroyNumstatI(NUMSTAT *nsPO)
{
    VALIDATE(nsPO,NUMSTAT_ID);
    CHECK_FILE(nsPO->in);
    CHECK_NFILE(nsPO->out,nsPO->outname);
    CHECK_NUMLIST(nsPO->vals);
    CHECK_NUMLIST(nsPO->svals);
    FREE(nsPO);
    return(TRUE);
}
/*************************************************************************/
int AddNumstatSecValsI(NUMSTAT *nsPO) 
{
    VALIDATE(nsPO,NUMSTAT_ID);
    if(!(nsPO->svals = CreateNumlistPO(IS_DOUB, NULL, 0)) ) {
        return(FALSE);
    }
    return(TRUE);
}
/*************************************************************************
*   Init structure
*/
void InitNumstat(NUMSTAT *nsPO)
{
    VALIDATE(nsPO,NUMSTAT_ID);
    INIT_S(nsPO->inname);
    nsPO->in = NULL;
    INIT_S(nsPO->outname);
    nsPO->out = NULL;
    nsPO->vals = nsPO->svals = NULL;
    nsPO->num = nsPO->lines = 0;
    nsPO->h_bin = nsPO->h_lo = nsPO->h_hi = BAD_D;
    nsPO->h_abin = nsPO->h_alo = nsPO->h_ahi = BAD_D;
    nsPO->h_maxbin = BOGUS;
    nsPO->h_mbin_set = FALSE;
    nsPO->h_pwide = HPLOT_SLEN;
    nsPO->do_hends = TRUE;
    nsPO->do_hist = FALSE;
    nsPO->do_hntb = FALSE;
    nsPO->do_dif = FALSE;
    nsPO->htb_xfold = HIS_BTR_XFOLD;
    strcpy(nsPO->h_pfmt,DEF_H_PTMF_S);
    strcpy(nsPO->hvsep,DEF_H_SEP_S);
    nsPO->do_perc = FALSE;  
    nsPO->col = nsPO->scol = BOGUS;
    nsPO->col_set = FALSE;
    INIT_S(nsPO->prls);
    nsPO->prsd = BAD_D;
    INIT_S(nsPO->echo);
    nsPO->do_efi = FALSE;
    nsPO->do_sk = FALSE;
    nsPO->do_sl = FALSE;
    nsPO->do_hplot = nsPO->do_ploti = FALSE;
    return;
}
/**************************************************************************
*   Handle loading values from column of data file
*   Allocates array space as needed
*/
int LoadColValsFromFileI(NUMSTAT *nsPO)
{
    int ok;
    char bufS[MAX_WIDTH+1];
    DOUB rD, sD;

    while(fgets(bufS,MAX_WIDTH,nsPO->in) != NULL)
    {
        nsPO->lines++;
        ok = GetLineDataValI(nsPO, bufS, &rD, &sD);
        if(IS_BOG(ok)) {
            PROBLINE;
            printf("Error reading data from line %d\n",nsPO->lines);
            fputs(bufS,stdout);
            return(FALSE);
        }
        if(ok) {
            ok = AppendNumlistDoubI(nsPO->vals, rD);
            if(nsPO->scol > 0) {
                ok = AppendNumlistDoubI(nsPO->svals, sD);
            }
            if(!ok) {
                PROBLINE;
                printf("Error saving data from line %d\n",nsPO->lines);
                fputs(bufS,stdout);
                return(FALSE);
            }
            nsPO->num++;
        }
    }
    /***
    *   Got anything?
    */
    if(nsPO->num < 1) {
        PROBLINE;
        printf("TOO FEW DATA FOR STATS\n");
        return(FALSE);
    }
    return(TRUE);
}
/*****************************************************************************/
int OpenNumstatFilesI(NUMSTAT *nsPO) 
{
    if(!(nsPO->in = OpenUFilePF(nsPO->inname,"r",NULL))) {      
        return(FALSE);  
    }
    if(!NO_S(nsPO->outname)) {
        if(!(nsPO->out = OpenUFilePF(nsPO->outname,"w",NULL))) {    
            return(FALSE); 
        }
    }
    HAND_NFILE(nsPO->out);
    return(TRUE);
}
/*****************************************************************************
*
*/
int CheckNumstatOptionsI(NUMSTAT *nsPO)
{
    if(!BAD_DOUB(nsPO->h_bin)) {    
        if(nsPO->h_bin < MIN_HIS_BIN) {
            printf("Histogram bin is too small; %f\n",nsPO->h_bin);
            return(FALSE);
        }
        nsPO->do_hist = TRUE; 
    }
    else if(nsPO->do_hplot || nsPO->do_ploti) { 
        nsPO->do_hist = TRUE;   
    }
    if( (!BAD_DOUB(nsPO->prsd)) || (!NO_S(nsPO->prls)) ) {  
        nsPO->do_perc = TRUE; 
    }
    if( nsPO->do_perc ) {
        nsPO->do_hist = nsPO->do_splot = FALSE;
    }
    else if(nsPO->do_hist || nsPO->do_splot) {
        nsPO->do_perc = FALSE;
    }
    /***
    *   Hist plot width
    */
    if( (nsPO->h_pwide < 1) || (nsPO->h_pwide > (DEF_BS-50)) )
    {
        printf("# Resetting histogram width %d to %d\n",
            nsPO->h_pwide,HPLOT_SLEN);
        nsPO->h_pwide = HPLOT_SLEN;
    }
    if(nsPO->h_maxbin > 0) {
        nsPO->h_mbin_set = TRUE;
    }
    else {
        nsPO->h_maxbin = DEF_HMB;
    }
    /***
    *   Possible dominant-bin truncation 
    */
    LIMIT_NUM(nsPO->htb_xfold, 1.0, TOO_BIG_D)
    /***
    *   Outputing lines with prefix? Get name; Add trailing tab
    */
    if(nsPO->do_efi) { 
        GetFilePartsI(nsPO->inname,NULL,nsPO->echo,NULL); 
    }
    if(!NO_S(nsPO->echo)) {
        strcat(nsPO->echo,"\t");
    }
    /***
    *   Input column set already?
    */
    if(nsPO->col > 0) {
        nsPO->col_set = TRUE;
    }
    /***
    *   Second data column
    */
    if(nsPO->scol > 0) {
        if(!AddNumstatSecValsI(nsPO)) {
            return(FALSE);
        }
    }
    /***
    *   Scatter plot
    */
    if( nsPO->do_splot && (!nsPO->svals) ) {
        printf("Scatter plot needs second column\n");
        return(FALSE);
    }
    /***
    *   Difference between cols
    */
    if( nsPO->do_dif && (!nsPO->svals) ) {
        printf("Difference needs second column\n");
        return(FALSE);
    }
    return(TRUE);
}
/*****************************************************************************/
int HandleHisMaxBinsI(NUMSTAT *nsPO)
{
    int n;

    n = GetNumlistLengthI(nsPO->vals);
    if( ((n/2) < nsPO->h_maxbin) && (!nsPO->h_mbin_set) ) {
        nsPO->h_maxbin = INT(n/2);
    }
    return(TRUE);
}
/*****************************************************************************
*   Get a value from data line
*/
int GetLineDataValI(NUMSTAT *nsPO, char *bufS, DOUB *rPD, DOUB *sPD)
{
    int j,ok,col;
    DOUB rD,sD;
    char *cPC;

    /***
    *   Comment? blank? Then ignore any set chars up front
    */
    if(COM_LINE(bufS)) {
        return(FALSE);
    }
    if(BlankStringI(bufS)) {
        return(FALSE);
    }
    cPC = bufS;
    /* Ignore chars up to? */
    for(j=0; j<nsPO->do_ic; j++)
    {
        if(!ISLINE(*cPC)) {
            break;
        }
        cPC++;
    }
    /***
    *   If data col not set, try to find first value on line
    */
    if(!nsPO->col_set) {
        col = 1;
        ok = GetWordDataValI(cPC, col, TRUE, NULL);
        while(ok == FALSE)
        {
            col++;
            ok = GetWordDataValI(cPC, col, FALSE, NULL);
        }
        if(ok == TRUE) {
            nsPO->col = col;
            nsPO->col_set = TRUE;
        }
        else {
            return(FALSE);
        }
    }
    /***
    *   Get first data col, maybe second too
    */
    ok = GetWordDataValI(cPC, nsPO->col, nsPO->do_sk, &rD);
    if(ok && sPD && (nsPO->scol >= 0)) {
        ok = GetWordDataValI(cPC, nsPO->scol, nsPO->do_sk, &sD);
    }
    if(ok) {
        *rPD = rD;
        if(sPD && (nsPO->scol >= 0)) {
            *sPD = sD;
        }
    }
    return(ok);
}
/**************************************************************************
*   Try to get value from word in col
*   skip makes this permissive of missing values, extra chars: ([,])
*/
int GetWordDataValI(char *linePC, int col, int skip, DOUB *vPD)
{
    DOUB vD;
    char wordS[DEF_BS],cwordS[DEF_BS];

    if(!GetNthWordI(linePC,col,wordS)) {
        if(skip) {
            return(FALSE);
        }
        return(BOGUS);
    }
    /***
    *   Possibly clean up extra chars in / around number
    */
    vD = BAD_D;
    if(skip) {
        RemoveChars("[(,)]", wordS, cwordS);
        sscanf(cwordS,"%lf",&vD);
    }
    else {
        sscanf(wordS,"%lf",&vD);
    }
    if(BAD_DOUB(vD)) {  
        if(skip) {
            return(FALSE);
        }
        return(BOGUS);
    }
    if(vPD) {
        *vPD = vD;
    }
    return(TRUE);
}
/**************************************************************************/
void NumstatOneLineOut(NUMSTAT *nsPO)
{
    char bufS[NSIZE];

    GetFilePartsI(nsPO->inname,NULL,bufS,NULL);
    fprintf(nsPO->out,"%s %d %4.4f %4.4f %4.4f %4.4f %4.4f\n", bufS, nsPO->num, 
    nsPO->min, nsPO->av - nsPO->sd, nsPO->av, nsPO->av + nsPO->sd, nsPO->max);
    return;
}
/**************************************************************************/
void NumstatHandleHeader(NUMSTAT *nsPO)
{
    char fS[DEF_BS], sS[DEF_BS];

    INIT_S(fS);
    if(nsPO->svals) {
        sprintf(fS,"Col_%d ",nsPO->col);
        sprintf(sS,"Col_%d ",nsPO->scol);
    }
    NumstatReportStats(nsPO, FALSE, fS, nsPO->out);
    /***
    *   If not doing diff, report second column story
    */
    if( (nsPO->svals) && (!nsPO->do_dif) ) {
        fprintf(nsPO->out,"%s#\n",nsPO->echo);
        NumstatReportStats(nsPO, TRUE, sS, nsPO->out);
    }
    return;
}
/*************************************************************************/
void NumstatReportStats(NUMSTAT *nsPO, int scol, char *exS, FILE *outPF)
{
    NUMLIST *valsPO;
    DOUB loD, hiD, avD, sdD, sumD, abinD, aloD, ahiD;
    char preS[DEF_BS];

    HAND_NFILE(outPF);
    if(scol) {
        valsPO = nsPO->svals;
    }   
    else {
        valsPO = nsPO->vals;
    }
    if(!valsPO) {
        return;
    }
    INIT_S(preS);
    if(exS && (!NO_S(exS))) {
        sprintf(preS,"%s",exS);
    }
    NumlistStatsI(valsPO, -1, -1, &loD, &hiD, &avD, &sdD);
    NumlistSumI(valsPO, -1, -1, &sumD);
    NumlistNaturalHistBinI(valsPO, nsPO->h_maxbin, TRUE, &abinD, &aloD, &ahiD);
    fprintf(nsPO->out,"%s#%s Number      %5d\n",nsPO->echo,preS, nsPO->num);
    fprintf(nsPO->out,"%s#%s Range       %10.4f to %8.4f = %8.4f\n",nsPO->echo,preS,
        loD, hiD, hiD-loD);
    fprintf(nsPO->out,"%s#%s Ave (SD)    %10.4f ( %4.4f )\n",nsPO->echo,preS,avD,sdD); 
    fprintf(nsPO->out,"%s#%s Sum       %12.4f\n", nsPO->echo,preS,sumD);
    fprintf(nsPO->out,"%s#%s BestBinning %10.4f from %10.4f to %10.4f\n", nsPO->echo,preS, 
        abinD, aloD, ahiD);
    return;
}
/**************************************************************************
*   Precentiles given as a comma-sep list
*/
void HandleListPercentiles(NUMSTAT *nsPO, DOUB *rvalsPD, int num, FILE *outPF)
{
    char *cPC;
    int p;
    DOUB pD,prsD;

    HAND_NFILE(outPF);
    /***
    *   Read list of comma-delimited percentiles
    */
    ReplaceChars(',',nsPO->prls,' ',nsPO->prls);
    fprintf(outPF,"%s# Percentiles: %s\n", nsPO->echo, nsPO->prls);
    if(nsPO->do_sl) {
        fprintf(outPF,"%sPercentiles",nsPO->echo);
    }
    cPC = nsPO->prls;
    PASS_BLANK(cPC);
    while(ISLINE(*cPC))
    {
        prsD = BAD_D; 
        sscanf(cPC,"%lf",&prsD);
        if(BAD_DOUB(prsD))
        {
            PROBLINE;
            printf("Bad percentile list format\n");
            printf("|%s|\n",nsPO->prls); 
            printf("Here: %s\n",cPC); 
            return;
        }
        if( (prsD<0.0) || (prsD>100.0) )
        {
            PROBLINE;
            printf("Bad percentile value listed\n");
            printf("|%s|\n",nsPO->prls); 
            printf("Bad: %f\n",prsD); 
            return;
        }
        pD = DNUM(nsPO->num) * prsD / 100.0;
        p = ROUND(pD);
        LIMIT_NUM(p, 0, (num - 1));
        ReportSingPercentile(nsPO, p, rvalsPD[p], num, outPF); 
        NEXT_WORD(cPC);
    }
    if(nsPO->do_sl) {
        fprintf(outPF,"\n");
    }
    return;
}
/**************************************************************************
*   Precentiles at given steps
*/
void HandleStepPercentiles(NUMSTAT *nsPO, DOUB *rvalsPD, int num, FILE *outPF)
{
    int p,last;
    DOUB pD,jD;

    if(BAD_DOUB(nsPO->prsd)) {
        return;
    }
    HAND_NFILE(outPF);
    /***
    *   Limit step
    */
    LIMIT_NUM(nsPO->prsd, 100.0/DNUM(num), 100.0);
    fprintf(outPF,"%s# Percentile steps of %4.2f (min to max)\n", nsPO->echo, 
        nsPO->prsd);
    if(nsPO->do_sl) {
        fprintf(outPF,"%sPercentiles",nsPO->prls);
    }
    /*** 
    *   Start at min == 0 percentile
    */ 
    jD = DNUM(num) * nsPO->prsd / 100.0;
    pD = 0.0;
    p = last = 0;
    while(p < num)
    {
        ReportSingPercentile(nsPO, p, rvalsPD[p], num, outPF); 
        pD += jD;
        last = p;
        p = ROUND(pD);
    }
    /*** 
    *   If not already reported, max == last percentile
    */ 
    if( (p != last) && (p <= num) ) {
        p = num - 1;
        ReportSingPercentile(nsPO, p, rvalsPD[p], num, outPF); 
    }
    if(nsPO->do_sl) {
        fprintf(outPF,"\n");
    }
    return;
}
/**************************************************************************
*   Precentiles given as a comma-sep list
*/
void ReportSingPercentile(NUMSTAT *nsPO, int p, DOUB vD, int num, FILE *outPF)
{
    HAND_NFILE(outPF);
    
    if(nsPO->do_sl) {
        fprintf(outPF, PERC_SL_FORM_S, vD);
    }
    else {
        fprintf(outPF,"%s",nsPO->echo);
        fprintf(outPF, PERCENTILE_FORM_S, ROUND(PERCENT_R(p,num)), 
                vD, PERCENT_R(p,num));
    }
    return;
}
/**************************************************************************
*   Replace values with differences between col 2 and 1
*/
int CalcColDifsI(NUMSTAT *nsPO)
{
    int i;
    DOUB fD, sD;

    /***
    *   Fill table values
    */
    for(i=0; i<nsPO->num; i++)
    {
        GetNumlistDoubI(nsPO->vals,i,&fD);
        GetNumlistDoubI(nsPO->svals,i,&sD);
        fD = sD - fD;
        SetNumlistDoubI(nsPO->vals,i,fD);
    }
    return(TRUE);
}
