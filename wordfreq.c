/*
* wordfreq.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <ctype.h>
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "table.h"
#include "wfutil.h"
#include "wordfreq.h"

#define DB_WF if(DB[99])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(WfUtilI(argc,argv),NULL) ); }
/**************************************************************************/
void WfUtilUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>   Sequence file\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -out XXX   Set output file to XXX\n");
    printf("   -size #    Set word size to # (def = %d)\n",DEF_WSIZE);
    printf("   -range # # Range of reported words # to # only\n");
    printf("   -norm      Normalize output (average word = 1.0)\n");
    printf("   -deg       Combine counts for +/- degenerate word (rev-comp)\n");
    printf("   -ilc       Ignore Lowercase case sequence\n");
    printf("   -iuc       Ignore Uppercase case sequence\n");
    printf("   -step #    Step between samples (default = 1)\n");
    printf("   -sran # #  Sequence range # to # only\n");
    printf("   -bran # #  Base range # to # only (i.e. sub-seqs)\n");
    printf("   -rre       Base range relative to end; i.e. Backwards\n");
    printf("   -pmat # #  Position-specific matrix from bases # to #\n");
    printf("   -sif XXX   Score input sequences via frequency tab XXX\n");
    printf("   -ssc       Score standard count values, e.g. first, second, etc\n");
    printf("   -scnl XXX  Score count list; comma delimited as \"1,2,5\"\n");
    printf("   -ssp       Score standard percentile values\n");
    printf("   -sprl XXX  Score percentile list; comma delimited as \"90,80,50\"\n");
    printf("   -ds        Dump (report) sequences appended as last column\n");
    printf("   -quiet     Suppress warnings about non-ACGT chars\n");
}
/**************************************************************************
*   Main program
*/
int WfUtilI(int argc, char **argv)
{
    int ok,tal,iraw,iseq,ifas,quiet,n;
    WF_UTIL *wfPO;

    wfPO = CreateWf_utilPO();
    iraw = iseq = ifas = quiet = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -iraw B -ifas B -siz I -deg B -bran I2 -rre B\
        -sran I2 -range I2 -norm B -sif S -pmat I2\
        -ilc B -iuc B -step I -quiet B -ds B -iseq B\
        -ssc B -ssp B -sprl S -scnl S",
        wfPO->inname, wfPO->outname, &iraw, &ifas, &wfPO->size, 
        &wfPO->do_deg, &wfPO->firstb,&wfPO->lastb, &wfPO->do_rre,
        &wfPO->firsts,&wfPO->lasts, &wfPO->min,&wfPO->max, &wfPO->do_norm,
        &wfPO->lisname, &wfPO->pmat_s,&wfPO->pmat_e,
        &wfPO->do_ilc, &wfPO->do_iuc, &wfPO->step, 
        &quiet, &wfPO->do_ds, &iseq,
        &wfPO->do_ssc, &wfPO->do_ssp, wfPO->sprl, wfPO->scnl,
        (int *)NULL))
    {
        WfUtilUse();
        CHECK_WF_UTIL(wfPO);
        return(FALSE);
    }
    /***
    *   Set input format 
    */
    wfPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas, wfPO->inname, TRUE);
    if(!wfPO->iform) {
        printf("Problem with input seq(s)\n");
        CHECK_WF_UTIL(wfPO);
        return(FALSE);
    }
    /***
    *   What output options?
    */
    if(!CheckWfuOptionsI(wfPO)) {
        ABORTLINE;
        CHECK_WF_UTIL(wfPO);
        return(FALSE);
    }
    if(!NO_S(wfPO->lisname)) {
        /* Have a sequence list to score, normalize */
        wfPO->do_norm = TRUE;
    }
    /***
    *   Open in/out files
    */
    if(!OpenWfuFilesI(wfPO)) {
        ABORTLINE;
        CHECK_WF_UTIL(wfPO);
        return(FALSE);
    }
    /***
    *   Set up auxillary stuff
    */
    if(!SetUpWfuAuxDataI(wfPO)) {
        ABORTLINE;
        CHECK_WF_UTIL(wfPO);
        return(FALSE);
    }
    SummaryHeader(wfPO, wfPO->out);
    /***
    *   Loop through the input sequence collection
    */
    wfPO->n = tal = n = 0;
    while(TRUE) {
        /***
        *   Parse sequence; FALSE = done
        */
        n++;
        ok = ParseSeqI(wfPO->in, wfPO->iform, n, SCLEAN_HI, !quiet, wfPO->fseq);
        if(ok==FALSE) {
            break;
        }
        if(ok!=TRUE) {
            continue;
        }
        /***
        *   Filter here BEFORE any manipulations
        *   Update count and qualify
        */
        if(!IsCurrentSeqOkI(wfPO, wfPO->fseq, wfPO->n)) {
            wfPO->n += 1;
            continue;
        }
        wfPO->n += 1;
        /* xxx Unconditional copy unneeded / slow??? */
        CopySeqI(wfPO->fseq, wfPO->seq, -1, -1);
        /***
        *   Trim, case-based-mask?
        */
        if(!HandleWfuSubseqI(wfPO, wfPO->seq, wfPO->fseq)) {
            ABORTLINE;
            break;
        }
        if( (wfPO->do_ilc) || (wfPO->do_iuc) ) {
            HandleWfuCaseMaskI(wfPO, wfPO->seq);
        }
        /***
        *   Processing current seq, or Tallying?
        */
        if(!NO_S(wfPO->lisname)) {
            HandleWfuSeqOutputI(wfPO, wfPO->seq, wfPO->fseq, wfPO->out);
        }
        else {
            tal = HandleWfTallyI(wfPO, wfPO->seq);
            if(tal<1) {
                break;
            }
            wfPO->tw += tal;
        }
    }
    /***
    *   Matrix or Stats?
    */
    if(wfPO->pmat_s > 0) {
        if(tal>0) {
            HandleWfuPosMat(wfPO, wfPO->out);
        }
    }
    else {
        HandleWfuStats(wfPO, wfPO->out);
    }
    /***
    *   All done
    */
    CHECK_WF_UTIL(wfPO);
    return(TRUE);
}
/*****************************************************************************
*   Create data struct
*/
WF_UTIL *CreateWf_utilPO()
{
    WF_UTIL *wfPO;

    if(! (wfPO = (WF_UTIL *)ALLOC(1,sizeof(WF_UTIL)) ) )
    {
        printf("# Failed to allocate working object\n");
        return(NULL);
    }
    wfPO->ID = WF_UTIL_ID;
    wfPO->seq = CreateSeqPO(0,NULL,NULL);
    wfPO->fseq = CreateSeqPO(0,NULL,NULL);
    InitWf_util(wfPO);
    return(wfPO);
}
/*****************************************************************************
*   Free datastructure and substructs
*/
int DestroyWf_utilI(WF_UTIL *wfPO)
{
    VALIDATE(wfPO,WF_UTIL_ID);
    CHECK_FILE(wfPO->in);
    CHECK_NFILE(wfPO->out,wfPO->outname);
    CHECK_SEQ(wfPO->seq);
    CHECK_SEQ(wfPO->fseq);
    CHECK_WORDFREQ(wfPO->wf);
    CHECK_TABLE(wfPO->pmat);
    FREE(wfPO);
    return(TRUE);
}
/*****************************************************************************
*   Set null / default values
*/
void InitWf_util(WF_UTIL *wfPO)
{
    VALIDATE(wfPO,WF_UTIL_ID);

    INIT_S(wfPO->inname);
    INIT_S(wfPO->outname);
    wfPO->in = NULL;
    wfPO->iform = BOGUS;
    wfPO->out = NULL;
    wfPO->owhat = BOGUS;
    wfPO->do_not = FALSE;
    INIT_S(wfPO->lisname);
    wfPO->do_ssc = wfPO->do_ssp = FALSE;
    INIT_S(wfPO->sprl);
    INIT_S(wfPO->scnl);
    wfPO->do_ds = FALSE;
    wfPO->wf = NULL;
    wfPO->size = DEF_WSIZE;
    wfPO->step = 1;
    wfPO->do_deg = FALSE;
    wfPO->do_rre = FALSE;
    wfPO->firsts = wfPO->firstb = -1; 
    wfPO->lasts = wfPO->lastb = TOO_BIG;
    wfPO->min = 0; 
    wfPO->max = TOO_BIG;
    wfPO->n = 0;
    wfPO->tw = 0;
    wfPO->pmat_s = wfPO->pmat_e = BOGUS;
    wfPO->do_ilc = FALSE;
    wfPO->do_iuc = FALSE;
}
/*************************************************************************
*   Check for option consistency
*/
int CheckWfuOptionsI(WF_UTIL *wfPO)
{
    /***
    *   Not scoring seqs, then size has to be OK
    */
    if(!NO_S(wfPO->lisname)) {
        if(wfPO->size<1) {
            PROBLINE;
            printf("Bogus size: %d\n",wfPO->size);
            return(FALSE);
        }
    }
    /***
    *   Position-specific matrix checks
    */
    if(wfPO->pmat_s >= 0) {
        if(wfPO->size > MAX_PSM_SIZE) {
            PROBLINE;
            printf("Max word size for Position-specific matrix is %d\n",
                MAX_PSM_SIZE);
            printf("  Sorry, %d is too big\n",wfPO->size);
            return(FALSE);
        }
        if( wfPO->pmat_s > wfPO->pmat_e ) {
            PROBLINE;
            printf("Bad range for Position-specific matrix:\n");
            printf("  Start=%d  End=%d\n", wfPO->pmat_s,wfPO->pmat_e);
            return(FALSE);
        }
    }
    /***
    *   Case ignorning options
    */
    if( (wfPO->do_ilc) && (wfPO->do_iuc) ) {
        PROBLINE;
        printf("Can't ignore both upper and lower case!\n"); 
        return(FALSE);
    }
    /***
    *   Step size has to be >=1 
    */
    if(wfPO->step < 1) {
        PROBLINE;
        printf("Step size has to be >=1\n"); 
        return(FALSE);
    }
    return(TRUE);
}
/***************************************************************************
*   Open any needed files or die
*/
int OpenWfuFilesI(WF_UTIL *wfPO)
{
    if(!(wfPO->in=OpenUFilePF(wfPO->inname,"r",NULL))) {
        return(FALSE);
    }
    if(!NO_S(wfPO->outname)) {
        if(!(wfPO->out=OpenUFilePF(wfPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    HAND_NFILE(wfPO->out);
    return(TRUE);
}
/***************************************************************************
*   Set up auxillary data structs / settings
*/
int SetUpWfuAuxDataI(WF_UTIL *wfPO)
{
    if(!NO_S(wfPO->lisname)) {
        if(!GetWordFreqsI(wfPO->lisname,&wfPO->wf)) {
            PROBLINE;
            printf("Failed to load word freq table\n");
            return(FALSE);
        }
        if(wfPO->do_norm) {
            NormalizeFrecs(wfPO->wf);
        }
    }
    else {
        wfPO->wf = CreateWordfreqPO(wfPO->size,ALPHDIM); 
        if(!wfPO->wf) {
            PROBLINE;
            printf("Failed to allocate space for wordsize %d\n",wfPO->size);
            return(FALSE);
        }
    }
    /***
    *   Position specific matrix; rows = N-mer, cols = base position
    */
    if(wfPO->pmat_s > 0) {
        wfPO->pmat_r = CalcInDimI(wfPO->size,ALPHDIM);
        wfPO->pmat_c = wfPO->pmat_e - wfPO->pmat_s + 1;
        if(! (wfPO->pmat=CreateTablePO(wfPO->pmat_r,wfPO->pmat_c)) ) {
            PROBLINE;
            printf("Failed to allocate postition-specific matrix\n");
            return(FALSE);
        }
        SetPosMatLablesI(wfPO);
    }
    /***
    *   Score percent or count list?
    */
    if(wfPO->do_ssc) {
        strcpy(wfPO->scnl, SSC_LIST_S);
    }
    if(wfPO->do_ssp) {
        strcpy(wfPO->sprl, SSP_LIST_S);
    }
    /***
    *   Replace any commas in list strings with space (tabs for header print)
    */
    if(!NO_S(wfPO->scnl)) {
        ReplaceChars(',', wfPO->scnl, '\t', wfPO->scnl);
    }
    if(!NO_S(wfPO->sprl)) {
        ReplaceChars(',', wfPO->sprl, '\t', wfPO->sprl);
    }
    return(TRUE);
}
/***************************************************************************
*   Set row and column lables for position matrix
*/
int SetPosMatLablesI(WF_UTIL *wfPO)
{
    int i;
    char nameS[NSIZE];
    
    /***
    *   Cols = base position relative to 3' or 5' end
    */
    for(i=0; i<wfPO->pmat_c; i++) {
        sprintf(nameS,"%d",wfPO->pmat_s + i);
        SetTableColLabI(wfPO->pmat,i,nameS);
    }
    /***
    *   Rows = N-mer
    */
    for(i=0; i<wfPO->pmat_r; i++) {
        FillWordfreqSeqStringI(wfPO->wf,i,nameS);
        SetTableRowLabI(wfPO->pmat,i,nameS);
    }
    /***
    *   Print format; <print> <sep> <preface>
    *   If not normalizing, simple counts are good
    */
    if(wfPO->do_norm) {
        SetTablePrintformI(wfPO->pmat,"%6.4f",NULL,"\t",NULL,NULL);
    }
    else {
        SetTablePrintformI(wfPO->pmat,"%6.0f",NULL,"\t",NULL,NULL);
    }
    return(TRUE);
}
/****************************************************************************
*   If limited base range, copy full seq fseqPO and modify (bound) seqPO
*   Also set fseqPO case to reflect any bounding 
*/
int HandleWfuSubseqI(WF_UTIL *wfPO, SEQ *seqPO, SEQ *fseqPO)
{
    int len, flen;

    if(wfPO->firstb > 0) {
        len = wfPO->lastb - wfPO->firstb + 1;
        SetCaseSeqSubseqI(fseqPO, FALSE, -1, -1);
        if(wfPO->do_rre) {
            NarrowSeqI(seqPO, wfPO->firstb-1, len, REVERSE, FALSE);
            flen = GetSeqLenI(fseqPO);
            SetCaseSeqSubseqI(fseqPO, TRUE, flen - wfPO->lastb, flen - wfPO->firstb + 1);
        }
        else {
            NarrowSeqI(seqPO, wfPO->firstb-1, len, FORWARD, FALSE);
            SetCaseSeqSubseqI(fseqPO, TRUE, wfPO->firstb-1, wfPO->lastb);
        }
    }
    return(TRUE);
}
/***************************************************************************
*   If masking lowercase, switch these to n and return the count 
*/
int HandleWfuCaseMaskI(WF_UTIL *wfPO, SEQ *seqPO)
{
    int i,n,len;
    char *seqPC;

    len = GetSeqLenI(seqPO);
    GetSeqSeqI(seqPO, &seqPC);
    n = 0;
    for(i=0;i<len;i++) {
        if( (wfPO->do_ilc) && (islower(INT(seqPC[i]))) ) {
            seqPC[i] = 'n';
            n++;
        }
        if( (wfPO->do_iuc) && (isupper(INT(seqPC[i]))) ) {
            seqPC[i] = 'n';
            n++;
        }
    }
    return(n);
}
/**************************************************************************
*   Screen current seq against filters
*/
int IsCurrentSeqOkI(WF_UTIL *wfPO, SEQ *seqPO, int n)
{
    int ok;

    ok = TRUE;
    if( (n < wfPO->firsts) || (n > wfPO->lasts) ) {
        ok = FALSE;
    }
    /***
    *   If not, invert qualification
    */
    if(wfPO->do_not) {
        ok = !ok;
    }
    return(ok);
}
/***************************************************************************/
void HandleWfuStats(WF_UTIL *wfPO,FILE *outPF)
{
    HAND_NFILE(outPF);
    if(wfPO->n < 1) {
        fprintf(outPF,"NO SEQS!\n");
    }
    if(NO_S(wfPO->lisname)) {
        if(wfPO->do_norm) {
            NormalizeFrecs(wfPO->wf);
        }
        DumpWordsI(wfPO->wf, wfPO->do_deg, wfPO->min, wfPO->max, outPF);
    }
}
/***************************************************************************
*   Score current seq against (already loaded) wordfreqs
*/
int HandleWfuSeqOutputI(WF_UTIL *wfPO, SEQ *seqPO, SEQ *fseqPO, FILE *outPF)
{
    int n,seqlen;
    char nameS[NSIZE], *seqPC;
    NUMLIST *numsPO;
    DOUB avD,minD,maxD;

    HAND_NFILE(outPF);
    GetSeqSeqI(seqPO,&seqPC);
    seqlen = GetSeqLenI(seqPO);
    /* Output name regardless */
    FillSeqNameStringI(seqPO, nameS, NSIZE);
    fprintf(outPF,"%s\t", nameS);
    /* List to collect counts */
    numsPO = CreateNumlistPO(IS_DOUB, NULL, seqlen);
    n = CollectWfuSeqValsI(wfPO, seqPO, numsPO);
    if(n<1) {
        CHECK_NUMLIST(numsPO);
        fprintf(outPF,"No words to count!\n");
        return(FALSE);
    }
    /***
    *   Outputing what numbers?
    */
    if(!NO_S(wfPO->scnl)) {
        if(!DumpListedNumlistValsI(numsPO, wfPO->scnl, FALSE, outPF)) {
            fprintf(outPF,"Problem score counts |%s|\n",wfPO->scnl);
        }
    }
    else if(!NO_S(wfPO->sprl)) {
        if(!DumpListedNumlistValsI(numsPO, wfPO->sprl, TRUE, outPF)) {
            fprintf(outPF,"Problem score percentiles |%s|\n",wfPO->sprl);
        }
    }
    else {
        NumlistStatsI(numsPO, -1,-1, &minD, &maxD, &avD, NULL);
        fprintf(outPF,"%5.4f\t%5.4f\t%5.4f",minD,avD,maxD);
    }
    /* Dump seq too? */
    if(wfPO->do_ds) {
        GetSeqSeqI(fseqPO, &seqPC);
        fprintf(outPF,"\t%s", seqPC);
    }
    fprintf(outPF,"\n");
/*
DumpNumlist(numsPO, -1, -1, "[%d]\t%5.4f\n", NULL);
*/
    CHECK_NUMLIST(numsPO);
    return(TRUE);
}
/***************************************************************************
*   Dump numlist values based on numbers in list string. 
*   If do_perc, treat as percentiles, else count (rank index)
*/
int DumpListedNumlistValsI(NUMLIST *valsPO, char *lisS, int do_per, FILE *outPF)
{
    int i,n,len;
    DOUB vD;
    char *cPC;
    
    HAND_NFILE(outPF);
    len = GetNumlistLengthI(valsPO);
    cPC = lisS;
    PASS_BLANK(cPC);
    n = 0;
    while(ISLINE(*cPC))
    {
        /* Parse number from string */
        vD = -1.0;
        sscanf(cPC,"%lf",&vD);
        if(vD < 0.0)
        {
            PROBLINE;
            printf("Bad number list |%s|\n", lisS); 
            printf("Here: %s\n",cPC); 
            return(FALSE);
        }
        /***
        *   Count index or percentile index; Percentile backwards
        */
        if(do_per) {
            i = len - INT(vD * DNUM(len) / 100.0);
        }
        else {
            i = INT(vD) - 1;
        }
        LIMIT_NUM(i, 0, len-1);
        /* Get value ... or don't */
        GetNumlistDoubI(valsPO, i, &vD); 
        if(n>0) {
            fprintf(outPF,"\t");
        }
        fprintf(outPF,"%5.4f",vD);
        n++;
        NEXT_WORD(cPC);
    }
    return(n);
}
/***************************************************************************
*   Set word freq count values for seq into numlist
*/
int CollectWfuSeqValsI(WF_UTIL *wfPO, SEQ *seqPO, NUMLIST *valsPO)
{
    int i,n,ind,seqlen;
    DOUB valD;
    char *seqPC;
    WORDFREQ *wordfPO;

    wordfPO = wfPO->wf;
    GetSeqSeqI(seqPO,&seqPC);
    seqlen = GetSeqLenI(seqPO);
    n = 0;
    for(i=0; i<=(seqlen - wordfPO->size); i += wfPO->step)
    {
        ind = IndexFromSeqI(&seqPC[i], wordfPO->size, wordfPO->n, wordfPO->ald);
        if( (ind >= wordfPO->n) || (ind<0) ) {
            printf("Bogus index: %d i=%d seq|%s|\n",ind,i,seqPC);
            ERR("CollectWfuSeqValsI","bad index");
            return(FALSE);
        }
        valD = wordfPO->freqs[ind].n;
        /*** 
        * Degenerate means need to add in compliment too
        */
        if(wfPO->do_deg) {
            ind = CompIndexI(ind, wordfPO->pmax, wordfPO->ald);
            if( (ind >= wordfPO->n) || (ind<0) ) {
                printf("Bogus comp index: %d i=%d seq|%s|\n",ind,i,seqPC);
                ERR("CollectWfuSeqValsI","bad index");
                return(FALSE);
            }
            valD += wordfPO->freqs[ind].n;
        }
        AddNumlistDoubI(valsPO, i, valD);
        n++;
    }
    /***
    *   Set length populated and sort decending
    */
    SetNumlistLengthI(valsPO, n);
    SortNumlistI(valsPO, -1);
    return(n);
}
/***************************************************************************
*   Report position-specific matrix
*/
void HandleWfuPosMat(WF_UTIL *wfPO, FILE *outPF)
{
    int r,c;
    DOUB vD;

    HAND_NFILE(outPF);
    fprintf(outPF,"# Position-specific word frequencies\n");
    fprintf(outPF,"# Words of size %d\n",wfPO->size);
    fprintf(outPF,"# Bases %d to %d",wfPO->pmat_s,wfPO->pmat_e);
    if(wfPO->do_rre) {   
        fprintf(outPF," from 3' end\n"); 
    }
    else {   
        fprintf(outPF," from 5' end\n"); 
    }
    /***
    *   Scale shams?
    */
    if(wfPO->do_norm) {
        for(r=0;r<wfPO->pmat_r;r++) {
            for(c=0;c<wfPO->pmat_c;c++) {
                GetTableValI(wfPO->pmat,r,c,&vD);
                vD = vD / DNUM(wfPO->n);
                SetTableValI(wfPO->pmat,r,c,vD);
            }
        }
    }   
    /***
    *   Dump table
    */
    DumpTable(wfPO->pmat,FALSE,FALSE,outPF);
}
/***************************************************************************
*   Tally words for current sequence
*/
int HandleWfTallyI(WF_UTIL *wfPO, SEQ *seqPO)
{
    int n;

    if(GetSeqLenI(seqPO) < wfPO->size) {
        return(0);
    }
    /***
    *   Position matrix or just word talley
    */
    if(wfPO->pmat_s > 0) {
        n = TallyPosMatWordsI(wfPO, seqPO);
    }
    else {
        n = TallyWordsI(seqPO, wfPO->wf, wfPO->step);
    }
    return(n);
}
/***************************************************************************
*   Update postition-specific matrix table with current seq
*   If problem, return FALSE
*/
int TallyPosMatWordsI(WF_UTIL *wfPO, SEQ *seqPO)
{
    int i,s,ind,len;
    char *seqPC;
    WORDFREQ *wordPO;

    wordPO = wfPO->wf;
    len = GetSeqLenI(seqPO);
    GetSeqSeqI(seqPO, &seqPC);
    /***
    *   For each position (col of matrix)
    */
    for(i=0; i < wfPO->pmat_c; i++) {
        /***
        *   Where to start sampling from?
        *   Given start coord is 1-based, so subtract 1
        */
        if(wfPO->do_rre) {
            s = len - (wfPO->pmat_s -1) - i - wfPO->size;
        }
        else {
            s = (wfPO->pmat_s -1) + i;
        }
        /***
        *   Out of sampling room for seq len?
        */
        if( (s<0) || (s>(len-wfPO->size)) ) {
            break;
        }
        /***
        *   Get word index (row of matrix) and update table
        */
        ind = IndexFromSeqI(&seqPC[s], wordPO->size, wordPO->n, wordPO->ald);
        if(!ModTableValI(wfPO->pmat,ind,i,1.0,MATH_ADD)) {
            printf("Row=ind=%d col=i=%d\n",ind,i);
            ERR("TallyPosMatWordsI","ModTableValI failed");
            return(FALSE);
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Write header info for wordfreq data
*/
void SummaryHeader(WF_UTIL *wfPO, FILE *outPF)
{
    HAND_NFILE(outPF);
    VersionSplash(outPF,VERSION_S,"#  ",TRUE);
    fprintf(outPF,"# Input        %s\n",wfPO->inname);
    if(wfPO->firstb > 0) {
        if(wfPO->do_rre) {
            fprintf(outPF,"#   Last %d to %d bases only\n",wfPO->firstb,wfPO->lastb);
        }
        else {
            fprintf(outPF,"#   First %d to %d bases only\n",wfPO->firstb,wfPO->lastb);
        }
    }
    if(!NO_S(wfPO->lisname)) {
        fprintf(outPF,"# Scoring sequences by normalized frequencies\n");
        fprintf(outPF,"# Word freq table: %s\n", wfPO->lisname);
        if(!NO_S(wfPO->scnl)) {
            fprintf(outPF,"# Score count (rank index) values\n");
            fprintf(outPF,"# <name>\t%s\n", wfPO->scnl);
        }
        else if(!NO_S(wfPO->sprl)) {
            fprintf(outPF,"# Score percentile values\n");
            fprintf(outPF,"# <name>\t%s\n", wfPO->sprl);
        }
        else {
            fprintf(outPF,"# Score simple stats\n");
            fprintf(outPF,"# <name>\t<min>\t<mean>\t<max>\n");
        }
    }
    else {
        WordfreqSummary(wfPO->wf, wfPO->min,wfPO->max, outPF);
        fprintf(outPF,"#\n");
    }
}
