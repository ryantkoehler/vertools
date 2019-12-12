/*
* tm_util.c
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
#include <math.h>
#include <ctype.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "table.h"
#include "fbound.h"
#include "tm_pars.h"
#include "dset_tm.h"
#include "competfb.h"
#include "tm_util.h"

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(TmUtilI(argc,argv), NULL) ); }
/**************************************************************************/
void TmUtilUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>   Sequences file to read in\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -out XXX   Set output to XXX\n");
    printf("   -par XXX   Use parameters in file XXX for Tm\n");
    printf("   -pdc #     Parameter \"default\" column set to # (1-base)\n");
    printf("   -con #     Set DNA concentration to # (M)\n");
    printf("   -scon #    Set second strand DNA conc. to # (M)\n");
    printf("   -tcon #    Set third strand DNA conc. to # (coaxial; M)\n");
    printf("   -sal #     Set Salt concentration to # (M)\n");
    printf("   -mg #      Set Mg concentration to # (M)\n");
    printf("   -tem #     Set temperature (for dG) to # C'\n");
    printf("   -tmsl      Use SantaLucia's algorithm\n");
    printf("   -tmpey     Use Peyret's (SantaLucia) algorithm\n");
    printf("   -tmoli     Use the algorithm as in \"Oligo\" program\n");
    printf("   -tmelt     Use the algorithm as in \"Melting\" program\n");
    printf("   -tmpna     Use the PNA algorithm\n");
    printf("   -tm24      Use the algorithm Tm = 2AT + 4GC\n");
    printf("   -bran # #  Base range cut to # to # (1-base coords)\n");
    printf("   -rre       Base range relative to end; i.e. backwards\n");
    printf("   -dpar      Dump (report) Tm parameters / settings used\n");
    printf("   -ds        Dump (report) sequences appended as last column\n");
    printf("   -the       Report thermodynamic quantities dG dH dS\n");
    printf("   -fds       Report fraction double strand (i.e. Fraction Bound, FB)\n");
    printf("   -tlen #    Report sequence length required for Tm (or FB) of #\n");
    printf("   -otls      Output sequences for -tlen option\n");
    printf("   -tlmin     Output seqs at least target value for -otls option\n");
    printf("   -tlp #     Set position for -tlen option to base seq[#]\n");
    printf("   -tlrev     Evaluate -tlen (from -tlp) IN REVERSE\n");
    printf("   -tmpro     Report Tm profile as function of length\n");
    printf("   -tab # #   Report table for lengths # to #\n");
    printf("   -flg # #   Flag sequences with Tm's in the range # to #\n");
    printf("   -tbj #     Table option jump between starting bases (i.e. step)\n");
    printf("   -not       Flagging logical not; out>-->in & in>-->out\n");
    printf("   -eraw      Extract flagged sequences in raw format (with extra info / line)\n");
    printf("   -emin      Extract min-length sequences (with -tab and -flg)\n");
    printf("   -emid      Extract mid-length sequences (with -tab and -flg)\n");
    printf("   -emax      Extract max-length sequences (with -tab and -flg)\n");
    printf("   -tes       Two explicit seqs (alternating lines), 5'-3' and 3'-5'\n");
    printf("   -btes      Two explicit seqs (e.g. Blast output), both 5'-3'\n");
    printf("   -dthe      Report Delta thermodynamic quantities (-tes -btes)\n");
/*
    printf("   -den       Dangling end seqs (Peyret algorithm)\n");
*/
    printf("   -cmb       Competitive miss-match binding (Peyret algorithm)\n");
    printf("\n");
    printf("NOTE: To extract seqs of lengths L1 to L2 with Tm vals T1 to T2:\n");
    printf("   -tab L1 L2 -flg T1 T2 -eraw (or -emin, -emid, -emax)\n");
    printf("\n");
}
/**************************************************************************
*   Main program
*/
int TmUtilI(int argc, char **argv)
{
    int tmoli,tmpna,tmelt,tm24,tmgb,tmsl,tmpey,slen,good,ngood,n;
    int iraw,ifas,iseq;
    TM_UTIL *tuPO;
    char *seqPC,nameS[NSIZE];

    tuPO = CreateTm_utilPO();
    tmoli = tmpna = tmelt = tm24 = tmgb = tmsl = tmpey = FALSE;
    iraw = ifas = iseq = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -par S -con D\
        -sal D -tmoli B -tmelt B -the B -dpar B\
        -tmpro B -tm24 B -tmgb B -mg D -quiet B\
        -flg D2 -not B -tlen D -tlp I\
        -tlrev B -tmpna B\
        -tmsl B -tmpey B -scon D -rsd B -tem D -fds B\
        -tlmin B -tes B -btes B -den B\
        -tcon D -eraw B -bran I2 -rre B -otls B\
        -iraw B -ifas B -dthe B -cmb B -pdc I -tab I2\
        -emin B -emid B -emax B -tbj I -ds B -iseq B -rc B",
        tuPO->inname, tuPO->outname, tuPO->parname, &tuPO->con1, 
        &tuPO->salt, &tmoli, &tmelt, &tuPO->do_therm, &tuPO->do_dpar,
        &tuPO->do_tmpro, &tm24, &tmgb, &tuPO->mg, &tuPO->quiet, 
        &tuPO->minv,&tuPO->maxv, &tuPO->do_not, &tuPO->tlen, &tuPO->tlp, 
        &tuPO->do_tlrev, &tmpna,
        &tmsl, &tmpey, &tuPO->con2, &tuPO->do_rsd, &tuPO->temp, &tuPO->do_fds, 
        &tuPO->do_tlmin, &tuPO->do_tes, &tuPO->do_btes,  
        &tuPO->do_den, &tuPO->con3, &tuPO->do_eraw,
        &tuPO->firstb,&tuPO->lastb, &tuPO->do_rre, &tuPO->do_otls,
        &iraw, &ifas,
        &tuPO->do_dtherm, &tuPO->do_cmb, 
        &tuPO->pdc,
        &tuPO->do_tab_lo,&tuPO->do_tab_hi,
        &tuPO->do_emin, &tuPO->do_emid, &tuPO->do_emax, 
        &tuPO->do_tab_j, &tuPO->do_ds, 
        &iseq, &tuPO->do_rc,
        (int *)NULL))
    {
        TmUtilUse();
        CHECK_TM_UTIL(tuPO);
        return(FALSE);
    }
    /***
    *   Open files
    */
    if(!OpenTmutilFilesI(tuPO)) {
        ABORTLINE;
        CHECK_TM_UTIL(tuPO);
        return(FALSE);
    }
    /***
    *   Set input format 
    */
    tuPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas, tuPO->inname, TRUE);
    if( ! tuPO->iform ) {
        printf("Problem with input seq(s)\n");
        CHECK_TM_UTIL(tuPO);
        return(FALSE);
    }
    /***
    *   Set algorithm then make sure options are consistent / workable
    */
    if(tm24) {
        tuPO->algo = ALGO_24;
    }
    else if(tmoli) {
        tuPO->algo = ALGO_OLIGO;
    }
    else if(tmelt) {
        tuPO->algo = ALGO_MELTING;
    }
    else if(tmpna) {
        tuPO->algo = ALGO_PNA;  
    }
    else if(tmpey) {
        tuPO->algo = ALGO_PEYRET;   
    }
    else if(tmsl) {
        tuPO->algo = ALGO_SANTA;    
    }
    else {
        tuPO->algo = DEF_TM_ALGO;   
    }
    if(!CheckTmutilOptionsI(tuPO))
    {
        ABORTLINE;
        CHECK_TM_UTIL(tuPO);
        return(FALSE);
    }
    /***
    *   Report settings / what will be output
    */
    ReportTmutilSettings(tuPO,tuPO->out);
    ReportTmutilOutput(tuPO,tuPO->out);
    /***
    *   Loop over input until no more
    */
    ngood = n = 0; 
    while(TRUE)
    {
        n++;
        slen = LoadTmutilSequenceI(tuPO, n);
        if(slen==BOGUS) {
            break;
        }
        /***
        *   Get name and check that can get a Tm
        *   If doing table / extract, don't cap len
        */
        FillSeqNameStringI(tuPO->seq,nameS,NSIZE);
        if( (!tuPO->do_tab_lo) &&
            ((slen<MIN_TM_LEN) || (slen>=MAX_TM_LEN)) ) {
            fprintf(tuPO->out,"# %s is too short/long (%d; %d to %d OK)\n",
                nameS,slen, MIN_TM_LEN, MAX_TM_LEN);
            continue;
        }
        if(GetSeqSnpCountI(tuPO->seq) > 0) {
            fprintf(tuPO->out,"# %s has SNP(s)\n",nameS);   
            continue;
        }
        if(AnySeqAmbigsI(tuPO->seq)) {
            fprintf(tuPO->out,"# %s has ambiguous bases\n",nameS);
            continue;
        }
        /***
        *   Thermo table or profile?
        */
        if( tuPO->do_tab_lo ) {
            HandleThermoSampleTableI(tuPO,tuPO->out);
            continue;
        }
        if(tuPO->do_tmpro) {
            fprintf(tuPO->out,"# %s\n",nameS);  
            HandleTmProfileI(tuPO,tuPO->out);
            continue;
        }
        /***
        *   Looking for target Tm lengths?
        */
        if(!BAD_REAL(tuPO->tlen)) {
            HandleTmLenI(tuPO,tuPO->out);
            continue;
        }
        GetSeqSeqI(tuPO->seq,&seqPC);
        slen = GetSeqLenI(tuPO->seq);
        /***
        *   Nic Peyret's two-explicit-strand case?
        *   CallDsetI sets calculated numbers into passed TM_UTIL struct
        */
        if( (tuPO->algo==ALGO_PEYRET) && ((tuPO->do_tes) || (tuPO->do_btes)) ) {
            if(!CallDsetTwoSeqsI(tuPO)) {
                ReportDsetProblemSeq(nameS,tuPO->out);
                if(tuPO->do_padbad) {
                    tuPO->Tm = tuPO->fds = tuPO->dG = tuPO->dH = tuPO->dS = -TOO_BIG_R;
                }
            }
        }
        else {
            SeqTmThermI(tuPO->tm, seqPC, slen, &tuPO->Tm, &tuPO->dG,
                &tuPO->dH, &tuPO->dS);
            /***
            *   Fraction double stranded == fraction bound
            *   SHAM; sometimes fds is out of range 0 to 100
            */
            if(tuPO->do_fds) {
                tuPO->fds = Fraction2D(tuPO->con1, tuPO->con2, tuPO->temp, tuPO->dG);
            }
            LIMIT_NUM(tuPO->fds, 0.0, 100.0);
            LIMIT_NUM(tuPO->fds2, 0.0, 100.0);
        }
        /***
        *   Qualify via any filters
        */
        good = SeqFlaggingValI(tuPO,NULL);
        ngood += good;
        HandleTmutilOutputI(tuPO,good,seqPC,slen,nameS,tuPO->out);
    }
    /***
    *   All done
    */
    CHECK_TM_UTIL(tuPO);
    return(TRUE);
}
/*****************************************************************************
*   Create data structure
*/
TM_UTIL *CreateTm_utilPO()
{
    TM_UTIL *tuPO;

    if(! (tuPO=(TM_UTIL *)ALLOC(1,sizeof(TM_UTIL))) ) {
        return(NULL);
    }
    tuPO->ID = TM_UTIL_ID;
    tuPO->tm = CreateTm_parsPO();
    tuPO->seq = CreateSeqPO(MAX_TM_LEN, NULL, NULL);
    tuPO->sseq = CreateSeqPO(MAX_TM_LEN, NULL, NULL);
    tuPO->fseq = CreateSeqPO(MAX_TM_LEN, NULL, NULL);
    tuPO->ltm_tab = CreateTablePO(DEF_LTM_ROWS, DEF_LTM_COLS);
    if( (!tuPO->tm) || (!tuPO->seq) || (!tuPO->sseq) || (!tuPO->fseq) || (!tuPO->ltm_tab) ) {
        CHECK_TM_UTIL(tuPO);
        return(NULL);
    }
    SetTablePrintformI(tuPO->ltm_tab, TM_UTIL_PFMT, NULL, "\t", NULL, NULL); 
    InitTm_utilI(tuPO);
    return(tuPO);
}
/*****************************************************************************
*   Clean up data structure and sub fields
*/
int DestroyTm_utilI(TM_UTIL *tuPO)
{
    VALIDATE(tuPO,TM_UTIL_ID);
    CHECK_TM_PARS(tuPO->tm);
    CHECK_SEQ(tuPO->seq);
    CHECK_SEQ(tuPO->sseq);
    CHECK_SEQ(tuPO->fseq);
    CHECK_TABLE(tuPO->ltm_tab);
    CHECK_FILE(tuPO->in);
    CHECK_NFILE(tuPO->out,tuPO->outname);
    CHECK_FREE(tuPO);
    return(TRUE);
}
/*****************************************************************************
*   Init data structure
*/
int InitTm_utilI(TM_UTIL *tuPO)
{
    VALIDATE(tuPO,TM_UTIL_ID);
    INIT_S(tuPO->inname);
    INIT_S(tuPO->outname);
    tuPO->in = tuPO->out = NULL;
    tuPO->iform = BOGUS;
    tuPO->oform = BOGUS;
    tuPO->con1 = tuPO->con2 = tuPO->con3 = DEF_NN_CONC;
    tuPO->salt = tuPO->mg = tuPO->temp = BAD_R;
    INIT_S(tuPO->parname);
    tuPO->pdc = 1;
    tuPO->do_rsd = FALSE;
    tuPO->quiet = FALSE;
    tuPO->do_dpar = FALSE;
    tuPO->do_ds = FALSE;
    tuPO->do_therm = tuPO->do_dtherm = FALSE;
    tuPO->do_cmb = FALSE;
    tuPO->do_fds = tuPO->do_den = FALSE;
    tuPO->do_tes = tuPO->do_btes = FALSE;
    tuPO->tlen = BAD_R;
    tuPO->tlp = 0;
    tuPO->do_tmpro = FALSE;
    tuPO->do_tab_lo = tuPO->do_tab_hi = 0;
    tuPO->do_tab_j = 1;
    tuPO->do_tabex = FALSE;
    tuPO->minv = BAD_R;
    tuPO->do_tlmin = tuPO->do_tlrev = FALSE;
    tuPO->do_flag = FALSE;
    tuPO->do_not = FALSE;
    tuPO->do_eraw = FALSE;
    tuPO->do_emin = tuPO->do_emid = tuPO->do_emax = FALSE;
    tuPO->do_padbad = TRUE;
    tuPO->firstb = tuPO->lastb = -1;
    tuPO->do_rre = FALSE;
    tuPO->do_rc = FALSE;
    ClearTm_utilThermoVals(tuPO);
    return(TRUE);
}
/*****************************************************************************
*   Clear per-seq, run-time thermo values 
*/
void ClearTm_utilThermoVals(TM_UTIL *tuPO)
{
    VALIDATE(tuPO,TM_UTIL_ID);
    tuPO->Tm = tuPO->dG = tuPO->dH = tuPO->dS = tuPO->fds = BAD_R;
    return;
}
/*****************************************************************************
*   Handle the -tlen options
*/
int HandleTmLenI(TM_UTIL *tuPO, FILE *outPF)
{
    int good,slen;
    DOUB tmD,utmD;
    char *seqPC, nameS[NSIZE];

    HAND_NFILE(outPF);
    FillSeqNameStringI(tuPO->seq,nameS,NSIZE);
    if(!GetSeqSeqI(tuPO->seq,&seqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(tuPO->seq);
    /***
    *   Check any offsets
    */
    if( (tuPO->tlp>slen) || (tuPO->tlp<0) ) {
        fprintf(outPF,"%s PROBLEM: BAD SUPPLIED OFFSET %d SeqLen=%d\n",
            nameS,tuPO->tlp,slen);
        return(FALSE);
    }
    /***
    *   If reverse, pass start of string but limit end
    */
    if(tuPO->do_tlrev) {
        good = SeqLenForThermI(tuPO->tm, seqPC, slen-tuPO->tlp, tuPO->tlen,
            REVERSE, &tmD, &utmD, tuPO->do_fds);
    }
    else {
        seqPC += tuPO->tlp;
        good = SeqLenForThermI(tuPO->tm, seqPC, slen-tuPO->tlp, tuPO->tlen,
            FORWARD, &tmD, &utmD, tuPO->do_fds);
    }
    /***
    *   Returned values are:
    *   good = length for over-shot thermo (if found)
    *   tmD = under-target Tm
    *   utmD = over-target Tm
    */
    if(good < 1) {
        fprintf(outPF,"%-12s PROBLEM: TOO SHORT SeqLen=%d Offset=%d\n",
            nameS,slen,tuPO->tlp);
        return(FALSE);
    }
    /**
    *   Output sequence; First pick one value from the two
    */
    if(tuPO->do_otls) {
        if(tuPO->do_tlmin) {
            tmD = utmD;
        }
        /***
        *   If lower is closest Tm to target, use this and shorten else use upper
        */
        else {
            if( (tuPO->tlen-tmD) < (utmD-tuPO->tlen) )
                good--;
            else
                tmD = utmD;
        }
        if(tuPO->do_tlrev) {
            NarrowSeqI(tuPO->seq,0,good,REVERSE,TRUE);
        }
        else {
            NarrowSeqI(tuPO->seq,0,good,FORWARD,TRUE);
        }
        fprintf(outPF,"# %-20s\t%d\t%6.3f\n",nameS,good,tmD);
        WriteSeq(tuPO->seq,SEQFM_RAW,outPF);
    }
    /***
    *   Report bounding lengths and thermo values
    */
    else {
        fprintf(outPF,"%-20s\t%d\t%d\t%6.3f\t%6.3f\n",nameS,good-1,good,tmD,utmD);
    }
    return(TRUE);
}
/**************************************************************************
*   Dump along-length thermo profile
*/
int HandleTmProfileI(TM_UTIL *tuPO, FILE *outPF)
{
    int j,slen;
    DOUB vD;
    char *seqPC, nameS[NSIZE];

    HAND_NFILE(outPF);
    FillSeqNameStringI(tuPO->seq,nameS,NSIZE);
    if(!GetSeqSeqI(tuPO->seq,&seqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(tuPO->seq);
    for(j=MIN_TM_LEN; j<=slen; j++)
    {
        SeqTmThermI(tuPO->tm, seqPC, j, &tuPO->Tm, &tuPO->dG,
            &tuPO->dH, &tuPO->dS);
        if(tuPO->do_fds) {
            vD = Fraction2D(tuPO->con1, tuPO->con2, tuPO->temp, tuPO->dG);
        }
        else {
            vD = tuPO->Tm;
        }
        fprintf(outPF,"%-20s\t%d\t%6.3f\n",nameS,j,vD); 
    }
    return(TRUE);
}
/**************************************************************************
*   Matrix of sequence offsets (starting pos) X min-to-max lengths
*/
int HandleThermoSampleTableI(TM_UTIL *tuPO, FILE *outPF)
{
    int i,j,slen,end,row,col,nok,fok,lok;
    DOUB vD;
    char *seqPC, nameS[NSIZE], rowS[NSIZE], seqS[MAX_TM_LEN + 1];

    HAND_NFILE(outPF);
    FillSeqNameStringI(tuPO->seq,nameS,NSIZE);
    if(!GetSeqSeqI(tuPO->seq,&seqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(tuPO->seq);
    end = slen - tuPO->do_tab_lo;
    if(end < 1) {
        return(FALSE);
    }
    fprintf(outPF,"# %s\n",nameS);
    /***
    *   Fill up table with thermo vals; First clear values and masks
    */
    InitTableValsI(tuPO->ltm_tab, BAD_D, FALSE);
    SetTableMasks(tuPO->ltm_tab, FALSE);
    row = 0;
    for(i=0; i<end; i += tuPO->do_tab_j) 
    {
        sprintf(rowS,"%s__%02d\t",nameS,i+1);
        SetTableRowLabI(tuPO->ltm_tab, row, rowS);
        col = 0;
        for(j = tuPO->do_tab_lo; j <= tuPO->do_tab_hi; j += tuPO->do_tab_j) 
        {
            ClearTm_utilThermoVals(tuPO);
            if(j<slen) {
                SeqTmThermI(tuPO->tm, seqPC, j, &tuPO->Tm, &tuPO->dG, &tuPO->dH, &tuPO->dS);
                if(tuPO->do_fds) {
                    tuPO->fds = Fraction2D(tuPO->con1, tuPO->con2, tuPO->temp, tuPO->dG);
                    vD = tuPO->fds;
                }
                else {
                    vD = tuPO->Tm;
                }
                AddTableValI(tuPO->ltm_tab, row, col, vD);
            }
            col++;
        }
        row++;
        seqPC++;
        slen--;
    }
    /***
    *   Not extracting seq = Dump the table
    */
    if(! tuPO->do_tabex) {
        /* Only actually dump used rows and columns; mask those */
        SetTableRowRangeMaskI(tuPO->ltm_tab, 0, row, TRUE);
        SetTableColRangeMaskI(tuPO->ltm_tab, 0, col, TRUE);
        /* Set name; Replace BAD_D with pretty-print number */
        SetTableNamesI(tuPO->ltm_tab, nameS, NULL, -1);
        InitTableMatchingValsI(tuPO->ltm_tab, BAD_D, -1.0, TRUE);
        /* Args: 0 = no header info; TRUE = use masking */
        DumpTable(tuPO->ltm_tab, 0, TRUE, outPF);
        return(TRUE);
    }
    /***
    *   Extracting seq(s) per row
    */ 
    GetSeqSeqI(tuPO->seq, &seqPC);
    row = 0;
    for(i=0; i<end; i += tuPO->do_tab_j) 
    {
        nok = GetTableRowFlagValIndicesI(tuPO, row, &fok, &lok);
        if(nok > 0) {
            /***
            *   Single special picks? If so, 'tighten loop', else given range 
            */
            if(tuPO->do_emin) {
                lok = fok+1;
            }
            else if(tuPO->do_emid) {
                fok = (fok + lok -1) / 2;
                lok = fok+1;
            }
            else if(tuPO->do_emax) {
                fok = lok-1;
            }
            for(col=fok; col<lok; col++) 
            {
                /***
                *   It's possible that values are bad; i.e. shorter seq has higher Tm
                */
                GetTableValI(tuPO->ltm_tab, row, col, &vD);
                if(BAD_DOUB(vD)) {
                    continue;
                }
                j = tuPO->do_tab_lo + col;
                strncpy(seqS, seqPC, j);
                seqS[j] = '\0';
                fprintf(outPF,"%s__%02d_%02d\t%s\t%3.3f\t%d\n",nameS,i+1,i+j+1,seqS,vD,j);
            }
        }
        row++;
        seqPC++;
    }
    return(TRUE);
}
/***************************************************************************
*   Find first, last and count of row elements with values in min-max bounds
*   Also sets out-of-bound values in table to bad
*/
int GetTableRowFlagValIndicesI(TM_UTIL *tuPO, int row, int *fPI, int *lPI)
{
    int n,col,first,last,ncol;
    DOUB vD;

    ncol = GetTableColsI(tuPO->ltm_tab, FALSE);
    n = 0;
    first = last = -1;
    for(col=0; col<ncol; col++)
    {
        GetTableValI(tuPO->ltm_tab, row, col, &vD);
        if( (vD >= tuPO->minv) && (vD <= tuPO->maxv) ){
            n++;
            if(n == 1){
                first = col;
            }
            last = col;
        }
        else {
            SetTableValI(tuPO->ltm_tab, row, col, BAD_D);
        }
    }
    if(n > 0) {
        *fPI = first;
        last = (last<0) ? first+1 : last+1;
        *lPI = last;
    }
    return(n);
}
/***************************************************************************
*   Returns length if got sequence
*   Returns BOGUS if ParseSeq is done
*/
int LoadTmutilSequenceI(TM_UTIL *tuPO, int n)
{
    int fok,sok;

    /***
    *   Two explicit strands = sequential seqs;
    *       Raw format, no cleaning, report errors
    */
    if( (tuPO->do_tes) || (tuPO->do_btes) ) {
        fok = ParseSeqI(tuPO->in, tuPO->iform, n, FALSE, TRUE, tuPO->seq);
        sok = ParseSeqI(tuPO->in, tuPO->iform, n, FALSE, TRUE, tuPO->sseq);
    }
    else {
        fok = sok = ParseSeqI(tuPO->in, tuPO->iform, n, FALSE, TRUE, tuPO->seq);
    }
    /***
    *   End of file, return BOGUS to break out of loop
    */
    if( (fok==FALSE) || (sok==FALSE) ) {
        return(BOGUS);
    }
    else if( (fok==BOGUS) || (sok==BOGUS) ) {
        printf("Problem loading sequence\n");
        return(FALSE);
    }
    /***
    *   Rev comp? Flip input here; Not allowed with two exp seqs
    */
    if(tuPO->do_rc) {
        ReverseCompSeqSeqI(tuPO->seq);
    }
    /***
    *   Length; If two explicit seqs, need to check this   
    */
    fok = sok = GetSeqLenI(tuPO->seq);
    HandleTmuSubseqI(tuPO, tuPO->seq, tuPO->fseq);
    if( (tuPO->do_tes) || (tuPO->do_btes) ) {
        sok = GetSeqLenI(tuPO->sseq);
        HandleTmuSubseqI(tuPO, tuPO->sseq, tuPO->fseq);
    }
    /*
    *   Return length (minimum) of seq
    */
    return(MIN_NUM(fok,sok));
}
/***************************************************************************
*   Call Nic Peyret's Dset thermodynamics code
*   Copy resulting numbers into fields in passed TM_UTIL struct
*/
int CallDsetTwoSeqsI(TM_UTIL *tuPO)
{
    int mode,slen;
    DOUB tmDA[DSET_EV_NUM],gDA[DSET_EV_NUM],hDA[DSET_EV_NUM],sDA[DSET_EV_NUM];
    DOUB fDA[DSET_EV_NUM];
    char *seqPC, *sseqPC, seqS[MAX_TM_LEN];

    if(!GetSeqSeqI(tuPO->seq,&seqPC)) {
        return(FALSE);
    }
    /***
    *   If delta energies, need single sequence alone first 
    */
    if(tuPO->do_dtherm) {
        if( ! SeqDsetEnergyI(tuPO->tm,seqPC,NULL,TRUE,DSET_IMP,tmDA,gDA,hDA,sDA,fDA)) {
            return(FALSE);
        }
        tuPO->Tm = tmDA[0];
        tuPO->dG = gDA[0];
        tuPO->dS = sDA[0];
        tuPO->dH = hDA[0];
        tuPO->fds = fDA[0];
    }
    /***
    *   Second sequence
    *   If either has "/", then coaxial case
    *   -btes option needs second sequence reverse complimented;
    *   Nic's code expects 5' >--> 3' + 3' >--> 5' for two seqs
    */
    if(!GetSeqSeqI(tuPO->sseq,&sseqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(tuPO->sseq);
    if( (strchr(seqPC,'/')) || (strchr(sseqPC,'/')) ) {
        mode = DSET_EXPCOAX;
    }
    else {
        mode = DSET_EXPNOCOAX;
    }
    if(tuPO->do_btes) {
        CompDNASeqI(sseqPC,slen,seqS);
    }
    else {
        sprintf(seqS,"%s",sseqPC);
    }
    seqS[slen] = '\0';
    /***
    *   Two explicit sequence call
    */
    if ( ! SeqDsetEnergyI(tuPO->tm,seqPC,seqS,TRUE,mode,tmDA,gDA,hDA,sDA,fDA)) {
        return(FALSE);
    }
    /***
    *   If delta energies, store both baseline and two-strand values
    */
    if(tuPO->do_dtherm) {
        tuPO->Tm2 = tmDA[0];
        tuPO->dG2 = gDA[0];
        tuPO->dS2 = sDA[0];
        tuPO->dH2 = hDA[0];
        tuPO->fds2 = fDA[0];
    }
    else {
        tuPO->Tm = tmDA[0];
        tuPO->dG = gDA[0];
        tuPO->dS = sDA[0];
        tuPO->dH = hDA[0];
        tuPO->fds = fDA[0];
    }
    return(TRUE);
}
/***************************************************************************
*   Open files or die
*/
int OpenTmutilFilesI(TM_UTIL *tuPO)
{
    if(!(tuPO->in=OpenUFilePF(tuPO->inname,"r",NULL))) {
        return(FALSE);
    }
    if(!NO_S(tuPO->outname)) {
        if(!(tuPO->out=OpenUFilePF(tuPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    HAND_NFILE(tuPO->out);
    return(TRUE);
}
/***************************************************************************
*   Set options to tm_util and check consistency 
*   TRUE if all looks OK
*   FALSE if any problem / inconsistency
*/
int CheckTmutilOptionsI(TM_UTIL *tuPO)
{
    char bufS[DEF_BS];

    if(!SetTmParOptionsI(tuPO)) {
        PROBLINE;
        printf("Problem setting options Thermo\n");
        return(FALSE);
    }
    /***
    *   Get the actually used values into local data structure
    */
    GetTmParConcI(tuPO->tm, &tuPO->con1, &tuPO->con2, &tuPO->con3);
    GetTmParTempI(tuPO->tm, &tuPO->temp);
    GetTmParSaltI(tuPO->tm, &tuPO->salt);
    GetTmParMgI(tuPO->tm, &tuPO->mg);
    /***
    *   Restricted bases?
    */
    if(tuPO->firstb > 0) {
        if( (tuPO->firstb<1) || (tuPO->firstb>=tuPO->lastb) ) {
            PROBLINE;
            printf("Bad base range: %d %d\n", tuPO->firstb, tuPO->lastb);
            return(FALSE);
        }
    }
    /***
    *   Sham? Delta values = full dump
    */
    if( tuPO->do_dtherm ) {
        tuPO->do_fds = TRUE;
    }
    /***
    *   Can't do competitive binding if not doing delta binding
    */
    if( (tuPO->do_cmb) && (!tuPO->do_dtherm) ) {
        PROBLINE;
        printf("Competitive miss-match binding only works with -dthe\n");
        return(FALSE);
    }
    /***
    *   Check for algoritm specific options
    */
    if( tuPO->do_fds && (ConcTermsForTmAlgoI(tuPO->tm)<2) ) {
        PROBLINE;
        FillTmAlgoString(tuPO->tm,bufS);
        printf("Fraction double stranded needs 2-conc algorithm\n"); 
        printf("  This algorithm won't work:\n  %s\n",bufS);
        return(FALSE);
    }
    /***
    *   Two explicit strand restrictions
    */
    if( (tuPO->do_tes) || (tuPO->do_btes) ) {
        if(tuPO->algo != ALGO_PEYRET) {
            PROBLINE;
            printf("Two explicit seqs only work's with -tmpey\n"); 
            FillTmAlgoString(tuPO->tm,bufS);
            printf("Current algorithm won't work:\n  %s\n",bufS);
            return(FALSE);
        }
        if(tuPO->do_tmpro) {
            PROBLINE;
            printf("Two explicit seqs don't work with Tm profile\n");
            return(FALSE);
        }
        if(!BAD_REAL(tuPO->tlen)) {
            PROBLINE;
            printf("Two explicit seqs don't work with Tm target length\n");
            return(FALSE);
        }
        if(tuPO->do_rc) {
            PROBLINE;
            printf("Two explicit seqs don't work with rev comp\n");
            return(FALSE);
        }
    }
    /***
    *   Dangling ends?
    */
    if( tuPO->do_den ) {
        if(tuPO->algo != ALGO_PEYRET) {
            PROBLINE;
            printf("Dangling ends only work's with -tmpey\n"); 
            FillTmAlgoString(tuPO->tm,bufS);
            printf("  This algorithm won't work:\n  %s\n",bufS);
            return(FALSE);
        }
    }
    /***
    *   Delta Therm?
    */
    if( (tuPO->do_dtherm) && ( ! ((tuPO->do_tes) || (tuPO->do_btes)) ) ) {
        PROBLINE;
        printf("Delta thermodynamics (-dthe) only works with -tes or -btes\n");
        return(FALSE);
    }
    /***
    *   Competitive missmatch binding
    */
    if( (tuPO->do_cmb ) && ( ! ((tuPO->do_tes) || (tuPO->do_btes)) ) ) {
        PROBLINE;
        printf(
        "Competitive missmatch binding (-cmb) only works with -tes or -btes\n");
        return(FALSE);
    }
    /***
    *   Target length with bogus position?
    */
    if(!BAD_REAL(tuPO->tlen)) {
        if(tuPO->tlp<0) {
            PROBLINE;
            printf("Bogus -tlp argument supplied = %d [0,1,2...]\n",tuPO->tlp);
            return(FALSE);
        }
        /***
        *    SHAM; should work (and extend to dG)
        */
        if(tuPO->algo==ALGO_PEYRET) {
            PROBLINE;
            printf("Sorry, -tlen won't work with -tmpey\n"); 
            return(FALSE);
        }
    }
    /***
    *   Flagging?
    */
    if(!BAD_REAL(tuPO->minv)) {
        if(tuPO->maxv < tuPO->minv) {
            PROBLINE;
            printf("Impossible flagging range given %3e to %3e\n",
                tuPO->minv, tuPO->maxv);
            return(FALSE);
        }
        tuPO->do_flag = TRUE;
    }
    /***
    *   Thermo table
    */
    if( tuPO->do_tab_lo ) {
        if( (tuPO->do_tab_lo < MIN_TM_LEN ) || (tuPO->do_tab_hi <= tuPO->do_tab_lo) ) {
            PROBLINE;
            printf("Table min/max won't work: %d and %d\n",
                tuPO->do_tab_hi, tuPO->do_tab_lo);
            return(FALSE);
        }
        if(tuPO->do_tab_hi > MAX_TM_LEN) {
            PROBLINE;
            printf("Table max is too long: %d Max Tm len = %d\n",
                tuPO->do_tab_hi, MAX_TM_LEN);
            return(FALSE);
        }
    }
    LIMIT_NUM(tuPO->do_tab_j, 1, 100);
    /***
    *   Special case with table-sampling + flagging + extract
    *   First, special case to make sure do_eraw is set if do_emin, mid, max ard set
    */
    if(tuPO->do_emin || tuPO->do_emid || tuPO->do_emax) {
        tuPO->do_eraw = TRUE;
        tuPO->do_not = FALSE;
    }
    if( (tuPO->do_tab_lo) && (! BAD_REAL(tuPO->minv)) && (tuPO->do_eraw) ) {
        tuPO->do_tabex = TRUE;
    }
    /***
    *   Incompatable options
    */
    if(tuPO->algo == ALGO_24) {
        tuPO->do_therm = FALSE;
    }
    return(TRUE);
}
/***************************************************************************
*   Set options to real TM structure and check consistency 
*   TRUE if all looks OK
*   FALSE if any problem / inconsistency
*/
int SetTmParOptionsI(TM_UTIL *tuPO)
{
    /***
    *   What algorithm / parameter defaults to use?
    */
    SetTmParAlgo(tuPO->tm,tuPO->algo,TRUE); 
    if(tuPO->algo != GetTmParAlgoI(tuPO->tm)) {
        PROBLINE;
        printf("Problem setting algorithm (%d)\n",tuPO->algo);
        return(FALSE);
    }
    /***
    *   Parameters to load from file; these overide any set above
    */
    if(!NO_S(tuPO->parname)) {
        if(!LoadTmParNNEnergiesI(tuPO->tm,tuPO->parname)) {
            return(FALSE);
        }
    }
    /***
    *   Set concentration / buffer terms?
    */
    SetTmParConc(tuPO->tm,tuPO->con1,tuPO->con2,tuPO->con3);    
    if(!BAD_REAL(tuPO->salt)) {
        SetTmParSalt(tuPO->tm,tuPO->salt);  
    }
    if(!BAD_REAL(tuPO->mg)) {
        SetTmParMg(tuPO->tm,tuPO->mg);  
    }
    if(!BAD_REAL(tuPO->temp)) {
        SetTmParTemp(tuPO->tm,tuPO->temp);  
    }
    if(tuPO->do_rsd) {
        SetTmStrandDir(tuPO->tm,REVERSE);   
    }
    if(tuPO->pdc > 0) {
        SetTmParDefCol(tuPO->tm, -1, tuPO->pdc);
    }
    return(TRUE);
}
/**************************************************************************
*   Dump settings
*/
void ReportTmutilSettings(TM_UTIL *tuPO, FILE *outPF)
{
    char hostS[DEF_BS],osS[DEF_BS],verS[DEF_BS];

    HAND_NFILE(outPF);
    fprintf(outPF,"# Output from %s\n",VERSION_S);
    fprintf(outPF,"#   %s %s %s\n",BD_S,__DATE__,__TIME__);
    GetSystemInfo(NULL,hostS,osS,verS,NULL);
    /* printf(outPF,"#   Run on %s, %s (%s)\n",hostS,osS,verS); */
    fprintf(outPF,"#   Run on %s, %s\n",hostS,osS);
    fprintf(outPF,"# Input file  %s\n",tuPO->inname);
    if(tuPO->do_rc) {
        fprintf(outPF,"#  Seqs converted to Reverse Complement\n");
    }
    ReportFlaggedOutputI(tuPO,outPF);
    if(tuPO->firstb > 0) {
        fprintf(outPF,"#   Base-Range:    %d to %d\t", tuPO->firstb, 
            tuPO->lastb);
        if(tuPO->do_rre) {
            fprintf(outPF,"(from end)\n");
        }
        else {
            fprintf(outPF,"(from start)\n");
        }
    }
    DumpTmPars(tuPO->tm, "# ", tuPO->do_dpar, outPF);
}
/**************************************************************************
*   Dump any flagging criteria
*/
void ReportFlaggedOutputI(TM_UTIL *tuPO, FILE *outPF)
{
    if(!tuPO->do_flag) {
        return;
    }
    if(tuPO->do_fds) {
        fprintf(outPF,"# %% Double strand %3.3f to %3.3f\n",tuPO->minv,tuPO->maxv);
    }
    else {
        fprintf(outPF,"# Tm limits   %3.3f to %3.3f\n",tuPO->minv,tuPO->maxv);
    }
    fprintf(outPF,"#\n");
}
/**************************************************************************
*   Dump column headings
*/
void ReportTmutilOutput(TM_UTIL *tuPO, FILE *outPF)
{
    int i;

    HAND_NFILE(outPF);
    if( tuPO->do_tab_lo ) {
        fprintf(outPF,"# Thermo values for lengths %d to %d\n", tuPO->do_tab_lo, tuPO->do_tab_hi);
        if( tuPO->do_tabex ) {
            fprintf(outPF,"# Extracting seqs for values %6.3f to %6.3f\n", tuPO->minv, tuPO->maxv);
        }
        else {
            fprintf(outPF,"Sequence\t");
            for(i=tuPO->do_tab_lo; i <= tuPO->do_tab_hi; i += tuPO->do_tab_j) {
                fprintf(outPF,"%d ",i);
            }
            fprintf(outPF,"\n");
        }
        return;
    }
    /***
    *   What will be dumped?
    */
    if(!BAD_REAL(tuPO->tlen)) {
        fprintf(outPF,"# %-20s\tLength\tTm\n","Name");
        return;
    }
    fprintf(outPF,"# %-20s\tTm","Name");
    if(tuPO->do_fds) {
        fprintf(outPF,"\tFB");  
    }
    if(tuPO->do_therm) {
        fprintf(outPF,"\tdG\tdH\tdS");
    }
    /***
    *   Delta values?
    */
    if(tuPO->do_dtherm) {
        if(!tuPO->do_therm) {
            fprintf(outPF,"\tdG");
        }
        fprintf(outPF,"\tTm2");
        fprintf(outPF,"\tFB2"); 
        fprintf(outPF,"\tdG2");
        fprintf(outPF,"\tdTm");
        fprintf(outPF,"\tdFB"); 
        fprintf(outPF,"\tddG");
    }
    if(tuPO->do_cmb) {
        fprintf(outPF,"\tcFB1");    
        fprintf(outPF,"\tcFB2");    
        fprintf(outPF,"\tcdFB");    
    }
    fprintf(outPF,"\n");    
}
/**************************************************************************/
int SeqFlaggingValI(TM_UTIL *tuPO, DOUB *valPD)
{
    int good;
    DOUB valD;

    if(!tuPO->do_flag) {
        return(TRUE);
    }
    if(tuPO->do_fds) {
        valD = tuPO->fds;
    }
    else {
        valD = tuPO->Tm;
    }
    good = ( (valD >= tuPO->minv)&&(valD <= tuPO->maxv) );
    if(tuPO->do_not) {
        good = !good;
    }
    if(valPD) {
        *valPD = valD;
    }
    return(good);
}
/**************************************************************************
*
*/
int HandleTmutilOutputI(TM_UTIL *tuPO, int good, char *seqS, int slen,
    char *nameS, FILE *outPF)
{
    char *fseqPC;

    HAND_NFILE(outPF);
    /***
    *   Sequence output?
    */
    if(tuPO->do_eraw) {
        if(good) {
            fprintf(outPF,"%-20s\t%s\n",nameS,seqS);
        }
        return(good);
    }
    if(tuPO->do_flag) {
        fprintf(outPF,"%d\t",good);
    }
    fprintf(outPF,"%-20s\t%6.3f",nameS,tuPO->Tm); 
    if(tuPO->do_fds) {
        fprintf(outPF,"\t%6.3f",tuPO->fds); 
    }
    if(tuPO->do_therm) {    
        fprintf(outPF,"\t%6.3f\t%6.3f\t%6.3f",tuPO->dG,tuPO->dH,tuPO->dS);  
    }
    /***
    *   Delta thermo?
    */
    if(tuPO->do_dtherm) {
        if(!tuPO->do_therm) {
            fprintf(outPF,"\t%6.3f",tuPO->dG);
        }
        /***
        *   Imperfect case
        */
        fprintf(outPF,"\t%6.3f",tuPO->Tm2);
        fprintf(outPF,"\t%6.3f",tuPO->fds2);
        fprintf(outPF,"\t%6.3f",tuPO->dG2);
        /***
        *   Delta baseline - imperfect case
        */
        fprintf(outPF,"\t%6.3f",tuPO->Tm - tuPO->Tm2);
        fprintf(outPF,"\t%6.3f",tuPO->fds - tuPO->fds2);
        fprintf(outPF,"\t%6.3f",tuPO->dG - tuPO->dG2);
    }
    /***
    *   Competitive miss-match, need to calculate here and report
    */
    if(tuPO->do_cmb) {
        if(!HandleCompetitiveFbOutputI(tuPO,outPF)) {
            ReportDsetProblemSeq(nameS,outPF);
        }
    }
    /*  Dump input seq? */
    if(tuPO->do_ds) {
        if(tuPO->firstb > 0) {
            GetSeqSeqI(tuPO->fseq, &fseqPC);
        }
        else {
            GetSeqSeqI(tuPO->seq, &fseqPC);
        }
        fprintf(outPF,"\t%s", fseqPC);
    }
    /***
    *   End of line, all done
    */
    fprintf(outPF,"\n");    
    return(good);
}
/***************************************************************************
*   Competitive calc and output
*/
int HandleCompetitiveFbOutputI(TM_UTIL *tuPO,FILE *outPF)
{
    DOUB fb1D,fb2D;
    
    HAND_NFILE(outPF);
    if( !CalcCompEqI(tuPO->con3, tuPO->con1, tuPO->con2, tuPO->dG, tuPO->dG2,
            tuPO->temp, &fb1D, &fb2D) ) {
        return(FALSE);
    }
    fb1D *= 100.0;
    fb2D *= 100.0;
    fprintf(outPF,"\t%6.3f",fb1D);
    fprintf(outPF,"\t%6.3f",fb2D);
    fprintf(outPF,"\t%6.3f",fb1D - fb2D);
    return(TRUE);
}
/***************************************************************************
*   Problem with Dset
*/
void ReportDsetProblemSeq(char *nameS,FILE *outPF)
{
    HAND_NFILE(outPF);
    fprintf(outPF,"# Problem with %s (dset failed)\n",nameS);
}
/***************************************************************************
*   Shrink current sequence to given base range (if any)
*/
int HandleTmuSubseqI(TM_UTIL *tuPO, SEQ *seqPO, SEQ *fseqPO)
{
    int len, flen;

    if(tuPO->firstb > 0) {
        len = tuPO->lastb - tuPO->firstb + 1;
        CopySeqI(seqPO, fseqPO, -1, -1);
        SetCaseSeqSubseqI(fseqPO, FALSE, -1, -1);
        if(tuPO->do_rre) {
            NarrowSeqI(seqPO,tuPO->firstb-1,len,REVERSE,TRUE);
            flen = GetSeqLenI(fseqPO);
            SetCaseSeqSubseqI(fseqPO, TRUE, flen - tuPO->lastb, flen - tuPO->firstb + 1);
        }
        else {
            NarrowSeqI(seqPO,tuPO->firstb-1,len,FORWARD,TRUE);
            SetCaseSeqSubseqI(fseqPO, TRUE, tuPO->firstb-1, tuPO->lastb);
        }
    }
    return(TRUE);
}
