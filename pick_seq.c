/*
* pick_seq.c
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
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "dna_pair.h"
#include "score.h"
#include "pick_seq.h"

#define DB_PS if(DB[182])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit(AllDoneI(PickSeqI(argc,argv),NULL)); }
/**************************************************************************/
void PickSeqUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>   Input sequences to pick from\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -out XXX   Send output to file XXX\n");
    printf("   -num #     Pick # sequences\n");
    printf("   -swap      Use the monte carlo swap picking algorithm\n");
    printf("   -seed #    Set random number seed to #\n");
    printf("   -cyc #     Set swap cycles to # (Default %d)\n", DEF_PS_CYC);
    printf("   -rep #     Set replacement tries to # (Default %d)\n", DEF_PS_REP);
    printf("   -sirf #    Save intermediate results every # cycle\n");
    printf("   -siru      Save intermediate results unique name / cycle\n");
    printf("   -smat # #  Similar & Comp total match score factor to #\n");
    printf("   -scon # #  Similar & Comp contig match score factor to #\n");
    printf("   -scb # #   Similar & Comp block match score factor to #\n");
    printf("   -mswm      Score simple match via weighted matching\n");
    printf("   -cswm      Score contig match via weighted matching\n");
    printf("   -ham       Score via Hamming distance (Similarity only)\n");
    printf("   -mwf XXX   Match weights from XXX (AA X, AC Y,... per line)\n");
    printf("   -fix XXX   Fix (keep) a subset of seqs listed in file XXX\n");
    printf("   -bwn #     Bail (stop) cycling if number of 'worst' is this or more\n");
    printf("   -bwf #     Bail (stop) cycling if fraction 'worst' is this or more\n");
    printf("   -quiet     Suppress run-time status reporting\n");
}
/**************************************************************************
*   main function
*/
int PickSeqI(int argc, char **argv)
{
    int swap,ifas,iraw,iseq;
    PICKSEQ *psPO;

    psPO = CreatePickseqPO();
    swap = ifas = iraw = iseq = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -seed I -num I -swap B -cyc I -rep I -sirf I -siru B\
        -smat R2 -scon R2 -scb R2 -mswm B -cswm B -mwf S -ifas B -iraw B\
        -quiet B -fix S -iseq B -ham B -bwn I -bwf R",
        psPO->seqname, psPO->outname, &psPO->rseed, &psPO->num, &swap, 
        &psPO->cyc, &psPO->rep, &psPO->sirf, &psPO->siru,
        &psPO->smat, &psPO->cmat, &psPO->scon, &psPO->ccon, 
        &psPO->scb, &psPO->ccb, &psPO->do_mswm, &psPO->do_cswm, 
        &psPO->mwfname, &ifas, &iraw, &psPO->verbose, &psPO->fixname,
        &iseq, &psPO->do_ham, &psPO->do_bwn, &psPO->do_bwf,
        (int *)NULL))
    {
        PickSeqUse();
        CHECK_PICKSEQ(psPO);
        return(FALSE);
    }
    /***
    *   Set input format 
    */
    psPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas,psPO->seqname,TRUE);
    if(!psPO->iform) {
        printf("Problem with input seq(s)\n");
        CHECK_PICKSEQ(psPO);
        return(FALSE);
    }
    /***
    *   Set command line options and check
    */
    if(swap) {
        psPO->alg = PS_ALG_SWAP;
    }
    if(!CheckPsGenericOptionsI(psPO)) {
        ABORTLINE;
        CHECK_PICKSEQ(psPO);
        return(FALSE);
    }
    /***
    *   Set up everything needed to go
    */
    if(!RealizePickseqI(psPO)) {
        ABORTLINE;
        CHECK_PICKSEQ(psPO);
        return(FALSE);
    }
    /***
    *   Check run-time options
    */
    if(!CheckPsSeqOptionsI(psPO)) {
        ABORTLINE;
        CHECK_PICKSEQ(psPO);
        return(FALSE);
    }
    /***
    *   Handle output
    */
    if(!SetPickseqOutputI(psPO)) {
        ABORTLINE;
        return(FALSE);
    }
    /***
    *   Do selection algrothim
    */
    switch(psPO->alg)
    {
        case PS_ALG_SWAP: 
            MonteCarloSwapPickI(psPO);
            break;
        default:
            PROBLINE;
            printf("Bogus algorithm code = %d\n",psPO->alg);
            printf("Not going to do anything!!!!\n");
            break;
    }
    /***
    *   All done
    */
    WriteSeqset(psPO->seqs,SEQFM_RAW,TRUE,psPO->out);
    CHECK_PICKSEQ(psPO);
    printf("# All done\n");
    return(TRUE);
}
/*****************************************************************************
*   Load sequences and set up memory 
*/
int RealizePickseqI(PICKSEQ *psPO)
{
    /***
    *   Load sequences
    */
    if(IS_BOG(psPO->iform)) {
        psPO->iform = GuessSeqFileTypeI(psPO->seqname,FALSE);
    }
    if(!ReadInSeqsetI(psPO->seqname,psPO->iform,0,&psPO->seqs,TRUE)) {
        printf("Failed to load seqs\n");
        return(FALSE);
    }
    psPO->nseq = GetSeqsetNumI(psPO->seqs);
    if(psPO->verbose) {
        printf("# Have %d sequences\n",psPO->nseq);
    }
    /***
    *   Random seed
    */
    Srand(psPO->rseed);
    /***
    *   Add seq-size space
    */
    psPO->tmask = (char *)ALLOC(psPO->nseq,sizeof(char));
    if(!psPO->tmask) {
        PROBLINE;
        printf("Failed to allocate temp space (%d seqs)\n",psPO->nseq);
        return(FALSE);
    }
    /***
    *   Any fixed subset to mask?
    */
    if(!NO_S(psPO->fixname)) {
        if(!SetFixedSubsetI(psPO)) {
            PROBLINE;
            printf("Failed to load fixed subset from file: %s\n",psPO->fixname);
            return(FALSE);
        }
    }
    psPO->rnum = psPO->num - psPO->nfix;
    /***
    *   Add num-size space
    */
    psPO->recs = CreateScorecsPO(psPO->num);
    psPO->tscores = (REAL *)ALLOC(psPO->num,sizeof(REAL));
    if( (!psPO->recs) || (!psPO->tscores) ) {
        PROBLINE;
        printf("Failed to allocate temp space (num %d)\n",psPO->num);
        return(FALSE);
    }
    /***
    *   Any match weighting file?
    */
    if(!NO_S(psPO->mwfname)) {
        if(!LoadPairingParsI(psPO->mwfname,psPO->pp)) {
            printf("Failed to load weights from %s\n",psPO->mwfname);
            return(FALSE);
        }
    }
    /***
    *   Hamming distance
    */
    if(psPO->do_ham){
        SetPparsHammingI(psPO->pp, psPO->do_ham);
    }
    /***
    *   Bailing number or fraction?
    */
    if(psPO->do_bwn > 0){
        psPO->do_bwf = RNUM(psPO->do_bwn) / RNUM(psPO->num);
    } 
    else if(psPO->do_bwf > 0.0){
        psPO->do_bwn = INT(psPO->do_bwf * RNUM(psPO->num));
    } 
    else{
        psPO->do_bwn = psPO->num;
    }
    return(TRUE);
}
/*****************************************************************************
*   Attempt to allocate fixed map and load it from named file
*/
int SetFixedSubsetI(PICKSEQ *psPO)
{
    int i;
    DOUB nD;
    char bufS[NSIZE],nameS[NSIZE];
    FILE *fPF;

    if(!(fPF = OpenUFilePF(psPO->fixname,"r",NULL))) {
        return(FALSE);
    }
    /***
    *   Allocate
    */
    psPO->fixmap = (char *)ALLOC(psPO->nseq,sizeof(char));
    if(!psPO->fixmap) {
        PROBLINE;
        printf("Failed to allocate fix subset space (%d seqs)\n",psPO->nseq);
        return(FALSE);
    }
    /***
    *   Mark each guy listed in file
    */
    while(fgets(bufS,NSIZE,fPF) != NULL)
    {
        if(COM_LINE(bufS)) {
            continue;
        }
        INIT_S(nameS);
        sscanf(bufS,"%s",nameS);
        if(NO_S(nameS)) {
            continue;
        }
        i = FindNamedSeqInSeqsetI(psPO->seqs,nameS,TRUE,NULL,NULL);
        if(IS_BOG(i)) {
            PROBLINE;
            printf("Failed to find listed name in sequence collection\n");
            printf("   Name: |%s|\n",nameS);
            CHECK_FILE(fPF);
            return(FALSE);
        }
        psPO->fixmap[i] = 1;
    }
    CHECK_FILE(fPF);
    /***
    *   Set count of masked seqs and check
    */
    ArraySum(psPO->fixmap,IS_CHAR,0,psPO->nseq,&nD);
    psPO->nfix = INT(nD);
    if(psPO->nfix<1) {
        return(FALSE);
    }
    if(psPO->verbose) {
        printf("# Set %d sequences fixed from %s\n",psPO->nfix,psPO->fixname); 
    }
    return(TRUE);
}
/*****************************************************************************
*   Check options that don't depend on number of sequences
*/
int CheckPsGenericOptionsI(PICKSEQ *psPO)
{
    if(IS_BOG(psPO->alg)) {
        PROBLINE;
        printf("Bad algorithm specified\n");
        return(FALSE);
    }
    if(psPO->num < 1) {
        printf("Number to pick is too small (%d)\n",psPO->num);
        return(FALSE);
    }
    return(TRUE);
}
/*****************************************************************************
*   Check that options are ok
*/
int CheckPsSeqOptionsI(PICKSEQ *psPO)
{
    if(psPO->num >= psPO->nseq) {
        printf("All sequences qualify; no picking needed\n");
        return(FALSE);
    }
    if(psPO->nfix >= psPO->num) {
        printf("Number of fixed seqs is too big\n");
        printf("  Number fixed:   %d\n",psPO->nfix);
        printf("  Number to pick: %d\n",psPO->num);
        return(FALSE);
    }
    return(TRUE);
}
/*****************************************************************************
*   Set up output file
*/
int SetPickseqOutputI(PICKSEQ *psPO)
{
    if(!NO_S(psPO->outname)) {
        if(!(psPO->out = OpenUFilePF(psPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    HAND_NFILE(psPO->out);
    WritePickSeqHeader(psPO,psPO->out);
    return(TRUE);
}
/*****************************************************************************
*   Pick subset by iteratively improving a random start set a number of cycles
*/
int MonteCarloSwapPickI(PICKSEQ *psPO)
{
    int i,n,lastie,rep,nsc;
    SEQSET *ssPO;
    SCOREC *recsPO;

    DB_PS DB_PrI(">> MonteCarloSwapPickI cyc=%d num=%d\n",psPO->cyc,psPO->num);
    ssPO = psPO->seqs;
    recsPO = psPO->recs;
    /***
    *   Start starting sequence masking; 
    *   seqs that are masked are those that are fixed and/or randomly chosen
    */
    SetStartingSetI(psPO,psPO->fixmap,ssPO->mask);
    /***
    *   Can't try more replacements than there are to try 
    */
    psPO->mtry = psPO->rep;
    if(psPO->mtry > (psPO->nseq - psPO->num)) {
        psPO->mtry = psPO->nseq - psPO->num;
    }
    DB_PS DB_PrI("+ rep=%d mtry=%d\n",psPO->rep,psPO->mtry);
    /***    
    *   Tell the story
    */
    SetInitialInternalScoresI(psPO);
    printf("# Evaluating %d sequences; %d cycles, %d replace attempts/bad/cyc\n",
        psPO->num,psPO->cyc,psPO->mtry);
    if(psPO->do_bwn < psPO->num) {
        printf("# Threshold to bail if num / frac worst-case records too high: %d / %0.2f\n",
            psPO->do_bwn, psPO->do_bwf);
    }
    fflush(stdout);
    /***
    *   Run through n cycles
    */
    n = lastie = 0;
    while(n < psPO->cyc)
    {
        DB_PS DB_PrI("+++++++ Cycle %d +++++++\n",n+1);
        /***
        *   Get best and worst in current score set
        */
        DB_PS DB_PrI("+ Tallying worsts\n");
        FindChosenBestWorstStatsI(psPO,n+1);
        /***
        *   If all selected seqs tie for worst, up the count of tie cases
        *   Only allow one cycle with tie case, then give up
        */
        if(psPO->tie == psPO->num) {
            printf("# All records have equivalent scores; DONE\n");
            break;
        }
        /***
        *   If the number of worst-scoring records meets threshold, bail
        */
        if(psPO->tie >= psPO->do_bwn) {
            printf("# Number / fraction of worst-case records too high (%d / %0.2f)\n",
                psPO->do_bwn, psPO->do_bwf);
            break;
        }
        /***
        *   Attempt to replace each worst-scoring record(s)
        */
        printf("#  Attempting to replace %3d worst records\n",psPO->tie);
        fflush(stdout);
        rep = nsc = 0;
        for(i=0;i<psPO->num;i++)
        {
            DB_PS DB_PrI("+ [%d] id=%d sc=%1.3f",i,recsPO[i].id,recsPO[i].sc);
            if(recsPO[i].sc < psPO->worst) {   
                DB_PS DB_PrI(" Not worst\n");
                continue;   
            }
            DB_PS DB_PrI(" WORST\n");
            /***
            *   Fixed?
            */
            if( (psPO->fixmap) && (psPO->fixmap[recsPO[i].id]) ) {
                DB_PS DB_PrI(" Fixed\n");
                continue;   
            }
            if(SwapThisGuyOutI(psPO, i)) {
                rep++;
            }
        }
        /***
        *   If failed to improve score (no replacements), then done
        */
        if(rep < 1) {
            printf("\n");
            printf("# Couldn't find any replacements; DONE\n");
            break;
        }
        /***
        *
        */
/**
        ravR /= RNUM(rep);
        printf("#  Av. Rep. Sc.  %f\n",ravR);
**/
        printf("#    %3d (%5.3f%%) replaced\n",rep,PERCENT_R(rep,psPO->tie));
        fflush(stdout);
        n++;
        /***
        *   Saving temp results?
        */
        if( (psPO->sirf>0) && ((n%psPO->sirf)==0) ) {
            DB_PS DB_PrI("+ Handling temp results\n");
            if(!HandleTempResultsI(psPO,n)) {
                ABORTLINE;
                break;
            }
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Try to swap out rep from chosen subset
*/
int SwapThisGuyOutI(PICKSEQ *psPO, int rep)
{
    int j,sub,sid,oid,try,nrep;
    REAL ravR,maxR,scR;
    SEQSET *ssPO;
    SCOREC *recsPO;

    ssPO = psPO->seqs;
    recsPO = psPO->recs;
    /***
    *   Clear temp mask of already attempted swapable seqs
    */
    for(j=0;j<psPO->nseq;j++)
    {   
        psPO->tmask[j] = 0; 
    }
    /***
    *   Up to max tries
    */
    ravR = 0.0;
    try = nrep = 0;
    while(try< psPO->mtry)
    {
        /***
        *   Random substitute as long as not currently in pool
        *   and not already tried as a sub
        */
        sub = RandI(psPO->nseq);
        if( ssPO->mask[sub] || psPO->tmask[sub] ) {     
            DB_PS DB_PrI("+  s=%d (in[%d] or tried[%d])\n",sub,
                ssPO->mask[sub],psPO->tmask[sub]);
            continue; 
        }
        /***
        *   Going to try sub
        */
        try++;
        psPO->tmask[sub] = 1;
        DB_PS DB_PrI("+ try%d s=%d\n",try,sub);
        /***
        *   Compare sub to all others in chosen set except the 
        *   one we're replacing
        */
        maxR = 0.0;
        for(j=0;j<psPO->num;j++)
        {
            sid = recsPO[j].id;
            if(sid == recsPO[rep].id) {
                continue;
            }
/***
            nsc++;
            if((nsc%UPDATE_NUM)==0)
            {
                printf("#    %d %d reps %4.2f%%\n",nsc,rep,
                    PERCENT_R(rep,psPO->tie));
                fflush(stdout);
            }
***/
            scR = PickSeqScoreR(psPO,sub,sid);
            psPO->tscores[j] = scR;
            maxR = MAX_NUM(scR,maxR);
        }
        DB_PS DB_PrI("+ try%d s=%d max=%1.3f",try,sub,maxR);
        /***
        *   If sub s is not better than worst i, no replace 
        */
        if(maxR >= psPO->worst) {
            DB_PS DB_PrI(" Not better (%d not accepted)\n",sub);
            continue;
        }
        /***
        *   Replace rep with sub 
        */
        oid = recsPO[rep].id;
        ssPO->mask[oid] = 0;
        recsPO[rep].id = sub;
        recsPO[rep].sc = maxR;
        ssPO->mask[sub] = 1;
        DB_PS DB_PrI(" REPLACE %d with %d\n",oid,sub);
        /***
        *   Replace all scores in current pool that have trial scores 
        *   greater than the pre-sub values
        */
        for(j=0;j<psPO->num;j++)
        {
            /***
            *   rep is now sub, the new guy with score set above
            */
            if(j==rep) {
                continue;
            }
            if(psPO->tscores[j] > recsPO[j].sc) {
                recsPO[j].sc = psPO->tscores[j];
            }
        }
        nrep++;
/**
???
*/
        ravR += maxR;
        break;
    }
    return(nrep);
}
/****************************************************************************
*   Find the best and worst scores in the current selected set
*/
int FindChosenBestWorstStatsI(PICKSEQ *psPO, int cyc)
{
    int i,n,tie,btie;
    REAL bestR, worstR, rworstR, avR;
    SCOREC *recsPO;

    worstR = rworstR = -TOO_BIG_R;
    bestR = TOO_BIG_R;
    n = 0;
    recsPO = psPO->recs;
    for(i=0;i<psPO->num;i++)
    {
        /***
        *   Fixed guys don't count here
        */
        rworstR = MAX_NUM(recsPO[i].sc,rworstR);
        if( (psPO->fixmap) && (psPO->fixmap[recsPO[i].id]) ) {
            continue;
        }
        n++;
        worstR = MAX_NUM(recsPO[i].sc,worstR);
        bestR = MIN_NUM(recsPO[i].sc,bestR);
    }
    /***
    *   Count ties and report the story
    */
    tie = btie = 0;
    avR = 0.0;
    for(i=0; i<psPO->num; i++)
    {
        DB_PS DB_PrI("+ Rec[%d] %d %1.3f\n", i, recsPO[i].id,recsPO[i].sc);
        if(recsPO[i].sc == worstR) {   
            tie++;  
        }
        if(recsPO[i].sc == bestR) {   
            btie++; 
        }
        avR += recsPO[i].sc;
    }
    printf("\n# Cycle %d\n",cyc);
    printf("#  Average score %8.4f\n",avR/RNUM(psPO->num));
    printf("#  Best score    %8.4f (%d)\n",bestR,btie);
    printf("#  Worst score   %8.4f (%d)\n",worstR,tie);
    fflush(stdout);
    /***
    *   Set values, return the number of real things considered
    */
    psPO->best = bestR;
    psPO->worst = worstR;
    psPO->rworst = rworstR;
    psPO->ave = avR;
    psPO->tie = tie;
    psPO->btie = btie;
    return(n);
}
/****************************************************************************
*   Set initial set
*/
int SetStartingSetI(PICKSEQ *psPO, char *fixmapPC, char *pickmapPC)
{
    int i,n,fdir,sdir;
    SCOREC *recsPO;

    /***
    *   Get random set 
    */
    MaskRandSubsetI(pickmapPC,psPO->nseq,RNUM(psPO->num)/RNUM(psPO->nseq));
/*
DumpArray(fixmapPC,psPO->nseq,"fix   %d\n",IS_CHAR,NULL);
DumpArray(pickmapPC,psPO->nseq,"rand %d\n",IS_CHAR,NULL);
*/
    /***
    *   If we have a fixed subset, mark these guys
    */
    if(fixmapPC) {
        /***
        *   Scan fixed set, looking for any not in set
        */
        for(i=0;i<psPO->nseq;i++)
        {
            if(pickmapPC[i]) {
                continue;
            }
            /***
            *   Scan mask for first non-fixed one to swap; 
            *   Randomly start in either direction
            */
            if(RandR(1.0)>0.5) {
                fdir = FORWARD;
                sdir = REVERSE;
            }
            else {
                fdir = REVERSE;
                sdir = FORWARD;
            }
            if(RandFixedSwapI(fixmapPC,pickmapPC,psPO->nseq,i,fdir)) {
                continue;
            }
            if(RandFixedSwapI(fixmapPC,pickmapPC,psPO->nseq,i,sdir)) {
                continue;
            }
            BOG_CHECK(TRUE);
        }
    }
/*
DumpArray(pickmapPC,psPO->nseq,"done %d\n",IS_CHAR,NULL);
*/
    /***
    *   Copy indices for masked subset into score-record datatypes
    */
    recsPO = psPO->recs;
    n = 0;
    for(i=0;i<psPO->nseq;i++)
    {   
        if(pickmapPC[i]) {
            DB_PS DB_PrI("+ pick[%d] = %d\n",n,i);
            recsPO[n].id = i;   
            recsPO[n].sc = -TOO_BIG_R;  
            n++;
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Set initial scores
*/
int SetInitialInternalScoresI(PICKSEQ *psPO)
{
    int i,j,fid,sid,ntot,nsc;
    REAL scR;
    SCOREC *recsPO;

    recsPO = psPO->recs;
    /***
    *   Set initial score for everybody in starting set
    *   Symmetric pairwise scoring, so only one triangle of matrix needed
    *   For each pair, the worst interaction (highest score) is set
    */
    ntot = psPO->num * (psPO->num -1) / 2;
    printf("# Doing initial pair-wise comparision (%d X %d = %d)\n",
        psPO->num, psPO->num, ntot);
    nsc = 0;
    for(i=0;i<(psPO->num-1);i++)
    {
        fid = recsPO[i].id;
        for(j=i+1;j<psPO->num;j++)
        {
            sid = recsPO[j].id;
            nsc++; 
            if((nsc%UPDATE_NUM)==0) {
                printf("#    %d %5.2f%%\n",nsc,PERCENT_R(nsc,ntot));
                fflush(stdout);
            }
            scR = PickSeqScoreR(psPO,fid,sid);
            recsPO[i].sc = MAX_NUM(scR,recsPO[i].sc);
            recsPO[j].sc = MAX_NUM(scR,recsPO[j].sc);
        }
    }
    printf("\n");
    return(TRUE);
}
/******************************************************************************
*   Try to find a non-fixed masked element to swap, moving in the indicated 
*   direction down the collection from the fixed element f
*/
int RandFixedSwapI(char *fixmapPC, char *pickmapPC, int n,int f,int dir)
{
    int j;

    if(dir==REVERSE) {
        for(j=f-1;j>=0;j--)
        {
            if( (!fixmapPC[j]) && (pickmapPC[j]) ) {
                pickmapPC[j] = FALSE;
                pickmapPC[f] = TRUE;
                return(TRUE);
            }
        }
    }
    else {
        for(j=f+1;j<n;j++)
        {
            if( (!fixmapPC[j]) && (pickmapPC[j]) ) {
                pickmapPC[j] = FALSE;
                pickmapPC[f] = TRUE;
                return(TRUE);
            }
        }
    }
    return(FALSE);
}
/****************************************************************************
*   Compute a sequence score 
*/
REAL PickSeqScoreR(PICKSEQ *psPO,int fid,int sid)
{
    char *fPC,*sPC;
    int flen,slen;
    REAL sR,scR;

    DB_PS DB_PrI(">> PickSeqScoreR %d %d\n",fid,sid);
    GetSeqsetSeqStringI(psPO->seqs,fid,&fPC);
    GetSeqsetSeqStringI(psPO->seqs,sid,&sPC);
    flen = GetSeqsetSeqLenI(psPO->seqs,fid);
    slen = GetSeqsetSeqLenI(psPO->seqs,sid);
    scR = 0.0;
    /***
    *   Tally score from different comparisons
    *   Total match 
    */
    if(psPO->smat > 0.0) {
        /***
        *   Similarity total match
        */
        SetPair_parIdentMatch(psPO->pp);
        SetPparsCtypeI(psPO->pp,PP_COM_SIM);
        SetPparsAlgI(psPO->pp,ALGO_MATCH);
        ScoreSeqCompareI(fPC,flen,sPC,slen,psPO->pp,&sR);
        scR += (sR * psPO->smat);
    }
    if(psPO->cmat > 0.0) {
        /***
        *   Compliment total match
        */
        SetPair_parCompMatch(psPO->pp);
        SetPparsCtypeI(psPO->pp,PP_COM_COM);
        SetPparsAlgI(psPO->pp,ALGO_MATCH);
        ScoreSeqCompareI(fPC,flen,sPC,slen,psPO->pp,&sR);
        scR += (sR * psPO->cmat);
    }
    /***
    *   Contig match
    */
    if(psPO->scon > 0.0) {
        /***
        *   Similarity contig match
        */
        SetPair_parIdentMatch(psPO->pp);
        SetPparsCtypeI(psPO->pp,PP_COM_SIM);
        SetPparsAlgI(psPO->pp,ALGO_CONT);
        ScoreSeqCompareI(fPC,flen,sPC,slen,psPO->pp,&sR);
        scR += (sR * psPO->scon);
    }
    if(psPO->ccon > 0.0) {
        /***
        *   Compliment contig match
        */
        SetPair_parCompMatch(psPO->pp);
        SetPparsCtypeI(psPO->pp,PP_COM_COM);
        SetPparsAlgI(psPO->pp,ALGO_CONT);
        ScoreSeqCompareI(fPC,flen,sPC,slen,psPO->pp,&sR);
        scR += (sR * psPO->ccon);
    }
    /***
    *   Block match
    */
    if(psPO->scb > 0.0) {
        /***
        *   Similarity block-weight match
        */
        SetPair_parIdentMatch(psPO->pp);
        SetPparsCtypeI(psPO->pp,PP_COM_SIM);
        SetPparsAlgI(psPO->pp,ALGO_BMATCH);
        ScoreSeqCompareI(fPC,flen,sPC,slen,psPO->pp,&sR);
        scR += (sR * psPO->scb);
    }
    if(psPO->ccb > 0.0) {
        /***
        *   Compliment block-weighted match
        */
        SetPair_parCompMatch(psPO->pp);
        SetPparsCtypeI(psPO->pp,PP_COM_COM);
        SetPparsAlgI(psPO->pp,ALGO_BMATCH);
        ScoreSeqCompareI(fPC,flen,sPC,slen,psPO->pp,&sR);
        scR += (sR * psPO->ccb);
    }
    return(scR);
}
/****************************************************************************
*   Save temp results (i.e. seqs) to file
*/
int HandleTempResultsI(PICKSEQ *psPO,int n)
{
    int i,id;
    char bufS[DEF_BS],outS[NSIZE];
    FILE *outPF;

    if(NO_S(psPO->outname)) {
        GetFilePartsI(psPO->seqname,NULL,bufS,NULL);
    }
    else {
        GetFilePartsI(psPO->outname,NULL,bufS,NULL);
    }
    if(psPO->siru) {
        sprintf(outS,"%s-pick_seq-%03d.temp",bufS,n);
    }
    else {
        sprintf(outS,"%s-pick_seq.temp",bufS);
    }
    if(!(outPF = OpenUFilePF(outS,"w",NULL))) {
        printf("Can't write temp results!!!\n");
        return(FALSE);
    }
    /***
    *   Normal header, disclaimer and seqs
    */
    WritePickSeqHeader(psPO,outPF);
    fflush(outPF);
    fprintf(outPF,"#\n");
    fprintf(outPF,"# INCOMPLETE RUN\n"); 
    fprintf(outPF,"# ONLY %d OF %d CYCLES EXECUTED\n",n,psPO->cyc); 
    fprintf(outPF,"#\n");
    fprintf(outPF,"# Average score %8.4f\n",psPO->ave/RNUM(psPO->num));
    fprintf(outPF,"# Best score    %8.4f (%d)\n",psPO->best,psPO->btie);
    fprintf(outPF,"# Worst score   %8.4f (%d)\n",psPO->worst,psPO->tie);
    fprintf(outPF,"#\n");
    /***
    *   Dump seqs and scores
    */
    for(i=0;i<psPO->num;i++)
    {
        id = psPO->recs[i].id;
        FillSeqsetSeqNameI(psPO->seqs,id,outS,20);
        FillSeqsetSeqStringI(psPO->seqs,id,bufS,DEF_BS);
        fprintf(outPF,"%-12s %-30s %5.3f\n",outS,bufS,psPO->recs[i].sc);
    }
    CHECK_NFILE(outPF,outS);
    return(TRUE);   
}
/****************************************************************************
*   Allocate pickseq data structure
*/
PICKSEQ *CreatePickseqPO()
{
    PICKSEQ *psPO;

    if(!(psPO = (PICKSEQ *)ALLOC(1,sizeof(PICKSEQ)))) {
        return(NULL);
    }
    psPO->ID = PICKSEQ_ID;
    if(!InitPickseqI(psPO)) {
        CHECK_PICKSEQ(psPO);
        return(NULL);
    }
    return(psPO);
}
/****************************************************************************
*   Free pickseq data structure and any associated sub-structs
*/
int DestroyPickseqI(PICKSEQ *psPO)
{
    VALIDATE(psPO,PICKSEQ_ID);
    CHECK_PPARS(psPO->pp);
    CHECK_SEQSET(psPO->seqs);
    CHECK_SCOREC(psPO->recs);
    CHECK_FREE(psPO->tmask);
    CHECK_FREE(psPO->tscores);
    CHECK_FREE(psPO->bias);
    CHECK_FREE(psPO->fixmap);
    CHECK_NFILE(psPO->out,psPO->outname);
    FREE(psPO);
    return(TRUE);
}
/****************************************************************************
*   Initialize pickseq data struct fields
*/
int InitPickseqI(PICKSEQ *psPO)
{
    VALIDATE(psPO,PICKSEQ_ID);
    INIT_S(psPO->seqname);
    psPO->iform = BOGUS;
    psPO->seqs = NULL;
    psPO->nseq = 0;
    INIT_S(psPO->mwfname);
    INIT_S(psPO->biasname);
    INIT_S(psPO->fixname);
    psPO->alg = DEF_PS_ALG;
    psPO->num = DEF_PS_NUM;
    psPO->rseed = BAD_I;
    psPO->cyc = DEF_PS_CYC;
    psPO->rep = DEF_PS_REP;
    psPO->smat = DEF_PS_SMAT;
    psPO->cmat = DEF_PS_CMAT;
    psPO->scon = DEF_PS_SCON;
    psPO->ccon = DEF_PS_CCON;
    psPO->scb = DEF_PS_SCB;
    psPO->ccb = DEF_PS_CCB;
    psPO->do_mswm = FALSE;
    psPO->do_cswm = FALSE;
    psPO->do_ham = FALSE;
    if( !(psPO->pp = CreatePparsPO())) {
        return(FALSE);
    }
    psPO->tmask = NULL;
    psPO->tscores = NULL;
    psPO->bias = NULL;
    psPO->fixmap = NULL;
    psPO->nfix = 0;
    psPO->verbose = TRUE;
    psPO->do_bwn = -1;
    psPO->do_bwf = -1.0;
    return(TRUE);
}
/****************************************************************************
*
*/
void WritePickSeqHeader(PICKSEQ *psPO,FILE *outPF)
{
    char bufS[DEF_BS];

    HAND_NFILE(outPF);
    fprintf(outPF,"# %s\n",VERSION_S);
    fprintf(outPF,"#     %s %s %s\n",BD_S,__DATE__,__TIME__);
    fprintf(outPF,"#     %s\n",RTK_S);
    fprintf(outPF,"# Input file  %s (%d seqs)\n",psPO->seqname,psPO->nseq);
    TimeStamp("# Date        ",outPF);
    FillPsAlgoStringI(psPO->alg,bufS);
    fprintf(outPF,"# Algorithm   %s\n",bufS);
    if(BAD_INT(psPO->rseed)) {
        fprintf(outPF,"# Random seed (clock)\n");
    }
    else {
        fprintf(outPF,"# Random seed %d\n",psPO->rseed);
    }
    if(!NO_S(psPO->fixname)) {
        fprintf(outPF,"# Fixed subset: %s\n",psPO->fixname);
        fprintf(outPF,"# Number fixed: %d\n",psPO->nfix);
    }
    fprintf(outPF,"#\n");
    fprintf(outPF,"# Sim Matches     %4.2f\n",psPO->smat);
    fprintf(outPF,"# Com Matches     %4.2f\n",psPO->cmat);
    if(psPO->do_mswm) {
        if(!NO_S(psPO->mwfname)) {
            fprintf(outPF,"#     Matches     weighted: %s\n",psPO->mwfname);
        }
        else {
            fprintf(outPF,"#     Matches     weighted (default)\n");
        }
    }
    else {
        fprintf(outPF,"#     Matches     simple count\n");
    }
    fprintf(outPF,"# Sim Contig      %4.2f\n",psPO->scon);
    fprintf(outPF,"# Com Contig      %4.2f\n",psPO->ccon);
    if(psPO->do_cswm) {
        if(!NO_S(psPO->mwfname)) {
            fprintf(outPF,"#     Contig      weighted: %s\n",psPO->mwfname);
        }
        else {
            fprintf(outPF,"#     Contig      weighted (default)\n");
        }
    }
    else {
        fprintf(outPF,"#     Contig      simple count\n");
    }
    /***
    *   Block match
    */
    fprintf(outPF,"# Sim Blk-Wt-Mat  %4.2f\n",psPO->scb);
    fprintf(outPF,"# Com Blk-Wt-Mat  %4.2f\n",psPO->ccb);
    /***
    *   Weight so dump weights (have to turn on here as well)
    */
    if( (psPO->do_mswm) || (psPO->do_cswm) ) {
        fprintf(outPF,"# Base pair weights (for simple and contig matches)\n");
        DumpPparsWeights(psPO->pp,outPF);
    }
    /***
    *   Bias file??? = SHAM?
    */
    if(psPO->bias) {
        fprintf(outPF,"# Bias file   %s\n",psPO->biasname);
    }
    fprintf(outPF,"#\n");
    fflush(outPF);
}
/****************************************************************************
*
*/
int FillPsAlgoStringI(int alg,char *bufS)
{
    switch(alg)
    {
        case PS_ALG_SWAP: 
            sprintf(bufS,"Monte Carlo Swapping");
            return(TRUE);
    }
    return(FALSE);
}
