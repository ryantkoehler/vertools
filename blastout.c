/*
* blastout.c
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
#define __MAIN__
#include "prim.h"
#include "blastout.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(BlastOutI(argc,argv),NULL) ); }
/*******************************************************************/
void BlastOutUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <blastfile> [...options]\n");
    printf("   <blastfile>  Blast search log file / Blastout dci file\n");
    printf("   -merg XXX    Merge together XXX and <blastlog>\n");
    printf("   -out XXX     Set output file to XXX\n");
    printf("   -dsum        Dump summary of N-hits and Max-hit per Query\n");
    printf("   -dci         Dump Condensed Info; hits, seqs, coords, scores\n");
    printf("   -dhit        Dump hit summaries\n");
    printf("   -dsco        Dump hit Scores\n");
    printf("   -dseq        Dump hit Sequences\n");
    printf("   -nsqh XXX    Name output seqs with query + hit separated by XXX\n");
    printf("   -smc -sml    Mark output seqs via Case: Matches = upper / lower\n");
    printf("   -dhc         Dump hit coordinates\n");
    printf("   -dfml #      Dump Fraction Match list for top # hits\n");
    printf("   -dffl        Dump fraction match full-length (not any base range limit)\n");
    printf("   -dhis #      Dump hit match-base histogram starting at #\n");
    printf("   -phis #      Pad (with 0) dumped match histogram up to #\n");
    printf("   -chis        Cumulative hit counts in histogram\n");
    printf("   -dmbc        Dump matching base counts (along query seqs)\n");
    printf("   -mid #       Set minimum match identities to #\n");
    printf("   -mif #       Set minimum match fraction, i.e. 12/16 = 0.75\n");
    printf("   -mfq         Min Fraction of Query (not of hit)\n");
    printf("   -mnot        Inverse (i.e. 'not') for minimum match qualifications\n");
    printf("   -con         Count contiguous matches, not total\n");
    printf("   -co3         Count contiguous matches starting from 3' end\n");
    printf("   -bran # #    Consider query bases only in range # to #\n");
    printf("   -rre         Base range relative to end; i.e. backwards\n");
    printf("   -qran # #    Restrict output to query range # to #\n"); 
    printf("   -hran # #    Restrict output to hits range # to #\n"); 
    printf("   -opq XXX     Output file per query with extension XXX\n");
    printf("   -mhit #      Set maximum number of hits to # (def = %d)\n",
        DEF_MAXHITS);
}
/**************************************************************************
*   top level function
*/
int BlastOutI(int argc, char **argv)
{
    char inS[DEF_BS],mergS[DEF_BS];
    BLASTOUT *bPO;
    int n,r1,r2,mhit,warn;

    warn = TRUE;
    mhit = DEF_MAXHITS;
    INIT_S(mergS);
    bPO = CreateBlastoutPO();
    if(!ParseArgsI(argc, argv,
        "S -out S -dhit B -mid I -mif R -dsum B -dhis I -phis I\
        -dseq B -con B -dhc B -dmbc B -bran I2 -rre B -dsco B -dci B\
        -qran I2 -opq S -chis B -mhit I -merg S -hran I2 -mfq B\
        -smc B -dfml I -mnot B -sml B -nsqh S -dffl B -co3 B",
        inS, bPO->outname, &bPO->dhit, &bPO->mid, &bPO->mif, 
        &bPO->dsum, &bPO->dhis, &bPO->phis, &bPO->dseq, 
        &bPO->do_con, &bPO->dhc, &bPO->dmbc, &bPO->firstb,&bPO->lastb, 
        &bPO->rre, &bPO->dsco, &bPO->do_dci, &bPO->firstq,&bPO->lastq, 
        &bPO->opq, &bPO->chis, &mhit, mergS, &bPO->firsth,&bPO->lasth,
        &bPO->do_mfq, &bPO->do_smc, &bPO->do_dfml, &bPO->do_mnot,
        &bPO->do_sml, bPO->nsqh, &bPO->do_dffl, &bPO->do_co3,
        (int *)NULL))
    {
        BlastOutUse();
        CHECK_BLASTOUT(bPO);
        return(FALSE);
    }
    /***
    *   Check options
    */
    if(!OkBlastoutOptsI(bPO,mhit)) {
        CHECK_BLASTOUT(bPO);
        ABORTLINE;
        return(FALSE);
    }
    /***
    *    Add space for anwsers & open input(s)
    */
    if(!AddHitSpaceI(mhit,&bPO->ans)) {
        CHECK_BLASTOUT(bPO);
        ABORTLINE;
        return(FALSE);
    }
    if(!SetBlastansFileI(inS,bPO->ans)) {
        CHECK_BLASTOUT(bPO);
        ABORTLINE;
        return(FALSE);
    }
    if(!NO_S(mergS)) {
        if(!AddHitSpaceI(mhit,&bPO->sans)) {
            CHECK_BLASTOUT(bPO);
            ABORTLINE;
            return(FALSE);
        }
        if(!SetBlastansFileI(mergS,bPO->sans)) {
            CHECK_BLASTOUT(bPO);
            ABORTLINE;
            return(FALSE);
        }
        if(!AddHitSpaceI(mhit,&bPO->mans)) {
            CHECK_BLASTOUT(bPO);
            ABORTLINE;
            return(FALSE);
        }
    }
    /***
    *   output 
    */
    if(!HandleBlastoutOutfileI(bPO,bPO->ans,0,TRUE)) {
        printf("Unrecognized input format\n");
        CHECK_BLASTOUT(bPO);
        ABORTLINE;
        return(FALSE);
    }
    WriteHeader(bPO,bPO->ans,bPO->sans,bPO->out);
    /***
    *   Party through log file(s)
    */
    n = 0;
    if(!NO_S(mergS)) {
        while(TRUE)
        {
            r1 = LoadBlastoutRecordI(bPO->ans,warn);
            r2 = LoadBlastoutRecordI(bPO->sans,warn);
            if( (!r1) || (!r2) ) {
                break;
            }
            if(!MergeHitListsI(bPO,bPO->ans,bPO->sans,bPO->mans)) {
                printf("Problem with records %d\n",n+1);
                ABORTLINE;
                break;
            }
            n++;
            HandleBlastoutOutfileI(bPO,bPO->mans,n,FALSE);
            DumpBlastout(bPO,bPO->mans,n,bPO->out);
        }
    }
    else {
        while(LoadBlastoutRecordI(bPO->ans,warn))
        {
            n++;
            if( (n<bPO->firstq) || (n>bPO->lastq) ) {
                continue;
            }
            if( (bPO->owhat==BOT_PERQ) || (bPO->owhat==BOT_PERHIT) ) {
                ProcHitList(bPO,bPO->ans);
                HandleBlastoutOutfileI(bPO,bPO->ans,n,FALSE);
                DumpBlastout(bPO,bPO->ans,n,bPO->out);
            }
        }
    }
    /***
    *   All done
    */
    printf("# Total queries in log: %d\n",n);
    CHECK_NFILE(bPO->out,bPO->outname);
    CHECK_BLASTOUT(bPO);
    return(TRUE);
}
/****************************************************************************
*   Check / set consistent options
*/
int OkBlastoutOptsI(BLASTOUT *bPO,int max)
{
    if(max<1) {
        PROBLINE;
        printf("Bad max hits: %d\n",max);
        return(FALSE);
    }
    if( (!NO_S(bPO->outname)) && (!NO_S(bPO->opq)) ) {
        PROBLINE;
        printf("Can't have a named output AND dump per query\n");
        return(FALSE);
    }
    if( (bPO->mif > 1.0) || (bPO->mif < 0.0) ) {
        PROBLINE;
        printf("Minimum fraction should be in the range 0.0 to 1.0\n");
        printf("    %3.4f won't work\n",bPO->mif);
        return(FALSE);
    } 
    if(bPO->firstb > 0) {
        if( (bPO->dhis>0) || (bPO->dmbc) ) {
            PROBLINE;
            printf("Base range restrictions aren't allowed here:\n");
            printf("   -dhis\n");
            printf("   -dmbc\n");
            return(FALSE);
        }
        if(bPO->lastb < bPO->firstb) {
            PROBLINE;
            printf("Base range %d to %d won't work\n",bPO->firstb,bPO->lastb);
            return(FALSE);
        }
    }
    /***
    *   Condensed info = summ of other options
    */
    if(bPO->do_dci) {
        bPO->dhit = bPO->dseq = bPO->dhc = bPO->dsco = TRUE;
    }
    /***
    *   Per query or per hit?
    */
    if( (bPO->dsum ) || (bPO->dhis > 0) || (bPO->dmbc) || (bPO->do_dfml > 0) ) {
        bPO->owhat = BOT_PERQ;
        LIMIT_NUM(bPO->phis,0,(HISTDIM-1));
        bPO->dhit = bPO->dseq = bPO->dhc = bPO->dsco = FALSE;
    }
    else if( (bPO->dhit ) || (bPO->dseq) || (bPO->dhc) || (bPO->dsco) ) {
        bPO->owhat = BOT_PERHIT;
        bPO->dsum = bPO->dhis = bPO->dmbc = FALSE;
    }
    return(TRUE);
}
/****************************************************************************
*   Write header story reporting what's comming out of the blast log
*/
void WriteHeader(BLASTOUT *bPO,BLASTANS *aPO,BLASTANS *mPO,FILE *outPF)
{
    char nameS[DEF_BS];

    HAND_NFILE(outPF);
    fprintf(outPF,"# %s, %s %s %s\n",VERSION_S,BD_S,__DATE__,__TIME__);
    fprintf(outPF,"#   %s\n",RTK_S);
    TimeStamp("# ",outPF);
    if(aPO) {
        fprintf(outPF,"# Input %s\n",aPO->input);
        FillBlastFormatString(aPO->itype,nameS);
        fprintf(outPF,"#   Format %s\n",nameS);
    }
    if(mPO) {
        fprintf(outPF,"# MergeInput %s\n",mPO->input);
        FillBlastFormatString(mPO->itype,nameS);
        fprintf(outPF,"#   Format %s\n",nameS);
    }
    if(bPO->owhat == 0) {
        return;
    }
    fprintf(outPF,"#\n");
    /***
    *   Restrictions
    */
    if(bPO->do_co3) {
        fprintf(outPF,"# 3' end contiguous matches considered\n");
    }
    else if(bPO->do_con) {
        fprintf(outPF,"# Contiguous matches considered\n");
    }
    else {
        fprintf(outPF,"# Total matches considered\n");
    }
    if(bPO->firstb > 0) {
        fprintf(outPF,"# Query bases restricted %d to %d",
            bPO->firstb,bPO->lastb);
        if(bPO->rre) {
            fprintf(outPF," (from end)\n");
        }
        else {
            fprintf(outPF," (from start)\n");
        }
    }
    else {
        fprintf(outPF,"# Query bases unrestricted (all considered)\n");
    }
    if(bPO->mid > 0) {
        fprintf(outPF,"# Minimum match count %d\n",bPO->mid);
    }
    else {
        fprintf(outPF,"# Minimum match count None\n");
    }
    if(bPO->mif > 0.0) {
        fprintf(outPF,"# Minimum match fraction %4.3f\n",bPO->mif);
        if(bPO->do_mfq) {
            fprintf(outPF,"# Minimum fraction of Query length\n");
        }
        else {
            fprintf(outPF,"# Minimum fraction of Hit length\n");
        }
    }
    else {
        fprintf(outPF,"# Minimum match fraction None\n");
    }
    if(bPO->do_mnot){
        fprintf(outPF,"# Minimum match tests Inverted (i.e. NOT)\n");
    }
    /***
    *   What output 
    */
    if(bPO->owhat == BOT_PERHIT) {
        fprintf(outPF,"# Reports for each hit\n");
    }
    else {
        fprintf(outPF,"# Reports for each query\n");
    }
    if(bPO->dhis > 0) {
        fprintf(outPF,"# Histogram starting at %d",bPO->dhis);
        if(bPO->chis) {
            fprintf(outPF," (cumulative)");
        }
        fprintf(outPF,"\n");
    }
    if(bPO->dmbc) {
        fprintf(outPF,"# Query & Matched base counts along query sequence\n");
    }
    if(bPO->dhc) {
        fprintf(outPF,"# Hit coordinates: Query then Subject\n");
    }
    if(bPO->dhit) {
        fprintf(outPF,"# HitStat: <name> <match> <match'> <mfrac> <qfrac>\n");
        fprintf(outPF,"#  <name>   = name of hit\n");
        fprintf(outPF,"#  <match>  = match reported by blast\n");
        fprintf(outPF,"#  <match'> = match as restricted here\n");
        fprintf(outPF,"#  <mfrac>  = matching base fraction of hit length\n");
        fprintf(outPF,"#  <qfrac>  = matching base fraction of query length\n");
    }
    return;
}
