/*
* seqtweak.c
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
#include "dna.h"
#include "seqtweak.h"

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); AllDoneI(SeqTweakI(argc,argv),NULL); exit(0); }
/**************************************************************************/
void SeqTweakUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>   Sequence file\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -out XXX   Set output basename to XXX\n");
    printf("   -bran # #  Limit base range # to # for modifications\n");
    printf("   -rre       Base range relative to end; i.e. backwards\n");
    printf("   -mis #     Randomly put in # mismatches\n");
    printf("   -del #     Randomly put in # deletions\n");
    printf("   -ins #     Randomly put in # insertions\n");
    printf("   -ids #     Insert/deletion word size (i.e. how many bases)\n");
    printf("   -mds #     Set Minimum Disruption Separation to # bases\n");
    printf("   -amb       All mismatch bases for each mismatch (i.e. 3/seq)\n");
    printf("   -smm       \"Smart\" mismatch placement\n");
    printf("   -cpm       Close-pack mismatch placement\n");
    printf("   -sh -fsh   Shuffle sequence; Full shuffle changes *all* bases\n");
    printf("   -seed #    Set random seed to #\n");
    printf("   -bname XXX Append output names with base XXX\n");
    printf("   -nsim      Use 'simple' naming (i.e. old)\n");
    printf("   -nre       Name relative to end (i.e. 3' == 1)\n");
    printf("   -quiet     No verbose description of tweaking; only seqs\n");
}
/**************************************************************************/
int SeqTweakI(int argc,char **argv)
{
    int ok,iraw,iseq,ifas,verb,nseq,n;
    SEQTWEAK *stPO;

    stPO = CreateSeqtweakPO();
    verb = TRUE;
    iraw = iseq = ifas = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -iraw B -ifas B -mis I -del I -ins I -ids I -bran I2\
        -seed I -mds I -quiet B -amb B -smm B -cpm B -bname S\
        -sh B -fsh B -rre B -nsim B -nre B -iseq B",
        stPO->inname, stPO->outname, &iraw, &ifas, 
        &stPO->mis, &stPO->del, &stPO->ins, &stPO->ids, 
        &stPO->firstb,&stPO->lastb, &stPO->seed, 
        &stPO->mds, &verb, &stPO->do_amb, &stPO->do_smm, &stPO->do_cpm,
        &stPO->bname,
        &stPO->do_sh, &stPO->do_fsh, &stPO->do_rre,
        &stPO->do_nsim, &stPO->do_nre, &iseq,
        (int *)NULL))
    {
        SeqTweakUse();
        CHECK_SEQTWEAK(stPO);
        return(FALSE);
    }
    /***
    *   Set input format
    */
    stPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas, stPO->inname, TRUE);
    if(!stPO->iform) {
        printf("Problem with input seq(s)\n");
        CHECK_SEQTWEAK(stPO);
        return(FALSE);
    }
    /***
    *   Check options
    */
    if(!OkSeqTweakOptsI(stPO)) {
        CHECK_SEQTWEAK(stPO);
        ABORTLINE;
        return(FALSE);
    }
    /***
    *   Set up for use
    */
    if(!SetUpSeqtweakI(stPO)) {
        CHECK_SEQTWEAK(stPO);
        ABORTLINE;
        return(FALSE);
    }
    /***
    *   Header whip
    */
    if(verb) {
        WriteSeqtweakHeader(stPO,stPO->out);
    }
    /***
    *   Party through file; getting seqs and messing them up
    */
    nseq = n = 0;
    while(TRUE) {
        /***
        *   Parse sequence; FALSE = done
        */
        ok = ParseSeqI(stPO->in, stPO->iform, n+1, FALSE, TRUE, stPO->seq);
        if(ok==FALSE) {
            break;
        }
        if(ok!=TRUE) {
            continue;
        }
        nseq++;
        if(!TweakableSeqI(stPO,stPO->seq,nseq)) {
            continue;
        }
        if( stPO->do_sh || stPO->do_fsh ) {
            n += TweakShuffleSeqI(stPO,verb,stPO->out);
        }
        else {
            n += TweakSeqBasesI(stPO,verb,stPO->out);
        }
    }
    if(verb)
    {
        printf("# Input sequences:  %d\n",nseq);
        printf("# Output sequences: %d\n",n);
    }
    /***
    *   All done
    */
    CHECK_SEQTWEAK(stPO);
    return(1);
}
/**************************************************************************
*   Allocate data structure
*/
SEQTWEAK *CreateSeqtweakPO()
{
    SEQTWEAK *stPO;

    if(! (stPO = (SEQTWEAK *)ALLOC(1,sizeof(SEQTWEAK)) ) )
    {
        printf("# Failed to allocate working object\n");
        return(NULL);
    }
    stPO->ID = SEQTWEAK_ID;
    stPO->seq = CreateSeqPO(0,NULL,NULL);
    if(!stPO->seq)
    {
        printf("# Failed to allocate seq space\n");
        CHECK_SEQTWEAK(stPO);
        return(NULL);
    }
    InitSeqtweak(stPO);
    return(stPO);
}
/**************************************************************************
*   Free data structrue and substructures
*/
int DestroySeqtweakI(SEQTWEAK *stPO)
{
    VALIDATE(stPO,SEQTWEAK_ID);
    CHECK_FILE(stPO->in);
    CHECK_NFILE(stPO->out,stPO->outname);
    CHECK_SEQ(stPO->seq);
    FREE(stPO);
    return(TRUE);
}
/**************************************************************************
*   Initialize data structure
*/
void InitSeqtweak(SEQTWEAK *stPO)
{
    VALIDATE(stPO,SEQTWEAK_ID);
    INIT_S(stPO->inname);
    INIT_S(stPO->outname);
    stPO->in = NULL;
    stPO->iform = BOGUS;
    stPO->out = NULL;
    stPO->firstb = stPO->lastb = BOGUS;
    stPO->do_rre = FALSE;
    stPO->mis = 0;
    stPO->ins = 0;
    stPO->del = 0;
    stPO->ids = 1;
    stPO->mds = 0;
    stPO->do_amb = FALSE;
    stPO->do_smm = FALSE;
    stPO->do_cpm = FALSE;
    stPO->do_sh = stPO->do_fsh = FALSE;
    stPO->seed = BOGUS;
    INIT_S(stPO->bname);
    stPO->do_nre = FALSE;
    stPO->do_nsim = FALSE;
}
/**************************************************************************
*   Check options
*/
int OkSeqTweakOptsI(SEQTWEAK *stPO)
{
    if( (stPO->mds>0) && ( (stPO->ins>0)||(stPO->del>0) ) ) {
        PROBLINE;
        printf(" -mds option doesn't work for In/Dels\n");
        printf("\n");
        return(FALSE);
    }
    if( (stPO->del + stPO->ins) > 0) {
        printf("# In/Del so setting to shotgun placement\n");
        stPO->do_smm = stPO->do_cpm = FALSE;
    }
    return(TRUE);
}
/**************************************************************************
*   Open files / init
*/
int SetUpSeqtweakI(SEQTWEAK *stPO)
{
    /***
    *   Input output 
    */
    if(!(stPO->in=OpenUFilePF(stPO->inname,"r",NULL))) {
        return(FALSE);
    }
    if(!NO_S(stPO->outname)) {
        if(!(stPO->out=OpenUFilePF(stPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    HAND_NFILE(stPO->out);
    Srand(stPO->seed);
    /***
    *   Only single mismatch, else simple (old) style naming
    */
    if( (stPO->mis > 1) || (stPO->ins > 0) || (stPO->del > 0) ) {
        stPO->do_nsim = TRUE;
    }
    return(TRUE);
}
/**************************************************************************
*   Dump settings
*/
void WriteSeqtweakHeader(SEQTWEAK *stPO,FILE *outPF)
{
    char sS[DEF_BS];

    HAND_NFILE(outPF);
    fprintf(outPF,"# %s, Build date %s %s\n",VERSION_S,__DATE__,__TIME__);
    fprintf(outPF,"# %s\n",RTK_S);
    TimeStamp("# ",outPF);
    fprintf(outPF,"#\n");
    fprintf(outPF,"# Input: %s\n",stPO->inname);
    FillRandSeedString(stPO->seed,sS);
    fprintf(outPF,"# Random seed: %s\n",sS);

    fprintf(outPF,"# Mismatchs: %d\n",stPO->mis);
    fprintf(outPF,"# Insertions: %d (%d bases)\n",stPO->ins,stPO->ids);
    fprintf(outPF,"# Deletions:  %d (%d bases)\n",stPO->del,stPO->ids);
    fprintf(outPF,"# Minimum disruption separation: %d\n",stPO->mds);
    if( ! IS_BOG(stPO->firstb) ) {
        fprintf(outPF,"# Base range: %d %d\n",stPO->firstb,stPO->lastb);
    }
    if(stPO->do_amb) {
        fprintf(outPF,"# All (3) mismatch bases / mismatch\n");
    }
    if(stPO->do_smm) {
        fprintf(outPF,"# \"Smart\" disruption placement\n");
    }
    else if(stPO->do_cpm) {
        fprintf(outPF,"# Close-pack disruption placement\n");
    }
    else {
        fprintf(outPF,"# Shotgun disruption placement\n");
    }
    fprintf(outPF,"#\n");
    fprintf(outPF,"\n");
}
/***************************************************************************
*   Ok to mess with?
*/
int TweakableSeqI(SEQTWEAK *stPO, SEQ *seqPO, int nseq)
{
    if(seqPO->len > MAX_TWEAKSEQ) {
        printf("# Sequence %d too long: %d\n",nseq,seqPO->len);
        printf("#   Max = %d\n",MAX_TWEAKSEQ);
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************
*   Too big and complicated!
*/
int TweakSeqBasesI(SEQTWEAK *stPO, int verb, FILE *outPF)
{
    int i,j,k,slen,tlen,tmod,nmod,trys,ok,firstb,lastb;
    char nameS[NSIZE],nnameS[NSIZE];
    char twkS[DEF_BS],bufS[DEF_BS],maskS[MAX_TWEAKSEQ];
    char seqS[MAX_TWEAKSEQ], sseqS[MAX_TWEAKSEQ], tseqS[MAX_TWEAKSEQ];
    char bC, *seqPC;
    SEQ *seqPO;

    HAND_NFILE(outPF);
    if( (stPO->mis<1) && (stPO->del<1) && (stPO->ins<1) ) {
        return(FALSE);
    }
    tmod = stPO->mis + (stPO->del * stPO->ids) + (stPO->ins * stPO->ids);
    nmod = stPO->mis + stPO->del + stPO->ins;
    /***
    *   Get starting seq & info
    */
    seqPO = stPO->seq;
    if(!GetSeqSeqI(seqPO,&seqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(seqPO);
    FillSeqNameStringI(seqPO,nameS,NSIZE);
    /***
    *   Bound region to be modified and check things fit
    */
    if( ! IS_BOG(stPO->firstb) ) {
        firstb = (stPO->do_rre) ? (slen - stPO->lastb + 1)  : stPO->firstb;
        lastb = (stPO->do_rre)  ? (slen - stPO->firstb + 1) : stPO->lastb;
    }
    else {
        firstb = 1;
        lastb = slen;
    }
    LIMIT_NUM(firstb,1,slen);
    LIMIT_NUM(lastb,0,slen);
    tlen = lastb - firstb + 1;
    if( (tlen - ((nmod-1) * stPO->mds) ) < tmod) {
        printf("# Too many modifications %d\n",tmod);
        printf("#   Can't work for length range %d\n",tlen);
        return(FALSE);
    }
    if(stPO->ids>MAX_IDSIZE) {
        printf("# In-del size %d too big (max %d)\n",stPO->ids,MAX_IDSIZE);
        return(FALSE);
    }
    /***
    *   Try to find a tweak pattern that fits in allowed range
    */
    ok = trys = 0;
    while(trys<MAX_TWK_STARTS) 
    {
        if(FillTweakMaskI(stPO, tlen, firstb-1, maskS, slen)) {
            ok++;
            break;
        }
        trys++;
    }
    if(!ok) {
        printf("# Failed to place in-del / mismatch positions\n");
        return(FALSE);
    }
    /***
    *   For single mismatch naming..
    */
    INIT_S(stPO->smm_sub);
    stPO->smm_pos = BOGUS;
    /***
    *   Do the subs indicated in mask
    */
    i = j = k = 0;
    for(i=0;i<slen;i++)
    {
        switch(maskS[i])
        {
            /* 
            *   Insertion 
            */
            case 'I':   
                seqS[k++] = DNAIndexBaseC(RandI(4));
                break;
            /* 
            *   Deletion
            */
            case 'D':
                j++;
                break;
            /* 
            *   Miss-match 
            */
            case 'M':
                if(stPO->do_amb)
                {
                    switch(seqPC[j])
                    {
                        case 'a': case 'A':
                            seqS[k] = 'C';
                            sseqS[k] = 'G';
                            tseqS[k] = 'T';
                            break;
                        case 'c': case 'C':
                            seqS[k] = 'A';
                            sseqS[k] = 'G';
                            tseqS[k] = 'T';
                            break;
                        case 'g': case 'G':
                            seqS[k] = 'A';
                            sseqS[k] = 'C';
                            tseqS[k] = 'T';
                            break;
                        case 't': case 'T':
                            seqS[k] = 'A';
                            sseqS[k] = 'C';
                            tseqS[k] = 'G';
                            break;
                    }
                    k++;
                }
                else {
                    bC = DNAIndexBaseC(RandI(4));
                    while(UPPER(bC)==UPPER(seqPC[j]))
                    {
                        bC = DNAIndexBaseC(RandI(4));
                    }
                    /* aving naming stuff */
                    stPO->smm_sub[0] = seqPC[j];
                    stPO->smm_sub[1] = ':';
                    stPO->smm_sub[2] = bC;
                    stPO->smm_pos = k;
                    /* ake sub and advance */
                    seqS[k] = bC;
                    k++;
                }
                j++;
                break;
            default:
                seqS[k] = sseqS[k] = tseqS[k] = seqPC[j++];
                k++;
                break;
        }
    }
    while(k<slen)
    {
        if(j<slen) {
            seqS[k++] = seqPC[j++];
        }
        else {
            seqS[k++] = DNAIndexBaseC(RandI(4));
        }
    }
    seqS[slen] = sseqS[slen] = tseqS[slen] = '\0';
    /***
    *   Report the disruptions and modify name
    */
    INIT_S(twkS);
    trys = 0;
    for(i=0;i<slen;i++)
    {
        if(maskS[i]=='.') {
            continue;
        }
        if( (maskS[i]=='I') || (maskS[i]=='D') ) {
            sprintf(bufS,"%c%d=%d",maskS[i],stPO->ids,i+1);
            i+= (stPO->ids-1);
        }
        else {
            sprintf(bufS,"%c=%d",maskS[i],i+1);
        }
        if(trys>0) {
            strcat(twkS," ");
        }
        strcat(twkS,bufS);
        trys++;
    }
    maskS[slen] = '\0';
    /***
    *   Tell the story
    */
    SetTweakedName(stPO,nameS,nnameS,slen);
    if(verb) {
        fprintf(outPF,"# Disrupts %d %s\n",(stPO->ins+stPO->del+stPO->mis),twkS);
        sprintf(bufS,"#%s",nameS);
        PadString(bufS,' ',strlen(nnameS));
        bufS[strlen(nnameS)] = '\0';
        fprintf(outPF,"%s\t%s\n",bufS,seqPC);
        sprintf(bufS,"#");
        PadString(bufS,' ',strlen(nnameS));
        bufS[strlen(nnameS)] = '\0';
        fprintf(outPF,"%s\t%s\n",bufS,maskS);
    }
    /***
    *   Special output names for All-Missmatch-Base case
    */
    if(stPO->do_amb) {
        fprintf(outPF,"%s_x\t%s\n",nameS,seqS);
        fprintf(outPF,"%s_y\t%s\n",nameS,sseqS);
        fprintf(outPF,"%s_z\t%s\n",nameS,tseqS);
    }
    else {
        fprintf(outPF,"%s \t%s\n",nnameS,seqS);
    }
    if(verb) {
        fprintf(outPF,"\n");
    }
    return(TRUE);
}
/**************************************************************************
*   Shuffle sequence 
*/
int TweakShuffleSeqI(SEQTWEAK *stPO,int verb, FILE *outPF)
{
    int shuffIA[MAX_TWEAKSEQ];
    int i, firstb, lastb, slen, mlen;
    char *seqPC, nameS[NSIZE], newS[MAX_TWEAKSEQ];
    SEQ *seqPO;

    HAND_NFILE(outPF);
    INIT_S(newS);
    if(stPO->do_fsh) {
        printf("Full shuffle not yet implemented\n");
        return(FALSE);
    }
    /***
    *   Get starting seq & info
    */
    seqPO = stPO->seq;
    if(!GetSeqSeqI(seqPO,&seqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(seqPO);
    FillSeqNameStringI(seqPO,nameS,NSIZE);
    if(slen > MAX_TWEAKSEQ) {
        printf("Sequence %s is too long; %d > max = %d\n",nameS,slen,MAX_TWEAKSEQ);
        return(FALSE);
    }
    /***
    *   Bound region to be modified and check things fit
    */
    firstb = IS_BOG(stPO->firstb) ? 0: stPO->firstb - 1;
    lastb = IS_BOG(stPO->lastb) ? slen: stPO->lastb;
    LIMIT_NUM(firstb,0,slen);
    LIMIT_NUM(lastb,0,slen);
    if( (lastb-firstb) < 2) {
        printf("Cannot shuffle less than 2 bases: %d to %d won't work\n",firstb+1,lastb);
        return(FALSE);
    }
    /***
    *   Copy starting seq so we don't mess with it
    */
    strcpy(newS, seqPC);
    mlen = lastb - firstb;
    /***
printf("first = %d, last = %d, mlen=%d, slen=%d\n",firstb,lastb,mlen,slen);
    *   Get shuffle index array; Only need the number we're shuffling 
    *   Shuffled indices should be unique (i.e. 0 to mlen)
    */
    ArrayRandSequenceI(shuffIA,mlen,mlen,TRUE);
/*
    DumpArray(shuffIA,mlen,NULL,IS_INT,NULL);
*/
    if(verb) {
        fprintf(outPF,"# Shuffling sequence bases %d to %d\n",firstb+1,lastb); 
    }
    for(i=0;i<mlen;i++) {
        if(verb) {
            printf("# Switching [%d %d]\n",i + firstb,shuffIA[i] + firstb);
        }
        newS[i + firstb] = seqPC[shuffIA[i] + firstb];
    }
    newS[slen] = '\0';
    fprintf(outPF,"%s-orig\t%s\n",nameS,seqPC);
    fprintf(outPF,"%s-shuf\t%s\n",nameS,newS);
    return(TRUE);
}
/**************************************************************************
*   Fill passed mask with positions for tweaks
*   tlen = working length
*   firstb = ofset into working length
*   maskS = mask to be filled
*   slen = length of mask
*/
int FillTweakMaskI(SEQTWEAK *stPO, int tlen, int firstb, char *maskS, int slen)
{
    int n,i;
    char *maskPC;

    /***
    *   Fill mask with what to do before doing anything
    *   Then get sub-seq to actually mess with
    */
    for(i=0;i<slen;i++)
    {
        maskS[i]='.';
    }
    maskPC = &maskS[firstb];
    /***
    *   Ins/dels are first, as these may be in blocks of size ids; 
    *       mismatches can fill in holes after
    */
    if( (stPO->ins + stPO->del) > 0) {
        n = FillTweakInDelMaskI(stPO, maskPC, tlen);
        if( n < (stPO->ins + stPO->del) ) {
            return(FALSE);
        }
    }
    /***
    *   Now mismatches
    */
    if(stPO->do_smm) {
        n = FillTweakMissSpeedMaskI(stPO, maskPC,tlen);
    }
    else if(stPO->do_cpm) {
        n = FillTweakMissClosePackMaskI(stPO, maskPC,tlen);
    }
    else {
        n = FillTweakMissShotgunMaskI(stPO, maskPC,tlen);
    }
    if(n<stPO->mis) {
        return(FALSE);
    }
    /***
    *   Clean temp markings in mask
    */
    for(i=0;i<slen;i++)
    {
        if(maskS[i] != toupper(maskS[i])) {
            maskS[i]='.';
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Fill passed mask with positions for InDel tweaks
*   Returns the number of InDels inserted
*/
int FillTweakInDelMaskI(SEQTWEAK *stPO, char *maskS, int slen)
{
    int i,n,ndel,nins,ran,rran,trys,mark;

    ndel = stPO->del;
    nins = stPO->ins;
    rran = slen - stPO->ids + 1;
    trys = n = 0;
    while( ((nins+ndel)>0) && (trys<MAX_TWK_TRYS) )
    {
        trys++;
        ran = RandI(rran); 
        mark = TRUE;
        for(i=0;i<stPO->ids;i++)
        {
            if(maskS[ran+i]!='.') {
                mark=FALSE;
                break;
            }
        }
        if(!mark) {
            continue;
        }
        if(nins>0) {
            for(i=0;i<stPO->ids;i++)
            {
                maskS[ran+i]='I';
            }
            nins--;
            n++;
        }
        else if(ndel>0) {
            for(i=0;i<stPO->ids;i++)
            {
                maskS[ran+i]='D';
            }
            ndel--;
            n++;
        }
    }
    return(n);
}
/**************************************************************************
*   Place masking in a smart (fast) way 
*/
int FillTweakMissShotgunMaskI(SEQTWEAK *stPO, char *maskS,int slen)
{
    int n,trys,ran,rran,nmis;

    nmis = stPO->mis;
    trys = 0;
    rran = slen;
    n = 0;
    while( (nmis>0) && (trys<MAX_TWK_TRYS) )
    {
        trys++;
        ran = RandI(rran);
        if(OkDisSpaceI(ran,maskS,slen,stPO->mds)) {
            MaskDisruption(ran,maskS,slen,stPO->mds);
            nmis--;
            n++;
        }
    }
    return(n);
}
/**************************************************************************
*   Place masking in a close packed way
*/
int FillTweakMissClosePackMaskI(SEQTWEAK *stPO, char *maskS, int slen)
{
    int i,n;

    i = n = 0;
    while( (n < stPO->mis) && (i<slen) )
    {
        MaskDisruption(i,maskS,slen,stPO->mds);
        i += (1 + stPO->mds);
        n++;
    }
    return(n);
}
/**************************************************************************
*   Place masking in a smart (fast) way 
*/
int FillTweakMissSpeedMaskI(SEQTWEAK *stPO, char *maskS,int slen)
{
    int n,m,k,len,min,ran,rran;

    m = stPO->mis;
    n = k = 0;
    while(m>0)
    {
        len = slen - k;
        /***
        *   Minimum space sequence, given spacing between disruptions
        */
        min = m + (m-1) * stPO->mds;
        /***
        *   Range to disrupt = full seq - minimum 
        */
        rran = len - min;
        ran = RandI(rran);
/*
printf("m=%d len=%d min=%d rran=%d k=%d ran=%d dis=%d\n", m,len,min,rran,k,ran,ran|k);
*/
        MaskDisruption(ran+k,maskS,slen,stPO->mds);
        k += (ran + stPO->mds + 1);
        n++;
        m--;
    }
/*
    printf("Smarty |");
    PrintString(maskS,slen,NULL);
    printf("|\n");
    n = FillTweakMissInDelMaskI(nmis, mds, maskS,slen);
*/
    return(n);
}
/**************************************************************************
*   Check that the spacing between things is ok
*/
int OkDisSpaceI(int ran,char *maskS,int slen,int mds)
{
    int i;

    if(isupper(INT(maskS[ran]))) {
        return(FALSE);
    }
    for(i=1;i<=mds;i++)
    {
        if((ran+i)>=slen) {
            break;
        }
        if(isupper(INT(maskS[ran+i]))) {
            return(FALSE);
        }
    }
    for(i=1;i<=mds;i++)
    {
        if((ran-i)<0) {
            break;
        }
        if(isupper(INT(maskS[ran-i]))) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Mask out adjacent parts of mask
*/
void MaskDisruption(int ran,char *maskS,int slen,int mds)
{
    int i;

    maskS[ran]='M';
    for(i=1;i<=mds;i++)
    {
        if((ran+i)>=slen) {
            break;
        }
        maskS[ran+i]='m';
    }
    for(i=1;i<=mds;i++)
    {
        if((ran-i)<0) {
            break;
        }
        maskS[ran-i]='m';
    }
}
/************************************************************************
*   Modify name of new guy
*/
void SetTweakedName(SEQTWEAK *stPO,char *nameS, char *newS, int len)
{
    strcpy(newS,nameS);
    if(!NO_S(stPO->bname)) {
        strcat(newS,"_");
        strcat(newS,stPO->bname);
    }
    else if (stPO->do_nsim) {
        SetTweakedSimName(stPO, nameS, newS);
    }
    else {
        SetTweakedNmodName(stPO, nameS, newS, len);
    }
    return;
}
/************************************************************************
*   Simple (old) naming 
*/
void SetTweakedSimName(SEQTWEAK *stPO,char *nameS, char *newS)
{
    char bufS[DEF_BS];

    if(stPO->ins>0) {
        sprintf(bufS,"_I%d=%d",stPO->ids,stPO->ins);
        strcat(newS,bufS);
    }
    if(stPO->del>0) {
        sprintf(bufS,"_D%d=%d",stPO->ids,stPO->del);
        strcat(newS,bufS);
    }
    if(stPO->mis>0) {
        sprintf(bufS,"_M=%d",stPO->mis);
        strcat(newS,bufS);
    }
    return;
}
/************************************************************************
*   Number-modification naming 
*/
void SetTweakedNmodName(SEQTWEAK *stPO,char *nameS, char *newS, int len)
{
    char bufS[DEF_BS];

    stPO->smm_sub[3] = '\0';
    if(stPO->do_nre) {
        sprintf(bufS,"_mmR:%d:%s", len - stPO->smm_pos, stPO->smm_sub);
    }
    else {
        sprintf(bufS,"_mmF:%d:%s", stPO->smm_pos + 1, stPO->smm_sub);
    }
    strcat(newS,bufS);
    return;
}
