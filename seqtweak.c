/*
* seqtweak.c
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
    printf("   -iuc -ilc  Ignore (don't modify) Upper / Lower case\n");
    printf("   -mis #     Randomly put in # mismatches\n");
    printf("   -mip #     Randomly put in # percent mismatches (of allowed positions)\n");
    printf("   -del #     Randomly put in # deletions\n");
    printf("   -ins #     Randomly put in # insertions\n");
    printf("   -ids #     Insert/deletion word size (i.e. how many bases)\n");
    printf("   -mds #     Set Minimum Disruption Separation to # bases\n");
    printf("   -amb       All mismatch bases for each mismatch (i.e. 3/seq)\n");
    printf("   -cpm       Close-pack mismatch placement\n");
    printf("   -sh        Shuffle sequence\n");
    printf("   -seed #    Set random seed to #\n");
    printf("   -num #     Number of different tweaks / input; default 1\n");
    printf("   -muc -mlc  Mark tweaks Upper / Lower case\n");
    printf("   -dmask     Mark masking (non-tweakable) positions\n");
    printf("   -bname XXX Append output names with base XXX\n");
    printf("   -nre       Name relative to end (i.e. 3' == 1)\n");
    printf("   -quiet     No verbose description of tweaking; only seqs\n");
}
/**************************************************************************/
int SeqTweakI(int argc,char **argv)
{
    int ok,iraw,iseq,ifas,verb,n;
    SEQTWEAK *stPO;

    stPO = CreateSeqtweakPO();
    verb = TRUE;
    iraw = iseq = ifas = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -iraw B -ifas B -mis I -del I -ins I -ids I -bran I2\
        -seed I -mds I -quiet B -amb B -cpm B -bname S\
        -sh B -rre B -nre B -iseq B -num I\
        -muc B -mlc B -iuc B -ilc B -dmask B -mip I",
        stPO->inname, stPO->outname, &iraw, &ifas, 
        &stPO->mis, &stPO->del, &stPO->ins, &stPO->ids, 
        &stPO->firstb,&stPO->lastb, &stPO->seed, 
        &stPO->mds, &verb, &stPO->do_amb, &stPO->do_cpm, &stPO->bname,
        &stPO->do_sh, &stPO->do_rre,
        &stPO->do_nre, &iseq, &stPO->num,
        &stPO->do_muc, &stPO->do_mlc, &stPO->do_iuc, &stPO->do_ilc, 
        &stPO->do_dmask, &stPO->mip,
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
    n = 0;
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
        stPO->seq_n += 1;
        if( (!UpdateSeqVarsI(stPO)) || (!TweakableSeqI(stPO)) ) {
            continue;
        }
        /***
        *   N number of tweaks
        */
        stPO->tweak_n = 0;
        while(stPO->tweak_n < stPO->num)
        {
            if(!InitTweakSeqMaskI(stPO)) {
                break;
            }
            stPO->tweak_n += 1;
            if(stPO->do_sh) {
                n += TweakShuffleSeqI(stPO,verb,stPO->out);
            }
            else {
                n += TweakSeqBasesI(stPO,verb,stPO->out);
            }
        }
    }
    if(verb) {
        printf("# Input sequences:  %d\n",stPO->seq_n);
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
    stPO->mip = BOGUS;
    stPO->ins = 0;
    stPO->del = 0;
    stPO->ids = 1;
    stPO->mds = 0;
    stPO->do_amb = FALSE;
    stPO->do_cpm = FALSE;
    stPO->do_sh = FALSE;
    stPO->do_muc = stPO->do_mlc = FALSE;
    stPO->do_iuc = stPO->do_ilc = FALSE;
    stPO->do_dmask = FALSE;
    stPO->seed = BOGUS;
    INIT_S(stPO->bname);
    stPO->do_nre = FALSE;
    stPO->num = 1;
    stPO->seq_n = stPO->tweak_n = 0;
    INIT_S(stPO->mask);
    INIT_S(stPO->tseq1);
    INIT_S(stPO->tseq2);
    INIT_S(stPO->tseq3);
    INIT_S(stPO->mbase);
}
/**************************************************************************
*   Check options
*/
int OkSeqTweakOptsI(SEQTWEAK *stPO)
{
    if(stPO->ids > MAX_IDSIZE) {
        printf("# In-del size %d too big (max %d)\n",stPO->ids,MAX_IDSIZE);
        return(FALSE);
    }
    if((!stPO->do_sh) && (stPO->mis<1) && (IS_BOG(stPO->mip)) && (stPO->del<1) && (stPO->ins<1)){
        printf("# No mismatch, No insert, No deletion, Not shuffle = nothing to do!\n");
        return(FALSE);
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
    *   Total modify base count
    */
    stPO->tmod = stPO->mis + (stPO->del * stPO->ids) + (stPO->ins * stPO->ids);
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
    fprintf(outPF,"# Mismatchs: %d (count)\n",stPO->mis);
    fprintf(outPF,"# Mismatchs: %d (percent)\n",stPO->mip);
    fprintf(outPF,"# Insertions: %d (%d bases each)\n",stPO->ins,stPO->ids);
    fprintf(outPF,"# Deletions:  %d (%d bases each)\n",stPO->del,stPO->ids);
    fprintf(outPF,"# Minimum disruption separation: %d\n",stPO->mds);
    if( ! IS_BOG(stPO->firstb) ) {
        fprintf(outPF,"# Base range: %d %d\n",stPO->firstb,stPO->lastb);
    }
    if(stPO->do_iuc) {
        fprintf(outPF,"# Ignore (don't change) Uppercase bases\n");
    }
    if(stPO->do_ilc) {
        fprintf(outPF,"# Ignore (don't change) lowercase bases\n");
    }
    if(stPO->do_amb) {
        fprintf(outPF,"# All (3) mismatch bases / mismatch\n");
    }
    if(stPO->do_cpm) {
        fprintf(outPF,"# Close-pack disruption placement\n");
    }
    else {
        fprintf(outPF,"# Shotgun disruption placement\n");
    }
    fprintf(outPF,"#\n");
    fprintf(outPF,"\n");
}
/***************************************************************************
*   Set run-time vars for seq into seqtweak
*/
int UpdateSeqVarsI(SEQTWEAK *stPO)
{
    stPO->seqlen = GetSeqLenI(stPO->seq);
    if(!FillSeqNameStringI(stPO->seq, stPO->seqname, NSIZE)) {
        return(FALSE);
    }
    if(!GetSeqSeqI(stPO->seq, &stPO->seqseq)) {
        return(FALSE);
    }
    return(TRUE);
}
/***************************************************************************
*   Check if ok to mess with
*   If yes, 
*/
int TweakableSeqI(SEQTWEAK *stPO)
{
    if( stPO->seqlen > MAX_TWEAKSEQ ) {
        printf("# Sequence %d too long: %d\n",stPO->seq_n, stPO->seqlen);
        printf("#   Max = %d\n",MAX_TWEAKSEQ);
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************
*   Shuffle sequence 
*/
int TweakShuffleSeqI(SEQTWEAK *stPO, int verb, FILE *outPF)
{
    int i, k, shuffIA[MAX_TWEAKSEQ];
    char nnameS[NSIZE], newS[MAX_TWEAKSEQ], nbaseS[MAX_TWEAKSEQ];

    HAND_NFILE(outPF);
    /***
    *   Get shuffle index array; Only need the number we're shuffling 
    *   Shuffled indices should be unique
    *   Collect tweakable bases (only; ignore non-tweak)
    */
    ArrayRandSequenceI(shuffIA, stPO->ntpos, stPO->ntpos, TRUE);
    k = 0;
    for(i=0; i<stPO->seqlen; i++)
    {
        if(stPO->mask[i] == MASK_YES) {
            nbaseS[k++] = stPO->seqseq[i];
        }
    }
/* 
    DumpArray(shuffIA, IS_INT, 0,stPO->ntpos, NULL, NULL);
*/
    if(verb) {
        fprintf(outPF,"# Shuffling sequence bases\n");
    }
    /***
    * Copy starting seq to new
    * For open positions, copy from rand[ind] to current pos; 
    * For closed positions, copy direct source > dest
    */
    k = 0;
    for(i=0; i<stPO->seqlen; i++)
    {
        if(stPO->mask[i] == MASK_YES) {
            newS[i] = nbaseS[shuffIA[k++]];
        }
        else {
            newS[i] = stPO->seqseq[i];
        }
    }
    newS[stPO->seqlen] = '\0';
    SetTweakedName(stPO, nnameS);
    if(verb) {
        ReportOldMaskOuts(stPO, nnameS, outPF);
    }
    fprintf(outPF,"%s\t%s\n", nnameS, newS);
    return(TRUE);
}
/**************************************************************************
*   Tweak bases for current seq
*/
int TweakSeqBasesI(SEQTWEAK *stPO, int verb, FILE *outPF)
{
    int i,j,k,trys,ok;
    char nnameS[NSIZE], twkS[DEF_BS];
    char bC;

    HAND_NFILE(outPF);
    /***
    *   Try to find a tweak pattern that fits in allowed range
    */
    ok = trys = 0;
    while(trys<MAX_TWK_STARTS) 
    {
        if(FillTweakMaskI(stPO)) {
            ok++;
            break;
        }
        trys++;
    }
    if(!ok) {
        printf("# Failed to place in-del / mismatch positions\n");
        printf("#   No tweaks for %s\n",stPO->seqname);
        return(FALSE);
    }
    /***
    *   Do the subs indicated in mask
    *   i scans mask, j = source (start) seq, k = dest (tweaked) seq
    */
    j = k = 0;
    for(i=0; i<stPO->seqlen; i++)
    {
        switch(stPO->mask[i])
        {
            case MASK_INSERT:   
                bC = DNAIndexBaseC(RandI(4));
                stPO->tseq1[k] = stPO->tseq2[k] = stPO->tseq3[k] = bC;
                stPO->mbase[k] = TRUE;
                k++;
                break;
            case MASK_DELETE:
                j++;
                break;
            case MASK_MISMATCH:
                if(stPO->do_amb)
                {
                    switch(stPO->seqseq[j])
                    {
                        case 'a': case 'A':
                            stPO->tseq1[k] = 'C';
                            stPO->tseq2[k] = 'G';
                            stPO->tseq3[k] = 'T';
                            break;
                        case 'c': case 'C':
                            stPO->tseq1[k] = 'A';
                            stPO->tseq2[k] = 'G';
                            stPO->tseq3[k] = 'T';
                            break;
                        case 'g': case 'G':
                            stPO->tseq1[k] = 'A';
                            stPO->tseq2[k] = 'C';
                            stPO->tseq3[k] = 'T';
                            break;
                        case 't': case 'T':
                            stPO->tseq1[k] = 'A';
                            stPO->tseq2[k] = 'C';
                            stPO->tseq3[k] = 'G';
                            break;
                    }
                }
                else {
                    bC = DNAIndexBaseC(RandI(4));
                    while(TOUPPER(bC) == TOUPPER(stPO->seqseq[j]))
                    {
                        bC = DNAIndexBaseC(RandI(4));
                    }
                    stPO->tseq1[k] = bC;
                }
                stPO->mbase[k] = TRUE;
                k++;
                j++;
                break;
            default:
                stPO->tseq1[k] = stPO->tseq2[k] = stPO->tseq3[k] = stPO->seqseq[j];
                j++;
                k++;
                break;
        }
    }
    while(k < stPO->seqlen)
    {
        if(j < stPO->seqlen) {
            bC = stPO->seqseq[j++];
            stPO->tseq1[k] = stPO->tseq2[k] = stPO->tseq3[k] = bC;
        }
        else {
            bC = DNAIndexBaseC(RandI(4));
            stPO->tseq1[k] = stPO->tseq2[k] = stPO->tseq3[k] = bC;
            stPO->mbase[k] = TRUE;
        }
        k++;
    }
    stPO->tseq1[k] = stPO->tseq2[k] = stPO->tseq3[k] = '\0';
    /***
    *   Tell the story
    */
    ChangeDisMaskSeqCase(stPO);
    FillDisruptString(stPO, twkS);
    SetTweakedName(stPO, nnameS);
    if(verb) {
        fprintf(outPF,"# Disrupts %d %s\n",(stPO->ins+stPO->del+stPO->mis),twkS);
        ReportOldMaskOuts(stPO, nnameS, outPF);
    }
    /***
    *   Special output names for All-Missmatch-Base case
    */
    if(stPO->do_amb) {
        fprintf(outPF,"%s_x\t%s\n", stPO->seqname, stPO->tseq1);
        fprintf(outPF,"%s_y\t%s\n", stPO->seqname, stPO->tseq2);
        fprintf(outPF,"%s_z\t%s\n", stPO->seqname, stPO->tseq3);
    }
    else {
        fprintf(outPF,"%s\t%s\n", nnameS, stPO->tseq1);
    }
    if(verb) {
        fprintf(outPF,"\n");
    }
    return(TRUE);
}
/**************************************************************************/
void ReportOldMaskOuts(SEQTWEAK *stPO, char *nnameS, FILE *outPF)
{
    int padsize;
    char bufS[NSIZE+1];

    HAND_NFILE(outPF);
    padsize = stPO->do_amb ? strlen(stPO->seqname)+2 : strlen(nnameS);
    sprintf(bufS,"#%s", stPO->seqname);
    PadString(bufS, ' ', padsize, bufS);
    fprintf(outPF,"%s\t%s\n",bufS, stPO->seqseq);
    sprintf(bufS,"#");
    PadString(bufS, ' ', padsize, bufS);
    stPO->mask[stPO->seqlen] = '\0';
    fprintf(outPF,"%s\t%s\n", bufS, stPO->mask);
}
/**************************************************************************
*   Change distruption masking and sequence case depending on settings
*/
void ChangeDisMaskSeqCase(SEQTWEAK *stPO) 
{
    int i;

    stPO->mask[stPO->seqlen] = '\0';
    if(stPO->do_muc) {
        Lowerize(stPO->tseq1);
        Lowerize(stPO->tseq2);
        Lowerize(stPO->tseq3);
    }
    else if(stPO->do_mlc) {
        Upperize(stPO->tseq1);
        Upperize(stPO->tseq2);
        Upperize(stPO->tseq3);
    }
    for(i=0; i<stPO->seqlen; i++)
    {
        /***
        * If open (YES) ignore; If masked (NO), set yes then ignore
        */
        if(stPO->mask[i]==MASK_YES) {
            continue;
        }
        if(stPO->mask[i]==MASK_NO) {
            if(! stPO->do_dmask) {
                stPO->mask[i] = MASK_YES; 
            }
            continue;
        }
        /***
        *   Marking changes with case? The mbase masking has modified bases
        */
        if(stPO->mbase[i]) {
            if(stPO->do_muc) {
                stPO->tseq1[i] = TOUPPER(stPO->tseq1[i]);
                stPO->tseq2[i] = TOUPPER(stPO->tseq2[i]);
                stPO->tseq3[i] = TOUPPER(stPO->tseq3[i]);
            }
            else if (stPO->do_mlc) {
                stPO->tseq1[i] = TOLOWER(stPO->tseq1[i]);
                stPO->tseq2[i] = TOLOWER(stPO->tseq2[i]);
                stPO->tseq3[i] = TOLOWER(stPO->tseq3[i]);
            }
        }
    }
    return;
}
/**************************************************************************
*   Fill disruption annotation string based on mask
*/
void FillDisruptString(SEQTWEAK *stPO, char *twkS) 
{
    int i,k;
    char bufS[DEF_BS];

    INIT_S(twkS);
    k = 0;
    for(i=0; i<stPO->seqlen; i++)
    {
        if( (stPO->mask[i]==MASK_YES) || (stPO->mask[i]==MASK_NO) ) {
            continue;
        }
        /***
        * Annotation of indel or single base collected, then appended
        */
        if( (stPO->mask[i]==MASK_INSERT) || (stPO->mask[i]==MASK_DELETE) ) {
            sprintf(bufS,"%c%d=%d", stPO->mask[i], stPO->ids, i+1);
            i+= (stPO->ids-1);
        }
        else if(stPO->mask[i]==MASK_MISMATCH) {
            sprintf(bufS,"%c=%d", stPO->mask[i], i+1);
        }
        /* Add space if we've already got a change */
        if(k>0) {
            strcat(twkS," ");
        }
        strcat(twkS,bufS);
        k++;
    }
    return;
}
/**************************************************************************
*   Initialize tweak mask for current seq
*   Return FALSE if tweaks can't fit in sequence
*/
int InitTweakSeqMaskI(SEQTWEAK *stPO)
{
    int firstb, lastb, i;

    /***
    * Initialize mask to non-modifable, seqs to NULL
    */
    InitArrayI(stPO->mask, IS_CHAR, 0, MAX_TWEAKSEQ, MASK_NO);
    InitArrayI(stPO->tseq1, IS_CHAR, 0, MAX_TWEAKSEQ, '\0');
    InitArrayI(stPO->tseq2, IS_CHAR, 0, MAX_TWEAKSEQ, '\0');
    InitArrayI(stPO->tseq3, IS_CHAR, 0, MAX_TWEAKSEQ, '\0');
    InitArrayI(stPO->mbase, IS_CHAR, 0, MAX_TWEAKSEQ, FALSE);
    /***
    *   Bound region to be modified? (if fisrtb has been set)
    */
    if( ! IS_BOG(stPO->firstb) ) {
        firstb = (stPO->do_rre) ? (stPO->seqlen - stPO->lastb + 1)  : stPO->firstb - 1;
        lastb = (stPO->do_rre)  ? (stPO->seqlen - stPO->firstb + 1) : stPO->lastb;
    }
    else {
        firstb = 0;
        lastb = stPO->seqlen;
    }
    LIMIT_NUM(firstb,0,stPO->seqlen);
    LIMIT_NUM(lastb,0,stPO->seqlen);
    /*** 
    *   Set modifiable part
    */
    stPO->sti = firstb;
    stPO->eni = lastb;
    stPO->ntpos = 0;
    for(i=stPO->sti; i<stPO->eni; i++)
    {
        /* Ignore uppercase? */
        if(stPO->do_iuc && ISUPPER(stPO->seqseq[i])) {
            continue;
        }
        /* Ignore lowercase? */
        if(stPO->do_ilc && ISLOWER(stPO->seqseq[i])) {
            continue;
        }
        stPO->mask[i] = MASK_YES;
        stPO->ntpos++;
    }
    /***
    *   If already not enough space, bail (doesn't count between spacing)
    */
    if(stPO->ntpos < stPO->tmod) {
        printf("# Too many modifications %d\n", stPO->tmod);
        printf("#   Only %d total positions available\n", stPO->ntpos);
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************
*   Fill mask with positions for tweaks
*/
int FillTweakMaskI(SEQTWEAK *stPO)
{
    int n;

    /***
    *   Ins/dels are first, as these may be bigger (blocks of size ids); 
    *   If don't get requested number of indels, bail false
    */
    if( (stPO->ins + stPO->del) > 0) {
        n = FillTweakInDelMaskI(stPO);
        if( n < (stPO->ins + stPO->del) ) {
            return(FALSE);
        }
    }
    /***
    *   If percentage mismatch, set actual number
    */
    if(!IS_BOG(stPO->mip)) {
        n = INT( DNUM(stPO->ntpos * stPO->mip) / 100.0);
        LIMIT_NUM(n,0,stPO->ntpos);
        stPO->mis = n;
    }
    /***
    *   Now mismatches; Close pack or random; Not enough, bail false
    */
    if(stPO->do_cpm) {
        n = FillTweakMissClosePackMaskI(stPO);
    }
    else {
        n = FillTweakMissShotgunMaskI(stPO);
    }
    if(n<stPO->mis) {
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************
*   Fill passed mask with positions for InDel tweaks
*   Updates (decrements) the number of tweakable positions remaining
*   Returns the number of InDels inserted
*/
int FillTweakInDelMaskI(SEQTWEAK *stPO)
{
    int i,n,ndel,nins,ran,rran,trys,dis;

    ndel = stPO->del;
    nins = stPO->ins;
    rran = stPO->seqlen - stPO->ids + 1;
    /***
    *   Random placements while number tweakable positions > indel size,
    *   number indels > 0, and not too many tries
    */
    trys = n = 0;
    while( (stPO->ntpos >= stPO->ids) && ((nins+ndel) > 0) && (trys < MAX_TWK_TRYS) )
    {
        trys++;
        /***
        *   Random start has to be followed by open (YES) space for indel
        */
        ran = RandI(rran); 
        if(stPO->mask[ran] != MASK_YES) {
            continue;
        }
        dis = ClosestDisruption(stPO, ran, stPO->ids);
        if(dis < stPO->mds) {
            continue;
        }
        /***
        *   Mark and up counts
        */
        if(nins>0) {
            for(i=0;i<stPO->ids;i++)
            {
                stPO->mask[ran+i] = MASK_INSERT;
                stPO->ntpos--;
            }
            nins--;
            n++;
        }
        else if(ndel>0) {
            for(i=0;i<stPO->ids;i++)
            {
                stPO->mask[ran+i] = MASK_DELETE;
                stPO->ntpos--;
            }
            ndel--;
            n++;
        }
        /* Mask any adjacents */
        MaskDisruption(stPO, ran, stPO->ids);
    }
    return(n);
}
/**************************************************************************
*   Place mismatches at random
*/
int FillTweakMissShotgunMaskI(SEQTWEAK *stPO)
{
    int n,trys,ran,rran,nmis,dis;

    nmis = stPO->mis;
    rran = stPO->seqlen;
    /***
    *   Random placements while number tweakable positions > 0,
    *   number mismatchs > 0, and not too many tries
    */
    trys = n = 0;
    while( (stPO->ntpos > 0) && (n < stPO->mis) && (trys < MAX_TWK_TRYS) )
    {
        trys++;
        ran = RandI(rran);
        if(stPO->mask[ran] == MASK_YES) {
            dis = ClosestDisruption(stPO, ran, 1);
            if(dis >= stPO->mds) {
                stPO->mask[ran] = MASK_MISMATCH;
                MaskDisruption(stPO, ran, 1);
                nmis--;
                n++;
            }
        }
    }
    return(n);
}
/**************************************************************************
*   Place mismatches in a close packed way
*/
int FillTweakMissClosePackMaskI(SEQTWEAK *stPO)
{
    int i,n,dis;

    i = n = 0;
    while( (stPO->ntpos > 0) && (n < stPO->mis) && (i < stPO->seqlen) )
    {
        if(stPO->mask[i] == MASK_YES) {
            dis = ClosestDisruption(stPO, i, 1);
            if(dis >= stPO->mds) {
                stPO->mask[i] = MASK_MISMATCH;
                stPO->ntpos--;
                MaskDisruption(stPO, i, 1);
                n++;
            }
        }
        i++;
    }
    return(n);
}
/**************************************************************************
*   Mask out adjacent parts of mask
*   Also decrements number of tweakable positions for any added mask
*/
void MaskDisruption(SEQTWEAK *stPO, int pos, int size)
{
    int i,up,down;

    for(i=0;i<stPO->mds;i++)
    {
        up = pos - i - 1;
        down = pos + i + size;
        if( up >= 0) {
            if(stPO->mask[up] == MASK_YES) {
                stPO->mask[up] = MASK_NO;
                stPO->ntpos--;
            }
        }
        if( down < stPO->seqlen) {
            if(stPO->mask[down] == MASK_YES) {
                stPO->mask[down] = MASK_NO;
                stPO->ntpos--;
            }
        }
    }
}
/**************************************************************************
*   Find the closest 'set' position in mask relative to pos
*/
int ClosestDisruption(SEQTWEAK *stPO, int pos, int size)
{
    int up,down,dis;

    dis = 0;
    up = pos - 1;
    down = pos + size;
    while( (up>=0) || (down<stPO->seqlen) )
    {
        /***
        *   YES and NO are not set = disruption (MIS, IN, DEL)
        *   (only if index is legit)
        */
        if(up>=0) {
            if( (stPO->mask[up] != MASK_YES) && (stPO->mask[up] != MASK_NO) ) {
                break;
            }
        }
        if(down<stPO->seqlen) {
            if( (stPO->mask[down] != MASK_YES) && (stPO->mask[down] != MASK_NO) ) {
                break;
            }
        }
        /* Up distance and expand window */
        dis++;
        up--;
        down++;
    }
    return(dis);
}
/************************************************************************
*   Modify name of new guy
*/
void SetTweakedName(SEQTWEAK *stPO, char *newS)
{
    char ntsufS[DEF_BS];

    strcpy(newS, stPO->seqname);
    if(!NO_S(stPO->bname)) {
        strcat(newS,"_");
        strcat(newS,stPO->bname);
    }
    else if(stPO->do_sh) {
        strcat(newS,"_shuf");
    }
    else {
        SetTweakedSimName(stPO, newS);
    }
    /***
    * Any num-tweak suffix?
    */
    if( stPO->num > 1 ){
        sprintf(ntsufS, "_twk%02d", stPO->tweak_n);
        strcat(newS, ntsufS);
    }
    return;
}
/************************************************************************
*   Simple (old) naming 
*/
void SetTweakedSimName(SEQTWEAK *stPO, char *newS)
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
