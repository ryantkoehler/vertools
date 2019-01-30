/*
* wfutil.c
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
#include <ctype.h>
#include <math.h>
#include "prim.h"
#include "score.h"
#include "wfutil.h"
#include "wordfreq.h"

#define DB_WF if(DB[99])

/****************************************************************************/
WORDFREQ *CreateWordfreqPO(int size,int ald)
{
    WORDFREQ *wfPO;

    wfPO = (WORDFREQ*)ALLOC(1,sizeof(WORDFREQ));
    if(!wfPO) {
        ERR("CreateWordfreqPO","Failed to allocate WORDFREQ");
        return(NULL);
    }
    wfPO->ID = WORDFREQ_ID;
    wfPO->size = size;
    wfPO->ald = ald;
    wfPO->n = CalcInDimI(size,ald);
    wfPO->pmax = wfPO->n / wfPO->ald;
    wfPO->freqs = (FREC *)ALLOC(wfPO->n,sizeof(FREC));
    if(!wfPO->freqs) {
        printf("Attempted n = %d\n",wfPO->n);
        CHECK_FREE(wfPO);
        ERR("CreateWordfreqPO","Failed to allocate freq records");
        return(NULL);
    }
    InitFrecIDs(wfPO->freqs,wfPO->n);
    return(wfPO);
}
/****************************************************************************/
int DestroyWordfreqI(WORDFREQ *wfPO)
{
    VALIDATE(wfPO,WORDFREQ_ID);
    CHECK_FREE(wfPO->freqs);
    CHECK_FREE(wfPO);
    return(TRUE);
}
/****************************************************************************/
void InitFrecIDs(FREC *frecsPO,int nw)
{
    int i;

    for(i=0;i<nw;i++)
    {   
        frecsPO[i].id = i;  
    }
}
/****************************************************************************/
int CalcInDimI(int size,int aldim)
{
    int i,pmax;

    pmax = 1;
    for(i=0;i<size;i++) 
    {   
        pmax *= aldim; 
    }
    return(pmax);
}
/****************************************************************************
*   fasPO = record with sequence
*/
int TallyWordsI(SEQ *seqPO, WORDFREQ *wfPO, int step)
{
    int i,ind;

    VALIDATE(seqPO,SEQ_ID);
    VALIDATE(wfPO,WORDFREQ_ID);
    /***
    *   The indices are constructed inverse so seqs are in alphabetic order
    */
    for(i=0;i<=(seqPO->len - wfPO->size);i += step) {
        ind = IndexFromSeqI(&seqPO->seq[i],wfPO->size,wfPO->n,wfPO->ald);
        if(ind >= 0) {
            wfPO->freqs[ind].n += 1.0;
        }
    }
    return(TRUE);
}
/****************************************************************************
*   fPO = frequence record array[nw] 
*   degen is a flag to dump degnerate seqs (seq & comp) together
*   lo and hi restrict range of counts reported
*/
int DumpWordsI(WORDFREQ *wfPO, int dgen, DOUB loD, DOUB hiD, FILE *outPF)
{
    int i,cind;
    char wordS[NSIZE],compS[NSIZE];
    DOUB wcD,ccD,bothD,totD,tsumD;

    VALIDATE(wfPO,WORDFREQ_ID);
    HAND_NFILE(outPF);
    totD = TotalWordsD(wfPO,loD,hiD);
    if(totD < 1.0) {
        printf("No word counts to report!\n");
        return(FALSE);
    }
    tsumD = TotalWordsD(wfPO, 0.0, TOO_BIG_D);
    fprintf(outPF,"#ALPHDIM %d\n",wfPO->ald);
    /***
    *   Header line (still commented so can reload
    */
    fprintf(outPF,"#\n");
    if(dgen) {
        fprintf(outPF,"# Word\tComp\tCount-N\tPercent\tWord-N\tComp-N\n");
    }
    else {
        fprintf(outPF,"# Word\tCount-N\tPercent\n");
    }
    /***
    *   Dump each word occuring lo to hi times
    */
    INIT_S(wordS);
    INIT_S(compS);
    for(i=0;i<wfPO->n;i++) {
        wcD = ccD = bothD = 0.0;
        /***
        *   Degenerate word case?
        */
        if(dgen) {
            /***
            *   Get seq & compliment (wordS and compS)
            *   If comp index < current don't report (inverse alphabet order)
            */
            SeqFromIndexI(wfPO->freqs[i].id,wfPO->pmax,wfPO->ald,wordS);
            cind = CompIndexI(wfPO->freqs[i].id,wfPO->pmax,wfPO->ald);
            SeqFromIndexI(cind,wfPO->pmax,wfPO->ald,compS);
            if(cind < wfPO->freqs[i].id) {
                continue;
            }
            wcD = wfPO->freqs[i].n;
            ccD = wfPO->freqs[cind].n;
            bothD = wcD + ccD;
        }
        else {
            SeqFromIndexI(wfPO->freqs[i].id,wfPO->pmax,wfPO->ald,wordS);
            bothD = wfPO->freqs[i].n;
        }
        /***
        *   If frequency is not in range, don't report
        */
        if( (bothD < loD) || (bothD > hiD) ) {
            continue;
        }
        /***
        *   Report the whip
        */
        if(dgen) {
            fprintf(outPF,"%s\t%s\t%5.0f\t%8.5f\t%5.0f\t%5.0f\n",wordS, compS,
                bothD, PERCENT_R(bothD,tsumD), wcD, ccD);
        }
        else {
            fprintf(outPF,"%s\t%8.0f\t%8.5f\n",wordS, bothD,
                 PERCENT_R(bothD,tsumD));
        }
    }
    return(TRUE);
}
/****************************************************************************/
void WordfreqSummary(WORDFREQ *wfPO, DOUB loD, DOUB hiD, FILE *outPF)
{
    int dif;
    DOUB totD,fracD;

    VALIDATE(wfPO,WORDFREQ_ID);
    HAND_NFILE(outPF);
    totD = TotalWordsD(wfPO,loD,hiD);
    dif = TotalDifWordsI(wfPO,loD,hiD);
    fracD = RNUM(dif) / RNUM(wfPO->n);
    fprintf(outPF,"# WordSize     %d\n",wfPO->size);
    if( (loD == 0.0) && (hiD == TOO_BIG_R) ) {
        fprintf(outPF,"# FreqCounts   %1.1f to %1.1f reported\n",loD,hiD);
    }
    else {
        fprintf(outPF,"# FreqCounts   All reported (%1.1f to %1.1f)\n",loD,hiD);
    }
    fprintf(outPF,"# WordsTotal   %1.0f\n",totD);
    fprintf(outPF,"# WordsDif     %d\n",dif);
    fprintf(outPF,"# WordsFrac    %f (%d of %d)\n",fracD,dif,wfPO->n);
}
/****************************************************************************
*   Returns the sum of words with counts in range low to high
*/
DOUB TotalWordsD(WORDFREQ *wfPO,DOUB loD,DOUB hiD)
{
    int i;
    DOUB totD;

    VALIDATE(wfPO,WORDFREQ_ID);
    totD = 0.0;
    for(i=0;i<wfPO->n;i++)
    {
        if( (wfPO->freqs[i].n >= loD) && (wfPO->freqs[i].n <= hiD) ) {
            totD += wfPO->freqs[i].n;
        }
    }
    return(totD);
}
/****************************************************************************
*   Returns the number of differing words
*/
int TotalDifWordsI(WORDFREQ *wfPO,DOUB loD,DOUB hiD)
{
    int i,tot;

    VALIDATE(wfPO,WORDFREQ_ID);
    tot = 0;
    for(i=0;i<wfPO->n;i++)
    {
        if( (wfPO->freqs[i].n >= loD) && (wfPO->freqs[i].n <= hiD) ) {
            tot++;
        }
    }
    return(tot);
}
/****************************************************************************
*   Given an index, fill the corresponding seq into wordS
*/
int SeqFromIndexI(int ind,int pmax,int aldim,char *wordS)
{
    int i,len,p;

    i = len = 0;
    p = pmax;
    while(p>=1)
    {
        switch((ind/p)%aldim)
        {
            case 0: wordS[i++] = 'A';   len++;  break;
            case 1: wordS[i++] = 'C';   len++;  break;
            case 2: wordS[i++] = 'G';   len++;  break;
            case 3: wordS[i++] = 'T';   len++;  break;
            default:    
                wordS[i++] = '?';   break;
        }
        p /= aldim;
    }
    wordS[len] = '\0';
    return(len);
}
/****************************************************************************
*   Take sequence (of length len) and return an index 
*   Weight first chars more for alphabetic sorting
*/
int IndexFromSeqI(char *seqS, int len, int pmax, int aldim)
{
    int j,ind,p;

    ind = 0;
    p = pmax;
    for(j=0;j<len;j++) {
        p /= aldim;
        switch(seqS[j]) {
            case 'a': case 'A':                 break; 
            case 'c': case 'C': ind += p;       break;
            case 'g': case 'G': ind += (2*p);   break;
            case 't': case 'T': ind += (3*p);   break;
            default:    return(BOGUS);
        }
    }
    return(ind);
}
/***************************************************************************
*   Full sequence string from index, given settings in WORDFREQ struct
*/
int FillWordfreqSeqStringI(WORDFREQ *wfPO, int ind, char *seqS)
{
    VALIDATE(wfPO,WORDFREQ_ID);
    INIT_S(seqS);
    if( (ind<0) || (ind>=wfPO->n) ) {
        return(FALSE);
    }
    SeqFromIndexI(ind,wfPO->pmax,wfPO->ald,seqS);
    return(TRUE);
} 
/***************************************************************************
*   Index of word's compliment
*/
int CompIndexI(int ind,int pmax,int aldim)
{
    int p,c,comp;

    DB_WF DB_PrI(">> CompIndexI pmax=%d ind=%d\n+ ",pmax,ind);
    comp = 0;
    c = 1;
    p = pmax;
    while(p>=1)
    {
        DB_WF DB_PrI(" (%d)=",(ind/p)%aldim);
        switch((ind/p)%aldim)
        {
            case 0: comp += (3*c);  break;
            case 1: comp += (2*c);  break;
            case 2: comp += c;      break;
            case 3:                 break;
            default:    
                DB_WF DB_PrI("<< CompIndexI BOGUS\n");
                return(BOGUS);
        }
        DB_WF DB_PrI("%d",comp);
        p /= aldim;
        c *= aldim;
    }
    DB_WF DB_PrI("\n<< CompIndexI %d\n",comp);
    return(comp);
}
/***************************************************************************
*   Collapse counts of degenerate sequences (self + comp) 
*   First-alphabetic seq of each pair gets count of first+second
*   Second-alphabetic seq of each pair gets zero
*/
void CollapseDegenRecs(WORDFREQ *wfPO)
{
    int i,cind;
    
    VALIDATE(wfPO,WORDFREQ_ID);
    for(i=0;i<wfPO->n;i++)
    {
        cind = CompIndexI(wfPO->freqs[i].id,wfPO->pmax,wfPO->ald);
        if(cind > wfPO->freqs[i].id) {
            wfPO->freqs[i].n += wfPO->freqs[cind].n;
            wfPO->freqs[cind].n = 0;
        }
    }
}
/*************************************************************************
*   Sort for output
*/
void SortWordFreqs(WORDFREQ *wfPO,int how)
{
    VALIDATE(wfPO,WORDFREQ_ID);
    switch(how)
    {
        case SORT_HILO:
            qsort(wfPO->freqs,wfPO->n,sizeof(FREC),qSortFrecsI);
            break;
        case SORT_IND:
            qsort(wfPO->freqs,wfPO->n,sizeof(FREC),qSortIndsI);
            break;
        default:
            printf("Bad code = %d\n",how);
            ERR("SortWordFreqs","Bogus sort code");
    }
}
/*************************************************************************/
int qSortFrecsI(const void *e1, const void *e2)
{
    /***
    *   On frequency
    */
    if( ((FREC *)e1)->n  > ((FREC *)e2)->n ) {
        return(-1);
    }
    if( ((FREC *)e1)->n  < ((FREC *)e2)->n ) {
        return(1);
    }
    /***
    *   Index for alphabetic output on freq tie
    */
    if( ((FREC *)e1)->id  < ((FREC *)e2)->id ) {
        return(-1);
    }
    if( ((FREC *)e1)->id  > ((FREC *)e2)->id ) {
        return(1);
    }
    return(0);
}
/*************************************************************************/
int qSortIndsI(const void *e1, const void *e2)
{
    if( ((FREC *)e1)->id  < ((FREC *)e2)->id ) {
        return(-1);
    }
    if( ((FREC *)e1)->id  > ((FREC *)e2)->id ) {
        return(1);
    }
    return(0);
}
/*************************************************************************
*   Attempt to load word freqs from file
*/
int GetWordFreqsI(char *nameS,WORDFREQ **wfPPO)
{
    int ok;
    FILE *fPF;
    WORDFREQ *wfPO;

    wfPO = *wfPPO = NULL;
    if(!(fPF = OpenUFilePF(nameS,"r",NULL))) {
        return(FALSE);
    }
    ok = LoadWordFreqDataI(fPF,&wfPO);
    FILECLOSE(fPF);
    if(ok) {
        strcpy(wfPO->name,nameS);
    }
    *wfPPO = wfPO;
    return(ok);
}
/*************************************************************************
*
*/
int LoadWordFreqDataI(FILE *inPF,WORDFREQ **wfPPO)
{
    char bufS[DEF_BS],seqS[NSIZE];
    int n,s,size,nw,ald;
    WORDFREQ *wfPO;
    DOUB fD;

    *wfPPO = NULL;
    /***
    *   Pass one = get size and count lines 
    */
    n = 0;
    size = s = ald = BOGUS;
    while(fgets(bufS,LINEGRAB,inPF))
    {
        if(COM_LINE(bufS)) {
            if(EQSTRING(bufS,"#ALPHDIM",8)) {
                sscanf(bufS,"%*s %d",&ald);
            }
            continue;
        }
        INIT_S(seqS);
        sscanf(bufS,"%s",seqS);
        if(size == BOGUS) {
            s = size = strlen(seqS);
        }
        else {
            s = strlen(seqS);
        }
        if(s!=size) {
            PROBLINE;
            printf("Inconsistent word size\n");
            printf("   Expecting %d\n",size);
            printf("   Found     %d |%s|\n",s,seqS);
            return(FALSE);
        }
        n++;
    }
    if(size < 1) {
        PROBLINE;
        printf("No sequence data read in\n");
        return(FALSE);
    }
    if(IS_BOG(ald)) {
        PROBLINE;
        printf("No ALPHDIM keyword read from seq data\n");
        return(FALSE);
    }
    /***
    *   Check line count
    */
    nw = CalcInDimI(size,ald);
    if(n > nw) {
        PROBLINE;
        printf("%d data lines; max for size %d = %d\n",n,size,nw);
        return(FALSE);
    }
    /***
    *   Allocate for the whip
    */
    wfPO = CreateWordfreqPO(size,ald);
    if(!wfPO) {
        PROBLINE;
        printf("Failed to allocate: %d words size %d (ald=%d)\n",nw,size,ald);
        return(FALSE);
    }
    /***
    *   Rewind file and load counts
    */
    rewind(inPF);
    while(fgets(bufS,LINEGRAB,inPF)) {
        if(COM_LINE(bufS)) {
            continue;
        }
        if(EQSTRING(bufS,"#ALPHDIM",8)) {
            continue;
        }
        INIT_S(seqS);   fD = BAD_R;
        sscanf(bufS,"%s %lf",seqS,&fD);
        s = IndexFromSeqI(seqS,wfPO->size,wfPO->n,wfPO->ald);
        if( (s<0) || (s>=wfPO->n) ) {
            PROBLINE;
            printf("Bad index calculated for word %s %d\n",seqS,s);
            CHECK_WORDFREQ(wfPO);
            return(FALSE);
        }
        wfPO->freqs[s].n = fD;
    }
    /***
    *   Set pointers
    */
    *wfPPO = wfPO;
    return(TRUE);
}
/****************************************************************************
*   Normalize so that array sums to nw; i.e. each possible word = 1 on average
*/
void NormalizeFrecs(WORDFREQ *wfPO)
{
    int i;
    DOUB denD,totD;

    totD = TotalWordsD(wfPO,0,TOO_BIG);
    denD = totD/DNUM(wfPO->n);
    for(i=0;i<wfPO->n;i++)
    {
        wfPO->freqs[i].n /= denD;
    }
}
/***************************************************************************/
void  LogFrecs(WORDFREQ *wfPO)
{
    int i;
    REAL logR;

    logR = LOG_10(TINY_R);
    for(i=0;i<wfPO->n;i++)
    {
        if(wfPO->freqs[i].n < TINY_R) {
            wfPO->freqs[i].n = logR;
        }
        else {
            wfPO->freqs[i].n = LOG_10(wfPO->freqs[i].n);
        }
    }
}
/*****************************************************************************
*   Check that 2 wordfreq configurations are compatible 
*/
int CompatWordfreqsI(WORDFREQ *fPO, WORDFREQ *sPO,int verb)
{
    if( (!fPO) || (!sPO) ) {
        return(FALSE);
    }
    if(fPO->size != sPO->size) {
        if(verb) {
            printf("Incompatable word sizes: %d %d\n",fPO->size,sPO->size);
        }
        return(FALSE);
    }
    if(fPO->ald != sPO->ald) {
        if(verb) {
            printf("Incompatable alphabet sizes: %d %d\n",fPO->ald,sPO->ald);
        }
        return(FALSE);
    }
    return(TRUE);
}
/*****************************************************************************
*   Combine two wordfreq sets if they are compatible
*/
int MergeWordfreqsI(WORDFREQ *fPO, WORDFREQ *sPO,int how, WORDFREQ *aPO)
{
    int i;

    if( (!CompatWordfreqsI(fPO,sPO,FALSE)) || (!CompatWordfreqsI(fPO,aPO,FALSE)) ) {
        return(FALSE);
    }
    for(i=0;i<aPO->n;i++)
    {
        switch(how)
        {
            case ADD_WF:
                aPO->freqs[i].n = fPO->freqs[i].n + sPO->freqs[i].n;
                break;
            case SUB_WF:
                aPO->freqs[i].n = fPO->freqs[i].n - sPO->freqs[i].n;
                break;
            case DSR_WF:
                if(fPO->freqs[i].n + sPO->freqs[i].n > 0.0)
                {
                    aPO->freqs[i].n = (fPO->freqs[i].n - sPO->freqs[i].n) /
                        (fPO->freqs[i].n + sPO->freqs[i].n);
                }
                else
                {
                    aPO->freqs[i].n = 0.0;
                }
                break;
            case FDSR_WF:
                if(fPO->freqs[i].n + sPO->freqs[i].n > 0.0)
                {
                    aPO->freqs[i].n = (fPO->freqs[i].n - sPO->freqs[i].n) /
                        sqrt(fPO->freqs[i].n + sPO->freqs[i].n);
                }
                else
                {
                    aPO->freqs[i].n = 0.0;
                }
                break;
            default:
                printf("Bogus merge code: %d\n",how);
                ERR("MergeWordfreqsI","Don't know how to merge?");
        }
    }
    return(TRUE);
}
