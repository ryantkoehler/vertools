/*
* dna_seqs.c
*
* Copyright 2016 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include "prim.h"
#include "dna.h"

#define DB_SSEQ if(DB[103])


/****************************************************************************
*   Sets a pointer to sequence in seq set ssPO[ind]
*/
int GetSeqsetSeqI(SEQSET *ssPO,int ind, SEQ **seqPPO)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(FALSE);
    }
    seqPO = ssPO->seqs[ind];
    *seqPPO = seqPO;
    return(TRUE);
}
/****************************************************************************
*   Sets a pointer to sequence in seq set ssPO[ind]
*/
int GetSeqsetSeqStringI(SEQSET *ssPO,int ind, char **seqPPC)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(FALSE);
    }
    seqPO = ssPO->seqs[ind];
    VALIDATE(seqPO,SEQ_ID);
    *seqPPC = seqPO->seq;
    return(TRUE);
}
/****************************************************************************
*   Copies sequence from seq ind in set ssPO into seqS up to max chars
*   returns BOGUS if bad index or sequence can't fit in passed string
*/
int FillSeqsetSeqStringI(SEQSET *ssPO,int ind,char *seqS,int max)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(BOGUS);
    }

    seqPO = ssPO->seqs[ind];
/* XX sham? */
    LIMIT_NUM(max,0,seqPO->len);
    max = MIN_NUM(max,seqPO->len);
    strncpy(seqS,seqPO->seq,max);
    seqS[max] = '\0';
    return(max);
}
/****************************************************************************
*   Set the name of sequence in seq set ssPO[ind] to name nS[n]
*
*   returns FALSE if index is bad
*/
int SetSeqsetSeqNameI(SEQSET *ssPO,int ind, char *nS, int nlen)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(FALSE);
    }
    LIMIT_NUM(nlen,0,NSIZE-1);
    seqPO = ssPO->seqs[ind];
    VALIDATE(seqPO,SEQ_ID);
    strncpy(seqPO->name, nS, nlen);
    return(TRUE);
}
/****************************************************************************
*   Sets a pointer to name of sequence in seq set ssPO[ind]
*/
int GetSeqsetSeqNameI(SEQSET *ssPO,int ind, char **namePPC)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(FALSE);
    }
    seqPO = ssPO->seqs[ind];
    VALIDATE(seqPO,SEQ_ID);
    *namePPC = seqPO->name;
    return(TRUE);
}
/****************************************************************************
*   Copies the name of seq in set ssPO[ind] into nS up to max chars
*
*   returns FALSE if index is bad
*/
int FillSeqsetSeqNameI(SEQSET *ssPO,int ind,char *nS,int max)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(FALSE);
    }
    seqPO = ssPO->seqs[ind];
    VALIDATE(seqPO,SEQ_ID);
    LIMIT_NUM(max,0,NSIZE-1);
    strncpy(nS,seqPO->name,max);
    nS[max] = '\0';
    return(TRUE);
}
/****************************************************************************
*   Returns the length of sequence in seq set ssPO[ind]
*
*   returns BOGUS is index is bad
*/
int GetSeqsetSeqLenI(SEQSET *ssPO,int ind)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(BOGUS);
    }
    seqPO = ssPO->seqs[ind];
    VALIDATE(seqPO,SEQ_ID);
    return(seqPO->len);
}
/****************************************************************************
*   Returns the number of SNPs in seq set ssPO[ind]
*   returns BOGUS if index is bad
*/
int GetSeqsetSeqSnpCountI(SEQSET *ssPO,int ind)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(BOGUS);
    }
    seqPO = ssPO->seqs[ind];
    VALIDATE(seqPO,SEQ_ID);
    return(seqPO->nsnp);
}
/****************************************************************************
*   Returns the number of sequences in seq set ssPO
*/
int GetSeqsetNumI(SEQSET *ssPO)
{
    VALIDATE(ssPO,SEQSET_ID);
    return(ssPO->n);
}
/****************************************************************************
*   Returns the size of the smallest sequence in seq set ssPO
*/
int SeqsetMinLenI(SEQSET *ssPO)
{
    int min;

    VALIDATE(ssPO,SEQSET_ID);
    SeqsetMinMaxLenI(ssPO,&min,NULL);
    return(min);
}
/****************************************************************************
*   Returns the size of the largest sequence in seq set ssPO
*/
int SeqsetMaxLenI(SEQSET *ssPO)
{
    int max;

    VALIDATE(ssPO,SEQSET_ID);
    SeqsetMinMaxLenI(ssPO,NULL,&max);
    return(max);
}
/****************************************************************************
*   Returns the size of the smallest sequence in seq set ssPO
*/
int SeqsetMinMaxLenI(SEQSET *ssPO,int *minPI,int *maxPI)
{
    int i,n,min,max;
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if(minPI)
    {   *minPI = BOGUS; }
    if(maxPI)
    {   *maxPI = BOGUS; }
    min = TOO_BIG;
    max = -TOO_BIG;
    n = 0;
    for(i=0;i<ssPO->n;i++)
    {
        if(ssPO->mask[i])
        {
            seqPO = ssPO->seqs[i];
            min = MIN_NUM(seqPO->len,min);
            max = MAX_NUM(seqPO->len,max);
            n++;
        }
    }
    if(minPI)
    {   *minPI = min; }
    if(maxPI)
    {   *maxPI = max; }
    return(n);
}
/*****************************************************************************
*   Set name of seqset to nameS
*/
void SetSeqsetName(SEQSET *ssPO,char *nameS)
{
    VALIDATE(ssPO,SEQSET_ID);
    strncpy(ssPO->name,nameS,NSIZE-1);
    ssPO->name[NSIZE-1] = '\0';
}
/*****************************************************************************
*   Get name of seqset to nameS
*/
int FillSeqsetNameStringI(SEQSET *ssPO,char *nameS,int max)
{
    VALIDATE(ssPO,SEQSET_ID);
    LIMIT_NUM(max,0,NSIZE-1);
    INIT_S(nameS);
    strncpy(nameS,ssPO->name,max);
    nameS[max] = '\0';
    return(max);
}
/*****************************************************************************
*   Set source of seqset to nameS
*/
void SetSeqsetSource(SEQSET *ssPO,char *nameS)
{
    VALIDATE(ssPO,SEQSET_ID);
    strncpy(ssPO->source,nameS,NSIZE-1);
    ssPO->source[NSIZE-1] = '\0';
}
/*****************************************************************************
*   Get source (filename) of seqset to nameS
*/
int FillSeqsetSourceStringI(SEQSET *ssPO,char *nameS,int max)
{
    VALIDATE(ssPO,SEQSET_ID);
    LIMIT_NUM(max,0,NSIZE-1);
    INIT_S(nameS);
    strncpy(nameS,ssPO->source,max);
    nameS[max] = '\0';
    return(max);
}
/*****************************************************************************
*   Look for sequence by name in set
*   If kc is true, case is kept in comparision
*   If tC is non-null, names will be truncated at this character 
*   If seqPPO is real, pointer to the (found) SEQ is set
*
*   Returns the index to found sequence, or BOGUS if not found
*/
int FindNamedSeqInSeqsetI(SEQSET *ssPO, char *nameS, int kc, char *tPC,
    SEQ **seqPPO)
{
    char snameS[NSIZE];
    int i,m;

    VALIDATE(ssPO,SEQSET_ID);
    m = BOGUS;
    for(i=0;i<ssPO->n;i++)
    {
        FillSeqsetSeqNameI(ssPO,i,snameS,NSIZE);
        if(tPC)
        {
            ReplaceChars(tPC[0],nameS,'\0',nameS);
        }
        if(kc)
        {
            if(!strcmp(nameS,snameS))
            {
                m = i;
                break;
            }
        }
        else
        {
            if(!strcasecmp(nameS,snameS))
            {
                m = i;
                break;
            }
        }
    }
    /***
    *   Set pointer if real 
    */
    if( (seqPPO) && (m>=0) )
    {
        if(!GetSeqsetSeqI(ssPO,m,seqPPO))
        {
            m = BOGUS;
        }
    }
    return(m);
}
/*****************************************************************************
*   Set sequence of seq to seqS 
*/
int SetSeqSequenceI(SEQ *seqPO,char *seqS,int len)
{
    VALIDATE(seqPO,SEQ_ID);
    if(len < 0) {
        len = strlen(seqS);
    }
    if(!AdjustSeqSizeI(seqPO, len+1, TRUE)) {
        return(FALSE);
    }
    strncpy(seqPO->seq,seqS,len);
    seqPO->len = len;
    seqPO->seq[seqPO->len] = '\0';
    return(TRUE);
}
/*****************************************************************************
*   Append sequence string to what's already in the SEQ struture
*/
int AppendSeqSequenceI(SEQ *seqPO,char *seqS,int len)
{
    int i;

    VALIDATE(seqPO,SEQ_ID);
    if(len < 0) {
        len = strlen(seqS);
    }
    if(!AdjustSeqSizeI(seqPO, seqPO->len + len, TRUE)) {
        return(FALSE);
    }
    for(i=0;i<len;i++)
    {
        BOG_CHECK( (seqPO->len + i) >= seqPO->ssize);
        seqPO->seq[seqPO->len + i] = seqS[i];
    }
    seqPO->len += len;
    seqPO->seq[seqPO->len] = '\0';
    return(TRUE);
}
/***************************************************************************
*   Append passed character to sequence, checking space as needed
*/
int AppendSeqCharI(SEQ *seqPO, char c, int error)
{
    VALIDATE(seqPO,SEQ_ID);
    if(seqPO->len >= seqPO->ssize) {
        if(!AdjustSeqSizeI(seqPO, seqPO->len + 1, error)) {
            return(FALSE);
        }
    }
    seqPO->seq[seqPO->len] = c;
    seqPO->len++;
    seqPO->seq[seqPO->len] = '\0';
    return(TRUE);
}
/*****************************************************************************
*   Set name of seq to nameS
*/
void SetSeqName(SEQ *seqPO,char *nameS)
{
    VALIDATE(seqPO,SEQ_ID);
    KillTrailStringBlanks(nameS);
    strncpy(seqPO->name,nameS,NSIZE-1);
    seqPO->name[NSIZE-1] = '\0';
}
/*****************************************************************************
*   Get length of seq
*/
int GetSeqLenI(SEQ *seqPO)
{
    VALIDATE(seqPO,SEQ_ID);
    return(seqPO->len);
}
/*****************************************************************************
*   Get number of SNPs in seq
*/
int GetSeqSnpCountI(SEQ *seqPO)
{
    VALIDATE(seqPO,SEQ_ID);
    return(seqPO->nsnp);
}
/*****************************************************************************
*   True if any ambiguous bases 
*/
int AnySeqAmbigsI(SEQ *seqPO)
{
    VALIDATE(seqPO,SEQ_ID);
    return(IS_SEQ_AMB(seqPO->flag));
}
/*****************************************************************************
*   Set pointer to sequence of seq
*/
int GetSeqSeqI(SEQ *seqPO, char **seqPPC)
{
    VALIDATE(seqPO,SEQ_ID);
    *seqPPC = seqPO->seq;
    return(TRUE);
}
/****************************************************************************
*   Fill passed string with sequence up to max
*/
int FillSeqSeqStringI(SEQ *seqPO,char *seqS,int max)
{
    VALIDATE(seqPO,SEQ_ID);
    if(max<0) {
        max = seqPO->len;
    }
    LIMIT_NUM(max,0,seqPO->len);
    max = MIN_NUM(max,seqPO->len);
    strncpy(seqS,seqPO->seq,max);
    seqS[max] = '\0';
    return(max);
}
/*****************************************************************************
*   Get name of seq to nameS
*/
int FillSeqNameStringI(SEQ *seqPO,char *nameS,int max)
{
    VALIDATE(seqPO,SEQ_ID);
    if(max<0) {
        max = strlen(seqPO->name);
    }
    LIMIT_NUM(max,0,NSIZE-1);
    if(nameS) {
        strncpy(nameS,seqPO->name,max);
        nameS[max] = '\0';
    }
    return(max);
}
/*****************************************************************************
*   Set SEQ sequnece to uppercase 
*/
int UppercaseSeqSeqI(SEQ *seqPO)
{
    VALIDATE(seqPO,SEQ_ID);
    Upperize(seqPO->seq);
    FinishSeqSettingsI(seqPO,FALSE,FALSE);
    return(TRUE);
}
/*****************************************************************************
*   Set SEQSEQ sequences to uppercase 
*/
int UppercaseSeqsetSeqsI(SEQSET *ssPO)
{
    int i;

    VALIDATE(ssPO,SEQSET_ID);
    for(i=0;i<ssPO->n;i++)
    {
        UppercaseSeqSeqI(ssPO->seqs[i]);
    }
    return(TRUE);
}
/*****************************************************************************
*   Replace SEQ sequnece with Reverse Complement
*/
int ReverseCompSeqSeqI(SEQ *seqPO)
{
    VALIDATE(seqPO,SEQ_ID);
    CompDNASeqI(seqPO->seq,seqPO->len,seqPO->seq);
    FinishSeqSettingsI(seqPO,FALSE,FALSE);
    return(TRUE);
}
/*****************************************************************************
*   Replace SEQSEQ sequences with Reverse Complements
*/
int ReverseCompSeqsetSeqsI(SEQSET *ssPO)
{
    int i;

    VALIDATE(ssPO,SEQSET_ID);
    for(i=0;i<ssPO->n;i++)
    {
        ReverseCompSeqSeqI(ssPO->seqs[i]);
    }
    return(TRUE);
}
