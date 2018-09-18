/*
* dna.c
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
#include <ctype.h>
#include <string.h>
#include "prim.h"
#include "dna.h"

#define DB_DNA      if(DB[100])

/***************************************************************************
*   Allocates a SEQ data structure (for single sequence)
*   If len is positive, will call to allocate this much space for sequence
*   If seqS is passed, the sequence is initialized with this sequence
*   If nameS is passed, the name is set
*   Seq space includes space for terminating NULL char (i.e. dim = [len+1] )
*/
SEQ *CreateSeqPO(int len, char *seqS, char *nameS)
{
    SEQ *seqPO;

    DB_DNA DB_PrI(">> CreateSeqPO len=%d seq=%p\n",len,seqS);
    if(! (seqPO = (SEQ *)ALLOC(1,sizeof(SEQ)))) {
        printf("Failed to allocate SEQ\n");
        return(NULL);
    }
    seqPO->ID = SEQ_ID;
    /***
    *   If no specified length and we've got a sequence, figure length
    */
    if( (len<0) && (seqS) ) {
        len = strlen(seqS);
    }
    if(len>0) {
        if(!AdjustSeqSizeI(seqPO, len+1, TRUE)) {
            CHECK_SEQ(seqPO);
            return(NULL);
        }
    }
    /***
    *   Add optionally passed initilizers
    */
    if(seqS) {
        SetSeqSequenceI(seqPO, seqS, len);   
    }
    if(nameS) {
        SetSeqName(seqPO, nameS);
    }
    DB_DNA DB_PrI("<< CreateSeqPO %p\n",seqPO);
    return(seqPO);
}
/***************************************************************************
*   Frees one or more SEQ data structures
*/
int DestroySeqI(SEQ *seqPO)
{
    VALIDATE(seqPO,SEQ_ID);
    CHECK_FREE(seqPO->seq);
    FREE(seqPO);
    return(TRUE);
}
/*****************************************************************************
*   If needed, allocate space for sequence up to size len + 1 (seq + null)
*/
int AdjustSeqSizeI(SEQ *seqPO, int len, int error)
{
    int n;

    VALIDATE(seqPO,SEQ_ID);
    DB_DNA DB_PrI(">> AdjustSeqSizeI len=%d ssize=%d\n",len,seqPO->ssize);
    /***
    *   If we don't have enough allocated size, we need to adjust 
    */
    if(len >= seqPO->ssize) {
        /***
        *   Allocate size = bigger of current + block or new length + 1
        */
        n = MAX_NUM(seqPO->ssize + SEQ_BLOCK, len + 1);
        if(seqPO->ssize==0) {
            DB_DNA DB_PrI("+ allocating size=%d\n",n);
            seqPO->seq = (char *)ALLOC(n,sizeof(char));
        }
        else {
            DB_DNA DB_PrI("+ re-allocating %p, size=%d\n",seqPO->seq,n);
            seqPO->seq = (char *)REALLOC(seqPO->seq,n,sizeof(char));
        }
        if(!seqPO->seq) {
            if(error) {
                printf("PROBLEM ALLOCATING FOR SEQ len %d (size %d)\n",len,n);
            }
            return(FALSE);
        }
        seqPO->ssize = n;
    }
    DB_DNA DB_PrI("<< AdjustSeqSizeI seq->seq=%p TRUE\n",seqPO->seq);
    return(TRUE);
}
/***************************************************************************
*   Report contents of SEQ data structure
*/
void DumpSeq(SEQ *seqPO,FILE *outPF)
{
    VALIDATE(seqPO,SEQ_ID);
    HAND_NFILE(outPF);
    fprintf(outPF, LINEBAR_S);
    fprintf(outPF,"Seq[%d] %s\n",seqPO->ind,seqPO->name);
    fprintf(outPF,"    Len=%d  \n",seqPO->len);
    fprintf(outPF,"    Ad=%p   Flag=%x\n",seqPO,seqPO->flag);
    fprintf(outPF,"    |%s|\n",seqPO->seq);
    fprintf(outPF, LINEBAR_S);
}
/***************************************************************************
*   Init contents of SEQ
*   If par, reference to parent structure is cleared
*   If mem, memory is cleared 
*/
void InitSeq(SEQ *seqPO, int par, int mem)
{
    VALIDATE(seqPO,SEQ_ID);
    INIT_S(seqPO->name);
    if(par) {
        seqPO->ind = 0;
        seqPO->par = NULL;
    }
    seqPO->len = 0;
    seqPO->nsnp = 0;
    seqPO->flag = 0;
    if(mem) {
        CHECK_FREE(seqPO->seq);
        seqPO->ssize = 0;
    }
}
/***************************************************************************
*   Copy contents of fseqPO into sseqPO, possibly calling to adjust it's size
*   Doesn't copy parent pointer / index, as these likely differ
*   If non-neg, only the sequence from st to length len will be copied
*/
int CopySeqI(SEQ *fseqPO, SEQ *sseqPO, int st, int len)
{
    DB_DNA DB_PrI(">> CopySeqI %p %p st=%d len=%d\n",fseqPO,sseqPO,st,len);
    VALIDATE(fseqPO,SEQ_ID);
    VALIDATE(sseqPO,SEQ_ID);
    /***
    *   Figure range of sequence to copy
    */
    if(st<0) {
        st = 0;
        len = fseqPO->len;
    }
    else {
        if( (len<1) || ((st+len) > fseqPO->len) ) {
            printf("Bad subseq dimensions: st=%d len=%d seqlen=%d\n",
                st,len,fseqPO->len);
            return(FALSE);
        }
    }
    /***
    *   Size adjustment; this handle's setting sseqPO->ssize
    */
    if(!AdjustSeqSizeI(sseqPO, fseqPO->ssize, TRUE)) {
        printf("Can't copy if can't allocate!\n");
        return(FALSE);
    }
    strcpy(sseqPO->name,fseqPO->name);
    strncpy(sseqPO->seq,&fseqPO->seq[st],len);
    sseqPO->seq[len] = '\0';
    sseqPO->len = len;
    sseqPO->nsnp = fseqPO->nsnp;
    sseqPO->flag = fseqPO->flag;
    DB_DNA DB_PrI("<< CopySeqI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Shrinks a sequence from start to len; 
*   dir = direction (end of string) from which to narrow 
*   fits = flag; if true, dims must fit withing passed seq, else fails
*/
int NarrowSeqI(SEQ *seqPO, int start, int len, int dir, int fits)
{
    int i,f,n;

    DB_DNA DB_PrI(">> NarrowSeqI %p start=%d len=%d dir=%d fits=%d\n",seqPO,
        start,len,dir,fits);
    VALIDATE(seqPO,SEQ_ID);
    if( (start<0) || (len<0) || (start>=seqPO->len) || (len>seqPO->len) ) {
        DB_DNA DB_PrI("+ out of bounds; seq-len=%d\n",seqPO->len);
        if(fits==TRUE) {
            DB_DNA DB_PrI("<< NarrowSeqI FALSE\n");
            return(FALSE);
        }
    }
    if(dir<0) {
        f = seqPO->len - len - start;
    }
    else {
        f = start;
    }
    DB_DNA DB_PrI("+ start=%d len=%d f=%d\n",start,len,f);
    n = 0;
    for(i=0; i<len; i++)
    {
        if( (f >=0 ) && (f < seqPO->len) ) {
            seqPO->seq[n++] = seqPO->seq[f];
        }
        f++;
    }
    seqPO->seq[n] = '\0';
    seqPO->len = n;
    DB_DNA DB_PrI("<< NarrowSeqI slen=%d TRUE\n",seqPO->len);
    return(TRUE);
}
/***************************************************************************
*   Allocates a SEQSET data structure
*   If len is positive, will allocate this much space for sequence
*/
SEQSET *CreateSeqsetPO(int n)
{
    SEQSET *ssPO;

    DB_DNA DB_PrI(">> CreateSeqsetPO n=%d\n",n);
    if(! (ssPO = (SEQSET *)ALLOC(1,sizeof(SEQSET)))) {
        printf("Failed to allocate SEQSET\n");
        return(NULL);
    }
    ssPO->ID = SEQSET_ID;
    if (n>0) {
        if(!AdjustSeqsetSizeI(ssPO,n)) {
            CHECK_SEQSET(ssPO);
            return(NULL);
        }
    }
    DB_DNA DB_PrI("<< CreateSeqsetPO %p\n",ssPO);
    return(ssPO);
}
/***************************************************************************
*   Frees SEQSET data structures
*/
int DestroySeqsetI(SEQSET *ssPO)
{
    int i;

    VALIDATE(ssPO,SEQSET_ID);
    for(i=0;i<ssPO->n;i++)
    {
        CHECK_SEQ(ssPO->seqs[i]);
    }
    CHECK_FREE(ssPO->seqs);
    CHECK_FREE(ssPO->mask);
    FREE(ssPO);
    return(TRUE);
}
/*****************************************************************************
*   If needed, allocate space to hold up to size
*   If this is the first call, initialize mask so every seq is ok
*/
int AdjustSeqsetSizeI(SEQSET *ssPO,int size)
{
    DB_DNA DB_PrI(">> AdjustSeqsetSizeI size=%d\n",size);
    VALIDATE(ssPO,SEQSET_ID);
    if (size >= ssPO->nsize) {
        size = ssPO->nsize + SSN_BLOCK;
        if (ssPO->nsize==0) {
            ssPO->seqs = (SEQ **)ALLOC(size,sizeof(SEQ *));
            ssPO->mask = (char *)ALLOC(size,sizeof(char));
            if(ssPO->mask) {
                InitArrayI(ssPO->mask,IS_CHAR,0,size,TRUE);
            }
        }
        else {
            ssPO->seqs = (SEQ **)REALLOC(ssPO->seqs,size,sizeof(SEQ *));
            ssPO->mask = (char *)REALLOC(ssPO->mask,size,sizeof(char));
        }
        if( (!ssPO->seqs) || (!ssPO->mask) ) {
            printf("PROBLEM ALLOCATING FOR SEQSET array (size %d)\n",size);
            return(FALSE);
        }
        ssPO->nsize = size;
    }
    DB_DNA DB_PrI("<< AdjustSeqsetSizeI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Dump SEQSET data structure
*/
void DumpSeqset(SEQSET *ssPO,int all,FILE *outPF)
{
    int i;
    
    HAND_NFILE(outPF);
    VALIDATE(ssPO,SEQSET_ID);
    fprintf(outPF, LINEBAR_S);
    fprintf(outPF,"SEQSET @  %p\n",ssPO);
    fprintf(outPF,"Name      %s\n",ssPO->name);
    fprintf(outPF,"Source    %s\n",ssPO->source);
    fprintf(outPF,"Type      %d\n",ssPO->type);
    fprintf(outPF,"N         %d\n",ssPO->n);
    fprintf(outPF,"Mask @    %p\n",ssPO->mask);
    if(all) {
        for(i=0;i<ssPO->n;i++)
        {
            DumpSeq(ssPO->seqs[i],outPF);
        }
    }
    fprintf(outPF, LINEBAR_S);
}
/*****************************************************************************
*   Add sequence specified by char string and name to seq set
*   Returns TRUE if OK, FALSE otherwise
*/
int AddNamedSequenceToSeqsetI(char *seqS, char *nameS, SEQSET *ssPO,
    int *indPI)
{
    SEQ *seqPO;

    VALIDATE(ssPO,SEQSET_ID);
    if(!(seqPO=CreateSeqPO(-1, seqS, nameS)) ) {
        return(FALSE);
    }
    /***
    *   Set datastructure bits and flags based on seqS sequence
    *   No "cleaning" of supplied seq and no error reporting
    */
    FinishSeqSettingsI(seqPO,FALSE,FALSE);
    /***
    *   Stick this guy into the set
    */
    if(!AddSeqToSeqsetI(seqPO,ssPO)) {
        CHECK_SEQ(seqPO);
        return(FALSE);
    }
    if(indPI) {
        *indPI = seqPO->ind;
    }
    return(TRUE);
}
/*****************************************************************************
*   Attempt to add sequence to sequence set
*/
int AddSeqToSeqsetI(SEQ *seqPO,SEQSET *ssPO)
{
    VALIDATE(ssPO,SEQSET_ID);
    VALIDATE(seqPO,SEQ_ID);
    /***
    *   Make sure there is room for this guy
    */
    if(!AdjustSeqsetSizeI(ssPO,ssPO->n)) {
        printf("Failed to adjust size for %d sequences\n",ssPO->n);
        return(FALSE);
    }
    /***
    *   Add to array, set mask and set pointers
    */
    ssPO->seqs[ssPO->n] = seqPO;
    ssPO->mask[ssPO->n] = TRUE;
    seqPO->par = ssPO;
    seqPO->ind = ssPO->n;
    ssPO->n += 1;
    return(TRUE);   
}
/***************************************************************************
*   Init contents of seqset. 
*/
void InitSeqset(SEQSET *ssPO)
{
    VALIDATE(ssPO,SEQSET_ID);
    INIT_S(ssPO->name);
}
/***************************************************************************
*   Sets dimension numbers of unmasked seqs;
*   nPI = number (not masked)
*   sPI = min (small) length (of not masked)
*   mPI = max length (of not masked)
*   snPI = max snp count (of not masked)
*/
void SeqsetUnmaskDims(SEQSET *ssPO,int *nPI,int *sPI,int *mPI, int *snPI)
{
    int i,n,s,m,sn,len,nsnp;

    VALIDATE(ssPO,SEQSET_ID);
    s = m = n = sn = 0;
    for(i=0;i<ssPO->n;i++)
    {
        if(!ssPO->mask[i]) {
            continue;
        }
        len = GetSeqsetSeqLenI(ssPO,i);
        m = MAX_NUM(m,len);
        s = MIN_NUM(m,len);
        nsnp = GetSeqsetSeqSnpCountI(ssPO,i);
        sn = MAX_NUM(sn,nsnp);
        n++;
    }
    if(nPI)
    { *nPI = n; }
    if(mPI)
    { *mPI = m; }
    if(sPI)
    { *sPI = s; }
    if(snPI)
    { *snPI = sn; }
}
/****************************************************************************
*   Sets bit flags indicative of what info / chars are in seq
*   Also applies "cleaning" to non-standard bases 
*
*   Returns:
*   TRUE if all good
*   FALSE if problem (warnings)
*   BOGUS if problem (error)
*/
int FinishSeqSettingsI(SEQ *seqPO, int clean, int error)
{
    int i,j,nsnp,ok,insnp;

    VALIDATE(seqPO,SEQ_ID);
    ok = TRUE;
    /***
    *   If no name, set default
    */
    if(FillSeqNameStringI(seqPO,NULL,-1) < 1) {
        SetSeqName(seqPO,DEF_SEQ_NAME_S);
    }
    /***
    *   Clear flag then set based on what we have
    */
    seqPO->flag = 0;
    /***
    *   Scan sequence and set flags based on what's found
    */
    insnp = nsnp = j = 0;
    for(i=0;i<seqPO->len;i++)
    {
        switch(seqPO->seq[i]) {
            /***
            *    Normal ACGT & U
            */
            case 'A': case 'a':
            case 'C': case 'c':
            case 'G': case 'g':
            case 'T': case 't':
                break;
            case 'U': case 'u':
                SET_SEQ_HASU(seqPO->flag);
                break;
            /***
            *   2-fold degen; i.e. S = G/C, W = A/T
            *   3-fold degen; i.e. B= not A = C/G/T
            */
            case 'S': case 's':
            case 'W': case 'w':
            case 'M': case 'm':
            case 'K': case 'k':
            case 'R': case 'r':
            case 'Y': case 'y':
            case 'B': case 'b':
            case 'D': case 'd':
            case 'H': case 'h':
            case 'V': case 'v':
                SET_SEQ_AMB(seqPO->flag);
                if (clean > SCLEAN_MID) {
                    if(error) {
                        printf("\n");
                        printf("# WARNING: Bad seq char |%c| @%d in %s\n",
                            seqPO->seq[i],i+1,seqPO->name);
                        printf("#   Setting to \"N\"\n");
                    }
                    seqPO->seq[i] = 'N';
                }
                else {
                    SET_SEQ_DEG(seqPO->flag);
                }
                break;
            /***
            *   4 = N
            */
            case 'N': case 'n':
                SET_SEQ_AMB(seqPO->flag);
                break;
            /***
            *   Special annotations; SNPs, deletions ...
            */
            case '[': 
                if(insnp) {
                    if(error) {
                        printf("\n");
                        printf("# WARNING: second \"[\" in SNP @%d in %s\n",
                            i+1,seqPO->name);
                    }
                    ok = FALSE;
                }
                insnp = TRUE; 
                break;
            case ']': 
                if(!insnp) {
                    if(error) {
                        printf("\n");
                        printf("# WARNING: \"]\" outside of SNP @%d in %s\n",
                            i+1,seqPO->name);
                    }
                    ok = FALSE;
                }
                insnp = FALSE; 
                nsnp++;
                break;
            case '/': 
                if(!insnp) {
                    SET_SEQ_COAX(seqPO->flag);
                }
                break;
            case '*':
                SET_SEQ_DEL(seqPO->flag);
                break;
            case '-':
                SET_SEQ_INS(seqPO->flag);
                break;
            /***
            *   Default = some bogus letter?
            */
            default:
                /* Not letter = sham */
                if( ! isalpha(INT(seqPO->seq[i])) ) {
                    if(error) {
                        printf("# ERROR: Bad seq char |%c| @%d in %s\n",
                            seqPO->seq[i],i+1,seqPO->name);
                    }
                    ok = BOGUS;
                    break;
                }
                /* Letter ... but not covered above */
                SET_SEQ_AMB(seqPO->flag);
                if(error) {
                    printf("\n");
                    printf("# WARNING: Bad seq char |%c| @%d in %s\n",
                        seqPO->seq[i],i+1,seqPO->name);
                    printf("#   Setting to \"N\"\n");
                }
                seqPO->seq[i] = 'N';
                if ( clean < SCLEAN_MID) {
                    SET_SEQ_NS(seqPO->flag);
                    break;
                }
                ok = FALSE;
                break;
        }
        /***
        *   Already a problem, so don't keep checking
        */
        if(ok != TRUE) {
            break;
        }
    }
    /***
    *   Remember any SNPs
    */
    seqPO->nsnp = nsnp;
    if(seqPO->nsnp>0) {
        SET_SEQ_SNP(seqPO->flag);
    }
    return(ok);
}
/**************************************************************************
*   Allocate SEQCOMP
*/
SEQCOMP *CreateSeqcompPO()
{
    SEQCOMP *scPO;

    if(! (scPO = (SEQCOMP *)ALLOC(1,sizeof(SEQCOMP)))) {
        printf("Failed to allocate SEQCOMP\n");
        return(NULL);
    }
    scPO->ID = SEQCOMP_ID;
    InitSeqcomp(scPO);
    return(scPO);
}
/***************************************************************************
*   Free one SEQCOMP data structure
*/
int DestroySeqcompI(SEQCOMP *scPO)
{
    VALIDATE(scPO,SEQCOMP_ID);
    FREE(scPO);
    return(TRUE);
}
/***************************************************************************
*   Init contents of SEQCOMP
*/
void InitSeqcomp(SEQCOMP *scPO)
{
    VALIDATE(scPO,SEQCOMP_ID);
    scPO->slen = scPO->nbase = 0;
    scPO->ra = scPO->rc = scPO->rg = scPO->rt = 0;
    scPO->rs = scPO->rw = scPO->rr = scPO->ry = scPO->rk = scPO->rm = 0;
    scPO->na = scPO->nc = scPO->ng = scPO->nt = 0;
    scPO->fa = scPO->fc = scPO->fg = scPO->ft = 0.0;
    scPO->n_dinuc = 0;
    InitArrayI(scPO->dinuc, IS_INT, 0, 16, 0);
}
/***************************************************************************
*   Init contents of SEQCOMP
*/
void DumpSeqcomp(SEQCOMP *scPO,FILE *outPF)
{
    int i;

    VALIDATE(scPO,SEQCOMP_ID);
    HAND_NFILE(outPF);
    fprintf(outPF, LINEBAR_S);
    fprintf(outPF,"SEQCOMP slen=%d\n",scPO->slen);
    fprintf(outPF,"  nA=%-5d nC=%-5d nG=%-5d nT=%-5d\n", 
        scPO->na,scPO->nc,scPO->ng,scPO->nt);
    fprintf(outPF,"  fA=%5.3f fC=%5.3f fG=%5.3f fT=%5.3f\n", 
        scPO->fa,scPO->fc,scPO->fg,scPO->ft);
    fprintf(outPF,"  rA=%-5d rC=%-5d rG=%-5d rT=%-5d\n", 
        scPO->ra,scPO->rc,scPO->rg,scPO->rt);
    fprintf(outPF,"  rS=%-5d rW=%-5d rR=%-5d rY=%-5d rK=%-5d rM=%-5d\n", 
        scPO->rs,scPO->rw,scPO->rr,scPO->ry,scPO->rk,scPO->rm);
    fprintf(outPF,"  n_dinuc = %d\n",scPO->n_dinuc);
    fprintf(outPF,"  Dinucs:");
    for(i=0;i<16;i++) { 
        fprintf(outPF,"\t%d",scPO->dinuc[i]);
    }
    fprintf(outPF,"\n");
    fprintf(outPF, LINEBAR_S);
}
/**************************************************************************
*   Allocate SEQTRIM
*/
SEQTRIM *CreateSeqtrimPO()
{
    SEQTRIM *stPO;

    if(! (stPO = (SEQTRIM *)ALLOC(1,sizeof(SEQTRIM)))) {
        printf("Failed to allocate SEQTRIM\n");
        return(NULL);
    }
    stPO->ID = SEQTRIM_ID;
    InitSeqtrim(stPO);
    return(stPO);
}
/***************************************************************************
*   Free one SEQCOMP data structure
*/
int DestroySeqtrimI(SEQTRIM *stPO)
{
    VALIDATE(stPO,SEQTRIM_ID);
    FREE(stPO);
    return(TRUE);
}
/***************************************************************************
*   Init SEQTRIM
*/
void InitSeqtrim(SEQTRIM *stPO)
{
    VALIDATE(stPO,SEQTRIM_ID);
    stPO->base_s = stPO->base_e = -1;
    stPO->basf_s = stPO->basf_e = -1.0;
    stPO->trim_s = stPO->trim_e = -1;
    stPO->wind = -1;
    stPO->winf = -1.0;
    stPO->wcent = -1;
    stPO->wcenf = -1.0;
    stPO->rre = FALSE;
    stPO->nmask = FALSE;
    stPO->umask = FALSE;
    stPO->lcase = FALSE;
    stPO->ucase = FALSE;
    stPO->one_base = TRUE;
    stPO->end_in = TRUE;
}
