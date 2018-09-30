/*
* dna_out.c
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
#include "prim.h"
#include "dna.h"

#define DB_DNA_IO   if(DB[105])

/**********************************************************************/
void WriteSeqset(SEQSET *ssPO,int oform,int head,FILE *outPF)
{
    int i;
    char oformS[DEF_BS];

    DB_DNA_IO DB_PrI(">> WriteSeqset oform=%d head=%d\n",oform,head);
    if(head)
    {
        FillSeqFtypeDescString(oform,oformS);
        WriteSeqsetHeader(ssPO,ssPO->name,oformS,outPF);
    }
    for(i=0;i<ssPO->n;i++)
    {
        if(!ssPO->mask[i])
        {
            continue;   
        }
        WriteSeq(ssPO->seqs[i],oform,outPF);
    }
    DB_DNA_IO DB_PrI("<< WriteSeqset\n");
}
/**********************************************************************/
void WriteSeqsetHeader(SEQSET *ssPO,char *snameS,char *formS,FILE *outPF)
{
    int n,min,max,msnp;

    HAND_NFILE(outPF);
    if(snameS)
    {   fprintf(outPF,"# %s\n",snameS); }
    fprintf(outPF,"# %s format sequences\n",formS);
    TimeStamp("# ",outPF);
    SeqsetUnmaskDims(ssPO,&n,&min,&max,&msnp);
    fprintf(outPF,"# %d sequences\n",n);
    fprintf(outPF,"#   Lengths %d to %d\n",min,max);
    if(msnp)
    {   fprintf(outPF,"# Sequences have SNP records (%d max)\n",msnp); }
    fprintf(outPF,"#\n");
}
/***************************************************************************/
int WriteSeqsetSeqI(SEQSET *ssPO,int ind,int oform,FILE *outPF)
{
    VALIDATE(ssPO,SEQSET_ID);
    if( (ind < 0) || (ind >= ssPO->n) )
    {
        return(FALSE);
    }
    WriteSeq(ssPO->seqs[ind],oform,outPF);
    return(TRUE);
}
/***************************************************************************/
void WriteSeq(SEQ *seqPO,int oform,FILE *outPF)
{
    char nameS[NSIZE];

    DB_DNA_IO DB_PrI(">> WriteSeq seqPO=%p oform=%d outPF=%p\n",seqPO,oform,outPF);
    VALIDATE(seqPO,SEQ_ID);
    HAND_NFILE(outPF);
    switch(oform)
    {
        case SEQFM_FASTA:
            fprintf(outPF,">%s\n",seqPO->name);
            WriteMultiLineSeq(seqPO->seq,seqPO->len,FALSE,-1,-1,outPF);
            break;
        case SEQFM_NFAS:
            fprintf(outPF,">%s\n",seqPO->name);
            WriteMultiLineSeq(seqPO->seq,seqPO->len,TRUE,-1,-1,outPF);
            break;
        case SEQFM_RAW:
        default:
            FillSeqNameStringI(seqPO,nameS,NSIZE-1);
            ReplaceChars(' ',nameS,'_',nameS);
            DB_DNA_IO DB_PrI("+ name |%s|\n",nameS);
            fprintf(outPF,RAW_PFORM_S,nameS,seqPO->seq);
            break;
    }
    DB_DNA_IO DB_PrI("<<  WriteSeq\n");
}
/***************************************************************************/
void WriteMultiLineSeq(char *seqS,int slen, int nice, int line, int block, FILE *outPF)
{
    int i,j,n;

    line = (line<=0) ? NICE_LLEN: line;
    block = (block<=0) ? NICE_BLEN: block;
    HAND_NFILE(outPF);
    n = i = 0;
    while(i<slen)
    {
        if(seqS[i]=='[') {
            fputc(seqS[i++],outPF);
            j = 0;
            while(seqS[i]!=']') 
            {
                if(i>=slen) {
                    break;
                }
                j++;
                fputc(seqS[i++],outPF);
            }
            n++;
        }
        else {
            n++;
        }
        fputc(seqS[i],outPF);
        if( (nice) && (n) && ( (n % line) == 0) ) {
            fputc('\n',outPF);
        }
        else if( nice && n && ( (n % block) == 0) ) {
            fputc(' ',outPF);
        }
        i++;
    }
    fputc('\n',outPF);
}
/***************************************************************************
*   Dumps out sequence, nicely spaced for reading
*
*   If preS is a real string, output lines will be prepended with this
*/
void WriteNiceSeqLine(char *preS,char *seqS,int len,FILE *outPF)
{
    int n,nl;
    char bufS[DEF_BS];

    HAND_NFILE(outPF);
    if(preS)
    {   
        strncpy(bufS,preS,DEF_BS-1);    
        bufS[DEF_BS-1] = '\0';
    }
    else
    {   
        sprintf(bufS,"    ");   
    }
    n=nl=0;
    while(n<len)
    {
        if(!isgraph(INT(seqS[n])))
        {
            break;
        }
        if((n%NICE_LLEN)==0)
        {
/* XX sham 10/29/10 
            fprintf(outPF,bufS);
*/
            fprintf(outPF,"%s",bufS);
        }
        fputc(seqS[n],outPF);
        if(((n+1)%NICE_LLEN)==0)
        {
            fputc('\n',outPF);
            nl=1;
        }
        else if( n && (((n+1)%NICE_BLEN)==0) )
        {
            fputc(' ',outPF);
        }
        else
        {   
            nl=0;
        }
        n++;
    }
    if(!nl)
    {   fputc('\n',outPF);  }
}
/***************************************************************************/
void WriteSeqsetNameList(SEQSET *ssPO,FILE *oPF)
{
    int i;
    char nameS[NSIZE];

    HAND_NFILE(oPF);
    for(i=0;i<ssPO->n;i++)
    {
        if(ssPO->mask)
        {
            if(!ssPO->mask[i])
            {   continue;   }
        }
        FillSeqsetSeqNameI(ssPO,i,nameS,NSIZE-1);
        fprintf(oPF,"%s\n",nameS);
    }
}
