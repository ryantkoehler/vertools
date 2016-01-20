/*
* dna_splc.c
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
#include <ctype.h>
#include "prim.h"
#include "dna.h"


#define DB_SPLC if(DB[107])

/*************************************************************************
*   Merge two seqs together
*   Returns length of final merged seq
*   Sequence firS is slid onto secS
*   If dir = forward, firS -> secS to yield firS-secS
*   If dir = reverse, secS <- firS to yield secS-firS
*   upS is the "upstream" part of secS (or null) 
*/
int SpliceTwoSeqsI(char *firS, char *secS, char *upS, int dir, 
    int over, int ext, char *inS, char *newS)
{
    int len,flen,slen,ulen,ov,done;
    char ninS[DEF_BS];

    flen = strlen(firS);
    slen = strlen(secS);
    ulen = 0;
    if(upS)
    {
        ulen = strlen(upS);
    }
    DB_SPLC 
    {
        DB_PrI(">> SpliceTwoSeqsI dir=%d over=%d exten=%d\n",dir,over,ext);
        DB_PrI("+ fir |%s| %d\n",firS,flen);
        DB_PrI("+ sec |%s| %d\n",secS,slen);
        if(upS)
        {
            DB_PrI("+ up  |%s| %d\n",upS,ulen);
        }
    }
    /***
    *   Init things
    */
    INIT_S(ninS);
    ov = 0;
    done = FALSE;
    /***
    *   Try to overlap, unless it would require extension or exceed max
    */
    if(over)
    {
        ov = ProbeTargOverlapI(firS,strlen(firS),secS,strlen(secS),dir);
        DB_SPLC DB_PrI("+ ov=%d\n",ov);
        if( (ov<1) || (ov>over) || (NeedExtensionI(firS,ov,upS,dir,ext)) )
        {
            ov = 0;
        }
        else
        {
            if(dir<0)
            {
                sprintf(newS,"%s%s",secS,&firS[ov]);
            }
            else
            {
                sprintf(newS,"%s%s",firS,&secS[ov]);
            }
            done = TRUE;
        }
    }
    /***
    *   If we didn't get a good overlap already check for extension and build
    */
    DB_SPLC DB_PrI("+ ov=%d done=%d\n",ov,done);
    if(!done)
    {
        if(NeedExtensionI(firS,ov,upS,dir,ext))
        {
            GetExtensionSeqI(upS,dir,ext,ninS);
            DB_SPLC DB_PrI("+ extension |%s|\n",ninS);
        }
        if(dir<0)
        {
            sprintf(newS,"%s%s%s",secS,ninS,firS);
        }
        else
        {
            sprintf(newS,"%s%s%s",firS,ninS,secS);
        }
    }
    /***
    *   Copy insert if passed and return length of final seq
    */
    len = strlen(newS);
    if(inS)
    {
        strcpy(inS,ninS);
    }
    DB_SPLC DB_PrI("<< SpliceTwoSeqsI |%s| %d\n",newS,len);
    return(len);    
}
/*************************************************************************
*   Check if base extension is needed (to avoid extra complementarity)
*   Compare seqS to upS, looking for identities 
*   If dir = FORWARD
*       SSSSSSSSSSSSS
*                  UUUUU
*                  ^^ = overlap = 2 = ov 
*   else if dir = REVERSE
*       SSSSSSSSSSSSS
*    UUUU
*       ^ = overlap = 1 = ov 
*/
int NeedExtensionI(char *seqS,int ov, char *upS, int dir, int ext)
{
    int i,j,d,u,flen,ulen;

    DB_SPLC DB_PrI(">> NeedExtensionI ov=%d dir=%d\n",ov,dir);
    /***
    *   Not doing extension, no upstream string, empty upstream string
    */
    if( (ext<1) || (!upS) || (NO_S(upS)) )
    {
        DB_SPLC DB_PrI("<< NeedExtensionI FALSE nothing to do\n");
        return(FALSE);
    }
    flen = strlen(seqS);
    ulen = strlen(upS);
    /***
    *   Direction specific comparision
    */
    if(dir<0)
    {
        i = ov;
        j = 0;
        d = 1;
    }
    else
    {
        i = flen-ov-1;
        j = ulen-1;
        d = -1;
    }
    DB_SPLC DB_PrI("+ f|%s|=%d u|%s|=%d i=%d j=%d d=%d\n",
        seqS,flen,upS,ulen,i,j,d);
    /***
    *   Compare up to max or until we run out of sequence
    */
    for(u=0;u<ext;u++)
    {
        DB_SPLC DB_PrI("+  %d |%c| |%c|\n",u,seqS[i],upS[j]);
        if( UPPER(seqS[i]) == UPPER(upS[j]) )
        {
            DB_SPLC DB_PrI("<< NeedExtensionI TRUE\n");
            return(TRUE);
        }
        i += d;
        j += d;
        if( (i<0) || (i>=flen) || (j<0) || (j>=ulen) )
        {
            break;
        }
    }
    DB_SPLC DB_PrI("<< NeedExtensionI FALSE\n");
    return(FALSE);  
}
/*************************************************************************
*   Fill in extension seq
*   Assumes that what is passed as upS is the "upstream" part of a sequence from
*   the merge site and that this is the sense strand (not reverse complement).
*   Returns what should be inserted to "disrupt" the matching strand
*/
int GetExtensionSeqI(char *upS,int dir,int ext,char *inS)
{
    int i,j,d,ulen;

    DB_SPLC DB_PrI(">> GetExtensionSeqI up|%s| dir=%d ext=%d\n",upS,dir,ext);
    ulen = strlen(upS);
    if(dir<0)
    {
        d = 1;
        j = 0;
    }
    else
    {
        d = -1;
        j = ulen -1;
    }
    for(i=0;i<ext;i++)
    {
        switch(upS[j])
        {
            /***
            *   A implies T, so this is C-T mismatch
            */
            case 'A': case 'a': 
                inS[i] = 'C';   
                break;
            /***
            *   C implies G, so this is A-G mismatch
            */
            case 'C': case 'c': 
                inS[i] = 'A';   
                break;
            /***
            *   G implies C, so this is C-C mismatch
            */
            case 'G': case 'g': 
                inS[i] = 'C';   
                break;
            /***
            *   T implies A, so this is C-A mismatch
            */
            case 'T': case 't': 
                inS[i] = 'C';   
                break;
            default: 
                inS[i] = 'C';   break;
        }
        j += d;
        if( (j<0) || (j>=ulen) )
        {
            break;
        }
    }
    inS[i] = '\0';
    DB_SPLC DB_PrI("<< GetExtensionSeqI in|%s|\n",inS);
    return(i);
}
/*****************************************************************************
*   Find max overlap of probe onto target
*   dir governs direction of comparison:
*       If REVERSE, probe is slid backward (appended suffix) onto target
*       If FORWARD, probe is slid forwards (prepended prefix) onto target
*/
int ProbeTargOverlapI(char *probS,int plen,char *targS,int tlen,int dir)
{
    int i,j,nc,ok,max;

    DB_SPLC DB_PrI(">> ProbeTargOverlapI p=%d t=%d dir=%d",plen,tlen,dir);
    nc = MIN_NUM(plen,tlen);
    DB_SPLC 
    { 
        DB_PrI(" nc=%d\n",nc);
        WriteNiceSeqLine("+ Probe  ",probS,plen,outfileGPF);
        WriteNiceSeqLine("+ Target ",targS,tlen,outfileGPF);
    }
    /***
    *   Slide probe onto target 
    *   If any difference in an alignment, can't overlap 
    */
    max =  0;
    for(i=1;i<=nc;i++)
    {
        DB_SPLC DB_PrI("+ [%2d]",i);
        ok = TRUE;
        for(j=0;j<i;j++)
        {
            switch(dir)
            {
                case REVERSE:
                    if(UPPER(probS[j]) != UPPER(targS[tlen-i+j]))
                    {
                        DB_SPLC DB_PrI(" P[%2d]%c != T[%2d]%c\n",j,probS[j],
                            tlen-i+j,targS[tlen-i+j]);
                        ok = FALSE;
                    }
                    break;
                case FORWARD:
                    if(UPPER(probS[plen-i+j]) != UPPER(targS[j]))
                    {
                        DB_SPLC DB_PrI(" P[%2d]%c != T[%2d]%c\n",plen-i+j,
                            probS[plen-i+j],j,targS[j]);
                        ok = FALSE; 
                    }
                    break;
                default:
                    printf("Bogus direction code = %d\n",dir);
                    ERR("ProbeTargOverlapI","bad direction code");
            }
            if(!ok)
            {
                break;
            }
        }
        if(ok)
        {
            max = MAX_NUM(max,i);
            DB_SPLC DB_PrI(" max=%d\n",max);
        }
    }
    DB_SPLC DB_PrI("<< ProbeTargOverlapI %d\n",max);
    return(max);
}
/*****************************************************************************
*   Merge sequences; addPO is appended or prepended to seqPO based on how
*/
int MergeSeqsI(SEQ *seqPO,SEQ *addPO,int how)
{
    int nsize;

    VALIDATE(seqPO,SEQ_ID);
    VALIDATE(addPO,SEQ_ID);
    nsize = seqPO->len + addPO->len;
    if(!AdjustSeqSizeI(seqPO,nsize+1, TRUE))
    {
        printf("Can't allocate to merge together\n");
        return(FALSE);
    }
    /***
    *   Append is easy
    */
    if(how >= 0)
    {
        strcat(seqPO->seq,addPO->seq);
        seqPO->len = nsize;
    }
    /***
    *   Prepend is tricky; reverse both, append, reverse again
    */
    else
    {   
        ReverseString(seqPO->seq,seqPO->len,seqPO->seq);
        ReverseString(addPO->seq,addPO->len,addPO->seq);
        strcat(seqPO->seq,addPO->seq);
        seqPO->len = nsize;
        ReverseString(seqPO->seq,seqPO->len,seqPO->seq);
        ReverseString(addPO->seq,addPO->len,addPO->seq);
    }
    /***
    *   Merge counts and flags
    */
    seqPO->flag |= addPO->flag;
    seqPO->nsnp += addPO->nsnp;
    return(TRUE);
}
