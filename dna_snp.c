/*
* dna_snp.c
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
#include <ctype.h>
#include <string.h>
#include "prim.h"
#include "dna.h"

#define DB_SNPS     if(DB[101])

/***************************************************************************
*   Expand to an explicit sequence a given SNP/allele from a sequence
*
*   snp = snp record in sequence, i.e. the #'th "[X/Y]" block
*   al = allele for snp; i.e. for above, al=0 is "X", al=1 is "Y"
*   seqS = place to put answer (may be NULL)
*   up = the number of chars upstream from the SNP to (try to) get
*   max = max length of seq to get; i.e. answer goes from up to max
*   start = offset for extracted sequence relative to source sequence
*   targ = SNP coordinate in extracted sequence 
*   alS = allele sequence; i.e. for above, "X" or "Y" 
*   asize = len of allele; i.e. for above, both "X" & "Y" = 1
*
*   If OK, returns the coord of the snp start in seqS, else BOGUS
*/
int ExpandSeqSnpI(SEQ *seqPO, int snp, int al, char *seqS, int up, int max,
    int *startPI, int *targPI, char *alS, int *asizePI)
{
    int i,j,sn,ss,an,as,sp,n,del,asize1,asize2;

    DB_SNPS 
    {
        DB_PrI(">> ExpandSeqSnpI snp=%d al=%d %d %d\n",snp,al,up,max);
        DB_PrI("+ seqS=%p alS=%p asize=%p\n",seqS,alS,asizePI);
    }
    /***
    *   Init passed strings
    */
    if(seqS)
    {       
        INIT_S(seqS);  
    }
    if(alS)
    {       
        INIT_S(alS);  
    }
    if(seqPO->nsnp == 0)
    {
        DB_SNPS DB_PrI("<< ExpandSeqSnpI nsnp=%d BOGUS\n",seqPO->nsnp);
        return(BOGUS);
    }
    if( (snp<0) || (snp>=seqPO->nsnp) )
    {
        printf("***************** PROBLEM *****************\n");
        printf("Problem snp spec: %d n=%d\n",snp,seqPO->nsnp);
        ERR("ExpandSeqSnpI","Bogus SNP number");
        return(BOGUS);
    }
    /***
    *   Find SNP start position
    */
    ss = sn = 0;
    for(i=0;i<seqPO->len;i++)
    {
        if(seqPO->seq[i]=='[') 
        {
            sn++;
        }
        if(sn==(snp+1))
        {
            ss = i;
            break;
        }
    }
    if(sn < 1)
    {
        printf("***************** PROBLEM *****************\n");
        printf("Problem finding SNP %d in seq; n=%d\n",snp,seqPO->nsnp);
        ERR("ExpandSeqSnpI","Missing SNP start?");
        return(BOGUS);
    }
    DB_SNPS DB_PrI("+ snp start=%d\n",i);
    /***
    *   Back up from current SNP start
    */
    n = 0;
    i = ss;
    while( (i>0) && (n<up) )
    {
        /***
        *   If upstream SNP record, it expands to ambigs for larger allele
        *   For example [A/T] = N, [GG/-] = NN, [AGCC/TT] NNNN, etc
        */
        if(seqPO->seq[i]==']')
        {
/* ham */
            i--;
            del = asize1 = asize2 = 0;
            while( (seqPO->seq[i]!='[') && (i>=0) )
            {
                if(seqPO->seq[i]=='/')
                {
                    del++;
                }
                else
                {
                    if(del==0)
                    {
                        asize1++;
                    }
                    else
                    {
                        asize2++;
                    }
                }
                i--;
            }
            DB_SNPS DB_PrI("+  SNP@ %d sizes %d %d\n",i,asize1,asize2);
            n += MAX_NUM(asize1,asize2);
        }
        else
        {
            n++;
        }
        i--;
    }
    DB_SNPS DB_PrI("+ n=%d i = start=%d\n",n,i);
    if(startPI)
    {
        *startPI = i;
    }
    /***
    *   Copy source into supplied seq
    */
    n = an = as = 0;
    sp = BOGUS;
    while( (n<max) && (i<seqPO->len) )
    {
        /***
        *   SNP block?
        */
        if(seqPO->seq[i]=='[')
        {
            /***
            *   If this is the targeted snp, i = ss
            */
            if(i==ss)
            {
                i++;    
                /***
                *   Find allele in block and copy it
                */
                an++;
                while(seqPO->seq[i]!=']')
                {
                    if(i>=seqPO->len)
                    {
                        break;
                    }
                    if(seqPO->seq[i]=='/')
                    {
                        an++;
                        i++;    
                        continue;
                    }
                    if(an==(al+1))
                    {
                        if(IS_BOG(sp))
                        {
                            DB_SNPS DB_PrI("+ targ %d\n",n);
                            sp = n;
                        }
                        if(isalpha(INT(seqPO->seq[i])))
                        {
                            if(seqS)
                            {
                                seqS[n] = seqPO->seq[i];
                            }
                            if(alS)
                            {
                                alS[as] = seqPO->seq[i];
                            }
                            n++;
                            as++;
                        }
                    }
                    i++;
                }
                DB_SNPS DB_PrI("+ allele size %d\n",as);
            }
            /***
            *   Non-target SNP block, pad with Ns to max allele size
            */
            else
            {
                DB_SNPS DB_PrI("+  SNP@ %d",i);
                i++;
                del = 0;
                asize1 = asize2 = 0;
                while(seqPO->seq[i]!=']') 
                {
                    if(i>=seqPO->len)
                    {
                        break;
                    }
                    if(seqPO->seq[i]=='/')
                    {
                        del++;
                    }
                    else
                    {
                        if(del==0)
                            asize1++;
                        else
                            asize2++;
                    }
                    i++;
                }
                /***
                *   Non-target snp gets single ambiguous marking; then pass
                */
                asize1 = MAX_NUM(asize1,asize2);
                DB_SNPS DB_PrI(" padding %d\n",asize1);
                for(j=0;j<asize1;j++)
                {
                    if(seqS)
                    {
                        seqS[n] = 'N';
                    }
                    n++;
                }
            }
        }
        /***
        *   Normal seq copy
        */
        else
        {
            if(seqS)
            {
                seqS[n] = seqPO->seq[i];
            }
            n++; 
        }
        i++;
    }
    /***
    *   Did we get the intended snp allele?
    */
    if(IS_BOG(sp))
    {
        DB_SNPS DB_PrI("+ FAILED to find snp %d\n",snp);
        n = 0;
    }
    if(seqS)
    {
        seqS[n] = '\0';
        DB_SNPS DB_PrI("+ seq |%s|\n",seqS);
    }
    if(alS)
    {
        alS[as] = '\0';
        DB_SNPS DB_PrI("+ alS |%s|\n",alS);
    }
    if(targPI)
    {
        *targPI = sp;
    }
    if(asizePI)
    {
        *asizePI = as;
    }
    DB_SNPS DB_PrI("<< ExpandSeqSnpI %d\n",n);
    return(n);
}
/***************************************************************************
*
*/
int GetSeqsetSnpWindowI(SEQSET *seqsPO, int ind, int snp, int up, int down, 
    int max, char *seqS, int *sposPI)
{
    SEQ *snpPO;

    if(!GetSeqsetSeqI(seqsPO,ind,&snpPO))
    {
        return(FALSE);
    }
    return(GetSeqSnpWindowI(snpPO,snp,up,down,max,seqS,sposPI));
}
/***************************************************************************
*   Isolate window around SNP (i.e. up/down from [ ]) from sequence
*   Returns the size of string, or BOGUS if problem
*   Vars up and down indicate how much "upstream" and "downstream" sequence
*       is to be collected *AROUND* the SNP (so actual collected sequence is
*       bigger than up + down, as it also includes the SNP [X/Y]). 
*   Var max limits the total size of collected window (for above reason);
*/
int GetSeqSnpWindowI(SEQ *seqPO, int snp, int up, int down, int max, 
    char *seqS, int *sposPI)
{
    int i,j,sn,start,end,pos;
    char *seqPC;

    DB_SNPS DB_PrI(">> GetSeqSnpWindowI snp=%d up=%d down=%d\n",snp,up,down);
    VALIDATE(seqPO,SEQ_ID);
    if( (snp<1) || (snp>seqPO->nsnp) )
    {
        DB_SNPS DB_PrI("<< GetSeqSnpWindowI Seq %d SNPs BOGUS\n",seqPO->nsnp);
        return(BOGUS);
    }
    if(seqS)
    {   
        INIT_S(seqS);   
    }
    /***
    *   Find SNP 
    */
    sn = j = 0;
    for(i=0;i<seqPO->len;i++)
    {
        if(seqPO->seq[i]=='[') 
        {
            sn++;
        }
        if(sn==snp)
        {
            j = i;
            while(seqPO->seq[j]!=']')
            {
                if(j>=seqPO->len)
                {
                    j = BOGUS;
                    break;
                }
                j++;
            }
            break;
        }
    }
    DB_SNPS DB_PrI("+ i=%d j=%d\n",i,j);
    if( (sn!=snp) || (j<=i) )
    {
        DB_SNPS DB_PrI("<< GetSeqSnpWindowI didn't find SNP = BOGUS\n");
        return(BOGUS);
    }
    /***
    *   Figure out start, end, and SNP position
    */
    start = i - up;
    end = j + down + 1;
    LIMIT_NUM(start,0,seqPO->len-1);
    LIMIT_NUM(end,0,seqPO->len-1);
    pos = i-start;
    DB_SNPS DB_PrI("+ start=%d end=%d pos=%d\n",start,end,start);
    /***
    *   Set what's real
    */
    sn = end - start + 1;
    if(seqS)
    {
        if(!GetSeqSeqI(seqPO,&seqPC))
        {
            ERR("GetSeqSnpWindowI","Failed to get sequence string");
            return(BOGUS);
        }
        LIMIT_NUM(sn,0,max-1);
        strncpy(seqS,&seqPC[start],sn);
        seqS[sn] = '\0';
        DB_SNPS DB_PrI(" n=%d |%s|\n",sn,seqS);
    }
    if(sposPI)
    {       
        *sposPI = pos;
    }
    DB_SNPS DB_PrI("<< GetSeqSnpWindowI %d\n",sn);
    return(sn);
}
/***************************************************************************
*   Isolate SNP block (i.e. between [ ]) from sequence
*   If seqS is real, the block is copied there up to max chars
*   Returns the size of block string, or BOGUS
*/
int GetSeqSnpBlockI(SEQ *seqPO, int snp, char *seqS, int max)
{
    int n;

    DB_SNPS DB_PrI(">> GetSeqSnpBlockI snp=%d max=%d\n",snp,max);
    n = GetSeqSnpWindowI(seqPO, snp, 0, 0, max, seqS, NULL);
    DB_SNPS DB_PrI("<< GetSeqSnpBlockI %d\n",n);
    return(n);
}
/***************************************************************************
*   Isolates allele for specific snp from record sequence
*   If seqS is real, the allele seq is copied there up to max
*   Returns the size of filled string, or BOGUS
*/
int GetSeqSnpAlleleI(SEQ *seqPO, int snp, int al, char *seqS, int max)
{
    int i,sn,an,n;

    DB_SNPS DB_PrI(">> GetSeqSnpAlleleI snp=%d al=%d %d\n",snp,al,max);
    if(seqS)
    {   INIT_S(seqS);   }
    if( (snp<1) || (snp>seqPO->nsnp) )
    {
        printf("***************** PROBLEM *****************\n");
        printf("Problem snp spec: %d n=%d\n",snp,seqPO->nsnp);
        ERR("GetSeqSnpAlleleI","Bogus SNP number");
        return(BOGUS);
    }
    /***
    *   Find SNP 
    */
    n = an = sn = 0;
    for(i=0;i<seqPO->len;i++)
    {
        if(seqPO->seq[i]=='[') 
            sn++;
        if(sn==snp)
        {
            /***
            *   Copy specified allele of this snp block
            */
            i++;    
            an++;
            while(seqPO->seq[i]!=']')
            {
                if( (i>=seqPO->len) || (n>=max) )
                    break;
                if(seqPO->seq[i]=='/')
                {
                    an++;
                    i++;    
                    continue;
                }
                if(an==al)
                {
                    if(isalpha(INT(seqPO->seq[i])))
                    {
                        if(seqS)
                            seqS[n] = seqPO->seq[i];
                        n++;
                    }
                }
                i++;
            }
            if(seqS)
                seqS[n] = '\0';
            DB_SNPS DB_PrI("<< GetSeqSnpAlleleI %d\n",n);
            return(n);
        }
    }
    DB_SNPS DB_PrI("<< GetSeqSnpAlleleI BOGUS\n");
    return(BOGUS);
}
/***************************************************************************
*   Sham function to "clean up" badly formatted SNPs in seq
*/
int CleanUpSnpseqLineI(SEQ *seqPO)
{
    int i,j,n,olen; 
    char *oseqPC;

    /***
    *   Look for snps where there is no explicit delimiter and count these
    */
    i = n = 0;
    while(i<seqPO->len)
    {
        if(seqPO->seq[i]=='[')
        {
            if(!SnpHasDelimI(&seqPO->seq[i],seqPO->len-i))
                n++;
        }
        i++;
    }
    /***
    *   No expansion needed
    */
    if(n==0)
    {
        return(TRUE);
    }
    /***
    *   Will need to reallocate sequence space
    *   TOXIC; shouldn't mess with memory here
    */
    oseqPC = seqPO->seq;
    olen = seqPO->len;
    if( (olen+n) > seqPO->ssize )
    {
printf("Toxic cheese\n");
        return(FALSE);
        if(! (seqPO->seq = (char *)ALLOC(olen+n,sizeof(char)) ))
        {
            seqPO->seq = oseqPC;
            seqPO->len = olen;
        }
    }
    seqPO->len = olen+n;
    /***
    *   Loop again, this time copying
    */
    i = n = 0;
    while(i<olen)
    {
        if(oseqPC[i]=='[')
        {
            seqPO->seq[n++] = oseqPC[i++];
            if(!SnpHasDelimI(&oseqPC[i],olen-i))
            {
                j = 0;
                while( (i<olen) && (oseqPC[i]!=']') )
                {
                    seqPO->seq[n++] = oseqPC[i++];
                    if(j==0)
                        seqPO->seq[n++] = '/';
                }
            }
            else
            {
                while( (i<olen) && (oseqPC[i]!=']') )
                {
                    seqPO->seq[n++] = oseqPC[i++];
                }
            }
        }
        else
        {
            seqPO->seq[n++] = oseqPC[i++];
        }
    }
    /***
    *   Kill old seq
    */
    FREE(oseqPC);
    return(TRUE);
}
/**************************************************************************/
int SnpHasDelimI(char *seqS,int max)
{
    int i,del;

    i = del = 0;
    while( (i<max) && (seqS[i]!=']') )
    {
        if(seqS[i]=='/')
        {
            del++;
            break;
        }
        i++;
    }
    return(del);
}
/***************************************************************************
*   Expand IUPAC degenerate single-base SNP codes into explicit [X/Y] format
*/
int ExpandSeqSingBaseSNPsI(SEQ *seqPO)
{
    int i,j,n,slen,nlen;
    char snpS[DEF_BS];

    VALIDATE(seqPO,SEQ_ID);
    /***
    *   Figure out expanded length
    */
    nlen = 0;
    for(i=0;i<seqPO->len;i++)
    {
        n = ExpandDegBaseI(seqPO->seq[i],snpS);
        switch(n)
        {
            /***
            *   snpS has contents of [] brackets, so add 2 for these too
            *   n = 2 means [X/Y]
            *   n = 3 means [X/Y/Z]
            */
            case 2:
            case 3:
                nlen += 2;
                nlen += strlen(snpS);
                break;
            /***
            *   Just normal, non-change base
            */
            default:
                nlen++;
        }
    }
    /***
    *   Nothing to do?
    */
    if(nlen == seqPO->len)
    {
        return(TRUE);
    }
    /***
    *   Make room for new expanded SNP chars
    */
    slen = seqPO->len;
    if(!AdjustSeqSizeI(seqPO,nlen+1,TRUE))
    {
        return(FALSE);
    }
    seqPO->len = nlen;
    seqPO->seq[nlen] = '\0';
    /***
    *   Copy expansions backwards so don't walk on chars upstream
    */
    j = nlen - 1;
    for(i=(slen-1);i>=0;i--)
    {
        n = ExpandDegBaseI(seqPO->seq[i],snpS);
        switch(n)
        {
            case 2:
            case 3:
                n = strlen(snpS);
                seqPO->seq[j-n-1] = '[';
                strncpy(&seqPO->seq[j-n],snpS,n);
                seqPO->seq[j] = ']';
                j -= (n+2);
                break;
            default:
                seqPO->seq[j] = seqPO->seq[i];
                j--;
        }
    }
    /***
    *   Update final settings with new sequence
    */
    FinishSeqSettingsI(seqPO,FALSE,FALSE);
    return(TRUE);
}
