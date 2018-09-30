/*
* dna_pair.c
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
#include <math.h>
#include "prim.h"
#include "dna.h"
#include "dna_pair.h"

#define DB_PAIR     if(DB[116])
#define DB_PAIRLO   if(DB[117])

/*****************************************************************************
*   Compare single sequence to all others in the sequence set
*
*   fPC is the sequence, of length flen, to be compared
*   ssPO holds sequences to compare to 
*   pparPO holds pair-wise comparision parameters
*   Scores are set into scoresPR array
*   If offset array offsPI is passed, offsets cooresponding to scores go here
*/
int CompareSeqToSeqsetI(SEQ *qseqPO, SEQSET *ssPO, PPARS *pparPO,
    REAL *scoresPR, int *offsPI)
{
    int j,off,flen,slen;
    char *fPC,*sPC;
    REAL scR;

    VALIDATE(qseqPO,SEQ_ID);
    VALIDATE(ssPO,SEQSET_ID);
    GetSeqSeqI(qseqPO,&fPC);
    flen = GetSeqLenI(qseqPO);
    DB_PAIR
    {
        DB_PrI(">> CompareSeqToSeqsetI scores=%p offs=%p\n",scoresPR,offsPI);
        DB_PrI("+ Query |%s| %d\n",fPC,flen);
    }
    /***
    *   Compare to each seq in set
    */
    for(j=0;j<ssPO->n;j++)
    {
        GetSeqsetSeqStringI(ssPO,j,&sPC);
        slen = GetSeqsetSeqLenI(ssPO,j);
        DB_PAIR 
        {   
            DB_PrI("+ [%d] |",j); PrintString(sPC,slen,NULL);
            DB_PrI("| %d\n",slen);
        }
        off = ScoreSeqCompareI(fPC,flen,sPC,slen,pparPO,&scR);
        scoresPR[j] = scR;
        if(offsPI) {
            offsPI[j] = off;
        }
        DB_PAIR DB_PrI("+ [%d]=%f off=%d\n",j,scoresPR[j],off);
    }
    DB_PAIR DB_PrI("<< CompareSeqToSeqsetI\n");
    return(ssPO->n);
}
/****************************************************************************
*   Compare first seq fS[flen] to second seq sS[slen] with pars in ppPO
*   If loop's are to be evaluated, sS is ignored (may be NULL)
*   Return value is offset (coord of second on first) of score for comparison
*   Score is returned via pointer scPR
*/
int ScoreSeqCompareI(char *fS, int flen, char *sS, int slen, PPARS *ppPO,
    REAL *scPR)
{
    REAL scR;
    int off;

    DB_PAIR 
    {   
        DB_PrI(">> ScoreSeqCompareI flen=%d slen=%d\n",flen,slen);
        DB_PrI("+ Compare type %d\n",ppPO->ctype);
        DB_PrI("+  fS|"); PrintString(fS,flen,NULL); DB_PrI("| %d\n",flen); 
        if(sS)
        {
            DB_PrI("+  sS|"); PrintString(sS,slen,NULL); DB_PrI("| %d\n",slen); 
        }
    }
    off = 0;
    scR = 0.0;
    /***
    * Hamming dist; Set both lengths to min of either; "minword" to full len
    */
    if(ppPO->do_ham) {
        flen = slen = MIN_NUM(flen,slen);
        SetPparsWordMinI(ppPO,flen);
    }
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
        case PP_COM_PCOM:
            off = ScoreSeqParallelI(fS,flen,sS,slen,ppPO,&scR);
            break;
        case PP_COM_COM:
            off = ScoreSeqAntiParallelI(fS,flen,sS,slen,ppPO,&scR);
            break;
        case PP_COM_LOOP:
            off = ScoreLoopOverlapI(fS,flen,ppPO,&scR);
            break;
        default:
            printf("Bogus compare code = %d\n",ppPO->ctype);
            ERR("ScoreSeqCompareI","Bad code");
            return(FALSE);
    }
    *scPR = scR;
    DB_PAIR DB_PrI("<< ScoreSeqCompareI %4.3f off\n",*scPR);
    return(off);
}
/****************************************************************************
*   Evaluate pairwise similarity of fS[flen] to sS[slen] via all alignments
*   Return offset (first onto second) for highest scoring alignment
*/
int ScoreSeqParallelI(char *fS,int flen,char *sS,int slen,PPARS *ppPO, REAL *scPR)
{
    int i,f,s,m,n,find,mind,minlen,mover,best;
    REAL scR,bestR;

    DB_PAIR DB_PrI(">> ScoreSeqParallelI\n");
    /***
    *   Figure out how many comparisons need to be made (m)
    *   find = index to figure starts for each alignment
    */
    minlen = MIN_NUM(flen,slen);
    m = flen + slen -(2 * ppPO->min_word) + 1;
    mover = m + ppPO->min_word - 1;
    find = flen - ppPO->min_word;
    DB_PAIR DB_PrI("+ min_word=%d m=%d minlen=%d find=%d\n",
        ppPO->min_word,m,minlen,find);
    /***
    *   Initialize based on what will be returned
    */
    switch(ppPO->rval)
    {
        case PP_R_MAX:  
            bestR = -TOO_BIG_R; 
            break;
        case PP_R_TOT:
        case PP_R_NUM:
            bestR = 0.0;    
            break;
        default:
            printf("Bad rval code = %d\n",ppPO->rval);
            ERR("ScoreSeqParallelI","Bogus rval code");
            return(FALSE);
    }
    best = 0;
    /***
    *   Try each of m comparisons, remembering best score and offset of this
    */
    for(i=0;i<m;i++)
    {
        f = MAX_NUM(0,find-i);
        s = MAX_NUM(0,i-find);
        mind = MIN_NUM(ppPO->min_word+i, mover-i);
        n = MIN_NUM(minlen,mind);
        DB_PAIR DB_PrI("\n+ M[%d] n=%d f=%d s=%d mind=%d\n",i,n,f,s,mind);
        /***
        *   If clamps, make sure these are ok before regular comparision
        */
        if( (ppPO->cl3>0) || (ppPO->alg==ALGO_CON3) )
        {
            if( (f+n) < flen) 
            {
                DB_PAIR DB_PrI("+  PAST 3' end, can't clamp = BREAK\n");
                break;
            }
            if( (flen - ppPO->cl3 - f) < 0) 
            {
                DB_PAIR DB_PrI("+  Too short 3', can't clamp %d = next\n",
                    ppPO->cl3);
                continue;
            }
            if(!OkAlignClampI(ppPO,fS,flen-ppPO->cl3, sS,s+n-ppPO->cl3,ppPO->cl3))
            {
                DB_PAIR DB_PrI("+  No 3' clamp = next\n");
                continue;
            }
        }
        if( (ppPO->cl5>0) || (ppPO->alg==ALGO_CON5) )
        {
            if(f > 0) 
            {
                DB_PAIR DB_PrI("+  past 5' end, can't clamp = next\n");
                continue;
            }
            if(n < ppPO->cl5) 
            {
                DB_PAIR DB_PrI("+  Too short 5' end, can't clamp = BREAK\n");
                continue;
            }
            if(!OkAlignClampI(ppPO,fS,0,sS,s,ppPO->cl5))
            {
                DB_PAIR DB_PrI("+  No 5' clamp = next\n");
                continue;
            }   
        }
        /***
        *   Compare region of n letters starting from fS[f] and sS[s]
        *   Always save the "best" offset, regardless of score returned
        */
        scR = AlignedSeqScoreR(ppPO,fS,f,sS,s,n);
        if(scR > bestR)
        {
            best = find-i;
        }
        /***
        *   What to do with score?
        */
        switch(ppPO->rval)
        {
            case PP_R_MAX:
                if(scR > bestR)
                {
                    bestR = scR;
                    DB_PAIR DB_PrI("+  Best score = %3.3f @%d\n",bestR,best);
                }   
                break;
            case PP_R_NUM:
                if(scR >= ppPO->thresh)
                {
                    bestR += 1.0;
                    DB_PAIR DB_PrI("+  Best count = %2.0f @%d\n",bestR,best);
                }
                break;
            case PP_R_TOT:
                if(scR >= ppPO->thresh)
                {
                    bestR += scR;
                    DB_PAIR DB_PrI("+  Best tot = %3.3f @%d\n",bestR,best);
                }
                break;
        }
    }
    *scPR = bestR;
    DB_PAIR DB_PrI("<< ScoreSeqParallelI %d (%4.2f)\n",best,bestR);
    return(best);
}
/****************************************************************************
*   Score sequence complimentarity; both fS[flen] and sS[slen] are 5'>-->3'
*/
int ScoreSeqAntiParallelI(char *fS,int flen,char *sS,int slen,PPARS *ppPO, 
    REAL *scPR)
{
    int i,f,s,m,n,find,mind,minlen,mover,best;
    REAL scR,bestR;

    DB_PAIR DB_PrI(">> ScoreSeqAntiParallelI\n");
    /***
    *   Figure out how many comparisons need to be made (m)
    *   find = index to figure starts for each alignment
    */
    m = flen + slen -(2 * ppPO->min_word) + 1;
    mover = m + ppPO->min_word - 1;
    minlen = MIN_NUM(flen,slen);
    find = flen - ppPO->min_word;
    DB_PAIR DB_PrI("+ min_word=%d m=%d minlen=%d find=%d\n",
        ppPO->min_word,m,minlen,find);
    /***
    *   Initialize based on what will be returned
    */
    switch(ppPO->rval)
    {
        case PP_R_MAX:  
            bestR = -TOO_BIG_R; 
            break;
        case PP_R_TOT:
        case PP_R_NUM:
            bestR = 0.0;    
            break;
        default:
            printf("Bad rval code = %d\n",ppPO->rval);
            ERR("ScoreSeqAntiParallelI","Bogus rval code");
            return(FALSE);
    }
    best = 0;
    /***
    *   Try each of m comparisons, remembering best score and offset of this
    */
    f = flen;
    s = slen;
    for(i=0;i<m;i++)
    {
        f--;
        s--;
        f = MAX_NUM(f,0);
        s = MAX_NUM(s,0);
        mind = MIN_NUM(ppPO->min_word+i, mover-i);
        n = MIN_NUM(minlen,mind);
        DB_PAIR DB_PrI("\n+ M[%d] n=%d f=%d s=%d mind=%d\n",i,n,f,s,mind);
        /***
        *   If clamps, make sure these are ok before regular comparision
        */
        if( (ppPO->cl3>0) || (ppPO->alg==ALGO_CON3) )
        {
            if( (f+n) < flen) 
            {
                DB_PAIR DB_PrI("+  PAST 3' end, can't clamp = BREAK\n");
                break;
            }
            if( (flen - ppPO->cl3 - f) < 0) 
            {
                DB_PAIR DB_PrI("+  Too short 3' end, can't clamp %d; next\n",
                    ppPO->cl3);
                continue;
            }
            if(!OkAlignClampI(ppPO,fS, flen-ppPO->cl3, sS,s,ppPO->cl3))
            {
                DB_PAIR DB_PrI("+  No 3' clamp = next\n");
                continue;
            }
        }
        if( (ppPO->cl5>0) || (ppPO->alg==ALGO_CON5) )
        {
            if(f > 0) 
            {
                DB_PAIR DB_PrI("+  past 5' end, can't clamp = next\n");
                continue;
            }
            if(n < ppPO->cl5) 
            {
                DB_PAIR DB_PrI("+  Too short 5' end, can't clamp = BREAK\n");
                continue;
            }
            if(!OkAlignClampI(ppPO,fS,0,sS,s,ppPO->cl5))
            {
                DB_PAIR DB_PrI("+  No 5' clamp = next\n");
                continue;
            }   
        }
        /***
        *   Compare region of n letters starting from fS[f] and sS[s]
        *   Always save the "best" offset, regardless of score returned
        */
        scR = AlignedSeqScoreR(ppPO,fS,f,sS,s,n);
        if(scR > bestR)
        {
            best = find-i;
        }
        /***
        *   What to do with score?
        */
        switch(ppPO->rval)
        {
            case PP_R_MAX:
                if(scR > bestR)
                {
                    bestR = scR;
                    DB_PAIR DB_PrI("+  Best score = %3.3f @%d\n",bestR,best);
                }   
                break;
            case PP_R_NUM:
                if(scR >= ppPO->thresh)
                {
                    bestR += 1.0;
                    DB_PAIR DB_PrI("+  Best count = %2.0f @%d\n",bestR,best);
                }
                break;
            case PP_R_TOT:
                if(scR >= ppPO->thresh)
                {
                    bestR += scR;
                    DB_PAIR DB_PrI("+  Best tot = %3.3f @%d\n",bestR,best);
                }
                break;
        }
    }
    *scPR = bestR;
    DB_PAIR DB_PrI("<< ScoreSeqAntiParallelI %d (%4.2f)\n",best,bestR);
    return(best);
}
/***************************************************************************
*   Check that specified end "clamp" is ok
*/
int OkAlignClampI(PPARS *ppPO, char *fS,int fi,char *sS,int si,int n)
{
    int i,ok,pind;

    DB_PAIRLO DB_PrI(">> OkAlignClampI n=%d\n",n);
    if(n<1)
    {
        return(TRUE);
    }
    DB_PAIRLO 
    {
        DB_PrI("+  F|"); PrintString(&fS[fi],n,NULL); DB_PrI("| fi=%d\n",fi); 
        DB_PrI("+  S|"); PrintString(&sS[si],n,NULL); DB_PrI("| si=%d\n",si); 
    }
    ok = TRUE;
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
        case PP_COM_PCOM:
            DB_PAIRLO DB_PrI("+ PP_COM_SIM / PP_COM_PCOM\n");
            for(i=0;i<n;i++)
            {
                pind = SeqBasePairIndexI(fS[fi+i],sS[si+i]);
                if(IS_BOG(pind))
                {
                    ok=FALSE;
                }
                else
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[si+i],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] <= 0.0)
                    {
                        ok = FALSE;
                    }
                }
                if(!ok)
                {
                    break;
                }
            }
            break;
        case PP_COM_COM:
            DB_PAIRLO DB_PrI("+ PP_COM_COM\n");
            for(i=0;i<n;i++)
            {
                pind = SeqBasePairIndexI(fS[fi+i],sS[si+n-i-1]);
                if(IS_BOG(pind))
                {
                    ok=FALSE;
                }
                else
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                            fS[fi+i],sS[si+n-i-1],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] <= 0.0)
                    {
                        ok = FALSE;
                    }
                }
                if(!ok)
                {
                    break;
                }
            }
            break;
        default:
            printf("Bogus compare code = %d\n",ppPO->ctype);
            ERR("OkAlignClampI","Bad code");
            return(FALSE);
    }
    DB_PAIRLO DB_PrI("<< OkAlignClampI %d\n",ok);
    return(ok);
}
/***************************************************************************
*   Compare n letters in aligned seqs fS and sS, starting at indices fi & si
*/
REAL AlignedSeqScoreR(PPARS *ppPO, char *fS,int fi,char *sS,int si,int n)
{
    REAL scR;

    DB_PAIRLO 
    {
        DB_PrI(">> AlignedSeqScoreR fi=%d si=%d n=%d alg=%d ctype=%d\n",
            fi, si, n, ppPO->alg, ppPO->ctype);
        DB_PrI("+  F|"); PrintString(&fS[fi],n,NULL); DB_PrI("|\n"); 
        DB_PrI("+  S|"); PrintString(&sS[si],n,NULL); DB_PrI("|\n"); 
    }
    switch(ppPO->alg)
    {
        case ALGO_MATCH: 
            scR = AlignedSeqMatchScoreR(ppPO,fS,fi,sS,si,n);
            break;
        case ALGO_CONT:
            scR = AlignedSeqContScoreR(ppPO,fS,fi,sS,si,n);
            break;
        case ALGO_CON5:
            scR = AlignedSeqCon5ScoreR(ppPO,fS,fi,sS,si,n);
            break;
        case ALGO_CON3:
            scR = AlignedSeqCon3ScoreR(ppPO,fS,fi,sS,si,n);
            break;
        case ALGO_BMATCH:
            scR = AlignedSeqBmatchScoreR(ppPO,fS,fi,sS,si,n);   
            break;
        default:
            printf("Bad pairing alg code = %d\n",ppPO->alg);
            ERR("AlignedSeqScoreR","Bogus algorithm code");
            scR = BAD_R;
    }
    DB_PAIRLO DB_PrI("<< AlignedSeqScoreR %3.3f\n",scR);
    return(scR);
}
/***************************************************************************
*   Total matching letters
*/
REAL AlignedSeqMatchScoreR(PPARS *ppPO, char *fS,int fi,char *sS,int si,int n)
{
    int i,j,pind;
    REAL scR;

    DB_PAIRLO DB_PrI(">> AlignedSeqMatchScoreR fi=%d si=%d n=%d\n",fi,si,n);
    scR = 0.0;
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
        case PP_COM_PCOM:
            DB_PAIRLO DB_PrI("+ PP_COM_SIM / PP_COM_PCOM\n");
            for(i=0;i<n;i++)
            {
                pind = SeqBasePairIndexI(fS[fi+i],sS[si+i]);
                if(!IS_BOG(pind))
                {
                    scR += ppPO->mscore[pind];
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f %3.2f\n",
                        fS[fi+i],sS[si+i],pind,ppPO->mscore[pind],scR);
                }
            }
            break;
        case PP_COM_COM:
            DB_PAIRLO DB_PrI("+ PP_COM_COM\n");
            j = si+n-1;
            for(i=0;i<n;i++)
            {
                pind = SeqBasePairIndexI(fS[fi+i],sS[j]);
                if(!IS_BOG(pind))
                {
                    scR += ppPO->mscore[pind];
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f %3.2f\n",
                        fS[fi+i],sS[j],pind,ppPO->mscore[pind],scR);
                }
                j--;
            }
            break;
        default:
            printf("Bad ctype code = %d\n",ppPO->ctype);
            ERR("AlignedSeqMatchScoreR","Bogus ctype code");
    }
    DB_PAIRLO DB_PrI("<< AlignedSeqMatchScoreR %4.2f\n",scR);
    return(scR);
}
/***************************************************************************
*   Max contiguous matching letters
*/
REAL AlignedSeqContScoreR(PPARS *ppPO, char *fS,int fi,char *sS,int si,int n)
{
    int i,j,c,pind;
    REAL scR,matR;

    DB_PAIRLO DB_PrI(">> AlignedSeqContScoreR fi=%d si=%d n=%d\n",fi,si,n);
    scR = matR = 0.0;
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
        case PP_COM_PCOM:
            for(i=0;i<n;i++)
            {
                c = 0;
                pind = SeqBasePairIndexI(fS[fi+i],sS[si+i]);
                if(!IS_BOG(pind))
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[si+i],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] > 0.0)
                    {
                        matR += ppPO->mscore[pind]; 
                        c++;
                    }
                }
                if(!c)
                {
                    scR = MAX_NUM(matR,scR);
                    matR = 0.0;
                }
            }
            scR = MAX_NUM(matR,scR);
            break;
        case PP_COM_COM:
            j = si+n-1;
            for(i=0;i<n;i++)
            {
                c = 0;
                pind = SeqBasePairIndexI(fS[fi+i],sS[j]);
                if(!IS_BOG(pind))
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[j],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] > 0.0)
                    {
                        matR += ppPO->mscore[pind]; 
                        c++;
                    }
                }
                if(!c)
                {
                    scR = MAX_NUM(matR,scR);
                    matR = 0.0;
                }
                j--;
            }
            scR = MAX_NUM(matR,scR);
            break;
        default:
            printf("Bad ctype code = %d\n",ppPO->ctype);
            ERR("AlignedSeqContScoreR","Bogus ctype code");
    }
    DB_PAIRLO DB_PrI("<< AlignedSeqContScoreR %4.2f\n",scR);
    return(scR);
}
/***************************************************************************
*   Max contiguous matching letters from 5' end
*/
REAL AlignedSeqCon5ScoreR(PPARS *ppPO, char *fS,int fi,char *sS,int si,int n)
{
    int i,j,pind,c;
    REAL scR;

    DB_PAIRLO DB_PrI(">> AlignedSeqCon5ScoreR fi=%d si=%d n=%d\n",fi,si,n);
    scR = 0;
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
        case PP_COM_PCOM:
            /***
            *   Front of fS to sS miss-match
            */
            for(i=0;i<n;i++)
            {
                c = 0;
                pind = SeqBasePairIndexI(fS[fi+i],sS[si+i]);
                if(!IS_BOG(pind))
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[si+i],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] > 0.0)
                    {
                        scR += ppPO->mscore[pind];
                        c++;
                    }
                }
                if(!c)
                {
                    break;
                }
            }
            break;
        case PP_COM_COM:
            /***
            *   Front of fS to end of sS until miss-match
            */
            j = si+n-1;
            for(i=0;i<n;i++)
            {
                c = 0;
                pind = SeqBasePairIndexI(fS[fi+i],sS[j]);
                if(!IS_BOG(pind))
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[j],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] > 0.0)
                    {
                        scR += ppPO->mscore[pind];
                        c++;
                    }
                }
                if(!c)
                {
                    break;
                }
                j--;
            }
            break;
        default:
            printf("Bad ctype code = %d\n",ppPO->ctype);
            ERR("AlignedSeqCon5ScoreR","Bogus ctype code");
    }
    DB_PAIRLO DB_PrI("<< AlignedSeqCon5ScoreR %2.0f\n",scR);
    return(scR);
}
/***************************************************************************
*   Max contiguous matching letters from 3' end of fS
*/
REAL AlignedSeqCon3ScoreR(PPARS *ppPO,char *fS,int fi,char *sS,int si,int n)
{
    int i,j,pind,c;
    REAL scR;

    DB_PAIRLO DB_PrI(">> AlignedSeqCon3ScoreR fi=%d si=%d n=%d\n",fi,si,n);
    scR = 0.0;
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
        case PP_COM_PCOM:
            /***
            *   Start at end of fS compared to sS
            */
            for(i=n-1;i>=0;i--)
            {
                c = 0;
                pind = SeqBasePairIndexI(fS[fi+i],sS[si+i]);
                if(!IS_BOG(pind))
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[si+i],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] > 0.0)
                    {
                        scR += ppPO->mscore[pind];
                        c++;
                    }
                }
                if(!c)
                {
                    break;
                }
            }
            break;
        case PP_COM_COM:
            /***
            *   Start at end of fS, compared to start of sS; break on miss-mat
            */
            j = si;
            for(i=n-1;i>=0;i--)
            {
                c = 0;
                pind = SeqBasePairIndexI(fS[fi+i],sS[j]);
                if(!IS_BOG(pind))
                {
                    DB_PAIRLO DB_PrI("+    %c--%c pind=%d %3.2f\n",
                        fS[fi+i],sS[j],pind,ppPO->mscore[pind]);
                    if(ppPO->mscore[pind] > 0.0)
                    {
                        scR += ppPO->mscore[pind];
                        c++;
                    }
                }
                if(!c)
                {
                    break;
                }
                j++;
            }
            break;
        default:
            printf("Bad ctype code = %d\n",ppPO->ctype);
            ERR("AlignedSeqCon3ScoreR","Bogus ctype code");
    }
    DB_PAIRLO DB_PrI("<< AlignedSeqCon3ScoreR %2.0f\n",scR);
    return(scR);
}
/**************************************************************************
*   Block-weighted match between letters
*/
REAL AlignedSeqBmatchScoreR(PPARS *ppPO,char *fS,int fi,char *sS,int si,int n)
{
    int i,j,mat,c;
    REAL scR;

    DB_PAIRLO DB_PrI(">> AlignedSeqBmatchScoreR fi=%d si=%d n=%d\n",fi,si,n);
    scR = 0.0;
    mat = 0;
    switch(ppPO->ctype)
    {
        case PP_COM_SIM:
            for(i=0;i<n;i++)
            {
                /***
                *   Self then one-back, one-up
                */
                if(UPPER(fS[fi+i]) != UPPER(sS[si+i]))
                {
                    continue;
                }
                mat += 2;
                if(i>0)
                {
                    if(UPPER(fS[fi+i-1]) == UPPER(sS[si+i-1]))
                        mat++;
                }
                if(i<(n-1))
                {
                    if(UPPER(fS[fi+i+1]) == UPPER(sS[si+i+1]))
                        mat++;
                }
            }
            break;
        case PP_COM_PCOM:
            for(i=0;i<n;i++)
            {
                c = 0;
                switch(UPPER(fS[fi+i]))
                {
                    case 'A': if(UPPER(sS[si+i])=='T') c++; break;
                    case 'C': if(UPPER(sS[si+i])=='G') c++; break;
                    case 'G': if(UPPER(sS[si+i])=='C') c++; break;
                    case 'T': if(UPPER(sS[si+i])=='A') c++; break;
                }
                if(!c)
                {
                    continue;
                }
                mat += 2;
                if(i>0)
                {
                    switch(UPPER(fS[fi+i-1]))
                    {
                        case 'A': if(UPPER(sS[si+i-1])=='T') mat++; break;
                        case 'C': if(UPPER(sS[si+i-1])=='G') mat++; break;
                        case 'G': if(UPPER(sS[si+i-1])=='C') mat++; break;
                        case 'T': if(UPPER(sS[si+i-1])=='A') mat++; break;
                    }
                }
                if(i<(n-1))
                {
                    switch(UPPER(fS[fi+i+1]))
                    {
                        case 'A': if(UPPER(sS[si+i+1])=='T') mat++; break;
                        case 'C': if(UPPER(sS[si+i+1])=='G') mat++; break;
                        case 'G': if(UPPER(sS[si+i+1])=='C') mat++; break;
                        case 'T': if(UPPER(sS[si+i+1])=='A') mat++; break;
                    }
                }
            }
            break;
        case PP_COM_COM:
            j = si+n-1;
            for(i=0;i<n;i++)
            {
                c = 0;
                switch(UPPER(fS[fi+i]))
                {
                    case 'A': if(UPPER(sS[j])=='T') c++;    break;
                    case 'C': if(UPPER(sS[j])=='G') c++;    break;
                    case 'G': if(UPPER(sS[j])=='C') c++;    break;
                    case 'T': if(UPPER(sS[j])=='A') c++;    break;
                }
                if(!c)
                {
                    j--;
                    continue;
                }
                mat += 2;
                if( (i>0) && (j<(n-1)) )
                {
                    switch(UPPER(fS[fi+i-1]))
                    {
                        case 'A': if(UPPER(sS[j+1])=='T') mat++;    break;
                        case 'C': if(UPPER(sS[j+1])=='G') mat++;    break;
                        case 'G': if(UPPER(sS[j+1])=='C') mat++;    break;
                        case 'T': if(UPPER(sS[j+1])=='A') mat++;    break;
                    }
                }
                if( (i<(n-1)) && (j>0) )
                {
                    switch(UPPER(fS[fi+i+1]))
                    {
                        case 'A': if(UPPER(sS[j-1])=='T') mat++;    break;
                        case 'C': if(UPPER(sS[j-1])=='G') mat++;    break;
                        case 'G': if(UPPER(sS[j-1])=='C') mat++;    break;
                        case 'T': if(UPPER(sS[j-1])=='A') mat++;    break;
                    }
                }
                j--;
            }
            break;
        default:
            printf("Bad ctype code = %d\n",ppPO->ctype);
            ERR("AlignedSeqWmatchScoreR","Bogus ctype code");
    }
    scR = RNUM(mat)/4.0;
    DB_PAIRLO DB_PrI("<< AlignedSeqBmatchScoreR %2.3f\n",scR);
    return(scR);
}
