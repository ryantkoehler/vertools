/*
* dna_loop.c
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
#include <math.h>
#include "prim.h"
#include "dna.h"
#include "dna_pair.h"


#define DB_PAIR     if(DB[116])
#define DB_PAIRLO   if(DB[117])
    
/****************************************************************************
*   Score all possible "loop" overlaps of sequence sS[len] and return the
*       "offset" corresponding to the highest scoring; 
*       The actual score is returned via *scPR.
*   The pairing-parameter data structure holds the value of the smallest
*       allowable loop and the minimum amount of overlap to consider.
*   Returned offset is the position of sequence (3') end relative to the (5')
*       begining; if these two align together, offset = 0.
*       For the sequence TAACGGCCGTTAGG, offset = -2 (with min_loop 4)
*       5'>   TAAC-GG
*             ||||  :
*       3'< GGATTG-CC
*/
int ScoreLoopOverlapI(char *sS,int len,PPARS *ppPO,REAL *scPR)
{
    int i,f,s,m,n,find,best;
    REAL scR,bestR;

    DB_PAIR DB_PrI(">> ScoreLoopOverlapI len=%d (loop=%d over=%d)\n",len,
        ppPO->min_loop,ppPO->min_word);
    /***
    *   m = number of alignments to compare
    */
    m = 2 * (len - ppPO->min_loop - (2 * ppPO->min_word)) + 1;
    find = len - ppPO->min_loop - (2 * ppPO->min_word);
    bestR = -TOO_BIG_R;
    best = 0;
    DB_PAIR DB_PrI("+ m=%d find=%d\n",m,find);
    for(i=0;i<m;i++)
    {
        f = MAX_NUM(0,find-i);
        s = MAX_NUM(0,i-find);
        n = (len - MAX_NUM(f,s) - ppPO->min_loop)/2;
        DB_PAIR DB_PrI("+ M[%d] n=%d [%d][%d]\n",i,n,f,len-s-1);
        scR = SeqLoopScoreR(sS, len, f, len-s-1, n, ppPO);
        if(scR > bestR)
        {
            bestR = scR;
            best = find-i;
            DB_PAIR DB_PrI("+  Best=%d score=%4.2f\n",best,bestR);
        }   
    }
    *scPR = bestR;
    DB_PAIR DB_PrI("<< ScoreLoopOverlapI %d (%4.2f)\n",best,bestR);
    return(best);
}
/***************************************************************************
*   Compare seq subset to itself to evaluate loop overlap
*   sS is sequence [len]
*   fi & si are first & second starting indices, for compare of n base stretch
*/
REAL SeqLoopScoreR(char *sS,int len,int fi,int si,int n,PPARS *ppPO)
{
    REAL scR;
    int i,mat,con,m;

    DB_PAIRLO DB_PrI(">> SeqLoopScoreR len=%d fi=%d si=%d n=%d\n",len,fi,si,n);
    scR = 0.0;
    mat = con = 0;
    switch(ppPO->alg)
    {   
        case ALGO_MATCH: 
            for(i=0;i<n;i++)
            {
                switch(sS[fi+i])
                {
                    case 'A': case 'a': 
                        if( (sS[si-i]=='T') || (sS[si-i]=='t') )
                            mat++;
                        break;
                    case 'C': case 'c': 
                        if( (sS[si-i]=='G') || (sS[si-i]=='g') )
                            mat++;
                        break;
                    case 'G': case 'g': 
                        if( (sS[si-i]=='C') || (sS[si-i]=='c') )
                            mat++;
                        break;
                    case 'T': case 't': 
                        if( (sS[si-i]=='A') || (sS[si-i]=='a') )
                            mat++;
                        break;
                    default:
                        printf("Bogus base [%d]|%c|\n",fi+i,sS[fi+i]);
                        ERR("SeqLoopScoreR","Bad seq"); return(-TOO_BIG_R);
                }
            }
            scR = RNUM(mat);
            break;
        case ALGO_CONT: 
            for(i=0;i<n;i++)
            {
                m = 0;
/**
                DB_PAIRLO DB_PrI("+ comp[%d][%d] ",fi+i,si-i);
                DB_PAIRLO DB_PrI(" %c %c",sS[fi+i],sS[si-i]);
**/
                switch(sS[fi+i])
                {
                    case 'A': case 'a': 
                        if( (sS[si-i]=='T') || (sS[si-i]=='t') )
                            m++;    
                        break;
                    case 'C': case 'c': 
                        if( (sS[si-i]=='G') || (sS[si-i]=='g') )
                            m++;    
                        break;
                    case 'G': case 'g': 
                        if( (sS[si-i]=='C') || (sS[si-i]=='c') )
                            m++;    
                        break;
                    case 'T': case 't': 
                        if( (sS[si-i]=='A') || (sS[si-i]=='a') )
                            m++;    
                        break;
                    default:
                        printf("Bogus base [%d]|%c|\n",fi+i,sS[fi+i]);
                        ERR("SeqLoopScoreR","Bad seq"); return(-TOO_BIG_R);
                }
                if(m)
                {
                    con++;
                }   
                else
                {
                    con = 0;
                }
                mat = MAX_NUM(mat,con);
/**
                DB_PAIRLO DB_PrI(" m=%d con=%d mat=%d\n",m,con,mat);
**/
            }
            scR = RNUM(mat);
            break;
        case ALGO_BMATCH: 
            for(i=0;i<n;i++)
            {
                /***
                *   Current position
                */
                con = FALSE;
                switch(toupper(INT(sS[fi+i])))
                {
                    case 'A': if(sS[si-i]=='T'){con++;} break;
                    case 'C': if(sS[si-i]=='G'){con++;} break;
                    case 'G': if(sS[si-i]=='C'){con++;} break;
                    case 'T': if(sS[si-i]=='A'){con++;} break;
                    default:
                        printf("Bogus base [%d]|%c|\n",fi+i,sS[fi+i]);
                        ERR("SeqLoopScoreR","Bad seq"); return(-TOO_BIG_R);
                }
                if(con == FALSE)
                {
                    continue;
                }
                mat += 2;
                /***
                *   Previous position; away from hairpin loop itself.
                *   Check if within seq dims
                */
                if( ((fi+i-1)>=0) && ((si-i+1)<len) )
                {
                    con = FALSE;
                    switch(toupper(INT(sS[fi+i-1])))
                    {
                        case 'A': if(sS[si-i+1]=='T'){con++;}   break;
                        case 'C': if(sS[si-i+1]=='G'){con++;}   break;
                        case 'G': if(sS[si-i+1]=='C'){con++;}   break;
                        case 'T': if(sS[si-i+1]=='A'){con++;}   break;
                        default:
                            printf("Bogus base [%d]|%c|\n",fi+i-1,sS[fi+i-1]);
                            ERR("SeqLoopScoreR","Bad seq"); return(-TOO_BIG_R);
                    
                    }
                    if(con)
                        mat++;
                }
                /***
                *   Next position is closer to the hairpin loop; ignore if in
                *       the loop via check against number to compare, n
                */
                if((i+1)<n)
                {
                    con = FALSE;
                    switch(toupper(INT(sS[fi+i+1])))
                    {
                        case 'A': if(sS[si-i-1]=='T'){con++;}   break;
                        case 'C': if(sS[si-i-1]=='G'){con++;}   break;
                        case 'G': if(sS[si-i-1]=='C'){con++;}   break;
                        case 'T': if(sS[si-i-1]=='A'){con++;}   break;
                        default:
                            printf("Bogus base [%d]|%c|\n",fi+i+1,sS[fi+i+1]);
                            ERR("SeqLoopScoreR","Bad seq"); return(-TOO_BIG_R);
                    
                    }
                    if(con)
                        mat++;
                }
            }
            scR = RNUM(mat)/4.0;
            break;
        default:
            scR = -TOO_BIG_R;
    }
    DB_PAIRLO DB_PrI("<< SeqLoopScoreR %4.2f\n",scR);
    return(scR);
}
