/*
* snp_char.c
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
#include "prim.h"
#include "dna.h"


/**********************************************************************
*   Convert explicit SNP string into one-base code
*/
int CollapseSnpStringI(char *snpS,char *scS)
{
    int n,a,c,g,t,d,s,x;
    char *cPC;

    /***
    *   Count (single!) bases in string
    */
    a = c = g = t = d = s = x = 0;
    cPC = snpS;
    n = 0;
    while(isgraph(INT(*cPC)))
    {
        switch(*cPC)
        {
            case 'A': case 'a': a++;    break;
            case 'C': case 'c': c++;    break;
            case 'G': case 'g': g++;    break;
            case 'T': case 't': t++;    break;
            case '*': case '-': d++;    break;
            case '/':           s++;    break;
            case '[':           break;
            case ']':           break;
            default:    x++;
        }
        cPC++;
    }
    n = a + c + g + t;
    /***
printf("X |%s| s=%d n=%d x=%d d=%d\n",snpS,s,n,x,d);
    *   As-is case =    'X'     where X is acgtx exactly
    */
    if( ((n+x)==1) && (d<1) && (s<1) ) {
        scS[0] = snpS[0];
        scS[1] = '\0';
        return(strlen(scS));
    }
    /***
    *   indel case =    'X/Y'   where X or Y is deletion exactly
    */
    if( (d==1) && (n>0) && (s>0) && (x<1) ) {
        sprintf(scS,"-");
        return(strlen(scS));
    }
    /***
    *   Odd case =      '?'     Anything else
    */
    if( (s!=1) || (n!=2) || (x>0) || (d>0) ) {
        sprintf(scS,"x");
        return(strlen(scS));
    }
    /***
    *   "Normal" case = 'X/Y'   where X and Y are acgt exactly
    *   Case by case base combinations to get code
    */
    sprintf(scS,"z");
    if(a) {
        if(c) {
            if(g) {
                if(t) {
                    /*  A + C + G + T = N */
                    sprintf(scS,"N");
                }
                else {
                    /*  A + C + G = !T = V */
                    sprintf(scS,"V");
                }
            }
            else if(t) {
                /*  A + C + T = !G = H */
                sprintf(scS,"H");
            }
            else {
                /*  A + C = M */
                sprintf(scS,"M");
            }
        }
        else if(g) {
            if(t) {
                /*  A + G + T = !C = D */
                sprintf(scS,"D");
            }
            else {
                /*  A + G = R */
                sprintf(scS,"R");
            }
        }
        else if(t) {
            /*  A + T = W */
            sprintf(scS,"W");
        }
        else {
            /*  Only A */
            sprintf(scS,"A");
        }
    }
    else if(c) {
        if(g) {
            if(t) {
                /*  C + G + T = !A = B */
                sprintf(scS,"B");
            }
            else {
                /*  C + G = S */
                sprintf(scS,"S");
            }
        }
        else if(t) {
            /*  C + T = Y */
            sprintf(scS,"Y");
        }
        else {
            /*  Only C */
            sprintf(scS,"C");
        }
    }
    else if(g) {
        if(t) {
            /*  G + T = K */
            sprintf(scS,"K");
        }
        else {
            /*  Only G */
            sprintf(scS,"G");
        }
    }
    else if(t) {
        /*  Only T */
        sprintf(scS,"T");
    }
    return(strlen(scS));
}
/*************************************************************************
*   Expand various degenerate base codes into explicit AGCT strings
*/
int ExpandDegBaseI(char bC,char *snpS)
{
    int n;
    char seqS[DEF_BS];

    n = 0;
    INIT_S(seqS);
    switch(bC)
    {
        case 'A':   case 'a':   sprintf(seqS,"A");      n=1;    break;
        case 'B':   case 'b':   sprintf(seqS,"C/G/T");  n=3;    break;
        case 'C':   case 'c':   sprintf(seqS,"C");      n=1;    break;
        case 'D':   case 'd':   sprintf(seqS,"A/G/T");  n=3;    break;
        case 'G':   case 'g':   sprintf(seqS,"G");      n=1;    break;
        case 'H':   case 'h':   sprintf(seqS,"A/C/T");  n=3;    break;
        case 'K':   case 'k':   sprintf(seqS,"G/T");    n=2;    break;
        case 'M':   case 'm':   sprintf(seqS,"A/C");    n=2;    break;
        case 'N':   case 'n':   sprintf(seqS,"A/C/G/T");n=4;    break;
        case 'R':   case 'r':   sprintf(seqS,"A/G");    n=2;    break;
        case 'S':   case 's':   sprintf(seqS,"C/G");    n=2;    break;
        case 'T':   case 't':   sprintf(seqS,"T");      n=1;    break;
        case 'V':   case 'v':   sprintf(seqS,"A/C/G");  n=3;    break;
        case 'W':   case 'w':   sprintf(seqS,"A/T");    n=2;    break;
        case 'Y':   case 'y':   sprintf(seqS,"C/T");    n=2;    break;
    }
    if(snpS)
    {
        strcpy(snpS,seqS);
    }
    return(n);
}
/*************************************************************************
*   Returns the number of different possible bases for a given code
*   Failure = 0
*/
int BaseDegeneracyI(char bC)
{
    int n;

    n = 0;
    switch(bC)
    {
        case 'A':   case 'a':   n=1;    break;
        case 'B':   case 'b':   n=3;    break;
        case 'C':   case 'c':   n=1;    break;
        case 'D':   case 'd':   n=3;    break;
        case 'G':   case 'g':   n=1;    break;
        case 'H':   case 'h':   n=3;    break;
        case 'K':   case 'k':   n=2;    break;
        case 'M':   case 'm':   n=2;    break;
        case 'N':   case 'n':   n=4;    break;
        case 'R':   case 'r':   n=2;    break;
        case 'S':   case 's':   n=2;    break;
        case 'T':   case 't':   n=1;    break;
        case 'U':   case 'u':   n=1;    break;
        case 'V':   case 'v':   n=3;    break;
        case 'W':   case 'w':   n=2;    break;
        case 'Y':   case 'y':   n=2;    break;
    }
    return(n);
}
/*************************************************************************
*
*/
int IsSnpSeqExpandedI(char *snpS)
{
    if(strlen(snpS)>1)
    {
        return(TRUE);
    }
    return(FALSE);
}
/**********************************************************************
*   Is non-degenerate IUPAC SNP base bC a member of degenerate IUPAC base dC
*/
int IsDegenBaseI(char bC, char dC)
{
    int is;

    is = FALSE;
    switch(bC)
    {
        case 'a': case 'A':
            switch(dC)
            {
                case 'a': case 'A': is++; break;
                case 'm': case 'M': is++; break;
                case 'r': case 'R': is++; break;
                case 'w': case 'W': is++; break;
                case 'b': case 'B': is++; break;
                case 'n': case 'N': is++; break;
            }
            break;
        case 'c': case 'C':
            switch(dC)
            {
                case 'c': case 'C': is++; break;
                case 'm': case 'M': is++; break;
                case 's': case 'S': is++; break;
                case 'y': case 'Y': is++; break;
                case 'd': case 'D': is++; break;
                case 'n': case 'N': is++; break;
            }
            break;
        case 'g': case 'G':
            switch(dC)
            {
                case 'g': case 'G': is++; break;
                case 'r': case 'R': is++; break;
                case 's': case 'S': is++; break;
                case 'k': case 'K': is++; break;
                case 'h': case 'H': is++; break;
                case 'n': case 'N': is++; break;
            }
            break;
        case 't': case 'T':
            switch(dC)
            {
                case 't': case 'T': is++; break;
                case 'k': case 'K': is++; break;
                case 'w': case 'W': is++; break;
                case 'y': case 'Y': is++; break;
                case 'v': case 'V': is++; break;
                case 'n': case 'N': is++; break;
            }
            break;
    }
    return(is);
}
/*************************************************************************
*   Returns how many ambigous bases between s and e in sS
*/
int CountSeqAmbigsI(char *sS,int s,int e)
{
    int i,n;

    if(s < 0) {
        s = 0;
    }
    if(e < 0) {
        e = strlen(sS);
    }
    n = 0;
    for(i=s;i<e;i++)
    {
        /***
        *   Ignore SNP site / indel stuff
        */
        if( (sS[i] == '[') || (sS[i] == '/') || (sS[i] == ']') || 
            (sS[i] == '*') || (sS[i] == '-') )
        {
            continue;
        }
        switch(BaseDegeneracyI(sS[i]))
        {
            /***
            *   A,G,C,T and U
            */
            case 1: 
                break;
            /***
            *   Any other letter
            */
            default:
                n++;    break;
        }
    }
    return(n);
}
/*************************************************************************
*   Returns how many ambigs of specified degeneracy between s and e in sS
*   Degeneracy = 1 for normal ACGT, 2 for SWRYMK, 3 for BDHV, 4 for N
*/
int CountSeqAmbigDegensI(char *sS,int s,int e,int type)
{
    int i,n;

    if(s < 0)
    {
        s = 0;
    }
    if(e < 0)
    {
        e = strlen(sS);
    }
    n = 0;
    for(i=s;i<e;i++)
    {
        /***
        *   Ignore SNP site stuff
        */
        if( (sS[i] == '[') || (sS[i] == '/') || (sS[i] == ']') || 
            (sS[i] == '*') )
        {
            continue;
        }
        if(BaseDegeneracyI(sS[i])==type)
        {
            n++;
        }
    }
    return(n);
}
/*************************************************************************
*   Returns how many SNP site annotations are between s and e in sS
*/
int CountSeqSnpSitesI(char *sS,int s,int e)
{
    int i,n,snp;

    if(s < 0)
    {
        s = 0;
    }
    if(e < 0)
    {
        e = strlen(sS);
    }
    n = 0;
    snp = FALSE;
    for(i=s;i<e;i++)
    {
        switch(sS[i])
        {
            case '[':       /* Start of SNP record, so count */
                snp = TRUE;
                n++;
                break;
            case ']':       /* End of SNP record, only count if not in */
                if(!snp)
                {
                    n++;
                }
                snp = FALSE;
                break;
            case '/':       /* SNP allele seperator; only count if not in */
                if(!snp)
                {
                    n++;
                }
                break;
            case '*':       /* SNP deletion indicator; only count if not in */
                if(!snp)
                {
                    n++;
                }
                break;
        }
    }
    return(n);
}
