/*
* seq_info.c
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
#include <math.h>
#include "prim.h"
#include "divinds.h"
#include "seq_info.h"

#define DB_INFO if(DB[88])

/**************************************************************************
*   Get score for number of words size min to max in seq relative to possible
*/
DOUB SeqWordFreqInfoD(char *seqPC,int len,int min,int max)
{
    int i,w,n,nd,p,m,mask[MAX_WFI];
    DOUB dD,dwD;
    char wordS[MAX_WFI+1],seqS[MAX_WFI+1],*sPC,*cPC;
    
    DB_INFO
    {
        DB_PrI(">> SeqWordFreqInfoD |");
        PrintString(seqS,len,stdout);
        DB_PrI("|\n");
    }
    if(len > MAX_WFI)
    {
        printf("SHAM; SeqWordFreqInfoD won't fly with seqs %d long\n",len);
        return(-1.0);
    }
    if(max < min)
    {
        printf("SHAM; SeqWordFreqInfoD; Bogus word sizes = %d to %d\n",min,max);
        return(-1.0);
    }
    /***
    *   Get seq into local, null-term string
    */
    strncpy(seqS,seqPC,len);
    seqS[len] = '\0';
    /***
    *   Scan for each word size
    */
    dD = 0.0;
    p = len - min + 1;
    for(w=min; w<=max; w++)
    {
        /***
        *   Clear masking
        */
        for(i=0;i<len;i++)
        {
            mask[i] = 0;
        }
        /***
        *   Scan along seq len
        */
        nd = 0;
        for(i=0; i<p; i++)
        {
            /***
            *   If this position is masked, ignore
            */
            if(mask[i])
            {
                continue;
            }
            nd++;
            /***
            *   Get current word; at least one occurrence 
            */
            sPC = &seqS[i];
            strncpy(wordS,sPC,w);
            wordS[w] = '\0';
            n = 1;
            DB_INFO DB_PrI("+ %s = %d ",wordS,i);
            /***
            *   Scan rest of seq for word
            */
            sPC++;
            while( (cPC=strstr(sPC,wordS)) != NULL)
            {
                m = cPC - seqS;
                mask[m] = TRUE;
                DB_INFO DB_PrI(" %d ",m);
                n++;
                sPC = cPC+1;
            }
            DB_INFO DB_PrI("= %d\n",n);
        }
        if(p > 1)
            dwD = (RNUM(nd-1)/RNUM(p)) / (1.0 - 1.0/RNUM(p));
        else
            dwD = 1.0;
        DB_INFO DB_PrI("+ %d diff size %d; p=%d dw=%f\n",nd,w,p,dwD);
        dD += dwD;
        p--;
    }
    dD /= RNUM(max-min+1);
    DB_INFO DB_PrI("<< SeqWordFreqInfoD %f\n",dD);
    return(dD);
}
/*********************************************************************
*   SHAM; should use wordfreq general primitives and not this limited hack
*/
int Fill123WordCountsI(char *seqS,int len,int *onePI,int *twoPI,int *thrPI)
{
    int i,one,two,thr,n;

    /***
    *   Initalize to zero
    */
    for(i=0;i<4;i++)
    {
        onePI[i] = 0;
    }
    for(i=0;i<16;i++)
    {
        twoPI[i] = 0;
    }
    for(i=0;i<64;i++)
    {
        thrPI[i] = 0;
    }
    /***
    *   Count the whip
    */
    n = 0;
    for(i=0;i<len;i++)
    {
        one = two = thr = BOGUS;
        switch(toupper(seqS[i]))
        {
            case 'A':   one=0; break;
            case 'C':   one=1; break;
            case 'G':   one=2; break;
            case 'T':   one=3; break;
        }
        if(i<(len-1))
        {
            switch(toupper(seqS[i+1]))
            {
                case 'A':   two=0; break;
                case 'C':   two=1; break;
                case 'G':   two=2; break;
                case 'T':   two=3; break;
            }
        }
        if(i<(len-2))
        {
            switch(toupper(seqS[i+2]))
            {
                case 'A':   thr=0; break;
                case 'C':   thr=1; break;
                case 'G':   thr=2; break;
                case 'T':   thr=3; break;
            }
        }
        if( !IS_BOG(one) )
        {
            onePI[one] += 1;
            n++;
        }
        if( (!IS_BOG(one)) && (!IS_BOG(two)) )
        {
            twoPI[(one*4)+two] += 1;
        }
        if( (!IS_BOG(one)) && (!IS_BOG(two)) && (!IS_BOG(thr)) )
        {
            thrPI[(one*16)+(two*4)+thr] += 1;
        }
    }
    return(n);
}
/*********************************************************************
*   Hack shannon index; 1-mer, 2-mer, 3-mer sham
*/
int Seq123ShannonInfoI(char *seqPC, int len, DOUB *onePD, DOUB *twoPD,
    DOUB *thrPD)
{
    int onePI[4],twoPI[16],thrPI[64];

    if(Fill123WordCountsI(seqPC,len,onePI,twoPI,thrPI)!=len)
    {
        return(FALSE);
    }
    *onePD = ShannonWeaverIndexD(len,onePI,4);
    *twoPD = ShannonWeaverIndexD(len-1,twoPI,16);
    *thrPD = ShannonWeaverIndexD(len-2,thrPI,64);
    return(TRUE);   
}
/*********************************************************************
*   Hack shannon even index; 1-mer, 2-mer, 3-mer sham
*/
int Seq123ShannonEvenInfoI(char *seqPC, int len, DOUB *onePD, DOUB *twoPD,
    DOUB *thrPD)
{
    int onePI[4],twoPI[16],thrPI[64];

    if(Fill123WordCountsI(seqPC,len,onePI,twoPI,thrPI)!=len)
    {
        return(FALSE);
    }
    *onePD = ShannonWeEvenIndexD(len,onePI,4);
    *twoPD = ShannonWeEvenIndexD(len-1,twoPI,16);
    *thrPD = ShannonWeEvenIndexD(len-2,thrPI,64);
    return(TRUE);   
}
/*********************************************************************
*   Hack brillouin index; 1-mer, 2-mer, 3-mer sham
*/
int Seq123BrillouinInfoI(char *seqPC, int len, DOUB *onePD, DOUB *twoPD,
    DOUB *thrPD)
{
    int onePI[4],twoPI[16],thrPI[64];

    if(Fill123WordCountsI(seqPC,len,onePI,twoPI,thrPI)!=len)
    {
        return(FALSE);
    }
    *onePD = BrillouinIndexD(len,onePI,4);
    *twoPD = BrillouinIndexD(len-1,twoPI,16);
    *thrPD = BrillouinIndexD(len-2,thrPI,64);
    return(TRUE);   
}
/*********************************************************************
*   Hack berger-parker-like index; 1-mer, 2-mer, 3-mer sham
*   I thought up the normalization scheme so 0-1
*/
int Seq123BergerParkerInfoI(char *seqPC, int len, DOUB *onePD, DOUB *twoPD,
    DOUB *thrPD)
{
    int i,max1,max2,max3,onePI[4],twoPI[16],thrPI[64];

    if(Fill123WordCountsI(seqPC,len,onePI,twoPI,thrPI)!=len)
    {
        return(FALSE);
    }
    /***
    *   Find max for each size
    */
    max1 = max2 = max3 = 0;
    for(i=0;i<4;i++)
    {
        max1 = MAX_NUM(onePI[i],max1);
    }
    for(i=0;i<16;i++)
    {
        max2 = MAX_NUM(twoPI[i],max2);
    }
    for(i=0;i<64;i++)
    {
        max3 = MAX_NUM(thrPI[i],max3);
    }
    *onePD = NormBergerParkerIndexD(max1,len,4);
    *twoPD = NormBergerParkerIndexD(max2,len-1,16);
    *thrPD = NormBergerParkerIndexD(max3,len-1,64);
    return(TRUE);   
}
/**************************************************************************
*   Normalized (my sham) Berger-parker Diversity Index
*/
DOUB NormBergerParkerIndexD(int max,int num,int adim)
{
    DOUB dD;

    if( (max<1) || (num<1) || (adim<2) )
    {
        return(BAD_R);
    }
    dD = (1.0 - DNUM(max)/DNUM(num)) * (DNUM(adim)/DNUM(adim-1));
    return(dD); 
}
