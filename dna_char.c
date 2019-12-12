/*
* dna_char.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
*
* The programs and source code of the vertools collection are free software.
* They are distributed in the hope that they will be useful,
* WITHOUT ANY WARRANTY OF FITNESS FOR ANY PARTICULAR PURPOSE.  
* 
* Permission is granted for research, educational, and possibly commercial use 
*   and modification as long as 1) Code and any derived works are not 
*   redistributed for any fee, and 2) Proper credit is given to the authors. 
*   If you wish to include this software in a product, or use it commercially,
*   please contact the authors.
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


/***************************************************************************
*   TRUE if in standard DNA alphabet
*   BOGUS if degenerate DNA code
*   FALSE otherwise
*/
int GoodDNABaseI(char c)
{
    switch(c)
    {
        /***
        *   Standard ones = TRUE;
        */
        case 'A':   case 'a':   
        case 'C':   case 'c':       
        case 'G':   case 'g':       
        case 'T':   case 't':       
            return(TRUE);
        /***
        *   Non-standard degenerate ones = BOGUS;
        */
        case 'N':   case 'n':   
        case 'B':   case 'b':   
        case 'D':   case 'd':   
        case 'H':   case 'h':   
        case 'V':   case 'v':   
        case 'S':   case 's':   
        case 'W':   case 'w':   
        case 'Y':   case 'y':   
        case 'K':   case 'k':   
        case 'M':   case 'm':   
        case 'R':   case 'r':   
            return(BOGUS);
    }
    return(FALSE);
}
/***************************************************************************
*   TRUE if c and s are complimentary bases in DNA alphabet
*/
int GoodDNACompBasesI(char c,char s)
{
    int g;

    g = FALSE;
    switch(c)
    {
        case 'A':   case 'a':   
            switch(s)
            {
                case 'T': case 't': g=TRUE;
            }
            break;
        case 'C':   case 'c':   
            switch(s)
            {
                case 'G': case 'g': g=TRUE;
            }
            break;
        case 'G':   case 'g':   
            switch(s)
            {
                case 'C': case 'c': g=TRUE;
            }
            break;
        case 'T':   case 't':   
            switch(s)
            {
                case 'A': case 'a': g=TRUE;
            }
            break;
    }
    return(g);
}
/***************************************************************************
*   Returns integer index for DNA alphabet character; BOGUS if not in libarary
*/
int DNABaseIndexI(char seqC)
{
    switch(seqC)
    {
        case 'A':   case 'a':   return(DNA_A);
        case 'C':   case 'c':   return(DNA_C);
        case 'G':   case 'g':   return(DNA_G);
        case 'T':   case 't':   return(DNA_T);
        case 'N':   case 'n':   return(DNA_N);
    }
    return(BOGUS);
}
/***************************************************************************
*   Returns char in DNA alphabet for integer index; NULL if out of bounds
*/
char DNAIndexBaseC(int ind)
{
    switch(ind)
    {
        case DNA_A: return('A');
        case DNA_C: return('C');
        case DNA_G: return('G');
        case DNA_T: return('T');
        case DNA_N: return('N');
    }
    return(0);
}
/***************************************************************************
*   Returns the sequence pair index; i.e. alphabetic order
*   If any unrecognized sequence, retruns BOGUS
*/
int SeqPairIndexI(char *seqS)
{
    return(SeqBasePairIndexI(seqS[0],seqS[1]));
}
/***************************************************************************
*   Returns the sequence pair index; i.e. alphabetic order
*   If any unrecognized sequence, retruns BOGUS
*/
int SeqBasePairIndexI(char c1, char c2)
{
    switch(c1)
    {
        case 'A': case 'a':
            switch(c2)
            {
                case 'A':   case 'a':   return(0);
                case 'C':   case 'c':   return(1);
                case 'G':   case 'g':   return(2);
                case 'T':   case 't':   return(3);
            }
            break;
        case 'C': case 'c':
            switch(c2)
            {
                case 'A':   case 'a':   return(4);
                case 'C':   case 'c':   return(5);
                case 'G':   case 'g':   return(6);
                case 'T':   case 't':   return(7);
            }
            break;
        case 'G': case 'g':
            switch(c2)
            {
                case 'A':   case 'a':   return(8);
                case 'C':   case 'c':   return(9);
                case 'G':   case 'g':   return(10);
                case 'T':   case 't':   return(11);
            }
            break;
        case 'T': case 't':
            switch(c2)
            {
                case 'A':   case 'a':   return(12);
                case 'C':   case 'c':   return(13);
                case 'G':   case 'g':   return(14);
                case 'T':   case 't':   return(15);
            }
            break;
    }
    return(BOGUS);
}
/***************************************************************************
*   Returns compliment for DNA alphabet; 0 if not in alphabet
*/
char CompDNABaseC(char c)
{
    switch(c)
    {
        case 'A':   return('T');
        case 'a':   return('t');
        case 'C':   return('G');
        case 'c':   return('g');
        case 'G':   return('C');
        case 'g':   return('c'); 
        case 'T':   return('A');
        case 't':   return('a');
        /***
        *   Degen and brackets
        */
        case 'N':   return('N');
        case 'n':   return('n');
        case '[':   return(']');
        case ']':   return('[');
        case 'S':   return('S');    /*  = [C/G] => [G/C] = S */
        case 's':   return('s');
        case 'W':   return('W');    /*  = [A/T] => [T/A] = W */
        case 'w':   return('w');
        case 'M':   return('K');    /*  = [A/C] => [T/G] = K */
        case 'm':   return('k');
        case 'K':   return('M');    /*  = [G/T] => [C/A] = M */
        case 'k':   return('m');
        case 'R':   return('Y');    /*  = [A/G] => [T/C] = Y */
        case 'r':   return('y');
        case 'Y':   return('R');    /*  = [C/T] => [G/A] = R */
        case 'y':   return('r');
        /***
        *   Triple degenerate 
        */
        case 'B':   return('V');    /*  = [C/G/T] => [G/C/A] = V */
        case 'b':   return('v');
        case 'D':   return('H');    /*  = [A/G/T] => [T/C/A] = H */
        case 'd':   return('h');    
        case 'H':   return('D');    /*  = [T/C/A] => [A/G/T] = D */
        case 'h':   return('d');    
        case 'V':   return('B');    /*  = [A/C/G] => [T/G/C] = B */
        case 'v':   return('b');    
        default:    return(c);
    }
}
/***************************************************************************
*   Returns complimentary DNA sequence 5'>-->3'
*   compS may be the same string as seqS
*/
int CompDNASeqI(char *seqS,int len,char *compS)
{
    int i;
    char uC,dC;

    for(i=0;i<(len/2);i++)
    {
        uC = seqS[i];
        dC = seqS[len-i-1];
        compS[i] = CompDNABaseC(dC);
        compS[len-i-1] = CompDNABaseC(uC);
    }
    /***
    *   If odd length, compliment middle base
    */
    if(len%2)
    {
        compS[len/2] = CompDNABaseC(seqS[len/2]);
    }
    return(i);
}
/***************************************************************************
*   Returns inverse DNA sequence of the one passed (both assumed 5'>-->3')
*/
int InverseDNASeqI(char *seqS,int len,char *compS)
{
    int i;

    for(i=0;i<len;i++)
    {
        switch(seqS[i])
        {
            case '[': 
            case ']': 
                compS[i] = seqS[i]; break;
            default: 
                compS[i] = CompDNABaseC(seqS[i]);
        }
    }
    return(i);
}
/***************************************************************************
*   Returns the reverse DNA seqeucne of the one passed (i.e. 3'>-->5')
*   Lenght of seqS is len, compS should be at least as big
*/
int ReverseDNASeqI(char *seqS,int len,char *compS)
{
    int i,j;
    char fC,sC;

    j = len-1;
    for(i=0;i<((len+1)/2);i++)
    {
        fC = seqS[i];
        sC = seqS[j];
        switch(sC)
        {
            case '[': compS[i] = ']';   break;
            case ']': compS[i] = '[';   break;
            default:
                compS[i] = sC;
        }
        switch(fC)
        {
            case '[': compS[j--] = ']'; break;
            case ']': compS[j--] = '['; break;
            default:
                compS[j--] = fC;
        }
    }
    return(len);
}
/***************************************************************************
*   Fills seqS[len] with a random sequence
*   If fracsPI is non-null, there should be 4 numbers corresponding to
*       percent (0-100) for ACGT (that should sum to 100)
*/
int RandomDNASeqI(char *seqS, int len, int *fracsPI)
{
    int i,j,r;
    static int set = FALSE;
    static char baseperS[100];

    /***
    *   Initialize the fraction array
    */
    if(set != TRUE)
    {
        if(fracsPI)
        {
            r = 0;
            for(i=0;i<4;i++)
            {   r += fracsPI[i];    }
            if( r != 100 ) 
            {
                printf("SHAM SNAPPED! r=%d\n",r);
                ERR("RandomDNASeqI","?");
            }
            r = 0;
            for(i=0;i<4;i++)
            {
                for(j=0;j<fracsPI[i];j++)
                { baseperS[r++] = DNAIndexBaseC(i); }
            }
        }
        else
        {
            for(i=0;i<100;i++)
            {
                switch(i/25)
                {
                    case 0: baseperS[i] = 'A';  break;
                    case 1: baseperS[i] = 'C';  break;
                    case 2: baseperS[i] = 'G';  break;
                    case 3: baseperS[i] = 'T';  break;
                    default:
                        printf("SHAM SNAPPED! i=%d i/25=%d\n",i,i/25);
                        ERR("RandomDNASeqI","?");
                }
            }
        }
        set = TRUE;
    }
/*
    for(i=0;i<100;i++)
    {
        printf("X %2d = %c\n",i,baseperS[i]);
    }
*/
    /**
    *   Now fill the sequence
    */
    for(i=0;i<len;i++)
    {
        r = RandI(100);
        seqS[i] = baseperS[r];
    }
    return(i);
}
/***************************************************************************
*   Randomly tweaks sequence of len with tweak changes
*/
int RandTweakDNASeqI(char *seqS,int len,int tweak)
{
    int i,r,b;
    REAL fR,rR;

    LIMIT_NUM(tweak,1,len);
    for(i=0;i<len;i++)
    {
        fR = RNUM(tweak)/RNUM(len-i); 
        rR = RandR(1.0);
        if(fR > rR)
        {
            b = DNABaseIndexI(seqS[i]);
            r = RandI(4);
            while(r==b)
            {
                r = RandI(4);
            }
            seqS[i] = DNAIndexBaseC(r);
            tweak--;
        }
        if(tweak==0)
            break;
    }
    return(i);
}
/****************************************************************************
*   Count the number of non-standard bases in the passed string
*/
int CountNonStandBasesI(char *seqS,int len)
{
    int i,n;
    
    n = 0;
    for(i=0;i<len;i++)
    {   
        switch(seqS[i])
        {
            case 'a':
            case 'A':
            case 'c':
            case 'C':
            case 'g':
            case 'G':
            case 't':
            case 'T':
                break;
            default:
                if(isgraph(INT(seqS[i])))
                    n++;
        }
    }
    return(n);
}
/****************************************************************************
*   Extract a subsequence from seq and put it into subS
*/
int ExtractForOrRevSubSeqI(char *seqS,int slen,int st,int en,int dir,char *subS)
{
    int i,n;
    char c;

    n = en-st;
    LIMIT_NUM(n,0,slen);
    for(i=0;i<n;i++)
    {
        if(dir==REVERSE)
            c = CompDNABaseC(seqS[en-i-1]);
        else
            c = seqS[st+i];
        subS[i] = c;
    }
    subS[n]='\0';   
    return(n);
}
/***********************************************************************
*   Reduce sequence to only bases "fir" to "las" 
*   If rre is true, reference relative the end rather than start (backwards)
*/
int GetReducedSeqI(char *seqS,int len,int fir,int las,int rre,char *newS)
{
    int i,j,start,end;

    if(rre) {
        start=len-las;
        end=len-fir+1;
    }
    else {
        start=fir-1;
        end=las;
    }
    LIMIT_NUM(start,0,len);
    LIMIT_NUM(end,0,len);
    i = start;
    j = 0;
    while(i<end)
    {
        newS[j++] = seqS[i++];
    }
    newS[j] = '\0';
    return(j);
}
/**************************************************************************
*   Converts any non ACGT to N and all to uppercase
*   If ols is TRUE, one-letter-SNP codes are used for SNPs
*   If mlc is TRUE, mask-lower-case i.e. set lowercase to N's
*/
int CleanUpSeqI(char *inS, int slen, char *outS, int ols, int mlc)
{
    int i,j,s,nn,mn;
    char snpS[DEF_BS+1];

    i = j = 0;
    while(i<slen)
    {
        /***
        *   Masking lowercase?
        */
        if( (mlc) && (islower(INT(inS[i]))) ) {
            inS[i] = 'N';
        }
        switch(inS[i])
        {
            case 'a': case 'A':
            case 'c': case 'C':
            case 'g': case 'G':
            case 't': case 'T':
                outS[j++] = TOUPPER(inS[i++]);
                break;
            case '[': 
                s = mn = nn = 0;
                while( (inS[i] != ']') && (i<slen) && (s<DEF_BS) )
                {
                    snpS[s++] = inS[i++];
                }
                snpS[s] = '\0';
                s = CollapseSnpStringI(snpS,snpS);
                i++;
                /***
                *   If simple single-base, two-allele SNP and one-letter
                */
                if( (ols) && (strlen(snpS)==1) ) {
                    outS[j++] = snpS[0];
                }
                else {
                    while(mn>0)
                    {
                        outS[j++] = 'N';
                        mn--;
                    }
                }
                break;
            case '/': 
            case '*': 
                i++;
                break;
            default:
                if(isalpha(INT(inS[i]))) {
                    outS[j++] = 'N';
                }
                else {
                    outS[j++] = '?';
                }
                i++;
        }
    }
    outS[j] = '\0';
    return(j);
}
