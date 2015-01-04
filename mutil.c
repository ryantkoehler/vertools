/*
* mutil.c
*
* Copyright 2015 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <math.h>
#include <time.h>
#include "prim.h"

#define DB_RMASK    if(DB[16])

/*********************************************************************** rrr
*   Returns the value of base to the exponent power
*/
double PowD(double baseR, double expR)
{
    double rD;

    if(baseR <= 0.0) {
        return(0.0);
    }
    rD = expR * log(baseR);
    rD = exp(rD);
    return(rD);
}
/************************************************************************
*   Seed random number; if negative seed, set with time
*/
void Srand(int seed)
{   
    if(seed<=0) {   
        SRAND(time(NULL));  
    }
    else {  
        SRAND(seed);        
    }
}
/************************************************************************/
int RandI(int range)
{
    int r;

    r = INT( RNUM(RAND()) / RNUM(RAND_DEN) * RNUM(range) );
    return(r);
}
/************************************************************************/
REAL RandR(REAL rangeR)
{
    REAL rR;

    rR = RNUM(RAND()) / RNUM(RAND_DEN) * rangeR;
    return(rR);
}
/**************************************************************************/
DOUB RandD(DOUB rangeD)
{
    DOUB rD;

    rD = DNUM(RAND()) / DNUM(RAND_DEN) * rangeD;
    return(rD);
}
/**************************************************************************/
void FillRandSeedString(int seed, char *sS)
{
    if(seed<=0) {
        sprintf(sS,"clock");
    }
    else {
        sprintf(sS,"%d",seed);
    }
}
/****************************************************************************
*   Fill passed array, seq[len], with random ints 0 to max
*   If unique is true, no two random numbers should be the same
*/
int SetRandSequenceI(int *seq,int len,int max,int unique)
{
    int i,u,d,ntake;

    if((unique)&&(max < len)) {
        printf("Bad unique sequence paremeters: len %d, max %d\n",len,max);
        ERR("SetRandSequenceI","impossible parameters");
    }
    ntake = 0;
    while(ntake < len)
    {
        d = RandI(max);
        if(unique) {
            u = TRUE;
            for(i=0;i<ntake;i++)
            {
                if(d==seq[i]) {
                    u = FALSE;  break;
                }
            }
            if(u==TRUE) {
                seq[ntake++] = d;
            }
        }
        else {
            seq[ntake++] = d;
        }
    }
    return(ntake);
}
/****************************************************************************/
#define MAX_REAL_FACT_NUM   150
/****************************************************************************/
DOUB LogFactorialD(int nI)
{
    double fD, dD;

    if(nI <= MAX_REAL_FACT_NUM) {    
        fD = FactorialD(nI);
        dD = log10(fD);
    }
    else {
        dD = SterlingLogFactorialD(nI);
    }
    return(dD);
}
/****************************************************************************/
DOUB FactorialD(int nI)
{
    double dD, fD;

    dD = DNUM(nI);
    if(nI <= 1) {
        dD = 1.0;
    }
    else {    
        fD = FactorialD(nI-1);
        dD *= fD;
    }
    return(dD);
}
/****************************************************************************/
DOUB SterlingLogFactorialD(int nI)
{
    double dD;

    dD = DNUM(nI) + 0.5;
    dD *= log10(DNUM(nI));
    dD -= (DNUM(nI) * 0.434294482);
    dD += 0.39909;
    return(dD);
}
/***********************************************************************
*   Gaussian random number generator;
*   mR = mean, sR = standard deviation
*
*   Based on code from boxmuller.c (on net)
*   Implements the Polar form of the Box-Muller Transformation
*   (c) Copyright 1994, Everett F. Carter Jr.
*   Permission is granted by the author to use this software for any 
*   application provided this copyright notice is preserved.
*/
DOUB RandGaussD(REAL mD,REAL sD)
{                       
    DOUB x1D, x2D, wD, y1D;
    static DOUB y2D;
    static int use_last = 0;

    if (use_last)               /* se value from previous call */
    {
        y1D = y2D;
        use_last = 0;
    }
    else
    {
        do 
        {
            x1D = DNUM(RAND()) / DNUM(RAND_DEN);
            x2D = DNUM(RAND()) / DNUM(RAND_DEN);
            x1D = (x1D * 2.0) - 1.0;
            x2D = (x2D * 2.0) - 1.0;
            wD = x1D * x1D + x2D * x2D;
        } while ( wD >= 1.0 );

        wD = sqrt( (-2.0 * log( wD ) ) / wD );
        y1D = x1D * wD;
        y2D = x2D * wD;
        use_last = 1;
    }
    return( mD + y1D * sD );
}
