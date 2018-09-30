/*
* divinds.c
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
#include <math.h>
#include "prim.h"
#include "divinds.h"

#define DB_STAT     if(DB[17])
#define DB_DIVIN    if(DB[18])

/****************************************************************************
*   nI = Total number of individuals, sI = number of classes 
*/
double BrillouinIndexD(int nI, int *nAI, int sI)
{
    int iI;
    double divD, fD;

    DB_DIVIN DB_PrI(">> BrillouinIndexD nI %d, sI %d\n",nI,sI);
    divD = LogFactorialD(nI);
    DB_DIVIN DB_PrI("+ div %f\n");
    for(iI = 0; iI < sI; iI++)
    {
        fD = LogFactorialD(nAI[iI]);
        divD -= fD;
        DB_DIVIN DB_PrI("+  - %9.5f => div %f\n",fD,divD);
    }
    divD /= DNUM(nI);
    DB_DIVIN DB_PrI("<< BrillouinIndexD %f\n",divD);
    return(divD);
}
/****************************************************************************/
double BrillouinEvenIndexD(int nI, int *nAI, int sI)
{
    double evD, hD, mD;

    DB_DIVIN DB_PrI(">> BrillouinEvenIndexD nI %d, sI %d\n",nI,sI);
    hD = BrillouinIndexD(nI, nAI, sI);
    mD = BrillouinMaxIndexD(nI, sI);
    DB_DIVIN DB_PrI("+ H %f Hmax %f\n",hD,mD);
    evD = hD/mD;
    DB_DIVIN DB_PrI("<< BrillouinEvenIndexD %f\n",evD);
    return(evD);
}
/****************************************************************************/
double BrillouinMaxIndexD(int nI, int sI)
{
    int rI, nsI;
    double divD, rD, d1D, d2D;

    DB_DIVIN DB_PrI(">> BrillouinMaxIndexD nI %d, sI %d\n",nI,sI);
    divD = LogFactorialD(nI);
    nsI = INT(nI/sI);
    rI = nI - (sI * nsI);
    DB_DIVIN DB_PrI("+ div %f, ns %d, r %d\n",divD,nsI,rI);
    rD = LogFactorialD(nsI);
    d1D = rD * DNUM(sI-rI);
    rD = LogFactorialD(nsI+1);
    d2D = rD * DNUM(rI);
    DB_DIVIN DB_PrI("+ d1 %f, d2 %f ",d1D,d2D);
    divD -= (d1D + d2D);
    DB_DIVIN DB_PrI("==> %f\n",divD);
    divD /= DNUM(nI);
    DB_DIVIN DB_PrI("<< BrillouinMaxIndexD  %f\n",divD);
    return(divD);
}
/****************************************************************************
*   Shannon index
*   nI = total number of items (array sum)
*   nAI = array of counts for each of [sI] items
*/
double ShannonWeaverIndexD(int n, int *nAI, int s)
{
    int i;
    double pD,divD;

    divD = 0.0;
    for(i = 0; i < s; i++)
    {
        if(nAI[i] > 0)
        {
            pD = DNUM(nAI[i])/DNUM(n);
            pD *= LOG_2(pD);
/*
printf("[%d] = %d p=%f lg2p=%f pD=%f\n",i,nAI[i], DNUM(nAI[i])/DNUM(n),
    LOG_2(DNUM(nAI[i])/DNUM(n)) ,pD);
*/
            divD += pD;
        }
    }
    if(divD<0.0)
    {
        divD = -divD;
    }
    return(divD);
}
/****************************************************************************/
double ShannonWeMaxIndexD(int sI)
{
    double mD;

    mD = LOG_2(sI);
    return(mD);
}
/****************************************************************************/
double ShannonWeEvenIndexD(int nI, int *nAI, int sI)
{
    double evD, hD, mD;
    
    hD = ShannonWeaverIndexD(nI, nAI, sI);
    mD = ShannonWeMaxIndexD(sI);
    evD = hD/mD;
    return(evD);
}
/****************************************************************************
*   Shannon index
*   Log 10 version
*   nI = total number of items (array sum)
*   nAI = array of counts for each of [sI] items
*/
double ShannonWeaver10IndexD(int n, int *nAI, int s)
{
    int i;
    double pD,divD;

    divD = 0.0;
    for(i = 0; i < s; i++)
    {
        if(nAI[i] > 0 )
        {
            pD = DNUM(nAI[i])/DNUM(n);
            pD *= LOG_10(pD);
            divD += pD;
        }
    }
    if(divD<0.0)
    {
        divD = -divD;
    }
    return(divD);
}
/****************************************************************************
*   Log 10 version
*/
double ShannonWe10MaxIndexD(int sI)
{
    double mD;

    mD = LOG_10(sI);
    return(mD);
}
/****************************************************************************
*   Log 10 version
*/
double ShannonWe10EvenIndexD(int nI, int *nAI, int sI)
{
    double evD, hD, mD;
    
    hD = ShannonWeaver10IndexD(nI, nAI, sI);
    mD = ShannonWe10MaxIndexD(sI);
    evD = hD/mD;
    return(evD);
}
/****************************************************************************
*   Simpson's diversity index 
*/
double SimpsonIndexD(int n, int *nAI, int s)
{
    int i;
    double divD, pD, numD, denD;

    DB_DIVIN DB_PrI(">> SimpsonIndexD n %d, s %d ",n,s);
    denD = DNUM(n)*DNUM(n-1);
    DB_DIVIN DB_PrI("den %3.3f\n",denD);
    divD = 0.0;
    for(i = 0; i < s; i++)
    {
        numD = DNUM(nAI[i])*DNUM(nAI[i]-1);
        pD = numD/denD;
        divD += pD;
    }
    DB_DIVIN DB_PrI("<< SimpsonIndexD %3.3f\n",divD);
    return(divD);
}
