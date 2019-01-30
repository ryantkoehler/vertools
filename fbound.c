/*
* fbound.c
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prim.h"
#include "fbound.h"


/*****************************************************************/ 
double Fraction1D(double conc,double conc2,double temp, double gD)
{   
    double aD,bD,cD,qD,eD,xD,climD,cexcD;
    int ns;

    if (conc>conc2)
    {
        climD=conc2;
        cexcD=conc;
    }   
    else
    {
        climD=conc;
        cexcD=conc2;
    }
    eD=exp( (gD* -1000.0) / (1.98722*(temp+273.15)) );
    aD=-eD*climD;
    bD=1.0+eD*climD+eD*cexcD;
    cD=-eD*cexcD;

    ns = NUM_SIGN(bD);
    qD=-0.5*(bD+RNUM(ns)*sqrt(pow(bD,2.0)-4.0*aD*cD));

    xD=qD/aD*100.0;
    return xD;
}
/*****************************************************************/
double Fraction2D(double conc,double conc2,double temp ,double gD)
    {   
    double aD,bD,cD,qD,eD,xD,climD,cexcD;
    int ns;

    if (conc>conc2)
        {climD=conc2;
        cexcD=conc;}   
    else
        {climD=conc;
        cexcD=conc2;}
    eD=exp(-(gD*1000.0)/(1.98722*(temp+273.15)));
    aD=-eD*climD;
    bD=1+eD*climD+eD*cexcD;
    cD=-eD*cexcD;

    ns = NUM_SIGN(bD);
    qD=-0.5*(bD+RNUM(ns)*sqrt(pow(bD,2.0)-4.0*aD*cD));

    xD=cD/qD*100.0;
    return xD;
    }
