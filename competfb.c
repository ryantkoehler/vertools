/*
* competfb.c
*
* Copyright 2014 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prim.h"
#include "competfb.h"


/****************************************************************************
*	The input consists of: C0A C0B1    C0B2    DGAB1T  DGAB2T  T
*	C0A is the initial concentration of the strand A that is going to 
*		be competed for (genomic DNA in case of SNP).
*	C0B1 is the initial concentration of the strand (B1) that is going
*		to compete with B2 to bind A (FAM probe in case of SNP).
*   C0B2 is the initial concentration of the strand (B2) that is going
*		to compete with B1 to bind A (VIC probe in case of SNP).
*	DGAB1T is the free energy at temperature T associated with the
*		reaction A +B1 =AB1.
*	DGAB2T is the free energy at temperature T associated with the
*		reaction A +B2 =AB2.
*	T is the temperature.
*
*	The output is set to f1PD and f1PD
*	f1PD is the fraction duplex of AB1
*	f2PD is the fraction duplex of AB2
*/
int CalcCompEqI(DOUB C0A, DOUB C0B1, DOUB C0B2, DOUB DGAB1T, DOUB DGAB2T,
	DOUB T, DOUB *f1PD, DOUB *f2PD)
{
	DOUB A,B1,B2,KAB1,KAB2,Aprev,B1prev,B2prev,fAB1,fAB2,conv;

		/*Calculate equilibrium constants*/
	
		KAB1=exp(-(DGAB1T*1000)/(1.98722*(T+273.15)));
		KAB2=exp(-(DGAB2T*1000)/(1.98722*(T+273.15)));
				
		/*Define starting conditions*/
		
		A=C0A;
		B1=C0B1;
		B2=C0B2;
		Aprev=0.0;
		B1prev=0.0;
		B2prev=0.0;
	
		/*Iterate until convergence reached*/
			
		conv=1.0e-15;

		while ((fabs(A-Aprev)>conv) || (fabs(B1-B1prev)>conv) || (fabs(B2-B2prev)>conv))
			{
			Aprev=A;
			B1prev=B1;
			B2prev=B2;
			A=C0A/(1+KAB1*B1+KAB2*B2);
			B1=C0B1/(1+KAB1*A);
			B2=C0B2/(1+KAB2*A);
			}
			
		/*Calculate fraction duplex*/
		
		fAB1=KAB1*A*B1/MIN_NUM(C0A,C0B1);
		fAB2=KAB2*A*B2/MIN_NUM(C0A,C0B2);

	*f1PD = fAB1;
	*f2PD = fAB2;
	return(TRUE);
}
