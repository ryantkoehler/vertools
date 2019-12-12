/*
* dset.h
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

# define lmax 1000
# define convergeD 1e-15

# define MIN_TM -273.15 /* Absolute zero (minimum possible temperature */


int thermoacI(double CMoD,double CMgD,double tempD,double CtopD,double CbotD,int seqtop[100], int leffI,double Gprop[4][4][4][4],double Gterm[5][5][5][5], double Hprop[4][4][4][4],double Hterm[5][5][5][5], double HinitAT,double GinitAT, double HinitGC,double GinitGC,double H[4], double Ssalt[4],double Gsalt[5],  double Tm[5],double xD[5]);
int thermocI(double CMoD,double CMgD,double tempD,double CtopD,double CbotD,double Cbot2D,int seqtop[lmax],int seqbotI[lmax], int leffI,double Gprop[4][4][4][4], double Gterm[5][5][5][5],double Hprop[4][4][4][4], double Hterm[5][5][5][5],double HinitAT,double GinitAT, double HinitGC,double GinitGC,double Gcoax[4][4][4][4],double Hcoax[4][4][4][4],int lbotI,int *,double GgapD[30],double GATclosepenD,double H[4], double Ssalt[4],double Gsalt[5],  double Tm[5],double xD[5]);
int thermocpI(double CMoD,double CMgD,double tempD,double CtopD,double CbotD,int seqtop[lmax],int seqbot[lmax], int leff,double Gprop[4][4][4][4], double Gterm[5][5][5][5],double Hprop[4][4][4][4], double Hterm[5][5][5][5],double HinitAT,double HinitGC, double GinitAT,double GinitGC, int *errorI,double H[4], double Ssalt[4],double Gsalt[5],  double Tm[5],double xD[5]);
