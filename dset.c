/*
* dset.c
*
* Copyright 2016 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include <time.h>
#include "dset.h"

/***
*   Conditional debug printing
*
*   In order to use dynamic debug, the global variable char DB[] must
*   be defined (as it is in prim.h). If prim.h is included (and the
*   variable is initialized safely, then __PRIMH__ is defined and DB[]
*   may be used.
*
*   If you don't want prim.h, comment it out and define DB_DSET to
*   be ON (i.e. if(1) = all printing) or OFF (i.e. if(0) = no printing)
*/
#include "prim.h"
#ifdef __PRIMH__

    #define DB_DSET if(DB[108])

#else

    #define DB_DSET if(1)

#endif


/*************************************************************/
/*thermoacI*/
/*************************************************************/
int thermoacI(double CMoD,double CMgD,double tempD,double CtopD,double CbotD ,int seqtop[lmax], int leff,double Gprop[4][4][4][4],double Gterm[5][5][5][5], double Hprop[4][4][4][4], double Hterm[5][5][5][5],double HinitAT,double HinitGC, double GinitAT,double GinitGC,double H[4], double Ssalt[4],double Gsalt[5],  double Tm[5],double xD[5])
    {
    int scomp;
    int initright, initleft;
    double G37, Gtemp,   S;
    double saltlen;
    double saltcorrG(double Gtemp,   double CMo, double CMg, double saltlen,double tempD);
    double saltcorrS(double H, double G37, double CMo, double CMg, double saltlen);
    double tmcalc(double Ctop,double Cbot,double Ssalt,double H,  int scomp);
    double hybtemp(double H, double S, double temp);
    double FractionD(double CtopD,double CbotD,double tempD ,double Gsalt);
    int selfcoacI( int seqtopI[lmax], int leffI);
    double inittermprop37acD(int seqtop[lmax],int initright, int initleft, double Eprop[4][4][4][4], double Eterm[5][5][5][5], double EinitAT,double EinitGC, int leff, int scomp);

    DB_DSET { printf(">> thermoacI\n"); fflush(stdout); }
    G37=0.0;
    H[0]=0.0;   
    scomp=selfcoacI(seqtop,leff);
    initright=leff-1;
    initleft=0;
    /*
    printf("init # HAT%f HGC%f GAT%f GGC%f",HinitAT,HinitGC,GinitAT,GinitGC);
    */
    H[0]=inittermprop37acD(seqtop, initright, initleft, Hprop, Hterm, HinitAT,HinitGC, leff, 0);
    G37=inittermprop37acD(seqtop,initright, initleft,Gprop, Gterm,GinitAT,GinitGC, leff, scomp);
    /*
    printf("ac after inittermprop37 Gtemp %f %f ",G37, H[0]);
    */
    S=((H[0]-G37)/310.15)*1000.0;
    Gtemp=hybtemp(H[0], S, tempD);
    /*
    printf("Gtemp %f  ",Gtemp);
    */
    saltlen=leff-1;
    Gsalt[0]=saltcorrG(Gtemp,CMoD,CMgD,saltlen,tempD);
    /*
    printf("Gsalt %f ",Gsalt[0]);
    */
    Ssalt[0]=saltcorrS( H[0],  G37, CMoD, CMgD, saltlen);
    /*
    printf("Ssalt %f",*Ssalt);
    */
    Tm[0]=tmcalc( CtopD, CbotD, Ssalt[0], H[0],scomp);
    /*
    printf("Tm %f",Tm[0]);
    */
    xD[0]=FractionD(CtopD,CbotD,tempD,Gsalt[0]);
    /**
    printf("&nbsp;&nbsp;&nbsp;%10.1f&nbsp;&nbsp;&nbsp;%10.1f&nbsp;&nbsp;&nbsp;%11.2f&nbsp;&nbsp;&nbsp;%10.1f&nbsp;&nbsp;&nbsp;%f\n",*H,*Ssalt,*Gsalt,*Tm,*xD);
    **/
    return 0;
    }
/*************************************************************/
/*thermocI*/
/*************************************************************/
int thermocI(double CMoD,double CMgD,double tempD,double CtopD,double CbotD,double Cbot2D,int seqtop[lmax],int seqbot[lmax], int leff,double Gprop[4][4][4][4], double Gterm[5][5][5][5],double Hprop[4][4][4][4], double Hterm[5][5][5][5],double HinitAT,double HinitGC, double GinitAT,double GinitGC,double Gcoax[4][4][4][4],double Hcoax[4][4][4][4], int lbotI,int *errorI,double GgapD[30],double GATclosepenD,double H[4], double Ssalt[4],double Gsalt[5],  double Tm[5],double xD[5])
    {
    int gaposI,gaplenI;
    int coaxposI,scomp=0;
    int initright, initleft;
    double G37, Gtemp,  S,GgapcorD,GglobalD,HglobalD,SglobalD,saltlen1D,saltlen2D;
    double saltlen,GcoaxcorD,HcoaxcorD;
    int seqbot35[lmax];
    double saltcorrG(double Gtemp,   double CMoD, double CMgD, double saltlen,double tempD);
    double saltcorrS(double H, double G37, double CMoD, double CMgD, double saltlen);
    double tmcalc(double CtopD,double CbotD,double Ssalt,double H,  int scomp);
    double hybtemp(double H, double S, double temp);
    int gapstack(int seqtop[lmax],int seqbot35[lmax], int *,int lbotI ,int *,int *,int *,int * );
    double coaxcalc(double EcoaxD[4][4][4][4],int seqtop[lmax] , int seqbot35[lmax],int *);
    double gapcalcD(int gaposI,int gaplenI,double EgapD[30],double ATclospenD,double Eprop[4][4][4][4],int seqtop[lmax], int seqbot35[lmax]);
    double FractionD(double CtopD,double CbotD,double tempD ,double Gsalt);
    int selfco(int seqtop[lmax], int seqbot[lmax], int leff);
    int penleft(int seqtop[lmax], int seqbot[lmax], double Gprop[4][4][4][4], double Gterm[5][5][5][5],double GinitAT,double GinitGC, int leff);
    int penright(int seqtop[lmax], int seqbot[lmax], double Gprop[4][4][4][4], double Gterm[5][5][5][5],double GinitAT,double GinitGC, int leff);
    int saltlength(int seqtop[lmax],int seqbot[lmax],int initright, int initleft, int leff);
    double inittermprop37(int seqtop[lmax], int seqbot[lmax],int initright, int initleft, double Eprop[4][4][4][4], double Eterm[5][5][5][5], double EinitAT,double EinitGC, int leff, int scomp);
    void bot35(int seqbot[lmax],int seqbot35[lmax],int leff);
    double GnetcalcD(double tempD,double HglobalD,double SglobalD,double CtopD,double CbotD,double Cbot2D,double Gsalt[5],double Ssalt[4],double H[4]); 
    double TmnetcalcD(double HglobalD,double SglobalD,double CtopD,double dCbotD,double Cbot2D,double Gsalt[5],double Ssalt[4],double H[4]);
    double minD(double aD,double bD,double cD);
    double dedup2GD,dedup2HD,dedup1GD,dedup1HD,transD;

    DB_DSET {printf(">> thermocI\n"); fflush(stdout);   }
    G37=0.0;
    *H=0.0;
    GcoaxcorD=0.0;
    HcoaxcorD=0.0;
    GgapcorD=0.0;
    coaxposI=0;
    *errorI=0;
    /**
    printf("ThermocI entry");
    printf("Cbot2D=%f<Br>",Cbot2D);
    printf("error   %d",*errorI);
    fflush(stdout);
    **/ 
    bot35(seqbot,seqbot35,lbotI);
    gapstack(seqtop,seqbot35,&leff,lbotI,errorI,&coaxposI,&gaposI,&gaplenI);
    /**
    printf("coaxposI%d",coaxposI);
    printf("leffI=%d",leff);
    **/
    if (*errorI==0)
        {
        /*printf("-------------------");*/
        fflush(stdout);
        /*----------if coaxial stacking*/
        if(coaxposI!=0)
            {
            /*----------Attribution of concentrations*/
            /*Rules:CtopD never contains a nick*/
            /* CbotD always nicked at 5' extremity*/
            /* CbotD always nicked at 3'extremity*/
            if (coaxposI<0)
                {
                transD=CbotD;
                CbotD=Cbot2D;
                Cbot2D=transD;  
                }
            /*----------*/
            /*----------Calculation of coaxial stacking contribution*/
            GcoaxcorD=coaxcalc( Gcoax,seqtop, seqbot35,&coaxposI);
            HcoaxcorD=coaxcalc( Hcoax,seqtop, seqbot35,&coaxposI);          
            /*----------*/
            /*----------Calculation of duplex 0*/
            /*Duplex 0 thermodynamics is thermodynamics of Duplex 0 binding alone */
            /*printf("<Br>-------------------------------------------------duplex 0<Br>");*/        
            initleft=penleft(seqtop,seqbot35,Gprop,Gterm,GinitAT,GinitGC,abs(coaxposI)+1);
            initright=abs(coaxposI);
            /*printf("initright %d initleft %d",initright,initleft);*/
            /*0.5 added to correct for extra DE*/
            saltlen1D=saltlength(seqtop,seqbot35,initright,initleft,initright+1)+0.5;
/* Moved up to from below RTK */
            saltlen=saltlength(seqtop,seqbot35,initright, initleft,leff);
            /*--------------------Determination of extra DE contribution as a functuion of top or bottom coax*/
            if (coaxposI>0)
                {
                dedup1GD=Gterm[seqtop[coaxposI]][seqtop[coaxposI+1]][4][seqbot35[coaxposI]];
                dedup1HD=Hterm[seqtop[coaxposI]][seqtop[coaxposI+1]][4][seqbot35[coaxposI]];
                }
            else
                {
                dedup1GD=Gterm[seqtop[abs(coaxposI)]][4][seqbot35[abs(coaxposI)+1]][seqbot35[abs(coaxposI)]];
                dedup1HD=Hterm[seqtop[abs(coaxposI)]][4][seqbot35[abs(coaxposI)+1]][seqbot35[abs(coaxposI)]];
                } 
            /*printf("<Br>dedup1GD %f,dedup1HD %f <Br>",dedup1GD,dedup1HD);*/
            /*fflush(stdout);*/ 
            /*--------------------*/
            
            /*----------Calculate gap contribution if gap*/
            if ((gaposI!=0)&&(gaposI<abs(coaxposI)))
                {
                GgapcorD=gapcalcD(gaposI,gaplenI,GgapD,GATclosepenD,Gprop,seqtop,seqbot35);
                G37=G37+GgapcorD;
                saltlen=saltlen-gaplenI;
                }
            /*----------*/
            H[0]=inittermprop37(seqtop, seqbot35,initright, initleft, Hprop, Hterm, HinitAT,HinitGC, initright+1, 0)+dedup1HD;
            G37=G37+inittermprop37(seqtop, seqbot35,initright, initleft,Gprop, Gterm,GinitAT,GinitGC, initright+1, 0)+dedup1GD;
            /*printf("<Br>G37=%f**",G37);*/
            S=((H[0]-G37)/310.15)*1000.0;
            Gtemp=hybtemp(H[0], S, tempD);  
            Gsalt[0]=saltcorrG( Gtemp,   CMoD, CMgD, saltlen1D,tempD);
            Ssalt[0]=saltcorrS( H[0],  G37, CMoD, CMgD, saltlen1D);
            Tm[0]=tmcalc( CtopD, CbotD, Ssalt[0], H[0],scomp);
            xD[0]=FractionD(CtopD,CbotD,tempD,Gsalt[0]);
            /**/
            /*printf("<Br>1 H S G T X %f %f %f %f %f<Br>",H[0],Ssalt[0],Gsalt[0],Tm[0],xD[0]);*/
            /**/
            /*----------*/
            /*----------calculation of duplex 1*/
            /*Duplex 1 thermodynamics is thermodynamics of Duplex 0 binding assuming other duplex is 100% bound */
            /*printf("<Br>-------------------------------------------------duplex 1<Br>");*/
            /*printf("<Br>GcoaxcorD %f, HcoaxcorD %f",GcoaxcorD,HcoaxcorD);*/
            /*No dangling end stabilization*/
            H[1]=H[0]+HcoaxcorD-dedup1HD;
            G37=G37+GcoaxcorD-dedup1GD;
            S=((H[1]-G37)/310.15)*1000.0;
            Gtemp=hybtemp(H[1], S, tempD);
            /*No dangling end stabilization => -0.5*/
            saltlen1D=saltlen1D-0.5;    
            Gsalt[1]=saltcorrG( Gtemp,   CMoD, CMgD, saltlen1D,tempD);
            Ssalt[1]=saltcorrS( H[1],  G37, CMoD, CMgD, saltlen1D);
            Tm[1]=tmcalc( CtopD, CbotD, Ssalt[1], H[1],scomp);
            xD[1]=FractionD(CtopD,CbotD,tempD,Gsalt[1]);
            /**
            printf("<Br>2 H S G T X %f %f %f %f %f<Br>",H[1],Ssalt[1],Gsalt[1],Tm[1],xD[1]);
            **/
            /*----------*/
            /*----------calculation of duplex 2*/
            /*Duplex 2 thermodynamics is thermodynamics of Duplex 2 binding alone */
            /*printf("<Br>-------------------------------------------------duplex 2<Br>");*/
            G37=0.0;
            initleft=abs(coaxposI)+1;   
            initright=penright(seqtop,seqbot35,Gprop,Gterm,GinitAT,GinitGC,leff);
            /*leff for this duplex in saltlength is initright-initleft+1 - initleft is 0 - initright is initright-initleft*/
            /*This is equivalent to placing the duplex so that is starts at position 0*/
            /*0.5 added to correct for extra DE*/
            saltlen2D=saltlength(seqtop,seqbot35,initright-initleft, 0,initright-initleft+1)+0.5;
            /*--------------------Determination of extra DE contribution as a functuion of top or bottom coax*/
            if (coaxposI>0)
                {
                dedup2GD=Gterm[seqtop[coaxposI]][seqtop[coaxposI+1]][seqbot35[coaxposI+1]][4];
                dedup2HD=Hterm[seqtop[coaxposI]][seqtop[coaxposI+1]][seqbot35[coaxposI+1]][4];
                }
            else
                {
                dedup2GD=Gterm[4][seqtop[abs(coaxposI)+1]][seqbot35[abs(coaxposI)+1]][seqbot35[abs(coaxposI)]];
                dedup2HD=Hterm[4][seqtop[abs(coaxposI)+1]][seqbot35[abs(coaxposI)+1]][seqbot35[abs(coaxposI)]];
                } 
            /*printf("<Br>dedup2GD %f,dedup2HD %f <Br>",dedup2GD,dedup2HD);*/   
            /*--------------------*/
            
            /*----------Calculate gap contribution if gap*/
            /*printf("<Br>initright=%d, initleft=%d",initright,initleft);*/
            if ((gaposI!=0)&&(gaposI>abs(coaxposI)))
                {
                GgapcorD=gapcalcD(gaposI,gaplenI,GgapD,GATclosepenD,Gprop,seqtop,seqbot35);
                G37=G37+GgapcorD;
                saltlen=saltlen-gaplenI;
                
                }
            /*----------*/
            /*printf("<Br>initright=%d, initleft=%d",initright,initleft);*/
            /*-1.001 is needed because as initleft<>0, terminal elts are added in inittermprop37*/ 
            /*as the terminal elts are in fact perfect matches they are found in Gterm and Hterm as 1.001 filling default*/
            H[2]=inittermprop37(seqtop, seqbot35,initright, initleft, Hprop, Hterm, HinitAT,HinitGC, leff, 0)-1.001+dedup2HD;
            G37=G37+inittermprop37(seqtop, seqbot35,initright, initleft,Gprop, Gterm,GinitAT,GinitGC, leff, 0)-1.001+dedup2GD;
            S=((H[2]-G37)/310.15)*1000.0;
            Gtemp=hybtemp(H[2], S, tempD);      
            Gsalt[2]=saltcorrG( Gtemp,   CMoD, CMgD, saltlen2D,tempD);
            Ssalt[2]=saltcorrS( H[2],  G37, CMoD, CMgD, saltlen2D);
            Tm[2]=tmcalc( CtopD, Cbot2D, Ssalt[2], H[2],scomp);
            xD[2]=FractionD(CtopD,Cbot2D,tempD,Gsalt[2]);
            /**
            printf("<Br>3 H S G T X %f %f %f %f %f<Br>",H[2],Ssalt[2],Gsalt[2],Tm[2],xD[2]);
            **/
            /*----------*/
            /*----------calculation of second duplex 3*/
            /*Duplex 3 thermodynamics is thermodynamics of Duplex 2 binding assuming other duplex is 100% bound */
            /*printf("<Br>-------------------------------------------------duplex 3<Br>");*/
            /*No dangling end stabilization*/
            H[3]=H[2]+HcoaxcorD-dedup2HD;
            G37=G37+GcoaxcorD-dedup2GD;
            S=((H[3]-G37)/310.15)*1000.0;
            Gtemp=hybtemp(H[3], S, tempD);
            /*No dangling end stabilization => -0.5*/
            saltlen2D=saltlen2D-0.5;    
            Gsalt[3]=saltcorrG( Gtemp,   CMoD, CMgD, saltlen2D,tempD);
            Ssalt[3]=saltcorrS( H[3],  G37, CMoD, CMgD, saltlen2D);
            Tm[3]=tmcalc( CtopD, Cbot2D, Ssalt[3], H[3],scomp);
            xD[3]=FractionD(CtopD,Cbot2D,tempD,Gsalt[3]);
            /**
            printf("<Br>4 H S G T X %f %f %f %f %f<Br>",H[3],Ssalt[3],Gsalt[3],Tm[3],xD[3]);
            **/
            /*----------*/
            /*----------calculation of global structure*/
            /*printf("<Br>-------------------------------------------------triplex <Br>");*/
            GglobalD=Gsalt[1]+Gsalt[2]-dedup2GD;
            HglobalD=H[1]+H[2]-dedup2HD;
            SglobalD=((HglobalD-GglobalD)/(tempD+273.15))*1000.0;
                    /**
            printf("Hglobal %f Sglobal %f Gglobal %f",HglobalD,SglobalD,GglobalD);
            **/
            /*-----------*/
            /*-----------calculation of net free energy of global structure*/
            xD[4]=GnetcalcD(tempD,HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);
            /**
            printf("X[4] %f",xD[4]);
            printf("fraction: %f",FractionD(CtopD,CbotD,tempD,Gsalt[4]));
            w=0;
            printf("fraction: %f",FractionD(CtopD,CbotD,tempD,w));
            **/
            /*-----------*/
            /*-----------calculation of Tm effective of global structure*/
            Tm[4]=TmnetcalcD(HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);
            /**
            printf("<Br>Tm[4] %f <Br>",Tm[4]);
            **/
            /*-----------*/
            }
        /*----------if no coaxial stacking*/
        else if (coaxposI==0)
            {
            DB_DSET {printf("+ no-coax case\n"); fflush(stdout);    }
            Gsalt[1]=698282.0;  
            scomp=selfco( seqtop, seqbot35,leff);   
            initleft=penleft(seqtop,seqbot35,Gprop,Gterm,GinitAT,GinitGC, leff);    
            initright=penright(seqtop,seqbot35,Gprop,Gterm,GinitAT,GinitGC,leff);
/* Moved up to from below RTK */
            saltlen=saltlength(seqtop,seqbot35,initright, initleft,leff);
            /**
            printf("main after pen %d %d", initleft,initright);
            fflush(stdout);
            **/
            if (initright<=initleft)
                {
                /*----------uses double mismatch CCCC values for calculation of H and G if no consecutive pairing*/
                *H=leff*Hprop[1][1][1][1];
                G37=leff*Gprop[1][1][1][1];
                /*----------*/
                }
            else
                {       
                /*
                printf("saltlen %f  ",saltlen);
                */
                *H=inittermprop37(seqtop, seqbot35,initright, initleft, Hprop, Hterm, HinitAT,HinitGC, leff, 0);
                G37=inittermprop37(seqtop, seqbot35,initright, initleft,Gprop, Gterm,GinitAT,GinitGC, leff, scomp);
                /**
                printf("main after inittermprop37 Gtemp %f %f ",G37, *H);
                fflush(stdout);
                **/
                /*----------Calculate gap contribution if gap*/
                if (gaposI!=0)
                    {
                    GgapcorD=gapcalcD(gaposI,gaplenI,GgapD,GATclosepenD,Gprop,seqtop,seqbot35);
                    G37=G37+GgapcorD;
                    saltlen=saltlen-gaplenI;
                    }
                DB_DSET {printf("+ gap calculated\n"); fflush(stdout); }
                /*----------*/
                }           
            S=((H[0]-G37)/310.15)*1000.0;
            Gtemp=hybtemp(H[0], S, tempD);
/* Wrong place? Moved up to top RTK
            saltlen=saltlength(seqtop,seqbot35,initright, initleft,leff);
*/
            /*
            printf("Gtemp %f  ",Gtemp);
            */
            Gsalt[0]=saltcorrG( Gtemp,   CMoD, CMgD, saltlen,tempD);
            /*
            printf("Gsalt %f ",*Gsalt);
            */
            Ssalt[0]=saltcorrS( H[0],  G37, CMoD, CMgD, saltlen);
            /*
            printf("Ssalt %f",Ssalt[1]);
            */
            *Tm=tmcalc( CtopD, CbotD, *Ssalt, H[0],scomp);
            if(*Tm<MIN_TM)
                {
                *Tm = MIN_TM;
                }
            /*
            printf("Tm %f",*Tm);
            */
            xD[0]=FractionD(CtopD,CbotD,tempD,Gsalt[0]);
                
            
            }
        }
    
    DB_DSET {printf("<< thermocI 0\n"); fflush(stdout); }
    return 1;
    }

/*************************************************************/
/*Calculates Equilibrium concentrations in case of coaxial stacking, GnetcalcD and fraction triplex*/
/*************************************************************/
/* he following equilibria are considered:*/
/* +B=AB   A+C=AC   A+B+C=ABC */
/*instant concentrations of A, B, and C are CAD,CBD,and CCD*/
/*initial concentrations of A, B, and C are CtopD, CbotD, and Cbot2D)*/
/*Ks for the above equilibria are respectively KABD, KACD, and KABCD*/
/*H[0], Ssalt[0] correspond to KABD*/
/*H[2], Ssalt[2] correspond to KACD*/
/*HglobalD, SglobalD correspond to KABCD*/
double GnetcalcD(double tempD,double HglobalD,double SglobalD,double CtopD,double CbotD,double Cbot2D,double Gsalt[5],double Ssalt[4],double H[4])
    {
    const double R=1.98722;
/* TK 12/29/12
    double ABD, ACD;
*/
    int i;
    double minD(double aD,double bD,double cD);
    double KABD,KACD,KABCD,CAD,CBD,CCD,ABCD,CAprevD,CBprevD,CCprevD,KnetD,xnetD;
    i=0;
    /*printf("<Br>Ssalt[0] %f H[0] %f tempD %f R %f",Ssalt[0],H[0],tempD,R);*/
    /*----------Calculates Ks at the given temp*/
    KABD=exp(-(1000*H[0]-(tempD+273.15)*Ssalt[0])/(R*(tempD+273.15)));
    /*printf("KABD %f",KABD);*/
    /*printf("<Br>Ssalt[2] %f H[2] %f tempD %f R %f",Ssalt[2],H[2],tempD,R);*/
    KACD=exp(-(1000*H[2]-(tempD+273.15)*Ssalt[2])/(R*(tempD+273.15)));
    KABCD=exp(-(1000*HglobalD-(tempD+273.15)*SglobalD)/(R*(tempD+273.15)));
    /*printf("<Br>//////////////////////////<Br>");*/
    /*printf("<Br> K %g %g %g",KABD,KACD,KABCD);*/
    /*printf("<Br>C %f %f %f",CtopD,CbotD,Cbot2D);*/
    /*----------*/
    /*----------initialize SS concentrations*/
    CAD=CtopD;
    CBD=CbotD;
    CCD=Cbot2D; 
    /*printf("tempD %f ",tempD);*/
    /*----------*/
    /*----------iterates SS concentrations until convergence is reached*/
    do 
        {
        i=i+1;
        CAprevD=CAD;
        CBprevD=CBD;
        CCprevD=CCD;        
        CAD=CtopD/(1.0+KABD*CBD+KACD*CCD+KABCD*CBD*CCD);
        CBD=CbotD/(1.0+KABD*CAD+KABCD*CAD*CCD);
        CCD=Cbot2D/(1.0+KACD*CAD+KABCD*CAD*CBD);
        /**                     
        printf("///////%g",fabs(CAprevD-CAD));  
        **/
        }
    while ((fabs(CAprevD-CAD)>convergeD) ||(fabs(CBprevD-CBD)>convergeD) ||(fabs(CCprevD-CCD)>convergeD));
    /*----------*/
    /*----------calculates equilibrium concentration of duplexes and triplex*/
    ABCD=KABCD*CAD*CBD*CCD; 
/* TK 12/29/12
    ABD=KABD*CAD*CBD;
    ACD=KACD*CAD*CCD;
*/
    /*printf("<Br>ABCD ABD ACD %g %g %g" ,ABCD,ABD,ACD);*/
    /*printf("<Br> CAD CBD CCD %g %g %g" ,CAD,CBD,CCD);*/
    /*printf("cons %g %g %g",ABD+ACD+ABCD+CAD,ABD+ABCD+CBD,ACD+ABCD+CCD);*/
    /*----------*/
    /*----------Calculates Knet and Gnet*/
    KnetD=ABCD/((CtopD-ABCD)*(CbotD-ABCD)*(Cbot2D-ABCD));
    /*printf("Knet %g",KnetD);*/
    Gsalt[4]=(-R*(tempD+273.15)*log(KnetD))/1000.0;
    /*printf("Gnet %g ",Gsalt[4]);*/
    /*---------*/
    /*----------calculates fraction triplex*/
    /*printf("min %f",minD(CtopD,CbotD,Cbot2D));*/
    xnetD=ABCD/minD(CtopD,CbotD,Cbot2D)*100.0;
    /*printf("<Br>i %d xnet %f",i,xnetD);*/
    /*printf("<Br>//////////////////////////<Br>");*/
    /*----------*/
    return xnetD;
    }
/*************************************************************/
/*Takes minimum of 3 doubles*/
/*************************************************************/
double minD(double aD,double bD,double cD)  
    {
    double miniD;
    if (aD>bD)
        {
        miniD=bD;
        }
    else
        {
        miniD=aD;
        }
    if (miniD>cD)
        {
        miniD=cD;
        }
    return miniD;
    }
/*************************************************************/
/*Calculates the TMnet in case of coaxial stacking*/
/*************************************************************/
double TmnetcalcD(double HglobalD,double SglobalD,double CtopD,double CbotD,double Cbot2D,double Gsalt[5],double Ssalt[4],double H[4])
    {
    double tlowD,thighD,tmidD,flowD,fhighD,fmidD,TmnetD;
    tlowD=0.0;
    thighD=100.0;
    /** 
    printf("<Br>*-*-*-*-*-*- %f %f %f %f %f %f %f %f<Br>",HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Ssalt[2],Ssalt[0],H[0]);
    **/
    flowD=GnetcalcD(tlowD,HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);
    fhighD=GnetcalcD(thighD,HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);
    
    /*printf("<Br>-------<Br>flow %f fhigh %f",flowD,fhighD);*/
    if ((fhighD<50.0) && (flowD>50.0))
        {
        do
            {
            tmidD=(tlowD+thighD)/2;
            fmidD=GnetcalcD(tmidD,HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);  
            if (fmidD>50.0)
                {
                tlowD=tmidD;
                /*fhighD=GnetcalcD(thighD,HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);*/
                }
            else    
                {
                thighD=tmidD;
                /*flowD=GnetcalcD(tlowD,HglobalD,SglobalD,CtopD,CbotD,Cbot2D,Gsalt,Ssalt,H);*/
                }
            /*printf("<Br>t values %f %f %f",tlowD, thighD,tmidD);*/
            }
        while (fabs(tlowD-thighD)>0.01);
        TmnetD=tlowD;
        }
    else
        {
        TmnetD=698282.0;
        }
    return TmnetD;
    }
/*************************************************************/
/*thermocpI*/
/*************************************************************/
int thermocpI(double CMoD,double CMgD,double tempD,double CtopD,double CbotD,int seqtop[lmax],int seqbot[lmax], int leff,double Gprop[4][4][4][4], double Gterm[5][5][5][5],double Hprop[4][4][4][4], double Hterm[5][5][5][5],double HinitAT,double HinitGC, double GinitAT,double GinitGC, int *errorI,double H[4], double Ssalt[4],double Gsalt[5],  double Tm[5],double xD[5])
    {
    int scomp=0;
    int initright, initleft;
    double G37, Gtemp,  S;
    double saltlen;
    int seqbot35[lmax];
    double saltcorrG(double Gtemp,   double CMoD, double CMgD, double saltlen,double tempD);
    double saltcorrS(double H, double G37, double CMoD, double CMgD, double saltlen);
    double tmcalc(double CtopD,double CbotD,double Ssalt,double H,  int scomp);
    double hybtemp(double H, double S, double temp);
    double FractionD(double CtopD,double CbotD,double tempD ,double Gsalt);
    int selfco(int seqtop[lmax], int seqbot[lmax], int leff);
    double inittermprop37(int seqtop[lmax], int seqbot[lmax],int initright, int initleft, double Eprop[4][4][4][4], double Eterm[5][5][5][5], double EinitAT,double EinitGC, int leff, int scomp);
    void bot35(int seqbot[lmax],int seqbot35[lmax],int leff);
    G37=0.0;
    H[1]=0.0;
    /**
    printf("leff %d",leff);
    **/ 
    bot35(seqbot,seqbot35,leff);    
    initleft=1; 
    initright=leff-2;
    /**
    printf("main after pen %d %d", initleft,initright);
    **/
    if (initright<=initleft)
        {
        H[0]=0.0;
        Ssalt[0]=0.0;
        Gsalt[0]=0.0;
        Tm[0]=0.0;
        xD[0]=0.0;
        }
    else
        {
        saltlen=leff-2;
        /**
        printf("saltlen %f  ",saltlen);
        **/
        H[0]=inittermprop37(seqtop, seqbot35,initright, initleft, Hprop, Hterm, HinitAT,HinitGC, leff, 0);
        G37=inittermprop37(seqtop, seqbot35,initright, initleft,Gprop, Gterm,GinitAT,GinitGC, leff, scomp);
        /**
        printf("main after inittermprop37 Gtemp %f %f ",G37, H);
        **/
        S=((H[0]-G37)/310.15)*1000.0;
        Gtemp=hybtemp(H[0], S, tempD);
        /**
        printf("Gtemp %f  ",Gtemp);
        **/
        Gsalt[0]=saltcorrG( Gtemp,   CMoD, CMgD, saltlen,tempD);
        /**
        printf("Gsalt %f ",Gsalt[0]);
        **/
        Ssalt[0]=saltcorrS( H[0],  G37, CMoD, CMgD, saltlen);
        /**
        printf("Ssalt %f",Ssalt[0]);
        **/
        Tm[0]=tmcalc( CtopD, CbotD, Ssalt[0], H[0],scomp);
        /**
        printf("Tm %f",Tm[0]);
        **/
        xD[0]=FractionD(CtopD,CbotD,tempD,Gsalt[0]);
        /**
        printf("xD %f",xD[0]); 
        **/
        }
    /**     
    printf( "<font color=000066>");
    printf("&nbsp;&nbsp;&nbsp;%10.1f&nbsp;&nbsp;&nbsp;%10.1f&nbsp;&nbsp;&nbsp;%11.2f&nbsp;&nbsp;&nbsp;%10.1f&nbsp;&nbsp;&nbsp;%f\n",*H,*Ssalt,*Gsalt,*Tm,*xD);
    printf("</font><Br>");
    **/
    return 0;
    }
/*************************************************************/
/*Takes care of gaps and coaxial stacking*/
/*************************************************************/
int gapstack(int seqtop[lmax],int seqbot35[lmax], int *leffI,int lbotI ,int *errorI,int *coaxposI,int *gaposI, int *gaplenI)
    {
    int i,k,safareaI=0;
    int ncoaxI;
    *gaposI=0;
    *gaplenI=0;
    k=0;
    ncoaxI=0;
    *coaxposI=0;
    /**
    printf("<Br>gapstack entry\n<Br>");
    printf("errorI %d",*errorI);
    fflush(stdout);
    printf("<Br>\nleffI %d, lbotI %d", *leffI,lbotI);
    **/
    /*----------if more than one stacking: error*/
    for (i=0;i<=lbotI-1; i++)
        {
        if (seqbot35[i]==6)
            {
            ncoaxI=ncoaxI+1;
            }
        }
        
    for (i=0;i<=*leffI-1; i++)
        {
        if (seqtop[i]==6)
            {
            ncoaxI=ncoaxI+1;
            }
        }
    if (ncoaxI>1)
        {
        *errorI=1;
        /**
        printf(" ncoax %d",ncoaxI); 
        fflush(stdout);
        **/
        }
    else
        {
        /*----------Finds stacking in top strand and excise stacking symbol 6 and return -coaxposI*/
        for (i=0;i<=*leffI-1; i++)
            {
            if (seqtop[i]==6)
                {
                ncoaxI=1;   
                *coaxposI=-(i-1);   
                for (k=i;k<=*leffI;k++)
                    {
                    seqtop[k]=seqtop[k+1];
                    }
                *leffI=*leffI-1;
                }
            }
        /**
        printf("<Br>coaxposI top %d",*coaxposI); 
        printf("errorI %d",*errorI);
        fflush(stdout);
        **/
        /*----------*/
        /*----------Finds stacking in bottom strand and excise stacking symbol 6 and return +coaxposI*/
        /**
        printf("<Br>");
        for(i=0;i<=*leffI-1;i++)
            {
            printf(" %d ",seqtop[i]);
            
            }
        printf("<Br>");
        for(i=0;i<=lbotI-1;i++)
            {
            printf(" %d ",seqbot35[i]);
            
            }
        printf("<Br>");
        **/     
        for (i=0;i<=lbotI-1; i++)
            {
            if (seqbot35[i]==6)
                {
                
                *coaxposI=i-1;
            
                for (k=i;k<=lbotI;k++)
                    {
                    seqbot35[k]=seqbot35[k+1];
                    }
                lbotI=lbotI-1;
                }
            }
        /**     
        for(i=0;i<=*leffI-1;i++)
            {
            printf(" %d ",seqtop[i]);
            
            }
        printf("<Br>");
        for(i=0;i<=lbotI-1;i++)
            {
            printf(" %d ",seqbot35[i]);
            
            }   
        
        **/
        /**
        printf("<Br>coaxposI bot %d",*coaxposI);
        printf("errorI %d",*errorI); 
        fflush(stdout);
        **/
        /*----------*/
        }
    /*----------check that there is no mismatch at the stacking interface - errorI=3*/
    /*printf("coaxpos %d ncoax %d",*coaxposI,ncoaxI);*/
    fflush(stdout);
    if (ncoaxI!=0)
        {
        if ((abs(*coaxposI)<3)||(abs(*coaxposI)>(*leffI-4)))
            {
            *errorI=2;
            }   
        if ((seqbot35[abs(*coaxposI)]+seqtop[abs(*coaxposI)])!=3)
            {
            *errorI=3;
            printf("------%d %d",seqbot35[abs(*coaxposI)],seqtop[abs(*coaxposI)]);
            }
        if ((seqbot35[abs(*coaxposI)+1]+seqtop[abs(*coaxposI)+1])!=3)
            {
            *errorI=3;
            printf("------%d %d",seqbot35[abs(*coaxposI)],seqtop[abs(*coaxposI)]);
            }
        }
    /**
    printf("*errorI_1 %d",*errorI);
    fflush(stdout);
    **/
    /*----------check that the two strands have similar length - errorI=4*/
    
    if ((*leffI!=lbotI)&&(*errorI==0))  
        {
        *errorI=4;
        }
        /*printf("<Br>\nleffI %d, lbotI %d", *leffI,lbotI);*/
    /**
    printf("*errorI_4 %d",*errorI);
    fflush(stdout);
    **/
    /*----------Finds gaps in top strand*/
    /*Only allows one gap; gaps in top and bottom strand: errorI=5; more than one gap in top strand: errorI=6*/ 
    if (*errorI==0)
        {
        for (i=0;i<=*leffI-1; i++)
            {
            if (seqtop[i]==7)
                {
                if (*gaposI==0)
                    {       
                    for (k=0;k<=*leffI-1; k++)
                        {
                        if (seqbot35[k]==7)
                            {
                            *errorI=5;
                            }
                        else
                    
                            {
                            *gaposI=i-1;
                            *gaplenI=1;
                            }
                        }
                    }
                else
                    {
/**
                    printf("i %d",i);
**/
                    if (i==(*gaposI+*gaplenI+1))
                        {
                        *gaplenI=*gaplenI+1;
                        }
                    else                
                        {
                        *errorI=6;
                        }
                    }
                }
            }
        }
    /**
    printf("*errorI_5_6 %d",*errorI);
    fflush(stdout);
    **/
    /*Finds gaps in bottom strand*/
    /*----------Only allows one gap;  more than one gap in bottom strand: errorI=7*/    
    if ((*gaposI==0)&&(*errorI==0))
        {
        for (i=0;i<=*leffI-1; i++)
            {
            /*printf("b%db",seqbot35[i]);*/
            if (seqbot35[i]==7)
                {
                if (*gaposI==0)
                    {       
                    *gaposI=i-1;
                    *gaplenI=1;
                    /*printf("jkhgfksjh",*gaposI);*/
                    }
                else
                    {
                    if (i==(*gaposI+*gaplenI+1))
                        {
                        *gaplenI=*gaplenI+1;
                        }
                    else
                        {
                        *errorI=7;
                        }
                    }
                }
            }
        }
    /**
    printf("*errorI_7 %d",*errorI);
    fflush(stdout);
    **/
        /*printf("gaposI %d  gaplenI %d",*gaposI,*gaplenI);*/
    /*----------Determine safety area size as a function of gap length*/
    if (*errorI==0)
        {
        if (*gaplenI==1)
            {
            safareaI=4;
            }
        else if ((*gaplenI>1)&&(*gaplenI<=10))
            {
            safareaI=5;
            }
        else if (*gaplenI>10)
            {
            safareaI=6;
            }
        }
    /*----------Check no mismatch or duplex extremity within safareaI bp of gap; if any: errorI=8*/
    if ((*errorI==0)&&(*gaplenI!=0))
        {
        if ((*gaposI<safareaI-1)||((*gaposI+*gaplenI)>(*leffI-safareaI-1)))
            {
            *errorI=8;
/**
            printf("test0");
**/
            }
        if (*errorI==0)
            {
            for (i=*gaposI+*gaplenI+1;i<=*gaposI+*gaplenI+safareaI;i++)
                {
                if ((seqbot35[i]+seqtop[i])!=3)
                    {
/*
                    printf("test1");
*/
                    *errorI=8;
                    break;
                    }
                }
            }
        if (*errorI==0)
            {
            for (i=*gaposI;i>=*gaposI-safareaI+1;i--)
                {
                if ((seqbot35[i]+seqtop[i])!=3)
                    {
/*
                    printf("test2");
*/
                    *errorI=8;
                    break;
                    }
                }   
            }
        }
    /**
    printf("*errorI_8 %d",*errorI);
    fflush(stdout);
    **/
    /*----------Check no stacking within safareaI bp of gap;if any: errorI=9*/
    if ((*errorI==0)&&(ncoaxI!=0)&&(*gaplenI!=0))
        {       
        if ((abs(*coaxposI)>=*gaposI-safareaI+1)&&(abs(*coaxposI)<=*gaposI+*gaplenI+safareaI))
            {
            *errorI=9;
            }           
        }
/* 
        printf("<Br>xxx ");
*/

    /**
    for(i=0;i<=9;i++)
        {
        printf(" %d ",seqtop[i]);
            
        }
    printf("<Br>");
    for(i=0;i<=9;i++)
        {
        printf(" %d ",seqbot35[i]);     
        }
        printf("<Br>");
    **/
    /**
    printf("*errorI_9 %d",*errorI);
    fflush(stdout);
    **/     
    /*----------generate new top and bottom strand*/
    if (*errorI==0)
        {
        for (i=*gaposI+1;i<=*leffI;i++)
            {
            seqtop[i]=seqtop[i+*gaplenI];
            seqbot35[i]=seqbot35[i+*gaplenI];
            /*leffI=leffI-1;*/          
            }
        *leffI=*leffI-*gaplenI;
        }
        
        /**
        for(i=0;i<=*leffI-1;i++)
            {
            printf(" %d ",seqtop[i]);
            
            }
        printf("<Br>");
        for(i=0;i<=*leffI-1;i++)
            {
            printf(" %d ",seqbot35[i]);
            
            }
        **/
    /*----------decrement coaxial stacking position if gap before stacking interface*/
    /*and if there is a stacking*/
    if ((*gaposI<abs(*coaxposI))&&(ncoaxI!=0))
        {
        *coaxposI=*coaxposI/abs(*coaxposI)*(abs(*coaxposI)-*gaplenI);
        }
    
    /*----------*/
    /**
    printf("<Br>coaxposI %d",*coaxposI);
    printf("<Br>gaposI: %d gaplenI: %d", *gaposI,*gaplenI);
    printf("<Br>gapstack exit\n<Br>");
    fflush(stdout);
    **/     
    return 0;
    }
/*************************************************************/
/*Calculates gap contribution*/
/*************************************************************/ 
double gapcalcD(int gaposI,int gaplenI,double GgapD[30],double GATclosepenD,double Eprop[4][4][4][4],int seqtop[lmax], int seqbot35[lmax])
    {
    double EgapcorrD;
    EgapcorrD=0.0;
    /*----------if bulge of one add bulge increment*/
    if (gaplenI==1)
        {
        EgapcorrD=GgapD[0];
        }
    /*----------*/
    /*----------Limit gap length correction to 30 nts*/
    if (gaplenI>30)
        {
        gaplenI=30;
        }
    /*----------*/
    /*----------if bulge of 2-30, subtract NN add bulge increment and 0.5 kcal/mol penalty per AT closing*/
    if ((gaplenI>1)&&(gaplenI<=30))
        {
        /*----------add bulge increment*/
        /*gaplenI-1 because gaps in array starts at 0 for gap of 1*/
        EgapcorrD=GgapD[gaplenI-1];
        /*----------*/
        /*----------subtract NN increment*/
        EgapcorrD=EgapcorrD-Eprop[seqtop[gaposI]][seqtop[gaposI+1]][seqbot35[gaposI+1]][seqbot35[gaposI]];
        /*----------*/      
        /**
        printf("%d %d %d %d %d %d",GgapD[0],GgapD[1],GgapD[gaplenI-1],GgapD[10],GgapD[gaplenI-1],GgapD[20]);
        printf("gap %f %f gap", GgapD[gaplenI-1],Eprop[seqtop[gaposI]][seqtop[gaposI+1]][seqbot35[gaposI+1]][seqbot35[gaposI]]);
        printf("AT closing %f",GATclosepenD);
        **/
        /*----------add AT closing penalties*/
        if (seqtop[gaposI]==0 || seqtop[gaposI]==3)
            {
            EgapcorrD=EgapcorrD+GATclosepenD;
            }
        if (seqtop[gaposI+1]==0 || seqtop[gaposI+1]==3)
            {
            EgapcorrD=EgapcorrD+GATclosepenD;
            }
        /*----------*/
        }
    /**
    printf("EgapcorrD %f",EgapcorrD);
    **/
    return EgapcorrD;
    }
/*************************************************************/
/*changes bottom strand from 5'->3' to 3'->5'*/
/*************************************************************/
void bot35(int seqbot[lmax],int seqbot35[lmax],int leff)
    {
    int i;
    for (i=0; i<=leff-1; i++)
        {
        seqbot35[leff-i-1]=seqbot[i];
        }
    /**
    for (i=0; i<=9;i++)
        {                       
        printf("$%d$",seqbot35[i]);
        }
    **/
    }


/*************************************************************/
/*Tests for self complementarity - bottom strand explicit case*/
/*************************************************************/
/*seqtop sequence of top strand*/
/*seqbot sequence of bottom strand*/
/*leff effective length of strands*/
/*scomp is 0 if non-self-complementary and 1 if self-complementary*/
int selfco( int seqtopI[lmax], int seqbotI[lmax], int leff)
    {
    int i,scompI;
    scompI=1;
    /*printf("%d\n",scomp);*/
    for (i=0; i<=(leff-1); i++)
        {  
        if (seqtopI[i] != seqbotI[leff-i-1])
            {
            scompI=0;
            break;
            }   
        }
    return scompI;
    }
/*************************************************************/
/*Tests for self complementarity - bottom strand implicit case*/
/*************************************************************/
int selfcoacI( int seqtopI[lmax], int leffI)
    {int i,scompI;
    scompI=1;
    /*printf("%d\n",scomp);*/
    for (i=0; i<=(leffI-1); i++)
        {  
        if ((seqtopI[i]+seqtopI[leffI-i-1])!=3)
            {
            scompI=0;
            break;
            }   
        }
    return scompI;
    }
/*************************************************************/
/*Treats left penultimate mismatches and terminal elements*/
/*************************************************************/
/*initleft start of propagation on left side*/
int penleft(int seqtop[lmax], int seqbot[lmax], double Gprop[4][4][4][4], double Gterm[5][5][5][5],double GinitAT,double GinitGC, int leff)
    {
    int i, initleft;
    double Gopen,Gclose;
    /*----------determines initleft*/
    DB_DSET {printf(">> penleft\n"); fflush(stdout);    }
    i=-1;
    initleft=-1;
    do
        {   
        i=i+1; 
        /*printf("i=%d %d %d\n",i,seqtop[i],seqbot[i]);*/
        if (seqtop[i]+seqbot[i]==3)
            {
             if (seqtop[i+1]+seqbot[i+1]==3)
                {
                initleft=i;
                }
            }
        }
    while ((initleft==-1)&&(i<(leff-1)));
    /*----------*/
    /*----------if there is an initleft try to calculate pen */
    if ((seqtop[i-1]+seqbot[i-1]!=3) && (seqtop[i-2]+seqbot[i-2]==3)&&(i>2)&&(initleft!=-1))
        {
        Gopen=Gterm[seqtop[i-1]][seqtop[i]][seqbot[i]][seqbot[i-1]];
        Gclose=Gprop[seqtop[i-1]][seqtop[i]][seqbot[i]][seqbot[i-1]]+Gprop[seqtop[i-2]][seqtop[i-1]][seqbot[i-1]][seqbot[i-2]];
        if (seqtop[i-2]==1 || seqtop [i-2]==4)
            {
            Gclose=Gclose+GinitAT;
            }
        else
            {
            Gclose=Gclose+GinitGC;
            }
        if (i+2!=leff-1)
            {
            Gclose=Gclose+Gterm[seqtop[i-3]][seqtop[i-2]][seqbot[i-2]][seqbot[i-3]];
            }
        if (Gopen>Gclose)
            {
            initleft=i;
            }
        else
            {
            initleft=i-2;
            }
        }
    
    DB_DSET {printf("<< penleft %d\n",initleft); fflush(stdout);    }
    return initleft;
    }
/*************************************************************/
/*Treats right penultimate mismatches and terminal elements*/
/*************************************************************/
/*initleft start of propagation on left side*/
/*initright start of propagation on right side*/
/*determines initright*/
int penright(int seqtop[lmax], int seqbot[lmax], double Gprop[4][4][4][4], double Gterm[5][5][5][5],double GinitAT,double GinitGC,int leff)
    {
    int i,initright;
    double Gopen,Gclose;

    DB_DSET {printf(">> penright\n"); fflush(stdout);   }
    initright=-1;
    i=leff;
    /*----------Determines initright*/
    do
        {   
        i=i-1; 
        /*printf("i%d*",i);*/
        if (seqtop[i]+seqbot[i]==3)
            {   
            if (seqtop[i-1]+seqbot[i-1]==3)
                {
                initright=i;
                /*printf("initright%d ",initright);*/
                DB_DSET {printf("-2-initright=%d\n",initright); fflush(stdout); }
                }
            }
        }
    while ((initright==-1)&&(i>=1));
    /*----------*/
    /*----------If there is an initright calculattes pen*/
    DB_DSET {printf("-3-initright=%d\n",initright); fflush(stdout); }
    DB_DSET {printf("+ pastwhile i=%d\n",i); fflush(stdout);    }
    if ((seqtop[i+1]+seqbot[i+1]!=3) && (seqtop[i+2]+seqbot[i+2]==3)&&(i!=-1)&&(i<(leff-2))&&(initright!=-1))
        {
            
        fflush(stdout);   
        DB_DSET {printf("+ 3&3 i=%d\n",i); fflush(stdout);  }
        Gopen=Gterm[seqtop[i]][seqtop[i+1]][seqbot[i+1]][seqbot[i]];
        Gclose=Gprop[seqtop[i]][seqtop[i+1]][seqbot[i+1]][seqbot[i]]+Gprop[seqtop[i+1]][seqtop[i+2]][seqbot[i+2]][seqbot[i+1]];
        if (seqtop[i+2]==1 || seqtop [i+2]==4)
            {
            Gclose=Gclose+GinitAT;
            }
        else
            {
            Gclose=Gclose+GinitGC;
            }
        if (i+2!=leff-1)
            {
            Gclose=Gclose+Gterm[seqtop[i+2]][seqtop[i+3]][seqbot[i+3]][seqbot[i+2]];
            }
            DB_DSET {printf("+ half %d\n",i); fflush(stdout);   }
            DB_DSET {printf("-5-initright=%d\n",initright); fflush(stdout); }
            
        if (Gopen>Gclose)
            {
            initright=i;
            }
        else
            {
            initright=i+2;
            }
            
        }
    /*----------*/
    DB_DSET {printf("<< penright %d\n",initright); fflush(stdout);  }
    return initright;
    }
/*************************************************************/
/*Determines H and G37 and S for initiation,propagation and*/
/*terminal elements */
/*************************************************************/
/*initright initiation contribution on right side*/
/*initleft initiation contribution on left side*/
/*initAT initiation contribution for termminal AT*/
/*initGC initiation contribution for terminal GC*/
/*termleft left terminal element contribution*/
/*termright right terminal element contribution*/ 
double inittermprop37(int seqtop[lmax], int seqbot[lmax],int initright, int initleft, double Eprop[4][4][4][4], double Eterm[5][5][5][5], double EinitAT,double EinitGC, int leff, int scomp)
    {
    int i;
    double Einitright=0, Einitleft=0, Etermright=0,Etermleft=0,E;

    DB_DSET {printf(">> inittermprop37 scomp=%d\n",scomp); fflush(stdout);  }
    E=0.00;
    /*----------sum propagation increments*/
    DB_DSET {printf("initright=%d, initleft=%d",initright, initleft); fflush(stdout);}
    for (i=initleft;i<initright;i++)
        {
        E=E+Eprop[seqtop[i]][seqtop[i+1]][seqbot[i+1]][seqbot[i]];
        /*printf("E  after prop %f\n",E);*/
        }
    DB_DSET {printf("test"); fflush(stdout);    }   
    /*----------*/
    /*----------add initiation increments*/
    if (seqtop[initright]==0 || seqtop [initright]==3)
        {
        Einitright=EinitAT;
        }
    else
        {
        Einitright=EinitGC;
        } 
    if (seqtop[initleft]==0 || seqtop [initleft]==3)
        {
        Einitleft=EinitAT;
        }
    else
        {
        Einitleft=EinitGC;
        } 
    DB_DSET {printf("Einitright,Einitleft %f%f\n",Einitright,Einitleft); fflush(stdout);    }
    /*printf("Einitright,Einitleft %f%f\n",Einitright,Einitleft);*/
    /*----------*/
    /*----------add terminal-element increments*/
    if (initleft !=0 )
        {
        Etermleft=Eterm[seqtop[initleft-1]][seqtop[initleft]][seqbot[initleft]][seqbot[initleft-1]];
        /*printf("- %d %d %d %d-",seqtop[initleft-1],seqtop[initleft],seqbot[initleft],seqbot[initleft-1]);*/
        }
    if (initright != leff-1)
        {
        Etermright=Eterm[seqtop[initright]][seqtop[initright+1]][seqbot[initright+1]][seqbot[initright]];
        /*printf("- %d %d %d %d-",seqtop[initright],seqtop[initright+1],seqbot[initright+1],seqbot[initright]);*/
        }
    DB_DSET {printf("Etermright,Etermleft %f%f\n",Etermright,Etermleft); fflush(stdout);    }   
    /*printf("Etermright,Etermleft %f%f\n",Etermright,Etermleft);*/
    /*----------*/
    /*calculation of G37, H and S*/
    E=E+Einitleft+Einitright+Etermleft+Etermright+scomp*0.4;
    /*printf("inittermprop37 completed %f\n",E);*/
    DB_DSET {printf("<< inittermprop37 %f\n",E); fflush(stdout);    }
    return E;
    }
/*--------------------------------------------------------------*/
double inittermprop37acD(int seqtop[lmax], int initright, int initleft, double Eprop[4][4][4][4], double Eterm[5][5][5][5], double EinitAT,double EinitGC, int leff, int scomp)
    {
    int i;
    double Einitright=0, Einitleft=0,E;
    E=0.00;
    /*
    printf(" %f %f ",EinitAT, EinitGC);
    */
    /*propagation*/
    /*printf("initright %d initleft %d",initright, initleft);*/
    for (i=initleft;i<initright;i++)
        {
        E=E+Eprop[seqtop[i]][seqtop[i+1]][3-seqtop[i+1]][3-seqtop[i]];
        /*printf("E  after prop %f\n",E);*/
        }
    /*initiation*/
    if (seqtop[initright]==0 || seqtop [initright]==3)
        {
        Einitright=EinitAT;
        }
    else
        {
        Einitright=EinitGC;
        } 
    if (seqtop[initleft]==0 || seqtop [initleft]==3)
        {
        Einitleft=EinitAT;
        }
    else
        {
        Einitleft=EinitGC;
        } 
        /*printf("Einitright,Einitleft %f%f\n",Einitright,Einitleft);*/
    /*terminal elements*/
    /*calculation of G37, H and S*/
    E=E+Einitleft+Einitright+scomp*0.4;
    /*printf("inittermprop37 completed %f\n",E);*/
    return E;
    }
/*************************************************************/
/*Calculates effective length for salt correction*/
/*************************************************************/
/*assumes no terminal phosphates*/
/*to use terminal phosphates use initright-initleft+1 instead of iniright-initleft*/
/*stars count of dangling ends involved*/
/*each DE accounts for 0.5*/
/*each TM accounts for 1.0*/ 
int saltlength(int seqtop[lmax],int seqbot[lmax],int initright, int initleft, int leff)
    {
    int stars;
    double saltlen;

    DB_DSET {printf(">> saltlength\n"); fflush(stdout); }
    saltlen=0;
    stars=0;
    /**
    printf ("<Br>function saltlength: initright %d, initleft %d leff %d <Br>",initright, initleft,leff);
    **/
    /*----------if no terminal element*/
    if ((initleft==0) && (initright==leff-1))
        { 
        saltlen=(initright-initleft);
        }
    /*----------*/
    /*----------if 1 terminal element at rigtht extremity*/
    else if ((initleft ==0) && (initright != leff-1)) 
        {
        /*----------if DE then 0.5 else TM then 1.0*/
        if ( (seqbot[initright+1]==4)||(seqtop[initright+1]==4))
            {
            saltlen=(initright-initleft)+0.5;
            }
        else
            {
            saltlen=(initright-initleft)+1;
            } 
        /*----------*/   
        }
    /*----------*/
    /*----------if 1 terminal element at left extremity*/
    else if((initleft!=0) && (initright==leff-1)) 
        {
        /*----------if DE at left extremity then 0.5 else TM then1.0*/
        if ((seqbot[initleft-1]==4)||(seqtop[initleft-1]==4))
            { 
            saltlen=(initright-initleft)+0.5;
            }
        else
            {
                    saltlen=(initright-initleft)+1.0;
                    }
                /*----------*/    
        }
    /*----------if 2 terminal elements*/  
    else if((initleft!=0) && (initright!=leff-1)) 
        {
        if (seqbot[initleft-1]==4)
            stars=1;       
        if (seqbot[initright+1]==4 )
            stars=stars+1;       
        if (seqtop[initleft-1]==4)
            stars=stars+1;      
        if (seqtop[initright+1]==4)
            stars=stars+1;     
        if (stars==0) 
            saltlen=(initright-initleft)+2.0;
        else if (stars==1) 
            saltlen=(initright-initleft)+0.5+1.0; 
        else if (stars==2) 
            saltlen=(initright-initleft)+1.0;      
        }
    /*----------*/
    /**
    printf("saltlength completed %f\n",saltlen);
    **/
    DB_DSET {printf("<< saltlength %f\n",saltlen); fflush(stdout);  }
    return saltlen;
    }
/*************************************************************/
/*Corrects Gtemp for salt*/
/*************************************************************/
/*Gsalt, Ssalt free energy and entropy corrected for salt*/
/*saltlen effective length for salt correction*/
/*CMo and CMg concentrations of monovalent cations and Mg2+*/
double saltcorrG( double Gtemp,  double CMo, double CMg, double saltlen,double tempD)
    {
    double Gsalt;
    Gsalt=Gtemp-0.114*(saltlen)*log(CMo+3.3*sqrt(CMg))*((tempD+273.15)/310.15);
    /**
    printf("<Br>function saltcorrG: CMo %f, CMg %f, saltlen %f, Gtemp %f, Gsalt %f<Br>", CMo, CMg, saltlen, Gtemp, Gsalt);
    **/
    return Gsalt;
    }
/*************************************************************/
/*Corrects S for salt*/
/*************************************************************/
/*Gsalt, Ssalt free energy and entropy corrected for salt*/
/*saltlen effective length for salt correction*/
/*CMo and CMg concentrations of monovalent cations and Mg2+*/
double saltcorrS(double H, double G37, double CMo, double CMg, double saltlen)
    {
    double Ssalt;
    Ssalt=((H-(G37-0.114*(saltlen)*log(CMo+3.3*sqrt(CMg))))/310.15)*1000.0;
    return Ssalt;
    }
/*************************************************************/
/*Corrects G for non cannonical hybridization temp*/
/*************************************************************/
/*G37 total free energy at 37*/
/*Gtemp total free energy at temp*/
/*H total enthalpy*/
/*S total entropy*/
double hybtemp(double H, double S, double temp)
    {
    double Gtemp;
    Gtemp=(H*1000.0-(temp+273.15)*S)/1000.0;
    /*printf("hybtemp completed %f\n",Gtemp);*/
    return Gtemp;
    }
/*************************************************************/
/*Calculates Tm*/
/*************************************************************/
/*Tm melting temperature*/
/*Ctop, Cbot concentrations of top and bottom strands*/
/* salt entropy corrected for salt*/
/*  perfect gas constant in cal/(mol*k)*/
double tmcalc(double Ctop,double Cbot,double Ssalt,double H,  int scomp)
    {
    const double R=1.98722;
    double Tm=0.0;
    /*----------if non self complementary*/
    if (scomp !=1)
        {   
        if (Ctop >= Cbot)
            {
            Tm=(H*1000.0F/(Ssalt+R*log(Ctop-0.5*Cbot)))-273.15;
            /*printf("Tm %f\n",*Tm)*/;
            }
        else
            {
            Tm=(H*1000/(Ssalt+R*log(Cbot-0.5*Ctop)))-273.15; 
            }
        }
    /*----------*/
    /*----------if self complementary*/
    else
        {
        Tm=(H*1000/(Ssalt+R*log(Ctop)))-273.15;
        /*printf("scomp TM calc");*/
        }
    /*----------*/
        /*printf("tmcalc completed %f\n",*Tm);*/
    return Tm;
    }
/*************************************************************/
/*Calculates coaxial stacking contribution*/
/*************************************************************/
/*Hcoax enthalpy for coaxial stacking*/
/*Gcoax free energy for coaxial stacking*/
/*Coaxpos positions of stacking interfaces*/
double coaxcalc(double EcoaxD[4][4][4][4],int seqtopI[lmax] , int seqbot35I[lmax],int *coaxposI)
    {   
    double EcoaxcorD;
    EcoaxcorD=0.0;
    /*----------if coaxial stacking on bottom strand*/
    if (*coaxposI>0)
        {
        /**
        printf("----------=+%d %d %d %d =+",seqtopI[*coaxposI-1],seqtopI[*coaxposI],seqbot35I[*coaxposI],seqbot35I[*coaxposI-1]);
        **/
        /*printf(" %d",*coaxposI);*/
        /**
        printf("=+%d %d %d %d =+",seqtopI[*coaxposI],seqtopI[*coaxposI+1],seqbot35I[*coaxposI+1],seqbot35I[*coaxposI]);
        **/
        EcoaxcorD=EcoaxcorD+EcoaxD[seqtopI[*coaxposI]][seqtopI[*coaxposI+1]][seqbot35I[*coaxposI+1]][seqbot35I[*coaxposI]];
        fflush(stdout);
        /**
        printf("Ecoax after coax %f",EcoaxcorD);
        **/ 
        }
    /*----------*/
    /*----------if coaxial stacking on top strand*/
    else if (*coaxposI<0)
        {
        *coaxposI=-*coaxposI;
        /*printf("jhgghg <0<0");*/
        EcoaxcorD=EcoaxcorD+EcoaxD[seqbot35I[*coaxposI+1]][seqbot35I[*coaxposI]][seqtopI[*coaxposI]][seqtopI[*coaxposI+1]];
        /*printf("=+%d %d %d %d =+",seqbot35I[*coaxposI+1],seqbot35I[*coaxposI],seqtopI[*coaxposI],seqtopI[*coaxposI+1]);*/
        *coaxposI=-*coaxposI;
        }
    /*-----------*/
    /**
    printf("Ecoax end of coax %f coaxposI %d",EcoaxcorD,*coaxposI);
    **/         
    return EcoaxcorD;
    }
/*************************************************************/
/*Calculates fraction duplex*/
/*************************************************************/
double FractionD(double CtopD,double CbotD,double tempD ,double Gsalt)
    {   
    double aD,bD,cD,qD,KABD,xD,x1D,climD,nsD,difD;

    if (CtopD>CbotD)
        {
        climD=CbotD;
        }   
    else
        {
        climD=CtopD;
        }
    /*----------Calculates K*/
    /*printf("Gsalt **%f--",Gsalt);*/
    KABD=exp(-(Gsalt*1000.0)/(1.98722*(tempD+273.15)));
    /*Calculates coefficients of quadratic equation a*AB^2+b*AB+c=0*/
    aD=KABD;
    bD=-KABD*CtopD-KABD*CbotD-1;
    cD=KABD*CtopD*CbotD;
    /*----------*/
    /*----------Solve quadratic equation -numerical solution from numerical recipes in fortran p78*/
    if (bD<0)
        {
        nsD=-1.00;
        }
        else
        {
        nsD=1.00;
        }
    difD = bD*bD-4.0*aD*cD;
    /*----------in case difD<0 to avoid NAN - no theoretical justification*/ 
    if (difD<0.0)
        {
        difD=sqrt(difD*difD);
        }
    /*----------*/
    qD = -0.5*(bD+nsD*sqrt(difD));
    x1D=cD/qD;
    if ((x1D>=0)&&(x1D<=climD))
        {
        xD=cD/qD;
        }
    else
        {
        xD=qD/aD;
        }
    /*----------calculates fraction bound*/
        xD=xD/climD*100.0;              
    /*printf("tempD %f K %f, a %f, b %f, c %f, q %f, X %f\n",tempD,KABD,aD,bD,cD,qD,xD);*/
    /*printf("%f ,%f ,%f ,%f,qD: %f ,xD: %f\n",eD,aD,bD,cD,qD,xD);*/
    /*----------*/

    return xD;
    }
