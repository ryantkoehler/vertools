/*
* dset_tm.c
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

/*---------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include "prim.h"
# include "dna.h"
# include "tm_pars.h"
# include "dset.h"
# include "dset_tm.h"

#define DB_DSET if(DB[123])

/*************************************************************/
/*Main interface between Ryan environment and DSET trunk*/
/*************************************************************
*	seqtopPC & seqbotPC are null-terminated strings that should be of the
*	same length (not counting / chars). Depending on mode, seqbotPC may
*	be null.
*/
int SeqDsetEnergyI(TM_PARS *tmPO, char *seqtopPC, char *seqbotPC, int verboseI,
	int modeI, DOUB *tPD, DOUB *gPD, DOUB *hPD, DOUB *sPD, DOUB *xPD)
{
	int k,l,m,n,seqtopI[lmax],seqbotI[lmax],leffI,lbotI,errorI,ok;
	double CtopD,CbotD,CMoD,CMgD,tempD,Cbot2D;
	double Gterm[5][5][5][5],Hterm[5][5][5][5];
	double Gprop[4][4][4][4],Hprop[4][4][4][4];
	double Gcoax[4][4][4][4],Hcoax[4][4][4][4];
	double GgapD[30];
	double H[4],Ssalt[4],Gsalt[5], Tm[5], xD[5];
	double HinitAT,GinitAT,HinitGC,GinitGC,GATclosepenD;

	DB_DSET DB_PrI(">> SeqDsetEnergyI mode=%d\n",modeI);
	VALIDATE(tmPO,TM_PARS_ID);
	if(!ValidDsetModeI(seqtopPC,seqbotPC,modeI)) {
		DB_DSET DB_PrI("<< SeqDsetEnergyI mode not valid FALSE\n");
		return(FALSE);
	}
	DB_DSET {
		DumpTmPars(tmPO, "+ ", TRUE, NULL);
	}
	errorI=0;
	ok = TRUE;

	/*----------loads thermodynamic parameters*/
        /*paramloaddna(Gprop,Gterm,Hprop, Hterm,&GinitAT,&GinitGC, &HinitAT,&HinitGC,Gcoax,Hcoax,GgapD,&GATclosepenD);*/
	/*----------*/
	/*----------initialization of local hybridization variables from structure*/
	CtopD=tmPO->conc;
	CbotD=tmPO->conc2;
	Cbot2D=tmPO->conc3;
	CMoD=tmPO->salt;
	CMgD=tmPO->mg;
	tempD=tmPO->tp;
	/*----------*/
	/*----------initialization of local thermodynamic variables from structure*/
	HinitAT=tmPO->HinitAT;
	GinitAT=tmPO->GinitAT;
	HinitGC=tmPO->HinitGC;
	GinitGC=tmPO->GinitGC;
	GATclosepenD=tmPO->GATclosepen;
	for (k=0;k<5;k++)
	{
		for (l=0;l<5;l++)
		{
			for (m=0;m<5;m++)
			{
				for (n=0;n<5;n++)
				{
        			Gterm[k][l][m][n]=tmPO->Gterm[k][l][m][n];
        			Hterm[k][l][m][n]=tmPO->Hterm[k][l][m][n];
        		}
			}
		}
	}
    for (k=0;k<4;k++)
	{
		for (l=0;l<4;l++)
		{
			for (m=0;m<4;m++)
			{
				for (n=0;n<4;n++)
				{
       				Gprop[k][l][m][n]=tmPO->Gprop[k][l][m][n];
       				Hprop[k][l][m][n]=tmPO->Hprop[k][l][m][n];
       			}
			}
		}
	}
    for (k=0;k<4;k++)
	{
		for (l=0;l<4;l++)
		{
			for (m=0;m<4;m++)
			{
				for (n=0;n<4;n++)
				{
       				Gcoax[k][l][m][n]=tmPO->Gcoax[k][l][m][n];
       				Hcoax[k][l][m][n]=tmPO->Hcoax[k][l][m][n];
				}
			}
		}
	}
    for (k=0;k<30;k++)	
	{
		GgapD[k]=tmPO->Ggap[k];
	}
	/*----------*/
	k=0;
	/*----------Case of 1 strand, the other is implicit*/
	if (modeI==DSET_IMP)
	{
		leffI = SeqToNumsI(seqtopPC,seqtopI,strlen(seqtopPC),modeI,verboseI);
		if(IS_BOG(leffI))
		{
			errorI = 10;
		}
		else
		{
			thermoacI(CMoD,CMgD,tempD,CtopD,CbotD,seqtopI,leffI,Gprop,
				Gterm,Hprop,Hterm, HinitAT,HinitGC, GinitAT,GinitGC,H,
				Ssalt,Gsalt,Tm,xD);
		}
	}
	/*----------case of 2 explicit strands*/
	else if ((modeI==DSET_EXPCOAX)||(modeI==DSET_EXPNOCOAX))
	{	
		leffI = SeqToNumsI(seqtopPC,seqtopI,strlen(seqtopPC),modeI,verboseI);
		lbotI = SeqToNumsI(seqbotPC,seqbotI,strlen(seqbotPC),modeI,verboseI);
		if( IS_BOG(leffI) || IS_BOG(lbotI) )
		{
			errorI = 10;
		}
		else
		{
			DB_DSET 
				DB_PrI("+ topD=%e botD=%e bot2D=%e tempD=%e cMoD=%e cMgD=%e\n",
					CtopD,CbotD,Cbot2D, tempD, CMoD,CMgD);

			ok = thermocI(CMoD,CMgD,tempD,CtopD,CbotD,Cbot2D,seqtopI,seqbotI,
				leffI,Gprop,Gterm,Hprop, Hterm,HinitAT,HinitGC, GinitAT,
				GinitGC,Gcoax,Hcoax,lbotI,&errorI,GgapD,GATclosepenD,H,
				Ssalt,Gsalt,Tm,xD);
		}
	}
	/*----------*/
	/*----------case of a primer*/		
	else if (modeI==DSET_PRIM)
	{
		/*non a,A,C,C,t,T,u,U,g,G characters are edited*/
		/* A=0,C=1,G=2,T=U=3*/
		leffI = SeqToNumsI(seqtopPC,seqtopI,strlen(seqtopPC),modeI,verboseI);
		if(IS_BOG(leffI))
		{
			errorI = 10;
		}
		else
		{
			/*----------Determines length and sets bottom extremities to DE*/
			seqbotI[0]=4;
			seqbotI[leffI-1]=4;
			/*----------*/
			/*----------Constructs bottom strand complementary to top strand 3' -> 5'*/
			for (l=1;l<=leffI-2;l++)
			{
				seqbotI[l]=3-seqtopI[leffI-l-1];
				/*rintf("seqbot %d",seqbotI[l]);*/	
			}
			/*----------*/
			/*----------executes function corresponding to mode*/
			thermocpI(CMoD,CMgD,tempD,CtopD,CbotD,seqtopI,seqbotI,leffI,
				Gprop,Gterm,Hprop,Hterm, HinitAT,HinitGC, GinitAT,GinitGC,
				&errorI,H,Ssalt,Gsalt,Tm,xD);
			/*----------*/
		}			
	}					 	
	/*----------*/
	/*----------Assigns the results to data structure*/
	/*----------print error messages if verbose on*/
	if(!ok)
	{
		if(verboseI)
			printf("# Error: Less than 2 consecutive base pairs.\n");
		return(FALSE);
	}
	switch(errorI)
	{
		case 1:
			if(verboseI)
				printf("# Error: more than one coaxial stacking interface.\n");
			return(FALSE);
		case 2:
			if(verboseI)
				printf("# Error: coaxial stacking interface too close to duplex extremities.\n");
			return(FALSE);
		case 3:
			if(verboseI)
				printf("# Error: mismatch at stacking interface.\n");
			return(FALSE);
		case 4:
			if(verboseI)
				printf("# Error: strands have different lengths.\n");
			return(FALSE);
		case 5:
			if(verboseI)
				printf("# Error: gaps in top and bottom strands.\n");
			return(FALSE);
		case 6:
			if(verboseI)
				printf("# Error: more than one gap in top strand.\n");
			return(FALSE);
		case 7:
			if(verboseI)
				printf("# Error: more than one gap in bottom strand.\n");
			return(FALSE);
		case 8:
			if(verboseI)
				printf("# Error: mismatch or strand end too close to gap.\n");
			return(FALSE);
		case 9:
			if(verboseI)
				printf("# Error: stacking interface too close to gap.\n");
			return(FALSE);
	}
	
	switch(modeI)
		{
		case DSET_EXPCOAX:
			if(gPD)
			{
/**
 DumpArray(Gsalt,5,NULL,IS_DOUB,NULL);
**/
				gPD[DSET_EV_5P]=Gsalt[0];
				gPD[DSET_EV_5P_ST]=Gsalt[1];
				gPD[DSET_EV_3P]=Gsalt[2];
				gPD[DSET_EV_3P_ST]=Gsalt[3];
				gPD[DSET_EV_FULL]=Gsalt[4];
			}
			if(hPD)
			{
				hPD[DSET_EV_5P]=H[0];
				hPD[DSET_EV_5P_ST]=H[1];
				hPD[DSET_EV_3P]=H[2];
				hPD[DSET_EV_3P_ST]=H[3];
				hPD[DSET_EV_FULL]=-TOO_BIG_R;
			}
			if(xPD)
			{
				xPD[DSET_EV_5P]=xD[0];
				xPD[DSET_EV_5P_ST]=xD[1];
				xPD[DSET_EV_3P]=xD[2];
				xPD[DSET_EV_3P_ST]=xD[3];
				xPD[DSET_EV_FULL]=xD[4];
			}
			if(sPD)
			{
				sPD[DSET_EV_5P]=Ssalt[0];
				sPD[DSET_EV_5P_ST]=Ssalt[1];
				sPD[DSET_EV_3P]=Ssalt[2];
				sPD[DSET_EV_3P_ST]=Ssalt[3];
				sPD[DSET_EV_FULL]=-TOO_BIG_R;
			}
			if(tPD)
			{
 /*
    DumpArray(Tm,5,NULL,IS_DOUB,NULL);
*/
    DumpArray(Tm, IS_DOUB, 0,5, NULL, NULL);
				tPD[DSET_EV_5P]=Tm[0];
				tPD[DSET_EV_5P_ST]=Tm[1];
				tPD[DSET_EV_3P]=Tm[2];
				tPD[DSET_EV_3P_ST]=Tm[3];
				tPD[DSET_EV_FULL]=Tm[4];
			}
			break;
		case DSET_EXPNOCOAX:
		case DSET_IMP:
		case DSET_PRIM:
			if(gPD)
			{
				gPD[DSET_EV_FULL]=Gsalt[0];
			}
			if(hPD)
			{
				hPD[DSET_EV_FULL]=H[0];
			}
			if(xPD)
			{
				xPD[DSET_EV_FULL]=xD[0];
			}
			if(sPD)
			{
				sPD[DSET_EV_FULL]=Ssalt[0];
			}
			if(tPD)
			{
				tPD[DSET_EV_FULL]=Tm[0];
			}
			break;
		default:
			printf("Bad mode = %d\n",modeI);
			ERR("SeqDsetEnergyI","Bogus mode given");
			return(FALSE);
	}
	/*----------*/
	DB_DSET DB_PrI("<< SeqDsetEnergyI TRUE\n");
	return(TRUE);	
}
/*************************************************************
*Checks if mode is in agreement with data passed to TmutilDsetI
*	seqbotPC may be NULL
*************************************************************/	
int ValidDsetModeI(char *seqtopPC,char *seqbotPC,int modeI)
{
	DB_DSET 
	{
		DB_PrI(">> ValidDsetModeI mode=%d\n",modeI);
		DB_PrI("+ top=%p bot=%p\n",seqtopPC,seqbotPC);
	}
	if(seqtopPC==NULL)
	{
		DB_DSET DB_PrI("<< ValidDsetModeI top==NULL FALSE\n");
		return(FALSE);
	}
	DB_DSET DB_PrI("+ top %d |%s|\n",strlen(seqtopPC),seqtopPC);
	switch(modeI)
	{
		case DSET_PRIM:
		case DSET_IMP:
			if(seqbotPC !=NULL)
			{
				DB_DSET DB_PrI("<< ValidDsetModeI bot!=NULL FALSE\n");
				return(FALSE);
			}
			break;
		case DSET_EXPNOCOAX:
			if(seqbotPC==NULL)
			{
				DB_DSET DB_PrI("<< ValidDsetModeI bot==NULL FALSE\n");
				return(FALSE);
			}
			DB_DSET DB_PrI("+ bot %d |%s|\n",strlen(seqbotPC),seqbotPC);
			if(strstr(seqbotPC,"/")!=NULL)
			{
				DB_DSET DB_PrI("<< ValidDsetModeI bot with / FALSE\n");
				return(FALSE);
			}
			if(strstr(seqtopPC,"/")!=NULL)
			{
				DB_DSET DB_PrI("<< ValidDsetModeI top with / FALSE\n");
				return(FALSE);
			}
			break;
		case DSET_EXPCOAX:
			if(seqbotPC==NULL)
			{
				DB_DSET DB_PrI("<< ValidDsetModeI bot==NULL FALSE\n");
				return(FALSE);
			}
			DB_DSET DB_PrI("+ bot %d |%s|\n",strlen(seqbotPC),seqbotPC);
			if(strstr(seqbotPC,"/")!=NULL)
				break;
			if(strstr(seqtopPC,"/")!=NULL)
				break;
			DB_DSET DB_PrI("<< ValidDsetModeI no / found FALSE\n");
			return(FALSE);
			break;
		default:
			DB_DSET DB_PrI("<< ValidDsetModeI bad mode FALSE\n");
			return(FALSE);
	}
	DB_DSET DB_PrI("<< ValidDsetModeI TRUE\n");
	return (TRUE);
}
/*************************************************************/
/*Returns the number of thermodynamic results to read depending on mode*/
/*************************************************************/	
int DsetNumValsI(int modeI)
	{
	int numberI;
	numberI=1;
	if (modeI==DSET_EXPCOAX)
		{
		numberI=5;	
		}
	return numberI;
	}
/*************************************************************
*   4/3/14 RTK; Update to catch return of fscanf ... No testing
*       Not sure this is even ever used? 
*       Nevertheless, try to keep original (hideous) formatting
*/
/*************************************************************/
/*Loads thermodynamic params into data struct from file*/
/*************************************************************/	
/*Load thermodynamic parameters and*/ 
/*hybridization parameters into tmPO*/
/*If problem, return FALSE, else TRUE*/

int LoadTmCompleteParsI(TM_PARS *tmPO,FILE *fPF)
	{
    int fok;
	int i,k,l,m,n;
	double Gtemp, Htemp,Gttemp, Httemp;
	double a1=1.001;
	char g[15],bufS[DEF_BS];

	DB_DSET DB_PrI(">> LoadTmCompleteParsI %p\n",fPF);
	/***
	*	Check for magic keyword
	*/
	n = fok = 0;
	while(fgets(bufS,LINEGRAB,fPF))
	{
		if(COM_LINE(bufS)) {
			continue;
        }
		Upperize(bufS);
		if(EQSTRING(bufS,"DSET_PARS",9)) {
			n++;
			break;
		}
	}
	if(!n) {
		PROBLINE;
		printf("File is missing magic keyword\n");
		DB_DSET DB_PrI("<< LoadTmCompleteParsI FALSE\n");
		return(FALSE);
	}
	/*parameter array initialization*/
	/*printf("test");*/
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
					for (n=0;n<=3;n++)										
					{
					tmPO->Gprop[k][l][m][n]=a1;
					tmPO->Hprop[k][l][m][n]=a1;
					}
	for (k=0; k<=4; k++)
		for (l=0; l<=4; l++)
			for (m=0;m<=4;m++)
				for (n=0;n<=4;n++)                                             
				  {
				  tmPO->Gterm[k][l][m][n]=a1;
				  tmPO->Hterm[k][l][m][n]=a1;
				  } 
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
					for (n=0;n<=3;n++)										
					{
					tmPO->Gcoax[k][l][m][n]=a1;
					tmPO->Hcoax[k][l][m][n]=a1;
					}
	/*printf ("testrtrt%.2f %.2f\n",Gprop[0][0][3][3],Gterm[3][4][3][0]);*/
	/*parameter array fill out by reading parameter file*/

	for (i=1; i<=16; ++i)
		{fok = fscanf(fPF,"%s%d%d%d%d%lf%lf",g,&k,&l,&m,&n,&Gtemp,&Htemp);
		 tmPO->Gprop[k][l][m][n]=Gtemp;
   		 tmPO->Hprop[k][l][m][n]=Htemp;}
	fok = fscanf(fPF,"%12c%lf%lf",g,&Gtemp,&Htemp);
	tmPO->GinitAT=Gtemp;
	tmPO->HinitAT=Htemp;
	fok = fscanf(fPF,"%12c%lf%lf",g,&Gtemp,&Htemp);
	tmPO->GinitGC=Gtemp;
	tmPO->HinitGC=Htemp;
	for (i=1;i<=48; ++i)
    	{
		fok = fscanf(fPF,"%s%d%d%d%d%lf%lf%lf%lf",g,&k,&l,&m,&n,&Gtemp,&Htemp,&Gttemp,&Httemp);    
		tmPO->Gprop[k][l][m][n]=Gtemp;
		tmPO->Hprop[k][l][m][n]=Htemp;
		tmPO->Gterm[k][l][m][n]=Gttemp;
		tmPO->Hterm[k][l][m][n]=Httemp;}
	for (i=1;i<=32; ++i)
		{fok = fscanf(fPF,"%s %d %d %d %d %lf %lf",g,&k,&l,&m,&n,&Gttemp,&Httemp);
		tmPO->Gterm[k][l][m][n]=Gttemp;
		tmPO->Hterm[k][l][m][n]=Httemp;}
	for (i=1;i<=16; ++i)
		{fok = fscanf(fPF,"%s %d %d %d %d %lf %lf",g,&k,&l,&m,&n,&Gtemp,&Htemp);
		tmPO->Gcoax[k][l][m][n]=Gtemp;
		tmPO->Hcoax[k][l][m][n]=Htemp;}	
	/*parameter array symetrization*/
	a1=1.001;
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
				for (n=0;n<=3;n++)
					/*printf("%d,%d,%d,%d--",k,l,m,n);*/
					if ((tmPO->Gprop[k][l][m][n])==a1) 
					{
					tmPO->Gprop[k][l][m][n]=tmPO->Gprop[m][n][k][l];
					tmPO->Hprop[k][l][m][n]=tmPO->Hprop[m][n][k][l];
					/*printf("rtes");*/
					}
	for (k=0; k<=4; k++)
		for (l=0; l<=4; l++)
			for (m=0;m<=4;m++)
				for (n=0;n<=4;n++)                             
					if (tmPO->Gterm[k][l][m][n]==a1) 
					{
					tmPO->Gterm[k][l][m][n]=tmPO->Gterm[m][n][k][l];
					tmPO->Hterm[k][l][m][n]=tmPO->Hterm[m][n][k][l];
					}
	/*----------Gaps*/
	for (i=0;i<=29; i++)
		{
		fok = fscanf(fPF,"%s %d %lf",g,&k,&Gtemp);
		tmPO->Ggap[k]=Gtemp;
		/**
		printf("%d",k);
		printf ("tes  %f\n",GgapD[k]);
		**/
		}
	fok = fscanf(fPF,"%s %lf",g,&Gtemp);
	tmPO->GATclosepen=Gtemp;
	/*----------*/	    
	/*printf("paramload completed\n");
	printf ("testrtrt%f %f\n",Gterm[3][0][3][4],Gterm[3][4][3][0]);
	scanf("%c");*/
	DB_DSET DB_PrI("<< LoadTmCompleteParsI fok=%d TRUE\n",fok);
	    return(TRUE);	
	}
/*************************************************************/
/*Loads HARDCODED thermodynamic params into data struct*/
/*************************************************************/		
/*initialize all of *your* vars*/
int InitTmCompleteParsI(TM_PARS *tmPO)
	/*array initialization*/
	{int k,l,m,n;
       double a1=1.001;
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
					for (n=0;n<=3;n++)										
					{
					tmPO->Gprop[k][l][m][n]=a1;
					tmPO->Hprop[k][l][m][n]=a1;
					}
	for (k=0; k<=4; k++)
		for (l=0; l<=4; l++)
			for (m=0;m<=4;m++)
				for (n=0;n<=4;n++)                                             
				  {
				  tmPO->Gterm[k][l][m][n]=a1;
				  tmPO->Hterm[k][l][m][n]=a1;
				  } 
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
					for (n=0;n<=3;n++)										
					{
					tmPO->Gcoax[k][l][m][n]=a1;
					tmPO->Hcoax[k][l][m][n]=a1;
					}
      /*array filling with parameters*/
	/*WC:AA/TT */     
	tmPO->Gprop[0][0][3][3]=-1.00;  
	tmPO->Hprop[0][0][3][3]=-7.9;
	/*WC:AT/TA */          
	tmPO->Gprop[0][3][0][3]=-0.88;  
	tmPO->Hprop[0][3][0][3]=-7.2;
	/*WC:TA/AT */          
	tmPO->Gprop[3][0][3][0]=-0.58;  
	tmPO->Hprop[3][0][3][0]=-7.2;
	/*WC:CA/GT */          
	tmPO->Gprop[1][0][3][2]=-1.45;  
	tmPO->Hprop[1][0][3][2]=-8.5;
	/*WC:GT/CA */          
	tmPO->Gprop[2][3][0][1]=-1.44;  
	tmPO->Hprop[2][3][0][1]=-8.4;
	/*WC:CT/GA */          
	tmPO->Gprop[1][3][0][2]=-1.28;  
	tmPO->Hprop[1][3][0][2]=-7.8;
	/*WC:GA/CT */         
	tmPO->Gprop[2][0][3][1]=-1.30;  
	tmPO->Hprop[2][0][3][1]=-8.2;
	/*WC:CG/GC */         
	tmPO->Gprop[1][2][1][2]=-2.17; 
	tmPO->Hprop[1][2][1][2]=-10.6;
	/*WC:GC/CG */          
	tmPO->Gprop[2][1][2][1]=-2.24;  
	tmPO->Hprop[2][1][2][1]=-9.8;
	/*WC:GG/CC */          
	tmPO->Gprop[2][2][1][1]=-1.84;  
	tmPO->Hprop[2][2][1][1]=-8.0;
	/*init:AT */                     
	tmPO->GinitAT=1.03;   
	tmPO->HinitAT=2.3;
	/*init:GC */                     
	tmPO->GinitGC=0.98;   
	tmPO->HinitGC=0.1;
	/*AA:AA/TA */          
	tmPO->Gprop[0][0][0][3]=0.61;   
	tmPO->Hprop[0][0][0][3]=1.2;  
	tmPO->Gterm[0][0][0][3]=-0.67;  
	tmPO->Hterm[0][0][0][3]=-3.1;
	/*AA:CA/GA */          
	tmPO->Gprop[1][0][0][2]=0.43;  
	tmPO->Hprop[1][0][0][2]=-0.9;  
	tmPO->Gterm[1][0][0][2]=-1.01;  
	tmPO->Hterm[1][0][0][2]=-4.3;
	/*AA:GA/CA */          
	tmPO->Gprop[2][0][0][1]=0.17;  
	tmPO->Hprop[2][0][0][1]=-2.9;  
	tmPO->Gterm[2][0][0][1]=-0.99;  
	tmPO->Hterm[2][0][0][1]=-8.0;
	/*AA:TA/AA */          
	tmPO->Gprop[3][0][0][0]=0.69;   
	tmPO->Hprop[3][0][0][0]=4.7;  
	tmPO->Gterm[3][0][0][0]=-0.58;  
	tmPO->Hterm[3][0][0][0]=-2.5;
	/*CC:AC/TC */          
	tmPO->Gprop[0][1][1][3]=1.33;   
	tmPO->Hprop[0][1][1][3]=0.0;  
	tmPO->Gterm[0][1][1][3]=-0.21;  
	tmPO->Hterm[0][1][1][3]=-0.1;
	/*CC:CC/GC */          
	tmPO->Gprop[1][1][1][2]=0.70;  
	tmPO->Hprop[1][1][1][2]=-1.5;  
	tmPO->Gterm[1][1][1][2]=-0.52;  
	tmPO->Hterm[1][1][1][2]=-4.0;
	/*CC:GC/CC */          
	tmPO->Gprop[2][1][1][1]=0.79;   
	tmPO->Hprop[2][1][1][1]=3.6;  
	tmPO->Gterm[2][1][1][1]=-0.62;  
	tmPO->Hterm[2][1][1][1]=-3.9;
	/*CC:TC/AC */          
	tmPO->Gprop[3][1][1][0]=1.05;   
	tmPO->Hprop[3][1][1][0]=6.1;  
	tmPO->Gterm[3][1][1][0]=-0.29;  
	tmPO->Hterm[3][1][1][0]=-0.7;
	/*GG:AG/TG */          
	tmPO->Gprop[0][2][2][3]=-0.13;  
	tmPO->Hprop[0][2][2][3]=-3.1;  
	tmPO->Gterm[0][2][2][3]=-0.42;  
	tmPO->Hterm[0][2][2][3]=-1.1;
	/*GG:CG/GG */          
	tmPO->Gprop[1][2][2][2]=-0.11;  
	tmPO->Hprop[1][2][2][2]=-4.9;  
	tmPO->Gterm[1][2][2][2]=-0.83; 
	tmPO->Hterm[1][2][2][2]=-3.8;
	/*GG:GG/CG */          
	tmPO->Gprop[2][2][2][1]=-1.11;  
	tmPO->Hprop[2][2][2][1]=-6.0;  
	tmPO->Gterm[2][2][2][1]=-0.96;  
	tmPO->Hterm[2][2][2][1]=-0.7;
	/*GG:TG/AG */          
	tmPO->Gprop[3][2][2][0]=0.44;   
	tmPO->Hprop[3][2][2][0]=1.6;  
	tmPO->Gterm[3][2][2][0]=-0.29;  
	tmPO->Hterm[3][2][2][0]=-1.1;
	/*TT:AT/TT */          
	tmPO->Gprop[0][3][3][3]=0.69;  
	tmPO->Hprop[0][3][3][3]=-2.7;  
	tmPO->Gterm[0][3][3][3]=-0.45;  
	tmPO->Hterm[0][3][3][3]=-2.4;
	/*TT:CT/GT */          
	tmPO->Gprop[1][3][3][2]=-0.12;  
	tmPO->Hprop[1][3][3][2]=-5.0;  
	tmPO->Gterm[1][3][3][2]=-0.87;  
	tmPO->Hterm[1][3][3][2]=-6.1;
	/*TT:GT/CT */          
	tmPO->Gprop[2][3][3][1]=0.45;  
	tmPO->Hprop[2][3][3][1]=-2.2;  
	tmPO->Gterm[2][3][3][1]=-0.86;  
	tmPO->Hterm[2][3][3][1]=-7.4;
	/*TT:TT/AT */          
	tmPO->Gprop[3][3][3][0]=0.68;   
	tmPO->Hprop[3][3][3][0]=0.2;  
	tmPO->Gterm[3][3][3][0]=-0.48;  
	tmPO->Hterm[3][3][3][0]=-3.2;
	/*CA:AA/TC */          
	tmPO->Gprop[0][0][1][3]=0.88;   
	tmPO->Hprop[0][0][1][3]=2.3;  
	tmPO->Gterm[0][0][1][3]=-0.35;  
	tmPO->Hterm[0][0][1][3]=-1.6;
	/*CA:AC/TA */          
	tmPO->Gprop[0][1][0][3]=0.77;   
	tmPO->Hprop[0][1][0][3]=5.3;  
	tmPO->Gterm[0][1][0][3]=-0.59;  
	tmPO->Hterm[0][1][0][3]=-1.8;
	/*CA:CA/GC */          
	tmPO->Gprop[1][0][1][2]=0.75;   
	tmPO->Hprop[1][0][1][2]=1.9;  
	tmPO->Gterm[1][0][1][2]=-0.76;  
	tmPO->Hterm[1][0][1][2]=-2.6;
	/*CA:CC/GA */          
	tmPO->Gprop[1][1][0][2]=0.79;   
	tmPO->Hprop[1][1][0][2]=0.6;  
	tmPO->Gterm[1][1][0][2]=-0.85;  
	tmPO->Hterm[1][1][0][2]=-2.7;
	/*CA:GA/CC */          
	tmPO->Gprop[2][0][1][1]=0.81;   
	tmPO->Hprop[2][0][1][1]=5.2;  
	tmPO->Gterm[2][0][1][1]=-0.71;  
	tmPO->Hterm[2][0][1][1]=-5.0;
	/*CA:GC/CA */          
	tmPO->Gprop[2][1][0][1]=0.47;  
	tmPO->Hprop[2][1][0][1]=-0.7;  
	tmPO->Gterm[2][1][0][1]=-1.01;  
	tmPO->Hterm[2][1][0][1]=-3.2;
	/*CA:TA/AC */          
	tmPO->Gprop[3][0][1][0]=0.92;   
	tmPO->Hprop[3][0][1][0]=3.4;  
	tmPO->Gterm[3][0][1][0]=-0.45;  
	tmPO->Hterm[3][0][1][0]=-2.3;
	/*CA:TC/AA */          
	tmPO->Gprop[3][1][0][0]=1.33;   
	tmPO->Hprop[3][1][0][0]=7.6;  
	tmPO->Gterm[3][1][0][0]=-0.55;  
	tmPO->Hterm[3][1][0][0]=-2.7;
	/*CT:AC/TT */          
	tmPO->Gprop[0][1][3][3]=0.64;   
	tmPO->Hprop[0][1][3][3]=0.7;  
	tmPO->Gterm[0][1][3][3]=-0.33;  
	tmPO->Hterm[0][1][3][3]=-0.9;
	/*CT:AT/TC */          
	tmPO->Gprop[0][3][1][3]=0.73;  
	tmPO->Hprop[0][3][1][3]=-1.2;  
	tmPO->Gterm[0][3][1][3]=-0.35;  
	tmPO->Hterm[0][3][1][3]=-2.3;
	/*CT:CC/GT */         
	tmPO->Gprop[1][1][3][2]=0.62;  
	tmPO->Hprop[1][1][3][2]=-0.8;  
	tmPO->Gterm[1][1][3][2]=-0.69;  
	tmPO->Hterm[1][1][3][2]=-3.2;
	/*CT:CT/GC */          
	tmPO->Gprop[1][3][1][2]=0.40;  
	tmPO->Hprop[1][3][1][2]=-1.5;  
	tmPO->Gterm[1][3][1][2]=-0.60;  
	tmPO->Hterm[1][3][1][2]=-3.9;
	/*CT:GC/CT */          
	tmPO->Gprop[2][1][3][1]=0.62;   
	tmPO->Hprop[2][1][3][1]=2.3;  
	tmPO->Gterm[2][1][3][1]=-0.72;  
	tmPO->Hterm[2][1][3][1]=-4.9;
	/*CT:GT/CC */          
	tmPO->Gprop[2][3][1][1]=0.98;   
	tmPO->Hprop[2][3][1][1]=5.2;  
	tmPO->Gterm[2][3][1][1]=-0.61;  
	tmPO->Hterm[2][3][1][1]=-3.0;
	/*
	Bug fix 6/13/03
	tmPO->Gterm[2][3][1][1]=-3.0;
	*/
	/*CT:TC/AT */          
	tmPO->Gprop[3][1][3][0]=0.97;   
	tmPO->Hprop[3][1][3][0]=1.2;  
	tmPO->Gterm[3][1][3][0]=-0.52;  
	tmPO->Hterm[3][1][3][0]=-2.5;
	/*CT:TT/AC */          
	tmPO->Gprop[3][3][1][0]=0.75;   
	tmPO->Hprop[3][3][1][0]=1.0;  
	tmPO->Gterm[3][3][1][0]=-0.34;  
	tmPO->Hterm[3][3][1][0]=-0.7;
	/*GA:AA/TG */          
	tmPO->Gprop[0][0][2][3]=0.14;  
	tmPO->Hprop[0][0][2][3]=-0.6;  
	tmPO->Gterm[0][0][2][3]=-0.52;  
	tmPO->Hterm[0][0][2][3]=-1.9;
	/*GA:AG/TA */          
	tmPO->Gprop[0][2][0][3]=0.02;  
	tmPO->Hprop[0][2][0][3]=-0.7;  
	tmPO->Gterm[0][2][0][3]=-0.65;  
	tmPO->Hterm[0][2][0][3]=-2.5;
	/*GA:CA/GG */          
	tmPO->Gprop[1][0][2][2]=0.03;  
	tmPO->Hprop[1][0][2][2]=-0.7;  
	tmPO->Gterm[1][0][2][2]=-0.88;  
	tmPO->Hterm[1][0][2][2]=-3.9;
	/*GA:CG/GA */          
	tmPO->Gprop[1][2][0][2]=0.11;  
	tmPO->Hprop[1][2][0][2]=-4.0;  
	tmPO->Gterm[1][2][0][2]=-1.23;  
	tmPO->Hterm[1][2][0][2]=-6.0;
	/*GA:GA/CG */          
	tmPO->Gprop[2][0][2][1]=-0.25;  
	tmPO->Hprop[2][0][2][1]=-0.6;  
	tmPO->Gterm[2][0][2][1]=-0.80;  
	tmPO->Hterm[2][0][2][1]=-4.3;
	/*GA:GG/CA */          
	tmPO->Gprop[2][2][0][1]=-0.52;   
	tmPO->Hprop[2][2][0][1]=0.5;  
	tmPO->Gterm[2][2][0][1]=-1.08;  
	tmPO->Hterm[2][2][0][1]=-4.6;
	/*GA:TA/AG */          
	tmPO->Gprop[3][0][2][0]=0.42;   
	tmPO->Hprop[3][0][2][0]=0.6;  
	tmPO->Gterm[3][0][2][0]=-0.53;  
	tmPO->Hterm[3][0][2][0]=-2.0;
	/*GA:TG/AA */          
	tmPO->Gprop[3][2][0][0]=0.74;   
	tmPO->Hprop[3][2][0][0]=3.0;  
	tmPO->Gterm[3][2][0][0]=-0.57;  
	tmPO->Hterm[3][2][0][0]=-2.4;
	/*GT:AG/TT */          
	tmPO->Gprop[0][2][3][3]=0.71;   
	tmPO->Hprop[0][2][3][3]=1.0;  
	tmPO->Gterm[0][2][3][3]=-0.45;  
	tmPO->Hterm[0][2][3][3]=-3.2;
	/*GT:AT/TG */          
	tmPO->Gprop[0][3][2][3]=0.07;  
	tmPO->Hprop[0][3][2][3]=-2.5;  
	tmPO->Gterm[0][3][2][3]=-0.54;  
	tmPO->Hterm[0][3][2][3]=-3.5;
	/*GT:CG/GT */          
	tmPO->Gprop[1][2][3][2]=-0.47;  
	tmPO->Hprop[1][2][3][2]=-4.1;  
	tmPO->Gterm[1][2][3][2]=-0.96;  
	tmPO->Hterm[1][2][3][2]=-3.8;
	/*GT:CT/GG */          
	tmPO->Gprop[1][3][2][2]=-0.32;  
	tmPO->Hprop[1][3][2][2]=-2.8;  
	tmPO->Gterm[1][3][2][2]=-0.81;  
	tmPO->Hterm[1][3][2][2]=-6.6;
	/*GT:GG/CT */          
	tmPO->Gprop[2][2][3][1]=0.08;   
	tmPO->Hprop[2][2][3][1]=3.3;  
	tmPO->Gterm[2][2][3][1]=-0.92;  
	tmPO->Hterm[2][2][3][1]=-5.9;
	/*GT:GT/CG */          
	tmPO->Gprop[2][3][2][1]=-0.59;  
	tmPO->Hprop[2][3][2][1]=-4.4;  
	tmPO->Gterm[2][3][2][1]=-0.40;  
	tmPO->Hterm[2][3][2][1]=-4.4;
	/*GT:TG/AT */          
	tmPO->Gprop[3][2][3][0]=0.43;  
	tmPO->Hprop[3][2][3][0]=-0.1;  
	tmPO->Gterm[3][2][3][0]=-0.59;  
	tmPO->Hterm[3][2][3][0]=-3.9;
	/*GT:TT/AG */          
	tmPO->Gprop[3][3][2][0]=0.34;  
	tmPO->Hprop[3][3][2][0]=-1.3;  
	tmPO->Gterm[3][3][2][0]=-0.59;  
	tmPO->Hterm[3][3][2][0]=-3.6;
	/*DE:AA/T* */          
	tmPO->Gterm[0][0][4][3]=-0.12;  
	tmPO->Hterm[0][0][4][3]=-0.5;
	/*DE:AC/T* */          
	tmPO->Gterm[0][1][4][3]=0.28;   
	tmPO->Hterm[0][1][4][3]=4.7;
	/*DE:AG/T* */          
	tmPO->Gterm[0][2][4][3]=-0.01;  
	tmPO->Hterm[0][2][4][3]=-4.1;
	/*DE:AT/T* */          
	tmPO->Gterm[0][3][4][3]=0.13;  
	tmPO->Hterm[0][3][4][3]=-3.8;
	/*DE:TA/A* */          
	tmPO->Gterm[3][0][4][0]=-0.48;  
	tmPO->Hterm[3][0][4][0]=-0.7;
	/*DE:TC/A* */          
	tmPO->Gterm[3][1][4][0]=-0.19;   
	tmPO->Hterm[3][1][4][0]=4.4;
	/*DE:TG/A* */          
	tmPO->Gterm[3][2][4][0]=-0.50;  
	tmPO->Hterm[3][2][4][0]=-1.6;
	/*DE:TT/A* */          
	tmPO->Gterm[3][3][4][0]=-0.29;   
	tmPO->Hterm[3][3][4][0]=2.9;
	/*DE:CA/G* */          
	tmPO->Gterm[1][0][4][2]=-0.82;  
	tmPO->Hterm[1][0][4][2]=-5.9;
	/*DE:CC/G* */          
	tmPO->Gterm[1][1][4][2]=-0.31;  
	tmPO->Hterm[1][1][4][2]=-2.6;
	/*DE:CG/G* */          
	tmPO->Gterm[1][2][4][2]=-0.01;  
	tmPO->Hterm[1][2][4][2]=-3.2;
	/*DE:CT/G* */          
	tmPO->Gterm[1][3][4][2]=-0.52;  
	tmPO->Hterm[1][3][4][2]=-5.2;
	/*DE:GA/C* */          
	tmPO->Gterm[2 ][0][4][1]=-0.92;  
	tmPO->Hterm[2 ][0][4][1]=-2.1;
	/*DE:GC/C* */          
	tmPO->Gterm[2][1][4][1]=-0.23;  
	tmPO->Hterm[2][1][4][1]=-0.2;
	/*DE:GG/C* */          
	tmPO->Gterm[2][2][4][1]=-0.44;  
	tmPO->Hterm[2][2][4][1]=-3.9;
	/*DE:GT/C* */          
	tmPO->Gterm[2][3][4][1]=-0.35;  
	tmPO->Hterm[2][3][4][1]=-4.4;
	// 	/*DE:AA/*T */          
	tmPO->Gterm[0][0][3][4]=-0.51;   
	tmPO->Hterm[0][0][3][4]=0.2;
	//	/*DE:CA/*T */          
	tmPO->Gterm[1][0][3][4]=-0.42;   
	tmPO->Hterm[1][0][3][4]=-0.6;
	//	/*DE:GA/*T */          
	tmPO->Gterm[2][0][3][4]=-0.62;  
	tmPO->Hterm[2][0][3][4]=-1.1;
	//	/*DE:TA/*T */          
	tmPO->Gterm[3][0][3][4]=-0.71;  
	tmPO->Hterm[3][0][3][4]=-6.9;
	//	/*DE:AT/*A */          
	tmPO->Gterm[0][3][0][4]=-0.50;  
	tmPO->Hterm[0][3][0][4]=-2.9;
	//	/*DE:CT/*A */          
	tmPO->Gterm[1][3][0][4]=-0.02;  
	tmPO->Hterm[1][3][0][4]=-4.1;
	//	/*DE:GT/*A */          
	tmPO->Gterm[2][3][0][4]=0.48;  
	tmPO->Hterm[2][3][0][4]=-4.2;
	//	/*DE:TT/*A */          
	tmPO->Gterm[3][3][0][4]=-0.10;  
	tmPO->Hterm[3][3][0][4]=-0.2;
	//	/*DE:AC/*G */          
	tmPO->Gterm[0][1][2][4]=-0.96;  
	tmPO->Hterm[0][1][2][4]=-6.3;
	//	/*DE:CC/*G */          
	tmPO->Gterm[1][1][2][4]=-0.52;  
	tmPO->Hterm[1][1][2][4]=-4.4;
	//	/*DE:GC/*G */          
	tmPO->Gterm[2][1][2][4]=-0.72;  
	tmPO->Hterm[2][1][2][4]=-5.1;
	//	/*DE:TC/*G */          
	tmPO->Gterm[0][2][1][4]=-0.58;  
	tmPO->Hterm[0][2][1][4]=-4.0;
	//	/*DE:AG/*C */          
	tmPO->Gterm[0][2][1][4]=-0.58;  
	tmPO->Hterm[0][2][1][4]=-3.7;
	//	/*DE:CG/*C */          
	tmPO->Gterm[1][2][1][4]=-0.34;  
	tmPO->Hterm[1][2][1][4]=-4.0;
	//	/*DE:GG/*C */          
	tmPO->Gterm[2][2][1][4]=-0.56;  
	tmPO->Hterm[2][2][1][4]=-3.9;
	//	/*DE:TG/*C */          
	tmPO->Gterm[3][2][1][4]=-0.61;  
	tmPO->Hterm[3][2][1][4]=-4.9;
	/*CS:A-A//T/T */      
	tmPO->Gcoax[0][0][3][3]=-1.42; 
	tmPO->Hcoax[0][0][3][3]=-14.6;
	/*CS:T-T//A/A */       
	tmPO->Gcoax[3][3][0][0]=-2.49; 
	tmPO->Hcoax[3][3][0][0]=-14.3;
	/*CS:A-T//T/A */       
	tmPO->Gcoax[0][3][0][3]=-2.40; 
	tmPO->Hcoax[0][3][0][3]=-15.1;
	/*CS:T-A//A/T */       
	tmPO->Gcoax[3][0][3][0]=-1.10;  
	tmPO->Hcoax[3][0][3][0]=-9.0;
	/*CS:C-A//G/T */       
	tmPO->Gcoax[1][0][3][2]=-0.94; 
	tmPO->Hcoax[1][0][3][2]=-10.0;
	/*CS:T-G//A/C */       
	tmPO->Gcoax[3][2][0][1]=-1.18;  
	tmPO->Hcoax[3][2][0][1]=-9.6;
	/*CS:G-T//C/A */       
	tmPO->Gcoax[2][3][0][1]=-2.73; 
	tmPO->Hcoax[2][3][0][1]=-15.1;
	/*CS:A-C//T/G */       
	tmPO->Gcoax[0][1][2][3]=-2.12; 
	tmPO->Hcoax[0][1][2][3]=-12.7;
	/*CS:C-T//G/A */       
	tmPO->Gcoax[1][3][0][2]=-1.97;  
	tmPO->Hcoax[1][3][0][2]=-7.6;
	/*CS:A-G//T/C */       
	tmPO->Gcoax[0][2][1][3]=-1.78; 
	tmPO->Hcoax[0][2][1][3]=-14.2;
	/*CS:G-A//C/T */       
	tmPO->Gcoax[2][0][3][1]=-2.60; 
	tmPO->Hcoax[2][0][3][1]=-17.7;
	/*CS:T-C//A/G */       
	tmPO->Gcoax[3][1][2][0]=-2.61;  
	tmPO->Hcoax[3][1][2][0]=-15.6;
	/*CS:C-G//G/C */       
	tmPO->Gcoax[1][2][1][2]=-0.75;  
	tmPO->Hcoax[1][2][1][2]=-5.0;
	/*CS:G-C//C/G */       
	tmPO->Gcoax[2][1][2][1]=-2.23; 
	tmPO->Hcoax[2][1][2][1]=-14.7;
	/*CS:C-C//G/G */       
	tmPO->Gcoax[1][1][2][2]=-2.60; 
	tmPO->Hcoax[1][1][2][2]=-18.7;
	/*CS:G-G//C/C */       
	tmPO->Gcoax[2][2][1][1]=-1.14;  
	tmPO->Hcoax[2][2][1][1]=-6.8;
	/*parameter array symetrization*/
	a1=1.001;
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
				for (n=0;n<=3;n++)
					/*printf("%d,%d,%d,%d--",k,l,m,n);*/
					if ((tmPO->Gprop[k][l][m][n])==a1) 
					{
					tmPO->Gprop[k][l][m][n]=tmPO->Gprop[m][n][k][l];
					tmPO->Hprop[k][l][m][n]=tmPO->Hprop[m][n][k][l];
					/*printf("rtes");*/
					}
	for (k=0; k<=4; k++)
		for (l=0; l<=4; l++)
			for (m=0;m<=4;m++)
				for (n=0;n<=4;n++)                             
					if (tmPO->Gterm[k][l][m][n]==a1) 
					{
					tmPO->Gterm[k][l][m][n]=tmPO->Gterm[m][n][k][l];
					tmPO->Hterm[k][l][m][n]=tmPO->Hterm[m][n][k][l];
	
							}
	/*----------Gaps*/
	tmPO->Ggap[0]=4.0;
	tmPO->Ggap[1]=2.9;
	tmPO->Ggap[2]=3.1;
	tmPO->Ggap[3]=3.2;
	tmPO->Ggap[4]=3.3;
	tmPO->Ggap[5]=3.5;
	tmPO->Ggap[6]=3.7;
	tmPO->Ggap[7]=3.8;
	tmPO->Ggap[8]=3.9;
	tmPO->Ggap[9]=4.0;
	tmPO->Ggap[10]=4.1;	      
	tmPO->Ggap[11]=4.2;
	tmPO->Ggap[12]=4.2;
	tmPO->Ggap[13]=4.3;
	tmPO->Ggap[14]=4.3;
	tmPO->Ggap[15]=4.4;
	tmPO->Ggap[16]=4.4;
	tmPO->Ggap[17]=4.5;
	tmPO->Ggap[18]=4.6;
	tmPO->Ggap[19]=4.6;
	tmPO->Ggap[20]=4.7;
	tmPO->Ggap[21]=4.7;
	tmPO->Ggap[22]=4.8;
	tmPO->Ggap[23]=4.8;
	tmPO->Ggap[24]=4.9;
	tmPO->Ggap[25]=4.9;
	tmPO->Ggap[26]=5.0;
	tmPO->Ggap[27]=5.1;
	tmPO->Ggap[28]=5.1;
	tmPO->Ggap[29]=5.2;
	tmPO->GATclosepen=0.5;   
	/*--------*/
	return(TRUE);
}
/*************************************************************/
/*Reports thermodynami & hybridization parameters*/
/*************************************************************/	
/*Dump out (write to file outPF) all the thermodynamic*/ 
/*parameters you ever care to report*/

void DumpCompleteTmPars(TM_PARS *tmPO,FILE *outPF)
	{
	int k,l,m,n;
	double Gtemp, Htemp,Gtemps, Htemps;

	HAND_NFILE(outPF);
	fprintf(outPF,"# Parameters in TM_PARS\n");
	fprintf(outPF,"# Hybridization parameters\n");
	fprintf(outPF,"# [top strand](M)/[bottom strand](M)/[cation+](M)/[Mg++](M)/temperature(C)\n");
	fprintf(outPF,"# %.2e/%.2e/%.2e/%.2e/%.2f\n",tmPO->conc,tmPO->conc2,tmPO->salt,tmPO->mg,tmPO->tp);
	fprintf(outPF,"# Nearest-neighbor parameters\n");
        fprintf(outPF,"# Watson-Crick, internal mismatches\n");
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
				for (n=0;n<=3;n++)										
					{
					Gtemp=tmPO->Gprop[k][l][m][n];
					Htemp=tmPO->Hprop[k][l][m][n];
					Gtemps=tmPO->Gprop[m][n][k][l];
					Htemps=tmPO->Hprop[m][n][k][l];

					if ((Gtemp!=1.001)||(Htemp !=1.001))
						fprintf(outPF,"# %d %d %d %d %.2f %.1f %.2f %.1f\n",k,l,m,n,Gtemp,Htemp,Gtemps,Htemps);
					}
	fprintf(outPF,"# initiation\n");
					Gtemp=tmPO->GinitAT;
					Htemp=tmPO->HinitAT;
					Gtemps=tmPO->GinitGC;
					Htemps=tmPO->HinitGC;
                    fprintf(outPF,"# AT %.2f %.1f GC %.2f %.1f\n",Gtemp,Htemp,Gtemps,Htemps);
	fprintf(outPF,"# Terminal mismatches and dandling ends\n");
	for (k=0; k<=4; k++)
		for (l=0; l<=4; l++)
			for (m=0;m<=4;m++)
				for (n=0;n<=4;n++)                                             
				  {
					Gtemp=tmPO->Gterm[k][l][m][n];
					Htemp=tmPO->Hterm[k][l][m][n];
					Gtemps=tmPO->Gterm[m][n][k][l];
					Htemps=tmPO->Hterm[m][n][k][l];
					if ((Gtemp!=1.001)||(Htemp !=1.001))
						fprintf(outPF,"# %d %d %d %d %.2f %.1f %.2f %.1f\n",k,l,m,n,Gtemp,Htemp, Gtemps, Htemps);
				  } 
	fprintf(outPF,"# Coaxial stacking \n");
	for (k=0; k<=3; k++)
		for (l=0; l<=3; l++)
			for (m=0;m<=3;m++)
				for (n=0;n<=3;n++)										
					{
					Gtemp=tmPO->Gcoax[k][l][m][n];
					Htemp=tmPO->Hcoax[k][l][m][n];
					if ((Gtemp!=1.001)||(Htemp !=1.001))
						fprintf(outPF,"# %d %d %d %d %.2f %.1f\n",k,l,m,n,Gtemp,Htemp);
	/*----------Gaps*/				}
	fprintf(outPF,"# Gaps \n");
	for (k=0; k<=29; k++)
		{
		Gtemp=tmPO->Ggap[k];
		fprintf(outPF,"# %d  %.2f \n",k+1,Gtemp);
		}	
	fprintf(outPF,"# AT closing penalty for gaps \n");
	Gtemp=tmPO->GATclosepen;
	/*----------*/  
	fprintf(outPF,"# %.2f \n",Gtemp);
}
/*******************************************************************/
int SeqToNumsI(char *seqS,int *numsPI,int max,int mode,int verb)
{
	int i,n;

	/*non a,A,C,C,t,T,u,U,g,G,s,S,- characters are edited*/
	/* A=0,C=1,G=2,T=U=3,S=6,-=7*/	
	n = 0;
	for (i=0;i<max;i++)
	{
		switch(seqS[i])
		{
			case 'A': case 'a':
				numsPI[n++] = 0;
				break;
			case 'C': case 'c':
				numsPI[n++] = 1;
				break;
			case 'G': case 'g':
				numsPI[n++] = 2;
				break;
			case 'T': case 't':
			case 'U': case 'u':
				numsPI[n++] = 3;
				break;
			case '*':
				switch(mode)
				{
					case DSET_EXPNOCOAX:
					case DSET_EXPCOAX:
						numsPI[n++] = 4;
						break;
					default:
						if(verb)
							printf("# Error: bad char '*' (mode %d)\n",mode);
						return(BOGUS);
				}
				break;
			case '-':
				switch(mode)
				{
					case DSET_EXPNOCOAX:
					case DSET_EXPCOAX:
						numsPI[n++] = 7;
						break;
					default:
						if(verb)
							printf("# Error: bad char '-' (mode %d)\n",mode);
						return(BOGUS);
				}
				break;
			case '/':
				switch(mode)
				{
					case DSET_EXPCOAX:
						numsPI[n++] = 6;
						break;
					default:
						if(verb)
							printf("# Error: bad char '/' (mode %d)\n",mode);
						return(BOGUS);
				}
				break;
			default:
				if(verb)
					printf("# Error: bad char '%c'\n",seqS[i]);
				return(BOGUS);
		}
	}
	return(n);
}
