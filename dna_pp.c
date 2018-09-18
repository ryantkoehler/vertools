/*
* dna_pp.c
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
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "prim.h"
#include "dna.h"
#include "dna_pair.h"


#define DB_PP   if(DB[115])


/************************************************************************ ppp
*   Create PPARS data structure 
*/
PPARS *CreatePparsPO()
{
    PPARS *pparPO;

    if(!(pparPO = (PPARS *)ALLOC(1,sizeof(PPARS))))
    {   
        return(NULL);
    }
    pparPO->ID = PPARS_ID;
    InitPpars(pparPO);
    return(pparPO); 
}
/************************************************************************
*   Frees up PPARS data structure
*/
int DestroyPparsI(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    FREE(pparPO);
    return(TRUE);
}
/***********************************************************************
*   Initialize PPARS data structure vars
*/
void InitPpars(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    INIT_S(pparPO->parfile);
    pparPO->alg = DEF_ALGO;
    pparPO->ctype = DEF_CTYPE;
    pparPO->rval = DEF_RVAL;
    pparPO->min_word = DEF_MIN_WORD;
    pparPO->min_loop = DEF_MIN_LOOP;
    pparPO->cl3 = DEF_CLAMP3;
    pparPO->cl5 = DEF_CLAMP5;
    pparPO->thresh = DEF_C_THRESH;
    INIT_S(pparPO->msfile);
    SetPair_parIdentMatch(pparPO);
    pparPO->do_ham = FALSE;
    return;
}
/***********************************************************************
*   Set pair-wise match weights for identity check
*/
void SetPair_parIdentMatch(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->mscore[0]   = 1.0;      /* AA */
    pparPO->mscore[1]   = 0.0;      /* AC */
    pparPO->mscore[2]   = 0.0;      /* AG */
    pparPO->mscore[3]   = 0.0;      /* AT */
    pparPO->mscore[4]   = 0.0;      /* CA */
    pparPO->mscore[5]   = 1.0;      /* CC */
    pparPO->mscore[6]   = 0.0;      /* CG */
    pparPO->mscore[7]   = 0.0;      /* CT */
    pparPO->mscore[8]   = 0.0;      /* GA */
    pparPO->mscore[9]   = 0.0;      /* GC */
    pparPO->mscore[10]  = 1.0;      /* GG */
    pparPO->mscore[11]  = 0.0;      /* GT */
    pparPO->mscore[12]  = 0.0;      /* TA */
    pparPO->mscore[13]  = 0.0;      /* TC */
    pparPO->mscore[14]  = 0.0;      /* TG */
    pparPO->mscore[15]  = 1.0;      /* TT */
    pparPO->parset = PP_PSET_ID;
    SetPair_parMaxScores(pparPO);
}
/***********************************************************************
*   Set pair-wise match / mismatch weights for default complement match
*/
void SetPair_parCompMatch(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->mscore[0]   = 0.0;      /* AA */
    pparPO->mscore[1]   = 0.0;      /* AC */
    pparPO->mscore[2]   = 0.0;      /* AG */
    pparPO->mscore[3]   = 1.0;      /* AT */
    pparPO->mscore[4]   = 0.0;      /* CA */
    pparPO->mscore[5]   = 0.0;      /* CC */
    pparPO->mscore[6]   = 1.0;      /* CG */
    pparPO->mscore[7]   = 0.0;      /* CT */
    pparPO->mscore[8]   = 0.0;      /* GA */
    pparPO->mscore[9]   = 1.0;      /* GC */
    pparPO->mscore[10]  = 0.0;      /* GG */
    pparPO->mscore[11]  = 0.0;      /* GT */
    pparPO->mscore[12]  = 1.0;      /* TA */
    pparPO->mscore[13]  = 0.0;      /* TC */
    pparPO->mscore[14]  = 0.0;      /* TG */
    pparPO->mscore[15]  = 0.0;      /* TT */
    pparPO->parset = PP_PSET_COM;
    SetPair_parMaxScores(pparPO);
}
/***********************************************************************
*   Set Thermodynamic pair-wise match / mismatch weights
*   1/20/03: Values are -dG from Nic Peyret's "averaging" of missmatchs
*/
void SetPair_parThermoWmatch(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->mscore[0]   = -0.5;     /* AA */
    pparPO->mscore[1]   = -0.85;    /* AC */
    pparPO->mscore[2]   = -0.1;     /* AG */
    pparPO->mscore[3]   = 1.1;      /* AT */
    pparPO->mscore[4]   = -0.85;    /* CA */
    pparPO->mscore[5]   = -0.95;    /* CC */
    pparPO->mscore[6]   = 1.7;      /* CG */
    pparPO->mscore[7]   = -0.7;     /* CT */
    pparPO->mscore[8]   = -0.1;     /* GA */
    pparPO->mscore[9]   = 1.7;      /* GC */
    pparPO->mscore[10]  = -0.15;    /* GG */
    pparPO->mscore[11]  = -0.05;    /* GT */
    pparPO->mscore[12]  = 1.1;      /* TA */
    pparPO->mscore[13]  = -0.7;     /* TC */
    pparPO->mscore[14]  = -0.05;    /* TG */
    pparPO->mscore[15]  = -0.45;    /* TT */
    pparPO->parset = PP_PSET_TH;
    SetPair_parMaxScores(pparPO);
}
/***********************************************************************
*   Set GC=3 AT=2 pair-wise match weights
*/
void SetPair_parS3W2Wmatch(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->mscore[0]   = 0.0;      /* AA */
    pparPO->mscore[1]   = 0.0;      /* AC */
    pparPO->mscore[2]   = 0.0;      /* AG */
    pparPO->mscore[3]   = 2.0;      /* AT */
    pparPO->mscore[4]   = 0.0;      /* CA */
    pparPO->mscore[5]   = 0.0;      /* CC */
    pparPO->mscore[6]   = 3.0;      /* CG */
    pparPO->mscore[7]   = 0.0;      /* CT */
    pparPO->mscore[8]   = 0.0;      /* GA */
    pparPO->mscore[9]   = 3.0;      /* GC */
    pparPO->mscore[10]  = 0.0;      /* GG */
    pparPO->mscore[11]  = 0.0;      /* GT */
    pparPO->mscore[12]  = 2.0;      /* TA */
    pparPO->mscore[13]  = 0.0;      /* TC */
    pparPO->mscore[14]  = 0.0;      /* TG */
    pparPO->mscore[15]  = 0.0;      /* TT */
    pparPO->parset = PP_PSET_TH;
    SetPair_parMaxScores(pparPO);
}
/***************************************************************************
*   Return match score for pair of chars
*/
int GetPair_parMatchScoreI(PPARS *pparPO, char f, char s, REAL *scPR)
{
    int pind;

    VALIDATE(pparPO,PPARS_ID);
    pind = SeqBasePairIndexI(f,s);
    if(IS_BOG(pind))
    {
        return(FALSE);
    }
    if(scPR)
    {
        *scPR = pparPO->mscore[pind];
    }
    return(TRUE);
}
/***************************************************************************
*   Load pair-wise pairing parameters 
*   Returns the TRUE if all parameters are loaded; FALSE if problem
*/
int LoadPairingParsI(char *inS,PPARS *pparPO)
{
    int i,prob;
    char bufS[DEF_BS],pairS[DEF_BS];
    FILE *fPF;
    REAL pR, pvalsR[16];
    
    DB_PP DB_PrI(">> LoadPairingParsI |%s|\n",inS);
    if(!(fPF=OpenUFilePF(inS,"r",NULL)))
    {
        DB_PP DB_PrI("<< LoadPairingParsI FALSE\n");
        return(FALSE);
    }
    /***
    *   Init parameters to bogus values
    */
    for(i=0;i<PP_IND_NUM;i++)
    {
        pvalsR[i] = BAD_R;
    }
    /***
    *   Load data values with assumed format XX Pairing
    *       where XX is the sequence pair & XX is a pairing score value
    */
    DB_PP DB_PrI("+ Loading data values\n");
    prob = 0;
    while(fgets(bufS,LINEGRAB,fPF) != NULL)
    {
        if(COM_LINE(bufS))
        {   continue;   }
        pR = BAD_R;
        INIT_S(pairS);
        sscanf(bufS,"%s %lf",pairS,&pR);
        if(NO_S(pairS))
        {   continue;   }
        Upperize(pairS);
        /***
        *   Either code-pair keys or specific init keywords
        */
        i = SeqPairIndexI(pairS); 
        if(BAD_REAL(pR) || IS_BOG(i))
        {
            PROBLINE;
            printf("Bad / unrecognized data line:\n");
            fputs(bufS,stdout);
            prob++; break;
        }
        pvalsR[i] = pR;
    }
    FILECLOSE(fPF);
    if(prob)
    {
        return(FALSE);
    }
    /***
    *   Check all parameters and complain if any missed 
    */
    for(i=0;i<PP_IND_NUM;i++)
    {
        if( BAD_REAL(pvalsR[i]) )
        {
            PROBLINE;
            pairS[0] = DNAIndexBaseC(i/4);
            pairS[1] = DNAIndexBaseC(i%4);
            pairS[2] = '\0';
            printf("Failed to get score value for pairing %s\n",pairS); 
            DB_PP DB_PrI("<< LoadPairingParsI FALSE\n");
            return(FALSE);
        }
    }
    /***
    *   Finally copy into data structure
    */
    for(i=0;i<PP_IND_NUM;i++)
    {
        pparPO->mscore[i] = pvalsR[i];
    }
    strcpy(pparPO->parfile,inS);
    SetPair_parMaxScores(pparPO);
    DB_PP DB_PrI("<< LoadPairingParsI TRUE\n");
    return(TRUE);
}
/***********************************************************************
*   Set max values of score for each base type
*/
void SetPair_parMaxScores(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    /* A-X */
    pparPO->mA = MAX_NUM(pparPO->mscore[0],pparPO->mscore[1]);
    pparPO->mA = MAX_NUM(pparPO->mscore[2],pparPO->mA);
    pparPO->mA = MAX_NUM(pparPO->mscore[3],pparPO->mA);
    /* C-X */
    pparPO->mC = MAX_NUM(pparPO->mscore[4],pparPO->mscore[5]);
    pparPO->mC = MAX_NUM(pparPO->mscore[6],pparPO->mC);
    pparPO->mC = MAX_NUM(pparPO->mscore[7],pparPO->mC);
    /* G-X */
    pparPO->mG = MAX_NUM(pparPO->mscore[8],pparPO->mscore[9]);
    pparPO->mG = MAX_NUM(pparPO->mscore[10],pparPO->mG);
    pparPO->mG = MAX_NUM(pparPO->mscore[11],pparPO->mG);
    /* T-X */
    pparPO->mT = MAX_NUM(pparPO->mscore[12],pparPO->mscore[13]);
    pparPO->mT = MAX_NUM(pparPO->mscore[14],pparPO->mT);
    pparPO->mT = MAX_NUM(pparPO->mscore[15],pparPO->mT);
}
/****************************************************************************
*   Set pair-wise comparision algorithm
*/
int SetPparsAlgI(PPARS *pparPO,int alg)
{
    VALIDATE(pparPO,PPARS_ID);
    switch(alg)
    {
        case ALGO_BMATCH:       
        case ALGO_CONT:     
        case ALGO_CON3:     
        case ALGO_CON5:     
        case ALGO_MATCH:    
            pparPO->alg = alg;
            return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************/
int GetPparsAlgI(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    return(pparPO->alg);
}
/***************************************************************************/
void FillPparsAlgTagString(int alg,char *bufS)
{
    switch(alg)
    {
        case ALGO_BMATCH:       sprintf(bufS,"SCB");    break;
        case ALGO_CONT:         sprintf(bufS,"SCON");   break;
        case ALGO_MATCH:        sprintf(bufS,"SMAT");   break;
        case ALGO_CON3:         sprintf(bufS,"SCO3");   break;
        case ALGO_CON5:         sprintf(bufS,"SCO5");   break;
        default:                sprintf(bufS,">>>BOGUS<<<");    break;
    }
}
/**************************************************************************/
void FillPparsAlgString(int alg,char *bufS)
{
    switch(alg)
    {
        case ALGO_MATCH:    
            sprintf(bufS,"Max matching");   break;
        case ALGO_BMATCH:       
            sprintf(bufS,"Block-weighted matching");    break;
        case ALGO_CONT:     
            sprintf(bufS,"Max contiguous stretch"); break;
        case ALGO_CON3:     
            sprintf(bufS,"Max contiguous stretch @ 3' end");    break;
        case ALGO_CON5:     
            sprintf(bufS,"Max contiguous stretch @ 5' end");    break;
        default:            
            sprintf(bufS,"??? DEFAULT ???");    break;
    }
}
/**************************************************************************
*   Set return-value type
*/
int SetPparsRvalI(PPARS *pparPO,int type)
{
    VALIDATE(pparPO,PPARS_ID);
    switch(type)
    {
        case PP_R_MAX:
        case PP_R_TOT:
        case PP_R_NUM:
        case PP_R_EFR:
            pparPO->rval = type;
            return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************/
int GetPparsRvalI(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    return(pparPO->rval);
}
/**************************************************************************/
void FillPparsRvalString(int rval,char *bufS)
{
    switch(rval)
    {
        case PP_R_MAX:  sprintf(bufS,"Maximum value");      break;
        case PP_R_TOT:  sprintf(bufS,"Total value (at least thresh)");  break;
        case PP_R_NUM:  sprintf(bufS,"Alignment count (at least thresh)");break;
        case PP_R_EFR:  sprintf(bufS,"Expect/Find");    break;
        default:        sprintf(bufS,">>>BOGUS<<<");    break;
    }
}
/**************************************************************************
*   Set comparision type
*/
int SetPparsCtypeI(PPARS *pparPO,int type)
{
    VALIDATE(pparPO,PPARS_ID);
    switch(type)
    {
        case PP_COM_SIM:
        case PP_COM_COM:
        case PP_COM_LOOP:
        case PP_COM_PCOM:
            pparPO->ctype = type;
            return(TRUE);
    }
    return(FALSE);
}
/**************************************************************************/
int GetPparsCtypeI(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    return(pparPO->ctype);
}
/**************************************************************************/
void FillPparsCtypeString(int ctype,char *bufS)
{
    switch(ctype)
    {
        case PP_COM_SIM:    sprintf(bufS,"Similarity");         break;
        case PP_COM_COM:    sprintf(bufS,"Complimentarity");    break;
        case PP_COM_LOOP:   sprintf(bufS,"Hairpin Loop");       break;
        case PP_COM_PCOM:   sprintf(bufS,"Parallel Comp");  break;
        default:        sprintf(bufS,">>>BOGUS<<<");    break;
    }
}
/*************************************************************************/
int SetPparsThreshI(PPARS *pparPO,REAL tR)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->thresh = tR;
    return(TRUE);
}
/*************************************************************************/
REAL GetPparsThreshR(PPARS *pparPO)
{
    VALIDATE(pparPO,PPARS_ID);
    return(pparPO->thresh);
}
/*************************************************************************/
REAL MaxSeqCompareScoreR(PPARS *pparPO,char *seqPC,int len)
{
    int i;
    REAL scR;

    VALIDATE(pparPO,PPARS_ID);
    scR = 0.0;
    for(i=0;i<len;i++)
    {
        switch(UPPER(seqPC[i]))
        {
            case 'A':   scR += pparPO->mA;  break;
            case 'C':   scR += pparPO->mC;  break;
            case 'G':   scR += pparPO->mG;  break;
            case 'T':   scR += pparPO->mT;  break;
        }
    }
    return(scR);
}
/****************************************************************************
*   Set pair-wise minimum word size for comparision
*/
int SetPparsWordMinI(PPARS *pparPO,int v)
{
    VALIDATE(pparPO,PPARS_ID);
    if(v > 0) {
        pparPO->min_word = v;
        return(TRUE);
    }
    return(FALSE);
}
/***************************************************************************/
int SetPparsHammingI(PPARS *pparPO,int v)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->do_ham = v;
    return(TRUE);
}
/****************************************************************************
*   Set clamp value
*/
int SetPparsClampsI(PPARS *pparPO,int v3,int v5)
{
    VALIDATE(pparPO,PPARS_ID);
    pparPO->cl3 = v3;
    pparPO->cl5 = v5;
    return(TRUE);
}
/****************************************************************************
*   Set pair-wise loop minimum for loop evaluations
*/
int SetPparsLoopMinI(PPARS *pparPO,int v)
{
    VALIDATE(pparPO,PPARS_ID);
    if(v > 0)
    {
        pparPO->min_loop = v;
        return(TRUE);
    }
    return(FALSE);
}
/****************************************************************************
*   Dump PPARS data structure vars
*/
void DumpPpars(PPARS *pparPO,FILE *oPF)
{
    char pairS[DEF_BS];

    VALIDATE(pparPO,PPARS_ID);
    HAND_NFILE(oPF);
    fprintf(oPF,"######################################################\n");
    fprintf(oPF,"#         Parameters for Sequence Comparison \n");
    fprintf(oPF,"######################################################\n");
    fprintf(oPF,"# Comparison:    ");   
    FillPparsCtypeString(pparPO->ctype,pairS);
    fprintf(oPF,"%s\n",pairS);  
    fprintf(oPF,"# Algorithm:     ");   
    FillPparsAlgString(pparPO->alg,pairS);
    fprintf(oPF,"%s\n",pairS);  
    if(pparPO->do_ham) {
        fprintf(oPF,"# Hamming distance (no align)\n"); 
    }     
    else {
        fprintf(oPF,"# Non-gap align (not Hamming)\n"); 
    }
    fprintf(oPF,"# Return value:  ");   
    FillPparsRvalString(pparPO->rval,pairS);
    fprintf(oPF,"%s\n",pairS);  
    fprintf(oPF,"# Threshold:     %3.2f\n",pparPO->thresh); 
    fprintf(oPF,"# Min word size: %d\n",pparPO->min_word);  
    fprintf(oPF,"# Min loop size: %d\n",pparPO->min_loop);  
    fprintf(oPF,"# Clamp 3':      %d\n",pparPO->cl3);   
    fprintf(oPF,"# Clamp 5':      %d\n",pparPO->cl5);   
    DumpPparsWeights(pparPO,oPF);
    fprintf(oPF,"#\n");
}
/****************************************************************************/
void DumpPparsWeights(PPARS *pparPO,FILE *oPF)
{
    int i;
    char pairS[DEF_BS];

    VALIDATE(pparPO,PPARS_ID);
    HAND_NFILE(oPF);
    fprintf(oPF,"# Pair weightings\n");
    if(!NO_S(pparPO->parfile))
    {
        fprintf(oPF,"#   From parameter file:    %s\n",pparPO->parfile);    
    }
    else
    {
        switch(pparPO->parset)
        {
            case PP_PSET_ID:
                fprintf(oPF,"#   Base identity\n"); 
                break;
            case PP_PSET_COM:
                fprintf(oPF,"#   Base complement\n");   
                break;
            case PP_PSET_TH:
                fprintf(oPF,"#   Thermodynamic defaults\n");    
                break;
            default:
                printf("Bogus parset = %d\n",pparPO->parset);
                ERR("DumpPparsWeights","Bad parset code");
            
        }
    }
    for(i=0;i<PP_IND_NUM;i++)
    {
        pairS[0] = DNAIndexBaseC(i/4);
        pairS[1] = DNAIndexBaseC(i%4);
        pairS[2] = '\0';
        fprintf(oPF,"#    %s = %4.3f\n",pairS,pparPO->mscore[i]);
    }
}
