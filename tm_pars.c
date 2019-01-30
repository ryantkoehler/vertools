/*
* tm_pars.c
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
#include <math.h>
#include "prim.h"
#include "dna.h"
#include "table.h"
#include "tm_pars.h"
/*
#include "dset_tm.h"
*/

#define DB_DNA_TM   if(DB[110])

/***************************************************************************
*   Allocate Tm data structure, and call init 
*/
TM_PARS *CreateTm_parsPO()
{
    TM_PARS *tmPO;

    DB_DNA_TM DB_PrI(">> CreateTm_parsPO\n");
    if(!(tmPO = (TM_PARS *)ALLOC(1,sizeof(TM_PARS))))
    {
        printf("Failed to allocate TM_PARS structure\n");
        return(NULL);
    }
    tmPO->ID = TM_PARS_ID;
    DB_DNA_TM DB_PrI("+ allocated %p size(%d)\n",tmPO,sizeof(TM_PARS));
    InitTmPars(tmPO);
    DB_DNA_TM DB_PrI("<< CreateTm_parsPO %p\n",tmPO);
    return(tmPO);
}
/***************************************************************************
*   Free Tm data structure
*/
int DestroyTm_parsI(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> DestroyTm_parsI %p\n",tmPO);
    VALIDATE(tmPO,TM_PARS_ID);
    CHECK_TABLE(tmPO->htab);
    CHECK_TABLE(tmPO->stab);
    FREE(tmPO);
    DB_DNA_TM DB_PrI("<< DestroyTm_parsI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Initialization Tm data structure to default values
*/
void InitTmPars(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> InitTmPars %p\n",tmPO);
    VALIDATE(tmPO,TM_PARS_ID);
    DB_DNA_TM DB_PrI("+ valid\n");
    /***
    *   Set default values
    */
    tmPO->algo = DEF_TM_ALGO;
    INIT_S(tmPO->source);
    tmPO->conc = tmPO->conc2 = tmPO->conc3 = DEF_NN_CONC;
    tmPO->salt = DEF_NN_SALT;
    tmPO->mg = DEF_NN_MG;
    tmPO->tp = DEF_NN_TEMP;
    DB_DNA_TM 
    {
        DB_PrI("+ Conc %e\n",tmPO->conc);
        DB_PrI("+ Salt %e\n",tmPO->salt);
        DB_PrI("+ Mg   %e\n",tmPO->mg);
    }
    tmPO->ambig_ok = DEF_AMBIG_OK;
    InitThermoPars(tmPO);
    DB_DNA_TM DB_PrI("<< InitTmPars\n");
}
/**************************************************************************
*   Initialize parameters for simple NN models
*/
void InitThermoPars(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> InitThermoPars algo=%d\n",tmPO->algo);
    ClearTmParNNTables(tmPO);
    SetTmParAlgo(tmPO,DEF_TM_ALGO,TRUE);
    UpdateTmParsNNTabVars(tmPO);
    DB_DNA_TM DB_PrI("<< InitThermoPars\n");
}
/**************************************************************************
*   Clear any thermo tables 
*/
void ClearTmParNNTables(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> ClearTmParNNTables\n");
    VALIDATE(tmPO,TM_PARS_ID);
    CHECK_TABLE(tmPO->htab);
    CHECK_TABLE(tmPO->stab);
    UpdateTmParsNNTabVars(tmPO);
    DB_DNA_TM DB_PrI("<< ClearTmParNNTables\n");
}
/**************************************************************************
*   Create thermo tables with specified row col dimensions
*/
int InitTmParNNTablesI(TM_PARS *tmPO, int r, int c)
{
    DB_DNA_TM DB_PrI(">> InitTmParNNTablesI r=%d c=%d\n",r,c);
    VALIDATE(tmPO,TM_PARS_ID);
    CHECK_TABLE(tmPO->htab);
    CHECK_TABLE(tmPO->stab);
    tmPO->htab = CreateTablePO(r,c);
    tmPO->stab = CreateTablePO(r,c);
    PrettyTmParTableForm(tmPO);
    DB_DNA_TM DB_PrI("<< InitTmParNNTablesI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Set TM_PARS values for thermo tables
*/
void UpdateTmParsNNTabVars(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> UpdateTmParsNNTabVars\n");
    VALIDATE(tmPO,TM_PARS_ID);
    tmPO->ht_col = 0;
    tmPO->st_col = 0;
    /***
    *   Number of columns and find any default columns
    */
    if(tmPO->htab) {
        DB_DNA_TM DB_PrI("+ H table\n");
        tmPO->ht_col = GetTableColsI(tmPO->htab,FALSE);
        SetTmParsTableColLabs(tmPO, H_TABLE);
/*
        n = FindTableNamedColI(tmPO->htab,-1,"D",FALSE,1);
        tmPO->ht_d = (n >= 0) ? n : 0;
*/
    } 
    if(tmPO->stab) {
        DB_DNA_TM DB_PrI("+ S table\n");
        tmPO->st_col = GetTableColsI(tmPO->stab,FALSE);
        SetTmParsTableColLabs(tmPO, S_TABLE);
    } 
    DB_DNA_TM {
        DB_DNA_TM DB_PrI("+ H col=%d d=%d\n",tmPO->ht_col,tmPO->ht_d);
        DB_DNA_TM DB_PrI("+ S col=%d d=%d\n",tmPO->st_col,tmPO->st_d);
        DB_DNA_TM DB_PrI("<< UpdateTmParsNNTabVars\n");
    }
}
/***************************************************************************
*   Set pretty formating for tables ...
*/
void PrettyTmParTableForm(TM_PARS *tmPO)
{
    char nameS[NSIZE];

    VALIDATE(tmPO,TM_PARS_ID);
    if(tmPO->htab) {
        VALIDATE(tmPO->htab,TABLE_ID);
        sprintf(nameS,"Enthalpy dH");
        SetTableNamesI(tmPO->htab,nameS,NULL,strlen(nameS));
        SetTablePrintformI(tmPO->htab,"%9.2f",NULL,"\t",NULL,"#\t");
    }   
    if(tmPO->stab) {
        VALIDATE(tmPO->stab,TABLE_ID);
        sprintf(nameS,"Entropy dS");
        SetTableNamesI(tmPO->stab,nameS,NULL,strlen(nameS));
        SetTablePrintformI(tmPO->stab,"%9.2f",NULL,"\t",NULL,"#\t");
    }
}
/***************************************************************************
*   Explicitly set enthalpy and entropy for SantaLucia's algorithm
*/
void SetSantaTmThermPars(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> SetSantaTmThermPars\n");
    VALIDATE(tmPO,TM_PARS_ID);
    InitTmParNNTablesI(tmPO, NN_IND_NUM, 1);
    /***
    *   Set H and S values
    */
    SetThermoTableNNValI(tmPO, H_TABLE, "AA", 0, -7900);
    SetThermoTableNNValI(tmPO, H_TABLE, "AC", 0, -8400);
    SetThermoTableNNValI(tmPO, H_TABLE, "AG", 0, -7800);
    SetThermoTableNNValI(tmPO, H_TABLE, "AT", 0, -7200);
    SetThermoTableNNValI(tmPO, H_TABLE, "CA", 0, -8500);
    SetThermoTableNNValI(tmPO, H_TABLE, "CC", 0, -8000);
    SetThermoTableNNValI(tmPO, H_TABLE, "CG", 0, -10600);
    SetThermoTableNNValI(tmPO, H_TABLE, "CT", 0, -7800);
    SetThermoTableNNValI(tmPO, H_TABLE, "GA", 0, -8200);
    SetThermoTableNNValI(tmPO, H_TABLE, "GC", 0, -9800);
    SetThermoTableNNValI(tmPO, H_TABLE, "GG", 0, -8000);
    SetThermoTableNNValI(tmPO, H_TABLE, "GT", 0, -8400);
    SetThermoTableNNValI(tmPO, H_TABLE, "TA", 0, -7200);
    SetThermoTableNNValI(tmPO, H_TABLE, "TC", 0, -8200);
    SetThermoTableNNValI(tmPO, H_TABLE, "TG", 0, -8500);
    SetThermoTableNNValI(tmPO, H_TABLE, "TT", 0, -7900);
    SetThermoTableNNValI(tmPO, H_TABLE, "InitAT", 0, 2300);
    SetThermoTableNNValI(tmPO, H_TABLE, "InitCG", 0, 100);
    SetThermoTableNNValI(tmPO, S_TABLE, "AA", 0, -22.2);
    SetThermoTableNNValI(tmPO, S_TABLE, "AC", 0, -22.4);
    SetThermoTableNNValI(tmPO, S_TABLE, "AG", 0, -21.0);
    SetThermoTableNNValI(tmPO, S_TABLE, "AT", 0, -20.4);
    SetThermoTableNNValI(tmPO, S_TABLE, "CA", 0, -22.7);
    SetThermoTableNNValI(tmPO, S_TABLE, "CC", 0, -19.9);
    SetThermoTableNNValI(tmPO, S_TABLE, "CG", 0, -27.2);
    SetThermoTableNNValI(tmPO, S_TABLE, "CT", 0, -21.0);
    SetThermoTableNNValI(tmPO, S_TABLE, "GA", 0, -22.2);
    SetThermoTableNNValI(tmPO, S_TABLE, "GC", 0, -24.4);
    SetThermoTableNNValI(tmPO, S_TABLE, "GG", 0, -19.9);
    SetThermoTableNNValI(tmPO, S_TABLE, "GT", 0, -22.4);
    SetThermoTableNNValI(tmPO, S_TABLE, "TA", 0, -21.3);
    SetThermoTableNNValI(tmPO, S_TABLE, "TC", 0, -22.2);
    SetThermoTableNNValI(tmPO, S_TABLE, "TG", 0, -22.7);
    SetThermoTableNNValI(tmPO, S_TABLE, "TT", 0, -22.2);
    SetThermoTableNNValI(tmPO, S_TABLE, "InitAT", 0, 4.1);
    SetThermoTableNNValI(tmPO, S_TABLE, "InitCG", 0, -2.8);
    sprintf(tmPO->source,"SantaLucia'98 (PNAS)"); 
    DB_DNA_TM DB_PrI("<< SetSantaTmThermPars\n");
}
/***************************************************************************
*   Explictly set enthalpy and entropy for Breslaur's ("Oligo") algorithm 
*/
void SetOligoTmThermPars(TM_PARS *tmPO)
{
    DB_DNA_TM DB_PrI(">> SetOligoTmThermPars\n");
    VALIDATE(tmPO,TM_PARS_ID);
    InitTmParNNTablesI(tmPO, NN_IND_NUM, 1);
    /***
    *   Set H and S values
    */
    SetThermoTableNNValI(tmPO, H_TABLE, "AA", 0, -9100);
    SetThermoTableNNValI(tmPO, H_TABLE, "AC", 0, -6500);
    SetThermoTableNNValI(tmPO, H_TABLE, "AG", 0, -7800);
    SetThermoTableNNValI(tmPO, H_TABLE, "AT", 0, -8600);
    SetThermoTableNNValI(tmPO, H_TABLE, "CA", 0, -5800);
    SetThermoTableNNValI(tmPO, H_TABLE, "CC", 0, -11000);
    SetThermoTableNNValI(tmPO, H_TABLE, "CG", 0, -11900);
    SetThermoTableNNValI(tmPO, H_TABLE, "CT", 0, -7800);
    SetThermoTableNNValI(tmPO, H_TABLE, "GA", 0, -5600);
    SetThermoTableNNValI(tmPO, H_TABLE, "GC", 0, -11100);
    SetThermoTableNNValI(tmPO, H_TABLE, "GG", 0, -11000);
    SetThermoTableNNValI(tmPO, H_TABLE, "GT", 0, -6500);
    SetThermoTableNNValI(tmPO, H_TABLE, "TA", 0, -6000);
    SetThermoTableNNValI(tmPO, H_TABLE, "TC", 0, -5600);
    SetThermoTableNNValI(tmPO, H_TABLE, "TG", 0, -5800);
    SetThermoTableNNValI(tmPO, H_TABLE, "TT", 0, -9100);
    SetThermoTableNNValI(tmPO, H_TABLE, "InitAT", 0, 0);
    SetThermoTableNNValI(tmPO, H_TABLE, "InitCG", 0, 0);
    SetThermoTableNNValI(tmPO, S_TABLE, "AA", 0, -24.0);
    SetThermoTableNNValI(tmPO, S_TABLE, "AC", 0, -17.3);
    SetThermoTableNNValI(tmPO, S_TABLE, "AG", 0, -20.8);
    SetThermoTableNNValI(tmPO, S_TABLE, "AT", 0, -23.9);
    SetThermoTableNNValI(tmPO, S_TABLE, "CA", 0, -12.9);
    SetThermoTableNNValI(tmPO, S_TABLE, "CC", 0, -26.6);
    SetThermoTableNNValI(tmPO, S_TABLE, "CG", 0, -27.8);
    SetThermoTableNNValI(tmPO, S_TABLE, "CT", 0, -20.8);
    SetThermoTableNNValI(tmPO, S_TABLE, "GA", 0, -13.5);
    SetThermoTableNNValI(tmPO, S_TABLE, "GC", 0, -26.7);
    SetThermoTableNNValI(tmPO, S_TABLE, "GG", 0, -26.6);
    SetThermoTableNNValI(tmPO, S_TABLE, "GT", 0, -17.3);
    SetThermoTableNNValI(tmPO, S_TABLE, "TA", 0, -16.9);
    SetThermoTableNNValI(tmPO, S_TABLE, "TC", 0, -13.5);
    SetThermoTableNNValI(tmPO, S_TABLE, "TG", 0, -12.9);
    SetThermoTableNNValI(tmPO, S_TABLE, "TT", 0, -24.0);
    SetThermoTableNNValI(tmPO, S_TABLE, "InitAT", 0, -15.1);
    SetThermoTableNNValI(tmPO, S_TABLE, "InitCG", 0, -15.1);
    sprintf(tmPO->source,"Oligo / Breslauer'86 (PNAS)"); 
    DB_DNA_TM DB_PrI("<< SetOligoTmThermPars\n");
}
/***************************************************************************
*   Set specific NN thermo value in specific table 
*/
int SetThermoTableNNValI(TM_PARS *tmPO, int which, char *nnS, int c, DOUB vD)
{
    int r;
    TABLE *tabPO;

    /***
    *   Get correct table
    */
    VALIDATE(tmPO,TM_PARS_ID);
    if(which == H_TABLE) {
        tabPO = tmPO->htab;
    }
    else if(which == S_TABLE) {
        tabPO = tmPO->stab;
    }
    else {
        printf("Table=%d=?\n",which);
        ERR("SetThermoTableNNValI","Bogus table code");
        return(FALSE);
    }
    VALIDATE(tabPO,TABLE_ID);
    /***
    *   Index and also set lable too
    */
    r = NNSeqPairIndexI(nnS);
    if(r<0) {
        return(FALSE);
    }
    c = (c<0) ? 0 : c;
    SetTableValI(tabPO,r,c,vD);
    SetTableRowLabI(tabPO,r,nnS);
    return(TRUE);
}
/***************************************************************************
*   NN index = NN + Init terms
*/
int NNSeqPairIndexI(char *inS)
{
    int i;
    char nnS[NSIZE];

    i = BOGUS;
    strcpy(nnS,inS);
    Upperize(nnS);
    if( EQSTRING(nnS,"INITAT",6) || EQSTRING(nnS,"INAT",4) ) {  
        i = INITAT_IND;  
    }
    else if( EQSTRING(nnS,"INITCG",6) || EQSTRING(nnS,"INGC",4) ) {     
        i = INITCG_IND;  
    }
    else {  
        i = SeqPairIndexI(nnS); 
    }
    return(i);
}
/***************************************************************************
*   Set algorithm flag and optionally appropriate parameters in data structure
*/
void SetTmParAlgo(TM_PARS *tmPO, int which, int init)
{
    DB_DNA_TM DB_PrI(">> SetTmParAlgo which=%d init=%d\n",which,init);
    VALIDATE(tmPO,TM_PARS_ID);
    switch(which)
    {
        case ALGO_24:
        case ALGO_PNA:
            break;
        case ALGO_MELTING:
        case ALGO_OLIGO:
            if(init) {
                SetOligoTmThermPars(tmPO); 
            }
            break;
        case ALGO_PEYRET:
            if(init) {
                InitTmCompleteParsI(tmPO);
            }
            break;
        case ALGO_SANTA:
            if(init) {
                SetSantaTmThermPars(tmPO);
            }
            break;
        default:
            printf("Bogus algorithm code=%d\n",tmPO->algo);
            ERR("SetTmParAlgo","Bogus alg code");
    }
    tmPO->algo = which;     
    DB_DNA_TM DB_PrI("<< SetTmParAlgo \n");
}
/*************************************************************************/
int GetTmParAlgoI(TM_PARS *tmPO)
{
    VALIDATE(tmPO,TM_PARS_ID);
    return(tmPO->algo);
}
/***************************************************************************/
int LegitTmParConcI(DOUB cD)
{
    if((cD<MIN_NN_CONC)||(cD>MAX_NN_CONC))
    {   return(FALSE);  }
    return(TRUE);
}
/***************************************************************************
*   Set DNA concentration value(s) in global structure
*/
void SetTmParConc(TM_PARS *tmPO,DOUB cD, DOUB secD, DOUB thirD)
{
    DB_DNA_TM DB_PrI(">> SetTmParConc %e %e %e\n",cD,secD,thirD);
    VALIDATE(tmPO,TM_PARS_ID);
    if(!LegitTmParConcI(cD)) {
        WARNLINE;
        printf("# Allowed concentration range: %e to %e M\n",
            MIN_NN_CONC,MAX_NN_CONC);
        printf("# Setting DNA concentration to %e not %e\n",DEF_NN_CONC,cD);
        printf("\n");
        tmPO->conc = DEF_NN_CONC;
        printf("\n");
    }
    else {
        tmPO->conc = cD;
    }
    /***
    *   Second concentration:
    *   If real, set it, otherwise set it to the first concentration
    */
    if(!BAD_REAL(secD)) {
        if(!LegitTmParConcI(secD)) {
            WARNLINE;
            printf("# Allowed concentration range: %e to %e M\n",
                MIN_NN_CONC,MAX_NN_CONC);
            printf("# Setting 2'nd DNA concentration to %e not %e\n",
                DEF_NN_CONC,cD);
            printf("\n");
            tmPO->conc2 = DEF_NN_CONC;
        }
        else {
            tmPO->conc2 = secD;
        }
    }
    else {
        tmPO->conc2 = tmPO->conc;
    }
    /***
    *   Third concentration:
    *   If real, set it, otherwise set it to the second concentration
    */
    if(!BAD_REAL(thirD)) {
        if(!LegitTmParConcI(thirD)) {
            WARNLINE;
            printf("# Allowed concentration range: %e to %e M\n",
                MIN_NN_CONC,MAX_NN_CONC);
            printf("# Setting 3'nd DNA concentration to %e not %e\n",
                DEF_NN_CONC,cD);
            printf("\n");
            tmPO->conc3 = DEF_NN_CONC;
        }
        else {
            tmPO->conc3 = thirD;
        }
    }
    else {
        tmPO->conc3 = tmPO->conc;
    }
    DB_DNA_TM DB_PrI("<< SetTmParConc\n");
}
/***************************************************************************/
int GetTmParConcI(TM_PARS *tmPO,DOUB *fPD, DOUB *sPD, DOUB *tPD)
{
    VALIDATE(tmPO,TM_PARS_ID);
    if(fPD) {
        *fPD = tmPO->conc;
    }
    if(sPD) {
        *sPD = tmPO->conc2;
    }
    if(tPD) {
        *tPD = tmPO->conc3;
    }
    return(TRUE);
}
/***************************************************************************/
int ConcTermsForTmAlgoI(TM_PARS *tmPO)
{
    int n;

    VALIDATE(tmPO,TM_PARS_ID);
    n = BOGUS;
    switch(tmPO->algo)
    {
        case ALGO_24:       
            n = 0;  break;
        case ALGO_PEYRET:   
            n = 3;  break;
        case ALGO_SANTA:    
            n = 2;  break;
        default:
            n = 1;  break;
    }
    return(n);
}
/***************************************************************************/
int LegitTmParSaltI(DOUB cD)
{
    if((cD<MIN_NN_SALT)||(cD>MAX_NN_SALT))
    {   return(FALSE);  }
    return(TRUE);
}
/***************************************************************************
*   Set salt value in global structure
*/
void SetTmParSalt(TM_PARS *tmPO,DOUB cD)
{
    DB_DNA_TM DB_PrI(">> SetTmParSalt %f\n",cD);
    VALIDATE(tmPO,TM_PARS_ID);
    if(!LegitTmParSaltI(cD)) {
        WARNLINE;
        printf("# Allowed concentration range: %e to %e M\n",MIN_NN_SALT,
            MAX_NN_SALT);
        printf("# Setting salt concentration to %e not %e\n",DEF_NN_SALT,cD);
        printf("\n");
        tmPO->salt = DEF_NN_SALT;
    }
    else {
        tmPO->salt = cD;
    }
    DB_DNA_TM DB_PrI("<< SetTmParSalt\n");
}
/***************************************************************************/
int GetTmParSaltI(TM_PARS *tmPO,DOUB *saltPD)
{
    VALIDATE(tmPO,TM_PARS_ID);
    if(saltPD) {
        *saltPD = tmPO->salt;
    }
    return(TRUE);
}
/***************************************************************************/
int LegitTmParTempI(DOUB tD)
{
    if( (tD<MIN_NN_TEMP) || (tD>MAX_NN_TEMP) )
    {   return(FALSE);  }
    return(TRUE);
}
/***************************************************************************/
void SetTmParTemp(TM_PARS *tmPO,DOUB cD)
{
    DB_DNA_TM DB_PrI(">> SetTmParTemp %f\n",cD);
    VALIDATE(tmPO,TM_PARS_ID);
    if(!LegitTmParTempI(cD)) {
        WARNLINE;
        printf("# Setting temperature to %e not %e\n",DEF_NN_TEMP,cD);
        printf("\n");
        tmPO->tp = DEF_NN_TEMP;
    }
    else {
        tmPO->tp = cD;
    }
    DB_DNA_TM DB_PrI("<< SetTmParTemp\n");
}
/***************************************************************************/
int GetTmParTempI(TM_PARS *tmPO,DOUB *tPD)
{
    VALIDATE(tmPO,TM_PARS_ID);
    if(tPD) {
        *tPD = tmPO->tp;
    }
    return(TRUE);
}
/***************************************************************************/
int LegitTmParMgI(DOUB cD)
{
    if((cD<MIN_NN_MG)||(cD>MAX_NN_MG))
    {   return(FALSE);  }
    return(TRUE);
}
/************************************************************************/
void SetTmParMg(TM_PARS *tmPO,DOUB cD)
{
    DB_DNA_TM DB_PrI(">> SetTmParMg %f\n",cD);
    VALIDATE(tmPO,TM_PARS_ID);
    if(!LegitTmParMgI(cD)) {
        WARNLINE;
        printf("# Allowed concentration range: %e to %e M\n",MIN_NN_MG,
            MAX_NN_MG);
        printf(" Setting Mg concentration to %e not %e\n",DEF_NN_MG,cD);
        printf("\n");
        tmPO->mg = DEF_NN_MG;
    }
    else {
        tmPO->mg = cD;
    }
    DB_DNA_TM DB_PrI("<< SetTmParMg\n");
}
/************************************************************************/
int GetTmParMgI(TM_PARS *tmPO,DOUB *tPD)
{
    VALIDATE(tmPO,TM_PARS_ID);
    if(tPD) {
        *tPD = tmPO->mg;
    }
    return(TRUE);
}
/***************************************************************************
*   Set base to substitute for ambigs in global structure
*/
void SetTmParAmbigBaseOK(TM_PARS *tmPO,int ok)
{
    VALIDATE(tmPO,TM_PARS_ID);
    tmPO->ambig_ok = ok;
}
/****************************************************************************/
void SetTmStrandDir(TM_PARS *tmPO,int dir)
{
    VALIDATE(tmPO,TM_PARS_ID);
    tmPO->dir = dir;
}
/***************************************************************************
*   Set algorithm flag in data structure
*/
void SetTmParDefCol(TM_PARS *tmPO,int which, int col)
{
    int h,s;

    DB_DNA_TM DB_PrI(">> SetTmParDefCol which=%d col=%d\n",which,col);
    VALIDATE(tmPO,TM_PARS_ID);
    UpdateTmParsNNTabVars(tmPO);
    h = s = 0;
    switch(which) 
    {
        case H_TABLE:   h++;    break;
        case S_TABLE:   s++;    break;
        default:
            if(which < 0) {
                h++;    
                s++;
                break;
            }
            printf("Table=%d=?\n",which);
            ERR("SetTmParDefCol","Bogus table code");
            return;
    }
    col -= 1;
    if(h) {
        LIMIT_NUM(col,0,tmPO->ht_col-1);
        tmPO->ht_d = col;
        SetTmParsTableColLabs(tmPO, H_TABLE);
    }
    if(s) {
        LIMIT_NUM(col,0,tmPO->st_col-1);
        tmPO->st_d = col;
        SetTmParsTableColLabs(tmPO, S_TABLE);
    }
    DB_DNA_TM DB_PrI("<< SetTmParDefCol\n");
}
/***************************************************************************
*   Set table col labels for base posistion...
*/
int SetTmParsTableColLabs(TM_PARS *tmPO, int which)
{
    int i,j,n,c;
    char labS[DEF_BS];
    TABLE *tabPO;

    VALIDATE(tmPO,TM_PARS_ID);
    /***
    *   Get correct table
    */
    if(which == H_TABLE) {
        tabPO = tmPO->htab;
        n = tmPO->ht_col;
        c = tmPO->ht_d;
    }
    else if(which == S_TABLE) {
        tabPO = tmPO->stab;
        n = tmPO->st_col;
        c = tmPO->st_d;
    }
    else {
        printf("Table=%d=?\n",which);
        ERR("SetTmParsTableColLabs","Bogus table code");
        return(FALSE);
    }
    VALIDATE(tabPO,TABLE_ID);
    /***
    *   If just one, no lable and bale
    */
    if(n<2) {
        HideTableColLabsI(tabPO);
        return(TRUE);
    }
    /***
    *   Index and also set lable too
    */
    SetTableColLabI(tabPO,c,"Def"); 
    for(i=0;i<c;i++)
    {
        sprintf(labS,"5'-%d",i+1);
        SetTableColLabI(tabPO,i,labS); 
    }
    j = 1;
    for(i=n-1;i>c;i--)
    {
        sprintf(labS,"3'-%d",j++);
        SetTableColLabI(tabPO,i,labS); 
    }
    return(TRUE);
}
