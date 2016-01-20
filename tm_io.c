/*
* tm_io.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prim.h"
#include "table.h"
#include "tm_pars.h"

#define DB_DNA_TM   if(DB[110])

/*************************************************************************
*   Dump out TmPars datastructure
*/
void DumpTmPars(TM_PARS *tmPO, char *preS, int full, FILE *outPF)
{
    char bufS[DEF_BS],upS[DEF_BS];

    VALIDATE(tmPO,TM_PARS_ID);
    HAND_NFILE(outPF);  
    if(preS) {
        sprintf(upS,"%s",preS);
    }
    else {
        INIT_S(upS);
    }
    fprintf(outPF,"######################################################\n");
    fprintf(outPF,"#    Parameters for DNA Tm / Thermo calculations\n");
    fprintf(outPF,"######################################################\n");
    FillTmAlgoString(tmPO,bufS);
    fprintf(outPF,"%sAlgorithm:  %s\n",upS,bufS);
    FillTmConcString(tmPO,bufS);
    fprintf(outPF,"%sConc.       %s\n",upS,bufS);
    FillTmSaltString(tmPO,bufS);
    fprintf(outPF,"%sSalt        %s\n",upS,bufS);
    FillTmMgString(tmPO,bufS);
    fprintf(outPF,"%sMg          %s\n",upS,bufS);
    FillTmTempString(tmPO,bufS);
    fprintf(outPF,"%sTemp        %s\n",upS,bufS);
    fprintf(outPF,"%s\n",upS);
    FillTmParSourceString(tmPO,bufS);
    fprintf(outPF,"%sParameters: %s\n",upS,bufS);
    if(full) {
        switch(tmPO->algo)
        {
            case ALGO_PEYRET:
                DumpCompleteTmPars(tmPO,outPF); 
                break;
            case ALGO_24:
                break;
            default:
                DumpSimpleTmPars(tmPO,preS,outPF); 
                break;
        }
    }
}
/*************************************************************************
*   Dump simple NN params  (i.e. only NN + inits)
*/
void DumpSimpleTmPars(TM_PARS *tmPO, char *preS,FILE *outPF)
{
    VALIDATE(tmPO,TM_PARS_ID);
    HAND_NFILE(outPF);
    DumpTable(tmPO->htab,FALSE,FALSE,outPF);
    fprintf(outPF,"#\n");
    DumpTable(tmPO->stab,FALSE,FALSE,outPF);
    fprintf(outPF,"#\n");
}
/**************************************************************************
*   Fills passed string with name of Tm algorithm specified in tmPO
*/
void FillTmAlgoString(TM_PARS *tmPO,char *bufS)
{
    VALIDATE(tmPO,TM_PARS_ID);
    if(tmPO->algo == ALGO_24) {
        sprintf(bufS,"Simple Tm = 2AT + 4GC");
        return;
    }
    sprintf(bufS,"Nearest neighbor; ");
    switch(tmPO->algo)
    {
        case ALGO_MELTING:  
            strcat(bufS,"as in \"Melting\" program"); 
            break;
        case ALGO_OLIGO:    
            strcat(bufS,"as in \"Oligo\" program"); 
            break; 
        case ALGO_PNA:      
            strcat(bufS,"for PNA"); 
            break;
        case ALGO_SANTA:    
            strcat(bufS,"SantaLucia'98 (with init terms)"); 
            break;
        case ALGO_PEYRET:   
            strcat(bufS,"Peyret (mismatch / terminal)"); 
            break;
        default:
            printf("Bogus algorithm code=%d\n",tmPO->algo);
            ERR("","Bogus alg code");
    }
}
/***************************************************************************/
void FillTmConcString(TM_PARS *tmPO,char *bufS)
{
    int n;
        
    VALIDATE(tmPO,TM_PARS_ID);
    n = ConcTermsForTmAlgoI(tmPO);
    switch(n)
    {
        case 0:         
            sprintf(bufS,"(not used)");     
            break; 
        case 1:     
            sprintf(bufS,"%1.3e (M)",tmPO->conc);
            break; 
        case 2:     
            sprintf(bufS,"%1.3e  %1.3e (M)", tmPO->conc, tmPO->conc2);
            break; 
        case 3:     
            sprintf(bufS,"%1.3e  %1.3e  %1.3e (M)", tmPO->conc, 
                tmPO->conc2, tmPO->conc3);
            break; 
        default:
            sprintf(bufS,"?????");
    }
}
/***************************************************************************/
void FillTmSaltString(TM_PARS *tmPO,char *bufS)
{
    VALIDATE(tmPO,TM_PARS_ID);
    switch(tmPO->algo)
    {
        case ALGO_24:       sprintf(bufS,"(not used)");     break; 
        default:
            sprintf(bufS,"%1.3e (M)",tmPO->salt);
    }
}
/***************************************************************************/
void FillTmTempString(TM_PARS *tmPO,char *bufS)
{
    VALIDATE(tmPO,TM_PARS_ID);
    switch(tmPO->algo)
    {
        case ALGO_24:       
            sprintf(bufS,"(not used)");     
            break; 
        default:
            sprintf(bufS,"%4.1f (C)",tmPO->tp);
    }
}
/***************************************************************************/
void FillTmMgString(TM_PARS *tmPO,char *bufS)
{
    VALIDATE(tmPO,TM_PARS_ID);
    switch(tmPO->algo)
    {
        case ALGO_OLIGO:        
        case ALGO_PEYRET:       
            sprintf(bufS,"%1.3e (M)",tmPO->mg); break; 
        default:
            sprintf(bufS,"(not used)");     
    }
}
/*************************************************************************/
void FillTmStrandDirString(TM_PARS *tmPO,char *bufS)
{
    VALIDATE(tmPO,TM_PARS_ID);
    switch(tmPO->algo)
    {
        default:
            sprintf(bufS,"Not used (symmetric)");
    }
}
/****************************************************************************
*   Parses Tm algorithm by name
*/
int ParseTmAlgoI(char *nameS)
{
    char upS[DEF_BS];
    int algo;

    strcpy(upS,nameS);
    Upperize(upS);
    if(EQSTRING(upS,"DEF",3)) {     
        algo = DEF_TM_ALGO;     
    }
    else if(EQSTRING(upS,"SANTA",5)) {  
        algo = ALGO_SANTA;  
    }
    else if(EQSTRING(upS,"OLI",3)) {    
        algo = ALGO_OLIGO;  
    }
    else if(EQSTRING(upS,"MELT",4)) {   
        algo = ALGO_MELTING;    
    }
    else if(EQSTRING(upS,"24",2)) {     
        algo = ALGO_24;     
    }
    else if(EQSTRING(upS,"PNA",3)) {    
        algo = ALGO_PNA;    
    }
    else if(EQSTRING(upS,"PEY",3)) {    
        algo = ALGO_PEYRET;     
    }
    else {  
        algo = BOGUS;   
    }
    return(algo);
}
/**************************************************************************
*   Fills passed string with name of Tm parameters specified in tmPO
*/
void FillTmParSourceString(TM_PARS *tmPO,char *bufS)
{
    VALIDATE(tmPO,TM_PARS_ID);
    switch(tmPO->algo)
    {
        case ALGO_24:       sprintf(bufS,"(not used)");     break; 
        default:
            sprintf(bufS,"%s",tmPO->source);
    }
}
/***************************************************************************
*   Load nearest-neighbor energy perameters from file. 
*   Returns the TRUE if parameters are loaded; FALSE if problem
*/
int LoadTmParNNEnergiesI(TM_PARS *tmPO,char *inS)
{
    int n,ok;
    char bufS[DEF_BS];
    FILE *inPF;
    TABLE *htabPO, *stabPO;
    
    DB_DNA_TM DB_PrI(">> LoadTmParNNEnergiesI |%s|\n",inS);
    VALIDATE(tmPO,TM_PARS_ID);
    if(!(inPF=OpenUFilePF(inS,"r",NULL)))
    {
        DB_DNA_TM DB_PrI("<< LoadTmParNNEnergiesI FALSE\n");
        return(FALSE);
    }
    /***
    *   Different input files depending on algorithm
    */
    switch(tmPO->algo)
    {
        case ALGO_PEYRET:
            ok = LoadTmCompleteParsI(tmPO,inPF);
            FILECLOSE(inPF);
            return(ok);
        default:
            break;
    }
    /***
    *   Enthalpy / entropy tables
    */
    ok = TRUE;
    htabPO = stabPO = NULL;
    if( (ok) && (!ParseTableI(inPF,TRUE,FALSE,FALSE,FALSE,&htabPO)) ) {
        ok = FALSE;
    }
    if( (ok) && (!ParseTableI(inPF,TRUE,FALSE,FALSE,FALSE,&stabPO)) ) {
        ok = FALSE;
    }
    FILECLOSE(inPF);
    /***
    *   Check that tables have right number of rows
    */
    if(ok) {
        n = GetTableRowsI(htabPO,FALSE);
        if(n != NN_IND_NUM) {
            PROBLINE;
            printf("H Table loaded has %d rows; expect %d\n",n,NN_IND_NUM);
            ok = FALSE;
        }
        n = GetTableRowsI(stabPO,FALSE);
        if(n != NN_IND_NUM) {
            PROBLINE;
            printf("S Table loaded has %d rows; expect %d\n",n,NN_IND_NUM);
            ok = FALSE;
        }
    }
    /***
    *   If good, attach tables...
    */
    if(ok) {
        CHECK_TABLE(tmPO->htab);
        tmPO->htab = htabPO;
        CHECK_TABLE(tmPO->stab);
        tmPO->stab = stabPO;
        PrettyTmParTableForm(tmPO);
        UpdateTmParsNNTabVars(tmPO);
        INIT_S(bufS);
        if( (tmPO->ht_col > 1) || (tmPO->st_col > 1) ) {
            sprintf(bufS," (%d dH cols; %d dS cols)",
                tmPO->ht_col, tmPO->st_col);
        }
        sprintf(tmPO->source,"File %s%s\n",inS,bufS);
    }
    /***
    *   Clean up?
    */
    if(!ok) {
        CHECK_TABLE(htabPO);
        CHECK_TABLE(stabPO);
    }
    return(ok);
}
