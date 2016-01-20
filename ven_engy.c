/*
* ven_engy.c
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
#include "dna.h"
/** Define for global vienna_ok **/
#define __VENERGY__ 
#include "venpipe.h"
/** Vienna headers **/
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"


#define DB_VLIB if(DB[140])


void read_parameter_file(const char fname[]);


/******************************************************************************
*   Dump out vienna global vars
*/
void DumpViennaGlobals(FILE *fPF)
{
    HAND_NFILE(fPF);
    fprintf(fPF,"+  VGV temperature=%2.2f\n",temperature);
    fprintf(fPF,"+  VGV do_backtrack=%d\n",do_backtrack);
    fprintf(fPF,"+  VGV noGU=%d\n",noGU);
    fprintf(fPF,"+  VGV no_closingGU=%d\n",no_closingGU);
    fprintf(fPF,"+  VGV noLonelyPairs=%d\n",noLonelyPairs);
    fprintf(fPF,"+  VGV tetra_loop=%d\n",tetra_loop);
    fprintf(fPF,"+  VGV energy_set=%d\n",energy_set);
    fprintf(fPF,"+  VGV fold_constrained=%d\n",fold_constrained);
    fprintf(fPF,"+  VGV dangles=%d\n",dangles);
    fprintf(fPF,"+  VGV logML=%d\n",logML);
    fprintf(fPF,"+  VGV pf_scale=%f\n",pf_scale);
}
/****************************************************************************
*   Set up Vienna environment variables based on setting in passed structure
*/
int SetUpViennaEnvI(VENERGY *vePO)
{
    VALIDATE(vePO,VENERGY_ID);
    DB_VLIB DB_PrI(">> SetUpViennaEnvI\n");
    /***
    *   Set vienna global vars 
    */
    temperature = vePO->temp;
    noLonelyPairs = 1;
    do_backtrack = 1;       /* eed for partition "structure" string */
    DB_VLIB DumpViennaGlobals(NULL);

    /***
    *   Set flag for intialized
    */
    vienna_okGI = TRUE;
    /***
    *   Parameter file?
    */
    if(!NO_S(vePO->vparfile)) {
        DB_VLIB DB_PrI("+ loading params from |%s|\n",vePO->vparfile);
        if(!SetViennaParametersI(vePO, vePO->vparfile)) {
            vienna_okGI = FALSE;
            DB_VLIB DB_PrI("<< SetUpViennaEnvI FALSE\n");
            return(FALSE);
        }
        DB_VLIB DB_PrI("+ params loaded\n");
    }
    else {
        DB_VLIB DB_PrI("+ no params to load\n");
    }
    DB_VLIB DB_PrI("<< SetUpViennaEnvI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Set Vienna temperature 
*/
int SetViennaTemperatureI(VENERGY *vePO, DOUB tempD)
{
    DB_VLIB DB_PrI(">> SetViennaTemperatureI %4.2lf\n",tempD);
    VALIDATE(vePO,VENERGY_ID);
    if(!vienna_okGI) {
        DB_VLIB DB_PrI("<< SetViennaTemperatureI no init FALSE\n");
        return(FALSE);
    }
    temperature = vePO->temp = tempD;
    DB_VLIB DB_PrI("+ calling update_fold_params\n");
    update_fold_params();
    DB_VLIB DB_PrI("<< SetViennaTemperatureI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Set Vienna salt 
*   7/7/14 RTK; Only set if not 1.0
*/
int SetViennaSaltI(VENERGY *vePO, DOUB saltD, int use, int keep)
{
    DB_VLIB DB_PrI(">> SetViennaSaltI %4.2lf use=%d keep=%d\n",saltD, use, keep);
    VALIDATE(vePO,VENERGY_ID);
    if(!vienna_okGI) {
        DB_VLIB DB_PrI("<< SetViennaSaltI no init FALSE\n");
        return(FALSE);
    }
    vePO->salt = saltD;
    /***
    *   Only do things if not 1.0
    */
    if (saltD != 1.0) {
        vePO->do_saltcorrect = use;
        vePO->do_saltkeep = keep;
    }
    else {
        vePO->do_saltcorrect = FALSE;
        vePO->do_saltkeep = TRUE;
    }
    if( (use) && (!NO_S(vePO->vparfile)) ) {
        DB_VLIB DB_PrI("+ Used and have parameter file so reparse\n");
        if(!SetViennaParametersI(vePO, vePO->vparfile)) {
            DB_VLIB DB_PrI("<< SetViennaSaltI FALSE\n");
            return(FALSE);
        }   
    }
    else {
        DB_VLIB DB_PrI("+ no-reparse: use=%d vpar=|%s|\n",use,vePO->vparfile);
    }
    DB_VLIB DB_PrI("<< SetViennaSaltI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Set Vienna parameters from file
*   Remembers file name, checks that it's real and calles vienna's parsing
*/
int SetViennaParametersI(VENERGY *vePO, char *parS)
{
    FILE *fPF;
    char fnameS[NSIZE],fixS[NSIZE],bufS[NSIZE+DEF_BS];

    DB_VLIB DB_PrI(">> SetViennaParametersI |%s|\n",parS);
    VALIDATE(vePO,VENERGY_ID);
    /***
    *   If empty or NULL string, ignore
    */
    if( (parS==NULL) || (NO_S(parS)) ) {
        DB_VLIB DB_PrI("<< SetViennaParametersI no parfile TRUE\n");
        return(TRUE);
    }
    strcpy(vePO->vparfile,parS);
    if(!vienna_okGI) {
        DB_VLIB DB_PrI("<< SetViennaParametersI no init FALSE\n");
        return(FALSE);
    }
    /***
    *   Check if can open file, and save expanded name (no environment vars)
    */
    if(!(fPF=OpenUFilePF(parS,"r",fnameS))) {
        DB_VLIB DB_PrI("<< SetViennaParametersI bad file FALSE\n");
        return(FALSE);
    }
    CHECK_FILE(fPF);
    /***
    *   If we have to correct salt, get fixed guy
    */
    if( (vePO->salt != 1.0) && (vePO->do_saltcorrect) ) {
        if(!SaltCorrectViennaParsI(fnameS,vePO->salt,fixS)) {
            DB_VLIB DB_PrI("<< SetViennaParametersI failed salt-cor file\n");
            return(FALSE);
        }
    }
    else {
        DB_VLIB DB_PrI("+ no salt correction; will use |%s|\n",fnameS);
        strcpy(fixS,fnameS);
    }
    /***
    *   Call vienna's parsing routine (no useful error handling)
    */
    DB_VLIB DB_PrI("+ loading params from |%s|\n",fixS);
    read_parameter_file(fixS);
    DB_VLIB DB_PrI("+ calling update_fold_params\n");
    update_fold_params();
    /***
    *   Clean up temp file?
    */
    if( (vePO->do_saltcorrect) && (!vePO->do_saltkeep) ) {
        sprintf(bufS,"rm -f %s",fixS);
        if ( system(bufS) ) {
            WARNLINE;
            printf("Failed to remove temp salt file\n|%s|\n",bufS);
        }
    }
    DB_VLIB DB_PrI("<< SetViennaParametersI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Attempt to clean up Vienna environment 
*/
int CleanViennaEnvI()
{
    DB_VLIB DB_PrI(">> CleanViennaEnvI\n");
    if(!vienna_okGI) {
        DB_VLIB DB_PrI("<< CleanViennaEnvI no init FALSE\n");
        return(FALSE);
    }
    vienna_okGI = FALSE;
    DB_VLIB DB_PrI("<< CleanViennaEnvI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Calculate vienna structure and corresponding energy
*   Input is the (DNA) sequence string rseqS[slen]
*   Ouput is a "structure" string put in ssS, if non-null and energy value
*       set to variable enPD
*/
int GetViennaStructI(char *rseqS,int slen, char *ssS, DOUB *enPD,
                        int do_pfold, char *pssS, DOUB *penPD)
{
    DOUB enD, penD;
    char seqS[MAX_VSLEN], structS[MAX_VSLEN], pstructS[MAX_VSLEN];

    DB_VLIB DB_PrI(">> GetViennaStructI slen=%d\n",slen);
    if(!vienna_okGI) {
        DB_VLIB DB_PrI("<< GetViennaStructI no init FALSE\n");
        return(FALSE);
    }
    /***
    *   Too long?
    */
    if(slen>=MAX_VSLEN) {
        printf("Max venpipe seqlen exceeded %d (max=%d)\n",slen,MAX_VSLEN);
        return(FALSE);
    }
    /***
    *   Clean up input seq, make sure uppercase and any T is set to U
    */
    if(!PrepViennaInputSeqI(rseqS,strlen(rseqS),seqS)) {
        DB_VLIB DB_PrI("<< GetViennaStructI seq prep failed FALSE\n");
        return(FALSE);
    }
    enD = penD = TOO_BIG_R;
    if(ssS) {
        INIT_S(ssS);
    }
    if(pssS) {
        INIT_S(pssS);
    }
    INIT_S(structS);
    INIT_S(pstructS);
    /***
    *   Fold, then maybe call partition fold too 
    *   Note: The free_*_arrays() functions are built into fold...
    */
    DB_VLIB DB_PrI("+ Folding\n");
    enD = fold(seqS,structS);
    DB_VLIB DB_PrI("+  |%s|\n",seqS); 
    DB_VLIB DB_PrI("+ f|%s|\n",structS); 
    DB_VLIB DB_PrI("+ fEnergy = %f\n",enD);
    if(do_pfold) {
        SetSeqViennaVarsI(seqS,enD);
        penD = pf_fold(seqS,pstructS);
        DB_VLIB DB_PrI("+ p|%s|\n",pstructS); 
        DB_VLIB DB_PrI("+ pEnergy = %f\n",penD);
    }
    /***
    *   Return what we can
    */
    if(enPD) {
        *enPD = enD;
    }
    if(ssS) {
        sprintf(ssS,"%s",structS);
    }
    if(enPD) {
        *penPD = penD;
    }
    if(pssS) {
        sprintf(pssS,"%s",pstructS);
    }
    DB_VLIB DB_PrI("<< GetViennaStructI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Need to set (global) pf_scale for partition folding
*/
int SetSeqViennaVarsI(char *seqS, DOUB enD) 
{
    DOUB kT;

    /***
    *   Variable for partition folding
    *   Comments and code take from "example.c" example program, 6/24/14 RTK
    *
    * for longer sequences one should also set a scaling factor for
    *   partition function folding, e.g: *
    *   kT = (temperature+273.15)*1.98717/1000.; * kT in kcal/mol *
    *   pf_scale = exp(-e1/kT/strlen(seq1));
    */
    kT = (temperature+273.15)*1.98717/1000.; 
    pf_scale = exp(-enD/kT/strlen(seqS));
/*
printf("# sss enD=%f\ttemp=%f\tkT=%f\tpf_scale=%f\n",enD,temperature,kT,pf_scale);
*/
    
    /***
    *   This should be called after temp change.... or here for each seq
    */
    update_pf_params(strlen(seqS));

    return(TRUE);
}
/***************************************************************************
*   Get the energy cooresponding to structure string ssS and sequence
*       in rseqS.
*/
int GetViennaStructEnergyI(char *rseqS,char *ssS,DOUB *enPD)
{
    DOUB eD;    
    char seqS[MAX_VSLEN];

    DB_VLIB DB_PrI(">> GetViennaStructEnergyI\n");
    if(!vienna_okGI)
    {
        DB_VLIB DB_PrI("<< GetViennaStructEnergyI no init FALSE\n");
        return(FALSE);
    }
    /***
    *   Clean up input seq, make sure uppercase and any T is set to U
    */
    if(!PrepViennaInputSeqI(rseqS,strlen(rseqS),seqS)) {
        DB_VLIB DB_PrI("<< GetViennaStructEnergyI seq preparation failed\n");
        return(FALSE);
    }
    /***
    *   Call to get energy
    */
/* depreciated
    eD = energy_of_struct(seqS,ssS);
*/
    DB_VLIB DB_PrI("+ calling energy_of_structure\n");
    eD = energy_of_structure(seqS,ssS,0);
    *enPD = eD;
    DB_VLIB DB_PrI("<< GetViennaStructEnergyI %f TRUE\n",eD);
    return(TRUE);
}
