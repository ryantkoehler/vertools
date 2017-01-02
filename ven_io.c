/*
* ven_io.c
*
* Copyright 2017 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <ctype.h>
#include "prim.h"
#include "venpipe.h"

#define DB_VENIO if(DB[142])

/*************************************************************************
*   Write the whip
*/
void WriteVenpipeHeader(VENPIPE *vPO, char *nameS, FILE *oPF)
{
    char bufS[DEF_BS];

    HAND_NFILE(oPF);
    if(!NO_S(nameS)) {
        fprintf(oPF,"# %s\n",nameS);
    }
    if(vPO->do_pfe) {
        sprintf(bufS,"MFE\tPFE");
    }
    else {
        sprintf(bufS,"MFE");
    }
    /***
    *   Main sham
    */
    VersionSplash(oPF,VEN_VERSION_S,"#  ",TRUE);
    fprintf(oPF,"# Input     %s\n",vPO->inname);
    fprintf(oPF,"# Reporting %s\n",bufS);
    fprintf(oPF,"#  MFE, Miniumum Free Energy\n");
    if(vPO->do_pfe) {
        fprintf(oPF,"#  PFE, Partition Free Energy\n");
    }
    if(vPO->do_dmb) {
        fprintf(oPF,"# Matching base counts\n");
        if(!BAD_INT(vPO->dmb_f)) {
            fprintf(oPF,"#  Positions %d to %d considered",vPO->dmb_f,vPO->dmb_l);
            if(vPO->mrre) {
                fprintf(oPF," (from end)\n");
            }
            fprintf(oPF,"\n");
        }
        else {
            fprintf(oPF,"#  All positions considered\n");
        }
    }
    /***
    *   Vienna settings 
    */
    DumpVenergy(vPO->ven, oPF);
    /***
    *   Columns 
    */
    fprintf(oPF,"# Name\tMFE");
    if(vPO->do_pfe) {
        fprintf(oPF,"\tPFE");
    }
    if(vPO->do_dmb) {
        fprintf(oPF,"\tLen\tHyb\tOpen");
        if(vPO->do_pfe) {
            fprintf(oPF,"\tpHySt\tpHyWk\tpOpen\tpHyb");
        }
    }
    fprintf(oPF,"\n");
}
/***********************************************************************/
void DumpVenergy(VENERGY *vePO, FILE *oPF)
{
    VALIDATE(vePO,VENERGY_ID);
    HAND_NFILE(oPF);
    fprintf(oPF,"#\n");
    fprintf(oPF,"################## Vienna Settings ##################\n");
    fprintf(oPF,"# Version        %s\n",VLIB_VERSION_S);
    if(NO_S(vePO->vparfile)) {
        fprintf(oPF,"# ViennaParams   none = defaults\n");
    }
    else {
        fprintf(oPF,"# ViennaParams   %s\n",vePO->vparfile);
        if(vePO->do_saltcorrect) {
            fprintf(oPF,"# CorrectedSalt  %4.2e M\n",vePO->salt);
        }
        else {
            fprintf(oPF,"# CorrectedSalt  FALSE\n");
        }
    }
    fprintf(oPF,"# ViennaTemp     %2.1f C\n",vePO->temp);
/*
    fprintf(oPF,"# ViennaBackTrac %d\n",vePO->do_backtrack);
    fprintf(oPF,"# noLonelyPairs  %d\n",vePO->nolonelypairs);
*/
    fprintf(oPF,"#\n");
}
/***********************************************************************/
void HandleVenpipeOut(VENPIPE *venPO, FILE *outPF)
{
    int len, ssc_array[MAX_VSLEN], pssc_array[MAX_VSLEN];
    char nameS[DEF_BS],enS[DEF_BS], matS[DEF_BS];
    char *ssPC, *pssPC, *seqPC, *fseqPC;

    HAND_NFILE(outPF);
    strcpy(nameS,venPO->tname);
    seqPC = venPO->tseq;
    ssPC = venPO->tss;
    pssPC = venPO->tss2;
    len = venPO->tlen;
    /***
    *   Energy story
    */
    if(venPO->do_pfe) {
        sprintf(enS,"%7.3f\t%7.3f", venPO->ten, venPO->ten2);
    }
    else {
        sprintf(enS,"%7.3f",venPO->ten);
    }
    /***
    *   Matching base story
    */
    INIT_S(matS);
    if(venPO->do_dmb) {
        StrucStringToArrayI(ssPC, len, ssc_array);
        if(venPO->do_pfe) {
            StrucStringToArrayI(pssPC, len, pssc_array);
        }
        FillSeqStructMatchOutI(venPO, ssc_array, pssc_array ,len, matS);
    }
    /***
    *   What to dump?
    */
    GetSeqSeqI(venPO->seq,&fseqPC);
    if( (venPO->do_ddb) || (venPO->do_mbtab) ) {
        fprintf(outPF,"# %-15s\t%s",nameS,enS);
        if(venPO->do_dmb) {
            fprintf(outPF,"\t%s",matS);
        }
        if(venPO->do_ds) {
            fprintf(outPF,"\t%s",fseqPC);
        }
        fprintf(outPF,"\n");
        fprintf(outPF,"# Sequence  %s\n",seqPC);
        if(venPO->do_ddb) {
            fprintf(outPF,"  MFE_Strc  %s\n",ssPC);
            if(venPO->do_pfe) {
                fprintf(outPF,"  PFE_Strc  %s\n",pssPC);
            }
            fprintf(outPF,"\n");
        }
        else {
            if(venPO->do_pfe) {
                DumpMatchBaseTablesI(venPO,nameS,len,ssc_array,pssc_array,outPF);
            }
            else {
                DumpMatchBaseTablesI(venPO,nameS,len,ssc_array,NULL,outPF);
            }
        }
    }
    else if(venPO->ofas) {
        fprintf(outPF,">%s  %s",nameS,enS);
        if(venPO->do_dmb) {
            fprintf(outPF,"\t%s",matS);
        }
        fprintf(outPF,"\n");
        fprintf(outPF,"%-15s\n",seqPC);
    }
    else if(venPO->oraw) {
        fprintf(outPF,"# %-15s\t%s",nameS,enS);
        if(venPO->do_dmb) {
            fprintf(outPF,"\t%s",matS);
        }
        fprintf(outPF,"\n");
        fprintf(outPF,RAW_PFORM_S,nameS,seqPC);
    }
    else {
        fprintf(outPF,"%-15s\t%s",nameS,enS);
        if(venPO->do_dmb) {
            fprintf(outPF,"\t%s",matS);
        }
        if(venPO->do_ds) {
            fprintf(outPF,"\t%s",fseqPC);
        }
        fprintf(outPF,"\n");
    }
    return;
}
/***********************************************************************/
int DumpMatchBaseTablesI(VENPIPE *venPO, char *nameS, int len, int *sscPI, 
    int *psscPI, FILE *outPF)
{
    char bufS[NSIZE];
    int ok;

    if(psscPI) {
        sprintf(bufS,"%s__MFE-tab ",nameS);
        ok = DumpOneMatchBaseTabI(venPO, bufS, len, sscPI, outPF);
        sprintf(bufS,"%s__PFE-tab ",nameS);
        ok = DumpOneMatchBaseTabI(venPO, bufS, len, psscPI, outPF);
    }
    else {
        sprintf(bufS,"%-15s",nameS);
        ok = DumpOneMatchBaseTabI(venPO, bufS, len, sscPI, outPF);
    }
    return(ok);
}
/***********************************************************************/
int DumpOneMatchBaseTabI(VENPIPE *venPO, char *nameS, int len, int *sscPI, 
    FILE *outPF)
{
    int i,first,last;

    if(!GetMatchBaseSeqCoordsI(venPO, len, &first, &last, NULL)) {
        return(FALSE);
    }
    HAND_NFILE(outPF);
    fprintf(outPF,"%s",nameS);
    for(i=first;i<last;i++) 
    {
        fprintf(outPF," %d",sscPI[i]);
    }
    fprintf(outPF,"\n");
    return(TRUE);
}
/***********************************************************************
*   Figure out start / end and range for given length based on settings
*/
int GetMatchBaseSeqCoordsI(VENPIPE *vpPO, int len, int *stPI, int *enPI, int *rPI)
{
    int first,last,range;

    first = 0;
    last = len;
    if(!BAD_INT(vpPO->dmb_f)) {
        if(vpPO->mrre) {
            first = len - vpPO->dmb_l;
            last = len - vpPO->dmb_f + 1;
        }
        else {
            first = vpPO->dmb_f - 1;
            last = vpPO->dmb_l;
        }
    }
    range = last - first;
    if(range < 1) {
        return(FALSE);
    }
    if(stPI) {
        *stPI = first;
    }
    if(enPI) {
        *enPI = last;
    }
    if(rPI) {
        *rPI = range;
    }
    return(TRUE);
}
/***********************************************************************
*   Fill string with formatted numbers for match base counts
*/
int FillSeqStructMatchOutI(VENPIPE *vpPO, int *sscPI, int *psscPI, int len, char *matS)
{
    int first,last,range,n1,n2;
    char bufS[DEF_BS];

    if(!GetMatchBaseSeqCoordsI(vpPO, len, &first, &last, &range)) {
        return(FALSE);
    }
    /***
    *   xxx Sham? Quietly bound without any warning ???
    */
    LIMIT_NUM(first,0,len);
    LIMIT_NUM(last,first,len);
    range = last - first;
    /***
    *   Get count(s) and build string to output
    */
    n1 = NumArrayValsI(sscPI, IS_INT, first, last, 1.0, 1000.0);
    sprintf(matS,"%3d\t%3d\t%3d",range,n1,range-n1);
    if(vpPO->do_pfe) {
        n1 = NumArrayValsI(psscPI, IS_INT, first, last, 1.0, 1.0);
        n2 = NumArrayValsI(psscPI, IS_INT, first, last, 2.0, 1000.0);
        sprintf(bufS,"\t%3d\t%3d\t%3d\t%3d",n2, n1, range-(n1+n2), n1+n2);
        strcat(matS,bufS);
    }
    return(TRUE);
}
/***********************************************************************
*   Take a "structure" string and convert this to numbers
*/
int StrucStringToArrayI(char *seqS, int len, int *aPI)
{
    int i,v;

    for(i=0;i<len;i++)
    {
        switch(seqS[i])
        {
            case '.':    v=0;     break;
            case ',':    v=1;     break;
            case '{':    v=1;     break;
            case '}':    v=1;     break;
            case '|':    v=2;     break;
            case '(':    v=2;     break;
            case ')':    v=2;     break;
            default:
                return(FALSE);
        }
        aPI[i] = v;
    }
    return(TRUE);
}
