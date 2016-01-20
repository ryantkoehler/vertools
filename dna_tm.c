/*
* dna_tm.c
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
#include "table.h"
#include "tm_pars.h"
#include "dset_tm.h"
#include "fbound.h"

#define DB_DNA_TM   if(DB[111])
#define DB_TM_LEN   if(DB[112])

/***************************************************************************
*   Calculates Tm for a sequence seqS[len]
*   Will set Tm, G, H, and S if var pointers are non-NULL
*
*   Returns TRUE if good calculation
*   FALSE if some problem
*/
int SeqTmThermI(TM_PARS *tmPO, char *seqS, int len, DOUB *tPR, DOUB *gPR,
    DOUB *hPR, DOUB *sPR)
{
    int ok;
    DOUB hiD,siD,hD,sD,tmD;
    DOUB tmDA[DSET_EV_NUM], dgDA[DSET_EV_NUM];
    DOUB dhDA[DSET_EV_NUM], dsDA[DSET_EV_NUM];
    char tS[MAX_TM_LEN+1];

    DB_DNA_TM
    { 
        DB_PrI(">> SeqTmThermI algo=%d len=%d\n",tmPO->algo,len);
        PrintString(seqS,len,outfileGPF); DB_PrI("\n");
    }
    /***
    *   Check and initialize passed in args
    */
    VALIDATE(tmPO,TM_PARS_ID);
    if(tPR)
    {   *tPR = BAD_R;   }
    if(gPR)
    {   *gPR = BAD_R;   }
    if(hPR)
    {   *hPR = BAD_R;   }
    if(sPR)
    {   *sPR = BAD_R;   }
    /***
    *   Check length
    */
    if(len>MAX_TM_LEN)
    {
        printf("Too long for Tm: %d (max=%d)\n",len,MAX_TM_LEN);
        ERR("SeqTmThermI","Too long, too long, too long");
        return(FALSE);
    }
    /***
    *   Different plans for different algorithms
    */
    ok = FALSE;
    switch(tmPO->algo)
    {
        case ALGO_24:
            DB_DNA_TM DB_PrI("+ simple composition algorithm; ALGO_24\n");
            tmD = TmFromGCContD(tmPO,seqS,len);
            if(tPR) {
                *tPR = tmD;
            }
            ok = TRUE;
            break;
        /***
        *   PNA 
        */
        case ALGO_PNA:
            tmD = CalcPNATmI(tmPO,seqS,len);
            if(tPR) {
                *tPR = tmD;
            }
            ok = TRUE;
            break;
        /***
        *   PEYRET
        */
        case ALGO_PEYRET:
            /***
            *   Nic Peyret's code call.
            *   This code sets values in passed arrays, so recover after
            *
            *   Here only single sequence, which must be copied to 
            *   it's own null-terminates string
            *
            *   SeqDsetEnergyI(TM_PARS *tmPO, char *seqtopPC,
            *       char *seqbotPC, int verboseI, int modeI, DOUB *tPD, 
            *       DOUB *gPD, DOUB *hPD, DOUB *sP, DOUB *xPD);
            */
            strncpy(tS,seqS,len);
            tS[len]='\0';
            DB_DNA_TM DB_PrI("+ Calling Peyret's Dset\n");
            ok = SeqDsetEnergyI(tmPO, tS, NULL, TRUE, DSET_IMP,
                tmDA, dgDA, dhDA, dsDA, NULL);
            if(!ok) {
                break;
            }
            if(tPR) {
                *tPR = tmDA[0];
            }
            if(gPR) {   
                *gPR = dgDA[0]; 
            }
            if(hPR) {   
                *hPR = dhDA[0]; 
            }
            if(sPR) {   
                *sPR = dsDA[0]; 
            }
            DB_DNA_TM DB_PrI("+ ok=%d Tm=%f dG=%f dH=%f dS=%f\n",ok,tmDA[0],dgDA[0],dhDA[0],dsDA[0]);
            break;
        /***
        *   "Normal" nearesst neighbor case(s)
        *   Locally coded SantaLucia algorithm ...
        *   Or the "Oligo" or "Melting" program variants
        */
        case ALGO_SANTA:
        case ALGO_OLIGO:
        case ALGO_MELTING:
            /***
            *   Get enthalpy and entropy for sequence
            *   Initialization thermo terms then nearest neighbors
            */
            if(!SetSeqTmInitTermsI(tmPO,seqS,len,&hiD,&siD)) {
                return(FALSE);
            }
            if(!TallySeqTmNNTermsI(tmPO,seqS,len,&hD,&sD)) {
                return(FALSE);
            }
            DB_DNA_TM {
                DB_PrI("+ H = %f (Tally %f + Init %f)\n",hD+hiD,hD,hiD);
                DB_PrI("+ S = %f (Tally %f + Init %f)\n",sD+siD,sD,siD);
            }
            hD += hiD;
            sD += siD;
            /***
            *   SantaLucia length-dependent salt correction term
            */
            if(tmPO->algo == ALGO_SANTA) {
                sD += (0.368 * log(tmPO->salt) * RNUM(len));
            }
            /***
            *   Calculate Tm and set other vars if pointers really passed
            */
            if(tPR) {   
                *tPR = TmFromThermoSumsD(tmPO,hD,sD,len);
                DB_DNA_TM DB_PrI("+ Tm=%f\n",*tPR);
            }
            if(gPR) {   
                *gPR = (hD - ((tmPO->tp+273.15) * sD) ) / 1000.0;
                DB_DNA_TM DB_PrI("+ dG=%f\n",*gPR);
            }
            if(hPR) {   
                *hPR = hD/1000.0;   
            }
            if(sPR) {   
                *sPR = sD;  
            }
            ok = TRUE;
            break;
        default:
            printf("Bogus algorithm code=%d\n",tmPO->algo);
            ERR("SeqTmThermI","Bogus alg code");
    }
    DB_DNA_TM DB_PrI("<< SeqTmThermI %d\n",ok);
    return(ok);
}
/*************************************************************************
*   Calculate Tm based on only GC / AT counts in sequence
*/
DOUB TmFromGCContD(TM_PARS *tmPO,char *seqS,int len)
{
    int i;
    DOUB tmD;

    DB_DNA_TM DB_PrI(">> TmFromGCContD\n");
    tmD = 0.0;
    for(i=0;i<len;i++)
    {
        switch(seqS[i])
        {
            case 'A':   case 'a':
            case 'T':   case 't':   
                tmD += 2.0; 
                break;
            case 'C':   case 'c':
            case 'G':   case 'g':   
                tmD += 4.0; 
                break;
            default:
                if(tmPO->ambig_ok)
                {
                    break;
                }
                printf("# Bogus DNA code letter = %c\n",seqS[i]);
                ERR("TmFromGCContD","Bad DNA letter");
                return(-TOO_BIG_R);
        }
    }
    DB_DNA_TM DB_PrI("<< TmFromGCContD %f\n",tmD);
    return(tmD);
}
/**************************************************************************
*   Set Tm "initialization" entropy / enthalpy terms for passed sequence
*   If ends of seq differ, set average values for both of end terms
*   Ambiguous bases (non ACGT) are ignored
*/
int SetSeqTmInitTermsI(TM_PARS *tmPO,char *seqS,int len,DOUB *hPD,DOUB *sPD)
{
    int term;
    DOUB hD, sD, phD, psD;

    DB_DNA_TM DB_PrI(">> SetSeqTmInitTermsI len=%d\n",len);
    hD = sD = 0.0;
    /***
    *   5' end
    */
    DB_DNA_TM DB_PrI("+ 5'|%c|\n",seqS[0]);
    term = SeqEndInitIndex(seqS[0]);
    if(term < 0) {
        return(FALSE);
    }
    if(!GetThermoTabTermsI(tmPO, term, 0, len, &hD, &sD)) {
        return(FALSE);
    }
    /***
    *   3' end
    */
    DB_DNA_TM DB_PrI("+ 3'|%c|\n",seqS[len-1]);
    term = SeqEndInitIndex(seqS[len-1]);
    if(term < 0) {
        return(FALSE);
    }
    /***
    *   We're looking into the NN table; last starts len-2
    */
    if(!GetThermoTabTermsI(tmPO, term, len-2, len, &phD, &psD)) {
        return(FALSE);
    }
    hD += phD;
    sD += psD;
    /***
    *   SantaLucia doesn't use the average
    */
    switch(tmPO->algo)
    {
        case ALGO_SANTA:        
        case ALGO_PEYRET: 
            break;
        default:
            hD /= 2.0;
            sD /= 2.0;
    }
    *hPD = hD;
    *sPD = sD;
    DB_DNA_TM DB_PrI("<< SetSeqTmInitTermsI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Get thermo dH and dS values for specific NN term and position from tables 
*/
int GetThermoTabTermsI(TM_PARS *tmPO, int term, int pos, int len, DOUB *hPD, DOUB *sPD)
{
    int h,s,bk;

    bk = len - pos - 1;
    /***
    *   Figure column index; 3' end or 5' end or default:  3' 3'' D 5''' 5'' 5'
    *   dH
    */
    if( bk <= (tmPO->ht_col - tmPO->ht_d) ) {
        h = tmPO->ht_col - bk;
    }
    else if( pos < tmPO->ht_d) {
        h = pos;
    }
    else {
        h = tmPO->ht_d;
    }
/*
printf("xxx term=%2d pos[%2d] h=%d (col=%d d=%d bk=%d)\t",term,pos,h,tmPO->ht_col,tmPO->ht_d,bk);
*/
    if(!GetTableValI(tmPO->htab,term,h,hPD) ) {
        return(FALSE);
    }
    /***
    *   dS
    */
    if( bk <= (tmPO->st_col - tmPO->st_d) ) {
        s = tmPO->st_col - bk;
    }
    else if( pos < tmPO->st_d) {
        s = pos;
    }
    else {
        s = tmPO->st_d;
    }
    if(!GetTableValI(tmPO->stab,term,s,sPD) ) {
        return(FALSE);
    }
/*
printf(" H=%5.2f  S=%5.2f\n",*hPD,*sPD);
*/
    return(TRUE);
}
/**************************************************************************
*   Get thermo dH and dS values for specific NN term and position from tables 
*/
int SeqEndInitIndex(char c)
{
    int i;

    i = BOGUS;
    switch(c)
    {
        case 'A': case 'a':
        case 'T': case 't':
            i = INITAT_IND;
            break;
        case 'C': case 'c':
        case 'G': case 'g':
            i = INITCG_IND;
            break;
    }
    return(i);
}
/****************************************************************************
*   Tally up enthalpy and entropy for each pair in sequence
*/
int TallySeqTmNNTermsI(TM_PARS *tmPO,char *seqS,int len,DOUB *hPD, DOUB *sPD)
{
    int i,ind;
    DOUB phD, psD;

    DB_DNA_TM DB_PrI(">> TallySeqTmNNTermsI len=%d\n",len);
    *hPD = *sPD = 0.0;
    for(i=0; i<(len-1); i++)
    {
        ind = SeqPairIndexI(&seqS[i]);
        if( IS_BOG(ind) && (!tmPO->ambig_ok) )
        {
            PROBLINE;
            printf("# Bogus DNA code letter(s) encountered = [%d]|%c|%c|\n",
                i,seqS[i],seqS[i+1]);
            return(FALSE);
        }
        /***
        *   Ignore bogus index (i.e. as from ambiguous bases) 
        */ 
        if(!IS_BOG(ind))
        {
            GetThermoTabTermsI(tmPO, ind, i, len, &phD, &psD);
            *hPD += phD;
            *sPD += psD;
        }
        DB_DNA_TM 
            DB_PrI("+  [%2d] %c%c=%-2d hR=%7.0f sR=%4.1f (%7.0f %4.1f)\n",
                i,seqS[i],seqS[i+1],ind,phD,psD,*hPD,*sPD);
    }
    DB_DNA_TM DB_PrI("<< TallySeqTmNNTermsI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Calculate's Tm from passed Enthalpy and Entropy
*   Different algorithm's may be used
*/
DOUB TmFromThermoSumsD(TM_PARS *tmPO,DOUB hD,DOUB sD,int len)
{
    DOUB tmD,saltD;

    DB_DNA_TM DB_PrI(">> TmFromThermoSumsD %f %f %d algo=%d\n",hD,sD,tmPO->algo,len);
    switch(tmPO->algo)
    {
        case ALGO_MELTING:
            saltD = tmPO->salt; 
            DB_DNA_TM DB_PrI("+ Melting; salt =%0.9e\n",saltD);
            saltD = log10(saltD/(1.0 + 0.7*saltD));
            DB_DNA_TM DB_PrI("+ Melting; salt' =%0.9e\n",saltD);
            tmD = hD/(sD+1.987*log(tmPO->conc)) + 16.6*saltD - 269.3;   
            break;
        case ALGO_OLIGO:
            saltD = tmPO->salt; 
            DB_DNA_TM DB_PrI("+ Oligo; Mg =%0.9e\n",tmPO->mg);
            saltD += (4.0 * sqrt(tmPO->mg));
            DB_DNA_TM DB_PrI("+ Oligo; salt =%0.9e\n",saltD);
            saltD = log10(saltD/(1.0 + 0.7*saltD));
            DB_DNA_TM DB_PrI("+ Oligo; salt =%0.9e\n",saltD);
            tmD = hD/(sD+1.987*log(tmPO->conc/4.0)) + 16.6*saltD - 273.15;
            break;
        case ALGO_SANTA:
            if(tmPO->conc >= tmPO->conc2) {
                tmD = hD/(sD+1.98722*log(tmPO->conc-0.5*tmPO->conc2));
                tmD -= 273.15;
            }
            else {
                tmD = hD/(sD+1.98722*log(tmPO->conc2-0.5*tmPO->conc));
                tmD -= 273.15;
            }
            break;
        default:
            printf("Bad alg=%d\n",tmPO->algo);
            ERR("TmFromThermoSumsD","Bogus alg");
            tmD = -TOO_BIG_R;
    }
    DB_DNA_TM DB_PrI("<< TmFromThermoSumsD %f\n",tmD);
    return(tmD);
}
/**************************************************************************
*   Find length of passed sequence that acheives target Themo value, targD
*       Assumes thermo value is low at short length (e.g. Tm or frac bound)
*       and increases with length.
*   The process can start at the beginning or end of the passed seq as:
*       if dir >= 0, the "direction" of seq traversal is 0 >--> end 
*       else, the "direction" is backward; i.e. end >--> 0      
*   NOTE: Thermo calculation is independent of sampling direction, so
*       if the thermo calculation is not symmetric (i.e. MGB) the Tm strand
*       direction has to be explictly set *before* this call.
*
*   Returns sequence length for *OVERSHOT* Thermo if OK; 
*   BOGUS if problem (i.e. too short to reach target)
*   If ltmPD / htmPD are non-NULL, actual thermo bounding targD are set to these
*/
int SeqLenForThermI(TM_PARS *tmPO, char *seqS, int len, DOUB targD, int dir,
    DOUB *ltmPD, DOUB *htmPD, int do_fds)
{
    int i,n;
    DOUB tmD,otmD,dgD;
    
    VALIDATE(tmPO,TM_PARS_ID);
    DB_TM_LEN {
        DB_PrI(">> SeqLenForThermI len=%d targ=%5.2f do_fds=%d\n",len,targD,do_fds);
        DB_PrI("+ |"); PrintString(seqS,len,outfileGPF); DB_PrI("|\n");
        DB_PrI("+ growth-dir=%d, Tm-strand=%d\n",dir,tmPO->dir);
    }
    n = BOGUS;
    otmD = - TOO_BIG_R;
    if(ltmPD) {
        *ltmPD = -TOO_BIG_R;
    }
    if(htmPD) {
        *htmPD = -TOO_BIG_R;
    }
    /***
    *   Sample increasing lengths until hit target Tm
    */
    for(i=MIN_TM_LEN;i<len;i++)
    {
        if(i>=MAX_TM_LEN) {
            break;
        }
        /***
        *   Backwards (upstream) from sequence end
        */
        if(dir < 0) {
            if(!SeqTmThermI(tmPO,&seqS[len-i],i,&tmD,&dgD,NULL,NULL)) {
                tmD = BAD_R;
            }
        }
        else {
            if(!SeqTmThermI(tmPO,seqS,i,&tmD,&dgD,NULL,NULL)) {
                tmD = BAD_R;
                break;
            }
        }
        DB_TM_LEN DB_PrI("+ len[%2d] tm=%f dG=%f",i,tmD,dgD);
        /***
        *   If fraction bound, calculate that to replace Tm
        */
        if(do_fds) {
            tmD = Fraction2D(tmPO->conc, tmPO->conc2, tmPO->tp, dgD);
            DB_TM_LEN DB_PrI(" fb%f",tmD);
        }
        /***
        *   We've passed the target; done sampling
        */
        if(tmD > targD) {
            DB_TM_LEN DB_PrI(" >targ DONE %f %f\n",otmD,tmD);
            n = i;
            break;
        }
        DB_TM_LEN DB_PrI("\n");
        otmD = tmD;
    }
    /***
    *   Make sure we've actually passed the target Tm
    */
    if(tmD < targD) {
        DB_TM_LEN DB_PrI("<< SeqLenForThermI TOO SHORT = BOGUS\n");
        return(BOGUS);
    }
    /***
    *   Set Tms and return length
    */
    if(ltmPD) {
        *ltmPD = otmD;
    }
    if(htmPD) {
        *htmPD = tmD;
    }
    DB_TM_LEN DB_PrI("<< SeqLenForThermI %d\n",n);
    return(n);
}
/**************************************************************************
*   Calculate Tm with PNA algorithm;
*/
DOUB CalcPNATmI(TM_PARS *tmPO,char *seqS,int len)
{
    int i,n_pyr;
    DOUB tmD,hD,sD,hsumD,ssumD;

    DB_DNA_TM DB_PrI(">> CalcPNATmI\n");
    /***
    *   Tally sequence pairs and pyridines
    */
    hsumD = ssumD = 0.0;
    n_pyr = 0;
    for(i=0;i<(len-1);i++)
    {
        /***
        *   Entropy / enthalpy terms
        */
        hD = sD = 0.0;
        switch(seqS[i])
        {
            case 'A': case 'a':
                switch(seqS[i+1])
                {
                    case 'A': case 'a': 
                        hD = -8400;    sD = -23.6; break;   
                    case 'C': case 'c': 
                        hD = -8600;    sD = -23.0; break;   
                    case 'G': case 'g': 
                        hD = -6100;    sD = -16.1; break;   
                    case 'T': case 't': 
                        hD = -6500;    sD = -18.8; break;   
                }
                break;
            case 'C':
                n_pyr++;
                switch(seqS[i+1])
                {
                    case 'A': case 'a': 
                        hD = -7400;    sD = -19.3; break; 
                    case 'C': case 'c': 
                        hD = -6700;    sD = -15.6; break;
                    case 'G': case 'g': 
                        hD = -10100;   sD = -25.5; break;
                    case 'T': case 't': 
                        hD = -6100;    sD = -16.1; break;
                }
                break;
            case 'G':
                switch(seqS[i+1])
                {
                    case 'A': case 'a': 
                        hD = -7700;    sD = -20.3; break; 
                    case 'C': case 'c': 
                        hD = -11100;   sD = -28.4; break;
                    case 'G': case 'g': 
                        hD = -6700;    sD = -15.6; break;
                    case 'T': case 't': 
                        hD = -8600;    sD = -23.0; break;
                }
                break;
            case 'T':
                n_pyr++;
                switch(seqS[i+1])
                {
                    case 'A': case 'a': 
                        hD = -6300;    sD = -18.5; break;
                    case 'C': case 'c': 
                        hD = -7700;    sD = -20.3; break;
                    case 'G': case 'g': 
                        hD = -7400;    sD = -19.3; break;
                    case 'T': case 't': 
                        hD = -8400;    sD = -23.6; break;
                }
                break;

        }
        hsumD += hD;
        ssumD += sD;
    }
    if( (seqS[len-1]=='C') || (seqS[len-1]=='c') ||
        (seqS[len-1]=='T') || (seqS[len-1]=='t') )
    {
        n_pyr++;
    }
    DB_DNA_TM DB_PrI("+ H=%f S=%f npyr=%d\n",hsumD,ssumD,n_pyr);
    /***
    *   Calculate Tm
    */
    tmD = hsumD / (ssumD + 1.987 * log(tmPO->conc/4.0)) - 273.15;
    tmD *= PNA_TM_SC_D;
    tmD += PNA_TM_CON_D;
    tmD += (PNA_TM_PYR_D * (DNUM(n_pyr)/DNUM(len)));
    tmD += (PNA_TM_LEN_D * DNUM(len));
    DB_DNA_TM DB_PrI("<< CalcPNATmI %f\n",tmD);
    return(tmD);
}
