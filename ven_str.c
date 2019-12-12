/*
* ven_str.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "prim.h"
#include "venpipe.h"

/****************************************************************************
*   Allocate venpipe data struct
*/
VENPIPE *CreateVenpipePO()
{
    VENPIPE *vpPO;

    vpPO = (VENPIPE *)ALLOC(1,sizeof(VENPIPE));
    if(!vpPO) {
        return(NULL);
    }
    vpPO->ID = VENPIPE_ID;
    vpPO->seq = CreateSeqPO(MAX_VSLEN, NULL, NULL);
    vpPO->ven = CreateVenergyPO();
    if( (!vpPO->seq) || (!vpPO->ven) ) {
        CHECK_VENPIPE(vpPO);
        return(NULL);
    }
    /***
    *   Initialize
    */
    if(!InitVenpipeI(vpPO,TRUE)) {
        CHECK_VENPIPE(vpPO);
        return(NULL);
    }
    return(vpPO);   
}
/****************************************************************************
*   Cleans up and frees venpipe data struct
*/
int DestroyVenpipeI(VENPIPE *vpPO)
{
    VALIDATE(vpPO,VENPIPE_ID);
    CHECK_SEQ(vpPO->seq);
    CHECK_VENERGY(vpPO->ven);
    CHECK_FILE(vpPO->in);
    CHECK_NFILE(vpPO->out,vpPO->outname);
    FREE(vpPO);
    return(TRUE);
}
/****************************************************************************
*   Init venpipe data struct settings to defaults
*   If full, fully clean things, else only per-seq info
*/
int InitVenpipeI(VENPIPE *vpPO, int full)
{
    if(full)
    {
        /***
        *   Files / path names
        */
        INIT_S(vpPO->inname);
        INIT_S(vpPO->outname);
        INIT_S(vpPO->vparfile);
        vpPO->do_saltcorrect = DEF_SALTCORRECT;
        vpPO->do_ksapar = FALSE;
        vpPO->iform = BOGUS;
        vpPO->cleanseq = SCLEAN_LOW;
        /***
        *   Global settings / flags
        */
        vpPO->do_pfe = FALSE;
        vpPO->firstb = vpPO->lastb = BAD_I;
        vpPO->rre = vpPO->mrre = FALSE;
        vpPO->do_dmb = FALSE;
        vpPO->dmb_f = vpPO->dmb_l = BAD_I;
        vpPO->mst = vpPO->men = BAD_R;
        vpPO->mj = DEF_MELT_JUMP_D;
        vpPO->do_mask = FALSE;
        vpPO->do_not = FALSE;
        vpPO->do_mbtab = FALSE;
        /***
        *   Thermo parameters
        */
        vpPO->temp = DEF_TEMP;
        vpPO->salt = DEF_SALT;
    }
    /***
    *   Sequence / structure / probes
    */
    INIT_S(vpPO->tname);
    INIT_S(vpPO->tseq);
    INIT_S(vpPO->tseq2);
    INIT_S(vpPO->tss);
    INIT_S(vpPO->tss2);
    vpPO->tlen = 0;
    vpPO->ten = vpPO->ten2 = TOO_BIG_R;
    return(TRUE);
}
/****************************************************************************
*   Allocate venergy data struct
*/
VENERGY *CreateVenergyPO()
{
    VENERGY *vePO;

    /***
    *   Only one of these at once; Check if global already set
    */
    if(vienna_okGI)
    {
        printf("Vienna already set up\n");
        printf("CreateVenergyPO returning NULL\n");
        return(NULL);
    }
    /***
    *   Normal create function
    */
    vePO = (VENERGY *)ALLOC(1,sizeof(VENERGY));
    if(!vePO) {
        return(NULL);
    }
    vePO->ID = VENERGY_ID;
    if(!InitVenergyI(vePO)) {
        CHECK_VENERGY(vePO);
        return(NULL);
    }
    return(vePO);   
}
/****************************************************************************
*   Clean up Vienna environment
*/
int DestroyVenergyI(VENERGY *vePO)
{
    VALIDATE(vePO,VENERGY_ID);
    CleanViennaEnvI();
    FREE(vePO);
    return(TRUE);
}
/****************************************************************************
*   Init structure values and set up Vienna environment
*/
int InitVenergyI(VENERGY *vePO)
{
    VALIDATE(vePO,VENERGY_ID);
    INIT_S(vePO->vparfile);
    vePO->do_saltcorrect = DEF_SALTCORRECT;
    vePO->do_saltkeep = FALSE;
    vePO->salt = DEF_SALT;
    vePO->temp = DEF_TEMP;

/*
    vePO->nogu = DEF_NOGU;
    vePO->noclosinggu = DEF_NOCLOSINGGU;
    vePO->nolonelypairs = DEF_NOLONLEYPAIRS; 
    vePO->tetra_loop = DEF_TETRA_LOOPS;
    vePO->energy_set = DEF_ENERGY_SET;
    vePO->fold_constrained = DEF_FOLD_CONSTRAINED;
    vePO->dangles = DEF_DANGLES;
    vePO->logml = DEF_LOGLM;
*/

    if(!SetUpViennaEnvI(vePO)) {
        return(FALSE);
    }
    return(TRUE);
}
/*************************************************************************
*   Attempt to instantiate venpipe structure for use
*/
int RealizeVenpipeI(VENPIPE *vpPO)
{
    int ok;

    ok = SetUpViennaEnvI(vpPO->ven);
    if(ok) {
        ok = SetViennaSaltI(vpPO->ven,vpPO->salt,vpPO->do_saltcorrect,
            vpPO->do_ksapar);
    }
    if(ok) {
        ok = SetViennaParametersI(vpPO->ven,vpPO->vparfile);
    }
    if(ok) {
        ok = SetViennaTemperatureI(vpPO->ven,vpPO->temp);
    }
    if( !ok ) {
        PROBLINE;
        printf("Problem setting up Vienna environment\n");
        return(FALSE);
    }
    /***
    *   Open files
    */
    if( ! (vpPO->in = OpenUFilePF(vpPO->inname,"r",NULL))) {
        PROBLINE;
        printf("Problem opening input file\n");
        return(FALSE);
    }
    if(!NO_S(vpPO->outname)) {
        if( ! (vpPO->out = OpenUFilePF(vpPO->outname,"w",NULL))) {
            PROBLINE;
            printf("Problem opening output log\n");
            return(FALSE);
        }
    }
    HAND_NFILE(vpPO->out);
    return(TRUE);
}
/****************************************************************************
*   Copy seq from SEQ structure into working buffer
*/
int CopyWorkingSeqI(VENPIPE *venPO, SEQ *seqPO)
{
    int i,n,ok,len;
    char nameS[NSIZE],*seqPC;

    VALIDATE(venPO,VENPIPE_ID);
    len = GetSeqLenI(seqPO);
    if(len > MAX_VSLEN) {
        printf("# %s is too long %d (%d max)\n",nameS,len,MAX_VSLEN);
        return(FALSE);
    }
    if(!GetSeqSeqI(seqPO,&seqPC)) {
        return(FALSE);
    }
    /***
    *   Copy into buffer, extracting only specified range if specified
    */
    n = 0;
    for(i=0;i<len;i++) {
        ok = TRUE;
        /***
        *   Restricting bases?
        */
        if(venPO->firstb > 0) {
            if(venPO->rre) {
                if( (len-i)<venPO->firstb ) {
                    ok = FALSE;
                }
                if( (len-i) > venPO->lastb ) {
                    ok = FALSE;
                }
            }
            else {
                if( (i+1) < venPO->firstb ) {
                    ok = FALSE;
                }
                if( (i+1) > venPO->lastb ) {
                    ok = FALSE;
                }
            }
        }
        /***
        *   Mask or ignore; Seq2 always get un-masked version
        */
        ok = (venPO->do_not) ? !ok : ok;
        if(!ok) {
            if(venPO->do_mask) {
                venPO->tseq[n] = 'n';
            }
            else {
                continue;
            }
        }
        else {
            venPO->tseq[n] = seqPC[i];
        }
        venPO->tseq2[n] = seqPC[i];
        n++;
    }
    venPO->tseq[n] = venPO->tseq2[n] = '\0';
    venPO->tlen = n;
    FillSeqNameStringI(seqPO, nameS, NSIZE);
    strcpy(venPO->tname, nameS);
    /***
    *   If subseq only and dumpign seq, set for display
    *   All lowercase, except explicit range
    */
    if( (venPO->do_ds) && (venPO->firstb > 0) ) {
        SetCaseSeqSubseqI(seqPO, FALSE, -1, -1);
        if(venPO->rre) {
            SetCaseSeqSubseqI(seqPO, TRUE, len - venPO->lastb, len - venPO->firstb + 1);
        }
        else {
            SetCaseSeqSubseqI(seqPO, TRUE, venPO->firstb-1, venPO->lastb);
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Revert U's to T's for display (as they T's need to be U for folding)
*/
void ChangeUtoT(char *seqS,int len)
{
    int i;

    for(i=0;i<len;i++) {
        switch(seqS[i]) {
            case 'U': seqS[i] = 'T'; break;
            case 'u': seqS[i] = 't'; break;
        }
    }
}
/****************************************************************************
*   Switch any T's to U's and set uppercase (needed for folding)
*   Return TRUE if all is well
*   Return FALSE if too long or bogus chars in sequence; IUB codes ignored 
*/
int PrepViennaInputSeqI(char *seqS, int len,char *newS)
{
    int i;

    if(len>=MAX_VSLEN) {
        PROBLINE;
        printf("Maximum length for venpipe folding exceeded\n");
        printf("  Max=%d, current=%d\n",MAX_VSLEN,len);
    }
    for(i=0;i<len;i++)
    {
        switch(seqS[i])
        {
            case 'T':
            case 't':
                newS[i] = 'U';
                break;
            default:
                if(BaseDegeneracyI(seqS[i]) == 0) {
                    PROBLINE;
                    printf("Bogus character in vienna seq = %c @ %d\n",
                        seqS[i],i);
                    return(FALSE);
                }
                newS[i] = TOUPPER(seqS[i]);
        }
    }
    newS[i] = '\0';
    return(TRUE);
}
/****************************************************************************
*   Call vienna interface code with settings from venpipe object
*/
int GetVenpipeStructureI(VENPIPE *vpPO)
{
    int ok;

    VALIDATE(vpPO,VENPIPE_ID);
    ok = GetViennaStructI(vpPO->tseq, vpPO->tlen, vpPO->tss, &vpPO->ten,
                        vpPO->do_pfe, vpPO->tss2, &vpPO->ten2);
    return(ok);
}

