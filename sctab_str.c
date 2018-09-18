/*
* sctab_str.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prim.h"
#include "score.h"
#include "table.h"
#include "stat.h"
#include "scoretab.h"

/*****************************************************************************
*   Create scoretab data type
*/
SCORETAB *CreateScoretabPO()
{
    SCORETAB *stPO;

    if(! (stPO = (SCORETAB *)ALLOC(1,sizeof(SCORETAB)) ) ) {
        printf("Failed to allocate working object\n");
        return(NULL);
    }
    stPO->ID = SCORETAB_ID;
    InitScoretab(stPO);
    return(stPO);
}
/*****************************************************************************
*   Clean up scoretab object
*/
int DestroyScoretabI(SCORETAB *stPO)
{
    int i;

    VALIDATE(stPO,SCORETAB_ID);
    CHECK_NFILE(stPO->out,stPO->outname);
    CHECK_TABLE(stPO->tab);
    CHECK_TABLE(stPO->stab);
    CHECK_NUMLIST(stPO->tvals1);
    CHECK_NUMLIST(stPO->tvals2);
    for(i=0;i<stPO->nscores;i++) {
        CHECK_SCFIELD(stPO->scores[i]);
    }
    CHECK_FREE(stPO->scores);
    FREE(stPO);
    return(TRUE);
}
/*****************************************************************************
*   Initialize scoretab object
*/
void InitScoretab(SCORETAB *stPO)
{
    VALIDATE(stPO,SCORETAB_ID);

    INIT_S(stPO->inname);
    INIT_S(stPO->outname);
    INIT_S(stPO->sdfname);
    INIT_S(stPO->rlisname);
    INIT_S(stPO->clisname);
    INIT_S(stPO->mergname);
    INIT_S(stPO->macpre);
    INIT_S(stPO->macsuf);
    stPO->do_scg = FALSE;
    stPO->out = NULL;
    stPO->owhat = BOGUS;
    stPO->tab = NULL;
    stPO->tvals1 = stPO->tvals2 = NULL;
    stPO->nscores = 0;
    stPO->scores = NULL;
    stPO->rlab = TRUE;
    stPO->clab = TRUE;
    stPO->corn = TRUE;
    stPO->do_skp = FALSE;
    stPO->mval = 1.0;
    stPO->sval = 0.0;
    stPO->blval = -TOO_BIG_R;
    stPO->bhval = TOO_BIG_R;
    stPO->exp = 1.0;
    stPO->do_srow = 0;
    stPO->do_scol = 0;
    strcpy(stPO->pvsep,DEF_PVSEP_S);
    stPO->usp_form = FALSE;
    stPO->quiet = FALSE;
    stPO->do_mrow = FALSE;
    stPO->do_ncv = FALSE;
    stPO->do_qcv = 0;
    stPO->minrow = stPO->mincol = 0;
    stPO->maxrow = stPO->maxcol = TOO_BIG;
    stPO->ccor = BOGUS;
    stPO->cinfo = BOGUS;
    stPO->inf_maxhb = INFO_DEF_MBIN;
    stPO->flglo = stPO->flghi = BAD_D;
    stPO->do_not = FALSE;
    stPO->do_wfm = FALSE;
    stPO->do_dwfm = FALSE;
    stPO->do_symu = FALSE;
    stPO->do_symd = FALSE;
    stPO->do_igd = FALSE;
    stPO->do_tran = FALSE;
    stPO->do_stru = stPO->do_strd = FALSE;
    stPO->do_abs = FALSE;
    stPO->psize = 0;
    stPO->pmin = 0.0;
    stPO->partalg = DEF_PALG;
    INIT_S(stPO->dumpbase);
    stPO->do_wsub = stPO->do_wst = stPO->do_kc = FALSE;
    stPO->seed = BOGUS;
    stPO->do_gax = FALSE;
    stPO->gaxc = stPO->gaxr = 0.0; 
    stPO->do_gam = FALSE;
    stPO->gamf = 0.0;
    stPO->gamg = 1.0; 
    return;
}
/************************************************************************/
void SetOutFormatting(SCORETAB *stPO,int oftw,int oftp)
{
    char formS[DEF_BS];
    int auto_rlab;

    if(stPO->do_ossv) {
        strcpy(stPO->pvsep," ");
    }
    else if(stPO->do_ocsv) {
        strcpy(stPO->pvsep,",");
    }
    else {
        strcpy(stPO->pvsep,"\t");
    }
    /***
    *   Value formatting; Either user-supplied or derive auto 
    */
    INIT_S(formS);
    if( oftp >= 0 ) {
        FloatFormatString(oftw,oftp,formS);
        stPO->usp_form = TRUE;
    }
    else {
        AutoTableOutFormattingI(stPO->tab, TRUE, FALSE);
    }
    SetTablePrintformI(stPO->tab,formS,NULL,stPO->pvsep,NULL,NULL);
    /***
    *   Set row lable width automatically (adds spaces), except if csv output 
    */
    auto_rlab = (stPO->do_ocsv) ? FALSE : TRUE;
    AutoTableOutFormattingI(stPO->tab, FALSE, auto_rlab);
    return;
}
/*************************************************************************
*   Check for output option consistency
*/
int CheckSctOutOptionsI(SCORETAB *stPO)
{
    if( stPO->minrow > stPO->maxrow ) {
        PROBLINE;
        printf("Bad row range: %d to %d\n",stPO->minrow,stPO->maxrow);
        return(FALSE);
    }
    if( stPO->mincol > stPO->maxcol ) {
        PROBLINE;
        printf("Bad column range: %d to %d\n",stPO->mincol,stPO->maxcol);
        return(FALSE);
    }
    if( stPO->bhval < stPO->blval ) {
        PROBLINE;
        printf("Bad value bounds: %f to %f\n",stPO->blval,stPO->bhval);
        return(FALSE);
    }
    /***
    *   Math
    */
    if(stPO->exp <= 0.0) {
        PROBLINE;
        printf("Bad exponent given: %f\n",stPO->exp);
        return(FALSE);
    }
    if(stPO->exp != 1.0) {
        stPO->do_abs = TRUE;
    }
    /***
    *   GA Stuff?
    */
    if( (stPO->gaxc > 0.0) || (stPO->gaxr > 0.0) ) {
        if(NO_S(stPO->mergname) ) {
            PROBLINE;
            printf("GA crossover needs a table to merge\n");
            return(FALSE);
        }
        stPO->do_gax = TRUE;
    }
    if(stPO->gamf > 0.0) {
        if (stPO->gamg < TINY_R) {
            PROBLINE;
            printf("GA guassian SD too small: %f\n",stPO->gamg);
            return(FALSE);
        }
        stPO->do_gam = TRUE;
    }
    Srand(stPO->seed);
    /***
    *   Output set?
    */
    switch(stPO->owhat) {
        case SCTO_FULL:
        case SCTO_RCL:
        case SCTO_CRL:
        case SCTO_RPRO:
        case SCTO_RSTAT:
        case SCTO_CSTAT:
        case SCTO_CINFO:
        case SCTO_CCOR:
        case SCTO_FCCM:
        case SCTO_FSTAT:
        case SCTO_FLAG:
        case SCTO_NONE:
            break;
        default:
            printf("Bogus output code = %d\n",stPO->owhat);
            ERR("CheckSctOutOptionsI","Bogus output code");
            return(FALSE);
    }
    return(TRUE);
}
/***************************************************************************
*   Extra space for temp arrays and load any scores if transforming
*/
int SetUpSctAuxDataI(SCORETAB *stPO)
{
    int n;
    FILE *fPF;
    SCFIELD **scPPO;

    /***
    *   Allocate temp arrays for rows/cols
    */
    if(!SetUpSctTableSpaceI(stPO)) {
        return(FALSE);
    }
    /***
    *   Loading any score transforms?
    */
    if(!NO_S(stPO->sdfname)) {
        if(!(fPF=OpenUFilePF(stPO->sdfname,"r",NULL))) {
            printf("Failed to open score definitions file\n");
            return(FALSE);
        }
        n = LoadScfieldArrayI(fPF,&scPPO);
        CHECK_FILE(fPF);
        if(n<1) {
            printf("Failed to load score definitions from %s\n",stPO->sdfname);
            return(FALSE);
        }
        stPO->scores = scPPO;
        stPO->nscores = n;
        if(!stPO->quiet) {
            printf("# Loaded %d scores from %s\n",stPO->nscores,stPO->sdfname);
            ReportScfieldArrayI(stPO->scores, stPO->nscores, "#\t", NULL);
        }
    }
    /***
    *   Weak first mapping?
    */
    if(stPO->do_wfm) {
        if(!CopyTableI(stPO->tab,TRUE,FALSE,&stPO->stab)) {
            PROBLINE;
            printf("Failed to create duplicate table for weak-first results\n");
            return(FALSE);
        }
    }
    return(TRUE);
}
/*************************************************************************
*   Allocate "extra" space associated with table
*   cols[nrow] and rows[ncol]
*/
int SetUpSctTableSpaceI(SCORETAB *stPO)
{
    /***
    *   Clean anything already here
    */
    CHECK_NUMLIST(stPO->tvals1);
    CHECK_NUMLIST(stPO->tvals2);
    /***
    *   New space; 
    *   Temp arrays allocated with zero space; Added if needed when needed
    */
    stPO->tvals1 = CreateNumlistPO(IS_DOUB, NULL, 0);
    stPO->tvals2 = CreateNumlistPO(IS_DOUB, NULL, 0);
    if( (!stPO->tvals1) || (!stPO->tvals2) ) {
        PROBLINE;
        printf("Failed to allocate temp arrays\n");
        return(FALSE);
    }
    return(TRUE);
}
/*************************************************************************
*   Mask out specific rows / cols
*/
int RestrictTableRangesI(SCORETAB *stPO, TABLE *tabPO)
{
    int r,c,ok,n,report;

    n = 0;
    /***
    *   Row or col range restrictions?
    */
    if(stPO->minrow != 0) {
        HandleMaskTableRowsI(stPO,tabPO);
        n++;
    }
    if(stPO->mincol != 0) {
        HandleMaskTableColsI(stPO,tabPO);
        n++;
    }
    /***
    *   Subset lists of row column names?
    */
    report = 0;
    if(!NO_S(stPO->rlisname)) {
        if ( !report ) {
            ReportWordMatchMods("# ", stPO->do_kc, stPO->do_wst, stPO->do_wsub, NULL);
            report++;
        }
        if(HandleTabRowSubsetI(stPO,tabPO) < 0) {
            return(BOGUS);
        }
        n++;
    }
    if(!NO_S(stPO->clisname)) {
        if ( !report ) {
            ReportWordMatchMods("# ", stPO->do_kc, stPO->do_wst, stPO->do_wsub, NULL);
            report++;
        }
        if(HandleTabColSubsetI(stPO,tabPO) < 0) {
            return(BOGUS);
        }
        n++;
    }
    /***
    *   Invert qualifications?
    */
    if(stPO->do_not) {
        for(r=0;r<tabPO->nrow;r++) {
            ok = GetTableRowMaskI(tabPO,r);
            SetTableRowMaskI(tabPO, r, !ok);
        }
        for(c=0;c<tabPO->ncol;c++) {
            ok = GetTableColMaskI(tabPO,c);
            SetTableColMaskI(tabPO, c, !ok);
        }
        n++;
    }
    /***
    *   Ignore diagonal by setting to zero?
    */
    if(stPO->do_igd) {
        SetTableDiagValI(tabPO,0.0);
    }
    return(n);
}
/******************************************************************************
*   Apply masking to rows outside allowed range
*/
int HandleMaskTableRowsI(SCORETAB *stPO, TABLE *tabPO)
{
    int r,n,ok;

    n = 0;
    for(r=0;r<tabPO->nrow;r++) {
        ok = TRUE;
        if( ((r+1)<stPO->minrow) || ((r+1)>stPO->maxrow) ) {
            ok = FALSE;
            n++;
        }
        SetTableRowMaskI(tabPO,r,ok);
    }
    return(n);
}
/******************************************************************************
*   Apply masking to colums outside allowed range
*/
int HandleMaskTableColsI(SCORETAB *stPO, TABLE *tabPO)
{
    int c,n,ok;

    n = 0;
    for(c=0;c<tabPO->ncol;c++) {
        ok = TRUE;
        if( ((c+1)<stPO->mincol) || ((c+1)>stPO->maxcol) ) {
            ok = FALSE;
            n++;
        }
        SetTableColMaskI(tabPO,c,ok);
    }
    return(n);
}
/***************************************************************************
*   Apply mask to subset of listed row names
*/
int HandleTabRowSubsetI(SCORETAB *stPO, TABLE *tabPO)
{
    int r,ok,m;
    char rowS[NSIZE];
    WORDLIST *wlisPO;

    if( ! (wlisPO = CreateWordlistPO(stPO->rlisname,NSIZE))) {
        PROBLINE;
        printf("Failed to get tokens from %s\n",stPO->rlisname);
        return(FALSE);
    }
    /***
    *   scan each row against name collection
    */
    m = 0;
    for(r=0; r<tabPO->nrow; r++) {
        GetTableRowLabI(tabPO,r,rowS,-1);
        ok = FALSE;
        if(WordInWordlistI(wlisPO, rowS, stPO->do_kc, stPO->do_wst, stPO->do_wsub, NULL)) {
            ok++;
            m++;
        }
        SetTableRowMaskI(tabPO,r,ok);
    }
    CHECK_WORDLIST(wlisPO);
    return(m);
}
/***************************************************************************
*   Apply mask to subset of listed row names
*/
int HandleTabColSubsetI(SCORETAB *stPO, TABLE *tabPO)
{
    int c,ok,m;
    char colS[NSIZE];
    WORDLIST *wlisPO;

    if( ! (wlisPO = CreateWordlistPO(stPO->clisname,NSIZE))) {
        PROBLINE;
        printf("Failed to get tokens from %s\n",stPO->clisname);
        return(FALSE);
    }
    /***
    *   scan each col against name collection
    */
    m = 0;
    for(c=0; c<tabPO->ncol; c++) {
        GetTableColLabI(tabPO,c,colS,-1);
        ok = FALSE;
        if(WordInWordlistI(wlisPO, colS, stPO->do_kc, stPO->do_wst, stPO->do_wsub, NULL)) {
            ok++;
            m++;
        }
        SetTableColMaskI(tabPO,c,ok);
    }
    CHECK_WORDLIST(wlisPO);
    return(m);
}
/******************************************************************************
*   Apply scaling / shift factor / bounds to all values 
*/
int AdjustTableValuesI(SCORETAB *stPO,TABLE *tabPO)
{
    int nr,nc,r,c;
    DOUB vD;
 
    /***
    *   GA mutation?
    */
    if(stPO->do_gam) {
        return(HandleSctGAMutationI(stPO));
    }
    /***
    *   Nothing to do?
    */
    if( (stPO->mval==1.0) && (stPO->sval==0.0) && (stPO->bhval==TOO_BIG_R) && 
        (!stPO->do_abs) && (stPO->exp==1.0) ) { 
        return(TRUE);
    }
    /***
    *   Update values; Handle masking in loops
    */
    nc = GetTableColsI(tabPO,FALSE);
    nr = GetTableRowsI(tabPO,FALSE);
    for(r=0;r<nr;r++) 
    {
        if(!GetTableRowMaskI(tabPO,r)) {
            continue;
        }
        for(c=0;c<nc;c++) 
        {
            if(!GetTableColMaskI(tabPO,c)) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            vD *= stPO->mval;
            vD += stPO->sval;
            LIMIT_NUM(vD,stPO->blval,stPO->bhval);
            if(stPO->do_abs) {
                vD = ABS_VAL(vD);
            }
            if(stPO->exp != 1.0) {
                vD = PowD(vD, stPO->exp);
            }
            SetTableValI(tabPO,r,c,vD);
        }
    }
    /***
    *   Tell the story
    */
    if(!stPO->quiet) {
        if(stPO->mval!=1.0) {
            printf("# Values multiplied by %f\n",stPO->mval);
        }
        if(stPO->sval!=0.0) {
            printf("# Values shifted by %f\n",stPO->sval);
        }
        if(stPO->bhval!=TOO_BIG_R) {
            printf("# Values bounded %f to %f\n",stPO->blval,stPO->bhval);
        }
        if(stPO->do_abs) {
            printf("# Values replaced with absolute value\n");
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Handle make-symmetric options
*/
void HandleSymmetrization(SCORETAB *stPO, TABLE *tabPO)
{
    int r,c;
    DOUB v1D,v2D;

    if( (!stPO->do_symu) && (!stPO->do_symd) ) {
        return;
    }
    for(r=0;r<tabPO->nrow;r++) {
        for(c=r+1;c<tabPO->ncol;c++) {
            GetTableValI(tabPO,r,c,&v1D);
            GetTableValI(tabPO,c,r,&v2D);
            if(stPO->do_symu) {
                v1D = MAX_NUM(v1D,v2D);
            }
            else {
                v1D = MIN_NUM(v1D,v2D);
            }
            SetTableValI(tabPO,r,c,v1D);
            SetTableValI(tabPO,c,r,v1D);
        }
    }
}
