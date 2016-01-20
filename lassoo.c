/*
* lassoo.c
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
#include <ctype.h>
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "score.h"
#include "bitpool.h"
#include "numlist.h"
#include "lassoo.h"


#define DB_LASP if(DB[95])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(LassooI(argc,argv),NULL) ); }
/**************************************************************************/
void LassooUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <elib>  [...options]\n");
    printf("   <elib>     External library file (bit strings)\n");
    printf("   -ilib XXX  Internal library file XXX (bit strings)\n");
    printf("   -glib XXX  Good reference library file XXX (bit strings)\n");
    printf("   -blib XXX  Bad reference library file XXX (bit strings)\n");
    printf("   -out XXX   Output file XXX\n");
    printf("   -pick #    Number of external records to select\n");
    printf("   -win #     Window size for score updates\n");
    printf("   -dep #     Distance function endpoint (Meaning depends on form)\n");
    printf("   -dmax #    Distance function max; Beyond this scores = 0.0\n");
    printf("   -dgau      Gaussian distance function: To 0.01%% at endpoint\n");
    printf("   -dinv      Inverse distance function: 1/X \n");
    printf("   -dlin      Linear distance function: 1 sloping to 0 past endpoint\n");
    printf("   -dst       Step distance function: 1 up to max\n");
    printf("   -ico #     Internal library weighting coefficent\n");
    printf("   -gco #     Good library weighting coefficent\n");
    printf("   -bco #     Bad library weighting coefficent\n");
    printf("   -stp       Strict parsing; Default is to skip input errors\n");
    printf("   -quiet     Minimal feedback\n");
}
/**************************************************************************/
int LassooI(int argc,char **argv)
{
    int cyc, do_inv, do_gau, do_lin, do_st;
    LASSOO *lasPO;

    do_inv = do_gau = do_lin = do_st = FALSE;
    lasPO = CreateLassooPO();
    if(!ParseArgsI(argc,argv,
        "S -out S -ilib S -glib S -blib S -pick I -win I -dep I -dmax I\
        -ico D -gco D -bco D -quiet B -stp B\
        -dgau B -dinv B -dlin B -dst B",
        lasPO->elibname, lasPO->outname, 
        lasPO->ilibname, lasPO->glibname, lasPO->blibname, 
        &lasPO->pick, &lasPO->win, &lasPO->dfep, &lasPO->dfmax,
        &lasPO->ilibco, &lasPO->glibco, &lasPO->blibco,
        &lasPO->verbose, &lasPO->istrict,
        &do_gau, &do_inv, &do_lin, &do_st,
        (int *)NULL))
    {
        LassooUse();
        CHECK_LASSOO(lasPO);
        return(FALSE);
    }
    /***
    *   Guess input format and override with command line
    */
    if(do_inv) {
        lasPO->dform = LASDF_INV;
    }
    else if(do_gau){
        lasPO->dform = LASDF_GAUS;
    }
    else if(do_lin){
        lasPO->dform = LASDF_LINE;
    }
    else if(do_st){
        lasPO->dform = LASDF_STEP;
    }
    if(!CheckLassooOptionsI(lasPO)) {
        ABORTLINE;
        CHECK_LASSOO(lasPO);
        return(FALSE);
    }
    WriteLassooBanner(lasPO,NULL);
    /***
    *   Set up for use
    */
    if(!RealizeLassooI(lasPO)) {
        ABORTLINE;
        CHECK_LASSOO(lasPO);
        return(FALSE);
    }
    /***
    *   Initialize static scores
    */
    if(!CalcLassooInitalScoresI(lasPO)) {
        ABORTLINE;
        CHECK_LASSOO(lasPO);
        return(FALSE);
    }
    /***
    *   Output story then party until we're done
    */
    WriteLassooHeader(lasPO,lasPO->out);
    cyc = 0;
    while( (lasPO->taken < lasPO->pick) && (lasPO->remain > 0) )
    {
        cyc++;
        CalcLassooScoresI(lasPO);
        PickBestExternalsI(lasPO,cyc);
        WriteLassooOutput(lasPO,cyc,lasPO->out);
        if( (lasPO->taken < lasPO->pick) && (lasPO->remain > 0) ) {
            UpdateInternalScoresI(lasPO,cyc);
        }
    }
    /***
    *   All done
    */
    CHECK_LASSOO(lasPO);
    return(TRUE);
}
/*****************************************************************************
*   Create main lassoo data struct
*/
LASSOO *CreateLassooPO()
{
    LASSOO *lasPO;

    if(! (lasPO = (LASSOO *)ALLOC(1,sizeof(LASSOO)) ) ) {
        printf("# Failed to allocate working object\n");
        return(NULL);
    }
    lasPO->ID = LASSOO_ID;
    InitLassoo(lasPO);
    return(lasPO);
}
/*****************************************************************************
*   Free datastructure and substructs
*/
int DestroyLassooI(LASSOO *lasPO)
{
    VALIDATE(lasPO,LASSOO_ID);
    CHECK_NFILE(lasPO->out,lasPO->outname);
    CHECK_BITPOOL(lasPO->glib);
    CHECK_BITPOOL(lasPO->blib);
    CHECK_BITPOOL(lasPO->elib);
    CHECK_BITPOOL(lasPO->ilib);
    CHECK_NUMLIST(lasPO->status);
    CHECK_NUMLIST(lasPO->prank);
    CHECK_NUMLIST(lasPO->iscores);
    CHECK_NUMLIST(lasPO->gscores);
    CHECK_NUMLIST(lasPO->bscores);
    CHECK_SCOREC(lasPO->srank);
    FREE(lasPO);
    return(TRUE);
}
/*****************************************************************************
*   Set null / default values
*/
void InitLassoo(LASSOO *lasPO)
{
    VALIDATE(lasPO,LASSOO_ID);
    INIT_S(lasPO->outname);
    INIT_S(lasPO->glibname);
    lasPO->glib = NULL;
    lasPO->nglib = 0;
    INIT_S(lasPO->blibname);
    lasPO->blib = NULL;
    lasPO->nblib = 0;
    INIT_S(lasPO->elibname);
    lasPO->elib = NULL;
    lasPO->nelib = 0;
    lasPO->remain = 0;
    lasPO->taken = 0;
    INIT_S(lasPO->ilibname);
    lasPO->ilib = NULL;
    lasPO->nilib = 0;
    lasPO->status = NULL;
    lasPO->prank = NULL;
    lasPO->iscores = NULL;
    lasPO->gscores = NULL;
    lasPO->bscores = NULL;
    lasPO->ilibco = DEF_ICOEF;
    lasPO->glibco = DEF_GCOEF;
    lasPO->blibco = DEF_BCOEF;
    lasPO->srank = NULL;
    lasPO->iform = DEF_IFORM;
    lasPO->istrict = FALSE;
    lasPO->pick = DEF_PICK;
    lasPO->win = DEF_WIN;
    lasPO->dform = DEF_DFUNC;
    InitArrayI(lasPO->dfunc, IS_DOUB, 0, DFMAX, 0.0);
    lasPO->dfep = DEF_DFEP;
    lasPO->dfmax = BAD_I;
    lasPO->verbose = DEF_VERB;
}
/*************************************************************************
*   Check for option consistency
*/
int CheckLassooOptionsI(LASSOO *lasPO)
{
    if( (lasPO->dfep < 1) || (lasPO->dfep > DFMAX) ) {
        PROBLINE;
        printf("Bad distance function end point: %d\n",lasPO->dfep);
        printf(" Allowed values = %d to %d\n",1,DFMAX);
        return(FALSE);
    }
    if(lasPO->dfmax < 1) {
        lasPO->dfmax = lasPO->dfep;
    } 
    if(lasPO->dfmax > DFMAX) {
        PROBLINE;
        printf("Bad distance function max: %d\n",lasPO->dfmax);
        printf(" Allowed values = %d to %d\n",1,DFMAX);
        return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*   Set up things for use
*/
int RealizeLassooI(LASSOO *lasPO)
{
    /***
    *   Load input files; 
    */
    if(!LoadLassooPoolsI(lasPO)) {
        return(FALSE);
    }
    lasPO->remain = lasPO->nelib;
    /***
    *   Check that all pools have same bit dimension
    */
    GetBitpoolDimsI(lasPO->elib, &lasPO->bdim, NULL);
    if(!CompatBitpoolBitDimI(lasPO->ilib,lasPO->bdim)) {
        return(FALSE);
    }
    if(!CompatBitpoolBitDimI(lasPO->glib,lasPO->bdim)) {
        return(FALSE);
    }
    if(!CompatBitpoolBitDimI(lasPO->blib,lasPO->bdim)) {
        return(FALSE);
    }
    /***
    *   Allocate number list and score space
    */
    lasPO->status = CreateNumlistPO(IS_INT, NULL, lasPO->nelib);
    lasPO->prank =  CreateNumlistPO(IS_INT, NULL, lasPO->nelib);
    lasPO->iscores = CreateNumlistPO(IS_DOUB, NULL, lasPO->nelib);
    lasPO->gscores = CreateNumlistPO(IS_DOUB, NULL, lasPO->nelib);
    lasPO->bscores = CreateNumlistPO(IS_DOUB, NULL, lasPO->nelib);
    lasPO->srank = CreateScorecsPO(lasPO->nelib);
    if( (!lasPO->iscores) || (!lasPO->gscores) || (!lasPO->bscores) ||
        (!lasPO->status) || (!lasPO->prank) || (!lasPO->srank) )
    {
        return(FALSE);
    }
    /***
    *   Set distance function
    */
    SetLassooDistFunctI(lasPO);
    /***
    *   Output file?
    */
    if(!NO_S(lasPO->outname)) {
        if(!(lasPO->out=OpenUFilePF(lasPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Set distance function array
*/
int SetLassooDistFunctI(LASSOO *lasPO)
{
    int ok;

    ok = FALSE;
    switch(lasPO->dform) 
    {
        case LASDF_GAUS:    
            ok = SetLasDistFuncGausI(lasPO); 
            break;
        case LASDF_STEP:    
            ok = SetLasDistFuncStepI(lasPO); 
            break;
        case LASDF_LINE:
            ok = SetLasDistFuncLineI(lasPO); 
            break;
        case LASDF_INV:
            ok = SetLasDistFuncInvI(lasPO); 
            break;
        default:
            printf("Bogus dform=%d\n",lasPO->dform);
            ERR("SetLassooDistFunctI","Bad dist func code");
    }
    return(ok);
}
/****************************************************************************/
int SetLasDistFuncGausI(LASSOO *lasPO)
{
    int i;
    DOUB vD;

    for(i=0;i<=lasPO->dfmax;i++)
    {
        vD = 2.1458 * DNUM(i) / DNUM(lasPO->dfep);
        vD = vD * vD;
        vD = exp(-vD);
        lasPO->dfunc[i] = vD;
    }
    return(i);
}
/****************************************************************************/
int SetLasDistFuncStepI(LASSOO *lasPO) 
{
    int ok;

    ok = InitArrayI(lasPO->dfunc, IS_DOUB, 0, (lasPO->dfep +1), 1.0);
    lasPO->dfmax = lasPO->dfep;
    return(ok);
}
/****************************************************************************/
int SetLasDistFuncLineI(LASSOO *lasPO) 
{
    int i, min;
    DOUB vD,dD;

    min = MIN_NUM(lasPO->dfep, lasPO->dfmax);
    dD = 1.0 / DNUM(lasPO->dfep + 1);
    vD = 1.0;
    for(i=0;i<=min;i++)
    {
        lasPO->dfunc[i] = vD;
        vD -= dD;
    }
    lasPO->dfmax = min;
    return(i);
}
/****************************************************************************/
int SetLasDistFuncInvI(LASSOO *lasPO) 
{
    int i;

    for(i=0;i<=lasPO->dfmax;i++)
    {
        lasPO->dfunc[i] = 1.0 / DNUM(i + 1);
    }
    return(i);
}
/****************************************************************************/
int FillLasDistFuncNameStringI(int df, char *nameS)
{
    switch(df) 
    {
        case LASDF_GAUS:    sprintf(nameS,"Gaussian");  break;
        case LASDF_STEP:    sprintf(nameS,"Step");      break;
        case LASDF_LINE:    sprintf(nameS,"Linear");    break;
        case LASDF_INV:     sprintf(nameS,"Inverse");   break;
        default:
            printf("Bogus dform=%d\n",df);
            ERR("FillLasDistFuncNameStringI","Bad dist func code");
            return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*   Write loading of bitstring pools
*/
int LoadLassooPoolsI(LASSOO *lasPO)
{
    /***
    *   External library (must have)
    */
    if(!LoadThisLassooPoolI(lasPO,LASP_EXT)) {
        return(FALSE);
    }
    /***
    *   Internal library (contitional)
    */
    if(!NO_S(lasPO->ilibname)) {
        if(!LoadThisLassooPoolI(lasPO,LASP_INT)) {
            return(FALSE);
        }
    }
    /***
    *   Good library (contitional)
    */
    if(!NO_S(lasPO->glibname)) { 
        if(!LoadThisLassooPoolI(lasPO,LASP_GOOD)) {
            return(FALSE);
        }
    }
    /***
    *   Bad library (contitional)
    */
    if(!NO_S(lasPO->blibname)) {
        if(!LoadThisLassooPoolI(lasPO,LASP_BAD)) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Fill name of pool type into passed string   
*/
int FillLassooPoolTypeStringI(int what,char *typeS)
{
    INIT_S(typeS);
    switch(what)
    {
        case LASP_EXT: sprintf(typeS,"External");   break;
        case LASP_INT: sprintf(typeS,"Internal");   break;
        case LASP_GOOD: sprintf(typeS,"Good");      break;
        case LASP_BAD: sprintf(typeS,"Bad");        break;
        default:
            printf("Bogus what=%d\n",what);
            ERR("FillLassooPoolTypeStringI","Bad pool indicator code");
            return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*   Attempt to load bitstring pool
*/
int LoadThisLassooPoolI(LASSOO *lasPO,int what)
{
    int n;
    char whatS[NSIZE];

    n = 0;
    switch(what)
    {
        case LASP_EXT:
            n = GetBitpoolI(lasPO->elibname, lasPO->iform, lasPO->istrict, &lasPO->elib);
            break;
        case LASP_INT:
            n = GetBitpoolI(lasPO->ilibname, lasPO->iform, lasPO->istrict, &lasPO->ilib);
            break;
        case LASP_GOOD:
            n = GetBitpoolI(lasPO->glibname, lasPO->iform, lasPO->istrict, &lasPO->glib);
            break;
        case LASP_BAD:
            n = GetBitpoolI(lasPO->blibname, lasPO->iform, lasPO->istrict, &lasPO->blib);
            break;
        default:
            printf("Bogus what=%d\n",what);
            ERR("LoadThisLassooPoolI","Bad pool indicator code");
            return(FALSE);
    }
    /***
    *   Did we get anything?
    */
    if(n<1) {
        return(FALSE);
    }
    /***
    *   Set what we got
    */
    switch(what)
    {
        case LASP_EXT: 
            lasPO->nelib = n; 
            break;
        case LASP_INT: 
            lasPO->nilib = n; 
            SetBitpoolCoef(lasPO->ilib,lasPO->ilibco);
            break;
        case LASP_GOOD: 
            lasPO->nglib = n; 
            SetBitpoolCoef(lasPO->glib,lasPO->glibco);
            break;
        case LASP_BAD: 
            lasPO->nblib = n; 
            SetBitpoolCoef(lasPO->blib,lasPO->blibco);
            break;
    }
    if(lasPO->verbose) {
        FillLassooPoolTypeStringI(what,whatS);
        printf("# Loaded %5d %s records\n",n,whatS);
    }
    return(TRUE);
}
/****************************************************************************
*   Check new bitstring pool against passed bitstring dimension
*/
int CompatBitpoolBitDimI(BITPOOL *bpPO,int dim)
{
    int b;

    if(!bpPO) {
        return(TRUE);
    }
    VALIDATE(bpPO,BITPOOL_ID);
    GetBitpoolDimsI(bpPO, &b, NULL);
    if(b != dim) {
        PROBLINE;
        printf("Incompatibile bitstring dimension for %s\n",bpPO->name);
        printf("  Expecting %d, Found %d\n",dim,b);
        return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*   Report program banner info
*/
void WriteLassooBanner(LASSOO *lasPO, FILE *outPF)
{
    HAND_NFILE(outPF);
    if(!lasPO->verbose)
    {
        return;
    }
    fprintf(outPF,"# %s, %s %s %s\n",VERSION_S,BD_S,__DATE__,__TIME__);
    fprintf(outPF,"# %s\n",RTK_S);
    TimeStamp("# ",outPF);
    fprintf(outPF,"#\n");
}
/****************************************************************************
*   Report program inputs and parameters
*/
void WriteLassooHeader(LASSOO *lasPO, FILE *outPF)
{
    HAND_NFILE(outPF);
    if(!lasPO->verbose)
    {
        return;
    }
    if(outPF != stdout)
    {
        WriteLassooBanner(lasPO, outPF);
    }
    WriteBitpoolInfo(lasPO->elib, "External library", outPF);
    WriteBitpoolInfo(lasPO->ilib, "Internal library", outPF);
    WriteBitpoolInfo(lasPO->glib, "Good reference library", outPF);
    WriteBitpoolInfo(lasPO->blib, "Bad reference library", outPF);
    WriteLassooPickVars(lasPO,outPF);
    WriteLassooDistFunc(lasPO,outPF);
    WriteLassooColLabs(lasPO, outPF);
}
/****************************************************************************
*   Report contents of single bitstring pool
*/
void WriteBitpoolInfo(BITPOOL *bpPO, char *whatS, FILE *outPF)
{
    if(!bpPO) {
        return;
    }
    VALIDATE(bpPO,BITPOOL_ID);
    if(bpPO->coef == 0.0) {
        return;
    }
    HAND_NFILE(outPF);
    fprintf(outPF,"##########################\n");
    fprintf(outPF,"# %s\n",whatS);
    fprintf(outPF,"#   Name: %s\n",bpPO->name);
    fprintf(outPF,"#   Size: %d\n",bpPO->num);
    fprintf(outPF,"#   Coef: %5.4f\n",bpPO->coef);
    return;
}
/****************************************************************************
*   Report picking variable parameters
*/
void WriteLassooPickVars(LASSOO *lasPO, FILE *outPF)
{
    HAND_NFILE(outPF);
    fprintf(outPF,"##########################\n");
    fprintf(outPF,"# Compound selection settings\n");
    fprintf(outPF,"#   Pick   = %d\n",lasPO->pick);
    fprintf(outPF,"#   Window = %d\n",lasPO->win);
}
/****************************************************************************
*   Report distance function variables
*/
void WriteLassooDistFunc(LASSOO *lasPO, FILE *outPF)
{
    int i;
    char nameS[DEF_BS];

    HAND_NFILE(outPF);
    fprintf(outPF,"##########################\n");
    fprintf(outPF,"# Distance function settings\n");
    FillLasDistFuncNameStringI(lasPO->dform,nameS);
    fprintf(outPF,"#   Functional form    = %s\n",nameS);
    fprintf(outPF,"#   End point          = %d bits\n",lasPO->dfep);
    fprintf(outPF,"#   Maximum (truncate) = %d bits\n",lasPO->dfmax);
    fprintf(outPF,"#   Values:\n");
    for(i=0; i<=(lasPO->dfmax+1); i++)
    {
        fprintf(outPF,"#    dBit_%02d = %f\n",i,lasPO->dfunc[i]);
    }
}
/****************************************************************************
*   Write output column lables
*/
void WriteLassooColLabs(LASSOO *lasPO, FILE *outPF)
{
    HAND_NFILE(outPF);
    if(!lasPO->verbose)
    {
        return;
    }
    /***
    *   Header line story
    */
    fprintf(outPF,"%-10s\t%s\t%s\t%s","#Name","Rank","Cyc","FScore");
    if(lasPO->ilib) {
        fprintf(outPF,"\t%s","IScore  ");
    } 
    if(lasPO->glib) {
        fprintf(outPF,"\t%s","GScore  ");
    }
    if(lasPO->blib) {
        fprintf(outPF,"\t%s","BScore  ");
    }
    fprintf(outPF,"\n");
}
/****************************************************************************
*   Dump current cycle of chosen records
*/
void WriteLassooOutput(LASSOO *lasPO, int cyc, FILE *outPF)
{
    int i,e,r,v;
    DOUB vD;
    char nameS[NSIZE];

    HAND_NFILE(outPF);
    if(lasPO->verbose) {
        fprintf(outPF,"# Cycle %d\n",cyc);
    }
    /***
    *   Each record with current cycle number
    */
    for(i=0;i<lasPO->nelib;i++)
    {
        e = lasPO->srank[i].id;
        /***
        *   Status = this cycle
        */
        GetNumlistIntI(lasPO->status,e,&v);
        if(v != cyc) {
            continue;
        }
        GetThisBitpoolNameI(lasPO->elib,e,nameS); 
        /***
        *   Minimal case
        */
        if(!lasPO->verbose) {
            fprintf(outPF,"%-10s\n",nameS);
            continue;
        }
        /***
        *   Normal output
        */
        GetNumlistIntI(lasPO->status,e,&v);
        GetNumlistIntI(lasPO->prank,e,&r);
        fprintf(outPF,"%-10s\t%4d %2d",nameS,r,v);
        fprintf(outPF,FSCO_PRINT,lasPO->srank[i].sc);
        if(lasPO->ilib) {
            GetNumlistDoubI(lasPO->iscores,e,&vD);
            fprintf(outPF,PSCO_PRINT,vD);
        }
        if(lasPO->glib) {
            GetNumlistDoubI(lasPO->gscores,e,&vD);
            fprintf(outPF,PSCO_PRINT,vD);
        }
        if(lasPO->blib) {
            GetNumlistDoubI(lasPO->bscores,e,&vD);
            fprintf(outPF,PSCO_PRINT,vD);
        }
        fprintf(outPF,"\n");
    }
    fflush(outPF);
}
/****************************************************************************
*   Update internal scores based on newly added members for single cycle
*/
int UpdateInternalScoresI(LASSOO *lasPO,int cyc)
{
    int i,j,d,n, *statusPI;
    DOUB aD,scD,*iscoresPD;
    BITPTR *bitsPI,*pbitsPI;

    /***
    *   Get raw array pointers (to speed up; makes obvious difference!)
    */
    if(!GetNumlistPtrDoubsI(lasPO->iscores, &iscoresPD, &n)) {
        printf("n=%d nelib=%d\n", n, lasPO->nelib);
        ERR("UpdateInternalScoresI","iscore numbers don't match");
        return(FALSE);
    }
    if(!GetNumlistPtrIntsI(lasPO->status, &statusPI, &n)) {
        printf("n=%d nelib=%d\n", n, lasPO->nelib);
        ERR("UpdateInternalScoresI","status numbers don't match");
        return(FALSE);
    }
    /***
    *   Count how many in current cycle
    */
    n = 0;
    for(i=0;i<lasPO->nelib;i++)
    {
        if(statusPI[i] == cyc) {
            n++;
        }
    }
    if(n<1) {
        return(0);
    }
    /***
    *   Update previously computed internal score components
    *   Old score gets scaled by ratio of old / new pool-size normalization
    */
    aD = DNUM(lasPO->nilib + lasPO->taken) / 
         DNUM(lasPO->nilib + lasPO->taken + n);
    for(i=0;i<lasPO->nelib;i++)
    {
        /***
        *   If not external ignore
        */
        if(statusPI[i] > 0) {
            continue;
        }
        iscoresPD[i] = aD * iscoresPD[i];
    }
    /***
    *   Get new score components for all new internals picked this cycle
    *   All remaining externals get updated 
    */
    for(i=0;i<lasPO->nelib;i++)
    {
        /***
        *   Not external
        */
        if(statusPI[i] > 0) {
            continue;
        }
        GetThisBitpoolPtrBitsI(lasPO->elib, i, &bitsPI, NULL);
        scD = 0.0;
        /***
        *   Each member in the current cycle
        */
        for(j=0;j<lasPO->nelib;j++)
        {
            if(statusPI[j] != cyc) {
                continue;
            }
            GetThisBitpoolPtrBitsI(lasPO->elib, j, &pbitsPI, NULL);
            /***
            *   Distance in bits then score
            */
            d = DifBitCountI(bitsPI,pbitsPI,lasPO->bdim);
            if(d <= lasPO->dfmax) {
                scD += (lasPO->dfunc[d] * lasPO->ilibco);
            }
        }
        /***
        *   Normalize score compoent by total internal library size
        */
        scD /= DNUM(lasPO->nilib + lasPO->taken + n);
        iscoresPD[i] += scD;
    }
    return(n);
}
/****************************************************************************
*   Calculate intial score components
*/
int CalcLassooInitalScoresI(LASSOO *lasPO)
{
    if(lasPO->verbose) {
        printf("#\n");
        printf("# Calculating static scores\n");
    }
    /***
    *   Internal reference pool
    */
    if(!NO_S(lasPO->ilibname)) {
        CalcScoresForBitpoolI(lasPO,lasPO->ilib,lasPO->iscores);
    }
    /***
    *   Good reference pool
    */
    if(!NO_S(lasPO->glibname)) {
        CalcScoresForBitpoolI(lasPO,lasPO->glib,lasPO->gscores);
    }
    /***
    *   Bad reference pool
    */
    if(!NO_S(lasPO->blibname)) {
        CalcScoresForBitpoolI(lasPO,lasPO->blib,lasPO->bscores);
    }
    return(TRUE);
}
/****************************************************************************
*   Calculate scores for externals against passed library pool
*/
int CalcScoresForBitpoolI(LASSOO *lasPO, BITPOOL *bpPO, NUMLIST *scoresPO)
{
    int i,j,d,n,t,v;
    BITPTR *bitsPI,*pbitsPI;
    DOUB scD;

    VALIDATE(lasPO,LASSOO_ID);
    VALIDATE(bpPO,BITPOOL_ID);
    /***
    *   Feedback reporting stories
    */
    n = 0;
    t = lasPO->remain * bpPO->num;
    if(lasPO->verbose) {
        printf("#  %d Externals X %d %s = %d total pair-wise scores\n",
            lasPO->remain, bpPO->num, bpPO->name, t);
    }
    /***
    *   Each external bitstring
    */
    for(i=0;i<lasPO->nelib;i++)
    {
        /***
        *   If not external anymore, will have a value in status so ignore
        */
        GetNumlistIntI(lasPO->status,i,&v);
        if(v > 0) {
            continue;
        }
        GetThisBitpoolPtrBitsI(lasPO->elib, i, &bitsPI, NULL);
        /***
        *   For each bitstring in passed pool get distance then score
        */
        scD = 0.0;
        for(j=0;j<bpPO->num;j++)
        {
            n++;
            if( (lasPO->verbose) && ((n%SCUPDATE)==0) ) {
                printf("#   %9d scores   %5.2f%% done\n",n,PERCENT_R(n,t));
                fflush(stdout);
            }
            GetThisBitpoolPtrBitsI(bpPO, j, &pbitsPI, NULL);
            d = DifBitCountI(bitsPI,pbitsPI,bpPO->bsize);
            if(d <= lasPO->dfmax) {
                scD += (lasPO->dfunc[d] * bpPO->coef);
            }
        }
        /***
        *   Normalize to pool size
        */
        scD /= DNUM(bpPO->num);
        SetNumlistDoubI(scoresPO,i,scD);
    }
    return(TRUE);
}
/****************************************************************************
*   Calculate current scores for each external record
*/
int CalcLassooScoresI(LASSOO *lasPO)
{
    int i,v;
    DOUB scD,goodD,badD,inD;

    /***
    *   Sort by ID; Each external rec
    */
    SortScorecIds(lasPO->srank,lasPO->nelib,SORT_ASCEND);
    for(i=0;i<lasPO->nelib;i++)
    {
        /***
        *   If not external anymore, will have non-zero in status so ignore
        */
        GetNumlistIntI(lasPO->status,i,&v);
        if(v>0) {
            continue;
        }
        /***
        *   Score = Good - Bad - Internal
        */
        GetNumlistDoubI(lasPO->gscores,i,&goodD);
        GetNumlistDoubI(lasPO->bscores,i,&badD);
        GetNumlistDoubI(lasPO->iscores,i,&inD);
        scD = goodD - badD - inD;
        lasPO->srank[i].sc = scD;
    }
    SortScorecVals(lasPO->srank,lasPO->nelib,SORT_DECEND);
    return(TRUE);
}
/****************************************************************************
*   Pick best externals into internal library 
*/
int PickBestExternalsI(LASSOO *lasPO,int cyc)
{
    int i,e,n,v,p;

    /***
    *   Sort by value; Each external rec
    */
    SortScorecVals(lasPO->srank,lasPO->nelib,SORT_DECEND);
    n = 0;
    for(i=0;i<lasPO->nelib;i++)
    {
        e = lasPO->srank[i].id;
        /***
        *   Status = 0 for external bitstrings, otherwise already picked
        */
        GetNumlistIntI(lasPO->status,e,&v);
        if(v>0) {
            continue;
        }
        /***
        *   Set this guy as internal for this cycle number
        */
        SetNumlistIntI(lasPO->status,e,cyc);
        p = lasPO->taken + 1;
        SetNumlistIntI(lasPO->prank,e,p);
        /***
        *   Update counts 
        */
        n++;
        lasPO->taken += 1;
        lasPO->remain -= 1;
        if( (n >= lasPO->win) || (lasPO->taken > lasPO->pick) ) {
            break;
        }
    }
    /***
    *   Update counts of chosen and remaining
    */
    return(n);
}
