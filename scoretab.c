/*
* scoretab.c
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
#define __MAIN__
#include "prim.h"
#include "score.h"
#include "stat.h"
#include "table.h"
#include "mut_info.h"
#include "scoretab.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit(AllDoneI(ScoreTabI(argc,argv),NULL)); }
/**************************************************************************/
void ScoreTabUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> [...options]\n");
    printf("   <infile>     Value table file\n");
    printf("   -nrlab       NO Rows labels in input\n");
    printf("   -nclab       NO Column labels in input\n");
    printf("   -ncorn       NO 'Corner' token between row / col labels\n");
    printf("   -sk          Skip missing / bogus rows; Default is abort\n");
    printf("   -out XXX     Set output to file XXX\n");
    printf("   -tsp -igd    Transpose / Ignore diagonal\n");
    printf("   -cran # #    Restrict Columns to range # to #\n");
    printf("   -rran # #    Restrict Rows to range # to #\n");
    printf("   -clis XXX    Restrict Column range to names listed in file XXX\n");
    printf("   -rlis XXX    Restrict Row range to names listed in file XXX\n");
    printf("   -kc          Restrict; Keep case for subset name comparisons\n");
    printf("   -wst         Restrict; Word start only needs to match (not full token)\n");
    printf("   -wsub        Restrict; Word substring only needs to match (not full token)\n");
    printf("   -flg # #     Flag values in the range # #\n");
    printf("   -not         Invert flagging / bounding / subset criteria\n");
    printf("   -sdf XXX     Score definition file to transform values\n");
    printf("   -scg         Score globally (not one per column)\n");
    printf("   -ncv         Normalize column values (i.e. 0 +/- 1 s.d.)\n");
    printf("   -qcv #       Quantile column values into # classes\n");
    printf("   -symu -symd  Make symmetric with Up (Max) or Down (Min) values\n");
    printf("   -mult #      Multiply values by #\n");
    printf("   -shif #      Shift values by #\n");
    printf("   -bval # #    Bound values range # #\n");
    printf("   -exp #       Raise values to # power (Abs|x| first)\n");
    printf("   -abs         Absolute value\n");
    printf("   -srow #      Smooth rows with +/- window # (i.e. 1 = -1,0,+1)\n");
    printf("   -scol #      Smooth cols with +/- window # (i.e. 1 = -1,0,+1)\n");
    printf("   -stru -strd  Sort rows Up / Down (ascend / decend)\n");
    printf("   -merg XXX    Merge in table XXX (add values)\n");
    printf("   -mlis        Merge in List of tables XXX (add values)\n");
    printf("   -mrow        Merge rows; New rows = f(Odd, Even); def = add\n");
    printf("   -msub        Merge by Subtracting values (First - Second)\n");
    printf("   -mmul        Merge by Multiplying values (First x Second)\n");
    printf("   -mdiv        Merge by Dividing values (First / Second)\n");
    printf("   -mmin -mmax  Merge via Min / Max values\n");
    printf("   -mapc        Merge via Appending Columns\n");
    printf("   -macp XXX    Merge appending columns with label prefix XXX\n");
    printf("   -macs XXX    Merge appending columns with label suffix XXX\n");
    printf("   -ofull       Output full table\n");
    printf("   -orc -ocr    Output list in Row-Col / Col-Row order\n");
    printf("   -ocsv -ossv  Output comma-seperated value / space-separated value\n");
    printf("   -pfm # #     Output print format # wide # precision (%%#.#f)\n");
    printf("   -rstat       Stats; Row statisticics (i.e. per row)\n");
    printf("   -cstat       Stats; Column statisticics (i.e. per col)\n");
    printf("   -fstat       Stats; Full table statisticics (whole thing)\n");
    printf("   -rpro        Row product\n");
    printf("   -fccm        Full column-column correlation matrix\n");
    printf("   -ccor #      Correlations between column # and remaining cols\n");
    printf("   -cinf #      Mutual information between column # and remaining cols\n");
    printf("   -ifmh #      Mutual information maximum histogram bins #\n");
    printf("   -wfm         Weak-first matching of rows to cols\n");
    printf("   -dwfm        Weak-first dump matching matrix\n");
    printf("   -part #      Weak-first partition into pools of size #\n");
    printf("   -pmin #      Weak-first partition minimum compatability score = #\n");
    printf("   -dmat XXX    Weak-first dump matrix files with base XXX (WF or GA)\n");
    printf("   -seed #      Set random number seed to #\n");
    printf("   -gaxr #      GA row cross-over fraction # between parents\n");
    printf("   -gaxc #      GA col cross-over fraction # between parents\n");
    printf("   -gamf #      GA mutation Fraction # of elements\n");
    printf("   -gamg #      GA mutation Gaussian with s.d. #\n");
    printf("   -quiet       Don't report status, processing, etc\n");
}
/**************************************************************************
*   Main program
*/
int ScoreTabI(int argc, char **argv)
{
    int r,c,ofull,rstat,cstat,fstat,rpro,fmtw,fmtp;
    int fccm,orc,ocr,poolsize;
    SCORETAB *stPO;

    stPO = CreateScoretabPO();
    ofull = rpro = FALSE;
    rstat = cstat = fstat = fccm = orc = ocr = FALSE;
    fmtw = fmtp = BAD_I;
    if(!ParseArgsI(argc,argv,
        "S -sdf S -out S -ofull B -nrlab B -nclab B -ncorn B -rpro B\
        -quiet B -cran I2 -mrow B -rran I2 -rlis S -clis S\
        -pfm I2 -cstat B -rstat B\
        -ccor I -fccm B -fstat B -igd B -mult D -shif D -wfm B -ncv B\
        -bval D2 -merg S -mlis B -msub B -mmul B\
        -mdiv B -orc B -ocr B -flg D2 -not B -dwfm B -scg B -symu B -symd B\
        -part I -pmin D -mmin B -mmax B -dmat S -wst B -kc B -wsub B\
        -tsp B -seed I -gaxr D -gaxc D -gamf D -gamg D -srow I -scol I\
        -stru B -strd B -ocsv B -ossv B -mapc B -sk B -abs B -exp D\
        -macp S -macs S -qcv I -cinf I -ifmh I"
        ,
        stPO->inname, stPO->sdfname, stPO->outname, &ofull,
        &stPO->rlab, &stPO->clab, &stPO->corn, &rpro, &stPO->quiet, 
        &stPO->mincol,&stPO->maxcol, &stPO->do_mrow,
        &stPO->minrow,&stPO->maxrow, &stPO->rlisname, &stPO->clisname, 
        &fmtw,&fmtp, &cstat, &rstat, &stPO->ccor, &fccm, &fstat, &stPO->do_igd,
        &stPO->mval, &stPO->sval, &stPO->do_wfm, &stPO->do_ncv,
        &stPO->blval, &stPO->bhval,
        &stPO->mergname, &stPO->do_mlis, &stPO->do_msub, &stPO->do_mmul,
        &stPO->do_mdiv, &orc, &ocr,
        &stPO->flglo, &stPO->flghi, &stPO->do_not, &stPO->do_dwfm, &stPO->do_scg,
        &stPO->do_symu, &stPO->do_symd, &stPO->psize, &stPO->pmin,
        &stPO->do_mmin, &stPO->do_mmax, stPO->dumpbase, &stPO->do_wst, 
        &stPO->do_kc, &stPO->do_wsub, &stPO->do_tran,
        &stPO->seed, &stPO->gaxr, &stPO->gaxc, &stPO->gamf, &stPO->gamg,
        &stPO->do_srow, &stPO->do_scol, &stPO->do_stru, &stPO->do_strd,
        &stPO->do_ocsv, &stPO->do_ossv, &stPO->do_mapc, &stPO->do_skp,
        &stPO->do_abs, &stPO->exp, stPO->macpre, stPO->macsuf, &stPO->do_qcv,
        &stPO->cinfo, &stPO->inf_maxhb,
        (int *)NULL))
    {
        ScoreTabUse();
        CHECK_SCORETAB(stPO);
        return(FALSE);
    }
    /***
    *   What output options?
    */
    if(fccm)
    {   stPO->owhat = SCTO_FCCM; }
    else if(!IS_BOG(stPO->ccor) )
    {   stPO->owhat = SCTO_CCOR; }
    else if(!IS_BOG(stPO->cinfo) )
    {   stPO->owhat = SCTO_CINFO; }
    else if(fstat)
    {   stPO->owhat = SCTO_FSTAT; }
    else if(cstat)
    {   stPO->owhat = SCTO_CSTAT; }
    else if(rstat)
    {   stPO->owhat = SCTO_RSTAT; }
    else if(ofull)
    {   stPO->owhat = SCTO_FULL; }
    else if(orc)
    {   stPO->owhat = SCTO_RCL; }
    else if(ocr)
    {   stPO->owhat = SCTO_CRL; }
    else if(rpro)
    {   stPO->owhat = SCTO_RPRO; }
    else if(!BAD_DOUB(stPO->flglo))
    {   stPO->owhat = SCTO_FLAG; }
    else
    {   stPO->owhat = SCTO_NONE; }
    if(!CheckSctOutOptionsI(stPO)) {
        ABORTLINE;
        CHECK_SCORETAB(stPO);
        return(FALSE);
    }
    /***
    *   Get the data file
    */
    if(!GetTableI(stPO->inname,stPO->rlab,stPO->clab,stPO->corn,stPO->do_skp,&stPO->tab)) {
        PROBLINE;
        printf("No table loaded!\n");
        ABORTLINE;
        CHECK_SCORETAB(stPO);
        return(FALSE);
    }
    c = GetTableColsI(stPO->tab,TRUE);
    r = GetTableRowsI(stPO->tab,TRUE);
    SetOutFormatting(stPO,fmtw,fmtp);
    if(!stPO->quiet) {
        printf("# Loaded table with %d rows, %d columns\n",r,c);
    }
    /***
    *   Transposing?    
    */
    if(stPO->do_tran) {
        printf("# Transposing...\n");
        HandleTransposeI(stPO->tab,&stPO->tab);
        c = GetTableColsI(stPO->tab,TRUE);
        r = GetTableRowsI(stPO->tab,TRUE);
        printf("# Transposed with %d rows, %d columns\n",r,c);
    }
    if(stPO->do_stru) {
        printf("# Sorting rows (ascending)...\n");
        HandleRowSortI(stPO->tab,1,FALSE);
    }
    if(stPO->do_strd) {
        printf("# Sorting rows (decending)...\n");
        HandleRowSortI(stPO->tab,-1,FALSE);
    }
    /***
    *   Set up extra stuff
    */
    if(!NO_S(stPO->outname)) {
        if(!(stPO->out = OpenUFilePF(stPO->outname,"w",NULL))) {
            ABORTLINE;
            CHECK_SCORETAB(stPO);
            return(FALSE);
        }
    }
    if(!SetUpSctAuxDataI(stPO)) {
        ABORTLINE;
        CHECK_SCORETAB(stPO);
        return(FALSE);
    }
    /***
    *   Masking out any cols/rows?
    */
    r = RestrictTableRangesI(stPO,stPO->tab);
    if(IS_BOG(r)) {
        ABORTLINE;
        CHECK_SCORETAB(stPO);
        return(FALSE);
    }
    else if(r>0) {
        if(!stPO->quiet) {
            r = GetTableRowsI(stPO->tab,TRUE);
            c = GetTableColsI(stPO->tab,TRUE);
            printf("# Restricted to %d total rows, %d total columns\n", r,c);
        }
    }
    /***
    *   Merging?     (GA cross over here too)
    */
    if(!NO_S(stPO->mergname)) {
        if(!HandleSctTableMergingI(stPO)) {
            ABORTLINE;
            CHECK_SCORETAB(stPO);
            return(FALSE);
        }
    }
    /***
    *   Scale / shift table values  (GA mutation here too)
    */
    AdjustTableValuesI(stPO,stPO->tab);
    HandleSmoothingI(stPO,stPO->tab);
    HandleNormalizationI(stPO,stPO->tab);
    /***
    *   Score transform of values?
    */
    if( (stPO->nscores>0) && (!HandleScoreTransformI(stPO,stPO->tab)) ) {
        ABORTLINE;
        CHECK_SCORETAB(stPO);
        return(FALSE);
    }
    /***
    *   Handle making symmetric, transpose
    */
    HandleSymmetrization(stPO,stPO->tab);
    /***
    *   Difference rows; new = odd-even
    */
    if(stPO->do_mrow) {
        if(!HandleScoreRowMergeI(stPO,stPO->tab)) {
            ABORTLINE;
            CHECK_SCORETAB(stPO);
            return(FALSE);
        }
    }
    /***
    *   Weak-first mapping?
    */
    if(stPO->do_wfm) {
        poolsize = GetTableColsI(stPO->tab,TRUE);
        if(!WeakFirstMappingI(stPO->tab,stPO->stab,TRUE,poolsize,20)) {
            printf("Weak-first mapping failed\n");
            ABORTLINE;
            CHECK_SCORETAB(stPO);
            return(FALSE);
        }
        if(stPO->do_dwfm) {
            AutoTableOutFormattingI(stPO->stab,TRUE,TRUE);
            DumpTable(stPO->stab,TRUE,TRUE,NULL);
        }
    }
    /***
    *   Partitioning?
    */
    if(stPO->psize > 0) {
        if(!HandleSetPartitioningI(stPO,stPO->tab)) {
            printf("Partitioning failed\n");
            ABORTLINE;
            CHECK_SCORETAB(stPO);
            return(FALSE);
        }
    }
    /***
    *   Report the story; If print format not user specified, auto set it
    */
    if(!stPO->usp_form) {
        AutoTableOutFormattingI(stPO->tab,TRUE,FALSE);
    }
    HandleSctOutput(stPO,stPO->out);
    /***
    *   All done
    */
    CHECK_SCORETAB(stPO);
    return(TRUE);
}
/**************************************************************************
*   Handle merging of odd-even rows
*/
int HandleScoreRowMergeI(SCORETAB *stPO,TABLE *tabPO)
{
    int r,c,mix;
    DOUB v1D,v2D,newD;
    char mixS[DEF_BS];

    mix = FigureMergeMixI(stPO, mixS);
    /***
    *   Go through cols, merging each pair of rows and masking out second val
    */
    for(c=0;c<tabPO->ncol;c++) {
        r = 0;
        while(r<tabPO->nrow) {
            GetTableValI(tabPO,r,c,&v1D);
            if( (r+1) < tabPO->nrow ) {
                GetTableValI(tabPO,r+1,c,&v2D);
                switch(mix)
                {
                    case MATH_ADD:  newD = v1D + v2D;   break;
                    case MATH_SUB:  newD = v1D - v2D;   break;
                    case MATH_MUL:  newD = v1D * v2D;   break;
                    case MATH_DIV:  
                        if(v2D == 0.0)
                        {   newD = 0.0;   }
                        else
                        {   newD = v1D / v2D;   }
                        break;
                    case MATH_MIN:  newD = MIN_NUM(v1D,v2D);    break;
                    case MATH_MAX:  newD = MAX_NUM(v1D,v2D);    break;
                    default:
                        printf("Bogus mix code: %d\n",mix);
                        ERR("HandleScoreRowMergeI","bad mixing code");
                        return(FALSE);
                }
                SetTableValI(tabPO,r,c,newD);
                tabPO->rmask[r+1] = FALSE;
            }
            r += 2;
        }
    }
    if(!stPO->quiet) {
        printf("# Rows merged: Odd and even %s\n",mixS);
    }
    return(TRUE);
}
/**************************************************************************
*   Handle normalization / quantation of table
*/
int HandleNormalizationI(SCORETAB *stPO,TABLE *tabPO)
{
    if(stPO->do_qcv > 1) {
        QuantileTableColsI(tabPO, stPO->do_qcv, TRUE, !stPO->quiet);
        if(!stPO->quiet) {
            printf("# Column values split into %d quantiles\n",stPO->do_qcv);
        }
    }
    else if(stPO->do_ncv) {
        NormTableI(tabPO,TRUE,TABLE_COL);
        if(!stPO->quiet) {
            printf("# Column values normalized\n");
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Handle smoothing of table
*/
int HandleSmoothingI(SCORETAB *stPO, TABLE *tabPO)
{
    if(stPO->do_srow > 0)
    {
        SmoothTableI(tabPO,TABLE_ROW,stPO->do_srow, FALSE);
        if(!stPO->quiet) {
            printf("# Row values smoothed with window +/- %d (%d wide)\n",
                stPO->do_srow, stPO->do_srow * 2 + 1);
        }
    }
    if(stPO->do_scol > 0)
    {
        SmoothTableI(tabPO,TABLE_COL,stPO->do_scol, FALSE);
        if(!stPO->quiet) {
            printf("# Col values smoothed with window +/- %d (%d wide)\n",
                stPO->do_scol, stPO->do_scol * 2 + 1);
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Handle output 
*/
void HandleSctOutput(SCORETAB *stPO,FILE *outPF)
{
    TABLE *tabPO,*stabPO;
    char bufS[DEF_BS];

    HAND_NFILE(outPF);
    tabPO = stPO->tab;
    stabPO = stPO->stab;
    switch(stPO->owhat)
    {
        case SCTO_FULL:
            fprintf(outPF,"# Full table\n");
            DumpTable(tabPO,FALSE,TRUE,outPF);
            break;
        case SCTO_RCL:
            fprintf(outPF,"# Full listing (row-col)\n");
            HandleSctOutList(stPO,tabPO,outPF);
            break;
        case SCTO_CRL:
            fprintf(outPF,"# Full listing (col-row)\n");
            HandleSctOutList(stPO,tabPO,outPF);
            break;
        case SCTO_RPRO:
            fprintf(outPF,"# Row product\n");
            HandleSctOutRows(stPO,tabPO,outPF);
            break;
        case SCTO_RSTAT:
            fprintf(outPF,"# Row statistics\n");
            fprintf(outPF,"# Row\tMin\tAve\tMax\tRange\tSd\tCV\tSum\n");
            HandleSctOutRows(stPO,tabPO,outPF);
            break;
        case SCTO_CSTAT:
            fprintf(outPF,"# Column statistics\n");
            fprintf(outPF,"# Col\tMin\tAve\tMax\tRange\tSd\tCV\tSum\n");
            HandleSctOutCols(stPO,tabPO,outPF);
            break;
        case SCTO_FSTAT:
            fprintf(outPF,"# Full table statistics\n");
            fprintf(outPF,"# Min\tAve\tMax\tRange\tSd\tCV\tSum\n");
            HandleSctOutStats(stPO,tabPO,TABLE_FULL,0,outPF);
            break;
        case SCTO_CINFO:
            if( (stPO->cinfo<1) || (stPO->cinfo>tabPO->ncol) ) {
                printf("Table has %d columns; Can't compare col %d to others\n",
                    tabPO->ncol,stPO->cinfo);
                return;
            }
            FillColOutputHeadLine(tabPO,stPO->cinfo,bufS);
            fprintf(outPF,"# Mutual information for %s\n",bufS);
            fprintf(outPF,"# Name\tM.I.\tN X N bins\tNorm M.I\n");
            HandleSctOutCols(stPO,tabPO,outPF);
            break;
        case SCTO_CCOR:
            if( (stPO->ccor<1) || (stPO->ccor>tabPO->ncol) ) {
                printf("Table has %d columns; Can't correlate col %d to others\n",
                    tabPO->ncol,stPO->ccor);
                return;
            }
            FillColOutputHeadLine(tabPO,stPO->ccor,bufS);
            fprintf(outPF,"# Correlations for %s\n",bufS);
            fprintf(outPF,"# Name\tPearson\tRank\n");
            HandleSctOutCols(stPO,tabPO,outPF);
            break;
        case SCTO_FCCM:
            fprintf(outPF,"# Full column correlation matrix\n");
            fprintf(outPF,"# Pearson = above diagonal (upper right)\n");
            fprintf(outPF,"# Rank    = below diagonal (lower left)\n");
            HandleSctColCors(stPO,tabPO,outPF);
            break;
        case SCTO_FLAG:
            if(stPO->do_not)
            {
                fprintf(outPF,"# Flagged values NOT in range %f to %f\n", 
                    stPO->flglo, stPO->flghi);
            }
            else
            {
                fprintf(outPF,"# Flagged values in range %f to %f\n", 
                    stPO->flglo, stPO->flghi);
            }
            HandleSctFlaggedOut(stPO,tabPO,outPF);
            break;
        case SCTO_NONE:
            break;
        default:
            printf("Uhhh....confused here, owhat = %d\n",stPO->owhat);
            ERR("HandleSctOutput","Bogus output code");
    }
    /***
    *   Weak first mapping done (and didn't dump matrix)?
    */
    if( (stPO->do_wfm) && (!stPO->do_dwfm) )
    {
        fprintf(outPF,"# Weak-first mapping results\n");
        HandleSctWfmOut(stPO,tabPO,stabPO,outPF);
    }
    return;
}
/**************************************************************************
*   
*/
void FillColOutputHeadLine(TABLE *tabPO, int col, char *bufS)
{
    char nameS[NSIZE];

    VALIDATE(tabPO,TABLE_ID);
    if(tabPO->clab) {
        GetTableColLabI(tabPO,col-1,nameS,-1);
        sprintf(bufS,"\"%s\" (column %d)",nameS,col);
    }
    else {
        sprintf(bufS,"column %d",col);
    }
    return;
}
/**************************************************************************
*   Handles list dump
*/
void HandleSctOutList(SCORETAB *stPO,TABLE *tabPO,FILE *outPF)
{
    int r,c;
    DOUB vD;
    char rnameS[NSIZE],cnameS[NSIZE];

    HAND_NFILE(outPF);
    if(stPO->owhat==SCTO_RCL) {
        for(r=0;r<tabPO->nrow;r++) {
            if(!tabPO->rmask[r])
            {   continue;   }
            if(tabPO->rlab)
            {   GetTableRowLabI(tabPO,r,rnameS,-1); }
            else
            {   sprintf(rnameS,"%05d",r); }
            for(c=0;c<tabPO->ncol;c++)
            {
                if(!tabPO->cmask[c])
                {   continue;   }
                if( (stPO->do_igd) && (c==r) )
                {   continue;   }
                if(tabPO->clab)
                {   GetTableColLabI(tabPO,c,cnameS,-1); }
                else
                {   sprintf(cnameS,"%05d",c); }
                GetTableValI(tabPO,r,c,&vD);
                PrintTabLabsValLineI(tabPO,rnameS,cnameS,vD,outPF);
            }
        }
    }
    else if(stPO->owhat==SCTO_CRL)
    {
        for(c=0;c<tabPO->ncol;c++)
        {
            if(!tabPO->cmask[c])
            {   continue;   }
            if(tabPO->clab)
            {   GetTableColLabI(tabPO,c,cnameS,-1); }
            else
            {   sprintf(cnameS,"%05d",c); }
            for(r=0;r<tabPO->nrow;r++)
            {
                if(!tabPO->rmask[r])
                {   continue;   }
                if( (stPO->do_igd) && (c==r) )
                {   continue;   }
                if(tabPO->rlab)
                {   GetTableRowLabI(tabPO,r,rnameS,-1); }
                else
                {   sprintf(rnameS,"%05d",r); }
                GetTableValI(tabPO,r,c,&vD);
                PrintTabLabsValLineI(tabPO,cnameS,rnameS,vD,outPF);
            }
        }
    }
    return;
}
/*************************************************************************/
void PrintTabLabsValLineI(TABLE *tabPO, char *fS, char *sS, DOUB vD, FILE *outPF)
{
    VALIDATE(tabPO,TABLE_ID);
    HAND_NFILE(outPF);
    if(fS) {
        fprintf(outPF,tabPO->prlform,fS);
    } 
    if(sS) {
        fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->prlform,sS);
    }
    fprintf(outPF,"%s",tabPO->pvsep); 
    fprintf(outPF,tabPO->pvform,vD); 
    fprintf(outPF,"\n");
    return;
}
/**************************************************************************
*   Handles output of row values for table
*/
void HandleSctOutRows(SCORETAB *stPO,TABLE *tabPO,FILE *outPF)
{
    int r,c;
    DOUB vD,tD;
    char nameS[NSIZE];

    HAND_NFILE(outPF);
    INIT_S(nameS);
    for(r=0;r<tabPO->nrow;r++) {
        if(!tabPO->rmask[r]) {
            continue;
        }
        if(tabPO->rlab) {
            GetTableRowLabI(tabPO,r,nameS,-1);
        }
        else {
            sprintf(nameS,"Row_%d",r+1);
        }
        fprintf(outPF,tabPO->prlform,nameS);
        if(stPO->owhat == SCTO_RSTAT) {
            HandleSctOutStats(stPO,tabPO,TABLE_ROW,r,outPF);
        }
        else if(stPO->owhat == SCTO_RPRO) {
            tD = 1.0;
            for(c=0;c<tabPO->ncol;c++) {
                if(!tabPO->cmask[c]) {
                    continue;
                }
                GetTableValI(tabPO,r,c,&vD);
                tD *= vD;
            }
            fprintf(outPF,"%s",tabPO->pvsep); 
            fprintf(outPF,tabPO->pvform,tD); 
            fprintf(outPF,"\n");
        }
        else {
            printf("Bogus process code: %d\n", stPO->owhat);
            ERR("HandleSctOutRows","Bad owhat code");
        }
    }
    return;
}
/**************************************************************************
*   Handles output of col values for table
*/
void HandleSctOutCols(SCORETAB *stPO,TABLE *tabPO,FILE *outPF)
{
    int c;
    char nameS[NSIZE];

    HAND_NFILE(outPF);
    INIT_S(nameS);
    for(c=0;c<tabPO->ncol;c++) {
        if(!tabPO->cmask[c]) {
            continue;
        }
        if(tabPO->clab) {
            GetTableColLabI(tabPO,c,nameS,-1);
        }
        else {
            sprintf(nameS,"Col_%d",c+1);
        }
        fprintf(outPF,"%s", nameS); 
        if(stPO->owhat == SCTO_CSTAT) {
            HandleSctOutStats(stPO,tabPO,TABLE_COL,c,outPF);
        }
        else if(stPO->owhat == SCTO_CCOR) {
            HandleSctSingColCors(stPO,tabPO,stPO->ccor-1,c,outPF);
        }
        else if(stPO->owhat == SCTO_CINFO) {
            HandleSctSingColInfo(stPO,tabPO,stPO->cinfo-1,c,outPF);
        }
        else {
            printf("Bogus process code: %d\n", stPO->owhat);
            ERR("HandleSctOutCols","Bad owhat code");
            return;
        }
    }
}
/**************************************************************************
*   Calc and report stats for rows, cols, or whole table
*/
void HandleSctOutStats(SCORETAB *stPO,TABLE *tabPO,int what,int which,
    FILE *outPF)
{
    int ok,n,mask;
    DOUB minD,maxD,avD,rgD,sdD,sumD,cvD;

    mask = TRUE;
    ok = FALSE;
    HAND_NFILE(outPF);
    if(what == TABLE_ROW) {
        ok = TableRowStatsI(tabPO,which,mask,&minD,&maxD,&avD,&sdD);
        /* Each row is n-col long; Need sep after row name */
        n = GetTableColsI(tabPO,mask);
        fprintf(outPF,"%s",tabPO->pvsep);
    }
    else if(what == TABLE_COL) {
        ok = TableColStatsI(tabPO,which,mask,&minD,&maxD,&avD,&sdD);
        /* Each col is n-row long; Need sep after col name */
        n = GetTableRowsI(tabPO,mask);
        fprintf(outPF,"%s",tabPO->pvsep);
    }
    else if(what == TABLE_FULL) {
        ok = TableStatsI(tabPO, 0, tabPO->nrow, 0, tabPO->ncol, mask,
            &minD,&maxD,&avD,&sdD);
        n = GetTableRowsI(tabPO,mask) * GetTableColsI(tabPO,mask);
    }
    else {
        printf("Sham value of what: %d\n",what);
        ERR("HandleSctOutStats","Bad col-vs-row code");
    }
    /***
    *   Print values or not
    */
    if(ok) {
        sumD = avD * RNUM(n);
        rgD = maxD - minD;
        if ( avD > 0 ) {
            cvD = sdD / avD;
        }
        else {
            cvD = 0.0;
        }
        /***
        *   Av, SD and CV use extended-precision, pform2
        */
        fprintf(outPF,tabPO->pvform,minD);   fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->pvform2,avD);   fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->pvform,maxD);   fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->pvform,rgD);    fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->pvform2,sdD);   fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->pvform2,cvD);   fprintf(outPF,"%s",tabPO->pvsep); 
        fprintf(outPF,tabPO->pvform,sumD);   
        fprintf(outPF,"\n");
    }
    else {
        fprintf(outPF,"NO STATS\n");
    }
}
/***********************************************************************
*   Mutual info for values in column c1 with c2
*/
void HandleSctSingColInfo(SCORETAB *stPO,TABLE *tabPO,int c1, int c2,
    FILE *outPF)
{
    DOUB vD;
    int nb1, nb2, ok;

    HAND_NFILE(outPF);
    /***
    *   sham; Should check these all work ok?
    */
    GetTableColValsI(tabPO,c1,stPO->tvals1,TRUE);   
    GetTableColValsI(tabPO,c2,stPO->tvals2,TRUE);   
    ok = NumlistMutInfoI(stPO->tvals1, stPO->tvals2, stPO->inf_maxhb, &vD, &nb1, &nb2);
    if(!ok) {
        vD = 0.0;
        nb1 = nb2 = 0;
    }
    fprintf(outPF,"%s",tabPO->pvsep); 
    fprintf(outPF,"%8.6f",vD); 
    fprintf(outPF,"%s",tabPO->pvsep); 
    fprintf(outPF,"%d %d",nb1,nb2);
    fprintf(outPF,"%s",tabPO->pvsep); 
    vD = (ok) ? vD / DNUM(nb2) : 0.0;
    fprintf(outPF,"%8.6f",vD); 
    fprintf(outPF,"\n"); 
}
/***********************************************************************
*   Correlate values in column c1 with c2
*/
void HandleSctSingColCors(SCORETAB *stPO,TABLE *tabPO,int c1, int c2,
    FILE *outPF)
{
    DOUB vD;

    HAND_NFILE(outPF);
    /***
    *   sham; Should check these all work ok?
    */
    GetTableColValsI(tabPO,c1,stPO->tvals1,TRUE);   
    GetTableColValsI(tabPO,c2,stPO->tvals2,TRUE);   
    /* Pearson */
    NumlistPairCorI(stPO->tvals1, stPO->tvals2, FALSE, &vD);
    fprintf(outPF,"%s",tabPO->pvsep); 
    fprintf(outPF,tabPO->pvform2,vD); 
    /* Rank correlation */
    NumlistPairCorI(stPO->tvals1, stPO->tvals2, TRUE, &vD);
    fprintf(outPF,"%s",tabPO->pvsep); 
    fprintf(outPF,tabPO->pvform2,vD); 
    fprintf(outPF,"\n"); 
}
/***********************************************************************
*   Do full col X col coorelations
*/
void HandleSctColCors(SCORETAB *stPO,TABLE *tabPO,FILE *outPF)
{
    int c1,c2,p,rank;
    DOUB vD;
    char nameS[NSIZE];

    HAND_NFILE(outPF);
    DumpTableColHeadings(tabPO,TRUE,outPF);
    for(c1=0; c1<tabPO->ncol; c1++)
    {
        if(!tabPO->cmask[c1]) {
            continue;
        }
        GetTableColValsI(tabPO,c1,stPO->tvals1,TRUE);   
        if(tabPO->clab) {
            GetTableColLabI(tabPO,c1,nameS,-1);
            fprintf(outPF,tabPO->prlform,nameS);
            fprintf(outPF,"%s",tabPO->pvsep); 
        }
        p = 0;
        for(c2=0; c2<tabPO->ncol; c2++) {
            if(!tabPO->cmask[c2]) {
                continue;
            }
            if(c1==c2) {
                vD = 1.0;
            }
            else {
                GetTableColValsI(tabPO,c2,stPO->tvals2,TRUE);   
                rank = (c2<c1) ? TRUE : FALSE;
                NumlistPairCorI(stPO->tvals1, stPO->tvals2, rank, &vD);
            }
            if(p>0) {
                fprintf(outPF,"%s",tabPO->pvsep); 
            }
            fprintf(outPF,tabPO->pvform2,vD); 
            p++;
        }
        fprintf(outPF,"\n"); 
    }
}
/***********************************************************************
*
*/
void HandleSctFlaggedOut(SCORETAB *stPO,TABLE *tabPO,FILE *outPF)
{
    int ok,r,c;
    DOUB vD;
    char rnameS[NSIZE],cnameS[NSIZE];

    /***
    *   Go through table transforming values
    */
    for(r=0;r<tabPO->nrow;r++)
    {
        if(!tabPO->rmask[r]) {
            continue;
        }
        for(c=0;c<tabPO->ncol;c++)
        {
            if(!tabPO->cmask[c]) {
                continue;
            }
            if( (stPO->do_igd) && (c==r) ) {    
                continue;   
            }
            GetTableValI(tabPO,r,c,&vD);
            ok = FALSE;
            if( (vD>=stPO->flglo) && (vD<=stPO->flghi) ) {
                ok = TRUE;
            }
            if(stPO->do_not) {
                ok = !ok;
            }
            if(!ok) {
                continue;
            }
            /***
            *   Get or concoct lables
            */
            if(tabPO->rlab) {   
                GetTableRowLabI(tabPO,r,rnameS,-1); 
            }
            else {  
                sprintf(rnameS,"Row-%03d",r); 
            }
            if(tabPO->clab) {   
                GetTableColLabI(tabPO,c,cnameS,-1); 
            }
            else {  
                sprintf(cnameS,"Col-%03d",c); 
            }
            PrintTabLabsValLineI(tabPO, rnameS, cnameS, vD, outPF);
        }
    }
}
/***********************************************************************
*
*/
void HandleSctWfmOut(SCORETAB *stPO,TABLE *vtabPO,TABLE *tabPO,FILE *outPF)
{
    int r,c,which,map;
    DOUB vD,mD;
    char rnameS[NSIZE], cnameS[NSIZE];

    HAND_NFILE(outPF);
    for(r=0;r<tabPO->nrow;r++)
    {
        if(!tabPO->rmask[r]) {
            continue;
        }
        GetTableRowLabI(tabPO,r,rnameS,-1);
        if(NO_S(rnameS)) {
            sprintf(rnameS,"[%d]",r);
        }
        map = FALSE;
        for(c=0;c<tabPO->ncol;c++)
        {
            if(!tabPO->cmask[c]) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            which = INT(vD);
            if(which < 0) {
                continue;
            }
            GetTableColLabI(tabPO,c,cnameS,-1);
            if(NO_S(cnameS)) {
                sprintf(cnameS,"[%d]",c);
            }
            /***
            *   Dump the mapping
            */
            map++;
            GetTableValI(vtabPO,r,c,&vD);
            GetTableValI(tabPO,r,c,&mD);
            fprintf(outPF,"# Row-Col %d-%d\t%1.0f\t%f\n",r+1,c+1,mD+1.0,vD);
            fprintf(outPF,"WFMapping %s\t%s\t%1.0f\n",rnameS,cnameS,mD+1.0);
        }
        /***
        *   Was anything found?
        */
        if(!map) {
            fprintf(outPF,"# NO MAPPING FOR ROW %d %s\n",r+1,rnameS);
        }
    }
}
