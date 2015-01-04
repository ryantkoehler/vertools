/*
* plotmat.c
*
* Copyright 2015 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include "numlist.h"
#include "wordlist.h"
#include "table.h"
#include "color.h"
#include "gd.h"     /* D graphics library; must preceed plot_lib.h */
#include "plot_lib.h"
#include "plotmat.h"

#define DB_TAB  if(DB[70])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(PlotMatI(argc,argv),NULL) ); }
/**************************************************************************/
void PlotMatUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage <in> ['-' for stdin] [...options]\n");
    printf("   <in>        Data file\n");
    printf("   -out XXX    Name of output (jpeg) file\n");
    printf("   -nrlab      Rows have NO lables (in input)\n");
    printf("   -nclab      Columns have NO lables (in input)\n");
    printf("   -ncorn      NO 'Corner' token between row / col labels\n");
    printf("   -sk         Skip missing / bogus rows; Default is abort\n");
    printf("   -pdim # #   Set plot X/Y dims to # # (matrix space)\n");
    printf("   -cdim # #   Set cell X/Y dims to # # (if < 1, auto fit)\n");
    printf("   -tm # -bm # Set top/bottom margins\n");
    printf("   -rm # -lm # Set right/left margins\n");
    printf("   -vcr # #    Set value color range from # to #\n");
    printf("   -acr        Auto color range from low to high\n");
    printf("   -sdc #      Set Standard Deviation color range +/- # sd\n");
    printf("   -nolab      No lables in output\n");
    printf("   -norl -nocl No row / col lables in output\n");
    printf("   -rlw #      Set row lable width to # (def = %d)\n",DEF_RLABW);
    printf("   -pcv        Print cell values\n");
    printf("   -pxv        Print only Extreme cell values (over / under)\n");
    printf("   -pfmt # #   Set print format # wide # precision (%%#.#f)\n");
    printf("   -fcol XXX   Use color XXX for printed cell value font\n");
    printf("   -igd        Ignore diagonal (leave blank)\n");
    printf("   -ngrid      No grid\n");
    printf("   -corcol     Use corelation color scheme (Blue White Red)\n");
    printf("   -corgbr     Use corelation color scheme (Green Black Red)\n");
    printf("   -ucol XXX   Use color XXX for under-range\n");
    printf("   -ocol XXX   Use color XXX for over-range\n");
    printf("   -colmen     List color menu\n");
    printf("   -dump       Dump out matrix\n");
}
/**************************************************************************
*   Main function
*/
int PlotMatI(int argc,char **argv)
{
    int nr,nc,rlab,clab,corn,nolab,dump,colmen,skip;
    char ocolS[DEF_BS],ucolS[DEF_BS];
    PLOTMAT *pmPO;

    pmPO = CreatePlotmatPO();
    INIT_S(ocolS); INIT_S(ucolS); 
    rlab = clab = corn = TRUE;
    nolab = skip = FALSE;
    dump = colmen = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -pdim I2 -dump B -vcr D2 -igd B -nrlab B -nclab B\
        -cdim I2 -tm I -bm I -rm I -lm I -corcol B -corgbr B\
        -nolab B -acr B -rlw I -ucol S -ocol S -norl B -nocl B -pcv B -pfmt I2\
        -ngrid B -colmen B -fcol S -sdc D -pxv B -sk B -ncorn B",
        pmPO->inname, pmPO->outname, &pmPO->xdim,&pmPO->ydim,
        &dump, &pmPO->lo,&pmPO->hi, &pmPO->do_igd, &rlab, &clab,
        &pmPO->xc,&pmPO->yc, &pmPO->tm, &pmPO->bm, &pmPO->rm, &pmPO->lm,
        &pmPO->do_corcol, &pmPO->do_corgbr, &nolab,
        &pmPO->do_acr, &pmPO->rowlabw, ucolS, ocolS,
        &pmPO->do_rowlab, &pmPO->do_collab,
        &pmPO->do_pcv, &pmPO->pftw,&pmPO->pftp, &pmPO->pgrid, &colmen,
        pmPO->fcolor, &pmPO->sd_col, &pmPO->do_pxv, &skip, &corn,
        (int *)NULL))
    {
        PlotMatUse();
        CHECK_PLOTMAT(pmPO);
        return(FALSE);
    }
    /***
    *   Only want color menu
    */
    if(colmen) {
        DefColorMenu();
        CHECK_PLOTMAT(pmPO);
        return(TRUE);
    }
    /***
    *   Load table 
    */
    if(!GetTableI(pmPO->inname,rlab,clab,corn,skip,&pmPO->tab)) {
        PROBLINE;
        printf("No table = no plot!\n");
        CHECK_PLOTMAT(pmPO);
        return(FALSE);
    }
    AutoTableOutFormattingI(pmPO->tab,TRUE,TRUE);
    nr = GetTableRowsI(pmPO->tab,FALSE);
    nc = GetTableColsI(pmPO->tab,FALSE);
    printf("# Have %d X %d table (R X C)\n",nr,nc);
    if(dump) {
        DumpTable(pmPO->tab,TRUE,FALSE,NULL);
    }
    /***
    *   Check / set output name. If stdin, name is "-" so change that
    */
    if(NO_S(pmPO->outname)) {
        GetFilePartsI(pmPO->inname,NULL,pmPO->outname,NULL);
        if(!strcmp(pmPO->inname,"-")) {
            sprintf(pmPO->outname,"%s",DEF_OUTNAME);
        }
        strcat(pmPO->outname,".jpg");
    }
    /***
    *   Check options and set vars
    */
    if(nolab) {
        pmPO->do_rowlab = pmPO->do_collab = FALSE;
    }
    if(!CheckPlotmatOptionsI(pmPO)) {
        PROBLINE;
        CHECK_PLOTMAT(pmPO);
        return(FALSE);
    }
    SetPlotmatFormatString(pmPO);
    /***
    *   Set up jpg pix and set up initial stuff
    */
    FigurePlotDims(pmPO,nr,nc);
    ReportPlotSettings(pmPO);
    if(!SetUpPlottingI(pmPO,ucolS,ocolS)) {
        PROBLINE;
        printf("Problem setting up plotting vars\n");
        CHECK_PLOTMAT(pmPO);
        return(FALSE);
    }
    /***
    *   Draw the picture and save
    */
    PlotTableOnPage(pmPO);
    SetImageplotOutQuality(pmPO->plot,95);
    SaveImageplotImageI(pmPO->plot);
    /***
    *   Close the whip
    */
    CHECK_PLOTMAT(pmPO);
    return(TRUE);
}
/*************************************************************************
*
*/
PLOTMAT *CreatePlotmatPO()
{
    PLOTMAT *pmPO;

    if( ! (pmPO = (PLOTMAT *)ALLOC(1,sizeof(PLOTMAT))) ) {
        return(NULL);
    }
    pmPO->ID = PLOTMAT_ID;
    InitPlotmat(pmPO);
    return(pmPO);
}
/*************************************************************************
*   Free structure
*/
int DestroyPlotmatI(PLOTMAT *pmPO)
{
    VALIDATE(pmPO,PLOTMAT_ID);
    CHECK_TABLE(pmPO->tab);
    CHECK_IMAGEPLOT(pmPO->plot);
    FREE(pmPO);
    return(TRUE);
}
/*************************************************************************
*   Initialize 
*/
void InitPlotmat(PLOTMAT *pmPO)
{
    VALIDATE(pmPO,PLOTMAT_ID);
    
    INIT_S(pmPO->inname);
    INIT_S(pmPO->outname);
    pmPO->hi = 1.0;
    pmPO->lo = 0.0;
    pmPO->xdim = DEF_DIMS;
    pmPO->ydim = DEF_DIMS;
    pmPO->tm = DEF_TMAR;
    pmPO->bm = DEF_BMAR;
    pmPO->lm = DEF_LMAR;
    pmPO->rm = DEF_RMAR;
    pmPO->do_rowlab = TRUE;
    pmPO->do_collab = TRUE;
    pmPO->do_pcv = pmPO->do_pxv = FALSE;
    pmPO->pftw = pmPO->pftp = BOGUS;
    pmPO->pgrid = TRUE;
    pmPO->do_acr = FALSE;
    pmPO->sd_col = 0.0;
    pmPO->do_corcol = FALSE;
    pmPO->do_corgbr = FALSE;
    pmPO->rowlabw = BOGUS;
    pmPO->collabh = DEF_CLABH;
    pmPO->ucolor = BOGUS;
    pmPO->ocolor = BOGUS;
    INIT_S(pmPO->fcolor);
    pmPO->pcvcolor = BOGUS;
    return;
}
/*************************************************************************
*
*/
void SetPlotmatFormatString(PLOTMAT *pmPO)
{
    int rw;
    TABLE *tabPO;

    tabPO = pmPO->tab;
    VALIDATE(tabPO,TABLE_ID);
    /***
    *   Print cell value format; If not specified, get auto from table
    */
    if(!IS_BOG(pmPO->pftw)) {
        FloatFormatString(pmPO->pftw,pmPO->pftp,pmPO->pnform);
    }
    else {
        GetTablePrintformI(tabPO,pmPO->pnform,NULL,NULL,NULL,NULL);
    }
    /***
    *   Row lable width
    */
    if(IS_BOG(pmPO->rowlabw)) {
        WordlistAutoFormatStringI(tabPO->rlabs, &rw, NULL);
        pmPO->rowlabw = rw * DEF_RLPPC; 
    }
    return;
}
/*************************************************************************
*
*/
int CheckPlotmatOptionsI(PLOTMAT *pmPO)
{
    return(TRUE);
}
/*************************************************************************
*   Get working cell dimensions
*/
void FigurePlotDims(PLOTMAT *pmPO,int nr,int nc)
{
    /***
    *   Explict X cell dimensions?
    */
    if(pmPO->xdim <= 0) {
        pmPO->xdim = DEF_DIMS;
    }
    if(pmPO->xc <= 0) {
        pmPO->xc = ROUND(pmPO->xdim / nc);
        LIMIT_NUM(pmPO->xc,1,TOO_BIG);
    }
    pmPO->xdim = pmPO->xc * nc + pmPO->lm + pmPO->rm;
    if(pmPO->do_rowlab) {
        pmPO->xdim += pmPO->rowlabw;
    }
    /***
    *   Explict Y cell dimensions?
    */
    if(pmPO->ydim <= 0) {
        pmPO->ydim = DEF_DIMS;
    }
    if(pmPO->yc <= 0) {
        pmPO->yc = ROUND((pmPO->ydim - pmPO->tm - pmPO->bm - DEF_MSPACE)/ nr);
        LIMIT_NUM(pmPO->yc,1,TOO_BIG);
    }
    pmPO->ydim = pmPO->yc * nr + pmPO->tm + pmPO->bm + DEF_MSPACE;
    if(pmPO->do_collab) {
        pmPO->ydim += pmPO->collabh;
    }
    return;
}
/************************************************************************/
void ReportPlotSettings(PLOTMAT *pmPO)
{
    printf("#  Dimensions: %4d by %4d\n",pmPO->xdim,pmPO->ydim);
    printf("#  Cells:      %4d by %4d\n",pmPO->xc,pmPO->yc);
    return;
}
/*************************************************************************
*   Create graphics object and set dims, colors, etc.
*/
int SetUpPlottingI(PLOTMAT *pmPO, char *ucolS, char *ocolS)
{
    DOUB loD,hiD,avD,sdD;

    if(! (pmPO->plot = CreateImageplotPO(pmPO->xdim,pmPO->ydim)) ) {
        PROBLINE;
        printf("Failed to allocate image object\n");
        return(FALSE);
    }
    SetImageplotFilename(pmPO->plot,pmPO->outname);
    /***
    *   If auto color range or SD coloring, need stats 
    */
    loD = pmPO->lo;
    hiD = pmPO->hi;
    if( pmPO->do_acr || (pmPO->sd_col > 0.0) ) {
        /***
        *   Call with no row or col restrictioins, no masking
        */
        TableStatsI(pmPO->tab, -1, -1, -1, -1, FALSE, &loD, &hiD, &avD, &sdD);
        /***
        *   Very small adjust (shrink) range so auto over-under values show 
        */
        if(pmPO->do_acr) {
            loD += ASCORE_OFF;
            hiD -= ASCORE_OFF;
        }
        if( pmPO->sd_col > 0.0) {
            loD = avD - (sdD * pmPO->sd_col);
            hiD = avD + (sdD * pmPO->sd_col);
            printf("# Plot range %.4f to %.4f, mean = %.4f, sd = %.4f\n",
                loD,hiD,avD,sdD);
        }
    }
    /***
    *   Save used values and set in plot
    */
    pmPO->lo = loD;
    pmPO->hi = hiD;
    SetImageplotLoHiScore(pmPO->plot,loD,hiD);
    /***
    *   Get Over / Under color indices or fail
    */
    if(!NO_S(ucolS)) {
        pmPO->ucolor = ParseColorNameI(ucolS,TRUE);
        if(IS_BOG(pmPO->ucolor)) {
            return(FALSE);
        }
    }
    if(!NO_S(ocolS)) {
        pmPO->ocolor = ParseColorNameI(ocolS,TRUE);
        if(IS_BOG(pmPO->ocolor)) {
            return(FALSE);
        }
    }
    if(!NO_S(pmPO->fcolor)) {
        pmPO->pcvcolor = ParseColorNameI(pmPO->fcolor,TRUE);
        if(IS_BOG(pmPO->pcvcolor)) {
            return(FALSE);
        }
    }
    SetLegendPosition(pmPO);
    SetupPlotPage(pmPO,pmPO->plot);
    return(TRUE);
}
/*************************************************************************
*   Set the postion of the color legend
*/
void SetLegendPosition(PLOTMAT *pmPO)
{
    int x,y,legw;

    /***
    *   Figure out where to put the legend
    */
    x = pmPO->lm + DEF_LEGTEX;
    y = pmPO->ydim - DEF_LEGY;
    if(pmPO->do_rowlab)
    {
        x += pmPO->rowlabw;
    }
    /***
    *   Legend width reduced if narrowed plot
    */
    legw = DEF_LEGW;
    if( (x+legw) >= (pmPO->xdim - 10) )
    {
        legw = pmPO->xdim - 10;
    }
/*
printf("XXX x=%d y=%d leg=%d\n",x,y,legw);
*/
    SetImageplotSpecLegDims(pmPO->plot,x,y,legw,DEF_LEGH,DEF_LEGTEX);
}
/*************************************************************************
*   Initialize page and draw spectrum
*/
void SetupPlotPage(PLOTMAT *pmPO, IMAGEPLOT *plotPO)
{
    ImageplotColorBoxI(plotPO,0,0,plotPO->x,plotPO->y,WHITE);
    if(pmPO->do_corcol)
    {
        SetCorrelationRWBSpectrumI(plotPO);
    }
    else if(pmPO->do_corgbr)
    {
        SetCorrelationRKGSpectrumI(plotPO);
    }
    else
    {
        SetRainbowRHighSpectrumI(plotPO);
    }
    /***
    *   Call here to walk on over / under colors set via above spectrums
    */
    SetImageplotSpecUnOvColorI(pmPO->plot,pmPO->ucolor,pmPO->ocolor);
    ImageplotSpectrumLegend(plotPO);
}
/*************************************************************************
*   Loop through matrix plotting the whip
*/
void PlotTableOnPage(PLOTMAT *pmPO)
{
    int r,c,w,h,x,y,color,p,nr,nc;
    DOUB rD;
    char nameS[NSIZE];
    
    w = pmPO->xc;
    h = pmPO->yc;
    nr = GetTableRowsI(pmPO->tab,FALSE);
    nc = GetTableColsI(pmPO->tab,FALSE);
    /***
    *   If col header
    */
    if(pmPO->do_collab) {
        y = pmPO->tm;
        x = pmPO->lm;
        if(pmPO->do_rowlab) {
            x += pmPO->rowlabw;
        }
        for(c=0;c<nc;c++) {
            GetTableColLabI(pmPO->tab,c,nameS,NSIZE-1);
            ImageplotTextStringI(pmPO->plot,x,y,nameS,GFONT_MEDBOLD,BLACK);
            x += w;
        }
    }
    /***
    *   Loop down the rows
    */
    for(r=0;r<nr;r++) {
        y = pmPO->tm + r * h;
        if(pmPO->do_collab) {
            y += pmPO->collabh;
        }
        x = pmPO->lm;
        if(pmPO->do_rowlab) {
            GetTableRowLabI(pmPO->tab,r,nameS,NSIZE-1);
            ImageplotTextStringI(pmPO->plot,x,y,nameS,GFONT_MEDBOLD,BLACK);
            x += pmPO->rowlabw;
        }
        /***
        *   Each cell in row
        */
        for(c=0;c<nc;c++) {
            GetTableValI(pmPO->tab,r,c,&rD);
            color = ImageplotColorForScoreI(pmPO->plot,rD);
            if( (r==c) && (pmPO->do_igd) ) {
                p = FALSE;
            }
            else {
                p = TRUE;
            }
            if(p) {
                if(pmPO->pgrid) {
                    if( (r==(nr-1)) || (c==(nc-1)) ) {
                        ImageplotColorBoxI(pmPO->plot,x-1,y-1,w+2,h+2,BLACK);
                    }
                    else {
                        ImageplotColorBoxI(pmPO->plot,x-1,y-1,w+1,h+1,BLACK);
                    }
                }
                ImageplotColorBoxI(pmPO->plot,x,y,w,h,color);
            }
            /***
            *   Printing cell values?
            */
            if(pmPO->do_pcv) {
                PlotCellValue(pmPO, x, y, rD, color);
            }
            else if(pmPO->do_pxv && ((rD < pmPO->lo) || (rD > pmPO->hi)) ) {
                PlotCellValue(pmPO, x, y, rD, color);
            }
            x += w;
        }
    }
}
/******************************************************************************
*   XXX SHAM with hard-coded color bounds!!!
*/
void PlotCellValue(PLOTMAT *pmPO, int x, int y, DOUB rD, int color)
{
    int fcol,fx,fy;
    char bufS[DEF_BS];

    sprintf(bufS,pmPO->pnform,rD);
    fcol = WHITE;
    /***
    *   Correlation Red-White-Blue coloring 
    */
    if(pmPO->do_corcol) {
        if( (color>=36) && (color<=44) ) {
            fcol = BLACK;
        }
    }
    /***
    *   Green-Black-Red coloring 
    */
    else if(pmPO->do_corgbr) {
        if( (color<=31) || (color>=47) ) {
            fcol = BLACK;
        }
    }
    /***
    *   Assume rainbow coloring
    */
    else {
        if( (color>=34) && (color<=41) ) {
            fcol = BLACK;
        }
    }
    /***
    *   If explictly specified, use that
    */
    if(!NO_S(pmPO->fcolor)) {
        fcol = pmPO->pcvcolor;
    }
    /***
    *   XXX val lable placment
    */
    fx = x + 3;
    fy = y + 3;
    /* mageplotTextStringI(pmPO->plot,x+3,y+10,bufS,GFONT_MEDBOLD,fcol); */
    ImageplotTextStringI(pmPO->plot,fx,fy,bufS,GFONT_MEDBOLD,fcol);
}
