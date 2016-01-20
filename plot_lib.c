/*
* plot_lib.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gd.h>     
#include <gdfontt.h>
#include <gdfonts.h>
#include <gdfontmb.h>
#include <gdfontl.h>
#include <gdfontg.h>
#include "prim.h"
#include "color.h"
#include "plot_lib.h"

#define DB_IMAGE    if(DB[125])
#define DB_DRAW if(DB[126])

/**************************************************************************
*   Creates gifplot data struct
*/
IMAGEPLOT *CreateImageplotPO(int x,int y)
{
    IMAGEPLOT *plotPO;

    DB_IMAGE DB_PrI(">> CreateImageplotPO\n");
    if(!(plotPO = (IMAGEPLOT *)ALLOC(1,sizeof(IMAGEPLOT))))
    {
        DB_IMAGE DB_PrI("<< CreateImageplotPO NULL\n");
        return(NULL);
    }
    plotPO->ID = IMAGEPLOT_ID;
    /***
    *   Init
    */
    InitImageplotPars(plotPO);
    if( (x>0) && (y>0) )
    {
        if(!InitImageplotImageI(plotPO,x,y))
        {
            CHECK_IMAGEPLOT(plotPO);
            return(NULL);
        }
    }
    DB_IMAGE DB_PrI("<< CreateImageplotPO %p\n",plotPO);
    return(plotPO);
}
/**************************************************************************
*   Frees up gifplot data structure 
*/
int DestroyImageplotI(IMAGEPLOT *plotPO)
{
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(plotPO->gd)
    {
        gdImageDestroy(plotPO->gd);
    }
    FREE(plotPO);
    return(TRUE);
}
/**************************************************************************
*   Initialize parameters for gifplot data struct
*/
void InitImageplotPars(IMAGEPLOT *plotPO)
{
    int i;

    DB_IMAGE DB_PrI(">> InitImageplotPars\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    plotPO->havegd = FALSE;
    plotPO->havedefcol = FALSE;
    plotPO->havespeccol = FALSE;
    /***
    *   file and dim settings
    */
    INIT_S(plotPO->outname);
    plotPO->oqual = -1;
    plotPO->x = DEF_XDIM;
    plotPO->y = DEF_YDIM;
    /***
    *   Colors
    */
    plotPO->bgcol = plotPO->fgcol = BOGUS;
    plotPO->ncol = 0;
    for(i=0;i<NUM_COL;i++)
    {
        plotPO->colors[i] = BOGUS;
    }
    plotPO->nccol = 0;
    plotPO->ccoff = BOGUS;
    plotPO->nscol = 0;
    plotPO->scoff = BOGUS;
    plotPO->sucol = plotPO->socol = BOGUS;
    /***
    *   Spectrum 
    */
    plotPO->specx = DEF_SPECX;
    plotPO->specy = DEF_SPECY;
    plotPO->specw = DEF_SPECW;
    plotPO->spech = DEF_SPECH;
    plotPO->spectex = BOGUS;
    SetImageplotLoHiScore(plotPO,DEF_SPECLO,DEF_SPECHI);
    strcpy(plotPO->pnform,DEF_IMPNFORM_S);
    DB_IMAGE DB_PrI("<< InitImageplotPars\n");
}
/**************************************************************************
*   Initializes a image drawing area
*/ 
int InitImageplotImageI(IMAGEPLOT *plotPO,int x,int y)
{
    DB_IMAGE DB_PrI("\n>> InitImageplotImageI DIMS %d %d\n",x,y);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(plotPO->havegd)
    {
        WARNLINE;
        printf("InitImageplotImageI called with image already existing\n");
        return(FALSE);
    }
    /***
    *   Create new structure to X Y dims
    */
    BOG_CHECK(plotPO->gd);
    LIMIT_NUM(x,1,TOO_BIG);
    LIMIT_NUM(y,1,TOO_BIG);
    plotPO->gd = gdImageCreate(x,y);
    if(!plotPO->gd)
    {
        PROBLINE;
        printf("Failed to allocate a image structure (%d by %d)\n",x,y);
        return(FALSE);
    }
    plotPO->havegd = TRUE;
    plotPO->x = x;
    plotPO->y = y;
    SetImageplotDefColorsI(plotPO);
    DB_IMAGE DB_PrI("<< InitImageplotImageI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Close drawing area, writing to file if we've got an image
*/
int SaveImageplotImageI(IMAGEPLOT *plotPO)
{
    FILE *fPF;

    DB_IMAGE DB_PrI(">> SaveImageplotImageI\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(NO_S(plotPO->outname))
    {
        DB_IMAGE DB_PrI("<< SaveImageplotImageI no filename FALSE\n");
        return(FALSE);
    }
    /***
    *   Open file for writing 
    */
    if(!(fPF = OpenUFilePF(plotPO->outname,"wb",NULL)))
    {
        PROBLINE;
        printf("Can't write to image file\n");
        DB_IMAGE DB_PrI("<< SaveImageplotImageI bad file FALSE\n");
        return(FALSE);
    }
    /***
    *   What format to write?
    */
    switch(plotPO->otype)
    {
        case PLOTOUT_PNG:
            DB_IMAGE DB_PrI("+ Have file, calling gdImageGif\n");
            gdImagePng(plotPO->gd,fPF);
            break;
        case PLOTOUT_JPG:
        default:
            DB_IMAGE DB_PrI("+ Have file, calling gdImageJpeg\n");
            gdImageJpeg(plotPO->gd,fPF,plotPO->oqual);
            break;
    }
    CHECK_NFILE(fPF,plotPO->outname);
    DB_IMAGE DB_PrI("<< SaveImageplotImageI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Is there a writable image yet? 
*/
int GetImpageplotStatusI(IMAGEPLOT *plotPO)
{
    VALIDATE(plotPO,IMAGEPLOT_ID);
    return(plotPO->havegd);
}
/**************************************************************************
*   Spectrum legend placement 
*/
void SetImageplotSpecLegDims(IMAGEPLOT *plotPO,int x,int y,int w,int h,int tex)
{
    DB_IMAGE DB_PrI(">> SetImageplotSpecLegDims x=%d y=%d w=%d h=%d tex=%d\n",
        x,y,w,h,tex);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    plotPO->specx = x + tex;
    plotPO->specy = y;
    plotPO->specw = w;
    plotPO->spech = h;
    plotPO->spectex = tex;
    DB_IMAGE DB_PrI("<< SetImageplotSpecLegDims\n");
}
/**************************************************************************
*   Sets output file name into structure
*/
void SetImageplotFilename(IMAGEPLOT *plotPO,char *fileS)
{
    DB_IMAGE DB_PrI(">> SetImageplotFilename |%s|\n",fileS);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    strcpy(plotPO->outname,fileS);  
    DB_IMAGE DB_PrI("<< SetImageplotFilename\n");
}
/**************************************************************************
*   Sets output file format
*/
void SetImageplotOutFormat(IMAGEPLOT *plotPO,int type)
{
    DB_IMAGE DB_PrI(">> SetImageplotOutFormat %d\n",type);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    switch(type)
    {
        case PLOTOUT_PNG:   
            plotPO->otype = PLOTOUT_PNG;    
            break;
        case PLOTOUT_JPG:   
        default:
            plotPO->otype = PLOTOUT_JPG;    
            break;
    }
    DB_IMAGE DB_PrI("<< SetImageplotOutFormat\n");
}
/**************************************************************************
*   Sets output file quality
*/
void SetImageplotOutQuality(IMAGEPLOT *plotPO,int qual)
{
    DB_IMAGE DB_PrI(">> SetImageplotOutQuality %d\n",qual);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    LIMIT_NUM(qual,-1,100);
    plotPO->oqual = qual;   
    DB_IMAGE DB_PrI("<< SetImageplotOutQuality\n");
}
/**************************************************************************
*   Set default colors
*/
int SetImageplotDefColorsI(IMAGEPLOT *plotPO)
{
    int i;

    DB_IMAGE DB_PrI(">> SetImageplotDefColorsI\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    /***
    *   Check if don't yet have gif structure or already set def colors
    */
    if( (!plotPO->havegd) || (plotPO->havedefcol) )
    {   
        DB_IMAGE DB_PrI("<< FALSE; Gif %d Defcol %d\n",
            plotPO->havegd, plotPO->havedefcol);
        return(FALSE);  
    }
    /***
    *   Bunch of defined color indices; Set color count to U_COLOR 
    */
    SetImageplotColorI(plotPO,  0.0,    0.0,    0.0,    BLACK);
    SetImageplotColorI(plotPO,  1.0,    1.0,    1.0,    WHITE);
    SetImageplotColorI(plotPO,  1.0,    0.0,    0.0,    RED);
    SetImageplotColorI(plotPO,  1.0,    0.25,   0.0,    REDORANGE);
    SetImageplotColorI(plotPO,  1.0,    0.5,    0.0,    ORANGE);
    SetImageplotColorI(plotPO,  1.0,    0.8,    0.0,    YELLOWORANGE);
    SetImageplotColorI(plotPO,  1.0,    1.0,    0.0,    YELLOW);
    SetImageplotColorI(plotPO,  0.5,    1.0,    0.0,    YELLOWGREEN);
    SetImageplotColorI(plotPO,  0.0,    1.0,    0.0,    GREEN);
    SetImageplotColorI(plotPO,  0.0,    0.9,    0.6,    CYANGREEN);
    SetImageplotColorI(plotPO,  0.0,    0.95,   1.0,    CYAN);
    SetImageplotColorI(plotPO,  0.0,    0.5,    1.0,    CYANBLUE);
    SetImageplotColorI(plotPO,  0.0,    0.0,    1.0,    BLUE);
    SetImageplotColorI(plotPO,  0.6,    0.0,    1.0,    PURPLE);
    SetImageplotColorI(plotPO,  1.0,    0.0,    1.0,    MAGENTA);
    SetImageplotColorI(plotPO,  1.0,    0.0,    0.5,    REDMAGENTA);
    SetImageplotColorI(plotPO,  0.8,    0.8,    0.8,    LGREY);
    SetImageplotColorI(plotPO,  0.5,    0.5,    0.5,    GREY);
    SetImageplotColorI(plotPO,  0.3,    0.3,    0.3,    DGREY);
    SetImageplotColorI(plotPO,  0.8,    0.6,    0.2,    LBROWN);
    SetImageplotColorI(plotPO,  0.6,    0.5,    0.2,    BROWN);
    SetImageplotColorI(plotPO,  0.4,    0.3,    0.2,    DBROWN);
    SetImageplotColorI(plotPO,  0.0,    0.7,    0.0,    DGREEN);
    SetImageplotColorI(plotPO,  0.6,    0.6,    1.0,    LBLUE);
    SetImageplotColorI(plotPO,  0.0,    0.9,    0.0,    SEQCOL_A);
    SetImageplotColorI(plotPO,  0.2,    0.2,    1.0,    SEQCOL_C);
    SetImageplotColorI(plotPO,  0.0,    0.0,    0.0,    SEQCOL_G);
    SetImageplotColorI(plotPO,  1.0,    0.0,    0.0,    SEQCOL_T);
    SetImageplotColorI(plotPO,  0.6,    0.6,    0.0,    SEQCOL_N);
    /***
    *   Check ok?
    */
    for(i=0;i<plotPO->ncol;i++)
    {
        DB_IMAGE DB_PrI("+ Col[%d] = %d\n",i,plotPO->colors[i]);
    }
    /***
    *   Set for/background and flag that defs are set
    */
    plotPO->fgcol = plotPO->colors[BLACK];
    plotPO->bgcol = plotPO->colors[WHITE];
    plotPO->havedefcol = TRUE;
    DB_IMAGE DB_PrI("<< SetImageplotDefColorsI ncol=%d TRUE\n",plotPO->ncol);
    return(TRUE);
}
/**************************************************************************
*   Set spectrum "under" and "over" colors
*/
int SetImageplotSpecUnOvColorI(IMAGEPLOT *plotPO,int uc,int oc)
{
    DB_IMAGE DB_PrI(">> SetImageplotSpecUnOvColorI uc=%d, oc=%d\n", uc,oc);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< SetImageplotSpecUnOvColorI no image FALSE\n");
        return(FALSE);  
    }
    if( (uc>=0) && (uc<plotPO->ncol) )
    {
        plotPO->sucol = uc;
        DB_IMAGE DB_PrI("+ ucol set to %d\n",plotPO->sucol);
    }
    if( (oc>=0) && (oc<plotPO->ncol) )
    {
        plotPO->socol = oc;
        DB_IMAGE DB_PrI("+ ocol set to %d\n",plotPO->socol);
    }
    CheckImageplotSpecUnOvCols(plotPO);
    DB_IMAGE DB_PrI("<< SetImageplotSpecUnOvColorI Un=%d Ov=%d\n",
            plotPO->sucol, plotPO->socol);
    return(TRUE);
}
/**************************************************************************/
void CheckImageplotSpecUnOvCols(IMAGEPLOT *plotPO)
{
    DB_IMAGE DB_PrI(">> CheckImageplotSpecUnOvCols\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< CheckImageplotSpecUnOvCols no image\n");
        return; 
    }
    if(IS_BOG(plotPO->sucol))
    {
        plotPO->sucol = plotPO->colors[WHITE];
        DB_IMAGE DB_PrI("+ ucol set to %d\n",plotPO->sucol);
    }
    if(IS_BOG(plotPO->socol))
    {
        plotPO->socol = plotPO->colors[WHITE];
        DB_IMAGE DB_PrI("+ ocol set to %d\n",plotPO->sucol);
    }
    DB_IMAGE DB_PrI("<< CheckImageplotSpecUnOvCols\n");
}
/**************************************************************************
*   Sets color into plot structure
*   If index is < 0, assign color to next available index and return this
*/
int SetImageplotColorI(IMAGEPLOT *plotPO,REAL rR,REAL gR,REAL bR,int ind)
{
    int color;

    DB_IMAGE DB_PrI(">> SetImageplotColorI rgb=%f %f %f, ind=%d\n",rR,gR,bR,ind);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< SetImageplotColorI FALSE\n");
        return(BOGUS);  
    }
    /***
    *   If index is too big, problem
    */
    if(ind>=NUM_COL)
    {
        PROBLINE;
        printf("Bad color index %d; can't allocate color\n",ind);
        return(BOGUS);
    }
    /***
    *   If negative, assign to next color
    */
    if(ind<0)
    {
        ind = plotPO->ncol;
    }
    /***
    *   Get color index into gd gif structure
    */
    LIMIT_NUM(rR,0.0,1.0); LIMIT_NUM(gR,0.0,1.0); LIMIT_NUM(bR,0.0,1.0);
    color = AllocGdColorI(plotPO,INT(rR*255.0),INT(gR*255.0),INT(bR*255.0));
    DB_IMAGE DB_PrI("+ RGB = %f %f %f = color %d\n", rR, gR, bR,color);
    if(color < 0)
    {
        PROBLINE;
        printf("Failed to allocate color %f %f %f\n", rR, gR, bR);
        return(BOGUS);
    } 
    plotPO->colors[ind] = color;
    plotPO->ncol += 1;
    DB_IMAGE DB_PrI("<< SetImageplotColorI ncol=%d [%d]=%d\n",plotPO->ncol,
        ind, plotPO->colors[ind]);
    return(ind);
}
/***************************************************************************
*
*/
int AllocGdColorI(IMAGEPLOT *plotPO,int r,int g,int b)
{
    int col,n;

    DB_IMAGE DB_PrI(">> AllocGdColorI %d %d %d ",r,g,b);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< AllocGdColorI BOGUS\n");
        return(BOGUS);  
    }
    col = BOGUS;
    n = 0;
    while(n<MAX_COL_TRIES)
    {
        col = gdImageColorAllocate(plotPO->gd,r,g,b);
        if(col>=0)
            break;
        n++;
    }
    DB_IMAGE DB_PrI("<< AllocGdColorI %d\n",col);
    return(col);
}
/**************************************************************************
*   Set low and hi score values for comp_seq spectrum 
*/
void SetImageplotLoHiScore(IMAGEPLOT *plotPO, REAL loR, REAL hiR)
{
    DB_IMAGE DB_PrI(">> SetImageplotLoHiScore lo=%f hi=%f\n",loR,hiR);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if((hiR - loR) <= 0.0)
    {
        PROBLINE;
        printf("Bad score range for spectrum: low=%f hi=%f\n",loR,hiR);
        printf("  No change\n");
        return;
    }
    plotPO->losc = loR;
    plotPO->hisc = hiR;
    plotPO->dsc = hiR - loR;
    DB_IMAGE DB_PrI("<< SetImageplotLoHiScore dif=%f\n",plotPO->dsc);
}
/*************************************************************************/
int GetImageplotLoHiScoreI(IMAGEPLOT *plotPO,REAL *loPR, REAL *hiPR)
{
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(loPR) {
        *loPR = plotPO->losc;
    }
    if(hiPR) {
        *hiPR = plotPO->hisc;
    }
    return(TRUE);
}
/*************************************************************************/
int GetImageplotForBackColorsI(IMAGEPLOT *plotPO,int *forPI,int *backPI)
{
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        return(FALSE);
    }
    if(forPI)
    {
        *forPI = plotPO->fgcol;
    }
    if(backPI)
    {
        *backPI = plotPO->bgcol;
    }
    return(TRUE);
}
/*************************************************************************
*   Set ncol colors in a spectrum from starting HSV and dH, dS, dV
*   Where HSV = hue (1=red,0.33=green,0.66=blue), saturation (1=full), 
*       value (1=full) and dH, dS, dV are delta HSV for each step 
*/
int SetImageplotHSVSpectrumI(IMAGEPLOT *plotPO, int ncol, REAL hR, REAL sR,
    REAL vR, REAL dhR, REAL dsR,REAL dvR,int new)
{
    int i,off;
    COLOR colO,*colPO;

    DB_IMAGE DB_PrI(">> SetImageplotHSVSpectrumI n=%d new=%d\n", ncol,new);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< SetImageplotHSVSpectrumI no image FALSE\n");
        return(FALSE);  
    }
    /***
    *   Local color sham
    */
    colPO = &colO; colPO->ID = COLOR_ID;
    colPO->h = hR;
    colPO->s = sR;
    colPO->v = vR;
    DB_IMAGE DB_PrI("+ HSV %f %f %f\n",hR,sR,vR);
    /***
    *   Walk on old guys or extend spectrum collection
    */
    if(new)
    {
        off = plotPO->ncol;
        plotPO->scoff = off;
        plotPO->nscol = ncol;
        DB_IMAGE DB_PrI("+ NEW set off to ncol = %d\n",off);
    }
    else 
    {
        off = plotPO->scoff + plotPO->nscol;
        plotPO->nscol += ncol;
        DB_IMAGE DB_PrI("+ off adjust to %d, nscol=%d\n",off,plotPO->nscol);
    }
    /***
    *   Ramp colors
    */
    for(i=0;i<ncol;i++)
    {
        BoundColor(colPO);
        HSV_to_RGB_I(colPO);
        if(!SetImageplotColorI(plotPO,colPO->r,colPO->g,colPO->b,i+off))
        {
            PROBLINE;
            printf("Failed to set color %d\n",i+off);
            return(FALSE);
        }
        colPO->h += dhR;
        colPO->s += dsR;
        colPO->v += dvR;
    }
    plotPO->havespeccol = TRUE;
    CheckImageplotSpecUnOvCols(plotPO);
    DB_IMAGE DB_PrI("<< SetImageplotHSVSpectrumI TRUE\n");
    return(TRUE);
}
/*****************************************************************************
*   Returns the color associated with score (based on spectrum settings)
*/
int ImageplotColorForScoreI(IMAGEPLOT *plotPO, REAL scR)
{
    int color;
    REAL cindR;

    DB_IMAGE {
        DB_PrI(">> ImageplotColorForScoreI sc=%f (lo=%f hi=%f)\n",
            scR,plotPO->losc,plotPO->hisc);
        DB_PrI("+ sucol=%d socol=%d\n", plotPO->sucol, plotPO->socol);
    }
    VALIDATE(plotPO,IMAGEPLOT_ID);
    /***
    *   Check for under / over case
    */
    if( (scR<plotPO->losc) && (!IS_BOG(plotPO->sucol)) )
    {
        color = plotPO->sucol;
        DB_IMAGE DB_PrI("<< %d (under score color)\n",color);
        return(color);
    }
    if( (scR>plotPO->hisc) && (!IS_BOG(plotPO->socol)) )
    {
        color = plotPO->socol;
        DB_IMAGE DB_PrI("<< %d (over score color)\n",color);
        return(color);
    }
    /***
    *   Bound score and get color from index
    */
    LIMIT_NUM(scR,plotPO->losc,plotPO->hisc);
    DB_IMAGE DB_PrI("+ limited to %f\n",scR);
    cindR = (scR - plotPO->losc) / plotPO->dsc;
    DB_IMAGE DB_PrI("+ frac = %f\n",cindR);
    color = ROUND(cindR * RNUM(plotPO->nscol));
    DB_IMAGE DB_PrI("+ colorR = %d\n",color);
    LIMIT_NUM(color,0,(plotPO->nscol-1));
    DB_IMAGE DB_PrI("+ limit to = %d\n",color);
    color += plotPO->scoff;
    DB_IMAGE DB_PrI("+ offset to = %d\n",color);
    color = plotPO->colors[color];
    DB_IMAGE DB_PrI("<< ImageplotColorForScoreI %d\n",color);
    return(color);
}
/****************************************************************************
*   Draw & lable spectrum
*/
void ImageplotSpectrumLegend(IMAGEPLOT *plotPO)
{
    char upS[DEF_BS],dnS[DEF_BS];

    DB_IMAGE DB_PrI(">> ImageplotSpectrumLegend\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< ImageplotSpectrumLegend no image\n");
        return;
    }
    sprintf(upS,plotPO->pnform,plotPO->losc);
    sprintf(dnS,plotPO->pnform,plotPO->hisc);
    ImageplotSpecLegendColor(plotPO);
    ImageplotSpecLegendTex(plotPO,upS,dnS);
    DB_IMAGE DB_PrI("<< ImageplotSpectrumLegend\n");
}
/****************************************************************************
*   Lable spectrum with supplied text strings
*/
void ImageplotSpecLegendTex(IMAGEPLOT *plotPO, char *sS,char *eS)
{
    int color,sx,sy;

    DB_IMAGE DB_PrI(">> ImageplotSpecLegendTex %p %p\n",sS,eS);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< ImageplotSpecLegendTex no image\n");
        return;
    }
    color = plotPO->colors[BLACK];
    sx = plotPO->specx - plotPO->spectex;
    sy = plotPO->specy;
    if(sS)
    {
        DB_IMAGE DB_PrI("+ S|%s| sx=%d sy=%d\n",sS,sx,sy);
        ImageplotTextStringI(plotPO,sx,sy,sS,GFONT_MEDBOLD,color,FALSE);
    }
    if(eS)
    {
        sx = sx + plotPO->specw + 2 * plotPO->spectex;
        DB_IMAGE DB_PrI("+ E|%s| sx=%d sy=%d\n",sS,sx,sy);
        ImageplotTextStringI(plotPO,sx,sy,eS,GFONT_MEDBOLD,color,FALSE);
    }
    DB_IMAGE DB_PrI("<< ImageplotSpecLegendTex\n");
}
/****************************************************************************
*   Draw color spectrum legend
*/
void ImageplotSpecLegendColor(IMAGEPLOT *plotPO)
{
    int i,sx,sy,w,h;
    REAL xR,dxR;

    DB_IMAGE DB_PrI(">> ImageplotSpecLegendColor\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< no image\n");
        return;
    }
    /***
    *   Start coords and dims
    */
    xR = RNUM(plotPO->specx);
    dxR = RNUM(plotPO->specw) / RNUM(plotPO->nscol);
    sx = INT(xR);
    sy = plotPO->specy;
    w = INT(dxR);
    h = plotPO->spech;
    DB_IMAGE DB_PrI("+ xR=%f dxR=%f x=%d y=%d w=%d h=%d\n",xR,dxR,sx,sy,w,h);
    /***
    *   First under-score color
    */
    if(!IS_BOG(plotPO->sucol))
    {
        ImageplotColorBoxI(plotPO,sx-1,sy-1,w+1,h+2,BLACK);
        ImageplotColorBoxI(plotPO,sx,sy,w,h,plotPO->sucol);
    }
    xR += dxR;
    /***
    *   Run through each color
    */
    for(i=0;i<plotPO->nscol;i++)
    {
        sx = ROUND(xR);
        /***
        *   Black bound for boxes add up to rim around whole
        */
        if(i==0)
        {
            ImageplotColorBoxI(plotPO,sx-1,sy-1,w+2,h+2,BLACK);
        }
        else
        {
            ImageplotColorBoxI(plotPO,sx,sy-1,w+2,h+2,BLACK);
        }
        ImageplotColorBoxI(plotPO,sx,sy,w,h,plotPO->colors[i+plotPO->scoff]);
        xR += dxR;
    }
    /***
    *   Last = over-score color
    */
    if(!IS_BOG(plotPO->socol))
    {
        sx = ROUND(xR);
        ImageplotColorBoxI(plotPO,sx-1,sy-1,w+2,h+2,BLACK);
        ImageplotColorBoxI(plotPO,sx,sy,w,h,plotPO->socol);
    }
    DB_IMAGE DB_PrI("<< ImageplotSpecLegendColor\n");
}
/****************************************************************************
*   Write text below spectrum 
*/
void ImageplotSubSpecTex(IMAGEPLOT *plotPO,char *sS)
{
    int sx,sy,color;

    DB_IMAGE DB_PrI(">> ImageplotSubSpecTex\n");
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) 
    {
        DB_IMAGE DB_PrI("<< no image\n");
        return;
    }
    color = plotPO->colors[BLACK];
    sx = plotPO->specx;
    sy = plotPO->specy + plotPO->spech + 2;
    ImageplotTextStringI(plotPO,sx,sy,sS,GFONT_MEDBOLD,color,FALSE);
}
/****************************************************************************
*   11/30/10 RTK; Had to sham up string char pointer to unsigned char for
*       library calls; Mac gcc issued pointer missmatch warning.
*   3/23/15 RTK; Add rotate argument and option to call regular or "Up" 
*/
int ImageplotTextStringI(IMAGEPLOT *plotPO, int x, int y, char *texPC, int font,
    int color, int rotate)
{
    unsigned char *texS;

    DB_DRAW DB_PrI(">> ImageplotTextStringI x=%d y=%d font=%d color=%d rotate=%d\n",
        x,y,font,color,rotate);
    texS = (unsigned char *) texPC;
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd) {
        DB_DRAW DB_PrI("<< ImageplotTextStringI no image\n");
        return(FALSE);
    }
    if(color<0) {
        color = plotPO->colors[BLACK];
        DB_DRAW DB_PrI("+ Black = col=%d\n",color);
    }
    switch(font)
    {
        case GFONT_TINY:
            if(rotate) {
                gdImageStringUp(plotPO->gd,gdFontTiny,x,y,texS,color);
            }
            else {
                gdImageString(plotPO->gd,gdFontTiny,x,y,texS,color);
            }
            break;
        case GFONT_SMALL:
            if(rotate) {
                gdImageStringUp(plotPO->gd,gdFontSmall,x,y,texS,color);
            }
            else {
                gdImageString(plotPO->gd,gdFontSmall,x,y,texS,color);
            }
            break;
        case GFONT_MEDBOLD:
            if(rotate) {
                gdImageStringUp(plotPO->gd,gdFontMediumBold,x,y,texS,color);
            }
            else {
                gdImageString(plotPO->gd,gdFontMediumBold,x,y,texS,color);
            }
            break;
        case GFONT_LARGE:
            if(rotate) {
                gdImageStringUp(plotPO->gd,gdFontLarge,x,y,texS,color);
            }
            else {
                gdImageString(plotPO->gd,gdFontLarge,x,y,texS,color);
            }
            break;
        case GFONT_GIANT:
            if(rotate) {
                gdImageStringUp(plotPO->gd,gdFontGiant,x,y,texS,color);
            }
            else {
                gdImageString(plotPO->gd,gdFontGiant,x,y,texS,color);
            }
            break;
        default:
            if(rotate) {
                gdImageStringUp(plotPO->gd,gdFontMediumBold,x,y,texS,color);
            }
            else {
                gdImageString(plotPO->gd,gdFontMediumBold,x,y,texS,color);
            }
            break;
    }
    DB_DRAW DB_PrI("<< ImageplotTextStringI TRUE\n");
    return(TRUE);
}
/************************************************************************
*
*/
int ImageplotColorBoxI(IMAGEPLOT *plotPO, int x,int y,int w,int h,int color)
{
    int cind;

    DB_DRAW DB_PrI(">> ImageplotColorBoxI x=%d y=%d w=%d h=%d col=%d\n",
        x,y,w,h,color);
    VALIDATE(plotPO,IMAGEPLOT_ID);
    if(!plotPO->havegd)
    {
        DB_DRAW DB_PrI("<< ImageplotColorBoxI no image FALSE\n");
        return(FALSE);
    }
    cind = plotPO->colors[color];
    if(IS_BOG(cind))
    {
        DB_DRAW DB_PrI("<< ImageplotColorBoxI bad color FALSE\n");
        return(FALSE);
    }
    gdImageFilledRectangle(plotPO->gd,x,y,x+w,y+h,color);
    DB_DRAW DB_PrI("<< ImageplotColorBoxI TRUE\n");
    return(TRUE);
}
