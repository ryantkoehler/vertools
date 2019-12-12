/*
* plot_spec.c
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <gd.h>
#include "prim.h"
#include "color.h"
#include "plot_lib.h"

#define DB_IMAGE    if(DB[125])

/*************************************************************************
*   Set spectrum for correlation; Red high, White middle, Blue low
*/
int SetCorrelationRWBSpectrumI(IMAGEPLOT *plotPO)
{
    int nc,u_col,o_col;
    REAL hR,sR,vR,dhR,dsR,dvR;

    DB_IMAGE DB_PrI(">> SetCorrelationRWBSpectrumI\n");
    if(!GetImpageplotStatusI(plotPO))
    {
        DB_IMAGE DB_PrI("<< SetCorrelationRWBSpectrumI FALSE\n");
        return(FALSE);
    }
    /***
    *   First ramp; low to middle
    *   Blue to white
    */
    hR = 0.666;
    sR = 1.0;
    vR = 1.0;
    dhR = 0.0;
    dsR = -0.1;
    dvR = 0.0;
    nc = 10;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, TRUE);
    /***
    *   Second ramp; middle to high
    *   White to Red
    */
    hR = 1.0;
    sR = 0.0;
    vR = 1.0;
    dhR = 0.0;
    dsR = 0.1;
    dvR = 0.0;
    nc = 10;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, FALSE);
    /***
    *   Set Under/Over colors (one at a time)
    *   Under = Cyan-ish
    *   Over = Magenta-ish
    */
    u_col = SetImageplotColorI(plotPO, 0.0, 0.5, 1.0, BOGUS);
    o_col = SetImageplotColorI(plotPO, 1.0, 0.1, 1.0, BOGUS);
    SetImageplotSpecUnOvColorI(plotPO,u_col,o_col);
    /***
    *   Draw it since we've just set it
    */
    ImageplotSpectrumLegend(plotPO);
    DB_IMAGE DB_PrI("<< SetCorrelationRWBSpectrumI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Set spectrum for correlation; 
*   Red high, Black middle, Green low
*/
int SetCorrelationRKGSpectrumI(IMAGEPLOT *plotPO)
{
    int nc,u_col,o_col;
    REAL hR,sR,vR,dhR,dsR,dvR;

    DB_IMAGE DB_PrI(">> SetCorrelationRKGSpectrumI\n");
    if(!GetImpageplotStatusI(plotPO))
    {
        DB_IMAGE DB_PrI("<< SetCorrelationRKGSpectrumI FALSE\n");
        return(FALSE);
    }
    /***
    *   First ramp; low to middle
    *   Green to black
    */
    hR = 0.333;
    sR = 1.0;
    vR = 1.0;
    dhR = 0.0;
    dsR = 0.0;
    dvR = -0.1;
    nc = 10;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, TRUE);
    /***
    *   Second ramp; middle to high
    *   Black to Red
    */
    hR = 1.0;
    sR = 1.0;
    vR = 0.0;
    dhR = 0.0;
    dsR = 0.0;
    dvR = 0.1;
    nc = 10;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, FALSE);
    /***
    *   Set Under/Over colors (one at a time)
    *   Over = Magenta
    *   Under = Cyan
    */
    u_col = SetImageplotColorI(plotPO, 0.0, 1.0, 1.0, BOGUS);
    o_col = SetImageplotColorI(plotPO, 1.0, 0.0, 1.0, BOGUS);
    SetImageplotSpecUnOvColorI(plotPO,u_col,o_col);
    /***
    *   Draw it since we've just set it
    */
    ImageplotSpectrumLegend(plotPO);
    DB_IMAGE DB_PrI("<< SetCorrelationRKGSpectrumI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Set rainbow spectrum;
*   Low = cool, high = hot; Purple, Blue, Green, Yellow, Red
*/
int SetRainbowRHighSpectrumI(IMAGEPLOT *plotPO)
{
    int nc,u_col,o_col;
    REAL hR,sR,vR,dhR,dsR,dvR;

    DB_IMAGE DB_PrI(">> SetRainbowRHighSpectrumI\n");
    if(!GetImpageplotStatusI(plotPO))
    {
        DB_IMAGE DB_PrI("<< SetRainbowRHighSpectrumI FALSE\n");
        return(FALSE);
    }
    /***
    *   First ramp; Low towards middle
    *   purple to green
    */
    hR = 0.76;
    sR = 0.44;
    vR = 1.0;
    dhR = -0.045;
    dsR = 0.09;
    dvR = 0.00;
    nc = 7;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, TRUE);
    /***
    *   Second ramp; further up
    *   Darker green to yellow
    */
    hR = hR + RNUM(nc)*dhR;
    sR = 0.8;
    vR = 0.9;
    dhR = -0.095;
    dsR = 0.1;
    dvR = 0.05;
    nc = 3;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, FALSE);
    /***
    *   Third ramp; middle to high
    *   yellow to red
    */
    hR = hR + RNUM(nc)*dhR;
    sR = 1.0;
    vR = 0.95;
    dhR = -0.028;
    dsR = 0.08;
    dvR = 0.05;
    nc = 6;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, FALSE);
    /***
    *   Set Under/Over colors (one at a time)
    *   Under = light grey
    *   Over = white
    */
    u_col = SetImageplotColorI(plotPO, 0.8, 0.8, 0.8, BOGUS);
    o_col = SetImageplotColorI(plotPO, 1.0, 1.0, 1.0, BOGUS);
    SetImageplotSpecUnOvColorI(plotPO,u_col,o_col);
    /***
    *   Draw it since we've just set it
    */
    ImageplotSpectrumLegend(plotPO);
    DB_IMAGE DB_PrI("<< SetRainbowRHighSpectrumI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Set spectrum for correlation; 
*   Orange high, grey middle, Cyan low
*/
int SetOrangeGreyCyanSpectrumI(IMAGEPLOT *plotPO)
{
    int nc,u_col,o_col;
    REAL hR,sR,vR,dhR,dsR,dvR;

    DB_IMAGE DB_PrI(">> SetOrangeGreyCyanSpectrumI\n");
    if(!GetImpageplotStatusI(plotPO))
    {
        DB_IMAGE DB_PrI("<< SetOrangeGreyCyanSpectrumI FALSE\n");
        return(FALSE);
    }
    /***
    *   Number of steps down and up
    */
    nc = 4;
    /***
    *   First ramp; low to middle
    *   Cyan to grey
    */
    hR = 0.56;
    sR = 1.0;
    vR = 1.0;
    dhR = -0.02;
    dsR = -0.35;
    dvR = -0.12;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, TRUE);
    /***
    *   Second ramp; middle to high
    *   grey to Orange
    */
    hR = 0.15;
    sR = sR + (dsR * RNUM(nc)); 
    vR = vR + (dvR * RNUM(nc-1)); 
    dhR = -0.03;
    dsR = -dsR;
    dvR = -dvR;
    SetImageplotHSVSpectrumI(plotPO, nc, hR, sR, vR, dhR, dsR, dvR, FALSE);
    /***
    *   Set Under/Over colors (one at a time)
    *   Over = Magenta
    *   Under = Light grey
    */
    o_col = SetImageplotColorI(plotPO, 1.0, 0.0, 1.0, BOGUS);
    u_col = SetImageplotColorI(plotPO, 0.8, 0.8, 0.8, BOGUS);
    SetImageplotSpecUnOvColorI(plotPO,u_col,o_col);
    /***
    *   Draw it since we've just set it
    */
    ImageplotSpectrumLegend(plotPO);
    DB_IMAGE DB_PrI("<< SetOrangeGreyCyanSpectrumI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Parses R,G,B string into numbers
*/
int ParseRGBStringI(char *rgbS,REAL *rPR, REAL *gPR, REAL *bPR)
{
    char bufS[DEF_BS], *cPC;
    REAL rR,gR,bR;

    /***
    *   Init
    */
    if(rPR)
    {   *rPR = BAD_R;   }
    if(gPR)
    {   *gPR = BAD_R;   }
    if(bPR)
    {   *bPR = BAD_R;   }
    rR = gR = bR = BAD_R;
    /***
    *   Transfer single-token R,G,B into R G B
    */
    ReplaceChars(',', rgbS, ' ', bufS);
    /***
    *   Party through collection 
    */
    cPC = bufS;
    sscanf(cPC,"%lf",&rR);
    NEXT_WORD(cPC);
    sscanf(cPC,"%lf",&gR);
    NEXT_WORD(cPC);
    sscanf(cPC,"%lf",&bR);
    if( BAD_REAL(rR) || BAD_REAL(gR) || BAD_REAL(bR) )
    {
        return(FALSE);
    }
    /***
    *   Set values
    */
    if(rPR)
    {   *rPR = rR;  }
    if(gPR)
    {   *gPR = gR;  }
    if(bPR)
    {   *bPR = bR;  }
    return(TRUE);
}
