/*
* color.c
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
#include <math.h>
#include <ctype.h>
#include "prim.h"
#include "color.h"

#define DB_COLOR if(DB[45])


/**************************************************************************/
COLOR *CreateColorPO()
{
    COLOR *colPO;
    
    colPO = (COLOR *)ALLOC(1,sizeof(COLOR));
    colPO->ID = COLOR_ID;
    return(colPO);
}
/**************************************************************************/
int DestroyColorI(COLOR *colPO)
{
    VALIDATE(colPO,COLOR_ID);
    FREE(colPO);
    return(TRUE);
}
/**************************************************************************/
void DumpColor(COLOR *colPO)
{
    VALIDATE(colPO,COLOR_ID);
    PrintI("RGB = %3.3f %3.3f %3.3f,  HSV = %3.3f %3.3f %3.3f\n",
        colPO->r, colPO->g, colPO->b, colPO->h, colPO->s, colPO->v);
}
/**************************************************************************/
void InitColor(COLOR *colPO)
{
    VALIDATE(colPO,COLOR_ID);
    colPO->r = colPO->g = colPO->b = 1.0;
    colPO->h = 0.0;
    colPO->s = colPO->v = 1.0;
}
/**************************************************************************
*   Make's sure HSV / RGB values are in range 0-1
*/
void BoundColor(COLOR *colPO)
{
    VALIDATE(colPO,COLOR_ID);
    if(colPO->h > 1.0)
    {   colPO->s -= FLOOR(colPO->s);    }
    if(colPO->h < 0.0)
    {   colPO->h -= FLOOR(colPO->h);    }
    LIMIT_NUM(colPO->s,0.0,1.0);
    LIMIT_NUM(colPO->v,0.0,1.0);
    LIMIT_NUM(colPO->r,0.0,1.0);
    LIMIT_NUM(colPO->g,0.0,1.0);
    LIMIT_NUM(colPO->b,0.0,1.0);
}
/**************************************************************************/
int HSV_to_RGB_I(COLOR *colPO)
{
    REAL fR,kR,mR,nR,hueR;
    int i;

    VALIDATE(colPO,COLOR_ID);
    DB_COLOR printf(">> HSV_to_RGB_I, HSV = %4.3f %4.3f %4.3f\n",
        colPO->h, colPO->s, colPO->v);
    /***
    *   check for the achromatic case 
    */
    if(EQUAL(colPO->s,0.0,MIN_COLOR))
    {
        colPO->r = colPO->g = colPO->b = colPO->v;
        DB_COLOR printf("+ sat too small, = value %4.3f\n",colPO->v);
    }
    else
    {
        hueR = colPO->h * 6.0;
        DB_COLOR printf("+ setting hue %4.3f\n",hueR);
        if(EQUAL(hueR,6.0,MIN_COLOR))
        {
            hueR = 0.0;
            DB_COLOR printf("+ capping hue %4.3f\n",hueR);
        }
        i = INT(hueR);
        fR = hueR - RNUM(i);
        mR = colPO->v * (1.0 - colPO->s);
        nR = colPO->v * (1.0 - colPO->s * fR);
        kR = colPO->v * (1.0 - colPO->s * (1.0 - fR));
        DB_COLOR printf("+ i = %d, f %4.3f, m %4.3f, n %4.3f, k %4.3f\n",i,
            fR,mR,nR,kR);
        switch(i)
        {
            case 0: SET_RGB(colPO,colPO->v,kR,mR);  break;
            case 1: SET_RGB(colPO,nR,colPO->v,mR);  break;
            case 2: SET_RGB(colPO,mR,colPO->v,kR);  break;
            case 3: SET_RGB(colPO,mR,nR,colPO->v);  break;
            case 4: SET_RGB(colPO,kR,mR,colPO->v);  break;
            case 5: SET_RGB(colPO,colPO->v,mR,nR);  break;
        }
    }
    DB_COLOR printf("<< HSV_to_RGB_I, RGB = %4.3f %4.3f %4.3f\n",
        colPO->r, colPO->g, colPO->b);
    return(TRUE);
}
/**************************************************************************/
int RGB_to_HSV_I(COLOR *colPO)
{
    REAL tempR,rR,gR,bR;

    VALIDATE(colPO,COLOR_ID);
    DB_COLOR printf(">> RGB_to_HSV_I, RGB = %4.3f %4.3f %4.3f\n",
        colPO->r,colPO->g,colPO->b);
    /***    
    *   determine the value 
    */
    colPO->v = MAX_NUM(colPO->r,colPO->g);
    colPO->v = MAX_NUM(colPO->b,colPO->v);
    DB_COLOR printf("+ value = %4.3f\n",colPO->v);
    /***
    *   determine saturation 
    */
    tempR = MIN_NUM(colPO->r,colPO->g);
    tempR = MIN_NUM(colPO->b,tempR);
    DB_COLOR printf("+ temp = %4.3f, ",tempR);
    if(EQUAL(colPO->v,0.0,MIN_COLOR))
    {
        colPO->s = 0.0;
        DB_COLOR printf("too small, sat = %4.3f\n",colPO->s);
    }
    else
    {
        colPO->s = (colPO->v - tempR) / colPO->v;
        DB_COLOR printf("sat = %4.3f\n",colPO->s);
    }
    /***
    *   determine the hue 
    */
    if(EQUAL(colPO->s,0.0,MIN_COLOR))
    {
        /* rbitrary sham ---> set to zero */
        colPO->h = 0.0;
        DB_COLOR printf("+ sat too small, hue = %4.3f\n",colPO->h);
    }
    else
    {
        rR = (colPO->v - colPO->r) / (colPO->v - tempR);
        gR = (colPO->v - colPO->g) / (colPO->v - tempR);
        bR = (colPO->v - colPO->b) / (colPO->v - tempR);
        DB_COLOR printf("+ sat ok, rgb = %4.3f, %4.3f, %4.3f\n",rR,gR,bR);
        /* he color is between yellow and magenta */
        if(colPO->r == colPO->v)
        {
            colPO->h = bR - gR;
            DB_COLOR printf("+ Y-M, hue = %4.3f\n",colPO->h);
        }
        /* he color is between cyan and yellow */
        else if(colPO->g == colPO->v)
        {
            colPO->h = 2.0 + rR - bR;
            DB_COLOR printf("+ C-Y, hue = %4.3f\n",colPO->h);
        }
        /* he color is between magenta and cyan */
        else
        {
            colPO->h = 4.0 + gR - rR;
            DB_COLOR printf("+ M-C, hue = %4.3f\n",colPO->h);
        }
        /* onvert to range 0-1 */
        colPO->h *= 0.166666;
        DB_COLOR printf("+ scaling %4.3f, ",colPO->h);
        /* revent negative value */
        if(colPO->h < 0.0)
        {
            DB_COLOR printf("+ bumping %4.3f, ",colPO->h);
            colPO->h += 1.0;
        }
    }
    DB_COLOR printf("<< RGB_to_HSV_I, HSV = %4.3f %4.3f %4.3f\n",
        colPO->h,colPO->s,colPO->v);
    return(TRUE);
}       
/**********************************************************************
*   Mix two input colors, fcol and scol, into another color, col, using
*   alpha to determine how much of scol to use: 0 = none, 1 = all
*
*   Does both RGB and HSV
*/
void MixColors(COLOR *fcolPO, COLOR *scolPO, REAL alR, COLOR *colPO)
{
    VALIDATE(fcolPO, COLOR_ID);
    VALIDATE(scolPO, COLOR_ID);
    VALIDATE(colPO, COLOR_ID);
    if(!OK_ALPHA(alR))
    {
        printf("BAD ALPHA: %f\n",alR);
        ERR("MixColors","Bad alpha value");
    }
    DB_COLOR 
    {
        DB_PrI(">> MixColors, alR %3.3f\n",alR);
        DB_PrI("+ fcol "); DumpColor(fcolPO);
        DB_PrI("+ scol "); DumpColor(scolPO);
    }
    colPO->r = fcolPO->r + ((scolPO->r - fcolPO->r) * alR);
    colPO->g = fcolPO->g + ((scolPO->g - fcolPO->g) * alR);
    colPO->b = fcolPO->b + ((scolPO->b - fcolPO->b) * alR);
    colPO->h = fcolPO->h + ((scolPO->h - fcolPO->h) * alR);
    colPO->s = fcolPO->s + ((scolPO->s - fcolPO->s) * alR);
    colPO->v = fcolPO->v + ((scolPO->v - fcolPO->v) * alR);
    DB_COLOR 
    {
        DB_PrI("+ mix "); DumpColor(colPO);
        DB_PrI("<< MixColors\n");
    }
}
/**********************************************************************
*
*/
void DefColorMenu()
{
    BARLINE;
    printf("Defined color list:\n");
    BARLINE;
    printf("    BK  Black\n");
    printf("    W   White\n");
    printf("    RE  Red\n");
    printf("    RO  Red Orange\n");
    printf("    O   Orange\n");
    printf("    YO  Yellow Orange\n");
    printf("    YE  Yellow\n");
    printf("    YG  Yellow Green\n");
    printf("    GR  Green\n");
    printf("    GD  Dark Green\n");
    printf("    CG  Cyan Green\n");
    printf("    CY  Cyan\n");
    printf("    CB  Cyan Blue\n");
    printf("    BL  Blue\n");
    printf("    BI  Light Blue\n");
    printf("    P   Purple\n");
    printf("    MA  Magenta\n");
    printf("    RM  Red Magenta\n");
    printf("    LG  Light Grey\n");
    printf("    MG  Medium grey\n");
    printf("    DG  Dark dgrey\n");
    printf("    LB  Light Brown\n");
    printf("    MB  Medium Brown\n");
    printf("    DB  Dark Brown\n");
}
/**********************************************************************
*   Parse color string
*/
int ParseColorNameI(char *colS, int warn)
{
    int col;

    col = BOGUS;
    switch(UPPER(colS[0]))
    {
        case 'B':
            switch(UPPER(colS[1]))
            {
                case 'K':   col = BLACK;    break;
                case 'L':   col = BLUE; break;
                case 'I':   col = LBLUE;    break;
            }
            break;
        case 'C':
            switch(UPPER(colS[1]))
            {
                case 'G':   col = CYANGREEN;    break;
                case 'Y':   col = CYAN; break;
                case 'B':   col = CYANBLUE; break;
            }
            break;
        case 'D':
            switch(UPPER(colS[1]))
            {
                case 'B':   col = DBROWN;   break;
                case 'G':   col = DGREY;    break;
            }
            break;
        case 'G':
            switch(UPPER(colS[1]))
            {
                case 'R':   col = GREEN;    break;
                case 'D':   col = DGREEN;   break;
            }
            break;
        case 'L':
            switch(UPPER(colS[1]))
            {
                case 'B':   col = LBROWN;   break;
                case 'G':   col = LGREY;    break;
            }
            break;
        case 'M':
            switch(UPPER(colS[1]))
            {
                case 'A':   col = MAGENTA;  break;
                case 'G':   col = GREY; break;
                case 'B':   col = BROWN;    break;
            }
            break;
        case 'O':
            col = ORANGE;
            break;
        case 'P':
            col = PURPLE;
            break;
        case 'R':
            switch(UPPER(colS[1]))
            {
                case 'E':   col = RED;  break;
                case 'M':   col = REDMAGENTA;   break;
                case 'O':   col = REDORANGE;    break;
            }
            break;
        case 'W':
            col = WHITE;
            break;
        case 'Y':
            switch(UPPER(colS[1]))
            {
                case 'O':   col = YELLOWORANGE; break;
                case 'E':   col = YELLOW;   break;
                case 'G':   col = YELLOWGREEN;  break;
            }
    }
    if( (IS_BOG(col)) && (warn) )
    {
        printf("Unrecognized color: %s\n",colS);
        DefColorMenu();
    }
    return(col);
}
