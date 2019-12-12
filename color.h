/*
* color.h
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


/*********** Color indices **************/
#define BLACK       0
#define WHITE       1
#define RED         2
#define REDORANGE   3
#define ORANGE      4
#define YELLOWORANGE 5
#define YELLOW      6
#define YELLOWGREEN 7
#define GREEN       8
#define CYANGREEN   9
#define CYAN        10
#define CYANBLUE    11
#define BLUE        12
#define PURPLE      13
#define MAGENTA     14
#define REDMAGENTA  15
#define LGREY       16
#define GREY        17
#define DGREY       18
#define LBROWN      19
#define BROWN       20
#define DBROWN      21
#define DGREEN      22
#define LBLUE       23
#define SEQCOL_A    24
#define SEQCOL_C    25
#define SEQCOL_G    26
#define SEQCOL_T    27
#define SEQCOL_N    28
#define U_COLOR     29


/************ Color structure ****************/
typedef struct COLOR
{
    int ID;
    REAL r,g,b;
    REAL h,s,v;
}COLOR;

#define MIN_COLOR           0.0025


/******************
*   Macros to set specific colors
*/
#define SET_RGB(co,x,y,z)   co->r = x; co->g = y; co->b = z
#define SET_WHITE(co)       SET_RGB(co,1.0,1.0,1.0);
#define SET_BLACK(co)       SET_RGB(co,0.0,0.0,0.0);
#define SET_GREY(co)        SET_RGB(co,0.5,0.5,0.5);
#define SET_RED(co)         SET_RGB(co,1.0,0.0,0.0);
#define SET_ORANGE(co)      SET_RGB(co,1.0,0.5,0.0);
#define SET_YELLOW(co)      SET_RGB(co,1.0,1.0,0.0);
#define SET_LIME(co)        SET_RGB(co,0.5,1.0,0.0);
#define SET_GREEN(co)       SET_RGB(co,0.0,1.0,0.0);
#define SET_CYAN(co)        SET_RGB(co,0.0,1.0,1.0);
#define SET_LBLUE(co)       SET_RGB(co,0.5,0.5,1.0);
#define SET_BLUE(co)        SET_RGB(co,0.0,0.0,1.0);
#define SET_PURPLE(co)      SET_RGB(co,0.5,0.0,1.0);
#define SET_MAGENTA(co)     SET_RGB(co,1.0,0.0,1.0);
#define SET_GREY1(co)       SET_RGB(co,0.1,0.1,0.1);
#define SET_GREY2(co)       SET_RGB(co,0.2,0.2,0.2);
#define SET_GREY3(co)       SET_RGB(co,0.3,0.3,0.3);
#define SET_GREY4(co)       SET_RGB(co,0.4,0.4,0.4);
#define SET_GREY5(co)       SET_RGB(co,0.5,0.5,0.5);
#define SET_GREY6(co)       SET_RGB(co,0.6,0.6,0.6);
#define SET_GREY7(co)       SET_RGB(co,0.7,0.7,0.7);
#define SET_GREY8(co)       SET_RGB(co,0.8,0.8,0.8);
#define SET_GREY9(co)       SET_RGB(co,0.9,0.9,0.9);

#define OK_ALPHA(al)    ( ((al) >= 0.0) && ((al) <= 1.0) )
#define BAD_RGB(co)     ( BAD_REAL((co)->r) && BAD_REAL((co)->g) && \
                            BAD_REAL((co)->b) )
#define BAD_HSV(co)     ( BAD_REAL((co)->h) && BAD_REAL((co)->s) && \
                            BAD_REAL((co)->v) )
#define BAD_COL(co)     (BAD_RGB(co) && BAD_HSV(co))


#define CHECK_COLOR(co) if(co){DestroyColorI(co); co=NULL;}

/*********************** ppp ********************
* C function listing generated by gen_prot
* Sun Jun  1 16:43:21 2003
*/
/****************************************************************
* color.c
*/
COLOR *CreateColorPO(void);
int DestroyColorI(COLOR *colPO);
void DumpColor(COLOR *colPO);
void InitColor(COLOR *colPO);
void BoundColor(COLOR *colPO);
int HSV_to_RGB_I(COLOR *colPO);
int RGB_to_HSV_I(COLOR *colPO);
void MixColors(COLOR *fcolPO, COLOR *scolPO, REAL alR, COLOR *colPO);
void DefColorMenu(void);
int ParseColorNameI(char *colS, int warn);

