/*
* num.c
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
#include <ctype.h>
#include <math.h>
#define __MAIN__
#include "prim.h"

#define VERSION_S   "NumberTool Version 1.24"

#define OP_ADD  1
#define OP_SUB  2
#define OP_MUL  3
#define OP_DIV  4
#define OP_MOD  5
#define OP_LOG  8
#define OP_LN   9
#define OP_LG2  10
#define OP_EXP  11
#define OP_ALOG 12
#define OP_PW2  13
#define OP_SIN  14
#define OP_COS  15
#define OP_TAN  16
#define OP_SQRT 17
#define OP_MAX  18
#define OP_MIN  19
#define OP_FMT  20

#define MIN_SIZE        1e-20
#define MIN_LOG         1e-99

int main(int argc, char **argv);
void NumUse(void);
int NumI(int argc, char **argv);
void HandleFormatString(REAL rR,FILE *outPF);
int ParseNumOpI(char *opS, int nargs);

/**************************************************************************/
int main(int argc, char **argv)
{   Init(argc,argv); exit( AllDoneI(NumI(argc,argv),NULL) );    }
/************************************************************************/
void NumUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: \n");
    printf("   -add X Y [z]   Report sum of X and Y with precision z\n");
    printf("   -sub X Y [z]   Report differenct of X and Y with precision z\n");
    printf("   -mul X Y [z]   Report product of X and Y with precision z\n");
    printf("   -div X Y [z]   Report quotient of X and Y with precision z\n");
    printf("   -mod X Y       Report X modulo Y\n");
    printf("   -sqr X [z]     Report square root of X with precision z\n");
    printf("   -log X [z]     Report log(x) (base 10) with precision z\n");
    printf("   -ln  X [z]     Report ln(x) (base e) with precision z\n");
    printf("   -lg2 X [z]     Report base 2 log(x) with precision z\n");
    printf("   -alog X [z]    Report anti-log (10^x) with precision z\n");
    printf("   -exp X [z]     Report exp (e^x) with precision z\n");
    printf("   -pw2 X [z]     Report power of 2 (2^x) with precision z\n");
    printf("   -sin X [z]     Report sin x (degrees) with precision z\n");
    printf("   -cos X [z]     Report cos x (degrees) with precision z\n");
    printf("   -tan X [z]     Report tan x (degrees) with precision z\n");
    printf("   -min X Y [z]   Report minimum of X and Y with precision z\n");
    printf("   -max X Y [z]   Report maximum of X and Y with precision z\n");
    printf("   -fmt X         Format string for value (i.e. range) X\n");
    printf("Note: to report e-based precision, put an \"e\" before [z]\n");
}
/**************************************************************************/
int NumI(int argc, char **argv)
{
    char bufS[100], chS[100], *cPC;
    int one,op,p,eform;
    double r1R, r2R, anR;
    
    if(argc < 2)
    {
        NumUse();   return(FALSE);
    }
    op = ParseNumOpI(argv[1],argc-2);
    if(op == 0)
    {
        NumUse();   return(FALSE);
    }
    eform = FALSE;
    /***
    *   Get arguments
    */
    anR = BAD_R;
    switch(op)
    {
        /***
        *   2 mand; 3 optional args
        */
        case OP_ADD:
        case OP_SUB:
        case OP_MUL:
        case OP_DIV:
        case OP_MOD:
        case OP_MIN:
        case OP_MAX:
            r1R = r2R = BAD_R;
            if(argc > 4)
            {
                p = BAD_I;
                sprintf(bufS,"%s %s",argv[2],argv[3]);
                sscanf(bufS,"%lf %lf",&r1R,&r2R);
                /***    
                *   Precision; print E notation?
                */
                sprintf(chS,"%s",argv[4]);
                cPC = chS;
                if(toupper(INT(*cPC))=='E')
                {
                    eform = TRUE;   
                    while((*cPC)&&(!isdigit(INT(*cPC))))
                    {   cPC++;  }
                    sscanf(cPC,"%d",&p);
                }
                else
                {
                    sscanf(argv[4],"%d",&p);
                }
                if(BAD_REAL(r1R) || BAD_REAL(r2R) || BAD_INT(p))
                {
                    printf("Bad numerical arguments: %s %s %s\n",
                        argv[2],argv[3],argv[4]);
                    NumUse();   return(FALSE);
                }
            }
            else
            {
                p = 0;
                sprintf(bufS,"%s %s",argv[2],argv[3]);
                sscanf(bufS,"%lf %lf",&r1R,&r2R);
                if(BAD_REAL(r1R) || BAD_REAL(r2R))
                {
                    printf("Bad numerical arguments: %s %s\n",
                        argv[2],argv[3]);
                    NumUse();   return(FALSE);
                }
            }
            if((op == OP_DIV)&&(ABS_VAL(r2R)<MIN_SIZE))
            {
                printf("Divisor too small: %s\n",argv[4]);
                NumUse();   return(FALSE);
            }
            if(op == OP_MOD)
            {
                p = 0;
            }
            break;
        /***
        *   1 mand; 2 optional args
        */
        case OP_LOG:
        case OP_LN:
        case OP_LG2:
        case OP_PW2:
        case OP_ALOG:
        case OP_EXP:
        case OP_SIN:
        case OP_COS:
        case OP_TAN:
        case OP_SQRT:
        case OP_FMT:
            r1R = BAD_R;
            if(argc > 3)
            {
                p = BAD_I;
                sprintf(bufS,"%s",argv[2]);
                sscanf(bufS,"%lf",&r1R);
                /***    
                *   Precision; print E notation?
                */
                sprintf(chS,"%s",argv[3]);
                cPC = chS;
                if(toupper(INT(*cPC)) == 'E')
                {
                    eform = TRUE;   
                    while((*cPC)&&(!isdigit(INT(*cPC))))
                    {   cPC++;  }
                    sscanf(cPC,"%d",&p);
                }
                else
                {
                    sscanf(argv[3],"%d",&p);
                }
                if(BAD_REAL(r1R) || BAD_INT(p))
                {
                    printf("Bad numerical arguments: %s %s\n",
                        argv[2],argv[3]);
                    NumUse();   return(FALSE);
                }
            }
            else
            {
                p = 0;
                sprintf(bufS,"%s",argv[2]);
                sscanf(bufS,"%lf",&r1R);
                if(BAD_REAL(r1R))
                {
                    printf("Bad numerical arguments: %s\n", argv[2]);
                    NumUse();   return(FALSE);
                }
            }
            if( ((op==OP_LN)||(op==OP_LOG)||(op==OP_LG2)) && (r1R<MIN_LOG) )
            {
                printf("Argument too small / negative: %s\n",argv[2]);
                NumUse();   return(FALSE);
            }
            break;
    }
    /****
    *   Do the work
    */
    one = FALSE;
    switch(op)
    {
        case OP_ADD:    anR = r1R + r2R;    one++;  break;
        case OP_SUB:    anR = r1R - r2R;    one++;  break;
        case OP_MUL:    anR = r1R * r2R;    one++;  break;
        case OP_DIV:    anR = r1R / r2R;    one++;  break;
        case OP_MOD:    anR = INT(r1R) % INT(r2R);  one++;  break;
        case OP_LOG:    anR = log(r1R) / LOG_E_TEN_R;   one++;  break;
        case OP_LN:     anR = log(r1R);     one++;  break;
        case OP_LG2:    anR = log(r1R) / LOG_E_TWO_R;   one++;  break;
        case OP_ALOG:   anR = exp(r1R * LOG_E_TEN_R);       one++;  break;
        case OP_EXP:    anR = exp(r1R);     one++;  break;
        case OP_PW2:    anR = exp(r1R * LOG_E_TWO_R);       one++;  break;
        case OP_SIN:    anR = SIN_R(r1R);   one++;  break;
        case OP_COS:    anR = COS_R(r1R);   one++;  break;
        case OP_TAN:    anR = TAN_R(r1R);   one++;  break;
        case OP_SQRT:   anR = SQRT_R(r1R);  one++;  break;
        case OP_MIN:    anR = MIN_NUM(r1R,r2R); one++;  break;
        case OP_MAX:    anR = MAX_NUM(r1R,r2R); one++;  break;
        case OP_FMT:    
            HandleFormatString(r1R,NULL);
            break;
    }
    if(one)
    {
        if(p > 0)
        {
            if(eform)
                sprintf(bufS,"%%1.%de\n",p);
            else
                sprintf(bufS,"%%1.%df\n",p);
            printf(bufS,anR);   
        }
        else
        {
            printf("%d\n",ROUND(anR));
        }
    }
    return(TRUE);
}
/************************************************************************
*   What print precision to use?
*/
void HandleFormatString(REAL rR,FILE *outPF)
{
    int p;

    HAND_NFILE(outPF);
    p = 0;
    while(rR<2.0)
    {
        p++;
        rR *= 10.0;
    }
    fprintf(outPF,"%%1.%df\n",p);
}
/************************************************************************
*   Parses numerical operator; i.e. command line option starting with "-"
*/
int ParseNumOpI(char *opS, int nargs)
{
    int op,minop;

    op = minop = 0;
    if(opS[0] != '-')
        return(FALSE);
    switch(toupper(opS[1]))
    {
        case 'A':   
            switch(toupper(opS[2]))
            {
                case 'D': op = OP_ADD;  minop = 2;  break;
                case 'L': op = OP_ALOG; minop = 1;  break;
            }
            break;
        case 'C':   
            switch(toupper(opS[2]))
            {
                case 'O':   op = OP_COS;    minop = 1;  break;
            }
            break;
        case 'D':   op = OP_DIV;    minop = 2;  break;
        case 'E':   op = OP_EXP;    minop = 1;  break;
        case 'F':   op = OP_FMT;    minop = 1;  break;
        case 'L':   
            switch(toupper(opS[2]))
            {
                case 'O': op = OP_LOG;  minop = 1;  break;
                case 'N': op = OP_LN;   minop = 1;  break;
                case 'G': op = OP_LG2;  minop = 1;  break;
            }
            break;
        case 'M':   
            switch(toupper(opS[2]))
            {
                case 'O': op = OP_MOD;  minop = 2;  break;
                case 'U': op = OP_MUL;  minop = 2;  break;
                case 'I': op = OP_MIN;  minop = 2;  break;
                case 'A': op = OP_MAX;  minop = 2;  break;
            }
            break;
        case 'P':   op = OP_PW2;    minop = 1;  break;
        case 'S':   
            switch(toupper(opS[2]))
            {
                case 'I': op = OP_SIN;  minop = 1;  break;
                case 'Q': op = OP_SQRT; minop = 1;  break;
                case 'U': op = OP_SUB;  minop = 2;  break;
            }
            break;
        case 'T':   op = OP_TAN;    minop = 1;  break;
    }
    if(op == 0)
    {
        printf("Unrecognized option; %s\n",opS);
        return(FALSE);
    }
    if(nargs < minop)
    {
        printf("Too few arguments; %d required\n",minop);
        return(FALSE);
    }
    return(op);
}
