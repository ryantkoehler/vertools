/*
* putil.c
*
* Copyright 2017 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <stdarg.h>
#include "prim.h"

#define DB_PARS     if(DB[5]) 

/** input parser **/
#define BOOL_TYPE   0
#define INT_TYPE    1
#define REAL_TYPE   2           /* Now treated as double; 11/10/04 RTK */
#define STR_TYPE    3
#define DOUB_TYPE   4
#define CHAR_TYPE   5

#define MAX_ARGS        200     /* Max supported command line args */
#define MAX_ARG_OPT     4       /* Max values associated with one arg key */

typedef struct V_ARG
{
     char key[DEF_BS];          /* Keyword to trigger this arg */
     int index;                 /* Index in collection (for debugging) */
     int type;                  /* Flag for type of values for this arg */
     int opt;                   /* Flag if this is an optional arg */
     int nvic;                  /* Number of associated values */
     int *vics[MAX_ARG_OPT];    /* Associated int values */
     REAL *rvics[MAX_ARG_OPT];  /* Associated REAL values */
     double *dvics[MAX_ARG_OPT];/* Associated double values */
     char *cvics[MAX_ARG_OPT];  /* Associated char values */
}V_ARG;

#define NEG_NUM(bu) ((bu[0]=='-')&&(isdigit(INT(bu[1]))))

PRIV_I GivenListedKeysMatchI(char *inkeyS, char *argkeyS);
PRIV_V ReportAmbigArgs(char *inkeyS, int n_mat, int *mats, V_ARG *varglisPO);
PRIV_V CleanV_argList(V_ARG *varglisPO, int num);
PRIV_I ParseSingVargArgI(V_ARG *vargPO,char *bufS,int j);

/*************************************************************************
*   Sets passed (pointers to) variables based on type/number specified in
*   a description string (arg 3).  If ok arguments are found returns TRUE
*
*   Optional arguments follow a specifer with a leading dash; e.g. -option
*
*   If the argument -noDB, or nodebug is found, debug is cleard
*   If the argument -dumpDB is given, debug flags are reported
*/
int ParseArgsI(int argc, char **argv, char *formPC, ...)
{
    int i, j, comarg, n_pvarg, n_argblock, n_mand, cur_mand, *vPI;
    int n_mat, matches[MAX_ARGS];
    char bufS[255], *cPC, *vPC;
    REAL *vPR;
    double *vPD;
    V_ARG v_argsOA[MAX_ARGS], *vargPO, *tvargPO;
    va_list ap;

    DB_PARS
    {
        DB_PrI(">> ParseArgsI, argc %d\n",argc);
        for(i=0;i<argc;i++)
        {
            DB_PrI("+ arg %2d |%s|\n",i,argv[i]);
        }
        DB_PrI("+ form: |%s|\n",formPC);
    }
    /***
    *   Initialize internal varg specifier and counts:
    *   n_argblock = number of argument blocks,
    *   n_pvarg = number of actual (variable length) arguments passed here
    */
    CleanV_argList(v_argsOA, MAX_ARGS);
    DB_PARS DB_PrI("+ varg structs initilized\n");
    n_argblock = n_pvarg = 0;
    va_start(ap,formPC);
    /***
    *   scan down the argument-specifying format string
    */
    cPC = formPC;
    while(*cPC != '\0')
    {
        if(n_argblock >= MAX_ARGS) {
            printf("Parsing problem with variable argument format\n");
            printf("|%s|\n",formPC);
            printf("TOO MANY ARGUMENTS, %d max\n",MAX_ARGS);
            ERR("ParseArgsI","bad parseable argument format");
            return(FALSE);
        }
        vargPO = &v_argsOA[n_argblock];
        /***
        *   Isolate current token
        */
        sscanf(cPC,"%s",bufS);
        DB_PARS DB_PrI("+  |%s| ",bufS);
        /***
        *   Optional argument starts with "-"
        */
        if(*bufS == '-') {
            DB_PARS DB_PrI("optional ");
            vargPO->opt = TRUE;
            sprintf(vargPO->key,"%s",bufS);
            DB_PARS DB_PrI("= key |%s|, ",vargPO->key);
            NEXT_WORD(cPC);
            if(*cPC == '\0')
            {
                printf("Parsing problem with variable argument format\n");
                printf("|%s|\n",formPC);
                printf("ARGUMENT %d has no type specifier\n",n_argblock);
                ERR("ParseArgsI","bad parseable argument format");
                return(FALSE);
            }
        }
        else {
            DB_PARS DB_PrI("manditory ");
        }
        /***
        *   Set the type of variable associated 
        */
        switch(*cPC)
        {
            case 'B':   vargPO->type = BOOL_TYPE;   break;
            case 'I':   vargPO->type = INT_TYPE;    break;
            case 'R':   vargPO->type = REAL_TYPE;   break;
            case 'D':   vargPO->type = DOUB_TYPE;   break;
            case 'C':   vargPO->type = CHAR_TYPE;   break;
            case 'S':   vargPO->type = STR_TYPE;    break;
            default:
                printf("Parsing problem with variable argument format\n");
                printf("|%s|\n",formPC);
                printf("ARGUMENT %d has bad type specifier |%c|\n",n_argblock,*cPC);
                ERR("ParseArgsI","bad parse-specifying argument format");
                return(FALSE);
        }
        /***
        *   If not boolean, check for number of associated arguments
        */
        if(vargPO->type != BOOL_TYPE) {
            cPC++;
            if(isdigit(INT(*cPC))) {
                sscanf(cPC,"%d",&vargPO->nvic);
            }
        }
        if(vargPO->nvic > MAX_ARG_OPT) {
            printf("Parsing problem with variable argument format\n");
            printf("|%s|\n",formPC);
            printf("ARGUMENT %d has lists %d values but %d max |%c|\n",
                n_argblock, vargPO->nvic, MAX_ARG_OPT, *cPC);
            ERR("ParseArgsI","bad parse-specifying argument var count");
            return(FALSE);
        }
        DB_PARS DB_PrI("type %d, nvic %d\n",vargPO->type,vargPO->nvic);
        /***
        *   Attach the actuall passed variable argument(s) to the key vargPO
        */
        for(i=0; i<vargPO->nvic; i++)
        {
            vPI = NULL; 
            vPR = NULL; 
            vPD = NULL; 
            vPC = NULL;
            switch(vargPO->type)
            {
                case REAL_TYPE: vPR = va_arg(ap,REAL*);     break;
                case DOUB_TYPE: vPD = va_arg(ap,double*);   break;
                case CHAR_TYPE: vPC = va_arg(ap,char*);     break;
                default:        vPI = va_arg(ap,int*);
            }
            DB_PARS DB_PrI("+    vi_arg...[%d], %d\n",i,n_pvarg);
            if( (vPI == NULL) && (vPR == NULL) && (vPD == NULL) && (vPC==NULL)) {
                printf("Bad variable argument for format\n");
                printf(" vic[%d], type=%d\n",i,vargPO->type);
                printf("|%s|\n",formPC);
                printf("ARGUMENT %d not passed\n",n_pvarg);
                ERR("ParseArgsI","bad variable argument passed");
                return(FALSE);
            }
            vargPO->dvics[i] = vPD;
            vargPO->rvics[i] = vPR;
            vargPO->cvics[i] = vPC;
            vargPO->vics[i] = vPI;
            n_pvarg++;
        }
        /***
        *   Next key in format string
        */
        NEXT_WORD(cPC);
        n_argblock++;   
    }
    /***
    *   Nothing read in?
    */
    if(n_argblock == 0) {
        printf("Parsing problem with variable argument format\n");
        printf("|%s|\n",formPC);
        printf("NO ARGUMENTS READ\n");
        ERR("ParseArgsI","bad parseable argument format");
    }
    /***
    *   Command line is now loaded into parts
    *   Count manditory argument/s; return if already too few
    */
    DB_PARS DB_PrI("+ counting manditory arguments\n");
    n_mand = j = 0;
    for(i = 0; i < n_argblock; i++)
    {
        vargPO = &v_argsOA[i];
        DB_PARS DB_PrI("+ [%d] |%s| opt=%d, nvic=%d\n",
            i,vargPO->key,vargPO->opt,vargPO->nvic);
        if(vargPO->opt) {
            continue;
        }
        n_mand++;
        j += vargPO->nvic;
    }
    DB_PARS DB_PrI("+ n_mand=%d  j=%d, argc = %d\n",n_mand,j,argc);
    if(j > (argc -1)) {
        printf("\n");
        printf("Missing manditory command line arguments;\n");
        printf("    %d non-optional argument(s) expected\n",n_mand);
        printf("\n");
        DB_PARS DB_PrI("<< ParseArgsI FALSE\n");
        return(FALSE);
    }
    /***
    *   Party with actual argument list: passed command vs option keys
    */
    DB_PARS DB_PrI("+ now parsing passed command line:\n");
    cur_mand = 0;   
    comarg = 1;
    while(comarg < argc)
    {
        sprintf(bufS,"%s", argv[comarg]);
        DB_PARS DB_PrI("+ [%2d] |%s|\n",comarg,bufS);
        comarg++;
        /***
        *   Special debug things (parsed elsewhere at init)
        */
        if( EQSTRING(bufS,"-noDB",5) || EQSTRING(bufS,"-nodb",5) ) {
            continue;
        }
        if( EQSTRING(bufS,"-dumpDB",7) ) {
            continue;
        }
        /***
        *   Search collection of arguments for matching keyword(s) if starts with '-'
        *   Don't consider '-' only (i.e. stdin) or neg number cases
        */
        n_mat = 0;
        if( (bufS[0] == '-') && (strlen(bufS) > 1) && (!NEG_NUM(bufS)) ) {
            for(j = 0; j < n_argblock; j++)
            {
                tvargPO = &v_argsOA[j];
                DB_PARS DB_PrI("+  checking [%2d] key|%s|\n",tvargPO->index,tvargPO->key);
                /***
                *   If not optional, has empty key, or is only '-', ignore 
                */
                if( (!tvargPO->opt) || NO_S(tvargPO->key) ) {
                    continue;
                }
                if(GivenListedKeysMatchI(bufS, tvargPO->key)) {
                    matches[n_mat] = j;
                    n_mat++;
                }
            }
        }
        /***
        *   Ambiguous match?
        */
        if(n_mat > 1) {
            ReportAmbigArgs(bufS, n_mat, matches, v_argsOA);
            return(FALSE);
        }
        /***
        *   if one argument match then this is an optional argument
        */
        if(n_mat == 1) {
            vargPO = &v_argsOA[matches[0]];
            DB_PARS DB_PrI("+  Found varg[%d]  key |%s| ",vargPO->index,vargPO->key);
            if(vargPO->type == BOOL_TYPE) {
                *vargPO->vics[0] = !(*vargPO->vics[0]);
                DB_PARS DB_PrI("= bool toggle (%p = %d)\n",
                    vargPO->vics[0],*vargPO->vics[0]);
            }
            else {
                DB_PARS DB_PrI("= eat %d args\n",vargPO->nvic);
                for(j= 0; j<vargPO->nvic; j++)
                {
                    if(comarg >= argc) {
                        printf("\n");
                        printf(
                            "Option \"%s\" is missing arguments; %d required\n",
                            vargPO->key, vargPO->nvic);
                        printf("\n");
                        return(FALSE);
                    }
                    sprintf(bufS,"%s",argv[comarg++]);
                    DB_PARS DB_PrI("+    %4d |%s|\n",j,bufS);
                    if(!ParseSingVargArgI(vargPO,bufS,j)) {
                        printf("    Problem with \"%s\"\n", bufS);
                        printf("\n");
                        return(FALSE);
                    }
                }
            }
        }
        /***
        *   non-optional argument; find next non_optional one and set it
        */
        else
        {
            if(cur_mand >= n_mand) {
                printf("\n");
                printf("Invalid or unrecognized option: \"%s\"\n",bufS);
                printf("\n");
                return(FALSE);
            }
            DB_PARS DB_PrI("Null key; looking for next mand (cur %d)\n",cur_mand);
            i = 0;
            vargPO = NULL;
            for(j = 0; j < n_argblock; j++)
            {
                tvargPO = &v_argsOA[j];
                if(tvargPO->opt) {
                    continue;
                }
                if(i==cur_mand) {
                    vargPO = tvargPO;
                    cur_mand++;
                    break;
                }
                i++;
            }
            if(vargPO == NULL) {
                DB_PARS DB_PrI("<< ParseArgsI FALSE, missing mand command\n");
                return(FALSE);
            }
            DB_PARS DB_PrI("+  mand %d, type %d = eat %d args\n",
                cur_mand,vargPO->type,vargPO->nvic);
            for(j= 0; j<vargPO->nvic; j++)
            {
                if(j > 0) {
                    if(comarg >= argc) {
                        printf("\n");
                        printf("Manditory input parameters are missing\n");
                        printf("\n");
                        return(FALSE);
                    }
                    sprintf(bufS,"%s",argv[comarg++]);
                }
                DB_PARS DB_PrI("+    %4d |%s|\n",j,bufS);
                if(!ParseSingVargArgI(vargPO,bufS,j)) {
                    printf("    Problem with \"%s\"\n", bufS);
                    printf("\n");
                    return(FALSE);
                }
            }
        }
    }
    /***
    *
    */
    if(cur_mand < n_mand) {
        printf("\n");
        printf("Missing manditory command line arguments;\n");
        printf("    %d non-optional argument(s) expected\n",n_mand);
        printf("\n");
        return(FALSE);
    }
    /***
    *   Clean and bail
    */
    va_end(ap);
    DB_PARS DB_PrI("<< ParseArgsI TRUE\n");
    return(TRUE);
}
/*************************************************************************/
PRIV_I GivenListedKeysMatchI(char *inkeyS, char *argkeyS)
{
    int len;

    len = MIN_NUM(strlen(inkeyS), strlen(argkeyS));
    return EQSTRING(inkeyS,argkeyS,len);
}
/*************************************************************************/
PRIV_V ReportAmbigArgs(char *inkeyS, int n_mat, int *mats, V_ARG *varglisPO)
{
    int i;
    V_ARG *vargPO;

    printf("\n");
    printf("Ambiguous keyword argument\n");
    printf("    %s matches %d options:\n", inkeyS, n_mat);
    for(i=0; i<n_mat; i++)
    {
        vargPO = &varglisPO[mats[i]];
        printf("    %s ?\n", vargPO->key);
    }
    printf("\n");
}
/*************************************************************************/
PRIV_V CleanV_argList(V_ARG *varglisPO, int num)
{
    int i;
    V_ARG *vargPO;

    for(i=0; i<num; i++)
    {
        vargPO = &varglisPO[i];
        vargPO->key[0] = '\0';
        vargPO->type = STR_TYPE;
        vargPO->opt = FALSE;
        vargPO->nvic = 1;
    }
}
/*************************************************************************/
PRIV_I ParseSingVargArgI(V_ARG *vargPO, char *bufS, int j)
{
    int ival;
    REAL rvalR;
    double dvalD;
    char cvalC;

    DB_PARS DB_PrI(">> ParseSingVargArgI |%s| setting in [%d][%d]\n",
        bufS,vargPO->index,j);
    if(j >= vargPO->nvic) {
        printf("\n nvic = %d, j = %d\n",vargPO->nvic,j);
        ERR("ParseSingVargArgI","bad argument index");
    }
    /***
    *   If starts with "-" can only be exactly that (for stdin) or neg number
    */
    if (bufS[0] == '-') {
        if ( (strcmp(bufS,"-")) && (!isdigit(INT(bufS[1]))) ) {
            printf("\n");
            printf("Invalid option syntax; options are as \"[-<key>] <val>\"\n");
            return(FALSE);
        }
    }
    switch(vargPO->type)
    {
        case INT_TYPE: 
            if( (!isdigit(INT(bufS[0]))) && (!NEG_NUM(bufS)) )
            {
                printf("\n missing int argument (%s # %d)\n",vargPO->key,j+1);
                return(FALSE);
            }
            sscanf(bufS,"%d",&ival);
            *vargPO->vics[j] = ival;
            break;
        case REAL_TYPE:
            if( (!isdigit(INT(bufS[0]))) && (!NEG_NUM(bufS)) )
            {
                printf("\n missing real argument (%s # %d)\n",vargPO->key,j+1);
                return(FALSE);
            }
            sscanf(bufS,"%lf",&rvalR);
            *vargPO->rvics[j] = (REAL)rvalR;
            break;
        case DOUB_TYPE:
            if( (!isdigit(INT(bufS[0]))) && (!NEG_NUM(bufS)) )
            {
                printf("\n missing double arg (%s # %d)\n",vargPO->key,j+1);
                return(FALSE);
            }
            sscanf(bufS,"%lf",&dvalD);
            *vargPO->dvics[j] = (double)dvalD;
            break;
        case CHAR_TYPE:
            sscanf(bufS,"%c",&cvalC);
            *vargPO->cvics[j] = (char)cvalC;
            break;
        case STR_TYPE: 
            strcpy( (char *)vargPO->vics[j], bufS);
/***    Only gets first non-space token from supplied string
            sscanf(bufS,"%s",(char *)vargPO->vics[j]);
*/
            break;
        default:
            printf("BAD V_ARG type %d\n",vargPO->type);
            ERR("ParseSingVargArgI","bad varg type");
    }
    return(TRUE);
}
