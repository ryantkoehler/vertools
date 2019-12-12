/*
* score_io.c
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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "prim.h"
#include "score.h"

#define DB_SCF if(DB[61])

/*************************************************************************
*   Load settings for SCFIELD from parameter file
*/
int LoadScfieldI(FILE *inPF, int error, SCFIELD **sfPPO)
{
    char bufS[BBUFF_SIZE],nameS[NSIZE];
    int n, s;
    DOUB xD,yD;
    SCFIELD *sfPO;

    DB_SCF DB_PrI(">> LoadScfieldI\n");
    sfPO = NULL;
    INIT_S(nameS); 
    n = s = 0;
    while(fgets(bufS,BLINEGRAB,inPF))
    {
        if(COM_LINE(bufS)) {   
            continue;
        }
        if(BlankStringI(bufS)) {
            continue;
        }
        /***
        *   Start indicator
        */
        if(strstr(bufS,SCFIELD_START_S)) {
            if(s) {
                if(error) {
                    ReportParseErrorLine(bufS,"LoadScfieldI","Second start token found");
                }
                DB_SCF DB_PrI("<< LoadScfieldI second start token FALSE\n");
                return(FALSE);
            }
            sscanf(bufS,"%*s %s",nameS);
            s = 1;
            continue;
        }
        /***
        *   Some non-blank, non-comment, non-start line but not yet start?
        */
        if(s < 1) {
            if(error) {
                ReportParseErrorLine(bufS,"LoadScfieldI","No start token found");
            }
            DB_SCF DB_PrI("<< LoadScfieldI no start token FALSE\n");
            return(FALSE);
        }
        /***
        *   End indicator
        */
        if(strstr(bufS,SCFIELD_END_S)) {
            s++;
            break;
        }
        /***
        *   If no SCFIELD yet, create on
        */
        if( !sfPO ) {
            sfPO=CreateScfieldPO(nameS);
        }
        if( !sfPO ) {
            DB_SCF DB_PrI("<< LoadScfieldI create failed FALSE\n");
            return(FALSE);
        }
        /***
        *   X Y value pair
        */
        xD = yD = BAD_D;
        sscanf(bufS,"%lf %lf",&xD,&yD);
        if( BAD_DOUB(xD) || BAD_DOUB(yD) ) {
            if(error) {
                ReportParseErrorLine(bufS,"LoadScfieldI","Failed to get X Y number pair");
            }
            DB_SCF DB_PrI("<< LoadScfieldI bad XY FALSE\n");
            return(FALSE);
        } 
        AppendNumlistDoubI(sfPO->xvals, xD);
        AppendNumlistDoubI(sfPO->yvals, yD);
        n++;
    }
    /***
    *   Check we got end token and enough values
    */
    if( (s != 2) || (n < 1) ) {
        if(error) {
            sprintf(bufS,"Failed to get end token (%d) or enough X Y pairs (%d) for %s",s,n,nameS);
            ReportParseErrorLine(NULL,"LoadScfieldI",bufS);
        }
        CHECK_SCFIELD(sfPO);
        DB_SCF DB_PrI("<< LoadScfieldI block and val count not good FALSE\n");
        return(FALSE);
    }
    /***
    *   Check all X values should be in increasing order 
    *   Set min / max values
    */
    if( !NumlistIsSortedI(sfPO->xvals, 1, MIN_SCF_XSTEP) ) {
        if(error) {
            sprintf(bufS,"X values must increase by at least (%f) for %s",MIN_SCF_XSTEP,nameS);
            DumpNumlist(sfPO->xvals, -1, -1, NULL, NULL);
            DumpNumlist(sfPO->yvals, -1, -1, NULL, NULL);
            ReportParseErrorLine(NULL,"LoadScfieldI",bufS);
        }
        CHECK_SCFIELD(sfPO);
        DB_SCF DB_PrI("<< LoadScfieldI block and val count not good FALSE\n");
        return(FALSE);
    }
    GetNumlistDoubI(sfPO->xvals, 0, &sfPO->minx);
    GetNumlistDoubI(sfPO->xvals, n-1, &sfPO->maxx);
    sfPO->n = n;
    /***
    *   Set pointer and return
    */
    *sfPPO = sfPO;
    DB_SCF DB_PrI("<< LoadScfieldI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Attempt to load an array of SCFIELDS;
*   These are returned via an array of pointers 
*/
int LoadScfieldArrayI(FILE *inPF, SCFIELD ***scPPPO)
{
    int n,i,fpos;
    char bufS[BBUFF_SIZE];
    SCFIELD *sfPO, *fsetPA[MAX_SCF_LIST], **scPPO;

    DB_SCF DB_PrI(">> LoadScfieldArrayI\n");
    *scPPPO = NULL;
    /***
    *   Go through file, looking for start blocks and trying to parse those
    */
    n = 0;
    fpos = ftell(inPF);
    while(fgets(bufS,BLINEGRAB,inPF)) {
        if(COM_LINE(bufS)) {   
            continue;
        }
        if(BlankStringI(bufS)) {
            continue;
        }
        /***
        *   Start indicator, try and parse and save pointer
        */
        if(strstr(bufS,SCFIELD_START_S)) {
            if(n >= MAX_SCF_LIST) {
                PROBLINE;
                printf("Too many scores indicated, max = %d; Previous:\n",MAX_SCF_LIST);
                ReportScfieldI(fsetPA[n-1], FALSE, NULL);
                return(FALSE);
            }
            fseek(inPF,fpos,0);
            if(! LoadScfieldI(inPF, TRUE, &sfPO) ) {
                PROBLINE;
                printf("Failed to load Scfield [%d]\n",n);
                return(FALSE);
            }
            fsetPA[n++] = sfPO;
        }
        fpos = ftell(inPF);
    }
    if(n<1) {
        DB_SCF DB_PrI("<< LoadScfieldArrayI (none) %d\n",n);
        return(n);
    }
    /***
    *   Allocate SCFLIELD pointer array and fill
    */
    scPPO = (SCFIELD **)ALLOC(n,sizeof(SCFIELD *));
    if(!scPPO) {
        printf("Failed to allocate for %d score fields\n",n);
        return(FALSE);
    }
    for(i=0;i<n;i++) {
        DB_SCF {
            DB_PrI("+ SCFIELD[%d] :\n",i);
            ReportScfieldI(fsetPA[i], "+ ", NULL);
        }
        scPPO[i] = fsetPA[i];
    }
    *scPPPO = scPPO;
    DB_SCF DB_PrI("<< LoadScfieldArrayI %d\n",n);
    return(n);
}
/****************************************************************************
*   Dump a scfield to file
*   The flag full dictates output behavior
*/
int ReportScfieldI(SCFIELD *sfPO, char *prefixS, FILE *oPF)
{
    int i;
    DOUB xD, yD;
    char preS[DEF_BS];

    VALIDATE(sfPO,SCFIELD_ID);
    HAND_NFILE(oPF);
    INIT_S(preS);
    if (prefixS != NULL) {
        sprintf(preS,"%s", prefixS);
    }
    /***
    *   Start keyword, X,Y values, End keyword
    */
    fprintf(oPF,"%s%s\t%s\n", preS, SCFIELD_START_S, sfPO->name);
    for(i=0; i< sfPO->n; i++) {
        GetNumlistDoubI(sfPO->xvals, i, &xD);
        GetNumlistDoubI(sfPO->yvals, i, &yD);
        fprintf(oPF,"%s%0.3f\t%0.3f\n", preS, xD, yD);
    }
    fprintf(oPF,"%s%s\n", preS, SCFIELD_END_S);
    return(sfPO->n);
}
/****************************************************************************
*   Report an array of SCFIELDs
*/
int ReportScfieldArrayI(SCFIELD **scPPO, int n, char *prefixS, FILE *oPF)
{
    int i;

    HAND_NFILE(oPF);
    for(i=0; i<n; i++) {
        if ( (prefixS != NULL) && (i>0) ) {
            fprintf(oPF,"%s\n", prefixS);
        }
        ReportScfieldI(scPPO[i], prefixS, oPF);
    }
    return(n);
}
