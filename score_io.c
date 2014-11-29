/*
* score_io.c
*
* Copyright 2014 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <string.h>
#include <math.h>
#include "prim.h"
#include "score.h"

#define DB_SCF if(DB[61])

/*************************************************************************
*	Load settings for SCFIELD from parameter file
*/
int LoadScfieldI(FILE *inPF, int error, SCFIELD **sfPPO)
{
	char bufS[BBUFF_SIZE],nameS[NSIZE];
    DOUB xD,yD,xPI[MAX_SC_VALS], yPI[MAX_SC_VALS], stepD, fxD, dyD;
    int ngive, i, n, s, g;
    SCFIELD *sfPO;

    DB_SCF DB_PrI(">> LoadScfieldI\n");
    InitArrayI(xPI, IS_DOUB, 0, MAX_SC_VALS, 0.0);
    InitArrayI(yPI, IS_DOUB, 0, MAX_SC_VALS, 0.0);
	/***
    *   Should start with 'SCFIELD START' 
	*	First real line should be 
	*/
    INIT_S(nameS); 
    s = ngive = 0;
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
        *   X Y value pair
        */
        if( ngive >= MAX_SC_VALS ) {
            if(error) {
                sprintf(nameS,"Too many X Y pairs (%d = max)",MAX_SC_VALS);
                ReportParseErrorLine(bufS,"LoadScfieldI",nameS);
            }
            DB_SCF DB_PrI("<< LoadScfieldI too many XY FALSE\n");
            return(FALSE);
        }
        xD = yD = BAD_D;
        sscanf(bufS,"%lf %lf",&xD,&yD);
        if( BAD_DOUB(xD) || BAD_DOUB(yD) ) {
            if(error) {
                ReportParseErrorLine(bufS,"LoadScfieldI","Failed to get X Y number pair");
            }
            DB_SCF DB_PrI("<< LoadScfieldI bad XY FALSE\n");
            return(FALSE);
        } 
        if( (ngive>0) && (xD <= xPI[ngive-1]) ) {
            if(error) {
                sprintf(nameS,"X values must increase in X Y pairs (%f -vs- %f)",xD,xPI[ngive-1]);
                ReportParseErrorLine(bufS,"LoadScfieldI",nameS);
            }
            DB_SCF DB_PrI("<< LoadScfieldI X vals not increasing FALSE\n");
            return(FALSE);
        }
        DB_SCF DB_PrI("+  saving [%d] x=%f y=%f\n",ngive,xD,yD);
        xPI[ngive] = xD;
        yPI[ngive] = yD;
        ngive++;
    }
    DB_SCF DB_PrI("+ after block, s=%d ngive=%d\n",s,ngive);
    /***
    *   Check we got end token and enough values
    */
    if( (s != 2) || (ngive < 2) ) {
        if(error) {
            sprintf(bufS,"Failed to get end token (%d) or enough X Y pairs (%d)",s,ngive);
            ReportParseErrorLine(NULL,"LoadScfieldI",bufS);
        }
        DB_SCF DB_PrI("<< LoadScfieldI block and val count not good FALSE\n");
        return(FALSE);
    }
    /***
    *   Find min value separation (=step) and calculate number of values to hold
    *   All X values should be in increasing order 
    */
    stepD = xPI[ngive-1] - xPI[0];
    for(i=1; i<ngive; i++) {
        stepD = MIN_NUM(stepD, (xPI[i] - xPI[i-1]) );
    }
    n = INT ( (xPI[ngive-1] - xPI[0]) / stepD) + 1;
    DB_SCF DB_PrI("+ step=%f, n=%d\n",stepD,n);
    /***
    *   Create new SCFIELD to hold data
    */
    if( ! (sfPO=CreateScfieldPO(ngive, n, nameS)) ) {
        DB_SCF DB_PrI("<< LoadScfieldI create failed FALSE\n");
        return(FALSE);
    }
    DB_SCF DB_PrI("+ Created with ngive=%d n=%d\n",ngive,n);
    sfPO->min = xPI[0];
    sfPO->max = xPI[ngive-1];
    sfPO->step = stepD;
    /***
    *   Set real given values, after intializing all to impossible
    */
    InitArrayI(sfPO->yval, IS_DOUB, 0, n, BAD_D);
    for(i=0; i<ngive; i++) {
        sfPO->gxv[i] = xPI[i];
        sfPO->gyv[i] = yPI[i];
        s = INT( (xPI[i] - xPI[0]) / stepD);
        BOG_CHECK(s >= n);
        DB_SCF DB_PrI("+  S val[%d] setting given[%d] =%f\n",s,i,yPI[i]);
        sfPO->yval[s] = yPI[i];
    }
    /***
    *   Set non-given values to intrpolated values
    *   Both first and last values ( [0] and [n-1] should be real ) 
    */
    xD = xPI[0];
    g = 0;
    for(i=0; i< n; i++) {
        DB_SCF DB_PrI("+  I [%d] x=%0.3f g=%d\n",i,xD,g);
        if(BAD_DOUB(sfPO->yval[i]) ) {
            fxD = (xD - xPI[g-1] ) / (xPI[g] - xPI[g-1] );
            dyD = yPI[g] - yPI[g-1];
            sfPO->yval[i] = yPI[g-1] + fxD * dyD;;
            DB_SCF DB_PrI("+    fx=%0.3f dy=%0.3f y[g-1]=%0.3f val[%d]=%0.3f\n",
                fxD, dyD, yPI[g-1], i, sfPO->yval[i]);
        }
        else {
            DB_SCF DB_PrI("+    already set\n",i);
            g++;
        }
        xD += stepD;
    }
	/***
	*	Set pointer and return
	*/
    *sfPPO = sfPO;
    DB_SCF DB_PrI("<< LoadScfieldI TRUE\n");
	return(TRUE);
}
/****************************************************************************
*	Attempt to load an array of SCFIELDS;
*	These are returned via an array of pointers 
*/
int LoadScfieldArrayI(FILE *inPF, SCFIELD ***scPPPO)
{
	int n,i,fpos;
    char bufS[BBUFF_SIZE];
	SCFIELD *sfPO, *fsetPA[MAX_SC_VALS], **scPPO;

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
            if(n >= MAX_SC_VALS) {
                PROBLINE;
                printf("Too many scores indicated, max = %d; Previous:\n",MAX_SC_VALS);
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
	*	Allocate SCFLIELD pointer array and fill
	*/
	scPPO = (SCFIELD **)ALLOC(n,sizeof(SCFIELD *));
	if(!scPPO) {
		printf("Failed to allocate for %d score fields\n",n);
		return(FALSE);
	}
	for(i=0;i<n;i++) {
        DB_SCF {
            DB_PrI("+ SCFIELD[%d] :\n",i);
            ReportScfieldI(fsetPA[i], TRUE, NULL);
        }
		scPPO[i] = fsetPA[i];
    }
	*scPPPO = scPPO;
    DB_SCF DB_PrI("<< LoadScfieldArrayI %d\n",n);
	return(n);
}
/****************************************************************************
*	Dump a scfield to file
*   The flag full dictates output behavior
*/
int ReportScfieldI(SCFIELD *sfPO, int full, FILE *oPF)
{
	int i;
    DOUB xD;
    char preS[DEF_BS];

	VALIDATE(sfPO,SCFIELD_ID);
	HAND_NFILE(oPF);
    INIT_S(preS);
    if (full) {
        sprintf(preS,"#\t");
    }
    /***
    *   Standard fields; full > 0 prepend with #
    */
    fprintf(oPF,"%s%s\t%s\n", preS, SCFIELD_START_S, sfPO->name);
    for(i=0; i< sfPO->ngiven; i++) {
        fprintf(oPF,"%s%0.3f\t%0.3f\n", preS, sfPO->gxv[i], sfPO->gyv[i]);
    }
    fprintf(oPF,"%s%s\n", preS, SCFIELD_END_S);
    /***
    *   full > 1 = all the rest
    */
    if (full > 1) {
        fprintf(oPF,"#\tFull %d interpolated values\n",sfPO->n);
        fprintf(oPF,"#\tMin  %0.3f\n",sfPO->min);
        fprintf(oPF,"#\tMax  %0.3f\n",sfPO->max);
        fprintf(oPF,"#\tStep %0.3f\n",sfPO->step);
        xD = sfPO->min;
        for(i=0; i< sfPO->n; i++) {
            fprintf(oPF,"#\t[%d]\t%0.3f\t%0.3f\n", i, xD, sfPO->yval[i]);
            xD += sfPO->step;
        }
    }
	return(sfPO->n);
}
/****************************************************************************
*	Report an array of SCFIELDs
*/
int ReportScfieldArrayI(SCFIELD **scPPO, int n, int full, FILE *oPF)
{
	int i;

	HAND_NFILE(oPF);
    for(i=0; i<n; i++) {
        if ((full > 1) && (i>0)) {
            fprintf(oPF,"#\n");
        }
        ReportScfieldI(scPPO[i], full, oPF);
    }
    return(n);
}
