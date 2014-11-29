/*
* score.c
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

#define DB_SCM if(DB[60])
#define DB_SCF if(DB[61])

PRIV_I qsort_scorec_id(const void *e1, const void *e2);
PRIV_I qsort_scorec_id_decend(const void *e1, const void *e2);
PRIV_I qsort_scorecs(const void *e1, const void *e2);
PRIV_I qsort_scorecs_decend(const void *e1, const void *e2);

/**************************************************************************
*	Allocate n SCOREC structures
*/
SCOREC *CreateScorecsPO(int n)
{
	SCOREC *recsPO;

	recsPO = (SCOREC *)ALLOC(n,sizeof(SCOREC));
	if(!recsPO) {
		printf("Failed to allocate for %d scores\n",n);
		return(NULL);
	}
	InitScorecs(recsPO,n);
	return(recsPO);
}	
/**************************************************************************
*	Free SCOREC structure(s)
*/
int DestroyScorecsI(SCOREC *recsPO)
{
	BOG_CHECK(!recsPO);
	FREE(recsPO);
	return(TRUE);
}
/**************************************************************************
*	Initialize SCOREC(s) with indices and values 
*/
void InitScorecs(SCOREC *recsPO,int n)
{
	int i;

	for(i=0;i<n;i++) {
		recsPO[i].id = i;
		recsPO[i].sc = 0.0;
	}
}
/**************************************************************************
*	Create an SCFIELD structure
*/
SCFIELD *CreateScfieldPO(int ngive, int nval, char *nameS)
{
	SCFIELD *sfPO;

	DB_SCM DB_PrI(">> CreateScfieldPO\n");
	if( ! (sfPO = (SCFIELD *)ALLOC(1,sizeof(SCFIELD)) ) ) {
		return(NULL);
	}
	sfPO->ID = SCFIELD_ID;
    if(ngive > 0) {
        sfPO->gxv = (DOUB *)ALLOC(ngive, sizeof(DOUB));
        sfPO->gyv = (DOUB *)ALLOC(ngive, sizeof(DOUB));
        if( (!sfPO->gxv) || (!sfPO->gyv) ) {
            CHECK_SCFIELD(sfPO);
            return(NULL);
        }
        sfPO->ngiven = ngive;
    }
    if(nval > 0) {
	    if( ! (sfPO->yval = (DOUB *)ALLOC(nval,sizeof(DOUB)) ) ) {
            CHECK_SCFIELD(sfPO);
            return(NULL);
        }
        sfPO->n = nval;
    }
    if(nameS) {
        SetScfieldNameI(sfPO,nameS,NSIZE);
    }
	InitScfield(sfPO);
	DB_SCM DB_PrI("<< CreateScfieldPO %p\n",sfPO);
	return(sfPO);
}
/**************************************************************************
*	Frees SCFIELD datatype 
*/
int DestroyScfieldI(SCFIELD *sfPO)
{
	VALIDATE(sfPO,SCFIELD_ID);
	CHECK_FREE(sfPO->yval);
	CHECK_FREE(sfPO->gxv);
	CHECK_FREE(sfPO->gyv);
	FREE(sfPO);
	return(TRUE);
}
/****************************************************************************
*	Initialize SCFIELD structure
*/
void InitScfield(SCFIELD *sfPO)
{
	VALIDATE(sfPO,SCFIELD_ID);
    InitArrayI(sfPO->yval, IS_DOUB, 0, sfPO->n, 0.0);
}
/****************************************************************************
*	Sort SCOREC list[n] by value
*/
void SortScorecVals(SCOREC *scPO,int n,int dir)
{
	if(dir < 0) {
		qsort(scPO,n,sizeof(SCOREC),qsort_scorecs_decend);
	}
	else {
		qsort(scPO,n,sizeof(SCOREC),qsort_scorecs);
	}
}
/***************************************************************************
*	Sort SCOREC list[n] by index
*/
void SortScorecIds(SCOREC *scPO,int n,int dir)
{
	if(dir < 0) {
		qsort(scPO,n,sizeof(SCOREC),qsort_scorec_id_decend);
	}
	else {
		qsort(scPO,n,sizeof(SCOREC),qsort_scorec_id);
	}
}
/****************************************************************************/
PRIV_I qsort_scorec_id(const void *e1, const void *e2)
{
	if( INT( ((SCOREC *)e1)->id ) > INT( ((SCOREC *)e2)->id) ) {
		return(1);
	}
	if( INT( ((SCOREC *)e1)->id ) < INT( ((SCOREC *)e2)->id) ) {
		return(-1);
	}
	return(0);
}
/****************************************************************************/
PRIV_I qsort_scorec_id_decend(const void *e1, const void *e2)
{
	if( INT( ((SCOREC *)e1)->id ) > INT( ((SCOREC *)e2)->id) ) {
		return(-1);
	}
	if( INT( ((SCOREC *)e1)->id ) < INT( ((SCOREC *)e2)->id) ) {
		return(1);
	}
	return(0);
}
/****************************************************************************/
PRIV_I qsort_scorecs(const void *e1, const void *e2)
{
	if( RNUM( ((SCOREC *)e1)->sc ) > RNUM( ((SCOREC *)e2)->sc) ) {
		return(1);
	}
	if( RNUM( ((SCOREC *)e1)->sc ) < RNUM( ((SCOREC *)e2)->sc) ) {
		return(-1);
	}
	return(0);
}
/****************************************************************************/
PRIV_I qsort_scorecs_decend(const void *e1, const void *e2)
{
	if( RNUM( ((SCOREC *)e1)->sc ) > RNUM( ((SCOREC *)e2)->sc) ) {
		return(-1);
	}
	if( RNUM( ((SCOREC *)e1)->sc ) < RNUM( ((SCOREC *)e2)->sc) ) {
		return(1);
	}
	return(0);
}
/**************************************************************************
*	Dump SCOREC values to file/stdout
*/
void DumpScorecsI(SCOREC *recsPO,int n,FILE *outPF)
{
	int i;

	HAND_NFILE(outPF);
	for(i=0;i<n;i++) {
		fprintf(outPF," [%4d]   Id %4d   Sc %8.5f\n",i,
			recsPO[i].id,recsPO[i].sc);
	}
}
/***************************************************************************
*	Set SCFIELD with nameS 
*/
int SetScfieldNameI(SCFIELD *sfPO,char *nameS,int max)
{
	VALIDATE(sfPO,SCFIELD_ID);
    max = (max < 0) ? NSIZE : max;
	LIMIT_NUM(max, 0, NSIZE);
	strncpy(sfPO->name,nameS,max);
	return(TRUE);
}
/***************************************************************************
*	Fills nameS with name of SCFIELD structure
*/
int FillScfieldNameStringI(SCFIELD *sfPO,char *nameS,int max)
{
	VALIDATE(sfPO,SCFIELD_ID);
    max = (max < 0) ? NSIZE : max;
	LIMIT_NUM(max, 0, NSIZE);
	strncpy(nameS,sfPO->name,max);
	return(TRUE);
}
/***********************************************************************
*	Evaluate score for value vR given settings in SCFIELD
*/
REAL EvalScfieldScoreR(SCFIELD *sfPO,REAL vR)
{
	DOUB scR,dxR,fxR,y1R,y2R,dyR;
	int x1,x2;

	DB_SCF DB_PrI(">> EvalScfieldScoreR %s (n=%d) %4.3f\n", sfPO->name,sfPO->n,vR);
	/***
	*	No transform case = default val
	*/
	if(sfPO->n < 2) {
		DB_SCF DB_PrI("<< EvalScfieldScoreR no func elements = 1\n");
		return(1.0);
	}
    /***
    *   Out of bounds or index cases
    */
	DB_SCF DB_PrI("+ min=%0.3f max=%0.3f step=%0.3f\n",
        sfPO->min, sfPO->max, sfPO->step);
    if(vR <= sfPO->min) {
	    DB_SCF DB_PrI("+ under [0]\n", sfPO->n,vR);
        scR = sfPO->yval[0];
    }
    else if(vR >= sfPO->max) {
	    DB_SCF DB_PrI("+ over [%d]\n", sfPO->n - 1);
        scR = sfPO->yval[sfPO->n - 1];
    }
    else {
        /***
        *   Interpolate 
        */
        dxR = (vR - sfPO->min) / sfPO->step;
		x1 = INT(dxR);
        x2 = x1 + 1;
	    DB_SCF DB_PrI("+ dx=%f X1=%d X2=%d\n",dxR,x1,x2);
        BOG_CHECK(x2 >= sfPO->n);
        y1R = sfPO->yval[x1];
        y2R = sfPO->yval[x2];
		dyR = y1R - y2R;
		fxR = dxR - RNUM(INT(dxR));
		DB_SCF DB_PrI("+ Y1=%f Y2=%f dY=%f fX=%f\n",y1R,y2R,dyR,fxR);
		scR = y1R - (dyR * fxR);
	}	
	DB_SCF DB_PrI("<< EvalScfieldScoreR %f\n",scR);
	return(scR);
}
