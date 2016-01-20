/*
* score.c
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
*   Allocate n SCOREC structures
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
*   Free SCOREC structure(s)
*/
int DestroyScorecsI(SCOREC *recsPO)
{
    BOG_CHECK(!recsPO);
    FREE(recsPO);
    return(TRUE);
}
/**************************************************************************
*   Initialize SCOREC(s) with indices and values 
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
*   Create an SCFIELD structure
*   Only substructs / name is initialized; No X,Y values
*/
SCFIELD *CreateScfieldPO(char *nameS)
{
    SCFIELD *sfPO;

    DB_SCM DB_PrI(">> CreateScfieldPO\n");
    if( ! (sfPO = (SCFIELD *)ALLOC(1,sizeof(SCFIELD)) ) ) {
        return(NULL);
    }
    sfPO->ID = SCFIELD_ID;
    sfPO->xvals = CreateNumlistPO(IS_DOUB,NULL,0);
    sfPO->yvals = CreateNumlistPO(IS_DOUB,NULL,0);
    SetNumlistNamesI(sfPO->xvals,"Scfield-Xvals",NULL,NSIZE);
    SetNumlistNamesI(sfPO->yvals,"Scfield-Yvals",NULL,NSIZE);
    sfPO->minx = sfPO->maxx = 0.0;
    if(nameS) {
        SetScfieldNameI(sfPO,nameS,NSIZE);
    }
    DB_SCM DB_PrI("<< CreateScfieldPO %p\n",sfPO);
    return(sfPO);
}
/**************************************************************************
*   Frees SCFIELD datatype 
*/
int DestroyScfieldI(SCFIELD *sfPO)
{
    VALIDATE(sfPO,SCFIELD_ID);
    CHECK_NUMLIST(sfPO->xvals);
    CHECK_NUMLIST(sfPO->yvals);
    FREE(sfPO);
    return(TRUE);
}
/****************************************************************************
*   Sort SCOREC list[n] by value
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
*   Sort SCOREC list[n] by index
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
*   Dump SCOREC values to file/stdout
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
*   Set SCFIELD with nameS 
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
*   Fills nameS with name of SCFIELD structure
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
*   Evaluate score for value (i.e. Y for X) given settings in SCFIELD
*   Assumes X values are increasing and sorted for correct behavior
*/
REAL EvalScfieldScoreR(SCFIELD *sfPO, REAL vR)
{
    DOUB xval1D, xval2D, yval1D, yval2D, yD;
    int xi;

    VALIDATE(sfPO,SCFIELD_ID);
    DB_SCF DB_PrI(">> EvalScfieldScoreR %s %4.3f\n", sfPO->name,vR);
    /***
    *   Out of bound X = min / max Y index value
    */
    if(vR <= sfPO->minx) {
        GetNumlistDoubI(sfPO->yvals, 0, &yD);
    }
    else if(vR >= sfPO->maxx) {
        GetNumlistDoubI(sfPO->yvals, sfPO->n - 1, &yD);
    }
    else {
        BOG_CHECK(sfPO->n < 2);
        /***
        *   Interpolate Y after finding bounding X indices
        */
        xi = 0;
        GetNumlistDoubI(sfPO->xvals, xi, &xval2D);
        while( xval2D < vR ) 
        {
            xi++;
            GetNumlistDoubI(sfPO->xvals, xi, &xval2D);
        }
        GetNumlistDoubI(sfPO->xvals, xi - 1, &xval1D);
        /* Y values */
        GetNumlistDoubI(sfPO->yvals, xi, &yval2D);
        GetNumlistDoubI(sfPO->yvals, xi - 1, &yval1D);
        DB_SCF DB_PrI("+ x1=%f x2=%f\ty1=%f y2=%f\n",xval1D,xval2D,yval1D,yval2D);
        yD = yval1D + (yval2D-yval1D) * (vR-xval1D) / (xval2D-xval1D);
    }
    DB_SCF DB_PrI("<< EvalScfieldScoreR %f\n",yD);
    return(yD);
}
