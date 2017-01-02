/*
* numlist.c
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
#include "prim.h"
#include "numlist.h"

#define DB_NLOW     if(DB[30])

/**************************************************************************/
NUMLIST *CreateNumlistPO(int type, char *fnameS, int n)
{
    int ok;
    NUMLIST *nlPO;

    DB_NLOW DB_PrI(">> CreateNumlistPO type=%d n=%d\n",type,n);
    if(!(nlPO=(NUMLIST *)ALLOC(1,sizeof(NUMLIST)))) {
        return(NULL);
    }
    nlPO->ID = NUMLIST_ID;
    InitNumlist(nlPO, type);
    if( ! NumlistTypeOkI(nlPO) ) {
        PROBLINE;
        printf("Numlist type not OK: called with %d\n",type);
        CHECK_NUMLIST(nlPO);
        return(NULL);
    }
    /***
    *   Getting values from file
    */
    if( fnameS && (!NO_S(fnameS))) {
        DB_NLOW DB_PrI("+ getting values from file |%s|...\n",fnameS);
        if(!NumlistLoadFromFileI(nlPO, type, fnameS, 1)) {
            CHECK_NUMLIST(nlPO);
            return(NULL);
        }
    }
    else if( n > 0 ) {
        DB_NLOW DB_PrI("+ not from file; Simply initializing to len %d\n",n);
        /***
        *   Setting length and initialize
        */
        ok = SetNumlistLengthI(nlPO, n);
        if(ok) {
            if(nlPO->type == IS_INT) {
                ok = SetNumlistRangeIntsI(nlPO, 0, n, 0);
            }
            else if (nlPO->type == IS_DOUB) {
                ok = SetNumlistRangeDoubsI(nlPO, 0, n, 0.0);
            }
        }
        if(!ok) {
            CHECK_NUMLIST(nlPO);
            return(NULL);
        }
    }
    DB_NLOW DB_PrI("<< CreateNumlistPO %p\n",nlPO);
    return(nlPO);
}
/*************************************************************************/
int DestroyNumlistI(NUMLIST *nlPO)
{
    DB_NLOW DB_PrI(">> DestroyNumlistI %p\n",nlPO);
    VALIDATE(nlPO,NUMLIST_ID);
    CHECK_FREE(nlPO->ivals);
    CHECK_FREE(nlPO->dvals);
    FREE(nlPO);
    DB_NLOW DB_PrI("<< DestroyNumlistI TRUE\n");
    return(TRUE);
}
/*************************************************************************/
void InitNumlist(NUMLIST *nlPO, int type)
{
    DB_NLOW DB_PrI(">> InitNumlist %p type=%d\n",nlPO,type);
    VALIDATE(nlPO,NUMLIST_ID);
    INIT_S(nlPO->fname);
    CHECK_FREE(nlPO->ivals);
    CHECK_FREE(nlPO->dvals);
    nlPO->n_vals = 0;
    nlPO->n = 0;
    nlPO->type = type;
    DB_NLOW DB_PrI("<< InitNumlist\n");
    return;
}
/*************************************************************************/
int NumlistTypeOkI(NUMLIST *nlPO)
{
    int ok;

    VALIDATE(nlPO,NUMLIST_ID);
    ok = FALSE;
    if( (nlPO->type == IS_INT) || (nlPO->type == IS_DOUB) ) {
        ok++;
    }
    return(ok);
}
/***************************************************************************
*   Sets length after making sure there's enough space for this
*/
int SetNumlistLengthI(NUMLIST *nlPO, int len)
{
    DB_NLOW DB_PrI(">> SetNumlistLengthI %p len=%d\n",nlPO,len);
    VALIDATE(nlPO,NUMLIST_ID);
    if( !HandleNumlistSpaceI(nlPO, len) ) {
        DB_NLOW DB_PrI("<< SetNumlistLengthI FALSE\n");
        return(FALSE);
    }
    nlPO->n = len;
    DB_NLOW DB_PrI("<< SetNumlistLengthI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Check if passed index will fit in allocated space 
*   Allocate space if needed
*/
int HandleNumlistSpaceI(NUMLIST *nlPO, int i)
{
    int nm;

    DB_NLOW DB_PrI(">> HandleNumlistSpaceI %p i=%d\n",nlPO,i);
    VALIDATE(nlPO,NUMLIST_ID);
    /***
    *   Index needs space so alloc
    */
    if(i >= nlPO->n_vals) {
        nm = MAX_NUM(i, (nlPO->n_vals + ALLOC_BLOCK));
        if(nlPO->type == IS_INT) {
            if(nlPO->n_vals == 0) {
                nlPO->ivals = (int *)ALLOC(nm ,sizeof(int));
                DB_NLOW DB_PrI("+ Allocated nm=%d ints %p\n",nm,nlPO->ivals);
            }
            else {
                nlPO->ivals = (int *)REALLOC(nlPO->ivals, nm ,sizeof(int));
                DB_NLOW DB_PrI("+ ReAllocated nm=%d ints %p\n",nm,nlPO->ivals);
            }
            if(!nlPO->ivals) {
                PROBLINE;
                printf("HandleNumlistSpaceI Failed allocate %d ints\n",nm);
                return(FALSE);
            }
        }
        else if(nlPO->type == IS_DOUB) {
            if(nlPO->n_vals == 0) {
                nlPO->dvals = (DOUB *)ALLOC(nm ,sizeof(DOUB));
                DB_NLOW DB_PrI("+ Allocated nm=%d doubs %p\n",nm,nlPO->dvals);
            }
            else {
                nlPO->dvals = (DOUB *)REALLOC(nlPO->dvals, nm ,sizeof(DOUB));
                DB_NLOW DB_PrI("+ ReAllocated nm=%d doubs %p\n",nm,nlPO->dvals);
            }
            if(!nlPO->dvals) {
                PROBLINE;
                printf("HandleNumlistSpaceI Failed allocate %d doubs\n",nm);
                return(FALSE);
            }
        }
        else {
            DB_NLOW DB_PrI("<< HandleNumlistSpaceI type=%d FALSE\n",nlPO->type);
            return(FALSE);
        }
        nlPO->n_vals = nm;
    }
    /***
    *   Number = high water mark = max of existing or new index
    */
    DB_NLOW DB_PrI("<< HandleNumlistSpaceI n_vals=%d TRUE\n",nlPO->n_vals);
    return(TRUE);
}
/*************************************************************************/
NUMLIST *DuplicateNumlistPO(NUMLIST *nlPO, int ntype)
{
    int i,v,len,type;
    char nameS[NSIZE],fnameS[NSIZE];
    DOUB vD;
    NUMLIST *newPO;

    DB_NLOW DB_PrI(">> DuplicateNumlistPO given %p ntype=%d\n",nlPO,ntype);
    VALIDATE(nlPO,NUMLIST_ID);
    len = GetNumlistLengthI(nlPO);
    type = GetNumlistTypeI(nlPO);
    DB_NLOW DB_PrI("+ len=%d type=%d\n",len,type);
    /***
    *   New one after first checking new type / keeping source type
    */
    if( (ntype != IS_INT) && (ntype != IS_DOUB) ) {
        ntype = type;
    }
    if(!(newPO = CreateNumlistPO(ntype,NULL,len))) {
        DB_NLOW DB_PrI("<< DuplicateNumlistPO NULL\n");
        return(NULL);
    }
    DB_NLOW DB_PrI("+ copying names and values\n");
    GetNumlistNamesI(nlPO, nameS, fnameS, NSIZE);
    SetNumlistNamesI(newPO, nameS, fnameS, NSIZE);
    if(type == IS_INT) {
        for(i=0;i<len;i++) 
        {
            GetNumlistIntI(nlPO,i,&v);
            if(ntype == IS_INT) {
                SetNumlistIntI(newPO,i,v);
            }
            else {
                vD = DNUM(v);
                SetNumlistDoubI(newPO,i,vD);
            }
        }
    }
    else if(type == IS_DOUB) {
        for(i=0;i<len;i++) 
        {
            GetNumlistDoubI(nlPO,i,&vD);
            if(ntype == IS_DOUB) {
                SetNumlistDoubI(newPO,i,vD);
            }
            else {
                v = INT(vD);
                SetNumlistIntI(newPO,i,v);
            }
        }
    }
    DB_NLOW DB_PrI("<< DuplicateNumlistPO %p\n",newPO);
    return(newPO);
}
/*************************************************************************/
void DumpNumlist(NUMLIST *nlPO, int st, int en, FILE *outPF)
{
    int n;

    VALIDATE(nlPO,NUMLIST_ID);
    HAND_NFILE(outPF);
    fprintf(outPF,"Numlist at %p\n",nlPO);
    fprintf(outPF,"Name  |%s|\n",nlPO->name);
    fprintf(outPF,"Fname |%s|\n",nlPO->fname);
    n = NumlistGoodStartEndI(nlPO, st, &st, en, &en);
    if(nlPO->type == IS_INT) {
        fprintf(outPF,"Contains %d integers\n",n);
        if(en > st) {
            fprintf(outPF,"Dumping %d to %d\n",st,en);
            DumpArray(nlPO->ivals, IS_INT, st, en, NULL, outPF);
        }
    }
    else if(nlPO->type == IS_DOUB) {
        fprintf(outPF,"Contains %d doubles\n",n);
        if(en > st) {
            fprintf(outPF,"Dumping %d to %d\n",st,en);
            DumpArray(nlPO->dvals, IS_DOUB, st, en, NULL, outPF);
        }
    }
    else {
        fprintf(outPF,"No values; type=%d\n",nlPO->type);
    }
    return;
}
/***************************************************************************/
int GetNumlistTypeI(NUMLIST *nlPO)
{
    VALIDATE(nlPO,NUMLIST_ID);
    return(nlPO->type);
}
/***************************************************************************/
int GetNumlistLengthI(NUMLIST *nlPO)
{
    VALIDATE(nlPO,NUMLIST_ID);
    return(nlPO->n);
}
/****************************************************************************
*   Add value to end of numlist, growing it by one
*/
int AppendNumlistIntI(NUMLIST *nlPO, int v)
{
    int i;
    
    i = GetNumlistLengthI(nlPO);
    return(AddNumlistIntI(nlPO, i, v));
}
/***************************************************************************/
int AppendNumlistDoubI(NUMLIST *nlPO, DOUB vD)
{
    int i;
    
    i = GetNumlistLengthI(nlPO);
    return(AddNumlistDoubI(nlPO, i, vD));
}
/***************************************************************************
*   Set usable start end values; Set bounded values and return number of values
*   Negative indices mean set to (existing) bounds 
*/
int NumlistGoodStartEndI(NUMLIST *nlPO, int st, int *stPI, int en, int *enPI)
{
    int n;

    n = GetNumlistLengthI(nlPO);
    st = (st<0) ? 0 : st;
    en = (en<0) ? n : en;
    LIMIT_NUM(st, 0, n);
    LIMIT_NUM(en, 0, n);
    if(stPI) {
        *stPI = st;
    }
    if(enPI) {
        *enPI = en;
    }
    return(n);
}
/**************************************************************************
*   Set range of numlist values from start to end to value v
*   Since calls AddNumlistInt, this can add space to list if needed
*/
int SetNumlistRangeIntsI(NUMLIST *nlPO, int st, int en, int v)
{
    int i, n;

    VALIDATE(nlPO,NUMLIST_ID);
    /**
    *   Note that length can be zero if not set previously 
    */
    n = NumlistGoodStartEndI(nlPO, st, &st, en, &en);
    if(en > n) {
        if(!SetNumlistLengthI(nlPO, en)) {
            return(FALSE);
        }
    }
    for(i=st; i<en; i++) {
        if(!SetNumlistIntI(nlPO, i, v)) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/**************************************************************************/
int SetNumlistRangeDoubsI(NUMLIST *nlPO, int st, int en, DOUB vD)
{
    int i, n;

    VALIDATE(nlPO,NUMLIST_ID);
    n = NumlistGoodStartEndI(nlPO, st, &st, en, &en);
    if(en > n) {
        if(!SetNumlistLengthI(nlPO, en)) {
            return(FALSE);
        }
    }
    for(i=st; i<en; i++) {
        if(!SetNumlistDoubI(nlPO, i, vD)) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/***************************************************************************
*   Add int to num list, updating length if needed
*/
int AddNumlistIntI(NUMLIST *nlPO, int i, int v)
{
    int n;

    n = GetNumlistLengthI(nlPO);
    if( i >= n ) {
        if( !SetNumlistLengthI(nlPO, i+1) ) {
            return(FALSE);
        }
    }
    SetNumlistIntI(nlPO, i, v);
    return(TRUE);
}
/**************************************************************************/
int AddNumlistDoubI(NUMLIST *nlPO, int i, DOUB vD)
{
    int n;

    n = GetNumlistLengthI(nlPO);
    if( i >= n ) {
        if( !SetNumlistLengthI(nlPO, i+1) ) {
            return(FALSE);
        }
    }
    SetNumlistDoubI(nlPO, i, vD);
    return(TRUE);
}
/***************************************************************************
*   Set value in numlist ONLY if it already fits in the existing length
*/
int SetNumlistIntI(NUMLIST *nlPO, int i, int v)
{
    VALIDATE(nlPO,NUMLIST_ID);
    if( (i < 0) || (i >= nlPO->n) || (nlPO->type != IS_INT) ) {
        return(FALSE);
    }
    nlPO->ivals[i] = v;
    return(TRUE);
}
/***************************************************************************/
int SetNumlistDoubI(NUMLIST *nlPO, int i, DOUB vD)
{
    VALIDATE(nlPO,NUMLIST_ID);
    if( (i < 0) || (i >= nlPO->n) || (nlPO->type != IS_DOUB) ) {
        return(FALSE);
    }
    nlPO->dvals[i] = vD;
    return(TRUE);
}
/***************************************************************************
*   Get value form numlist; If int pointer is supplied, value is set here
*/
int GetNumlistIntI(NUMLIST *nlPO, int i, int *vPI)
{
    VALIDATE(nlPO,NUMLIST_ID);
    if( (i < 0) || (i >= nlPO->n) || (nlPO->type != IS_INT) ) {
        return(FALSE);
    }
    if(vPI) {
        *vPI = nlPO->ivals[i];
    }
    return(TRUE);
}
/***************************************************************************/
int GetNumlistDoubI(NUMLIST *nlPO, int i, DOUB *vPD)
{
    VALIDATE(nlPO,NUMLIST_ID);
    if( (i < 0) || (i >= nlPO->n) || (nlPO->type != IS_DOUB) ) {
        return(FALSE);
    }
    if(vPD) {
        *vPD = nlPO->dvals[i];
    }
    return(TRUE);
}
/***************************************************************************
*   Get numlist int and doub, regardless of actual type
*/
int GetNumlistIntDoubI(NUMLIST *nlPO, int i, int *vPI, DOUB *vPD)
{
    DOUB vD;
    int v;

    VALIDATE(nlPO,NUMLIST_ID);
    if( (i < 0) || (i >= nlPO->n) ) {
        return(FALSE);
    }
    if(nlPO->type == IS_DOUB) {
        vD = nlPO->dvals[i];
        v = INT(vD);
    }
    else {
        v = nlPO->ivals[i];
        vD = DNUM(v);
    }
    if(vPI) {
        *vPI = v;
    }
    if(vPD) {
        *vPD = vD;
    }
    return(TRUE);
}
/***************************************************************************
*   Get a pointer to the ints contained in numlist; Can also pt to len 
*/
int GetNumlistPtrIntsI(NUMLIST *nlPO, int **vPPI, int *nPI)
{
    VALIDATE(nlPO,NUMLIST_ID);
    if( nlPO->type != IS_INT ) {
        return(FALSE);
    }
    if(vPPI) {
        *vPPI = nlPO->ivals;
    }
    if(nPI) {
        *nPI = GetNumlistLengthI(nlPO);
    }
    return(TRUE);
}
/***************************************************************************/
int GetNumlistPtrDoubsI(NUMLIST *nlPO, DOUB **vPPD, int *nPI)
{
    VALIDATE(nlPO,NUMLIST_ID);
    if( nlPO->type != IS_DOUB ) {
        return(FALSE);
    }
    if(vPPD) {
        *vPPD = nlPO->dvals;
    }
    if(nPI) {
        *nPI = GetNumlistLengthI(nlPO);
    }
    return(TRUE);
}
/*************************************************************************/
int SetNumlistNamesI(NUMLIST *nlPO, char *nameS, char *fnameS, int max)
{
    int n;

    VALIDATE(nlPO,NUMLIST_ID);
    if(nameS) {
        n = (max < 0) ? strlen(nameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(nlPO->name,nameS,n);
        nlPO->name[n] = '\0';
    }
    if(fnameS) {
        n = (max < 0) ? strlen(fnameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(nlPO->fname,fnameS,n);
        nlPO->fname[n] = '\0';
    }
    return(TRUE);
}
/*************************************************************************/
int GetNumlistNamesI(NUMLIST *nlPO, char *nameS, char *fnameS, int max)
{
    int n;

    VALIDATE(nlPO,NUMLIST_ID);
    if(nameS) {
        n = (max < 0) ? strlen(nameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(nameS,nlPO->name,n);
        nameS[n] = '\0';
    }
    if(fnameS) {
        n = (max < 0) ? strlen(fnameS) : max;
        LIMIT_NUM(n,0,NSIZE);
        strncpy(fnameS,nlPO->fname,n);
        fnameS[n] = '\0';
    }
    return(TRUE);
}
/*************************************************************************/
int NumlistAutoFormatStringI(NUMLIST *nlPO, int *pPI, int *wPI, char *formS)
{
    DOUB lohiD[2], *valsPD;
    int num,type;

    type = GetNumlistTypeI(nlPO);
    if(type == IS_INT) {
        if(!NumlistStatsI(nlPO,-1,-1, &lohiD[0], &lohiD[1], NULL,NULL)) {
            return(FALSE);
        }
        valsPD = lohiD;
        num = 2;
    }
    else {
        if(!GetNumlistPtrDoubsI(nlPO, &valsPD, &num)) {
            return(FALSE);
        }
    }
    DoubArrayPrecisionI(valsPD, num, pPI, wPI, formS);
    return(TRUE);
}
/*************************************************************************/
int NumlistSameTypeI(NUMLIST *flPO, NUMLIST *slPO)
{
    int same;

    VALIDATE(flPO,NUMLIST_ID);
    VALIDATE(slPO,NUMLIST_ID);
    same = TRUE;
    if( GetNumlistTypeI(flPO) != GetNumlistTypeI(slPO) ) {
        same = FALSE;
    }
    return(same);
}
/*************************************************************************/
int NumlistSameLenI(NUMLIST *flPO, NUMLIST *slPO)
{
    int same;

    VALIDATE(flPO,NUMLIST_ID);
    VALIDATE(slPO,NUMLIST_ID);
    same = TRUE;
    if( GetNumlistLengthI(flPO) != GetNumlistLengthI(slPO) ) {
        same = FALSE;
    }
    return(same);
}
/*************************************************************************/
int NumlistAreSameI(NUMLIST *flPO, NUMLIST *slPO)
{
    return( NumlistSameTypeI(flPO, slPO) && NumlistSameLenI(flPO, slPO) );
}
/**************************************************************************/
int SortNumlistI(NUMLIST *nlPO, int dir)
{
    int ok, n, *vPI;
    DOUB *vPD;

    VALIDATE(nlPO,NUMLIST_ID);
    ok = FALSE;
    if(nlPO->type == IS_INT) {
        if(GetNumlistPtrIntsI(nlPO, &vPI, &n)) {
            SortArray(vPI,IS_INT,n,dir);
            ok++;
        }
    }
    else if(nlPO->type == IS_DOUB) {
        if(GetNumlistPtrDoubsI(nlPO, &vPD, &n)) {
            SortArray(vPD,IS_DOUB,n,dir);
            ok++;
        }
    } 
    return(ok);
}
/***************************************************************************
*   Check if is sorted in given direction; Consecutive values must differ by mindif 
*/
int NumlistIsSortedI(NUMLIST *nlPO, int dir, DOUB mindifD)
{
    int i;
    DOUB v1D, v2D;

    VALIDATE(nlPO,NUMLIST_ID);
    if( GetNumlistLengthI(nlPO) < 2) {
        return(TRUE);
    }
    GetNumlistIntDoubI(nlPO, 0, NULL, &v1D);
    i = 1; 
    while( GetNumlistIntDoubI(nlPO, i++, NULL, &v2D) )
    {
        if(dir < 0) {
            if( (v1D - v2D) < mindifD) {
                return(FALSE);
            }
        }
        else {
            if( (v2D - v1D) < mindifD) {
                return(FALSE);
            }
        }
        v1D = v2D;
    }
    return(TRUE);
}
/***************************************************************************
*   Smoothing (sliding window) of first numlist into second numlist
*/
int SmoothNumlistI(NUMLIST *flPO, NUMLIST *slPO, int win)
{
    int ok, n, m, *fPI, *sPI;
    DOUB *fPD, *sPD;
    
    VALIDATE(flPO,NUMLIST_ID);
    VALIDATE(slPO,NUMLIST_ID);
    if( !NumlistSameTypeI(flPO, slPO) ) {
        return(FALSE);
    }
    /***
    *   Make sure second (smoothed) list is same size as first; 
    *   Setting length with allocate space if needed
    */
    n = GetNumlistLengthI(flPO);
    m = GetNumlistLengthI(slPO);
    if( m != n ) {
        if( !SetNumlistLengthI(slPO,n) ) {
            return(FALSE);
        }
    }
    /***
    *   Get pointers to raw arrays then smooth with these
    */
    ok = FALSE;
    if(flPO->type == IS_INT) {
        GetNumlistPtrIntsI(flPO, &fPI, &n);
        GetNumlistPtrIntsI(slPO, &sPI, NULL);
        SmoothArray(fPI,sPI,IS_INT,n,win);
        ok++;
    }
    else if(flPO->type == IS_DOUB) {
        GetNumlistPtrDoubsI(flPO, &fPD, &n);
        GetNumlistPtrDoubsI(slPO, &sPD, NULL);
        SmoothArray(fPD,sPD,IS_DOUB,n,win);
        ok++;
    } 
    return(ok);
}
/***************************************************************************/
int NumlistLoadFromFileI(NUMLIST *nlPO, int type, char *fnameS, int col)
{
    DB_NLOW DB_PrI(">> NumlistLoadFromFileI type=%d |%s| col=%d\n",type,fnameS,col);
    NOT_YET;
    DB_NLOW DB_PrI("<< NumlistLoadFromFileI FALSE\n");
    return(FALSE);
}
/***************************************************************************/
int NumlistSumI(NUMLIST *nlPO, int st, int en, DOUB *sumPD)
{
    int i,n;
    DOUB vD,sumD;

    VALIDATE(nlPO,NUMLIST_ID);
    NumlistGoodStartEndI(nlPO, st, &st, en, &en);
    sumD = 0.0;
    n = 0;
    for(i=st; i<en; i++)
    {
        GetNumlistIntDoubI(nlPO, i, NULL, &vD);
        sumD += vD;
        n++;
    }
    if(sumPD) {
        *sumPD = sumD;
    }
    return(n);
}
/***************************************************************************
*   Numlist interface to raw stat functions
*/
int NumlistStatsI(NUMLIST *nlPO, int st, int en, DOUB *loD, DOUB *hiD,
    DOUB *avD, DOUB *sdD)
{
    int ok,type;

    NumlistGoodStartEndI(nlPO, st, &st, en, &en);
    type = GetNumlistTypeI(nlPO);
    ok = FALSE;
    if(type == IS_INT) {
        ok = ArrayStatsI(nlPO->ivals,IS_INT,st,en,loD,hiD,avD,sdD);
    }
    else if(type == IS_DOUB) {
        ok = ArrayStatsI(nlPO->dvals,IS_DOUB,st,en,loD,hiD,avD,sdD);
    }
    return(ok);
}

