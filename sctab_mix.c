/*
* sctab_mix.c
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
* See https://www.verdascend.com/ for more
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "prim.h"
#include "numlist.h"
#include "wordlist.h"
#include "score.h"
#include "table.h"
#include "scoretab.h"


/*******************************************************************************
*   Handle merging options
*/
int HandleSctTableMergingI(SCORETAB *stPO)
{
    int ok;
    char bufS[DEF_BS],tabS[DEF_BS];
    FILE *fPF;
    TABLE *tabPO;

    /***
    *   GS cross over?
    */
    if(stPO->do_gax) {
        return(HandleSctGACrossOverI(stPO));
    }
    /***
    *   Table merging
    */
    ok = TRUE;
    tabPO = NULL;
    if(stPO->do_mlis) {
        if(!(fPF= OpenUFilePF(stPO->mergname,"r",NULL))) {
            return(FALSE);
        }
        while((fgets(bufS,LINEGRAB,fPF)) != NULL) {
            if(COM_LINE(bufS)) {
                continue;
            }
            INIT_S(tabS);
            sscanf(bufS,"%s",tabS);
            if(NO_S(tabS)) {
                continue;
            }
            CHECK_TABLE(tabPO);
            if(!GetTableI(tabS,stPO->rlab,stPO->clab,stPO->corn,stPO->do_skp,&tabPO)) {
                printf("Problem loading table: %s\n",tabS);
                ok = FALSE;
                break;
            }
            if(!HandleTableMergeI(stPO,tabPO,tabS)) {
                ok = FALSE;
                break;
            }   
        }
        CHECK_FILE(fPF);
    }
    else {
        if(!GetTableI(stPO->mergname,stPO->rlab,stPO->clab,stPO->corn,stPO->do_skp,&tabPO)) {
            printf("Problem loading table: %s\n",stPO->mergname);
            ok = FALSE;
        }
        else {
            if(!HandleTableMergeI(stPO,tabPO,stPO->mergname)) {
                ok = FALSE;
            }   
        }
    }
    CHECK_TABLE(tabPO);
    return(ok);
}
/*****************************************************************************
*
*/
int HandleTableMergeI(SCORETAB *stPO, TABLE *tabPO, char *tabS)
{
    int mix, rowmatch;
    char mixS[DEF_BS];
    TABLE *newPO;

    /* check row names match */
    rowmatch = TRUE;
    /***    
    *   Append?
    */
    if(stPO->do_mapc) {
        if(!HandleTableClabPreSufI(tabPO, stPO->macpre, stPO->macsuf)) {
            printf("Problem adding prefix / suffix to merge table...?\n");
            return(FALSE);
        }
        newPO = NewPastedTablePO(stPO->tab,tabPO,TRUE,rowmatch,stPO->do_kc);
        if(!newPO) {
            printf("Problem appending tables together\n");
            return(FALSE);
        }
        /* Replace current table with new one */
        CHECK_TABLE(stPO->tab);
        stPO->tab = newPO;
        SetUpSctTableSpaceI(stPO);
    }
    /***
    *   Simply mix values
    */
    else {
        mix = FigureMergeMixI(stPO, mixS);
        if(!MixTablesI(stPO->tab,tabPO,stPO->tab,TRUE,mix,TRUE)) {
            return(FALSE);
        }
        if(!stPO->quiet) {
            printf("# Table %s %s input table\n",tabS,mixS);
        }
    }
    return(TRUE);
}
/****************************************************************************
*
*/
int HandleTableClabPreSufI(TABLE *tabPO, char *preS, char *sufS)
{
    int c,nc,mod;
    char nameS[NSIZE], newS[NSIZE], pS[NSIZE], sS[NSIZE];

    /***
    *   Anything to do?
    */
    mod = 0;
    INIT_S(pS); 
    if((preS) && (!NO_S(preS))) {
        sprintf(pS,"%s",preS);
        mod++;
    }
    INIT_S(sS);
    if((sufS) && (!NO_S(sufS))) {
        sprintf(sS,"%s",sufS);
        mod++;
    }
    if(!mod) {
        return(TRUE);
    }
    /***
    *   Update each label
    */
    nc = GetTableColsI(tabPO,FALSE);
    for(c=0;c<nc;c++) 
    {
        GetTableColLabI(tabPO,c,nameS,-1);
        sprintf(newS,"%s%s%s",pS,nameS,sS);
        SetTableColLabI(tabPO,c,newS);
    }
    return(TRUE);
}
/*****************************************************************************
*   Decide what merging option is appropriate and tell story if given string
*/
int FigureMergeMixI(SCORETAB *stPO, char *whatS)
{
    int mix;

    mix = MATH_ADD;
    if(whatS)
        sprintf(whatS,"added");
    if(stPO->do_msub) {
        if(whatS)
            sprintf(whatS,"subtracted");
        mix = MATH_SUB;
    }
    else if(stPO->do_mmul) {
        if(whatS)
            sprintf(whatS,"multipled");
        mix = MATH_MUL;
    }
    else if(stPO->do_mdiv) {
        if(whatS)
            sprintf(whatS,"divided");
        mix = MATH_DIV;
    }
    else if(stPO->do_mmin) {
        if(whatS)
            sprintf(whatS,"Min values");
        mix = MATH_MIN;
    }
    else if(stPO->do_mmax) {
        if(whatS)
            sprintf(whatS,"Max values");
        mix = MATH_MAX;
    }
    return(mix);
}
/*************************************************************************
*   Copy to new transposed one and kill original
*/
int HandleTransposeI(TABLE *tabPO, TABLE **tranPPO)
{
    int ok;

    /***
    *   Copy with masking = TRUE, transpose = TRUE
    */
    ok = CopyTableI(tabPO,TRUE,TRUE,tranPPO);
    CHECK_TABLE(tabPO);
    return(ok);
}
/*************************************************************************
*   Sort row values in table; Same row lable, just columns permuted
*/
int HandleRowSortI(TABLE *tabPO, int sdir, int mask)
{
    int r,nc,nr;
    NUMLIST *nlPO;
    DOUB *valsPD;

    /***    
    *   Local array (via numlist) for row
    */
    nc = GetTableColsI(tabPO,mask);
    if ( nc < 2 ) {
        return(nc);
    }
    nlPO = CreateNumlistPO(IS_DOUB, NULL, nc);
    if ( ! nlPO ) {
        return(FALSE);
    }
    /***
    *   Each row
    */
    nr = GetTableRowsI(tabPO,mask);
    for(r=0; r<nr; r++) {
        if( mask && (!tabPO->rmask[r]) ) {
            continue;
        }
        /***
        *   Values from table to numlist, then get raw array pointer from that,
        *    sort the array, and set values back into table
        */
        GetTableRowValsI(tabPO, r, nlPO, mask);
        GetNumlistPtrDoubsI(nlPO, &valsPD, NULL);
        SortArray(valsPD, IS_DOUB, nc, sdir);
        SetTableRowValsI(tabPO, r, nlPO, nc, mask);
    }
    CHECK_NUMLIST(nlPO);
    return(nc);
}
