/*
* table_io.c
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
#include <string.h>
#include <ctype.h>
#include "prim.h"
#include "table.h"

#define DB_TAB  if(DB[70])

/****************************************************************************
*   Attempt to read in a table from file named nameS
*   Flags rlab and clab indicate number of rows, cols to treat as labels
*/
int GetTableI(char *nameS, int rlab, int clab, int corn, int skp, TABLE **tPPO)
{
    int ok;
    FILE *fPF;
    TABLE *tabPO;

    DB_TAB DB_PrI(">> GetTableI |%s|\n",nameS);
    if(!(fPF=OpenUFilePF(nameS,"r",NULL))) {
        printf("Can't open table file: %s\n",nameS);
        return(FALSE);
    }
    ok = ParseTableI(fPF,rlab,clab,corn,skp,&tabPO);
    CHECK_FILE(fPF);
    /***
    *   Set input filename and name (if not set already)
    */
    if(ok) {
        sprintf(tabPO->fname,"%s",nameS);
        if(NO_S(tabPO->name)) {
            GetFilePartsI(nameS,NULL,tabPO->name,NULL);
        }
        *tPPO = tabPO;
    }
    DB_TAB DB_PrI("<< GetTableI %d\n",ok);
    return(ok);
}
/****************************************************************************
*   Attempt to load table from file
*   Flags rlab and clab specify if there are row and/or col lables
*   If skp then simply skip problem lines rather than bail out
*/
int ParseTableI(FILE *inPF, int rlab, int clab, int corn, int skp, TABLE **tPPO)
{
    char bufS[TABBUFF_SIZE+1], nameS[NSIZE], cornS[NSIZE];
    int line,nrow,ncol,n,tot,ok;
    TABLE *tabPO;
    NUMLIST *valsPO;
    WORDLIST *rlabsPO, *clabsPO;

    DB_TAB DB_PrI(">> ParseTableI rlab=%d clab=%d\n",rlab,clab);
    tabPO = NULL;
    valsPO = NULL;
    rlabsPO = clabsPO = NULL;
    INIT_S(cornS);
    /***
    *   Create number and word lists to collect values and labels
    */
    valsPO = CreateNumlistPO(IS_DOUB, NULL, 0);
    if(rlab) {
        rlabsPO = CreateWordlistPO(NULL, -1);
    }
    if(clab) {
        clabsPO = CreateWordlistPO(NULL, -1);
    }
    if( (!valsPO) || ( clab && (!clabsPO)) || ( rlab && (!rlabsPO)) ) {
        DB_TAB DB_PrI("<< ParseTableI vals=%p cols=%p rows=%p FALSE\n",valsPO,clabsPO,rlabsPO);
        CHECK_NUMLIST(valsPO);
        CHECK_WORDLIST(rlabsPO);
        CHECK_WORDLIST(clabsPO);
        return(FALSE); 
    }
    /***
    *   Get to first real line of file
    */
    DB_TAB DB_PrI("+ Getting to first row... ");
    line = 0;
    while(fgets(bufS,TABBUFF_SIZE,inPF)) {
        line++;
        if( (COM_LINE(bufS)) || BlankStringI(bufS) ) {
            continue;
        }
        break;
    }
    ok = (NO_S(bufS)) ? FALSE : TRUE;
    DB_TAB DB_PrI(" ok = %d\n",ok);
    /***
    *   Now collect labels or values from first row and count columns
    */
    ncol = nrow = tot = 0;
    if(ok) {
        if(clab) {
            ncol = ParseLabelRowFromBuffI(bufS, line, corn, cornS, clabsPO);
        }
        else {
            ncol = ParseDataRowFromBuffI(bufS, line, nrow++, rlabsPO, valsPO);
            tot = ncol;
        }
    }
    ok = (ncol > 0) ? TRUE : FALSE;
    /***
    *   Rest of lines or until specal END keyword
    */
    while( ok && (fgets(bufS,TABBUFF_SIZE,inPF)) ) {
        line++;
        if( (COM_LINE(bufS)) || BlankStringI(bufS) ) {
            continue;
        }
        /***
        *   Check for end/stop keywords
        */
        INIT_S(nameS);
        sscanf(bufS,"%s",nameS);
        Upperize(nameS);
        if(EQSTRING(nameS,"ENDTAB",6)) {
            break;
        }
        /***
        *   Legit row?
        */
        n = ParseDataRowFromBuffI(bufS, line, nrow, rlabsPO, valsPO);
        if( n != ncol ) {
            if(skp) {
                WARNLINE;
                Chomp(bufS);
                printf("# Got %d values from line %d when expecting %d: |%s|\n",n,line,ncol,bufS);
                /***
                *   Ignore = reset already-collected values / maybe row lable
                */
printf("xxx Setting numlist len to %d, rlabs %d\n",tot,nrow);
                SetNumlistLengthI(valsPO,tot);
                if(rlabsPO) {
                    SetWordlistLengthI(rlabsPO,nrow);
                }
            }
            else {
                PROBLINE;
                printf("Got %d values from line %d when expecting %d: |%s|\n",n,line,ncol,bufS);
                nrow = 0;
                break;
            }
        }
        else {
            nrow++;
            tot += n;
        }
    }
    ok = (tot > 0) ? TRUE : FALSE;
    /***
    *   If OK, allocate (empty) table then attach values and labels to it
    */
    if(ok) {
        DB_TAB DB_PrI("+ Allocating empty table structure\n");
        if( !(tabPO = CreateTablePO(-1,-1)) ) {
            PROBLINE;
            printf("Failed to allocate for table %d row X %d col\n",nrow,ncol);
            ok = FALSE;
        }
    }
    if(ok) {
        DB_TAB DB_PrI("+ Setting values into table structure\n");
        tabPO->nrow = nrow;
        tabPO->ncol = ncol;
        tabPO->vals = valsPO;       /* et values to what we've collected */
        tabPO->rlab = rlab;
        tabPO->clab = clab;
        tabPO->rlabs = rlabsPO;
        tabPO->clabs = clabsPO;
        sprintf(tabPO->corner, "%s", cornS);
        DB_TAB DB_PrI("+ Finishing off auxillary stuff and initilizing\n");
        ok = TableAllocAuxI(tabPO);
        InitTableI(tabPO, FALSE);
    }
    /***
    *   Set pointer to new table or clean up if problem
    */
    if(ok) {
        DB_TAB DB_PrI("+ Good so setting pointer to table: %p\n",tabPO);
        *tPPO = tabPO;
    }
    else {
        DB_TAB DB_PrI("+ Problems, cleaning up\n");
        CHECK_NUMLIST(valsPO);
        CHECK_WORDLIST(rlabsPO);
        CHECK_WORDLIST(clabsPO);
        CHECK_TABLE(tabPO);
    }
    DB_TAB DB_PrI("<< ParseTableI %d\n",ok);
    return(ok);
}
/*************************************************************************
*   Loads column lables from passed string; If rlab, put first token in cornS
*   Returns the number of col lables
*/
int ParseLabelRowFromBuffI(char *bufS, int line, int corn, char *cornS, WORDLIST *labsPO)
{
    char *cPC, nameS[NSIZE];
    int n, ok, tok;

    DB_TAB DB_PrI(">> ParseLabelRowFromBuffI line=%d corn=%d\n",line,corn);
    cPC = bufS;
    ok = TRUE;
    n = tok = 0;
    while( ISLINE(*cPC) && ok ) {
        PASS_BLANK(cPC);
        INIT_S(nameS);
        sscanf(cPC,"%s",nameS);
        if(NO_S(nameS)) {
            break;
        }
        /***
        *   If rows have labels, first token is "corner" case; Keep if cornS
        */
        if( corn && (tok==0) ) {
            if( cornS ) {
                sprintf(cornS,"%s",nameS);    
            }
        }
        else {
            DB_TAB DB_PrI("+ Adding [%d] |%s|\n", n, nameS);
            ok = AddWordlistWordI(labsPO, n++, nameS, TRUE);
        }
        PASS_WORD(cPC);
        tok++;
    }
    if(!ok) {
        PROBLINE;
        printf("Parsing table col labels at line %d |%s|\n",line,bufS);
        return(FALSE);
    }
    DB_TAB DB_PrI("<< ParseLabelRowFromBuffI %d\n",n);
    return(n);
}
/*************************************************************************
*   Loads values and optionally row lable (if labsPO) from passed string 
*   Returns the number of data columns
*/
int ParseDataRowFromBuffI(char *bufS, int line, int row, WORDLIST *labsPO, NUMLIST *valsPO)
{
    char *cPC, nameS[NSIZE];
    int n, ok, tok;
    DOUB vD;

    DB_TAB DB_PrI(">> ParseDataRowFromBuffI line=%d row=%d labs=%p\n",line,row,labsPO);
    cPC = bufS;
    ok = TRUE;
    n = tok = 0;
    while( ISLINE(*cPC) && ok ) {
        PASS_BLANK(cPC);
        /***
        *   First token is a label?
        */
        if( (labsPO) && (tok==0) ) {
            INIT_S(nameS);
            sscanf(cPC,"%s",nameS);
            DB_TAB DB_PrI("+ Adding lable [%d] |%s|\n",row,nameS);
            ok = AddWordlistWordI(labsPO, row, nameS, TRUE);
        }
        else {
            vD = BAD_R;
            sscanf(cPC,"%lf",&vD);
            if(BAD_REAL(vD)) {
                break;
            }
            DB_TAB DB_PrI("+ Adding val n=%d vD=%lf\n",n,vD);
            ok = AppendNumlistDoubI(valsPO, vD);
            n++;
        }
        PASS_WORD(cPC);
        tok++;
    }
    if(!ok) {
        PROBLINE;
        printf("Parsing table values at line %d, row %d |%s|\n",line,row+1,bufS);
        return(FALSE);
    }
    DB_TAB DB_PrI("<< ParseDataRowFromBuffI %d\n",n);
    return(n);
}
/**************************************************************************
*   Write out table contents to named file 
*/
int WriteTableFileI(TABLE *tabPO, int head, int mask, char *nameS)
{
    FILE *outPF;

    if(!(outPF=OpenUFilePF(nameS,"w",NULL))) {
        return(FALSE);
    }
    DumpTable(tabPO,head,mask,outPF);
    CHECK_NFILE(outPF,nameS);
    return(TRUE);
}
/**************************************************************************
*   Write out table contents to file outPF
*/
void DumpTable(TABLE *tabPO,int settings,int mask,FILE *outPF)
{
    int r,c;
    char nameS[NSIZE];
    DOUB vD;

    VALIDATE(tabPO,TABLE_ID);
    HAND_NFILE(outPF);
    DumpTableDescription(tabPO,mask,"# ",outPF);
    if(settings) {
        DumpTableSettings(tabPO,mask,outPF);
    }
    DumpTableColHeadings(tabPO,mask,outPF);
    for(r=0;r<tabPO->nrow;r++)
    {
        if( mask && (!tabPO->rmask[r]) ) {
            continue;
        }
        fprintf(outPF,"%s",tabPO->prefix); 
        if(tabPO->rlab) {
            GetTableRowLabI(tabPO,r,nameS,NSIZE-1);
            fprintf(outPF,tabPO->prlform,nameS); 
            fprintf(outPF,"%s",tabPO->pvsep); 
        }
        for(c=0; c < tabPO->ncol; c++)
        {
            if( mask && (!tabPO->cmask[c]) ) {
                continue;
            }
            GetTableValI(tabPO,r,c,&vD);
            fprintf(outPF,tabPO->pvform,vD); 
            if(c < (tabPO->ncol -1) ) {
                fprintf(outPF,"%s",tabPO->pvsep); 
            }
        }
        fprintf(outPF,"\n"); 
    }
}
/**************************************************************************
*   Column headings
*/
void DumpTableColHeadings(TABLE *tabPO,int mask,FILE *outPF)
{
    int c,w;
    char nameS[NSIZE],cformS[DEF_BS];

    VALIDATE(tabPO,TABLE_ID);
    if(!tabPO->clab) {
        return;
    }
    HAND_NFILE(outPF);
    /***
    *   If not CSV, make lable width same as number format
    */
    sprintf(cformS,"%%s");
    if(!IsTableCSVFormatI(tabPO)) {
        if(sscanf(tabPO->pvform,"%%%d.",&w)) {
            sprintf(cformS,"%%-%ds",w);
        }
    }
    fprintf(outPF,"%s",tabPO->prefix); 
    if(tabPO->rlab) {
        if(!NO_S(tabPO->corner)) {
            fprintf(outPF,tabPO->prlform,tabPO->corner); 
        }
        else {
            fprintf(outPF,tabPO->prlform,"RowNames"); 
        }
        fprintf(outPF,"%s",tabPO->pvsep); 
    }
    for(c=0;c<tabPO->ncol;c++)
    {
        if( mask && (!tabPO->cmask[c]) ) {
            continue;
        }
        GetTableColLabI(tabPO,c,nameS,NSIZE-1);
        fprintf(outPF,cformS,nameS); 
        fprintf(outPF,"%s",tabPO->pvsep); 
    }
    fprintf(outPF,"\n"); 
}
/**************************************************************************
*   Summary story about table
*/
void DumpTableSettings(TABLE *tabPO,int mask,FILE *outPF)
{
    VALIDATE(tabPO,TABLE_ID);
    HAND_NFILE(outPF);
    fprintf(outPF,"# Prec.  =  %d\n",tabPO->prenum);
    fprintf(outPF,"# RLabs  =  %d\n",tabPO->rlab);
    fprintf(outPF,"# CLabs  =  %d\n",tabPO->clab);
    fprintf(outPF,"# PvForm = |%s|\n",tabPO->pvform);
    fprintf(outPF,"# PrForm = |%s|\n",tabPO->prlform);
    fprintf(outPF,"# Prefix = |%s|\n",tabPO->prefix);
}
/**************************************************************************
*   Description about table
*/
void DumpTableDescription(TABLE *tabPO,int mask,char *preS,FILE *outPF)
{
    int r,c;
    char nameS[NSIZE],fnameS[NSIZE];

    VALIDATE(tabPO,TABLE_ID);
    HAND_NFILE(outPF);
    fprintf(outPF,"%s",preS); 
    fprintf(outPF,"Table at %p\n",tabPO); 
    GetTableNamesI(tabPO,nameS,fnameS,NSIZE-1);
    if(!NO_S(nameS)) {
        fprintf(outPF,"%s",preS); 
        fprintf(outPF,"Table name: %s\n",nameS);
    }
    if(!NO_S(fnameS)) {
        fprintf(outPF,"%s",preS); 
        fprintf(outPF,"Table file: %s\n",fnameS);
    }
    GetTableDimsI(tabPO,&r,&c,mask);
    fprintf(outPF,"%s",preS); 
    fprintf(outPF,"%d Rows X %d Cols\n",r,c);
}
