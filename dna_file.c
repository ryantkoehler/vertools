/*
* dna_file.c
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
#include "prim.h"
#include "dna.h"

#define DB_DNA_IO if(DB[105])


/***************************************************************************
*   Figure out sequence type given options and input name
*   Return format code or false if error
*/
int FigureSeqFileTypeI(int iraw, int iseq, int ifas, char *fnameS, int error)
{
    int form;

    if(iraw) {
        form = SEQFM_RAW;
    }
    else if(ifas) {
        form = SEQFM_FASTA;
    }
    else if(iseq) {
        form = SEQFM_SEQ;
    }
    else {
        if( STDIN_STR(fnameS) ) {
            form = SEQFM_RAW;
        }
        else {
            form = GuessSeqFileTypeI(fnameS, error);
        }
    }
    if(!OkSeqInFormatI(form,fnameS,error)) {
        return(FALSE);
    }
    return(form);
}
/***************************************************************************
*   Guess a file's type based on its extension
*   fnameS is the file name
*/
int GuessSeqFileTypeI(char *fnameS, int error)
{
    int type;
    char extS[DEF_BS];

    DB_DNA_IO DB_PrI(">> GuessSeqFileTypeI fname=|%s|\n",fnameS);
    INIT_S(extS);
    GetFilePartsI(fnameS,NULL,NULL,extS);
    Upperize(extS);
    type = ParseSeqTypeI(extS, TRUE);
    if(IS_BOG(type) && error) {
        PROBLINE;
        printf("Unrecognized sequence file format extension\n");
        printf("    File: %s\n",fnameS);
    }
    DB_DNA_IO DB_PrI("<< GuessSeqFileTypeI %d\n",type);
    return(type);
}
/***************************************************************************
*   Parse file extension and return associated filetype if recognized
*   If exact then full correct extension is needed, else just first
*       unambiguous characters of the extension
*
*   Returns BOGUS for unknown / incomplete file extensions
*/
int ParseSeqTypeI(char *extS, int exact)
{
    int type;

    DB_DNA_IO DB_PrI(">> ParseSeqTypeI |%s| exact=%d\n",extS,exact);
    type = BOGUS;
    if(!(isgraph(INT(extS[0])))) {
        DB_DNA_IO DB_PrI("<< ParseSeqTypeI %d !isgraph\n",type);
        return(type);
    }
    if(exact) {
        if(EQSTRING(extS,"DNA",3)) { 
            type = SEQFM_RAW;
        }
        else if(EQSTRING(extS,"FAS",3)) {
            type = SEQFM_FASTA;
        }
        else if(EQSTRING(extS,"SEQ",3)) {
            type = SEQFM_SEQ;
        }
    }
    else {
        switch(toupper(INT(extS[0])))
        {
            case 'D':   type = SEQFM_RAW;   break;
            case 'F':   type = SEQFM_FASTA; break;
            case 'S':   type = SEQFM_SEQ;   break;
        }
    }
    DB_DNA_IO DB_PrI("<< ParseSeqTypeI %d\n",type);
    return(type);
}
/*************************************************************************
*   Fill seq file format extension string
*/
void FillSeqFtypeExtString(int type,char *typeS)
{
    switch(type)
    {
        case SEQFM_RAW:     sprintf(typeS,"raw");   break;
        case SEQFM_SEQ:     sprintf(typeS,"seq");   break;
        case SEQFM_FASTA:   sprintf(typeS,"fas");   break;
        default:            sprintf(typeS,"unk");   break;
    }
}
/*************************************************************************
*   Fill seq file format description string
*/
void FillSeqFtypeDescString(int type,char *typeS)
{
    switch(type)
    {
        case SEQFM_RAW:     sprintf(typeS,"raw / dna"); break;
        case SEQFM_SEQ:     sprintf(typeS,"seq"); break;
        case SEQFM_FASTA:   sprintf(typeS,"fasta"); break;
        default:            sprintf(typeS,"Unknown??? (%d)",type);  break;
    }
}
/*************************************************************************
*   Is format code ok?
*/
int OkSeqInFormatI(int type,char *nameS,int error)
{
    switch(type) {
        case SEQFM_RAW:
        case SEQFM_SEQ:
        case SEQFM_FASTA:
            return(TRUE);
    }
    if(nameS) {
        if( (!NO_S(nameS)) && error ) {
            printf("Unknown format (type code: %d) for file\n",type);
            printf("  |%s| = ????\n",nameS);
        }
    }
    return(FALSE);
}
