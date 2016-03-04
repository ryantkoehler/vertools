/*
* dna_in.c
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
#include "dna.h"

#define DB_DNA_IO   if(DB[105])     

/***************************************************************************
*   Guess format based on extension and attempt to load sequence set
*   If error then complain if problems, else keep silent
*/
int GuessAndGetSeqsetI(char *fnameS, SEQSET **ssPPO,int clean,int error)
{
    int type,got;
    SEQSET *ssPO;

    DB_DNA_IO DB_PrI(">> GuessAndGetSeqsetI fname=|%s|\n",fnameS);
    *ssPPO=NULL;
    type = GuessSeqFileTypeI(fnameS,error);
    if(IS_BOG(type)) { 
        if(error) {
            PROBLINE;
            printf("Unrecognized file extension; don't know format type\n"); 
        }
        DB_DNA_IO DB_PrI("<< GuessAndGetSeqsetI FALSE\n");
        return(FALSE); 
    }
    got = ReadInSeqsetI(fnameS, type, clean, &ssPO, error);
    if(!got) { 
        DB_DNA_IO DB_PrI("<< GuessAndGetSeqsetI FALSE\n");
        return(FALSE); 
    }
    *ssPPO = ssPO;
    DB_DNA_IO DB_PrI("<< GuessAndGetSeqsetI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Attempt to read in sequence set from file fnameS with format type
*   Returns TRUE is successful and sets SEQSET pointer to new object`
*
*   If error flag is TRUE, then complain about problems, else keep silent
*/
int ReadInSeqsetI(char *fnameS, int type, int clean, SEQSET **ssPPO, int error)
{
    SEQSET *ssPO;
    FILE *fPF;

    DB_DNA_IO DB_PrI(">> ReadInSeqsetI type=%d name=|%s|\n",type,fnameS);
    *ssPPO = NULL;
    if(!(fPF=OpenUFilePF(fnameS, "r", NULL))) {
        DB_DNA_IO DB_PrI("<< ReadInSeqsetI FILE FAIL = FALSE\n");
        return(FALSE);
    }
    ssPO = GetSeqsetPO(fPF, type, clean, error);
    FILECLOSE(fPF);
    if(ssPO == NULL) {
        DB_DNA_IO DB_PrI("<< ReadInSeqsetI failed to get seqset = FALSE\n");
        return(FALSE);
    }
    SetSeqsetSource(ssPO, fnameS);
    GetFilePartsI(fnameS, NULL, ssPO->name, NULL);
    *ssPPO = ssPO;
    DB_DNA_IO DB_PrI("<< ReadInSeqsetI %p TRUE\n",ssPO);
    return(TRUE);
}
/****************************************************************************
*   Returns new SEQSET datatype or NULL on failure
*/
SEQSET *GetSeqsetPO(FILE *fPF,int type,int clean,int error)
{
    int n,ok;
    SEQSET *ssPO;
    SEQ *seqPO;

    DB_DNA_IO DB_PrI(">> GetSeqsetPO type=%d\n",type);
    /***
    *   Allocate empty set object
    */
    if(!(ssPO=CreateSeqsetPO(0))) {
        printf("Problem allocating sequence set object\n");
        return(NULL);
    }
    DB_DNA_IO DB_PrI("+ have SEQSET shell allocated\n");
    /***
    *   For each seq, allocate, parse, and keep
    */
    n = 0;
    while(TRUE)
    {
        if(!(seqPO=CreateSeqPO(0, NULL, NULL))) {
            printf("Problem allocating sequence to load\n");
            CHECK_SEQSET(ssPO);
            return(NULL);
        }
        /***
        *   Parse sequence
        */
        ok = ParseSeqI(fPF, type, n+1, clean, error, seqPO);
        if(ok!=TRUE) {
            CHECK_SEQ(seqPO);
            break;
        }
        DB_DNA_IO DB_PrI("+ ok read seq %d\n",n);
        /***
        *   Add to collection
        */
        if(!AddSeqToSeqsetI(seqPO, ssPO)) {
            CHECK_SEQ(seqPO);
            CHECK_SEQSET(ssPO);
            return(NULL);
        }
        n++;
    }
    DB_DNA_IO DB_PrI("+ n=%d\n",n);
    /***
    *   If nothing actually read in, kill the holder set and return NULL
    */
    if (n < 1) {
        CHECK_SEQSET(ssPO);
        DB_DNA_IO DB_PrI("<< GetSeqsetPO NULL\n");
        return(NULL);
    }
    DB_DNA_IO DB_PrI("<< GetSeqsetPO %p\n",ssPO);
    return(ssPO);
}
/****************************************************************************
*   Attempts to parse sequence of format iform from file
*
*   Loaded ok = TRUE
*   End of file (done) = FALSE
*   Errors = BOGUS
*/
int ParseSeqI(FILE *inPF, int iform, int num, int clean, int error, SEQ *seqPO)
{
    int ok;
    char nameS[NSIZE];

    DB_DNA_IO DB_PrI(">> ParseSeqI iform=%d, num=%d\n",iform,num);
    VALIDATE(seqPO,SEQ_ID);
    InitSeq(seqPO,FALSE,FALSE);
    /***
    *   Get name first
    */
    switch(iform)
    {
        case SEQFM_RAW:     
            ok = ParseOneLineSeqI(inPF, TRUE, error, seqPO);
            break;
        case SEQFM_SEQ:     
            ok = ParseOneLineSeqI(inPF, FALSE, error, seqPO);
            sprintf(nameS, SEQ_NAME_S, num);
            SetSeqName(seqPO, nameS);
            break;
        case SEQFM_FASTA:   
            ok = ParseFastaSeqI(inPF,error,seqPO);
            break;
        default:
            printf("Bogus input format code=%d\n",iform);
            ERR("ParseSeqI","bogus format code");
            return(FALSE);
    }
    if(ok != TRUE)
    {
        DB_DNA_IO DB_PrI("<< ParseSeqI PROB_seq %d\n",ok);
        return(ok);
    }
    ok = FinishSeqSettingsI(seqPO,clean,error);
    DB_DNA_IO DB_PrI("<< ParseSeqI GOOD %d\n",ok);
    return(ok);
}
/*************************************************************************
*   Single line sequence input.
*       If has_name is TRUE, then line should include <name> <seq>
*       Else, only <sequence>
*   Returns TRUE if got sequence
*   Returns FALSE if end of file
*   Returns BOGUS if problem
*/
int ParseOneLineSeqI(FILE *fPF, int has_name, int error, SEQ *seqPO)
{
    int c,t,n,s;
    char nameS[NSIZE];
    char seqS[BBUFF_SIZE+1];

    VALIDATE(seqPO,SEQ_ID);
    InitSeq(seqPO,FALSE,FALSE);
    INIT_S(seqS);
    /***
    *   Depending on expected name, set name and spacer-token to true at start
    */
    s = 0;
    if( has_name ) {
        INIT_S(nameS);
        n = 0;
    }
    else {
        sprintf(nameS, "%s", "no-name");
        n = 1;
    }
    t = n;
    while((c = fgetc(fPF)) != EOF) {
        /***
        *   End of line. If we've got any seq, done
        */
        if( !ISLINE(c) ) {
            if(s) {
                break;
            }
            continue;
        }
        /***
        *   '#' = comment, so ignore rest of line
        */
        if( c=='#' ) {
            EatOneLineI(fPF);
            continue;
        }
        /***
        *   Space; Next token, but ignore
        */
        if( !isgraph(INT(c)) ) {
            if(n) {
                t++;
            }
            continue;
        }
        /***
        *   Name = first token; If first char is '#' it's a comment line
        */
        if( (!t) && has_name ) {
            if(n < NSIZE) {
                nameS[n++] = c;
            }
        }
        /***    
        *   Sequence chars
        */
        else {
            /***
            *   Don't collect numbers
            */
            if( isdigit(INT(c)) ) {
                continue;
            }
            /***
            *   If we've filled up the temp buffer, add this to seq and reset
            */
            if(s >= BBUFF_SIZE) {
                if(!AppendSeqSequenceI(seqPO,seqS,s)) {
                    return(BOGUS);
                }
                s = 0;
            }
            seqS[s++] = c;
            seqS[s] = '\0';
        }
    }
    /***    
    *   Nothing = end of file = FALSE
    */
    if( !s ) {
        return(FALSE);
    }
    /***
    *   Finish name (if has one) & sequence
    */
    if(has_name) {
        nameS[n] = '\0';
    }
    if(s) {
        if(!AppendSeqSequenceI(seqPO, seqS, s)) {
            return(BOGUS);
        }
    }
    /***    
    *   No sequence = error; The append above sets length
    */
    if(seqPO->len < 1) {
        if(error) {
            PROBLINE;
            printf("PROBLEM parsing RAW format line. Starts with |%s|\n",nameS);
        }
        return(BOGUS);
    }
    SetSeqName(seqPO, nameS);
    return(TRUE);
}
/*************************************************************************
*   Fasta format = line starting with > then lines until next > or end
*   Returns TRUE if got sequence
*   Returns FALSE if end of file
*   Returns BOGUS if problem
*/
int ParseFastaSeqI(FILE *fPF,int error,SEQ *seqPO)
{
    int c,t,n,s;
    char nameS[NSIZE];
    char seqS[BBUFF_SIZE+1];

    VALIDATE(seqPO,SEQ_ID);
    InitSeq(seqPO,FALSE,FALSE);
    t = n = s = 0;
    /***
    *   Find head line, starts with '>' in first position
    */
    while((c = fgetc(fPF)) != EOF) {
        if(c == '>') {
            t++;
            break;
        }
        EatOneLineI(fPF);
    }
    /***
    *   Nothing to start with so we're done
    */
    if(!t) {
        return(FALSE);
    }
    /***
    *   Collect name
    */
    while((c = fgetc(fPF)) != EOF) {
        if(!ISLINE(c)) {
            break;
        }
        if(n < NSIZE) {
            nameS[n++] = c;
        }
    }
    /***
    *   Now collect sequence until end or next record
    */
    t = 0;
    while((c = fgetc(fPF)) != EOF) {
        if( (c == '>') && (t == 0) ) {
            ungetc(c,fPF);
            break;
        }
        if(!ISLINE(c)) {
            t = 0;
            continue;
        }
        /***
        *   Don't collect spaces or numbers
        */
        if( (!isgraph(INT(c))) || (isdigit(INT(c))) ) {
            continue;
        }
        if( (t==0) && (c=='#') ) {
            EatOneLineI(fPF);
            continue;
        }
        t++;
        if(s >= BBUFF_SIZE) {
            if(!AppendSeqSequenceI(seqPO,seqS,s)) {
                return(BOGUS);
            }
            s = 0;
        }
        seqS[s++] = c;
        seqS[s] = '\0';
    }
    /***
    *   Finish name & sequence
    */
    if(n) {
        nameS[n] = '\0';
    }
    else {
        strcpy(nameS,DEF_SEQ_NAME_S);
    }
    if(s) {
        if(!AppendSeqSequenceI(seqPO,seqS,s)) {
            return(BOGUS);
        }
    }
    /***    
    *   No sequence = error; The append above sets length
    */
    TrimSeqTrailingCharsI(seqPO, TRUE);
    if(seqPO->len < 1) {
        if(error) {
            PROBLINE;
            printf("PROBLEM parsing FASTA format. Starts with |%s|\n",nameS);
        }
        return(BOGUS);
    }
    SetSeqName(seqPO, nameS);
    return(TRUE);
}
/*************************************************************************
*   Trims any chars from end of sequence
*   Default = non-alphabetic, as might be at end of fasta record (e.g. "//")
*/
int TrimSeqTrailingCharsI(SEQ *seqPO, int what)
{
    int i,n;

    VALIDATE(seqPO,SEQ_ID);
    i = seqPO->len -1;
    n = 0;
    while(i>=0) {
        if( isalpha(INT(seqPO->seq[i])) ) {
            break;
        }
        seqPO->seq[i--] = '\0';
        seqPO->len--;
        n++;
    }
    return(n);
}
