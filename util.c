/*
* util.c
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
#include <stdarg.h>
#include "prim.h"

#define DB_ALLOC    if(DB[1])
#define DB_FREE     if(DB[2])
#define DB_INQFILE  if(DB[3])
#define DB_ENVVAR   if(DB[4])

/***
*   Flag to exit program on all done? 
*/
#define EXIT_ON_DONE    TRUE
#define ERR_TO_STDERR   FALSE

/***
*   Name of debug file to look for at start up 
*/
#define DB_FILE_S "debug"



void DBComlineCheck(int argc, char **argv);

/********************************************************************* iii
*   First call
*   Initializes globals and debug
*/
void Init(int argc, char **argv)
{
    InitGvars();
    DBInit();
    if(!argv) {
        argc = 0;
    }
    DBComlineCheck(argc,argv);
    return;
}
/*************************************************************************
*   Init globals
*/
void InitGvars(void)
{
    filecountGI = numblockGI = 0;
    outfileGPF = stdout;
    return;
}
/*************************************************************************
*   Last call
*   Reports problems if fail or any dropped memory or file pointers
*   TRUE and FALSE are legitimate-non-problem status codes
*/
int AllDoneI(int status,char *progS)
{
    int ok=0;

    if(CheckGvarsI(TRUE)) {
        ok = NOT_OK;
        fprintf(stdout,"--- ABNORMAL TERMINATION CONDITIONS -------------\n");
    }
    switch(status)
    {
        case TRUE:
            /* ok = OK; */
            ok = EXIT_SUCCESS;
            break;
        case FALSE:
            /* ok = NOT_OK; */
            ok = EXIT_FAILURE;
            break;
        default:
            /* ok = NOT_OK; */
            ok = EXIT_FAILURE;
            fprintf(stdout," TERMINATING with error status: %d\n",status);
    }
    /***
    *   Report program if given
    */
    if(progS && (ok != OK)) {
        fprintf(stdout,"\n%s  Status=%d\n\n",progS,ok);
    }
    fflush(stdout);
    /***
    *   If really all done here, exit
    */
    if(EXIT_ON_DONE) {  
        exit(ok);
    }
    return(ok);
}
/**************************************************************************
*   Check global accounting variables and complain if suspect
*/
int CheckGvarsI(int err)
{
    int i;

    i = 0;
    if(filecountGI != 0) {
        if(err) {
            PrintI(" FILE COUNT DISCREPENCY: %d open\n",filecountGI);
        }
        i++;
    }
    if(numblockGI != 0) {
        if(err) {
            PrintI(" MEMORY HAS BEEN DROPPED\n");
            PrintI("       BLOCKS ALLOCATED: %d\n",numblockGI);
        }
        i++;
    }
    return(i);
}
/**************************************************************************
*   Initializes debug structure flags
*   Reads flags (numbers) from file, if present, and sets these
*/
void DBInit()
{
    int d, n, rep, com;
    FILE *fPF;
    char bufS[LINEGRAB+1];

    DBClear();
    fPF = fopen(DB_FILE_S,"r");
    if(fPF == NULL) {   
        return; 
    }
    fclose(fPF);
    fPF = FileOpenPF(DB_FILE_S,"r",FALSE);
    n = rep = com = 0;
    while(fgets(bufS,LINEGRAB,fPF) != NULL)
    {
        if(COM_LINE(bufS)) {    
            continue;   
        }
        if(strstr(bufS,"/*")) { 
            com = TRUE; 
        }
        if(strstr(bufS,"*/")) { 
            com = FALSE;    
        }
        if(com) {   
            continue;   
        }
        if(strstr(bufS,"break")) {  
            break;  
        }
        if(strstr(bufS,"Report")) {
            rep=1;  continue;
        }
        sscanf(bufS,"%d",&d);
        if((d < 0)||(d >= NUM_DB_FLAGS)) {  
            continue;   
        }
        dbflagsGC[d] = 1;
        n++;
    }
    FILECLOSE(fPF);
    if(rep) {   
        DBReport(); 
    }
    return;
}
/**************************************************************************
*   Look at command line args to see if debug should be disabled / dumped
*/
void DBComlineCheck(int argc, char **argv)
{
    int i,c,r;

    c = r = FALSE;
    for(i=1;i<argc;i++)
    {
        if( EQSTRING(argv[i],"-noDB",5) || EQSTRING(argv[i],"-nodb",5) ) {
            c++;
        }
        if( EQSTRING(argv[i],"-dumpDB",7) || EQSTRING(argv[i],"-dumpdb",7) ) {
            r++;
        }
    }
    if(r) {
        DBReport(); 
    }
    if(c) {
        DBClear();  
    }
    return;
}
/**************************************************************************/
void DBClear()
{
    int i;

    for(i=0; i<NUM_DB_FLAGS; i++)
    {   
        dbflagsGC[i] = 0;   
    }
    return;
}
/*************************************************************************
*   Lists active debug flags
*/
void DBReport()
{
    int i, n;

    n = 0;
    DB_PrI("----- REPORTING ALL DB FLAGS ON: -----\n");
    for(i=0; i<NUM_DB_FLAGS; i++)
    {
        if(dbflagsGC[i]) {
            DB_PrI("   DB: %4d\n",i);   n++;
        }
    }
    if(n==0) {  
        DB_PrI("NONE\n");   
    }
    else {  
        DB_PrI("TOTAL: %d\n",n);    
    }
    DB_PrI("\n");
    return;
}
/*************************************************************************
*   Prints error message with the following info:
*   Function name (funS)
*   Source code line number (line)
*   Source code filename (fileS)
*   Some descriptive message (mesS)
*/
void ErrorMsg(char *funS, int line, char *fileS, char *mesS)
{
    FILE *outPF;

    if(ERR_TO_STDERR) {
        outPF = stderr;
    }
    else {
        outPF = stdout;
    }
    fflush(stdout);
    fflush(stderr);

    fprintf(outPF,
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(outPF,
    "!!!!!!!!!!!!!!!! ERROR ENCOUNTERED -Sham has snapped !!!!!!!!!!!!!!\n");
    fprintf(outPF,"Function: %s\n",funS);
    fprintf(outPF,"File:     %s\n",fileS);
    fprintf(outPF,"Line:     %d\n",line);
    fprintf(outPF,"Problem:  %s\n",mesS);
    fprintf(outPF,
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(outPF,
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fflush(outPF);

    /***
    *   End the show here
    */
    AllDoneI(BOGUS,NULL);
    return;
}
/****************************************************************************
*   Prints error message banner and sets flags to prevent redundant banners
*/ 
void ErrorBanner(int level,FILE *oPF)
{
    HAND_NFILE(oPF);
    switch(level)
    {
        case WARN_LEV:
            fprintf(oPF,
                    "#***************************** WARNING ******************************\n");
            break;
        case PROB_LEV:
            if(level > errlevelGI) {
                fprintf(oPF,"\n");
                fprintf(oPF,
                    "#***************************** PROBLEM ******************************\n");
                fprintf(oPF,"\n");
            }
            break;
        case ABORT_LEV:
            if(level > errlevelGI) {
                fprintf(oPF,"\n");
                fprintf(oPF,
                    "#***************************** ABORTING *****************************\n");
                fprintf(oPF,"\n");
            }
            break;
    }
    fflush(oPF);
    if(level > errlevelGI) {
        errlevelGI = level; 
    }
    return;
}
/**************************************************************************
*   Compares an object pointer to the type expected of it and gives an
*   error message if they are different
*/
void ValidObj(OBJPTR obPO, int is, int type, int line, char *fileS)
{
    FILE *outPF;

    if(ERR_TO_STDERR) {
        outPF = stderr;
    }
    else {
        outPF = stdout;
    }
    fflush(stdout);
    fflush(stderr);
    if(obPO == NULL) {  
        ErrorMsg("ValidObj",line,fileS,"NULL OBJ"); 
    }
    if(is != type) {
        fprintf(outPF,"ACTUAL ID:    X%x\t%d\n",is,is);
        fprintf(outPF,"EXPECTED ID:  X%x\t%d\n",type,type);
        fprintf(outPF,"ADDRESS:      %p\n",obPO);
        fflush(outPF);
        ErrorMsg("ValidObj",line,fileS,"BAD DATATYPE");
    }
    return;
}
/******************************************************************** ppp
*   Prints to outfileGPF
*/
int PrintI(char *formPS, ...)
{
    va_list args;
    int cnt;

    va_start(args,formPS);
    cnt = vfprintf(outfileGPF,formPS,args);
    va_end(args);
    return(cnt);
}
/********************************************************************
*   Prints to standard out (for debug)
*/
int DB_PrI(char *formPS, ...)
{
    va_list args;
    int cnt;

    va_start(args,formPS);
    cnt = vfprintf(stdout,formPS,args);
    va_end(args);
    fflush(stdout);
    return(cnt);
}
/********************************************************************/
void TimeStamp(char *preS, FILE *fPF)
{
    char bufS[DEF_BS];

    if(!fPF) {
        fPF = stdout;
    }
    if(preS) {
        fprintf(fPF,"%s",preS);
    }
    FillDateString(bufS);
    fprintf(fPF,"%s\n",bufS);
    return;
}
/********************************************************************
*   Attempts to open a file of passed name and mode
*   if successful, increments global active file count
*/
FILE *FileOpenPF(char *nameS, char *modeS, int error)
{
    FILE *fPF;
    
    if( (!strcmp(nameS,"-")) || (!strcmp(nameS,"stdin")) ) {
        fPF = stdin;
    }
    else {
        fPF = fopen(nameS,modeS);
    }
    if(fPF != NULL) {   
        filecountGI++;  
    }
    else {  
        if(error) {
            PROBLINE;
            PrintI("FILE OPEN FAILED: |%s|\n",nameS);   
            PrintI("            MODE: |%s|\n",modeS);   
            SEPLINE;
        }
    }
    return(fPF);    
}
/************************************************************************/
int IsFileStdinI(FILE *fPF)
{
    if(fPF) {
        if(fPF == stdin) {
            return(TRUE);
        }
    }
    return(FALSE);
}
/*************************************************************************
*   Closes the passed file and decrements the global active file count
*/
void FileClose(FILE *fPF)
{
    BOG_CHECK(fPF == NULL);
    if(! IsFileStdinI(fPF) ){
        fclose(fPF);    
    }
    filecountGI--;
    return;
}
/*****************************************************************************
*   Close file (if real) and report name
*/
void NewFileClose(FILE *fPF, char *nameS)
{
    if(fPF) {
        FileClose(fPF); 
        printf("NEW FILE: %s\n",nameS); 
    }
    return;
}
/************************************************************************
*   Tries to open file using literal name first, if fails trys to
*   open again after expanding environment variables
*   -If passed fname isn't NULL, write actual file-open name to it
*
*   Extend mode so that "quite" may be specified for failure,
*       e.g. if mode = "r q" then open with "r" and don't report failure
*/
FILE *OpenUFilePF(char *nameS, char *mS, char *fnameS)
{
    FILE *fPF;
    char bufS[NSIZE],modeS[100];
    int warn;

    DB_INQFILE DB_PrI(">> OpenUFilePF |%s| |%s|\n",nameS,mS);
    if(fnameS != NULL) {
        INIT_S(fnameS);
    }
    /***
    *   Get mode arg to really use, and any flags to interpret here
    */
    warn = TRUE;
    INIT_S(modeS); INIT_S(bufS);
    sscanf(mS,"%s %s",modeS,bufS);
    if( (bufS[0]=='Q') || (bufS[0]=='q') ) {
        DB_INQFILE DB_PrI("+ mode has Q; setting warn off\n");
        warn = FALSE;
    }
    /***
    *   Open as is and retry after expanding environment vars if fails
    */
    INIT_S(bufS);
    fPF = FileOpenPF(nameS,modeS,FALSE);
    if(fPF == NULL) {
        ExpandEnvVarI(nameS,bufS);
        fPF = FileOpenPF(bufS,modeS,FALSE);
    }
    /***
    *   Failed so... warn?
    */
    if(fPF == NULL) {
        if(warn) {
            PROBLINE;
            printf("FILE OPEN FAILED: |%s|\n",nameS);
            if(!NO_S(bufS))
                printf("   Expanded name: |%s|\n",bufS);
            printf("            MODE: |%s|\n",modeS);
            SEPLINE;
        }
    }
    /***
    *   Opened file. Place to store used name?
    */
    if(fPF && fnameS) {
        if(NO_S(bufS)) {
            strcpy(fnameS,nameS); 
        }
        else {
            strcpy(fnameS,bufS); 
        }
    }
    DB_INQFILE DB_PrI("<< OpenUFilePF %p\n",fPF);
    return(fPF);
}
/************************************************************************
*   Dereferences any environment variables in the first passed string
*   and fills the second passed string
*   Environment vars are detected by '$' in name
*/
int ExpandEnvVarI(char *inS, char *outS)
{
    int i;
    char *pS,*eS,envS[BBUFF_SIZE],tempS[BBUFF_SIZE];

    DB_ENVVAR {
        DB_PrI(">> ExpandEnvVarI\n");
        DB_PrI("+ |%s|\n",inS);
    }
    /***
    *   Save input, as output may be the same string
    */
    sprintf(tempS,"%s",inS);    
    if(!(pS = (char *)strstr(inS,"$"))) {
        DB_ENVVAR {DB_PrI("<< ExpandEnvVarI FALSE no '$' found\n");}
        sprintf(outS,"%s",inS);
        return(FALSE);
    }
    /***
    *   Starting after $, seperate environment var from rest of string
    */
    pS++;
    DB_ENVVAR DB_PrI("+ have $, looking for '/', ' ', '}'\n");
    sprintf(envS,"%s",pS);
    i = 0;
    while(ISLINE(envS[i])) {
        if(envS[i]=='/') {
            break;
        }
        i++;
    }
    envS[i] = '\0';
    DB_ENVVAR DB_PrI("+ have env var |%s|\n",envS);
    eS = getenv(envS);  
    if(eS == NULL) {
        DB_ENVVAR {DB_PrI("<< ExpandEnvVarI FALSE |%s| not found\n",envS);}
        return(FALSE);
    }
    pS = &tempS[i]; pS++;
    DB_ENVVAR DB_PrI("+ have filename |%s|\n",pS);
    sprintf(outS,"%s%s",eS,pS);
    DB_ENVVAR { 
        DB_PrI("+ expanded |%s|\n",outS);
        DB_PrI("<< ExpandEnvVarI TRUE\n");
    }
    return(TRUE);   
}
/********************************************************************* mmm
*   Allocates block of memory and increments the global memory counts
*   12/99 RTK kill size encoding, as toxic for doubles (8 byte blocks!)
*   8/05 Simplify away all fancy size accounting stuff....simple sham now!
*/
void *AllocPO(int num, int size,int line, char *fileS)
{
    void *newPI;

    DB_ALLOC DB_PrI(">> AllocPO num=%d, size=%d\n",num,size);
    if( (num<1) || (size<1) ) {
        PROBLINE;
        printf("num=%d\tsize=%d\n",num,size);
        ERR("AllocPO","Number or size less than 1");
        return(NULL);
    }
    /***
    *   Allocate size + 1 UNLG to hold size value
    */
    newPI = calloc(num,size);
    if(newPI == NULL) {
        fprintf(stderr,"Asking for memory %d of size %d = %d total\n",num,
            size,num*size);
        fflush(stderr);
        ERR("AllocPO","calloc failed");
        return(NULL);
    }
    numblockGI++;
    DB_ALLOC {
        DB_PrI("+ Called from %s %d\n",fileS,line);
        DB_PrI("<< AllocPO %p blocks %d\n",newPI,numblockGI);
    }
    return(newPI);  
}
/**************************************************************************
*   12/99 RTK kill size encoding, as toxic for doubles (8 byte blocks!)
*   8/05 Simplify away all fancy size accounting stuff....simple sham now!
*/
void *ReAllocPO(void *obPO, int num, int size, int line, char *fileS)
{
    void *newPI;

    DB_ALLOC DB_PrI(">> ReAllocPO ob=%p num=%d, size=%d\n",obPO,num,size);
    if(obPO == NULL) {
        ErrorMsg("ReAllocPO",line,fileS,"Passed NULL OBJECT\n");
        return(NULL);
    }
    if( (num<1) || (size<1) ) {
        PROBLINE;
        printf("num=%d\tsize=%d\n",num,size);
        ERR("AllocPO","Number or size less than 1");
        return(NULL);
    }
    /***
    *   Reallocte to new size + extra UNLG
    */
    newPI = realloc(obPO, num * size);
    if(newPI == NULL) {
        fprintf(stderr,"reAsking for memory %d of size %d = %d total\n",num,
            size,num*size);
        fflush(stderr);
        ERR("AllocPO","realloc failed");
        return(NULL);
    }
    DB_ALLOC {
        DB_PrI("+ Called from %s %d\n",fileS,line);
        DB_PrI("<< ReAllocPO %p blocks %d\n", newPI, numblockGI);
    }
    return(newPI);  
}
/*************************************************************************
*   Frees passed block of memory and decrements global memory counts
*   12/99 RTK kill size encoding, as toxic for doubles (8 byte blocks!)
*/
int FreeI(void *obPO, int line, char *fileS)
{
    if(obPO == NULL) {
        ErrorMsg("FreeI",line,fileS,"Passed NULL OBJECT\n");
    }
    DB_FREE DB_PrI(">> FreeI %p, size $d\n",obPO,sizeof(*obPO));
    numblockGI--;
    free(obPO);
    DB_FREE 
    {
        DB_PrI("+ Called from %s %d\n",fileS,line);
        DB_PrI("<< FreeI, blocks %d\n",numblockGI);
    }
    return(TRUE);
}
/**************************************************************************/
void ShowMemoryStatus(char *whatS)
{
    PrintI("------------- Memory staus: -------------------------\n");
    if(whatS) {
        PrintI("%s\n",whatS);
    }
    PrintI("       BLOCKS ALLOCATED: %d\n\n",numblockGI);
} 
