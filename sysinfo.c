/*
* sysinfo.c
*
* Copyright 2015 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <errno.h>
#include "prim.h"

#define DB_SYS  if(DB[6])


/**************************************************************************
*   Standardized version splash report
*/
void VersionSplash(FILE *outPF, char *verS, char *preS, int bars)
{
    char timeS[DEF_BS],pfxS[DEF_BS];
    
    HAND_NFILE(outPF);
    INIT_S(pfxS);
    if(preS) {  
        sprintf(pfxS,"%s",preS);    
    }
    if(bars) {  
        fprintf(outPF, LINEBAR_S); 
    }
    fprintf(outPF,"%s%s, %s %s %s\n",pfxS,verS,BD_S,__DATE__,__TIME__);
    fprintf(outPF,"%s%s\n",pfxS,RTK_S);
    PrintSysInfo(outPF,pfxS);
    FillDateString(timeS);
    fprintf(outPF,"%sRun date %s\n",pfxS,timeS);
    if(bars) {  
        fprintf(outPF, LINEBAR_S); 
    }
}
/*****************************************************************************
*   Fill passed strings with sysinfo content
*   If passed strings are non-null, assumes big enough to hold answer
*/
void GetSystemInfo(char *userS, char *hostS, char *osS, char *verS, 
    char *archS)
{
    int stat;
    char *cPC;
    struct utsname utsO;
    
    DB_SYS DB_PrI(">> GetSystemInfo\n");
    /***
    *   Initalize passed strings
    */
    if(userS)   INIT_S(userS);
    if(hostS)   INIT_S(hostS);
    if(osS)     INIT_S(osS);
    if(verS)    INIT_S(verS);
    if(archS)   INIT_S(archS);
    /***
    *   User name
    */
    if(userS) {
        cPC = getlogin(); 
        if(cPC) {
            strcpy(userS,cPC);
        }
    }
    /***
    *   Uname; non-neg means success
    */
    stat = uname(&utsO);
    if( (hostS) && (stat>=0) ) {
        strcpy(hostS,utsO.nodename);
    }
    if( (osS) && (stat>=0) ) {
        strcpy(osS,utsO.sysname);
    }
    if( (verS) && (stat>=0) ) {
        strcpy(verS,utsO.release);
        strcat(verS,"; ");
        strcat(verS,utsO.version);
    }
    if( (archS) && (stat>=0) ) {
        strcpy(archS,utsO.machine);
    }
    DB_SYS DB_PrI("<< GetSystemInfo\n");
}
/**************************************************************************
*   Prints (subset of) systeminfo to output file (or stdout on NULL)
*   If preS is passed, this is printed as a prefix
*/
void PrintSysInfo(FILE *outPF,char *preS)
{
    char userS[NSIZE],hostS[NSIZE],osS[NSIZE],archS[NSIZE];

    HAND_NFILE(outPF);
    GetSystemInfo(userS, hostS, osS, NULL, archS);
    if(preS) {
        fprintf(outPF,"%s",preS);
    }
    fprintf(outPF,"Running on %s (%s %s) by %s\n",hostS,osS,archS,userS);
}
/***************************************************************************
*   Dereference a symbolic link
*   If dereferenced, return TRUE
*   If not a symbolic link, return FALSE
*   If error, return BOGUS
*/
int DeRefSymLinkI(char *lnS, char *drS, int max)
{
    int ok;

    DB_SYS DB_PrI(">> DeRefSymLinkI |%s|\n",lnS);
    ok = readlink(lnS,drS,max);
    if(ok<0) {
        DB_SYS DB_PrI("+ ok = %d, errno = %d \n",ok,errno);
        if(errno == EINVAL) {
            DB_SYS DB_PrI("<< DeRefSymLinkI EINVAL (not link) FALSE\n",errno);
            return(FALSE);
        }
        DB_SYS DB_PrI("<< DeRefSymLinkI BOGUS\n");
        return(BOGUS);
    }
    drS[ok] = '\0';
    DB_SYS DB_PrI("+ |%s| ok=%d\n",drS,ok);
    DB_SYS DB_PrI("<< DeRefSymLinkI TRUE\n");
    return(TRUE);
}
