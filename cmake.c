/*
* cmake.c
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
#define __MAIN__
#include "prim.h"

#define VERSION_S   "CMake Version 2.55"

#define DEBUG if(DB[8])

/***
*   Option values
*/
#define OPT_CF_S    "-O2"
#define WARN_CF_S   "-Wall"
#define WARNCC_CF_S "-fullwarn"
#define STD99_CF_S  "-std=c99"
#define LM_SYS_S    "-lm"   
#define LP_CF_S     "-pg"
#define LP_SYS_S    "-pg"   

#define LGD_SYS_S   "-lgd -lpng -lz -ljpeg" 

#define LXML_SYS_S  "-L/usr/local/lib -R/usr/local/lib -lxml2 -lz -lm -lsocket -lnsl"
#define LXML_INC_S  "-I/usr/local/include/libxml2"

#define LXGD_SYS_S  "-L/usr/local/lib -L/usr/local/glib/lib -R/usr/local/lib\
                    -lgdome -lglib -lxml2 -lz -lm -lsocket -lnsl"
#define LXGD_INC_S  "-I/usr/local/include\
                    -I/usr/local/include/libgdome\
                    -I/usr/local/glib/include/glib-1.2\
                    -I/usr/local/glib/lib/glib/include\
                    -I/usr/local/include/libxml2"

#define LOPT_SYS_S  "-L/opt/local/lib"
#define LOPT_INC_S  "-I/opt/local/include"

/***
*   Big file directives
xxx check on other systems?
*
#define BIG_FILE_S  "LFS_CFLAGS    := $(shell getconf LFS_CFLAGS 2>/dev/null) -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE\
LFS_LDFLAGS   := $(shell getconf LFS_LDFLAGS 2>/dev/null)\
LFS_LIBS      := $(shell getconf LFS_LIBS 2>/dev/null)\
CFLAGS        += $(LFS_CFLAGS)\
LDFLAGS       = $(LFS_LDFLAGS)\
LDLIBS        = $(LFS_LIBS)"
*/
#define BIG_CF_S    "-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64"




/***
*   Default names, buff sizes, prototypes
*/
#define LIST_FILE_S "cfile.lis"
#define MAKE_FILE_S "makefile"
#define M_F 250
#define M_S 50
#define IS_C_FILE(tc) ((tc)=='c')
#define IS_A_FILE(tc) ((tc)=='a')

int main(int argc, char **argv);
void Gen_makeUse(void);
int Gen_makefileI(int argc, char **argv);
int GetSourceListI(FILE *fPF, char files[M_F][M_S], char *typeS);
int CheckSourceListI(char fils[M_F][M_S],char *typeS,int tl,int nf);
int FindProgNameI(char filesS[M_F][M_S],char *typeS,int nf,int tl,char *nameS);
int FindHeaderFileI(char *bufS, char *hfileS);

/**************************************************************************/
int main(int argc, char **argv) 
{   Init(argc,argv); exit( AllDoneI(Gen_makefileI(argc,argv),NULL) ); }
/**************************************************************************/
void Gen_makeUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("   -i XXX  Set input file listing (1 source/line) as XXX\n");
    printf("   -o XXX  Set output makefile to name XXX\n");
    printf("   -p XXX  Set target program or library name to XXX\n");
    printf("   -lib    Make a library (rather than program)\n");
    printf("   -mak    Output \"<prog>.make\"\n");
    printf("   -nw     No checks / warnings for CFLAGS\n");
    printf("   -nm     No math library linking\n");
    printf("   -nop    No optimization flags (normally %s)\n",OPT_CF_S);
    printf("   -stc99  Add standard flags: %s\n",STD99_CF_S);
    printf("   -lp     Link profiler (-pg for CFLAGS and SYSLIBS)\n");
    printf("   -lgd    Link with graphics libs required by GD lib\n");
    printf("   -lxml   Link with xmlib2 (xml2) XML libraries\n");
    printf("   -lxgd   Link with gdome XML libraries\n");
    printf("   -lopt   Link (and include) with path /opt/local/* e.g. Mac ports)\n");
    printf("   -cc     Use 'cc' rather than 'gcc' for compile command\n");
    printf("   -nlarg  No large file directives (default is LARGE)\n");
    BARLINE;
    printf("   Note:   Default input is from \"%s\"\n",LIST_FILE_S);
    printf("           Listed C files up to a blank are processed, with the\n");
    printf("               last C file name used to name the target program\n");
    printf("           Default output is \"%s\"\n",MAKE_FILE_S);
}
/**************************************************************************
*
*/
int Gen_makefileI(int argc, char **argv)
{
    int i,files,hfiles;
    int link_math,do_nw,do_opt,do_cc,targ_lib,do_mak,do_lp,do_no,do_large;
    int do_lxml,do_lxgd,do_lgd,do_std,do_lopt;
    char bufS[DEF_BS], typeS[M_F]; 
    char sourceS[M_F][M_S]; 
    char inS[DEF_BS], nameS[DEF_BS], outS[DEF_BS], hfileS[DEF_BS];
    FILE *inPF, *outPF, *cPF;

    INIT_S(outS); INIT_S(nameS); INIT_S(inS); 
    do_cc = do_nw = targ_lib = do_mak = do_lp = do_no = do_lgd = FALSE;
    do_lxml = do_lxgd = do_std = do_lopt = FALSE;
    do_opt = link_math = do_large = TRUE;
    if(!ParseArgsI(argc,argv,
        "-p S -i S -o S -nw B -nm B -nop B -lib B -mak B -cc B -lp B\
        -lgd B -lxml B -lxgd B -stc99 B -nlarg B -lopt B",
        nameS,inS,outS,&do_nw,&link_math,&do_opt,&targ_lib,&do_mak,&do_cc,
        &do_lp, &do_lgd, &do_lxml, &do_lxgd, &do_std, &do_large, &do_lopt,
        (int *)NULL))
    {
        Gen_makeUse();
        return(FALSE);
    }
    if(NO_S(inS))
    {
        sprintf(inS,"%s",LIST_FILE_S);
    }
    if(!(inPF = OpenUFilePF(inS,"r",NULL)))
    {
        return(FALSE);
    }
    /***
    *   Story of the world
    */
    printf("\n");
    printf("%s, %s %s %s\n",VERSION_S,BD_S,__DATE__,__TIME__); 
    printf("%s\n",RTK_S);
    printf("\n");
    printf("Reading source file list from %s\n",inS);
    printf("\n");
    /***
    *   Load file list 
    */
    if(!(files = GetSourceListI(inPF,sourceS,typeS)))
    {
        printf("Problem with listed source files\n");
        printf("cmake ABORTING\n");
        FILECLOSE(inPF); return(FALSE);
    }
    FILECLOSE(inPF);
    if(!CheckSourceListI(sourceS,typeS,files,targ_lib))
    {
        printf("Problem with listed source files\n");
        printf("cmake ABORTING\n");
        return(FALSE);
    }
    if(NO_S(nameS))
    {
        if(!FindProgNameI(sourceS,typeS,files,targ_lib,nameS))
        {
            printf("Failed to find an output name from input: %s\n",inS);
            printf("cmake ABORTING\n");
            FILECLOSE(inPF); return(FALSE);
        }
    }
    /***
    *   output file
    */
    if(do_mak)
    {
        sprintf(outS,"%s.make",nameS);
    }
    else if(NO_S(outS))
    {
        sprintf(outS,"%s",MAKE_FILE_S);
    }
    if(!(outPF = FileOpenPF(outS,"w",TRUE)))
    {
        FILECLOSE(inPF); return(FALSE);
    }
    if(targ_lib)
        printf("   Target Library....... %s\n",nameS);
    else
        printf("   Target Program....... %s\n",nameS);
    if(do_cc)
        printf("   Compiler............. cc\n");
    else
        printf("   Compiler............. gcc\n");
    printf("   Total source files... %d\n\n   ",files);
    /***
    *   Header whip
    */
    fprintf(outPF,"# Makefile for %s\n",nameS);
    fprintf(outPF,"# Generated by %s\n",VERSION_S);
    fprintf(outPF,"#    %s %s %s\n",BD_S,__DATE__,__TIME__);
    fprintf(outPF,"#    %s\n",RTK_S);
    TimeStamp("# ",outPF);
    fprintf(outPF,"#\n");
    fprintf(outPF,"\n");
    fprintf(outPF,"TARGET = %s\n",nameS);
    /***
    *   Do we need math?
    */
    if(do_lgd || do_lxml || do_lxgd) {
        link_math = TRUE;
    }
    /***
    *   gcc -vs- cc
    */
    if(do_cc) {
        fprintf(outPF,"CC = cc\n");
    }
    else {
        fprintf(outPF,"CC = gcc\n");
    }
    /***
    *   CFLAGS
    */
    fprintf(outPF,"CFLAGS =");
    if(!do_nw) {
        if(do_cc) {
            fprintf(outPF," %s",WARNCC_CF_S);
        }
        else {
            fprintf(outPF," %s",WARN_CF_S);
        }
    }
    if(do_std) {
        fprintf(outPF," %s",STD99_CF_S);
    }
    if(do_opt) {
        fprintf(outPF," %s",OPT_CF_S);
    }
    if(do_lp) {
        fprintf(outPF," %s",LP_CF_S);
    }
    if(do_large) {
        fprintf(outPF," %s",BIG_CF_S);
    }
    if(do_lxml) {
        fprintf(outPF," %s",LXML_INC_S);
    }
    if(do_lxgd) {
        fprintf(outPF," %s",LXGD_INC_S);
    }
    if(do_lopt) {
        fprintf(outPF," %s",LOPT_INC_S);
    }
    fprintf(outPF,"\n");    
    /***
    *   SYSLIBS 
    */
    fprintf(outPF,"SYSLIBS ="); 
    if(link_math) {
        fprintf(outPF," %s",LM_SYS_S);
    }
    if(do_lp) {
        fprintf(outPF," %s",LP_SYS_S);
    }
    if(do_lgd) {
        fprintf(outPF," %s",LGD_SYS_S);
    }
    if(do_lxml) {
        fprintf(outPF," %s",LXML_SYS_S);
    }
    if(do_lxgd) {
        fprintf(outPF," %s",LXGD_SYS_S);
    }
    if(do_lopt) {
        fprintf(outPF," %s",LOPT_SYS_S);
    }
    fprintf(outPF,"\n");    
    /***
    *   Object collection
    */
    fprintf(outPF,"OBJS =    ");
    for(i=0;i<files;i++)
    {
        DEBUG DB_PrI("+ file %2d %s\n",i,sourceS[i]);
        /***
        *   If target is a lib and last file is .a then ignore
        */
        if( (targ_lib) && (i == (files-1)) && IS_A_FILE(typeS[i]) ) {
            DEBUG DB_PrI("+  last = .a so ignoring\n");
            continue;
        }
        GetFilePartsI(sourceS[i],NULL,nameS,NULL);
        if(IS_A_FILE(typeS[i])) {
            strcat(nameS,".a");
        }
        else {
            strcat(nameS,".o");
        }
        if( (i>0) && ((i%5)==0) ) {
            fprintf(outPF," \\\n          ");
        }
        sprintf(bufS,"%-10s",nameS);
        fprintf(outPF,"%s ",bufS);
    }
    /***
    *   Target is a library -vs- program
    */
    fprintf(outPF,"\n\n$(TARGET) : $(OBJS)\n");
    if(targ_lib) {
        fprintf(outPF,"\tar r $(TARGET) $(OBJS)\n");
    }
    else {
        fprintf(outPF,"\t$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}\n");
    }
    /***
    *   Header dependencies of each object
    */
    DEBUG DB_PrI("+ expanding source header depend\n");
    for(i=0;i<files;i++)
    {
        DEBUG DB_PrI("+  file %2d %s\n",i,sourceS[i]);
        if(IS_A_FILE(typeS[i])) {
            continue;
        }
        if(!(cPF = FileOpenPF(sourceS[i],"r",FALSE))) {
            continue;
        }
        printf(".");
        GetFilePartsI(sourceS[i],NULL,nameS,NULL);
        strcat(nameS,".o");
        sprintf(bufS,"%s: ",nameS);
        fprintf(outPF,"\n%-10s",bufS);
        /***
        *   treat source as first "header"; hfiles = 1  
        */
        sprintf(bufS,"%-10s",sourceS[i]);
        fprintf(outPF,"%s ",bufS);
        hfiles = 1;
        DEBUG DB_PrI("+  looking for headers\n");
        while (fgets(bufS,LINEGRAB,cPF) != NULL)
        {
/** picky version
            if(!EQSTRING(bufS,"#include",8))
                continue;
**/
            if( (bufS[0]!='#') || (!strstr(bufS,"#include")) ) {
                continue;
            }
            if(FindHeaderFileI(bufS,hfileS)) {
                DEBUG DB_PrI("+   header[%d] %s\n",hfiles,hfileS);
                if( (hfiles%5)==0) {
                    fprintf(outPF," \\\n          ");
                }
                sprintf(bufS,"%-10s",hfileS);
                fprintf(outPF,"%s ",bufS);
                hfiles++;
            }   
        }
        fprintf(outPF,"\n");
        FILECLOSE(cPF);
    }   
    CHECK_NFILE(outPF,outS);
    return(TRUE);
}
/***************************************************************************/
int GetSourceListI(FILE *fPF, char filesS[M_F][M_S], char *typeS)
{
    char bufS[DEF_BS], nameS[DEF_BS], exS[10];
    int i,n,com;

    DEBUG DB_PrI(">> GetSourceListI\n");
    com = FALSE;
    n = i = 0;
    while(fgets(bufS,LINEGRAB,fPF))
    {
        i++;
        if(COM_LINE(bufS)) {
            continue;
        }
        if((!com) && (strstr(bufS,"/*"))) {
            com = TRUE;
        }
        if((com) && (strstr(bufS,"*/"))) {
            com = FALSE;
        }
        if(com) {
            continue;
        }
        if(strstr(bufS,"break")) {
            break;
        }
        INIT_S(nameS);
        sscanf(bufS,"%s",nameS);
        DEBUG DB_PrI(" %3d (line %3d) |%s|\n",n,i,nameS);
        if(!isgraph(INT(*nameS))) {
            break;
        }
        GetFilePartsI(nameS,NULL,bufS,exS);
        if( (strlen(bufS)>8) || (strlen(exS)>3) ) {
            printf("WARNING: %s is incompatable with dos\n",nameS);
        }
        if(strlen(nameS) >= M_S) {
            printf("ERROR:\n");
            printf("File %s is too long (max %d)\n",nameS,M_S);
            return(FALSE);
        }
/** Don't mess with case
        Lowerize(nameS);
**/
        DEBUG DB_PrI("  |%s| ex |%c|\n",nameS,*exS);
        sprintf(filesS[n],"%s",nameS);
        typeS[n] = tolower(*exS);
        if( (!IS_C_FILE(typeS[n])) && (!IS_A_FILE(typeS[n])) ) {
            printf("ERROR:\n");
            printf("File %s bogus extension ('.c' '.a' recognized)\n",
                filesS[n]);
            return(FALSE);
        }
        n++;
    }
    DEBUG DB_PrI("<< GetSourceListI %d\n",n);
    return(n);
}
/***************************************************************************/
int CheckSourceListI(char fileS[M_F][M_S],char *typeS,int nf,int tl)
{
    int i;
    FILE *tPF;

    DEBUG DB_PrI(">> CheckSourceListI\n");
    for(i=0;i<nf;i++)
    {
        DEBUG DB_PrI("+ %2d |%s|\n",i,fileS[i]);
        if((tPF = FileOpenPF(fileS[i],"r",FALSE))) {
            FILECLOSE(tPF);     
            continue;
        }
        /***
        *   If last source file is a library, treat as target and don't check
        */
        if( (tl) && (i == (nf-1)) && IS_A_FILE(typeS[i]) ) {
            DEBUG DB_PrI("+  last = .a so ignoring\n");
            continue;
        }
        printf("ERROR\n");
        printf("Can't open: %s\n",fileS[i]);
        DEBUG DB_PrI("<< CheckSourceListI FALSE\n");
        return(FALSE);
    }
    DEBUG DB_PrI("<< CheckSourceListI TRUE\n");
    return(TRUE);   
}
/*************************************************************************
*   Look for program / library name as last legit file in list before break
*/
int FindProgNameI(char fileS[M_F][M_S],char *typeS,int nf, int tl,char *nameS)
{
    int i;

    DEBUG DB_PrI(">> FindProgName, tl = %d\n",tl);
    INIT_S(nameS);
    for(i=0;i<nf;i++)
    {
        /***
        *   If not a library target, ignore library files for names
        */
        if(IS_A_FILE(typeS[i])) {
            if(!tl) {
                continue;
            }
            else {
                sprintf(nameS,"%s",fileS[i]);
            }
        }
        else
            GetFilePartsI(fileS[i],NULL,nameS,NULL);
    }
    if(NO_S(nameS)) {
        DEBUG DB_PrI("<< FindProgName FALSE; no name\n");
        return(FALSE);  
    }
    DEBUG DB_PrI("<< FindProgName TRUE |%s|\n",nameS);
    return(TRUE);
}
/*************************************************************************/
int FindHeaderFileI(char *bufS, char *hfileS)
{
    int i,h;

    DEBUG DB_PrI(">> FindHeaderFileI\n");
    INIT_S(hfileS);
    bufS[strlen(bufS)-1] = '\0';
    DEBUG DB_PrI("+ |%s|\n",bufS);
    h = i = 0;
    while((bufS[i] != '"') && (bufS[i] != '\0'))
    {
        i++;
    }
    if(bufS[i] != '"') {
        return(FALSE);
    }
    i++;
    while((bufS[i] != '"') && (bufS[i] != '\0'))
    {
        hfileS[h++] = bufS[i++];
    }
    hfileS[h] = '\0';
    DEBUG DB_PrI("<< FindHeaderFileI TRUE |%s|\n",hfileS);
    return(TRUE);
}
