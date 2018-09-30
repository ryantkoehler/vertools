/*
* dfutil.c
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
#include <ctype.h>
#include <math.h>
#include "prim.h"

#define DB_DFU  if(DB[13])

/****************************************************************************
*   Count cols and rows in data table file
*/
int DataTableDimsI(FILE *inPF,int hr,int *ncolPI,int *nrowPI)
{
    char *cPC, bufS[BBUFF_SIZE+1];
    int line,ncol,nc,nrow;

    *ncolPI = *nrowPI = 0;
    line = ncol = nrow = 0;
    while(fgets(bufS,BBUFF_SIZE,inPF) != NULL)
    {
        if(line < hr) {
            line++;
            continue;
        }
        line++;
        /****
        *   Count the columns; first row is special as others compare to this
        */
        cPC = bufS;
        PASS_BLANK(cPC);
        nc = 0;
        while(isgraph(INT(*cPC)))
        {
            nc++;
            NEXT_WORD(cPC);
        }
        if(nrow == 0) {
            ncol = nc;
        }
        else if(nc != ncol) {
            printf("Missing table data on line %d\n",line);
            printf("Expecting %d values, %d found\n",ncol,nc);
            return(FALSE);
        }
        nrow++;
    }
    *ncolPI = ncol;
    *nrowPI = nrow;
    return(TRUE);
}
/*************************************************************************
*   Wrappers for FileStatsI and FileLinesI
*   ... someday nice to have functions determine if given file or name????
*/
int FilenameStatsI(char *fnameS, int *linePI, int *minPI, int *maxPI, 
    int *comPI, int *blankPI)
{
    int n;
    FILE *fPF;

    if(!(fPF=OpenUFilePF(fnameS, "r", NULL))) {
        return(-1);
    }
    n = FileStatsI(fPF, linePI, minPI, maxPI, comPI, blankPI);
    FILECLOSE(fPF);
    return(n);
}
/************************************************************************/
int FilenameLinesI(char *fnameS, int ig_com, int ig_blank)
{
    int n;
    FILE *fPF;

    if(!(fPF=OpenUFilePF(fnameS, "r", NULL))) {
        return(-1);
    }
    n = FileLinesI(fPF, ig_com, ig_blank);
    FILECLOSE(fPF);
    return(n);
}
/*************************************************************************
*   Reads through a file collecting stats on line number and size as well
*       as number of comment and blank lines
*   Returns total number of characters read in file
*/
int FileStatsI(FILE *fPF, int *linePI, int *minPI, int *maxPI, int *comPI, 
    int *blankPI)
{
    int c,n,p,g,line,min,max,com,blank;
    off_t fpos;

    if(!fPF) {   
        return(0);  
    }
    line = n = 0;
    fpos = ftell(fPF);
    max = -TOO_BIG;
    min = TOO_BIG;
    com = blank = p = g = 0;
    /***
    *   Char at a time
    */
    while( (c=fgetc(fPF)) != EOF)
    {
        n++;
        p++;
        /***
        *   Comment line
        */
        if( (p==1) && (c=='#') ) {
            com++;
        }
        /***
        *   Graphic?
        */
        if(isgraph(INT(c))) {
            g++;
        }
        /***
        *   End of line?
        */
        if(!ISLINE(c)) {
            line++;
            if(g==0) {
                blank++;
            }
            min = MIN_NUM(min,p);
            max = MAX_NUM(max,p);
            g = p = 0;
        }
    }
    /***
    *   Set values to any real pointers
    */
    if(linePI) {
        *linePI = line;
    }
    if(minPI) {
        *minPI = min;
    }
    if(maxPI) {
        *maxPI = max;
    }
    if(comPI) {
        *comPI = com;
    }
    if(blankPI) {
        *blankPI = blank;
    }
    fseek(fPF,fpos,0);
    return(n);
}
/*************************************************************************
*   Returns the number of lines in a file
*   If ig_com is true, comment lines are removed from the count
*   If ig_blank is true, blank lines are removed from the count
*/
int FileLinesI(FILE *fPF, int ig_com, int ig_blank)
{
    int n,nc,nb;

    if(!fPF) {   
        return(0); 
    }
    FileStatsI(fPF,&n,NULL,NULL,&nc,&nb);
    if(ig_com) {
        n -= nc;
    }
    if(ig_blank) {
        n -= nb;
    }
    return(n);
}
/*************************************************************************
*   Eats up to first non-print character in file (i.e. this line)
*   Returns number of chars read (i.e. line length)
*   MOVES FILE AHEAD
*/
int EatOneLineI(FILE *fPF)
{
    int n,c;

    n = 0;
    if(!fPF)
    {   return(n);  }
    n = 0;
    while( (c=fgetc(fPF)) != EOF)
    {
        if(!ISLINE(c)) {
            break;
        }
        n++;
    }
    return(n);
}
