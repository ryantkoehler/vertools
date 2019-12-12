/*
* shuffle.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
*
* The programs and source code of the vertools collection are free software.
* They are distributed in the hope that they will be useful,
* WITHOUT ANY WARRANTY OF FITNESS FOR ANY PARTICULAR PURPOSE.  
* 
* Permission is granted for research, educational, and possibly commercial use 
*   and modification as long as 1) Code and any derived works are not 
*   redistributed for any fee, and 2) Proper credit is given to the authors. 
*   If you wish to include this software in a product, or use it commercially,
*   please contact the authors.
*
* See https://www.verdascend.com/ for more
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define __MAIN__
#include "prim.h"

#define VERSION_S "Shuffle Version 0.23"

typedef struct LINEREC
{
    int id;
    off_t fpos;
}LINEREC;

#define BLOCK_SIZE 5000

int main(int argc, char **argv);
void ShuffleUse();
int ShuffleI(int argc, char **argv);
int qsort_lrec(const void *e1, const void *e2);


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(ShuffleI(argc,argv),NULL) );}
/*******************************************************************/
void ShuffleUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> [...options]\n");
    printf("   <infile>  A text file with lines to be shuffled\n");
    printf("   -out XXX  Output to file XXX\n");
    printf("   -seed #   Set seed for random subset qualification to #\n");
    printf("\n");
}
/**************************************************************************/
int ShuffleI(int argc, char **argv)
{
    int i,n,nline,seed;
    off_t fpos;
    LINEREC *linesPO;
    FILE *inPF,*outPF;
    char inS[NSIZE],outS[NSIZE],bufS[BBUFF_SIZE+1];

    INIT_S(outS);   
    seed = BAD_I;
    if(!ParseArgsI(argc, argv,
        "S -out S -seed I",
        inS,outS,&seed,
        (int *)NULL))
    {
        ShuffleUse();
        return(FALSE);
    }
    if(!(inPF= OpenUFilePF(inS,"r",NULL)))
    {   return(FALSE);  }
    nline = BLOCK_SIZE;
    if(!(linesPO = (LINEREC *) ALLOC(nline,sizeof(LINEREC)))) {
        printf("Failed to allocate for %d lines\n",nline);
        FILECLOSE(inPF);
        return(FALSE);
    }
    outPF = NULL;
    if(!NO_S(outS)) {
        if(!(outPF= OpenUFilePF(outS,"w",NULL))) {  
            FILECLOSE(inPF); CHECK_FREE(linesPO);
            return(FALSE);  
        }
    }
    HAND_NFILE(outPF);
    Srand(seed);
    /***
    *   Load lines
    */
    n = 0;
    fpos = ftell(inPF);
    while(fgets(bufS,BBUFF_SIZE,inPF) != NULL)
    {
        if(COM_LINE(bufS)) {
            fpos = ftell(inPF);
            continue;
        }
        if(n >= nline) {
            nline += BLOCK_SIZE;
            linesPO = (LINEREC *) REALLOC(linesPO,nline,sizeof(LINEREC));
            if(!linesPO)
            {
                printf("Failed to allocate for %d lines\n",nline);
                FILECLOSE(inPF);
                return(FALSE);
            }
        }
        linesPO[n++].fpos = fpos;
        fpos = ftell(inPF);
    }
    /***
    *   Assign random numbers to eash line and sort on these
    */
    for(i=0;i<n;i++) {
        linesPO[i].id = RandI(n*10);
    }
    qsort(linesPO,n,sizeof(LINEREC),qsort_lrec);
    /***
    *   Rewind and spit out
    */
    rewind(inPF);
    for(i=0;i<n;i++) {
        fseek(inPF,linesPO[i].fpos,0);
        if( !fgets(bufS,BBUFF_SIZE,inPF)){ BOG_CHECK(TRUE); }
        fputs(bufS,outPF);
    }
    CHECK_FREE(linesPO);
    FILECLOSE(inPF);
    CHECK_NFILE(outPF,outS);
    return(TRUE);
}
/****************************************************************************/
int qsort_lrec(const void *e1, const void *e2)
{
    if( INT( ((LINEREC *)e1)->id ) > INT( ((LINEREC *)e2)->id) )
        return(1);
    if( INT( ((LINEREC *)e1)->id ) < INT( ((LINEREC *)e2)->id) )
        return(-1);
    return(0);
}
