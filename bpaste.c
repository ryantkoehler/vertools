/*
* bpaste.c
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
#include <string.h>
#include <ctype.h>
#define __MAIN__
#include "prim.h"

int main(int argc, char **argv);
void BPasteUse(void);
int BPasteI(int argc, char **argv);

#define VERSION_S "BPaste version 0.22"

#define LINE_SIZE   100000
#define DEF_DEL_S   " "

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(BPasteI(argc,argv),NULL) );}
/*******************************************************************/
void BPasteUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <ffile> <sfile> [...options]\n");
    printf("   <ffile>   First file to get lines from\n");
    printf("   <sfile>   Second file to get lines from\n");
    printf("   -out XXX  Output to file XXX\n");
    printf("   -dt       Use tab for Delimiter\n");
    printf("   -del XXX  Use XXX for Delimiter\n");
    printf("   -icom     Ignore commented lines\n");
    printf("\n");
}
/**************************************************************************/
int BPasteI(int argc, char **argv)
{
    FILE *finPF,*sinPF,*outPF;
    char finS[NSIZE], sinS[NSIZE], outS[NSIZE], delS[100];
    char fbufS[LINE_SIZE+1], sbufS[LINE_SIZE+1];
    int ok,icom,delt;

    icom = delt = FALSE;
    strcpy(delS,DEF_DEL_S);
    INIT_S(outS);   
    if(!ParseArgsI(argc, argv,
        "S S -out S -del S -dt B -icom B",
        finS, sinS, outS, delS, &delt, &icom,
        (int *)NULL))
    {
        BPasteUse();
        return(FALSE);
    }
    finPF = OpenUFilePF(finS,"r",NULL);
    sinPF = OpenUFilePF(sinS,"r",NULL);
    if( (!finPF) || (!sinPF) )
    {
        CHECK_FILE(finPF);
        CHECK_FILE(sinPF);
        return(FALSE);  
    }
    outPF = NULL;
    if(!NO_S(outS))
    {
        if(!(outPF= OpenUFilePF(outS,"w",NULL)))
        {   
            FILECLOSE(sinPF);   
            FILECLOSE(finPF);   
            return(FALSE); 
        }
    }
    HAND_NFILE(outPF);
    if(delt)
    {
        sprintf(delS,"\t");
    }
    ok = TRUE;
    while(ok)
    {
        if(!fgets(fbufS,LINE_SIZE,finPF))
        {
            break;
        }
        if(!fgets(sbufS,LINE_SIZE,sinPF))
        {
            break;
        }
        if(icom)
        {
            while(COM_LINE(fbufS))
            {
                if(!fgets(fbufS,LINE_SIZE,finPF))
                {
                    ok = FALSE;
                    break;
                }
            }
            while(COM_LINE(sbufS))
            {
                if(!fgets(sbufS,LINE_SIZE,sinPF))
                {
                    ok = FALSE;
                    break;
                }
            }
        }
        if(!ok)
        {
            break;
        }
        fbufS[strlen(fbufS)-1] = '\0';
        sbufS[strlen(sbufS)-1] = '\0';
        fputs(fbufS,outPF);
        fputs(delS,outPF);
        fputs(sbufS,outPF);
        fputs("\n",outPF);
    }
    CHECK_NFILE(outPF,outS);
    FILECLOSE(sinPF);   
    FILECLOSE(finPF);   
    return(TRUE);
}
