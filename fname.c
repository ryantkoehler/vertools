/*
* fname.c
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
#define __MAIN__
#include "prim.h"

#define VERSION_S "Fname Version 0.32"

int main(int argc, char **argv);
void FnameUse(void);
int FnameI(int argc, char **argv);

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(FnameI(argc,argv),NULL) ); }
/**************************************************************************/
void FnameUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <filename> [...options]\n");
    printf("    <infile>   Is a filename\n");
    printf("    -p         Print path\n");
    printf("    -b         Print base\n");
    printf("    -e         Print extension\n");
    printf("    -l         Length of filename component\n");
    printf("Options may be combined\n");
}
/**************************************************************************/
int FnameI(int argc, char **argv)
{
    char pathS[BBUFF_SIZE], baseS[BBUFF_SIZE], extS[BBUFF_SIZE];
    char inS[BBUFF_SIZE], outS[BBUFF_SIZE];
    int i,p,b,e,l;

    p=b=e=l=FALSE;
    if(!ParseArgsI(argc,argv,"S -p B -b B -e B -l B",
        inS,&p,&b,&e,&l,
        (int *)NULL))
    {
        FnameUse();
        return(FALSE);
    }
    /***
    *   If path and extension, need base too
    */
    if(p && e)
    {
        b = TRUE;
    }
    /***
    *   Get parts into substrings then put together the parts
    */
    GetFilePartsI(inS,pathS,baseS,extS);
/*
printf("path |%s|\n",pathS);
printf("base |%s|\n",baseS);
printf("ext  |%s|\n",extS);
*/
    INIT_S(outS);
    i = 0;
    if(p)
    {
        strcat(outS,pathS); i++;
    }
    if( (b) && (!NO_S(baseS)) )
    {
        if( (p) && (!NO_S(pathS)) )
        {   
            if( pathS[strlen(pathS)-1]!='/' ) 
                strcat(outS,"/");   
        }
        strcat(outS,baseS); i++;
    }
    if( (e) && (!NO_S(extS)) )
    {
        if( (b) && (!NO_S(baseS)) )
        {       
            strcat(outS,".");   
        }
        strcat(outS,extS); i++;
    }
    if(l)
    {
        printf("%d\n",INT(strlen(outS)));   
    }
    else if(i)
    {   
        printf("%s\n",outS);    
    }
    return(TRUE);
}
