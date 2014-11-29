/*
* ccase.c
*
* Copyright 2014 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <ctype.h>
#define __MAIN__
#include "prim.h"

#define VERSION_S "CCase Version 0.2"

int main(int argc, char **argv);
void CCaseUse(void);
int CCaseI(int argc, char **argv);

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(CCaseI(argc,argv),NULL) ); }
/**************************************************************************/
void CCaseUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
	printf("Usage: <input> [...options]\n");
	printf("  -u   Change to Upper case\n");
	printf("  -l   Change to Lower case\n");
	printf("  -f   Treat input as a file and lowerize contents\n");
}
/**************************************************************************/
int CCaseI(int argc, char **argv)
{
	int file,up,low,c;
	char inS[DEF_BS*2];
	FILE *fPF;	

	file = up = low = FALSE;
    if(!ParseArgsI(argc,argv,"S -f B -u B -l B",
		inS, &file, &up, &low,
		(int *)NULL))
    {
        CCaseUse();
        return(FALSE);
    }
/* printf("inS |%s|\n",inS); */
	if(file)
	{
		if(!(fPF=OpenUFilePF(inS,"r",NULL)))
		{
        	return(FALSE);
		}
		while((c = fgetc(fPF)) != EOF)
		{
			if(up)
			{
				c = toupper(c);
			}
			else if(low)
			{
				c = tolower(c);
			}
			fputc(c,stdout);
		}
		FILECLOSE(fPF);
	}
	else
	{
		if(low)
		{
			Lowerize(inS);
		}
		else if(up)
		{
			Upperize(inS);
		}
		printf("%s\n",inS);
	}
	return(TRUE);
}
