/*
* ven_salt.c
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
#include <math.h>
#include "prim.h"
#include "venpipe.h"


#define DB_VLIB if(DB[140])


/****************************************************************************
*   Takes vienna parameter file inS (assumed salt = 1M) and generates
*   an on-the-fly new one that has salt/temperature adjusted energy
*   parameters. The new filename is returned via newS.
*/
int SaltCorrectViennaParsI(char *inS, DOUB saltD, char *newS)
{
    int change;
    char bufS[DEF_BS],wordS[DEF_BS],*cPC;
    FILE *inPF,*outPF;
    DOUB oD,nD;

    DB_VLIB DB_PrI(">> SaltCorrectViennaParsI salt=%f\n",saltD);
    /***
    *   Open input 
    */
    if(!(inPF=OpenUFilePF(inS,"r",NULL))) {
        DB_VLIB DB_PrI("<< SaltCorrectViennaParsI FALSE\n");
        return(FALSE);
    }
    /***
    *   Open output
    */
    if(!OpenSaltCorOutputFileI(inS,newS,&outPF)) {
        PROBLINE;
        printf("Couldn't open temp salt-corrected vienna parameter file\n");
        printf("    Last try: |%s|\n",newS);
        DB_VLIB DB_PrI("<< SaltCorrectViennaParsI FALSE\n");
        return(FALSE);
    }
    /***
    *   First line = vienna header, followed by correction story
    */
    if(!fgets(bufS,LINEGRAB,inPF)) {
        DB_VLIB DB_PrI("<< SaltCorrectViennaParsI FALSE header line\n");
        return(FALSE);
    }
    fputs(bufS,outPF);
    fprintf(outPF,"\n");
    fprintf(outPF,"/* ALT CORRECTED VIENNA PARAMETER FILE */\n");
    fprintf(outPF,"/* ILENAME       %s */\n",newS);
    fprintf(outPF,"/* NPUT          %s */\n",inS);
    fprintf(outPF,"/* FFECTIVE SALT %4.4e */\n",saltD);
    fprintf(outPF,"\n");
    /***
    *   Party through the file, modifying parts as needed
    */
    change = 0;
    while(fgets(bufS,LINEGRAB,inPF)) {
        /***
        *   If blank just echo it
        */
        if(BlankStringI(bufS)) {
            fputs(bufS,outPF);
            continue;
        }
        /***
        *   Block id or just comment
        */
        if(bufS[0] == '#') {
            INIT_S(wordS);
            sscanf(bufS,"# %s",wordS);
            if(!NO_S(wordS))
            {
                change = ParseViennaBlockTokenI(wordS);
                if(IS_BOG(change))
                {
                    PROBLINE;
                    printf("Unrecognized vienna parameter token:\n");
                    fputs(bufS,stdout);
                    CHECK_FILE(inPF);
                    CHECK_FILE(outPF);
                    return(FALSE);
                }
            }
            fputs(bufS,outPF);
/**
            fprintf(outPF,"* Salt correction changes = %d *\n",change);
**/
            continue;
        }
        /***
        *   No change 
        */
        if(!change) {
            fputs(bufS,outPF);
            continue;
        }
        /***
        *   Now get and modify each number
        */
        cPC = bufS;
        PASS_BLANK(cPC);
        while(ISLINE(*cPC)) {
            INIT_S(wordS);
            sscanf(cPC,"%s",wordS);
            if(NO_S(wordS)) {
                break;
            }
            /***
            *   If not a number, dump token as is
            */
            oD = BAD_R;
            sscanf(wordS,"%lf",&oD);
            if( (!strncasecmp(wordS,"INF",3)) || (BAD_REAL(oD)) ) {
                fprintf(outPF," %5s",wordS);
            }
            else {
                nD = SaltCorrectVenNumberD(oD,change,saltD);
                fprintf(outPF," %5.0f",nD);
            }
            NEXT_WORD(cPC);
        }   
        fprintf(outPF,"\n");
    }
    /***
    *   Close 
    */
    CHECK_FILE(inPF);
    CHECK_FILE(outPF);
    DB_VLIB DB_PrI("<< SaltCorrectViennaParsI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Open output file to write (temp) salt-corrected parameters
*/
int OpenSaltCorOutputFileI(char *inS, char *newS, FILE **outPPF)
{
    int n,ok;
    char baseS[NSIZE];
    FILE *outPF;

    DB_VLIB DB_PrI(">> OpenSaltCorOutputFileI |%s|\n",inS);
    INIT_S(newS);
    outPF = NULL;
    /***
    *   Start with default name
    */
    GetFilePartsI(inS,NULL,baseS,NULL);
    strcat(baseS,"_saltcor");
    /***
    *   Loop until file doesn't already exist (can't open to read)
    */
    n = ok = 0;
    while(n<SALTCOR_MAXFILE) {
        sprintf(newS,"%s%d.par",baseS,n);
        outPF = FileOpenPF(newS,"r",FALSE);
        DB_VLIB DB_PrI("+ %d |%s| %p ",n,newS,outPF);
        if(!outPF) {
            DB_VLIB DB_PrI("Doesn't exist\n",n,newS);
            ok++;
            break;
        }
        DB_VLIB DB_PrI("exists\n",n,newS);
        FILECLOSE(outPF);
        n++;
    }
    /***
    *   If !ok, fail
    */
    if(!ok) {
        DB_VLIB DB_PrI("<< OpenSaltCorOutputFileI FALSE\n");
        return(FALSE);
    }
    /***
    *   Now try to open for writing
    */
    if(!(outPF=OpenUFilePF(newS,"w",NULL))) {
        DB_VLIB DB_PrI("<< SaltCorrectViennaParsI FALSE\n");
        return(FALSE);
    }
    *outPPF = outPF;
    DB_VLIB DB_PrI("<< OpenSaltCorOutputFileI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Parses a token from vienna parameter file
*   Returns:
*       BOGUS if unrecognized
*       0 if numbers are not to be modified
*       1 if numbers are to be modified at half case (dangles)
*       2 if numbers are to be modified 
*/
int ParseViennaBlockTokenI(char *wordS)
{
    int stat;

    stat = BOGUS;
    if(!strcmp(wordS,"stack_energies"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"stack_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_hairpin"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_interior"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"mismatch_multi"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"dangle5"))
    {   stat = 1;   }
    else if(!strcmp(wordS,"dangle3"))
    {   stat = 1;   }
    else if(!strcmp(wordS,"dangle5_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"dangle3_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"hairpin"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"bulge"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"internal_loop"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"ML_params"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"NINIO"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"Tetraloops"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"Triloops"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"END"))
    {   stat = 0;   }
    /***
    *   RNA stuff
    */
    else if(!strcmp(wordS,"int11_energies"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"int11_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"int12_energies"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"int12_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"int21_energies"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"int21_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"int22_energies"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"int22_enthalpies"))
    {   stat = 0;   }
    /***
    *   New in parmeter format v2.0; 
    *   Set these analgous to above (RTK 1/7/13), with enthalpies = 0 and the others
    *   assumed to be salt-dependent
    */
    else if(!strcmp(wordS,"Misc"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"Hexaloops"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"ML_params"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"bulge_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"hairpin_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"int11"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"int21"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"int22"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"interior"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"interior_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_interior_1n"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"mismatch_interior_1n_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_interior_23"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"mismatch_interior_23_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_interior_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_exterior"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"mismatch_exterior_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_multi"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"mismatch_multi_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"stack"))
    {   stat = 2;   }
    else if(!strcmp(wordS,"stack_enthalpies"))
    {   stat = 0;   }
    else if(!strcmp(wordS,"mismatch_hairpin_enthalpies"))
    {   stat = 0;   }
    return(stat);
}
/****************************************************************************
*   Salt corrected dG
*   As per Nic Peyret
*/
DOUB SaltCorrectVenNumberD(DOUB oD,int type,DOUB saltD)
{
    DOUB newD;

    if(oD==0.0) {
        return(0.0);
    }
    newD = oD;
    oD /= 100.0;
    switch(type)
    {
        /***
        *   Dangle end case, half the correction term
        */
        case 1:
            newD = oD - 0.057 * log(saltD);
            break;
        /***
        *   Normal case correction term
        */
        case 2:
            newD = oD - 0.114 * log(saltD);
            break;
    }
    newD *= 100.0;
    return(newD);
}
