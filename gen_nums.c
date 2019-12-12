/*
* gen_nums.c
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
#define __MAIN__
#include "prim.h"
#include "gen_nums.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); AllDoneI(Gen_numsI(argc,argv),NULL); return(TRUE); }
/************************************************************************/
void Gen_numsUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: -num # [options...]\n");
    printf("   -num #     Is the number of items (lines) to output\n");
    printf("   -st #      Start at number # (default = 1)\n");
    printf("   -end #     End value # (For step increment)\n");
    printf("   -inc #     Increment by # (default = 1)\n");
    printf("   -max #     After greater than #, reset to start\n");
    printf("   -min #     After less than #, reset to start\n");
    printf("   -let XY    Cycle letters from X to Y; e.g. AZ or az or Az\n");
    printf("   -out XXX   Set output to XXX\n");
    printf("   -pdp #     Print # decimal places\n");
    printf("   -pfmt # #  Print with format # wide # precision (%%#.#f)\n");
    printf("   -zp        Zero pad output\n");
    printf("   -npl #     Numbers per line (default = 1)\n");
    printf("   -tab       Separate multiple items per line with tabs\n");
    printf("   -cs XXX    Print character string XXX rather than number\n");
    printf("   -pre XXX   Output as prefix to lines in file XXX\n");
    printf("   -suf XXX   Output as suffix to lines in file XXX\n");
    printf("   -ran # #   Random numbers in range # to #\n");
    printf("   -gaus # #  Gaussian distribution centered at # with S.D. #\n");
    printf("   -bits #    Bits with # percent on / line (e.g. 33.3 = 1/3)\n");
    printf("   -seed #    Set random number seed\n");
}
/**************************************************************************
*
*/
int Gen_numsI(int argc, char **argv)
{
    char bufS[DEF_BS],outS[NSIZE],stS[DEF_BS], *maskPC;
    char preS[NSIZE], sufS[NSIZE], lineS[BBUFF_SIZE];
    char letS[DEF_BS], fC,lC,cC;
    int i,j,num,seed;
    int p_wid,p_pre,do_zp,do_tab,numpl,p_pdp;
    DOUB vD,startD,endD,incD,mrandD,randD,drandD,gmeanD,gsdD,maxD,minD;
    REAL bitsR;
    FILE *outPF, *inPF;

    INIT_S(stS); INIT_S(outS);
    INIT_S(preS); INIT_S(sufS);
    do_zp = do_tab = FALSE;
    numpl = 1;
    num = BOGUS;
    mrandD = randD = gmeanD = gsdD = BAD_R;
    bitsR = BAD_R;
    seed = BAD_I;
    p_wid = DEF_FORM_W;
    p_pre = DEF_FORM_P;
    p_pdp = BAD_I;
    startD = 1.0;
    endD = BAD_R;
    incD = 1.0;
    minD = -TOO_BIG_D;
    maxD = TOO_BIG_D;
    INIT_S(letS);
    if(!ParseArgsI(argc,argv,
        "-num I -st D -zp B -cs S -seed I -ran D2 -inc D -out S -pfmt I2\
        -gaus D2 -bits R -min D -max D -npl I -tab B -en D -pdp I -let S\
        -pre S -suf S",
        &num, &startD, &do_zp, stS, &seed, &mrandD,&randD, &incD, outS,
        &p_wid,&p_pre, &gmeanD,&gsdD, &bitsR, &minD, &maxD, &numpl,
        &do_tab, &endD, &p_pdp, &letS, preS, sufS,
        (int *)NULL))
    {
        Gen_numsUse();
        return(FALSE);
    }
    if(IS_BOG(num)) {
        Gen_numsUse();
        return(FALSE);
    }
    /***
    *   If prefix / suffix, num = number of lines
    */
    inPF = NULL;
    if(!NO_S(preS)) {
        if(!(inPF=OpenUFilePF(preS, "r", NULL))) {
            return(FALSE);
        }
    }
    else if(!NO_S(sufS)) {
        if(!(inPF=OpenUFilePF(sufS, "r", NULL))) {
            return(FALSE);
        }
    }
    /***
    *   Check if can do anything before starting
    */
    if( (num < 1) && (!inPF) ) {
        /*
        PROBLINE;
        printf(" Can't enumerate %d items!\n",num);     
        ABORTLINE;
        return(FALSE);
        */
        CHECK_FILE(inPF);
        return(TRUE);
    }
    /***
    *   Initialize
    */
    outPF = NULL;
    maskPC = NULL;
    Srand(seed);
    fC = lC = cC = '+';
    /***
    *   Check options, set up, etc
    */
    if(numpl<1) {
        PROBLINE;
        printf("Number of items per line too small: %d\n",numpl);
        CHECK_FILE(inPF);
        return(FALSE);
    }
    if( (!BAD_REAL(endD)) && (num > 1) ) {
        incD = (endD - startD) / RNUM(num - 1); 
    }
    if(!BAD_INT(p_pdp)) {
        p_pre = p_pdp;
    }
    if(!BAD_REAL(bitsR)) {
        if( (bitsR<0.0) || (bitsR>100.0) ) {
            PROBLINE;
            printf("Bit fraction out of bounds (0-100): %f\n",bitsR);
            ABORTLINE;
            CHECK_FILE(inPF);
            return(FALSE);
        }
        if(! (maskPC= (char *)ALLOC(numpl,sizeof(char)) ) ) {
            PROBLINE;
            printf("Failed to allocate for bit mask\n");
            ABORTLINE;
            CHECK_FILE(inPF);
            return(FALSE);
        }
        INIT_S(stS);
    }
    else if(!BAD_REAL(gmeanD)) {
        if(gsdD < TINY_R) {
            PROBLINE;
            printf("SD too small: %f\n",gsdD);
            ABORTLINE;
            CHECK_FILE(inPF);
            return(FALSE);
        }
        INIT_S(stS);
    }
    else if(!BAD_REAL(mrandD)) {
        drandD = randD - mrandD;
        if(drandD < TINY_R) {
            PROBLINE;
            printf(" Random range %f to %f is too narrow\n",mrandD,randD);      
            ABORTLINE;
            CHECK_FILE(inPF);
            return(FALSE);
        }
        INIT_S(stS);
    }
    else if(!NO_S(letS)) {
        if(!ParseLetterArgI(letS, &fC, &lC)) {
            PROBLINE;
            printf(" Didn't parse letters XY from |%s|\n",letS);        
            ABORTLINE;
            CHECK_FILE(inPF);
            return(FALSE);
        }
        cC = fC;
    }
    /***
    *   Output file?
    */
    if(!NO_S(outS)) {
        if(!(outPF = OpenUFilePF(outS,"w",NULL))) {
            ABORTLINE;
            CHECK_FILE(inPF);
            return(FALSE);
        }
    }
    HAND_NFILE(outPF);
    /***
    *   Format string
    */
    if(!NO_S(stS)) {
        sprintf(bufS,"%%s");
    }
    else if(!NO_S(letS)) {
        sprintf(bufS,"%%c");
    }
    else if(do_zp) {
        sprintf(bufS,"%%0%d.%df",p_wid,p_pre);
    }
    else {
        sprintf(bufS,"%%%d.%df",p_wid,p_pre);
    }
    /***
    *   Party until done 
    */
    vD = startD;    
    i = 0;
    while(TRUE)
    {
        /***
        *   Input line or next number; Break if done
        */
        if(inPF) {
            if(! fgets(lineS,LINEGRAB,inPF)) {
                break;
            }
            ReplaceChars('\n',lineS,'\0',lineS);
            /* If number suffix, print line here */
            if(!NO_S(sufS)) {
                fprintf(outPF,"%s\t",lineS);
            }
        }
        else {
            i++;
            if(i>num) {
                break;
            }
        }
        /* Bit mask sham */
        if(maskPC) {
            MaskRandSubsetI(maskPC,numpl,bitsR/100.0);
        }
        /***
        *   For each item / line
        */
        for(j=0;j<numpl;j++)
        {
            if(j>0) {
                if(do_tab) {
                    fprintf(outPF,"\t");
                }
                else {
                    fprintf(outPF," ");
                }
            }
            if(!NO_S(stS)) {
                fprintf(outPF,bufS,stS);
            }
            else if(!NO_S(letS)) {
                fprintf(outPF,bufS,cC);
                GetNextLetterI(cC, fC, lC, &cC);
            }
            else {
                if(maskPC) {
                    vD = DNUM(maskPC[j]); 
                }
                else if(!BAD_REAL(gmeanD)) {
                    vD = RandGaussD(gmeanD,gsdD); 
                }
                else if(!BAD_REAL(randD)) {     
                    vD = RandR(randD) + mrandD; 
                }
                fprintf(outPF,bufS,vD);
            }
            /***
            *   Default behavior is to increment and check
            */
            vD += incD; 
            if( (vD>maxD) || (vD<minD) ) {
                vD = startD;
            }
        }
        /* If number is prefix, print line here */
        if(!NO_S(preS)) {
            fprintf(outPF,"\t%s",lineS);
        }
        fprintf(outPF,"\n");
    }
    /**
    *   All done
    */
    CHECK_FREE(maskPC);
    CHECK_FILE(inPF);
    CHECK_NFILE(outPF,outS);
    return(TRUE);
}
/***************************************************************************/
int ParseLetterArgI(char *letS, char *fPC, char *lPC)
{
    if(strlen(letS) < 2) {
        PROBLINE;        
        printf("Letter argument too short\n");
        return(FALSE);
    }
    if( (!isalpha(INT(letS[0]))) || (!isalpha(INT(letS[1]))) ) {
        printf("Letter argument not alphabetic\n");
        return(FALSE);
    }
    if(letS[0] >= letS[1]) {
        printf("Letter args not asscending: %c > %c\n", letS[0], letS[1]);
        return(FALSE);
    }
    *fPC = letS[0];
    *lPC = letS[1];
    return(TRUE);
}
/***************************************************************************/
int GetNextLetterI(char cC, char fC, char lC, char *nPC)
{
    while(cC <= lC) 
    {
        cC++;
        if(isalpha(INT(cC))) {
            break;
        }
    }
    if(cC > lC) {
        cC = fC;
    }
    *nPC = cC;
    return(TRUE);
}
