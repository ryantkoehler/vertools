/*
* chardump.c
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
#include <ctype.h>
#define __MAIN__
#include "prim.h"

#define VERSION_S   "CharDump version 0.41"

#define CCONT_SIZE  260

/***
*   Character type codes
*/
#define PRINT_C     0
#define PRINT_S     1
#define NP_NEWL     2
#define NP_TAB      3
#define NP_WIERD    10

int main(int argc, char **argv);
void Char_dumpUse();
int Char_dumpI(int argc, char **argv);

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(Char_dumpI(argc,argv),NULL) ); }
/**************************************************************************/
void Char_dumpUse()
{
    VersionSplash(NULL,VERSION_S,"#   ",TRUE);
    printf("Usage: <file> ['-' for stdin] [...options]\n");
    printf("    <file>   is any input file to check\n");
    printf("    -da      Dump all character\n");
    printf("    -dn      Dump non-print character info\n");
    printf("    -dw      Dump weird non-print character info\n");
    printf("    -det     Report details of non-print occurance\n");
    printf("    -tw      Terminate on a weird non-print character\n");
    printf("    -out XXX Write out new file named XXX\n");
    printf("    -clean   Clean output by stripping out wierd chars\n");
    printf("    -clnl    Clean by replacing wierd chars with `new line'\n");
    printf("    -clsp    Clean by replacing wierd chars with space chars\n");
    printf("    -rp X    Replace weird chars with X\n");
    printf("    -rs      Replace weird chars Singly (one / block)\n");
    printf("    -stat    Stats only, no output printed\n");
    printf("NOTE: Characters are printed in OCTAL\n");
    printf("      Non-weird, non-print = space, newline, tab\n");
    printf("\n");
}
/**************************************************************************/
int Char_dumpI(int argc, char **argv)
{
    int i,ok,num,got,ctype,nope,ch,line,ccountPI[CCONT_SIZE];
    int dall,dnp,dwierd,twierd,dump,to255,stat,detail,clean,clnl,clsp;
    int do_rs, prev_rs;
    char inS[NSIZE], outS[NSIZE], rpS[DEF_BS];
    char c,pc,fc;
    FILE *inPF, *outPF;

    dall = dnp = dwierd = stat = clean = clnl = clsp = FALSE;
    twierd = detail = do_rs = FALSE;
    INIT_S(outS);   
    INIT_S(rpS);    
    if(!ParseArgsI(argc,argv,
        "S -da B -dn B -dw B -tw B -out S -stat B -det B -clean B -clnl B \
        -clsp B -rp S -rs B",
        inS, &dall, &dnp, &dwierd, &twierd, outS, &stat, &detail, &clean,
        &clnl, &clsp, rpS, &do_rs,
        (int *)NULL))
    {
        Char_dumpUse();
        return(FALSE);
    }
    inPF= OpenUFilePF(inS,"rb",NULL);
    if(inPF == NULL) {
        return(FALSE);
    }
    outPF = NULL;
    if(!NO_S(outS)) {
        outPF= OpenUFilePF(outS,"w",NULL);
        if(outPF == NULL) {     
            FILECLOSE(inPF);    
            return(FALSE);
        }
    }
    HAND_NFILE(outPF);
    /***
    *   Options
    */
    if(dall) {
        dnp = dwierd = TRUE;
    }
    if(dnp) {
        dwierd = TRUE;
    }
    if( dall || dnp || dwierd ) {
        dump = FALSE;
    }
    else {
        dump = TRUE;
    }
    if(detail) {
        for(i=0;i<256;i++) {
            ccountPI[i]=0;
        }
    }
    /***
    *   For each character read in
    */
    i = nope = line = ch = to255 = 0;
    got = pc = 0;
    prev_rs = FALSE;
    while( (c=fgetc(inPF)) != EOF ) {
        if(got==0) {
            pc = c;
        }
        got++;
        /***
        *   Classify this char
        */
        ctype = -1;
        if( !isprint(INT(c)) ) {
            if(c == '\n') {
                ctype = NP_NEWL;
            }
            else if(c == '\t') {
                ctype = NP_TAB;
            }
            else {
                ctype = NP_WIERD;
            }
        }
        else if(ctype<0) {
            if(isgraph(INT(c))) {
                ctype = PRINT_C;
            }
            else {
                ctype = PRINT_S;
            }
        }
        ch++;
        /***
        *    case analysis on char and (maybe) dump out
        */
        ok = TRUE;
        switch(ctype)
        {
            case PRINT_C:
            case PRINT_S:
                if((dall)&&(!stat)) {
                    printf("%8d  |%c| = %o\n",got,c,c);
                }
                if((dall)&&(detail)) {
                    ccountPI[INT(c)] += 1;
                }
                break;
            case NP_NEWL:
                line++;
                ch = 0;
                if((dnp)&&(!stat)) {
                    printf("%8d  %o = newline (%d)\n",got,c,line);
                }
                if(detail) {
                    ccountPI[INT(c)] += 1;
                }
                break;
            case NP_TAB:
                if((dnp)&&(!stat)) {
                    printf("%8d  %o = tab (line %d, char %d)\n",
                         got,c,line+1,ch);
                }
                if(detail) {
                    ccountPI[INT(c)] += 1;
                }
                break;
            case NP_WIERD:
                nope++;
                if( (dwierd) && (!stat) ) {
                    printf("%8d  %o = nonprint (line %d, char %d)",
                        got,c,line+1,ch);
                    if(isprint(INT(pc))) {
                        printf(" pc=|%c|",pc);
                    }
                    fc = fgetc(inPF);
                    if(isprint(INT(fc))) {
                        printf(" fc=|%c|",fc);
                    }
                    printf("\n");
                    ungetc(fc,inPF);
                }
                if(detail) {
                    ccountPI[INT(c)] += 1;
                }
                if(clean) {
                    ok = FALSE;
                }
                if(clnl) {
                    c = '\n';
                }
                if(clsp) {
                    c = ' ';
                }
                if(!NO_S(rpS)) {
                    c = rpS[0];
                    if(do_rs && prev_rs) {
                        ok = FALSE;
                    }
                    prev_rs++;
                }
                break;
        }
        if( (ctype == NP_WIERD) && twierd ) 
        {
            printf("Wierd character termination; char %d = %o\n",i,c);
            break;
        }
        if( (dump) && (ok) && (!stat) ) {
            fputc(c,outPF);
        }
        /***
        *   Have replace char and current char isn't weird
        */
        if( (!NO_S(rpS)) && (ctype != NP_WIERD) ) {
            prev_rs = FALSE;
        }
        pc = c;
    }
    if(stat) {
        printf("# ---- Summary ---------------------------------\n");
        printf("# Characters:  %10d  (%d kb)\n",got,ROUND(got/1000));
        printf("# Lines:       %10d\n",line);
        printf("# Wierd:       %10d\n",nope);
        if(to255>0) {
            printf("# Over255:     %10d\n",nope);
        }
        printf("# NOTE: Characters are printed in OCTAL\n");
        printf("#      Non-weird, non-print = space, newline, tab\n");
        if(detail) {
            printf("# Details: (char)  (count)\n");
            num = 0;
            for(i=0;i<256;i++) {
                if(ccountPI[i]) {
                    num++;
                    printf("#          %o %8d\n",i,ccountPI[i]);
                }
            }
            if(num == 0) {
                printf("# NONE\n");
            }
        }
    }
    CHECK_FILE(inPF);
    CHECK_NFILE(outPF,outS);
    return(TRUE);
}
