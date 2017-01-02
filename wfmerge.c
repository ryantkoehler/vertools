/*
* wfmerge.c
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
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "wfutil.h"
#include "wfmerge.h"

#define DB_WF if(DB[99])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); AllDoneI(WFMergeI(argc,argv),NULL); exit(0); }
/**************************************************************************/
void WFMergeUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Use <in> [...options]\n");
    printf("   <in>        Frequency data file *listing*\n");
    printf("   -sec XXX    Second Freq data file; treat first as freq dat\n");
    printf("   -out XXX    Set output file XXX\n");
    printf("   -range # #  Limit reported output values # to #\n");
    printf("   -sub        Subtract second from first (def=add)\n");
    printf("   -dsr        Difference sum ratio\n");
    printf("   -fdsr       Frequency weighted difference sum ratio\n");
    printf("   -sf         Sort output on frequency (high to low)\n");
    printf("   -deg        Combine counts for +/- degenerate seqs\n");
}
/**************************************************************************/
int WFMergeI(int argc,char **argv)
{
    int nlis,dgen,fsort,sub,dsr,fdsr;
    DOUB loD,hiD;
    WORDFREQ *fPO,*sPO;
    char fS[DEF_BS],sS[DEF_BS],outS[DEF_BS],bufS[DEF_BS],wordS[DEF_BS];
    FILE *inPF,*outPF;

    sub = fsort = dgen = dsr = fdsr = FALSE;
    loD = 0.0;
    hiD = TOO_BIG_D;
    INIT_S(sS);
    INIT_S(outS);
    if(!ParseArgsI(argc,argv,
        "S -sec S -out S -sub B -sf B -deg B -dsr B -fdsr B\
        -ran D2",
        fS, sS, outS, &sub, &fsort, &dgen, &dsr, &fdsr, &loD,&hiD,
        (int *)NULL))
    {
        WFMergeUse();
        return(FALSE);
    }
    fPO = sPO = NULL;
    inPF = outPF = NULL;
    nlis = 0;
    /***
    *   Consistent options
    */
    if( dsr || fdsr )
    {
        if(NO_S(sS))
        {
            PROBLINE;
            printf(" -dsr & -fdsr only work with -sec <seq>\n");
            return(FALSE);
        }
    }
    if(fsort && dgen) {
        PROBLINE;
        printf("Sort and Degen don't work together\n");
        return(FALSE);
    }
    /***
    *   Input = two real things or file list?
    */
    if(!NO_S(sS)) {
        if(!GetWordFreqsI(fS,&fPO)) {
            printf("Couldn't load word freqs from %s\n",fS);
            return(FALSE);
        }
        if(!GetWordFreqsI(sS,&sPO)) {
            CHECK_WORDFREQ(fPO);
            printf("Couldn't load word freqs from %s\n",sS);
            return(FALSE);
        }
        if( dsr || fdsr ) {
            NormalizeFrecs(fPO);
            NormalizeFrecs(sPO);
        }
        if(fdsr) {
            MergeWordfreqsI(fPO,sPO,FDSR_WF,fPO);
        }
        else if(dsr) {
            MergeWordfreqsI(fPO,sPO,DSR_WF,fPO);
        }
        else if(sub) {
            MergeWordfreqsI(fPO,sPO,SUB_WF,fPO);
        }
        else {
            MergeWordfreqsI(fPO,sPO,ADD_WF,fPO);
        }
    }
    else {
        if(!(inPF=OpenUFilePF(fS,"r",NULL))) {
            return(FALSE);
        }
        /***
        *   Get the first one on the list
        */
        INIT_S(wordS);
        while(fgets(bufS,LINEGRAB,inPF)) {
            if(COM_LINE(bufS)) {
                continue;
            }
            if(BlankStringI(bufS)) {
                continue;
            }
            sscanf(bufS,"%s",wordS);
            break;
        }
        if(!GetWordFreqsI(wordS,&fPO)) {
            FILECLOSE(inPF);
            printf("Couldn't load word freqs from %s\n",wordS);
            return(FALSE);
        }
        printf("# Have first: %s\n",fPO->name);
        nlis++;
        /***
        *   Now start partying with rest in list
        */
        while(fgets(bufS,LINEGRAB,inPF)) {
            if(COM_LINE(bufS)) {
                continue;
            }
            if(BlankStringI(bufS)) {
                continue;
            }
            INIT_S(wordS);
            sscanf(bufS,"%s",wordS);
            if(!GetWordFreqsI(wordS,&sPO)) {
                PROBLINE;
                printf("Couldn't load word freqs from %s\n",wordS);
                CHECK_FILE(inPF);
                CHECK_WORDFREQ(fPO); CHECK_WORDFREQ(sPO);
                ABORTLINE;
                return(FALSE);
            }
            if(!MergeWordfreqsI(fPO,sPO,ADD_WF,fPO)) {
                PROBLINE;
                printf("Couldn't merge words from %s\n",wordS);
                CHECK_FILE(inPF);
                CHECK_WORDFREQ(fPO); CHECK_WORDFREQ(sPO);
                ABORTLINE;
                return(FALSE);
            }
            nlis++;
            CHECK_WORDFREQ(sPO);
        }
    }
    /***
    *   Output file
    */
    outPF = NULL;
    if(!NO_S(outS)) {
        if(!(outPF=OpenUFilePF(outS,"w",NULL)))
        {
            CHECK_FILE(inPF);
            CHECK_WORDFREQ(fPO); CHECK_WORDFREQ(sPO);
            return(FALSE);
        }
    }
    HAND_NFILE(outPF);
    /***
    *   dump the whip
    */
    if(fsort) {
        SortWordFreqs(fPO,SORT_HILO);
    }
    sprintf(fPO->name,"%s",fS);
    WordfreqSummary(fPO,loD,hiD,outPF);
    DumpWordsI(fPO,dgen,loD,hiD,outPF);
    /***
    *   All done
    */
    CHECK_FILE(inPF);
    CHECK_NFILE(outPF,outS);
    CHECK_WORDFREQ(fPO); CHECK_WORDFREQ(sPO);
    return(TRUE);
}
