/*
* filter.c
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
#include <ctype.h>
#define __MAIN__
#include "prim.h"
#include "wordlist.h"
#include "filter.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(Filter_numsI(argc,argv),NULL) );}
/*******************************************************************/
void Filter_numsUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>  File to be filtered\n");
    printf("   -out XXX  Output to file XXX\n");
    printf("   -col #    Take values from column # (def = 1)\n");
    printf("   -sc #     Skip (ignore) characters up to position # / line\n");
    printf("   -rg # #   Qualify line if value is in range from # to #\n");
    printf("   -gt #     Qualify line if value is greater than (or equal) #\n");
    printf("   -lt #     Qualify line if value is less than (or equal) #\n");
    printf("   -lrg # #  Qualify line number range # to #     [NOTE: Only data lines count]\n");
    printf("   -wlis XXX Qualify line with words (first token) listed in XXX\n");
    printf("   -kc       Keep case for token comparison (default ignore)\n");
    printf("   -wst      Word start only needs to match line (not full token)\n");
    printf("   -wsub     Word substring only needs to match line (not full token)\n");
    printf("   -rann #   Qualify random number of lines #     [NOTE: Exact, NO stdin]\n");
    printf("   -ranf #   Qualify random fraction (0 - 1) #    [NOTE: Exact, NO stdin]\n");
    printf("   -ranp #   Qualify random probability (0 -1) #  [NOTE: Approximate]\n");
    printf("   -seed #   Set random seed to #\n");
    printf("   -A #      Report # lines After qualifiying lines (like grep -A)\n");
    printf("   -B #      Report # lines Before qualifiying lines [NOTE: NO stdin]\n");
    printf("   -icbn     Ignore chars before numbers (i.e. strip leading chars)\n");
    printf("   -not      Invert line qualification test(s)\n");
    printf("   -flag     Preceed lines with 1/0 for good/bad\n");
    printf("   -pln      Preceed lines with line number\n");
    printf("   -stat     Report only stats about values\n");
    printf("   -quiet    No summary report\n");
    printf("\n");
    printf("NOTE: Filtered file lines limited to %d chars wide\n",FILTBUF_SIZE);
    printf("\n");
}
/**************************************************************************/
int Filter_numsI(int argc, char **argv)
{
    char bufS[FILTBUF_SIZE+1], *cPC;
    int ok,nok,line,pout,extra;
    FILTER *filtPO;
    
    filtPO = CreateFilterPO();
    if(!ParseArgsI(argc, argv,
        "S -not B -stat B -rg D2 -gt D -lt D -col I -out S\
        -seed I -ranf D -ranp D\
        -rann I -sc I -flag B -icbn B\
        -lrg I2 -wlis S -kc B -wst B -wsub B\
        -pln B -A I -B I -qu B",
        filtPO->inname, &filtPO->do_not, &filtPO->do_stat, 
        &filtPO->min,&filtPO->max, &filtPO->min, &filtPO->max, 
        &filtPO->col, filtPO->outname, 
        &filtPO->seed, &filtPO->ranf, &filtPO->ranp, 
        &filtPO->rann, &filtPO->skipc, &filtPO->do_flag, 
        &filtPO->do_icbn,
        &filtPO->firstl,&filtPO->lastl, filtPO->wlisname, 
        &filtPO->do_kc, &filtPO->do_wst, &filtPO->do_wsub, 
        &filtPO->do_pln, &filtPO->do_A, &filtPO->do_B, 
        &filtPO->do_quiet, 
        (int *)NULL))
    {
        Filter_numsUse();
        CHECK_FILTER(filtPO);
        return(FALSE);
    }
    /***
    *   Check options and set things up
    */
    if(!OpenFilterFilesI(filtPO)) {
        CHECK_FILTER(filtPO);
        return(FALSE);
    }
    if(!CheckFilterOptionsI(filtPO)) {
        CHECK_FILTER(filtPO);
        return(FALSE);
    }
    /***
    *   Process lines 
    */
    line = nok = extra = pout = 0;
    while(fgets(bufS,FILTBUF_SIZE,filtPO->in) != NULL) 
    {
        if( SkipThisLineI(filtPO, bufS) ) {
            continue;
        }
        line++;
        cPC = GetLineStartPC(filtPO, bufS);
        ok = IsFiltLineOkI(filtPO, line, cPC, TRUE); 
        /*   Count ok, and set for any -After extra  */
        if(ok) {
            nok++;
            if(filtPO->do_A > 0) {
                extra = filtPO->do_A + 1;
            }
        }
        /*  Not only stats = output */
        if( !filtPO->do_stat) {
            pout = FALSE;
            if( (filtPO->do_flag) || (ok) || (extra>0) ) {
                pout++;
            }
            if(pout) {
                if(filtPO->do_pln) {
                    fprintf(filtPO->out,"%d\t",line);
                }
                if(filtPO->do_flag) {
                    fprintf(filtPO->out,"%d\t",ok);
                }
                fputs(bufS,filtPO->out);
            }
        }
        extra--;
    }
    ReportFilterStats(filtPO, line, nok); 
    CHECK_FILTER(filtPO);
    return(TRUE);
}
/*************************************************************************/
int SkipThisLineI(FILTER *filtPO, char *bufS)
{
    /*   comment, blank? */
    if( COM_LINE(bufS) || BlankStringI(bufS) ) {        
        return(TRUE);   
    }
    /*  Not a data line; Nothing to check? */
    if( ! GetLineStartPC(filtPO, bufS) ) {
        return(TRUE);
    }
    return(FALSE);
}
/*************************************************************************/
char *GetLineStartPC(FILTER *filtPO, char *bufS)
{
    int j;
    char *cPC;

    cPC = bufS;
    for(j=0;j<filtPO->skipc;j++)
    {
        if(!ISLINE(*cPC)) {
            break;
        }
        cPC++;
    }
    return(cPC);
}
/*************************************************************************/
int SetBeforeMaskingI(FILTER *filtPO) 
{
    char bufS[FILTBUF_SIZE+1], *cPC;
    int line, keep, rand;

    /***
    *   Screen all lines, setting mask for qualified ones
    */
    rand = DoingRandFilterI(filtPO, TRUE);
    line = 0;
    while(fgets(bufS,FILTBUF_SIZE,filtPO->in) != NULL) 
    {
        if( SkipThisLineI(filtPO, bufS) ) {
            continue;
        }
        line ++;
        cPC = GetLineStartPC(filtPO, bufS);
        if ( IsFiltLineOkI(filtPO, line, cPC, rand) ) {
            filtPO->mask[line-1] = 1;
        }
    }
    rewind(filtPO->in);
    /***
    *   Mark all mask lines before qualified ones
    */
    keep = 0;
    while(line > 0) 
    {
        if(filtPO->mask[line-1]) {
            keep = filtPO->do_B + 1;
        }
        if(keep > 0) {
            filtPO->mask[line-1] = 1;
        }
        keep--;
        line--;
    }
    return(TRUE);
}
/*************************************************************************/
FILTER *CreateFilterPO()
{
    FILTER *fpPO;

    if(!(fpPO=(FILTER *)ALLOC(1,sizeof(FILTER)))) {
        return(NULL);
    }
    fpPO->ID = FILTER_ID;
    InitFilter(fpPO);
    return(fpPO);
}
/*************************************************************************/
int DestroyFilterI(FILTER *fpPO)
{
    VALIDATE(fpPO,FILTER_ID);
    CHECK_FILE(fpPO->in);
    CHECK_NFILE(fpPO->out,fpPO->outname);
    CHECK_FREE(fpPO->mask);
    CHECK_WORDLIST(fpPO->wlis);
    FREE(fpPO);
    return(TRUE);
}
/*************************************************************************
*   Init structure
*/
void InitFilter(FILTER *fpPO)
{
    VALIDATE(fpPO,FILTER_ID);
    INIT_S(fpPO->inname);
    fpPO->in = NULL;
    fpPO->mask = NULL;
    fpPO->n_mask = 0;
    INIT_S(fpPO->outname);
    fpPO->out = NULL;
    INIT_S(fpPO->wlisname);
    fpPO->wlis = NULL;
    fpPO->col = DEF_COL;
    fpPO->min = -TOO_BIG_D;
    fpPO->max = TOO_BIG_D;
    fpPO->firstl = fpPO->lastl = BAD_I;
    fpPO->seed = BAD_I;
    fpPO->rann = -1;
    fpPO->ranf = fpPO->ranp = -1.0;
    fpPO->do_not = FALSE;
    fpPO->do_stat = FALSE;
    fpPO->do_kc = fpPO->do_wst = fpPO->do_wsub = FALSE;
    fpPO->do_icbn = FALSE;
    fpPO->skipc = 0;
    fpPO->do_flag = fpPO->do_pln = FALSE;
    fpPO->do_A = fpPO->do_B = 0;
    fpPO->do_quiet = FALSE;
}
/*************************************************************************/
int CheckFilterOptionsI(FILTER *filtPO)
{
    Srand(filtPO->seed);
    /*** 
    *   Allocate and set mask if random exact
    */
    if(!HandleFilterMaskingI(filtPO)) {
        return(FALSE);
    }
    /*  Set before masking */
    if(filtPO->do_B > 0) {
        SetBeforeMaskingI(filtPO);
    }
    if(!HandleFilterListSetsI(filtPO)) {
        return(FALSE);
    }
    return(TRUE);
}
/*************************************************************************
*   Check options and allocate / set mask if needed
*/
int HandleFilterMaskingI(FILTER *filtPO)
{
    DOUB fracD;

    if( NeedFilterMaskingI(filtPO)) {
        if( IsFileStdinI(filtPO->in) ) {
            PROBLINE;
            printf("Options not compatible with stdin\n\n");
            return(FALSE);
        }
        /***
        *   Last two args mean ignore comments and blanks
        *   Then rewind file and allocate space for lines
        */
        filtPO->n_mask = FileLinesI(filtPO->in,TRUE,TRUE);
        rewind(filtPO->in);
        filtPO->mask = (char *)ALLOC(filtPO->n_mask,sizeof(char));
        if(!filtPO->mask) {
            PROBLINE;
            printf("Can't allocate mask for %d input lines\n",filtPO->n_mask);
            return(FALSE);
        }
        /***
        *   If doing random fraction, calc frac and set
        */
        fracD = -1.0;
        if( filtPO->rann > 0) {
            fracD = RNUM(filtPO->rann) / RNUM(filtPO->n_mask);
        }
        else if( filtPO->ranf >= 0.0) {
            fracD = filtPO->ranf;
        }
        if(fracD >= 0.0) {
            LIMIT_NUM(fracD, 0.0, 1.0);
            MaskRandSubsetI(filtPO->mask,filtPO->n_mask,RNUM(fracD));
        }
    }
    return(TRUE);
}
/************************************************************************/
int DoingRandFilterI(FILTER *filtPO, int do_ranp)
{
    int rand;

    rand = FALSE;
    rand = ( filtPO->rann > 0 ) ? TRUE : rand;
    rand = ( filtPO->ranf > 0.0 ) ? TRUE : rand;
    if(do_ranp) {
        rand = ( filtPO->ranp > 0.0 ) ? TRUE : rand;
    }
    return(rand);
}
/************************************************************************/
int NeedFilterMaskingI(FILTER *filtPO)
{
    int need;

    /***
    *   Need If random number or fraction (exact; not probability)
    */
    need = DoingRandFilterI(filtPO, FALSE);
    /*  Need if reporting do_Before lines */
    if( filtPO->do_B > 0 ) {
        /* Random probability but do_Before too; Set fractional prob */
        if( filtPO->ranp >= 0.0 ) {
            filtPO->ranf = filtPO->ranp;
        }
        need = TRUE;
    }
    return(need);
}
/*************************************************************************
*   Loading words from files?
*/
int HandleFilterListSetsI(FILTER *filtPO)
{
    if(!NO_S(filtPO->wlisname)) {
        if( ! (filtPO->wlis = CreateWordlistPO(filtPO->wlisname,NSIZE))) {
            PROBLINE;
            printf("Failed to get tokens from %s\n",filtPO->wlisname);
            return(FALSE);
        }
    }
    return(TRUE);
}
/*************************************************************************
*   Open in / out files
*/
int OpenFilterFilesI(FILTER *filtPO)
{
    if(!(filtPO->in = OpenUFilePF(filtPO->inname,"r",NULL))) {  
        return(FALSE);  
    }
    if(!NO_S(filtPO->outname)) {
        if(!(filtPO->out = OpenUFilePF(filtPO->outname,"w",NULL)))
        {   
            return(FALSE); 
        }
    }
    HAND_NFILE(filtPO->out);
    return(TRUE);
}
/**************************************************************************/
void ReportFilterStats(FILTER *filtPO, int num, int nok) 
{
    if(!filtPO->do_quiet) {
        printf("#  %d data lines\n",num);
        printf("#  %d qualified data lines (%5.2f%%)\n",
            nok, PERCENT_R(nok,num));   
    }
}
/*************************************************************************
*   Get number value from word and set it to pointer
*/
int FiltGetWordNumValI(FILTER *filtPO, char *wordS, DOUB *valPD)
{
    int w;
    char *wPC;
    DOUB rD;

    /***
    *   Ignore chars before digits?
    */
    w = 0;
    if(filtPO->do_icbn) {
        while( (isgraph(INT(wordS[w]))) && (!isdigit(INT(wordS[w]))) ) {
            w++;
        }
        if( (w>0) && (wordS[w-1]=='-') ) {
            w--;
        }
    }
    wPC = &wordS[w];
    /***
    *   Get number 
    */
    rD = BAD_D;
    sscanf(wPC,"%lf",&rD);
    if(BAD_DOUB(rD)) {
        return(FALSE);
    }
    *valPD = rD;
    return(TRUE);
}
/*************************************************************************
*   Is the line pointed to by cPC OK based on filter criteria?
*/
int IsFiltLineOkI(FILTER *filtPO, int line, char *cPC, int use_mask) 
{
    char wordS[NSIZE];
    DOUB rD;
    int ok;

    /***
    *   Just filtering on line number
    */
    ok = TRUE;
    if( !BAD_INT(filtPO->firstl) ) {
        if( (line < filtPO->firstl) || (line > filtPO->lastl) ) {
            ok = FALSE;
        }
    }
    /*  Mask?  */
    else if( (filtPO->mask) && (use_mask) ) {
        /***
        *   If mask is smaller than number of lines, disqualify
        */
        if(line > filtPO->n_mask) {
            ok = FALSE;
        }
        else {
            ok = INT(filtPO->mask[line-1]);
        }
    }
    /*  Random probability? */
    else if(filtPO->ranp >= 0.0) {
        rD = RandD(1.0);
        ok = (rD < filtPO->ranp);
    }
    /***
    *   Get column entry then check 
    *   If can't get word or number from column, simply set not ok
    */
    else {
        ok = GetNthWordI(cPC,filtPO->col,wordS);
        /***    
        *   Get number from word (if not token or number list)
        */
        if( ok && (!filtPO->wlis) ) {
            ok = FiltGetWordNumValI(filtPO, wordS, &rD);
            /***
            *   Too big or small?
            */
            if( ok && ( (rD < filtPO->min) || (rD > filtPO->max)) ) {   
                ok = FALSE; 
            }
        }
        /***
        *   In number or token list?
        */
        if( ok && (filtPO->wlis) ) {
            ok = WordInWordlistI(filtPO->wlis, wordS, filtPO->do_kc, filtPO->do_wst, 
                                filtPO->do_wsub, NULL);
        }
    }
    /* Inverting qualificaiton? */
    if(filtPO->do_not) {
        ok = !ok;
    }
    return(ok);
}
