/*
* filter.c
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
    printf("   -rg # #   Qualify line if value is in range from # to # inclusive\n");
    printf("   -gt #     Qualify line if value is greater than or equal to #\n");
    printf("   -lt #     Qualify line if value is less than or equal to #\n");
    printf("   -vex      Value exclusive filtering; e.g. make -rg, -lt, -gt NOT inclusive\n");
    printf("   -abs      Use absolute value\n");
    printf("   -lrg # #  Qualify line number range # to #     [NOTE: Only data lines count]\n");
    printf("   -brg # #  Qualify block number range # to #    [NOTE: Only data lines count]\n");
    printf("   -wlis XXX Qualify line with words (first token) listed in XXX\n");
    printf("   -kc       Keep case for token comparison (default ignore)\n");
    printf("   -wst      Word start only needs to match line (not full token)\n");
    printf("   -wsub     Word substring only needs to match line (not full token)\n");
    printf("   -rann #   Qualify random number of lines #     [NOTE: Exact, NO stdin]\n");
    printf("   -ranf #   Qualify random fraction (0 - 1) #    [NOTE: Exact, NO stdin]\n");
    printf("   -ranp #   Qualify random probability (0 -1) #  [NOTE: Approximate]\n");
    printf("   -seed #   Set random seed to #\n");
    printf("   -maxout # Limit output lines to #\n");
    printf("   -A #      Report # lines After qualifiying lines (like grep -A)\n");
    printf("   -B #      Report # lines Before qualifiying lines [NOTE: NO stdin]\n");
    printf("   -all      All lines considered; Default ignores blank and # comments\n");
    printf("   -blk #    Process lines in 'blocks' of size #; default = 1\n");
    printf("   -bof #    Start blocks offset # lines; default = 0 = skip no lines\n");
    printf("   -blis XXX Comma seperated list of within-block line numbers; first = 1\n");
    printf("   -bnot     Invert within-block line list)\n");
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
    int ok,nok,prev_ok,b_line,line,outline,pout,extra;
    FILTER *filtPO;
    
    filtPO = CreateFilterPO();
    if(!ParseArgsI(argc, argv,
        "S -not B -stat B -rg D2 -gt D -lt D -col I -out S\
        -seed I -ranf D -ranp D\
        -rann I -sc I -flag B -icbn B\
        -lrg I2 -wlis S -kc B -wst B -wsub B\
        -pln B -A I -B I -qu B -vex B -all B -bof I -blk I\
        -brg I2 -bnot B -blis S -maxout I -abs B",
        filtPO->inname, &filtPO->do_not, &filtPO->do_stat, 
        &filtPO->min,&filtPO->max, &filtPO->min, &filtPO->max, 
        &filtPO->col, filtPO->outname, 
        &filtPO->seed, &filtPO->ranf, &filtPO->ranp, 
        &filtPO->rann, &filtPO->skipc, &filtPO->do_flag, 
        &filtPO->do_icbn,
        &filtPO->firstl,&filtPO->lastl, filtPO->wlisname, 
        &filtPO->do_kc, &filtPO->do_wst, &filtPO->do_wsub, 
        &filtPO->do_pln, &filtPO->do_A, &filtPO->do_B, 
        &filtPO->do_quiet, &filtPO->do_vex, &filtPO->do_all, 
        &filtPO->l_bof, &filtPO->l_blk, &filtPO->firstb,&filtPO->lastb,
        &filtPO->do_bm_not, filtPO->blk_mlis, &filtPO->maxout,
        &filtPO->do_abs,
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
    line = outline = nok = extra = pout = 0;
    prev_ok = TRUE;
    while(fgets(bufS,FILTBUF_SIZE,filtPO->in) != NULL) 
    {
        if( SkipThisLineI(filtPO, bufS) ) {
            continue;
        }
        line++;
        /***
        *   If start-of-block, calc and save flag
        */
        b_line = LineInFiltBlockI(filtPO, line);
        if( b_line == 1 ) {
            cPC = GetLineStartPC(filtPO, bufS);
            ok = IsFiltLineOkI(filtPO, line, cPC, TRUE); 
            prev_ok = ok;
        }
        /*  Before blocks start */
        else if ( b_line < 1 ) {
            ok = FALSE;
        }
        /*  In block, but not start; Use previous value, and possibly mask if ok */
        else {
            ok = prev_ok;
            if(ok) {
                ok = LineInBlockOkI(filtPO, b_line);
            }
        }
        /*   Count ok, and set for any -After extra  */
        if(ok) {
            nok++;
            if(filtPO->do_A > 0) {
                extra = filtPO->do_A + 1;
            }
        }
        /*  Not only stats = per-line output */
        if( !filtPO->do_stat) {
            pout = FALSE;
            if( (filtPO->do_flag) || (ok) || (extra>0) ) {
                pout++;
                outline++;
            }
            /* If maxout is set, unset print flag if too many lines */
            if( (filtPO->maxout > 0) && (outline > filtPO->maxout) ) {
                pout = FALSE;
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
    int skip;

    skip = 0;
    /*   comment, blank? */
    if( COM_LINE(bufS) || BlankStringI(bufS) ) {        
        if( ! filtPO->do_all ) {
            skip++;
        }
    }
    /*  Not a data line; Nothing to check? */
    if( ! GetLineStartPC(filtPO, bufS) ) {
        skip++;
    }
    return(skip);
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
        if(SkipThisLineI(filtPO, bufS)) {
            continue;
        }
        line ++;
        cPC = GetLineStartPC(filtPO, bufS);
        if(IsFiltLineOkI(filtPO, line, cPC, rand)) {
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
    CHECK_FREE(fpPO->blk_mask);
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
    fpPO->n_mask = fpPO->n_block = 0;
    INIT_S(fpPO->outname);
    fpPO->out = NULL;
    INIT_S(fpPO->wlisname);
    fpPO->wlis = NULL;
    fpPO->col = DEF_COL;
    fpPO->min = -TOO_BIG_D;
    fpPO->max = TOO_BIG_D;
    fpPO->firstl = fpPO->lastl = BAD_I;
    fpPO->firstb = fpPO->lastb = BAD_I;
    INIT_S(fpPO->blk_mlis);
    fpPO->blk_mask = NULL;
    fpPO->do_bm_not = FALSE;
    fpPO->l_bof = 0;
    fpPO->l_blk = 1;
    fpPO->seed = BAD_I;
    fpPO->rann = -1;
    fpPO->ranf = fpPO->ranp = -1.0;
    fpPO->maxout = -1;
    fpPO->do_not = FALSE;
    fpPO->do_stat = FALSE;
    fpPO->do_kc = fpPO->do_wst = fpPO->do_wsub = FALSE;
    fpPO->do_icbn = FALSE;
    fpPO->skipc = 0;
    fpPO->do_flag = fpPO->do_pln = FALSE;
    fpPO->do_A = fpPO->do_B = 0;
    fpPO->do_quiet = FALSE;
    fpPO->do_vex = FALSE;
    fpPO->do_all = FALSE;
    fpPO->do_abs = FALSE;
}
/*************************************************************************/
int CheckFilterOptionsI(FILTER *filtPO)
{
    Srand(filtPO->seed);
    /*  Set up line block story */
    if(!HandleFilterBlocksI(filtPO)) {
        return(FALSE);
    }
    /*  Set up masking */
    if(!HandleFilterMaskingI(filtPO)) {
        return(FALSE);
    }
    /*  Set qual-line-before masking */
    if(filtPO->do_B > 0) {
        SetBeforeMaskingI(filtPO);
    }
    if(!HandleFilterListSetsI(filtPO)) {
        return(FALSE);
    }
    return(TRUE);
}
/*************************************************************************/
int HandleFilterBlocksI(FILTER *filtPO)
{
    int i;
    char *cPC;

    if( (filtPO->l_blk < 1) || (filtPO->l_blk > MAX_BLK_SIZE) || (filtPO->l_bof < 0) ) {
        printf("Line block size (%d) and offset (%d) won't work\n", filtPO->l_blk, filtPO->l_bof);
        return(FALSE);
    }
    /***
    *   Within-block masking?
    */
    if( (filtPO->l_blk > 1) && (!NO_S(filtPO->blk_mlis)) ) {
        /*  Allocate mask then parse list and set; 1-based so add 1 extra */
        filtPO->blk_mask = (char *)ALLOC(filtPO->l_blk + 1, sizeof(char));
        ReplaceChars(',', filtPO->blk_mlis, ' ', filtPO->blk_mlis);
        cPC = filtPO->blk_mlis;
        while(ISLINE(*cPC))
        {
            i = -1;
            sscanf(cPC,"%d",&i);
            if(i < 1) {
                printf("Problem with block mask list:\n|%s| |%s| %d\n", filtPO->blk_mlis, cPC, i);
                return(FALSE);
            }
            if(i <= filtPO->l_blk) {
                filtPO->blk_mask[i] = 1;
            }
            NEXT_WORD(cPC);
        }
        /* Invert? */
        if(filtPO->do_bm_not) {
            InvertMask(filtPO->blk_mask, filtPO->l_blk + 1);
        }
    }
    return(TRUE);
}
/*************************************************************************
*   Check options and allocate / set mask if needed
*/
int HandleFilterMaskingI(FILTER *filtPO)
{
    int block;
    DOUB fracD;

    if( NeedFilterMaskingI(filtPO) ) {
        if( IsFileStdinI(filtPO->in) ) {
            PROBLINE;
            printf("Options not compatible with stdin\n\n");
            return(FALSE);
        }
        /***
        *   Last two args mean ignore comments and blanks
        *   Then rewind file and allocate space for line mask
        */
        filtPO->n_mask = FileLinesI(filtPO->in,TRUE,TRUE);
        filtPO->n_block = INT((filtPO->n_mask - filtPO->l_bof) / filtPO->l_blk);
        rewind(filtPO->in);
        filtPO->mask = (char *)ALLOC(filtPO->n_mask, sizeof(char));
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
            if( filtPO->l_blk > 1) {
                block = BlockForFiltLineI(filtPO, filtPO->n_mask);
                fracD = RNUM(filtPO->rann) / RNUM(block);
            }
            else {
                fracD = RNUM(filtPO->rann) / RNUM(filtPO->n_mask);
            }
        }
        else if( filtPO->ranf >= 0.0) {
            fracD = filtPO->ranf;
        }
        if(fracD >= 0.0) {
            SetRandFilterMaskingI(filtPO, fracD);
        }
/*
DumpArray(filtPO->mask, IS_CHAR, 0, filtPO->n_mask, " M=%d\n", NULL);
*/
    }
    return(TRUE);
}
/***********************************************************************
*   Returns 1-based block number for line number (i.e. first = 1)
*   Lines assumed 1-based (i.e. first = 1); 
*/
int BlockForFiltLineI(FILTER *filtPO, int line)
{
    int block;

    if( line <= filtPO->l_bof ) {
        block = 0;
    }
    else {
        block = 1 + INT((line - 1 - filtPO->l_bof) / filtPO->l_blk);
    }
    return(block);
}
/************************************************************************
*   Return the 1-based within-block line number for line number; 1 = first
*   Lines assumed 1-based (i.e. first = 1)
*/
int LineInFiltBlockI(FILTER *filtPO, int line)
{
    int bline;

    if( line <= filtPO->l_bof ) {
        bline = 0;
    }
    else {
        bline = 1 + (line - 1 - filtPO->l_bof) % filtPO->l_blk;
    }
    return(bline);
}
/************************************************************************/
int LineInBlockOkI(FILTER *filtPO, int line)
{
    int ok;

    ok = TRUE;
    if( filtPO->blk_mask && (line <= filtPO->l_blk) ) {
        ok = filtPO->blk_mask[line];
    }
    return(ok);
}
/************************************************************************
*   Set random mask; Handles lines and blocks
*/
int SetRandFilterMaskingI(FILTER *filtPO, DOUB fracD)
{
    int i, line;
    char *tmaskPC;

    LIMIT_NUM(fracD, 0.0, 1.0);
    if( filtPO->l_blk > 1) {
        tmaskPC = (char *)ALLOC(filtPO->n_block, sizeof(char));
        MaskRandSubsetI(tmaskPC, filtPO->n_block, RNUM(fracD));
        for(i=0; i<filtPO->n_block; i++)
        {
            line = filtPO->l_bof + (filtPO->l_blk * i);
            BOG_CHECK(line >= filtPO->n_mask);
            if(tmaskPC[i]) {
                filtPO->mask[line] = 1;
            }
        }
        CHECK_FREE(tmaskPC);
    }
    else {
        MaskRandSubsetI(filtPO->mask, filtPO->n_mask, RNUM(fracD));
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
    char *wPC, rwordS[DEF_BS], cwordS[DEF_BS];
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
    *   Copy to clean up extra chars, then try to get number
    */
    sscanf(wPC,"%s",rwordS);
    RemoveChars("[(,)]", rwordS, cwordS);
    rD = BAD_D;
    sscanf(cwordS,"%lf",&rD);
    if(BAD_DOUB(rD)) {
        return(FALSE);
    }
    *valPD = rD;
    return(TRUE);
}
/*************************************************************************
*   Is the line pointed to by cPC OK based on all criteria?
*   line = 1-based line number
*/
int IsFiltLineOkI(FILTER *filtPO, int line, char *cPC, int use_mask) 
{
    char wordS[NSIZE];
    DOUB rD;
    int ok, test, block;

    ok = TRUE;
    test = 0;
    /*  Line number ok? */
    if( !BAD_INT(filtPO->firstl) ) {
        if( (line < filtPO->firstl) || (line > filtPO->lastl) ) {
            ok = FALSE;
        }
        test++;
    }
    /*  Block number ok? */
    if( ok && (!BAD_INT(filtPO->firstb)) ) {
        block = BlockForFiltLineI(filtPO, line);
        if( (block < filtPO->firstb) || (block > filtPO->lastb) ) {
            ok = FALSE;
        }
        test++;
    }
    /*  Mask?  */
    if( ok && (filtPO->mask) && (use_mask) ) {
        /***
        *   If mask is smaller than number of lines, disqualify
        */
        if(line > filtPO->n_mask) {
            ok = FALSE;
        }
        else {
            ok = INT(filtPO->mask[line-1]);
        }
        test++;
    }
    /*  Random probability? */
    if( !test && (filtPO->ranp >= 0.0) ) {
        rD = RandD(1.0);
        ok = (rD < filtPO->ranp);
        test++;
    }
    /***
    *   Get column entry then check 
    *   If can't get word or number from column, simply set not ok
    */
    if( !test ) {
        ok = GetNthWordI(cPC,filtPO->col,wordS);
        /***    
        *   Get number from word (if not token or number list)
        */
        if( ok && (!filtPO->wlis) ) {
            ok = FiltGetWordNumValI(filtPO, wordS, &rD);
            /***
            *   Absolute value?
            */
            if( ok && filtPO->do_abs && (rD < 0.0) ) {
                rD = -rD;
            }
            /***
            *   Too big or small? Inclusive or exclusive bounds differ
            */
            if( ok ) {
                if (( filtPO->do_vex ) && ((rD <= filtPO->min) || (rD >= filtPO->max)) ) {   
                    ok = FALSE; 
                }
                else if ((rD < filtPO->min) || (rD > filtPO->max)) {   
                    ok = FALSE; 
                }
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
