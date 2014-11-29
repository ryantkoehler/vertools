/*
* filter.c
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
	printf("   -lrg # #  Qualify line number range # to #\n");
	printf("   -wlis XXX Qualify line with words (first token) listed in XXX\n");
	printf("   -kc       Keep case for token comparison (default ignore)\n");
	printf("   -wst      Word start only needs to match line (not full token)\n");
	printf("   -wsub     Word substring only needs to match line (not full token)\n");
	printf("   -rann #   Qualify random number of lines #     [NOTE: exact, NO stdin]\n");
	printf("   -ranf #   Qualify random fraction (0 - 1) #    [NOTE: exact, NO stdin]\n");
	printf("   -ranm #   Set random mask dimension to # (rather than count input lines)\n");
	printf("   -ranp #   Qualify random probability (0 -1) #  [NOTE: approximate]\n");
	printf("   -seed #   Set random seed to #\n");
	printf("   -icbn     Ignore chars before numbers (i.e. strip leading chars)\n");
	printf("   -not      Invert line qualification test(s)\n");
	printf("   -flag     Preceed lines with 1/0 for good/bad\n");
	printf("   -pln      Preceed lines with line number\n");
	printf("   -stat     Report only stats about values\n");
	printf("\n");
	printf("NOTE: Filtered file lines limited to %d chars wide\n",FILTBUF_SIZE);
	printf("      Also, -rann and -ranf fractions limited by -ranm 'M' to at most 'M'\n");
	printf("\n");
}
/**************************************************************************/
int Filter_numsI(int argc, char **argv)
{
	char bufS[FILTBUF_SIZE+1], *cPC;
	int j,ok,num,nok,line,tline;
    FILTER *filtPO;
    
    filtPO = CreateFilterPO();
	if(!ParseArgsI(argc, argv,
		"S -not B -stat B -rg D2 -gt D -lt D -col I -out S\
		-seed I -ranf D -ranp D\
		-rann I -sc I -flag B -icbn B\
		-lrg I2 -wlis S -kc B -tst B -tsub B\
        -ranm I -pln B",
		filtPO->inname, &filtPO->do_not, &filtPO->do_stat, 
        &filtPO->min,&filtPO->max, &filtPO->min, &filtPO->max, 
        &filtPO->col, filtPO->outname, 
        &filtPO->seed, &filtPO->ranf, &filtPO->ranp, 
        &filtPO->rann, &filtPO->skipc, &filtPO->do_flag, 
        &filtPO->do_icbn,
		&filtPO->firstl,&filtPO->lastl, filtPO->wlisname, 
        &filtPO->do_kc, &filtPO->do_wst, &filtPO->do_wsub, 
        &filtPO->ranm, &filtPO->do_pln,
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
	*	Filter lines 
	*/
	tline = line = num = nok = 0;
	while(fgets(bufS,FILTBUF_SIZE,filtPO->in) != NULL) 
    {
		tline++;
		/***
		*	Only non-comment, non-blank lines count 
		*/
		if( COM_LINE(bufS) || BlankStringI(bufS) ) {		
			continue;	
		}
		line++;
		/***
		*	Jump past any up front chars to ingore 
		*/
		cPC = bufS;
        for(j=0;j<filtPO->skipc;j++)
        {
            if(!ISLINE(*cPC)) {
                break;
			}
            cPC++;
        }
        if(!ISLINE(*cPC)) {
			continue;
		}
        /***
        *   Qualified?
        */
        ok = IsFiltLineOkI(filtPO, line, cPC); 
        num++;
		/***
		*	Invert status
		*/
		if(filtPO->do_not) {
			ok = !ok;
		}
		nok += ok;
        /***
        *   Now output
        */
        if(filtPO->do_stat) {
            continue;
        }
		if(filtPO->do_flag) {
			if(filtPO->do_pln) {
                fprintf(filtPO->out,"%d\t",line);
            }
			fprintf(filtPO->out,"%d\t",ok);
		    fputs(bufS,filtPO->out);
		}
        else if(ok) {
			if(filtPO->do_pln) {
                fprintf(filtPO->out,"%d\t",line);
            }
		    fputs(bufS,filtPO->out);
        }
	}
    ReportFilterStats(filtPO, num, nok); 
    CHECK_FILTER(filtPO);
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
    fpPO->rann = fpPO->ranm = BAD_I;
    fpPO->ranf = fpPO->ranp = BAD_D;
    fpPO->do_not = FALSE;
    fpPO->do_stat = FALSE;
    fpPO->do_kc = fpPO->do_wst = fpPO->do_wsub = FALSE;
    fpPO->do_icbn = FALSE;
    fpPO->skipc = 0;
    fpPO->do_flag = fpPO->do_pln = FALSE;
}
/*************************************************************************
*
*/
int CheckFilterOptionsI(FILTER *filtPO)
{
    /***
    *   Random set and masks
    */
	Srand(filtPO->seed);
    if(!HandleFilterMaskingI(filtPO)) {
        return(FALSE);
    }
    if(!HandleFilterListSetsI(filtPO)) {
        return(FALSE);
    }
    return(TRUE);
}
/*************************************************************************
*   Random number or fraction cases need a per-line mask
*/
int HandleFilterMaskingI(FILTER *filtPO)
{
    int n;
    DOUB fracD;

    n = 0;
    /***
    *   If need mask, bound or count for dimension
    */
    if(!BAD_INT(filtPO->ranm)) {
        n = filtPO->ranm;
        LIMIT_NUM(n, 1, TOO_BIG);
    }
	else if ( (!BAD_INT(filtPO->rann)) || (!BAD_DOUB(filtPO->ranf)) ) {
        /*
	    *	Last two args mean ignore comments and blanks
        *   Then rewind file and allocate space for lines
        */
	    n = FileLinesI(filtPO->in,TRUE,TRUE);
	    rewind(filtPO->in);
    }
    /***
    *   Anything to allocate and set?
    */
    if(n>0) {
        filtPO->n_mask = n;
	    filtPO->mask = (char *)ALLOC(filtPO->n_mask,sizeof(char));
        if(!filtPO->mask) {
		    PROBLINE;
		    printf("Can't allocate mask for %d input lines\n",filtPO->n_mask);
            return(FALSE);
        }
        /***
        *   Calculate / use fraction to mask subset
        */
	    if(!BAD_INT(filtPO->rann)) {
		    fracD = RNUM(filtPO->rann) / RNUM(filtPO->n_mask);
        }
        else {
		    fracD = filtPO->ranf;
        }
        LIMIT_NUM(fracD, 0.0, 1.0);
		MaskRandSubsetI(filtPO->mask,filtPO->n_mask,RNUM(fracD));
    }
    return(TRUE);
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
/*************************************************************************
*
*/
void ReportFilterStats(FILTER *filtPO, int num, int nok) 
{
    printf("#  %d data lines\n",num);
    printf("#  %d qualified data lines (%5.2f%%)\n",
        nok, PERCENT_R(nok,num));	
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
	*	Ignore chars before digits?
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
	*	Get number 
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
int IsFiltLineOkI(FILTER *filtPO, int line, char *cPC) 
{
    char wordS[NSIZE];
    DOUB rD;
    int ok;

	/***
	*	Just filtering on line number
	*/
	ok = TRUE;
    if( !BAD_INT(filtPO->firstl) ) {
		if( (line < filtPO->firstl) || (line > filtPO->lastl) ) {
		    ok = FALSE;
		}
    }
	/***
	*	Random fractions
    *   If we've got a mask (i.e. random fractions), use that
	*/
    else if(filtPO->mask) {
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
	else if(!BAD_INT(filtPO->ranp)) {
        rD = RandD(1.0);
        ok = (rD < filtPO->ranp);
    }
	/***
	*	Get column entry then check 
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
    		*	Too big or small?
    		*/
    		if( ok && ( (rD < filtPO->min) || (rD > filtPO->max)) ) {	
    			ok = FALSE;	
    		}
	    }
	    /***
		*	In number or token list?
		*/
	    if( ok && (filtPO->wlis) ) {
            ok = WordInWordlistI(filtPO->wlis, wordS, filtPO->do_kc, filtPO->do_wst, 
                                filtPO->do_wsub, NULL);
        }
	}
    return(ok);
}
