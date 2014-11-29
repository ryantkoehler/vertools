/*
* wordfreq.c
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
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "table.h"
#include "wfutil.h"
#include "wordfreq.h"

#define DB_WF if(DB[99])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(WfUtilI(argc,argv),NULL) ); }
/**************************************************************************/
void WfUtilUse(void)
{
	VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
	printf("   <infile>   Sequence file\n");
	printf("   -out XXX   Set output file to XXX\n");
    printf("   -iraw -ifas Treat input as \"raw\" / fasta format\n");
	printf("   -size #    Set word size to # (def = %d)\n",DEF_WSIZE);
	printf("   -range # # Range of reported words # to # only\n");
	printf("   -norm      Normalize output (average word = 1.0)\n");
	printf("   -deg       Combine counts for +/- degenerate word (rev-comp)\n");
	printf("   -ilc       Ignore Lowercase case sequence\n");
	printf("   -iuc       Ignore Uppercase case sequence\n");
	printf("   -step #    Step between samples (default = 1)\n");
	printf("   -sran # #  Sequence range # to # only\n");
	printf("   -bran # #  Base range # to # only (i.e. sub-seqs)\n");
    printf("   -rre       Base range relative to end; i.e. Backwards\n");
    printf("   -pmat # #  Position-specific matrix from bases # to #\n");
	printf("   -sif XXX   Score input sequences via frequency tab XXX\n");
	printf("   -quiet     Suppress warnings about non-ACGT chars\n");
}
/**************************************************************************
*	Main program
*/
int WfUtilI(int argc, char **argv)
{
	int ok,tal,iraw,ifas,quiet;
	WF_UTIL *wfPO;

	wfPO = CreateWf_utilPO();
	iraw = ifas = quiet = FALSE;
	if(!ParseArgsI(argc,argv,
        "S -out S -iraw B -ifas B -siz I -deg B -bran I2 -rre B\
		-sran I2 -ran I2 -norm B -sif S -pmat I2\
		-ilc B -iuc B -step I -quiet B",
		wfPO->inname, wfPO->outname, &iraw, &ifas, &wfPO->size, 
		&wfPO->do_deg, &wfPO->firstb,&wfPO->lastb, &wfPO->do_rre,
		&wfPO->firsts,&wfPO->lasts, &wfPO->min,&wfPO->max, &wfPO->do_norm,
		&wfPO->lisname, &wfPO->pmat_s,&wfPO->pmat_e,
		&wfPO->do_ilc, &wfPO->do_iuc, &wfPO->step,
        &quiet,
        (int *)NULL))
    {
        WfUtilUse();
		CHECK_WF_UTIL(wfPO);
		return(FALSE);
	}
    /***
    *   Set input format 
    */
    wfPO->iform = FigureSeqFileTypeI(iraw,ifas,wfPO->inname,TRUE);
    if(!wfPO->iform) {
        printf("Problem with input seq(s)\n");
		CHECK_WF_UTIL(wfPO);
        return(FALSE);
    }
	/***
	*	What output options?
	*/
	if(!CheckWfuOptionsI(wfPO)) {
		ABORTLINE;
		CHECK_WF_UTIL(wfPO);
		return(FALSE);
	}
	if(!NO_S(wfPO->lisname)) {
        /* If have a sequence list to score, normalize */
		wfPO->do_norm = TRUE;
	}
	/***
	*	Open in/out files
	*/
	if(!OpenWfuFilesI(wfPO)) {
		ABORTLINE;
		CHECK_WF_UTIL(wfPO);
		return(FALSE);
	}
	/***
	*	Set up auxillary stuff
	*/
	if(!SetUpWfuAuxDataI(wfPO)) {
		ABORTLINE;
		CHECK_WF_UTIL(wfPO);
		return(FALSE);
	}
	SummaryHeader(wfPO, wfPO->out);
	/***
	*	Loop through the input sequence collection
	*/
	wfPO->n = tal = 0;
	while(TRUE) {
		/***
		*	Parse sequence; FALSE = done
		*/
		ok = ParseSeqI(wfPO->in, wfPO->iform, SCLEAN_HI, !quiet, wfPO->seq);
		if(ok==FALSE) {
			break;
		}
		if(ok!=TRUE) {
			if(!wfPO->igprob) {
				break;
			}
			continue;
		}
		/***
		*	Filter here BEFORE any manipulations
		*	Update count and qualify
		*/
		if(!IsCurrentSeqOkI(wfPO,wfPO->n)) {
			wfPO->n += 1;
			continue;
		}
		wfPO->n += 1;
		/***
		*	trim, case-based-mask?
		*/
		if(!HandleWfuSubseqI(wfPO)) {
			ABORTLINE;
			break;
		}
		if( (wfPO->do_ilc) || (wfPO->do_iuc) ) {
			HandleWfuCaseMaskI(wfPO);
		}
		/***
		*	Processing current seq, or Tallying?
		*/
		if(!NO_S(wfPO->lisname)) {
			HandleWfuSeqOutputI(wfPO,wfPO->out);
		}
		else {
			tal = HandleWfTallyI(wfPO);
			if(tal<1) {
				break;
			}
			wfPO->tw += tal;
		}
	}
	/***
	*	Matrix or Stats?
	*/
	if(wfPO->pmat_s > 0) {
 		if(tal>0) {
			HandleWfuPosMat(wfPO,wfPO->out);
		}
	}
	else {
		HandleWfuStats(wfPO,wfPO->out);
	}
	/***
	*	All done
	*/
	CHECK_WF_UTIL(wfPO);
	return(TRUE);
}
/*****************************************************************************
*	Create data struct
*/
WF_UTIL *CreateWf_utilPO()
{
	WF_UTIL *wfPO;

	if(! (wfPO = (WF_UTIL *)ALLOC(1,sizeof(WF_UTIL)) ) )
	{
		printf("# Failed to allocate working object\n");
		return(NULL);
	}
	wfPO->ID = WF_UTIL_ID;
	wfPO->seq = CreateSeqPO(0,NULL,NULL);
	if(!wfPO->seq)
	{
		printf("# Failed to allocate seq space\n");
		CHECK_WF_UTIL(wfPO);
		return(NULL);
	}
	InitWf_util(wfPO);
	return(wfPO);
}
/*****************************************************************************
*	Free datastructure and substructs
*/
int DestroyWf_utilI(WF_UTIL *wfPO)
{
	VALIDATE(wfPO,WF_UTIL_ID);
	CHECK_FILE(wfPO->in);
	CHECK_NFILE(wfPO->out,wfPO->outname);
	CHECK_SEQ(wfPO->seq);
	CHECK_WORDFREQ(wfPO->wf);
	CHECK_TABLE(wfPO->pmat);
	FREE(wfPO);
	return(TRUE);
}
/*****************************************************************************
*	Set null / default values
*/
void InitWf_util(WF_UTIL *wfPO)
{
	VALIDATE(wfPO,WF_UTIL_ID);

	INIT_S(wfPO->inname);
	INIT_S(wfPO->outname);
	wfPO->in = NULL;
	wfPO->iform = BOGUS;
	wfPO->out = NULL;
	wfPO->owhat = BOGUS;
	wfPO->igprob = FALSE;
	wfPO->do_not = FALSE;
	INIT_S(wfPO->lisname);
	wfPO->wf = NULL;
	wfPO->size = DEF_WSIZE;
	wfPO->step = 1;
	wfPO->do_deg = FALSE;
	wfPO->do_rre = FALSE;
	wfPO->firsts = wfPO->firstb = -TOO_BIG;	
	wfPO->lasts = wfPO->lastb = TOO_BIG;
	wfPO->min = 0; 
	wfPO->max = TOO_BIG;
	wfPO->n = 0;
	wfPO->tw = 0;
	wfPO->pmat_s = wfPO->pmat_e = BOGUS;
	wfPO->do_ilc = FALSE;
	wfPO->do_iuc = FALSE;
}
/*************************************************************************
*	Check for option consistency
*/
int CheckWfuOptionsI(WF_UTIL *wfPO)
{
    /***
    *   Not scoring seqs, then size has to be OK
    */
	if(!NO_S(wfPO->lisname)) {
		if(wfPO->size<1) {
			PROBLINE;
			printf("Bogus size: %d\n",wfPO->size);
			return(FALSE);
		}
	}
	/***
	*	Position-specific matrix checks
	*/
	if(wfPO->pmat_s >= 0) {
		if(wfPO->size > MAX_PSM_SIZE) {
			PROBLINE;
			printf("Max word size for Position-specific matrix is %d\n",
				MAX_PSM_SIZE);
			printf("  Sorry, %d is too big\n",wfPO->size);
			return(FALSE);
		}
		if( wfPO->pmat_s > wfPO->pmat_e ) {
			PROBLINE;
			printf("Bad range for Position-specific matrix:\n");
			printf("  Start=%d  End=%d\n", wfPO->pmat_s,wfPO->pmat_e);
			return(FALSE);
		}
	}
	/***
	*	Case ignorning options
	*/
	if( (wfPO->do_ilc) && (wfPO->do_iuc) ) {
		PROBLINE;
		printf("Can't ignore both upper and lower case!\n"); 
		return(FALSE);
	}
	/***
	*	Step size has to be >=1 
	*/
	if(wfPO->step < 1) {
		PROBLINE;
		printf("Step size has to be >=1\n"); 
		return(FALSE);
	}
	return(TRUE);
}
/***************************************************************************
*	Open any needed files or die
*/
int OpenWfuFilesI(WF_UTIL *wfPO)
{
	if(!(wfPO->in=OpenUFilePF(wfPO->inname,"r",NULL))) {
		return(FALSE);
	}
	if(!NO_S(wfPO->outname)) {
		if(!(wfPO->out=OpenUFilePF(wfPO->outname,"w",NULL))) {
			return(FALSE);
		}
	}
	HAND_NFILE(wfPO->out);
	return(TRUE);
}
/***************************************************************************
*	Set up auxillary data structs
*/
int SetUpWfuAuxDataI(WF_UTIL *wfPO)
{
	if(!NO_S(wfPO->lisname)) {
		if(!GetWordFreqsI(wfPO->lisname,&wfPO->wf)) {
			PROBLINE;
			printf("Failed to load word freq table\n");
			return(FALSE);
		}
		if(wfPO->do_norm) {
			NormalizeFrecs(wfPO->wf);
		}
	}
	else {
		wfPO->wf = CreateWordfreqPO(wfPO->size,ALPHDIM); 
        if(!wfPO->wf) {
			PROBLINE;
			printf("Failed to allocate space for wordsize %d\n",wfPO->size);
			return(FALSE);
		}
	}
	/***
	*	Position speicif matrix; rows = N-mer, cols = base position
	*/
	if(wfPO->pmat_s > 0) {
		wfPO->pmat_r = CalcInDimI(wfPO->size,ALPHDIM);
		wfPO->pmat_c = wfPO->pmat_e - wfPO->pmat_s + 1;
		if(! (wfPO->pmat=CreateTablePO(wfPO->pmat_r,wfPO->pmat_c)) ) {
			PROBLINE;
			printf("Failed to allocate postition-specific matrix\n");
			return(FALSE);
		}
		SetPosMatLablesI(wfPO);
	}
	return(TRUE);
}
/***************************************************************************
*	Set row and column lables for position matrix
*/
int SetPosMatLablesI(WF_UTIL *wfPO)
{
	int i;
	char nameS[NSIZE];
	
	/***
	*	Cols = base position relative to 3' or 5' end
	*/
	for(i=0; i<wfPO->pmat_c; i++) {
		sprintf(nameS,"%d",wfPO->pmat_s + i);
		SetTableColLabI(wfPO->pmat,i,nameS);
	}
	/***
	*	Rows = N-mer
	*/
	for(i=0; i<wfPO->pmat_r; i++) {
		FillWordfreqSeqStringI(wfPO->wf,i,nameS);
		SetTableRowLabI(wfPO->pmat,i,nameS);
	}
	/***
	*	Print format; <print> <sep> <preface>
    *   If not normalizing, simple counts are good
	*/
    if(wfPO->do_norm) {
	    SetTablePrintformI(wfPO->pmat,"%6.4f",NULL,"\t",NULL,NULL);
    }
    else {
	    SetTablePrintformI(wfPO->pmat,"%6.0f",NULL,"\t",NULL,NULL);
    }
	return(TRUE);
}
/***************************************************************************
*
*/
int HandleWfuSubseqI(WF_UTIL *wfPO)
{
	int len;

	if(wfPO->firstb >= 0) {
		len = wfPO->lastb - wfPO->firstb + 1;
		if(wfPO->do_rre) {
			NarrowSeqI(wfPO->seq,wfPO->firstb-1,len,REVERSE,FALSE);
		}
		else {
			NarrowSeqI(wfPO->seq,wfPO->firstb-1,len,FORWARD,FALSE);
		}
	}
	return(TRUE);
}
/***************************************************************************
*	If masking lowercase, switch these to n and return the count 
*/
int HandleWfuCaseMaskI(WF_UTIL *wfPO)
{
	int i,n,len;
	char *seqPC;

	len = GetSeqLenI(wfPO->seq);
	GetSeqSeqI(wfPO->seq,&seqPC);
	n = 0;
	for(i=0;i<len;i++) {
		if( (wfPO->do_ilc) && (islower(INT(seqPC[i]))) ) {
			seqPC[i] = 'n';
			n++;
		}
		if( (wfPO->do_iuc) && (isupper(INT(seqPC[i]))) ) {
			seqPC[i] = 'n';
			n++;
		}
	}
	return(n);
}
/**************************************************************************
*	Screen current seq against filters
*/
int IsCurrentSeqOkI(WF_UTIL *wfPO,int n)
{
	int ok;

	ok = TRUE;
	if( (n < wfPO->firsts) || (n > wfPO->lasts) ) {
		ok = FALSE;
	}
	/***
	*	If not, invert qualification
	*/
	if(wfPO->do_not) {
		ok = !ok;
	}
	return(ok);
}
/**************************************************************************
*	Sequence "cleaning" options
*/
void HandleWfuClean(WF_UTIL *wfPO)
{
	SEQ *seqPO;
	int slen;

	seqPO = wfPO->seq;
	slen = CleanUpSeqI(seqPO->seq,seqPO->len,seqPO->seq,FALSE,FALSE);
	seqPO->len = slen;
}
/**************************************************************************
*
*/
void HandleWfuStats(WF_UTIL *wfPO,FILE *outPF)
{
	HAND_NFILE(outPF);
	if(wfPO->n < 1) {
		fprintf(outPF,"NO SEQS!\n");
	}
	if(NO_S(wfPO->lisname)) {
		if(wfPO->do_norm) {
			NormalizeFrecs(wfPO->wf);
		}
		DumpWordsI(wfPO->wf, wfPO->do_deg, wfPO->min, wfPO->max, outPF);
	}
}
/***************************************************************************
*	Score current seq against 
*/
int HandleWfuSeqOutputI(WF_UTIL *wfPO,FILE *outPF)
{
	int i,n,ind;
	char nameS[NSIZE];
	WORDFREQ *wordfPO;
	SEQ *seqPO;
	DOUB sD,minD,maxD;

	HAND_NFILE(outPF);
	seqPO = wfPO->seq;
	wordfPO = wfPO->wf;
    minD = TOO_BIG_R;
    maxD = -TOO_BIG_R;
	sD = 0.0;
	n = 0;
	for(i=0; i<=(seqPO->len - wordfPO->size); i += wfPO->step)
	{
		ind = IndexFromSeqI(&seqPO->seq[i],wordfPO->size,
			wordfPO->n, wordfPO->ald);
		if( (ind >= wordfPO->n) || (ind<0) ) {
			printf("Bogus index: %d\n",ind);
			ERR("HandleWfuSeqOutputI","bad index");
		}
		if(ind < 0) {
			continue;
		}
		maxD = MAX_NUM(maxD, wordfPO->freqs[ind].n);
		minD = MIN_NUM(minD, wordfPO->freqs[ind].n);
		sD += wordfPO->freqs[ind].n;
		n++;
        /***
        *   If degenerate words, need to add compliment too
        */
        if(wfPO->do_deg) {
            ind = CompIndexI(ind, wordfPO->pmax, wordfPO->ald);
		    if( (ind >= wordfPO->n) || (ind<0) ) {
			    printf("Bogus index: %d\n",ind);
			    ERR("HandleWfuSeqOutputI","bad index");
		    }
		    maxD = MAX_NUM(maxD, wordfPO->freqs[ind].n);
		    minD = MIN_NUM(minD, wordfPO->freqs[ind].n);
		    sD += wordfPO->freqs[ind].n;
		    n++;
        }
	}
	FillSeqNameStringI(seqPO,nameS,NSIZE);
	if(n<1) {
	    fprintf(outPF,"# %-15s\tNo words to count!\n",nameS);
        return(FALSE);
	}
	sD /= DNUM(n);
	fprintf(outPF,"%-15s\t%5.4f\t%5.4f\t%5.4f\n",nameS,minD,sD,maxD);
	return(TRUE);
}
/***************************************************************************
*	Report position-specific matrix
*/
void HandleWfuPosMat(WF_UTIL *wfPO, FILE *outPF)
{
	int r,c;
	DOUB vD;

	HAND_NFILE(outPF);
	fprintf(outPF,"# Position-specific word frequencies\n");
	fprintf(outPF,"# Words of size %d\n",wfPO->size);
	fprintf(outPF,"# Bases %d to %d",wfPO->pmat_s,wfPO->pmat_e);
	if(wfPO->do_rre)
	{ 	fprintf(outPF," from 3' end\n"); }
	else
	{ 	fprintf(outPF," from 5' end\n"); }
	/***
	*	Scale shams?
	*/
	if(wfPO->do_norm) {
	    for(r=0;r<wfPO->pmat_r;r++) {
		    for(c=0;c<wfPO->pmat_c;c++) {
			    GetTableValI(wfPO->pmat,r,c,&vD);
			    vD = vD / DNUM(wfPO->n);
			    SetTableValI(wfPO->pmat,r,c,vD);
            }
		}
	}	
	/***
	*	Dump table
	*/
	DumpTable(wfPO->pmat,FALSE,FALSE,outPF);
}
/***************************************************************************
*	Tally words for current sequence
*/
int HandleWfTallyI(WF_UTIL *wfPO)
{
	int n;
	SEQ *seqPO;

	seqPO = wfPO->seq;
	if(seqPO->len < wfPO->size) {
		return(0);
	}
	/***
	*	Position matrix or just word talley
	*/
	if(wfPO->pmat_s > 0) {
		n = TallyPosMatWordsI(wfPO,seqPO);
	}
	else {
		n = TallyWordsI(seqPO,wfPO->wf,wfPO->step);
	}
	return(n);
}
/***************************************************************************
*	Update postition-specific matrix table with current seq
*	If problem, return FALSE
*/
int TallyPosMatWordsI(WF_UTIL *wfPO, SEQ *seqPO)
{
	int i,s,ind,len;
	char *seqPC;
	WORDFREQ *wordPO;

	wordPO = wfPO->wf;
	len = GetSeqLenI(seqPO);
	GetSeqSeqI(seqPO,&seqPC);
	/***
	*	For each position (col of matrix)
	*/
	for(i=0; i < wfPO->pmat_c; i++) {
		/***
		*	Where to start sampling from?
        *   Given start coord is 1-based, so subtract 1
		*/
		if(wfPO->do_rre) {
			s = len - (wfPO->pmat_s -1) - i - wfPO->size;
		}
		else {
			s = (wfPO->pmat_s -1) + i;
		}
		/***
		*	Out of sampling room for seq len?
		*/
		if( (s<0) || (s>(len-wfPO->size)) ) {
			break;
		}
		/***
		*	Get word index (row of matrix) and update table
		*/
		ind = IndexFromSeqI(&seqPC[s], wordPO->size, wordPO->n, wordPO->ald);
		if(!ModTableValI(wfPO->pmat,ind,i,1.0,MATH_ADD)) {
			printf("Row=ind=%d col=i=%d\n",ind,i);
			ERR("TallyPosMatWordsI","ModTableValI failed");
			return(FALSE);
		}
	}
	return(TRUE);
}
/**************************************************************************
*	Write header info for wordfreq data
*/
void SummaryHeader(WF_UTIL *wfPO, FILE *outPF)
{
	HAND_NFILE(outPF);
	VersionSplash(outPF,VERSION_S,"#  ",TRUE);
	fprintf(outPF,"# Input        %s\n",wfPO->inname);
	if(wfPO->firstb > 0) {
		if(wfPO->do_rre) {
			fprintf(outPF,"#   Last %d to %d bases only\n",wfPO->firstb,wfPO->lastb);
		}
		else {
			fprintf(outPF,"#   First %d to %d bases only\n",wfPO->firstb,wfPO->lastb);
		}
	}
	if(!NO_S(wfPO->lisname)) {
	    fprintf(outPF,"# Scoring sequences by normalized frequencies\n");
	    fprintf(outPF,"# Word freq table: %s\n",wfPO->lisname);
	    fprintf(outPF,"# <name>\t<min>\t<mean>\t<max>\n");
    }
    else {
        WordfreqSummary(wfPO->wf, wfPO->min,wfPO->max, outPF);
	    fprintf(outPF,"#\n");
    }
}
