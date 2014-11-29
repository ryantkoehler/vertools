/*
* comp_seq.c
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


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "dna_pair.h"
#include "comp_seq.h"

#define DB_COMSEQ 	if(DB[114])

#define MAX_LINE 70

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit(AllDoneI(Comp_seqI(argc,argv),NULL)); }
/**************************************************************************/
void Comp_seqUse(void)
{
	VersionSplash(NULL,VERSION_S,"#  ",TRUE);
	printf("Use: <infile> ['-' for stdin] [-sim|-com|-ham|-self|-psel|-loop] [...options]\n");
	printf("   <infile>     Is a sequence file (\"raw\" format)\n");
	printf("   -iraw -ifas  Treat input as raw / fasta format\n");
	printf("   -out XXX     Set output to XXX\n");
	printf("   -self -psel  Evaluate for self complimentarity / parallel\n");
	printf("   -loop #      Evaluate hairpin with min \"loop\" size #\n");
	printf("   -sim         Evaluate similarity between seqs (best no-gap align)\n");
	printf("   -com         Evaluate complementarity between seqs (best no-gap align)\n");
	printf("   -ham         Evaluate \"Hamming distance\" (Similarity & No align)\n");
	printf("   -smat -scon  Score via Simple match / Contiguous match\n");
	printf("   -sco3 -sco5  Score via contiguous match at 3' / 5' end\n");
	printf("   -scb	        Score via Block-weighted match\n");
	printf("   -swm         Score via thermo weighted matching (Only complement)\n");
	printf("   -sw32        Score via GC=3 AT=2 matching (Only complement)\n");
	printf("   -mwf XXX     Match weights from XXX (AA x, AC x, one pair/line)\n");
	printf("   -rset XXX    Compare infile seqs to reference set seqs XXX\n");
	printf("   -crs         Compare rset in serial order (not full set)\n");
	printf("   -cl3 -cl5 #  Force 3' / 5' end \"clamps\" in matches\n");
	printf("   -mwd # -th # Set minimum word size / threshold to #\n");
	printf("   -tot -num    Report over-thresh Total / Number as score\n");
/*
	printf("   -ostat       Report ONLY stats (low av high scores)\n");
*/
	printf("   -fmat        Report full pair-wise score matrix\n");
	printf("   -rm          Report matches (high score with position offset)\n");
	printf("   -sa          Show alignments (of high score match)\n");
	printf("   -alf #       Set threhold for full alignment display\n");
	printf("   -alp #       Set threhold for partial alignment display\n");
	printf("   -norm        Normalize scores by maximum possible\n");
	printf("   -flg # #     Flag sequences with scores in the range # to #\n");
	printf("   -eraw        Extract flagged sequences in raw format\n");
	printf("   -not         Invert qualification for flagging / extraction\n");
    printf("\n");
}
/**************************************************************************
*	Main program
*/
int Comp_seqI(int argc, char **argv)
{
	int scon,swm,scb,sw32,smat,cl3,cl5,sco3,sco5,iraw,ifas;
	int r_max,r_num,r_tot,mword,quiet,nsco;
	int i,flen,slen;
	char smfS[DEF_BS], *fPC, *sPC, fnameS[DEF_BS];
	COMPSEQ *csPO;
	SEQ *seqPO, *sseqPO;
	REAL scR,thR,partmatR;

	csPO = CreateCompseqPO();
	scon = swm = sw32 = scb = smat = FALSE;
	sco3 = sco5 = r_tot = r_max = r_num = quiet = iraw = ifas = FALSE;
	cl3 = cl5 = 0;
	thR = 0.0;
	partmatR = BAD_R;
	mword = 1;
	INIT_S(smfS);
	if(!ParseArgsI(argc,argv,
        "S -out S -sa B -sim B -com B -self B\
		-scon B -rm B -swm B -rset S -mwf S -flg R2 -scb B\
		-smat B -loop I -norm B -cl3 I -cl5 I -sco3 B -sco5 B -num B\
		-tot B -th R -mwd I -fmat B -psel B -crs B -not B -eraw B \
		-iraw B -ifas B -alf R -alp R -sw32 B -ham B",
		csPO->inname, csPO->outname, &csPO->do_sa, &csPO->do_sim, &csPO->do_com,
		&csPO->do_self, &scon, &csPO->do_rm, &swm, csPO->rinname,
		smfS, &csPO->flo,&csPO->fhi, 
		&scb, &smat, &csPO->do_loop,&csPO->do_norm,
		&cl3, &cl5, &sco3, &sco5, &r_num, &r_tot,
		&thR, &mword, &csPO->do_fmat, &csPO->do_pself, &csPO->do_crs, 
		&csPO->do_not, &csPO->do_eraw, &iraw, &ifas,
		&csPO->fullmat, &partmatR, &sw32, &csPO->do_ham,
        (int *)NULL))
    {
        Comp_seqUse();
		CHECK_COMPSEQ(csPO);
		return(TRUE);
	}
	/***
	*	How are we comparing / what output?
	*/
	if(csPO->do_loop>0)
	{ 	csPO->owhat = CSO_LOOP; }
	else if(csPO->do_pself)
	{ 	csPO->owhat = CSO_PSELF; }
	else if(csPO->do_self)
	{ 	csPO->owhat = CSO_SELF; }
	else if(csPO->do_crs)
	{ 	csPO->owhat = CSO_RSER; }
	else 
	{ 	csPO->owhat = CSO_FULL; }

    if(csPO->do_ham) {
        SetPparsHammingI(csPO->pp, csPO->do_ham);
	    csPO->do_sim = TRUE;	
    }
	/***
	*	Load / set comparison algorithms; 
	*/
	if(csPO->do_sim) { 	
        SetPair_parIdentMatch(csPO->pp); 
    }
	else { 	
        SetPair_parCompMatch(csPO->pp); 
    }
	if( (swm) || (sw32) || (!NO_S(smfS)) ) {	
		if(!NO_S(smfS)) {
			if(!LoadPairingParsI(smfS,csPO->pp)) {
				CHECK_COMPSEQ(csPO);
				return(FALSE);
			}
		}
		else { 	
            if(sw32) { 
                SetPair_parS3W2Wmatch(csPO->pp); 
            }
            else {
			    SetPair_parThermoWmatch(csPO->pp); 
			    if(BAD_REAL(partmatR)) {
				    csPO->partmat = -0.05;
			    }
            }
		}
	}
	if(!BAD_REAL(partmatR)) {
		csPO->partmat = partmatR;
	}
	/***
	*	Default algorithm case = simple match; arg "smat" is this case
	*/
	if(scb) {	
        SetPparsAlgI(csPO->pp,ALGO_BMATCH);	
    }
	else if(scon) {	
        SetPparsAlgI(csPO->pp,ALGO_CONT);	
    }
	else if(sco3) {	
        SetPparsAlgI(csPO->pp,ALGO_CON3);	
    }
	else if(sco5) {	
        SetPparsAlgI(csPO->pp,ALGO_CON5);	
    }
	else if(smat || TRUE) {	
        SetPparsAlgI(csPO->pp,ALGO_MATCH);	
    }
	SetPparsClampsI(csPO->pp,cl3,cl5);
	SetPparsWordMinI(csPO->pp,mword);
	/***
	*	Output-specific settings
	*/
	switch(csPO->owhat)
	{
		case CSO_LOOP:
			SetPparsLoopMinI(csPO->pp,csPO->do_loop);
			SetPparsCtypeI(csPO->pp,PP_COM_LOOP);
			csPO->do_self = TRUE;
			break;
		case CSO_PSELF:
			SetPparsCtypeI(csPO->pp,PP_COM_PCOM);
			csPO->do_self = TRUE;
			break;
		case CSO_SELF:
			SetPparsCtypeI(csPO->pp,PP_COM_COM);
			csPO->do_self = TRUE;
			break;
		case CSO_RSER:
		case CSO_FULL:
			if(csPO->do_com)
			{ 	SetPparsCtypeI(csPO->pp,PP_COM_COM); }
			else if(csPO->do_sim)
			{ 	SetPparsCtypeI(csPO->pp,PP_COM_SIM); }
			csPO->do_self = FALSE;
			break;
	}
	/***
	*	What score to report?	
	*	r_max == default
	*/
	if(r_tot) { 	
        SetPparsRvalI(csPO->pp,PP_R_TOT); 
    }
	else if(r_num) { 	
        SetPparsRvalI(csPO->pp,PP_R_NUM); 
    }
	else if(r_max || TRUE)	/* 	xxx SHAM? */ { 	
        SetPparsRvalI(csPO->pp,PP_R_MAX); 
    }
	SetPparsThreshI(csPO->pp,thR);
	/***
	*	Get input format; guess then override with command line
	*/
    csPO->iform = FigureSeqFileTypeI(iraw,ifas,csPO->inname,TRUE);
	if(!csPO->iform) {
		CHECK_COMPSEQ(csPO);
		ABORTLINE;
		return(FALSE);
	}
	/***
	*	Get the sequences / input file; and allocate space for results
	*/
	if(!GetCompSeqSeqsI(csPO,quiet)) {
		CHECK_COMPSEQ(csPO);
		ABORTLINE;
		return(FALSE);
	}
	if(!AddCompseqScoreSpaceI(csPO)) {
		CHECK_COMPSEQ(csPO);
		ABORTLINE;
		return(FALSE);
	}
	/***
	*	Check options now that we have sequences
	*/
	if( (!csPO->do_sa) && (!csPO->do_rm) ) {
		csPO->ostat = TRUE;
	}
	if(!CheckCompseqOptionsI(csPO)) {
		CHECK_COMPSEQ(csPO);
		ABORTLINE;
		return(FALSE);
	}
	/***
	*	Output files?
	*/
	if(!OpenCompseqOutputI(csPO)) {
		CHECK_COMPSEQ(csPO);
		ABORTLINE;
		return(FALSE);
	}
	/***
	*	Write headers
	*/
	WriteCompseqHeader(csPO,csPO->out); 
	/***
	*	Run the comparisions
	*/
	DB_COMSEQ DB_PrI("+ looking at nseqs %d\n",csPO->nseqs);
	i = 0;
	while(TRUE)
	{
		if(!GetNextSeqI(csPO,i,&seqPO)) {
			break;
		}
		flen = GetSeqLenI(seqPO);
/* xxx
        if(csPO->do_ham) {
            SetPparsWordMinI(csPO->pp,flen);
        }
*/
		GetSeqSeqI(seqPO,&fPC);
		FillSeqNameStringI(seqPO,fnameS,30);
		/***
		*	Get pair-wise scores
		*/
		DB_COMSEQ DB_PrI("+ seq[%d]",i);
		switch(csPO->owhat)
		{
			case CSO_LOOP:
				DB_COMSEQ DB_PrI(" loop\n");
				csPO->offs[0] = ScoreLoopOverlapI(fPC,flen,csPO->pp,&scR);
				csPO->scores[0] = scR;
				break;
			case CSO_SELF:
			case CSO_PSELF:
				DB_COMSEQ DB_PrI(" [p]self\n");
				csPO->offs[0] = ScoreSeqCompareI(fPC,flen,fPC,flen,
					csPO->pp,&scR);
				csPO->scores[0] = scR;
				break;
			case CSO_RSER:
				DB_COMSEQ DB_PrI(" rser\n");
				GetSeqsetSeqStringI(csPO->target,i,&sPC);				
				slen = GetSeqsetSeqLenI(csPO->target,i);
				csPO->offs[0] = ScoreSeqCompareI(fPC,flen,sPC,slen,
					csPO->pp,&scR);
				csPO->scores[0] = scR;
				break;
			case CSO_FULL:
				DB_COMSEQ DB_PrI(" full\n");
				CompareSeqToSeqsetI(seqPO, csPO->target, csPO->pp, csPO->scores, 
					csPO->offs);
				break;
		}
		/***
		*	Process score collection and different outputs
		*/
		nsco = ProcessScoresI(csPO,i,fPC,flen);
		if(HandleSeqFlaggingI(csPO,seqPO,csPO->out))
		{ 	
			i++;
			continue; 	
		}
		if(HandleMatrixOutI(csPO,fnameS,csPO->out))
		{ 	
			i++;
			continue; 	
		}
		switch(csPO->owhat)
		{
			case CSO_LOOP:
				HandleSeqLoopOut(csPO,seqPO,0);
				break;
			case CSO_PSELF:
				HandleSeqPairOut(csPO,seqPO,NULL,0,FALSE,TRUE);
				break;
			case CSO_SELF:
				HandleSeqPairOut(csPO,seqPO,NULL,0,TRUE,FALSE);
				break;
			case CSO_RSER:
                GetSeqsetSeqI(csPO->target,i,&sseqPO);
				if(csPO->do_com) 
				{
					HandleSeqPairOut(csPO,seqPO,sseqPO,0,TRUE,FALSE);
				}
				else
				{
					HandleSeqPairOut(csPO,seqPO,sseqPO,0,FALSE,FALSE);
				}
				break;
			case CSO_FULL:
				HandleFullScoreOut(csPO,nsco,i,seqPO);
		}
		/***
		*	Next seq
		*/
		i++;
	}
	/***
	*	All done
	*/
	CHECK_COMPSEQ(csPO);
	return(TRUE);
}
/*****************************************************************************
*	Return the next seq to be processed, or FALSE if at end of list
*/
int GetNextSeqI(COMPSEQ *csPO,int s, SEQ **seqPPO)
{
	SEQ *seqPO;

	if(csPO->do_self) {
		if(!ParseSeqI(csPO->in,csPO->iform,TRUE,TRUE,csPO->seq)) {
			return(FALSE);
		}
		seqPO = csPO->seq;
	}
	else {
		if(!GetSeqsetSeqI(csPO->seqs,s,&seqPO)) {
			return(FALSE);
		}
	}
	*seqPPO = seqPO;
	return(TRUE);
}
/*****************************************************************************/
void HandleNoScoreOut(COMPSEQ *csPO,SEQ *seqPO)
{
	char fnameS[NSIZE];

	FillSeqNameStringI(seqPO,fnameS,25);
	fprintf(csPO->out,"%-10s = 0 (NO SCORE)\n",fnameS); 
}
/*****************************************************************************
*	Output sequence pair info
*/
void HandleSeqPairOut(COMPSEQ *csPO, SEQ *fseqPO, SEQ *sseqPO, int sc, int
	comp, int antip)
{
	char fseqS[MAX_CSL],sseqS[MAX_CSL],pS[MAX_CSL],compS[MAX_CSL];
	char *fPC,*sPC,fnameS[NSIZE],snameS[NSIZE],formS[DEF_BS];
	int i,n,fi,si,p,flen,slen,off;
	REAL scR,alscR;

	DB_COMSEQ DB_PrI(">> HandleSeqPairOut sc=%d com=%d ap=%d\n",sc,comp,antip);
	HAND_NFILE(csPO->out);
	off = csPO->offs[sc];
	scR = csPO->scores[sc];
	DB_COMSEQ DB_PrI("+ scR=%f off=%d\n",scR,off);
	/***
	*	Get names, lens, seq pointers
	*/
	FillSeqNameStringI(fseqPO,fnameS,30);
	if(csPO->ostat)
	{
		fprintf(csPO->out,"%-14s %7.3f\n",fnameS,scR);
		return;
	}
	flen = GetSeqLenI(fseqPO);
	GetSeqSeqI(fseqPO,&fPC);
	if(sseqPO!=NULL) {
		FillSeqNameStringI(sseqPO,snameS,30);
		slen = GetSeqLenI(sseqPO);
		GetSeqSeqI(sseqPO,&sPC);
	}
	else {
		FillSeqNameStringI(fseqPO,snameS,30);
		slen = GetSeqLenI(fseqPO);
		GetSeqSeqI(fseqPO,&sPC);
	}
	/***
	*	If failing / bogus score then don't draw anything
	*/
	if(scR < 0.0)
	{
		if(csPO->do_rm) {
			fprintf(csPO->out,"%-10s -vs- %-10s = 0 (NO SCORE)\n",
				fnameS,snameS); 
		}
		return;
	}
	/***
	*	What will be printed as comparison? actual or inverse (3' to 5')
	*/
	if(comp) { 	
		ReverseDNASeqI(sPC,slen,compS); 
	}
	else { 	
        strncpy(compS,sPC,slen);
        compS[slen] = '\0';
	}
	/***
	*	Build up front part of match strings
	*/
	si = fi = p = n = 0;
	for(i=0;i<ABS_VAL(off);i++)
	{
		if(off < 0) {
			fseqS[p] = ' ';
			sseqS[p] = compS[si++];
		}
		else {
			fseqS[p] = fPC[fi++];
			sseqS[p] = ' ';
		}
		pS[p] = ' ';
		p++;
	}
	/***
	*	Build matching part, counting observed matchs
	*/
	while(TRUE)
	{
		if( (fi >= flen) || (si >= slen) ) {
			break;
		}
		fseqS[p] = fPC[fi];
		sseqS[p] = compS[si];
		if(GetPair_parMatchScoreI(csPO->pp,fseqS[p],sseqS[p],&alscR))
		{
			if(alscR >= csPO->fullmat)
			{
				pS[p] = '|';
				n++;
			}
			else if(alscR >= csPO->partmat)
			{
				pS[p] = ':';
				n++;
			}
			else
			{
				pS[p] = ' ';
			}
		}
		else
		{
			pS[p] = 'X';
		}
		fi++;
		si++;
		p++;
	}
	/***
	*	Anything left?
	*/
	while( (fi < flen) || (si < slen) )
	{
		fseqS[p] = sseqS[p] = ' ';
		if(fi < flen)
		{
			fseqS[p] = fPC[fi];
		}
		if(si < slen)
		{
			sseqS[p] = compS[si];
		}
		pS[p] = ' ';
		fi++;
		si++;
		p++;
	}
	fseqS[p] = sseqS[p] = pS[p] = '\0';
	/***
	*	Summary statement?
	*/
	if(csPO->do_rm)
	{
		fprintf(csPO->out,"%-10s -vs- %-10s = %6.4f",fnameS,snameS,scR); 
		fprintf(csPO->out," @[%d]",off);
		fprintf(csPO->out," %5d/%d\n",n,flen);
	}
	/***
	*	If not actually showing alignment, bail
	*/
	if(!csPO->do_sa)
	{
		return;
	}
	/***
	*	Ajdustable format string for different name lengths
	*/
	n = MAX_NUM(strlen(fnameS),strlen(snameS));
	LIMIT_NUM(n,10,40);
	sprintf(formS,"%%-%ds",n);
	/***
	*	Report the match
	*/
	if(comp)
	{
		fprintf(csPO->out,formS,fnameS);
		fprintf(csPO->out," >5' %s\n",fseqS);
		fprintf(csPO->out,formS,"matching:");
		fprintf(csPO->out,"     %s\n",pS);
		fprintf(csPO->out,formS,snameS);
		fprintf(csPO->out," <3' %s\n",sseqS);
	}
	else
	{
		fprintf(csPO->out,formS,fnameS);
		fprintf(csPO->out," >5' %s\n",fseqS);
		fprintf(csPO->out,formS,"matching:");
		fprintf(csPO->out,"     %s\n",pS);
		fprintf(csPO->out,formS,snameS);
		fprintf(csPO->out," >5' %s\n",sseqS);
	}
	fprintf(csPO->out,"\n");
}
/****************************************************************************
*	Make pretty picture of a loop'd sequence
*	off = offset of 3'end relative to 5' end.
*/
void HandleSeqLoopOut(COMPSEQ *csPO, SEQ *seqPO, int sc)
{
	int i,p,n,ls,fi,si,len,off;
	char startS[MAX_CSL],endS[MAX_CSL],pS[MAX_CSL];
	char nameS[DEF_BS], *sPC;
	REAL scR; 

	off = csPO->offs[sc];
	scR = csPO->scores[sc];
	HAND_NFILE(csPO->out);
	/***
	*	Sequence pointer, len, name
	*/
	len = GetSeqLenI(seqPO);
	GetSeqSeqI(seqPO,&sPC);
	FillSeqNameStringI(seqPO,nameS,30);
	/***
	*	If only stats, report name and score
	*/
	if(csPO->ostat)
	{
		fprintf(csPO->out,"%-14s %7.3f\n",nameS,scR);
		return;
	}
	/***
	*	Fill in any overhaning parts
	*/
	si = len -1;
	p = fi = 0;
	for(i=0;i<ABS_VAL(off);i++)
	{
		if(off < 0)
		{
			startS[p] = ' ';
			endS[p] = sPC[si--];
		}
		else
		{
			startS[p] = sPC[fi++];
			endS[p] = ' ';
		}
		pS[p] = ' ';
		p++;
	}
	/***
	*	Fill in matching "stem" parts
	*/
	n = 0;
	ls = (len - csPO->pp->min_loop - ABS_VAL(off))/2;
	for(i=0;i<ls;i++)
	{
		startS[p] = sPC[fi++];
		endS[p] = sPC[si--];
		if(GoodDNACompBasesI(startS[p],endS[p]))
		{
			pS[p] = '|';
			n++;
		}
		else
		{
			pS[p] = ' ';
		}
		p++;
	}
	/***
	*	End loop stuff
	*/
	startS[p] = endS[p] = '-';
	pS[p++] = ' ';
	ls = (csPO->pp->min_loop +1)/2;
	for(i=0;i<ls; i++)
	{
		startS[p] = sPC[fi++];
		if((ABS_VAL(off + len) % 2)>0)
		{
			if(i==0)
				endS[p] = '-';
			endS[p+1] = sPC[si--];
		}
		else
		{
			endS[p] = sPC[si--];
		}
		pS[p] = ' ';
		p++;
	}
	pS[p-1] = ':';
	/***
	*	Final output
	*/
	startS[p] = endS[p] = pS[p] = '\0';
	p = 0;
	if(csPO->do_rm)
	{
		fprintf(csPO->out,"%s LoopScore = %6.4f",nameS,scR); 
/** SHAM
		fprintf(csPO->out," %4d/%d",n,len);
**/
		fprintf(csPO->out," @[%d]\n",off);
		p++;
	}
	if(csPO->do_sa)
	{
		fprintf(csPO->out,"%-12s 5'> %s\n",nameS,startS);
		fprintf(csPO->out,"%-12s     %s\n",nameS,pS);
		fprintf(csPO->out,"%-12s 3'< %s\n",nameS,endS);
		p++;
	}
	if(p)
	{
		fprintf(csPO->out,"\n");
	}
}
/***************************************************************************
*	Load sequences 
*	If self-compare, only open the input file to eat one at time
*/
int GetCompSeqSeqsI(COMPSEQ *csPO,int quiet)
{
	SEQSET *ssPO, *rssPO;
	int max;

	DB_COMSEQ DB_PrI(">> GetCompSeqSeqsI\n");
	/***
	*	Special case for self-compare; only open input file, don't load
	*/
	if(csPO->do_self)
	{
		if(!(csPO->in=OpenUFilePF(csPO->inname,"r",NULL)))
		{
			return(FALSE);
		}
		csPO->nseqs = 1;
		return(TRUE);
	}
	/***
	*	Load sequence sets
	*/
	ssPO = rssPO = NULL;
	if(!ReadInSeqsetI(csPO->inname,csPO->iform,TRUE,&ssPO,TRUE))
	{
		PROBLINE;
		printf("Didn't get sequences from %s\n",csPO->inname);
		return(FALSE);
	}
	csPO->seqs = ssPO;
	if(!quiet)
	{
		printf("# Have %d sequences loaded\n",ssPO->n);
	}
	max = SeqsetMaxLenI(ssPO);
	if(max >= MAX_CSL)
	{
		PROBLINE;
		printf("Max sequence length exceeded (%d)\n",max);
		printf("May only compare to length %d\n",MAX_CSL);
		return(FALSE);
	}
	/***
	*	Reference set?
	*/
	if(!NO_S(csPO->rinname))
	{
		if(!ReadInSeqsetI(csPO->rinname,csPO->iform,TRUE,&rssPO,TRUE))
		{
			printf("Didn't get reference seqs from %s\n",csPO->rinname);
			return(FALSE);
		}
		csPO->rseqs = rssPO;
		if(!quiet)
		{
			printf("# Have %d reference sequences loaded\n",rssPO->n);
		}
		max = SeqsetMinLenI(rssPO);
		if(max >= MAX_CSL)
		{
			PROBLINE;
			printf("Maximum seq length exceeded (ref set %d)\n",max);
			printf("May only compare to length %d\n",MAX_CSL);
			return(FALSE);
		}
	}
	/***
	*	Set pointer to actual sequences; Which one is target?
	*	If a ref seqset is real, use that, else use the first set
	*/
	csPO->nseqs = ssPO->n;
	if(rssPO)
	{	
		csPO->target = rssPO;
		csPO->ntargseqs = rssPO->n;
	}
	else
	{	
		csPO->target = ssPO;
		csPO->ntargseqs = ssPO->n;
	}
	DB_COMSEQ DB_PrI("<< GetCompSeqSeqsI ntargseqs=%d TRUE\n",csPO->ntargseqs);
	return(TRUE);
}
/***************************************************************************
*	Header info to dump with answers
*/
void WriteCompseqHeader(COMPSEQ *csPO,FILE *outPF)
{
	char bufS[DEF_BS];

	HAND_NFILE(outPF);
	fprintf(outPF,"# %s, %s %s %s\n",VERSION_S,BD_S,__DATE__,__TIME__);
	fprintf(outPF,"# Input........ %s\n",csPO->inname);
	if(!NO_S(csPO->rinname)) {
		fprintf(outPF,"# Compared to.. %s\n",csPO->rinname);
	}
	FillCompseqOutWhatString(csPO->owhat,bufS);
	fprintf(outPF,"# Reporting.... %s\n",bufS);
    /* xxx ? better place .... or just replace whole sham! */
    if(csPO->do_ham) {
	    fprintf(outPF,"# No alignment; Hamming distance\n");
    }
    else {
	    fprintf(outPF,"# Gap-free alignment max score\n");
    }
	fprintf(outPF,"# \n");
	DumpPpars(csPO->pp,outPF);	
	if(csPO->do_fmat) { 	
		WriteHeaderForFmat(csPO,outPF); 
		return;
	}
	if(csPO->ostat) {
		if(csPO->do_flag) {
			if(csPO->do_eraw) {
				fprintf(outPF,"# Name      Sequence\n");
            }
			else {
				fprintf(outPF,"# Flag  Name                    Score\n");
            }
		}
		else if( (csPO->do_self) || (csPO->owhat==CSO_RSER) ) {
			fprintf(outPF,"# Name         Score\n");
		}
		else {
			/* SHAM formatting should be dynamic? */
			fprintf(outPF,"# %-18s\t%s\t%s\t%s\t%s\n",
				"Name", "Min", "Ave", "Max", "N-Max");
		}
	}
    return;
}
/***************************************************************************/
void FillCompseqOutWhatString(int owhat,char *bufS)
{
	switch(owhat)
	{
		case CSO_LOOP: 	sprintf(bufS,"Hairpin loop stems"); break;
		case CSO_PSELF:	sprintf(bufS,"Parallel self complement"); break;
		case CSO_SELF:	sprintf(bufS,"Self complement"); break;
		case CSO_RSER:	sprintf(bufS,"Serially listed reference"); break;
		case CSO_FULL:	sprintf(bufS,"Full pair-wise comparison"); break;
		default:
			INIT_S(bufS);
	}
    return;
}
/***************************************************************************
*	Column labes for matrix dump
*/
void WriteHeaderForFmat(COMPSEQ *csPO,FILE *outPF)
{
	int i,j;
	char snameS[NSIZE];

	HAND_NFILE(outPF);
	fprintf(outPF,"Name");
	for(i=0;i<csPO->ntarg;i++)
	{
		if(i>=csPO->ntargseqs) {
			j = i - csPO->ntargseqs;
		}
		else {
			j = i;
		}
		FillSeqsetSeqNameI(csPO->target,j,snameS,NSIZE); 
		fprintf(outPF,"\t%s",snameS);
	}
	fprintf(outPF,"\n");
    return;
}
/***************************************************************************
*	Add space to hold scores
*/
int AddCompseqScoreSpaceI(COMPSEQ *csPO)
{
	DB_COMSEQ DB_PrI(">> AddCompseqScoreSpaceI sim=%d com=%d\n",
		csPO->do_sim,csPO->do_com);
	if(csPO->do_sim || csPO->do_com) { 	
		csPO->ntarg = csPO->ntargseqs;  
	}
	else {
		csPO->ntarg = 1;
	}
	DB_COMSEQ DB_PrI("+ ntarg = %d\n",csPO->ntarg);
	csPO->scores = (REAL *)ALLOC(csPO->ntarg,sizeof(REAL));
	csPO->offs = (int *)ALLOC(csPO->ntarg,sizeof(int));
	if( (!csPO->scores) || (!csPO->offs) ) {
		PROBLINE;
		printf("Failed to allocate for %d scores/offsets\n",csPO->ntarg);
		return(FALSE);
	}
	DB_COMSEQ DB_PrI("<< AddCompseqScoreSpaceI TRUE\n");
	return(TRUE);
}
/***************************************************************************/
int OpenCompseqOutputI(COMPSEQ *csPO)
{
	if(!NO_S(csPO->outname)) {
		if(!(csPO->out=OpenUFilePF(csPO->outname,"w",NULL))) {
			return(FALSE);
		}
	}
	HAND_NFILE(csPO->out);
	return(TRUE);
}
/***************************************************************************
*	Alloc and init compseq structure
*/
COMPSEQ *CreateCompseqPO()
{
	COMPSEQ *csPO;

	if(!(csPO = (COMPSEQ *)ALLOC(1,sizeof(COMPSEQ)))) {
		return(NULL);
	}
	csPO->ID = COMPSEQ_ID;
	csPO->pp = CreatePparsPO();
	csPO->seq = CreateSeqPO(0,NULL,NULL);
	if( (!csPO->pp) || (!csPO->seq) ) {
		CHECK_COMPSEQ(csPO);
		return(NULL);
	}
	InitCompseq(csPO);
	return(csPO);
}
/***************************************************************************
*	Free up structure and sub-structs
*/
int DestroyCompseqI(COMPSEQ *csPO)
{
	VALIDATE(csPO,COMPSEQ_ID);
	CHECK_FREE(csPO->offs);
	CHECK_FREE(csPO->scores);
	CHECK_PPARS(csPO->pp);
	CHECK_SEQ(csPO->seq);
	CHECK_SEQSET(csPO->seqs);
	CHECK_SEQSET(csPO->rseqs);
	CHECK_FILE(csPO->in);
	CHECK_NFILE(csPO->out,csPO->outname);
	CHECK_FREE(csPO);
	return(TRUE);
}
/****************************************************************************
*	Initialize global structure
*/
void InitCompseq(COMPSEQ *csPO)
{
	INIT_S(csPO->inname);
	INIT_S(csPO->rinname);
	INIT_S(csPO->outname);
	INIT_S(csPO->com);
	csPO->do_sa = FALSE;
	csPO->do_rm = FALSE;
	csPO->do_norm = FALSE;
	csPO->do_flag = csPO->do_eraw = FALSE;
	csPO->flo = csPO->fhi = BAD_R;
	csPO->do_sim = csPO->do_com = csPO->do_self = FALSE;
	csPO->do_loop = -1;
	csPO->fullmat = 1.0;
	csPO->partmat = 1.0;
}
/***************************************************************************
*	Check options 
*/
int CheckCompseqOptionsI(COMPSEQ *csPO)
{
	/***
	*	Anything to really do?
	*/
	if( (!csPO->do_sim) && (!csPO->do_com) && (!csPO->do_self) ) {
        Comp_seqUse();
		printf("\n");
		printf("   You must specify -sim, -comp, -self, -psel, or -loop!\n");
		printf("\n");
		return(FALSE);
	}
	if( (csPO->do_sim) && (csPO->do_com) ) {
		PROBLINE;
		printf("\n");
		printf("   Both -sim and -com are not allowed together\n");
		printf("\n");
		return(FALSE);
	}
	if( (csPO->owhat==CSO_RSER) && NO_S(csPO->rinname) ) {
		PROBLINE;
		printf("Serial comparision (-crs) only works with ref set (-rset)\n");
		return(FALSE);
	}
	if(!BAD_REAL(csPO->flo)) {
		if( csPO->flo > csPO->fhi ) {
			PROBLINE;
			printf("Bad flagging values: %4.2f to %4.2f\n",csPO->flo,csPO->fhi);
			return(FALSE);
		}
		csPO->do_flag = TRUE;
	}
	if( (csPO->do_fmat) && (csPO->do_self) )
	{
		PROBLINE;
		printf("-fmat option won't work with self-comparision\n");
		return(FALSE);
	}
	if( (csPO->do_eraw) && (!csPO->do_flag) )
	{
		PROBLINE;
		printf("-eraw option only works with flagging\n");
		return(FALSE);
	}
	return(TRUE);
}
/****************************************************************************
*	Process scores
*/
int ProcessScoresI(COMPSEQ *csPO,int cur,char *seqPC,int len)
{
	int j,nave;
	REAL maxR;

	/***
	*	Normalize by length?
	*/
	if(csPO->do_norm)
	{
		maxR = MaxSeqCompareScoreR(csPO->pp,seqPC,len);
		for(j=0;j<csPO->ntarg;j++)
		{
			if(csPO->scores[j] > 0.0)
			{
				csPO->scores[j] /= RNUM(maxR);
			}
		}
	}
	/***
	*	Tally all collected scores for current sequence to find 
	*	low, high, average
	*/
	csPO->max_sc = -TOO_BIG_R;
	csPO->min_sc = TOO_BIG_R;
	csPO->av_sc = 0.0;
	nave = 0;
	DB_COMSEQ DB_PrI("+ Tallying score array for [%d]\n",cur);
	for(j=0; j < csPO->ntarg; j++)
	{
		DB_COMSEQ DB_PrI("+ [%d]score[%d]=%6.2f ",cur,j,csPO->scores[j]);
		/***
		*	Don't include self-similarity (ident) in score tally
		*/
		if( (cur==j) && (csPO->do_sim) && (csPO->rseqs==NULL) )
		{
			DB_COMSEQ DB_PrI("ignore self\n");
			continue;
		}
		if(csPO->scores[j] < 0.0)
		{
			continue;
		}
		csPO->max_sc = MAX_NUM(csPO->scores[j],csPO->max_sc);
		csPO->min_sc = MIN_NUM(csPO->scores[j],csPO->min_sc);
		csPO->av_sc += csPO->scores[j];
		DB_COMSEQ 
			DB_PrI("min=%6.2f max=%6.2f av=%6.2f\n",csPO->min_sc,
				csPO->max_sc,csPO->av_sc);
		nave++;
	}
	if(nave>0)
	{ 	
		csPO->av_sc /= RNUM(nave); 
	}
	return(nave);
}
/****************************************************************************
*	Mark current seq if flagging is set
*/
int HandleSeqFlaggingI(COMPSEQ *csPO,SEQ *seqPO,FILE *outPF)
{
	int ok;
	char fnameS[NSIZE];

	if(!csPO->do_flag) 
	{
		return(FALSE);
	}
	/***
	*	Max score in range?
	*/
	HAND_NFILE(outPF);
	ok = TRUE;
	if( (csPO->max_sc<csPO->flo) || (csPO->max_sc>csPO->fhi) ) 
	{
		ok = FALSE;
	}
	if(csPO->do_not)
	{
		ok = !ok;
	}
	/***
	*	Extracting raw seq or simply flagging?
	*/
	FillSeqNameStringI(seqPO,fnameS,25);
	if(csPO->do_eraw)
	{
		if(ok)
		{
			fprintf(outPF,"# %-20s\t%3.2f\n",fnameS,csPO->max_sc);
			WriteSeq(seqPO,SEQFM_RAW,outPF);
		}
	}
	else
	{
		fprintf(outPF,"%d\t%-20s\t%3.2f\n",ok,fnameS,csPO->max_sc);
	}
	return(TRUE);
}
/****************************************************************************
*	Dump score matrix row if set
*	Return TRUE if any output, false otherwise
*/
int HandleMatrixOutI(COMPSEQ *csPO,char *fnameS,FILE *outPF)
{
	int j;
    REAL scR;

	HAND_NFILE(outPF);
	if(csPO->do_fmat)
	{
		fprintf(outPF,"%s",fnameS);
		for(j=0;j<csPO->ntarg;j++)
		{
		    scR = (csPO->scores[j] < 0.0) ? 0.0 : csPO->scores[j];
			fprintf(outPF,"\t%4.2f",scR);
		}
		fprintf(outPF,"\n");
		return(TRUE);
	}
	return(FALSE);
}
/****************************************************************************
*	Handle possibly multiple-outputs for current sequence
*/
int HandleFullScoreOut(COMPSEQ *csPO,int nsco,int cur,SEQ *seqPO)
{
	int j,nns,nok,ok;
	char fnameS[NSIZE];
	SEQ *sseqPO;

	/***
	*	If already no there are no scores, handle here and bail
	*/
	if(nsco<1)
	{
		HandleNoScoreOut(csPO,seqPO);
		return(TRUE);
	}
	/***
	*	Pass two over targets, considering each record score
	*/
	nns = nok = 0;
	for(j=0;j<csPO->ntarg;j++)
	{
		/***
		*	Ignore self similarity
		*	Set value to average in case computing stats at end
		*/
		if( (j==cur) && (csPO->do_sim) && (csPO->rseqs==NULL) )
		{
			csPO->scores[j] = csPO->av_sc;
			continue;
		}
		/***
		*	Check against any score range specified
		*/
		ok = FALSE;
		if(!BAD_REAL(csPO->flo))
		{
			if( (csPO->scores[j]>=csPO->flo) && (csPO->scores[j]<=csPO->fhi) )
			{
				ok++;	
			}
		}
		else if(csPO->scores[j]==csPO->max_sc)
		{	
			ok++;	
		}
		if(ok)
		{
			nok += ok;
		}
		/***
		*	If ok handle the match
		*/
		if( ok && (!csPO->ostat) )
		{ 		
			if(csPO->max_sc<0.0)
			{
				nns++;
			}
			if(j>=csPO->ntargseqs)
			{
				GetSeqsetSeqI(csPO->target,j - csPO->ntargseqs, &sseqPO);
			}
			else
			{
				GetSeqsetSeqI(csPO->target, j, &sseqPO);
			}
			if( (csPO->do_com) && (!nns) )	
			{
				HandleSeqPairOut(csPO,seqPO,sseqPO,j,TRUE,FALSE);
			}
			else if(!nns)
			{
				HandleSeqPairOut(csPO,seqPO,sseqPO,j,FALSE,FALSE);
			}
		}
	}
	/***
	*	End of comparisions for current record
	*/	
	if(csPO->ostat)
	{
		FillSeqsetSeqNameI(csPO->seqs,cur,fnameS,25);
        /*  SHAM should be dynamic formatting */
		fprintf(csPO->out,"%-20s\t%6.2f\t%6.2f\t%6.2f\t%d\n", 
            fnameS, csPO->min_sc, csPO->av_sc, csPO->max_sc, nok);
	}
	else if( (csPO->do_sa || csPO->do_rm) && (nok>0) )
	{
		fprintf(csPO->out,"\n");
	}
	return(TRUE);
}
