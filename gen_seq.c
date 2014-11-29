/*
* gen_seq.c
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
#include <math.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "dna_cons.h"
#include "gen_seq.h"

#define DB_GS if(DB[200])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit(AllDoneI(GenSeqI(argc,argv),"gen_seq")); }
/**************************************************************************/
void GenSeqUse(void)
{
	VersionSplash(NULL,VERSION_S,"#  ",TRUE);
	printf("Useage: [-filt <infile> or '-' for stdin] [...options]\n");
	printf("   -par XXX   Parameters (intrinsic constraints) from file XXX\n");
	printf("   -dpar      Dump parameter constraints used\n");
	printf("   -filt XXX  Read seqs from XXX and filter on parameters\n");
    printf("   -iraw      Treat input sequences as 'raw' format\n");
    printf("   -ifas      Treat input sequences as fasta format\n");
	printf("   -ran       Random sequences\n");
	printf("   -base XXX  Random base percentages; A,C,G,T in format #,#,#,#\n");
	printf("   -maxt #    Random seq maximum tries before quitting\n");
	printf("   -num #     Random seq number of seqs to generate\n");
	printf("   -enu       Exhaustive enumaration of all seqs (to len)\n");
	printf("   -len #     Generate sequences of length #\n");
	printf("   -not       Invert constraints; BAD = ok, GOOD = not ok\n");
	printf("   -both      Check both strands (enumeration algorithm)\n");
	printf("   -out XXX   Output to file XXX\n");
	printf("   -oraw      Output sequences int 'raw' format\n");
	printf("   -ofas      Output sequences int fasta format\n");
	printf("   -dall      Dump all sequences generated; ok or not\n");
	printf("   -dgb       Dump good/bad status with seqs\n");
	printf("   -bname XXX Set seq base name to XXX\n");
	printf("   -stat      Report stats of generation\n");
	printf("   -ostat     Only report stats of generation; nothing else\n");
	printf("   -seed #    Set random seed to #\n");
	printf("   -quiet     Quiet mode\n");
}
/**************************************************************************
*	Main gen_seq program
*/
int GenSeqI(int argc, char **argv)
{
	int nt,ng,slen,mtry,dump,quiet,ok;
	char bnameS[NSIZE];
	GENSEQ *gsPO;
	SEQ *seqPO;

	INIT_S(bnameS); 
	dump = quiet = FALSE;
	mtry = TOO_BIG;
 	gsPO = CreateGenseqPO();
	if(!gsPO)
	{
		printf("Failed to allocate genseq structure for settings\n");
		ABORTLINE;
		return(FALSE);
	}
	if(!ParseArgsI(argc,argv,
        "-par S -dpar B -dall B -stat B -oraw B\
		-ofas B -out S -bname S\
		-len I -both B -enu B -seed I -num I -dgb B\
		-ran B -ostat B -filt S -not B -maxt I -quiet B\
		-base S -iraw B -ifas B",
		gsPO->parname, &dump, &gsPO->dump_all, &gsPO->do_stat,&gsPO->out_raw,
		&gsPO->out_fasta, gsPO->outfile, bnameS,
		&gsPO->len, &gsPO->do_both, &gsPO->do_enu, &gsPO->rseed, &gsPO->num, 
		&gsPO->dump_ok, &gsPO->do_ran, &gsPO->do_ostat, 
		gsPO->filtname, &gsPO->do_not, &mtry, &quiet,
		&gsPO->baseper, 
		&gsPO->iraw, &gsPO->ifas, 
        (int *)NULL))
    {
        GenSeqUse();
		CHECK_GENSEQ(gsPO);
		return(FALSE);
	}
	if(!CheckGenSeqOptionsI(gsPO))
	{
		CHECK_GENSEQ(gsPO);
		return(FALSE);
	}
	if(!SetUpGenSeqI(gsPO))
	{
		ABORTLINE;
		CHECK_GENSEQ(gsPO);
		return(FALSE);
	}
	/***
	*	More convienent pointer and update size to fit
	*/
	seqPO = gsPO->seq;
	AdjustSeqSizeI(seqPO,gsPO->len,TRUE);
	/***
	*	Sequence base name?
	*/
	if(!NO_S(bnameS)) {
		sprintf(gsPO->bname,"%s_",bnameS);
	}
	/***
	*	Dump seqs / parameters?
	*/
	if(dump || gsPO->out_raw) {
		WriteGenseqHeader(gsPO,gsPO->out);
	}
	/***
	*	Filter preexisting set?
	*/
	ng = nt = 0;
	if(gsPO->filt) {
		while(TRUE)
		{
			/***
			*	Parse sequence; FALSE = done
			*/
			ok = ParseSeqI(gsPO->filt,gsPO->fform,SCLEAN_HI,TRUE,seqPO);
			if(ok==FALSE) {
				break;
			}
			if(ok!=TRUE) {
				if(!gsPO->igprob) {
					break;
				}
				continue;
			}
			slen = GetSeqLenI(seqPO);
			FillSeqNameStringI(seqPO,bnameS,NSIZE-1);
			gsPO->ok = FullSeqInConsOkI(seqPO->seq,slen,TRUE,gsPO->cons);
			if(gsPO->do_not) {
				gsPO->ok = !gsPO->ok;
			}
			if(gsPO->dump_all) {
				HandleCompleteSeqI(seqPO->seq,slen,bnameS,gsPO->cons,gsPO);
				if(nt >= gsPO->num) {
					break;
				}
			}
			else {
				if(gsPO->ok) {
					HandleCompleteSeqI(seqPO->seq,slen,bnameS,gsPO->cons,gsPO);
					ng++;
				}
			}
			nt++;
		}
	}
	/***
	*	Random sequence generation?
	*/
	else if(gsPO->do_ran) {
		if(!quiet) {
			PrintI("# Attempting to generate %d random sequences\n",gsPO->num);
        }
		while(ng < gsPO->num)
		{
			RandomDNASeqI(seqPO->seq,gsPO->len,gsPO->basepnum);
			gsPO->ok = FullSeqInConsOkI(seqPO->seq,gsPO->len,TRUE,gsPO->cons);
			if(gsPO->do_not) { 	
                gsPO->ok = !gsPO->ok; 
            }
			if(gsPO->ok || gsPO->dump_all) {
				gsPO->n_good += 1;
				HandleCompleteSeqI(seqPO->seq,gsPO->len,NULL,gsPO->cons,gsPO);
				if(gsPO->ok)
					ng++;
			}
			nt++;
			if((nt%GENSEQ_UDF)==0) {	
				PrintI("# %d of %d trys (%5.4f%%)\n",ng,nt,
				PERCENT_R(ng,nt));	
			}
			if(nt > mtry) {
				break;
			}
		}
	}
	/***
	*	Top call to recursive seq enumeration function; 
	*/
	else if(gsPO->do_enu) {
		/***
		*	If both strands, simply call to start at all 4 bases
		*	Otherwise start here with A and C only; G and T are compliments
		*/
		gsPO->n_max = NumSkippedD(0,gsPO->len);
		if(gsPO->do_both) {
			printf("# Evaluating %1.4E sequences (both strands)\n",gsPO->n_max);
			ExtendSeqI(seqPO->seq,0,gsPO->cons,gsPO);
		}
		else {
			printf("# Evaluating %1.4E sequences (one strand)\n",
				gsPO->n_max/2.0);
			/***
			*	Recur with A 
			*/
			seqPO->seq[0] = 'A';
			ExtendSeqI(seqPO->seq,1,gsPO->cons,gsPO);
			/***
			*	Recur with C 
			*/
			seqPO->seq[0] = 'C';
			ExtendSeqI(seqPO->seq,1,gsPO->cons,gsPO);
		}
	}
	/***
	*	Summary?
	*/
	if( gsPO->do_ostat || gsPO->do_stat)
	{
		if(gsPO->out_raw || gsPO->out_fasta) {
			printf("# ");
        }
		if(gsPO->do_enu) {
			printf("%1.6E of %1.6E seqs %d long OK (%1.4f%%)\n",gsPO->n_good,
				gsPO->n_max, gsPO->len, PERCENT_R(gsPO->n_good,gsPO->n_max));
        }
		else {
			printf("%d of %d seqs OK (%1.4f%%)\n",ng,nt,PERCENT_R(ng,nt));
        }
	}
	CHECK_GENSEQ(gsPO);
	return(TRUE);
}
/***************************************************************************
*	Allocate GENSEQ structure
*/
GENSEQ *CreateGenseqPO()
{
	GENSEQ *gsPO;

	if(!(gsPO = (GENSEQ *)ALLOC(1,sizeof(GENSEQ)))) {
		return(NULL);
	}
	gsPO->ID = GENSEQ_ID;
	gsPO->cons = CreateInConsPO();
	gsPO->seq = CreateSeqPO(DEF_LEN,NULL,NULL);
	if( (!gsPO->cons) || (!gsPO->seq) ) {
		CHECK_GENSEQ(gsPO);
		return(NULL);
	}
	InitGenseq(gsPO);
	return(gsPO);
}
/***************************************************************************
*	Free GENSEQ struct
*/
int DestroyGenseqI(GENSEQ *gsPO)
{
	VALIDATE(gsPO,GENSEQ_ID);
	CHECK_SEQ(gsPO->seq);
	CHECK_IN_CONS(gsPO->cons);
	CHECK_FILE(gsPO->filt);
	CHECK_NFILE(gsPO->out,gsPO->outfile);
	CHECK_FREE(gsPO);
	return(TRUE);
}
/*****************************************************************************
*	Initialize GENSEQ data structure vars
*/
void InitGenseq(GENSEQ *gsPO)
{
	VALIDATE(gsPO,GENSEQ_ID);
	INIT_S(gsPO->parname);
	INIT_S(gsPO->baseper);
	InitArrayI(gsPO->basepnum,IS_INT,0,4,25);
	INIT_S(gsPO->outfile);
	gsPO->out = NULL;
	INIT_S(gsPO->filtname);
	gsPO->filt = NULL;
	gsPO->fform = BOGUS;
    gsPO->iraw = gsPO->ifas = FALSE;
	sprintf(gsPO->bname,"%s_",DEF_BNAME_S);
	gsPO->n_tot = 0;
	gsPO->n_max = 0;
	gsPO->n_good = 0;
	gsPO->nexttime = GENSEQ_UDF;
	gsPO->verb = FALSE;
	gsPO->do_both = FALSE;
	gsPO->do_stat = FALSE;
	gsPO->do_ostat = FALSE;
	gsPO->do_not = FALSE;
	gsPO->dump_all = FALSE;
	gsPO->dump_ok = FALSE;
	gsPO->out_raw = FALSE;
	gsPO->out_fasta = FALSE;
	gsPO->rseed = BOGUS;
	gsPO->do_ran = FALSE;
	gsPO->do_enu = FALSE;
	gsPO->len = DEF_LEN;
	gsPO->num = DEF_NUM;
}
/*****************************************************************************
*	Allocate and set objects, open files
*/
int SetUpGenSeqI(GENSEQ *gsPO)
{
	/***
	*	Filter set input file?
	*/
	if(!NO_S(gsPO->filtname))
	{
		if( ! (gsPO->filt=OpenUFilePF(gsPO->filtname,"r",NULL)) ) {
			printf("Can't open input to filter!\n");
			return(FALSE);
		}
        gsPO->fform = FigureSeqFileTypeI(gsPO->iraw,gsPO->ifas,gsPO->filtname,TRUE);
        if(!gsPO->fform) {
            printf("Problem with input seq(s)\n");
            return(FALSE);
        }
	}
	/***
	*	Load input parameters 
	*/
	if(!NO_S(gsPO->parname)) {
		if(!LoadInConsParsI(gsPO->parname,gsPO->cons)) {
			PROBLINE;
			printf("Problem parsing parameter file\n");
			return(FALSE);
		}
		printf("# Parameters loaded from %s\n",gsPO->parname);
	}
	if(gsPO->filt) {
		SetInConsMaxlenI(gsPO->cons,IN_CONS_MAX-1);
	}
	else {
		SetInConsMaxlenI(gsPO->cons,gsPO->len);
	}
	/***
	*	Set up run-time values then make sure they can yield seqs
	*/
	if(!PrepareInConsI(gsPO->cons)) {
		PROBLINE;
		printf("Problem preparing constraints\n");
		return(FALSE);
	}
	if(!HandleBasePerSetupI(gsPO)) {
		PROBLINE;
		printf("Problem with base percentiles\n");
		return(FALSE);
	}
	if(!ConsistInConsI(gsPO->cons)) {
		PROBLINE;
		printf("Inconsistent / impossible constraints\n");
		return(FALSE);
	}
	/***
	*	Random seed
	*/
	Srand(gsPO->rseed);
	/***
	*	Output file?
	*/
	gsPO->out = NULL;
	if(!NO_S(gsPO->outfile)) {
		if(!(gsPO->out = OpenUFilePF(gsPO->outfile,"w",NULL))) {
			CHECK_GENSEQ(gsPO);
			ABORTLINE;
			return(FALSE);
		}
	}
	HAND_NFILE(gsPO->out);
	return(TRUE);
}
/**************************************************************************
*	Have to have some options or not worth running
*/
int CheckGenSeqOptionsI(GENSEQ *gsPO)
{
	VALIDATE(gsPO,GENSEQ_ID);
/*  Old sham, pre v1.0
	if( (!gsPO->do_ran) && (!gsPO->do_enu) && (NO_S(gsPO->filtname)) )
	{
		GenSeqUse();
		PROBLINE;
		printf("You must specify -ran, -enu, or -filt!\n"); 
		return(FALSE);
	}
*/
    if( (!gsPO->do_enu) && (NO_S(gsPO->filtname)) ) {
	    gsPO->do_ran = TRUE;
    }
	/***
	*	What kind of output?
	*/
	if( (!gsPO->out_raw) && (!gsPO->out_fasta) ) {
		gsPO->out_raw = TRUE;
	}
	if( gsPO->do_ostat ) {
		gsPO->out_raw = gsPO->out_fasta = FALSE;
	}
	return(TRUE);
}
/************************************************************************
*	Parse base perencile into array percents (if exists)
*/
int HandleBasePerSetupI(GENSEQ *gsPO)
{
	int a,c,g,t;

	if(NO_S(gsPO->baseper)) {
		return(TRUE);
	}
	a = c = g = t = BOGUS;
	sscanf(gsPO->baseper, "%d,%d,%d,%d",&a,&c,&g,&t);
	if(IS_BOG(a) || IS_BOG(c) || IS_BOG(g) || IS_BOG(t) )
	{
		PROBLINE;
		printf("Failed to parse 4 percentiles from |%s|\n",gsPO->baseper);
		printf("Format expected to be #,#,#,# as e.g. \"10,20,30,40\"\n");
		return(FALSE);
	}
	if( (a+c+g+t) != 100)
	{
		PROBLINE;
		printf("Percentiles for ACGT don't add to 100\n");
		printf("Base percentiles:\tA=%d\tC=%d\tG=%d\tT=%d\t=%d%%\n",a,c,g,t,a+c+g+t);
		return(FALSE);
	}
	gsPO->basepnum[0] = a;
	gsPO->basepnum[1] = c;
	gsPO->basepnum[2] = g;
	gsPO->basepnum[3] = t;
	return(TRUE);	
}
/************************************************************************/
void WriteGenseqHeader(GENSEQ *gsPO,FILE *outPF)
{
	int ran;
    char sS[DEF_BS];

	VALIDATE(gsPO,GENSEQ_ID);
	HAND_NFILE(outPF);
	fprintf(outPF,"# %s\n",gsPO->outfile);
	fprintf(outPF,"# Sequences output from %s\n",VERSION_S);
	fprintf(outPF,"#    %s %s %s\n",BD_S,__DATE__,__TIME__);
	fprintf(outPF,"#    %s\n",RTK_S);
	TimeStamp("# ",outPF);
	fprintf(outPF,"#\n");
	ran = FALSE;
	if(gsPO->filt)
	{
		fprintf(outPF,"# Filtering sequences from %s (%d)\n",
			gsPO->filtname,gsPO->num);
	}
	else if(gsPO->do_ran)
	{
		fprintf(outPF,"# Random sequences (len %d)\n",gsPO->len);
	    if(! NO_S(gsPO->baseper))
        {
		    fprintf(outPF,"# Base percentiles (A,C,G,T): %s\n",gsPO->baseper);
        }
		ran++;
	}
	else if(gsPO->do_enu)
	{
		fprintf(outPF,"# Enumerating sequences (len %d)\n",gsPO->len);
	}
	if(ran)
	{
        FillRandSeedString(gsPO->rseed,sS);
        fprintf(outPF,"# Random seed: %s\n",sS);
	}
	DumpInCons(gsPO->cons,outPF);
	if(gsPO->do_not) {
		fprintf(outPF,"# NOT = TRUE so constraint criteria are inverted\n");
    }
}
/****************************************************************************
*	Recursive sequence enumeration function
*/
int ExtendSeqI(char *seqS, int len, IN_CONS *iconPO, GENSEQ *gsPO)
{
	char compS[IN_CONS_MAX+1];
	int n;

	/***
	*	If not full length check if ok; if no seq at all, ok to extend
	*/
	if((len) && (len<iconPO->maxlen))
	{
		gsPO->ok = SeqInConsOkI(seqS,len,FALSE,iconPO);
	}
	else if(!len)
	{
		gsPO->ok = TRUE;
	}
	/***
	*	If long enough already, handle here
	*	If looking at both sequences, don't consider the complementary seq;
	*		otherwise look at both seq and compliment
	*/
	if(len == iconPO->maxlen)
	{
		/***
		*	Provide feedback?
		*/
		if(gsPO->n_tot >= gsPO->nexttime)
		{
			printf("  %1.4E seqs\n",gsPO->n_tot);
			fflush(stdout);
			gsPO->nexttime += GENSEQ_UDF;
		}
		/***
		*	Check on finished seq constraints
		*/
		n = 0;
		gsPO->ok = FullSeqInConsOkI(seqS,len,FALSE,iconPO);
		HandleCompleteSeqI(seqS,len,NULL,iconPO,gsPO);
		gsPO->n_tot += 1;
		gsPO->n_good += gsPO->ok;
		n += gsPO->ok;
		/***
		*	If not explicitly enumerating both strands, check compliment here
		*/
		if(!gsPO->do_both)
		{
			InverseDNASeqI(seqS,len,compS);
			gsPO->ok = SeqInConsOkI(compS,len,FALSE,iconPO);
			HandleCompleteSeqI(compS,len,NULL,iconPO,gsPO);
			gsPO->n_tot += 1;
			gsPO->n_good += gsPO->ok;
			n += gsPO->ok;
		}
		return(n);
	}
	/***
	*	If not keeping everything and sequence is bad, bail early
	*/
	if((!gsPO->ok)&&(!gsPO->dump_all))
	{
		return(FALSE);
	}
	/***
	*	Recur extended with each letter
	*/
	seqS[len] = 'A'; 
	ExtendSeqI(seqS,len+1,iconPO,gsPO);
	seqS[len] = 'C'; 
	ExtendSeqI(seqS,len+1,iconPO,gsPO);
	seqS[len] = 'G'; 
	ExtendSeqI(seqS,len+1,iconPO,gsPO);
	seqS[len] = 'T'; 
	ExtendSeqI(seqS,len+1,iconPO,gsPO);
	return(TRUE);
}
/****************************************************************************
*	How many seqs skipped if enumeration tree truncated early?
*/
double NumSkippedD(int len,int max)
{
	int i,l;
	double nD;

	l = max-len;
	nD=1.0;
	for(i=0;i<l;i++)
		nD*=4.0;
	return(nD);
}
/****************************************************************************
*	Sequence is completely enumerated so handle it based on output options
*/
int HandleCompleteSeqI(char *seqS,int len,char *nS,IN_CONS *iconPO,
	GENSEQ *gsPO)
{
	/***
	*	If only stats, done
	*/
	if(gsPO->do_ostat)
	{	
		return(TRUE);	
	}
	/***
	*	Dump everything else only if ok
	*/
	if(gsPO->dump_all || gsPO->ok)
	{
		if(gsPO->out_fasta)
		{
			if(nS)
				fprintf(gsPO->out,"> %s",nS);
			else
				fprintf(gsPO->out,"> %s%06.0f",gsPO->bname,gsPO->n_good);
			if(gsPO->dump_ok)
			{
				if(gsPO->ok)
					fprintf(gsPO->out," GOOD");
				else
					fprintf(gsPO->out," BAD ");
			}
			fprintf(gsPO->out,"\n");
			fprintf(gsPO->out,"\t");
			PrintString(seqS,len,gsPO->out);
			fprintf(gsPO->out,"\n");
		}
		else if(gsPO->out_raw)
		{
			if(nS)
				fprintf(gsPO->out,"%-12s",nS);
			else
				fprintf(gsPO->out,"%s%06.0f",gsPO->bname,gsPO->n_good);
			if(gsPO->dump_ok)
			{
				if(gsPO->ok)
					fprintf(gsPO->out," GOOD");
				else
					fprintf(gsPO->out," BAD ");
			}
			fprintf(gsPO->out," ");
			PrintString(seqS,len,gsPO->out);
			fprintf(gsPO->out,"\n");
		}
	}
	return(TRUE);
}
