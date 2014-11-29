/*
* blasthit.c
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
#include "prim.h"
#include "blastout.h"

#define DB_BIO  if(DB[133])

/****************************************************************************
*	Process raw hit list and possibly cull some records
*/
void ProcHitList(BLASTOUT *bPO, BLASTANS *aPO)
{
	int i,n,ok,num,max,den,alen;
	REAL fracR;

    DB_BIO DB_PrI(">> ProcHitList\n");
	n = max = 0;
	for(i=0;i<aPO->nhits;i++) {
		if( ((i+1) < bPO->firsth) || ((i+1) > bPO->lasth) ) {
			continue;
		}
		num = HitMatchCountI(aPO,i,bPO->firstb,bPO->lastb,bPO->rre,
			bPO->do_con,&den);
		max = MAX_NUM(num,max);
        /***
        *   Qualify by max match 
        */
        ok = TRUE;
		if(num < bPO->mid) {
			ok = FALSE;
		}
		alen = AdjustHitLenI(bPO,aPO,i,den);
		fracR = RNUM(num)/RNUM(alen);
		if(fracR < bPO->mif) {
			ok = FALSE;
		}
        ok = bPO->do_mnot ? (!ok) : ok;
        if (!ok) {
			continue;
        }
		/***
		*	Current record is cool; shrink list if any have been skipped
		*/
		if(i != n) {
			strcpy(&aPO->scos[n*BLBSIZE],&aPO->scos[i*BLBSIZE]);
			strcpy(&aPO->idens[n*BLBSIZE],&aPO->idens[i*BLBSIZE]);
			strcpy(&aPO->hits[n*BLBSIZE],&aPO->hits[i*BLBSIZE]);
			strcpy(&aPO->qseqs[n*BLBSIZE],&aPO->qseqs[i*BLBSIZE]);
			strcpy(&aPO->sseqs[n*BLBSIZE],&aPO->sseqs[i*BLBSIZE]);
		}
		aPO->hnums[n] = num;
		aPO->hlens[n] = alen;
		aPO->hfmb[n] = fracR;
		n++;
	}
	aPO->nhits = n;
	aPO->maxhit = max;
    DB_BIO DB_PrI("<< ProcHitList\n");
}
/****************************************************************************
*	Merge hits from two BLASTANS into a third
*/
int MergeHitListsI(BLASTOUT *bPO, BLASTANS *a1PO,BLASTANS *a2PO,BLASTANS *a3PO)
{
	int i,n,a,a1,a2,num,max,den,alen;
	REAL fracR;
	BLASTANS *aPO;

	if(!SameBlastQueryI(a1PO,a2PO))
	{
		return(FALSE);
	}
	/***
	*	Merge hits into a3PO
	*/
	n = max = a = a1 = a2 = 0;
	aPO = NULL;
	for(i=0;i<a3PO->asize;i++)
	{
		/***
		*	Figure out which source to copy from
		*/
		if( (a1>=a1PO->nhits) && (a2>=a2PO->nhits) )
		{
			break;
		}
		else if( (a1<a1PO->nhits) && (a2<a2PO->nhits) )
		{
			if(ODD_NUM(n))
			{
				aPO = a2PO;
				a = a2++;
			}
			else
			{
				aPO = a1PO;
				a = a1++;
			}
		}
		else if(a1<a1PO->nhits)
		{
			aPO = a1PO;
			a = a1++;
		}
		else if(a2<a2PO->nhits)
		{
			aPO = a2PO;
			a = a2++;
		}
		num = HitMatchCountI(aPO,a,bPO->firstb,bPO->lastb,bPO->rre,
			bPO->do_con,&den);
		max = MAX_NUM(num,max);
		if(num < bPO->mid)
		{
			continue;
		}
		alen = AdjustHitLenI(bPO,aPO,a,den);
		fracR = RNUM(num)/RNUM(alen);
		if(fracR < bPO->mif)
		{
			continue;
		}
		/***
		*	Copy current record into answer set
		*/
		strcpy(&a3PO->scos[n*BLBSIZE],&aPO->scos[a*BLBSIZE]);
		strcpy(&a3PO->idens[n*BLBSIZE],&aPO->idens[a*BLBSIZE]);
		strcpy(&a3PO->hits[n*BLBSIZE],&aPO->hits[a*BLBSIZE]);
		strcpy(&a3PO->qseqs[n*BLBSIZE],&aPO->qseqs[a*BLBSIZE]);
		strcpy(&a3PO->sseqs[n*BLBSIZE],&aPO->sseqs[a*BLBSIZE]);
		a3PO->hnums[n] = num;
		a3PO->hlens[n] = alen;
		a3PO->hfmb[n] = fracR;
		a3PO->qhsc[n] = aPO->qhsc[a];
		a3PO->qhec[n] = aPO->qhec[a];
		a3PO->shsc[n] = aPO->shsc[a];
		a3PO->shec[n] = aPO->shec[a];
		n++;
	}
	a3PO->nhits = n;
	a3PO->maxhit = max;
	strcpy(a3PO->query,a1PO->query);
	a3PO->qlen = a1PO->qlen;
	return(TRUE);
}
/****************************************************************************
*	Fill histogram with matching base counts
*/
int FillHitHistI(BLASTANS *aPO,int first,int last,int rre,int con)
{
	int i,num,max;

	for(i=0;i<HISTDIM;i++)
	{
		aPO->ihist[i] = 0;
	}
	max = 0;
/**
printf("nhits=%d\n",aPO->nhits);
*/
	for(i=0;i<aPO->nhits;i++)
	{
		num = HitMatchCountI(aPO,i,first,last,rre,con,NULL);
		max = MAX_NUM(num,max);
		if(num >= HISTDIM)
			continue;
		aPO->ihist[num] += 1;
/**
printf("i=%d num=%d ihis=%d max=%d\n",i,num,aPO->ihist[num],max);
**/
	}
	return(max);
}
/****************************************************************************
*	Fill histogram with matching base counts
*/
void IntegrateHist(BLASTANS *aPO)
{
	int i,num;

	num = 0;
	for(i=HISTDIM-1;i>0;i--)
	{
		num += aPO->ihist[i];
		aPO->ihist[i] = num;
	}
}
/****************************************************************************
*	Returns the count of (qualified) matching bases 
*	Also sets length of (qualified) alignment if passed third arg
*/
int HitMatchCountI(BLASTANS *aPO, int hit, int first, int last, int rre,
	int con, int *alenPI)
{
	int qs,qe,gap,ok,b,j,c,alen,max,num;
	char *qPC,*sPC;

	qs = aPO->qhsc[hit]; 
	qe = aPO->qhec[hit];
/* xxx 
printf("xxx hit=%d\tqs=%d\tqe=%d\n",hit,qs,qe);
*/
	j = gap = num = c = max = alen = 0;
	/***
	*	Scan alignment until out of chars
	*/
	qPC = &aPO->qseqs[hit * BLBSIZE];
	sPC = &aPO->sseqs[hit * BLBSIZE];
	while(isgraph(INT(qPC[j])))
	{
		/***
		*	Restricted base range?
		*/
		ok = TRUE;
		if(qs<qe)
		{
			b = qs+j-1-gap;
		}
		else
		{
			b = qs-j-1+gap;
		}
		if(first>0)
		{
			if(rre)
			{
				if( ((aPO->qlen-b)<first) || 
					((aPO->qlen-b)>last) )
					ok = FALSE;
			}
			else
			{
				if( ((b+1)<first) || (b>=last) )
					ok = FALSE;
			}
		}
		if(!ok)
		{
			j++;
			continue;
		}
		alen++;
		/***
		*	Missmatch
		* SHAM; Could check if degenerate IUB codes match?
		*/
		if(qPC[j] != sPC[j])
		{
			if(qPC[j] == '-')
			{
				gap++;
			}
			max = MAX_NUM(c,max);
			c = 0;
			j++;
			continue;
		}
		/***
		*	Matching case, up counts
		*/
		c++;
		num++;
		j++;
	}
	/***
	*	Final counts and what to return?
	*/
	max = MAX_NUM(c,max);
	if(con)
	{
		num = max;
	}
	if(alenPI)
	{
		*alenPI = alen;
	}
	return(num);
}
/**************************************************************************/
int SameBlastQueryI(BLASTANS *fPO, BLASTANS *sPO)
{
	if(strcmp(fPO->query,sPO->query))
	{
		PROBLINE;
		printf("Different query names\n");
		printf(" From %s\n",fPO->input);
		printf(" Name %s\n",fPO->query);
		printf(" From %s\n",sPO->input);
		printf(" Name %s\n",sPO->query);
		return(FALSE);
	}
	return(TRUE);
}
