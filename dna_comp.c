/*
* dna_comp.c
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
#include "prim.h"
#include "dna.h"

#define DB_CONS if(DB[115])

/***************************************************************************
*   Calculates base composition and runs of like bases, including IUB 
*   If clean is true, any non-ACGT will fail
*/
int GetSeqCompositionI(char *seqS, int slen, int clean, SEQCOMP *scPO)
{
    int i,na,nc,ng,nt,tot;
    int ra,rc,rg,rt,rs,rw,rr,ry,rk,rm;
    int rowa,rowc,rowg,rowt,rows,roww,rowr,rowy,rowk,rowm;
    int prev, cur, din;

    DB_CONS DB_PrI(">> GetSeqCompositionI slen=%d\n",slen);
    /***
    *   Initialize seq base comp structure
    */
    InitSeqcomp(scPO);
    scPO->slen = slen;
    /***
    *   Initialize counts to 0
    */
    na = nc = ng = nt = 0;
    rowa = rowc = rowg = rowt = ra = rc = rg = rt = 0;
    rows = roww = rowr = rowy = rowm = rowk = 0;
    rs = rw = rr = ry = rm = rk = 0;
    prev = -1;
    /***
    *   Screen the bases
    */
    for(i=0;i<slen;i++)
    {
        switch(seqS[i])
        {
            /***
            *   A = W,R,M; not S,Y,K
            */
            case 'A': case 'a':
                na++;
                ra++;   
                rw++; rr++; rm++;
                rc = rg = rt = 0;
                rs = ry = rk = 0;
                cur = 0;
                break;
            /***
            *   C = S,Y,M; not W,R,K
            */
            case 'C': case 'c':
                nc++;
                rc++;
                rs++; ry++; rm++;
                ra = rg = rt = 0;
                rw = rr = rk = 0;
                cur = 1;
                break;
            /***
            *   G = S,R,K; not W,Y,M
            */
            case 'G': case 'g':
                ng++;
                rg++;
                rs++; rr++; rk++;
                ra = rc = rt = 0;
                rw = ry = rm = 0;
                cur = 2;
                break;
            /***
            *   T = W,Y,K; not S,R,M
            */
            case 'T': case 't':
                nt++;
                rt++;
                rw++; ry++; rk++;
                ra = rc = rg = 0;
                rs = rr = rm = 0;
                cur = 3;
                break;
            default:
                if(clean) {
                    DB_CONS DB_PrI("<< GetSeqCompositionI clean |%c| FALSE\n",
                        seqS[i]);
                    return(FALSE);
                }
                cur = -1;
        }
        if(ra>rowa)
        {   rowa = ra; }
        if(rc>rowc)
        {   rowc = rc; }
        if(rg>rowg)
        {   rowg = rg; }
        if(rt>rowt)
        {   rowt = rt; }
        if(rs>rows)
        {   rows = rs; }
        if(rw>roww)
        {   roww = rw; }
        if(rr>rowr)
        {   rowr = rr; }
        if(ry>rowy)
        {   rowy = ry; }
        if(rm>rowm)
        {   rowm = rm; }
        if(rk>rowk)
        {   rowk = rk; }
        /***
        *   Dinucleotide counts
        */
        if (prev >= 0) {
            if (cur >= 0) {
                din = (prev * 4) + cur;
                scPO->dinuc[din] += 1;
                scPO->n_dinuc += 1;
            }
        }
        prev = cur;
    }
    /***
    *   Total bases looked at
    */
    tot = na + nc + ng + nt;
    scPO->nbase = tot;
    /***
    *   Set values
    */
    DB_CONS 
    {
        DB_PrI("+ coma=%d comc=%d comg=%d comt=%d\n",na,nc,ng,nt);
        DB_PrI("+ rowa=%d rowc=%d rowg=%d rowt=%d\n",rowa,rowc,rowg,rowt);
        DB_PrI("+ rows=%d roww=%d rowr=%d rowy=%d rowm=%d rowk=%d\n",
            rows,roww,rowr,rowy,rowm,rowk);
    }
    scPO->na = na;
    scPO->nc = nc;
    scPO->ng = ng;
    scPO->nt = nt;
    if(tot > 0) {
        scPO->fa = RNUM(na)/RNUM(tot);
        scPO->fc = RNUM(nc)/RNUM(tot);
        scPO->fg = RNUM(ng)/RNUM(tot);
        scPO->ft = RNUM(nt)/RNUM(tot);
    }
    scPO->ra = rowa;
    scPO->rc = rowc;
    scPO->rg = rowg;
    scPO->rt = rowt;
    scPO->rs = rows;
    scPO->rw = roww;
    scPO->rr = rowr;
    scPO->ry = rowy;
    scPO->rm = rowm;
    scPO->rk = rowk;
    DB_CONS DB_PrI("<< GetSeqCompositionI TRUE\n");
    return(TRUE);
}
