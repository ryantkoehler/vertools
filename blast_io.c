/*
* blast_io.c
*
* Copyright 2015 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
*   Load current record hit info 
*/
int LoadBlastoutRecordI(BLASTANS *aPO,int warn)
{
    int ok;

    InitBlastans(aPO,FALSE);
    switch(aPO->itype) {
        case BOT_NCBI:  
        case BOT_WU:    
        case BOT_SMWU:  
            ok = LoadBlastRecordI(aPO,aPO->in,warn);
            break;
        case BOT_BO:    
            ok = LoadBdciRecordI(aPO,aPO->in,warn);
            break;
        default:
            ok = FALSE;
    }
    return(ok);
}
/****************************************************************************
*   Load info from a blast search output for one record into data struct
*   Returns TRUE if seems to parse ok;
*           FALSE if problems
*/
int LoadBlastRecordI(BLASTANS *aPO,FILE *inPF,int warn)
{
    int in,hit,n,fpos,ok;
    char bufS[BBUFF_SIZE];

    DB_BIO DB_PrI(">> LoadBlastRecordI\n");
    in = hit = n = 0;
    /***   
    *   Pass whatever header to get to start of query
    */ 
    while(fgets(bufS,BLINEGRAB,inPF)) {
        if(EQSTRING(bufS,"Query=",6)) {
/*
            if( ! LoadBlastQueryLineData(aPO, bufS, inPF) ) {
                return(FALSE);
            }
*/
            ok = LoadBlastQueryLineData(aPO, bufS, inPF);
            if( ok != TRUE ) {
                if( (warn) && (ok < 0) ) {
                    PROBLINE;
                    printf("Failed to load header block for query\n|%s|\n",bufS);
                }
                return(FALSE);
            }
            in++;
            break;
        }
    }
    if(!in) {
        DB_BIO DB_PrI("<< LoadBlastRecordI no Query= found = FALSE\n");
        return (FALSE);
    }
    /***
    *   Now each hit until next record block... at which point we back up and bail
    */
    fpos = ftell(inPF);
    while(fgets(bufS,BLINEGRAB,inPF)) {
        /***
        *   Next record block so back up and break out
        */
        if(EQSTRING(bufS,"Query=",6)) {
            DB_BIO DB_PrI("+ Next record so back up and bail\n");
            fseek(inPF,fpos,0);
            break;
        }
        /***
        *   Load hit (subject) block (follows > line)
        */
        if(bufS[0] == '>') {
            n = AddBlastHitSetI(bufS,inPF,aPO,hit,warn);
            if(n<1) {
                DB_BIO DB_PrI("+ AddHit =%d break\n",n);
                break;
            }
            hit += n;
        }
    }
    aPO->nhits = hit;
    DB_BIO DB_PrI("<< LoadBlastRecordI hit=%d TRUE\n",hit);
    return(TRUE);
}
/****************************************************************************
*   Load info for single query data from passed input line (and maybe next)
*   Update for NCBI case with "Length=" 7/5/14 RTK
*/
int LoadBlastQueryLineData(BLASTANS *aPO, char *inbufS, FILE *inPF)
{
    char bufS[BBUFF_SIZE];

    DB_BIO DB_PrI(">> LoadBlastQueryLineData\n");
    sscanf(inbufS,"%*s %s",aPO->query);
    if( aPO->itype == BOT_SMWU ) {
        /***
        *   Query= <query> (<len> letters)
        */
        ReplaceChars('(',inbufS,' ',bufS);
        sscanf(bufS,"%*s %*s %d",&aPO->qlen);
    }
    else if( aPO->itype == BOT_NCBI) {
        /***
        *   Query= <query>
        *                           <-- Next line 1
        *   Length=<len>            <-- Next line 2
        */
        if(!fgets(bufS,BLINEGRAB,inPF) ){
            return(BOGUS);
        }
        if(!fgets(bufS,BLINEGRAB,inPF) ){
            return(BOGUS);
        }
        if(!EQSTRING(bufS,"Length=",7)) {
            return(BOGUS);
        }
        ReplaceChars('=',bufS,' ',bufS);
        sscanf(bufS,"%*s %d",&aPO->qlen);
    }
    else {
        /***
        *   Query= <query>
        *       (<len> letters)     <-- Need to get next line
        */
        if(!fgets(bufS,BLINEGRAB,inPF) ){
            return(BOGUS);
        }
        ReplaceChars('(',bufS,' ',bufS);
        sscanf(bufS,"%d",&aPO->qlen);
    }
    DB_BIO DB_PrI("<< LoadBlastQueryLineData qlen=%d TRUE\n",aPO->qlen);
    return(TRUE);
}
/****************************************************************************
*   Check if line indicates end of record
*/
int BlastRecordEndLine(BLASTANS *aPO, char *bufS)
{
    if( (aPO->itype==BOT_WU) || (aPO->itype==BOT_SMWU) ) {
        if( strstr(bufS,"NONE") || EQSTRING(bufS,"Parameters",10) ) {
            return(TRUE);
        }
    }
    else {
        if( strstr(bufS,"No hits found") || EQSTRING(bufS,"Matrix:",7) ) {
            return(TRUE);
        }
    }
    return(FALSE);
}
/****************************************************************************
*   Load info for single query bloack (can be multiple alignment hits)
*   Return the number of actual aligment hits
*/
int AddBlastHitSetI(char *headS, FILE *inPF, BLASTANS *aPO, int hit, int warn)
{
    int thit, fpos, match, v1, v2;
    int coordsIA[4];
    char bufS[BBUFF_SIZE], tbufS[BBUFF_SIZE], *cPC;
    char hitS[BLBSIZE+1], scoreS[BLBSIZE+1], identS[BLBSIZE+1];
    char qseqS[BLBSIZE+1], sseqS[BLBSIZE+1];

    DB_BIO DB_PrI(">> AddBlastHitSetI hit=%d\n",hit);
    /***
    *   Start of hit block starts like:
    *   ">name"
    *   Get and save the "name"
    */
    if( headS[0]!='>') {
        return(FALSE);
    }
    INIT_S(hitS);
    sscanf(headS,">%s",hitS);
    /***
    *   asize = array size = max hit number to load info
    *   If too many, bail (maybe warn first?)
    */
    if(hit >= aPO->asize) {
        if(warn) {
            printf("# Warning: truncating hit set to %d\n",hit);
        }
        DB_BIO DB_PrI("<< AddBlastHitSetI FALSE asize=%d hit=%d\n",
            aPO->asize,hit);
        return(FALSE);
    }
    /***
    *   Remember file position to restore if we get to next query
    */
    fpos = ftell(inPF);
    /***
    *   Init count and record fields then eat lines until we run into nex query or out of file
    */
    thit = 0;
    INIT_S(scoreS); INIT_S(identS); INIT_S(qseqS); INIT_S(sseqS);
    coordsIA[0] = coordsIA[2] = TOO_BIG;
    coordsIA[1] = coordsIA[3] = 0;
    while(fgets(bufS,BLINEGRAB,inPF)) {
        DB_BIO {
            DB_PrI("+ "); fputs(bufS,stdout);
        }
        /***
        *   If we get to the next hit or record, bail
        */
        if( (bufS[0]=='>') || EQSTRING(bufS,"Query=",6) ) 
        {
            DB_BIO DB_PrI("+ Next record; Restoring file to %d (break)\n",fpos);
            fseek(inPF,fpos,0);
            break;
        }
        /***
        *   Score = starts single hit alignment block
        */
        if( (cPC = strstr(bufS,"Score =")) != NULL)
        {
            /***
            *    If we already have info, save that to a hit and reset everything
            */
            if( !NO_S(scoreS) ) {
                if( ! SaveOneHitRecI(aPO, hit + thit, hitS, scoreS, identS, qseqS, sseqS, coordsIA)) {
                    break;
                }
                thit++;
                INIT_S(scoreS); INIT_S(identS); INIT_S(qseqS); INIT_S(sseqS);
                coordsIA[0] = coordsIA[2] = TOO_BIG;
                coordsIA[1] = coordsIA[3] = 0;
            }
            ReplaceChars('\n',cPC,'\0',scoreS);
        }
        /***
        *   Identities line
        */
        if( (cPC = strstr(bufS,"Identities =")) != NULL)
        {
            NEXT_WORD(cPC); NEXT_WORD(cPC);
            sscanf(cPC,"%s",identS);
        }
        /***
        *   Query sequence line as:
        *   Query: <s_coord> <seq> <e_coord>
        *   7/5/14 RTK; Newer NCBI Blast+ output does not have the trailing colon;
xxx SHAM; Old version not working; No queries parse ... not sure if this is problem???
        */
        match = FALSE;
        if( aPO->itype == BOT_NCBOLD) {
            match = (strstr(bufS,"Query:")) ? TRUE : FALSE;
        }
        else {
            match = (strstr(bufS,"Query")) ? TRUE : FALSE;
        }
        if( match ) {
            INIT_S(tbufS);
            sscanf(bufS, "%*s %d %s %d", &v1, tbufS, &v2);
            coordsIA[0] = MIN_NUM(coordsIA[0],v1);
            coordsIA[1] = MAX_NUM(coordsIA[1],v2);
            if( (strlen(qseqS) + strlen(tbufS)) >= BLBSIZE) {
                if(warn) {
                        printf("# Warning: Query alignment length > %d ignored\n",BLBSIZE);
                }
            }
            else {
                strcat(qseqS,tbufS);
            }
        }
        /***
        *   Subject format same as Query
        */
        match = FALSE;
        if( aPO->itype == BOT_NCBOLD) {
            match = (strstr(bufS,"Sbjct:")) ? TRUE : FALSE;
        }
        else {
            match = (strstr(bufS,"Sbjct")) ? TRUE : FALSE;
        }
        if( match ) {
            INIT_S(tbufS);
            sscanf(bufS, "%*s %d %s %d", &v1, tbufS, &v2);
            coordsIA[2] = MIN_NUM(coordsIA[2],v1);
            coordsIA[3] = MAX_NUM(coordsIA[3],v2);
            if( (strlen(sseqS) + strlen(tbufS)) >= BLBSIZE) {
                if(warn) {
                        printf("# Warning: Subject alignment length > %d ignored\n",BLBSIZE);
                }
            }
            else {
                strcat(sseqS,tbufS);
            }
        }
        /***
        *   Remember position before eating next line
        */
        fpos = ftell(inPF);
    }
    /***
    *   Last or only record?
    */
    if( !NO_S(scoreS) ) {
        if( SaveOneHitRecI(aPO, hit + thit, hitS, scoreS, identS, qseqS, sseqS, coordsIA)) {
            thit++;
        }
    }

    DB_BIO DB_PrI("<< AddBlastHitSetI thit=%d\n",thit);
    return(thit);
}

/****************************************************************************
*   Load hits for record in blastout dci format
*/
int SaveOneHitRecI(BLASTANS *aPO, int hit, char *hitS, char *scoreS, char *identS,
    char *qseqS, char *sseqS, int *coordPI)
{
    if(hit >= aPO->asize) {
        return(FALSE);
    }
    strncpy(&aPO->hits[BLBSIZE*hit], hitS, BLBSIZE);
    strncpy(&aPO->scos[BLBSIZE*hit], scoreS, BLBSIZE);
    strncpy(&aPO->idens[BLBSIZE*hit], identS, BLBSIZE);
    strncpy(&aPO->qseqs[BLBSIZE*hit], qseqS, BLBSIZE);
    strncpy(&aPO->sseqs[BLBSIZE*hit], sseqS, BLBSIZE);
    aPO->qhsc[hit] = coordPI[0];
    aPO->qhec[hit] = coordPI[1];
    aPO->shsc[hit] = coordPI[2]; 
    aPO->shec[hit] = coordPI[3];;
    return(TRUE);
}

/****************************************************************************
*   Load hits for record in blastout dci format
*/
int LoadBdciRecordI(BLASTANS *aPO,FILE *inPF,int warn)
{
    int in,hit,n,line,expect;
    char bufS[BBUFF_SIZE];

    in = line = 0;
    expect = BAD_I;
    while(fgets(bufS,BLINEGRAB,inPF))
    {
        line++;
        if(line>HEAD_CHECK) { 
            break;
        }
        if( (bufS[0]=='#') && (strstr(bufS,"QueryNumber")) ) {
            in++;
        }
        if(!in) {
            continue;
        }
        INIT_S(bufS);
        if(!fgets(bufS,BLINEGRAB,inPF) ) {
            return(FALSE);
        }
        sscanf(bufS,"%*s %s %*s %d",aPO->query,&expect);
        if(!fgets(bufS,BLINEGRAB,inPF) ) {
            return(FALSE);
        }
        sscanf(bufS,"%*s %*s %*s %d",&aPO->qlen);
        break;
    }
/**
printf("query |%s| %d in=%d\n",aPO->query,aPO->qlen,in);
**/
    if(!in) {
        return(FALSE);
    }
    /***
    *   Load up to max records
    */
    hit = 0;
    while(hit<aPO->asize)
    {
        n = AddBlastDciRecI(inPF,aPO,hit,warn);
/*
printf(" hit=%d n=%d\n",hit,n);
*/
        if(n<1) {
            if(IS_BOG(n)) {
                return(FALSE);
            }
            break;
        }
        hit+=n;
    }
    /***
    *   Compare actual to alleged 
    */
    if(hit != expect) {
        PROBLINE;
        printf("Differing actual / alleged hit counts for query\n");
        printf("  Indicated: %d\n",expect);
        printf("  Found:     %d\n",hit);
        return(BOGUS);
    }
    aPO->nhits = hit;
    return(TRUE);
}
/***************************************************************************
*   TRUE for added hit
*   FALSE for end; nothing added
*   BOGUS for problem
*/
int AddBlastDciRecI(FILE *inPF,BLASTANS *aPO,int hit,int warn)
{
    int line,f,n,t;
    char bufS[BBUFF_SIZE];

    /***
    *   Get hit name from buffer then continue to eat lines till done
    */
    f = line = 0;
    while(fgets(bufS,BLINEGRAB,inPF))
    {
        line++;
        if(line>HEAD_CHECK) {
            break;
        }
        if(bufS[0]=='#') {
            if(strstr(bufS,"EndQuery")) {
                break;
            }
            if(strstr(bufS,"HitNumber")) {
                f++;
                continue;
            }
        }
        switch(f)
        {
            case 1:     /* Name & idents */
                if( !strstr(bufS,"HitStat") )
                {
                    PROBLINE;
                    printf("Derailed reading blastout hit stats\n");
                    fputs(bufS,stdout);
                    return(BOGUS);
                }
                n = t = BAD_I;
                sscanf(bufS,"%*s %*s %s (%d/%d)",&aPO->hits[BLBSIZE*hit],&n,&t);
                sprintf(&aPO->idens[BLBSIZE*hit],"%d/%d",n,t);
                f++;
                break;
            case 2:     /* Score line */
                ReplaceChars('\n',&bufS[2],'\0',&aPO->scos[BLBSIZE*hit]);
                f++;
                break;
            case 3:     /* Hit coords */
                if( !strstr(bufS,"HitCoords") )
                {
                    PROBLINE;
                    printf("Derailed reading blastout hit coords\n");
                    fputs(bufS,stdout);
                    return(BOGUS);
                }
                sscanf(bufS,"%*s %*s %d %d %d %d",&aPO->qhsc[hit], 
                    &aPO->qhec[hit], &aPO->shsc[hit],&aPO->shec[hit]);
                f++;
                break;
            case 4:     /* Query (sub)sequence */
/*
                if(bufS[0]!='Q')
                {
                    PROBLINE;
                    printf("Derailed reading blastout query seq\n");
                    fputs(bufS,stdout);
                    return(BOGUS);
                }
*/
                sscanf(bufS,"%*s %s",&aPO->qseqs[BLBSIZE*hit]);
                f++;
                break;
            case 5:     /* Subject (sub)sequence */
/*
                if(bufS[0]!='H')
                {
                    PROBLINE;
                    printf("Derailed reading blastout hit seq\n");
                    fputs(bufS,stdout);
                    return(BOGUS);
                }
                sscanf(&bufS[2],"%s %s",&aPO->hits[BLBSIZE*hit],
                    &aPO->sseqs[BLBSIZE*hit]);
*/
                sscanf(bufS,"%s %s",&aPO->hits[BLBSIZE*hit],
                    &aPO->sseqs[BLBSIZE*hit]);
                f++;
                break;
        }
        if(f>5)
            break;
    }
    /***
    *   there should have been six fields loaded after "HitNumber"
    */
    if(f==0)
    {
        return(FALSE);
    }
    else if(f!=6)
    {
        PROBLINE;
        printf("Derailed reading blastout hit record; got only %d fields\n",f);
        return(BOGUS);
    }
    return(TRUE);
}
/****************************************************************************
*   Report hit data held in bPO structure
*/
void DumpBlastout(BLASTOUT *bPO, BLASTANS *aPO,int qn, FILE *outPF)
{
    int i,max, mod_case;
    static int head = FALSE;
    char sS[BLBSIZE],bufS[DEF_BS];

    HAND_NFILE(outPF);
    /***
    *   Hit histogram, then bail
    */
    if(bPO->dhis > 0) {
        max = FillHitHistI(aPO,bPO->firstb,bPO->lastb,bPO->rre,bPO->do_con);
        if(bPO->chis) {
            IntegrateHist(aPO);
        }
        max = MAX_NUM(max,(bPO->phis));
        /***
        *   Header row?
        */
        if(!head) {
            fprintf(outPF,"%-15s","QueryName");
            for(i=bPO->dhis; i<=max; i++) {
                fprintf(outPF," %3d",i);
            }
            fprintf(outPF,"\n");
            head++;
        }
        /***
        *   Lable and each count
        */
        fprintf(outPF,"%-15s",aPO->query);
        for(i=bPO->dhis; i<=max; i++) {
            fprintf(outPF," %3d",aPO->ihist[i]);
        }
        fprintf(outPF,"\n");
        return;
    }
    /***
    *   Per query reports, then bail
    */
    if(bPO->dsum) {
        fprintf(outPF,"%12s %4d %4d\n",aPO->query,aPO->nhits,aPO->maxhit);
        return;
    }
    if(bPO->do_dfml > 0) {
        DumpBlastFracMatchList(bPO,aPO,outPF);
        return;
    }
    if( bPO->dmbc ) {
        DumpBlastMatchBaseCounts(bPO,aPO,outPF);
        return;
    }
    /***
    *   Dumping hit seqs/coords
    */
    if( bPO->dseq || bPO->dhc || bPO->dmbc ) {
        fprintf(outPF,"#\n");
        fprintf(outPF,"# QueryNumber %d\n",qn);
        fprintf(outPF,"# %s HitTotal %d\n",aPO->query,aPO->nhits);
        fprintf(outPF,"# %s QueryLen %d\n",aPO->query,aPO->qlen);
    }
    else {
        fprintf(outPF,"%12s %4d\n",aPO->query,aPO->nhits);
    }
    /***
    *   Dump hit story
    */
    for(i=0;i<aPO->nhits;i++) {
        FillBlastRecHitName(aPO, i, sS);
        if(bPO->dhit) {
            fprintf(outPF,"# HitNumber %d\n",i+1);
            if(bPO->dseq)
                fprintf(outPF,"# ");
            FillBlastRecHitStatLineI(aPO, i, bufS);
            fprintf(outPF,"HitStat %s\t%s\n",sS,bufS);
        }
        if(bPO->dsco) {
            if(bPO->dseq)
                fprintf(outPF,"# ");
            fprintf(outPF,"%s\n",&aPO->scos[i*BLBSIZE]);
        }
        if(bPO->dhc) {
            if(bPO->dseq)
                fprintf(outPF,"# ");
            FillBlastRecHitCoordsLineI(aPO, i, bufS);
            fprintf(outPF,"HitCoords\t%d\t%s\t", i+1, bufS);
            FillBlastRecHitStatLineI(aPO, i, bufS);
            fprintf(outPF,"%s\n",bufS);
        }
        if(bPO->dseq) {
            if( (bPO->do_smc) || (bPO->do_sml) ) {
                mod_case = (bPO->do_smc) ? 1 : -1;
            }
            else {
                mod_case = 0;
            }
            DumpBlastoutHitSeqs(aPO, i, mod_case, bPO->nsqh, outPF);
        }
    }
    if(bPO->do_dci) {
        fprintf(outPF,"# EndQuery %d\n",qn);
    }
    fflush(outPF);
}
/**************************************************************************
*   Dump sequences 
*/
void DumpBlastoutHitSeqs(BLASTANS *aPO, int hit, int mod_case, char *nsqhS, FILE *outPF)
{
    char qnameS[NSIZE],snameS[NSIZE],qseqS[BLBSIZE+1],sseqS[BLBSIZE+1];
    int i;

    HAND_NFILE(outPF);
    sprintf(qnameS,"%s", aPO->query);
    FillBlastRecHitName(aPO, hit, snameS);
    sprintf(qseqS,"%s",&aPO->qseqs[hit*BLBSIZE]);
    sprintf(sseqS,"%s",&aPO->sseqs[hit*BLBSIZE]);
    /***
    *   Sequence marking via case?
    */
    if (mod_case) {
        if ( mod_case < 0 ) {
            Lowerize(qseqS);
            Lowerize(sseqS);
        }
        else {
            Upperize(qseqS);
            Upperize(sseqS);
        }
        i = 0;
        while(isgraph(qseqS[i])) {
            if (sseqS[i] != qseqS[i]) {
                if ( mod_case < 0 ) {
                    sseqS[i] = toupper(sseqS[i]);
                }
                else {
                    sseqS[i] = tolower(sseqS[i]);
                }
            }
            i++;
        }
    }
    /***
    *   
    */
    fprintf(outPF,"%s\t%s\n",qnameS,qseqS);
    if ( (nsqhS) && (! NO_S(nsqhS)) ) {
        fprintf(outPF,"%s%s%s\t%s\n",qnameS,nsqhS,snameS,sseqS);
    }
    else {
        fprintf(outPF,"%s\t%s\n",snameS,sseqS);
    }
}
/**************************************************************************
*   Fill name for blasthit record
*/
int FillBlastRecHitName(BLASTANS *aPO, int r, char *bufS)
{
    if( (r < 0) || ( r >= aPO->nhits) )
    {
        return(FALSE);
    }
    strncpy(bufS,&aPO->hits[r*BLBSIZE],NAME_MAX);
    bufS[NAME_MAX] = '\0';
    return(TRUE);
}
/**************************************************************************
*   Fill coords for match
*/
int FillBlastRecHitCoordsLineI(BLASTANS *aPO, int r, char *bufS)
{
    int qs,qe,ss,se,t;
    char qS[NAME_MAX+1], sS[NAME_MAX+1], strandS[DEF_BS];

    if( (r < 0) || ( r >= aPO->nhits) )
    {
        return(FALSE);
    }
    /*
    *   Query part
    */
    qs = aPO->qhsc[r]; 
    qe = aPO->qhec[r];
    strncpy(qS,aPO->query,NAME_MAX);
    qS[NAME_MAX] = '\0';
    /*
    *   Subject (hit) part
    */
    ss = aPO->shsc[r]; 
    se = aPO->shec[r];
    FillBlastRecHitName(aPO, r, sS);
    if ( ss < se ) 
    {
        sprintf(strandS,"Pos");
    }
    else
    {
        sprintf(strandS,"Neg");
        t = ss;
        ss = se;
        se = t;
    }
    sprintf(bufS,"%s\t%d\t%d\t%s\t%d\t%d\t%s",
        qS,qs,qe, sS,ss,se,strandS);
    return(TRUE);
}
/**************************************************************************
*   Fill stat report for match
*/
int FillBlastRecHitStatLineI(BLASTANS *aPO, int r, char *bufS)
{
    if( (r < 0) || ( r >= aPO->nhits) ) {
        return(FALSE);
    }
    sprintf(bufS,"(%s)\t%2d/%2d  %4.3f  %4.3f",
        &aPO->idens[r*BLBSIZE], 
        aPO->hnums[r], aPO->hlens[r], aPO->hfmb[r],
        RNUM(aPO->hnums[r]) / RNUM(aPO->qlen) );
    return(TRUE);
}
/**************************************************************************
*   Dump (sorted) list of hit match fractions
*/
void DumpBlastFracMatchList(BLASTOUT *bPO,BLASTANS *aPO,FILE *outPF)
{
    int i, num, den;
    DOUB fR, sortR[DEF_MAXHITS], normD;

    /***
    *   Normalization factor; If full length or no len restriction...
    */
    if( (bPO->do_dffl) || (bPO->firstb<1) ) {
        normD = DNUM(aPO->qlen);
    }
    else {
        normD = DNUM(bPO->lastb - bPO->firstb + 1);
    } 
    /***
    *   Collect fractions and sort
    */
    InitArrayI(sortR, IS_DOUB, 0, DEF_MAXHITS, 0.0);
    for(i=0; i<bPO->do_dfml; i++) {
        if(i<aPO->nhits){
            num = HitMatchCountI(aPO, i, bPO->firstb, bPO->lastb, bPO->rre, 
                bPO->do_con, &den);
            fR = RNUM(num)/normD;
        }
        else {
            fR = 0.0;
        }
        if(i<DEF_MAXHITS) {
            sortR[i] = fR;
        }
    }
    /***
    *   Now dump the list
    */
    i = MIN_NUM(i,DEF_MAXHITS);
    SortArray(sortR, IS_DOUB, i, SORT_DECEND);
    fprintf(outPF,"%-15s\t",aPO->query);
    for(i=0; i<bPO->do_dfml; i++) {
        fprintf(outPF,"%2.2f ",sortR[i]);
    }
    fprintf(outPF,"\n");
    return;
}
/**************************************************************************
*   Fill and dump base-matching counts along query sequence
*/
void DumpBlastMatchBaseCounts(BLASTOUT *bPO,BLASTANS *aPO,FILE *outPF)
{
    int i,j,b,hist[BLBSIZE],gap,qs,qe;
    char *qPC,*sPC;

    if(aPO->qlen < 1) {
        printf("SHAM: %s len=%d\n",aPO->query,aPO->qlen);
    }
    for(i=0;i<aPO->qlen;i++) {
        hist[i] = 0;
    }
    for(i=0;i<aPO->nhits;i++) {
        qs = aPO->qhsc[i]; 
        qe = aPO->qhec[i];
        j = gap = 0;
        qPC = &aPO->qseqs[BLBSIZE*i];
        sPC = &aPO->sseqs[BLBSIZE*i];
        while(isgraph(INT(qPC[j]))) {
            if(qPC[j] != sPC[j]) {
                if(qPC[j] == '-') {
                    gap++;
                }
                j++;
                continue;
            }
            if(qs<qe)
                b = qs+j-1-gap;
            else
                b = qs-j-1+gap;
            if( (b<0) || (b>=aPO->qlen) ) {
                printf("Bogus base bin b=%d i=%d j=%d gap=%d qlen=%d\n",
                    b,i,j,gap,aPO->qlen);
                ERR("DumpBlastMatchBaseCounts","Sham with matching bin");
                return;
            }
            hist[b] += 1;
            j++;
            if(!isgraph(INT(sPC[j])))
                break;
        }
    }
    /***
    *   Dump the counts
    */
    fprintf(outPF,"%-15s\t",aPO->query);
    i = 0;
    while(i<aPO->qlen) {
        b = TRUE;
        /***
        *   Restricted bases?
        */
        if(bPO->firstb>0) {
            if(bPO->rre) {
                if( ((aPO->qlen-i)<bPO->firstb) || 
                    ((aPO->qlen-i)>=bPO->lastb) )
                    b = FALSE;
            }
            else {
                if( ((i+1)<bPO->firstb) || (i>=bPO->lastb) )
                    b = FALSE;
            }
        }
        if(b) {
            fprintf(outPF," %d",hist[i]);
        }
        i++;
    }
    fprintf(outPF,"\n");
}
/****************************************************************************
*
*/
int GuessBlastInputTypeI(FILE *inPF)
{
    int line,fpos,itype;
    char bufS[BBUFF_SIZE];

    DB_BIO DB_PrI(">> GuessBlastInputTypeI\n");
    fpos = ftell(inPF);
    itype = BOGUS;
    line = 0;
    while(fgets(bufS,BLINEGRAB,inPF))
    {
        line++;
        if(bufS[0]=='#') {
            Upperize(bufS);
            if(strstr(bufS,"BLASTOUT")) {
                itype = BOT_BO;
            }
        }
        else if(EQSTRING(bufS,"BLAST",5)) {
            if(strstr(bufS,"WashU")) {
                if(strstr(bufS,"linux")) {
                    itype = BOT_SMWU;
                }
                else {
                    itype = BOT_WU;
                }
            }
            else {
                if(!strstr(bufS,"+")) {
                    itype = BOT_NCBOLD;
                }
                else {
                    itype = BOT_NCBI;
                }
            }
            break;
        }
        if(line>HEAD_CHECK)
        {
            break;
        }
    }
    fseek(inPF,fpos,0);
    DB_BIO DB_PrI("<< GuessBlastInputTypeI %d\n",itype);
    return(itype);
}
/**************************************************************************/
void FillBlastFormatString(int type,char *nameS)
{
    switch(type)
    {
        case BOT_NCBI:      sprintf(nameS,"NCBI (Blast+)"); break;
        case BOT_NCBOLD:    sprintf(nameS,"NCBI (old Blast)");  break;
        case BOT_WU:        sprintf(nameS,"Wu (normal)");   break;
        case BOT_SMWU:      sprintf(nameS,"Wu (\"smart\")");    break;
        case BOT_BO:        sprintf(nameS,"BlastOut (dci)");    break;
        default:            sprintf(nameS,"???");
    }
}
/**************************************************************************/
int HandleBlastoutOutfileI(BLASTOUT *bPO,BLASTANS *aPO,int q,int head)
{
    /***
    *   Output per query
    */
    if(!NO_S(bPO->opq))
    {
        CHECK_NFILE(bPO->out,bPO->outname);
        if(q>0)
        {
            sprintf(bPO->outname,"%s.%s",aPO->query,bPO->opq);
            if(!(bPO->out = OpenUFilePF(bPO->outname,"w",NULL)))
            {
                HAND_NFILE(bPO->out);
                return(FALSE);
            }
        }
    }
    /***
    *   Normal output file, once
    */
    else
    {
        if( (!NO_S(bPO->outname)) && (bPO->out==NULL) )
        {
            if(!(bPO->out = OpenUFilePF(bPO->outname,"w",NULL)))
            {
                HAND_NFILE(bPO->out);
                return(FALSE);
            }
        }
    }
    /***
    *   Set to default
    */
    HAND_NFILE(bPO->out);
    return(TRUE);
}
/****************************************************************************
*
*/
int SetBlastansFileI(char *inS,BLASTANS *aPO)
{
    VALIDATE(aPO,BLASTANS_ID);
    if(!(aPO->in = OpenUFilePF(inS,"r",NULL)))
    {
        return(FALSE);
    }
    aPO->itype = GuessBlastInputTypeI(aPO->in);
    if(IS_BOG(aPO->itype))
    {
        PROBLINE;
        printf("Unrecognized input format\n");
        return(FALSE);
    }
    strcpy(aPO->input,inS);
    return(TRUE);
}
