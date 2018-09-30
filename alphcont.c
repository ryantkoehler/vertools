/*
* alphcont.c
*
* Copyright 2018 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
* See https://www.verdascend.com/ for more
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "alphcont.h"

#define DB_ALFC if(DB[99])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(AlphContI(argc,argv),NULL) ); }
/**************************************************************************/
void AlphContUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>   Sequence file\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -out XXX   Set output to XXX\n");
    printf("   -bran # #  Consider base range # to # (i.e. sub-seqs)\n");
    printf("   -rre       Base range relative to end; i.e. Backwards\n");
    printf("   -row XXX   Report (max) row of string XXX\n");
    printf("   -con XXX   Report count of string XXX\n");
    printf("   -skew XXX  Report skew of string XXX\n");
    printf("   -cent XXX  Report centroid of string XXX\n");
    printf("   -btab      Base table output: ACGTSWRYMK content\n");
    printf("   -rtab      Rows table output: ACGTSWRYMK max-length rows\n");
    printf("   -dtab      Dinucleotide table output: AA AC AG ... TC TG TT\n");
    printf("   -ecc       Explicit char counts (e.g. \"S\"==S, not G or C)\n");
    printf("   -dnum      Dump count of words (number); Default is fraction\n");
    printf("   -ign       Ignore N's when computing fraction matching\n");
    printf("   -cwin #    Evaluate counts in windows of # bases\n");
    printf("   -dcw       Dump (cwin) values along length\n");
    printf("   -dcc       Dump cumulative counts along sequence\n");
    printf("   -dfg       Dump flagging with matches along sequence length\n");
    printf("   -dfs       Dump flagging and sequence\n");
    printf("   -flg # #   Flag sequences with values # to #\n");
    printf("   -ds        Dump (report) sequences appended as last column\n");
    printf("   -eraw      Extract flagged sequences in raw format\n");
    printf("   -not       Invert flagging criteria\n");
    printf("   -oval      Output only values (no seqs, e.g. with -eraw)\n");
    printf("   -stat      Output only stats with -flg option\n");
    printf("   -win # #   Sample sequence windows of size # stepping by #\n");
    printf("   -wst #     Sample windows starting base\n");
    printf("   -wen #     Sample windows ending base\n");
    printf("\n");
}
/**************************************************************************
*
*/
int AlphContI(int argc,char **argv)
{
    int ifas, iraw, iseq, dcc, dfs, dfg, dcw, ok, oval, eraw, start, slen, n;
    char snameS[NSIZE],winS[NSIZE];
    ALPHCONT *alPO;

    alPO = CreateAlphcontPO();
    iraw = ifas = iseq = dcc = dfs = dfg = dcw = oval = eraw = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -igp B -bran I2 -flg R2 -row S -con S -not B -dnum B\
        -dfs B -dcc B -iraw B -ifas B -cwin I -dfg B -dcw B -rre B\
        -eraw B -oval B -stat B -ecc B -btab B -skew S -cent S\
        -win I2 -wst I -wen I -ign B -rtab B -dtab B -ds B -iseq B",
        alPO->inname, alPO->outname, &alPO->igprob, 
        &alPO->firstb,&alPO->lastb,
        &alPO->min,&alPO->max, alPO->row, alPO->con,
        &alPO->do_not, &alPO->do_dfr, &dfs, &dcc, &iraw, &ifas, &alPO->cwin,
        &dfg, &dcw, &alPO->do_rre, &eraw, &oval,
        &alPO->do_stat, &alPO->do_ecc, &alPO->do_btab,
        &alPO->skew, &alPO->cent,
        &alPO->winsize,&alPO->winstep, &alPO->winst, &alPO->winen, 
        &alPO->do_ign, &alPO->do_rtab, &alPO->do_dtab, &alPO->do_ds,
        &iseq,
        (int *)NULL))
    {
        AlphContUse();
        CHECK_ALPHCONT(alPO);
        return(TRUE);
    }
    /***
    *   Set input format to explitly chosen, else guess (if not stdin)
    */
    alPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas, alPO->inname, TRUE);
    if(!alPO->iform) {
        printf("Problem with input seq(s)\n");
        CHECK_ALPHCONT(alPO);
        return(FALSE);
    }
    /***
    *   output options?
    */
    if(!BAD_REAL(alPO->min)) {
        alPO->do_qual = TRUE;
    }
    else {
        alPO->do_qual = FALSE;
        alPO->min = -TOO_BIG_R;
        alPO->max = TOO_BIG_R;
    }
    if(dcc) {
        alPO->owhat = ALCO_DCC;
    } 
    else if(dcw) {
        alPO->owhat = ALCO_DCW;
    }
    else if(dfg) {
        alPO->owhat = ALCO_DFG;
    }
    else if(dfs) {
        alPO->owhat = ALCO_DFS;
    }
    else if(eraw) {
        alPO->owhat = ALCO_SEQ;
    }
    else {
        alPO->owhat = ALCO_VAL;
    }
    if(oval) {
        alPO->owhat = ALCO_VAL;
    }
    if(!CheckAlcoOptionsI(alPO)) {
        printf("Problem with command options\n");
        ABORTLINE;
        CHECK_ALPHCONT(alPO);
        return(FALSE);
    }
    /***
    *   Open in/out files
    */
    if(!OpenAlcoFilesI(alPO)) {
        printf("Problem with alco files\n");
        ABORTLINE;
        CHECK_ALPHCONT(alPO);
        return(FALSE);
    }
    HandleAlcoHeader(alPO,alPO->out);
    /***
    *   Loop through the collection
    */
    n = 0;
    while(TRUE) {
        /***
        *   Parse (full) sequence; FALSE = done
        */
        n++;
        ok = ParseSeqI(alPO->in, alPO->iform, n, alPO->iclean, TRUE, alPO->fseq);
        if(ok==FALSE) {
            break;
        }
        if(ok!=TRUE) {
            continue;
        }
        /***
        *   Single input or window sampling?
        */
        FillSeqNameStringI(alPO->fseq, snameS, NSIZE-1);
        if(alPO->winsize < 1) {
            CopySeqI(alPO->fseq, alPO->seq, -1, -1);
            if(!AlphconProcSeqI(alPO, alPO->seq, alPO->fseq, snameS)) {
                printf("Problem processing seq |%s|\n",snameS);
                ABORTLINE;
                return(FALSE);
            }
        }
        else {
            slen = GetSeqLenI(alPO->fseq);
            if(!IS_BOG(alPO->winen)) {
                slen = MIN_NUM(slen,alPO->winen);
            }
            start = alPO->winst - 1;
            while( (start + alPO->winsize) < slen) {
                CopySeqI(alPO->fseq, alPO->seq, start, alPO->winsize);
                sprintf(winS,"%s_Win-%03d_%03d",snameS, 
                    start+1, start+1+alPO->winsize);
                if(!AlphconProcSeqI(alPO, alPO->seq, alPO->fseq, winS)) {
                    printf("Problem processing seq, window size %d |%s|\n", alPO->winsize,snameS);
                    ABORTLINE;
                    return(FALSE);
                }
                start += alPO->winstep;
            }
        }
    }
    /***
    *   All done
    */
    CHECK_ALPHCONT(alPO);
    return(TRUE);
}
/*****************************************************************************
*   Process passed sequence
*/
int AlphconProcSeqI(ALPHCONT *alPO, SEQ *seqPO, SEQ *fseqPO, char *nameS)
{
    int ok;

    /***
    *   Clean & init for seq counting
    */
    HandleAlcoSeqCleaningI(seqPO, alPO->do_ecc);
    /***
    *   Trim sequence?
    */
    if(!HandleAlcoSubseqI(alPO, seqPO, fseqPO)) {
        return(FALSE);
    }
    /***
    *   Get counts 
    */
    if(!HandleSeqCountsI(alPO, seqPO)) {
        return(FALSE);
    }
    /***
    *   Output
    */
    ok = IsAlphcontSeqOkI(alPO);
    HandleAlcoOutputI(alPO, seqPO, nameS, ok, alPO->out);
    return(TRUE);
}
/*****************************************************************************
*   Create data struct
*/
ALPHCONT *CreateAlphcontPO()
{
    ALPHCONT *alPO;

    if(! (alPO = (ALPHCONT *)ALLOC(1,sizeof(ALPHCONT)) ) ) {
        printf("# Failed to allocate working object\n");
        return(NULL);
    }
    alPO->ID = ALPHCONT_ID;
    alPO->seq = CreateSeqPO(0,NULL,NULL);
    alPO->fseq = CreateSeqPO(0,NULL,NULL);
    alPO->comp = CreateSeqcompPO();
    if( (!alPO->seq) || (!alPO->fseq) || (!alPO->comp ) ) {
        printf("# Failed to allocate seq space\n");
        CHECK_ALPHCONT(alPO);
        return(NULL);
    }
    InitAlphcont(alPO);
    return(alPO);
}
/*****************************************************************************
*   Free datastructure and substructs
*/
int DestroyAlphcontI(ALPHCONT *alPO)
{
    VALIDATE(alPO,ALPHCONT_ID);
    CHECK_FILE(alPO->in);
    CHECK_NFILE(alPO->out,alPO->outname);
    CHECK_SEQ(alPO->seq);
    CHECK_SEQ(alPO->fseq);
    CHECK_SEQCOMP(alPO->comp);
    CHECK_FREE(alPO->smask);
    FREE(alPO);
    return(TRUE);
}
/*****************************************************************************
*   Set null / default values
*/
void InitAlphcont(ALPHCONT *alPO)
{
    VALIDATE(alPO,ALPHCONT_ID);
    INIT_S(alPO->inname);
    INIT_S(alPO->outname);
    alPO->in = NULL;
    alPO->iform = BOGUS;
    alPO->iclean = SCLEAN_MID;
    alPO->out = NULL;
    alPO->igprob = FALSE;
    alPO->do_stat = FALSE;
    alPO->quiet = FALSE;
    alPO->do_not = FALSE;
    alPO->firstb = -1;    
    alPO->lastb = -1;
    alPO->do_rre = FALSE;
    alPO->cwin = 1;
    alPO->do_dfr = TRUE;
    alPO->min = alPO->max = BAD_R;
    alPO->do_qual = FALSE;
    alPO->do_btab = alPO->do_rtab = FALSE;
    alPO->do_dtab = FALSE;
    sprintf(alPO->con,DEF_CON_S);
    INIT_S(alPO->row);
    INIT_S(alPO->skew);
    INIT_S(alPO->cent);
    alPO->do_uc = FALSE;
    alPO->do_lc = FALSE;
    alPO->winsize = -1;
    alPO->winstep = -1;
    alPO->winst = 1;
    alPO->winen = BOGUS;
    alPO->smask = NULL;
    alPO->n_smask = 0;
}
/*************************************************************************
*   Check for option consistency
*   Set flags for options
*/
int CheckAlcoOptionsI(ALPHCONT *alPO)
{
    alPO->rowlen = strlen(alPO->row);
    alPO->conlen = strlen(alPO->con);
    alPO->skewlen = strlen(alPO->skew);
    alPO->centlen = strlen(alPO->cent);
    Upperize(alPO->row);
    Upperize(alPO->con);
    Upperize(alPO->skew);
    Upperize(alPO->cent);
    /***    
    *   If row table output, set for this
    */
    if ( alPO->do_rtab ) {
        if (alPO->do_dtab) {
            PROBLINE;
            printf("Rows don't work with dinucleotides (-dtab)\n");
            return(FALSE);
        }
        alPO->do_dfr = FALSE;
    }
    /***
    *   Skew is reported as a fraction
    */
    if (alPO->skewlen > 0) {
        alPO->do_dfr = TRUE;
    }
    /***
    *   Row / skew restrictions
    */
    if (alPO->rowlen > 0) {
        if( alPO->owhat==ALCO_DCC) {
            PROBLINE;
            printf("Row doesn't work with -dcc\n");
            return(FALSE);
        }
    }
    if (alPO->skewlen > 0) {
        if( (alPO->owhat==ALCO_DFG) || (alPO->owhat==ALCO_DCC) ||
            (alPO->owhat==ALCO_DFS) )
        {
            PROBLINE;
            printf("Skew doesn't work with -dcc, -dfg, or -dfs\n");
            return(FALSE);
        }
    }
    /***
    *   Window on counts?
    */
    if( (alPO->cwin>1) || (alPO->owhat == ALCO_DCW) ) {
        if(alPO->rowlen > 0) {
            PROBLINE;
            printf("Window doesn't work with row counts\n");
            return(FALSE);
        }
        else if(alPO->cwin < 1) {
            PROBLINE;
            printf("Window can't be less than one; %d won't work\n",alPO->cwin);
            return(FALSE);
        }
    }
    if( (alPO->firstb > 0) && (alPO->owhat == ALCO_DCW) ) {
        PROBLINE;
        printf("Can't dump window counts with restricted base range\n");
        return(FALSE);
    }
    /***
    *   Output
    */ 
    if(alPO->rowlen > 0) {
        alPO->do_dfr = FALSE;
    }
    /***
    *   Qualifying is only cool with seq / number output
    */
    if( alPO->do_qual ) {
        switch(alPO->owhat) {
            case ALCO_VAL:
            case ALCO_SEQ:
            case ALCO_DCW:
                break;
            default:
                PROBLINE;
                printf("Flagging only works with simple output\n");
                printf("  (any sequence flagging won't work here)\n");
                return(FALSE);
        }
    }
    if(alPO->winsize > 0) {
        if( (alPO->winstep < 1) || (alPO->winst < 1) || 
            (!IS_BOG(alPO->winen)&&(alPO->winen < alPO->winst)) )
        {
            PROBLINE;
            printf("Window sampling values won't work!\n");
            printf("size=%d\nstep=%d\nstart=%d\nend=%d\n",
                alPO->winsize,alPO->winstep,alPO->winst,alPO->winen);
            return(FALSE);
        }
    }
    return(TRUE);
}
/***************************************************************************
*   Open any needed files or die
*/
int OpenAlcoFilesI(ALPHCONT *alPO)
{
    if(!(alPO->in=OpenUFilePF(alPO->inname,"r",NULL))) {
        return(FALSE);
    }
    if(!NO_S(alPO->outname)) {
        if(!(alPO->out=OpenUFilePF(alPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    HAND_NFILE(alPO->out);
    return(TRUE);
}
/****************************************************************************
*   Report what's getting counted / dumped
*/
void HandleAlcoHeader(ALPHCONT *alPO, FILE *outPF)
{
    char bufS[DEF_BS];

    HAND_NFILE(outPF);
    if(alPO->quiet)
    {
        return;
    }
    VersionSplash(outPF,VERSION_S,"# ",TRUE);
    fprintf(outPF,"#\n");
    fprintf(outPF,"# Input: %s\n",alPO->inname);
    fprintf(outPF,"# Base range:");
    if(alPO->firstb > 0) {
        fprintf(outPF," %d to %d",alPO->firstb, alPO->lastb);
        if(alPO->do_rre) {
            fprintf(outPF," from end\n");
        }
        else {
            fprintf(outPF," from start\n");
        }
    }
    else {
        fprintf(outPF," all\n");
    }
    /***
    *   What's reported?
    */
    if( alPO->do_btab || alPO->do_dtab || alPO->do_rtab ) {
        if(alPO->do_rtab) {
            sprintf(bufS,"A C G T Max1 S W R Y M K Max2");
            ReplaceChars(' ', bufS, '\t', bufS);
            fprintf(outPF,"# Base row table\n");
            fprintf(outPF,"#\tSeq\t%s\n",bufS);
        }
        else {
            if(alPO->do_dtab) {
                sprintf(bufS,"AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT");
                ReplaceChars(' ', bufS, '\t', bufS);
                fprintf(outPF,"# Dinucleotide table\n");
                fprintf(outPF,"#\tSeq\t%s\n",bufS);
            }
            else {
                sprintf(bufS,"A C G T S W R Y M K");
                ReplaceChars(' ', bufS, '\t', bufS);
                fprintf(outPF,"# Base count table\n");
                fprintf(outPF,"#\tSeq\t%s\n",bufS);
            }
        }
    }
    else if(!NO_S(alPO->skew)) {
        fprintf(outPF,"# Evaluating skew of: %s\n",alPO->skew);
        fprintf(outPF,"#    Col-2 = (Up - Down) / (Side length)\n");
        fprintf(outPF,"#    Col-3 = (Up - Down) / (Total matches)\n");
        fprintf(outPF,"#    Col-4 = Side length\n");
        fprintf(outPF,"#    Col-5 = ABS| (Up - Down) / (Side length)   |\n");
        fprintf(outPF,"#    Col-6 = ABS| (Up - Down) / (Total matches) |\n");
    }
    else if(!NO_S(alPO->cent)) {
        fprintf(outPF,"# Evaluating centroid of: %s\n",alPO->cent);
    }
    else if(!NO_S(alPO->row)) {
        fprintf(outPF,"# Counting max row of: %s\n",alPO->row);
    }
    else {
        fprintf(outPF,"# Counting content of: %s\n",alPO->con);
    }
    /***
    *   Explicit or degenerate
    */
    if( (!alPO->do_btab) && (!alPO->do_dtab) && (!alPO->do_rtab) ) {
        if(alPO->do_ecc) {
            fprintf(outPF,"# Explicit char matching (e.g. 'S' = 'S')\n");
        }
        else {
            fprintf(outPF,"# Degen char matching (e.g. 'S' = 'C' or 'G')\n");
        }
    }
    if(alPO->winsize > 0) {
        fprintf(outPF,"# Window sampling size %d, step %d, starting at %d to",
            alPO->winsize, alPO->winstep, alPO->winst);
        if(!IS_BOG(alPO->winen)) {
            fprintf(outPF," %d\n",alPO->winen);
        }
        else {
            fprintf(outPF," end\n");
        }
    }
    /***
    *   How's it reported?
    */
    switch(alPO->owhat)
    {
        case ALCO_DCC:
            fprintf(outPF,"# Cumulative counts\n");
            break;
        case ALCO_DFS:
            fprintf(outPF,"# Flagged sequences\n");
            break;
        case ALCO_DFG:
            fprintf(outPF,"# Flagging for sequences\n");
            break;
        case ALCO_DCW:
            fprintf(outPF,"# Window (%d) count values\n",alPO->cwin);
            break;
        case ALCO_SEQ:
            fprintf(outPF,"# Sequence subset\n");
            if(alPO->do_not)
            {
                fprintf(outPF,"#   Limtited to values NOT %4.3f to %4.3f\n",
                    alPO->min,alPO->max);
            }
            else
            {
                fprintf(outPF,"#   Limited to values %4.3f to %4.3f\n",
                    alPO->min,alPO->max);
            }
            break;
    }
}
/****************************************************************************
*   Clean up non-standard chars and init mask
*/
int HandleAlcoSeqCleaningI(SEQ *seqPO, int eec)
{
    VALIDATE(seqPO,SEQ_ID);
    /***
    *   If not looking at explict characters, clean string to DNA alphabet
    *   Last two FALSE args = don't switch SNP to IUPAC, don't mask lowercase
    */
    if(eec) {
        seqPO->len = CleanUpSeqI(seqPO->seq,seqPO->len,seqPO->seq,FALSE,FALSE);
    }
    return(TRUE);
}
/***************************************************************************
*   Shrink to range if needed
*/
int HandleAlcoSubseqI(ALPHCONT *alPO, SEQ *seqPO, SEQ *fseqPO)
{
    int len, flen;

    if(alPO->firstb >= 0) {
        len = alPO->lastb - alPO->firstb + 1;
        SetCaseSeqSubseqI(fseqPO, FALSE, -1, -1);
        if(alPO->do_rre) {
            NarrowSeqI(seqPO,alPO->firstb-1,len,REVERSE,FALSE);
            flen = GetSeqLenI(fseqPO);
            SetCaseSeqSubseqI(fseqPO, TRUE, flen - alPO->lastb, flen - alPO->firstb + 1);
        }
        else {
            NarrowSeqI(seqPO,alPO->firstb-1,len,FORWARD,FALSE);
            SetCaseSeqSubseqI(fseqPO, TRUE, alPO->firstb-1, alPO->lastb);
        }
    }
    return(TRUE);
}
/****************************************************************************
*   What to calucluate
*/
int HandleSeqCountsI(ALPHCONT *alPO, SEQ *seqPO)
{
    int ok;

    VALIDATE(seqPO,SEQ_ID);
    alPO->val = alPO->val2 = 0.0;
    /***
    *   Base table, else skew, else centroid, else row or content
    *   Args: Seq, len, only ACGT, structure for results
    */
    if( alPO->do_btab || alPO->do_dtab || alPO->do_rtab ) {
        ok = GetSeqCompositionI(seqPO->seq, seqPO->len, FALSE, alPO->comp);
    }
    else if(alPO->skewlen>0) {
        ok = GetSeqSkewI(alPO,seqPO);
    }
    else if(alPO->centlen>0) {
        ok = GetSeqCentI(alPO,seqPO);
    }
    else {
        ok = ScanSeqCountsI(alPO, seqPO);
        if (ok) {
            HandleFractionValueI(alPO, seqPO);
        }
    }
    return(ok);
}
/****************************************************************************
*   Skew of query word in sequence up/down parts
*/
int GetSeqSkewI(ALPHCONT *alPO, SEQ *seqPO)
{
    int i,n,m,t,umat,dmat;

    /***
    *   Number of positions to screen and position of mid-point
    */
    n = INT(seqPO->len / 2) - alPO->skewlen;
    if(n<1) {
        return(TRUE);
    }
    m = seqPO->len / 2;
    if(ODD_NUM(seqPO->len)) {
        m++;
    }
    /***
    *   Evaluate number of instances in upstream and downstream parts
    */
    umat = dmat = t = 0;
    for(i=0;i<=n;i++)
    {
        t++;
        if(QSMatchI(alPO->skew,&seqPO->seq[i],alPO->skewlen,alPO->do_ecc)) {
            umat++;
        }
        if(QSMatchI(alPO->skew,&seqPO->seq[m+i],alPO->skewlen,alPO->do_ecc)) {
            dmat++;
        }
    }
    /***
    *   Two values
    *   Skew = up-down / total samples
    *   Skew2 = up-down / total matches
    */
    n = umat + dmat;
    alPO->val = DNUM(umat - dmat)/DNUM(t);
    alPO->val2 = DNUM(umat - dmat)/DNUM(n);
    return(TRUE);
}
/****************************************************************************
*   Centroid of query word in sequence 
*/
int GetSeqCentI(ALPHCONT *alPO, SEQ *seqPO)
{
    int i,n,pos,mat;

    /***
    *   Number of positions to screen 
    */
    n = seqPO->len - alPO->conlen;
    if(n<1) {
        return(TRUE);
    }
    /***
    *   Tally counts and positions 
    */
    mat = 0;
    pos = 0;
    for(i=0;i<=n;i++) {
        if(QSMatchI(alPO->cent,&seqPO->seq[i],alPO->centlen,alPO->do_ecc)) {
            mat++;
            pos += i;
        }
    }
    /***
    *   Centroid simply = "average" position
    */
    if(mat>0) {
        alPO->val = DNUM(pos)/DNUM(mat);
        if(!CheckAlphcontSeqmask(alPO,seqPO,FALSE)) {
            return(FALSE);
        } 
        alPO->smask[INT(alPO->val)] = 1;
    }
    else {
        alPO->val = 0.0;
    }
    return(TRUE);
}
/****************************************************************************
*   Check to see if seq mask is big enough;
*   If init_zero is passed, it is also initialized to zero
*/
int CheckAlphcontSeqmask(ALPHCONT *alPO, SEQ *seqPO, int init_zero)
{
    int sl;

    sl = GetSeqLenI(seqPO);
    if ( sl > alPO->n_smask ) {
        if ( ! alPO->smask ) {
            alPO->smask = ALLOC(sl, sizeof(char));
        }
        else {
            alPO->smask = REALLOC(alPO->smask, sl, sizeof(char));
        }
    }
    if ( ! alPO->smask) {
        PROBLINE;
        return(FALSE);
    }
    alPO->n_smask = sl;
    if ( init_zero ) {
        InitArrayI(alPO->smask, IS_CHAR, 0, alPO->n_smask, 0);
    }
    return(TRUE);
}
/****************************************************************************
*   Get max row or simple count of word in sequence
*/
int ScanSeqCountsI(ALPHCONT *alPO, SEQ *seqPO)
{
    int i,j,n,mat,max;
    char *tmaskS = NULL;

    /***
    *   Max row scan
    */
    if(alPO->rowlen>0) {
        n = seqPO->len - alPO->rowlen;
        max = mat = 0;
        for(i=0;i<=n;i++) {
            j = mat = 0;
            while(QSMatchI(alPO->row, &seqPO->seq[i+j], alPO->rowlen,
                alPO->do_ecc) )
            {
                mat++;
                j += alPO->rowlen;
                if( (i+j+alPO->rowlen) > seqPO->len )
                    break;
            }
            max = MAX_NUM(max,mat);
        }
        max = MAX_NUM(max,mat);
        alPO->val = RNUM(max);
    }
    /***
    *   Composition / masking scan
    */
    else {
        n = seqPO->len - alPO->conlen;
        if( ! CheckAlphcontSeqmask(alPO, seqPO, TRUE) ) {
            return(FALSE);
        } 
        mat = 0; 
        for(i=0;i<=n;i++) {
            if ( QSMatchI(alPO->con, &seqPO->seq[i], alPO->conlen, alPO->do_ecc)) {
                mat++;
                for(j=0;j<alPO->conlen;j++) {
                    alPO->smask[i+j] = 1;
                }
            }
        }
        /***
        *   If window sampling sham, do that
        */
        if(alPO->cwin>1) {
            mat = max = 0;
            n = seqPO->len - alPO->cwin;
            for(i=0;i<=n;i++) {
                mat = MaskedInWindowI(alPO->smask, i, alPO->cwin);
                max = MAX_NUM(max, mat);
            }
            max = MAX_NUM(max, mat);
            alPO->val = RNUM(max) / RNUM(alPO->conlen);
        }
        else {
            alPO->val = RNUM(mat);
        }
    }
    /***
    *   Some output options need second pass for rows / cwin 
    *   Allocate temp mask here as well
    */
    n = FALSE;
    switch(alPO->owhat)
    {
        case ALCO_DCC:
        case ALCO_DFS: 
        case ALCO_DFG: 
            n++;
    }
    if(n) {
        if( ! CheckAlphcontSeqmask(alPO, seqPO, FALSE) ) {
            return(FALSE);
        } 
        if( ! (tmaskS = ALLOC(GetSeqLenI(seqPO), sizeof(char))) ) {
            return(FALSE);
        }
        /***
        *   Rows
        */
        if(alPO->rowlen>0) {
            /***
            *   If any hits, check each starting position against max 
            */
            if(INT(alPO->val) > 0) {
                n = seqPO->len - alPO->rowlen;
                for(i=0;i<=n;i++) {
                    j = mat = 0;
                    while(QSMatchI(alPO->row, &seqPO->seq[i+j], alPO->rowlen,
                        alPO->do_ecc) )
                    {
                        mat++;
                        j += alPO->rowlen;
                        if( (i+j+alPO->rowlen) > seqPO->len )
                            break;
                    }
                    if(mat == INT(alPO->val)) {
                        mat *= alPO->rowlen;
                        for(j=0;j<mat;j++) {
                            alPO->smask[i+j] = 1;
                        }
                    }
                }
            }
        }
        /***
        *   Window sampling
        */
        else if(alPO->cwin>1) {
            CleanString(tmaskS,seqPO->len);
            /***
            *   Check each starting position against max 
            */
            if(INT(alPO->val) > 0) {
                n = seqPO->len - alPO->cwin;
                for(i=0;i<=n;i++) {
                    mat = MaskedInWindowI(alPO->smask,i,alPO->cwin);
                    mat /= alPO->conlen;
                    if(mat == INT(alPO->val)) {
                        for(j=0;j<alPO->cwin;j++) {
                            tmaskS[i+j] = 1;
                        }
                    }
                }
            }
            for(j=0; j< seqPO->len; j++) {
                alPO->smask[j] = tmaskS[j];
            }
        }
    }
    CHECK_FREE(tmaskS);
    return(TRUE);
}
/**************************************************************************
*   How many matches in window?
*/
int HandleFractionValueI(ALPHCONT *alPO, SEQ *seqPO)
{
    int i,n,numn;

    /***
    *   No fraction if windows
    */
    if ( alPO->cwin > 1 ) {
        return(TRUE);
    } 
    /***
    *   If ingorning N's count them
    */
    numn = 0;
    if(alPO->do_ign) {
        for(i=0;i<seqPO->len;i++) {
            if( (seqPO->seq[i]=='N') || (seqPO->seq[i]=='n') ) {
                numn++;
            }
        }
    }
    /***
    *   If reporting fractional value, adjust here
    */
    if(alPO->do_dfr) {
        n = seqPO->len - alPO->conlen + 1;
        n -= numn;
        if(n>0) {
            alPO->val /=  RNUM(n);
        }
        else {
            alPO->val = 0.0;
        }
    }
    return(TRUE);
}
/**************************************************************************
*   How many matches in window?
*/
int MaskedInWindowI(char *maskS,int s,int n)
{
    int j,mat;

    mat = 0;
    for(j=0;j<n;j++)
    {
        if(maskS[s+j])
            mat++;
    }
    return(mat);
}
/**************************************************************************
*   Check current seq values against score constraints
*/
int IsAlphcontSeqOkI(ALPHCONT *alPO)
{
    int ok;

    ok = TRUE;
    if( (alPO->val<alPO->min) || (alPO->val>alPO->max) )
    {
        ok = FALSE;
    }
    /***
    *   If not, invert qualification
    */
    if(alPO->do_not)
    {
        ok = !ok;
    }
    return(ok);
}
/**************************************************************************
*   Handle output for current seq
*/
int HandleAlcoOutputI(ALPHCONT *alPO, SEQ *seqPO, char *nameS, int sok,
    FILE *outPF)
{
    int slen,ok,n,i,mat;
    char forS[DEF_BS],*seqPC, *fseqPC;

    HAND_NFILE(outPF);
    /***
    *   alPO->seq = full sequence, in case seqPO has been trimmed
    */
    GetSeqSeqI(seqPO,&seqPC);
    slen = GetSeqLenI(seqPO);
    GetSeqSeqI(alPO->fseq,&fseqPC);
    /***
    *   If simple base table, handle and return
    */
    if( alPO->do_btab || alPO->do_dtab || alPO->do_rtab ) {
        ReportAlphcontBaseTables(alPO,nameS,seqPO,alPO->comp,outPF);
        return(TRUE);
    }
    /***
    *   If flagging seqs, this goes in front of other outputs
    */
    INIT_S(forS);
    if(alPO->do_qual) {
        if(alPO->owhat!=ALCO_SEQ) {
            if(sok)
            {   sprintf(forS,"1 ");     }
            else
            {   sprintf(forS,"0 ");     }
        }
    }
    /***
    *   What remaining output
    */
    switch(alPO->owhat)
    {
        case ALCO_VAL:
            fprintf(outPF,"%s%-15s\t",forS,nameS);
            if(alPO->do_dfr) {
                fprintf(outPF,"%6.3f",alPO->val);
                if(alPO->skewlen>0) {
                    fprintf(outPF,"\t%6.3f",alPO->val2);
                    /***
                    *   If dumping skew, also add length and ABS of both
                    */
                    fprintf(outPF,"\t%2d",slen);
                    fprintf(outPF,"\t%6.3f",ABS_VAL(alPO->val));
                    fprintf(outPF,"\t%6.3f",ABS_VAL(alPO->val2));
                }
            }
            else {
                fprintf(outPF,"%2d",INT(alPO->val));
            }
            if(alPO->do_ds) {
                fprintf(outPF,"\t%s",fseqPC);
            }
            fprintf(outPF,"\n");
            break;
        case ALCO_SEQ:
            ok = TRUE;
            if( (alPO->val<alPO->min) || (alPO->val>alPO->max) ) {   
                ok = FALSE; 
            }
            if(alPO->do_not) {   
                ok = !ok;   
            }
            if(ok) {
                fprintf(outPF,"# %-15s\t",nameS);
                if(alPO->do_dfr) {
                    fprintf(outPF,"%6.3f",alPO->val);
                    if(alPO->skewlen>0) {
                        fprintf(outPF,"\t%6.3f",alPO->val2);
                    }
                    fprintf(outPF,"\n");
                }
                else {
                    fprintf(outPF,"%2d\n",INT(alPO->val));
                }
                fprintf(outPF,"%s%-15s\t%s\n\n",forS,nameS,fseqPC);
            }
            break;
        case ALCO_DCC:
            fprintf(outPF,"%-15s \t",nameS);
            n = 0;
            for(i=0;i<slen;i++)
            {
                if(alPO->smask[i]) {   
                    n++;    
                }
                fprintf(outPF,"%d ",n);
            }
            fprintf(outPF,"\n");
            break;
        case ALCO_DFG:
            fprintf(outPF,"%-15s\t",nameS);
            DumpAlphconMask(alPO, slen, outPF);
            break;
        case ALCO_DFS:
            fprintf(outPF,"\n");
            fprintf(outPF,"# %-15s\t",nameS);
            DumpAlphconMask(alPO, slen, outPF);
            fprintf(outPF,"%-15s  \t%s\n",nameS,seqPC);
            break;
        case ALCO_DCW:
            fprintf(outPF,"%s%-15s\t",forS,nameS);
            n = slen - alPO->cwin;
            for(i=0;i<=n;i++)
            {
                mat = MaskedInWindowI(alPO->smask,i,alPO->cwin);
                fprintf(outPF," %d",mat);
            }
            fprintf(outPF,"\n");
            break;
    }
    return(TRUE);
}
/**************************************************************************
*   Dump X and . based on mask to passed output file
*/
void DumpAlphconMask(ALPHCONT *alPO, int len, FILE *outPF)
{
    int i;
    
    HAND_NFILE(outPF);
    for(i=0;i<len;i++)
    {
        if(alPO->smask[i])
        {   fprintf(outPF,"X"); }
        else
        {   fprintf(outPF,"."); }
    }
    fprintf(outPF,"\n");
}
/**************************************************************************
*   Dump table output
*/
void ReportAlphcontBaseTables(ALPHCONT *alPO, char *nameS, SEQ *seqPO, SEQCOMP *compPO, FILE *outPF)
{
    char *fseqPC;

    VALIDATE(compPO,SEQCOMP_ID);
    HAND_NFILE(outPF);

    fprintf(outPF,"%s",nameS);
    if(alPO->do_dtab) {
        ReportAlphcontDinucCounts(compPO, alPO->do_dfr, outPF);
        return;
    }
    fprintf(outPF,"\t");
    /***
    *   Fraction bases 
    */
    if(alPO->do_dfr) {
        fprintf(outPF,"%6.3f\t%6.3f\t%6.3f\t%6.3f",
            compPO->fa, compPO->fc, compPO->fg, compPO->ft);
        fprintf(outPF,"\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f",
            (compPO->fg + compPO->fc), (compPO->fa + compPO->ft),
            (compPO->fg + compPO->fa), (compPO->fc + compPO->ft),
            (compPO->fa + compPO->fc), (compPO->fg + compPO->ft) );
    }
    /***
    *   Rows or counts 
    */
    else {
        if(alPO->do_rtab) {
            fprintf(outPF,"%3d\t%3d\t%3d\t%3d\t%3d",
                compPO->ra, compPO->rc, compPO->rg, compPO->rt, compPO->rmax1);
            fprintf(outPF,"\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d",
                compPO->rs, compPO->rs, compPO->rr, compPO->ry, compPO->rm, compPO->rk, compPO->rmax2);
        }
        else {
            fprintf(outPF,"%3d\t%3d\t%3d\t%3d",
                compPO->na, compPO->nc, compPO->ng, compPO->nt);
            fprintf(outPF,"\t%3d\t%3d\t%3d\t%3d\t%3d\t%d",
                (compPO->ng + compPO->nc), (compPO->na + compPO->nt), 
                (compPO->ng + compPO->na), (compPO->nc + compPO->nt), 
                (compPO->na + compPO->nc), (compPO->ng + compPO->nt) );
        }
    }
    GetSeqSeqI(alPO->fseq,&fseqPC);
    if(alPO->do_ds) {
        fprintf(outPF,"\t%s",fseqPC);
    }
    fprintf(outPF,"\n");
}
/**************************************************************************
*
*/
void ReportAlphcontDinucCounts(SEQCOMP *compPO, int frac, FILE *outPF)
{
    int i,n;

    VALIDATE(compPO,SEQCOMP_ID);
    HAND_NFILE(outPF);
    n = compPO->n_dinuc;
    if (n < 1) {
        fprintf(outPF," TOO SHORT, no dinucs\n");
        return;
    }
    for(i=0;i<16;i++) {
        if (frac) {
            fprintf(outPF,"\t%6.3f",RNUM(compPO->dinuc[i])/RNUM(n));
        }
        else {
            fprintf(outPF,"\t%d",compPO->dinuc[i]);
        }
    }
    fprintf(outPF,"\n");
}
/**************************************************************************
*   String comparision function; This is the main work function
*/
int QSMatchI(char *qS,char *sS,int len,int exp)
{
    int i,ok;
    
    DB_ALFC 
    {
        DB_PrI(">> QSMatchI |%s| |",qS); 
        PrintString(sS,len,NULL); 
        DB_PrI("| exp=%d\n",exp);
    }
    /***
    *   If explicit, then only true stingmatch is a match
    */
    if(exp) {
        if(!strncasecmp(qS,sS,len)) {
            DB_ALFC DB_PrI("<< QSMatchI exp TRUE\n"); 
            return(TRUE);
        }
        DB_ALFC DB_PrI("<< QSMatchI exp FALSE\n"); 
        return(FALSE);
    }
    /***
    *   Letter-by-letter comparison
    */
    ok = TRUE;
    for(i=0;i<len;i++)
    {
        switch(sS[i])
        {
            case 'A': case 'a':
                switch(qS[i])
                {
                    case 'A': case 'a':
                    case 'M': case 'm':
                    case 'R': case 'r':
                    case 'W': case 'w':
                    case 'D': case 'd':
                    case 'H': case 'h': 
                    case 'V': case 'v':
                    case 'N': case 'n':
                        break;
                    default:
                        ok = FALSE;
                }
                break;
            case 'C': case 'c':
                switch(qS[i])
                {
                    case 'C': case 'c':
                    case 'M': case 'm':
                    case 'S': case 's':
                    case 'Y': case 'y':
                    case 'B': case 'b':
                    case 'H': case 'h':
                    case 'V': case 'v':
                    case 'N': case 'n':
                        break;
                    default:
                        ok = FALSE;
                }
                break;
            case 'G': case 'g':
                switch(qS[i])
                {
                    case 'G': case 'g':
                    case 'R': case 'r':
                    case 'S': case 's':
                    case 'K': case 'k':
                    case 'B': case 'b':
                    case 'D': case 'd':
                    case 'V': case 'v':
                    case 'N': case 'n':
                        break;
                    default:
                        ok = FALSE;
                }
                break;
            case 'T': case 't':
                switch(qS[i])
                {
                    case 'T': case 't':
                    case 'W': case 'w':
                    case 'Y': case 'y':
                    case 'K': case 'k':
                    case 'B': case 'b':
                    case 'D': case 'd':
                    case 'H': case 'h':
                    case 'N': case 'n':
                        break;
                    default:
                        ok = FALSE;
                }
                break;
            case 'N': case 'n':
                switch(qS[i])
                {
                    case 'N': case 'n':
                        break;
                    default:
                        ok = FALSE;
                }
                break;
            default:
                ok = FALSE;
        }
        if(!ok) {
            break;
        }
    }
    DB_ALFC DB_PrI("<< QSMatchI %d\n",ok); 
    return(ok);
}
