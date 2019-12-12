/*
* dna_util.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
*
* The programs and source code of the vertools collection are free software.
* They are distributed in the hope that they will be useful,
* WITHOUT ANY WARRANTY OF FITNESS FOR ANY PARTICULAR PURPOSE.  
* 
* Permission is granted for research, educational, and possibly commercial use 
*   and modification as long as 1) Code and any derived works are not 
*   redistributed for any fee, and 2) Proper credit is given to the authors. 
*   If you wish to include this software in a product, or use it commercially,
*   please contact the authors.
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
#include "seq_info.h"
#include "wordlist.h"
#include "dna_util.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(DnaUtilI(argc,argv),NULL) ); }
/**************************************************************************/
void DnaUtilUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <infile>   Sequences file to read in\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -out XXX   Set output to XXX\n");
    printf("   -sran # #  Sequence range restricted # to # (1-base count)\n");
    printf("   -com       Reverse Complement: e.g. CCCAT >-> ATGGG \n");
    printf("   -inv       Inverse Compliment: e.g. CCCAT >-> GGGTA\n");
    printf("   -rev       Reverse direction;  e.g. CCCAT >-> TACCC\n");
    printf("   -nan       No automatic added naming for -com -inv -rev\n");
    printf("   -cln -cll  Clean non-ACGT / lowercase chars to N\n");
    printf("   -isc       Invert sequence case; e.g. AcGt >-> aCgT\n");
    printf("   -cls -exi  Clean SNP [X/Y] to IUB codes / Expand to [X/Y]\n");
    printf("   -cstat     Case-base stat report (UPPER lower ambigN)\n");
    printf("   -csep      Case-base separation of subseqs\n");
    printf("   -tnb       Truncate name on blank (first token only)\n");
    printf("   -b2u       Replace blanks in name with underscores\n");
    printf("   -oraw      Output \"raw\" format\n");
    printf("   -ofas      Output fasta format\n");
    printf("   -nfas      Output 'nice' fasta format, lines and blocks\n");
    printf("   -nfl #     Nice fasta output in lines of #\n");
    printf("   -nfb #     Nice fasta output in blocks of #\n");
    printf("   -olis      Output list of names\n");
    printf("   -flg       Flag output status; i.e. 1=in, 0=out\n");
    printf("   -ostat     Output sequences stats (Length Amb N LC [SNP])\n");
    printf("   -stat      Report overall stats for input file\n");
    printf("   -ds -di    Dump modified / input sequences appended as last column\n");
    printf("   -bran # #  Base range # to # (1-base coords)\n");
    printf("   -braf # #  Fractional base range # to # (0-1)\n");
    printf("   -rre       Base range relative to end; i.e. backwards\n");
    printf("   -trim # #  Trim # and # from ends\n");
    printf("   -win #     Window of # bases; Default centered\n");
    printf("   -wif #     Fraction window of # (0-1)\n");
    printf("   -wcp #     Window center postion\n");
    printf("   -wcf #     Fractional window center postion\n");
    printf("   -msk -imsk Mask / Invert mask selected range with N\n");
    printf("   -muc -mlc  Mark selected range Upper / Lower case\n");
    printf("   -wlis XXX  Qualify with words (first token) listed in XXX\n");
    printf("   -kc        Keep case for token comparison (default ignore)\n");
    printf("   -wst       Word start only needs to match line (not full token)\n");
    printf("   -wsub      Word substring only needs to match line (not full token)\n");
    printf("   -flen # #  Flag sequences with lengths of # to #\n");
    printf("   -famb -fsnp Flag sequences with ambigous bases / SNP sites\n");
    printf("   -not       Flagging logical not; out>-->in & in>-->out\n");
    printf("   -pol #     Probes of length #\n");
    printf("   -pml # #   Probes multi length # to #\n");
    printf("   -psj # #   Probes starting at # and jumping #\n");
    printf("   -plp XXX   Probes listed in XXX; <len> <start> / line\n");
    printf("   -pew # #   Probes \"extra window\" bases up/down\n");
    printf("   -pat #     Probes allowed to be truncated to len # (seq end)\n");
    printf("   -pco       Probes \"clean\" only; no ambigs or SNPs\n");
    printf("   -inwf # #  Info-index based on word sizes # to #\n");
    printf("   -insh -inbr -inbp\tInfo-index: Shannon / Brillouin / Berger-Parker\n");
    printf("   -igp       Ignore non-fatal problems and just keep going!\n");
    printf("\n");
}
/**************************************************************************
*   Main program
*/
int DnaUtilI(int argc, char **argv)
{
    int n,nok,ok,oraw,ofas,nfas,olis,ostat,stat,iraw,iseq,ifas,do_cln,do_cls,do_cll;
    DNA_UTIL *duPO;
    SEQTRIM *stPO;

    duPO = CreateDna_utilPO();
    stPO = duPO->subseq;
    oraw = ofas = nfas = ostat = stat = olis = iraw = ifas = iseq = FALSE;
    do_cln = do_cls = do_cll = FALSE;
    if(!ParseArgsI(argc,argv,
        "S -out S -oraw B -com B -inv B -rev B -ofas B -nfas B\
        -stat B -flen I2 -famb B -not B\
        -bran I2 -braf D2 -rre B -trim I2\
        -win I -wif D -wcp I -wcf D\
        -msk B -imsk B -muc B -mlc B\
        -sran I2 -wlis S -kc B -olis B -fsnp B\
        -iraw B -ifas B -ostat B -plp S -pew I2 -cln B -cls B -pol I -igp B\
        -flg B -psj I2 -inwf I2 -insh B -inbr B -inbp B\
        -cll B -pat I -exi B -tnb B -b2u B\
        -pco B -pml I2\
        -cstat B -csep B -wst B -tsub B -nan B -nfl I -nfb I -ds B -di B\
        -iseq B -isc B",
        duPO->inname, duPO->outname, &oraw, &duPO->do_comp, &duPO->do_inv, 
        &duPO->do_rev, &ofas, &nfas, &stat, &duPO->min_len,&duPO->max_len, 
        &duPO->do_famb, &duPO->do_not,  
        &stPO->base_s,&stPO->base_e, &stPO->basf_s,&stPO->basf_e, &stPO->rre, 
        &stPO->trim_s,&stPO->trim_e, 
        &stPO->wind, &stPO->winf, &stPO->wcent, &stPO->wcenf,
        &stPO->nmask, &stPO->umask, &stPO->ucase, &stPO->lcase,
        &duPO->first,&duPO->last, 
        duPO->wlisname, &duPO->do_kc, &olis, 
        &duPO->do_fsnp, &iraw, &ifas, &ostat, duPO->olpname,
        &duPO->olp_up,&duPO->olp_dn, &do_cln, &do_cls, &duPO->do_pol,
        &duPO->igprob, &duPO->do_flg, &duPO->opl_st,&duPO->opl_j,
        &duPO->ifwmin,&duPO->ifwmax,
        &duPO->do_insh, &duPO->do_inbr, &duPO->do_inbp,
        &do_cll, &duPO->do_pat, &duPO->do_exi,
        &duPO->do_tnb, &duPO->do_b2u, &duPO->do_pco, &duPO->do_pol,&duPO->do_pml,
        &duPO->do_cstat, &duPO->do_csep, 
        &duPO->do_wst, &duPO->do_wsub, &duPO->do_nan,
        &duPO->nfline, &duPO->nfblock, &duPO->do_ds, &duPO->do_di,
        &iseq, &duPO->do_isc,
        (int *)NULL))
    {
        DnaUtilUse();
        CHECK_DNA_UTIL(duPO);
        return(FALSE);
    }
    /***
    *   Set input format 
    */
    duPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas, duPO->inname, TRUE);
    if(!duPO->iform) {
        printf("Problem with input seq(s)\n");
        CHECK_DNA_UTIL(duPO);
        return(FALSE);
    }
    /***
    *   What output options? (precedence increasing down list)
    */
    if(stat) {
        duPO->owhat = DNUO_OSTAT;
    }
    else if(olis) {
        duPO->owhat = DNUO_LIS;
    }
    else if(ostat) {
        duPO->owhat = DNUO_STAT;
    }
    else if(AnyInfoSettingsI(duPO)) {
        duPO->owhat = DNUO_INFO;
    }
    else if( (duPO->do_cstat) || (duPO->do_csep) ) {
        duPO->owhat = DNUO_CSEP;
    }
    else {
        duPO->owhat = DNUO_SEQ;
    }
    /***
    *   Output formats?
    */
    duPO->oform = SEQFM_RAW;
    if(oraw) {
        duPO->oform = SEQFM_RAW;
    }
    else if(ofas) {
        duPO->oform = SEQFM_FASTA;
    }
    else if(nfas) {
        duPO->oform = SEQFM_NFAS;
    }
    /***
    *   Sequence "cleaning" options
    */
    if(duPO->do_exi) {
        duPO->do_clean = 0;
    }
    else if(do_cls) {
        duPO->do_clean = CLEAN_SNP;
    }
    else if(do_cll) {
        duPO->do_clean = CLEAN_CASE;
    }
    else if(do_cln) {
        duPO->do_clean = CLEAN_ALL;
    }
    if(!CheckDnuOptionsI(duPO)) {
        ABORTLINE;
        CHECK_DNA_UTIL(duPO);
        return(FALSE);
    }
    /***
    *   Open in/out files
    */
    if(!OpenDnuFilesI(duPO)) {
        ABORTLINE;
        CHECK_DNA_UTIL(duPO);
        return(FALSE);
    }
    /***
    *   Output header 
    */
    SetUpOptionFlags(duPO);
    HandleDnuFlagHeader(duPO,duPO->out);
    HandleDnuInfoHeader(duPO,duPO->out);
    HandleDnuOstatHeader(duPO,duPO->out);
    /***
    *   Loop through the input collection
    */
    n = nok = 0;
    while(TRUE) {
        /***
        *   Parse sequence; FALSE = done
        */
        ok = ParseSeqI(duPO->in, duPO->iform, n+1, duPO->iclean, TRUE, duPO->fseq);
        if(ok==FALSE) {
            break;
        }
        if(ok!=TRUE) {
            if(!duPO->igprob) {
                ABORTLINE;
                break;
            }
            continue;
        }
        /***
        *   Filter here BEFORE any manipulations
        *   Update count and qualify
        */
        n++;
        ok = IsCurrentSeqOkI(duPO, duPO->fseq, n);
        if(ok) {
            nok++;
            HandleDuCleanExpI(duPO, duPO->fseq);
            /* Copy full to working seq */
            CopySeqI(duPO->fseq, duPO->seq, -1, -1);
            /***
            *   Do any masking, trim subseqs, sequence flips
            */
            if(!DnaTrimSeqI(duPO->subseq, duPO->seq)) {
                ABORTLINE;
                break;
            }
            HandleDuSeqFlipsI(duPO, duPO->seq);
            if(duPO->do_isc) {
                InvertSeqCase(duPO->seq);
            }
            /***
            *   If outputting listed / to-length probes, handle here
            */
            if (duPO->owhat==DNUO_PROBE) {
                if(!HandleDuProbesOutI(duPO, duPO->seq, duPO->out)) {
                    if(!duPO->igprob) {
                        ABORTLINE;
                        break;
                    }   
                }
                continue;
            }
        }
        /***
        *   Output story
        */
        HandleDuOutputI(duPO, duPO->seq, ok, duPO->fseq, duPO->out);
    }
    /***
    *   Stats?
    */
    if (duPO->owhat==DNUO_OSTAT) {
        HandleDuStats(duPO,n,nok,duPO->out);
    }
    /***
    *   All done
    */
    CHECK_DNA_UTIL(duPO);
    return(TRUE);
}
/*****************************************************************************
*   Create data struct
*/
DNA_UTIL *CreateDna_utilPO()
{
    DNA_UTIL *duPO;

    if(! (duPO = (DNA_UTIL *)ALLOC(1,sizeof(DNA_UTIL)) ) )
    {
        printf("# Failed to allocate working object\n");
        return(NULL);
    }
    duPO->ID = DNA_UTIL_ID;
    duPO->seq = CreateSeqPO(0,NULL,NULL);
    duPO->fseq = CreateSeqPO(0,NULL,NULL);
    duPO->subseq = CreateSeqtrimPO();
    InitDna_util(duPO);
    return(duPO);
}
/*****************************************************************************
*   Free datastructure and substructs
*/
int DestroyDna_utilI(DNA_UTIL *duPO)
{
    VALIDATE(duPO,DNA_UTIL_ID);
    CHECK_FILE(duPO->in);
    CHECK_WORDLIST(duPO->wlis);
    CHECK_FILE(duPO->olp);
    CHECK_NFILE(duPO->out,duPO->outname);
    CHECK_SEQ(duPO->seq);
    CHECK_SEQ(duPO->fseq);
    CHECK_SEQTRIM(duPO->subseq);
    FREE(duPO);
    return(TRUE);
}
/*****************************************************************************
*   Set null / default values
*/
void InitDna_util(DNA_UTIL *duPO)
{
    VALIDATE(duPO,DNA_UTIL_ID);

    INIT_S(duPO->inname);
    INIT_S(duPO->outname);
    duPO->in = NULL;
    duPO->iform = BOGUS;
    duPO->iclean = SCLEAN_MID;
    duPO->out = NULL;
    duPO->owhat = BOGUS;
    duPO->nfblock = duPO->nfline = BOGUS;
    duPO->igprob = FALSE;
    duPO->do_kc = duPO->do_wst = duPO->do_wsub = FALSE;
    duPO->do_tnb = duPO->do_b2u = FALSE;
    duPO->do_comp = duPO->do_rev = duPO->do_inv = FALSE;
    duPO->do_nan = FALSE;
    duPO->do_not = FALSE;
    duPO->do_clean = FALSE;
    duPO->do_exi = FALSE;
    duPO->do_cstat = FALSE;
    duPO->do_csep = FALSE;
    duPO->do_flg = FALSE;
    duPO->do_famb = FALSE;
    duPO->do_fsnp = FALSE;
    duPO->first = BOGUS;
    duPO->last = BOGUS;
    duPO->n_flags = 0;
    duPO->min_len = BOGUS; 
    duPO->max_len = BOGUS;
    duPO->len_hi = duPO->snp_hi = -TOO_BIG;
    duPO->len_lo = duPO->snp_lo = TOO_BIG;
    duPO->snp_c = duPO->amb_c = 0;
    duPO->do_pol = duPO->do_pml = 0;
    duPO->do_pat = 0;
    duPO->do_pco = FALSE;
    duPO->opl_st = 0;
    duPO->opl_j = 0;
    INIT_S(duPO->olpname);
    duPO->olp = NULL;
    duPO->olp_up = duPO->olp_dn = 0;
    duPO->do_insh = duPO->do_inbr = duPO->do_inbp = FALSE;
    duPO->ifwmin = duPO->ifwmax = 0;
    duPO->do_ds = duPO->do_di = FALSE;
}
/**************************************************************************
*   Information indices set?
*/
int AnyInfoSettingsI(DNA_UTIL *duPO)
{
    if( duPO->do_insh || duPO->do_inbr || duPO->do_inbp ) {
        return(TRUE);
    }
    if( (duPO->ifwmin > 0) || (duPO->ifwmax > 0 ) ) {
        return(TRUE);
    }
    return(FALSE);
}
/*************************************************************************
*   Check for (some) option consistency
*/
int CheckDnuOptionsI(DNA_UTIL *duPO)
{
    if(!CheckDnuProbeOptionsI(duPO)) {
        return(TRUE);
    }
    if( (duPO->ifwmin>0) || (duPO->ifwmax>0 ) ) {
        if( duPO->ifwmin > duPO->ifwmax ) {
            PROBLINE;
            printf("Bad word sizes for -inwf: %d %d\n",
                duPO->ifwmin, duPO->ifwmax);
            return(FALSE);
        }
    }
    return(TRUE);
}
/***************************************************************************
*   Open any needed files or die
*/
int OpenDnuFilesI(DNA_UTIL *duPO)
{
    if(!(duPO->in=OpenUFilePF(duPO->inname,"r",NULL))) {
        return(FALSE);
    }
    if(!NO_S(duPO->wlisname)) {
        if( ! (duPO->wlis = CreateWordlistPO(duPO->wlisname,NSIZE))) {
            PROBLINE;
            printf("Failed to get tokens from %s\n",duPO->wlisname);
            return(FALSE);
        }
    }
    if(!NO_S(duPO->olpname)) {
        if(!(duPO->olp =OpenUFilePF(duPO->olpname,"r",NULL))) {
            return(FALSE);
        }
    }
    if(!NO_S(duPO->outname)) {
        if(!(duPO->out=OpenUFilePF(duPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    HAND_NFILE(duPO->out);
    return(TRUE);
}
/***************************************************************************
*   Check if any filters are active (to be reported)
*/
int SetUpOptionFlags(DNA_UTIL *duPO)
{
    int n;

    n = 0;
    if(IS_BOG(duPO->first)) {
        duPO->first = duPO->last = -1;
    }
    else {
        n++;
    }
    if(IS_BOG(duPO->min_len)) {
        duPO->min_len = duPO->max_len = -1;
    }
    else {
        n++;
    }
    if( duPO->do_famb ) {
        n++;
    }
    if( duPO->do_fsnp ) {
        n++;
    }
    if(!NO_S(duPO->wlisname)) {
        n++;
    }
    duPO->n_flags = n;
    return(n);
}
/***************************************************************************
*   Flip around seq and maybe rename them
*/
int HandleDuSeqFlipsI(DNA_UTIL *duPO, SEQ *seqPO)
{
    char *seqPC, nameS[NSIZE], addS[DEF_BS];
    int len, any;

    seqPC = NULL;
    any = FALSE;
    INIT_S(addS);
    len = 0;    /* make compiler happy */
    if(seqPO) {
        GetSeqSeqI(seqPO, &seqPC);
        len = GetSeqLenI(seqPO);
    }
    /*** 
    *   Options that involve flipping / tweaking 
    */
    if(duPO->do_comp) {
        any++;
        if(seqPO) {
            CompDNASeqI(seqPC,len,seqPC);
            sprintf(addS,"_rc");
        }
    }
    else if(duPO->do_inv) {
        any++;
        if(seqPO) {
            InverseDNASeqI(seqPC,len,seqPC);
            sprintf(addS,"_inv");
        }
    }
    else if(duPO->do_rev) {
        any++;
        if(seqPO) {
            ReverseDNASeqI(seqPC,len,seqPC);
            sprintf(addS,"_rev");
        }
    }
    if( (seqPO) && (any) && (!duPO->do_nan) ) {
        FillSeqNameStringI(seqPO,nameS,-1);
        strcat(nameS,addS);
        SetSeqName(seqPO,nameS);
    }
    return(any);
}
/**************************************************************************
*   Screen current seq against filters
*/
int IsCurrentSeqOkI(DNA_UTIL *duPO, SEQ *seqPO, int n)
{
    int ok;
    char nameS[NSIZE];

    ok = TRUE;
    if(duPO->first > 0) {
        if( (n < duPO->first) || (n > duPO->last) ) {
            ok = FALSE;
        }
    }
    if(duPO->min_len > 0) {
        if( (seqPO->len < duPO->min_len) || (seqPO->len > duPO->max_len) ) {
            ok = FALSE;
        }
    }
    if( duPO->do_famb && IS_SEQ_AMB(seqPO->flag) ) {
        ok = FALSE;
    }
    if( duPO->do_fsnp && IS_SEQ_SNP(seqPO->flag) ) {
        ok = FALSE;
    }
    if( ok && duPO->wlis ) {
        FillSeqNameStringI(seqPO,nameS,-1);
        ok = WordInWordlistI(duPO->wlis, nameS, duPO->do_kc, duPO->do_wst, 
                            duPO->do_wsub, NULL);
    } 
    /***
    *   If not, invert qualification
    */
    if(duPO->do_not) {
        ok = !ok;
    }
    return(ok);
}
/**************************************************************************
*   Merge seqPO with current seq 
*/
int HandleMergeSeqsI(SEQ *fPO, SEQ *sPO, int dir, int over, SEQ **newPPO)
{
    int nlen;
    SEQ *newPO;
    char *fPC, *sPC, *nPC, nameS[NSIZE];

    /***
    *   Create new seq big enough for both parts
    */
    nlen = GetSeqLenI(fPO) + GetSeqLenI(sPO);;
    if(!(newPO=CreateSeqPO(nlen,NULL,NULL)) ) {
        return(FALSE);
    }
    if( (!GetSeqSeqI(fPO,&fPC)) || (!GetSeqSeqI(sPO,&sPC)) || 
        (!GetSeqSeqI(newPO,&nPC)) ) {
        return(FALSE);
    }
    if(!SpliceTwoSeqsI(fPC,sPC,NULL,dir,over,FALSE,NULL,nPC)) {
        CHECK_SEQ(newPO);
        return(FALSE);
    }
    /***
    *   Finish new guy, copying name from start
    */
    SetSeqSequenceI(newPO,nPC,strlen(nPC));
    FillSeqNameStringI(fPO,nameS,NSIZE-1);
    SetSeqName(newPO,nameS);
    *newPPO = newPO;
    return(TRUE);
}
/**************************************************************************
*   Sequence "cleaning" options
*/
int HandleDuCleanExpI(DNA_UTIL *duPO, SEQ *seqPO)
{
    int len, slen;
    char nameS[NSIZE],tokS[NSIZE], *seqPC;

    VALIDATE(seqPO,SEQ_ID);
    GetSeqSeqI(seqPO, &seqPC);
    len = GetSeqLenI(seqPO);
    /***
    *   Sequence per se
    */
    switch(duPO->do_clean)
    {
        case CLEAN_CASE:
            slen = CleanUpSeqI(seqPC, len, seqPC, FALSE,TRUE);
            seqPO->len = slen;
            break;
        case CLEAN_ALL:
            slen = CleanUpSeqI(seqPC, len, seqPC, FALSE,FALSE);
            seqPO->len = slen;
            break;
        case CLEAN_SNP:
            slen = CleanUpSeqI(seqPC, len, seqPC, TRUE,FALSE);
            seqPO->len = slen;
            break;
    }
    if(duPO->do_exi) {
        ExpandSeqSingBaseSNPsI(seqPO);
    }
    /***
    *   Name first token only
    */
    if(duPO->do_tnb) {
        FillSeqNameStringI(seqPO, nameS, NSIZE-1);
        sscanf(nameS,"%s",tokS);
        SetSeqName(seqPO,tokS); 
    }
    else if(duPO->do_b2u) {
        FillSeqNameStringI(seqPO, nameS, NSIZE-1);
        ReplaceChars(' ', nameS, '_', nameS);
        SetSeqName(seqPO,nameS); 
    }
    return(TRUE);
}
/**************************************************************************
*   Handle output for current seq
*/
int HandleDuOutputI(DNA_UTIL *duPO, SEQ *seqPO, int ok, SEQ *fseqPO, FILE *outPF)
{
    int slen,snps;
    char nameS[NSIZE],*seqPC;

    HAND_NFILE(outPF);
    FillSeqNameStringI(seqPO, nameS, NSIZE-1);
    GetSeqSeqI(seqPO, &seqPC);
    slen = GetSeqLenI(seqPO);
    snps = GetSeqSnpCountI(seqPO);
    switch(duPO->owhat)
    {
        case DNUO_SEQ:  
            if(duPO->do_flg) {
                if( duPO->oform == SEQFM_FASTA) {
                    fprintf(outPF,"# ");
                }
                if(ok) {
                    fprintf(outPF,"1 ");
                }
                else {
                    fprintf(outPF,"0 ");
                }
                WriteSeq(seqPO,duPO->oform,outPF); 
            }
            else if(ok) {
                if(slen < 1) {
                    fprintf(outPF,"# %s zero length\n",nameS);
                }
                else {
                    WriteSeq(seqPO, duPO->oform, outPF); 
                }
            }
            break;
        case DNUO_LIS:
            if(duPO->do_flg) {
                if(ok) {
                    fprintf(outPF,"1 ");
                }
                else {
                    fprintf(outPF,"0 ");
                }
                fprintf(outPF,"%s\n",nameS);
            }
            else if(ok) {
                fprintf(outPF,"%s\n",nameS);
            }
            break;
        /***
        *   Length, number of ambigs, N's, lowercase, SNPs
        */
        case DNUO_STAT:
            if(ok) {
                HandleDnaOstatI(duPO, nameS, seqPO, fseqPO, outPF);
            }
            break;
        case DNUO_INFO:
            if(ok) {
                if(snps > 0) {
                    fprintf(outPF,"%s\tSNPS so NO INFO\n",nameS);
                }
                else {
                    HandleDnaInfoOutI(duPO, duPO->seq, duPO->fseq, outPF);
                }
            }
            break;
        case DNUO_CSEP:
            if(ok) {
                HandleCaSepOutputI(duPO, seqPO, outPF);
            }
            break;
        /***
        *   If not some explicit output, tally stats for final report
        */
        default:
            duPO->len_hi = MAX_NUM(slen,duPO->len_hi);
            duPO->len_lo = MIN_NUM(slen,duPO->len_lo);
            duPO->snp_hi = MAX_NUM(snps,duPO->snp_hi);
            duPO->snp_lo = MIN_NUM(snps,duPO->snp_lo);
            if(snps>0) {
                duPO->snp_c += 1;
            }
            if(AnySeqAmbigsI(seqPO)) {
                duPO->amb_c += 1;
            }
    }
    return(TRUE);
}
/**************************************************************************
*   Report stats for given sequence
*/
int HandleDnaOstatI(DNA_UTIL *duPO, char *nameS, SEQ *seqPO, SEQ *fseqPO, FILE *outPF)
{
    int na,nn,lc,snp, slen;
    char *seqPC;

    if(!GetSeqSeqI(seqPO, &seqPC)) {
        return(FALSE);
    }
    slen = GetSeqLenI(seqPO);
    HAND_NFILE(outPF);
    na = CountSeqAmbigsI(seqPC,0,slen);
    nn = CountSeqAmbigDegensI(seqPC,0,slen,4);
    CountStringCaseI(seqPC,slen,&lc,NULL);
    snp = CountSeqSnpSitesI(seqPC,0,slen);
    fprintf(outPF,"%s\t%d\t%d\t%d\t%d\t%d",nameS,slen,na,nn,lc,snp);
    if(duPO->do_ds) {
        GetSeqSeqI(seqPO, &seqPC);
        fprintf(outPF,"\t%s",seqPC);
    }
    if(duPO->do_di) {
        GetSeqSeqI(fseqPO, &seqPC);
        fprintf(outPF,"\t%s",seqPC);
    }
    fprintf(outPF,"\n");
    return(TRUE);
}
/**************************************************************************
*   Report stats for collection of seqs looked at
*/
void HandleDuStats(DNA_UTIL *duPO,int n,int nok,FILE *outPF)
{
    char bufS[DEF_BS];

    HAND_NFILE(outPF);
    fprintf(outPF,"Name     %s\n",duPO->inname);
    fprintf(outPF,"Number   %d\n",n);
    FillSeqFtypeDescString(duPO->iform,bufS);
    fprintf(outPF,"Format   %s\n",bufS);
    if(n<1) {
        return;
    }
    fprintf(outPF,"NumOk    %d  %5.2f%%\n",nok,PERCENT_R(nok,n));
    if(nok<1) {
        return;
    }
    fprintf(outPF,"MinLen   %d\n",duPO->len_lo);
    fprintf(outPF,"MaxLen   %d\n",duPO->len_hi);
    if(duPO->snp_hi>0) {
        fprintf(outPF,"MinSNP   %d\n",duPO->snp_lo);
        fprintf(outPF,"MaxSNP   %d\n",duPO->snp_hi);
        fprintf(outPF,"WithSNPs %d  %5.2f%%\n",duPO->snp_c,
            PERCENT_R(duPO->snp_c,n));
    }
    else {
        fprintf(outPF,"WithSNPs 0\n");
    }
    if(duPO->amb_c>0) {
        fprintf(outPF,"WithAmbs %d  %5.2f%%\n",duPO->amb_c,
            PERCENT_R(duPO->amb_c,n));
    }
    else {
        fprintf(outPF,"WithAmbs 0\n");
    }
}
/**************************************************************************
*   Output header describing stat columns
*/
void HandleDnuFlagHeader(DNA_UTIL *duPO,FILE *outPF)
{
    char okS[DEF_BS];
    if( (!duPO->n_flags) && (!duPO->do_flg) ) {
        return;
    }
    HAND_NFILE(outPF);
    fprintf(outPF,"# Filtering criteria:");
    if(duPO->do_not) {
        fprintf(outPF," (Inverted)");
        sprintf(okS,"Filtered out");
    }
    else {
        sprintf(okS,"Selected");
    }
    fprintf(outPF,"\n");
    if( duPO->min_len > 0) {
        fprintf(outPF,"#  Length: %d to %d %s\n", duPO->min_len, duPO->max_len, okS);
    }
    if( duPO->first > 0) {
        fprintf(outPF,"#  Inputs: %d to %d %s\n", duPO->first, duPO->last, okS);
    }
    if( duPO->do_famb ) {
        fprintf(outPF,"#  Ambigs: %s\n", okS);
    }
    if( duPO->do_fsnp ) {
        fprintf(outPF,"#  SNPs:   %s\n", okS);
    }
    if(!NO_S(duPO->wlisname)) { 
        fprintf(outPF,"#  Listed (%s) %s\n", duPO->wlisname,okS);
    }
}
/**************************************************************************
*   Output header describing stat columns
*/
void HandleDnuOstatHeader(DNA_UTIL *duPO,FILE *outPF)
{
    if(duPO->owhat != DNUO_STAT)
    {
        return;
    }
    HAND_NFILE(outPF);
    fprintf(outPF,"# Reporting sequence stats\n");
    fprintf(outPF,"#   Col1 = Name\n");
    fprintf(outPF,"#   Col2 = Length\n");
    fprintf(outPF,"#   Col3 = Ambiguous base count (including N)\n");
    fprintf(outPF,"#   Col4 = N (masked) base count\n");
    fprintf(outPF,"#   Col5 = Lowercase base count\n");
    fprintf(outPF,"#   Col6 = SNP annoation count; i.e. as \"...[X/Y]...\"\n");
}
/**************************************************************************
*   Output header describing info indices for sequences
*/
void HandleDnuInfoHeader(DNA_UTIL *duPO,FILE *outPF)
{
    if(duPO->owhat != DNUO_INFO)
    {
        return;
    }
    /***
    *   Tell story
    */
    HAND_NFILE(outPF);
    if( (duPO->ifwmin>0) || (duPO->ifwmax>0 ) )
    {
        fprintf(outPF,"# Reporting normalized \"word-count info\" scores\n");
        fprintf(outPF,"#   These are indicative of sequence complexity\n");
        fprintf(outPF,"#   They range from 0 (minimal) to 1 (maximal)\n");
        fprintf(outPF,"#   A 0 indicates all subsequence words are the same\n");
        fprintf(outPF,"#   A 1 indicates all subsequence words differ\n");
        fprintf(outPF,"#   Words (i.e. n-mers) of size %d to %d used here\n",
            duPO->ifwmin,duPO->ifwmax);
        fprintf(outPF,"#\n");
        return;
    }
    else if(duPO->do_insh)
    {
        fprintf(outPF,"# Reporting Shannon diversity and evenness indices\n");
        fprintf(outPF,"#   Words (i.e. n-mers) of size 1 to 3 used here\n");
        fprintf(outPF,"#   Columns 2-4 have diverity for 1- 2- and 3-mers\n");
        fprintf(outPF,"#   Columns 5-7 have evenness for 1- 2- and 3-mers\n");
        fprintf(outPF,"#\n");
        return;
    }
    else if(duPO->do_inbr)
    {
        fprintf(outPF,"# Reporting Brillouin diversity indices\n");
        fprintf(outPF,"#   Words (i.e. n-mers) of size 1 to 3 used here\n");
        fprintf(outPF,"#   Columns 2-4 correspond to 1- 2- and 3-mers\n");
        fprintf(outPF,"#\n");
    }
    else if(duPO->do_inbp)
    {
        fprintf(outPF,"# Reporting normalized Berger-Park diversity indices\n");
        fprintf(outPF,"#   Words (i.e. n-mers) of size 1 to 3 used here\n");
        fprintf(outPF,"#   Columns 2-4 correspond to 1- 2- and 3-mers\n");
    }
}
/**************************************************************************
*   Output info indices for sequences
*/
int HandleDnaInfoOutI(DNA_UTIL *duPO, SEQ *seqPO, SEQ *fseqPO, FILE *outPF)
{
    int slen;
    char *seqPC,nameS[NSIZE];
    DOUB dD,d2D,d3D;

    HAND_NFILE(outPF);
    GetSeqSeqI(seqPO,&seqPC);
    slen = GetSeqLenI(seqPO);
    FillSeqNameStringI(seqPO,nameS,NSIZE-1);
    if( (duPO->ifwmin>0) || (duPO->ifwmax>0 ) ) {
        dD = SeqWordFreqInfoD(seqPC, slen, duPO->ifwmin, duPO->ifwmax);
        fprintf(outPF,"%s\t%5.4f",nameS,dD);
    }
    /***
    *   Shannon div & evenness indices 
    */
    else if(duPO->do_insh) {
        if(!Seq123ShannonInfoI(seqPC, slen, &dD, &d2D, &d3D))
        {
            fprintf(outPF,"%s\tPROBLEM WITH SHANNON INDEX\n",nameS);
            return(FALSE);
        }
        fprintf(outPF,"%s\t%4.3f\t%4.3f\t%4.3f",nameS,dD,d2D,d3D);
        Seq123ShannonEvenInfoI(seqPC, slen, &dD, &d2D, &d3D);
        fprintf(outPF,"\t%4.3f\t%4.3f\t%4.3f",dD,d2D,d3D);
    }
    /***
    *   Various info indices 
    */
    else if(duPO->do_inbr) {
        if(!Seq123BrillouinInfoI(seqPC, slen, &dD, &d2D, &d3D)) {
            fprintf(outPF,"%s\tPROBLEM WITH BRILLOUIN INDEX\n",nameS);
            return(FALSE);
        }
        fprintf(outPF,"%s\t%4.3f\t%4.3f\t%4.3f",nameS,dD,d2D,d3D);
    }
    else if(duPO->do_inbp) {
        if(!Seq123BergerParkerInfoI(seqPC, slen, &dD, &d2D, &d3D)) {
            fprintf(outPF,"%s\tPROBLEM WITH BERGER-PARKER INDEX\n",nameS);
            return(FALSE);
        }
        fprintf(outPF,"%s\t%4.3f\t%4.3f\t%4.3f",nameS,dD,d2D,d3D);
    }
    else {
        ERR("HandleDnaInfoOut","sham with div index flag");
        return(FALSE);
    }
    if(duPO->do_ds) {
        GetSeqSeqI(seqPO,&seqPC);
        fprintf(outPF,"\t%s",seqPC);
    }
    if(duPO->do_di) {
        GetSeqSeqI(fseqPO,&seqPC);
        fprintf(outPF,"\t%s",seqPC);
    }
    fprintf(outPF,"\n");
    return(TRUE);
}
/*************************************************************************
*
*/
int HandleCaSepOutputI(DNA_UTIL *duPO, SEQ *seqPO, FILE *outPF)
{
    int i,n,s,e,t,u, slen;
    static int countsIA[4];
    char baseS[DEF_BS], nameS[NSIZE], *seqPC;

    VALIDATE(duPO,DNA_UTIL_ID);
    HAND_NFILE(outPF);
    FillSeqNameStringI(seqPO,nameS,NSIZE-1);
    GetSeqSeqI(seqPO, &seqPC);
    slen = GetSeqLenI(seqPO);
    if(slen<1) {
        return(0);
    }
    n = countsIA[1] = countsIA[2] = countsIA[3] = 0;
    /***
    *   Through the sequence, each case-class noted
    */
    t = CaSepCharClass(seqPC[0]);
    s = e = 1;
    for(i=1;i<slen;i++) 
    {
        u = CaSepCharClass(seqPC[i]);
        if(u != t) {
            e = i;
            n += UpSepCharCountStory(t, countsIA, baseS);
            /***
            *   Separate seqs or just story
            */
            HandCaSepOneSeqOut(duPO, seqPC, s, e, nameS, baseS, outPF);
            /***
            *   Next char segment initialize...
            */
            t = u;
            s = e + 1;
        }
    }
    n += UpSepCharCountStory(t, countsIA, baseS);
    HandCaSepOneSeqOut(duPO, seqPC, s, i, nameS, baseS, outPF);
    return(n);
}
/*************************************************************************
*   Returns class code for given charater (i.e. "base")
*/
int CaSepCharClass(char c)
{
    int i;

    /***
    *   ACGT = TRUE >--> 1 UPPER case, 2 lower case
    *   IUB = BOGUS >--> 3
    *   Unknown = FALSE >--> 4
    */
    i = GoodDNABaseI(c);
    if(i == TRUE) {
        i = isupper(c) ? 1 : 2;
    }
    else {
        i = IS_BOG(i) ? 3: 4;
    }
    return(i);
}
/*************************************************************************
*   Update counts based on class code, also set sep-fragment story string
*/
int UpSepCharCountStory(int t, int *cPI, char *nS)
{
    char wordS[] = "_ULA";

    if( (t<1) || (t>3) ) {
        return(FALSE);
    }
    cPI[t] += 1;
    sprintf(nS,"%c%d",wordS[t],cPI[t]);
    return(TRUE);
}
/*************************************************************************
*   Report output for single case-sep sub-seq
*/
void HandCaSepOneSeqOut(DNA_UTIL *duPO, char *seqS, int s, int e, char *nameS, char *baseS, 
    FILE *outPF)
{
    SEQ *seqPO;
    char onameS[NSIZE];

    VALIDATE(duPO,DNA_UTIL_ID);
    HAND_NFILE(outPF);
    /***
    *   If sep seqs or (only) stats
    */
    if(duPO->do_csep) {
        if(duPO->do_cstat) {
            sprintf(onameS,"%s_%s__%d_%d",nameS,baseS,s,e);
        }
        else {
            sprintf(onameS,"%s_%s",nameS,baseS);
        }
        seqPO = CreateSeqPO(e-s+1, &seqS[s-1], onameS);
        WriteSeq(seqPO,duPO->oform,outPF); 
        CHECK_SEQ(seqPO);
    }
    else if(duPO->do_cstat) {
        fprintf(outPF,"%s\t%s\t%d\t%d\n",nameS,baseS,s,e);
    }
}
/*************************************************************************
*   
*/
void InvertSeqCase(SEQ *seqPO)
{
    int i,len;
    char *seqPC;

    len = GetSeqLenI(seqPO);
    GetSeqSeqI(seqPO, &seqPC);
    for(i=0; i<len; i++)
    {
        if( ISUPPER(seqPC[i]) ) {
            seqPC[i] = TOLOWER(seqPC[i]);
        }
        else if( ISLOWER(seqPC[i]) ) {
            seqPC[i] = TOUPPER(seqPC[i]);
        }
    }
}
