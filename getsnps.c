/*
* getsnps.c
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#define __MAIN__
#include "prim.h"
#include "dna.h"
#include "getsnps.h"

#define DB_GS   if(DB[107])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(GetSnpsI(argc,argv),NULL) ); }
/**************************************************************************/
void GetSnpsIUse(void)
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Use <infile>  [...options]\n");
    printf("   <infile>   SNP fasta file from dbSNP (i.e. rs.fasta)\n");
    printf("   -iraw      Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq      Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas      Treat input as fasta format\n");
    printf("   -isdb      Treat input as SNPdb-format fasta file\n");
    printf("   -out XXX   Set output to XXX\n");
    printf("   -olis      Output list of qualified records\n");
    printf("   -osum      Output summary (dump)\n");
    printf("   -ofas      Output in fasta format\n");
    printf("   -oraw      Output in \"raw\" format\n");
    printf("   -oexs      Output explicit sequences for each SNP\n");
    printf("   -oexd      Output explicit sequences named with dashs \"-\"\n");
    printf("   -cs        Collapse SNPs to single letter IUPAC codes\n");
    printf("   -win # #   Extract windows # up and # down around SNP sites\n");
    printf("   -awin      Extract all windows, even if can't fully extend\n");
    printf("   -twin      Trim windows of any ambiguous bases\n");
    printf("   -nwin #    Set max ambiguous bases (i.e. N) in windows to #\n");
    printf("   -iub XXX   Limit to SNP IUPAC code(s) XXX (multiple ok)\n");
    printf("   -nal #     Limit to number of alleles #\n");
    printf("   -slen #    Limit to supplied seqs of length # to #\n");
    printf("   -class #   Limit to class #\n");
/*
    printf("   -msnp #    Set max \"SNP\" size to # (default = %d)\n",DEF_MSNP);
    printf("   -vtype XXX Limit to VarType descripion XXX\n");
    printf("   -cstat XXX Limit to CurationStatus descripion XXX\n");
    printf("   -gr XXX    Limit to GeneRegion description XXX\n");
*/
    printf("   -not       Invert qualification tests\n");
    printf("   -quiet     No report of progress while processing\n");
    printf("   -noh       No header info reported\n");
    printf("   -all       Process all records; don't bail on problems\n");
}
/**************************************************************************/
int GetSnpsI(int argc,char **argv)
{
    int uwin,dwin,min,max,not,verb,ok,n,nok,class,isdb;
    int ofas,oraw,header,iraw,iseq,ifas,all,c,osum,olis;
    char inS[NSIZE],outS[NSIZE],iubS[NSIZE];
    char vtypeS[DEF_BS], cstatS[DEF_BS], generS[DEF_BS];
    FILE *inPF, *outPF;
    SNPREC snO, *snPO;
    SEQ *seqPO;

    snPO = &snO;
    InitSnprec(snPO,TRUE);
    class = BOGUS;
    uwin = dwin = -1;
    min = 0; max = TOO_BIG;
    isdb = ofas = oraw = iraw = iseq = ifas = all = osum = olis = FALSE;
    verb = header = TRUE;
    INIT_S(outS); INIT_S(iubS); 
    INIT_S(vtypeS); INIT_S(cstatS); INIT_S(generS); 
    if(!ParseArgsI(argc,argv,
        "S -out S -win I2 -iub S -class I -slen I2 -not B -quiet B\
        -isdb B -vtype S -cstat S -gr S -cs B -oexs B\
        -ofas B -oraw B -noh B -iraw B -awin B -nwin I -ifas B\
        -all B -nal I -twin B -oexd B -osum B -olis B -iseq B",
        inS, outS, &uwin,&dwin, iubS, &class, &min,&max, &not, &verb,
        &isdb, vtypeS, cstatS, generS, &snPO->do_cs, &snPO->do_oexs,
        &ofas, &oraw, &header, &iraw, &snPO->do_awin, &snPO->do_nwin, &ifas,
        &all, &snPO->nal,
        &snPO->do_twin, &snPO->do_oexd, &osum, &olis, &iseq,
        (int *)NULL))
    {
        GetSnpsIUse();
        return(FALSE);
    }
    /***
    *   Set input format to explitly chosen, else guess (if not stdin)
    */
    if(isdb) {
        snPO->in_form = INF_SNPDB;
    }
    else {
        snPO->in_form =  FigureSeqFileTypeI(iraw,iseq,ifas,inS,TRUE);
    }
    if(!snPO->in_form) {
        GetSnpsIUse();
        printf("\n");
        printf(" >>> An input format must be specified! <<<\n");
        printf("\n");
        return(FALSE);
    }
    /***    
    *   Output
    */
    if(ofas) {
        snPO->out_form = SEQFM_FASTA;
    }
    else if(oraw) {
        snPO->out_form = SEQFM_RAW;
    }
    Upperize(iubS);
    /***
    *   Input and output files?
    */
    if(!(inPF=OpenUFilePF(inS,"r",NULL))) {
        return(FALSE);
    }
    outPF = NULL;
    if(!NO_S(outS)) {
        if(!(outPF=OpenUFilePF(outS,"w",NULL))) {
            FILECLOSE(inPF);
            return(FALSE);
        }
    }
    HAND_NFILE(outPF);
    if(header) {
        OutputHeader(snPO,inS,iubS,class,min,max,not,outS,outPF);
        if( (uwin>0) || (dwin>0) ) {
            OutputHeader(snPO,inS,iubS,class,min,max,not,outS,outPF);
            fprintf(outPF,"# SNP windows %d up / %d down\n",uwin,dwin);
            if(snPO->do_twin) {
                fprintf(outPF,"# SNP windows trimmed of ambiguous bases\n");
            }
        }
        fprintf(outPF,"#\n");
    }
    /***
    *   Party through file; getting info from fasta header lines
    */
    seqPO = CreateSeqPO(0,NULL,NULL);
    n = nok = 0;
    while(TRUE)
    {
        /***
        *   Load record in various formats
        */
        ok = FALSE;
        switch(snPO->in_form)
        {
            case INF_SNPDB:
                ok = LoadSNPdbFastaRecI(inPF, snPO);
                break;
            case SEQFM_RAW:
            case SEQFM_FASTA:
                ok = ParseSeqI(inPF, snPO->in_form, n+1, FALSE, TRUE, seqPO);
                if(ok > 0) {
                    ok = ExtractSNPInfoFromSeqI(seqPO,snPO);
                }
                else {
                    ok = FALSE;     /* Parse can return TRUE, FALSE or BOGUS */
                }
                break;
        }
        if(!ok) {
            if(!all) {
                break;
            }
            if( (c=fgetc(inPF)) == EOF) {
                break;
            }
            ungetc(c,inPF);
            continue;
        }
        /***
        *   Process SNP info
        */
        ParseSnpAllelesI(snPO);
        if(snPO->do_cs) {
            strcpy(snPO->snp,snPO->iubc);
        }
        /***
        *   Qualification tests 
        */
        ok = OkSnpRecI(snPO,iubS,class,min,max,vtypeS,cstatS,generS);
        if(not) {
            ok = !ok;
        }
        n++;
        if(ok) {
            nok++;
        }
        if( verb && ((n%10000) == 0) ) {
            printf("# %8d; %8d ok\n",n,nok);
            fflush(stdout);
        }
        if(!ok) {
            continue;
        }
        /***
        *   Output what
        */
        if(olis) {
            HAND_NFILE(outPF);
            fprintf(outPF,"%s\n",snPO->name);
        }
        else if(osum) {
            DumpSnprec(snPO, NULL);
        }
        else if(uwin > -1) {
            ok = ExtractSnpWinI(snPO,uwin,dwin,outPF);
            if(verb && (!ok) ) {
                printf("# FAILED to extract window %s (bases: up=%d down=%d)\n",
                    snPO->name, INT(strlen(snPO->sequp)), INT(strlen(snPO->seqdn)) );
            }
        }
        else if(ofas || oraw) {
            ExtractSnpWinI(snPO,-1,-1,outPF);
        }
    }
    if(header)
    {
        printf("# Total records: %d\n",n);
        printf("#    Ok records: %d\n",nok);
    }
    /***
    *   Clean up; all done
    */
    CHECK_SEQ(seqPO);
    CHECK_FILE(inPF);
    CHECK_NFILE(outPF,outS);
    return(1);
}
/**********************************************************************
*   Init structure vars, run-time or everything
*/
void InitSnprec(SNPREC *snPO,int full)
{
    INIT_S(snPO->head);
    INIT_S(snPO->name);
    INIT_S(snPO->seq);
    INIT_S(snPO->sequp);
    INIT_S(snPO->seqdn);
    INIT_S(snPO->snp);
    INIT_S(snPO->iubc);
    INIT_S(snPO->vtype);
    INIT_S(snPO->cstat);
    INIT_S(snPO->freq);
    INIT_S(snPO->genen);
    INIT_S(snPO->gener);
    snPO->taxid = snPO->class = snPO->apos = BOGUS;
    snPO->slen = snPO->snplen = 0;
    snPO->uplen = snPO->dnlen = 0;
    snPO->nales = 0;
    if(full)
    {
        /***
        *   Various run-time global options
        */
        snPO->do_cs = FALSE;
        snPO->do_oexs = FALSE;
        snPO->do_awin = FALSE;
        snPO->do_twin = FALSE;
        snPO->do_nwin = -1;
        snPO->in_form = 0;
        snPO->out_form = SEQFM_RAW;
        snPO->nal = BOGUS;
        snPO->do_oexd = FALSE;
    }
}
/************************************************************************
*
*/
void DumpSnprec(SNPREC *snPO, FILE *outPF)
{
    int i;

    HAND_NFILE(outPF);
    fprintf(outPF,"#####################################################\n");
    fprintf(outPF,"HEAD:   |%s|\n",snPO->head);
    fprintf(outPF,"NAME:   |%s|\n",snPO->name);
    fprintf(outPF,"TAXID:  %d\n",snPO->taxid);
    fprintf(outPF,"CLASS:  %d\n",snPO->class);
    fprintf(outPF,"VTYPE:  |%s|\n",snPO->vtype);
    fprintf(outPF,"CSTAT:  |%s|\n",snPO->cstat);
    fprintf(outPF,"GENEN:  |%s|\n",snPO->genen);
    fprintf(outPF,"GENER:  |%s|\n",snPO->gener);
    fprintf(outPF,"ALLELE: |%s| %d\n",snPO->snp, snPO->snplen);
    for(i=0; i<snPO->nales; i++)
    {
        fprintf(outPF,"ALLELE_%d:|%s|\n",i+1,snPO->ales[i]);
    }
    fprintf(outPF,"IUBC:   |%s|\n",snPO->iubc);
    fprintf(outPF,"FREQS:  |%s|\n",snPO->freq);
    fprintf(outPF,"APOS:   %d\n",snPO->apos);
    fprintf(outPF,"SLEN:   %d\n",snPO->slen);
    fprintf(outPF,"SEQUP:  |%s|\n",snPO->sequp);
    fprintf(outPF,"SEQDN:  |%s|\n",snPO->seqdn);
    fprintf(outPF,"SEQ:    |%s|\n",snPO->seq);
}
/**********************************************************************
*   Load snp record in "SNPdb" fasta format
*/
int LoadSNPdbFastaRecI(FILE *inPF, SNPREC *snPO)
{
    int c, i, u, d, b;
    char bufS[BBUFF_SIZE];

    DB_GS DB_PrI(">> eoadSNPdbFastaRecI\n");
    InitSnprec(snPO,FALSE);
    /***
    *   First get header line
    */
    while(fgets(bufS,BLINEGRAB,inPF))
    {
        if(COM_LINE(bufS))
            continue;
        if(bufS[0] == '>')
        {
            ReplaceChars('\n',bufS,'\0',snPO->head);
            break;
        }
    }
    if(NO_S(snPO->head)) {
        DB_GS DB_PrI("<< LoadSNPdbFastaRecI no head FALSE\n");
        return(FALSE);
    }
    if(! ParseSNPdbHeadFields(snPO, snPO->head) ) {
        PROBLINE;
        printf("Failed to parse SNP header line\n|%s|\n",snPO->head);
        return(FALSE);
    }
    /***
    *   Now get up, down and snp sequence
    *   Not sure of lines / spaces so char at a time...
    */
    i = u = d = b = 0;
    while((c = fgetc(inPF)) != EOF)
    {
        /*
        *   Start of next fasta ...or comment 
        */
        if( (c=='>') || (c=='#') ){
            ungetc(c,inPF);
            break;
        }
        /***
        *   More than one non-print (e.g. new lines) = done
        */
        if( (!isgraph(c)) && (c!=' ') ) {
            b++;
            if(b>1) {
                break;
            }
            continue;
        }
        b = 0;
        if(c==' ') {
            continue;
        }
        /***
        *   SNP position; Position is 1-based 
        */
        i++;
        if(i < snPO->apos) {
            snPO->sequp[u++] = c;
        }
        else if (i > snPO->apos) {
            snPO->seqdn[d++] = c;
        }
    }
    snPO->sequp[u] = '\0';
    snPO->seqdn[d] = '\0';
    snPO->uplen = u;
    snPO->dnlen = d;
    if( NO_S(snPO->sequp) || NO_S(snPO->sequp) ) {
        PROBLINE;
        printf("Failed to get up/down sequence; Header:\n|%s|\n",snPO->head);
        return(FALSE);
    }
    DB_GS DB_PrI("+ up |%s|\n+ sn |%s|\n+ dn |%s|\n", snPO->sequp, snPO->snp, snPO->seqdn);
    sprintf(snPO->seq,"%s[%s]%s",snPO->sequp,snPO->snp,snPO->seqdn);
    /***
    *   All done; set values into passed datastruc 
    */
    DB_GS DB_PrI("<< LoadSNPdbFastaRecI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Parse fields from dbSNP header line string
*/
int ParseSNPdbHeadFields(SNPREC *snPO, char *headS)
{
    int n,m;
    char *cPC, bufS[BBUFF_SIZE], pairS[NSIZE], tokS[NSIZE], valS[NSIZE];

    /***
    *   Assume looks like below. First break into tokens then process each...
    *   >gnl|dbSNP|rs3087742 rs=3087742|pos=512|len=524|taxid=9606|mol="genomic"|class=1|alleles="C/T"|build=132
    */
    ReplaceSomeChars('|',headS,' ',bufS,BBUFF_SIZE-1);
    cPC = bufS;
    PASS_BLANK(cPC);
    n = m = 0;
    while(ISLINE(*cPC)) {
        INIT_S(pairS);
        INIT_S(tokS);
        INIT_S(valS);
        sscanf(cPC,"%s",pairS);
        if(n==2) {
            strcpy(snPO->name,pairS);
            m++;
        }
        else if(strstr(pairS,"=")) {
            ReplaceChars('=',pairS,' ',pairS);
            sscanf(pairS,"%s %s",tokS,valS);
            if( NO_S(tokS) || NO_S(valS) ) {
                printf("Problem parsing |%s| from line |%s|\n",pairS,headS);
                return(FALSE);
            }
            Upperize(tokS);
            if(EQSTRING(tokS,"LEN",3)) {
                sscanf(valS,"%d",&snPO->slen);
                m++;
            }
            if(EQSTRING(tokS,"POS",3)) {
                sscanf(valS,"%d",&snPO->apos);
                m++;
            }
            if(EQSTRING(tokS,"TAX",3)) {
                sscanf(valS,"%d",&snPO->taxid);
                m++;
            }
            if(EQSTRING(tokS,"CLASS",5)) {
                sscanf(valS,"%d",&snPO->class);
                m++;
            }
            if(EQSTRING(tokS,"ALLELES",7)) {
                /***
                *   SNP may be wrapped in sing / double quotes "A/G" 'a/g'
                */
                RemoveChars("\"'", valS, valS);
                sscanf(valS,"%s",snPO->snp);
                snPO->snplen = strlen(snPO->snp);
                m++;
            }
        }
        n++;
        NEXT_WORD(cPC);
    }
    return(m);
}
/*************************************************************************
*   Qualify current SNP 
*/
int OkSnpRecI(SNPREC *snPO, char *iubS, int class, int min, int max,
    char *vtypeS, char *cstatS, char *generS)
{
    int ok;

    if( (class > 0) && (class != snPO->class) ) {
        return(FALSE);
    }
    if(min > 0) {
        if( (snPO->slen<min) || (snPO->slen>max) ) {
            return(FALSE);
        }
    }
    /***
    *   Base degeneracy
    */
    if(!IS_BOG(snPO->nal)) {
        if(snPO->nales  != snPO->nal) {
            return(FALSE);
        }
    }
    /***
    *   Any of (possibly listed) IUPAC codes 
    */
    if(!NO_S(iubS)) {
        ok = FALSE;
        while(isgraph(INT(*iubS)))
        {
            if(snPO->iubc[0] == *iubS) {
                ok++;
                break;
            }
            iubS++;
        }
        if(!ok) {
            return(FALSE);
        }
    }
    if(!NO_S(vtypeS)) {
        if(!strstr(snPO->vtype,vtypeS)) {
            return(FALSE);
        }
    }
    if(!NO_S(cstatS)) {
        if(!strstr(snPO->cstat,cstatS)) {
            return(FALSE);
        }
    }
    if(!NO_S(generS)) {
        if(!strstr(snPO->gener,generS)) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Write out SNP window from uwin to dwin 
*/
int ExtractSnpWinI(SNPREC *snPO, int uwin, int dwin, FILE *outPF)
{
    char nameS[NSIZE], seqS[SEQMAX],snpS[SNPMAX],upS[SEQMAX],dnS[SEQMAX], sepS[10];
    int e,s,n;

    DB_GS DB_PrI(">> ExtractSnpWinI uwin=%d dwin=%d\n",uwin,dwin);
    HAND_NFILE(outPF);
    /***
    *   Figure / check extraction dims
    */
    if(!SetWindowStartEndI(snPO,uwin,dwin,&s,&e)) {
        DB_GS DB_PrI("<< ExtractSnpWinI window start/end FALSE\n");
        return(FALSE);
    }
    /***
    *   Get sequence parts into local strings
    */
    strcpy(snpS,snPO->snp);
    strcpy(upS,&snPO->sequp[s]);
    strncpy(dnS,snPO->seqdn,e);
    dnS[e] = '\0';
    /***
    *   Check for ambigs in the window
    */
    if(snPO->do_nwin > -1) {
        n = CountSeqAmbigsI(upS,-1,-1);
        n += CountSeqAmbigsI(dnS,-1,-1);
        if(n > snPO->do_nwin) {
            DB_GS DB_PrI("<< ExtractSnpWinI ambigs FALSE\n");
            return(FALSE);
        }
    }
    /***
    *   Whole thing only
    */
    if(! snPO->do_oexs) {
        sprintf(seqS,"%s[%s]%s",upS,snpS,dnS);
        WriteOutSequence(snPO->name, seqS, snPO->out_form, outPF);
        DB_GS DB_PrI("<< ExtractSnpWinI TRUE\n");
        return(TRUE);
    }
    /***
    *   Separator character for names
    */
    if(snPO->do_oexd) {
        sprintf(sepS,"-");
    }
    else {
        sprintf(sepS,"_");
    }
    /***
    *   Each allele
    */
    n = 0;
    while(n<snPO->nales) {
        RemoveChars("-*",snPO->ales[n],snpS);
        if(NO_S(snpS)) {
            sprintf(nameS,"%s%s%d%s%s",snPO->name,sepS,n+1,sepS,"del");
        }
        else {
            sprintf(nameS,"%s%s%d%s%s",snPO->name,sepS,n+1,sepS,snpS);
        } 
        sprintf(seqS,"%s%s%s",upS,snpS,dnS);
        WriteOutSequence(nameS, seqS, snPO->out_form, outPF);
        n++;
    }
    DB_GS DB_PrI("<< ExtractSnpWinI TRUE\n");
    return(TRUE);
}
/***************************************************************************
*   Handle limiting / adjusting of SNP window start & end output indices
*   For negative window size, whole thing 
*/
int SetWindowStartEndI(SNPREC *snPO,int uwin,int dwin,int *sPI,int *ePI)
{
    int i,s,e,ulen,dlen;

    DB_GS DB_PrI(">> SetWindowStartEndI uwin=%d dwin=%d\n",uwin,dwin);
    ulen = snPO->uplen;
    dlen = snPO->dnlen;
    /***
    *   Window start / end
    */  
    if(uwin < 0) {
        s = 0;
    }   
    else {
        s = ulen - uwin;
    }   
    if(dwin < 0) {
        e = dlen;
    }   
    else {
        e = dwin;
    }   
    DB_GS DB_PrI("+ ulen=%d dlen=%d s=%d e=%d\n",ulen,dlen,s,e);
    if( (s<0) || (e>dlen) ) {
        if(!snPO->do_awin) {
            DB_GS DB_PrI("<< SetWindowStartEndI can't fit FALSE\n");
            return(FALSE);
        }
        LIMIT_NUM(s,0,ulen);
        LIMIT_NUM(e,0,dlen);
    }
    DB_GS DB_PrI("+ limited s=%d e=%d\n",s,e);
    /***
    *   Trimming any ambigs?
    */
    if(snPO->do_twin)
    {
        for(i=(ulen-1);i>=s;i--) 
        {
            if( BaseDegeneracyI(snPO->sequp[i]) != 1) {
                s = i + 1;
                break;
            }
        }
        for(i=0;i<e;i++) 
        {
            if( BaseDegeneracyI(snPO->seqdn[i]) != 1) {
                e = i;
                break;
            }
        }
    }
    *sPI = s;
    *ePI = e;
    DB_GS DB_PrI("<< SetWindowStartEndI s=%d e=%d TRUE\n",s,e);
    return(TRUE);
}
/***************************************************************************
*
*/
void WriteOutSequence(char *nameS, char *seqS, int oform, FILE *outPF)
{
    HAND_NFILE(outPF);
    if(oform == SEQFM_FASTA) {
        fprintf(outPF,">%s\n",nameS);
    }
    else {
        fprintf(outPF,"%-12s\t",nameS);
    }
    fprintf(outPF,"%s\n",seqS);
}
/*************************************************************************/
void OutputHeader(SNPREC *snPO, char *inS, char *snpS, int class, int min, 
    int max, int not, char *outS,FILE *fPF)
{
    int line;
    HAND_NFILE(fPF);

    fprintf(fPF,"# %s\n",outS);
    fprintf(fPF,"# Output by %s, %s %s %s\n",VERSION_S,BD_S,__DATE__,__TIME__);
    fprintf(fPF,"# %s\n",RTK_S);
    fprintf(fPF,"# Input file: %s\n",inS);
    TimeStamp("# ",fPF);
    line = 0;
    if(!NO_S(snpS))
    {
        if(!line)
            fprintf(fPF,"#\n");
        line++;
        fprintf(fPF,"# SNP type restriction: %s\n",snpS);
    }
    if(snPO->nal > 0)
    {
        if(!line)
            fprintf(fPF,"#\n");
        line++;
        fprintf(fPF,"# SNP restricted to %d alleles\n",snPO->nal);
    }
    if(class > 0)
    {
        if(!line)
            fprintf(fPF,"#\n");
        line++;
        fprintf(fPF,"# SNP class restriction: %d\n",class);
    }
    if(min > 0)
    {
        if(!line)
            fprintf(fPF,"#\n");
        line++;
        fprintf(fPF,"# SNP SeqLen restriction: %d %d\n",min,max);
    }
    if( (not) && (line>0) )
    {
        fprintf(fPF,"# All RESTRICTIONS INVERTED\n");
    }
    fprintf(fPF,"#\n");
}
/***************************************************************************
*
*/
int ParseSnpAllelesI(SNPREC *snPO)
{
    char *cPC,snpS[SNPMAX],alS[SNPMAX];
    int n;

    ReplaceSomeChars('/',snPO->snp,' ',snpS,SNPMAX-1);
    cPC = snpS;
    PASS_BLANK(cPC);
    n = 0;
    while( (ISLINE(*cPC)) && (n<N_SNP) ) {
        INIT_S(alS);
        sscanf(cPC,"%s",alS);
        sprintf(snPO->ales[n],"%s",alS);
        n++;
        NEXT_WORD(cPC);
    }
    snPO->nales = n;
/*
    CollapseSnpString(snPO->snp,snPO->iubc);
*/
    CollapseSnpStringI(snPO->snp,snPO->iubc);
    Upperize(snPO->iubc);
    return(n);
}
/***************************************************************************
*
*/
int ExtractSNPInfoFromSeqI(SEQ *seqPO, SNPREC *snPO)
{
    int i,u,d,s,b;

    InitSnprec(snPO,FALSE);
    FillSeqNameStringI(seqPO, snPO->name, -1);
    FillSeqSeqStringI(seqPO, snPO->seq, SEQMAX-1);
    i = u = d = s = b = 0;
    while(snPO->seq[i] != '\0')
    {
        if( (snPO->seq[i] == '[') || ( snPO->seq[i] == ']') ) {
            b++;
        }
        else
        {
            switch(b)
            {
                case 0:     snPO->sequp[u++] = snPO->seq[i];    break;
                case 1:     snPO->snp[s++] = snPO->seq[i];      break;
                case 2:     snPO->seqdn[d++] = snPO->seq[i];    break;
                default:
                    PROBLINE;
                    printf("Sequence has bad [x/y] SNP:\n|%s|\n",snPO->seq);
                    return(FALSE);
            }
        }
        i++;
    }
    snPO->sequp[u] = '\0';
    snPO->snp[s] = '\0';
    snPO->seqdn[d] = '\0';
    snPO->uplen = u;
    snPO->dnlen = d;
    return(TRUE);
}
