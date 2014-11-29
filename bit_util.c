/*
* bit_util.c
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
#include "bitpool.h"
#include "table.h"
#include "bit_util.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(Bit_utilI(argc,argv),NULL) ); }
/**************************************************************************/
void Bit_utilUse(void)
{
	VersionSplash(NULL,VERSION_S,"#  ",TRUE);
	printf("Usage: <infile> ['-' for stdin] [...options]\n");
	printf("   <infile>   Bit record file to read in\n");
    printf("   -ibit      Input bit format: <name> <01011...> per line\n");
    printf("   -itab      Input table format\n");
    printf("   -stp       Strict parsing; Default is to skip input errors\n");
	printf("   -out XXX   Set output to XXX\n");
    printf("   -obit      Output bit format\n");
    printf("   -osbit     Output spaced bit format\n");
    printf("   -otab      Output table\n");
    printf("   -dump      Dump bitstring data\n");
	printf("   -rran # #  Record range restricted # to # (1-base count)\n");
	printf("   -ostat     Output bitstring stats\n");
	printf("   -stat      Report overall stats for input file\n");
	printf("   -fmat      Full matrix\n");
	printf("   -rap       Report all pairs, even if redundant\n");
	printf("   -sbit XXX  Second set of bit records from file XXX\n");
    printf("   -band      Bitwise AND between <in> and <sbin> records\n");
    printf("   -bor       Bitwise OR  between <in> and <sbin> records\n");
    printf("   -bxor      Bitwise XOR between <in> and <sbin> records\n");
    printf("   -b3        All 3 bitwise operators between <in> and <in'>\n");
    printf("   -bnot      Bitwise NOT for records\n");
	printf("   -twn #     Tweak number # of bits\n");
	printf("   -twf #     Tweak faction # of bits\n");
	printf("   -seed #    Set random seed for tweaks\n");
    printf("\n");
}
/**************************************************************************
*	Main program
*/
int Bit_utilI(int argc, char **argv)
{
    int itab, ibit, otab, obit, osbit;
	BIT_UTIL *buPO;

	buPO = CreateBit_utilPO();
    itab = ibit = FALSE;
    otab = obit = osbit = FALSE;
	if(!ParseArgsI(argc,argv,
        "S -out S -rran I2 -dump B\
        -stat B -ostat B -fmat B\
        -ibit B -itab B -obit B -osbit B -otab B\
        -twn I -twf D -sbit S -band B -bor B -bxor B -bnot B\
        -seed I -b3 B -stp B -rap B",
		buPO->inname, buPO->outname, 
        &buPO->firstr,&buPO->lastr, &buPO->do_dump,
        &buPO->do_stat, &buPO->do_ostat, &buPO->do_fmat,
        &ibit, &itab, &obit, &osbit, &otab,
        &buPO->tw_num, &buPO->tw_frac, &buPO->sbname,
        &buPO->do_band, &buPO->do_bor, &buPO->do_bxor, &buPO->do_bnot,
        &buPO->seed, &buPO->do_b3, &buPO->istrict, &buPO->do_rap,
        (int *)NULL))
    {
        Bit_utilUse();
		CHECK_BIT_UTIL(buPO);
		return(FALSE);
	}
	/***
    *   Set input / output format 
	*/
    buPO->iform = FigureBitFileTypeI(ibit,ibit,FALSE,itab,buPO->inname,TRUE);
	if(!buPO->iform) {
        printf("Problem with input seq(s)\n");
		CHECK_BIT_UTIL(buPO);
		return(FALSE);
	}
	buPO->oform = buPO->iform;
    if( obit || osbit || otab ) {
        buPO->oform = FigureBitFileTypeI(obit,osbit,FALSE,otab,buPO->outname,TRUE);
    }
	/***
	*	Set up / load
	*/
	if(!SetupBtuI(buPO)) {
		CHECK_BIT_UTIL(buPO);
		return(FALSE);
	}
    if(!CheckBtuOptionsI(buPO)) {
		CHECK_BIT_UTIL(buPO);
		return(FALSE);
    }
	HandleBtuHeader(buPO,buPO->out);
	/***
	*   Do whatever and output
	*/
    HandleBtuModsI(buPO);
    HandleBtuFiltersI(buPO);
    HandleBtuOutput(buPO,buPO->out);
	/***
	*	All done
	*/
	CHECK_BIT_UTIL(buPO);
	return(TRUE);
}
/*****************************************************************************
*	Create data struct
*/
BIT_UTIL *CreateBit_utilPO()
{
	BIT_UTIL *buPO;

	if(! (buPO = (BIT_UTIL *)ALLOC(1,sizeof(BIT_UTIL)) ) )
	{
		printf("# Failed to allocate working object\n");
		return(NULL);
	}
	buPO->ID = BIT_UTIL_ID;
	InitBit_util(buPO);
	return(buPO);
}
/*****************************************************************************
*	Free datastructure and substructs
*/
int DestroyBit_utilI(BIT_UTIL *buPO)
{
	VALIDATE(buPO,BIT_UTIL_ID);
	CHECK_NFILE(buPO->out,buPO->outname);
	CHECK_BITPOOL(buPO->bits);
	CHECK_FREE(buPO->mmask);
	CHECK_FREE(buPO->bmask);
	CHECK_BITPOOL(buPO->sbits);
	CHECK_BITPOOL(buPO->pretwk);
	FREE(buPO);
	return(TRUE);
}
/*****************************************************************************
*	Set null / default values
*/
void InitBit_util(BIT_UTIL *buPO)
{
	VALIDATE(buPO,BIT_UTIL_ID);

	INIT_S(buPO->inname);
	buPO->iform = BOGUS;
	INIT_S(buPO->outname);
	buPO->out = NULL;
	buPO->oform = BOGUS;
    buPO->outbits = buPO->outcomp = FALSE;
    buPO->bits = buPO->sbits = NULL;
    buPO->bsize = 0;
    buPO->firstr = buPO->lastr = BOGUS;
	INIT_S(buPO->sbname);
    buPO->do_stat = buPO->do_ostat = FALSE;
    buPO->do_fmat = buPO->do_rap = FALSE;
    buPO->do_band = buPO->do_bor = buPO->do_bxor = FALSE;
    buPO->do_b3 = buPO->do_bnot = FALSE;
    buPO->tw_num = 0;
    buPO->tw_frac = 0.0;
    buPO->seed = BAD_I;
    return;
}
/**************************************************************************/
int CheckBtuOptionsI(BIT_UTIL *buPO)
{
    if( buPO->do_band || buPO->do_bor || buPO->do_bxor || buPO->do_b3 ||
        buPO->sbits ) {
        buPO->outcomp = TRUE;
    }
    if(buPO->do_stat || buPO->do_ostat || buPO->do_fmat || buPO->outcomp) {
        buPO->outbits = FALSE;
    }
    else {
        buPO->outbits = TRUE;
    }
	return(TRUE);
}
/**************************************************************************/
int SetupBtuI(BIT_UTIL *buPO)
{
    if(!OpenBtuFilesI(buPO)) {
        return(FALSE);
    }
    if(!LoadBtuBitpoolsI(buPO)) {
        PROBLINE;
        printf("Failed to load bit pools\n");
        return(FALSE);
    }
    if(!SetUpBtuAuxDataI(buPO)) {
        PROBLINE;
        printf("Failed to set up aux data\n");
        return(FALSE);
    }
    if(!SetupTweaksI(buPO)) {
        PROBLINE;
        printf("Failed to set up for random tweaking\n");
        return(FALSE);
    }
    return(TRUE);
}
/**************************************************************************/
int SetupTweaksI(BIT_UTIL *buPO)
{
    int bsize,num;

    GetBitpoolDimsI(buPO->bits, &bsize, &num);
    if(buPO->tw_frac > 0.0) {
        buPO->tw_num = INT(buPO->tw_frac * DNUM(bsize));
    }
    LIMIT_NUM(buPO->tw_num, 0, bsize);
    if(buPO->tw_num) {
        Srand(buPO->seed);
        buPO->bmask = (char *)ALLOC(bsize,sizeof(char));
        buPO->pretwk = CreateBitpoolPO(bsize,num);
        if((!buPO->bmask) || (!buPO->pretwk)) {
            printf("Failed to allocate for tweaks num=%d size=%d\n",num,bsize);
            return(FALSE);
        }
    }
    return(TRUE);
}
/**************************************************************************/
int OpenBtuFilesI(BIT_UTIL *buPO)
{
	if(!NO_S(buPO->outname)) {
		if(!(buPO->out=OpenUFilePF(buPO->outname,"w",NULL))) {
			return(FALSE);
		}
	}
	HAND_NFILE(buPO->out);
	return(TRUE);
}
/**************************************************************************/
int LoadBtuBitpoolsI(BIT_UTIL *buPO)
{
    int ok,fsize,ssize;

    ok = LoadOneNamedBitpoolI(buPO->inname, buPO->iform, buPO->istrict, &buPO->bits);
    if( ok && (!NO_S(buPO->sbname)) ) {
        ok = LoadOneNamedBitpoolI(buPO->sbname, buPO->iform, buPO->istrict, &buPO->sbits);
        if(ok) {
            GetBitpoolDimsI(buPO->bits, &fsize, NULL);
            GetBitpoolDimsI(buPO->sbits, &ssize, NULL);
            if(fsize != ssize) {
                PROBLINE;
                printf("Bit size dimensions differ:\n");
                printf("%s = %d,   %s = %d\n",buPO->inname,fsize,buPO->sbname,ssize);
                ok = FALSE;
            }
        }
    }
    return(ok);
}
/**************************************************************************/
int LoadOneNamedBitpoolI(char *nameS, int iform, int error, BITPOOL **bpPPO)
{
    int ok;
    BITPOOL *bpPO;
    TABLE *tabPO;

    ok = FALSE;
    if( (iform == BPF_BITS) || (iform == BPF_SBITS) || (iform == BPF_HEX) ) {
        ok = GetBitpoolI(nameS, iform, error, &bpPO);
    }
    else if(iform == BPF_TAB) {
        tabPO = NULL;
        ok = GetTableI(nameS,TRUE,TRUE,TRUE,FALSE,&tabPO);
        if(ok) {
            ok = BitpoolFromTableI(tabPO, &bpPO);
        }
        CHECK_TABLE(tabPO);
    }
    if(ok) {
        *bpPPO = bpPO;
    }
    return(ok);
}
/**************************************************************************/
int SetUpBtuAuxDataI(BIT_UTIL *buPO)
{
    GetBitpoolDimsI(buPO->bits, &buPO->bsize, &buPO->num);
    if( (buPO->bsize * buPO->num) < 1 ){
        PROBLINE;
        printf("Bit pool has zero bits (%d) or members (%d)\n",buPO->bsize,buPO->num); 
        return(FALSE);
    }
    if( !(buPO->mmask = (char*)ALLOC(buPO->num,sizeof(char))) ){
        return(FALSE);
    }
    InitArrayI(buPO->mmask,IS_CHAR,0,buPO->num,TRUE);
    return(TRUE);
}
/**************************************************************************/
int HandleBtuModsI(BIT_UTIL *buPO) 
{
    if(buPO->do_bnot) {
        ModThisBitpoolI(buPO->bits,-1,-1,BIT_NOT);
    }
    if(buPO->tw_num) {
        HandleBtuTweaksI(buPO);
    }
    return(TRUE);
}
/**************************************************************************/
int HandleBtuTweaksI(BIT_UTIL *buPO)
{
    int m,b,v;
    DOUB fracD;

    BOG_CHECK(!buPO->bmask);
    fracD = (buPO->tw_frac > 0.0) ? buPO->tw_frac : RNUM(buPO->tw_num)/RNUM(buPO->bsize);
    for(m=0;m<buPO->num;m++) {
        CopyThisBitpoolMemberI(buPO->bits,m, buPO->pretwk,m);
        MaskRandSubsetI(buPO->bmask,buPO->bsize,fracD);
        for(b=0;b<buPO->bsize;b++) {
            if(buPO->bmask[b]) {
                GetThisBitpoolBitI(buPO->bits,m,b,&v);
                v = (v) ? FALSE : TRUE;
                SetThisBitpoolBitI(buPO->bits,m,b,v);
            }
        }
    }
    return(TRUE);
}
/**************************************************************************/
int HandleBtuFiltersI(BIT_UTIL *buPO) 
{
    int i,ok;

    for(i=0;i<buPO->num;i++) {
        ok = IsBitpoolRecOkI(buPO, buPO->bits, i);
        buPO->mmask[i] = ok;
    }
    return(TRUE);
}
/**************************************************************************/
int IsBitpoolRecOkI(BIT_UTIL *buPO, BITPOOL *bpPO, int r)
{
	int ok;

    ok = TRUE;
    /***
    *   Record in range?
    */
    if(buPO->firstr > 0) {
	    if( (r < (buPO->firstr -1)) || (r >= buPO->lastr) ) {
            ok = FALSE;
        }
    }
	return(ok);
}
/*************************************************************************/
void HandleBtuHeader(BIT_UTIL *buPO, FILE *outPF)
{
    int bsize,n;
    char bufS[DEF_BS];

    HAND_NFILE(outPF);
    GetBitpoolDimsI(buPO->bits,&bsize,&n);
    fprintf(outPF,"# Loaded %d records of %d bits each from %s\n",n,bsize,buPO->inname);
    if(buPO->sbits) {
        GetBitpoolDimsI(buPO->sbits,&bsize,&n);
        fprintf(outPF,"# Loaded %d records of %d bits each from %s\n",n,bsize,buPO->sbname);
    }
    if(buPO->outbits && buPO->tw_num) {
        fprintf(outPF,"# Tweaking %d bits",buPO->tw_num);
        if(buPO->tw_frac > 0.0) {
            fprintf(outPF," (%5.3f fraction of %d)",buPO->tw_frac,bsize);
        }
        fprintf(outPF,"\n");
        FillRandSeedString(buPO->seed,bufS);
        fprintf(outPF,"# Random seed: %s\n",bufS);
    }
    return;
}
/*************************************************************************/
void HandleBtuOutput(BIT_UTIL *buPO, FILE *outPF)
{
    int i,b,n,min,max;
    DOUB sumD;
    char nameS[NSIZE], pformS[DEF_BS];

	HAND_NFILE(outPF);
    /***
    *   Set auto formatting and get format string 
    */
    AutoBitpoolOutFormattingI(buPO->bits);
    GetBitpoolPrintFormI(buPO->bits,pformS);
    /***
    *   Special case outputs
    */
    if(buPO->do_dump){
        DumpBitpool(buPO->bits,-1,-1,NULL,outPF);
        return;
    }
    if( buPO->do_fmat || buPO->outcomp ) {
        HandleBtuFmat(buPO, outPF);
        return;
    }
    if( buPO->outbits && (buPO->oform == BPF_TAB)){
        HandleBtuOutTable(buPO, outPF);
        return;
    }
    /***
    *   Init for stats then loop through collection
    */
    min = TOO_BIG;
    max = -TOO_BIG;
    n = 0;
    sumD = 0.0;
    if( (buPO->do_ostat) && (!buPO->outbits) && (!buPO->outcomp) ){
        fprintf(outPF,"# Name\tOn\tOff\n");
    }
    for(i=0;i<buPO->num;i++)
    {
        if(!buPO->mmask[i]) {
            continue;
        }
        n++;
        if(buPO->outbits) {
            if(buPO->pretwk) {
                DumpThisBitpoolMemberI(buPO->pretwk,i,buPO->oform,"# ",outPF);
            }
            DumpThisBitpoolMemberI(buPO->bits,i,buPO->oform,NULL,outPF);
            continue;
        }
        GetThisBitpoolNameI(buPO->bits,i,nameS);
        GetThisBitpoolOnCountI(buPO->bits, i, &b);
        sumD += DNUM(b);
        min = MIN_NUM(min,b);
        max = MAX_NUM(max,b);
        if( buPO->do_ostat) {
            fprintf(outPF,pformS,nameS);
            fprintf(outPF,"\t%d\t%d\n", b, buPO->bsize-b);
        }
    }
    /***
    *   Final story
    */
    if( buPO->do_stat) {
        fprintf(outPF,"OK records:   %d\n",n);
        if(n>0) {
            fprintf(outPF,"Min set bits: %d\n",min);
            fprintf(outPF,"Max set bits: %d\n",max);
            fprintf(outPF,"Ave set bits: %5.2f\n", sumD / DNUM(n));
        }
    }
    return;
}
/************************************************************************/
void HandleBtuOutTable(BIT_UTIL *buPO, FILE *outPF)
{
    TABLE *tabPO;

    TableFromBitpoolI(buPO->bits, buPO->mmask, &tabPO);
    DumpTable(tabPO,TRUE,FALSE,outPF);
    CHECK_TABLE(tabPO);
    return;
}
/************************************************************************
*   Full pairwise matrix of comparisons
*/
void HandleBtuFmat(BIT_UTIL *buPO, FILE *outPF)
{
    int r,c,ncol,nrow,v,op;
    char nameS[NSIZE], cnameS[NSIZE];
    BITPOOL *rbitsPO, *cbitsPO;

	HAND_NFILE(outPF);
    WriteBtuFmatHeader(buPO, outPF);
    /***
    *   Get row and col bitpools and dimensions
    *   Rows are primary members; Columns secondary members (diff or same set)
    */
    rbitsPO = buPO->bits;
    cbitsPO = (buPO->sbits) ? buPO->sbits : buPO->bits;
    GetBitpoolDimsI(rbitsPO, NULL, &nrow);
    GetBitpoolDimsI(cbitsPO, NULL, &ncol);
    /***
    *   For each primary bit set member compared to second bit set member
    *   Only primary members can be masked (i.e. mmask)
    */
    op = BOGUS;
    for(r=0;r<nrow;r++) 
    {
        if(!buPO->mmask[r]) {
            continue;
        }
        GetThisBitpoolNameI(rbitsPO,r,nameS);
        if(buPO->do_fmat) {
            fprintf(outPF,"%s",nameS);
        }
        for(c=0;c<ncol;c++) 
        {
            if( (!buPO->sbits) && (!buPO->mmask[c])) {
                continue;
            }
            if(buPO->do_band || buPO->do_b3) {
                op = BIT_AND;
            }
            else if(buPO->do_bor) {
                op = BIT_OR;
            }
            else if(buPO->do_bxor || buPO->sbits) {
                op = BIT_XOR;
            }
            else if(buPO->do_fmat) {
                op = (c<r) ? BIT_OR : BIT_AND;
            }
            if(!BitpoolBitwiseOpCountI(rbitsPO,r,cbitsPO,c,op,&v)) {
                printf("SHAM r=%d c=%d op=%d\n",r,c,op);
                return;
            }
            /***
            *   Full matrix just add number to output table; Else second name to line
            */
            if(buPO->do_fmat) {
                fprintf(outPF,"\t%d",v);
            }
            else {
                /***
                *   If single-source, only report member pairs once 
                */
                if((!buPO->sbits) && (c<r) && (!buPO->do_rap)) {
                    continue;
                }
                GetThisBitpoolNameI(cbitsPO,c,cnameS);
                fprintf(outPF,"%s\t%s\t%d",nameS,cnameS,v);
                if(buPO->do_b3) {
                    BitpoolBitwiseOpCountI(rbitsPO,r,cbitsPO,c,BIT_OR,&v);
                    fprintf(outPF,"\t%d",v);
                    BitpoolBitwiseOpCountI(rbitsPO,r,cbitsPO,c,BIT_XOR,&v);
                    fprintf(outPF,"\t%d",v);
                }
                fprintf(outPF,"\n");    
            }
        }
        if(buPO->do_fmat) {
            fprintf(outPF,"\n");    
        }
    }
    return;
}
/************************************************************************/
void WriteBtuFmatHeader(BIT_UTIL *buPO, FILE *outPF)
{
    int c,num;
    char nameS[NSIZE];
    BITPOOL *bitsPO;

	HAND_NFILE(outPF);
    fprintf(outPF,"# Full pairwise comparison matrix\n");
    if(buPO->outcomp) {
        if(buPO->do_b3) {
            fprintf(outPF,"#  Member1 Member2  AND  OR  XOR\n");
        }
        else if(buPO->do_band) {
            fprintf(outPF,"#  Member1 Member2  AND\n");
        }
        else if(buPO->do_bor) {
            fprintf(outPF,"#  Member1 Member2  OR\n");
        }
        else {
            fprintf(outPF,"#  Member1 Member3  XOR\n");
        }
    }
    /***
    *   Full matrix table header row; First or second bit pool
    */
    if(buPO->do_fmat) {
        if(!buPO->outcomp){
            fprintf(outPF,"#    AND = above diagonal (upper right)\n");
            fprintf(outPF,"#    OR  = below diagonal (lower left)\n");
        }
        fprintf(outPF,"Names");
        bitsPO = (buPO->sbits) ? buPO->sbits : buPO->bits;
        GetBitpoolDimsI(bitsPO, NULL, &num);
        for(c=0;c<num;c++) 
        {
            if( (!buPO->sbits) && (!buPO->mmask[c])) {
                continue;
            }
            GetThisBitpoolNameI(bitsPO,c,nameS);
            fprintf(outPF,"\t%s",nameS);
        }
        fprintf(outPF,"\n");
    }
    return;
}
/*************************************************************************
*   Take values in table and create a new bitpool with these
*/
int BitpoolFromTableI(TABLE *tabPO, BITPOOL **bpPPO)
{
    int r,c,nrow,ncol,v;
    char nameS[NSIZE], fnameS[NSIZE];
    DOUB vD;
    BITPOOL *bpPO;

    nrow = ncol = 0;
    GetTableDimsI(tabPO,&nrow,&ncol,FALSE);
    if((nrow * ncol) < 1) {
        return(FALSE);
    }
    if( !(bpPO = CreateBitpoolPO(ncol,nrow)) ) {
        PROBLINE;
        printf("Failed to create Bitpool row=%d col=%d\n",nrow,ncol);
        return(FALSE);
    }
    /***
    *   Copy names then party through the table
    */
    GetTableNamesI(tabPO,nameS,fnameS,NSIZE-1);
    SetBitpoolNamesI(bpPO,nameS,fnameS,-1);
    for(r=0;r<nrow;r++) 
    {
        GetTableRowLabI(tabPO,r,nameS,-1);
        SetThisBitpoolNameI(bpPO,r,nameS);
        for(c=0;c<ncol;c++) 
        {
            GetTableValI(tabPO,r,c,&vD);
            v = (vD > 0.0) ? 1 : 0;
            SetThisBitpoolBitI(bpPO,r,c,v);
        }
    }
    *bpPPO = bpPO;
    return(TRUE);
}
/*************************************************************************/
int TableFromBitpoolI(BITPOOL *bpPO, char *maskPC, TABLE **tabPPO)
{
    int r,c,m,bsize,num,nok,v;
    char nameS[NSIZE], fnameS[NSIZE];
    TABLE *tabPO;

    bsize = num = nok = 0;
    GetBitpoolDimsI(bpPO,&bsize,&num);
    nok = num;
    if(maskPC) {
        nok = NumArrayValsI(maskPC,IS_CHAR,0,num,0.1,100.0);
    }
    if((bsize * nok) < 1) {
        return(FALSE);
    }
    if( !(tabPO = CreateTablePO(num,bsize)) ) {
        PROBLINE;
        printf("Failed to create Table bsize=%d num=%d\n",bsize,num);
        return(FALSE);
    }
    /***
    *   Set names and col labels == bit number
    */
    GetBitpoolNamesI(bpPO,nameS,fnameS,-1);
    SetTableNamesI(tabPO,nameS,fnameS,-1);
    for(c=0;c<bsize;c++) 
    {
        sprintf(nameS,"Bit_%02d",c+1);
        SetTableColLabI(tabPO,c,nameS);
    }
    m = 0;
    for(r=0;r<num;r++) 
    {
        if( maskPC && (!maskPC[r]) ){
            continue;
        }
        GetThisBitpoolNameI(bpPO,r,nameS);
        SetTableRowLabI(tabPO,m,nameS);
        for(c=0;c<bsize;c++) 
        {
            GetThisBitpoolBitI(bpPO,r,c,&v);
            SetTableValI(tabPO,m,c,DNUM(v));
        }
        m++;
    }
    *tabPPO = tabPO;
    return(TRUE);
}
