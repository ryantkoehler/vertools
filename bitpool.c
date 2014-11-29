/*
* bitpool.c
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
#include "prim.h"
#include "bitpool.h"

#define DB_BITS if(DB[24])      /* Bit operations */
#define DB_BITP if(DB[25])      /* Bitpool struct */


/*****************************************************************************
*   Globals to handle bits
*/
unsigned long bitmaskGI[32] =
    {
            1 ,  1<<1,  1<<2,  1<<3,  1<<4 ,  1<<5,  1<<6,  1<<7,
         1<<8 ,  1<<9, 1<<10, 1<<11, 1<<12 , 1<<13, 1<<14, 1<<15,
        1<<16 , 1<<17, 1<<18, 1<<19, 1<<20 , 1<<21, 1<<22, 1<<23,
        1<<24 , 1<<25, 1<<26, 1<<27, 1<<28 , 1<<29, 1<<30, 2147483648
    };
short wbcountGI[256] =
    {
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
    };

/*****************************************************************************
*	Create bit-pool data struct
*/
BITPOOL *CreateBitpoolPO(int bsize, int num)
{
	BITPOOL *bpPO;

	DB_BITP DB_PrI(">> CreateBitpoolPO bsize=%d num=%d\n",bsize,num);
	if(! (bpPO = (BITPOOL *)ALLOC(1,sizeof(BITPOOL)) ) )
	{
		printf("# Failed to allocate bit-pool object\n");
	    DB_BITP DB_PrI("<< CreateBitpoolPO NULL\n");
		return(NULL);
	}
	bpPO->ID = BITPOOL_ID;
    InitBitpoolI(bpPO);
	/***
    *   Name list (empty) and call to set dims (which allocates if needed)
	*/
	bpPO->mnames = CreateWordlistPO(NULL,BP_NSIZE);
    if( !SetBitpoolDimsI(bpPO, bsize, num)) {
        CHECK_BITPOOL(bpPO);
	    DB_BITP DB_PrI("<< CreateBitpoolPO NULL\n");
        return(NULL);
    }
	DB_BITP DB_PrI("<< CreateBitpoolPO %p\n",bpPO);
	return(bpPO);
}
/*****************************************************************************
*	Free datastructure and substructs
*/
int DestroyBitpoolI(BITPOOL *bpPO)
{
	DB_BITP DB_PrI(">> DestroyBitpoolI %p\n",bpPO);
	VALIDATE(bpPO,BITPOOL_ID);
    InitBitpoolI(bpPO);
	FREE(bpPO);
	DB_BITP DB_PrI("<< DestroyBitpoolI TRUE\n");
	return(TRUE);
}
/****************************************************************************/
int InitBitpoolI(BITPOOL *bpPO)
{
	VALIDATE(bpPO,BITPOOL_ID);
    INIT_S(bpPO->name);
    INIT_S(bpPO->fname);
    bpPO->num = bpPO->n_blocks = bpPO->bsize = bpPO->bkprec = 0;
    CHECK_FREE(bpPO->blocks);
    CHECK_WORDLIST(bpPO->mnames);
    strcpy(bpPO->prlform,DEF_BP_PRLFORM_S);
	return(TRUE);
}
/*****************************************************************************
*   Set dimensions for bitpool; Number of bits per record, number of records
*   If the number doesn't fit, call to make space
*/
int SetBitpoolDimsI(BITPOOL *bpPO, int bsize, int num)
{
	VALIDATE(bpPO,BITPOOL_ID);
    
	DB_BITP DB_PrI(">> SetBitpoolDimsI %p bsize=%d num=%d\n",bpPO,bsize,num);
    /***
    *   Cannot change bit dimension if already have records
    */
    if(bsize > 0) {
        if(bpPO->num) {
            printf("Sham: bsize=%d, num already=%d\n",bsize,bpPO->num);
            ERR("SetBitpoolDimsI","Can't readjust bitsize if already have number set");
            return(FALSE);
        }
        bpPO->bsize = bsize;
        bpPO->bkprec = INT(bsize/BU_WORDLEN);
	    if( (bpPO->bkprec * BU_WORDLEN) < bsize) { 
	        bpPO->bkprec += 1;
	    }
	    DB_BITP DB_PrI("+ Set bsize=%d bkprec=%d\n", bpPO->bsize, bpPO->bkprec);
    }
    /***
    *   Number of member records; maybe have to add space
    */
    if(num > 0 ) {
        if( !HandleBitpoolSpaceI(bpPO, num)) {
            return(FALSE);
        }
	    DB_BITP DB_PrI("+ Set num=%d\n", bpPO->num);
        bpPO->num = num;
    }
	DB_BITP DB_PrI("<< SetBitpoolDimsI TRUE\n");
    return(TRUE);
}
/****************************************************************************/
int GetBitpoolDimsI(BITPOOL *bpPO, int *bsizePI, int *numPI)
{
	VALIDATE(bpPO,BITPOOL_ID);
    if(bsizePI) {
        *bsizePI = bpPO->bsize;
    }
    if(numPI) {
        *numPI = bpPO->num;
    }
    return(TRUE);
}
/*****************************************************************************
*   Check if passed record number fits in allocated space
*   Allocate if needed
*/
int HandleBitpoolSpaceI(BITPOOL *bpPO, int num)
{
    int nb;

	DB_BITP DB_PrI(">> HandleBitpoolSpaceI %p num=%d\n",bpPO,num);
	VALIDATE(bpPO,BITPOOL_ID);
    /***
    *   If index X blocks-per-rec doesn't fit, expand
    */
    if ( (bpPO->bkprec * num) >= bpPO->n_blocks ) {
        nb = MAX_NUM( (bpPO->bkprec * num), (bpPO->n_blocks + ALLOC_BLOCK) );
	    DB_BITP DB_PrI("+ Expanding; nb=%d n_blocks=%d\n",nb,bpPO->n_blocks);
        if( bpPO->n_blocks == 0 ) {
            bpPO->blocks = (BITPTR *)ALLOC(nb ,sizeof(BITPTR));
        }
        else {
            bpPO->blocks = (BITPTR *)REALLOC(bpPO->blocks, nb ,sizeof(BITPTR));
        }
        if( !bpPO->blocks ) {
            PROBLINE;
            printf("HandleBitpoolSpaceI Failed allocate %d BITPTR\n",nb);
            return(FALSE);
        }
        bpPO->n_blocks = nb;
	    DB_BITP DB_PrI("+ blocks=%p n_blocks=%d\n",bpPO->blocks, bpPO->n_blocks);
    }
	DB_BITP DB_PrI("<< HandleBitpoolSpaceI TRUE\n");
    return(TRUE);
}
/*****************************************************************************
*	Set coefficient for pool
*/
void SetBitpoolCoef(BITPOOL *bpPO,DOUB coD)
{
	VALIDATE(bpPO,BITPOOL_ID);
	bpPO->coef = coD;
}
/****************************************************************************/
int SetBitpoolNamesI(BITPOOL *bpPO, char *nameS, char *fnameS, int max)
{
    int n;

    VALIDATE(bpPO,BITPOOL_ID);
    LIMIT_NUM(max,0,NSIZE);
    if(nameS) {
        n = (max < 0) ? strlen(nameS) : max;
        strncpy(bpPO->name,nameS,n);
        bpPO->name[n] = '\0';
    }
    if(fnameS) {
        n = (max < 0) ? strlen(fnameS) : max;
        strncpy(bpPO->fname,fnameS,n);
        bpPO->fname[n] = '\0';
    }
    return(TRUE);
}
/****************************************************************************/
int GetBitpoolNamesI(BITPOOL *bpPO, char *nameS, char *fnameS, int max)
{
    int n;

    VALIDATE(bpPO,BITPOOL_ID);
    LIMIT_NUM(max,0,NSIZE);
    if(nameS) {
        n = (max < 0) ? strlen(nameS) : max;
        strncpy(nameS,bpPO->name,n);
        nameS[n] = '\0';
    }
    if(fnameS) {
        n = (max < 0) ? strlen(fnameS) : max;
        strncpy(fnameS,bpPO->fname,n);
        fnameS[n] = '\0';
    }
    return(TRUE);
}
/****************************************************************************/
int AutoBitpoolOutFormattingI(BITPOOL *bpPO)
{
    char formS[DEF_BS];

    if(!WordlistAutoFormatStringI(bpPO->mnames, NULL, formS)){
        return(FALSE);
    }
    SetBitpoolPrintFormI(bpPO, formS);
    return(TRUE);
}
/****************************************************************************/
int SetBitpoolPrintFormI(BITPOOL *bpPO, char *prlformS)
{
    VALIDATE(bpPO,BITPOOL_ID);
    strcpy(bpPO->prlform,prlformS);
    return(TRUE);
}
/****************************************************************************/
int GetBitpoolPrintFormI(BITPOOL *bpPO, char *prlformS)
{
    VALIDATE(bpPO,BITPOOL_ID);
    strcpy(prlformS,bpPO->prlform);
    return(TRUE);
}
/***************************************************************************/
int CopyThisBitpoolMemberI(BITPOOL *fbPO, int fm, BITPOOL *sbPO, int sm)
{
    int ok,i,n;
    BITPTR *fbitsPI, *sbitsPI;
    char nameS[NSIZE];

    /***
    *   Same bsize, and get pointers to first and sec bit members
    */
    if(!SameBitpoolDimsI(fbPO,sbPO,TRUE,FALSE)) {
        return(FALSE);
    }
    ok = GetThisBitpoolPtrBitsI(fbPO, fm, &fbitsPI, &n);
    ok *= GetThisBitpoolPtrBitsI(sbPO, sm, &sbitsPI, NULL);
    if(!ok){
        return(FALSE);
    }
    /***
    *   Copy bit blocks then name
    */
    for(i=0;i<n;i++) 
    {
        sbitsPI[i] = fbitsPI[i];
    }
    ok = GetThisBitpoolNameI(fbPO, fm, nameS);
    ok *= SetThisBitpoolNameI(sbPO, sm, nameS);
    return(ok);
}
/****************************************************************************
*   Add name and bit values to m'th record; 
*   Extends bitpool length (and calls to add space) if needed
*/
int AddThisBitpoolMemberI(BITPOOL *bpPO, int m, char *bitS, char *nameS)
{
    int ok;

	DB_BITP DB_PrI(">> AddThisBitpoolMemberI %p m=%d\n",bpPO,m);
	VALIDATE(bpPO,BITPOOL_ID);
    ok = TRUE;
    if( m >= bpPO->num ) {
	    DB_BITP DB_PrI("+ %d >= %d so setting dims %d\n",m,bpPO->num,m+1);
        ok = SetBitpoolDimsI(bpPO, -1, m+1);
    }
    if(ok) {
	    DB_BITP DB_PrI("+ Setting member m=%d\n",m);
        ok = SetThisBitpoolMemberI(bpPO, m, bitS, nameS);
    }
	DB_BITP DB_PrI("<< AddThisBitpoolMemberI %d\n",ok);
    return(ok);
}
/****************************************************************************
*	Set name and bit values for m'th record in pool
*/
int SetThisBitpoolMemberI(BITPOOL *bpPO, int m, char *bitS, char *nameS)
{
	DB_BITP DB_PrI(">> SetThisBitpoolMemberI %p m=%d\n",bpPO,m);
	VALIDATE(bpPO,BITPOOL_ID);
    if(nameS) {
	    DB_BITP DB_PrI("+ Setting name m=%d |%s|\n",m,nameS);
        if(!SetThisBitpoolNameI(bpPO, m, nameS)){
            printf("Failed to set %p [%d] name |%s|\n",bpPO,m,nameS);
            return(FALSE);
        }
    }
    if(bitS) {
	    DB_BITP DB_PrI("+ Setting bits m=%d |%s|\n",m,bitS);
        if(!SetThisBitpoolBitstringI(bpPO, m, bitS)){
            printf("Failed to set %p [%d] bits |%s|\n",bpPO,m,bitS);
            return(FALSE);
        }
    }
	DB_BITP DB_PrI("<< SetThisBitpoolMemberI TRUE\n");
	return(TRUE);
}
/****************************************************************************/
int SetThisBitpoolNameI(BITPOOL *bpPO, int m, char *nameS)
{
    int ok;

	VALIDATE(bpPO,BITPOOL_ID);
    ok = TRUE;
    if(!NO_S(nameS)){
        ok = AddWordlistWordI(bpPO->mnames,m,nameS,BP_NSIZE);
    }
    return(ok);
}
/****************************************************************************/
int GetThisBitpoolNameI(BITPOOL *bpPO, int m, char *nameS)
{
	VALIDATE(bpPO,BITPOOL_ID);
    return(GetWordlistWordI(bpPO->mnames,m,nameS,BP_NSIZE));
}
/****************************************************************************
*   Set string of 01 to actual bits for member of bpool
*/
int SetThisBitpoolBitstringI(BITPOOL *bpPO, int m, char *bitS)
{
	int b;
	BITPTR *bitsPI;
    char *cPC;

	DB_BITP DB_PrI(">> SetThisBitpoolBitstringI %p m=%d\n",bpPO,m);
    if(!GetThisBitpoolPtrBitsI(bpPO, m, &bitsPI, NULL)) {
	    DB_BITP DB_PrI("<< SetThisBitpoolBitstringI no bit Ptr FALSE\n");
        return(FALSE);
    }
	b = 0;
    cPC = bitS;
    while( ISLINE(*cPC) && (b<bpPO->bsize) )
	{
        if(*cPC =='0'){
			SetThisBitI(bitsPI,b++,FALSE);	
        }
        else if(*cPC =='1'){
			SetThisBitI(bitsPI,b++,TRUE);	
        }
        else if( isgraph(INT(*cPC)) ) {
            PROBLINE;
            printf("Setting mem[%d] bit[%d] not 0 or 1 in |%s|\n",m,b+1,bitS);
            return(FALSE);
        }
        cPC++;
    }
	DB_BITP DB_PrI("<< SetThisBitpoolBitstringI %d\n",b);
    return(b);
}
/****************************************************************************
*   Set specific bit [b] in member [m] of bitpool to value of on
*/
int SetThisBitpoolBitI(BITPOOL *bpPO, int m, int b, int on)
{
    int ok,bsize,num;
	BITPTR *bitsPI;

    ok = FALSE;
    if( GetBitpoolDimsI(bpPO,&bsize,&num)) {
        if( (m>=0) && (m<num) && (b>=0) && (b<bsize) ) {
            if(GetThisBitpoolPtrBitsI(bpPO, m, &bitsPI,NULL)) {
                SetThisBitI(bitsPI,b,on);
                ok++;
            }
        }
    }
    return(ok);
}
/****************************************************************************
*	Set pointer to bitstring for member m of pool
*/
int GetThisBitpoolPtrBitsI(BITPOOL *bpPO, int m, BITPTR **bPPI, int *nPI)
{
    int ok;

	VALIDATE(bpPO,BITPOOL_ID);
    ok = 0;
    if((m >= 0) && ( m < bpPO->num)){
        ok++;
        if(bPPI){
	        *bPPI = &bpPO->blocks[m * bpPO->bkprec];
        }
    }
    if(ok && nPI) {
        *nPI = bpPO->bkprec;
    }
	return(ok);
}
/****************************************************************************/
int GetThisBitpoolBitI(BITPOOL *bpPO, int m, int b, int *onPI)
{
    int ok,bsize,num;
	BITPTR *bitsPI;

    ok = FALSE;
    if( GetBitpoolDimsI(bpPO,&bsize,&num)) {
        if( (m>=0) && (m<num) && (b>=0) && (b<bsize) ) {
            if(GetThisBitpoolPtrBitsI(bpPO, m, &bitsPI, NULL)) {
                *onPI = GetThisBitI(bitsPI,b);
                ok++;
            }
        }
    }
    return(ok);
}
/***************************************************************************/
int GetThisBitpoolOnCountI(BITPOOL *bpPO, int m, int *onPI)
{
    int bsize;
	BITPTR *bitsPI;

    if( GetBitpoolDimsI(bpPO, &bsize, NULL)) {
        if(GetThisBitpoolPtrBitsI(bpPO, m, &bitsPI, NULL)) {
            *onPI = OnBitCountI(bitsPI, bsize);
            return(TRUE);
        }
    }
    return(FALSE);
}
/***************************************************************************/
int BitpoolBitwiseOpCountI(BITPOOL *fbPO, int fm, BITPOOL *sbPO, int sm, int op, 
    int *onPI)
{
    int fsize, ssize, on;
    BITPTR *fbitsPI, *sbitsPI;

    fsize = ssize = 0;
    GetBitpoolDimsI(fbPO, &fsize, NULL);
    GetBitpoolDimsI(sbPO, &ssize, NULL);
    /***
    *   Bitsizes must be non-zero and the same
    */
    if( ((fsize * ssize)<1) || (fsize != ssize) ) {
printf("xxx fsize=%d ssize=%d\n",fsize,ssize);
        return(FALSE);
    }
    if(!GetThisBitpoolPtrBitsI(fbPO, fm, &fbitsPI, NULL)) {
printf("xxx no ptr for %d\n",fm);
        return(FALSE);
    }
    if(!GetThisBitpoolPtrBitsI(sbPO, sm, &sbitsPI, NULL)) {
printf("xxx no ptr for %d\n",sm);
        return(FALSE);
    }
    on = BitwiseOpCountI(fbitsPI, sbitsPI, fsize, op); 
    if(on < 0) {
printf("xxx on = %d\n",on);
        return(FALSE);
    }
    *onPI = on;
    return(TRUE);
}
/***************************************************************************/
int ModThisBitpoolI(BITPOOL *bpPO, int st, int en, int op) 
{
    int i,num;

    GetBitpoolDimsI(bpPO, NULL, &num);
    for(i=0;i<num;i++) {
        if(!ModThisBitpoolBitsI(bpPO, i, st, en, op)) {
            return(FALSE);
        } 
    }
    return(TRUE);
}
/***************************************************************************/
int ModThisBitpoolBitsI(BITPOOL *bpPO, int m, int st, int en, int op) 
{
    int i,v,bsize;

    GetBitpoolDimsI(bpPO, &bsize, NULL);
    st = (st<0) ? 0 : st;
    en = (en<0) ? bsize: en;
    LIMIT_NUM(st,0,bsize);
    LIMIT_NUM(en,0,bsize);
    for(i=st;i<en;i++) 
    {
        if(!GetThisBitpoolBitI(bpPO,m,i,&v)) {
            return(FALSE);
        }
        switch(op)
        {
            case BIT_NOT:   v = (v) ? FALSE : TRUE;   break;
            case BIT_ON:    v = TRUE;  break;
            case BIT_OFF:   v = FALSE; break;
            default:
                printf("Bad opcode=%d\n",op);
                ERR("ModThisBitpoolBitsI","Bad op code");
                return(FALSE);
        }
        SetThisBitpoolBitI(bpPO,m,i,v);
    }
    return(TRUE);
}
/***************************************************************************/
int SameBitpoolDimsI(BITPOOL *fbPO, BITPOOL *sbPO, int size, int num)
{
    int fsize, ssize, fnum, snum;

    GetBitpoolDimsI(fbPO, &fsize, &fnum);
    GetBitpoolDimsI(sbPO, &ssize, &snum);
    if( size && (fsize != ssize)) {
        return(FALSE);
    }
    if( num && (fnum != snum)) {
        return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*	Dump pool members of bitstring
*/
void DumpBitpool(BITPOOL *bpPO, int st, int en, char *preS, FILE *outPF)
{
    DumpBitpoolDescription(bpPO, preS, outPF);
    DumpBitpoolMembers(bpPO, st, en, preS, outPF);
}
/**************************************************************************/
void DumpBitpoolDescription(BITPOOL *bpPO, char *pS, FILE *outPF)
{
    int bsize,num;
    char nameS[NSIZE], fnameS[NSIZE], preS[DEF_BS];

    VALIDATE(bpPO,BITPOOL_ID);
    HAND_NFILE(outPF);
    INIT_S(preS);
    if(pS) {
        sprintf(preS,"%s",pS);
    }
    fprintf(outPF,"%s",preS);
    fprintf(outPF,"Bitpool at %p\n",bpPO);
    GetBitpoolNamesI(bpPO,nameS,fnameS,NSIZE-1);
    if(!NO_S(nameS)) {
        fprintf(outPF,"%s",preS);
        fprintf(outPF,"Bitpool name: %s\n",nameS);
    }
    if(!NO_S(fnameS)) {
        fprintf(outPF,"%s",preS);
        fprintf(outPF,"Bitpool file: %s\n",nameS);
    }
    GetBitpoolDimsI(bpPO, &bsize, &num);
    fprintf(outPF,"%s",preS);
    fprintf(outPF,"%d Members with %d bits each\n",num,bsize);
}
/**************************************************************************/
void DumpBitpoolMembers(BITPOOL *bpPO, int st, int en, char *pS, FILE *outPF)
{
	BITPTR *bitsPI;
	int i,n,num,bsize;
    char bufS[NSIZE], preS[DEF_BS];

	VALIDATE(bpPO,BITPOOL_ID);
    HAND_NFILE(outPF);
    INIT_S(preS);
    if(pS) {
        sprintf(preS,"%s",pS);
    }
    GetBitpoolDimsI(bpPO, &bsize, &num);
    st = (st<0) ? 0 : st;
    en = (en<0) ? num : en;
    LIMIT_NUM(st,0,num);
    LIMIT_NUM(en,0,num);
	for(i=st;i<en;i++)
	{
        GetThisBitpoolPtrBitsI(bpPO, i, &bitsPI, NULL);
		n = OnBitCountI(bitsPI,bpPO->bsize);
        sprintf(bufS,"%s[%d]\t(%d)\t",preS,i,n);
        DumpThisBitpoolMemberI(bpPO, i, BPF_BITS, bufS, outPF);
	}
    return;
}
/**************************************************************************/
int DumpThisBitpoolMemberI(BITPOOL *bpPO, int m, int oform, char *pS, FILE *outPF)
{
	BITPTR *bitsPI;
	int i,on,nb;
    char nameS[NSIZE], preS[DEF_BS];

	VALIDATE(bpPO,BITPOOL_ID);
    HAND_NFILE(outPF);
    INIT_S(preS);
    if(pS) {
        sprintf(preS,"%s",pS);
    }
    if(!GetThisBitpoolNameI(bpPO, m, nameS)) {
        return(FALSE);
    }
    if(!GetThisBitpoolPtrBitsI(bpPO, m, &bitsPI, &nb)) {
        return(FALSE);
    }
    switch(oform) 
    {
        case BPF_BITS:
        case BPF_SBITS:
        case BPF_HEX:
		    fprintf(outPF,"%s",preS);
		    fprintf(outPF,bpPO->prlform,nameS);
            if(oform == BPF_HEX) {
                for(i=0;i<nb;i++) 
                {
		            fprintf(outPF,"\t%x",(unsigned int)bitsPI[i]);
                }
            }
            else {
		        fprintf(outPF,"\t");
                for(i=0;i<bpPO->bsize;i++) 
                {
                    GetThisBitpoolBitI(bpPO, m, i, &on);
                    if(on) {
		                fprintf(outPF,"1");
                    }
                    else {
		                fprintf(outPF,"0");
                    }
                    if(oform == BPF_SBITS) {
		                fprintf(outPF," ");
                    }
                }
            }
		    fprintf(outPF,"\n");
            break;
        default:
            printf("Bad oform=%d\n",oform);
            ERR("DumpThisBitpoolMemberI","Bad out format code");
            return(FALSE);
    }
    return(TRUE);
}
/***************************************************************************
*   Figure out bit file type given options and input name
*   Return format code or false if error
*/
int FigureBitFileTypeI(int ibit, int isbit, int ihex, int itab, char *fnameS, int error)
{
    int form;

    if(ibit) {
        form = BPF_BITS;
    }
    else if(isbit) {
        form = BPF_SBITS;
    }
    else if(ihex) {
        form = BPF_HEX;
    }
    else if(itab) {
        form = BPF_TAB;
    }
    else {
        if( STDIN_STR(fnameS) ) {
            form = BPF_BITS;
        }
        else {
            form = GuessBitFileTypeI(fnameS, error);
        }
    }
    if(!OkBitFileFormatI(form,fnameS,error)) {
        return(FALSE);
    }
    return(form);
}
/***************************************************************************/
int GuessBitFileTypeI(char *nameS, int error)
{
    int form;
    char extS[DEF_BS];

    GetFilePartsI(nameS,NULL,NULL,extS);
    form = 0;
    switch(extS[0]) {
        case 'b': case 'B': form = BPF_BITS;    break;
        case 'h': case 'H': form = BPF_HEX;     break;
        case 't': case 'T': form = BPF_TAB;     break;
    }
    return(form);
}
/***************************************************************************/
int OkBitFileFormatI(int form, char *nameS, int error)
{
    int ok;

    ok = FALSE;
    switch(form) {
        case BPF_BITS:
        case BPF_SBITS:
        case BPF_HEX:
        case BPF_TAB:
            ok++;   break;
        default:
            if(error) {
                PROBLINE;
                printf("Unrecognized file format (%d)",form);
                if(nameS) {
                    printf(" |%s|",nameS);
                }
                printf("\n");
            }
    }
    return(ok);
}
/***************************************************************************/
int OkBitwiseOpI(int op, char *opS)
{
    int ok;
    char bufS[DEF_BS];

    INIT_S(bufS);
    ok = FALSE;
    switch(op) 
    {
        case BIT_AND:   ok++;   sprintf(bufS,"bitwise AND");    break;
        case BIT_OR:    ok++;   sprintf(bufS,"bitwise OR");    break;
        case BIT_XOR:   ok++;   sprintf(bufS,"bitwise XOR");    break;
    }
    if(ok && opS) {
        sprintf(opS,"%s",bufS);
    }
    return(ok);
}
/****************************************************************************
*	Attempt to parse bitpool from named file.
*/
int GetBitpoolI(char *fileS, int iform, int error, BITPOOL **bpPPO)
{
    int n;
	FILE *fPF;
    BITPOOL *newPO;
    char nameS[NSIZE];

	DB_BITP DB_PrI(">> GetBitpoolI |%s| iform=%d\n",fileS,iform);
	if(!(fPF=OpenUFilePF(fileS,"r",NULL))) {
	    DB_BITP DB_PrI("<< GetBitpoolI FALSE\n");
		return(FALSE);
	}
    n = ParseBitpoolI(fPF, iform, error, &newPO);
    if (n) {
	    GetFilePartsI(fileS,NULL,nameS,NULL);
        SetBitpoolNamesI(newPO, nameS, fileS, NSIZE-1);
        *bpPPO = newPO;
    }
	CHECK_FILE(fPF);
	DB_BITP DB_PrI("<< GetBitpoolI %d\n",n);
    return(n);
}
/***************************************************************************
*   Parse bit pool from given file
*/
int ParseBitpoolI(FILE *inPF, int iform, int error, BITPOOL **bpPPO)
{
    int ok;

    ok = FALSE;
    switch(iform) {
		case BPF_BITS:
		case BPF_SBITS:
            ok = ParseBitpoolOneZeroI(inPF,error,bpPPO);
            break;
		default:
			printf("Bad format code=%d\n",iform);
			ERR("ParseBitpoolI","Bad input format code");
            break;
    }
    return(ok);
}
/***************************************************************************/
int ParseBitpoolOneZeroI(FILE *inPF, int error, BITPOOL **bpPPO)
{
	int s,b,m,line;
	char bufS[BBUFF_SIZE], nameS[NSIZE], *cPC;
	BITPOOL *bpPO;

	DB_BITP DB_PrI(">> ParseBitpoolOneZeroI\n");
	/***
	*	Create empty structure
	*/
	if(!(bpPO=CreateBitpoolPO(0,0))) {
		return(FALSE);
	}
	/***
	*	Load records
	*/
	b = s = m = line = 0;
	while(fgets(bufS,BLINEGRAB,inPF) != NULL)
	{
        line++;
        if( (COM_LINE(bufS)) || BlankStringI(bufS) ) {
			continue;
		}
        cPC = bufS;
        PASS_BLANK(cPC);            /* Pass any spaces before name */
		INIT_S(nameS);
		sscanf(bufS,"%s",nameS);
        PASS_WORD(cPC);             /* Pass name */
        /***
        *   First record deterines bitstring dimension, else make sure same
        *   Otherwise 
        */
        if(!BitStringOneZeroLineDimI(cPC, &s)) {
            m = 0;
            break;
        }
        if( m==0 ) {
            SetBitpoolDimsI(bpPO, s, -1);
            b = s;
        }
        else if (s != b) {
            if(error) {
                PROBLINE;
                printf("# Expecting %d bits, found %d; Line %d |%s|\n",b,s,line,bufS);
                m = 0;
                break;
            }
            else {
                WARNLINE;
                Chomp(bufS);
                printf("# Expecting %d bits, found %d; Ignoring line %d |%s|\n",b,s,line,bufS);
                continue;
            }
        }
        /***
        *   Add to collection
        */
        if( !AddThisBitpoolMemberI(bpPO, m, cPC, nameS)) {
            PROBLINE;
            printf("Failed to add [%d] bitstring. Line %d |%s|\n",m,line,bufS);
            m = 0;
            break;
        }
		m++;
	}
	/***
	*	Clean up or set pointer and all done
	*/
    if( m < 1 ) {
        CHECK_BITPOOL(bpPO);    
    }
    else {
	    *bpPPO = bpPO;
    }
	DB_BITP DB_PrI("<< ParseBitpoolI %d\n",m);
	return(m);	
}
/***************************************************************************/
int BitStringOneZeroLineDimI(char *lineS, int *bsizePI)
{
    int n;
    char *cPC;

    cPC = lineS;
    n = 0;
    while( ISLINE(*cPC) ) {
        if( (*cPC =='0') || (*cPC =='1') ){
            n++;
        }
        else if( isgraph(INT(*cPC)) ) {
            PROBLINE;
            printf("Char not 0 or 1; %d of |%s|\n",n+1,lineS);
            return(FALSE);
        }
        cPC++;
    }
    if(bsizePI) {
        *bsizePI = n;
    }
    return(TRUE);
}
/****************************************************************************
*	Set bit value in bitstring bitsPI
*	No dimension checking so be carefull what you pass these guy
*/
int SetThisBitI(BITPTR *bitsPI,int bit,int on)
{
	int w,b;

	w = bit / BU_WORDLEN;
	b = bit % BU_WORDLEN;
	if(on) {
		bitsPI[w] |= bitmaskGI[b];
	}
	else {
		bitsPI[w] &= ~bitmaskGI[b];
	}
    return(TRUE);
}
/*************************************************************************/
int GetThisBitI(BITPTR *bitsPI,int bit)
{
	int w,b,on;

	w = bit / BU_WORDLEN;
	b = bit % BU_WORDLEN;
	on = (bitsPI[w] & bitmaskGI[b]) ? 1 : 0;
    return(on);
}
/*************************************************************************/
int ArrayToBitsI(int *aPI, int len, BITPTR *bitsPI)
{
	int i,n,w,b;

	DB_BITS DB_PrI(">> ArrayToBitsI %p n=%d %p\n",aPI,len,bitsPI);
	n = 0;
	for(i=0;i<len;i++)
	{
		w = i / BU_WORDLEN;
		b = i % BU_WORDLEN;
		DB_BITS DB_PrI("+  W[%d]B[%d] ",w,b); 
		if(aPI[i]) {
			bitsPI[w] |= bitmaskGI[b];
			DB_BITS DB_PrI("ON\n");
			n++;
		}
		else {
			bitsPI[w] &= ~bitmaskGI[b];
			DB_BITS DB_PrI("off\n");
		}
	}
	DB_BITS DB_PrI("<< ArrayToBitsI %d\n",n);
	return(n);
}
/*******************************************************************/
int BitsToArrayI(BITPTR *bitsPI, int len, int clean, int *aPI)
{
	int i,n,w,b;

	DB_BITS DB_PrI(">> BitsToArrayI %p n=%d %d %p\n",bitsPI,len,clean,aPI);
	n = 0;
	for(i=0;i<len;i++)
	{
		w = i / BU_WORDLEN;
		b = i % BU_WORDLEN;
		DB_BITS 
		{
			if(b==0) {
				DB_PrI("+  W = %x\n",bitsPI[w]);
            }
			DB_PrI("+  W[%d]B[%d] ",w,b);
		}
		if(bitsPI[w] & bitmaskGI[b]) {
			DB_BITS DB_PrI("ON\n");
			if(clean) {
				aPI[i] = 1;
            }
			else {
				aPI[i] += 1;
            }
			n++;
		}
		else {
			if(clean) {
				aPI[i] = 0;
            }
			DB_BITS DB_PrI("off\n");
		}
	}
	DB_BITS DB_PrI("<< BitsToArrayI %d\n",n);
	return(n);
}
/*******************************************************************/
int OnBitCountI(BITPTR *bitsPI, int len)
{
	int i,j,n,nb;
	unsigned char *tbitsPC;

	DB_BITS DB_PrI(">> OnBitCountI %p len=%d\n",bitsPI,len);
	n = 0;
	/***
	*	Block at a time
	*/
	nb = (len / BU_BLK_SIZE);
	DB_BITS DB_PrI("+ %d blocks\n",nb);
	tbitsPC = (unsigned char *)bitsPI;
	for(i=0;i<nb;i++)
	{
		n += wbcountGI[INT(tbitsPC[i])];
		DB_BITS DB_PrI("+ [%d] %3d, u = %d, n = %d\n",
			i,INT(tbitsPC[i]),wbcountGI[INT(tbitsPC[i])],n);
	}
	/***
	*	One at a time
	*/
	nb = len % BU_BLK_SIZE;
	DB_BITS DB_PrI("+ %d remaining bits\n",nb);
	for(j=0;j<nb;j++)
	{
		DB_BITS DB_PrI("+  j = %d ",j);
		if(tbitsPC[i] & bitmaskGI[j]) {
			n++;
			DB_BITS DB_PrI("ON = %d",n);
		}
		DB_BITS DB_PrI("\n");
	}
	DB_BITS DB_PrI("<< OnBitCountI %d\n",n);
	return(n);
}
/*******************************************************************/
int DifBitCountI(BITPTR *bitsPI, BITPTR *sbitsPI, int len)
{
    return(BitwiseOpCountI(bitsPI, sbitsPI, len, BIT_XOR));
}
/********************************************************************/
int BitwiseOpCountI(BITPTR *bitsPI, BITPTR *sbitsPI, int len, int op)
{
	int i,n,nb,no;
	unsigned char *tbitsPC, *stbitsPC;
	BITPTR last;
	
	tbitsPC = (unsigned char *)bitsPI;
	stbitsPC = (unsigned char *)sbitsPI;
	n = 0;
	/***
	*	Block (8-bit) at a time for first bits
	*/
	no = len / BU_WORDLEN;
	nb = no * BU_WORDLEN / BU_BLK_SIZE;
	for(i=0;i<nb;i++)
	{
        switch(op) {
            case BIT_AND:
		        n += wbcountGI[INT(tbitsPC[i] & stbitsPC[i])];
                break;
            case BIT_OR:
		        n += wbcountGI[INT(tbitsPC[i] | stbitsPC[i])];
                break;
            case BIT_XOR:
		        n += wbcountGI[INT(tbitsPC[i] ^ stbitsPC[i])];
                break;
            default:
                printf("Bad bitwise op code=%d\n",op);
                ERR("BitLogicCountI","Bad op code");
                return(BOGUS);
        }
	}
	/***
	*	One bit at a time for last bits
	*/
    last = 0;
    switch(op) {
        case BIT_AND:
	        last = bitsPI[no] & sbitsPI[no];    break;
        case BIT_OR:
	        last = bitsPI[no] | sbitsPI[no];    break;
        case BIT_XOR:
	        last = bitsPI[no] ^ sbitsPI[no];    break;
    }
	nb = len - (nb * BU_BLK_SIZE);
	for(i=0;i<nb;i++)
	{
		if(last & bitmaskGI[i]) {
			n++;
        }
	}
	return(n);
}
/*******************************************************************
*/
int DifByMaxBitsI(BITPTR *bitsPI, BITPTR *sbitsPI, int len, int max)
{
	int i,n,nb;
	unsigned char *tbitsPC, *stbitsPC, lastC;
	
	DB_BITS DB_PrI(">> DifByMaxBitsI len %d max %d\n",len,max);
	n = 0;
	/***
	*	Block at a time
	*/
	nb = len / BU_BLK_SIZE;
	tbitsPC = (unsigned char *)bitsPI;
	stbitsPC = (unsigned char *)sbitsPI;
	for(i=0;i<nb;i++)
	{
		DB_BITS DB_PrI("+ [%d] ",i);
		n += wbcountGI[INT(tbitsPC[i] ^ stbitsPC[i])];
		DB_BITS DB_PrI("%d\n",n);
	}
	if(n>max) {
		DB_BITS DB_PrI("<< DifByMaxBitsI TRUE\n");
		return(TRUE);
	}
	/***
	*	One at a time
	*/
	lastC = tbitsPC[i] ^ stbitsPC[i];
	nb = len % BU_BLK_SIZE;
	for(i=0;i<nb;i++)
	{
		DB_BITS DB_PrI("+ [%d] ",i);
		n += (lastC & bitmaskGI[i]);
		DB_BITS DB_PrI("%d\n",n);
	}
	if(n>max) {
		DB_BITS DB_PrI("<< DifByMaxBitsI TRUE\n");
		return(TRUE);
	}
	DB_BITS DB_PrI("<< DifByMaxBitsI FALSE\n");
	return(FALSE);
}
/************************************************************************
*	Returns Tanamoto similarity coefficient bewteen two bit strings up to len
*/
REAL BitTanCoefR(BITPTR *bitsPI, BITPTR *sbitsPI, int len)
{
	REAL coR;

	int i,j,nf,ns,nc,b,nb,den;
	unsigned char *tbitsPC, *stbitsPC;

	DB_BITS DB_PrI(">> BitTanCoefR %p %p %d\n",bitsPI,sbitsPI,len);
	/***
	*	Block at a time
	*/
	nb = len / BU_BLK_SIZE;
	DB_BITS DB_PrI("+ %d blocks\n",nb);
	tbitsPC = (unsigned char *)bitsPI;
	stbitsPC = (unsigned char *)sbitsPI;
	nf = ns = nc = 0;
	for(i=0;i<nb;i++)
	{
		nf += wbcountGI[INT(tbitsPC[i])];
		ns += wbcountGI[INT(stbitsPC[i])];
		nc += wbcountGI[INT(tbitsPC[i] & stbitsPC[i])];
	}
	/***
	*	One at a time
	*/
	nb = len % BU_BLK_SIZE;
	DB_BITS DB_PrI("+ %d remaining bits\n",nb);
	for(j=0;j<nb;j++)
	{
		b = j % BU_WORDLEN;
		if(tbitsPC[i] & bitmaskGI[b]) {
			nf++;
        }
		if(stbitsPC[i] & bitmaskGI[b]) {
			ns++;
        }
		if((tbitsPC[i] & stbitsPC[i]) & bitmaskGI[b]) {
			nc++;
        }
	}
	DB_BITS DB_PrI("+ nf = %d, ns = %d, nc = %d\n",nf,ns,nc);
	den = nf + ns - nc;
	if(den > 0) {
		coR = RNUM(nc)/RNUM(den);
    }
	else {
		coR = 1.0;
    }
	DB_BITS DB_PrI("<< BitTanCoefR %f\n",coR);
	return(coR);
}
