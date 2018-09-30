/*
* sctab_part.c
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
#include <math.h>
#include "prim.h"
#include "score.h"
#include "table.h"
#include "scoretab.h"

#define DB_PART if(DB[75])


/**************************************************************************
*   Handle 
*/
int HandleSetPartitioningI(SCORETAB *stPO, TABLE *tabPO)
{
    int p,npools;
    TABLE *poolsPO, *dupPO;
    char dnameS[NSIZE];

    DB_PART DB_PrI(">> HandleSetPartitioningI\n");
    /***
    *   Symmetric?
    */
    if(tabPO->nrow != tabPO->ncol) {
        PROBLINE;
        printf("Partitioning only works on symmetric matrices\n");
        DB_PART DB_PrI("<< HandleSetPartitioningI FALSE\n");
        return(FALSE);
    }
    /***
    *   How many pools?
    */
    npools = tabPO->nrow / stPO->psize;
    if( (npools * stPO->psize) < tabPO->nrow ) {
        npools++;
    }
    WritePartitionHeader(stPO,tabPO,npools,stPO->out);
    /***
    *   Get working partition table and intialize
    *   Each row = pool partition
    */
    if(! (poolsPO = CreateTablePO(npools,stPO->psize))) {
        return(FALSE);
    }
    SetStartingPoolsI(stPO,tabPO,poolsPO,npools,stPO->psize);
    DB_PART DumpTable(poolsPO,TRUE,TRUE,NULL);
    /***
    *   If dumping matrices, create copy to mess with and write out starting
    */
    dupPO = NULL;
    if(!NO_S(stPO->dumpbase))
    {
        if(!CopyTableI(tabPO,TRUE,FALSE,&dupPO)) {
            CHECK_TABLE(poolsPO);
            return(FALSE);
        }
        MaskPoolMemValsI(poolsPO,tabPO,dupPO,FALSE);
        sprintf(dnameS,"%s-beg-part.mat",stPO->dumpbase);
        WriteTableFileI(dupPO,TRUE,TRUE,dnameS);
        MaskPoolMemValsI(poolsPO,tabPO,dupPO,TRUE);
        sprintf(dnameS,"%s-beg-pool.mat",stPO->dumpbase);
        WriteTableFileI(dupPO,TRUE,TRUE,dnameS);
    }
    /***
    *   
    */
    for(p=0;p<npools;p++)
    {
        ScreenPoolI(poolsPO,p,tabPO,stPO->pmin,stPO->out);
    }
    for(p=0;p<npools;p++)
    {
        ScanShamsI(poolsPO,p,tabPO,stPO->pmin);
        DumpPoolMembersI(poolsPO,p,tabPO,stPO->out);
    }
    DB_PART DumpTable(poolsPO,TRUE,TRUE,NULL);
    /***
    *   Dumping results?
    */
    if( ! NO_S(stPO->dumpbase) ) {
        MaskPoolMemValsI(poolsPO,tabPO,dupPO,FALSE);
        sprintf(dnameS,"%s-end-part.mat",stPO->dumpbase);
        WriteTableFileI(dupPO,TRUE,TRUE,dnameS);
        MaskPoolMemValsI(poolsPO,tabPO,dupPO,TRUE);
        sprintf(dnameS,"%s-end-pool.mat",stPO->dumpbase);
        WriteTableFileI(dupPO,TRUE,TRUE,dnameS);
    }
    CHECK_TABLE(dupPO);
    CHECK_TABLE(poolsPO);
    DB_PART DB_PrI("<< HandleSetPartitioningI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Set intial pool membership table    
*/
int SetStartingPoolsI(SCORETAB *stPO,TABLE *tabPO, TABLE *poolsPO, int npools,
    int psize)
{
    int r,p,m;
    DOUB vD;

    DB_PART DB_PrI(">> SetStartingPoolsI\n");
    InitTableValsI(poolsPO,-1.0,TRUE);
    /***
    *   Simply in order
    */
    r = 0;
    for(p=0;p<npools;p++)
    {
        for(m=0;m<psize;m++)
        {
            vD = RNUM(r);
            SetTableValI(poolsPO,p,m,vD);
            r++;
            if( r >= tabPO->nrow) {
                break;
            }
        }
    }
    DB_PART DB_PrI("<< SetStartingPoolsI\n");
    return(TRUE);
}
/**************************************************************************
*   Get table index for pool member
*/
int GetPoolMemTabIndexI(TABLE *poolsPO, int p1, int m1, int *indPI)
{
    int i;
    DOUB vD;

    if( ! GetTableValI(poolsPO,p1,m1,&vD) ) {
        return(FALSE);
    }
    i = INT(vD);
    if( i<0 ) {
        return(FALSE);
    }
    *indPI = i;
    return(TRUE);
}
/**************************************************************************
*   Get value associated with members of pools
*/
int GetPoolPairValueI(TABLE *poolsPO, int p1, int m1, int p2, int m2,
    TABLE *tabPO, DOUB *vPD)
{
    int i,j;
    DOUB vD;

    DB_PART DB_PrI(">> GetPoolPairValueI p1=%d m1=$d p2=%d m2=%d\n",p1,m1,p2,m2);
    if(!GetPoolMemTabIndexI(poolsPO, p1, m1, &i)) {
        DB_PART DB_PrI("<< GetPoolPairValueI FALSE\n");
        return(FALSE);
    }
    if(!GetPoolMemTabIndexI(poolsPO, p2, m2, &j)) {
        DB_PART DB_PrI("<< GetPoolPairValueI FALSE\n");
        return(FALSE);
    }
    if(!GetTableValI(tabPO,i,j,&vD)) {
        return(FALSE);
    }
    *vPD = vD;
    DB_PART DB_PrI("<< GetPoolPairValueI vD=%f TRUE\n",vD);
    return(TRUE);
}
/**************************************************************************
*   Screen pool pairwise members
*/
int ScreenPoolI(TABLE *poolsPO,int pool,TABLE *tabPO,DOUB minD, FILE *outPF)
{
    int i,j,m1,m2;
    DOUB vD;

    DB_PART DB_PrI(">> ScreenPoolI pool=%d min=%f\n",pool,minD);
    /***
    *   Each pair i-j
    */
    for(i=0;i<poolsPO->ncol;i++)
    {
        for(j=i+1;j<poolsPO->ncol;j++)
        {
            if( ! GetPoolPairValueI(poolsPO,pool,i,pool,j,tabPO,&vD) ) {
                continue;
            }
            DB_PART DB_PrI("+ pair[%d][%d] = %f\n",m1,m2,vD);
            /***
            *   Happy pair so ignore
            */
            if(vD > minD) {
                continue;
            }
            /***
            *   Have a forbidden conflict so must swap or kill
            */
            GetPoolMemTabIndexI(poolsPO, pool, i, &m1);
            GetPoolMemTabIndexI(poolsPO, pool, j, &m2);
            ReportConflict(tabPO,pool,m1,m2,vD,outPF);

            if(!CanSwapOutI(poolsPO,pool,i,j,tabPO,minD)) {
                printf("# No swap found pool %d [%d][%d] = %f\n",pool,m1,m2,vD);
                KillThisGuyI(j,poolsPO,pool,tabPO,outPF);
            }
        }
    }
    DB_PART DB_PrI("<< ScreenPoolI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Swap out i or j from pool into another one
*/
int CanSwapOutI(TABLE *poolsPO,int pool,int i,int j,TABLE *tabPO, DOUB minD)
{
    int q,k,iok,jok;

    DB_PART DB_PrI(">> CanSwapOutI pool=%d i=%d j=%d\n",pool,i,j);
    /***
    *   Screen i and j against members of all other pools
    */
    for(q=0;q<poolsPO->nrow;q++)
    {
        if( q==pool ) {
            continue;
        }
        /***
        *   Check if i or j can go into pool q
        */
        iok = OkInPoolI(poolsPO,pool,i,q,tabPO,minD);
        jok = OkInPoolI(poolsPO,pool,j,q,tabPO,minD);
        DB_PART DB_PrI("+ newpool q=%d iok=%d jok=%d\n",q,iok,jok);
        if( (!iok) && (!jok) ) {
            continue;
        }
        /***
        *   Check if any member k of pool q can come into pool
        */
        for(k=0;k<poolsPO->ncol;k++)
        {
            if(!OkInPoolI(poolsPO,q,k,pool,tabPO,minD)) {
                DB_PART DB_PrI("+  k=%d nope!\n",k);
                continue;
            }
            DB_PART DB_PrI("+  k=%d OK!\n",k);
            /***
            *   swap i into q and k into pool
            */
            if(iok) {
                HandlePoolSwapI(poolsPO,pool,i,q,k);
                DB_PART DB_PrI("<< CanSwapOutI i&k TRUE\n");
                return(TRUE);
            }
            else if(jok) {
                HandlePoolSwapI(poolsPO,pool,j,q,k);
                DB_PART DB_PrI("<< CanSwapOutI j&k TRUE\n");
                return(TRUE);
            }
        }
    }
    DB_PART DB_PrI("<< CanSwapOutI FALSE\n");
    return(FALSE);
}
/**************************************************************************
*   Swap index values between pools and members
*/
int HandlePoolSwapI(TABLE *poolsPO,int pool1,int m1,int pool2,int m2)
{
    DOUB v1D,v2D;

    GetTableValI(poolsPO,pool1,m1,&v1D);
    GetTableValI(poolsPO,pool2,m2,&v2D);
    SetTableValI(poolsPO,pool1,m1,v2D);
    SetTableValI(poolsPO,pool2,m2,v1D);
    printf("#  Swapped pool %d [%d] with %d [%d]\n",pool1,m1,pool2,m2);
    return(TRUE);
}
/**************************************************************************
*   Check if old pool member is ok in new pool
*/
int OkInPoolI(TABLE *poolsPO, int old, int mem, int new, TABLE *tabPO, DOUB minD)
{
    int i;
    DOUB vD;

    DB_PART DB_PrI(">> OkInPoolI old=%d mem=%d new=%d\n",old,mem,new);
    /***
    *   For each member of new pool check 
    */
    for(i=0;i<poolsPO->ncol;i++)
    {
        /***
        *   Pair-wise value old-mem new-i is forbidden?
        */
        if( GetPoolPairValueI(poolsPO,old,mem,new,i,tabPO,&vD) && (vD<=minD) ) {
            DB_PART DB_PrI("<< OkInPoolI i=%d v=%f FALSE\n",i,vD);
            return(FALSE);
        }
    }
    DB_PART DB_PrI("<< OkInPoolI TRUE\n");
    return(TRUE);
}
/**************************************************************************
*   Kill member pi from pool 
*/
int KillThisGuyI(int pi,TABLE *poolsPO,int pool,TABLE *tabPO, FILE *outPF)
{
    int m1;
    char nameS[NSIZE];
    DOUB vD;

    HAND_NFILE(outPF);
    GetTableValI(poolsPO,pool,pi,&vD);
    m1 = INT(vD);
/*
printf("pool=%d pi=%d vD=%f\n",pool,pi,vD);
*/
    BOG_CHECK(m1<0)
    BOG_CHECK(m1>=tabPO->nrow)
    GetTableRowLabI(tabPO,m1,nameS,NSIZE-1);
    fprintf(outPF,"POOL_FAIL\t%3d\t%s\n",m1,nameS);
    SetTableValI(poolsPO,pool,pi,-1.0);
    return(TRUE);
}
/**************************************************************************
*   Dump pool membership info to file
*/
int DumpPoolMembersI(TABLE *poolsPO, int pool, TABLE *tabPO, FILE *outPF)
{
    int i,m1;
    char nameS[NSIZE];
    DOUB vD;

    HAND_NFILE(outPF);
    for(i=0;i<poolsPO->ncol;i++)
    {
        GetTableValI(poolsPO,pool,i,&vD);
        m1 = INT(vD);
        if( m1<0 ) {
            continue;
        }
        GetTableRowLabI(tabPO,m1,nameS,NSIZE-1);
        fprintf(outPF,"POOL_%02d\t%3d\t%s\n",pool+1,m1,nameS);
    }
    return(FALSE);
}
/**************************************************************************
*   Report conflict between two pool elements
*/
void ReportConflict(TABLE *tabPO,int pool, int m1,int m2,DOUB vD,FILE *outPF)
{
    char name1S[NSIZE], name2S[NSIZE];
    
    HAND_NFILE(outPF);
    GetTableRowLabI(tabPO,m1,name1S,NSIZE-1);
    GetTableRowLabI(tabPO,m2,name2S,NSIZE-1);
    fprintf(outPF,"#\n");
    fprintf(outPF,"# Conflict in pool [%d] pair [%d][%d] %s %s %f\n",
        pool,m1,m2,name1S,name2S,vD);
}
/**************************************************************************
*
*/
int MaskPoolMemValsI(TABLE *poolsPO,TABLE *tabPO,TABLE *dupPO, int pool)
{
    int q,i,j,r,c;
    DOUB vD;

    InitTableValsI(dupPO,0.0,FALSE);
    /***
    *   For each pool
    */
    for(q=0;q<poolsPO->nrow;q++)
    {
        /***
        *   For each member pair i j
        */
        for(i=0;i<poolsPO->ncol;i++)
        {
            if( ! GetPoolMemTabIndexI(poolsPO, q, i, &r) ) {
                continue;
            }
            for(j=i+1;j<poolsPO->ncol;j++)
            {
                if( ! GetPoolMemTabIndexI(poolsPO, q, j, &c) ) {
                    continue;
                }
                if(pool) {
                    vD = DNUM(q+1);
                }
                else {
                    GetTableValI(tabPO,r,c,&vD);
                }
                SetTableValI(dupPO,r,c,vD);
                SetTableValI(dupPO,c,r,vD);
            }
        }
    }
    return(TRUE);
}
/**************************************************************************
*   Check members in pool all have values above min
*/
int ScanShamsI(TABLE *poolsPO,int pool, TABLE *tabPO, DOUB minD)
{
    int i,j,m1,m2;
    DOUB vD;

    /***
    *   Each member i
    */
    for(i=0;i<poolsPO->ncol;i++)
    {
        GetTableValI(poolsPO,pool,i,&vD);
        m1 = INT(vD);
        if(m1<0)
        {
            continue;
        }
        /***
        *   Each member j
        */
        for(j=i+1;j<poolsPO->ncol;j++)
        {
            GetTableValI(poolsPO,pool,j,&vD);
            m2 = INT(vD);
            if(m2<0)
            {
                continue;
            }
            GetTableValI(tabPO,m1,m2,&vD);
            if(vD <= minD)
            {
    printf("SHAM X pool %d [%d][%d] = %f\n",pool,m1,m2,vD);
            }
        }
    }
    return(TRUE);
}
/**************************************************************************
*
*/
void WritePartitionHeader(SCORETAB *stPO, TABLE *tabPO, int npools, FILE *outPF)
{
    char bufS[NSIZE],fnameS[NSIZE];

    HAND_NFILE(outPF);
    fprintf(outPF,"#\n");
    fprintf(outPF,"# Output from  %s, %s %s %s\n",VERSION_S,BD_S,__DATE__,__TIME__);
    TimeStamp("# Date         ",outPF);
    GetTableNamesI(tabPO,bufS,fnameS,NSIZE-1);
    fprintf(outPF,"# Input table: %s\n",bufS);
    fprintf(outPF,"# Input file:  %s\n",fnameS);
    FillPartitionAlgoStringI(stPO->partalg,bufS);
    fprintf(outPF,"# Algorithm:   %s\n",bufS);
    fprintf(outPF,"# Partitioning into %d pools of size %d, min = %3.4f\n",
        npools,stPO->psize,stPO->pmin);
}
/**************************************************************************
*
*/
int FillPartitionAlgoStringI(int alg, char *nameS)
{
    int ok;

    ok = FALSE;
    switch(alg)
    {
        case PALG_BTS:  
            sprintf(nameS,"Block-then-swap"); 
            ok++;   
            break;
        default:
            sprintf(nameS,"????"); 
            break;
    }
    return(ok);
}
