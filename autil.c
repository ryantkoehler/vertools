/*
* autil.c
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
#include <math.h>
#include "prim.h"

#define DB_ARRAY    if(DB[15])
#define DB_RMASK    if(DB[16])

/******************************************************************** aaa
*   Initializes passed array with passed value 
*   For loop from start < end
*/
int InitArrayI(void *aPO, int vt, int start, int end, DOUB valD)
{
    DOUB *aPD;
    int *aPI, i;
    short int *fPH;
    char *aPC;

    DB_ARRAY DB_PrI(">> InitArray %d %d vt=%d",start,end,vt);
    if( (start<0) || (end<=start) )
    {
        return(BOGUS);
    }
    switch(vt)
    {
        case IS_DOUB:
            aPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   aPD[i] = valD;  }
            break;
        case IS_CHAR:
            aPC = (char *)aPO;
            for(i=start; i<end; i++)
            {   aPC[i] = CHAR(valD);    }
            break;
        case IS_INT:
            aPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {   aPI[i] = INT(valD); }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {   fPH[i] = SHORT(valD);   }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("InitArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< InitArray TRUE\n");
    return(TRUE);
}
/******************************************************************** aaa
*   Scales the passed array with passed value (multiplication)
*/
int ScaleArrayI(void *aPO,int vt, int start, int end, DOUB valD)
{
    DOUB *aPD;
    int i, *aPI;
    short int *fPH;

    DB_ARRAY DB_PrI(">> ScaleArray vR=%f vt=%d\n",valD,vt);
    if( (start<0) || (end<=start) )
    {
        return(FALSE);
    }
    switch(vt)
    {
        case IS_DOUB:
            aPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   aPD[i] *= valD; }
            break;
        case IS_INT:
            aPI = (int *)aPO;
            for(i=start; i<end; i++)
            {   aPI[i] = INT(DNUM(aPI[i]) * valD);  }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {   fPH[i] = SHORT(DNUM(fPH[i]) * valD);    }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("ScaleArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< ScaleArray TRUE\n");
    return(TRUE);
}
/******************************************************************** aaa
*   Shifts the passed array with passed value (addition subtraction)
*/
int ShiftArrayI(void *aPO,int vt, int start, int end, DOUB valD)
{
    DOUB *aPD;
    int i, *aPI;
    short int *fPH;

    DB_ARRAY DB_PrI(">> ShiftArray %d %d vt=%d\n",start,end,vt);
    if( (start<0) || (end<=start) )
    {
        DB_ARRAY DB_PrI("<< ShiftArray bad start-end FALSE\n");
        return(FALSE);
    }
    switch(vt)
    {
        case IS_DOUB:
            aPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   aPD[i] += valD; }
            break;
        case IS_INT:
            aPI = (int *)aPO;
            for(i=start; i<end; i++)
            {   aPI[i] = INT( DNUM(aPI[i]) + valD); }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {   fPH[i] = SHORT(DNUM(fPH[i]) + valD);    }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("ShiftArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< ShiftArray TRUE\n");
    return(TRUE);
}
/******************************************************************** aaa
*   Normalizes the passed array 
*       range 0 to 1 for REAL, DOUB
*       range 0 to 100 for int
*/
int NormArrayI(void *aPO,int vt, int start, int end)
{
    DOUB *aPD, hiD,loD,delD;
    int *aPI, i;
    short int *fPH;

    DB_ARRAY DB_PrI(">> NormArray vt=%d\n",vt);
    if(!ArrayStatsI(aPO, vt, start, end, &loD, &hiD, NULL, NULL)) {
        DB_ARRAY DB_PrI("<< NormArray FALSE\n");
        return(FALSE);
    }
    delD = hiD - loD;
    LIMIT_NUM(delD,TINY_R,TOO_BIG_R);
    switch(vt)
    {   
        case IS_DOUB:
            aPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   
                aPD[i] = (aPD[i] - loD) / delD; 
            }
            break;
        case IS_INT:
            aPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {
                aPI[i] = ROUND( (100.0 * RNUM(aPI[i]) - loD) / delD);
            }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {
                fPH[i] = ROUND( (100.0 * RNUM(fPH[i]) - loD) / delD);
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("NormArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< NormArrayI TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Bound array values 
*/
int BoundArrayI(void *aPO, int vt, int start, int end, DOUB lowD, DOUB hiD)
{
    DOUB *aPD;
    int *aPI, i, hi, low;
    short int *aPH;

    DB_ARRAY DB_PrI(">> BoundArray vt=%d\n",vt);
    if( (start<0) || (end<=start) ) {
        DB_ARRAY DB_PrI("<< BoundArray bad start-end FALSE\n");
        return(FALSE);
    }
    switch(vt)
    {
        case IS_DOUB:
            aPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   
                LIMIT_NUM(aPD[i],lowD,hiD); 
            }
            break;
        case IS_INT:
            aPI = (int *)aPO;   
            hi = INT(hiD);  
            low = INT(lowD);    
            for(i=start; i<end; i++)
            {   
                LIMIT_NUM(aPI[i],low,hi); 
            }
            break;
        case IS_SHORT:
            aPH = (short int *)aPO; 
            hi = SHORT(hiD);    
            low = SHORT(lowD);  
            for(i=start; i<end; i++)
            {   
                LIMIT_NUM(aPH[i],low,hi); 
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("BoundArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< BoundArray TRUE\n");
    return(TRUE);
}
/*************************************************************************
*   Count values in range low to hi
*/
int NumArrayValsI(void *aPO, int vt, int start, int end, DOUB lowD, DOUB hiD)
{
    DOUB *aPD;
    int *aPI, i, n;
    short int *aPH;
    char *cPC;

    DB_ARRAY DB_PrI(">> NumArrayValsI vt=%d %f %f\n",vt,lowD,hiD);
    if( (start<0) || (end<=start) ) {
        DB_ARRAY DB_PrI("<< BoundArray bad start-end BOGUS\n");
        return(BOGUS);
    }
    n = 0;
    switch(vt)
    {
        case IS_DOUB:
            aPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   
                if( (DNUM(aPD[i])>=lowD) && (DNUM(aPD[i])<=hiD) )
                {   n++;    }
            }
            break;
        case IS_INT:
            aPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {   
                if( (DNUM(aPI[i])>=lowD) && (DNUM(aPI[i])<=hiD) )
                {   n++;    }
            }
            break;
        case IS_SHORT:
            aPH = (short int *)aPO;     
            for(i=start; i<end; i++)
            {   
                if( (DNUM(aPH[i])>=lowD) && (DNUM(aPH[i])<=hiD) )
                {   n++;    }
            }
            break;
        case IS_CHAR:
            cPC = (char *)aPO;  
            for(i=start; i<end; i++)
            {   
                if( (DNUM(cPC[i])>=lowD) && (DNUM(cPC[i])<=hiD) )
                {   n++;    }
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("NumArrayValsI","Bogus code");
    }
    DB_ARRAY DB_PrI("<< NumArrayValsI %d\n",n);
    return(n);
}
/*************************************************************************
*   Copy one array into another
*/
void CopyArray(void *fPO, void *sPO, int vt, int start, int end)
{
    DOUB *fPD, *sPD;
    int i, j, *fPI, *sPI;

    DB_ARRAY DB_PrI(">> CopyArray vt=%d %d to %d\n",vt,start,end);
    j = 0;
    switch(vt)
    {
        case IS_DOUB:
            fPD = (DOUB *)fPO;  
            sPD = (DOUB *)sPO;  
            for(i=start; i<end; i++)
            {   
                sPD[j++] = fPD[i];  
            }
            break;
        case IS_INT:
            fPI = (int *)fPO;   
            sPI = (int *)sPO;   
            for(i=start; i<end; i++)
            {   
                sPI[j++] = fPI[i];  
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("CopyArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< CopyArray\n");
}
/*************************************************************************
*   Smooth array 
*/
void SmoothArray(void *fPO, void *sPO, int vt, int n, int win)
{
    DOUB *fPD, *sPD, avD;
    int i, j, s, e, w;

    DB_ARRAY DB_PrI(">> SmoothArray vt=%d win=%d\n",vt,win);
    LIMIT_NUM(win,1,n);
    switch(vt)
    {
        case IS_DOUB:
            fPD = (DOUB *)fPO;  
            sPD = (DOUB *)sPO;  
            for(i=0; i<n; i++)
            {   
                s = i - win;
                e = i + win;
                LIMIT_NUM(s,0,n);
                LIMIT_NUM(e,0,n);
                w = 0;
                avD = 0.0;
                for(j=s; j<e; j++) {
                    avD += fPD[j];
                    w++;
                }
                avD = avD / RNUM(w);
                sPD[i] = avD;   
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("SmoothArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< SmoothArray\n");
}
/*************************************************************************
*   Mix first array (fPO) with second (sPO) into target (tPO)
*   All should be n long
*   Combination operator is mix
*   NO CHECK FOR /= 0
*/
void MixArrays(void *fPO, void*sPO, void *tPO, int vt, int n, int mix)
{
    DOUB *fPD, *sPD, *tPD;
    int i, *fPI, *sPI, *tPI;

    DB_ARRAY DB_PrI(">> MixArrays vt=%d\n",vt);
    fPD = sPD = tPD = NULL;
    fPI = sPI = tPI = NULL;
    switch(vt)
    {
        case IS_INT:
            fPI = (int *)fPO;   sPI = (int *)sPO;   tPI = (int *)tPO;   
            for(i=0; i<n; i++)
            {
                switch(mix)
                {
                    case MATH_ADD:  tPI[i] = fPI[i] + sPI[i];       break;
                    case MATH_SUB:  tPI[i] = fPI[i] - sPI[i];       break;
                    case MATH_MUL:  tPI[i] = fPI[i] * sPI[i];       break;
                    case MATH_DIV:  tPI[i] = fPI[i] / sPI[i];       break;
                    default:    
                        DB_PrI("BAD MIX CODE: %d\n",mix);
                        ERR("MixArrays","Unknown mixing operation code");
                        return;
                }
            }
            break;
        case IS_DOUB:
            fPD = (DOUB *)fPO;  sPD = (DOUB *)sPO;  tPD = (DOUB *)tPO;  
            for(i=0; i<n; i++)
            {
                switch(mix)
                {
                    case MATH_ADD:  tPD[i] = fPD[i] + sPD[i];       break;
                    case MATH_SUB:  tPD[i] = fPD[i] - sPD[i];       break;
                    case MATH_MUL:  tPD[i] = fPD[i] * sPD[i];       break;
                    case MATH_DIV:  tPD[i] = fPD[i] / sPD[i];       break;
                    default:    
                        DB_PrI("BAD MIX CODE: %d\n",mix);
                        ERR("MixArrays","Unknown mixing operation code");
                        return;
                }
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("MixArrays","Bogus code");
            return;
    }
    DB_ARRAY DB_PrI("<< MixArrays");
}
/*************************************************************************
*   Dumps an int or real array using the format passed
*/
void DumpArray(void *aPO, int vt, int start, int end, char *formS, FILE *oPF)
{
    int i, *fPI;
    DOUB *fPD;
    char *aPC;
    char pformS[DEF_BS];

    HAND_NFILE(oPF);
    DB_ARRAY DB_PrI(">> DumpArray vt=%d\n",vt);
    if(formS) {
        strcpy(pformS,formS);
    }
    switch(vt) {
        case IS_DOUB:   
            if(!formS) {
                sprintf(pformS,"[%%d]\t%%4.2f\n");  break;
            }
            break;
        case IS_INT:    
            if(!formS) {
                sprintf(pformS,"[%%d]\t%%4d\n");    break;
            }
            break;
        case IS_CHAR:   
            if(!formS) {
                sprintf(pformS,"[%%d]\t%%2d\n");    break;
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("DumpArray","Bogus code");
    }
    for(i=start; i<end; i++)
    {
        switch(vt)
        {
            case IS_DOUB:
                fPD = (DOUB *)aPO;  
                fprintf(oPF,pformS,i,fPD[i]);   
                break;
            case IS_INT:
                fPI = (int *)aPO;   
                fprintf(oPF,pformS,i,fPI[i]);   
                break;
            case IS_CHAR:
                aPC = (char *)aPO;
                fprintf(oPF,pformS,aPC[i]); 
                break;
        }
    }
    DB_ARRAY DB_PrI("<< DumpArray\n");
} 
/*************************************************************************
*   Return summed value of all array elements
*/
void ArraySum(void *aPO, int vt, int start, int end, DOUB *sumPD)
{
    int i, *fPI;
    DOUB *fPD,sumD;
    char *aPC;

    sumD = 0.0;
    switch(vt)
    {
        case IS_DOUB:
            fPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {   
                sumD += fPD[i]; 
            }
            break;
        case IS_INT:
            fPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {   
                sumD += DNUM(fPI[i]);   
            }
            break;
        case IS_CHAR:
            aPC = (char *)aPO;
            for(i=start; i<end; i++)
            {   
                sumD += DNUM(aPC[i]);   
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("SumArray","Bogus code");
    }
    *sumPD = sumD;
}
/*************************************************************************
*   calculates stats (low, high, average) 
*   SD is only calculated if a real variable address is passed
*/
int ArrayStatsI(void *aPO, int vt, int start, int end, DOUB *loPD, DOUB *hiPD, 
    DOUB *avPD, DOUB *sdPD)
{
    int *fPI,i;
    short int *fPH;
    DOUB loD,hiD,sumD,avD,*fPD;
    char *aPC;
    
    DB_ARRAY DB_PrI(">> ArrayStats, %d %d, ct=%d\n",start,end,vt);
    if( (start<0) || (end<=start) )
    {
        return(FALSE);
    }
    loD = TOO_BIG_R; hiD = -TOO_BIG_R;
    sumD = 0.0;
    switch(vt)
    {
        case IS_DOUB:
            fPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {
                sumD += DNUM(fPD[i]);
                loD = MIN_NUM(fPD[i],loD);
                hiD = MAX_NUM(fPD[i],hiD);
            }
            break;
        case IS_INT:
            fPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {
                sumD += DNUM(fPI[i]);
                loD = MIN_NUM(DNUM(fPI[i]),loD);
                hiD = MAX_NUM(DNUM(fPI[i]),hiD);
            }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {
                sumD += DNUM(fPH[i]);
                loD = MIN_NUM(DNUM(fPH[i]),loD);
                hiD = MAX_NUM(DNUM(fPH[i]),hiD);
            }
            break;
        case IS_CHAR:
            aPC = (char *)aPO;
            for(i=start; i<end; i++)
            {
                sumD += DNUM(aPC[i]);
                loD = MIN_NUM(DNUM(aPC[i]),loD);
                hiD = MAX_NUM(DNUM(aPC[i]),hiD);
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("ArrayStats","Bogus code");

    }
    avD = sumD / DNUM(end - start);
    /***
    *   Set low, hi, average numbers
    */
    if(loPD) {
        *loPD = loD;
    }
    if(hiPD) {
        *hiPD = hiD;
    }
    if(avPD) {
        *avPD = avD;
    }
    DB_ARRAY DB_PrI("+ ArrayStatsI lo=%f hi=%f sum=%f av=%f\n",loD,hiD,sumD,avD);
    /***
    *   If not doing standard deviation, we're done
    */
    if(!sdPD) {
        DB_ARRAY DB_PrI("<< ArrayStatsI no SD so done TRUE\n");
        return(TRUE);
    }
    /***
    *   Pass two to get sd
    */
    sumD = 0.0;
    switch(vt)
    {
        case IS_DOUB:
            fPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {
                hiD = fPD[i] - avD;
                sumD += (hiD * hiD);
            }
            break;
        case IS_INT:
            fPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {
                hiD = DNUM(fPI[i]) - avD;
                sumD += (hiD * hiD);
            }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {
                hiD = DNUM(fPH[i]) - avD;
                sumD += (hiD * hiD);
            }
            break;
        case IS_CHAR:
            aPC = (char *)aPO;
            for(i=start; i<end; i++)
            {
                hiD = DNUM(aPC[i]) - avD;
                sumD += (hiD * hiD);
            }
            break;
    }
    sumD /= DNUM(end-start);
    *sdPD = sqrt(sumD);
    return(TRUE);
}
/************************************************************************
*   Allocates and fills a histogram count array for the n elements of the
*       passed array aPO.  
*   aPO = array pointer
*   vt = array variable type code
*   start = array start
*   end = array end
*   bsD = Bucket size 
*   lvD = Low value, i.e. where hist starts
*   hvD = High value, i.e. where hist ends
*   ncPI = number of counted elements; value set 
*   hPPI = allocated histogram of counts
*
*   Returns number of bins in histogram (i.e. array dim)
*/
int ArrayHistI(void *aPO, int vt, int start, int end, DOUB bsD, 
    DOUB lvD, DOUB hvD, int *ncPI, DOUB **hPPD)
{
    int *fPI,i,buck,nb,nc;
    short int *fPH;
    DOUB *fPD,*hisPD,hiD,loD;
    char *aPC;

    DB_ARRAY DB_PrI(">> ArrayHistI %d %d bs=%f lv=%f hv=%f\n",
        start,end,bsD,lvD,hvD);
    *hPPD = NULL;
    if( (start<0) || (end<=start) ) {
        DB_ARRAY DB_PrI("<< ArrayHistI start-end bad FALSE\n");
        return(FALSE);
    }
    /***
    *   Count range of values
    */
    if(!ArrayStatsI(aPO, vt, start, end, &loD, &hiD, NULL, NULL)) {
        DB_ARRAY DB_PrI("<< ArrayHistI stats bad FALSE\n");
        return(FALSE);
    }
    DB_ARRAY DB_PrI("+ lo=%f hi=%f   start=%d end=%d\n",loD,hiD,start,end);
    /***
    *   How many buckets? Allocate for them
    */
    nb = INT( (hiD-lvD) / bsD ) + 1;
    if(nb <= 0) {   
        DB_ARRAY DB_PrI("<< ArrayHistI no bins FALSE\n");
        return(FALSE);  
    }
    if(!(hisPD = (DOUB *)ALLOC(nb,sizeof(DOUB)))) { 
        DB_ARRAY DB_PrI("<< ArrayHistI Alloc failed FALSE\n");
        return(FALSE);  
    }
    DB_ARRAY DB_PrI("+ ArrayHistI allocated %d buckets\n",nb);
    /***
    *   Pass two, tally values
    */
    nc = 0;
    switch(vt)
    {
        case IS_DOUB:
            fPD = (DOUB *)aPO;  
            for(i=start; i<end; i++)
            {
                if(fPD[i] < lvD) {
                    continue;
                }
                buck = INT( (fPD[i] - lvD) / bsD );
                hisPD[buck] += 1;
                nc++;
            }
            break;
        case IS_INT:
            fPI = (int *)aPO;   
            for(i=start; i<end; i++)
            {
                if(fPI[i] < INT(lvD)) {
                    continue;
                }
                buck = INT( (DNUM(fPI[i]) - lvD) / bsD );
                hisPD[buck] += 1;
                nc++;
            }
            break;
        case IS_SHORT:
            fPH = (short int *)aPO; 
            for(i=start; i<end; i++)
            {
                if(fPH[i] < INT(lvD)) {
                    continue;
                }
                buck = INT( (DNUM(fPH[i]) - lvD) / bsD );
                hisPD[buck] += 1;
                nc++;
            }
            break;
        case IS_CHAR:
            aPC = (char *)aPO;
            for(i=start; i<end; i++)
            {
                if(aPC[i] < INT(lvD)) {
                    continue;
                }
                buck = INT( (DNUM(aPC[i]) - lvD) / bsD );
                hisPD[buck] += 1;
                nc++;
            }
            break;
    }
    DB_ARRAY DB_PrI("+ ArrayHistI whip\n");
    *hPPD = hisPD;
    if(ncPI != NULL) {  
        *ncPI = nc; 
    }
    DB_ARRAY DB_PrI("<< ArrayHistI %d\n",nb);
    return(nb);
}
/****************************************************************************
*   Fill passed array, seq[len], with random ints 0 to max
*   If unique is true, no two random numbers should be the same
*/
int ArrayRandSequenceI(int *seq, int len, int max, int unique)
{
    int i,u,d,ntake;

    if( (unique) && (max<len) ) {
        printf("Bad unique sequence parameters: unique=%d len=%d, max=%d\n",unique,len,max);
        ERR("ArrayRandSequenceI","impossible parameters");
    }
    ntake = 0;
    while(ntake < len)
    {
        d = RandI(max);
        if(unique) {
            u = TRUE;
            for(i=0;i<ntake;i++)
            {
                if(d==seq[i]) {
                    u = FALSE;  break;
                }
            }
            if(u==TRUE) {
                seq[ntake++] = d;
            }
        }
        else {
            seq[ntake++] = d;
        }
    }
    return(ntake);
}
/****************************************************************************/
void InvertMask(char *maskPC,int len)
{
    int i;

    for(i=0;i<len;i++)
    {
        if(maskPC[i]) {
            maskPC[i] = FALSE;
        }
        else {
            maskPC[i] = TRUE;
        }
    }
}
/****************************************************************************/
int MaskCountI(char *maskPC,int len)
{
    int i,n;

    n = 0;
    for(i=0;i<len;i++)
    {
        if(maskPC[i]) {
            n++;
        }
    }
    return(n);
}
/**************************************************************************
*   Fill maskPC[len] with a random subset of fracR fraction of 1's
*/
int MaskRandSubsetI(char *maskPC, int len, DOUB fracR)
{
    int i,n,rval,need;
    char v;

    DB_RMASK DB_PrI(">> MaskRandSubsetI len %d, frac %f\n",len,fracR);
    if(len<=0)
    {
        DB_RMASK DB_PrI("<< MaskRandSubsetI len<=0 0\n");
        return(0);
    }
    /***
    *   initialize all unset or all set
    */
    LIMIT_NUM(fracR,0.0,1.0);
    if(fracR < 0.5)
    {
        v = 0;
        need = ROUND(RNUM(len) * fracR);
        DB_RMASK DB_PrI("+ v = %d, need = %d\n",v,need);
    }
    else
    {
        v = 1;
        need = ROUND(RNUM(len) * (1.0 - fracR));
        DB_RMASK DB_PrI("+ v = %d, need = %d\n",v,need);
    }
    for(i=0;i<len;i++)
    {   
        maskPC[i] = v;  
    }
    /***
    *   Set until all set
    */
    DB_RMASK DB_PrI("+ mask initialized, setting %d\n",need);
    n = 0;
    while(n < need)
    {
        for(i=0;i<len;i++)
        {
            if(maskPC[i] != v)
                continue;
            rval = RandI(len);
            DB_RMASK DB_PrI("+  %d vs %d ",rval,need);
            if(rval < need)
            {
                maskPC[i] = !v;
                n++;
                DB_RMASK DB_PrI("n = %d\n",n);
            }
            else DB_RMASK DB_PrI("nope\n");
            if(n == need)
                break;
        }
    }
    if(v == 1)
    { n = len - n; }
    DB_RMASK DB_PrI("<< MaskRandSubsetI %d\n",n);
    return(n);
}
/******************************************************************** 
*   Qsort functions
*/
PRIV_I qsort_doub(const void *e1, const void *e2);
PRIV_I qsort_doub_decend(const void *e1, const void *e2);
PRIV_I qsort_real(const void *e1, const void *e2);
PRIV_I qsort_real_decend(const void *e1, const void *e2);
PRIV_I qsort_int(const void *e1, const void *e2);
PRIV_I qsort_int_decend(const void *e1, const void *e2);
PRIV_I qsort_short(const void *e1, const void *e2);
PRIV_I qsort_short_decend(const void *e1, const void *e2);
PRIV_I qsort_char(const void *e1, const void *e2);
PRIV_I qsort_doub_decend(const void *e1, const void *e2);
/**********************************************************/
PRIV_I qsort_doub(const void *e1, const void *e2)
{
    DOUB fD,sD;
    
    fD = *(DOUB *)e1;   sD = *(DOUB *)e2;
    if(fD>sD)   {return(1);}
    if(fD<sD)   {return(-1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_doub_decend(const void *e1, const void *e2)
{
    DOUB fD,sD;
    
    fD = *(DOUB *)e1;   sD = *(DOUB *)e2;
    if(fD>sD)   {return(-1);}
    if(fD<sD)   {return(1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_int(const void *e1, const void *e2)
{
    int f,s;
    
    f = *(int *)e1; s = *(int *)e2;
    if(f>s) {return(1);}
    if(f<s) {return(-1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_int_decend(const void *e1, const void *e2)
{
    int f,s;
    
    f = *(int *)e1; s = *(int *)e2;
    if(f>s) {return(-1);}
    if(f<s) {return(1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_short(const void *e1, const void *e2)
{
    int f,s;
    
    f = *(short int *)e1;   s = *(short int *)e2;
    if(f>s) {return(1);}
    if(f<s) {return(-1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_short_decend(const void *e1, const void *e2)
{
    int f,s;
    
    f = *(short int *)e1;   s = *(short int *)e2;
    if(f>s) {return(-1);}
    if(f<s) {return(1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_char(const void *e1, const void *e2)
{
    char fC,sC;
    
    fC = *(char *)e1;   sC = *(char *)e2;
    if(fC>sC)   {return(1);}
    if(fC<sC)   {return(-1);}
    return(0);
}
/**********************************************************/
PRIV_I qsort_char_decend(const void *e1, const void *e2)
{
    char fC,sC;
    
    fC = *(char *)e1;   sC = *(char *)e2;
    if(fC>sC)   {return(-1);}
    if(fC<sC)   {return(1);}
    return(0);
}
/******************************************************************** 
*   Sort passed array 
*/
void SortArray(void *aPO, int vt, int n, int dir)
{
    DB_ARRAY DB_PrI(">> SortArray %d dir=%d vt=%d",n,dir,vt);
    switch(vt)
    {
        case IS_DOUB:
            if(dir<0) {
                qsort(aPO,n,sizeof(DOUB),qsort_doub_decend);
            }
            else {
                qsort(aPO,n,sizeof(DOUB),qsort_doub);
            }
            break;
        case IS_INT:
            if(dir<0) {
                qsort(aPO,n,sizeof(int),qsort_int_decend);
            }
            else {
                qsort(aPO,n,sizeof(int),qsort_int);
            }
            break;
        case IS_SHORT:
            if(dir<0) {
                qsort(aPO,n,sizeof(short int),qsort_short_decend);
            }
            else {
                qsort(aPO,n,sizeof(short int),qsort_short);
            }
            break;
        case IS_CHAR:
            if(dir<0) {
                qsort(aPO,n,sizeof(char),qsort_char_decend);
            }
            else {
                qsort(aPO,n,sizeof(char),qsort_char);
            }
            break;
        default:
            printf("Bogus value-type code = %d\n",vt);
            ERR("SortArray","Bogus code");
    }
    DB_ARRAY DB_PrI("<< SortArray\n");
}
