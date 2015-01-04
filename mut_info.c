/*
* mut_info.c
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
#include <math.h>
#include "prim.h"
#include "numlist.h"
#include "stat.h"
#include "table.h"
#include "histogram.h"
#include "mut_info.h"

#define DB_MINF     if(DB[38])

/**************************************************************************/
int NumlistMutInfoI(NUMLIST *nums1PO, NUMLIST *nums2PO, int maxb, DOUB *minfPD, 
    int *nb1PI, int *nb2PI)
{
    int i,r,c,n,nrow,ncol,dis_bin;
    DOUB vD, bin1D, bin2D, st1D, st2D, en1D, en2D, infoD;
    HISTOGRAM *his1PO, *his2PO;
    TABLE *tabPO;

    DB_MINF DB_PrI(">> NumlistMutInfoI maxb=%d\n",maxb);
    VALIDATE(nums1PO,NUMLIST_ID);
    VALIDATE(nums2PO,NUMLIST_ID);
    if( !NumlistSameLenI(nums1PO,nums2PO) ) {
        DB_MINF DB_PrI("<< NumlistMutInfoI not same len FALSE\n");
        return(FALSE);
    }
    n = GetNumlistLengthI(nums1PO);
    DB_MINF DB_PrI("+ number lists %d long\n",n);
    /***
    *   Get histograms from number lists
    */
    dis_bin = FALSE;
    NumlistNaturalHistBinI(nums1PO, maxb, dis_bin, &bin1D, &st1D, &en1D);
    NumlistNaturalHistBinI(nums2PO, maxb, dis_bin, &bin2D, &st2D, &en2D);
    DB_MINF DB_PrI("+ Hist1 bin=%f st=%f en=%f\n",bin1D,st1D,en1D);
    DB_MINF DB_PrI("+ Hist2 bin=%f st=%f en=%f\n",bin2D,st2D,en2D);
    his1PO = his2PO = NULL;
    NumlistToHistogramI(nums1PO, bin1D, st1D, en1D, &his1PO);
    NumlistToHistogramI(nums2PO, bin2D, st2D, en2D, &his2PO);
    if( (!his1PO) || (!his2PO) ) {
        CHECK_HISTOGRAM(his1PO);
        CHECK_HISTOGRAM(his2PO);
        DB_MINF DB_PrI("<< NumlistMutInfoI nums to histograms %p %p FALSE\n",his1PO,his2PO);
        return(FALSE);
    }
    nrow = GetHistogramNumBinsI(his1PO);
    ncol = GetHistogramNumBinsI(his2PO);
    DB_MINF DB_PrI("+ histogram bins nrow=%d ncol=%d\n",nrow,ncol);
    if( (nrow < 2) || (ncol < 2) ) {
        CHECK_HISTOGRAM(his1PO);
        CHECK_HISTOGRAM(his2PO);
        DB_MINF DB_PrI("<< NumlistMutInfoI not enough bins FALSE\n");
        return(FALSE);
    }
    /***
    *   Allocate space
    */
    tabPO = CreateTablePO(nrow, ncol);
    if( !tabPO ) {
        CHECK_HISTOGRAM(his1PO);
        CHECK_HISTOGRAM(his2PO);
        DB_MINF DB_PrI("<< NumlistMutInfoI failed table space FALSE\n");
        return(FALSE);
    }
    /***
    *   Fill joint probabilities from pairwise screen of input data
    */
    for(i=0;i<n;i++)
    {
        GetNumlistIntDoubI(nums1PO,i,NULL,&vD);
        HistogramBinForValueI(his1PO, vD, &r, NULL,NULL);
        GetNumlistIntDoubI(nums2PO,i,NULL,&vD);
        HistogramBinForValueI(his2PO, vD, &c, NULL,NULL);
        GetTableValI(tabPO,r,c,&vD);
        vD += 1.0;
        SetTableValI(tabPO,r,c,vD);
    }

    DB_MINF {
        printf("\n xxxxxxxxxxxxxxxxxxxxxx\n");
        DumpNumlist(his1PO->bins,-1,-1,NULL);
        DumpNumlist(his2PO->bins,-1,-1,NULL);
        AutoTableOutFormattingI(tabPO,TRUE,TRUE);
        DumpTable(tabPO,TRUE,FALSE,NULL);
    }

    /***
    *   Now the info
    */
    infoD = 0.0;
    for(r=0;r<nrow;r++)
    {
        GetNumlistIntDoubI(his1PO->bins,r,NULL,&bin1D);
        bin1D = bin1D / DNUM(n);
        for(c=0;c<ncol;c++)
        {
            GetTableValI(tabPO,r,c,&vD);
            if(vD == 0.0) {
                continue;
            }
            vD = vD / DNUM(n);

            GetNumlistIntDoubI(his2PO->bins,c,NULL,&bin2D);
            bin2D = bin2D / DNUM(n);

/*
printf("[%d][%d]\tbin1=%f\tbin2=%f\tv=%f",r,c,bin1D,bin2D,vD);
*/
            vD = vD * LOG_2(vD / (bin1D * bin2D));
            infoD += vD;
/*
printf("\tV=%f\tinfo=%f\n",vD,infoD);
*/

            SetTableValI(tabPO,r,c,vD);
        }
    }

    DB_MINF {
        AutoTableOutFormattingI(tabPO,TRUE,TRUE);
        DumpTable(tabPO,TRUE,FALSE,NULL);
    }

    DB_MINF DB_PrI("+ info = %f\n",infoD);
    /***
    *   Set values and clean up
    */
    if(nb1PI) {
        *nb1PI = nrow;
    }
    if(nb2PI) {
        *nb2PI = ncol;
    }
    if(minfPD) {
        *minfPD = infoD;
    }
    CHECK_TABLE(tabPO);
    CHECK_HISTOGRAM(his1PO);
    CHECK_HISTOGRAM(his2PO);
    DB_MINF DB_PrI("<< NumlistMutInfoI TRUE\n");
    return(TRUE);
}
