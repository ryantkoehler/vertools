/*
* dna_cons.c
*
* Copyright 2019 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <math.h>
#include "prim.h"
#include "dna.h"
#include "tm_pars.h"
#include "dna_cons.h"

#define DB_CONS if(DB[120])

/*****************************************************************************
*   Allocate IN_CONS data structure
*/
IN_CONS *CreateInConsPO()
{
    IN_CONS *iconPO;

    if(!(iconPO = (IN_CONS *)ALLOC(1,sizeof(IN_CONS))))
    {   return(NULL);   }
    iconPO->ID = IN_CONS_ID;
    InitInCons(iconPO);
    return(iconPO);
}
/*****************************************************************************
*   Free IN_CONS data structure 
*/
int DestroyInConsI(IN_CONS *iconPO)
{
    VALIDATE(iconPO,IN_CONS_ID);
    FREE(iconPO);
    return(TRUE);
}
/*****************************************************************************
*   Initialize IN_CONS data structure vars
*/
void InitInCons(IN_CONS *iconPO)
{
    VALIDATE(iconPO,IN_CONS_ID);
    INIT_S(iconPO->parname);
    /***
    *   Default intrinsic constraints = no limits or bogus values 
    */
    iconPO->maxlen = 0;
    iconPO->minAc = iconPO->minCc = iconPO->minGc = iconPO->minTc = 0.0;
    iconPO->maxAc = iconPO->maxCc = iconPO->maxGc = iconPO->maxTc = 1.0;
    iconPO->minSc = iconPO->minRc = iconPO->minKc = 0.0;
    iconPO->minWc = iconPO->minYc = iconPO->minMc = 0.0;
    iconPO->maxSc = iconPO->maxRc = iconPO->maxKc = 1.0;
    iconPO->maxWc = iconPO->maxYc = iconPO->maxMc = 1.0;
    iconPO->max_rowA = iconPO->max_rowC = BOGUS;
    iconPO->max_rowG = iconPO->max_rowT = BOGUS;
    iconPO->max_rowS = iconPO->max_rowW = BOGUS;
    iconPO->max_rowY = iconPO->max_rowR = BOGUS;
    iconPO->max_rowM = iconPO->max_rowK = BOGUS;
    iconPO->max_repRY = iconPO->max_repYR = BOGUS;
    iconPO->max_repMK = iconPO->max_repKM = BOGUS;
}
/*****************************************************************************
*   Load parameters 
*/
int LoadInConsParsI(char *fnameS, IN_CONS *pxPO)
{
    FILE *fPF;
    int ok;

    if(!(fPF=OpenUFilePF(fnameS,"r",NULL)))
    {
        return(FALSE);
    }
    ok = ParseInConsParsI(fPF,pxPO);
    if(ok)
    {
        strcpy(pxPO->parname,fnameS);
    }
    FILECLOSE(fPF);
    return(ok);
}

/*****************************************************************************
*   Parse gen_seq parameters from keyword-value file
*/
int ParseInConsParsI(FILE *fPF,IN_CONS *iconPO)
{
    int l,com;
    char bufS[DEF_BS];

    VALIDATE(iconPO,IN_CONS_ID);
    com=FALSE;
    l = 0;
    while(fgets(bufS,LINEGRAB,fPF) != NULL)
    {
        l++;
        if(EQSTRING(bufS,"break",5))
            break;
        if(COM_LINE(bufS))
            continue;
        if(BlankStringI(bufS))
            continue;
        if(strstr(bufS,"/*"))
        {   
            com=TRUE;
            continue;   
        }
        if(strstr(bufS,"*/"))
        {
            com=FALSE;
            continue;
        }
        if(com)
        {   continue;   }
        Upperize(bufS);
        /***
        *   A content 
        */
        if(EQSTRING(bufS,"MIN_AC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minAc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_AC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxAc); 
            continue;
        }
        /***
        *   C content 
        */
        if(EQSTRING(bufS,"MIN_CC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minCc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_CC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxCc); 
            continue;
        }
        /***
        *   G content 
        */
        if(EQSTRING(bufS,"MIN_GC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minGc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_GC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxGc); 
            continue;
        }
        /***
        *   T content 
        */
        if(EQSTRING(bufS,"MIN_TC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minTc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_TC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxTc); 
            continue;
        }
        /***
        *   S & W content 
        */
        if(EQSTRING(bufS,"MIN_SC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minSc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_SC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxSc); 
            continue;
        }
        if(EQSTRING(bufS,"MIN_WC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minWc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_WC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxWc); 
            continue;
        }
        /***
        *   R & Y content 
        */
        if(EQSTRING(bufS,"MIN_RC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minRc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_RC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxRc); 
            continue;
        }
        if(EQSTRING(bufS,"MIN_YC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minYc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_YC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxYc); 
            continue;
        }
        /***
        *   M & K content 
        */
        if(EQSTRING(bufS,"MIN_MC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minMc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_MC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxMc); 
            continue;
        }
        if(EQSTRING(bufS,"MIN_KC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->minKc); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_KC",6))
        { 
            sscanf(bufS,"%*s %lf",&iconPO->maxKc); 
            continue;
        }
        /***
        *   Max bases in a row
        */
        if(EQSTRING(bufS,"MAX_AROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowA); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_CROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowC); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_GROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowG); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_TROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowT); 
            continue;
        }
        /***
        *   Max dengerate bases in a row
        */
        if(EQSTRING(bufS,"MAX_SROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowS); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_WROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowW); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_YROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowY); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_RROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowR); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_MROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowM); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_KROW",8))
        {   
            sscanf(bufS,"%*s %d",&iconPO->max_rowK); 
            continue;
        }
        /*** SHAM ?
        *   Max pYridine / puRine repeats
        */
        if(EQSTRING(bufS,"MAX_RYREP",9))
        { 
            sscanf(bufS,"%*s %d",&iconPO->max_repRY); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_YRREP",9))
        { 
            sscanf(bufS,"%*s %d",&iconPO->max_repYR); 
            continue;
        }
        /***
        *   Max aMine / Ketone repeats
        */
        if(EQSTRING(bufS,"MAX_MKREP",9))
        { 
            sscanf(bufS,"%*s %d",&iconPO->max_repMK); 
            continue;
        }
        if(EQSTRING(bufS,"MAX_KMREP",9))
        { 
            sscanf(bufS,"%*s %d",&iconPO->max_repKM); 
            continue;
        }
        /***
        *   Default = what's this?
        */
        PROBLINE;
        printf("Unrecognized keyword in line %d\n",l);
        fputs(bufS,stdout);
        return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*   Prepare some run-time settings based on parameters
*/
int PrepareInConsI(IN_CONS *iconPO)
{
    int i;

    VALIDATE(iconPO,IN_CONS_ID);
    /***
    *   Initialize to no constraints 
    */
    iconPO->cont_con = 0;
    iconPO->row_con = 0;
    /***
    *   Any min content contraints?
    */
    if(iconPO->minAc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minCc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minGc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minTc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minSc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minWc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minRc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minYc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minMc>0.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->minKc>0.0)
    {   iconPO->cont_con += 1;  }
    /***
    *   Any max content contraints?
    */
    if(iconPO->maxAc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxCc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxGc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxTc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxSc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxWc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxRc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxYc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxKc<1.0)
    {   iconPO->cont_con += 1;  }
    if(iconPO->maxMc<1.0)
    {   iconPO->cont_con += 1;  }
    /***
    *   Set max values and length dependent limits on content
    *   Can never have more than max in sequence 
    */
    iconPO->max_A = FLOOR(RNUM(iconPO->maxlen * iconPO->maxAc));
    iconPO->max_C = FLOOR(RNUM(iconPO->maxlen * iconPO->maxCc));
    iconPO->max_G = FLOOR(RNUM(iconPO->maxlen * iconPO->maxGc));
    iconPO->max_T = FLOOR(RNUM(iconPO->maxlen * iconPO->maxTc));
    iconPO->max_S = FLOOR(RNUM(iconPO->maxlen * iconPO->maxSc));
    iconPO->max_W = FLOOR(RNUM(iconPO->maxlen * iconPO->maxWc));
    iconPO->max_R = FLOOR(RNUM(iconPO->maxlen * iconPO->maxRc));
    iconPO->max_Y = FLOOR(RNUM(iconPO->maxlen * iconPO->maxYc));
    iconPO->max_M = FLOOR(RNUM(iconPO->maxlen * iconPO->maxMc));
    iconPO->max_K = FLOOR(RNUM(iconPO->maxlen * iconPO->maxKc));
    iconPO->min_A[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minAc));
    iconPO->min_C[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minCc));
    iconPO->min_G[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minGc));
    iconPO->min_T[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minTc));
    iconPO->min_S[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minSc));
    iconPO->min_W[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minWc));
    iconPO->min_R[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minRc));
    iconPO->min_Y[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minYc));
    iconPO->min_M[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minMc));
    iconPO->min_K[iconPO->maxlen-1]=CEIL(RNUM(iconPO->maxlen * iconPO->minKc));
    for(i=(iconPO->maxlen-2);i>=0;i--)
    {
        iconPO->min_A[i] = iconPO->min_A[i+1] - 1;
        iconPO->min_C[i] = iconPO->min_C[i+1] - 1;
        iconPO->min_G[i] = iconPO->min_G[i+1] - 1;
        iconPO->min_T[i] = iconPO->min_T[i+1] - 1;
        iconPO->min_S[i] = iconPO->min_S[i+1] - 1;
        iconPO->min_W[i] = iconPO->min_W[i+1] - 1;
        iconPO->min_R[i] = iconPO->min_R[i+1] - 1;
        iconPO->min_Y[i] = iconPO->min_Y[i+1] - 1;
        iconPO->min_M[i] = iconPO->min_M[i+1] - 1;
        iconPO->min_K[i] = iconPO->min_K[i+1] - 1;
    }
    /***
    *   Any row constraints?
    */
    if(IS_BOG(iconPO->max_rowA)) 
        iconPO->max_rowA = iconPO->maxlen;  
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowC)) 
        iconPO->max_rowC = iconPO->maxlen;  
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowG)) 
        iconPO->max_rowG = iconPO->maxlen;  
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowT)) 
        iconPO->max_rowT = iconPO->maxlen; 
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowS)) 
        iconPO->max_rowS = iconPO->maxlen;
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowW)) 
        iconPO->max_rowW = iconPO->maxlen;
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowR)) 
        iconPO->max_rowR = iconPO->maxlen;  
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowY)) 
        iconPO->max_rowY = iconPO->maxlen; 
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowK)) 
        iconPO->max_rowK = iconPO->maxlen;
    else 
        iconPO->row_con = TRUE;
    if(IS_BOG(iconPO->max_rowM)) 
        iconPO->max_rowM = iconPO->maxlen;
    else 
        iconPO->row_con = TRUE;
    return(TRUE);
}
/****************************************************************************
*   Check that seqs can actually be made
*/
int ConsistInConsI(IN_CONS *iconPO)
{
    VALIDATE(iconPO,IN_CONS_ID);
    if((iconPO->maxlen < 1)||(iconPO->maxlen > IN_CONS_MAX))
    {
        PROBLINE;
        printf("Intrinsic constraint evalation limited to %d\n",IN_CONS_MAX);
        printf("Length %d won't work\n",iconPO->maxlen);
        return(FALSE);
    }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minAc,iconPO->maxAc,
        "A content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minCc,iconPO->maxCc,
        "C content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minGc,iconPO->maxGc,
        "G content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minTc,iconPO->maxTc,
        "T content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minSc,iconPO->maxSc,
        "S content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minWc,iconPO->maxWc,
        "W content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minRc,iconPO->maxRc,
        "R content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minYc,iconPO->maxYc,
        "Y content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minKc,iconPO->maxKc,
        "K content"))
    {   return(FALSE);  }
    if(!ConsistPercentParsI(iconPO->maxlen,iconPO->minMc,iconPO->maxMc,
        "M content"))
    {   return(FALSE);  }
    return(TRUE);
}
/****************************************************************************
*   Check that the percents specified can actually be acheived
*/
int ConsistPercentParsI(int len,REAL minR,REAL maxR,char *whatS)
{
    int min,max;

    if((minR<0.0)||(maxR>1.0)||(minR>maxR))
    {
        PROBLINE;
        printf("Inconsistent range specified for %s\n",whatS);
        printf("Min %f  Max %f\n",minR,maxR);
        return(FALSE);
    }
    max = FLOOR(RNUM(RNUM(len) * maxR));
    min = CEIL(RNUM(RNUM(len) * minR));
    if(max < min)
    {
        PROBLINE;
        printf("Impossible range specified for %s\n",whatS);
        printf("Min %f  Max %f\n",minR,maxR);
        printf("No whole numbers obtain this percentage with length %d\n",len);
        return(FALSE);
    }
    return(TRUE);
}
/****************************************************************************
*   Dumps intrinsic constraint settings file
*/
void DumpInCons(IN_CONS *iconPO,FILE *oPF)
{
    int m;

    VALIDATE(iconPO,IN_CONS_ID);
    HAND_NFILE(oPF);
    fprintf(oPF,"#\n");
    fprintf(oPF,"# Intrinsic Sequence Constraints\n");
    if(!NO_S(iconPO->parname))
        fprintf(oPF,"# Parameter file: %s\n",iconPO->parname);
    m = 0;
    m += ReportMinMaxCompI("#  A%        ",iconPO->minAc,iconPO->maxAc,oPF);
    m += ReportMinMaxCompI("#  C%        ",iconPO->minCc,iconPO->maxCc,oPF);
    m += ReportMinMaxCompI("#  G%        ",iconPO->minGc,iconPO->maxGc,oPF);
    m += ReportMinMaxCompI("#  T%        ",iconPO->minTc,iconPO->maxTc,oPF);
    m += ReportMinMaxCompI("#  S%        ",iconPO->minSc,iconPO->maxSc,oPF);
    m += ReportMinMaxCompI("#  W%        ",iconPO->minWc,iconPO->maxWc,oPF);
    m += ReportMinMaxCompI("#  R%        ",iconPO->minRc,iconPO->maxRc,oPF);
    m += ReportMinMaxCompI("#  Y%        ",iconPO->minYc,iconPO->maxYc,oPF);
    m += ReportMinMaxCompI("#  M%        ",iconPO->minMc,iconPO->maxMc,oPF);
    m += ReportMinMaxCompI("#  K%        ",iconPO->minKc,iconPO->maxKc,oPF);
    m += ReportMaxRowI("#  Max A row ",iconPO->max_rowA,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max C row ",iconPO->max_rowC,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max G row ",iconPO->max_rowG,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max T row ",iconPO->max_rowT,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max S row ",iconPO->max_rowS,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max W row ",iconPO->max_rowW,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max R row ",iconPO->max_rowR,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max Y row ",iconPO->max_rowY,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max M row ",iconPO->max_rowM,iconPO->maxlen,oPF);
    m += ReportMaxRowI("#  Max K row ",iconPO->max_rowK,iconPO->maxlen,oPF);
    if(!m)
        fprintf(oPF,"# No constraints specified; anything goes!\n");
    fprintf(oPF,"# \n");
}
/************************************************************************/
int ReportMinMaxCompI(char *nameS,REAL minR,REAL maxR,FILE *oPF)
{
    HAND_NFILE(oPF);
    if( (minR>0.0) || (maxR<1.0) )
    {
        fprintf(oPF,"%s %5.2f to %5.2f\n",nameS,minR,maxR);
        return(TRUE);
    }
    return(FALSE);
}
/************************************************************************/
int ReportMaxRowI(char *nameS,int max,int mlen,FILE *oPF)
{
    HAND_NFILE(oPF);
    if(max<mlen)
    {
        fprintf(oPF,"%s %d\n",nameS,max);
        return(TRUE);
    }
    return(FALSE);
}
/********************************************************************** ooo
*   Checks current sequence against constraints
*/
int SeqInConsOkI(char *seqS, int len, int full, IN_CONS *iconPO)
{
    /***
    *   Check content constraints?
    */
    if(iconPO->cont_con)
    {
        if(!SeqContConsOkI(seqS,len,full,iconPO))
        {   
            DB_CONS DB_PrI("<< SeqInConsOkI FALSE; !ContOK\n");
            return(FALSE);  
        }
    }
    /***
    *   Check for rows of base type constraints?
    */
    if(iconPO->row_con)
    {
        if(!SeqRowConsOkI(seqS,len,iconPO))
        {   
            DB_CONS DB_PrI("<< SeqInConsOkI FALSE; !RowsOK\n");
            return(FALSE);  
        }
    }
    return(TRUE);
}
/****************************************************************************
*   Checks current sequence against constraints for full-sequence
*/
int FullSeqInConsOkI(char *seqS,int len,int full,IN_CONS *iconPO)
{
    DB_CONS 
    { 
        DB_PrI(">> FullSeqInConsOkI |");
        PrintString(seqS,len,NULL);
        DB_PrI("|\n");  
    } 
    if(!SeqInConsOkI(seqS,len,full,iconPO))
    {
        DB_CONS DB_PrI("<< FullSeqInConsOkI FALSE\n");
        return(FALSE);
    }
    DB_CONS DB_PrI("<< FullSeqInConsOkI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Checks current sequence against content constraints
*   If full is TRUE, use absolute limits
*/
int SeqContConsOkI(char *seqS, int len, int full, IN_CONS *iconPO)
{
    int i,s,w,r,y,k,m,con[DNA_NUM];
    REAL ilenR;

    DB_CONS
    {
        DB_PrI(">> SeqContConsOkI len=%d ",len);
        PrintString(seqS,len,NULL);DB_PrI("\n");    
    }   
    /***
    *   Tally the sham
    */
    con[0]=con[1]=con[2]=con[3]=0;
    for(i=0;i<len;i++)
    {
        switch(seqS[i])
        {
            case 'A':   case 'a':   con[DNA_A] += 1;    break;
            case 'C':   case 'c':   con[DNA_C] += 1;    break;
            case 'G':   case 'g':   con[DNA_G] += 1;    break;
            case 'T':   case 't':   con[DNA_T] += 1;    break;
        }
    }
    /***
    *   Any simple comp exceeded?
    */
    ilenR = 1.0 / RNUM(len);
    if(full)
    {
        if( (RNUM(con[DNA_A])*ilenR) > iconPO->maxAc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_C])*ilenR) > iconPO->maxCc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_G])*ilenR) > iconPO->maxGc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_T])*ilenR) > iconPO->maxTc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_A])*ilenR) < iconPO->minAc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_C])*ilenR) < iconPO->minCc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_G])*ilenR) < iconPO->minGc)
        {   return(FALSE);  }
        if( (RNUM(con[DNA_T])*ilenR) < iconPO->minTc)
        {   return(FALSE);  }
    }
    else
    {
        if(con[DNA_A] > iconPO->max_A)
        {   return(FALSE);  }
        if(con[DNA_C] > iconPO->max_C)
        {   return(FALSE);  }
        if(con[DNA_G] > iconPO->max_G)
        {   return(FALSE);  }
        if(con[DNA_T] > iconPO->max_T)
        {   return(FALSE);  }
        len -= 1;
        if(con[DNA_A] < iconPO->min_A[len])
        {   return(FALSE);  }
        if(con[DNA_C] < iconPO->min_C[len])
        {   return(FALSE);  }
        if(con[DNA_G] < iconPO->min_G[len])
        {   return(FALSE);  }
        if(con[DNA_T] < iconPO->min_T[len])
        { return(FALSE);    }
    }
    /***
    *   Any composite max exceeded?
    */
    s = con[DNA_C] + con[DNA_G];
    w = con[DNA_A] + con[DNA_T];
    r = con[DNA_A] + con[DNA_G];
    y = con[DNA_C] + con[DNA_T];
    k = con[DNA_G] + con[DNA_T];
    m = con[DNA_A] + con[DNA_C];
    /***
    *   Any composite min exceeded?
    */
    if(full)
    {
        if( (RNUM(s)*ilenR) > iconPO->maxSc)
        {   return(FALSE);  }
        if( (RNUM(w)*ilenR) > iconPO->maxWc)
        {   return(FALSE);  }
        if( (RNUM(r)*ilenR) > iconPO->maxRc)
        {   return(FALSE);  }
        if( (RNUM(y)*ilenR) > iconPO->maxYc)
        {   return(FALSE);  }
        if( (RNUM(k)*ilenR) > iconPO->maxKc)
        {   return(FALSE);  }
        if( (RNUM(m)*ilenR) > iconPO->maxMc)
        {   return(FALSE);  }
        if( (RNUM(s)*ilenR) < iconPO->minSc)
        {   return(FALSE);  }
        if( (RNUM(w)*ilenR) < iconPO->minWc)
        {   return(FALSE);  }
        if( (RNUM(r)*ilenR) < iconPO->minRc)
        {   return(FALSE);  }
        if( (RNUM(y)*ilenR) < iconPO->minYc)
        {   return(FALSE);  }
        if( (RNUM(k)*ilenR) < iconPO->minKc)
        {   return(FALSE);  }
        if( (RNUM(m)*ilenR) < iconPO->minMc)
        {   return(FALSE);  }
    }
    else
    {
        if(s > iconPO->max_S)
        {   return(FALSE);  }
        if(w > iconPO->max_W)
        {   return(FALSE);  }
        if(r > iconPO->max_R)
        {   return(FALSE);  }
        if(y > iconPO->max_Y)
        {   return(FALSE);  }
        if(k > iconPO->max_K)
        {   return(FALSE);  }
        if(m > iconPO->max_M)
        {   return(FALSE);  }
        if(s < iconPO->min_S[len])
        {   return(FALSE);  }
        if(w < iconPO->min_W[len])
        {   return(FALSE);  }
        if(r < iconPO->min_R[len])
        {   return(FALSE);  }
        if(y < iconPO->min_Y[len])
        {   return(FALSE);  }
        if(k < iconPO->min_K[len])
        {   return(FALSE);  }
        if(m < iconPO->min_M[len])
        {   return(FALSE);  }
    }
    /***
    *   Composition ok
    */
    return(TRUE);
}
/****************************************************************************
*   Checks current sequence against row-of-base-type constraints
*/
int SeqRowConsOkI(char *seqS, int len, IN_CONS *iconPO)
{
    int i;
    int rowA,mrowA,rowC,mrowC,rowG,mrowG,rowT,mrowT;
    int rowS,mrowS,rowW,mrowW,rowM,mrowM,rowK,mrowK,rowY,mrowY,rowR,mrowR;

    DB_CONS 
    {
        DB_PrI(">> SeqRowConsOkI len=%d ",len);
        PrintString(seqS,len,NULL);DB_PrI("\n");    
    }
    rowA=mrowA=rowC=mrowC=rowG=mrowG=rowT=mrowT=0;
    rowS=mrowS=rowW=mrowW=rowM=mrowM=rowK=mrowK=rowY=mrowY=rowR=mrowR=0;
    for(i=0;i<len;i++)
    {
        switch(seqS[i])
        {
            case 'A': case 'a':
                rowA++; rowW++; rowR++; rowM++;
                /***
                *   Set all non-A to zero after keeping the max
                */
                mrowC = MAX_NUM(rowC,mrowC);
                mrowG = MAX_NUM(rowG,mrowG);
                mrowT = MAX_NUM(rowT,mrowT);
                rowC = rowG = rowT = 0;
                mrowS = MAX_NUM(rowS,mrowS);
                mrowY = MAX_NUM(rowY,mrowY);
                mrowK = MAX_NUM(rowK,mrowK);
                rowS = rowY = rowK = 0;
                break;
            case 'C': case 'c':
                rowC++; rowS++; rowY++; rowM++;
                /***
                *   Set all non-C to zero after keeping the max
                */
                mrowA = MAX_NUM(rowA,mrowA);
                mrowG = MAX_NUM(rowG,mrowG);
                mrowT = MAX_NUM(rowT,mrowT);
                rowA = rowG = rowT = 0;
                mrowW = MAX_NUM(rowW,mrowW);
                mrowR = MAX_NUM(rowR,mrowR);
                mrowK = MAX_NUM(rowK,mrowK);
                rowW = rowR = rowK = 0;
                break;
            case 'G': case 'g':
                rowG++; rowS++; rowR++; rowK++;
                /***
                *   Set all non-A to zero after keeping the max
                */
                mrowA = MAX_NUM(rowA,mrowA);
                mrowC = MAX_NUM(rowC,mrowC);
                mrowT = MAX_NUM(rowT,mrowT);
                rowA = rowC = rowT = 0;
                mrowW = MAX_NUM(rowW,mrowW);
                mrowY = MAX_NUM(rowY,mrowY);
                mrowM = MAX_NUM(rowM,mrowM);
                rowW = rowY = rowM = 0;
                break;
            case 'T': case 't':
                rowT++; rowW++; rowY++; rowK++;
                /***
                *   Set all non-A to zero after keeping the max
                */
                mrowA = MAX_NUM(rowA,mrowA);
                mrowC = MAX_NUM(rowC,mrowC);
                mrowG = MAX_NUM(rowG,mrowG);
                rowA = rowC = rowG = 0;
                mrowS = MAX_NUM(rowS,mrowS);
                mrowR = MAX_NUM(rowR,mrowR);
                mrowM = MAX_NUM(rowM,mrowM);
                rowS = rowR = rowM = 0;
                break;
        }
    }
    /***
    *   Update maxs for everything considered
    */
    mrowA = MAX_NUM(rowA,mrowA);
    mrowC = MAX_NUM(rowC,mrowC);
    mrowG = MAX_NUM(rowG,mrowG);
    mrowT = MAX_NUM(rowT,mrowT);
    mrowS = MAX_NUM(rowS,mrowS);
    mrowW = MAX_NUM(rowW,mrowW);
    mrowR = MAX_NUM(rowR,mrowR);
    mrowY = MAX_NUM(rowY,mrowY);
    mrowM = MAX_NUM(rowM,mrowM);
    mrowK = MAX_NUM(rowK,mrowK);
    DB_CONS DB_PrI(" S=%d W=%d R=%d Y=%d M=%d K=%d\n",
        mrowS, mrowW, mrowR, mrowY, mrowM, mrowK);
    /***
    *   Check counts to constraints
    */
    if(mrowA > iconPO->max_rowA)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE A\n");
        return(FALSE);
    }
    if(mrowC > iconPO->max_rowC)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE C\n");
        return(FALSE);
    }
    if(mrowG > iconPO->max_rowG)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE G\n");
        return(FALSE);
    }
    if(mrowT > iconPO->max_rowT)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE T\n");
        return(FALSE);
    }
    if(mrowS > iconPO->max_rowS)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE S\n");
        return(FALSE);
    }
    if(mrowW > iconPO->max_rowW)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE W\n");
        return(FALSE);
    }
    if(mrowR > iconPO->max_rowR)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE R\n");
        return(FALSE);
    }
    if(mrowY > iconPO->max_rowY)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE Y\n");
        return(FALSE);
    }
    if(mrowM > iconPO->max_rowM)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE M\n");
        return(FALSE);
    }
    if(mrowK > iconPO->max_rowK)
    {
        DB_CONS DB_PrI("<< SeqRowConsOkI FALSE K\n");
        return(FALSE);
    }
    DB_CONS DB_PrI("<< SeqRowConsOkI TRUE\n");
    return(TRUE);
}
/*****************************************************************************
*   Set maximum length
*/
int SetInConsMaxlenI(IN_CONS *iconPO, int len)
{
    VALIDATE(iconPO,IN_CONS_ID);
    if( (len<1) || (len>=IN_CONS_MAX) )
        return(FALSE);
    iconPO->maxlen = len;
    return(TRUE);
}
/*****************************************************************************
*   Get maximum length
*/
int GetInConsMaxlenI(IN_CONS *iconPO)
{
    VALIDATE(iconPO,IN_CONS_ID);
    return(iconPO->maxlen);
}
