/*
* venpipe.c
*
* Copyright 2016 Ryan Koehler, VerdAscend Sciences, ryan@verdascend.com
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
#include <ctype.h>
#include <string.h>
#define __MAIN__
#include "prim.h"
#include "venpipe.h"

#define DB_VENG if(DB[142])

/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit(AllDoneI(VenPipeI(argc,argv),NULL));    }
/**************************************************************************/
void VenPipeUse()
{
    VersionSplash(NULL,VEN_VERSION_S,"#  ",TRUE);
    printf("# %s\n",VLIB_VERSION_S);
    printf("#\n");
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <input>     Input sequence file (format guessed from .ext)\n");
    printf("   -iraw       Treat input as \"raw\" format; <name> <seq> / line\n");
    printf("   -iseq       Treat input as simmple sequence; <seq> / line\n");
    printf("   -ifas       Treat input as fasta format\n");
    printf("   -vpar XXX   Parameter file (Vienna energy table)\n");
    printf("   -nosc       No salt correction for vienna parameter file\n");
    printf("   -out XXX    Set output to = XXX\n");
    printf("   -ofas       Output fasta file\n");
    printf("   -oraw       Output \"raw\" file\n");
    printf("   -bran # #   Restrict evaluation to bases # to #\n");
    printf("   -rre        Base range relative to end; i.e. backwards\n");
    printf("   -mask       Mask non-selected range with N (default = cut)\n");
    printf("   -not        Invert mask; selected range masked with N\n");
    printf("   -temp #     Set temperature to # (C')\n");
    printf("   -sal #      Set salt to # (Molar)\n");
    printf("   -ksapar     Keep salt-corrected parameter file\n");
    printf("   -pfe        Partition Free Energy; Default is Min Free Eng\n");
    printf("   -dss        Dump sequence structure\n");
    printf("   -dmb        Dump matching (paired) base counts\n");
    printf("   -mbtab      Dump matching base table (0 = free; >0 = matching)\n");
    printf("   -mbr # #    Matching base count range from # to #\n");
    printf("   -mrre       Matching base count range relative to end\n");
    printf("   -dseq       Dump (report) sequences appended as last column\n");
    printf("   -melt # #   Dump melting profile from temp # to #\n");
    printf("   -mstep #    Set melting profile step to # (default %4.2f)\n",
        DEF_MELT_JUMP_D);
}
/**************************************************************************
*   top level function
*/
int VenPipeI(int argc, char **argv)
{
    int ok,nprob,nseq,iraw,iseq,ifas,com;
    VENPIPE *vpPO;

    if(!(vpPO = CreateVenpipePO()))
    {
        printf("Failed to allocate venpipe data structure\n");
        ABORTLINE; return(FALSE);
    }
    iraw = iseq = ifas = com = FALSE;
    if(!ParseArgsI(argc, argv,
        "S -vpar S -out S -temp D -com B -dss B -pfe B -pcon D -tcon D -sal D\
        -ofas B -oraw B -dmb B -bran I2\
        -rre B -mbr I2 -mrre B -melt D2 -mstep D -iraw B -ifas B\
        -nosc B -ksapar B -mask B -not B -mbtab B -dseq B -iseq B",
        vpPO->inname, vpPO->vparfile, vpPO->outname, &vpPO->temp, &com, 
        &vpPO->do_ss, &vpPO->do_pfe, &vpPO->pcon, &vpPO->tcon, &vpPO->salt,
        &vpPO->ofas, &vpPO->oraw, &vpPO->do_dmb, &vpPO->firstb,&vpPO->lastb,
        &vpPO->rre, &vpPO->dmb_f,&vpPO->dmb_l, &vpPO->mrre, 
        &vpPO->mst,&vpPO->men, &vpPO->mj, &iraw, &ifas, &vpPO->do_saltcorrect,
        &vpPO->do_ksapar, &vpPO->do_mask, &vpPO->do_not, &vpPO->do_mbtab,
        &vpPO->do_ds, &iseq,
        (int *)NULL))
    {
        VenPipeUse();
        CHECK_VENPIPE(vpPO);
        return(FALSE);
    }
    /***
    *   Set intput format
    */
    vpPO->iform = FigureSeqFileTypeI(iraw, iseq, ifas, vpPO->inname, TRUE);
    if(!vpPO->iform) {
        printf("Problem with input seq(s)\n");
        CHECK_VENPIPE(vpPO);
        return(FALSE);
    }
    /***
    *   Check options aand set everything up to use
    */
    if(!OkVenpipeOptionsI(vpPO)) {
        ABORTLINE;
        CHECK_VENPIPE(vpPO);
        return(FALSE);
    }
    if(!RealizeVenpipeI(vpPO)) {
        ABORTLINE;
        CHECK_VENPIPE(vpPO);
        return(FALSE);
    }
    WriteVenpipeHeader(vpPO,vpPO->outname,vpPO->out);
    /***
    *   Process each record listed in input list or single record
    */
    nprob = nseq = 0;
    while(TRUE)
    {
        /***
        *   Parse sequence; FALSE = done
        */
        ok = ParseSeqI(vpPO->in, vpPO->iform, nseq+1, vpPO->cleanseq, TRUE, vpPO->seq);
        if(ok==FALSE) {
            break;
        }
        nseq++;
        /***
        *   Init vienna struct (not fully), and copy sequence 
        */
        if(ok==TRUE) {
            InitVenpipeI(vpPO, FALSE);
            ok = CopyWorkingSeqI(vpPO, vpPO->seq);
        }
        if(!ok) {
            nprob++;
            continue;
        }
        /***
        *   If melting, only one sequence is done
        */
        if(!BAD_REAL(vpPO->mst)) {
            HandleVenpipeMeltingI(vpPO, vpPO->out);
            break;
        }
        else {
            if(!GetVenpipeStructureI(vpPO)) {
                PROBLINE;
                printf("Failed to seq structure for seq # %d\n",nseq);
                break;
            }
            HandleVenpipeOut(vpPO, vpPO->out);
        }
    }
    /***
    *   All done
    */
    CHECK_VENPIPE(vpPO);
    return(TRUE);
}
/***************************************************************************
*   Check options
*/
int OkVenpipeOptionsI(VENPIPE *vpPO)
{
    if(NO_S(vpPO->vparfile)) {
        WARNLINE;
        printf("Vienna energy file should probably be used... Defaults = ???\n");
    }
    if( (!BAD_INT(vpPO->dmb_f)) || (vpPO->do_mbtab) ){
        vpPO->do_dmb = TRUE;
    }
    if(!BAD_INT(vpPO->firstb)) {
        if( (vpPO->firstb<1) || (vpPO->firstb>=vpPO->lastb) )
        {
            printf("Bad base range given: %d to %d\n",vpPO->firstb,vpPO->lastb);
            return(FALSE);
        }
    }
    if(!BAD_INT(vpPO->dmb_f)) {
        if( (vpPO->dmb_f < 1) || (vpPO->dmb_f > vpPO->dmb_l) )
        {
            PROBLINE;
            printf("Bad dump base range given: %d to %d\n",vpPO->dmb_f,vpPO->dmb_l);
            return(FALSE);
        }
    }
    /***
    *   Melting 
    */
    if(!BAD_REAL(vpPO->mst)) {
        if( (vpPO->mst < 0.0) || (vpPO->mst > 100.0) || 
            (vpPO->men < 0.0) || (vpPO->men > 100.0) ||
            (ABS_VAL(vpPO->mj) < MIN_MELT_JUMP_D) ) {
            PROBLINE;
            printf("Bad melting temperature numbers (0-100, min jump %2.2f)\n",
                MIN_MELT_JUMP_D);
            return(FALSE);
        }
    }
    return(TRUE);
}
/***************************************************************************
*   Dump a matrix "melting" the current structure
*/
int HandleVenpipeMeltingI(VENPIPE *venPO,FILE *outPF)
{
    int i,n,m,nsham;
    DOUB tD,fD,jumpD;

    HAND_NFILE(outPF);
    nsham = 6;
    /***
    *   Header story 
    */
    fprintf(outPF,"# Melting profile for  %s\n",venPO->tname);
    fprintf(outPF,"# Rows = temperature, columns = bases\n");
    fprintf(outPF,"# Starting temperature %5.2f\n",venPO->mst);
    fprintf(outPF,"# Ending temperature   %5.2f\n",venPO->men);
    fprintf(outPF,"# Temperature step     %5.2f\n",venPO->mj);
    fprintf(outPF,"#\n");
    fprintf(outPF,"# NOTE: The last %d cols show the fraction of paired bases\n",
        nsham);
    fprintf(outPF,"#       This is duplicated to show up on plots (sham!)\n");
    /***
    *   Legend
    */
    fprintf(outPF,"Temp ");
    for(i=0;i<venPO->tlen;i++)
    {
        fprintf(outPF,"%d ",i+1);
    }
    fprintf(outPF,". Frac ");
    for(i=1;i<nsham;i++)
    {
        fprintf(outPF,". ");
    }
    fprintf(outPF,"\n");
    /***
    *   Temperature ramp
    */
    tD = venPO->mst;
    jumpD = venPO->mj;
    LIMIT_NUM(tD,0.0,100.0);
    while(TRUE)
    {
        SetViennaTemperatureI(venPO->ven,tD);
        if(!GetVenpipeStructureI(venPO)) {
            PROBLINE;
            printf("Failed to get structure (temp = %f)\n",tD);
            return(FALSE);
        }
        fprintf(outPF,"%6.2f\t",tD);
        n = strlen(venPO->tss);
        m = 0;
        for(i=0;i<n;i++)
        {
            switch(venPO->tss[i])
            {
                case '.':   
                    fprintf(outPF," 0");    
                    break;
                case '(':   
                case ')':   
                    fprintf(outPF," 1");    
                    m++;
                    break;
                default:    
                    fprintf(outPF," ?");    
            }
        }
        /***
        *   Sham with multiple percents to shows up on plots
        */
        fD = RNUM(m)/RNUM(n);
        fprintf(outPF," 0");
        for(i=0;i<nsham;i++)
        {
            fprintf(outPF," %4.3f",fD);
        }
        fprintf(outPF,"\n");
        /***
        *   Increment and check if still in bounds
        */
        tD += jumpD;
        if( (jumpD<0.0)&&(tD<venPO->men) )
        {
            break;
        }
        if( (jumpD>0.0)&&(tD>venPO->men) )
        {
            break;
        }
    }
    return(TRUE);
}
