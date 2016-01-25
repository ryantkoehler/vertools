/*
* dna_prob.c
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
#include <string.h>
#include <ctype.h>
#include "prim.h"
#include "dna.h"
#include "seq_info.h"
#include "dna_util.h"


/**************************************************************************
*   Check probe-replated options
*/
int CheckDnuProbeOptionsI(DNA_UTIL *duPO)
{
    if( !NO_S(duPO->olpname) ) {
        duPO->owhat = DNUO_PROBE;
        if( (duPO->olp_up<0) || (duPO->olp_dn<0) ) {
            PROBLINE;
            printf("Bad olp windows specified: %d up and %d down\n",
                duPO->olp_up,duPO->olp_dn);
            return(FALSE);
        }
    }
    else if(duPO->do_pol>0) {
        duPO->owhat = DNUO_PROBE;
        if( (duPO->opl_st<0) || (duPO->opl_j<0) ) {
            PROBLINE;
            printf("Bad probe step %d and/or jump %d\n",
                duPO->opl_st,duPO->opl_j);
            return(FALSE);
        }
        if(duPO->do_pml) { 
            if(duPO->do_pml < duPO->do_pol) {
                PROBLINE;
                printf("Bad multi-length probe range %d to %d\n",
                    duPO->do_pol, duPO->do_pml);
                return(FALSE);
            }
        }
        LIMIT_NUM(duPO->opl_st,1,TOO_BIG);
        LIMIT_NUM(duPO->opl_j,1,TOO_BIG);
    }
    return(TRUE);
}
/**************************************************************************
*   Output probes to length
*/
int HandleDuProbesOutI(DNA_UTIL *duPO, SEQ *seqPO, FILE *outPF)
{
    int len,plen,plen2,minlen,pos,fpos,flen;
    char bufS[DEF_BS],*seqPC,nameS[NSIZE],pnameS[NSIZE];

    HAND_NFILE(outPF);
    FillSeqNameStringI(seqPO,nameS,NSIZE-1);
    len = GetSeqLenI(seqPO);
    if(!GetSeqSeqI(seqPO,&seqPC)) {
        return(FALSE);
    }
    pos = duPO->opl_st -1;
    plen = duPO->do_pol;
    plen2 = (duPO->do_pml > 0) ? duPO->do_pml : plen;
    minlen = (duPO->do_pat > 0) ? duPO->do_pat : plen - 1;
    /***
    *   Tell story.
    *   Listed or number-based probes
    */
    if( !NO_S(duPO->olpname) ) {
        fprintf(outPF,"# Listed probes from file: %s\n",duPO->olpname);
        /* Make sure to rewind input file */
        BOG_CHECK(! duPO->olp);
        rewind(duPO->olp);
    }
    else {
        if( (plen<1) || (len<minlen) ) {
            printf("# Sequence %s has only %d bases\n",nameS,len);
            printf("# Can't extract probes of length %d!!!\n",duPO->do_pol);    
            return(FALSE);
        }
        if(duPO->do_pml > 0) {
            fprintf(outPF,"# Sampled probes %d to %d long every %d from %d\n",
                plen, plen2, duPO->opl_j, pos + 1);
        }
        else {
            fprintf(outPF,"# Sampled probes %d long every %d from %d\n",
                plen, duPO->opl_j, pos + 1);
        }
    }
    if( duPO->olp_up ) {
        fprintf(outPF,"#    Extra window of %d up and %d down added\n", 
            duPO->olp_up, duPO->olp_dn);
    }
    /***    
    *   Walk along the sham
    */
    while(TRUE)
    {
        /***
        *   Reading from file?
        */
        if(duPO->olp) {
            if(!fgets(bufS,LINEGRAB,duPO->olp)) {
                break;
            }
            if(COM_LINE(bufS)) {
                continue;
            }
            plen = pos = 0;
            sscanf(bufS,"%d %d",&plen,&pos);
            if( (plen<1) || (pos<1) ) {
                continue;
            }
            pos -= 1;
            minlen = plen;
        }
        /***
        *   Adjust coords by windows?
        */
        fpos = pos - duPO->olp_up;
        LIMIT_NUM(fpos,0,len-1);
        flen = plen + duPO->olp_up + duPO->olp_dn;
        flen = MIN_NUM(flen,len-fpos);
        plen2 = (duPO->do_pml > 0) ? duPO->do_pml : flen;
/*
printf("xxx up=%d dn=%d fpos=%d flen=%d\n", duPO->olp_up,duPO->olp_dn,fpos,flen);
*/
        /***
        *   Supplied file coords 
        */
        if(duPO->olp) {
            if( (fpos < 0) || (flen < plen) || ((fpos+flen) > len) ) {
                printf("# Sequence %s is length %d\n",nameS,len);
                printf("# Can't extract probe of %d bases starting at %d\n",
                    plen, fpos+1);
                printf("#   Input <len> <start> :\t"); fputs(bufS,stdout);
                continue;
            }
            sprintf(pnameS,"%s__%d_%d",nameS, fpos + 1, fpos + flen);
            HandleDuOneProbeOutI(duPO, &seqPC[fpos], flen, pnameS, outPF);
        }
        /***
        *   Stepping coords
        */
        else 
        {
/*
printf("xxx flen=%d plen2=%d fpos+flen=%d\n", flen,plen2,fpos+flen);
*/
            /* as "flen < minlen" for V0.64 */
            if(flen <= minlen) {
                break;
            }
            /***
            *   Dump probe(s)
            */
            while( (flen <= plen2) && ((fpos + flen) <= len) ){
                sprintf(pnameS,"%s__%d_%d",nameS, fpos + 1, fpos + flen);
                HandleDuOneProbeOutI(duPO, &seqPC[fpos], flen, pnameS, outPF);
                flen++;
            }
        }
        pos += duPO->opl_j;
    }
    return(TRUE);
}
/**************************************************************************
*   Dump one probe
*/
int HandleDuOneProbeOutI(DNA_UTIL *duPO, char *seqPC, int len, char *pnameS, FILE *outPF)
{
    HAND_NFILE(outPF);
    /***
    *   Check ambigs / SNPs
    */
    if( duPO->do_pco && (CountSeqAmbigsI(seqPC,0,len)>0) ) {
        printf("# Ambiguous bases in %s; Skipping\n",pnameS);
    }
    else if( duPO->do_pco && (CountSeqSnpSitesI(seqPC,0,len)>0) ) {
        printf("# SNPs in %s; Skipping\n",pnameS);
    } 
    else {
        if(duPO->oform == SEQFM_FASTA) {
            fprintf(outPF,"> %s\n",pnameS);
        }   
        else {
            ReplaceChars(' ',pnameS,'_',pnameS);
            fprintf(outPF,"%s\t",pnameS);
        }
        PrintString(seqPC,len,outPF);
        fprintf(outPF,"\n");
    }
    return(TRUE);
}
