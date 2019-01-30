/*
* dna_trim.c
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
#include <ctype.h>
#include "prim.h"
#include "dna.h"


/***************************************************************************
*   Handle various trimming / modifications of sequence.
*/
int DnaTrimSeqI(SEQTRIM *stPO, SEQ *seqPO)
{
    int slen, wlen, wcent, start, end, t;

    VALIDATE(stPO, SEQTRIM_ID);
    if( !AnySeqTrimmingI(stPO) ) {
        return(TRUE);
    }
    slen = GetSeqLenI(seqPO);
    start = 0;
    end = slen;
    /***
    *   Base range
    */
    if( (stPO->basf_s >= 0.0) || (stPO->basf_e >= 0.0) ){
        start = ROUND(stPO->basf_s * DNUM(slen));
        end = ROUND(stPO->basf_e * DNUM(slen));
    }
    else if( (stPO->base_s >= 0) || (stPO->base_e >= 0) ){
        start = stPO->base_s;
        end = stPO->base_e;
    }
    /***
    *   Trimming ends 
    */
    else if( (stPO->trim_s >= 0.0) || (stPO->trim_e >= 0.0) ){
        start = stPO->trim_s + 1;
        end = slen - stPO->trim_e;
    }
    /***
    *   Window
    */
    else if( (stPO->wind > 0) || (stPO->winf > 0.0) ){
        if (stPO->winf > 0.0) {
            wlen = ROUND(stPO->winf * DNUM(slen));
        }
        else {
            wlen = stPO->wind;
        }
        /* Center position */
        if (stPO->wcenf >= 0) {
            wcent = ROUND(stPO->wcenf * DNUM(slen));
        }
        else if (stPO->wcent > 0) {
            wcent = stPO->wcent;
        }
        else {
            wcent = slen / 2;
        }
        /***
        *   Start / end
        */
        start = wcent - wlen / 2 + 1;
        end = start + wlen - 1;
    }
    /***
    *   If relative to end, invert
    */
    if(stPO->rre) {
        t = start;
        start = slen - end + 1;
        end = slen - t + 1;
    }
    /***
    *   Possibly adjust coords and limit the sham to within bounds
    */
    if(stPO->one_base) {
        start--;
        end--;
    }
    if(stPO->end_in) {
        end++;
    }
    LIMIT_NUM(start,0,slen);
    LIMIT_NUM(end,0,slen);
    /***
    *   Mask?
    */
    if (stPO->nmask) {
        SetMaskSeqSubseqI(seqPO, start, end);
    }
    else if (stPO->umask) {
        SetMaskSeqSubseqI(seqPO, 0, start);
        SetMaskSeqSubseqI(seqPO, end, -1);
    }
    /* Change case */
    else if (stPO->ucase) {
        SetCaseSeqSubseqI(seqPO, FALSE, -1, -1);
        SetCaseSeqSubseqI(seqPO, TRUE, start, end);
    }
    else if (stPO->lcase) {
        SetCaseSeqSubseqI(seqPO, TRUE, -1, -1);
        SetCaseSeqSubseqI(seqPO, FALSE, start, end);
    }
    /* Trim */
    else {
        slen = end - start;
        NarrowSeqI(seqPO, start, slen, FORWARD, FALSE);
    }
    return(TRUE);
}
/***************************************************************************
*   Are any trimseq settings non-default (i.e. active?)
*/
int AnySeqTrimmingI(SEQTRIM *stPO)
{
    VALIDATE(stPO, SEQTRIM_ID);
    if( (stPO->base_s >= 0) || (stPO->base_s >= 0) ||
        (stPO->basf_s >= 0.0) || (stPO->basf_s >= 0.0) ||
        (stPO->trim_s >= 0) || (stPO->trim_e >= 0) ||
        (stPO->wind >= 0) || (stPO->winf >= 0.0) ) 
    {
            return(TRUE);
    }
    return(FALSE);
}
