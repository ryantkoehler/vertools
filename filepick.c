/*
* filepick.c
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
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define __MAIN__
#include "prim.h"
#include "wordlist.h"
#include "filepick.h"


/**************************************************************************/
int main(int argc, char **argv)
{ Init(argc,argv); exit( AllDoneI(FilePickI(argc,argv),NULL) ); }
/*******************************************************************/
void FilePickUse()
{
    VersionSplash(NULL,VERSION_S,"#  ",TRUE);
    printf("Usage: <infile> ['-' for stdin] [...options]\n");
    printf("   <file>     A keyword delimited file\n");
    printf("   -out XXX   Set output file to XXX\n");
    printf("   -rst XXX   Set record START token to XXX (def = \"%s\")\n", DEF_STOK);
    printf("   -rsany     Record start token may be anywhere (strstr)\n");
    printf("   -ret XXX   Set record END token to XXX\n");
    printf("   -ner       No end record output\n");
    printf("   -snl       Separate records with new line\n");
    printf("   -srn       Separate records with record number line\n");
    printf("   -rg # #    Select subset of records # to #\n");
    printf("   -wlis XXX  Select records matching words in XXX (first token / line)\n");
    printf("   -kc        Keep case for word comparison (default ignore)\n");
    printf("   -wst       Word start only needs to match line (not full token)\n");
    printf("   -wsub      Word substring only needs to match line (not full token)\n");
    printf("   -not       Invert selection criteria\n");
    printf("   -dh        Dump input header (before START token)\n");
    printf("   -oh        Only dump input header (before START token)\n");
    printf("   -esf       Extract records to separate files\n");
    printf("   -esb XXX   Set extract file base name to XXX\n");
    printf("   -esx XXX   Set extract file extension to XXX\n");
    printf("   -stat      Just report the numbers\n");
    printf("\n");
}
/**************************************************************************
*   top level function
*/
int FilePickI(int argc, char **argv)
{
    int stat, record, nokrec, in_rec, ok_rec, out_line, rec_line;
    FILEPICK *fpPO;
    char bufS[BBUFF_SIZE], tempS[DEF_BS], curecS[NSIZE];

    stat = FALSE;
    fpPO=CreateFilepickPO();
    if(!ParseArgsI(argc, argv,
        "S -out S -wlis S -rst S -ret S -rg I2 -not B -kc B -dh B -stat B\
        -rsany B -esf B -esb S -esx S -oh B -snl B -srn B -tst B -tsub B\
        -ner B",
        fpPO->inname, fpPO->outname, fpPO->wlisname, fpPO->stok, fpPO->etok, 
        &fpPO->first,&fpPO->last, &fpPO->do_not, &fpPO->do_kc, &fpPO->do_dh, 
        &stat, &fpPO->do_rsany, &fpPO->do_esf, fpPO->esbase, fpPO->esext,
        &fpPO->do_oh, &fpPO->do_snl, &fpPO->do_srn, 
        &fpPO->do_wst, &fpPO->do_wsub, &fpPO->do_ner,
        (int *)NULL))
    {
        FilePickUse();
        CHECK_FILEPICK(fpPO);
        return(FALSE);
    }
    /***
    *   Set up
    */
    if(!(fpPO->in = OpenUFilePF(fpPO->inname,"r",NULL)))
    {
        ABORTLINE;
        CHECK_FILEPICK(fpPO);
        return(FALSE);
    }
    if(!CheckFilepickOptions(fpPO)) {
        ABORTLINE;
        CHECK_FILEPICK(fpPO);
        return(FALSE);
    }
    /***
    *   Party through input file
    */
    record = nokrec = 0;
    in_rec = ok_rec = 0;
    rec_line = 0;
    while(fgets(bufS,BLINEGRAB,fpPO->in))
    {
        if(LineMatchRecTokenI(fpPO,bufS,fpPO->stok,curecS)) {
            record++;
            in_rec++;
            ok_rec = IsRecordOkI(fpPO,record,curecS);
            if(ok_rec) {
                nokrec++;
            }
            rec_line = 0;
        }
        /***
        *   only header and have record
        */
        if( (fpPO->do_oh) && (record) ) {
            break;
        }
        /***
        *   Match end token? If so, still may print this line
        */
        if(LineMatchRecTokenI(fpPO,bufS,fpPO->etok,NULL)) {
            if(fpPO->do_ner) {
                out_line = FALSE;
            }
            else {
                out_line = FpOutputThisLineI(fpPO,record,in_rec,ok_rec);
            }
            in_rec = 0;
        }
        else {
            out_line = FpOutputThisLineI(fpPO,record,in_rec,ok_rec);
        }
        /***
        *   Output if not only reporting stats
        */
        if( (out_line) && (!stat) ) {
            /***
            *   First line of record, check output and extra stuff to print 
            */
            if(rec_line == 0) {
                if(!HandleFilepickOutfileI(fpPO,record,curecS)) {
                    ABORTLINE;
                    break;
                }
                if(fpPO->do_snl) {
                    FpOutOneLineI(fpPO, "\n", fpPO->out);
                }
                if(fpPO->do_srn) {
                    sprintf(tempS, "# %s Filepick Record %d %s\n", 
                        REC_SEP_S,record,REC_SEP_S);
                    FpOutOneLineI(fpPO, tempS, fpPO->out);
                }
            }
            FpOutOneLineI(fpPO, bufS, fpPO->out);
            rec_line++;
        }
    }
    if(stat)
    {
        printf("# Total %d records\n",record);
        printf("#       %d qualify\n",nokrec);
    }
    /***
    *   All Done
    */
    CHECK_FILEPICK(fpPO);
    return(TRUE);
}
/*************************************************************************/
FILEPICK *CreateFilepickPO()
{
    FILEPICK *fpPO;

    if(!(fpPO=(FILEPICK *)ALLOC(1,sizeof(FILEPICK)))) {
        return(NULL);
    }
    fpPO->ID = FILEPICK_ID;
    InitFilepick(fpPO);
    return(fpPO);
}
/*************************************************************************/
int DestroyFilepickI(FILEPICK *fpPO)
{
    VALIDATE(fpPO,FILEPICK_ID);
    CHECK_WORDLIST(fpPO->wlis);
    CHECK_FILE(fpPO->in);
    CHECK_NFILE(fpPO->out,fpPO->outname);
    FREE(fpPO);
    return(TRUE);
}
/*************************************************************************
*   Init structure
*/
void InitFilepick(FILEPICK *fpPO)
{
    VALIDATE(fpPO,FILEPICK_ID);
    INIT_S(fpPO->inname);
    fpPO->in = NULL;
    INIT_S(fpPO->outname);
    fpPO->out = NULL;
    strcpy(fpPO->stok,DEF_STOK);
    INIT_S(fpPO->etok);
    INIT_S(fpPO->wlisname);
    fpPO->wlis = NULL;
    fpPO->first = fpPO->last = BOGUS;
    fpPO->do_not = FALSE;
    fpPO->do_ner = FALSE;
    fpPO->do_kc = FALSE;
    fpPO->do_dh = FALSE;
    fpPO->do_oh = FALSE;
    fpPO->do_rsany = FALSE;
    fpPO->do_esf = FALSE;
    INIT_S(fpPO->esbase);
    INIT_S(fpPO->esext);
}
/*************************************************************************
*
*/
int CheckFilepickOptions(FILEPICK *fpPO) 
{
    if(fpPO->do_oh) {
        fpPO->do_dh = TRUE;
    }
    if(!NO_S(fpPO->wlisname)) {
        if( ! (fpPO->wlis = CreateWordlistPO(fpPO->wlisname,NSIZE))) {
            PROBLINE;
            printf("Failed to get tokens from %s\n",fpPO->wlisname);
            return(FALSE);
        }
    }
    return(TRUE);
}
/*************************************************************************
*   True if passed line matches based on token
*   Sets record value string (i.e. name) into passed string
*/
int LineMatchRecTokenI(FILEPICK *fpPO, char *bufS, char *tokS, char *nameS)
{
    int n;
    char *cPC, wordS[DEF_BS];

    /***
    *   No token to look for, so no match!
    */
    if( NO_S(tokS) || (BlankStringI(bufS)) ) {
        return(FALSE);
    }
    /***
    *   Look anywhere or at first word on line
    */
    n = 0;
    if (nameS) {
        INIT_S(nameS);
    }
    if(fpPO->do_rsany) {
        cPC = strstr(bufS,tokS);
        if(cPC) {
            n++;
            if (nameS) {
                sscanf(cPC,"%*s %s",nameS);
            }
        }
    }
    else {
        sscanf(bufS,"%s",wordS);
        if(WordStringMatchI(tokS,wordS,fpPO->do_kc,FALSE,FALSE)) {
            n++;
            if (nameS) {
                sscanf(bufS,"%*s %s",nameS);
            }
        }
    }
    return(n);
}
/*************************************************************************
*   Check record status and return if it's qualified or not 
*/
int IsRecordOkI(FILEPICK *fpPO,int record,char *curecS)
{
    int ok;

    ok = TRUE;
    /***
    *   If range set, in range?
    */
    if(!IS_BOG(fpPO->first)) {
        ok = ((record >= fpPO->first) && (record <= fpPO->last)) ? TRUE : FALSE; 
    }
    /***
    *   Subset of names to check?
    */
    if( (ok) && (fpPO->wlis) ) {
        ok = WordInWordlistI(fpPO->wlis, curecS, fpPO->do_kc, fpPO->do_wst, 
                            fpPO->do_wsub, NULL);
    }
    /***
    *   Invert qualification?
    */
    if(fpPO->do_not) {
        ok = !ok;
    }
    return(ok);
}
/*************************************************************************
*   Check option flags and decide if this line should be output or not
*/
int FpOutputThisLineI(FILEPICK *fpPO,int record,int in_rec,int ok_rec)
{
    int ok;
 
    /***
    *   Header = before first record
    */
    if(record < 1) {
        ok = (fpPO->do_dh) ? TRUE : FALSE;
    }
    else {
        ok = (in_rec && ok_rec) ? TRUE : FALSE;
    }
    return(ok);
}
/*************************************************************************
*   Make sure output file is set right
*/
int HandleFilepickOutfileI(FILEPICK *fpPO,int record,char *recS)
{
    char baseS[NSIZE],extS[NSIZE];

    /***
    *   Separate outputs for each file, construct name after checking for old
    */
    if(fpPO->do_esf) {
        CHECK_NFILE(fpPO->out,fpPO->outname);
        GetFilePartsI(fpPO->inname,NULL,baseS,extS);
        if(!NO_S(fpPO->esbase)) {
            strcpy(baseS,fpPO->esbase);
        }
        if(!NO_S(fpPO->esext)) {
            strcpy(extS,fpPO->esext);
        }
        sprintf(fpPO->outname,"%s_part%02d.%s",baseS,record,extS);
    }
    /***
    *   Try to open?
    */
    if(!NO_S(fpPO->outname)) {
        if(!(fpPO->out = OpenUFilePF(fpPO->outname,"w",NULL))) {
            return(FALSE);
        }
    }
    return(TRUE);
}
/*************************************************************************
*   Handle output line
*/
int FpOutOneLineI(FILEPICK *fpPO, char *lineS, FILE *outPF)
{
    HAND_NFILE(outPF);
    fprintf(outPF,"%s",lineS);
    return(TRUE);
}
