/*
* sutil.c
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
#include <time.h>
#include "prim.h"

/***
*   Dos = \
*   unix = /
*/
#define DIR_SEP_CHAR '/'

#define DB_STRNG    if(DB[10])
#define DB_FNAME    if(DB[11])

/******************************************************************** sss
*   Writes the date into the passed string
*/
void FillDateString(char *dateS)
{
    int i;
    time_t tt;

    (void) time(&tt);
    sprintf(dateS,"%s",asctime (localtime (&tt)));
    /***
    *   Sham to kill \n at end of date
    */
    i = 0;
    while(ISLINE(dateS[i]))
    {   
        if(dateS[i] == '\n') {
            break;
        }
        i++;    
    }
    dateS[i] = '\0';
    return;
}
/******************************************************************** sss
*   Writes the number-format date into the passed string
*/
void FillNumDateString(char *dateS)
{
    int y,m,d;
    char tmpS[DEF_BS],wS[DEF_BS];

    FillDateString(tmpS);
    y = m = d = 0;
    /***
    *   Format assumed as:
    *   Tue Oct 21 11:25:59 2003
    */
    sscanf(tmpS,"%*s %s %d %*s %d",wS,&d,&y);
    switch(UPPER(wS[0]))
    {
        case 'J':
            if(UPPER(wS[1])=='A')
                m = 1;
            else if(UPPER(wS[2])=='N')
                m = 6;
            else if(UPPER(wS[2])=='L')
                m = 7;
            break;
        case 'F':   m = 2;  break;
        case 'M':   
            switch(UPPER(wS[2]))
            {
                case 'R':   m = 3;  break;
                case 'Y':   m = 5;  break;
            }
            break;
        case 'A':   
            switch(UPPER(wS[1]))
            {
                case 'P':   m = 4;  break;
                case 'U':   m = 8;  break;
            }
            break;  /* XXX was a bug! */
        case 'S':   m = 9;  break;
        case 'O':   m = 10; break;
        case 'N':   m = 11; break;
        case 'D':   m = 12; break;
    }
    sscanf(tmpS,"%*s %*s %*s %s",wS);
    sprintf(dateS,"%d-%02d-%02dT%s",y,m,d,wS);
    return;
}
/*************************************************************************
*   Print time from start to end; Time strings expected as from 
*       FillDateString() 
*/
int PrintTimeStoryI(char *stimeS,char *etimeS,int ve, char *preS)
{
    int stime,etime,dtime,d,h,m,s;
    char prS[DEF_BS];

    INIT_S(prS);
    if(preS) {  
        sprintf(prS,"%s",preS); 
    }
    stime = GetDateI(stimeS);
    etime = GetDateI(etimeS);
    if( IS_BOG(stime) || IS_BOG(etime) ) {
        return(FALSE);
    }
    dtime = etime - stime;
    d = DaysI(dtime);
    h = HoursI(dtime);
    m = MinutesI(dtime);
    s = SecondsI(dtime);
    if(ve) {
        printf("\n");
        printf("%sRun time report\n",prS);
        printf("%s  Start: %s\n",prS,stimeS);
        printf("%s  End:   %s\n",prS,etimeS);
    }
    printf("%sDelta time: %d dy, %02d:%02d:%02d\n",prS,d,h,m,s); 
    if(ve) {
        printf("\n");
    }
    return(TRUE);
}
/*************************************************************************/
int GetDateI(char *bufS)
{
    int d,h,m,s;

    h = m = s = -1;
    sscanf(bufS,"%*s %*s %d %d:%d:%d",&d,&h,&m,&s);
    if ( (d < 0) || (h < 0) || (m < 0) || (s < 0) ) {
        return(-1);
    }
    return(d*86400 + h*3600 + m*60 + s);
}
/*************************************************************************/
int DaysI(int sec)
{
    return(sec/86400);
}
/*************************************************************************/
int HoursI(int sec)
{
    return((sec%86400) / 3600);
}
/*************************************************************************/
int MinutesI(int sec)
{
    return((sec%3600) / 60);
}
/*************************************************************************/
int SecondsI(int sec)
{
    return(sec%60);
}
/*************************************************************************
*   Fills passed string with '\0' characters to the passed length 
*/
void CleanString(char *stringS,int num)
{   
    int i;

    for(i=0; i<num; i++)
    {   stringS[i] = '\0';  }
    return;
}
/***********************************************************************
*   True if no "graph" chars on a line; it's blank
*/
int BlankStringI(char *lS)
{
    while(ISLINE(*lS))
    {
        if(isgraph(INT(*lS))) {
            return(FALSE);
        }
        lS++;   
    }
    return(TRUE);
}
/***********************************************************************
*   Replaces occurances of sC (in sS) with rC (in rS)
*/
void ReplaceChars(char sC, char *sS, char rC, char *rS)
{ 
    ReplaceSomeChars(sC,sS,rC,rS,strlen(sS)); 
    rS[strlen(sS)] = '\0'; 
    return;
}
/***********************************************************************
*   Replace chars only up to limit
*/
void ReplaceSomeChars(char sC, char *sS, char rC, char *rS, int lim)
{
    int i;

    i = 0;
    while(i < lim)
    {
        if(sS[i] == sC) {
            rS[i] = rC;
        }
        else {
            rS[i] = sS[i];
        }
        i++;
    }
    return;
}
/****************************************************************************
*   removes occurance of all listed chars in listS from source string sS
*   and places the result in rS
*/
void RemoveChars(char *listS, char *sS, char *rS)
{
    /***
    *   First process sS into rS then continue with rS
    */
    RemoveThisCharI(*listS,sS,rS);
    listS++;
    while(*listS != '\0')
    {
        RemoveThisCharI(*listS,rS,rS);
        listS++;
    }
    return;
}
/***********************************************************************
*   Removes occurances of sC in source string sS into replacment string rS
*/
int RemoveThisCharI(char sC, char *sS, char *rS)
{   
    int j;

    j = RemoveSomeCharsI(sC,sS,rS,strlen(sS)); 
    rS[j] = '\0';
    return(j);
}
/****************************************************************************
*   Remove chars (i.e. don't copy them) only up to limit
*/
int RemoveSomeCharsI(char sC, char *sS, char *rS, int lim)
{
    int i,j;

    i = j = 0;
    while(i < lim)
    {
        if(sS[i] != sC) {
            rS[j] = sS[i];
            j++;
        }
        i++;
    }
    return(j);
}
/************************************************************************
*   Pads non-graph/print characters in a string up to len 
*   DOES NOT CAP THE STRING WITH \0 at the end
*/
void PadString(char *sS,char pC,int len)
{
    int i,all;

    all = FALSE;
    for(i=0; i<len; i++)
    {
        if(!ISLINE(sS[i])) {
            all++;
        }
        if((!isgraph(INT(sS[i])))||(all)) {
            sS[i] = pC;
        }
    }
    return;
}
/*************************************************************************
*   Truncates passed string on char
*/
int TruncateStringI(char *bufS, char tC, char *newS)
{
    int j;

    j = 0;
    while(ISLINE(bufS[j]))
    {
        if(bufS[j] == tC) {
            break;
        }
        newS[j] = bufS[j];
        j++;
    }
    newS[j] = '\0';
    return(j);
}
/*************************************************************************
*   Print sourceS up to len to output file outPF
*/
void PrintString(char *sourceS, int len, FILE *outPF)
{
    int i;

    HAND_NFILE(outPF);
    for(i=0;i<len;i++)
    {   
        fputc(sourceS[i],outPF);    
    }
    return;
}
/*************************************************************************
*   Changes all characters in passed string to uppercase
*/
void Upperize(char *stringS)
{
    UpperizeToLen(stringS, strlen(stringS));
    return;
}
/************************************************************************/
void UpperizeToLen(char *stringS, int n)
{
    int i;

    for(i=0;i<n;i++) {
        stringS[i] = toupper(INT(stringS[i]));  
    }
}
/*************************************************************************
*   Changes all characters in passed string to lowercase
*/
void Lowerize(char *stringS)
{
    LowerizeToLen(stringS, strlen(stringS));
}
/************************************************************************/
void LowerizeToLen(char *stringS, int n) 
{
    int i;

    for(i=0;i<n;i++) {
        stringS[i] = tolower(INT(stringS[i]));  
    }
}
/*************************************************************************
*   Count upper / lowercase chars in string
*   Count's upper and lower cases in passed string
*   Returns the number of chars not either
*/
int CountStringCaseI(char *sS, int len, int *lowPI, int *upPI)
{
    int i,nl,nu,nn;

    nl = nu = nn = 0;
    for(i=0;i<len;i++)
    {
        if(isupper(INT(sS[i])))
        {   nu++; }
        else if(islower(INT(sS[i])))
        {   nl++; }
        else
        {   nn++; }
    }
    if(lowPI) {
        *lowPI = nl;
    }
    if(upPI) {
        *upPI = nu;
    }
    return(nn);
}
/*************************************************************************
*   Kills trailing blanks in passed string
*/
void KillTrailStringBlanks(char *sS)
{
    int i;

    if(!sS) {
        return;
    }
    i = strlen(sS);
    if(i==0) {
        return;
    }
    while( (i>=0) && (!isgraph(INT(sS[i]))) ) {
        i--;
    }
    sS[i+1] = '\0';
    return;
}
/*************************************************************************/
void Chomp(char *sS)
{
    int i;

    if(!sS) {
        return;
    }
    i = strlen(sS);
    if(i==0) {
        return;
    }
    while(i>=0) {
        if(sS[i] == '\n') {
            sS[i] = '\0';
            break;
        }
        i--;
    }
    return;
}
/**************************************************************************
*   Removes path from file name; walks on passed string
*/
void StripFilePath(char *fnameS)
{
    char *cPC = NULL;
    int i;

    i = strlen(fnameS);
    while(i > 0)
    {
        if(fnameS[i] == '/') {
            cPC = &fnameS[i+1];
            break;
        }
        i--;
    }
    if(i > 0) {
        sprintf(fnameS,"%s",cPC);
    }
    return;
}
/**************************************************************************
*   Gets components from passed filename: path, base, extension; any of which
*   may be NULL
*/
int GetFilePartsI(char *fnameS, char *pathS, char *baseS, char *extS)
{
    int i,n,nb,ne,b,e,len,plen,blen,elen;

    DB_FNAME DB_PrI(">> GetFilePartsI |%s| %p %p %p\n",fnameS,pathS,baseS,extS);
    if(extS) {
        INIT_S(extS);
    }
    if(pathS) {
        INIT_S(pathS);
    }
    if(baseS) {
        INIT_S(baseS);
    }
    /***
    *   look for last . and separator 
    */
    len = i = b = e = nb = ne = 0;
    while(ISLINE(fnameS[i]))
    {
        if(fnameS[i]=='.') {
            ne++;
            e = i;
        }
        if(fnameS[i]==DIR_SEP_CHAR) {
            nb++;
            b = i;
        }
        i++;
        len++;
    }
    DB_FNAME DB_PrI("+ len=%d, ne=%d e=%d, nb=%d b=%d\n",len,ne,e,nb,b);
    /***
    *   Special case limits 
    */
    if(nb==0) {
        b = -1;
        DB_FNAME DB_PrI("+ no path\n");
    }
    if(ne==0) {
        e = -1;
        DB_FNAME DB_PrI("+ no exten\n");
    }
    if( b==(len-1) ) {
        b = len;
        DB_FNAME DB_PrI("+ sss/ case (slash at end of path)\n");
    }
    if( e==(len-1) ) {
        e = len;
        DB_FNAME DB_PrI("+ sss. case (dot at end of base)\n");
    }
    if( (e<b) && (ne>0) ) {
        b = len;
        DB_FNAME DB_PrI("+ ./ case (only path)\n");
    }
    if( (b==0) && (nb>0) ) {
        b = len;
        DB_FNAME DB_PrI("+ /sss case (only path)\n");
    }
    DB_FNAME DB_PrI("+ After checks e=%d b=%d\n",e,b);
    /***
    *   Now forward, filling the whip
    */
    i = len = plen = blen = elen = 0;
    while(ISLINE(fnameS[i]))
    {
        if( (i<b) && (nb>0) ) {
            if(pathS) {
                pathS[plen++] = fnameS[i];
            }
        }
        else if( (i>e) && (ne>0) ) {
            if(extS) {
                extS[elen++] = fnameS[i];
            }
        }
        else if( (i!=e) && (i!=b) ) {
            if( baseS ) {
                baseS[blen++] = fnameS[i];
            }
        }
        i++;
        len++;
    }
    DB_FNAME DB_PrI("+ lens: %d, p=%d b=%d e=%d\n",len,plen,blen,elen);
    /***
    *   Count fields and null-terminate / clean up substrings
    */
    n = 0;
    if(elen>0) {
        if(extS) {
            extS[elen] = '\0';
        }
        n++;
    }
    if(plen>0) {
        if(pathS) {
            pathS[plen] = '\0';
        }
        n++;
    }
    if(blen>0) {
        if(baseS) {
            baseS[blen] = '\0';
        }
        n++;
    }
    DB_FNAME DB_PrI("<< GetFilePartsI %d\n",n);
    return(n);
}
/***************************************************************************/
int FnameFromPartsI(char *pathS, char *baseS, char *extS, char *fnameS)
{
    int ok;

    INIT_S(fnameS);
    ok = FALSE;
    if(pathS != NULL) {
        if(isgraph(INT(*pathS))) {
            strcat(fnameS,pathS);
            strcat(fnameS,"/");
        }
    }
    if(baseS != NULL) {
        if(isgraph(INT(*baseS))) {
            strcat(fnameS,baseS);
            ok++;
        }
    }
    if(extS != NULL) {
        if(isgraph(INT(extS[0]))) {
            strcat(fnameS,".");
            strcat(fnameS,extS);
        }
    }
    if( (NO_S(fnameS)) || (!ok) )
    {   return(FALSE);  }
    return(TRUE);
}
/***************************************************************************
*   Deliminates a substring bounded by start and end characters sC and eC
*   If subS is not NULL, the delimited part is copied to subS.
*   Substring length is returned 
*/
int DelimSubStringI(char *inS, char sC, char eC, char *subS)
{
    int i,j,lev;

    DB_STRNG DB_PrI(">> DelimSubstringI '%c' '%c' %x input:\n|%s|\n", sC,eC,subS,inS);
    lev = i = j = 0;
    while(ISLINE(inS[i]))
    {
        DB_STRNG DB_PrI("+ %4d |%c| ",i,inS[i]);
        if(inS[i] == sC) {
            lev++;
            if(lev == 1) {
                DB_STRNG DB_PrI("\n");
                i++;    continue;
            }
        }
        if(lev == 0) {
            DB_STRNG DB_PrI("not yet\n");
            i++;    continue;
        }
        if(inS[i] == eC) {
            lev--;
        }
        if(lev < 1) {
            DB_STRNG DB_PrI("end\n");
            break;
        }
        if(subS) {
            subS[j] = inS[i];
        }   
        DB_STRNG DB_PrI("%d j %4d\n",lev,j);
        i++;    j++;
    }
    if(lev != 0) {
        return(-1);
    }
    DB_STRNG DB_PrI("<< DelimSubstringI %d\n",j);
    return(j);
}
/**************************************************************************/
int GetNthWordI(char *bufS,int n,char *wordS)
{
    return(CopyNthWordI(bufS,n,wordS,TOO_BIG));
}
/**************************************************************************/
int CopyNthWordI(char *bufS,int n,char *wordS,int max)
{
    int i;
    char *cPC;
    
    if(wordS)
    {   INIT_S(wordS);  }
    cPC = bufS;
    PASS_BLANK(cPC);
    i = 1;
    while(i<n)
    {
        if(!isgraph(INT(*cPC))) {
            break;
        }
        NEXT_WORD(cPC);
        i++;
    }
    if(!isgraph(INT(*cPC))) {
        return(FALSE);  
    }
    if(wordS) {
        sscanf(cPC,"%s",wordS);
    }
    return(TRUE);
}
/**************************************************************************
*   Take input string and reverse it
*/
void ReverseString(char *inS,int len,char *outS)
{
    int i;
    char uC,dC;

    for(i=0;i<(len/2);i++)
    {
        uC = inS[i];
        dC = inS[len-i-1];
        outS[i] = dC;
        outS[len-i-1] = uC;
    }
    if(len%2) {
        outS[len/2] = inS[len/2];
    }
    return;
}
/***********************************************************************
*   Parse True/False or Yes/No
*   Returns TRUE, FALSE, or BOGUS is unrecognized
*/
int ParseTrueFalseI(char *wordS,int warn)
{
    int a;

    if( (UPPER(wordS[0])=='T') || (UPPER(wordS[0])=='Y') )
    {   a=TRUE;     }
    else if( (UPPER(wordS[0])=='F') || (UPPER(wordS[0])=='N') )
    {   a=FALSE; }
    else
    {
        if(warn) {
            printf("Problem with keyword value (not True/False Yes/No)\n");
        }
        a = BOGUS;
    }
    return(a);
}
/************************************************************************
*   Print True/False into string based on v
*/
void FillTrueFalseString(int v,char *sS)
{
    if(v)
    {   sprintf(sS,"True"); }
    else
    {   sprintf(sS,"False"); }
    return;
}
/*************************************************************************
*   Set name string for optional file 
*/
int FillOptionalNameStringI(char *nameS,char *deS)
{
    int ok;

    ok = FALSE;
    if(NO_S(nameS))
    {   
        if(deS)
        {   sprintf(deS,"None");    }
    }
    else
    {   
        if(deS)
        {   strcpy(deS,nameS);  }
        ok++; 
    }
    return(ok);
}
/*************************************************************************
*   Reports a parsing error Set name string for optional file 
*/
void ReportParseErrorLine(char *bufS, char *funcS, char *storyS)
{
    PROBLINE;
    printf("Problem parsing line");
    if(funcS) {
        printf(" (%s)",funcS);
    }
    printf("\n"); 
    if(storyS) {
        printf("%s\n",storyS);
    } 
    if(bufS) {
        printf("Line: |%s|\n",bufS);
    } 
    return;
}
/*************************************************************************/
int SameStringsI(char *fS, char *sS, int kc) 
{
    return(WordStringMatchI(fS, sS, kc, FALSE, FALSE));
}
/****************************************************************************
*   Report the status of various word-string match modifier options
*/
void ReportWordMatchMods(char *prefS, int kc, int st, int sub, FILE *outPF)
{
    char pS[DEF_BS];

    HAND_NFILE(outPF);
    INIT_S(pS);
    if(prefS) {
        sprintf(pS,"%s",prefS);
    }
    fprintf(outPF,"%sComparing strings / words is...",pS);
    fprintf(outPF,"%sCase ",pS);
    if(kc) {
        fprintf(outPF,"sensitive (i.e. case matters)\n");
    }
    else {
        fprintf(outPF,"insensitive (i.e. case ignored)\n");
    }
    /***
    *   Subset or Start or Full match
    */
    if(sub) {
        fprintf(outPF,"%sSubstring (Word1 can be a substring anywhere in Word2)",pS);
    }
    else if(st) {
        fprintf(outPF,"%sStart (Word1 can be a substring at the start of Word2)",pS);
    }
    else {
        fprintf(outPF,"%sFull match (Word1 must match fully Word2)",pS);
    }
}
/*************************************************************************
*   See if first word "matches" second with various rules:
*   If kc, then keep case, else ignore
*   If st, only start of fS needs to match sS; sS can be longer than fS
*   If sub, fS only must be a substring of sS; sS can be longer than fS
*/
int WordStringMatchI(char *fS, char *sS, int kc, int st, int sub) 
{
    int ok;
    char ufS[BBUFF_SIZE],usS[BBUFF_SIZE];

    /***
    *   Match start, substring or whole thing
    */
    if(st) {
        ok = (kc) ? (!strncmp(fS,sS,strlen(fS))) : (!strncasecmp(fS,sS,strlen(fS))); 
    }
    else if(sub) {
        ok = FALSE;
        if(kc) {
            if(strstr(sS,fS)) {
                ok = TRUE;
            }
        }
        else {
            strncpy(ufS,fS,BBUFF_SIZE-1);
            Upperize(ufS);
            strncpy(usS,sS,BBUFF_SIZE-1);
            Upperize(usS);
            if(strstr(usS,ufS)) {
                ok = TRUE;
            }
        }
        /*
        ok = (kc) ? (strstr(sS,fS)) : (strcasestr(sS,fS));
        */
    }
    else {
        ok = (kc) ? (!strcmp(fS,sS)) : (!strcasecmp(fS,sS));
    }
    return(ok);
}
/************************************************************************/
void FloatFormatString(int w,int p,char *fmtS)
{
    sprintf(fmtS,"%%%d.%dlf",w,p);
    return;
}
/**********************************************************************
*   Figure out (print) precision for number
*   SHAM Tried doing this with math via INT(x) but couldn't get sham working!
*
*   Sets nid to number of digits for integer part
*   Sets nfd to number of digits for fractional part 
*   Sets pw  to print width (i.e. int + frac with possible '-' and '.')
*/
int NumPrecisionI(DOUB vD, int *nidPI, int *nfdPI, int *pwPI, char *pformS)
{
    int i,nid,nfd,fp,pw,frac;
    char bufS[DEF_BS], formS[DEF_BS];

    sprintf(formS,"%%.%df",MAX_PRT_PREC);
    sprintf(bufS,formS,vD);
    i = frac = nid = nfd = fp = 0;
    while(ISLINE(bufS[i])) 
    {
        if(bufS[i] == '.') {
            frac++;
        }
        /***
        *   Digit: If fraction, save last non-zero one
        *          Non-fraction, simply count
        */
        if(isdigit(INT(bufS[i]))) {
            if(frac) { 
                fp++;
                if(bufS[i] != '0') {
                    nfd  = fp;
                }
            }
            else {
                nid++;
            }
        }
        i++;
    }
    /***
    *   Print width will be w + [neg] + [. frac] 
    */
    pw = nid;
    if(vD < 0.0) {
        pw++;
    }
    if(nfd>0) {
        pw++;
        pw += nfd;
    }
    /***
    *   Set what we've got
    */
    if(nidPI) {
        *nidPI = nid;
    }
    if(nfdPI) {
        *nfdPI = nfd;
    }
    if(pwPI) {
        *pwPI = pw;
    }
    if(pformS) {
        FloatFormatString(pw,nfd,pformS);
    }
    return(TRUE);
}
/*************************************************************************/
int DoubArrayPrecisionI(DOUB *valsPD, int n, int *pPI, int *wPI, char *formS)
{
    int i,w,p,wmax,pmax,neg;

    wmax = pmax = neg = 0;
    for(i=0;i<n;i++)
    {
        NumPrecisionI(valsPD[i], &w, &p, NULL, NULL);
        wmax = MAX_NUM(w,wmax);
        pmax = MAX_NUM(p,pmax);
        if(valsPD[i] < 0.0) {
            neg++;
        }
    }
    w = wmax;
    p = pmax;
    /***
    *   Actual print width = int + [neg] + [.frac]
    */
    if(neg) {
        w++;
    }
    if(p) {
        w++;
        w += p;
    }
    /***
    *   Return things
    */
    if(pPI) {
        *pPI = p;
    }
    if(wPI) {
        *wPI = w;
    }
    if(formS) {
        FloatFormatString(w,p,formS);
    }
    return(TRUE);
}
