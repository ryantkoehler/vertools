/*
* wordlist.c
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
#include "wordlist.h"

#define DB_DFU  if(DB[13])

/**************************************************************************/
WORDLIST *CreateWordlistPO(char *fnameS, int wsize)
{
    WORDLIST *wlPO;

    if(!(wlPO=(WORDLIST *)ALLOC(1,sizeof(WORDLIST)))) {
        return(NULL);
    }
    wlPO->ID = WORDLIST_ID;
    InitWordlist(wlPO, wsize);
    if( fnameS && (!NO_S(fnameS)) ) {
        if(!WordlistLoadFromFileI(wlPO, fnameS, wsize)) {
            CHECK_WORDLIST(wlPO);
            return(NULL);
        }
    }
    return(wlPO);
}
/*************************************************************************/
int DestroyWordlistI(WORDLIST *wlPO)
{
    VALIDATE(wlPO,WORDLIST_ID);
    CHECK_FREE(wlPO->words);
    FREE(wlPO);
    return(TRUE);
}
/*************************************************************************/
void InitWordlist(WORDLIST *wlPO, int wsize)
{
    VALIDATE(wlPO,WORDLIST_ID);
    INIT_S(wlPO->fname);
    CHECK_FREE(wlPO->words);
    wlPO->n = 0;
    if(wsize > 0) {
        wlPO->wsize = wsize;
    }
    else {
        wlPO->wsize = DEF_WLWSIZE;
    }
}
/*************************************************************************
*   Check if passed index will fit in allocated space and allocate if needed
*/
int HandleWordlistSpaceI(WORDLIST *wlPO, int i)
{
    int nm;

    DB_DFU DB_PrI(">> HandleWordlistSpaceI %p i=%d\n",wlPO,i);
    VALIDATE(wlPO,WORDLIST_ID);
    if(wlPO->wsize < 1) {
        DB_DFU DB_PrI("<< HandleWordlistSpaceI (wsize < 1) FALSE\n");
        return(FALSE);
    }
    /***
    *   If index won't fit in allocated space, alloc
    */
    if(i >= wlPO->n_words) {
        nm = MAX_NUM(i, (wlPO->n_words + ALLOC_BLOCK));
        if(wlPO->n_words == 0) {
            wlPO->words = (char *)ALLOC(nm * wlPO->wsize,sizeof(char));
            DB_DFU DB_PrI("+ Allocated nm=%d %p\n",nm,wlPO->words);
        }
        else {
            wlPO->words = (char *)REALLOC(wlPO->words, nm * wlPO->wsize,sizeof(char));
            DB_DFU DB_PrI("+ ReAllocated nm=%d %p\n",nm,wlPO->words);
        }
        if(!wlPO->words) {
            PROBLINE;
            printf("HandleWordlistSpaceI Failed allocate %d words, %d\n",nm,wlPO->wsize);
            return(FALSE);
        }
        wlPO->n_words = nm;
    }
    /***
    *   Update "number" == high water mark == max of existing or index
    */
    wlPO->n = MAX_NUM(wlPO->n, i+1);
    DB_DFU DB_PrI("<< HandleWordlistSpaceI n=%d n_words=%d TRUE\n", wlPO->n, wlPO->n_words);
    return(TRUE);
}
/*************************************************************************/
void DumpWordlist(WORDLIST *wlPO, int st, int en, FILE *outPF)
{
    int i,n;
    char tS[NSIZE];

    VALIDATE(wlPO,WORDLIST_ID);
    HAND_NFILE(outPF);
    fprintf(outPF,"Wordlist at %p\n",wlPO);
    fprintf(outPF,"Name  |%s|\n",wlPO->name);
    fprintf(outPF,"Fname |%s|\n",wlPO->fname);
    fprintf(outPF,"%d words, max size %d\n",wlPO->n,wlPO->wsize);
    n = GetWordlistLengthI(wlPO);
    st = (st<0) ? 0 : st;
    en = (en<0) ? n : en;
    LIMIT_NUM(st,0,n);
    LIMIT_NUM(en,0,n);
    for(i=st;i<en;i++) {
        GetWordlistWordI(wlPO, i, tS, -1);
        fprintf(outPF,"w[%d] = |%s|\n",i,tS);
    }
}
/***************************************************************************/
int SetWordlistWsizeI(WORDLIST *wlPO, int wsize)
{
    VALIDATE(wlPO,WORDLIST_ID);
    if(wlPO->words) {
        PROBLINE;
        printf("Cannot set wordlist wordsize if already have words\n");
        return(FALSE);
    }
    wlPO->wsize = wsize;
    return(TRUE);
}
/***************************************************************************/
int SetWordlistLengthI(WORDLIST *wlPO, int len)
{
    VALIDATE(wlPO,WORDLIST_ID);
    if(!HandleWordlistSpaceI(wlPO, len)) {
        return(FALSE);
    }
    wlPO->n = len;
    return(TRUE);
}
/***************************************************************************/
int GetWordlistLengthI(WORDLIST *wlPO)
{
    VALIDATE(wlPO,WORDLIST_ID);
    return(wlPO->n);
}
/***************************************************************************/
int AppendWdordlistWordI(WORDLIST *wlPO, char *wS, int fit)
{
    int i;

    i = GetWordlistLengthI(wlPO);
    return(AddWordlistWordI(wlPO, i+1, wS, fit));
}
/****************************************************************************
*   Add past word into worklist, allocating space if needed
*   If fit is true, the word must fully fit
*/
int AddWordlistWordI(WORDLIST *wlPO, int w, char *wS, int fit)
{
    DB_DFU DB_PrI(">> AddWordlistWordI w=%d |%s|\n",w,wS);
    VALIDATE(wlPO,WORDLIST_ID);
    if(HandleWordlistSpaceI(wlPO, w)) {
        if(SetWordlistWordI(wlPO, w, wS, fit)) {
            DB_DFU DB_PrI("<< AddWordlistWordI TRUE\n");
            return(TRUE);
        }
    }
    DB_DFU DB_PrI("<< AddWordlistWordI FALSE\n");
    return(FALSE);
}
/****************************************************************************
*   Set passed word into wordlist; If fit is true, it must fully fit
*/
int SetWordlistWordI(WORDLIST *wlPO, int w, char *wS, int fit)
{
    int max;

    VALIDATE(wlPO,WORDLIST_ID);
    DB_DFU DB_PrI(">> SetWordlistWordI w=%d |%s|\n",w,wS);
    if( (w < 0) || (w >= wlPO->n) ) {
        DB_DFU DB_PrI("<< SetWordlistWordI n=%d FALSE\n",wlPO->n);
        return(FALSE);
    }
    max = strlen(wS);
    if(max >= wlPO->wsize) {
        if(fit) {
            DB_DFU DB_PrI("<< SetWordlistWordI wont fit %d %d FALSE\n",max,wlPO->wsize);
            return(FALSE);
        }
        max = wlPO->wsize - 1;
    }
    /***
    * Copy and cap
    */
    strncpy(&wlPO->words[w * wlPO->wsize], wS, max);
    wlPO->words[w * wlPO->wsize + max] = '\0';

    DB_DFU DB_PrI("<< SetWordlistWordI TRUE\n");
    return(TRUE);
}
/****************************************************************************
*   Get null-terminated string from word list
*/
int GetWordlistWordI(WORDLIST *wlPO, int w, char *wS, int max)
{
    VALIDATE(wlPO,WORDLIST_ID);
    if( (w < 0) || (w >= wlPO->n) ) {
        return(FALSE);
    }
    if( (max<0) || (max >= wlPO->wsize) ) {
        max = wlPO->wsize - 1;
    }
    strncpy(wS, &wlPO->words[w * wlPO->wsize], max);
    wS[max] = '\0';
    return(TRUE);
}
/***************************************************************************
*   Get stats for word lengths in wordlist
*/
int WordlistLenStatsI(WORDLIST *wlPO, int *minPI, int *maxPI, DOUB *avPD)
{
    int i,min,max,len;
    char wordS[BBUFF_SIZE];
    DOUB avD;

    min = BBUFF_SIZE;
    max = 0;
    avD = 0.0;
    for(i=0;i<wlPO->n;i++) {
        GetWordlistWordI(wlPO, i, wordS, -1);
        len = strlen(wordS);
        min = MIN_NUM(len,min);
        max = MAX_NUM(len,max);
        avD += DNUM(len);
    }
    if(minPI) {
        *minPI = min;
    }
    if(maxPI) {
        *maxPI = max;
    }
    if(avPD && (wlPO->n > 0)) {
        *avPD = avD / DNUM(wlPO->n);
    }
    return(i);
}
/***************************************************************************
*   Figure formatting string from word lengths; So all words "fit" the same
*/
int WordlistAutoFormatStringI(WORDLIST *wlPO, int *wPI, char *formS)
{
    int max;
    char bufS[DEF_BS];

    if(!WordlistLenStatsI(wlPO, NULL, &max, NULL)){
        return(FALSE);
    }
    sprintf(bufS,"%%-%ds",max);
    if(formS) {
        strcpy(formS,bufS);
    }
    if(wPI) {
        *wPI = max;
    }
    return(TRUE);
}
/****************************************************************************
*   Check if word is in wordlist
*   Word is specific (i.e. whole) while list items more general (i.e. partial)
*   Returns TRUE if found, sets which element matches if whichPI is passed 
*/
int WordInWordlistI(WORDLIST *wlPO, char *wordS, int kc, int st, int sub, int *whichPI)
{
    int i,in;
    char tS[NSIZE];

    VALIDATE(wlPO,WORDLIST_ID);
    in = FALSE;
    for(i=0;i<wlPO->n;i++) {
        GetWordlistWordI(wlPO, i, tS, -1);
        /***
        *   Compare list items as queries compared to specific word
        */
        if(WordStringMatchI(tS,wordS,kc,st,sub)) {
            if(whichPI) {
                *whichPI = i;
            }
            in++;
            break;
        }
    }
    return(in);
}
/*************************************************************************
*   Load words from passed filename, using max word size wsize 
*/
int WordlistLoadFromFileI(WORDLIST *wlPO, char *fnameS, int wsize)
{
    FILE *inPF;

    VALIDATE(wlPO,WORDLIST_ID);
    if(!(inPF= OpenUFilePF(fnameS,"r",NULL))) {
        return(FALSE);
    }
    wlPO->n = LoadWordStringFromFileI(inPF, wsize, &wlPO->words);
    FILECLOSE(inPF);
    if (wlPO->n < 1) {
        return(FALSE);
    }
    /***
    *   Keep name and word size 
    */
    strcpy(wlPO->fname,fnameS);
    wlPO->wsize = wsize;
    return(TRUE);
}
/****************************************************************************
*   Loads a string array from file. 
*   First token (up to max size) taken from each non-blank, non-comment line
*   Returns number of items collected
*   ALLOCATES space and returns via passed pointer to string 
*   MOVES FILE AHEAD
*/
int LoadWordStringFromFileI(FILE *fPF, int max, char **arrayPPC)
{
    int c,i,p,w,n,nm;
    char bufS[BBUFF_SIZE],*arrayPC;

    DB_DFU DB_PrI(">> LoadWordStringFromFileI f=%p max=%d\n",fPF,max);
    arrayPC = *arrayPPC = NULL;
    if( (max>=BBUFF_SIZE) || (max<1) ) {
        printf("Impossible word size=%d\n",max);
        ERR("LoadWordStringFromFileI","bad word size");
        return(FALSE);
    }
    /***
    *   Load tokens into buffer then add to growing array if ok
    */
    p = n = nm = w = 0;
    while((c=fgetc(fPF))!=EOF)
    {
        /***
        *   Comment line or up front blank space
        */
        if( (p==0) && (c=='#') ) {
            DB_DFU DB_PrI("+ Eating comment line\n");
            EatOneLineI(fPF);
            continue;
        }
        p++;
        if( (!isgraph(c)) && (w==0) ) {
            continue;
        }
        /***
        *   Not done with word, fill in buffer and continue
        */
        if( isgraph(c) && ((w+1)<max) ) {
            bufS[w++] = c;
            continue;
        }
        /***
        *   Finish line if not done and if didn't get any word, continue
        */
        if(ISLINE(c)) {
            EatOneLineI(fPF);
        }
        p = 0;
        if(w<1) {
            DB_DFU DB_PrI("+ Blank line\n");
            continue;
        }
        /***
        *   Make space for word if needed
        */
        bufS[w]='\0';
        DB_DFU DB_PrI("+ [%d]word |%s|\n",n,bufS);
        if((n+1)>nm)
        {
            if(nm == 0) {
                nm = ALLOC_BLOCK;
                arrayPC = (char *)ALLOC(nm * max,sizeof(char));
                DB_DFU DB_PrI("+ Allocated nm=%d %p\n",nm,arrayPC);
            }
            else {
                nm += ALLOC_BLOCK;
                arrayPC = (char *)REALLOC(arrayPC,nm * max,sizeof(char));
                DB_DFU DB_PrI("+ ReAllocated nm=%d %p\n",nm,arrayPC);
            }
            if(!arrayPC) {
                PROBLINE;
                printf("Failed to allocate to size %d (words %d)\n",nm,max);
                return(BOGUS);
            }
        }
        /***
        *   Copy word/nulls into space for current word
        */
        for(i=0;i<max;i++)
        {
            if(isgraph(INT(bufS[i]))) {
                arrayPC[n*max+i] = bufS[i];
            }
            else {
                arrayPC[n*max+i] = '\0';
            }
        }
        arrayPC[n*max+i] = '\0';
        w = p = 0;
        n++;    
    }
    *arrayPPC = arrayPC;
    DB_DFU DB_PrI("<< LoadWordStringFromFileI %p n=%d\n",arrayPC,n);
    return(n);
}
