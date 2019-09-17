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
*       w = 0-based index; Negative = backwards from top of list
*       wS = word destination
*/
int GetWordlistWordI(WORDLIST *wlPO, int w, char *wS, int max)
{
    VALIDATE(wlPO,WORDLIST_ID);
    /* Negative word index = backwards from end of list */
    if(w < 0) {
        w = wlPO->n + w;
    }
    /* Out of range */
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
*   Load words from passed string
*   If sepC is supplied, use this as word separator, else white space
*       Does NOT load any empty 'words'
*
*   Returns number of words added; BOGUS if error
*/
int LoadWordlistFromStringI(WORDLIST *wlPO, char *bufS, char sepC)
{
    int n,sep,i;
    char *cPC, wordS[BBUFF_SIZE+1];

    VALIDATE(wlPO,WORDLIST_ID);
    cPC = bufS;
    INIT_S(wordS);
    n = i = 0;
    while(ISLINE(*cPC)) {
        /***
        *   End-of-word separater?
        */
        sep = 0;
        if(sepC) {
            if(*cPC == sepC) {
                sep++;
            }
        }
        else {
            if( !isgraph(INT(*cPC)) ) {
                sep++;
            }
        }
        /***
        *   If sep, keep word and reset, else collect char to word
        */
        if(sep) {
            if(i>0) {
                wordS[i] = '\0';
                if(! AddWordlistWordI(wlPO, n++, wordS, TRUE) ) {
                    return(BOGUS);
                }
            }
            i=0;
        }
        else {
            if(i>=BBUFF_SIZE) {
                return(BOGUS);
            }
            wordS[i++] = *cPC;
        }
        cPC++;
    }
    /* Last one? */
    if(i>0) {
        wordS[i] = '\0';
        if(! AddWordlistWordI(wlPO, n++, wordS, TRUE) ) {
            return(BOGUS);
        }
    }
    return(n);
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
/**************************************************************************
*   STRINGWORDS object
*   Intended to allow easy indexing of 'words' in string parsed in various
*       ways (e.g. different 'word' seperator) *without copying* the strings
*       as with wordlist
*
*   n = max words (size of allocated index arrays); Default if < 1
*   sepC = word-separator (see set function)
*/
STRINGWORDS *CreateStringwordsPO(int n, char sepC)
{
    STRINGWORDS *swPO;

    if(!(swPO=(STRINGWORDS *)ALLOC(1,sizeof(STRINGWORDS)))) {
        return(NULL);
    }
    swPO->ID = STRINGWORDS_ID;
    if(n < 1) {
        n = DEF_SWNUM;
    }
    swPO->n_pos = n;
    swPO->st_pos = (int *)ALLOC(n, sizeof(int));
    swPO->en_pos = (int *)ALLOC(n, sizeof(int));
    if( (!(swPO->st_pos)) || (!(swPO->en_pos)) ) {
        CHECK_STRINGWORDS(swPO);
        return(NULL);
    }
    /* Init string-specific vars and separater char */
    InitStringwords(swPO);
    SetStringwordsSep(swPO, sepC);
    return(swPO);
}
/*************************************************************************
*   Clean up / de-allocate structure
*/
int DestroyStringwordsI(STRINGWORDS *swPO)
{
    VALIDATE(swPO,STRINGWORDS_ID);
    CHECK_FREE(swPO->st_pos);
    CHECK_FREE(swPO->en_pos);
    FREE(swPO);
    return(TRUE);
}
/************************************************************************
*   Set string-specific vars to initial values
*/
void InitStringwords(STRINGWORDS *swPO)
{
    VALIDATE(swPO,STRINGWORDS_ID);
    swPO->string = NULL;
    swPO->slen = 0;
    swPO->num = 0;
}
/************************************************************************
*   Dump stringwords structure; st and en are start and end word[w]
*/
void DumpStringwords(STRINGWORDS *swPO, int st, int en, FILE *outPF)
{
    int i,n,spos,epos;
    char tS[NSIZE];

    VALIDATE(swPO,STRINGWORDS_ID);
    HAND_NFILE(outPF);
    fprintf(outPF,"Stringwords at %p\n",swPO);
    fprintf(outPF,"n_pos: %d\n",swPO->n_pos);
    n = GetStringwordsNumI(swPO);
    if(swPO->string) {
        fprintf(outPF,"string: |%s|\n",swPO->string);
    }
    else {
        fprintf(outPF,"string: <none>\n");
    }
    if(swPO->sep) {
        fprintf(outPF,"sep:    |%c|\n",swPO->sep);
    }
    else {
        fprintf(outPF,"sep:    <whitespace>\n");
    }
    fprintf(outPF,"num:   %d\n",swPO->num);
    st = (st<0) ? 0 : st;
    en = (en<0) ? n : en;
    LIMIT_NUM(st,0,n);
    LIMIT_NUM(en,0,n);
    for(i=st;i<en;i++) {
        GetStringwordsCoordsI(swPO, i, &spos, &epos);
        GetStringwordsWordI(swPO, i, tS, -1);
        fprintf(outPF,"w[%d] = %d %d |%s|\n",i,spos,epos,tS);
    }
}
/************************************************************************
*   Set word-separator char
*   sepC = word-separator char
*       If char given, only that will be used as separator
*       If zero (0), then any whitespace is used
*/
void SetStringwordsSep(STRINGWORDS *swPO, char sepC)
{
    VALIDATE(swPO,STRINGWORDS_ID);
    swPO->sep = sepC;
}
/*************************************************************************
*   Set up stringwords for passed string; i.e. parse string and remember
*   Passed string is bufS
*   Length of string to parse is slen; If < 0, call strlen(bufS)
*
*   Consequtive separators are ignored; Only non-sep = words
*   
*   Returns number of words
*/
int LoadStringwordsI(STRINGWORDS *swPO, char *bufS, int slen)
{
    int i,w,is_sep,inword;

    VALIDATE(swPO,STRINGWORDS_ID);
    InitStringwords(swPO);
    if(slen < 0) {
        slen = strlen(bufS);
    }
    swPO->slen = slen;
    /* Point to string ... though don't know who owns space */
    swPO->string = bufS;
    /* Each char, saving start / end of words */
    inword = FALSE;
    i = w = 0;
    while(i < swPO->slen) {
        /* Sep char? */
        if(swPO->sep) {
            is_sep = (bufS[i] == swPO->sep) ? TRUE : FALSE;
        }
        else {
            is_sep = (isspace(bufS[i])) ? TRUE : FALSE;
        }
        /* Now handle token / word accounting */
        if(is_sep) {
            if(inword) {
                BOG_CHECK(w >= swPO->n_pos);
                swPO->en_pos[w] = i;
                w++;
            }
            inword = FALSE;
        }
        else {
            if(!inword) {
                swPO->st_pos[w] = i;
            }
            inword = TRUE;
        }
        i++;
    }
    /* Final end position */
    if(inword) {
        BOG_CHECK(w >= swPO->n_pos);
        swPO->en_pos[w] = i;
        w++;
    }
    swPO->num = w;
    return(w);
}
/************************************************************************
*   Get number of words
*/
int GetStringwordsNumI(STRINGWORDS *swPO)
{
    VALIDATE(swPO,STRINGWORDS_ID);
    return(swPO->num);
}
/************************************************************************
*   Get (start, end) coords for word[w]; Negative w is from list end
*   Sets values in passed pointers if word is good
*/
int GetStringwordsCoordsI(STRINGWORDS *swPO, int w, int *sPI, int *ePI)
{
    VALIDATE(swPO,STRINGWORDS_ID);
    /* No words == no answers */
    if(swPO->num < 1) {
        return(FALSE);
    }
    /* Negative index = backwards offset */
    w = (w < 0) ? swPO->num + w : w;
    /* Out of bounds == no answers */
    if( (w < 0) || (w >= swPO->num) ) {
        return(FALSE);
    }
    /* Set pointer values only if real */
    if(sPI) {
        *sPI = swPO->st_pos[w];
    }
    if(ePI) {
        *ePI = swPO->en_pos[w];
    }
    return(TRUE);
}
/************************************************************************
*   Get string (i.e. word) for word[w]; Negative w is from list end
*   Copies into passed wordS (up to max)
*/
int GetStringwordsWordI(STRINGWORDS *swPO, int w, char *wordS, int max)
{
    int wlen;

    VALIDATE(swPO,STRINGWORDS_ID);
    INIT_S(wordS);
    /* No words == no answers */
    if(swPO->num < 1) {
        return(FALSE);
    }
    /* Negative index = backwards offset */
    w = (w < 0) ? swPO->num + w + 1 : w;
    /* Out of bounds == no answers */
    if( (w < 0) || (w >= swPO->num) ) {
        return(FALSE);
    }
    /* Should have string to copy from */
    BOG_CHECK(! swPO->string);
    /* Unless max not set (neg), set to smaller of max and word len */
    wlen = swPO->en_pos[w] - swPO->st_pos[w];
    if(max > 0) {
        wlen = MIN_NUM(max, wlen);
    }
/*
printf("www %d wlen = %d (%d %d)\n",w,wlen,swPO->st_pos[w],swPO->en_pos[w]);
*/
    strncpy(wordS, &swPO->string[swPO->st_pos[w]], wlen);
    wordS[wlen] = '\0';
    return(TRUE);
}
