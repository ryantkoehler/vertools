/*
* gen_nums.h
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

#define VERSION_S "Gen_Nums version 1.5"

#define DEF_FORM_W  1
#define DEF_FORM_P  0
#define DEF_SEED    1234

/*********************** ppp ********************
* C function listing generated by gen_prot
* Thu Apr  5 06:58:35 2018
*/
/****************************************************************
* gen_nums.c
*/
int main(int argc, char **argv);
void Gen_numsUse(void);
int Gen_numsI(int argc, char **argv);
int OpenInFileI(char *inS, FILE **inPPF);
int ParseLetterArgI(char *letS, char *fPC, char *lPC);
int GetNextLetterI(char cC, char fC, char lC, char *nPC);

