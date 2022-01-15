/* ====================================================================
 * Copyright (c) 2002 Carnegie Mellon University. All rights reserved 
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The names "BLMT" and "Carnegie Mellon" must not be used to
 *    endorse or promote products derived from this software without
 *    prior written permission. To obtain permission, contact 
 *    blmt@cs.cmu.edu.
 *
 * 4. Products derived from this software may not be called "BLMT"
 *    nor may "BLMT" appear in their names without prior written
 *    permission of Carnegie Mellon University. To obtain permission,
 *    contact blmt@cs.cmu.edu.
 *
 * 5. Redistributions of any form whatsoever must retain the following
 *    acknowledgment:
 *    "This product includes software developed by Carnegie
 *    Mellon University (http://www.cs.cmu.edu/~blmt)."
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 */

/**************************************************************************
 * Program: ngrams
 *
 * Description:
 *
 *	Reads Sorted Suffix Array and LCP array and finds the longest 
 *      repeating subsequence in the genome.
 *
 * Usage:
 *
 *	./ngrams -help
 *
 * References:
 *      [1].M. Ganapathiraju, J. Klein-Seetharaman, et. al.
 *          Comparative n-gram analysis of whole-genome protein sequences,
 *          Proceedings of HLT2002: Human Language Technologies Conference,
 *          San Diego, USA, 2002.
 *
 * Date: FEB-26-2002
 *
 * Author: Madhavi K. Ganapathiraju (madhavi@cs.cmu.edu)
 *
 * For Further Information Contact: 
 *
 *       Judith Klein-Seetharaman (judithks@cs.cmu.edu) or
 *       Madhavi K. Ganapathiraju (madhavi@cs.cmu.edu)
 *
 ***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <malloc.h>
#include "mylib.h"

#define INFO

void parse_arguments(int argC, char **argV);
void display_parameters();
void display_usage(char *program);
char string[1000],srtfilename[256], lcpfilename[256], ngramsfilename[256];
unsigned short flag_printngram = 0, flag_sortc=0, ngram_length=0, topn=0;

/*********************************************/

int main (int argc, char **argv){
  unsigned int superstrlen, no_of_ngrams;
  unsigned int *ngcounts, *indices;
  char **ngrams;
  unsigned int i,j,k, n;
  char **ptr;
  char *superstring, **pheader, *str, protein[MAXLEN],tmpStr[MAXLEN],printbuffer[MAXLEN];
  FILE *fsrt=NULL, *flcp=NULL, *fout=NULL, *fngrams=NULL;
  unsigned int *lcp=NULL;
  
  superstring = (char *) malloc (sizeof(char) * SUPERLEN);
  ptr = (char **) malloc (sizeof(char *) * SUPERLEN);
  pheader = (char **) malloc (sizeof(char *) * MAXPROTEINS);
  str = (char *) malloc (sizeof(char) * 400);
  parse_arguments(argc, argv);
  if (!strlen(srtfilename)){
    printf("Enter Sorted-Suffix-Array filename:");
    fflush(stdout);
    while (!scanf("%s",srtfilename));
  }
  if (!strlen(lcpfilename)){
    printf("Enter LCP filename:");
    fflush(stdout);
    while (!scanf("%s",lcpfilename));
  }
  if (!ngram_length){
    printf("Enter N-gram length:");
    fflush(stdout);
    while (!scanf("%d",&ngram_length));
  } 
  if (!strlen(ngramsfilename)){
    sprintf(ngramsfilename, "%s.%dgrams", srtfilename, ngram_length);
  }
  if (topn)
    flag_sortc = 1;
#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  display_parameters();
  printf ("..........................................................\n");
#endif

  if ((fsrt = fopen(srtfilename,"rb")) == NULL){
    sprintf(string, "ERROR: Could not open input file %s for writing;",srtfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }

  if ((flcp = fopen(lcpfilename,"rb")) == NULL){
    sprintf(string, "ERROR: Could not open input file %s for writing;",lcpfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }
  if (ngram_length <= 0){
    printf ("ERROR: N-gram length not positive. Exiting.\n");
    printf ("INFO: Try with another n [ >0 ]. \n");
    return (FAILED);
  }
  if ((fngrams = fopen(ngramsfilename,"w")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",ngramsfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }
  if ((superstrlen = readFromFile(fsrt, superstring, ptr)) <= 0){
    printf ("ERROR: 0 elements read from Sorted-Suffix-Array. \n");
    printf ("INFO: Exiting.\n");
    return (FAILED);
  }
  lcp = (unsigned int *) malloc (sizeof(unsigned int) * superstrlen);
  i=readLCPFromFile(flcp, lcp);
  ngrams = (char **) malloc (sizeof (char *) * superstrlen);
  ngcounts = (unsigned int *) malloc (sizeof (unsigned int) * superstrlen);
  indices = (unsigned int *) malloc (sizeof (unsigned int) * superstrlen);
  printf ("INFO: Allocated pointers\n");
#ifdef INFO
  printf("INFO: Length of Genome: %d\n",superstrlen);
#endif
  no_of_ngrams = find_ngrams(superstrlen, superstring, lcp, ptr, ngram_length, ngcounts, ngrams);
  if (!topn)
    topn = no_of_ngrams;
  for (i=0; i< superstrlen; i++){
    indices[i] = i;
  }
  if (flag_sortc){
    inplaceBinarySortUInteger (no_of_ngrams, ngcounts, indices);
  }
  printf ("\nINFO: Sorted by Count.\n");
  string[0]='\0';
  for (i=no_of_ngrams-1, j=0; i >= 0 ; i--){
    if (flag_sortc == 0) 
      k = no_of_ngrams-1 - i;  /* if printing alphabetically, dont print in revers */
    else 
      k = i;
    if (j<topn){
      if (flag_printngram){
	strncpy(string, ngrams[indices[k]], ngram_length);
	fprintf (fngrams, "%s  ", string);
      }
      fprintf (fngrams, "%d\n",ngcounts[indices[k]]);
      j++;
    }
    else break;
  }
#ifdef INFO
  printf ("INFO: No of %d-grams = %d\n", ngram_length, no_of_ngrams);
  printf ("INFO: Wrote n-gram counts to file.\n");
  printf ("\n----------------------------------------------------------\n");
#endif 
  if (fsrt)
    fclose(fsrt);
  if (flcp)
    fclose (flcp);
  if (fngrams)
    fclose (fngrams);
  return(0);
}




/*********************************************/
void parse_arguments(int argC, char **argV){
 
  int i,j,strNum=1, charNum;
  char *Switch[]={"fsrt","flcp","fngrams",  "printall", "sortc", "n", "top", "help", "##END##"};
  typedef enum { f_srt,f_lcp, f_ngrams, f_printngram, f_sortc, nlen, f_topn, help, f_end
  }Params;
  srtfilename[0]='\0';
  while(strNum < argC){
    if (argV[strNum][0] == '-'){
      for (i=0; i<f_end;i++)
        if(!strcmp(argV[strNum]+1,Switch[i]))
          break;
      if (i==f_end)  {            /* The word doesnt match switch */
        fprintf(stderr,"\nInapproppriate Switch %s\tIgnored\n",argV[strNum]);
        strNum++;
        continue;
      }
      switch(i){
      case help:
	display_usage(argV[0]);
	exit(0);
      case f_srt:
        strcpy(srtfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case f_lcp:
        strcpy(lcpfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case f_ngrams:
        strcpy(ngramsfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case nlen:
	ngram_length = atoi (argV[strNum+1]);
	strNum += 2;
	continue;
      case f_topn:
	topn = atoi (argV[strNum+1]);
	strNum += 2;
	continue;
      case f_printngram:
	flag_printngram = 1;
	strNum +=1 ;
	continue;
      case f_sortc:
	flag_sortc = 1;
	strNum +=1 ;
	continue;
      }
    }
    strNum++; // This statement is reached only if the '-' didnt match above
  }
}
/*********************************************/

void display_parameters(){
  printf("Input Parameters:\n");
  printf("Sorted-Suffix Array (Input) File  : %s\n", srtfilename); 
  printf("LCP Array           (Input) File  : %s\n", lcpfilename); 
  printf("N-gram             (Output) File  : %s\n", ngramsfilename); 
  printf("n-gram length                     : %d\n", ngram_length);
  if (topn) printf("No of N-grams to print            : %d\n", topn);
  printf("Print ngram in Output file?       : %s\n", flag_printngram?"Yes": "No");
  printf("Sorting order                     : %s\n", flag_sortc?"By count":"By n-gram");
}
/*********************************************/
void display_usage(char *program){
  printf("\n\nUsage:\n%s\n\n", program);
  printf("\t-fsrt        <Sorted Suffix Array (Input)filename>\n");
  printf("\t-flcp        <LCP array file (Input)>\n");
  printf("\t-fngrams     <n-gram count (Output) filename>\n");
  printf("\t-n           <n-gram length>\n");
  printf("\t-top N       <print only top N n-grams>\n");
  printf("\t             (with this option, n-grams are sorted by count)\n");
  printf("\t             (Also, if N is 0 or -top option not givenm all n-grams are printed)\n");
  printf("\t[-printall]  (default: ON. Give this flag if you want n-gram \n");
  printf("\t              to be printed besides the count of n-gram)\n");
  printf("\t[-sortc]     (default: OFF. Give this flag to sort n-grams\n");
  printf("\t              by count instead alphabetically)\n");
  printf("\t-help        display this help message\n\n");
}

/*********************************************/












