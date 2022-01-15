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
 * Program: proteinNGram
 *
 * Description:
 *
 *	Reads Sorted Suffix Array and LCP array and finds the longest 
 *      repeating subsequence in the genome.
 *
 * Usage:
 *
 *	./proteinNGram -help
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

char string[1000],faafilename[26], srtfilename[26], lcpfilename[26], ngramsfilename[20];
char statsfilename[256];
unsigned short flag_printheaders=0, flag_nosort=0, ngram_length=4;

/*********************************************/

int main (int argc, char **argv){
  unsigned int superstrlen, no_of_proteins;
  unsigned int pposition[MAXPROTEINS], plengths[MAXPROTEINS], indices[MAXPROTEINS];
  unsigned int i,j,k, n;
  char **ptr;
  char *superstring, **pheader, *str, protein[MAXLEN],tmpStr[MAXLEN],printbuffer[MAXLEN];
  FILE *fsrt=NULL, *flcp=NULL, *fout=NULL, *fngrams=NULL, *fstats=NULL;

  unsigned int no_of_ngrams, ng;
  unsigned int *lcp=NULL;
  char longstr[MAXLEN],tmp[MAXLEN];
  int no_of_matched_proteins;

  superstring = (char *) malloc (sizeof(char) * SUPERLEN);
  ptr = (char **) malloc (sizeof(char *) * SUPERLEN);
  pheader = (char **) malloc (sizeof(char *) * MAXPROTEINS);
  str = (char *) malloc (sizeof(char) * MAXLEN);
  parse_arguments(argc, argv);
  if (!strlen(srtfilename) || !(strlen(lcpfilename)) || !(strlen(statsfilename))){
    display_usage(argv[0]);
    return (FAILED);
  }
#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  display_parameters();
  printf ("..........................................................\n");
#endif

  if ((fsrt = fopen(srtfilename,"r")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",srtfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }

  if ((flcp = fopen(lcpfilename,"r")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",lcpfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }

  if ((fstats = fopen(statsfilename,"w")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",statsfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }

  if (strlen(ngramsfilename)){
    if ((fngrams = fopen(ngramsfilename,"r")) == NULL){
      sprintf(string, "ERROR: Could not open input file %s for reading;",ngramsfilename);
      sprintf(string, "INFO: Will take user input\n");
      perror (string);
    }
  }


  if (fngrams != NULL){
    no_of_ngrams = 0;
    do {
      string[0]='\0';
      fgets(string,MAXLEN,fngrams);
      if (strlen(string)){
	if (string[strlen(string)-1] == '\n')
	  string[strlen(string)-1] = '\0';
	strcpy (str,string);
	no_of_ngrams++;
        break;
      }
    } while (!feof (fngrams));
  }
  else {
    printf("Enter the Protein Sequence: ");
    fflush(stdout);
    while (scanf("%s", tmpStr)>=1){
      if (tmpStr[strlen(tmpStr)] == '\n')
	tmpStr[strlen(tmpStr)] = '\0';
      if (!strcmp(tmpStr,"X")  || !strcmp(tmpStr,"x"))
	break;
      else if (strlen(tmpStr)){
	strcpy(str,tmpStr);
	no_of_ngrams++;
	break;
      }
      printf("Enter the Protein Sequence: ");
      fflush(stdout);
    }
  }
  if (no_of_ngrams < 1){
    fprintf(stderr,"ERROR: No Protein Sequence given, exiting\n");
    exit(0);
  }
    
  superstrlen = readFromFile(fsrt, superstring, ptr);
  if (superstrlen <=0){
    printf ("ERROR: Could not Read Sorted-Suffix-Array of Genome\n");
    printf ("INFO: Exiting\n");
    return (FAILED);
  }
  lcp = (unsigned int *) malloc (sizeof(unsigned int) * superstrlen);
  i = readLCPFromFile(flcp, lcp);
  if (superstrlen <= 0 || i <=0 || i!= superstrlen){
    printf ("ERROR: LCP length and Sorted-Suffix-Array length not matching.\n");
    printf ("No of elements read from LCP = %d\n", i);
    printf ("INFO: Exiting\n");
    return (FAILED);
  }
#ifdef INFO
  printf("INFO: Length of Genome: %d\n",superstrlen);
#endif
  strcpy(longstr, str);
#ifdef INFO
  printf("INFO: Computing Global n-gram statistics for the sequence %s\n",longstr);
#endif
  
  no_of_ngrams = strlen(longstr)-ngram_length+1;
  for (i=0; i<strlen(longstr); i++){
    longstr[i] = toupper(longstr[i]);
  }
  for (ng=0; ng<no_of_ngrams; ng++){
    str[0] = '\0';
    strncpy(str, longstr+ng, ngram_length);
    str[ngram_length] = '\0';
    no_of_matched_proteins = 0;
    no_of_matched_proteins = countOccurancesFromLCP(str,  superstrlen, superstring, ptr, lcp, strlen(str));
    fprintf(fstats, "%s %d\n", str, no_of_matched_proteins);
  }
  fprintf(fstats, "\n");
#ifdef INFO
  printf ("INFO: Wrote statistics to output file\n");
  printf ("INFO: Done\n");
  printf ("----------------------------------------------------------\n"); 
#endif
  if (fsrt)
    fclose(fsrt);
  if (flcp)
    fclose(flcp);
  if (fstats)
    fclose(fstats);
  if(fngrams)
    fclose(fngrams);
  return(0);
}

/*********************************************/
void parse_arguments(int argC, char **argV){
 
  int i,j,strNum=1, charNum;
  char *Switch[]={"fsrt","flcp","fprot", "fstats", "n","help", "##END##"};
  typedef enum {f_srt,f_lcp,f_prot, f_stats, f_nlen, help, f_end
  }Params;
  faafilename[0]='\0';
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
      case f_stats:
        strcpy(statsfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case f_prot:
        strcpy(ngramsfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case f_nlen:
	ngram_length = (unsigned short) atoi(argV[strNum+1]);
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
  printf("Suffix Array File           : %s\n", srtfilename);
  printf("LCP Array File              : %s\n", lcpfilename);
  printf("Protein Sequence File       : %s\n", ngramsfilename);
  printf("Statistics (Output) File    : %s\n", statsfilename);
}
/*********************************************/
void display_usage(char *program){
  printf("\n\nUsage:\n%s\n\n", program);
  printf("\t-fsrt        <Sorted Suffix Array (Input)filename>\n");
  printf("\t-flcp        <LCP filename>\n");
  printf("\t-fprot       <Input Protein Sequence (for n-gram analysis)>\n");
  printf("\t-fstats      <Output file to write n-gram statitics>\n");
  printf("\t-n           <n-gram length; Default: 4>\n");
  printf("\t-help        display this help message\n\n");
}

/*********************************************/





