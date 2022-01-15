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
 * Program: srt2lcp
 *
 * Description:
 *
 *	Reads Sorted Suffix Array of Genome and creates Rank Array and 
 *       Least Common Prefix (LCP) arrays. These two arrays are written
 *       out to two different files. 
 *
 *       The format of the input file is as follows:
 *                  1 unsigned int ;(length of Genome string, say N)
 *                  N unsigned int ;indices of suffixes after sorting.
 *                  N char         ; Genome-string
 *                  
 *       The format of the output files is as follows:
 *       
 *       Rank Array: 1 unsigned int ; (length of the Rank Array =  N);
 *                   N unsigned int ; (Ranks)
 *
 *       LCP Array:  1 unsigned int ; (length of the LCP Array = N);
 *                   N unsigned int ; (LCPs)
 *
 * LCP and Rank arrays are computed using the linear-time algorithm 
 * given in [1].
 *           
 * Reference:
 *	[1]. Kasai, T., et al. 
 *            Linear-Time Longest-Common-Prefix computation in Suffix Arrays 
 *            and Its applications,
 *            Annual Symposium on Combinatorial Pattern Matching CPM-2001,
 *            Jerusalem, Israel
 *       [2]. M. Ganapathiraju, J. Klein-Seetharaman, et. al.
 *            Comparative n-gram analysis of whole-genome protein sequences,
 *            Proceedings of HLT2002: Human Language Technologies Conference,
 *            San Diego, USA, 2002.
 * 
 * Usage:
 *
 *	./srt2lcp -help
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
#include <time.h>
#include "mylib.h"

#define INFO

void parse_arguments(int argC, char **argV);
void display_usage(char *program);
void display_parameters();

char string[1000], srtfilename[256], lcpfilename[256], rnkfilename[256];

int main (int argc, char **argv){
  unsigned int i,j, superstrlen;
  unsigned int *lcp; 
  char **ptr;
  char *superstring;
  
  FILE *fin=NULL, *fout=NULL, *flcp=NULL, *frnk=NULL;
  
  superstring = (char *) malloc (sizeof(char) * SUPERLEN);
  ptr = (char **) malloc (sizeof(char *) * SUPERLEN);
  
  parse_arguments(argc, argv);
  if (!strlen(srtfilename)){
    printf("Enter Sorted-Suffix-Array of Genome filename:");
    fflush(stdout);
    while (!scanf("%s",srtfilename));
    sprintf(lcpfilename,"%s.lcp",srtfilename);
    sprintf(rnkfilename,"%s.rnk",srtfilename);
  }

#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  printf ("INFO: Input (Suffix-Array-of-Genome) file: %s\n", srtfilename);
  printf ("INFO: Output        (LCP)            file: %s\n", lcpfilename);
  printf ("INFO: Output        (Rank)           file: %s\n", rnkfilename);
#endif

  if ((fin = fopen(srtfilename,"rb")) == NULL){
    sprintf(string, "ERROR: Could not open input file %s for reading;",srtfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }
  
  if ((flcp = fopen(lcpfilename,"wb")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",lcpfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }
  
  if ((frnk = fopen(rnkfilename,"wb")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",rnkfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }

#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  printf ("\nINFO: Beginning to read Suffix Array...\n");
#endif
  
  superstring[0] = '\0';
  superstrlen = readFromFile(fin, superstring, ptr);
  
  
#ifdef INFO
  printf("INFO: Read Sorted-Suffix Array from file.\nINFO: Genome Length = %u\n",superstrlen);
#endif
  
  
  lcp = (unsigned int *) malloc (sizeof(unsigned int) * superstrlen);
  
  
#ifdef INFO	
  printf("INFO: Now finding LCP's...\n");
#endif
  
  
  findFastLCPs(superstring, ptr, superstrlen, lcp, flcp, frnk);
  
  if (fin) 
    fclose(fin);
  if (flcp){
    fflush(flcp);
    fclose(flcp);
  }
  if (frnk){
    fflush(frnk);
    fclose(frnk);
  }

#ifdef INFO
  printf ("INFO: Wrote LCP's to file.Done.\n");
  printf ("\n----------------------------------------------------------\n\n");
#endif
  return (0);
}

/*********************************************/
void parse_arguments(int argC, char **argV){
  int i,j,strNum=1, charNum;
  char *Switch[]={"fsrt","flcp","frnk", "help", "##END##"};
  typedef enum {f_srt,f_lcp,f_rnk,help, f_end
  }Params;
  
  srtfilename[0]='\0';
  while(strNum < argC){
    if (argV[strNum][0] == '-'){
      for (i=0; i<f_end;i++)
        if(!strcmp(argV[strNum]+1,Switch[i]))
          break;
      if (i==f_end)  {            /* The word doesnt match switch */
#ifdef INFO
        printf ("INFO: Inapproppriate Switch %s\tIgnored\n",argV[strNum]);
#endif
        strNum++;
        continue;
      }
      switch(i){
      case help:
	display_usage(argV[0]);
	exit(0);
      case f_srt:
        strcpy(srtfilename,argV[strNum+1]);
	sprintf(lcpfilename,"%s.lcp",srtfilename);
	sprintf(rnkfilename,"%s.rnk",srtfilename);
        strNum+=2;
	continue;
      case f_lcp:
        strcpy(lcpfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case f_rnk:
        strcpy(rnkfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      }
    }
    strNum++; // This statement is reached only if the '-' didnt match above
  }
}
/*********************************************/

void display_parameters(){
  printf("Sorted-Suffixes-Genome (Input) File:\t\t%s\n", srtfilename);
  printf("LCP (Output) File:\t\%s\n",lcpfilename);
  printf("Rank (Output) File:\t\%s\n",rnkfilename);
  
}

/*********************************************/
void display_usage(char *program){
  printf("\n\nUsage:\n%s\n", program);
  printf("\t-fsrt    <Sorted-Suffix-Array of Genome (Input) filename  (.faa.srt)>\n");
  printf("\t-flcp    <LCP (Output) filename  (default: .srt.lcp)>\n");
  printf("\t-frnk    <Rank (Output) filename (default: .srt.rnk)>\n");
  printf("\t-help    display this help message\n\n");
}

/*********************************************/













