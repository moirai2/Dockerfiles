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
* Program: faa2srt
*
* Description:
*
*	Reads Genome file and creates a Sorted Suffix-Array.
*       Input file is a fastaa format Genome file. Proteins are read into
*       a single character array. A blank space seperates consecutive
*       proteins. 
*
*       inplace-binary sort with a 3-character radix is used for sorting
*       the large genome-suffix array.
*       
* Reference:
*
*      [1].Manber, U. and G. Myers, 
*          Suffix arrays: A New Method for On-Line String Searches. 
*          SIAM Journal on Computing, 1993. 22(5): p. 935-948.
*      [2].M. Ganapathiraju, J. Klein-Seetharaman, et. al.
*          Comparative n-gram analysis of whole-genome protein sequences,
*          Proceedings of HLT2002: Human Language Technologies Conference,
*          San Diego, USA, 2002.
*
* Usage:
*
*	./faa2srt -help
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

char string[1000], faafilename[256], srtfilename[256];

int main (int argc, char **argv){
  unsigned int superstrlen;
  unsigned int Time_start, Time_end, Time; 
  char **ptr;
  char *superstring;
  
  FILE *fin=NULL, *fsrt=NULL;
  superstring = (char *) malloc (sizeof(char) * SUPERLEN);
  ptr = (char **) malloc (sizeof(char *) * SUPERLEN);
  
  parse_arguments(argc, argv);
  if (!strlen(faafilename)){
    printf("Enter Genome filename:");
    fflush(stdout);
    while (!scanf("%s",faafilename));
    sprintf(srtfilename,"%s.srt",faafilename);
  }

#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  display_parameters();
#endif


  if ((fin = fopen(faafilename,"r")) == NULL){
    sprintf(string, "ERROR: Could not open input file %s for reading;",faafilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }
  if ((fsrt = fopen(srtfilename,"wb")) == NULL){
    sprintf(string, "ERROR: Could not open output file %s for writing;",srtfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }
  Time_start = time(NULL);
#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  printf ("\nINFO: Beginning to read Genome String ...\n");
#endif
  superstrlen = fileToString(fin, superstring);
  Time_end = time(NULL);
  Time = Time_end - Time_start;

#ifdef INFO
  printf ("INFO: Read Sorted-Suffix Array from file.\n");
  printf ("INFO: Genome Length = %u\n",superstrlen);
  printf ("INFO: Time taken: %d %s.\n", Time, Time==1 ? "second": "seconds");
  printf ("INFO: Beginning to Sort Suffixes...\n");
#endif

  Time_start = time(NULL);
  inplaceBinaryFasterSort(superstring, superstrlen, ptr);

#ifdef INFO
  printf ("INFO: Completed Sorting Suffix Array.\n");
  Time_end = time(NULL);
  Time = Time_end - Time_start;
  printf ("INFO: Time taken: %d %s.\n", Time, Time==1 ? "second": "seconds");
  printf ("INFO: Writing Suffix Array to file...\n");
#endif

  if (!writeSorted(fsrt, ptr, superstring, superstrlen)){
    fprintf(stderr, "ERROR: Could not write Suffix Array to file;\n");
    fprintf(stderr, "Press Any key to continue...\n");
    getchar();
    exit(1);
  }

#ifdef INFO
  printf ("INFO: Wrote Suffix Array to file.\n");
  printf ("INFO: Done.\n");
#endif
 /* 
  if (fin) 
    fclose(fin);
  if (fsrt){
    fflush(fsrt);
    fclose(fsrt);
  }
*/
  exit(0);
}


/*********************************************/
void parse_arguments(int argC, char **argV){
  int i,j,strNum=1, charNum;
  char *Switch[]={"ffaa","fsrt","help", "##END##"};
  typedef enum {f_faa,f_srt,help, f_end
  }Params;
  
  faafilename[0]='\0';
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
      case f_faa:
        strcpy(faafilename,argV[strNum+1]);
	sprintf(srtfilename,"%s.srt",faafilename);
        strNum+=2;
	continue;
      case f_srt:
        strcpy(srtfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      }
    }
    strNum++; // This statement is reached only if the '-' didnt match above
  }
}
/*********************************************/

void display_parameters(){
  printf("INFO: Genome                  (Input) File:\t\t%s\n",faafilename);
  printf("INFO: Sorted-Suffixes-Genome (Output) File:\t\t%s\n", srtfilename);
}

/*********************************************/
void display_usage(char *program){
  printf("\n\nUsage:\n%s\n", program);
  printf("\t-ffaa    <Genome (Input) filename  (.faa)>\n");
  printf("\t-fsrt    <Sorted-Suffix-Array of Genome (Output) filename  (default: .faa.srt)>\n");
  printf("\t-help    display this help message\n\n");
}

/*********************************************/

