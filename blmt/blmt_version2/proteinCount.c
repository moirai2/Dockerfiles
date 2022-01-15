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
 * Program: proteinCount
 *
 * Description:
 *
 *	Reads Genome file and prints
 *          1. Total number of proteins in the Genome
 *          2. Average length of proteins in the Genome
 *          3. Lengths of proteins (with an option to print in sorted order)
 *          4. Optionally, the header of the protein. 
 *       
 * Usage:
 *
 *	./proteinCount -help
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
 *
 * UPDATES:
 * FEB-09-2003: Madhavi Ganapathiraju: Corrected order of preference
 * in using srtfilename and faafilename when both are given.
 * fixed a bug so that the correct number of proteins are output.
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
char string[1000],faafilename[256], srtfilename[256],protfilename[256];
unsigned short flag_printheaders = 0, flag_nosort=0;

/*********************************************/

int main (int argc, char **argv){
  unsigned int superstrlen, no_of_proteins;
  unsigned int pposition[MAXPROTEINS], plengths[MAXPROTEINS], indices[MAXPROTEINS];
  unsigned int i,j,k, n;
  char **ptr;
  char *superstring, **pheader, *str, protein[MAXLEN],tmpStr[MAXLEN],printbuffer[MAXLEN];
  FILE *fin=NULL, *fsrt=NULL, *fprot=NULL, *flcp=NULL;
  
  superstring = (char *) malloc (sizeof(char) * SUPERLEN);
  ptr = (char **) malloc (sizeof(char *) * SUPERLEN);
  pheader = (char **) malloc (sizeof(char *) * MAXPROTEINS);
  parse_arguments(argc, argv);
  if (strlen (srtfilename)){
#ifdef INFO
    printf("INFO: Srtfilename used. Faafile if mentioned is ignored\n");
#endif
    if ((fin = fopen(srtfilename,"r")) == NULL){
      sprintf(string, "ERROR: Could not open input file %s for reading;",srtfilename);
#ifdef INFO
      printf ("%s\n", string);
#endif
      perror (string);
      return (FAILED);
    }
    else {
      fsrt = fin;
    }
  }
  else if (strlen (faafilename)){
    if ((fin = fopen(faafilename,"r")) == NULL){
      sprintf(string, "ERROR: Could not open input file %s for writing;",faafilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
      perror (string);
      return (FAILED);
    }
  }
  else {
    display_usage(argv[0]);
    return (FAILED);
  }
#ifdef INFO
  printf ("\n----------------------------------------------------------\n");
  display_parameters();
  printf ("..........................................................\n");
#endif

  if ((fprot = fopen(protfilename,"w")) == NULL){
    sprintf(string, "ERROR: Could not open input file %s for reading;\n",protfilename);
#ifdef INFO
    printf ("%s\n", string);
#endif
    perror (string);
    return (FAILED);
  }

  if (strlen(srtfilename)){
    superstrlen = readFromFile(fsrt, superstring, ptr);
    free (ptr);
    no_of_proteins = readProteinLengths(superstring, superstrlen, pposition);
#ifdef INFO
    printf("INFO: Total number of Proteins = %d\n",no_of_proteins);
#endif
    superstrlen=0;
    for (i=0; i< no_of_proteins; i++){
      indices[i] = i;
      superstrlen+=(pposition[i+1]-pposition[i]);
      fprintf (fprot, "%d\n", pposition[i+1]-pposition[i]);
    }
  }
  else {
    no_of_proteins = readProteinHeaders(fin, pheader, pposition);
#ifdef INFO
    printf ("INFO: Completed Reading Protein Headers\n");
    printf ("INFO: Total number of Proteins = %d\n",no_of_proteins);
#endif
    superstrlen=0;
    for (i=0; i< no_of_proteins; i++){
      indices[i] = i;
      plengths[i] = pposition[i+1]-pposition[i]; 
      superstrlen+= plengths[i];
    }
  }
  if (!flag_nosort){
    inplaceBinarySortUInteger (no_of_proteins, plengths, indices);
  }
#ifdef INFO
  printf("INFO: Average Protein Length = %.1f\n", (float)superstrlen/no_of_proteins);
  printf("INFO: Done.\n");
  printf ("----------------------------------------------------------\n\n");
#endif
  for (i=0; i < no_of_proteins ; i++){
    fprintf (fprot, "%d %s\n", plengths[indices[i]], (flag_printheaders?pheader[indices[i]]: ""));
  }
  fflush(fprot);
  if (fin)
    fclose(fin);
  if (fsrt)
    fclose(fsrt);
  fclose(fprot);
  return(0);
}




/*********************************************/
void parse_arguments(int argC, char **argV){
 
  int i,j,strNum=1, charNum;
  char *Switch[MAX_ARGS]={"ffaa","fsrt","fprot", "printall", "nosort", "help", "##END##"};
  typedef enum {f_faa, f_srt,f_prot, f_printheaders, f_nosort, help, f_end
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
      case f_faa:
        strcpy(faafilename,argV[strNum+1]);
	sprintf(protfilename,"%s.prot",faafilename);
        strNum+=2;
	continue;
      case f_srt:
        strcpy(srtfilename,argV[strNum+1]);
	sprintf(protfilename,"%s.prot",faafilename);
        strNum += 2;
        continue;
      case f_prot:
        strcpy(protfilename,argV[strNum+1]);
        strNum += 2;
        continue;
      case f_printheaders:
	flag_printheaders = 1;
	strNum +=1 ;
	continue;
      case f_nosort:
	flag_nosort = 1;
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
  printf("Genome String (Input)  File : %s\n", faafilename);
  printf("Protein Count (Output) File : %s\n",protfilename);
  printf("Sort by Protein Length?     : %s\n", flag_nosort?"No": "Yes");
  printf("Print Headers to file?      : %s\n", flag_printheaders?"Yes": "No");
}
/*********************************************/
void display_usage(char *program){
  printf("\n\nUsage:\n%s\n\n", program);
  printf("\t-ffaa        <Genome (.faa) (Input) filename>\n");
  printf("\t-fsrt        <Sorted Suffix Array (Input)filename>\n");
  printf("\t-fprot       <Protein count (Output) filename>\n");
  printf("\t[-printall]  (default: OFF. Give this flag if you want proteins Headers to be printed\n");
  printf("\t[-nosort]    (default: OFF. Give this flag to NOT sort proteins by length\n");
  printf("\t-help        display this help message\n\n");
}

/*********************************************/





