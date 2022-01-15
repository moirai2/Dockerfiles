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
* Collection of functions
*
* Reference:
*
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
*
* Madhavi Ganapathiraju: Corrected readProteinLengths, which was returning 
*   one number greater than the actual number of proteins.
*
***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "mylib.h"

#define INFO
/***********************************************************
 * int fileToString (FILE **fp, char *superstring)
 * reads the Gene protein sequences from a fastaa format file 
 * into a char array. The char array is limited
 * in size by SUPERLEN. The header lines of each
 * protein sequence, that start with a '>' are
 * ignored. The end of each sequence is preserved
 * by inserting a SPACE in the superstring
 *
 * Arguments: FILE **fp: input-file pointer
 *			 char *superstring: array into which 
 *                                           string is read
 * Return value: length of superstring
 ************************************************************/
unsigned int fileToString (FILE *fp, char *superstring){
  char string[MAXLEN],buffer[500000]; 
  superstring[0]='\0';
  buffer[0]='\0';
  unsigned int i=0,l=0;
  do {
    string[0]='\0';
    fgets(string,MAXLEN,fp);
    if (string[0] == '>'){
      strcat (buffer+l, " ");
      l++;
      continue;
    }
    i=strlen(string);
    if (i){
      if (string[i-1] == '\n')
	string[--i] = '\0';
      strcat (buffer+l,string);
    }
    l+=i;
    if (l> 500000-MAXLEN-1){
      strcat(superstring, buffer);
      buffer[0]='\0';
      l=0;
      printf("."); fflush(stdout);
    }
  } while (!feof (fp));
  fclose (fp);
  strcat (buffer, " ");
  strcat(superstring, buffer);
  
  return (strlen(superstring));
}

/****************************************************************
 * read protein headers
 * This function may be called after seperately from 
 * reading the genesequence itself.. It is possible to
 * integrate ths function readProteinHeaders with fileToString 
 * so that the file read is done only once.
 * An additional field ENDOFPROTEINS is appended in string 
 * 	position following the last protein. 
 * Arguments: File *fp: pointer to opened input file
 *   char **pheader: pointer to an array of strings
 *			into which the proteins headers are readin
 *   unsigned int *pposition: Pointer to an integer array, into which
 *			the positions of the beginning of the proteins
 *			are read in.
 * Return value: Integer giving the number of proteins read in.
 *****************************************************************/

unsigned int readProteinHeaders (FILE *fp, char **pheader, unsigned int *pposition){
  char string[MAXLEN];
  unsigned int num=0, curr_position = 0;
  do {
    string[0]='\0';
    fgets(string,MAXLEN,fp);
    if (strlen(string)){
      if (string[strlen(string)-1] == '\n')
	string[strlen(string)-1] = '\0';
      if (isalpha(string[0])){
	curr_position += strlen(string);
	continue;
      }
      else if (string[0] == '>'){
	curr_position += 1;
	pheader[num] = (char *) malloc (sizeof(char) * (strlen(string) + 1));
	strcpy(pheader[num], string);
	pposition[num] = curr_position;
	num++;
	continue;
		  }
    }
  } while (!feof (fp));
  pposition[num] = curr_position;
  pheader[num] = (char *) malloc (sizeof(char) * strlen("ENDOFPROTEINS"));
  strcpy(pheader[num],"ENDOFPROTEINS");
  fclose (fp);
  return(num);
}
/*******************************************************************/

unsigned int readProteinLengths (char *ss,unsigned int ssl,  unsigned int *pposition){
  unsigned int i, pctr;
  pposition[0] = 0;
  for (i=1, pctr=1; i < ssl; i++){
    if (ss[i] == ' '){
      pposition[pctr] = i+1;
      pctr++;
    }
  }
  return (pctr);
}
/*********************************************************************
 * makeindextable 
 * indextable is an optional three letter index that may be created
 * for faster access across the genome
 **********************************************************************/

void makeindextable(unsigned int superstrlen,char **ptr, unsigned int ***indextable,unsigned int N){
  unsigned int i,j,k,n;
  char *substr;
  substr = (char *) malloc (sizeof(char) * (N+1));
  strcpy(substr,"   ");
  for (n=0; n<superstrlen; n++){
    if ((ptr[n][0] >=  'A' && ptr[n][0] <=  'Z') && 
	(ptr[n][1] >=  'A' && ptr[n][1] <=  'Z') &&
	(ptr[n][2] >=  'A' && ptr[n][2] <=  'Z')){
      if (strncmp(substr,ptr[n],N)){
	strncpy(substr,ptr[n],N);
	indextable[substr[0]-'A'][substr[1]-'A'][substr[2]-'A'] = n;
      }
    }
  }
}
/*******************************************************/

int locateProteins(char *str, unsigned int *num, unsigned int superstrlen, char *superstring,char **ptr, unsigned int ***indextable,int N){
  char tmp[80];
  unsigned int base_index = 0, n;
  int i;
  base_index = find_position(ptr, 0, superstrlen-1, str);
  for (n=base_index; n<superstrlen;n++){
    for (i=0;i<N; i++){
      if (isspace(ptr[n][i])){
	base_index++;
	break;
      }
    }
    if(i==N)
      break;
  }
  for (n=base_index,i=0; i<MAX_COUNT && n<superstrlen;n++){
    num[i] = 0;
    if (strncmp(str, ptr[n],N) < 0)
      break;
    else if (strncmp(str, ptr[n],N) == 0){
      num[i] = ptr[n]-superstring;
      i++;
    }
    else if (strncmp(str, ptr[n],N) < 0)
      continue;
  }
  printf("No of occurances of %s is %d\n", str, i);
  return (i);
}

/*********************************************************************
 * countOccurances: 
 * Arguments:
 * char * str: n-gram to be searched and counted
 * unsigned int superstrlen: Length of the genome
 * char *superstring: Genome
 * char **ptr: Suffix Array (pointers in memory to the sorted suffixes)
 * int N: Length of the n-gram (the n in the n-gram)
 **********************************************************************/
int countOccurances(char *str, unsigned int superstrlen, char *superstring,char **ptr, int N){
  unsigned int base_index = 0, n;
  int i;
  base_index = find_position(ptr, 0, superstrlen-1, str);
  for (n=base_index; n<superstrlen;n++){
    for (i=0;i<N; i++){
      if (isspace(ptr[n][i])){
	base_index++;
	break;
      }
    }
    if(i==N)
      break;
  }
  for (n=base_index,i=0; n<superstrlen;n++){
    if (strncmp(str, ptr[n],N) < 0)
      break;
    else if (strncmp(str, ptr[n],N) == 0)
      i++;
  }
  return (i);
}
/*********************************************************************
 * countOccurancesFromLCP:
 * Arguments:
 * char * str: n-gram to be searched and counted
 * unsigned int superstrlen: Length of the genome
 * char *superstring: Genome
 * char **ptr: Suffix Array (pointers in memory to the sorted suffixes)
 * int N: Length of the n-gram (the n in the n-gram)
 **********************************************************************/
int countOccurancesFromLCP(char *str, unsigned int superstrlen, char *superstring,char **ptr, unsigned int *lcp, int N){
  unsigned int base_index = 0, n;
  int i;
  base_index = find_position(ptr, 0, superstrlen-1, str);
  for (n=base_index; n<superstrlen;n++){
    if (strncmp(str, ptr[n],N) == 0)
      break;
    else if (strncmp(str, ptr[n],N) < 0)
      return (0);
    else base_index++;
  }
  base_index++; // You should start counting from the LCP of the next suffix;
  for (n=base_index,i=1; n<superstrlen;n++){
    if (lcp[n] >= N)
      i++;
    else break;
  }
  return (i);
}
/**********************************************************************
* findLCPs(char *ptr[SUPERLEN], unsigned int ssl, unsigned int *lcp);
* finds Least Common Prefix lengths in consecutive suffix
* arrays,  for 0,1,..., N number of mismatches where N is 
* MAXMM defined in the headers.
* Arguments: ptr: ptrs to suffix arrays
*			 ssl: Total number of suffix arrays (length of superstring)
*			 lcp: NxM int integer array that holds length of prefix 
*					matches. N is the number of suffix arrays and M is 
*					the number of mismatches allowed plus 1
* Returns: None
************************************************************************/
void findLCPs(char *ptr[SUPERLEN], unsigned int ssl, unsigned int *lcp, FILE *fmsg){
	unsigned int i,j;
	
	char *p1, *p2;
	for (i=1; i < ssl; i++){
	  p1 = ptr[i-1];
	  p2 = ptr[i];
	  
	  for (j=0; i+j < ssl; j++){
	    if (p1[j] == ' ' ||  p2[j] == ' ' || p1[j] != p2[j]){
	      lcp[i] = j;
	      break;  
	    }
	  }
	}
	printf("Done with LCP's\n");
	if (fwrite(lcp, sizeof(unsigned int), ssl, fmsg) < ssl){
	  printf("ERROR in Witing LCP\n");
	}
	fflush(fmsg);
}

/***********************************************************
 * int writeSorted(fout, ptr, superstring, ssl));
 * write the positions of the sorted suffix arrays as they 
 * were in the original string. For example, for a word APPLE,
 * the output that will be written out will be 0,4,3,2,1
 * Output format is:
 * 1 unsigned int :  Length of the char array including the 
 spaces that seperate proteins from one another
 * N unsigned ints that give the position of consecutive
 suffix arrays after sorting, as they were in the 
 original char array (as explained above);
 N is the length of the char array.
 * N bytes of original char array
 * Everything is written in binary format.
 * 
 * Return value: SUCCEED on success, FAILED on failure
 ***********************************************************/


int writeSorted(FILE *fout, char *ptr[SUPERLEN], char superstring[SUPERLEN], unsigned int ssl){
  unsigned int i,l,k, disp;
  if (!fwrite (&ssl, sizeof(unsigned int), 1, fout)){
    perror("Could not write ssl");
    return (FAILED);
  }
  
  printf("Wrote the strlen %d to file\n",ssl);fflush(stdout);
  for (i=0;i<ssl;i++){
   disp = ptr[i]-superstring;
   if ((k = fwrite ((void *)&disp, sizeof(unsigned int), 1, fout)) < 1){
        fprintf(stderr,"Wrote only %u elements of disp\n",i);
	perror("Could not write disp");
	return (FAILED);
   }
  }
  printf("Wrote PTRs\n");fflush(stdout);
  superstring[ssl]='\0';
  if ((k = fwrite (superstring, sizeof(char), ssl, fout)) < ssl){
      fprintf(stderr,"Wrote only %u elements of Genome-String\n",k);
      perror("Could not write superstring");
      return (FAILED);
  }
  printf("Wrote Genome-String to Srt file\n");
  fflush(fout);
  fclose (fout);
  return (SUCCESS);
}


/**************************************************************/
unsigned int readFromFile(FILE *fin, char superstring[SUPERLEN], char *ptr[SUPERLEN]){
  unsigned int i, j,ssl, disp;
  if (!fin){
#ifdef INFO
    printf ("INFO: readFromFile:  Null File pointer. Cant read from file\n");
#endif
    fprintf(stderr,"INFO: readFromFile: Null File pointer. Cant read from file\n");
    return (FAILED);
  }

  if (!fread(&ssl, sizeof(unsigned int), 1, fin)){
#ifdef INFO
    printf ("ERROR: could not read superstrlen from file");
#endif
    perror ("ERROR: could not read superstrlen from file");
    return (FAILED);
  }
  for(j=0; j<ssl; j++){ 
    if ((i = fread((void *)&disp, sizeof(unsigned int), 1, fin)) < 1){
      fprintf(stderr,"Read only %u elements of disp\n", i);
#ifdef INFO
      printf ("ERROR: Could not read pointers to sorted suffexes from input file\n");
#endif
      perror ("ERROR: could not read displacements of suffix arrays");
      return (FAILED);
    }
    ptr[j] = superstring + disp;
  }
#ifdef INFO
  printf ("INFO: Read Pointers to sorted suffixes.\n");
#endif
  if ((i = fread((void *)superstring, sizeof(char), ssl, fin)) < ssl) {
    
#ifdef INFO
    printf ("ERROR: Could not read Genome-String from input file\n");
#endif
    perror("ERROR: could not read chars of suffix arrays");
    return (FAILED);
  }
#ifdef INFO
  printf ("INFO: Finished reading the entire suffix array from file.\n");
#endif 
  
  fclose (fin);
  return (ssl);
}
/**************************************************************/
unsigned int readLCPFromFile(FILE *fin, unsigned int lcp[]){
  unsigned int i, ssl;
  if (!fin){
#ifdef INFO
    printf ("INFO: Null File pointer. Cant read from file\n");
#endif
    fprintf(stderr,"INFO: Null File pointer. Cant read from file\n");
    return (FAILED);
  }

  if (!fread(&ssl, sizeof(unsigned int), 1, fin)){
#ifdef INFO
    printf ("ERROR: could not read superstrlen from file");
#endif
    perror ("ERROR: could not read superstrlen from file");
    return (FAILED);
  }
  if ((i = fread((void *)lcp, sizeof(unsigned int), ssl, fin)) < ssl){
    fprintf(stderr,"ERROR: Read only %u elements of lcp\n", i);
#ifdef INFO
      printf ("ERROR: Could not read LCP from file\n");
#endif
      perror("");
      return (FAILED);
  }
#ifdef INFO
  printf ("INFO: Read LCP Array from file.\n");
#endif
  fclose (fin);
  return (ssl);
}

/**************************************************************/
unsigned int readRankFromFile(FILE *fin, unsigned int rnk[]){
  unsigned int i, ssl;
  if (!fin){
#ifdef INFO
    printf ("INFO: Null File pointer. Cant read from file\n");
#endif
    fprintf(stderr,"INFO: Null File pointer. Cant read from file\n");
    return (FAILED);
  }

  if (!fread(&ssl, sizeof(unsigned int), 1, fin)){
#ifdef INFO
    printf ("ERROR: could not read superstrlen from file");
#endif
    perror ("ERROR: could not read superstrlen from file");
    return (FAILED);
  }
  if ((i = fread((void *)rnk, sizeof(unsigned int), ssl, fin)) < ssl){
    fprintf(stderr,"ERROR: Read only %u elements of rnk\n", i);
#ifdef INFO
      printf ("ERROR: Could not read Rank from file\n");
#endif
      perror("");
      return (FAILED);
  }
#ifdef INFO
  printf ("INFO: Read Rank Array from file.\n");
#endif
  fclose (fin);
  return (ssl);
}

/**********************************************************************
 * indices[i], i=0,1,...,26 stores the position at which a suffix arrays 
 * beginning with a particular letter begin, during sorting. This speeds
 * up the sorting speed since if a suffix array beginning N is encountered
 * its position can be searched beginning at indices['N'-'A'+1]. indices[0]
 * corresponds to the suffix arrays that begin with a space
 * Arguments: index--- the position of the first letter of the suffix array
 *  in the alphabet. A=1, B=2 and so on. ' ' = 0
 * indices: array that gives the positions of the suffix arrays
 * as described above
 * position: the position at which the current first char occured
 * Returns: None
 ************************************************************************/
void setIndex (int index1, int index2, int index3, unsigned int indices1[27], unsigned int indices2[27][27], unsigned int indices3[27][27][27], int position){
  int i,j,k;
  if (indices1[index1] == 0)
    indices1[index1] = position;
  if (indices2[index1][index2] == 0)
    indices2[index1][index2] = position;
  if (indices3[index1][index2][index3] == 0)
    indices3[index1][index2][index3] = position;
  
  for (i=index1+1; i < 27; i++){
		if (indices1[i] > 0)
		  indices1[i]++;
  }
  for (i = index1; i<27;i++){
    for (j = index2+1; j<27; j++){
      if (indices2[i][j] > 0)
	indices2[i][j]++;
    }
  }
  for (i=index1; i < 27; i++){
    for (j=index2; j<27; j++){
      for (k=index3+1; j<27; j++){
	if (indices3[i][j][k] > 0)
	  (indices3[i][j][k])++;
      }
    }
  }
}
/*****************************************************************************************
 * The three indices give position of a suffix array based on its first char, first two chars and first three
 * chars respectively.
 * Arguments are : first char, second char, third char displacements from 'A', 
 * and the three arrays that give positions based on 1st char, 1st 2 chars, 1st 3chars respectively
 * Return value: Position from which to compare
 *************************************************************************************************************/
unsigned int getIndex(unsigned int index1, unsigned int index2, unsigned int index3,
		       unsigned int  indices1[27], unsigned int indices2[27][27], unsigned int indices3[27][27][27]){
  unsigned int i, j, k, max;
  i = indices1[index1];
  j = indices2[index1][index2];
  k = indices3[index1][index2][index3];
  
  max = (i>j?i:j);
  max = (max > k? max:k);
  return (max);
}


/*********************************************************
 * find_position:
 * binary search method to locate the first suffix beginning with a 
 * substring pattern. 
 * Arguments:
 * char **ptr: Suffix Array (pointers in memory to the sorted suffixes)
 * unsigned int l, r: optional left and right  positions within which
 *	the search is to be restricted. If the enire genome is to be searched
 *	these values should be 0 and superstrlen-1;
 * char *subs: substring
 *
 * Return value: unisgned int, giving the position of the first occurance
 *	of substring in the suffix array 
 **********************************************************/

unsigned int find_position(char *ptr[SUPERLEN],unsigned int l, unsigned int r, char *subs){
  unsigned int m;
  int res;
  m = (l+r)/2;
  while (r-l>1){
    res = strcmp(subs,ptr[m]);
    if(res < 0)
      r = m;
    else
      l = m;
    m = (l+r)/2;
  }
  if (strcmp(subs, ptr[r]) > 0)
    return(r);
  if (strcmp(subs, ptr[l]) > 0 && strcmp(subs,ptr[r]) <=0 )
    return(r);
  else 
    return(l);
}

/****************************************************************
 * inplaceBinaryFasterSort:
 * uses binary sort with a radix based on the first three letters
 ****************************************************************/

void inplaceBinaryFasterSort(char *ss, unsigned int ssl, char *ptr[SUPERLEN]){
  unsigned int i,j, k, left, Indices[ALPHABET_SIZE][ALPHABET_SIZE][ALPHABET_SIZE][3];
  unsigned int subssl;
  unsigned int p;
  char *subsuperstring;
  strcat (ss, "AA");
  makeRadix(ss, ssl, ptr, Indices);
  for (i=1; i<ssl;i++){
    subsuperstring = ptr[i];
    subssl = i;
    left = getRadix(ptr[i], Indices);
    p = find_position(ptr, left, subssl, subsuperstring);
    for (j=subssl; j> p; j--)
      ptr[j] = ptr[j-1];
    ptr[p] = subsuperstring;
  }	
  printf("INFO: Sorted suffixes\n");
}

/***************************************************************************
 * findFastLCPs(char *ptr[SUPERLEN], unsigned int ssl, unsigned int *lcp);
 * finds Least Common Prefix lengths in consecutive suffix
 * arrays,  for 0,1,..., N number of mismatches where N is 
 * MAXMM defined in the headers.
 * Arguments: ptr: ptrs to suffix arrays
 *			 ssl: Total number of suffix arrays (length of superstring)
 *			 lcp: NxM int integer array that holds length of prefix 
 *				matches. N is the number of suffix arrays and M is 
 *				the number of mismatches allowed plus 1
 * Returns: None
 ************************************************************************/
void findFastLCPs(char *superstring, char *ptr[SUPERLEN], unsigned int ssl, unsigned int *lcp, FILE *flcp, FILE *frank){
  unsigned int i,j;
  unsigned int h;
  unsigned int *rank;
  unsigned long Time_start, Time_end, Time;
  rank = (unsigned int *) malloc (sizeof(unsigned int) * ssl);
  printf("Allocated mem to ranks\n");fflush(stdout);
  Time_start = time(NULL);
  for (i=0; i<ssl; i++){
   // printf("%d %d %d %d\n",i,ptr[i],ptr[i]-superstring, ssl);fflush(stdout);
    rank[ptr[i]-superstring] = i;
  }
  //printf("Found mem to ranks\n");fflush(stdout);
  fwrite (&ssl, sizeof(unsigned int), 1, frank);
  if (fwrite(rank, sizeof(unsigned int), ssl, frank) < ssl){
      perror("ERROR in Writing Rank Array\n");
      fflush(stderr);
  }
  //printf("Wrote mem to ranks\n");fflush(stdout);
  fflush (frank);

  h = 0;
  lcp[0] = 0;
  for (i=1; i < ssl; i++){
    if (rank[i] > 1){
      j = ptr[rank[i]-1]  - superstring;
      while (superstring[i+h] == superstring[j+h]){
	h++;
      }
      lcp[rank[i]] = h;
      if (h > 0)
	h--;
    }
  }
  printf("Found LCP to ranks\n");fflush(stdout);
  Time_end = time(NULL);
  Time = Time_end-Time_start;
#ifdef INFO
  printf ("INFO: Calculated LCPs.\n");
  printf ("INFO: Time taken to find LCP's: %d %s.\n", Time, ((Time == 1)?"second":"seconds"));
  fflush (stdout);
#endif
  fwrite (&ssl, sizeof(unsigned int), 1, flcp);
  if (fwrite(lcp, sizeof(unsigned int), ssl, flcp) < ssl){
      perror("ERROR in Witing LCP\n");
      fflush(stderr);
  }
  printf("WRote LCP to ranks\n");fflush(stdout);
  fflush (flcp);
}


/****************************************************************
 * inplaceBinaryFasterSort:
 * uses binary sort with a radix based on the first three letters
 ****************************************************************/

void inplaceBinarySortUInteger (unsigned int length, unsigned int numarray[], unsigned int *indices){
  unsigned int i,j, num, p, *tmp;
  unsigned int **ptr;
  ptr = (unsigned int**) malloc (sizeof(unsigned int *) * length);
  for (i=0; i<length; i++){
    ptr[i] = numarray+i;
  }

  for (i=1; i<length; i++){
    tmp = ptr[i];
    p = find_position_uint(ptr,0,i, *ptr[i]);
    for (j=i; j > p; j--)
      ptr[j] = ptr[j-1];
    ptr[p] = tmp;
  }
  for (i=0; i<length; i++){
    indices[i] = ptr[i]-numarray;
  }
}


/*********************************************************
 * find_position:
 * binary search method to locate the first suffix beginning with a 
 * substring pattern. 
 * Arguments:
 * char **ptr: Suffix Array (pointers in memory to the sorted suffixes)
 * unsigned int l, r: optional left and right  positions within which
 *	the search is to be restricted. If the enire genome is to be searched
 *	these values should be 0 and superstrlen-1;
 * char *subs: substring
 *
 * Return value: unisgned int, giving the position of the first occurance
 *	of substring in the suffix array 
 **********************************************************/

unsigned int find_position_uint(unsigned int *ptr[], unsigned int l, unsigned int r, unsigned int num){
  unsigned int m;
  m = (l+r)/2;
  while (r-l>1){
    if(num < *ptr[m])
      r = m;
    else
      l = m;
    m = (l+r)/2;
  }
  if (num > *ptr[l] > 0 && num <= *ptr[r])
    return(r);
  else 
    return(l);
}


unsigned int find_ngrams(unsigned int ssl,char *ss, unsigned int *lcp,char **ptr,unsigned short N, unsigned int *ngcount, char **ngrams){
  unsigned int i, j, k, p=1, no_of_ngrams=0;
  unsigned short ctr=0;
  char tmp[1000];
  while (lcp[p] >0) p++; // to skip the suffixes starting with a space
  ctr=0;
  j=0;
  while (p < ssl){
    if (lcp[p] < N){
      k=0;
      if (j==N){
	if (ctr ==1){
	  no_of_ngrams--;
	}
	else{
	  ngcount[no_of_ngrams-1] = ctr;
	}
      }
      for (j=0; j< N; j++)
	if (ptr[p][j] == ' ')
	  break;
      if (j == N){
	ngrams[no_of_ngrams] = ptr[p];
	no_of_ngrams++;
      }
      ctr=0;
    }
    ctr++;
    p++;
  }
  if (j==N){
    if (ctr == 1)
      no_of_ngrams--;
    else
      ngcount[no_of_ngrams-1] = ctr;
  }
  return (no_of_ngrams);
  
}

/***********************************************************/
void makeRadix(char *ss, unsigned int ssl, char **ptr, unsigned int Indices[ALPHABET_SIZE][ALPHABET_SIZE][ALPHABET_SIZE][3]){
  unsigned int ind[ALPHABET_SIZE][ALPHABET_SIZE][ALPHABET_SIZE][3], i,j, k,l;
  char c1,c2, c3;
  unsigned int m1, m2, m3;
  for (i=0; i<ALPHABET_SIZE; i++){
    for (j=0; j<ALPHABET_SIZE; j++){
      for (k=0; k<ALPHABET_SIZE; k++){
	Indices[i][j][k][0]=0;
	Indices[i][j][k][1]=0;
	Indices[i][j][k][2]=0;
	ind[i][j][k][0]=0;
	ind[i][j][k][1]=0;
	ind[i][j][k][2]=0;
      }
    }
  }
  printf("Allocated memory\n"); 
  for (i=0; i<ssl; i++){
    c1 = ss[i];
    c2 = ss[i+1];
    c3 = ss[i+2];

    m1 = isupper(c1)? c1-'A'+1:0;
    m2 = isupper(c2)? c2-'A'+1:0;
    m3 = isupper(c3)? c3-'A'+1:0;
    
    (Indices[m1][m2][m3][0])++;
  }
  for(i=0, j=0; i<ALPHABET_SIZE; i++){
    for (k=0; k< ALPHABET_SIZE; k++){
      for (l=0; l<ALPHABET_SIZE; l++){
	Indices[i][k][l][1] = j;
	Indices[i][k][l][2] = j + (Indices[i][k][l][0] ? Indices[i][k][l][0]:1) - 1;
	j += Indices[i][k][l][0];
	//printf ("%d %d %d %d %d\n", i,k,l, Indices[i][k][l][1], Indices[i][k][l][2]);
      }
    }
  }
  printf("Calculated left and right indices positions \n");fflush(stdout);
  for (i=0; i<ssl;i++){
    c1 = ss[i];
    c2 = ss[i+1];
    c3 = ss[i+2];
    
    m1 = isupper(c1)? c1-'A'+1:0;
    m2 = isupper(c2)? c2-'A'+1:0;
    m3 = isupper(c3)? c3-'A'+1:0;
    
    j = Indices[m1][m2][m3][1] + ind[m1][m2][m3][0];
    (ind[m1][m2][m3][0])++;
    ptr[j] = ss+i;
  }
}
  


/***********************************************************/
unsigned int getRadix(char *ss, unsigned int Indices[ALPHABET_SIZE][ALPHABET_SIZE][ALPHABET_SIZE][3]){
  char c1,c2, c3;
  unsigned int m1, m2, m3;
  
  c1 = ss[0];
  c2 = ss[1];
  c3 = ss[2];
  
  m1 = isupper(c1)? c1-'A'+1:0;
  m2 = isupper(c2)? c2-'A'+1:0;
  m3 = isupper(c3)? c3-'A'+1:0;
  
  return(Indices[m1][m2][m3][1]);
  
}
  


