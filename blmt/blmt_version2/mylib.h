#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#define SUCCESS 1
#define FAILED  0

#define MAXMM 10
#define MAXLEN 80000 /* Increase this if the "input protein" is longer */
#define SUPERLEN 1200000000

#define MAXPROTEINS 35000
#define MAX_PHEADER_LEN 200

#define ALPHABET_SIZE 27
#define MAX_COUNT  300000

#define MAXNGRAMS 700
#define MAX_ARGS 7

extern unsigned int fileToString (FILE *fp, char *superstring);
extern int openFiles(FILE *fin, FILE *fout, FILE *fmsg, int argc, char **argv);
extern void inplaceSort(char *ss, unsigned int ssl, char *ptr[SUPERLEN]);
extern void findLCPs(char *ptr[SUPERLEN], unsigned int ssl, unsigned int *lcp, FILE *fmsg);
extern int writeSorted(FILE *fout, char *ptr[SUPERLEN], char superstring[SUPERLEN], unsigned int ssl);
extern unsigned int readFromFile(FILE *fin, char superstring[SUPERLEN], char *ptr[SUPERLEN]);
extern unsigned int readLCPFromFile(FILE *fin, unsigned int lcp[]);
extern unsigned int readRankFromFile(FILE *fin, unsigned int rnk[]);
extern void setIndex (int index1, int index2, int index3, unsigned int indices1[27], unsigned int indices2[27][27], unsigned int indices3[27][27][27], int position);
extern unsigned int getIndex(unsigned int index1, unsigned int index2, unsigned int index3,
		 unsigned int  indices1[27], unsigned int indices2[27][27], unsigned int indices3[27][27][27]);
extern unsigned int readProteinHeaders (FILE *fp, char **pheader, unsigned int *pposition);
extern void makeindextable(unsigned int superstrlen,char **ptr, unsigned int ***indextable,unsigned int N);
extern int locateProteins(char *str, unsigned int *num, unsigned int superstrlen, char *superstring,char **ptr, unsigned int ***indextable,int N);
extern unsigned int find_position(char *ptr[SUPERLEN],unsigned int l, unsigned int r, char *subs);
extern unsigned int find_position_uint(unsigned int *ptr[], unsigned int l, unsigned int r, unsigned int num);
extern void inplaceBinarySort(char *ss, unsigned int ssl, char *ptr[SUPERLEN]);
extern void inplaceBinaryFasterSort(char *ss, unsigned int ssl, char *ptr[SUPERLEN]);
extern void inplaceBinarySortUInteger (unsigned int length, unsigned int *numarray, unsigned int *indices);
extern void findFastLCPs(char *superstrlen, char *ptr[SUPERLEN], unsigned int ssl, unsigned int *lcp, FILE *flcp, FILE *frank);
extern int countOccurances(char *str, unsigned int superstrlen, char *superstring,char **ptr, int N);
extern int countOccurancesFromLCP(char *str, unsigned int superstrlen, char *superstring,char **ptr, unsigned int *lcp, int N);
extern unsigned int readProteinLengths (char *ss,unsigned int ssl,  unsigned int *pposition);
extern unsigned int find_ngrams(unsigned int ssl,char *ss, unsigned int *lcp,char **ptr,unsigned short N, unsigned int *ngcount, char **ngrams);
extern void makeRadix(char *ss, unsigned int ssl, char **ptr, unsigned int Indices[ALPHABET_SIZE][ALPHABET_SIZE][ALPHABET_SIZE][3]);
extern unsigned int getRadix(char *ss, unsigned int Indices[ALPHABET_SIZE][ALPHABET_SIZE][ALPHABET_SIZE][3]);
