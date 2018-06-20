#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>
#define MAX_FN_LEN (511)
#define MAX_ID_LEN (511)
#define MAX_SEQ_LEN (300000000)
#define MAX_FQ_LEN (2047)

/* Structure for a single sequence, like from a fasta file */
typedef struct chr {
  char id[MAX_ID_LEN+1];
  char seq[MAX_SEQ_LEN];
  size_t len;
} Chr;
typedef struct chr* ChrP;

/* Structure for a read pair, like from paired end sequencing */
typedef struct sqp {
  char fid [ MAX_ID_LEN + 1];
  char rid [ MAX_ID_LEN + 1];
  char fseq[  MAX_FQ_LEN + 1 ];
  char rseq[  MAX_FQ_LEN + 1 ];
  char fqual[ MAX_FQ_LEN + 1 ];
  char rqual[ MAX_FQ_LEN + 1 ];
  size_t flen;
  size_t rlen;
} Sqp;
typedef struct sqp* SQP;
  
ChrP newSeq( void );

/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]);

/* read_next_fasta
   args 1. pointer to file to be read
        2. ChrP to put the sequence
   returns: FALSE (0) if no problem
            non-zero if there is a problem
*/
int read_next_fasta( FILE* genome_file, ChrP chr );
int gz_read_next_fasta( gzFile genome_file, ChrP chr );
int is_gz( const char* fq_fn );
int read_next_fastqs( FILE* ffq, FILE* rfq, SQP fqpair );
int gz_read_next_fastqs( gzFile gzffq, gzFile gzrfq, SQP fqpair );
int read_fastq( FILE* fastq, char id[],
		char seq[], char qual[], size_t* len );
int gzread_fastq( gzFile fastq, char id[],
		  char seq[], char qual[], size_t* len );

