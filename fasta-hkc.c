#include <stdlib.h>                                                          
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "kmer.h"
#include "file_io.h"

#define NUM_TESTS (1000000)
#define MAX_COUNTS (511)
#define HKC_LEN (1027);
typedef struct kmer_data {
  size_t count;
} kd;

void increment_kmer_count( klnP kmer );
size_t get_kmer_count( klnP kmer );
void print_hist( KSP kmers );
void search_kmer_tree( ktnP tree_node,
		       size_t depth, size_t* hist );
void add_kmers_from_seq( ChrP seq, KSP kmers );
void find_and_write_hkcs( ChrP seq, KSP kmers );

void help( void ) {
  printf( "fasta-hkc -f <fasta file> -k <kmer length>\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;

  int ich, i;
  
  size_t k;
  char* kmer_str;
  KSP kmers;
  klnP kmer;
  kd* data;
  gzFile fp_gz;
  FILE* fp;
  int gzipped     = 0;
  int read_status = 0;
  int total_read  = 0;
  ChrP seq;
  char fn[MAX_FN_LEN+1];
    
  if( argc == 1 ) {
    help();
  }
  while( (ich=getopt( argc, argv, "k:f:" )) != -1 ) {
    switch(ich) {
    case 'k' :
      k = (size_t)atoi( optarg );
      break;
    case 'h' :
      help();
    case 'f' :
      strcpy( fn, optarg );
      break;
    default :
      help();
    }
  }

  fprintf( stderr, "[Initializing data structures]\n" );
  kmer_str = (char*)malloc(sizeof(char) * k);
  kmers = init_KSP( k );
  seq = newSeq();
  
  /* Get a handle on the input file to pass to parser */
  if ( is_gz( fn ) ) {
    gzipped = 1;
    fprintf( stderr, "Opening gzipped file: %s\n", fn );
    fp_gz = gzopen( fn, "r" );
    if ( fp_gz == NULL ) {
      fprintf( stderr,
	       "ERROR: Problem reading compressed fasta file.\n" );
      exit( 1 );
    }
  }
  else {
    fp = fileOpen( fn, "r" );
    fprintf( stderr, "Opening file: %s\n", fn );
    if ( fp == NULL ) {
      fprintf( stderr,
	       "ERROR: Problem reading fasta file.\n" );
      exit( 1 );
    }
  }

  fprintf( stderr, "[Reading fasta sequences:" );
  while( read_status == 0 ) {
    if ( gzipped ) {
      read_status = gz_read_next_fasta( fp_gz, seq );
    }
    else {
      read_status = read_next_fasta( fp, seq );
    }
    if ( read_status == 0 ) {
      /* Got a good sequence in seq. Find the kmers */
      add_kmers_from_seq( seq, kmers );
      total_read++;
      if ( total_read%10 == 0 ) {
	fprintf( stderr, "." );
      }
    }
  }
  fprintf( stderr, "]\n" );

  if ( is_gz(fn) ) {
    gzclose( fp_gz );
  }
  else {
    fclose( fp );
  }

  /* Get a handle on the input file to pass to parser */
  if ( is_gz( fn ) ) {
    gzipped = 1;
    fprintf( stderr, "Opening gzipped file: %s\n", fn );
    fp_gz = gzopen( fn, "r" );
    if ( fp_gz == NULL ) {
      fprintf( stderr,
	       "ERROR: Problem reading compressed fasta file.\n" );
      exit( 1 );
    }
  }
  else {
    fp = fileOpen( fn, "r" );
    fprintf( stderr, "Opening file: %s\n", fn );
    if ( fp == NULL ) {
      fprintf( stderr,
	       "ERROR: Problem reading fasta file.\n" );
      exit( 1 );
    }
  }

  printf( "%lu\t.\n", kmers->k );
  fprintf( stderr, "[Reading fasta sequences for hkcs]\n" );
  read_status = 0;
  while( read_status == 0 ) {
    if ( gzipped ) {
      read_status = gz_read_next_fasta( fp_gz, seq );
    }
    else {
      read_status = read_next_fasta( fp, seq );
    }
    if ( read_status == 0 ) {
      find_and_write_hkcs( seq, kmers );
    }
  }
  /* Like Elsa says, "Let it go!" */
  free(seq);

}

void add_kmers_from_seq( const ChrP seq, KSP kmers ) {
  size_t i;
  kd* data;
  klnP kmer;
  
  for( i = 0; i < seq->len - kmers->k + 1; i++ ) {
    kmer = get_canonical_kmer( &seq->seq[i], kmers );
    if ( kmer == NULL ) { // new kmer or bad kmer
      kmer = add_canonical_kmer( &seq->seq[i], kmers );
      if ( kmer != NULL ) { //bad kmer?
	data = (kd*)malloc(sizeof(kd));
	data->count = 1;
	kmer->data = data;
      }
    }
    else {
      increment_kmer_count( kmer );
    }
  }
}

void find_and_write_hkcs( const ChrP seq, const KSP kmers ) {
  size_t i, hkc_start, hkc_end, cov;
  char* HKC_seq; // place to copy the HKC for printing
  int in_hkc = 0; // boolean - are we in an HKC at this position?
  size_t max_hkc_len = HKC_LEN;
  kd* data;
  klnP kmer;

  /* Initialize HKC_seq. We'll grow it later if necessary */
  HKC_seq = (char*)malloc(sizeof(char) * (max_hkc_len+1));

  for( i = 0; i < seq->len - kmers->k + 1; i++ ) {
    kmer = get_canonical_kmer( &seq->seq[i], kmers );
    cov = get_kmer_count( kmer );
    if ( cov == 1 ) { // this is an HKC position
      if ( in_hkc ) { // continuing HKC already started
	; // keep going
      }
      else { // starting a new HKC
	hkc_start = i;
	in_hkc    = 1;
      }
    }
    else { // not HKC
      if ( in_hkc ) { // must have just ended a HKC;
	// Therefore, i is first position after end of HKC
	hkc_end = i + kmers->k - 1; // set hkc_end to first position
	// outside of kmer, i.e. 0-index open ended coordinate
	// minus one because i was the first position *outside*
	// the HKC
	if ( (hkc_end - hkc_start) > max_hkc_len ) {
	  free( HKC_seq );
	  HKC_seq = (char*)malloc(sizeof(char) *
				  (hkc_end - hkc_start + 1));
	}
	strncpy( HKC_seq, &seq->seq[hkc_start], hkc_end-hkc_start );
	HKC_seq[ hkc_end - hkc_start ] = '\0';
	printf( "%s\t.\t.\t.\t.\t.\n", HKC_seq );
      }
      in_hkc = 0;
    }
  }
  free( HKC_seq );
}

void print_hist( KSP kmers ) {
  size_t i, len;
  size_t* hist;
  hist = (size_t*)malloc(sizeof(size_t)*MAX_COUNTS);
  for( i = 0; i < MAX_COUNTS; i++ ) {
    hist[i] = 0;
  }
  len = 1<<(K_AR_SIZE*2);
  /* Look in each position of the array part */
  for( i = 0; i < len; i++ ) {
    if ( kmers->ka[i] != NULL ) {
      search_kmer_tree( kmers->ka[i],
			(kmers->k - kmers->k_ar_size),
			hist );
    }
  }
  for( i = 0; i < MAX_COUNTS; i++ ) {
    printf( "%lu %lu\n", i, hist[i] );
  }
  free( hist );
}

void search_kmer_tree( ktnP tree_node,
		       size_t depth, size_t* hist ) {
  klnP leaf_node;
  kd* data;
  if ( depth > 1 ) {
    if ( tree_node->Ap != NULL ) {
      search_kmer_tree( tree_node->Ap, depth - 1, hist );
    }
    if ( tree_node->Cp != NULL ) {
      search_kmer_tree( tree_node->Cp, depth - 1, hist );
    }
    if ( tree_node->Gp != NULL ) {
      search_kmer_tree( tree_node->Gp, depth - 1, hist );
    }
    if ( tree_node->Tp != NULL ) {
      search_kmer_tree( tree_node->Tp, depth - 1, hist );
    }
    return;
  }
  else {
    if ( tree_node->Ap != NULL ) {
      leaf_node = tree_node->Ap;
      data = (kd*)leaf_node->data;
      if ( data->count <= MAX_COUNTS ) {
	hist[data->count]++;
      }
    }
    if ( tree_node->Cp != NULL ) {
      leaf_node = tree_node->Cp;
      data = (kd*)leaf_node->data;
      if ( data->count <= MAX_COUNTS ) {
	hist[data->count]++;
      }
    }
    if ( tree_node->Gp != NULL ) {
      leaf_node = tree_node->Gp;
      data = (kd*)leaf_node->data;
      if ( data->count <= MAX_COUNTS ) {
	hist[data->count]++;
      }
    }
    if ( tree_node->Tp != NULL ) {
      leaf_node = tree_node->Tp;
      data = (kd*)leaf_node->data;
      if ( data->count <= MAX_COUNTS ) {
	hist[data->count]++;
      }
    }
  }
}

size_t get_kmer_count( klnP kmer ) {
  kd* data;
  if ( kmer == NULL ) {
    return 0;
  }
  data = (kd*)kmer->data;
  return data->count;
}
 
void increment_kmer_count( klnP kmer ) {
  kd* data;
  data = (kd*)kmer->data;
  if ( data->count < MAX_COUNTS ) {
    data->count++;
  }
}

void make_random_kmer( char* kmer_str, size_t k ) {
  int i, r;
  for( i = 0; i < k; i++ ) {
    r = rand()%4;
    switch(r) {
    case 0:
      kmer_str[i] = 'A';
      break;
    case 1:
      kmer_str[i] = 'C';
      break;
    case 2:
      kmer_str[i] = 'G';
      break;
    case 3:
      kmer_str[i] = 'T';
      break;
    }
  }
}

      
 
