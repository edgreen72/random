#include <stdlib.h>                                                          
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "kmer.h"
#include "file_io.h"

#define NUM_TESTS (1000000)
#define MAX_COUNTS (511)
typedef struct kmer_data {
  size_t count;
} kd;

void increment_kmer_count( klnP kmer );
void print_hist( KSP kmers );
void search_kmer_tree( ktnP tree_node,
		       size_t depth, size_t* hist );
void add_kmers_from_seq( ChrP seq, KSP kmers );

void help( void ) {
  printf( "fasta-kmer-spectrum -f <fasta file> -k <kmer length>\n" );
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
    }
  }
  
  /* Like Elsa says, "Let it go!" */
  free(seq);
    
  print_hist( kmers );
  
}

void add_kmers_from_seq( const ChrP seq, KSP kmers ) {
  size_t i;
  kd* data;
  klnP kmer;
  
  for( i = 0; i < seq->len - kmers->k; i++ ) {
    kmer = get_canonical_kmer( &seq->seq[i], kmers );
    if ( kmer == NULL ) {
      kmer = add_canonical_kmer( &seq->seq[i], kmers );
      data = (kd*)malloc(sizeof(kd));
      data->count = 1;
      kmer->data = data;
    }
    else {
      increment_kmer_count( kmer );
    }
  }
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

void increment_kmer_count( klnP kmer ) {
  kd* data;
  data = (kd*)kmer->data;
  data->count++;
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

      
 
