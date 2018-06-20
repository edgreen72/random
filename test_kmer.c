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

void make_random_kmer( char* kmer_str, size_t k );
void increment_kmer_count( klnP kmer );
void print_hist( KSP kmers );
void search_kmer_tree( ktnP tree_node,
		       size_t depth, size_t* hist );

void help( void ) {
  printf( "test_kmer -k <kmer length> -c [canonical kmers]\n" );
  printf( "Tests the kmer.o object\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;

  int ich, i;
  int canonical_kmer = 0;
  size_t k;
  char* kmer_str;
  KSP kmers;
  klnP kmer;
  kd* data;
  
  if( argc == 1 ) {
    help();
  }
  while( (ich=getopt( argc, argv, "k:hc" )) != -1 ) {
    switch(ich) {
    case 'k' :
      k = (size_t)atoi( optarg );
      break;
    case 'h' :
      help();
    case 'c' :
      canonical_kmer = 1;
      break;
    default :
      help();
    }
  }

  kmer_str = (char*)malloc(sizeof(char) * k);
  kmers = init_KSP( k );
  
  for( i = 0; i < NUM_TESTS; i++ ) {
    make_random_kmer( kmer_str, k );

    if ( canonical_kmer ) {
      kmer = get_canonical_kmer( kmer_str, kmers );
    }
    else {
      kmer = get_kmer( kmer_str, kmers );
    }
    
    if ( kmer == NULL ) {
      if ( canonical_kmer ) {
	kmer = add_canonical_kmer( kmer_str, kmers );
      }
      else {
	kmer = add_kmer( kmer_str, kmers );
      }
      data = (kd*)malloc(sizeof(kd));
      data->count = 0;
      kmer->data = data;
    }
    increment_kmer_count( kmer );
    
    revcom_kmer( kmer_str, k );
    if ( canonical_kmer ) {
      kmer = get_canonical_kmer( kmer_str, kmers );
    }
    else {
      kmer = get_kmer( kmer_str, kmers );
    }

    if ( kmer == NULL ) {
      if ( canonical_kmer ) {
	kmer = add_canonical_kmer( kmer_str, kmers );
      }
      else {
	kmer = add_kmer( kmer_str, kmers );
      }
      data = (kd*)malloc(sizeof(kd));
      data->count = 0;
      kmer->data = data;
    }
    increment_kmer_count( kmer );
  }
  
  print_hist( kmers );
  
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

      
 
