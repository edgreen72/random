#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#define MAX_K (35) // biggest K we can deal with
#define K_AR_SIZE (12) // default length of the array size of the kmer structure
#define MAX_HKC_PER_K (64) // biggest number of HKC that a k can point to

typedef struct kmer_tree_node {
  void* Ap;
  void* Cp;
  void* Gp;
  void* Tp;
} ktn;
typedef struct kmer_tree_node* ktnP;

/* Structure for leaf nodes in k-mer suffix tree */
typedef struct kmer_leaf_node {
  void* data;
} kln;
typedef struct kmer_leaf_node* klnP;

typedef struct kmers {
  size_t k; // length of kmers
  size_t k_ar_size ; // length of kmer part that we'll handle in the array
                 // and not the tree
  ktnP* ka; // the array part;
} Kmers;
typedef struct kmers* KSP;

/* Function prototypes */
KSP init_KSP( int k );
klnP add_kmer( const char* kmer, KSP ks );
klnP get_kmer( const char* kmer, KSP ks );
int remove_kmer( const char* kmer, KSP ks );
int kmer2inx( const char* kmer,
              const size_t kmer_len,
              size_t* inx );
void revcom_kmer( char* kmer, size_t k );
char revcom_base( char c );
ktnP init_ktn( void );
klnP init_kln( void );
void free_kln_data( klnP kp );
