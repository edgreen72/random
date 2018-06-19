#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "kmer.h"

KSP init_KSP( int k ) {
  KSP ks;
  size_t i, len;
  len = 1<<(K_AR_SIZE*2); // Length of the array for the array part
  ks = (KSP)malloc(sizeof(Kmers));
  
  ks->k_ar_size = K_AR_SIZE;
  ks->k = k;
  ks->ka = (ktnP*)malloc(sizeof(ktnP) * len);
  for( i = 0; i < len; i++ ) {
    ks->ka[i] = NULL;
  }
  return ks;
}

/* add_kmer
   This function takes a kmer as input and returns the data
   associated with that kmer.
   Args: const char* kmer - pointer to the kmer
         KSP ks - pointer the kmer structure
   Returns: klnP - pointer to the new leave node that was added
            NULL if the kmer is bad or couldn't be added
*/
klnP add_kmer( const char* kmer, KSP ks ) {
  size_t inx;
  size_t i;
  ktnP curr_node;
  ktnP next_node;
  klnP last_node;
  //ktnP curr_node, next_node;
  size_t kmer_pos;
  char next_base;

  /* Find what the index position is for the first K_AR_SIZE
     bases in this kmer */
  if ( kmer2inx( kmer, ks->k_ar_size, &inx ) ) {
  
    /* Start at the index position of the first ks->k_ar_size bases */
    curr_node = (ktnP)ks->ka[inx];
    
    /* If we've never seen that before, then initialize it */
    if ( curr_node == NULL ) {
      curr_node = (ktnP)init_ktn();
      ks->ka[inx] = curr_node;
    }
    
    kmer_pos = ks->k_ar_size;
    
    while( kmer_pos < ks->k - 1 ) {
      next_base = toupper(kmer[kmer_pos]);
      switch( next_base ) {
      case 'A' :
        if ( curr_node->Ap == NULL ) {
          curr_node->Ap =  (ktnP)init_ktn();
        }
        next_node = curr_node->Ap;
        break;
      case 'C' :
        if ( curr_node->Cp == NULL ) {
          curr_node->Cp = (ktnP)init_ktn();
        }
        next_node = curr_node->Cp;
        break;
      case 'G' :
        if ( curr_node->Gp == NULL ) {
          curr_node->Gp = (ktnP)init_ktn();
        }
        next_node = curr_node->Gp;
        break;
      case 'T' :
        if ( curr_node->Tp == NULL ) {
          curr_node->Tp = (ktnP)init_ktn();
        }
        next_node = curr_node->Tp;
        break;
      default :
        // not a good base, not a good kmer, we're done
        return NULL;
      }
      curr_node = next_node;
      kmer_pos++;
    }
    /* Now, the final base => leaf node */
    next_base = toupper(kmer[kmer_pos]);
    switch(next_base) {
    case 'A':
      if ( curr_node->Ap == NULL ) {
        curr_node->Ap = (klnP)init_kln();
      }
      last_node = curr_node->Ap;
      break;
    case 'C' :
      if ( curr_node->Cp == NULL ) {
        curr_node->Cp = (klnP)init_ktn();
      }
      last_node = curr_node->Cp;
      break;
    case 'G' :
      if ( curr_node->Gp == NULL ) {
        curr_node->Gp = (klnP)init_ktn();
      }
      last_node = curr_node->Gp;
      break;
    case 'T' :
      if ( curr_node->Tp == NULL ) {
        curr_node->Tp = (klnP)init_ktn();
      }
      last_node = curr_node->Tp;
      break;
    default : // not a good base
      return NULL;
    }
    return last_node;
  }
  else {
    return NULL;
  }
}

/* Just like add_kmer, except...
   Tests if the revcom form of the kmer string is
   topologically less than the form given. If so, 
   adds that instead. This is a cheap way to handle
   the fact that reverse complement kmers can (in
   some contexts) be considered the same kmer as
   the not reverse complemented form.
*/
klnP add_canonical_kmer( const char* kmer, KSP ks ) {
  char revcom_kmer_str[ MAX_K + 1];
  strncpy( revcom_kmer_str, kmer, ks->k );
  revcom_kmer( revcom_kmer_str, ks->k );

  if ( strcmp( kmer, revcom_kmer_str ) < 0 ) {
    return add_kmer( kmer, ks );
  }
  else {
    return add_kmer( revcom_kmer_str, ks );
  }
}



/* get_kmer
   This function returns the klnP (leaf node pointer) associated with
   this kmer, if it is present in the data structure. NULL if it is
   a bad kmer or never seen
*/
klnP get_kmer( const char* kmer, KSP ks ) {
  size_t inx;
  size_t i;
  ktnP curr_node;
  ktnP next_node;
  klnP last_node;
  //ktnP curr_node, next_node;
  size_t kmer_pos;
  char next_base;

  /* Find what the index position is for the first K_AR_SIZE
     bases in this kmer */
  if ( kmer2inx( kmer, ks->k_ar_size, &inx ) ) {
  
    /* Start at the index position of the first ks->k_ar_size bases */
    curr_node = (ktnP)ks->ka[inx];
    
    /* If we've never seen that before, then this kmer is not
       present => return NULL, we're done */
    if ( curr_node == NULL ) {
      return NULL;
    }
    kmer_pos = ks->k_ar_size;

    while( kmer_pos < ks->k - 1 ) {
      next_base = toupper(kmer[kmer_pos]);
      switch( next_base ) {
      case 'A' :
        if ( curr_node->Ap == NULL ) {
	  return NULL;
        }
        next_node = curr_node->Ap;
        break;
      case 'C' :
        if ( curr_node->Cp == NULL ) {
	  return NULL;
        }
        next_node = curr_node->Cp;
        break;
      case 'G' :
        if ( curr_node->Gp == NULL ) {
	  return NULL;
        }
        next_node = curr_node->Gp;
        break;
      case 'T' :
        if ( curr_node->Tp == NULL ) {
	  return NULL;
        }
        next_node = curr_node->Tp;
        break;
      default :
        // not a good base, not a good kmer, we're done
        return NULL;
      }
      curr_node = next_node;
      kmer_pos++;
    }
    /* Now, the final base => leaf node */
    next_base = toupper(kmer[kmer_pos]);
    switch(next_base) {
    case 'A':
      return curr_node->Ap;
    case 'C' :
      return curr_node->Cp;
    case 'G' :
      return curr_node->Gp;
    case 'T' :
      return curr_node->Tp;
    default :
      return NULL;
    }
  }
  else { // didn't even see the beginning part of the kmer
    return NULL;
  }
}

klnP get_canonical_kmer( const char* kmer, KSP ks ) {
  char revcom_kmer_str[ MAX_K + 1];
  strncpy( revcom_kmer_str, kmer, ks->k );
  revcom_kmer( revcom_kmer_str, ks->k );

  if ( strcmp( kmer, revcom_kmer_str ) < 0 ) {
    return get_kmer( kmer, ks );
  }
  else {
    return get_kmer( revcom_kmer_str, ks );
  }
}


/* Returns 0 => was present, now it's gone
           1 => was never there!
*/
int remove_kmer( const char* kmer, KSP ks ) {
  size_t inx;
  size_t i;
  ktnP curr_node;
  ktnP next_node;
  klnP last_node;
  //ktnP curr_node, next_node;
  size_t kmer_pos;
  char next_base;

  /* Find what the index position is for the first K_AR_SIZE
     bases in this kmer */
  if ( kmer2inx( kmer, ks->k_ar_size, &inx ) ) {
  
    /* Start at the index position of the first ks->k_ar_size bases */
    curr_node = (ktnP)ks->ka[inx];
    
    /* If we've never seen that before, then we're done */
    if ( curr_node == NULL ) {
      return 1;
    }
    kmer_pos = ks->k_ar_size;

    while( kmer_pos < ks->k - 1 ) {
      next_base = toupper(kmer[kmer_pos]);
      switch( next_base ) {
      case 'A' :
        if ( curr_node->Ap == NULL ) {
	  return 1;
        }
        next_node = curr_node->Ap;
        break;
      case 'C' :
        if ( curr_node->Cp == NULL ) {
	  return 1;
        }
        next_node = curr_node->Cp;
        break;
      case 'G' :
        if ( curr_node->Gp == NULL ) {
	  return 1;
        }
        next_node = curr_node->Gp;
        break;
      case 'T' :
        if ( curr_node->Tp == NULL ) {
	  return 1;
        }
        next_node = curr_node->Tp;
        break;
      default :
        // not a good base, not a good kmer, we're done
        return 1;
      }
      curr_node = next_node;
      kmer_pos++;
    }
    /* Now, the final base => leaf node */
    next_base = toupper(kmer[kmer_pos]);
    switch(next_base) {
    case 'A':
      last_node = curr_node->Ap;
      curr_node->Ap = NULL;
      break;
    case 'C' :
      last_node = curr_node->Cp;
      curr_node->Cp = NULL;
      break;
    case 'G' :
      last_node = curr_node->Gp;
      curr_node->Gp = NULL;
      break;
    case 'T' :
      last_node = curr_node->Tp;
      curr_node->Tp = NULL;
      break;
    default :
      last_node = NULL;
    }
    if ( last_node == NULL ) {
      return 1;
    }
    else {
      free_kln_data( last_node );
      free( last_node );
    }
  }
  else { // didn't even see the beginning part of the kmer
    return 1;
  }
}

void free_kln_data( klnP kp ) {
  if ( kp == NULL ) {
    return;
  }
  if ( kp->data != NULL ) {
    free( kp->data );
  }
}
    
/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
             might not be null-terminated
         (2) length of the kmer
         (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error.
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index.
*/
int kmer2inx( const char* kmer,
              const size_t kmer_len,
              size_t* inx ) {
  size_t l_inx  = 0;
  int i = 0;
  char curr_char;

  while( i < kmer_len ) {
    l_inx = l_inx << 2;
    curr_char = toupper(kmer[i]); // Upper case it in case it is not
    switch( curr_char ) {
    case 'A' :
      l_inx += 0;
      break;
    case 'C' :
      l_inx += 1;
      break;
    case 'G' :
      l_inx += 2;
      break;
    case 'T' :
      l_inx += 3;
      break;
    default :
      return 0; // not valid!
    }
    i++;
  }
  *inx = l_inx;
  return 1; // valid!
}

void revcom_kmer( char* kmer, size_t k ) {
  char tmp;
  size_t i;
  for( i = 0; i < k/2; i++ ) {
    /* Put the beginning bases in the ending positions, 
       complemented and vice-versa */
    tmp = kmer[i];
    kmer[i] = revcom_base( kmer[k-(i+1)] );
    kmer[k-(i+1)] = revcom_base(tmp);
  }
  /* Even length string? We're done */
  if ( k%2 == 0 ) {
    return;
  }
  else { /* Odd length. Gotta revcom the middle base */
    kmer[i] = revcom_base( kmer[i] );
  }
}

char revcom_base( char c ) {
  switch (c) {
  case 'A' :
    return 'T';
    
  case 'C' :
    return 'G';
    
  case 'G' :
    return 'C';

  case 'T' :
    return 'A';

  case 'a' :
    return 't';
    
  case 'c' :
    return 'g';
    
  case 'g' :
    return 'c';

  case 't' :
    return 'a';

  default :
    return c;
  }
}



ktnP init_ktn( void ) {
  ktnP new_ktn;
  new_ktn = (ktnP)malloc(sizeof(ktn));
  new_ktn->Ap = NULL;
  new_ktn->Cp = NULL;
  new_ktn->Gp = NULL;
  new_ktn->Tp = NULL;
  return new_ktn;
}

klnP init_kln( void ) {
  klnP new_kln;
  new_kln = (klnP)malloc(sizeof(kln));
  new_kln->data = NULL;
  return new_kln;
}
