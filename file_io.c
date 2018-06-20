#include "file_io.h"

ChrP newSeq( void ) {
  ChrP seq;
  seq = (ChrP)malloc(sizeof(Chr));
  seq->len = 0;
  return seq;
}

/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}

/* gz_read_next_fasta
   args 1. gzFile pointer to file to be read
        2. ChrP to put the sequence
   returns: FALSE (0) if no problem
            non-zero if there is a problem
*/
int gz_read_next_fasta( gzFile genome_file, ChrP chr ) {
  char c;
  size_t i;
  c = gzgetc( genome_file );
    if ( c == EOF ) return -1;
    if ( c != '>' ) return -1;

  // get id
  i = 0;
  while( (!isspace( c=gzgetc( genome_file ) ) &&
	  (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    chr->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      //id is too long - truncate it
      chr->id[i] = '\0';
    }
  }
  chr->id[i] = '\0';

  // everything else on this line is description, if there is anything
  if ( c == '\n' ) {
    ;
  }
  else { // not end of line, so skip past the stupid whitespace...
    while( (c != '\n') &&
	   (isspace(c)) ) {
      c = gzgetc( genome_file );
    }
    ///...everthing else is description
    i = 0;
    gzungetc( c, genome_file );
    while ( c != '\n' ) {
      c = gzgetc( genome_file );
    }
  }

  // read sequence
  i = 0;
  c = gzgetc( genome_file );
  while ( ( c != '>' ) &&
          ( c != EOF ) &&
	  (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      chr->seq[i++] = c;
    }
    c = gzgetc( genome_file );
  }
  chr->seq[i] = '\0';

  chr->len = i;

  if ( c == '>' ) {
    gzungetc( '>', genome_file );
    return 0;
  }

  /* Run up against the sequence length limit so truncate it here,
     wind through the genome_file filehandle, and return this guy */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '>') &&
	    (c != EOF) ) {
      c = gzgetc( genome_file );
    }
    if ( c == '>' ) {
      gzungetc( '>', genome_file );
    }
    fprintf( stderr, "%s is longer than allowed length: %d\n",
	     chr->id, MAX_SEQ_LEN );
    return 0;
  }

  return 0;
}



/* read_next_fasta
   args 1. pointer to file to be read
        2. ChrP to put the sequence
   returns: FALSE (0) if no problem
            non-zero if there is a problem
*/
int read_next_fasta( FILE* genome_file, ChrP chr ) {
  char c;
  size_t i;
  c = fgetc( genome_file );
  if ( c == EOF ) return -1;
  if ( c != '>' ) return -1;

  // get id
  i = 0;
  while( (!isspace( c=fgetc( genome_file ) ) &&
	  (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return -1;
    }
    chr->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      //id is too long - truncate it
      chr->id[i] = '\0';
    }
  }
  chr->id[i] = '\0';

  // everything else on this line is description, if there is anything
  if ( c == '\n' ) {
    ;
  }
  else { // not end of line, so skip past the stupid whitespace...
    while( (c != '\n') &&
	   (isspace(c)) ) {
      c = fgetc( genome_file );
    }
    ///...everthing else is description
    i = 0;
    ungetc( c, genome_file );
    while ( c != '\n' ) {
      c = fgetc( genome_file );
    }
  }

  // read sequence
  i = 0;
  c = fgetc( genome_file );
  while ( ( c != '>' ) &&
          ( c != EOF ) &&
	  (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      chr->seq[i++] = c;
    }
    c = fgetc( genome_file );
  }
  chr->seq[i] = '\0';

  chr->len = i;

  if ( c == '>' ) {
    ungetc( '>', genome_file );
    return 0;
  }

  /* Run up against the sequence length limit so truncate it here,
     wind through the genome_file filehandle, and return this guy */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '>') &&
	    (c != EOF) ) {
      c = fgetc( genome_file );
    }
    if ( c == '>' ) {
      ungetc( '>', genome_file );
    }
    fprintf( stderr, "%s is longer than allowed length: %d\n",
	     chr->id, MAX_SEQ_LEN );
    return 0;
  }

  return 0;
}

/* is_gz
   Returns true if the filename argument ends in .gz
*/
int is_gz( const char* fq_fn ) {
  size_t fn_len;
  fn_len = strlen( fq_fn );
  if ( (fq_fn[fn_len-3] == '.') &&
       (fq_fn[fn_len-2] == 'g') &&
       (fq_fn[fn_len-1] == 'z') ) {
    return 1;
  }
  return 0;
}

/* Args: FILE* ffq - file pointer to forward fastq file
         FILE* rfq - file pointer to reverse fastq file
         SQP fqpair - pointer to struct sqp (SQP) for holding read pair
   Returns: 0 (false) if everything got loaded, i.e., no problem
            non-zero if EOF or a problem
   Reads the next fastq sequences in the input files whose
   file pointers are passed. Stores their info in the SQP passed in
*/
int read_next_fastqs( FILE* ffq, FILE* rfq, SQP fqpair ) {
  int frs, rrs;

  frs = read_fastq( ffq, fqpair->fid, fqpair->fseq, fqpair->fqual,
		    &fqpair->flen );
  rrs = read_fastq( rfq, fqpair->rid, fqpair->rseq, fqpair->rqual,
		    &fqpair->flen );

  if ( frs || rrs ) {
    return 1;
  }
  return 0;
}


/* Args: gzFile ffq - gzFile pointer to forward fastq file
         gzFile rfq - gzFile pointer to reverse fastq file
         SQP fqpair - pointer to struct sqp (SQP) for holding read pair
   Returns: 0 (false) if everything got loaded, i.e., no problem
            non-zero if EOF or a problem
   Reads the next fastq sequences in the input files whose
   file pointers are passed. Stores their info in the SQP passed in
*/
int gz_read_next_fastqs( gzFile gzffq, gzFile gzrfq, SQP fqpair ) {
  int frs, rrs;

  frs = gzread_fastq( gzffq, fqpair->fid, fqpair->fseq, fqpair->fqual,
		      &fqpair->flen );
  rrs = gzread_fastq( gzrfq, fqpair->rid, fqpair->rseq, fqpair->rqual,
		      &fqpair->rlen );

  if ( frs || rrs ) {
    return -1;
  }
  return 0;
}



/* Args: FILE* fastq - file pointer to a fastq file ready to read 
                       the next sequence
         char id[] - pointer to character array for the ID
         char seq[] - pointer to the character array for the seq
         char qual[] - pointer to the char array for the qscores
   Returns: 0 - everything is copacetic
            non-zero if EOF or other problem
*/
int read_fastq( FILE* fastq, char id[], char seq[], char qual[], size_t* len ) {
  char c;
  size_t i;
  c = fgetc( fastq );
  if ( c == EOF ) return -1;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return -1;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
	  (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return -1;
    }
    id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      id[i] = '\0';
    }
  }
  id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
	  (c != EOF) ) {
    c = fgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
  while ( (c != '\n') &&
	  (c != EOF) &&
	  (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      seq[i++] = c;
    }
    c = fgetc( fastq );
  }
  seq[i] = '\0';
  *len = i;
  
  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
	 (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
  i = 0;
  while( (c != '\n') &&
	 (c != EOF) &&
	 (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      qual[i++] = c;
    }
    c = fgetc( fastq );
  }
  qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return -1;
  }
  return 0;
}

/* Args: gzFile fastq - file pointer to a fastq file ready to read 
                       the next sequence
         char id[] - pointer to character array for the ID
         char seq[] - pointer to the character array for the seq
         char qual[] - pointer to the char array for the qscores
   Returns: 0 - everything is copacetic
            non-zero if EOF or other problem
*/
int gzread_fastq( gzFile fastq, char id[], char seq[], char qual[], size_t* len ) {
  char c;
  size_t i;
  c = gzgetc( fastq );
  if ( c == EOF ) return -1;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return -1;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=gzgetc( fastq ) ) &&
	  (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return -1;
    }
    id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      id[i] = '\0';
    }
  }
  id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
	  (c != EOF) ) {
    c = gzgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = gzgetc( fastq );
  while ( (c != '\n') &&
	  (c != EOF) &&
	  (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      seq[i++] = c;
    }
    c = gzgetc( fastq );
  }
  seq[i] = '\0';
  *len = i;
  
  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = gzgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = gzgetc( fastq );
  while( (c != '\n') &&
	 (c != EOF) ) {
    c = gzgetc( fastq );
  }

  /* Now, get the quality score line */
  c = gzgetc( fastq );
  i = 0;
  while( (c != '\n') &&
	 (c != EOF) &&
	 (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      qual[i++] = c;
    }
    c = gzgetc( fastq );
  }
  qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
	    (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return -1;
  }
  return 0;
}
