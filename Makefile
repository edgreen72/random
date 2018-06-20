CC=gcc
CFLAGS=-gdwarf-2 -g

file_io.o : file_io.h file_io.c
	echo "Making file_io.o ..."
	$(CC) $(CFLAGS) file_io.c -c -lz -o file_io.o 

kmer.o : kmer.h kmer.c
	echo "Making kmer.o ..."
	$(CC) $(CFLAGS) -c -o kmer.o kmer.c

test_kmer : kmer.o test_kmer.c
	echo "Making test_kmer ..."
	$(CC) $(CFLAGS) kmer.o -o test_kmer test_kmer.c

fasta-kmer-spectrum : fasta-kmer-spectrum.c kmer.o file_io.o
	echo "Making fasta-kmer-spectrum..."
	$(CC) $(CFLAGS) kmer.o file_io.o -lz -o fasta-kmer-spectrum fasta-kmer-spectrum.c
