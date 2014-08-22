#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <pthread.h>
//#include <time.h>

#include "Reads.h"
#include "KmerCode.h"
#include "Store.h"
#include "ErrorCorrection.h"

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;
bool agressiveCorrection ;
int MAX_FIX_PER_100 ;

void PrintHelp()
{
	printf( "Usage: ./a.out [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
		"\t-c jf_dump: the kmer counts dumped by JellyFish\n"
		"\t-k kmer_length\n"
		"Other parameters:\n"
		"\t-od output_file_directory (default: ./)\n"
		"\t-t number of threads to use (default: 1)\n"
		//"\t-trim allow trimming (default: false)\n"
		//"\t-all: output all the reads including those unfixable (default: false)\n"
		"\t-maxcor: the maximum number of correction every 100bp (default: 8)\n" ) ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int kmerLength = -1 ;
	Reads reads ;
	Store kmers ;
	FILE *fpJellyFishDump = NULL ;
	char buffer[100] ;

	// variables for threads
	int numOfThreads ;
	pthread_attr_t pthreadAttr ;
	pthread_t *threads ;
	pthread_mutex_t mutexErrorCorrection ;

	if ( argc == 1 )
	{
		PrintHelp() ;
		exit( 0 ) ;
	}

	numOfThreads = 1 ;
	agressiveCorrection = false ;
	MAX_FIX_PER_100 = 8 ;

	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			//reads.AddReadFile( argv[i + 1] ) ;
			++i ;
			continue ;
		}
		else if ( !strcmp( "-od", argv[i] ) )
		{
			reads.SetOutputDirectory( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-c", argv[i] ) )
		{
			fpJellyFishDump = fopen( argv[i + 1], "r" ) ;
			if ( fpJellyFishDump == NULL )
			{
				printf( "Could not open file %s\n", argv[i + 1]) ;
				exit( 1 ) ;
			}
			++i ;
		}
		else if ( !strcmp( "-k", argv[i] ) )
		{
			kmerLength = atoi( argv[i + 1] ) ;
			i += 1 ;
		}
		else if ( !strcmp( "-t", argv[i] ) )
		{
			numOfThreads = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-maxcor", argv[i] ) )
		{
			MAX_FIX_PER_100 = atoi( argv[i + 1] ) ;
			++i ;
		}
		/*else if ( !strcmp( "--agressive", argv[i] ) )
		{
			agressiveCorrection = true ;	
		}*/
		else if ( !strcmp( "-h", argv[i] ) )
		{
			PrintHelp() ;
			exit( 0 ) ;
		}
		else
		{
			printf( "Unknown argument: %s\n", argv[i] ) ;
			exit( 0 ) ;
		}
	}
	
	// Go the second round to get the reads files
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1 ] ) ;
			++i;
		}
	}
	
	KmerCode kcode( kmerLength ) ; 
	
	if ( numOfThreads > 1 )
	{
		// Initialized pthread variables
		pthread_attr_init( &pthreadAttr ) ;
		pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_JOINABLE ) ;
		threads = ( pthread_t * )malloc( sizeof( pthread_t ) * numOfThreads ) ;
		pthread_mutex_init( &mutexErrorCorrection, NULL ) ;
	}

	// Read in the kmers from the dump of JellyFish
	while ( fscanf( fpJellyFishDump, "%s", buffer ) != EOF )
	{
		int count = atoi( &buffer[1] ) ;
		fscanf( fpJellyFishDump, "%s", buffer ) ;
		if ( count <= 1 )
			continue ;

		kcode.Restart() ;
		for ( i = 0 ; buffer[i] ; ++i )
			kcode.Append( buffer[i] ) ;
		kmers.Put( kcode, count ) ;
		kcode.ReverseComplement() ;
		kmers.Put( kcode, count ) ;
	}
	// Scan the reads to get the information of the kmers.
	reads.Rewind() ;
	if ( numOfThreads == 1 )
	{
		while ( reads.Next() )
		{
			char *readId = reads.id ;
			char *seq = reads.seq ;
			char *qual ;
			if ( reads.HasQuality() )
				qual = reads.qual ;
			else
				qual = NULL ;

			int ecResult ;
			ecResult = ErrorCorrection( seq, qual, kcode, kmers ) ;	
			//if ( !strcmp( seq, "GGACTTTGAAAAGAGAGTCAAAGAGTGCTTGAAATTGTCGGGAGGGAAGGGGATGGGGGCCGGGGATGGGGCGGG" ) )
			//	exit( 1 ) ;
			/*printf( "%d\n", ecResult ) ;
			  if ( ecResult <= 0 )
			  printf( "%s\n%s\n", readId, seq ) ;*/
			//printf( "%d\n", ecResult ) ;
			reads.Output( ecResult, 0, 0, false ) ;
		}
	}
	else
	{
		int maxBatchSize = 512 * numOfThreads ;
		int batchSize ;
		struct _ErrorCorrectionThreadArg arg ;
		void *pthreadStatus ;

		struct _Read *readBatch = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;

		arg.kmerLength = kmerLength ;
		arg.kmers = &kmers ;
		arg.readBatch = readBatch ;
		arg.lock = &mutexErrorCorrection ;
		
		while ( 1 )
		{
			batchSize = reads.GetBatch( readBatch, maxBatchSize, true, true ) ;
			if ( batchSize == 0 )
				break ; 
			//printf( "batchSize=%d\n", batchSize ) ;
			arg.batchSize = batchSize ;
			arg.batchUsed = 0 ;
			for ( i = 0 ; i < numOfThreads ; ++i )
				pthread_create( &threads[i], &pthreadAttr, ErrorCorrection_Thread, (void *)&arg ) ;	

			for ( i = 0 ; i < numOfThreads ; ++i )
				pthread_join( threads[i], &pthreadStatus ) ;

			reads.OutputBatch( readBatch, batchSize, false ) ;
		}

		free( readBatch ) ;	


	}
	
	fclose( fpJellyFishDump ) ;
	return 0 ;
}
