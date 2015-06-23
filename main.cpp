#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
//#include <pthread.h>
//#include <time.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

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
int MAX_FIX_PER_K ;
double ERROR_RATE ;

struct _summary
{
	uint64_t totalCorrections ;
	uint64_t totalReads ;
} ;


int CompDouble( const void *p1, const void *p2 ) 
{
	double d = ( *(double *)p1 ) - ( *(double *)p2 ) ;
	if ( d > 0 )
		return 1 ;
	else if ( d < 0 )
		return -1 ;
	else
		return 0 ;
}

void PrintHelp()
{
	printf( "Usage: ./a.out [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
		"\t-p seq_file_left seq_file_right: the paths to the paired-end data set. Can use multiple -p to specifiy multiple sequence files\n"
		"\t-c jf_dump: the kmer counts dumped by JellyFish\n"
		"\t-k kmer_length\n"
		"Other parameters:\n"
		"\t-od output_file_directory (default: ./)\n"
		"\t-t number of threads to use (default: 1)\n"
		//"\t-trim allow trimming (default: false)\n"
		//"\t-all: output all the reads including those unfixable (default: false)\n"
		"\t-maxcor INT: the maximum number of correction every 100bp (default: 8)\n" 
		"\t-maxcorK INT: the maximum number of correction within k-bp window (default: 4)\n"
		) ;
}

void UpdateSummary( int corCnt, struct _summary &summary )
{
	++summary.totalReads ;
	if ( corCnt < 0 )
		return ;
	summary.totalCorrections += corCnt ;
}

void PrintSummary( const struct _summary &summary )
{
	printf( "Processed %" PRIu64 " reads\n"
		"\tCorrected %" PRIu64 " bases.\n",
		summary.totalReads, summary.totalCorrections ) ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int kmerLength = 23 ;
	Reads reads ;
	Reads pairedReads ;

	struct _summary summary ;

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
	MAX_FIX_PER_K = 4 ;

	summary.totalCorrections = 0 ;
	summary.totalReads = 0 ;

	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			//reads.AddReadFile( argv[i + 1] ) ;
			++i ;
			continue ;
		}
		else if ( !strcmp( "-p", argv[i] ) )
		{
			i += 2 ;
			continue ;
		}
		else if ( !strcmp( "-od", argv[i] ) )
		{
			mkdir( argv[i + 1], 0700 ) ;
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
		else if ( !strcmp( "-maxcor", argv[i] ) )
		{
			MAX_FIX_PER_K = atoi( argv[i + 1] ) ;
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
		else if ( !strcmp( "-p", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1], true ) ;
			pairedReads.AddReadFile( argv[i + 2], true ) ;
			i += 2 ;
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

	// Test how many kmers are in the jellyfish file
	/*uint64_t kmerCnt = 0 ;
	while ( fscanf( fpJellyFishDump, "%s", buffer ) != EOF )
	{
		int count = atoi( &buffer[1] ) ;
		fscanf( fpJellyFishDump, "%s", buffer ) ;
		if ( count <= 1 )
			continue ;
		++kmerCnt ;
	}
	kmers.Allocate( kmerCnt ) ;

	// Read in the kmers from the dump of JellyFish
	rewind( fpJellyFishDump ) ;*/
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
	}

	srand( 17 ) ;

	// Read in the first abundant 100K kmers from the dump of jellyfish to estimate 
	// the error rate.
	rewind( fpJellyFishDump ) ;
	k = 0 ;
	const int rateSize = 100000 ;
	double *rates = (double *)malloc( sizeof( double ) * rateSize ) ;
	rates[0] = 0 ;
	while ( fscanf( fpJellyFishDump, "%s", buffer ) != EOF && k < rateSize )
	{
		int count = atoi( &buffer[1] ) ;
		fscanf( fpJellyFishDump, "%s", buffer ) ;
		if ( !kcode.IsValid() )
			continue ;
	
		kcode.Restart() ;
		for ( i = 0 ; buffer[i] ; ++i )
			kcode.Append( buffer[i] ) ;
		count = 0 ;
		int max = 0 ;
		int secondMax = 0 ;
		for ( i = 0 ; i < 4 ; ++i )
		{
			kcode.ShiftRight( 1 ) ;
			kcode.Append( numToNuc[i] ) ;
			count = kmers.GetCount( kcode ) ;
			//printf( "%d ", count ) ;
			if ( count > max )
			{
				secondMax = max ;
				max = count ;
			}
			else if ( count > secondMax )
				secondMax = count ;	
		}
		if ( max < 1000 )
			continue ;
		rates[k] = (double)secondMax / (double)max ;
		//printf( "%d %d %lf\n", max, secondMax, rates[k] ) ;
		++k ;
	}
		
	qsort( rates, k, sizeof( rates[0] ), CompDouble ) ;
	ERROR_RATE = rates[(int)( k * 0.95 )] ;
	if ( ERROR_RATE == 0 || k < 100 )
		ERROR_RATE = 0.01 ;
	free( rates ) ;
	//printf( "%lf\n", ERROR_RATE ) ;
	//exit ( 1 ) ;

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

			char *readId2, *seq2, *qual2 ;
			int t = -1, t1 = -1, t2 = -1 ;
			if ( reads.IsPaired() )
			{
				pairedReads.Next() ;
				readId2 = pairedReads.id ;
				seq2 = pairedReads.seq ;
				
				if ( pairedReads.HasQuality() )
					qual2 = pairedReads.qual ;
				else
					qual2 = NULL ;

				t1 = GetStrongTrustedThreshold( seq, qual, kcode, kmers ) ;
				t2 = GetStrongTrustedThreshold( seq2, qual2, kcode, kmers ) ;
				t = ( t1 < t2 ) ? t1 : t2 ;
			}

			int ecResult ;
			ecResult = ErrorCorrection( seq, qual, t, kcode, kmers ) ;	
			
			UpdateSummary( ecResult, summary ) ;
			//if ( !strcmp( seq, "GGACTTTGAAAAGAGAGTCAAAGAGTGCTTGAAATTGTCGGGAGGGAAGGGGATGGGGGCCGGGGATGGGGCGGG" ) )
			//	exit( 1 ) ;
			/*printf( "%d\n", ecResult ) ;
			  if ( ecResult <= 0 )
			  printf( "%s\n%s\n", readId, seq ) ;*/
			//printf( "%d\n", ecResult ) ;
			reads.Output( ecResult, 0, 0, false ) ;

			if ( reads.IsPaired() )
			{
				ecResult = ErrorCorrection( seq2, qual2, t, kcode, kmers ) ;
				pairedReads.Output( ecResult, 0, 0, false ) ;
				
				UpdateSummary( ecResult, summary ) ;
			}

		}
	}
	else
	{
		int maxBatchSize = 512 * numOfThreads ;
		int batchSize ;
		struct _ErrorCorrectionThreadArg arg ;
		void *pthreadStatus ;

		struct _Read *readBatch = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		struct _Read *readBatch2 = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;

		arg.kmerLength = kmerLength ;
		arg.kmers = &kmers ;
		arg.readBatch = readBatch ;
		arg.readBatch2 = NULL ;
		if ( reads.IsPaired() )
		{
			arg.readBatch2 = readBatch2 ;
		}
		arg.lock = &mutexErrorCorrection ;
		
		while ( 1 )
		{
			batchSize = reads.GetBatch( readBatch, maxBatchSize, true, true ) ;
			if ( reads.IsPaired() )
			{
				int tmp = pairedReads.GetBatch( readBatch2, maxBatchSize, true, true ) ;
				if ( tmp != batchSize )
				{
					printf( "ERROR: The files are not paired!\n" ) ; 
					exit ( 1 ) ;
				}
			}
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
			for ( i = 0 ; i < batchSize ; ++i )
				UpdateSummary( readBatch[i].correction, summary ) ;

			if ( reads.IsPaired() )
			{
				pairedReads.OutputBatch( readBatch2, batchSize, false ) ;
				for ( i = 0 ; i < batchSize ; ++i )
					UpdateSummary( readBatch2[i].correction, summary ) ;
			}
			
		}

		free( readBatch ) ;
		free( readBatch2 ) ;
	}
	
	fclose( fpJellyFishDump ) ;

	PrintSummary( summary ) ;
	return 0 ;
}
