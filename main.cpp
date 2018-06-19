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
char badQualityThreshold ; // quality <= this is bad
bool zlibVersionChecked = false ;
bool outputStdout = false ;
bool VERBOSE = false ;

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
	fprintf( stderr, "Usage: ./rcorrector [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
		"\t-p seq_file_left seq_file_right: the paths to the paired-end data set. Can use multiple -p to specifiy multiple sequence files\n"
		"\t-i seq_file: seq_file is the path to the interleaved mate-pair sequence file. Can use multiple -i\n"
		"\t-c jf_dump: the kmer counts dumped by JellyFish\n"
		"\t-k kmer_length\n"
		"Other parameters:\n"
		"\t-od output_file_directory (default: ./)\n"
		"\t-t number of threads to use (default: 1)\n"
		//"\t-trim allow trimming (default: false)\n"
		//"\t-all: output all the reads including those unfixable (default: false)\n"
		"\t-maxcor INT: the maximum number of correction every 100bp (default: 8)\n" 
		"\t-maxcorK INT: the maximum number of correction within k-bp window (default: 4)\n"
		"\t-wk FLOAT: the proportion of kmers that are used to estimate weak kmer count threshold (default: 0.95)\n"
		"\t-stdout: output the corrected sequences to stdout (default: not used)\n"
		"\t-verbose: output some correction information to stdout (default: not used)\n"
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
	fprintf( stderr, "Processed %" PRIu64 " reads\n"
		"\tCorrected %" PRIu64 " bases.\n",
		summary.totalReads, summary.totalCorrections ) ;
}

char GetBadQuality( Reads &reads )
{
	int i ;
	int qualHisto[300], firstQualHisto[300] ;
	int totalCnt, cnt ;
	int t1, t2 ;
	//Reads reads( readFile ) ;
	if ( !reads.HasQuality() )
		return 0 ;

	memset( qualHisto, 0, sizeof( qualHisto ) ) ;
	memset( firstQualHisto, 0, sizeof( firstQualHisto )) ;
	for ( i = 0 ; i < 1000000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		++qualHisto[ (int)reads.qual[ strlen( reads.seq ) - 1 ] ] ;
		++firstQualHisto[ (int)reads.qual[0] ] ;
	}

	totalCnt = i ;
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += firstQualHisto[i] ;
		if ( cnt > totalCnt * 0.05 )
			break ;
	}
	t1 = i - 1 ;
	
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += qualHisto[i] ;
		if ( cnt > totalCnt * 0.05 )
			break ;
	}
	t2 = i ;
	//printf( "%d %d\n", t1, t2 ) ;
	return (char)( t2 < t1 ? t2 : t1 ) ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int kmerLength = 23 ;
	double errorRateKmerPortion ;
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
	errorRateKmerPortion = 0.95 ;

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
		else if ( !strcmp( "-i", argv[i] ) )
		{
			++i ;
			continue ;
		}
		else if ( !strcmp( "-od", argv[i] ) )
		{
			mkdir( argv[i + 1], 0700 ) ;
			reads.SetOutputDirectory( argv[i + 1] ) ;
			pairedReads.SetOutputDirectory( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-c", argv[i] ) )
		{
			fpJellyFishDump = fopen( argv[i + 1], "r" ) ;
			if ( fpJellyFishDump == NULL )
			{
				fprintf( stderr, "Could not open file %s\n", argv[i + 1]) ;
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
		else if ( !strcmp( "-maxcorK", argv[i] ) )
		{
			MAX_FIX_PER_K = atoi( argv[i + 1] ) ;
			++i ;
		}
		/*else if ( !strcmp( "--agressive", argv[i] ) )
		{
			agressiveCorrection = true ;	
		}*/
		else if ( !strcmp( "-wk", argv[i] ) )
		{
			errorRateKmerPortion = atof( argv[i + 1 ] )  ;
			++i ;
		}
		else if ( !strcmp( "-stdout", argv[i] ) )
		{
			outputStdout = true ;
		}
		else if ( !strcmp( "-verbose", argv[i] ) )
		{
			VERBOSE = true ;
		}
		else if ( !strcmp( "-h", argv[i] ) )
		{
			PrintHelp() ;
			exit( 0 ) ;
		}
		else
		{
			fprintf( stderr, "Unknown argument: %s\n", argv[i] ) ;
			exit( 0 ) ;
		}
	}
	
	// Go the second round to get the reads files
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1 ], false, false ) ;
			++i;
		}
		else if ( !strcmp( "-p", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1], true, false ) ;
			pairedReads.AddReadFile( argv[i + 2], true, false ) ;
			i += 2 ;
		}
		else if ( !strcmp( "-i", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1], false, true ) ;
			++i ;
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
	int kmerCount = 0 ;
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
		++kmerCount ;
	}
	fprintf( stderr, "Stored %d kmers\n", kmerCount ) ;

	srand( 17 ) ;

	// Read in the first abundant 100K kmers from the dump of jellyfish to estimate 
	// the error rate.
	rewind( fpJellyFishDump ) ;
	k = 0 ;
	const int rateSize = 100000 ;
	double *rates = (double *)malloc( sizeof( double ) * ( rateSize + 1 ) ) ;
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
	rates[k] = rates[k - 1] ;
	ERROR_RATE = rates[(int)( k * errorRateKmerPortion )] ;
	if ( ERROR_RATE == 0 || k < 100 )
		ERROR_RATE = 0.01 ;
	free( rates ) ;
	fprintf( stderr, "Weak kmer threshold rate: %lf (estimated from %.3lf/1 of the chosen kmers)\n", ERROR_RATE, errorRateKmerPortion ) ;
	//exit ( 1 ) ;

	// Get the bad quality screo
	reads.Rewind() ;
	badQualityThreshold = GetBadQuality( reads ) ;
	fprintf( stderr, "Bad quality threshold is %d\n", badQualityThreshold ) ;

	// Scan the reads to get the information of the kmers.
	reads.Rewind() ;
	if ( numOfThreads == 1 )
	{
		struct _Read readBuffer[2] ;
		while ( reads.NextWithBuffer( readBuffer[0].id, readBuffer[0].seq, readBuffer[0].qual ) )
		{
			char *readId = readBuffer[0].id ;
			char *seq = readBuffer[0].seq ;
			char *qual ;
			if ( reads.HasQuality() )
				qual = readBuffer[0].qual ;
			else
				qual = NULL ;

			char *readId2, *seq2, *qual2 ;
			int t = -1, t1 = -1, t2 = -1 ;
			if ( reads.IsPaired() || reads.IsInterleaved() )
			{
				if ( reads.IsPaired() )
				{
					pairedReads.NextWithBuffer( readBuffer[1].id, readBuffer[1].seq, readBuffer[1].qual ) ;
				}
				else if ( reads.IsInterleaved() )
				{
					reads.NextWithBuffer( readBuffer[1].id, readBuffer[1].seq, readBuffer[1].qual ) ;	
				}
				readId2 = readBuffer[1].id ;
				seq2 = readBuffer[1].seq ;

				if ( qual != NULL )
					qual2 = readBuffer[1].qual ;
				else 
					qual2 = NULL ;

				t1 = GetStrongTrustedThreshold( seq, qual, kcode, kmers ) ;
				t2 = GetStrongTrustedThreshold( seq2, qual2, kcode, kmers ) ;
				t = ( t1 < t2 ) ? t1 : t2 ;
			}

			int ecResult ;
			ecResult = ErrorCorrection( readId, seq, qual, t, kcode, kmers ) ;
				
			UpdateSummary( ecResult, summary ) ;
			//if ( !strcmp( seq, "GGACTTTGAAAAGAGAGTCAAAGAGTGCTTGAAATTGTCGGGAGGGAAGGGGATGGGGGCCGGGGATGGGGCGGG" ) )
			//	exit( 1 ) ;
			/*printf( "%d\n", ecResult ) ;
			  if ( ecResult <= 0 )
			  printf( "%s\n%s\n", readId, seq ) ;*/
			//printf( "%d\n", ecResult ) ;
			readBuffer[0].correction = ecResult ;
			GetKmerInformation( seq, kmerLength, kmers, readBuffer[0].l, readBuffer[0].m, readBuffer[0].h ) ;
			reads.OutputBatch( &readBuffer[0], 1, false ) ;

			if ( reads.IsPaired() )
			{
				ecResult = ErrorCorrection( readId2, seq2, qual2, t, kcode, kmers ) ;
				readBuffer[1].correction = ecResult ; 
				GetKmerInformation( seq2, kmerLength, kmers, readBuffer[1].l, readBuffer[1].m, readBuffer[1].h ) ;
				pairedReads.OutputBatch( &readBuffer[1], 1, false ) ;
				UpdateSummary( ecResult, summary ) ;
			}
			if ( reads.IsInterleaved() )
			{
				ecResult = ErrorCorrection( readId2, seq2, qual2, t, kcode, kmers ) ;
				readBuffer[1].correction = ecResult ; 
				GetKmerInformation( seq2, kmerLength, kmers, readBuffer[1].l, readBuffer[1].m, readBuffer[1].h ) ;
				reads.OutputBatch( &readBuffer[1], 1, false) ;
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
		int fileInd1, fileInd2 ;

		arg.kmerLength = kmerLength ;
		arg.kmers = &kmers ;
		arg.readBatch = readBatch ;
		//arg.readBatch2 = NULL ;
		arg.readBatch2 = readBatch2 ;
		arg.lock = &mutexErrorCorrection ;
		
		while ( 1 )
		{
			batchSize = reads.GetBatch( readBatch, maxBatchSize, fileInd1, true, true ) ;
			if ( reads.IsPaired() )
			{
				int tmp = pairedReads.GetBatch( readBatch2, maxBatchSize, fileInd2, true, true ) ;
				if ( tmp != batchSize )
				{
					fprintf( stderr, "ERROR: The files are not paired!\n" ) ; 
					exit ( 1 ) ;
				}
				arg.readBatch2 = readBatch2 ;
			}
			else
				arg.readBatch2 = NULL ;

			if ( batchSize == 0 )
				break ; 
			//printf( "batchSize=%d\n", batchSize ) ;
			arg.batchSize = batchSize ;
			arg.batchUsed = 0 ;
			arg.interleaved = reads.IsInterleaved() ;
			for ( i = 0 ; i < numOfThreads ; ++i )
				pthread_create( &threads[i], &pthreadAttr, ErrorCorrection_Thread, (void *)&arg ) ;	

			for ( i = 0 ; i < numOfThreads ; ++i )
				pthread_join( threads[i], &pthreadStatus ) ;
			
			
			if ( outputStdout )
			{
				if ( reads.IsPaired() )
				{
					for ( i = 0 ; i < batchSize ; ++i )
					{
						reads.OutputBatch( readBatch + i, 1, false, fileInd1 ) ;
						UpdateSummary( readBatch[i].correction, summary ) ;
						pairedReads.OutputBatch( readBatch + i, 1, false, fileInd2 ) ;
						UpdateSummary( readBatch2[i].correction, summary ) ;
					}
				}
				else
				{
					reads.OutputBatch( readBatch, batchSize, false, fileInd1 ) ;
					for ( i = 0 ; i < batchSize ; ++i )
						UpdateSummary( readBatch[i].correction, summary ) ;
				}
			}	
			else
			{
				reads.OutputBatch( readBatch, batchSize, false, fileInd1 ) ;
				for ( i = 0 ; i < batchSize ; ++i )
					UpdateSummary( readBatch[i].correction, summary ) ;

				if ( reads.IsPaired() )
				{
					pairedReads.OutputBatch( readBatch2, batchSize, false, fileInd2 ) ;
					for ( i = 0 ; i < batchSize ; ++i )
						UpdateSummary( readBatch2[i].correction, summary ) ;
				}
			}
		}

		free( readBatch ) ;
		free( readBatch2 ) ;
	}

	fclose( fpJellyFishDump ) ;

	PrintSummary( summary ) ;
	return 0 ;
}
