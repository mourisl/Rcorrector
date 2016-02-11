#include "ErrorCorrection.h"

#include <vector>

//#define DEBUG

#define MAX_TRIAL 1025

extern char nucToNum[26] ;
extern char numToNuc[26] ;

extern int MAX_FIX_PER_100 ;
extern int MAX_FIX_PER_K ;
extern double ERROR_RATE ;
extern bool VERBOSE ;

extern char badQualityThreshold ;

// Collect the information for fixing starting from a specific position
struct _fix
{
	int fix[MAX_READ_LENGTH] ;
	int bottleNeck ;
	int fixCnt ;
} ;

struct _segmentInfo
{
	int fixCnt ;
	int bestFixCnt ;
	int top2FixBottleNeck[2] ;
	int lanchor, ranchor ;
	int from, to ;
} ;


int CompInt( const void *p1, const void *p2 )
{
	return *(int *)p1 - *(int *)p2 ;
}

void PrintKmer( KmerCode &kcode )
{
	int kl = kcode.GetKmerLength() ;
	int i ;
	for ( i = kl - 1 ; i >= 0 ; --i )
	{
		printf( "%c", numToNuc[ ( kcode.GetCode() ) >> (uint64_t)( i * 2 ) & 3ull ] ) ;
	}
	printf( "\n" ) ;
}

bool IsPolyA( char *buffer, int length, int threshold = 2 )
{
	int i ;
	int cnt = 0 ;
	for ( i = 0 ; i < length ; ++i )
		if ( buffer[i] == 'A' )
			++cnt ;
	if ( cnt >= length - threshold )
		return true ;

	cnt = 0 ;
	for ( i = 0 ; i < length ; ++i )
		if ( buffer[i] == 'T' )
			++cnt ;
	if ( cnt >= length - threshold )
		return true ;

	return false ;
}

void *ErrorCorrection_Thread( void *arg )
{
	int i, ind ;
	int correction ;//, badPrefix, badSuffix ;
	struct _ErrorCorrectionThreadArg *myArg = ( struct _ErrorCorrectionThreadArg *)arg ; 	
	
	KmerCode kcode( myArg->kmerLength ) ;
	int increment = 1 ;
	if ( myArg->interleaved )
	{
		increment = 2 ;
	}
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		ind = myArg->batchUsed ;
		myArg->batchUsed += increment ;
		pthread_mutex_unlock( myArg->lock ) ;
		//printf( "%d %d\n", ind, myArg->batchSize ) ;
		if ( ind >= myArg->batchSize )
			break ;
		//correction = ErrorCorrection_Wrapper( myArg->readBatch[ind].seq, kmerCode, myArg->kmers, 
		//		badPrefix, badSuffix ) ;
		int t = -1, t1 = -1, t2 = -1 ;
		if ( myArg->readBatch2 != NULL || myArg->interleaved )
		{
			t1 = GetStrongTrustedThreshold( myArg->readBatch[ind].seq, myArg->readBatch[ind].qual, kcode, *myArg->kmers ) ;
			
			if ( myArg->readBatch2 != NULL )
				t2 = GetStrongTrustedThreshold( myArg->readBatch2[ind].seq, myArg->readBatch2[ind].qual, kcode, *myArg->kmers ) ;
			else if ( myArg->interleaved )
				t2 = GetStrongTrustedThreshold( myArg->readBatch[ind + 1].seq, myArg->readBatch[ind + 1].qual, kcode, *myArg->kmers ) ;
			t = ( t1 < t2 ) ? t1 : t2 ;
		}

		correction = ErrorCorrection( myArg->readBatch[ind].id, myArg->readBatch[ind].seq, myArg->readBatch[ind].qual, t, kcode, *myArg->kmers ) ;
		myArg->readBatch[ind].correction = correction ;
		myArg->readBatch[ind].badPrefix = 0 ;
		myArg->readBatch[ind].badSuffix = 0 ;
		GetKmerInformation( myArg->readBatch[ind].seq, myArg->kmerLength, *myArg->kmers, 
			myArg->readBatch[ind].l, myArg->readBatch[ind].m, myArg->readBatch[ind].h ) ;

		if ( myArg->readBatch2 != NULL )
		{
			correction = ErrorCorrection( myArg->readBatch2[ind].id, myArg->readBatch2[ind].seq, myArg->readBatch2[ind].qual, t, kcode, *myArg->kmers ) ;
			myArg->readBatch2[ind].correction = correction ;
			myArg->readBatch2[ind].badPrefix = 0 ;
			myArg->readBatch2[ind].badSuffix = 0 ;
			GetKmerInformation( myArg->readBatch2[ind].seq, myArg->kmerLength, *myArg->kmers, 
					myArg->readBatch2[ind].l, myArg->readBatch2[ind].m, myArg->readBatch2[ind].h ) ;
		}
		if ( myArg->interleaved )
		{
			correction = ErrorCorrection( myArg->readBatch[ind + 1].id, myArg->readBatch[ind + 1].seq, myArg->readBatch[ind + 1].qual, t, kcode, *myArg->kmers ) ;
			myArg->readBatch[ind + 1].correction = correction ;
			myArg->readBatch[ind + 1].badPrefix = 0 ;
			myArg->readBatch[ind + 1].badSuffix = 0 ;
			GetKmerInformation( myArg->readBatch[ind + 1].seq, myArg->kmerLength, *myArg->kmers, 
					myArg->readBatch[ind + 1].l, myArg->readBatch[ind + 1].m, myArg->readBatch[ind + 1].h ) ;
		}
	}

	pthread_exit( NULL ) ;
}


inline double GetBound( int c )
{
	return c * ERROR_RATE  + 6.0 * sqrt( c  * ERROR_RATE ) + 1 ;
}

int InferPosThreshold( KmerCode &kcode, Store &kmers, int direction, int upperBound ) // -1:left, 1:right
{
	int i ;
	int ret ;
	KmerCode tmpKcode = kcode ;
	int maxCnt = 0 ; 
	int tmp ;

	for ( i = 0 ; i < 4 ; ++i )
	{
		tmpKcode = kcode ;
		if ( direction == -1 )
			tmpKcode.Prepend( numToNuc[i] ) ;
		else 
			tmpKcode.Append( numToNuc[i] ) ;
		tmp = kmers.GetCount( tmpKcode ) ; 
		if ( tmp > maxCnt )
			maxCnt = tmp ;
	}
	
	ret = GetBound( maxCnt ) ;
	if ( ret < 1 )
		ret = 1 ;
	//return ret ;
	//tmp = GetBound( upperBound ) ;
	if ( upperBound > ret || upperBound <= 0 )
		return ret ;
	else 
		return upperBound ;
}

bool TestExtend( int t, KmerCode kcode, int depth, int direction, Store &kmers )
{
	int i ;
	KmerCode tmpKcode( kcode ) ;
	if ( depth <= 0 )
		return true ;

	for ( i = 0 ; i < 4 ; ++i )
	{
		tmpKcode = kcode ;
		if ( direction == 1 )
			tmpKcode.Append( numToNuc[i] ) ;
		else
			tmpKcode.Prepend( numToNuc[i] ) ;
#ifdef DEBUG
		printf( "%s %d %d: %d\n", __func__, i, depth, kmers.GetCount( tmpKcode ) ) ;
#endif
		if ( kmers.GetCount( tmpKcode ) < t )
			continue ;	
		if ( TestExtend( t, tmpKcode, depth - 1 , direction, kmers ) )
			return true ;
	}
	return false ;
}

// Return how many possible candidates
void SearchPaths_Right( int start, int to, int pos, int t, bool isPaired, char *seq, int fixCnt, int *fix, int fixBottleNeck,
	int &maxFixCnt, int *bestFix, int &bestFixCnt, int &bestFixBottleNeck, int top2FixBottleNeck[], bool isStrongTrusted[], bool isPolyAKmer[],
	KmerCode &kcode, Store &kmers, int &trialCnt )
{
	int i, cnt, tmpBottleNeck ;
	int extension = 0 ;
	int kmerLength = kcode.GetKmerLength() ;
	int threshold ;

	// Give up fixing this read
	if ( trialCnt > MAX_TRIAL )
	{
		if ( maxFixCnt > 2 )
		{
			--maxFixCnt ;
			bestFixCnt = 0 ; // Remove the result before
			bestFixBottleNeck = -1 ;
			trialCnt = 0 ;
		}
		else
			return ;
	}
	if ( fixCnt > maxFixCnt )
		return ;
	//if ( !isPaired && fixCnt == maxFixCnt && fixBottleNeck < bestFixBottleNeck )
	//	return ;

	if ( pos >= to )
	{
		if ( seq[pos] == '\0' )
		{
			/*int teThreshold = ( t < fixBottleNeck ? t : fixBottleNeck ) ;
			int tmp = teThreshold ;
			teThreshold = teThreshold - 2 * sqrt( teThreshold ) ;
			if ( teThreshold < tmp / 2 )
				teThreshold = tmp / 2  ;
			if ( teThreshold < 1 )
				teThreshold = 1 ;*/
			//if ( !TestExtend( 1, kcode, 1, 1, kmers ) ) 
			//	return ;
		}
		
		if ( fixBottleNeck < t )
			++fixCnt ;

		if ( fixCnt < maxFixCnt )
		{
			top2FixBottleNeck[0] = fixBottleNeck ;
			top2FixBottleNeck[1] = -1 ;
		}
		else if ( fixCnt == maxFixCnt )
		{
			if ( fixBottleNeck > top2FixBottleNeck[0] )
			{
				top2FixBottleNeck[1] = top2FixBottleNeck[0] ;
				top2FixBottleNeck[0] = fixBottleNeck ;
			}
			else if ( fixBottleNeck > top2FixBottleNeck[1] )
			{
				top2FixBottleNeck[1] = fixBottleNeck ;
			}
		}
		
		//printf( "%d %d %d %d\n", fixCnt, fixBottleNeck, maxFixCnt, bestFixBottleNeck ) ;
		if ( fixCnt < maxFixCnt 
		     || /*( !isPaired && fixCnt == maxFixCnt 
		     		&& (double)rand() / RAND_MAX <= (double)fixBottleNeck / ( fixBottleNeck + bestFixBottleNeck ) ) 
		     || ( isPaired &&*/( fixCnt == maxFixCnt && fixBottleNeck > bestFixBottleNeck ) )
		{
			if ( fixCnt < maxFixCnt )
				trialCnt = -( maxFixCnt - fixCnt + 1 ) * MAX_TRIAL ;
			for ( i = start ; i < pos ; ++i )
				bestFix[i] = fix[i] ;
			maxFixCnt = fixCnt ;
			//if ( fixCnt == maxFixCnt )
			//	bestFixBottleNeck = 0 ;
			bestFixBottleNeck = fixBottleNeck ;
			bestFixCnt = 1 ;
		}
		else if ( fixCnt == maxFixCnt && fixBottleNeck == bestFixBottleNeck )
		{
			bestFixCnt += 1 ;
		}
		return ;
	}

	threshold = InferPosThreshold( kcode, kmers, 1, t ) ;

	KmerCode tmpKcode = kcode ;
	/*for ( i = 0 ; seq[i] ; ++i )
		printf( "%d ", fix[i] ) ;
	printf( "\n" ) ;*/

	if ( nucToNum[ seq[pos] - 'A' ] != -1 )
	{
		tmpKcode.Append( seq[pos] ) ;	
		cnt = kmers.GetCount( tmpKcode ) ; 
#ifdef DEBUG
		printf( "right: %d %d %c=>%c (%d, %d)\n", fixCnt, pos, seq[pos], seq[pos], cnt, bestFixBottleNeck ) ;
		//PrintKmer( tmpKcode ) ;
#endif
		if ( cnt >= threshold ) 
		{
			fix[pos] = -1 ; 
			tmpBottleNeck = fixBottleNeck ;
			if ( cnt < tmpBottleNeck )
				tmpBottleNeck = cnt ;
			++extension ;
			SearchPaths_Right( start, to, pos + 1, t, isPaired, seq, fixCnt, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
					trialCnt ) ;		
		}
		else if ( threshold == 1 && t <= 2 ) 
		{
			// See if it is accidentally below the threshold, by jumping
			int j, k = 0 ;
			for ( i = pos  ; cnt < threshold && k < kmerLength ; ++k )	
			{
				++i ;
				if ( i >= to )
					break ;
				tmpKcode.Append( seq[i] ) ;
				cnt = kmers.GetCount( tmpKcode ) ;
				//PrintKmer( tmpKcode ) ;
			}
		
			if ( k < kmerLength && i < to )
			{
				for ( j = pos ; j <= i ; ++j )
					fix[j] = -1 ;
				++extension ;
				tmpBottleNeck = fixBottleNeck ;
				//tmpKcode.Prepend( 'A' ) ;
				SearchPaths_Right( start, to, i + 1, t, isPaired, seq, fixCnt + 1, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
					trialCnt ) ;
			}
		}
	}

	char c[5] = "ACGT" ;

	if ( isStrongTrusted[pos] == false && !isPolyAKmer[ pos - kmerLength + 1 ] )
	{
		for ( i = 0 ; i < 4 ; ++i )
		{
			if ( c[i] == seq[pos] )
				continue ;
			tmpKcode = kcode ;
			tmpKcode.Append( c[i] ) ;
			cnt = kmers.GetCount( tmpKcode ) ; 
			if ( cnt >= threshold ) 
			{
				int tmpFixCnt = fixCnt ;
				fix[pos] = i ;

				++trialCnt ;
				++extension ;
#ifdef DEBUG
				printf( "right: %d %d %c=>%c (%d, %d)\n", fixCnt, pos, seq[pos], c[i], cnt, bestFixBottleNeck ) ;
#endif
				if ( nucToNum[ seq[pos] - 'A' ] != -1 )
					++tmpFixCnt ;
				tmpBottleNeck = fixBottleNeck ;
				if ( cnt < tmpBottleNeck )
					tmpBottleNeck = cnt ;
				SearchPaths_Right( start, to, pos + 1, t, isPaired, seq, tmpFixCnt, fix, tmpBottleNeck,
						maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, 
						isPolyAKmer, tmpKcode, kmers, trialCnt ) ;
			}
		}
	}
	/*else if ( isPolyAKmer[ pos - kmerLength + 1 ] )
	{
		// Jump through the PolyA region.
		tmpKcode = kcode ;
		int dest = pos + kmerLength - 1 ;
		for ( i = pos ; seq[i] && i < dest ; ++i )
		{
			tmpKcode.Append( seq[i] ) ;
			if ( nucToNum[ seq[i] - 'A' ] < 0 )
			{
				dest += kmerLength ;
			}
		}
		SearchPaths_Right( start, dest + 1, t, seq, fixCnt, fix, fixBottleNeck,
				maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
				trialCnt ) ;

	}*/

	// Jump through some unfixable region
	if ( extension == 0 )
	{
		int j, k ;
		tmpKcode = kcode ;
		for ( i = pos ; i < to ; ++i )
		{
			threshold = InferPosThreshold( tmpKcode, kmers, 1, t ) ;

			tmpKcode.Append( seq[i] ) ;
			cnt = kmers.GetCount( tmpKcode ) ;
			fix[i] = -1 ;
			if ( cnt >= threshold )
				break ;
		}
		tmpBottleNeck = fixBottleNeck ;
		/*tmpBottleNeck = 0 ;
		for ( j = 0 ; j < 4 ; ++j )
		{
			KmerCode enumKmerCode = kcode ;
			enumKmerCode.Append( numToNuc[j] ) ;
			tmpBottleNeck += kmers.GetCount( enumKmerCode ) ;
		}
		if ( tmpBottleNeck > fixBottleNeck )
			tmpBottleNeck = fixBottleNeck ;
		else if ( tmpBottleNeck <= 0 )
			tmpBottleNeck = 1 ;*/

		if ( seq[i] )
		{
			k = i - pos - kmerLength + 1 ;
		}
		else
		{
			k = ( i - pos ) / 2 ;
			//tmpBottleNeck = 0 ;
		}
		if ( k <= 0 )
			k = 1 ;
#ifdef DEBUG
		printf( "jump right: %d=>%d(%d)\n", pos, i, k ) ;
#endif
		fixCnt += k ;
		if ( i >= to )
			i -= 1 ;
		//tmpKcode.Prepend( 'A' ) ;
		SearchPaths_Right( start, to, i + 1, t, isPaired, seq, fixCnt, fix, tmpBottleNeck,
				maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
				trialCnt ) ;
	}
}

void SearchPaths_Left( int start, int to, int pos, int t, bool isPaired, char *seq, int fixCnt, int *fix, int fixBottleNeck, 
	int &maxFixCnt, int *bestFix, int &bestFixCnt, int &bestFixBottleNeck, int top2FixBottleNeck[], bool isStrongTrusted[], bool isPolyAKmer[], 
	KmerCode &kcode, Store &kmers, int &trialCnt )
{
	int i, cnt, tmpBottleNeck ;
	int extension = 0 ;
	int kmerLength = kcode.GetKmerLength() ;
	int threshold = t ;

	if ( trialCnt > MAX_TRIAL )
	{
		if ( maxFixCnt > 2 )
		{
			--maxFixCnt ;
			bestFixCnt = 0 ; // Remove the result before
			bestFixBottleNeck = -1 ;
			trialCnt = 0 ;
		}
		else
			return ;
	}	
	
	if ( fixCnt > maxFixCnt ) // TODO: it is possible that a wrong path has smaller fix point but with much lower bottleneck
		return ;
	if ( pos < to )
	{
		if ( pos <= -1 )
		{
			/*int teThreshold = ( t < fixBottleNeck ? t : fixBottleNeck ) ;
			int tmp = teThreshold ;
			teThreshold = teThreshold - 2 * sqrt( teThreshold ) ;
			if ( teThreshold < tmp / 2 )
				teThreshold = tmp / 2  ;
			if ( teThreshold < 1 )
				teThreshold = 1 ;*/
			//if ( !TestExtend( 1, kcode, 3, -1, kmers ) ) 
			//	return ;
		}
		
		if ( fixBottleNeck < t )
			++fixCnt ;
		
		if ( fixCnt < maxFixCnt )
		{
			top2FixBottleNeck[0] = fixBottleNeck ;
			top2FixBottleNeck[1] = -1 ;
		}
		else if ( fixCnt == maxFixCnt )
		{
			if ( fixBottleNeck > top2FixBottleNeck[0] )
			{
				top2FixBottleNeck[1] = top2FixBottleNeck[0] ;
				top2FixBottleNeck[0] = fixBottleNeck ;
			}
			else if ( fixBottleNeck > top2FixBottleNeck[1] )
			{
				top2FixBottleNeck[1] = fixBottleNeck ;
			}
		}
		//printf( "%d %d %d %d\n", fixCnt, fixBottleNeck, maxFixCnt, bestFixBottleNeck ) ;
		if ( fixCnt < maxFixCnt 
		     || /*( !isPaired && fixCnt == maxFixCnt 
		     		&& (double)rand() / RAND_MAX <= (double)fixBottleNeck / ( fixBottleNeck + bestFixBottleNeck ) ) 
		     || ( isPaired &&*/ ( fixCnt == maxFixCnt && fixBottleNeck > bestFixBottleNeck ) )
		{
			if ( fixCnt < maxFixCnt )
				trialCnt = -( maxFixCnt - fixCnt + 1 ) * MAX_TRIAL ;
			for ( i = start ; i > pos ; --i )
				bestFix[i] = fix[i] ;
			maxFixCnt = fixCnt ;
			//if ( fixCnt == maxFixCnt )
			//	bestFixBottleNeck = 0 ;
			bestFixBottleNeck = fixBottleNeck ;
			bestFixCnt = 1 ;
		}
		else if ( fixCnt == maxFixCnt && fixBottleNeck == bestFixBottleNeck )
		{
			bestFixCnt += 1 ;
		}
		return ;
	}
	threshold = InferPosThreshold( kcode, kmers, -1, t ) ;

	KmerCode tmpKcode = kcode ;
	/*for ( i = 0 ; seq[i] ; ++i )
		printf( "%d ", fix[i] ) ;
	printf( "\n" ) ;*/

	if ( nucToNum[ seq[pos] - 'A' ] != -1 )
	{
		tmpKcode.Prepend( seq[pos] ) ;	
		cnt = kmers.GetCount( tmpKcode ) ; 
#ifdef DEBUG
		printf( "left: %d %d %c=>%c (%d[%d], %d)\n", fixCnt, pos, seq[pos], seq[pos], cnt, threshold, bestFixBottleNeck ) ;
#endif
		if ( cnt >= threshold ) 
		{
			fix[pos] = -1 ; 
			tmpBottleNeck = fixBottleNeck ;
			if ( cnt < tmpBottleNeck )
				tmpBottleNeck = cnt ;
			++extension ;
			SearchPaths_Left( start, to, pos - 1, threshold, isPaired, seq, fixCnt, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers, 
					trialCnt ) ;		
		}
		else if ( threshold == 1 && t <= 2 ) 
		{
			// See if it is accidentally below the threshold
			int j, k = 0 ;
			for ( i = pos ; cnt < threshold && k < kmerLength ; ++k )	
			{
				--i ;
				if ( i < to )
					break ;
				tmpKcode.Prepend( seq[i] ) ;
				cnt = kmers.GetCount( tmpKcode ) ;
			}
		
			if ( k < kmerLength && i >= to )
			{
				for ( j = i ; j <= pos ; ++j )
					fix[j] = -1 ;
				tmpBottleNeck = fixBottleNeck ;
				++extension ;
				SearchPaths_Left( start, to, i - 1, threshold, isPaired, seq, fixCnt + 1, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, 
					isPolyAKmer, tmpKcode, kmers, trialCnt ) ;
			}
		}
	}

	char c[5] = "ACGT" ;
	if ( isStrongTrusted[pos] == false && !isPolyAKmer[pos] ) 
	{
		for ( i = 0 ; i < 4 ; ++i )
		{
			if ( c[i] == seq[pos] )
				continue ;
			tmpKcode = kcode ;

			tmpKcode.Prepend( c[i] ) ;
			cnt = kmers.GetCount( tmpKcode ) ; 
			//printf( "%d: %c %d %lld %s\n", pos, c[i], cnt, tmpKcode.GetCode(), &seq[pos] ) ;
			if ( cnt >= threshold ) 
			{
				int tmpFixCnt = fixCnt ;
				fix[pos] = i ;
				++trialCnt ;

#ifdef DEBUG
				printf( "left: %d %d %c=>%c (%d[%d], %d)\n", fixCnt, pos, seq[pos], c[i], cnt, threshold, bestFixBottleNeck ) ;
#endif
				if ( nucToNum[ seq[pos] - 'A' ] != -1 )
					++tmpFixCnt ;
				tmpBottleNeck = fixBottleNeck ;
				if ( cnt < tmpBottleNeck )
					tmpBottleNeck = cnt ;
				++extension ;
				SearchPaths_Left( start, to, pos - 1, threshold, isPaired, seq, tmpFixCnt, fix, tmpBottleNeck,
						maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, 
						isStrongTrusted, isPolyAKmer, tmpKcode, kmers, trialCnt ) ;
			}
		}
	}
	/*else if ( isPolyAKmer[ pos] )
	{
		// Jump through the PolyA region.
		tmpKcode = kcode ;
		int dest = pos - kmerLength + 1 ;
		for ( i = pos ; i >= 0 && i > dest ; --i )
		{
			tmpKcode.Prepend( seq[i] ) ;
			if ( nucToNum[ seq[i] - 'A'] < 0 )
			{
				dest -= kmerLength ;
			}
		}
		SearchPaths_Left( start, dest - 1, t, seq, fixCnt, fix, fixBottleNeck,
				maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
				trialCnt ) ;

	}*/
	
	// Jump through some unfixable region
	if ( extension == 0 )
	{
		int k ;
		tmpKcode = kcode ;
		for ( i = pos ; i >= to ; --i )
		{
			threshold = InferPosThreshold( tmpKcode, kmers, -1, t ) ;
			tmpKcode.Prepend( seq[i] ) ;
			cnt = kmers.GetCount( tmpKcode ) ;
			fix[i] = -1 ;
			if ( cnt >= threshold )
				break ;
		}
		
		tmpBottleNeck = fixBottleNeck ;
		/*tmpBottleNeck = 0 ;
		int j ;
		for ( j = 0 ; j < 4 ; ++j )
		{
			KmerCode enumKmerCode = kcode ;
			enumKmerCode.Prepend( numToNuc[j] ) ;
			tmpBottleNeck += kmers.GetCount( enumKmerCode ) ;
		}
		if ( tmpBottleNeck > fixBottleNeck )
			tmpBottleNeck = fixBottleNeck ;
		else if ( tmpBottleNeck <= 0 )
			tmpBottleNeck = 1 ;*/
		
		if ( i >= 0 )
		{
			k = pos - i - kmerLength + 1 ;
		}
		else
		{
			k = pos - 1 ;
			//tmpBottleNeck = 0 ;
		}
		if ( k <= 0 )
			k = 1 ;
#ifdef DEBUG
		printf( "jump left: %d=>%d(%d)\n", pos, i, k ) ;
#endif
		fixCnt += k ;	
		if ( i <= to )
			++i ;
		SearchPaths_Left( start, to, i - 1, threshold, isPaired, seq, fixCnt, fix, tmpBottleNeck,
				maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, top2FixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
				trialCnt ) ;
	}
}



int ErrorCorrection( char *id, char *seq, char *qual, int pairStrongTrustThreshold, KmerCode &kcode, Store &kmers )
{
	//printf( "Correct %s\n", seq ) ; fflush( stdout ) ;

	if ( VERBOSE )
	{
		printf( "%s\n", id ) ;
	}
	int i, j, k ;
	int kcnt = 0 ;
	int kmerLength = kcode.GetKmerLength() ;
	int readLength ; 
	int counts[MAX_READ_LENGTH] ;
	int iBuffer[MAX_READ_LENGTH], fix[MAX_READ_LENGTH] ;
	double dBuffer[MAX_READ_LENGTH] ;
	bool isStrongTrusted[MAX_READ_LENGTH] ; 
	bool isPolyAKmer[MAX_READ_LENGTH] ;
	int tstart, tend, longestTrustedLength, trustThreshold ;
	//int searchResult ;
	int totalFix= 0, allowedFix ;
	int badSegmentCnt = 0 ;
	int ret ;
	int fixBottleNeck, bestFixBottleNeck ;
	int strongTrustThreshold ;
	bool unfixable, forceNextRound ;
	bool isPaired ;
	bool flag ;

	struct _segmentInfo segmentInfo[MAX_READ_LENGTH], trustedIsland[MAX_READ_LENGTH] ;
	int segmentInfoCnt = 0, trustedIslandCnt = 0 ;

	kcode.Restart() ;
	for ( i = 0 ; i < kmerLength - 1 ; ++i )
		kcode.Append( seq[i] ) ;
	for ( ; seq[i] ; ++i, ++kcnt )
	{
		kcode.Append( seq[i] ) ;
		counts[kcnt] = kmers.GetCount( kcode ) ;
	}

	//printf( "\n",seq ) ;
	readLength = i ;
	allowedFix = readLength ;//* MAX_FIX_PER_100 / 100 ;	
	unfixable = false ;
	
	isPaired = pairStrongTrustThreshold == -1 ? false : true ;
	//isPaired = true ;

	//printf( "%d\n", allowedFix ) ;
	// Test whether there are too many Ns
	k = 0 ; 
	for ( i = 0 ; i < readLength ; ++i )
		if ( seq[i] == 'N' )
			++k ;
	if ( k > 5 )
		return -1 ;

	// Test whether there are too many A
	k = 0 ;
	for ( i = 0 ; i < readLength ; ++i )
		if ( seq[i] == 'A' )
			++k ;
	if ( k > readLength - kmerLength )
		return -1 ;
	
	k = 0 ;
	for ( i = 0 ; i < readLength ; ++i )
		if ( seq[i] == 'T' )
			++k ;
	if ( k > readLength - kmerLength )
		return -1 ;


	//printf( "%s\n", seq ) ;
	if ( VERBOSE )
	{
		printf( "Before correction:\n%s\n", seq ) ;
		for ( i = 0 ; i < kcnt ; ++i )
		{
			if ( counts[i] != 0 )
				printf( "%d ", counts[i] ) ;
			else
				printf( "1 " ) ;
		}
		printf( "\n" ) ;
	}

	trustThreshold = 2 ;

	for ( i = 0 ; i < kcnt ; ++i )
	{
		int t = 7 ;
		if ( kmerLength / 2 > t )
			t = kmerLength / 2 ;
		if ( IsPolyA( seq + i, kmerLength, t ) )
			iBuffer[i] = -1 ;
		else
			iBuffer[i] = counts[i] ;
	}
	qsort( iBuffer, kcnt, sizeof( int ), CompInt ) ;


	for ( i = kcnt - 1 ; i >= 1  ; --i )
	{
		if ( iBuffer[i] > 2 * iBuffer[i - 1]  && iBuffer[i] > 10 ) 
			break ;
	}

	flag = false ;
	if ( i >= 1 )
	//if ( i < kcnt )
	{
		//trustThreshold = iBuffer[( kcnt - 1 + i ) / 2] / 20 + 1 ;
		trustThreshold = GetBound( iBuffer[i] ) ;
		strongTrustThreshold = iBuffer[i] ;
		if ( strongTrustThreshold >= 20 && iBuffer[i - 1] == 2 && trustThreshold < 3 )
		{
			flag = true ;
			trustThreshold = 3 ;
		}
		//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[i] ) ;
	}
	else
	{
		//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[kcnt / 2 ] ) ;

		for ( i = 0 ; i < kcnt ; ++i )
			if ( iBuffer[i] > 0 )
				break ;
		strongTrustThreshold = iBuffer[ ( i + kcnt - 1 ) / 2 ] ;
		trustThreshold = GetBound( strongTrustThreshold ) ;
	}

	if ( pairStrongTrustThreshold >= 1 && strongTrustThreshold > pairStrongTrustThreshold )
	{
		if ( flag == false || pairStrongTrustThreshold < 20 )
			trustThreshold = GetBound( pairStrongTrustThreshold ) ;
		strongTrustThreshold = pairStrongTrustThreshold ;
	}
	for ( i = 0 ; i < kcnt ; ++i )
	{
		if ( IsPolyA( seq + i, kmerLength ) )
			isPolyAKmer[i] = true ;
		else
			isPolyAKmer[i] = false ;
	}

/*#ifdef DEBUG
	printf( "isPolyAKmer:\n") ;
	for ( i = 0 ; i < kcnt ; ++i )
	{
		printf( "%d", isPolyAKmer[i] ) ;
	}
	printf( "\n" ) ;
#endif*/
	//trustThreshold /= 2 ;
	if ( trustThreshold < 2 )
		trustThreshold = 2 ;
	//if ( trustThreshold > 10 )
	//	trustThreshold = 10 ;

	//if ( trustThreshold >= 10 )
	//	allowedFix -= trustThreshold / 10 ;
	//if ( allowedFix < 0 )
	//	allowedFix = 0 ;


	// Try different thresholds
	int iter = 0 ;
	while ( 1 )
	{
		if ( VERBOSE )
			printf( "strong trust threshold=%d threshold=%d\n", strongTrustThreshold, trustThreshold ) ;
		
		allowedFix = readLength ;//* MAX_FIX_PER_100 / 100 ;
		totalFix = 0 ;
		unfixable = false ;
		forceNextRound = false ;
		// Find the longest trusted portion
		trustedIslandCnt = 0 ;
		longestTrustedLength = -1 ;
		j = 0 ;
		memset( isStrongTrusted, false, sizeof( bool ) * readLength ) ;
		for ( i = 0 ; i < readLength ; ++i )
			iBuffer[i] = -1 ;
		for ( i = 0 ; i < kcnt ; ++i )
		{
			//int t = 7 ;
			//if ( kmerLength / 2 > t )
			//	t = kmerLength / 2 ;
			if ( counts[i] >= strongTrustThreshold && !IsPolyA( seq + i, kmerLength ) )
			{
				//printf( "(%d %d)\n", counts[i], strongTrustThreshold ) ;
				++j ;
			}
			else
			{
				if ( j > longestTrustedLength )
				{
					longestTrustedLength = j ;
					tstart = i - longestTrustedLength ;
					tend = i - 1 ;
				}

				if ( j >= 2 )
				{
					//printf( "%d %d: %d %d\n", i, j, i - j, i - 1 + kmerLength - 1 ) ;
					/*for ( k = i - j ; k <= i - 1 + kmerLength - 1 ; ++k )
					{
						if ( iBuffer[k] == -1 )
							isStrongTrusted[k] = true ;
						else
							isStrongTrusted[k] = false ;
						//iBuffer[k] = i - j ; 
					}*/
					trustedIsland[ trustedIslandCnt ].from = i - j ;
					trustedIsland[ trustedIslandCnt ].to = i - 1 ;
					++trustedIslandCnt ;
				}
				/*if ( IsPolyA( seq + i, kmerLength ) )
				{
					for ( k = i ; k < i + kmerLength ; ++i )
						isStrongTrusted[k] = true ;
				}*/
				j = 0 ;
			}
		}
		if ( j > longestTrustedLength )
		{
			longestTrustedLength = j ;
			tstart = i - longestTrustedLength ;
			tend = i - 1 ;
		}
		if ( j >= 2 )
		{
			/*for ( k = i - j ; k <= i - 1 + kmerLength - 1 ; ++k )
			{
				if ( iBuffer[k] == -1 )
					isStrongTrusted[k] = true ;
				else
					isStrongTrusted[k] = false ;
				//iBuffer[k] = i - j ; 
			}*/
			trustedIsland[ trustedIslandCnt ].from = i - j ;
			trustedIsland[ trustedIslandCnt ].to = i - 1 ;
			++trustedIslandCnt ;
		}

		// Adjust the boundary of trusted kmer islands
		for ( i = 1 ; i < trustedIslandCnt ; ++i )
		{
			if ( trustedIsland[i].from <= trustedIsland[i - 1].to + kmerLength - 1 + 1 )	
			{
				int len1 = trustedIsland[i - 1].to - trustedIsland[i - 1].from ;
				int len2 = trustedIsland[i].to - trustedIsland[i].from ;

				int overlap = trustedIsland[i - 1].to + kmerLength - 1 + 1 - trustedIsland[i].from ;
				//printf( "%d\n", overlap ) ;
				// check the count between them
				for ( j = trustedIsland[i - 1].to + 1 ; j < trustedIsland[i].from ; ++j )
				{
					if ( counts[j] <= 2 && counts[j] < trustThreshold )
						break ;
				}

				if ( j >= trustedIsland[i].from )
					continue ;
				if ( overlap > 3 )
					continue ;

				if ( len1 < len2 )
				{
					// change the i-1 island
					trustedIsland[i - 1].to -= ( overlap + 1 ) ;
				}
				else
				{
					trustedIsland[i].from += ( overlap + 1 ) ;
				}
			}
		}

		// Set trusted positions
		for ( i = 0 ; i < trustedIslandCnt ; ++i )
		{
			if ( trustedIsland[i].from > trustedIsland[i].to )
				continue ;
			for ( j = trustedIsland[i].from ; j <= trustedIsland[i].to + kmerLength - 1 ; ++j )
				isStrongTrusted[j] = true ;
		}
		
	
		trustedIslandCnt = 0 ;
		j = -1 ;
		for ( i = 0 ; seq[i] ; ++i )
		{
			if ( j == -1 && isStrongTrusted[i] == true )
			{
				j = i ;
			}
			if ( j != -1 && isStrongTrusted[i] == false && isStrongTrusted[i - 1] == true )
			{
				trustedIsland[ trustedIslandCnt ].from = j ;
				trustedIsland[ trustedIslandCnt ].to = i - 1 ;
				++trustedIslandCnt ;
				j = -1 ;
			}
		}

		if ( j != -1 )
		{
			trustedIsland[ trustedIslandCnt ].from = j ;
			trustedIsland[ trustedIslandCnt ].to = i - 1 ;
			++trustedIslandCnt ;
			j = -1 ;
		}

		if ( trustedIslandCnt == 0 )
		{
			trustedIsland[0].from = tstart ;
			trustedIsland[0].to = tend + kmerLength - 1 ;
			trustedIslandCnt = 1 ;
		}

		segmentInfoCnt = 0 ;
		if ( trustedIsland[0].from > 0 )
		{
			segmentInfo[ segmentInfoCnt ].from = 0;
			segmentInfo[ segmentInfoCnt ].to = trustedIsland[0].from - 1 ; 
			segmentInfo[ segmentInfoCnt ].fixCnt = -1 ;
			segmentInfo[ segmentInfoCnt ].bestFixCnt = -1 ;
			segmentInfo[ segmentInfoCnt ].top2FixBottleNeck[0] = -1 ;
			segmentInfo[ segmentInfoCnt ].top2FixBottleNeck[1] = -1 ;
			segmentInfo[ segmentInfoCnt ].lanchor = 0 ;
			segmentInfo[ segmentInfoCnt ].ranchor = trustedIsland[0].to - trustedIsland[0].from + 1 ;
			++segmentInfoCnt ;
		}

		for ( i = 0 ; i < trustedIslandCnt - 1 ; ++i )
		{
			segmentInfo[ segmentInfoCnt ].from = trustedIsland[i].to + 1 ;
			segmentInfo[ segmentInfoCnt ].to = trustedIsland[i + 1].from - 1 ; 
			segmentInfo[ segmentInfoCnt ].fixCnt = -1 ;
			segmentInfo[ segmentInfoCnt ].bestFixCnt = -1 ;
			segmentInfo[ segmentInfoCnt ].top2FixBottleNeck[0] = -1 ;
			segmentInfo[ segmentInfoCnt ].top2FixBottleNeck[1] = -1 ;
			segmentInfo[ segmentInfoCnt ].lanchor = trustedIsland[i].to - trustedIsland[i].from + 1 ;
			segmentInfo[ segmentInfoCnt ].ranchor = trustedIsland[i + 1].to - trustedIsland[i + 1].from + 1 ;
			++segmentInfoCnt ;
		}
		if ( trustedIsland[i].to < readLength - 1 )
		{
			segmentInfo[ segmentInfoCnt ].from = trustedIsland[i].to + 1 ;
			segmentInfo[ segmentInfoCnt ].to = readLength ; // the range is actually [from, to)
			segmentInfo[ segmentInfoCnt ].fixCnt = -1 ;
			segmentInfo[ segmentInfoCnt ].bestFixCnt = -1 ;
			segmentInfo[ segmentInfoCnt ].top2FixBottleNeck[0] = -1 ;
			segmentInfo[ segmentInfoCnt ].top2FixBottleNeck[1] = -1 ;
			segmentInfo[ segmentInfoCnt ].lanchor = trustedIsland[i].to - trustedIsland[i].from + 1 ;
			segmentInfo[ segmentInfoCnt ].ranchor = 0 ;
			++segmentInfoCnt ;
		}

		// Remove some of the strong trusted positions if they are next to untrusted regions shorter than k
		// We don't need to adjust boundary cases
		/*for ( i = 0 ; i < readLength - 1 && isStrongTrusted[i] && !isStrongTrusted[i + 1] ; ++i )
			;
		i = i - kmerLength + 2 ;
		for (  ; i < kcnt ; ++i )
		{
			if ( counts[i] < trustThreshold )
			{
				++j ;
			}
			else
			{
				if ( j < kmerLength && j > 0 )
				{
					for ( k = i - 1 - kmerLength + 1  ; k < i - j ; ++k )
						isStrongTrusted[k] = false ;
					for ( k = i ; k < i + ( kmerLength - j ) && k < readLength ; ++k )
						isStrongTrusted[k] = false ;
				}
				j = 0 ;
			}
		}

		// Adjust the longest trusted region
		if ( !isStrongTrusted[ tstart ] )
		{	
			while ( tstart < kcnt && !isStrongTrusted[tstart] )
				++tstart ;
		}

		if ( !isStrongTrusted[ tend + kmerLength - 1 ] )
		{
			while ( tend >= 0 && !isStrongTrusted[ tend + kmerLength - 1 ] )
				++tend ;
		}

		longestTrustedLength = tend - tstart + 1 ;
		if ( longestTrustedLength < 0 )
			longestTrustedLength = -1 ;*/
		if ( VERBOSE )
		{
			printf( "Is corresponding base strong trusted?\n" ) ;
			for ( i = 0 ; i < readLength ; ++i )
				printf( "%d", isStrongTrusted[i] ) ;
			printf( "\n" ) ;
		}

#ifdef DEBUG
		printf( "trustedIslandCnt=%d segmentInfoCnt=%d\n", trustedIslandCnt, segmentInfoCnt ) ;
#endif
		//printf( "trusted region:%d %d: %c %d\n", tstart, tend, seq[tstart], longestTrustedLength ) ;


		//if ( strongTrustThreshold >= 20 && longestTrustedLength >= 5 )
		//	printf( "AS events?\n" ) ;

		//return 0 ;

		if ( longestTrustedLength == -1 )
			return -1 ;

		if ( longestTrustedLength == kcnt )
			return 0 ;

		for ( i = 0 ; i < readLength ; ++i )
			fix[i] = -1 ;

		// Scan towards right
		badSegmentCnt = 0 ;	
		if ( segmentInfoCnt > 0 )
		{
			int bestFixCnt = -1 ;
			int trialCnt = 0 ;
			int maxFixCnt = allowedFix ;
			
			bestFixBottleNeck = INF ; 
			for ( k = 0 ; k < segmentInfoCnt ; ++k )
			{
				trialCnt = 0 ;
				maxFixCnt = ( segmentInfo[k].to - segmentInfo[k].from + 1 ) * MAX_FIX_PER_K / kmerLength * 2 + 1 ; //allowedFix ;
				//printf( "%d %d %d %d\n", segmentInfo[k].from, segmentInfo[k].to, MAX_FIX_PER_K, maxFixCnt ) ;
				if ( maxFixCnt < MAX_FIX_PER_K )
					maxFixCnt = MAX_FIX_PER_K ;

				fixBottleNeck = 1000000000 ;
				int tmpBestFixBottleNeck = -1 ;

				if ( segmentInfo[k].lanchor >= segmentInfo[k].ranchor )
				{
					int extend = ( segmentInfo[k].to == readLength  ? 0 : ( kmerLength - 1 ) ) ; 
					tend = segmentInfo[k].from - kmerLength ;
					kcode.Restart() ;
					for ( i = tend ; i < tend + kmerLength ; ++i )
						kcode.Append( seq[i] ) ;
					SearchPaths_Right( tend + kmerLength, segmentInfo[k].to + extend, tend + kmerLength, 
							trustThreshold, isPaired, seq, 0, iBuffer, fixBottleNeck, maxFixCnt, fix, bestFixCnt, 
							tmpBestFixBottleNeck, segmentInfo[k].top2FixBottleNeck, isStrongTrusted, 
							isPolyAKmer, kcode, kmers, trialCnt ) ;
				}
				else // Search left
				{
					int extend = ( segmentInfo[k].from == 0 ? 0 : ( kmerLength - 1 ) ) ; 
					tstart = segmentInfo[k].to + 1 ;
					kcode.Restart() ;
					for ( i = tstart ; i < tstart + kmerLength ; ++i )
						kcode.Append( seq[i] ) ;

					SearchPaths_Left( tstart - 1, segmentInfo[k].from - extend, tstart - 1, 
							trustThreshold, isPaired, seq, 0, iBuffer, fixBottleNeck, maxFixCnt, fix, bestFixCnt, 
							tmpBestFixBottleNeck, segmentInfo[k].top2FixBottleNeck, isStrongTrusted, 
							isPolyAKmer, kcode, kmers, trialCnt ) ;
				}
				//printf( "%d\n", tmpBestFixBottleNeck ) ;
				if ( tmpBestFixBottleNeck == -1 )
				{
					++badSegmentCnt ;
					continue ; // ignore a bad segment.
				}

				if ( tmpBestFixBottleNeck < bestFixBottleNeck )
					bestFixBottleNeck = tmpBestFixBottleNeck ;
				if ( bestFixBottleNeck == -1 )
					break ;
				if ( trialCnt > MAX_TRIAL )
					return -1 ;
				totalFix += maxFixCnt ;
			}
			// Compute bestFixCnt ;
			//printf( "%d\n", bestFixBottleNeck ) ;
			if ( bestFixBottleNeck != -1 )
			{
				bestFixCnt = 1 ;
				for ( i = 0 ; i < segmentInfoCnt ; ++i )
				{
					//printf( "%d %d\n", segmentInfo[i].top2FixBottleNeck[0], segmentInfo[i].top2FixBottleNeck[1] ) ;
					if ( segmentInfo[i].top2FixBottleNeck[1] >= bestFixBottleNeck )
						bestFixCnt *= 2 ;
					else if ( segmentInfo[i].top2FixBottleNeck[0] >= bestFixBottleNeck )
					{
						bestFixCnt *= 1 ;
					}
					// the two if statements above should always be true
				}
			}
			//printf( "%d %d\n", bestFixCnt, maxFixCnt ) ;

			if ( bestFixBottleNeck != -1 && iter == 0 && bestFixBottleNeck < GetBound( strongTrustThreshold ) )
			{
				/*int t = bestFixBottleNeck ;
				fixBottleNeck = 1000000000 ;
				bestFixBottleNeck = -1 ; 
				for ( i = tend + kmerLength ; i < readLength ; ++i )
					fix[i] = -1 ;
				bestFixCnt = 0 ;
				trialCnt = 0 ;
				maxFixCnt = allowedFix ;

				SearchPaths_Right( tend + kmerLength, tend + kmerLength, t, seq, 0, iBuffer, fixBottleNeck, 
						maxFixCnt, fix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, kcode, kmers, trialCnt ) ;*/
				//bestFixCnt = 0 ;
				//return -1 ;
				forceNextRound = true ;
			}

			//printf( "%d %lf %d\n", bestFixBottleNeck, GetBound( strongTrustThreshold ), forceNextRound  ) ;
			/*if ( maxFixCnt == 0 && bestFixBottleNeck < strongTrustThreshold && bestFixBottleNeck >= GetBound( strongTrustThreshold )  
				&& readLength - tend - kmerLength > kmerLength )
			{
				forceNextRound = true ;
				//strongTrustThreshold = bestFixBottleNeck ;
			}*/

			/*printf( "bestFixCnt=%d\n", bestFixCnt ) ;
			for ( i = 0 ; i < readLength ; ++i )
				printf( "%d ", fix[i] ) ;
			printf( "\n" ) ;*/
			//printf( "%d\n", trialCnt ) ;
			if ( bestFixCnt >= 2 )
				return -1 ;
			else if ( bestFixCnt <= 0 )
				unfixable = true ;
		}
		if ( totalFix == 0 && forceNextRound )
			return 0 ;
		//printf( "totalFix=%d\n", totalFix ) ;
		if ( totalFix > allowedFix )
			unfixable = true ;
		if ( unfixable == false && !forceNextRound )
			break ;

		// Maybe due to alternative splicing, a highly expressed 
		// txpt may shadow a less expressed txpt.
		if ( trustThreshold < 10 && !forceNextRound )
			return -1 ;
		//return -1 ;
		// And also the gap between

		//return -1 ;	
		for ( i = 0 ; i < kcnt ; ++i )
		{
			int t = 7 ;
			if ( kmerLength / 2 > t )
				t = kmerLength / 2 ;
			if ( IsPolyA( seq + i, kmerLength, t ) )
				iBuffer[i] = -1 ;
			else
				iBuffer[i] = counts[i] ;
		}
		qsort( iBuffer, kcnt, sizeof( int ), CompInt ) ;
		bool hasDrop = false ;
		for ( i = kcnt - 1 ; i >= 1  ; --i )
		{
			if ( iBuffer[i] > strongTrustThreshold )
				continue ;
			if ( iBuffer[i] > 2 * iBuffer[i - 1]  && iBuffer[i] > 10 ) 
			{
				hasDrop = true ;
				if ( iBuffer[i] < strongTrustThreshold )
					break ;
			}
			else if ( iBuffer[i - 1] == 0 && iBuffer[i] >= 5 )
			{
				hasDrop = true ;
				if ( iBuffer[i] < strongTrustThreshold )
					break ;
			}
		}
		if ( hasDrop )
		{
			//printf( "Iteration %d: ", iter ) ;
			++iter ;
			//printf( "%d", trustThreshold ) ;
			//trustThreshold = iBuffer[( kcnt - 1 + i ) / 2] / 20 + 1 ;
			trustThreshold = GetBound( iBuffer[i] ) ;
			strongTrustThreshold = iBuffer[i] ;
			//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[i] ) ;
			//printf( "=>%d\n", trustThreshold ) ;
			//fflush( stdout ) ;
		}
		else
			break ;
			//return -1 ;
	}
	
	//printf( "totalFix=%d\n", totalFix ) ;
	// If the fix is due to the repeats, then we remove it.
	//bool undo[MAX_READ_LENGTH] ;
	int cnt = 0 ;
	for ( i = 0 ; i < readLength ; ++i )
	{
		if ( seq[i] == 'N' || fix[i] == -1 )
			continue ;
		iBuffer[ cnt ] = i ;
		++cnt ;
	}
		
	//printf( "\n" ) ;
#ifdef DEBUG
	for ( i = 0 ; i < readLength ; ++i )
		printf( "%d", fix[i] == -1 ? 5 : fix[i] ) ;
	printf( "\n" ) ;
#endif

	// If multiple errors within a kmer window did not decrease the count much,
	// then, they are not errors.
	for ( i = 1 ; i < cnt ; ++i )
	{
		if ( qual[0] != '\0' && ( qual[ iBuffer[i] ] <= badQualityThreshold && qual[ iBuffer[i - 1] ] <= badQualityThreshold ) )
			continue ;
		if ( iBuffer[i] - iBuffer[i - 1] + 1 <= kmerLength )
		{
			//printf( "%d %d\n", counts[ iBuffer[i - 1]], counts[ iBuffer[i] ] ) ;
			j = iBuffer[i - 1] - kmerLength + 1 ;
			int minSingle = INF ;
			int minDouble = INF ;
			int taga = -1, tagb = -1 ;

			if ( j < 0 )
				j = 0 ;
			for ( ; j < kcnt ; ++j )
			{
				if ( i >= 2 && j <= iBuffer[i - 2] )
					continue ;
				if ( i < cnt - 1 && j + kmerLength - 1 >= iBuffer[i + 1] )
					break ;

				if ( j + kmerLength - 1 >= iBuffer[i] )
					break ;
				if ( counts[j] < minSingle )
				{
					minSingle = counts[j] ;
					taga = j ;
				}
			}

			for ( ; j < kcnt ; ++j )
			{
				if ( i < cnt - 1 && j + kmerLength - 1 >= iBuffer[i + 1] )
					break ;
				if ( j > iBuffer[i - 1] )
					break ;
				if ( counts[j] < minDouble )
				{
					minDouble = counts[j] ;
					tagb = j ;
				}
			}
		
			for ( ; j < kcnt ; ++j )
			{
				if ( i < cnt - 1 && j + kmerLength - 1 >= iBuffer[i + 1] )
					break ;

				if ( j > iBuffer[i] )
					break ;
				if ( counts[j] < minSingle )
				{
					minSingle = counts[j] ;
					taga = j ;
				}
			}


#ifdef DEBUG
			printf( "%d(%d) %d(%d): %d %d\n", taga, iBuffer[i - 1], tagb, iBuffer[i], counts[taga], counts[tagb] ) ;
#endif
			if ( 	minSingle != INF && minDouble != INF &&
				counts[ taga ] > 1 &&
				counts[ tagb ] > 1 &&
				counts[ taga ] > counts[ tagb ] / 2 &&
				counts[ taga ] < 2 * counts[ tagb ] ) 
			{
				fix[ iBuffer[i] ] = -1 ;
				fix[ iBuffer[i - 1] ] = -1 ;

				j = i - 2 ;
				while ( j >= 0 && iBuffer[j + 1] - iBuffer[j] + 1 <= kmerLength )
				{
					fix[ iBuffer[j] ] = -1 ;
					--j ;
				}

				while ( i + 1 < cnt && iBuffer[i + 1] - iBuffer[i] + 1 <= kmerLength )
				{
					fix[ iBuffer[i + 1] ] = -1 ;	
					++i ;
				}
			}
		} 
	}
#ifdef DEBUG
	for ( i = 0 ; i < readLength ; ++i )
		printf( "%d", fix[i] == -1 ? 5 : fix[i] ) ;
	printf( "\n" ) ;
#endif


	// If the fix are at the boundary and fixed all the position, then we should remove these fixes
	if ( totalFix > 3 && readLength > 10 )
	{
		int tmp = 0 ;
		for ( i = 0 ; i < 10 ; ++i )
			if ( fix[i] != -1 && seq[i] != 'N' && qual[i] > badQualityThreshold )
				++tmp ; 
		if ( tmp >= 2 )
			for ( i = 0 ; i < 10 ; ++i )
			{
				if ( seq[i] != 'N' )
					fix[i] = -1 ;
			}

		tmp = 0 ;
		for ( i = readLength - 10 ; i < readLength ; ++i )
			if ( fix[i] != -1 && seq[i] != 'N' && qual[i] > badQualityThreshold )
				++tmp ; 
		if ( tmp >= 3 )
			for ( i = readLength - 10 ; i < readLength ; ++i )
			{
				if ( seq[i] != 'N' )
					fix[i] = -1 ;
			}
	}	

	if ( totalFix >= MAX_FIX_PER_K )
	{
		dBuffer[0] = 0 ;
		for ( i = 0 ; i < readLength ; ++i )
			if ( seq[i] != 'N' && fix[i] != -1 )
			{
				if ( qual[i] > badQualityThreshold )
					dBuffer[i + 1] = dBuffer[i] + 1 ;
				else
					dBuffer[i + 1] = dBuffer[i] + 0.5 ;
			}
			else
				dBuffer[i + 1] = dBuffer[i] ;

		for ( i = 0 ; i < kcnt ; ++i )
		{
			int threshold = MAX_FIX_PER_K ;
			// Allowing slightly more on the 3' end
			//if ( i > 0.7 * readLength )
			//	++threshold ;

			if ( dBuffer[i + kmerLength] - dBuffer[i] > threshold )
				return -1 ;
		}

		/*if ( readLength > 10 )
		{
			if ( iBuffer[10] - iBuffer[0] > 3 || 
				iBuffer[readLength] - iBuffer[readLength - 11] > 3 )
			{
				// undo the fix

			}
		}*/
	}
	
	ret = 0 ;
	for ( i = 0 ; seq[i] ; ++i )
	{
		if ( fix[i] != -1 )
		{
			seq[i] = numToNuc[ fix[i] ] ;
			++ret ;
		}
	}
	if ( ret == 0 && badSegmentCnt > 0 )
		return -1 ;
	return ret ;
}

int GetStrongTrustedThreshold( char *seq, char *qual, KmerCode &kcode, Store &kmers ) 
{
	int i, k ;
	int kcnt = 0 ;
	int counts[MAX_READ_LENGTH] ;
	int iBuffer[MAX_READ_LENGTH] ;
	int readLength ;
	int kmerLength = kcode.GetKmerLength() ;

	kcode.Restart() ;
	for ( i = 0 ; i < kmerLength - 1 ; ++i )
		kcode.Append( seq[i] ) ;
	for ( ; seq[i] ; ++i, ++kcnt )
	{
		kcode.Append( seq[i] ) ;
		counts[kcnt] = kmers.GetCount( kcode ) ;
	}

	//printf( "\n",seq ) ;
	readLength = i ;
	//printf( "%d\n", allowedFix ) ;
	// Test whether there are too many Ns
	k = 0 ; 
	for ( i = 0 ; i < readLength ; ++i )
		if ( seq[i] == 'N' )
			++k ;
	if ( k > 5 )
		return -1 ;

	// Test whether there are too many A
	k = 0 ;
	for ( i = 0 ; i < readLength ; ++i )
		if ( seq[i] == 'A' )
			++k ;
	if ( k > readLength - kmerLength )
		return -1 ;

	k = 0 ;
	for ( i = 0 ; i < readLength ; ++i )
		if ( seq[i] == 'T' )
			++k ;
	if ( k > readLength - kmerLength )
		return -1 ;


	for ( i = 0 ; i < kcnt ; ++i )
	{
		int t = 7 ;
		if ( kmerLength / 2 > t )
			t = kmerLength / 2 ;
		if ( IsPolyA( seq + i, kmerLength, t ) )
			iBuffer[i] = -1 ;
		else
			iBuffer[i] = counts[i] ;
	}
	qsort( iBuffer, kcnt, sizeof( int ), CompInt ) ;


	for ( i = kcnt - 1 ; i >= 1  ; --i )
	{
		if ( iBuffer[i] > 2 * iBuffer[i - 1]  && iBuffer[i] > 10 ) 
			break ;
	}
	if ( i >= 1 )
		//if ( i < kcnt )
	{
		//trustThreshold = iBuffer[( kcnt - 1 + i ) / 2] / 20 + 1 ;
		//trustThreshold = GetBound( iBuffer[i] ) ;
		return iBuffer[i] ;
		//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[i] ) ;
	}
	else
	{
		//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[kcnt / 2 ] ) ;
		for ( i = 0 ; i < kcnt ; ++i )
			if ( iBuffer[i] > 0 )
				break ;
		return iBuffer[ ( i + kcnt - 1 ) / 2] ;
	}

}

void GetKmerInformation( char *seq, int kmerLength, Store &kmers, int &l, int &m, int &h )
{
	int kmerCount[MAX_READ_LENGTH] ;
	int i ;
	int k = 0 ;
	KmerCode kcode( kmerLength ) ;
	for ( i = 0 ; seq[i] && i < kmerLength - 1 ; ++i )
		kcode.Append( seq[i] ) ;
	l = m = h = 0 ;
	for ( ; seq[i] ; ++i )		
	{
		kcode.Append( seq[i] ) ;
		if ( kcode.IsValid() )
		{
			kmerCount[k] = kmers.GetCount( kcode ) ;
			if ( kmerCount[k] == 0 )
				kmerCount[k] = 1 ;
			++k ;
		}
	}
	if ( k == 0 )
		return ;

	if ( VERBOSE )	
	{
		printf( "After coorrection:\n" ) ;
		printf( "%s\n", seq ) ;
		for ( i = 0 ; i < k ; ++i )
			printf( "%d ", kmerCount[i] ) ;
		printf( "\n" ) ;
	}
	qsort( kmerCount, k, sizeof( int ), CompInt ) ;
	l = kmerCount[0] ;
	m = kmerCount[ k / 2 ] ;
	h = kmerCount[ k - 1 ] ;
}
