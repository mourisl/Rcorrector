#include "ErrorCorrection.h"

//#define DEBUG

#define MAX_TRIAL 10000

extern char nucToNum[26] ;
extern char numToNuc[26] ;

extern int MAX_FIX_PER_100 ;

// Collect the information for fixing starting from a specific position
/*struct _fix
{
	int fix[MAX_READ_LENGTH] ;
} ;*/


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
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		ind = myArg->batchUsed ;
		++myArg->batchUsed ;
		pthread_mutex_unlock( myArg->lock ) ;
		//printf( "%d %d\n", ind, myArg->batchSize ) ;
		if ( ind >= myArg->batchSize )
			break ;
		//correction = ErrorCorrection_Wrapper( myArg->readBatch[ind].seq, kmerCode, myArg->kmers, 
		//		badPrefix, badSuffix ) ;

		correction = ErrorCorrection( myArg->readBatch[ind].seq, myArg->readBatch[ind].qual, kcode, *myArg->kmers ) ;
		myArg->readBatch[ind].correction = correction ;
		myArg->readBatch[ind].badPrefix = 0 ;
		myArg->readBatch[ind].badSuffix = 0 ;
	}

	pthread_exit( NULL ) ;
}


inline double GetBound( int c )
{
	return c / 300.0 + 8.0 * sqrt( c / 300.0 ) + 1 ;
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
	return ret ;
	/*tmp = GetBound( upperBound ) ;
	if ( tmp > ret )
		return ret ;
	else
		return tmp ;*/
}

// Return how many possible candidates
void SearchPaths_Right( int start, int pos, int t, char *seq, int fixCnt, int *fix, int fixBottleNeck,
	int &maxFixCnt, int *bestFix, int &bestFixCnt, int &bestFixBottleNeck, bool isStrongTrusted[], bool isPolyAKmer[],
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
	if ( fixCnt == maxFixCnt && fixBottleNeck < bestFixBottleNeck )
		return ;
	if ( !seq[pos] )
	{
		if ( fixCnt < maxFixCnt || ( fixCnt == maxFixCnt && fixBottleNeck > bestFixBottleNeck ) )
		{
			//printf( "%d %d %d %d\n", fixCnt, fixBottleNeck, maxFixCnt, bestFixBottleNeck ) ;
			if ( fixCnt < maxFixCnt )
				trialCnt = -( maxFixCnt - fixCnt + 1 ) * MAX_TRIAL ;
			for ( i = start ; i < pos ; ++i )
				bestFix[i] = fix[i] ;
			maxFixCnt = fixCnt ;
			bestFixBottleNeck = fixBottleNeck ;
			bestFixCnt = 1 ;
		}
		else if ( fixCnt == maxFixCnt && fixBottleNeck == bestFixBottleNeck )
		{
			bestFixCnt += 1 ;
		}
		return ;
	}

	if ( t == -1 )
		threshold = InferPosThreshold( kcode, kmers, 1, fixBottleNeck ) ;
	else
		threshold = t ;

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
			SearchPaths_Right( start, pos + 1, t, seq, fixCnt, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
					trialCnt ) ;		
		}
		else if ( threshold == 1 ) 
		{
			// See if it is accidentally below the threshold, by jumping
			int j, k = 0 ;
			for ( i = pos  ; cnt < threshold && k < kmerLength ; ++k )	
			{
				++i ;
				if ( !seq[i] )
					break ;
				tmpKcode.Append( seq[i] ) ;
				cnt = kmers.GetCount( tmpKcode ) ;
				//PrintKmer( tmpKcode ) ;
			}
		
			if ( k < kmerLength && seq[i] )
			{
				for ( j = pos ; j <= i ; ++j )
					fix[j] = -1 ;
				++extension ;
				tmpBottleNeck = fixBottleNeck ;
				//tmpKcode.Prepend( 'A' ) ;
				SearchPaths_Right( start, i + 1, t, seq, fixCnt, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
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
				SearchPaths_Right( start, pos + 1, t, seq, tmpFixCnt, fix, tmpBottleNeck,
						maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
						trialCnt ) ;
			}
		}
	}
	// Jump through some unfixable region
	if ( extension == 0 )
	{
		int j, k ;
		tmpKcode = kcode ;
		for ( i = pos ; seq[i] ; ++i )
		{
			if ( t == -1 )
				threshold = InferPosThreshold( tmpKcode, kmers, 1, fixBottleNeck ) ;

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
		if ( !seq[i] )
			i -= 1 ;
		//tmpKcode.Prepend( 'A' ) ;
		SearchPaths_Right( start, i + 1, t, seq, fixCnt, fix, tmpBottleNeck,
				maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
				trialCnt ) ;
	}
}

void SearchPaths_Left( int start, int pos, int t, char *seq, int fixCnt, int *fix, int fixBottleNeck, 
	int &maxFixCnt, int *bestFix, int &bestFixCnt, int &bestFixBottleNeck, bool isStrongTrusted[], bool isPolyAKmer[], 
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
	if ( pos < 0 )
	{
		//printf( "%d %d %d %d\n", fixCnt, fixBottleNeck, maxFixCnt, bestFixBottleNeck ) ;
		if ( fixCnt < maxFixCnt || ( fixCnt == maxFixCnt && fixBottleNeck > bestFixBottleNeck ))
		{
			if ( fixCnt < maxFixCnt )
				trialCnt = -( maxFixCnt - fixCnt + 1 ) * MAX_TRIAL ;
			for ( i = start ; i > pos ; --i )
				bestFix[i] = fix[i] ;
			maxFixCnt = fixCnt ;
			bestFixBottleNeck = fixBottleNeck ;
			bestFixCnt = 1 ;
		}
		else if ( fixCnt == maxFixCnt && fixBottleNeck == bestFixBottleNeck )
		{
			bestFixCnt += 1 ;
		}
		return ;
	}
	if ( threshold == -1 )
		threshold = InferPosThreshold( kcode, kmers, -1, fixBottleNeck ) ;

	KmerCode tmpKcode = kcode ;
	/*for ( i = 0 ; seq[i] ; ++i )
		printf( "%d ", fix[i] ) ;
	printf( "\n" ) ;*/

	if ( nucToNum[ seq[pos] - 'A' ] != -1 )
	{
		tmpKcode.Prepend( seq[pos] ) ;	
		cnt = kmers.GetCount( tmpKcode ) ; 
#ifdef DEBUG
		printf( "left: %d %d %c=>%c (%d, %d)\n", fixCnt, pos, seq[pos], seq[pos], cnt, bestFixBottleNeck ) ;
#endif
		if ( cnt >= threshold ) 
		{
			fix[pos] = -1 ; 
			tmpBottleNeck = fixBottleNeck ;
			if ( cnt < tmpBottleNeck )
				tmpBottleNeck = cnt ;
			++extension ;
			SearchPaths_Left( start, pos - 1, threshold, seq, fixCnt, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers, 
					trialCnt ) ;		
		}
		else if ( threshold == 1 ) 
		{
			// See if it is accidentally below the threshold
			int j, k = 0 ;
			for ( i = pos ; cnt < threshold && k < kmerLength ; ++k )	
			{
				--i ;
				if ( i < 0 )
					break ;
				tmpKcode.Prepend( seq[i] ) ;
				cnt = kmers.GetCount( tmpKcode ) ;
			}
		
			if ( k < kmerLength && i >= 0 )
			{
				for ( j = i ; j <= pos ; ++j )
					fix[j] = -1 ;
				tmpBottleNeck = fixBottleNeck ;
				++extension ;
				SearchPaths_Left( start, i - 1, threshold, seq, fixCnt, fix, tmpBottleNeck,
					maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
					trialCnt ) ;
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
				printf( "left: %d %d %c=>%c (%d, %d)\n", fixCnt, pos, seq[pos], c[i], cnt, bestFixBottleNeck ) ;
#endif
				if ( nucToNum[ seq[pos] - 'A' ] != -1 )
					++tmpFixCnt ;
				tmpBottleNeck = fixBottleNeck ;
				if ( cnt < tmpBottleNeck )
					tmpBottleNeck = cnt ;
				++extension ;
				SearchPaths_Left( start, pos - 1, threshold, seq, tmpFixCnt, fix, tmpBottleNeck,
						maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
						trialCnt ) ;
			}
		}
	}
	
	// Jump through some unfixable region
	if ( extension == 0 )
	{
		int k ;
		tmpKcode = kcode ;
		for ( i = pos ; i >= 0 ; --i )
		{
			if ( t == -1 )
				threshold = InferPosThreshold( tmpKcode, kmers, -1, fixBottleNeck ) ;
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
		if ( i <= 0 )
			++i ;
		SearchPaths_Left( start, i - 1, threshold, seq, fixCnt, fix, tmpBottleNeck,
				maxFixCnt, bestFix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, tmpKcode, kmers,
				trialCnt ) ;
	}
}



int ErrorCorrection( char *seq, char *qual, KmerCode &kcode, Store &kmers )
{
	//printf( "Correct %s\n", seq ) ; fflush( stdout ) ;
	int i, j, k ;
	int kcnt = 0 ;
	int kmerLength = kcode.GetKmerLength() ;
	int readLength ; 
	int counts[MAX_READ_LENGTH] ;
	int iBuffer[MAX_READ_LENGTH], fix[MAX_READ_LENGTH] ;
	bool isStrongTrusted[MAX_READ_LENGTH] ; 
	bool isPolyAKmer[MAX_READ_LENGTH] ;
	int tstart, tend, longestTrustedLength, trustThreshold ;
	int searchResult ;
	int totalFix= 0, allowedFix ;
	int ret ;
	int fixBottleNeck, bestFixBottleNeck ;
	int strongTrustThreshold ;
	bool unfixable, forceNextRound ;

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
	allowedFix = readLength * MAX_FIX_PER_100 / 100 ;	
	unfixable = false ;

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


#ifdef DEBUG
	printf( "%s\n", seq ) ;
	for ( i = 0 ; i < kcnt ; ++i )
	{
		printf( "%d ", counts[i] ) ;
	}
	printf( "\n" ) ;
#endif

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
	if ( i >= 1 )
	//if ( i < kcnt )
	{
		//trustThreshold = iBuffer[( kcnt - 1 + i ) / 2] / 20 + 1 ;
		trustThreshold = GetBound( iBuffer[i] ) ;
		strongTrustThreshold = iBuffer[i] ;
		//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[i] ) ;
	}
	else
	{
		trustThreshold = GetBound( iBuffer[kcnt / 2] ) ;
		//trustThreshold = iBuffer[kcnt / 2] - 4 * sqrt( iBuffer[kcnt / 2 ] ) ;
		strongTrustThreshold = iBuffer[kcnt / 2] ;
	}


	for ( i = 0 ; i < kcnt ; ++i )
	{
		if ( IsPolyA( seq + i, kmerLength ) )
			isPolyAKmer[i] = true ;
		else
			isPolyAKmer[i] = false ;
	}

#ifdef DEBUG
	printf( "strong trust threshold=%d threshold=%d\n", strongTrustThreshold, trustThreshold ) ; fflush( stdout ) ;
#endif
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
		allowedFix = readLength * MAX_FIX_PER_100 / 100 ;
		totalFix = 0 ;
		unfixable = false ;
		forceNextRound = false ;
		// Find the longest trusted portion
		longestTrustedLength = -1 ;
		j = 0 ;
		memset( isStrongTrusted, false, sizeof( bool ) * readLength ) ;
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
					for ( k = i - j ; k <= i - 1 + kmerLength - 1 ; ++k )
					{
						isStrongTrusted[k] = true ;
					}
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
			for ( k = i - j ; k <= i - 1 + kmerLength - 1 ; ++k )
			{
				isStrongTrusted[k] = true ;
			}
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

#ifdef DEBUG
		for ( i = 0 ; i < readLength ; ++i )
			printf( "%d", isStrongTrusted[i] ) ;
		printf( "\n" ) ;
#endif
		//printf( "trusted region:%d %d: %c %d\n", tstart, tend, seq[tstart], longestTrustedLength ) ;


		//if ( strongTrustThreshold >= 20 && longestTrustedLength >= 5 )
		//	printf( "AS events?\n" ) ;

		//return 0 ;

		if ( longestTrustedLength  == -1 )
			return -1 ;

		if ( longestTrustedLength == kcnt )
			return 0 ;

		for ( i = 0 ; i < readLength ; ++i )
			fix[i] = -1 ;

		// Scan towards right
		kcode.Restart() ;
		if ( tend < kcnt - 1 )
		{
			int bestFixCnt = 0 ;
			int trialCnt = 0 ;
			int maxFixCnt = allowedFix ;
			for ( i = tend ; i < tend + kmerLength ; ++i )
				kcode.Append( seq[i] ) ;
			fixBottleNeck = 1000000000 ;
			bestFixBottleNeck = -1 ; 
			SearchPaths_Right( tend + kmerLength, tend + kmerLength, iter==0?-1:trustThreshold, seq, 0, iBuffer, fixBottleNeck, 
					maxFixCnt, fix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, kcode, kmers, trialCnt ) ;

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
			if ( trialCnt > MAX_TRIAL )
			{
				return -1 ;
			}
			if ( bestFixCnt >= 2 )
				return -1 ;
			else if ( bestFixCnt <= 0 )
				unfixable = true ;
			totalFix += maxFixCnt ;
		}
		//printf( "%d\n", totalFix ) ;	
		kcode.Restart() ;
		if ( tstart > 0 && !unfixable )
		{
			int bestFixCnt = 0 ;
			int maxFixCnt = allowedFix - totalFix ;
			int trialCnt = 0 ;
			for ( i = tstart ; i < tstart + kmerLength ; ++i )
				kcode.Append( seq[i] ) ;
			fixBottleNeck = 1000000000 ;
			bestFixBottleNeck = -1 ; 
			SearchPaths_Left( tstart - 1, tstart - 1, iter==0?-1:trustThreshold, seq, 0, iBuffer, fixBottleNeck, 
					maxFixCnt, fix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, kcode, kmers, trialCnt ) ;
			if ( bestFixBottleNeck != -1 && iter == 0 && bestFixBottleNeck < GetBound( strongTrustThreshold ) )
			{
				/*int t = bestFixBottleNeck ;
				fixBottleNeck = 1000000000 ;
				bestFixBottleNeck = -1 ; 
				for ( i = 0 ; i < tstart ; ++i )
					fix[i] = -1 ;
				bestFixCnt = 0 ;
				trialCnt = 0 ;
				maxFixCnt = allowedFix - totalFix ;
				SearchPaths_Left( tstart - 1, tstart - 1, t, seq, 0, iBuffer, fixBottleNeck, 
						maxFixCnt, fix, bestFixCnt, bestFixBottleNeck, isStrongTrusted, isPolyAKmer, kcode, kmers, trialCnt ) ;
				//bestFixCnt = 0 ;
				//return -1 ;*/
				forceNextRound = true ;
			}
			
			/*if ( maxFixCnt == 0 && bestFixBottleNeck < strongTrustThreshold && bestFixBottleNeck >= GetBound( strongTrustThreshold )
				&& tstart > kmerLength )
			{
				forceNextRound = true ;
				//strongTrustThreshold = bestFixBottleNeck ;
			}*/
			/*printf( "bestFixCnt=%d\n", bestFixCnt ) ;
			for ( i = 0 ; i < readLength ; ++i )
				printf( "%d ", fix[i] ) ;
			printf( "\n" ) ;*/
			if ( trialCnt > MAX_TRIAL )
				return -1 ;
			if ( bestFixCnt >= 2 )
				return -1 ;
			else if ( bestFixCnt <= 0 )
				unfixable = true ;
			totalFix += maxFixCnt ;

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
			//break ;
			return -1 ;
	}

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

	// If multiple errors within a kmer window did not decrease the count much,
	// then, they are not errors.
	for ( i = 1 ; i < cnt ; ++i )
	{
		if ( iBuffer[i] - iBuffer[i - 1] + 1 <= kmerLength )
		{
			if ( counts[ iBuffer[i - 1] ] > 1 &&
				counts[ iBuffer[i] ] > 1 &&
				counts[ iBuffer[i - 1] ] > counts[ iBuffer[i] ] / 2 &&
				counts[ iBuffer[i - 1] ] < 2 * counts[ iBuffer[i] ] )
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


	// If the fix are at the boundary and fixed all the position, then we should remove these fixes
	if ( totalFix > 3 && readLength > 10 )
	{
		int tmp = 0 ;
		for ( i = 0 ; i < 10 ; ++i )
			if ( fix[i] != -1 && seq[i] != 'N' )
				++tmp ; 
		if ( tmp >= 2 )
			for ( i = 0 ; i < 10 ; ++i )
			{
				if ( seq[i] != 'N' )
					fix[i] = -1 ;
			}

		tmp = 0 ;
		for ( i = readLength - 10 ; i < readLength ; ++i )
			if ( fix[i] == -1 && seq[i] != 'N' )
				++tmp ; 
		if ( tmp >= 3 )
			for ( i = readLength - 10 ; i < readLength ; ++i )
			{
				if ( seq[i] != 'N' )
					fix[i] = -1 ;
			}
	}	

	if ( totalFix > 3 )
	{
		iBuffer[0] = 0 ;
		for ( i = 0 ; i < readLength ; ++i )
			if ( seq[i] != 'N' && fix[i] != -1 )
				iBuffer[i + 1] = iBuffer[i] + 1 ;
			else
				iBuffer[i + 1] = iBuffer[i] ;

		for ( i = 0 ; i < kcnt ; ++i )
		{
			int threshold = 4 ;
			if ( i > 0.7 * readLength )
				threshold = 5 ;
			if ( iBuffer[i + kmerLength] - iBuffer[i] > threshold )
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
	return ret ;
}
