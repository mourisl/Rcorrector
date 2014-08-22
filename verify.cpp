// A verification program to test the correction result from mason simulator.
// Usage: ./a.out *.fa/fq [OPTIONS]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char buffer[2048] ;	
char id[2048] ;
char oriSeq[2048] ;
char seq[2048] ;
char qual[2048] ;
char origSeq[2048] ;
char cigar[2048] ;
int align[2048] ;
int lcs[501][501] ;
int choose[501][501] ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

void ReverseComplement( char *s )
{
	char buffer[2048] ;
	int len = strlen( s ), i ;
	for ( i = 0 ; i < len ; ++i )
		buffer[i] = numToNuc[ 3 - nucToNum[ s[len - i - 1] - 'A' ] ] ;
	buffer[i] = '\0' ;
	strcpy( s, buffer ) ;
}

bool StrCompWithTrim( char *ref, char *s )
{
	int i ;
	for ( i = 0 ; s[i] && ref[i] ; ++i )
		if ( s[i] != ref[i] )
			break ;
	return s[i] != '\0' ;
}

char *FindIdColumn( char *id, const char *tag ) 
{
	char *p = strstr( id, tag ) ;
	if ( p == NULL )
		return p ;
	while ( *p && *p != '=' )
		++p ;
	++p ;
	//return ( p + strlen( tag ) + 1 ) ;
	return  p ;
}

int Lcs( int i, int j, char *a, char *b  )
{
	if ( i < 0 || j < 0 )
		return 0 ;
	if ( lcs[i][j] != -1 )
		return lcs[i][j] ;
	int tmp ;
	int max = -1 ;
	tmp = ( a[i] == b[j] ) + Lcs( i - 1, j - 1, a, b ) ;
	if ( tmp > max )
	{
		max = tmp ;
		choose[i][j] = 1 ;
	}

	if ( j > i )
	{
		tmp = Lcs( i, j - 1, a, b ) ;
		if ( tmp > max )
		{
			max = tmp ;
			choose[i][j] = 0 ;
		}
	}
	
	if ( i > j )
	{
		tmp = Lcs( i - 1, j, a, b ) ;
		if ( tmp > max )
		{
			max = tmp ;
			choose[i][j] = 2 ;
		}
	}
	return lcs[i][j] = max ;
}

// Align b to a
void Alignment( int lena, char *a, int lenb, char *b  )
{
	int i, j ;
	if ( lena == lenb )
	{
		for ( i = 0 ; i < lenb ; ++i )
			align[i] = i ;
		return ;
	}
	memset( lcs, -1, sizeof( lcs ) ) ;	
	int tmp = Lcs( lena - 1, lenb - 1, a, b ) ;
	for ( j = lenb - 1, i = lena - 1 ; j >= 0 && i >= 0 ; )
	{
		if ( choose[i][j] == 1 )
		{
			align[j] = i ;
			--i ; --j ;
		}
		else if ( choose[i][j] == 0 )
		{
			align[j] = -1 ;
			--j ;
		}
		else if ( choose[i][j] == 2 )
		{
			align[j] = i - 1 ;
			--i ;
		}
	}
	while ( j >= 0 )
	{
		align[j] = -1 ;
		--j ;
	}
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int len ;
	char *p ;
	FILE *fp ;
	int FILE_TYPE ; // 0-fasta, 1-fastq
	int correctCount = 0, errorCount = 0 ;
	int sameCount = 0 ;
	int trimCount = 0, trimSum = 0 ;
	int exp ; //0-low,1-med,2-high,3-unknown
	bool verbose = false ;
	bool baseVerbose = false ;
	int baseTP[4] = {0,0,0,0}, baseFP[4] = {0,0,0,0}, baseFN[4] = {0,0,0,0} ;
	int readTP[4] = {0,0,0,0}, readFP[4] = {0,0,0,0}, readFN[4] = {0,0,0,0} ;
	bool useExp = false ;
	bool ignoreIndel = false ;
	for ( i = 2 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i],"-v" ) )
			verbose = true ;
		else if ( !strcmp( argv[i], "-bv" ) )
			baseVerbose = true ;
		else if ( !strcmp( argv[i], "-exp" ) )
			useExp = true ;
		else if ( !strcmp( argv[i], "-noindel" ) )
			ignoreIndel = true ;
	}

	// Decide whether it is FASTQ or FASTA.
	fp = fopen( argv[1], "r" ) ;
	fscanf( fp, "%s", buffer ) ;
	if ( buffer[0] == '>' )
		FILE_TYPE = 0 ;
	else
		FILE_TYPE = 1 ;
	fclose( fp ) ;

	fp = fopen( argv[1], "r" ) ;
	while ( fgets( id, sizeof( id ), fp ) != NULL )
	{
		if ( FILE_TYPE == 0 )
		{
			fgets( seq, sizeof( seq ), fp ) ;
		}
		else if ( FILE_TYPE == 1 )
		{
			fgets( seq, sizeof( seq ), fp ) ;
			fgets( buffer, sizeof( buffer ), fp ) ;
			fgets( qual, sizeof( qual ), fp ) ;
		}
		//printf( "%s%s%s", id, seq,qual ) ;
		// Clean the return symbol	
		len = strlen( id ) ;		
		if ( id[len - 1] == '\n')
			id[len - 1] = '\0' ;
		len = strlen( seq ) ;
		if ( seq[len - 1] == '\n' )
			seq[len - 1] = '\0' ;
		if ( qual[len - 1] == '\n' )
			qual[len - 1] = '\0' ;
		
		// Parse the id field
		p = FindIdColumn( id, "haplotype_infix" ) ;
		sscanf( p, "%s", origSeq ) ;
	
		if ( ignoreIndel && strlen( seq ) != strlen( origSeq ) )
			continue ;
		if ( strlen( seq ) != strlen( origSeq ) )
		{
			printf( "%s\n%s\n", id, seq ) ;
		}
		p = FindIdColumn( id, "edit_string" ) ;
		sscanf( p, "%s", cigar ) ;
		
		p = FindIdColumn( id, "strand=reverse" ) ;
		if ( p != NULL )
		{
			ReverseComplement( origSeq ) ;
		}
		
		p = FindIdColumn( id, "exp" ) ;
		if ( p != NULL )
		{
			sscanf( p, "%s", buffer ) ;
			if ( !strcmp( buffer, "high" ) )
				exp = 2 ;
			else if ( !strcmp( buffer, "medium" ) )
				exp = 1 ;
			else if ( !strcmp( buffer, "low" ) )
				exp = 0 ;
			else
				exp = 3 ;
		}
		else
			exp = 3 ;
		if ( verbose || baseVerbose )
			printf( "%s\n", id ) ;
		//printf( "%s %s\n", seq, origSeq ) ;
		if ( StrCompWithTrim( origSeq, seq ) )
		{
			/*if ( verbose )
			{
				printf( "Diff:\n%s\n%s\n", id, seq ) ;
			}*/
			++errorCount ;
			
			for ( i = 0 ; cigar[i] ; ++i )
			{
				if ( cigar[i] != 'M' )
					break ;
			}
			if ( cigar[i] )
			{
				if ( verbose )
					printf( "FN\n" ) ;
				++readFN[exp] ;
			}
			else
			{
				if ( verbose )
					printf( "FP\n" ) ;
				++readFP[exp] ;
			}
		}
		else
		{
			/*if ( verbose )
			{
				printf( "Same:\n%s\n%s\n", id, seq ) ;
			}*/
			//printf( "S\n" ) ;
			++correctCount ;
			
			for ( i = 0 ; cigar[i] ; ++i )
			{
				if ( cigar[i] != 'M' )
					break ;
			}
			if ( cigar[i] )
			{
				if ( verbose )
					printf( "TP\n" ) ;
				++readTP[exp] ;
			}
		}

		for ( i = 0 ; cigar[i] ; ++i )
		{
			if ( cigar[i] != 'M' )
				break ;
		}
		if ( !cigar[i] )
			++sameCount ;


		p = FindIdColumn( id, "trim" ) ;
		if ( p != NULL )
		{
			int tmp = atoi( p ) ;
			//printf( "%d %s\n", tmp, p ) ;
			++trimCount ;
			trimSum += tmp ;
		}

		// Collect information of TP, FP, FN for base level
		int verboseType = 0 ;
		int lena = strlen( origSeq ) ;
		int lenb = strlen( seq ) ;
		bool visited[2048] ;
		memset( visited, false, sizeof( bool ) * lena ) ;
		Alignment( lena, origSeq, lenb, seq ) ; 
		for ( i = 0 ; seq[i] ; ++i )
		{
			//printf( "(%d, %d) ", align[i], baseFP[exp] ) ;
			if ( align[i] == -1 )
			{
				//if ( i > 2 && i < lenb - 2 )
				//	printf( "%s\n%s\n", id, seq ) ;
				++baseFP[exp] ;
				continue ;
			}
			visited[ align[i] ] = true ;
			/*if ( i == 0 )
				baseFP[exp] += align[i] ;
			else 
				baseFP[exp] += ( align[i] - align[i - 1] - 1 ) ;*/

			if ( cigar[ align[i] ] == 'M' )
			{
				if ( seq[i] != origSeq[ align[i]] )
				{
					verboseType = 1 ;
					++baseFP[exp] ;
				}
			}
			else if ( cigar[align[i]] == 'E' )
			{
				if ( seq[i] == origSeq[align[i]] )
				{
					if ( verboseType == 0 )
						verboseType = 2 ;
					//printf( "%d %d\n", i, align[i] ) ;
					++baseTP[exp] ;
				}
				else
				{
					if ( verboseType == 0 || verboseType == 2 )
						verboseType = 3 ;
					++baseFN[exp] ;
				}
			}
		}
		
		for ( k = lena - 1 ; k >= 0 ; --k )
			if ( visited[k] )
				break ;

		for ( i = 0 ; i < k + 1 ; ++i )
			if ( visited[i] == false )
			{
				if ( cigar[i] == 'M' )	
					++baseFP[exp] ;
				else
					++baseFN[exp] ;
			}

		for ( i = 0 ; ; ++i )
			if ( align[i] == -1 )
				--baseFP[exp] ;
			else
				break ;
		for ( i = lenb - 1 ; ; --i )
			if ( align[i] == -1 )	
				--baseFP[exp] ;
			else
				break ;
		//printf( "\n" ) ;
		/*int verboseType = 0 ;
		for ( i = 0 ; seq[i] ; ++i )
		{
			if ( cigar[i] == 'M' )
			{
				if ( seq[i] != origSeq[i] )
				{
					verboseType = 1 ;
					++baseFP[exp] ;
				}
			}
			else if ( cigar[i] == 'E' )
			{
				if ( seq[i] == origSeq[i] )
				{
					if ( verboseType == 0 )
						verboseType = 2 ;
					++baseTP[exp] ;
				}
				else
				{
					if ( verboseType == 0 || verboseType == 2 )
						verboseType = 3 ;
					++baseFN[exp] ;
				}
			}
		}*/

		if ( baseVerbose )
		{
			if ( verboseType == 1 )
				printf( "FP\n" ) ;
			else if ( verboseType == 2 )
				printf( "TP\n" ) ;
			else if ( verboseType == 3 )
				printf( "FN\n" ) ;
		}
 	}

	int TP, FP, FN ;
	printf( "correct #: %d\n"
		"error #: %d\n", correctCount, errorCount ) ;
	printf( "Original Correct Reads Count: %d\n", sameCount ) ;
	printf( "Trimmed Reads Count: %d. Average trim length: %lf\n", trimCount, (double)trimSum / trimCount ) ;

	printf( "Overall:\n") ;
	TP = FP = FN = 0 ;
	for ( i = 0  ; i < 4 ; ++i )
	{
		TP += baseTP[i] ;
		FP += baseFP[i] ;
		FN += baseFN[i] ;
	}
	printf( "\nBase level:\n" ) ;
	printf( "TP: %d\nFP: %d\nFN: %d\n", TP, FP, FN ) ;
	double recall = ( double )TP/(TP+FN) ;
	double precision = (double)TP/(TP+FP) ;
	printf( "Recall: %lf\n"
		"Precision: %lf\n"
		"F-score: %lf\n"
		"Gain: %lf\n", recall, precision, 2*recall*precision / ( recall + precision ),
				(double)(TP-FP)/(TP+FN) ) ;
	
	TP = FP = FN = 0 ;
	for ( i = 0  ; i < 4 ; ++i )
	{
		TP += readTP[i] ;
		FP += readFP[i] ;
		FN += readFN[i] ;
	}
	printf( "\nRead level:\n" ) ;
	printf( "TP: %d\nFP: %d\nFN: %d\n", TP, FP, FN ) ;
	recall = ( double )TP/(TP+FN) ;
	precision = (double)TP/(TP+FP) ;
	printf( "Recall: %lf\n"
			"Precision: %lf\n"
			"F-score: %lf\n"
			"Gain: %lf\n", recall, precision, 2*recall*precision / ( recall + precision ),
			(double)(TP-FP)/(TP+FN) ) ;

	if ( useExp )
	{
		for ( i = 0 ; i < 3 ; ++i )
		{
			TP = baseTP[i] ;
			FP = baseFP[i] ;
			FN = baseFN[i] ;
			printf( "\nExpress level: %d\n", i ) ;
			printf( "Base level:\n" ) ;
			printf( "TP: %d\nFP: %d\nFN: %d\n", TP, FP, FN ) ;
			double recall = ( double )TP/(TP+FN) ;
			double precision = (double)TP/(TP+FP) ;
			printf( "Recall: %lf\n"
					"Precision: %lf\n"
					"F-score: %lf\n"
					"Gain: %lf\n", recall, precision, 2*recall*precision / ( recall + precision ),
					(double)(TP-FP)/(TP+FN) ) ;

			TP = readTP[i] ;
			FP = readFP[i] ;
			FN = readFN[i] ;
			printf( "\nRead level:\n" ) ;
			printf( "TP: %d\nFP: %d\nFN: %d\n", TP, FP, FN ) ;
			recall = ( double )TP/(TP+FN) ;
			precision = (double)TP/(TP+FP) ;
			printf( "Recall: %lf\n"
					"Precision: %lf\n"
					"F-score: %lf\n"
					"Gain: %lf\n", recall, precision, 2*recall*precision / ( recall + precision ),
					(double)(TP-FP)/(TP+FN) ) ;

		}
	}
	return 0 ;
}
