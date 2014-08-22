#ifndef _MOURISL_KMERINFO
#define _MOURISL_KMERINFO

#include <map>
#include <vector>
#include <stdlib.h>

#include "KmerCode.h"
#include "utils.h"

struct _kmerInRead
{
	std::size_t readIdx ;
	int kmerPos ;
	bool rc ; // need reverse and complementary
} ;

class KmerInfo
{
private:
	std::map<uint64_t, std::vector<struct _kmerInRead> > kmerInfo ;
	std::map<uint64_t, std::size_t> kmerCnt ;
	int kmerLength ;
public:	
	KmerInfo( int kl ):kmerLength( kl ) 
	{
	} 

	~KmerInfo()
	{
	}

	int GetKmerLength()
	{
		return kmerLength ;
	}

	void AddKmerCount( KmerCode &kmer )
	{
		uint64_t code = kmer.GetCanonicalCode() ;
		if ( kmerCnt.count( code ) == 0 )
			kmerCnt[ code ] = 1 ;
		else
		{
			++kmerCnt[ code ] ;
		}
	}

	void AddReadInfo( KmerCode &kmer, std::size_t readIdx, int kmerPos )
	{
		struct _kmerInRead info ;
		uint64_t code = kmer.GetCanonicalCode() ;
		info.readIdx = readIdx ;
		info.kmerPos = kmerPos ;
		info.rc = ( code != kmer.GetCode() ) ;

		kmerInfo[ code ].push_back( info ) ;
	}

	
	/*int GetKmerInfo( KmerCode &kmer, struct _kmerInRead *buffer, int bufferSize )
	{
	}*/
	int GetReadsInfo( KmerCode &kmer, struct _kmerInRead **buffer )
	{
		uint64_t code = kmer.GetCanonicalCode() ;
		int i, size ;
		
		std::vector<struct _kmerInRead> &infoList = kmerInfo[code] ;
		size = infoList.size() ;
		*buffer = ( struct _kmerInRead * )malloc( sizeof( struct _kmerInRead ) * size ) ;
		for ( i = 0 ; i < size ; ++i )
			(*buffer)[i] = infoList[i] ;
		return size ;
	}

	int GetKmerCount( KmerCode &kmer ) 
	{
		uint64_t code = kmer.GetCanonicalCode() ;
		if ( kmerCnt.count( code ) == 0 )
			return 0 ;
		else
			return kmerCnt[code] ;
	}
} ;

#endif
