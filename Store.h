#ifndef _MOURISL_CLASS_STORE
#define _MOURISL_CLASS_STORE
/**
  The wrapper for kmers' storage.
*/
#include <stdio.h>
#include <stdint.h>

#if __cplusplus<=199711L
#include <tr1/unordered_map>
#else
#include <unordered_map>
#endif

#include "KmerCode.h"

class Store
{
private:
	uint64_t size ;
	//std::map<uint64_t, int> hash ;
#if __cplusplus<=199711L
	std::tr1::unordered_map<uint64_t, int> hash ;
#else
	std::unordered_map<uint64_t, int> hash ;
#endif

public:
	Store()
	{
	}

	Store( uint64_t cnt )
	{
#if __cplusplus>199711L
		hash.reserve( cnt ) ;
#endif
	}

	~Store() 
	{
	}

	void Allocate( uint64_t cnt )
	{
#if __cplusplus>199711L
		hash.reserve( cnt ) ;
#endif
	}

	int Put( KmerCode &code, int cnt  ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		hash[ code.GetCanonicalKmerCode() ] = cnt ;
		return 0 ;
	}

	int GetCount( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		if ( hash.count( code.GetCanonicalKmerCode() ) == 0 )
			return 0 ;
		return hash[ code.GetCanonicalKmerCode() ] ;
	}


	uint64_t GetCanonicalKmerCode( uint64_t code, int k ) 
	{
		int i ;
		uint64_t crCode = 0ull ; // complementary code
		for ( i = 0 ; i < k ; ++i )
		{
			//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
			//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;

			uint64_t tmp = ( code >> ( 2ull * i ) ) & 3ull ;
			crCode = ( crCode << 2ull ) | ( 3ull - tmp ) ;
		}
		return crCode < code ? crCode : code ;
	}

	int Clear() 
	{
		return 0 ;
	}
} ;
#endif
