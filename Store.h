#ifndef _MOURISL_CLASS_STORE
#define _MOURISL_CLASS_STORE
/**
  The wrapper for kmers' storage.
*/
#include <stdio.h>
#include <stdint.h>
#include <tr1/unordered_map>

#include "KmerCode.h"

class Store
{
private:
	uint64_t size ;
	//std::map<uint64_t, int> hash ;
	std::tr1::unordered_map<uint64_t, int> hash ;
public:
	Store()
	{
	} 

	~Store() 
	{
	}

	
	int Put( KmerCode &code, int cnt  ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		hash[ code.GetCode() ] = cnt ;
		return 0 ;
	}

	int GetCount( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		if ( hash.count( code.GetCode() ) == 0 )
			return 0 ;
		return hash[ code.GetCode() ] ;
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
