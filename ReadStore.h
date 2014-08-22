#ifndef _MOURISL_READSTORE
#define _MOURISL_READSTORE

// The interface to store the reads information

#include <vector>
#include <string>

class ReadStore
{
private:
	std::vector<std::string> store ;
public:
	ReadStore()
	{
	}

	~ReadStore() 
	{
	}

	// @return: the index for this read
	int AddRead( char *seq )
	{
		std::string s(seq) ;
		store.push_back(s) ;
		
		return store.size() - 1 ;
	}

	const char *GetRead( int idx )
	{
		return store[idx].c_str() ;
	}
} ;

#endif
