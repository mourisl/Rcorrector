#ifndef _RCORRECTOR_ERRORCORRECTION_MOURISL
#define _RCORRECTOR_ERRORCORRECTION_MOURISL

#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#include "Store.h"
#include "Reads.h"
#include "utils.h"

struct _ErrorCorrectionThreadArg
{
	int kmerLength ;
	Store *kmers ;
	struct _Read *readBatch ;
	int batchSize ;
	int batchUsed ;
	
	pthread_mutex_t *lock ;
} ;

void *ErrorCorrection_Thread( void *arg )  ;
int ErrorCorrection( char *seq, char *qual, KmerCode &kcode, Store &kmers ) ;

#endif
