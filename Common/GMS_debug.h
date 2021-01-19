
#ifndef __GMS_DEBUG_H__
#define __GMS_DEBUG_H__

namespace file_info {

	const unsigned int gGMS_DEBUG_MAJOR = 1U;

	const unsigned int gGMS_DEBUG_MINOR = 0U;

	const unsigned int gGMS_DEBUG_MICRO = 0U;

	const unsigned int gGMS_DEBUG_FULLVER = 
		1000U * gGMS_DEBUG_MAJOR + 100U * gGMS_DEBUG_MINOR + 10U * gGMS_DEBUG_MICRO;

	const char * const pgGMS_DEBUG_CREATE_DATE = "27-09-2019 20:30 +00200 (FRI 27 SEP 2019 GMT+2)";

	const char * const pgGMS_DEBUG_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_DEBUG_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_DEBUG_SYNOPSIS = "Debug check memory allocation functions.";
}

#include <cstdint>

#include "GMS_config.h"

namespace gms {
	namespace common {

#if (GMS_DEBUG_ON) == 1
        

		

	

		
                                            

	  	
		double * gms_dmallocu_dbg(const size_t);			 

		float * gms_fmallocu_dbg(const size_t);

		int32_t * gms_imallocu_dbg(const size_t);

		double * gms_dmalloca_dbg(const size_t,const size_t);

		float  * gms_fmalloca_dbg(const size_t,const size_t);

		int32_t * gms_imalloca_dbg(const size_t,const size_t);

		double * gms_edmalloca_dbg(const size_t,const size_t);

		float * gms_efmalloca_dbg(const size_t,const size_t);

		int32_t * gms_eimalloca_dbg(const size_t,const size_t);
					 
		
       
		


#else
#error "Include this file only in debug mode!!"
#endif
	}
}





#endif /*__LAM_DEBUG_H__*/
