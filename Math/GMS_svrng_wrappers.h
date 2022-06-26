
#ifndef __GMS_SVRNG_WRAPPERS_H__
#define __GMS_SVRNG_WRAPPERS_H__


// Compilation of this translation unit require an Intel Compiler
#if !defined __ICC || !defined __INTEL_COMPILER
    #error "Compilation of this translation unit requires an Intel Compiler and implementation of SVRNG library."
#endif

namespace file_info {

     const unsigned int gGMS_SVRNG_WRAPPERS_MAJOR = 1;
     const unsigned int gGMS_SVRNG_WRAPPERS_MINOR = 0;
     const unsigned int gGMS_SVRNG_WRAPPERS_MICRO = 0;
     const unsigned int gGMS_SVRNG_WRAPPERS_FULLVER = 1000U*gGMS_SVRNG_WRAPPERS_MAJOR+100U*gGMS_SVRNG_WRAPPERS_MINOR+
                                                      10U*gGMS_SVRNG_WRAPPERS_MICRO;
     const char * const pgGMS_SVRNG_WRAPPERS_CREATION_DATE = "29-12-2019 11:55 +00200 (SUN 29 DEC 2019 GMT+2)";
     const char * const pgGMS_SVRNG_WRAPPERS_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const pgGMS_SVRNG_WRAPPERS_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
     const char * const pgGMS_SVRNG_WRAPPERS_SYNOPSYS      = "GMS wrappers around Intel vector random number generators library SVRNG.";
}

#include <cstdint>
#include "GMS_config.h"
#include "GMS_avxvecf32.h"
#include "GMS_avx512c4f32.h"
#include "GMS_avxc8f32.h"
#include "GMS_avxc4f64.h"

namespace gms{

        namespace math {

                  namespace stat {

                      // mt19937 engine float8
		      void svrng_wrapper_mt19937_init_float8(float * __restrict __ATTR_ALIGN__(64),
		                                             const int64_t ,
							     const float,
							     const float,
							     const int32_t,
							     int32_t &  ) __ATTR_COLD__ __ATTR_ALIGN__(32);
		      // mt19937 engine double4
		      void svrng_wrapper_mt19937_init_double4(double * __restrict __ATTR_ALIGN__(64),
		                                              const int64_t,
							      const double,
							      const double,
							      const int32_t,
							      int32_t & ) __ATTR_COLD__ __ATTR_ALIGN__(32);

		      // mt19937 engine AVXVec8
		      void svrng_wrapper_mt19937_init_avxvec8(AVXVec8 * __restrict __ATTR_ALIGN__(64),
		                                              const int64_t,
							      const float,
							      const float,
							      const int32_t,
							      int32_t & ) __ATTR_COLD__ __ATTR_ALIGN__(32);

		      // mt19937 engine AVX512c4f32
		      void svnrg_wrapper_mt19937_init_avx512c4f32(AVX512c4f32 * __restrict __ATTR_ALIGN__(64),
		                                                  const int64_t,
								  const float,
								  const float,
								  const float,
								  const float,
								  const int32_t,
								  int32_t & ) __ATTR_COLD__ __ATTR_ALIGN__(32);
		      // mt19937 engine AVXc8f32
		      void svrng_wrapper_mt19937_init_avxc8f32(AVXc8f32 * __restrict __ATTR_ALIGN__(64),
							       const int64_t,
							       const float,
							       const float,
							       const float,
							       const float,
							       const int32_t,
							       int32_t & ) __ATTR_COLD__ __ATTR_ALIGN__(32);

		      // mt19937 engine AVXc4f64
		      void svrng_wrapper_mt19937_init_avxc4f64(AVXc4f64 * __restrict __ATTR_ALIGN__(64),
							       const int64_t,
							       const double,
							       const double,
							       const double,
							       const double,
							       const int32_t,
							       int32_t & ) __ATTR_COLD__ __ATTR_ALIGN__(32);
 
	 }  // stat
    }  // math
 
} // gms








#endif /*__GMS_SVRNG_WRAPPERS_H__*/
