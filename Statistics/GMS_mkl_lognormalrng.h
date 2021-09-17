
#ifndef __GMS_MKL_LOGNORMALRNG_H__
#define __GMS_MKL_LOGNORMALRNG_H__ 20420181155

namespace file_info {


  
const unsigned int gGMS_MKL_LOGNORMALRNG_MAJOR = 1;

const unsigned int gGMS_MKL_LOGNORMALRNG_MINOR = 1;

const unsigned int gGMS_MKL_LOGNORMALRNG_MICRO = 0;

const unsigned int gGMS_MKL_LOGNORMALRNG_FULLVER = 
	1000U*gGMS_MKL_LOGNORMALRNG_MAJOR + 100U*gGMS_MKL_LOGNORMALRNG_MINOR + 10U*gGMS_MKL_LOGNORMALRNG_MICRO;

const char * const pgGMS_MKL_LOGNORMALRNG_CREATE_DATE = "25-04-2018 11:55 +00200 (WED 25 APR 2018 GMT+2)";

const char * const pgGMS_MKL_LOGNORMALRNG_BUILD_DATE = __DATE__ ":" __TIME__;

const char * const pgGMS_MKL_LOGNORMALRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_LOGNORMALRNG_DESCRIPT = "C++ wrapper around Intel MKL vdRngLognormal procedure.";

}




#include <iostream>

#include "GMS_config.h"

#if (USE_MKL) == 1
    #include <mkl_vsl.h>
#else
    #error Intel MKL library is required for the compilation
#endif

namespace gms {
	namespace math {
		namespace stat {
				
				//
				//	C++ wrapper around Intel MKL vdRngLognormal procedure.
				//

                    __attribute__((align(64)))   struct MKLLRNGData {


                                                        double * __restrict m_rvec;

							double   m_a;

				                        double   m_sigma;

				                        double   m_b;

				                        double   m_beta;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO_ALIGNED(8,1, 8)
#endif

			 	                       MKL_INT m_nvalues;

				                       MKL_UINT m_brng;

				                       MKL_INT  m_seed;

			                               MKL_INT  m_error;
			};


                        __attribute__((align(64))) struct MKLLognormalRNG {

				  
					
					
					MKLLRNGData data;

					//
					//	Construction and destruction.
					//

					MKLLognormalRNG();

					MKLLognormalRNG(const MKL_INT,
						        const MKL_UINT,
						        const MKL_INT,
						        const double,
						        const double,
						        const double,
						        const double);

					MKLLognormalRNG(const MKLLognormalRNG &);

					MKLLognormalRNG(MKLLognormalRNG &&)noexcept(true);

					~MKLLognormalRNG()noexcept(true);

					MKLLognormalRNG & operator=(const MKLLognormalRNG &);

					MKLLognormalRNG & operator=(KLLognormalRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);

			};

			std::ostream &
			operator<<(std::ostream &,
			           const MKLLognormalRNG &);
		}
	}
}




#endif /*__GMS_MKL_LOGNORMALRNG_H__*/
