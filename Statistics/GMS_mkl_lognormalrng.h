
#ifndef __GMS_MKL_LOGNORMALRNG_H__
#define __GMS_MKL_LOGNORMALRNG_H__

namespace file_info {

#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif
  
const unsigned int gGMS_MKL_LOGNORMALRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_LOGNORMALRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_MKL_LOGNORMALRNG_MICRO = gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_LOGNORMALRNG_FULLVER = 
	1000U*gGMS_MKL_LOGNORMALRNG_MAJOR + 100U*gGMS_MKL_LOGNORMALRNG_MINOR + 10U*gGMS_MKL_LOGNORMALRNG_MICRO;

const char * const pgGMS_MKL_LOGNORMALRNG_CREATE_DATE = "25-04-2018 11:55 +00200 (WED 25 APR 2018 GMT+2)";

const char * const pgGMS_MKL_LOGNORMALRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_LOGNORMALRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_LOGNORMALRNG_DESCRIPT = "C++ wrapper around Intel MKL vdRngLognormal procedure.";

}




#include <iostream>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#else
    #error Unsupported Operating System (Linux and Win64 are currently supported)
#endif
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
#if defined _WIN64
		    __declspec(align(64))	 struct MKLLRNGData {
#elif defined __linux
                    __attribute__((align(64)))   struct MKLLRNGData {
#endif
#if defined _WIN64
				_Field_size_(m_nvalues) double * __restrict m_rvec;
#elif defined __linux
                                                        double * __restrict m_rvec;
#endif
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

#if defined _WIN64				
			__declspec(align(64)) struct  MKLLognormalRNG{
#elif defined __linux
                        __attribute__((align(64))) struct MKLLognormalRNG {
#endif
				  
					
					
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
