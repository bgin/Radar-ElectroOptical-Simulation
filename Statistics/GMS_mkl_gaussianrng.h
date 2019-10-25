
#ifndef __GMS_MKL_GAUSSIANRNG_H__
#define __GMS_MKL_GAUSSIANRNG_H__

namespace file_info {

#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

const unsigned int gGMS_MKL_GAUSSIANRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_GAUSSIANRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_MKL_GAUSSIANRNG_MICRO = gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_GAUSSIANRNG_FULLVER = 
	1000U*gGMS_MKL_GAUSSIANRNG_MAJOR + 100U*gGMS_MKL_GAUSSIANRNG_MINOR + 10U*gGMS_MKL_GAUSSIANRNG_MICRO;

const char * const pgGMS_MKL_GAUSSIANRNG_CREATE_DATE = "22-04-2018 12:16 +00200 (SUN 22 APR 2018 GMT+2)";

const char * const pgGMS_MKL_GAUSSIANRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_GAUSSIANRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_GAUSSIANRNG_DESCRIPT = "C++ wrappers around Intel MKL random number generator procedures.";

}


#include <iostream>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#else
    #error Unsupported Operating System (Linux and Win64 -- only supported currently)
#endif
#if (USE_MKL) == 1
#include <mkl_vsl.h>
#else
#error Intel MKL library required for this compilation
#endif

namespace gms {
	namespace math {
		namespace stat {
			
				// Wrapper around MKL vsrnggaussian( method, stream, n, r, a, sigma )
#if defined _WIN64
			__declspec(align(64))struct MKLGaussianRNG{
#elif defined __linux
			__attribute__((align(64))) struct MKLGaussianRNG{
#endif
#if defined _WIN64			  
				    _Field_size_(m_nvalues) double * __restrict m_rvec;
#elif defined __linux
			                double * __restrict m_rvec;
#endif
					double m_a;

					double m_sigma;

#if (USE_STRUCT_PADDING) == 1
					PAD_TO(1,8)
#endif

					MKL_INT m_nvalues;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(2,4)
#endif
					MKL_UINT m_brng;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(3,4)
#endif
					MKL_INT  m_seed;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(4,4)
#endif
					MKL_INT  m_error;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(5,4)
#endif
					
					
					

					//
					//	Construction and Destruction
					//
					MKLGaussianRNG();

					MKLGaussianRNG(const MKL_INT, 
						       const MKL_UINT,
						       const MKL_INT,
						       const double,
						       const double );

					MKLGaussianRNG(const MKLGaussianRNG &);

					MKLGaussianRNG(MKLGaussianRNG &&) noexcept(true);

					~MKLGaussianRNG() noexcept(true);

					MKLGaussianRNG & operator=(const MKLGaussianRNG &);

					MKLGaussianRNG & operator=(const MKLGaussianRNG & );

					MKLGaussianRNG & operator=(MKLGaussianRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT); // Should be an interface in case of static polymorphism.

					void compute_rand_distribution(VSLStreamStatePtr, const MKL_INT);
			};

			std::ostream & operator<<(std::ostream &,
						   const MKLGaussianRNG &);




		}
	}
}


#endif /*__GMS_MKL_RANDOM_H__*/
