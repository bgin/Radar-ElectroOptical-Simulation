
#ifndef __GMS_MKL_EXPONENTIALRNG_H__
#define __GMS_MKL_EXPONENTIALRNG_H__

namespace file_info {

#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

const unsigned int gGMS_MKL_EXPONENTIALRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_EXPONENTIALRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_MKL_EXPONENTIALRNG_MICRO = gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_EXPONENTIALRNG_FULLVER = 
	1000U*gGMS_MKL_EXPONENTIALRNG_MAJOR + 100U*gGMS_MKL_EXPONENTIALRNG_MINOR + 10U*gGMS_MKL_EXPONENTIALRNG_MICRO;

const char * const pgGMS_MKL_EXPONENTIALRNG_CREATE_DATE = "24-04-2018 15:25 +00200 (TUE 24 APR 2018 GMT+2)";

const char * const pgGMS_MKL_EXPONENTIALRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_EXPONENTIALRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_EXPONENTIALRNG_DESCRIPT = "C++ wrapper around Intel MKL vdRngExponential procedure.";

}



#include <iostream>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#endif
#if (USE_MKL) == 1
#include <mkl_vsl.h>
#endif

namespace gms {
	namespace math {
		namespace stat {

			// C++ wrapper around Intel MKL vdRngExponential procedure.
#if defined _WIN64
			__declspec(align(64)) struct MKLExponentialRNG {
#elif defined __linux
			__attribute__((align(64))) struct MKLExponentialRNG {
#endif
#if defined _WIN64
				    _Field_size_(m_nvalues) double * __restrict m_rvec;
#elif defined __linux
			                double * __restrict m_rvec;
#endif
					double   m_a;

					double   m_beta;
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
					//	Constructor and destructor
					//
					MKLExponentialRNG();

					MKLExponentialRNG(const MKL_INT,
							  const MKL_UINT,
							  const MKL_INT,
							  const double,
							  const double);

					MKLExponentialRNG(const MKLExponentialRNG &);

					MKLExponentialRNG(MKLExponentialRNG &&)noexcept(true);

					~MKLExponentialRNG()noexcept(true);

					MKLExponentialRNG & operator=(const MKLExponentialRNG &);

					MKLExponentialRNG & operator=(MKLExponentialRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);

			};

			std::ostream &
			operator<<(std::ostream &,
				   const MKLExponentialRNG &);
		}
	}
}


#endif /*__GMS_MKL_EXPONENTIALRNG_H__*/
