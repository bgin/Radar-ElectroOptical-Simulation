
#ifndef __GMS_MKL_WIEBULLRNG_H__
#define __GMS_MKL_WEIBULLRNG_H__

namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif
  
const unsigned int gGMS_MKL_WEIBULLRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_WEIBULLRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_MKL_WEIBULLRNG_MICRO = gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_WEIBULLRNG_FULLER = 
	1000U*gGMS_MKL_WEIBULLRNG_MAJOR + 100U*gGMS_MKL_WEIBULLRNG_MINOR + 10U*gGMS_MKL_WEIBULLRNG_MICRO;

const char * const pgGMS_MKL_WEIBULLRNG_CREATE_DATE = "26-04-2018 14:17 +00200 (THR 26 APR 2018 GMT+2)";

const char * const pgGMS_MKL_WEIBULLRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_WEIBULLRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_WEIBULLRNG_DESCRIPT = "C++ wrapper for Intel MKL vdRngWeibull procedure.";

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
    #error Intel MKL library is required for this translation unit compilation
#endif

namespace gms {
	namespace math {
		namespace stat {

			//
			//	C++ wrapper for Intel MKL vdRngWeibull procedure
			//
#if defined _WIN64
			__declspec(align(64)) struct MKLWRNGData {
#elif defined __linux
                        __attribute__(align((64))) struct MKLWRNGData {
#endif
#if defined _WIN64
				_Field_size_(m_nvalues)double * __restrict m_rvec;
#elif defined __linux
                                   double * __restrict m_rvec;
#endif
				   double m_alpha;
				   double m_a;
				   double m_beta;
				   MKL_INT m_nvalues;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(1, 4)
#endif
				   MKL_INT m_brng;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(2, 4)
#endif
			       MKL_INT m_seed;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(3, 4)
#endif
				   MKL_INT m_error;
#if (USE_STRUCT_PADDING) == 1
				   PAD_TO(4, 4)
#endif
			};

			

				   

				  
					
#if defined _WIN64
			__declspec(align(64)) struct MKLWeibullRNG {
#elif defined
                        __attribute__((align(64))) struct MKLWeibullRNG {
#endif
				MKLWRNGData data;
					//
					//	Constuction and destruction
					//

					MKLWeibullRNG();

					MKLWeibullRNG(const MKL_INT,
					              const MKL_UINT,
						      const MKL_INT,
						      const double,
						      const double,
					              const double);

					MKLWeibullRNG(const MKLWeibullRNG &);

					MKLWeibullRNG(MKLWeibullRNG &&)noexcept(true);

					~MKLWeibullRNG()noexcept(true);

					MKLWeibullRNG & operator=(const MKLWeibullRNG &);

					MKLWeibullRNG & operator=(MKLWeibullRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,
			           const MKLWeibullRNG &);


		}
	}
}



#endif /*__GMS_WEIBULLRNG_H__*/
