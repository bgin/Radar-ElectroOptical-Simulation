
#ifndef __GMS_MKL_UNIFORMRNG_H__
#define __GMS_MKL_UNIFORMRNG_H__

namespace file_info {

#if defined _WIN64
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif
  
const unsigned int gGMS_MKL_UNIFORMRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_UNIFORMRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_MKL_UNIFORMRNG_MICRO = gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_UNIFORMRNG_FULLVER = 
	1000U*gGMS_MKL_UNIFORMRNG_MAJOR + 100U*gGMS_MKL_UNIFORMRNG_MINOR + 10U*gGMS_MKL_UNIFORMRNG_MICRO;

const char * const pgGMS_MKL_UNIFORMRNG_CREATE_DATE = "24-04-2018 08:42 +00200 (TUE 24 APR 2018 GMT+2)";

const char * const pgGMS_MKL_UNIFORMRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_UNIFORMRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_UNIFORMRNG_DESCRIPT = "C++ wrapper over MKL vdRngUniform procedure.";

}

#include <iostream>
#if defined _WIN64
    #include "../GMS_config.h"
#elif defined __linux
    #include "GMS_config.h"
#else
    #error UNsupported Operating System (Linux and Win64 are currently supported)
#endif
#if (USE_MKL) == 1
    #include <mkl_vsl.h>
#else
    #error Intel MKL library is required for compilation of this translation unit
#endif

namespace gms {
	namespace math {
		namespace stat {

			//
			//	C++ wrapper around Intel MKL vdRngUniform procedure.
			//
#if defined _WIN64
			__declspec(align(64))  struct MKLUniformRNG{
#elif defined __linux
                        __attribute__((align(64))) struct MKLUniformRNG {
#endif
#if defined _WIN64
				    _Field_size_(m_nvalues) double * __restrict m_rvec; // vector of random values
#elif defined __linux
                                    double * __restrict m_rvec;
#endif
				    double   m_a; // param 'a'

				    double   m_b;   // param 'b'

#if (USE_STRUCT_PADDING) == 1
					PAD_TO(1,8)
#endif
					
					MKL_INT m_nvalues; // number of random values in vector m_rvec
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(2,4)
#endif
					MKL_UINT m_brng;  // Basic random number generator
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(3,4)
#endif
					MKL_INT  m_seed;  // initial seed 
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(4,4)
#endif
					MKL_INT  m_error; // error returned by the MKL
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(5,4)
#endif
					

					//
					//	Constructor and destructor
					//

					MKLUniformRNG();

					MKLUniformRNG(const MKL_INT,
						      const MKL_UINT,
						      const MKL_INT,
						      const double,
					              const double);

					MKLUniformRNG(const MKLUniformRNG &);

					MKLUniformRNG(MKLUniformRNG &&)noexcept(true);

					~MKLUniformRNG()noexcept(true);

					MKLUniformRNG & operator=(const MKLUniformRNG &);

					MKLUniformRNG & operator=(MKLUniformRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,const MKLUniformRNG &);
		}
	}
}

#endif /*__GMS_MKL_UNIFORMRNG_H__*/
