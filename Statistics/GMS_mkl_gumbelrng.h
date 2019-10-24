
#ifndef __GMS_MKL_GUMBELRNG_H__
#define __GMS_MKL_GUMBELRNG_H__

namespace file_info {

#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif

const unsigned int gGMS_MKL_GUMBELRNG_MAJOR =   gms::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_GUMBELRNG_MINOR =   gms::common::gVersionInfo.m_VersionMinor;
  
const unsigned int gGMS_MKL_GUMBELRNG_MICRO =   gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_GUMBELRNG_FULLVER = 
	1000U * gGMS_MKL_GUMBELRNG_MAJOR + 100U * gGMS_MKL_GUMBELRNG_MINOR + 10U * gGMS_MKL_GUMBELRNG_MICRO;

const char * const pgGMS_MKL_GUMBELRNG_CREATE_DATE = "26-04-2018 10:47 +00200 (THR 26 APR 2018 GMT+2)";

const char * const pgGMS_MKL_GUMBELRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_GUMBELRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_GUMBELRNG_DESCRIPT = "C++ wrapper around Intel MKL vdRngGumbel procedure.";

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

			//
			//	C++ wrapper for Intel MKL vdRngGumbel procedure.
			//
#if defined _WIN64
			__declspec(align(64)) struct  MKLGumbelRNG {
#elif defined __linux
			__attribute__((align(64))) struct MKLGumbelRNG {
#endif
#if defined _WIN64
				       _Field_size_(m_nvalues)double * __restrict m_rvec;
#elif defined __linux
                                        double * __restrict m_rvec;
#endif
					double      m_a;

					double      m_beta;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(1,8)
#endif

					MKL_INT		m_nvalues;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(2,4)
#endif
					MKL_UINT    m_brng;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(3,4)
#endif
					MKL_INT     m_seed;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(4,4)
#endif
					MKL_INT     m_error;
#if (USE_STRUCT_PADDING) == 1
					PAD_TO(5,4)
#endif
					

				

					//
					//	Construction and destruction
					//

					MKLGumbelRNG();

					MKLGumbelRNG(const MKL_INT,
						     const MKL_UINT,
						     const MKL_INT,
						     const double,
					             const double);

					MKLGumbelRNG(const MKLGumbelRNG &);

					MKLGumbelRNG(MKLGumbelRNG &&)noexcept(true);

					~MKLGumbelRNG()noexcept(true);

					MKLGumbelRNG & operator=(const MKLGumbelRNG &);

					MKLGumbelRNG & operator=(MKLGumbelRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr, const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,
			          const MKLGumbelRNG &);
		}
	}
}



#endif /*__GMS_MKL_GUMBELRNG_H__*/
