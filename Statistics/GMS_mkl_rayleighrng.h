
#ifndef __GMS_MKL_RAYLEIGHRNG_H__
#define __GMS_MKL_RAYLEIGHRNG_H__

namespace file_info {

#if defined _WIN64
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif
	const unsigned int gGMS_MKL_RAYLEIGHRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_MKL_RAYLEIGHRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_MKL_RAYLEIGHRNG_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_MKL_RAYLEIGHRNG_FULLVER = 
	     1000U * gGMS_MKL_RAYLEIGHRNG_MAJOR + 100U*gGMS_MKL_RAYLEIGHRNG_MINOR + 10U*gGMS_MKL_RAYLEIGHRNG_MICRO;


	const char * const pgGMS_MKL_RAYLEIGHRNG_CREATE_DATE = "27-04-2018 09:31 +00200 (FRI 27 APR 2018 GMT+2)";

	const char * const pgGMS_MKL_RAYLEIGHRNG_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_MKL_RAYLEIGHRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

	const char * const pgGMS_MKL_RAYLEIGHRNG_DESCRIPT = "C++ wrapper for Intel MKL vdRngRayleigh procedure.";

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
    #error Intel MKL library is requaired for this compilation
#endif

namespace gms {
	namespace math {
		namespace stat {

			//
			//	C++ wrapper for Intel MKL vdRngRayleigh procedure.
			//
#if defined _WIN64
			__declspec(align(64)) struct MKLRNGData {
#elif defined __linux
                        __attribute__((align(64))) struct MKLRNGData {
#endif
#if defined _WIN64
				_Field_size_(m_nvalues)double * __restrict m_rvec;
#elif defined __linux
                                double * __restrict m_rvec;
#endif
				double      m_a;

				double      m_beta;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1, 8)
#endif

					MKL_INT		m_nvalues;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2, 4)
#endif
					MKL_UINT    m_brng;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3, 4)
#endif
					MKL_INT     m_seed;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(4, 4)
#endif
					MKL_INT     m_error;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(5, 4)
#endif
				
			};
#if defined _WIN64
			__declspec(align(64))  struct MKLRayleighRNG{

#elif defined __linux
                        __attribute__((align(64))) struct MKLRayleighRNG{
#endif
					
					MKLRNGData datum;
					

					//
					//	Construction and destruction
					//

					MKLRayleighRNG();

					MKLRayleighRNG(const MKL_INT,
						       const MKL_UINT,
						       const MKL_INT,
						       const double,
						       const double);

					MKLRayleighRNG(const MKLRayleighRNG &);

					MKLRayleighRNG(MKLRayleighRNG &&)noexcept(true);

					~MKLRayleighRNG()noexcept(true);

					MKLRayleighRNG & operator=(const MKLRayleighRNG &);

					MKLRayleighRNG & operator=(MKLRayleighRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);
			};

			std::ostream &
			operator<<(std::ostream &,
				    const MKLRayleighRNG &);
		}
	}
}



#endif /*__GMS_MKL_RAYLEIGHRNG_H__*/
