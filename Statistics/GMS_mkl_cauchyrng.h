
#ifndef __GMS_MKL_CAUCHYRNG_H__
#define __GMS_MKL_CAYCHYRNG_H__

namespace file_info {
#if defined _WIN64
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

  const unsigned int gGMS_MKL_CAUCHYRNG_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

  const unsigned int gGMS_MKL_CAUCHYRNG_MINOR = gms::common::gVersionInfo.m_VersionMinor;

  const unsigned int gGMS_MKL_CAUCHYRNG_MICRO =   gms::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_CAUCHYRNG_FULLVER = 
	1000U * gGMS_MKL_CAUCHYRNG_MAJOR + 100U*gGMS_MKL_CAUCHYRNG_MINOR + 10U*gGMS_MKL_CAUCHYRNG_MICRO;

const char * const pgGMS_MKL_CAUCHYRNG_CREATE_DATE = "23-04-2018 09:16 +00200 (MON 23 APR 2018 GMT+2)";

const char * const pgGMS_MKL_CAUCHYRNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_CAUCHYRNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_CAUCHYRNG_DESCRIPT = "C++ wrapper over MKL vdRngCauchy procedure.";

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

			// Wrapper around Mkl_Api(int,vdRngCauchy,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))
#if defined _WIN64
			   __declspec(align(64))struct MKLCauchyRNG {
#elif defined __linux
			     __attribute__((align(64))) MKLCauchyRNG {
#endif
#if defined _WIN64
				_Field_size_(m_nvalues) double * __restrict m_rvec; // vector of random values
#elif defined __linux
				double * __restrict m_rvec;
#endif
				double m_a;       // param 'a'

				double m_beta;   // param 'beta'
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,8)
#endif

				MKL_INT m_nvalues;  // number of random values in vector m_rvec
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2,4)
#endif
				MKL_UINT m_brng;    // Basic random number generator
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3,4)
#endif
				MKL_INT m_seed;     // initial seed
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(4,4)
#endif
				MKL_INT m_error;   // error returned by the MKL
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(5,4)
#endif
				

				

				//
				//	Construction and destruction
				//
				MKLCauchyRNG();

				MKLCauchyRNG(const MKL_INT,
					     const MKL_UINT,
					     const MKL_INT,
					     const double,
					     const double);

				MKLCauchyRNG(const MKLCauchyRNG &);

				MKLCauchyRNG(MKLCauchyRNG &&)noexcept(true);

				~MKLCauchyRNG()noexcept(true);

				MKLCauchyRNG & operator=(const MKLCauchyRNG &);

				MKLCauchyRNG & operator=(MKLCauchyRNG &&)noexcept(true);

				void compute_rand_distribution(const MKL_INT);

				void compute_rand_distribution(VSLStreamStatePtr,
							       const MKL_INT);


			};

			std::ostream &
			operator<<(std::ostream &,const MKLCauchyRNG &);
		}
	}
}



#endif /*__GMS_MKL_CAUCHYRNG_H__*/
