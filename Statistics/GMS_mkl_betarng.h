
#ifndef __GMS_MKL_BETARNG_H__
#define __GMS_MKL_BETARNG_H__

namespace file_info {
 // Include master version file
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif
const unsigned int gGMS_MKL_BETARNG_MAJOR = lam::common::gVersionInfo.m_VersionMajor;

const unsigned int gGMS_MKL_BETARNG_MINOR = lam::common::gVersionInfo.m_VersionMinor;

const unsigned int gGMS_MKL_BETARNG_MICRO = lam::common::gVersionInfo.m_VersionMicro;

const unsigned int gGMS_MKL_BETARNG_FULLVER = 
	1000U*gGMS_MKL_BETARNG_MAJOR + 100U*gGMS_MKL_BETARNG_MINOR + 10U*gGMS_MKL_BETARNG_MICRO;



const char * const pgGMS_MKL_BETARNG_CREATE_DATE = "25-04-2018 15:37 +00200 (WED 25 APR 2018 GMT+2)";

const char * const pgGMS_MKL_BETARNG_BUILD_DATE = "00-00-0000 00:00";

const char * const pgGMS_MKL_BETARNG_AUTHOR = "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";

const char * const pgGMS_MKL_BETARNG_DESCRIPT = "C++ wrapper around Intel MKL vdRngBeta procedure.";

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
				//	C++ wrapper for Intel MKL vdRngBeta procedure.
				//
#if defined _WIN64
			__declspec(align(64)) struct MKLBRNGData{
#elif defined __linux
			__attribute__((align(64))) struct MKLBRNGData {
#endif
			         // POD and trivial type.
				// Payload aligned on 8-byte and starting at 64-byte boundary.
#if defined _WIN64
				_Field_size_(m_nvalues) double * __restrict m_rvec;
#elif defined __linux
			        double * __restrict m_rvec;
#endif
				double m_p;
				double m_q;
				double m_a;
				double m_beta;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,8)
#endif
				MKL_INT m_nvalues;
				MKL_INT m_brng;
				MKL_INT m_seed;
				MKL_INT m_error;
			};
#if defined _WIN64 
			__declspec(align(64)) struct  MKLBetaRNG { // Non POD type and non trivial type.
#elif defined __linux
			  __attribute__((align(64))) struct MKLBetaRNG {
#endif
				  

					
					MKLBRNGData datum;
					

					//
					//	Construction and destruction
					//

					MKLBetaRNG();

					MKLBetaRNG(const MKL_INT,
						   const MKL_UINT,
						   const MKL_INT,
					           const double,
						   const double,
						   const double,
						   const double);

					MKLBetaRNG(const MKLBetaRNG &);

					MKLBetaRNG(MKLBetaRNG &&)noexcept(true);

					~MKLBetaRNG()noexcept(true);

					MKLBetaRNG & operator=(const MKLBetaRNG &);

					MKLBetaRNG & operator=(MKLBetaRNG &&)noexcept(true);

					void compute_rand_distribution(const MKL_INT);

					void compute_rand_distribution(VSLStreamStatePtr,const MKL_INT);

			};

			std::ostream &
			operator<<(std::ostream &,
				   const MKLBetaRNG &);
		}
	}
}


#endif /*__GMS_MKL_BETARNG_H_*/
