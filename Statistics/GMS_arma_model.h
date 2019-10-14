
#ifndef __GMS_ARMA_MODEL_H__
#define __GMS_ARMA_MODEL_H__



namespace file_info {
#if defined _WIN64  
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

	const unsigned int gGMS_ARMA_MODEL_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_ARMA_MODEL_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_ARMA_MODEL_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_ARMA_MODEL_FULLVER = 
	 1000U*gGMS_ARMA_MODEL_MAJOR+100U*gGMS_ARMA_MODEL_MINOR+10U*gGMS_ARMA_MODEL_MICRO;

	const char * const pgGMS_ARMA_MODEL_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_ARMA_MODEL_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_ARMA_MODEL_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_ARMA_MODEL_SYNOPSIS = "AUTOMATIC AR-MA MODEL FITTING -- SCALAR CASE"
}

#if defined _WIN64
  #include "../GMS_config.h"
#elif defined __linux
  #include "GMS_config.h"
#endif
#if defined (__INTEL_COMPILER)
#include <../perf_headers/c++/valarray>
#else
#include <valarray>
#endif

#include <cstdint>

namespace gms {
	namespace math {
		namespace stat {
#if defined _WIN64
		  __declspec(align(64))	  struct ArmaModel {
#elif defined __linux
		  __attribute__((align(64))) struct ArmaModel {
#endif
				using VAf64 = std::valarray<double>;
				using VAi32 = std::valarray<int32_t>;
			
				/*
				Low-level computational code -- members are public.
				*/
			
				VAi32 *   __restrict    m_iq;

				VAf64 *   __restrict    m_b2;

				VAi32 *   __restrict    m_ip;

				VAf64 *   __restrict    m_a2;

				VAf64 *   __restrict    m_std;

				VAf64 *   __restrict    m_cxx2;

				VAf64 *   __restrict    m_g;

				VAf64 *   __restrict    m_saic;

				double   m_aicm;

#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,32)
#endif

				int32_t  m_newn;

				int32_t  m_nmax;

				int32_t  m_mmax;

				int32_t  m_kq;

				int32_t  m_kp;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2,4)
#endif
				

				

				ArmaModel();

				ArmaModel(const int32_t,
					  const int32_t,
					  const int32_t,
					  const double);

				ArmaModel(const ArmaModel &);

				ArmaModel(ArmaModel &&) noexcept(true);

				~ArmaModel() = default;

				ArmaModel & operator=(const ArmaModel &);

				ArmaModel & operator=(ArmaModel &&) noexcept(true);

				void compute_arma_model(int32_t,
							int32_t,
							VAf64 &,
							int32_t,
							VAi32 &,
							VAf64 &,
							VAi32 &,
							VAf64 &,
							int32_t,
							int32_t[6] );


			};
		}
	}
}


#endif  /*__GMS_ARMA_MODEL_H__*/
