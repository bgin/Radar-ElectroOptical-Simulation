
#ifndef __GMS_BAYESIAN_ESTIMATE_H__
#define __GMS_BAYESIAN_ESTIMATE_H__



namespace file_info {
#if defined _WIN64  
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

	const unsigned int gGMS_BASYESIAN_ESTIMATE_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_BAYESIAN_ESTIMATE_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_BAYESIAN_ESTIMATE_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_BAYESIAN_ESTIMATE_FULLVER = 
	 1000U*gGMS_BAYESIAN_ESTIMATE_MAJOR+100U*gGMS_BAYESIAN_ESTIMATE_MINOR+10U*gGMS_BAYESIAN_ESTIMATE_MICRO;

	const char * const pgGMS_BAYESIAN_ESTIMATE_CREATE_DATE = "06-10-2019 14:31 +00200 (SUN 06 OCT 2019 GMT+2)";

	const char * const pgGMS_BAYESIAN_ESTIMATE_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_BAYESIAN_ESTIMATE_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_BAYESIAN_ESTIMATE_SYNOPSIS =  " BAYESIAN ESTIMATES OF TIME SERIES MODELS ";
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
		  __declspec(align(64))	struct BayesEstimate {
#elif defined __linux
		  __attribute__((align(64))) struct BayesianEstimate {
#endif
				using VAf64 = std::valarray<double>;
				using VAi32 = std::valarray<int32_t>;

				/*
					Low-level computational code -- members are public.
				*/

				int32_t m_k;

				int32_t m_n;

				int32_t m_il;

				int32_t m_mj2;

				int32_t m_m;

				double m_zmean;

				double m_sum;

				double m_aicm;

				double m_sdm;

				double m_aicb;

				double m_sdb;

				double m_ek;

				double m_oeic;

				double m_omean;

				double m_om;

				VAi32  m_ind;

				VAf64  m_a1;

				VAf64  m_sd;

				VAf64  m_aic;

				VAf64  m_dic;

				VAf64  m_a2;

				VAf64  m_c;

				VAf64  m_c1;

				VAf64  m_c2;

				VAf64  m_b;

				VAf64  m_esum;

				VAf64  m_e;

				VAf64  m_emean;

				VAf64  m_vari;

				VAf64  m_skew;

				VAf64  m_peak;

				VAf64  m_cov;

				VAf64  m_sxx;

				BayesEstimate() = default;

				BayesEstimate(const int32_t ,
					      const int32_t ,
					      const int32_t ,
					      const int32_t ,
					      const int32_t ,
					      const double );

				BayesEstimate(const BayesEstimate &);

				BayesEstimate(BayesEstimate &&) noexcept(true);

				~BayesEstimate() = default;

				BayesEstimate & operator=(const BayesEstimate &);

				BayesEstimate & operator=(BayesEstimate &&) noexcept(true);

				void compute(VAf64 &,
					     int32_t,
					     int32_t,
					     VAi32 &,
					     VAi32 &,
					     int32_t[6]);


			};
		}
	}
}



#endif /*__GMS_BAYESIAN_ESTIMATE_H__*/
