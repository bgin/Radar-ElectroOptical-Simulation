
#ifndef __GMS_MARKOV_MODEL_H__
#define __GMS_MARKOV_MODEL_H__





namespace file_info {
#if defined _WIN64  
  #include "../GMS_version.h"
#elif defined __linux
  #include "GMS_version.h"
#endif

	const unsigned int gGMS_MARKOV_MODEL_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_MARKOV_MODEL_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_MARKOV_MODEL_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_MARKOV_MODEL_FULLVER = 
	 1000U*gGMS_MARKOV_MODEL_MAJOR+100U*gGMS_MARKOV_MODEL_MINOR+10U*gGMS_MARKOV_MODEL_MICRO;

	const char * const pgGMS_MARKOV_MODEL_CREATE_DATE = "15-10-2019 10:51 +00200 (TUE 15 OCT 2019 GMT+2)";

	const char * const pgGMS_MARKOV_MODEL_BUILD_DATE = "00-00-0000 00:00";

	const char * const pgGMS_MARKOV_MODEL_AUTHOR = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";

	const char * const pgGMS_MARKOV_MODEL_SYNOPSIS = "MAXIMUM LIKELIHOOD COMPUTATION OF MARKOVIAN MODEL (malloc aligned C-pointer version)"; 
}

#include <cstdint>

namespace gms {
	namespace math {
		namespace stat {
#if defined _WIN64					
			__declspec(align(64)) struct MarkovModel {
#elif defined __linux
			__attribute__((align(64))) struct MarkovModel {
#endif
				/*
					Public member access low-level computational
					code.
				*/

				int32_t	  m_k;

				int32_t   m_id;

				int32_t   m_mj4;

				int32_t   m_mj6;

				int32_t   m_mj7;

				int32_t   m_ipq;

				int32_t   m_iqm;

				double    m_aicd;
#if defined _WIN64

				_Field_size_(m_k)	      int32_t * __restrict m_idd;

				_Field_size_(m_k)	      int32_t * __restrict m_ir;

				_Field_size_(m_k)	      int32_t * __restrict m_ij;
				
				_Field_size_(m_k)	      int32_t * __restrict m_ik;

				_Field_size_(m_mj4)           double  * __restrict m_g;

				_Field_size_(m_k*m_k)         double  * __restrict m_a1;

				_Field_size_(m_k*m_k)         double  * __restrict m_a;

				_Field_size_(m_k*m_id)        double  * __restrict m_b;

				_Field_size_(m_mj4*m_mj4)     double  * __restrict m_vd;

				_Field_size_(m_id*m_id*m_mj6) double  * __restrict m_bm;

				_Field_size_(m_id*m_id*m_mj7) double  * __restrict m_au;

				_Field_size_(m_id*m_id*m_mj7) double  * __restrict m_zz;

				_Field_size_(m_id*m_id)       double  * __restrict m_c0;
#elif defined __linux
			     	 int32_t * __restrict m_idd;

				 int32_t * __restrict m_ir;

				 int32_t * __restrict m_ij;
				
				 int32_t * __restrict m_ik;

			         double  * __restrict m_g;

			         double  * __restrict m_a1;

			         double  * __restrict m_a;

			         double  * __restrict m_b;

			         double  * __restrict m_vd;

			         double  * __restrict m_bm;

			         double  * __restrict m_au;

			         double  * __restrict m_zz;

			         double  * __restrict m_c0;
#endif

				MarkovModel();

				MarkovModel(const int32_t,
					    const int32_t,
					    const int32_t,
					    const int32_t,
					    const int32_t);

				MarkovModel(const MarkovModel &);

				MarkovModel(MarkovModel &&) noexcept(true);

				~MarkovModel() noexcept(true);

				MarkovModel & operator=(const MarkovModel &);

				MarkovModel & operator=(MarkovModel &&) noexcept(true);

				void compute_markov_model( int32_t,
							   int32_t,
							   double * __restrict,
							   int32_t * __restrict,
							   int32_t,
							   double * __restrict,
							   double * __restrict,
							   int32_t,
							   int32_t ,
							   int32_t[6]);
			};
		}
	}
}


#endif /*__GMS_MARKOV_MODEL_H__*/
