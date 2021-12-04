
#ifndef __GMS_ROTATION_KERNELS_AVX512_OMP_H__
#define __GMS_ROTATION_KERNELS_AVX512_OMP_H__ 041220210910



namespace file_info {

const unsigned int gGMS_ROTATION_KERNELS_AVX512_OMP_MAJOR = 1U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_OMP_MINOR = 0U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_OMP_MICRO = 1U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_OMP_FULLVER =
       1000U*gGMS_ROTATION_KERNELS_AVX512_OMP_MAJOR+
       100U*gGMS_ROTATION_KERNELS_AVX512_OMP_MINOR +
       10U*gGMS_ROTATION_KERNELS_AVX512_OMP_MICRO;
const char * const pgGMS_ROTATION_KERNELS_AVX512_OMP_CREATION_DATE = "04-12-2021 09:10 AM +00200 (SAT 04 DEC 2021 GMT+2)";
const char * const pgGMS_ROTATION_KERNELS_AVX512_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_ROTATION_KERNELS_AVX512_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_ROTATION_KERNELS_AVX512_OMP_DESCRIPTION   = "AVX512 vectorized basic rotation operations OpenMP parallelized (coarse level).";
}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_dcm_avx512.hpp"

namespace  gms {

          namespace  math {


	          void
		  q4x16_to_rmat9x16_zmm16r4_omp(const __m512 * __restrict,
		                                const __m512 * __restrict,
						const __m512 * __restrict,
						const __m512 * __restrict,
						DCM9x16 * __restrict,
						const int32_t) __ATTR_HOT__
						               __ATTR_ALIGN__(32);


		  void
		  q4x8_to_rmat9x8_zmm8r8_omp(const __m512d * __restrict,
		                             const __m512d * __restrict,
					     const __m512d * __restrict,
					     const __m512d * __restrict,
					     DCM9x8 * __restrict,
					     const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


                  void
		  rand_sphere_rm9x16_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
						 const __m512 * __restrict,
						 DCM9x16 * __restrict,
						 const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		   void
		   rand_sphere_rm9x8_zmm8r8_omp(const __m512d * __restrict,
		                                const __m512d * __restrict,
					        const __m512d * __restrict,
					        DCM9x8 * __restrict,
					        const int32_t) __ATTR_HOT__
						               __ATTR_ALIGN__(32);


		   void
		   urand_q4x16_a_zmm16r4_omp(const __m512 * __restrict,
		                             const __m512 * __restrict,
					     const __m512 * __restrict,
					     const __m512 * __restrict,
					     float * __restrict,
					     float * __restrict,
					     float * __restrict,
					     float * __restrict,
					     const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


		   void
		   urand_q4x16_u_zmm16r4_omp(const __m512 * __restrict,
		                             const __m512 * __restrict,
					     const __m512 * __restrict,
					     const __m512 * __restrict,
					     float * __restrict,
					     float * __restrict,
					     float * __restrict,
					     float * __restrict,
					     const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


		   void
		   urand_q4x8_a_zmm8r8_omp(  const __m512d * __restrict,
		                             const __m512d * __restrict,
					     const __m512d * __restrict,
					     const __m512d * __restrict,
					     double * __restrict,
					     double * __restrict,
					     double * __restrict,
					     double * __restrict,
					     const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


		   void
		   urand_q4x8_u_zmm8r8_omp(  const __m512d * __restrict,
		                             const __m512d * __restrict,
					     const __m512d * __restrict,
					     const __m512d * __restrict,
					     double * __restrict,
					     double * __restrict,
					     double * __restrict,
					     double * __restrict,
					     const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


		   void
		   q4x16_to_ea3x16_a_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
					         const __m512 * __restrict,
					         const __m512 * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					     	 const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		   void
		   q4x16_to_ea3x16_u_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
					         const __m512 * __restrict,
					         const __m512 * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					     	 const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		   void
		   q4x8_to_ea3x8_a_zmm8r8_omp(const __m512d * __restrict,
		                              const __m512d * __restrict,
					      const __m512d * __restrict,
					      const __m512d * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      const int32_t) __ATTR_HOT__
						             __ATTR_ALIGN__(32);


		   void
		   q4x8_to_ea3x8_u_zmm8r8_omp(const __m512d * __restrict,
		                              const __m512d * __restrict,
					      const __m512d * __restrict,
					      const __m512d * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      const int32_t) __ATTR_HOT__
						             __ATTR_ALIGN__(32);


		   void
		   q4x16_to_ax4x16_a_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
					         const __m512 * __restrict,
					         const __m512 * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		   void
		   q4x16_to_ax4x16_u_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
					         const __m512 * __restrict,
					         const __m512 * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		   void
		   q4x8_to_ax4x8_a_zmm8r8_omp(const __m512d * __restrict,
		                              const __m512d * __restrict,
					      const __m512d * __restrict,
					      const __m512d * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


		   void
		   q4x8_to_ax4x8_u_zmm8r8_omp(const __m512d * __restrict,
		                              const __m512d * __restrict,
					      const __m512d * __restrict,
					      const __m512d * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      const int32_t) __ATTR_HOT__
						            __ATTR_ALIGN__(32);


		   void
		   q4x16_to_rv4x16_a_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
					         const __m512 * __restrict,
					         const __m512 * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		   void
		   q4x16_to_rv4x16_u_zmm16r4_omp(const __m512 * __restrict,
		                                 const __m512 * __restrict,
					         const __m512 * __restrict,
					         const __m512 * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         float * __restrict,
					         const int32_t) __ATTR_HOT__
						                __ATTR_ALIGN__(32);


		  void
		  q4x8_to_rv4x8_a_zmm8r8_omp(const __m512d * __restrict,
		                              const __m512d * __restrict,
					      const __m512d * __restrict,
					      const __m512d * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      const int32_t) __ATTR_HOT__
						             __ATTR_ALIGN__(32);


		  void
		  q4x8_to_rv4x8_u_zmm8r8_omp( const __m512d * __restrict,
		                              const __m512d * __restrict,
					      const __m512d * __restrict,
					      const __m512d * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      double * __restrict,
					      const int32_t) __ATTR_HOT__
						             __ATTR_ALIGN__(32);


		  void
		  rmat9x16_to_ea3x16_a_zmm16r4_omp(const DCM9x16 * __restrict,
		                                   float * __restrict,
					           float * __restrict,
					           float * __restrict,
					           const int32_t) __ATTR_HOT__
						                  __ATTR_ALIGN__(32);


		  void
		  rmat9x16_to_ea3x16_u_zmm16r4_omp(const DCM9x16 * __restrict,
		                                   float * __restrict,
					           float * __restrict,
					           float * __restrict,
					           const int32_t) __ATTR_HOT__
						                  __ATTR_ALIGN__(32);


		  void
		  rmat9x8_to_ea3x8_a_zmm8r8_omp(const DCM9x8 * __restrict,
		                                double * __restrict,
					        double * __restrict,
					        double * __restrict,
					        const int32_t) __ATTR_HOT__
						               __ATTR_ALIGN__(32);


		  void
		  rmat9x8_to_ea3x8_u_zmm8r8_omp(const DCM9x8 * __restrict,
		                                double * __restrict,
					        double * __restrict,
					        double * __restrict,
					        const int32_t) __ATTR_HOT__
						               __ATTR_ALIGN__(32);


		  


		  



		   


							     

     }


}






#endif /*__GMS_ROTATION_KERNELS_AVX512_OMP_H__*/
