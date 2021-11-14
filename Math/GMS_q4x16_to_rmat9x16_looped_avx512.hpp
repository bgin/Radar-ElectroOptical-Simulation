

#ifndef __GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_HPP__
#define __GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_HPP__ 141120211425


namespace file_info {

   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MAJOR = 1;
   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MINOR = 0;
   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MICRO = 0;
   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_FULLVER =
    1000*GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MAJOR+
    100*GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MINOR+
    10*GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MICRO;
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_CREATE_DATE = "14-11-2021 14:25 +00200 (SUN 14 NOV 2021 GMT+2)";
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_BUILD_DATE  = __DATE__":"__TIME__;
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_DESCRIPTION = "Optimized kernels (AVX512) for basic rotation operations";
}


#include <immintrin.h>
#include <cstdint>
#include "GMS_cephes.h"
#include "GMS_config.h"
#include "GMS_rotations_avx512_helpers.hpp"



namespace gms {

          namespace  math {


                  
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      void
                      q4x16_rm9x16_looped_u_nonunroll_zmm16r4(float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      const float * __restrict __ATTR_ALIGN__(64) q_x,
							      const float * __restrict __ATTR_ALIGN__(64) q_y,
							      const float * __restrict __ATTR_ALIGN__(64) q_z,
							      const float * __restrict __ATTR_ALIGN__(64) q_w,
							      const int32_t n) {

                               if(__builtin_expect(0<=n,0)) {return;}
			       int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif                         for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
                                   _mm_prefetch((const char*)&q_x[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+16],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_loadu_ps(&q_x[i]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_loadu_ps(&q_y[i]);
				   register const __m512 vqw = _mm512_loadu_ps(&q_w[i]);
				   _mm_prefetch((const char*)&q_z[i+16],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_loadu_ps(&q_z[i]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_storeu_ps(&row5[i],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row9[i],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_storeu_ps(&row2[i],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row7[i],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_storeu_ps(&row6[i],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_storeu_ps(&row3[i],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_storeu_ps(&row1[i],_mm512_fmadd_ps(v16_2,t4,t3));
				   
                               }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#pragma novector
#endif
                                    for(; i != n; ++i) {
                                        const double qx = q_x[i];
					const double t0 = qx*qx;
					const double qy = q_y[i];
					const double qz = q_z[i];
					const double qw = q_w[i];
					const double t1 = (qy+qz+qw)*(qy+qz+qw);
					const double q2 = t0-t1;
					row1[i] = q2+2.0F*qy*qy;
					row5[i] = q2+2.0F*qz*qz;
					row9[i] = q2+2.0F*qw*qw;
					
					row2[i] = 2.0F*(qy*qz-qx*qw);
					row6[i] = 2.0F*(qz*qw-qx*qy);
					row7[i] = 2.0F*(qw*qy-qx*qz);
					row4[i] = 2.0F*(qz*qy+qx*qw);
					row8[i] = 2.0F*(qw*qz+qx*qy);
					row3[i] = 2.0F*(qy*qw+qx*qz);
				    }
			       
		      }
							  



                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      void
                      q4x16_rm9x16_looped_a_nonunroll_zmm16r4(float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      float * __restrict __ATTR_ALIGN__(64) q_x,
							      float * __restrict __ATTR_ALIGN__(64) q_y,
							      float * __restrict __ATTR_ALIGN__(64) q_z,
							      float * __restrict __ATTR_ALIGN__(64) q_w,
							      const int32_t n) {

                               if(__builtin_expect(0<=n,0)) {return;}
			       int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
                               __assume_aligned(row1,64);
			       __assume_aligned(row2,64);
			       __assume_aligned(row3,64);
			       __assume_aligned(row4,64);
			       __assume_aligned(row5,64);
			       __assume_aligned(row6,64);
			       __assume_aligned(row7,64);
			       __assume_aligned(row8,64);
			       __assume_aligned(row9,64);
			       __assume_aligned(q_x,64);
			       __assume_aligned(q_y,64);
			       __assume_aligned(q_z,64);
			       __assume_aligned(q_w,64);
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC)
                               row1 = (float*)__builtin_assume_aligned(row1,64);
			       row2 = (float*)__builtin_assume_aligned(row2,64);
			       row3 = (float*)__builtin_assume_aligned(row3,64);
			       row4 = (float*)__builtin_assume_aligned(row4,64);
			       row5 = (float*)__builtin_assume_aligned(row5,64);
			       row6 = (float*)__builtin_assume_aligned(row6,64);
			       row7 = (float*)__builtin_assume_aligned(row7,64);
			       row8 = (float*)__builtin_assume_aligned(row8,64);
			       row9 = (float*)__builtin_assume_aligned(row9,64);
			       q_x = (float*)__builtin_assume_aligned(q_x,64);
			       q_y = (float*)__builtin_assume_aligned(q_y,64);
			       q_z = (float*)__builtin_assume_aligned(q_z,64);
			       q_w = (float*)__builtin_assume_aligned(q_w,64);
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif                         for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
                                   _mm_prefetch((const char*)&q_x[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+16],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_load_ps(&q_x[i]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_load_ps(&q_y[i]);
				   register const __m512 vqw = _mm512_load_ps(&q_w[i]);
				   _mm_prefetch((const char*)&q_z[i+16],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_load_ps(&q_z[i]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_store_ps(&row5[i],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row9[i],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_store_ps(&row2[i],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row4[i],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row7[i],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_store_ps(&row6[i],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row8[i],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_store_ps(&row3[i],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_store_ps(&row1[i],_mm512_fmadd_ps(v16_2,t4,t3));
				   
                               }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#pragma novector
#endif
                                    for(; i != n; ++i) {
                                        const double qx = q_x[i];
					const double t0 = qx*qx;
					const double qy = q_y[i];
					const double qz = q_z[i];
					const double qw = q_w[i];
					const double t1 = (qy+qz+qw)*(qy+qz+qw);
					const double q2 = t0-t1;
					row1[i] = q2+2.0F*qy*qy;
					row5[i] = q2+2.0F*qz*qz;
					row9[i] = q2+2.0F*qw*qw;
					
					row2[i] = 2.0F*(qy*qz-qx*qw);
					row6[i] = 2.0F*(qz*qw-qx*qy);
					row7[i] = 2.0F*(qw*qy-qx*qz);
					row4[i] = 2.0F*(qz*qy+qx*qw);
					row8[i] = 2.0F*(qw*qz+qx*qy);
					row3[i] = 2.0F*(qy*qw+qx*qz);
				    }
			       
		      }
							  









	  
     }

}










#endif /*__GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_HPP__*/
