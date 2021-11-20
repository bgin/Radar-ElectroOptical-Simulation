

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



                     
		      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      void
                      q4x16_rm9x16_looped_u_unroll4x_zmm16r4(float * __restrict __ATTR_ALIGN__(64) row1,
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
#endif
			      for(i = 0; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_x[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+64],_MM_HINT_T0);
				   register const __m512 vqx  = _mm512_loadu_ps(&q_x[i+0]);
				   register const __m512 vqx2 = _mm512_loadu_ps(&q_x[i+16]);
				   register const __m512 vqx3 = _mm512_loadu_ps(&q_x[i+32]);
				   register const __m512 vqx4 = _mm512_loadu_ps(&q_x[i+48]);
				   register const __m512 t0   = _mm512_mul_ps(vqx,vqx);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 t03  = _mm512_mul_ps(vqx3,vqx3);
				   register const __m512 t04  = _mm512_mul_ps(vqx4,vqx4);
				   register const __m512 vqy  = _mm512_loadu_ps(&q_y[i+0]);
				   register const __m512 vqy2 = _mm512_loadu_ps(&q_y[i+16]);
				   register const __m512 vqy3 = _mm512_loadu_ps(&q_y[i+32]);
				   register const __m512 vqy4 = _mm512_loadu_ps(&q_y[i+48]);
				   register const __m512 vqw  = _mm512_loadu_ps(&q_w[i+0]);
				   register const __m512 vqw2 = _mm512_loadu_ps(&q_w[i+16]);
				   register const __m512 vqw3 = _mm512_loadu_ps(&q_w[i+32]);
				   register const __m512 vqw4 = _mm512_loadu_ps(&q_w[i+48]);
				   _mm_prefetch((const char*)&q_z[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_z[i+64],_MM_HINT_T0);
				   register const __m512 vqz  = _mm512_loadu_ps(&q_z[i+0]);
				   register const __m512 vqz2 = _mm512_loadu_ps(&q_z[i+16]);
				   register const __m512 vqz3 = _mm512_loadu_ps(&q_z[i+32]);
				   register const __m512 vqz4 = _mm512_loadu_ps(&q_z[i+48]);
				   register const __m512 t1   = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t13  = _mm512_add_ps(_mm512_add_ps(vqx3,vqy3),vqw3);
				   register const __m512 t14  = _mm512_add_ps(_mm512_add_ps(vqx4,vqy4),vqw4);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t23  = _mm512_mul_ps(t13,t13);
				   register const __m512 t24  = _mm512_mul_ps(t14,t14);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   register const __m512 t33  = _mm512_sub_ps(t03,t23);
				   register const __m512 t34  = _mm512_sub_ps(t04,t24);
				   _mm512_storeu_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_storeu_ps(&row5[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz3,vqz3),t33));
				   _mm512_storeu_ps(&row5[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz4,vqz4),t34));
				   _mm512_storeu_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   _mm512_storeu_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   _mm512_storeu_ps(&row9[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vq3,vq3),t33));
				   _mm512_storeu_ps(&row9[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw4,vqw4),t34));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t52 = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t53 = _mm512_mul_ps(vqy3,vqz3);
				   register const __m512 t54 = _mm512_mul_ps(vqy4,vqz4);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t63  = _mm512_mul_ps(vqz3,vqw3);
				   register const __m512 t64  = _mm512_mul_ps(vqz4,vqw4);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t73  = _mm512_mul_ps(vqw3,vqy3);
				   register const __m512 t74  = _mm512_mul_ps(vqw4,vqy4);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   register const __m512 t83  = _mm512_mul_ps(vqx3,vqw3);
				   register const __m512 t84  = _mm512_mul_ps(vqx4,vqw4);
				   _mm512_storeu_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_storeu_ps(&row2[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t53,t83)));
				   _mm512_storeu_ps(&row2[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t54,t84)));
				   _mm512_storeu_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_storeu_ps(&row4[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t53,t83)));
				   _mm512_storeu_ps(&row4[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t54,t84)));
				   _mm512_storeu_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   _mm512_storeu_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   _mm512_storeu_ps(&row7[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t73,t83)));
				   _mm512_storeu_ps(&row7[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t74,t84)));
				   register const __m512 t9   = _mm512_mul_ps(vqx,vqy);
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   register const __m512 t93  = _mm512_mul_ps(vqx3,vqy3);
				   register const __m512 t94  = _mm512_mul_ps(vqx4,vqy4);
				   _mm512_storeu_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_storeu_ps(&row6[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t63,t93)));
				   _mm512_storeu_ps(&row6[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t64,t94)));
				   _mm512_storeu_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   _mm512_storeu_ps(&row8[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t63,t93)));
				   _mm512_storeu_ps(&row8[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t64,t94)));
				   register const __m512 t10  = _mm512_mul_ps(vqx,vqz);
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   register const __m512 t103 = _mm512_mul_ps(vqx3,vqz3);
				   register const __m512 t104 = _mm512_mul_ps(vqx4,vqz4);
				   _mm512_storeu_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   _mm512_storeu_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   _mm512_storeu_ps(&row3[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t73,t103)));
				   _mm512_storeu_ps(&row3[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t74,t104)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   register const __m512 t42 = _mm512_mul_ps(vqy2,vqy2);
				   register const __m512 t43 = _mm512_mul_ps(vqy3,vqy3);
				   register const __m512 t44 = _mm512_mul_ps(vqy4,vqy4);
				   _mm512_storeu_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   _mm512_storeu_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
				   _mm512_storeu_ps(&row1[i+32],_mm512_fmadd_ps(v16_2,t43,t33));
				   _mm512_storeu_ps(&row1[i+48],_mm512_fmadd_ps(v16_2,t44,t34));
#else
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_loadu_ps(&q_x[i+0]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_loadu_ps(&q_y[i+0]);
				   register const __m512 vqw = _mm512_loadu_ps(&q_w[i+0]);
				   _mm_prefetch((const char*)&q_z[i+32],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_loadu_ps(&q_z[i+0]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_storeu_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_storeu_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_storeu_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_storeu_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_storeu_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   _mm_prefetch((const char*)&q_x[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+64],_MM_HINT_T0);
				   register const __m512 vqx2 = _mm512_loadu_ps(&q_x[i+16]);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 vqy2 = _mm512_loadu_ps(&q_y[i+16]);
				   register const __m512 vqw2 = _mm512_loadu_ps(&q_w[i+16]);
				   _mm_prefetch((const char*)&q_z[i+64],_MM_HINT_T0);
				   register const __m512 vqz2 = _mm512_loadu_ps(&q_z[i+16]);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   _mm512_storeu_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_storeu_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   register const __m512 t52  = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   _mm512_storeu_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_storeu_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_storeu_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   _mm512_storeu_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_storeu_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   _mm512_storeu_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   register const __m512 t42  = _mm512_mul_ps(vqy2,vqy2);
				   _mm512_storeu_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
				   register const __m512 vqx3 = _mm512_loadu_ps(&q_x[i+32]);
				   register const __m512 t03  = _mm512_mul_ps(vqx3,vqx3);
				   register const __m512 vqy3 = _mm512_loadu_ps(&q_y[i+32]);
				   register const __m512 vqw3 = _mm512_loadu_ps(&q_w[i+32]);
				   register const __m512 vqz3 = _mm512_loadu_ps(&q_z[i+32]);
				   register const __m512 t13  = _mm512_add_ps(_mm512_add_ps(vqx3,vqy3),vqw3);
				   register const __m512 t23  = _mm512_mul_ps(t13,t13);
				   register const __m512 t33  = _mm512_sub_ps(t03,t23);
				   _mm512_storeu_ps(&row5[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz3,vqz3),t33));
				   _mm512_storeu_ps(&row9[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw3,vqw3),t33));
				   register const __m512 t53  = _mm512_mul_ps(vqy3,vqz3);
				   register const __m512 t63  = _mm512_mul_ps(vqz3,vqw3);
				   register const __m512 t73  = _mm512_mul_ps(vqw3,vqy3);
				   register const __m512 t83  = _mm512_mul_ps(vqx3,vqw3);
				   _mm512_storeu_ps(&row2[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t53,t83)));
				   _mm512_storeu_ps(&row4[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t53,t83)));
				   _mm512_storeu_ps(&row7[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t73,t83)));
				   register const __m512 t93  = _mm512_mul_ps(vqx3,vqy3);
				   _mm512_storeu_ps(&row6[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t63,t93)));
				   _mm512_storeu_ps(&row8[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t63,t93)));
				   register const __m512 t103 = _mm512_mul_ps(vqx3,vqz3);
				   _mm512_storeu_ps(&row3[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t73,t103)));
				   register const __m512 t43  = _mm512_mul_ps(vqy3,vqy3);
				   _mm512_storeu_ps(&row1[i+32],_mm512_fmadd_ps(v16_2,t43,t33));
				   _mm_prefetch((const char*)&q_x[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+64],_MM_HINT_T0);
				   register const __m512 vqx4 = _mm512_loadu_ps(&q_x[i+48]);
				   register const __m512 t04  = _mm512_mul_ps(vqx4,vqx4);
				   register const __m512 vqy4 = _mm512_loadu_ps(&q_y[i+48]);
				   register const __m512 vqw4 = _mm512_loadu_ps(&q_w[i+48]);
				   _mm_prefetch((const char*)&q_z[i+64],_MM_HINT_T0);
				   register const __m512 vqz4 = _mm512_loadu_ps(&q_z[i]);
				   register const __m512 t14  = _mm512_add_ps(_mm512_add_ps(vqx4,vqy4),vqw4);
				   register const __m512 t24  = _mm512_mul_ps(t14,t14);
				   register const __m512 t34  = _mm512_sub_ps(t04,t24);
				   _mm512_storeu_ps(&row5[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz4,vqz4),t34));
				   _mm512_storeu_ps(&row9[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw4,vqw4),t34));
				   register const __m512 t54  = _mm512_mul_ps(vqy4,vqz4);
				   register const __m512 t64  = _mm512_mul_ps(vqz4,vqw4);
				   register const __m512 t74  = _mm512_mul_ps(vqw4,vqy4);
				   register const __m512 t84  = _mm512_mul_ps(vqx4,vqw4);
				   _mm512_storeu_ps(&row2[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t54,t84)));
				   _mm512_storeu_ps(&row4[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t54,t84)));
				   _mm512_storeu_ps(&row7[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t74,t84)));
				   register const __m512 t94  = _mm512_mul_ps(vqx4,vqy4);
				   _mm512_storeu_ps(&row6[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t64,t94)));
				   _mm512_storeu_ps(&row8[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t64,t94)));
				   register const __m512 t104 = _mm512_mul_ps(vqx4,vqz4);
				   _mm512_storeu_ps(&row3[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t74,t104)));
				   register const __m512 t44  = _mm512_mul_ps(vqy4,vqy4);
				   _mm512_storeu_ps(&row1[i+48],_mm512_fmadd_ps(v16_2,t44,t34));
#endif
			       }

			       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   register const __m512 vqx  = _mm512_loadu_ps(&q_x[i+0]);
				   register const __m512 vqx2 = _mm512_loadu_ps(&q_x[i+16]);
				   register const __m512 t0   = _mm512_mul_ps(vqx,vqx);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 vqy  = _mm512_loadu_ps(&q_y[i+0]);
				   register const __m512 vqy2 = _mm512_loadu_ps(&q_y[i+16]);
				   register const __m512 vqw  = _mm512_loadu_ps(&q_w[i+0]);
				   register const __m512 vqw2 = _mm512_loadu_ps(&q_w[i+16]);
				   register const __m512 vqz  = _mm512_loadu_ps(&q_z[i+0]);
				   register const __m512 vqz2 = _mm512_loadu_ps(&q_z[i+16]);
				   register const __m512 t1   = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   _mm512_storeu_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_storeu_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   _mm512_storeu_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t52 = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   _mm512_storeu_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_storeu_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_storeu_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   _mm512_storeu_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   register const __m512 t9   = _mm512_mul_ps(vqx,vqy);
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   _mm512_storeu_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_storeu_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   register const __m512 t10  = _mm512_mul_ps(vqx,vqz);
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   _mm512_storeu_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   _mm512_storeu_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   register const __m512 t42 = _mm512_mul_ps(vqy2,vqy2);
				   _mm512_storeu_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   _mm512_storeu_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
				  
#else
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_loadu_ps(&q_x[i+0]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_loadu_ps(&q_y[i+0]);
				   register const __m512 vqw = _mm512_loadu_ps(&q_w[i+0]);
				   _mm_prefetch((const char*)&q_z[i+32],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_loadu_ps(&q_z[i+0]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_storeu_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_storeu_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_storeu_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_storeu_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_storeu_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   register const __m512 vqx2 = _mm512_loadu_ps(&q_x[i+16]);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 vqy2 = _mm512_loadu_ps(&q_y[i+16]);
				   register const __m512 vqw2 = _mm512_loadu_ps(&q_w[i+16]);
				   register const __m512 vqz2 = _mm512_loadu_ps(&q_z[i+16]);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   _mm512_storeu_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_storeu_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   register const __m512 t52  = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   _mm512_storeu_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_storeu_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_storeu_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   _mm512_storeu_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_storeu_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   _mm512_storeu_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   register const __m512 t42  = _mm512_mul_ps(vqy2,vqy2);
				   _mm512_storeu_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
		
#endif
				 
                                }

			      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   _mm_prefetch((const char*)&q_x[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+16],_MM_HINT_T0);
				   register const __m512 vqx  = _mm512_loadu_ps(&q_x[i+0]);
				   register const __m512 t0   = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy  = _mm512_loadu_ps(&q_y[i+0]);
				   register const __m512 vqw  = _mm512_loadu_ps(&q_w[i+0]);
				   register const __m512 vqz  = _mm512_loadu_ps(&q_z[i+0]);
				   register const __m512 t1   = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_storeu_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_storeu_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9   = _mm512_mul_ps(vqx,vqy);
				   _mm512_storeu_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10  = _mm512_mul_ps(vqx,vqz);
				   _mm512_storeu_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_storeu_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				 
#else
                                   _mm_prefetch((const char*)&q_x[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+16],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_loadu_ps(&q_x[i+0]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_loadu_ps(&q_y[i+0]);
				   register const __m512 vqw = _mm512_loadu_ps(&q_w[i+0]);
				   _mm_prefetch((const char*)&q_z[i+16],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_loadu_ps(&q_z[i+0]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_storeu_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_storeu_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_storeu_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_storeu_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_storeu_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_storeu_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_storeu_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_storeu_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
#endif

				}

                                for(; (i+0) < n; i += 1) {
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
                      q4x16_rm9x16_looped_a_unroll4x_zmm16r4( float * __restrict __ATTR_ALIGN__(64) row1,
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
#endif
			      for(i = 0; (i+63) < n; i += 64) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_x[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+64],_MM_HINT_T0);
				   register const __m512 vqx  = _mm512_load_ps(&q_x[i+0]);
				   register const __m512 vqx2 = _mm512_load_ps(&q_x[i+16]);
				   register const __m512 vqx3 = _mm512_load_ps(&q_x[i+32]);
				   register const __m512 vqx4 = _mm512_load_ps(&q_x[i+48]);
				   register const __m512 t0   = _mm512_mul_ps(vqx,vqx);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 t03  = _mm512_mul_ps(vqx3,vqx3);
				   register const __m512 t04  = _mm512_mul_ps(vqx4,vqx4);
				   register const __m512 vqy  = _mm512_load_ps(&q_y[i+0]);
				   register const __m512 vqy2 = _mm512_load_ps(&q_y[i+16]);
				   register const __m512 vqy3 = _mm512_load_ps(&q_y[i+32]);
				   register const __m512 vqy4 = _mm512_load_ps(&q_y[i+48]);
				   register const __m512 vqw  = _mm512_load_ps(&q_w[i+0]);
				   register const __m512 vqw2 = _mm512_load_ps(&q_w[i+16]);
				   register const __m512 vqw3 = _mm512_load_ps(&q_w[i+32]);
				   register const __m512 vqw4 = _mm512_load_ps(&q_w[i+48]);
				   _mm_prefetch((const char*)&q_z[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_z[i+64],_MM_HINT_T0);
				   register const __m512 vqz  = _mm512_load_ps(&q_z[i+0]);
				   register const __m512 vqz2 = _mm512_load_ps(&q_z[i+16]);
				   register const __m512 vqz3 = _mm512_load_ps(&q_z[i+32]);
				   register const __m512 vqz4 = _mm512_load_ps(&q_z[i+48]);
				   register const __m512 t1   = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t13  = _mm512_add_ps(_mm512_add_ps(vqx3,vqy3),vqw3);
				   register const __m512 t14  = _mm512_add_ps(_mm512_add_ps(vqx4,vqy4),vqw4);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t23  = _mm512_mul_ps(t13,t13);
				   register const __m512 t24  = _mm512_mul_ps(t14,t14);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   register const __m512 t33  = _mm512_sub_ps(t03,t23);
				   register const __m512 t34  = _mm512_sub_ps(t04,t24);
				   _mm512_store_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_store_ps(&row5[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz3,vqz3),t33));
				   _mm512_store_ps(&row5[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz4,vqz4),t34));
				   _mm512_store_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   _mm512_store_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   _mm512_store_ps(&row9[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vq3,vq3),t33));
				   _mm512_store_ps(&row9[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw4,vqw4),t34));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t52 = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t53 = _mm512_mul_ps(vqy3,vqz3);
				   register const __m512 t54 = _mm512_mul_ps(vqy4,vqz4);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t63  = _mm512_mul_ps(vqz3,vqw3);
				   register const __m512 t64  = _mm512_mul_ps(vqz4,vqw4);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t73  = _mm512_mul_ps(vqw3,vqy3);
				   register const __m512 t74  = _mm512_mul_ps(vqw4,vqy4);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   register const __m512 t83  = _mm512_mul_ps(vqx3,vqw3);
				   register const __m512 t84  = _mm512_mul_ps(vqx4,vqw4);
				   _mm512_store_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_store_ps(&row2[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t53,t83)));
				   _mm512_store_ps(&row2[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t54,t84)));
				   _mm512_store_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_store_ps(&row4[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t53,t83)));
				   _mm512_store_ps(&row4[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t54,t84)));
				   _mm512_store_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   _mm512_store_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   _mm512_store_ps(&row7[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t73,t83)));
				   _mm512_store_ps(&row7[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t74,t84)));
				   register const __m512 t9   = _mm512_mul_ps(vqx,vqy);
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   register const __m512 t93  = _mm512_mul_ps(vqx3,vqy3);
				   register const __m512 t94  = _mm512_mul_ps(vqx4,vqy4);
				   _mm512_store_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_store_ps(&row6[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t63,t93)));
				   _mm512_store_ps(&row6[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t64,t94)));
				   _mm512_store_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   _mm512_store_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   _mm512_store_ps(&row8[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t63,t93)));
				   _mm512_store_ps(&row8[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t64,t94)));
				   register const __m512 t10  = _mm512_mul_ps(vqx,vqz);
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   register const __m512 t103 = _mm512_mul_ps(vqx3,vqz3);
				   register const __m512 t104 = _mm512_mul_ps(vqx4,vqz4);
				   _mm512_store_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   _mm512_store_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   _mm512_store_ps(&row3[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t73,t103)));
				   _mm512_store_ps(&row3[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t74,t104)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   register const __m512 t42 = _mm512_mul_ps(vqy2,vqy2);
				   register const __m512 t43 = _mm512_mul_ps(vqy3,vqy3);
				   register const __m512 t44 = _mm512_mul_ps(vqy4,vqy4);
				   _mm512_storeu_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   _mm512_storeu_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
				   _mm512_storeu_ps(&row1[i+32],_mm512_fmadd_ps(v16_2,t43,t33));
				   _mm512_storeu_ps(&row1[i+48],_mm512_fmadd_ps(v16_2,t44,t34));
#else
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_load_ps(&q_x[i+0]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_load_ps(&q_y[i+0]);
				   register const __m512 vqw = _mm512_load_ps(&q_w[i+0]);
				   _mm_prefetch((const char*)&q_z[i+32],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_load_ps(&q_z[i+0]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_store_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_store_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_store_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_store_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_store_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   _mm_prefetch((const char*)&q_x[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+64],_MM_HINT_T0);
				   register const __m512 vqx2 = _mm512_load_ps(&q_x[i+16]);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 vqy2 = _mm512_load_ps(&q_y[i+16]);
				   register const __m512 vqw2 = _mm512_load_ps(&q_w[i+16]);
				   _mm_prefetch((const char*)&q_z[i+64],_MM_HINT_T0);
				   register const __m512 vqz2 = _mm512_load_ps(&q_z[i+16]);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   _mm512_store_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_store_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   register const __m512 t52  = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   _mm512_store_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_store_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_store_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   _mm512_store_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_store_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   _mm512_store_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   register const __m512 t42  = _mm512_mul_ps(vqy2,vqy2);
				   _mm512_store_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
				   register const __m512 vqx3 = _mm512_load_ps(&q_x[i+32]);
				   register const __m512 t03  = _mm512_mul_ps(vqx3,vqx3);
				   register const __m512 vqy3 = _mm512_load_ps(&q_y[i+32]);
				   register const __m512 vqw3 = _mm512_load_ps(&q_w[i+32]);
				   register const __m512 vqz3 = _mm512_load_ps(&q_z[i+32]);
				   register const __m512 t13  = _mm512_add_ps(_mm512_add_ps(vqx3,vqy3),vqw3);
				   register const __m512 t23  = _mm512_mul_ps(t13,t13);
				   register const __m512 t33  = _mm512_sub_ps(t03,t23);
				   _mm512_store_ps(&row5[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz3,vqz3),t33));
				   _mm512_store_ps(&row9[i+32],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw3,vqw3),t33));
				   register const __m512 t53  = _mm512_mul_ps(vqy3,vqz3);
				   register const __m512 t63  = _mm512_mul_ps(vqz3,vqw3);
				   register const __m512 t73  = _mm512_mul_ps(vqw3,vqy3);
				   register const __m512 t83  = _mm512_mul_ps(vqx3,vqw3);
				   _mm512_store_ps(&row2[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t53,t83)));
				   _mm512_store_ps(&row4[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t53,t83)));
				   _mm512_store_ps(&row7[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t73,t83)));
				   register const __m512 t93  = _mm512_mul_ps(vqx3,vqy3);
				   _mm512_store_ps(&row6[i+32],_mm512_mul_ps(v16_2,_mm512_sub_ps(t63,t93)));
				   _mm512_store_ps(&row8[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t63,t93)));
				   register const __m512 t103 = _mm512_mul_ps(vqx3,vqz3);
				   _mm512_store_ps(&row3[i+32],_mm512_mul_ps(v16_2,_mm512_add_ps(t73,t103)));
				   register const __m512 t43  = _mm512_mul_ps(vqy3,vqy3);
				   _mm512_store_ps(&row1[i+32],_mm512_fmadd_ps(v16_2,t43,t33));
				   _mm_prefetch((const char*)&q_x[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+64],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+64],_MM_HINT_T0);
				   register const __m512 vqx4 = _mm512_load_ps(&q_x[i+48]);
				   register const __m512 t04  = _mm512_mul_ps(vqx4,vqx4);
				   register const __m512 vqy4 = _mm512_load_ps(&q_y[i+48]);
				   register const __m512 vqw4 = _mm512_load_ps(&q_w[i+48]);
				   _mm_prefetch((const char*)&q_z[i+64],_MM_HINT_T0);
				   register const __m512 vqz4 = _mm512_load_ps(&q_z[i]);
				   register const __m512 t14  = _mm512_add_ps(_mm512_add_ps(vqx4,vqy4),vqw4);
				   register const __m512 t24  = _mm512_mul_ps(t14,t14);
				   register const __m512 t34  = _mm512_sub_ps(t04,t24);
				   _mm512_storeu_ps(&row5[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz4,vqz4),t34));
				   _mm512_storeu_ps(&row9[i+48],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw4,vqw4),t34));
				   register const __m512 t54  = _mm512_mul_ps(vqy4,vqz4);
				   register const __m512 t64  = _mm512_mul_ps(vqz4,vqw4);
				   register const __m512 t74  = _mm512_mul_ps(vqw4,vqy4);
				   register const __m512 t84  = _mm512_mul_ps(vqx4,vqw4);
				   _mm512_store_ps(&row2[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t54,t84)));
				   _mm512_store_ps(&row4[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t54,t84)));
				   _mm512_store_ps(&row7[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t74,t84)));
				   register const __m512 t94  = _mm512_mul_ps(vqx4,vqy4);
				   _mm512_store_ps(&row6[i+48],_mm512_mul_ps(v16_2,_mm512_sub_ps(t64,t94)));
				   _mm512_store_ps(&row8[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t64,t94)));
				   register const __m512 t104 = _mm512_mul_ps(vqx4,vqz4);
				   _mm512_store_ps(&row3[i+48],_mm512_mul_ps(v16_2,_mm512_add_ps(t74,t104)));
				   register const __m512 t44  = _mm512_mul_ps(vqy4,vqy4);
				   _mm512_store_ps(&row1[i+48],_mm512_fmadd_ps(v16_2,t44,t34));
#endif
			       }

			       for(; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   register const __m512 vqx  = _mm512_load_ps(&q_x[i+0]);
				   register const __m512 vqx2 = _mm512_load_ps(&q_x[i+16]);
				   register const __m512 t0   = _mm512_mul_ps(vqx,vqx);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 vqy  = _mm512_load_ps(&q_y[i+0]);
				   register const __m512 vqy2 = _mm512_load_ps(&q_y[i+16]);
				   register const __m512 vqw  = _mm512_load_ps(&q_w[i+0]);
				   register const __m512 vqw2 = _mm512_load_ps(&q_w[i+16]);
				   register const __m512 vqz  = _mm512_load_ps(&q_z[i+0]);
				   register const __m512 vqz2 = _mm512_load_ps(&q_z[i+16]);
				   register const __m512 t1   = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   _mm512_store_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_store_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   _mm512_store_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t52 = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   _mm512_store_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_store_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_store_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   _mm512_store_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   register const __m512 t9   = _mm512_mul_ps(vqx,vqy);
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   _mm512_store_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_store_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   _mm512_store_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   register const __m512 t10  = _mm512_mul_ps(vqx,vqz);
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   _mm512_store_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   _mm512_store_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   register const __m512 t42 = _mm512_mul_ps(vqy2,vqy2);
				   _mm512_store_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   _mm512_store_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
				  
#else
                                   _mm_prefetch((const char*)&q_x[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+32],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+32],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_load_ps(&q_x[i+0]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_load_ps(&q_y[i+0]);
				   register const __m512 vqw = _mm512_load_ps(&q_w[i+0]);
				   _mm_prefetch((const char*)&q_z[i+32],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_load_ps(&q_z[i+0]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_store_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_store_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_store_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_store_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_store_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				   register const __m512 vqx2 = _mm512_load_ps(&q_x[i+16]);
				   register const __m512 t02  = _mm512_mul_ps(vqx2,vqx2);
				   register const __m512 vqy2 = _mm512_load_ps(&q_y[i+16]);
				   register const __m512 vqw2 = _mm512_load_ps(&q_w[i+16]);
				   register const __m512 vqz2 = _mm512_load_ps(&q_z[i+16]);
				   register const __m512 t12  = _mm512_add_ps(_mm512_add_ps(vqx2,vqy2),vqw2);
				   register const __m512 t22  = _mm512_mul_ps(t12,t12);
				   register const __m512 t32  = _mm512_sub_ps(t02,t22);
				   _mm512_storeu_ps(&row5[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz2,vqz2),t32));
				   _mm512_storeu_ps(&row9[i+16],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw2,vqw2),t32));
				   register const __m512 t52  = _mm512_mul_ps(vqy2,vqz2);
				   register const __m512 t62  = _mm512_mul_ps(vqz2,vqw2);
				   register const __m512 t72  = _mm512_mul_ps(vqw2,vqy2);
				   register const __m512 t82  = _mm512_mul_ps(vqx2,vqw2);
				   _mm512_store_ps(&row2[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t52,t82)));
				   _mm512_store_ps(&row4[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t52,t82)));
				   _mm512_store_ps(&row7[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t72,t82)));
				   register const __m512 t92  = _mm512_mul_ps(vqx2,vqy2);
				   _mm512_store_ps(&row6[i+16],_mm512_mul_ps(v16_2,_mm512_sub_ps(t62,t92)));
				   _mm512_store_ps(&row8[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t62,t92)));
				   register const __m512 t102 = _mm512_mul_ps(vqx2,vqz2);
				   _mm512_store_ps(&row3[i+16],_mm512_mul_ps(v16_2,_mm512_add_ps(t72,t102)));
				   register const __m512 t42  = _mm512_mul_ps(vqy2,vqy2);
				   _mm512_store_ps(&row1[i+16],_mm512_fmadd_ps(v16_2,t42,t32));
		
#endif
				 
                                }

			      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                   _mm_prefetch((const char*)&q_x[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+16],_MM_HINT_T0);
				   register const __m512 vqx  = _mm512_load_ps(&q_x[i+0]);
				   register const __m512 t0   = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy  = _mm512_load_ps(&q_y[i+0]);
				   register const __m512 vqw  = _mm512_load_ps(&q_w[i+0]);
				   register const __m512 vqz  = _mm512_load_ps(&q_z[i+0]);
				   register const __m512 t1   = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_store_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_store_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9   = _mm512_mul_ps(vqx,vqy);
				   _mm512_store_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10  = _mm512_mul_ps(vqx,vqz);
				   _mm512_store_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				   _mm512_store_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
				 
#else
                                   _mm_prefetch((const char*)&q_x[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_y[i+16],_MM_HINT_T0);
				   _mm_prefetch((const char*)&q_w[i+16],_MM_HINT_T0);
				   register const __m512 vqx = _mm512_load_ps(&q_x[i+0]);
				   register const __m512 t0  = _mm512_mul_ps(vqx,vqx);
				   register const __m512 vqy = _mm512_load_ps(&q_y[i+0]);
				   register const __m512 vqw = _mm512_load_ps(&q_w[i+0]);
				   _mm_prefetch((const char*)&q_z[i+16],_MM_HINT_T0);
				   register const __m512 vqz = _mm512_load_ps(&q_z[i+0]);
				   register const __m512 t1  = _mm512_add_ps(_mm512_add_ps(vqx,vqy),vqw);
				   register const __m512 t2  = _mm512_mul_ps(t1,t1);
				   register const __m512 t3  = _mm512_sub_ps(t0,t2);
				   _mm512_store_ps(&row5[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqz,vqz),t3));
				   _mm512_store_ps(&row9[i+0],_mm512_fmadd_ps(v16_2,_mm512_mul_ps(vqw,vqw),t3));
				   register const __m512 t5  = _mm512_mul_ps(vqy,vqz);
				   register const __m512 t6  = _mm512_mul_ps(vqz,vqw);
				   register const __m512 t7  = _mm512_mul_ps(vqw,vqy);
				   register const __m512 t8  = _mm512_mul_ps(vqx,vqw);
				   _mm512_store_ps(&row2[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t5,t8)));
				   _mm512_store_ps(&row4[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t5,t8)));
				   _mm512_store_ps(&row7[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t7,t8)));
				   register const __m512 t9  = _mm512_mul_ps(vqx,vqy);
				   _mm512_store_ps(&row6[i+0],_mm512_mul_ps(v16_2,_mm512_sub_ps(t6,t9)));
				   _mm512_store_ps(&row8[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t6,t9)));
				   register const __m512 t10 = _mm512_mul_ps(vqx,vqz);
				   _mm512_store_ps(&row3[i+0],_mm512_mul_ps(v16_2,_mm512_add_ps(t7,t10)));
				   register const __m512 t4  = _mm512_mul_ps(vqy,vqy);
				  _mm512_store_ps(&row1[i+0],_mm512_fmadd_ps(v16_2,t4,t3));
#endif

				}

                                for(; (i+0) < n; i += 1) {
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
