
#ifndef __GMS_RAND_RMAT9X16_LOOPED_AVX512_HPP__
#define __GMS_RAND_RMAT9X16_LOOPED_AVX512_HPP__ 271120211007




namespace file_info {

   const unsigned int GMS_RAND_RMAT9X16_LOOPED_AVX512_MAJOR = 1;
   const unsigned int GMS_RAND_RMAT9X16_LOOPED_AVX512_MINOR = 0;
   const unsigned int GMS_RAND_RMAT9X16_LOOPED_AVX512_MICRO = 0;
   const unsigned int GMS_RAND_RMAT9X16_LOOPED_AVX512_FULLVER =
    1000*GMS_RAND_RMAT9X16_LOOPED_AVX512_MAJOR+
    100*GMS_RAND_RMAT9X16_LOOPED_AVX512_MINOR+
    10*GMS_RAND_RMAT9X16_LOOPED_AVX512_MICRO;
   const char * const GMS_RAND_RMAT9X16_LOOPED_AVX512_CREATE_DATE = "27-11-2021 10:07 +00200 (SAT 27 NOV 2021 GMT+2)";
   const char * const GMS_RAND_RMAT9X16_LOOPED_AVX512_BUILD_DATE  = __DATE__":"__TIME__;
   const char * const GMS_RAND_RMAT9X16_LOOPED_AVX512_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const GMS_RAND_RMAT9X16_LOOPED_AVX512_DESCRIPTION = "Random Rotation Matrix 9x16 (of sphere) kernel AVX512.";
}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_cephes.h"
#include "GMS_rotations_avx512_helpers.hpp"




namespace  gms {


           namespace   math {


	   
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      void
                      rand_rm9x16_looped_u_nonunroll_zmm16r4( float * __restrict row1,
		                                              float * __restrict row2,
							      float * __restrict row3,
							      float * __restrict row4,
							      float * __restrict row5,
							      float * __restrict  row6,
							      float * __restrict  row7,
							      float * __restrict  row8,
							      float * __restrict  row9,
							      const float * __restrict  rx, // ramdomly normally distributed vector [0,1]
							      const float * __restrict  ry, // ramdomly normally distributed vector [0,1]
							      const float * __restrict  rz, // ramdomly normally distributed vector [0,1]
							      const int32_t n) {

                               if(__builtin_expect(0<=n,0)) {return;}
			       int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif                         for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
                                     _mm_prefetch((const char*)&rx[i+16],_MM_HINT_T0);
				     register const __m512 vrx   = _mm512_loadu_ps(&rx[i]);
				     register const __m512 theta = _mm512_mul_ps(vrx,v16_2pi);
				     _mm_prefetch((const char*)&ry[i+16],_MM_HINT_T0);
				     register const __m512 vry   = _mm512_loadu_ps(&ry[i]);
				     register const __m512 phi   = _mm512_mul_ps(vry,v16_2pi);
				     _mm_prefetch((const char*)&rz[i+16],_MM_HINT_T0);
				     register const __m512 vrz   = _mm512_loadu_ps(&rz[i]);
				     register const __m512 z     = _mm512_mul_ps(vrz,vrz);
				     _mm512_storeu_ps(&row9[i], _mm512_sub_ps(v16_1,z));
				     const register __m512 r     = _mm512_sqrt_ps(z);
				     const register __m512 vx    = _mm512_mul_ps(r,_mm512_sin_ps(phi));
				     const register __m512 vy    = _mm512_mul_ps(r,_mm512_cos_ps(phi));
				     const register __m512 vz    = _mm512_sqrt_ps(_mm512_sub_ps(v16_2,z));
				     _mm512_storeu_ps(row6[i], _mm512_mul_ps(vy,vz));
				     _mm512_storeu_ps(row3[i], _mm512_mul_ps(vx,vz));
				     const register __m512 st    = _mm512_sin_ps(theta);
			             const register __m512 ct    = _mm512_cos_ps(theta);
			             const register __m512 sx    = _mm512_fmsub_ps(vx,ct,_mm512_mul_ps(vy,st));
				     _mm512_storeu_ps(row7[i], _mm512_mul_ps(vz,sx));
				     _mm512_storeu_ps(row1[i], _mm512_fmsub_ps(vx,sx,ct));
				     const register __m512 sy    = _mm512_fmadd_ps(vx,st,_mm512_mul_ps(vy,ct));
				     _mm512_storeu_ps(&row8[i], _mm512_mul_ps(vz,sy));
				     _mm512_storeu_ps(&row2[i], _mm512_fmsub_ps(vx,sy,st));
				     _mm512_storeu_ps(&row4[i], _mm512_fmadd_ps(vy,sx,st));
				     _mm512_storeu_ps(&row5[i], _mm512_fmsub_ps(vy,sy,ct));
                               }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#pragma novector
#endif
                                    for(; i != n; ++i) {
                                        const float x = rx[i];
					const float theta = x*6.2831853071795864769253F;
					const float st = cephes_sinf(theta);
					const float ct = cephes_cosf(theta);
					const float y = ry[i];
					const float phi = y*6.2831853071795864769253;
					const float z = rz[i];
					const float zz = z+z;
					row9[i]        = 1.0F-z;
					const float r  = cephes_sqrtf(zz);
					const float vx = cephes_sinf(phi) * r;
					const float vy = cephes_cosf(phi) * r;
					const float vz = cephes_sqrtf(2.0F-zz);
				        row3[i]        = vx * vz;
					row6[i]        = vy * vz;
					const float st = cephes_sinf(theta);
					const float ct = cephes_cosf(theta);
					const float sx = vx * ct - vy * st;
					row1[i]        = vx * sx - ct;
					const float sy = vx * st + vy * ct;
					row8[i]        = vz * sy;
					row2[i]        = vx * sy - st;
					row4[i]        = vy * sx + st;
					row5[i]        = vy * sy - ct;
				    }
		      }




		      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      void
                      rand_rm9x16_looped_a_nonunroll_zmm16r4( float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      const float * __restrict __ATTR_ALIGN__(64) rx, // ramdomly normally distributed vector [0,1]
							      const float * __restrict __ATTR_ALIGN__(64) ry, // ramdomly normally distributed vector [0,1]
							      const float * __restrict __ATTR_ALIGN__(64) rz, // ramdomly normally distributed vector [0,1]
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
			       __assume_aligned(rx,64);
			       __assume_aligned(ry,64);
			       __assume_aligned(rz,64);
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
			       rx = (float*)__builtin_assume_aligned(rx,64);
			       ry = (float*)__builtin_assume_aligned(ry,64);
			       rz = (float*)__builtin_assume_aligned(rz,64);
			       
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif                         for(i = 0; i != ROUND_TO_SIXTEEN(n,16); i += 16) {
                                     _mm_prefetch((const char*)&rx[i+16],_MM_HINT_T0);
				     register const __m512 vrx   = _mm512_loadu_ps(&rx[i]);
				     register const __m512 theta = _mm512_mul_ps(vrx,v16_2pi);
				     _mm_prefetch((const char*)&ry[i+16],_MM_HINT_T0);
				     register const __m512 vry   = _mm512_loadu_ps(&ry[i]);
				     register const __m512 phi   = _mm512_mul_ps(vry,v16_2pi);
				     _mm_prefetch((const char*)&rz[i+16],_MM_HINT_T0);
				     register const __m512 vrz   = _mm512_loadu_ps(&rz[i]);
				     register const __m512 z     = _mm512_mul_ps(vrz,vrz);
				     _mm512_storeu_ps(&row9[i], _mm512_sub_ps(v16_1,z));
				     const register __m512 r     = _mm512_sqrt_ps(z);
				     const register __m512 vx    = _mm512_mul_ps(r,_mm512_sin_ps(phi));
				     const register __m512 vy    = _mm512_mul_ps(r,_mm512_cos_ps(phi));
				     const register __m512 vz    = _mm512_sqrt_ps(_mm512_sub_ps(v16_2,z));
				     _mm512_storeu_ps(row6[i], _mm512_mul_ps(vy,vz));
				     _mm512_storeu_ps(row3[i], _mm512_mul_ps(vx,vz));
				     const register __m512 st    = _mm512_sin_ps(theta);
			             const register __m512 ct    = _mm512_cos_ps(theta);
			             const register __m512 sx    = _mm512_fmsub_ps(vx,ct,_mm512_mul_ps(vy,st));
				     _mm512_storeu_ps(row7[i], _mm512_mul_ps(vz,sx));
				     _mm512_storeu_ps(row1[i], _mm512_fmsub_ps(vx,sx,ct));
				     const register __m512 sy    = _mm512_fmadd_ps(vx,st,_mm512_mul_ps(vy,ct));
				     _mm512_storeu_ps(&row8[i], _mm512_mul_ps(vz,sy));
				     _mm512_storeu_ps(&row2[i], _mm512_fmsub_ps(vx,sy,st));
				     _mm512_storeu_ps(&row4[i], _mm512_fmadd_ps(vy,sx,st));
				     _mm512_storeu_ps(&row5[i], _mm512_fmsub_ps(vy,sy,ct));
                               }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(16)
#pragma novector
#endif
                                    for(; i != n; ++i) {
                                        const float x = rx[i];
					const float theta = x*6.2831853071795864769253F;
					const float st = cephes_sinf(theta);
					const float ct = cephes_cosf(theta);
					const float y = ry[i];
					const float phi = y*6.2831853071795864769253;
					const float z = rz[i];
					const float zz = z+z;
					row9[i]        = 1.0F-z;
					const float r  = cephes_sqrtf(zz);
					const float vx = cephes_sinf(phi) * r;
					const float vy = cephes_cosf(phi) * r;
					const float vz = cephes_sqrtf(2.0F-zz);
				        row3[i]        = vx * vz;
					row6[i]        = vy * vz;
					const float st = cephes_sinf(theta);
					const float ct = cephes_cosf(theta);
					const float sx = vx * ct - vy * st;
					row1[i]        = vx * sx - ct;
					const float sy = vx * st + vy * ct;
					row8[i]        = vz * sy;
					row2[i]        = vx * sy - st;
					row4[i]        = vy * sx + st;
					row5[i]        = vy * sy - ct;
				    }
		      }




		      
     }


}











#endif /*__GMS_RAND_RMAT9X16_LOOPED_AVX512_HPP__*/
