
#ifndef __GMS_RAND_RMAT9X8_LOOPED_AVX512_HPP__
#define __GMS_RAND_RMAT9X8_LOOPED_AVX512_HPP__ 271120211007




namespace file_info {

   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_MAJOR = 1;
   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_MINOR = 0;
   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_MICRO = 0;
   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_FULLVER =
    1000*GMS_RAND_RMAT9X8_LOOPED_AVX512_MAJOR+
    100*GMS_RAND_RMAT9X8_LOOPED_AVX512_MINOR+
    10*GMS_RAND_RMAT9X8_LOOPED_AVX512_MICRO;
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_CREATE_DATE = "27-11-2021 10:07 +00200 (SAT 27 NOV 2021 GMT+2)";
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_BUILD_DATE  = __DATE__":"__TIME__;
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_DESCRIPTION = "Random Rotation Matrix 9x8 (of sphere) kernel AVX512.";
}


#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include "GMS_cephes.h"
#include "GMS_rotations_avx512_helpers.hpp"




namespace  gms {


           namespace   math {


	   
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      void
                      rand_rm9x8_looped_u_nonunroll_zmm8r8(   double * __restrict row1,
		                                              double * __restrict row2,
							      double * __restrict row3,
							      double * __restrict row4,
							      double * __restrict row5,
							      double * __restrict  row6,
							      double * __restrict  row7,
							      double * __restrict  row8,
							      double * __restrict  row9,
							      const double * __restrict  rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  rz, // ramdomly normally distributed vector [0,1]
							      const int32_t n) {

                               if(__builtin_expect(0<=n,0)) {return;}
			       int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif                         for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                     _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0);
				     register const __m512d vrx   = _mm512_loadu_pd(&rx[i]);
				     register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				     _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0);
				     register const __m512d vry   = _mm512_loadu_pd(&ry[i]);
				     register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				     _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0);
				     register const __m512d vrz   = _mm512_loadu_pd(&rz[i]);
				     register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				     _mm512_storeu_pd(&row9[i], _mm512_sub_ps(v8_1,z));
				     const register __m512d r     = _mm512_sqrt_pd(z);
				     const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				     const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				     const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
				     _mm512_storeu_pd(row6[i], _mm512_mul_pd(vy,vz));
				     _mm512_storeu_pd(row3[i], _mm512_mul_pd(vx,vz));
				     const register __m512d st    = _mm512_sin_pd(theta);
			             const register __m512d ct    = _mm512_cos_pd(theta);
			             const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_ps(vy,st));
				     _mm512_storeu_pd(row7[i], _mm512_mul_pd(vz,sx));
				     _mm512_storeu_pd(row1[i], _mm512_fmsub_pd(vx,sx,ct));
				     const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				     _mm512_storeu_pd(&row8[i], _mm512_mul_pd(vz,sy));
				     _mm512_storeu_pd(&row2[i], _mm512_fmsub_pd(vx,sy,st));
				     _mm512_storeu_pd(&row4[i], _mm512_fmadd_pd(vy,sx,st));
				     _mm512_storeu_pd(&row5[i], _mm512_fmsub_pd(vy,sy,ct));
                               }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#pragma novector
#endif
                                    for(; i != n; ++i) {
                                        const double x = rx[i];
					const double theta = x*6.2831853071795864769253;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double y = ry[i];
					const double phi = y*6.2831853071795864769253;
					const double z = rz[i];
					const double zz = z+z;
					row9[i]        = 1.0-z;
					const double r  = sqrt(zz);
					const double vx = sin(phi) * r;
					const double vy = cos(phi) * r;
					const double vz = sqrt(2.0-zz);
				        row3[i]        = vx * vz;
					row6[i]        = vy * vz;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double sx = vx * ct - vy * st;
					row1[i]        = vx * sx - ct;
					const double sy = vx * st + vy * ct;
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
			rand_rm9x8_looped_a_nonunroll_zmm8r8( double * __restrict __ATTR_ALIGN__(64) row1,
		                                              double * __restrict __ATTR_ALIGN__(64) row2,
							      double * __restrict __ATTR_ALIGN__(64) row3,
							      double * __restrict __ATTR_ALIGN__(64) row4,
							      double * __restrict __ATTR_ALIGN__(64) row5,
							      double * __restrict __ATTR_ALIGN__(64) row6,
							      double * __restrict __ATTR_ALIGN__(64) row7,
							      double * __restrict __ATTR_ALIGN__(64) row8,
							      double * __restrict __ATTR_ALIGN__(64) row9,
							      const double * __restrict  __ATTR_ALIGN__(64) rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  __ATTR_ALIGN__(64) ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  __ATTR_ALIGN__(64) rz, // ramdomly normally distributed vector [0,1]
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
                               row1 = (double*)__builtin_assume_aligned(row1,64);
			       row2 = (double*)__builtin_assume_aligned(row2,64);
			       row3 = (double*)__builtin_assume_aligned(row3,64);
			       row4 = (double*)__builtin_assume_aligned(row4,64);
			       row5 = (double*)__builtin_assume_aligned(row5,64);
			       row6 = (double*)__builtin_assume_aligned(row6,64);
			       row7 = (double*)__builtin_assume_aligned(row7,64);
			       row8 = (double*)__builtin_assume_aligned(row8,64);
			       row9 = (double*)__builtin_assume_aligned(row9,64);
			       rx = (double*)__builtin_assume_aligned(rx,64);
			       ry = (double*)__builtin_assume_aligned(ry,64);
			       rz = (double*)__builtin_assume_aligned(rz,64);
			       
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif                         for(i = 0; i != ROUND_TO_EIGHT(n,8); i += 8) {
                                     _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0);
				     register const __m512d vrx   = _mm512_load_pd(&rx[i]);
				     register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				     _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0);
				     register const __m512d vry   = _mm512_load_pd(&ry[i]);
				     register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				     _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0);
				     register const __m512d vrz   = _mm512_load_pd(&rz[i]);
				     register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				     _mm512_store_pd(&row9[i], _mm512_sub_ps(v8_1,z));
				     const register __m512d r     = _mm512_sqrt_pd(z);
				     const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				     const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				     const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
				     _mm512_store_pd(row6[i], _mm512_mul_pd(vy,vz));
				     _mm512_store_pd(row3[i], _mm512_mul_pd(vx,vz));
				     const register __m512d st    = _mm512_sin_pd(theta);
			             const register __m512d ct    = _mm512_cos_pd(theta);
			             const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_ps(vy,st));
				     _mm512_store_pd(row7[i], _mm512_mul_pd(vz,sx));
				     _mm512_store_pd(row1[i], _mm512_fmsub_pd(vx,sx,ct));
				     const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				     _mm512_store_pd(&row8[i], _mm512_mul_pd(vz,sy));
				     _mm512_store_pd(&row2[i], _mm512_fmsub_pd(vx,sy,st));
				     _mm512_store_pd(&row4[i], _mm512_fmadd_pd(vy,sx,st));
				     _mm512_store_pd(&row5[i], _mm512_fmsub_pd(vy,sy,ct));
                               }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#pragma novector
#endif
                                    for(; i != n; ++i) {
                                        const double x = rx[i];
					const double theta = x*6.2831853071795864769253;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double y = ry[i];
					const double phi = y*6.2831853071795864769253;
					const double z = rz[i];
					const double zz = z+z;
					row9[i]        = 1.0-z;
					const double r  = sqrt(zz);
					const double vx = sin(phi) * r;
					const double vy = cos(phi) * r;
					const double vz = sqrt(2.0-zz);
				        row3[i]        = vx * vz;
					row6[i]        = vy * vz;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double sx = vx * ct - vy * st;
					row1[i]        = vx * sx - ct;
					const double sy = vx * st + vy * ct;
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
                      rand_rm9x8_looped_u_unroll4x_zmm8r8(    double * __restrict row1,
		                                              double * __restrict row2,
							      double * __restrict row3,
							      double * __restrict row4,
							      double * __restrict row5,
							      double * __restrict  row6,
							      double * __restrict  row7,
							      double * __restrict  row8,
							      double * __restrict  row9,
							      const double * __restrict  rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  rz, // ramdomly normally distributed vector [0,1]
							      const int32_t n) {

                               if(__builtin_expect(0<=n,0)) {return;}
			       int32_t i;
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			      for(i = 0; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&rx[i+0],   _MM_HINT_T0);
				  _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0)
                                  _mm_prefetch((const char*)&rx[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&rx[i+24],_MM_HINT_T0);
                                  register const __m512d vrx0   = _mm512_loadu_pd(&rx[i+0]);
				  register const __m512d vrx1   = _mm512_loadu_pd(&rx[i+8]);
				  register const __m512d vrx2   = _mm512_loadu_pd(&rx[i+16]);
				  register const __m512d vrx3   = _mm512_loadu_pd(&rx[i+24]);
				  register const __m512d theta0 = _mm512_mul_pd(vrx0,v8_2pi);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  register const __m512d theta2 = _mm512_mul_pd(vrx2,v8_2pi);
				  register const __m512d theta3 = _mm512_mul_pd(vrx3,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],   _MM_HINT_T0);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0)
                                  _mm_prefetch((const char*)&ry[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&ry[i+24],_MM_HINT_T0);
				  register const __m512d vry0   = _mm512_loadu_pd(&ry[i+0]);
				  register const __m512d vry1    = _mm512_loadu_pd(&ry[i+8]);
				  register const __m512d vry2   = _mm512_loadu_pd(&ry[i+16]);
				  register const __m512d vry3   = _mm512_loadu_pd(&ry[i+24]);
				  register const __m512d phi0   = _mm512_mul_ps(vry0,v8_2pi);
				  register const __m512d phi1   = _mm512_mul_ps(vry1,v8_2pi);
				  register const __m512d phi2   = _mm512_mul_ps(vry2,v8_2pi);
				  register const __m512d phi3   = _mm512_mul_ps(vry3,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],   _MM_HINT_T0);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0)
                                  _mm_prefetch((const char*)&rz[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&rz[i+24],_MM_HINT_T0);
				  register const __m512d vrz0   = _mm512_loadu_pd(&rz[i+0]);
				  register const __m512d vrz1   = _mm512_loadu_pd(&rz[i+8]);
				  register const __m512d vrz2   = _mm512_loadu_pd(&rz[i+16]);
				  register const __m512d vrz3   = _mm512_loadu_pd(&rz[i+24]);
				  register const __m512d z0     = _mm512_mul_pd(vrz0,vrz0);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  register const __m512d z2     = _mm512_mul_pd(vrz2,vrz2);
				  register const __m512d z3     = _mm512_mul_pd(vrz3,vrz3);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z0);
				  _mm512_storeu_pd(&row9[i+8],_mm512_sub_pd(v8_1,z1);
				  _mm512_storeu_pd(&row9[i+16],_mm512_sub_pd(v8_1,z2);
				  _mm512_storeu_pd(&row9[i+24],_mm512_sub_pd(v8_1,z3);
				  register const __m512d r0     = _mm512_sqrt_pd(z0);
				  register const __m512d r1     = _mm512_sqrt_pd(z1);
				  register const __m512d r2     = _mm512_sqrt_pd(z2);
				  register const __m512d r3     = _mm512_sqrt_pd(z3);
				  register const __m512d vx0    = _mm512_mul_pd(r0,_mm512_sin_pd(phi0));
				  register const __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  register const __m512d vx2    = _mm512_mul_pd(r2,_mm512_sin_pd(phi2));
				  register const __m512d vx3    = _mm512_mul_pd(r3,_mm512_sin_pd(phi3));
				  register const __m512d vy0    = _mm512_mul_pd(r0,_mm512_cos_pd(phi0));
				  register const __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  register const __m512d vy2    = _mm512_mul_pd(r2,_mm512_cos_pd(phi2));
				  register const __m512d vy3    = _mm512_mul_pd(r3,_mm512_cos_pd(phi3));
				  register const __m512d vz0    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z0));
				  register const __m512d vz1    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z1));
				  register const __m512d vz2    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z2));
				  register const __m512d vz3    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z3));
				  _mm512_storeu_pd(&row6[i+0],  _mm512_mul_pd(vy0,vz0));
				  _mm512_storeu_pd(&row6[i+8],  _mm512_mul_pd(vy1,vz1));
				  _mm512_storeu_pd(&row6[i+16], _mm512_mul_pd(vy2,vz2));
				  _mm512_storeu_pd(&row6[i+24], _mm512_mul_pd(vy3,vz3));
				  _mm512_storeu_pd(&row3[i+0],  _mm512_mul_pd(vx0,vz0));
				  _mm512_storeu_pd(&row3[i+8], _mm512_mul_pd(vx1,vz1));
				  _mm512_storeu_pd(&row3[i+16], _mm512_mul_pd(vx2,vz2));
				  _mm512_storeu_pd(&row3[i+24], _mm512_mul_pd(vx3,vz3));
				  register const __m512d st0     = _mm512_sin_pd(theta0);
				  register const __m512d st1     = _mm512_sin_pd(theta1);
				  register const __m512d st2     = _mm512_sin_pd(theta2);
				  register const __m512d st3     = _mm512_sin_pd(theta3);
				  register const __m512d ct0     = _mm512_cos_pd(theta0);
				  register const __m512d ct1     = _mm512_cos_pd(theta1);
				  register const __m512d ct2     = _mm512_cos_pd(theta2);
				  register const __m512d ct3     = _mm512_cos_pd(theta3);
				  register const __m512d sx0     = _mm512_fmsub_pd(vx0,ct0,_mm512_mul_pd(vy0,st0));
				  register const __m512d sx1     = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  register const __m512d sx2     = _mm512_fmsub_pd(vx2,ct2,_mm512_mul_pd(vy2,st2));
				  register const __m512d sx3     = _mm512_fmsub_pd(vx3,ct3,_mm512_mul_pd(vy3,st3));
				  _mm512_storeu_pd(&row7[i+0],  _mm512_mul_pd(vz0,sx0));
				  _mm512_storeu_pd(&row7[i+8 ], _mm512_mul_pd(vz1,sx1));
				  _mm512_storeu_pd(&row7[i+16], _mm512_mul_pd(vz2,sx2));
				  _mm512_storeu_pd(&row7[i+24], _mm512_mul_pd(vz3,sx3));
				  _mm512_storeu_pd(&row1[i+0],  _mm512_fmsub_pd(vx0,sx0,ct0));
				  _mm512_storeu_pd(&row1[i+8 ], _mm512_fmsub_pd(vx1,sx1,ct1));
				  _mm512_storeu_pd(&row1[i+16], _mm512_fmsub_pd(vx2,sx2,ct2));
				  _mm512_storeu_pd(&row1[i+24], _mm512_fmsub_pd(vx3,sx3,ct3));
				  register const __m512d sy0     = _mm512_fmadd_pd(vx0,st0,_mm512_mul_pd(vy0,ct0));
				  register const __m512d sy1     = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  register const __m512d sy2     = _mm512_fmadd_pd(vx2,st2,_mm512_mul_pd(vy2,ct2));
				  register const __m512d sy3     = _mm512_fmadd_pd(vx3,st3,_mm512_mul_pd(vy3,ct3));
				  _mm512_storeu_pd(&row8[i+0],  _mm512_mul_pd(vz0,sy0));
				  _mm512_storeu_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_storeu_pd(&row8[i+16], _mm512_mul_pd(vz2,sy2));
				  _mm512_storeu_pd(&row8[i+24], _mm512_mul_pd(vz3,sy3));
				  _mm512_storeu_pd(&row2[i+0],  _mm512_fmsub_pd(vx0,sy0,st0));
				  _mm512_storeu_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_storeu_pd(&row2[i+16], _mm512_fmsub_pd(vx2,sy2,st2));
				  _mm512_storeu_pd(&row2[i+24], _mm512_fmsub_pd(vx3,sy3,st3));
				  _mm512_storeu_pd(&row4[i+0],  _mm512_fmadd_pd(vy0,sx0,st0));
				  _mm512_storeu_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_storeu_pd(&row4[i+16], _mm512_fmadd_pd(vy2,sx2,st2));
				  _mm512_storeu_pd(&row4[i+24], _mm512_fmadd_pd(vy3,sx3,st3));
				  _mm512_storeu_pd(&row5[i+0],  _mm512_fmsub_pd(vy0,sy0,ct0));
				  _mm512_storeu_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
				  _mm512_storeu_pd(&row5[i+16], _mm512_fmsub_pd(vy2,sy2,ct2));
				  _mm512_storeu_pd(&row5[i+24], _mm512_fmsub_pd(vy3,sy3,ct3));
				  
#else
                                  _mm_prefetch((const char*)&rx[i+0],_MM_HINT_T0);
				  register const __m512d vrx   = _mm512_loadu_pd(&rx[i]);
				  register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],_MM_HINT_T0);
				  register const __m512d vry   = _mm512_loadu_pd(&ry[i]);
				  register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],_MM_HINT_T0);
				  register const __m512d vrz   = _mm512_loadu_pd(&rz[i]);
				  register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z));
				  const register __m512d r     = _mm512_sqrt_pd(z);
				  const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				  const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				  const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
				  _mm512_storeu_pd(row6[i+0], _mm512_mul_pd(vy,vz));
				  _mm512_storeu_pd(row3[i+0], _mm512_mul_pd(vx,vz));
				  const register __m512d st    = _mm512_sin_pd(theta);
			          const register __m512d ct    = _mm512_cos_pd(theta);
			          const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
				  _mm512_storeu_pd(row7[i+0], _mm512_mul_pd(vz,sx));
				  _mm512_storeu_pd(row1[i+0], _mm512_fmsub_pd(vx,sx,ct));
				  const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				  _mm512_storeu_pd(&row8[i+0], _mm512_mul_pd(vz,sy));
				  _mm512_storeu_pd(&row2[i+0], _mm512_fmsub_pd(vx,sy,st));
				  _mm512_storeu_pd(&row4[i+0], _mm512_fmadd_pd(vy,sx,st));
				  _mm512_storeu_pd(&row5[i+0], _mm512_fmsub_pd(vy,sy,ct));
				  _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0);
				  register const __m512d vrx1   = _mm512_loadu_pd(&rx[i+8]);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0);
				  register const __m512d vry1   = _mm512_loadu_pd(&ry[i+8]);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0);
				  register const __m512d vrz1   = _mm512_loadu_pd(&rz[i+8]);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  _mm512_storeu_pd(&row9[i+8], _mm512_sub_pd(v8_1,z1));
				  const register __m512d r1     = _mm512_sqrt_pd(z1);
				  const register __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  const register __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  const register __m512d vz1    = _mm512_sqrt_pd(_mm512_sub_ps(v8_2,z1));
				  _mm512_storeu_pd(row6[i+8], _mm512_mul_pd(vy1,vz1));
				  _mm512_storeu_pd(row3[i+8], _mm512_mul_pd(vx1,vz1));
				  const register __m512d st1    = _mm512_sin_pd(theta1);
			          const register __m512d ct1    = _mm512_cos_pd(theta1);
			          const register __m512d sx1    = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  _mm512_storeu_pd(row7[i+8], _mm512_mul_pd(vz1,sx1));
				  _mm512_storeu_pd(row1[i+8], _mm512_fmsub_pd(vx1,sx1,ct1));
				  const register __m512d sy1    = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  _mm512_storeu_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_storeu_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_storeu_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_storeu_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
				  _mm_prefetch((const char*)&rx[i+16],_MM_HINT_T0);
				  register const __m512d vrx2   = _mm512_loadu_pd(&rx[i+16]);
				  register const __m512d theta2 = _mm512_mul_pd(vrx2,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+16],_MM_HINT_T0);
				  register const __m512d vry2   = _mm512_loadu_pd(&ry[i+16]);
				  register const __m512d phi2   = _mm512_mul_pd(vry2,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+16],_MM_HINT_T0);
				  register const __m512d vrz2   = _mm512_loadu_pd(&rz[i+16]);
				  register const __m512d z2     = _mm512_mul_pd(vrz2,vrz2);
				  _mm512_storeu_pd(&row9[i+16], _mm512_sub_pd(v8_1,z2));
				  const register __m512d r2     = _mm512_sqrt_pd(z2);
				  const register __m512d vx2    = _mm512_mul_pd(r2,_mm512_sin_pd(phi2));
				  const register __m512d vy2    = _mm512_mul_pd(r2,_mm512_cos_pd(phi2));
				  const register __m512d vz2    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z2));
				  _mm512_storeu_pd(row6[i+16], _mm512_mul_pd(vy2,vz2));
				  _mm512_storeu_pd(row3[i+16], _mm512_mul_pd(vx2,vz2));
				  const register __m512d st2    = _mm512_sin_pd(theta2);
			          const register __m512d ct2    = _mm512_cos_pd(theta2);
			          const register __m512d sx2    = _mm512_fmsub_pd(vx2,ct2,_mm512_mul_pd(vy2,st2));
				  _mm512_storeu_pd(row7[i+16], _mm512_mul_pd(vz2,sx2));
				  _mm512_storeu_pd(row1[i+16], _mm512_fmsub_pd(vx2,sx2,ct2));
				  const register __m512d sy2    = _mm512_fmadd_pd(vx2,st2,_mm512_mul_pd(vy2,ct2));
				  _mm512_storeu_pd(&row8[i+16], _mm512_mul_pd(vz2,sy2));
				  _mm512_storeu_pd(&row2[i+16], _mm512_fmsub_pd(vx2,sy2,st2));
				  _mm512_storeu_pd(&row4[i+16], _mm512_fmadd_pd(vy2,sx2,st2));
				  _mm512_storeu_pd(&row5[i+16], _mm512_fmsub_pd(vy2,sy2,ct2));
				  _mm_prefetch((const char*)&rx[i+24],_MM_HINT_T0);
				  register const __m512d vrx3   = _mm512_loadu_pd(&rx[i+24]);
				  register const __m512d theta3 = _mm512_mul_pd(vrx3,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+24],_MM_HINT_T0);
				  register const __m512d vry3   = _mm512_loadu_pd(&ry[i+48]);
				  register const __m512d phi3   = _mm512_mul_pd(vry3,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+24],_MM_HINT_T0);
				  register const __m512d vrz3   = _mm512_loadu_pd(&rz[i+24]);
				  register const __m512d z3     = _mm512_mul_pd(vrz3,vrz3);
				  _mm512_storeu_pd(&row9[i+24], _mm512_sub_pd(v8_1,z3));
				  const register __m512d r3     = _mm512_sqrt_pd(z3);
				  const register __m512d vx3    = _mm512_mul_pd(r3,_mm512_sin_pd(phi3));
				  const register __m512d vy3    = _mm512_mul_pd(r3,_mm512_cos_pd(phi3));
				  const register __m512d vz3    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z1));
				  _mm512_storeu_pd(row6[i+24], _mm512_mul_pd(vy3,vz3));
				  _mm512_storeu_pd(row3[i+24], _mm512_mul_pd(vx3,vz3));
				  const register __m512d st3    = _mm512_sin_pd(theta3);
			          const register __m512d ct3    = _mm512_cos_pd(theta3);
			          const register __m512d sx3    = _mm512_fmsub_pd(vx3,ct3,_mm512_mul_pd(vy3,st3));
				  _mm512_storeu_ps(row7[i+24], _mm512_mul_pd(vz3,sx3));
				  _mm512_storeu_ps(row1[i+24], _mm512_fmsub_pd(vx3,sx3,ct3));
				  const register __m512d sy3    = _mm512_fmadd_pd(vx3,st3,_mm512_mul_pd(vy3,ct3));
				  _mm512_storeu_ps(&row8[i+24], _mm512_mul_pd(vz3,sy3));
				  _mm512_storeu_ps(&row2[i+24], _mm512_fmsub_pd(vx3,sy3,st3));
				  _mm512_storeu_ps(&row4[i+24], _mm512_fmadd_pd(vy3,sx3,st3));
				  _mm512_storeu_ps(&row5[i+24], _mm512_fmsub_pd(vy3,sy3,ct3));
#endif
                              }

			      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&rx[i+0], _MM_HINT_T0);
				  _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0)
                                  register const __m512d vrx0   = _mm512_loadu_pd(&rx[i+0]);
				  register const __m512d vrx1   = _mm512_loadu_pd(&rx[i+8]);
				  register const __m512d theta0 = _mm512_mul_pd(vrx0,v8_2pi);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0], _MM_HINT_T0);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0)
                                  register const __m512d vry0   = _mm512_loadu_pd(&ry[i+0]);
				  register const __m512d vry1   = _mm512_loadu_pd(&ry[i+8]);
				  register const __m512d phi0   = _mm512_mul_pd(vry0,v8_2pi);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0], _MM_HINT_T0);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0)
                                  register const __m512d vrz0   = _mm512_loadu_pd(&rz[i+0]);
				  register const __m512d vrz1   = _mm512_loadu_pd(&rz[i+8]);
				  register const __m512d z0     = _mm512_mul_pd(vrz0,vrz0);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z0);
				  _mm512_storeu_pd(&row9[i+8],_mm512_sub_pd(v8_1,z1);
				  register const __m512d r0     = _mm512_sqrt_pd(z0);
				  register const __m512d r1     = _mm512_sqrt_pd(z1);
				  register const __m512d vx0    = _mm512_mul_pd(r0,_mm512_sin_pd(phi0));
				  register const __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  register const __m512d vy0    = _mm512_mul_pd(r0,_mm512_cos_pd(phi0));
				  register const __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  register const __m512d vz0    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z0));
				  register const __m512d vz1    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z1));
				  _mm512_storeu_pd(&row6[i+0],  _mm512_mul_pd(vy0,vz0));
				  _mm512_storeu_pd(&row6[i+8], _mm512_mul_pd(vy1,vz1));
				  _mm512_storeu_pd(&row3[i+0],  _mm512_mul_pd(vx0,vz0));
				  _mm512_storeu_pd(&row3[i+8], _mm512_mul_pd(vx1,vz1));
				  register const __m512d st0     = _mm512_sin_pd(theta0);
				  register const __m512d st1     = _mm512_sin_pd(theta1);
				  register const __m512d ct0     = _mm512_cos_pd(theta0);
				  register const __m512d ct1     = _mm512_cos_pd(theta1);
				  register const __m512d sx0     = _mm512_fmsub_pd(vx0,ct0,_mm512_mul_pd(vy0,st0));
				  register const __m512d sx1     = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  _mm512_storeu_pd(&row7[i+0],  _mm512_mul_pd(vz0,sx0));
				  _mm512_storeu_pd(&row7[i+8], _mm512_mul_pd(vz1,sx1));
				  _mm512_storeu_pd(&row1[i+0],  _mm512_fmsub_pd(vx0,sx0,ct0));
				  _mm512_storeu_pd(&row1[i+8], _mm512_fmsub_pd(vx1,sx1,ct1));
				  register const __m512d sy0     = _mm512_fmadd_pd(vx0,st0,_mm512_mul_pd(vy0,ct0));
				  register const __m512d sy1     = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  _mm512_storeu_pd(&row8[i+0],  _mm512_mul_pd(vz0,sy0));
				  _mm512_storeu_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_storeu_pd(&row2[i+0],  _mm512_fmsub_pd(vx0,sy0,st0));
				  _mm512_storeu_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_storeu_pd(&row4[i+0],  _mm512_fmadd_pd(vy0,sx0,st0));
				  _mm512_storeu_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_storeu_pd(&row5[i+0],  _mm512_fmsub_pd(vy0,sy0,ct0));
				  _mm512_storeu_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
				  
#else
                                  _mm_prefetch((const char*)&rx[i+0],_MM_HINT_T0);
				  register const __m512d vrx   = _mm512_loadu_pd(&rx[i]);
				  register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],_MM_HINT_T0);
				  register const __m512d vry   = _mm512_loadu_pd(&ry[i]);
				  register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],_MM_HINT_T0);
				  register const __m512d vrz   = _mm512_loadu_pd(&rz[i]);
				  register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z));
				  const register __m512d r     = _mm512_sqrt_pd(z);
				  const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				  const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				  const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
				  _mm512_storeu_pd(row6[i+0], _mm512_mul_pd(vy,vz));
				  _mm512_storeu_pd(row3[i+0], _mm512_mul_pd(vx,vz));
				  const register __m512d st    = _mm512_sin_pd(theta);
			          const register __m512d ct    = _mm512_cos_pd(theta);
			          const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
				  _mm512_storeu_pd(row7[i+0], _mm512_mul_pd(vz,sx));
				  _mm512_storeu_pd(row1[i+0], _mm512_fmsub_pd(vx,sx,ct));
				  const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				  _mm512_storeu_pd(&row8[i+0], _mm512_mul_pd(vz,sy));
				  _mm512_storeu_pd(&row2[i+0], _mm512_fmsub_pd(vx,sy,st));
				  _mm512_storeu_pd(&row4[i+0], _mm512_fmadd_pd(vy,sx,st));
				  _mm512_storeu_pd(&row5[i+0], _mm512_fmsub_pd(vy,sy,ct));
				  register const __m512d vrx1   = _mm512_loadu_pd(&rx[i+8]);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0);
				  register const __m512d vry1   = _mm512_loadu_pd(&ry[i+8]);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0);
				  register const __m512d vrz1   = _mm512_loadu_pd(&rz[i+8]);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  _mm512_storeu_pd(&row9[i+8], _mm512_sub_pd(v16_1,z1));
				  const register __m512d r1     = _mm512_sqrt_pd(z1);
				  const register __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  const register __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  const register __m512d vz1    = _mm512_sqrt_pd(_mm512_sub_ps(v8_2,z1));
				  _mm512_storeu_pd(row6[i+8], _mm512_mul_pd(vy1,vz1));
				  _mm512_storeu_pd(row3[i+8], _mm512_mul_pd(vx1,vz1));
				  const register __m512d st1    = _mm512_sin_pd(theta1);
			          const register __m512d ct1    = _mm512_cos_pd(theta1);
			          const register __m512d sx1    = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  _mm512_storeu_pd(row7[i+8], _mm512_mul_pd(vz1,sx1));
				  _mm512_storeu_pd(row1[i+8], _mm512_fmsub_pd(vx1,sx1,ct1));
				  const register __m512d sy1    = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  _mm512_storeu_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_storeu_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_storeu_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_storeu_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
#endif
			      }

			      for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&rx[i+0], _MM_HINT_T0);
				  register const __m512d vrx0   = _mm512_loadu_pd(&rx[i+0]);
				  register const __m512d theta0 = _mm512_mul_pd(vrx0,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0], _MM_HINT_T0);
				  register const __m512d vry0   = _mm512_loadu_pd(&ry[i+0]);
				  register const __m512d phi0   = _mm512_mul_pd(vry0,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0], _MM_HINT_T0);
				  register const __m512d vrz0   = _mm512_loadu_pd(&rz[i+0]);
				  register const __m512d z0     = _mm512_mul_pd(vrz0,vrz0);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z0);
				  register const __m512d r0     = _mm512_sqrt_pd(z0);
				  register const __m512d vx0    = _mm512_mul_pd(r0,_mm512_sin_pd(phi0));
				  register const __m512d vy0    = _mm512_mul_pd(r0,_mm512_cos_pd(phi0));
				  register const __m512d vz0    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z0));
				  _mm512_storeu_pd(&row6[i+0],  _mm512_mul_pd(vy0,vz0));
				  _mm512_storeu_pd(&row3[i+0],  _mm512_mul_pd(vx0,vz0));
				  register const __m512d st0     = _mm512_sin_pd(theta0);
				  register const __m512d ct0     = _mm512_cos_pd(theta0);
				  register const __m512d sx0     = _mm512_fmsub_pd(vx0,ct0,_mm512_mul_pd(vy0,st0));
				  _mm512_storeu_pd(&row7[i+0],  _mm512_mul_pd(vz0,sx0));
				  _mm512_storeu_pd(&row1[i+0],  _mm512_fmsub_pd(vx0,sx0,ct0));
				  register const __m512d sy0     = _mm512_fmadd_pd(vx0,st0,_mm512_mul_pd(vy0,ct0));
				  _mm512_storeu_pd(&row8[i+0],  _mm512_mul_pd(vz0,sy0));
				  _mm512_storeu_pd(&row2[i+0],  _mm512_fmsub_pd(vx0,sy0,st0));
				  _mm512_storeu_pd(&row4[i+0],  _mm512_fmadd_pd(vy0,sx0,st0));
				  _mm512_storeu_pd(&row5[i+0],  _mm512_fmsub_pd(vy0,sy0,ct0));
				 
#else
                                   _mm_prefetch((const char*)&rx[i+0],_MM_HINT_T0);
				  register const __m512d vrx   = _mm512_loadu_pd(&rx[i]);
				  register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],_MM_HINT_T0);
				  register const __m512d vry   = _mm512_loadu_pd(&ry[i]);
				  register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],_MM_HINT_T0);
				  register const __m512d vrz   = _mm512_loadu_pd(&rz[i]);
				  register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z));
				  const register __m512d r     = _mm512_sqrt_pd(z);
				  const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				  const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				  const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_ps(v8_2,z));
				  _mm512_storeu_pd(row6[i+0], _mm512_mul_pd(vy,vz));
				  _mm512_storeu_pd(row3[i+0], _mm512_mul_pd(vx,vz));
				  const register __m512d st    = _mm512_sin_pd(theta);
			          const register __m512d ct    = _mm512_cos_pd(theta);
			          const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
				  _mm512_storeu_pd(row7[i+0], _mm512_mul_pd(vz,sx));
				  _mm512_storeu_pd(row1[i+0], _mm512_fmsub_pd(vx,sx,ct));
				  const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				  _mm512_storeu_pd(&row8[i+0], _mm512_mul_pd(vz,sy));
				  _mm512_storeu_pd(&row2[i+0], _mm512_fmsub_pd(vx,sy,st));
				  _mm512_storeu_pd(&row4[i+0], _mm512_fmadd_pd(vy,sx,st));
				  _mm512_storeu_pd(&row5[i+0], _mm512_fmsub_pd(vy,sy,ct));
#endif
			      }

#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#pragma novector
#endif
                                for(; (i+0) < n; i += 1) {
                                         const double x = rx[i];
					const double theta = x*6.2831853071795864769253;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double y = ry[i];
					const double phi = y*6.2831853071795864769253;
					const double z = rz[i];
					const double zz = z+z;
					row9[i]        = 1.0-z;
					const double r  = sqrt(zz);
					const double vx = sin(phi) * r;
					const double vy = cos(phi) * r;
					const double vz = sqrt(2.0-zz);
				        row3[i]        = vx * vz;
					row6[i]        = vy * vz;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double sx = vx * ct - vy * st;
					row1[i]        = vx * sx - ct;
					const double sy = vx * st + vy * ct;
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
                      rand_rm9x8_looped_a_unroll4x_zmm8r8(    double * __restrict __ATTR_ALIGN__(64) row1,
		                                              double * __restrict __ATTR_ALIGN__(64) row2,
							      double * __restrict __ATTR_ALIGN__(64) row3,
							      double * __restrict __ATTR_ALIGN__(64) row4,
							      double * __restrict __ATTR_ALIGN__(64) row5,
							      double * __restrict __ATTR_ALIGN__(64) row6,
							      double * __restrict __ATTR_ALIGN__(64) row7,
							      double * __restrict __ATTR_ALIGN__(64) row8,
							      double * __restrict __ATTR_ALIGN__(64) row9,
							      const double * __restrict __ATTR_ALIGN__(64) rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict __ATTR_ALIGN__(64)  ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict __ATTR_ALIGN__(64) rz, // ramdomly normally distributed vector [0,1]
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
                               row1 = (double*)__builtin_assume_aligned(row1,64);
			       row2 = (double*)__builtin_assume_aligned(row2,64);
			       row3 = (double*)__builtin_assume_aligned(row3,64);
			       row4 = (double*)__builtin_assume_aligned(row4,64);
			       row5 = (double*)__builtin_assume_aligned(row5,64);
			       row6 = (double*)__builtin_assume_aligned(row6,64);
			       row7 = (double*)__builtin_assume_aligned(row7,64);
			       row8 = (double*)__builtin_assume_aligned(row8,64);
			       row9 = (double*)__builtin_assume_aligned(row9,64);
			       rx = (double*)__builtin_assume_aligned(rx,64);
			       ry = (double*)__builtin_assume_aligned(ry,64);
			       rz = (double*)__builtin_assume_aligned(rz,64);
			       
#endif
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma code_align(32)
#endif
			      for(i = 0; (i+31) < n; i += 32) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&rx[i+0],   _MM_HINT_T0);
				  _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0)
                                  _mm_prefetch((const char*)&rx[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&rx[i+24],_MM_HINT_T0);
                                  register const __m512d vrx0   = _mm512_load_pd(&rx[i+0]);
				  register const __m512d vrx1   = _mm512_load_pd(&rx[i+8]);
				  register const __m512d vrx2   = _mm512_load_pd(&rx[i+16]);
				  register const __m512d vrx3   = _mm512_load_pd(&rx[i+24]);
				  register const __m512d theta0 = _mm512_mul_pd(vrx0,v8_2pi);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  register const __m512d theta2 = _mm512_mul_pd(vrx2,v8_2pi);
				  register const __m512d theta3 = _mm512_mul_pd(vrx3,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],   _MM_HINT_T0);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0)
                                  _mm_prefetch((const char*)&ry[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&ry[i+24],_MM_HINT_T0);
				  register const __m512d vry0   = _mm512_load_pd(&ry[i+0]);
				  register const __m512d vry1    = _mm512_load_pd(&ry[i+8]);
				  register const __m512d vry2   = _mm512_load_pd(&ry[i+16]);
				  register const __m512d vry3   = _mm512_load_pd(&ry[i+24]);
				  register const __m512d phi0   = _mm512_mul_pd(vry0,v8_2pi);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  register const __m512d phi2   = _mm512_mul_pd(vry2,v8_2pi);
				  register const __m512d phi3   = _mm512_mul_pd(vry3,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],   _MM_HINT_T0);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0)
                                  _mm_prefetch((const char*)&rz[i+16],_MM_HINT_T0);
				  _mm_prefetch((const char*)&rz[i+24],_MM_HINT_T0);
				  register const __m512d vrz0   = _mm512_load_pd(&rz[i+0]);
				  register const __m512d vrz1   = _mm512_load_pd(&rz[i+8]);
				  register const __m512d vrz2   = _mm512_load_pd(&rz[i+16]);
				  register const __m512d vrz3   = _mm512_load_pd(&rz[i+24]);
				  register const __m512d z0     = _mm512_mul_pd(vrz0,vrz0);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  register const __m512d z2     = _mm512_mul_pd(vrz2,vrz2);
				  register const __m512d z3     = _mm512_mul_pd(vrz3,vrz3);
				  _mm512_store_pd(&row9[i+0], _mm512_sub_pd(v8_1,z0);
				  _mm512_store_pd(&row9[i+8],_mm512_sub_pd(v8_1,z1);
				  _mm512_store_pd(&row9[i+16],_mm512_sub_pd(v8_1,z2);
				  _mm512_store_pd(&row9[i+24],_mm512_sub_pd(v8_1,z3);
				  register const __m512d r0     = _mm512_sqrt_pd(z0);
				  register const __m512d r1     = _mm512_sqrt_pd(z1);
				  register const __m512d r2     = _mm512_sqrt_pd(z2);
				  register const __m512d r3     = _mm512_sqrt_pd(z3);
				  register const __m512d vx0    = _mm512_mul_pd(r0,_mm512_sin_pd(phi0));
				  register const __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  register const __m512d vx2    = _mm512_mul_pd(r2,_mm512_sin_pd(phi2));
				  register const __m512d vx3    = _mm512_mul_pd(r3,_mm512_sin_pd(phi3));
				  register const __m512d vy0    = _mm512_mul_pd(r0,_mm512_cos_pd(phi0));
				  register const __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  register const __m512d vy2    = _mm512_mul_pd(r2,_mm512_cos_pd(phi2));
				  register const __m512d vy3    = _mm512_mul_pd(r3,_mm512_cos_pd(phi3));
				  register const __m512d vz0    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z0));
				  register const __m512d vz1    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z1));
				  register const __m512d vz2    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z2));
				  register const __m512d vz3    = _mm512_sqrt_ps(_mm512_sub_pd(v8_2,z3));
				  _mm512_store_pd(&row6[i+0],  _mm512_mul_pd(vy0,vz0));
				  _mm512_store_pd(&row6[i+8],  _mm512_mul_pd(vy1,vz1));
				  _mm512_store_pd(&row6[i+16], _mm512_mul_pd(vy2,vz2));
				  _mm512_store_pd(&row6[i+24], _mm512_mul_pd(vy3,vz3));
				  _mm512_store_pd(&row3[i+0],  _mm512_mul_pd(vx0,vz0));
				  _mm512_store_pd(&row3[i+8], _mm512_mul_pd(vx1,vz1));
				  _mm512_store_pd(&row3[i+16], _mm512_mul_pd(vx2,vz2));
				  _mm512_store_pd(&row3[i+24], _mm512_mul_pd(vx3,vz3));
				  register const __m512d st0     = _mm512_sin_pd(theta0);
				  register const __m512d st1     = _mm512_sin_pd(theta1);
				  register const __m512d st2     = _mm512_sin_pd(theta2);
				  register const __m512d st3     = _mm512_sin_pd(theta3);
				  register const __m512d ct0     = _mm512_cos_pd(theta0);
				  register const __m512d ct1     = _mm512_cos_pd(theta1);
				  register const __m512d ct2     = _mm512_cos_pd(theta2);
				  register const __m512d ct3     = _mm512_cos_pd(theta3);
				  register const __m512d sx0     = _mm512_fmsub_pd(vx0,ct0,_mm512_mul_pd(vy0,st0));
				  register const __m512d sx1     = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  register const __m512d sx2     = _mm512_fmsub_pd(vx2,ct2,_mm512_mul_pd(vy2,st2));
				  register const __m512d sx3     = _mm512_fmsub_pd(vx3,ct3,_mm512_mul_pd(vy3,st3));
				  _mm512_storeu_pd(&row7[i+0],  _mm512_mul_pd(vz0,sx0));
				  _mm512_storeu_pd(&row7[i+8 ], _mm512_mul_pd(vz1,sx1));
				  _mm512_storeu_pd(&row7[i+16], _mm512_mul_pd(vz2,sx2));
				  _mm512_storeu_pd(&row7[i+24], _mm512_mul_pd(vz3,sx3));
				  _mm512_storeu_pd(&row1[i+0],  _mm512_fmsub_pd(vx0,sx0,ct0));
				  _mm512_storeu_pd(&row1[i+8 ], _mm512_fmsub_pd(vx1,sx1,ct1));
				  _mm512_storeu_pd(&row1[i+16], _mm512_fmsub_pd(vx2,sx2,ct2));
				  _mm512_storeu_pd(&row1[i+24], _mm512_fmsub_pd(vx3,sx3,ct3));
				  register const __m512d sy0     = _mm512_fmadd_pd(vx0,st0,_mm512_mul_pd(vy0,ct0));
				  register const __m512d sy1     = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  register const __m512d sy2     = _mm512_fmadd_pd(vx2,st2,_mm512_mul_pd(vy2,ct2));
				  register const __m512d sy3     = _mm512_fmadd_pd(vx3,st3,_mm512_mul_pd(vy3,ct3));
				  _mm512_store_pd(&row8[i+0],  _mm512_mul_pd(vz0,sy0));
				  _mm512_store_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_store_pd(&row8[i+16], _mm512_mul_pd(vz2,sy2));
				  _mm512_store_pd(&row8[i+24], _mm512_mul_pd(vz3,sy3));
				  _mm512_store_pd(&row2[i+0],  _mm512_fmsub_pd(vx0,sy0,st0));
				  _mm512_store_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_store_pd(&row2[i+16], _mm512_fmsub_pd(vx2,sy2,st2));
				  _mm512_store_pd(&row2[i+24], _mm512_fmsub_pd(vx3,sy3,st3));
				  _mm512_store_pd(&row4[i+0],  _mm512_fmadd_pd(vy0,sx0,st0));
				  _mm512_store_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_store_pd(&row4[i+16], _mm512_fmadd_pd(vy2,sx2,st2));
				  _mm512_store_pd(&row4[i+24], _mm512_fmadd_pd(vy3,sx3,st3));
				  _mm512_store_pd(&row5[i+0],  _mm512_fmsub_pd(vy0,sy0,ct0));
				  _mm512_storeu_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
				  _mm512_store_pd(&row5[i+16], _mm512_fmsub_pd(vy2,sy2,ct2));
				  _mm512_store_pd(&row5[i+24], _mm512_fmsub_pd(vy3,sy3,ct3));
				  
#else
                                  _mm_prefetch((const char*)&rx[i+0],_MM_HINT_T0);
				  register const __m512d vrx   = _mm512_load_pd(&rx[i]);
				  register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],_MM_HINT_T0);
				  register const __m512d vry   = _mm512_load_pd(&ry[i]);
				  register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],_MM_HINT_T0);
				  register const __m512d vrz   = _mm512_load_pd(&rz[i]);
				  register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				  _mm512_store_pd(&row9[i+0], _mm512_sub_pd(v8_1,z));
				  const register __m512d r     = _mm512_sqrt_pd(z);
				  const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				  const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				  const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
				  _mm512_store_pd(row6[i+0], _mm512_mul_pd(vy,vz));
				  _mm512_store_pd(row3[i+0], _mm512_mul_pd(vx,vz));
				  const register __m512d st    = _mm512_sin_pd(theta);
			          const register __m512d ct    = _mm512_cos_pd(theta);
			          const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
				  _mm512_store_pd(row7[i+0], _mm512_mul_pd(vz,sx));
				  _mm512_store_pd(row1[i+0], _mm512_fmsub_pd(vx,sx,ct));
				  const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				  _mm512_store_pd(&row8[i+0], _mm512_mul_pd(vz,sy));
				  _mm512_store_pd(&row2[i+0], _mm512_fmsub_pd(vx,sy,st));
				  _mm512_store_pd(&row4[i+0], _mm512_fmadd_pd(vy,sx,st));
				  _mm512_store_pd(&row5[i+0], _mm512_fmsub_pd(vy,sy,ct));
				  _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0);
				  register const __m512d vrx1   = _mm512_load_pd(&rx[i+8]);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0);
				  register const __m512d vry1   = _mm512_load_pd(&ry[i+8]);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0);
				  register const __m512d vrz1   = _mm512_load_pd(&rz[i+8]);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  _mm512_store_pd(&row9[i+8], _mm512_sub_pd(v8_1,z1));
				  const register __m512d r1     = _mm512_sqrt_pd(z1);
				  const register __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  const register __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  const register __m512d vz1    = _mm512_sqrt_pd(_mm512_sub_ps(v8_2,z1));
				  _mm512_store_pd(row6[i+8], _mm512_mul_pd(vy1,vz1));
				  _mm512_store_pd(row3[i+8], _mm512_mul_pd(vx1,vz1));
				  const register __m512d st1    = _mm512_sin_pd(theta1);
			          const register __m512d ct1    = _mm512_cos_pd(theta1);
			          const register __m512d sx1    = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  _mm512_store_pd(row7[i+8], _mm512_mul_pd(vz1,sx1));
				  _mm512_store_pd(row1[i+8], _mm512_fmsub_pd(vx1,sx1,ct1));
				  const register __m512d sy1    = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  _mm512_store_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_store_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_store_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_store_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
				  _mm_prefetch((const char*)&rx[i+16],_MM_HINT_T0);
				  register const __m512d vrx2   = _mm512_load_pd(&rx[i+16]);
				  register const __m512d theta2 = _mm512_mul_pd(vrx2,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+16],_MM_HINT_T0);
				  register const __m512d vry2   = _mm512_load_pd(&ry[i+16]);
				  register const __m512d phi2   = _mm512_mul_pd(vry2,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+16],_MM_HINT_T0);
				  register const __m512d vrz2   = _mm512_load_pd(&rz[i+16]);
				  register const __m512d z2     = _mm512_mul_pd(vrz2,vrz2);
				  _mm512_store_pd(&row9[i+16], _mm512_sub_pd(v8_1,z2));
				  const register __m512d r2     = _mm512_sqrt_pd(z2);
				  const register __m512d vx2    = _mm512_mul_pd(r2,_mm512_sin_pd(phi2));
				  const register __m512d vy2    = _mm512_mul_pd(r2,_mm512_cos_pd(phi2));
				  const register __m512d vz2    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z2));
				  _mm512_store_pd(row6[i+16], _mm512_mul_pd(vy2,vz2));
				  _mm512_store_pd(row3[i+16], _mm512_mul_pd(vx2,vz2));
				  const register __m512d st2    = _mm512_sin_pd(theta2);
			          const register __m512d ct2    = _mm512_cos_pd(theta2);
			          const register __m512d sx2    = _mm512_fmsub_pd(vx2,ct2,_mm512_mul_pd(vy2,st2));
				  _mm512_store_pd(row7[i+16], _mm512_mul_pd(vz2,sx2));
				  _mm512_store_pd(row1[i+16], _mm512_fmsub_pd(vx2,sx2,ct2));
				  const register __m512d sy2    = _mm512_fmadd_pd(vx2,st2,_mm512_mul_pd(vy2,ct2));
				  _mm512_store_pd(&row8[i+16], _mm512_mul_pd(vz2,sy2));
				  _mm512_store_pd(&row2[i+16], _mm512_fmsub_pd(vx2,sy2,st2));
				  _mm512_store_pd(&row4[i+16], _mm512_fmadd_pd(vy2,sx2,st2));
				  _mm512_store_pd(&row5[i+16], _mm512_fmsub_pd(vy2,sy2,ct2));
				  _mm_prefetch((const char*)&rx[i+24],_MM_HINT_T0);
				  register const __m512d vrx3   = _mm512_load_pd(&rx[i+24]);
				  register const __m512d theta3 = _mm512_mul_pd(vrx3,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+24],_MM_HINT_T0);
				  register const __m512d vry3   = _mm512_load_pd(&ry[i+48]);
				  register const __m512d phi3   = _mm512_mul_pd(vry3,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+24],_MM_HINT_T0);
				  register const __m512d vrz3   = _mm512_load_pd(&rz[i+24]);
				  register const __m512d z3     = _mm512_mul_pd(vrz3,vrz3);
				  _mm512_store_pd(&row9[i+24], _mm512_sub_pd(v8_1,z3));
				  const register __m512d r3     = _mm512_sqrt_pd(z3);
				  const register __m512d vx3    = _mm512_mul_pd(r3,_mm512_sin_pd(phi3));
				  const register __m512d vy3    = _mm512_mul_pd(r3,_mm512_cos_pd(phi3));
				  const register __m512d vz3    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z1));
				  _mm512_store_pd(row6[i+24], _mm512_mul_pd(vy3,vz3));
				  _mm512_store_pd(row3[i+24], _mm512_mul_pd(vx3,vz3));
				  const register __m512d st3    = _mm512_sin_pd(theta3);
			          const register __m512d ct3    = _mm512_cos_pd(theta3);
			          const register __m512d sx3    = _mm512_fmsub_pd(vx3,ct3,_mm512_mul_pd(vy3,st3));
				  _mm512_store_ps(row7[i+24], _mm512_mul_pd(vz3,sx3));
				  _mm512_store_ps(row1[i+24], _mm512_fmsub_pd(vx3,sx3,ct3));
				  const register __m512d sy3    = _mm512_fmadd_pd(vx3,st3,_mm512_mul_pd(vy3,ct3));
				  _mm512_store_ps(&row8[i+24], _mm512_mul_pd(vz3,sy3));
				  _mm512_store_ps(&row2[i+24], _mm512_fmsub_pd(vx3,sy3,st3));
				  _mm512_store_ps(&row4[i+24], _mm512_fmadd_pd(vy3,sx3,st3));
				  _mm512_store_ps(&row5[i+24], _mm512_fmsub_pd(vy3,sy3,ct3));
#endif
                              }

			      for(; (i+15) < n; i += 16) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&rx[i+0], _MM_HINT_T0);
				  _mm_prefetch((const char*)&rx[i+8],_MM_HINT_T0)
                                  register const __m512d vrx0   = _mm512_load_pd(&rx[i+0]);
				  register const __m512d vrx1   = _mm512_load_pd(&rx[i+8]);
				  register const __m512d theta0 = _mm512_mul_pd(vrx0,v8_2pi);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0], _MM_HINT_T0);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0)
                                  register const __m512d vry0   = _mm512_load_pd(&ry[i+0]);
				  register const __m512d vry1   = _mm512_load_pd(&ry[i+8]);
				  register const __m512d phi0   = _mm512_mul_pd(vry0,v8_2pi);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0], _MM_HINT_T0);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0)
                                  register const __m512d vrz0   = _mm512_load_pd(&rz[i+0]);
				  register const __m512d vrz1   = _mm512_load_pd(&rz[i+8]);
				  register const __m512d z0     = _mm512_mul_pd(vrz0,vrz0);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  _mm512_store_pd(&row9[i+0], _mm512_sub_pd(v8_1,z0);
				  _mm512_store_pd(&row9[i+8],_mm512_sub_pd(v8_1,z1);
				  register const __m512d r0     = _mm512_sqrt_pd(z0);
				  register const __m512d r1     = _mm512_sqrt_pd(z1);
				  register const __m512d vx0    = _mm512_mul_pd(r0,_mm512_sin_pd(phi0));
				  register const __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  register const __m512d vy0    = _mm512_mul_pd(r0,_mm512_cos_pd(phi0));
				  register const __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  register const __m512d vz0    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z0));
				  register const __m512d vz1    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z1));
				  _mm512_store_pd(&row6[i+0],  _mm512_mul_pd(vy0,vz0));
				  _mm512_store_pd(&row6[i+8], _mm512_mul_pd(vy1,vz1));
				  _mm512_store_pd(&row3[i+0],  _mm512_mul_pd(vx0,vz0));
				  _mm512_store_pd(&row3[i+8], _mm512_mul_pd(vx1,vz1));
				  register const __m512d st0     = _mm512_sin_pd(theta0);
				  register const __m512d st1     = _mm512_sin_pd(theta1);
				  register const __m512d ct0     = _mm512_cos_pd(theta0);
				  register const __m512d ct1     = _mm512_cos_pd(theta1);
				  register const __m512d sx0     = _mm512_fmsub_pd(vx0,ct0,_mm512_mul_pd(vy0,st0));
				  register const __m512d sx1     = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  _mm512_store_pd(&row7[i+0],  _mm512_mul_pd(vz0,sx0));
				  _mm512_store_pd(&row7[i+8], _mm512_mul_pd(vz1,sx1));
				  _mm512_storeu_pd(&row1[i+0],  _mm512_fmsub_pd(vx0,sx0,ct0));
				  _mm512_storeu_pd(&row1[i+8], _mm512_fmsub_pd(vx1,sx1,ct1));
				  register const __m512d sy0     = _mm512_fmadd_pd(vx0,st0,_mm512_mul_pd(vy0,ct0));
				  register const __m512d sy1     = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  _mm512_store_pd(&row8[i+0],  _mm512_mul_pd(vz0,sy0));
				  _mm512_store_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_store_pd(&row2[i+0],  _mm512_fmsub_pd(vx0,sy0,st0));
				  _mm512_store_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_store_pd(&row4[i+0],  _mm512_fmadd_pd(vy0,sx0,st0));
				  _mm512_store_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_store_pd(&row5[i+0],  _mm512_fmsub_pd(vy0,sy0,ct0));
				  _mm512_store_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
				  
#else
                                  _mm_prefetch((const char*)&rx[i+0],_MM_HINT_T0);
				  register const __m512d vrx   = _mm512_load_pd(&rx[i]);
				  register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],_MM_HINT_T0);
				  register const __m512d vry   = _mm512_load_pd(&ry[i]);
				  register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],_MM_HINT_T0);
				  register const __m512d vrz   = _mm512_load_pd(&rz[i]);
				  register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				  _mm512_store_pd(&row9[i+0], _mm512_sub_pd(v8_1,z));
				  const register __m512d r     = _mm512_sqrt_pd(z);
				  const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				  const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				  const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z));
				  _mm512_store_pd(row6[i+0], _mm512_mul_pd(vy,vz));
				  _mm512_store_pd(row3[i+0], _mm512_mul_pd(vx,vz));
				  const register __m512d st    = _mm512_sin_pd(theta);
			          const register __m512d ct    = _mm512_cos_pd(theta);
			          const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
				  _mm512_store_pd(row7[i+0], _mm512_mul_pd(vz,sx));
				  _mm512_store_pd(row1[i+0], _mm512_fmsub_pd(vx,sx,ct));
				  const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				  _mm512_store_pd(&row8[i+0], _mm512_mul_pd(vz,sy));
				  _mm512_store_pd(&row2[i+0], _mm512_fmsub_pd(vx,sy,st));
				  _mm512_store_pd(&row4[i+0], _mm512_fmadd_pd(vy,sx,st));
				  _mm512_store_pd(&row5[i+0], _mm512_fmsub_pd(vy,sy,ct));
				  register const __m512d vrx1   = _mm512_load_pd(&rx[i+8]);
				  register const __m512d theta1 = _mm512_mul_pd(vrx1,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+8],_MM_HINT_T0);
				  register const __m512d vry1   = _mm512_load_pd(&ry[i+8]);
				  register const __m512d phi1   = _mm512_mul_pd(vry1,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+8],_MM_HINT_T0);
				  register const __m512d vrz1   = _mm512_load_pd(&rz[i+8]);
				  register const __m512d z1     = _mm512_mul_pd(vrz1,vrz1);
				  _mm512_store_pd(&row9[i+8], _mm512_sub_pd(v16_1,z1));
				  const register __m512d r1     = _mm512_sqrt_pd(z1);
				  const register __m512d vx1    = _mm512_mul_pd(r1,_mm512_sin_pd(phi1));
				  const register __m512d vy1    = _mm512_mul_pd(r1,_mm512_cos_pd(phi1));
				  const register __m512d vz1    = _mm512_sqrt_pd(_mm512_sub_ps(v8_2,z1));
				  _mm512_store_pd(row6[i+8], _mm512_mul_pd(vy1,vz1));
				  _mm512_store_pd(row3[i+8], _mm512_mul_pd(vx1,vz1));
				  const register __m512d st1    = _mm512_sin_pd(theta1);
			          const register __m512d ct1    = _mm512_cos_pd(theta1);
			          const register __m512d sx1    = _mm512_fmsub_pd(vx1,ct1,_mm512_mul_pd(vy1,st1));
				  _mm512_store_pd(row7[i+8], _mm512_mul_pd(vz1,sx1));
				  _mm512_store_pd(row1[i+8], _mm512_fmsub_pd(vx1,sx1,ct1));
				  const register __m512d sy1    = _mm512_fmadd_pd(vx1,st1,_mm512_mul_pd(vy1,ct1));
				  _mm512_store_pd(&row8[i+8], _mm512_mul_pd(vz1,sy1));
				  _mm512_store_pd(&row2[i+8], _mm512_fmsub_pd(vx1,sy1,st1));
				  _mm512_store_pd(&row4[i+8], _mm512_fmadd_pd(vy1,sx1,st1));
				  _mm512_store_pd(&row5[i+8], _mm512_fmsub_pd(vy1,sy1,ct1));
#endif
			      }

			      for(; (i+7) < n; i += 8) {
#if (GMS_INTERLEAVE_SIMD_OPS_SCHEDULE) == 1
                                  _mm_prefetch((const char*)&rx[i+0], _MM_HINT_T0);
				  register const __m512d vrx0   = _mm512_load_pd(&rx[i+0]);
				  register const __m512d theta0 = _mm512_mul_pd(vrx0,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0], _MM_HINT_T0);
				  register const __m512d vry0   = _mm512_load_pd(&ry[i+0]);
				  register const __m512d phi0   = _mm512_mul_pd(vry0,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0], _MM_HINT_T0);
				  register const __m512d vrz0   = _mm512_load_pd(&rz[i+0]);
				  register const __m512d z0     = _mm512_mul_pd(vrz0,vrz0);
				  _mm512_storeu_pd(&row9[i+0], _mm512_sub_pd(v8_1,z0);
				  register const __m512d r0     = _mm512_sqrt_pd(z0);
				  register const __m512d vx0    = _mm512_mul_pd(r0,_mm512_sin_pd(phi0));
				  register const __m512d vy0    = _mm512_mul_pd(r0,_mm512_cos_pd(phi0));
				  register const __m512d vz0    = _mm512_sqrt_pd(_mm512_sub_pd(v8_2,z0));
				  _mm512_store_pd(&row6[i+0],  _mm512_mul_pd(vy0,vz0));
				  _mm512_store_pd(&row3[i+0],  _mm512_mul_pd(vx0,vz0));
				  register const __m512d st0     = _mm512_sin_pd(theta0);
				  register const __m512d ct0     = _mm512_cos_pd(theta0);
				  register const __m512d sx0     = _mm512_fmsub_pd(vx0,ct0,_mm512_mul_pd(vy0,st0));
				  _mm512_store_pd(&row7[i+0],  _mm512_mul_pd(vz0,sx0));
				  _mm512_store_pd(&row1[i+0],  _mm512_fmsub_pd(vx0,sx0,ct0));
				  register const __m512d sy0     = _mm512_fmadd_pd(vx0,st0,_mm512_mul_pd(vy0,ct0));
				  _mm512_store_pd(&row8[i+0],  _mm512_mul_pd(vz0,sy0));
				  _mm512_store_pd(&row2[i+0],  _mm512_fmsub_pd(vx0,sy0,st0));
				  _mm512_store_pd(&row4[i+0],  _mm512_fmadd_pd(vy0,sx0,st0));
				  _mm512_store_pd(&row5[i+0],  _mm512_fmsub_pd(vy0,sy0,ct0));
				 
#else
                                   _mm_prefetch((const char*)&rx[i+0],_MM_HINT_T0);
				  register const __m512d vrx   = _mm512_load_pd(&rx[i]);
				  register const __m512d theta = _mm512_mul_pd(vrx,v8_2pi);
				  _mm_prefetch((const char*)&ry[i+0],_MM_HINT_T0);
				  register const __m512d vry   = _mm512_load_pd(&ry[i]);
				  register const __m512d phi   = _mm512_mul_pd(vry,v8_2pi);
				  _mm_prefetch((const char*)&rz[i+0],_MM_HINT_T0);
				  register const __m512d vrz   = _mm512_load_pd(&rz[i]);
				  register const __m512d z     = _mm512_mul_pd(vrz,vrz);
				  _mm512_store_pd(&row9[i+0], _mm512_sub_pd(v8_1,z));
				  const register __m512d r     = _mm512_sqrt_pd(z);
				  const register __m512d vx    = _mm512_mul_pd(r,_mm512_sin_pd(phi));
				  const register __m512d vy    = _mm512_mul_pd(r,_mm512_cos_pd(phi));
				  const register __m512d vz    = _mm512_sqrt_pd(_mm512_sub_ps(v8_2,z));
				  _mm512_store_pd(row6[i+0], _mm512_mul_pd(vy,vz));
				  _mm512_store_pd(row3[i+0], _mm512_mul_pd(vx,vz));
				  const register __m512d st    = _mm512_sin_pd(theta);
			          const register __m512d ct    = _mm512_cos_pd(theta);
			          const register __m512d sx    = _mm512_fmsub_pd(vx,ct,_mm512_mul_pd(vy,st));
				  _mm512_store_pd(row7[i+0], _mm512_mul_pd(vz,sx));
				  _mm512_store_pd(row1[i+0], _mm512_fmsub_pd(vx,sx,ct));
				  const register __m512d sy    = _mm512_fmadd_pd(vx,st,_mm512_mul_pd(vy,ct));
				  _mm512_store_pd(&row8[i+0], _mm512_mul_pd(vz,sy));
				  _mm512_store_pd(&row2[i+0], _mm512_fmsub_pd(vx,sy,st));
				  _mm512_store_pd(&row4[i+0], _mm512_fmadd_pd(vy,sx,st));
				  _mm512_store_pd(&row5[i+0], _mm512_fmsub_pd(vy,sy,ct));
#endif
			      }

#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(8)
#pragma novector
#endif
                                for(; (i+0) < n; i += 1) {
                                         const double x = rx[i];
					const double theta = x*6.2831853071795864769253;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double y = ry[i];
					const double phi = y*6.2831853071795864769253;
					const double z = rz[i];
					const double zz = z+z;
					row9[i]        = 1.0-z;
					const double r  = sqrt(zz);
					const double vx = sin(phi) * r;
					const double vy = cos(phi) * r;
					const double vz = sqrt(2.0-zz);
				        row3[i]        = vx * vz;
					row6[i]        = vy * vz;
					const double st = sin(theta);
					const double ct = cos(theta);
					const double sx = vx * ct - vy * st;
					row1[i]        = vx * sx - ct;
					const double sy = vx * st + vy * ct;
					row8[i]        = vz * sy;
					row2[i]        = vx * sy - st;
					row4[i]        = vy * sx + st;
					row5[i]        = vy * sy - ct;
				}
		      }



		    

		      
     }


}











#endif /*__GMS_RAND_RMAT9X16_LOOPED_AVX512_HPP__*/
