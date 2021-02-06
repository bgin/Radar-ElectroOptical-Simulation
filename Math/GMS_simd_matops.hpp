

#ifndef __GMS_SIMD_MATOPS_HPP__
#define __GMS_SIMD_MATOPS_HPP__



namespace file_info {

     const unsigned int gGMS_SIMD_MATOPS_MAJOR = 1U;
     const unsigned int gGMS_SIMD_MATOPS_MINOR = 0U;
     const unsigned int gGMS_SIMD_MATOPS_MICRO = 0U;
     const unsigned int gGMS_SIMD_MATOPS_FULLVER =
        1000U*gGMS_SIMD_MATOPS_MAJOR+100U*gGMS_SIMD_MATOPS_MINOR+
	10U*gGMS_SIMD_MATOPS_MICRO;
     const char * const pgGMS_SIMD_MATOPS_CREATE_DATE = "06-02-2021 2:00PM +00200 (SAT 02 FEB 2021 GMT+2)";
     const char * const pgGMS_SIMD_MATOPS_BUILD_DATE  = __DATE__ ":" __TIME__;
     const char * const pgGMS_SIMD_MATOPS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
}

#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

        namespace math {

	         __ATTR_ALWAYS_INLINE__
                 __ATTR_HOT__
		 __ATTR_ALIGN__(32)
		 static inline
                 void det_8mat4x4_sse_ps_u(float * __restrict  pSrcBuf, // size of arrays = 8*4*4 floats (8 4x4 matrices)
		                           float * __restrict  pDstBuf) {
                      __m128 zero = _mm_setzero_ps();
		      __m128 a11,a12,a13,a14;
		      __m128 a21,a22,a23,a24;
		      __m128 a31,a32,a33,a34;
		      __m128 a41,a42,a43,a44;
		      __m128 xmm0,xmm1,xmm2,xmm3,
		             xmm4,xmm5,xmm6,xmm7;
		      a11 = zero,a12 = zero,
		      a13 = zero,a14 = zero;
		      a21 = zero,a22 = zero,
		      a23 = zero,a24 = zero;
		      a31 = zero,a32 = zero,
		      a33 = zero,a34 = zero;
		      a41 = zero,a42 = zero,
		      a43 = zero,a44 = zero;
		      xmm0 = zero,xmm1 = zero,xmm2 = zero,
		      xmm3 = zero,xmm4 = zero,xmm5 = zero,
		      xmm6 = zero,xmm7 = zero;
                      float * __restrict  pSrc   = (float*)pSrcBuff;
		      float * __restrict  pDst   = (float*)pDstBuff;
		      constexpr int32_t n_iter   = 2;
		      constexpr int32_t n_stride = 32;
		      constexpr int32_t mat_size = 16;
		      constexpr int32_t n_mat    = 8;
		      for(int32_t i = 0;i != n_iter; ++i,pSrc += n_stride) {
                          xmm0 = _mm_loadu_ps( pSrc+0);
                          xmm1 = _mm_loadu_ps( pSrc+16);
			  xmm4 = _mm_unpacklo_ps(xmm0,xmm1);
			  xmm0 = _mm_unpackhi_ps(xmm0,xmm1);
			  xmm2 = _mm_loadu_ps( pSrc+32);
			  xmm3 = _mm_loadu_ps( pSrc+48);
			  xmm1 = _mm_unpacklo_ps(xmm2,xmm3);
			  xmm2 = _mm_unpackhi_ps(xmm2,xmm3);
                          xmm3 = _mm_shuffle_ps(xmm4,xmm1,0x4e);
			  a11  = _mm_blend_ps(xmm4,xmm3,0xc);
			  a12  = _mm_blend_ps(xmm1,xmm3,0x3);
			  xmm3 = _mm_shuffle_ps(xmm0,xmm2,0x4e);
			  a13  = _mm_blend_ps(xmm0,xmm3,0xc);
			  a14  = _mm_blend_ps(xmm2,xmm3,0x3);

			  xmm0 = _mm_loadu_ps(pSrc+4);
			  xmm1 = _mm_loadu_ps(pSrc+20);
			  xmm4 = _mm_unpacklo_ps(xmm0,xmm1);
			  xmm0 = _mm_unpackhi_ps(xmm0,xmm1);
			  xmm2 = _mm_loadu_ps(pSrc+36);
			  xmm3 = _mm_loadu_ps(pSrc+52);
			  xmm1 = _mm_unpacklo_ps(xmm2,xmm3);
			  xmm2 = _mm_unpackhi_ps(xmm2,xmm3);
			  xmm3 = _mm_shuffle_ps(xmm4,xmm1,0x4e);
			  a21  = _mm_blend_ps(xmm4,xmm3,0xc);
			  a22  = _mm_blend_ps(xmm1,xmm3,0x3);
			  xmm3 = _mm_shuffle_ps(xmm0,xmm2,0x4e);
			  a23  = _mm_blend_ps(xmm0,xmm3,0xc);
			  a24  = _mm_blend_ps(xmm2,xmm3,0x3);
			  
		      }
		 }

    } // math


} // gms






#endif /*__GMS_SIMD_MATOPS_HPP__*/
