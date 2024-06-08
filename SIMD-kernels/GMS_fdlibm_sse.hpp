
#ifndef __GMS_FDLIBM_SSE_HPP__
#define __GMS_FDLIBM_SSE_HPP__

/*
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 */
 
 /*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_info {

 const unsigned int GMS_FDLIBM_SSE_MAJOR = 1U;
 const unsigned int GMS_FDLIBM_SSE_MINOR = 0U;
 const unsigned int GMS_FDLIBM_SSE_MICRO = 0U;
 const unsigned int GMS_FDLIBM_SSE_FULLVER =
  1000U*GMS_FDLIBM_SSE_MAJOR+100U*GMS_FDLIBM_SSE_MINOR+10U*GMS_FDLIBM_SSE_MICRO;
 const char * const GMS_FDLIBM_SSE_CREATION_DATE  = "08-06-2024 08:22AM +00200 (SAT 08 JUN 2024 08:22AM GMT+2)";
 const char * const GMS_PDF_CDF_AVX_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const GMS_PDF_CDF_AVX_AUTHOR        = "Author: Sun Microsystems, Inc, modified by Bernard Gingold, contact: beniekg@gmail.com";
 const char * const GMS_PDF_CDF_AVX_SYNOPSIS      = "Manually vectorized (version: SSE) FDLIBM library content."


}

#include <immintrin.h>
#include <limits>
#include "GMS_config.h"
//#include "GMS_simd_utils.hpp"

#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x

#define RODATA_STORAGE_CONSTANTS 0

namespace gms {


         namespace math {
         
/* __ieee754_acos(x)
 * Method :                  
 *	acos(x)  = pi/2 - asin(x)
 *	acos(-x) = pi/2 + asin(x)
 * For |x|<=0.5
 *	acos(x) = pi/2 - (x + x*x^2*R(x^2))	(see asin.c)
 * For x>0.5
 * 	acos(x) = pi/2 - (pi/2 - 2asin(sqrt((1-x)/2)))
 *		= 2asin(sqrt((1-x)/2))  
 *		= 2s + 2s*z*R(z) 	...z=(1-x)/2, s=sqrt(z)
 *		= 2f + (2c + 2s*z*R(z))
 *     where f=hi part of s, and c = (z-f*f)/(s+f) is the correction term
 *     for f so that f+c ~ sqrt(z).
 * For x<-0.5
 *	acos(x) = pi - 2asin(sqrt((1-|x|)/2))
 *		= pi - 0.5*(s+s*z*R(z)), where z=(1-|x|)/2,s=sqrt(z)
 *
 * Special cases:
 *	if x is NaN, return x itself;
 *	if |x|>1, return NaN with invalid signal.
 *
 * Function needed: sqrt
 */
     
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d e_acos_xmm2r8(const __m128d x) {

#if (RODATA_STORAGE_CONSTANTS) == 1
                            static  __m128d one     = _mm_set1_pd(1.00000000000000000000e+00); //0x3FF00000, 0x00000000
                            static  __m128d pi      = _mm_set1_pd(3.14159265358979311600e+00); /* 0x400921FB, 0x54442D18 */
                            static  __m128d pio2_hi = _mm_set1_pd(1.57079632679489655800e+00); /* 0x3FF921FB, 0x54442D18 */
                            static  __m128d pio2_lo = _mm_set1_pd(6.12323399573676603587e-17); /* 0x3C91A626, 0x33145C07 */
                            static  __m128d pS0     = _mm_set1_pd(1.66666666666666657415e-01); /* 0x3FC55555, 0x55555555 */
                            static  __m128d pS1     = _mm_set1_pd(-3.25565818622400915405e-01); /* 0xBFD4D612, 0x03EB6F7D */
                            static  __m128d pS2     = _mm_set1_pd(2.01212532134862925881e-01); /* 0x3FC9C155, 0x0E884455 */
                            static  __m128d pS3     = _mm_set1_pd(-4.00555345006794114027e-02); /* 0xBFA48228, 0xB5688F3B */
                            static  __m128d pS4     = _mm_set1_pd(7.91534994289814532176e-04); /* 0x3F49EFE0, 0x7501B288 */
                            static  __m128d pS5     = _mm_set1_pd(3.47933107596021167570e-05); /* 0x3F023DE1, 0x0DFDF709 */
                            static  __m128d qS1     = _mm_set1_pd(-2.40339491173441421878e+00); /* 0xC0033A27, 0x1C8A2D4B */
                            static  __m128d qS2     = _mm_set1_pd(2.02094576023350569471e+00); /* 0x40002AE5, 0x9C598AC8 */
                            static  __m128d qS3     = _mm_set1_pd(-6.88283971605453293030e-01); /* 0xBFE6066C, 0x1B8D0159 */
                            static  __m128d qS4     = _mm_set1_pd(7.70381505559019352791e-02); /* 0x3FB3B8C5, 0xB12E9282 */
#else		              
		              __m128d one     = _mm_set1_pd(1.00000000000000000000e+00); //0x3FF00000, 0x00000000
                              __m128d pi      = _mm_set1_pd(3.14159265358979311600e+00); /* 0x400921FB, 0x54442D18 */
                              __m128d pio2_hi = _mm_set1_pd(1.57079632679489655800e+00); /* 0x3FF921FB, 0x54442D18 */
                              __m128d pio2_lo = _mm_set1_pd(6.12323399573676603587e-17); /* 0x3C91A626, 0x33145C07 */
                              __m128d pS0     = _mm_set1_pd(1.66666666666666657415e-01); /* 0x3FC55555, 0x55555555 */
                              __m128d pS1     = _mm_set1_pd(-3.25565818622400915405e-01); /* 0xBFD4D612, 0x03EB6F7D */
                              __m128d pS2     = _mm_set1_pd(2.01212532134862925881e-01); /* 0x3FC9C155, 0x0E884455 */
                              __m128d pS3     = _mm_set1_pd(-4.00555345006794114027e-02); /* 0xBFA48228, 0xB5688F3B */
                              __m128d pS4     = _mm_set1_pd(7.91534994289814532176e-04); /* 0x3F49EFE0, 0x7501B288 */
                              __m128d pS5     = _mm_set1_pd(3.47933107596021167570e-05); /* 0x3F023DE1, 0x0DFDF709 */
                              __m128d qS1     = _mm_set1_pd(-2.40339491173441421878e+00); /* 0xC0033A27, 0x1C8A2D4B */
                              __m128d qS2     = _mm_set1_pd(2.02094576023350569471e+00); /* 0x40002AE5, 0x9C598AC8 */
                              __m128d qS3     = _mm_set1_pd(-6.88283971605453293030e-01); /* 0xBFE6066C, 0x1B8D0159 */
                              __m128d qS4     = _mm_set1_pd(7.70381505559019352791e-02); /* 0x3FB3B8C5, 0xB12E9282 */
#endif     
                              __m128d zero    = _mm_setzero_pd();  
                              __m128d NAN     = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());   
                              __m128d z,p,q,r;
                              __m128d w,s,c,df;
                              __m128d ret;
                              __m128i vhx,vix,vlo;
                              double * __restrict__ dptr = nullptr;
                              int hx1,hx2;  
                              int lx1,lx2; 
                              __mmask8 m0; 
                              dptr = (double*)&x;
                              hx1  = __HI(dptr[0]);
                              hx2  = __HI(dptr[1]);
                              vhx  = _mm_setr_epi32(hx1,hx2,0,0);
                              vix  = _mm_and_si128(vhx,_mm_set1_epi32(0x7fffffff));
                              if(_mm_cmp_epi32_mask(vix,
                                              _mm_set1_epi32(0x3ff00000),_MM_CMPINT_NLE)) {
                                   lx1 = __LO(dptr[0]);
                                   lx2 = __LO(dptr[1]);
                                   vlo = _mm_setr_epi32(lx1,lx2,0,0);
                                   __m128i t0 = _mm_or_si128(_mm_sub_epi32(vix,_mm_set1_epi32(0x3ff00000)),vlo);
                                   if(_mm_cmp_epi32_mask(t0,
                                              _mm_setzero_si128(),_MM_CMPINT_EQ)){
                                          m0 = _mm_cmp_epi32_mask(vhx,_mm_setzero_si128(),_MM_CMPINT_NLT);
                                          ret= _mm_mask_blend_pd(m0,_mm_add_pd(pi,_mm_add_pd(pio2_lo,pio2_lo)),zero);
                                          return (ret);
                                   } 
                                   return (NAN);       
                              }
                              if(_mm_cmp_epi32_mask(vix,
                                              _mm_set1_epi32(0x3fe00000),_MM_CMPINT_LT)) {
                                      if(_mm_cmp_epi32_mask(vix,
                                                      _mm_set1_epi32(0x3c600000),_MM_CMPINT_LE)) 
                                                                                return (_mm_add_pd(pio2_hi,pio2_lo));
                                      z   = _mm_add_pd(x,x);
                                      p   = _mm_mul_pd(z,_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(pS5,z,pS4),z,pS3),z,pS2),z,pS1),z,pS0));
                                      q   = _mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(qS4,z,qS3),z,qS2),z,qS1),z,one);
                                      r   = _mm_div_pd(p,q);
                                      ret = _mm_sub_pd(pio2_hi,_mm_sub_pd(x,_mm_sub_pd(pio2_lo,_mm_mul_pd(x,r))));
                                      return (ret);
                              }
                              else if(_mm_cmp_epi32_mask(vhx,
                                                   _mm_setzero_si128(),_MM_CMPINT_LT)) {
                                      z   = _mm_mul_pd(_mm_add_pd(one,x),_mm_set1_pd(0.5));
                                      p   = _mm_mul_pd(z,_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(pS5,z,pS4),z,pS3),z,pS2),z,pS1),z,pS0));
                                      q   = _mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(qS4,z,qS3),z,qS2),z,qS1),z,one);
                                      s   = _mm_sqrt_pd(z);
                                      r   = _mm_div_pd(p,q);
                                      w   = _mm_fmsub_pd(r,s,pio2_lo);
                                      ret = _mm_sub_pd(pi,_mm_mul_pd(_mm_set1_pd(2.0),_mm_add_pd(s,w)));
                                      return (ret);                  
                              }
                              else {
                                      z   = _mm_mul_pd(_mm_sub_pd(one,x),_mm_set1_pd(0.5));
                                      s   = _mm_sqrt_pd(z);
                                      df  = s;
                                      double * __restrict ptr = (double*)&df;
                                      p   = _mm_mul_pd(z,_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(pS5,z,pS4),z,pS3),z,pS2),z,pS1),z,pS0));
                                      __LO(ptr[0]) = 0;
                                      __LO(ptr[1]) = 0;
                                      q   = _mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(_mm_fmadd_pd(qS4,z,qS3),z,qS2),z,qS1),z,one);
                                      r   = _mm_div_pd(p,q);
                                      df  = _mm_load_pd(&ptr[0]);
                                      c   = _mm_div_pd(_mm_sub_pd(z,_mm_mul_pd(df,df)),
                                                       _mm_add_pd(s,df));
                                      w   = _mm_fmadd_pd(r,s,c);
                                      ret = _mm_mul_pd(_mm_set1_pd(2.0),_mm_add_pd(df,w));
                                      return (ret);
                              }
		      }
     
        }// math


} // gms

























#endif /*__GMS_FDLIBM_SSE_HPP__*/
