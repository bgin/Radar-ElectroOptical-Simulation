
#ifndef __GMS_PDF_CDF_AVX512_HPP__
#define __GMS_PDF_CDF_AVX512_HPP__ 290520221332

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

 const unsigned int gGMS_PDF_CDF_AVX512_MAJOR = 1U;
 const unsigned int gGMS_PDF_CDF_AVX512_MINOR = 0U;
 const unsigned int gGMS_PDF_CDF_AVX512_MICRO = 0U;
 const unsigned int gGMS_PDF_CDF_AVX512_FULLVER =
  1000U*gGMS_PDF_CDF_AVX512_MAJOR+100U*gGMS_PDF_CDF_AVX512_MINOR+10U*gGMS_PDF_CDF_AVX512_MICRO;
 const char * const pgGMS_PDF_CDF_AVX512_CREATION_DATE = "29-05-2022 13:32 +00200 (SUN 29 MAY 2022 13:32 GMT+2)";
 const char * const pgGMS_PDF_CDF_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_PDF_CDF_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_PDF_CDF_AVX512_SYNOPSIS      = "Manually vectorized [AVX512] PDF,CDF functions"


}

#include <immintrin.h>
#include <limits>
#include "GMS_config.h"
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif
#include "GMS_simd_utils.hpp"


namespace gms {

        namespace math {

/*
!*****************************************************************************80
!
!! ANGLIT_CDF evaluates the Anglit CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
*/
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_cdf_zmm8r8(const __m512d x) {

                           const __m512d pi    = _mm512_set1_pd(3.14159265358979323846264338328);
			   const __m512d pi2   = _mm512_set1_pd(1.57079632679489661923132169164);
			   const __m512d _0_5  = _mm512_set1_pd(0.5);
			   const __m512d pi4   = _mm512_set1_pd(0.78539816339744830961566084582);
			   const __m512d _2    = _mm512_set1_pd(2.0);
			   const __m512d _1    = _mm512_set1_pd(1.0);
			   __m512d cdf,tmp;
			   __mmask8 m1,m2;
			   m1  = _mm512_cmp_pd_mask(x,pi4,_CMP_LT_OQ);
#if (USE_SLEEF_LIB) == 1
			   tmp = xcos(_mm512_fmadd_pd(_2,x,pi2));
#else
                           tmp = _mm512_cos_pd(_mm512_fmadd_pd(_2,x,pi2));
#endif
                           cdf = _mm512_mask_blend_pd(m1,_1,
			                          _mm512_sub_pd(_0_5,
						            _mm512_mul_pd(_0_5,tmp)));
			   m2  = _mm512_cmp_pd_mask(x,zmm8r8_negate(pi4),_CMP_LT_OQ);
			   cdf = _mm512_mask_blend_pd(m2,cdf,_mm512_setzero_pd());
			   return (cdf);
		  }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 anglit_cdf_zmm16r4(const __m512 x) {

                           const __m512 pi    = _mm512_set1_pd(3.14159265358979323846264338328f);
			   const __m512 pi2   = _mm512_set1_pd(1.57079632679489661923132169164f);
			   const __m512 _0_5  = _mm512_set1_pd(0.5f);
			   const __m512 pi4   = _mm512_set1_pd(0.78539816339744830961566084582f);
			   const __m512 _2    = _mm512_set1_pd(2.0f);
			   const __m512 _1    = _mm512_set1_pd(1.0f);
			   __m512 cdf,tmp;
			   __mmask16 m1,m2;
			   m1  = _mm512_cmp_ps_mask(x,pi4,_CMP_LT_OQ);
#if (USE_SLEEF_LIB) == 1
			   tmp = xcosf(_mm512_fmadd_ps(_2,x,pi2));
#else
                           tmp = _mm512_cos_ps(_mm512_fmadd_ps(_2,x,pi2));
#endif
                           cdf = _mm512_mask_blend_ps(m1,_1,
			                          _mm512_sub_ps(_0_5,
						            _mm512_mul_ps(_0_5,tmp)));
			   m2  = _mm512_cmp_ps_mask(x,zmm16r4_negate(pi4),_CMP_LT_OQ);
			   cdf = _mm512_mask_blend_ps(m2,cdf,_mm512_setzero_ps());
			   return (cdf);
		  }


/*
!*****************************************************************************80
!
!! ANGLIT_CDF_INV inverts the Anglit CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_cdf_inv_zmm8r8(const __m512d cdf) {

                           const __m512d _0    = _mm512_setzero_pd();
			   const __m512d _1    = _mm512_set1_pd(1.0);
			   const __m512d pi2   = _mm512_set1_pd(1.57079632679489661923132169164);
			   const __m512d _0_5  = _mm512_set1_pd(0.5);
			   const __m512d _2    = _mm512_set1_pd(2.0);
			   __m512d x,tmp;
			   __mmask8 m1,m2;
			   m1  = _mm512_cmp_pd_mask(cdf,_0,_CMP_LT_OQ);
			   m2  = _mm512_cmp_pd_mask(cdf,_1,_CMP_GT_OQ);
			   if(m1 || m2) {
                              x = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			      return (x);
			   }
#if (USE_SLEEF_LIB) == 1
                             tmp = xacos(_mm512_sub_pd(_1,
			                           _mm512_mul_pd(_2,cdf)));
#else
                             tmp = _mm512_acos_pd(_mm512_sub_pd(_1,
			                           _mm512_mul_pd(_2,cdf)));
#endif
                             x   = _mm512_fmsub_pd(_0_5,tmp,pi2);
			     return (x);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 anglit_cdf_inv_zmm16r4(const __m512 cdf) {

                           const __m512 _0    = _mm512_setzero_ps();
			   const __m512 _1    = _mm512_set1_ps(1.0f);
			   const __m512 pi2   = _mm512_set1_ps(1.57079632679489661923132169164f);
			   const __m512 _0_5  = _mm512_set1_ps(0.5f);
			   const __m512 _2    = _mm512_set1_ps(2.0f);
			   __m512 x,tmp;
			   __mmask16 m1,m2;
			   m1  = _mm512_cmp_ps_mask(cdf,_0,_CMP_LT_OQ);
			   m2  = _mm512_cmp_ps_mask(cdf,_1,_CMP_GT_OQ);
			   if(m1 || m2) {
                              x = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
			      return (x);
			   }
#if (USE_SLEEF_LIB) == 1
                             tmp = xacosf(_mm512_sub_ps(_1,
			                           _mm512_mul_ps(_2,cdf)));
#else
                             tmp = _mm512_acos_ps(_mm512_sub_ps(_1,
			                           _mm512_mul_ps(_2,cdf)));
#endif
                             x   = _mm512_fmsub_ps(_0_5,tmp,pi2);
			     return (x);
		    }

/*
 !*****************************************************************************80
!
!! ANGLIT_MEAN returns the mean of the Anglit PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_mean_zmm8r8() {

		             return (_mm512_setzero_pd());
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 anglit_mean_zmm16r4() {

		             return (_mm512_setzero_ps());
		     }

/*
!*****************************************************************************80
!
!! ANGLIT_PDF evaluates the Anglit PDF.
!
!  Discussion:
!
!    PDF(X) = sin ( 2 * X + PI / 2 ) for -PI/4 <= X <= PI/4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_pdf_zmm8r8(const __m512d x) {

                         const __m512d pi2  = _mm512_set1_pd(1.57079632679489661923132169164);
			 const __m512d _2   = _mm512_set1_pd(2.0);
			 const __m512d _0   = _mm512_setzero_pd();
			 __m512d pdf;
			 __mmask8 m,m1,m2;
			 m1  = _mm512_cmp_pd_mask(zmm8r8_negate(pi2),x,_CMP_LT_OQ);
			 m2  = _mm512_cmp_pd_mask(pi2,x,_CMP_LE_OQ);
			 m   = m1||m2;
#if (USE_SLEEF_LIB) == 1
			 tmp = xsin(_mm512_fmadd_pd(_2,x,pi2));
#else
                         tmp = _mm512_sin_pd(_mm512_fmadd_pd(_2,x,pi2));
#endif
			 pdf = _mm512_mask_blend_pd(m,tmp,_0);
			 return (pdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 anglit_pdf_zmm16r4(const __m512 x) {

                         const __m512 pi2  = _mm512_set1_ps(1.57079632679489661923132169164f);
			 const __m512 _2   = _mm512_set1_ps(2.0f);
			 const __m512 _0   = _mm512_setzero_ps();
			 __m512 pdf;
			 __mmask16 m,m1,m2;
			 m1  = _mm512_cmp_ps_mask(zmm16r4_negate(pi2),x,_CMP_LT_OQ);
			 m2  = _mm512_cmp_ps_mask(pi2,x,_CMP_LE_OQ);
			 m   = m1||m2;
#if (USE_SLEEF_LIB) == 1
			 tmp = xsinf(_mm512_fmadd_ps(_2,x,pi2));
#else
                         tmp = _mm512_sin_ps(_mm512_fmadd_ps(_2,x,pi2));
#endif
			 pdf = _mm512_mask_blend_ps(m,tmp,_0);
			 return (pdf);
		    }

/*
!*****************************************************************************80
!
!! ANGLIT_SAMPLE samples the Anglit PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
*/

#if defined(__ICC) || defined(__INTEL_COMPILER)
#include <svrng.h>
#else
#error 'Required Intel Compiler distribution'
#endif

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_sample_zmm8r8() {

                         __m512d cdf;
			 svrng_engine_t engine;
			 svrng_distribution_t uniform;
			 uint32_t seed    = 0U;
			 int32_t result   = -9999;
			 int32_t err      = -9999;
			 result           = _rdrand32_step(&seed);
			 if(!result) seed = 1563548129U;
			 engine           = svrng_new_mt19937_engine(seed);
			 err              = svrng_get_status();
			 if(err!=SVRNG_STATUS_OK) {
                            const __m512d nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_uniform_distribution_double(0.0,1.0);
			 const double * __restrict ptr = (const double*)(&svrng_generate8_double(engine,uniform));
			 cdf              = anglit_cdf_inv_zmm8r8(_mm512_loadu_pd(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_sample_zmm8r8(const __m512 cdf) {

                            return (anglit_cdf_inv_zmm8r8(cdf));
		      }
		    
		    

     }

}














#endif /*__GMS_PDF_CDF_AVX512_HPP__*/
