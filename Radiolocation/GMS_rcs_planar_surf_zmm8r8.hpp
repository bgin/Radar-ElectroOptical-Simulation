

#ifndef __GMS_RCS_PLANAR_SURF_ZMM8R8_HPP__
#define __GMS_RCS_PLANAR_SURF_ZMM8R8_HPP__ 180420231543



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

namespace file_version {

    const unsigned int GMS_RCS_PLANAR_SURF_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM8R8_FULLVER =
      1000U*GMS_RCS_PLANAR_SURF_ZMM8R8_MAJOR+
      100U*GMS_RCS_PLANAR_SURF_ZMM8R8_MINOR+
      10U*GMS_RCS_PLANAR_SURF_ZMM8R8_MICRO;
    const char * const GMS_RCS_PLANAR_SURF_ZMM8R8_CREATION_DATE = "18-04-2023 15:43 PM +00200 (TUE 18 APR 2023 GMT+2)";
    const char * const GMS_RCS_PLANAR_SURF_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_PLANAR_SURF_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_PLANAR_SURF_ZMM8R8_DESCRIPTION   = "AVX512 optimized Planar Surfaces Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimddp.hpp"
#include "GMS_complex_zmm8r8.hpp"
#include "GMS_simd_utils.hpp"


namespace  gms {


         namespace radiolocation {


                 /*
                       Complex impedances.
                       Formula 7.1-6
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm8r8(const __m512d tht,
                                        const __m512d mur,
                                        const __m512d mui,
                                        const __m512d epsr,
                                        const __m512d epsi,
                                        __m512d * __restrict zr,
                                        __m512d * __restrict zi) {

                        register __m512d cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcos(tht);
                        wrkc = _mm512_setzero_pd();
                        cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_pd(cost);
                        csqrt_zmm8r8(divr,divi,&wrkc,&csqr,&csqi);
                        *zr = _mm512_mul_pd(invc,csqr);
                        *zi = _mm512_mul_pd(invc,csqi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                        const double * __restrict __ATTR_ALIGN__(64) pmur,
                                        const double * __restrict __ATTR_ALIGN__(64) pmui,
                                        const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                        const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                        double * __restrict __ATTR_ALIGN__(64) zr,
                                        double * __restrict __ATTR_ALIGN__(64)  zi) {

                        register __m512d tht  = _mm512_load_pd(&ptht[0]);
                        register __m512d mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d mui  = _mm512_load_pd(&pmui[0]);
                        register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcos(tht);
                        wrkc = _mm512_setzero_pd();
                        cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_pd(cost);
                        csqrt_zmm8r8(divr,divi,&wrkc,&csqr,&csqi);
                        _mm512_store_pd(&zr[0] ,_mm512_mul_pd(invc,csqr));
                        _mm512_store_pd(&zi[0] ,_mm512_mul_pd(invc,csqi));
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm8r8_u(const double * __restrict  ptht,
                                        const double * __restrict  pmur,
                                        const double * __restrict  pmui,
                                        const double * __restrict  pepsr,
                                        const double * __restrict  pepsi,
                                        double * __restrict  zr,
                                        double * __restrict   zi) {

                        register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                        register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                        register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcos(tht);
                        wrkc = _mm512_setzero_pd();
                        cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_pd(cost);
                        csqrt_zmm8r8(divr,divi,&wrkc,&csqr,&csqi);
                        _mm512_storeu_pd(&zr[0] ,_mm512_mul_pd(invc,csqr));
                        _mm512_storeu_pd(&zi[0] ,_mm512_mul_pd(invc,csqi));
                }


                 /*
                          Equivalent complex impedances.
                          Formula 7.1-4
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f714_zmm8r8( const __m512d tht1,
                                        const __m512d mur1,
                                        const __m512d mui1,
                                        const __m512d epsr1,
                                        const __m512d epsi1,
                                        const __m512d tht2,
                                        const __m512d mur2,
                                        const __m512d mui2,
                                        const __m512d epsr2,
                                        const __m512d epsi2,
                                        __m512d * __restrict Rr,
                                        __m512d * __restrict Ri) {

                     using namespace gms::math;
                     register __m512d z1r,z1i,z2r,z2i;
                     register __m512d t0r,t0i,t1r,t1i;
                     zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     zi_f716_zmm8r8(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     t0r = _mm512_sub_pd(z1r,z2r);
                     t1r = _mm512_add_pd(z1r,z2r);
                     t0i = _mm512_sub_pd(z1i,z2i);
                     t1i = _mm512_add_pd(z1i,z2i);
                     t0r = negate_zmm8r8(t0r);
                     t0i = negate_zmm8r8(t0i);
                     cdiv_zmm8r8(t0r,t0i,t1r,t1i,*Rr,*Ri);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f714_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                        const double * __restrict __ATTR_ALIGN__(64) pmur1,
                                        const double * __restrict __ATTR_ALIGN__(64) pmui1,
                                        const double * __restrict __ATTR_ALIGN__(64) pepsr1,
                                        const double * __restrict __ATTR_ALIGN__(64) pepsi1,
                                        const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                        const double * __restrict __ATTR_ALIGN__(64) pmur2,
                                        const double * __restrict __ATTR_ALIGN__(64) pmui2,
                                        const double * __restrict __ATTR_ALIGN__(64) pepsr2,
                                        const double * __restrict __ATTR_ALIGN__(64) pepsi2,
                                        double * __restrict __ATTR_ALIGN__(64) Rr,
                                        double * __restrict __ATTR_ALIGN__(64) Ri) {

                     using namespace gms::math;
                     register __m512d tht1  = _mm512_load_pd(&ptht1[0]);
                     register __m512d mur1  = _mm512_load_pd(&pmur1[0]);
                     register __m512d mui1  = _mm512_load_pd(&pmui1[0]);
                     register __m512d epsr1 = _mm512_load_pd(&pepsr1[0]);
                     register __m512d epsi1 = _mm512_load_pd(&pepsi1[0]);
                     register __m512d tht2  = _mm512_load_pd(&ptht2[0]);
                     register __m512d mur2  = _mm512_load_pd(&pmur2[0]);
                     register __m512d mui2  = _mm512_load_pd(&pmui2[0]);
                     register __m512d epsr2 = _mm512_load_pd(&pepsr2[0]);
                     register __m512d epsi2 = _mm512_load_pd(&pepsi2[0]);
                     register __m512d z1r,z1i,z2r,z2i;
                     register __m512d t0r,t0i,t1r,t1i,resr,resi;
                     zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     zi_f716_zmm8r8(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     t0r = _mm512_sub_pd(z1r,z2r);
                     t1r = _mm512_add_pd(z1r,z2r);
                     t0i = _mm512_sub_pd(z1i,z2i);
                     t1i = _mm512_add_pd(z1i,z2i);
                     t0r = negate_zmm8r8(t0r);
                     t0i = negate_zmm8r8(t0i);
                     cdiv_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                     _mm512_store_pd(&Rr[0], resr);
                     _mm512_store_pd(&Ri[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f714_zmm8r8_u( const double * __restrict  ptht1,
                                        const double * __restrict  pmur1,
                                        const double * __restrict  pmui1,
                                        const double * __restrict  pepsr1,
                                        const double * __restrict  pepsi1,
                                        const double * __restrict  ptht2,
                                        const double * __restrict  pmur2,
                                        const double * __restrict  pmui2,
                                        const double * __restrict  pepsr2,
                                        const double * __restrict  pepsi2,
                                        double * __restrict  Rr,
                                        double * __restrict  Ri) {

                     using namespace gms::math;
                     register __m512d tht1  = _mm512_loadu_pd(&ptht1[0]);
                     register __m512d mur1  = _mm512_loadu_pd(&pmur1[0]);
                     register __m512d mui1  = _mm512_loadu_pd(&pmui1[0]);
                     register __m512d epsr1 = _mm512_loadu_pd(&pepsr1[0]);
                     register __m512d epsi1 = _mm512_loadu_pd(&pepsi1[0]);
                     register __m512d tht2  = _mm512_loadu_pd(&ptht2[0]);
                     register __m512d mur2  = _mm512_loadu_pd(&pmur2[0]);
                     register __m512d mui2  = _mm512_loadu_pd(&pmui2[0]);
                     register __m512d epsr2 = _mm512_loadu_pd(&pepsr2[0]);
                     register __m512d epsi2 = _mm512_loadu_pd(&pepsi2[0]);
                     register __m512d z1r,z1i,z2r,z2i;
                     register __m512d t0r,t0i,t1r,t1i,resr,resi;
                     zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     zi_f716_zmm8r8(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     t0r = _mm512_sub_pd(z1r,z2r);
                     t1r = _mm512_add_pd(z1r,z2r);
                     t0i = _mm512_sub_pd(z1i,z2i);
                     t1i = _mm512_add_pd(z1i,z2i);
                     t0r = negate_zmm8r8(t0r);
                     t0i = negate_zmm8r8(t0i);
                     cdiv_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                     _mm512_storeu_pd(&Rr[0], resr);
                     _mm512_storeu_pd(&Ri[0], resi);
                }


                  /*
                        Transmission coefficient components.
                        Formula 7.1-5
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void T_f715_zmm8r8( const __m512d tht1,
                                        const __m512d mur1,
                                        const __m512d mui1,
                                        const __m512d epsr1,
                                        const __m512d epsi1,
                                        const __m512d tht2,
                                        const __m512d mur2,
                                        const __m512d mui2,
                                        const __m512d epsr2,
                                        const __m512d epsi2,
                                        __m512d * __restrict Tr,
                                        __m512d * __restrict Ti) {

                     const __m512d _2 = _mm512_set1_pd(2.0f);
                     register __m512d z1r,z1i,z2r,z2i;
                     register __m512d t0r,t0i,t1r,t1i;
                     zi_f716_zmm8r8(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     t0r = _mm512_mul_pd(_2,z2r);
                     t1r = _mm512_add_pd(z1r,z2r);
                     t0i = _mm512_mul_pd(_2,z2i);
                     t1i = _mm512_add_pd(z1i,z2i);
                     cdiv_zmm8r8(t0r,t0i,t1r,t1i,*Tr,*Ti);
             }


                /*
                        Reflection coefficient special cases:
                        1) k1<k2, eps1,eps2 (real), mu1 = m2 = mu0
                        Formula 7.1-17
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7117_zmm8r8(const __m512d tht,
                                          const __m512d eps1,
                                          const __m512d eps2) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d e1e2,sqr1,sqr2,num,den,R;
                          register __m512d cost,sint,x0,x1;
                          e1e2 = _mm512_div_pd(eps1,eps2);
                          cost = xcos(tht);
                          sqr1 = _mm512_sqrt_pd(e1e2);
                          sint = xsin(tht);
                          x0   = _mm512_mul_pd(sqr1,cost);
                          x1   = _mm512_mul_pd(_mm512_sub_pd(_1,e1e2),
                                                     _mm512_mul_pd(sint,sint));
                          sqr2 = _mm512_sqrt_pd(x1);
                          num  = _mm512_sub_pd(x0,x1);
                          den  = _mm512_add_pd(x0,x1);
                          R    = _mm512_div_pd(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7117_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                            const double * __restrict __ATTR_ALIGN__(64) peps1,
                                            const double * __restrict __ATTR_ALIGN__(64) peps2) {

                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps2= _mm512_load_pd(&peps2[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d e1e2,sqr1,sqr2,num,den,R;
                          register __m512d cost,sint,x0,x1;
                          e1e2 = _mm512_div_pd(eps1,eps2);
                          cost = xcos(tht);
                          sqr1 = _mm512_sqrt_pd(e1e2);
                          sint = xsin(tht);
                          x0   = _mm512_mul_pd(sqr1,cost);
                          x1   = _mm512_mul_pd(_mm512_sub_pd(_1,e1e2),
                                                     _mm512_mul_pd(sint,sint));
                          sqr2 = _mm512_sqrt_pd(x1);
                          num  = _mm512_sub_pd(x0,x1);
                          den  = _mm512_add_pd(x0,x1);
                          R    = _mm512_div_pd(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7117_zmm8r8_u(const double * __restrict  ptht,
                                            const double * __restrict  peps1,
                                            const double * __restrict  peps2) {

                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps2= _mm512_loadu_pd(&peps2[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d e1e2,sqr1,sqr2,num,den,R;
                          register __m512d cost,sint,x0,x1;
                          e1e2 = _mm512_div_pd(eps1,eps2);
                          cost = xcos(tht);
                          sqr1 = _mm512_sqrt_pd(e1e2);
                          sint = xsin(tht);
                          x0   = _mm512_mul_pd(sqr1,cost);
                          x1   = _mm512_mul_pd(_mm512_sub_pd(_1,e1e2),
                                                     _mm512_mul_pd(sint,sint));
                          sqr2 = _mm512_sqrt_pd(x1);
                          num  = _mm512_sub_pd(x0,x1);
                          den  = _mm512_add_pd(x0,x1);
                          R    = _mm512_div_pd(num,den);
                          return (R);
                }


                   /*
                        Reflection coefficient special cases:
                        1) k1<k2, eps1,eps2 (real), mu1 = m2 = mu0
                        Formula 7.1-18
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7118_zmm8r8(const __m512d tht,
                                          const __m512d eps1,
                                          const __m512d eps2) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d e1e2,sqr1,sqr2,num,den,R;
                          register __m512d cost,sint,x0,x1;
                          e1e2 = _mm512_div_pd(eps2,eps1);
                          cost = xcos(tht);
                          sqr1 = _mm512_sqrt_pd(e1e2);
                          sint = xsin(tht);
                          x0   = _mm512_mul_pd(sqr1,cost);
                          x1   = _mm512_mul_pd(_mm512_sub_pd(_1,e1e2),
                                                     _mm512_mul_pd(sint,sint));
                          sqr2 = _mm512_sqrt_pd(x1);
                          num  = _mm512_sub_pd(x0,x1);
                          den  = _mm512_add_pd(x0,x1);
                          R    = _mm512_div_pd(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7118_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                            const double * __restrict __ATTR_ALIGN__(64) peps1,
                                            const double * __restrict __ATTR_ALIGN__(64) peps2) {

                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps2= _mm512_load_pd(&peps2[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d e1e2,sqr1,sqr2,num,den,R;
                          register __m512d cost,sint,x0,x1;
                          e1e2 = _mm512_div_pd(eps2,eps1);
                          cost = xcos(tht);
                          sqr1 = _mm512_sqrt_pd(e1e2);
                          sint = xsin(tht);
                          x0   = _mm512_mul_pd(sqr1,cost);
                          x1   = _mm512_mul_pd(_mm512_sub_pd(_1,e1e2),
                                                     _mm512_mul_pd(sint,sint));
                          sqr2 = _mm512_sqrt_pd(x1);
                          num  = _mm512_sub_pd(x0,x1);
                          den  = _mm512_add_pd(x0,x1);
                          R    = _mm512_div_pd(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7118_zmm8r8_u(const double * __restrict  ptht,
                                            const double * __restrict  peps1,
                                            const double * __restrict  peps2) {

                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps2= _mm512_loadu_pd(&peps2[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d e1e2,sqr1,sqr2,num,den,R;
                          register __m512d cost,sint,x0,x1;
                          e1e2 = _mm512_div_pd(eps2,eps1);
                          cost = xcos(tht);
                          sqr1 = _mm512_sqrt_pd(e1e2);
                          sint = xsin(tht);
                          x0   = _mm512_mul_pd(sqr1,cost);
                          x1   = _mm512_mul_pd(_mm512_sub_pd(_1,e1e2),
                                                     _mm512_mul_pd(sint,sint));
                          sqr2 = _mm512_sqrt_pd(x1);
                          num  = _mm512_sub_pd(x0,x1);
                          den  = _mm512_add_pd(x0,x1);
                          R    = _mm512_div_pd(num,den);
                          return (R);
                }



                    /*
                        Reflection coefficient special cases:
                        2) k2<k1, eps1,eps2 (real), mu1 = mu2 = mu0
                        Formula 7.1-23
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7123_zmm8r8(const __m512d tht,
                                        const __m512d eps2,
                                        const __m512d eps1,
                                        __m512d * __restrict Rr,
                                        __m512d * __restrict Ri) {

                        const __m512d n2 = _mm512_set1_pd(-2.0f);
                        register __m512d ear,eai,cer,cei;
                        register __m512d sint,cost,e2e1,rat,arg,atarg;
                        register __m512d x0,x1;
                        e2e1 = _mm512_div_pd(eps2,eps1);
                        cost = xcos(tht);
                        ear  = _mm512_setzero_pd();
                        sint = xsin(tht);
                        x0   = _mm512_mul_pd(cost,cost);
                        x1   = _mm512_sub_pd(_mm512_mul_pd(sint,sint),e2e1);
                        rat  = _mm512_div_pd(x1,x0);
                        arg  = _mm512_sqrt_pd(rat);
                        atarg= xatan(arg);
                        eai  = _mm512_mul_pd(n2,atarg);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        *Rr = cer;
                        *Ri = cei;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7123_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) ptht,
                                            const double * __restrict __ATTR_ALIGN__(64) peps1,
                                            const double * __restrict __ATTR_ALIGN__(64) peps2
                                            double * __restrict __ATTR_ALIGN__(64) Rr,
                                            double * __restrict __ATTR_ALIGN__(64) Ri) {

                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d eps1= _mm512_load_pd(&peps1[0]);
                        register __m512d eps2= _mm512_load_pd(&peps2[0]);
                        const __m512d n2 = _mm512_set1_pd(-2.0f);
                        register __m512d ear,eai,cer,cei;
                        register __m512d sint,cost,e2e1,rat,arg,atarg;
                        register __m512d x0,x1;
                        e2e1 = _mm512_div_pd(eps2,eps1);
                        cost = xcos(tht);
                        ear  = _mm512_setzero_pd();
                        sint = xsin(tht);
                        x0   = _mm512_mul_pd(cost,cost);
                        x1   = _mm512_sub_pd(_mm512_mul_pd(sint,sint),e2e1);
                        rat  = _mm512_div_pd(x1,x0);
                        arg  = _mm512_sqrt_pd(rat);
                        atarg= xatan(arg);
                        eai  = _mm512_mul_pd(n2,atarg);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        _mm512_store_pd(&Rr[0] ,cer);
                        _mm512_store_pd(&Ri[0] ,cei);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7123_zmm8r8_u(  const double * __restrict  ptht,
                                            const double * __restrict  peps1,
                                            const double * __restrict  peps2
                                            double * __restrict  Rr,
                                            double * __restrict  Ri) {

                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                        register __m512d eps2= _mm512_loadu_pd(&peps2[0]);
                        const __m512d n2 = _mm512_set1_pd(-2.0f);
                        register __m512d ear,eai,cer,cei;
                        register __m512d sint,cost,e2e1,rat,arg,atarg;
                        register __m512d x0,x1;
                        e2e1 = _mm512_div_pd(eps2,eps1);
                        cost = xcos(tht);
                        ear  = _mm512_setzero_pd();
                        sint = xsin(tht);
                        x0   = _mm512_mul_pd(cost,cost);
                        x1   = _mm512_sub_pd(_mm512_mul_pd(sint,sint),e2e1);
                        rat  = _mm512_div_pd(x1,x0);
                        arg  = _mm512_sqrt_pd(rat);
                        atarg= xatan(arg);
                        eai  = _mm512_mul_pd(n2,atarg);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        _mm512_storeu_pd(&Rr[0] ,cer);
                        _mm512_storeu_pd(&Ri[0] ,cei);
                }


                    /*
                        Reflection coefficient special cases:
                        2) k2<k1, eps1,eps2 (real), mu1 = mu2 = mu0
                        Formula 7.1-24
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7124_zmm8r8(const __m512d tht,
                                        const __m512d eps2,
                                        const __m512d eps1,
                                        __m512d * __restrict Rr,
                                        __m512d * __restrict Ri) {

                        const __m512d n2 = _mm512_set1_pd(-2.0f);
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei;
                        register __m512d sint,cost,e2e1,e1e2,rat,arg,atarg;
                        register __m512d x0,x1;
                        e2e1 = _mm512_div_pd(eps2,eps1);
                        sint = xsin(tht);
                        ear  = _mm512_setzero_pd();
                        cost = xcos(tht);
                        x0   = _mm512_mul_pd(e2e1,_mm512_fmsub_pd(sint,sint,_1));
                        e1e2 = _mm512_div_pd(eps1,eps2);
                        x1   = _mm512_mul_pd(e1e2,_mm512_mul_pd(cost,cost));
                        rat  = _mm512_div_pd(x0,x1);
                        arg  = _mm512_sqrt_pd(rat);
                        atarg= xatan(arg);
                        eai  = _mm512_mul_pd(n2,atarg);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        *Rr  = cer;
                        *Ri  = cei;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7124_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) ptht,
                                            const double * __restrict __ATTR_ALIGN__(64) peps1,
                                            const double * __restrict __ATTR_ALIGN__(64) peps2
                                            double * __restrict __ATTR_ALIGN__(64) Rr,
                                            double * __restrict __ATTR_ALIGN__(64) Ri) {

                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d eps1= _mm512_load_pd(&peps1[0]);
                        register __m512d eps2= _mm512_load_pd(&peps2[0]);
                        const __m512d n2 = _mm512_set1_pd(-2.0f);
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei;
                        register __m512d sint,cost,e2e1,e1e2,rat,arg,atarg;
                        register __m512d x0,x1;
                        e2e1 = _mm512_div_pd(eps2,eps1);
                        sint = xsin(tht);
                        ear  = _mm512_setzero_pd();
                        cost = xcos(tht);
                        x0   = _mm512_mul_pd(e2e1,_mm512_fmsub_pd(sint,sint,_1));
                        e1e2 = _mm512_div_pd(eps1,eps2);
                        x1   = _mm512_mul_pd(e1e2,_mm512_mul_pd(cost,cost));
                        rat  = _mm512_div_pd(x0,x1);
                        arg  = _mm512_sqrt_pd(rat);
                        atarg= xatan(arg);
                        eai  = _mm512_mul_pd(n2,atarg);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        _mm512_store_pd(&Rr[0] ,cer);
                        _mm512_store_pd(&Ri[0] ,cei);
                }



                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7124_zmm8r8_u(  const double * __restrict  ptht,
                                            const double * __restrict  peps1,
                                            const double * __restrict  peps2
                                            double * __restrict  Rr,
                                            double * __restrict  Ri) {

                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                        register __m512d eps2= _mm512_loadu_pd(&peps2[0]);
                        const __m512d n2 = _mm512_set1_pd(-2.0f);
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei;
                        register __m512d sint,cost,e2e1,e1e2,rat,arg,atarg;
                        register __m512d x0,x1;
                        e2e1 = _mm512_div_pd(eps2,eps1);
                        sint = xsin(tht);
                        ear  = _mm512_setzero_pd();
                        cost = xcos(tht);
                        x0   = _mm512_mul_pd(e2e1,_mm512_fmsub_pd(sint,sint,_1));
                        e1e2 = _mm512_div_pd(eps1,eps2);
                        x1   = _mm512_mul_pd(e1e2,_mm512_mul_pd(cost,cost));
                        rat  = _mm512_div_pd(x0,x1);
                        arg  = _mm512_sqrt_pd(rat);
                        atarg= xatan(arg);
                        eai  = _mm512_mul_pd(n2,atarg);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        _mm512_storeu_pd(&Rr[0] ,cer);
                        _mm512_storeu_pd(&Ri[0] ,cei);
                }


                  /*
                       Lateral displacement of the incident ray.
                       Formula 7.1-27
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d D_f7127_zmm8r8(    const __m512d gam0,
                                              const __m512d tht,
                                              const __m512d eps2,
                                              const __m512d eps1) {

                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
                         register __m512d g0pi,ttht,sint,e2e1,sqr,rat;
                         register __m512d D;
                         g0pi = _mm512_mul_pd(gam0,invpi);
                         e2e1 = _mm512_div_pd(eps2,eps1);
                         ttht = xtan(tht);
                         sint = _mm512_fmsub_pd(sint,sint,e2e1);
                         sqr  = _mm512_sqrt_pd(sint);
                         rat  = _mm512_div_pd(ttht,sqr);
                         D    = _mm512_mul_pd(g0pi,rat);
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d D_f7127_zmm8r8_a(    const double * __restrict __ATTR_ALIGN__(64)  pgam0,
                                                const double * __restrict __ATTR_ALIGN__(64)  ptht,
                                                const double * __restrict __ATTR_ALIGN__(64)  peps2,
                                                const double * __restrict __ATTR_ALIGN__(64)  peps1) {

                         register __m512d gam0= _mm512_load_pd(&pgam0[0]);
                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         register __m512d eps1= _mm512_load_pd(&peps1[0]);
                         register __m512d eps2= _mm512_load_pd(&peps2[0]);
                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
                         register __m512d g0pi,ttht,sint,e2e1,sqr,rat;
                         register __m512d D;
                         g0pi = _mm512_mul_pd(gam0,invpi);
                         e2e1 = _mm512_div_pd(eps2,eps1);
                         ttht = xtan(tht);
                         sint = _mm512_fmsub_pd(sint,sint,e2e1);
                         sqr  = _mm512_sqrt_pd(sint);
                         rat  = _mm512_div_pd(ttht,sqr);
                         D    = _mm512_mul_pd(g0pi,rat);
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d D_f7127_zmm8r8_u(    const double * __restrict   pgam0,
                                                const double * __restrict   ptht,
                                                const double * __restrict   peps2,
                                                const double * __restrict   peps1) {

                         register __m512d gam0= _mm512_loadu_pd(&pgam0[0]);
                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                         register __m512d eps2= _mm512_loadu_pd(&peps2[0]);
                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
                         register __m512d g0pi,ttht,sint,e2e1,sqr,rat;
                         register __m512d D;
                         g0pi = _mm512_mul_pd(gam0,invpi);
                         e2e1 = _mm512_div_pd(eps2,eps1);
                         ttht = xtan(tht);
                         sint = _mm512_fmsub_pd(sint,sint,e2e1);
                         sqr  = _mm512_sqrt_pd(sint);
                         rat  = _mm512_div_pd(ttht,sqr);
                         D    = _mm512_mul_pd(g0pi,rat);
                         return (D);
                }


                    /*
                       Lateral displacement of the incident ray.
                       Formula 7.1-28
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d D_f7128_zmm8r8(    const __m512d gam0,
                                              const __m512d tht,
                                              const __m512d eps2,
                                              const __m512d eps1) {

                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
                         register __m512d g0pi,ttht,sint,e2e1,e1e2,sqr,rat;
                         register __m512d D;
                         g0pi = _mm512_mul_pd(gam0,invpi);
                         e2e1 = _mm512_div_pd(eps2,eps1);
                         ttht = xtan(tht);
                         sint = _mm512_fmsub_pd(sint,sint,e2e1);
                         e1e2 = _mm512_div_pd(eps1,eps2);
                         sqr  = _mm512_sqrt_pd(sint);
                         rat  = _mm512_div_pd(ttht,sqr);
                         D    = _mm512_mul_pd(g0pi,_mm512_mul_pd(e1e2,rat));
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d D_f7128_zmm8r8_a(    const double * __restrict __ATTR_ALIGN__(64)  pgam0,
                                                const double * __restrict __ATTR_ALIGN__(64)  ptht,
                                                const double * __restrict __ATTR_ALIGN__(64)  peps2,
                                                const double * __restrict __ATTR_ALIGN__(64)  peps1) {

                         register __m512d gam0= _mm512_load_pd(&pgam0[0]);
                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         register __m512d eps1= _mm512_load_pd(&peps1[0]);
                         register __m512d eps2= _mm512_load_pd(&peps2[0]);
                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
                         register __m512d g0pi,ttht,sint,e2e1,e1e2,sqr,rat;
                         register __m512d D;
                         g0pi = _mm512_mul_pd(gam0,invpi);
                         e2e1 = _mm512_div_pd(eps2,eps1);
                         ttht = xtan(tht);
                         sint = _mm512_fmsub_pd(sint,sint,e2e1);
                         e1e2 = _mm512_div_pd(eps1,eps2);
                         sqr  = _mm512_sqrt_pd(sint);
                         rat  = _mm512_div_pd(ttht,sqr);
                         D    = _mm512_mul_pd(g0pi,_mm512_mul_pd(e1e2,rat));
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d D_f7128_zmm8r8_u(    const double * __restrict   pgam0,
                                                const double * __restrict   ptht,
                                                const double * __restrict   peps2,
                                                const double * __restrict   peps1) {

                         register __m512d gam0= _mm512_loadu_pd(&pgam0[0]);
                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                         register __m512d eps2= _mm512_loadu_pd(&peps2[0]);
                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
                         register __m512d g0pi,ttht,sint,e2e1,e1e2,sqr,rat;
                         register __m512d D;
                         g0pi = _mm512_mul_pd(gam0,invpi);
                         e2e1 = _mm512_div_pd(eps2,eps1);
                         ttht = xtan(tht);
                         sint = _mm512_fmsub_pd(sint,sint,e2e1);
                         e1e2 = _mm512_div_pd(eps1,eps2);
                         sqr  = _mm512_sqrt_pd(sint);
                         rat  = _mm512_div_pd(ttht,sqr);
                         D    = _mm512_mul_pd(g0pi,_mm512_mul_pd(e1e2,rat));
                         return (D);
                }


                      /*
                             For (k1/k2)^2*sin^2(theta)<<1 (Simplification
                             of formulae 7.1-9 and 7.1-10).
                             Formula 7.1-29
                        */


                    
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7129_zmm8r8(const __m512d tht,
                                        const __m512d mur1,
                                        const __m512d mui1,
                                        const __m512d epsr1,
                                        const __m512d epsi1,
                                        const __m512d mur2,
                                        const __m512d mui2,
                                        const __m512d epsr2,
                                        const __m512d epsi2,
                                        __m512d * __restrict Rr,
                                        __m512d * __restrict Ri) {

                       register __m512d z1r,z1i,z2r,z2i;
                       register __m512d cost,numr,numi,denr,deni;
                       register __m512d resr,resi;
                       zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcos(tht);
                       zi_f716_zmm8r8(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_fmsub_pd(z2r,cost,z1r);
                       denr = _mm512_fmadd_pd(z2r,cost,z2r);
                       numi = _mm512_fmsub_pd(z2i,cost,z1i);
                       deni = _mm512_fmadd_pd(z2i,cost,z2i);
                       cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7129_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur1,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui1,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr1,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi1,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur2,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui2,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr2,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi2,
                                          double * __restrict __ATTR_ALIGN__(64) Rr,
                                          double * __restrict __ATTR_ALIGN__(64) Ri) {

                       register __m512d tht  = _mm512_load_pd(&ptht[0]);
                       register __m512d mur1  = _mm512_load_pd(&pmur1[0]);
                       register __m512d mui1  = _mm512_load_pd(&pmui1[0]);
                       register __m512d epsr1 = _mm512_load_pd(&pepsr1[0]);
                       register __m512d epsi1 = _mm512_load_pd(&pepsi1[0]);
                       register __m512d mur2  = _mm512_load_pd(&pmur2[0]);
                       register __m512d mui2  = _mm512_load_pd(&pmui2[0]);
                       register __m512d epsr2 = _mm512_load_pd(&pepsr2[0]);
                       register __m512d epsi2 = _mm512_load_pd(&pepsi2[0]);
                       register __m512d z1r,z1i,z2r,z2i;
                       register __m512d cost,numr,numi,denr,deni;
                       register __m512d resr,resi;
                       zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcos(tht);
                       zi_f716_zmm8r8(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_fmsub_pd(z2r,cost,z1r);
                       denr = _mm512_fmadd_pd(z2r,cost,z2r);
                       numi = _mm512_fmsub_pd(z2i,cost,z1i);
                       deni = _mm512_fmadd_pd(z2i,cost,z2i);
                       cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7129_zmm8r8_u(const double * __restrict  ptht,
                                          const double * __restrict  pmur1,
                                          const double * __restrict  pmui1,
                                          const double * __restrict  pepsr1,
                                          const double * __restrict  pepsi1,
                                          const double * __restrict  pmur2,
                                          const double * __restrict  pmui2,
                                          const double * __restrict  pepsr2,
                                          const double * __restrict  pepsi2,
                                          double * __restrict  Rr,
                                          double * __restrict  Ri) {

                       register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                       register __m512d mur1  = _mm512_loadu_pd(&pmur1[0]);
                       register __m512d mui1  = _mm512_loadu_pd(&pmui1[0]);
                       register __m512d epsr1 = _mm512_loadu_pd(&pepsr1[0]);
                       register __m512d epsi1 = _mm512_loadu_pd(&pepsi1[0]);
                       register __m512d mur2  = _mm512_loadu_pd(&pmur2[0]);
                       register __m512d mui2  = _mm512_loadu_pd(&pmui2[0]);
                       register __m512d epsr2 = _mm512_loadu_pd(&pepsr2[0]);
                       register __m512d epsi2 = _mm512_loadu_pd(&pepsi2[0]);
                       register __m512d z1r,z1i,z2r,z2i;
                       register __m512d cost,numr,numi,denr,deni;
                       register __m512d resr,resi;
                       zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcos(tht);
                       zi_f716_zmm8r8(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_fmsub_pd(z2r,cost,z1r);
                       denr = _mm512_fmadd_pd(z2r,cost,z2r);
                       numi = _mm512_fmsub_pd(z2i,cost,z1i);
                       deni = _mm512_fmadd_pd(z2i,cost,z2i);
                       cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
                }


                  /*
                             For (k1/k2)^2*sin^2(theta)<<1 (Simplification
                             of formulae 7.1-9 and 7.1-10).
                             Formula 7.1-30

                     */

   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7130_zmm8r8(const __m512d tht,
                                        const __m512d mur1,
                                        const __m512d mui1,
                                        const __m512d epsr1,
                                        const __m512d epsi1,
                                        const __m512d mur2,
                                        const __m512d mui2,
                                        const __m512d epsr2,
                                        const __m512d epsi2,
                                        __m512d * __restrict Rr,
                                        __m512d * __restrict Ri) {

                       register __m512d z1r,z1i,z2r,z2i;
                       register __m512d cost,numr,numi,denr,deni;
                       register __m512d resr,resi,t0r,t0i;
                       zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcos(tht);
                       t0r  = _mm512_mul_pd(z1r,cost);
                       t0i  = _mm512_mul_pd(z1i,cost);
                       zi_f716_zmm8r8(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_sub_pd(z2r,t0r);
                       denr = _mm512_add_pd(z2r,t0r);
                       numi = _mm512_sub_pd(z2i,t0i);
                       deni = _mm512_add_pd(z2i,t0i);
                       cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7130_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur1,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui1,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr1,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi1,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur2,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui2,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr2,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi2,
                                          double * __restrict __ATTR_ALIGN__(64) Rr,
                                          double * __restrict __ATTR_ALIGN__(64) Ri) {

                       register __m512d tht   = _mm512_load_pd(&ptht[0]);
                       register __m512d mur1  = _mm512_load_pd(&pmur1[0]);
                       register __m512d mui1  = _mm512_load_pd(&pmui1[0]);
                       register __m512d epsr1 = _mm512_load_pd(&pepsr1[0]);
                       register __m512d epsi1 = _mm512_load_pd(&pepsi1[0]);
                       register __m512d mur2  = _mm512_load_pd(&pmur2[0]);
                       register __m512d mui2  = _mm512_load_pd(&pmui2[0]);
                       register __m512d epsr2 = _mm512_load_pd(&pepsr2[0]);
                       register __m512d epsi2 = _mm512_load_pd(&pepsi2[0]);
                       register __m512d z1r,z1i,z2r,z2i;
                       register __m512d cost,numr,numi,denr,deni;
                       register __m512d resr,resi,t0r,t0i;
                       zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcos(tht);
                       t0r  = _mm512_mul_pd(z1r,cost);
                       t0i  = _mm512_mul_pd(z1i,cost);
                       zi_f716_zmm8r8(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_sub_pd(z2r,t0r);
                       denr = _mm512_add_pd(z2r,t0r);
                       numi = _mm512_sub_pd(z2i,t0i);
                       deni = _mm512_add_pd(z2i,t0i);
                       cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                       _mm512_store_pd(&Rr[0] ,resr);
                       _mm512_store_pd(&Ri[0] ,resi);
               }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7130_zmm8r8_u(const double * __restrict  ptht,
                                          const double * __restrict  pmur1,
                                          const double * __restrict  pmui1,
                                          const double * __restrict  pepsr1,
                                          const double * __restrict  pepsi1,
                                          const double * __restrict  pmur2,
                                          const double * __restrict  pmui2,
                                          const double * __restrict  pepsr2,
                                          const double * __restrict  pepsi2,
                                          double * __restrict  Rr,
                                          double * __restrict  Ri) {

                       register __m512d tht   = _mm512_loadu_pd(&ptht[0]);
                       register __m512d mur1  = _mm512_loadu_pd(&pmur1[0]);
                       register __m512d mui1  = _mm512_loadu_pd(&pmui1[0]);
                       register __m512d epsr1 = _mm512_loadu_pd(&pepsr1[0]);
                       register __m512d epsi1 = _mm512_loadu_pd(&pepsi1[0]);
                       register __m512d mur2  = _mm512_loadu_pd(&pmur2[0]);
                       register __m512d mui2  = _mm512_loadu_pd(&pmui2[0]);
                       register __m512d epsr2 = _mm512_loadu_pd(&pepsr2[0]);
                       register __m512d epsi2 = _mm512_loadu_pd(&pepsi2[0]);
                       register __m512d z1r,z1i,z2r,z2i;
                       register __m512d cost,numr,numi,denr,deni;
                       register __m512d resr,resi,t0r,t0i;
                       zi_f716_zmm8r8(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcos(tht);
                       t0r  = _mm512_mul_pd(z1r,cost);
                       t0i  = _mm512_mul_pd(z1i,cost);
                       zi_f716_zmm8r8(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_sub_pd(z2r,t0r);
                       denr = _mm512_add_pd(z2r,t0r);
                       numi = _mm512_sub_pd(z2i,t0i);
                       deni = _mm512_add_pd(z2i,t0i);
                       cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                       _mm512_storeu_pd(&Rr[0] ,resr);
                       _mm512_storeu_pd(&Ri[0] ,resi);
               }


                 /*
                       Reflection coefficients for (alpha<cos^2(theta)).
                       Formula 7.2-15
                  */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7215_f7216_zmm8r8(const __m512d d,
                                                const __m512d k0,
                                                const __m512d alp,
                                                const __m512d tht) {

                         const __m512d vmin = _mm512_set1_pd(1.17549e-38);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __mmask16 m = _mm512_cmp_pd_mask(d,vmin,_CMP_EQ_OQ);
                         register __m512d pid2,cost,cos2t,num,den,x0,x1,x2;
                         register __m512d k,k01a,sin2t,k02k,sqr;
                         register __m512d R;
                         if(!m) {
                            pid2 = _mm512_mul_pd(pi,_mm512_mul_pd(d,hlf));
                            cost = xcos(tht);
                            cos2t= _mm512_fmadd_pd(cost,cost,alp);
                            x0   = _mm512_sqrt_pd(cos2t);
                            x1   = _mm512_fmsub_pd(pid2,cost,x0);
                            num  = _mm512_sinh_pd(x1);
                            x2   = _mm512_fmadd_pd(pid2,cost,x0);
                            den  = _mm512_sinh_pd(x2);
                            R    = _mm512_div_pd(num,den);
                            return (R);
                        }
                        else {
                            const __m512d _1 = _mm512_sqrt_pd(1.0f);
                            k    = _mm512_sqrt_pd(_mm512_sub_pd(_1,alp));
                            cost = xcos(tht);
                            k02k = _mm512_div_pd(_mm512_mul_pd(k0,k0),k);
                            sint = xsin(tht);
                            sin2t= _mm512_mul_pd(sint,sint);
                            x0   = _mm512_sub_pd(_1,_mm512_mul_pd(k02k,sin2t));
                            sqr  = _mm512_sqrt_pd(x0);
                            x1   = _mm512_mul_pd(k,sqr);
                            num  = _mm512_fmsub_pd(k0,cost,x1);
                            den  = _mm512_fmadd_pd(k0,cost,x1);
                            R    = _mm512_div_pd(num,den);
                            return (R);
                        }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7215_f7216_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pd,
                                                  const double * __restrict __ATTR_ALIGN__(64) pk0,
                                                  const double * __restrict __ATTR_ALIGN__(64) palp,
                                                  const double * __restrict __ATTR_ALIGN__(64) ptht) {

                         register __m512d d   = _mm512_load_pd(&pd[0]);
                         register __m512d k0  = _mm512_load_pd(&pk0[0]);
                         register __m512d alp = _mm512_load_pd(&palp[0]);
                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         const __m512d vmin = _mm512_set1_pd(1.17549e-38);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __mmask16 m = _mm512_cmp_pd_mask(d,vmin,_CMP_EQ_OQ);
                         register __m512d pid2,cost,cos2t,num,den,x0,x1,x2;
                         register __m512d k,k01a,sin2t,k02k,sqr;
                         register __m512d R;
                         if(!m) {
                            pid2 = _mm512_mul_pd(pi,_mm512_mul_pd(d,hlf));
                            cost = xcos(tht);
                            cos2t= _mm512_fmadd_pd(cost,cost,alp);
                            x0   = _mm512_sqrt_pd(cos2t);
                            x1   = _mm512_fmsub_pd(pid2,cost,x0);
                            num  = _mm512_sinh_pd(x1);
                            x2   = _mm512_fmadd_pd(pid2,cost,x0);
                            den  = _mm512_sinh_pd(x2);
                            R    = _mm512_div_pd(num,den);
                            return (R);
                        }
                        else {
                            const __m512d _1 = _mm512_sqrt_pd(1.0f);
                            k    = _mm512_sqrt_pd(_mm512_sub_pd(_1,alp));
                            cost = xcos(tht);
                            k02k = _mm512_div_pd(_mm512_mul_pd(k0,k0),k);
                            sint = xsin(tht);
                            sin2t= _mm512_mul_pd(sint,sint);
                            x0   = _mm512_sub_pd(_1,_mm512_mul_pd(k02k,sin2t));
                            sqr  = _mm512_sqrt_pd(x0);
                            x1   = _mm512_mul_pd(k,sqr);
                            num  = _mm512_fmsub_pd(k0,cost,x1);
                            den  = _mm512_fmadd_pd(k0,cost,x1);
                            R    = _mm512_div_pd(num,den);
                            return (R);
                        }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d R_f7215_f7216_zmm8r8_u(const double * __restrict  pd,
                                                  const double * __restrict  pk0,
                                                  const double * __restrict  palp,
                                                  const double * __restrict  ptht) {

                         register __m512d d   = _mm512_loadu_pd(&pd[0]);
                         register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                         register __m512d alp = _mm512_loadu_pd(&palp[0]);
                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         const __m512d vmin = _mm512_set1_pd(1.17549e-38);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __mmask16 m = _mm512_cmp_pd_mask(d,vmin,_CMP_EQ_OQ);
                         register __m512d pid2,cost,cos2t,num,den,x0,x1,x2;
                         register __m512d k,k01a,sin2t,k02k,sqr;
                         register __m512d R;
                         if(!m) {
                            pid2 = _mm512_mul_pd(pi,_mm512_mul_pd(d,hlf));
                            cost = xcos(tht);
                            cos2t= _mm512_fmadd_pd(cost,cost,alp);
                            x0   = _mm512_sqrt_pd(cos2t);
                            x1   = _mm512_fmsub_pd(pid2,cost,x0);
                            num  = _mm512_sinh_pd(x1);
                            x2   = _mm512_fmadd_pd(pid2,cost,x0);
                            den  = _mm512_sinh_pd(x2);
                            R    = _mm512_div_pd(num,den);
                            return (R);
                        }
                        else {
                            const __m512d _1 = _mm512_sqrt_pd(1.0f);
                            k    = _mm512_sqrt_pd(_mm512_sub_pd(_1,alp));
                            cost = xcos(tht);
                            k02k = _mm512_div_pd(_mm512_mul_pd(k0,k0),k);
                            sint = xsin(tht);
                            sin2t= _mm512_mul_pd(sint,sint);
                            x0   = _mm512_sub_pd(_1,_mm512_mul_pd(k02k,sin2t));
                            sqr  = _mm512_sqrt_pd(x0);
                            x1   = _mm512_mul_pd(k,sqr);
                            num  = _mm512_fmsub_pd(k0,cost,x1);
                            den  = _mm512_fmadd_pd(k0,cost,x1);
                            R    = _mm512_div_pd(num,den);
                            return (R);
                        }
                  }


                    /*
                            Infinite strips, low frequency region.
                            E-field (scattered) along 'z'.
                            Formula 7.4-1
                        */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f741_zmm8r8(const __m512d k0,
                                         const __m512d r,
                                         const __m512d a,
                                         const __m512d tht,
                                         const __m512d Eir,
                                         const __m512d Eii,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {

                        const __m512d _1   = _mm512_set1_pd(1.0f);
                        const __m512d gam  = _mm512_set1_pd(1.7811f);
                        const __m512d qtr  = _mm512_set1_pd(0.25f);
                        const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                        const __m512d pi2  = _mm512_set1_pd(0.5f*3.14159265358979323846264338328);
                        const __m512d _4   = _mm512_set1_pd(4.0f);
                        const __m512d pi4  = _mm512_set1_pd(0.25f*3.14159265358979323846264338328);
                        register __m512d denr,deni,ear,eai,cer,cei,arg;
                        register __m512d t0r,t0i,num,k02,a2,k0a,k0r,cost,trm;
                        register __m512d x0,x1,t1r,t1i;
                        deni = pi2;
                        cost = xcos(tht);
                        k0r  = _mm512_mul_pd(k0,r);
                        a2   = _mm512_mul_pd(a,a);
                        ear  = _mm512_setzero_pd();
                        k0a  = _mm512_mul_pd(k0,a);
                        eai  = _mm512_add_pd(k0r,pi4);
                        k02  = _mm512_mul_pd(k0,k0);
                        x0   = _mm512_div_pd(pi,_mm512_add_pd(k0r,k0r));
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        arg  = _mm512_div_pd(_4,_mm512_mul_pd(gam,k0a));
                        trm  = _mm512_sqrt_pd(x0);
                        x1   = _mm512_mul_pd(cost,cost);
                        denr = xlog(arg);
                        x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                        num  = _mm512_mul_pd(x0,x1);
                        t0r  = _mm512_mul_pd(_mm512_div_pd(num,denr),trm);
                        t0i  = _mm512_mul_pd(_mm512_div_pd(num,deni),trm);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        cmul_zmm8r8(Eir,Eii,t1r,t1i,*Esr,*Esi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f741_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) pEir,
                                           const double * __restrict __ATTR_ALIGN__(64) pEii,
                                           double * __restrict __ATTR_ALIGN__(64) Esr,
                                           double * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d a   = _mm512_load_pd(&pa[0]);
                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d Eir = _mm512_load_pd(&pEir[0]); 
                        register __m512d Eii = _mm512_load_pd(&pEii[0]);
                        const __m512d _1   = _mm512_set1_pd(1.0f);
                        const __m512d gam  = _mm512_set1_pd(1.7811f);
                        const __m512d qtr  = _mm512_set1_pd(0.25f);
                        const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                        const __m512d pi2  = _mm512_set1_pd(0.5f*3.14159265358979323846264338328);
                        const __m512d _4   = _mm512_set1_pd(4.0f);
                        const __m512d pi4  = _mm512_set1_pd(0.25f*3.14159265358979323846264338328);
                        register __m512d denr,deni,ear,eai,cer,cei,arg;
                        register __m512d t0r,t0i,num,k02,a2,k0a,k0r,cost,trm;
                        register __m512d x0,x1,t1r,t1i,resr,resi;
                        deni = pi2;
                        cost = xcos(tht);
                        k0r  = _mm512_mul_pd(k0,r);
                        a2   = _mm512_mul_pd(a,a);
                        ear  = _mm512_setzero_pd();
                        k0a  = _mm512_mul_pd(k0,a);
                        eai  = _mm512_add_pd(k0r,pi4);
                        k02  = _mm512_mul_pd(k0,k0);
                        x0   = _mm512_div_pd(pi,_mm512_add_pd(k0r,k0r));
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        arg  = _mm512_div_pd(_4,_mm512_mul_pd(gam,k0a));
                        trm  = _mm512_sqrt_pd(x0);
                        x1   = _mm512_mul_pd(cost,cost);
                        denr = xlog(arg);
                        x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                        num  = _mm512_mul_pd(x0,x1);
                        t0r  = _mm512_mul_pd(_mm512_div_pd(num,denr),trm);
                        t0i  = _mm512_mul_pd(_mm512_div_pd(num,deni),trm);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        cmul_zmm8r8(Eir,Eii,t1r,t1i,&resr,&resi);
                        _mm512_store_pd(&Esr[0], resr);
                        _mm512_store_pd(&Esi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f741_zmm8r8_u(const double * __restrict  pk0,
                                           const double * __restrict  pr,
                                           const double * __restrict  pa,
                                           const double * __restrict  ptht,
                                           const double * __restrict  pEir,
                                           const double * __restrict  pEii,
                                           double * __restrict  Esr,
                                           double * __restrict  Esi) {

                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d a   = _mm512_loadu_pd(&pa[0]);
                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d Eir = _mm512_loadu_pd(&pEir[0]); 
                        register __m512d Eii = _mm512_loadu_pd(&pEii[0]);
                        const __m512d _1   = _mm512_set1_pd(1.0f);
                        const __m512d gam  = _mm512_set1_pd(1.7811f);
                        const __m512d qtr  = _mm512_set1_pd(0.25f);
                        const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                        const __m512d pi2  = _mm512_set1_pd(0.5f*3.14159265358979323846264338328);
                        const __m512d _4   = _mm512_set1_pd(4.0f);
                        const __m512d pi4  = _mm512_set1_pd(0.25f*3.14159265358979323846264338328);
                        register __m512d denr,deni,ear,eai,cer,cei,arg;
                        register __m512d t0r,t0i,num,k02,a2,k0a,k0r,cost,trm;
                        register __m512d x0,x1,t1r,t1i,resr,resi;
                        deni = pi2;
                        cost = xcos(tht);
                        k0r  = _mm512_mul_pd(k0,r);
                        a2   = _mm512_mul_pd(a,a);
                        ear  = _mm512_setzero_pd();
                        k0a  = _mm512_mul_pd(k0,a);
                        eai  = _mm512_add_pd(k0r,pi4);
                        k02  = _mm512_mul_pd(k0,k0);
                        x0   = _mm512_div_pd(pi,_mm512_add_pd(k0r,k0r));
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        arg  = _mm512_div_pd(_4,_mm512_mul_pd(gam,k0a));
                        trm  = _mm512_sqrt_pd(x0);
                        x1   = _mm512_mul_pd(cost,cost);
                        denr = xlog(arg);
                        x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                        num  = _mm512_mul_pd(x0,x1);
                        t0r  = _mm512_mul_pd(_mm512_div_pd(num,denr),trm);
                        t0i  = _mm512_mul_pd(_mm512_div_pd(num,deni),trm);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        cmul_zmm8r8(Eir,Eii,t1r,t1i,&resr,&resi);
                        _mm512_storeu_pd(&Esr[0], resr);
                        _mm512_storeu_pd(&Esi[0], resi);
                 }


                     
                    /*
                            Infinite strips, low frequency region.
                            H-field (scattered) along 'z'.
                            Formula 7.4-2
                        */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f742_zmm8r8(const __m512d k0a,
                                         const __m512d k0r,
                                         const __m512d tht,
                                         const __m512d Hir,
                                         const __m512d Hii,
                                         __m512d * __restrict Hsr,
                                         __m512d * __restrict Hsi) {

                         const __m512d _1o8 = _mm512_set1_pd(0.125f);
                         const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d pi4  = _mm512_set1_pd(0.25f*3.14159265358979323846264338328);
                         register __m512d ear,eai,cer,cei;
                         register __m512d trm,t0r,t0i,num,cost,x0,x1;
                         ear  = _mm512_setzero_pd();
                         cost = xcos(tht);
                         eai  = _mm512_add_pd(k0r,pi4);
                         x0   = _mm512_div_pd(pi,_mm512_add_pd(k0r,k0r));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         trm  = _mm512_sqrt_pd(x0);
                         x1   = _mm512_add_pd(k0a,k0a);
                         x0   = _mm512_mul_pd(cost,cost);
                         num  = _mm512_mul_pd(_mm512_mul_pd(x1,x1),x0);
                         t0r  = _mm512_mul_pd(trm,cer);
                         t0i  = _mm512_mul_pd(trm,cei);
                         num  = _mm512_mul_pd(_1o8,num);
                         x0   = _mm512_mul_pd(Hsr,num);
                         x1   = _mm512_mul_pd(Hsi,num);
                         cmul_zmm8r8(x0,x1,t0r,t0i,*Hsr,*Hsi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f742_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) pHir,
                                           const double * __restrict __ATTR_ALIGN__(64) pHii,
                                           double * __restrict __ATTR_ALIGN__(64) Hsr,
                                           double * __restrict __ATTR_ALIGN__(64) Hsi) {

                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         register __m512d k0r = _mm512_load_pd(&pk0r[0]);
                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         register __m512d Hir = _mm512_load_pd(&pHir[0]);
                         register __m512d Hii = _mm512_load_pd(&pHii[0]);
                         const __m512d _1o8 = _mm512_set1_pd(0.125f);
                         const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d pi4  = _mm512_set1_pd(0.25f*3.14159265358979323846264338328);
                         register __m512d ear,eai,cer,cei;
                         register __m512d trm,t0r,t0i,num,cost,x0,x1,resr,resi;
                         ear  = _mm512_setzero_pd();
                         cost = xcos(tht);
                         eai  = _mm512_add_pd(k0r,pi4);
                         x0   = _mm512_div_pd(pi,_mm512_add_pd(k0r,k0r));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         trm  = _mm512_sqrt_pd(x0);
                         x1   = _mm512_add_pd(k0a,k0a);
                         x0   = _mm512_mul_pd(cost,cost);
                         num  = _mm512_mul_pd(_mm512_mul_pd(x1,x1),x0);
                         t0r  = _mm512_mul_pd(trm,cer);
                         t0i  = _mm512_mul_pd(trm,cei);
                         num  = _mm512_mul_pd(_1o8,num);
                         x0   = _mm512_mul_pd(Hsr,num);
                         x1   = _mm512_mul_pd(Hsi,num);
                         cmul_zmm8r8(x0,x1,t0r,t0i,&resr,&resi);
                         _mm512_store_pd(&Hsr[0], resr);
                         _mm512_store_pd(&Hsi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f742_zmm8r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  pk0r,
                                           const double * __restrict  ptht,
                                           const double * __restrict  pHir,
                                           const double * __restrict  pHii,
                                           double * __restrict  Hsr,
                                           double * __restrict  Hsi) {

                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d k0r = _mm512_loadu_pd(&pk0r[0]);
                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         register __m512d Hir = _mm512_loadu_pd(&pHir[0]);
                         register __m512d Hii = _mm512_loadu_pd(&pHii[0]);
                         const __m512d _1o8 = _mm512_set1_pd(0.125f);
                         const __m512d pi   = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d pi4  = _mm512_set1_pd(0.25f*3.14159265358979323846264338328);
                         register __m512d ear,eai,cer,cei,resr,resi;
                         register __m512d trm,t0r,t0i,num,cost,x0,x1;
                         ear  = _mm512_setzero_pd();
                         cost = xcos(tht);
                         eai  = _mm512_add_pd(k0r,pi4);
                         x0   = _mm512_div_pd(pi,_mm512_add_pd(k0r,k0r));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         trm  = _mm512_sqrt_pd(x0);
                         x1   = _mm512_add_pd(k0a,k0a);
                         x0   = _mm512_mul_pd(cost,cost);
                         num  = _mm512_mul_pd(_mm512_mul_pd(x1,x1),x0);
                         t0r  = _mm512_mul_pd(trm,cer);
                         t0i  = _mm512_mul_pd(trm,cei);
                         num  = _mm512_mul_pd(_1o8,num);
                         x0   = _mm512_mul_pd(Hsr,num);
                         x1   = _mm512_mul_pd(Hsi,num);
                         cmul_zmm8r8(x0,x1,t0r,t0i,&resr,&resi);
                         _mm512_storeu_pd(&Hsr[0], resr);
                         _mm512_storeu_pd(&Hsi[0], resi);
                }


                  /*
                       The resultant backscatter RCS of perpendicular
                       polarization.
                       Formula 7.4-3
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f743_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht) {

                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d pisq= _mm512_set1_pd(0.25f*9.869604401089358618834490999876);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d c0  = _mm512_set1_pd(4.48f);
                          register __m512d fac,num,den,cost,k02,a2,k0a;
                          register __m512d rcs,x0,x1,arg,larg,rat;
                          k0a  = _mm512_mul_pd(k0,a);
                          fac  = _mm512_div_pd(pis,k0);
                          k02  = _mm512_mul_pd(k0,k0);
                          a2   = _mm512_mul_pd(a,a);
                          cost = xcos(tht);
                          x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                          arg  = _mm512_div_pd(c0,_mm512_add_pd(k0a,k0a));
                          x1   = _mm512_mul_pd(x0,_mm512_mul_pd(cost,cost));
                          larg = xlog(arg);
                          num  = _mm512_mul_pd(x1,x1);
                          den  = _mm512_fmadd_pd(larg,larg,pisq);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512mul_pd(fac,rat);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f743_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0 = _mm512_load_pd(&pk0[0]);
                          register __m512d a  = _mm512_load_pd(&pa[0]);
                          register __m512d tht= _mm512_load_pd(&ptht[0]);
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d pisq= _mm512_set1_pd(0.25f*9.869604401089358618834490999876);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d c0  = _mm512_set1_pd(4.48f);
                          register __m512d fac,num,den,cost,k02,a2,k0a;
                          register __m512d rcs,x0,x1,arg,larg,rat;
                          k0a  = _mm512_mul_pd(k0,a);
                          fac  = _mm512_div_pd(pis,k0);
                          k02  = _mm512_mul_pd(k0,k0);
                          a2   = _mm512_mul_pd(a,a);
                          cost = xcos(tht);
                          x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                          arg  = _mm512_div_pd(c0,_mm512_add_pd(k0a,k0a));
                          x1   = _mm512_mul_pd(x0,_mm512_mul_pd(cost,cost));
                          larg = xlog(arg);
                          num  = _mm512_mul_pd(x1,x1);
                          den  = _mm512_fmadd_pd(larg,larg,pisq);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512_mul_pd(fac,rat);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f743_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa,
                                             const double * __restrict  ptht) {

                          register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a  = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht= _mm512_loadu_pd(&ptht[0]);
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d pisq= _mm512_set1_pd(0.25f*9.869604401089358618834490999876);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d c0  = _mm512_set1_pd(4.48f);
                          register __m512d fac,num,den,cost,k02,a2,k0a;
                          register __m512d rcs,x0,x1,arg,larg,rat;
                          k0a  = _mm512_mul_pd(k0,a);
                          fac  = _mm512_div_pd(pis,k0);
                          k02  = _mm512_mul_pd(k0,k0);
                          a2   = _mm512_mul_pd(a,a);
                          cost = xcos(tht);
                          x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                          arg  = _mm512_div_pd(c0,_mm512_add_pd(k0a,k0a));
                          x1   = _mm512_mul_pd(x0,_mm512_mul_pd(cost,cost));
                          larg = xlog(arg);
                          num  = _mm512_mul_pd(x1,x1);
                          den  = _mm512_fmadd_pd(larg,larg,pisq);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512_mul_pd(fac,rat);
                          return (rcs);
                }


                  /*
                        
                       The resultant backscatter RCS of parallel
                       polarization.
                       Formula 7.4-4  

                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f744_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht) {

                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876); 
                          const __m512d 1o64= _mm512_set1_pd(0.015625f);
                          register __m512d k0a,k0a2,cost,cos2t,fac;
                          register __m512d rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          k0a2 = _mm512_add_pd(k0a,k0a);
                          fac  = _mm512_div_pd(pis,k0);
                          x2   = _mm512_mul_pd(cos2t,cos2t);
                          x0   = _mm512_mul_pd(k0a2,k0a2);
                          x1   = _mm512_mul_pd(x0,x0);
                          num  = _mm512_mul_pd(_mm512_mul_pd(x1,x2),1o64);
                          rcs  = _mm512_mul_pd(fac,num);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f744_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876); 
                          const __m512d 1o64= _mm512_set1_pd(0.015625);
                          register __m512d k0a,k0a2,cost,cos2t,fac;
                          register __m512d rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          k0a2 = _mm512_add_pd(k0a,k0a);
                          fac  = _mm512_div_pd(pis,k0);
                          x2   = _mm512_mul_pd(cos2t,cos2t);
                          x0   = _mm512_mul_pd(k0a2,k0a2);
                          x1   = _mm512_mul_pd(x0,x0);
                          num  = _mm512_mul_pd(_mm512_mul_pd(x1,x2),1o64);
                          rcs  = _mm512_mul_pd(fac,num);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f744_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa,
                                             const double * __restrict  ptht) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876); 
                          const __m512d 1o64= _mm512_set1_pd(0.015625);
                          register __m512d k0a,k0a2,cost,cos2t,fac;
                          register __m512d rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          k0a2 = _mm512_add_pd(k0a,k0a);
                          fac  = _mm512_div_pd(pis,k0);
                          x2   = _mm512_mul_pd(cos2t,cos2t);
                          x0   = _mm512_mul_pd(k0a2,k0a2);
                          x1   = _mm512_mul_pd(x0,x0);
                          num  = _mm512_mul_pd(_mm512_mul_pd(x1,x2),1o64);
                          rcs  = _mm512_mul_pd(fac,num);
                          return (rcs);
                }


                  /*
                          General bistatic case.
                          The Rayleigh scattering results.
                          Plane-perpendicular.
                          Formula 7.4-5
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f745_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht,
                                           const __m512d tht2) {

                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d pisq= _mm512_set1_pd(0.25f*9.869604401089358618834490999876);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d c0  = _mm512_set1_pd(4.48f);
                          register __m512d fac,num,den,cost,k02,a2,k0a;
                          register __m512d rcs,x0,x1,arg,larg,rat,cost2;
                          k0a  = _mm512_mul_pd(k0,a);
                          fac  = _mm512_div_pd(pis,k0);
                          k02  = _mm512_mul_pd(k0,k0);
                          a2   = _mm512_mul_pd(a,a);
                          cost = xcos(tht);
                          x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                          cost2= xcos(tht2);
                          arg  = _mm512_div_pd(c0,_mm512_add_pd(k0a,k0a));
                          x1   = _mm512_mul_pd(x0,_mm512_mul_pd(cost,cost2));
                          larg = xlog(arg);
                          num  = _mm512_mul_pd(x1,x1);
                          den  = _mm512_fmadd_pd(larg,larg,pisq);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512mul_pd(fac,rat);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f745_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht2) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d tht2= _mm512_load_pd(&ptht2[0]);  
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d pisq= _mm512_set1_pd(0.25f*9.869604401089358618834490999876);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d c0  = _mm512_set1_pd(4.48f);
                          register __m512d fac,num,den,cost,k02,a2,k0a;
                          register __m512d rcs,x0,x1,arg,larg,rat,cost2;
                          k0a  = _mm512_mul_pd(k0,a);
                          fac  = _mm512_div_pd(pis,k0);
                          k02  = _mm512_mul_pd(k0,k0);
                          a2   = _mm512_mul_pd(a,a);
                          cost = xcos(tht);
                          x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                          cost2= xcos(tht2);
                          arg  = _mm512_div_pd(c0,_mm512_add_pd(k0a,k0a));
                          x1   = _mm512_mul_pd(x0,_mm512_mul_pd(cost,cost2));
                          larg = xlog(arg);
                          num  = _mm512_mul_pd(x1,x1);
                          den  = _mm512_fmadd_pd(larg,larg,pisq);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512mul_pd(fac,rat);
                          return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f745_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict pa,
                                             const double * __restrict  ptht,
                                             const double * __restrict  ptht2) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register __m512d tht2= _mm512_loadu_pd(&ptht2[0]);  
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d pisq= _mm512_set1_pd(0.25f*9.869604401089358618834490999876);
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d c0  = _mm512_set1_pd(4.48f);
                          register __m512d fac,num,den,cost,k02,a2,k0a;
                          register __m512d rcs,x0,x1,arg,larg,rat,cost2;
                          k0a  = _mm512_mul_pd(k0,a);
                          fac  = _mm512_div_pd(pis,k0);
                          k02  = _mm512_mul_pd(k0,k0);
                          a2   = _mm512_mul_pd(a,a);
                          cost = xcos(tht);
                          x0   = _mm512_fmadd_pd(k02,_mm512_mul_pd(a2,qtr),_1);
                          cost2= xcos(tht2);
                          arg  = _mm512_div_pd(c0,_mm512_add_pd(k0a,k0a));
                          x1   = _mm512_mul_pd(x0,_mm512_mul_pd(cost,cost2));
                          larg = xlog(arg);
                          num  = _mm512_mul_pd(x1,x1);
                          den  = _mm512_fmadd_pd(larg,larg,pisq);
                          rat  = _mm512_div_pd(num,den);
                          rcs  = _mm512mul_pd(fac,rat);
                          return (rcs);
                }


                   /*
                          General bistatic case.
                          The Rayleigh scattering results.
                          Plane-parallel.
                          Formula 7.4-6
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f746_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht,
                                           const __m512d tht2) {

                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876); 
                          const __m512d 1o64= _mm512_set1_pd(0.015625);
                          register __m512d k0a,k0a2,cost,cos2t,fac,cost2;
                          register __m512d rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          k0a2 = _mm512_add_pd(k0a,k0a);
                          cost2= xcos(tht2);
                          fac  = _mm512_div_pd(pis,k0);
                          x2   = _mm512_mul_pd(_mm512_mul_pd(cost2.cost2),
                                 _mm512_mul_pd(cos2t,cos2t));
                          x0   = _mm512_mul_pd(k0a2,k0a2);
                          x1   = _mm512_mul_pd(x0,x0);
                          num  = _mm512_mul_pd(_mm512_mul_pd(x1,x2),1o64);
                          rcs  = _mm512_mul_pd(fac,num);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f746_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht2) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d tht2= _mm512_load_pd(&ptht2[0]);
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876); 
                          const __m512d 1o64= _mm512_set1_pd(0.015625);
                          register __m512d k0a,k0a2,cost,cos2t,fac,cost2;
                          register __m512d rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          k0a2 = _mm512_add_pd(k0a,k0a);
                          cost2= xcos(tht2);
                          fac  = _mm512_div_pd(pis,k0);
                          x2   = _mm512_mul_pd(_mm512_mul_pd(cost2.cost2),
                                 _mm512_mul_pd(cos2t,cos2t));
                          x0   = _mm512_mul_pd(k0a2,k0a2);
                          x1   = _mm512_mul_pd(x0,x0);
                          num  = _mm512_mul_pd(_mm512_mul_pd(x1,x2),1o64);
                          rcs  = _mm512_mul_pd(fac,num);
                          return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f746_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa,
                                             const double * __restrict  ptht,
                                             const double * __restrict  ptht2) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register __m512d tht2= _mm512_loadu_pd(&ptht2[0]);
                          const __m512d pis = _mm512_set1_pd(9.869604401089358618834490999876); 
                          const __m512d 1o64= _mm512_set1_pd(0.015625);
                          register __m512d k0a,k0a2,cost,cos2t,fac,cost2;
                          register __m512d rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          k0a2 = _mm512_add_pd(k0a,k0a);
                          cost2= xcos(tht2);
                          fac  = _mm512_div_pd(pis,k0);
                          x2   = _mm512_mul_pd(_mm512_mul_pd(cost2.cost2),
                                 _mm512_mul_pd(cos2t,cos2t));
                          x0   = _mm512_mul_pd(k0a2,k0a2);
                          x1   = _mm512_mul_pd(x0,x0);
                          num  = _mm512_mul_pd(_mm512_mul_pd(x1,x2),1o64);
                          rcs  = _mm512_mul_pd(fac,num);
                          return (rcs);
                }


                  /*
                         High Frequency Region.
                         For k0a>>1, PO solution of backscatter RCS.
                         Formula 7.4-7
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f747_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht) {

                          register __m512d invk0,k0a,cost,sint,arg,sarg;
                          register __m512d rcs,num,sqr,x0;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          invk0= _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          arg  = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          sarg = xsin(arg);
                          num  = _mm512_mul_pd(cost,sarg);
                          x0   = _mm512_div_pd(num,sint);
                          sqr  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_mul_pd(invk0,sqr);
                          return (rcs); 
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f747_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register __m512d invk0,k0a,cost,sint,arg,sarg;
                          register __m512d rcs,num,sqr,x0;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          invk0= _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          arg  = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          sarg = xsin(arg);
                          num  = _mm512_mul_pd(cost,sarg);
                          x0   = _mm512_div_pd(num,sint);
                          sqr  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_mul_pd(invk0,sqr);
                          return (rcs); 
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f747_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa,
                                             const double * __restrict  ptht) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register __m512d invk0,k0a,cost,sint,arg,sarg;
                          register __m512d rcs,num,sqr,x0;
                          k0a  = _mm512_mul_pd(k0,a);
                          cost = xcos(tht);
                          invk0= _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          arg  = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          sarg = xsin(arg);
                          num  = _mm512_mul_pd(cost,sarg);
                          x0   = _mm512_div_pd(num,sint);
                          sqr  = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_mul_pd(invk0,sqr);
                          return (rcs); 
                }


                 /*
                       Backscattered fields from the edges of strips.
                       Helper function for the formula 7.4-9
                       Electric-field (over z).
                       Formula 7.4-15
                  */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefg12_f7415_zmm8r8(  const __m512d k0a,
                                                const __m512d tht,
                                                __m512d * __restrict gamm1,
                                                __m512d * __restrict gamm2){
                                                 
                          const __m512d C0318309886183790671537767526745 = _mm512_set1_pd(0.318309886183790671537767526745);
                          const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d C05   = _mm512_set1_pd(0.5f);
                          register __m512d thth,arg1,arg2,carg1,carg2,sqr,x0;
                          x0   = _mm512_add_pd(k0a,k0a);
                          thth = _mm512_mul_pd(C05,tht);
                          sqr  = _mm512_sqrt_pd(_mm512_mul_pd(x0,C0318309886183790671537767526745));
                          arg1 = _mm512_add_pd(C078539816339744830961566084582,thth);
                          carg1= xcos(arg1);
                          x0   = _mm512_add_pd(sqr,sqr);
                          arg2 = _mm512_sub_pd(C078539816339744830961566084582,thth);
                          carg2= xcos(arg2);
                          *gamm1 = _mm512_mul_pd(x0,_mm512_abs_pd(carg1));
                          *gamm2 = _mm512_mul_pd(x0,_mm512_abs_pd(carg2));
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefg12_f7415_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                                const double * __restrict __ATTR_ALIGN__(64) ptht,
                                                double * __restrict __ATTR_ALIGN__(64) gamm1,
                                                double * __restrict __ATTR_ALIGN__(64) gamm2){
                                
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);                 
                          const __m512d C0318309886183790671537767526745 = _mm512_set1_pd(0.318309886183790671537767526745);
                          const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d C05   = _mm512_set1_pd(0.5f);
                          register __m512d thth,arg1,arg2,carg1,carg2,sqr,x0;
                          x0   = _mm512_add_pd(k0a,k0a);
                          thth = _mm512_mul_pd(C05,tht);
                          sqr  = _mm512_sqrt_pd(_mm512_mul_pd(x0,C0318309886183790671537767526745));
                          arg1 = _mm512_add_pd(C078539816339744830961566084582,thth);
                          carg1= xcos(arg1);
                          x0   = _mm512_add_pd(sqr,sqr);
                          arg2 = _mm512_sub_pd(C078539816339744830961566084582,thth);
                          carg2= xcos(arg2);
                          _mm512_store_pd(&gamm1[0] ,_mm512_mul_pd(x0,_mm512_abs_pd(carg1)));
                          _mm512_store_pd(&gamm2[0] ,_mm512_mul_pd(x0,_mm512_abs_pd(carg2)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefg12_f7415_zmm8r8_u(  const double * __restrict  pk0a,
                                                const double * __restrict  ptht,
                                                double * __restrict  gamm1,
                                                double * __restrict  gamm2){
                                
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);                 
                          const __m512d C0318309886183790671537767526745 = _mm512_set1_pd(0.318309886183790671537767526745);
                          const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d C05   = _mm512_set1_pd(0.5f);
                          register __m512d thth,arg1,arg2,carg1,carg2,sqr,x0;
                          x0   = _mm512_add_pd(k0a,k0a);
                          thth = _mm512_mul_pd(C05,tht);
                          sqr  = _mm512_sqrt_pd(_mm512_mul_pd(x0,C0318309886183790671537767526745));
                          arg1 = _mm512_add_pd(C078539816339744830961566084582,thth);
                          carg1= xcos(arg1);
                          x0   = _mm512_add_pd(sqr,sqr);
                          arg2 = _mm512_sub_pd(C078539816339744830961566084582,thth);
                          carg2= xcos(arg2);
                          _mm512_storeu_pd(&gamm1[0] ,_mm512_mul_pd(x0,_mm512_abs_pd(carg1)));
                          _mm512_storeu_pd(&gamm2[0] ,_mm512_mul_pd(x0,_mm512_abs_pd(carg2)));
                }



                      /*
                       Backscattered fields from the edges of strips.
                       Helper function for the formula 7.4-9
                       Electric-field (over z).
                       Formula 7.4-13
                  */

#include "GMS_rcs_common_zmm8r8.hpp"

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefA12_f7413_zmm8r8(const __m512d k0a,
                                              const __m512d tht,
                                              __m512d * __restrict A1r,
                                              __m512d * __restrict A1i,
                                              __m512d * __restrict A2r,
                                              __m512d * __restrict A2i) {

                        const __m512d C078539816339744830961566084582  = _mm512_set1_pd(-0.78539816339744830961566084582);
                        const __m512d C141421356237309504880168872421  = _mm512_set1_pd(1.41421356237309504880168872421);
                        register __m512d ear,eai,cer,cei,Cr1,Si1,Cr2,Si2;
                        register __m512d gam1,gam2;
                        ear = _mm512_setzero_pd();
                        eai = C078539816339744830961566084582;
                        coefg12_f7415_zmm8r8(k0a,tht,&gam1,&gam2); 
                        Cr1 = fresnel_C_zmm8r8(gam1);
                        Si1 = fresnel_S_zmm8r8(gam1);
                        Cr2 = fresnel_C_zmm8r8(gam2);
                        Si2 = fresnel_S_zmm8r8(gam2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cer = _mm512_mul_pd(C141421356237309504880168872421,cer);
                        cei = _mm512_mul_pd(C141421356237309504880168872421,cei);
                        cmul_zmm8r8(cer,cei,Cr1,Si1,*A1r,*A1i);
                        cmul_zmm8r8(cer,cei,Cr2,Si2,*A2r,*A2i);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefA12_f7413_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              double * __restrict __ATTR_ALIGN__(64) A1r,
                                              double * __restrict __ATTR_ALIGN__(64) A1i,
                                              double * __restrict __ATTR_ALIGN__(64) A2r,
                                              double * __restrict __ATTR_ALIGN__(64) A2i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        const __m512d C078539816339744830961566084582  = _mm512_set1_pd(-0.78539816339744830961566084582);
                        const __m512d C141421356237309504880168872421  = _mm512_set1_pd(1.41421356237309504880168872421);
                        register __m512d ear,eai,cer,cei,Cr1,Si1,Cr2,Si2;
                        register __m512d gam1,gam2,res1r,res1i,res2r,res2i;
                        ear = _mm512_setzero_pd();
                        eai = C078539816339744830961566084582;
                        coefg12_f7415_zmm8r8(k0a,tht,&gam1,&gam2); 
                        Cr1 = fresnel_C_zmm8r8(gam1);
                        Si1 = fresnel_S_zmm8r8(gam1);
                        Cr2 = fresnel_C_zmm8r8(gam2);
                        Si2 = fresnel_S_zmm8r8(gam2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cer = _mm512_mul_pd(C141421356237309504880168872421,cer);
                        cei = _mm512_mul_pd(C141421356237309504880168872421,cei);
                        cmul_zmm8r8(cer,cei,Cr1,Si1,&res1r,&res1i);
                        __m512d_store_pd(&A1r[0], res1r);
                        __m512d_store_pd(&A1i[0], res1i);
                        cmul_zmm8r8(cer,cei,Cr2,Si2,&res2r,&res2i);
                        __m512d_store_pd(&A2r[0], res2r);
                        __m512d_store_pd(&A2i[0], res2i);
               }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefA12_f7413_zmm8r8_u(const double * __restrict  pk0a,
                                                const double * __restrict  ptht,
                                                double * __restrict  A1r,
                                                double * __restrict  A1i,
                                                double * __restrict  A2r,
                                                double * __restrict  A2i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        const __m512d C078539816339744830961566084582  = _mm512_set1_pd(-0.78539816339744830961566084582);
                        const __m512d C141421356237309504880168872421  = _mm512_set1_pd(1.41421356237309504880168872421);
                        register __m512d ear,eai,cer,cei,Cr1,Si1,Cr2,Si2;
                        register __m512d gam1,gam2,res1r,res1i,res2r,res2i;
                        ear = _mm512_setzero_pd();
                        eai = C078539816339744830961566084582;
                        coefg12_f7415_zmm8r8(k0a,tht,&gam1,&gam2); 
                        Cr1 = fresnel_C_zmm8r8(gam1);
                        Si1 = fresnel_S_zmm8r8(gam1);
                        Cr2 = fresnel_C_zmm8r8(gam2);
                        Si2 = fresnel_S_zmm8r8(gam2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cer = _mm512_mul_pd(C141421356237309504880168872421,cer);
                        cei = _mm512_mul_pd(C141421356237309504880168872421,cei);
                        cmul_zmm8r8(cer,cei,Cr1,Si1,&res1r,&res1i);
                        __m512d_storeu_pd(&A1r[0], res1r);
                        __m512d_storeu_pd(&A1i[0], res1i);
                        cmul_zmm8r8(cer,cei,Cr2,Si2,&res2r,&res2i);
                        __m512d_storeu_pd(&A2r[0], res2r);
                        __m512d_storeu_pd(&A2i[0], res2i);
               }


                   /*
                       Backscattered fields from the edges of strips.
                       Helper function for the formula 7.4-9
                       Electric-field (over z).
                       Formula 7.4-14
                  */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefB12_f7414_zmm8r8(const __m512d k0a,
                                              const __m512d tht,
                                              __m512d * __restrict B1r,
                                              __m512d * __restrict B1i,
                                              __m512d * __restrict B2r,
                                              __m512d * __restrict B2i) {

                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(-0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C05                              = _mm512_set1_pd(0.5f);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear1,eai1,eai2,cer1,cei1,cer2,cei2;
                         register __m512d ir,ii,sint,htht,carg1,carg2,arg1,arg2,abs1,abs2;
                         register __m512d x0,x1,x2,x3,x4,x5,t0r,t0i,t1r,t1i,k0a2;
                         register __m512d A1r,A1i,A2r,A2i;
                         ir   = _mm512_setzero_pd();
                         x0   = _mm512_sqrt_pd(_mm512_mul_pd(C314159265358979323846264338328,k0a));
                         coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         htht = _mm512_mul_pd(C05,tht);
                         ear1 = ir;
                         ii   = _mm512_add_pd(x0,x0);
                         k0a2 = _mm512_add_pd(k0a,k0a);
                         sint = xsin(tht);
                         x0   = _mm512_add_pd(C10,sint);
                         eai2 = _mm512_fmsub_pd(k0a2,x0,C078539816339744830961566084582);
                         cexp_zmm8r8(ear1,eai2,&cer1,&cei1);
                         arg1 = _mm512_sub_pd(C078539816339744830961566084582,htht);
                         carg1= xcos(arg1);
                         x1   = _mm512_sub_pd(C10,sint);
                         eai1 = _mm512_fmsub_pd(k0a2,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear1,eai1,&cer2,&cei2);
                         abs1 = _mm512_abs_pd(carg1);
                         t0r  = _mm512_div_pd(cer2,abs1);
                         arg2 = _mm512_add_pd(C078539816339744830961566084582,htht);
                         t0i  = _mm512_div_pd(cei2,abs1);
                         carg2= xcos(arg2);
                         abs2 = _mm512_abs_pd(carg2);
                         t1r  = _mm512_div_pd(cer1,abs2);
                         cmul_zmm8r8(ir,ii,t0r,t0i,&x2,&x3);
                         t1i  = _mm512_div_pd(cei1,abs2);
                         cmul_zmm8r8(ir,ii,t1r,t1i,&x4,&x5);
                         *B1r  = _mm512_add_pd(A1r,x2);
                         *B2r  = _mm512_add_pd(A2r,x4);
                         *B1i  = _mm512_add_pd(A1i,x3);
                         *B2i  = _mm512_add_pd(A2i,x5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefB12_f7414_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0a,
                                                const double * __restrict __ATTR_ALIGN__(64)  ptht,
                                                double * __restrict __ATTR_ALIGN__(64) B1r,
                                                double * __restrict __ATTR_ALIGN__(64) B1i,
                                                double * __restrict __ATTR_ALIGN__(64) B2r,
                                                double * __restrict __ATTR_ALIGN__(64) B2i) {

                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(-0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C05                              = _mm512_set1_pd(0.5f);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear1,eai1,eai2,cer1,cei1,cer2,cei2;
                         register __m512d ir,ii,sint,htht,carg1,carg2,arg1,arg2,abs1,abs2;
                         register __m512d x0,x1,x2,x3,x4,x5,t0r,t0i,t1r,t1i,k0a2;
                         register __m512d A1r,A1i,A2r,A2i;
                         ir   = _mm512_setzero_pd();
                         x0   = _mm512_sqrt_pd(_mm512_mul_pd(C314159265358979323846264338328,k0a));
                         coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         htht = _mm512_mul_pd(C05,tht);
                         ear1 = ir;
                         ii   = _mm512_add_pd(x0,x0);
                         k0a2 = _mm512_add_pd(k0a,k0a);
                         sint = xsin(tht);
                         x0   = _mm512_add_pd(C10,sint);
                         eai2 = _mm512_fmsub_pd(k0a2,x0,C078539816339744830961566084582);
                         cexp_zmm8r8(ear1,eai2,&cer1,&cei1);
                         arg1 = _mm512_sub_pd(C078539816339744830961566084582,htht);
                         carg1= xcos(arg1);
                         x1   = _mm512_sub_pd(C10,sint);
                         eai1 = _mm512_fmsub_pd(k0a2,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear1,eai1,&cer2,&cei2);
                         abs1 = _mm512_abs_pd(carg1);
                         t0r  = _mm512_div_pd(cer2,abs1);
                         arg2 = _mm512_add_pd(C078539816339744830961566084582,htht);
                         t0i  = _mm512_div_pd(cei2,abs1);
                         carg2= xcos(arg2);
                         abs2 = _mm512_abs_pd(carg2);
                         t1r  = _mm512_div_pd(cer1,abs2);
                         cmul_zmm8r8(ir,ii,t0r,t0i,&x2,&x3);
                         t1i  = _mm512_div_pd(cei1,abs2);
                         cmul_zmm8r8(ir,ii,t1r,t1i,&x4,&x5);
                         _mm512_store_pd(&B1r[0], _mm512_add_pd(A1r,x2));
                         _mm512_store_pd(&B2r[0], _mm512_add_pd(A2r,x4));
                         _mm512_store_pd(&B1i[0], _mm512_add_pd(A1i,x3));
                         _mm512_store_pd(&B2i[0], _mm512_add_pd(A2i,x5));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefB12_f7414_zmm8r8_u(const double * __restrict   pk0a,
                                                const double * __restrict   ptht,
                                                double * __restrict  B1r,
                                                double * __restrict  B1i,
                                                double * __restrict  B2r,
                                                double * __restrict  B2i) {

                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(-0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C05                              = _mm512_set1_pd(0.5f);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear1,eai1,eai2,cer1,cei1,cer2,cei2;
                         register __m512d ir,ii,sint,htht,carg1,carg2,arg1,arg2,abs1,abs2;
                         register __m512d x0,x1,x2,x3,x4,x5,t0r,t0i,t1r,t1i,k0a2;
                         register __m512d A1r,A1i,A2r,A2i;
                         ir   = _mm512_setzero_pd();
                         x0   = _mm512_sqrt_pd(_mm512_mul_pd(C314159265358979323846264338328,k0a));
                         coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         htht = _mm512_mul_pd(C05,tht);
                         ear1 = ir;
                         ii   = _mm512_add_pd(x0,x0);
                         k0a2 = _mm512_add_pd(k0a,k0a);
                         sint = xsin(tht);
                         x0   = _mm512_add_pd(C10,sint);
                         eai2 = _mm512_fmsub_pd(k0a2,x0,C078539816339744830961566084582);
                         cexp_zmm8r8(ear1,eai2,&cer1,&cei1);
                         arg1 = _mm512_sub_pd(C078539816339744830961566084582,htht);
                         carg1= xcos(arg1);
                         x1   = _mm512_sub_pd(C10,sint);
                         eai1 = _mm512_fmsub_pd(k0a2,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear1,eai1,&cer2,&cei2);
                         abs1 = _mm512_abs_pd(carg1);
                         t0r  = _mm512_div_pd(cer2,abs1);
                         arg2 = _mm512_add_pd(C078539816339744830961566084582,htht);
                         t0i  = _mm512_div_pd(cei2,abs1);
                         carg2= xcos(arg2);
                         abs2 = _mm512_abs_pd(carg2);
                         t1r  = _mm512_div_pd(cer1,abs2);
                         cmul_zmm8r8(ir,ii,t0r,t0i,&x2,&x3);
                         t1i  = _mm512_div_pd(cei1,abs2);
                         cmul_zmm8r8(ir,ii,t1r,t1i,&x4,&x5);
                         _mm512_storeu_pd(&B1r[0], _mm512_add_pd(A1r,x2));
                         _mm512_storeu_pd(&B2r[0], _mm512_add_pd(A2r,x4));
                         _mm512_storeu_pd(&B1i[0], _mm512_add_pd(A1i,x3));
                         _mm512_storeu_pd(&B2i[0], _mm512_add_pd(A2i,x5));
               }


                  /*
                       Very Important!!
                       Backscattered fields from the edges of strips.
                       Ufimtsev derivation.
                       Electric-field (over z).
                       Formula 7.4-9
                 */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f749_zmm8r8(const __m512d tht,
                                         const __m512d k0a,
                                         const __m512d k0r,
                                         const __m512d Eir,
                                         const __m512d Eii,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {

                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C6283185307179586476925286766559 = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512d B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                         register __m512d ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         register __m512d sqr,x0,x1,k0a2,y0;
                         ear   = _mm512_setzero_pd();
                         sqr   = _mm512_sqrt_pd(_mm512_mul_pd(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_pd(k0a,k0a);
                         eai   = _mm512_add_pd(k0r,C078539816339744830961566084582);
                         sint  = xsin(tht);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_pd(sint,sint);
                         coefB12_f7414_zmm8r8(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                         sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                         cmul_zmm8r8(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                         ear2  = ear;
                         sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                         y0    = _mm512_mul_pd(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_pd(sin1,cer2);
                         cei2  = _mm512_mul_pd(sin1,cei2);
                         ear3  = ear;
                         cmul_zmm8r8(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                         eai3  = gms::math::negate_zmm8r8(y0);
                         cexp_zmm8r8(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_pd(sin2,cer3);
                         cei3  = _mm512_mul_pd(sin2,cei3);
                         cmul_zmm8r8(B2rs,B2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm8r8(B1rs,B1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_pd(t0r,t1r);
                         t2i = _mm512_sub_pd(t0i,t1i);
                         cmul_zmm8r8(Eir,Eii,t2r,t2i,&x0,&x1);
                         cmul_zmm8r8(x0,x1,cer,cei,*Esr,*Esi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f749_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                           const double * __restrict __ATTR_ALIGN__(64) pEir,
                                           const double * __restrict __ATTR_ALIGN__(64) pEii,
                                           double * __restrict __ATTR_ALIGN__(64)  Esr,
                                           double * __restrict __ATTR_ALIGN__(64)  Esi) {

                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         register __m512d k0r = _mm512_load_pd(&pk0r[0]);
                         register __m512d Eir = _mm512_load_pd(&pEir[0]);
                         register __m512d Eii = _mm512_load_pd(&pEii[0]);
                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C6283185307179586476925286766559 = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         __m512d ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512d B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                         register __m512d ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512d sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_pd();
                         sqr   = _mm512_sqrt_pd(_mm512_mul_pd(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_pd(k0a,k0a);
                         eai   = _mm512_add_pd(k0r,C078539816339744830961566084582);
                         sint  = xsin(tht);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_pd(sint,sint);
                         coefB12_f7414_zmm8r8(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                         sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                         cmul_zmm8r8(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                         ear2  = ear;
                         sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                         y0    = _mm512_mul_pd(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_pd(sin1,cer2);
                         cei2  = _mm512_mul_pd(sin1,cei2);
                         ear3  = ear;
                         cmul_zmm8r8(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                         eai3  = gms::math::negate_zmm8r8(y0);
                         cexp_zmm8r8(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_pd(sin2,cer3);
                         cei3  = _mm512_mul_pd(sin2,cei3);
                         cmul_zmm8r8(B2rs,B2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm8r8(B1rs,B1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_pd(t0r,t1r);
                         t2i = _mm512_sub_pd(t0i,t1i);
                         cmul_zmm8r8(Eir,Eii,t2r,t2i,&x0,&x1);
                         cmul_zmm8r8(x0,x1,cer,cei,&resr,&resi);
                         _mm512_store_pd(&Esr[0], resr);
                         _mm512_store_pd(&Esi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f749_zmm8r8_u(const double * __restrict  ptht,
                                           const double * __restrict  pk0a,
                                           const double * __restrict  pk0r,
                                           const double * __restrict  pEir,
                                           const double * __restrict  pEii,
                                           double * __restrict   Esr,
                                           double * __restrict   Esi) {

                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d k0r = _mm512_loadu_pd(&pk0r[0]);
                         register __m512d Eir = _mm512_loadu_pd(&pEir[0]);
                         register __m512d Eii = _mm512_loadu_pd(&pEii[0]);
                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C6283185307179586476925286766559 = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         __m512d ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512d B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                         register __m512d ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512d sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_pd();
                         sqr   = _mm512_sqrt_pd(_mm512_mul_pd(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_pd(k0a,k0a);
                         eai   = _mm512_add_pd(k0r,C078539816339744830961566084582);
                         sint  = xsin(tht);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_pd(sint,sint);
                         coefB12_f7414_zmm8r8(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                         sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                         cmul_zmm8r8(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                         ear2  = ear;
                         sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                         y0    = _mm512_mul_pd(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_pd(sin1,cer2);
                         cei2  = _mm512_mul_pd(sin1,cei2);
                         ear3  = ear;
                         cmul_zmm8r8(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                         eai3  = gms::math::negate_zmm8r8(y0);
                         cexp_zmm8r8(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_pd(sin2,cer3);
                         cei3  = _mm512_mul_pd(sin2,cei3);
                         cmul_zmm8r8(B2rs,B2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm8r8(B1rs,B1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_pd(t0r,t1r);
                         t2i = _mm512_sub_pd(t0i,t1i);
                         cmul_zmm8r8(Eir,Eii,t2r,t2i,&x0,&x1);
                         cmul_zmm8r8(x0,x1,cer,cei,&resr,&resi);
                         _mm512_storeu_pd(&Esr[0], resr);
                         _mm512_storeu_pd(&Esi[0], resi);
                 }


                    /*
                       Very Important!!
                       Backscattered fields from the edges of strips.
                       Ufimtsev derivation.
                       Electric-field (over z).
                       Formula 7.4-10
                 */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f749_zmm8r8(const __m512d tht,
                                         const __m512d k0a,
                                         const __m512d k0r,
                                         const __m512d Hir,
                                         const __m512d Hii,
                                         __m512d * __restrict Hsr,
                                         __m512d * __restrict Hsi) {

                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C6283185307179586476925286766559 = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512d A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                         register __m512d ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         register __m512d sqr,x0,x1,k0a2,y0;
                         ear   = _mm512_setzero_pd();
                         sqr   = _mm512_sqrt_pd(_mm512_mul_pd(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_pd(k0a,k0a);
                         eai   = _mm512_add_pd(k0r,C078539816339744830961566084582);
                         sint  = xsin(tht);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_pd(sint,sint);
                         coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                         cmul_zmm8r8(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                         ear2  = ear;
                         sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                         y0    = _mm512_mul_pd(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_pd(sin2,cer2);
                         cei2  = _mm512_mul_pd(sin2,cei2);
                         ear3  = ear;
                         cmul_zmm8r8(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                         eai3  = gms::math::negate_zmm8r8(y0);
                         cexp_zmm8r8(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_pd(sin1,cer3);
                         cei3  = _mm512_mul_pd(sin1,cei3);
                         cmul_zmm8r8(A2rs,A2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm8r8(A1rs,A1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_pd(t0r,t1r);
                         t2i = _mm512_sub_pd(t0i,t1i);
                         cmul_zmm8r8(Hir,Hii,t2r,t2i,&x0,&x1);
                         cmul_zmm8r8(x0,x1,cer,cei,*Hsr,*Hsi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f749_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                           const double * __restrict __ATTR_ALIGN__(64) pHir,
                                           const double * __restrict __ATTR_ALIGN__(64) pHii,
                                           double * __restrict __ATTR_ALIGN__(64)  Hsr,
                                           double * __restrict __ATTR_ALIGN__(64)  Hsi) {

                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         register __m512d k0r = _mm512_load_pd(&pk0r[0]);
                         register __m512d Eir = _mm512_load_pd(&pHir[0]);
                         register __m512d Eii = _mm512_load_pd(&pHii[0]);
                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C6283185307179586476925286766559 = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512d A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                         register __m512d ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512d sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_pd();
                         sqr   = _mm512_sqrt_pd(_mm512_mul_pd(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_pd(k0a,k0a);
                         eai   = _mm512_add_pd(k0r,C078539816339744830961566084582);
                         sint  = xsin(tht);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_pd(sint,sint);
                         coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                         cmul_zmm8r8(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                         ear2  = ear;
                         sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                         y0    = _mm512_mul_pd(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_pd(sin2,cer2);
                         cei2  = _mm512_mul_pd(sin2,cei2);
                         ear3  = ear;
                         cmul_zmm8r8(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                         eai3  = gms::math::negate_zmm8r8(y0);
                         cexp_zmm8r8(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_pd(sin1,cer3);
                         cei3  = _mm512_mul_pd(sin1,cei3);
                         cmul_zmm8r8(A2rs,A2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm8r8(A1rs,A1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_pd(t0r,t1r);
                         t2i = _mm512_sub_pd(t0i,t1i);
                         cmul_zmm8r8(Hir,Hii,t2r,t2i,&x0,&x1);
                         cmul_zmm8r8(x0,x1,cer,cei,&resr,*resi);
                         _mm512_store_pd(&Hsr[0], resr);
                         _mm512_store_pd(&Hsi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f749_zmm8r8_u(const double * __restrict  ptht,
                                           const double * __restrict  pk0a,
                                           const double * __restrict  pk0r,
                                           const double * __restrict  pHir,
                                           const double * __restrict  pHii,
                                           double * __restrict __ATTR_ALIGN__(64)  Hsr,
                                           double * __restrict __ATTR_ALIGN__(64)  Hsi) {

                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d k0r = _mm512_loadu_pd(&pk0r[0]);
                         register __m512d Eir = _mm512_loadu_pd(&pHir[0]);
                         register __m512d Eii = _mm512_loadu_pd(&pHii[0]);
                         const __m512d C078539816339744830961566084582  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                         const __m512d C6283185307179586476925286766559 = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d C10                              = _mm512_set1_pd(1.0f);
                         register __m512d ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512d A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                         register __m512d ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512d sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_pd();
                         sqr   = _mm512_sqrt_pd(_mm512_mul_pd(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_pd(k0a,k0a);
                         eai   = _mm512_add_pd(k0r,C078539816339744830961566084582);
                         sint  = xsin(tht);
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_pd(sint,sint);
                         coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                         cmul_zmm8r8(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                         ear2  = ear;
                         sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                         y0    = _mm512_mul_pd(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_pd(sin2,cer2);
                         cei2  = _mm512_mul_pd(sin2,cei2);
                         ear3  = ear;
                         cmul_zmm8r8(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                         eai3  = gms::math::negate_zmm8r8(y0);
                         cexp_zmm8r8(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_pd(sin1,cer3);
                         cei3  = _mm512_mul_pd(sin1,cei3);
                         cmul_zmm8r8(A2rs,A2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm8r8(A1rs,A1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_pd(t0r,t1r);
                         t2i = _mm512_sub_pd(t0i,t1i);
                         cmul_zmm8r8(Hir,Hii,t2r,t2i,&x0,&x1);
                         cmul_zmm8r8(x0,x1,cer,cei,&resr,*resi);
                         _mm512_storeu_pd(&Hsr[0], resr);
                         _mm512_storeu_pd(&Hsi[0], resi);
                 }


                       /*
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Perpendicular RCS.
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-11
                 */


                 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7411_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d tht) {

                          const __m512d C10   = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512d B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                          register __m512d ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512d t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          coefB12_f7414_zmm8r8(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                          sin2t = _mm512_add_pd(sint,sint);
                          ear1  = _mm512_setzero_pd();
                          x0    = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm8r8(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm8r8(x0);
                          cmul_zmm8r8(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                          sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                          cer1  = _mm512_mul_pd(sin1,cer1);
                          cei1  = _mm512_mul_pd(sin1,cei1);
                          cmul_zmm8r8(B2rs,B2is,cer1,cei1,&t0r,&t0i);
                          sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_pd(sin2,cer2);
                          cei2  = _mm512_mul_pd(sin2,cei2);
                          cmul_zmm8r8(B1rs,B1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_pd(t0r,t1r);
                          x1    = _mm512_sub_pd(t0i,t1i);
                          cabs  = cabs_zmm8r8(x0,x1);
                          rcs   = _mm512_mul_pd(ink0,cabs);
                          return (rcs);
                }


 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7411_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0 = _mm512_load_pd(&pk0[0]);
                          register __m512d a  = _mm512_load_pd(&pa[0]);
                          register __m512d tht= _mm512_load_pd(&ptht[0]);
                          const __m512d C10   = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512d B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                          register __m512d ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          __m512d t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          coefB12_f7414_zmm8r8(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                          sin2t = _mm512_add_pd(sint,sint);
                          ear1  = _mm512_setzero_pd();
                          x0    = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm8r8(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm8r8(x0);
                          cmul_zmm8r8(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                          sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                          cer1  = _mm512_mul_pd(sin1,cer1);
                          cei1  = _mm512_mul_pd(sin1,cei1);
                          cmul_zmm8r8(B2rs,B2is,cer1,cei1,&t0r,&t0i);
                          sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_pd(sin2,cer2);
                          cei2  = _mm512_mul_pd(sin2,cei2);
                          cmul_zmm8r8(B1rs,B1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_pd(t0r,t1r);
                          x1    = _mm512_sub_pd(t0i,t1i);
                          cabs  = cabs_zmm8r8(x0,x1);
                          rcs   = _mm512_mul_pd(ink0,cabs);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7411_zmm8r8_u(const double * __restrict pk0,
                                            const double * __restrict  pa,
                                            const double * __restrict  ptht) {

                          register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a  = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht= _mm512_loadu_pd(&ptht[0]);
                          const __m512d C10   = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512d B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                          register __m512d ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          __m512d t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          coefB12_f7414_zmm8r8(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                          sin2t = _mm512_add_pd(sint,sint);
                          ear1  = _mm512_setzero_pd();
                          x0    = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm8r8(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm8r8(x0);
                          cmul_zmm8r8(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                          sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                          cer1  = _mm512_mul_pd(sin1,cer1);
                          cei1  = _mm512_mul_pd(sin1,cei1);
                          cmul_zmm8r8(B2rs,B2is,cer1,cei1,&t0r,&t0i);
                          sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_pd(sin2,cer2);
                          cei2  = _mm512_mul_pd(sin2,cei2);
                          cmul_zmm8r8(B1rs,B1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_pd(t0r,t1r);
                          x1    = _mm512_sub_pd(t0i,t1i);
                          cabs  = cabs_zmm8r8(x0,x1);
                          rcs   = _mm512_mul_pd(ink0,cabs);
                          return (rcs);
                }


                   
                       /*
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Parallel RCS.
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-12
                 */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7412_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d tht) {

                          const __m512d C10   = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512d A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                          register __m512d ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512d t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                          sin2t = _mm512_add_pd(sint,sint);
                          ear1  = _mm512_setzero_pd();
                          x0    = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm8r8(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm8r8(x0);
                          cmul_zmm8r8(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                          sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                          sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                          cer1  = _mm512_mul_pd(sin2,cer1);
                          cei1  = _mm512_mul_pd(sin2,cei1);
                          cmul_zmm8r8(A2rs,A2is,cer1,cei1,&t0r,&t0i);
                          cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_pd(sin1,cer2);
                          cei2  = _mm512_mul_pd(sin1,cei2);
                          cmul_zmm8r8(A1rs,A1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_pd(t0r,t1r);
                          x1    = _mm512_sub_pd(t0i,t1i);
                          cabs  = cabs_zmm8r8(x0,x1);
                          rcs   = _mm512_mul_pd(ink0,cabs);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7412_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0 = _mm512_load_pd(&pk0[0]);
                          register __m512d a  = _mm512_load_pd(&pa[0]);
                          register __m512d tht= _mm512_load_pd(&ptht[0]);
                          const __m512d C10   = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512d A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                          register __m512d ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512d t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                          sin2t = _mm512_add_pd(sint,sint);
                          ear1  = _mm512_setzero_pd();
                          x0    = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm8r8(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm8r8(x0);
                          cmul_zmm8r8(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                          sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                          sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                          cer1  = _mm512_mul_pd(sin2,cer1);
                          cei1  = _mm512_mul_pd(sin2,cei1);
                          cmul_zmm8r8(A2rs,A2is,cer1,cei1,&t0r,&t0i);
                          cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_pd(sin1,cer2);
                          cei2  = _mm512_mul_pd(sin1,cei2);
                          cmul_zmm8r8(A1rs,A1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_pd(t0r,t1r);
                          x1    = _mm512_sub_pd(t0i,t1i);
                          cabs  = cabs_zmm8r8(x0,x1);
                          rcs   = _mm512_mul_pd(ink0,cabs);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7412_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  ptht) {

                          register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a  = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht= _mm512_loadu_pd(&ptht[0]);
                          const __m512d C10   = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512d A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                          register __m512d ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512d t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          sint = xsin(tht);
                          coefA12_f7413_zmm8r8(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                          sin2t = _mm512_add_pd(sint,sint);
                          ear1  = _mm512_setzero_pd();
                          x0    = _mm512_mul_pd(_mm512_add_pd(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm8r8(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm8r8(x0);
                          cmul_zmm8r8(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                          sin1  = _mm512_div_pd(_mm512_sub_pd(C10,sint),sin2t);
                          cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                          sin2  = _mm512_div_pd(_mm512_add_pd(C10,sint),sin2t);
                          cer1  = _mm512_mul_pd(sin2,cer1);
                          cei1  = _mm512_mul_pd(sin2,cei1);
                          cmul_zmm8r8(A2rs,A2is,cer1,cei1,&t0r,&t0i);
                          cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_pd(sin1,cer2);
                          cei2  = _mm512_mul_pd(sin1,cei2);
                          cmul_zmm8r8(A1rs,A1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_pd(t0r,t1r);
                          x1    = _mm512_sub_pd(t0i,t1i);
                          cabs  = cabs_zmm8r8(x0,x1);
                          rcs   = _mm512_mul_pd(ink0,cabs);
                          return (rcs);
                }


                   /*
                          
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Incident angle very close to zero.
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-16, 7.4-17
              
                     */

                    
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7416_f7417_zmm8r8(const __m512d k0,
                                                  const __m512d a,
                                                  const __m512d tht) {

                          const __m512d C0001  = _mm512_set1_pd(0.0001); //very close value to zero incidence (rad)
                          const __m512d C00    = _mm512_setzero_pd();
                          const __m512d C10    = _mm512_set1_pd(1.0f);
                          const __mmask16 m1  = _mm512_cmp_pd_mask(tht,C0001,_CMP_LE_OQ);
                          const __mmask16 m2  = _mm512_cmp_pd_mask(tht,C00,_CMPEQ_OQ);
                          register __m512d k0a,ink0,rcs;
                          
                          if(m1) {
                              register __m512d k0a2,sint,arg,sarg,carg,sints,x0,x1,rat;
                              k0a  = _mm512_mul_pd(k0,a);
                              sint = xsin(tht);
                              k0a2 = _mm512_add_pd(k0a,k0a);
                              arg  = _mm512_mul_pd(k0a2,sint);
                              ink0 = _mm512_rcp14_pd(k0);
                              sarg = xsin(arg);
                              sints= _mm512_mul_pd(sint,sint);
                              x0   = _mm512_mul_pd(carg,carg);
                              x1   = _mm512_mul_pd(sarg,sarg);
                              rat  = _mm512_div_pd(x1,sints);
                              rcs  = _mm512_mul_pd(ink0, _mm512_add_pd(rat,carg));
                              return (rcs);
                          }
                          else if(m2) {
                              register __m512d x0,k0a2;
                              k0a  = _mm512_mul_pd(k0,a);
                              ink0 = _mm512_rcp14_pd(k0);
                              k0a2 = _mm512_add_pd(k0a,k0a);
                              x0   = _mm512_fmadd_pd(k0a2,k0a2,C10);
                              rcs  = _mm512_mul_pd(ink0,x0);
                              return (rcs);
                          }
                          
                  }
                  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7416_f7417_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                                    const double * __restrict __ATTR_ALIGN__(64) pa,
                                                    const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          const __m512d C0001  = _mm512_set1_pd(0.0001); //very close value to zero incidence (rad)
                          const __m512d C00    = _mm512_setzero_pd();
                          const __m512d C10    = _mm512_set1_pd(1.0f);
                          const __mmask16 m1  = _mm512_cmp_pd_mask(tht,C0001,_CMP_LE_OQ);
                          const __mmask16 m2  = _mm512_cmp_pd_mask(tht,C00,_CMPEQ_OQ);
                          register __m512d k0a,ink0,rcs;
                          
                          if(m1) {
                              register __m512d k0a2,sint,arg,sarg,carg,sints,x0,x1,rat;
                              k0a  = _mm512_mul_pd(k0,a);
                              sint = xsin(tht);
                              k0a2 = _mm512_add_pd(k0a,k0a);
                              arg  = _mm512_mul_pd(k0a2,sint);
                              ink0 = _mm512_rcp14_pd(k0);
                              sarg = xsin(arg);
                              sints= _mm512_mul_pd(sint,sint);
                              x0   = _mm512_mul_pd(carg,carg);
                              x1   = _mm512_mul_pd(sarg,sarg);
                              rat  = _mm512_div_pd(x1,sints);
                              rcs  = _mm512_mul_pd(ink0, _mm512_add_pd(rat,carg));
                              return (rcs);
                          }
                          else if(m2) {
                              register __m512d x0,k0a2;
                              k0a  = _mm512_mul_pd(k0,a);
                              ink0 = _mm512_rcp14_pd(k0);
                              k0a2 = _mm512_add_pd(k0a,k0a);
                              x0   = _mm512_fmadd_pd(k0a2,k0a2,C10);
                              rcs  = _mm512_mul_pd(ink0,x0);
                              return (rcs);
                          }
                         
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7416_f7417_zmm8r8_u(const double * __restrict  pk0,
                                                    const double * __restrict  pa,
                                                    const double * __restrict  ptht) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          const __m512d C0001  = _mm512_set1_pd(0.0001); //very close value to zero incidence (rad)
                          const __m512d C00    = _mm512_setzero_pd();
                          const __m512d C10    = _mm512_set1_pd(1.0f);
                          const __mmask16 m1  = _mm512_cmp_pd_mask(tht,C0001,_CMP_LE_OQ);
                          const __mmask16 m2  = _mm512_cmp_pd_mask(tht,C00,_CMPEQ_OQ);
                          register __m512d k0a,ink0,rcs;
                          
                          if(m1) {
                              register __m512d k0a2,sint,arg,sarg,carg,sints,x0,x1,rat;
                              k0a  = _mm512_mul_pd(k0,a);
                              sint = xsin(tht);
                              k0a2 = _mm512_add_pd(k0a,k0a);
                              arg  = _mm512_mul_pd(k0a2,sint);
                              ink0 = _mm512_rcp14_pd(k0);
                              sarg = xsin(arg);
                              sints= _mm512_mul_pd(sint,sint);
                              x0   = _mm512_mul_pd(carg,carg);
                              x1   = _mm512_mul_pd(sarg,sarg);
                              rat  = _mm512_div_pd(x1,sints);
                              rcs  = _mm512_mul_pd(ink0, _mm512_add_pd(rat,carg));
                              return (rcs);
                          }
                          else if(m2) {
                              register __m512d x0,k0a2;
                              k0a  = _mm512_mul_pd(k0,a);
                              ink0 = _mm512_rcp14_pd(k0);
                              k0a2 = _mm512_add_pd(k0a,k0a);
                              x0   = _mm512_fmadd_pd(k0a2,k0a2,C10);
                              rcs  = _mm512_mul_pd(ink0,x0);
                              return (rcs);
                          }
                         
                  }


                       /*
                          
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Incident angle at theta = PI/2
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-18 (perpendicula)
              
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7418_zmm8r8(const __m512d k0,
                                            const __m512d a) {
                                              
                          const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328f);
                          const __m512d C160                             = _mm512_set1_pd(16.0f);
                          const __m512d C80                              = _mm512_set1_pd(8.0f);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,k0a8,sin,k0a16,rat;
                          register __m512d rcs;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          k0a8 = _mm512_mul_pd(C80,k0a);
                          k0a16= _mm512_mul_pd(C160,k0a);
                          sin  = xsin(k0a8);
                          rat  = _mm512_div_pd(sin,k0a16);
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(C10,rat));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7418_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa) {
                             
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);                 
                          const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328f);
                          const __m512d C160                             = _mm512_set1_pd(16.0f);
                          const __m512d C80                              = _mm512_set1_pd(8.0f);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,k0a8,sin,k0a16,rat;
                          register __m512d rcs;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          k0a8 = _mm512_mul_pd(C80,k0a);
                          k0a16= _mm512_mul_pd(C160,k0a);
                          sin  = xsin(k0a8);
                          rat  = _mm512_div_pd(sin,k0a16);
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(C10,rat));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7418_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa) {
                             
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);                 
                          const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328f);
                          const __m512d C160                             = _mm512_set1_pd(16.0f);
                          const __m512d C80                              = _mm512_set1_pd(8.0f);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          register __m512d ink0,k0a,k0a8,sin,k0a16,rat;
                          register __m512d rcs;
                          k0a  = _mm512_mul_pd(k0,a);
                          ink0 = _mm512_rcp14_pd(k0);
                          k0a8 = _mm512_mul_pd(C80,k0a);
                          k0a16= _mm512_mul_pd(C160,k0a);
                          sin  = xsin(k0a8);
                          rat  = _mm512_div_pd(sin,k0a16);
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(C10,rat));
                          return (rcs);
                 }


                    /*
                          
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Incident angle at theta = PI/2
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-19 (parallel)
              
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7419_zmm8r8() { return _mm512_setzero_pd();}


                    /*
                          
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Incident angle at theta near PI/2
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-20 (perpendicular)
              
                     */

#include <limits>
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7420_zmm8r8(const __m512d k0,
                                            const __m512d tht) {

                          const __m512d C0001                            = _mm512_set1_pd(0.0001); 
                          const __m512d C157079632679489661923132169164  = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          const __m512d diff  = _mm512_sub_pd(C157079632679489661923132169164,
                                                             _mm512_abs_pd(tht));
                          const __mmask16 m1 = _mm512_cmp_pd_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                             register __m512d ink0,sint,sin2t,abs,sqr,x0;
                             register __m512d rcs;
                             ink0   = _mm512_rcp14_pd(k0);
                             sint   = xsin(tht);
                             sint2t = _mm512_mul_pd(sint,sint);
                             x0     = _mm512_add_pd(C10,_mm512_abs_pd(sint));
                             sqr    = _mm512_div_pd(x0,sin2t);
                             rcs    = _mm512_mul_pd(ink0,_mm512_mul_pd(sqr,sqr));
                             return (rcs);
                         }
                          else {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                         }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7420_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0 = _mm512_load_pd(&pk0[0]);
                          register __m512d tht= _mm512_load_pd(&ptht[0]); 
                          const __m512d C0001                            = _mm512_set1_pd(0.0001); 
                          const __m512d C157079632679489661923132169164  = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          const __m512d diff  = _mm512_sub_pd(C157079632679489661923132169164,
                                                             _mm512_abs_pd(tht));
                          const __mmask16 m1 = _mm512_cmp_pd_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                             register __m512d ink0,sint,sin2t,abs,sqr,x0;
                             register __m512d rcs;
                             ink0   = _mm512_rcp14_pd(k0);
                             sint   = xsin(tht);
                             sint2t = _mm512_mul_pd(sint,sint);
                             x0     = _mm512_add_pd(C10,_mm512_abs_pd(sint));
                             sqr    = _mm512_div_pd(x0,sin2t);
                             rcs    = _mm512_mul_pd(ink0,_mm512_mul_pd(sqr,sqr));
                             return (rcs);
                         }
                          else {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                         }
                  }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7420_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  ptht) {

                          register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                          register __m512d tht= _mm512_loadu_pd(&ptht[0]); 
                          const __m512d C0001                            = _mm512_set1_pd(0.0001); 
                          const __m512d C157079632679489661923132169164  = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          const __m512d diff  = _mm512_sub_pd(C157079632679489661923132169164,
                                                             _mm512_abs_pd(tht));
                          const __mmask16 m1 = _mm512_cmp_pd_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                             register __m512d ink0,sint,sin2t,abs,sqr,x0;
                             register __m512d rcs;
                             ink0   = _mm512_rcp14_pd(k0);
                             sint   = xsin(tht);
                             sint2t = _mm512_mul_pd(sint,sint);
                             x0     = _mm512_add_pd(C10,_mm512_abs_pd(sint));
                             sqr    = _mm512_div_pd(x0,sin2t);
                             rcs    = _mm512_mul_pd(ink0,_mm512_mul_pd(sqr,sqr));
                             return (rcs);
                         }
                          else {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                         }
                  }


                      /*
                          
                       Very Important!!
                       Resultant RCS of backscattered fields from the edges of strips.
                       Incident angle at theta near PI/2
                       See: formula 7.4-10, 7.4-9   
                       Formula: 7.4-21 (parallel)
              
                     */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7421_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d tht) {

                          const __m512d C0001                            = _mm512_set1_pd(0.0001); 
                          const __m512d C0318309886183790671537767526745 = _mm512_set1_pd(0.318309886183790671537767526745);
                          const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          const __m512d diff  = _mm512_sub_pd(C157079632679489661923132169164,
                                                             _mm512_abs_pd(tht));
                          const __mmask16 m1 = _mm512_cmp_pd_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                              register __m512d ink0,k0a,k0a2,sint,sin2t,x0,x1;
                              register __m512d rcs,;
                              k0a  = _mm512_mul_pd(k0,a);
                              sint = xsin(tht);
                              k0a2 = _mm512_mul_pd(k0a,k0a);
                              sin2t= xsin(sint,sint);
                              ink0 = _mm512_rcp14_pd(k0);
                              x0   = _mm512_mul_pd(_mm512_mul_pd(k0a2,k0a2),
                                                          C0318309886183790671537767526745);
                              x1   = _mm512_div_pd(_mm512_sub_pd(C10,sin2t),sint);
                              x0   = _mm512_mul_pd(x0,x1);
                              rcs  = _mm512_mul_pd(ink0,_mm512_mul_pd(x0,x0));
                              return (rcs);
                       }
                       else {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                      }

                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7421_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          const __m512d C0001                            = _mm512_set1_pd(0.0001f); 
                          const __m512d C0318309886183790671537767526745 = _mm512_set1_pd(0.318309886183790671537767526745);
                          const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          const __m512d diff  = _mm512_sub_pd(C157079632679489661923132169164,
                                                             _mm512_abs_pd(tht));
                          const __mmask16 m1 = _mm512_cmp_pd_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                              register __m512d ink0,k0a,k0a2,sint,sin2t,x0,x1;
                              register __m512d rcs,;
                              k0a  = _mm512_mul_pd(k0,a);
                              sint = xsin(tht);
                              k0a2 = _mm512_mul_pd(k0a,k0a);
                              sin2t= xsin(sint,sint);
                              ink0 = _mm512_rcp14_pd(k0);
                              x0   = _mm512_mul_pd(_mm512_mul_pd(k0a2,k0a2),
                                                          C0318309886183790671537767526745);
                              x1   = _mm512_div_pd(_mm512_sub_pd(C10,sin2t),sint);
                              x0   = _mm512_mul_pd(x0,x1);
                              rcs  = _mm512_mul_pd(ink0,_mm512_mul_pd(x0,x0));
                              return (rcs);
                       }
                       else {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                      }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7421_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  ptht) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          const __m512d C0001                            = _mm512_set1_pd(0.0001f); 
                          const __m512d C0318309886183790671537767526745 = _mm512_set1_pd(0.318309886183790671537767526745);
                          const __m512d C314159265358979323846264338328  = _mm512_set1_pd(3.14159265358979323846264338328);
                          const __m512d C10                              = _mm512_set1_pd(1.0f);
                          const __m512d diff  = _mm512_sub_pd(C157079632679489661923132169164,
                                                             _mm512_abs_pd(tht));
                          const __mmask16 m1 = _mm512_cmp_pd_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                              register __m512d ink0,k0a,k0a2,sint,sin2t,x0,x1;
                              register __m512d rcs,;
                              k0a  = _mm512_mul_pd(k0,a);
                              sint = xsin(tht);
                              k0a2 = _mm512_mul_pd(k0a,k0a);
                              sin2t= xsin(sint,sint);
                              ink0 = _mm512_rcp14_pd(k0);
                              x0   = _mm512_mul_pd(_mm512_mul_pd(k0a2,k0a2),
                                                          C0318309886183790671537767526745);
                              x1   = _mm512_div_pd(_mm512_sub_pd(C10,sin2t),sint);
                              x0   = _mm512_mul_pd(x0,x1);
                              rcs  = _mm512_mul_pd(ink0,_mm512_mul_pd(x0,x0));
                              return (rcs);
                       }
                       else {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                      }

                 }


                   /*
                          Bistatic RCS at high frequencies.
                          Approximation by Sommerfeld-MacDonald technique.
                          Case of parallel RCS.
                          Formula 7.4-22

                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_zmm8r8(const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d k0,
                                          const __m512d a,
                                          __m512d * __restrict rcs1,
                                          __m512d * __restrict rcs2) {

                        const __m512d C05 = _mm512_set1_pd(0.5f);
                        register __m512d A1r1,A1i1,A1r2,A1i2;
                        register __m512d A2r1,A2i1,A2r2,A2i2;
                        register __m512d ear1,eai1,ear2,ea2;
                        register __m512d cer1,cei1,cer2,cei2;
                        register __m512d t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512d ink0,sints,costm,sint1,sint2;
                        register __m512d cost12,rat,k0a,tp,tm;
                        register __m512d x0,x1,x2,x3;
                        k0a   = _mm512_mul_pd(k0,a);
                        sint1 = xsin(tht1);
                        ear1  = _mm512_setzero_pd();
                        tp    = _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                        coefA12_f7413_zmm8r8(k0a,tht1,&A1r1,&A1i1,&A2r1,&A2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                        sint2 = xsin(tht2);
                        sints = _mm512_add_pd(sint1,sint2);
                        eai1  = _mm512_mul_pd(k0a,sints);
                        cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_pd(k0);
                        eai2  = gms::math::negate_zmm8r8(eai1);
                        cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                        sintp = xsin(tp);
                        costm = xcos(tm);
                        x0    = _mm512_div_pd(_mm512_add_pd(sintp,costm),sints);
                        coefA12_f7413_zmm8r8(k0a,tht2,&A1r2,&A1i2,&A2r2,&A2i2);
                        cer1  = _mm512_mul_pd(cer1,x0);
                        cer2  = _mm512_mul_pd(cer2,x0);
                        cei1  = _mm512_mul_pd(cei1,x0);
                        cei2  = _mm512_mul_pd(cei2,x0);
                        cmul_zmm8r8(A2r1,A2i1,A2r2,A2i2,&t0r,&t0i);
                        cmul_zmm8r8(A1r1,A1i1,A1r2,A1i2,&t1r,&t1i);
                        cmul_zmm8r8(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm8r8(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_pd(res1r,res2r);
                        ear2 = _mm512_add_pd(res1r,res2r);
                        eai1 = _mm512_sub_pd(resi1,res2i);
                        eai2 = _mm512_add_pd(res1i,res2i);
                        *rcs1= _mm512_mul_pd(ink0,cabs_zmm8r8(ear1,eai1));
                        *rcs2= _mm512_mul_pd(ink0,cabs_zmm8r8(ear2,eai2));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                            const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            double * __restrict __ATTR_ALIGN__(64) rcs1,
                                            double * __restrict __ATTR_ALIGN__(64) rcs2) {

                        register __m512d  tht1 = _mm512_load_pd(&ptht1[0]);
                        register __m512d  tht2 = _mm512_load_pd(&ptht2[0]);
                        register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d  a    = _mm512_load_pd(&pa[0]);
                        const __m512d C05 = _mm512_set1_pd(0.5f);
                        register __m512d A1r1,A1i1,A1r2,A1i2;
                        register __m512d A2r1,A2i1,A2r2,A2i2;
                        register __m512d ear1,eai1,ear2,ea2;
                        register __m512d cer1,cei1,cer2,cei2;
                        register __m512d t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512d ink0,sints,costm,sint1,sint2;
                        register __m512d cost12,rat,k0a,tp,tm;
                        __m512d x0,x1,x2,x3;
                        k0a   = _mm512_mul_pd(k0,a);
                        sint1 = xsin(tht1);
                        ear1  = _mm512_setzero_pd();
                        tp    = _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                        coefA12_f7413_zmm8r8(k0a,tht1,&A1r1,&A1i1,&A2r1,&A2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                        sint2 = xsin(tht2);
                        sints = _mm512_add_pd(sint1,sint2);
                        eai1  = _mm512_mul_pd(k0a,sints);
                        cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_pd(k0);
                        eai2  = gms::math::negate_zmm8r8(eai1);
                        cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                        sintp = xsin(tp);
                        costm = xcos(tm);
                        x0    = _mm512_div_pd(_mm512_add_pd(sintp,costm),sints);
                        coefA12_f7413_zmm8r8(k0a,tht2,&A1r2,&A1i2,&A2r2,&A2i2);
                        cer1  = _mm512_mul_pd(cer1,x0);
                        cer2  = _mm512_mul_pd(cer2,x0);
                        cei1  = _mm512_mul_pd(cei1,x0);
                        cei2  = _mm512_mul_pd(cei2,x0);
                        cmul_zmm8r8(A2r1,A2i1,A2r2,A2i2,&t0r,&t0i);
                        cmul_zmm8r8(A1r1,A1i1,A1r2,A1i2,&t1r,&t1i);
                        cmul_zmm8r8(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm8r8(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_pd(res1r,res2r);
                        ear2 = _mm512_add_pd(res1r,res2r);
                        eai1 = _mm512_sub_pd(resi1,res2i);
                        eai2 = _mm512_add_pd(res1i,res2i);
                        _mm512_store_pd(&rcs1[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear1,eai1)));
                        _mm512_store_pd(&rcs2[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear2,eai2)));
                }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_zmm8r8_u(const double * __restrict  ptht1,
                                            const double * __restrict  ptht2,
                                            const double * __restrict  pk0,
                                            const double * __restrict pa,
                                            double * __restrict  rcs1,
                                            double * __restrict  rcs2) {

                        register __m512d  tht1 = _mm512_loadu_pd(&ptht1[0]);
                        register __m512d  tht2 = _mm512_loadu_pd(&ptht2[0]);
                        register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d  a    = _mm512_loadu_pd(&pa[0]);
                        const __m512d C05 = _mm512_set1_pd(0.5f);
                        register __m512d A1r1,A1i1,A1r2,A1i2;
                        register __m512d A2r1,A2i1,A2r2,A2i2;
                        register __m512d ear1,eai1,ear2,ea2;
                        register __m512d cer1,cei1,cer2,cei2;
                        register __m512d t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512d ink0,sints,costm,sint1,sint2;
                        register __m512d cost12,rat,k0a,tp,tm;
                        __m512d x0,x1,x2,x3;
                        k0a   = _mm512_mul_pd(k0,a);
                        sint1 = xsin(tht1);
                        ear1  = _mm512_setzero_pd();
                        tp    = _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                        coefA12_f7413_zmm8r8(k0a,tht1,&A1r1,&A1i1,&A2r1,&A2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                        sint2 = xsin(tht2);
                        sints = _mm512_add_pd(sint1,sint2);
                        eai1  = _mm512_mul_pd(k0a,sints);
                        cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_pd(k0);
                        eai2  = gms::math::negate_zmm8r8(eai1);
                        cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                        sintp = xsin(tp);
                        costm = xcos(tm);
                        x0    = _mm512_div_pd(_mm512_add_pd(sintp,costm),sints);
                        coefA12_f7413_zmm8r8(k0a,tht2,&A1r2,&A1i2,&A2r2,&A2i2);
                        cer1  = _mm512_mul_pd(cer1,x0);
                        cer2  = _mm512_mul_pd(cer2,x0);
                        cei1  = _mm512_mul_pd(cei1,x0);
                        cei2  = _mm512_mul_pd(cei2,x0);
                        cmul_zmm8r8(A2r1,A2i1,A2r2,A2i2,&t0r,&t0i);
                        cmul_zmm8r8(A1r1,A1i1,A1r2,A1i2,&t1r,&t1i);
                        cmul_zmm8r8(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm8r8(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_pd(res1r,res2r);
                        ear2 = _mm512_add_pd(res1r,res2r);
                        eai1 = _mm512_sub_pd(resi1,res2i);
                        eai2 = _mm512_add_pd(res1i,res2i);
                        _mm512_storeu_pd(&rcs1[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear1,eai1)));
                        _mm512_storeu_pd(&rcs2[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear2,eai2)));
                }



                  
                   /*
                          Bistatic RCS at high frequencies.
                          Approximation by Sommerfeld-MacDonald technique.
                          Case of perpendicular RCS.
                          Formula 7.4-23

                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7423_zmm8r8(const __m512d tht1,
                                          const __m512d tht2,
                                          const __m512d k0,
                                          const __m512d a,
                                          __m512d * __restrict rcs1,
                                          __m512d * __restrict rcs2) {

                        const __m512d C05 = _mm512_set1_pd(0.5f);
                        register __m512d B1r1,B1i1,B1r2,B1i2;
                        register __m512d B2r1,B2i1,B2r2,B2i2;
                        register __m512d ear1,eai1,ear2,ea2;
                        register __m512d cer1,cei1,cer2,cei2;
                        register __m512d t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512d ink0,sints,costm,sint1,sint2;
                        register __m512d cost12,rat,k0a,tp,tm;
                        register __m512d x0,x1,x2,x3;
                        k0a   = _mm512_mul_pd(k0,a);
                        sint1 = xsin(tht1);
                        ear1  = _mm512_setzero_pd();
                        tp    = _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                        coefB12_f7414_zmm8r8(k0a,tht1,&B1r1,&B1i1,&B2r1,&B2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                        sint2 = xsin(tht2);
                        sints = _mm512_add_pd(sint1,sint2);
                        eai1  = _mm512_mul_pd(k0a,sints);
                        cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_pd(k0);
                        eai2  = gms::math::negate_zmm8r8(eai1);
                        cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                        sintp = xsin(tp);
                        costm = xcos(tm);
                        x0    = _mm512_div_pd(_mm512_add_pd(sintp,costm),sints);
                        coefB12_f7414_zmm8r8(k0a,tht2,&B1r2,&B1i2,&B2r2,&B2i2);
                        cer1  = _mm512_mul_pd(cer1,x0);
                        cer2  = _mm512_mul_pd(cer2,x0);
                        cei1  = _mm512_mul_pd(cei1,x0);
                        cei2  = _mm512_mul_pd(cei2,x0);
                        cmul_zmm8r8(B2r1,B2i1,B2r2,B2i2,&t0r,&t0i);
                        cmul_zmm8r8(B1r1,B1i1,B1r2,B1i2,&t1r,&t1i);
                        cmul_zmm8r8(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm8r8(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_pd(res1r,res2r);
                        ear2 = _mm512_add_pd(res1r,res2r);
                        eai1 = _mm512_sub_pd(resi1,res2i);
                        eai2 = _mm512_add_pd(res1i,res2i);
                        *rcs1= _mm512_mul_pd(ink0,cabs_zmm8r8(ear1,eai1));
                        *rcs2= _mm512_mul_pd(ink0,cabs_zmm8r8(ear2,eai2));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7423_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                            const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            double * __restrict __ATTR_ALIGN__(64) rcs1,
                                            double * __restrict __ATTR_ALIGN__(64) rcs2) {

                        register __m512d  tht1 = _mm512_load_pd(&ptht1[0]);
                        register __m512d  tht2 = _mm512_load_pd(&ptht2[0]);
                        register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d  a    = _mm512_load_pd(&pa[0]);

                        const __m512d C05 = _mm512_set1_pd(0.5f);
                        register __m512d B1r1,B1i1,B1r2,B1i2;
                        register __m512d B2r1,B2i1,B2r2,B2i2;
                        register __m512d ear1,eai1,ear2,ea2;
                        register __m512d cer1,cei1,cer2,cei2;
                        register __m512d t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512d ink0,sints,costm,sint1,sint2;
                        register __m512d cost12,rat,k0a,tp,tm;
                        register __m512d x0,x1,x2,x3;
                        k0a   = _mm512_mul_pd(k0,a);
                        sint1 = xsin(tht1);
                        ear1  = _mm512_setzero_pd();
                        tp    = _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                        coefB12_f7414_zmm8r8(k0a,tht1,&B1r1,&B1i1,&B2r1,&B2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                        sint2 = xsin(tht2);
                        sints = _mm512_add_pd(sint1,sint2);
                        eai1  = _mm512_mul_pd(k0a,sints);
                        cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_pd(k0);
                        eai2  = gms::math::negate_zmm8r8(eai1);
                        cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                        sintp = xsin(tp);
                        costm = xcos(tm);
                        x0    = _mm512_div_pd(_mm512_add_pd(sintp,costm),sints);
                        coefB12_f7414_zmm8r8(k0a,tht2,&B1r2,&B1i2,&B2r2,&B2i2);
                        cer1  = _mm512_mul_pd(cer1,x0);
                        cer2  = _mm512_mul_pd(cer2,x0);
                        cei1  = _mm512_mul_pd(cei1,x0);
                        cei2  = _mm512_mul_pd(cei2,x0);
                        cmul_zmm8r8(B2r1,B2i1,B2r2,B2i2,&t0r,&t0i);
                        cmul_zmm8r8(B1r1,B1i1,B1r2,B1i2,&t1r,&t1i);
                        cmul_zmm8r8(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm8r8(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_pd(res1r,res2r);
                        ear2 = _mm512_add_pd(res1r,res2r);
                        eai1 = _mm512_sub_pd(resi1,res2i);
                        eai2 = _mm512_add_pd(res1i,res2i);
                        _mm512_store_pd(&rcs1[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear1,eai1)));
                        _mm512_store_pd(&rcs2[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear2,eai2)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7423_zmm8r8_u(const double * __restrict  ptht1,
                                            const double * __restrict  ptht2,
                                            const double * __restrict  pk0,
                                            const double * __restrict  pa,
                                            double * __restrict  rcs1,
                                            double * __restrict  rcs2) {

                        register __m512d  tht1 = _mm512_loadu_pd(&ptht1[0]);
                        register __m512d  tht2 = _mm512_loadu_pd(&ptht2[0]);
                        register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d  a    = _mm512_loadu_pd(&pa[0]);

                        const __m512d C05 = _mm512_set1_pd(0.5f);
                        register __m512d B1r1,B1i1,B1r2,B1i2;
                        register __m512d B2r1,B2i1,B2r2,B2i2;
                        register __m512d ear1,eai1,ear2,ea2;
                        register __m512d cer1,cei1,cer2,cei2;
                        register __m512d t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512d ink0,sints,costm,sint1,sint2;
                        register __m512d cost12,rat,k0a,tp,tm;
                        register __m512d x0,x1,x2,x3;
                        k0a   = _mm512_mul_pd(k0,a);
                        sint1 = xsin(tht1);
                        ear1  = _mm512_setzero_pd();
                        tp    = _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                        coefB12_f7414_zmm8r8(k0a,tht1,&B1r1,&B1i1,&B2r1,&B2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                        sint2 = xsin(tht2);
                        sints = _mm512_add_pd(sint1,sint2);
                        eai1  = _mm512_mul_pd(k0a,sints);
                        cexp_zmm8r8(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_pd(k0);
                        eai2  = gms::math::negate_zmm8r8(eai1);
                        cexp_zmm8r8(ear2,eai2,&cer2,&cei2);
                        sintp = xsin(tp);
                        costm = xcos(tm);
                        x0    = _mm512_div_pd(_mm512_add_pd(sintp,costm),sints);
                        coefB12_f7414_zmm8r8(k0a,tht2,&B1r2,&B1i2,&B2r2,&B2i2);
                        cer1  = _mm512_mul_pd(cer1,x0);
                        cer2  = _mm512_mul_pd(cer2,x0);
                        cei1  = _mm512_mul_pd(cei1,x0);
                        cei2  = _mm512_mul_pd(cei2,x0);
                        cmul_zmm8r8(B2r1,B2i1,B2r2,B2i2,&t0r,&t0i);
                        cmul_zmm8r8(B1r1,B1i1,B1r2,B1i2,&t1r,&t1i);
                        cmul_zmm8r8(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm8r8(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_pd(res1r,res2r);
                        ear2 = _mm512_add_pd(res1r,res2r);
                        eai1 = _mm512_sub_pd(resi1,res2i);
                        eai2 = _mm512_add_pd(res1i,res2i);
                        _mm512_storeu_pd(&rcs1[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear1,eai1)));
                        _mm512_storeu_pd(&rcs2[0], _mm512_mul_pd(ink0,cabs_zmm8r8(ear2,eai2)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_f7423_zmm8r8(const __m512d tht1,
                                                const __m512d tht2,
                                                const __m512d k0,
                                                const __m512d a,
                                                __m512d * __restrict rcsx1,
                                                __m512d * __restrict rcsx2,
                                                __m512d * __restrict rcsy1,
                                                __m512d * __restrict rcsy2) {

                        rcs_f7422_zmm8r8(tht1,tht2,k0,a,*rcsx1,*rcsx2);
                        rcs_f7423_zmm8r8(tht1,tht2,k0,a,*rcsy1,*rcsy2);
               }


                    /*
                            RCS forward direction scattering.
                            Perpendicular and parallel.
                            Formula 7.4-24
                       */

                   /*
                                 Argument checking removed!!
                         */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7424_zmm8r8(const __m512d tht1,
                                            const __m512d tht2,
                                            const __m512d k0,
                                            const __m512d a) {

                          const __m512d C05 = _mm512_set1_pd(0.5f);
                          register __m512d thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512d rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_pd(k0,a);
                          thtp =  _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                          sint1=  xsin(tht1);
                          ink0 =  _mm512_rcp14_pd(k0);
                          sint2=  xsin(tht2);
                          arg  = _mm512_mul_pd(k0a,_mm512_add_pd(sint1,sint2));
                          sintp= xsin(thtp);
                          thtm =  _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                          sarg =  xsin(arg);
                          sqr1 = _mm512_mul_pd(sarg,sarg);
                          costm=  xcos(thtm);
                          thtp = _mm512_div_pd(sqr1,_mm512_mul_pd(costm,costm));
                          carg =  xcos(arg);
                          sqr2 = _mm512_mul_pd(carg,carg);
                          thtm = _mm512_div_pd(sqr2,_mm512_mul_pd(sintp,sintp));
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(thtp,thtm));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7424_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512d tht1  = _mm512_load_pd(&ptht1[0]);
                          register __m512d tht2  = _mm512_load_pd(&ptht2[0]);
                          register __m512d k0    = _mm512_load_pd(&pk0[0]);
                          register __m512d a     = _mm512_load_pd(&pa[0]);
                          const __m512d C05 = _mm512_set1_pd(0.5f);
                          register __m512d thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512d rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_pd(k0,a);
                          thtp =  _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                          sint1=  xsin(tht1);
                          ink0 =  _mm512_rcp14_pd(k0);
                          sint2=  xsin(tht2);
                          arg  = _mm512_mul_pd(k0a,_mm512_add_pd(sint1,sint2));
                          sintp= xsin(thtp);
                          thtm =  _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                          sarg =  xsin(arg);
                          sqr1 = _mm512_mul_pd(sarg,sarg);
                          costm=  xcos(thtm);
                          thtp = _mm512_div_pd(sqr1,_mm512_mul_pd(costm,costm));
                          carg =  xcos(arg);
                          sqr2 = _mm512_mul_pd(carg,carg);
                          thtm = _mm512_div_pd(sqr2,_mm512_mul_pd(sintp,sintp));
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(thtp,thtm));
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7424_zmm8r8_u(const double * __restrict  ptht1,
                                              const double * __restrict  ptht2,
                                              const double * __restrict  pk0,
                                              const double * __restrict  pa) {

                          register __m512d tht1  = _mm512_loadu_pd(&ptht1[0]);
                          register __m512d tht2  = _mm512_loadu_pd(&ptht2[0]);
                          register __m512d k0    = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a     = _mm512_loadu_pd(&pa[0]);
                          const __m512d C05 = _mm512_set1_pd(0.5f);
                          register __m512d thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512d rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_pd(k0,a);
                          thtp =  _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                          sint1=  xsin(tht1);
                          ink0 =  _mm512_rcp14_pd(k0);
                          sint2=  xsin(tht2);
                          arg  = _mm512_mul_pd(k0a,_mm512_add_pd(sint1,sint2));
                          sintp= xsin(thtp);
                          thtm =  _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                          sarg =  xsin(arg);
                          sqr1 = _mm512_mul_pd(sarg,sarg);
                          costm=  xcos(thtm);
                          thtp = _mm512_div_pd(sqr1,_mm512_mul_pd(costm,costm));
                          carg =  xcos(arg);
                          sqr2 = _mm512_mul_pd(carg,carg);
                          thtm = _mm512_div_pd(sqr2,_mm512_mul_pd(sintp,sintp));
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(thtp,thtm));
                          return (rcs);
                 }


                   /*
                           Backward hemisphere RCS. for |theta2| < PI/2
                           theta1 != PI/2 (edges interaction neglected).
                           Formula 7.4-27
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7427_zmm8r8(const __m512d tht1,
                                            const __m512d tht2,
                                            const __m512d k0,
                                            const __m512d a) {

                          const __m512d C05 = _mm512_set1_pd(0.5f);
                          register __m512d thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512d rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_pd(k0,a);
                          thtp =  _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                          sint1=  xsin(tht1);
                          ink0 =  _mm512_rcp14_pd(k0);
                          sint2=  xsin(tht2);
                          arg  = _mm512_mul_pd(k0a,_mm512_add_pd(sint1,sint2));
                          sintp= xsin(thtp);
                          thtm =  _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                          sarg =  xsin(arg);
                          sqr1 = _mm512_mul_pd(sarg,sarg);
                          costm=  xcos(thtm);
                          thtp = _mm512_div_pd(sqr1,_mm512_mul_pd(sintp,sintp));
                          carg =  xcos(arg);
                          sqr2 = _mm512_mul_pd(carg,carg);
                          thtm = _mm512_div_pd(sqr2,_mm512_mul_pd(costm,costm));
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(thtp,thtm));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7427_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) ptht1,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht2,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512d tht1  = _mm512_load_pd(&ptht1[0]);
                          register __m512d tht2  = _mm512_load_pd(&ptht2[0]);
                          register __m512d k0    = _mm512_load_pd(&pk0[0]);
                          register __m512d a     = _mm512_load_pd(&pa[0]);

                          const __m512d C05 = _mm512_set1_pd(0.5f);
                          register __m512d thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512d rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_pd(k0,a);
                          thtp =  _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                          sint1=  xsin(tht1);
                          ink0 =  _mm512_rcp14_pd(k0);
                          sint2=  xsin(tht2);
                          arg  = _mm512_mul_pd(k0a,_mm512_add_pd(sint1,sint2));
                          sintp= xsin(thtp);
                          thtm =  _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                          sarg =  xsin(arg);
                          sqr1 = _mm512_mul_pd(sarg,sarg);
                          costm=  xcos(thtm);
                          thtp = _mm512_div_pd(sqr1,_mm512_mul_pd(sintp,sintp));
                          carg =  xcos(arg);
                          sqr2 = _mm512_mul_pd(carg,carg);
                          thtm = _mm512_div_pd(sqr2,_mm512_mul_pd(costm,costm));
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(thtp,thtm));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f7427_zmm8r8_u(  const double * __restrict  ptht1,
                                                const double * __restrict  ptht2,
                                                const double * __restrict  pk0,
                                                const double * __restrict  pa) {

                          register __m512d tht1  = _mm512_loadu_pd(&ptht1[0]);
                          register __m512d tht2  = _mm512_loadu_pd(&ptht2[0]);
                          register __m512d k0    = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a     = _mm512_loadu_pd(&pa[0]);

                          const __m512d C05 = _mm512_set1_pd(0.5f);
                          register __m512d thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512d rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_pd(k0,a);
                          thtp =  _mm512_mul_pd(_mm512_add_pd(tht1,tht2),C05);
                          sint1=  xsin(tht1);
                          ink0 =  _mm512_rcp14_pd(k0);
                          sint2=  xsin(tht2);
                          arg  = _mm512_mul_pd(k0a,_mm512_add_pd(sint1,sint2));
                          sintp= xsin(thtp);
                          thtm =  _mm512_mul_pd(_mm512_sub_pd(tht1,tht2),C05);
                          sarg =  xsin(arg);
                          sqr1 = _mm512_mul_pd(sarg,sarg);
                          costm=  xcos(thtm);
                          thtp = _mm512_div_pd(sqr1,_mm512_mul_pd(sintp,sintp));
                          carg =  xcos(arg);
                          sqr2 = _mm512_mul_pd(carg,carg);
                          thtm = _mm512_div_pd(sqr2,_mm512_mul_pd(costm,costm));
                          rcs  = _mm512_mul_pd(ink0,_mm512_add_pd(thtp,thtm));
                          return (rcs);
                 }


                     /*
                             Flat plates approximations.
                             Low-frequency approximations for k0a<<1,
                             plane wave unit amplitude xz-plane propagating
                             at angle 'theta'.
                             Perpendicular polarization.
                             Theta component.
                             Formula 7.5-1
                        */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f751_zmm8r8(const __m512d tht,
                                         const __m512d phi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d a,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {

                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_pd(k0,r);
                        sinp = xsin(phi);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(sinp,cost));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        *Esr = _mm512_mul_pd(x1,cer);
                        cei  = _mm512_mul_pd(x0,cei);
                        *Esi = _mm512_mul_pd(x1,cei);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f751_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         const double * __restrict __ATTR_ALIGN__(64) pr,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         double * __restrict __ATTR_ALIGN__(64) Esr,
                                         double * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d phi = _mm512_load_pd(&pphi[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d a   = _mm512_load_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_pd(k0,r);
                        sinp = xsin(phi);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(sinp,cost));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_store_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_store_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));
                 }



                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f751_zmm8r8_u(const double * __restrict ptht,
                                         const double * __restrict  pphi,
                                         const double * __restrict  pk0,
                                         const double * __restrict  pr,
                                         const double * __restrict  pa,
                                         double * __restrict  Esr,
                                         double * __restrict  Esi) {

                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d a   = _mm512_loadu_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_pd(k0,r);
                        sinp = xsin(phi);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(sinp,cost));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_storeu_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_storeu_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));
                 }


                   /*
                             Flat plates approximations.
                             Low-frequency approximations for k0a<<1,
                             plane wave unit amplitude xz-plane propagating
                             at angle 'theta'.
                             Perpendicular polarization.
                             Phi component.
                             Formula 7.5-1
                      */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f751_zmm8r8(const __m512d tht,
                                         const __m512d tht2
                                         const __m512d phi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d a,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {

                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        const __m512d C0212206590789193781025178351163 = 
                                              _mm512_set1_pd(0.212206590789193781025178351163);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cosp,sint,x0,x1,inr,sint2;
                        k0r  = _mm512_mul_pd(k0,r);
                        sint = xsin(tht);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cosp = xcos(phi);
                        k02  = _mm512_mul_pd(k0,k0);
                        sint2= xsin(tht2);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_fmadd_pd(_mm512_mul_pd(C0212206590789193781025178351163,
                                         sint),sint2,_mm512_mul_pd(C0424413181578387562050356702327,cosp));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        *Esr = _mm512_mul_pd(x1,cer);
                        cei  = _mm512_mul_pd(x0,cei);
                        *Esi = _mm512_mul_pd(x1,cei);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f751_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                         const double * __restrict __ATTR_ALIGN__(64) ptht2
                                         const double * __restrict __ATTR_ALIGN__(64) pphi,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         const double * __restrict __ATTR_ALIGN__(64) pr,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         double * __restrict __ATTR_ALIGN__(64) Esr,
                                         double * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d tht2= _mm512_load_pd(&ptht2[0]);
                        register __m512d phi = _mm512_load_pd(&pphi[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d a   = _mm512_load_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        const __m512d C0212206590789193781025178351163 = 
                                              _mm512_set1_pd(0.212206590789193781025178351163);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cosp,sint,x0,x1,inr,sint2;
                        k0r  = _mm512_mul_pd(k0,r);
                        sint = xsin(tht);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cosp = xcos(phi);
                        k02  = _mm512_mul_pd(k0,k0);
                        sint2= xsin(tht2);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_fmadd_pd(_mm512_mul_pd(C0212206590789193781025178351163,
                                         sint),sint2,_mm512_mul_pd(C0424413181578387562050356702327,cosp));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_store_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_store_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f751_zmm8r8_u(const double * __restrict  ptht,
                                         const double * __restrict  ptht2
                                         const double * __restrict  pphi,
                                         const double * __restrict  pk0,
                                         const double * __restrict  pr,
                                         const double * __restrict  pa,
                                         double * __restrict  Esr,
                                         double * __restrict  Esi) {

                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d tht2= _mm512_loadu_pd(&ptht2[0]);
                        register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d a   = _mm512_loadu_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        const __m512d C0212206590789193781025178351163 = 
                                              _mm512_set1_pd(0.212206590789193781025178351163);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cosp,sint,x0,x1,inr,sint2;
                        k0r  = _mm512_mul_pd(k0,r);
                        sint = xsin(tht);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cosp = xcos(phi);
                        k02  = _mm512_mul_pd(k0,k0);
                        sint2= xsin(tht2);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_fmadd_pd(_mm512_mul_pd(C0212206590789193781025178351163,
                                         sint),sint2,_mm512_mul_pd(C0424413181578387562050356702327,cosp));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_storeu_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_storeu_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));
                 }   
                 
                 
                    /*
                             Flat plates approximations.
                             Low-frequency approximations for k0a<<1,
                             plane wave unit amplitude xz-plane propagating
                             at angle 'theta'.
                             Parallel polarization.
                             Theta component.
                             Formula 7.5-2
                        */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f752_zmm8r8(const __m512d tht,
                                         const __m512d tht2,
                                         const __m512d phi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d a,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {

                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,cosp,cost2,x0,x1,inr,t0;
                        k0r  = _mm512_mul_pd(k0,r);
                        cosp = xcos(phi);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        cost2= xcos(tht2);
                        t0   = _mm512_mul_pd(cost,cost2);
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(cosp,t0));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        *Esr = _mm512_mul_pd(x1,cer);
                        cei  = _mm512_mul_pd(x0,cei);
                        *Esi = _mm512_mul_pd(x1,cei);
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f752_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                         const double * __restrict __ATTR_ALIGN__(64) ptht2
                                         const double * __restrict __ATTR_ALIGN__(64) pphi,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         const double * __restrict __ATTR_ALIGN__(64) pr,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         double * __restrict __ATTR_ALIGN__(64) Esr,
                                         double * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d tht2= _mm512_load_pd(&ptht2[0]);
                        register __m512d phi = _mm512_load_pd(&pphi[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d a   = _mm512_load_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,cosp,cost2,x0,x1,inr,t0;
                        k0r  = _mm512_mul_pd(k0,r);
                        cosp = xcos(phi);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        cost2= xcos(tht2);
                        t0   = _mm512_mul_pd(cost,cost2);
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(cosp,t0));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_store_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_store_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f752_zmm8r8_u(const double * __restrict ptht,
                                         const double * __restrict ptht2
                                         const double * __restrict  pphi,
                                         const double * __restrict  pk0,
                                         const double * __restrict  pr,
                                         const double * __restrict  pa,
                                         double * __restrict Esr,
                                         double * __restrict  Esi) {

                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d tht2= _mm512_loadu_pd(&ptht2[0]);
                        register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d a   = _mm512_loadu_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,cosp,cost2,x0,x1,inr,t0;
                        k0r  = _mm512_mul_pd(k0,r);
                        cosp = xcos(phi);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        cost2= xcos(tht2);
                        t0   = _mm512_mul_pd(cost,cost2);
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(cosp,t0));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_storeu_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_storeu_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));
                 }
                 
                 
                 /*
                             Flat plates approximations.
                             Low-frequency approximations for k0a<<1,
                             plane wave unit amplitude xz-plane propagating
                             at angle 'phi'.
                             Parallel polarization.
                             Theta component.
                             Formula 7.5-2
                        */
                        
                        
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f752_zmm8r8(const __m512d tht,
                                         const __m512d phi2,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d a,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {
                        
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_pd(k0,r);
                        sinp = xsin(phi2);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(sinp,cost));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        *Esr = _mm512_mul_pd(x1,cer);
                        cei  = _mm512_mul_pd(x0,cei);
                        *Esi = _mm512_mul_pd(x1,cei);       
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f752_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           double * __restrict __ATTR_ALIGN__(64) Esr,
                                           double * __restrict __ATTR_ALIGN__(64) Esi) {
                        
                        register __m512d tht  = _mm512_load_pd(&ptht[0]);
                        register __m512d phi2 = _mm512_load_pd(&pphi2[0]);
                        register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d r    = _mm512_load_pd(&pr[0]);
                        register __m512d a    = _mm512_load_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_pd(k0,r);
                        sinp = xsin(phi2);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(sinp,cost));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_store_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_store_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));       
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f752_zmm8r8_u(const double * __restrict  ptht,
                                           const double * __restrict  pphi2,
                                           const double * __restrict  pk0,
                                           const double * __restrict  pr,
                                           const double * __restrict pa,
                                           double * __restrict Esr,
                                           double * __restrict Esi) {
                        
                        
                        register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                        register __m512d phi2 = _mm512_loadu_pd(&pphi2[0]);
                        register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d r    = _mm512_loadu_pd(&pr[0]);
                        register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        const __m512d C0424413181578387562050356702327 = 
                                              _mm512_set1_pd(0.424413181578387562050356702327);
                        register __m512d ear,eai,cer,cei,k0r,k02,a3;
                        register __m512d cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_pd(k0,r);
                        sinp = xsin(phi2);
                        ear  = _mm512_setzero_pd();
                        eai  = k0r;
                        cost = xcos(tht);
                        k02  = _mm512_mul_pd(k0,k0);
                        inr  = _mm512_rcp14_pd(r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                        k02  = gms::math::negate_zmm8r8(k02);
                        x0   = _mm512_mul_pd(C0424413181578387562050356702327,
                                                         _mm512_mul_pd(sinp,cost));
                        x1   = _mm512_mul_pd(k02,a3);
                        cer  = _mm512_mul_pd(x0,cer);
                        _mm512_storeu_pd(&Esr[0] ,_mm512_mul_pd(x1,cer));
                        cei  = _mm512_mul_pd(x0,cei);
                        _mm512_storeu_pd(&Esi[0] ,_mm512_mul_pd(x1,cei));       
                 }
                 
                 
                 /*
                         Backscatter RCS (of 7.5-,7.5-2)
                         Formula: 7.5-3 (perpendicular)
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f753_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht) {
                  
                          const __m512d C0565884242104516749400475603102 = 
                                        _mm512_set1_pd(0.565884242104516749400475603102);
                          const __m512d C20 = _mm512_set1_pd(2.0f);
                          register __m512d k02,k04,a2,a6,sint,sqr,x0,x1;
                          register __m512d rcs;
                          k02   = _mm512_mul_pd(k0,k0);
                          a2    = _mm512_mul_pd(a,a);
                          sint  = xsin(tht);
                          k04   = _mm512_mul_pd(k02,k02);
                          x0    = _mm512_fmadd_pd(sint,sint,C20);
                          a6    = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          sqr   = _mm512_mul_pd(x0,x0);
                          x1    = _mm512_mul_pd(k04,a6);
                          rcs   = _mm512_mul_pd(C0565884242104516749400475603102,
                                                               _mm512_mul_pd(x1,sqr));
                          return (rcs);                  
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f753_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht) {
                  
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          const __m512d C0565884242104516749400475603102 = 
                                        _mm512_set1_pd(0.565884242104516749400475603102);
                          const __m512d C20 = _mm512_set1_pd(2.0f);
                          register __m512d k02,k04,a2,a6,sint,sqr,x0,x1;
                          register __m512d rcs;
                          k02   = _mm512_mul_pd(k0,k0);
                          a2    = _mm512_mul_pd(a,a);
                          sint  = xsin(tht);
                          k04   = _mm512_mul_pd(k02,k02);
                          x0    = _mm512_fmadd_pd(sint,sint,C20);
                          a6    = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          sqr   = _mm512_mul_pd(x0,x0);
                          x1    = _mm512_mul_pd(k04,a6);
                          rcs   = _mm512_mul_pd(C0565884242104516749400475603102,
                                                               _mm512_mul_pd(x1,sqr));
                          return (rcs);                  
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f753_zmm8r8_u(const double * __restrict pk0,
                                             const double * __restrict pa,
                                             const double * __restrict  ptht) {
                  
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          const __m512d C0565884242104516749400475603102 = 
                                        _mm512_set1_pd(0.565884242104516749400475603102);
                          const __m512d C20 = _mm512_set1_pd(2.0f);
                          register __m512d k02,k04,a2,a6,sint,sqr,x0,x1;
                          register __m512d rcs;
                          k02   = _mm512_mul_pd(k0,k0);
                          a2    = _mm512_mul_pd(a,a);
                          sint  = xsin(tht);
                          k04   = _mm512_mul_pd(k02,k02);
                          x0    = _mm512_fmadd_pd(sint,sint,C20);
                          a6    = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          sqr   = _mm512_mul_pd(x0,x0);
                          x1    = _mm512_mul_pd(k04,a6);
                          rcs   = _mm512_mul_pd(C0565884242104516749400475603102,
                                                               _mm512_mul_pd(x1,sqr));
                          return (rcs);                  
                 }
                 
                 
                   /*
                         Backscatter RCS (of 7.5-,7.5-2)
                         Formula: 7.5-4 (paralell)
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f754_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht) {
                                           
                          const __m512d C2263536968418066997601902412409 = 
                                          _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d k02,k04,a2,a6,cost,cost2t,cost4t,x0;
                          register __m512d rcs;
                          k02  = _mm512_mul_pd(k0,k0);
                          cost = xcos(tht);
                          a2   = _mm512_mul_pd(a,a);
                          cos2t= _mm512_mul_pd(cost,cost);
                          k04  = _mm512_mul_pd(k02,k02);
                          a6   = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));  
                          x0   = _mm512_mul_pd(k04,a6);
                          cos4t= _mm512_mul_pd(cos2t,cos2t);
                          rcs  = _mm512_mul_pd(C2263536968418066997601902412409,
                                                                 _mm512_mul_pd(x0,cos4t));
                          return (rcs);                   
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f754_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht) {
                  
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                                           
                          const __m512d C2263536968418066997601902412409 = 
                                          _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d k02,k04,a2,a6,cost,cost2t,cost4t,x0;
                          register __m512d rcs;
                          k02  = _mm512_mul_pd(k0,k0);
                          cost = xcos(tht);
                          a2   = _mm512_mul_pd(a,a);
                          cos2t= _mm512_mul_pd(cost,cost);
                          k04  = _mm512_mul_pd(k02,k02);
                          a6   = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));  
                          x0   = _mm512_mul_pd(k04,a6);
                          cos4t= _mm512_mul_pd(cos2t,cos2t);
                          rcs  = _mm512_mul_pd(C2263536968418066997601902412409,
                                                                 _mm512_mul_pd(x0,cos4t));
                          return (rcs);                   
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f754_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa,
                                             const double * __restrict  ptht) {
                  
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                                           
                          const __m512d C2263536968418066997601902412409 = 
                                          _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d k02,k04,a2,a6,cost,cost2t,cost4t,x0;
                          register __m512d rcs;
                          k02  = _mm512_mul_pd(k0,k0);
                          cost = xcos(tht);
                          a2   = _mm512_mul_pd(a,a);
                          cos2t= _mm512_mul_pd(cost,cost);
                          k04  = _mm512_mul_pd(k02,k02);
                          a6   = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));  
                          x0   = _mm512_mul_pd(k04,a6);
                          cos4t= _mm512_mul_pd(cos2t,cos2t);
                          rcs  = _mm512_mul_pd(C2263536968418066997601902412409,
                                                                 _mm512_mul_pd(x0,cos4t));
                          return (rcs);                   
               }
               
               
               /*
                    Normal incidence backscatter RCS.
                    Formula: 7.5-5
               */
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f755_zmm8r8(const __m512d k0,
                                           const __m512d a) {
                                           
                         const __m512d C2263536968418066997601902412409 = 
                                          _mm512_set1_pd(2.263536968418066997601902412409);
                         register __m512d rcs,k02,k04,a2,a6;
                         k02 = _mm512_mul_pd(k0,k0);
                         a2  = _mm512_mul_pd(a,a);
                         k04 = _mm512_mul_pd(k02,k02);
                         a6  = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                         rcs = _mm512_mul_pd( C2263536968418066997601902412409,
                                                                _mm512_mul_pd(k04,a6));
                         return (rcs);                         
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f755_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa) {
                                         
                         register __m512d k0  = _mm512_load_pd(&pk0[0]);
                         register __m512d a   = _mm512_load_pd(&pa[0]);
                               
                         const __m512d C2263536968418066997601902412409 = 
                                          _mm512_set1_pd(2.263536968418066997601902412409);
                         register __m512d rcs,k02,k04,a2,a6;
                         k02 = _mm512_mul_pd(k0,k0);
                         a2  = _mm512_mul_pd(a,a);
                         k04 = _mm512_mul_pd(k02,k02);
                         a6  = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                         rcs = _mm512_mul_pd( C2263536968418066997601902412409,
                                                                _mm512_mul_pd(k04,a6));
                         return (rcs);                         
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f755_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict pa) {
                                         
                         register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                         register __m512d a   = _mm512_loadu_pd(&pa[0]);
                               
                         const __m512d C2263536968418066997601902412409 = 
                                          _mm512_set1_pd(2.263536968418066997601902412409);
                         register __m512d rcs,k02,k04,a2,a6;
                         k02 = _mm512_mul_pd(k0,k0);
                         a2  = _mm512_mul_pd(a,a);
                         k04 = _mm512_mul_pd(k02,k02);
                         a6  = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                         rcs = _mm512_mul_pd( C2263536968418066997601902412409,
                                                                _mm512_mul_pd(k04,a6));
                         return (rcs);                         
                 }
                 
                 
                 /*
                       Wedge-on incidence.
                       RCS (perpendicular)
                       Formula: 7.5-6
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f756_zmm8r8(const __m512d k0,
                                           const __m512d a) {
                                           
                          const __m512d C509295817894065074460428042792 = 
                                                      _mm512_set1_pd(5.09295817894065074460428042792);
                          register __m512d rcs,k02,k04,a2,a6;
                          k02 = _mm512_mul_pd(k0,k0);
                          a2  = _mm512_mul_pd(a,a);
                          k04 = _mm512_mul_pd(k02,k02);
                          a6  = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          rcs = _mm512_mul_pd(C509295817894065074460428042792,
                                                               _mm512_mul_pd(k04,a6));
                          return (rcs);                        
                  }
                  
                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f756_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa) {
                                         
                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          const __m512d C509295817894065074460428042792 = 
                                                      _mm512_set1_pd(5.09295817894065074460428042792);
                          register __m512d rcs,k02,k04,a2,a6;
                          k02 = _mm512_mul_pd(k0,k0);
                          a2  = _mm512_mul_pd(a,a);
                          k04 = _mm512_mul_pd(k02,k02);
                          a6  = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          rcs = _mm512_mul_pd(C509295817894065074460428042792,
                                                               _mm512_mul_pd(k04,a6));
                          return (rcs);                        
                  }
                  

                    
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f756_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa) {
                                         
                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          const __m512d C509295817894065074460428042792 = 
                                                      _mm512_set1_pd(5.09295817894065074460428042792);
                          register __m512d rcs,k02,k04,a2,a6;
                          k02 = _mm512_mul_pd(k0,k0);
                          a2  = _mm512_mul_pd(a,a);
                          k04 = _mm512_mul_pd(k02,k02);
                          a6  = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          rcs = _mm512_mul_pd(C509295817894065074460428042792,
                                                               _mm512_mul_pd(k04,a6));
                          return (rcs);                        
                  }
                  
                  
                    /*
                       Wedge-on incidence.
                       RCS (parallel)
                       Formula: 7.5-7
                  */
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f757_zmm8r8() { return (_mm512_setzero_pd());}
                   
                   
                  /*
                          Low-frequency forward scatter RCS for theta2=PI-theta1
                          phi2=PI, perpendicular
                          Formula 7.5-8
                  */
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f758_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht2) {
                                           
                          return (rcs_f753_zmm8r8(k0,a,tht2));                         
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f758_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                             const double * __restrict __ATTR_ALIGN__(64)  pa,
                                             const double * __restrict __ATTR_ALIGN__(64)  ptht2) {
                                           
                          return (rcs_f753_zmm8r8_a(pk0,pa,ptht2));                         
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f758_zmm8r8_u(const double * __restrict  pk0,
                                             const double * __restrict  pa,
                                             const double * __restrict  ptht2) {
                                           
                          return (rcs_f753_zmm8r8_u(pk0,pa,ptht2));                         
                 }
                 
                 
                 
                   /*
                          Low-frequency forward scatter RCS for theta2=PI-theta1
                          phi2=PI, parallel
                          Formula 7.5-9
                  */ 
                  
                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f759_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d tht2) {
                                           
                          return (rcs_f754_zmm8r8(k0,a,tht2));                         
                 }
                 
                 
                 
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f759_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                             const double * __restrict __ATTR_ALIGN__(64)  pa,
                                             const double * __restrict __ATTR_ALIGN__(64)  ptht2) {
                                           
                          return (rcs_f754_zmm8r8_a(pk0,pa,ptht2));                         
                 }
                 
                 
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512d rcs_f759_zmm8r8_u(const double * __restrict   pk0,
                                             const double * __restrict   pa,
                                             const double * __restrict   ptht2) {
                                           
                          return (rcs_f754_zmm8r8_u(pk0,pa,ptht2));                         
                 }
                 
                 
                 /*
                       High frequency approximations.
                       Perfectly conducting disc.
                       Geometrical Diifraction backscatter RCS,
                       for 0<|theta|<PI.
                       !!NOT VALID FOR: theta==0, theta==PI/2!!
                       Formula: 7.5-13
                 */
                 
                 
                  // Input argument surpessed for this kernel.
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7513_zmm8r8(const __m512d k0a,
                                            const __m512d tht) {
                              
                                        
                         /* const __m512d C157079632679489661923132169164 = 
                                                _mm512_set1_pd(1.57079632679489661923132169164f);
                          const __m512d C00                             = _mm512_setzero_pd();
                          const __mmask16 m1 = _mm512_cmp_pd_mask(tht,C00,_CMP_EQ_OQ);
                          const __mmask16 m2 = _mm512_cmp_pd_mask(_mm512_abs_pd(tht),
                                                                      C157079632679489661923132169164,_CMP_EQ_OQ);
                          if(m1 || m2) {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                          }   */ 
                         const __m512d C39478417604357434475337963999505 = 
                                                              _mm512_set1_pd(39.478417604357434475337963999505);
                         register __m512d k0a2,sint,arg,carg,sarg,sin2t,x0,c2,s2;
                         register __m512d rcs;
                         k0a2 = _mm512_add_pd(k0a,k0a);
                         sint = xsin(tht);
                         arg  = _mm512_mul_pd(k0a2,sint);
                         x0   = _mm512_div_pd(k0a,
                                              _mm512_mul_pd(C39478417604357434475337963999505,sint));
                         sin2t= _mm512_mul_pd(sint,sint);
                         carg = xcos(arg);
                         sarg = xsin(arg);
                         c2   = _mm512_mul_pd(carg,carg);
                         s2   = _mm512_div_pd(_mm512_mul_pd(sarg,sarg),sint2t);
                         rcs  = _mm512_mul_pd(x0,_mm512_add_pd(c2,s2));
                         return (rcs);
                 }
                 
                 
                  
                 
                 
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7513_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht) {
                              
                                        
                         /* const __m512d C157079632679489661923132169164 = 
                                                _mm512_set1_pd(1.57079632679489661923132169164f);
                          const __m512d C00                             = _mm512_setzero_pd();
                          const __mmask16 m1 = _mm512_cmp_pd_mask(tht,C00,_CMP_EQ_OQ);
                          const __mmask16 m2 = _mm512_cmp_pd_mask(_mm512_abs_pd(tht),
                                                                      C157079632679489661923132169164,_CMP_EQ_OQ);
                          if(m1 || m2) {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                          }   */ 
                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         register __m512d tht = _mm512_load_pd(&ptht[0]);
                         const __m512d C39478417604357434475337963999505 = 
                                                              _mm512_set1_pd(39.478417604357434475337963999505);
                         register __m512d k0a2,sint,arg,carg,sarg,sin2t,x0,c2,s2;
                         register __m512d rcs;
                         k0a2 = _mm512_add_pd(k0a,k0a);
                         sint = xsin(tht);
                         arg  = _mm512_mul_pd(k0a2,sint);
                         x0   = _mm512_div_pd(k0a,
                                              _mm512_mul_pd(C39478417604357434475337963999505,sint));
                         sin2t= _mm512_mul_pd(sint,sint);
                         carg = xcos(arg);
                         sarg = xsin(arg);
                         c2   = _mm512_mul_pd(carg,carg);
                         s2   = _mm512_div_pd(_mm512_mul_pd(sarg,sarg),sint2t);
                         rcs  = _mm512_mul_pd(x0,_mm512_add_pd(c2,s2));
                         return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7513_zmm8r8_u(const double * __restrict  pk0a,
                                              const double * __restrict  ptht) {
                              
                                        
                         /* const __m512d C157079632679489661923132169164 = 
                                                _mm512_set1_pd(1.57079632679489661923132169164f);
                          const __m512d C00                             = _mm512_setzero_pd();
                          const __mmask16 m1 = _mm512_cmp_pd_mask(tht,C00,_CMP_EQ_OQ);
                          const __mmask16 m2 = _mm512_cmp_pd_mask(_mm512_abs_pd(tht),
                                                                      C157079632679489661923132169164,_CMP_EQ_OQ);
                          if(m1 || m2) {
                             __m512d NAN = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
                             return (NAN);
                          }   */ 
                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                         const __m512d C39478417604357434475337963999505 = 
                                                              _mm512_set1_pd(39.478417604357434475337963999505);
                         register __m512d k0a2,sint,arg,carg,sarg,sin2t,x0,c2,s2;
                         register __m512d rcs;
                         k0a2 = _mm512_add_pd(k0a,k0a);
                         sint = xsin(tht);
                         arg  = _mm512_mul_pd(k0a2,sint);
                         x0   = _mm512_div_pd(k0a,
                                              _mm512_mul_pd(C39478417604357434475337963999505,sint));
                         sin2t= _mm512_mul_pd(sint,sint);
                         carg = xcos(arg);
                         sarg = xsin(arg);
                         c2   = _mm512_mul_pd(carg,carg);
                         s2   = _mm512_div_pd(_mm512_mul_pd(sarg,sarg),sint2t);
                         rcs  = _mm512_mul_pd(x0,_mm512_add_pd(c2,s2));
                         return (rcs);
                 }
                 
                 
                 /*
                     Rectangular plates.
                     Backscatter RCS of perfectly conducting square.
                     Normal incidence.
                     Formula: 7.5-31
                 */
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7531_zmm8r8(const __m512d a, // side length of plate
                                            const __m512d gam0) {
                                            
                          const __m512d C1140 = _mm512_set1_pd(114.0f);
                          register __m512d a2,a6,gam2,gam4;
                          register __m512d rcs;
                          a2   = _mm512_mul_pd(a,a);
                          gam2 = _mm512_mul_pd(gam0,gam0);
                          a6   = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          gam4 = _mm512_mul_pd(gam2,gam2);
                          rcs  = _mm512_mul_pd(C1140,_mm512_mul_pd(a6,gam4));
                          return (rcs);                  
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7531_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa, // side length of plate
                                              const double * __restrict __ATTR_ALIGN__(64)  pgam0) {
                                     
                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d gam0 = _mm512_load_pd(&pgam0[0]);      
                          const __m512d C1140 = _mm512_set1_pd(114.0f);
                          register __m512d a2,a6,gam2,gam4;
                          register __m512d rcs;
                          a2   = _mm512_mul_pd(a,a);
                          gam2 = _mm512_mul_pd(gam0,gam0);
                          a6   = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          gam4 = _mm512_mul_pd(gam2,gam2);
                          rcs  = _mm512_mul_pd(C1140,_mm512_mul_pd(a6,gam4));
                          return (rcs);                  
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7531_zmm8r8_u(const double * __restrict  pa, // side length of plate
                                              const double * __restrict  pgam0) {
                                     
                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d gam0 = _mm512_loadu_pd(&pgam0[0]);      
                          const __m512d C1140 = _mm512_set1_pd(114.0f);
                          register __m512d a2,a6,gam2,gam4;
                          register __m512d rcs;
                          a2   = _mm512_mul_pd(a,a);
                          gam2 = _mm512_mul_pd(gam0,gam0);
                          a6   = _mm512_mul_pd(a2,_mm512_mul_pd(a2,a2));
                          gam4 = _mm512_mul_pd(gam2,gam2);
                          rcs  = _mm512_mul_pd(C1140,_mm512_mul_pd(a6,gam4));
                          return (rcs);                  
                 }
                 
                 
                 
                 /*
                     
                     Arbitrary shape (no-circular) plates high frequency.
                     Backscatter RCS.
                     Normal incidence.
                     Formula: 7.5-32
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7532_zmm8r8(const __m512d A,
                                            const __m512d gam0) {
                                            
                          const __m512d C12566370614359172953850573533118  = 
                                                          _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,A2,gam02,rat;
                          A2    = _mm512_mul_pd(A,A);
                          gam02 = _mm512_mul_pd(gam0,
                          rat   = _mm512_div_pd(A2,gam02);
                          rcs   = _mm512_mul_pd(C12566370614359172953850573533118,rat);
                          return (rcs);              
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7532_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pA,
                                            const  double * __restrict __ATTR_ALIGN__(64)  gam0) {
                                      
                          register __m512d A    = _mm512_load_pd(&pA[0]);
                          register __m512d gam0 = _mm512_load_pd(&pgam0[0]);     
                          const __m512d C12566370614359172953850573533118  = 
                                                          _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,A2,gam02,rat;
                          A2    = _mm512_mul_pd(A,A);
                          gam02 = _mm512_mul_pd(gam0,
                          rat   = _mm512_div_pd(A2,gam02);
                          rcs   = _mm512_mul_pd(C12566370614359172953850573533118,rat);
                          return (rcs);              
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7532_zmm8r8_u(const double * __restrict  pA,
                                            const  double * __restrict  gam0) {
                                      
                          register __m512d A    = _mm512_loadu_pd(&pA[0]);
                          register __m512d gam0 = _mm512_loadu_pd(&pgam0[0]);     
                          const __m512d C12566370614359172953850573533118  = 
                                                          _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d rcs,A2,gam02,rat;
                          A2    = _mm512_mul_pd(A,A);
                          gam02 = _mm512_mul_pd(gam0,
                          rat   = _mm512_div_pd(A2,gam02);
                          rcs   = _mm512_mul_pd(C12566370614359172953850573533118,rat);
                          return (rcs);              
                }
                
                
                /*
                       Rectangular plate perfectly conducting.
                       Physical Optics backascatter RCS at theta angle.
                       Formula: 7.5-33
                
                */
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7533_zmm8r8(const __m512d A,
                                            const __m512d k0a,
                                            const __m512d k0b,
                                            const __m512d tht,
                                            const __m512d phi,
                                            const __m512d gam0) {
                                            
                          const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d A2,gam02,sint,cost,cosp,sinp,fac;
                          register __m512d rcs,trm1,trm2,arg1,arg2,sarg1,sarg2,x0,x1;
                          A2   = _mm512_mul_pd(A,A);
                          cost = xcos(tht);
                          gam02= _mm512_mul_pd(gam0,gam0);
                          sint = xsin(tht);
                          fac  = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),
                                                                                             gam02);  
                          cosp = xcos(phi);
                          arg1 = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                          sarg1= xsin(arg1);
                          trm1 = _mm512_div_pd(sarg1,arg1);
                          x0   = _mm512_mul_pd(trm1,trm1);
                          sinp = xsin(phi);
                          fac  = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                          arg2 = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                          sarg2= xsin(arg2);
                          trm2 = _mm512_div_pd(sarg2,arg2);
                          x1   = _mm512_mul_pd(trm2,trm2);
                          rcs  = _mm512_mul_pd(fac,_mm512_mul_pd(x0,x1));
                          return (rcs);                 
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7533_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pA,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0b,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi,
                                              const double * __restrict __ATTR_ALIGN__(64) pgam0) {
                             
                          register __m512d A    = _mm512_load_pd(&pA[0]);
                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          register __m512d k0b  = _mm512_load_pd(&pk0b[0]);
                          register __m512d tht  = _mm512_load_pd(&ptht[0]);
                          register __m512d phi  = _mm512_load_pd(&pphi[0]);
                          register __m512d gam0 = _mm512_load_pd(&pgam0[0]);             
                          const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d A2,gam02,sint,cost,cosp,sinp,fac;
                          register __m512d rcs,trm1,trm2,arg1,arg2,sarg1,sarg2,x0,x1;
                          A2   = _mm512_mul_pd(A,A);
                          cost = xcos(tht);
                          gam02= _mm512_mul_pd(gam0,gam0);
                          sint = xsin(tht);
                          fac  = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),
                                                                                             gam02);  
                          cosp = xcos(phi);
                          arg1 = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                          sarg1= xsin(arg1);
                          trm1 = _mm512_div_pd(sarg1,arg1);
                          x0   = _mm512_mul_pd(trm1,trm1);
                          sinp = xsin(phi);
                          fac  = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                          arg2 = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                          sarg2= xsin(arg2);
                          trm2 = _mm512_div_pd(sarg2,arg2);
                          x1   = _mm512_mul_pd(trm2,trm2);
                          rcs  = _mm512_mul_pd(fac,_mm512_mul_pd(x0,x1));
                          return (rcs);                 
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7533_zmm8r8_u(const double * __restrict  pA,
                                              const double * __restrict pk0a,
                                              const double * __restrict pk0b,
                                              const double * __restrict  ptht,
                                              const double * __restrict  pphi,
                                              const double * __restrict  pgam0) {
                             
                          register __m512d A    = _mm512_loadu_pd(&pA[0]);
                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d k0b  = _mm512_loadu_pd(&pk0b[0]);
                          register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                          register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                          register __m512d gam0 = _mm512_loadu_pd(&pgam0[0]);             
                          const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d A2,gam02,sint,cost,cosp,sinp,fac;
                          register __m512d rcs,trm1,trm2,arg1,arg2,sarg1,sarg2,x0,x1;
                          A2   = _mm512_mul_pd(A,A);
                          cost = xcos(tht);
                          gam02= _mm512_mul_pd(gam0,gam0);
                          sint = xsin(tht);
                          fac  = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),
                                                                                             gam02);  
                          cosp = xcos(phi);
                          arg1 = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                          sarg1= xsin(arg1);
                          trm1 = _mm512_div_pd(sarg1,arg1);
                          x0   = _mm512_mul_pd(trm1,trm1);
                          sinp = xsin(phi);
                          fac  = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                          arg2 = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                          sarg2= xsin(arg2);
                          trm2 = _mm512_div_pd(sarg2,arg2);
                          x1   = _mm512_mul_pd(trm2,trm2);
                          rcs  = _mm512_mul_pd(fac,_mm512_mul_pd(x0,x1));
                          return (rcs);                 
                }
                
                
                /*
                      Triangular plates.
                      High-frequency backscatter RCS by methods of 
                      Physical Optics.
                      Formula 7.5-59
                */
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7559_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d gam0,
                                            const __m512d k0a,
                                            const __m512d tht,
                                            const __m512d k0b,
                                            const __m512d phi) {
                                            
                         const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                         const __m512d C05                               = _mm512_set1_pd(0.5f);
                         const __m512d C025                              = _mm512_set1_pd(0.25f);
                         register __m512d A2,alp,beta,sint,beta2,sinp,cosp,sin2a,x0,x1,cost;
                         register __m512d rcs,sinb,sina,sinb2,alp2,fac,gam02,aab,den,trm1,trm2,trm3;
                         x0    = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                         sint  = xsin(tht);
                         A2    = _mm512_mul_pd(x0,x0);
                         cosp  = xcos(phi);
                         gam02 = _mm512_mul_pd(gam0,gam0);
                         alp   = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                         fac   = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                         sinp  = xsin(phi);
                         beta  = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                         cost  = xcos(tht);
                         sinb  = xsin(beta);
                         beta2 = _mm512_mul_pd(beta,C05);
                         fac   = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                         sina  = xsin(alp);
                         aab   = _mm512_div_pd(_mm512_add_pd(a,a),b);
                         sinb2 = xsin(beta2);
                         trm1  = _mm512_fmsub_pd(sina,sina,
                                             _mm512_mul_pd(sinb2,sinb2));
                         x1    = xsin(_mm512_add_pd(alp,alp));
                         den   = _mm512_fmsub_pd(alp,alp,
                                             _mm512_mul_pd(beta2,beta2));
                         trm3  = _mm512_fmsub_pd(cosp,sinb,
                                             _mm512_mul_pd(sinp,x1));
                         trm2  = _mm512_mul_pd(C025,_mm512_mul_pd(sinp,sinp));
                         trm3  = _mm512_mul_pd(aab,trm3);
                         x0    = _mm512_mul_pd(trm1,trm1);
                         x1    = _mm512_mul_pd(trm3,trm3);
                         alp   = _mm512_fmadd_pd(trm2,x1,x0);
                         beta  = _mm512_div_pd(alp,den);
                         rcs   = _mm512_mul_pd(fac,beta);
                         return (rcs);
                  }
                  
                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7559_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0b,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {
                              
                         register __m512d a    = _mm512_load_pd(&pa[0]);
                         register __m512d b    = _mm512_load_pd(&pb[0]);
                         register __m512d gam0 = _mm512_load_pd(&pgam0[0]);
                         register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                         register __m512d tht  = _mm512_load_pd(&ptht[0]);
                         register __m512d k0b  = _mm512_load_pd(&pk0b[0]);
                         register __m512d phi  = _mm512_load_pd(&pphi[0]);             
                         const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                         const __m512d C05                               = _mm512_set1_pd(0.5f);
                         const __m512d C025                              = _mm512_set1_pd(0.25f);
                         register __m512d A2,alp,beta,sint,beta2,sinp,cosp,sin2a,x0,x1,cost;
                         register __m512d rcs,sinb,sina,sinb2,alp2,fac,gam02,aab,den,trm1,trm2,trm3;
                         x0    = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                         sint  = xsin(tht);
                         A2    = _mm512_mul_pd(x0,x0);
                         cosp  = xcos(phi);
                         gam02 = _mm512_mul_pd(gam0,gam0);
                         alp   = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                         fac   = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                         sinp  = xsin(phi);
                         beta  = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                         cost  = xcos(tht);
                         sinb  = xsin(beta);
                         beta2 = _mm512_mul_pd(beta,C05);
                         fac   = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                         sina  = xsin(alp);
                         aab   = _mm512_div_pd(_mm512_add_pd(a,a),b);
                         sinb2 = xsin(beta2);
                         trm1  = _mm512_fmsub_pd(sina,sina,
                                             _mm512_mul_pd(sinb2,sinb2));
                         x1    = xsin(_mm512_add_pd(alp,alp));
                         den   = _mm512_fmsub_pd(alp,alp,
                                             _mm512_mul_pd(beta2,beta2));
                         trm3  = _mm512_fmsub_pd(cosp,sinb,
                                             _mm512_mul_pd(sinp,x1));
                         trm2  = _mm512_mul_pd(C025,_mm512_mul_pd(sinp,sinp));
                         trm3  = _mm512_mul_pd(aab,trm3);
                         x0    = _mm512_mul_pd(trm1,trm1);
                         x1    = _mm512_mul_pd(trm3,trm3);
                         alp   = _mm512_fmadd_pd(trm2,x1,x0);
                         beta  = _mm512_div_pd(alp,den);
                         rcs   = _mm512_mul_pd(fac,beta);
                         return (rcs);
                  }
                  
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7559_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pgam0,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  ptht,
                                              const double * __restrict  pk0b,
                                              const double * __restrict  pphi) {
                              
                         register __m512d a    = _mm512_loadu_pd(&pa[0]);
                         register __m512d b    = _mm512_loadu_pd(&pb[0]);
                         register __m512d gam0 = _mm512_loadu_pd(&pgam0[0]);
                         register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                         register __m512d k0b  = _mm512_loadu_pd(&pk0b[0]);
                         register __m512d phi  = _mm512_loadu_pd(&pphi[0]);             
                         const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                         const __m512d C05                               = _mm512_set1_pd(0.5f);
                         const __m512d C025                              = _mm512_set1_pd(0.25f);
                         register __m512d A2,alp,beta,sint,beta2,sinp,cosp,sin2a,x0,x1,cost;
                         register __m512d rcs,sinb,sina,sinb2,alp2,fac,gam02,aab,den,trm1,trm2,trm3;
                         x0    = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                         sint  = xsin(tht);
                         A2    = _mm512_mul_pd(x0,x0);
                         cosp  = xcos(phi);
                         gam02 = _mm512_mul_pd(gam0,gam0);
                         alp   = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                         fac   = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                         sinp  = xsin(phi);
                         beta  = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                         cost  = xcos(tht);
                         sinb  = xsin(beta);
                         beta2 = _mm512_mul_pd(beta,C05);
                         fac   = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                         sina  = xsin(alp);
                         aab   = _mm512_div_pd(_mm512_add_pd(a,a),b);
                         sinb2 = xsin(beta2);
                         trm1  = _mm512_fmsub_pd(sina,sina,
                                             _mm512_mul_pd(sinb2,sinb2));
                         x1    = xsin(_mm512_add_pd(alp,alp));
                         den   = _mm512_fmsub_pd(alp,alp,
                                             _mm512_mul_pd(beta2,beta2));
                         trm3  = _mm512_fmsub_pd(cosp,sinb,
                                             _mm512_mul_pd(sinp,x1));
                         trm2  = _mm512_mul_pd(C025,_mm512_mul_pd(sinp,sinp));
                         trm3  = _mm512_mul_pd(aab,trm3);
                         x0    = _mm512_mul_pd(trm1,trm1);
                         x1    = _mm512_mul_pd(trm3,trm3);
                         alp   = _mm512_fmadd_pd(trm2,x1,x0);
                         beta  = _mm512_div_pd(alp,den);
                         rcs   = _mm512_mul_pd(fac,beta);
                         return (rcs);
                  }
                  
                  
                   /*
                      Triangular plates.
                      High-frequency backscatter RCS by methods of 
                      Physical Optics (wave incident in plane phi = 0).
                      Formula 7.5-60
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7560_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d gam0,
                                            const __m512d tht,
                                            const __m512d phi,
                                            const __m512d k0a) {
                                            
                         const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                         const __m512d C05                               = _mm512_set1_pd(0.5f);
                         const __m512d C40                               = _mm512_set1_pd(4.0f);
                         register __m512d A2,gam02,cost,alp,fac,sina,sina4,alp4;
                         register __m512d rcs,sin2a,alp2,falp4,x0,x1;
                         x0   = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                         cost = xcos(tht);
                         gam02= _mm512_mul_pd(gam0,gam0);
                         A2   = _mm512_mul_pd(x0,x0);
                         sint = xsin(tht);
                         cosp = xcos(phi);
                         fac  = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                         alp  = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                         x1   = _mm512_mul_pd(alp,alp);
                         sina = xsin(alp);
                         alp4 = _mm512_mul_pd(x1,x1);
                         falp4= _mm512_mul_pd(C40,alp4);
                         alp2 = _mm512_add_pd(alp,alp);
                         sin2a= xsin(alp2);
                         x0   = _mm512_sub_pd(sin2a,alp2);
                         fac  = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                         x1   = _mm512_mul_pd(x0,x0);
                         A2   = _mm512_mul_pd(sina,sina);
                         gam02= _mm512_div_pd(x1,falp4);
                         x0   = _mm512_mul_pd(A2,A2);
                         x1   = _mm512_div_pd(x0,alp4);
                         rcs  = _mm512_mul_pd(fac,_mm512_add_pd(x1,gam02));
                         return (rcs);                        
                }
                
                
                 __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7560_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {
                              
                         register __m512d a    = _mm512_load_pd(&pa[0]);
                         register __m512d b    = _mm512_load_pd(&pb[0]);
                         register __m512d gam0 = _mm512_load_pd(&pgam0[0]);
                         register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                         register __m512d tht  = _mm512_load_pd(&ptht[0]);
                         register __m512d phi  = _mm512_load_pd(&pphi[0]); 
                                            
                         const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                         const __m512d C05                               = _mm512_set1_pd(0.5f);
                         const __m512d C40                               = _mm512_set1_pd(4.0f);
                         register __m512d A2,gam02,cost,alp,fac,sina,sina4,alp4;
                         register __m512d rcs,sin2a,alp2,falp4,x0,x1;
                         x0   = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                         cost = xcos(tht);
                         gam02= _mm512_mul_pd(gam0,gam0);
                         A2   = _mm512_mul_pd(x0,x0);
                         sint = xsin(tht);
                         cosp = xcos(phi);
                         fac  = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                         alp  = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                         x1   = _mm512_mul_pd(alp,alp);
                         sina = xsin(alp);
                         alp4 = _mm512_mul_pd(x1,x1);
                         falp4= _mm512_mul_pd(C40,alp4);
                         alp2 = _mm512_add_pd(alp,alp);
                         sin2a= xsin(alp2);
                         x0   = _mm512_sub_pd(sin2a,alp2);
                         fac  = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                         x1   = _mm512_mul_pd(x0,x0);
                         A2   = _mm512_mul_pd(sina,sina);
                         gam02= _mm512_div_pd(x1,falp4);
                         x0   = _mm512_mul_pd(A2,A2);
                         x1   = _mm512_div_pd(x0,alp4);
                         rcs  = _mm512_mul_pd(fac,_mm512_add_pd(x1,gam02));
                         return (rcs);                        
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7560_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pgam0,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  ptht,
                                              const double * __restrict  pphi) {
                              
                         register __m512d a    = _mm512_loadu_pd(&pa[0]);
                         register __m512d b    = _mm512_loadu_pd(&pb[0]);
                         register __m512d gam0 = _mm512_loadu_pd(&pgam0[0]);
                         register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                         register __m512d phi  = _mm512_loadu_pd(&pphi[0]); 
                                            
                         const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                         const __m512d C05                               = _mm512_set1_pd(0.5f);
                         const __m512d C40                               = _mm512_set1_pd(4.0f);
                         register __m512d A2,gam02,cost,alp,fac,sina,sina4,alp4;
                         register __m512d rcs,sin2a,alp2,falp4,x0,x1;
                         x0   = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                         cost = xcos(tht);
                         gam02= _mm512_mul_pd(gam0,gam0);
                         A2   = _mm512_mul_pd(x0,x0);
                         sint = xsin(tht);
                         cosp = xcos(phi);
                         fac  = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                         alp  = _mm512_mul_pd(k0a,_mm512_mul_pd(sint,cosp));
                         x1   = _mm512_mul_pd(alp,alp);
                         sina = xsin(alp);
                         alp4 = _mm512_mul_pd(x1,x1);
                         falp4= _mm512_mul_pd(C40,alp4);
                         alp2 = _mm512_add_pd(alp,alp);
                         sin2a= xsin(alp2);
                         x0   = _mm512_sub_pd(sin2a,alp2);
                         fac  = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                         x1   = _mm512_mul_pd(x0,x0);
                         A2   = _mm512_mul_pd(sina,sina);
                         gam02= _mm512_div_pd(x1,falp4);
                         x0   = _mm512_mul_pd(A2,A2);
                         x1   = _mm512_div_pd(x0,alp4);
                         rcs  = _mm512_mul_pd(fac,_mm512_add_pd(x1,gam02));
                         return (rcs);                        
                }
                
                
                /*
                      Triangular plates.
                      High-frequency backscatter RCS by methods of 
                      Physical Optics (wave incident in plane phi = PI/2).
                      Formula 7.5-61 
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7561_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d k0b,
                                            const __m512d tht,
                                            const __m512d phi,
                                            const __m512d gam0) {
                                            
                        const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                        
                        register __m512d A2,gam02,fac,cost,x0,x1,sint,sinp,beta2,sbeta,sbeta4,beta4;
                        register __m512d rcs;
                        x0    = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                        cost  = xcos(tht);
                        gam02 = _mm512_mul_pd(gam0,gam0);
                        sint  = xsin(tht);
                        A2    = _mm512_mul_pd(x0,x0);
                        sinp  = xsin(phi);
                        fac   = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                        x0    = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                        fac   = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                        beta2 = _mm512_mul_pd(x0,C05);
                        sbeta = xsin(beta2);
                        x0    = _mm512_mul_pd(beta2,beta2);
                        x1    = _mm512_mul_pd(sbeta,sbeta);
                        beta4 = _mm512_mul_pd(x0,x0);   
                        sbeta4= _mm512_mul_pd(x1,x1);
                        sint  = _mm512_div_pd(sbeta4,beta4);
                        rcs   = _mm512_mul_pd(fac,sint);
                        return (rcs);       
               }
               
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7561_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0b,
                                              const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {
                              
                        register __m512d a    = _mm512_load_pd(&pa[0]);
                        register __m512d b    = _mm512_load_pd(&pb[0]);
                        register __m512d gam0 = _mm512_load_pd(&pgam0[0]);
                        register __m512d k0a  = _mm512_load_pd(&pk0b[0]);
                        register __m512d tht  = _mm512_load_pd(&ptht[0]);
                        register __m512d phi  = _mm512_load_pd(&pphi[0]);
                        const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                        
                        register __m512d A2,gam02,fac,cost,x0,x1,sint,sinp,beta2,sbeta,sbeta4,beta4;
                        register __m512d rcs;
                        x0    = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                        cost  = xcos(tht);
                        gam02 = _mm512_mul_pd(gam0,gam0);
                        sint  = xsin(tht);
                        A2    = _mm512_mul_pd(x0,x0);
                        sinp  = xsin(phi);
                        fac   = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                        x0    = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                        fac   = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                        beta2 = _mm512_mul_pd(x0,C05);
                        sbeta = xsin(beta2);
                        x0    = _mm512_mul_pd(beta2,beta2);
                        x1    = _mm512_mul_pd(sbeta,sbeta);
                        beta4 = _mm512_mul_pd(x0,x0);   
                        sbeta4= _mm512_mul_pd(x1,x1);
                        sint  = _mm512_div_pd(sbeta4,beta4);
                        rcs   = _mm512_mul_pd(fac,sint);
                        return (rcs);       
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512d rcs_f7561_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pgam0,
                                              const double * __restrict  pk0b,
                                              const double * __restrict  ptht,
                                              const double * __restrict  pphi) {
                              
                        register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d gam0 = _mm512_loadu_pd(&pgam0[0]);
                        register __m512d k0a  = _mm512_loadu_pd(&pk0b[0]);
                        register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                        register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                        const __m512d C12566370614359172953850573533118 = 
                                                       _mm512_set1_pd(12.566370614359172953850573533118);
                        
                        register __m512d A2,gam02,fac,cost,x0,x1,sint,sinp,beta2,sbeta,sbeta4,beta4;
                        register __m512d rcs;
                        x0    = _mm512_mul_pd(_mm512_mul_pd(a,b),C05);
                        cost  = xcos(tht);
                        gam02 = _mm512_mul_pd(gam0,gam0);
                        sint  = xsin(tht);
                        A2    = _mm512_mul_pd(x0,x0);
                        sinp  = xsin(phi);
                        fac   = _mm512_div_pd(_mm512_mul_pd(C12566370614359172953850573533118,A2),gam02);
                        x0    = _mm512_mul_pd(k0b,_mm512_mul_pd(sint,sinp));
                        fac   = _mm512_mul_pd(fac,_mm512_mul_pd(cost,cost));
                        beta2 = _mm512_mul_pd(x0,C05);
                        sbeta = xsin(beta2);
                        x0    = _mm512_mul_pd(beta2,beta2);
                        x1    = _mm512_mul_pd(sbeta,sbeta);
                        beta4 = _mm512_mul_pd(x0,x0);   
                        sbeta4= _mm512_mul_pd(x1,x1);
                        sint  = _mm512_div_pd(sbeta4,beta4);
                        rcs   = _mm512_mul_pd(fac,sint);
                        return (rcs);       
               }
               
                                            
                  
                

      } // radiolocation


} // gms


























#endif /*__GMS_RCS_PLANAR_SURF_ZMM8R8_HPP__*/
