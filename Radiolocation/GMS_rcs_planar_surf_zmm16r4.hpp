

#ifndef __GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__
#define __GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__ 180420231543



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

    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_PLANAR_SURF_ZMM16R4_FULLVER =
      1000U*GMS_RCS_PLANAR_SURF_ZMM16R4_MAJOR+
      100U*GMS_RCS_PLANAR_SURF_ZMM16R4_MINOR+
      10U*GMS_RCS_PLANAR_SURF_ZMM16R4_MICRO;
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_CREATION_DATE = "18-04-2023 15:43 PM +00200 (TUE 18 APR 2023 GMT+2)";
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_PLANAR_SURF_ZMM16R4_DESCRIPTION   = "AVX512 optimized Planar Surfaces Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
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
                   void zi_f716_zmm16r4(const __m512 tht,
                                        const __m512 mur,
                                        const __m512 mui,
                                        const __m512 epsr,
                                        const __m512 epsi,
                                        __m512 * __restrict zr,
                                        __m512 * __restrict zi) {

                        register __m512 cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcosf(tht);
                        wrkc = _mm512_setzero_ps();
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_ps(cost);
                        csqrt_zmm16r4(divr,divi,&wrkc,&csqr,&csqi);
                        *zr = _mm512_mul_ps(invc,csqr);
                        *zi = _mm512_mul_ps(invc,csqi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                        const float * __restrict __ATTR_ALIGN__(64) pmur,
                                        const float * __restrict __ATTR_ALIGN__(64) pmui,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsr,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsi,
                                        float * __restrict __ATTR_ALIGN__(64) zr,
                                        float * __restrict __ATTR_ALIGN__(64)  zi) {

                        register __m512 tht  = _mm512_load_ps(&ptht[0]);
                        register __m512 mur  = _mm512_load_ps(&pmur[0]);
                        register __m512 mui  = _mm512_load_ps(&pmui[0]);
                        register __m512 epsr = _mm512_load_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_load_ps(&pepsi[0]);
                        register __m512 cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcosf(tht);
                        wrkc = _mm512_setzero_ps();
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_ps(cost);
                        csqrt_zmm16r4(divr,divi,&wrkc,&csqr,&csqi);
                        _mm512_store_ps(&zr[0] ,_mm512_mul_ps(invc,csqr));
                        _mm512_store_ps(&zi[0] ,_mm512_mul_ps(invc,csqi));
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void zi_f716_zmm16r4_u(const float * __restrict  ptht,
                                        const float * __restrict  pmur,
                                        const float * __restrict  pmui,
                                        const float * __restrict  pepsr,
                                        const float * __restrict  pepsi,
                                        float * __restrict  zr,
                                        float * __restrict   zi) {

                        register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                        register __m512 mur  = _mm512_loadu_ps(&pmur[0]);
                        register __m512 mui  = _mm512_loadu_ps(&pmui[0]);
                        register __m512 epsr = _mm512_loadu_ps(&pepsr[0]);
                        register __m512 epsi = _mm512_loadu_ps(&pepsi[0]);
                        register __m512 cost,invc,divr,divi,csqr,csqi,wrkc;
                        cost = xcosf(tht);
                        wrkc = _mm512_setzero_ps();
                        cdiv_zmm16r4(mur,mui,epsr,epsi,&divr,&divi);
                        invc = _mm512_rcp14_ps(cost);
                        csqrt_zmm16r4(divr,divi,&wrkc,&csqr,&csqi);
                        _mm512_storeu_ps(&zr[0] ,_mm512_mul_ps(invc,csqr));
                        _mm512_storeu_ps(&zi[0] ,_mm512_mul_ps(invc,csqi));
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
                   void R_f714_zmm16r4( const __m512 tht1,
                                        const __m512 mur1,
                                        const __m512 mui1,
                                        const __m512 epsr1,
                                        const __m512 epsi1,
                                        const __m512 tht2,
                                        const __m512 mur2,
                                        const __m512 mui2,
                                        const __m512 epsr2,
                                        const __m512 epsi2,
                                        __m512 * __restrict Rr,
                                        __m512 * __restrict Ri) {

                     using namespace gms::math;
                     register __m512 z1r,z1i,z2r,z2i;
                     register __m512 t0r,t0i,t1r,t1i;
                     zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     zi_f716_zmm16r4(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     t0r = _mm512_sub_ps(z1r,z2r);
                     t1r = _mm512_add_ps(z1r,z2r);
                     t0i = _mm512_sub_ps(z1i,z2i);
                     t1i = _mm512_add_ps(z1i,z2i);
                     t0r = negate_zmm16r4(t0r);
                     t0i = negate_zmm16r4(t0i);
                     cdiv_zmm16r4(t0r,t0i,t1r,t1i,*Rr,*Ri);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f714_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                        const float * __restrict __ATTR_ALIGN__(64) pmur1,
                                        const float * __restrict __ATTR_ALIGN__(64) pmui1,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsr1,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsi1,
                                        const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                        const float * __restrict __ATTR_ALIGN__(64) pmur2,
                                        const float * __restrict __ATTR_ALIGN__(64) pmui2,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsr2,
                                        const float * __restrict __ATTR_ALIGN__(64) pepsi2,
                                        float * __restrict __ATTR_ALIGN__(64) Rr,
                                        float * __restrict __ATTR_ALIGN__(64) Ri) {

                     using namespace gms::math;
                     register __m512 tht1  = _mm512_load_ps(&ptht1[0]);
                     register __m512 mur1  = _mm512_load_ps(&pmur1[0]);
                     register __m512 mui1  = _mm512_load_ps(&pmui1[0]);
                     register __m512 epsr1 = _mm512_load_ps(&pepsr1[0]);
                     register __m512 epsi1 = _mm512_load_ps(&pepsi1[0]);
                     register __m512 tht2  = _mm512_load_ps(&ptht2[0]);
                     register __m512 mur2  = _mm512_load_ps(&pmur2[0]);
                     register __m512 mui2  = _mm512_load_ps(&pmui2[0]);
                     register __m512 epsr2 = _mm512_load_ps(&pepsr2[0]);
                     register __m512 epsi2 = _mm512_load_ps(&pepsi2[0]);
                     register __m512 z1r,z1i,z2r,z2i;
                     register __m512 t0r,t0i,t1r,t1i,resr,resi;
                     zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     zi_f716_zmm16r4(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     t0r = _mm512_sub_ps(z1r,z2r);
                     t1r = _mm512_add_ps(z1r,z2r);
                     t0i = _mm512_sub_ps(z1i,z2i);
                     t1i = _mm512_add_ps(z1i,z2i);
                     t0r = negate_zmm16r4(t0r);
                     t0i = negate_zmm16r4(t0i);
                     cdiv_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                     _mm512_store_ps(&Rr[0], resr);
                     _mm512_store_ps(&Ri[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f714_zmm16r4_u( const float * __restrict  ptht1,
                                        const float * __restrict  pmur1,
                                        const float * __restrict  pmui1,
                                        const float * __restrict  pepsr1,
                                        const float * __restrict  pepsi1,
                                        const float * __restrict  ptht2,
                                        const float * __restrict  pmur2,
                                        const float * __restrict  pmui2,
                                        const float * __restrict  pepsr2,
                                        const float * __restrict  pepsi2,
                                        float * __restrict  Rr,
                                        float * __restrict  Ri) {

                     using namespace gms::math;
                     register __m512 tht1  = _mm512_loadu_ps(&ptht1[0]);
                     register __m512 mur1  = _mm512_loadu_ps(&pmur1[0]);
                     register __m512 mui1  = _mm512_loadu_ps(&pmui1[0]);
                     register __m512 epsr1 = _mm512_loadu_ps(&pepsr1[0]);
                     register __m512 epsi1 = _mm512_loadu_ps(&pepsi1[0]);
                     register __m512 tht2  = _mm512_loadu_ps(&ptht2[0]);
                     register __m512 mur2  = _mm512_loadu_ps(&pmur2[0]);
                     register __m512 mui2  = _mm512_loadu_ps(&pmui2[0]);
                     register __m512 epsr2 = _mm512_loadu_ps(&pepsr2[0]);
                     register __m512 epsi2 = _mm512_loadu_ps(&pepsi2[0]);
                     register __m512 z1r,z1i,z2r,z2i;
                     register __m512 t0r,t0i,t1r,t1i,resr,resi;
                     zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     zi_f716_zmm16r4(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     t0r = _mm512_sub_ps(z1r,z2r);
                     t1r = _mm512_add_ps(z1r,z2r);
                     t0i = _mm512_sub_ps(z1i,z2i);
                     t1i = _mm512_add_ps(z1i,z2i);
                     t0r = negate_zmm16r4(t0r);
                     t0i = negate_zmm16r4(t0i);
                     cdiv_zmm16r4(t0r,t0i,t1r,t1i,&resr,&resi);
                     _mm512_storeu_ps(&Rr[0], resr);
                     _mm512_storeu_ps(&Ri[0], resi);
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
                   void T_f715_zmm16r4( const __m512 tht1,
                                        const __m512 mur1,
                                        const __m512 mui1,
                                        const __m512 epsr1,
                                        const __m512 epsi1,
                                        const __m512 tht2,
                                        const __m512 mur2,
                                        const __m512 mui2,
                                        const __m512 epsr2,
                                        const __m512 epsi2,
                                        __m512 * __restrict Tr,
                                        __m512 * __restrict Ti) {

                     const __m512 _2 = _mm512_set1_ps(2.0f);
                     register __m512 z1r,z1i,z2r,z2i;
                     register __m512 t0r,t0i,t1r,t1i;
                     zi_f716_zmm16r4(tht2,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                     zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                     t0r = _mm512_mul_ps(_2,z2r);
                     t1r = _mm512_add_ps(z1r,z2r);
                     t0i = _mm512_mul_ps(_2,z2i);
                     t1i = _mm512_add_ps(z1i,z2i);
                     cdiv_zmm16r4(t0r,t0i,t1r,t1i,*Tr,*Ti);
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
                   __m512 R_f7117_zmm16r4(const __m512 tht,
                                          const __m512 eps1,
                                          const __m512 eps2) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 e1e2,sqr1,sqr2,num,den,R;
                          register __m512 cost,sint,x0,x1;
                          e1e2 = _mm512_div_ps(eps1,eps2);
                          cost = xcosf(tht);
                          sqr1 = _mm512_sqrt_ps(e1e2);
                          sint = xsinf(tht);
                          x0   = _mm512_mul_ps(sqr1,cost);
                          x1   = _mm512_mul_ps(_mm512_sub_ps(_1,e1e2),
                                                     _mm512_mul_ps(sint,sint));
                          sqr2 = _mm512_sqrt_ps(x1);
                          num  = _mm512_sub_ps(x0,x1);
                          den  = _mm512_add_ps(x0,x1);
                          R    = _mm512_div_ps(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 R_f7117_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps2) {

                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps2= _mm512_load_ps(&peps2[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 e1e2,sqr1,sqr2,num,den,R;
                          register __m512 cost,sint,x0,x1;
                          e1e2 = _mm512_div_ps(eps1,eps2);
                          cost = xcosf(tht);
                          sqr1 = _mm512_sqrt_ps(e1e2);
                          sint = xsinf(tht);
                          x0   = _mm512_mul_ps(sqr1,cost);
                          x1   = _mm512_mul_ps(_mm512_sub_ps(_1,e1e2),
                                                     _mm512_mul_ps(sint,sint));
                          sqr2 = _mm512_sqrt_ps(x1);
                          num  = _mm512_sub_ps(x0,x1);
                          den  = _mm512_add_ps(x0,x1);
                          R    = _mm512_div_ps(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 R_f7117_zmm16r4_u(const float * __restrict  ptht,
                                            const float * __restrict  peps1,
                                            const float * __restrict  peps2) {

                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps2= _mm512_loadu_ps(&peps2[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 e1e2,sqr1,sqr2,num,den,R;
                          register __m512 cost,sint,x0,x1;
                          e1e2 = _mm512_div_ps(eps1,eps2);
                          cost = xcosf(tht);
                          sqr1 = _mm512_sqrt_ps(e1e2);
                          sint = xsinf(tht);
                          x0   = _mm512_mul_ps(sqr1,cost);
                          x1   = _mm512_mul_ps(_mm512_sub_ps(_1,e1e2),
                                                     _mm512_mul_ps(sint,sint));
                          sqr2 = _mm512_sqrt_ps(x1);
                          num  = _mm512_sub_ps(x0,x1);
                          den  = _mm512_add_ps(x0,x1);
                          R    = _mm512_div_ps(num,den);
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
                   __m512 R_f7118_zmm16r4(const __m512 tht,
                                          const __m512 eps1,
                                          const __m512 eps2) {

                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 e1e2,sqr1,sqr2,num,den,R;
                          register __m512 cost,sint,x0,x1;
                          e1e2 = _mm512_div_ps(eps2,eps1);
                          cost = xcosf(tht);
                          sqr1 = _mm512_sqrt_ps(e1e2);
                          sint = xsinf(tht);
                          x0   = _mm512_mul_ps(sqr1,cost);
                          x1   = _mm512_mul_ps(_mm512_sub_ps(_1,e1e2),
                                                     _mm512_mul_ps(sint,sint));
                          sqr2 = _mm512_sqrt_ps(x1);
                          num  = _mm512_sub_ps(x0,x1);
                          den  = _mm512_add_ps(x0,x1);
                          R    = _mm512_div_ps(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 R_f7118_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps2) {

                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps2= _mm512_load_ps(&peps2[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 e1e2,sqr1,sqr2,num,den,R;
                          register __m512 cost,sint,x0,x1;
                          e1e2 = _mm512_div_ps(eps2,eps1);
                          cost = xcosf(tht);
                          sqr1 = _mm512_sqrt_ps(e1e2);
                          sint = xsinf(tht);
                          x0   = _mm512_mul_ps(sqr1,cost);
                          x1   = _mm512_mul_ps(_mm512_sub_ps(_1,e1e2),
                                                     _mm512_mul_ps(sint,sint));
                          sqr2 = _mm512_sqrt_ps(x1);
                          num  = _mm512_sub_ps(x0,x1);
                          den  = _mm512_add_ps(x0,x1);
                          R    = _mm512_div_ps(num,den);
                          return (R);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 R_f7118_zmm16r4_u(const float * __restrict  ptht,
                                            const float * __restrict  peps1,
                                            const float * __restrict  peps2) {

                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                          register __m512 eps2= _mm512_loadu_ps(&peps2[0]);
                          const __m512 _1 = _mm512_set1_ps(1.0f);
                          register __m512 e1e2,sqr1,sqr2,num,den,R;
                          register __m512 cost,sint,x0,x1;
                          e1e2 = _mm512_div_ps(eps2,eps1);
                          cost = xcosf(tht);
                          sqr1 = _mm512_sqrt_ps(e1e2);
                          sint = xsinf(tht);
                          x0   = _mm512_mul_ps(sqr1,cost);
                          x1   = _mm512_mul_ps(_mm512_sub_ps(_1,e1e2),
                                                     _mm512_mul_ps(sint,sint));
                          sqr2 = _mm512_sqrt_ps(x1);
                          num  = _mm512_sub_ps(x0,x1);
                          den  = _mm512_add_ps(x0,x1);
                          R    = _mm512_div_ps(num,den);
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
                   void R_f7123_zmm16r4(const __m512 tht,
                                        const __m512 eps2,
                                        const __m512 eps1,
                                        __m512 * __restrict Rr,
                                        __m512 * __restrict Ri) {

                        const __m512 n2 = _mm512_set1_ps(-2.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 sint,cost,e2e1,rat,arg,atarg;
                        register __m512 x0,x1;
                        e2e1 = _mm512_div_ps(eps2,eps1);
                        cost = xcosf(tht);
                        ear  = _mm512_setzero_ps();
                        sint = xsinf(tht);
                        x0   = _mm512_mul_ps(cost,cost);
                        x1   = _mm512_sub_ps(_mm512_mul_ps(sint,sint),e2e1);
                        rat  = _mm512_div_ps(x1,x0);
                        arg  = _mm512_sqrt_ps(rat);
                        atarg= xatanf(arg);
                        eai  = _mm512_mul_ps(n2,atarg);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        *Rr = cer;
                        *Ri = cei;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7123_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) ptht,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps2
                                            float * __restrict __ATTR_ALIGN__(64) Rr,
                                            float * __restrict __ATTR_ALIGN__(64) Ri) {

                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 eps1= _mm512_load_ps(&peps1[0]);
                        register __m512 eps2= _mm512_load_ps(&peps2[0]);
                        const __m512 n2 = _mm512_set1_ps(-2.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 sint,cost,e2e1,rat,arg,atarg;
                        register __m512 x0,x1;
                        e2e1 = _mm512_div_ps(eps2,eps1);
                        cost = xcosf(tht);
                        ear  = _mm512_setzero_ps();
                        sint = xsinf(tht);
                        x0   = _mm512_mul_ps(cost,cost);
                        x1   = _mm512_sub_ps(_mm512_mul_ps(sint,sint),e2e1);
                        rat  = _mm512_div_ps(x1,x0);
                        arg  = _mm512_sqrt_ps(rat);
                        atarg= xatanf(arg);
                        eai  = _mm512_mul_ps(n2,atarg);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        _mm512_store_ps(&Rr[0] ,cer);
                        _mm512_store_ps(&Ri[0] ,cei);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7123_zmm16r4_u(  const float * __restrict  ptht,
                                            const float * __restrict  peps1,
                                            const float * __restrict  peps2
                                            float * __restrict  Rr,
                                            float * __restrict  Ri) {

                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                        register __m512 eps2= _mm512_loadu_ps(&peps2[0]);
                        const __m512 n2 = _mm512_set1_ps(-2.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 sint,cost,e2e1,rat,arg,atarg;
                        register __m512 x0,x1;
                        e2e1 = _mm512_div_ps(eps2,eps1);
                        cost = xcosf(tht);
                        ear  = _mm512_setzero_ps();
                        sint = xsinf(tht);
                        x0   = _mm512_mul_ps(cost,cost);
                        x1   = _mm512_sub_ps(_mm512_mul_ps(sint,sint),e2e1);
                        rat  = _mm512_div_ps(x1,x0);
                        arg  = _mm512_sqrt_ps(rat);
                        atarg= xatanf(arg);
                        eai  = _mm512_mul_ps(n2,atarg);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        _mm512_storeu_ps(&Rr[0] ,cer);
                        _mm512_storeu_ps(&Ri[0] ,cei);
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
                   void R_f7124_zmm16r4(const __m512 tht,
                                        const __m512 eps2,
                                        const __m512 eps1,
                                        __m512 * __restrict Rr,
                                        __m512 * __restrict Ri) {

                        const __m512 n2 = _mm512_set1_ps(-2.0f);
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 sint,cost,e2e1,e1e2,rat,arg,atarg;
                        register __m512 x0,x1;
                        e2e1 = _mm512_div_ps(eps2,eps1);
                        sint = xsinf(tht);
                        ear  = _mm512_setzero_ps();
                        cost = xcosf(tht);
                        x0   = _mm512_mul_ps(e2e1,_mm512_fmsub_ps(sint,sint,_1));
                        e1e2 = _mm512_div_ps(eps1,eps2);
                        x1   = _mm512_mul_ps(e1e2,_mm512_mul_ps(cost,cost));
                        rat  = _mm512_div_ps(x0,x1);
                        arg  = _mm512_sqrt_ps(rat);
                        atarg= xatanf(arg);
                        eai  = _mm512_mul_ps(n2,atarg);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        *Rr  = cer;
                        *Ri  = cei;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7124_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) ptht,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps2
                                            float * __restrict __ATTR_ALIGN__(64) Rr,
                                            float * __restrict __ATTR_ALIGN__(64) Ri) {

                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 eps1= _mm512_load_ps(&peps1[0]);
                        register __m512 eps2= _mm512_load_ps(&peps2[0]);
                        const __m512 n2 = _mm512_set1_ps(-2.0f);
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 sint,cost,e2e1,e1e2,rat,arg,atarg;
                        register __m512 x0,x1;
                        e2e1 = _mm512_div_ps(eps2,eps1);
                        sint = xsinf(tht);
                        ear  = _mm512_setzero_ps();
                        cost = xcosf(tht);
                        x0   = _mm512_mul_ps(e2e1,_mm512_fmsub_ps(sint,sint,_1));
                        e1e2 = _mm512_div_ps(eps1,eps2);
                        x1   = _mm512_mul_ps(e1e2,_mm512_mul_ps(cost,cost));
                        rat  = _mm512_div_ps(x0,x1);
                        arg  = _mm512_sqrt_ps(rat);
                        atarg= xatanf(arg);
                        eai  = _mm512_mul_ps(n2,atarg);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        _mm512_store_ps(&Rr[0] ,cer);
                        _mm512_store_ps(&Ri[0] ,cei);
                }



                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7124_zmm16r4_u(  const float * __restrict  ptht,
                                            const float * __restrict  peps1,
                                            const float * __restrict  peps2
                                            float * __restrict  Rr,
                                            float * __restrict  Ri) {

                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                        register __m512 eps2= _mm512_loadu_ps(&peps2[0]);
                        const __m512 n2 = _mm512_set1_ps(-2.0f);
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 sint,cost,e2e1,e1e2,rat,arg,atarg;
                        register __m512 x0,x1;
                        e2e1 = _mm512_div_ps(eps2,eps1);
                        sint = xsinf(tht);
                        ear  = _mm512_setzero_ps();
                        cost = xcosf(tht);
                        x0   = _mm512_mul_ps(e2e1,_mm512_fmsub_ps(sint,sint,_1));
                        e1e2 = _mm512_div_ps(eps1,eps2);
                        x1   = _mm512_mul_ps(e1e2,_mm512_mul_ps(cost,cost));
                        rat  = _mm512_div_ps(x0,x1);
                        arg  = _mm512_sqrt_ps(rat);
                        atarg= xatanf(arg);
                        eai  = _mm512_mul_ps(n2,atarg);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        _mm512_storeu_ps(&Rr[0] ,cer);
                        _mm512_storeu_ps(&Ri[0] ,cei);
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
                   __m512 D_f7127_zmm16r4(    const __m512 gam0,
                                              const __m512 tht,
                                              const __m512 eps2,
                                              const __m512 eps1) {

                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
                         register __m512 g0pi,ttht,sint,e2e1,sqr,rat;
                         register __m512 D;
                         g0pi = _mm512_mul_ps(gam0,invpi);
                         e2e1 = _mm512_div_ps(eps2,eps1);
                         ttht = xtanf(tht);
                         sint = _mm512_fmsub_ps(sint,sint,e2e1);
                         sqr  = _mm512_sqrt_ps(sint);
                         rat  = _mm512_div_ps(ttht,sqr);
                         D    = _mm512_mul_ps(g0pi,rat);
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 D_f7127_zmm16r4_a(    const float * __restrict __ATTR_ALIGN__(64)  pgam0,
                                                const float * __restrict __ATTR_ALIGN__(64)  ptht,
                                                const float * __restrict __ATTR_ALIGN__(64)  peps2,
                                                const float * __restrict __ATTR_ALIGN__(64)  peps1) {

                         register __m512 gam0= _mm512_load_ps(&pgam0[0]);
                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         register __m512 eps1= _mm512_load_ps(&peps1[0]);
                         register __m512 eps2= _mm512_load_ps(&peps2[0]);
                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
                         register __m512 g0pi,ttht,sint,e2e1,sqr,rat;
                         register __m512 D;
                         g0pi = _mm512_mul_ps(gam0,invpi);
                         e2e1 = _mm512_div_ps(eps2,eps1);
                         ttht = xtanf(tht);
                         sint = _mm512_fmsub_ps(sint,sint,e2e1);
                         sqr  = _mm512_sqrt_ps(sint);
                         rat  = _mm512_div_ps(ttht,sqr);
                         D    = _mm512_mul_ps(g0pi,rat);
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 D_f7127_zmm16r4_u(    const float * __restrict   pgam0,
                                                const float * __restrict   ptht,
                                                const float * __restrict   peps2,
                                                const float * __restrict   peps1) {

                         register __m512 gam0= _mm512_loadu_ps(&pgam0[0]);
                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                         register __m512 eps2= _mm512_loadu_ps(&peps2[0]);
                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
                         register __m512 g0pi,ttht,sint,e2e1,sqr,rat;
                         register __m512 D;
                         g0pi = _mm512_mul_ps(gam0,invpi);
                         e2e1 = _mm512_div_ps(eps2,eps1);
                         ttht = xtanf(tht);
                         sint = _mm512_fmsub_ps(sint,sint,e2e1);
                         sqr  = _mm512_sqrt_ps(sint);
                         rat  = _mm512_div_ps(ttht,sqr);
                         D    = _mm512_mul_ps(g0pi,rat);
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
                   __m512 D_f7128_zmm16r4(    const __m512 gam0,
                                              const __m512 tht,
                                              const __m512 eps2,
                                              const __m512 eps1) {

                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
                         register __m512 g0pi,ttht,sint,e2e1,e1e2,sqr,rat;
                         register __m512 D;
                         g0pi = _mm512_mul_ps(gam0,invpi);
                         e2e1 = _mm512_div_ps(eps2,eps1);
                         ttht = xtanf(tht);
                         sint = _mm512_fmsub_ps(sint,sint,e2e1);
                         e1e2 = _mm512_div_ps(eps1,eps2);
                         sqr  = _mm512_sqrt_ps(sint);
                         rat  = _mm512_div_ps(ttht,sqr);
                         D    = _mm512_mul_ps(g0pi,_mm512_mul_ps(e1e2,rat));
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 D_f7128_zmm16r4_a(    const float * __restrict __ATTR_ALIGN__(64)  pgam0,
                                                const float * __restrict __ATTR_ALIGN__(64)  ptht,
                                                const float * __restrict __ATTR_ALIGN__(64)  peps2,
                                                const float * __restrict __ATTR_ALIGN__(64)  peps1) {

                         register __m512 gam0= _mm512_load_ps(&pgam0[0]);
                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         register __m512 eps1= _mm512_load_ps(&peps1[0]);
                         register __m512 eps2= _mm512_load_ps(&peps2[0]);
                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
                         register __m512 g0pi,ttht,sint,e2e1,e1e2,sqr,rat;
                         register __m512 D;
                         g0pi = _mm512_mul_ps(gam0,invpi);
                         e2e1 = _mm512_div_ps(eps2,eps1);
                         ttht = xtanf(tht);
                         sint = _mm512_fmsub_ps(sint,sint,e2e1);
                         e1e2 = _mm512_div_ps(eps1,eps2);
                         sqr  = _mm512_sqrt_ps(sint);
                         rat  = _mm512_div_ps(ttht,sqr);
                         D    = _mm512_mul_ps(g0pi,_mm512_mul_ps(e1e2,rat));
                         return (D);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 D_f7128_zmm16r4_u(    const float * __restrict   pgam0,
                                                const float * __restrict   ptht,
                                                const float * __restrict   peps2,
                                                const float * __restrict   peps1) {

                         register __m512 gam0= _mm512_loadu_ps(&pgam0[0]);
                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                         register __m512 eps2= _mm512_loadu_ps(&peps2[0]);
                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
                         register __m512 g0pi,ttht,sint,e2e1,e1e2,sqr,rat;
                         register __m512 D;
                         g0pi = _mm512_mul_ps(gam0,invpi);
                         e2e1 = _mm512_div_ps(eps2,eps1);
                         ttht = xtanf(tht);
                         sint = _mm512_fmsub_ps(sint,sint,e2e1);
                         e1e2 = _mm512_div_ps(eps1,eps2);
                         sqr  = _mm512_sqrt_ps(sint);
                         rat  = _mm512_div_ps(ttht,sqr);
                         D    = _mm512_mul_ps(g0pi,_mm512_mul_ps(e1e2,rat));
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
                   void R_f7129_zmm16r4(const __m512 tht,
                                        const __m512 mur1,
                                        const __m512 mui1,
                                        const __m512 epsr1,
                                        const __m512 epsi1,
                                        const __m512 mur2,
                                        const __m512 mui2,
                                        const __m512 epsr2,
                                        const __m512 epsi2,
                                        __m512 * __restrict Rr,
                                        __m512 * __restrict Ri) {

                       register __m512 z1r,z1i,z2r,z2i;
                       register __m512 cost,numr,numi,denr,deni;
                       register __m512 resr,resi;
                       zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcosf(tht);
                       zi_f716_zmm16r4(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_fmsub_ps(z2r,cost,z1r);
                       denr = _mm512_fmadd_ps(z2r,cost,z2r);
                       numi = _mm512_fmsub_ps(z2i,cost,z1i);
                       deni = _mm512_fmadd_ps(z2i,cost,z2i);
                       cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7129_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur1,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui1,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr1,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi1,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur2,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui2,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr2,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi2,
                                          float * __restrict __ATTR_ALIGN__(64) Rr,
                                          float * __restrict __ATTR_ALIGN__(64) Ri) {

                       register __m512 tht  = _mm512_load_ps(&ptht[0]);
                       register __m512 mur1  = _mm512_load_ps(&pmur1[0]);
                       register __m512 mui1  = _mm512_load_ps(&pmui1[0]);
                       register __m512 epsr1 = _mm512_load_ps(&pepsr1[0]);
                       register __m512 epsi1 = _mm512_load_ps(&pepsi1[0]);
                       register __m512 mur2  = _mm512_load_ps(&pmur2[0]);
                       register __m512 mui2  = _mm512_load_ps(&pmui2[0]);
                       register __m512 epsr2 = _mm512_load_ps(&pepsr2[0]);
                       register __m512 epsi2 = _mm512_load_ps(&pepsi2[0]);
                       register __m512 z1r,z1i,z2r,z2i;
                       register __m512 cost,numr,numi,denr,deni;
                       register __m512 resr,resi;
                       zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcosf(tht);
                       zi_f716_zmm16r4(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_fmsub_ps(z2r,cost,z1r);
                       denr = _mm512_fmadd_ps(z2r,cost,z2r);
                       numi = _mm512_fmsub_ps(z2i,cost,z1i);
                       deni = _mm512_fmadd_ps(z2i,cost,z2i);
                       cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7129_zmm16r4_u(const float * __restrict  ptht,
                                          const float * __restrict  pmur1,
                                          const float * __restrict  pmui1,
                                          const float * __restrict  pepsr1,
                                          const float * __restrict  pepsi1,
                                          const float * __restrict  pmur2,
                                          const float * __restrict  pmui2,
                                          const float * __restrict  pepsr2,
                                          const float * __restrict  pepsi2,
                                          float * __restrict  Rr,
                                          float * __restrict  Ri) {

                       register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                       register __m512 mur1  = _mm512_loadu_ps(&pmur1[0]);
                       register __m512 mui1  = _mm512_loadu_ps(&pmui1[0]);
                       register __m512 epsr1 = _mm512_loadu_ps(&pepsr1[0]);
                       register __m512 epsi1 = _mm512_loadu_ps(&pepsi1[0]);
                       register __m512 mur2  = _mm512_loadu_ps(&pmur2[0]);
                       register __m512 mui2  = _mm512_loadu_ps(&pmui2[0]);
                       register __m512 epsr2 = _mm512_loadu_ps(&pepsr2[0]);
                       register __m512 epsi2 = _mm512_loadu_ps(&pepsi2[0]);
                       register __m512 z1r,z1i,z2r,z2i;
                       register __m512 cost,numr,numi,denr,deni;
                       register __m512 resr,resi;
                       zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcosf(tht);
                       zi_f716_zmm16r4(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_fmsub_ps(z2r,cost,z1r);
                       denr = _mm512_fmadd_ps(z2r,cost,z2r);
                       numi = _mm512_fmsub_ps(z2i,cost,z1i);
                       deni = _mm512_fmadd_ps(z2i,cost,z2i);
                       cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
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
                   void R_f7130_zmm16r4(const __m512 tht,
                                        const __m512 mur1,
                                        const __m512 mui1,
                                        const __m512 epsr1,
                                        const __m512 epsi1,
                                        const __m512 mur2,
                                        const __m512 mui2,
                                        const __m512 epsr2,
                                        const __m512 epsi2,
                                        __m512 * __restrict Rr,
                                        __m512 * __restrict Ri) {

                       register __m512 z1r,z1i,z2r,z2i;
                       register __m512 cost,numr,numi,denr,deni;
                       register __m512 resr,resi,t0r,t0i;
                       zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcosf(tht);
                       t0r  = _mm512_mul_ps(z1r,cost);
                       t0i  = _mm512_mul_ps(z1i,cost);
                       zi_f716_zmm16r4(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_sub_ps(z2r,t0r);
                       denr = _mm512_add_ps(z2r,t0r);
                       numi = _mm512_sub_ps(z2i,t0i);
                       deni = _mm512_add_ps(z2i,t0i);
                       cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                       *Rr = resr;
                       *Ri = resi;
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7130_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur1,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui1,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr1,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi1,
                                          const float * __restrict __ATTR_ALIGN__(64) pmur2,
                                          const float * __restrict __ATTR_ALIGN__(64) pmui2,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsr2,
                                          const float * __restrict __ATTR_ALIGN__(64) pepsi2,
                                          float * __restrict __ATTR_ALIGN__(64) Rr,
                                          float * __restrict __ATTR_ALIGN__(64) Ri) {

                       register __m512 tht   = _mm512_load_ps(&ptht[0]);
                       register __m512 mur1  = _mm512_load_ps(&pmur1[0]);
                       register __m512 mui1  = _mm512_load_ps(&pmui1[0]);
                       register __m512 epsr1 = _mm512_load_ps(&pepsr1[0]);
                       register __m512 epsi1 = _mm512_load_ps(&pepsi1[0]);
                       register __m512 mur2  = _mm512_load_ps(&pmur2[0]);
                       register __m512 mui2  = _mm512_load_ps(&pmui2[0]);
                       register __m512 epsr2 = _mm512_load_ps(&pepsr2[0]);
                       register __m512 epsi2 = _mm512_load_ps(&pepsi2[0]);
                       register __m512 z1r,z1i,z2r,z2i;
                       register __m512 cost,numr,numi,denr,deni;
                       register __m512 resr,resi,t0r,t0i;
                       zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcosf(tht);
                       t0r  = _mm512_mul_ps(z1r,cost);
                       t0i  = _mm512_mul_ps(z1i,cost);
                       zi_f716_zmm16r4(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_sub_ps(z2r,t0r);
                       denr = _mm512_add_ps(z2r,t0r);
                       numi = _mm512_sub_ps(z2i,t0i);
                       deni = _mm512_add_ps(z2i,t0i);
                       cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                       _mm512_store_ps(&Rr[0] ,resr);
                       _mm512_store_ps(&Ri[0] ,resi);
               }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void R_f7130_zmm16r4_u(const float * __restrict  ptht,
                                          const float * __restrict  pmur1,
                                          const float * __restrict  pmui1,
                                          const float * __restrict  pepsr1,
                                          const float * __restrict  pepsi1,
                                          const float * __restrict  pmur2,
                                          const float * __restrict  pmui2,
                                          const float * __restrict  pepsr2,
                                          const float * __restrict  pepsi2,
                                          float * __restrict  Rr,
                                          float * __restrict  Ri) {

                       register __m512 tht   = _mm512_loadu_ps(&ptht[0]);
                       register __m512 mur1  = _mm512_loadu_ps(&pmur1[0]);
                       register __m512 mui1  = _mm512_loadu_ps(&pmui1[0]);
                       register __m512 epsr1 = _mm512_loadu_ps(&pepsr1[0]);
                       register __m512 epsi1 = _mm512_loadu_ps(&pepsi1[0]);
                       register __m512 mur2  = _mm512_loadu_ps(&pmur2[0]);
                       register __m512 mui2  = _mm512_loadu_ps(&pmui2[0]);
                       register __m512 epsr2 = _mm512_loadu_ps(&pepsr2[0]);
                       register __m512 epsi2 = _mm512_loadu_ps(&pepsi2[0]);
                       register __m512 z1r,z1i,z2r,z2i;
                       register __m512 cost,numr,numi,denr,deni;
                       register __m512 resr,resi,t0r,t0i;
                       zi_f716_zmm16r4(tht1,mur1,mui1,epsr1,epsi1,&z1r,&z1i);
                       cost = xcosf(tht);
                       t0r  = _mm512_mul_ps(z1r,cost);
                       t0i  = _mm512_mul_ps(z1i,cost);
                       zi_f716_zmm16r4(tht1,mur2,mui2,epsr2,epsi2,&z2r,&z2i);
                       numr = _mm512_sub_ps(z2r,t0r);
                       denr = _mm512_add_ps(z2r,t0r);
                       numi = _mm512_sub_ps(z2i,t0i);
                       deni = _mm512_add_ps(z2i,t0i);
                       cdiv_zmm16r4(numr,numi,denr,deni,&resr,&resi);
                       _mm512_storeu_ps(&Rr[0] ,resr);
                       _mm512_storeu_ps(&Ri[0] ,resi);
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
                   __m512 R_f7215_f7216_zmm16r4(const __m512 d,
                                                const __m512 k0,
                                                const __m512 alp,
                                                const __m512 tht) {

                         const __m512 vmin = _mm512_set1_ps(1.17549e-38);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __mmask16 m = _mm512_cmp_ps_mask(d,vmin,_CMP_EQ_OQ);
                         register __m512 pid2,cost,cos2t,num,den,x0,x1,x2;
                         register __m512 k,k01a,sin2t,k02k,sqr;
                         register __m512 R;
                         if(!m) {
                            pid2 = _mm512_mul_ps(pi,_mm512_mul_ps(d,hlf));
                            cost = xcosf(tht);
                            cos2t= _mm512_fmadd_ps(cost,cost,alp);
                            x0   = _mm512_sqrt_ps(cos2t);
                            x1   = _mm512_fmsub_ps(pid2,cost,x0);
                            num  = _mm512_sinh_ps(x1);
                            x2   = _mm512_fmadd_ps(pid2,cost,x0);
                            den  = _mm512_sinh_ps(x2);
                            R    = _mm512_div_ps(num,den);
                            return (R);
                        }
                        else {
                            const __m512 _1 = _mm512_sqrt_ps(1.0f);
                            k    = _mm512_sqrt_ps(_mm512_sub_ps(_1,alp));
                            cost = xcosf(tht);
                            k02k = _mm512_div_ps(_mm512_mul_ps(k0,k0),k);
                            sint = xsinf(tht);
                            sin2t= _mm512_mul_ps(sint,sint);
                            x0   = _mm512_sub_ps(_1,_mm512_mul_ps(k02k,sin2t));
                            sqr  = _mm512_sqrt_ps(x0);
                            x1   = _mm512_mul_ps(k,sqr);
                            num  = _mm512_fmsub_ps(k0,cost,x1);
                            den  = _mm512_fmadd_ps(k0,cost,x1);
                            R    = _mm512_div_ps(num,den);
                            return (R);
                        }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 R_f7215_f7216_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pd,
                                                  const float * __restrict __ATTR_ALIGN__(64) pk0,
                                                  const float * __restrict __ATTR_ALIGN__(64) palp,
                                                  const float * __restrict __ATTR_ALIGN__(64) ptht) {

                         register __m512 d   = _mm512_load_ps(&pd[0]);
                         register __m512 k0  = _mm512_load_ps(&pk0[0]);
                         register __m512 alp = _mm512_load_ps(&palp[0]);
                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         const __m512 vmin = _mm512_set1_ps(1.17549e-38);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __mmask16 m = _mm512_cmp_ps_mask(d,vmin,_CMP_EQ_OQ);
                         register __m512 pid2,cost,cos2t,num,den,x0,x1,x2;
                         register __m512 k,k01a,sin2t,k02k,sqr;
                         register __m512 R;
                         if(!m) {
                            pid2 = _mm512_mul_ps(pi,_mm512_mul_ps(d,hlf));
                            cost = xcosf(tht);
                            cos2t= _mm512_fmadd_ps(cost,cost,alp);
                            x0   = _mm512_sqrt_ps(cos2t);
                            x1   = _mm512_fmsub_ps(pid2,cost,x0);
                            num  = _mm512_sinh_ps(x1);
                            x2   = _mm512_fmadd_ps(pid2,cost,x0);
                            den  = _mm512_sinh_ps(x2);
                            R    = _mm512_div_ps(num,den);
                            return (R);
                        }
                        else {
                            const __m512 _1 = _mm512_sqrt_ps(1.0f);
                            k    = _mm512_sqrt_ps(_mm512_sub_ps(_1,alp));
                            cost = xcosf(tht);
                            k02k = _mm512_div_ps(_mm512_mul_ps(k0,k0),k);
                            sint = xsinf(tht);
                            sin2t= _mm512_mul_ps(sint,sint);
                            x0   = _mm512_sub_ps(_1,_mm512_mul_ps(k02k,sin2t));
                            sqr  = _mm512_sqrt_ps(x0);
                            x1   = _mm512_mul_ps(k,sqr);
                            num  = _mm512_fmsub_ps(k0,cost,x1);
                            den  = _mm512_fmadd_ps(k0,cost,x1);
                            R    = _mm512_div_ps(num,den);
                            return (R);
                        }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 R_f7215_f7216_zmm16r4_u(const float * __restrict  pd,
                                                  const float * __restrict  pk0,
                                                  const float * __restrict  palp,
                                                  const float * __restrict  ptht) {

                         register __m512 d   = _mm512_loadu_ps(&pd[0]);
                         register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                         register __m512 alp = _mm512_loadu_ps(&palp[0]);
                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         const __m512 vmin = _mm512_set1_ps(1.17549e-38);
                         const __m512 hlf  = _mm512_set1_ps(0.5f);
                         const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __mmask16 m = _mm512_cmp_ps_mask(d,vmin,_CMP_EQ_OQ);
                         register __m512 pid2,cost,cos2t,num,den,x0,x1,x2;
                         register __m512 k,k01a,sin2t,k02k,sqr;
                         register __m512 R;
                         if(!m) {
                            pid2 = _mm512_mul_ps(pi,_mm512_mul_ps(d,hlf));
                            cost = xcosf(tht);
                            cos2t= _mm512_fmadd_ps(cost,cost,alp);
                            x0   = _mm512_sqrt_ps(cos2t);
                            x1   = _mm512_fmsub_ps(pid2,cost,x0);
                            num  = _mm512_sinh_ps(x1);
                            x2   = _mm512_fmadd_ps(pid2,cost,x0);
                            den  = _mm512_sinh_ps(x2);
                            R    = _mm512_div_ps(num,den);
                            return (R);
                        }
                        else {
                            const __m512 _1 = _mm512_sqrt_ps(1.0f);
                            k    = _mm512_sqrt_ps(_mm512_sub_ps(_1,alp));
                            cost = xcosf(tht);
                            k02k = _mm512_div_ps(_mm512_mul_ps(k0,k0),k);
                            sint = xsinf(tht);
                            sin2t= _mm512_mul_ps(sint,sint);
                            x0   = _mm512_sub_ps(_1,_mm512_mul_ps(k02k,sin2t));
                            sqr  = _mm512_sqrt_ps(x0);
                            x1   = _mm512_mul_ps(k,sqr);
                            num  = _mm512_fmsub_ps(k0,cost,x1);
                            den  = _mm512_fmadd_ps(k0,cost,x1);
                            R    = _mm512_div_ps(num,den);
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
                   void Esz_f741_zmm16r4(const __m512 k0,
                                         const __m512 r,
                                         const __m512 a,
                                         const __m512 tht,
                                         const __m512 Eir,
                                         const __m512 Eii,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                        const __m512 _1   = _mm512_set1_ps(1.0f);
                        const __m512 gam  = _mm512_set1_ps(1.7811f);
                        const __m512 qtr  = _mm512_set1_ps(0.25f);
                        const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 pi2  = _mm512_set1_ps(0.5f*3.14159265358979323846264338328f);
                        const __m512 _4   = _mm512_set1_ps(4.0f);
                        const __m512 pi4  = _mm512_set1_ps(0.25f*3.14159265358979323846264338328f);
                        register __m512 denr,deni,ear,eai,cer,cei,arg;
                        register __m512 t0r,t0i,num,k02,a2,k0a,k0r,cost,trm;
                        register __m512 x0,x1,t1r,t1i;
                        deni = pi2;
                        cost = xcosf(tht);
                        k0r  = _mm512_mul_ps(k0,r);
                        a2   = _mm512_mul_ps(a,a);
                        ear  = _mm512_setzero_ps();
                        k0a  = _mm512_mul_ps(k0,a);
                        eai  = _mm512_add_ps(k0r,pi4);
                        k02  = _mm512_mul_ps(k0,k0);
                        x0   = _mm512_div_ps(pi,_mm512_add_ps(k0r,k0r));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        arg  = _mm512_div_ps(_4,_mm512_mul_ps(gam,k0a));
                        trm  = _mm512_sqrt_ps(x0);
                        x1   = _mm512_mul_ps(cost,cost);
                        denr = xlogf(arg);
                        x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                        num  = _mm512_mul_ps(x0,x1);
                        t0r  = _mm512_mul_ps(_mm512_div_ps(num,denr),trm);
                        t0i  = _mm512_mul_ps(_mm512_div_ps(num,deni),trm);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        cmul_zmm16r4(Eir,Eii,t1r,t1i,*Esr,*Esi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f741_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pEir,
                                           const float * __restrict __ATTR_ALIGN__(64) pEii,
                                           float * __restrict __ATTR_ALIGN__(64) Esr,
                                           float * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 a   = _mm512_load_ps(&pa[0]);
                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 Eir = _mm512_load_ps(&pEir[0]); 
                        register __m512 Eii = _mm512_load_ps(&pEii[0]);
                        const __m512 _1   = _mm512_set1_ps(1.0f);
                        const __m512 gam  = _mm512_set1_ps(1.7811f);
                        const __m512 qtr  = _mm512_set1_ps(0.25f);
                        const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 pi2  = _mm512_set1_ps(0.5f*3.14159265358979323846264338328f);
                        const __m512 _4   = _mm512_set1_ps(4.0f);
                        const __m512 pi4  = _mm512_set1_ps(0.25f*3.14159265358979323846264338328f);
                        register __m512 denr,deni,ear,eai,cer,cei,arg;
                        register __m512 t0r,t0i,num,k02,a2,k0a,k0r,cost,trm;
                        register __m512 x0,x1,t1r,t1i,resr,resi;
                        deni = pi2;
                        cost = xcosf(tht);
                        k0r  = _mm512_mul_ps(k0,r);
                        a2   = _mm512_mul_ps(a,a);
                        ear  = _mm512_setzero_ps();
                        k0a  = _mm512_mul_ps(k0,a);
                        eai  = _mm512_add_ps(k0r,pi4);
                        k02  = _mm512_mul_ps(k0,k0);
                        x0   = _mm512_div_ps(pi,_mm512_add_ps(k0r,k0r));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        arg  = _mm512_div_ps(_4,_mm512_mul_ps(gam,k0a));
                        trm  = _mm512_sqrt_ps(x0);
                        x1   = _mm512_mul_ps(cost,cost);
                        denr = xlogf(arg);
                        x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                        num  = _mm512_mul_ps(x0,x1);
                        t0r  = _mm512_mul_ps(_mm512_div_ps(num,denr),trm);
                        t0i  = _mm512_mul_ps(_mm512_div_ps(num,deni),trm);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        cmul_zmm16r4(Eir,Eii,t1r,t1i,&resr,&resi);
                        _mm512_store_ps(&Esr[0], resr);
                        _mm512_store_ps(&Esi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f741_zmm16r4_u(const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict  pa,
                                           const float * __restrict  ptht,
                                           const float * __restrict  pEir,
                                           const float * __restrict  pEii,
                                           float * __restrict  Esr,
                                           float * __restrict  Esi) {

                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 a   = _mm512_loadu_ps(&pa[0]);
                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 Eir = _mm512_loadu_ps(&pEir[0]); 
                        register __m512 Eii = _mm512_loadu_ps(&pEii[0]);
                        const __m512 _1   = _mm512_set1_ps(1.0f);
                        const __m512 gam  = _mm512_set1_ps(1.7811f);
                        const __m512 qtr  = _mm512_set1_ps(0.25f);
                        const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 pi2  = _mm512_set1_ps(0.5f*3.14159265358979323846264338328f);
                        const __m512 _4   = _mm512_set1_ps(4.0f);
                        const __m512 pi4  = _mm512_set1_ps(0.25f*3.14159265358979323846264338328f);
                        register __m512 denr,deni,ear,eai,cer,cei,arg;
                        register __m512 t0r,t0i,num,k02,a2,k0a,k0r,cost,trm;
                        register __m512 x0,x1,t1r,t1i,resr,resi;
                        deni = pi2;
                        cost = xcosf(tht);
                        k0r  = _mm512_mul_ps(k0,r);
                        a2   = _mm512_mul_ps(a,a);
                        ear  = _mm512_setzero_ps();
                        k0a  = _mm512_mul_ps(k0,a);
                        eai  = _mm512_add_ps(k0r,pi4);
                        k02  = _mm512_mul_ps(k0,k0);
                        x0   = _mm512_div_ps(pi,_mm512_add_ps(k0r,k0r));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        arg  = _mm512_div_ps(_4,_mm512_mul_ps(gam,k0a));
                        trm  = _mm512_sqrt_ps(x0);
                        x1   = _mm512_mul_ps(cost,cost);
                        denr = xlogf(arg);
                        x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                        num  = _mm512_mul_ps(x0,x1);
                        t0r  = _mm512_mul_ps(_mm512_div_ps(num,denr),trm);
                        t0i  = _mm512_mul_ps(_mm512_div_ps(num,deni),trm);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        cmul_zmm16r4(Eir,Eii,t1r,t1i,&resr,&resi);
                        _mm512_storeu_ps(&Esr[0], resr);
                        _mm512_storeu_ps(&Esi[0], resi);
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
                   void Hsz_f742_zmm16r4(const __m512 k0a,
                                         const __m512 k0r,
                                         const __m512 tht,
                                         const __m512 Hir,
                                         const __m512 Hii,
                                         __m512 * __restrict Hsr,
                                         __m512 * __restrict Hsi) {

                         const __m512 _1o8 = _mm512_set1_ps(0.125f);
                         const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 pi4  = _mm512_set1_ps(0.25f*3.14159265358979323846264338328f);
                         register __m512 ear,eai,cer,cei;
                         register __m512 trm,t0r,t0i,num,cost,x0,x1;
                         ear  = _mm512_setzero_ps();
                         cost = xcosf(tht);
                         eai  = _mm512_add_ps(k0r,pi4);
                         x0   = _mm512_div_ps(pi,_mm512_add_ps(k0r,k0r));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         trm  = _mm512_sqrt_ps(x0);
                         x1   = _mm512_add_ps(k0a,k0a);
                         x0   = _mm512_mul_ps(cost,cost);
                         num  = _mm512_mul_ps(_mm512_mul_ps(x1,x1),x0);
                         t0r  = _mm512_mul_ps(trm,cer);
                         t0i  = _mm512_mul_ps(trm,cei);
                         num  = _mm512_mul_ps(_1o8,num);
                         x0   = _mm512_mul_ps(Hsr,num);
                         x1   = _mm512_mul_ps(Hsi,num);
                         cmul_zmm16r4(x0,x1,t0r,t0i,*Hsr,*Hsi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f742_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pHir,
                                           const float * __restrict __ATTR_ALIGN__(64) pHii,
                                           float * __restrict __ATTR_ALIGN__(64) Hsr,
                                           float * __restrict __ATTR_ALIGN__(64) Hsi) {

                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         register __m512 k0r = _mm512_load_ps(&pk0r[0]);
                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         register __m512 Hir = _mm512_load_ps(&pHir[0]);
                         register __m512 Hii = _mm512_load_ps(&pHii[0]);
                         const __m512 _1o8 = _mm512_set1_ps(0.125f);
                         const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 pi4  = _mm512_set1_ps(0.25f*3.14159265358979323846264338328f);
                         register __m512 ear,eai,cer,cei;
                         register __m512 trm,t0r,t0i,num,cost,x0,x1,resr,resi;
                         ear  = _mm512_setzero_ps();
                         cost = xcosf(tht);
                         eai  = _mm512_add_ps(k0r,pi4);
                         x0   = _mm512_div_ps(pi,_mm512_add_ps(k0r,k0r));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         trm  = _mm512_sqrt_ps(x0);
                         x1   = _mm512_add_ps(k0a,k0a);
                         x0   = _mm512_mul_ps(cost,cost);
                         num  = _mm512_mul_ps(_mm512_mul_ps(x1,x1),x0);
                         t0r  = _mm512_mul_ps(trm,cer);
                         t0i  = _mm512_mul_ps(trm,cei);
                         num  = _mm512_mul_ps(_1o8,num);
                         x0   = _mm512_mul_ps(Hsr,num);
                         x1   = _mm512_mul_ps(Hsi,num);
                         cmul_zmm16r4(x0,x1,t0r,t0i,&resr,&resi);
                         _mm512_store_ps(&Hsr[0], resr);
                         _mm512_store_ps(&Hsi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f742_zmm16r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  pk0r,
                                           const float * __restrict  ptht,
                                           const float * __restrict  pHir,
                                           const float * __restrict  pHii,
                                           float * __restrict  Hsr,
                                           float * __restrict  Hsi) {

                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 k0r = _mm512_loadu_ps(&pk0r[0]);
                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         register __m512 Hir = _mm512_loadu_ps(&pHir[0]);
                         register __m512 Hii = _mm512_loadu_ps(&pHii[0]);
                         const __m512 _1o8 = _mm512_set1_ps(0.125f);
                         const __m512 pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 pi4  = _mm512_set1_ps(0.25f*3.14159265358979323846264338328f);
                         register __m512 ear,eai,cer,cei,resr,resi;
                         register __m512 trm,t0r,t0i,num,cost,x0,x1;
                         ear  = _mm512_setzero_ps();
                         cost = xcosf(tht);
                         eai  = _mm512_add_ps(k0r,pi4);
                         x0   = _mm512_div_ps(pi,_mm512_add_ps(k0r,k0r));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         trm  = _mm512_sqrt_ps(x0);
                         x1   = _mm512_add_ps(k0a,k0a);
                         x0   = _mm512_mul_ps(cost,cost);
                         num  = _mm512_mul_ps(_mm512_mul_ps(x1,x1),x0);
                         t0r  = _mm512_mul_ps(trm,cer);
                         t0i  = _mm512_mul_ps(trm,cei);
                         num  = _mm512_mul_ps(_1o8,num);
                         x0   = _mm512_mul_ps(Hsr,num);
                         x1   = _mm512_mul_ps(Hsi,num);
                         cmul_zmm16r4(x0,x1,t0r,t0i,&resr,&resi);
                         _mm512_storeu_ps(&Hsr[0], resr);
                         _mm512_storeu_ps(&Hsi[0], resi);
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
                   __m512 rcs_f743_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht) {

                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 pisq= _mm512_set1_ps(0.25f*9.869604401089358618834490999876f);
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 c0  = _mm512_set1_ps(4.48f);
                          register __m512 fac,num,den,cost,k02,a2,k0a;
                          register __m512 rcs,x0,x1,arg,larg,rat;
                          k0a  = _mm512_mul_ps(k0,a);
                          fac  = _mm512_div_ps(pis,k0);
                          k02  = _mm512_mul_ps(k0,k0);
                          a2   = _mm512_mul_ps(a,a);
                          cost = xcosf(tht);
                          x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                          arg  = _mm512_div_ps(c0,_mm512_add_ps(k0a,k0a));
                          x1   = _mm512_mul_ps(x0,_mm512_mul_ps(cost,cost));
                          larg = xlogf(arg);
                          num  = _mm512_mul_ps(x1,x1);
                          den  = _mm512_fmadd_ps(larg,larg,pisq);
                          rat  = _mm512_div_ps(num,den);
                          rcs  = _mm512mul_ps(fac,rat);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f743_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0 = _mm512_load_ps(&pk0[0]);
                          register __m512 a  = _mm512_load_ps(&pa[0]);
                          register __m512 tht= _mm512_load_ps(&ptht[0]);
                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 pisq= _mm512_set1_ps(0.25f*9.869604401089358618834490999876f);
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 c0  = _mm512_set1_ps(4.48f);
                          register __m512 fac,num,den,cost,k02,a2,k0a;
                          register __m512 rcs,x0,x1,arg,larg,rat;
                          k0a  = _mm512_mul_ps(k0,a);
                          fac  = _mm512_div_ps(pis,k0);
                          k02  = _mm512_mul_ps(k0,k0);
                          a2   = _mm512_mul_ps(a,a);
                          cost = xcosf(tht);
                          x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                          arg  = _mm512_div_ps(c0,_mm512_add_ps(k0a,k0a));
                          x1   = _mm512_mul_ps(x0,_mm512_mul_ps(cost,cost));
                          larg = xlogf(arg);
                          num  = _mm512_mul_ps(x1,x1);
                          den  = _mm512_fmadd_ps(larg,larg,pisq);
                          rat  = _mm512_div_ps(num,den);
                          rcs  = _mm512_mul_ps(fac,rat);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f743_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ptht) {

                          register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a  = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht= _mm512_loadu_ps(&ptht[0]);
                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 pisq= _mm512_set1_ps(0.25f*9.869604401089358618834490999876f);
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 c0  = _mm512_set1_ps(4.48f);
                          register __m512 fac,num,den,cost,k02,a2,k0a;
                          register __m512 rcs,x0,x1,arg,larg,rat;
                          k0a  = _mm512_mul_ps(k0,a);
                          fac  = _mm512_div_ps(pis,k0);
                          k02  = _mm512_mul_ps(k0,k0);
                          a2   = _mm512_mul_ps(a,a);
                          cost = xcosf(tht);
                          x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                          arg  = _mm512_div_ps(c0,_mm512_add_ps(k0a,k0a));
                          x1   = _mm512_mul_ps(x0,_mm512_mul_ps(cost,cost));
                          larg = xlogf(arg);
                          num  = _mm512_mul_ps(x1,x1);
                          den  = _mm512_fmadd_ps(larg,larg,pisq);
                          rat  = _mm512_div_ps(num,den);
                          rcs  = _mm512_mul_ps(fac,rat);
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
                   __m512 rcs_f743_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht) {

                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f); 
                          const __m512 1o64= _mm512_set1_ps(0.015625f);
                          register __m512 k0a,k0a2,cost,cos2t,fac;
                          register __m512 rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          k0a2 = _mm512_add_ps(k0a,k0a);
                          fac  = _mm512_div_ps(pis,k0);
                          x2   = _mm512_mul_ps(cos2t,cos2t);
                          x0   = _mm512_mul_ps(k0a2,k0a2);
                          x1   = _mm512_mul_ps(x0,x0);
                          num  = _mm512_mul_ps(_mm512_mul_ps(x1,x2),1o64);
                          rcs  = _mm512_mul_ps(fac,num);
                          return (rcs);
                }



      } // radiolocation


} // gms


























#endif /*__GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__*/