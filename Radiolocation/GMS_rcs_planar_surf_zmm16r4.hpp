

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
                   __m512 rcs_f744_zmm16r4(const __m512 k0,
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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f744_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f744_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ptht) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
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
                   __m512 rcs_f745_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht,
                                           const __m512 tht2) {

                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 pisq= _mm512_set1_ps(0.25f*9.869604401089358618834490999876f);
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 c0  = _mm512_set1_ps(4.48f);
                          register __m512 fac,num,den,cost,k02,a2,k0a;
                          register __m512 rcs,x0,x1,arg,larg,rat,cost2;
                          k0a  = _mm512_mul_ps(k0,a);
                          fac  = _mm512_div_ps(pis,k0);
                          k02  = _mm512_mul_ps(k0,k0);
                          a2   = _mm512_mul_ps(a,a);
                          cost = xcosf(tht);
                          x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                          cost2= xcosf(tht2);
                          arg  = _mm512_div_ps(c0,_mm512_add_ps(k0a,k0a));
                          x1   = _mm512_mul_ps(x0,_mm512_mul_ps(cost,cost2));
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
                   __m512 rcs_f745_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht2) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          register __m512 tht2= _mm512_load_ps(&ptht2[0]);  
                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 pisq= _mm512_set1_ps(0.25f*9.869604401089358618834490999876f);
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 c0  = _mm512_set1_ps(4.48f);
                          register __m512 fac,num,den,cost,k02,a2,k0a;
                          register __m512 rcs,x0,x1,arg,larg,rat,cost2;
                          k0a  = _mm512_mul_ps(k0,a);
                          fac  = _mm512_div_ps(pis,k0);
                          k02  = _mm512_mul_ps(k0,k0);
                          a2   = _mm512_mul_ps(a,a);
                          cost = xcosf(tht);
                          x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                          cost2= xcosf(tht2);
                          arg  = _mm512_div_ps(c0,_mm512_add_ps(k0a,k0a));
                          x1   = _mm512_mul_ps(x0,_mm512_mul_ps(cost,cost2));
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
                   __m512 rcs_f745_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict pa,
                                             const float * __restrict  ptht,
                                             const float * __restrict  ptht2) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          register __m512 tht2= _mm512_loadu_ps(&ptht2[0]);  
                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f);
                          const __m512 pisq= _mm512_set1_ps(0.25f*9.869604401089358618834490999876f);
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          const __m512 qtr = _mm512_set1_ps(0.25f);
                          const __m512 c0  = _mm512_set1_ps(4.48f);
                          register __m512 fac,num,den,cost,k02,a2,k0a;
                          register __m512 rcs,x0,x1,arg,larg,rat,cost2;
                          k0a  = _mm512_mul_ps(k0,a);
                          fac  = _mm512_div_ps(pis,k0);
                          k02  = _mm512_mul_ps(k0,k0);
                          a2   = _mm512_mul_ps(a,a);
                          cost = xcosf(tht);
                          x0   = _mm512_fmadd_ps(k02,_mm512_mul_ps(a2,qtr),_1);
                          cost2= xcosf(tht2);
                          arg  = _mm512_div_ps(c0,_mm512_add_ps(k0a,k0a));
                          x1   = _mm512_mul_ps(x0,_mm512_mul_ps(cost,cost2));
                          larg = xlogf(arg);
                          num  = _mm512_mul_ps(x1,x1);
                          den  = _mm512_fmadd_ps(larg,larg,pisq);
                          rat  = _mm512_div_ps(num,den);
                          rcs  = _mm512mul_ps(fac,rat);
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
                   __m512 rcs_f746_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht,
                                           const __m512 tht2) {

                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f); 
                          const __m512 1o64= _mm512_set1_ps(0.015625f);
                          register __m512 k0a,k0a2,cost,cos2t,fac,cost2;
                          register __m512 rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          k0a2 = _mm512_add_ps(k0a,k0a);
                          cost2= xcosf(tht2);
                          fac  = _mm512_div_ps(pis,k0);
                          x2   = _mm512_mul_ps(_mm512_mul_ps(cost2.cost2),
                                 _mm512_mul_ps(cos2t,cos2t));
                          x0   = _mm512_mul_ps(k0a2,k0a2);
                          x1   = _mm512_mul_ps(x0,x0);
                          num  = _mm512_mul_ps(_mm512_mul_ps(x1,x2),1o64);
                          rcs  = _mm512_mul_ps(fac,num);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f746_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht2) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          register __m512 tht2= _mm512_load_ps(&ptht2[0]);
                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f); 
                          const __m512 1o64= _mm512_set1_ps(0.015625f);
                          register __m512 k0a,k0a2,cost,cos2t,fac,cost2;
                          register __m512 rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          k0a2 = _mm512_add_ps(k0a,k0a);
                          cost2= xcosf(tht2);
                          fac  = _mm512_div_ps(pis,k0);
                          x2   = _mm512_mul_ps(_mm512_mul_ps(cost2.cost2),
                                 _mm512_mul_ps(cos2t,cos2t));
                          x0   = _mm512_mul_ps(k0a2,k0a2);
                          x1   = _mm512_mul_ps(x0,x0);
                          num  = _mm512_mul_ps(_mm512_mul_ps(x1,x2),1o64);
                          rcs  = _mm512_mul_ps(fac,num);
                          return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f746_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ptht,
                                             const float * __restrict  ptht2) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          register __m512 tht2= _mm512_loadu_ps(&ptht2[0]);
                          const __m512 pis = _mm512_set1_ps(9.869604401089358618834490999876f); 
                          const __m512 1o64= _mm512_set1_ps(0.015625f);
                          register __m512 k0a,k0a2,cost,cos2t,fac,cost2;
                          register __m512 rcs,x0,x1,num,x2;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          k0a2 = _mm512_add_ps(k0a,k0a);
                          cost2= xcosf(tht2);
                          fac  = _mm512_div_ps(pis,k0);
                          x2   = _mm512_mul_ps(_mm512_mul_ps(cost2.cost2),
                                 _mm512_mul_ps(cos2t,cos2t));
                          x0   = _mm512_mul_ps(k0a2,k0a2);
                          x1   = _mm512_mul_ps(x0,x0);
                          num  = _mm512_mul_ps(_mm512_mul_ps(x1,x2),1o64);
                          rcs  = _mm512_mul_ps(fac,num);
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
                   __m512 rcs_f747_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht) {

                          register __m512 invk0,k0a,cost,sint,arg,sarg;
                          register __m512 rcs,num,sqr,x0;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          invk0= _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          arg  = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          sarg = xsinf(arg);
                          num  = _mm512_mul_ps(cost,sarg);
                          x0   = _mm512_div_ps(num,sint);
                          sqr  = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(invk0,sqr);
                          return (rcs); 
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f747_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          register __m512 invk0,k0a,cost,sint,arg,sarg;
                          register __m512 rcs,num,sqr,x0;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          invk0= _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          arg  = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          sarg = xsinf(arg);
                          num  = _mm512_mul_ps(cost,sarg);
                          x0   = _mm512_div_ps(num,sint);
                          sqr  = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(invk0,sqr);
                          return (rcs); 
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f747_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ptht) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          register __m512 invk0,k0a,cost,sint,arg,sarg;
                          register __m512 rcs,num,sqr,x0;
                          k0a  = _mm512_mul_ps(k0,a);
                          cost = xcosf(tht);
                          invk0= _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          arg  = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          sarg = xsinf(arg);
                          num  = _mm512_mul_ps(cost,sarg);
                          x0   = _mm512_div_ps(num,sint);
                          sqr  = _mm512_mul_ps(x0,x0);
                          rcs  = _mm512_mul_ps(invk0,sqr);
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
                   void coefg12_f7415_zmm16r4(  const __m512 k0a,
                                                const __m512 tht,
                                                __m512 * __restrict gamm1,
                                                __m512 * __restrict gamm2){
                                                 
                          const __m512 C0318309886183790671537767526745 = _mm512_set1_ps(0.318309886183790671537767526745f);
                          const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 C05   = _mm512_set1_ps(0.5f);
                          register __m512 thth,arg1,arg2,carg1,carg2,sqr,x0;
                          x0   = _mm512_add_ps(k0a,k0a);
                          thth = _mm512_mul_ps(C05,tht);
                          sqr  = _mm512_sqrt_ps(_mm512_mul_ps(x0,C0318309886183790671537767526745));
                          arg1 = _mm512_add_ps(C078539816339744830961566084582,thth);
                          carg1= xcosf(arg1);
                          x0   = _mm512_add_ps(sqr,sqr);
                          arg2 = _mm512_sub_ps(C078539816339744830961566084582,thth);
                          carg2= xcosf(arg2);
                          *gamm1 = _mm512_mul_ps(x0,_mm512_abs_ps(carg1));
                          *gamm2 = _mm512_mul_ps(x0,_mm512_abs_ps(carg2));
                }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefg12_f7415_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                                const float * __restrict __ATTR_ALIGN__(64) ptht,
                                                float * __restrict __ATTR_ALIGN__(64) gamm1,
                                                float * __restrict __ATTR_ALIGN__(64) gamm2){
                                
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);                 
                          const __m512 C0318309886183790671537767526745 = _mm512_set1_ps(0.318309886183790671537767526745f);
                          const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 C05   = _mm512_set1_ps(0.5f);
                          register __m512 thth,arg1,arg2,carg1,carg2,sqr,x0;
                          x0   = _mm512_add_ps(k0a,k0a);
                          thth = _mm512_mul_ps(C05,tht);
                          sqr  = _mm512_sqrt_ps(_mm512_mul_ps(x0,C0318309886183790671537767526745));
                          arg1 = _mm512_add_ps(C078539816339744830961566084582,thth);
                          carg1= xcosf(arg1);
                          x0   = _mm512_add_ps(sqr,sqr);
                          arg2 = _mm512_sub_ps(C078539816339744830961566084582,thth);
                          carg2= xcosf(arg2);
                          _mm512_store_ps(&gamm1[0] ,_mm512_mul_ps(x0,_mm512_abs_ps(carg1)));
                          _mm512_store_ps(&gamm2[0] ,_mm512_mul_ps(x0,_mm512_abs_ps(carg2)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefg12_f7415_zmm16r4_u(  const float * __restrict  pk0a,
                                                const float * __restrict  ptht,
                                                float * __restrict  gamm1,
                                                float * __restrict  gamm2){
                                
                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);                 
                          const __m512 C0318309886183790671537767526745 = _mm512_set1_ps(0.318309886183790671537767526745f);
                          const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                          const __m512 C05   = _mm512_set1_ps(0.5f);
                          register __m512 thth,arg1,arg2,carg1,carg2,sqr,x0;
                          x0   = _mm512_add_ps(k0a,k0a);
                          thth = _mm512_mul_ps(C05,tht);
                          sqr  = _mm512_sqrt_ps(_mm512_mul_ps(x0,C0318309886183790671537767526745));
                          arg1 = _mm512_add_ps(C078539816339744830961566084582,thth);
                          carg1= xcosf(arg1);
                          x0   = _mm512_add_ps(sqr,sqr);
                          arg2 = _mm512_sub_ps(C078539816339744830961566084582,thth);
                          carg2= xcosf(arg2);
                          _mm512_storeu_ps(&gamm1[0] ,_mm512_mul_ps(x0,_mm512_abs_ps(carg1)));
                          _mm512_storeu_ps(&gamm2[0] ,_mm512_mul_ps(x0,_mm512_abs_ps(carg2)));
                }



                      /*
                       Backscattered fields from the edges of strips.
                       Helper function for the formula 7.4-9
                       Electric-field (over z).
                       Formula 7.4-13
                  */

#include "GMS_rcs_common_zmm16r4.hpp"

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefA12_f7413_zmm16r4(const __m512 k0a,
                                              const __m512 tht,
                                              __m512 * __restrict A1r,
                                              __m512 * __restrict A1i,
                                              __m512 * __restrict A2r,
                                              __m512 * __restrict A2i) {

                        const __m512 C078539816339744830961566084582  = _mm512_set1_ps(-0.78539816339744830961566084582f);
                        const __m512 C141421356237309504880168872421  = _mm512_set1_ps(1.41421356237309504880168872421f);
                        register __m512 ear,eai,cer,cei,Cr1,Si1,Cr2,Si2;
                        register __m512 gam1,gam2;
                        ear = _mm512_setzero_ps();
                        eai = C078539816339744830961566084582;
                        coefg12_f7415_zmm16r4(k0a,tht,&gam1,&gam2); 
                        Cr1 = fresnel_C_zmm16r4(gam1);
                        Si1 = fresnel_S_zmm16r4(gam1);
                        Cr2 = fresnel_C_zmm16r4(gam2);
                        Si2 = fresnel_S_zmm16r4(gam2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cer = _mm512_mul_ps(C141421356237309504880168872421,cer);
                        cei = _mm512_mul_ps(C141421356237309504880168872421,cei);
                        cmul_zmm16r4(cer,cei,Cr1,Si1,*A1r,*A1i);
                        cmul_zmm16r4(cer,cei,Cr2,Si2,*A2r,*A2i);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefA12_f7413_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht,
                                              float * __restrict __ATTR_ALIGN__(64) A1r,
                                              float * __restrict __ATTR_ALIGN__(64) A1i,
                                              float * __restrict __ATTR_ALIGN__(64) A2r,
                                              float * __restrict __ATTR_ALIGN__(64) A2i) {

                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        const __m512 C078539816339744830961566084582  = _mm512_set1_ps(-0.78539816339744830961566084582f);
                        const __m512 C141421356237309504880168872421  = _mm512_set1_ps(1.41421356237309504880168872421f);
                        register __m512 ear,eai,cer,cei,Cr1,Si1,Cr2,Si2;
                        register __m512 gam1,gam2,res1r,res1i,res2r,res2i;
                        ear = _mm512_setzero_ps();
                        eai = C078539816339744830961566084582;
                        coefg12_f7415_zmm16r4(k0a,tht,&gam1,&gam2); 
                        Cr1 = fresnel_C_zmm16r4(gam1);
                        Si1 = fresnel_S_zmm16r4(gam1);
                        Cr2 = fresnel_C_zmm16r4(gam2);
                        Si2 = fresnel_S_zmm16r4(gam2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cer = _mm512_mul_ps(C141421356237309504880168872421,cer);
                        cei = _mm512_mul_ps(C141421356237309504880168872421,cei);
                        cmul_zmm16r4(cer,cei,Cr1,Si1,&res1r,&res1i);
                        __m512_store_ps(&A1r[0], res1r);
                        __m512_store_ps(&A1i[0], res1i);
                        cmul_zmm16r4(cer,cei,Cr2,Si2,&res2r,&res2i);
                        __m512_store_ps(&A2r[0], res2r);
                        __m512_store_ps(&A2i[0], res2i);
               }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefA12_f7413_zmm16r4_u(const float * __restrict  pk0a,
                                                const float * __restrict  ptht,
                                                float * __restrict  A1r,
                                                float * __restrict  A1i,
                                                float * __restrict  A2r,
                                                float * __restrict  A2i) {

                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        const __m512 C078539816339744830961566084582  = _mm512_set1_ps(-0.78539816339744830961566084582f);
                        const __m512 C141421356237309504880168872421  = _mm512_set1_ps(1.41421356237309504880168872421f);
                        register __m512 ear,eai,cer,cei,Cr1,Si1,Cr2,Si2;
                        register __m512 gam1,gam2,res1r,res1i,res2r,res2i;
                        ear = _mm512_setzero_ps();
                        eai = C078539816339744830961566084582;
                        coefg12_f7415_zmm16r4(k0a,tht,&gam1,&gam2); 
                        Cr1 = fresnel_C_zmm16r4(gam1);
                        Si1 = fresnel_S_zmm16r4(gam1);
                        Cr2 = fresnel_C_zmm16r4(gam2);
                        Si2 = fresnel_S_zmm16r4(gam2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cer = _mm512_mul_ps(C141421356237309504880168872421,cer);
                        cei = _mm512_mul_ps(C141421356237309504880168872421,cei);
                        cmul_zmm16r4(cer,cei,Cr1,Si1,&res1r,&res1i);
                        __m512_storeu_ps(&A1r[0], res1r);
                        __m512_storeu_ps(&A1i[0], res1i);
                        cmul_zmm16r4(cer,cei,Cr2,Si2,&res2r,&res2i);
                        __m512_storeu_ps(&A2r[0], res2r);
                        __m512_storeu_ps(&A2i[0], res2i);
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
                   void coefB12_f7414_zmm16r4(const __m512 k0a,
                                              const __m512 tht,
                                              __m512 * __restrict B1r,
                                              __m512 * __restrict B1i,
                                              __m512 * __restrict B2r,
                                              __m512 * __restrict B2i) {

                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(-0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C05                              = _mm512_set1_ps(0.5f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear1,eai1,eai2,cer1,cei1,cer2,cei2;
                         register __m512 ir,ii,sint,htht,carg1,carg2,arg1,arg2,abs1,abs2;
                         register __m512 x0,x1,x2,x3,x4,x5,t0r,t0i,t1r,t1i,k0a2;
                         register __m512 A1r,A1i,A2r,A2i;
                         ir   = _mm512_setzero_ps();
                         x0   = _mm512_sqrt_ps(_mm512_mul_ps(C314159265358979323846264338328,k0a));
                         coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         htht = _mm512_mul_ps(C05,tht);
                         ear1 = ir;
                         ii   = _mm512_add_ps(x0,x0);
                         k0a2 = _mm512_add_ps(k0a,k0a);
                         sint = xsinf(tht);
                         x0   = _mm512_add_ps(C10,sint);
                         eai2 = _mm512_fmsub_ps(k0a2,x0,C078539816339744830961566084582);
                         cexp_zmm16r4(ear1,eai2,&cer1,&cei1);
                         arg1 = _mm512_sub_ps(C078539816339744830961566084582,htht);
                         carg1= xcosf(arg1);
                         x1   = _mm512_sub_ps(C10,sint);
                         eai1 = _mm512_fmsub_ps(k0a2,x1,C078539816339744830961566084582);
                         cexp_zmm16r4(ear1,eai1,&cer2,&cei2);
                         abs1 = _mm512_abs_ps(carg1);
                         t0r  = _mm512_div_ps(cer2,abs1);
                         arg2 = _mm512_add_ps(C078539816339744830961566084582,htht);
                         t0i  = _mm512_div_ps(cei2,abs1);
                         carg2= xcosf(arg2);
                         abs2 = _mm512_abs_ps(carg2);
                         t1r  = _mm512_div_ps(cer1,abs2);
                         cmul_zmm16r4(ir,ii,t0r,t0i,&x2,&x3);
                         t1i  = _mm512_div_ps(cei1,abs2);
                         cmul_zmm16r4(ir,ii,t1r,t1i,&x4,&x5);
                         *B1r  = _mm512_add_ps(A1r,x2);
                         *B2r  = _mm512_add_ps(A2r,x4);
                         *B1i  = _mm512_add_ps(A1i,x3);
                         *B2i  = _mm512_add_ps(A2i,x5);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefB12_f7414_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                                const float * __restrict __ATTR_ALIGN__(64)  ptht,
                                                float * __restrict __ATTR_ALIGN__(64) B1r,
                                                float * __restrict __ATTR_ALIGN__(64) B1i,
                                                float * __restrict __ATTR_ALIGN__(64) B2r,
                                                float * __restrict __ATTR_ALIGN__(64) B2i) {

                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(-0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C05                              = _mm512_set1_ps(0.5f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear1,eai1,eai2,cer1,cei1,cer2,cei2;
                         register __m512 ir,ii,sint,htht,carg1,carg2,arg1,arg2,abs1,abs2;
                         register __m512 x0,x1,x2,x3,x4,x5,t0r,t0i,t1r,t1i,k0a2;
                         register __m512 A1r,A1i,A2r,A2i;
                         ir   = _mm512_setzero_ps();
                         x0   = _mm512_sqrt_ps(_mm512_mul_ps(C314159265358979323846264338328,k0a));
                         coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         htht = _mm512_mul_ps(C05,tht);
                         ear1 = ir;
                         ii   = _mm512_add_ps(x0,x0);
                         k0a2 = _mm512_add_ps(k0a,k0a);
                         sint = xsinf(tht);
                         x0   = _mm512_add_ps(C10,sint);
                         eai2 = _mm512_fmsub_ps(k0a2,x0,C078539816339744830961566084582);
                         cexp_zmm16r4(ear1,eai2,&cer1,&cei1);
                         arg1 = _mm512_sub_ps(C078539816339744830961566084582,htht);
                         carg1= xcosf(arg1);
                         x1   = _mm512_sub_ps(C10,sint);
                         eai1 = _mm512_fmsub_ps(k0a2,x1,C078539816339744830961566084582);
                         cexp_zmm16r4(ear1,eai1,&cer2,&cei2);
                         abs1 = _mm512_abs_ps(carg1);
                         t0r  = _mm512_div_ps(cer2,abs1);
                         arg2 = _mm512_add_ps(C078539816339744830961566084582,htht);
                         t0i  = _mm512_div_ps(cei2,abs1);
                         carg2= xcosf(arg2);
                         abs2 = _mm512_abs_ps(carg2);
                         t1r  = _mm512_div_ps(cer1,abs2);
                         cmul_zmm16r4(ir,ii,t0r,t0i,&x2,&x3);
                         t1i  = _mm512_div_ps(cei1,abs2);
                         cmul_zmm16r4(ir,ii,t1r,t1i,&x4,&x5);
                         _mm512_store_ps(&B1r[0], _mm512_add_ps(A1r,x2));
                         _mm512_store_ps(&B2r[0], _mm512_add_ps(A2r,x4));
                         _mm512_store_ps(&B1i[0], _mm512_add_ps(A1i,x3));
                         _mm512_store_ps(&B2i[0], _mm512_add_ps(A2i,x5));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void coefB12_f7414_zmm16r4_u(const float * __restrict   pk0a,
                                                const float * __restrict   ptht,
                                                float * __restrict  B1r,
                                                float * __restrict  B1i,
                                                float * __restrict  B2r,
                                                float * __restrict  B2i) {

                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(-0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C05                              = _mm512_set1_ps(0.5f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear1,eai1,eai2,cer1,cei1,cer2,cei2;
                         register __m512 ir,ii,sint,htht,carg1,carg2,arg1,arg2,abs1,abs2;
                         register __m512 x0,x1,x2,x3,x4,x5,t0r,t0i,t1r,t1i,k0a2;
                         register __m512 A1r,A1i,A2r,A2i;
                         ir   = _mm512_setzero_ps();
                         x0   = _mm512_sqrt_ps(_mm512_mul_ps(C314159265358979323846264338328,k0a));
                         coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         htht = _mm512_mul_ps(C05,tht);
                         ear1 = ir;
                         ii   = _mm512_add_ps(x0,x0);
                         k0a2 = _mm512_add_ps(k0a,k0a);
                         sint = xsinf(tht);
                         x0   = _mm512_add_ps(C10,sint);
                         eai2 = _mm512_fmsub_ps(k0a2,x0,C078539816339744830961566084582);
                         cexp_zmm16r4(ear1,eai2,&cer1,&cei1);
                         arg1 = _mm512_sub_ps(C078539816339744830961566084582,htht);
                         carg1= xcosf(arg1);
                         x1   = _mm512_sub_ps(C10,sint);
                         eai1 = _mm512_fmsub_ps(k0a2,x1,C078539816339744830961566084582);
                         cexp_zmm16r4(ear1,eai1,&cer2,&cei2);
                         abs1 = _mm512_abs_ps(carg1);
                         t0r  = _mm512_div_ps(cer2,abs1);
                         arg2 = _mm512_add_ps(C078539816339744830961566084582,htht);
                         t0i  = _mm512_div_ps(cei2,abs1);
                         carg2= xcosf(arg2);
                         abs2 = _mm512_abs_ps(carg2);
                         t1r  = _mm512_div_ps(cer1,abs2);
                         cmul_zmm16r4(ir,ii,t0r,t0i,&x2,&x3);
                         t1i  = _mm512_div_ps(cei1,abs2);
                         cmul_zmm16r4(ir,ii,t1r,t1i,&x4,&x5);
                         _mm512_storeu_ps(&B1r[0], _mm512_add_ps(A1r,x2));
                         _mm512_storeu_ps(&B2r[0], _mm512_add_ps(A2r,x4));
                         _mm512_storeu_ps(&B1i[0], _mm512_add_ps(A1i,x3));
                         _mm512_storeu_ps(&B2i[0], _mm512_add_ps(A2i,x5));
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
                   void Esz_f749_zmm16r4(const __m512 tht,
                                         const __m512 k0a,
                                         const __m512 k0r,
                                         const __m512 Eir,
                                         const __m512 Eii,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512 B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                         register __m512 ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         register __m512 sqr,x0,x1,k0a2,y0;
                         ear   = _mm512_setzero_ps();
                         sqr   = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_ps(k0a,k0a);
                         eai   = _mm512_add_ps(k0r,C078539816339744830961566084582);
                         sint  = xsinf(tht);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_ps(sint,sint);
                         coefB12_f7414_zmm16r4(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                         sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                         cmul_zmm16r4(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                         ear2  = ear;
                         sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                         y0    = _mm512_mul_ps(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_ps(sin1,cer2);
                         cei2  = _mm512_mul_ps(sin1,cei2);
                         ear3  = ear;
                         cmul_zmm16r4(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                         eai3  = gms::math::negate_zmm16r4(y0);
                         cexp_zmm16r4(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_ps(sin2,cer3);
                         cei3  = _mm512_mul_ps(sin2,cei3);
                         cmul_zmm16r4(B2rs,B2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm16r4(B1rs,B1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_ps(t0r,t1r);
                         t2i = _mm512_sub_ps(t0i,t1i);
                         cmul_zmm16r4(Eir,Eii,t2r,t2i,&x0,&x1);
                         cmul_zmm16r4(x0,x1,cer,cei,*Esr,*Esi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f749_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                           const float * __restrict __ATTR_ALIGN__(64) pEir,
                                           const float * __restrict __ATTR_ALIGN__(64) pEii,
                                           float * __restrict __ATTR_ALIGN__(64)  Esr,
                                           float * __restrict __ATTR_ALIGN__(64)  Esi) {

                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         register __m512 k0r = _mm512_load_ps(&pk0r[0]);
                         register __m512 Eir = _mm512_load_ps(&pEir[0]);
                         register __m512 Eii = _mm512_load_ps(&pEii[0]);
                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         __m512 ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512 B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                         register __m512 ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512 sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_ps();
                         sqr   = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_ps(k0a,k0a);
                         eai   = _mm512_add_ps(k0r,C078539816339744830961566084582);
                         sint  = xsinf(tht);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_ps(sint,sint);
                         coefB12_f7414_zmm16r4(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                         sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                         cmul_zmm16r4(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                         ear2  = ear;
                         sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                         y0    = _mm512_mul_ps(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_ps(sin1,cer2);
                         cei2  = _mm512_mul_ps(sin1,cei2);
                         ear3  = ear;
                         cmul_zmm16r4(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                         eai3  = gms::math::negate_zmm16r4(y0);
                         cexp_zmm16r4(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_ps(sin2,cer3);
                         cei3  = _mm512_mul_ps(sin2,cei3);
                         cmul_zmm16r4(B2rs,B2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm16r4(B1rs,B1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_ps(t0r,t1r);
                         t2i = _mm512_sub_ps(t0i,t1i);
                         cmul_zmm16r4(Eir,Eii,t2r,t2i,&x0,&x1);
                         cmul_zmm16r4(x0,x1,cer,cei,&resr,&resi);
                         _mm512_store_ps(&Esr[0], resr);
                         _mm512_store_ps(&Esi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Esz_f749_zmm16r4_u(const float * __restrict  ptht,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pk0r,
                                           const float * __restrict  pEir,
                                           const float * __restrict  pEii,
                                           float * __restrict   Esr,
                                           float * __restrict   Esi) {

                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 k0r = _mm512_loadu_ps(&pk0r[0]);
                         register __m512 Eir = _mm512_loadu_ps(&pEir[0]);
                         register __m512 Eii = _mm512_loadu_ps(&pEii[0]);
                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         __m512 ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512 B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                         register __m512 ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512 sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_ps();
                         sqr   = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_ps(k0a,k0a);
                         eai   = _mm512_add_ps(k0r,C078539816339744830961566084582);
                         sint  = xsinf(tht);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_ps(sint,sint);
                         coefB12_f7414_zmm16r4(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                         sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                         cmul_zmm16r4(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                         ear2  = ear;
                         sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                         y0    = _mm512_mul_ps(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_ps(sin1,cer2);
                         cei2  = _mm512_mul_ps(sin1,cei2);
                         ear3  = ear;
                         cmul_zmm16r4(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                         eai3  = gms::math::negate_zmm16r4(y0);
                         cexp_zmm16r4(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_ps(sin2,cer3);
                         cei3  = _mm512_mul_ps(sin2,cei3);
                         cmul_zmm16r4(B2rs,B2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm16r4(B1rs,B1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_ps(t0r,t1r);
                         t2i = _mm512_sub_ps(t0i,t1i);
                         cmul_zmm16r4(Eir,Eii,t2r,t2i,&x0,&x1);
                         cmul_zmm16r4(x0,x1,cer,cei,&resr,&resi);
                         _mm512_storeu_ps(&Esr[0], resr);
                         _mm512_storeu_ps(&Esi[0], resi);
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
                   void Hsz_f749_zmm16r4(const __m512 tht,
                                         const __m512 k0a,
                                         const __m512 k0r,
                                         const __m512 Hir,
                                         const __m512 Hii,
                                         __m512 * __restrict Hsr,
                                         __m512 * __restrict Hsi) {

                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512 A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                         register __m512 ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         register __m512 sqr,x0,x1,k0a2,y0;
                         ear   = _mm512_setzero_ps();
                         sqr   = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_ps(k0a,k0a);
                         eai   = _mm512_add_ps(k0r,C078539816339744830961566084582);
                         sint  = xsinf(tht);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_ps(sint,sint);
                         coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                         cmul_zmm16r4(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                         ear2  = ear;
                         sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                         y0    = _mm512_mul_ps(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_ps(sin2,cer2);
                         cei2  = _mm512_mul_ps(sin2,cei2);
                         ear3  = ear;
                         cmul_zmm16r4(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                         eai3  = gms::math::negate_zmm16r4(y0);
                         cexp_zmm16r4(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_ps(sin1,cer3);
                         cei3  = _mm512_mul_ps(sin1,cei3);
                         cmul_zmm16r4(A2rs,A2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm16r4(A1rs,A1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_ps(t0r,t1r);
                         t2i = _mm512_sub_ps(t0i,t1i);
                         cmul_zmm16r4(Hir,Hii,t2r,t2i,&x0,&x1);
                         cmul_zmm16r4(x0,x1,cer,cei,*Hsr,*Hsi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f749_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0r,
                                           const float * __restrict __ATTR_ALIGN__(64) pHir,
                                           const float * __restrict __ATTR_ALIGN__(64) pHii,
                                           float * __restrict __ATTR_ALIGN__(64)  Hsr,
                                           float * __restrict __ATTR_ALIGN__(64)  Hsi) {

                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         register __m512 k0r = _mm512_load_ps(&pk0r[0]);
                         register __m512 Eir = _mm512_load_ps(&pHir[0]);
                         register __m512 Eii = _mm512_load_ps(&pHii[0]);
                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512 A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                         register __m512 ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512 sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_ps();
                         sqr   = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_ps(k0a,k0a);
                         eai   = _mm512_add_ps(k0r,C078539816339744830961566084582);
                         sint  = xsinf(tht);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_ps(sint,sint);
                         coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                         cmul_zmm16r4(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                         ear2  = ear;
                         sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                         y0    = _mm512_mul_ps(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_ps(sin2,cer2);
                         cei2  = _mm512_mul_ps(sin2,cei2);
                         ear3  = ear;
                         cmul_zmm16r4(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                         eai3  = gms::math::negate_zmm16r4(y0);
                         cexp_zmm16r4(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_ps(sin1,cer3);
                         cei3  = _mm512_mul_ps(sin1,cei3);
                         cmul_zmm16r4(A2rs,A2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm16r4(A1rs,A1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_ps(t0r,t1r);
                         t2i = _mm512_sub_ps(t0i,t1i);
                         cmul_zmm16r4(Hir,Hii,t2r,t2i,&x0,&x1);
                         cmul_zmm16r4(x0,x1,cer,cei,&resr,*resi);
                         _mm512_store_ps(&Hsr[0], resr);
                         _mm512_store_ps(&Hsi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   void Hsz_f749_zmm16r4_u(const float * __restrict  ptht,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pk0r,
                                           const float * __restrict  pHir,
                                           const float * __restrict  pHii,
                                           float * __restrict __ATTR_ALIGN__(64)  Hsr,
                                           float * __restrict __ATTR_ALIGN__(64)  Hsi) {

                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 k0r = _mm512_loadu_ps(&pk0r[0]);
                         register __m512 Eir = _mm512_loadu_ps(&pHir[0]);
                         register __m512 Eii = _mm512_loadu_ps(&pHii[0]);
                         const __m512 C078539816339744830961566084582  = _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                         const __m512 C6283185307179586476925286766559 = _mm512_set1_ps(6.283185307179586476925286766559f);
                         const __m512 C10                              = _mm512_set1_ps(1.0f);
                         register __m512 ear,eai,cer,cei,sint,sin2t,sin1,sin2;
                         register __m512 A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                         register __m512 ear2,eai2,ear3,eai3,cer2,cei2,cer3,cei3,t0r,t0i,t1r,t1i,t2r,t2i;
                         __m512 sqr,x0,x1,k0a2,y0,resr,resi;
                         ear   = _mm512_setzero_ps();
                         sqr   = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0r));
                         k0a2  = _mm512_add_ps(k0a,k0a);
                         eai   = _mm512_add_ps(k0r,C078539816339744830961566084582);
                         sint  = xsinf(tht);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         sin2t = _mm512_add_ps(sint,sint);
                         coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                         sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                         cmul_zmm16r4(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                         ear2  = ear;
                         sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                         y0    = _mm512_mul_ps(k0a2,sint);
                         eai2  = y0;
                         cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                         cer2  = _mm512_mul_ps(sin2,cer2);
                         cei2  = _mm512_mul_ps(sin2,cei2);
                         ear3  = ear;
                         cmul_zmm16r4(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                         eai3  = gms::math::negate_zmm16r4(y0);
                         cexp_zmm16r4(ear3,eai3,&cer3,&cei3);
                         cer3  = _mm512_mul_ps(sin1,cer3);
                         cei3  = _mm512_mul_ps(sin1,cei3);
                         cmul_zmm16r4(A2rs,A2is,cer2,cei2,&t0r,&t0i);
                         cmul_zmm16r4(A1rs,A1is,cer3,cei3,&t1r,&t1i);
                         t2r = _mm512_sub_ps(t0r,t1r);
                         t2i = _mm512_sub_ps(t0i,t1i);
                         cmul_zmm16r4(Hir,Hii,t2r,t2i,&x0,&x1);
                         cmul_zmm16r4(x0,x1,cer,cei,&resr,*resi);
                         _mm512_storeu_ps(&Hsr[0], resr);
                         _mm512_storeu_ps(&Hsi[0], resi);
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
                   __m512 rcs_f7411_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 tht) {

                          const __m512 C10   = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512 B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          coefB12_f7414_zmm16r4(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                          sin2t = _mm512_add_ps(sint,sint);
                          ear1  = _mm512_setzero_ps();
                          x0    = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm16r4(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm16r4(x0);
                          cmul_zmm16r4(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                          sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          cer1  = _mm512_mul_ps(sin1,cer1);
                          cei1  = _mm512_mul_ps(sin1,cei1);
                          cmul_zmm16r4(B2rs,B2is,cer1,cei1,&t0r,&t0i);
                          sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_ps(sin2,cer2);
                          cei2  = _mm512_mul_ps(sin2,cei2);
                          cmul_zmm16r4(B1rs,B1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_ps(t0r,t1r);
                          x1    = _mm512_sub_ps(t0i,t1i);
                          cabs  = cabs_zmm16r4(x0,x1);
                          rcs   = _mm512_mul_ps(ink0,cabs);
                          return (rcs);
                }


 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7411_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                            const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0 = _mm512_load_ps(&pk0[0]);
                          register __m512 a  = _mm512_load_ps(&pa[0]);
                          register __m512 tht= _mm512_load_ps(&ptht[0]);
                          const __m512 C10   = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512 B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          __m512 t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          coefB12_f7414_zmm16r4(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                          sin2t = _mm512_add_ps(sint,sint);
                          ear1  = _mm512_setzero_ps();
                          x0    = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm16r4(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm16r4(x0);
                          cmul_zmm16r4(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                          sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          cer1  = _mm512_mul_ps(sin1,cer1);
                          cei1  = _mm512_mul_ps(sin1,cei1);
                          cmul_zmm16r4(B2rs,B2is,cer1,cei1,&t0r,&t0i);
                          sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_ps(sin2,cer2);
                          cei2  = _mm512_mul_ps(sin2,cei2);
                          cmul_zmm16r4(B1rs,B1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_ps(t0r,t1r);
                          x1    = _mm512_sub_ps(t0i,t1i);
                          cabs  = cabs_zmm16r4(x0,x1);
                          rcs   = _mm512_mul_ps(ink0,cabs);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7411_zmm16r4_u(const float * __restrict pk0,
                                            const float * __restrict  pa,
                                            const float * __restrict  ptht) {

                          register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a  = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht= _mm512_loadu_ps(&ptht[0]);
                          const __m512 C10   = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512 B1r,B1i,B2r,B2i,B1rs,B1is,B2rs,B2is;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          __m512 t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          coefB12_f7414_zmm16r4(k0a,tht,&B1r,&B1i,&B2r,&B2i);
                          sin2t = _mm512_add_ps(sint,sint);
                          ear1  = _mm512_setzero_ps();
                          x0    = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm16r4(B1r,B1i,B1r,B1i,&B1rs,&B1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm16r4(x0);
                          cmul_zmm16r4(B2r,B2i,B2r,B2i,&B2rs,&B2is);
                          sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          cer1  = _mm512_mul_ps(sin1,cer1);
                          cei1  = _mm512_mul_ps(sin1,cei1);
                          cmul_zmm16r4(B2rs,B2is,cer1,cei1,&t0r,&t0i);
                          sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_ps(sin2,cer2);
                          cei2  = _mm512_mul_ps(sin2,cei2);
                          cmul_zmm16r4(B1rs,B1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_ps(t0r,t1r);
                          x1    = _mm512_sub_ps(t0i,t1i);
                          cabs  = cabs_zmm16r4(x0,x1);
                          rcs   = _mm512_mul_ps(ink0,cabs);
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
                   __m512 rcs_f7412_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 tht) {

                          const __m512 C10   = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512 A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                          sin2t = _mm512_add_ps(sint,sint);
                          ear1  = _mm512_setzero_ps();
                          x0    = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm16r4(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm16r4(x0);
                          cmul_zmm16r4(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                          sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                          cer1  = _mm512_mul_ps(sin2,cer1);
                          cei1  = _mm512_mul_ps(sin2,cei1);
                          cmul_zmm16r4(A2rs,A2is,cer1,cei1,&t0r,&t0i);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_ps(sin1,cer2);
                          cei2  = _mm512_mul_ps(sin1,cei2);
                          cmul_zmm16r4(A1rs,A1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_ps(t0r,t1r);
                          x1    = _mm512_sub_ps(t0i,t1i);
                          cabs  = cabs_zmm16r4(x0,x1);
                          rcs   = _mm512_mul_ps(ink0,cabs);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7412_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0 = _mm512_load_ps(&pk0[0]);
                          register __m512 a  = _mm512_load_ps(&pa[0]);
                          register __m512 tht= _mm512_load_ps(&ptht[0]);
                          const __m512 C10   = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512 A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                          sin2t = _mm512_add_ps(sint,sint);
                          ear1  = _mm512_setzero_ps();
                          x0    = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm16r4(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm16r4(x0);
                          cmul_zmm16r4(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                          sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                          cer1  = _mm512_mul_ps(sin2,cer1);
                          cei1  = _mm512_mul_ps(sin2,cei1);
                          cmul_zmm16r4(A2rs,A2is,cer1,cei1,&t0r,&t0i);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_ps(sin1,cer2);
                          cei2  = _mm512_mul_ps(sin1,cei2);
                          cmul_zmm16r4(A1rs,A1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_ps(t0r,t1r);
                          x1    = _mm512_sub_ps(t0i,t1i);
                          cabs  = cabs_zmm16r4(x0,x1);
                          rcs   = _mm512_mul_ps(ink0,cabs);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7412_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  ptht) {

                          register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a  = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht= _mm512_loadu_ps(&ptht[0]);
                          const __m512 C10   = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,sint,sin2t,sin1,sin2;
                          register __m512 A1r,A1i,A2r,A2i,A1rs,A1is,A2rs,A2is;
                          register __m512 ear1,eai1,cer1,cei1,ear2,eai2,cer2,cei2;
                          register __m512 t0r,t0i,t1r,t1i,cabs,rcs,x0,x1;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          sint = xsinf(tht);
                          coefA12_f7413_zmm16r4(k0a,tht,&A1r,&A1i,&A2r,&A2i);
                          sin2t = _mm512_add_ps(sint,sint);
                          ear1  = _mm512_setzero_ps();
                          x0    = _mm512_mul_ps(_mm512_add_ps(k0a,k0a),sint);
                          eai1  = x0;
                          cmul_zmm16r4(A1r,A1i,A1r,A1i,&A1rs,&A1is);
                          ear2  = ear1;
                          eai2  = gms::math::negate_zmm16r4(x0);
                          cmul_zmm16r4(A2r,A2i,A2r,A2i,&A2rs,&A2is);
                          sin1  = _mm512_div_ps(_mm512_sub_ps(C10,sint),sin2t);
                          cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                          sin2  = _mm512_div_ps(_mm512_add_ps(C10,sint),sin2t);
                          cer1  = _mm512_mul_ps(sin2,cer1);
                          cei1  = _mm512_mul_ps(sin2,cei1);
                          cmul_zmm16r4(A2rs,A2is,cer1,cei1,&t0r,&t0i);
                          cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                          cer2  = _mm512_mul_ps(sin1,cer2);
                          cei2  = _mm512_mul_ps(sin1,cei2);
                          cmul_zmm16r4(A1rs,A1is,cer2,cei2,&t1r,&t1i);
                          x0    = _mm512_sub_ps(t0r,t1r);
                          x1    = _mm512_sub_ps(t0i,t1i);
                          cabs  = cabs_zmm16r4(x0,x1);
                          rcs   = _mm512_mul_ps(ink0,cabs);
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
                   __m512 rcs_f7416_f7417_zmm16r4(const __m512 k0,
                                                  const __m512 a,
                                                  const __m512 tht) {

                          const __m512 C0001  = _mm512_set1_ps(0.0001f); //very close value to zero incidence (rad)
                          const __m512 C00    = _mm512_setzero_ps();
                          const __m512 C10    = _mm512_set1_ps(1.0f);
                          const __mmask16 m1  = _mm512_cmp_ps_mask(tht,C0001,_CMP_LE_OQ);
                          const __mmask16 m2  = _mm512_cmp_ps_mask(tht,C00,_CMPEQ_OQ);
                          register __m512 k0a,ink0,rcs;
                          
                          if(m1) {
                              register __m512 k0a2,sint,arg,sarg,carg,sints,x0,x1,rat;
                              k0a  = _mm512_mul_ps(k0,a);
                              sint = xsinf(tht);
                              k0a2 = _mm512_add_ps(k0a,k0a);
                              arg  = _mm512_mul_ps(k0a2,sint);
                              ink0 = _mm512_rcp14_ps(k0);
                              sarg = xsinf(arg);
                              sints= _mm512_mul_ps(sint,sint);
                              x0   = _mm512_mul_ps(carg,carg);
                              x1   = _mm512_mul_ps(sarg,sarg);
                              rat  = _mm512_div_ps(x1,sints);
                              rcs  = _mm512_mul_ps(ink0, _mm512_add_ps(rat,carg));
                              return (rcs);
                          }
                          else if(m2) {
                              register __m512 x0,k0a2;
                              k0a  = _mm512_mul_ps(k0,a);
                              ink0 = _mm512_rcp14_ps(k0);
                              k0a2 = _mm512_add_ps(k0a,k0a);
                              x0   = _mm512_fmadd_ps(k0a2,k0a2,C10);
                              rcs  = _mm512_mul_ps(ink0,x0);
                              return (rcs);
                          }
                          
                  }
                  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7416_f7417_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                                    const float * __restrict __ATTR_ALIGN__(64) pa,
                                                    const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          const __m512 C0001  = _mm512_set1_ps(0.0001f); //very close value to zero incidence (rad)
                          const __m512 C00    = _mm512_setzero_ps();
                          const __m512 C10    = _mm512_set1_ps(1.0f);
                          const __mmask16 m1  = _mm512_cmp_ps_mask(tht,C0001,_CMP_LE_OQ);
                          const __mmask16 m2  = _mm512_cmp_ps_mask(tht,C00,_CMPEQ_OQ);
                          register __m512 k0a,ink0,rcs;
                          
                          if(m1) {
                              register __m512 k0a2,sint,arg,sarg,carg,sints,x0,x1,rat;
                              k0a  = _mm512_mul_ps(k0,a);
                              sint = xsinf(tht);
                              k0a2 = _mm512_add_ps(k0a,k0a);
                              arg  = _mm512_mul_ps(k0a2,sint);
                              ink0 = _mm512_rcp14_ps(k0);
                              sarg = xsinf(arg);
                              sints= _mm512_mul_ps(sint,sint);
                              x0   = _mm512_mul_ps(carg,carg);
                              x1   = _mm512_mul_ps(sarg,sarg);
                              rat  = _mm512_div_ps(x1,sints);
                              rcs  = _mm512_mul_ps(ink0, _mm512_add_ps(rat,carg));
                              return (rcs);
                          }
                          else if(m2) {
                              register __m512 x0,k0a2;
                              k0a  = _mm512_mul_ps(k0,a);
                              ink0 = _mm512_rcp14_ps(k0);
                              k0a2 = _mm512_add_ps(k0a,k0a);
                              x0   = _mm512_fmadd_ps(k0a2,k0a2,C10);
                              rcs  = _mm512_mul_ps(ink0,x0);
                              return (rcs);
                          }
                         
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7416_f7417_zmm16r4_u(const float * __restrict  pk0,
                                                    const float * __restrict  pa,
                                                    const float * __restrict  ptht) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          const __m512 C0001  = _mm512_set1_ps(0.0001f); //very close value to zero incidence (rad)
                          const __m512 C00    = _mm512_setzero_ps();
                          const __m512 C10    = _mm512_set1_ps(1.0f);
                          const __mmask16 m1  = _mm512_cmp_ps_mask(tht,C0001,_CMP_LE_OQ);
                          const __mmask16 m2  = _mm512_cmp_ps_mask(tht,C00,_CMPEQ_OQ);
                          register __m512 k0a,ink0,rcs;
                          
                          if(m1) {
                              register __m512 k0a2,sint,arg,sarg,carg,sints,x0,x1,rat;
                              k0a  = _mm512_mul_ps(k0,a);
                              sint = xsinf(tht);
                              k0a2 = _mm512_add_ps(k0a,k0a);
                              arg  = _mm512_mul_ps(k0a2,sint);
                              ink0 = _mm512_rcp14_ps(k0);
                              sarg = xsinf(arg);
                              sints= _mm512_mul_ps(sint,sint);
                              x0   = _mm512_mul_ps(carg,carg);
                              x1   = _mm512_mul_ps(sarg,sarg);
                              rat  = _mm512_div_ps(x1,sints);
                              rcs  = _mm512_mul_ps(ink0, _mm512_add_ps(rat,carg));
                              return (rcs);
                          }
                          else if(m2) {
                              register __m512 x0,k0a2;
                              k0a  = _mm512_mul_ps(k0,a);
                              ink0 = _mm512_rcp14_ps(k0);
                              k0a2 = _mm512_add_ps(k0a,k0a);
                              x0   = _mm512_fmadd_ps(k0a2,k0a2,C10);
                              rcs  = _mm512_mul_ps(ink0,x0);
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
                   __m512 rcs_f7418_zmm16r4(const __m512 k0,
                                            const __m512 a) {
                                              
                          const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 C160                             = _mm512_set1_ps(16.0f);
                          const __m512 C80                              = _mm512_set1_ps(8.0f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,k0a8,sin,k0a16,rat;
                          register __m512 rcs;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          k0a8 = _mm512_mul_ps(C80,k0a);
                          k0a16= _mm512_mul_ps(C160,k0a);
                          sin  = xsinf(k0a8);
                          rat  = _mm512_div_ps(sin,k0a16);
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(C10,rat));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7418_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa) {
                             
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);                 
                          const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 C160                             = _mm512_set1_ps(16.0f);
                          const __m512 C80                              = _mm512_set1_ps(8.0f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,k0a8,sin,k0a16,rat;
                          register __m512 rcs;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          k0a8 = _mm512_mul_ps(C80,k0a);
                          k0a16= _mm512_mul_ps(C160,k0a);
                          sin  = xsinf(k0a8);
                          rat  = _mm512_div_ps(sin,k0a16);
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(C10,rat));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7418_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa) {
                             
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);                 
                          const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 C160                             = _mm512_set1_ps(16.0f);
                          const __m512 C80                              = _mm512_set1_ps(8.0f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          register __m512 ink0,k0a,k0a8,sin,k0a16,rat;
                          register __m512 rcs;
                          k0a  = _mm512_mul_ps(k0,a);
                          ink0 = _mm512_rcp14_ps(k0);
                          k0a8 = _mm512_mul_ps(C80,k0a);
                          k0a16= _mm512_mul_ps(C160,k0a);
                          sin  = xsinf(k0a8);
                          rat  = _mm512_div_ps(sin,k0a16);
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(C10,rat));
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
                   __m512 rcs_f7419_zmm16r4() { return _mm512_setzero_ps();}


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
                   __m512 rcs_f7420_zmm16r4(const __m512 k0,
                                            const __m512 tht) {

                          const __m512 C0001                            = _mm512_set1_ps(0.0001f); 
                          const __m512 C157079632679489661923132169164  = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          const __m512 diff  = _mm512_sub_ps(C157079632679489661923132169164,
                                                             _mm512_abs_ps(tht));
                          const __mmask16 m1 = _mm512_cmp_ps_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                             register __m512 ink0,sint,sin2t,abs,sqr,x0;
                             register __m512 rcs;
                             ink0   = _mm512_rcp14_ps(k0);
                             sint   = xsinf(tht);
                             sint2t = _mm512_mul_ps(sint,sint);
                             x0     = _mm512_add_ps(C10,_mm512_abs_ps(sint));
                             sqr    = _mm512_div_ps(x0,sin2t);
                             rcs    = _mm512_mul_ps(ink0,_mm512_mul_ps(sqr,sqr));
                             return (rcs);
                         }
                          else {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                         }
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7420_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0 = _mm512_load_ps(&pk0[0]);
                          register __m512 tht= _mm512_load_ps(&ptht[0]); 
                          const __m512 C0001                            = _mm512_set1_ps(0.0001f); 
                          const __m512 C157079632679489661923132169164  = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          const __m512 diff  = _mm512_sub_ps(C157079632679489661923132169164,
                                                             _mm512_abs_ps(tht));
                          const __mmask16 m1 = _mm512_cmp_ps_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                             register __m512 ink0,sint,sin2t,abs,sqr,x0;
                             register __m512 rcs;
                             ink0   = _mm512_rcp14_ps(k0);
                             sint   = xsinf(tht);
                             sint2t = _mm512_mul_ps(sint,sint);
                             x0     = _mm512_add_ps(C10,_mm512_abs_ps(sint));
                             sqr    = _mm512_div_ps(x0,sin2t);
                             rcs    = _mm512_mul_ps(ink0,_mm512_mul_ps(sqr,sqr));
                             return (rcs);
                         }
                          else {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                         }
                  }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7420_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  ptht) {

                          register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                          register __m512 tht= _mm512_loadu_ps(&ptht[0]); 
                          const __m512 C0001                            = _mm512_set1_ps(0.0001f); 
                          const __m512 C157079632679489661923132169164  = _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          const __m512 diff  = _mm512_sub_ps(C157079632679489661923132169164,
                                                             _mm512_abs_ps(tht));
                          const __mmask16 m1 = _mm512_cmp_ps_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                             register __m512 ink0,sint,sin2t,abs,sqr,x0;
                             register __m512 rcs;
                             ink0   = _mm512_rcp14_ps(k0);
                             sint   = xsinf(tht);
                             sint2t = _mm512_mul_ps(sint,sint);
                             x0     = _mm512_add_ps(C10,_mm512_abs_ps(sint));
                             sqr    = _mm512_div_ps(x0,sin2t);
                             rcs    = _mm512_mul_ps(ink0,_mm512_mul_ps(sqr,sqr));
                             return (rcs);
                         }
                          else {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
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
                   __m512 rcs_f7421_zmm16r4(const __m512 k0,
                                            const __m512 a,
                                            const __m512 tht) {

                          const __m512 C0001                            = _mm512_set1_ps(0.0001f); 
                          const __m512 C0318309886183790671537767526745 = _mm512_set1_ps(0.318309886183790671537767526745f);
                          const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          const __m512 diff  = _mm512_sub_ps(C157079632679489661923132169164,
                                                             _mm512_abs_ps(tht));
                          const __mmask16 m1 = _mm512_cmp_ps_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                              register __m512 ink0,k0a,k0a2,sint,sin2t,x0,x1;
                              register __m512 rcs,;
                              k0a  = _mm512_mul_ps(k0,a);
                              sint = xsinf(tht);
                              k0a2 = _mm512_mul_ps(k0a,k0a);
                              sin2t= xsinf(sint,sint);
                              ink0 = _mm512_rcp14_ps(k0);
                              x0   = _mm512_mul_ps(_mm512_mul_ps(k0a2,k0a2),
                                                          C0318309886183790671537767526745);
                              x1   = _mm512_div_ps(_mm512_sub_ps(C10,sin2t),sint);
                              x0   = _mm512_mul_ps(x0,x1);
                              rcs  = _mm512_mul_ps(ink0,_mm512_mul_ps(x0,x0));
                              return (rcs);
                       }
                       else {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                      }

                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7421_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht) {

                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          const __m512 C0001                            = _mm512_set1_ps(0.0001f); 
                          const __m512 C0318309886183790671537767526745 = _mm512_set1_ps(0.318309886183790671537767526745f);
                          const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          const __m512 diff  = _mm512_sub_ps(C157079632679489661923132169164,
                                                             _mm512_abs_ps(tht));
                          const __mmask16 m1 = _mm512_cmp_ps_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                              register __m512 ink0,k0a,k0a2,sint,sin2t,x0,x1;
                              register __m512 rcs,;
                              k0a  = _mm512_mul_ps(k0,a);
                              sint = xsinf(tht);
                              k0a2 = _mm512_mul_ps(k0a,k0a);
                              sin2t= xsinf(sint,sint);
                              ink0 = _mm512_rcp14_ps(k0);
                              x0   = _mm512_mul_ps(_mm512_mul_ps(k0a2,k0a2),
                                                          C0318309886183790671537767526745);
                              x1   = _mm512_div_ps(_mm512_sub_ps(C10,sin2t),sint);
                              x0   = _mm512_mul_ps(x0,x1);
                              rcs  = _mm512_mul_ps(ink0,_mm512_mul_ps(x0,x0));
                              return (rcs);
                       }
                       else {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                      }

                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7421_zmm16r4_u(const float * __restrict  pk0,
                                              const float * __restrict  pa,
                                              const float * __restrict  ptht) {

                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          const __m512 C0001                            = _mm512_set1_ps(0.0001f); 
                          const __m512 C0318309886183790671537767526745 = _mm512_set1_ps(0.318309886183790671537767526745f);
                          const __m512 C314159265358979323846264338328  = _mm512_set1_ps(3.14159265358979323846264338328f);
                          const __m512 C10                              = _mm512_set1_ps(1.0f);
                          const __m512 diff  = _mm512_sub_ps(C157079632679489661923132169164,
                                                             _mm512_abs_ps(tht));
                          const __mmask16 m1 = _mm512_cmp_ps_mask(C0001,diff,_CMP_LE_OQ);
                          if(m1) {
                              register __m512 ink0,k0a,k0a2,sint,sin2t,x0,x1;
                              register __m512 rcs,;
                              k0a  = _mm512_mul_ps(k0,a);
                              sint = xsinf(tht);
                              k0a2 = _mm512_mul_ps(k0a,k0a);
                              sin2t= xsinf(sint,sint);
                              ink0 = _mm512_rcp14_ps(k0);
                              x0   = _mm512_mul_ps(_mm512_mul_ps(k0a2,k0a2),
                                                          C0318309886183790671537767526745);
                              x1   = _mm512_div_ps(_mm512_sub_ps(C10,sin2t),sint);
                              x0   = _mm512_mul_ps(x0,x1);
                              rcs  = _mm512_mul_ps(ink0,_mm512_mul_ps(x0,x0));
                              return (rcs);
                       }
                       else {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
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
                   void rcs_f7422_zmm16r4(const __m512 tht1,
                                          const __m512 tht2,
                                          const __m512 k0,
                                          const __m512 a,
                                          __m512 * __restrict rcs1,
                                          __m512 * __restrict rcs2) {

                        const __m512 C05 = _mm512_set1_ps(0.5f);
                        register __m512 A1r1,A1i1,A1r2,A1i2;
                        register __m512 A2r1,A2i1,A2r2,A2i2;
                        register __m512 ear1,eai1,ear2,ea2;
                        register __m512 cer1,cei1,cer2,cei2;
                        register __m512 t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512 ink0,sints,costm,sint1,sint2;
                        register __m512 cost12,rat,k0a,tp,tm;
                        register __m512 x0,x1,x2,x3;
                        k0a   = _mm512_mul_ps(k0,a);
                        sint1 = xsinf(tht1);
                        ear1  = _mm512_setzero_ps();
                        tp    = _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                        coefA12_f7413_zmm16r4(k0a,tht1,&A1r1,&A1i1,&A2r1,&A2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                        sint2 = xsinf(tht2);
                        sints = _mm512_add_ps(sint1,sint2);
                        eai1  = _mm512_mul_ps(k0a,sints);
                        cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_ps(k0);
                        eai2  = gms::math::negate_zmm16r4(eai1);
                        cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                        sintp = xsinf(tp);
                        costm = xcosf(tm);
                        x0    = _mm512_div_ps(_mm512_add_ps(sintp,costm),sints);
                        coefA12_f7413_zmm16r4(k0a,tht2,&A1r2,&A1i2,&A2r2,&A2i2);
                        cer1  = _mm512_mul_ps(cer1,x0);
                        cer2  = _mm512_mul_ps(cer2,x0);
                        cei1  = _mm512_mul_ps(cei1,x0);
                        cei2  = _mm512_mul_ps(cei2,x0);
                        cmul_zmm16r4(A2r1,A2i1,A2r2,A2i2,&t0r,&t0i);
                        cmul_zmm16r4(A1r1,A1i1,A1r2,A1i2,&t1r,&t1i);
                        cmul_zmm16r4(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm16r4(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_ps(res1r,res2r);
                        ear2 = _mm512_add_ps(res1r,res2r);
                        eai1 = _mm512_sub_ps(resi1,res2i);
                        eai2 = _mm512_add_ps(res1i,res2i);
                        *rcs1= _mm512_mul_ps(ink0,cabs_zmm16r4(ear1,eai1));
                        *rcs2= _mm512_mul_ps(ink0,cabs_zmm16r4(ear2,eai2));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                            const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0,
                                            const float * __restrict __ATTR_ALIGN__(64) pa,
                                            float * __restrict __ATTR_ALIGN__(64) rcs1,
                                            float * __restrict __ATTR_ALIGN__(64) rcs2) {

                        register __m512  tht1 = _mm512_load_ps(&ptht1[0]);
                        register __m512  tht2 = _mm512_load_ps(&ptht2[0]);
                        register __m512  k0   = _mm512_load_ps(&pk0[0]);
                        register __m512  a    = _mm512_load_ps(&pa[0]);
                        const __m512 C05 = _mm512_set1_ps(0.5f);
                        register __m512 A1r1,A1i1,A1r2,A1i2;
                        register __m512 A2r1,A2i1,A2r2,A2i2;
                        register __m512 ear1,eai1,ear2,ea2;
                        register __m512 cer1,cei1,cer2,cei2;
                        register __m512 t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512 ink0,sints,costm,sint1,sint2;
                        register __m512 cost12,rat,k0a,tp,tm;
                        __m512 x0,x1,x2,x3;
                        k0a   = _mm512_mul_ps(k0,a);
                        sint1 = xsinf(tht1);
                        ear1  = _mm512_setzero_ps();
                        tp    = _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                        coefA12_f7413_zmm16r4(k0a,tht1,&A1r1,&A1i1,&A2r1,&A2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                        sint2 = xsinf(tht2);
                        sints = _mm512_add_ps(sint1,sint2);
                        eai1  = _mm512_mul_ps(k0a,sints);
                        cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_ps(k0);
                        eai2  = gms::math::negate_zmm16r4(eai1);
                        cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                        sintp = xsinf(tp);
                        costm = xcosf(tm);
                        x0    = _mm512_div_ps(_mm512_add_ps(sintp,costm),sints);
                        coefA12_f7413_zmm16r4(k0a,tht2,&A1r2,&A1i2,&A2r2,&A2i2);
                        cer1  = _mm512_mul_ps(cer1,x0);
                        cer2  = _mm512_mul_ps(cer2,x0);
                        cei1  = _mm512_mul_ps(cei1,x0);
                        cei2  = _mm512_mul_ps(cei2,x0);
                        cmul_zmm16r4(A2r1,A2i1,A2r2,A2i2,&t0r,&t0i);
                        cmul_zmm16r4(A1r1,A1i1,A1r2,A1i2,&t1r,&t1i);
                        cmul_zmm16r4(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm16r4(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_ps(res1r,res2r);
                        ear2 = _mm512_add_ps(res1r,res2r);
                        eai1 = _mm512_sub_ps(resi1,res2i);
                        eai2 = _mm512_add_ps(res1i,res2i);
                        _mm512_store_ps(&rcs1[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear1,eai1)));
                        _mm512_store_ps(&rcs2[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear2,eai2)));
                }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_zmm16r4_u(const float * __restrict  ptht1,
                                            const float * __restrict  ptht2,
                                            const float * __restrict  pk0,
                                            const float * __restrict pa,
                                            float * __restrict  rcs1,
                                            float * __restrict  rcs2) {

                        register __m512  tht1 = _mm512_loadu_ps(&ptht1[0]);
                        register __m512  tht2 = _mm512_loadu_ps(&ptht2[0]);
                        register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512  a    = _mm512_loadu_ps(&pa[0]);
                        const __m512 C05 = _mm512_set1_ps(0.5f);
                        register __m512 A1r1,A1i1,A1r2,A1i2;
                        register __m512 A2r1,A2i1,A2r2,A2i2;
                        register __m512 ear1,eai1,ear2,ea2;
                        register __m512 cer1,cei1,cer2,cei2;
                        register __m512 t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512 ink0,sints,costm,sint1,sint2;
                        register __m512 cost12,rat,k0a,tp,tm;
                        __m512 x0,x1,x2,x3;
                        k0a   = _mm512_mul_ps(k0,a);
                        sint1 = xsinf(tht1);
                        ear1  = _mm512_setzero_ps();
                        tp    = _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                        coefA12_f7413_zmm16r4(k0a,tht1,&A1r1,&A1i1,&A2r1,&A2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                        sint2 = xsinf(tht2);
                        sints = _mm512_add_ps(sint1,sint2);
                        eai1  = _mm512_mul_ps(k0a,sints);
                        cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_ps(k0);
                        eai2  = gms::math::negate_zmm16r4(eai1);
                        cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                        sintp = xsinf(tp);
                        costm = xcosf(tm);
                        x0    = _mm512_div_ps(_mm512_add_ps(sintp,costm),sints);
                        coefA12_f7413_zmm16r4(k0a,tht2,&A1r2,&A1i2,&A2r2,&A2i2);
                        cer1  = _mm512_mul_ps(cer1,x0);
                        cer2  = _mm512_mul_ps(cer2,x0);
                        cei1  = _mm512_mul_ps(cei1,x0);
                        cei2  = _mm512_mul_ps(cei2,x0);
                        cmul_zmm16r4(A2r1,A2i1,A2r2,A2i2,&t0r,&t0i);
                        cmul_zmm16r4(A1r1,A1i1,A1r2,A1i2,&t1r,&t1i);
                        cmul_zmm16r4(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm16r4(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_ps(res1r,res2r);
                        ear2 = _mm512_add_ps(res1r,res2r);
                        eai1 = _mm512_sub_ps(resi1,res2i);
                        eai2 = _mm512_add_ps(res1i,res2i);
                        _mm512_storeu_ps(&rcs1[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear1,eai1)));
                        _mm512_storeu_ps(&rcs2[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear2,eai2)));
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
                   void rcs_f7423_zmm16r4(const __m512 tht1,
                                          const __m512 tht2,
                                          const __m512 k0,
                                          const __m512 a,
                                          __m512 * __restrict rcs1,
                                          __m512 * __restrict rcs2) {

                        const __m512 C05 = _mm512_set1_ps(0.5f);
                        register __m512 B1r1,B1i1,B1r2,B1i2;
                        register __m512 B2r1,B2i1,B2r2,B2i2;
                        register __m512 ear1,eai1,ear2,ea2;
                        register __m512 cer1,cei1,cer2,cei2;
                        register __m512 t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512 ink0,sints,costm,sint1,sint2;
                        register __m512 cost12,rat,k0a,tp,tm;
                        register __m512 x0,x1,x2,x3;
                        k0a   = _mm512_mul_ps(k0,a);
                        sint1 = xsinf(tht1);
                        ear1  = _mm512_setzero_ps();
                        tp    = _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                        coefB12_f7414_zmm16r4(k0a,tht1,&B1r1,&B1i1,&B2r1,&B2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                        sint2 = xsinf(tht2);
                        sints = _mm512_add_ps(sint1,sint2);
                        eai1  = _mm512_mul_ps(k0a,sints);
                        cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_ps(k0);
                        eai2  = gms::math::negate_zmm16r4(eai1);
                        cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                        sintp = xsinf(tp);
                        costm = xcosf(tm);
                        x0    = _mm512_div_ps(_mm512_add_ps(sintp,costm),sints);
                        coefB12_f7414_zmm16r4(k0a,tht2,&B1r2,&B1i2,&B2r2,&B2i2);
                        cer1  = _mm512_mul_ps(cer1,x0);
                        cer2  = _mm512_mul_ps(cer2,x0);
                        cei1  = _mm512_mul_ps(cei1,x0);
                        cei2  = _mm512_mul_ps(cei2,x0);
                        cmul_zmm16r4(B2r1,B2i1,B2r2,B2i2,&t0r,&t0i);
                        cmul_zmm16r4(B1r1,B1i1,B1r2,B1i2,&t1r,&t1i);
                        cmul_zmm16r4(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm16r4(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_ps(res1r,res2r);
                        ear2 = _mm512_add_ps(res1r,res2r);
                        eai1 = _mm512_sub_ps(resi1,res2i);
                        eai2 = _mm512_add_ps(res1i,res2i);
                        *rcs1= _mm512_mul_ps(ink0,cabs_zmm16r4(ear1,eai1));
                        *rcs2= _mm512_mul_ps(ink0,cabs_zmm16r4(ear2,eai2));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7423_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                            const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0,
                                            const float * __restrict __ATTR_ALIGN__(64) pa,
                                            float * __restrict __ATTR_ALIGN__(64) rcs1,
                                            float * __restrict __ATTR_ALIGN__(64) rcs2) {

                        register __m512  tht1 = _mm512_load_ps(&ptht1[0]);
                        register __m512  tht2 = _mm512_load_ps(&ptht2[0]);
                        register __m512  k0   = _mm512_load_ps(&pk0[0]);
                        register __m512  a    = _mm512_load_ps(&pa[0]);

                        const __m512 C05 = _mm512_set1_ps(0.5f);
                        register __m512 B1r1,B1i1,B1r2,B1i2;
                        register __m512 B2r1,B2i1,B2r2,B2i2;
                        register __m512 ear1,eai1,ear2,ea2;
                        register __m512 cer1,cei1,cer2,cei2;
                        register __m512 t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512 ink0,sints,costm,sint1,sint2;
                        register __m512 cost12,rat,k0a,tp,tm;
                        register __m512 x0,x1,x2,x3;
                        k0a   = _mm512_mul_ps(k0,a);
                        sint1 = xsinf(tht1);
                        ear1  = _mm512_setzero_ps();
                        tp    = _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                        coefB12_f7414_zmm16r4(k0a,tht1,&B1r1,&B1i1,&B2r1,&B2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                        sint2 = xsinf(tht2);
                        sints = _mm512_add_ps(sint1,sint2);
                        eai1  = _mm512_mul_ps(k0a,sints);
                        cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_ps(k0);
                        eai2  = gms::math::negate_zmm16r4(eai1);
                        cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                        sintp = xsinf(tp);
                        costm = xcosf(tm);
                        x0    = _mm512_div_ps(_mm512_add_ps(sintp,costm),sints);
                        coefB12_f7414_zmm16r4(k0a,tht2,&B1r2,&B1i2,&B2r2,&B2i2);
                        cer1  = _mm512_mul_ps(cer1,x0);
                        cer2  = _mm512_mul_ps(cer2,x0);
                        cei1  = _mm512_mul_ps(cei1,x0);
                        cei2  = _mm512_mul_ps(cei2,x0);
                        cmul_zmm16r4(B2r1,B2i1,B2r2,B2i2,&t0r,&t0i);
                        cmul_zmm16r4(B1r1,B1i1,B1r2,B1i2,&t1r,&t1i);
                        cmul_zmm16r4(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm16r4(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_ps(res1r,res2r);
                        ear2 = _mm512_add_ps(res1r,res2r);
                        eai1 = _mm512_sub_ps(resi1,res2i);
                        eai2 = _mm512_add_ps(res1i,res2i);
                        _mm512_store_ps(&rcs1[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear1,eai1)));
                        _mm512_store_ps(&rcs2[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear2,eai2)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7423_zmm16r4_u(const float * __restrict  ptht1,
                                            const float * __restrict  ptht2,
                                            const float * __restrict  pk0,
                                            const float * __restrict  pa,
                                            float * __restrict  rcs1,
                                            float * __restrict  rcs2) {

                        register __m512  tht1 = _mm512_loadu_ps(&ptht1[0]);
                        register __m512  tht2 = _mm512_loadu_ps(&ptht2[0]);
                        register __m512  k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512  a    = _mm512_loadu_ps(&pa[0]);

                        const __m512 C05 = _mm512_set1_ps(0.5f);
                        register __m512 B1r1,B1i1,B1r2,B1i2;
                        register __m512 B2r1,B2i1,B2r2,B2i2;
                        register __m512 ear1,eai1,ear2,ea2;
                        register __m512 cer1,cei1,cer2,cei2;
                        register __m512 t0r,t0i,t1r,t1i,resr1,resi1,resr2,resi2;
                        register __m512 ink0,sints,costm,sint1,sint2;
                        register __m512 cost12,rat,k0a,tp,tm;
                        register __m512 x0,x1,x2,x3;
                        k0a   = _mm512_mul_ps(k0,a);
                        sint1 = xsinf(tht1);
                        ear1  = _mm512_setzero_ps();
                        tp    = _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                        coefB12_f7414_zmm16r4(k0a,tht1,&B1r1,&B1i1,&B2r1,&B2i1);
                        ear2  = ear1;
                        tm    = _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                        sint2 = xsinf(tht2);
                        sints = _mm512_add_ps(sint1,sint2);
                        eai1  = _mm512_mul_ps(k0a,sints);
                        cexp_zmm16r4(ear1,eai1,&cer1,&cei1);
                        ink0  = _mm512_rcp14_ps(k0);
                        eai2  = gms::math::negate_zmm16r4(eai1);
                        cexp_zmm16r4(ear2,eai2,&cer2,&cei2);
                        sintp = xsinf(tp);
                        costm = xcosf(tm);
                        x0    = _mm512_div_ps(_mm512_add_ps(sintp,costm),sints);
                        coefB12_f7414_zmm16r4(k0a,tht2,&B1r2,&B1i2,&B2r2,&B2i2);
                        cer1  = _mm512_mul_ps(cer1,x0);
                        cer2  = _mm512_mul_ps(cer2,x0);
                        cei1  = _mm512_mul_ps(cei1,x0);
                        cei2  = _mm512_mul_ps(cei2,x0);
                        cmul_zmm16r4(B2r1,B2i1,B2r2,B2i2,&t0r,&t0i);
                        cmul_zmm16r4(B1r1,B1i1,B1r2,B1i2,&t1r,&t1i);
                        cmul_zmm16r4(t0r,t0i,cer1,cei1,&resr1,&resi1);
                        cmul_zmm16r4(t1r,t1i,cer2,cei2,&resr2,&resi2);
                        ear1 = _mm512_sub_ps(res1r,res2r);
                        ear2 = _mm512_add_ps(res1r,res2r);
                        eai1 = _mm512_sub_ps(resi1,res2i);
                        eai2 = _mm512_add_ps(res1i,res2i);
                        _mm512_storeu_ps(&rcs1[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear1,eai1)));
                        _mm512_storeu_ps(&rcs2[0], _mm512_mul_ps(ink0,cabs_zmm16r4(ear2,eai2)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void rcs_f7422_f7423_zmm16r4(const __m512 tht1,
                                                const __m512 tht2,
                                                const __m512 k0,
                                                const __m512 a,
                                                __m512 * __restrict rcsx1,
                                                __m512 * __restrict rcsx2,
                                                __m512 * __restrict rcsy1,
                                                __m512 * __restrict rcsy2) {

                        rcs_f7422_zmm16r4(tht1,tht2,k0,a,*rcsx1,*rcsx2);
                        rcs_f7423_zmm16r4(tht1,tht2,k0,a,*rcsy1,*rcsy2);
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
                   __m512 rcs_f7424_zmm16r4(const __m512 tht1,
                                            const __m512 tht2,
                                            const __m512 k0,
                                            const __m512 a) {

                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          register __m512 thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512 rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_ps(k0,a);
                          thtp =  _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                          sint1=  xsinf(tht1);
                          ink0 =  _mm512_rcp14_ps(k0);
                          sint2=  xsinf(tht2);
                          arg  = _mm512_mul_ps(k0a,_mm512_add_ps(sint1,sint2));
                          sintp= xsinf(thtp);
                          thtm =  _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                          sarg =  xsinf(arg);
                          sqr1 = _mm512_mul_ps(sarg,sarg);
                          costm=  xcosf(thtm);
                          thtp = _mm512_div_ps(sqr1,_mm512_mul_ps(costm,costm));
                          carg =  xcosf(arg);
                          sqr2 = _mm512_mul_ps(carg,carg);
                          thtm = _mm512_div_ps(sqr2,_mm512_mul_ps(sintp,sintp));
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(thtp,thtm));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7424_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512 tht1  = _mm512_load_ps(&ptht1[0]);
                          register __m512 tht2  = _mm512_load_ps(&ptht2[0]);
                          register __m512 k0    = _mm512_load_ps(&pk0[0]);
                          register __m512 a     = _mm512_load_ps(&pa[0]);
                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          register __m512 thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512 rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_ps(k0,a);
                          thtp =  _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                          sint1=  xsinf(tht1);
                          ink0 =  _mm512_rcp14_ps(k0);
                          sint2=  xsinf(tht2);
                          arg  = _mm512_mul_ps(k0a,_mm512_add_ps(sint1,sint2));
                          sintp= xsinf(thtp);
                          thtm =  _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                          sarg =  xsinf(arg);
                          sqr1 = _mm512_mul_ps(sarg,sarg);
                          costm=  xcosf(thtm);
                          thtp = _mm512_div_ps(sqr1,_mm512_mul_ps(costm,costm));
                          carg =  xcosf(arg);
                          sqr2 = _mm512_mul_ps(carg,carg);
                          thtm = _mm512_div_ps(sqr2,_mm512_mul_ps(sintp,sintp));
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(thtp,thtm));
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7424_zmm16r4_u(const float * __restrict  ptht1,
                                              const float * __restrict  ptht2,
                                              const float * __restrict  pk0,
                                              const float * __restrict  pa) {

                          register __m512 tht1  = _mm512_loadu_ps(&ptht1[0]);
                          register __m512 tht2  = _mm512_loadu_ps(&ptht2[0]);
                          register __m512 k0    = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a     = _mm512_loadu_ps(&pa[0]);
                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          register __m512 thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512 rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_ps(k0,a);
                          thtp =  _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                          sint1=  xsinf(tht1);
                          ink0 =  _mm512_rcp14_ps(k0);
                          sint2=  xsinf(tht2);
                          arg  = _mm512_mul_ps(k0a,_mm512_add_ps(sint1,sint2));
                          sintp= xsinf(thtp);
                          thtm =  _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                          sarg =  xsinf(arg);
                          sqr1 = _mm512_mul_ps(sarg,sarg);
                          costm=  xcosf(thtm);
                          thtp = _mm512_div_ps(sqr1,_mm512_mul_ps(costm,costm));
                          carg =  xcosf(arg);
                          sqr2 = _mm512_mul_ps(carg,carg);
                          thtm = _mm512_div_ps(sqr2,_mm512_mul_ps(sintp,sintp));
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(thtp,thtm));
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
                   __m512 rcs_f7427_zmm16r4(const __m512 tht1,
                                            const __m512 tht2,
                                            const __m512 k0,
                                            const __m512 a) {

                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          register __m512 thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512 rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_ps(k0,a);
                          thtp =  _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                          sint1=  xsinf(tht1);
                          ink0 =  _mm512_rcp14_ps(k0);
                          sint2=  xsinf(tht2);
                          arg  = _mm512_mul_ps(k0a,_mm512_add_ps(sint1,sint2));
                          sintp= xsinf(thtp);
                          thtm =  _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                          sarg =  xsinf(arg);
                          sqr1 = _mm512_mul_ps(sarg,sarg);
                          costm=  xcosf(thtm);
                          thtp = _mm512_div_ps(sqr1,_mm512_mul_ps(sintp,sintp));
                          carg =  xcosf(arg);
                          sqr2 = _mm512_mul_ps(carg,carg);
                          thtm = _mm512_div_ps(sqr2,_mm512_mul_ps(costm,costm));
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(thtp,thtm));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7427_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) ptht1,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht2,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512 tht1  = _mm512_load_ps(&ptht1[0]);
                          register __m512 tht2  = _mm512_load_ps(&ptht2[0]);
                          register __m512 k0    = _mm512_load_ps(&pk0[0]);
                          register __m512 a     = _mm512_load_ps(&pa[0]);

                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          register __m512 thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512 rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_ps(k0,a);
                          thtp =  _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                          sint1=  xsinf(tht1);
                          ink0 =  _mm512_rcp14_ps(k0);
                          sint2=  xsinf(tht2);
                          arg  = _mm512_mul_ps(k0a,_mm512_add_ps(sint1,sint2));
                          sintp= xsinf(thtp);
                          thtm =  _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                          sarg =  xsinf(arg);
                          sqr1 = _mm512_mul_ps(sarg,sarg);
                          costm=  xcosf(thtm);
                          thtp = _mm512_div_ps(sqr1,_mm512_mul_ps(sintp,sintp));
                          carg =  xcosf(arg);
                          sqr2 = _mm512_mul_ps(carg,carg);
                          thtm = _mm512_div_ps(sqr2,_mm512_mul_ps(costm,costm));
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(thtp,thtm));
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f7427_zmm16r4_u(  const float * __restrict  ptht1,
                                                const float * __restrict  ptht2,
                                                const float * __restrict  pk0,
                                                const float * __restrict  pa) {

                          register __m512 tht1  = _mm512_loadu_ps(&ptht1[0]);
                          register __m512 tht2  = _mm512_loadu_ps(&ptht2[0]);
                          register __m512 k0    = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a     = _mm512_loadu_ps(&pa[0]);

                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          register __m512 thtp,thtm,sint1,sint2,k0a,ink0;
                          register __m512 rcs,sintp,costm,arg,sarg,carg,sqr1,sqr2;
                          k0a  =  _mm512_mul_ps(k0,a);
                          thtp =  _mm512_mul_ps(_mm512_add_ps(tht1,tht2),C05);
                          sint1=  xsinf(tht1);
                          ink0 =  _mm512_rcp14_ps(k0);
                          sint2=  xsinf(tht2);
                          arg  = _mm512_mul_ps(k0a,_mm512_add_ps(sint1,sint2));
                          sintp= xsinf(thtp);
                          thtm =  _mm512_mul_ps(_mm512_sub_ps(tht1,tht2),C05);
                          sarg =  xsinf(arg);
                          sqr1 = _mm512_mul_ps(sarg,sarg);
                          costm=  xcosf(thtm);
                          thtp = _mm512_div_ps(sqr1,_mm512_mul_ps(sintp,sintp));
                          carg =  xcosf(arg);
                          sqr2 = _mm512_mul_ps(carg,carg);
                          thtm = _mm512_div_ps(sqr2,_mm512_mul_ps(costm,costm));
                          rcs  = _mm512_mul_ps(ink0,_mm512_add_ps(thtp,thtm));
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
                   void Eth_f751_zmm16r4(const __m512 tht,
                                         const __m512 phi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 a,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_ps(k0,r);
                        sinp = xsinf(phi);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(sinp,cost));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        *Esr = _mm512_mul_ps(x1,cer);
                        cei  = _mm512_mul_ps(x0,cei);
                        *Esi = _mm512_mul_ps(x1,cei);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f751_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                         const float * __restrict __ATTR_ALIGN__(64) pphi,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         const float * __restrict __ATTR_ALIGN__(64) pr,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         float * __restrict __ATTR_ALIGN__(64) Esr,
                                         float * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 a   = _mm512_load_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_ps(k0,r);
                        sinp = xsinf(phi);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(sinp,cost));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_store_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_store_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));
                 }



                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f751_zmm16r4_u(const float * __restrict ptht,
                                         const float * __restrict  pphi,
                                         const float * __restrict  pk0,
                                         const float * __restrict  pr,
                                         const float * __restrict  pa,
                                         float * __restrict  Esr,
                                         float * __restrict  Esi) {

                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 a   = _mm512_loadu_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_ps(k0,r);
                        sinp = xsinf(phi);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(sinp,cost));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_storeu_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_storeu_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));
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
                   void Eph_f751_zmm16r4(const __m512 tht,
                                         const __m512 tht2
                                         const __m512 phi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 a,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        const __m512 C0212206590789193781025178351163 = 
                                              _mm512_set1_ps(0.212206590789193781025178351163f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cosp,sint,x0,x1,inr,sint2;
                        k0r  = _mm512_mul_ps(k0,r);
                        sint = xsinf(tht);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cosp = xcosf(phi);
                        k02  = _mm512_mul_ps(k0,k0);
                        sint2= xsinf(tht2);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_fmadd_ps(_mm512_mul_ps(C0212206590789193781025178351163,
                                         sint),sint2,_mm512_mul_ps(C0424413181578387562050356702327,cosp));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        *Esr = _mm512_mul_ps(x1,cer);
                        cei  = _mm512_mul_ps(x0,cei);
                        *Esi = _mm512_mul_ps(x1,cei);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f751_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                         const float * __restrict __ATTR_ALIGN__(64) ptht2
                                         const float * __restrict __ATTR_ALIGN__(64) pphi,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         const float * __restrict __ATTR_ALIGN__(64) pr,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         float * __restrict __ATTR_ALIGN__(64) Esr,
                                         float * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 tht2= _mm512_load_ps(&ptht2[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 a   = _mm512_load_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        const __m512 C0212206590789193781025178351163 = 
                                              _mm512_set1_ps(0.212206590789193781025178351163f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cosp,sint,x0,x1,inr,sint2;
                        k0r  = _mm512_mul_ps(k0,r);
                        sint = xsinf(tht);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cosp = xcosf(phi);
                        k02  = _mm512_mul_ps(k0,k0);
                        sint2= xsinf(tht2);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_fmadd_ps(_mm512_mul_ps(C0212206590789193781025178351163,
                                         sint),sint2,_mm512_mul_ps(C0424413181578387562050356702327,cosp));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_store_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_store_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f751_zmm16r4_u(const float * __restrict  ptht,
                                         const float * __restrict  ptht2
                                         const float * __restrict  pphi,
                                         const float * __restrict  pk0,
                                         const float * __restrict  pr,
                                         const float * __restrict  pa,
                                         float * __restrict  Esr,
                                         float * __restrict  Esi) {

                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 tht2= _mm512_loadu_ps(&ptht2[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 a   = _mm512_loadu_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        const __m512 C0212206590789193781025178351163 = 
                                              _mm512_set1_ps(0.212206590789193781025178351163f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cosp,sint,x0,x1,inr,sint2;
                        k0r  = _mm512_mul_ps(k0,r);
                        sint = xsinf(tht);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cosp = xcosf(phi);
                        k02  = _mm512_mul_ps(k0,k0);
                        sint2= xsinf(tht2);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_fmadd_ps(_mm512_mul_ps(C0212206590789193781025178351163,
                                         sint),sint2,_mm512_mul_ps(C0424413181578387562050356702327,cosp));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_storeu_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_storeu_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));
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
                   void Eth_f752_zmm16r4(const __m512 tht,
                                         const __m512 tht2,
                                         const __m512 phi,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 a,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,cosp,cost2,x0,x1,inr,t0;
                        k0r  = _mm512_mul_ps(k0,r);
                        cosp = xcosf(phi);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        cost2= xcosf(tht2);
                        t0   = _mm512_mul_ps(cost,cost2);
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(cosp,t0));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        *Esr = _mm512_mul_ps(x1,cer);
                        cei  = _mm512_mul_ps(x0,cei);
                        *Esi = _mm512_mul_ps(x1,cei);
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f752_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                         const float * __restrict __ATTR_ALIGN__(64) ptht2
                                         const float * __restrict __ATTR_ALIGN__(64) pphi,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         const float * __restrict __ATTR_ALIGN__(64) pr,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         float * __restrict __ATTR_ALIGN__(64) Esr,
                                         float * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 tht2= _mm512_load_ps(&ptht2[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 a   = _mm512_load_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,cosp,cost2,x0,x1,inr,t0;
                        k0r  = _mm512_mul_ps(k0,r);
                        cosp = xcosf(phi);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        cost2= xcosf(tht2);
                        t0   = _mm512_mul_ps(cost,cost2);
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(cosp,t0));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_store_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_store_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eth_f752_zmm16r4_u(const float * __restrict ptht,
                                         const float * __restrict ptht2
                                         const float * __restrict  pphi,
                                         const float * __restrict  pk0,
                                         const float * __restrict  pr,
                                         const float * __restrict  pa,
                                         float * __restrict Esr,
                                         float * __restrict  Esi) {

                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 tht2= _mm512_loadu_ps(&ptht2[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 a   = _mm512_loadu_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,cosp,cost2,x0,x1,inr,t0;
                        k0r  = _mm512_mul_ps(k0,r);
                        cosp = xcosf(phi);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        cost2= xcosf(tht2);
                        t0   = _mm512_mul_ps(cost,cost2);
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(cosp,t0));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_storeu_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_storeu_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));
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
                   void Eph_f752_zmm16r4(const __m512 tht,
                                         const __m512 phi2,
                                         const __m512 k0,
                                         const __m512 r,
                                         const __m512 a,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {
                        
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_ps(k0,r);
                        sinp = xsinf(phi2);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(sinp,cost));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        *Esr = _mm512_mul_ps(x1,cer);
                        cei  = _mm512_mul_ps(x0,cei);
                        *Esi = _mm512_mul_ps(x1,cei);       
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f752_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           float * __restrict __ATTR_ALIGN__(64) Esr,
                                           float * __restrict __ATTR_ALIGN__(64) Esi) {
                        
                        register __m512 tht  = _mm512_load_ps(&ptht[0]);
                        register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        register __m512 r    = _mm512_load_ps(&pr[0]);
                        register __m512 a    = _mm512_load_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_ps(k0,r);
                        sinp = xsinf(phi2);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(sinp,cost));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_store_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_store_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));       
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   void Eph_f752_zmm16r4_u(const float * __restrict  ptht,
                                           const float * __restrict  pphi2,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pr,
                                           const float * __restrict pa,
                                           float * __restrict Esr,
                                           float * __restrict Esi) {
                        
                        
                        register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                        register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const __m512 C0424413181578387562050356702327 = 
                                              _mm512_set1_ps(0.424413181578387562050356702327f);
                        register __m512 ear,eai,cer,cei,k0r,k02,a3;
                        register __m512 cost,sinp,x0,x1,inr;
                        k0r  = _mm512_mul_ps(k0,r);
                        sinp = xsinf(phi2);
                        ear  = _mm512_setzero_ps();
                        eai  = k0r;
                        cost = xcosf(tht);
                        k02  = _mm512_mul_ps(k0,k0);
                        inr  = _mm512_rcp14_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        a3   = _mm512_mul_ps(a,_mm512_mul_ps(a,a));
                        k02  = gms::math::negate_zmm16r4(k02);
                        x0   = _mm512_mul_ps(C0424413181578387562050356702327,
                                                         _mm512_mul_ps(sinp,cost));
                        x1   = _mm512_mul_ps(k02,a3);
                        cer  = _mm512_mul_ps(x0,cer);
                        _mm512_storeu_ps(&Esr[0] ,_mm512_mul_ps(x1,cer));
                        cei  = _mm512_mul_ps(x0,cei);
                        _mm512_storeu_ps(&Esi[0] ,_mm512_mul_ps(x1,cei));       
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
                   __m512 rcs_f753_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht) {
                  
                          const __m512 C0565884242104516749400475603102 = 
                                        _mm512_set1_ps(0.565884242104516749400475603102f);
                          const __m512 C20 = _mm512_set1_ps(2.0f);
                          register __m512 k02,k04,a2,a6,sint,sqr,x0,x1;
                          register __m512 rcs;
                          k02   = _mm512_mul_ps(k0,k0);
                          a2    = _mm512_mul_ps(a,a);
                          sint  = xsinf(tht);
                          k04   = _mm512_mul_ps(k02,k02);
                          x0    = _mm512_fmadd_ps(sint,sint,C20);
                          a6    = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          sqr   = _mm512_mul_ps(x0,x0);
                          x1    = _mm512_mul_ps(k04,a6);
                          rcs   = _mm512_mul_ps(C0565884242104516749400475603102,
                                                               _mm512_mul_ps(x1,sqr));
                          return (rcs);                  
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f753_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht) {
                  
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                          const __m512 C0565884242104516749400475603102 = 
                                        _mm512_set1_ps(0.565884242104516749400475603102f);
                          const __m512 C20 = _mm512_set1_ps(2.0f);
                          register __m512 k02,k04,a2,a6,sint,sqr,x0,x1;
                          register __m512 rcs;
                          k02   = _mm512_mul_ps(k0,k0);
                          a2    = _mm512_mul_ps(a,a);
                          sint  = xsinf(tht);
                          k04   = _mm512_mul_ps(k02,k02);
                          x0    = _mm512_fmadd_ps(sint,sint,C20);
                          a6    = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          sqr   = _mm512_mul_ps(x0,x0);
                          x1    = _mm512_mul_ps(k04,a6);
                          rcs   = _mm512_mul_ps(C0565884242104516749400475603102,
                                                               _mm512_mul_ps(x1,sqr));
                          return (rcs);                  
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f753_zmm16r4_u(const float * __restrict pk0,
                                             const float * __restrict pa,
                                             const float * __restrict  ptht) {
                  
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                          const __m512 C0565884242104516749400475603102 = 
                                        _mm512_set1_ps(0.565884242104516749400475603102f);
                          const __m512 C20 = _mm512_set1_ps(2.0f);
                          register __m512 k02,k04,a2,a6,sint,sqr,x0,x1;
                          register __m512 rcs;
                          k02   = _mm512_mul_ps(k0,k0);
                          a2    = _mm512_mul_ps(a,a);
                          sint  = xsinf(tht);
                          k04   = _mm512_mul_ps(k02,k02);
                          x0    = _mm512_fmadd_ps(sint,sint,C20);
                          a6    = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          sqr   = _mm512_mul_ps(x0,x0);
                          x1    = _mm512_mul_ps(k04,a6);
                          rcs   = _mm512_mul_ps(C0565884242104516749400475603102,
                                                               _mm512_mul_ps(x1,sqr));
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
                   __m512 rcs_f754_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht) {
                                           
                          const __m512 C2263536968418066997601902412409 = 
                                          _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 k02,k04,a2,a6,cost,cost2t,cost4t,x0;
                          register __m512 rcs;
                          k02  = _mm512_mul_ps(k0,k0);
                          cost = xcosf(tht);
                          a2   = _mm512_mul_ps(a,a);
                          cos2t= _mm512_mul_ps(cost,cost);
                          k04  = _mm512_mul_ps(k02,k02);
                          a6   = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));  
                          x0   = _mm512_mul_ps(k04,a6);
                          cos4t= _mm512_mul_ps(cos2t,cos2t);
                          rcs  = _mm512_mul_ps(C2263536968418066997601902412409,
                                                                 _mm512_mul_ps(x0,cos4t));
                          return (rcs);                   
               }
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f754_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht) {
                  
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 tht = _mm512_load_ps(&ptht[0]);
                                           
                          const __m512 C2263536968418066997601902412409 = 
                                          _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 k02,k04,a2,a6,cost,cost2t,cost4t,x0;
                          register __m512 rcs;
                          k02  = _mm512_mul_ps(k0,k0);
                          cost = xcosf(tht);
                          a2   = _mm512_mul_ps(a,a);
                          cos2t= _mm512_mul_ps(cost,cost);
                          k04  = _mm512_mul_ps(k02,k02);
                          a6   = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));  
                          x0   = _mm512_mul_ps(k04,a6);
                          cos4t= _mm512_mul_ps(cos2t,cos2t);
                          rcs  = _mm512_mul_ps(C2263536968418066997601902412409,
                                                                 _mm512_mul_ps(x0,cos4t));
                          return (rcs);                   
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f754_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ptht) {
                  
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                                           
                          const __m512 C2263536968418066997601902412409 = 
                                          _mm512_set1_ps(2.263536968418066997601902412409f);
                          register __m512 k02,k04,a2,a6,cost,cost2t,cost4t,x0;
                          register __m512 rcs;
                          k02  = _mm512_mul_ps(k0,k0);
                          cost = xcosf(tht);
                          a2   = _mm512_mul_ps(a,a);
                          cos2t= _mm512_mul_ps(cost,cost);
                          k04  = _mm512_mul_ps(k02,k02);
                          a6   = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));  
                          x0   = _mm512_mul_ps(k04,a6);
                          cos4t= _mm512_mul_ps(cos2t,cos2t);
                          rcs  = _mm512_mul_ps(C2263536968418066997601902412409,
                                                                 _mm512_mul_ps(x0,cos4t));
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
                   __m512 rcs_f755_zmm16r4(const __m512 k0,
                                           const __m512 a) {
                                           
                         const __m512 C2263536968418066997601902412409 = 
                                          _mm512_set1_ps(2.263536968418066997601902412409f);
                         register __m512 rcs,k02,k04,a2,a6;
                         k02 = _mm512_mul_ps(k0,k0);
                         a2  = _mm512_mul_ps(a,a);
                         k04 = _mm512_mul_ps(k02,k02);
                         a6  = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                         rcs = _mm512_mul_ps( C2263536968418066997601902412409,
                                                                _mm512_mul_ps(k04,a6));
                         return (rcs);                         
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f755_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa) {
                                         
                         register __m512 k0  = _mm512_load_ps(&pk0[0]);
                         register __m512 a   = _mm512_load_ps(&pa[0]);
                               
                         const __m512 C2263536968418066997601902412409 = 
                                          _mm512_set1_ps(2.263536968418066997601902412409f);
                         register __m512 rcs,k02,k04,a2,a6;
                         k02 = _mm512_mul_ps(k0,k0);
                         a2  = _mm512_mul_ps(a,a);
                         k04 = _mm512_mul_ps(k02,k02);
                         a6  = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                         rcs = _mm512_mul_ps( C2263536968418066997601902412409,
                                                                _mm512_mul_ps(k04,a6));
                         return (rcs);                         
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f755_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict pa) {
                                         
                         register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                         register __m512 a   = _mm512_loadu_ps(&pa[0]);
                               
                         const __m512 C2263536968418066997601902412409 = 
                                          _mm512_set1_ps(2.263536968418066997601902412409f);
                         register __m512 rcs,k02,k04,a2,a6;
                         k02 = _mm512_mul_ps(k0,k0);
                         a2  = _mm512_mul_ps(a,a);
                         k04 = _mm512_mul_ps(k02,k02);
                         a6  = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                         rcs = _mm512_mul_ps( C2263536968418066997601902412409,
                                                                _mm512_mul_ps(k04,a6));
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
                   __m512 rcs_f756_zmm16r4(const __m512 k0,
                                           const __m512 a) {
                                           
                          const __m512 C509295817894065074460428042792 = 
                                                      _mm512_set1_ps(5.09295817894065074460428042792f);
                          register __m512 rcs,k02,k04,a2,a6;
                          k02 = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          k04 = _mm512_mul_ps(k02,k02);
                          a6  = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          rcs = _mm512_mul_ps(C509295817894065074460428042792,
                                                               _mm512_mul_ps(k04,a6));
                          return (rcs);                        
                  }
                  
                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f756_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa) {
                                         
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          const __m512 C509295817894065074460428042792 = 
                                                      _mm512_set1_ps(5.09295817894065074460428042792f);
                          register __m512 rcs,k02,k04,a2,a6;
                          k02 = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          k04 = _mm512_mul_ps(k02,k02);
                          a6  = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          rcs = _mm512_mul_ps(C509295817894065074460428042792,
                                                               _mm512_mul_ps(k04,a6));
                          return (rcs);                        
                  }
                  

                    
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f756_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa) {
                                         
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const __m512 C509295817894065074460428042792 = 
                                                      _mm512_set1_ps(5.09295817894065074460428042792f);
                          register __m512 rcs,k02,k04,a2,a6;
                          k02 = _mm512_mul_ps(k0,k0);
                          a2  = _mm512_mul_ps(a,a);
                          k04 = _mm512_mul_ps(k02,k02);
                          a6  = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          rcs = _mm512_mul_ps(C509295817894065074460428042792,
                                                               _mm512_mul_ps(k04,a6));
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
                   __m512 rcs_f757_zmm16r4() { return (_mm512_setzero_ps());}
                   
                   
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
                   __m512 rcs_f758_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht2) {
                                           
                          return (rcs_f753_zmm16r4(k0,a,tht2));                         
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f758_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                             const float * __restrict __ATTR_ALIGN__(64)  pa,
                                             const float * __restrict __ATTR_ALIGN__(64)  ptht2) {
                                           
                          return (rcs_f753_zmm16r4_a(pk0,pa,ptht2));                         
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f758_zmm16r4_u(const float * __restrict  pk0,
                                             const float * __restrict  pa,
                                             const float * __restrict  ptht2) {
                                           
                          return (rcs_f753_zmm16r4_u(pk0,pa,ptht2));                         
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
                   __m512 rcs_f759_zmm16r4(const __m512 k0,
                                           const __m512 a,
                                           const __m512 tht2) {
                                           
                          return (rcs_f754_zmm16r4(k0,a,tht2));                         
                 }
                 
                 
                 
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f759_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                             const float * __restrict __ATTR_ALIGN__(64)  pa,
                                             const float * __restrict __ATTR_ALIGN__(64)  ptht2) {
                                           
                          return (rcs_f754_zmm16r4_a(pk0,pa,ptht2));                         
                 }
                 
                 
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline   
                   __m512 rcs_f759_zmm16r4_u(const float * __restrict   pk0,
                                             const float * __restrict   pa,
                                             const float * __restrict   ptht2) {
                                           
                          return (rcs_f754_zmm16r4_u(pk0,pa,ptht2));                         
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
                   __m512 rcs_f7513_zmm16r4(const __m512 k0a,
                                            const __m512 tht) {
                              
                                        
                         /* const __m512 C157079632679489661923132169164 = 
                                                _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 C00                             = _mm512_setzero_ps();
                          const __mmask16 m1 = _mm512_cmp_ps_mask(tht,C00,_CMP_EQ_OQ);
                          const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(tht),
                                                                      C157079632679489661923132169164,_CMP_EQ_OQ);
                          if(m1 || m2) {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                          }   */ 
                         const __m512 C39478417604357434475337963999505 = 
                                                              _mm512_set1_ps(39.478417604357434475337963999505f);
                         register __m512 k0a2,sint,arg,carg,sarg,sin2t,x0,c2,s2;
                         register __m512 rcs;
                         k0a2 = _mm512_add_ps(k0a,k0a);
                         sint = xsinf(tht);
                         arg  = _mm512_mul_ps(k0a2,sint);
                         x0   = _mm512_div_ps(k0a,
                                              _mm512_mul_ps(C39478417604357434475337963999505,sint));
                         sin2t= _mm512_mul_ps(sint,sint);
                         carg = xcosf(arg);
                         sarg = xsinf(arg);
                         c2   = _mm512_mul_ps(carg,carg);
                         s2   = _mm512_div_ps(_mm512_mul_ps(sarg,sarg),sint2t);
                         rcs  = _mm512_mul_ps(x0,_mm512_add_ps(c2,s2));
                         return (rcs);
                 }
                 
                 
                  
                 
                 
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7513_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht) {
                              
                                        
                         /* const __m512 C157079632679489661923132169164 = 
                                                _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 C00                             = _mm512_setzero_ps();
                          const __mmask16 m1 = _mm512_cmp_ps_mask(tht,C00,_CMP_EQ_OQ);
                          const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(tht),
                                                                      C157079632679489661923132169164,_CMP_EQ_OQ);
                          if(m1 || m2) {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                          }   */ 
                         register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                         register __m512 tht = _mm512_load_ps(&ptht[0]);
                         const __m512 C39478417604357434475337963999505 = 
                                                              _mm512_set1_ps(39.478417604357434475337963999505f);
                         register __m512 k0a2,sint,arg,carg,sarg,sin2t,x0,c2,s2;
                         register __m512 rcs;
                         k0a2 = _mm512_add_ps(k0a,k0a);
                         sint = xsinf(tht);
                         arg  = _mm512_mul_ps(k0a2,sint);
                         x0   = _mm512_div_ps(k0a,
                                              _mm512_mul_ps(C39478417604357434475337963999505,sint));
                         sin2t= _mm512_mul_ps(sint,sint);
                         carg = xcosf(arg);
                         sarg = xsinf(arg);
                         c2   = _mm512_mul_ps(carg,carg);
                         s2   = _mm512_div_ps(_mm512_mul_ps(sarg,sarg),sint2t);
                         rcs  = _mm512_mul_ps(x0,_mm512_add_ps(c2,s2));
                         return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7513_zmm16r4_u(const float * __restrict  pk0a,
                                              const float * __restrict  ptht) {
                              
                                        
                         /* const __m512 C157079632679489661923132169164 = 
                                                _mm512_set1_ps(1.57079632679489661923132169164f);
                          const __m512 C00                             = _mm512_setzero_ps();
                          const __mmask16 m1 = _mm512_cmp_ps_mask(tht,C00,_CMP_EQ_OQ);
                          const __mmask16 m2 = _mm512_cmp_ps_mask(_mm512_abs_ps(tht),
                                                                      C157079632679489661923132169164,_CMP_EQ_OQ);
                          if(m1 || m2) {
                             __m512 NAN = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
                             return (NAN);
                          }   */ 
                         register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                         const __m512 C39478417604357434475337963999505 = 
                                                              _mm512_set1_ps(39.478417604357434475337963999505f);
                         register __m512 k0a2,sint,arg,carg,sarg,sin2t,x0,c2,s2;
                         register __m512 rcs;
                         k0a2 = _mm512_add_ps(k0a,k0a);
                         sint = xsinf(tht);
                         arg  = _mm512_mul_ps(k0a2,sint);
                         x0   = _mm512_div_ps(k0a,
                                              _mm512_mul_ps(C39478417604357434475337963999505,sint));
                         sin2t= _mm512_mul_ps(sint,sint);
                         carg = xcosf(arg);
                         sarg = xsinf(arg);
                         c2   = _mm512_mul_ps(carg,carg);
                         s2   = _mm512_div_ps(_mm512_mul_ps(sarg,sarg),sint2t);
                         rcs  = _mm512_mul_ps(x0,_mm512_add_ps(c2,s2));
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
                   __m512 rcs_f7531_zmm16r4(const __m512 a, // side length of plate
                                            const __m512 gam0) {
                                            
                          const __m512 C1140 = _mm512_set1_ps(114.0f);
                          register __m512 a2,a6,gam2,gam4;
                          register __m512 rcs;
                          a2   = _mm512_mul_ps(a,a);
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          a6   = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          gam4 = _mm512_mul_ps(gam2,gam2);
                          rcs  = _mm512_mul_ps(C1140,_mm512_mul_ps(a6,gam4));
                          return (rcs);                  
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7531_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa, // side length of plate
                                              const float * __restrict __ATTR_ALIGN__(64)  pgam0) {
                                     
                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 gam0 = _mm512_load_ps(&pgam0[0]);      
                          const __m512 C1140 = _mm512_set1_ps(114.0f);
                          register __m512 a2,a6,gam2,gam4;
                          register __m512 rcs;
                          a2   = _mm512_mul_ps(a,a);
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          a6   = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          gam4 = _mm512_mul_ps(gam2,gam2);
                          rcs  = _mm512_mul_ps(C1140,_mm512_mul_ps(a6,gam4));
                          return (rcs);                  
                 }
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7531_zmm16r4_u(const float * __restrict  pa, // side length of plate
                                              const float * __restrict  pgam0) {
                                     
                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);      
                          const __m512 C1140 = _mm512_set1_ps(114.0f);
                          register __m512 a2,a6,gam2,gam4;
                          register __m512 rcs;
                          a2   = _mm512_mul_ps(a,a);
                          gam2 = _mm512_mul_ps(gam0,gam0);
                          a6   = _mm512_mul_ps(a2,_mm512_mul_ps(a2,a2));
                          gam4 = _mm512_mul_ps(gam2,gam2);
                          rcs  = _mm512_mul_ps(C1140,_mm512_mul_ps(a6,gam4));
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
                   __m512 rcs_f7532_zmm16r4(const __m512 A,
                                            const __m512 gam0) {
                                            
                          const __m512 C12566370614359172953850573533118  = 
                                                          _mm512_set1_ps(12.566370614359172953850573533118f);
                          register __m512 rcs,A2,gam02,rat;
                          A2    = _mm512_mul_ps(A,A);
                          gam02 = _mm512_mul_ps(gam0,
                          rat   = _mm512_div_ps(A2,gam02);
                          rcs   = _mm512_mul_ps(C12566370614359172953850573533118,rat);
                          return (rcs);              
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7532_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pA,
                                            const  float * __restrict __ATTR_ALIGN__(64)  gam0) {
                                      
                          register __m512 A    = _mm512_load_ps(&pA[0]);
                          register __m512 gam0 = _mm512_load_ps(&pgam0[0]);     
                          const __m512 C12566370614359172953850573533118  = 
                                                          _mm512_set1_ps(12.566370614359172953850573533118f);
                          register __m512 rcs,A2,gam02,rat;
                          A2    = _mm512_mul_ps(A,A);
                          gam02 = _mm512_mul_ps(gam0,
                          rat   = _mm512_div_ps(A2,gam02);
                          rcs   = _mm512_mul_ps(C12566370614359172953850573533118,rat);
                          return (rcs);              
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7532_zmm16r4_u(const float * __restrict  pA,
                                            const  float * __restrict  gam0) {
                                      
                          register __m512 A    = _mm512_loadu_ps(&pA[0]);
                          register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);     
                          const __m512 C12566370614359172953850573533118  = 
                                                          _mm512_set1_ps(12.566370614359172953850573533118f);
                          register __m512 rcs,A2,gam02,rat;
                          A2    = _mm512_mul_ps(A,A);
                          gam02 = _mm512_mul_ps(gam0,
                          rat   = _mm512_div_ps(A2,gam02);
                          rcs   = _mm512_mul_ps(C12566370614359172953850573533118,rat);
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
                   __m512 rcs_f7533_zmm16r4(const __m512 A,
                                            const __m512 k0a,
                                            const __m512 k0b,
                                            const __m512 tht,
                                            const __m512 phi,
                                            const __m512 gam0) {
                                            
                          const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                          register __m512 A2,gam02,sint,cost,cosp,sinp,fac;
                          register __m512 rcs,trm1,trm2,arg1,arg2,sarg1,sarg2,x0,x1;
                          A2   = _mm512_mul_ps(A,A);
                          cost = xcosf(tht);
                          gam02= _mm512_mul_ps(gam0,gam0);
                          sint = xsinf(tht);
                          fac  = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),
                                                                                             gam02);  
                          cosp = xcosf(phi);
                          arg1 = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                          sarg1= xsinf(arg1);
                          trm1 = _mm512_div_ps(sarg1,arg1);
                          x0   = _mm512_mul_ps(trm1,trm1);
                          sinp = xsinf(phi);
                          fac  = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                          arg2 = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                          sarg2= xsinf(arg2);
                          trm2 = _mm512_div_ps(sarg2,arg2);
                          x1   = _mm512_mul_ps(trm2,trm2);
                          rcs  = _mm512_mul_ps(fac,_mm512_mul_ps(x0,x1));
                          return (rcs);                 
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7533_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pA,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0b,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) pgam0) {
                             
                          register __m512 A    = _mm512_load_ps(&pA[0]);
                          register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          register __m512 k0b  = _mm512_load_ps(&pk0b[0]);
                          register __m512 tht  = _mm512_load_ps(&ptht[0]);
                          register __m512 phi  = _mm512_load_ps(&pphi[0]);
                          register __m512 gam0 = _mm512_load_ps(&pgam0[0]);             
                          const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                          register __m512 A2,gam02,sint,cost,cosp,sinp,fac;
                          register __m512 rcs,trm1,trm2,arg1,arg2,sarg1,sarg2,x0,x1;
                          A2   = _mm512_mul_ps(A,A);
                          cost = xcosf(tht);
                          gam02= _mm512_mul_ps(gam0,gam0);
                          sint = xsinf(tht);
                          fac  = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),
                                                                                             gam02);  
                          cosp = xcosf(phi);
                          arg1 = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                          sarg1= xsinf(arg1);
                          trm1 = _mm512_div_ps(sarg1,arg1);
                          x0   = _mm512_mul_ps(trm1,trm1);
                          sinp = xsinf(phi);
                          fac  = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                          arg2 = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                          sarg2= xsinf(arg2);
                          trm2 = _mm512_div_ps(sarg2,arg2);
                          x1   = _mm512_mul_ps(trm2,trm2);
                          rcs  = _mm512_mul_ps(fac,_mm512_mul_ps(x0,x1));
                          return (rcs);                 
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7533_zmm16r4_u(const float * __restrict  pA,
                                              const float * __restrict pk0a,
                                              const float * __restrict pk0b,
                                              const float * __restrict  ptht,
                                              const float * __restrict  pphi,
                                              const float * __restrict  pgam0) {
                             
                          register __m512 A    = _mm512_loadu_ps(&pA[0]);
                          register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 k0b  = _mm512_loadu_ps(&pk0b[0]);
                          register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                          register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                          register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);             
                          const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                          register __m512 A2,gam02,sint,cost,cosp,sinp,fac;
                          register __m512 rcs,trm1,trm2,arg1,arg2,sarg1,sarg2,x0,x1;
                          A2   = _mm512_mul_ps(A,A);
                          cost = xcosf(tht);
                          gam02= _mm512_mul_ps(gam0,gam0);
                          sint = xsinf(tht);
                          fac  = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),
                                                                                             gam02);  
                          cosp = xcosf(phi);
                          arg1 = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                          sarg1= xsinf(arg1);
                          trm1 = _mm512_div_ps(sarg1,arg1);
                          x0   = _mm512_mul_ps(trm1,trm1);
                          sinp = xsinf(phi);
                          fac  = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                          arg2 = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                          sarg2= xsinf(arg2);
                          trm2 = _mm512_div_ps(sarg2,arg2);
                          x1   = _mm512_mul_ps(trm2,trm2);
                          rcs  = _mm512_mul_ps(fac,_mm512_mul_ps(x0,x1));
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
                   __m512 rcs_f7559_zmm16r4(const __m512 a,
                                            const __m512 b,
                                            const __m512 gam0,
                                            const __m512 k0a,
                                            const __m512 tht,
                                            const __m512 k0b,
                                            const __m512 phi) {
                                            
                         const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                         const __m512 C05                               = _mm512_set1_ps(0.5f);
                         const __m512 C025                              = _mm512_set1_ps(0.25f);
                         register __m512 A2,alp,beta,sint,beta2,sinp,cosp,sin2a,x0,x1,cost;
                         register __m512 rcs,sinb,sina,sinb2,alp2,fac,gam02,aab,den,trm1,trm2,trm3;
                         x0    = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                         sint  = xsinf(tht);
                         A2    = _mm512_mul_ps(x0,x0);
                         cosp  = xcosf(phi);
                         gam02 = _mm512_mul_ps(gam0,gam0);
                         alp   = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                         fac   = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                         sinp  = xsinf(phi);
                         beta  = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                         cost  = xcosf(tht);
                         sinb  = xsinf(beta);
                         beta2 = _mm512_mul_ps(beta,C05);
                         fac   = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                         sina  = xsinf(alp);
                         aab   = _mm512_div_ps(_mm512_add_ps(a,a),b);
                         sinb2 = xsinf(beta2);
                         trm1  = _mm512_fmsub_ps(sina,sina,
                                             _mm512_mul_ps(sinb2,sinb2));
                         x1    = xsinf(_mm512_add_ps(alp,alp));
                         den   = _mm512_fmsub_ps(alp,alp,
                                             _mm512_mul_ps(beta2,beta2));
                         trm3  = _mm512_fmsub_ps(cosp,sinb,
                                             _mm512_mul_ps(sinp,x1));
                         trm2  = _mm512_mul_ps(C025,_mm512_mul_ps(sinp,sinp));
                         trm3  = _mm512_mul_ps(aab,trm3);
                         x0    = _mm512_mul_ps(trm1,trm1);
                         x1    = _mm512_mul_ps(trm3,trm3);
                         alp   = _mm512_fmadd_ps(trm2,x1,x0);
                         beta  = _mm512_div_ps(alp,den);
                         rcs   = _mm512_mul_ps(fac,beta);
                         return (rcs);
                  }
                  
                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7559_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0b,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {
                              
                         register __m512 a    = _mm512_load_ps(&pa[0]);
                         register __m512 b    = _mm512_load_ps(&pb[0]);
                         register __m512 gam0 = _mm512_load_ps(&pgam0[0]);
                         register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                         register __m512 tht  = _mm512_load_ps(&ptht[0]);
                         register __m512 k0b  = _mm512_load_ps(&pk0b[0]);
                         register __m512 phi  = _mm512_load_ps(&pphi[0]);             
                         const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                         const __m512 C05                               = _mm512_set1_ps(0.5f);
                         const __m512 C025                              = _mm512_set1_ps(0.25f);
                         register __m512 A2,alp,beta,sint,beta2,sinp,cosp,sin2a,x0,x1,cost;
                         register __m512 rcs,sinb,sina,sinb2,alp2,fac,gam02,aab,den,trm1,trm2,trm3;
                         x0    = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                         sint  = xsinf(tht);
                         A2    = _mm512_mul_ps(x0,x0);
                         cosp  = xcosf(phi);
                         gam02 = _mm512_mul_ps(gam0,gam0);
                         alp   = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                         fac   = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                         sinp  = xsinf(phi);
                         beta  = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                         cost  = xcosf(tht);
                         sinb  = xsinf(beta);
                         beta2 = _mm512_mul_ps(beta,C05);
                         fac   = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                         sina  = xsinf(alp);
                         aab   = _mm512_div_ps(_mm512_add_ps(a,a),b);
                         sinb2 = xsinf(beta2);
                         trm1  = _mm512_fmsub_ps(sina,sina,
                                             _mm512_mul_ps(sinb2,sinb2));
                         x1    = xsinf(_mm512_add_ps(alp,alp));
                         den   = _mm512_fmsub_ps(alp,alp,
                                             _mm512_mul_ps(beta2,beta2));
                         trm3  = _mm512_fmsub_ps(cosp,sinb,
                                             _mm512_mul_ps(sinp,x1));
                         trm2  = _mm512_mul_ps(C025,_mm512_mul_ps(sinp,sinp));
                         trm3  = _mm512_mul_ps(aab,trm3);
                         x0    = _mm512_mul_ps(trm1,trm1);
                         x1    = _mm512_mul_ps(trm3,trm3);
                         alp   = _mm512_fmadd_ps(trm2,x1,x0);
                         beta  = _mm512_div_ps(alp,den);
                         rcs   = _mm512_mul_ps(fac,beta);
                         return (rcs);
                  }
                  
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7559_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pgam0,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  ptht,
                                              const float * __restrict  pk0b,
                                              const float * __restrict  pphi) {
                              
                         register __m512 a    = _mm512_loadu_ps(&pa[0]);
                         register __m512 b    = _mm512_loadu_ps(&pb[0]);
                         register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);
                         register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                         register __m512 k0b  = _mm512_loadu_ps(&pk0b[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]);             
                         const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                         const __m512 C05                               = _mm512_set1_ps(0.5f);
                         const __m512 C025                              = _mm512_set1_ps(0.25f);
                         register __m512 A2,alp,beta,sint,beta2,sinp,cosp,sin2a,x0,x1,cost;
                         register __m512 rcs,sinb,sina,sinb2,alp2,fac,gam02,aab,den,trm1,trm2,trm3;
                         x0    = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                         sint  = xsinf(tht);
                         A2    = _mm512_mul_ps(x0,x0);
                         cosp  = xcosf(phi);
                         gam02 = _mm512_mul_ps(gam0,gam0);
                         alp   = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                         fac   = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                         sinp  = xsinf(phi);
                         beta  = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                         cost  = xcosf(tht);
                         sinb  = xsinf(beta);
                         beta2 = _mm512_mul_ps(beta,C05);
                         fac   = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                         sina  = xsinf(alp);
                         aab   = _mm512_div_ps(_mm512_add_ps(a,a),b);
                         sinb2 = xsinf(beta2);
                         trm1  = _mm512_fmsub_ps(sina,sina,
                                             _mm512_mul_ps(sinb2,sinb2));
                         x1    = xsinf(_mm512_add_ps(alp,alp));
                         den   = _mm512_fmsub_ps(alp,alp,
                                             _mm512_mul_ps(beta2,beta2));
                         trm3  = _mm512_fmsub_ps(cosp,sinb,
                                             _mm512_mul_ps(sinp,x1));
                         trm2  = _mm512_mul_ps(C025,_mm512_mul_ps(sinp,sinp));
                         trm3  = _mm512_mul_ps(aab,trm3);
                         x0    = _mm512_mul_ps(trm1,trm1);
                         x1    = _mm512_mul_ps(trm3,trm3);
                         alp   = _mm512_fmadd_ps(trm2,x1,x0);
                         beta  = _mm512_div_ps(alp,den);
                         rcs   = _mm512_mul_ps(fac,beta);
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
                   __m512 rcs_f7560_zmm16r4(const __m512 a,
                                            const __m512 b,
                                            const __m512 gam0,
                                            const __m512 tht,
                                            const __m512 phi,
                                            const __m512 k0a) {
                                            
                         const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                         const __m512 C05                               = _mm512_set1_ps(0.5f);
                         const __m512 C40                               = _mm512_set1_ps(4.0f);
                         register __m512 A2,gam02,cost,alp,fac,sina,sina4,alp4;
                         register __m512 rcs,sin2a,alp2,falp4,x0,x1;
                         x0   = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                         cost = xcosf(tht);
                         gam02= _mm512_mul_ps(gam0,gam0);
                         A2   = _mm512_mul_ps(x0,x0);
                         sint = xsinf(tht);
                         cosp = xcosf(phi);
                         fac  = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                         alp  = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                         x1   = _mm512_mul_ps(alp,alp);
                         sina = xsinf(alp);
                         alp4 = _mm512_mul_ps(x1,x1);
                         falp4= _mm512_mul_ps(C40,alp4);
                         alp2 = _mm512_add_ps(alp,alp);
                         sin2a= xsinf(alp2);
                         x0   = _mm512_sub_ps(sin2a,alp2);
                         fac  = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                         x1   = _mm512_mul_ps(x0,x0);
                         A2   = _mm512_mul_ps(sina,sina);
                         gam02= _mm512_div_ps(x1,falp4);
                         x0   = _mm512_mul_ps(A2,A2);
                         x1   = _mm512_div_ps(x0,alp4);
                         rcs  = _mm512_mul_ps(fac,_mm512_add_ps(x1,gam02));
                         return (rcs);                        
                }
                
                
                 __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7560_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {
                              
                         register __m512 a    = _mm512_load_ps(&pa[0]);
                         register __m512 b    = _mm512_load_ps(&pb[0]);
                         register __m512 gam0 = _mm512_load_ps(&pgam0[0]);
                         register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                         register __m512 tht  = _mm512_load_ps(&ptht[0]);
                         register __m512 phi  = _mm512_load_ps(&pphi[0]); 
                                            
                         const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                         const __m512 C05                               = _mm512_set1_ps(0.5f);
                         const __m512 C40                               = _mm512_set1_ps(4.0f);
                         register __m512 A2,gam02,cost,alp,fac,sina,sina4,alp4;
                         register __m512 rcs,sin2a,alp2,falp4,x0,x1;
                         x0   = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                         cost = xcosf(tht);
                         gam02= _mm512_mul_ps(gam0,gam0);
                         A2   = _mm512_mul_ps(x0,x0);
                         sint = xsinf(tht);
                         cosp = xcosf(phi);
                         fac  = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                         alp  = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                         x1   = _mm512_mul_ps(alp,alp);
                         sina = xsinf(alp);
                         alp4 = _mm512_mul_ps(x1,x1);
                         falp4= _mm512_mul_ps(C40,alp4);
                         alp2 = _mm512_add_ps(alp,alp);
                         sin2a= xsinf(alp2);
                         x0   = _mm512_sub_ps(sin2a,alp2);
                         fac  = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                         x1   = _mm512_mul_ps(x0,x0);
                         A2   = _mm512_mul_ps(sina,sina);
                         gam02= _mm512_div_ps(x1,falp4);
                         x0   = _mm512_mul_ps(A2,A2);
                         x1   = _mm512_div_ps(x0,alp4);
                         rcs  = _mm512_mul_ps(fac,_mm512_add_ps(x1,gam02));
                         return (rcs);                        
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7560_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pgam0,
                                              const float * __restrict  pk0a,
                                              const float * __restrict  ptht,
                                              const float * __restrict  pphi) {
                              
                         register __m512 a    = _mm512_loadu_ps(&pa[0]);
                         register __m512 b    = _mm512_loadu_ps(&pb[0]);
                         register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);
                         register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                         register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                         register __m512 phi  = _mm512_loadu_ps(&pphi[0]); 
                                            
                         const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                         const __m512 C05                               = _mm512_set1_ps(0.5f);
                         const __m512 C40                               = _mm512_set1_ps(4.0f);
                         register __m512 A2,gam02,cost,alp,fac,sina,sina4,alp4;
                         register __m512 rcs,sin2a,alp2,falp4,x0,x1;
                         x0   = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                         cost = xcosf(tht);
                         gam02= _mm512_mul_ps(gam0,gam0);
                         A2   = _mm512_mul_ps(x0,x0);
                         sint = xsinf(tht);
                         cosp = xcosf(phi);
                         fac  = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                         alp  = _mm512_mul_ps(k0a,_mm512_mul_ps(sint,cosp));
                         x1   = _mm512_mul_ps(alp,alp);
                         sina = xsinf(alp);
                         alp4 = _mm512_mul_ps(x1,x1);
                         falp4= _mm512_mul_ps(C40,alp4);
                         alp2 = _mm512_add_ps(alp,alp);
                         sin2a= xsinf(alp2);
                         x0   = _mm512_sub_ps(sin2a,alp2);
                         fac  = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                         x1   = _mm512_mul_ps(x0,x0);
                         A2   = _mm512_mul_ps(sina,sina);
                         gam02= _mm512_div_ps(x1,falp4);
                         x0   = _mm512_mul_ps(A2,A2);
                         x1   = _mm512_div_ps(x0,alp4);
                         rcs  = _mm512_mul_ps(fac,_mm512_add_ps(x1,gam02));
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
                   __m512 rcs_f7561_zmm16r4(const __m512 a,
                                            const __m512 b,
                                            const __m512 k0b,
                                            const __m512 tht,
                                            const __m512 phi,
                                            const __m512 gam0) {
                                            
                        const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                        
                        register __m512 A2,gam02,fac,cost,x0,x1,sint,sinp,beta2,sbeta,sbeta4,beta4;
                        register __m512 rcs;
                        x0    = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                        cost  = xcosf(tht);
                        gam02 = _mm512_mul_ps(gam0,gam0);
                        sint  = xsinf(tht);
                        A2    = _mm512_mul_ps(x0,x0);
                        sinp  = xsinf(phi);
                        fac   = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                        x0    = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                        fac   = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                        beta2 = _mm512_mul_ps(x0,C05);
                        sbeta = xsinf(beta2);
                        x0    = _mm512_mul_ps(beta2,beta2);
                        x1    = _mm512_mul_ps(sbeta,sbeta);
                        beta4 = _mm512_mul_ps(x0,x0);   
                        sbeta4= _mm512_mul_ps(x1,x1);
                        sint  = _mm512_div_ps(sbeta4,beta4);
                        rcs   = _mm512_mul_ps(fac,sint);
                        return (rcs);       
               }
               
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7561_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pgam0,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0b,
                                              const float * __restrict __ATTR_ALIGN__(64) ptht,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi) {
                              
                        register __m512 a    = _mm512_load_ps(&pa[0]);
                        register __m512 b    = _mm512_load_ps(&pb[0]);
                        register __m512 gam0 = _mm512_load_ps(&pgam0[0]);
                        register __m512 k0a  = _mm512_load_ps(&pk0b[0]);
                        register __m512 tht  = _mm512_load_ps(&ptht[0]);
                        register __m512 phi  = _mm512_load_ps(&pphi[0]);
                        const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                        
                        register __m512 A2,gam02,fac,cost,x0,x1,sint,sinp,beta2,sbeta,sbeta4,beta4;
                        register __m512 rcs;
                        x0    = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                        cost  = xcosf(tht);
                        gam02 = _mm512_mul_ps(gam0,gam0);
                        sint  = xsinf(tht);
                        A2    = _mm512_mul_ps(x0,x0);
                        sinp  = xsinf(phi);
                        fac   = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                        x0    = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                        fac   = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                        beta2 = _mm512_mul_ps(x0,C05);
                        sbeta = xsinf(beta2);
                        x0    = _mm512_mul_ps(beta2,beta2);
                        x1    = _mm512_mul_ps(sbeta,sbeta);
                        beta4 = _mm512_mul_ps(x0,x0);   
                        sbeta4= _mm512_mul_ps(x1,x1);
                        sint  = _mm512_div_ps(sbeta4,beta4);
                        rcs   = _mm512_mul_ps(fac,sint);
                        return (rcs);       
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline 
                   __m512 rcs_f7561_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pb,
                                              const float * __restrict  pgam0,
                                              const float * __restrict  pk0b,
                                              const float * __restrict  ptht,
                                              const float * __restrict  pphi) {
                              
                        register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        register __m512 b    = _mm512_loadu_ps(&pb[0]);
                        register __m512 gam0 = _mm512_loadu_ps(&pgam0[0]);
                        register __m512 k0a  = _mm512_loadu_ps(&pk0b[0]);
                        register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                        register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                        const __m512 C12566370614359172953850573533118 = 
                                                       _mm512_set1_ps(12.566370614359172953850573533118f);
                        
                        register __m512 A2,gam02,fac,cost,x0,x1,sint,sinp,beta2,sbeta,sbeta4,beta4;
                        register __m512 rcs;
                        x0    = _mm512_mul_ps(_mm512_mul_ps(a,b),C05);
                        cost  = xcosf(tht);
                        gam02 = _mm512_mul_ps(gam0,gam0);
                        sint  = xsinf(tht);
                        A2    = _mm512_mul_ps(x0,x0);
                        sinp  = xsinf(phi);
                        fac   = _mm512_div_ps(_mm512_mul_ps(C12566370614359172953850573533118,A2),gam02);
                        x0    = _mm512_mul_ps(k0b,_mm512_mul_ps(sint,sinp));
                        fac   = _mm512_mul_ps(fac,_mm512_mul_ps(cost,cost));
                        beta2 = _mm512_mul_ps(x0,C05);
                        sbeta = xsinf(beta2);
                        x0    = _mm512_mul_ps(beta2,beta2);
                        x1    = _mm512_mul_ps(sbeta,sbeta);
                        beta4 = _mm512_mul_ps(x0,x0);   
                        sbeta4= _mm512_mul_ps(x1,x1);
                        sint  = _mm512_div_ps(sbeta4,beta4);
                        rcs   = _mm512_mul_ps(fac,sint);
                        return (rcs);       
               }
               
                                            
                  
                

      } // radiolocation


} // gms


























#endif /*__GMS_RCS_PLANAR_SURF_ZMM16R4_HPP__*/
