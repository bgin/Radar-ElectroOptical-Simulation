
#ifndef __GMS_RCS_CYLINDER_ZMM8R8_HPP__
#define __GMS_RCS_CYLINDER_ZMM8R8_HPP__ 200120231636


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

    const unsigned int GMS_RCS_CYLINDER_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_RCS_CYLINDER_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM8R8_FULLVER =
      1000U*GMS_RCS_CYLINDER_ZMM8R8_MAJOR+
      100U*GMS_RCS_CYLINDER_ZMM8R8_MINOR+
      10U*GMS_RCS_CYLINDER_ZMM8R8_MICRO;
    const char * const GMS_RCS_CYLINDER_ZMM8R8_CREATION_DATE = "28-02-2023 13:04 PM +00200 (TUE 28 FEB 2023 GMT+2)";
    const char * const GMS_RCS_CYLINDER_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_CYLINDER_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_CYLINDER_ZMM8R8_DESCRIPTION   = "AVX512 optimized Cylinder Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimddp.hpp"


namespace gms {


          namespace radiolocation {




               namespace {
                   const  static __m512d Ir  = _mm512_setzero_pd();
                   const  static __m512d Ii  = _mm512_set1_pd(1.0);
                   const  static __m512d nIr = _mm512_set1_pd(-0.0);
                   const  static __m512d nIi = _mm512_set1_pd(-1.0);
                   const  static __m512d PI  = _mm512_set1_pd(3.14159265358979323846264338328);

               }


                   /* 
                         Low frequency scattering widths (k0a << 1).
                         Backscatter scattering width for E-field 
                         cylinder-parallel,formula 4.1-19
                    */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f419_zmm8r8(const __m512d a,
                                           const __m512d k0a) {

                          const register __m512d num = _mm512_mul_pd(a, 
                                                           _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d pi4 = _mm512_set1_pd(2.467401100272339654708622749969);
                          const register __m512d c0  = _mm512_set1_pd(0.8905);
                          const register __m512d arg = _mm512_mul_pd(k0a,c0);
                          __m512d ln,ln2,rcs,den;
                          ln = logk(arg);
                          ln2= _mm512_mul_pd(ln,ln);
                          den= _mm512_fmadd_pd(k0a,ln2,pi4);
                          rcs= _mm512_div_pd(num,den);
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f419_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512d a   = _mm512_load_pd(&pa[0]);
                          const register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const register __m512d num = _mm512_mul_pd(a, 
                                                           _mm512_set1_pd(9.869604401089358618834490999876f));
                          const register __m512d pi4 = _mm512_set1_pd(2.467401100272339654708622749969f);
                          const register __m512d c0  = _mm512_set1_pd(0.8905f);
                          const register __m512d arg = _mm512_mul_pd(k0a,c0);
                          __m512d ln,ln2,rcs,den;
                          ln = logk(arg);
                          ln2= _mm512_mul_pd(ln,ln);
                          den= _mm512_fmadd_pd(k0a,ln2,pi4);
                          rcs= _mm512_div_pd(num,den);
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f419_zmm8r8_u(const double * __restrict  pa,
                                             const double * __restrict  pk0a) {

                          const register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          const register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const register __m512d num = _mm512_mul_pd(a, 
                                                           _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d pi4 = _mm512_set1_pd(2.467401100272339654708622749969);
                          const register __m512d c0  = _mm512_set1_pd(0.8905f);
                          const register __m512d arg = _mm512_mul_pd(k0a,c0);
                          __m512d ln,ln2,rcs,den;
                          ln = logk(arg);
                          ln2= _mm512_mul_pd(ln,ln);
                          den= _mm512_fmadd_pd(k0a,ln2,pi4);
                          rcs= _mm512_div_pd(num,den);
                          return (rcs);
              }


                /* 
                         Low frequency scattering widths (k0a << 1).
                         Backscatter scattering width for H-field 
                         cylinder-parallel,formula 4.1-20
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4120_zmm8r8(const __m512d a,
                                            const __m512d k0a) {

                          const register __m512d pi2a = _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d c0   = _mm512_set1_pd(2.25);
                          const register __m512d k0a3 = _mm512_mul_pd(k0a,
                                                                 _mm512_mul_pd(k0a,k0a));
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(pi2a,_mm512_mul_pd(c0,k0a3));
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4120_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512d a    = _mm512_load_pd(&pa[0]);
                          const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const register __m512d pi2a = _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d c0   = _mm512_set1_pd(2.25);
                          const register __m512d k0a3 = _mm512_mul_pd(k0a,
                                                                 _mm512_mul_pd(k0a,k0a));
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(pi2a,_mm512_mul_pd(c0,k0a3));
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4120_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict  pk0a) {

                          const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const register __m512d pi2a = _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d c0   = _mm512_set1_pd(2.25f);
                          const register __m512d k0a3 = _mm512_mul_pd(k0a,
                                                                 _mm512_mul_pd(k0a,k0a));
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(pi2a,_mm512_mul_pd(c0,k0a3));
                          return (rcs);
              }


                /*
                        Bistatic scattering widths, E-field cylinder axis-parallel
                        Formula 4.1-21
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4121_zmm8r8(const __m512d a,
                                            const __m512d k0a) {

                          return (rcs_f4120_zmm8r8(a,k0a));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4121_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          return (rcs_f4120_zmm8r8_a(pa,pk0a));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4121_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a) {

                          return (rcs_f4120_zmm8r8_u(pa,pk0a));
              }



                 /*
                        Bistatic scattering widths, H-field cylinder axis-parallel
                        Formula 4.1-22
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4122_zmm8r8(const __m512d phi,
                                            const __m512d a,
                                            const __m512d k0a) {

                          const register __m512d pi2a= _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_pd(hlf,cosph);
                          sqr   = _mm512_mul_pd(frac,frac);
                          rcs   = _mm512_mul_pd(pi2a,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4122_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphi,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512d phi = _mm512_load_pd(&pphi[0]);
                          const register __m512d a   = _mm512_load_pd(&pa[0]);
                          const register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const register __m512d pi2a= _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_pd(hlf,cosph);
                          sqr   = _mm512_mul_pd(frac,frac);
                          rcs   = _mm512_mul_pd(pi2a,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4122_zmm8r8_u(const double * __restrict  pphi,
                                              const double * __restrict  pa,
                                              const double * __restrict  pk0a) {

                          const register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                          const register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          const register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const register __m512d pi2a= _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_pd(hlf,cosph);
                          sqr   = _mm512_mul_pd(frac,frac);
                          rcs   = _mm512_mul_pd(pi2a,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                 }


                   /*
                       Forward scattering widths, E-field.
                       Formula 4.1-23
                   */
 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4123_zmm8r8(const __m512d a,
                                            const __m512d k0a) {

                          return (rcs_f4120_zmm8r8(a,k0a));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4123_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          return (rcs_f4120_zmm8r8_a(pa,pk0a));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4123_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a) {

                          return (rcs_f4120_zmm8r8_u(pa,pk0a));
              }


                  /*
                       Forward scattering widths, H-field.
                       Formula 4.1-24
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4124_zmm8r8(const __m512d a,
                                            const __m512d k0a) {

                          const register __m512d pi2a= _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const register __m512d qtr = _mm512_set1_pd(0.25);
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(pi2a,_mm512_mul_pd(k0a3,qtr));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4124_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512d a   = _mm512_load_pd(&pa[0]);
                          const register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const register __m512d pi2a= _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const register __m512d qtr = _mm512_set1_pd(0.25);
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(pi2a,_mm512_mul_pd(k0a3,qtr));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4124_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a) {

                          const register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          const register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const register __m512d pi2a= _mm512_mul_pd(a, 
                                                    _mm512_set1_pd(9.869604401089358618834490999876));
                          const register __m512d k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const register __m512d qtr = _mm512_set1_pd(0.25);
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(pi2a,_mm512_mul_pd(k0a3,qtr));
                          return (rcs);
                  }


                    /*
                          Surface currents (k0a << 1), for long cylinder (wire).
                          E-field cylinder axis parallel.
                          Formula 4.1-25
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm8r8(const __m512d eps0,
                                         const __m512d mu0,
                                         const __m512d Er,
                                         const __m512d Ei,
                                         const __m512d k0a,
                                         __m512d * __restrict Kzr,
                                         __m512d * __restrict Kzi) {

                        const register __m512d lna = _mm512_mul_pd(k0a,
                                                     _mm512_set1_pd(0.8905));
                        const register __m512d ln  = logkf(lna);
                        const __m512d sqr = _mm512_sqrt_pd(_mm512_div_pd(eps0,mu0));
                        const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                        __m512d t0r,t0i,t1r,t1i,divr,divi;
                        t0r = nIr;
                        t0i = _mm512_mul_pd(nIi,sqr);
                        t1r = _mm512_mul_pd(k0a,ln);
                        t1i = _mm512_mul_pd(nIi,pi2);
                        cdiv_zmm8r8(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm8r8(t0r,t0i,divr,divi,&Kzr,&Kzi);
                      
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) peps0,
                                           const  double * __restrict __ATTR_ALIGN__(64) pmu0,
                                           const   double * __restrict __ATTR_ALIGN__(64) pEr,
                                           const   double * __restrict __ATTR_ALIGN__(64) pEi,
                                           const   double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           double * __restrict __ATTR_ALIGN__(64) Kzr,
                                           double * __restrict __ATTR_ALIGN__(64) Kzi) {

                        const register __m512d eps0 = _mm512_load_pd(&peps0[0]);
                        const register __m512d mu0  = _mm512_load_pd(&pmu0[0]);
                        const register __m512d Er   = _mm512_load_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_load_pd(&pEi[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d lna = _mm512_mul_pd(k0a,
                                                     _mm512_set1_pd(0.8905));
                        const register __m512d ln  = logkf(lna);
                        const __m512d sqr = _mm512_sqrt_pd(_mm512_div_pd(eps0,mu0));
                        const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                        __m512d t0r,t0i,t1r,t1i,divr,divi,resr,resi;
                        t0r = nIr;
                        t0i = _mm512_mul_pd(nIi,sqr);
                        t1r = _mm512_mul_pd(k0a,ln);
                        t1i = _mm512_mul_pd(nIi,pi2);
                        cdiv_zmm8r8(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_pd(&Kzr[0], resr);
                        _mm512_store_pd(&Kzi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm8r8_u(const  double * __restrict  peps0,
                                           const  double * __restrict  pmu0,
                                           const   double * __restrict  pEr,
                                           const   double * __restrict  pEi,
                                           const   double * __restrict  pk0a,
                                           double * __restrict  Kzr,
                                           double * __restrict  Kzi) {

                        const register __m512d eps0 = _mm512_loadu_pd(&peps0[0]);
                        const register __m512d mu0  = _mm512_loadu_pd(&pmu0[0]);
                        const register __m512d Er   = _mm512_loadu_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_loadu_pd(&pEi[0]);
                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d lna = _mm512_mul_pd(k0a,
                                                     _mm512_set1_pd(0.8905));
                        const register __m512d ln  = logkf(lna);
                        const __m512d sqr = _mm512_sqrt_pd(_mm512_div_pd(eps0,mu0));
                        const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                        __m512d t0r,t0i,t1r,t1i,divr,divi,resr,resi;
                        t0r = nIr;
                        t0i = _mm512_mul_pd(nIi,sqr);
                        t1r = _mm512_mul_pd(k0a,ln);
                        t1i = _mm512_mul_pd(nIi,pi2);
                        cdiv_zmm8r8(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_storeu_pd(&Kzr[0], resr);
                        _mm512_storeu_pd(&Kzi[0], resi);
                 }


                  /*
                          Surface currents (k0a << 1), for long cylinder (wire).
                          H-field cylinder axis parallel.
                          Formula 4.1-26
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm8r8(const __m512d Hr,
                                          const __m512d Hi,
                                          __m512d * __restrict Kphr,
                                          __m512d * __restrict Kphi) {

                        *Kphr = _mm512_mul_pd(nIi,Hr);
                        *Kphi = _mm512_mul_pd(nIi,Hi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) Hr,
                                            const double * __restrict __ATTR_ALIGN__(64) Hi,
                                           double * __restrict __ATTR_ALIGN__(64) Kphr,
                                          double * __restrict __ATTR_ALIGN__(64) Kphi) {

                        _mm512_store_pd(&Kphr[0] ,_mm512_mul_pd(nIi,_mm512_load_pd(&Hr[0]));
                        _mm512_store_pd(&Kphi[0] ,_mm512_mul_pd(nIi,_mm512_load_pd(&Hi[0]));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm8r8_u(const double * __restrict  Hr,
                                            const double * __restrict  Hi,
                                           double * __restrict  Kphr,
                                          double * __restrict Kphi) {

                        _mm512_storeu_pd(&Kphr[0] ,_mm512_mul_pd(nIi,_mm512_loadu_pd(&Hr[0]));
                        _mm512_storeu_pd(&Kphi[0] ,_mm512_mul_pd(nIi,_mm512_loadu_pd(&Hi[0]));
                 }


                   /*
                        The toal current along the wire.
                        Formula 4.1-27 

                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm8r8(const __m512d eps0,
                                         const __m512d mu0,
                                         const __m512d Er,
                                         const __m512d Ei,
                                         const __m512d k0a,
                                         const __m512d k0,
                                         __m512d * __restrict Izr,
                                         __m512d * __restrict Izi) {

                        const register __m512d lna = _mm512_mul_pd(k0a,
                                                     _mm512_set1_pd(0.8905));
                        const register __m512d _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
                        const register __m512d ln  = logkf(lna);
                        const __m512d sqr = _mm512_sqrt_pd(_mm512_div_pd(eps0,mu0));
                        const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                        __m512d t0r,t0i,t1r,t1i,divr,divi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_pd(nIi,sqr);
                        t2r = _mm512_mul_pd(_2pi,Er);
                        t1r = _mm512_mul_pd(k0,ln);
                        t2i = _mm512_mul_pd(_2pi,Ei);
                        t1i = _mm512_mul_pd(nIi,pi2);
                        cdiv_zmm8r8(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm8r8(t0r,t0i,divr,divi,&Izr,&Izi);
                      
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) peps0,
                                           const  double * __restrict __ATTR_ALIGN__(64) pmu0,
                                           const   double * __restrict __ATTR_ALIGN__(64) pEr,
                                           const   double * __restrict __ATTR_ALIGN__(64) pEi,
                                           const   double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const   double * __restrict __ATTR_ALIGN__(64) pk0,
                                           double * __restrict __ATTR_ALIGN__(64) Izr,
                                           double * __restrict __ATTR_ALIGN__(64) Izi) {

                        const register __m512d eps0 = _mm512_load_pd(&peps0[0]);
                        const register __m512d _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
                        const register __m512d mu0  = _mm512_load_pd(&pmu0[0]);
                        const register __m512d sqr = _mm512_sqrt_pd(_mm512_div_pd(eps0,mu0));
                        const register __m512d Er   = _mm512_load_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_load_pd(&pEi[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d lna = _mm512_mul_pd(k0a,
                                                     _mm512_set1_pd(0.8905f));
                        const register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        const register __m512d ln  = logkf(lna);
                        const register __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                        __m512d t0r,t0i,t1r,t1i,divr,divi,resr,resi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_pd(nIi,sqr);
                        t1r = _mm512_mul_pd(k0,ln);
                        t2r = _mm512_mul_pd(_2pi,Er);
                        t1i = _mm512_mul_pd(nIi,pi2);
                        t2i = _mm512_mul_pd(_2pi,Ei);
                        cdiv_zmm8r8(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_pd(&Izr[0], resr);
                        _mm512_store_pd(&Izi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm8r8_u(const  double * __restrict  peps0,
                                           const  double * __restrict  pmu0,
                                           const   double * __restrict  pEr,
                                           const   double * __restrict  pEi,
                                           const   double * __restrict  pk0a,
                                           const   double * __restrict  pk0,
                                           double * __restrict  Izr,
                                           double * __restrict  Izi) {

                        const register __m512d eps0 = _mm512_load_pd(&peps0[0]);
                        const register __m512d _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
                        const register __m512d mu0  = _mm512_load_pd(&pmu0[0]);
                        const register __m512d sqr = _mm512_sqrt_pd(_mm512_div_pd(eps0,mu0));
                        const register __m512d Er   = _mm512_load_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_load_pd(&pEi[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d lna = _mm512_mul_pd(k0a,
                                                     _mm512_set1_pd(0.8905f));
                        const register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        const register __m512d ln  = logkf(lna);
                        const register __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                        __m512d t0r,t0i,t1r,t1i,divr,divi,resr,resi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_pd(nIi,sqr);
                        t1r = _mm512_mul_pd(k0,ln);
                        t2r = _mm512_mul_pd(_2pi,Er);
                        t1i = _mm512_mul_pd(nIi,pi2);
                        t2i = _mm512_mul_pd(_2pi,Ei);
                        cdiv_zmm8r8(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_pd(&Izr[0], resr);
                        _mm512_store_pd(&Izi[0], resi);
                 }


                   /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        (0<<phi<pi/2, k0a > 2)
                        Electric-field.
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm8r8(const __m512d phi2,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d k0a,
                                         const __m512d Er,
                                         const __m512d Ei,
                                         __m512d * __restrict EOr,
                                         __m512d * __restrict EOi) {

                        register __m512d t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512d t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512d k0as,fac,earg,_2a,cer,cei,exr,exi;
                        register __m512d t2r,t2i;
                        const register __m512d c0 = _mm512_set1_pd(0.375);
                        cosf2 = xcosf(phi2);
                        const register __m512d c1 = _mm512_set1_pd(0.1171875);
                        cos2f2 = _mm512_mul_pd(cosf2,cosf2);
                        const register __m512d c2 = _mm512_set1_pd(4.0);
                        cos4f2 = _mm512_mul_pd(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d c3 = _mm512_set1_pd(8.0);
                        _2r    = _mm512_add_pd(r,r);
                        _2a    = _mm512_add_pd(a,a);
                        const register __m512d c4 = _mm512_set1_pd(33.0);
                        k0as   = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c5 = _mm512_set1_pd(5.0);
                        t0     = _mm512_mul_pd(a,cosf2);
                        const register __m512d c6 = _mm512_set1_pd(1.0);
                        t1     = _mm512_div_pd(t0,_2r);
                        fac    = _mm512_sqrt_pd(t1);
                        earg   = _mm512_mul_pd(k0,
                                          _mm512_sub_pd(r,_mm512_mul_pd(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_pd(Ii,earg);
                        cexp_zmm8r8(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_pd(cos2f2);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_pd(t3,c0);//taken t3
                        t0     = _mm512_mul_pd(c2,_mm512_mul_pd(k0as,cos2f2));
                        t4     = _mm512_rcp14_pd(t0);//taken t4
                        t1     = _mm512_mul_pd(c3,cos2f2);
                        t2     = _mm512_sub_pd(c1,_mm512_div_pd(c4,t1)); // t2 taken
                        t0     = _mm512_div_pd(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_pd(Ii,_mm512_mul_pd(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_pd(c6,t2i);
                        cmul_zmm8r8(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_pd(t4,_mm512_add_pd(t2,t0));//taken t5
                        t2r    = _mm512_add_pd(t2r,t5);
                        t2i    = _mm512_add_pd(t2i,t5);
                        cmul_zmm8r8(exr,exi,t2r,t2i,&t0r,&t0i);
                        *EOr = t0r;
                        *EOi = t0i;
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pEr,
                                           const double * __restrict __ATTR_ALIGN__(64) pEi,
                                           double * __restrict __ATTR_ALIGN__(64) EOr,
                                           double * __restrict __ATTR_ALIGN__(64) EOi) {

                        const register __m512d phi2 = _mm512_load_pd(&pphi2[0]);
                        const register __m512d a    = _mm512_load_pd(&pa[0]);
                        const register __m512d r    = _mm512_load_pd(&pr[0]);
                        const register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d Er   = _mm512_load_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_load_pd(&pEi[0]);
                        register __m512d t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512d t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512d k0as,fac,earg,_2a,cer,cei,exr,exi;;
                        register __m512d t2r,t2i;
                        const register __m512d c0 = _mm512_set1_pd(0.375f;
                        cosf2 = xcosf(phi2);
                        const register __m512d c1 = _mm512_set1_pd(0.1171875);
                        cos2f2 = _mm512_mul_pd(cosf2,cosf2);
                        const register __m512d c2 = _mm512_set1_pd(4.0);
                        cos4f2 = _mm512_mul_pd(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d c3 = _mm512_set1_pd(8.0);
                        _2r    = _mm512_add_pd(r,r);
                        _2a    = _mm512_add_pd(a,a);
                        const register __m512d c4 = _mm512_set1_pd(33.0);
                        k0as   = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c5 = _mm512_set1_pd(5.0);
                        t0     = _mm512_mul_pd(a,cosf2);
                        const register __m512d c6 = _mm512_set1_pd(1.0);
                        t1     = _mm512_div_pd(t0,_2r);
                        fac    = _mm512_sqrt_pd(t1);
                        earg   = _mm512_mul_pd(k0,
                                          _mm512_sub_pd(r,_mm512_mul_pd(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_pd(Ii,earg);
                        cexp_zmm8r8(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_pd(cos2f2);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_pd(t3,c0);//taken t3
                        t0     = _mm512_mul_pd(c2,_mm512_mul_pd(k0as,cos2f2));
                        t4     = _mm512_rcp14_pd(t0);//taken t4
                        t1     = _mm512_mul_pd(c3,cos2f2);
                        t2     = _mm512_sub_pd(c1,_mm512_div_pd(c4,t1)); // t2 taken
                        t0     = _mm512_div_pd(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_pd(Ii,_mm512_mul_pd(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_pd(c6,t2i);
                        cmul_zmm8r8(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_pd(t4,_mm512_add_pd(t2,t0));//taken t5
                        t2r    = _mm512_add_pd(t2r,t5);
                        t2i    = _mm512_add_pd(t2i,t5);
                        cmul_zmm8r8(exr,exi,t2r,t2i,&t0r,&t0i);
                        _mm512_store_pd(&EOr[0], t0r);
                        _mm512_store_pd(&EOi[0], t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm8r8_u(const double * __restrict  pphi2,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0,
                                           const double * __restrict  pk0a,
                                           const double * __restrict  pEr,
                                           const double * __restrict  pEi,
                                           double * __restrict  EOr,
                                           double * __restrict  EOi) {

                        const register __m512d phi2 = _mm512_loadu_pd(&pphi2[0]);
                        const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        const register __m512d r    = _mm512_loadu_pd(&pr[0]);
                        const register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d Er   = _mm512_loadu_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_loadu_pd(&pEi[0]);
                        register __m512d t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512d t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512d k0as,fac,earg,_2a,cer,cei,exr,exi;;
                        register __m512d t2r,t2i;
                        const register __m512d c0 = _mm512_set1_pd(0.375);
                        cosf2 = xcosf(phi2);
                        const register __m512d c1 = _mm512_set1_pd(0.1171875);
                        cos2f2 = _mm512_mul_pd(cosf2,cosf2);
                        const register __m512d c2 = _mm512_set1_pd(4.0);
                        cos4f2 = _mm512_mul_pd(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d c3 = _mm512_set1_pd(8.0);
                        _2r    = _mm512_add_pd(r,r);
                        _2a    = _mm512_add_pd(a,a);
                        const register __m512d c4 = _mm512_set1_pd(33.0);
                        k0as   = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c5 = _mm512_set1_pd(5.0);
                        t0     = _mm512_mul_pd(a,cosf2);
                        const register __m512d c6 = _mm512_set1_pd(1.0);
                        t1     = _mm512_div_pd(t0,_2r);
                        fac    = _mm512_sqrt_pd(t1);
                        earg   = _mm512_mul_pd(k0,
                                          _mm512_sub_pd(r,_mm512_mul_pd(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_pd(Ii,earg);
                        cexp_zmm8r8(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_pd(cos2f2);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_pd(t3,c0);//taken t3
                        t0     = _mm512_mul_pd(c2,_mm512_mul_pd(k0as,cos2f2));
                        t4     = _mm512_rcp14_pd(t0);//taken t4
                        t1     = _mm512_mul_pd(c3,cos2f2);
                        t2     = _mm512_sub_pd(c1,_mm512_div_pd(c4,t1)); // t2 taken
                        t0     = _mm512_div_pd(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_pd(Ii,_mm512_mul_pd(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_pd(c6,t2i);
                        cmul_zmm8r8(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_pd(t4,_mm512_add_pd(t2,t0));//taken t5
                        t2r    = _mm512_add_pd(t2r,t5);
                        t2i    = _mm512_add_pd(t2i,t5);
                        cmul_zmm8r8(exr,exi,t2r,t2i,&t0r,&t0i);
                        _mm512_storeu_pd(&EOr[0], t0r);
                        _mm512_storeu_pd(&EOi[0], t0i);
                 }


                     /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        (0<<phi<pi/2, k0a > 2)
                        Magnetic-field.
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm8r8(const __m512d phi2,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d k0a,
                                         const __m512d Hr,
                                         const __m512d Hi,
                                         __m512d * __restrict HOr,
                                         __m512d * __restrict HOi) {

                        register __m512d t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512d t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512d k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512d t2r,t2i;
                        const register __m512d c0 = _mm512_set1_pd(0.375);
                        cosf2 = xcosf(phi2);
                        const register __m512d c1 = _mm512_set1_pd(0.1171875);
                        cos2f2 = _mm512_mul_pd(cosf2,cosf2);
                        const register __m512d c2 = _mm512_set1_pd(4.0);
                        cos4f2 = _mm512_mul_pd(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d c3 = _mm512_set1_pd(8.0);
                        _2r    = _mm512_add_pd(r,r);
                        _2a    = _mm512_add_pd(a,a);
                        const register __m512d c4 = _mm512_set1_pd(33.0);
                        k0as   = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c5 = _mm512_set1_pd(7.0);
                        t0     = _mm512_mul_pd(a,cosf2);
                        const register __m512d c6 = _mm512_set1_pd(1.0);
                        t1     = _mm512_div_pd(t0,_2r);
                        fac    = _mm512_sqrt_pd(t1);
                        earg   = _mm512_mul_pd(k0,
                                          _mm512_sub_pd(r,_mm512_mul_pd(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_pd(Ii,earg);
                        cexp_zmm8r8(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_pd(cos2f2);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_pd(t3,c0);//taken t3
                        t0     = _mm512_mul_pd(c2,_mm512_mul_pd(k0as,cos2f2));
                        t4     = _mm512_rcp14_pd(t0);//taken t4
                        t1     = _mm512_mul_pd(c3,cos2f2);
                        t2     = _mm512_add_pd(c1,_mm512_div_pd(c4,t1)); // t2 taken
                        t0     = _mm512_div_pd(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_pd(Ii,_mm512_mul_pd(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_pd(c6,t2i);
                        cmul_zmm8r8(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_pd(t4,_mm512_sub_pd(t2,t0));//taken t5
                        t2r    = _mm512_add_pd(t2r,t5);
                        t2i    = _mm512_add_pd(t2i,t5);
                        cmul_zmm8r8(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        *HOr = t0r;
                        *HOi = t0i;
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pHr,
                                           const double * __restrict __ATTR_ALIGN__(64) pHi,
                                           double * __restrict __ATTR_ALIGN__(64) HOr,
                                           double * __restrict __ATTR_ALIGN__(64) HOi) {

                        const register __m512d phi2 = _mm512_load_pd(&pphi2[0]);
                        const register __m512d a    = _mm512_load_pd(&pa[0]);
                        const register __m512d r    = _mm512_load_pd(&pr[0]);
                        const register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d Hr   = _mm512_load_pd(&pHr[0]);
                        const register __m512d Hi   = _mm512_load_pd(&pHi[0]);
                        register __m512d t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512d t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512d k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512d t2r,t2i;
                        const register __m512d c0 = _mm512_set1_pd(0.375);
                        cosf2 = xcosf(phi2);
                        const register __m512d c1 = _mm512_set1_pd(0.1171875);
                        cos2f2 = _mm512_mul_pd(cosf2,cosf2);
                        const register __m512d c2 = _mm512_set1_pd(4.0);
                        cos4f2 = _mm512_mul_pd(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d c3 = _mm512_set1_pd(8.0);
                        _2r    = _mm512_add_pd(r,r);
                        _2a    = _mm512_add_pd(a,a);
                        const register __m512d c4 = _mm512_set1_pd(33.0);
                        k0as   = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c5 = _mm512_set1_pd(7.0);
                        t0     = _mm512_mul_pd(a,cosf2);
                        const register __m512d c6 = _mm512_set1_pd(1.0);
                        t1     = _mm512_div_pd(t0,_2r);
                        fac    = _mm512_sqrt_pd(t1);
                        earg   = _mm512_mul_pd(k0,
                                          _mm512_sub_pd(r,_mm512_mul_pd(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_pd(Ii,earg);
                        cexp_zmm8r8(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_pd(cos2f2);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_pd(t3,c0);//taken t3
                        t0     = _mm512_mul_pd(c2,_mm512_mul_pd(k0as,cos2f2));
                        t4     = _mm512_rcp14_pd(t0);//taken t4
                        t1     = _mm512_mul_pd(c3,cos2f2);
                        t2     = _mm512_add_pd(c1,_mm512_div_pd(c4,t1)); // t2 taken
                        t0     = _mm512_div_pd(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_pd(Ii,_mm512_mul_pd(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_pd(c6,t2i);
                        cmul_zmm8r8(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_pd(t4,_mm512_sub_pd(t2,t0));//taken t5
                        t2r    = _mm512_add_pd(t2r,t5);
                        t2i    = _mm512_add_pd(t2i,t5);
                        cmul_zmm8r8(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        _mm512_store_pd(&HOr[0], t0r);
                        _mm512_store_pd(&HOi[0], t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm8r8_u(const double * __restrict  pphi2,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0,
                                           const double * __restrict  pk0a,
                                           const double * __restrict  pHr,
                                           const double * __restrict  pHi,
                                           double * __restrict  HOr,
                                           double * __restrict  HOi) {

                        const register __m512d phi2 = _mm512_loadu_pd(&pphi2[0]);
                        const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        const register __m512d r    = _mm512_loadu_pd(&pr[0]);
                        const register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d Hr   = _mm512_loadu_pd(&pHr[0]);
                        const register __m512d Hi   = _mm512_loadu_pd(&pHi[0]);
                        register __m512d t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512d t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512d k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512d t2r,t2i;
                        const register __m512d c0 = _mm512_set1_pd(0.375);
                        cosf2 = xcosf(phi2);
                        const register __m512d c1 = _mm512_set1_pd(0.1171875);
                        cos2f2 = _mm512_mul_pd(cosf2,cosf2);
                        const register __m512d c2 = _mm512_set1_pd(4.0);
                        cos4f2 = _mm512_mul_pd(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d c3 = _mm512_set1_pd(8.0);
                        _2r    = _mm512_add_pd(r,r);
                        _2a    = _mm512_add_pd(a,a);
                        const register __m512d c4 = _mm512_set1_pd(33.0);
                        k0as   = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c5 = _mm512_set1_pd(7.0);
                        t0     = _mm512_mul_pd(a,cosf2);
                        const register __m512d c6 = _mm512_set1_pd(1.0);
                        t1     = _mm512_div_pd(t0,_2r);
                        fac    = _mm512_sqrt_pd(t1);
                        earg   = _mm512_mul_pd(k0,
                                          _mm512_sub_pd(r,_mm512_mul_pd(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_pd(Ii,earg);
                        cexp_zmm8r8(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_pd(cos2f2);
                        cmul_zmm8r8(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_pd(t3,c0);//taken t3
                        t0     = _mm512_mul_pd(c2,_mm512_mul_pd(k0as,cos2f2));
                        t4     = _mm512_rcp14_pd(t0);//taken t4
                        t1     = _mm512_mul_pd(c3,cos2f2);
                        t2     = _mm512_add_pd(c1,_mm512_div_pd(c4,t1)); // t2 taken
                        t0     = _mm512_div_pd(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_pd(Ii,_mm512_mul_pd(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_pd(c6,t2i);
                        cmul_zmm8r8(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_pd(t4,_mm512_sub_pd(t2,t0));//taken t5
                        t2r    = _mm512_add_pd(t2r,t5);
                        t2i    = _mm512_add_pd(t2i,t5);
                        cmul_zmm8r8(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        _mm512_storeu_pd(&HOr[0], t0r);
                        _mm512_storeu_pd(&HOi[0], t0i);
                 }


                 
                 /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        (0<<phi<pi/2, k0a > 2)
                        Electric-field.
                        Formula 4.1-30
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm8r8(const __m512d Er,
                                         const __m512d Ei,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d k0a,
                                         const __m512d phi,
                                         __m512d * __restrict ECr,
                                         __m512d * __restrict ECi) {

                        register __m512d e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512d ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512d Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        const  __m512d k0ai16 = _mm512_pow_pd(k0a,
                                                         _mm512_set1_pd(0.166666666666666666666666666667));
                        const register __m512d k0apaphi = _mm512_fmadd_pd(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm8r8(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_pd(k0ai16);
                        const register __m512d k0apsphi = _mm512_fmsub_pd(k0a,PI,phi);
                        const  __m512d c0   = _mm512_set1_pd(0.910721);
                        const  register __m512d k0rp12 = _mm512_fmadd_pd(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm8r8(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512d c0r  = _mm512_set1_pd(0.9358135);
                        sqr    = _mm512_div_pd(a,_mm512_add_pd(r,r);
                        const  __m512d c0i  = _mm512_set1_pd(1.607129);
                        sqr    = _mm512_sqrt_pd(sqr);
                        const  __m512d c1r  = _mm512_set1_pd(0.057397);
                        Etr    = _mm512_mul_pd(Er,sqr);// first complex term
                        Eti    = _mm512_mul_pd(Ei,sqr);// first complex term
                        const  __m512d c1i  = _mm512_set1_pd(-0.0994145f);
                        const  __m512d _1   = _mm512_set1_pd(1.0f);
                        const  __m512d k0an23 = _mm512_pow_pd(k0a,
                                                          _mm512_set1_pd(0.666666666666666666666666666667));
                        k0an23 = _mm512_rcp14_pd(k0an23);
                        const __m512d k0an43= _mm512_pow_pd(k0a,
                                                        _mm512_set1_pd(1.333333333333333333333333333333));
                        k0an43 = _mm512_rcp14_pd(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_pd(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_pd(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_pd(c1r,k0an43);
                        t1i = _mm512_mul_pd(c1i,k0an43);
                        tmp1r = _mm512_sub_pd(t0r,t1r);
                        tmp1i = _mm512_sub_pd(t0i,t1i);
                        cmul_zmm8r8(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm8r8(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_pd();
                        t0i = _mm512_setzero_pd();
                        cmul_zmm8r8(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        *ECr = _mm512_add_pd(t0r,tmp3r);
                        *ECi = _mm512_add_pd(t0i,tmp3i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pEr,
                                           const double * __restrict __ATTR_ALIGN__(64) pEi,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pphi,
                                           double * __restrict __ATTR_ALIGN__(64) ECr,
                                           double * __restrict __ATTR_ALIGN__(64) ECi) {

                        const register __m512d phi  = _mm512_load_pd(&pphi[0]);
                        const register __m512d a    = _mm512_load_pd(&pa[0]);
                        const register __m512d r    = _mm512_load_pd(&pr[0]);
                        const register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d Er   = _mm512_load_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_load_pd(&pEi[0]);
                        register __m512d e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512d ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512d Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        const  __m512d k0ai16 = _mm512_pow_pd(k0a,
                                                         _mm512_set1_pd(0.166666666666666666666666666667));
                        const register __m512d k0apaphi = _mm512_fmadd_pd(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm8r8(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_pd(k0ai16);
                        const register __m512d k0apsphi = _mm512_fmsub_pd(k0a,PI,phi);
                        const  __m512d c0   = _mm512_set1_pd(0.910721);
                        const  register __m512d k0rp12 = _mm512_fmadd_pd(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm8r8(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512d c0r  = _mm512_set1_pd(0.9358135);
                        sqr    = _mm512_div_pd(a,_mm512_add_pd(r,r);
                        const  __m512d c0i  = _mm512_set1_pd(1.607129);
                        sqr    = _mm512_sqrt_pd(sqr);
                        const  __m512d c1r  = _mm512_set1_pd(0.057397);
                        Etr    = _mm512_mul_pd(Er,sqr);// first complex term
                        Eti    = _mm512_mul_pd(Ei,sqr);// first complex term
                        const  __m512d c1i  = _mm512_set1_pd(-0.0994145);
                        const  __m512d _1   = _mm512_set1_pd(1.0);
                        const  __m512d k0an23 = _mm512_pow_pd(k0a,
                                                          _mm512_set1_pd(0.666666666666666666666666666667));
                        k0an23 = _mm512_rcp14_pd(k0an23);
                        const __m512d k0an43= _mm512_pow_pd(k0a,
                                                        _mm512_set1_pd(1.333333333333333333333333333333));
                        k0an43 = _mm512_rcp14_pd(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_pd(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_pd(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_pd(c1r,k0an43);
                        t1i = _mm512_mul_pd(c1i,k0an43);
                        tmp1r = _mm512_sub_pd(t0r,t1r);
                        tmp1i = _mm512_sub_pd(t0i,t1i);
                        cmul_zmm8r8(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm8r8(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_pd();
                        t0i = _mm512_setzero_pd();
                        cmul_zmm8r8(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_store_pd(&ECr[0] ,_mm512_add_pd(t0r,tmp3r));
                        _mm512_store_pd(&ECi[0] ,_mm512_add_pd(t0i,tmp3i));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm8r8_a(const double * __restrict  pEr,
                                           const double * __restrict  pEi,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0,
                                           const double * __restrict  pk0a,
                                           const double * __restrict  pphi,
                                           double * __restrict  ECr,
                                           double * __restrict ECi) {

                        const register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                        const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        const register __m512d r    = _mm512_loadu_pd(&pr[0]);
                        const register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d Er   = _mm512_loadu_pd(&pEr[0]);
                        const register __m512d Ei   = _mm512_loadu_pd(&pEi[0]);
                        register __m512d e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512d ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512d Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        const  __m512d k0ai16 = _mm512_pow_pd(k0a,
                                                         _mm512_set1_pd(0.166666666666666666666666666667));
                        const register __m512d k0apaphi = _mm512_fmadd_pd(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm8r8(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_pd(k0ai16);
                        const register __m512d k0apsphi = _mm512_fmsub_pd(k0a,PI,phi);
                        const  __m512d c0   = _mm512_set1_pd(0.910721);
                        const  register __m512d k0rp12 = _mm512_fmadd_pd(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm8r8(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512d c0r  = _mm512_set1_pd(0.9358135);
                        sqr    = _mm512_div_pd(a,_mm512_add_pd(r,r);
                        const  __m512d c0i  = _mm512_set1_pd(1.607129);
                        sqr    = _mm512_sqrt_pd(sqr);
                        const  __m512d c1r  = _mm512_set1_pd(0.057397);
                        Etr    = _mm512_mul_pd(Er,sqr);// first complex term
                        Eti    = _mm512_mul_pd(Ei,sqr);// first complex term
                        const  __m512d c1i  = _mm512_set1_pd(-0.0994145);
                        const  __m512d _1   = _mm512_set1_pd(1.0f);
                        const  __m512d k0an23 = _mm512_pow_pd(k0a,
                                                          _mm512_set1_pd(0.666666666666666666666666666667));
                        k0an23 = _mm512_rcp14_pd(k0an23);
                        const __m512d k0an43= _mm512_pow_pd(k0a,
                                                        _mm512_set1_pd(1.333333333333333333333333333333));
                        k0an43 = _mm512_rcp14_pd(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_pd(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_pd(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_pd(c1r,k0an43);
                        t1i = _mm512_mul_pd(c1i,k0an43);
                        tmp1r = _mm512_sub_pd(t0r,t1r);
                        tmp1i = _mm512_sub_pd(t0i,t1i);
                        cmul_zmm8r8(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm8r8(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_pd();
                        t0i = _mm512_setzero_pd();
                        cmul_zmm8r8(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_storeu_pd(&ECr[0] ,_mm512_add_pd(t0r,tmp3r));
                        _mm512_storeu_pd(&ECi[0] ,_mm512_add_pd(t0i,tmp3i));
                 }


                    /*
                        Approximation for upper-middle and high-frequency region
                        (k0a > 2).
                        Bistatic creeping wave approximation for resonance region
                        valid only for (0<<phi<pi/2, k0a > 2)
                        Magnetic-field.
                        Formula 4.1-32
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm8r8(const __m512d Hr,
                                         const __m512d Hi,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d k0a,
                                         const __m512d phi,
                                         __m512d * __restrict HCr,
                                         __m512d * __restrict HCi) {

                        register __m512d e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512d ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512d Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        const  __m512d k0ai16 = _mm512_pow_pd(k0a,
                                                         _mm512_set1_pd(0.166666666666666666666666666667f));
                        const register __m512d k0apaphi = _mm512_fmadd_pd(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm8r8(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_pd(k0ai16);
                        const register __m512d k0apsphi = _mm512_fmsub_pd(k0a,PI,phi);
                        const  __m512d c0   = _mm512_set1_pd(1.531915f);
                        const  register __m512d k0rp12 = _mm512_fmadd_pd(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm8r8(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512d c0r  = _mm512_set1_pd(0.404308f);
                        sqr    = _mm512_div_pd(a,_mm512_add_pd(r,r);
                        const  __m512d c0i  = _mm512_set1_pd(0.70028f);
                        sqr    = _mm512_sqrt_pd(sqr);
                        const  __m512d c1r  = _mm512_set1_pd(0.072732f);
                        Etr    = _mm512_mul_pd(Er,sqr);// first complex term
                        Eti    = _mm512_mul_pd(Ei,sqr);// first complex term
                        const  __m512d c1i  = _mm512_set1_pd(-0.1259755f);
                        const  __m512d _1   = _mm512_set1_pd(1.0f);
                        const  __m512d k0an23 = _mm512_pow_pd(k0a,
                                                          _mm512_set1_pd(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_pd(k0an23);
                        const __m512d k0an43= _mm512_pow_pd(k0a,
                                                        _mm512_set1_pd(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_pd(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_pd(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_pd(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_pd(c1r,k0an43);
                        t1i = _mm512_mul_pd(c1i,k0an43);
                        tmp1r = _mm512_sub_pd(t0r,t1r);
                        tmp1i = _mm512_sub_pd(t0i,t1i);
                        cmul_zmm8r8(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm8r8(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_pd();
                        t0i = _mm512_setzero_pd();
                        cmul_zmm8r8(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        *HCr = _mm512_add_pd(t0r,tmp3r);
                        *HCi = _mm512_add_pd(t0i,tmp3i);
                 }


                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64)  pHr,
                                           const  double * __restrict __ATTR_ALIGN__(64)  pHi,
                                           const  double * __restrict __ATTR_ALIGN__(64)  pa,
                                           const  double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const  double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const  double * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const  double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           double * __restrict __ATTR_ALIGN__(64)  HCr,
                                           double * __restrict __ATTR_ALIGN__(64)  HCi) {

                        const register __m512d phi  = _mm512_load_pd(&pphi[0]);
                        const register __m512d a    = _mm512_load_pd(&pa[0]);
                        const register __m512d r    = _mm512_load_pd(&pr[0]);
                        const register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d Hr   = _mm512_load_pd(&pHr[0]);
                        const register __m512d Hi   = _mm512_load_pd(&pHi[0]);
                        register __m512d e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512d ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512d Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        const  __m512d k0ai16 = _mm512_pow_pd(k0a,
                                                         _mm512_set1_pd(0.166666666666666666666666666667f));
                        const register __m512d k0apaphi = _mm512_fmadd_pd(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm8r8(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_pd(k0ai16);
                        const register __m512d k0apsphi = _mm512_fmsub_pd(k0a,PI,phi);
                        const  __m512d c0   = _mm512_set1_pd(1.531915f);
                        const  register __m512d k0rp12 = _mm512_fmadd_pd(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm8r8(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512d c0r  = _mm512_set1_pd(0.404308f);
                        sqr    = _mm512_div_pd(a,_mm512_add_pd(r,r);
                        const  __m512d c0i  = _mm512_set1_pd(0.70028f);
                        sqr    = _mm512_sqrt_pd(sqr);
                        const  __m512d c1r  = _mm512_set1_pd(0.072732f);
                        Etr    = _mm512_mul_pd(Er,sqr);// first complex term
                        Eti    = _mm512_mul_pd(Ei,sqr);// first complex term
                        const  __m512d c1i  = _mm512_set1_pd(-0.1259755f);
                        const  __m512d _1   = _mm512_set1_pd(1.0f);
                        const  __m512d k0an23 = _mm512_pow_pd(k0a,
                                                          _mm512_set1_pd(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_pd(k0an23);
                        const __m512d k0an43= _mm512_pow_pd(k0a,
                                                        _mm512_set1_pd(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_pd(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_pd(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_pd(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_pd(c1r,k0an43);
                        t1i = _mm512_mul_pd(c1i,k0an43);
                        tmp1r = _mm512_sub_pd(t0r,t1r);
                        tmp1i = _mm512_sub_pd(t0i,t1i);
                        cmul_zmm8r8(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm8r8(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_pd();
                        t0i = _mm512_setzero_pd();
                        cmul_zmm8r8(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_store_pd(&HCr[0] ,_mm512_add_pd(t0r,tmp3r));
                        _mm512_store_pd(&HCi[0] ,_mm512_add_pd(t0i,tmp3i));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm8r8_u(const  double * __restrict  pHr,
                                           const  double * __restrict  pHi,
                                           const  double * __restrict  pa,
                                           const  double * __restrict  pr,
                                           const  double * __restrict  pk0,
                                           const  double * __restrict  pk0a,
                                           const  double * __restrict  pphi,
                                           double * __restrict   HCr,
                                           double * __restrict   HCi) {

                        const register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                        const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        const register __m512d r    = _mm512_loadu_pd(&pr[0]);
                        const register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d Hr   = _mm512_loadu_pd(&pHr[0]);
                        const register __m512d Hi   = _mm512_loadu_pd(&pHi[0]);
                        register __m512d e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512d ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512d Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        const  __m512d k0ai16 = _mm512_pow_pd(k0a,
                                                         _mm512_set1_pd(0.166666666666666666666666666667));
                        const register __m512d k0apaphi = _mm512_fmadd_pd(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm8r8(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_pd(k0ai16);
                        const register __m512d k0apsphi = _mm512_fmsub_pd(k0a,PI,phi);
                        const  __m512d c0   = _mm512_set1_pd(1.531915);
                        const  register __m512d k0rp12 = _mm512_fmadd_pd(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm8r8(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512d c0r  = _mm512_set1_pd(0.404308);
                        sqr    = _mm512_div_pd(a,_mm512_add_pd(r,r);
                        const  __m512d c0i  = _mm512_set1_pd(0.70028);
                        sqr    = _mm512_sqrt_pd(sqr);
                        const  __m512d c1r  = _mm512_set1_pd(0.072732f);
                        Etr    = _mm512_mul_pd(Er,sqr);// first complex term
                        Eti    = _mm512_mul_pd(Ei,sqr);// first complex term
                        const  __m512d c1i  = _mm512_set1_pd(-0.1259755f);
                        const  __m512d _1   = _mm512_set1_pd(1.0f);
                        const  __m512d k0an23 = _mm512_pow_pd(k0a,
                                                          _mm512_set1_pd(0.666666666666666666666666666667));
                        k0an23 = _mm512_rcp14_pd(k0an23);
                        const __m512d k0an43= _mm512_pow_pd(k0a,
                                                        _mm512_set1_pd(1.333333333333333333333333333333));
                        k0an43 = _mm512_rcp14_pd(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_pd(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_pd(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_pd(c1r,k0an43);
                        t1i = _mm512_mul_pd(c1i,k0an43);
                        tmp1r = _mm512_sub_pd(t0r,t1r);
                        tmp1i = _mm512_sub_pd(t0i,t1i);
                        cmul_zmm8r8(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm8r8(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_pd();
                        t0i = _mm512_setzero_pd();
                        cmul_zmm8r8(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_storeu_pd(&HCr[0] ,_mm512_add_pd(t0r,tmp3r));
                        _mm512_storeu_pd(&HCi[0] ,_mm512_add_pd(t0i,tmp3i));
                 }


                   /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Optical wave component e-field, formula 4.1-33
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm8r8(const __m512d Er,
                                         const __m512d Ei,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0a,
                                         __m512d * __restrict EOr,
                                         __m512d * __restrict EOi) {

                        register __m512d _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512d ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        const __m512d _16= _mm512_set1_pd(16.0);
                        _2r             = _mm512_add_pd(r,r);
                        _2k0a           = _mm512_add_pd(k0a,k0a);
                        const __m512d _5 = _mm512_set1_pd(5.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d c0 = _mm512_set1_pd(127.0);
                        t0              = _mm512_sub_pd(k0r,_2k0a);
                        ear             = Ir;
                        const __m512d c1 = _mm512_set1_pd(512.0);
                        eai             = t0;
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        t1              = _mm512_div_pd(a,_2r);
                        t2              = _mm512_sqrt_pd(t1);
                        facr            = _mm512_mul_pd(Er,t2);
                        faci            = _mm512_mul_pd(Ei,t2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(_16,k0a);
                        t2              = _mm512_mul_pd(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_pd(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_pd(_1,t0i);
                        t3              = _mm512_div_pd(c0,t2);
                        t0r             = _mm512_add_pd(t3,t0r);
                        t0i             = _mm512_add_pd(t3,t0i);
                        cmul_zmm8r8(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm8r8(cer1,cei1,t0r,t0i,*EOr,*EOi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pEr,
                                           const double * __restrict __ATTR_ALIGN__(64) pEi,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           double * __restrict __ATTR_ALIGN__(64) EOr,
                                           double * __restrict __ATTR_ALIGN__(64) EOi) {

                        register __m512d Er = _mm512_load_pd(&pEr[0]);
                        register __m512d Ei = _mm512_load_pd(&pEi[0]);
                        register __m512d a  = _mm512_load_pd(&pa[0]);
                        register __m512d r  = _mm512_load_pd(&pr[0]);
                        register __m512d k0a= _mm512_load_pd(&pk0a[0]);
                        register __m512d _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512d ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        register __m512d resr,resi;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        const __m512d _16= _mm512_set1_pd(16.0);
                        _2r             = _mm512_add_pd(r,r);
                        _2k0a           = _mm512_add_pd(k0a,k0a);
                        const __m512d _5 = _mm512_set1_pd(5.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d c0 = _mm512_set1_pd(127.0);
                        t0              = _mm512_sub_pd(k0r,_2k0a);
                        ear             = Ir;
                        const __m512d c1 = _mm512_set1_pd(512.0);
                        eai             = t0;
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        t1              = _mm512_div_pd(a,_2r);
                        t2              = _mm512_sqrt_pd(t1);
                        facr            = _mm512_mul_pd(Er,t2);
                        faci            = _mm512_mul_pd(Ei,t2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(_16,k0a);
                        t2              = _mm512_mul_pd(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_pd(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_pd(_1,t0i);
                        t3              = _mm512_div_pd(c0,t2);
                        t0r             = _mm512_add_pd(t3,t0r);
                        t0i             = _mm512_add_pd(t3,t0i);
                        cmul_zmm8r8(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm8r8(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_store_pd(&EOr[0],resr);
                        _mm512_store_pd(&EOi[0],resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm8r8_u(const double * __restrict  pEr,
                                           const double * __restrict  pEi,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0a,
                                           double * __restrict  EOr,
                                           double * __restrict  EOi) {

                        register __m512d Er = _mm512_loadu_pd(&pEr[0]);
                        register __m512d Ei = _mm512_loadu_pd(&pEi[0]);
                        register __m512d a  = _mm512_loadu_pd(&pa[0]);
                        register __m512d r  = _mm512_loadu_pd(&pr[0]);
                        register __m512d k0a= _mm512_loadu_pd(&pk0a[0]);
                        register __m512d _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512d ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        register __m512d resr,resi;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        const __m512d _16= _mm512_set1_pd(16.0);
                        _2r             = _mm512_add_pd(r,r);
                        _2k0a           = _mm512_add_pd(k0a,k0a);
                        const __m512d _5 = _mm512_set1_pd(5.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d c0 = _mm512_set1_pd(127.0);
                        t0              = _mm512_sub_pd(k0r,_2k0a);
                        ear             = Ir;
                        const __m512d c1 = _mm512_set1_pd(512.0);
                        eai             = t0;
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        t1              = _mm512_div_pd(a,_2r);
                        t2              = _mm512_sqrt_pd(t1);
                        facr            = _mm512_mul_pd(Er,t2);
                        faci            = _mm512_mul_pd(Ei,t2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(_16,k0a);
                        t2              = _mm512_mul_pd(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_pd(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_pd(_1,t0i);
                        t3              = _mm512_div_pd(c0,t2);
                        t0r             = _mm512_add_pd(t3,t0r);
                        t0i             = _mm512_add_pd(t3,t0i);
                        cmul_zmm8r8(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm8r8(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_storeu_pd(&EOr[0],resr);
                        _mm512_storeu_pd(&EOi[0],resi);
                 }


                     /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Optical wave component h-field, formula 4.1-35
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm8r8(const __m512d Hr,
                                         const __m512d Hi,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0a,
                                         __m512d * __restrict HOr,
                                         __m512d * __restrict HOi) {

                        register __m512d _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512d ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        const __m512d _16= _mm512_set1_pd(16.0);
                        _2r             = _mm512_add_pd(r,r);
                        _2k0a           = _mm512_add_pd(k0a,k0a);
                        const __m512d _11 = _mm512_set1_pd(11.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d c0 = _mm512_set1_pd(353.0);
                        t0              = _mm512_sub_pd(k0r,_2k0a);
                        ear             = Ir;
                        const __m512d c1 = _mm512_set1_pd(512.0);
                        eai             = t0;
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        t1              = _mm512_div_pd(a,_2r);
                        t2              = _mm512_sqrt_pd(t1);
                        facr            = _mm512_mul_pd(Er,t2);
                        faci            = _mm512_mul_pd(Ei,t2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(_16,k0a);
                        t2              = _mm512_mul_pd(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_pd(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_pd(_1,t0i);
                        t3              = _mm512_div_pd(c0,t2);
                        t0r             = _mm512_add_pd(t3,t0r);
                        t0i             = _mm512_add_pd(t3,t0i);
                        cmul_zmm8r8(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm8r8(cer1,cei1,t0r,t0i,*HOr,*HOi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pHr,
                                           const double * __restrict __ATTR_ALIGN__(64) pHi,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           double * __restrict __ATTR_ALIGN__(64) HOr,
                                           double * __restrict __ATTR_ALIGN__(64) HOi) {

                        register __m512d Hr = _mm512_load_pd(&pHr[0]);
                        register __m512d Hi = _mm512_load_pd(&pHi[0]);
                        register __m512d a  = _mm512_load_pd(&pa[0]);
                        register __m512d r  = _mm512_load_pd(&pr[0]);
                        register __m512d k0a= _mm512_load_pd(&pk0a[0]);
                        register __m512d _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i,resr,resi;
                        register __m512d ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        const __m512d _16= _mm512_set1_pd(16.0);
                        _2r             = _mm512_add_pd(r,r);
                        _2k0a           = _mm512_add_pd(k0a,k0a);
                        const __m512d _11 = _mm512_set1_pd(11.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d c0 = _mm512_set1_pd(353.0);
                        t0              = _mm512_sub_pd(k0r,_2k0a);
                        ear             = Ir;
                        const __m512d c1 = _mm512_set1_pd(512.0);
                        eai             = t0;
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        t1              = _mm512_div_pd(a,_2r);
                        t2              = _mm512_sqrt_pd(t1);
                        facr            = _mm512_mul_pd(Er,t2);
                        faci            = _mm512_mul_pd(Ei,t2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(_16,k0a);
                        t2              = _mm512_mul_pd(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_pd(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_pd(_1,t0i);
                        t3              = _mm512_div_pd(c0,t2);
                        t0r             = _mm512_add_pd(t3,t0r);
                        t0i             = _mm512_add_pd(t3,t0i);
                        cmul_zmm8r8(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm8r8(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_store_pd(&HOr[0], resr);
                        _mm512_store_pd(&HOi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm8r8_u(const double * __restrict  pHr,
                                           const double * __restrict  pHi,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0a,
                                           double * __restrict  HOr,
                                           double * __restrict  HOi) {

                        register __m512d Hr = _mm512_loadu_pd(&pHr[0]);
                        register __m512d Hi = _mm512_loadu_pd(&pHi[0]);
                        register __m512d a  = _mm512_loadu_pd(&pa[0]);
                        register __m512d r  = _mm512_loadu_pd(&pr[0]);
                        register __m512d k0a= _mm512_loadu_pd(&pk0a[0]);
                        register __m512d _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i,resr,resi;
                        register __m512d ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        const __m512d _16= _mm512_set1_pd(16.0);
                        _2r             = _mm512_add_pd(r,r);
                        _2k0a           = _mm512_add_pd(k0a,k0a);
                        const __m512d _11 = _mm512_set1_pd(11.0f;
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d c0 = _mm512_set1_pd(353.0);
                        t0              = _mm512_sub_pd(k0r,_2k0a);
                        ear             = Ir;
                        const __m512d c1 = _mm512_set1_pd(512.0);
                        eai             = t0;
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        t1              = _mm512_div_pd(a,_2r);
                        t2              = _mm512_sqrt_pd(t1);
                        facr            = _mm512_mul_pd(Er,t2);
                        faci            = _mm512_mul_pd(Ei,t2);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(_16,k0a);
                        t2              = _mm512_mul_pd(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_pd(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_pd(_1,t0i);
                        t3              = _mm512_div_pd(c0,t2);
                        t0r             = _mm512_add_pd(t3,t0r);
                        t0i             = _mm512_add_pd(t3,t0i);
                        cmul_zmm8r8(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm8r8(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_storeu_pd(&HOr[0], resr);
                        _mm512_storeu_pd(&HOi[0], resi);
                 }


                   /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Creeping wave component e-field, formula 4.1-34
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm8r8(const __m512d Er,
                                         const __m512d Ei,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0,
                                         __m512d * __restrict ECr,
                                         __m512d * __restrict ECi) {

                        register __m512d k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512d fracr,fraci,_2r,t0,t1,t2;
                        register __m512d e1ar,e1ai,exar;
                        register __m512d ce1r,ce1i,rex,t0r,t0i;
                        const __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        k0r   = _mm512_mul_pd(k0,r);
                        const __m512d c0   = _mm512_set1_pd(2.939945);
                        k0a   = _mm512_mul_pd(k0,a);
                        const __m512d c1   = _mm512_set1_pd(0.180318); 
                        k0a13 = _mm512_pow_pd(k0a,
                                         _mm512_set1_pd(0.333333333333333333333333333333333333));
                        const __m512d c2   = _mm512_set1_pd(1.821442);
                        _2r   = _mm512_add_pd(r,r);
                        const __m512d c3   = _mm512_set1_pd(-5.048945);
                        k0an13= _mm512_rcp14_pd(k0a13);
                        const __m512d c4   = _mm512_set1_pd(0.312320);
                        t0    = _mm512_div_pd(a,_2r);
                        t1    = _mm512_pow_pd(k0a,
                                          _mm512_set1_pd(0.166666666666666666666666666667));
                        k0an16= _mm512_rcp14_pd(t1);
                        t2    = _mm512_sqrt_pd(t0);
                        fracr = _mm512_mul_pd(Er,t2);
                        fraci = _mm512_mul_pd(Ei,t2);
                        t0    = _mm512_fmsub_pd(c0,k0a13,
                                            _mm512_mul_pd(c1,k0an13));
                        t1    = _mm512_fmadd_pd(k0a,PI,_mm512_add_pd(pi12,t0));
                        t1    = _mm512_add_pd(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_pd(c3,k0a13,
                                            _mm512_mul_pd(c4,k0an13));
                        t1    = _mm512_mul_pd(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_pd(t2);
                        rex   = _mm512_mul_pd(rex,t1);
                        cmul_zmm8r8(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        *ECr = _mm512_mul_pd(t0r,rex);
                        *ECi = _mm512_mul_pd(t0i,rex);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pEr,
                                         const double * __restrict __ATTR_ALIGN__(64) pEi,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         const double * __restrict __ATTR_ALIGN__(64) pr,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         double * __restrict __ATTR_ALIGN__(64) ECr,
                                         double * __restrict __ATTR_ALIGN__(64) ECi) {

                        const register __m512d Er = _mm512_load_pd(&pEr[0]);
                        const register __m512d Ei = _mm512_load_pd(&pEi[0]);
                        const register __m512d a  = _mm512_load_pd(&pa[0]);
                        const register __m512d r  = _mm512_load_pd(&pr[0]);
                        const register __m512d k0 = _mm512_load_pd(&pk0[0]);
                        register __m512d k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512d fracr,fraci,_2r,t0,t1,t2;
                        register __m512d e1ar,e1ai,exar;
                        register __m512d ce1r,ce1i,rex,t0r,t0i;
                        const __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        k0r   = _mm512_mul_pd(k0,r);
                        const __m512d c0   = _mm512_set1_pd(2.939945);
                        k0a   = _mm512_mul_pd(k0,a);
                        const __m512d c1   = _mm512_set1_pd(0.180318); 
                        k0a13 = _mm512_pow_pd(k0a,
                                         _mm512_set1_pd(0.333333333333333333333333333333333333));
                        const __m512d c2   = _mm512_set1_pd(1.821442);
                        _2r   = _mm512_add_pd(r,r);
                        const __m512d c3   = _mm512_set1_pd(-5.048945);
                        k0an13= _mm512_rcp14_pd(k0a13);
                        const __m512d c4   = _mm512_set1_pd(0.312320);
                        t0    = _mm512_div_pd(a,_2r);
                        t1    = _mm512_pow_pd(k0a,
                                          _mm512_set1_pd(0.166666666666666666666666666667));
                        k0an16= _mm512_rcp14_pd(t1);
                        t2    = _mm512_sqrt_pd(t0);
                        fracr = _mm512_mul_pd(Er,t2);
                        fraci = _mm512_mul_pd(Ei,t2);
                        t0    = _mm512_fmsub_pd(c0,k0a13,
                                            _mm512_mul_pd(c1,k0an13));
                        t1    = _mm512_fmadd_pd(k0a,PI,_mm512_add_pd(pi12,t0));
                        t1    = _mm512_add_pd(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_pd(c3,k0a13,
                                            _mm512_mul_pd(c4,k0an13));
                        t1    = _mm512_mul_pd(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_pd(t2);
                        rex   = _mm512_mul_pd(rex,t1);
                        cmul_zmm8r8(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_store_pd(&ECr[0] ,_mm512_mul_pd(t0r,rex));
                        _mm512_store_pd(&ECi[0] ,_mm512_mul_pd(t0i,rex));
                }


                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm8r8_u(const double * __restrict  pEr,
                                           const double * __restrict  pEi,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0,
                                           double * __restrict  ECr,
                                           double * __restrict  ECi) {

                        const register __m512d Er = _mm512_loadu_pd(&pEr[0]);
                        const register __m512d Ei = _mm512_loadu_pd(&pEi[0]);
                        const register __m512d a  = _mm512_loadu_pd(&pa[0]);
                        const register __m512d r  = _mm512_loadu_pd(&pr[0]);
                        const register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                        register __m512d k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512d fracr,fraci,_2r,t0,t1,t2;
                        register __m512d e1ar,e1ai,exar;
                        register __m512d ce1r,ce1i,rex,t0r,t0i;
                        const __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        k0r   = _mm512_mul_pd(k0,r);
                        const __m512d c0   = _mm512_set1_pd(2.939945);
                        k0a   = _mm512_mul_pd(k0,a);
                        const __m512d c1   = _mm512_set1_pd(0.180318); 
                        k0a13 = _mm512_pow_pd(k0a,
                                         _mm512_set1_pd(0.333333333333333333333333333333333333));
                        const __m512d c2   = _mm512_set1_pd(1.821442);
                        _2r   = _mm512_add_pd(r,r);
                        const __m512d c3   = _mm512_set1_pd(-5.048945);
                        k0an13= _mm512_rcp14_pd(k0a13);
                        const __m512d c4   = _mm512_set1_pd(0.312320);
                        t0    = _mm512_div_pd(a,_2r);
                        t1    = _mm512_pow_pd(k0a,
                                          _mm512_set1_pd(0.166666666666666666666666666667));
                        k0an16= _mm512_rcp14_pd(t1);
                        t2    = _mm512_sqrt_pd(t0);
                        fracr = _mm512_mul_pd(Er,t2);
                        fraci = _mm512_mul_pd(Ei,t2);
                        t0    = _mm512_fmsub_pd(c0,k0a13,
                                            _mm512_mul_pd(c1,k0an13));
                        t1    = _mm512_fmadd_pd(k0a,PI,_mm512_add_pd(pi12,t0));
                        t1    = _mm512_add_pd(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_pd(c3,k0a13,
                                            _mm512_mul_pd(c4,k0an13));
                        t1    = _mm512_mul_pd(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_pd(t2);
                        rex   = _mm512_mul_pd(rex,t1);
                        cmul_zmm8r8(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_storeu_pd(&ECr[0] ,_mm512_mul_pd(t0r,rex));
                        _mm512_storeu_pd(&ECi[0] ,_mm512_mul_pd(t0i,rex));
                }



                   /*

                       Backscattering creeping-wave approximation for resonance region
                       (phi == 0, k0a > 2).
                       Creeping wave component h-field, formula 4.1-36
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm8r8(const __m512d Hr,
                                         const __m512d Hi,
                                         const __m512d a,
                                         const __m512d r,
                                         const __m512d k0,
                                         __m512d * __restrict HCr,
                                         __m512d * __restrict HCi) {

                        register __m512d k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512d fracr,fraci,_2r,t0,t1,t2;
                        register __m512d e1ar,e1ai,exar;
                        register __m512d ce1r,ce1i,rex,t0r,t0i;
                        const __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        k0r   = _mm512_mul_pd(k0,r);
                        const __m512d c0   = _mm512_set1_pd(1.2701695);
                        k0a   = _mm512_mul_pd(k0,a);
                        const __m512d c1   = _mm512_set1_pd(0.2284945); 
                        k0a13 = _mm512_pow_pd(k0a,
                                         _mm512_set1_pd(0.333333333333333333333333333333333333));
                        const __m512d c2   = _mm512_set1_pd(3.063830);
                        _2r   = _mm512_add_pd(r,r);
                        const __m512d c3   = _mm512_set1_pd(-2.200000);
                        k0an13= _mm512_rcp14_pd(k0a13);
                        const __m512d c4   = _mm512_set1_pd(0.3957635);
                        t0    = _mm512_div_pd(a,_2r);
                        t1    = _mm512_pow_pd(k0a,
                                          _mm512_set1_pd(0.166666666666666666666666666667));
                        k0an16= _mm512_rcp14_pd(t1);
                        t2    = _mm512_sqrt_pd(t0);
                        fracr = _mm512_mul_pd(Hr,t2);
                        fraci = _mm512_mul_pd(Hi,t2);
                        t0    = _mm512_fmsub_pd(c0,k0a13,
                                            _mm512_mul_pd(c1,k0an13));
                        t1    = _mm512_fmadd_pd(k0a,PI,_mm512_add_pd(pi12,t0));
                        t1    = _mm512_add_pd(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_pd(c3,k0a13,
                                            _mm512_mul_pd(c4,k0an13));
                        t1    = _mm512_mul_pd(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_pd(t2);
                        rex   = _mm512_mul_pd(rex,t1);
                        cmul_zmm8r8(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        *HCr = _mm512_mul_pd(t0r,rex);
                        *HCi = _mm512_mul_pd(t0i,rex);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pHr,
                                           const double * __restrict __ATTR_ALIGN__(64) pHi,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           double * __restrict __ATTR_ALIGN__(64) HCr,
                                           double * __restrict __ATTR_ALIGN__(64) HCi) {

                        const register __m512d Hr = _mm512_load_pd(&pHr[0]);
                        const register __m512d Hi = _mm512_load_pd(&pHi[0]);
                        const register __m512d a  = _mm512_load_pd(&pa[0]);
                        const register __m512d r  = _mm512_load_pd(&pr[0]);
                        const register __m512d k0 = _mm512_load_pd(&pk0[0]);
                        register __m512d k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512d fracr,fraci,_2r,t0,t1,t2;
                        register __m512d e1ar,e1ai,exar;
                        register __m512d ce1r,ce1i,rex,t0r,t0i;
                        const __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        k0r   = _mm512_mul_pd(k0,r);
                        const __m512d c0   = _mm512_set1_pd(1.2701695);
                        k0a   = _mm512_mul_pd(k0,a);
                        const __m512d c1   = _mm512_set1_pd(0.2284945); 
                        k0a13 = _mm512_pow_pd(k0a,
                                         _mm512_set1_pd(0.333333333333333333333333333333333333));
                        const __m512d c2   = _mm512_set1_pd(3.063830);
                        _2r   = _mm512_add_pd(r,r);
                        const __m512d c3   = _mm512_set1_pd(-2.200000);
                        k0an13= _mm512_rcp14_pd(k0a13);
                        const __m512d c4   = _mm512_set1_pd(0.3957635);
                        t0    = _mm512_div_pd(a,_2r);
                        t1    = _mm512_pow_pd(k0a,
                                          _mm512_set1_pd(0.166666666666666666666666666667));
                        k0an16= _mm512_rcp14_pd(t1);
                        t2    = _mm512_sqrt_pd(t0);
                        fracr = _mm512_mul_pd(Hr,t2);
                        fraci = _mm512_mul_pd(Hi,t2);
                        t0    = _mm512_fmsub_pd(c0,k0a13,
                                            _mm512_mul_pd(c1,k0an13));
                        t1    = _mm512_fmadd_pd(k0a,PI,_mm512_add_pd(pi12,t0));
                        t1    = _mm512_add_pd(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_pd(c3,k0a13,
                                            _mm512_mul_pd(c4,k0an13));
                        t1    = _mm512_mul_pd(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_pd(t2);
                        rex   = _mm512_mul_pd(rex,t1);
                        cmul_zmm8r8(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_store_pd(&HCr[0] ,_mm512_mul_pd(t0r,rex));
                        _mm512_store_pd(&HCi[0] ,_mm512_mul_pd(t0i,rex));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm8r8_u(const double * __restrict  pHr,
                                           const double * __restrict  pHi,
                                           const double * __restrict  pa,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0,
                                           double * __restrict  HCr,
                                           double * __restrict  HCi) {

                        const register __m512d Hr = _mm512_loadu_pd(&pHr[0]);
                        const register __m512d Hi = _mm512_loadu_pd(&pHi[0]);
                        const register __m512d a  = _mm512_loadu_pd(&pa[0]);
                        const register __m512d r  = _mm512_loadu_pd(&pr[0]);
                        const register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                        register __m512d k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512d fracr,fraci,_2r,t0,t1,t2;
                        register __m512d e1ar,e1ai,exar;
                        register __m512d ce1r,ce1i,rex,t0r,t0i;
                        const __m512d pi12 = _mm512_set1_pd(0.261799387799149436538553615273);
                        k0r   = _mm512_mul_pd(k0,r);
                        const __m512d c0   = _mm512_set1_pd(1.2701695);
                        k0a   = _mm512_mul_pd(k0,a);
                        const __m512d c1   = _mm512_set1_pd(0.2284945f; 
                        k0a13 = _mm512_pow_pd(k0a,
                                         _mm512_set1_pd(0.333333333333333333333333333333333333));
                        const __m512d c2   = _mm512_set1_pd(3.063830);
                        _2r   = _mm512_add_pd(r,r);
                        const __m512d c3   = _mm512_set1_pd(-2.200000);
                        k0an13= _mm512_rcp14_pd(k0a13);
                        const __m512d c4   = _mm512_set1_pd(0.3957635);
                        t0    = _mm512_div_pd(a,_2r);
                        t1    = _mm512_pow_pd(k0a,
                                          _mm512_set1_pd(0.166666666666666666666666666667));
                        k0an16= _mm512_rcp14_pd(t1);
                        t2    = _mm512_sqrt_pd(t0);
                        fracr = _mm512_mul_pd(Hr,t2);
                        fraci = _mm512_mul_pd(Hi,t2);
                        t0    = _mm512_fmsub_pd(c0,k0a13,
                                            _mm512_mul_pd(c1,k0an13));
                        t1    = _mm512_fmadd_pd(k0a,PI,_mm512_add_pd(pi12,t0));
                        t1    = _mm512_add_pd(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm8r8(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_pd(c3,k0a13,
                                            _mm512_mul_pd(c4,k0an13));
                        t1    = _mm512_mul_pd(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_pd(t2);
                        rex   = _mm512_mul_pd(rex,t1);
                        cmul_zmm8r8(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_storeu_pd(&HCr[0] ,_mm512_mul_pd(t0r,rex));
                        _mm512_storeu_pd(&HCi[0] ,_mm512_mul_pd(t0i,rex));
                }


                  /*
                        Bistatic scattering width in high frequency limit (k0a > 20)
                        for |PI-phi| > k0a^0.3
                        Formula 4.1-37
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4137_zmm8r8(const __m512d a,
                                            const __m512d phi2) {

                          register __m512d rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_pd(PI,_mm512_mul_pd(a,cosp2));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4137_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi2) {

                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d phi2 = _mm512_load_pd(&pphi2[0]);
                          register __m512d rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_pd(PI,_mm512_mul_pd(a,cosp2));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4137_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pphi2) {

                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d phi2 = _mm512_loadu_pd(&pphi2[0]);
                          register __m512d rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_pd(PI,_mm512_mul_pd(a,cosp2));
                          return (rcs);
                 }


                    /*
                          Backscattering Width in High-Frequency Limit (k0a > 20)
                          Formula 4.1-38
                     */

                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4138_zmm8r8(const __m512d a) {

                          return (__m512d_mul_pd(a,PI));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4138_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512d a = _mm512_load_pd(&pa[0]);
                          return (__m512d_mul_pd(a,PI));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4138_zmm8r8_u(const double * __restrict  pa) {

                          register __m512d a = _mm512_loadu_pd(&pa[0]);
                          return (__m512d_mul_pd(a,PI));
                  }


                   /*
                         Forward scattering widths and pattern in high-frequency limit
                         (k0a>20.0)
                         Formula 4.1-40, RCS.
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4140_zmm8r8(const __m512d k0a,
                                            const __m512d alpha) {

                          register __m512d sinc,k0alp,k0as,t0;
                          register __m512d rcs;
                          const __m512d _4       = _mm512_set1_pd(4.0);
                          k0alp = _mm512_mul_pd(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_pd(t0,k0alp); 
                          k0as  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a)); 
                          rcs   = _mm512_mul_pd(k0as,_mm512_mul_pd(sinc,sinc));
                          return (rcs);        
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4140_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) palpha) {

                          register __m512d k0a   = _mm512_load_pd(&pk0a[0]);
                          register __m512d alpha = _mm512_load_pd(&palpha[0]);
                          register __m512d sinc,k0alp,k0as,t0;
                          register __m512d rcs;
                          const __m512d _4       = _mm512_set1_pd(4.0);
                          k0alp = _mm512_mul_pd(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_pd(t0,k0alp); 
                          k0as  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a)); 
                          rcs   = _mm512_mul_pd(k0as,_mm512_mul_pd(sinc,sinc));
                          return (rcs);        
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4140_zmm8r8_u(const double * __restrict  pk0a,
                                              const double * __restrict  palpha) {

                          register __m512d k0a   = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d alpha = _mm512_loadu_pd(&palpha[0]);
                          register __m512d sinc,k0alp,k0as,t0;
                          register __m512d rcs;
                          const __m512d _4       = _mm512_set1_pd(4.0);
                          k0alp = _mm512_mul_pd(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_pd(t0,k0alp); 
                          k0as  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a)); 
                          rcs   = _mm512_mul_pd(k0as,_mm512_mul_pd(sinc,sinc));
                          return (rcs);        
                  }


                     /*
                         Forward scattering widths and pattern in high-frequency limit
                         (k0a>20.0), forward scattered (diffracted) e-field
                         Formula 4.1-39.

                       */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm8r8(const __m512d Er,
                                         const __m512d Ei,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d alp
                                         const __m512d k0a,
                                         __m512d * __restrict Esr,
                                         __m512d * __restrict Esi) {

                        register __m512d _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512d facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512d cer,cei;
                        const __m512d _2  = _mm512_set1_pd(2.0);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        _2k0as = _mm512_add_pd(_2,_mm512_mul_pd(k0a,k0a));
                        k0r    = _mm512_mul_pd(k0,r);
                        k0alp  = _mm512_mul_pd(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_pd(k0r,pir);
                        sinc   = _mm512_div_pd(t0,k0alp);
                        cexp_zmm8r8(arr,ari,&cer,&cei);
                        div    = _mm512_div_pd(_2k0as,pi4);
                        t0     = _mm512_sqrt_pd(div);
                        facr   = _mm512_mul_pd(Er,t0);
                        t0r    = _mm512_mul_pd(cer,sinc);
                        faci   = _mm512_mul_pd(Ei,t0);
                        t0i    = _mm512_mul_pd(cei,sinc);
                        cmul_zmm8r8(facr,faci,t0r,t0i,*Esr,*Esi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pEr,
                                           const double * __restrict __ATTR_ALIGN__(64) pEi,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) palp
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           double * __restrict __ATTR_ALIGN__(64) Esr,
                                           double * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512d Er = _mm512_load_pd(&pEr[0]);
                        register __m512d Ei = _mm512_load_pd(&pEi[0]);
                        register __m512d r  = _mm512_load_pd(&pr[0]);
                        register __m512d k0 = _mm512_load_pd(&pk0[0]);
                        register __m512d alp= _mm512_load_pd(&palp[0]);
                        register __m512d k0a= _mm512_load_pd(&pk0a[0]);
                        register __m512d _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512d facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512d cer,cei,resr,resi;
                        const __m512d _2  = _mm512_set1_pd(2.0);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        _2k0as = _mm512_add_pd(_2,_mm512_mul_pd(k0a,k0a));
                        k0r    = _mm512_mul_pd(k0,r);
                        k0alp  = _mm512_mul_pd(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_pd(k0r,pir);
                        sinc   = _mm512_div_pd(t0,k0alp);
                        cexp_zmm8r8(arr,ari,&cer,&cei);
                        div    = _mm512_div_pd(_2k0as,pi4);
                        t0     = _mm512_sqrt_pd(div);
                        facr   = _mm512_mul_pd(Er,t0);
                        t0r    = _mm512_mul_pd(cer,sinc);
                        faci   = _mm512_mul_pd(Ei,t0);
                        t0i    = _mm512_mul_pd(cei,sinc);
                        cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                        _mm512_store_pd(&Esr[0], resr);
                        _mm512_store_pd(&Esi[0], resi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm8r8_u(const double * __restrict  pEr,
                                           const double * __restrict  pEi,
                                           const double * __restrict  pr,
                                           const double * __restrict  pk0,
                                           const double * __restrict  palp
                                           const double * __restrict  pk0a,
                                           double * __restrict  Esr,
                                           double * __restrict  Esi) {

                        register __m512d Er = _mm512_loadu_pd(&pEr[0]);
                        register __m512d Ei = _mm512_loadu_pd(&pEi[0]);
                        register __m512d r  = _mm512_loadu_pd(&pr[0]);
                        register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
                        register __m512d alp= _mm512_loadu_pd(&palp[0]);
                        register __m512d k0a= _mm512_loadu_pd(&pk0a[0]);
                        register __m512d _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512d facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512d cer,cei,resr,resi;
                        const __m512d _2  = _mm512_set1_pd(2.0);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        _2k0as = _mm512_add_pd(_2,_mm512_mul_pd(k0a,k0a));
                        k0r    = _mm512_mul_pd(k0,r);
                        k0alp  = _mm512_mul_pd(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_pd(k0r,pir);
                        sinc   = _mm512_div_pd(t0,k0alp);
                        cexp_zmm8r8(arr,ari,&cer,&cei);
                        div    = _mm512_div_pd(_2k0as,pi4);
                        t0     = _mm512_sqrt_pd(div);
                        facr   = _mm512_mul_pd(Er,t0);
                        t0r    = _mm512_mul_pd(cer,sinc);
                        faci   = _mm512_mul_pd(Ei,t0);
                        t0i    = _mm512_mul_pd(cei,sinc);
                        cmul_zmm8r8(facr,faci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_pd(&Esr[0], resr);
                        _mm512_storeu_pd(&Esi[0], resi);
              }


                  /*
                         Forward scattering widths and pattern in high-frequency limit
                         (k0a>20.0), constant angle (alpha=0)
                         Formula 4.1-41, RCS.
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4141_zmm8r8(const __m512d k0a) {

                          const __m512d _4 = _mm512_set1_pd(4.0);
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a));
                          return (rcs);
                 }

                

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4141_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0a) {

                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const __m512d _4 = _mm512_set1_pd(4.0);
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4141_zmm8r8_u(const double * __restrict   pk0a) {

                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d _4 = _mm512_set1_pd(4.0);
                          register __m512d rcs;
                          rcs = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,k0a));
                          return (rcs);
                 }


                   /*
                        Approximations for the low frequency region (k0a<<1,k1a<<1)
                        Scattered far-zone e-field, formula 4.1-45
                    */

                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm8r8(const __m512d EIr,
                                         const __m512d EIi,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d k0a,
                                         const __m512d phi,
                                         const __m512d eps0,
                                         const __m512d eps1,
                                         const __m512d mu0,
                                         const __m512d mu1,
                                         __m512d * __restrict ESr,
                                         __m512d * __restrict ESi) {

                        register __m512d k0r,k0as,fracr,fraci,k0as2;
                        register __m512d ear,eai,cer,cei,t0r,t0i;
                        register __m512d t0,t1,cosp,t2,sk0r,t3,mul;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d _2 = _mm512_set1_pd(2.0);
                        sk0r            = _mm512_sqrt_pd(k0r);
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        const __m512d hlf= _mm512_set1_pd(0.5);
                        k0as2           = _mm512_mul_pd(hlf,k0as);
                        const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                        cosp            = xcosf(phi);
                        const __m512d pi2= _mm512_set1_pd(1.253314137315500251207882642406);
                        fracr           = _mm512_mul_pd(EIr,pi2);
                        t0              = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                        fraci           = _mm512_mul_pd(EIi,pi2);
                        t1              = _mm512_div_pd(_mm512_sub_pd(mu1,mu0),
                                                        _mm512_add_pd(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_pd(t1,t1);
                        eai             = _mm512_sub_pd(k0r,pi4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(t2,cosp);
                        cer = _mm512_div_pd(cer,sk0r);
                        t3  = _mm512_sub_pd(t0,t1);
                        cei = _mm512_div_pd(cei,sk0r);
                        mul = _mm512_mul_pd(k0as2,t3);
                        t0r = _mm512_mul_pd(cer,mul);
                        t0i = _mm512_mul_pd(cei,mul);
                        cmul_zmm8r8(fracr,fraci,t0r,t0i,*ESr,*ESi);
               }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pEIr,
                                           const double * __restrict __ATTR_ALIGN__(64) pEIi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const double * __restrict __ATTR_ALIGN__(64)  peps0,
                                           const double * __restrict __ATTR_ALIGN__(64)  peps1,
                                           const double * __restrict __ATTR_ALIGN__(64)  pmu0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pmu1,
                                           double * __restrict __ATTR_ALIGN__(64)  ESr,
                                           double * __restrict __ATTR_ALIGN__(64)  ESi) {

                        register __m512d EIr = _mm512_load_pd(&pEIr[0]);
                        register __m512d EIi = _mm512_load_pd(&pEIi[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d phi = _mm512_load_pd(&pphi[0]);
                        register __m512d eps0= _mm512_load_pd(&peps0[0]);
                        register __m512d eps1= _mm512_load_pd(&peps1[0]);
                        register __m512d mu0 = _mm512_load_pd(&pmu0[0]);
                        register __m512d mu1 = _mm512_load_pd(&pmu1[0]);
                        register __m512d k0r,k0as,fracr,fraci,k0as2;
                        register __m512d ear,eai,cer,cei,t0r,t0i;
                        register __m512d t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d _2 = _mm512_set1_pd(2.0);
                        sk0r            = _mm512_sqrt_pd(k0r);
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        const __m512d hlf= _mm512_set1_pd(0.5);
                        k0as2           = _mm512_mul_pd(hlf,k0as);
                        const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                        cosp            = xcosf(phi);
                        const __m512d pi2= _mm512_set1_pd(1.253314137315500251207882642406);
                        fracr           = _mm512_mul_pd(EIr,pi2);
                        t0              = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                        fraci           = _mm512_mul_pd(EIi,pi2);
                        t1              = _mm512_div_pd(_mm512_sub_pd(mu1,mu0),
                                                        _mm512_add_pd(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_pd(t1,t1);
                        eai             = _mm512_sub_pd(k0r,pi4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(t2,cosp);
                        cer = _mm512_div_pd(cer,sk0r);
                        t3  = _mm512_sub_pd(t0,t1);
                        cei = _mm512_div_pd(cei,sk0r);
                        mul = _mm512_mul_pd(k0as2,t3);
                        t0r = _mm512_mul_pd(cer,mul);
                        t0i = _mm512_mul_pd(cei,mul);
                        cmul_zmm8r8(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_store_pd(&ESr[0], resr);
                        _mm512_store_pd(&ESi[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm8r8_u(const double * __restrict  pEIr,
                                           const double * __restrict  pEIi,
                                           const double * __restrict  pr,
                                           const double * __restrict   pk0,
                                           const double * __restrict  pk0a,
                                           const double * __restrict   pphi,
                                           const double * __restrict   peps0,
                                           const double * __restrict  peps1,
                                           const double * __restrict   pmu0,
                                           const double * __restrict   pmu1,
                                           double * __restrict   ESr,
                                           double * __restrict   ESi) {

                        register __m512d EIr = _mm512_loadu_pd(&pEIr[0]);
                        register __m512d EIi = _mm512_loadu_pd(&pEIi[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                        register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                        register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                        register __m512d mu0 = _mm512_loadu_pd(&pmu0[0]);
                        register __m512d mu1 = _mm512_loadu_pd(&pmu1[0]);
                        register __m512d k0r,k0as,fracr,fraci,k0as2;
                        register __m512d ear,eai,cer,cei,t0r,t0i;
                        register __m512d t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d _2 = _mm512_set1_pd(2.0);
                        sk0r            = _mm512_sqrt_pd(k0r);
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        const __m512d hlf= _mm512_set1_pd(0.5);
                        k0as2           = _mm512_mul_pd(hlf,k0as);
                        const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                        cosp            = xcosf(phi);
                        const __m512d pi2= _mm512_set1_pd(1.253314137315500251207882642406);
                        fracr           = _mm512_mul_pd(EIr,pi2);
                        t0              = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                        fraci           = _mm512_mul_pd(EIi,pi2);
                        t1              = _mm512_div_pd(_mm512_sub_pd(mu1,mu0),
                                                        _mm512_add_pd(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_pd(t1,t1);
                        eai             = _mm512_sub_pd(k0r,pi4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(t2,cosp);
                        cer = _mm512_div_pd(cer,sk0r);
                        t3  = _mm512_sub_pd(t0,t1);
                        cei = _mm512_div_pd(cei,sk0r);
                        mul = _mm512_mul_pd(k0as2,t3);
                        t0r = _mm512_mul_pd(cer,mul);
                        t0i = _mm512_mul_pd(cei,mul);
                        cmul_zmm8r8(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_pd(&ESr[0], resr);
                        _mm512_storeu_pd(&ESi[0], resi);
               }


                  /*
                        Approximations for the low frequency region (k0a<<1,k1a<<1)
                        Scattered far-zone h-field, formula 4.1-46
                    */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm8r8(const __m512d HIr,
                                         const __m512d HIi,
                                         const __m512d r,
                                         const __m512d k0,
                                         const __m512d k0a,
                                         const __m512d phi,
                                         const __m512d eps0,
                                         const __m512d eps1,
                                         const __m512d mu0,
                                         const __m512d mu1,
                                         __m512d * __restrict HSr,
                                         __m512d * __restrict HSi) {

                        register __m512d k0r,k0as,fracr,fraci,k0as2;
                        register __m512d ear,eai,cer,cei,t0r,t0i;
                        register __m512d t0,t1,cosp,t2,sk0r,t3,mul;
                        const __m512d _1 = _mm512_set1_pd(1.0);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d _2 = _mm512_set1_pd(2.0);
                        sk0r            = _mm512_sqrt_pd(k0r);
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        const __m512d hlf= _mm512_set1_pd(0.5);
                        k0as2           = _mm512_mul_pd(hlf,k0as);
                        const __m512d pi4= _mm512_set1_pd(0.785398163397448309615660845820);
                        cosp            = xcosf(phi);
                        const __m512d pi2= _mm512_set1_pd(1.253314137315500251207882642406);
                        fracr           = _mm512_mul_pd(HIr,pi2);
                        t0              = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                        fraci           = _mm512_mul_pd(HIi,pi2);
                        t1              = _mm512_div_pd(_mm512_sub_pd(eps1,eps0),
                                                        _mm512_add_pd(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_pd(t1,t1);
                        eai             = _mm512_sub_pd(k0r,pi4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(t2,cosp);
                        cer = _mm512_div_pd(cer,sk0r);
                        t3  = _mm512_sub_pd(t0,t1);
                        cei = _mm512_div_pd(cei,sk0r);
                        mul = _mm512_mul_pd(k0as2,t3);
                        t0r = _mm512_mul_pd(cer,mul);
                        t0i = _mm512_mul_pd(cei,mul);
                        cmul_zmm8r8(fracr,fraci,t0r,t0i,*HSr,*HSi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pHIr,
                                           const double * __restrict __ATTR_ALIGN__(64) pHIi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const double * __restrict __ATTR_ALIGN__(64)  peps0,
                                           const double * __restrict __ATTR_ALIGN__(64)  peps1,
                                           const double * __restrict __ATTR_ALIGN__(64)  pmu0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pmu1,
                                           double * __restrict __ATTR_ALIGN__(64)  HSr,
                                           double * __restrict __ATTR_ALIGN__(64)  HSi) {

                        register __m512d HIr = _mm512_load_pd(&pHIr[0]);
                        register __m512d HIi = _mm512_load_pd(&pHIi[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d phi = _mm512_load_pd(&pphi[0]);
                        register __m512d eps0= _mm512_load_pd(&peps0[0]);
                        register __m512d eps1= _mm512_load_pd(&peps1[0]);
                        register __m512d mu0 = _mm512_load_pd(&pmu0[0]);
                        register __m512d mu1 = _mm512_load_pd(&pmu1[0]);
                        register __m512d k0r,k0as,fracr,fraci,k0as2;
                        register __m512d ear,eai,cer,cei,t0r,t0i;
                        register __m512d t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d _2 = _mm512_set1_pd(2.0f);
                        sk0r            = _mm512_sqrt_pd(k0r);
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        const __m512d hlf= _mm512_set1_pd(0.5f);
                        k0as2           = _mm512_mul_pd(hlf,k0as);
                        const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                        cosp            = xcosf(phi);
                        const __m512d pi2= _mm512_set1_pd(1.253314137315500251207882642406);
                        fracr           = _mm512_mul_pd(HIr,pi2);
                        t0              = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                        fraci           = _mm512_mul_pd(HIi,pi2);
                        t1              = _mm512_div_pd(_mm512_sub_pd(eps1,eps0),
                                                        _mm512_add_pd(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_pd(t1,t1);
                        eai             = _mm512_sub_pd(k0r,pi4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(t2,cosp);
                        cer = _mm512_div_pd(cer,sk0r);
                        t3  = _mm512_sub_pd(t0,t1);
                        cei = _mm512_div_pd(cei,sk0r);
                        mul = _mm512_mul_pd(k0as2,t3);
                        t0r = _mm512_mul_pd(cer,mul);
                        t0i = _mm512_mul_pd(cei,mul);
                        cmul_zmm8r8(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_store_pd(&HSr[0], resr);
                        _mm512_store_pd(&HSi[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm8r8_u(const double * __restrict   pHIr,
                                           const double * __restrict   pHIi,
                                           const double * __restrict   pr,
                                           const double * __restrict   pk0,
                                           const double * __restrict   pk0a,
                                           const double * __restrict   pphi,
                                           const double * __restrict   peps0,
                                           const double * __restrict   peps1,
                                           const double * __restrict   pmu0,
                                           const double * __restrict   pmu1,
                                           double * __restrict  HSr,
                                           double * __restrict   HSi) {

                        register __m512d HIr = _mm512_loadu_pd(&pHIr[0]);
                        register __m512d HIi = _mm512_loadu_pd(&pHIi[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                        register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                        register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                        register __m512d mu0 = _mm512_loadu_pd(&pmu0[0]);
                        register __m512d mu1 = _mm512_loadu_pd(&pmu1[0]);
                        register __m512d k0r,k0as,fracr,fraci,k0as2;
                        register __m512d ear,eai,cer,cei,t0r,t0i;
                        register __m512d t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        k0r             = _mm512_mul_pd(k0,r);
                        const __m512d _2 = _mm512_set1_pd(2.0f);
                        sk0r            = _mm512_sqrt_pd(k0r);
                        k0as            = _mm512_mul_pd(k0a,k0a);
                        const __m512d hlf= _mm512_set1_pd(0.5f);
                        k0as2           = _mm512_mul_pd(hlf,k0as);
                        const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                        cosp            = xcosf(phi);
                        const __m512d pi2= _mm512_set1_pd(1.253314137315500251207882642406);
                        fracr           = _mm512_mul_pd(HIr,pi2);
                        t0              = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                        fraci           = _mm512_mul_pd(HIi,pi2);
                        t1              = _mm512_div_pd(_mm512_sub_pd(eps1,eps0),
                                                        _mm512_add_pd(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_pd(t1,t1);
                        eai             = _mm512_sub_pd(k0r,pi4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_pd(t2,cosp);
                        cer = _mm512_div_pd(cer,sk0r);
                        t3  = _mm512_sub_pd(t0,t1);
                        cei = _mm512_div_pd(cei,sk0r);
                        mul = _mm512_mul_pd(k0as2,t3);
                        t0r = _mm512_mul_pd(cer,mul);
                        t0i = _mm512_mul_pd(cei,mul);
                        cmul_zmm8r8(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_pd(&HSr[0], resr);
                        _mm512_storeu_pd(&HSi[0], resi);
               }


                 /*
                      Bistatic scattering width (k0a<<1, k1a<<1) at the angle 'phi'
                      Formula 4.1-47

                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4147_zmm8r8(const __m512d a,
                                            const __m512d k0a,
                                            const __m512d phi,
                                            const __m512d eps1,
                                            const __m512d eps0,
                                            const __m512d mu1,
                                            const __m512d mu0) {

                          register __m512d t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512d rcs;
                          cosp             = xcosf(phi);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,_mm512_mul_pd(mut,cosp));
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                }


                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4147_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) pphi,
                                            const double * __restrict __ATTR_ALIGN__(64) peps1,
                                            const double * __restrict __ATTR_ALIGN__(64) peps0,
                                            const double * __restrict __ATTR_ALIGN__(64) pmu1,
                                            const double * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512d a = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d phi = _mm512_load_pd(&pphi[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps0= _mm512_load_pd(&peps0[0]);
                          register __m512d pmu1= _mm512_load_pd(&pmu1[0]);
                          register __m512d pmu0= _mm512_load_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512d rcs;
                          cosp             = xcosf(phi);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,_mm512_mul_pd(mut,cosp));
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4147_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict  pk0a,
                                            const double * __restrict  pphi,
                                            const double * __restrict  peps1,
                                            const double * __restrict  peps0,
                                            const double * __restrict  pmu1,
                                            const double * __restrict  pmu0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                          register __m512d pmu1= _mm512_loadu_pd(&pmu1[0]);
                          register __m512d pmu0= _mm512_loadu_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512d rcs;
                          cosp             = xcosf(phi);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,_mm512_mul_pd(mut,cosp));
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                }   


                    /*
                      Bistatic scattering width (k0a<<1, k1a<<1) at the angle 'phi'
                      Formula 4.1-48

                   */   

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4148_zmm8r8(const __m512d a,
                                            const __m512d k0a,
                                            const __m512d phi,
                                            const __m512d eps1,
                                            const __m512d eps0,
                                            const __m512d mu1,
                                            const __m512d mu0) {

                          register __m512d t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512d rcs;
                          cosp             = xcosf(phi);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,_mm512_mul_pd(mut,cosp));
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                }  


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4148_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) pphi,
                                            const double * __restrict __ATTR_ALIGN__(64) peps1,
                                            const double * __restrict __ATTR_ALIGN__(64) peps0,
                                            const double * __restrict __ATTR_ALIGN__(64) pmu1,
                                            const double * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d phi = _mm512_load_pd(&pphi[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps0= _mm512_load_pd(&peps0[0]);
                          register __m512d pmu1= _mm512_load_pd(&pmu1[0]);
                          register __m512d pmu0= _mm512_load_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512d rcs;
                          cosp             = xcosf(phi);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,_mm512_mul_pd(mut,cosp));
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                }  


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4148_zmm8r8_u(const double * __restrict  pa,
                                            const double * __restrict pk0a,
                                            const double * __restrict  pphi,
                                            const double * __restrict  peps1,
                                            const double * __restrict  peps0,
                                            const double * __restrict  pmu1,
                                            const double * __restrict  pmu0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                          register __m512d pmu1= _mm512_loadu_pd(&pmu1[0]);
                          register __m512d pmu0= _mm512_loadu_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512d rcs;
                          cosp             = xcosf(phi);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,_mm512_mul_pd(mut,cosp));
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);
                }  


                   /*
                         Backscattering width (k0a<<1,k1a<<1), when phi = 0
                         Formula 4.1-49
                    */
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4149_zmm8r8(const __m512d a,
                                            const __m512d k0a,
                                            const __m512d eps1,
                                            const __m512d eps0,
                                            const __m512d mu1,
                                            const __m512d mu0) {

                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4149_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) peps1,
                                              const double * __restrict __ATTR_ALIGN__(64) peps0,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps0= _mm512_load_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_load_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_load_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4149_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  peps1,
                                              const double * __restrict  peps0,
                                              const double * __restrict  pmu1,
                                              const double * __restrict  pmu0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_loadu_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_loadu_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                     /*
                         Backscattering width (k0a<<1,k1a<<1), when phi = 0
                         Formula 4.1-50
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4150_zmm8r8(const __m512d a,
                                            const __m512d k0a,
                                            const __m512d eps1,
                                            const __m512d eps0,
                                            const __m512d mu1,
                                            const __m512d mu0) {

                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4150_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) peps1,
                                              const double * __restrict __ATTR_ALIGN__(64) peps0,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps0= _mm512_load_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_load_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_load_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4150_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  peps1,
                                              const double * __restrict  peps0,
                                              const double * __restrict  pmu1,
                                              const double * __restrict  pmu0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_loadu_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_loadu_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_sub_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                    /*
                         Forward scattering width (k0a<<1, k1a<<1), phi = pi
                         Formula 4.1-51
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4151_zmm8r8(const __m512d a,
                                            const __m512d k0a,
                                            const __m512d eps1,
                                            const __m512d eps0,
                                            const __m512d mu1,
                                            const __m512d mu0) {

                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_add_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4151_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) peps1,
                                              const double * __restrict __ATTR_ALIGN__(64) peps0,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps0= _mm512_load_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_load_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_load_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_add_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4151_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  peps1,
                                              const double * __restrict  peps0,
                                              const double * __restrict  pmu1,
                                              const double * __restrict  pmu0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_loadu_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_loadu_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(eps1,eps0),_1);
                          t1               = _mm512_sub_pd(mu1,mu0);
                          t2               = _mm512_add_pd(mu1,mu0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_add_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                      /*
                         Forward scattering width (k0a<<1, k1a<<1), phi = pi
                         Formula 4.1-52
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4152_zmm8r8(const __m512d a,
                                            const __m512d k0a,
                                            const __m512d eps1,
                                            const __m512d eps0,
                                            const __m512d mu1,
                                            const __m512d mu0) {

                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_add_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4152_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) peps1,
                                              const double * __restrict __ATTR_ALIGN__(64) peps0,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu1,
                                              const double * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_load_pd(&peps1[0]);
                          register __m512d eps0= _mm512_load_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_load_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_load_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_add_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4152_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  peps1,
                                              const double * __restrict  peps0,
                                              const double * __restrict  pmu1,
                                              const double * __restrict  pmu0) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a3[0]);
                          register __m512d eps1= _mm512_loadu_pd(&peps1[0]);
                          register __m512d eps0= _mm512_loadu_pd(&peps0[0]);
                          register __m512d mu1 = _mm512_loadu_pd(&pmu1[0]);
                          register __m512d mu0 = _mm512_loadu_pd(&pmu0[0]);
                          register __m512d t0,t1,k0a3,epst,mut,sqr,t2,diff;
                          register __m512d rcs;
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          t0               = _mm512_mul_pd(pi4,_mm512_mul_pd(PI,a));
                          const __m512d _1  = _mm512_set1_pd(1.0f);
                          k0a3             = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          const __m512d _2  = _mm512_set1_pd(2.0f);
                          epst             = _mm512_sub_pd(_mm512_div_pd(mu1,mu0),_1);
                          t1               = _mm512_sub_pd(eps1,eps0);
                          t2               = _mm512_add_pd(eps1,eps0);
                          mut              = _mm512_mul_pd(_2,_mm512_div_pd(t1,t2));
                          diff             = _mm512_add_pd(epst,mut);
                          sqr              = _mm512_mul_pd(diff,diff);
                          rcs              = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,sqr));
                          return (rcs);

                   }


                     /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-72
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4172_zmm8r8(const __m512d mur,
                                          const __m512d mui,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d psi,
                                          __m512d * __restrict Tinr,
                                          __m512d * __restrict Tini) {

                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d sin2p,cosp,divr,divi,t1;
                        register __m512d sqr1,sqi1,sqr2,sqi2,t0;
                        register __m512d mulr,muli,t0r,t0i,t1r,t1i;
                        register __m512d t2r,t2i,t3r,t3i;
                        cosp = xcosf(psi);
                        t0   = xsinf(psi);
                        sin2p= _mm512_mul_pd(t0,t0);
                        cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                        t1   = _mm512_sub_pd(_1,sin2p);
                        cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                        csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                        t0r = _mm512_div_pd(t1,mulr);
                        t0i = _mm512_div_pd(t1,muli);
                        csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                        t2r = _mm512_add_pd(sqr1,sqr1);
                        t2i = _mm512_add_pd(sqi1,sqi1);
                        cmul_zmm8r8(t2r,t2i,sqr2,sqi2,&t1r,&t1i);//numerator
                        cmul_zmm8r8(sqr1,sqi1,sqr2,sqi2,&t3r,&t3i); // denum
                        t3r = _mm512_add_pd(cosp,t3r);
                        t3i = _mm512_add_pd(cosp,t3i);
                        cdiv_zmm8r8(t1r,t1i,t3r,t3i,*Tinr,*Tini);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4172_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pmur,
                                            const double * __restrict __ATTR_ALIGN__(64) pmui,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                            double * __restrict __ATTR_ALIGN__(64) Tinr,
                                            double * __restrict __ATTR_ALIGN__(64) Tini) {

                        register __m512d mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d mui  = _mm512_load_pd(&pmui[0]);
                        register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d sin2p,cosp,divr,divi,t1;
                        register __m512d sqr1,sqi1,sqr2,sqi2,t0;
                        register __m512d mulr,muli,t0r,t0i,t1r,t1i;
                        register __m512d t2r,t2i,t3r,t3i,resr,resi;
                        cosp = xcosf(psi);
                        t0   = xsinf(psi);
                        sin2p= _mm512_mul_pd(t0,t0);
                        cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                        t1   = _mm512_sub_pd(_1,sin2p);
                        cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                        csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                        t0r = _mm512_div_pd(t1,mulr);
                        t0i = _mm512_div_pd(t1,muli);
                        csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                        t2r = _mm512_add_pd(sqr1,sqr1);
                        t2i = _mm512_add_pd(sqi1,sqi1);
                        cmul_zmm8r8(t2r,t2i,sqr2,sqi2,&t1r,&t1i);//numerator
                        cmul_zmm8r8(sqr1,sqi1,sqr2,sqi2,&t3r,&t3i); // denum
                        t3r = _mm512_add_pd(cosp,t3r);
                        t3i = _mm512_add_pd(cosp,t3i);
                        cdiv_zmm8r8(t1r,t1i,t3r,t3i,&resr,&resi);
                        _mm512_store_pd(&Tinr[0], resr);
                        _mm512_store_pd(&Tini[0], resi);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4172_zmm8r8_u(const double * __restrict  pmur,
                                            const double * __restrict  pmui,
                                            const double * __restrict  pepsr,
                                            const double * __restrict  pepsi,
                                            const double * __restrict  ppsi,
                                            double * __restrict  Tinr,
                                            double * __restrict  Tini) {

                        register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                        register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                        const __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d sin2p,cosp,divr,divi,t1;
                        register __m512d sqr1,sqi1,sqr2,sqi2,t0;
                        register __m512d mulr,muli,t0r,t0i,t1r,t1i;
                        register __m512d t2r,t2i,t3r,t3i,resr,resi;
                        cosp = xcosf(psi);
                        t0   = xsinf(psi);
                        sin2p= _mm512_mul_pd(t0,t0);
                        cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                        t1   = _mm512_sub_pd(_1,sin2p);
                        cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                        csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                        t0r = _mm512_div_pd(t1,mulr);
                        t0i = _mm512_div_pd(t1,muli);
                        csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                        t2r = _mm512_add_pd(sqr1,sqr1);
                        t2i = _mm512_add_pd(sqi1,sqi1);
                        cmul_zmm8r8(t2r,t2i,sqr2,sqi2,&t1r,&t1i);//numerator
                        cmul_zmm8r8(sqr1,sqi1,sqr2,sqi2,&t3r,&t3i); // denum
                        t3r = _mm512_add_pd(cosp,t3r);
                        t3i = _mm512_add_pd(cosp,t3i);
                        cdiv_zmm8r8(t1r,t1i,t3r,t3i,&resr,&resi);
                        _mm512_storeu_pd(&Tinr[0], resr);
                        _mm512_storeu_pd(&Tini[0], resi);
                  }


                     /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-73
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4173_zmm8r8(const __m512d mur,
                                          const __m512d mui,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d psi,
                                          __m512d * __restrict Tinr,
                                          __m512d * __restrict Tini) {

                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         register __m512d cosp,_2cosp,divr,divi;
                         register __m512d sqr1,sqi1,sqr2,sqi2;
                         register __m512d sinp,sin2p,mulr,muli;
                         register __m512d t0r,t0i,_1msp;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         _2cosp = _mm512_add_pd(cosp,cosp);
                         sin2p  = _mm512_mul_pd(sinp,sinp);
                         _1msp  = _mm512_sub_pd(_1,sin2p);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&divi);
                         t0r = _mm512_div_pd(_1msp,mulr);
                         t0i = _mm512_div_pd(_1msp,muli);
                         csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                         csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                         *Tinr = _mm512_fmadd_pd(sqr1,sqr2,cosp);
                         *Tini = _mm512_fmadd_pd(sqi1,sqi2,cosp);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4173_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pmur,
                                            const double * __restrict __ATTR_ALIGN__(64) pmui,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                            double * __restrict __ATTR_ALIGN__(64) Tinr,
                                            double * __restrict __ATTR_ALIGN__(64) Tini) {

                         register __m512d mur  = _mm512_load_pd(&pmur[0]);
                         register __m512d mui  = _mm512_load_pd(&pmui[0]);
                         register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         register __m512d cosp,_2cosp,divr,divi;
                         register __m512d sqr1,sqi1,sqr2,sqi2;
                         register __m512d sinp,sin2p,mulr,muli;
                         register __m512d t0r,t0i,_1msp;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         _2cosp = _mm512_add_pd(cosp,cosp);
                         sin2p  = _mm512_mul_pd(sinp,sinp);
                         _1msp  = _mm512_sub_pd(_1,sin2p);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&divi);
                         t0r = _mm512_div_pd(_1msp,mulr);
                         t0i = _mm512_div_pd(_1msp,muli);
                         csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                         csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                         _mm512_store_pd(&Tinr[0] ,_mm512_fmadd_pd(sqr1,sqr2,cosp));
                         _mm512_store_pd(&Tini[0] ,_mm512_fmadd_pd(sqi1,sqi2,cosp));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4173_zmm8r8_u(const double * __restrict  pmur,
                                            const double * __restrict  pmui,
                                            const double * __restrict  pepsr,
                                            const double * __restrict  pepsi,
                                            const double * __restrict  ppsi,
                                            double * __restrict  Tinr,
                                            double * __restrict  Tini) {

                         register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                         register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         register __m512d cosp,_2cosp,divr,divi;
                         register __m512d sqr1,sqi1,sqr2,sqi2;
                         register __m512d sinp,sin2p,mulr,muli;
                         register __m512d t0r,t0i,_1msp;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         _2cosp = _mm512_add_pd(cosp,cosp);
                         sin2p  = _mm512_mul_pd(sinp,sinp);
                         _1msp  = _mm512_sub_pd(_1,sin2p);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&divi);
                         t0r = _mm512_div_pd(_1msp,mulr);
                         t0i = _mm512_div_pd(_1msp,muli);
                         csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                         csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                         _mm512_storeu_pd(&Tinr[0] ,_mm512_fmadd_pd(sqr1,sqr2,cosp));
                         _mm512_storeu_pd(&Tini[0] ,_mm512_fmadd_pd(sqi1,sqi2,cosp));
                 }


                    /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-74
                       */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4174_zmm8r8(const __m512d mur,
                                          const __m512d mui,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d psi,
                                          __m512d * __restrict Toutr,
                                          __m512d * __restrict Touti) {

                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d divr,divi,sqr1,sqi1;
                         register __m512d mulr,muli,sqr2,sqi2;
                         register __m512d cosp,sinp,sin2p,t0r,t0i;
                         register __m512d _2sqr1,_2sqi1,t1r,t1i;
                         register __m512d numr,numi,denr,deni;
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&div);
                         cosp = xcosf(psi);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mulr,&muli);
                         sinp = xsinf(psi);
                         csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                         sin2p= _mm512_mul_pd(sinp,sinp);
                         //_2sqr1 = _mm512_mul_pd(_2,sqr1);
                         t0r    = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         //_2sqi1 = _mm512_mul_pd(_2,sqi1);
                         t0i    = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                         cmul_zmm8r8(sqr1,sqi1,t0r,t0i,&t1r,&t1i);
                         numr = _mm512_mul_pd(_2,t1r);
                         denr = _mm512_add_pd(cosp,t1r);
                         numi = _mm512_mul_pd(_2,t1i);
                         deni = _mm512_add_pd(cosp,t1i);
                         cdiv_zmm8r8(numr,numi,denr,deni,*Toutr,*Touti);
                 }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4174_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pmur,
                                            const double * __restrict __ATTR_ALIGN__(64) pmui,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                            double * __restrict __ATTR_ALIGN__(64) Toutr,
                                            double * __restrict __ATTR_ALIGN__(64) Touti) {

                         register __m512d mur  = _mm512_load_pd(&pmur[0]);
                         register __m512d mui  = _mm512_load_pd(&pmui[0]);
                         register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d divr,divi,sqr1,sqi1;
                         register __m512d mulr,muli,sqr2,sqi2;
                         register __m512d cosp,sinp,sin2p,t0r,t0i;
                         register __m512d _2sqr1,_2sqi1,t1r,t1i;
                         register __m512d numr,numi,denr,deni,resr,resi;
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&div);
                         cosp = xcosf(psi);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mulr,&muli);
                         sinp = xsinf(psi);
                         csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                         sin2p= _mm512_mul_pd(sinp,sinp);
                         //_2sqr1 = _mm512_mul_pd(_2,sqr1);
                         t0r    = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         //_2sqi1 = _mm512_mul_pd(_2,sqi1);
                         t0i    = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                         cmul_zmm8r8(sqr1,sqi1,t0r,t0i,&t1r,&t1i);
                         numr = _mm512_mul_pd(_2,t1r);
                         denr = _mm512_add_pd(cosp,t1r);
                         numi = _mm512_mul_pd(_2,t1i);
                         deni = _mm512_add_pd(cosp,t1i);
                         cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                         _mm512_store_pd(&Toutr[0], resr);
                         _mm512_store_pd(&Touti[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4174_zmm8r8_u(const double * __restrict pmur,
                                            const double * __restrict  pmui,
                                            const double * __restrict  pepsr,
                                            const double * __restrict  pepsi,
                                            const double * __restrict  ppsi,
                                            double * __restrict  Toutr,
                                            double * __restrict  Touti) {

                         register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                         register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d divr,divi,sqr1,sqi1;
                         register __m512d mulr,muli,sqr2,sqi2;
                         register __m512d cosp,sinp,sin2p,t0r,t0i;
                         register __m512d _2sqr1,_2sqi1,t1r,t1i;
                         register __m512d numr,numi,denr,deni,resr,resi;
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&div);
                         cosp = xcosf(psi);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mulr,&muli);
                         sinp = xsinf(psi);
                         csqrt_zmm8r8(divr,divi,&sqr1,&sqi1);
                         sin2p= _mm512_mul_pd(sinp,sinp);
                         //_2sqr1 = _mm512_mul_pd(_2,sqr1);
                         t0r    = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         //_2sqi1 = _mm512_mul_pd(_2,sqi1);
                         t0i    = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr2,&sqi2);
                         cmul_zmm8r8(sqr1,sqi1,t0r,t0i,&t1r,&t1i);
                         numr = _mm512_mul_pd(_2,t1r);
                         denr = _mm512_add_pd(cosp,t1r);
                         numi = _mm512_mul_pd(_2,t1i);
                         deni = _mm512_add_pd(cosp,t1i);
                         cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                         _mm512_storeu_pd(&Toutr[0], resr);
                         _mm512_storeu_pd(&Touti[0], resi);
                 }


                   /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-75
                       */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4175_zmm8r8(const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           const __m512d psi,
                                           __m512d * __restrict Toutr,
                                           __m512d * __restrict Touti) {

                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         register __m512d cosp,sinp,sin2p,sqr1,sqi1;
                         register __m512d sqr2,sqi2,_2cosp,divr,divi;
                         register __m512d mulr,muli,denr,deni,t0r,t0i;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_mul_pd(sinp,sinp);
                         _2cosp= _mm512_add_pd(cosp,cosp);
                         cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         t0i = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         cmul_zmm8r8(sqr1,sqi1,sqr2,sqi2,&denr,&deni);
                         denr = _mm512_add_pd(cosp,denr);
                         deni = _mm512_add_pd(cosp,deni);
                         *Toutr = _mm512_div_pd(_2cosp,denr);
                         *Touti = _mm512_div_pd(_2cosp,deni);
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4175_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pmur,
                                            const double * __restrict __ATTR_ALIGN__(64) pmui,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                            const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                            const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                            double * __restrict __ATTR_ALIGN__(64) Toutr,
                                            double * __restrict __ATTR_ALIGN__(64) Touti) {

                         register __m512d mur  = _mm512_load_pd(&pmur[0]);
                         register __m512d mui  = _mm512_load_pd(&pmui[0]);
                         register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         register __m512d cosp,sinp,sin2p,sqr1,sqi1;
                         register __m512d sqr2,sqi2,_2cosp,divr,divi;
                         register __m512d mulr,muli,denr,deni,t0r,t0i;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_mul_pd(sinp,sinp);
                         _2cosp= _mm512_add_pd(cosp,cosp);
                         cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         t0i = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         cmul_zmm8r8(sqr1,sqi1,sqr2,sqi2,&denr,&deni);
                         denr = _mm512_add_pd(cosp,denr);
                         deni = _mm512_add_pd(cosp,deni);
                         _mm512_store_pd(&Toutr[0] ,_mm512_div_pd(_2cosp,denr));
                         _mm512_store_pd(&Touti[0] ,_mm512_div_pd(_2cosp,deni));
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4175_zmm8r8_u(const double * __restrict  pmur,
                                            const double * __restrict  pmui,
                                            const double * __restrict  pepsr,
                                            const double * __restrict  pepsi,
                                            const double * __restrict  ppsi,
                                            double * __restrict  Toutr,
                                            double * __restrict  Touti) {

                         register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                         register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         register __m512d cosp,sinp,sin2p,sqr1,sqi1;
                         register __m512d sqr2,sqi2,_2cosp,divr,divi;
                         register __m512d mulr,muli,denr,deni,t0r,t0i;
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_mul_pd(sinp,sinp);
                         _2cosp= _mm512_add_pd(cosp,cosp);
                         cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         t0i = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         cmul_zmm8r8(sqr1,sqi1,sqr2,sqi2,&denr,&deni);
                         denr = _mm512_add_pd(cosp,denr);
                         deni = _mm512_add_pd(cosp,deni);
                         _mm512_storeu_pd(&Toutr[0] ,_mm512_div_pd(_2cosp,denr));
                         _mm512_storeu_pd(&Touti[0] ,_mm512_div_pd(_2cosp,deni));
                 }


                   /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-76
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4176_zmm8r8( const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           const __m512d psi,
                                           __m512d * __restrict Rinr,
                                           __m512d * __restrict Rini) {

                         register __m512d cosp,sinp,sin2p,divr,divi;
                         register __m512d mulr,muli,denr,deni,numr,numi;
                         register __m512d sqr1,sqi1,sqr2,sqi2,t0r,t0i;
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_pd(sinp,sinp);
                         cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_div_pd(sin2p,mulr));
                         t0i = _mm512_sub_pd(_1,_mm512_div_pd(sin2p,muli));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_pd(cosp,sqr2);
                         sqi2 = _mm512_mul_pd(cosp,sqi2);
                         numr = _mm512_sub_pd(sqr2,sqr1);
                         denr = _mm512_add_pd(sqr2,sqr1);
                         numi = _mm512_sub_pd(sqi2,sqi1);
                         deni = _mm512_add_pd(sqi2,sqi1);
                         cdiv_zmm8r8(numr,numi,denr,deni,*Rinr,*Rini);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4176_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                             double * __restrict __ATTR_ALIGN__(64) Rinr,
                                             double * __restrict __ATTR_ALIGN__(64) Rini) {

                         register __m512d mur  = _mm512_load_pd(&pmur[0]);
                         register __m512d mui  = _mm512_load_pd(&pmui[0]);
                         register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                         register __m512d cosp,sinp,sin2p,divr,divi;
                         register __m512d mulr,muli,denr,deni,numr,numi;
                         register __m512d sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_pd(sinp,sinp);
                         cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_div_pd(sin2p,mulr));
                         t0i = _mm512_sub_pd(_1,_mm512_div_pd(sin2p,muli));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_pd(cosp,sqr2);
                         sqi2 = _mm512_mul_pd(cosp,sqi2);
                         numr = _mm512_sub_pd(sqr2,sqr1);
                         denr = _mm512_add_pd(sqr2,sqr1);
                         numi = _mm512_sub_pd(sqi2,sqi1);
                         deni = _mm512_add_pd(sqi2,sqi1);
                         cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                         _mm512_store_pd(&Rinr[0], resr);
                         _mm512_store_pd(&Rini[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4176_zmm8r8_u( const double * __restrict  pmur,
                                             const double * __restrict  pmui,
                                             const double * __restrict  pepsr,
                                             const double * __restrict  pepsi,
                                             const double * __restrict  ppsi,
                                             double * __restrict  Rinr,
                                             double * __restrict  Rini) {

                         register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                         register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d cosp,sinp,sin2p,divr,divi;
                         register __m512d mulr,muli,denr,deni,numr,numi;
                         register __m512d sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_pd(sinp,sinp);
                         cdiv_zmm8r8(mur,mui,epsr,epsi,&divr,&divi);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_div_pd(sin2p,mulr));
                         t0i = _mm512_sub_pd(_1,_mm512_div_pd(sin2p,muli));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_pd(cosp,sqr2);
                         sqi2 = _mm512_mul_pd(cosp,sqi2);
                         numr = _mm512_sub_pd(sqr2,sqr1);
                         denr = _mm512_add_pd(sqr2,sqr1);
                         numi = _mm512_sub_pd(sqi2,sqi1);
                         deni = _mm512_add_pd(sqi2,sqi1);
                         cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                         _mm512_storeu_pd(&Rinr[0], resr);
                         _mm512_storeu_pd(&Rini[0], resi);
                 }


                    /*
                           Fresnel reflection and transmission coefficients
                           Formula 4.1-77
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4177_zmm8r8( const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           const __m512d psi,
                                           __m512d * __restrict Rinr,
                                           __m512d * __restrict Rini) {

                         register __m512d cosp,sinp,sin2p,divr,divi;
                         register __m512d mulr,muli,denr,deni,numr,numi;
                         register __m512d sqr1,sqi1,sqr2,sqi2,t0r,t0i;
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_pd(sinp,sinp);
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&divi);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         t0i = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_pd(sqr2,cosp);
                         sqi2 = _mm512_mul_pd(sqi2,cosp);
                         numr = _mm512_sub_pd(sqr2,sqr1);
                         denr = _mm512_add_pd(sqr2,sqr1);
                         numi = _mm512_sub_pd(sqi2,sqi1);
                         deni = _mm512_add_pd(sqi2,sqi1);
                         cdiv_zmm8r8(numr,numi,denr,deni,*Rinr,*Rini);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4177_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                             double * __restrict __ATTR_ALIGN__(64) Rinr,
                                             double * __restrict __ATTR_ALIGN__(64) Rini ) {

                         register __m512d mur  = _mm512_load_pd(&pmur[0]);
                         register __m512d mui  = _mm512_load_pd(&pmui[0]);
                         register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                         register __m512d cosp,sinp,sin2p,divr,divi;
                         register __m512d mulr,muli,denr,deni,numr,numi;
                         register __m512d sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_pd(sinp,sinp);
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&divi);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         t0i = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_pd(sqr2,cosp);
                         sqi2 = _mm512_mul_pd(sqi2,cosp);
                         numr = _mm512_sub_pd(sqr2,sqr1);
                         denr = _mm512_add_pd(sqr2,sqr1);
                         numi = _mm512_sub_pd(sqi2,sqi1);
                         deni = _mm512_add_pd(sqi2,sqi1);
                         cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                         _mm512_store_pd(&Rinr[0], resr);
                         _mm512_store_pd(&Rini[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rin_f4177_zmm8r8_u( const double * __restrict  pmur,
                                             const double * __restrict  pmui,
                                             const double * __restrict  pepsr,
                                             const double * __restrict  pepsi,
                                             const double * __restrict  ppsi,
                                             double * __restrict  Rinr,
                                             double * __restrict  Rini ) {

                         register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                         register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d cosp,sinp,sin2p,divr,divi;
                         register __m512d mulr,muli,denr,deni,numr,numi;
                         register __m512d sqr1,sqi1,sqr2,sqi2,t0r,t0i,resr,resi;
                         const __m512d _1 = _mm512_set1_pd(1.0f);
                         cosp = xcosf(psi);
                         sinp = xsinf(psi);
                         sin2p= _mm512_add_pd(sinp,sinp);
                         cdiv_zmm8r8(epsr,epsi,mur,mui,&divr,&divi);
                         cmul_zmm8r8(mur,mui,epsr,epsi,&mulr,&muli);
                         t0r = _mm512_sub_pd(_1,_mm512_mul_pd(mulr,sin2p));
                         t0i = _mm512_sub_pd(_1,_mm512_mul_pd(muli,sin2p));
                         csqrt_zmm8r8(t0r,t0i,&sqr1,&sqi1);
                         csqrt_zmm8r8(divr,divi,&sqr2,&sqi2);
                         sqr2 = _mm512_mul_pd(sqr2,cosp);
                         sqi2 = _mm512_mul_pd(sqi2,cosp);
                         numr = _mm512_sub_pd(sqr2,sqr1);
                         denr = _mm512_add_pd(sqr2,sqr1);
                         numi = _mm512_sub_pd(sqi2,sqi1);
                         deni = _mm512_add_pd(sqi2,sqi1);
                         cdiv_zmm8r8(numr,numi,denr,deni,&resr,&resi);
                         _mm512_storeu_pd(&Rinr[0], resr);
                         _mm512_storeu_pd(&Rini[0], resi);
                }


                  /*
                          Specular rays reflection
                          Formula 4.1-64
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rext_f4164_zmm8r8(const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           __m512d * __restrict Rexr,
                                           __m512d * __restrict Rexi) {

                        register __m512d sqr1,sqi1,sqr2,sqi2;
                        register __m512d difr,difi,sumr,sumi;
                        csqrt_zmm8r8(mur,mui,&sqr1,sqi1);
                        csqrt_zmm8r8(epsr,epsi,&sqr2,&sqi2);
                        difr = _mm512_sub_pd(sqr1,sqr2);
                        sumr = _mm512_add_pd(sqr1,sqr2);
                        difi = _mm512_sub_pd(sqi1,sqi2);
                        sumi = _mm512_add_pd(sqi1,sqi2);
                        cdiv_zmm8r8(difr,difi,sumr,sumi,*Rexr,*Rexi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rext_f4164_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             double * __restrict __ATTR_ALIGN__(64) Rexr,
                                             double * __restrict __ATTR_ALIGN__(64) Rexi) {

                        register __m512d mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d mui  = _mm512_load_pd(&pmui[0]);
                        register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d sqr1,sqi1,sqr2,sqi2;
                        register __m512d difr,difi,sumr,sumi,resr,resi;
                        csqrt_zmm8r8(mur,mui,&sqr1,sqi1);
                        csqrt_zmm8r8(epsr,epsi,&sqr2,&sqi2);
                        difr = _mm512_sub_pd(sqr1,sqr2);
                        sumr = _mm512_add_pd(sqr1,sqr2);
                        difi = _mm512_sub_pd(sqi1,sqi2);
                        sumi = _mm512_add_pd(sqi1,sqi2);
                        cdiv_zmm8r8(difr,difi,sumr,sumi,&resr,&resi);
                        _mm512_store_pd(&Rexr[0], resr);
                        _mm512_store_pd(&Rexi[0], resi);
                }

                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rext_f4164_zmm8r8_u(const double * __restrict  pmur,
                                             const double * __restrict  pmui,
                                             const double * __restrict  pepsr,
                                             const double * __restrict  pepsi,
                                             double * __restrict  Rexr,
                                             double * __restrict  Rexi) {

                        register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                        register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d sqr1,sqi1,sqr2,sqi2;
                        register __m512d difr,difi,sumr,sumi,resr,resi;
                        csqrt_zmm8r8(mur,mui,&sqr1,sqi1);
                        csqrt_zmm8r8(epsr,epsi,&sqr2,&sqi2);
                        difr = _mm512_sub_pd(sqr1,sqr2);
                        sumr = _mm512_add_pd(sqr1,sqr2);
                        difi = _mm512_sub_pd(sqi1,sqi2);
                        sumi = _mm512_add_pd(sqi1,sqi2);
                        cdiv_zmm8r8(difr,difi,sumr,sumi,&resr,&resi);
                        _mm512_storeu_pd(&Rexr[0], resr);
                        _mm512_storeu_pd(&Rexi[0], resi);
                }


                  /*

                         Axial rays, when phi = 0
                         Formula 4.1-67
                    */

                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4167_zmm8r8( const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           __m512d * __restrict Tinr,
                                           __m512d * __restrict Tini) {

                          register __m512d sqr1,sqi1,sqr2,sqi2;
                          register __m512d sumr,sumi,mu2r,mu2i;
                          csqrt_zmm8r8(mur,mui,&sqr1,&sqi1);
                          mu2r = _mm512_add_pd(sqr1,sqr1);
                          mu2i = _mm512_add_pd(sqi1,sqi1);
                          csqrt_zmm8r8(epsr,epsi,&sqr2,&sqi2);
                          sumr = _mm512_add_pd(sqr1,sqr2);
                          sumi = _mm512_add_pd(sqi1,sqi2);
                          cdiv_zmm8r8(mu2r,mu2i,sumr,sumi,*Tinr,*Tini);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4167_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             double * __restrict __ATTR_ALIGN__(64) Tinr,
                                             double * __restrict __ATTR_ALIGN__(64) Tini ) {

                          register __m512d mur  = _mm512_load_pd(&pmur[0]);
                          register __m512d mui  = _mm512_load_pd(&pmui[0]);
                          register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                          register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                          register __m512d sqr1,sqi1,sqr2,sqi2,resr,resi;
                          register __m512d sumr,sumi,mu2r,mu2i;
                          csqrt_zmm8r8(mur,mui,&sqr1,&sqi1);
                          mu2r = _mm512_add_pd(sqr1,sqr1);
                          mu2i = _mm512_add_pd(sqi1,sqi1);
                          csqrt_zmm8r8(epsr,epsi,&sqr2,&sqi2);
                          sumr = _mm512_add_pd(sqr1,sqr2);
                          sumi = _mm512_add_pd(sqi1,sqi2);
                          cdiv_zmm8r8(mu2r,mu2i,sumr,sumi,&resr,&resi);
                          _mm512_store_pd(&Tinr[0], resr);
                          _mm512_store_pd(&Tini[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tin_f4167_zmm8r8_u( const double * __restrict  pmur,
                                             const double * __restrict  pmui,
                                             const double * __restrict  pepsr,
                                             const double * __restrict  pepsi,
                                             double * __restrict  Tinr,
                                             double * __restrict  Tini ) {

                          register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                          register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                          register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                          register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                          register __m512d sqr1,sqi1,sqr2,sqi2,resr,resi;
                          register __m512d sumr,sumi,mu2r,mu2i;
                          csqrt_zmm8r8(mur,mui,&sqr1,&sqi1);
                          mu2r = _mm512_add_pd(sqr1,sqr1);
                          mu2i = _mm512_add_pd(sqi1,sqi1);
                          csqrt_zmm8r8(epsr,epsi,&sqr2,&sqi2);
                          sumr = _mm512_add_pd(sqr1,sqr2);
                          sumi = _mm512_add_pd(sqi1,sqi2);
                          cdiv_zmm8r8(mu2r,mu2i,sumr,sumi,&resr,&resi);
                          _mm512_storeu_pd(&Tinr[0], resr);
                          _mm512_storeu_pd(&Tini[0], resi);
                }


                  /*
                          Axial rays, when phi = 0
                          Formula 4.1-68
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4168_zmm8r8( const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           __m512d * __restrict Toutr,
                                           __m512d * __restrict Touti) {

                          register __m512d sqr1,sqi1,sqr2,sqi2;
                          register __m512d sumr,sumi,eps2r,eps2i;
                          csqrt_zmm8r8(epsr,epsi,&sqr1,&sqi1);
                          eps2r = _mm512_add_pd(sqr1,sqr1);
                          eps2i = _mm512_add_pd(sqi1,sqi1);
                          csqrt_zmm8r8(mur,mui,&sqr2,&sqi2);
                          sumr = _mm512_add_pd(sqr1,sqr2);
                          sumi = _mm512_add_pd(sqi1,sqi2);
                          cdiv_zmm8r8(eps2r,eps2i,sumr,sumi,*Toutr,*Touti);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4168_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             double * __restrict __ATTR_ALIGN__(64) Toutr,
                                             double * __restrict __ATTR_ALIGN__(64) Touti ) {

                          register __m512d mur  = _mm512_load_pd(&pmur[0]);
                          register __m512d mui  = _mm512_load_pd(&pmui[0]);
                          register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                          register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                          register __m512d sqr1,sqi1,sqr2,sqi2;
                          register __m512d sumr,sumi,eps2r,eps2i,resr,resi;
                          csqrt_zmm8r8(epsr,epsi,&sqr1,&sqi1);
                          eps2r = _mm512_add_pd(sqr1,sqr1);
                          eps2i = _mm512_add_pd(sqi1,sqi1);
                          csqrt_zmm8r8(mur,mui,&sqr2,&sqi2);
                          sumr = _mm512_add_pd(sqr1,sqr2);
                          sumi = _mm512_add_pd(sqi1,sqi2);
                          cdiv_zmm8r8(eps2r,eps2i,sumr,sumi,&resr,&resi);
                          _mm512_store_pd(&Toutr[0], resr);
                          _mm512_store_pd(&Touti[0], resi);
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Tout_f4168_zmm8r8_u(  const double * __restrict  pmur,
                                             const double * __restrict  pmui,
                                             const double * __restrict  pepsr,
                                             const double * __restrict  pepsi,
                                             double * __restrict  Toutr,
                                             double * __restrict Touti ) {

                          register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                          register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                          register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                          register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                          register __m512d sqr1,sqi1,sqr2,sqi2;
                          register __m512d sumr,sumi,eps2r,eps2i,resr,resi;
                          csqrt_zmm8r8(epsr,epsi,&sqr1,&sqi1);
                          eps2r = _mm512_add_pd(sqr1,sqr1);
                          eps2i = _mm512_add_pd(sqi1,sqi1);
                          csqrt_zmm8r8(mur,mui,&sqr2,&sqi2);
                          sumr = _mm512_add_pd(sqr1,sqr2);
                          sumi = _mm512_add_pd(sqi1,sqi2);
                          cdiv_zmm8r8(eps2r,eps2i,sumr,sumi,&resr,&resi);
                          _mm512_storeu_pd(&Toutr[0], resr);
                          _mm512_storeu_pd(&Touti[0], resi);
                }


                  /*
                          Axial rays, when phi = 0
                          Formula 4.1-69
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rint_f4169_zmm8r8(const __m512d mur,
                                           const __m512d mui,
                                           const __m512d epsr,
                                           const __m512d epsi,
                                           __m512d * __restrict Rintr,
                                           __m512d * __restrict Rinti) {
                        
                        register __m512d t0r,t0i;
                        const __m512d n1 = _mm512_mul_pd(-1.0f);
                        Rext_f4164_zmm8r8(mur,mui,epsr,epsi,&t0r,&t0i);
                        *Rintr = _mm512_mul_pd(n1,t0r);
                        *Rinti = _mm512_mul_pd(n1,t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rint_f4169_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pmur,
                                             const double * __restrict __ATTR_ALIGN__(64) pmui,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                             const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                             double * __restrict __ATTR_ALIGN__(64) Rintr,
                                             double * __restrict __ATTR_ALIGN__(64) Rinti) {
                        
                        register __m512d mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d mui  = _mm512_load_pd(&pmui[0]);
                        register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d t0r,t0i;
                        const __m512d n1 = _mm512_mul_pd(-1.0f);
                        Rext_f4164_zmm8r8(mur,mui,epsr,epsi,&t0r,&t0i);
                        _mm512_store_pd(&Rintr[0] ,_mm512_mul_pd(n1,t0r));
                        _mm512_store_pd(&Rinti[0] ,_mm512_mul_pd(n1,t0i));
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Rint_f4169_zmm8r8_u( const double * __restrict pmur,
                                              const double * __restrict  pmui,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              double * __restrict  Rintr,
                                              double * __restrict  Rinti) {
                        
                        register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                        register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d t0r,t0i;
                        const __m512d n1 = _mm512_mul_pd(-1.0f);
                        Rext_f4164_zmm8r8(mur,mui,epsr,epsi,&t0r,&t0i);
                        _mm512_storeu_pd(&Rintr[0] ,_mm512_mul_pd(n1,t0r));
                        _mm512_storeu_pd(&Rinti[0] ,_mm512_mul_pd(n1,t0i));
                 }


                   /*
                       Backscatter widths in high-frequency limit.
                       Phi = 0, formula 4.1-91,for k1a>5.
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4191_zmm8r8(const __m512d a,
                                            const __m512d mur,
                                            const __m512d mui,
                                            const __m512d epsr,
                                            const __m512d epsi ) {

                          register __m512d t0r,t0i;
                          register __m512d cabs,rcs;
                          Rext_f4164_zmm8r8(mur,mui,epsr,epsi,&t0r,&t0i);
                          cabs = cabs_zmm8r8(t0r,t0i);
                          rcs  = _mm512_mul_pd(cabs,_mm512_mul_pd(PI,a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4191_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi, ) {

                          register __m512d mur  = _mm512_load_pd(&pmur[0]);
                          register __m512d mui  = _mm512_load_pd(&pmui[0]);
                          register __m512d epsr = _mm512_load_pd(&pepsr[0]);
                          register __m512d epsi = _mm512_load_pd(&pepsi[0]);
                          register __m512d t0r,t0i;
                          register __m512d cabs,rcs;
                          Rext_f4164_zmm8r8(mur,mui,epsr,epsi,&t0r,&t0i);
                          cabs = cabs_zmm8r8(t0r,t0i);
                          rcs  = _mm512_mul_pd(cabs,_mm512_mul_pd(PI,a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4191_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi, ) {

                          register __m512d mur  = _mm512_loadu_pd(&pmur[0]);
                          register __m512d mui  = _mm512_loadu_pd(&pmui[0]);
                          register __m512d epsr = _mm512_loadu_pd(&pepsr[0]);
                          register __m512d epsi = _mm512_loadu_pd(&pepsi[0]);
                          register __m512d t0r,t0i;
                          register __m512d cabs,rcs;
                          Rext_f4164_zmm8r8(mur,mui,epsr,epsi,&t0r,&t0i);
                          cabs = cabs_zmm8r8(t0r,t0i);
                          rcs  = _mm512_mul_pd(cabs,_mm512_mul_pd(PI,a));
                          return (rcs);
                 }


                    /*
                         Bistatic scattering width (k0a0<<1, k1a0<<1), function of phi angle.
                         Formula 4.1-104
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41104_zmm8r8(const __m512d a0,
                                             const __m512d a1,
                                             const __m512d k0a0,
                                             const __m512d phi,
                                             const __m512d mu1r,
                                             const __m512d mu1i,
                                             const __m512d mu0r,
                                             const __m512d mu0i,
                                             const __m512d eps1r,
                                             const __m512d eps1i,
                                             const __m512d eps0r,
                                             const __m512d eps0i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d cosp,divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          cosp  = xcosf(phi);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                          div2i = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                          t1r   = _mm512_sub_pd(t0r,div2r);
                          t1i   = _mm512_sub_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41104_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa0,
                                               const double * __restrict __ATTR_ALIGN__(64) pa1,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                               const double * __restrict __ATTR_ALIGN__(64) pphi,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu1r,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu1i,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu0r,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu0i,
                                               const double * __restrict __ATTR_ALIGN__(64) peps1r,
                                               const double * __restrict __ATTR_ALIGN__(64) peps1i,
                                               const double * __restrict __ATTR_ALIGN__(64) peps0r,
                                               const double * __restrict __ATTR_ALIGN__(64) peps0i) {

                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d phi   = _mm512_load_pd(&pphi[0]);
                          register __m512d mu1r  = _mm512_load_pd(&pmu1r[0]);
                          register __m512d mu1i  = _mm512_load_pd(&pmu1i[0]);
                          register __m512d mu0r  = _mm512_load_pd(&pmu0r[0]);
                          register __m512d mu0i  = _mm512_load_pd(&pmu0i[0]);
                          register __m512d eps1r = _mm512_load_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_load_pd(&peps1i[0]);
                          register __m512d eps0r = _mm512_load_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_load_pd(&peps1r[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d cosp,divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          cosp  = xcosf(phi);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                          div2i = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                          t1r   = _mm512_sub_pd(t0r,div2r);
                          t1i   = _mm512_sub_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41104_zmm8r8_u(const double * __restrict  pa0,
                                               const double * __restrict  pa1,
                                               const double * __restrict  pk0a0,
                                               const double * __restrict  pphi,
                                               const double * __restrict  pmu1r,
                                               const double * __restrict  pmu1i,
                                               const double * __restrict  pmu0r,
                                               const double * __restrict  pmu0i,
                                               const double * __restrict  peps1r,
                                               const double * __restrict  peps1i,
                                               const double * __restrict  peps0r,
                                               const double * __restrict  peps0i) {

                          register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                          register __m512d a1    = _mm512_loadu_pd(&pa1[0]);
                          register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                          register __m512d phi   = _mm512_loadu_pd(&pphi[0]);
                          register __m512d mu1r  = _mm512_loadu_pd(&pmu1r[0]);
                          register __m512d mu1i  = _mm512_loadu_pd(&pmu1i[0]);
                          register __m512d mu0r  = _mm512_loadu_pd(&pmu0r[0]);
                          register __m512d mu0i  = _mm512_loadu_pd(&pmu0i[0]);
                          register __m512d eps1r = _mm512_loadu_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_loadu_pd(&peps1i[0]);
                          register __m512d eps0r = _mm512_loadu_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_loadu_pd(&peps1r[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d cosp,divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          cosp  = xcosf(phi);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                          div2i = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                          t1r   = _mm512_sub_pd(t0r,div2r);
                          t1i   = _mm512_sub_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                /*
                         Backscattering  width (k0a0<<1, k1a0<<1), phi = 0
                         Formula 4.1-105
                  */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41105_zmm8r8(const __m512d a0,
                                             const __m512d a1,
                                             const __m512d k0a0,
                                             const __m512d mu1r,
                                             const __m512d mu1i,
                                             const __m512d mu0r,
                                             const __m512d mu0i,
                                             const __m512d eps1r,
                                             const __m512d eps1i,
                                             const __m512d eps0r,
                                             const __m512d eps0i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,div2r);
                          div2i = _mm512_mul_pd(_2,div2r);
                          t1r   = _mm512_sub_pd(t0r,div2r);
                          t1i   = _mm512_sub_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41105_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa0,
                                               const double * __restrict __ATTR_ALIGN__(64) pa1,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu1r,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu1i,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu0r,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu0i,
                                               const double * __restrict __ATTR_ALIGN__(64) peps1r,
                                               const double * __restrict __ATTR_ALIGN__(64) peps1i,
                                               const double * __restrict __ATTR_ALIGN__(64) peps0r,
                                               const double * __restrict __ATTR_ALIGN__(64) peps0i) {

                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d mu1r  = _mm512_load_pd(&pmu1r[0]);
                          register __m512d mu1i  = _mm512_load_pd(&pmu1i[0]);
                          register __m512d mu0r  = _mm512_load_pd(&pmu0r[0]);
                          register __m512d mu0i  = _mm512_load_pd(&pmu0i[0]);
                          register __m512d eps1r = _mm512_load_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_load_pd(&peps1i[0]);
                          register __m512d eps0r = _mm512_load_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_load_pd(&peps1r[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,div2r);
                          div2i = _mm512_mul_pd(_2,div2r);
                          t1r   = _mm512_sub_pd(t0r,div2r);
                          t1i   = _mm512_sub_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41105_zmm8r8_u(const double * __restrict  pa0,
                                               const double * __restrict  pa1,
                                               const double * __restrict  pk0a0,
                                               const double * __restrict  pmu1r,
                                               const double * __restrict  pmu1i,
                                               const double * __restrict  pmu0r,
                                               const double * __restrict  pmu0i,
                                               const double * __restrict  peps1r,
                                               const double * __restrict  peps1i,
                                               const double * __restrict  peps0r,
                                               const double * __restrict  peps0i) {

                          register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                          register __m512d a1    = _mm512_loadu_pd(&pa1[0]);
                          register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                          register __m512d mu1r  = _mm512_loadu_pd(&pmu1r[0]);
                          register __m512d mu1i  = _mm512_loadu_pd(&pmu1i[0]);
                          register __m512d mu0r  = _mm512_loadu_pd(&pmu0r[0]);
                          register __m512d mu0i  = _mm512_loadu_pd(&pmu0i[0]);
                          register __m512d eps1r = _mm512_loadu_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_loadu_pd(&peps1i[0]);
                          register __m512d eps0r = _mm512_loadu_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_loadu_pd(&peps1r[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,div2r);
                          div2i = _mm512_mul_pd(_2,div2r);
                          t1r   = _mm512_sub_pd(t0r,div2r);
                          t1i   = _mm512_sub_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                /*
                      Forward scattering width (k0a0<<1, k1a0<<1), phi = pi.
                      Formula 4.1-106
                 */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41106_zmm8r8(const __m512d a0,
                                             const __m512d a1,
                                             const __m512d k0a0,
                                             const __m512d mu1r,
                                             const __m512d mu1i,
                                             const __m512d mu0r,
                                             const __m512d mu0i,
                                             const __m512d eps1r,
                                             const __m512d eps1i,
                                             const __m512d eps0r,
                                             const __m512d eps0i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,div2r);
                          div2i = _mm512_mul_pd(_2,div2r);
                          t1r   = _mm512_add_pd(t0r,div2r);
                          t1i   = _mm512_add_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41106_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa0,
                                               const double * __restrict __ATTR_ALIGN__(64) pa1,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu1r,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu1i,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu0r,
                                               const double * __restrict __ATTR_ALIGN__(64) pmu0i,
                                               const double * __restrict __ATTR_ALIGN__(64) peps1r,
                                               const double * __restrict __ATTR_ALIGN__(64) peps1i,
                                               const double * __restrict __ATTR_ALIGN__(64) peps0r,
                                               const double * __restrict __ATTR_ALIGN__(64) peps0i) {

                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d mu1r  = _mm512_load_pd(&pmu1r[0]);
                          register __m512d mu1i  = _mm512_load_pd(&pmu1i[0]);
                          register __m512d mu0r  = _mm512_load_pd(&pmu0r[0]);
                          register __m512d mu0i  = _mm512_load_pd(&pmu0i[0]);
                          register __m512d eps1r = _mm512_load_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_load_pd(&peps1i[0]);
                          register __m512d eps0r = _mm512_load_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_load_pd(&peps1r[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d pia,k0a03,rcs,a1a0,_1pa,_1ma,cabs;
                          register __m512d divr,divi,e1mr,e1mi,t1r,t1i;
                          register __m512d e0mr,e0mi,frac,a1a0s,t0r,t0i;
                          register __m512d div2r,div2i,numr,numi,denr,deni;
                          pia = _mm512_mul_pd(PI,a0);
                          k0a03 = _mm512_mul_pd(k0a0,
                                            _mm512_mul_pd(k0a0,k0a0));
                          frac  = _mm512_mul_pd(pia,_mm512_mul_pd(pi4,k0a03));
                          a1a0  = _mm512_div_pd(a1,a0);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          _1pa  = _mm512_add_pd(_1,a1a0s);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          t0r   = _mm512_sub_pd(_mm512_mul_pd(divr,_1ma),_1);
                          t0i   = _mm512_sub_pd(_mm512_mul_pd(divi,_1ma),_1);
                          e1mr  = _mm512_mul_pd(eps1r,_1pa);
                          e0mr  = _mm512_mul_pd(eps0r,_1ma);
                          e1mi  = _mm512_mul_pd(eps1i,_1pa);
                          e0mi  = _mm512_mul_pd(eps0i,_1ma);
                          numr  = _mm512_sub_pd(e1mr,e0mr);
                          numi  = _mm512_sub_pd(e1mi,e0mi);
                          denr  = _mm512_add_pd(e1mr,e0mr);
                          deni  = _mm512_add_pd(e1mi,e0mi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&div2r,&div2i);
                          div2r = _mm512_mul_pd(_2,div2r);
                          div2i = _mm512_mul_pd(_2,div2r);
                          t1r   = _mm512_add_pd(t0r,div2r);
                          t1i   = _mm512_add_pd(t0i,div2i);
                          cabs  = cabs_zmm8r8(t1r,t1i);
                          rcs   = _mm512_mul_pd(frac,cabs);
                          return (rcs);
               }



                 /*
                       Hollow cylindrical shell.
                       Approximations for the low frequency region
                       (k0a0<<1, k1a0<<1).
                       Formula 4.1-124
                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41124_zmm8r8(const __m512d a1,
                                          const __m512d a0,
                                          const __m512d k0a0,
                                          const __m512d eps1r,
                                          const __m512d eps1i,
                                          const __m512d eps0r,
                                          const __m512d eps0i,
                                          __m512d * __restrict A0r,
                                          __m512d * __restrict A0i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512d t0r,t0i;
                          a1a0  = _mm512_div_pd(a1,a0);
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          fraci = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          cdiv_zmm8r8(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          divr = _mm512_sub_pd(divr,_1);
                          divi = _mm512_sub_pd(divi,_1);
                          t0r  = _mm512_mul_pd(divr,_1ma);
                          t0i  = _mm512_mul_pd(divi,_1ma);
                          cmul_zmm8r8(fracr,fraci,t0r,t0i,*A0r,*A0i);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41124_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  double * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps1r,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps1i,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps0r,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps0i,
                                            double * __restrict __ATTR_ALIGN__(64) A0r,
                                            double * __restrict __ATTR_ALIGN__(64) A0i) {

                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d eps1r = _mm512_load_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_load_pd(&peps1i[0]); 
                          register __m512d eps0r = _mm512_load_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_load_pd(&peps0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512d t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_pd(a1,a0);
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          fraci = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          cdiv_zmm8r8(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          divr = _mm512_sub_pd(divr,_1);
                          divi = _mm512_sub_pd(divi,_1);
                          t0r  = _mm512_mul_pd(divr,_1ma);
                          t0i  = _mm512_mul_pd(divi,_1ma);
                          cmul_zmm8r8(fracr,fraci,t0r,t0i,&resr,&resi);
                          _mm512_store_pd(&A0r[0], resr);
                          _mm512_store_pd(&A0i[0], resi);

               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41124_zmm8r8_u(const  double * __restrict  pa1,
                                          const  double * __restrict pa0,
                                          const  double * __restrict  pk0a0,
                                          const  double * __restrict  peps1r,
                                          const  double * __restrict  peps1i,
                                          const  double * __restrict  peps0r,
                                          const  double * __restrict  peps0i,
                                          double * __restrict  A0r,
                                          double * __restrict  A0i) {

                          register __m512d a1    = _mm512_loadu_pd(&pa1[0]);
                          register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                          register __m512d eps1r = _mm512_loadu_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_loadu_pd(&peps1i[0]); 
                          register __m512d eps0r = _mm512_loadu_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_loadu_pd(&peps0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512d t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_pd(a1,a0);
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          fraci = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          cdiv_zmm8r8(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          divr = _mm512_sub_pd(divr,_1);
                          divi = _mm512_sub_pd(divi,_1);
                          t0r  = _mm512_mul_pd(divr,_1ma);
                          t0i  = _mm512_mul_pd(divi,_1ma);
                          cmul_zmm8r8(fracr,fraci,t0r,t0i,&resr,&resi);
                          _mm512_storeu_pd(&A0r[0], resr);
                          _mm512_storeu_pd(&A0i[0], resi);

               }


                  /*

                       Hollow cylindrical shell.
                       Approximations for the low frequency region
                       (k0a0<<1, k1a0<<1).
                       Formula 4.1-126
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41126_zmm8r8(const __m512d a1,
                                          const __m512d a0,
                                          const __m512d k0a0,
                                          const __m512d mu1r,
                                          const __m512d mu1i,
                                          const __m512d mu0r,
                                          const __m512d mu0i,
                                          __m512d * __restrict B0r,
                                          __m512d * __restrict B0i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512d t0r,t0i;
                          a1a0  = _mm512_div_pd(a1,a0);
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          fraci = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          divr = _mm512_sub_pd(divr,_1);
                          divi = _mm512_sub_pd(divi,_1);
                          t0r  = _mm512_mul_pd(divr,_1ma);
                          t0i  = _mm512_mul_pd(divi,_1ma);
                          cmul_zmm8r8(fracr,fraci,t0r,t0i,*B0r,*B0i);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41126_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  double * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu1r,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu1i,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu0r,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu0i,
                                            double * __restrict __ATTR_ALIGN__(64) B0r,
                                            double * __restrict __ATTR_ALIGN__(64) B0i) {

                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d mu1r = _mm512_load_pd(&pmu1r[0]);
                          register __m512d mu1i = _mm512_load_pd(&pmu1i[0]); 
                          register __m512d mu0r = _mm512_load_pd(&pmu0r[0]);
                          register __m512d mu0i = _mm512_load_pd(&pmu0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512d t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_pd(a1,a0);
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          fraci = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          divr = _mm512_sub_pd(divr,_1);
                          divi = _mm512_sub_pd(divi,_1);
                          t0r  = _mm512_mul_pd(divr,_1ma);
                          t0i  = _mm512_mul_pd(divi,_1ma);
                          cmul_zmm8r8(fracr,fraci,t0r,t0i,*resr,*resi);
                          _mm512_store_pd(&B0r[0], resr);
                          _mm512_store_pd(&B0i[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41126_zmm8r8_u(const  double * __restrict  pa1,
                                            const  double * __restrict  pa0,
                                            const  double * __restrict  pk0a0,
                                            const  double * __restrict  mu1r,
                                            const  double * __restrict  mu1i,
                                            const  double * __restrict  mu0r,
                                            const  double * __restrict  mu0i,
                                            double * __restrict  B0r,
                                            double * __restrict  B0i) {

                          register __m512d a1    = _mm512_loadu_pd(&pa1[0]);
                          register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                          register __m512d mu1r = _mm512_loadu_pd(&pmu1r[0]);
                          register __m512d mu1i = _mm512_loadu_pd(&pmu1i[0]); 
                          register __m512d mu0r = _mm512_loadu_pd(&pmu0r[0]);
                          register __m512d mu0i = _mm512_loadu_pd(&pmu0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,fracr,fraci,divr,divi,a1a0,a1a0s,_1ma;
                          register __m512d t0r,t0i,resr,resi;
                          a1a0  = _mm512_div_pd(a1,a0);
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          fracr = Ir;
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          fraci = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          divr = _mm512_sub_pd(divr,_1);
                          divi = _mm512_sub_pd(divi,_1);
                          t0r  = _mm512_mul_pd(divr,_1ma);
                          t0i  = _mm512_mul_pd(divi,_1ma);
                          cmul_zmm8r8(fracr,fraci,t0r,t0i,*resr,*resi);
                          _mm512_storeu_pd(&B0r[0], resr);
                          _mm512_storeu_pd(&B0i[0], resi);
               }


                  /*

                          Hollow cylindrical shell.
                          Approximations for the low frequency region
                          (k0a0<<1, k1a0<<1).
                           Formula 4.1-125
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41125_zmm8r8(const __m512d a1,
                                          const __m512d a0,
                                          const __m512d k0a0,
                                          const __m512d mu1r,
                                          const __m512d mu1i,
                                          const __m512d mu0r,
                                          const __m512d mu0i,
                                          __m512d * __restrict A1r,
                                          __m512d * __restrict A1i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512d divr,divi,divrs,divis,t1r,t1i;
                          register __m512d sqpr,sqpi,sqmr,sqmi,t0r,t0i;
                          register __m512d numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_pd(a1,a0);
                          fracr = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          cmul_zmm8r8(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_pd(divrs,_1);
                          divis = _mm512_sub_pd(divis,_1);
                          numr  = _mm512_mul_pd(divrs,_1ma);
                          numi  = _mm512_mul_pd(divis,_1ma);
                          t0r   = _mm512_add_pd(divr,_1);
                          t0i   = _mm512_add_pd(divi,_1);
                          cmul_zmm8r8(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_pd(divr,_1);
                          t1i   = _mm512_sub_pd(divi,_1);
                          cmul_zmm8r8(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_pd(sqmr,a1a02);
                          sqmi = _mm512_mul_pd(sqmi,a1a02);
                          denr = _mm512_sub_pd(sqpr,sqmr);
                          deni = _mm512_sub_pd(sqpi,sqmi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm8r8(fracr,fraci,ratr,rati,*A1r,*A1i);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41125_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  double * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu1r,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu1i,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu0r,
                                            const  double * __restrict __ATTR_ALIGN__(64) pmu0i,
                                            double * __restrict __ATTR_ALIGN__(64) A1r,
                                            double * __restrict __ATTR_ALIGN__(64) A1i) {

                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d mu1r = _mm512_load_pd(&pmu1r[0]);
                          register __m512d mu1i = _mm512_load_pd(&pmu1i[0]); 
                          register __m512d mu0r = _mm512_load_pd(&pmu0r[0]);
                          register __m512d mu0i = _mm512_load_pd(&pmu0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512d divr,divi,divrs,divis,t1r,t1i;
                          register __m512d sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512d numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_pd(a1,a0);
                          fracr = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          cmul_zmm8r8(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_pd(divrs,_1);
                          divis = _mm512_sub_pd(divis,_1);
                          numr  = _mm512_mul_pd(divrs,_1ma);
                          numi  = _mm512_mul_pd(divis,_1ma);
                          t0r   = _mm512_add_pd(divr,_1);
                          t0i   = _mm512_add_pd(divi,_1);
                          cmul_zmm8r8(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_pd(divr,_1);
                          t1i   = _mm512_sub_pd(divi,_1);
                          cmul_zmm8r8(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_pd(sqmr,a1a02);
                          sqmi = _mm512_mul_pd(sqmi,a1a02);
                          denr = _mm512_sub_pd(sqpr,sqmr);
                          deni = _mm512_sub_pd(sqpi,sqmi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm8r8(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_store_pd(&A1r[0], resr);
                          _mm512_store_pd(&A1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41125_zmm8r8_u(const  double * __restrict  pa1,
                                            const  double * __restrict  pa0,
                                            const  double * __restrict  pk0a0,
                                            const  double * __restrict  pmu1r,
                                            const  double * __restrict  pmu1i,
                                            const  double * __restrict  pmu0r,
                                            const  double * __restrict  pmu0i,
                                            double * __restrict  A1r,
                                            double * __restrict A1i) {

                          register __m512d a1    = _mm512_loadu_pd(&pa1[0]);
                          register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                          register __m512d mu1r = _mm512_loadu_pd(&pmu1r[0]);
                          register __m512d mu1i = _mm512_loadu_pd(&pmu1i[0]); 
                          register __m512d mu0r = _mm512_loadu_pd(&pmu0r[0]);
                          register __m512d mu0i = _mm512_loadu_pd(&pmu0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512d divr,divi,divrs,divis,t1r,t1i;
                          register __m512d sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512d numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_pd(a1,a0);
                          fracr = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          cdiv_zmm8r8(mu1r,mu1i,mu0r,mu0i,&divr,&divi);
                          cmul_zmm8r8(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_pd(divrs,_1);
                          divis = _mm512_sub_pd(divis,_1);
                          numr  = _mm512_mul_pd(divrs,_1ma);
                          numi  = _mm512_mul_pd(divis,_1ma);
                          t0r   = _mm512_add_pd(divr,_1);
                          t0i   = _mm512_add_pd(divi,_1);
                          cmul_zmm8r8(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_pd(divr,_1);
                          t1i   = _mm512_sub_pd(divi,_1);
                          cmul_zmm8r8(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_pd(sqmr,a1a02);
                          sqmi = _mm512_mul_pd(sqmi,a1a02);
                          denr = _mm512_sub_pd(sqpr,sqmr);
                          deni = _mm512_sub_pd(sqpi,sqmi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm8r8(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_storeu_pd(&A1r[0], resr);
                          _mm512_storeu_pd(&A1i[0], resi);
                }


                 /*

                          Hollow cylindrical shell.
                          Approximations for the low frequency region
                          (k0a0<<1, k1a0<<1).
                           Formula 4.1-127
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41127_zmm8r8(const __m512d a1,
                                          const __m512d a0,
                                          const __m512d k0a0,
                                          const __m512d eps1r,
                                          const __m512d eps1i,
                                          const __m512d eps0r,
                                          const __m512d eps0i,
                                          __m512d * __restrict B1r,
                                          __m512d * __restrict B1i) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512d divr,divi,divrs,divis,t1r,t1i;
                          register __m512d sqpr,sqpi,sqmr,sqmi,t0r,t0i;
                          register __m512d numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_pd(a1,a0);
                          fracr = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          cdiv_zmm8r8(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          cmul_zmm8r8(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_pd(divrs,_1);
                          divis = _mm512_sub_pd(divis,_1);
                          numr  = _mm512_mul_pd(divrs,_1ma);
                          numi  = _mm512_mul_pd(divis,_1ma);
                          t0r   = _mm512_add_pd(divr,_1);
                          t0i   = _mm512_add_pd(divi,_1);
                          cmul_zmm8r8(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_pd(divr,_1);
                          t1i   = _mm512_sub_pd(divi,_1);
                          cmul_zmm8r8(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_pd(sqmr,a1a02);
                          sqmi = _mm512_mul_pd(sqmi,a1a02);
                          denr = _mm512_sub_pd(sqpr,sqmr);
                          deni = _mm512_sub_pd(sqpi,sqmi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm8r8(fracr,fraci,ratr,rati,*B1r,*B1i);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41127_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) pa1,
                                            const  double * __restrict __ATTR_ALIGN__(64) pa0,
                                            const  double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps1r,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps1i,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps0r,
                                            const  double * __restrict __ATTR_ALIGN__(64) peps0i,
                                            double * __restrict __ATTR_ALIGN__(64) B1r,
                                            double * __restrict __ATTR_ALIGN__(64) B1i) {

                          register __m512d a1    = _mm512_load_pd(&pa1[0]);
                          register __m512d a0    = _mm512_load_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                          register __m512d eps1r = _mm512_load_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_load_pd(&peps1i[0]); 
                          register __m512d eps0r = _mm512_load_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_load_pd(&peps0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512d divr,divi,divrs,divis,t1r,t1i;
                          register __m512d sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512d numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_pd(a1,a0);
                          fracr = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          cdiv_zmm8r8(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          cmul_zmm8r8(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_pd(divrs,_1);
                          divis = _mm512_sub_pd(divis,_1);
                          numr  = _mm512_mul_pd(divrs,_1ma);
                          numi  = _mm512_mul_pd(divis,_1ma);
                          t0r   = _mm512_add_pd(divr,_1);
                          t0i   = _mm512_add_pd(divi,_1);
                          cmul_zmm8r8(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_pd(divr,_1);
                          t1i   = _mm512_sub_pd(divi,_1);
                          cmul_zmm8r8(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_pd(sqmr,a1a02);
                          sqmi = _mm512_mul_pd(sqmi,a1a02);
                          denr = _mm512_sub_pd(sqpr,sqmr);
                          deni = _mm512_sub_pd(sqpi,sqmi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm8r8(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_store_pd(&B1r[0], resr);
                          _mm512_store_pd(&B1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41127_zmm8r8_u(const  double * __restrict  pa1,
                                            const  double * __restrict  pa0,
                                            const  double * __restrict  pk0a0,
                                            const  double * __restrict  peps1r,
                                            const  double * __restrict  peps1i,
                                            const  double * __restrict  peps0r,
                                            const  double * __restrict  peps0i,
                                            double * __restrict  B1r,
                                            double * __restrict  B1i) {

                          register __m512d a1    = _mm512_loadu_pd(&pa1[0]);
                          register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                          register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                          register __m512d eps1r = _mm512_loadu_pd(&peps1r[0]);
                          register __m512d eps1i = _mm512_loadu_pd(&peps1i[0]); 
                          register __m512d eps0r = _mm512_loadu_pd(&peps0r[0]);
                          register __m512d eps0i = _mm512_loadu_pd(&peps0i[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          const __m512d pi4= _mm512_set1_pd(0.78539816339744830961566084582);
                          register __m512d k0a02,a1a0,a1a0s,_1ma,ratr,rati;
                          register __m512d divr,divi,divrs,divis,t1r,t1i;
                          register __m512d sqpr,sqpi,sqmr,sqmi,t0r,t0i,resr,resi;
                          register __m512d numr,numi,denr,deni,fracr,fraci;
                          k0a02 = _mm512_mul_pd(k0a,k0a);
                          fraci = Ir;
                          a1a0  = _mm512_div_pd(a1,a0);
                          fracr = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a02));
                          a1a0s = _mm512_mul_pd(a1a0,a1a0);
                          _1ma  = _mm512_sub_pd(_1,a1a0s);
                          cdiv_zmm8r8(eps1r,eps1i,eps0r,eps0i,&divr,&divi);
                          cmul_zmm8r8(divr,divi,divr,divi,&divrs,&divis);
                          divrs = _mm512_sub_pd(divrs,_1);
                          divis = _mm512_sub_pd(divis,_1);
                          numr  = _mm512_mul_pd(divrs,_1ma);
                          numi  = _mm512_mul_pd(divis,_1ma);
                          t0r   = _mm512_add_pd(divr,_1);
                          t0i   = _mm512_add_pd(divi,_1);
                          cmul_zmm8r8(t0r,t0i,t0r,t0i,&sqpr,&sqpi);
                          t1r   = _mm512_sub_pd(divr,_1);
                          t1i   = _mm512_sub_pd(divi,_1);
                          cmul_zmm8r8(t1r,t1i,t1r,t1i,&sqmr,&sqmi);
                          sqmr = _mm512_mul_pd(sqmr,a1a02);
                          sqmi = _mm512_mul_pd(sqmi,a1a02);
                          denr = _mm512_sub_pd(sqpr,sqmr);
                          deni = _mm512_sub_pd(sqpi,sqmi);
                          cdiv_zmm8r8(numr,numi,denr,deni,&ratr,&rati);
                          cmul_zmm8r8(fracr,fraci,ratr,rati,&resr,&resi);
                          _mm512_storeu_pd(&B1r[0], resr);
                          _mm512_storeu_pd(&B1i[0], resi);
                }


                    /*

                          Low-frequncy approximations (k0a<0.2)
                          Cylindrical Luneberg lens (k0a<0.2).
                          Formula 4.1-162
                     */

                    

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41162_zmm8r8(const __m512d k0a,
                                          __m512d * __restrict A0r,
                                          __m512d * __restrict A0i) {

                         const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf = _mm512_set1_pd(0.5f);
                         register __m512d k0a2,k0ah;
                         k0a2 = _mm512_mul_pd(k0a,k0a);
                         *A0r = Ir;
                         k0ah = _mm512_mul_pd(hlf,k0a2);
                         *A0i = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0ah));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41162_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            double * __restrict __ATTR_ALIGN__(64) A0r,
                                            double * __restrict __ATTR_ALIGN__(64) A0i) {

                   
                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf = _mm512_set1_pd(0.5f);
                         register __m512d k0a2,k0ah;
                         k0a2 = _mm512_mul_pd(k0a,k0a);
                         _mm512_store_pd(&A0r[0] ,Ir);
                         k0ah = _mm512_mul_pd(hlf,k0a2);
                         _mm512_store_pd(&A0i[0] ,_mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0ah)));
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41162_zmm8r8_u(const double * __restrict  pk0a,
                                            double * __restrict  A0r,
                                            double * __restrict  A0i) {

                   
                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf = _mm512_set1_pd(0.5f);
                         register __m512d k0a2,k0ah;
                         k0a2 = _mm512_mul_pd(k0a,k0a);
                         _mm512_storeu_pd(&A0r[0] ,Ir);
                         k0ah = _mm512_mul_pd(hlf,k0a2);
                         _mm512_storeu_pd(&A0i[0] ,_mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0ah)));
                }


                 

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41162_zmm8r8(__m512d * __restrict B0r,
                                          __m512d * __restrict B0i) {

                        *B0r = _mm512_setzero_pd();
                        *B0i = _mm512_setzero_pd();
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41162_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64) B0r,
                                          double * __restrict __ATTR_ALIGN__(64) B0i) {

                        _mm512_store_pd(&B0r[0] ,_mm512_setzero_pd());
                        _mm512_store_pd(&B0i[0] ,_mm512_setzero_pd());
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41162_zmm8r8_u(double * __restrict  B0r,
                                          double * __restrict  B0i) {

                        _mm512_storeu_pd(&B0r[0] ,_mm512_setzero_pd());
                        _mm512_storeu_pd(&B0i[0] ,_mm512_setzero_pd());
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41162_zmm8r8(__m512d * __restrict A1r,
                                          __m512d * __restrict A1i) {

                        *A1r = _mm512_setzero_pd();
                        *A1i = _mm512_setzero_pd();
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41162_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64) A1r,
                                          double * __restrict __ATTR_ALIGN__(64) A1i) {

                        _mm512_store_pd(&A1r[0] ,_mm512_setzero_pd());
                        _mm512_store_pd(&A1i[0] ,_mm512_setzero_pd());
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41162_zmm8r8_u(double * __restrict  A1r,
                                            double * __restrict  A1i) {

                        _mm512_storeu_pd(&A1r[0] ,_mm512_setzero_pd());
                        _mm512_storeu_pd(&A1i[0] ,_mm512_setzero_pd());
               }


                      __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41162_zmm8r8(const __m512d k0a,
                                          __m512d * __restrict B1r,
                                          __m512d * __restrict B1i) {

                         const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d c0 = _mm512_set1_pd(1.8992f);
                         register __m512d k0a2,k0ah;
                         k0a2 = _mm512_mul_pd(k0a,k0a);
                         *B1r = Ir;
                         k0ah = _mm512_mul_pd(c0,k0a2);
                         *B1i = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0ah));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41162_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            double * __restrict __ATTR_ALIGN__(64) B1r,
                                            double * __restrict __ATTR_ALIGN__(64) B1i) {

                   
                         register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                         const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d c0 = _mm512_set1_pd(1.8992f);
                         register __m512d k0a2,k0ah;
                         k0a2 = _mm512_mul_pd(k0a,k0a);
                         _mm512_store_pd(&B1r[0] ,Ir);
                         k0ah = _mm512_mul_pd(c0,k0a2);
                         _mm512_store_pd(&B1i[0] ,_mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0ah)));
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41162_zmm8r8_u(const double * __restrict  pk0a,
                                            double * __restrict  B1r,
                                            double * __restrict  B1i) {

                   
                         register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                         const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d c0 = _mm512_set1_pd(1.8992f);
                         register __m512d k0a2,k0ah;
                         k0a2 = _mm512_mul_pd(k0a,k0a);
                         _mm512_storeu_pd(&B1r[0] ,Ir);
                         k0ah = _mm512_mul_pd(c0,k0a2);
                         _mm512_storeu_pd(&B1i[0] ,_mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0ah)));
                }


                   /*
                          Low-frequncy approximations (k0a<0.2)
                          Cylindrical Luneberg lens (k0a<0.2).  
                          Scattering widths.
                          Formula 4.1-163
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41163_zmm8r8(const __m512d a,
                                             const __m512d k0a) {

                          const __m512d c0   = _mm512_set1_pd(0.0625f);
                          const __m512d pipi = _mm512_set1_pd( 9.869604401089358618834490999876);
                          register __m512d rcs,t0,k0a3;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          t0   = _mm512_mul_pd(pipi,a);
                          rcs  = _mm512_mul_pd(k0a3,t0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41163_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]); 
                          const __m512d c0   = _mm512_set1_pd(0.0625f);
                          const __m512d pipi = _mm512_set1_pd( 9.869604401089358618834490999876);
                          register __m512d rcs,t0,k0a3;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          t0   = _mm512_mul_pd(pipi,a);
                          rcs  = _mm512_mul_pd(k0a3,t0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41163_zmm8r8_u(const double * __restrict  pa,
                                             const double * __restrict  pk0a) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]); 
                          const __m512d c0   = _mm512_set1_pd(0.0625f);
                          const __m512d pipi = _mm512_set1_pd( 9.869604401089358618834490999876);
                          register __m512d rcs,t0,k0a3;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          t0   = _mm512_mul_pd(pipi,a);
                          rcs  = _mm512_mul_pd(k0a3,t0);
                          return (rcs);
                 }


                    /*
                          Low-frequncy approximations (k0a<0.2)
                          Cylindrical Luneberg lens (k0a<0.2).  
                          Scattering widths.
                          Formula 4.1-164
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41164_zmm8r8(const __m512d a,
                                             const __m512d k0a,
                                             const __m512d phi) {

                          const __m512d c0   = _mm512_set1_pd(0.03607f);
                          const __m512d pipi = _mm512_set1_pd( 9.869604401089358618834490999876);
                          register __m512d rcs,t0,cosp,k0a3;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(c0,_mm512_mul_pd(pipi,a));
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          rcs   = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,cos2p));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41164_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                               const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]); 
                          register __m512d phi = _mm512_load_pd(&pphi[0]);
                          const __m512d c0   = _mm512_set1_pd(0.03607f);
                          const __m512d pipi = _mm512_set1_pd( 9.869604401089358618834490999876);
                          register __m512d rcs,t0,cosp,k0a3;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(c0,_mm512_mul_pd(pipi,a));
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          rcs   = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,cos2p));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f41164_zmm8r8_u(const double * __restrict  pa,
                                               const double * __restrict  pk0a,
                                               const double * __restrict  pphi) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]); 
                          register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                          const __m512d c0   = _mm512_set1_pd(0.03607f);
                          const __m512d pipi = _mm512_set1_pd( 9.869604401089358618834490999876);
                          register __m512d rcs,t0,cosp,k0a3;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(c0,_mm512_mul_pd(pipi,a));
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          rcs   = _mm512_mul_pd(t0,_mm512_mul_pd(k0a3,cos2p));
                          return (rcs);
                 }


                  /*

                      Cylindrical Eaton-Lippman Lens, (k0a<0.2)
                      Formulae 4.1-165
                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41165_zmm8r8(const __m512d k0a,
                                          __m512d * __restrict A0r,
                                          __m512d * __restrict A0i) {

                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a2;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        *A0r = Ir;
                        *A0i = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a2));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41165_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            double * __restrict __ATTR_ALIGN__(64) A0r,
                                            double * __restrict __ATTR_ALIGN__(64) A0i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a2;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        _mm512_store_pd(&A0r[0], Ir);
                        _mm512_store_pd(&A0i[0], _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a2)));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A0_f41165_zmm8r8_u(const double * __restrict  pk0a,
                                            double * __restrict  A0r,
                                            double * __restrict  A0i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a2;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        _mm512_storeu_pd(&A0r[0], Ir);
                        _mm512_storeu_pd(&A0i[0], _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,k0a2)));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41165_zmm8r8(__m512d * __restrict A0r,
                                          __m512d * __restrict A0i) {

                        *A0r = Ir;
                        *A0i = Ir;
                }


                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41165_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64) A0r,
                                            double * __restrict __ATTR_ALIGN__(64) A0i) {

                        _mm512_store_pd(&A0r[0], Ir);
                        _mm512_store_pd(&A0i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A1_f41165_zmm8r8_u(double * __restrict  A0r,
                                            double * __restrict  A0i) {

                        _mm512_storeu_pd(&A0r[0], Ir);
                        _mm512_storeu_pd(&A0i[0], Ir);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41165_zmm8r8(__m512d * __restrict B0r,
                                          __m512d * __restrict B0i) {

                        *B0r = Ir;
                        *B0i = Ir;
                }


                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41165_zmm8r8_a(double * __restrict __ATTR_ALIGN__(64) B0r,
                                            double * __restrict __ATTR_ALIGN__(64) B0i) {

                        _mm512_store_pd(&B0r[0], Ir);
                        _mm512_store_pd(&B0i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B0_f41165_zmm8r8_u(double * __restrict  B0r,
                                            double * __restrict  B0i) {

                        _mm512_storeu_pd(&B0r[0], Ir);
                        _mm512_storeu_pd(&B0i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41165_zmm8r8(const __m512d k0a,
                                          __m512d * __restrict B1r,
                                          __m512d * __restrict B1i) {

                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d c0  = _mm512_set1_pd(0.43616f);
                        register __m512d k0a2,t0;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        t0   = _mm512_mul_pd(c0,k0a2);
                        *B1r = Ir;
                        *B1i = _mm512_mul_pd(Ii,_mm512_mul_pd(pi4,t0));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41165_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                          double * __restrict __ATTR_ALIGN__(64) B1r,
                                          double * __restrict __ATTR_ALIGN__(64) B1i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d c0  = _mm512_set1_pd(0.43616f);
                        register __m512d k0a2,t0;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        t0   = _mm512_mul_pd(c0,k0a2);
                        _mm512_store_pd(&B1r[0] ,Ir);
                        _mm512_store_pd(&B1i[0] ,_mm512_mul_pd(Ii,_mm512_mul_pd(pi4,t0)));
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B1_f41165_zmm8r8_u(const double * __restrict  pk0a,
                                          double * __restrict  B1r,
                                          double * __restrict  B1i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d c0  = _mm512_set1_pd(0.43616f);
                        register __m512d k0a2,t0;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        t0   = _mm512_mul_pd(c0,k0a2);
                        _mm512_storeu_pd(&B1r[0] ,Ir);
                        _mm512_storeu_pd(&B1i[0] ,_mm512_mul_pd(Ii,_mm512_mul_pd(pi4,t0)));
                }


                  /*

                       Cylindrical Eaton-Lippman Lens, (k0a<0.2) 
                       Scattering widths.
                       Formula: 1.4-166
                   */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f14166_zmm8r8(const __m512d a,
                                             const __m512d k0a) {

                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d pip = _mm512_set1_pd(9.869604401089358618834490999876);
                          register __m512d a4,k0a3,rcs;
                          a4   = _mm512_mul_pd(a,qtr);
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          rcs  = _mm512_mul_pd(k0a3,_mm512_mul_pd(pip,a4));
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f14166_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d pip = _mm512_set1_pd(9.869604401089358618834490999876);
                          register __m512d a4,k0a3,rcs;
                          a4   = _mm512_mul_pd(a,qtr);
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          rcs  = _mm512_mul_pd(k0a3,_mm512_mul_pd(pip,a4));
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f14166_zmm8r8_u(const double * __restrict pa,
                                               const double * __restrict  pk0a) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d qtr = _mm512_set1_pd(0.25f);
                          const __m512d pip = _mm512_set1_pd(9.869604401089358618834490999876);
                          register __m512d a4,k0a3,rcs;
                          a4   = _mm512_mul_pd(a,qtr);
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          rcs  = _mm512_mul_pd(k0a3,_mm512_mul_pd(pip,a4));
                          return (rcs);
               }


                  /*

                       Cylindrical Eaton-Lippman Lens, (k0a<0.2) 
                       Scattering widths.
                       Formula: 1.4-167
                   */
                 
 
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f14167_zmm8r8(const __m512d a,
                                             const __m512d k0a,
                                             const __m512d phi) {

                          const __m512d c0  = _mm512_set1_pd(0.19024f);
                          const __m512d pip = _mm512_set1_pd(9.869604401089358618834490999876);
                          register __m512d cosp,cos2p,k0a3,rcs,t1;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          t1   = _mm512_mul_pd(_mm512_mul_pd(c0,pip),
                                               _mm512_mul_pd(a,k0a3));
                          cosp = xcosf(phi);
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          rcs   = _mm512_mul_pd(t1,cos2p);         
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f14167_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                               const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d phi = _mm512_load_pd(&pphi[0]);
                          const __m512d c0  = _mm512_set1_pd(0.19024f);
                          const __m512d pip = _mm512_set1_pd(9.869604401089358618834490999876);
                          register __m512d cosp,cos2p,k0a3,rcs,t1;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          t1   = _mm512_mul_pd(_mm512_mul_pd(c0,pip),
                                               _mm512_mul_pd(a,k0a3));
                          cosp = xcosf(phi);
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          rcs   = _mm512_mul_pd(t1,cos2p);         
                          return (rcs);
               }

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f14167_zmm8r8_u(const double * __restrict  pa,
                                               const double * __restrict  pk0a,
                                               const double * __restrict  pphi) {

                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d phi = _mm512_loadu_pd(&pphi[0]);
                          const __m512d c0  = _mm512_set1_pd(0.19024f);
                          const __m512d pip = _mm512_set1_pd(9.869604401089358618834490999876);
                          register __m512d cosp,cos2p,k0a3,rcs,t1;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          t1   = _mm512_mul_pd(_mm512_mul_pd(c0,pip),
                                               _mm512_mul_pd(a,k0a3));
                          cosp = xcosf(phi);
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          rcs   = _mm512_mul_pd(t1,cos2p);         
                          return (rcs);
               }


                  /*

                        Infinitely long cylinder.
                        Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                        TM-incident E-field.
                        Formula 4.2-48
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4248_zmm8r8(const __m512d E0r,
                                         const __m512d E0i,
                                         const __m512d psi,
                                         const __m512d phi,
                                         const __m512d k0,
                                         const __m512d z,
                                         const __m512d r,
                                         const __m512d a0,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui,
                                         __m512d * __restrict Ezr,
                                         __m512d * __restrict Ezi) {

                        const __m512d spi2 = _mm512_set1_pd(0.886226925452758013649083741671);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d _2  = _mm512_set1_pd(2.0f);
                        register __m512d k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512d epsrp1,epsip1,epsrm1,epsim1;
                        register __m512d murp1,muip1,murm1,muim1,k0a02;
                        register __m512d mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512d ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512d frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_pd(k0,r);
                        k0z  = _mm512_mul_pd(k0,z);
                        k0a0 = _mm512_mul_pd(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_pd(hlf,_mm512_mul_pd(k0a0,k0a0));
                        epsrp1 = _mm512_add_pd(epsr,_1);
                        ear    = Ir;
                        epsip1 = _mm512_add_pd(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        scosps = _mm512_sqrt_pd(cosps);
                        epsim1 = _mm512_sub_pd(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        cos2ps= _mm512_mul_pd(cosps,cosps);
                        murm1  = _mm512_sub_pd(mur,_1);
                        cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_pd(mui,_1);
                        sin2ps= _mm512_mul_pd(sinps,sinps);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_pd(k0r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_pd(E0r,scosps);
                        fraci = _mm512_mul_pd(E0i,scosps);
                        cmul_zmm8r8(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_pd(spi2,_mm512_div_pd(frer,t1));
                        cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_pd(spi2,_mm512_div_pd(frei,t1));
                        cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_pd(epsrm1,cos2ps);
                        t0i = _mm512_mul_pd(epsim1,cos2ps);
                        numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                        cdiv_zmm8r8(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                        t1i = _mm512_mul_pd(_2,_mm512_mul_pd(div2i,cosp));
                        t2r = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0r,t1r));
                        t2i = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0i,t1i));
                        cmul_zmm8r8(div1r,div1i,t2r,t2i,*Ezr,*Ezi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4248_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pE0r,
                                           const double * __restrict __ATTR_ALIGN__(64) pE0i,
                                           const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                           const double * __restrict __ATTR_ALIGN__(64) pphi,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pz,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pa0,
                                           const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                           const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                           const double * __restrict __ATTR_ALIGN__(64) pmur,
                                           const double * __restrict __ATTR_ALIGN__(64) pmui,
                                           double * __restrict __ATTR_ALIGN__(64) Ezr,
                                           double * __restrict __ATTR_ALIGN__(64) Ezi) {

                        register __m512d E0r = _mm512_load_pd(&pE0r[0]);
                        register __m512d E0i = _mm512_load_pd(&pE0i[0]);
                        register __m512d psi = _mm512_load_pd(&ppsi[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d z   = _mm512_load_pd(&pz[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d a0  = _mm512_load_pd(&pa0[0]);
                        register __m512d epsr= _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi= _mm512_load_pd(&pepsi[0]);
                        register __m512d pmur= _mm512_load_pd(&pmur[0]);
                        register __m512d pmui= _mm512_load_pd(&pmui[0]);
                        const __m512d spi2 = _mm512_set1_pd(0.886226925452758013649083741671);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d _2  = _mm512_set1_pd(2.0f);
                        register __m512d k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512d epsrp1,epsip1,epsrm1,epsim1;
                        register __m512d murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512d mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512d ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512d frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_pd(k0,r);
                        k0z  = _mm512_mul_pd(k0,z);
                        k0a0 = _mm512_mul_pd(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_pd(hlf,_mm512_mul_pd(k0a0,k0a0));
                        epsrp1 = _mm512_add_pd(epsr,_1);
                        ear    = Ir;
                        epsip1 = _mm512_add_pd(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        scosps = _mm512_sqrt_pd(cosps);
                        epsim1 = _mm512_sub_pd(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        cos2ps= _mm512_mul_pd(cosps,cosps);
                        murm1  = _mm512_sub_pd(mur,_1);
                        cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_pd(mui,_1);
                        sin2ps= _mm512_mul_pd(sinps,sinps);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_pd(k0r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_pd(E0r,scosps);
                        fraci = _mm512_mul_pd(E0i,scosps);
                        cmul_zmm8r8(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_pd(spi2,_mm512_div_pd(frer,t1));
                        cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_pd(spi2,_mm512_div_pd(frei,t1));
                        cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_pd(epsrm1,cos2ps);
                        t0i = _mm512_mul_pd(epsim1,cos2ps);
                        numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                        cdiv_zmm8r8(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                        t1i = _mm512_mul_pd(_2,_mm512_mul_pd(div2i,cosp));
                        t2r = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0r,t1r));
                        t2i = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0i,t1i));
                        cmul_zmm8r8(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_store_pd(&Ezr[0], resr);
                        _mm512_store_pd(&Ezi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4248_zmm8r8_u(const double * __restrict  pE0r,
                                           const double * __restrict  pE0i,
                                           const double * __restrict  ppsi,
                                           const double * __restrict  pphi,
                                           const double * __restrict  pk0,
                                           const double * __restrict  pz,
                                           const double * __restrict  pr,
                                           const double * __restrict  pa0,
                                           const double * __restrict  pepsr,
                                           const double * __restrict  pepsi,
                                           const double * __restrict  pmur,
                                           const double * __restrict  pmui,
                                           double * __restrict  Ezr,
                                           double * __restrict  Ezi) {

                        register __m512d E0r = _mm512_loadu_pd(&pE0r[0]);
                        register __m512d E0i = _mm512_loadu_pd(&pE0i[0]);
                        register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d z   = _mm512_loadu_pd(&pz[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d a0  = _mm512_loadu_pd(&pa0[0]);
                        register __m512d epsr= _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi= _mm512_loadu_pd(&pepsi[0]);
                        register __m512d pmur= _mm512_loadu_pd(&pmur[0]);
                        register __m512d pmui= _mm512_loadu_pd(&pmui[0]);
                        const __m512d spi2 = _mm512_set1_pd(0.886226925452758013649083741671);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d _2  = _mm512_set1_pd(2.0f);
                        register __m512d k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512d epsrp1,epsip1,epsrm1,epsim1;
                        register __m512d murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512d mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512d ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512d frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_pd(k0,r);
                        k0z  = _mm512_mul_pd(k0,z);
                        k0a0 = _mm512_mul_pd(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_pd(hlf,_mm512_mul_pd(k0a0,k0a0));
                        epsrp1 = _mm512_add_pd(epsr,_1);
                        ear    = Ir;
                        epsip1 = _mm512_add_pd(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        scosps = _mm512_sqrt_pd(cosps);
                        epsim1 = _mm512_sub_pd(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        cos2ps= _mm512_mul_pd(cosps,cosps);
                        murm1  = _mm512_sub_pd(mur,_1);
                        cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_pd(mui,_1);
                        sin2ps= _mm512_mul_pd(sinps,sinps);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_pd(k0r);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_pd(E0r,scosps);
                        fraci = _mm512_mul_pd(E0i,scosps);
                        cmul_zmm8r8(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_pd(spi2,_mm512_div_pd(frer,t1));
                        cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_pd(spi2,_mm512_div_pd(frei,t1));
                        cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_pd(epsrm1,cos2ps);
                        t0i = _mm512_mul_pd(epsim1,cos2ps);
                        numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                        cdiv_zmm8r8(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                        t1i = _mm512_mul_pd(_2,_mm512_mul_pd(div2i,cosp));
                        t2r = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0r,t1r));
                        t2i = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0i,t1i));
                        cmul_zmm8r8(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_storeu_pd(&Ezr[0], resr);
                        _mm512_storeu_pd(&Ezi[0], resi);
                }


               

                       /*

                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TM-incident H-field.
                         Formula 4.2-51
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hp_f4251_zmm8r8(const __m512d E0r,
                                         const __m512d E0i,
                                         const __m512d psi,
                                         const __m512d phi,
                                         const __m512d k0,
                                         const __m512d z,
                                         const __m512d r,
                                         const __m512d a0,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui,
                                         __m512d * __restrict Hpr,
                                         __m512d * __restrict Hpi) {

                        const __m512d e0u0 = _mm512_set1_pd(0.000007036193308495678572187302);
                        const __m512d spi2 = _mm512_set1_pd(0.886226925452758013649083741671);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d _2  = _mm512_set1_pd(2.0f);
                        register __m512d k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512d epsrp1,epsip1,epsrm1,epsim1,rat,spirat;
                        register __m512d murp1,muip1,murm1,muim1,k0a02;
                        register __m512d mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512d ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512d frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_pd(k0,r);
                        k0z  = _mm512_mul_pd(k0,z);
                        k0a0 = _mm512_mul_pd(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_pd(hlf,_mm512_mul_pd(k0a0,k0a0));
                        epsrp1 = _mm512_add_pd(epsr,_1);
                        ear    = Ir;
                        rat    = _mm512_sqrt_pd(e0u0);
                        epsip1 = _mm512_add_pd(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        scosps = _mm512_sqrt_pd(cosps);
                        epsim1 = _mm512_sub_pd(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        cos2ps= _mm512_mul_pd(cosps,cosps);
                        murm1  = _mm512_sub_pd(mur,_1);
                        cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_pd(mui,_1);
                        spirat = _mm512_mul_pd(spi2,rat);
                        sin2ps= _mm512_mul_pd(sinps,sinps);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_pd(E0r,scosps);
                        fraci = _mm512_mul_pd(E0i,scosps);
                        cmul_zmm8r8(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_pd(spirat,_mm512_div_pd(frer,t1));
                        cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_pd(spirat,_mm512_div_pd(frei,t1));
                        cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_pd(epsrm1,cos2ps);
                        t0i = _mm512_mul_pd(epsim1,cos2ps);
                        numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                        cdiv_zmm8r8(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                        t1i = _mm512_mul_pd(_2,_mm512_mul_pd(div2i,cosp));
                        t2r = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0r,t1r));
                        t2i = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0i,t1i));
                        div1r = _mm512_mul_pd(nIi,div1r);
                        div1i = _mm512_mul_pd(nIi,div1i);
                        cmul_zmm8r8(div1r,div1i,t2r,t2i,*Hpr,*Hpi);
                        
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hp_f4251_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pE0r,
                                           const double * __restrict __ATTR_ALIGN__(64) pE0i,
                                           const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                           const double * __restrict __ATTR_ALIGN__(64) pphi,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0,
                                           const double * __restrict __ATTR_ALIGN__(64) pz,
                                           const double * __restrict __ATTR_ALIGN__(64) pr,
                                           const double * __restrict __ATTR_ALIGN__(64) pa0,
                                           const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                           const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                           const double * __restrict __ATTR_ALIGN__(64) pmur,
                                           const double * __restrict __ATTR_ALIGN__(64) pmui,
                                           double * __restrict __ATTR_ALIGN__(64) Hpr,
                                           double * __restrict __ATTR_ALIGN__(64) Hpi) {


                        register __m512d E0r = _mm512_load_pd(&pE0r[0]);
                        register __m512d E0i = _mm512_load_pd(&pE0i[0]);
                        register __m512d psi = _mm512_load_pd(&ppsi[0]);
                        register __m512d k0  = _mm512_load_pd(&pk0[0]);
                        register __m512d z   = _mm512_load_pd(&pz[0]);
                        register __m512d r   = _mm512_load_pd(&pr[0]);
                        register __m512d a0  = _mm512_load_pd(&pa0[0]);
                        register __m512d epsr= _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi= _mm512_load_pd(&pepsi[0]);
                        register __m512d pmur= _mm512_load_pd(&pmur[0]);
                        register __m512d pmui= _mm512_load_pd(&pmui[0]);
                        const __m512d e0u0 = _mm512_set1_pd(0.000007036193308495678572187302);
                        const __m512d spi2 = _mm512_set1_pd(0.886226925452758013649083741671);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d _2  = _mm512_set1_pd(2.0f);
                        register __m512d k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512d epsrp1,epsip1,epsrm1,epsim1,rat,spirat;
                        register __m512d murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512d mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512d ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512d frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_pd(k0,r);
                        k0z  = _mm512_mul_pd(k0,z);
                        rat  = _mm512_div_pd(eps0,mu0);
                        k0a0 = _mm512_mul_pd(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_pd(hlf,_mm512_mul_pd(k0a0,k0a0));
                        epsrp1 = _mm512_add_pd(epsr,_1);
                        ear    = Ir;
                        rat    = _mm512_sqrt_pd(e0u0);
                        epsip1 = _mm512_add_pd(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        scosps = _mm512_sqrt_pd(cosps);
                        epsim1 = _mm512_sub_pd(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        cos2ps= _mm512_mul_pd(cosps,cosps);
                        murm1  = _mm512_sub_pd(mur,_1);
                        cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_pd(mui,_1);
                        spirat = _mm512_mul_pd(spi2,rat);
                        sin2ps= _mm512_mul_pd(sinps,sinps);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_pd(E0r,scosps);
                        fraci = _mm512_mul_pd(E0i,scosps);
                        cmul_zmm8r8(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_pd(spirat,_mm512_div_pd(frer,t1));
                        cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_pd(spirat,_mm512_div_pd(frei,t1));
                        cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_pd(epsrm1,cos2ps);
                        t0i = _mm512_mul_pd(epsim1,cos2ps);
                        numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                        cdiv_zmm8r8(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                        t1i = _mm512_mul_pd(_2,_mm512_mul_pd(div2i,cosp));
                        t2r = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0r,t1r));
                        t2i = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0i,t1i));
                        div1r = _mm512_mul_pd(nIi,div1r);
                        div1i = _mm512_mul_pd(nIi,div1i);
                        cmul_zmm8r8(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_store_pd(&Hpr[0], resr);
                        _mm512_store_pd(&Hpi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hp_f4251_zmm8r8_u(  const double * __restrict  pE0r,
                                           const double * __restrict  pE0i,
                                           const double * __restrict  ppsi,
                                           const double * __restrict pphi,
                                           const double * __restrict  pk0,
                                           const double * __restrict  pz,
                                           const double * __restrict  pr,
                                           const double * __restrict  pa0,
                                           const double * __restrict  pepsr,
                                           const double * __restrict  pepsi,
                                           const double * __restrict  pmur,
                                           const double * __restrict  pmui,
                                           double * __restrict  Hpr,
                                           double * __restrict  Hpi) {


                        register __m512d E0r = _mm512_loadu_pd(&pE0r[0]);
                        register __m512d E0i = _mm512_loadu_pd(&pE0i[0]);
                        register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                        register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                        register __m512d z   = _mm512_loadu_pd(&pz[0]);
                        register __m512d r   = _mm512_loadu_pd(&pr[0]);
                        register __m512d a0  = _mm512_loadu_pd(&pa0[0]);
                        register __m512d epsr= _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi= _mm512_loadu_pd(&pepsi[0]);
                        register __m512d pmur= _mm512_loadu_pd(&pmur[0]);
                        register __m512d pmui= _mm512_loadu_pd(&pmui[0]);
                        const __m512d e0u0 = _mm512_set1_pd(0.000007036193308495678572187302);
                        const __m512d spi2 = _mm512_set1_pd(0.886226925452758013649083741671);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d _2  = _mm512_set1_pd(2.0f);
                        register __m512d k0r,k0z,k0a0,cosp,cosps,cos2ps,sinps,sin2ps;
                        register __m512d epsrp1,epsip1,epsrm1,epsim1,rat,spirat;
                        register __m512d murp1,muip1,murm1,muim1,k0a02,resr,resi;
                        register __m512d mul1r,mul1i,mul2r,mul2i,mul3r,mul3i,t0r,t0i,t2r,t2i;
                        register __m512d ear,eai,t0,t1,cer,cei,fracr,fraci,scosps;
                        register __m512d frer,frei,div1r,div1i,div2r,div2i,numr,numi,t1r,t1i;
                        k0r  = _mm512_mul_pd(k0,r);
                        k0z  = _mm512_mul_pd(k0,z);
                        rat  = _mm512_div_pd(eps0,mu0);
                        k0a0 = _mm512_mul_pd(k0,a0);
                        cosp = xcosf(phi);
                        k0a02= _mm512_mul_pd(hlf,_mm512_mul_pd(k0a0,k0a0));
                        epsrp1 = _mm512_add_pd(epsr,_1);
                        ear    = Ir;
                        rat    = _mm512_sqrt_pd(e0u0);
                        epsip1 = _mm512_add_pd(epsi,_1)
                        cosps= xcosf(psi);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        scosps = _mm512_sqrt_pd(cosps);
                        epsim1 = _mm512_sub_pd(epsi,_1)
                        sinps= xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        cos2ps= _mm512_mul_pd(cosps,cosps);
                        murm1  = _mm512_sub_pd(mur,_1);
                        cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                        muim1  = _mm512_sub_pd(mui,_1);
                        spirat = _mm512_mul_pd(spi2,rat);
                        sin2ps= _mm512_mul_pd(sinps,sinps);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                        eai    = t0;
                        t1     = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        fracr = _mm512_mul_pd(E0r,scosps);
                        fraci = _mm512_mul_pd(E0i,scosps);
                        cmul_zmm8r8(fracr,fraci,cer,cei,&frer,&frei);
                        div1r = _mm512_mul_pd(spirat,_mm512_div_pd(frer,t1));
                        cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                        div1i = _mm512_mul_pd(spirat,_mm512_div_pd(frei,t1));
                        cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                        t0r = _mm512_mul_pd(epsrm1,cos2ps);
                        t0i = _mm512_mul_pd(epsim1,cos2ps);
                        numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                        numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                        cdiv_zmm8r8(numr,numi,mul3r,mul3i,&div2r,&div2i);
                        t1r = _mm512_mul_pd(_2,_mm512_mul_pd(div2r,cosp));
                        t1i = _mm512_mul_pd(_2,_mm512_mul_pd(div2i,cosp));
                        t2r = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0r,t1r));
                        t2i = _mm512_mul_pd(k0a02,_mm512_sub_pd(t0i,t1i));
                        div1r = _mm512_mul_pd(nIi,div1r);
                        div1i = _mm512_mul_pd(nIi,div1i);
                        cmul_zmm8r8(div1r,div1i,t2r,t2i,&resr,&resi);
                        _mm512_storeu_pd(&Hpr[0], resr);
                        _mm512_storeu_pd(&Hpi[0], resi);
                }


               
                     /*

                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TM-incident E-field.
                         Formula 4.2-49
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4249_zmm8r8(const __m512d E0r,
                                          const __m512d E0i,
                                          const __m512d k0z,
                                          const __m512d k0r,
                                          const __m512d k0a0,
                                          const __m512d psi,
                                          const __m512d phi,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          __m512d * __restrict Ephr,
                                          __m512d * __restrict Ephi) {

                        const __m512d spi2 = _mm512_set1_pd(2.506628274631000502415765284811);
                        const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i;
                        register __m512d cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512d emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_pd(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        t0    = _mm512_mul_pd(k0r,cosp);
                        ear   = Ir;
                        cmul_zmm8r8(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_pd(t0);
                        emum1r = _mm512_sub_pd(emum1r,_1);
                        emum1i = _mm512_sub_pd(emum1i,_1);
                        epsp1r = _mm512_add_pd(epsr,_1);
                        epsp1i = _mm512_add_pd(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosp,pi4));
                        sinpsp = _mm512_mul_pd(sinps,sinp);
                        eai    = t0;
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cmul_zmm8r8(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_pd(fracr,den);
                        t0r   = _mm512_mul_pd(spi2,_mm512_mul_pd(fracr,k0a02));
                        fraci = _mm512_div_pd(fraci,den);
                        t0i   = _mm512_mul_pd(spi2,_mm512_mul_pd(fraci,k0a02));
                        cmul_zmm8r8(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm8r8(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_pd(sinpsp,t1r);
                        t1i = _mm512_mul_pd(sinpsp,t1i);
                        cmul_zmm8r8(fracr,fraci,t1r,t1i,*Ephr,*Ephi);
               }
                   


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4249_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pE0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pE0i,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          double * __restrict __ATTR_ALIGN__(64) Ephr,
                                          double * __restrict __ATTR_ALIGN__(64) Ephi) {

                        register __m512d E0r   = _mm512_load_pd(&pE0r[0]);
                        register __m512d E0i   = _mm512_load_pd(&pE0i[0]);
                        register __m512d k0z   = _mm512_load_pd(&pk0z[0]);
                        register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                        register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                        register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                        register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_load_pd(&pmur[0]);
                        register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                        const __m512d spi2 = _mm512_set1_pd(2.506628274631000502415765284811);
                        const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,resr,resi;
                        register __m512d cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512d emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_pd(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        t0    = _mm512_mul_pd(k0r,cosp);
                        ear   = Ir;
                        cmul_zmm8r8(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_pd(t0);
                        emum1r = _mm512_sub_pd(emum1r,_1);
                        emum1i = _mm512_sub_pd(emum1i,_1);
                        epsp1r = _mm512_add_pd(epsr,_1);
                        epsp1i = _mm512_add_pd(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosp,pi4));
                        sinpsp = _mm512_mul_pd(sinps,sinp);
                        eai    = t0;
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cmul_zmm8r8(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_pd(fracr,den);
                        t0r   = _mm512_mul_pd(spi2,_mm512_mul_pd(fracr,k0a02));
                        fraci = _mm512_div_pd(fraci,den);
                        t0i   = _mm512_mul_pd(spi2,_mm512_mul_pd(fraci,k0a02));
                        cmul_zmm8r8(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm8r8(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_pd(sinpsp,t1r);
                        t1i = _mm512_mul_pd(sinpsp,t1i);
                        cmul_zmm8r8(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_store_pd(&Ephr[0],resr);
                        _mm512_store_pd(&Ephi[0],resi);
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4249_zmm8r8_u(const double * __restrict pE0r,
                                          const double * __restrict  pE0i,
                                          const double * __restrict pk0z,
                                          const double * __restrict  pk0r,
                                          const double * __restrict  pk0a0,
                                          const double * __restrict  ppsi,
                                          const double * __restrict  pphi,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          double * __restrict Ephr,
                                          double * __restrict  Ephi) {

                        register __m512d E0r   = _mm512_loadu_pd(&pE0r[0]);
                        register __m512d E0i   = _mm512_loadu_pd(&pE0i[0]);
                        register __m512d k0z   = _mm512_loadu_pd(&pk0z[0]);
                        register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                        register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                        register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                        register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                        const __m512d spi2 = _mm512_set1_pd(2.506628274631000502415765284811);
                        const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,resr,resi;
                        register __m512d cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512d emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_pd(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        t0    = _mm512_mul_pd(k0r,cosp);
                        ear   = Ir;
                        cmul_zmm8r8(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_pd(t0);
                        emum1r = _mm512_sub_pd(emum1r,_1);
                        emum1i = _mm512_sub_pd(emum1i,_1);
                        epsp1r = _mm512_add_pd(epsr,_1);
                        epsp1i = _mm512_add_pd(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosp,pi4));
                        sinpsp = _mm512_mul_pd(sinps,sinp);
                        eai    = t0;
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cmul_zmm8r8(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_pd(fracr,den);
                        t0r   = _mm512_mul_pd(spi2,_mm512_mul_pd(fracr,k0a02));
                        fraci = _mm512_div_pd(fraci,den);
                        t0i   = _mm512_mul_pd(spi2,_mm512_mul_pd(fraci,k0a02));
                        cmul_zmm8r8(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm8r8(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_pd(sinpsp,t1r);
                        t1i = _mm512_mul_pd(sinpsp,t1i);
                        cmul_zmm8r8(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_storeu_pd(&Ephr[0],resr);
                        _mm512_storeu_pd(&Ephi[0],resi);
               }



                 /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TM-incident H-field.
                         Formula 4.2-50

                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4250_zmm8r8( const __m512d E0r,
                                          const __m512d E0i,
                                          const __m512d k0z,
                                          const __m512d k0r,
                                          const __m512d k0a0,
                                          const __m512d psi,
                                          const __m512d phi,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          __m512d * __restrict Hzr,
                                          __m512d * __restrict Hzi) {

                        const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                        const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,scosp;
                        register __m512d cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli;
                        register __m512d emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_pd(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        ear   = Ir;
                        cmul_zmm8r8(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_pd(k0r);
                        emum1r = _mm512_sub_pd(emum1r,_1);
                        emum1i = _mm512_sub_pd(emum1i,_1);
                        epsp1r = _mm512_add_pd(epsr,_1);
                        epsp1i = _mm512_add_pd(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosp,pi4));
                        scosp  = _mm512_sqrt_pd(cosp);
                        sinpsp = _mm512_mul_pd(sinps,sinp);
                        eai    = t0;
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cer = _mm512_mul_pd(scosp,cer);
                        cei = _mm512_mul_pd(scosp,cei);
                        cmul_zmm8r8(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_pd(fracr,den);
                        t0r   = _mm512_mul_pd(e0u0,_mm512_mul_pd(fracr,k0a02));
                        fraci = _mm512_div_pd(fraci,den);
                        t0i   = _mm512_mul_pd(e0u0,_mm512_mul_pd(fraci,k0a02));
                        cmul_zmm8r8(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm8r8(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_pd(sinpsp,t1r);
                        t1i = _mm512_mul_pd(sinpsp,t1i);
                        cmul_zmm8r8(fracr,fraci,t1r,t1i,*Hzr,*Hzi);
               }


                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4250_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pE0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pE0i,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          double * __restrict __ATTR_ALIGN__(64) Hzr,
                                          double * __restrict __ATTR_ALIGN__(64) Hzi ) {

                        register __m512d E0r   = _mm512_load_pd(&pE0r[0]);
                        register __m512d E0i   = _mm512_load_pd(&pE0i[0]);
                        register __m512d k0z   = _mm512_load_pd(&pk0z[0]);
                        register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                        register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                        register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                        register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_load_pd(&pmur[0]);
                        register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                        const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                        const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,scosp;
                        register __m512d cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli,resr,resi;
                        register __m512d emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_pd(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        ear   = Ir;
                        cmul_zmm8r8(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_pd(k0r);
                        emum1r = _mm512_sub_pd(emum1r,_1);
                        emum1i = _mm512_sub_pd(emum1i,_1);
                        epsp1r = _mm512_add_pd(epsr,_1);
                        epsp1i = _mm512_add_pd(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosp,pi4));
                        scosp  = _mm512_sqrt_pd(cosp);
                        sinpsp = _mm512_mul_pd(sinps,sinp);
                        eai    = t0;
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cer = _mm512_mul_pd(scosp,cer);
                        cei = _mm512_mul_pd(scosp,cei);
                        cmul_zmm8r8(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_pd(fracr,den);
                        t0r   = _mm512_mul_pd(e0u0,_mm512_mul_pd(fracr,k0a02));
                        fraci = _mm512_div_pd(fraci,den);
                        t0i   = _mm512_mul_pd(e0u0,_mm512_mul_pd(fraci,k0a02));
                        cmul_zmm8r8(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm8r8(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_pd(sinpsp,t1r);
                        t1i = _mm512_mul_pd(sinpsp,t1i);
                        cmul_zmm8r8(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_store_pd(&Hzr[0], resr);
                        _mm512_store_pd(&Hzi[0], resi);
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4250_zmm8r8_u(const double * __restrict  pE0r,
                                          const double * __restrict pE0i,
                                          const double * __restrict  pk0z,
                                          const double * __restrict  pk0r,
                                          const double * __restrict  pk0a0,
                                          const double * __restrict  ppsi,
                                          const double * __restrict  pphi,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          double * __restrict  Hzr,
                                          double * __restrict  Hzi ) {

                        register __m512d E0r   = _mm512_loadu_pd(&pE0r[0]);
                        register __m512d E0i   = _mm512_loadu_pd(&pE0i[0]);
                        register __m512d k0z   = _mm512_loadu_pd(&pk0z[0]);
                        register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                        register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                        register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                        register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                        const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                        const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d ear,eai,cer,cei,den,t0,fracr,fraci,t0r,t0i,scosp;
                        register __m512d cosp,sinps,sinp,k0a02,cosps,sinpsp,mulr,muli,resr,resi;
                        register __m512d emum1r,emum1i,epsp1r,epsp1i,murp1,muip1,t1r,t1i;
                        k0a02 = _mm512_mul_pd(k0a0,k0a0);
                        cosp  = xcosf(phi);
                        ear   = Ir;
                        cmul_zmm8r8(epsr,epsi,mur,mui,&emum1r,&emum1i);
                        den   = _mm512_sqrt_pd(k0r);
                        emum1r = _mm512_sub_pd(emum1r,_1);
                        emum1i = _mm512_sub_pd(emum1i,_1);
                        epsp1r = _mm512_add_pd(epsr,_1);
                        epsp1i = _mm512_add_pd(epsi,_1);
                        sinps  = xsinf(psi);
                        murp1  = _mm512_add_pd(mur,_1);
                        muip1  = _mm512_add_pd(mui,_1);
                        sinp   = xsinf(psi);
                        cosps  = xcosf(psi);
                        t0     = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosp,pi4));
                        scosp  = _mm512_sqrt_pd(cosp);
                        sinpsp = _mm512_mul_pd(sinps,sinp);
                        eai    = t0;
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        cer = _mm512_mul_pd(scosp,cer);
                        cei = _mm512_mul_pd(scosp,cei);
                        cmul_zmm8r8(E0r,E0i,cer,cei,&fracr,&fraci);
                        fracr = _mm512_div_pd(fracr,den);
                        t0r   = _mm512_mul_pd(e0u0,_mm512_mul_pd(fracr,k0a02));
                        fraci = _mm512_div_pd(fraci,den);
                        t0i   = _mm512_mul_pd(e0u0,_mm512_mul_pd(fraci,k0a02));
                        cmul_zmm8r8(epsp1r,epsp1i,murp1,muip1,&mulr,&muli);
                        cdiv_zmm8r8(emum1r,emum1i,mulr,muli,&t1r,&t1i);
                        t1r = _mm512_mul_pd(sinpsp,t1r);
                        t1i = _mm512_mul_pd(sinpsp,t1i);
                        cmul_zmm8r8(fracr,fraci,t1r,t1i,&resr,&resi);
                        _mm512_storeu_pd(&Hzr[0], resr);
                        _mm512_storeu_pd(&Hzi[0], resi);
               }


                 /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident H-field.
                         Formula 4.2-52

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4252_zmm8r8(const __m512d H0r,
                                         const __m512d H0i,
                                         const __m512d psi,
                                         const __m512d phi,
                                         const __m512d k0r,
                                         const __m512d k0z,
                                         const __m512d k0a0,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui,
                                         __m512d * __restrict Hzr,
                                         __m512d * __restrict Hzi) {

                         const __m512d spi2 = _mm512_set1_pd(1.253314137315500251207882642406);
                         const __m512d _2   = _mm512_set1_pd(2.0f);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         register __m512d scosps,sinps,cosps,k0a02,sk0r;
                         register __m512d cos2ps,sin2ps,cosp,mul1r,mul1i;
                         register __m512d mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512d murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512d murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(hlf,_mm512_mul_pd(k0a,k0a));
                         scosps= _mm512_sqrt_pd(cosps);
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         t0    = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_pd(k0r);
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_pd(mur,_1);
                         fracr = _mm512_mul_pd(H0r,scosps);
                         sin2ps= _mm512_mul_pd(sinps,sinps);
                         muim1 = _mm512_sub_pd(mui,_1);
                         fraci = _mm512_mul_pd(H0i,scosps);
                         epsrm1= _mm512_sub_pd(epsr,_1);
                         epsim1= _mm512_sub_pd(epsi,_1);
                         murp1 = _mm512_add_pd(mur,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_pd(murm1,cos2ps);
                         t0r = _mm512_mul_pd(spi2,_mm512_div_pd(t0r,sk0r));
                         mucsi = _mm512_mul_pd(muim1,cos2ps);
                         t0i = _mm512_mul_pd(spi2,_mm512_div_pd(t0i,sk0r));
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_pd(t0r,k0a02);
                         t0i = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_pd(_2,_mm512_mul_pd(t1r,cosp));
                         t1i = _mm512_mul_pd(_2,_mm512_mul_pd(t1i,cosp));
                         mucsr = _mm512_sub_pd(mucsr,t1r);
                         mucsi = _mm512_sub_pd(mucsi,t1i);
                         cmul_zmm8r8(t0r,t0i,mucsr,mucsi,*Hzr,*Hzi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4252_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          double * __restrict __ATTR_ALIGN__(64) Hzr,
                                          double * __restrict __ATTR_ALIGN__(64) Hzi) {

                         register __m512d H0r   = _mm512_load_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_load_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_load_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d spi2 = _mm512_set1_pd(1.253314137315500251207882642406);
                         const __m512d _2   = _mm512_set1_pd(2.0f);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         register __m512d scosps,sinps,cosps,k0a02,sk0r;
                         register __m512d cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512d mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512d murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512d murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(hlf,_mm512_mul_pd(k0a,k0a));
                         scosps= _mm512_sqrt_pd(cosps);
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         t0    = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_pd(k0r);
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_pd(mur,_1);
                         fracr = _mm512_mul_pd(H0r,scosps);
                         sin2ps= _mm512_mul_pd(sinps,sinps);
                         muim1 = _mm512_sub_pd(mui,_1);
                         fraci = _mm512_mul_pd(H0i,scosps);
                         epsrm1= _mm512_sub_pd(epsr,_1);
                         epsim1= _mm512_sub_pd(epsi,_1);
                         murp1 = _mm512_add_pd(mur,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_pd(murm1,cos2ps);
                         t0r = _mm512_mul_pd(spi2,_mm512_div_pd(t0r,sk0r));
                         mucsi = _mm512_mul_pd(muim1,cos2ps);
                         t0i = _mm512_mul_pd(spi2,_mm512_div_pd(t0i,sk0r));
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_pd(t0r,k0a02);
                         t0i = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_pd(_2,_mm512_mul_pd(t1r,cosp));
                         t1i = _mm512_mul_pd(_2,_mm512_mul_pd(t1i,cosp));
                         mucsr = _mm512_sub_pd(mucsr,t1r);
                         mucsi = _mm512_sub_pd(mucsi,t1i);
                         cmul_zmm8r8(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_store_pd(&Hzr[0], resr);
                         _mm512_store_pd(&Hzi[0], resi);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hz_f4252_zmm8r8_u(const double * __restrict  pH0r,
                                          const double * __restrict  pH0i,
                                          const double * __restrict  pk0z,
                                          const double * __restrict  pk0r,
                                          const double * __restrict  pk0a0,
                                          const double * __restrict  ppsi,
                                          const double * __restrict  pphi,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict pmur,
                                          const double * __restrict  pmui,
                                          double * __restrict  Hzr,
                                          double * __restrict Hzi) {

                         register __m512d H0r   = _mm512_loadu_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_loadu_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_loadu_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d spi2 = _mm512_set1_pd(1.253314137315500251207882642406);
                         const __m512d _2   = _mm512_set1_pd(2.0f);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         register __m512d scosps,sinps,cosps,k0a02,sk0r;
                         register __m512d cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512d mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512d murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512d murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(hlf,_mm512_mul_pd(k0a,k0a));
                         scosps= _mm512_sqrt_pd(cosps);
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         t0    = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_pd(k0r);
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_pd(mur,_1);
                         fracr = _mm512_mul_pd(H0r,scosps);
                         sin2ps= _mm512_mul_pd(sinps,sinps);
                         muim1 = _mm512_sub_pd(mui,_1);
                         fraci = _mm512_mul_pd(H0i,scosps);
                         epsrm1= _mm512_sub_pd(epsr,_1);
                         epsim1= _mm512_sub_pd(epsi,_1);
                         murp1 = _mm512_add_pd(mur,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_pd(murm1,cos2ps);
                         t0r = _mm512_mul_pd(spi2,_mm512_div_pd(t0r,sk0r));
                         mucsi = _mm512_mul_pd(muim1,cos2ps);
                         t0i = _mm512_mul_pd(spi2,_mm512_div_pd(t0i,sk0r));
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_pd(t0r,k0a02);
                         t0i = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_pd(_2,_mm512_mul_pd(t1r,cosp));
                         t1i = _mm512_mul_pd(_2,_mm512_mul_pd(t1i,cosp));
                         mucsr = _mm512_sub_pd(mucsr,t1r);
                         mucsi = _mm512_sub_pd(mucsi,t1i);
                         cmul_zmm8r8(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_storeu_pd(&Hzr[0], resr);
                         _mm512_storeu_pd(&Hzi[0], resi);
                }


                    /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident E-field.
                         Formula 4.2-55

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4255_zmm8r8(const __m512d H0r,
                                         const __m512d H0i,
                                         const __m512d psi,
                                         const __m512d phi,
                                         const __m512d k0r,
                                         const __m512d k0z,
                                         const __m512d k0a0,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui,
                                         __m512d * __restrict Epr,
                                         __m512d * __restrict Epi) {

                         const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                         const __m512d _2   = _mm512_set1_pd(2.0f);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         register __m512d sinps,cosps,k0a02,sk0r;
                         register __m512d cos2ps,sin2ps,cosp,mul1r,mul1i;
                         register __m512d mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512d murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512d murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(hlf,_mm512_mul_pd(k0a,k0a));
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         t0    = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_pd(mur,_1);
                         fracr = H0r;
                         sin2ps= _mm512_mul_pd(sinps,sinps);
                         muim1 = _mm512_sub_pd(mui,_1);
                         fraci = H0i;
                         epsrm1= _mm512_sub_pd(epsr,_1);
                         epsim1= _mm512_sub_pd(epsi,_1);
                         murp1 = _mm512_add_pd(mur,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_pd(murm1,cos2ps);
                         t0r = _mm512_mul_pd(e0u0,_mm512_div_pd(t0r,sk0r));
                         mucsi = _mm512_mul_pd(muim1,cos2ps);
                         t0i = _mm512_mul_pd(e0u0,_mm512_div_pd(t0i,sk0r));
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_pd(t0r,k0a02);
                         t0i = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_pd(_2,_mm512_mul_pd(t1r,cosp));
                         t1i = _mm512_mul_pd(_2,_mm512_mul_pd(t1i,cosp));
                         mucsr = _mm512_sub_pd(mucsr,t1r);
                         mucsi = _mm512_sub_pd(mucsi,t1i);
                         cmul_zmm8r8(t0r,t0i,mucsr,mucsi,*Epr,*Epi);
                }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4255_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          double * __restrict __ATTR_ALIGN__(64) Epr,
                                          double * __restrict __ATTR_ALIGN__(64) Epi) {

                         register __m512d H0r   = _mm512_load_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_load_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_load_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                         const __m512d _2   = _mm512_set1_pd(2.0f);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         register __m512d sinps,cosps,k0a02,sk0r;
                         register __m512d cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512d mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512d murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512d murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(hlf,_mm512_mul_pd(k0a,k0a));
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         t0    = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_pd(mur,_1);
                         fracr = H0r;
                         sin2ps= _mm512_mul_pd(sinps,sinps);
                         muim1 = _mm512_sub_pd(mui,_1);
                         fraci = H0i;
                         epsrm1= _mm512_sub_pd(epsr,_1);
                         epsim1= _mm512_sub_pd(epsi,_1);
                         murp1 = _mm512_add_pd(mur,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_pd(murm1,cos2ps);
                         t0r = _mm512_mul_pd(e0u0,_mm512_div_pd(t0r,sk0r));
                         mucsi = _mm512_mul_pd(muim1,cos2ps);
                         t0i = _mm512_mul_pd(e0u0,_mm512_div_pd(t0i,sk0r));
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_pd(t0r,k0a02);
                         t0i = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_pd(_2,_mm512_mul_pd(t1r,cosp));
                         t1i = _mm512_mul_pd(_2,_mm512_mul_pd(t1i,cosp));
                         mucsr = _mm512_sub_pd(mucsr,t1r);
                         mucsi = _mm512_sub_pd(mucsi,t1i);
                         cmul_zmm8r8(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_store_pd(&Epr[0], resr);
                         _mm512_store_pd(&Epi[0], resi);
                }


                      __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Eph_f4255_zmm8r8_u(const double * __restrict pH0r,
                                          const double * __restrict  pH0i,
                                          const double * __restrict  pk0z,
                                          const double * __restrict  pk0r,
                                          const double * __restrict  pk0a0,
                                          const double * __restrict  ppsi,
                                          const double * __restrict  pphi,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          double * __restrict  Epr,
                                          double * __restrict  Epi) {

                         register __m512d H0r   = _mm512_loadu_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_loadu_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_loadu_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                         const __m512d _2   = _mm512_set1_pd(2.0f);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         register __m512d sinps,cosps,k0a02,sk0r;
                         register __m512d cos2ps,sin2ps,cosp,mul1r,mul1i,resr,resi;
                         register __m512d mul2r,mul2i,mul3r,mul3i,t0,t1r,t1i;
                         register __m512d murm1,muim1,epsrm1,epsim1,numr,numi;
                         register __m512d murp1,muip1,epsrp1,epsip1,mucsr,mucsi;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(hlf,_mm512_mul_pd(k0a,k0a));
                         sinps = xsinf(psi);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         t0    = _mm512_fmadd_pd(k0z,sinps,_mm512_fmadd_pd(k0r,cosps,pi4));
                         sk0r  = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                         ear   = Ir;
                         eai   = t0;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         cosp  = xcosf(phi);
                         murm1 = _mm512_sub_pd(mur,_1);
                         fracr = H0r;
                         sin2ps= _mm512_mul_pd(sinps,sinps);
                         muim1 = _mm512_sub_pd(mui,_1);
                         fraci = H0i;
                         epsrm1= _mm512_sub_pd(epsr,_1);
                         epsim1= _mm512_sub_pd(epsi,_1);
                         murp1 = _mm512_add_pd(mur,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(fracr,fraci,cer,cei,&t0r,&t0i);
                         mucsr = _mm512_mul_pd(murm1,cos2ps);
                         t0r = _mm512_mul_pd(e0u0,_mm512_div_pd(t0r,sk0r));
                         mucsi = _mm512_mul_pd(muim1,cos2ps);
                         t0i = _mm512_mul_pd(e0u0,_mm512_div_pd(t0i,sk0r));
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul1r,&mul1i);
                         t0r = _mm512_mul_pd(t0r,k0a02);
                         t0i = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul2r,mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&t1r,&t1i);
                         t1r = _mm512_mul_pd(_2,_mm512_mul_pd(t1r,cosp));
                         t1i = _mm512_mul_pd(_2,_mm512_mul_pd(t1i,cosp));
                         mucsr = _mm512_sub_pd(mucsr,t1r);
                         mucsi = _mm512_sub_pd(mucsi,t1i);
                         cmul_zmm8r8(t0r,t0i,mucsr,mucsi,&resr,&resi);
                         _mm512_storeu_pd(&Epr[0], resr);
                         _mm512_storeu_pd(&Epi[0], resi);
                }


                    /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident H-field.
                         Formula 4.2-53

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hph_f4253_zmm8r8(const __m512d H0r,
                                          const __m512d H0i,
                                          const __m512d k0z,
                                          const __m512d k0r,
                                          const __m512d psi,
                                          const __m512d phi,
                                          const __m512d k0a0,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          __m512d * __restrict Hpr,
                                          __m512d * __restrict Hpi) {

                         const __m512d s2pi = _mm512_set1_pd(2.506628274631000502415765284811);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d sinps,cosps,scpk0r,sinph;
                         register __m512d k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512d spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosps,pi4));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                         cmul_zmm8r8(H0r,H0i,cer,cei,&fracr,&fraci);
                         spsph = _mm512_mul_pd(sinps,sinph);
                         murp1 = _mm512_add_pd(mur,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         t0r   = _mm512_mul_pd(s2pi,_mm512_div_pd(fracr,scpk0r));
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         t0i   = _mm512_mul_pd(s2pi,_mm512_div_pd(fraci,scpk0r));
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,spsph);
                         t0r  = _mm512_mul_pd(t0r,k0a02);
                         divi = _mm512_mul_pd(divi,spsph);
                         t0i  = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(t0r,t0i,divr,divi,*Hpr,*Hpi);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hph_f4253_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          double * __restrict __ATTR_ALIGN__(64) Hpr,
                                          double * __restrict __ATTR_ALIGN__(64) Hpi) {

                         register __m512d H0r   = _mm512_load_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_load_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_load_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d s2pi = _mm512_set1_pd(2.506628274631000502415765284811);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512d k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512d spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosps,pi4));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                         cmul_zmm8r8(H0r,H0i,cer,cei,&fracr,&fraci);
                         spsph = _mm512_mul_pd(sinps,sinph);
                         murp1 = _mm512_add_pd(mur,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         t0r   = _mm512_mul_pd(s2pi,_mm512_div_pd(fracr,scpk0r));
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         t0i   = _mm512_mul_pd(s2pi,_mm512_div_pd(fraci,scpk0r));
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,spsph);
                         t0r  = _mm512_mul_pd(t0r,k0a02);
                         divi = _mm512_mul_pd(divi,spsph);
                         t0i  = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_store_pd(&Hpr[0], resr);
                         _mm512_store_pd(&Hpi[0], resi);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hph_f4253_zmm8r8_u(const double * __restrict  pH0r,
                                          const double * __restrict pH0i,
                                          const double * __restrict  pk0z,
                                          const double * __restrict  pk0r,
                                          const double * __restrict  pk0a0,
                                          const double * __restrict  ppsi,
                                          const double * __restrict  pphi,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          double * __restrict  Hpr,
                                          double * __restrict  Hpi) {

                         register __m512d H0r   = _mm512_loadu_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_loadu_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_loadu_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d s2pi = _mm512_set1_pd(2.506628274631000502415765284811);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512d k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512d spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosps,pi4));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_pd(_mm512_mul_pd(k0r,cosps));
                         cmul_zmm8r8(H0r,H0i,cer,cei,&fracr,&fraci);
                         spsph = _mm512_mul_pd(sinps,sinph);
                         murp1 = _mm512_add_pd(mur,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         t0r   = _mm512_mul_pd(s2pi,_mm512_div_pd(fracr,scpk0r));
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         t0i   = _mm512_mul_pd(s2pi,_mm512_div_pd(fraci,scpk0r));
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,spsph);
                         t0r  = _mm512_mul_pd(t0r,k0a02);
                         divi = _mm512_mul_pd(divi,spsph);
                         t0i  = _mm512_mul_pd(t0i,k0a02);
                         cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_storeu_pd(&Hpr[0], resr);
                         _mm512_storeu_pd(&Hpi[0], resi);
                 }


                   /*
                         Infinitely long cylinder.
                         Scattered fields (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                         TE-incident E-field.
                         Formula 4.2-54

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4254_zmm8r8( const __m512d H0r,
                                          const __m512d H0i,
                                          const __m512d k0z,
                                          const __m512d k0r,
                                          const __m512d psi,
                                          const __m512d phi,
                                          const __m512d k0a0,
                                          const __m512d epsr,
                                          const __m512d epsi,
                                          const __m512d mur,
                                          const __m512d mui,
                                          __m512d * __restrict Ezr,
                                          __m512d * __restrict Ezi) {
                         
                         const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d sinps,cosps,scpk0r,sinph;
                         register __m512d k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512d spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosps,pi4));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_pd(k0r);
                         cmul_zmm8r8(H0r,H0i,cer,cei,&fracr,&fraci);
                         fracr = _mm512_mul_pd(fracr,cosps);
                         spsph = _mm512_mul_pd(sinps,sinph);
                         fraci = _mm512_mul_pd(fraci,cosps);
                         murp1 = _mm512_add_pd(mur,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         t0r   = _mm512_mul_pd(e0u0,_mm512_div_pd(fracr,scpk0r));
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         t0i   = _mm512_mul_pd(e0u0,_mm512_div_pd(fraci,scpk0r));
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,spsph);
                         t0r  = _mm512_mul_pd(Ii,_mm512_mul_pd(t0r,k0a02));
                         divi = _mm512_mul_pd(divi,spsph);
                         t0i  = _mm512_mul_pd(Ii,_mm512_mul_pd(t0i,k0a02));
                         cmul_zmm8r8(t0r,t0i,divr,divi,*Ezr,*Ezi);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4254_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pH0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pH0i,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0z,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0r,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                          const double * __restrict __ATTR_ALIGN__(64) pmur,
                                          const double * __restrict __ATTR_ALIGN__(64) pmui,
                                          double * __restrict __ATTR_ALIGN__(64) Hzr,
                                          double * __restrict __ATTR_ALIGN__(64) Hzi ) {
                         
                         register __m512d H0r   = _mm512_load_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_load_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_load_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512d k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512d spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosps,pi4));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_pd(k0r);
                         cmul_zmm8r8(H0r,H0i,cer,cei,&fracr,&fraci);
                         fracr = _mm512_mul_pd(fracr,cosps);
                         spsph = _mm512_mul_pd(sinps,sinph);
                         fraci = _mm512_mul_pd(fraci,cosps);
                         murp1 = _mm512_add_pd(mur,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         t0r   = _mm512_mul_pd(e0u0,_mm512_div_pd(fracr,scpk0r));
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         t0i   = _mm512_mul_pd(e0u0,_mm512_div_pd(fraci,scpk0r));
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,spsph);
                         t0r  = _mm512_mul_pd(Ii,_mm512_mul_pd(t0r,k0a02));
                         divi = _mm512_mul_pd(divi,spsph);
                         t0i  = _mm512_mul_pd(Ii,_mm512_mul_pd(t0i,k0a02));
                         cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_store_pd(&Hzr[0], resr);
                         _mm512_store_pd(&Hzi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Ez_f4254_zmm8r8_u(const double * __restrict  pH0r,
                                          const double * __restrict pH0i,
                                          const double * __restrict  pk0z,
                                          const double * __restrict  pk0r,
                                          const double * __restrict  pk0a0,
                                          const double * __restrict  ppsi,
                                          const double * __restrict  pphi,
                                          const double * __restrict  pepsr,
                                          const double * __restrict  pepsi,
                                          const double * __restrict  pmur,
                                          const double * __restrict  pmui,
                                          double * __restrict  Hzr,
                                          double * __restrict  Hzi ) {
                         
                         register __m512d H0r   = _mm512_loadu_pd(&pH0r[0]);
                         register __m512d H0i   = _mm512_loadu_pd(&pH0i[0]);
                         register __m512d k0z   = _mm512_loadu_pd(&pk0z[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d e0u0 = _mm512_set1_pd(0.00001763712109284471382861586);
                         const __m512d pi4  = _mm512_set1_pd(0.78539816339744830961566084582);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d sinps,cosps,scpk0r,sinph,resr,resi;
                         register __m512d k0a02,mul1r,mul1i,mul2r,mul2i,divr,divi;
                         register __m512d spsph,epsrp1,epsip1,murp1,muip1;
                         register __m512d fracr,fraci,ear,eai,cer,cei,t0r,t0i;
                         sinps = xsinf(psi);
                         ear   = Ir;
                         cosps = xcosf(psi);
                         k0a02 = _mm512_mul_pd(k0a,k0a);
                         sinph = xsinf(phi);
                         eai   = _mm512_fmadd_pd(k0z,sinps,_mm512_fmsub_pd(k0r,cosps,pi4));
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         scpk0r = _mm512_sqrt_pd(k0r);
                         cmul_zmm8r8(H0r,H0i,cer,cei,&fracr,&fraci);
                         fracr = _mm512_mul_pd(fracr,cosps);
                         spsph = _mm512_mul_pd(sinps,sinph);
                         fraci = _mm512_mul_pd(fraci,cosps);
                         murp1 = _mm512_add_pd(mur,_1);
                         epsrp1= _mm512_add_pd(epsr,_1);
                         muip1 = _mm512_add_pd(mui,_1);
                         epsip1= _mm512_add_pd(epsi,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         t0r   = _mm512_mul_pd(e0u0,_mm512_div_pd(fracr,scpk0r));
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         t0i   = _mm512_mul_pd(e0u0,_mm512_div_pd(fraci,scpk0r));
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,spsph);
                         t0r  = _mm512_mul_pd(Ii,_mm512_mul_pd(t0r,k0a02));
                         divi = _mm512_mul_pd(divi,spsph);
                         t0i  = _mm512_mul_pd(Ii,_mm512_mul_pd(t0i,k0a02));
                         cmul_zmm8r8(t0r,t0i,divr,divi,&resr,&resi);
                         _mm512_storeu_pd(&Hzr[0], resr);
                         _mm512_storeu_pd(&Hzi[0], resi);
                 }


                  /*
                     Bistatic scattering width for (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                     Infinitely long cylinder.
                     TM-incident.
                     Formula 4.2-56
                 */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4256_zmm8r8(const __m512d a0,
                                            const __m512d k0a0,
                                            const __m512d psi,
                                            const __m512d phi,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {

                         const __m512d _4 = _mm512_set1_pd(4.0f);
                         const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512d epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512d epsrcps,epsicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512d t1r,t1i,cabs;
                         k0a03  = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrm1 = _mm512_sub_pd(epsr,_1);
                         spia   = _mm512_mul_pd(spi,a0);
                         epsim1 = _mm512_sub_pd(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_pd(epsr,_1);
                         cos2ps = _mm512_mul_pd(cosps,cosps);
                         epsip1 = _mm512_add_pd(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_pd(mur,_1);
                         sin2ps = _mm512_mul_pd(sinps,sinps);
                         muim1  = _mm512_sub_pd(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_pd(mur,_1);
                         t0     = _mm512_mul_pd(_4,cos2ps);
                         muip1  = _mm512_add_pd(mui,_1);
                         frac   = _mm512_div_pd(spia,t0);
                         epsrcps= _mm512_mul_pd(epsrm1,cos2ps);
                         frac   = _mm512_mul_pd(frac,k0a03);
                         epsicps= _mm512_mul_pd(epsim1,cos2ps);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_pd(_2,_mm512_mul_pd(divr,cosp));
                         t1r  = _mm512_sub_pd(epsrcps,t0r);
                         t0i  = _mm512_mul_pd(_2,_mm512_mul_pd(divi,cosp));
                         t1i  = _mm512_sub_pd(epsicps,t0i);
                         cabs = cabs_zmm8r8(t1r,t1i);
                         rcs  = _mm512_mul_pd(cabs,frac);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4256_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {

                         register __m512d a0   = _mm512_load_pd(&pa0[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d _4 = _mm512_set1_pd(4.0f);
                         const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512d epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512d epsrcps,epsicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512d t1r,t1i,cabs;
                         k0a03  = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrm1 = _mm512_sub_pd(epsr,_1);
                         spia   = _mm512_mul_pd(spi,a0);
                         epsim1 = _mm512_sub_pd(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_pd(epsr,_1);
                         cos2ps = _mm512_mul_pd(cosps,cosps);
                         epsip1 = _mm512_add_pd(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_pd(mur,_1);
                         sin2ps = _mm512_mul_pd(sinps,sinps);
                         muim1  = _mm512_sub_pd(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_pd(mur,_1);
                         t0     = _mm512_mul_pd(_4,cos2ps);
                         muip1  = _mm512_add_pd(mui,_1);
                         frac   = _mm512_div_pd(spia,t0);
                         epsrcps= _mm512_mul_pd(epsrm1,cos2ps);
                         frac   = _mm512_mul_pd(frac,k0a03);
                         epsicps= _mm512_mul_pd(epsim1,cos2ps);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_pd(_2,_mm512_mul_pd(divr,cosp));
                         t1r  = _mm512_sub_pd(epsrcps,t0r);
                         t0i  = _mm512_mul_pd(_2,_mm512_mul_pd(divi,cosp));
                         t1i  = _mm512_sub_pd(epsicps,t0i);
                         cabs = cabs_zmm8r8(t1r,t1i);
                         rcs  = _mm512_mul_pd(cabs,frac);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4256_zmm8r8_u(const double * __restrict  pa0,
                                              const double * __restrict  pk0a0,
                                              const double * __restrict  ppsi,
                                              const double * __restrict  pphi,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {

                         register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d _4 = _mm512_set1_pd(4.0f);
                         const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512d epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512d epsrcps,epsicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512d t1r,t1i,cabs;
                         k0a03  = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrm1 = _mm512_sub_pd(epsr,_1);
                         spia   = _mm512_mul_pd(spi,a0);
                         epsim1 = _mm512_sub_pd(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_pd(epsr,_1);
                         cos2ps = _mm512_mul_pd(cosps,cosps);
                         epsip1 = _mm512_add_pd(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_pd(mur,_1);
                         sin2ps = _mm512_mul_pd(sinps,sinps);
                         muim1  = _mm512_sub_pd(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_pd(mur,_1);
                         t0     = _mm512_mul_pd(_4,cos2ps);
                         muip1  = _mm512_add_pd(mui,_1);
                         frac   = _mm512_div_pd(spia,t0);
                         epsrcps= _mm512_mul_pd(epsrm1,cos2ps);
                         frac   = _mm512_mul_pd(frac,k0a03);
                         epsicps= _mm512_mul_pd(epsim1,cos2ps);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_pd(mul2r,sin2ps,mul1r);
                         numi = _mm512_fmadd_pd(mul2i,sin2ps,mul1i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_pd(_2,_mm512_mul_pd(divr,cosp));
                         t1r  = _mm512_sub_pd(epsrcps,t0r);
                         t0i  = _mm512_mul_pd(_2,_mm512_mul_pd(divi,cosp));
                         t1i  = _mm512_sub_pd(epsicps,t0i);
                         cabs = cabs_zmm8r8(t1r,t1i);
                         rcs  = _mm512_mul_pd(cabs,frac);
                         return (rcs);
                 }


                  /*
                     Bistatic scattering width for (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                     Infinitely long cylinder.
                     TE-incident.
                     Formula 4.2-58
                 */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4258_zmm8r8(const __m512d a0,
                                            const __m512d k0a0,
                                            const __m512d psi,
                                            const __m512d phi,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {

                         const __m512d _4 = _mm512_set1_pd(4.0f);
                         const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512d epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512d murcps,muicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512d t1r,t1i,cabs;
                         k0a03  = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrm1 = _mm512_sub_pd(epsr,_1);
                         spia   = _mm512_mul_pd(spi,a0);
                         epsim1 = _mm512_sub_pd(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_pd(epsr,_1);
                         cos2ps = _mm512_mul_pd(cosps,cosps);
                         epsip1 = _mm512_add_pd(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_pd(mur,_1);
                         sin2ps = _mm512_mul_pd(sinps,sinps);
                         muim1  = _mm512_sub_pd(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_pd(mur,_1);
                         t0     = _mm512_mul_pd(_4,cos2ps);
                         muip1  = _mm512_add_pd(mui,_1);
                         frac   = _mm512_div_pd(spia,t0);
                         murcps= _mm512_mul_pd(murm1,cos2ps);
                         frac   = _mm512_mul_pd(frac,k0a03);
                         muicps= _mm512_mul_pd(muim1,cos2ps);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_pd(mul1r,sin2ps,mul2r);
                         numi = _mm512_fmadd_pd(mul1i,sin2ps,mul2i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_pd(_2,_mm512_mul_pd(divr,cosp));
                         t1r  = _mm512_sub_pd(murcps,t0r);
                         t0i  = _mm512_mul_pd(_2,_mm512_mul_pd(divi,cosp));
                         t1i  = _mm512_sub_pd(muicps,t0i);
                         cabs = cabs_zmm8r8(t1r,t1i);
                         rcs  = _mm512_mul_pd(cabs,frac);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4258_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {

                         register __m512d a0    = _mm512_load_pd(&pa0[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d _4 = _mm512_set1_pd(4.0f);
                         const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512d epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512d murcps,muicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512d t1r,t1i,cabs;
                         k0a03  = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrm1 = _mm512_sub_pd(epsr,_1);
                         spia   = _mm512_mul_pd(spi,a0);
                         epsim1 = _mm512_sub_pd(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_pd(epsr,_1);
                         cos2ps = _mm512_mul_pd(cosps,cosps);
                         epsip1 = _mm512_add_pd(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_pd(mur,_1);
                         sin2ps = _mm512_mul_pd(sinps,sinps);
                         muim1  = _mm512_sub_pd(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_pd(mur,_1);
                         t0     = _mm512_mul_pd(_4,cos2ps);
                         muip1  = _mm512_add_pd(mui,_1);
                         frac   = _mm512_div_pd(spia,t0);
                         murcps= _mm512_mul_pd(murm1,cos2ps);
                         frac   = _mm512_mul_pd(frac,k0a03);
                         muicps= _mm512_mul_pd(muim1,cos2ps);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_pd(mul1r,sin2ps,mul2r);
                         numi = _mm512_fmadd_pd(mul1i,sin2ps,mul2i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_pd(_2,_mm512_mul_pd(divr,cosp));
                         t1r  = _mm512_sub_pd(murcps,t0r);
                         t0i  = _mm512_mul_pd(_2,_mm512_mul_pd(divi,cosp));
                         t1i  = _mm512_sub_pd(muicps,t0i);
                         cabs = cabs_zmm8r8(t1r,t1i);
                         rcs  = _mm512_mul_pd(cabs,frac);
                         return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4258_zmm8r8_u(const double * __restrict  pa0,
                                              const double * __restrict  pk0a0,
                                              const double * __restrict  ppsi,
                                              const double * __restrict  pphi,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {

                         register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d _4 = _mm512_set1_pd(4.0f);
                         const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2 = _mm512_set1_pd(2.0f);
                         register __m512d rcs,k0a03,frac,cosp,cos2ps,cosps,sinps,sin2ps,spia,t0,t0r,t0i;
                         register __m512d epsrm1,epsim1,epsrp1,epsip1,murm1,muim1,murp1,muip1,numr,numi;
                         register __m512d murcps,muicps,divr,divi,mul1r,mul1i,mul2r,mul2i,mul3r,mul3i;
                         register __m512d t1r,t1i,cabs;
                         k0a03  = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrm1 = _mm512_sub_pd(epsr,_1);
                         spia   = _mm512_mul_pd(spi,a0);
                         epsim1 = _mm512_sub_pd(epsi,_1);
                         cosps  = xcosf(psi);
                         epsrp1 = _mm512_add_pd(epsr,_1);
                         cos2ps = _mm512_mul_pd(cosps,cosps);
                         epsip1 = _mm512_add_pd(epsi,_1);
                         sinps  = xsinf(ps);
                         murm1  = _mm512_sub_pd(mur,_1);
                         sin2ps = _mm512_mul_pd(sinps,sinps);
                         muim1  = _mm512_sub_pd(mui,_1);
                         cosp   = xcosf(phi);
                         murp1  = _mm512_add_pd(mur,_1);
                         t0     = _mm512_mul_pd(_4,cos2ps);
                         muip1  = _mm512_add_pd(mui,_1);
                         frac   = _mm512_div_pd(spia,t0);
                         murcps= _mm512_mul_pd(murm1,cos2ps);
                         frac   = _mm512_mul_pd(frac,k0a03);
                         muicps= _mm512_mul_pd(muim1,cos2ps);
                         cmul_zmm8r8(epsrp1,epsip1,murm1,muim1,&mul1r,&mul1i);
                         cmul_zmm8r8(epsrm1,epsim1,murp1,muip1,&mul2r,&mul2i);
                         numr = _mm512_fmadd_pd(mul1r,sin2ps,mul2r);
                         numi = _mm512_fmadd_pd(mul1i,sin2ps,mul2i);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul3r,&mul3i);
                         cdiv_zmm8r8(numr,numi,mul3r,mul3i,&divr,&divi);
                         t0r  = _mm512_mul_pd(_2,_mm512_mul_pd(divr,cosp));
                         t1r  = _mm512_sub_pd(murcps,t0r);
                         t0i  = _mm512_mul_pd(_2,_mm512_mul_pd(divi,cosp));
                         t1i  = _mm512_sub_pd(muicps,t0i);
                         cabs = cabs_zmm8r8(t1r,t1i);
                         rcs  = _mm512_mul_pd(cabs,frac);
                         return (rcs);
                 }


                   /*

                           Bistatic scattering width for (k0a0 sqrt(epsr*mur-sin^2(Psi) < 0.5)
                           Infinitely long cylinder.
                           TM-incident.
                           Formula 4.2-57
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4257_zmm8r8(const __m512d a0,
                                            const __m512d k0a0,
                                            const __m512d psi,
                                            const __m512d phi,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {

                         const __m512d _4  = _mm512_set1_pd(4.0f);
                         const __m512d spi = _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2  = _mm512_set1_pd(2.0f);
                         register __m512d spi4,rcs,cos2ps,sinps,sinp,k0a03;
                         register __m512d frac,divr,divi,mul1r,mul1i,mul2r,mul2i;
                         register __m512d epsrp1,epsip1,murp1,muip1,t0,cabs;
                         k0a03 = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrp1= _mm512_add_pd(epsr,_1);
                         cosps = xcosf(psi);
                         epsip1= _mm512_add_pd(epsi,_1);
                         spi4  = _mm512_mul_pd(spi,_mm512_mul_pd(a0,a0));
                         spi4  = _mm512_mul_pd(_4,spi4);
                         murp1 = _mm512_add_pd(mur,_1);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         muip1 = _mm512_add_pd(mui,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         sinps = xsinf(psi);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         sinp  = xsinf(phi);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         frac  = _mm512_mul_pd(k0a03,_mm512_div_pd(spi4,cos2ps));
                         t0    = _mm512_mul_pd(sinps,sinp);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,t0);
                         divi = _mm512_mul_pd(divi,t0);
                         cabs = cabs_zmm8r8(divr,divi);
                         rcs  = _mm512_mul_pd(frac,cabs);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4257_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa0,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a0,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {

                         register __m512d a0    = _mm512_load_pd(&pa0[0]);
                         register __m512d k0a0  = _mm512_load_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_load_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_load_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_load_pd(&pmur[0]);
                         register __m512d mui   = _mm512_load_pd(&pmui[0]); 
                         const __m512d _4  = _mm512_set1_pd(4.0f);
                         const __m512d spi = _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2  = _mm512_set1_pd(2.0f);
                         register __m512d spi4,rcs,cos2ps,sinps,sinp,k0a03;
                         register __m512d frac,divr,divi,mul1r,mul1i,mul2r,mul2i;
                         register __m512d epsrp1,epsip1,murp1,muip1,t0,cabs;
                         k0a03 = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrp1= _mm512_add_pd(epsr,_1);
                         cosps = xcosf(psi);
                         epsip1= _mm512_add_pd(epsi,_1);
                         spi4  = _mm512_mul_pd(spi,_mm512_mul_pd(a0,a0));
                         spi4  = _mm512_mul_pd(_4,spi4);
                         murp1 = _mm512_add_pd(mur,_1);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         muip1 = _mm512_add_pd(mui,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         sinps = xsinf(psi);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         sinp  = xsinf(phi);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         frac  = _mm512_mul_pd(k0a03,_mm512_div_pd(spi4,cos2ps));
                         t0    = _mm512_mul_pd(sinps,sinp);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,t0);
                         divi = _mm512_mul_pd(divi,t0);
                         cabs = cabs_zmm8r8(divr,divi);
                         rcs  = _mm512_mul_pd(frac,cabs);
                         return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4257_zmm8r8_u(const double * __restrict  pa0,
                                              const double * __restrict  pk0a0,
                                              const double * __restrict  ppsi,
                                              const double * __restrict  pphi,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {

                         register __m512d a0    = _mm512_loadu_pd(&pa0[0]);
                         register __m512d k0a0  = _mm512_loadu_pd(&pk0a0[0]);
                         register __m512d psi   = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d pphi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                         register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                         register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                         register __m512d mui   = _mm512_loadu_pd(&pmui[0]); 
                         const __m512d _4  = _mm512_set1_pd(4.0f);
                         const __m512d spi = _mm512_set1_pd(9.869604401089358618834490999876);
                         const __m512d _2  = _mm512_set1_pd(2.0f);
                         register __m512d spi4,rcs,cos2ps,sinps,sinp,k0a03;
                         register __m512d frac,divr,divi,mul1r,mul1i,mul2r,mul2i;
                         register __m512d epsrp1,epsip1,murp1,muip1,t0,cabs;
                         k0a03 = _mm512_mul_pd(k0a0,_mm512_mul_pd(k0a0,k0a0));
                         epsrp1= _mm512_add_pd(epsr,_1);
                         cosps = xcosf(psi);
                         epsip1= _mm512_add_pd(epsi,_1);
                         spi4  = _mm512_mul_pd(spi,_mm512_mul_pd(a0,a0));
                         spi4  = _mm512_mul_pd(_4,spi4);
                         murp1 = _mm512_add_pd(mur,_1);
                         cos2ps= _mm512_mul_pd(cosps,cosps);
                         muip1 = _mm512_add_pd(mui,_1);
                         cmul_zmm8r8(epsr,epsi,mur,mui,&mul1r,&mul1i);
                         sinps = xsinf(psi);
                         mul1r = _mm512_sub_pd(mul1r,_1);
                         mul1i = _mm512_sub_pd(mul1i,_1);
                         sinp  = xsinf(phi);
                         cmul_zmm8r8(epsrp1,epsip1,murp1,muip1,&mul2r,&mul2i);
                         frac  = _mm512_mul_pd(k0a03,_mm512_div_pd(spi4,cos2ps));
                         t0    = _mm512_mul_pd(sinps,sinp);
                         cdiv_zmm8r8(mul1r,mul1i,mul2r,mul2i,&divr,&divi);
                         divr = _mm512_mul_pd(divr,t0);
                         divi = _mm512_mul_pd(divi,t0);
                         cabs = cabs_zmm8r8(divr,divi);
                         rcs  = _mm512_mul_pd(frac,cabs);
                         return (rcs);
                 }


                   /*
                       Circular cylinders of finite length.
                       Cylinder radius small (k0a<1.0)
                       Wire limit of cylinder (h>>a).
                       E-field
                       Formula 4.3-9
                   */

                    
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f439_zmm8r8(const __m512d EIr,
                                        const __m512d EIi,
                                        const __m512d r,
                                        const __m512d k0,
                                        const __m512d psii,
                                        const __m512d psis,
                                        const __m512d h,
                                        const __m512d ln4ha,
                                        __m512d * __restrict ESr,
                                        __m512d * __restrict ESi) {

                       const __m512d thrd = _mm512_set1_pd(0.333333333333333333333333333333f);
                       const __m512d _1   = _mm512_set1_pd(1.0f);
                       register __m512d ir,k02,ear,eai,cer,cei,h3,cpsii,cpsis;
                       register __m512d num,rat,den,t0r,t0i,mulr,muli;
                       cpsii = xcosf(psii);
                       k02   = _mm512_mul_pd(thrd,_mm512_mul_pd(k0,k0));
                       ir    = _mm512_rcp14_pd(r);
                       ear   = Ir;
                       cpsis = xcosf(psis);
                       eai   = _mm512_mul_pd(k0,r);
                       den   = _mm512_sub_pd(ln4ha,_1);
                       h3    = _mm512_mul_pd(h,_mm512_mul_pd(h,h));
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_pd(cer,ir);
                       num   = _mm512_mul_pd(h3,_mm512_mul_pd(cpsis,cpsii));
                       cei   = _mm512_mul_pd(cei,ir);
                       rat   = _mm512_div_pd(num,den);
                       t0r   = _mm512_mul_pd(EIr,rat);
                       t0i   = _mm512_mul_pd(EIi,rat);
                       cmul_zmm8r8(cer,cei,t0r,t0i,&mulr,&muli);
                       *ESr = _mm512_mul_pd(mulr,k02);
                       *ESi = _mm512_mul_pd(muli,k02);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f439_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pEIr,
                                          const double * __restrict __ATTR_ALIGN__(64) pEIi,
                                          const double * __restrict __ATTR_ALIGN__(64) pr,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                          const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                          const double * __restrict __ATTR_ALIGN__(64) ph,
                                          const double * __restrict __ATTR_ALIGN__(64) pln4ha,
                                          double * __restrict __ATTR_ALIGN__(64) ESr,
                                          double * __restrict __ATTR_ALIGN__(64) ESi) {

                       register __m512d EIr  = _mm512_load_pd(&pEIr[0]);
                       register __m512d EIi  = _mm512_load_pd(&pEIi[0]);
                       register __m512d r    = _mm512_load_pd(&pr[0]);
                       register __m512d k0   = _mm512_load_pd(&pk0[0]);
                       register __m512d psii = _mm512_load_pd(&ppsii[0]);
                       register __m512d psis = _mm512_load_pd(&ppsis[0]);
                       register __m512d h    = _mm512_load_pd(&ph[0]);
                       register __m512d ln4ha= _mm512_load_pd(&pln4ha[0]);
                       const __m512d thrd = _mm512_set1_pd(0.333333333333333333333333333333f);
                       const __m512d _1   = _mm512_set1_pd(1.0f);
                       register __m512d ir,k02,ear,eai,cer,cei,h3,cpsii,cpsis;
                       register __m512d num,rat,den,t0r,t0i,mulr,muli;
                       cpsii = xcosf(psii);
                       k02   = _mm512_mul_pd(thrd,_mm512_mul_pd(k0,k0));
                       ir    = _mm512_rcp14_pd(r);
                       ear   = Ir;
                       cpsis = xcosf(psis);
                       eai   = _mm512_mul_pd(k0,r);
                       den   = _mm512_sub_pd(ln4ha,_1);
                       h3    = _mm512_mul_pd(h,_mm512_mul_pd(h,h));
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_pd(cer,ir);
                       num   = _mm512_mul_pd(h3,_mm512_mul_pd(cpsis,cpsii));
                       cei   = _mm512_mul_pd(cei,ir);
                       rat   = _mm512_div_pd(num,den);
                       t0r   = _mm512_mul_pd(EIr,rat);
                       t0i   = _mm512_mul_pd(EIi,rat);
                       cmul_zmm8r8(cer,cei,t0r,t0i,&mulr,&muli);
                       _mm512_store_pd(&ESr[0] ,_mm512_mul_pd(mulr,k02));
                       _mm512_store_pd(&ESi[0] ,_mm512_mul_pd(muli,k02));
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f439_zmm8r8_u(const double * __restrict  pEIr,
                                          const double * __restrict  pEIi,
                                          const double * __restrict  pr,
                                          const double * __restrict  pk0,
                                          const double * __restrict  ppsii,
                                          const double * __restrict  ppsis,
                                          const double * __restrict  ph,
                                          const double * __restrict  pln4ha,
                                          double * __restrict  ESr,
                                          double * __restrict  ESi) {

                       register __m512d EIr  = _mm512_loadu_pd(&pEIr[0]);
                       register __m512d EIi  = _mm512_loadu_pd(&pEIi[0]);
                       register __m512d r    = _mm512_loadu_pd(&pr[0]);
                       register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                       register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                       register __m512d psis = _mm512_loadu_pd(&ppsis[0]);
                       register __m512d h    = _mm512_loadu_pd(&ph[0]);
                       register __m512d ln4ha= _mm512_loadu_pd(&pln4ha[0]);
                       const __m512d thrd = _mm512_set1_pd(0.333333333333333333333333333333f);
                       const __m512d _1   = _mm512_set1_pd(1.0f);
                       register __m512d ir,k02,ear,eai,cer,cei,h3,cpsii,cpsis;
                       register __m512d num,rat,den,t0r,t0i,mulr,muli;
                       cpsii = xcosf(psii);
                       k02   = _mm512_mul_pd(thrd,_mm512_mul_pd(k0,k0));
                       ir    = _mm512_rcp14_pd(r);
                       ear   = Ir;
                       cpsis = xcosf(psis);
                       eai   = _mm512_mul_pd(k0,r);
                       den   = _mm512_sub_pd(ln4ha,_1);
                       h3    = _mm512_mul_pd(h,_mm512_mul_pd(h,h));
                       cexp_zmm8r8(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_pd(cer,ir);
                       num   = _mm512_mul_pd(h3,_mm512_mul_pd(cpsis,cpsii));
                       cei   = _mm512_mul_pd(cei,ir);
                       rat   = _mm512_div_pd(num,den);
                       t0r   = _mm512_mul_pd(EIr,rat);
                       t0i   = _mm512_mul_pd(EIi,rat);
                       cmul_zmm8r8(cer,cei,t0r,t0i,&mulr,&muli);
                       _mm512_storeu_pd(&ESr[0] ,_mm512_mul_pd(mulr,k02));
                       _mm512_storeu_pd(&ESi[0] ,_mm512_mul_pd(muli,k02));
                 }


                  /*
                       Circular cylinders of finite length.
                       Cylinder radius small (k0a<1.0)
                       Wire limit of cylinder (h>>a).
                       RCS.
                       Formula 4.3-10

                    */
                  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4310_zmm8r8(const __m512d k0,
                                            const __m512d h,
                                            const __m512d psii,
                                            const __m512d psis,
                                            const __m512d ln4h) {

                          const __m512d _4pi9 = _mm512_set1_pd(1.396263401595463661538952614791);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d cpsii,cpsis,c2psii,c2psis,den,num,t0;
                          register __m512d k04,h6,rcs,t1,h2,rat,frac;
                          h2     = _mm512_mul_pd(h,h);
                          k04    = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
                                              _mm512_mul_pd(k0,k0));
                          t0     = _mm512_sub_pd(ln4h,_1);
                          cpsii  = xcosf(psii);
                          c2psii = _mm512_mul_pd(cpsii,cpsii);
                          den    = _mm512_mul_pd(t0,t0); 
                          t1     = _mm512_mul_pd(h,h2);
                          h6     = _mm512_mul_pd(t1,h2);
                          cpsis  = xcosf(psis);
                          c2psis = _mm512_mul_pd(cpsis,cpsis);
                          num    = _mm512_mul_pd(c2psis,c2psii);
                          frac   = _mm512_mul_pd(_4pi9,_mm512_mul_pd(k04,h6));
                          rat    = _mm512_div_pd(num,den);
                          rcs    = _mm512_mul_pd(frac,rat);
                          return (rcs);
                 } 


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4310_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const double * __restrict __ATTR_ALIGN__(64)  ph,
                                              const double * __restrict __ATTR_ALIGN__(64)  ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64)  ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64)  pln4h) {

                          register __m512d k0   = _mm512_load_pd(&pk0[0]);
                          register __m512d h    = _mm512_load_pd(&ph[0]);
                          register __m512d psii = _mm512_load_pd(&ppsii[0]);
                          register __m512d psis = _mm512_load_pd(&ppsis[0]);
                          register __m512d pln4h= _mm512_load_pd(&pln4h[0]);
                          const __m512d _4pi9 = _mm512_set1_pd(1.396263401595463661538952614791);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d cpsii,cpsis,c2psii,c2psis,den,num,t0;
                          register __m512d k04,h6,rcs,t1,h2,rat,frac;
                          h2     = _mm512_mul_pd(h,h);
                          k04    = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
                                              _mm512_mul_pd(k0,k0));
                          t0     = _mm512_sub_pd(ln4h,_1);
                          cpsii  = xcosf(psii);
                          c2psii = _mm512_mul_pd(cpsii,cpsii);
                          den    = _mm512_mul_pd(t0,t0); 
                          t1     = _mm512_mul_pd(h,h2);
                          h6     = _mm512_mul_pd(t1,h2);
                          cpsis  = xcosf(psis);
                          c2psis = _mm512_mul_pd(cpsis,cpsis);
                          num    = _mm512_mul_pd(c2psis,c2psii);
                          frac   = _mm512_mul_pd(_4pi9,_mm512_mul_pd(k04,h6));
                          rat    = _mm512_div_pd(num,den);
                          rcs    = _mm512_mul_pd(frac,rat);
                          return (rcs);
                 } 


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4310_zmm8r8_u(const double * __restrict   pk0,
                                              const double * __restrict   ph,
                                              const double * __restrict   ppsii,
                                              const double * __restrict   ppsis,
                                              const double * __restrict   pln4h) {

                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                          register __m512d h    = _mm512_loadu_pd(&ph[0]);
                          register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d psis = _mm512_loadu_pd(&ppsis[0]);
                          register __m512d pln4h= _mm512_loadu_pd(&pln4h[0]);
                          const __m512d _4pi9 = _mm512_set1_pd(1.396263401595463661538952614791);
                          const __m512d _1    = _mm512_set1_pd(1.0f);
                          register __m512d cpsii,cpsis,c2psii,c2psis,den,num,t0;
                          register __m512d k04,h6,rcs,t1,h2,rat,frac;
                          h2     = _mm512_mul_pd(h,h);
                          k04    = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
                                              _mm512_mul_pd(k0,k0));
                          t0     = _mm512_sub_pd(ln4h,_1);
                          cpsii  = xcosf(psii);
                          c2psii = _mm512_mul_pd(cpsii,cpsii);
                          den    = _mm512_mul_pd(t0,t0); 
                          t1     = _mm512_mul_pd(h,h2);
                          h6     = _mm512_mul_pd(t1,h2);
                          cpsis  = xcosf(psis);
                          c2psis = _mm512_mul_pd(cpsis,cpsis);
                          num    = _mm512_mul_pd(c2psis,c2psii);
                          frac   = _mm512_mul_pd(_4pi9,_mm512_mul_pd(k04,h6));
                          rat    = _mm512_div_pd(num,den);
                          rcs    = _mm512_mul_pd(frac,rat);
                          return (rcs);
                 } 


                  /*
                         The average dipole scattering RCS when the incidence
                         and scattered polarization direction coincide.
                         Formula 4.3-11
                    */

                   
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4311_zmm8r8(const __m512d k0,
                                            const __m512d h,
                                            const __m512d ln4h) {

                          const __m512d _4pi45 = _mm512_set1_pd(0.279252680319092732307790522958);
                          const __m512d _1     = _mm512_set1_pd(1.0f);
                          register __m512d rcs,den,inv,k04,h6,h2,t0,t1;
                          h2     = _mm512_mul_pd(h,h);
                          k04    = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
                                              _mm512_mul_pd(k0,k0));
                          t0     = _mm512_sub_pd(ln4h,_1);
                          t1     = _mm512_mul_pd(h,h2);
                          h6     = _mm512_mul_pd(t1,h2);
                          den    = _mm512_mul_pd(t0,t0);
                          inv    = _mm512_div_pd(_1,den);
                          t0     = _mm512_mul_pd(_4pi45,_mm512_mul_pd(k04,h6));
                          rcs    = _mm512_mul_pd(t0,inv);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4311_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) ph,
                                              const double * __restrict __ATTR_ALIGN__(64) pln4h) {

                          register __m512d k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d h   = _mm512_load_pd(&ph[0]);
                          register __m512d ln4h= _mm512_load_pd(&pln4h[0]);
                          const __m512d _4pi45 = _mm512_set1_pd(0.279252680319092732307790522958);
                          const __m512d _1     = _mm512_set1_pd(1.0f);
                          register __m512d rcs,den,inv,k04,h6,h2,t0,t1;
                          h2     = _mm512_mul_pd(h,h);
                          k04    = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
                                              _mm512_mul_pd(k0,k0));
                          t0     = _mm512_sub_pd(ln4h,_1);
                          t1     = _mm512_mul_pd(h,h2);
                          h6     = _mm512_mul_pd(t1,h2);
                          den    = _mm512_mul_pd(t0,t0);
                          inv    = _mm512_div_pd(_1,den);
                          t0     = _mm512_mul_pd(_4pi45,_mm512_mul_pd(k04,h6));
                          rcs    = _mm512_mul_pd(t0,inv);
                          return (rcs);
               }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4311_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  ph,
                                              const double * __restrict  pln4h) {

                          register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d h   = _mm512_loadu_pd(&ph[0]);
                          register __m512d ln4h= _mm512_loadu_pd(&pln4h[0]);
                          const __m512d _4pi45 = _mm512_set1_pd(0.279252680319092732307790522958);
                          const __m512d _1     = _mm512_set1_pd(1.0f);
                          register __m512d rcs,den,inv,k04,h6,h2,t0,t1;
                          h2     = _mm512_mul_pd(h,h);
                          k04    = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
                                              _mm512_mul_pd(k0,k0));
                          t0     = _mm512_sub_pd(ln4h,_1);
                          t1     = _mm512_mul_pd(h,h2);
                          h6     = _mm512_mul_pd(t1,h2);
                          den    = _mm512_mul_pd(t0,t0);
                          inv    = _mm512_div_pd(_1,den);
                          t0     = _mm512_mul_pd(_4pi45,_mm512_mul_pd(k04,h6));
                          rcs    = _mm512_mul_pd(t0,inv);
                          return (rcs);
               }


                  /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-18
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4318_zmm8r8(const __m512d EIr,
                                         const __m512d EIi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d psii,
                                         const __m512d psis,
                                         const __m512d phi,
                                         const __m512d a,
                                         __m512d * __restrict ESr,
                                         __m512d * __restrict ESi) {

                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         register __m512d ir,a3,k02,cosp,spsii,spsis,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         spsis= xsinf(psis);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_mul_pd(spsii,_mm512_mul_pd(spsis,cosp));
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4318_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const double * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const double * __restrict __ATTR_ALIGN__(64)  ppsii,
                                           const double * __restrict __ATTR_ALIGN__(64)  ppsis,
                                           const double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pa,
                                           double * __restrict __ATTR_ALIGN__(64) ESr,
                                           double * __restrict __ATTR_ALIGN__(64)  ESi) {

                         register __m512d EIr  = _mm512_load_pd(&pEIr[0]);
                         register __m512d EIi  = _mm512_load_pd(&pEIi[0]);
                         register __m512d k0   = _mm512_load_pd(&pk0[0]);
                         register __m512d r    = _mm512_load_pd(&pr[0]);
                         register __m512d psii = _mm512_load_pd(&ppsii[0]);
                         register __m512d psis = _mm512_load_pd(&ppsis[0]);
                         register __m512d phi  = _mm512_load_pd(&pphi[0]);
                         register __m512d a    = _mm512_load_pd(&a[0]);
                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         register __m512d ir,a3,k02,cosp,spsii,spsis,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         spsis= xsinf(psis);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_mul_pd(spsii,_mm512_mul_pd(spsis,cosp));
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_store_pd(&ESr[0], resr);
                         _mm512_store_pd(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4318_zmm8r8_u(const double * __restrict   pEIr,
                                           const double * __restrict   pEIi,
                                           const double * __restrict   pk0,
                                           const double * __restrict   pr,
                                           const double * __restrict   ppsii,
                                           const double * __restrict   ppsis,
                                           const double * __restrict   pphi,
                                           const double * __restrict   pa,
                                           double * __restrict  ESr,
                                           double * __restrict   ESi) {

                         register __m512d EIr  = _mm512_loadu_pd(&pEIr[0]);
                         register __m512d EIi  = _mm512_loadu_pd(&pEIi[0]);
                         register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                         register __m512d r    = _mm512_loadu_pd(&pr[0]);
                         register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                         register __m512d psis = _mm512_loadu_pd(&ppsis[0]);
                         register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d a    = _mm512_loadu_pd(&a[0]);
                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         register __m512d ir,a3,k02,cosp,spsii,spsis,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         spsis= xsinf(psis);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_mul_pd(spsii,_mm512_mul_pd(spsis,cosp));
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_storeu_pd(&ESr[0], resr);
                         _mm512_storeu_pd(&ESi[0], resi);
                }


                   /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-19
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4319_zmm8r8(const __m512d EIr,
                                         const __m512d EIi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d psii,
                                         const __m512d phi,
                                         const __m512d a,
                                         __m512d * __restrict ESr,
                                         __m512d * __restrict ESi) {

                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         register __m512d ir,a3,k02,sinp,spsii,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         sinp = xsinf(phi);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_mul_pd(spsii,sinp);
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,*ESr,*ESi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4319_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const double * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const double * __restrict __ATTR_ALIGN__(64)  ppsii,
                                           const double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pa,
                                           double * __restrict __ATTR_ALIGN__(64) ESr,
                                           double * __restrict __ATTR_ALIGN__(64)  ESi) {

                         register __m512d EIr  = _mm512_load_pd(&pEIr[0]);
                         register __m512d EIi  = _mm512_load_pd(&pEIi[0]);
                         register __m512d k0   = _mm512_load_pd(&pk0[0]);
                         register __m512d r    = _mm512_load_pd(&pr[0]);
                         register __m512d psii = _mm512_load_pd(&ppsii[0]);
                         register __m512d phi  = _mm512_load_pd(&pphi[0]);
                         register __m512d a    = _mm512_load_pd(&a[0]);
                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         register __m512d ir,a3,k02,sinp,spsii,t0,t1,resr,resi;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         sinp = xsinf(phi);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_mul_pd(spsii,sinp);
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_store_pd(&ESr[0], resr);
                         _mm512_store_pd(&ESi[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4319_zmm8r8_u(const double * __restrict   pEIr,
                                           const double * __restrict   pEIi,
                                           const double * __restrict   pk0,
                                           const double * __restrict   pr,
                                           const double * __restrict   ppsii,
                                           const double * __restrict   pphi,
                                           const double * __restrict   pa,
                                           double * __restrict  ESr,
                                           double * __restrict  ESi) {

                         register __m512d EIr  = _mm512_loadu_pd(&pEIr[0]);
                         register __m512d EIi  = _mm512_loadu_pd(&pEIi[0]);
                         register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                         register __m512d r    = _mm512_loadu_pd(&pr[0]);
                         register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                         register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d a    = _mm512_loadu_pd(&a[0]);
                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         register __m512d ir,a3,k02,sinp,spsii,t0,t1,resr,resi;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         sinp = xsinf(phi);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         spsii= xsinf(psii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_mul_pd(spsii,sinp);
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_storeu_pd(&ESr[0], resr);
                         _mm512_storeu_pd(&ESi[0], resi);
                }


                 /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-20
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4320_zmm8r8(const __m512d EIr,
                                         const __m512d EIi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d psis,
                                         const __m512d phi,
                                         const __m512d a,
                                         __m512d * __restrict ESr,
                                         __m512d * __restrict ESi) {
                      
                        ES_f4319_zmm8r8(EIr,EIi,k0,r,psis,phi,a,&ESr,&ESi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4320_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const double * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const double * __restrict __ATTR_ALIGN__(64)  ppsis,
                                           const double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pa,
                                           double * __restrict __ATTR_ALIGN__(64) ESr,
                                           double * __restrict __ATTR_ALIGN__(64)  ESi) {

                      ES_f4319_zmm8r8_a(pEIr,pEIi,pk0,pr,ppsis,pphi,ESr,ESi);
              }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4320_zmm8r8_u(const double * __restrict   pEIr,
                                           const double * __restrict   pEIi,
                                           const double * __restrict   pk0,
                                           const double * __restrict   pr,
                                           const double * __restrict   ppsis,
                                           const double * __restrict   pphi,
                                           const double * __restrict   pa,
                                           double * __restrict  ESr,
                                           double * __restrict  ESi) {

                      ES_f4319_zmm8r8_u(pEIr,pEIi,pk0,pr,ppsis,pphi,ESr,ESi);
              }


                   /*
                           Disc limit of cylinder (h<<a).
                           Scattered fields from the cylinder in the disc limit
                           Formula 4.3-21
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4321_zmm8r8(const __m512d EIr,
                                         const __m512d EIi,
                                         const __m512d k0,
                                         const __m512d r,
                                         const __m512d psii,
                                         const __m512d psis,
                                         const __m512d phi,
                                         const __m512d a,
                                         __m512d * __restrict ESr,
                                         __m512d * __restrict ESi) {

                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         const __m512d hlf   = _mm512_set1_pd(0.5f);
                         register __m512d ir,a3,k02,cosp,cpsii,cpsis,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         cpsis= xcosf(psis);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         cpsii= xcosf(psii);
                         cpsii= _mm512_mul_pd(hlf,cpsii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_fmadd_pd(cpsii,cpsis,cosp);
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,*ESr,*ESi);
              }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4321_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pEIr,
                                           const double * __restrict __ATTR_ALIGN__(64)  pEIi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const double * __restrict __ATTR_ALIGN__(64)  pr,
                                           const double * __restrict __ATTR_ALIGN__(64)  ppsii,
                                           const double * __restrict __ATTR_ALIGN__(64)  ppsis,
                                           const double * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const double * __restrict __ATTR_ALIGN__(64)  pa,
                                           double * __restrict __ATTR_ALIGN__(64) ESr,
                                           double * __restrict __ATTR_ALIGN__(64)  ESi) {

                         register __m512d EIr  = _mm512_load_pd(&pEIr[0]);
                         register __m512d EIi  = _mm512_load_pd(&pEIi[0]);
                         register __m512d k0   = _mm512_load_pd(&pk0[0]);
                         register __m512d r    = _mm512_load_pd(&pr[0]);
                         register __m512d psii = _mm512_load_pd(&ppsii[0]);
                         register __m512d psis = _mm512_load_pd(&ppsis[0]);
                         register __m512d phi  = _mm512_load_pd(&pphi[0]);
                         register __m512d a    = _mm512_load_pd(&a[0]);
                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         const __m512d hlf   = _mm512_set1_pd(0.5f);
                         register __m512d ir,a3,k02,cosp,cpsii,cpsis,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         cpsis= xcosf(psis);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         cpsii= xcosf(psii);
                         cpsii= _mm512_mul_pd(hlf,cpsii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_fmadd_pd(cpsii,cpsis,cosp);
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_store_pd(&ESr[0], resr);
                         _mm512_store_pd(&ESi[0], resi);
              }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void ES_f4321_zmm8r8_u(const double * __restrict   pEIr,
                                           const double * __restrict   pEIi,
                                           const double * __restrict   pk0,
                                           const double * __restrict   pr,
                                           const double * __restrict   ppsii,
                                           const double * __restrict   ppsis,
                                           const double * __restrict   pphi,
                                           const double * __restrict   pa,
                                           double * __restrict  ESr,
                                           double * __restrict  ESi) {

                         register __m512d EIr  = _mm512_loadu_pd(&pEIr[0]);
                         register __m512d EIi  = _mm512_loadu_pd(&pEIi[0]);
                         register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                         register __m512d r    = _mm512_loadu_pd(&pr[0]);
                         register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                         register __m512d psis = _mm512_loadu_pd(&ppsis[0]);
                         register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                         register __m512d a    = _mm512_loadu_pd(&a[0]);
                         const __m512d _43pi = _mm512_set1_pd(0.424413181578387562050356702327);
                         const __m512d hlf   = _mm512_set1_pd(0.5f);
                         register __m512d ir,a3,k02,cosp,cpsii,cpsis,t0,t1;
                         register __m512d ear,eai,cer,cei,t0r,t0i,t1r,t1i,resr,resi;
                         a3   = _mm512_mul_pd(a,_mm512_mul_pd(a,a));
                         ir   = _mm512_rcp14_pd(r);
                         cpsis= xcosf(psis);
                         k02  = _mm512_mul_pd(k0,k0);
                         ear  = _mm512_mul_pd(k0,r);
                         cosp = xcosf(phi); 
                         eai  = Ir;
                         cexp_zmm8r8(ear,eai,&cer,&cei);
                         t0   = _mm512_mul_pd(_43pi,k02);
                         cpsii= xcosf(psii);
                         cpsii= _mm512_mul_pd(hlf,cpsii);
                         t0r  = _mm512_mul_pd(t0,_mm512_mul_pd(cer,ir));
                         t0i  = _mm512_mul_pd(t0,_mm512_mul_pd(cei,ir)); 
                         t1   = _mm512_fmadd_pd(cpsii,cpsis,cosp);
                         t0   = _mm512_mul_pd(a3,t1);
                         t1r  = _mm512_mul_pd(EIr,t0);
                         t1i  = _mm512_mul_pd(EIi,t0);
                         cmul_zmm8r8(t0r,t0i,t1r,t1i,&resr,&resi);
                         _mm512_storeu_pd(&ESr[0], resr);
                         _mm512_storeu_pd(&ESi[0], resi);
              }


                 /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-22
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4322_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d psii,
                                            const __m512d psis,
                                            const __m512d phi) {

                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d rcs,k04,a6,t0,t1,spsii,spsis,cosp;
                          register __m512d s2psii,s2psis,cos2p,t2;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          t1    = _mm512_mul_pd(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          s2psii= _mm512_mul_pd(spsii,spsii);
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          spsis = xsinf(psis);
                          s2psis= _mm512_mul_pd(psis,psis);
                          t3    = _mm512_mul_pd(s2psii,_mm512_mul_pd(s2psis,cosp));
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4322_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d  k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d  a   = _mm512_load_pd(&pa[0]);
                          register __m512d  psii= _mm512_load_pd(&ppsii[0]);
                          register __m512d  psis= _mm512_load_pd(&ppsis[0]);
                          register __m512d  phi = _mm512_load_pd(&pphi[0]);
                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d rcs,k04,a6,t0,t1,spsii,spsis,cosp;
                          register __m512d s2psii,s2psis,cos2p,t2;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          t1    = _mm512_mul_pd(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          s2psii= _mm512_mul_pd(spsii,spsii);
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          spsis = xsinf(psis);
                          s2psis= _mm512_mul_pd(psis,psis);
                          t3    = _mm512_mul_pd(s2psii,_mm512_mul_pd(s2psis,cosp));
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4322_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  ppsii,
                                              const double * __restrict  ppsis,
                                              const double * __restrict  pphi) {

                          register __m512d  k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d  a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d  psii= _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  psis= _mm512_loadu_pd(&ppsis[0]);
                          register __m512d  phi = _mm512_loadu_pd(&pphi[0]);
                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d rcs,k04,a6,t0,t1,spsii,spsis,cosp;
                          register __m512d s2psii,s2psis,cos2p,t2;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          cos2p = _mm512_mul_pd(cosp,cosp);
                          t1    = _mm512_mul_pd(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          s2psii= _mm512_mul_pd(spsii,spsii);
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          spsis = xsinf(psis);
                          s2psis= _mm512_mul_pd(psis,psis);
                          t3    = _mm512_mul_pd(s2psii,_mm512_mul_pd(s2psis,cosp));
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,t3));
                          return (rcs);
                }


                    /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-23
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4323_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d psii,
                                            const __m512d phi) {

                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d rcs,k04,a6,t0,t1,spsii,sinp;
                          register __m512d s2psii,sin2p,t2;
                          sinp  = xsinf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          sin2p = _mm512_mul_pd(sinp,sinp);
                          t1    = _mm512_mul_pd(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          s2psii= _mm512_mul_pd(spsii,spsii);
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          t3    = _mm512_mul_pd(s2psii,sin2p);
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4323_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d  k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d  a   = _mm512_load_pd(&pa[0]);
                          register __m512d  psii= _mm512_load_pd(&ppsii[0]);
                          register __m512d  phi = _mm512_load_pd(&pphi[0]);
                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d rcs,k04,a6,t0,t1,spsii,sinp;
                          register __m512d s2psii,sin2p,t2;
                          sinp  = xsinf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          sin2p = _mm512_mul_pd(sinp,sinp);
                          t1    = _mm512_mul_pd(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          s2psii= _mm512_mul_pd(spsii,spsii);
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          t3    = _mm512_mul_pd(s2psii,sin2p);
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,t3));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4323_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  ppsii,
                                              const double * __restrict  pphi) {

                          register __m512d  k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d  a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d  psii= _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  phi = _mm512_loadu_pd(&pphi[0]);
                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          register __m512d rcs,k04,a6,t0,t1,spsii,sinp;
                          register __m512d s2psii,sin2p,t2;
                          sinp  = xsinf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          sin2p = _mm512_mul_pd(sinp,sinp);
                          t1    = _mm512_mul_pd(a,a);
                          spsii = xsinf(psii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          s2psii= _mm512_mul_pd(spsii,spsii);
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          t3    = _mm512_mul_pd(s2psii,sin2p);
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,t3));
                          return (rcs);
                }


                  /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-24
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4324_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d psis,
                                            const __m512d phi) {

                          register __m512d rcs;
                          rcs = rcs_f4323_zmm8r8(k0,a,psis,phi);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4324_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {

                           register __m512d rcs;
                           rcs = rcs_f4323_zmm8r8_a(pk0,pa,ppsis,pphi);
                           return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4324_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  ppsis,
                                              const double * __restrict  pphi) {

                           register __m512d rcs;
                           rcs = rcs_f4323_zmm8r8_u(pk0,pa,ppsis,pphi);
                           return (rcs);
                }


                  /*
                           Disc limit of cylinder (h<<a).
                           Bistatic scattering RCS for cylinder in the disc limit
                           Formula 4.3-25
                   */

  
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4325_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d psii,
                                            const __m512d psis,
                                            const __m512d phi) {

                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          const __m512d hlf    = _mm512_set1_pd(0.5f);
                          register __m512d rcs,k04,a6,t0,t1,cpsii,cpsis,cosp;
                          register __m512d t2,term;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          t1    = _mm512_mul_pd(a,a);
                          cpsii = xcosf(psii);
                          cpsii = _mm512_mul_pd(hlf,cpsii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          cpsis = xcosf(psis);
                          term  = _mm512_fmadd_pd(cpsis,cpsii,cosp);
                          term  = _mm512_mul_pd(term,term);
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,term));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4325_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d  k0  = _mm512_load_pd(&pk0[0]);
                          register __m512d  a   = _mm512_load_pd(&pa[0]);
                          register __m512d  psis= _mm512_load_pd(&ppsis[0]);
                          register __m512d  psii= _mm512_load_pd(&ppsii[0]);
                          register __m512d  phi = _mm512_load_pd(&pphi[0]);
                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          const __m512d hlf    = _mm512_set1_pd(0.5f);
                          register __m512d rcs,k04,a6,t0,t1,cpsii,cpsis,cosp;
                          register __m512d t2,term;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          t1    = _mm512_mul_pd(a,a);
                          cpsii = xcosf(psii);
                          cpsii = _mm512_mul_pd(hlf,cpsii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          cpsis = xcosf(psis);
                          term  = _mm512_fmadd_pd(cpsis,cpsii,cosp);
                          term  = _mm512_mul_pd(term,term);
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,term));
                          return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4325_zmm8r8_u(  const double * __restrict  pk0,
                                                const double * __restrict  pa,
                                                const double * __restrict  ppsis,
                                                const double * __restrict  ppsii,
                                                const double * __restrict  pphi) {

                          register __m512d  k0  = _mm512_loadu_pd(&pk0[0]);
                          register __m512d  a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d  psis= _mm512_loadu_pd(&ppsis[0]);
                          register __m512d  psii= _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  phi = _mm512_loadu_pd(&pphi[0]);
                          const __m512d _64pi9 = _mm512_set1_pd(2.263536968418066997601902412409);
                          const __m512d hlf    = _mm512_set1_pd(0.5f);
                          register __m512d rcs,k04,a6,t0,t1,cpsii,cpsis,cosp;
                          register __m512d t2,term;
                          cosp  = xcosf(phi);
                          t0    = _mm512_mul_pd(k0,k0);
                          t1    = _mm512_mul_pd(a,a);
                          cpsii = xcosf(psii);
                          cpsii = _mm512_mul_pd(hlf,cpsii);
                          k04   = _mm512_mul_pd(t0,t0);
                          t2    = _mm512_mul_pd(_64pi9,_mm512_mul_pd(k04,k04));
                          a6    = _mm512_mul_pd(t1,_mm512_mul_pd(t1,t1));
                          cpsis = xcosf(psis);
                          term  = _mm512_fmadd_pd(cpsis,cpsii,cosp);
                          term  = _mm512_mul_pd(term,term);
                          rcs   = _mm512_mul_pd(t2,_mm512_mul_pd(a6,term));
                          return (rcs);
                }


                   /*
                          Backscattering RCS for perfectly conducting wire.
                          (2*h>gamma/4)
                          Formula 4.3-29

                     */

                     /*
                          Parameter a1,a2,a3 of equation 4.3-29
                          Formula 4.3-30
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a1_f4330_zmm8r8(const __m512d k0h,
                                           const __m512d psi) {

                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d a1,_2k0h,spsi,arg,_2spsi,sarg;
                          _2k0h = _mm512_add_pd(k0h,k0h);
                          spsi  = xsinf(psi);
                          arg  = _mm512_mul_pd(_2k0h,spsi);
                          _2spsi = _mm512_add_pd(spsi,spsi); 
                          sarg   = xsinf(arg);
                          a1 = _mm512_div_pd(sarg,_2spsi);
                          return (a1); 
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a1_f4330_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d k0h = _mm512_load_pd(&pk0h[0]);
                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d a1,_2k0h,spsi,arg,_2spsi,sarg;
                          _2k0h = _mm512_add_pd(k0h,k0h);
                          spsi  = xsinf(psi);
                          arg  = _mm512_mul_pd(_2k0h,spsi);
                          _2spsi = _mm512_add_pd(spsi,spsi); 
                          sarg   = xsinf(arg);
                          a1 = _mm512_div_pd(sarg,_2spsi);
                          return (a1); 
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a1_f4330_zmm8r8_u(const double * __restrict  pk0h,
                                             const double * __restrict  ppsi) {

                          register __m512d k0h = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d a1,_2k0h,spsi,arg,_2spsi,sarg;
                          _2k0h = _mm512_add_pd(k0h,k0h);
                          spsi  = xsinf(psi);
                          arg  = _mm512_mul_pd(_2k0h,spsi);
                          _2spsi = _mm512_add_pd(spsi,spsi); 
                          sarg   = xsinf(arg);
                          a1 = _mm512_div_pd(sarg,_2spsi);
                          return (a1); 
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a2_f4330_zmm8r8(const __m512d k0h,
                                           const __m512d psi) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_sub_pd(_1,spsi);
                          arg  = _mm512_mul_pd(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_pd(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a2_f4330_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const double * __restrict __ATTR_ALIGN__(64) ppsi) {
   
                          register __m512d k0h = _mm512_load_pd(&pk0h[0]);
                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_sub_pd(_1,spsi);
                          arg  = _mm512_mul_pd(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_pd(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a2_f4330_zmm8r8_u(const double * __restrict  pk0h,
                                             const double * __restrict  ppsi) {
   
                          register __m512d k0h = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_sub_pd(_1,spsi);
                          arg  = _mm512_mul_pd(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_pd(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a3_f4330_zmm8r8(const __m512d k0h,
                                           const __m512d psi) {

                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_add_pd(_1,spsi);
                          arg  = _mm512_mul_pd(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_pd(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a3_f4330_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const double * __restrict __ATTR_ALIGN__(64) ppsi) {
   
                          register __m512d k0h = _mm512_load_pd(&pk0h[0]);
                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_add_pd(_1,spsi);
                          arg  = _mm512_mul_pd(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_pd(sarg,_1map);
                          return (a2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d a3_f4330_zmm8r8_u(const double * __restrict  pk0h,
                                             const double * __restrict  ppsi) {
   
                          register __m512d k0h = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d a2,spsi,_1msp,arg,sarg;
                          spsi = xsinf(psi);
                          _1msp= _mm512_add_pd(_1,spsi);
                          arg  = _mm512_mul_pd(k0h,_1msp);
                          sarg = xsinf(arg);
                          a2   = _mm512_div_pd(sarg,_1map);
                          return (a2);
                }


                    /*
                          Parameter F1,F2 of equation 4.3-29
                          Formula 4.3-31
                      */

                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F1_f4331_zmm8r8(const __m512d k0a) {

                          const __m512d c0 = _mm512_set1_pd(0.8905);
                          const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d n2 = _mm512_set1_pd(-2.0f);
                          register __m512d F1,om,om2,arg,larg;
                          arg = _mm512_mul_pd(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_pd(n2,larg);
                          om2 = _mm512_mul_pd(om,om);
                          F1  = _mm512_div_pd(om,_mm512_add_pd(om2,spi));
                          return (F1);
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F1_f4331_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const __m512d c0 = _mm512_set1_pd(0.8905);
                          const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d n2 = _mm512_set1_pd(-2.0f);
                          register __m512d F1,om,om2,arg,larg;
                          arg = _mm512_mul_pd(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_pd(n2,larg);
                          om2 = _mm512_mul_pd(om,om);
                          F1  = _mm512_div_pd(om,_mm512_add_pd(om2,spi));
                          return (F1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F1_f4331_zmm8r8_u(const double * __restrict  pk0a) {

                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d c0 = _mm512_set1_pd(0.8905);
                          const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d n2 = _mm512_set1_pd(-2.0f);
                          register __m512d F1,om,om2,arg,larg;
                          arg = _mm512_mul_pd(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_pd(n2,larg);
                          om2 = _mm512_mul_pd(om,om);
                          F1  = _mm512_div_pd(om,_mm512_add_pd(om2,spi));
                          return (F1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F2_f4331_zmm8r8(const __m512d k0a) {

                          const __m512d c0 = _mm512_set1_pd(0.8905);
                          const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d n2 = _mm512_set1_pd(-2.0f);
                          register __m512d F1,om,om2,arg,larg;
                          arg = _mm512_mul_pd(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_pd(n2,larg);
                          om2 = _mm512_mul_pd(om,om);
                          F1  = _mm512_div_pd(PI,_mm512_add_pd(om2,spi));
                          return (F1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F2_f4331_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const __m512d c0 = _mm512_set1_pd(0.8905);
                          const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d n2 = _mm512_set1_pd(-2.0f);
                          register __m512d F1,om,om2,arg,larg;
                          arg = _mm512_mul_pd(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_pd(n2,larg);
                          om2 = _mm512_mul_pd(om,om);
                          F1  = _mm512_div_pd(PI,_mm512_add_pd(om2,spi));
                          return (F1);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F2_f4331_zmm8r8_u(const double * __restrict  pk0a) {

                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d c0 = _mm512_set1_pd(0.8905);
                          const __m512d spi= _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d n2 = _mm512_set1_pd(-2.0f);
                          register __m512d F1,om,om2,arg,larg;
                          arg = _mm512_mul_pd(k0a,c0);
                          larg= xlogf(arg);
                          om  = _mm512_mul_pd(n2,larg);
                          om2 = _mm512_mul_pd(om,om);
                          F1  = _mm512_div_pd(PI,_mm512_add_pd(om2,spi));
                          return (F1);
                }


                     /*
                          Parameter (helper) Lambda of equation 4.3-29
                          Formula 4.3-34
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d L_f4334_zmm8r8(const __m512d k0h,
                                          const __m512d k0a) {

                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d c0  = _mm512_set1_pd(4.11);
                          const __m512d hlf = _mm512_set1_pd(-0.5f);
                          const __m512d n2  = _mm512_set1_pd(-2.0f);
                          const __m512d c1  = _mm512_set1_pd(0.8905);
                          const __m512d _0  = _mm512_setzero_pd();
                          register __m512d L,om,del,ck0h,sk0h,t0;
                          register __m512d ar1,ar2,lar1,lar2;
                          ar1  = _mm512_mul_pd(k0a,c1);
                          lar1 =xlogf(ar1);
                          ar2  = _mm512_div_pd(k0h,c0);
                          lar2 = xlogf(ar2);
                          om   = _mm512_mul_pd(n2,lar1);
                          del  = _mm512_mul_pd(hlf,lar2);
                          ck0h = xcosf(k0h);
                          t0   = _mm512_sub_pd(_0,_mm512sub_pd(om,del));
                          sk0h = xsinf(k0h);
                          L    = _mm512_fmadd_pd(pi4,sk0h,_mm512_mul_pd(ck0h,t0));
                          return (L);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d L_f4334_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d k0h = _mm512_load_pd(&pk0h[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d c0  = _mm512_set1_pd(4.11);
                          const __m512d hlf = _mm512_set1_pd(-0.5f);
                          const __m512d n2  = _mm512_set1_pd(-2.0f);
                          const __m512d c1  = _mm512_set1_pd(0.8905);
                          const __m512d _0  = _mm512_setzero_pd();
                          register __m512d L,om,del,ck0h,sk0h,t0;
                          register __m512d ar1,ar2,lar1,lar2;
                          ar1  = _mm512_mul_pd(k0a,c1);
                          lar1 =xlogf(ar1);
                          ar2  = _mm512_div_pd(k0h,c0);
                          lar2 = xlogf(ar2);
                          om   = _mm512_mul_pd(n2,lar1);
                          del  = _mm512_mul_pd(hlf,lar2);
                          ck0h = xcosf(k0h);
                          t0   = _mm512_sub_pd(_0,_mm512sub_pd(om,del));
                          sk0h = xsinf(k0h);
                          L    = _mm512_fmadd_pd(pi4,sk0h,_mm512_mul_pd(ck0h,t0));
                          return (L);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d L_f4334_zmm8r8_u(const double * __restrict pk0h,
                                            const double * __restrict pk0a) {

                          register __m512d k0h = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d c0  = _mm512_set1_pd(4.11);
                          const __m512d hlf = _mm512_set1_pd(-0.5f);
                          const __m512d n2  = _mm512_set1_pd(-2.0f);
                          const __m512d c1  = _mm512_set1_pd(0.8905);
                          const __m512d _0  = _mm512_setzero_pd();
                          register __m512d L,om,del,ck0h,sk0h,t0;
                          register __m512d ar1,ar2,lar1,lar2;
                          ar1  = _mm512_mul_pd(k0a,c1);
                          lar1 =xlogf(ar1);
                          ar2  = _mm512_div_pd(k0h,c0);
                          lar2 = xlogf(ar2);
                          om   = _mm512_mul_pd(n2,lar1);
                          del  = _mm512_mul_pd(hlf,lar2);
                          ck0h = xcosf(k0h);
                          t0   = _mm512_sub_pd(_0,_mm512sub_pd(om,del));
                          sk0h = xsinf(k0h);
                          L    = _mm512_fmadd_pd(pi4,sk0h,_mm512_mul_pd(ck0h,t0));
                          return (L);
                }


                   /*
                          Parameter (helper) Sigma of equation 4.3-29
                          Formula 4.3-35
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d S_f4335_zmm8r8(const __m512d k0a,
                                          const __m512d k0h) {

                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(7.12);
                          register __m512d ar,lar,sk0h,ck0h;
                          register __m512d S,t0;
                          ar  = _mm512_mul_pd(c0,k0a);
                          lar = xlogf(ar);
                          sk0h= xsinf(k0h);
                          ck0h= xcosf(k0h);
                          t0  = _mm512_mul_pd(hlf,lar);
                          S   = _mm512_fmsub_pd(t0,sk0h,
                                            _mm512_mul_pd(pi4,ck0h));
                          return (S);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d S_f4335_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0h) {

                          register __m512d k0h = _mm512_load_pd(&pk0h[0]);
                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(7.12);
                          register __m512d ar,lar,sk0h,ck0h;
                          register __m512d S,t0;
                          ar  = _mm512_mul_pd(c0,k0a);
                          lar = xlogf(ar);
                          sk0h= xsinf(k0h);
                          ck0h= xcosf(k0h);
                          t0  = _mm512_mul_pd(hlf,lar);
                          S   = _mm512_fmsub_pd(t0,sk0h,
                                            _mm512_mul_pd(pi4,ck0h));
                          return (S);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d S_f4335_zmm8r8_u(const double * __restrict  pk0a,
                                            const double * __restrict  pk0h) {

                          register __m512d k0h = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(7.12);
                          register __m512d ar,lar,sk0h,ck0h;
                          register __m512d S,t0;
                          ar  = _mm512_mul_pd(c0,k0a);
                          lar = xlogf(ar);
                          sk0h= xsinf(k0h);
                          ck0h= xcosf(k0h);
                          t0  = _mm512_mul_pd(hlf,lar);
                          S   = _mm512_fmsub_pd(t0,sk0h,
                                            _mm512_mul_pd(pi4,ck0h));
                          return (S);
                }


                  /*

                           Parameter G1,G2 of equation 4.3-29
                           Formula 4.3-32
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G2_f4332_zmm8r8(const __m512d k0h,
                                           const __m512d k0a) {

                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d G2,L,S,num,den;
                          L = L_f4334_zmm8r8(k0h,k0a);
                          S = S_f4335_zmm8r8(k0a,k0h);
                          num = _mm512_mul_pd(hlf,S);
                          den = _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          G2  = _mm512_div_pd(num,den);
                          return (G2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G2_f4332_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d G2,L,S,num,den;
                          L = L_f4334_zmm8r8_a(pk0h,pk0a);
                          S = S_f4335_zmm8r8_a(pk0a,pk0h);
                          num = _mm512_mul_pd(hlf,S);
                          den = _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          G2  = _mm512_div_pd(num,den);
                          return (G2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G2_f4332_zmm8r8_u(const double * __restrict  pk0h,
                                             const double * __restrict  pk0a) {

                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d G2,L,S,num,den;
                          L = L_f4334_zmm8r8_a(pk0h,pk0a);
                          S = S_f4335_zmm8r8_a(pk0a,pk0h);
                          num = _mm512_mul_pd(hlf,S);
                          den = _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          G2  = _mm512_div_pd(num,den);
                          return (G2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G1_f4332_zmm8r8(const __m512d k0h,
                                           const __m512d k0a) {

                         const __m512d hlf = _mm512_set1_pd(0.5f);
                         const __m512d n2  = _mm512_set1_pd(-2.0f);
                         const __m512d c0  = _mm512_set1_pd(0.8905);
                         register __m512d G1,L,S,om,G1,ln,num,den,om2,t0,rat;
                         L = L_f4334_zmm8r8(k0h,k0a);
                         S = S_f4335_zmm8r8(k0a,k0h);
                         ln= xlogf(_mm512_mul_pd(k0a,c0));
                         om= _mm512_mul_pd(m2,ln); 
                         G2= G2_f4332_zmm8r8(k0h,k0a);
                         om2= _mm512_add_pd(om,om);
                         num= _mm512_mul_pd(hlf,L);
                         t0 = _mm512_mul_pd(PI,G2);
                         ln = _mm512_div_pd(t0,om2);
                         den= _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                         rat= _mm512_div_pd(num,den);
                         G1 = _mm512_sub_pd(rat,ln);
                         return (G1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G1_f4332_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                         register __m512d k0h  = _mm512_load_pd(&pk0h[0]);
                         register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                         const __m512d hlf = _mm512_set1_pd(0.5f);
                         const __m512d n2  = _mm512_set1_pd(-2.0f);
                         const __m512d c0  = _mm512_set1_pd(0.8905);
                         register __m512d G1,L,S,om,G1,ln,num,den,om2,t0,rat;
                         L = L_f4334_zmm8r8(k0h,k0a);
                         S = S_f4335_zmm8r8(k0a,k0h);
                         ln= xlogf(_mm512_mul_pd(k0a,c0));
                         om= _mm512_mul_pd(m2,ln); 
                         G2= G2_f4332_zmm8r8(k0h,k0a);
                         om2= _mm512_add_pd(om,om);
                         num= _mm512_mul_pd(hlf,L);
                         t0 = _mm512_mul_pd(PI,G2);
                         ln = _mm512_div_pd(t0,om2);
                         den= _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                         rat= _mm512_div_pd(num,den);
                         G1 = _mm512_sub_pd(rat,ln);
                         return (G1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G1_f4332_zmm8r8_u(const double * __restrict  pk0h,
                                             const double * __restrict  pk0a) {

                         register __m512d k0h  = _mm512_loadu_pd(&pk0h[0]);
                         register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                         const __m512d hlf = _mm512_set1_pd(0.5f);
                         const __m512d n2  = _mm512_set1_pd(-2.0f);
                         const __m512d c0  = _mm512_set1_pd(0.8905);
                         register __m512d G1,L,S,om,G1,ln,num,den,om2,t0,rat;
                         L = L_f4334_zmm8r8(k0h,k0a);
                         S = S_f4335_zmm8r8(k0a,k0h);
                         ln= xlogf(_mm512_mul_pd(k0a,c0));
                         om= _mm512_mul_pd(m2,ln); 
                         G2= G2_f4332_zmm8r8(k0h,k0a);
                         om2= _mm512_add_pd(om,om);
                         num= _mm512_mul_pd(hlf,L);
                         t0 = _mm512_mul_pd(PI,G2);
                         ln = _mm512_div_pd(t0,om2);
                         den= _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                         rat= _mm512_div_pd(num,den);
                         G1 = _mm512_sub_pd(rat,ln);
                         return (G1);
                 }


                     /*

                           Parameter H1,H2 of equation 4.3-29
                           Formula 4.3-33
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d H2_f4333_zmm8r8(const __m512d k0h,
                                           const __m512d k0a) {

                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d H2,L,S,num,den,arg;
                          arg  = _mm512_mul_pd(pi2,k0h);
                          L    = L_f4334_zmm8r8(k0h,k0a);
                          S    = S_f4335_zmm8r8(k0a,k0h);
                          num  = _mm512_mul_pd(hlf,S);
                          den  = _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          H2  = _mm512_div_pd(num,den);
                          return (H2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d H2_f4333_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                             const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          register __m512d k0h  = _mm512_load_pd(&pk0h[0]);
                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d H2,L,S,num,den,arg;
                          arg  = _mm512_mul_pd(pi2,k0h);
                          L    = L_f4334_zmm8r8(k0h,k0a);
                          S    = S_f4335_zmm8r8(k0a,k0h);
                          num  = _mm512_mul_pd(hlf,S);
                          den  = _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          H2  = _mm512_div_pd(num,den);
                          return (H2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d H2_f4333_zmm8r8_u(const double * __restrict  pk0h,
                                             const double * __restrict  pk0a) {

                          register __m512d k0h  = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          register __m512d H2,L,S,num,den,arg;
                          arg  = _mm512_mul_pd(pi2,k0h);
                          L    = L_f4334_zmm8r8(k0h,k0a);
                          S    = S_f4335_zmm8r8(k0a,k0h);
                          num  = _mm512_mul_pd(hlf,S);
                          den  = _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          H2  = _mm512_div_pd(num,den);
                          return (H2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d H1_f4333_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(64) pk0h) {

                          register __m512d k0h  = _mm512_load_pd(&pk0h[0]);
                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d n2  = _mm512_set1_pd(-2.0f);
                          const __m512d c0  = _mm512_set1_pd(0.8905f);
                          register __m512d H1,H2,om,ar,lar,L,S,num,den;
                          register __m512d om2,t0,arg;
                          ar = _mm512_mul_pd(k0a,c0);
                          arg= _mm512_mul_pd(k0h,pi2);
                          lar= xlogf(ar);
                          om = _mm512_mul_pd(n2,lar);
                          L  = L_f4334_zmm8r8(k0h,k0a);
                          om2= _mm512_add_pd(om,om);
                          S  = S_f4335_zmm8r8(k0a,k0h);
                          H1 = H1_f4333_zmm8r8(k0h,k0a);
                          num= _mm512_mul_pd(hlf,L);
                          t0 = _mm512_div_pd(_mm512_mul_pd(PI,H1),om2);
                          den= _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          ar = _mm512_div_pd(num,den);
                          H1 = _mm512_sub_pd(ar,t0);
                          return (H1);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d H1_f4333_zmm8r8_u(const double * __restrict  pk0a,
                                             const double * __restrict  pk0h) {

                          register __m512d k0h  = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d n2  = _mm512_set1_pd(-2.0f);
                          const __m512d c0  = _mm512_set1_pd(0.8905f);
                          register __m512d H1,H2,om,ar,lar,L,S,num,den;
                          register __m512d om2,t0,arg;
                          ar = _mm512_mul_pd(k0a,c0);
                          arg= _mm512_mul_pd(k0h,pi2);
                          lar= xlogf(ar);
                          om = _mm512_mul_pd(n2,lar);
                          L  = L_f4334_zmm8r8(k0h,k0a);
                          om2= _mm512_add_pd(om,om);
                          S  = S_f4335_zmm8r8(k0a,k0h);
                          H1 = H1_f4333_zmm8r8(k0h,k0a);
                          num= _mm512_mul_pd(hlf,L);
                          t0 = _mm512_div_pd(_mm512_mul_pd(PI,H1),om2);
                          den= _mm512_fmadd_pd(L,L,_mm512_mul_pd(S,S));
                          ar = _mm512_div_pd(num,den);
                          H1 = _mm512_sub_pd(ar,t0);
                          return (H1);
               }


                 /*
                          Backscattering RCS for perfectly conducting wire.
                          (2*h>gamma/4)
                          Formula 4.3-29

                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4329_zmm8r8(const __m512d k0,
                                            const __m512d gami,
                                            const __m512d gams,
                                            const __m512d k0h,
                                            const __m512d k0a,
                                            const __m512d psi) {

                          const __m512d _16pi = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _2    = _mm512_set1_pd(2.0f);
                          register __m512d rcs,a1,a2,a3,F1,F2,G1,G2,H1,H2,first;
                          register __m512d cgami,cgams,c2gami,c2gams,sinps;
                          register __m512d arg,sarg,carg,t0,t1,t2,t3,t4,x0,x1,t5,b0;
                          register __m512d a1s,F1F2,G1G2,a2pa3,a2ma3,H1H2,a2sma3s;
                          register __m512d GHGH,_2a1,FGFG,FHFH,tmp1,tmp2,tmp3;
                          b0     = _mm512_div_pd(_16pi,_mm512_mul_pd(k0,k0));
                          a1     = a1_f4330_zmm8r8(k0h,psi);
                          _2a1   = _mm512_add_pd(a1,a1);
                          cgami  = xcosf(gami);
                          F1     = F1_f4331_zmm8r8(k0a);
                          c2gami = _mm512_mul_pd(cgami,cgami);
                          F2     = F2_f4331_zmm8r8(k0a);
                          cgams  = xcosf(gams);
                          G1     = G1_f4332_zmm8r8(k0h,k0a);
                          c2gams = _mm512_mul_pd(cgams,cgams);
                          a2     = a2_f4330_zmm8r8(k0h,psi);
                          first  = _mm512_mul_pd(b0,_mm512_mul_pd(c2gami,c2gams));
                          G2     = G1_f4332_zmm8r8(k0h,k0a);
                          sinps  = xsinf(psi);
                          a3     = a3_f4330_zmm8r8(k0h,psi);
                          H1     = H1_f4333_zmm8r8(k0h,k0a);
                          arg    = _mm512_mul_pd(k0h,sinps);
                          H2     = H2_f4333_zmm8r8(k0h,k0a);
                          sarg   = xsinf(arg);
                          a1s    = _mm512_mul_pd(a1,a1);
                          carg   = xcosf(arg);
                          x0     = _mm512_add_pd(a2,a3);
                          a2pa3  = _mm512_mul_pd(x0,x0);
                          F1F2   = _mm512_fmadd_pd(F1,F1,_mm512_mul_pd(F2,F2));
                          x1     = _mm512_sub_pd(a2,a3);
                          t0     = _mm512_mul_pd(a1s,F1F2);
                          a2ma3  = _mm512_mul_pd(x1,x1);
                          G1G2   = _mm512_fmadd_pd(G1,G1,_mm512_mul_pd(G2,G2));
                          t1     = _mm512_mul_pd(a2pa3,_mm512_mul_pd(G1G2,carg));
                          x0     = _mm512_mul_pd(sarg,sarg);
                          H1H2   = _mm512_fmadd_pd(H1,H1,_mm512_mul_pd(H2,H2));
                          t2     = _mm512_mul_pd(a2ma3,_mm512_mul_pd(H1H2,x0));
                          a2sma3s= _mm512_mul_pd(_2,_mm512_fmsub_pd(a2,a2,
                                                                _mm512_mul_pd(a3,a3)));
                          GHGH   = _mm512_fmadd_pd(G1,H1,_mm512_mul_pd(G2,H2));
                          x1     = _mm512_mul_pd(carg,sarg);
                          t3     = _mm512_mul_pd(a2sma3s,_mm512_mul_pd(GHGH,x1));
                          x0     = _mm512_mul_pd(_2a1,a2pa3);
                          FGFG   = _mm512_fmadd_pd(F1,G1,_mm512_mul_pd(F2,G2));
                          t4     = _mm512_mul_pd(x0,_mm512_mul_pd(FGFG,carg);
                          x1     = _mm512_mul_pd(_2a1,a2ma3);
                          FHFH   = _mm512_fmadd_pd(F1,H1,_mm512_mul_pd(F2,H2));
                          t5     = _mm512_mul_pd(x1,_mm512_mul_pd(FHFH,sarg));
                          tmp1   = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
                          tmp2   = _mm512_sub_pd(_mm512_add_pd(t3,t4),t5);
                          tmp3   = _mm512_sub_pd(tmp1,tmp2);
                          rcs    = _mm512_mul_pd(first,tmp3);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4329_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pgami,
                                              const double * __restrict __ATTR_ALIGN__(64) pgams,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0h,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsi) {


                          register __m512d k0     = _mm512_load_pd(&pk0[0]);
                          register __m512d gami   = _mm512_load_pd(&pgami[0]);
                          register __m512d gams   = _mm512_load_pd(&pgams[0]);
                          register __m512d k0h    = _mm512_load_pd(&pk0h[0]);
                          register __m512d k0a    = _mm512_load_pd(&pk0a[0]);
                          register __m512d psi    = _mm512_load_pd(&ppsi[0]);  
                          const __m512d _16pi = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _2    = _mm512_set1_pd(2.0f);
                          register __m512d rcs,a1,a2,a3,F1,F2,G1,G2,H1,H2,first;
                          register __m512d cgami,cgams,c2gami,c2gams,sinps;
                          register __m512d arg,sarg,carg,t0,t1,t2,t3,t4,x0,x1,t5,b0;
                          register __m512d a1s,F1F2,G1G2,a2pa3,a2ma3,H1H2,a2sma3s;
                          register __m512d GHGH,_2a1,FGFG,FHFH,tmp1,tmp2,tmp3;
                          b0     = _mm512_div_pd(_16pi,_mm512_mul_pd(k0,k0));
                          a1     = a1_f4330_zmm8r8(k0h,psi);
                          _2a1   = _mm512_add_pd(a1,a1);
                          cgami  = xcosf(gami);
                          F1     = F1_f4331_zmm8r8(k0a);
                          c2gami = _mm512_mul_pd(cgami,cgami);
                          F2     = F2_f4331_zmm8r8(k0a);
                          cgams  = xcosf(gams);
                          G1     = G1_f4332_zmm8r8(k0h,k0a);
                          c2gams = _mm512_mul_pd(cgams,cgams);
                          a2     = a2_f4330_zmm8r8(k0h,psi);
                          first  = _mm512_mul_pd(b0,_mm512_mul_pd(c2gami,c2gams));
                          G2     = G1_f4332_zmm8r8(k0h,k0a);
                          sinps  = xsinf(psi);
                          a3     = a3_f4330_zmm8r8(k0h,psi);
                          H1     = H1_f4333_zmm8r8(k0h,k0a);
                          arg    = _mm512_mul_pd(k0h,sinps);
                          H2     = H2_f4333_zmm8r8(k0h,k0a);
                          sarg   = xsinf(arg);
                          a1s    = _mm512_mul_pd(a1,a1);
                          carg   = xcosf(arg);
                          x0     = _mm512_add_pd(a2,a3);
                          a2pa3  = _mm512_mul_pd(x0,x0);
                          F1F2   = _mm512_fmadd_pd(F1,F1,_mm512_mul_pd(F2,F2));
                          x1     = _mm512_sub_pd(a2,a3);
                          t0     = _mm512_mul_pd(a1s,F1F2);
                          a2ma3  = _mm512_mul_pd(x1,x1);
                          G1G2   = _mm512_fmadd_pd(G1,G1,_mm512_mul_pd(G2,G2));
                          t1     = _mm512_mul_pd(a2pa3,_mm512_mul_pd(G1G2,carg));
                          x0     = _mm512_mul_pd(sarg,sarg);
                          H1H2   = _mm512_fmadd_pd(H1,H1,_mm512_mul_pd(H2,H2));
                          t2     = _mm512_mul_pd(a2ma3,_mm512_mul_pd(H1H2,x0));
                          a2sma3s= _mm512_mul_pd(_2,_mm512_fmsub_pd(a2,a2,
                                                                _mm512_mul_pd(a3,a3)));
                          GHGH   = _mm512_fmadd_pd(G1,H1,_mm512_mul_pd(G2,H2));
                          x1     = _mm512_mul_pd(carg,sarg);
                          t3     = _mm512_mul_pd(a2sma3s,_mm512_mul_pd(GHGH,x1));
                          x0     = _mm512_mul_pd(_2a1,a2pa3);
                          FGFG   = _mm512_fmadd_pd(F1,G1,_mm512_mul_pd(F2,G2));
                          t4     = _mm512_mul_pd(x0,_mm512_mul_pd(FGFG,carg);
                          x1     = _mm512_mul_pd(_2a1,a2ma3);
                          FHFH   = _mm512_fmadd_pd(F1,H1,_mm512_mul_pd(F2,H2));
                          t5     = _mm512_mul_pd(x1,_mm512_mul_pd(FHFH,sarg));
                          tmp1   = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
                          tmp2   = _mm512_sub_pd(_mm512_add_pd(t3,t4),t5);
                          tmp3   = _mm512_sub_pd(tmp1,tmp2);
                          rcs    = _mm512_mul_pd(first,tmp3);
                          return (rcs);
               }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4329_zmm8r8_u(const double * __restrict pk0,
                                              const double * __restrict  pgami,
                                              const double * __restrict  pgams,
                                              const double * __restrict  pk0h,
                                              const double * __restrict  pk0a,
                                              const double * __restrict  ppsi) {


                          register __m512d k0     = _mm512_loadu_pd(&pk0[0]);
                          register __m512d gami   = _mm512_loadu_pd(&pgami[0]);
                          register __m512d gams   = _mm512_loadu_pd(&pgams[0]);
                          register __m512d k0h    = _mm512_loadu_pd(&pk0h[0]);
                          register __m512d k0a    = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d psi    = _mm512_loadu_pd(&ppsi[0]);  
                          const __m512d _16pi = _mm512_set1_pd(50.265482457436691815402294132472);
                          const __m512d _2    = _mm512_set1_pd(2.0f);
                          register __m512d rcs,a1,a2,a3,F1,F2,G1,G2,H1,H2,first;
                          register __m512d cgami,cgams,c2gami,c2gams,sinps;
                          register __m512d arg,sarg,carg,t0,t1,t2,t3,t4,x0,x1,t5,b0;
                          register __m512d a1s,F1F2,G1G2,a2pa3,a2ma3,H1H2,a2sma3s;
                          register __m512d GHGH,_2a1,FGFG,FHFH,tmp1,tmp2,tmp3;
                          b0     = _mm512_div_pd(_16pi,_mm512_mul_pd(k0,k0));
                          a1     = a1_f4330_zmm8r8(k0h,psi);
                          _2a1   = _mm512_add_pd(a1,a1);
                          cgami  = xcosf(gami);
                          F1     = F1_f4331_zmm8r8(k0a);
                          c2gami = _mm512_mul_pd(cgami,cgami);
                          F2     = F2_f4331_zmm8r8(k0a);
                          cgams  = xcosf(gams);
                          G1     = G1_f4332_zmm8r8(k0h,k0a);
                          c2gams = _mm512_mul_pd(cgams,cgams);
                          a2     = a2_f4330_zmm8r8(k0h,psi);
                          first  = _mm512_mul_pd(b0,_mm512_mul_pd(c2gami,c2gams));
                          G2     = G1_f4332_zmm8r8(k0h,k0a);
                          sinps  = xsinf(psi);
                          a3     = a3_f4330_zmm8r8(k0h,psi);
                          H1     = H1_f4333_zmm8r8(k0h,k0a);
                          arg    = _mm512_mul_pd(k0h,sinps);
                          H2     = H2_f4333_zmm8r8(k0h,k0a);
                          sarg   = xsinf(arg);
                          a1s    = _mm512_mul_pd(a1,a1);
                          carg   = xcosf(arg);
                          x0     = _mm512_add_pd(a2,a3);
                          a2pa3  = _mm512_mul_pd(x0,x0);
                          F1F2   = _mm512_fmadd_pd(F1,F1,_mm512_mul_pd(F2,F2));
                          x1     = _mm512_sub_pd(a2,a3);
                          t0     = _mm512_mul_pd(a1s,F1F2);
                          a2ma3  = _mm512_mul_pd(x1,x1);
                          G1G2   = _mm512_fmadd_pd(G1,G1,_mm512_mul_pd(G2,G2));
                          t1     = _mm512_mul_pd(a2pa3,_mm512_mul_pd(G1G2,carg));
                          x0     = _mm512_mul_pd(sarg,sarg);
                          H1H2   = _mm512_fmadd_pd(H1,H1,_mm512_mul_pd(H2,H2));
                          t2     = _mm512_mul_pd(a2ma3,_mm512_mul_pd(H1H2,x0));
                          a2sma3s= _mm512_mul_pd(_2,_mm512_fmsub_pd(a2,a2,
                                                                _mm512_mul_pd(a3,a3)));
                          GHGH   = _mm512_fmadd_pd(G1,H1,_mm512_mul_pd(G2,H2));
                          x1     = _mm512_mul_pd(carg,sarg);
                          t3     = _mm512_mul_pd(a2sma3s,_mm512_mul_pd(GHGH,x1));
                          x0     = _mm512_mul_pd(_2a1,a2pa3);
                          FGFG   = _mm512_fmadd_pd(F1,G1,_mm512_mul_pd(F2,G2));
                          t4     = _mm512_mul_pd(x0,_mm512_mul_pd(FGFG,carg);
                          x1     = _mm512_mul_pd(_2a1,a2ma3);
                          FHFH   = _mm512_fmadd_pd(F1,H1,_mm512_mul_pd(F2,H2));
                          t5     = _mm512_mul_pd(x1,_mm512_mul_pd(FHFH,sarg));
                          tmp1   = _mm512_add_pd(t0,_mm512_add_pd(t1,t2));
                          tmp2   = _mm512_sub_pd(_mm512_add_pd(t3,t4),t5);
                          tmp3   = _mm512_sub_pd(tmp1,tmp2);
                          rcs    = _mm512_mul_pd(first,tmp3);
                          return (rcs);
               }


                  /*

                         Simplified back and bistatic scattering RCS for
                         half and full-wave dipole (2*h == gam0/2, and gam0)
                         gam0 -- wavelength.
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4337_zmm8r8(const __m512d gammi,
                                            const __m512d gamms,
                                            const __m512d psii,
                                            const __m512d psis,
                                            const __m512d g0 )  {//wavelength coeff

                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          register __m512d rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512d spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_pd(pi2,spsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_pd(pi2,spsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_pd(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_pd(cgami,cgami);
                          t0    = _mm512_mul_pd(g0,_mm512_mul_pd(c2gami,c2gams));
                          c1    = xcosf(carg1);
                          rat1  = _mm512_div_pd(c1,cpsii);
                          tmp0  = _mm512_mul_pd(rat1,rat1);
                          c2    = xcosf(carg2);
                          rat2  = _mm512_div_pd(c2,cpsis);
                          tmp1  = _mm512_mul_pd(rat2,rat2);
                          t1    = _mm512_mul_pd(tmp0,tmp1);
                          rcs   = _mm512_mul_pd(t0,t1);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4337_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgammi,
                                              const double * __restrict __ATTR_ALIGN__(64) pgamms,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64) pg0 )  {//wavelength coeff

                          register __m512d  gammi = _mm512_load_pd(&pgammi[0]);
                          register __m512d  gamms = _mm512_load_pd(&pgamms[0]);
                          register __m512d  psii  = _mm512_load_pd(&ppsii[0]);
                          register __m512d  psis  = _mm512_load_pd(&ppsis[0]);
                          register __m512d  g0    = _mm512_load_pd(&pg0[0]);
                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          register __m512d rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512d spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_pd(pi2,spsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_pd(pi2,spsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_pd(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_pd(cgami,cgami);
                          t0    = _mm512_mul_pd(g0,_mm512_mul_pd(c2gami,c2gams));
                          c1    = xcosf(carg1);
                          rat1  = _mm512_div_pd(c1,cpsii);
                          tmp0  = _mm512_mul_pd(rat1,rat1);
                          c2    = xcosf(carg2);
                          rat2  = _mm512_div_pd(c2,cpsis);
                          tmp1  = _mm512_mul_pd(rat2,rat2);
                          t1    = _mm512_mul_pd(tmp0,tmp1);
                          rcs   = _mm512_mul_pd(t0,t1);
                          return (rcs);
                 }
      


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4337_zmm8r8_u(const double * __restrict  pgammi,
                                              const double * __restrict  pgamms,
                                              const double * __restrict  ppsii,
                                              const double * __restrict  ppsis,
                                              const double * __restrict  pg0 )  {//wavelength coeff

                          register __m512d  gammi = _mm512_loadu_pd(&pgammi[0]);
                          register __m512d  gamms = _mm512_loadu_pd(&pgamms[0]);
                          register __m512d  psii  = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  psis  = _mm512_loadu_pd(&ppsis[0]);
                          register __m512d  g0    = _mm512_loadu_pd(&pg0[0]);
                          const __m512d pi2 = _mm512_set1_pd(1.57079632679489661923132169164);
                          register __m512d rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512d spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_pd(pi2,spsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_pd(pi2,spsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_pd(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_pd(cgami,cgami);
                          t0    = _mm512_mul_pd(g0,_mm512_mul_pd(c2gami,c2gams));
                          c1    = xcosf(carg1);
                          rat1  = _mm512_div_pd(c1,cpsii);
                          tmp0  = _mm512_mul_pd(rat1,rat1);
                          c2    = xcosf(carg2);
                          rat2  = _mm512_div_pd(c2,cpsis);
                          tmp1  = _mm512_mul_pd(rat2,rat2);
                          t1    = _mm512_mul_pd(tmp0,tmp1);
                          rcs   = _mm512_mul_pd(t0,t1);
                          return (rcs);
                 }


                   /*

                         Simplified back and bistatic scattering RCS for
                         Full-wave dipole (2*h == gam0)
                         gam0 -- wavelength.
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4340_zmm8r8(const __m512d gammi,
                                            const __m512d gamms,
                                            const __m512d psii,
                                            const __m512d psis,
                                            const __m512d g0 )  {//wavelength coeff

                         
                          register __m512d rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512d spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_pd(PI,spsii);
                          cpsii = _mm512_mul_pd(cpsii,cpsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_pd(PI,spsis);
                          cpsis = _mm512_mul_pd(cpsis,cpsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_pd(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_pd(cgami,cgami);
                          t0    = _mm512_mul_pd(g0,_mm512_mul_pd(c2gami,c2gams));
                          c1    = xsinf(carg1);
                          rat1  = _mm512_div_pd(c1,cpsii);
                          tmp0  = _mm512_mul_pd(rat1,rat1);
                          c2    = xsinf(carg2);
                          rat2  = _mm512_div_pd(c2,cpsis);
                          tmp1  = _mm512_mul_pd(rat2,rat2);
                          t1    = _mm512_mul_pd(tmp0,tmp1);
                          rcs   = _mm512_mul_pd(t0,t1);
                          return (rcs);
                 }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4340_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgammi,
                                              const double * __restrict __ATTR_ALIGN__(64) pgamms,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64) pg0 )  {//wavelength coeff

                         
                          register __m512d  gammi = _mm512_load_pd(&pgammi[0]);
                          register __m512d  gamms = _mm512_load_pd(&pgamms[0]);
                          register __m512d  psii  = _mm512_load_pd(&ppsii[0]);
                          register __m512d  psis  = _mm512_load_pd(&ppsis[0]);
                          register __m512d  g0    = _mm512_load_pd(&pg0[0]);
                          register __m512d rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512d spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_pd(PI,spsii);
                          cpsii = _mm512_mul_pd(cpsii,cpsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_pd(PI,spsis);
                          cpsis = _mm512_mul_pd(cpsis,cpsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_pd(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_pd(cgami,cgami);
                          t0    = _mm512_mul_pd(g0,_mm512_mul_pd(c2gami,c2gams));
                          c1    = xsinf(carg1);
                          rat1  = _mm512_div_pd(c1,cpsii);
                          tmp0  = _mm512_mul_pd(rat1,rat1);
                          c2    = xsinf(carg2);
                          rat2  = _mm512_div_pd(c2,cpsis);
                          tmp1  = _mm512_mul_pd(rat2,rat2);
                          t1    = _mm512_mul_pd(tmp0,tmp1);
                          rcs   = _mm512_mul_pd(t0,t1);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4340_zmm8r8_u(const double * __restrict  pgammi,
                                              const double * __restrict  pgamms,
                                              const double * __restrict  ppsii,
                                              const double * __restrict  ppsis,
                                              const double * __restrict  pg0 )  {//wavelength coeff

                         
                          register __m512d  gammi = _mm512_loadu_pd(&pgammi[0]);
                          register __m512d  gamms = _mm512_loadu_pd(&pgamms[0]);
                          register __m512d  psii  = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  psis  = _mm512_loadu_pd(&ppsis[0]);
                          register __m512d  g0    = _mm512_loadu_pd(&pg0[0]);
                          register __m512d rcs,cgami,cgams,c2gami,c2gams,t0,carg1,carg2;
                          register __m512d spsii,spsis,cpsii,cpsis,rat1,rat2,t1,c1,c2,tmp0,tmp1;
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          cpsii = xcosf(psii);
                          carg1 = _mm512_mul_pd(PI,spsii);
                          cpsii = _mm512_mul_pd(cpsii,cpsii);
                          cpsis = xcosf(psis);
                          carg2 = _mm512_mul_pd(PI,spsis);
                          cpsis = _mm512_mul_pd(cpsis,cpsis);
                          cgams = xcosf(gamms);
                          c2gams= _mm512_mul_pd(cgams,cgams);
                          cgami = xcosf(gammi);
                          c2gami= _mm512_mul_pd(cgami,cgami);
                          t0    = _mm512_mul_pd(g0,_mm512_mul_pd(c2gami,c2gams));
                          c1    = xsinf(carg1);
                          rat1  = _mm512_div_pd(c1,cpsii);
                          tmp0  = _mm512_mul_pd(rat1,rat1);
                          c2    = xsinf(carg2);
                          rat2  = _mm512_div_pd(c2,cpsis);
                          tmp1  = _mm512_mul_pd(rat2,rat2);
                          t1    = _mm512_mul_pd(tmp0,tmp1);
                          rcs   = _mm512_mul_pd(t0,t1);
                          return (rcs);
                 }


                     /*
                           Cylinder length much greater then wavelength (h>>gamma).
                           Biscattering RCS, formula 4.3-43
                      */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4343_zmm8r8(const __m512d rcs_inf, // rcs of inifnitely long cylinder (section 4.2)
                                            const __m512d k0,
                                            const __m512d h,
                                            const __m512d psis,
                                            const __m512d psii) {

                           const __m512d _4 = _mm512_set1_pd(4.0f);
                           register __m512d k0h,x0,term1,cpsis,c2psis,rcs;
                           register __m512d term2,spsii,spsis,arg,sarg,rat;
                           k0h   = _mm512_mul_pd(_4,_mm512_mul_pd(k0,h));
                           x0    = _mm512_mul_pd(k0h,k0h);
                           cpsis = xcosf(psis);
                           term1 = _mm512_div_pd(x0,PI);
                           c2psis= _mm512_mul_pd(cpsis,cpsis);
                           term1 = _mm512_mul_pd(term1,_mm512_mul_pd(c2psis,rcs_inf));
                           spsis = xsinf(psis);
                           spsii = xsinf(psii);
                           x0    = _mm512_add_pd(spsis,spsii);
                           arg   = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                           sarg  = xsinf(arg);
                           rat   = _mm512_div_pd(sarg,arg);
                           term2 = _mm512_mul_pd(rat,rat);
                           rcs   = _mm512_mul_pd(term1,term2);
                           return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4343_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) prcs_inf, // rcs of inifnitely long cylinder (section 4.2)
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) ph,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii) {

                           register __m512d prcs_inf  = _mm512_load_pd(&prcs_inf[0]);
                           register __m512d k0        = _mm512_load_pd(&pk0[0]);
                           register __m512d ph        = _mm512_load_pd(&ph[0]);
                           register __m512d psis      = _mm512_load_pd(&ppsis[0]);
                           register __m512d ppsis     = _mm512_load_pd(&ppsii[0]);
                           const __m512d _4 = _mm512_set1_pd(4.0f);
                           register __m512d k0h,x0,term1,cpsis,c2psis,rcs;
                           register __m512d term2,spsii,spsis,arg,sarg,rat;
                           k0h   = _mm512_mul_pd(_4,_mm512_mul_pd(k0,h));
                           x0    = _mm512_mul_pd(k0h,k0h);
                           cpsis = xcosf(psis);
                           term1 = _mm512_div_pd(x0,PI);
                           c2psis= _mm512_mul_pd(cpsis,cpsis);
                           term1 = _mm512_mul_pd(term1,_mm512_mul_pd(c2psis,rcs_inf));
                           spsis = xsinf(psis);
                           spsii = xsinf(psii);
                           x0    = _mm512_add_pd(spsis,spsii);
                           arg   = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                           sarg  = xsinf(arg);
                           rat   = _mm512_div_pd(sarg,arg);
                           term2 = _mm512_mul_pd(rat,rat);
                           rcs   = _mm512_mul_pd(term1,term2);
                           return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4343_zmm8r8_u(const double * __restrict  prcs_inf, // rcs of inifnitely long cylinder (section 4.2)
                                              const double * __restrict  pk0,
                                              const double * __restrict  ph,
                                              const double * __restrict  ppsis,
                                              const double * __restrict  ppsii) {

                           register __m512d prcs_inf  = _mm512_loadu_pd(&prcs_inf[0]);
                           register __m512d k0        = _mm512_loadu_pd(&pk0[0]);
                           register __m512d ph        = _mm512_loadu_pd(&ph[0]);
                           register __m512d psis      = _mm512_loadu_pd(&ppsis[0]);
                           register __m512d ppsis     = _mm512_loadu_pd(&ppsii[0]);
                           const __m512d _4 = _mm512_set1_pd(4.0f);
                           register __m512d k0h,x0,term1,cpsis,c2psis,rcs;
                           register __m512d term2,spsii,spsis,arg,sarg,rat;
                           k0h   = _mm512_mul_pd(_4,_mm512_mul_pd(k0,h));
                           x0    = _mm512_mul_pd(k0h,k0h);
                           cpsis = xcosf(psis);
                           term1 = _mm512_div_pd(x0,PI);
                           c2psis= _mm512_mul_pd(cpsis,cpsis);
                           term1 = _mm512_mul_pd(term1,_mm512_mul_pd(c2psis,rcs_inf));
                           spsis = xsinf(psis);
                           spsii = xsinf(psii);
                           x0    = _mm512_add_pd(spsis,spsii);
                           arg   = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                           sarg  = xsinf(arg);
                           rat   = _mm512_div_pd(sarg,arg);
                           term2 = _mm512_mul_pd(rat,rat);
                           rcs   = _mm512_mul_pd(term1,term2);
                           return (rcs);
                }


                  /*
                         General bistatic scattering RCS from long thin wire.
                         Formula 4.3-44
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4344_zmm8r8(const __m512d h,
                                            const __m512d k0,
                                            const __m512d k0a,
                                            const __m512d psii,
                                            const __m512d psis,
                                            const __m512d gams,
                                            const __m512d gami) {

                          const __m512d c0 = _mm512_set1_pd(12.566370614359172953850573533118);
                          const __m512d c1 = _mm512_set1_pd(2.467401100272339654708622749969);
                          const __m512d c2 = _mm512_set1_pd(0.8905f);
                          register __m512d term1,term2,term3,cgami,cgams,c2gami,c2gams;
                          register __m512d rcs,inv,arg,sarg,rat1,rat2,x0,x1,arg2,larg;
                          register __m512d cpsii,cpsis,fac,c2psii,c2psis,spsii,spsis;
                          fac    = _mm512_mul_pd(c0,_mm512_mul_pd(h,h));
                          arg2   = _mm512_mul_pd(k0a,c2);
                          cpsii  = xcosf(psii);
                          cpsis  = xcosf(psis);
                          c2psii = _mm512_mul_pd(cpsii,cpsii);
                          c2psis = _mm512_mul_pd(cpsis,cpsis);
                          arg2   = _mm512_mul_pd(cpsii,arg2);
                          rat1   = _mm512_div_pd(c2psis,c2psii);
                          cgami  = xcosf(gami);
                          c2gami = _mm512_mul_pd(cgami,cgami);
                          cgams  = xcosf(gams);
                          c2gams = _mm512_mul_pd(cgams,cgams);
                          x0     = _mm512_mul_pd(c2gams,c2gami);
                          term1  = _mm512_mul_pd(fac,_mm512_mul_pd(rat1,x0));
                          larg   = xlogf(arg2);
                          spsii  = xsinf(psii);
                          x1     = _mm512_fmadd_pd(larg,larg,c1);
                          inv    = _mm512_rcp14_pd(x1);
                          spsis  = xsinf(psis);
                          x0     = _mm512_add_pd(spsii,spsis);
                          arg    = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                          sarg   = xsinf(arg);
                          rat2   = _mm512_div_pd(sarg,arg);
                          term2  = _mm512_mul_pd(rat2,rat2);
                          rcs    = _mm512_mul_pd(term1,_mm512_mul_pd(inv,term2));
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4344_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  ph,
                                              const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const double * __restrict __ATTR_ALIGN__(64)  pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64)  ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64)  ppsis,
                                              const double * __restrict __ATTR_ALIGN__(64)  pgams,
                                              const double * __restrict __ATTR_ALIGN__(64)  pgami) {

                          register __m512d h    = _mm512_load_pd(&ph[0]);
                          register __m512d k0   = _mm512_load_pd(&pk0[0]); 
                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          register __m512d psii = _mm512_load_pd(&ppsii[0]);
                          register __m512d psis = _mm512_load_pd(&ppsis[0]);
                          register __m512d gams = _mm512_load_pd(&pgams[0]);
                          register __m512d gami = _mm512_load_pd(&pgami[0]);
                          const __m512d c0 = _mm512_set1_pd(12.566370614359172953850573533118);
                          const __m512d c1 = _mm512_set1_pd(2.467401100272339654708622749969);
                          const __m512d c2 = _mm512_set1_pd(0.8905);
                          register __m512d term1,term2,term3,cgami,cgams,c2gami,c2gams;
                          register __m512d rcs,inv,arg,sarg,rat1,rat2,x0,x1,arg2,larg;
                          register __m512d cpsii,cpsis,fac,c2psii,c2psis,spsii,spsis;
                          fac    = _mm512_mul_pd(c0,_mm512_mul_pd(h,h));
                          arg2   = _mm512_mul_pd(k0a,c2);
                          cpsii  = xcosf(psii);
                          cpsis  = xcosf(psis);
                          c2psii = _mm512_mul_pd(cpsii,cpsii);
                          c2psis = _mm512_mul_pd(cpsis,cpsis);
                          arg2   = _mm512_mul_pd(cpsii,arg2);
                          rat1   = _mm512_div_pd(c2psis,c2psii);
                          cgami  = xcosf(gami);
                          c2gami = _mm512_mul_pd(cgami,cgami);
                          cgams  = xcosf(gams);
                          c2gams = _mm512_mul_pd(cgams,cgams);
                          x0     = _mm512_mul_pd(c2gams,c2gami);
                          term1  = _mm512_mul_pd(fac,_mm512_mul_pd(rat1,x0));
                          larg   = xlogf(arg2);
                          spsii  = xsinf(psii);
                          x1     = _mm512_fmadd_pd(larg,larg,c1);
                          inv    = _mm512_rcp14_pd(x1);
                          spsis  = xsinf(psis);
                          x0     = _mm512_add_pd(spsii,spsis);
                          arg    = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                          sarg   = xsinf(arg);
                          rat2   = _mm512_div_pd(sarg,arg);
                          term2  = _mm512_mul_pd(rat2,rat2);
                          rcs    = _mm512_mul_pd(term1,_mm512_mul_pd(inv,term2));
                          return (rcs);
                }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4344_zmm8r8_u(const double * __restrict   ph,
                                              const double * __restrict   pk0,
                                              const double * __restrict   pk0a,
                                              const double * __restrict   ppsii,
                                              const double * __restrict   ppsis,
                                              const double * __restrict   pgams,
                                              const double * __restrict   pgami) {

                          register __m512d h    = _mm512_loadu_pd(&ph[0]);
                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]); 
                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d psis = _mm512_loadu_pd(&ppsis[0]);
                          register __m512d gams = _mm512_loadu_pd(&pgams[0]);
                          register __m512d gami = _mm512_loadu_pd(&pgami[0]);
                          const __m512d c0 = _mm512_set1_pd(12.566370614359172953850573533118);
                          const __m512d c1 = _mm512_set1_pd(2.467401100272339654708622749969);
                          const __m512d c2 = _mm512_set1_pd(0.8905);
                          register __m512d term1,term2,term3,cgami,cgams,c2gami,c2gams;
                          register __m512d rcs,inv,arg,sarg,rat1,rat2,x0,x1,arg2,larg;
                          register __m512d cpsii,cpsis,fac,c2psii,c2psis,spsii,spsis;
                          fac    = _mm512_mul_pd(c0,_mm512_mul_pd(h,h));
                          arg2   = _mm512_mul_pd(k0a,c2);
                          cpsii  = xcosf(psii);
                          cpsis  = xcosf(psis);
                          c2psii = _mm512_mul_pd(cpsii,cpsii);
                          c2psis = _mm512_mul_pd(cpsis,cpsis);
                          arg2   = _mm512_mul_pd(cpsii,arg2);
                          rat1   = _mm512_div_pd(c2psis,c2psii);
                          cgami  = xcosf(gami);
                          c2gami = _mm512_mul_pd(cgami,cgami);
                          cgams  = xcosf(gams);
                          c2gams = _mm512_mul_pd(cgams,cgams);
                          x0     = _mm512_mul_pd(c2gams,c2gami);
                          term1  = _mm512_mul_pd(fac,_mm512_mul_pd(rat1,x0));
                          larg   = xlogf(arg2);
                          spsii  = xsinf(psii);
                          x1     = _mm512_fmadd_pd(larg,larg,c1);
                          inv    = _mm512_rcp14_pd(x1);
                          spsis  = xsinf(psis);
                          x0     = _mm512_add_pd(spsii,spsis);
                          arg    = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                          sarg   = xsinf(arg);
                          rat2   = _mm512_div_pd(sarg,arg);
                          term2  = _mm512_mul_pd(rat2,rat2);
                          rcs    = _mm512_mul_pd(term1,_mm512_mul_pd(inv,term2));
                          return (rcs);
                }


                    /*

                          General backscatter (only) scattering RCS from long thin wire.
                          Formula 4.3-45
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4345_zmm8r8(const __m512d psi,
                                            const __m512d k0a,
                                            const __m512d gami,
                                            const __m512d gams,
                                            const __m512d k0,
                                            const __m512d h) {

                         const __m512d pi24 = _mm512_set1_pd(2.467401100272339654708622749969);
                         const __m512d _2pi = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d c0   = _mm512_set1_pd(0.8905);
                         register __m512d rat1,arg,sarg,arg2,larg2,k0h,t0,rat;
                         register __m512d rcs,cpsi,cgami,cgams,c2gami,c2gams,spsi;
                         register __m512d x0,x1;
                         k0h   = _mm512_mul_pd(k0,h);
                         t0    = _mm512_mul_pd(_2pi,_mm512_mul_pd(h,h));
                         x0    = _mm512_add_pd(k0h,k0h);
                         spsi  = xsinf(psi);
                         arg   = _mm512_mul_pd(x0,spsi);
                         cpsi  = xcosf(psi);
                         arg2  = _mm512_mul_pd(cpsi,_mm512_mul_pd(k0a,c0));
                         larg  = _mm512_fmadd_pd(arg2,arg2,pi24);
                         sarg  = xsinf(arg);
                         cgams = xcosf(gams);
                         rat   = _mm512_div_pd(sarg,arg);
                         cgami = xcosf(gami);
                         x1    = _mm512_mul_pd(rat,rat);
                         c2gams= _mm512_mul_pd(cgams,cgams);
                         c2gami= _mm512_mul_pd(cgami,cgami);
                         x0    = _mm512_mul_pd(t0,_mm512_mul_pd(c2gams,c2gami));
                         rat1  = _mm512_div_pd(x0,larg);
                         rcs   = _mm512_mul_pd(rat1,x1);
                         return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4345_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  ph,
                                              const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const double * __restrict __ATTR_ALIGN__(64)  pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64)  pgams,
                                              const double * __restrict __ATTR_ALIGN__(64)  pgami,
                                              const double * __restrict __ATTR_ALIGN__(64)  ppsi) {

                         register __m512d h    = _mm512_load_pd(&ph[0]);
                         register __m512d k0   = _mm512_load_pd(&pk0[0]); 
                         register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                         register __m512d psi  = _mm512_load_pd(&ppsi[0]);
                         register __m512d gams = _mm512_load_pd(&pgams[0]);
                         register __m512d gami = _mm512_load_pd(&pgami[0]);
                         const __m512d pi24 = _mm512_set1_pd(2.467401100272339654708622749969);
                         const __m512d _2pi = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d c0   = _mm512_set1_pd(0.8905);
                         register __m512d rat1,arg,sarg,arg2,larg2,k0h,t0,rat;
                         register __m512d rcs,cpsi,cgami,cgams,c2gami,c2gams,spsi;
                         register __m512d x0,x1;
                         k0h   = _mm512_mul_pd(k0,h);
                         t0    = _mm512_mul_pd(_2pi,_mm512_mul_pd(h,h));
                         x0    = _mm512_add_pd(k0h,k0h);
                         spsi  = xsinf(psi);
                         arg   = _mm512_mul_pd(x0,spsi);
                         cpsi  = xcosf(psi);
                         arg2  = _mm512_mul_pd(cpsi,_mm512_mul_pd(k0a,c0));
                         larg  = _mm512_fmadd_pd(arg2,arg2,pi24);
                         sarg  = xsinf(arg);
                         cgams = xcosf(gams);
                         rat   = _mm512_div_pd(sarg,arg);
                         cgami = xcosf(gami);
                         x1    = _mm512_mul_pd(rat,rat);
                         c2gams= _mm512_mul_pd(cgams,cgams);
                         c2gami= _mm512_mul_pd(cgami,cgami);
                         x0    = _mm512_mul_pd(t0,_mm512_mul_pd(c2gams,c2gami));
                         rat1  = _mm512_div_pd(x0,larg);
                         rcs   = _mm512_mul_pd(rat1,x1);
                         return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4345_zmm8r8_u(const double * __restrict   ph,
                                              const double * __restrict   pk0,
                                              const double * __restrict   pk0a,
                                              const double * __restrict   pgams,
                                              const double * __restrict   pgami,
                                              const double * __restrict   ppsi) {

                         register __m512d h    = _mm512_loadu_pd(&ph[0]);
                         register __m512d k0   = _mm512_loadu_pd(&pk0[0]); 
                         register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                         register __m512d psi  = _mm512_loadu_pd(&ppsi[0]);
                         register __m512d gams = _mm512_loadu_pd(&pgams[0]);
                         register __m512d gami = _mm512_loadu_pd(&pgami[0]);
                         const __m512d pi24 = _mm512_set1_pd(2.467401100272339654708622749969);
                         const __m512d _2pi = _mm512_set1_pd(6.283185307179586476925286766559);
                         const __m512d c0   = _mm512_set1_pd(0.8905);
                         register __m512d rat1,arg,sarg,arg2,larg2,k0h,t0,rat;
                         register __m512d rcs,cpsi,cgami,cgams,c2gami,c2gams,spsi;
                         register __m512d x0,x1;
                         k0h   = _mm512_mul_pd(k0,h);
                         t0    = _mm512_mul_pd(_2pi,_mm512_mul_pd(h,h));
                         x0    = _mm512_add_pd(k0h,k0h);
                         spsi  = xsinf(psi);
                         arg   = _mm512_mul_pd(x0,spsi);
                         cpsi  = xcosf(psi);
                         arg2  = _mm512_mul_pd(cpsi,_mm512_mul_pd(k0a,c0));
                         larg  = _mm512_fmadd_pd(arg2,arg2,pi24);
                         sarg  = xsinf(arg);
                         cgams = xcosf(gams);
                         rat   = _mm512_div_pd(sarg,arg);
                         cgami = xcosf(gami);
                         x1    = _mm512_mul_pd(rat,rat);
                         c2gams= _mm512_mul_pd(cgams,cgams);
                         c2gami= _mm512_mul_pd(cgami,cgami);
                         x0    = _mm512_mul_pd(t0,_mm512_mul_pd(c2gams,c2gami));
                         rat1  = _mm512_div_pd(x0,larg);
                         rcs   = _mm512_mul_pd(rat1,x1);
                         return (rcs);
                }


                  /*
                        Backscattering From a Perfectly Conducting Cylinder With Flat Ends.
                        Helper functions, M1,M2 for the main formula 4.3-48
                        Formula 4.3-50

                   */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d M1_f4350_zmm8r8(const __m512d psi) {

                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,M1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          M1   = _mm512_mul_pd(c2,_mm512_add_pd(x0,inv2));
                          return (M1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d M1_f4350_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,M1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          M1   = _mm512_mul_pd(c2,_mm512_add_pd(x0,inv2));
                          return (M1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d M1_f4350_zmm8r8_u(const double * __restrict  ppsi) {

                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,M1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          M1   = _mm512_mul_pd(c2,_mm512_add_pd(x0,inv2));
                          return (M1);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d M2_f4350_zmm8r8(const __m512d psi) {

                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,M2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          M2   = _mm512_mul_pd(c2,_mm512_sub_pd(x0,inv2));
                          return (M2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d M2_f4350_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,M2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          M2   = _mm512_mul_pd(c2,_mm512_sub_pd(x0,inv2));
                          return (M2);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d M2_f4350_zmm8r8_u(const double * __restrict  ppsi) {

                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,M2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          M2   = _mm512_mul_pd(c2,_mm512_sub_pd(x0,inv2));
                          return (M2);
                 }


                    /*
                        Backscattering From a Perfectly Conducting Cylinder With Flat Ends.
                        Helper functions, M1,M2 for the main formula 4.3-48
                        Formula 4.3-51

                   */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d N1_f4351_zmm8r8(const __m512d psi) {

                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n4 = _mm512_set1_pd(-4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,N1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          N1   = _mm512_mul_pd(c2,_mm512_sub_pd(n4,_mm512_sub_pd(x0,inv2)));
                          return (N1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d N1_f4351_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n4 = _mm512_set1_pd(-4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,N1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          N1   = _mm512_mul_pd(c2,_mm512_sub_pd(n4,_mm512_sub_pd(x0,inv2)));
                          return (N1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d N1_f4351_zmm8r8_u(const double * __restrict  ppsi) {

                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n4 = _mm512_set1_pd(-4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,N1,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          N1   = _mm512_mul_pd(c2,_mm512_sub_pd(n4,_mm512_sub_pd(x0,inv2)));
                          return (N1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d N2_f4351_zmm8r8(const __m512d psi) {

                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n4 = _mm512_set1_pd(-4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,N2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          N2   = _mm512_mul_pd(c2,_mm512_add_pd(n4,_mm512_add_pd(x0,inv2)));
                          return (N2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d N2_f4351_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n4 = _mm512_set1_pd(-4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,N2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          N2   = _mm512_mul_pd(c2,_mm512_add_pd(n4,_mm512_add_pd(x0,inv2)));
                          return (N2);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d N2_f4351_zmm8r8_u(const double * __restrict  ppsi) {

                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d n4 = _mm512_set1_pd(-4.0f);
                          const __m512d n1 = _mm512_set1_pd(-1.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          const __m512d c3 = _mm512_set1_pd(0.666666666666666666666666666667);
                          const __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d inv1,inv2,N2,arg1,arg2,carg1,carg2,x0,x1;
                          arg1 = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg1= xcosf(arg1);
                          x0   = _mm512_fmadd_pd(_2,psi,PI);
                          carg1= _mm512_add_pd(c1,carg1);
                          arg2 = _mm512_mul_pd(c3,x0);
                          inv1 = _mm512_rcp14_pd(carg1);
                          carg2= xcosf(arg2);
                          x1   = _mm512_add_pd(c1,carg2);
                          inv2 = _mm512_rcp14_pd(x1);
                          x0   = _mm512_mul_pd(n1,inv1);
                          N2   = _mm512_mul_pd(c2,_mm512_add_pd(n4,_mm512_add_pd(x0,inv2)));
                          return (N2);
                }


                   /*
                        Backscattering From a Perfectly Conducting Cylinder With Flat Ends.
                        Helper functions, M1,M2 for the main formula 4.3-48
                        Formula 4.3-52

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G_f4352_zmm8r8(const __m512d psi) {

                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d _2 = _mm512_set1_pd(-2.0f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          register __m512d G,inv,arg,carg,x0;
                          arg = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_pd(c1,carg);
                          inv = _mm512_rcp14_pd(x0);
                          G   = _mm512_mul_pd(c2,_mm512_sub_pd(n2,inv));
                          return (G);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G_f4352_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d _2 = _mm512_set1_pd(-2.0f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          register __m512d G,inv,arg,carg,x0;
                          arg = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_pd(c1,carg);
                          inv = _mm512_rcp14_pd(x0);
                          G   = _mm512_mul_pd(c2,_mm512_sub_pd(n2,inv));
                          return (G);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d G_f4352_zmm8r8_u(const double * __restrict  ppsi) {

                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d _2 = _mm512_set1_pd(-2.0f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          register __m512d G,inv,arg,carg,x0;
                          arg = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_pd(c1,carg);
                          inv = _mm512_rcp14_pd(x0);
                          G   = _mm512_mul_pd(c2,_mm512_sub_pd(n2,inv));
                          return (G);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F_f4352_zmm8r8(const __m512d psi) {

                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d _2 = _mm512_set1_pd(-2.0f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          register __m512d F,inv,arg,carg,x0;
                          arg = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_pd(c1,carg);
                          inv = _mm512_rcp14_pd(x0);
                          G   = _mm512_mul_pd(c2,_mm512_add_pd(n2,inv));
                          return (F);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F_f4352_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ppsi) {

                          register __m512d psi = _mm512_load_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d _2 = _mm512_set1_pd(-2.0f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          register __m512d F,inv,arg,carg,x0;
                          arg = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_pd(c1,carg);
                          inv = _mm512_rcp14_pd(x0);
                          G   = _mm512_mul_pd(c2,_mm512_add_pd(n2,inv));
                          return (F);
                  }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F_f4352_zmm8r8_u(const double * __restrict ppsi) {

                          register __m512d psi = _mm512_loadu_pd(&ppsi[0]);
                          const __m512d c0 = _mm512_set1_pd(0.333333333333333333333333333333333333333333f);
                          const __m512d _2 = _mm512_set1_pd(-2.0f);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          const __m512d c2 = _mm512_set1_pd(0.577350269189625764509148780502);
                          register __m512d F,inv,arg,carg,x0;
                          arg = _mm512_mul_pd(_mm512_mul_pd(_4,psi),c0);
                          carg= xcosf(arg);
                          x0  = _mm512_add_pd(c1,carg);
                          inv = _mm512_rcp14_pd(x0);
                          G   = _mm512_mul_pd(c2,_mm512_add_pd(n2,inv));
                          return (F);
                  }


                    /*
                           Scattering From Cylinder Near the Specular Direction.
                           Formula 4.3-53
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4353_zmm8r8(const __m512d k0a,
                                            const __m512d k0,
                                            const __m512d h,
                                            const __m512d phi,
                                            const __m512d psii,
                                            const __m512d psis) {

                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,trm2,trm3;
                          register __m512d cphi,cpsis,c2psis,cpsii,c2psii;
                          register __m512d spsii,spsis,arg,sarg,x0,x1;
                          x0    = _mm512_mul_pd(h,h);
                          x1    = _mm512_mul_pd(c1,phi);
                          trm1  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,x0));
                          cpsii = xcosf(psi);
                          cphi  = xcosf(x1);
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          x0    = _mm512_add_pd(spsii,spsis);
                          c2psis= _mm512_mul_pd(cpsis,cpsis);
                          arg   = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                          x1    = _mm512_mul_pd(c2psis,cphi);
                          sarg  = xsinf(arg);
                          trm2  = _mm512_div_pd(x1,cpsii);
                          trm3  = _mm512_div_pd(sarg,arg);
                          x1    = _mm512_mul_pd(trm1,trm2);
                          x0    = _mm512_mul_pd(trm3,trm3)
                          rcs   = _mm512_mul_pd(x1,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4353_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) ph,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsis) {

                          register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          register __m512d k0   = _mm512_load_pd(&pk0[0]);
                          register __m512d h    = _mm512_load_pd(&ph[0]);
                          register __m512d phi  = _mm512_load_pd(&pphi[0]);
                          register __m512d psii = _mm512_load_pd(&ppsii[0]);
                          register __m512d psis = _mm512_load_pd(&ppsis[0]);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,trm2,trm3;
                          register __m512d cphi,cpsis,c2psis,cpsii,c2psii;
                          register __m512d spsii,spsis,arg,sarg,x0,x1;
                          x0    = _mm512_mul_pd(h,h);
                          x1    = _mm512_mul_pd(c1,phi);
                          trm1  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,x0));
                          cpsii = xcosf(psi);
                          cphi  = xcosf(x1);
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          x0    = _mm512_add_pd(spsii,spsis);
                          c2psis= _mm512_mul_pd(cpsis,cpsis);
                          arg   = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                          x1    = _mm512_mul_pd(c2psis,cphi);
                          sarg  = xsinf(arg);
                          trm2  = _mm512_div_pd(x1,cpsii);
                          trm3  = _mm512_div_pd(sarg,arg);
                          x1    = _mm512_mul_pd(trm1,trm2);
                          x0    = _mm512_mul_pd(trm3,trm3)
                          rcs   = _mm512_mul_pd(x1,x0);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4353_zmm8r8_u(const double * __restrict  pk0a,
                                              const double * __restrict  pk0,
                                              const double * __restrict  ph,
                                              const double * __restrict  pphi,
                                              const double * __restrict  ppsii,
                                              const double * __restrict  ppsis) {

                          register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                          register __m512d h    = _mm512_loadu_pd(&ph[0]);
                          register __m512d phi  = _mm512_loadu_pd(&pphi[0]);
                          register __m512d psii = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d psis = _mm512_loadu_pd(&ppsis[0]);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,trm2,trm3;
                          register __m512d cphi,cpsis,c2psis,cpsii,c2psii;
                          register __m512d spsii,spsis,arg,sarg,x0,x1;
                          x0    = _mm512_mul_pd(h,h);
                          x1    = _mm512_mul_pd(c1,phi);
                          trm1  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,x0));
                          cpsii = xcosf(psi);
                          cphi  = xcosf(x1);
                          spsii = xsinf(psii);
                          spsis = xsinf(psis);
                          x0    = _mm512_add_pd(spsii,spsis);
                          c2psis= _mm512_mul_pd(cpsis,cpsis);
                          arg   = _mm512_mul_pd(k0,_mm512_mul_pd(x0,h));
                          x1    = _mm512_mul_pd(c2psis,cphi);
                          sarg  = xsinf(arg);
                          trm2  = _mm512_div_pd(x1,cpsii);
                          trm3  = _mm512_div_pd(sarg,arg);
                          x1    = _mm512_mul_pd(trm1,trm2);
                          x0    = _mm512_mul_pd(trm3,trm3)
                          rcs   = _mm512_mul_pd(x1,x0);
                          return (rcs);
                 }


                   /*

                            Specular direction -- RCS.
                            Formula 4.3-54
                       */

                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4354_zmm8r8(const __m512d k0a,
                                            const __m512d h,
                                            const __m512d psii,
                                            const __m512d phi) {

                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,phi2;
                          register __m512d h2,cpsii,cphi,x0;
                          h2    = _mm512_mul_pd(h,h);
                          phi2  = _mm512_mul_pd(c1,phi);
                          trm1  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          cpsii = xcosf(psii);
                          x0    = _mm512_mul_pd(trm1,cpsii);
                          cphi  = xcosf(phi2);
                          rcs   = _mm512_mul_pd(x0,cphi);
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4354_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ph,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d  k0a  = _mm512_load_pd(&pk0a[0]);
                          register __m512d  h    = _mm512_load_pd(&ph[0]);
                          register __m512d  psii = _mm512_load_pd(&ppsii[0]);
                          register __m512d  phi  = _mm512_load_pd(&pphi[0]);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,phi2;
                          register __m512d h2,cpsii,cphi,x0;
                          h2    = _mm512_mul_pd(h,h);
                          phi2  = _mm512_mul_pd(c1,phi);
                          trm1  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          cpsii = xcosf(psii);
                          x0    = _mm512_mul_pd(trm1,cpsii);
                          cphi  = xcosf(phi2);
                          rcs   = _mm512_mul_pd(x0,cphi);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4354_zmm8r8_u(const double * __restrict  pk0a,
                                              const double * __restrict  ph,
                                              const double * __restrict  ppsii,
                                              const double * __restrict  pphi) {

                          register __m512d  k0a  = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d  h    = _mm512_loadu_pd(&ph[0]);
                          register __m512d  psii = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  phi  = _mm512_loadu_pd(&pphi[0]);
                          const __m512d c1 = _mm512_set1_pd(0.5f);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,phi2;
                          register __m512d h2,cpsii,cphi,x0;
                          h2    = _mm512_mul_pd(h,h);
                          phi2  = _mm512_mul_pd(c1,phi);
                          trm1  = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          cpsii = xcosf(psii);
                          x0    = _mm512_mul_pd(trm1,cpsii);
                          cphi  = xcosf(phi2);
                          rcs   = _mm512_mul_pd(x0,cphi);
                          return (rcs);
                 }


                  /*

                         Backscattering direction -- RCS for incidence angles
                         near broadside.
                         Formula 4.3-54
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4354_zmm8r8(const __m512d k0a,
                                            const __m512d h,
                                            const __m512d k0,
                                            const __m512d psii) {

                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,trm2,cpsii,spsii;
                          register __m512d x0,x1,k0h,h2,arg,sarg;
                          k0h  = _mm512_mul_pd(k0,h);
                          h2   = _mm512_mul_pd(h,h);
                          x0   = _mm512_add_pd(k0h,k0h); 
                          x1   = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          cpsii= xcosf(psi);
                          spsii= xsinf(psi);
                          trm1 = _mm512_mul_pd(x1,cpsii);
                          arg  = _mm512_mul_pd(x0,spsii);
                          sarg = xsinf(arg);
                          x0   = _mm512_div_pd(sarg,arg);
                          trm2 = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_mul_pd(trm1,trm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4354_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) ph,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) ppsii) {

                          register __m512d  k0a  = _mm512_load_pd(&pk0a[0]);
                          register __m512d  h    = _mm512_load_pd(&ph[0]);
                          register __m512d  psii = _mm512_load_pd(&ppsii[0]);
                          register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,trm2,cpsii,spsii;
                          register __m512d x0,x1,k0h,h2,arg,sarg;
                          k0h  = _mm512_mul_pd(k0,h);
                          h2   = _mm512_mul_pd(h,h);
                          x0   = _mm512_add_pd(k0h,k0h); 
                          x1   = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          cpsii= xcosf(psi);
                          spsii= xsinf(psi);
                          trm1 = _mm512_mul_pd(x1,cpsii);
                          arg  = _mm512_mul_pd(x0,spsii);
                          sarg = xsinf(arg);
                          x0   = _mm512_div_pd(sarg,arg);
                          trm2 = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_mul_pd(trm1,trm2);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4354_zmm8r8_u(const double * __restrict  pk0a,
                                              const double * __restrict  ph,
                                              const double * __restrict  pk0,
                                              const double * __restrict  ppsii) {

                          register __m512d  k0a  = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d  h    = _mm512_loadu_pd(&ph[0]);
                          register __m512d  psii = _mm512_loadu_pd(&ppsii[0]);
                          register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,trm1,trm2,cpsii,spsii;
                          register __m512d x0,x1,k0h,h2,arg,sarg;
                          k0h  = _mm512_mul_pd(k0,h);
                          h2   = _mm512_mul_pd(h,h);
                          x0   = _mm512_add_pd(k0h,k0h); 
                          x1   = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          cpsii= xcosf(psi);
                          spsii= xsinf(psi);
                          trm1 = _mm512_mul_pd(x1,cpsii);
                          arg  = _mm512_mul_pd(x0,spsii);
                          sarg = xsinf(arg);
                          x0   = _mm512_div_pd(sarg,arg);
                          trm2 = _mm512_mul_pd(x0,x0);
                          rcs  = _mm512_mul_pd(trm1,trm2);
                          return (rcs);
                }


                 /*

                        Broadside (psi == 0) RCS.
                        Formula 4.3-56
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4356_zmm8r8(const __m512d k0a,
                                            const __m512d h) {

                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,h2;
                          h2 = _mm512_mul_pd(h,h);
                          rcs = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          return (rcs); 
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4356_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) ph) {

                          register __m512d  k0a  = _mm512_load_pd(&pk0a[0]);
                          register __m512d  h    = _mm512_load_pd(&ph[0]);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,h2;
                          h2 = _mm512_mul_pd(h,h);
                          rcs = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          return (rcs); 
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4356_zmm8r8_u(const double * __restrict  pk0a,
                                            const double * __restrict  ph) {

                          register __m512d  k0a  = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d  h    = _mm512_loadu_pd(&ph[0]);
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,h2;
                          h2 = _mm512_mul_pd(h,h);
                          rcs = _mm512_mul_pd(_4,_mm512_mul_pd(k0a,h2));
                          return (rcs); 
                }


                  /*
                       Elliptical cylinders.
                   */


                   /*
                         Low-frequency approximations (k0a<0.5, k0b<0.5)
                         TM-case,formula 4.4-11
                    */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4411_zmm8r8(const __m512d a,
                                         const __m512d b,
                                         const __m512d k0,
                                         __m512d * __restrict TMr,
                                         __m512d * __restrict TMi) {

                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         const __m512d imn  = _mm512_set1_pd(-1.57079632679489661923132169164);
                         const __m512d imp  = _mm512_set1_pd(1.57079632679489661923132169164);
                         const __m512d c0   = _mm512_set1_pd(0.8905);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d ab2,c0k0,arg,larg;
                         register __m512d invr,invi;
                         ab2  = _mm512_mul_pd(_mm512_add_pd(a,b),hlf);
                         c0k0 = _mm512_mul_pd(c0,k0);
                         arg  = _mm512_mul_pd(ab2,c0k0);
                         larg = xlogf(arg);
                         cdiv_zmm8r8(_1,_1,larg,imn,&invr,&invi);
                         *TMr = _mm512_mul_pd(imp,invr);
                         *TMi = _mm512_mul_pd(imp,invi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4411_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                         const double * __restrict __ATTR_ALIGN__(64) pb,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         double * __restrict __ATTR_ALIGN__(64) TMr,
                                         double * __restrict __ATTR_ALIGN__(64) TMi) {

                         register __m512d a = _mm512_load_pd(&pa[0]);
                         register __m512d b = _mm512_load_pd(&pb[0]);
                         register __m512d k0= _mm512_load_pd(&pk0[0]);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         const __m512d imn  = _mm512_set1_pd(-1.57079632679489661923132169164);
                         const __m512d imp  = _mm512_set1_pd(1.57079632679489661923132169164);
                         const __m512d c0   = _mm512_set1_pd(0.8905);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d ab2,c0k0,arg,larg;
                         register __m512d invr,invi;
                         ab2  = _mm512_mul_pd(_mm512_add_pd(a,b),hlf);
                         c0k0 = _mm512_mul_pd(c0,k0);
                         arg  = _mm512_mul_pd(ab2,c0k0);
                         larg = xlogf(arg);
                         cdiv_zmm8r8(_1,_1,larg,imn,&invr,&invi);
                         _mm512_store_pd(&TMr[0], _mm512_mul_pd(imp,invr));
                         _mm512_store_pd(&TMi[0], _mm512_mul_pd(imp,invi));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4411_zmm8r8_u(const double * __restrict  pa,
                                           const double * __restrict  pb,
                                           const double * __restrict  pk0,
                                           double * __restrict  TMr,
                                           double * __restrict  TMi) {

                         register __m512d a = _mm512_loadu_pd(&pa[0]);
                         register __m512d b = _mm512_loadu_pd(&pb[0]);
                         register __m512d k0= _mm512_loadu_pd(&pk0[0]);
                         const __m512d hlf  = _mm512_set1_pd(0.5f);
                         const __m512d imn  = _mm512_set1_pd(-1.57079632679489661923132169164);
                         const __m512d imp  = _mm512_set1_pd(1.57079632679489661923132169164);
                         const __m512d c0   = _mm512_set1_pd(0.8905);
                         const __m512d _1   = _mm512_set1_pd(1.0f);
                         register __m512d ab2,c0k0,arg,larg;
                         register __m512d invr,invi;
                         ab2  = _mm512_mul_pd(_mm512_add_pd(a,b),hlf);
                         c0k0 = _mm512_mul_pd(c0,k0);
                         arg  = _mm512_mul_pd(ab2,c0k0);
                         larg = xlogf(arg);
                         cdiv_zmm8r8(_1,_1,larg,imn,&invr,&invi);
                         _mm512_storeu_pd(&TMr[0], _mm512_mul_pd(imp,invr));
                         _mm512_storeu_pd(&TMi[0], _mm512_mul_pd(imp,invi));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4412_zmm8r8(const __m512d k0a,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d phi1,
                                         const __m512d phi2,
                                         __m512d * __restrict TEr,
                                         __m512d * __restrict TEi) {

                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d k0a2,ba,cphi1,sphi1,trm1,trm2,_1ba,x0,x1;
                        register __m512d cphi2,sphi2;
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        ba    = _mm512_div_pd(b,a);
                        cphi1 = xcosf(phi1);
                        _1ba  = _mm512_add_pd(_1,ba);
                        sphi1 = xsinf(phi1);
                        x0    = _mm512_mul_pd(pi4,k0a2);
                        cphi2 = xcosf(phi2);
                        x1    = _mm512_add_pd(ba,_1ba);
                        sphi2 = xsinf(phi2);
                        trm1  = _mm512_mul_pd(x0,x1);
                        x0    = _mm512_fmadd_pd(cphi2,cphi1,_mm512_mul_pd(sphi2,sphi1));
                        trm2  = _mm512_mul_pd(ba,x0);
                        x1    = _mm512_mul_pd(trm1,trm2);
                        *TEr  = nIi;
                        *TEi  = _mm512_mul_pd(nIi,x1);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4412_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pb,
                                           const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                           const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                           double * __restrict __ATTR_ALIGN__(64) TEr,
                                           double * __restrict __ATTR_ALIGN__(64) TEi) {

                        register __m512d a    = _mm512_load_pd(&pa[0]);
                        register __m512d b    = _mm512_load_pd(&pb[0]);
                        register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d phi1 = _mm512_load_pd(&phi1[0]);
                        register __m512d phi2 = _mm512_load_pd(&phi2[0]);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d k0a2,ba,cphi1,sphi1,trm1,trm2,_1ba,x0,x1;
                        register __m512d cphi2,sphi2;
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        ba    = _mm512_div_pd(b,a);
                        cphi1 = xcosf(phi1);
                        _1ba  = _mm512_add_pd(_1,ba);
                        sphi1 = xsinf(phi1);
                        x0    = _mm512_mul_pd(pi4,k0a2);
                        cphi2 = xcosf(phi2);
                        x1    = _mm512_add_pd(ba,_1ba);
                        sphi2 = xsinf(phi2);
                        trm1  = _mm512_mul_pd(x0,x1);
                        x0    = _mm512_fmadd_pd(cphi2,cphi1,_mm512_mul_pd(sphi2,sphi1));
                        trm2  = _mm512_mul_pd(ba,x0);
                        x1    = _mm512_mul_pd(trm1,trm2);
                        _mm512_store_pd(&TEr[0] ,nIi);
                        _mm512_store_pd(&TEi[0] ,_mm512_mul_pd(nIi,x1));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4412_zmm8r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  pa,
                                           const double * __restrict  pb,
                                           const double * __restrict  pphi1,
                                           const double * __restrict  pphi2,
                                           __m512d * __restrict TEr,
                                           __m512d * __restrict TEi) {

                        register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d phi1 = _mm512_loadu_pd(&phi1[0]);
                        register __m512d phi2 = _mm512_loadu_pd(&phi2[0]);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        register __m512d k0a2,ba,cphi1,sphi1,trm1,trm2,_1ba,x0,x1;
                        register __m512d cphi2,sphi2;
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        ba    = _mm512_div_pd(b,a);
                        cphi1 = xcosf(phi1);
                        _1ba  = _mm512_add_pd(_1,ba);
                        sphi1 = xsinf(phi1);
                        x0    = _mm512_mul_pd(pi4,k0a2);
                        cphi2 = xcosf(phi2);
                        x1    = _mm512_add_pd(ba,_1ba);
                        sphi2 = xsinf(phi2);
                        trm1  = _mm512_mul_pd(x0,x1);
                        x0    = _mm512_fmadd_pd(cphi2,cphi1,_mm512_mul_pd(sphi2,sphi1));
                        trm2  = _mm512_mul_pd(ba,x0);
                        x1    = _mm512_mul_pd(trm1,trm2);
                        _mm512_store_pd(&TEr[0], nIi);
                        _mm512_store_pd(&TEi[0], _mm512_mul_pd(nIi,x1));
                }


                 /*
                       TM-case, RCS.
                       Formula 4.4-13
                  */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4413_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d k0) {

                          const __m512d pi2 = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(0.8905);
                          const __m512d pi24= _mm512_set1_pd(2.467401100272339654708622749969);
                          register __m512d rcs,abh,k0abh,num,sqr1,sqr2,c0k0,arg,larg,x0,den,x1;
                          abh = _mm512_mul_pd(_mm512_add_pd(a,b),hlf);
                          c0k0= _mm512_mul_pd(c0,k0);
                          num = _mm512_mul_pd(pi2,abh);
                          arg = _mm512_mul_pd(c0k0,abh);
                          larg= xlogf(arg);
                          x0  = _mm512_fmadd_pd(larg,larg,pi24);
                          sqr1= _mm512_sqrt_pd(_mm512_mul_pd(k0,abh));
                          sqr2= _mm512_sqrt_pd(x0);
                          den = _mm512_mul_pd(sqr1,sqr2);
                          x1  = _mm512_mul_pd(den,den);
                          rcs = _mm512_div_pd(num,x1);
                          return (rcs);
                }
                                            


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4413_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pa,
                                              const double * __restrict __ATTR_ALIGN__(64)  pb,
                                              const double * __restrict __ATTR_ALIGN__(64)  pk0) {

                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d b    = _mm512_load_pd(&pb[0]);
                          register __m512d k0   = _mm512_load_pd(&pk0[0]);
                          const __m512d pi2 = _mm512_set1_pd(9.869604401089358618834490999876f);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(0.8905f);
                          const __m512d pi24= _mm512_set1_pd(2.467401100272339654708622749969f);
                          register __m512d rcs,abh,k0abh,num,sqr1,sqr2,c0k0,arg,larg,x0,den,x1;
                          abh = _mm512_mul_pd(_mm512_add_pd(a,b),hlf);
                          c0k0= _mm512_mul_pd(c0,k0);
                          num = _mm512_mul_pd(pi2,abh);
                          arg = _mm512_mul_pd(c0k0,abh);
                          larg= xlogf(arg);
                          x0  = _mm512_fmadd_pd(larg,larg,pi24);
                          sqr1= _mm512_sqrt_pd(_mm512_mul_pd(k0,abh));
                          sqr2= _mm512_sqrt_pd(x0);
                          den = _mm512_mul_pd(sqr1,sqr2);
                          x1  = _mm512_mul_pd(den,den);
                          rcs = _mm512_div_pd(num,x1);
                          return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4413_zmm8r8_u(const double * __restrict   pa,
                                              const double * __restrict   pb,
                                              const double * __restrict   pk0) {

                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d b    = _mm512_loadu_pd(&pb[0]);
                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                          const __m512d pi2 = _mm512_set1_pd(9.869604401089358618834490999876);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(0.8905);
                          const __m512d pi24= _mm512_set1_pd(2.467401100272339654708622749969);
                          register __m512d rcs,abh,k0abh,num,sqr1,sqr2,c0k0,arg,larg,x0,den,x1;
                          abh = _mm512_mul_pd(_mm512_add_pd(a,b),hlf);
                          c0k0= _mm512_mul_pd(c0,k0);
                          num = _mm512_mul_pd(pi2,abh);
                          arg = _mm512_mul_pd(c0k0,abh);
                          larg= xlogf(arg);
                          x0  = _mm512_fmadd_pd(larg,larg,pi24);
                          sqr1= _mm512_sqrt_pd(_mm512_mul_pd(k0,abh));
                          sqr2= _mm512_sqrt_pd(x0);
                          den = _mm512_mul_pd(sqr1,sqr2);
                          x1  = _mm512_mul_pd(den,den);
                          rcs = _mm512_div_pd(num,x1);
                          return (rcs);
                }


                    /*
                         High frequency approximations (k0a>5, k0b>5)
                         TM-case, formula 4.4-15
                      */

#include "GMS_simd_utils.hpp"


                    /*
                        Helper function for testing the condition of high-frequency limit.
                        Page. 322.

                     */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __mmask16 
                   TM_f4415_helper_zmm8r8(const __m512d k0,
                                           const __m512d a,
                                           const __m512d phi1,
                                           const __m512d phi2,
                                           const __m512d b) {
                      const __m512d c0 = _mm512_set1_pd(0.166666666666666666666666666667f);
                      register __m512d a2,b2,sphi1,cphi1,trm1,trm2,root6;
                      register __m512d k02,absp,sphi1s,cphi1s,k0a2,k0b2,x0;
                      __mmask16 m;
                      k02  = _mm512_mul_pd(k0,k0);
                      a2   = _mm512_mul_pd(a,a);
                      k0a2 = _mm512_mul_pd(k02,a2);
                      b2   = _mm512_mul_pd(b,b);
                      k0b2 = _mm512_mul_pd(k02,b2);
                      cphi1= xcosf(phi1);
                      absp = _mm512_abs_pd(_mm512_sub_pd(phi2,phi1));
                      cphi1s = _mm512_mul_pd(cphi1,cphi1)
                      sphi1= xsinf(phi1);
                      trm1 = _mm512_sub_pd(PI,absp);
                      sphi1s = _mm512_mul_pd(sphi1,sphi1);
                      trm2 = _mm512_fmadd_pd(k02a2,sphi1s,_mm512_mul_pd(k02b2,cphi1s));
                      x0   = _mm512_pow_pd(trm2,c0);
                      root6= _mm512_rcp14_pd(x0);
                      m    = _mm512_cmp_mask_pd(trm1,root6,_CMP_GT_OQ);
                      return (m);
                }

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4415_zmm8r8(const __m512d phi1,
                                         const __m512d phi2,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d k0,
                                         __m512d * __restrict TMr,
                                         __m512d * __restrict TMi,
                                         bool & status) {

                        using namespace gms::math;
                        __mmask16 m = TM_f4415_helper_zmm8r8(k0,a,phi1,phi2,b);
                        if(!m) {
                           status = false;
                           return;
                        }
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d ip4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d arg1,arg2,carg1,carg2,sarg2,sqr1,ear,eai,cer,cei,trm1;
                        register __m512d f,rho,a2b2,a2,b2,b2a2,k0a,cphi2,cphi1,sphi1,sphi2,frat;
                        register __m512d cphis,sphis,rhod,rhorat,x0,x1,tmp1,tmp2,b2a2s,carg2s,sarg2s;
                        arg1  = _mm512_mul_pd(_mm512_sub_pd(phi2,phi1),hlf);
                        a2    = _mm512_mul_pd(a,a);
                        b2    = _mm512_mul_pd(b,b);
                        k0a   = _mm512_mul_pd(k0,a);
                        arg2  = _mm512_mul_pd(_mm512_add_pd(phi2,phi1),hlf);
                        carg1 = xcosf(arg1);
                        a2b2  = _mm512_mul_pd(a2,b2);
                        b2a2  = _mm512_div_pd(b2,a2);
                        cphi1 = xcosf(phi1);
                        sphi1 = xsinf(phi1);
                        trm1  = _mm512_sqrt_pd(_mm512_mul_pd(PI,carg1));
                        cphi2 = xcosf(phi2);
                        sphi2 = xsinf(phi2);
                        cphis = _mm512_add_pd(cphi1,cphi2);
                        carg2 = xcosf(arg2);
                        sphis = _mm512_add_pd(sphi1,sphi2);
                        sarg2 = xsinf(arg2);
                        x0    = _mm512_mul_pd(carg2,carg2);
                        x1    = _mm512_mul_pd(sarg2,sarg2);
                        rhod  = _mm512_fmadd_pd(a2,x0,_mm512_mul_pd(b2,x1));
                        b2a2s = _mm512_mul_pd(b2a2,sphis);
                        tmp1  = _mm512_pow_pd(rhod,_mm512_set1_pd(1.5f));
                        rhorat= _mm512_div_pd(a2b2,tmp1);
                        x0    = _mm512_fmadd_pd(sarg2,b2a2s,carg2);
                        carg2s= _mm512_mul_pd(carg2,carg2)
                        tmp2  = _mm512_mul_pd(cphis,x0);
                        sarg2s= _mm512_mul_pd(sarg2,sarg2);
                        x1    = _mm512_fmadd_pd(b2a2,sarg2s,carg2s);
                        tmp1  = _mm512_sqrt_pd(x1);
                        frat  = _mm512_div_pd(tmp2,tmp1);
                        trm1  = zmm8r8_negate(trm1);
                        ear   = _mm512_add_pd(nIr,ip4);
                        x0    = _mm512_mul_pd(_mm512_sqrt_pd(_mm512_mul_pd(k0,rhorat)),hlf);
                        eai   = _mm512_mul_pd(nIi,_mm512_mul_pd(k0a,frat));
                        eai   = _mm512_add_pd(eai,ip4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        x1    = _mm512_mul_pd(trm1,x0);
                        *TMr = _mm512_mul_pd(x1,cer);
                        *TMi = _mm512_mul_pd(x1,cei);
                        status = true;
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4415_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         const double * __restrict __ATTR_ALIGN__(64) pb,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         double * __restrict __ATTR_ALIGN__(64) TMr,
                                         double * __restrict __ATTR_ALIGN__(64) TMi,
                                         bool & status) {

                        using namespace gms::math;

                        register __m512d phi1 = _mm512_load_pd(&phi1[0]);
                        register __m512d phi2 = _mm512_load_pd(&phi2[0]);
                        register __m512d a    = _mm512_load_pd(&pa[0]);
                        register __m512d b    = _mm512_load_pd(&pb[0]);
                        register __m512d k0   = _mm512_load_pd(&pk0[0]);
                        __mmask16 m = TM_f4415_helper_zmm8r8(k0,a,phi1,phi2,b);
                        if(!m) {
                           status = false;
                           return;
                        }
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d ip4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d arg1,arg2,carg1,carg2,sarg2,sqr1,ear,eai,cer,cei,trm1;
                        register __m512d f,rho,a2b2,a2,b2,b2a2,k0a,cphi2,cphi1,sphi1,sphi2,frat;
                        register __m512d cphis,sphis,rhod,rhorat,x0,x1,tmp1,tmp2,b2a2s,carg2s,sarg2s;
                        arg1  = _mm512_mul_pd(_mm512_sub_pd(phi2,phi1),hlf);
                        a2    = _mm512_mul_pd(a,a);
                        b2    = _mm512_mul_pd(b,b);
                        k0a   = _mm512_mul_pd(k0,a);
                        arg2  = _mm512_mul_pd(_mm512_add_pd(phi2,phi1),hlf);
                        carg1 = xcosf(arg1);
                        a2b2  = _mm512_mul_pd(a2,b2);
                        b2a2  = _mm512_div_pd(b2,a2);
                        cphi1 = xcosf(phi1);
                        sphi1 = xsinf(phi1);
                        trm1  = _mm512_sqrt_pd(_mm512_mul_pd(PI,carg1));
                        cphi2 = xcosf(phi2);
                        sphi2 = xsinf(phi2);
                        cphis = _mm512_add_pd(cphi1,cphi2);
                        carg2 = xcosf(arg2);
                        sphis = _mm512_add_pd(sphi1,sphi2);
                        sarg2 = xsinf(arg2);
                        x0    = _mm512_mul_pd(carg2,carg2);
                        x1    = _mm512_mul_pd(sarg2,sarg2);
                        rhod  = _mm512_fmadd_pd(a2,x0,_mm512_mul_pd(b2,x1));
                        b2a2s = _mm512_mul_pd(b2a2,sphis);
                        tmp1  = _mm512_pow_pd(rhod,_mm512_set1_pd(1.5f));
                        rhorat= _mm512_div_pd(a2b2,tmp1);
                        x0    = _mm512_fmadd_pd(sarg2,b2a2s,carg2);
                        carg2s= _mm512_mul_pd(carg2,carg2)
                        tmp2  = _mm512_mul_pd(cphis,x0);
                        sarg2s= _mm512_mul_pd(sarg2,sarg2);
                        x1    = _mm512_fmadd_pd(b2a2,sarg2s,carg2s);
                        tmp1  = _mm512_sqrt_pd(x1);
                        frat  = _mm512_div_pd(tmp2,tmp1);
                        trm1  = zmm8r8_negate(trm1);
                        ear   = _mm512_add_pd(nIr,ip4);
                        x0    = _mm512_mul_pd(_mm512_sqrt_pd(_mm512_mul_pd(k0,rhorat)),hlf);
                        eai   = _mm512_mul_pd(nIi,_mm512_mul_pd(k0a,frat));
                        eai   = _mm512_add_pd(eai,ip4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        x1    = _mm512_mul_pd(trm1,x0);
                        _mm512_store_pd(&TMr[0] ,_mm512_mul_pd(x1,cer));
                        _mm512_store_pd(&TMi[0] ,_mm512_mul_pd(x1,cei));
                        status = true;
                 }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4415_zmm8r8_u(const double * __restrict  pphi1,
                                         const double * __restrict  pphi2,
                                         const double * __restrict  pa,
                                         const double * __restrict  pb,
                                         const double * __restrict  pk0,
                                         double * __restrict  TMr,
                                         double * __restrict  TMi,
                                         bool & status) {

                        using namespace gms::math;
                        register __m512d phi1 = _mm512_loadu_pd(&phi1[0]);
                        register __m512d phi2 = _mm512_loadu_pd(&phi2[0]);
                        register __m512d a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                        __mmask16 m = TM_f4415_helper_zmm8r8(k0,a,phi1,phi2,b);
                        if(!m) {
                           status = false;
                           return;
                        }
                        const __m512d hlf = _mm512_set1_pd(0.5f);
                        const __m512d ip4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d arg1,arg2,carg1,carg2,sarg2,sqr1,ear,eai,cer,cei,trm1;
                        register __m512d f,rho,a2b2,a2,b2,b2a2,k0a,cphi2,cphi1,sphi1,sphi2,frat;
                        register __m512d cphis,sphis,rhod,rhorat,x0,x1,tmp1,tmp2,b2a2s,carg2s,sarg2s;
                        arg1  = _mm512_mul_pd(_mm512_sub_pd(phi2,phi1),hlf);
                        a2    = _mm512_mul_pd(a,a);
                        b2    = _mm512_mul_pd(b,b);
                        k0a   = _mm512_mul_pd(k0,a);
                        arg2  = _mm512_mul_pd(_mm512_add_pd(phi2,phi1),hlf);
                        carg1 = xcosf(arg1);
                        a2b2  = _mm512_mul_pd(a2,b2);
                        b2a2  = _mm512_div_pd(b2,a2);
                        cphi1 = xcosf(phi1);
                        sphi1 = xsinf(phi1);
                        trm1  = _mm512_sqrt_pd(_mm512_mul_pd(PI,carg1));
                        cphi2 = xcosf(phi2);
                        sphi2 = xsinf(phi2);
                        cphis = _mm512_add_pd(cphi1,cphi2);
                        carg2 = xcosf(arg2);
                        sphis = _mm512_add_pd(sphi1,sphi2);
                        sarg2 = xsinf(arg2);
                        x0    = _mm512_mul_pd(carg2,carg2);
                        x1    = _mm512_mul_pd(sarg2,sarg2);
                        rhod  = _mm512_fmadd_pd(a2,x0,_mm512_mul_pd(b2,x1));
                        b2a2s = _mm512_mul_pd(b2a2,sphis);
                        tmp1  = _mm512_pow_pd(rhod,_mm512_set1_pd(1.5f));
                        rhorat= _mm512_div_pd(a2b2,tmp1);
                        x0    = _mm512_fmadd_pd(sarg2,b2a2s,carg2);
                        carg2s= _mm512_mul_pd(carg2,carg2)
                        tmp2  = _mm512_mul_pd(cphis,x0);
                        sarg2s= _mm512_mul_pd(sarg2,sarg2);
                        x1    = _mm512_fmadd_pd(b2a2,sarg2s,carg2s);
                        tmp1  = _mm512_sqrt_pd(x1);
                        frat  = _mm512_div_pd(tmp2,tmp1);
                        trm1  = zmm8r8_negate(trm1);
                        ear   = _mm512_add_pd(nIr,ip4);
                        x0    = _mm512_mul_pd(_mm512_sqrt_pd(_mm512_mul_pd(k0,rhorat)),hlf);
                        eai   = _mm512_mul_pd(nIi,_mm512_mul_pd(k0a,frat));
                        eai   = _mm512_add_pd(eai,ip4);
                        cexp_zmm8r8(ear,eai,&cer,&cei);
                        x1    = _mm512_mul_pd(trm1,x0);
                        _mm512_storeu_pd(&TMr[0] ,_mm512_mul_pd(x1,cer));
                        _mm512_storeu_pd(&TMi[0] ,_mm512_mul_pd(x1,cei));
                        status = true;
                 }


                    /*
                         High frequency approximations (k0a>5, k0b>5)
                         TE-case, formula 4.4-16
                      */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4416_zmm8r8(const __m512d phi1,
                                         const __m512d phi2,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d k0,
                                         __m512d * __restrict TEr,
                                         __m512d * __restrict TEi,
                                         bool & status) {

                        __m512d resr,resi;
                        TM_f4415_zmm8r8(phi1,phi2,a,b,k0,&resr,&resi,status);
                        if(!status) {
                           return;
                        }
                         else {
                           *TEr = _mm512_mul_pd(nIi,resr);
                           *TEi = _mm512_mul_pd(nIi,resi);
                        }
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4416_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         const double * __restrict __ATTR_ALIGN__(64) pb,
                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         double * __restrict __ATTR_ALIGN__(64) TEr,
                                         double * __restrict __ATTR_ALIGN__(64) TEi,
                                         bool & status) {

                        
                        TM_f4415_zmm8r8_a(phi1,phi2,a,b,k0,TEr,TEi,status);
                        if(!status) {
                           return;
                        }
                         else {
                           _mm512_store_pd(&TEr[0] ,_mm512_mul_pd(nIi,_mm512_load_pd(&TEr[0])));
                           _mm512_store_pd(&TEi[0] ,_mm512_mul_pd(nIi,_mm512_load_pd(&TEi[0])));
                        }
               }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4416_zmm8r8_u(const double * __restrict  pphi1,
                                         const double * __restrict  pphi2,
                                         const double * __restrict  pa,
                                         const double * __restrict  pb,
                                         const double * __restrict  pk0,
                                         double * __restrict  TEr,
                                         double * __restrict  TEi,
                                         bool & status) {

                       
                        TM_f4415_zmm8r8_u(phi1,phi2,a,b,k0,TEr,TEi,status);
                        if(!status) {
                           return;
                        }
                         else {
                           _mm512_storeu_pd(&TEr[0] ,_mm512_mul_pd(nIi,_mm512_loadu_pd(&TEr[0])));
                           _mm512_storeu_pd(&TEi[0] ,_mm512_mul_pd(nIi,_mm512_loadu_pd(&TEr[0])));
                        }
               }


                 /*

                        Bistatic scattering width.
                        Formula 4.4-19
                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4419_zmm8r8(const __m512d phi1,
                                            const __m512d phi2,
                                            const __m512d a,
                                            const __m512d b) {

                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(1.5f);
                          register __m512d rcs,a2,b2,a2b2,num;
                          register __m512d arg,carg,carg2,sarg,sarg2;
                          register __m512d pow32,x0;
                          a2   = _mm512_mul_pd(a,a);
                          arg  = _mm512_mul_pd(_mm512_add_pd(phi2,phi1),hlf);
                          b2   = _mm512_mul_pd(b,b);
                          carg = xcosf(arg);
                          num  = _mm512_mul_pd(PI,_mm512_mul_pd(a2,b2));
                          sarg = xsinf(arg);
                          carg2= _mm512_mul_pd(carg,carg);
                          sarg2= _mm512_mul_pd(sarg,sarg);
                          x0   = _mm512_fmadd_pd(a2,carg2,_mm512_mul_pd(b2,sarg2));
                          pow32= _mm512_pow_pd(x0,c0);
                          rcs  = _mm512_div_pd(num,pow32);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4419_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb) {

                          register __m512d phi1 = _mm512_load_pd(&pphi1[0]);
                          register __m512d phi2 = _mm512_load_pd(&pphi2[0]);
                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d b    = _mm512_load_pd(&pb[0]);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(1.5f);
                          register __m512d rcs,a2,b2,a2b2,num;
                          register __m512d arg,carg,carg2,sarg,sarg2;
                          register __m512d pow32,x0;
                          a2   = _mm512_mul_pd(a,a);
                          arg  = _mm512_mul_pd(_mm512_add_pd(phi2,phi1),hlf);
                          b2   = _mm512_mul_pd(b,b);
                          carg = xcosf(arg);
                          num  = _mm512_mul_pd(PI,_mm512_mul_pd(a2,b2));
                          sarg = xsinf(arg);
                          carg2= _mm512_mul_pd(carg,carg);
                          sarg2= _mm512_mul_pd(sarg,sarg);
                          x0   = _mm512_fmadd_pd(a2,carg2,_mm512_mul_pd(b2,sarg2));
                          pow32= _mm512_pow_pd(x0,c0);
                          rcs  = _mm512_div_pd(num,pow32);
                          return (rcs);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4419_zmm8r8_u(const double * __restrict  pphi1,
                                              const double * __restrict  pphi2,
                                              const double * __restrict  pa,
                                              const double * __restrict  pb) {

                          register __m512d phi1 = _mm512_loadu_pd(&pphi1[0]);
                          register __m512d phi2 = _mm512_loadu_pd(&pphi2[0]);
                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d b    = _mm512_loadu_pd(&pb[0]);
                          const __m512d hlf = _mm512_set1_pd(0.5f);
                          const __m512d c0  = _mm512_set1_pd(1.5f);
                          register __m512d rcs,a2,b2,a2b2,num;
                          register __m512d arg,carg,carg2,sarg,sarg2;
                          register __m512d pow32,x0;
                          a2   = _mm512_mul_pd(a,a);
                          arg  = _mm512_mul_pd(_mm512_add_pd(phi2,phi1),hlf);
                          b2   = _mm512_mul_pd(b,b);
                          carg = xcosf(arg);
                          num  = _mm512_mul_pd(PI,_mm512_mul_pd(a2,b2));
                          sarg = xsinf(arg);
                          carg2= _mm512_mul_pd(carg,carg);
                          sarg2= _mm512_mul_pd(sarg,sarg);
                          x0   = _mm512_fmadd_pd(a2,carg2,_mm512_mul_pd(b2,sarg2));
                          pow32= _mm512_pow_pd(x0,c0);
                          rcs  = _mm512_div_pd(num,pow32);
                          return (rcs);
                 }


                   /*

                          Backscattering width, for phi2 == phi1.
                          Formula 4.4-20
                      */


                     
                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4420_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d phi) {

                          const __m512d c0  = _mm512_set1_pd(1.5f);
                          register __m512d rcs,a2,b2,a2b2,num;
                          register __m512d carg,carg2,sarg,sarg2;
                          register __m512d pow32,x0;
                          a2   = _mm512_mul_pd(a,a);
                          carg = xcosf(phi);
                          b2   = _mm512_mul_pd(b,b);
                          carg2= _mm512_mul_pd(carg,carg);
                          num  = _mm512_mul_pd(PI,_mm512_mul_pd(a2,b2));
                          sarg = xsinf(phi);
                          sarg2= _mm512_mul_pd(sarg,sarg);
                          x0   = _mm512_fmadd_pd(a2,carg2,_mm512_mul_pd(b2,sarg2));
                          pow32= _mm512_pow_pd(x0,c0);
                          rcs  = _mm512_div_pd(num,pow32);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4420_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi) {

                          register __m512d phi2 = _mm512_load_pd(&pphi[0]);
                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d b    = _mm512_load_pd(&pb[0]);
                          const __m512d c0  = _mm512_set1_pd(1.5f);
                          register __m512d rcs,a2,b2,a2b2,num;
                          register __m512d carg,carg2,sarg,sarg2;
                          register __m512d pow32,x0;
                          a2   = _mm512_mul_pd(a,a);
                          carg = xcosf(phi);
                          b2   = _mm512_mul_pd(b,b);
                          carg2= _mm512_mul_pd(carg,carg);
                          num  = _mm512_mul_pd(PI,_mm512_mul_pd(a2,b2));
                          sarg = xsinf(phi);
                          sarg2= _mm512_mul_pd(sarg,sarg);
                          x0   = _mm512_fmadd_pd(a2,carg2,_mm512_mul_pd(b2,sarg2));
                          pow32= _mm512_pow_pd(x0,c0);
                          rcs  = _mm512_div_pd(num,pow32);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4420_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pphi) {

                          register __m512d phi2 = _mm512_loadu_pd(&pphi[0]);
                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d b    = _mm512_loadu_pd(&pb[0]);
                          const __m512d c0  = _mm512_set1_pd(1.5f);
                          register __m512d rcs,a2,b2,a2b2,num;
                          register __m512d carg,carg2,sarg,sarg2;
                          register __m512d pow32,x0;
                          a2   = _mm512_mul_pd(a,a);
                          carg = xcosf(phi);
                          b2   = _mm512_mul_pd(b,b);
                          carg2= _mm512_mul_pd(carg,carg);
                          num  = _mm512_mul_pd(PI,_mm512_mul_pd(a2,b2));
                          sarg = xsinf(phi);
                          sarg2= _mm512_mul_pd(sarg,sarg);
                          x0   = _mm512_fmadd_pd(a2,carg2,_mm512_mul_pd(b2,sarg2));
                          pow32= _mm512_pow_pd(x0,c0);
                          rcs  = _mm512_div_pd(num,pow32);
                          return (rcs);
                 }


                   /*
                        Forward scattering pattern and width.
                        Formula 4.4-23 a scattering amplitude

                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __mmask16 
                   T_f4423_helper_zmm8r8( const __m512d k0,
                                           const __m512d a,
                                           const __m512d phi1,
                                           const __m512d phi2,
                                           const __m512d b) {
                      const __m512d c0 = _mm512_set1_pd(0.166666666666666666666666666667f);
                      register __m512d a2,b2,sphi1,cphi1,trm1,trm2,root6;
                      register __m512d k02,absp,sphi1s,cphi1s,k0a2,k0b2,x0;
                      __mmask16 m;
                      k02  = _mm512_mul_pd(k0,k0);
                      a2   = _mm512_mul_pd(a,a);
                      k0a2 = _mm512_mul_pd(k02,a2);
                      b2   = _mm512_mul_pd(b,b);
                      k0b2 = _mm512_mul_pd(k02,b2);
                      cphi1= xcosf(phi1);
                      trm1 = _mm512_sub_pd(phi1,phi2);
                      cphi1s = _mm512_add_pd(PI,_mm512_mul_pd(cphi1,cphi1));
                      sphi1= xsinf(phi1);
                      sphi1s = _mm512_mul_pd(sphi1,sphi1);
                      trm2 = _mm512_fmadd_pd(k02a2,sphi1s,_mm512_mul_pd(k02b2,cphi1s));
                      x0   = _mm512_pow_pd(trm2,c0);
                      root6= _mm512_rcp14_pd(x0);
                      m    = _mm512_cmp_mask_pd(_mm512_abs_pd(trm1),root6,_CMP_LT_OQ);
                      return (m);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d T_f4423_zmm8r8(const __m512d a,
                                          const __m512d b,
                                          const __m512d phi1,
                                          const __m512d phi2,
                                          const __m512d k0,
                                          bool & status) {

                          using namespace gms::math;
                          __mmask16 m = TM_f4423_helper_zmm8r8(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          register __m512d T,k0c,c,alp,a2,b2;
                          register __m512d sphi,sphi2,cphi,cphi2;
                          register __m512d arg,sarg,rat,x0;
                          a2   = _mm512_mul_pd(a,a);
                          alp  = _mm512_add_pd(PI,_mm512_sub_pd(phi2,phi1));
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          x0   = _mm512_fmadd_pd(a2,cphi2,_mm512_mul_pd(b2,sphi2));
                          c    = _mm512_sqrt_pd(x0);
                          k0c  = _mm512_mul_pd(k0,c);
                          arg  = _mm512_mul_pd(k0c,alp);
                          sarg = xsinf(arg);
                          k0c  = negate_zmm8r8(k0c);
                          rat  = _mm512_div_pd(sarg,arg);
                          T    = _mm512_mul_pd(k0c,rat);
                          status = true;
                          return (T);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d T_f4423_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                          const double * __restrict __ATTR_ALIGN__(64) pb,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                          const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
                                          bool & status) {

                          using namespace gms::math;
                          register __m512d phi1 = _mm512_load_pd(&phi1[0]);
                          register __m512d phi2 = _mm512_load_pd(&phi2[0]);
                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d b    = _mm512_load_pd(&pb[0]);
                          register __m512d k0   = _mm512_load_pd(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm8r8(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          register __m512d T,k0c,c,alp,a2,b2;
                          register __m512d sphi,sphi2,cphi,cphi2;
                          register __m512d arg,sarg,rat,x0;
                          a2   = _mm512_mul_pd(a,a);
                          alp  = _mm512_add_pd(PI,_mm512_sub_pd(phi2,phi1));
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          x0   = _mm512_fmadd_pd(a2,cphi2,_mm512_mul_pd(b2,sphi2));
                          c    = _mm512_sqrt_pd(x0);
                          k0c  = _mm512_mul_pd(k0,c);
                          arg  = _mm512_mul_pd(k0c,alp);
                          sarg = xsinf(arg);
                          k0c  = negate_zmm8r8(k0c);
                          rat  = _mm512_div_pd(sarg,arg);
                          T    = _mm512_mul_pd(k0c,rat);
                          status = true;
                          return (T);
                 }


                     __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d T_f4423_zmm8r8_u(const double * __restrict pa,
                                          const double * __restrict  pb,
                                          const double * __restrict  pphi1,
                                          const double * __restrict  pphi2,
                                          const double * __restrict  pk0,
                                          bool & status) {

                          using namespace gms::math;
                          register __m512d phi1 = _mm512_loadu_pd(&phi1[0]);
                          register __m512d phi2 = _mm512_loadu_pd(&phi2[0]);
                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d b    = _mm512_loadu_pd(&pb[0]);
                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm8r8(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          register __m512d T,k0c,c,alp,a2,b2;
                          register __m512d sphi,sphi2,cphi,cphi2;
                          register __m512d arg,sarg,rat,x0;
                          a2   = _mm512_mul_pd(a,a);
                          alp  = _mm512_add_pd(PI,_mm512_sub_pd(phi2,phi1));
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          x0   = _mm512_fmadd_pd(a2,cphi2,_mm512_mul_pd(b2,sphi2));
                          c    = _mm512_sqrt_pd(x0);
                          k0c  = _mm512_mul_pd(k0,c);
                          arg  = _mm512_mul_pd(k0c,alp);
                          sarg = xsinf(arg);
                          k0c  = negate_zmm8r8(k0c);
                          rat  = _mm512_div_pd(sarg,arg);
                          T    = _mm512_mul_pd(k0c,rat);
                          status = true;
                          return (T);
                 }


                   /*
                          Scattering width near the forward direction.
                          Formula 4.4-24

                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4424_zmm8r8(const __m512d a,
                                          const __m512d b,
                                          const __m512d phi1,
                                          const __m512d phi2,
                                          const __m512d k0,
                                          bool & status) {

                         
                          __mmask16 m = TM_f4423_helper_zmm8r8(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,k0c,c,alp,a2,b2;
                          register __m512d sphi,sphi2,cphi,cphi2;
                          register __m512d arg,sarg,rat,x0,x1,x2;
                          a2   = _mm512_mul_pd(a,a);
                          alp  = _mm512_add_pd(PI,_mm512_sub_pd(phi2,phi1));
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          x0   = _mm512_fmadd_pd(a2,cphi2,_mm512_mul_pd(b2,sphi2));
                          c    = _mm512_sqrt_pd(x0);
                          k0c  = _mm512_mul_pd(k0,c);
                          arg  = _mm512_mul_pd(k0c,alp);
                          sarg = xsinf(arg);
                          x1   = _mm512_mul_pd(_4,_mm512_mul_pd(k0c,k0c));
                          rat  = _mm512_div_pd(sarg,arg);
                          x2   = _mm512_mul_pd(rat,rat);
                          rcs  = _mm512_mul_pd(k0c,rat);
                          status = true;
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4424_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              bool & status) {

                          
                          register __m512d phi1 = _mm512_load_pd(&phi1[0]);
                          register __m512d phi2 = _mm512_load_pd(&phi2[0]);
                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d b    = _mm512_load_pd(&pb[0]);
                          register __m512d k0   = _mm512_load_pd(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm8r8(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,k0c,c,alp,a2,b2;
                          register __m512d sphi,sphi2,cphi,cphi2;
                          register __m512d arg,sarg,rat,x0,x1,x2;
                          a2   = _mm512_mul_pd(a,a);
                          alp  = _mm512_add_pd(PI,_mm512_sub_pd(phi2,phi1));
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          x0   = _mm512_fmadd_pd(a2,cphi2,_mm512_mul_pd(b2,sphi2));
                          c    = _mm512_sqrt_pd(x0);
                          k0c  = _mm512_mul_pd(k0,c);
                          arg  = _mm512_mul_pd(k0c,alp);
                          sarg = xsinf(arg);
                          x1   = _mm512_mul_pd(_4,_mm512_mul_pd(k0c,k0c));
                          rat  = _mm512_div_pd(sarg,arg);
                          x2   = _mm512_mul_pd(rat,rat);
                          rcs  = _mm512_mul_pd(k0c,rat);
                          status = true;
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4424_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pphi1,
                                              const double * __restrict  pphi2,
                                              const double * __restrict  pk0,
                                              bool & status) {

                          
                          register __m512d phi1 = _mm512_loadu_pd(&phi1[0]);
                          register __m512d phi2 = _mm512_loadu_pd(&phi2[0]);
                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d b    = _mm512_loadu_pd(&pb[0]);
                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]);
                          __mmask16 m = TM_f4423_helper_zmm8r8(k0,a,phi1,phi2,b);
                          if(!m) {
                             status = false;
                             return;
                          }
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,k0c,c,alp,a2,b2;
                          register __m512d sphi,sphi2,cphi,cphi2;
                          register __m512d arg,sarg,rat,x0,x1,x2;
                          a2   = _mm512_mul_pd(a,a);
                          alp  = _mm512_add_pd(PI,_mm512_sub_pd(phi2,phi1));
                          b2   = _mm512_mul_pd(b,b);
                          sphi = xsinf(phi1);
                          cphi = xcosf(phi1);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          x0   = _mm512_fmadd_pd(a2,cphi2,_mm512_mul_pd(b2,sphi2));
                          c    = _mm512_sqrt_pd(x0);
                          k0c  = _mm512_mul_pd(k0,c);
                          arg  = _mm512_mul_pd(k0c,alp);
                          sarg = xsinf(arg);
                          x1   = _mm512_mul_pd(_4,_mm512_mul_pd(k0c,k0c));
                          rat  = _mm512_div_pd(sarg,arg);
                          x2   = _mm512_mul_pd(rat,rat);
                          rcs  = _mm512_mul_pd(k0c,rat);
                          status = true;
                          return (rcs);
                 }


                   /*
                         Scattering width in the exact forward direction (alpha == 0).
                         Formula 4.4-25
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4425_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                            const __m512d phi) {

                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,_4k0,a2,b2,sphi,sphi2;
                          register __m512d cphi,cphi2;
                          a2   = _mm512_mul_pd(a,a);
                          sphi = xsinf(phi);
                          b2   = _mm512_mul_pd(b,b);
                          cphi = xcosf(phi);
                          _4k0 = _mm512_mul_pd(_4,k0);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          x0   = _mm512_fmadd_pd(a2,sphi2,_mm512_mul_pd(b2,cphi2));
                          rcs  = _mm512_mul_pd(_4k0,x0);
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4425_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0,
                                              const double * __restrict __ATTR_ALIGN__(64)  pa,
                                              const double * __restrict __ATTR_ALIGN__(64)  pb,
                                              const double * __restrict __ATTR_ALIGN__(64)  pphi) {

                          register __m512d phi2 = _mm512_load_pd(&phi2[0]);
                          register __m512d a    = _mm512_load_pd(&pa[0]);
                          register __m512d b    = _mm512_load_pd(&pb[0]);
                          register __m512d k0   = _mm512_load_pd(&pk0[0]); 
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,_4k0,a2,b2,sphi,sphi2;
                          register __m512d cphi,cphi2;
                          a2   = _mm512_mul_pd(a,a);
                          sphi = xsinf(phi);
                          b2   = _mm512_mul_pd(b,b);
                          cphi = xcosf(phi);
                          _4k0 = _mm512_mul_pd(_4,k0);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          x0   = _mm512_fmadd_pd(a2,sphi2,_mm512_mul_pd(b2,cphi2));
                          rcs  = _mm512_mul_pd(_4k0,x0);
                          return (rcs);
                 }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4425_zmm8r8_u(const double * __restrict   pk0,
                                              const double * __restrict   pa,
                                              const double * __restrict   pb,
                                              const double * __restrict  pphi) {

                          register __m512d phi2 = _mm512_loadu_pd(&phi2[0]);
                          register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          register __m512d b    = _mm512_loadu_pd(&pb[0]);
                          register __m512d k0   = _mm512_loadu_pd(&pk0[0]); 
                          const __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d rcs,_4k0,a2,b2,sphi,sphi2;
                          register __m512d cphi,cphi2;
                          a2   = _mm512_mul_pd(a,a);
                          sphi = xsinf(phi);
                          b2   = _mm512_mul_pd(b,b);
                          cphi = xcosf(phi);
                          _4k0 = _mm512_mul_pd(_4,k0);
                          cphi2= _mm512_mul_pd(cphi,cphi);
                          sphi2= _mm512_mul_pd(sphi,sphi);
                          x0   = _mm512_fmadd_pd(a2,sphi2,_mm512_mul_pd(b2,cphi2));
                          rcs  = _mm512_mul_pd(_4k0,x0);
                          return (rcs);
                 }


                    /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          TM-case, formula 4.4-26
                     */

                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4426_zmm8r8(const __m512d k0,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d phi1,
                                         const __m512d phi2,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui,
                                         __m512d * __restrict TMr,
                                         __m512d * __restrict TMi) {

                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,murpba,muipba;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512d facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        faci   = _mm512_mul_pd(pi4,_mm512_mul_pd(k0a2,ba));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        *TMr   = _mm512_mul_pd(facr,tmpr);
                        *TMi   = _mm512_mul_pd(faci,tmpi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4426_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         const double * __restrict __ATTR_ALIGN__(64) pb,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                         const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                         const double * __restrict __ATTR_ALIGN__(64) pmur,
                                         const double * __restrict __ATTR_ALIGN__(64) pmui,
                                         double * __restrict __ATTR_ALIGN__(64) TMr,
                                         double * __restrict __ATTR_ALIGN__(64) TMi) {

                        register __m512d k0    = _mm512_load_pd(&pk0[0]);
                        register __m512d a     = _mm512_load_pd(&pa[0]);
                        register __m512d b     = _mm512_load_pd(&pb[0]);
                        register __m512d phi1  = _mm512_load_pd(&pphi1[0]);
                        register __m512d phi2  = _mm512_load_pd(&pphi2[0]);
                        register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_load_pd(&pmur[0]);
                        register __m512d mui   = _mm512_load_pd(&pmui[0]);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,murpba,muipba;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512d facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        faci   = _mm512_mul_pd(pi4,_mm512_mul_pd(k0a2,ba));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_store_pd(&TMr[0]  ,_mm512_mul_pd(facr,tmpr));
                        _mm512_store_pd(&TMi[0]  ,_mm512_mul_pd(faci,tmpi));
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TM_f4426_zmm8r8_u(const double * __restrict  pk0,
                                         const double * __restrict  pa,
                                         const double * __restrict  pb,
                                         const double * __restrict  pphi1,
                                         const double * __restrict pphi2,
                                         const double * __restrict  pepsr,
                                         const double * __restrict  pepsi,
                                         const double * __restrict  pmur,
                                         const double * __restrict  pmui,
                                         double * __restrict  TMr,
                                         double * __restrict  TMi) {

                        register __m512d k0    = _mm512_loadu_pd(&pk0[0]);
                        register __m512d a     = _mm512_loadu_pd(&pa[0]);
                        register __m512d b     = _mm512_loadu_pd(&pb[0]);
                        register __m512d phi1  = _mm512_loadu_pd(&pphi1[0]);
                        register __m512d phi2  = _mm512_loadu_pd(&pphi2[0]);
                        register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui   = _mm512_loadu_pd(&pmui[0]);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,murpba,muipba;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512d facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        faci   = _mm512_mul_pd(pi4,_mm512_mul_pd(k0a2,ba));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_storeu_pd(&TMr[0]  ,_mm512_mul_pd(facr,tmpr));
                        _mm512_storeu_pd(&TMi[0]  ,_mm512_mul_pd(faci,tmpi));
                }


                   /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          TE-case, formula 4.4-27
                     */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4427_zmm8r8(const __m512d k0,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d phi1,
                                         const __m512d phi2,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui,
                                         __m512d * __restrict TEr,
                                         __m512d * __restrict TEi) {

                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512d facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        faci   = _mm512_mul_pd(pi4,_mm512_mul_pd(k0a2,ba));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        *TEr   = _mm512_mul_pd(facr,tmpr);
                        *TEi   = _mm512_mul_pd(faci,tmpi);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4427_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                         const double * __restrict __ATTR_ALIGN__(64) pa,
                                         const double * __restrict __ATTR_ALIGN__(64) pb,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                         const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                         const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                         const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                         const double * __restrict __ATTR_ALIGN__(64) pmur,
                                         const double * __restrict __ATTR_ALIGN__(64) pmui,
                                         double * __restrict __ATTR_ALIGN__(64) TEr,
                                         double * __restrict __ATTR_ALIGN__(64) TEi) {

                        register __m512d k0    = _mm512_load_pd(&pk0[0]);
                        register __m512d a     = _mm512_load_pd(&pa[0]);
                        register __m512d b     = _mm512_load_pd(&pb[0]);
                        register __m512d phi1  = _mm512_load_pd(&pphi1[0]);
                        register __m512d phi2  = _mm512_load_pd(&pphi2[0]);
                        register __m512d epsr  = _mm512_load_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_load_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_load_pd(&pmur[0]);
                        register __m512d mui   = _mm512_load_pd(&pmui[0]);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512d facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        faci   = _mm512_mul_pd(pi4,_mm512_mul_pd(k0a2,ba));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_store_pd(&TEr[0] ,_mm512_mul_pd(facr,tmpr));
                        _mm512_store_pd(&TEi[0] ,_mm512_mul_pd(faci,tmpi));
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void TE_f4427_zmm8r8_u(const double * __restrict  pk0,
                                         const double * __restrict  pa,
                                         const double * __restrict  pb,
                                         const double * __restrict  pphi1,
                                         const double * __restrict pphi2,
                                         const double * __restrict  pepsr,
                                         const double * __restrict  pepsi,
                                         const double * __restrict pmur,
                                         const double * __restrict  pmui,
                                         double * __restrict  TEr,
                                         double * __restrict  TEi) {

                        register __m512d k0    = _mm512_loadu_pd(&pk0[0]);
                        register __m512d a     = _mm512_loadu_pd(&pa[0]);
                        register __m512d b     = _mm512_loadu_pd(&pb[0]);
                        register __m512d phi1  = _mm512_loadu_pd(&pphi1[0]);
                        register __m512d phi2  = _mm512_loadu_pd(&pphi2[0]);
                        register __m512d epsr  = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d epsi  = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d mur   = _mm512_loadu_pd(&pmur[0]);
                        register __m512d mui   = _mm512_loadu_pd(&pmui[0]);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i;
                        register __m512d facr,faci,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi;
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        facr   = Ir;
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        faci   = _mm512_mul_pd(pi4,_mm512_mul_pd(k0a2,ba));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        _mm512_storeu_pd(&TEr[0] ,_mm512_mul_pd(facr,tmpr));
                        _mm512_storeu_pd(&TEi[0] ,_mm512_mul_pd(faci,tmpi));
                }


                  /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Bistatic scattering width (RCS).
                          TM-case.
                          Formula 4.4-28
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4428_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                            const __m512d phi1,
                                            const __m512d phi2,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {
                                        
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512d cphi2,cphi1,sphi2,sphi1,murpba,muipba,b2,a2,pia;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512d fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_pd(b,b);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        pia    = _mm512_mul_pd(PI,a);
                        a2     = _mm512_mul_pd(a,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_pd(b2,a2)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        k0a3   = _mm512_mul_pd(k0a2,k0a);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4428_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {
                                        
                        register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d  a    = _mm512_load_pd(&pa[0]);
                        register __m512d  b    = _mm512_load_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_load_pd(&pphi1[0]);
                        register __m512d  phi2 = _mm512_load_pd(&pphi2[0]);
                        register __m512d  epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_load_pd(&pmui[0]);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512d cphi2,cphi1,sphi2,sphi1,murpba,muipba,b2,a2,pia;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512d fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_pd(b,b);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        pia    = _mm512_mul_pd(PI,a);
                        a2     = _mm512_mul_pd(a,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_pd(b2,a2)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        k0a3   = _mm512_mul_pd(k0a2,k0a);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4428_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pphi1,
                                              const double * __restrict  pphi2,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {
                                        
                        register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d  a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d  b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_loadu_pd(&pphi1[0]);
                        register __m512d  phi2 = _mm512_loadu_pd(&pphi2[0]);
                        register __m512d  epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_loadu_pd(&pmui[0]);
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512d cphi2,cphi1,sphi2,sphi1,murpba,muipba,b2,a2,pia;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512d fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_pd(b,b);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        pia    = _mm512_mul_pd(PI,a);
                        a2     = _mm512_mul_pd(a,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_pd(b2,a2)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        k0a3   = _mm512_mul_pd(k0a2,k0a);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                  /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Bistatic scattering width (RCS).
                          TE-case.
                          Formula 4.4-29

                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4429_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                         const __m512d phi1,
                                         const __m512d phi2,
                                         const __m512d epsr,
                                         const __m512d epsi,
                                         const __m512d mur,
                                         const __m512d mui) {
                                         
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba,cabs;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512d fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_pd(a,a);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_pd(b,b);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        b2a2   = _mm512_div_pd(b2,a2);
                        pia    = _mm512_mul_pd(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        k0a3   = _mm512_mul_pd(k0a,k0a2);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                    __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4429_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {
                              
                        register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d  a    = _mm512_load_pd(&pa[0]);
                        register __m512d  b    = _mm512_load_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_load_pd(&pphi1[0]);
                        register __m512d  phi2 = _mm512_load_pd(&pphi2[0]);
                        register __m512d  epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_load_pd(&pmui[0]);           
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba,cabs;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512d fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_pd(a,a);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_pd(b,b);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        b2a2   = _mm512_div_pd(b2,a2);
                        pia    = _mm512_mul_pd(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        k0a3   = _mm512_mul_pd(k0a,k0a2);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }



                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4429_zmm8r8_u(const double * __restrict pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict pb,
                                              const double * __restrict pphi1,
                                              const double * __restrict  pphi2,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {
                              
                        register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d  a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d  b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_loadu_pd(&pphi1[0]);
                        register __m512d  phi2 = _mm512_loadu_pd(&pphi2[0]);
                        register __m512d  epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_loadu_pd(&pmui[0]);           
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi2,cphi1,sphi2,sphi1,epsrpba,epsipba,cabs;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512d fac,_1ba,cphit,sphit,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_pd(a,a);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_pd(b,b);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        b2a2   = _mm512_div_pd(b2,a2);
                        pia    = _mm512_mul_pd(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        k0a3   = _mm512_mul_pd(k0a,k0a2);
                        cphi2  = xcosf(phi2);
                        cphit  = _mm512_mul_pd(cphi2,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi2  = xsinf(phi2);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        sphit  = _mm512_mul_pd(sphi2,sphi1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphit,murpba);
                        t1i    = _mm512_div_pd(cphit,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphit,murmba);
                        t2i    = _mm512_div_pd(sphit,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                    /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Backscattering  width (RCS).
                          TM-case.
                          Formula 4.4-30
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4430_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                            const __m512d phi1,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {
                                        
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512d cphi1,sphi1,cphi1s,sphi1s,murpba,muipba,b2,a2,pia;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512d fac,_1ba,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_pd(b,b);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        pia    = _mm512_mul_pd(PI,a);
                        a2     = _mm512_mul_pd(a,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_pd(b2,a2)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi1s = _mm512_mul_pd(cphi1,cphi1);
                        k0a3   = _mm512_mul_pd(k0a2,k0a);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi1s  = _mm512_mul_pd(sphi1,sphi1);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphi1s,murpba);
                        t1i    = _mm512_div_pd(cphi1s,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphi1s,murmba);
                        t2i    = _mm512_div_pd(sphi1s,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4430_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {
                                 
                                 
                        register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d  a    = _mm512_load_pd(&pa[0]);
                        register __m512d  b    = _mm512_load_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_load_pd(&pphi1[0]);
                        register __m512d  epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_load_pd(&pmui[0]);              
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512d cphi1,sphi1,cphi1s,sphi1s,murpba,muipba,b2,a2,pia;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512d fac,_1ba,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_pd(b,b);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        pia    = _mm512_mul_pd(PI,a);
                        a2     = _mm512_mul_pd(a,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_pd(b2,a2)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi1s = _mm512_mul_pd(cphi1,cphi1);
                        k0a3   = _mm512_mul_pd(k0a2,k0a);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi1s  = _mm512_mul_pd(sphi1,sphi1);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphi1s,murpba);
                        t1i    = _mm512_div_pd(cphi1s,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphi1s,murmba);
                        t2i    = _mm512_div_pd(sphi1s,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                  __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4430_zmm8r8_u(const double * __restrict pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict pphi1,
                                              const double * __restrict pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {
                                 
                                 
                        register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d  a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d  b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_loadu_pd(&pphi1[0]);
                        register __m512d  epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_loadu_pd(&pmui[0]);              
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,ba,epsrm1,epsim1,murm1,muim1,k0a3;
                        register __m512d cphi1,sphi1,cphi1s,sphi1s,murpba,muipba,b2,a2,pia;
                        register __m512d murmba,muimba,t0r,t0i,t1r,t1i,t2r,t2i,b2a2;
                        register __m512d fac,_1ba,t3r,t3i,tmpr,tmpi,cabs;
                        b2     = _mm512_mul_pd(b,b);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        ba     = _mm512_div_pd(b,a);
                        pia    = _mm512_mul_pd(PI,a);
                        a2     = _mm512_mul_pd(a,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        sphi1  = xsinf(phi1)
                        b2a2   = _mm512_div_pd(b2,a2)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        cphi1s = _mm512_mul_pd(cphi1,cphi1);
                        k0a3   = _mm512_mul_pd(k0a2,k0a);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(epsrm1,murm1);
                        sphi1s  = _mm512_mul_pd(sphi1,sphi1);
                        t0i    = _mm512_sub_pd(epsim1,muim1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        murpba = _mm512_add_pd(mur,ba);
                        muipba = _mm512_add_pd(mui,ba);
                        t1r    = _mm512_div_pd(cphi1s,murpba);
                        t1i    = _mm512_div_pd(cphi1s,muipba);
                        murmba = _mm512_fmadd_pd(mur,ba,_1);
                        muimba = _mm512_fmadd_pd(mui,ba,_1);
                        t2r    = _mm512_div_pd(sphi1s,murmba);
                        t2i    = _mm512_div_pd(sphi1s,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                  /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Backscattering  width (RCS).
                          TE-case.
                          Formula 4.4-31
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4431_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                            const __m512d phi1,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {
                                         
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi1,sphi1,cphi1s,sphi1s,epsrpba,epsipba,cabs;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512d fac,_1ba,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_pd(a,a);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_pd(b,b);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        b2a2   = _mm512_div_pd(b2,a2);
                        pia    = _mm512_mul_pd(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        k0a3   = _mm512_mul_pd(k0a,k0a2);
                        cphi1s = _mm512_mul_pd(cphi1,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi1s = _mm512_mul_pd(sphi1,sphi1);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphi1s,murpba);
                        t1i    = _mm512_div_pd(cphi1s,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphi1s,murmba);
                        t2i    = _mm512_div_pd(sphi1s,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4431_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {
                                
                        register __m512d  k0   = _mm512_load_pd(&pk0[0]);
                        register __m512d  a    = _mm512_load_pd(&pa[0]);
                        register __m512d  b    = _mm512_load_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_load_pd(&pphi1[0]);
                        register __m512d  epsr = _mm512_load_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_load_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_load_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_load_pd(&pmui[0]);           
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi1,sphi1,cphi1s,sphi1s,epsrpba,epsipba,cabs;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512d fac,_1ba,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_pd(a,a);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_pd(b,b);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        b2a2   = _mm512_div_pd(b2,a2);
                        pia    = _mm512_mul_pd(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        k0a3   = _mm512_mul_pd(k0a,k0a2);
                        cphi1s = _mm512_mul_pd(cphi1,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi1s = _mm512_mul_pd(sphi1,sphi1);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphi1s,murpba);
                        t1i    = _mm512_div_pd(cphi1s,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphi1s,murmba);
                        t2i    = _mm512_div_pd(sphi1s,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4431_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pphi1,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {
                                
                        register __m512d  k0   = _mm512_loadu_pd(&pk0[0]);
                        register __m512d  a    = _mm512_loadu_pd(&pa[0]);
                        register __m512d  b    = _mm512_loadu_pd(&pb[0]);
                        register __m512d  phi1 = _mm512_loadu_pd(&pphi1[0]);
                        register __m512d  epsr = _mm512_loadu_pd(&pepsr[0]);
                        register __m512d  epsi = _mm512_loadu_pd(&pepsi[0]);
                        register __m512d  mur  = _mm512_loadu_pd(&pmur[0]);
                        register __m512d  mui  = _mm512_loadu_pd(&pmui[0]);           
                        const __m512d _1  = _mm512_set1_pd(1.0f);
                        const __m512d pi4 = _mm512_set1_pd(0.78539816339744830961566084582);
                        register __m512d rcs,k0a,k0a2,k0a3,ba,epsrm1,epsim1,murm1,muim1;
                        register __m512d cphi1,sphi1,cphi1s,sphi1s,epsrpba,epsipba,cabs;
                        register __m512d epsrmba,epsimba,t0r,t0i,t1r,t1i,t2r,t2i,b2,a2;
                        register __m512d fac,_1ba,t3r,t3i,tmpr,tmpi,b2a2,pia;
                        a2     = _mm512_mul_pd(a,a);
                        k0a    = _mm512_mul_pd(k0,a);
                        cphi1  = xcosf(phi1);
                        b2     = _mm512_mul_pd(b,b);
                        ba     = _mm512_div_pd(b,a);
                        epsrm1 = _mm512_sub_pd(epsr,_1);
                        b2a2   = _mm512_div_pd(b2,a2);
                        pia    = _mm512_mul_pd(PI,a);
                        sphi1  = xsinf(phi1)
                        k0a2   = _mm512_mul_pd(k0a,k0a);
                        epsim1 = _mm512_sub_pd(epsi,_1);
                        k0a3   = _mm512_mul_pd(k0a,k0a2);
                        cphi1s = _mm512_mul_pd(cphi1,cphi1);
                        murm1  = _mm512_sub_pd(mur,_1);
                        muim1  = _mm512_sub_pd(mui,_1);
                        _1ba   = _mm512_add_pd(_1,ba);
                        t0r    = _mm512_sub_pd(murm1,epsrm1);
                        sphi1s = _mm512_mul_pd(sphi1,sphi1);
                        t0i    = _mm512_sub_pd(muim1,epsim1);
                        fac    = _mm512_mul_pd(_mm512_mul_pd(pia,pi4),
                                               _mm512_mul_pd(k0a3,b2a2));
                        epsrpba= _mm512_add_pd(epsr,ba);
                        epsipba= _mm512_add_pd(epsi,ba);
                        t1r    = _mm512_div_pd(cphi1s,murpba);
                        t1i    = _mm512_div_pd(cphi1s,muipba);
                        epsrmba= _mm512_fmadd_pd(epsr,ba,_1);
                        epsimba= _mm512_fmadd_pd(epsi,ba,_1);
                        t2r    = _mm512_div_pd(sphi1s,murmba);
                        t2i    = _mm512_div_pd(sphi1s,muimba);
                        t3r    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1r,t2r));
                        t3i    = _mm512_mul_pd(_1ba,_mm512_add_pd(t1i,t2i));
                        cmul_zmm8r8(t0r,t0i,t3r,t3i,&tmpr,&tmpi);
                        cabs   = cabs_zmm8r8(tmpr,tmpi);
                        rcs    = _mm512_mul_pd(fac,cabs);
                        return (rcs);
                }


                 /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Forward scattering (phi2 = pi+phi1)  width (RCS).
                          TM-case.
                          Formula 4.4-32
                    */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4432_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                            const __m512d phi1,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {

                           return (rcs_f4430_zmm8r8(k0,a,b,phi1,epsr,epsi,mur,mui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4432_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {

                           return (rcs_f4430_zmm8r8_a(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4432_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pphi1,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {

                           return (rcs_f4430_zmm8r8_u(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               } 


                 /*
                          Infinitely long homogenous cylinder at normal
                          incidence.
                          Low frequency approximation (k0a<0.5,k0b<0.5,k1a<0.5,k1b<0.5)
                          Forward scattering (phi2 = pi+phi1)  width (RCS).
                          TE-case.
                          Formula 4.4-33

                   */


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4433_zmm8r8(const __m512d k0,
                                            const __m512d a,
                                            const __m512d b,
                                            const __m512d phi1,
                                            const __m512d epsr,
                                            const __m512d epsi,
                                            const __m512d mur,
                                            const __m512d mui) {

                           return (rcs_f4431_zmm8r8(k0,a,b,phi1,epsr,epsi,mur,mui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4433_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsr,
                                              const double * __restrict __ATTR_ALIGN__(64) pepsi,
                                              const double * __restrict __ATTR_ALIGN__(64) pmur,
                                              const double * __restrict __ATTR_ALIGN__(64) pmui) {

                           return (rcs_f4431_zmm8r8_a(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               }


                   __ATTR_ALWAYS_INLINE__
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f4433_zmm8r8_u(const double * __restrict  pk0,
                                              const double * __restrict  pa,
                                              const double * __restrict  pb,
                                              const double * __restrict  pphi1,
                                              const double * __restrict  pepsr,
                                              const double * __restrict  pepsi,
                                              const double * __restrict  pmur,
                                              const double * __restrict  pmui) {

                           return (rcs_f4431_zmm8r8_u(pk0,pa,pb,pphi1,
                                                       pepsr,pepsi,pmur,pmui));
               }


                








                  

 



                   






                   

 




      } // radiolocation


} // gms









#endif /*__GMS_RCS_CYLINDER_ZMM8R8_HPP__*/
