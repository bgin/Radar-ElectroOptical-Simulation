
#ifndef __GMS_RCS_CYLINDER_ZMM16R4_HPP__
#define __GMS_RCS_CYLINDER_ZMM16R4_HPP__ 200120231636


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

    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_FULLVER =
      1000U*GMS_RCS_CYLINDER_ZMM16R4_MAJOR+
      100U*GMS_RCS_CYLINDER_ZMM16R4_MINOR+
      10U*GMS_RCS_CYLINDER_ZMM16R4_MICRO;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_CREATION_DATE = "20-01-2023 16:36 PM +00200 (FRI 20 JAN 2023 GMT+2)";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_DESCRIPTION   = "AVX512 optimized Cylinder Radar Cross Section (analytic) functionality."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"


namespace gms {


          namespace radiolocation {




               namespace {
                   const  static __m512 Ir  = _mm512_setzero_ps();
                   const  static __m512 Ii  = _mm512_set1_ps(1.0f);
                   const  static __m512 nIr = _mm512_set1_ps(-0.0f);
                   const  static __m512 nIi = _mm512_set1_ps(-1.0f);
                   const  static __m512 PI  = _mm512_set1_ps(3.14159265358979323846264338328f);

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
                   __m512 rcs_f419_zmm16r4(const __m512 a,
                                           const __m512 k0a) {

                          const register __m512 num = _mm512_mul_ps(a, 
                                                           _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 pi4 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const register __m512 c0  = _mm512_set1_ps(0.8905f);
                          const register __m512 arg = _mm512_mul_ps(k0a,c0);
                          __m512 ln,ln2,rcs,den;
                          ln = logkf(arg);
                          ln2= _mm512_mul_ps(ln,ln);
                          den= _mm512_fmadd_ps(k0a,ln2,pi4);
                          rcs= _mm512_div_ps(num,den);
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f419_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 a   = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const register __m512 num = _mm512_mul_ps(a, 
                                                           _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 pi4 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const register __m512 c0  = _mm512_set1_ps(0.8905f);
                          const register __m512 arg = _mm512_mul_ps(k0a,c0);
                          __m512 ln,ln2,rcs,den;
                          ln = logkf(arg);
                          ln2= _mm512_mul_ps(ln,ln);
                          den= _mm512_fmadd_ps(k0a,ln2,pi4);
                          rcs= _mm512_div_ps(num,den);
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f419_zmm16r4_u(const float * __restrict  pa,
                                             const float * __restrict  pk0a) {

                          const register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 num = _mm512_mul_ps(a, 
                                                           _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 pi4 = _mm512_set1_ps(2.467401100272339654708622749969f);
                          const register __m512 c0  = _mm512_set1_ps(0.8905f);
                          const register __m512 arg = _mm512_mul_ps(k0a,c0);
                          __m512 ln,ln2,rcs,den;
                          ln = logkf(arg);
                          ln2= _mm512_mul_ps(ln,ln);
                          den= _mm512_fmadd_ps(k0a,ln2,pi4);
                          rcs= _mm512_div_ps(num,den);
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
                   __m512 rcs_f4120_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          const register __m512 pi2a = _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 c0   = _mm512_set1_ps(2.25f);
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,
                                                                 _mm512_mul_ps(k0a,k0a));
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(c0,k0a3));
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4120_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 a    = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                          const register __m512 pi2a = _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 c0   = _mm512_set1_ps(2.25f);
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,
                                                                 _mm512_mul_ps(k0a,k0a));
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(c0,k0a3));
                          return (rcs);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4120_zmm16r4_u(const float * __restrict  pa,
                                            const float * __restrict  pk0a) {

                          const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 pi2a = _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 c0   = _mm512_set1_ps(2.25f);
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,
                                                                 _mm512_mul_ps(k0a,k0a));
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(c0,k0a3));
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
                   __m512 rcs_f4121_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          return (rcs_f4120_zmm16r4(a,k0a));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4121_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          return (rcs_f4120_zmm16r4_a(pa,pk0a));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4121_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          return (rcs_f4120_zmm16r4_u(pa,pk0a));
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
                   __m512 rcs_f4122_zmm16r4(const __m512 phi,
                                            const __m512 a,
                                            const __m512 k0a) {

                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_ps(hlf,cosph);
                          sqr   = _mm512_mul_ps(frac,frac);
                          rcs   = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4122_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 phi = _mm512_load_ps(&pphi[0]);
                          const register __m512 a   = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_ps(hlf,cosph);
                          sqr   = _mm512_mul_ps(frac,frac);
                          rcs   = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4122_zmm16r4_u(const float * __restrict  pphi,
                                              const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          const register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                          const register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 hlf = _mm512_set1_ps(0.5f);
                          register __m512 cosph,k0a3,frac,sqr,rcs; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cosph = xcosf(phi);
                          frac  = _mm512_add_ps(hlf,cosph);
                          sqr   = _mm512_mul_ps(frac,frac);
                          rcs   = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,sqr));
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
                   __m512 rcs_f4123_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          return (rcs_f4120_zmm16r4(a,k0a));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4123_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          return (rcs_f4120_zmm16r4_a(pa,pk0a));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4123_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          return (rcs_f4120_zmm16r4_u(pa,pk0a));
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
                   __m512 rcs_f4124_zmm16r4(const __m512 a,
                                            const __m512 k0a) {

                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const register __m512 qtr = _mm512_set1_ps(0.25f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,qtr));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4124_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512 a   = _mm512_load_ps(&pa[0]);
                          const register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const register __m512 qtr = _mm512_set1_ps(0.25f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,qtr));
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4124_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          const register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          const register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const register __m512 pi2a= _mm512_mul_ps(a, 
                                                    _mm512_set1_ps(9.869604401089358618834490999876f));
                          const register __m512 k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const register __m512 qtr = _mm512_set1_ps(0.25f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(pi2a,_mm512_mul_ps(k0a3,qtr));
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
                   void Kz_f4125_zmm16r4(const __m512 eps0,
                                         const __m512 mu0,
                                         const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 k0a,
                                         __m512 * __restrict Kzr,
                                         __m512 * __restrict Kzi) {

                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0a,ln);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&Kzr,&Kzi);
                      
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) peps0,
                                           const  float * __restrict __ATTR_ALIGN__(64) pmu0,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const   float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) Kzr,
                                           float * __restrict __ATTR_ALIGN__(64) Kzi) {

                        const register __m512 eps0 = _mm512_load_ps(&peps0[0]);
                        const register __m512 mu0  = _mm512_load_ps(&pmu0[0]);
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0a,ln);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_ps(&Kzr[0], resr);
                        _mm512_store_ps(&Kzi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kz_f4125_zmm16r4_u(const  float * __restrict  peps0,
                                           const  float * __restrict  pmu0,
                                           const   float * __restrict  pEr,
                                           const   float * __restrict  pEi,
                                           const   float * __restrict  pk0a,
                                           float * __restrict  Kzr,
                                           float * __restrict  Kzi) {

                        const register __m512 eps0 = _mm512_loadu_ps(&peps0[0]);
                        const register __m512 mu0  = _mm512_loadu_ps(&pmu0[0]);
                        const register __m512 Er   = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_loadu_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0a,ln);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(Er,Ei,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_storeu_ps(&Kzr[0], resr);
                        _mm512_storeu_ps(&Kzi[0], resi);
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
                   void Kph_f4126_zmm16r4(const __m512 Hr,
                                          const __m512 Hi,
                                          __m512 * __restrict Kphr,
                                          __m512 * __restrict Kphi) {

                        *Kphr = _mm512_mul_ps(nIi,Hr);
                        *Kphi = _mm512_mul_ps(nIi,Hi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) Hr,
                                            const float * __restrict __ATTR_ALIGN__(64) Hi,
                                           float * __restrict __ATTR_ALIGN__(64) Kphr,
                                          float * __restrict __ATTR_ALIGN__(64) Kphi) {

                        _mm512_store_ps(&Kphr[0] ,_mm512_mul_ps(nIi,_mm512_load_ps(&Hr[0]));
                        _mm512_store_ps(&Kphi[0] ,_mm512_mul_ps(nIi,_mm512_load_ps(&Hi[0]));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Kph_f4126_zmm16r4_u(const float * __restrict  Hr,
                                            const float * __restrict  Hi,
                                           float * __restrict  Kphr,
                                          float * __restrict Kphi) {

                        _mm512_storeu_ps(&Kphr[0] ,_mm512_mul_ps(nIi,_mm512_loadu_ps(&Hr[0]));
                        _mm512_storeu_ps(&Kphi[0] ,_mm512_mul_ps(nIi,_mm512_loadu_ps(&Hi[0]));
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
                   void Iz_f4127_zmm16r4(const __m512 eps0,
                                         const __m512 mu0,
                                         const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 k0a,
                                         const __m512 k0,
                                         __m512 * __restrict Izr,
                                         __m512 * __restrict Izi) {

                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 _2pi= _mm512_set1_ps(6.283185307179586476925286766559f);
                        const register __m512 ln  = logkf(lna);
                        const __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t2r = _mm512_mul_ps(_2pi,Er);
                        t1r = _mm512_mul_ps(k0,ln);
                        t2i = _mm512_mul_ps(_2pi,Ei);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        cdiv_zmm16r4(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&Izr,&Izi);
                      
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) peps0,
                                           const  float * __restrict __ATTR_ALIGN__(64) pmu0,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const   float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const   float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const   float * __restrict __ATTR_ALIGN__(64) pk0,
                                           float * __restrict __ATTR_ALIGN__(64) Izr,
                                           float * __restrict __ATTR_ALIGN__(64) Izi) {

                        const register __m512 eps0 = _mm512_load_ps(&peps0[0]);
                        const register __m512 _2pi= _mm512_set1_ps(6.283185307179586476925286766559f);
                        const register __m512 mu0  = _mm512_load_ps(&pmu0[0]);
                        const register __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 ln  = logkf(lna);
                        const register __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0,ln);
                        t2r = _mm512_mul_ps(_2pi,Er);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        t2i = _mm512_mul_ps(_2pi,Ei);
                        cdiv_zmm16r4(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_ps(&Izr[0], resr);
                        _mm512_store_ps(&Izi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Iz_f4127_zmm16r4_u(const  float * __restrict  peps0,
                                           const  float * __restrict  pmu0,
                                           const   float * __restrict  pEr,
                                           const   float * __restrict  pEi,
                                           const   float * __restrict  pk0a,
                                           const   float * __restrict  pk0,
                                           float * __restrict  Izr,
                                           float * __restrict  Izi) {

                        const register __m512 eps0 = _mm512_load_ps(&peps0[0]);
                        const register __m512 _2pi= _mm512_set1_ps(6.283185307179586476925286766559f);
                        const register __m512 mu0  = _mm512_load_ps(&pmu0[0]);
                        const register __m512 sqr = _mm512_sqrt_ps(_mm512_div_ps(eps0,mu0));
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 lna = _mm512_mul_ps(k0a,
                                                     _mm512_set1_ps(0.8905f));
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 ln  = logkf(lna);
                        const register __m512 pi2 = _mm512_set1_ps(1.57079632679489661923132169164f);
                        __m512 t0r,t0i,t1r,t1i,divr,divi,resr,resi,t2r,t2i;
                        t0r = nIr;
                        t0i = _mm512_mul_ps(nIi,sqr);
                        t1r = _mm512_mul_ps(k0,ln);
                        t2r = _mm512_mul_ps(_2pi,Er);
                        t1i = _mm512_mul_ps(nIi,pi2);
                        t2i = _mm512_mul_ps(_2pi,Ei);
                        cdiv_zmm16r4(t2r,t2i,t1r,t1i,&divr,&divi);
                        cmul_zmm16r4(t0r,t0i,divr,divi,&resr,&resi);
                        _mm512_store_ps(&Izr[0], resr);
                        _mm512_store_ps(&Izi[0], resi);
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
                   void EO_f4129_zmm16r4(const __m512 phi2,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 Er,
                                         const __m512 Ei,
                                         __m512 * __restrict EOr,
                                         __m512 * __restrict EOi) {

                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,exr,exi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(5.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_sub_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_ps(c6,t2i);
                        cmul_zmm16r4(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_ps(t4,_mm512_add_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(exr,exi,t2r,t2i,&t0r,&t0i);
                        *EOr = t0r;
                        *EOi = t0i;
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           float * __restrict __ATTR_ALIGN__(64) EOr,
                                           float * __restrict __ATTR_ALIGN__(64) EOi) {

                        const register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,exr,exi;;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(5.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_sub_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_ps(c6,t2i);
                        cmul_zmm16r4(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_ps(t4,_mm512_add_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(exr,exi,t2r,t2i,&t0r,&t0i);
                        _mm512_store_ps(&EOr[0], t0r);
                        _mm512_store_ps(&EOi[0], t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4129_zmm16r4_u(const float * __restrict  pphi2,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           float * __restrict  EOr,
                                           float * __restrict  EOi) {

                        const register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_loadu_ps(&pEi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,exr,exi;;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(5.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_sub_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_add_ps(c6,t2i);
                        cmul_zmm16r4(Er,Ei,t1r,t1i,&exr,&exi);
                        t5     = _mm512_mul_ps(t4,_mm512_add_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(exr,exi,t2r,t2i,&t0r,&t0i);
                        _mm512_storeu_ps(&EOr[0], t0r);
                        _mm512_storeu_ps(&EOi[0], t0i);
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
                   void HO_f4131_zmm16r4(const __m512 phi2,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 Hr,
                                         const __m512 Hi,
                                         __m512 * __restrict HOr,
                                         __m512 * __restrict HOi) {

                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(7.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_add_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_ps(c6,t2i);
                        cmul_zmm16r4(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_ps(t4,_mm512_sub_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        *HOr = t0r;
                        *HOi = t0i;
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pHr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHi,
                                           float * __restrict __ATTR_ALIGN__(64) HOr,
                                           float * __restrict __ATTR_ALIGN__(64) HOi) {

                        const register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_load_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_load_ps(&pHi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(7.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_add_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_ps(c6,t2i);
                        cmul_zmm16r4(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_ps(t4,_mm512_sub_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        _mm512_store_ps(&HOr[0], t0r);
                        _mm512_store_ps(&HOi[0], t0i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4131_zmm16r4_u(const float * __restrict  pphi2,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pHr,
                                           const float * __restrict  pHi,
                                           float * __restrict  HOr,
                                           float * __restrict  HOi) {

                        const register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_loadu_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_loadu_ps(&pHi[0]);
                        register __m512 t0,t1,t2,cosf2,cos2f2,t3,t4,t5;
                        register __m512 t0r,t0i,t1r,t1i,_2k0a,_2r,cos4f2;
                        register __m512 k0as,fac,earg,_2a,cer,cei,hxr,hxi;
                        register __m512 t2r,t2i;
                        const register __m512 c0 = _mm512_set1_ps(0.375f);
                        cosf2 = xcosf(phi2);
                        const register __m512 c1 = _mm512_set1_ps(0.1171875f);
                        cos2f2 = _mm512_mul_ps(cosf2,cosf2);
                        const register __m512 c2 = _mm512_set1_ps(4.0f);
                        cos4f2 = _mm512_mul_ps(cos2f2,cos2f2);
                        _2k0a  = _mm512_add_ps(k0a,k0a);
                        const register __m512 c3 = _mm512_set1_ps(8.0f);
                        _2r    = _mm512_add_ps(r,r);
                        _2a    = _mm512_add_ps(a,a);
                        const register __m512 c4 = _mm512_set1_ps(33.0f);
                        k0as   = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c5 = _mm512_set1_ps(7.0f);
                        t0     = _mm512_mul_ps(a,cosf2);
                        const register __m512 c6 = _mm512_set1_ps(1.0f);
                        t1     = _mm512_div_ps(t0,_2r);
                        fac    = _mm512_sqrt_ps(t1);
                        earg   = _mm512_mul_ps(k0,
                                          _mm512_sub_ps(r,_mm512_mul_ps(_2a,cosf2)));
                        t0r    = Ir;
                        t0i    = _mm512_mul_ps(Ii,earg);
                        cexp_zmm16r4(t0r,t0i,&cer,&cei);
                        t3     = _mm512_rcp14_ps(cos2f2);
                        cmul_zmm16r4(t0r,t0i,cer,cei,&t1r,&t1i);
                        t3     = _mm512_sub_ps(t3,c0);//taken t3
                        t0     = _mm512_mul_ps(c2,_mm512_mul_ps(k0as,cos2f2));
                        t4     = _mm512_rcp14_ps(t0);//taken t4
                        t1     = _mm512_mul_ps(c3,cos2f2);
                        t2     = _mm512_add_ps(c1,_mm512_div_ps(c4,t1)); // t2 taken
                        t0     = _mm512_div_ps(c5,cos4f2);// t0 taken
                        t2r    = Ir;
                        t2i    = _mm512_div_ps(Ii,_mm512_mul_ps(_2k0a,cosf2));
                        t2r    = c6;
                        t2i    = _mm512_sub_ps(c6,t2i);
                        cmul_zmm16r4(Hr,Hi,t1r,t1i,&hxr,&hxi);
                        t5     = _mm512_mul_ps(t4,_mm512_sub_ps(t2,t0));//taken t5
                        t2r    = _mm512_add_ps(t2r,t5);
                        t2i    = _mm512_add_ps(t2i,t5);
                        cmul_zmm16r4(hxr,hxi,t2r,t2i,&t0r,&t0i);
                        _mm512_storeu_ps(&HOr[0], t0r);
                        _mm512_storeu_ps(&HOi[0], t0i);
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
                   void EC_f4130_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 a,
                                         const

                   




                    



      } // radiolocation


} // gms









#endif /*__GMS_RCS_CYLINDER_ZMM16R4_HPP__*/
