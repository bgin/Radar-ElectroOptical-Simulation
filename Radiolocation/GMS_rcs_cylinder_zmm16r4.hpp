
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
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         __m512 * __restrict ECr,
                                         __m512 * __restrict ECi) {

                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(0.910721f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.9358135f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(1.607129f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.057397f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.0994145f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        *ECr = _mm512_add_ps(t0r,tmp3r);
                        *ECi = _mm512_add_ps(t0i,tmp3i);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) pphi,
                                           float * __restrict __ATTR_ALIGN__(64) ECr,
                                           float * __restrict __ATTR_ALIGN__(64) ECi) {

                        const register __m512 phi  = _mm512_load_ps(&pphi[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_load_ps(&pEi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(0.910721f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.9358135f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(1.607129f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.057397f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.0994145f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_store_ps(&ECr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_store_ps(&ECi[0] ,_mm512_add_ps(t0i,tmp3i));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4130_zmm16r4_a(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict  pphi,
                                           float * __restrict  ECr,
                                           float * __restrict ECi) {

                        const register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Er   = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei   = _mm512_loadu_ps(&pEi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(0.910721f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.9358135f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(1.607129f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.057397f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.0994145f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_storeu_ps(&ECr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_storeu_ps(&ECi[0] ,_mm512_add_ps(t0i,tmp3i));
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
                   void HC_f4132_zmm16r4(const __m512 Hr,
                                         const __m512 Hi,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         __m512 * __restrict HCr,
                                         __m512 * __restrict HCi) {

                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(1.531915f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.404308f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(0.70028f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.072732f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.1259755f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        *HCr = _mm512_add_ps(t0r,tmp3r);
                        *HCi = _mm512_add_ps(t0i,tmp3i);
                 }


                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64)  pHr,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pHi,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pa,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const  float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           float * __restrict __ATTR_ALIGN__(64)  HCr,
                                           float * __restrict __ATTR_ALIGN__(64)  HCi) {

                        const register __m512 phi  = _mm512_load_ps(&pphi[0]);
                        const register __m512 a    = _mm512_load_ps(&pa[0]);
                        const register __m512 r    = _mm512_load_ps(&pr[0]);
                        const register __m512 k0   = _mm512_load_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_load_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_load_ps(&pHi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(1.531915f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.404308f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(0.70028f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.072732f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.1259755f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_store_ps(&HCr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_store_ps(&HCi[0] ,_mm512_add_ps(t0i,tmp3i));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4132_zmm16r4_u(const  float * __restrict  pHr,
                                           const  float * __restrict  pHi,
                                           const  float * __restrict  pa,
                                           const  float * __restrict  pr,
                                           const  float * __restrict  pk0,
                                           const  float * __restrict  pk0a,
                                           const  float * __restrict  pphi,
                                           float * __restrict   HCr,
                                           float * __restrict   HCi) {

                        const register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                        const register __m512 a    = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r    = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0   = _mm512_loadu_ps(&pk0[0]);
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 Hr   = _mm512_loadu_ps(&pHr[0]);
                        const register __m512 Hi   = _mm512_loadu_ps(&pHi[0]);
                        register __m512 e1ar,e1ai,ce1r,ce1i,e2ar,e2ai,e3ar,e3ai;
                        register __m512 ce2r,ce2i,ce3r,ce3i,sqr,tmp2r,tmp2i;
                        register __m512 Etr,Eti,tmp1r,tmp1i,t0r,t0i,t1r,t1i,tmp3r,tmp3i;
                        const  __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        const  __m512 k0ai16 = _mm512_pow_ps(k0a,
                                                         _mm512_set1_ps(0.166666666666666666666666666667f));
                        const register __m512 k0apaphi = _mm512_fmadd_ps(k0a,PI,phi);
                        e2ar   = Ir; 
                        e2ai   = k0apaphi;
                        cexp_zmm16r4(e2ar,e2ai,&ce2r,&ce2i); // second cexp
                        k0ai16 = _mm512_rcp14_ps(k0ai16);
                        const register __m512 k0apsphi = _mm512_fmsub_ps(k0a,PI,phi);
                        const  __m512 c0   = _mm512_set1_ps(1.531915f);
                        const  register __m512 k0rp12 = _mm512_fmadd_ps(k0,r,pi12);
                        e3ar   = Ir;
                        e3ai   = k0apsphi;
                        cexp_zmm16r4(e3ar,e3ai,&ce3r,&ce3i); // third cexp
                        const  __m512 c0r  = _mm512_set1_ps(0.404308f);
                        sqr    = _mm512_div_ps(a,_mm512_add_ps(r,r);
                        const  __m512 c0i  = _mm512_set1_ps(0.70028f);
                        sqr    = _mm512_sqrt_ps(sqr);
                        const  __m512 c1r  = _mm512_set1_ps(0.072732f);
                        Etr    = _mm512_mul_ps(Er,sqr);// first complex term
                        Eti    = _mm512_mul_ps(Ei,sqr);// first complex term
                        const  __m512 c1i  = _mm512_set1_ps(-0.1259755f);
                        const  __m512 _1   = _mm512_set1_ps(1.0f);
                        const  __m512 k0an23 = _mm512_pow_ps(k0a,
                                                          _mm512_set1_ps(0.666666666666666666666666666667f));
                        k0an23 = _mm512_rcp14_ps(k0an23);
                        const __m512 k0an43= _mm512_pow_ps(k0a,
                                                        _mm512_set1_ps(1.333333333333333333333333333333f));
                        k0an43 = _mm512_rcp14_ps(k0an43);
                        e1ar = Ir;
                        e1ai = k0rp12;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i); // first cexp
                        t0r = _mm512_fmadd_ps(_1,c0r,k0an23);// 1st complex term
                        t0i = _mm512_fmadd_ps(_1,c0i,k0an23);// //-//-//-
                        t1r = _mm512_mul_ps(c1r,k0an43);
                        t1i = _mm512_mul_ps(c1i,k0an43);
                        tmp1r = _mm512_sub_ps(t0r,t1r);
                        tmp1i = _mm512_sub_ps(t0i,t1i);
                        cmul_zmm16r4(ce2r,ce2i,tmp1r,tmp1i,&tmp2r,&tmp2i);
                        cmul_zmm16r4(ce3r,ce3i,tmp1r,tmp1i,&tmp3r,&tmp3i);
                        t0r = _mm512_setzero_ps();
                        t0i = _mm512_setzero_ps();
                        cmul_zmm16r4(Etr,Eti,tmp2r,tmp2i,&t0r,&t0i)
                        _mm512_storeu_ps(&HCr[0] ,_mm512_add_ps(t0r,tmp3r));
                        _mm512_storeu_ps(&HCi[0] ,_mm512_add_ps(t0i,tmp3i));
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
                   void EO_f4133_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0a,
                                         __m512 * __restrict EOr,
                                         __m512 * __restrict EOi) {

                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _5 = _mm512_set1_ps(5.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(127.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,*EOr,*EOi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) EOr,
                                           float * __restrict __ATTR_ALIGN__(64) EOi) {

                        register __m512 Er = _mm512_load_ps(&pEr[0]);
                        register __m512 Ei = _mm512_load_ps(&pEi[0]);
                        register __m512 a  = _mm512_load_ps(&pa[0]);
                        register __m512 r  = _mm512_load_ps(&pr[0]);
                        register __m512 k0a= _mm512_load_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        register __m512 resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _5 = _mm512_set1_ps(5.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(127.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&EOr[0],resr);
                        _mm512_store_ps(&EOi[0],resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EO_f4133_zmm16r4_u(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0a,
                                           float * __restrict  EOr,
                                           float * __restrict  EOi) {

                        register __m512 Er = _mm512_loadu_ps(&pEr[0]);
                        register __m512 Ei = _mm512_loadu_ps(&pEi[0]);
                        register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0a= _mm512_loadu_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        register __m512 resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _5 = _mm512_set1_ps(5.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(127.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_5,t1);
                        t0r             = _1;
                        t0i             = _mm512_add_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&EOr[0],resr);
                        _mm512_storeu_ps(&EOi[0],resi);
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
                   void HO_f4135_zmm16r4(const __m512 Hr,
                                         const __m512 Hi,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0a,
                                         __m512 * __restrict HOr,
                                         __m512 * __restrict HOi) {

                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _11 = _mm512_set1_ps(11.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(353.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,*HOr,*HOi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pHr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) HOr,
                                           float * __restrict __ATTR_ALIGN__(64) HOi) {

                        register __m512 Hr = _mm512_load_ps(&pHr[0]);
                        register __m512 Hi = _mm512_load_ps(&pHi[0]);
                        register __m512 a  = _mm512_load_ps(&pa[0]);
                        register __m512 r  = _mm512_load_ps(&pr[0]);
                        register __m512 k0a= _mm512_load_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i,resr,resi;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _11 = _mm512_set1_ps(11.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(353.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&HOr[0], resr);
                        _mm512_store_ps(&HOi[0], resi);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HO_f4135_zmm16r4_u(const float * __restrict  pHr,
                                           const float * __restrict  pHi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0a,
                                           float * __restrict  HOr,
                                           float * __restrict  HOi) {

                        register __m512 Hr = _mm512_loadu_ps(&pHr[0]);
                        register __m512 Hi = _mm512_loadu_ps(&pHi[0]);
                        register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0a= _mm512_loadu_ps(&pk0a[0]);
                        register __m512 _2r,_2k0a,k0as,k0r,t0,t1,t2,t3,t1r,t1i,resr,resi;
                        register __m512 ear,eai,cer,cei,facr,faci,t0r,t0i,cer1,cei1;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        const __m512 _16= _mm512_set1_ps(16.0f);
                        _2r             = _mm512_add_ps(r,r);
                        _2k0a           = _mm512_add_ps(k0a,k0a);
                        const __m512 _11 = _mm512_set1_ps(11.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 c0 = _mm512_set1_ps(353.0f);
                        t0              = _mm512_sub_ps(k0r,_2k0a);
                        ear             = Ir;
                        const __m512 c1 = _mm512_set1_ps(512.0f);
                        eai             = t0;
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        t1              = _mm512_div_ps(a,_2r);
                        t2              = _mm512_sqrt_ps(t1);
                        facr            = _mm512_mul_ps(Er,t2);
                        faci            = _mm512_mul_ps(Ei,t2);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(_16,k0a);
                        t2              = _mm512_mul_ps(c1,k0as);
                        t0r             = Ir;
                        t0i             = _mm512_div_ps(_11,t1);
                        t0r             = _1;
                        t0i             = _mm512_sub_ps(_1,t0i);
                        t3              = _mm512_div_ps(c0,t2);
                        t0r             = _mm512_add_ps(t3,t0r);
                        t0i             = _mm512_add_ps(t3,t0i);
                        cmul_zmm16r4(facr,faci,cer,cei,&cer1,&cei1);
                        cmul_zmm16r4(cer1,cei1,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&HOr[0], resr);
                        _mm512_storeu_ps(&HOi[0], resi);
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
                   void EC_f4134_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         __m512 * __restrict ECr,
                                         __m512 * __restrict ECi) {

                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(2.939945f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.180318f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(1.821442f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-5.048945f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.312320f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Er,t2);
                        fraci = _mm512_mul_ps(Ei,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        *ECr = _mm512_mul_ps(t0r,rex);
                        *ECi = _mm512_mul_ps(t0i,rex);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                         const float * __restrict __ATTR_ALIGN__(64) pEi,
                                         const float * __restrict __ATTR_ALIGN__(64) pa,
                                         const float * __restrict __ATTR_ALIGN__(64) pr,
                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
                                         float * __restrict __ATTR_ALIGN__(64) ECr,
                                         float * __restrict __ATTR_ALIGN__(64) ECi) {

                        const register __m512 Er = _mm512_load_ps(&pEr[0]);
                        const register __m512 Ei = _mm512_load_ps(&pEi[0]);
                        const register __m512 a  = _mm512_load_ps(&pa[0]);
                        const register __m512 r  = _mm512_load_ps(&pr[0]);
                        const register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(2.939945f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.180318f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(1.821442f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-5.048945f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.312320f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Er,t2);
                        fraci = _mm512_mul_ps(Ei,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_store_ps(&ECr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_store_ps(&ECi[0] ,_mm512_mul_ps(t0i,rex));
                }


                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EC_f4134_zmm16r4_u(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           float * __restrict  ECr,
                                           float * __restrict  ECi) {

                        const register __m512 Er = _mm512_loadu_ps(&pEr[0]);
                        const register __m512 Ei = _mm512_loadu_ps(&pEi[0]);
                        const register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(2.939945f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.180318f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(1.821442f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-5.048945f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.312320f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Er,t2);
                        fraci = _mm512_mul_ps(Ei,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_storeu_ps(&ECr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_storeu_ps(&ECi[0] ,_mm512_mul_ps(t0i,rex));
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
                   void HC_f4136_zmm16r4(const __m512 Hr,
                                         const __m512 Hi,
                                         const __m512 a,
                                         const __m512 r,
                                         const __m512 k0,
                                         __m512 * __restrict HCr,
                                         __m512 * __restrict HCi) {

                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(1.2701695f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.2284945f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(3.063830f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-2.200000f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.3957635f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Hr,t2);
                        fraci = _mm512_mul_ps(Hi,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        *HCr = _mm512_mul_ps(t0r,rex);
                        *HCi = _mm512_mul_ps(t0i,rex);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pHr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHi,
                                           const float * __restrict __ATTR_ALIGN__(64) pa,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           float * __restrict __ATTR_ALIGN__(64) HCr,
                                           float * __restrict __ATTR_ALIGN__(64) HCi) {

                        const register __m512 Hr = _mm512_load_ps(&pHr[0]);
                        const register __m512 Hi = _mm512_load_ps(&pHi[0]);
                        const register __m512 a  = _mm512_load_ps(&pa[0]);
                        const register __m512 r  = _mm512_load_ps(&pr[0]);
                        const register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(1.2701695f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.2284945f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(3.063830f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-2.200000f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.3957635f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Hr,t2);
                        fraci = _mm512_mul_ps(Hi,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_store_ps(&HCr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_store_ps(&HCi[0] ,_mm512_mul_ps(t0i,rex));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void HC_f4136_zmm16r4_u(const float * __restrict  pHr,
                                           const float * __restrict  pHi,
                                           const float * __restrict  pa,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           float * __restrict  HCr,
                                           float * __restrict  HCi) {

                        const register __m512 Hr = _mm512_loadu_ps(&pHr[0]);
                        const register __m512 Hi = _mm512_loadu_ps(&pHi[0]);
                        const register __m512 a  = _mm512_loadu_ps(&pa[0]);
                        const register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        const register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0r,k0a,k0a13,k0an13,k0an16;
                        register __m512 fracr,fraci,_2r,t0,t1,t2;
                        register __m512 e1ar,e1ai,exar;
                        register __m512 ce1r,ce1i,rex,t0r,t0i;
                        const __m512 pi12 = _mm512_set1_ps(0.261799387799149436538553615273f);
                        k0r   = _mm512_mul_ps(k0,r);
                        const __m512 c0   = _mm512_set1_ps(1.2701695f);
                        k0a   = _mm512_mul_ps(k0,a);
                        const __m512 c1   = _mm512_set1_ps(0.2284945f); 
                        k0a13 = _mm512_pow_ps(k0a,
                                         _mm512_set1_ps(0.333333333333333333333333333333333333f));
                        const __m512 c2   = _mm512_set1_ps(3.063830f);
                        _2r   = _mm512_add_ps(r,r);
                        const __m512 c3   = _mm512_set1_ps(-2.200000f);
                        k0an13= _mm512_rcp14_ps(k0a13);
                        const __m512 c4   = _mm512_set1_ps(0.3957635f);
                        t0    = _mm512_div_ps(a,_2r);
                        t1    = _mm512_pow_ps(k0a,
                                          _mm512_set1_ps(0.166666666666666666666666666667f));
                        k0an16= _mm512_rcp14_ps(t1);
                        t2    = _mm512_sqrt_ps(t0);
                        fracr = _mm512_mul_ps(Hr,t2);
                        fraci = _mm512_mul_ps(Hi,t2);
                        t0    = _mm512_fmsub_ps(c0,k0a13,
                                            _mm512_mul_ps(c1,k0an13));
                        t1    = _mm512_fmadd_ps(k0a,PI,_mm512_add_ps(pi12,t0));
                        t1    = _mm512_add_ps(k0r,t1);
                        e1ar  = Ir;
                        e1ai  = t1;
                        cexp_zmm16r4(e1ar,e1ai,&ce1r,&ce1i);
                        exar  = _mm512_fmsub_ps(c3,k0a13,
                                            _mm512_mul_ps(c4,k0an13));
                        t1    = _mm512_mul_ps(c2,k0an16);
                        t2    = xexpf(exar);
                        rex   = _mm512_rcp14_ps(t2);
                        rex   = _mm512_mul_ps(rex,t1);
                        cmul_zmm16r4(fracr,fraci,ce1r,ce1i,&t0r,&t0i);
                        _mm512_storeu_ps(&HCr[0] ,_mm512_mul_ps(t0r,rex));
                        _mm512_storeu_ps(&HCi[0] ,_mm512_mul_ps(t0i,rex));
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
                   __m512 rcs_f4137_zmm16r4(const __m512 a,
                                            const __m512 phi2) {

                          register __m512 rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_ps(PI,_mm512_mul_ps(a,cosp2));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4137_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pphi2) {

                          register __m512 a    = _mm512_load_ps(&pa[0]);
                          register __m512 phi2 = _mm512_load_ps(&pphi2[0]);
                          register __m512 rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_ps(PI,_mm512_mul_ps(a,cosp2));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4137_zmm16r4_u(const float * __restrict  pa,
                                              const float * __restrict  pphi2) {

                          register __m512 a    = _mm512_loadu_ps(&pa[0]);
                          register __m512 phi2 = _mm512_loadu_ps(&pphi2[0]);
                          register __m512 rcs,cosp2;
                          cosp2 = xcosf(phi2);
                          rcs   = _mm512_mul_ps(PI,_mm512_mul_ps(a,cosp2));
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
                   __m512 rcs_f4138_zmm16r4(const __m512 a) {

                          return (__m512_mul_ps(a,PI));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4138_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa) {

                          register __m512 a = _mm512_load_ps(&pa[0]);
                          return (__m512_mul_ps(a,PI));
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4138_zmm16r4_u(const float * __restrict  pa) {

                          register __m512 a = _mm512_loadu_ps(&pa[0]);
                          return (__m512_mul_ps(a,PI));
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
                   __m512 rcs_f4140_zmm16r4(const __m512 k0a,
                                            const __m512 alpha) {

                          register __m512 sinc,k0alp,k0as,t0;
                          register __m512 rcs;
                          const __m512 _4       = _mm512_set1_ps(4.0f);
                          k0alp = _mm512_mul_ps(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_ps(t0,k0alp); 
                          k0as  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a)); 
                          rcs   = _mm512_mul_ps(k0as,_mm512_mul_ps(sinc,sinc));
                          return (rcs);        
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4140_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(64) palpha) {

                          register __m512 k0a   = _mm512_load_ps(&pk0a[0]);
                          register __m512 alpha = _mm512_load_ps(&palpha[0]);
                          register __m512 sinc,k0alp,k0as,t0;
                          register __m512 rcs;
                          const __m512 _4       = _mm512_set1_ps(4.0f);
                          k0alp = _mm512_mul_ps(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_ps(t0,k0alp); 
                          k0as  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a)); 
                          rcs   = _mm512_mul_ps(k0as,_mm512_mul_ps(sinc,sinc));
                          return (rcs);        
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4140_zmm16r4_u(const float * __restrict  pk0a,
                                              const float * __restrict  palpha) {

                          register __m512 k0a   = _mm512_loadu_ps(&pk0a[0]);
                          register __m512 alpha = _mm512_loadu_ps(&palpha[0]);
                          register __m512 sinc,k0alp,k0as,t0;
                          register __m512 rcs;
                          const __m512 _4       = _mm512_set1_ps(4.0f);
                          k0alp = _mm512_mul_ps(k0a,alpha);
                          t0    = xsinf(k0alp);
                          sinc  = _mm512_div_ps(t0,k0alp); 
                          k0as  = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a)); 
                          rcs   = _mm512_mul_ps(k0as,_mm512_mul_ps(sinc,sinc));
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
                   void Es_f4139_zmm16r4(const __m512 Er,
                                         const __m512 Ei,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 alp
                                         const __m512 k0a,
                                         __m512 * __restrict Esr,
                                         __m512 * __restrict Esi) {

                        register __m512 _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512 facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512 cer,cei;
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        _2k0as = _mm512_add_ps(_2,_mm512_mul_ps(k0a,k0a));
                        k0r    = _mm512_mul_ps(k0,r);
                        k0alp  = _mm512_mul_ps(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_ps(k0r,pir);
                        sinc   = _mm512_div_ps(t0,k0alp);
                        cexp_zmm16r4(arr,ari,&cer,&cei);
                        div    = _mm512_div_ps(_2k0as,pi4);
                        t0     = _mm512_sqrt_ps(div);
                        facr   = _mm512_mul_ps(Er,t0);
                        t0r    = _mm512_mul_ps(cer,sinc);
                        faci   = _mm512_mul_ps(Ei,t0);
                        t0i    = _mm512_mul_ps(cei,sinc);
                        cmul_zmm16r4(facr,faci,t0r,t0i,*Esr,*Esi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEi,
                                           const float * __restrict __ATTR_ALIGN__(64) pr,
                                           const float * __restrict __ATTR_ALIGN__(64) pk0,
                                           const float * __restrict __ATTR_ALIGN__(64) palp
                                           const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           float * __restrict __ATTR_ALIGN__(64) Esr,
                                           float * __restrict __ATTR_ALIGN__(64) Esi) {

                        register __m512 Er = _mm512_load_ps(&pEr[0]);
                        register __m512 Ei = _mm512_load_ps(&pEi[0]);
                        register __m512 r  = _mm512_load_ps(&pr[0]);
                        register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 alp= _mm512_load_ps(&palp[0]);
                        register __m512 k0a= _mm512_load_ps(&pk0a[0]);
                        register __m512 _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512 facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512 cer,cei,resr,resi;
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        _2k0as = _mm512_add_ps(_2,_mm512_mul_ps(k0a,k0a));
                        k0r    = _mm512_mul_ps(k0,r);
                        k0alp  = _mm512_mul_ps(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_ps(k0r,pir);
                        sinc   = _mm512_div_ps(t0,k0alp);
                        cexp_zmm16r4(arr,ari,&cer,&cei);
                        div    = _mm512_div_ps(_2k0as,pi4);
                        t0     = _mm512_sqrt_ps(div);
                        facr   = _mm512_mul_ps(Er,t0);
                        t0r    = _mm512_mul_ps(cer,sinc);
                        faci   = _mm512_mul_ps(Ei,t0);
                        t0i    = _mm512_mul_ps(cei,sinc);
                        cmul_zmm16r4(facr,faci,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&Esr[0], resr);
                        _mm512_store_ps(&Esi[0], resi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4139_zmm16r4_u(const float * __restrict  pEr,
                                           const float * __restrict  pEi,
                                           const float * __restrict  pr,
                                           const float * __restrict  pk0,
                                           const float * __restrict  palp
                                           const float * __restrict  pk0a,
                                           float * __restrict  Esr,
                                           float * __restrict  Esi) {

                        register __m512 Er = _mm512_loadu_ps(&pEr[0]);
                        register __m512 Ei = _mm512_loadu_ps(&pEi[0]);
                        register __m512 r  = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 alp= _mm512_loadu_ps(&palp[0]);
                        register __m512 k0a= _mm512_loadu_ps(&pk0a[0]);
                        register __m512 _2k0as,k0r,k0alp,sinc,pir,div;
                        register __m512 facr,faci,arr,ari,t0r,t0i,t0;
                        register __m512 cer,cei,resr,resi;
                        const __m512 _2  = _mm512_set1_ps(2.0f);
                        const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                        _2k0as = _mm512_add_ps(_2,_mm512_mul_ps(k0a,k0a));
                        k0r    = _mm512_mul_ps(k0,r);
                        k0alp  = _mm512_mul_ps(k0a,alp);
                        t0     = xsinf(k0alp);
                        arr    = Ir;
                        ari    = _mm512_sub_ps(k0r,pir);
                        sinc   = _mm512_div_ps(t0,k0alp);
                        cexp_zmm16r4(arr,ari,&cer,&cei);
                        div    = _mm512_div_ps(_2k0as,pi4);
                        t0     = _mm512_sqrt_ps(div);
                        facr   = _mm512_mul_ps(Er,t0);
                        t0r    = _mm512_mul_ps(cer,sinc);
                        faci   = _mm512_mul_ps(Ei,t0);
                        t0i    = _mm512_mul_ps(cei,sinc);
                        cmul_zmm16r4(facr,faci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&Esr[0], resr);
                        _mm512_storeu_ps(&Esi[0], resi);
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
                   __m512 rcs_f4141_zmm16r4(const __m512 k0a) {

                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a));
                          return (rcs);
                 }

                

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4141_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64)  pk0a) {

                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4141_zmm16r4_u(const float * __restrict   pk0a) {

                          register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                          const __m512 _4 = _mm512_set1_ps(4.0f);
                          register __m512 rcs;
                          rcs = _mm512_mul_ps(_4,_mm512_mul_ps(k0a,k0a));
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
                   void Es_f4145_zmm16r4(const __m512 EIr,
                                         const __m512 EIi,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         const __m512 eps0,
                                         const __m512 eps1,
                                         const __m512 mu0,
                                         const __m512 mu1,
                                         __m512 * __restrict ESr,
                                         __m512 * __restrict ESi) {

                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(EIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                        fraci           = _mm512_mul_ps(EIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(mu1,mu0),
                                                        _mm512_add_ps(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,*ESr,*ESi);
               }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pEIr,
                                           const float * __restrict __ATTR_ALIGN__(64) pEIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps0,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps1,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu1,
                                           float * __restrict __ATTR_ALIGN__(64)  ESr,
                                           float * __restrict __ATTR_ALIGN__(64)  ESi) {

                        register __m512 EIr = _mm512_load_ps(&pEIr[0]);
                        register __m512 EIi = _mm512_load_ps(&pEIi[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 eps0= _mm512_load_ps(&peps0[0]);
                        register __m512 eps1= _mm512_load_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(EIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                        fraci           = _mm512_mul_ps(EIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(mu1,mu0),
                                                        _mm512_add_ps(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&ESr[0], resr);
                        _mm512_store_ps(&ESi[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Es_f4145_zmm16r4_u(const float * __restrict  pEIr,
                                           const float * __restrict  pEIi,
                                           const float * __restrict  pr,
                                           const float * __restrict   pk0,
                                           const float * __restrict  pk0a,
                                           const float * __restrict   pphi,
                                           const float * __restrict   peps0,
                                           const float * __restrict  peps1,
                                           const float * __restrict   pmu0,
                                           const float * __restrict   pmu1,
                                           float * __restrict   ESr,
                                           float * __restrict   ESi) {

                        register __m512 EIr = _mm512_loadu_ps(&pEIr[0]);
                        register __m512 EIi = _mm512_loadu_ps(&pEIi[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                        register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(EIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                        fraci           = _mm512_mul_ps(EIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(mu1,mu0),
                                                        _mm512_add_ps(mu1,mu0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&ESr[0], resr);
                        _mm512_storeu_ps(&ESi[0], resi);
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
                   void Hs_f4146_zmm16r4(const __m512 HIr,
                                         const __m512 HIi,
                                         const __m512 r,
                                         const __m512 k0,
                                         const __m512 k0a,
                                         const __m512 phi,
                                         const __m512 eps0,
                                         const __m512 eps1,
                                         const __m512 mu0,
                                         const __m512 mu1,
                                         __m512 * __restrict HSr,
                                         __m512 * __restrict HSi) {

                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(HIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                        fraci           = _mm512_mul_ps(HIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(eps1,eps0),
                                                        _mm512_add_ps(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,*HSr,*HSi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pHIr,
                                           const float * __restrict __ATTR_ALIGN__(64) pHIi,
                                           const float * __restrict __ATTR_ALIGN__(64)  pr,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64)  pphi,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps0,
                                           const float * __restrict __ATTR_ALIGN__(64)  peps1,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu0,
                                           const float * __restrict __ATTR_ALIGN__(64)  pmu1,
                                           float * __restrict __ATTR_ALIGN__(64)  HSr,
                                           float * __restrict __ATTR_ALIGN__(64)  HSi) {

                        register __m512 HIr = _mm512_load_ps(&pHIr[0]);
                        register __m512 HIi = _mm512_load_ps(&pHIi[0]);
                        register __m512 r   = _mm512_load_ps(&pr[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);
                        register __m512 eps0= _mm512_load_ps(&peps0[0]);
                        register __m512 eps1= _mm512_load_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_load_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_load_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(HIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                        fraci           = _mm512_mul_ps(HIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(eps1,eps0),
                                                        _mm512_add_ps(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_store_ps(&HSr[0], resr);
                        _mm512_store_ps(&HSi[0], resi);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void Hs_f4146_zmm16r4_u(const float * __restrict   pHIr,
                                           const float * __restrict   pHIi,
                                           const float * __restrict   pr,
                                           const float * __restrict   pk0,
                                           const float * __restrict   pk0a,
                                           const float * __restrict   pphi,
                                           const float * __restrict   peps0,
                                           const float * __restrict   peps1,
                                           const float * __restrict   pmu0,
                                           const float * __restrict   pmu1,
                                           float * __restrict  HSr,
                                           float * __restrict   HSi) {

                        register __m512 HIr = _mm512_loadu_ps(&pHIr[0]);
                        register __m512 HIi = _mm512_loadu_ps(&pHIi[0]);
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);
                        register __m512 eps0= _mm512_loadu_ps(&peps0[0]);
                        register __m512 eps1= _mm512_loadu_ps(&peps1[0]);
                        register __m512 mu0 = _mm512_loadu_ps(&pmu0[0]);
                        register __m512 mu1 = _mm512_loadu_ps(&pmu1[0]);
                        register __m512 k0r,k0as,fracr,fraci,k0as2;
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 t0,t1,cosp,t2,sk0r,t3,mul,resr,resi;
                        const __m512 _1 = _mm512_set1_ps(1.0f);
                        k0r             = _mm512_mul_ps(k0,r);
                        const __m512 _2 = _mm512_set1_ps(2.0f);
                        sk0r            = _mm512_sqrt_ps(k0r);
                        k0as            = _mm512_mul_ps(k0a,k0a);
                        const __m512 hlf= _mm512_set1_ps(0.5f);
                        k0as2           = _mm512_mul_ps(hlf,k0as);
                        const __m512 pi4= _mm512_set1_ps(0.78539816339744830961566084582f);
                        cosp            = xcosf(phi);
                        const __m512 pi2= _mm512_set1_ps(1.253314137315500251207882642406f);
                        fracr           = _mm512_mul_ps(HIr,pi2);
                        t0              = _mm512_sub_ps(_mm512_div_ps(mu1,mu0),_1);
                        fraci           = _mm512_mul_ps(HIi,pi2);
                        t1              = _mm512_div_ps(_mm512_sub_ps(eps1,eps0),
                                                        _mm512_add_ps(eps1,eps0));
                        ear             = Ir;
                        t2              = _mm512_add_ps(t1,t1);
                        eai             = _mm512_sub_ps(k0r,pi4);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        t1              = _mm512_mul_ps(t2,cosp);
                        cer = _mm512_div_ps(cer,sk0r);
                        t3  = _mm512_sub_ps(t0,t1);
                        cei = _mm512_div_ps(cei,sk0r);
                        mul = _mm512_mul_ps(k0as2,t3);
                        t0r = _mm512_mul_ps(cer,mul);
                        t0i = _mm512_mul_ps(cei,mul);
                        cmul_zmm16r4(fracr,fraci,t0r,t0i,&resr,&resi);
                        _mm512_storeu_ps(&HSr[0], resr);
                        _mm512_storeu_ps(&HSi[0], resi);
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
                   __m512 rcs_f4147_zmm16r4(const __m512 a,
                                            const __m512 k0a,
                                            const __m512 phi,
                                            const __m512 eps1,
                                            const __m512 eps0,
                                            const __m512 mu1,
                                            const __m512 mu0) {

                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }


                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f4147_zmm16r4(const float * __restrict __ATTR_ALIGN__(64) pa,
                                            const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(64) pphi,
                                            const float * __restrict __ATTR_ALIGN__(64) peps1,
                                            const float * __restrict __ATTR_ALIGN__(64) peps0,
                                            const float * __restrict __ATTR_ALIGN__(64) pmu1,
                                            const float * __restrict __ATTR_ALIGN__(64) pmu0) {

                          register __m512 a = _mm512_load_ps(&pa[0]);
                          register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                          register __m512 phi = _mm512_load_ps(&pphi[0]);
                          register __m512 eps1= _mm512_load_ps(&peps1[0]);
                          register __m512 eps0= _mm512_load_ps(&peps0[0]);
                          register __m512 pmu1= _mm512_load_ps(&pmu1[0]);
                          register __m512 pmu0= _mm512_load_ps(&pmu0[0]);
                          register __m512 t0,t1,k0a3,epst,mut,cosp,sqr,t2,diff;
                          register __m512 rcs;
                          cosp             = xcosf(phi);
                          const __m512 pi4 = _mm512_set1_ps(0.78539816339744830961566084582f);
                          t0               = _mm512_mul_ps(pi4,_mm512_mul_ps(PI,a));
                          const __m512 _1  = _mm512_set1_ps(1.0f);
                          k0a3             = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          const __m512 _2  = _mm512_set1_ps(2.0f);
                          epst             = _mm512_sub_ps(_mm512_div_ps(eps1,eps0),_1);
                          t1               = _mm512_sub_ps(mu1,mu0);
                          t2               = _mm512_add_ps(mu1,mu0);
                          mut              = _mm512_mul_ps(_2,_mm512_div_ps(t1,t2));
                          diff             = _mm512_sub_ps(epst,_mm512_mul_ps(mut,cosp));
                          sqr              = _mm512_mul_ps(diff,diff);
                          rcs              = _mm512_mul_ps(t0,_mm512_mul_ps(k0a3,sqr));
                          return (rcs);
                }                    


                    



      } // radiolocation


} // gms









#endif /*__GMS_RCS_CYLINDER_ZMM16R4_HPP__*/
