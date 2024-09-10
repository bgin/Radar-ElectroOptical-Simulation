
#ifndef __GMS_RCS_SPHERE_YMM8R4_HPP__
#define __GMS_RCS_SPHERE_YMM8R4_HPP__  100920240643


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

    const unsigned int GMS_RCS_SPHERE_YMM8R4_MAJOR = 1U;
    const unsigned int GMS_RCS_SPHERE_YMM8R4_MINOR = 0U;
    const unsigned int GMS_RCS_SPHERE_YMM8R4_MICRO = 0U;
    const unsigned int GMS_RCS_SPHERE_YMM8R4_FULLVER =
      1000U*GMS_RCS_SPHERE_YMM8R4_MAJOR+
      100U*GMS_RCS_SPHERE_YMM8R4_MINOR+
      10U*GMS_RCS_SPHERE_YMM8R4_MICRO;
    const char * const GMS_RCS_SPHERE_YMM8R4_CREATION_DATE = "10-09-2024 06:43 PM +00200 (TUE 10 SEP 2024 GMT+2)";
    const char * const GMS_RCS_SPHERE_YMM8R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_SPHERE_YMM8R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_SPHERE_YMM8R4_DESCRIPTION   = "AVX optimized Sphere Radar Cross Section (analytic) functionality."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

         namespace radiolocation {

                  /*
                       Radar Cross Section Handbook 1, page 147, formula 3.2-4
                       Backscattering function ,resonance region 0.4 .le. k0a .le. 20.0
                       Theta = 0, far-field
                       Valid for k0a < 1 only!!
                   */

                namespace {
                   const  static __m256 Ir  = _mm256_setzero_ps();
                   const  static __m256 Ii  = _mm256_set1_ps(1.0f);
                   const  static __m256 nIr = _mm256_set1_ps(-0.0f);
                   const  static __m256 nIi = _mm256_set1_ps(-1.0f);
                   const  static __m256 PI  = _mm256_set1_ps(3.14159265358979323846264338328f);

               }

                   __ATTR_ALWAYS_INLINE__
	           static inline
                   void bsc_func_324_ymm8r4(const __m256 k0a, // size of sphere expressed in wavenumber units
                                               __m256 * __restrict F0r, // the results
                                               __m256 * __restrict F0i) { // the results

                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        const  __m256 c0 = _mm256_set1_ps(1.5f);
                        const  __m256 c1 = _mm256_set1_ps(0.092592592592592592592592592593f);
                        const  __m256 c2 = _mm256_set1_ps(0.018888888888888888888888888889f);
                        const  __m256 c3 = _mm256_set1_ps(0.558656504577139497774418409339f);
                        const  __m256 c4 = _mm256_set1_ps(1.2f);
                        const  __m256 c5 = _mm256_set1_ps(0.5f);
                        __m256 k0a2,k0a3,k0a4,k0a6;
                        __m256 t0,t1,t2,t3;
                        k0a2 = _mm256_mul_ps(k0a,k0a);
                        t1   = _mm256_sub_ps(_1,_mm256_mul_ps(c1,k0a2));
                        k0a3 = _mm256_mul_ps(k0a2,k0a);
                        t0   = _mm256_mul_ps(c0,k0a3);
                        k0a4 = _mm256_mul_ps(k0a3,k0a);
                        k0a6 = _mm256_mul_ps(k0a3,k0a2);
                        t2   = _mm256_fmsub_ps(c2,k0a4,_mm256_mul_ps(c3,k0a6));
                        t3   = _mm256_mul_ps(t0,_mm256_add_ps(t1,t1));
                        *F0r = t3;
                        t0   = _mm256_mul_ps(c5,k0a6);
                        t1   = _mm256_add_ps(_1,_mm256_mul_ps(c4,k0a2));
                        *F0i = _mm256_mul_ps(t0,t1);
                }


                   __ATTR_ALWAYS_INLINE__
	           static inline
                   void bsc_func_324_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a, // size of sphere expressed in wavenumber units
                                               __m256 * __restrict F0r, // the results
                                               __m256 * __restrict F0i) { // the results

                         __m256 k0a = _mm256_load_ps(&pk0a[0]);
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        const  __m256 c0 = _mm256_set1_ps(1.5f);
                        const  __m256 c1 = _mm256_set1_ps(0.092592592592592592592592592593f);
                        const  __m256 c2 = _mm256_set1_ps(0.018888888888888888888888888889f);
                        const  __m256 c3 = _mm256_set1_ps(0.558656504577139497774418409339f);
                        const  __m256 c4 = _mm256_set1_ps(1.2f);
                        const  __m256 c5 = _mm256_set1_ps(0.5f);
                        __m256 k0a2,k0a3,k0a4,k0a6;
                        __m256 t0,t1,t2,t3;
                        k0a2 = _mm256_mul_ps(k0a,k0a);
                        t1   = _mm256_sub_ps(_1,_mm256_mul_ps(c1,k0a2));
                        k0a3 = _mm256_mul_ps(k0a2,k0a);
                        t0   = _mm256_mul_ps(c0,k0a3);
                        k0a4 = _mm256_mul_ps(k0a3,k0a);
                        k0a6 = _mm256_mul_ps(k0a3,k0a2);
                        t2   = _mm256_fmsub_ps(c2,k0a4,_mm256_mul_ps(c3,k0a6));
                        t3   = _mm256_mul_ps(t0,_mm256_add_ps(t1,t1));
                        *F0r = t3;
                        t0   = _mm256_mul_ps(c5,k0a6);
                        t1   = _mm256_add_ps(_1,_mm256_mul_ps(c4,k0a2));
                        *F0i = _mm256_mul_ps(t0,t1);
                }


                   __ATTR_ALWAYS_INLINE__
	           static inline
                   void bsc_func_324_ymm8r4_u(const float * __restrict  pk0a, // size of sphere expressed in wavenumber units
                                               __m256 * __restrict F0r, // the results
                                               __m256 * __restrict F0i) { // the results

                         __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        const  __m256 c0 = _mm256_set1_ps(1.5f);
                        const  __m256 c1 = _mm256_set1_ps(0.092592592592592592592592592593f);
                        const  __m256 c2 = _mm256_set1_ps(0.018888888888888888888888888889f);
                        const  __m256 c3 = _mm256_set1_ps(0.558656504577139497774418409339f);
                        const  __m256 c4 = _mm256_set1_ps(1.2f);
                        const  __m256 c5 = _mm256_set1_ps(0.5f);
                        __m256 k0a2,k0a3,k0a4,k0a6;
                        __m256 t0,t1,t2,t3;
                        k0a2 = _mm256_mul_ps(k0a,k0a);
                        t1   = _mm256_sub_ps(_1,_mm256_mul_ps(c1,k0a2));
                        k0a3 = _mm256_mul_ps(k0a2,k0a);
                        t0   = _mm256_mul_ps(c0,k0a3);
                        k0a4 = _mm256_mul_ps(k0a3,k0a);
                        k0a6 = _mm256_mul_ps(k0a3,k0a2);
                        t2   = _mm256_fmsub_ps(c2,k0a4,_mm256_mul_ps(c3,k0a6));
                        t3   = _mm256_mul_ps(t0,_mm256_add_ps(t1,t1));
                        *F0r = t3;
                        t0   = _mm256_mul_ps(c5,k0a6);
                        t1   = _mm256_add_ps(_1,_mm256_mul_ps(c4,k0a2));
                        *F0i = _mm256_mul_ps(t0,t1);
                }


                  /*
                        Radar Cross Section Handbook 1, page 147, formula 3.2-5
                        Backscattering cross section
                        
                    */
                   __ATTR_ALWAYS_INLINE__
	            static inline
                   __m256 rcs_f325_ymm8r4(const __m256 k0,
                                           const __m256 a ) {
                                        
                         __m256 a2,k0a,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,sigma; 
                        a2  = _mm256_mul_ps(a,a);
                        const  __m256 pi9 = _mm256_set1_ps(28.274333882308139146163790449516);
                        k0a = _mm256_mul_ps(k0,a); 
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        k0a2= _mm256_mul_ps(k0a,k0a);
                        const  __m256 c0 = _mm256_set1_ps(0.185185185185185185185185185185f);
                        k0a4= _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c1 = _mm256_set1_ps(0.04635116598079561042524005487f);
                        k0a6= _mm256_mul_ps(k0a4,_mm256_mul_ps(k0a,k0a));
                        const  __m256 c2 = _mm256_set1_ps(1.00969984042999916015789031662f);
                        t0  = _mm256_mul_ps(k0a4,_mm256_mul_ps(pi9,a2));
                        t1  = _mm256_sub_ps(_1,_mm256_mul_ps(c0,k0a2));
                        t2  = _mm256_fmsub_ps(c1,k0a4,_mm256_mul_ps(c2,k0a6));
                        sigma = _mm256_fmadd_ps(t1,t2,t0);
                        return (sigma);
                }


                   __ATTR_ALWAYS_INLINE__
	          static inline
                   __m256 rcs_f325_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0,
                                             const float * __restrict __ATTR_ALIGN__(32) pa ) {
                        
                         __m256 k0= _mm256_load_ps(&pk0[0]); 
                         __m256 a = _mm256_load_ps(&pa[0]);               
                         __m256 a2,k0a,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,sigma; 
                        a2  = _mm256_mul_ps(a,a);
                        const  __m256 pi9 = _mm256_set1_ps(28.274333882308139146163790449516);
                        k0a = _mm256_mul_ps(k0,a); 
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        k0a2= _mm256_mul_ps(k0a,k0a);
                        const  __m256 c0 = _mm256_set1_ps(0.185185185185185185185185185185f);
                        k0a4= _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c1 = _mm256_set1_ps(0.04635116598079561042524005487f);
                        k0a6= _mm256_mul_ps(k0a4,_mm256_mul_ps(k0a,k0a));
                        const  __m256 c2 = _mm256_set1_ps(1.00969984042999916015789031662f);
                        t0  = _mm256_mul_ps(k0a4,_mm256_mul_ps(pi9,a2));
                        t1  = _mm256_sub_ps(_1,_mm256_mul_ps(c0,k0a2));
                        t2  = _mm256_fmsub_ps(c1,k0a4,_mm256_mul_ps(c2,k0a6));
                        sigma = _mm256_fmadd_ps(t1,t2,t0);
                        return (sigma);
                }


                   __ATTR_ALWAYS_INLINE__
	         static inline
                   __m256 rcs_f325_ymm8r4_u(const float * __restrict pk0,
                                             const float * __restrict pa ) {
                        
                         __m256 k0= _mm256_loadu_ps(&pk0[0]); 
                         __m256 a = _mm256_loadu_ps(&pa[0]);               
                         __m256 a2,k0a,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,sigma; 
                        a2  = _mm256_mul_ps(a,a);
                        const  __m256 pi9 = _mm256_set1_ps(28.274333882308139146163790449516);
                        k0a = _mm256_mul_ps(k0,a); 
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        k0a2= _mm256_mul_ps(k0a,k0a);
                        const  __m256 c0 = _mm256_set1_ps(0.185185185185185185185185185185f);
                        k0a4= _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c1 = _mm256_set1_ps(0.04635116598079561042524005487f);
                        k0a6= _mm256_mul_ps(k0a4,_mm256_mul_ps(k0a,k0a));
                        const  __m256 c2 = _mm256_set1_ps(1.00969984042999916015789031662f);
                        t0  = _mm256_mul_ps(k0a4,_mm256_mul_ps(pi9,a2));
                        t1  = _mm256_sub_ps(_1,_mm256_mul_ps(c0,k0a2));
                        t2  = _mm256_fmsub_ps(c1,k0a4,_mm256_mul_ps(c2,k0a6));
                        sigma = _mm256_fmadd_ps(t1,t2,t0);
                        return (sigma);
                }



#include "GMS_complex_ymm8r4.hpp"


                  /*
                        Creeping wave term, F c(0) at the upper end of the resonance region.
                        Formula 3.2-8
                    */
                   __ATTR_ALWAYS_INLINE__
	          static inline
                   void bsc_f328_ymm8r4(const __m256 x,//k0a
                                         __m256 * __restrict Fc0r,
                                         __m256 * __restrict Fc0i) {

                       __m256 p43,pn16,pn13,pn23,p13;
                       __m256 t0r,t0i,t1r,t1i,t2r,t2i;  
                       __m256 t3r,t3i,t4r,t4i,t5r,t5i;
                       __m256 t6r,t6i,e1r,e1i,e2r,e2i;
                       __m256 e3r,e3i,t6r,t6i,t7r,t7i;
                       __m256 argr,argi,resr,resi;              
                      const  __m256 c0  = _mm256_set1_ps(1.357588f);
                      p43                       = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(1.3333333333333333333333333333333f));
                      t0r                       = _mm256_mul_ps(nIr,p43); // - / x 4/3
                      t0i                       = _mm256_mul_ps(nIi,p43); // - / x 4/3
                      const  __m256 c0r = _mm256_set1_ps(0.741196f);
                      const  __m256 c0i = _mm256_set1_ps(1.283788f);
                      pn16                      = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(0.166666666666666666666666666667f));
                      pn16                      = _mm256_rcp14_ps(pn16);
                      const  __m256 c1r = _mm256_set1_ps(2.200002f);
                      const  __m256 c1i = _mm256_set1_ps(-1.270172f);
                      pn23                      = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(0.666666666666666666666666666667f));
                      pn23                      = _mm256_rcp14_ps(pn23);
                      const  __m256 c2r = _mm256_set1_ps(0.445396f);
                      const  __m256 c2i = _mm256_set1_ps(0.257150f);
                      p13                       = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(-0.3333333333333333333333333333333333f));
                      const  __m256 c3r = _mm256_set1_ps(0.964654f);
                      const  __m256 c3i = _mm256_set1_ps(1.670829f);
                      pn13                      = _mm256_rcp14_ps(p13);
                      const  __m256 c4r = _mm256_set1_ps(7.014224f);
                      const  __m256 c4i = _mm256_set1_ps(-4.049663f);
                      t1r                       = _mm256_set1_ps(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm256_mul_ps(pn16,_mm256_set1_ps(-3.14159265358979323846264338328f)); //-//-//
                      const  __m256 c5r = _mm256_set1_ps(0.444477f);
                      const  __m256 c5i = _mm256_set1_ps(0.256619f);
                      t2r                       = _mm256_fmadd_ps(c0r,pn23,c0);
                      t2i                       = _mm256_fmadd_ps(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const  __m256 c6r = _mm256_set1_ps(0.798821f);
                      const  __m256 c6i = _mm256_set1_ps(1.383598f);
                      t3r                       = _mm256_mul_ps(pn13,c1r);
                      t3i                       = _mm256_mul_ps(pn13,c1i);
                      const  __m256 c7r = _mm256_set1_ps(5.048956f);
                      const  __m256 c7i = _mm256_set1_ps(-2.915016f);
                      t4r                       = _mm256_mul_ps(pn13,c2r);
                      t4i                       = _mm256_mul_ps(pn13,c2i);
                      const  __m256 c8r = _mm256_set1_ps(0.312321f);
                      const  __m256 c8i = _mm256_set1_ps(0.180319f);
                      t5r                       = _mm256_add_ps(t3r,t4r);
                      t5i                       = _mm256_add_ps(t3i,t4i);
                      cexp_ymm8c4(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm256_fmadd_ps(c3r,pn23,_mm256_set1_ps(0.695864f)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm256_fmadd_ps(c3i,pn23,_mm256_set1_ps(0.695864f)); //--//--//--//
                      t4r                       = _mm256_mul_ps(p13,c4r);
                      t4i                       = _mm256_mul_ps(p13,c4i);
                      t5r                       = _mm256_mul_ps(pn13,c5r);
                      t5i                       = _mm256_mul_ps(pn13,c5i);
                      t6r                       = _mm256_sub_ps(t4r,t5r);
                      t6i                       = _mm256_sub_ps(t4i,t5i);
                      cexp_ymm8c4(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm256_fmadd_ps(c4r,pn23,_mm256_set1_ps(0.807104f)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm256_fmadd_ps(c4i,pn23,_mm256_set1_ps(0.807104f)); //--//--//
                      t5r                       = _mm256_mul_ps(p13,c7r);
                      t5i                       = _mm256_mul_ps(p13,c7i);
                      t6r                       = _mm256_mul_ps(pn13,c8r);
                      t6i                       = _mm256_mul_ps(pn13,c8i);
                      t7r                       = _mm256_sub_ps(t5r,t6r);
                      t7i                       = _mm256_sub_ps(t5i,t6i);
                      cexp_ymm8c4(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm256_setzero_ps();
                      t5i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm256_setzero_ps();
                      t6i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm256_setzero_ps();
                      t7i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm256_add_ps(t5r,_mm256_sub_ps(t6r,t7r));
                      argi                      = _mm256_add_ps(t5i,_mm256_sub_ps(t6i,t7i));
                      t4r                       = _mm256_setzero_ps();
                      t4i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_ymm8c4(t4r,t4i,&resr,&resi);
                      cmul_ymm8c4(t04,t0i,resr,resi,&Fc0r,&Fc0i);                    
                }


                   __ATTR_ALWAYS_INLINE__
	          static inline
                   void bsc_f328_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) px,//k0a
                                           float * __restrict __ATTR_ALIGN__(32) Fc0r,
                                           float * __restrict __ATTR_ALIGN__(32) Fc0i) {

                       __m256 x = _mm256_load_ps(&px[0]);
                       __m256 p43,pn16,pn13,pn23,p13;
                       __m256 t0r,t0i,t1r,t1i,t2r,t2i;  
                       __m256 t3r,t3i,t4r,t4i,t5r,t5i;
                       __m256 t6r,t6i,e1r,e1i,e2r,e2i;
                       __m256 e3r,e3i,t6r,t6i,t7r,t7i;
                       __m256 argr,argi,resr,resi; 
                       __m256 fcr,fci;             
                      const  __m256 c0  = _mm256_set1_ps(1.357588f);
                      p43                       = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(1.3333333333333333333333333333333f));
                      t0r                       = _mm256_mul_ps(nIr,p43); // - / x 4/3
                      t0i                       = _mm256_mul_ps(nIi,p43); // - / x 4/3
                      const  __m256 c0r = _mm256_set1_ps(0.741196f);
                      const  __m256 c0i = _mm256_set1_ps(1.283788f);
                      pn16                      = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(0.166666666666666666666666666667f));
                      pn16                      = _mm256_rcp14_ps(pn16);
                      const  __m256 c1r = _mm256_set1_ps(2.200002f);
                      const  __m256 c1i = _mm256_set1_ps(-1.270172f);
                      pn23                      = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(0.666666666666666666666666666667f));
                      pn23                      = _mm256_rcp14_ps(pn23);
                      const  __m256 c2r = _mm256_set1_ps(0.445396f);
                      const  __m256 c2i = _mm256_set1_ps(0.257150f);
                      p13                       = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(-0.3333333333333333333333333333333333f));
                      const  __m256 c3r = _mm256_set1_ps(0.964654f);
                      const  __m256 c3i = _mm256_set1_ps(1.670829f);
                      pn13                      = _mm256_rcp14_ps(p13);
                      const  __m256 c4r = _mm256_set1_ps(7.014224f);
                      const  __m256 c4i = _mm256_set1_ps(-4.049663f);
                      t1r                       = _mm256_set1_ps(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm256_mul_ps(pn16,_mm256_set1_ps(-3.14159265358979323846264338328f)); //-//-//
                      const  __m256 c5r = _mm256_set1_ps(0.444477f);
                      const  __m256 c5i = _mm256_set1_ps(0.256619f);
                      t2r                       = _mm256_fmadd_ps(c0r,pn23,c0);
                      t2i                       = _mm256_fmadd_ps(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const  __m256 c6r = _mm256_set1_ps(0.798821f);
                      const  __m256 c6i = _mm256_set1_ps(1.383598f);
                      t3r                       = _mm256_mul_ps(pn13,c1r);
                      t3i                       = _mm256_mul_ps(pn13,c1i);
                      const  __m256 c7r = _mm256_set1_ps(5.048956f);
                      const  __m256 c7i = _mm256_set1_ps(-2.915016f);
                      t4r                       = _mm256_mul_ps(pn13,c2r);
                      t4i                       = _mm256_mul_ps(pn13,c2i);
                      const  __m256 c8r = _mm256_set1_ps(0.312321f);
                      const  __m256 c8i = _mm256_set1_ps(0.180319f);
                      t5r                       = _mm256_add_ps(t3r,t4r);
                      t5i                       = _mm256_add_ps(t3i,t4i);
                      cexp_ymm8c4(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm256_fmadd_ps(c3r,pn23,_mm256_set1_ps(0.695864f)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm256_fmadd_ps(c3i,pn23,_mm256_set1_ps(0.695864f)); //--//--//--//
                      t4r                       = _mm256_mul_ps(p13,c4r);
                      t4i                       = _mm256_mul_ps(p13,c4i);
                      t5r                       = _mm256_mul_ps(pn13,c5r);
                      t5i                       = _mm256_mul_ps(pn13,c5i);
                      t6r                       = _mm256_sub_ps(t4r,t5r);
                      t6i                       = _mm256_sub_ps(t4i,t5i);
                      cexp_ymm8c4(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm256_fmadd_ps(c4r,pn23,_mm256_set1_ps(0.807104f)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm256_fmadd_ps(c4i,pn23,_mm256_set1_ps(0.807104f)); //--//--//
                      t5r                       = _mm256_mul_ps(p13,c7r);
                      t5i                       = _mm256_mul_ps(p13,c7i);
                      t6r                       = _mm256_mul_ps(pn13,c8r);
                      t6i                       = _mm256_mul_ps(pn13,c8i);
                      t7r                       = _mm256_sub_ps(t5r,t6r);
                      t7i                       = _mm256_sub_ps(t5i,t6i);
                      cexp_ymm8c4(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm256_setzero_ps();
                      t5i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm256_setzero_ps();
                      t6i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm256_setzero_ps();
                      t7i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm256_add_ps(t5r,_mm256_sub_ps(t6r,t7r));
                      argi                      = _mm256_add_ps(t5i,_mm256_sub_ps(t6i,t7i));
                      t4r                       = _mm256_setzero_ps();
                      t4i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_ymm8c4(t4r,t4i,&resr,&resi);
                      cmul_ymm8c4(t04,t0i,resr,resi,&fcr,&fci);    
                      __m256_store_ps(&Fc0r[0], fcr);
                      __m256_store_ps(&Fc0i[0], fci);                
                }



                   __ATTR_ALWAYS_INLINE__
	           static inline
                   void bsc_f328_ymm8r4_u(const float * __restrict px,//k0a
                                           float * __restrict  Fc0r,
                                           float * __restrict  Fc0i) {

                       __m256 x = _mm256_loadu_ps(&px[0]);
                       __m256 p43,pn16,pn13,pn23,p13;
                       __m256 t0r,t0i,t1r,t1i,t2r,t2i;  
                       __m256 t3r,t3i,t4r,t4i,t5r,t5i;
                       __m256 t6r,t6i,e1r,e1i,e2r,e2i;
                       __m256 e3r,e3i,t6r,t6i,t7r,t7i;
                       __m256 argr,argi,resr,resi; 
                       __m256 fcr,fci;             
                      const  __m256 c0  = _mm256_set1_ps(1.357588f);
                      p43                       = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(1.3333333333333333333333333333333f));
                      t0r                       = _mm256_mul_ps(nIr,p43); // - / x 4/3
                      t0i                       = _mm256_mul_ps(nIi,p43); // - / x 4/3
                      const  __m256 c0r = _mm256_set1_ps(0.741196f);
                      const  __m256 c0i = _mm256_set1_ps(1.283788f);
                      pn16                      = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(0.166666666666666666666666666667f));
                      pn16                      = _mm256_rcp14_ps(pn16);
                      const  __m256 c1r = _mm256_set1_ps(2.200002f);
                      const  __m256 c1i = _mm256_set1_ps(-1.270172f);
                      pn23                      = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(0.666666666666666666666666666667f));
                      pn23                      = _mm256_rcp14_ps(pn23);
                      const  __m256 c2r = _mm256_set1_ps(0.445396f);
                      const  __m256 c2i = _mm256_set1_ps(0.257150f);
                      p13                       = _mm256_pow_ps(x,
                                                        _mm256_set1_ps(-0.3333333333333333333333333333333333f));
                      const  __m256 c3r = _mm256_set1_ps(0.964654f);
                      const  __m256 c3i = _mm256_set1_ps(1.670829f);
                      pn13                      = _mm256_rcp14_ps(p13);
                      const  __m256 c4r = _mm256_set1_ps(7.014224f);
                      const  __m256 c4i = _mm256_set1_ps(-4.049663f);
                      t1r                       = _mm256_set1_ps(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm256_mul_ps(pn16,_mm256_set1_ps(-3.14159265358979323846264338328f)); //-//-//
                      const  __m256 c5r = _mm256_set1_ps(0.444477f);
                      const  __m256 c5i = _mm256_set1_ps(0.256619f);
                      t2r                       = _mm256_fmadd_ps(c0r,pn23,c0);
                      t2i                       = _mm256_fmadd_ps(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const  __m256 c6r = _mm256_set1_ps(0.798821f);
                      const  __m256 c6i = _mm256_set1_ps(1.383598f);
                      t3r                       = _mm256_mul_ps(pn13,c1r);
                      t3i                       = _mm256_mul_ps(pn13,c1i);
                      const  __m256 c7r = _mm256_set1_ps(5.048956f);
                      const  __m256 c7i = _mm256_set1_ps(-2.915016f);
                      t4r                       = _mm256_mul_ps(pn13,c2r);
                      t4i                       = _mm256_mul_ps(pn13,c2i);
                      const  __m256 c8r = _mm256_set1_ps(0.312321f);
                      const  __m256 c8i = _mm256_set1_ps(0.180319f);
                      t5r                       = _mm256_add_ps(t3r,t4r);
                      t5i                       = _mm256_add_ps(t3i,t4i);
                      cexp_ymm8c4(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm256_fmadd_ps(c3r,pn23,_mm256_set1_ps(0.695864f)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm256_fmadd_ps(c3i,pn23,_mm256_set1_ps(0.695864f)); //--//--//--//
                      t4r                       = _mm256_mul_ps(p13,c4r);
                      t4i                       = _mm256_mul_ps(p13,c4i);
                      t5r                       = _mm256_mul_ps(pn13,c5r);
                      t5i                       = _mm256_mul_ps(pn13,c5i);
                      t6r                       = _mm256_sub_ps(t4r,t5r);
                      t6i                       = _mm256_sub_ps(t4i,t5i);
                      cexp_ymm8c4(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm256_fmadd_ps(c4r,pn23,_mm256_set1_ps(0.807104f)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm256_fmadd_ps(c4i,pn23,_mm256_set1_ps(0.807104f)); //--//--//
                      t5r                       = _mm256_mul_ps(p13,c7r);
                      t5i                       = _mm256_mul_ps(p13,c7i);
                      t6r                       = _mm256_mul_ps(pn13,c8r);
                      t6i                       = _mm256_mul_ps(pn13,c8i);
                      t7r                       = _mm256_sub_ps(t5r,t6r);
                      t7i                       = _mm256_sub_ps(t5i,t6i);
                      cexp_ymm8c4(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm256_setzero_ps();
                      t5i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm256_setzero_ps();
                      t6i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm256_setzero_ps();
                      t7i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm256_add_ps(t5r,_mm256_sub_ps(t6r,t7r));
                      argi                      = _mm256_add_ps(t5i,_mm256_sub_ps(t6i,t7i));
                      t4r                       = _mm256_setzero_ps();
                      t4i                       = _mm256_setzero_ps();
                      cmul_ymm8c4(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_ymm8c4(t4r,t4i,&resr,&resi);
                      cmul_ymm8c4(t04,t0i,resr,resi,&fcr,&fci);    
                      __m256_storeu_ps(&Fc0r[0], fcr);
                      __m256_storeu_ps(&Fc0i[0], fci);                
                }




                  



                  /*
                       The complex scattering amplitudes near the lower end of the resonance
                       region (say, 0.4 < k^a < 1) 
                       E-plane approximation
                       These equations are not valid at all for k0a > 1. They are
                       valid for all theta angles.
                   */

                   __ATTR_ALWAYS_INLINE__
	         static inline
                   void S1_f3213_zmm146r4( const __m256 k0a,
                                           const __m256 tht,
                                           __m256 * __restrict S1r,
                                           __m256 * __restrict S1i) {

                         __m256 cost,cos2t,cos3t;
                         __m256 k0a3,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,t4,t5,t6;
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                        cost                       = _mm256_cos_ps(tht);
                        const  __m256 _2   = _mm256_set1_ps(2.0f);
                        cos2t                      = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        cos3t                      = _mm256_cos_ps(_mm256_add_ps(tht,_mm256_add_ps(tht,tht)));
                        const  __m256 c0   = _mm256_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 c1   = _mm256_set1_ps(0.244444444444444444444444444444f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c2   = _mm256_set1_ps(0.083333333333333333333333333333f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c3   = _mm256_set1_ps(1.0f/120.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c4   = _mm256_set1_ps(1807.0f/70.0f);
                        t0                         = _mm256_add_ps(half,cost);  // 0.5+cost
                        const  __m256 c5   = _mm256_set1_ps(2531.0f/105.0f);
                        t1                         = _mm256_sub_ps(c0,_mm256_fmadd_ps(c1,cost,
                                                                            _mm256_mul_ps(c2,cos2t)));
                        t1                         = _mm256_mul_ps(t1,k0a2);
                        const  __m256 c6   = _mm256_set1_ps(57.0f/42.0f);
                        t2                         = _mm256_sub_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm256_fmadd_ps(c6,cos2t,
                                                                                          _mm256_mul_ps(c0,cos3t))));
                        t3                         = _mm256_mul_ps(c3,_mm256_mul_ps(t2,k0a4));
                        const  __m256 c7   = _mm256_set1_ps(0.166666666666666666666666666667f);
                        const  __m256 c8   = _mm256_set1_ps(0.2f);
                        t4                         = _mm256_mul_ps(c7,_mm256_fmsub_ps(_4,cost,_1));
                        t5                         = _mm256_mul_ps(c8,_mm256_fmadd_ps(_2,cost,_1));
                        t5                         = _mm256_mul_ps(t5,k0a2);
                        t6                         = _mm256_mul_ps(k0a6,_mm256_add_ps(t4,t5)); // imaginary part
                        *S1i                       = t6;
                        t4                         = _mm256_sub_ps(t0,_mm256_add_ps(t2,t3));
                        *S1r                       = _mm256_mul_ps(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	          static inline
                   void S1_f3213_zmm146r4_a( const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                             const float * __restrict __ATTR_ALIGN__(32) ptht,
                                             float * __restrict __ATTR_ALIGN__(32) S1r,
                                             float * __restrict __ATTR_ALIGN__(32) S1i) {

                         __m256 k0a = _mm256_load_ps(&pk0a[0]);
                         __m256 tht = _mm256_load_ps(&ptht[0]);
                         __m256 cost,cos2t,cos3t;
                         __m256 k0a3,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,t4,t5,t6;
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                        cost                       = _mm256_cos_ps(tht);
                        const  __m256 _2   = _mm256_set1_ps(2.0f);
                        cos2t                      = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        cos3t                      = _mm256_cos_ps(_mm256_add_ps(tht,_mm256_add_ps(tht,tht)));
                        const  __m256 c0   = _mm256_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 c1   = _mm256_set1_ps(0.244444444444444444444444444444f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c2   = _mm256_set1_ps(0.083333333333333333333333333333f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c3   = _mm256_set1_ps(1.0f/120.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c4   = _mm256_set1_ps(1807.0f/70.0f);
                        t0                         = _mm256_add_ps(half,cost);  // 0.5+cost
                        const  __m256 c5   = _mm256_set1_ps(2531.0f/105.0f);
                        t1                         = _mm256_sub_ps(c0,_mm256_fmadd_ps(c1,cost,
                                                                            _mm256_mul_ps(c2,cos2t)));
                        t1                         = _mm256_mul_ps(t1,k0a2);
                        const  __m256 c6   = _mm256_set1_ps(57.0f/42.0f);
                        t2                         = _mm256_sub_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm256_fmadd_ps(c6,cos2t,
                                                                                          _mm256_mul_ps(c0,cos3t))));
                        t3                         = _mm256_mul_ps(c3,_mm256_mul_ps(t2,k0a4));
                        const  __m256 c7   = _mm256_set1_ps(0.166666666666666666666666666667f);
                        const  __m256 c8   = _mm256_set1_ps(0.2f);
                        t4                         = _mm256_mul_ps(c7,_mm256_fmsub_ps(_4,cost,_1));
                        t5                         = _mm256_mul_ps(c8,_mm256_fmadd_ps(_2,cost,_1));
                        t5                         = _mm256_mul_ps(t5,k0a2);
                        t6                         = _mm256_mul_ps(k0a6,_mm256_add_ps(t4,t5)); // imaginary part
                        _mm256_store_ps(&S1i[0],t6);
                        
                        t4                         = _mm256_sub_ps(t0,_mm256_add_ps(t2,t3));
                        _mm256_store_ps(&S1r[0],_mm256_mul_ps(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	          static inline
                   void S1_f3213_zmm146r4_u( const float * __restrict  pk0a,
                                             const float * __restrict  ptht,
                                             float * __restrict S1r,
                                             float * __restrict S1i) {

                         __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                         __m256 tht = _mm256_loadu_ps(&ptht[0]);
                         __m256 cost,cos2t,cos3t;
                         __m256 k0a3,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,t4,t5,t6;
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                        cost                       = _mm256_cos_ps(tht);
                        const  __m256 _2   = _mm256_set1_ps(2.0f);
                        cos2t                      = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        cos3t                      = _mm256_cos_ps(_mm256_add_ps(tht,_mm256_add_ps(tht,tht)));
                        const  __m256 c0   = _mm256_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 c1   = _mm256_set1_ps(0.244444444444444444444444444444f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c2   = _mm256_set1_ps(0.083333333333333333333333333333f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c3   = _mm256_set1_ps(1.0f/120.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c4   = _mm256_set1_ps(1807.0f/70.0f);
                        t0                         = _mm256_add_ps(half,cost);  // 0.5+cost
                        const  __m256 c5   = _mm256_set1_ps(2531.0f/105.0f);
                        t1                         = _mm256_sub_ps(c0,_mm256_fmadd_ps(c1,cost,
                                                                            _mm256_mul_ps(c2,cos2t)));
                        t1                         = _mm256_mul_ps(t1,k0a2);
                        const  __m256 c6   = _mm256_set1_ps(57.0f/42.0f);
                        t2                         = _mm256_sub_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm256_fmadd_ps(c6,cos2t,
                                                                                          _mm256_mul_ps(c0,cos3t))));
                        t3                         = _mm256_mul_ps(c3,_mm256_mul_ps(t2,k0a4));
                        const  __m256 c7   = _mm256_set1_ps(0.166666666666666666666666666667f);
                        const  __m256 c8   = _mm256_set1_ps(0.2f);
                        t4                         = _mm256_mul_ps(c7,_mm256_fmsub_ps(_4,cost,_1));
                        t5                         = _mm256_mul_ps(c8,_mm256_fmadd_ps(_2,cost,_1));
                        t5                         = _mm256_mul_ps(t5,k0a2);
                        t6                         = _mm256_mul_ps(k0a6,_mm256_add_ps(t4,t5)); // imaginary part
                        _mm256_storeu_ps(&S1i[0],t6);
                        
                        t4                         = _mm256_sub_ps(t0,_mm256_add_ps(t2,t3));
                        _mm256_storeu_ps(&S1r[0],_mm256_mul_ps(k0a3,t4);
                }


                   /*
                       The complex scattering amplitudes near the lower end of the resonance
                       region (say, 0.4 < k^a < 1) 
                       H-plane approximation
                       These equations are not valid at all for k0a > 1. They are
                       valid for all theta angles.
                   */

                   __ATTR_ALWAYS_INLINE__
	           static inline
                   void S2_f3214_zmm146r4( const __m256 k0a,
                                           const __m256 tht,
                                           __m256 * __restrict S2r,
                                           __m256 * __restrict S2i) {

                         __m256 cost,cos2t,cos3t;
                         __m256 k0a3,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,t4,t5,t6;
                        const  __m256 c9   = _mm256_set1_ps(0.125f);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                        cost                       = _mm256_cos_ps(tht);
                        const  __m256 _2   = _mm256_set1_ps(2.0f);
                        cos2t                      = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        cos3t                      = _mm256_cos_ps(_mm256_add_ps(tht,_mm256_add_ps(tht,tht)));
                        const  __m256 c0   = _mm256_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 c1   = _mm256_set1_ps(23.0f/60.0f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c2   = _mm256_set1_ps(0.055555555555555555555555555556f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c3   = _mm256_set1_ps(1.0f/60.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c4   = _mm256_set1_ps(-1343.0f/150.0f);
                        t0                         = _mm256_fmadd_ps(half,cost,_1);  // 0.5*cost+1
                        const  __m256 c5   = _mm256_set1_ps(3769.0f/280.0f);
                        t1                         = _mm256_sub_ps(c0,_mm256_fmsub_ps(c1,cost,
                                                                            _mm256_mul_ps(c2,cos2t)));
                        t1                         = _mm256_mul_ps(t1,k0a2);
                        const  __m256 c6   = _mm256_set1_ps(57.0f/63.0f);
                        t2                         = _mm256_add_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm256_fmadd_ps(c6,cos2t,
                                                                                          _mm256_mul_ps(c9,cos3t))));
                        t3                         = _mm256_mul_ps(c3,_mm256_mul_ps(t2,k0a4));
                        const  __m256 c7   = _mm256_set1_ps(0.166666666666666666666666666667f);
                        const  __m256 c8   = _mm256_set1_ps(0.2f);
                        t4                         = _mm256_mul_ps(c7,_mm256_sub_ps(_4,cost));
                        t5                         = _mm256_mul_ps(c8,_mm256_add_ps(_2,cost));
                        t5                         = _mm256_mul_ps(t5,k0a2);
                        t6                         = _mm256_mul_ps(k0a6,_mm256_add_ps(t4,t5)); // imaginary part
                        *S2i                       = t6;
                        t4                         = _mm256_add_ps(t0,_mm256_add_ps(t2,t3));
                        *S2r                       = _mm256_mul_ps(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	          static inline
                   void S2_f3214_zmm146r4_a( const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                             const float * __restrict __ATTR_ALIGN__(32) ptht,
                                             float * __restrict __ATTR_ALIGN__(32) S2r,
                                             float * __restrict __ATTR_ALIGN__(32) S2i) {

                         __m256 k0a = _mm256_load_ps(&pk0a[0]);
                         __m256 tht = _mm256_load_ps(&ptht[0]);
                         __m256 cost,cos2t,cos3t;
                         __m256 k0a3,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,t4,t5,t6;
                        const  __m256 c9   = _mm256_set1_ps(0.125f);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                        cost                       = _mm256_cos_ps(tht);
                        const  __m256 _2   = _mm256_set1_ps(2.0f);
                        cos2t                      = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        cos3t                      = _mm256_cos_ps(_mm256_add_ps(tht,_mm256_add_ps(tht,tht)));
                        const  __m256 c0   = _mm256_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 c1   = _mm256_set1_ps(23.0f/60.0f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c2   = _mm256_set1_ps(0.055555555555555555555555555556f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c3   = _mm256_set1_ps(1.0f/60.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c4   = _mm256_set1_ps(-1343.0f/150.0f);
                        t0                         = _mm256_fmadd_ps(half,cost,_1);  // 0.5*cost+1
                        const  __m256 c5   = _mm256_set1_ps(3769.0f/280.0f);
                        t1                         = _mm256_sub_ps(c0,_mm256_fmsub_ps(c1,cost,
                                                                            _mm256_mul_ps(c2,cos2t)));
                        t1                         = _mm256_mul_ps(t1,k0a2);
                        const  __m256 c6   = _mm256_set1_ps(57.0f/63.0f);
                        t2                         = _mm256_add_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm256_fmadd_ps(c6,cos2t,
                                                                                          _mm256_mul_ps(c9,cos3t))));
                        t3                         = _mm256_mul_ps(c3,_mm256_mul_ps(t2,k0a4));
                        const  __m256 c7   = _mm256_set1_ps(0.166666666666666666666666666667f);
                        const  __m256 c8   = _mm256_set1_ps(0.2f);
                        t4                         = _mm256_mul_ps(c7,_mm256_sub_ps(_4,cost));
                        t5                         = _mm256_mul_ps(c8,_mm256_add_ps(_2,cost));
                        t5                         = _mm256_mul_ps(t5,k0a2);
                        t6                         = _mm256_mul_ps(k0a6,_mm256_add_ps(t4,t5)); // imaginary part
                        _mm256_store_ps(&S2i[0],   t6);
                        t4                         = _mm256_add_ps(t0,_mm256_add_ps(t2,t3));
                        _mm256_store_ps(&S2r[0],    _mm256_mul_ps(k0a3,t4);
                }



                   __ATTR_ALWAYS_INLINE__
	          static inline
                   void S2_f3214_zmm146r4_u( const float * __restrict  pk0a,
                                             const float * __restrict  ptht,
                                             float * __restrict  S2r,
                                             float * __restrict  S2i) {

                         __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                         __m256 tht = _mm256_loadu_ps(&ptht[0]);
                         __m256 cost,cos2t,cos3t;
                         __m256 k0a3,k0a2,k0a4,k0a6;
                         __m256 t0,t1,t2,t3,t4,t5,t6;
                        const  __m256 c9   = _mm256_set1_ps(0.125f);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                        cost                       = _mm256_cos_ps(tht);
                        const  __m256 _2   = _mm256_set1_ps(2.0f);
                        cos2t                      = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        cos3t                      = _mm256_cos_ps(_mm256_add_ps(tht,_mm256_add_ps(tht,tht)));
                        const  __m256 c0   = _mm256_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 c1   = _mm256_set1_ps(23.0f/60.0f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c2   = _mm256_set1_ps(0.055555555555555555555555555556f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c3   = _mm256_set1_ps(1.0f/60.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c4   = _mm256_set1_ps(-1343.0f/150.0f);
                        t0                         = _mm256_fmadd_ps(half,cost,_1);  // 0.5*cost+1
                        const  __m256 c5   = _mm256_set1_ps(3769.0f/280.0f);
                        t1                         = _mm256_sub_ps(c0,_mm256_fmsub_ps(c1,cost,
                                                                            _mm256_mul_ps(c2,cos2t)));
                        t1                         = _mm256_mul_ps(t1,k0a2);
                        const  __m256 c6   = _mm256_set1_ps(57.0f/63.0f);
                        t2                         = _mm256_add_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm256_fmadd_ps(c6,cos2t,
                                                                                          _mm256_mul_ps(c9,cos3t))));
                        t3                         = _mm256_mul_ps(c3,_mm256_mul_ps(t2,k0a4));
                        const  __m256 c7   = _mm256_set1_ps(0.166666666666666666666666666667f);
                        const  __m256 c8   = _mm256_set1_ps(0.2f);
                        t4                         = _mm256_mul_ps(c7,_mm256_sub_ps(_4,cost));
                        t5                         = _mm256_mul_ps(c8,_mm256_add_ps(_2,cost));
                        t5                         = _mm256_mul_ps(t5,k0a2);
                        t6                         = _mm256_mul_ps(k0a6,_mm256_add_ps(t4,t5)); // imaginary part
                        _mm256_storeu_ps(&S2i[0],   t6);
                        t4                         = _mm256_add_ps(t0,_mm256_add_ps(t2,t3));
                        _mm256_storeu_ps(&S2r[0],    _mm256_mul_ps(k0a3,t4);
                }


                  
                  /*
                       Formula 3.2-16, optics contribution at upper end of resonance region
                   */
                   __ATTR_ALWAYS_INLINE__
	           static inline
                   void S1_f3216_ymm8r4(const __m256 k0a,
                                         const __m256 tht,
                                         __m256 * __restrict S1r,
                                         __m256 * __restrict S1i) {

                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        const  __m256 _7   = _mm256_set1_ps(7.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                         __m256 k0ah,k0a2,k0aa,cexpr,cexpi;
                         __m256 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                         __m256 t0r,t0i,t1r,t1i,t2,t3,cos6t;
                        k0ah  = _mm256_mul_ps(half,k0a);
                        htht  = _mm256_mul_ps(half,tht);
                        cost  = _mm256_cos_ps(tht);
                        k0aa  = _mm256_add_ps(k0a,k0a);
                        k0a2  = _mm256_mul_ps(k0a,k0a);
                        sint  = _mm256_sin_ps(htht);
                        sin2t = _mm256_mul_ps(sint,sint);
                        carr  = _mm256_set1_ps(-0.0f);
                        cari  = _mm256_mul_ps(k0aa,htht);
                        cexp_ymm8c4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm256_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm256_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm256_mul_ps(cost,cost);
                        cos3t = _mm256_mul_ps(cost,cos2t);
                        cos6t = _mm256_mul_ps(cos3t,cos2t);
                        t3    = _mm256_div_ps(sin2t,cos6t);
                        t0r   = _mm256_setzero_ps();
                        t0i   = _mm256_div_ps(Ii,_mm256_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm256_sub_ps(_1,t0i); // second term
                        t2    = _mm256_div_ps(_7,_mm256_mul_ps(_4,k0a2));
                        cerr  = _mm256_mul_ps(t2,t3);
                        t1r   = _mm256_sub_ps(t0r,cerr);
                        t1i   = _mm256_sub_ps(t0i,cerr);
                        cmul_ymm8c4(cexpr,cexpi,t1r,t1i,&S1r,&S1i);
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3216_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(32) ptht,
                                           float * __restrict __ATTR_ALIGN__(32) S1r,
                                           float * __restrict __ATTR_ALIGN__(32) S1i) {
     
                        const  __m256 k0a  = _mm256_load_ps(&pk0a[0]);
                        const  __m256 tht  = _mm256_load_ps(&ptht[0]);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        const  __m256 _7   = _mm256_set1_ps(7.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                         __m256 k0ah,k0a2,k0aa,cexpr,cexpi;
                         __m256 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                         __m256 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm256_mul_ps(half,k0a);
                        htht  = _mm256_mul_ps(half,tht);
                        cost  = _mm256_cos_ps(tht);
                        k0aa  = _mm256_add_ps(k0a,k0a);
                        k0a2  = _mm256_mul_ps(k0a,k0a);
                        sint  = _mm256_sin_ps(htht);
                        sin2t = _mm256_mul_ps(sint,sint);
                        carr  = _mm256_set1_ps(-0.0f);
                        cari  = _mm256_mul_ps(k0aa,htht);
                        cexp_ymm8c4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm256_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm256_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm256_mul_ps(cost,cost);
                        cos3t = _mm256_mul_ps(cost,cos2t);
                        cos6t = _mm256_mul_ps(cos3t,cos2t);
                        t3    = _mm256_div_ps(sin2t,cos6t);
                        t0r   = _mm256_setzero_ps();
                        t0i   = _mm256_div_ps(Ii,_mm256_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm256_sub_ps(_1,t0i); // second term
                        t2    = _mm256_div_ps(_7,_mm256_mul_ps(_4,k0a2));
                        cerr  = _mm256_mul_ps(t2,t3);
                        t1r   = _mm256_sub_ps(t0r,cerr);
                        t1i   = _mm256_sub_ps(t0i,cerr);
                        cmul_ymm8c4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm256_store_ps(&S1r[0], resr);
                        _mm256_store_ps(&S1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3216_ymm8r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  ptht,
                                           float * __restrict  S1r,
                                           float * __restrict  S1i) {
     
                        const  __m256 k0a  = _mm256_loadu_ps(&pk0a[0]);
                        const  __m256 tht  = _mm256_loadu_ps(&ptht[0]);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 _4   = _mm256_set1_ps(4.0f);
                        const  __m256 _7   = _mm256_set1_ps(7.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                         __m256 k0ah,k0a2,k0aa,cexpr,cexpi;
                         __m256 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                         __m256 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm256_mul_ps(half,k0a);
                        htht  = _mm256_mul_ps(half,tht);
                        cost  = _mm256_cos_ps(tht);
                        k0aa  = _mm256_add_ps(k0a,k0a);
                        k0a2  = _mm256_mul_ps(k0a,k0a);
                        sint  = _mm256_sin_ps(htht);
                        sin2t = _mm256_mul_ps(sint,sint);
                        carr  = _mm256_set1_ps(-0.0f);
                        cari  = _mm256_mul_ps(k0aa,htht);
                        cexp_ymm8c4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm256_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm256_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm256_mul_ps(cost,cost);
                        cos3t = _mm256_mul_ps(cost,cos2t);
                        cos6t = _mm256_mul_ps(cos3t,cos2t);
                        t3    = _mm256_div_ps(sin2t,cos6t);
                        t0r   = _mm256_setzero_ps();
                        t0i   = _mm256_div_ps(Ii,_mm256_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm256_sub_ps(_1,t0i); // second term
                        t2    = _mm256_div_ps(_7,_mm256_mul_ps(_4,k0a2));
                        cerr  = _mm256_mul_ps(t2,t3);
                        t1r   = _mm256_sub_ps(t0r,cerr);
                        t1i   = _mm256_sub_ps(t0i,cerr);
                        cmul_ymm8c4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm256_storeu_ps(&S1r[0], resr);
                        _mm256_storeu_ps(&S1i[0], resi);
                }


                   /*
                       Formula 3.2-17, optics contribution at upper end of resonance region
                   */
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f3217_ymm8r4(const __m256 k0a,
                                         const __m256 tht,
                                         __m256 * __restrict S2r,
                                         __m256 * __restrict S2i) {

                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 quat = _mm256_set1_ps(0.25f);
                        const  __m256 _6   = _mm256_set1_ps(6.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                         __m256 k0ah,k0a2,k0aa,cexpr,cexpi;
                         __m256 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                         __m256 t0r,t0i,t1r,t1i,t2,t3,cos6t;
                        k0ah  = _mm256_mul_ps(half,k0a);
                        htht  = _mm256_mul_ps(half,tht);
                        cost  = _mm256_cos_ps(tht);
                        k0aa  = _mm256_add_ps(k0a,k0a);
                        k0a2  = _mm256_mul_ps(k0a,k0a);
                        sint  = _mm256_sin_ps(htht);
                        sin2t = _mm256_mul_ps(sint,sint);
                        carr  = _mm256_set1_ps(-0.0f);
                        cari  = _mm256_mul_ps(k0aa,htht);
                        cexp_ymm8c4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm256_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm256_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm256_mul_ps(cost,cost);
                        cos3t = _mm256_mul_ps(cost,cos2t);
                        cos6t = _mm256_mul_ps(cos3t,cos2t);
                        t3    = _mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_6,cost),sin2t),cos6t);
                        t0r   = _mm256_setzero_ps();
                        t0i   = _mm256_div_ps(_mm256_mul_ps(Ii,cost),_mm256_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm256_sub_ps(_1,t0i); // second term
                        t2    = _mm256_mul_ps(quat,k0a2);
                        cerr  = _mm256_mul_ps(t2,t3);
                        t1r   = _mm256_sub_ps(t0r,cerr);
                        t1i   = _mm256_sub_ps(t0i,cerr);
                        cmul_ymm8c4(cexpr,cexpi,t1r,t1i,&S2r,&S2i);
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f3217_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(32) ptht,
                                           float * __restrict __ATTR_ALIGN__(32) S2r,
                                           float * __restrict __ATTR_ALIGN__(32) S2i) {
                        
                         __m256 k0a        = _mm256_load_ps(&pk0a[0]);
                         __m256 tht        = _mm256_load_ps(&tht[0]);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 quat = _mm256_set1_ps(0.25f);
                        const  __m256 _6   = _mm256_set1_ps(6.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                         __m256 k0ah,k0a2,k0aa,cexpr,cexpi;
                         __m256 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                         __m256 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm256_mul_ps(half,k0a);
                        htht  = _mm256_mul_ps(half,tht);
                        cost  = _mm256_cos_ps(tht);
                        k0aa  = _mm256_add_ps(k0a,k0a);
                        k0a2  = _mm256_mul_ps(k0a,k0a);
                        sint  = _mm256_sin_ps(htht);
                        sin2t = _mm256_mul_ps(sint,sint);
                        carr  = _mm256_set1_ps(-0.0f);
                        cari  = _mm256_mul_ps(k0aa,htht);
                        cexp_ymm8c4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm256_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm256_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm256_mul_ps(cost,cost);
                        cos3t = _mm256_mul_ps(cost,cos2t);
                        cos6t = _mm256_mul_ps(cos3t,cos2t);
                        t3    = _mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_6,cost),sin2t),cos6t);
                        t0r   = _mm256_setzero_ps();
                        t0i   = _mm256_div_ps(_mm256_mul_ps(Ii,cost),_mm256_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm256_sub_ps(_1,t0i); // second term
                        t2    = _mm256_mul_ps(quat,k0a2);
                        cerr  = _mm256_mul_ps(t2,t3);
                        t1r   = _mm256_sub_ps(t0r,cerr);
                        t1i   = _mm256_sub_ps(t0i,cerr);
                        cmul_ymm8c4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm256_store_ps(&S2r[0],resr);
                        _mm256_store_ps(&S2i[0],resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f3217_ymm8r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  ptht,
                                           float * __restrict  S2r,
                                           float * __restrict  S2i) {
                        
                         __m256 k0a        = _mm256_loadu_ps(&pk0a[0]);
                         __m256 tht        = _mm256_loadu_ps(&tht[0]);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        const  __m256 quat = _mm256_set1_ps(0.25f);
                        const  __m256 _6   = _mm256_set1_ps(6.0f);
                        const  __m256 half = _mm256_set1_ps(0.5f);
                         __m256 k0ah,k0a2,k0aa,cexpr,cexpi;
                         __m256 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                         __m256 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm256_mul_ps(half,k0a);
                        htht  = _mm256_mul_ps(half,tht);
                        cost  = _mm256_cos_ps(tht);
                        k0aa  = _mm256_add_ps(k0a,k0a);
                        k0a2  = _mm256_mul_ps(k0a,k0a);
                        sint  = _mm256_sin_ps(htht);
                        sin2t = _mm256_mul_ps(sint,sint);
                        carr  = _mm256_set1_ps(-0.0f);
                        cari  = _mm256_mul_ps(k0aa,htht);
                        cexp_ymm8c4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm256_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm256_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm256_mul_ps(cost,cost);
                        cos3t = _mm256_mul_ps(cost,cos2t);
                        cos6t = _mm256_mul_ps(cos3t,cos2t);
                        t3    = _mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_6,cost),sin2t),cos6t);
                        t0r   = _mm256_setzero_ps();
                        t0i   = _mm256_div_ps(_mm256_mul_ps(Ii,cost),_mm256_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm256_sub_ps(_1,t0i); // second term
                        t2    = _mm256_mul_ps(quat,k0a2);
                        cerr  = _mm256_mul_ps(t2,t3);
                        t1r   = _mm256_sub_ps(t0r,cerr);
                        t1i   = _mm256_sub_ps(t0i,cerr);
                        cmul_ymm8c4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm256_storeu_ps(&S2r[0],resr);
                        _mm256_storeu_ps(&S2i[0],resi);
                }


                 /*
                        Low frequency region (k0a < 0.4) i.e. Rayleigh-region complex scattering
                        amplitudes.
                        Formula 3.2-20
                    */
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 S1_f3220_ymm8r4(const __m256 k0a,
                                           const __m256 tht) {

                           const __m256 half = _mm256_set1_ps(0.5f);
                           __m256 S1;
                           __m256 k0a3,cost;
                          k0a3 = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          cost = _mm256_cos_ps(tht);
                          S1   = _mm256_mul_ps(k0a3,_mm256_add_ps(half,cost));
                          return (S1);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3220_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(32) ptht,
                                           float * __restrict __ATTR_ALIGN__(32) S1) {

                           __m256 k0a = _mm256_load_ps(&pk0a[0]);
                           __m256 tht = _mm256_load_ps(&ptht[0]);
                           const __m256 half = _mm256_set1_ps(0.5f);
                           __m256 k0a3,cost;
                          k0a3 = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          cost = _mm256_cos_ps(tht);
                          _mm256_store_ps(&S1[0],_mm256_mul_ps(k0a3,_mm256_add_ps(half,cost)));
               
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3220_ymm8r4_u(const float * __restrict pk0a,
                                           const float * __restrict ptht,
                                           float * __restrict S1) {

                           __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                           __m256 tht = _mm256_loadu_ps(&ptht[0]);
                           const __m256 half = _mm256_set1_ps(0.5f);
                           __m256 k0a3,cost;
                          k0a3 = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          cost = _mm256_cos_ps(tht);
                          _mm256_storeu_ps(&S1[0],_mm256_mul_ps(k0a3,_mm256_add_ps(half,cost)));
               
               }



                   /*
                        Low frequency region (k0a < 0.4) i.e. Rayleigh-region complex scattering
                        amplitudes.
                        Formula 3.2-21
                    */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 S2_f3221_ymm8r4(const __m256 k0a,
                                           const __m256 tht) {

                           const __m256 half = _mm256_set1_ps(0.5f);
                           const __m256 _1   = _mm256_set1_ps(1.0f);
                           __m256 S2;
                           __m256 k0a3,cost,t0; 
                          k0a3  = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          cost  = _mm256_cos_ps(tht);
                          t0    = _mm256_fmadd_ps(half,cost,_1);
                          S2    = _mm256_mul_ps(k0a3,t0);
                          return (S2);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f3221_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(32) ptht,
                                           float * __restrict S2) {

                           __m256 k0a = _mm256_load_ps(&pk0a[0]);
                           __m256 tht = _mm256_load_ps(&ptht[0]);
                           const __m256 half = _mm256_set1_ps(0.5f);
                           const __m256 _1   = _mm256_set1_ps(1.0f);
                           __m256 k0a3,cost,t0; 
                          k0a3  = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          cost  = _mm256_cos_ps(tht);
                          t0    = _mm256_fmadd_ps(half,cost,_1);
                          _mm256_store_ps(&S2[0], _mm256_mul_ps(k0a3,t0));
                          
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f3221_ymm8r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  ptht,
                                           float * __restrict S2) {

                           __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                           __m256 tht = _mm256_loadu_ps(&ptht[0]);
                           const __m256 half = _mm256_set1_ps(0.5f);
                           const __m256 _1   = _mm256_set1_ps(1.0f);
                           __m256 k0a3,cost,t0; 
                          k0a3  = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          cost  = _mm256_cos_ps(tht);
                          t0    = _mm256_fmadd_ps(half,cost,_1);
                          _mm256_storeu_ps(&S2[0], _mm256_mul_ps(k0a3,t0));
                          
               }


                 /*
                         The E-plane and H-plane scattering cross sections.
                         Formula 3.2-22
                    */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3222_ymm8r4(const __m256 k0a,
                                            const __m256 a,
                                            const __m256 theta) {

                           const  __m256 _4   = _mm256_set1_ps(4.0f);
                           const  __m256 half = _mm256_set1_ps(0.5f);
                            __m256 k0a4,a2,sqr,cost,t0,t1,t2; 
                            __m256 rcs;
                           cost = _mm256_cos_ps(theta);
                           k0a4 = _mm256_mul_ps(k0a,k0a,_mm256_mul_ps(k0a,k0a));
                           a2   = _mm256_mul_ps(a,a);
                           t0   = _mm256_add_ps(half,cost);
                           sqr  = _mm256_mul_ps(t0,t0);
                           t1   = _mm256_mul_ps(pi,a2);
                           t2   = _mm256_mul_ps(_4,k0a4);
                           rcs  = _mm256_mul_ps(t1,_mm256_mul_ps(t2,sqr));
                           return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3222_ymm8r4_a(  const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(32) pa,
                                              const float * __restrict __ATTR_ALIGN__(32) ptheta,
                                              float * __restrict __ATTR_ALIGN__(32) rcs ) {

                            __m256  k0a       = _mm256_load_ps(&pk0a[0]);
                           const  __m256 _4   = _mm256_set1_ps(4.0f);
                            __m256  a         = _mm256_load_ps(&pa[0]);
                           const  __m256 half = _mm256_set1_ps(0.5f);
                            __m256 theta      = _mm256_load_ps(&ptheta[0]);
                            __m256 k0a4,a2,sqr,cost,t0,t1,t2; 
                            __m256 rcs;
                           cost = _mm256_cos_ps(theta);
                           k0a4 = _mm256_mul_ps(k0a,k0a,_mm256_mul_ps(k0a,k0a));
                           a2   = _mm256_mul_ps(a,a);
                           t0   = _mm256_add_ps(half,cost);
                           sqr  = _mm256_mul_ps(t0,t0);
                           t1   = _mm256_mul_ps(pi,a2);
                           t2   = _mm256_mul_ps(_4,k0a4);
                           _mm256_store_ps(&rcs[0], _mm256_mul_ps(t1,_mm256_mul_ps(t2,sqr)));
              }


                  __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3222_ymm8r4_u(  const float * __restrict  pk0a,
                                              const float * __restrict  pa,
                                              const float * __restrict ptheta,
                                              float * __restrict  rcs ) {

                            __m256  k0a       = _mm256_loadu_ps(&pk0a[0]);
                           const  __m256 _4   = _mm256_set1_ps(4.0f);
                            __m256  a         = _mm256_loadu_ps(&pa[0]);
                           const  __m256 half = _mm256_set1_ps(0.5f);
                            __m256 theta      = _mm256_loadu_ps(&ptheta[0]);
                            __m256 k0a4,a2,sqr,cost,t0,t1,t2; 
                            __m256 rcs;
                           cost = _mm256_cos_ps(theta);
                           k0a4 = _mm256_mul_ps(k0a,k0a,_mm256_mul_ps(k0a,k0a));
                           a2   = _mm256_mul_ps(a,a);
                           t0   = _mm256_add_ps(half,cost);
                           sqr  = _mm256_mul_ps(t0,t0);
                           t1   = _mm256_mul_ps(pi,a2);
                           t2   = _mm256_mul_ps(_4,k0a4);
                           _mm256_storeu_ps(&rcs[0], _mm256_mul_ps(t1,_mm256_mul_ps(t2,sqr)));
              }


                 /*
                         The E-plane and H-plane scattering cross sections.
                         Formula 3.2-23
                    */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3223_ymm8r4(const __m256 k0a,
                                            const __m256 a,
                                            const __m256 theta) {

                           const  __m256 _4   = _mm256_set1_ps(4.0f);
                           const  __m256 half = _mm256_set1_ps(0.5f);
                           const  __m256 _1   = _mm256_set1_ps(1.0f);
                            __m256 k0a4,a2,sqr,cost,t0,t1,t2; 
                            __m256 rcs;
                           cost = _mm256_cos_ps(theta);
                           k0a4 = _mm256_mul_ps(k0a,k0a,_mm256_mul_ps(k0a,k0a));
                           a2   = _mm256_mul_ps(a,a);
                           t0   = _mm256_fmadd_ps(half,cost,_1);
                           sqr  = _mm256_mul_ps(t0,t0);
                           t1   = _mm256_mul_ps(pi,a2);
                           t2   = _mm256_mul_ps(_4,k0a4);
                           rcs  = _mm256_mul_ps(t1,_mm256_mul_ps(t2,sqr));
                           return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3223_ymm8r4_a(  const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                              const float * __restrict __ATTR_ALIGN__(32) pa,
                                              const float * __restrict __ATTR_ALIGN__(32) ptheta,
                                              float * __restrict __ATTR_ALIGN__(32) rcs  ) {
       
                            __m256       k0a  = _mm256_load_ps(&pk0a[0]);
                           const  __m256 _4   = _mm256_set1_ps(4.0f);
                            __m256       a    = _mm256_load_ps(&pa[0]);
                           const  __m256 half = _mm256_set1_ps(0.5f);
                            __m256       theta= _mm256_load_ps(&ptheta[0]);
                           const  __m256 _1   = _mm256_set1_ps(1.0f);
                            __m256 k0a4,a2,sqr,cost,t0,t1,t2; 
                            __m256 rcs;
                           cost = _mm256_cos_ps(theta);
                           k0a4 = _mm256_mul_ps(k0a,k0a,_mm256_mul_ps(k0a,k0a));
                           a2   = _mm256_mul_ps(a,a);
                           t0   = _mm256_fmadd_ps(half,cost,_1);
                           sqr  = _mm256_mul_ps(t0,t0);
                           t1   = _mm256_mul_ps(pi,a2);
                           t2   = _mm256_mul_ps(_4,k0a4);
                           _mm256_store_ps(&rcs[0] ,_mm256_mul_ps(t1,_mm256_mul_ps(t2,sqr)));
                           
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3223_ymm8r4_u(  const float * __restrict  pk0a,
                                              const float * __restrict  pa,
                                              const float * __restrict  ptheta,
                                              float * __restrict  rcs  ) {
       
                            __m256       k0a  = _mm256_loadu_ps(&pk0a[0]);
                           const  __m256 _4   = _mm256_set1_ps(4.0f);
                            __m256       a    = _mm256_loadu_ps(&pa[0]);
                           const  __m256 half = _mm256_set1_ps(0.5f);
                            __m256       theta= _mm256_loadu_ps(&ptheta[0]);
                           const  __m256 _1   = _mm256_set1_ps(1.0f);
                            __m256 k0a4,a2,sqr,cost,t0,t1,t2; 
                            __m256 rcs;
                           cost = _mm256_cos_ps(theta);
                           k0a4 = _mm256_mul_ps(k0a,k0a,_mm256_mul_ps(k0a,k0a));
                           a2   = _mm256_mul_ps(a,a);
                           t0   = _mm256_fmadd_ps(half,cost,_1);
                           sqr  = _mm256_mul_ps(t0,t0);
                           t1   = _mm256_mul_ps(pi,a2);
                           t2   = _mm256_mul_ps(_4,k0a4);
                           _mm256_storeu_ps(&rcs[0] ,_mm256_mul_ps(t1,_mm256_mul_ps(t2,sqr)));
                           
               }


                 /*
                        High frequency region (k0a > 20).
                        Complex scattering amplitudes.
                        Formula 3.2-24
                   */


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S12_f3224_ymm8r4( const __m256 k0a,
                                           const __m256 tht,
                                           __m256 * __restrict S12r,
                                           __m256 * __restrict S12i) {

                        const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                        const  __m256 htht = _mm256_mul_ps(_mm256_set1_ps(0.5f),tht);
                         __m256 cosht,hk0a,_2k0a,carr,cari,cexr,cexi;
                        cosht = _mm256_cos_ps(htht);          
                        hk0a  = _mm256_mul_ps(nhlf,k0a);
                        _2k0a = _mm256_add_ps(k0a,k0a);
                        carr  = nIr;
                        cari  = _mm256_mul_ps(nIi,_mm256_mul_ps(_2k0a,cosht));
                        cexp_ymm8c4(carr,cari,&cexr,cexi);
                        *S12r  = _mm256_mul_ps(hk0a,cexr);
                        *S12i  = _mm256_mul_ps(hk0a,cexi);
              }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S12_f3224_ymm8r4_a( const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                             const float * __restrict __ATTR_ALIGN__(32) ptht,
                                             float * __restrict __ATTR_ALIGN__(32) S12r,
                                             float * __restrict __ATTR_ALIGN__(32) S12i) {

                        const  __m256 k0a  = _mm256_load_ps(&pk0a[0]);
                        const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                        const  __m256 tht  = _mm256_load_ps(&ptht[0]);
                        const  __m256 htht = _mm256_mul_ps(_mm256_set1_ps(0.5f),tht);
                         __m256 cosht,hk0a,_2k0a,carr,cari,cexr,cexi;
                        cosht = _mm256_cos_ps(htht);          
                        hk0a  = _mm256_mul_ps(nhlf,k0a);
                        _2k0a = _mm256_add_ps(k0a,k0a);
                        carr  = nIr;
                        cari  = _mm256_mul_ps(nIi,_mm256_mul_ps(_2k0a,cosht));
                        cexp_ymm8c4(carr,cari,&cexr,cexi);
                        _mm256_store_ps(&S12r[0], _mm256_mul_ps(hk0a,cexr));
                        _mm256_store_ps(&S12i[0], _mm256_mul_ps(hk0a,cexi));
              }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S12_f3224_ymm8r4_u( const float * __restrict  pk0a,
                                             const float * __restrict  ptht,
                                             float * __restrict  S12r,
                                             float * __restrict  S12i) {

                        const  __m256 k0a  = _mm256_loadu_ps(&pk0a[0]);
                        const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                        const  __m256 tht  = _mm256_loadu_ps(&ptht[0]);
                        const  __m256 htht = _mm256_mul_ps(_mm256_set1_ps(0.5f),tht);
                         __m256 cosht,hk0a,_2k0a,carr,cari,cexr,cexi;
                        cosht = _mm256_cos_ps(htht);          
                        hk0a  = _mm256_mul_ps(nhlf,k0a);
                        _2k0a = _mm256_add_ps(k0a,k0a);
                        carr  = nIr;
                        cari  = _mm256_mul_ps(nIi,_mm256_mul_ps(_2k0a,cosht));
                        cexp_ymm8c4(carr,cari,&cexr,cexi);
                        _mm256_storeu_ps(&S12r[0], _mm256_mul_ps(hk0a,cexr));
                        _mm256_storeu_ps(&S12i[0], _mm256_mul_ps(hk0a,cexi));
              }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3225_ymm8r4(const __m256 a) {
                         
                          __m256 rcs;
                         rcs = _mm256_mul_ps(pi,_mm256_mul_ps(a,a));
                         return (rcs);
              }


                /*
                       Resonance region (0.4 < k0a < 20.0), equations are valid only for k0a < 1.0.
                       Complex scattering amplitude represented as a scattering function -- formula 3.2-26
                 */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S12_f3226_ymm8r4(const __m256 k0a,
                                          __m256 * __restrict S12r,
                                          __m256 * __restrict S12i) {

                        __m256 k0a2,k0a3,k0a4,k0a6;
                        __m256 t0,t1,t2,t3;
                        const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c0   = _mm256_set1_ps(113.0f/90.0f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c1   = _mm256_set1_ps(1783.0f/2100.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c2   = _mm256_set1_ps(670057.0f/396900.0f);
                        t0                         = _mm256_mul_ps(nhlf,k0a3)
                        const  __m256 c3   = _mm256_set1_ps(0.833333333333333333333333333333f);
                        const  __m256 c4   = _mm256_set1_ps(0.24f);
                        t1                         = _mm256_fmsub_ps(_mm256_add_ps(_1,c0),k0a2,
                                                                 _mm256_fmsub_ps(c1,k0a4,
                                                                             _mm256_mul_ps(c2,k0a6)));
                        *S12r                      = _mm256_mul_ps(t0,t1);
                        t2                         = _mm256_mul_ps(Ini,_mm256_mul_ps(c3,k0a6));
                        t3                         = _mm256_fmadd_ps(c4,k0a3,_1);
                        *S12i                      = _mm256_mul_ps(t2,t3);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S12_f3226_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                            float * __restrict __ATTR_ALIGN__(32) S12r,
                                            float * __restrict __ATTR_ALIGN__(32) S12i) {

                        
                        __m256 k0a2,k0a3,k0a4,k0a6;
                        __m256 t0,t1,t2,t3;
                         __m256       k0a  = _mm256_load_ps(&pk0a[0]);
                        const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c0   = _mm256_set1_ps(113.0f/90.0f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c1   = _mm256_set1_ps(1783.0f/2100.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c2   = _mm256_set1_ps(670057.0f/396900.0f);
                        t0                         = _mm256_mul_ps(nhlf,k0a3)
                        const  __m256 c3   = _mm256_set1_ps(0.833333333333333333333333333333f);
                        const  __m256 c4   = _mm256_set1_ps(0.24f);
                        t1                         = _mm256_fmsub_ps(_mm256_add_ps(_1,c0),k0a2,
                                                                 _mm256_fmsub_ps(c1,k0a4,
                                                                             _mm256_mul_ps(c2,k0a6)));
                        _mm256_store_ps(&S12r[0],  _mm256_mul_ps(t0,t1));
                        t2                         = _mm256_mul_ps(Ini,_mm256_mul_ps(c3,k0a6));
                        t3                         = _mm256_fmadd_ps(c4,k0a3,_1);
                        _mm256_store_ps(&S12i[0],  _mm256_mul_ps(t2,t3));
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S12_f3226_ymm8r4_u(const float * __restrict  pk0a,
                                            float * __restrict  S12r,
                                            float * __restrict  S12i) {

                        
                        __m256 k0a2,k0a3,k0a4,k0a6;
                        __m256 t0,t1,t2,t3;
                         __m256       k0a  = _mm256_loadu_ps(&pk0a[0]);
                        const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                        k0a2                       = _mm256_mul_ps(k0a,k0a);
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        k0a3                       = _mm256_mul_ps(k0a2,k0a);
                        const  __m256 c0   = _mm256_set1_ps(113.0f/90.0f);
                        k0a4                       = _mm256_mul_ps(k0a2,k0a2);
                        const  __m256 c1   = _mm256_set1_ps(1783.0f/2100.0f);
                        k0a6                       = _mm256_mul_ps(k0a3,k0a2);
                        const  __m256 c2   = _mm256_set1_ps(670057.0f/396900.0f);
                        t0                         = _mm256_mul_ps(nhlf,k0a3)
                        const  __m256 c3   = _mm256_set1_ps(0.833333333333333333333333333333f);
                        const  __m256 c4   = _mm256_set1_ps(0.24f);
                        t1                         = _mm256_fmsub_ps(_mm256_add_ps(_1,c0),k0a2,
                                                                 _mm256_fmsub_ps(c1,k0a4,
                                                                             _mm256_mul_ps(c2,k0a6)));
                        _mm256_storeu_ps(&S12r[0],  _mm256_mul_ps(t0,t1));
                        t2                         = _mm256_mul_ps(Ini,_mm256_mul_ps(c3,k0a6));
                        t3                         = _mm256_fmadd_ps(c4,k0a3,_1);
                        _mm256_storeu_ps(&S12i[0],  _mm256_mul_ps(t2,t3));
               }


                 /*
                           Resonance region (0.4 < k0a < 20.0), equations are valid only for k0a < 1.0.
                           Radar cross-section, formula 3.2-27
                     */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3227_ymm8r4(const __m256 k0a,
                                            const __m256 a) {

                           __m256 a2,k0a2,k0a4,k0a6,t0;
                           __m256 rcs,t1;
                          a2                       = _mm256_mul_ps(a,a);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          k0a2                     = _mm256_mul_ps(k0a,k0a);
                          const  __m256 c0 = _mm256_set1_ps(113.0f/45.0f);
                          k0a4                     = _mm256_mul_ps(k0a2,k0a2);
                          const  __m256 c1 = _mm256_set1_ps(6899.0f/56700.0f);
                          k0a6                     = _mm256_mul_ps(k0a4,_mm256_mul_ps(k0a,k0a));
                          const  __m256 c2 = _mm256_set1_ps(5419129.0f/1984500.0f);
                          t0                       = _mm256_mul_ps(pi,_mm256_mul_ps(a2,k0a4));
                          t1                       = _mm256_fmsub_ps(_mm256_add_ps(_1,c0),k0a2,
                                                                             _mm256_fmsub_ps(c1,k0a4,
                                                                                         _mm256_mul_ps(c2,k0a6)));
                          rcs                      = _mm256_mul_ps(t0,t1);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3227_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(32) pa,
                                            float * __restrict __ATTR_ALIGN__(32) rcs) {

                           __m256 k0a      = _mm256_load_ps(&pk0a[0]);
                           __m256  a       = _mm256_load_ps(&pa[0]);
                           __m256 a2,k0a2,k0a4,k0a6,t0;
                           __m256 t1;
                          a2                       = _mm256_mul_ps(a,a);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          k0a2                     = _mm256_mul_ps(k0a,k0a);
                          const  __m256 c0 = _mm256_set1_ps(113.0f/45.0f);
                          k0a4                     = _mm256_mul_ps(k0a2,k0a2);
                          const  __m256 c1 = _mm256_set1_ps(6899.0f/56700.0f);
                          k0a6                     = _mm256_mul_ps(k0a4,_mm256_mul_ps(k0a,k0a));
                          const  __m256 c2 = _mm256_set1_ps(5419129.0f/1984500.0f);
                          t0                       = _mm256_mul_ps(pi,_mm256_mul_ps(a2,k0a4));
                          t1                       = _mm256_fmsub_ps(_mm256_add_ps(_1,c0),k0a2,
                                                                             _mm256_fmsub_ps(c1,k0a4,
                                                                                         _mm256_mul_ps(c2,k0a6)));
                          _mm256_store_ps(&rcs[0],  _mm256_mul_ps(t0,t1));
                         
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3227_ymm8r4_u(const float * __restrict  pk0a,
                                            const float * __restrict  pa,
                                            float * __restrict  rcs) {

                           __m256 k0a      = _mm256_loadu_ps(&pk0a[0]);
                           __m256  a       = _mm256_loadu_ps(&pa[0]);
                           __m256 a2,k0a2,k0a4,k0a6,t0;
                           __m256 t1;
                          a2                       = _mm256_mul_ps(a,a);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          k0a2                     = _mm256_mul_ps(k0a,k0a);
                          const  __m256 c0 = _mm256_set1_ps(113.0f/45.0f);
                          k0a4                     = _mm256_mul_ps(k0a2,k0a2);
                          const  __m256 c1 = _mm256_set1_ps(6899.0f/56700.0f);
                          k0a6                     = _mm256_mul_ps(k0a4,_mm256_mul_ps(k0a,k0a));
                          const  __m256 c2 = _mm256_set1_ps(5419129.0f/1984500.0f);
                          t0                       = _mm256_mul_ps(pi,_mm256_mul_ps(a2,k0a4));
                          t1                       = _mm256_fmsub_ps(_mm256_add_ps(_1,c0),k0a2,
                                                                             _mm256_fmsub_ps(c1,k0a4,
                                                                                         _mm256_mul_ps(c2,k0a6)));
                          _mm256_storeu_ps(&rcs[0],  _mm256_mul_ps(t0,t1));
                         
               }


                 /*
                       Scattering functions equation at upper end of resonance region.
                       Optics wave term and creeping wave term.
                       Optics wave term, formula 3.2-28
                   */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void FO_f3228_ymm8r4(const __m256 k0a,
                                         __m256 * __restrict FOr,
                                         __m256 * __restrict FOi) {

                        const  __m256 hlf = _mm256_set1_ps(0.5f);
                        const  __m256 c0  = _mm256_set1_ps(0.916666666666666666666666666667f);
                         __m256 k0a2,t0;
                        k0a2                      = _mm256_mul_ps(k0a,k0a);
                        t0                        = _mm256_sub_ps(k0a2,c0);
                        *FOr                      = Inr;
                        *FOi                      = _mm256_mul_ps(Ini,_mm256_mul_ps(hlf,t0));
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void FO_f3228_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                           float * __restrict __ATTR_ALIGN__(32) FOr,
                                           float * __restrict __ATTR_ALIGN__(32) FOi) {

                         __m256       k0a = _mm256_load_ps(&pk0a[0]);
                        const  __m256 hlf = _mm256_set1_ps(0.5f);
                        const  __m256 c0  = _mm256_set1_ps(0.916666666666666666666666666667f);
                         __m256 k0a2,t0;
                        k0a2                      = _mm256_mul_ps(k0a,k0a);
                        t0                        = _mm256_sub_ps(k0a2,c0);
                        _mm256_store_ps(&FOr[0],   nIr);
                        _mm256_store_ps(&FOi[0],   _mm256_mul_ps(nIi,_mm256_mul_ps(hlf,t0)));
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void FO_f3228_ymm8r4_u(const float * __restrict pk0a,
                                           float * __restrict FOr,
                                           float * __restrict FOi) {

                         __m256       k0a = _mm256_loadu_ps(&pk0a[0]);
                        const  __m256 hlf = _mm256_set1_ps(0.5f);
                        const  __m256 c0  = _mm256_set1_ps(0.916666666666666666666666666667f);
                         __m256 k0a2,t0;
                        k0a2                      = _mm256_mul_ps(k0a,k0a);
                        t0                        = _mm256_sub_ps(k0a2,c0);
                        _mm256_storeu_ps(&FOr[0],   Inr);
                        _mm256_storeu_ps(&FOi[0],   _mm256_mul_ps(Ini,_mm256_mul_ps(hlf,t0)));
                }


                    /*
                       Scattering functions equation at upper end of resonance region.
                       Optics wave term and creeping wave term.
                       Creeping wave term, formula 3.2-29
                   */
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void FC_f3229_ymm8r4(const __m256 x,
                                         __m256 * __restrict FCr,
                                         __m256 * __restrict FCi) {
                         
                         const __m256   x16        = _mm256_add_ps(x,
                                                               _mm256_set1_ps(0.666666666666666666666666666667f));
                         const __m256          _2pi= _mm256_set1_ps(6.283185307179586476925286766559f);
                         const  __m256 c0  = _mm256_set1_ps(1.357588f);
                         const __m256          px43= _mm256_pow_ps(x,
                                                               _mm256_set1_ps(1.333333333333333333333333333333f));
                         const  __m256 c1  = _mm256_set1_ps(0.807104f);
                         const __m256        pnx23 = _mm256_rcp14_ps(_mm256_pow_ps(x,
                                                               _mm256_set1_ps(0.666666666666666666666666666667f)));
                         const  __m256 c0r = _mm256_set1_ps(0.032927f);
                         const __m256        pnx43 = _mm256_rcp14_ps(px43);
                         const  __m256 c0i = _mm256_set1_ps(0.057154f);
                         const __m256        px13  = _mm256_pow_ps(x,
                                                               _mm256_set1_ps(0.333333333333333333333333333333333f));
                         const  __m256 c1r = _mm256_set1_ps(0.242679f);
                         const  __m256 c1i = _mm256_set1_ps(-0.710672f);
                         const __m256        pnx13 = _mm256_rcp14_ps(px13);
                         const  __m256 c2r = _mm256_set1_ps(0.001846f);
                         const __m256        npx13 = _mm256_sub_ps(_mm256_setzero_ps(),px13);
                         const  __m256 c2i = _mm256_set1_ps(0.008027f);
                         const  __m256 c3r = _mm256_set1_ps(0.741196f);
                         const  __m256 c3i = _mm256_set1_ps(1.283788f);
                         const  __m256 c4r = _mm256_set1_ps(4.400004f);
                         const  __m256 c4i = _mm256_set1_ps(-2.540343f);
                         const  __m256 c5r = _mm256_set1_ps(0.890792f);
                         const  __m256 c5i = _mm256_set1_ps(0.514299f);
                         const  __m256 c6r = _mm256_set1_ps(0.798821f);
                         const  __m256 c6i = _mm256_set1_ps(1.383598f);
                         const  __m256 c7r = _mm256_set1_ps(10.097912f);
                         const  __m256 c7i = _mm256_set1_ps(-5.830032f);
                         const  __m256 c8r = _mm256_set1_ps(0.624641f);
                         const  __m256 c8i = _mm256_set1_ps(0.360637f);
                         __m256 E1r,E1i,H1r,H1i,ix43r,ix43i,ear,eai,t3r,t3i;
                         __m256 t0r,t0i,t1r,t1i,t2r,t2i,cearr,ceari,tr4,t4i;
                         __m256 exp1r,exp1i,exp2r,exp2i,t5r,t5i;
                         ix43r = Inr;
                         ix43i = _mm256_mul_ps(Ini,px43);
                         ear = Ir;
                         eai = _mm256_add_ps(_mm256_mul_ps(Ii,_2pi),x16);
                         cexp_ymm8c4(ear,eai,&t3r,&t3i); 
                         cmul_ymm8c4(ix43r,ix43i,t3r,t3i,&cearr,&ceari);// ei2pi(x+1/6)
                         t0r = _mm256_fmadd_ps(pnx23,c1r,c0r);
                         t0i = _mm256_fmadd_ps(pnx23,c1i,c0i);
                         t1r = _mm256_sub_ps(t0r,_mm256_mul_ps(pnx43,c2r));
                         t1i = _mm256_sub_ps(t0i,_mm256_mul_ps(pnx43,c2i));
                         t2r = _mm256_mul_ps(ix43r,t1r); // first term re
                         t2i = _mm256_mul_ps(ix43i,t1i); // first term im
                         t0r = _mm256_fmadd_ps(xpn23,c3r,c0);
                         t0i = _mm256_fmadd_ps(xpn23,c3i,c0);
                         t4r = _mm256_fmadd_ps(npx13,c4r,
                                                    _mm256_mul_ps(pnx13,c5r));
                         t4i = _mm256_fmadd_ps(npx13,c4i,
                                                    _mm256_mul_ps(pnx13,c5i));
                         cexp_ymm8c4(t4r,t4i,&exp1r,&exp1i);
                         cmul_ymm8c4(t0r,t0i,exp1r,exp1i,&E1r,&E1i);//t3r,t3i hold second exp term
                         t1r = _mm256_fmadd_ps(pnx23,c6r,c1);
                         t1i = _mm256_fmadd_ps(pnx23,c6i,c1);
                         t5r = _mm256_fmsub_ps(npx13,c7r,
                                                     _mm256_mul_ps(pnx13,c8r));
                         t5i = _mm256_fmsub_ps(npx13,c7i,
                                                     _mm256_mul_ps(pnx13,c8i));
                         cexp_ymm8c4(t5r,t5i,&exp2r,&exp2i);
                         cmul_ymm8c4(t1r,t1i,exp2r,exp2i,&H1r,&H1i);
                         t3r = _mm256_add_ps(E1r,H1r);
                         t3i = _mm256_add_ps(E1i,H1i);
                         cmul_ymm8c4(cearr,ceari,t3r,t3i,&t1r,&t1i);
                         *FCr = _mm256_sub_ps(t2r,t1r);
                         *FCi = _mm256_sub_ps(t2i,t1i);
               } 


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void FC_f3229_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) px,
                                         float * __restrict __ATTR_ALIGN__(32) FCr,
                                         float * __restrict __ATTR_ALIGN__(32) FCi) {
                         
                          __m256 x         = _mm256_load_ps(&px[0]);
                         const __m256   x16        = _mm256_add_ps(x,
                                                               _mm256_set1_ps(0.666666666666666666666666666667f));
                         const __m256          _2pi= _mm256_set1_ps(6.283185307179586476925286766559f);
                         const  __m256 c0  = _mm256_set1_ps(1.357588f);
                         const __m256          px43= _mm256_pow_ps(x,
                                                               _mm256_set1_ps(1.333333333333333333333333333333f));
                         const  __m256 c1  = _mm256_set1_ps(0.807104f);
                         const __m256        pnx23 = _mm256_rcp14_ps(_mm256_pow_ps(x,
                                                               _mm256_set1_ps(0.666666666666666666666666666667f)));
                         const  __m256 c0r = _mm256_set1_ps(0.032927f);
                         const __m256        pnx43 = _mm256_rcp14_ps(px43);
                         const  __m256 c0i = _mm256_set1_ps(0.057154f);
                         const __m256        px13  = _mm256_pow_ps(x,
                                                               _mm256_set1_ps(0.333333333333333333333333333333333f));
                         const  __m256 c1r = _mm256_set1_ps(0.242679f);
                         const  __m256 c1i = _mm256_set1_ps(-0.710672f);
                         const __m256        pnx13 = _mm256_rcp14_ps(px13);
                         const  __m256 c2r = _mm256_set1_ps(0.001846f);
                         const __m256        npx13 = _mm256_sub_ps(_mm256_setzero_ps(),px13);
                         const  __m256 c2i = _mm256_set1_ps(0.008027f);
                         const  __m256 c3r = _mm256_set1_ps(0.741196f);
                         const  __m256 c3i = _mm256_set1_ps(1.283788f);
                         const  __m256 c4r = _mm256_set1_ps(4.400004f);
                         const  __m256 c4i = _mm256_set1_ps(-2.540343f);
                         const  __m256 c5r = _mm256_set1_ps(0.890792f);
                         const  __m256 c5i = _mm256_set1_ps(0.514299f);
                         const  __m256 c6r = _mm256_set1_ps(0.798821f);
                         const  __m256 c6i = _mm256_set1_ps(1.383598f);
                         const  __m256 c7r = _mm256_set1_ps(10.097912f);
                         const  __m256 c7i = _mm256_set1_ps(-5.830032f);
                         const  __m256 c8r = _mm256_set1_ps(0.624641f);
                         const  __m256 c8i = _mm256_set1_ps(0.360637f);
                         __m256 E1r,E1i,H1r,H1i,ix43r,ix43i,ear,eai,t3r,t3i;
                         __m256 t0r,t0i,t1r,t1i,t2r,t2i,cearr,ceari,tr4,t4i;
                         __m256 exp1r,exp1i,exp2r,exp2i,t5r,t5i;
                         ix43r = Inr;
                         ix43i = _mm256_mul_ps(Ini,px43);
                         ear = Ir;
                         eai = _mm256_add_ps(_mm256_mul_ps(Ii,_2pi),x16);
                         cexp_ymm8c4(ear,eai,&t3r,&t3i); 
                         cmul_ymm8c4(ix43r,ix43i,t3r,t3i,&cearr,&ceari);// ei2pi(x+1/6)
                         t0r = _mm256_fmadd_ps(pnx23,c1r,c0r);
                         t0i = _mm256_fmadd_ps(pnx23,c1i,c0i);
                         t1r = _mm256_sub_ps(t0r,_mm256_mul_ps(pnx43,c2r));
                         t1i = _mm256_sub_ps(t0i,_mm256_mul_ps(pnx43,c2i));
                         t2r = _mm256_mul_ps(ix43r,t1r); // first term re
                         t2i = _mm256_mul_ps(ix43i,t1i); // first term im
                         t0r = _mm256_fmadd_ps(xpn23,c3r,c0);
                         t0i = _mm256_fmadd_ps(xpn23,c3i,c0);
                         t4r = _mm256_fmadd_ps(npx13,c4r,
                                                    _mm256_mul_ps(pnx13,c5r));
                         t4i = _mm256_fmadd_ps(npx13,c4i,
                                                    _mm256_mul_ps(pnx13,c5i));
                         cexp_ymm8c4(t4r,t4i,&exp1r,&exp1i);
                         cmul_ymm8c4(t0r,t0i,exp1r,exp1i,&E1r,&E1i);//t3r,t3i hold second exp term
                         t1r = _mm256_fmadd_ps(pnx23,c6r,c1);
                         t1i = _mm256_fmadd_ps(pnx23,c6i,c1);
                         t5r = _mm256_fmsub_ps(npx13,c7r,
                                                     _mm256_mul_ps(pnx13,c8r));
                         t5i = _mm256_fmsub_ps(npx13,c7i,
                                                     _mm256_mul_ps(pnx13,c8i));
                         cexp_ymm8c4(t5r,t5i,&exp2r,&exp2i);
                         cmul_ymm8c4(t1r,t1i,exp2r,exp2i,&H1r,&H1i);
                         t3r = _mm256_add_ps(E1r,H1r);
                         t3i = _mm256_add_ps(E1i,H1i);
                         cmul_ymm8c4(cearr,ceari,t3r,t3i,&t1r,&t1i);
                         _mm256_store_ps(&FCr[0] ,_mm256_sub_ps(t2r,t1r));
                         _mm256_store_ps(&FCi[0] ,_mm256_sub_ps(t2i,t1i));
               } 

                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void FC_f3229_ymm8r4_u(const float * __restrict  px,
                                           float * __restrict  FCr,
                                           float * __restrict  FCi) {
                         
                          __m256 x         = _mm256_loadu_ps(&px[0]);
                         const __m256   x16        = _mm256_add_ps(x,
                                                               _mm256_set1_ps(0.666666666666666666666666666667f));
                         const __m256          _2pi= _mm256_set1_ps(6.283185307179586476925286766559f);
                         const  __m256 c0  = _mm256_set1_ps(1.357588f);
                         const __m256          px43= _mm256_pow_ps(x,
                                                               _mm256_set1_ps(1.333333333333333333333333333333f));
                         const  __m256 c1  = _mm256_set1_ps(0.807104f);
                         const __m256        pnx23 = _mm256_rcp14_ps(_mm256_pow_ps(x,
                                                               _mm256_set1_ps(0.666666666666666666666666666667f)));
                         const  __m256 c0r = _mm256_set1_ps(0.032927f);
                         const __m256        pnx43 = _mm256_rcp14_ps(px43);
                         const  __m256 c0i = _mm256_set1_ps(0.057154f);
                         const __m256        px13  = _mm256_pow_ps(x,
                                                               _mm256_set1_ps(0.333333333333333333333333333333333f));
                         const  __m256 c1r = _mm256_set1_ps(0.242679f);
                         const  __m256 c1i = _mm256_set1_ps(-0.710672f);
                         const __m256        pnx13 = _mm256_rcp14_ps(px13);
                         const  __m256 c2r = _mm256_set1_ps(0.001846f);
                         const __m256        npx13 = _mm256_sub_ps(_mm256_setzero_ps(),px13);
                         const  __m256 c2i = _mm256_set1_ps(0.008027f);
                         const  __m256 c3r = _mm256_set1_ps(0.741196f);
                         const  __m256 c3i = _mm256_set1_ps(1.283788f);
                         const  __m256 c4r = _mm256_set1_ps(4.400004f);
                         const  __m256 c4i = _mm256_set1_ps(-2.540343f);
                         const  __m256 c5r = _mm256_set1_ps(0.890792f);
                         const  __m256 c5i = _mm256_set1_ps(0.514299f);
                         const  __m256 c6r = _mm256_set1_ps(0.798821f);
                         const  __m256 c6i = _mm256_set1_ps(1.383598f);
                         const  __m256 c7r = _mm256_set1_ps(10.097912f);
                         const  __m256 c7i = _mm256_set1_ps(-5.830032f);
                         const  __m256 c8r = _mm256_set1_ps(0.624641f);
                         const  __m256 c8i = _mm256_set1_ps(0.360637f);
                         __m256 E1r,E1i,H1r,H1i,ix43r,ix43i,ear,eai,t3r,t3i;
                         __m256 t0r,t0i,t1r,t1i,t2r,t2i,cearr,ceari,tr4,t4i;
                         __m256 exp1r,exp1i,exp2r,exp2i,t5r,t5i;
                         ix43r = Inr;
                         ix43i = _mm256_mul_ps(Ini,px43);
                         ear = Ir;
                         eai = _mm256_add_ps(_mm256_mul_ps(Ii,_2pi),x16);
                         cexp_ymm8c4(ear,eai,&t3r,&t3i); 
                         cmul_ymm8c4(ix43r,ix43i,t3r,t3i,&cearr,&ceari);// ei2pi(x+1/6)
                         t0r = _mm256_fmadd_ps(pnx23,c1r,c0r);
                         t0i = _mm256_fmadd_ps(pnx23,c1i,c0i);
                         t1r = _mm256_sub_ps(t0r,_mm256_mul_ps(pnx43,c2r));
                         t1i = _mm256_sub_ps(t0i,_mm256_mul_ps(pnx43,c2i));
                         t2r = _mm256_mul_ps(ix43r,t1r); // first term re
                         t2i = _mm256_mul_ps(ix43i,t1i); // first term im
                         t0r = _mm256_fmadd_ps(xpn23,c3r,c0);
                         t0i = _mm256_fmadd_ps(xpn23,c3i,c0);
                         t4r = _mm256_fmadd_ps(npx13,c4r,
                                                    _mm256_mul_ps(pnx13,c5r));
                         t4i = _mm256_fmadd_ps(npx13,c4i,
                                                    _mm256_mul_ps(pnx13,c5i));
                         cexp_ymm8c4(t4r,t4i,&exp1r,&exp1i);
                         cmul_ymm8c4(t0r,t0i,exp1r,exp1i,&E1r,&E1i);//t3r,t3i hold second exp term
                         t1r = _mm256_fmadd_ps(pnx23,c6r,c1);
                         t1i = _mm256_fmadd_ps(pnx23,c6i,c1);
                         t5r = _mm256_fmsub_ps(npx13,c7r,
                                                     _mm256_mul_ps(pnx13,c8r));
                         t5i = _mm256_fmsub_ps(npx13,c7i,
                                                     _mm256_mul_ps(pnx13,c8i));
                         cexp_ymm8c4(t5r,t5i,&exp2r,&exp2i);
                         cmul_ymm8c4(t1r,t1i,exp2r,exp2i,&H1r,&H1i);
                         t3r = _mm256_add_ps(E1r,H1r);
                         t3i = _mm256_add_ps(E1i,H1i);
                         cmul_ymm8c4(cearr,ceari,t3r,t3i,&t1r,&t1i);
                         _mm256_storeu_ps(&FCr[0] ,_mm256_sub_ps(t2r,t1r));
                         _mm256_storeu_ps(&FCi[0] ,_mm256_sub_ps(t2i,t1i));
               } 


                 /*
                         Low frquency region (k0a < 0.4), Rayleigh approximation
                         for forward scattering function and cross-section.
                         Formulae: 3.2-31, 3.2-32
                         Forward scattering function.
                   */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 F_f3231_ymm8r4(const __m256 k0a) {

                          const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                           __m256 k0a3,Fpi;
                          k0a3 = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          Fpi  = _mm256_mul_ps(nhlf,k0a3);
                          return (Fpi); 
                  }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void F_f3231_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                          float * __restrict __ATTR_ALIGN__(32) Fpi) {

                          const  __m256 k0a  = _mm256_load_ps(&pk0a[0]);
                          const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                           __m256 k0a3,Fpi;
                          k0a3 = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          _mm256_store_ps(&Fpi[0] ,_mm256_mul_ps(nhlf,k0a3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void F_f3231_ymm8r4_u(const float * __restrict  pk0a,
                                          float * __restrict  Fpi) {

                          const  __m256 k0a  = _mm256_loadu_ps(&pk0a[0]);
                          const  __m256 nhlf = _mm256_set1_ps(-0.5f);
                           __m256 k0a3,Fpi;
                          k0a3 = _mm256_mul_ps(k0a,_mm256_mul_ps(k0a,k0a));
                          _mm256_storeu_ps(&Fpi[0] ,_mm256_mul_ps(nhlf,k0a3));
                 }


                    /*
                         Low frquency region (k0a < 0.4), Rayleigh approximation
                         for forward scattering function and cross-section.
                         Formulae: 3.2-31, 3.2-32
                         RCS.
                   */


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3232_ymm8r4(const __m256 k0a,
                                            const __m256 a) {

                           __m256 a2,ka04,rcs,t0;
                          a2   = _mm256_mul_ps(a,a);
                          t0   = _mm256_mul_ps(ka0,ka0);
                          ka04 = _mm256_mul_ps(t0,t0);
                          rcs  = _mm256_mul_ps(pi,_mm256_mul_ps(a2,ka04));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3232_ymm8r4_a(  const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(32) pa,
                                            float * __restrict __ATTR_ALIGN__(32) rcs) {

                           __m256 k0a = _mm256_load_ps(&pk0a[0]);
                           __m256 a   = _mm256_load_ps(&pa[0]);
                           __m256 a2,ka04,t0;
                          a2   = _mm256_mul_ps(a,a);
                          t0   = _mm256_mul_ps(ka0,ka0);
                          ka04 = _mm256_mul_ps(t0,t0);
                          _mm256_store_ps(&rcs[0] ,_mm256_mul_ps(
                                                     pi,_mm256_mul_ps(a2,ka04)));
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3232_ymm8r4_u(  const float * __restrict  pk0a,
                                              const float * __restrict    pa,
                                              float * __restrict  rcs) {

                           __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                           __m256 a   = _mm256_loadu_ps(&pa[0]);
                           __m256 a2,ka04,t0;
                          a2   = _mm256_mul_ps(a,a);
                          t0   = _mm256_mul_ps(ka0,ka0);
                          ka04 = _mm256_mul_ps(t0,t0);
                          _mm256_storeu_ps(&rcs[0] ,_mm256_mul_ps(
                                                     pi,_mm256_mul_ps(a2,ka04)));
                 }


                   /*
                         High-frequency region -- forward scattering function and 
                         cross-section (k0a > 20)
                         Formula 3.2-33 (Forward scattering function).
                     */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void F_f3233_ymm8r4(const __m256 k0a,
                                        __m256 * __restrict Fr,
                                        __m256 * __restrict Fi) {

                        const  __m256 hlf = _mm256_set1_ps(0.5f);
                         __m256 k0a2;
                        k0a2 = _mm256_mul_ps(k0a,k0a);
                        *Fr  = nIr;
                        *Fi  = _mm256_mul_ps(nIi,_mm256_mul_ps(hlf,k0a2));
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void F_f3233_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                          float * __restrict __ATTR_ALIGN__(32) Fr,
                                          float * __restrict __ATTR_ALIGN__(32) Fi) {

                         __m256 k0a       = _mm256_load_ps(&pk0a[0]);
                        const  __m256 hlf = _mm256_set1_ps(0.5f);
                         __m256 k0a2;
                        k0a2 = _mm256_mul_ps(k0a,k0a);
                        _mm256_store_ps(&Fr[0]  ,nIr);
                        _mm256_store_ps(&Fi[0]  ,_mm256_mul_ps(nIi,_mm256_mul_ps(hlf,k0a2)));
                 }


                    __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void F_f3233_ymm8r4_u(const float * __restrict  pk0a,
                                          float * __restrict  Fr,
                                          float * __restrict  Fi) {

                         __m256 k0a       = _mm256_loadu_ps(&pk0a[0]);
                        const  __m256 hlf = _mm256_set1_ps(0.5f);
                         __m256 k0a2;
                        k0a2 = _mm256_mul_ps(k0a,k0a);
                        _mm256_storeu_ps(&Fr[0]  ,nIr);
                        _mm256_storeu_ps(&Fi[0]  ,_mm256_mul_ps(nIi,_mm256_mul_ps(hlf,k0a2)));
                 }

                   /*
                         High-frequency region -- forward scattering function and 
                         cross-section (k0a > 20)
                         Formula 3.2-34 (RCS).
                     */

 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3234_ymm8r4(const __m256 k0a,
                                            const __m256 a) {

                           __m256 a2,ka02,rcs;
                          a2   = _mm256_mul_ps(a,a);
                          ka02   = _mm256_mul_ps(ka0,ka0);
                          rcs  = _mm256_mul_ps(pi,_mm256_mul_ps(a2,ka02));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3234_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                            const float * __restrict __ATTR_ALIGN__(32) pa,
                                            float * __restrict __ATTR_ALIGN__(32) rcs) {

                           __m256 k0a = _mm256_load_ps(&pk0a[0]);
                           __m256 a   = _mm256_load_ps(&pa[0]);
                           __m256 a2,ka02,rcs;
                          a2   = _mm256_mul_ps(a,a);
                          ka02   = _mm256_mul_ps(ka0,ka0);
                          _mm256_store_ps(&rcs[0], _mm256_mul_ps(pi,_mm256_mul_ps(a2,ka02)));
                          
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f3234_ymm8r4_u(const float * __restrict  pk0a,
                                            const float * __restrict  pa,
                                            float * __restrict  rcs) {

                           __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                           __m256 a   = _mm256_loadu_ps(&pa[0]);
                           __m256 a2,ka02,rcs;
                          a2   = _mm256_mul_ps(a,a);
                          ka02   = _mm256_mul_ps(ka0,ka0);
                          _mm256_storeu_ps(&rcs[0], _mm256_mul_ps(pi,_mm256_mul_ps(a2,ka02)));
                          
                 }


                  /*
                          Low-frequency region (k1a < 0.8).
                          Expansion by two series terms i.e. A0,A1 and B0,B1.
                          Formula 3.3-5
                    */
                   
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void A_coffs_f335_ymm8r4(const __m256 k0a5,
                                         const __m256 m1r,
                                         const __m256 m1i,
                                         __m256 * __restrict A1r,
                                         __m256 * __restrict A1i,
                                         __m256 * __restrict A2r,
                                         __m256 * __restrict A2i) {

                        const  __m256 c0 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                         __m256 c0r,c0i;
                        c0r = _mm256_mul_ps(_mm256_sub_ps(m1r,_1),k0a5);
                        c0i = _mm256_mul_ps(_mm256_sub_ps(m1i,_1),k0a5);
                        *A1r = _mm256_mul_ps(c0r,c0);
                        *A1i = _mm256_mul_ps(c0i,c0);
                        *A2r = Ir;
                        *A2i = Ir
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void A_coffs_f335_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pk0a5,
                                               const float * __restrict __ATTR_ALIGN__(32) pm1r,
                                               const float * __restrict __ATTR_ALIGN__(32) pm1i,
                                               float * __restrict __ATTR_ALIGN__(32) A1r,
                                               float * __restrict __ATTR_ALIGN__(32) A1i,
                                               float * __restrict __ATTR_ALIGN__(32) A2r,
                                               float * __restrict __ATTR_ALIGN__(32) A2i) {

                         __m256 k0a5     = _mm256_load_ps(&pk0a5[0]);
                         __m256 m1r      = _mm256_load_ps(&pm1r[0]); 
                         __m256 m1i      = _mm256_load_ps(&pm1i[0]);
                        const  __m256 c0 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                         __m256 c0r,c0i;
                        c0r = _mm256_mul_ps(_mm256_sub_ps(m1r,_1),k0a5);
                        c0i = _mm256_mul_ps(_mm256_sub_ps(m1i,_1),k0a5);
                        _mm256_store_ps(&A1r[0], _mm256_mul_ps(c0r,c0));
                        _mm256_store_ps(&A1i[0], _mm256_mul_ps(c0i,c0));
                        _mm256_store_ps(&A2r[0], Ir);
                        _mm256_store_ps(&A2i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void A_coffs_f335_ymm8r4_u(const float * __restrict  pk0a5,
                                               const float * __restrict  pm1r,
                                               const float * __restrict  pm1i,
                                               float * __restrict  A1r,
                                               float * __restrict  A1i,
                                               float * __restrict  A2r,
                                               float * __restrict  A2i) {

                         __m256 k0a5     = _mm256_loadu_ps(&pk0a5[0]);
                         __m256 m1r      = _mm256_loadu_ps(&pm1r[0]); 
                         __m256 m1i      = _mm256_loadu_ps(&pm1i[0]);
                        const  __m256 c0 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                         __m256 c0r,c0i;
                        c0r = _mm256_mul_ps(_mm256_sub_ps(m1r,_1),k0a5);
                        c0i = _mm256_mul_ps(_mm256_sub_ps(m1i,_1),k0a5);
                        _mm256_storeu_ps(&A1r[0], _mm256_mul_ps(c0r,c0));
                        _mm256_storeu_ps(&A1i[0], _mm256_mul_ps(c0i,c0));
                        _mm256_storeu_ps(&A2r[0], Ir);
                        _mm256_storeu_ps(&A2i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void B_coffs_f335_ymm8r4(const __m256 k0a3,
                                             const __m256 k0a5,
                                             const __m256 m1r,
                                             const __m256 m1i,
                                             __m256 * __restrict B1r,
                                             __m256 * __restrict B1i,
                                             __m256 * __restrict B2r,
                                             __m256 * __restrict B2i) {
                 
                         __m256 mm1r,mm1i,m1s1r,m1s1i,m1s2r,m1s2i;
                         __m256 m1a2r,m1a2i,_2m1r,_2m1i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mm1r                     = _mm256_mul_ps(m1r,m1r);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mm1i                     = _mm256_mul_ps(m1i,m1i);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        m1s1r                    = _mm256_sub_ps(mm1r,_1);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        m1s1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 c1 = _mm256_set1_ps(0.055555555555555555555555555556f);
                        m1s2r                    = _mm256_sub_ps(mm1r,_2);
                         __m256 t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
                        m1s2i                    = _mm256_sub_ps(mm1i,_2);
                        m1a2r                    = _mm256_add_ps(mm1r,_2);
                        m1a2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(m1a2r,m1a2i,m1a2r,m1a2i,&t3r,&t3i);
                        _2m1r                    = _mm256_fmadd_ps(_2,mm1r,_3);
                        _2m1i                    = _mm256_fmadd_ps(_2,mm1i,_3);
                        cdiv_ymm8c4(m1s1r,m1s1i,m1a2r,m1a2i,&t0r,&t0i);
                        t0r = _mm256_mul_ps(t0r,k0a3);
                        t0i = _mm256_mul_ps(t0i,k0a3);
                        cmul_ymm8c4(m1s1r,m1s1i,m1s2r,m1s2i,&t1r,&t1i);
                        cdiv_ymm8c4(t1r,t1i,t3r,t3i,&t2r,&t2i);
                        t2r = _mm256_mul_ps(c0,_mm256_mul_ps(t2r,k0a5));
                        t2i = _mm256_mul_ps(c0,_mm256_mul_ps(t2i,k0a5));
                        *B1r = _mm256_add_ps(t0r,t2r);
                        *B1i = _mm256_add_ps(t0i,t2i);
                        cdiv_ymm8c4(m1s1r,m1s1i,_2m1r,_2m1i,&t0r,&t0i);
                        *B2r = _mm256_mul_ps(c1,_mm256_mul_ps(t0r,k0a5));
                        *B2i = _mm256_mul_ps(c1,_mm256_mul_ps(t0i,k0a5));
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void B_coffs_f335_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32)  pk0a3,
                                               const float * __restrict __ATTR_ALIGN__(32)  pk0a5,
                                               const float * __restrict __ATTR_ALIGN__(32)  pm1r,
                                               const float * __restrict __ATTR_ALIGN__(32)  pm1i,
                                               float * __restrict __ATTR_ALIGN__(32) B1r,
                                               float * __restrict __ATTR_ALIGN__(32) B1i,
                                               float * __restrict __ATTR_ALIGN__(32) B2r,
                                               float * __restrict __ATTR_ALIGN__(32) B2i) {
                 
                         __m256 k0a3 = _mm256_load_ps(&pk0a3[0]);
                         __m256 k0a5 = _mm256_load_ps(&pk0a5[0]);
                         __m256 m1r  = _mm256_load_ps(&pm1r[0]);
                         __m256 m1i  = _mm256_load_ps(&pm1i[0]);
                         __m256 mm1r,mm1i,m1s1r,m1s1i,m1s2r,m1s2i;
                         __m256 m1a2r,m1a2i,_2m1r,_2m1i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mm1r                     = _mm256_mul_ps(m1r,m1r);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mm1i                     = _mm256_mul_ps(m1i,m1i);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        m1s1r                    = _mm256_sub_ps(mm1r,_1);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        m1s1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 c1 = _mm256_set1_ps(0.055555555555555555555555555556f);
                        m1s2r                    = _mm256_sub_ps(mm1r,_2);
                         __m256 t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
                        m1s2i                    = _mm256_sub_ps(mm1i,_2);
                        m1a2r                    = _mm256_add_ps(mm1r,_2);
                        m1a2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(m1a2r,m1a2i,m1a2r,m1a2i,&t3r,&t3i);
                        _2m1r                    = _mm256_fmadd_ps(_2,mm1r,_3);
                        _2m1i                    = _mm256_fmadd_ps(_2,mm1i,_3);
                        cdiv_ymm8c4(m1s1r,m1s1i,m1a2r,m1a2i,&t0r,&t0i);
                        t0r = _mm256_mul_ps(t0r,k0a3);
                        t0i = _mm256_mul_ps(t0i,k0a3);
                        cmul_ymm8c4(m1s1r,m1s1i,m1s2r,m1s2i,&t1r,&t1i);
                        cdiv_ymm8c4(t1r,t1i,t3r,t3i,&t2r,&t2i);
                        t2r = _mm256_mul_ps(c0,_mm256_mul_ps(t2r,k0a5));
                        t2i = _mm256_mul_ps(c0,_mm256_mul_ps(t2i,k0a5));
                        _mm256_store_ps(&B1r[0], _mm256_add_ps(t0r,t2r));
                        _mm256_store_ps(&B1i[0], _mm256_add_ps(t0i,t2i));
                        cdiv_ymm8c4(m1s1r,m1s1i,_2m1r,_2m1i,&t0r,&t0i);
                        _mm256_store_ps(&B2r[0], _mm256_mul_ps(c1,_mm256_mul_ps(t0r,k0a5)));
                        _mm256_store_ps(&B2i[0], _mm256_mul_ps(c1,_mm256_mul_ps(t0i,k0a5)));
                }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void B_coffs_f335_ymm8r4_u(const float * __restrict   pk0a3,
                                               const float * __restrict   pk0a5,
                                               const float * __restrict   pm1r,
                                               const float * __restrict   pm1i,
                                               float * __restrict  B1r,
                                               float * __restrict  B1i,
                                               float * __restrict  B2r,
                                               float * __restrict  B2i) {
                 
                         __m256 k0a3 = _mm256_loadu_ps(&pk0a3[0]);
                         __m256 k0a5 = _mm256_loadu_ps(&pk0a5[0]);
                         __m256 m1r  = _mm256_loadu_ps(&pm1r[0]);
                         __m256 m1i  = _mm256_loadu_ps(&pm1i[0]);
                         __m256 mm1r,mm1i,m1s1r,m1s1i,m1s2r,m1s2i;
                         __m256 m1a2r,m1a2i,_2m1r,_2m1i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mm1r                     = _mm256_mul_ps(m1r,m1r);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mm1i                     = _mm256_mul_ps(m1i,m1i);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        m1s1r                    = _mm256_sub_ps(mm1r,_1);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        m1s1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 c1 = _mm256_set1_ps(0.055555555555555555555555555556f);
                        m1s2r                    = _mm256_sub_ps(mm1r,_2);
                         __m256 t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
                        m1s2i                    = _mm256_sub_ps(mm1i,_2);
                        m1a2r                    = _mm256_add_ps(mm1r,_2);
                        m1a2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(m1a2r,m1a2i,m1a2r,m1a2i,&t3r,&t3i);
                        _2m1r                    = _mm256_fmadd_ps(_2,mm1r,_3);
                        _2m1i                    = _mm256_fmadd_ps(_2,mm1i,_3);
                        cdiv_ymm8c4(m1s1r,m1s1i,m1a2r,m1a2i,&t0r,&t0i);
                        t0r = _mm256_mul_ps(t0r,k0a3);
                        t0i = _mm256_mul_ps(t0i,k0a3);
                        cmul_ymm8c4(m1s1r,m1s1i,m1s2r,m1s2i,&t1r,&t1i);
                        cdiv_ymm8c4(t1r,t1i,t3r,t3i,&t2r,&t2i);
                        t2r = _mm256_mul_ps(c0,_mm256_mul_ps(t2r,k0a5));
                        t2i = _mm256_mul_ps(c0,_mm256_mul_ps(t2i,k0a5));
                        _mm256_storeu_ps(&B1r[0], _mm256_add_ps(t0r,t2r));
                        _mm256_storeu_ps(&B1i[0], _mm256_add_ps(t0i,t2i));
                        cdiv_ymm8c4(m1s1r,m1s1i,_2m1r,_2m1i,&t0r,&t0i);
                        _mm256_storeu_ps(&B2r[0], _mm256_mul_ps(c1,_mm256_mul_ps(t0r,k0a5)));
                        _mm256_storeu_ps(&B2i[0], _mm256_mul_ps(c1,_mm256_mul_ps(t0i,k0a5)));
                }


                   /*
                         Rayleigh backscattering RCS for dielectric spheres at angle 0.
                         Formula 3.3-7
                     */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f337_ymm8r4(const __m256 a,
                                           const __m256 k0a4,
                                           const __m256 m1r,
                                           const __m256 m1i) {

                           __m256 aa       = _mm256_mul_ps(a,a);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 _2 = _mm256_set1_ps(2.0f);
                          const  __m256 _4 = _mm256_set1_ps(4.0f);
                           __m256 fac      = _mm256_mul_ps(_4,
                                                               _mm256_mul_ps(PI,aa));
                           __m256 mm1r,mm1i,mma2r,mma2i,mms1r,mms1i,t0r,t0i;
                           __m256 cabs,rcs;
                          cmul_ymm8c4(m1r,m1i,m1r,m1i,&mm1r,&mm1i);
                          mms1r = _mm256_sub_ps(mm1r,_1);
                          mms1i = _mm256_sub_ps(mm1i,_1);
                          mma2r = _mm256_add_ps(mm1r,_2);
                          mma2i = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mm1sr,mm1si,mma2r,mma2i,&t0r,&t0i);
                          cabs =  cabs_ymm8c4(t0r,t0i);
                          rcs  = _mm256_mul_ps(fac,_mm256_mul_ps(cabs,k0a4));
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f337_ymm8r4_a(  const float * __restrict __ATTR_ALIGN__(32) pa,
                                           const float * __restrict __ATTR_ALIGN__(32) pk0a4,
                                           const float * __restrict __ATTR_ALIGN__(32) pm1r,
                                           const float * __restrict __ATTR_ALIGN__(32) pm1i,
                                           float * __restrict __ATTR_ALIGN__(32) rcs ) {

                           __m256 a        = _mm256_load_ps(&pa[0]);
                           __m256 aa       = _mm256_mul_ps(a,a);
                           __m256 k0a4     = _mm256_load_ps(&pk0a4[0]);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                           __m256 m1r      = _mm256_load_ps(&pm1r[0]);
                          const  __m256 _2 = _mm256_set1_ps(2.0f);
                           __m256 m1i      = _mm256_load_ps(&pm1i[0]);
                          const  __m256 _4 = _mm256_set1_ps(4.0f);
                           __m256 fac      = _mm256_mul_ps(_4,
                                                               _mm256_mul_ps(PI,aa));
                           __m256 mm1r,mm1i,mma2r,mma2i,mms1r,mms1i,t0r,t0i;
                           __m256 cabs,rcs;
                          cmul_ymm8c4(m1r,m1i,m1r,m1i,&mm1r,&mm1i);
                          mms1r = _mm256_sub_ps(mm1r,_1);
                          mms1i = _mm256_sub_ps(mm1i,_1);
                          mma2r = _mm256_add_ps(mm1r,_2);
                          mma2i = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mm1sr,mm1si,mma2r,mma2i,&t0r,&t0i);
                          cabs =  cabs_ymm8c4(t0r,t0i);
                          _mm256_store_ps(&rcs[0], _mm256_mul_ps(fac,_mm256_mul_ps(cabs,k0a4)));
                    
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void rcs_f337_ymm8r4_u(const float * __restrict  pa,
                                           const float * __restrict  pk0a4,
                                           const float * __restrict  pm1r,
                                           const float * __restrict  pm1i,
                                           float * __restrict  rcs ) {

                           __m256 a        = _mm256_loadu_ps(&pa[0]);
                           __m256 aa       = _mm256_mul_ps(a,a);
                           __m256 k0a4     = _mm256_loadu_ps(&pk0a4[0]);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                           __m256 m1r      = _mm256_loadu_ps(&pm1r[0]);
                          const  __m256 _2 = _mm256_set1_ps(2.0f);
                           __m256 m1i      = _mm256_loadu_ps(&pm1i[0]);
                          const  __m256 _4 = _mm256_set1_ps(4.0f);
                           __m256 fac      = _mm256_mul_ps(_4,
                                                               _mm256_mul_ps(PI,aa));
                           __m256 mm1r,mm1i,mma2r,mma2i,mms1r,mms1i,t0r,t0i;
                           __m256 cabs,rcs;
                          cmul_ymm8c4(m1r,m1i,m1r,m1i,&mm1r,&mm1i);
                          mms1r = _mm256_sub_ps(mm1r,_1);
                          mms1i = _mm256_sub_ps(mm1i,_1);
                          mma2r = _mm256_add_ps(mm1r,_2);
                          mma2i = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mm1sr,mm1si,mma2r,mma2i,&t0r,&t0i);
                          cabs =  cabs_ymm8c4(t0r,t0i);
                          _mm256_storeu_ps(&rcs[0], _mm256_mul_ps(fac,_mm256_mul_ps(cabs,k0a4)));
                    
               }


                 /*
                        Low-frequency bi-static scattering
                  */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f338_ymm8r4(const __m256 ka03,
                                        const __m256 ka05,
                                        const __m256 tht,
                                        const __m256 mm1r,
                                        const __m256 mm1i,
                                        __m256 * __restrict S1r,
                                        __m256 * __restrict S1i) {

                         __m256 mms1r,mms1i,mms2r,mms2i;
                         __m256 mma2r,mma2i,mma3r,mma3i;
                         __m256 div1r,div1i,div2r,div2i;
                         __m256 div3r,div3i,mulr,muli,t1r,t1i;
                         __m256 cost,cos2t,msqr,msqi,t0r,t0i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mms1r                    = _mm256_sub_ps(mm1r,_1);
                        mms1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mms2r                    = _mm256_sub_ps(mm1r,_2);
                        mms2i                    = _mm256_sub_ps(mm1i,_2);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        mma2r                    = _mm256_add_ps(mm1r,_2);
                        mma2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        mma3r                    = _mm256_add_ps(mm1r,_3);
                        mma3i                    = _mm256_add_ps(mm1i,_3);
                        const  __m256 c1 = _mm256_set1_ps(0.166666666666666666666666666667f);
                        cost                     = _mm256_cos_ps(tht);
                        const  __m256 c2 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        cos2t                    = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm256_mul_ps(div1r,_mm256_mul_ps(cost,ka03));
                        cmul_ymm8c4(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm256_mul_ps(div1i,_mm256_mul_ps(cost,ka03));
                        cdiv_ymm8c4(mulr,muli,msqr,msqi,&div2r,&div2i);
                        div2r = _mm256_mul_ps(div2r,cost);
                        div2i = _mm256_mul_ps(div2i,cost);
                        t0r   = _mm256_mul_ps(c2,mms1r);
                        t0i   = _mm256_mul_ps(c2,mms1i);
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm256_mul_ps(c1,_mm256_mul_ps(div3r,cos2t));
                        div3i = _mm256_mul_ps(c1,_mm256_mul_ps(div3i,cos2t));
                        t1r   = _mm256_sub_ps(div2r,_mm256_sub_ps(div3r,t0r));
                        t1i   = _mm256_sub_ps(div2i,_mm256_sub_ps(div3i,t0i));
                        *S1r  = _mm256_fmadd_ps(t1r,k0a5,div1r);
                        *S1i  = _mm256_fmadd_ps(t1i,k0a5,div1i);
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f338_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pka03,
                                        const float * __restrict __ATTR_ALIGN__(32) pka05,
                                        const float * __restrict __ATTR_ALIGN__(32) ptht,
                                        const float * __restrict __ATTR_ALIGN__(32) pmm1r,
                                        const float * __restrict __ATTR_ALIGN__(32) pmm1i,
                                        float * __restrict __ATTR_ALIGN__(32) S1r,
                                        float * __restrict __ATTR_ALIGN__(32) S1i) {

                         __m256 ka03     = _mm256_load_ps(&pka03[0]);
                         __m256 ka05     = _mm256_load_ps(&pka05[0]);
                         __m256 tht      = _mm256_load_ps(&ptht[0]);
                         __m256 mm1r     = _mm256_load_ps(&pmm1r[0]);
                         __m256 mm1i     = _mm256_load_ps(&pmm1i[0]);
                         __m256 mms1r,mms1i,mms2r,mms2i;
                         __m256 mma2r,mma2i,mma3r,mma3i;
                         __m256 div1r,div1i,div2r,div2i;
                         __m256 div3r,div3i,mulr,muli,t1r,t1i;
                         __m256 cost,cos2t,msqr,msqi,t0r,t0i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mms1r                    = _mm256_sub_ps(mm1r,_1);
                        mms1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mms2r                    = _mm256_sub_ps(mm1r,_2);
                        mms2i                    = _mm256_sub_ps(mm1i,_2);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        mma2r                    = _mm256_add_ps(mm1r,_2);
                        mma2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        mma3r                    = _mm256_add_ps(mm1r,_3);
                        mma3i                    = _mm256_add_ps(mm1i,_3);
                        const  __m256 c1 = _mm256_set1_ps(0.166666666666666666666666666667f);
                        cost                     = _mm256_cos_ps(tht);
                        const  __m256 c2 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        cos2t                    = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm256_mul_ps(div1r,_mm256_mul_ps(cost,ka03));
                        cmul_ymm8c4(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm256_mul_ps(div1i,_mm256_mul_ps(cost,ka03));
                        cdiv_ymm8c4(mulr,muli,msqr,msqi,&div2r,&div2i);
                        div2r = _mm256_mul_ps(div2r,cost);
                        div2i = _mm256_mul_ps(div2i,cost);
                        t0r   = _mm256_mul_ps(c2,mms1r);
                        t0i   = _mm256_mul_ps(c2,mms1i);
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm256_mul_ps(c1,_mm256_mul_ps(div3r,cos2t));
                        div3i = _mm256_mul_ps(c1,_mm256_mul_ps(div3i,cos2t));
                        t1r   = _mm256_sub_ps(div2r,_mm256_sub_ps(div3r,t0r));
                        t1i   = _mm256_sub_ps(div2i,_mm256_sub_ps(div3i,t0i));
                        _mm256_store_ps(&S1r[0] ,_mm256_fmadd_ps(t1r,k0a5,div1r));
                        _mm256_store_ps(&S1i[0] ,_mm256_fmadd_ps(t1i,k0a5,div1i));
               }


                  __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f338_ymm8r4_u(const float * __restrict  pka03,
                                        const float * __restrict  pka05,
                                        const float * __restrict  ptht,
                                        const float * __restrict  pmm1r,
                                        const float * __restrict  pmm1i,
                                        float * __restrict  S1r,
                                        float * __restrict  S1i) {

                         __m256 ka03     = _mm256_loadu_ps(&pka03[0]);
                         __m256 ka05     = _mm256_loadu_ps(&pka05[0]);
                         __m256 tht      = _mm256_loadu_ps(&ptht[0]);
                         __m256 mm1r     = _mm256_loadu_ps(&pmm1r[0]);
                         __m256 mm1i     = _mm256_loadu_ps(&pmm1i[0]);
                         __m256 mms1r,mms1i,mms2r,mms2i;
                         __m256 mma2r,mma2i,mma3r,mma3i;
                         __m256 div1r,div1i,div2r,div2i;
                         __m256 div3r,div3i,mulr,muli,t1r,t1i;
                         __m256 cost,cos2t,msqr,msqi,t0r,t0i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mms1r                    = _mm256_sub_ps(mm1r,_1);
                        mms1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mms2r                    = _mm256_sub_ps(mm1r,_2);
                        mms2i                    = _mm256_sub_ps(mm1i,_2);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        mma2r                    = _mm256_add_ps(mm1r,_2);
                        mma2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        mma3r                    = _mm256_add_ps(mm1r,_3);
                        mma3i                    = _mm256_add_ps(mm1i,_3);
                        const  __m256 c1 = _mm256_set1_ps(0.166666666666666666666666666667f);
                        cost                     = _mm256_cos_ps(tht);
                        const  __m256 c2 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        cos2t                    = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm256_mul_ps(div1r,_mm256_mul_ps(cost,ka03));
                        cmul_ymm8c4(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm256_mul_ps(div1i,_mm256_mul_ps(cost,ka03));
                        cdiv_ymm8c4(mulr,muli,msqr,msqi,&div2r,&div2i);
                        div2r = _mm256_mul_ps(div2r,cost);
                        div2i = _mm256_mul_ps(div2i,cost);
                        t0r   = _mm256_mul_ps(c2,mms1r);
                        t0i   = _mm256_mul_ps(c2,mms1i);
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm256_mul_ps(c1,_mm256_mul_ps(div3r,cos2t));
                        div3i = _mm256_mul_ps(c1,_mm256_mul_ps(div3i,cos2t));
                        t1r   = _mm256_sub_ps(div2r,_mm256_sub_ps(div3r,t0r));
                        t1i   = _mm256_sub_ps(div2i,_mm256_sub_ps(div3i,t0i));
                        _mm256_storeu_ps(&S1r[0] ,_mm256_fmadd_ps(t1r,k0a5,div1r));
                        _mm256_storeu_ps(&S1i[0] ,_mm256_fmadd_ps(t1i,k0a5,div1i));
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f338_ymm8r4(const __m256 ka03,
                                        const __m256 ka05,
                                        const __m256 tht,
                                        const __m256 mm1r,
                                        const __m256 mm1i,
                                        __m256 * __restrict S2r,
                                        __m256 * __restrict S2i) {

                         __m256 mms1r,mms1i,mms2r,mms2i;
                         __m256 mma2r,mma2i,mma3r,mma3i;
                         __m256 div1r,div1i,div2r,div2i;
                         __m256 div3r,div3i,mulr,muli,t1r,t1i;
                         __m256 cost,cos2t,msqr,msqi,t0r,t0i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mms1r                    = _mm256_sub_ps(mm1r,_1);
                        mms1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mms2r                    = _mm256_sub_ps(mm1r,_2);
                        mms2i                    = _mm256_sub_ps(mm1i,_2);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        mma2r                    = _mm256_add_ps(mm1r,_2);
                        mma2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        mma3r                    = _mm256_add_ps(mm1r,_3);
                        mma3i                    = _mm256_add_ps(mm1i,_3);
                        const  __m256 c1 = _mm256_set1_ps(0.166666666666666666666666666667f);
                        cost                     = _mm256_cos_ps(tht);
                        const  __m256 c2 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        cos2t                    = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm256_mul_ps(div1r,ka03);
                        cmul_ymm8c4(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm256_mul_ps(div1i,ka03);
                        cdiv_ymm8c4(mulr,muli,msqr,msqi,&div2r,&div2i);
                        t0r   = _mm256_mul_ps(c2,mms1r);
                        t0i   = _mm256_mul_ps(c2,mms1i);
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm256_mul_ps(c1,_mm256_mul_ps(div3r,cost));
                        div3i = _mm256_mul_ps(c1,_mm256_mul_ps(div2i,cost));
                        t1r   = _mm256_sub_ps(div2r,_mm256_sub_ps(div3r,t0r));
                        t1i   = _mm256_sub_ps(div2i,_mm256_sub_ps(div3i,t0i));
                        *S2r  = _mm256_fmadd_ps(t1r,cost,div1r);
                        *S2i  = _mm256_fmadd_ps(t1i,cost,div1i);
               }


                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f338_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pka03,
                                          const float * __restrict __ATTR_ALIGN__(32) pka05,
                                          const float * __restrict __ATTR_ALIGN__(32) ptht,
                                          const float * __restrict __ATTR_ALIGN__(32) pmm1r,
                                          const float * __restrict __ATTR_ALIGN__(32) pmm1i,
                                          float * __restrict __ATTR_ALIGN__(32) S2r,
                                          float * __restrict __ATTR_ALIGN__(32) S2i) {

                         __m256 ka03     = _mm256_load_ps(&pka03[0]);
                         __m256 ka05     = _mm256_load_ps(&pka05[0]);
                         __m256 tht      = _mm256_load_ps(&ptht[0]);
                         __m256 mm1r     = _mm256_load_ps(&pmm1r[0]);
                         __m256 mm1i     = _mm256_load_ps(&pmm1i[0]);
                         __m256 mms1r,mms1i,mms2r,mms2i;
                         __m256 mma2r,mma2i,mma3r,mma3i;
                         __m256 div1r,div1i,div2r,div2i;
                         __m256 div3r,div3i,mulr,muli,t1r,t1i;
                         __m256 cost,cos2t,msqr,msqi,t0r,t0i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mms1r                    = _mm256_sub_ps(mm1r,_1);
                        mms1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mms2r                    = _mm256_sub_ps(mm1r,_2);
                        mms2i                    = _mm256_sub_ps(mm1i,_2);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        mma2r                    = _mm256_add_ps(mm1r,_2);
                        mma2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        mma3r                    = _mm256_add_ps(mm1r,_3);
                        mma3i                    = _mm256_add_ps(mm1i,_3);
                        const  __m256 c1 = _mm256_set1_ps(0.166666666666666666666666666667f);
                        cost                     = _mm256_cos_ps(tht);
                        const  __m256 c2 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        cos2t                    = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm256_mul_ps(div1r,ka03);
                        cmul_ymm8c4(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm256_mul_ps(div1i,ka03);
                        cdiv_ymm8c4(mulr,muli,msqr,msqi,&div2r,&div2i);
                        t0r   = _mm256_mul_ps(c2,mms1r);
                        t0i   = _mm256_mul_ps(c2,mms1i);
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm256_mul_ps(c1,_mm256_mul_ps(div3r,cost));
                        div3i = _mm256_mul_ps(c1,_mm256_mul_ps(div2i,cost));
                        t1r   = _mm256_sub_ps(div2r,_mm256_sub_ps(div3r,t0r));
                        t1i   = _mm256_sub_ps(div2i,_mm256_sub_ps(div3i,t0i));
                        _mm256_store_ps(&S2r[0] ,_mm256_fmadd_ps(t1r,cost,div1r));
                        _mm256_store_ps(&S2i[0] ,_mm256_fmadd_ps(t1i,cost,div1i));
               }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S2_f338_ymm8r4_u(const float * __restrict  pka03,
                                          const float * __restrict  pka05,
                                          const float * __restrict  ptht,
                                          const float * __restrict  pmm1r,
                                          const float * __restrict  pmm1i,
                                          float * __restrict  S2r,
                                          float * __restrict  S2i) {

                         __m256 ka03     = _mm256_loadu_ps(&pka03[0]);
                         __m256 ka05     = _mm256_loadu_ps(&pka05[0]);
                         __m256 tht      = _mm256_loadu_ps(&ptht[0]);
                         __m256 mm1r     = _mm256_loadu_ps(&pmm1r[0]);
                         __m256 mm1i     = _mm256_loadu_ps(&pmm1i[0]);
                         __m256 mms1r,mms1i,mms2r,mms2i;
                         __m256 mma2r,mma2i,mma3r,mma3i;
                         __m256 div1r,div1i,div2r,div2i;
                         __m256 div3r,div3i,mulr,muli,t1r,t1i;
                         __m256 cost,cos2t,msqr,msqi,t0r,t0i;
                        const  __m256 _1 = _mm256_set1_ps(1.0f);
                        mms1r                    = _mm256_sub_ps(mm1r,_1);
                        mms1i                    = _mm256_sub_ps(mm1i,_1);
                        const  __m256 _2 = _mm256_set1_ps(2.0f);
                        mms2r                    = _mm256_sub_ps(mm1r,_2);
                        mms2i                    = _mm256_sub_ps(mm1i,_2);
                        const  __m256 _3 = _mm256_set1_ps(3.0f);
                        mma2r                    = _mm256_add_ps(mm1r,_2);
                        mma2i                    = _mm256_add_ps(mm1i,_2);
                        cmul_ymm8c4(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const  __m256 c0 = _mm256_set1_ps(0.6f);
                        mma3r                    = _mm256_add_ps(mm1r,_3);
                        mma3i                    = _mm256_add_ps(mm1i,_3);
                        const  __m256 c1 = _mm256_set1_ps(0.166666666666666666666666666667f);
                        cost                     = _mm256_cos_ps(tht);
                        const  __m256 c2 = _mm256_set1_ps(0.033333333333333333333333333333f);
                        cos2t                    = _mm256_cos_ps(_mm256_add_ps(tht,tht));
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm256_mul_ps(div1r,ka03);
                        cmul_ymm8c4(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm256_mul_ps(div1i,ka03);
                        cdiv_ymm8c4(mulr,muli,msqr,msqi,&div2r,&div2i);
                        t0r   = _mm256_mul_ps(c2,mms1r);
                        t0i   = _mm256_mul_ps(c2,mms1i);
                        cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm256_mul_ps(c1,_mm256_mul_ps(div3r,cost));
                        div3i = _mm256_mul_ps(c1,_mm256_mul_ps(div2i,cost));
                        t1r   = _mm256_sub_ps(div2r,_mm256_sub_ps(div3r,t0r));
                        t1i   = _mm256_sub_ps(div2i,_mm256_sub_ps(div3i,t0i));
                        _mm256_storeu_ps(&S2r[0] ,_mm256_fmadd_ps(t1r,cost,div1r));
                        _mm256_storeu_ps(&S2i[0] ,_mm256_fmadd_ps(t1i,cost,div1i));
               }

                   /*
                         E-plane and H-plane RCS.
                         Formulae: 3.3-10,3.3-11
                     */
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3310_ymm8r4(const __m256 tht,
                                            const __m256 a,
                                            const __m256 ka04,
                                            const __m256 mm1r,
                                            const __m256 mm1i) {

                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 aa = _mm256_mul_ps(a,a);
                          const  __m256 _2 = _mm256_set1_ps(2.0f); 
                          const  __m256 cost = _mm256_cos_ps(tht);
                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 cos2t,mms1r,mms1i,mma2r,mma2i;
                           __m256 divr,divi,rcs,frac,cabs;
                          cos2t   = _mm256_mul_ps(cost,cost);
                          mms1r   = _mm256_sub_ps(mm1r,_1);
                          mma2r   = _mm256_add_ps(mm1r,_2);
                          mms1i   = _mm256_sub_ps(mm1i,_1);
                          mma2i   = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm256_mul_ps(_4pi,aa);
                          cabs    =  cabs_ymm8c4(divr,divi);
                          rcs     = _mm256_mul_ps(_mm256_mul_ps(frac,cabs),
                                                  _mm256_mul_ps(cos2t,ka04));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3310_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) ptht,
                                            const float * __restrict __ATTR_ALIGN__(32) pa,
                                            const float * __restrict __ATTR_ALIGN__(32) pka04,
                                            const float * __restrict __ATTR_ALIGN__(32) pmm1r,
                                            const float * __restrict __ATTR_ALIGN__(32) pmm1i) {

                          const  __m256 tht = _mm256_load_ps(&ptht[0]);
                          const  __m256 a   = _mm256_load_ps(&pa[0]);
                          const  __m256 ka04= _mm256_load_ps(&pka04[0]);
                          const  __m256 mm1r= _mm256_load_ps(&pmm1r[0]);
                          const  __m256 mm1i= _mm256_load_ps(&pmm1i[0]);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 aa = _mm256_mul_ps(a,a);
                          const  __m256 _2 = _mm256_set1_ps(2.0f); 
                          const  __m256 cost = _mm256_cos_ps(tht);
                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 cos2t,mms1r,mms1i,mma2r,mma2i;
                           __m256 divr,divi,rcs,frac,cabs;
                          cos2t   = _mm256_mul_ps(cost,cost);
                          mms1r   = _mm256_sub_ps(mm1r,_1);
                          mma2r   = _mm256_add_ps(mm1r,_2);
                          mms1i   = _mm256_sub_ps(mm1i,_1);
                          mma2i   = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm256_mul_ps(_4pi,aa);
                          cabs    =  cabs_ymm8c4(divr,divi);
                          rcs     = _mm256_mul_ps(_mm256_mul_ps(frac,cabs),
                                                  _mm256_mul_ps(cos2t,ka04));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3310_ymm8r4_u(const float * __restrict  ptht,
                                              const float * __restrict  pa,
                                              const float * __restrict  pka04,
                                              const float * __restrict  pmm1r,
                                              const float * __restrict  pmm1i) {

                          const  __m256 tht = _mm256_loadu_ps(&ptht[0]);
                          const  __m256 a   = _mm256_loadu_ps(&pa[0]);
                          const  __m256 ka04= _mm256_loadu_ps(&pka04[0]);
                          const  __m256 mm1r= _mm256_loadu_ps(&pmm1r[0]);
                          const  __m256 mm1i= _mm256_loadu_ps(&pmm1i[0]);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 aa = _mm256_mul_ps(a,a);
                          const  __m256 _2 = _mm256_set1_ps(2.0f); 
                          const  __m256 cost = _mm256_cos_ps(tht);
                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 cos2t,mms1r,mms1i,mma2r,mma2i;
                           __m256 divr,divi,rcs,frac,cabs;
                          cos2t   = _mm256_mul_ps(cost,cost);
                          mms1r   = _mm256_sub_ps(mm1r,_1);
                          mma2r   = _mm256_add_ps(mm1r,_2);
                          mms1i   = _mm256_sub_ps(mm1i,_1);
                          mma2i   = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm256_mul_ps(_4pi,aa);
                          cabs    =  cabs_ymm8c4(divr,divi);
                          rcs     = _mm256_mul_ps(_mm256_mul_ps(frac,cabs),
                                                  _mm256_mul_ps(cos2t,ka04));
                          return (rcs); 
                 }


                    __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3311_ymm8r4(const __m256 tht,
                                            const __m256 a,
                                            const __m256 ka04,
                                            const __m256 mm1r,
                                            const __m256 mm1i) {

                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 aa = _mm256_mul_ps(a,a);
                          const  __m256 _2 = _mm256_set1_ps(2.0f); 
                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 mms1r,mms1i,mma2r,mma2i;
                           __m256 divr,divi,rcs,frac,cabs;
                          mms1r   = _mm256_sub_ps(mm1r,_1);
                          mma2r   = _mm256_add_ps(mm1r,_2);
                          mms1i   = _mm256_sub_ps(mm1i,_1);
                          mma2i   = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm256_mul_ps(_4pi,aa);
                          cabs    =  cabs_ymm8c4(divr,divi);
                          rcs     = _mm256_mul_ps(ka04,_mm256_mul_ps(frac,cabs));
                                                  
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3311_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) ptht,
                                              const float * __restrict __ATTR_ALIGN__(32) pa,
                                              const float * __restrict __ATTR_ALIGN__(32) pka04,
                                              const float * __restrict __ATTR_ALIGN__(32) pmm1r,
                                              const float * __restrict __ATTR_ALIGN__(32) pmm1i) {

                          const  __m256 tht = _mm256_load_ps(&ptht[0]);
                          const  __m256 a   = _mm256_load_ps(&pa[0]);
                          const  __m256 ka04= _mm256_load_ps(&pka04[0]);
                          const  __m256 mm1r= _mm256_load_ps(&pmm1r[0]);
                          const  __m256 mm1i= _mm256_load_ps(&pmm1i[0]);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 aa = _mm256_mul_ps(a,a);
                          const  __m256 _2 = _mm256_set1_ps(2.0f); 
                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 mms1r,mms1i,mma2r,mma2i;
                           __m256 divr,divi,rcs,frac,cabs;
                          mms1r   = _mm256_sub_ps(mm1r,_1);
                          mma2r   = _mm256_add_ps(mm1r,_2);
                          mms1i   = _mm256_sub_ps(mm1i,_1);
                          mma2i   = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm256_mul_ps(_4pi,aa);
                          cabs    =  cabs_ymm8c4(divr,divi);
                          rcs     = _mm256_mul_ps(ka04,_mm256_mul_ps(frac,cabs));
                                                  
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3311_ymm8r4_u(const float * __restrict  ptht,
                                              const float * __restrict  pa,
                                              const float * __restrict  pka04,
                                              const float * __restrict  pmm1r,
                                              const float * __restrict  pmm1i) {

                          const  __m256 tht = _mm256_loadu_ps(&ptht[0]);
                          const  __m256 a   = _mm256_loadu_ps(&pa[0]);
                          const  __m256 ka04= _mm256_loadu_ps(&pka04[0]);
                          const  __m256 mm1r= _mm256_loadu_ps(&pmm1r[0]);
                          const  __m256 mm1i= _mm256_loadu_ps(&pmm1i[0]);
                          const  __m256 _1 = _mm256_set1_ps(1.0f);
                          const  __m256 aa = _mm256_mul_ps(a,a);
                          const  __m256 _2 = _mm256_set1_ps(2.0f); 
                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 mms1r,mms1i,mma2r,mma2i;
                           __m256 divr,divi,rcs,frac,cabs;
                          mms1r   = _mm256_sub_ps(mm1r,_1);
                          mma2r   = _mm256_add_ps(mm1r,_2);
                          mms1i   = _mm256_sub_ps(mm1i,_1);
                          mma2i   = _mm256_add_ps(mm1i,_2);
                          cdiv_ymm8c4(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm256_mul_ps(_4pi,aa);
                          cabs    =  cabs_ymm8c4(divr,divi);
                          rcs     = _mm256_mul_ps(ka04,_mm256_mul_ps(frac,cabs));
                                                  
                          return (rcs); 
                 }


                    /*
                          Bistatic Geometric Optics Rays.
                          The RCS of sphere included N-rays.
                          E-plane or H-plane, for backscattering
                          and forward-scattering E and H RCS are 
                          identical
                     */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3314_ymm8r4(const __m256 * __restrict Snthr,
                                            const __m256 * __restrict Snthi,
                                            const __m256 * __restrict cjphr,
                                            const __m256 * __restrict cjphi,
                                            __m256 * __restrict wrkr,
                                            __m256 * __restrict wrki,
                                            const __m256 k02,
                                            const int32_t N) {

                          const  __m256 _4pi = _mm256_set1_ps(12.566370614359172953850573533118f);
                           __m256 tmpre,tmpim;
                           __m256 rcs,cabs,frac;
                          frac = _mm256_div_ps(_4pi,k02);
                          accre = _mm256_setzero_ps();
                          accim = _mm256_setzero_ps();
                          for(int32_t i = 0; i != N; ++i) {
                              cmul_ymm8c4(Snthr[i],Snthi[i],cjphr[i],cjphi[i], 
                                           &wrkr[i],&wrki[i]);
                                                   
                          }
                          for(int32_t i = 0; i != N-1; ++i) {
                               tmpre = _mm256_add_ps(wrkr[i],wrkr[i+1]);
                               accre = _mm256_add_ps(accre,tmpre); 
                               tmpim = _mm256_add_ps(wrki[i],wrki[i+1]);
                               accim = _mm256_add_ps(accim,tmpim);
                          }
                          cabs =  cabs_ymm8c4(accre,accim);
                          rcs  = _mm256_mul_ps(frac,cabs);
                          return (rcs);
                }


                  /*
                         Large sphere limit, k0a > 1.15/m1 (reflective region).
                         Backscattering RCS, formula 3.3-17
                    */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3317_ymm8r4(const __m256 m1r,
                                            const __m256 m1i,
                                            const __m256 a) {

                          const  __m256 _1   = _mm256_set1_ps(1.0f);
                          const  __m256 frac = _mm256_mul_ps(PI,
                                                                 _mm256_mul_ps(a,a));
                           __m256 divr,divi,m1s1r,m1s1i;
                           __m256 m1a1r,m1a1i,cabs,rcs;
                          m1s1r = _mm256_sub_ps(m1r,_1);
                          m1a1r = _mm256_add_ps(m1r,_1);
                          m1s1i = _mm256_sub_ps(m1i,_1);
                          m1a1i = _mm256_add_ps(m1i,_1);
                          cdiv_ymm8c4(m1s1r,m1s1i,m1a1r,m1a1i,&divr,&divi); 
                          cabs = cabs_ymm8c4(divr,divi);
                          rcs  = _mm256_mul_ps(frac,cabs);
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3317_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pm1r,
                                              const float * __restrict __ATTR_ALIGN__(32) pm1i,
                                              const float * __restrict __ATTR_ALIGN__(32) pa) {

                          const  __m256 m1r  = _mm256_load_ps(&pm1r[0]);
                          const  __m256 m1i  = _mm256_load_ps(&pm1i[0]);
                          const  __m256 a    = _mm256_load_ps(&pa[0]);
                          const  __m256 _1   = _mm256_set1_ps(1.0f);
                          const  __m256 frac = _mm256_mul_ps(PI,
                                                                 _mm256_mul_ps(a,a));
                           __m256 divr,divi,m1s1r,m1s1i;
                           __m256 m1a1r,m1a1i,cabs,rcs;
                          m1s1r = _mm256_sub_ps(m1r,_1);
                          m1a1r = _mm256_add_ps(m1r,_1);
                          m1s1i = _mm256_sub_ps(m1i,_1);
                          m1a1i = _mm256_add_ps(m1i,_1);
                          cdiv_ymm8c4(m1s1r,m1s1i,m1a1r,m1a1i,&divr,&divi); 
                          cabs = cabs_ymm8c4(divr,divi);
                          rcs  = _mm256_mul_ps(frac,cabs);
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3317_ymm8r4_u(const float * __restrict  pm1r,
                                              const float * __restrict  pm1i,
                                              const float * __restrict  pa) {

                          const  __m256 m1r  = _mm256_loadu_ps(&pm1r[0]);
                          const  __m256 m1i  = _mm256_loadu_ps(&pm1i[0]);
                          const  __m256 a    = _mm256_loadu_ps(&pa[0]);
                          const  __m256 _1   = _mm256_set1_ps(1.0f);
                          const  __m256 frac = _mm256_mul_ps(PI,
                                                                 _mm256_mul_ps(a,a));
                           __m256 divr,divi,m1s1r,m1s1i;
                           __m256 m1a1r,m1a1i,cabs,rcs;
                          m1s1r = _mm256_sub_ps(m1r,_1);
                          m1a1r = _mm256_add_ps(m1r,_1);
                          m1s1i = _mm256_sub_ps(m1i,_1);
                          m1a1i = _mm256_add_ps(m1i,_1);
                          cdiv_ymm8c4(m1s1r,m1s1i,m1a1r,m1a1i,&divr,&divi); 
                          cabs = cabs_ymm8c4(divr,divi);
                          rcs  = _mm256_mul_ps(frac,cabs);
                          return (rcs);
                  }


                    /*
                       Forward scattering RCS.
                       Formula 3.3-19
                         */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3319_ymm8r4(const __m256 a,
                                            const __m256 k0a) {

                          const  __m256 aa   = _mm256_mul_ps(a,a);
                          const  __m256 c0   = _mm256_set1_ps(2.25f);
                          const  __m256 k0a23= _mm256_pow_ps(k0a,
                                                                 _mm256_set1_ps(0.666666666666666666666666666667f);
                           __m256 rcs,fac;
                          fac = _mm256_mul_ps(PI,aa);
                          rcs = _mm256_mul_ps(fac,_mm256_mul_ps(c0,k0a23));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3319_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pa,
                                            const float * __restrict __ATTR_ALIGN__(32) pk0a) {

                          const  __m256 a    = _mm256_load_ps(&pa[0]);
                          const  __m256 aa   = _mm256_mul_ps(a,a);
                          const  __m256 k0a  = _mm256_load_ps(&pk0a[0]);
                          const  __m256 c0   = _mm256_set1_ps(2.25f);
                          const  __m256 k0a23= _mm256_pow_ps(k0a,
                                                                 _mm256_set1_ps(0.666666666666666666666666666667f);
                           __m256 rcs,fac;
                          fac = _mm256_mul_ps(PI,aa);
                          rcs = _mm256_mul_ps(fac,_mm256_mul_ps(c0,k0a23));
                          return (rcs);
                 }


                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   __m256 rcs_f3319_ymm8r4_u(const float * __restrict  pa,
                                              const float * __restrict  pk0a) {

                          const  __m256 a    = _mm256_loadu_ps(&pa[0]);
                          const  __m256 aa   = _mm256_mul_ps(a,a);
                          const  __m256 k0a  = _mm256_loadu_ps(&pk0a[0]);
                          const  __m256 c0   = _mm256_set1_ps(2.25f);
                          const  __m256 k0a23= _mm256_pow_ps(k0a,
                                                                 _mm256_set1_ps(0.666666666666666666666666666667f);
                           __m256 rcs,fac;
                          fac = _mm256_mul_ps(PI,aa);
                          rcs = _mm256_mul_ps(fac,_mm256_mul_ps(c0,k0a23));
                          return (rcs);
                 }


                    /*
                         Approximate solutions for far-field region (Rayleigh-Gans)
                         (abs(m1-1) << 1,2*k0a abs(m1-1) << 1)
                         Bistatic scattering formula 3.3-22
                     */

                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3322_ymm8r4(const __m256 m1r,
                                         const __m256 m1i,
                                         const __m256 tht,
                                         const __m256 k0a,
                                         __m256 * __restrict S1r,
                                         __m256 * __restrict S1i) {
                       
                         __m256 cosht,m1s1r,m1s1i,htht,cost,t0,t1;
                         __m256 st,ct,carg,carg2,facr,faci,sinc;
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        cost  = _mm256_cos_ps(tht);
                        const  __m256 k0a3 = _mm256_mul_ps(k0a,
                                                           _mm256_mul_ps(k0a,k0a));
                        m1s1r = _mm256_sub_ps(m1r,_1);
                        m1s1i = _mm256_sub_ps(m1i,_1);
                        const  __m256 _2k0a  = _mm256_add_ps(k0a,k0a);
                        const  __m256 _2k0a3 = _mm256_add_ps(k0a3,k0a3)
                        htht = _mm256_mul_ps(_mm256_set1_ps(0.5f),tht);
                        cosht= _mm256_cos_ps(htht);
                        facr = _mm256_mul_ps(_2k0a3,m1s1r);
                        faci = _mm256_mul_ps(_2k0a3,m1s1i);
                        carg = _mm256_mul_ps(_2k0a,cosht);
                        carg2= _mm256_mul_ps(carg,carg);
                        st   = _mm256_sin_ps(carg);
                        sinc = _mm256_div_ps(st,carg);
                        ct   = _mm256_cos_ps(carg);
                        t0   = _mm256_mul_ps(_mm256_sub_ps(st,ct),cost);
                        t1   = _mm256_div_ps(t0,carg2);
                        *S1r = _mm256_mul_ps(facr,t1);
                        *S1i = _mm256_mul_ps(faci,t1);
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3322_ymm8r4_a(const float * __restrict __ATTR_ALIGN__(32) pm1r,
                                           const float * __restrict __ATTR_ALIGN__(32) pm1i,
                                           const float * __restrict __ATTR_ALIGN__(32) ptht,
                                           const float * __restrict __ATTR_ALIGN__(32) pk0a,
                                           float * __restrict __ATTR_ALIGN__(32) S1r,
                                           float * __restrict __ATTR_ALIGN__(32) S1i) {
                       
                        const  __m256 m1r = _mm256_load_ps(&pm1r[0]);
                        const  __m256 m1i = _mm256_load_ps(&pm1i[0]);
                        const  __m256 tht = _mm256_load_ps(&ptht[0]);
                        const  __m256 k0a = _mm256_load_ps(&pk0a[0]);
                         __m256 cosht,m1s1r,m1s1i,htht,cost,t0,t1;
                         __m256 st,ct,carg,carg2,facr,faci,sinc;
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        cost  = _mm256_cos_ps(tht);
                        const  __m256 k0a3 = _mm256_mul_ps(k0a,
                                                           _mm256_mul_ps(k0a,k0a));
                        m1s1r = _mm256_sub_ps(m1r,_1);
                        m1s1i = _mm256_sub_ps(m1i,_1);
                        const  __m256 _2k0a  = _mm256_add_ps(k0a,k0a);
                        const  __m256 _2k0a3 = _mm256_add_ps(k0a3,k0a3)
                        htht = _mm256_mul_ps(_mm256_set1_ps(0.5f),tht);
                        cosht= _mm256_cos_ps(htht);
                        facr = _mm256_mul_ps(_2k0a3,m1s1r);
                        faci = _mm256_mul_ps(_2k0a3,m1s1i);
                        carg = _mm256_mul_ps(_2k0a,cosht);
                        carg2= _mm256_mul_ps(carg,carg);
                        st   = _mm256_sin_ps(carg);
                        sinc = _mm256_div_ps(st,carg);
                        ct   = _mm256_cos_ps(carg);
                        t0   = _mm256_mul_ps(_mm256_sub_ps(st,ct),cost);
                        t1   = _mm256_div_ps(t0,carg2);
                        _mm256_store_ps(&S1r[0], _mm256_mul_ps(facr,t1));
                        _mm256_store_ps(&S1i[0], _mm256_mul_ps(faci,t1));
                 }


                   __ATTR_ALWAYS_INLINE__
	           
	          
                  
	           static inline
                   void S1_f3322_ymm8r4_u(const float * __restrict  pm1r,
                                           const float * __restrict  pm1i,
                                           const float * __restrict  ptht,
                                           const float * __restrict  pk0a,
                                           float * __restrict S1r,
                                           float * __restrict  S1i) {
                       
                        const  __m256 m1r = _mm256_loadu_ps(&pm1r[0]);
                        const  __m256 m1i = _mm256_loadu_ps(&pm1i[0]);
                        const  __m256 tht = _mm256_loadu_ps(&ptht[0]);
                        const  __m256 k0a = _mm256_loadu_ps(&pk0a[0]);
                         __m256 cosht,m1s1r,m1s1i,htht,cost,t0,t1;
                         __m256 st,ct,carg,carg2,facr,faci,sinc;
                        const  __m256 _1   = _mm256_set1_ps(1.0f);
                        cost  = _mm256_cos_ps(tht);
                        const  __m256 k0a3 = _mm256_mul_ps(k0a,
                                                           _mm256_mul_ps(k0a,k0a));
                        m1s1r = _mm256_sub_ps(m1r,_1);
                        m1s1i = _mm256_sub_ps(m1i,_1);
                        const  __m256 _2k0a  = _mm256_add_ps(k0a,k0a);
                        const  __m256 _2k0a3 = _mm256_add_ps(k0a3,k0a3)
                        htht = _mm256_mul_ps(_mm256_set1_ps(0.5f),tht);
                        cosht= _mm256_cos_ps(htht);
                        facr = _mm256_mul_ps(_2k0a3,m1s1r);
                        faci = _mm256_mul_ps(_2k0a3,m1s1i);
                        carg = _mm256_mul_ps(_2k0a,cosht);
                        carg2= _mm256_mul_ps(carg,carg);
                        st   = _mm256_sin_ps(carg);
                        sinc = _mm256_div_ps(st,carg);
                        ct   = _mm256_cos_ps(carg);
                        t0   = _mm256_mul_ps(_mm256_sub_ps(st,ct),cost);
                        t1   = _mm256_div_ps(t0,carg2);
                        _mm256_storeu_ps(&S1r[0], _mm256_mul_ps(facr,t1));
                        _mm256_storeu_ps(&S1i[0], _mm256_mul_ps(faci,t1));
                 }


                  

     } // radiolocation

} // gms











#endif /*__GMS_RCS_SPHERE_YMM8R4_HPP__*/
