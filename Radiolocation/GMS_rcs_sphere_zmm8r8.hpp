
#ifndef __GMS_RCS_SPHERE_ZMM8R8_HPP__
#define __GMS_RCS_SPHERE_ZMM8R8_HPP__


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

    const unsigned int GMS_RCS_SPHERE_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_RCS_SPHERE_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_RCS_SPHERE_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_RCS_SPHERE_ZMM8R8_FULLVER =
      1000U*GMS_RCS_SPHERE_ZMM8R8_MAJOR+
      100U*GMS_RCS_SPHERE_ZMM8R8_MINOR+
      10U*GMS_RCS_SPHERE_ZMM8R8_MICRO;
    const char * const GMS_RCS_SPHERE_ZMM8R8_CREATION_DATE = "04-01-2023 12:45 AM +00200 (WED 04 JAN 2023 GMT+2)";
    const char * const GMS_RCS_SPHERE_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_SPHERE_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_SPHERE_ZMM8R8_DESCRIPTION   = "AVX512 optimized Sphere Radar Cross Section (analytic) functionality."

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
                   const  static __m512d Ir  = _mm512_setzero_pd();
                   const  static __m512d Ii  = _mm512_set1_pd(1.0f);
                   const  static __m512d nIr = _mm512_set1_pd(-0.0f);
                   const  static __m512d nIi = _mm512_set1_pd(-1.0f);
                   const  static __m512d PI  = _mm512_set1_pd(3.14159265358979323846264338328);

               }

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm8r8(const __m512d k0a, // size of sphere expressed in wavenumber units
                                               __m512d * __restrict F0r, // the results
                                               __m512d * __restrict F0i) { // the results

                        const register __m512d _1 = _mm512_set1_pd(1.0);
                        const register __m512d c0 = _mm512_set1_pd(1.5);
                        const register __m512d c1 = _mm512_set1_pd(0.092592592592592592592592592593);
                        const register __m512d c2 = _mm512_set1_pd(0.018888888888888888888888888889);
                        const register __m512d c3 = _mm512_set1_pd(0.558656504577139497774418409339);
                        const register __m512d c4 = _mm512_set1_pd(1.2);
                        const register __m512d c5 = _mm512_set1_pd(0.5);
                        __m512d k0a2,k0a3,k0a4,k0a6;
                        __m512d t0,t1,t2,t3;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        t1   = _mm512_sub_pd(_1,_mm512_mul_pd(c1,k0a2));
                        k0a3 = _mm512_mul_pd(k0a2,k0a);
                        t0   = _mm512_mul_pd(c0,k0a3);
                        k0a4 = _mm512_mul_pd(k0a3,k0a);
                        k0a6 = _mm512_mul_pd(k0a3,k0a2);
                        t2   = _mm512_fmsub_pd(c2,k0a4,_mm512_mul_pd(c3,k0a6));
                        t3   = _mm512_mul_pd(t0,_mm512_add_pd(t1,t1));
                        *F0r = t3;
                        t0   = _mm512_mul_pd(c5,k0a6);
                        t1   = _mm512_add_pd(_1,_mm512_mul_pd(c4,k0a2));
                        *F0i = _mm512_mul_pd(t0,t1);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a, // size of sphere expressed in wavenumber units
                                               __m512d * __restrict F0r, // the results
                                               __m512d * __restrict F0i) { // the results

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        const register __m512d _1 = _mm512_set1_pd(1.0);
                        const register __m512d c0 = _mm512_set1_pd(1.5);
                        const register __m512d c1 = _mm512_set1_pd(0.092592592592592592592592592593);
                        const register __m512d c2 = _mm512_set1_pd(0.018888888888888888888888888889);
                        const register __m512d c3 = _mm512_set1_pd(0.558656504577139497774418409339);
                        const register __m512d c4 = _mm512_set1_pd(1.2);
                        const register __m512d c5 = _mm512_set1_pd(0.5f);
                        __m512d k0a2,k0a3,k0a4,k0a6;
                        __m512d t0,t1,t2,t3;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        t1   = _mm512_sub_pd(_1,_mm512_mul_pd(c1,k0a2));
                        k0a3 = _mm512_mul_pd(k0a2,k0a);
                        t0   = _mm512_mul_pd(c0,k0a3);
                        k0a4 = _mm512_mul_pd(k0a3,k0a);
                        k0a6 = _mm512_mul_pd(k0a3,k0a2);
                        t2   = _mm512_fmsub_pd(c2,k0a4,_mm512_mul_pd(c3,k0a6));
                        t3   = _mm512_mul_pd(t0,_mm512_add_pd(t1,t1));
                        *F0r = t3;
                        t0   = _mm512_mul_pd(c5,k0a6);
                        t1   = _mm512_add_pd(_1,_mm512_mul_pd(c4,k0a2));
                        *F0i = _mm512_mul_pd(t0,t1);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm8r8_u(const double * __restrict  pk0a, // size of sphere expressed in wavenumber units
                                               __m512d * __restrict F0r, // the results
                                               __m512d * __restrict F0i) { // the results

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        const register __m512d c0 = _mm512_set1_pd(1.5f);
                        const register __m512d c1 = _mm512_set1_pd(0.092592592592592592592592592593);
                        const register __m512d c2 = _mm512_set1_pd(0.018888888888888888888888888889);
                        const register __m512d c3 = _mm512_set1_pd(0.558656504577139497774418409339);
                        const register __m512d c4 = _mm512_set1_pd(1.2f);
                        const register __m512d c5 = _mm512_set1_pd(0.5f);
                        __m512d k0a2,k0a3,k0a4,k0a6;
                        __m512d t0,t1,t2,t3;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        t1   = _mm512_sub_pd(_1,_mm512_mul_pd(c1,k0a2));
                        k0a3 = _mm512_mul_pd(k0a2,k0a);
                        t0   = _mm512_mul_pd(c0,k0a3);
                        k0a4 = _mm512_mul_pd(k0a3,k0a);
                        k0a6 = _mm512_mul_pd(k0a3,k0a2);
                        t2   = _mm512_fmsub_pd(c2,k0a4,_mm512_mul_pd(c3,k0a6));
                        t3   = _mm512_mul_pd(t0,_mm512_add_pd(t1,t1));
                        *F0r = t3;
                        t0   = _mm512_mul_pd(c5,k0a6);
                        t1   = _mm512_add_pd(_1,_mm512_mul_pd(c4,k0a2));
                        *F0i = _mm512_mul_pd(t0,t1);
                }


                  /*
                        Radar Cross Section Handbook 1, page 147, formula 3.2-5
                        Backscattering cross section
                        
                    */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f325_zmm8r8(const __m512d k0,
                                           const __m512d a ) {
                                        
                        register __m512d a2,k0a,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,sigma; 
                        a2  = _mm512_mul_pd(a,a);
                        const register __m512d pi9 = _mm512_set1_pd(28.274333882308139146163790449516);
                        k0a = _mm512_mul_pd(k0,a); 
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        k0a2= _mm512_mul_pd(k0a,k0a);
                        const register __m512d c0 = _mm512_set1_pd(0.185185185185185185185185185185);
                        k0a4= _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c1 = _mm512_set1_pd(0.04635116598079561042524005487);
                        k0a6= _mm512_mul_pd(k0a4,_mm512_mul_pd(k0a,k0a));
                        const register __m512d c2 = _mm512_set1_pd(1.00969984042999916015789031662);
                        t0  = _mm512_mul_pd(k0a4,_mm512_mul_pd(pi9,a2));
                        t1  = _mm512_sub_pd(_1,_mm512_mul_pd(c0,k0a2));
                        t2  = _mm512_fmsub_pd(c1,k0a4,_mm512_mul_pd(c2,k0a6));
                        sigma = _mm512_fmadd_pd(t1,t2,t0);
                        return (sigma);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f325_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pa ) {
                        
                        register __m512d k0= _mm512_load_pd(&pk0[0]); 
                        register __m512d a = _mm512_load_pd(&pa[0]);               
                        register __m512d a2,k0a,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,sigma; 
                        a2  = _mm512_mul_pd(a,a);
                        const register __m512d pi9 = _mm512_set1_pd(28.274333882308139146163790449516);
                        k0a = _mm512_mul_pd(k0,a); 
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        k0a2= _mm512_mul_pd(k0a,k0a);
                        const register __m512d c0 = _mm512_set1_pd(0.185185185185185185185185185185);
                        k0a4= _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c1 = _mm512_set1_pd(0.04635116598079561042524005487);
                        k0a6= _mm512_mul_pd(k0a4,_mm512_mul_pd(k0a,k0a));
                        const register __m512d c2 = _mm512_set1_pd(1.00969984042999916015789031662);
                        t0  = _mm512_mul_pd(k0a4,_mm512_mul_pd(pi9,a2));
                        t1  = _mm512_sub_pd(_1,_mm512_mul_pd(c0,k0a2));
                        t2  = _mm512_fmsub_pd(c1,k0a4,_mm512_mul_pd(c2,k0a6));
                        sigma = _mm512_fmadd_pd(t1,t2,t0);
                        return (sigma);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f325_zmm8r8_u(const double * __restrict pk0,
                                             const double * __restrict pa ) {
                        
                        register __m512d k0= _mm512_loadu_pd(&pk0[0]); 
                        register __m512d a = _mm512_loadu_pd(&pa[0]);               
                        register __m512d a2,k0a,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,sigma; 
                        a2  = _mm512_mul_pd(a,a);
                        const register __m512d pi9 = _mm512_set1_pd(28.274333882308139146163790449516);
                        k0a = _mm512_mul_pd(k0,a); 
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        k0a2= _mm512_mul_pd(k0a,k0a);
                        const register __m512d c0 = _mm512_set1_pd(0.185185185185185185185185185185);
                        k0a4= _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c1 = _mm512_set1_pd(0.04635116598079561042524005487);
                        k0a6= _mm512_mul_pd(k0a4,_mm512_mul_pd(k0a,k0a));
                        const register __m512d c2 = _mm512_set1_pd(1.00969984042999916015789031662);
                        t0  = _mm512_mul_pd(k0a4,_mm512_mul_pd(pi9,a2));
                        t1  = _mm512_sub_pd(_1,_mm512_mul_pd(c0,k0a2));
                        t2  = _mm512_fmsub_pd(c1,k0a4,_mm512_mul_pd(c2,k0a6));
                        sigma = _mm512_fmadd_pd(t1,t2,t0);
                        return (sigma);
                }



#include "GMS_complex_zmm8r8.hpp"


                  /*
                        Creeping wave term, F c(0) at the upper end of the resonance region.
                        Formula 3.2-8
                    */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_f328_zmm8r8(const __m512d x,//k0a
                                         __m512d * __restrict Fc0r,
                                         __m512d * __restrict Fc0i) {

                      register __m512d p43,pn16,pn13,pn23,p13;
                      register __m512d t0r,t0i,t1r,t1i,t2r,t2i;  
                      register __m512d t3r,t3i,t4r,t4i,t5r,t5i;
                      register __m512d t6r,t6i,e1r,e1i,e2r,e2i;
                      register __m512d e3r,e3i,t6r,t6i,t7r,t7i;
                      register __m512d argr,argi,resr,resi;              
                      const register __m512d c0  = _mm512_set1_pd(1.357588f);
                      p43                       = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(1.3333333333333333333333333333333));
                      t0r                       = _mm512_mul_pd(nIr,p43); // - / x 4/3
                      t0i                       = _mm512_mul_pd(nIi,p43); // - / x 4/3
                      const register __m512d c0r = _mm512_set1_pd(0.741196);
                      const register __m512d c0i = _mm512_set1_pd(1.283788);
                      pn16                      = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(0.166666666666666666666666666667));
                      pn16                      = _mm512_rcp14_pd(pn16);
                      const register __m512d c1r = _mm512_set1_pd(2.200002);
                      const register __m512d c1i = _mm512_set1_pd(-1.270172);
                      pn23                      = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(0.666666666666666666666666666667));
                      pn23                      = _mm512_rcp14_pd(pn23);
                      const register __m512d c2r = _mm512_set1_pd(0.445396);
                      const register __m512d c2i = _mm512_set1_pd(0.257150);
                      p13                       = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(-0.3333333333333333333333333333333333));
                      const register __m512d c3r = _mm512_set1_pd(0.964654);
                      const register __m512d c3i = _mm512_set1_pd(1.670829);
                      pn13                      = _mm512_rcp14_pd(p13);
                      const register __m512d c4r = _mm512_set1_pd(7.014224);
                      const register __m512d c4i = _mm512_set1_pd(-4.049663);
                      t1r                       = _mm512_set1_pd(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_pd(pn16,_mm512_set1_pd(-3.14159265358979323846264338328)); //-//-//
                      const register __m512d c5r = _mm512_set1_pd(0.444477);
                      const register __m512d c5i = _mm512_set1_pd(0.256619);
                      t2r                       = _mm512_fmadd_pd(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_pd(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512d c6r = _mm512_set1_pd(0.798821);
                      const register __m512d c6i = _mm512_set1_pd(1.383598);
                      t3r                       = _mm512_mul_pd(pn13,c1r);
                      t3i                       = _mm512_mul_pd(pn13,c1i);
                      const register __m512d c7r = _mm512_set1_pd(5.048956);
                      const register __m512d c7i = _mm512_set1_pd(-2.915016);
                      t4r                       = _mm512_mul_pd(pn13,c2r);
                      t4i                       = _mm512_mul_pd(pn13,c2i);
                      const register __m512d c8r = _mm512_set1_pd(0.312321);
                      const register __m512d c8i = _mm512_set1_pd(0.180319);
                      t5r                       = _mm512_add_pd(t3r,t4r);
                      t5i                       = _mm512_add_pd(t3i,t4i);
                      cexp_zmm8r8(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm512_fmadd_pd(c3r,pn23,_mm512_set1_pd(0.695864)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm512_fmadd_pd(c3i,pn23,_mm512_set1_pd(0.695864)); //--//--//--//
                      t4r                       = _mm512_mul_pd(p13,c4r);
                      t4i                       = _mm512_mul_pd(p13,c4i);
                      t5r                       = _mm512_mul_pd(pn13,c5r);
                      t5i                       = _mm512_mul_pd(pn13,c5i);
                      t6r                       = _mm512_sub_pd(t4r,t5r);
                      t6i                       = _mm512_sub_pd(t4i,t5i);
                      cexp_zmm8r8(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm512_fmadd_pd(c4r,pn23,_mm512_set1_pd(0.807104)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm512_fmadd_pd(c4i,pn23,_mm512_set1_pd(0.807104)); //--//--//
                      t5r                       = _mm512_mul_pd(p13,c7r);
                      t5i                       = _mm512_mul_pd(p13,c7i);
                      t6r                       = _mm512_mul_pd(pn13,c8r);
                      t6i                       = _mm512_mul_pd(pn13,c8i);
                      t7r                       = _mm512_sub_pd(t5r,t6r);
                      t7i                       = _mm512_sub_pd(t5i,t6i);
                      cexp_zmm8r8(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm512_setzero_pd();
                      t5i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm512_setzero_pd();
                      t6i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm512_setzero_pd();
                      t7i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm512_add_pd(t5r,_mm512_sub_pd(t6r,t7r));
                      argi                      = _mm512_add_pd(t5i,_mm512_sub_pd(t6i,t7i));
                      t4r                       = _mm512_setzero_pd();
                      t4i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_zmm8r8(t4r,t4i,&resr,&resi);
                      cmul_zmm8r8(t04,t0i,resr,resi,&Fc0r,&Fc0i);                    
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_f328_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px,//k0a
                                           double * __restrict __ATTR_ALIGN__(64) Fc0r,
                                           double * __restrict __ATTR_ALIGN__(64) Fc0i) {

                      register __m512d x = _mm512_load_pd(&px[0]);
                      register __m512d p43,pn16,pn13,pn23,p13;
                      register __m512d t0r,t0i,t1r,t1i,t2r,t2i;  
                      register __m512d t3r,t3i,t4r,t4i,t5r,t5i;
                      register __m512d t6r,t6i,e1r,e1i,e2r,e2i;
                      register __m512d e3r,e3i,t6r,t6i,t7r,t7i;
                      register __m512d argr,argi,resr,resi; 
                      register __m512d fcr,fci;             
                      const register __m512d c0  = _mm512_set1_pd(1.357588);
                      p43                       = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(1.3333333333333333333333333333333));
                      t0r                       = _mm512_mul_pd(nIr,p43); // - / x 4/3
                      t0i                       = _mm512_mul_pd(nIi,p43); // - / x 4/3
                      const register __m512d c0r = _mm512_set1_pd(0.741196);
                      const register __m512d c0i = _mm512_set1_pd(1.283788);
                      pn16                      = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(0.166666666666666666666666666667));
                      pn16                      = _mm512_rcp14_pd(pn16);
                      const register __m512d c1r = _mm512_set1_pd(2.200002);
                      const register __m512d c1i = _mm512_set1_pd(-1.270172);
                      pn23                      = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(0.666666666666666666666666666667));
                      pn23                      = _mm512_rcp14_pd(pn23);
                      const register __m512d c2r = _mm512_set1_pd(0.445396);
                      const register __m512d c2i = _mm512_set1_pd(0.257150);
                      p13                       = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(-0.3333333333333333333333333333333333));
                      const register __m512d c3r = _mm512_set1_pd(0.964654);
                      const register __m512d c3i = _mm512_set1_pd(1.670829);
                      pn13                      = _mm512_rcp14_pd(p13);
                      const register __m512d c4r = _mm512_set1_pd(7.014224);
                      const register __m512d c4i = _mm512_set1_pd(-4.049663);
                      t1r                       = _mm512_set1_pd(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_pd(pn16,_mm512_set1_pd(-3.14159265358979323846264338328)); //-//-//
                      const register __m512d c5r = _mm512_set1_pd(0.444477);
                      const register __m512d c5i = _mm512_set1_pd(0.256619);
                      t2r                       = _mm512_fmadd_pd(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_pd(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512d c6r = _mm512_set1_pd(0.798821);
                      const register __m512d c6i = _mm512_set1_pd(1.383598);
                      t3r                       = _mm512_mul_pd(pn13,c1r);
                      t3i                       = _mm512_mul_pd(pn13,c1i);
                      const register __m512d c7r = _mm512_set1_pd(5.048956);
                      const register __m512d c7i = _mm512_set1_pd(-2.915016);
                      t4r                       = _mm512_mul_pd(pn13,c2r);
                      t4i                       = _mm512_mul_pd(pn13,c2i);
                      const register __m512d c8r = _mm512_set1_pd(0.312321);
                      const register __m512d c8i = _mm512_set1_pd(0.180319);
                      t5r                       = _mm512_add_pd(t3r,t4r);
                      t5i                       = _mm512_add_pd(t3i,t4i);
                      cexp_zmm8r8(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm512_fmadd_pd(c3r,pn23,_mm512_set1_pd(0.695864)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm512_fmadd_pd(c3i,pn23,_mm512_set1_pd(0.695864)); //--//--//--//
                      t4r                       = _mm512_mul_pd(p13,c4r);
                      t4i                       = _mm512_mul_pd(p13,c4i);
                      t5r                       = _mm512_mul_pd(pn13,c5r);
                      t5i                       = _mm512_mul_pd(pn13,c5i);
                      t6r                       = _mm512_sub_pd(t4r,t5r);
                      t6i                       = _mm512_sub_pd(t4i,t5i);
                      cexp_zmm8r8(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm512_fmadd_pd(c4r,pn23,_mm512_set1_pd(0.807104)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm512_fmadd_pd(c4i,pn23,_mm512_set1_pd(0.807104)); //--//--//
                      t5r                       = _mm512_mul_pd(p13,c7r);
                      t5i                       = _mm512_mul_pd(p13,c7i);
                      t6r                       = _mm512_mul_pd(pn13,c8r);
                      t6i                       = _mm512_mul_pd(pn13,c8i);
                      t7r                       = _mm512_sub_pd(t5r,t6r);
                      t7i                       = _mm512_sub_pd(t5i,t6i);
                      cexp_zmm8r8(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm512_setzero_pd();
                      t5i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm512_setzero_pd();
                      t6i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm512_setzero_pd();
                      t7i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm512_add_pd(t5r,_mm512_sub_pd(t6r,t7r));
                      argi                      = _mm512_add_pd(t5i,_mm512_sub_pd(t6i,t7i));
                      t4r                       = _mm512_setzero_pd();
                      t4i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_zmm8r8(t4r,t4i,&resr,&resi);
                      cmul_zmm8r8(t04,t0i,resr,resi,&fcr,&fci);    
                      __m512d_store_pd(&Fc0r[0], fcr);
                      __m512d_store_pd(&Fc0i[0], fci);                
                }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_f328_zmm8r8_u(const double * __restrict px,//k0a
                                           double * __restrict  Fc0r,
                                           double * __restrict  Fc0i) {

                      register __m512d x = _mm512_loadu_pd(&px[0]);
                      register __m512d p43,pn16,pn13,pn23,p13;
                      register __m512d t0r,t0i,t1r,t1i,t2r,t2i;  
                      register __m512d t3r,t3i,t4r,t4i,t5r,t5i;
                      register __m512d t6r,t6i,e1r,e1i,e2r,e2i;
                      register __m512d e3r,e3i,t6r,t6i,t7r,t7i;
                      register __m512d argr,argi,resr,resi; 
                      register __m512d fcr,fci;             
                      const register __m512d c0  = _mm512_set1_pd(1.357588);
                      p43                       = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(1.3333333333333333333333333333333));
                      t0r                       = _mm512_mul_pd(nIr,p43); // - / x 4/3
                      t0i                       = _mm512_mul_pd(nIi,p43); // - / x 4/3
                      const register __m512d c0r = _mm512_set1_pd(0.741196);
                      const register __m512d c0i = _mm512_set1_pd(1.283788);
                      pn16                      = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(0.166666666666666666666666666667));
                      pn16                      = _mm512_rcp14_pd(pn16);
                      const register __m512d c1r = _mm512_set1_pd(2.200002);
                      const register __m512d c1i = _mm512_set1_pd(-1.270172);
                      pn23                      = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(0.666666666666666666666666666667));
                      pn23                      = _mm512_rcp14_pd(pn23);
                      const register __m512d c2r = _mm512_set1_pd(0.445396);
                      const register __m512d c2i = _mm512_set1_pd(0.257150);
                      p13                       = _mm512_pow_pd(x,
                                                        _mm512_set1_pd(-0.3333333333333333333333333333333333));
                      const register __m512d c3r = _mm512_set1_pd(0.964654);
                      const register __m512d c3i = _mm512_set1_pd(1.670829);
                      pn13                      = _mm512_rcp14_pd(p13);
                      const register __m512d c4r = _mm512_set1_pd(7.014224);
                      const register __m512d c4i = _mm512_set1_pd(-4.049663);
                      t1r                       = _mm512_set1_pd(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_pd(pn16,_mm512_set1_pd(-3.14159265358979323846264338328)); //-//-//
                      const register __m512d c5r = _mm512_set1_pd(0.444477);
                      const register __m512d c5i = _mm512_set1_pd(0.256619);
                      t2r                       = _mm512_fmadd_pd(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_pd(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512d c6r = _mm512_set1_pd(0.798821);
                      const register __m512d c6i = _mm512_set1_pd(1.383598);
                      t3r                       = _mm512_mul_pd(pn13,c1r);
                      t3i                       = _mm512_mul_pd(pn13,c1i);
                      const register __m512d c7r = _mm512_set1_pd(5.048956);
                      const register __m512d c7i = _mm512_set1_pd(-2.915016);
                      t4r                       = _mm512_mul_pd(pn13,c2r);
                      t4i                       = _mm512_mul_pd(pn13,c2i);
                      const register __m512d c8r = _mm512_set1_pd(0.312321);
                      const register __m512d c8i = _mm512_set1_pd(0.180319);
                      t5r                       = _mm512_add_pd(t3r,t4r);
                      t5i                       = _mm512_add_pd(t3i,t4i);
                      cexp_zmm8r8(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm512_fmadd_pd(c3r,pn23,_mm512_set1_pd(0.695864)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm512_fmadd_pd(c3i,pn23,_mm512_set1_pd(0.695864)); //--//--//--//
                      t4r                       = _mm512_mul_pd(p13,c4r);
                      t4i                       = _mm512_mul_pd(p13,c4i);
                      t5r                       = _mm512_mul_pd(pn13,c5r);
                      t5i                       = _mm512_mul_pd(pn13,c5i);
                      t6r                       = _mm512_sub_pd(t4r,t5r);
                      t6i                       = _mm512_sub_pd(t4i,t5i);
                      cexp_zmm8r8(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm512_fmadd_pd(c4r,pn23,_mm512_set1_pd(0.807104)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm512_fmadd_pd(c4i,pn23,_mm512_set1_pd(0.807104)); //--//--//
                      t5r                       = _mm512_mul_pd(p13,c7r);
                      t5i                       = _mm512_mul_pd(p13,c7i);
                      t6r                       = _mm512_mul_pd(pn13,c8r);
                      t6i                       = _mm512_mul_pd(pn13,c8i);
                      t7r                       = _mm512_sub_pd(t5r,t6r);
                      t7i                       = _mm512_sub_pd(t5i,t6i);
                      cexp_zmm8r8(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm512_setzero_pd();
                      t5i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm512_setzero_pd();
                      t6i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm512_setzero_pd();
                      t7i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm512_add_pd(t5r,_mm512_sub_pd(t6r,t7r));
                      argi                      = _mm512_add_pd(t5i,_mm512_sub_pd(t6i,t7i));
                      t4r                       = _mm512_setzero_pd();
                      t4i                       = _mm512_setzero_pd();
                      cmul_zmm8r8(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_zmm8r8(t4r,t4i,&resr,&resi);
                      cmul_zmm8r8(t04,t0i,resr,resi,&fcr,&fci);    
                      __m512d_storeu_pd(&Fc0r[0], fcr);
                      __m512d_storeu_pd(&Fc0i[0], fci);                
                }




                  

#include "GMS_sleefsimddp.hpp"

                  /*
                       The complex scattering amplitudes near the lower end of the resonance
                       region (say, 0.4 < k^a < 1) 
                       E-plane approximation
                       These equations are not valid at all for k0a > 1. They are
                       valid for all theta angles.
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3213_zmm146r4( const __m512d k0a,
                                           const __m512d tht,
                                           __m512d * __restrict S1r,
                                           __m512d * __restrict S1i) {

                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0f);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(0.244444444444444444444444444444);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.083333333333333333333333333333);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/120.0);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(1807.0f/70.0f);
                        t0                         = _mm512_add_pd(half,cost);  // 0.5+cost
                        const register __m512d c5   = _mm512_set1_pd(2531.0f/105.0f);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmadd_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0f/42.0f);
                        t2                         = _mm512_sub_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c0,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2);
                        t4                         = _mm512_mul_pd(c7,_mm512_fmsub_pd(_4,cost,_1));
                        t5                         = _mm512_mul_pd(c8,_mm512_fmadd_pd(_2,cost,_1));
                        t5                         = _mm512_mul_pd(t5,k0a2);
                        t6                         = _mm512_mul_pd(k0a6,_mm512_add_pd(t4,t5)); // imaginary part
                        *S1i                       = t6;
                        t4                         = _mm512_sub_pd(t0,_mm512_add_pd(t2,t3));
                        *S1r                       = _mm512_mul_pd(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3213_zmm146r4_a( const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht,
                                             double * __restrict __ATTR_ALIGN__(64) S1r,
                                             double * __restrict __ATTR_ALIGN__(64) S1i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0f);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(0.244444444444444444444444444444);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.083333333333333333333333333333);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0f/120.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(1807.0f/70.0f);
                        t0                         = _mm512_add_pd(half,cost);  // 0.5+cost
                        const register __m512d c5   = _mm512_set1_pd(2531.0f/105.0f);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmadd_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0f/42.0f);
                        t2                         = _mm512_sub_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c0,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2f);
                        t4                         = _mm512_mul_pd(c7,_mm512_fmsub_pd(_4,cost,_1));
                        t5                         = _mm512_mul_pd(c8,_mm512_fmadd_pd(_2,cost,_1));
                        t5                         = _mm512_mul_pd(t5,k0a2);
                        t6                         = _mm512_mul_pd(k0a6,_mm512_add_pd(t4,t5)); // imaginary part
                        _mm512_store_pd(&S1i[0],t6);
                        
                        t4                         = _mm512_sub_pd(t0,_mm512_add_pd(t2,t3));
                        _mm512_store_pd(&S1r[0],_mm512_mul_pd(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3213_zmm146r4_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict S1r,
                                             double * __restrict S1i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0f);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(0.244444444444444444444444444444);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.083333333333333333333333333333);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0f/120.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(1807.0f/70.0f);
                        t0                         = _mm512_add_pd(half,cost);  // 0.5+cost
                        const register __m512d c5   = _mm512_set1_pd(2531.0f/105.0f);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmadd_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0f/42.0f);
                        t2                         = _mm512_sub_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c0,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2f);
                        t4                         = _mm512_mul_pd(c7,_mm512_fmsub_pd(_4,cost,_1));
                        t5                         = _mm512_mul_pd(c8,_mm512_fmadd_pd(_2,cost,_1));
                        t5                         = _mm512_mul_pd(t5,k0a2);
                        t6                         = _mm512_mul_pd(k0a6,_mm512_add_pd(t4,t5)); // imaginary part
                        _mm512_storeu_pd(&S1i[0],t6);
                        
                        t4                         = _mm512_sub_pd(t0,_mm512_add_pd(t2,t3));
                        _mm512_storeu_pd(&S1r[0],_mm512_mul_pd(k0a3,t4);
                }


                   /*
                       The complex scattering amplitudes near the lower end of the resonance
                       region (say, 0.4 < k^a < 1) 
                       H-plane approximation
                       These equations are not valid at all for k0a > 1. They are
                       valid for all theta angles.
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3214_zmm146r4( const __m512d k0a,
                                           const __m512d tht,
                                           __m512d * __restrict S2r,
                                           __m512d * __restrict S2i) {

                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d c9   = _mm512_set1_pd(0.125);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0f);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(23.0f/60.0f);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.055555555555555555555555555556);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0f/60.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(-1343.0f/150.0f);
                        t0                         = _mm512_fmadd_pd(half,cost,_1);  // 0.5*cost+1
                        const register __m512d c5   = _mm512_set1_pd(3769.0f/280.0f);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmsub_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0f/63.0f);
                        t2                         = _mm512_add_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c9,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2f);
                        t4                         = _mm512_mul_pd(c7,_mm512_sub_pd(_4,cost));
                        t5                         = _mm512_mul_pd(c8,_mm512_add_pd(_2,cost));
                        t5                         = _mm512_mul_pd(t5,k0a2);
                        t6                         = _mm512_mul_pd(k0a6,_mm512_add_pd(t4,t5)); // imaginary part
                        *S2i                       = t6;
                        t4                         = _mm512_add_pd(t0,_mm512_add_pd(t2,t3));
                        *S2r                       = _mm512_mul_pd(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3214_zmm146r4_a( const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht,
                                             double * __restrict __ATTR_ALIGN__(64) S2r,
                                             double * __restrict __ATTR_ALIGN__(64) S2i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d c9   = _mm512_set1_pd(0.125);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0f);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(23.0f/60.0f);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.055555555555555555555555555556);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0f/60.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(-1343.0f/150.0f);
                        t0                         = _mm512_fmadd_pd(half,cost,_1);  // 0.5*cost+1
                        const register __m512d c5   = _mm512_set1_pd(3769.0f/280.0f);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmsub_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0f/63.0f);
                        t2                         = _mm512_add_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c9,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2f);
                        t4                         = _mm512_mul_pd(c7,_mm512_sub_pd(_4,cost));
                        t5                         = _mm512_mul_pd(c8,_mm512_add_pd(_2,cost));
                        t5                         = _mm512_mul_pd(t5,k0a2);
                        t6                         = _mm512_mul_pd(k0a6,_mm512_add_pd(t4,t5)); // imaginary part
                        _mm512_store_pd(&S2i[0],   t6);
                        t4                         = _mm512_add_pd(t0,_mm512_add_pd(t2,t3));
                        _mm512_store_pd(&S2r[0],    _mm512_mul_pd(k0a3,t4);
                }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3214_zmm146r4_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict  S2r,
                                             double * __restrict  S2i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d c9   = _mm512_set1_pd(0.125);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0f);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(23.0f/60.0f);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.055555555555555555555555555556);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0f/60.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(-1343.0f/150.0f);
                        t0                         = _mm512_fmadd_pd(half,cost,_1);  // 0.5*cost+1
                        const register __m512d c5   = _mm512_set1_pd(3769.0f/280.0f);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmsub_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0f/63.0f);
                        t2                         = _mm512_add_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c9,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2f);
                        t4                         = _mm512_mul_pd(c7,_mm512_sub_pd(_4,cost));
                        t5                         = _mm512_mul_pd(c8,_mm512_add_pd(_2,cost));
                        t5                         = _mm512_mul_pd(t5,k0a2);
                        t6                         = _mm512_mul_pd(k0a6,_mm512_add_pd(t4,t5)); // imaginary part
                        _mm512_storeu_pd(&S2i[0],   t6);
                        t4                         = _mm512_add_pd(t0,_mm512_add_pd(t2,t3));
                        _mm512_storeu_pd(&S2r[0],    _mm512_mul_pd(k0a3,t4);
                }


                  
                  /*
                       Formula 3.2-16, optics contribution at upper end of resonance region
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3216_zmm8r8(const __m512d k0a,
                                         const __m512d tht,
                                         __m512d * __restrict S1r,
                                         __m512d * __restrict S1i) {

                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        const register __m512d _7   = _mm512_set1_pd(7.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        register __m512d k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512d sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512d t0r,t0i,t1r,t1i,t2,t3,cos6t;
                        k0ah  = _mm512_mul_pd(half,k0a);
                        htht  = _mm512_mul_pd(half,tht);
                        cost  = xcos(tht);
                        k0aa  = _mm512_add_pd(k0a,k0a);
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        sint  = xsin(htht);
                        sin2t = _mm512_mul_pd(sint,sint);
                        carr  = _mm512_set1_pd(-0.0f);
                        cari  = _mm512_mul_pd(k0aa,htht);
                        cexp_zmm8r8(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_pd(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_pd(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_pd(cost,cost);
                        cos3t = _mm512_mul_pd(cost,cos2t);
                        cos6t = _mm512_mul_pd(cos3t,cos2t);
                        t3    = _mm512_div_pd(sin2t,cos6t);
                        t0r   = _mm512_setzero_pd();
                        t0i   = _mm512_div_pd(Ii,_mm512_mul_pd(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_pd(_1,t0i); // second term
                        t2    = _mm512_div_pd(_7,_mm512_mul_pd(_4,k0a2));
                        cerr  = _mm512_mul_pd(t2,t3);
                        t1r   = _mm512_sub_pd(t0r,cerr);
                        t1i   = _mm512_sub_pd(t0i,cerr);
                        cmul_zmm8r8(cexpr,cexpi,t1r,t1i,&S1r,&S1i);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3216_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           double * __restrict __ATTR_ALIGN__(64) S1r,
                                           double * __restrict __ATTR_ALIGN__(64) S1i) {
     
                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d tht  = _mm512_load_pd(&ptht[0]);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        const register __m512d _7   = _mm512_set1_pd(7.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        register __m512d k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512d sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512d t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_pd(half,k0a);
                        htht  = _mm512_mul_pd(half,tht);
                        cost  = xcos(tht);
                        k0aa  = _mm512_add_pd(k0a,k0a);
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        sint  = xsin(htht);
                        sin2t = _mm512_mul_pd(sint,sint);
                        carr  = _mm512_set1_pd(-0.0f);
                        cari  = _mm512_mul_pd(k0aa,htht);
                        cexp_zmm8r8(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_pd(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_pd(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_pd(cost,cost);
                        cos3t = _mm512_mul_pd(cost,cos2t);
                        cos6t = _mm512_mul_pd(cos3t,cos2t);
                        t3    = _mm512_div_pd(sin2t,cos6t);
                        t0r   = _mm512_setzero_pd();
                        t0i   = _mm512_div_pd(Ii,_mm512_mul_pd(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_pd(_1,t0i); // second term
                        t2    = _mm512_div_pd(_7,_mm512_mul_pd(_4,k0a2));
                        cerr  = _mm512_mul_pd(t2,t3);
                        t1r   = _mm512_sub_pd(t0r,cerr);
                        t1i   = _mm512_sub_pd(t0i,cerr);
                        cmul_zmm8r8(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_store_pd(&S1r[0], resr);
                        _mm512_store_pd(&S1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3216_zmm8r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  ptht,
                                           double * __restrict  S1r,
                                           double * __restrict  S1i) {
     
                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d _4   = _mm512_set1_pd(4.0f);
                        const register __m512d _7   = _mm512_set1_pd(7.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        register __m512d k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512d sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512d t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_pd(half,k0a);
                        htht  = _mm512_mul_pd(half,tht);
                        cost  = xcos(tht);
                        k0aa  = _mm512_add_pd(k0a,k0a);
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        sint  = xsin(htht);
                        sin2t = _mm512_mul_pd(sint,sint);
                        carr  = _mm512_set1_pd(-0.0f);
                        cari  = _mm512_mul_pd(k0aa,htht);
                        cexp_zmm8r8(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_pd(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_pd(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_pd(cost,cost);
                        cos3t = _mm512_mul_pd(cost,cos2t);
                        cos6t = _mm512_mul_pd(cos3t,cos2t);
                        t3    = _mm512_div_pd(sin2t,cos6t);
                        t0r   = _mm512_setzero_pd();
                        t0i   = _mm512_div_pd(Ii,_mm512_mul_pd(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_pd(_1,t0i); // second term
                        t2    = _mm512_div_pd(_7,_mm512_mul_pd(_4,k0a2));
                        cerr  = _mm512_mul_pd(t2,t3);
                        t1r   = _mm512_sub_pd(t0r,cerr);
                        t1i   = _mm512_sub_pd(t0i,cerr);
                        cmul_zmm8r8(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_storeu_pd(&S1r[0], resr);
                        _mm512_storeu_pd(&S1i[0], resi);
                }


                   /*
                       Formula 3.2-17, optics contribution at upper end of resonance region
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3217_zmm8r8(const __m512d k0a,
                                         const __m512d tht,
                                         __m512d * __restrict S2r,
                                         __m512d * __restrict S2i) {

                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d quat = _mm512_set1_pd(0.25);
                        const register __m512d _6   = _mm512_set1_pd(6.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        register __m512d k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512d sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512d t0r,t0i,t1r,t1i,t2,t3,cos6t;
                        k0ah  = _mm512_mul_pd(half,k0a);
                        htht  = _mm512_mul_pd(half,tht);
                        cost  = xcos(tht);
                        k0aa  = _mm512_add_pd(k0a,k0a);
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        sint  = xsin(htht);
                        sin2t = _mm512_mul_pd(sint,sint);
                        carr  = _mm512_set1_pd(-0.0f);
                        cari  = _mm512_mul_pd(k0aa,htht);
                        cexp_zmm8r8(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_pd(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_pd(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_pd(cost,cost);
                        cos3t = _mm512_mul_pd(cost,cos2t);
                        cos6t = _mm512_mul_pd(cos3t,cos2t);
                        t3    = _mm512_div_pd(_mm512_mul_pd(_mm512_add_pd(_6,cost),sin2t),cos6t);
                        t0r   = _mm512_setzero_pd();
                        t0i   = _mm512_div_pd(_mm512_mul_pd(Ii,cost),_mm512_mul_pd(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_pd(_1,t0i); // second term
                        t2    = _mm512_mul_pd(quat,k0a2);
                        cerr  = _mm512_mul_pd(t2,t3);
                        t1r   = _mm512_sub_pd(t0r,cerr);
                        t1i   = _mm512_sub_pd(t0i,cerr);
                        cmul_zmm8r8(cexpr,cexpi,t1r,t1i,&S2r,&S2i);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3217_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           double * __restrict __ATTR_ALIGN__(64) S2r,
                                           double * __restrict __ATTR_ALIGN__(64) S2i) {
                        
                        register __m512d k0a        = _mm512_load_pd(&pk0a[0]);
                        register __m512d tht        = _mm512_load_pd(&tht[0]);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d quat = _mm512_set1_pd(0.25);
                        const register __m512d _6   = _mm512_set1_pd(6.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        register __m512d k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512d sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512d t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_pd(half,k0a);
                        htht  = _mm512_mul_pd(half,tht);
                        cost  = xcos(tht);
                        k0aa  = _mm512_add_pd(k0a,k0a);
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        sint  = xsin(htht);
                        sin2t = _mm512_mul_pd(sint,sint);
                        carr  = _mm512_set1_pd(-0.0f);
                        cari  = _mm512_mul_pd(k0aa,htht);
                        cexp_zmm8r8(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_pd(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_pd(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_pd(cost,cost);
                        cos3t = _mm512_mul_pd(cost,cos2t);
                        cos6t = _mm512_mul_pd(cos3t,cos2t);
                        t3    = _mm512_div_pd(_mm512_mul_pd(_mm512_add_pd(_6,cost),sin2t),cos6t);
                        t0r   = _mm512_setzero_pd();
                        t0i   = _mm512_div_pd(_mm512_mul_pd(Ii,cost),_mm512_mul_pd(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_pd(_1,t0i); // second term
                        t2    = _mm512_mul_pd(quat,k0a2);
                        cerr  = _mm512_mul_pd(t2,t3);
                        t1r   = _mm512_sub_pd(t0r,cerr);
                        t1i   = _mm512_sub_pd(t0i,cerr);
                        cmul_zmm8r8(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_store_pd(&S2r[0],resr);
                        _mm512_store_pd(&S2i[0],resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3217_zmm8r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  ptht,
                                           double * __restrict  S2r,
                                           double * __restrict  S2i) {
                        
                        register __m512d k0a        = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d tht        = _mm512_loadu_pd(&tht[0]);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        const register __m512d quat = _mm512_set1_pd(0.25);
                        const register __m512d _6   = _mm512_set1_pd(6.0f);
                        const register __m512d half = _mm512_set1_pd(0.5f);
                        register __m512d k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512d sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512d t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_pd(half,k0a);
                        htht  = _mm512_mul_pd(half,tht);
                        cost  = xcos(tht);
                        k0aa  = _mm512_add_pd(k0a,k0a);
                        k0a2  = _mm512_mul_pd(k0a,k0a);
                        sint  = xsin(htht);
                        sin2t = _mm512_mul_pd(sint,sint);
                        carr  = _mm512_set1_pd(-0.0f);
                        cari  = _mm512_mul_pd(k0aa,htht);
                        cexp_zmm8r8(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_pd(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_pd(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_pd(cost,cost);
                        cos3t = _mm512_mul_pd(cost,cos2t);
                        cos6t = _mm512_mul_pd(cos3t,cos2t);
                        t3    = _mm512_div_pd(_mm512_mul_pd(_mm512_add_pd(_6,cost),sin2t),cos6t);
                        t0r   = _mm512_setzero_pd();
                        t0i   = _mm512_div_pd(_mm512_mul_pd(Ii,cost),_mm512_mul_pd(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_pd(_1,t0i); // second term
                        t2    = _mm512_mul_pd(quat,k0a2);
                        cerr  = _mm512_mul_pd(t2,t3);
                        t1r   = _mm512_sub_pd(t0r,cerr);
                        t1i   = _mm512_sub_pd(t0i,cerr);
                        cmul_zmm8r8(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_storeu_pd(&S2r[0],resr);
                        _mm512_storeu_pd(&S2i[0],resi);
                }


                 /*
                        Low frequency region (k0a < 0.4) i.e. Rayleigh-region complex scattering
                        amplitudes.
                        Formula 3.2-20
                    */
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d S1_f3220_zmm8r8(const __m512d k0a,
                                           const __m512d tht) {

                          register const __m512d half = _mm512_set1_pd(0.5f);
                          register __m512d S1;
                          register __m512d k0a3,cost;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cost = xcos(tht);
                          S1   = _mm512_mul_pd(k0a3,_mm512_add_pd(half,cost));
                          return (S1);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3220_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           double * __restrict __ATTR_ALIGN__(64) S1) {

                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register const __m512d half = _mm512_set1_pd(0.5f);
                          register __m512d k0a3,cost;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cost = xcos(tht);
                          _mm512_store_pd(&S1[0],_mm512_mul_pd(k0a3,_mm512_add_pd(half,cost)));
               
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3220_zmm8r8_u(const double * __restrict pk0a,
                                           const double * __restrict ptht,
                                           double * __restrict S1) {

                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register const __m512d half = _mm512_set1_pd(0.5f);
                          register __m512d k0a3,cost;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cost = xcos(tht);
                          _mm512_storeu_pd(&S1[0],_mm512_mul_pd(k0a3,_mm512_add_pd(half,cost)));
               
               }


                   /*
                        Low frequency region (k0a < 0.4) i.e. Rayleigh-region complex scattering
                        amplitudes.
                        Formula 3.2-21
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d S2_f3221_zmm8r8(const __m512d k0a,
                                           const __m512d tht) {

                          register const __m512d half = _mm512_set1_pd(0.5f);
                          register const __m512d _1   = _mm512_set1_pd(1.0f);
                          register __m512d S2;
                          register __m512d k0a3,cost,t0; 
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cost  = xcos(tht);
                          t0    = _mm512_fmadd_pd(half,cost,_1);
                          S2    = _mm512_mul_pd(k0a3,t0);
                          return (S2);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3221_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           double * __restrict S2) {

                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d tht = _mm512_load_pd(&ptht[0]);
                          register const __m512d half = _mm512_set1_pd(0.5f);
                          register const __m512d _1   = _mm512_set1_pd(1.0f);
                          register __m512d k0a3,cost,t0; 
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cost  = xcos(tht);
                          t0    = _mm512_fmadd_pd(half,cost,_1);
                          _mm512_store_pd(&S2[0], _mm512_mul_pd(k0a3,t0));
                          
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3221_zmm8r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  ptht,
                                           double * __restrict S2) {

                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          register const __m512d half = _mm512_set1_pd(0.5f);
                          register const __m512d _1   = _mm512_set1_pd(1.0f);
                          register __m512d k0a3,cost,t0; 
                          k0a3  = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          cost  = xcos(tht);
                          t0    = _mm512_fmadd_pd(half,cost,_1);
                          _mm512_storeu_pd(&S2[0], _mm512_mul_pd(k0a3,t0));
                          
               }


                 /*
                         The E-plane and H-plane scattering cross sections.
                         Formula 3.2-22
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3222_zmm8r8(const __m512d k0a,
                                            const __m512d a,
                                            const __m512d theta) {

                           const register __m512d _4   = _mm512_set1_pd(4.0f);
                           const register __m512d half = _mm512_set1_pd(0.5f);
                           register __m512d k0a4,a2,sqr,cost,t0,t1,t2; 
                           register __m512d rcs;
                           cost = xcos(theta);
                           k0a4 = _mm512_mul_pd(k0a,k0a,_mm512_mul_pd(k0a,k0a));
                           a2   = _mm512_mul_pd(a,a);
                           t0   = _mm512_add_pd(half,cost);
                           sqr  = _mm512_mul_pd(t0,t0);
                           t1   = _mm512_mul_pd(pi,a2);
                           t2   = _mm512_mul_pd(_4,k0a4);
                           rcs  = _mm512_mul_pd(t1,_mm512_mul_pd(t2,sqr));
                           return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3222_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ptheta,
                                              double * __restrict __ATTR_ALIGN__(64) rcs ) {

                           register __m512d  k0a       = _mm512_load_pd(&pk0a[0]);
                           const register __m512d _4   = _mm512_set1_pd(4.0f);
                           register __m512d  a         = _mm512_load_pd(&pa[0]);
                           const register __m512d half = _mm512_set1_pd(0.5f);
                           register __m512d theta      = _mm512_load_pd(&ptheta[0]);
                           register __m512d k0a4,a2,sqr,cost,t0,t1,t2; 
                           register __m512d rcs;
                           cost = xcos(theta);
                           k0a4 = _mm512_mul_pd(k0a,k0a,_mm512_mul_pd(k0a,k0a));
                           a2   = _mm512_mul_pd(a,a);
                           t0   = _mm512_add_pd(half,cost);
                           sqr  = _mm512_mul_pd(t0,t0);
                           t1   = _mm512_mul_pd(pi,a2);
                           t2   = _mm512_mul_pd(_4,k0a4);
                           _mm512_store_pd(&rcs[0], _mm512_mul_pd(t1,_mm512_mul_pd(t2,sqr)));
              }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3222_zmm8r8_u(  const double * __restrict  pk0a,
                                              const double * __restrict  pa,
                                              const double * __restrict ptheta,
                                              double * __restrict  rcs ) {

                           register __m512d  k0a       = _mm512_loadu_pd(&pk0a[0]);
                           const register __m512d _4   = _mm512_set1_pd(4.0f);
                           register __m512d  a         = _mm512_loadu_pd(&pa[0]);
                           const register __m512d half = _mm512_set1_pd(0.5f);
                           register __m512d theta      = _mm512_loadu_pd(&ptheta[0]);
                           register __m512d k0a4,a2,sqr,cost,t0,t1,t2; 
                           register __m512d rcs;
                           cost = xcos(theta);
                           k0a4 = _mm512_mul_pd(k0a,k0a,_mm512_mul_pd(k0a,k0a));
                           a2   = _mm512_mul_pd(a,a);
                           t0   = _mm512_add_pd(half,cost);
                           sqr  = _mm512_mul_pd(t0,t0);
                           t1   = _mm512_mul_pd(pi,a2);
                           t2   = _mm512_mul_pd(_4,k0a4);
                           _mm512_storeu_pd(&rcs[0], _mm512_mul_pd(t1,_mm512_mul_pd(t2,sqr)));
              }


                 /*
                         The E-plane and H-plane scattering cross sections.
                         Formula 3.2-23
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3223_zmm8r8(const __m512d k0a,
                                            const __m512d a,
                                            const __m512d theta) {

                           const register __m512d _4   = _mm512_set1_pd(4.0f);
                           const register __m512d half = _mm512_set1_pd(0.5f);
                           const register __m512d _1   = _mm512_set1_pd(1.0f);
                           register __m512d k0a4,a2,sqr,cost,t0,t1,t2; 
                           register __m512d rcs;
                           cost = xcos(theta);
                           k0a4 = _mm512_mul_pd(k0a,k0a,_mm512_mul_pd(k0a,k0a));
                           a2   = _mm512_mul_pd(a,a);
                           t0   = _mm512_fmadd_pd(half,cost,_1);
                           sqr  = _mm512_mul_pd(t0,t0);
                           t1   = _mm512_mul_pd(pi,a2);
                           t2   = _mm512_mul_pd(_4,k0a4);
                           rcs  = _mm512_mul_pd(t1,_mm512_mul_pd(t2,sqr));
                           return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3223_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) ptheta,
                                              double * __restrict __ATTR_ALIGN__(64) rcs  ) {
       
                           register __m512d       k0a  = _mm512_load_pd(&pk0a[0]);
                           const register __m512d _4   = _mm512_set1_pd(4.0f);
                           register __m512d       a    = _mm512_load_pd(&pa[0]);
                           const register __m512d half = _mm512_set1_pd(0.5f);
                           register __m512d       theta= _mm512_load_pd(&ptheta[0]);
                           const register __m512d _1   = _mm512_set1_pd(1.0f);
                           register __m512d k0a4,a2,sqr,cost,t0,t1,t2; 
                           register __m512d rcs;
                           cost = xcos(theta);
                           k0a4 = _mm512_mul_pd(k0a,k0a,_mm512_mul_pd(k0a,k0a));
                           a2   = _mm512_mul_pd(a,a);
                           t0   = _mm512_fmadd_pd(half,cost,_1);
                           sqr  = _mm512_mul_pd(t0,t0);
                           t1   = _mm512_mul_pd(pi,a2);
                           t2   = _mm512_mul_pd(_4,k0a4);
                           _mm512_store_pd(&rcs[0] ,_mm512_mul_pd(t1,_mm512_mul_pd(t2,sqr)));
                           
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3223_zmm8r8_u(  const double * __restrict  pk0a,
                                              const double * __restrict  pa,
                                              const double * __restrict  ptheta,
                                              double * __restrict  rcs  ) {
       
                           register __m512d       k0a  = _mm512_loadu_pd(&pk0a[0]);
                           const register __m512d _4   = _mm512_set1_pd(4.0f);
                           register __m512d       a    = _mm512_loadu_pd(&pa[0]);
                           const register __m512d half = _mm512_set1_pd(0.5f);
                           register __m512d       theta= _mm512_loadu_pd(&ptheta[0]);
                           const register __m512d _1   = _mm512_set1_pd(1.0f);
                           register __m512d k0a4,a2,sqr,cost,t0,t1,t2; 
                           register __m512d rcs;
                           cost = xcos(theta);
                           k0a4 = _mm512_mul_pd(k0a,k0a,_mm512_mul_pd(k0a,k0a));
                           a2   = _mm512_mul_pd(a,a);
                           t0   = _mm512_fmadd_pd(half,cost,_1);
                           sqr  = _mm512_mul_pd(t0,t0);
                           t1   = _mm512_mul_pd(pi,a2);
                           t2   = _mm512_mul_pd(_4,k0a4);
                           _mm512_storeu_pd(&rcs[0] ,_mm512_mul_pd(t1,_mm512_mul_pd(t2,sqr)));
                           
               }


                 /*
                        High frequency region (k0a > 20).
                        Complex scattering amplitudes.
                        Formula 3.2-24
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S12_f3224_zmm8r8( const __m512d k0a,
                                           const __m512d tht,
                                           __m512d * __restrict S12r,
                                           __m512d * __restrict S12i) {

                        const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                        const register __m512d htht = _mm512_mul_pd(_mm512_set1_pd(0.5f),tht);
                        register __m512d cosht,hk0a,_2k0a,carr,cari,cexr,cexi;
                        cosht = xcos(htht);          
                        hk0a  = _mm512_mul_pd(nhlf,k0a);
                        _2k0a = _mm512_add_pd(k0a,k0a);
                        carr  = nIr;
                        cari  = _mm512_mul_pd(nIi,_mm512_mul_pd(_2k0a,cosht));
                        cexp_zmm8r8(carr,cari,&cexr,cexi);
                        *S12r  = _mm512_mul_pd(hk0a,cexr);
                        *S12i  = _mm512_mul_pd(hk0a,cexi);
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S12_f3224_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht,
                                             double * __restrict __ATTR_ALIGN__(64) S12r,
                                             double * __restrict __ATTR_ALIGN__(64) S12i) {

                        const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                        const register __m512d tht  = _mm512_load_pd(&ptht[0]);
                        const register __m512d htht = _mm512_mul_pd(_mm512_set1_pd(0.5f),tht);
                        register __m512d cosht,hk0a,_2k0a,carr,cari,cexr,cexi;
                        cosht = xcos(htht);          
                        hk0a  = _mm512_mul_pd(nhlf,k0a);
                        _2k0a = _mm512_add_pd(k0a,k0a);
                        carr  = nIr;
                        cari  = _mm512_mul_pd(nIi,_mm512_mul_pd(_2k0a,cosht));
                        cexp_zmm8r8(carr,cari,&cexr,cexi);
                        _mm512_store_pd(&S12r[0], _mm512_mul_pd(hk0a,cexr));
                        _mm512_store_pd(&S12i[0], _mm512_mul_pd(hk0a,cexi));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S12_f3224_zmm8r8_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict  S12r,
                                             double * __restrict  S12i) {

                        const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                        const register __m512d tht  = _mm512_loadu_pd(&ptht[0]);
                        const register __m512d htht = _mm512_mul_pd(_mm512_set1_pd(0.5f),tht);
                        register __m512d cosht,hk0a,_2k0a,carr,cari,cexr,cexi;
                        cosht = xcos(htht);          
                        hk0a  = _mm512_mul_pd(nhlf,k0a);
                        _2k0a = _mm512_add_pd(k0a,k0a);
                        carr  = nIr;
                        cari  = _mm512_mul_pd(nIi,_mm512_mul_pd(_2k0a,cosht));
                        cexp_zmm8r8(carr,cari,&cexr,cexi);
                        _mm512_storeu_pd(&S12r[0], _mm512_mul_pd(hk0a,cexr));
                        _mm512_storeu_pd(&S12i[0], _mm512_mul_pd(hk0a,cexi));
              }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3225_zmm8r8(const __m512d a) {
                         
                         register __m512d rcs;
                         rcs = _mm512_mul_pd(pi,_mm512_mul_pd(a,a));
                         return (rcs);
              }


                /*
                       Resonance region (0.4 < k0a < 20.0), equations are valid only for k0a < 1.0.
                       Complex scattering amplitude represented as a scattering function -- formula 3.2-26
                 */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S12_f3226_zmm8r8(const __m512d k0a,
                                          __m512d * __restrict S12r,
                                          __m512d * __restrict S12i) {

                        __m512d k0a2,k0a3,k0a4,k0a6;
                        __m512d t0,t1,t2,t3;
                        const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c0   = _mm512_set1_pd(113.0f/90.0f);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c1   = _mm512_set1_pd(1783.0f/2100.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c2   = _mm512_set1_pd(670057.0f/396900.0f);
                        t0                         = _mm512_mul_pd(nhlf,k0a3)
                        const register __m512d c3   = _mm512_set1_pd(0.833333333333333333333333333333);
                        const register __m512d c4   = _mm512_set1_pd(0.24);
                        t1                         = _mm512_fmsub_pd(_mm512_add_pd(_1,c0),k0a2,
                                                                 _mm512_fmsub_pd(c1,k0a4,
                                                                             _mm512_mul_pd(c2,k0a6)));
                        *S12r                      = _mm512_mul_pd(t0,t1);
                        t2                         = _mm512_mul_pd(Ini,_mm512_mul_pd(c3,k0a6));
                        t3                         = _mm512_fmadd_pd(c4,k0a3,_1);
                        *S12i                      = _mm512_mul_pd(t2,t3);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S12_f3226_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            double * __restrict __ATTR_ALIGN__(64) S12r,
                                            double * __restrict __ATTR_ALIGN__(64) S12i) {

                        
                        __m512d k0a2,k0a3,k0a4,k0a6;
                        __m512d t0,t1,t2,t3;
                        register __m512d       k0a  = _mm512_load_pd(&pk0a[0]);
                        const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c0   = _mm512_set1_pd(113.0f/90.0f);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c1   = _mm512_set1_pd(1783.0f/2100.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c2   = _mm512_set1_pd(670057.0f/396900.0f);
                        t0                         = _mm512_mul_pd(nhlf,k0a3)
                        const register __m512d c3   = _mm512_set1_pd(0.833333333333333333333333333333);
                        const register __m512d c4   = _mm512_set1_pd(0.24);
                        t1                         = _mm512_fmsub_pd(_mm512_add_pd(_1,c0),k0a2,
                                                                 _mm512_fmsub_pd(c1,k0a4,
                                                                             _mm512_mul_pd(c2,k0a6)));
                        _mm512_store_pd(&S12r[0],  _mm512_mul_pd(t0,t1));
                        t2                         = _mm512_mul_pd(Ini,_mm512_mul_pd(c3,k0a6));
                        t3                         = _mm512_fmadd_pd(c4,k0a3,_1);
                        _mm512_store_pd(&S12i[0],  _mm512_mul_pd(t2,t3));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S12_f3226_zmm8r8_u(const double * __restrict  pk0a,
                                            double * __restrict  S12r,
                                            double * __restrict  S12i) {

                        
                        __m512d k0a2,k0a3,k0a4,k0a6;
                        __m512d t0,t1,t2,t3;
                        register __m512d       k0a  = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c0   = _mm512_set1_pd(113.0f/90.0f);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c1   = _mm512_set1_pd(1783.0f/2100.0f);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c2   = _mm512_set1_pd(670057.0f/396900.0f);
                        t0                         = _mm512_mul_pd(nhlf,k0a3)
                        const register __m512d c3   = _mm512_set1_pd(0.833333333333333333333333333333);
                        const register __m512d c4   = _mm512_set1_pd(0.24);
                        t1                         = _mm512_fmsub_pd(_mm512_add_pd(_1,c0),k0a2,
                                                                 _mm512_fmsub_pd(c1,k0a4,
                                                                             _mm512_mul_pd(c2,k0a6)));
                        _mm512_storeu_pd(&S12r[0],  _mm512_mul_pd(t0,t1));
                        t2                         = _mm512_mul_pd(Ini,_mm512_mul_pd(c3,k0a6));
                        t3                         = _mm512_fmadd_pd(c4,k0a3,_1);
                        _mm512_storeu_pd(&S12i[0],  _mm512_mul_pd(t2,t3));
               }


                 /*
                           Resonance region (0.4 < k0a < 20.0), equations are valid only for k0a < 1.0.
                           Radar cross-section, formula 3.2-27
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3227_zmm8r8(const __m512d k0a,
                                            const __m512d a) {

                          register __m512d a2,k0a2,k0a4,k0a6,t0;
                          register __m512d rcs,t1;
                          a2                       = _mm512_mul_pd(a,a);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          k0a2                     = _mm512_mul_pd(k0a,k0a);
                          const register __m512d c0 = _mm512_set1_pd(113.0f/45.0f);
                          k0a4                     = _mm512_mul_pd(k0a2,k0a2);
                          const register __m512d c1 = _mm512_set1_pd(6899.0f/56700.0f);
                          k0a6                     = _mm512_mul_pd(k0a4,_mm512_mul_pd(k0a,k0a));
                          const register __m512d c2 = _mm512_set1_pd(5419129.0f/1984500.0f);
                          t0                       = _mm512_mul_pd(pi,_mm512_mul_pd(a2,k0a4));
                          t1                       = _mm512_fmsub_pd(_mm512_add_pd(_1,c0),k0a2,
                                                                             _mm512_fmsub_pd(c1,k0a4,
                                                                                         _mm512_mul_pd(c2,k0a6)));
                          rcs                      = _mm512_mul_pd(t0,t1);
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3227_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            double * __restrict __ATTR_ALIGN__(64) rcs) {

                          register __m512d k0a      = _mm512_load_pd(&pk0a[0]);
                          register __m512d  a       = _mm512_load_pd(&pa[0]);
                          register __m512d a2,k0a2,k0a4,k0a6,t0;
                          register __m512d t1;
                          a2                       = _mm512_mul_pd(a,a);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          k0a2                     = _mm512_mul_pd(k0a,k0a);
                          const register __m512d c0 = _mm512_set1_pd(113.0f/45.0f);
                          k0a4                     = _mm512_mul_pd(k0a2,k0a2);
                          const register __m512d c1 = _mm512_set1_pd(6899.0f/56700.0f);
                          k0a6                     = _mm512_mul_pd(k0a4,_mm512_mul_pd(k0a,k0a));
                          const register __m512d c2 = _mm512_set1_pd(5419129.0f/1984500.0f);
                          t0                       = _mm512_mul_pd(pi,_mm512_mul_pd(a2,k0a4));
                          t1                       = _mm512_fmsub_pd(_mm512_add_pd(_1,c0),k0a2,
                                                                             _mm512_fmsub_pd(c1,k0a4,
                                                                                         _mm512_mul_pd(c2,k0a6)));
                          _mm512_store_pd(&rcs[0],  _mm512_mul_pd(t0,t1));
                         
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3227_zmm8r8_u(const double * __restrict  pk0a,
                                            const double * __restrict  pa,
                                            double * __restrict  rcs) {

                          register __m512d k0a      = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d  a       = _mm512_loadu_pd(&pa[0]);
                          register __m512d a2,k0a2,k0a4,k0a6,t0;
                          register __m512d t1;
                          a2                       = _mm512_mul_pd(a,a);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          k0a2                     = _mm512_mul_pd(k0a,k0a);
                          const register __m512d c0 = _mm512_set1_pd(113.0f/45.0f);
                          k0a4                     = _mm512_mul_pd(k0a2,k0a2);
                          const register __m512d c1 = _mm512_set1_pd(6899.0f/56700.0f);
                          k0a6                     = _mm512_mul_pd(k0a4,_mm512_mul_pd(k0a,k0a));
                          const register __m512d c2 = _mm512_set1_pd(5419129.0f/1984500.0f);
                          t0                       = _mm512_mul_pd(pi,_mm512_mul_pd(a2,k0a4));
                          t1                       = _mm512_fmsub_pd(_mm512_add_pd(_1,c0),k0a2,
                                                                             _mm512_fmsub_pd(c1,k0a4,
                                                                                         _mm512_mul_pd(c2,k0a6)));
                          _mm512_storeu_pd(&rcs[0],  _mm512_mul_pd(t0,t1));
                         
               }


                 /*
                       Scattering functions equation at upper end of resonance region.
                       Optics wave term and creeping wave term.
                       Optics wave term, formula 3.2-28
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void FO_f3228_zmm8r8(const __m512d k0a,
                                         __m512d * __restrict FOr,
                                         __m512d * __restrict FOi) {

                        const register __m512d hlf = _mm512_set1_pd(0.5f);
                        const register __m512d c0  = _mm512_set1_pd(0.916666666666666666666666666667);
                        register __m512d k0a2,t0;
                        k0a2                      = _mm512_mul_pd(k0a,k0a);
                        t0                        = _mm512_sub_pd(k0a2,c0);
                        *FOr                      = Inr;
                        *FOi                      = _mm512_mul_pd(Ini,_mm512_mul_pd(hlf,t0));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void FO_f3228_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           double * __restrict __ATTR_ALIGN__(64) FOr,
                                           double * __restrict __ATTR_ALIGN__(64) FOi) {

                        register __m512d       k0a = _mm512_load_pd(&pk0a[0]);
                        const register __m512d hlf = _mm512_set1_pd(0.5f);
                        const register __m512d c0  = _mm512_set1_pd(0.916666666666666666666666666667);
                        register __m512d k0a2,t0;
                        k0a2                      = _mm512_mul_pd(k0a,k0a);
                        t0                        = _mm512_sub_pd(k0a2,c0);
                        _mm512_store_pd(&FOr[0],   nIr);
                        _mm512_store_pd(&FOi[0],   _mm512_mul_pd(nIi,_mm512_mul_pd(hlf,t0)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void FO_f3228_zmm8r8_u(const double * __restrict pk0a,
                                           double * __restrict FOr,
                                           double * __restrict FOi) {

                        register __m512d       k0a = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d hlf = _mm512_set1_pd(0.5f);
                        const register __m512d c0  = _mm512_set1_pd(0.916666666666666666666666666667);
                        register __m512d k0a2,t0;
                        k0a2                      = _mm512_mul_pd(k0a,k0a);
                        t0                        = _mm512_sub_pd(k0a2,c0);
                        _mm512_storeu_pd(&FOr[0],   Inr);
                        _mm512_storeu_pd(&FOi[0],   _mm512_mul_pd(Ini,_mm512_mul_pd(hlf,t0)));
                }


                    /*
                       Scattering functions equation at upper end of resonance region.
                       Optics wave term and creeping wave term.
                       Creeping wave term, formula 3.2-29
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void FC_f3229_zmm8r8(const __m512d x,
                                         __m512d * __restrict FCr,
                                         __m512d * __restrict FCi) {
                         
                         const __m512d   x16        = _mm512_add_pd(x,
                                                               _mm512_set1_pd(0.666666666666666666666666666667));
                         const __m512d          _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
                         const register __m512d c0  = _mm512_set1_pd(1.357588);
                         const __m512d          px43= _mm512_pow_pd(x,
                                                               _mm512_set1_pd(1.333333333333333333333333333333));
                         const register __m512d c1  = _mm512_set1_pd(0.807104);
                         const __m512d        pnx23 = _mm512_rcp14_pd(_mm512_pow_pd(x,
                                                               _mm512_set1_pd(0.666666666666666666666666666667)));
                         const register __m512d c0r = _mm512_set1_pd(0.032927);
                         const __m512d        pnx43 = _mm512_rcp14_pd(px43);
                         const register __m512d c0i = _mm512_set1_pd(0.057154);
                         const __m512d        px13  = _mm512_pow_pd(x,
                                                               _mm512_set1_pd(0.333333333333333333333333333333333));
                         const register __m512d c1r = _mm512_set1_pd(0.242679);
                         const register __m512d c1i = _mm512_set1_pd(-0.710672);
                         const __m512d        pnx13 = _mm512_rcp14_pd(px13);
                         const register __m512d c2r = _mm512_set1_pd(0.001846);
                         const __m512d        npx13 = _mm512_sub_pd(_mm512_setzero_pd(),px13);
                         const register __m512d c2i = _mm512_set1_pd(0.008027);
                         const register __m512d c3r = _mm512_set1_pd(0.741196);
                         const register __m512d c3i = _mm512_set1_pd(1.283788);
                         const register __m512d c4r = _mm512_set1_pd(4.400004);
                         const register __m512d c4i = _mm512_set1_pd(-2.540343);
                         const register __m512d c5r = _mm512_set1_pd(0.890792);
                         const register __m512d c5i = _mm512_set1_pd(0.514299);
                         const register __m512d c6r = _mm512_set1_pd(0.798821);
                         const register __m512d c6i = _mm512_set1_pd(1.383598);
                         const register __m512d c7r = _mm512_set1_pd(10.097912);
                         const register __m512d c7i = _mm512_set1_pd(-5.830032);
                         const register __m512d c8r = _mm512_set1_pd(0.624641);
                         const register __m512d c8i = _mm512_set1_pd(0.360637);
                         __m512d E1r,E1i,H1r,H1i,ix43r,ix43i,ear,eai,t3r,t3i;
                         __m512d t0r,t0i,t1r,t1i,t2r,t2i,cearr,ceari,tr4,t4i;
                         __m512d exp1r,exp1i,exp2r,exp2i,t5r,t5i;
                         ix43r = Inr;
                         ix43i = _mm512_mul_pd(Ini,px43);
                         ear = Ir;
                         eai = _mm512_add_pd(_mm512_mul_pd(Ii,_2pi),x16);
                         cexp_zmm8r8(ear,eai,&t3r,&t3i); 
                         cmul_zmm8r8(ix43r,ix43i,t3r,t3i,&cearr,&ceari);// ei2pi(x+1/6)
                         t0r = _mm512_fmadd_pd(pnx23,c1r,c0r);
                         t0i = _mm512_fmadd_pd(pnx23,c1i,c0i);
                         t1r = _mm512_sub_pd(t0r,_mm512_mul_pd(pnx43,c2r));
                         t1i = _mm512_sub_pd(t0i,_mm512_mul_pd(pnx43,c2i));
                         t2r = _mm512_mul_pd(ix43r,t1r); // first term re
                         t2i = _mm512_mul_pd(ix43i,t1i); // first term im
                         t0r = _mm512_fmadd_pd(xpn23,c3r,c0);
                         t0i = _mm512_fmadd_pd(xpn23,c3i,c0);
                         t4r = _mm512_fmadd_pd(npx13,c4r,
                                                    _mm512_mul_pd(pnx13,c5r));
                         t4i = _mm512_fmadd_pd(npx13,c4i,
                                                    _mm512_mul_pd(pnx13,c5i));
                         cexp_zmm8r8(t4r,t4i,&exp1r,&exp1i);
                         cmul_zmm8r8(t0r,t0i,exp1r,exp1i,&E1r,&E1i);//t3r,t3i hold second exp term
                         t1r = _mm512_fmadd_pd(pnx23,c6r,c1);
                         t1i = _mm512_fmadd_pd(pnx23,c6i,c1);
                         t5r = _mm512_fmsub_pd(npx13,c7r,
                                                     _mm512_mul_pd(pnx13,c8r));
                         t5i = _mm512_fmsub_pd(npx13,c7i,
                                                     _mm512_mul_pd(pnx13,c8i));
                         cexp_zmm8r8(t5r,t5i,&exp2r,&exp2i);
                         cmul_zmm8r8(t1r,t1i,exp2r,exp2i,&H1r,&H1i);
                         t3r = _mm512_add_pd(E1r,H1r);
                         t3i = _mm512_add_pd(E1i,H1i);
                         cmul_zmm8r8(cearr,ceari,t3r,t3i,&t1r,&t1i);
                         *FCr = _mm512_sub_pd(t2r,t1r);
                         *FCi = _mm512_sub_pd(t2i,t1i);
               } 


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void FC_f3229_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) px,
                                         double * __restrict __ATTR_ALIGN__(64) FCr,
                                         double * __restrict __ATTR_ALIGN__(64) FCi) {
                         
                         register __m512d x         = _mm512_load_pd(&px[0]);
                         const __m512d   x16        = _mm512_add_pd(x,
                                                               _mm512_set1_pd(0.666666666666666666666666666667));
                         const __m512d          _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
                         const register __m512d c0  = _mm512_set1_pd(1.357588);
                         const __m512d          px43= _mm512_pow_pd(x,
                                                               _mm512_set1_pd(1.333333333333333333333333333333));
                         const register __m512d c1  = _mm512_set1_pd(0.807104);
                         const __m512d        pnx23 = _mm512_rcp14_pd(_mm512_pow_pd(x,
                                                               _mm512_set1_pd(0.666666666666666666666666666667)));
                         const register __m512d c0r = _mm512_set1_pd(0.032927);
                         const __m512d        pnx43 = _mm512_rcp14_pd(px43);
                         const register __m512d c0i = _mm512_set1_pd(0.057154);
                         const __m512d        px13  = _mm512_pow_pd(x,
                                                               _mm512_set1_pd(0.333333333333333333333333333333333));
                         const register __m512d c1r = _mm512_set1_pd(0.242679);
                         const register __m512d c1i = _mm512_set1_pd(-0.710672);
                         const __m512d        pnx13 = _mm512_rcp14_pd(px13);
                         const register __m512d c2r = _mm512_set1_pd(0.001846);
                         const __m512d        npx13 = _mm512_sub_pd(_mm512_setzero_pd(),px13);
                         const register __m512d c2i = _mm512_set1_pd(0.008027);
                         const register __m512d c3r = _mm512_set1_pd(0.741196);
                         const register __m512d c3i = _mm512_set1_pd(1.283788);
                         const register __m512d c4r = _mm512_set1_pd(4.400004);
                         const register __m512d c4i = _mm512_set1_pd(-2.540343);
                         const register __m512d c5r = _mm512_set1_pd(0.890792);
                         const register __m512d c5i = _mm512_set1_pd(0.514299);
                         const register __m512d c6r = _mm512_set1_pd(0.798821);
                         const register __m512d c6i = _mm512_set1_pd(1.383598);
                         const register __m512d c7r = _mm512_set1_pd(10.097912);
                         const register __m512d c7i = _mm512_set1_pd(-5.830032);
                         const register __m512d c8r = _mm512_set1_pd(0.624641);
                         const register __m512d c8i = _mm512_set1_pd(0.360637);
                         __m512d E1r,E1i,H1r,H1i,ix43r,ix43i,ear,eai,t3r,t3i;
                         __m512d t0r,t0i,t1r,t1i,t2r,t2i,cearr,ceari,tr4,t4i;
                         __m512d exp1r,exp1i,exp2r,exp2i,t5r,t5i;
                         ix43r = Inr;
                         ix43i = _mm512_mul_pd(Ini,px43);
                         ear = Ir;
                         eai = _mm512_add_pd(_mm512_mul_pd(Ii,_2pi),x16);
                         cexp_zmm8r8(ear,eai,&t3r,&t3i); 
                         cmul_zmm8r8(ix43r,ix43i,t3r,t3i,&cearr,&ceari);// ei2pi(x+1/6)
                         t0r = _mm512_fmadd_pd(pnx23,c1r,c0r);
                         t0i = _mm512_fmadd_pd(pnx23,c1i,c0i);
                         t1r = _mm512_sub_pd(t0r,_mm512_mul_pd(pnx43,c2r));
                         t1i = _mm512_sub_pd(t0i,_mm512_mul_pd(pnx43,c2i));
                         t2r = _mm512_mul_pd(ix43r,t1r); // first term re
                         t2i = _mm512_mul_pd(ix43i,t1i); // first term im
                         t0r = _mm512_fmadd_pd(xpn23,c3r,c0);
                         t0i = _mm512_fmadd_pd(xpn23,c3i,c0);
                         t4r = _mm512_fmadd_pd(npx13,c4r,
                                                    _mm512_mul_pd(pnx13,c5r));
                         t4i = _mm512_fmadd_pd(npx13,c4i,
                                                    _mm512_mul_pd(pnx13,c5i));
                         cexp_zmm8r8(t4r,t4i,&exp1r,&exp1i);
                         cmul_zmm8r8(t0r,t0i,exp1r,exp1i,&E1r,&E1i);//t3r,t3i hold second exp term
                         t1r = _mm512_fmadd_pd(pnx23,c6r,c1);
                         t1i = _mm512_fmadd_pd(pnx23,c6i,c1);
                         t5r = _mm512_fmsub_pd(npx13,c7r,
                                                     _mm512_mul_pd(pnx13,c8r));
                         t5i = _mm512_fmsub_pd(npx13,c7i,
                                                     _mm512_mul_pd(pnx13,c8i));
                         cexp_zmm8r8(t5r,t5i,&exp2r,&exp2i);
                         cmul_zmm8r8(t1r,t1i,exp2r,exp2i,&H1r,&H1i);
                         t3r = _mm512_add_pd(E1r,H1r);
                         t3i = _mm512_add_pd(E1i,H1i);
                         cmul_zmm8r8(cearr,ceari,t3r,t3i,&t1r,&t1i);
                         _mm512_store_pd(&FCr[0] ,_mm512_sub_pd(t2r,t1r));
                         _mm512_store_pd(&FCi[0] ,_mm512_sub_pd(t2i,t1i));
               } 

                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void FC_f3229_zmm8r8_u(const double * __restrict  px,
                                           double * __restrict  FCr,
                                           double * __restrict  FCi) {
                         
                         register __m512d x         = _mm512_loadu_pd(&px[0]);
                         const __m512d   x16        = _mm512_add_pd(x,
                                                               _mm512_set1_pd(0.666666666666666666666666666667));
                         const __m512d          _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
                         const register __m512d c0  = _mm512_set1_pd(1.357588);
                         const __m512d          px43= _mm512_pow_pd(x,
                                                               _mm512_set1_pd(1.333333333333333333333333333333));
                         const register __m512d c1  = _mm512_set1_pd(0.807104);
                         const __m512d        pnx23 = _mm512_rcp14_pd(_mm512_pow_pd(x,
                                                               _mm512_set1_pd(0.666666666666666666666666666667)));
                         const register __m512d c0r = _mm512_set1_pd(0.032927);
                         const __m512d        pnx43 = _mm512_rcp14_pd(px43);
                         const register __m512d c0i = _mm512_set1_pd(0.057154);
                         const __m512d        px13  = _mm512_pow_pd(x,
                                                               _mm512_set1_pd(0.333333333333333333333333333333333));
                         const register __m512d c1r = _mm512_set1_pd(0.242679);
                         const register __m512d c1i = _mm512_set1_pd(-0.710672);
                         const __m512d        pnx13 = _mm512_rcp14_pd(px13);
                         const register __m512d c2r = _mm512_set1_pd(0.001846);
                         const __m512d        npx13 = _mm512_sub_pd(_mm512_setzero_pd(),px13);
                         const register __m512d c2i = _mm512_set1_pd(0.008027);
                         const register __m512d c3r = _mm512_set1_pd(0.741196);
                         const register __m512d c3i = _mm512_set1_pd(1.283788);
                         const register __m512d c4r = _mm512_set1_pd(4.400004);
                         const register __m512d c4i = _mm512_set1_pd(-2.540343);
                         const register __m512d c5r = _mm512_set1_pd(0.890792);
                         const register __m512d c5i = _mm512_set1_pd(0.514299);
                         const register __m512d c6r = _mm512_set1_pd(0.798821);
                         const register __m512d c6i = _mm512_set1_pd(1.383598);
                         const register __m512d c7r = _mm512_set1_pd(10.097912);
                         const register __m512d c7i = _mm512_set1_pd(-5.830032);
                         const register __m512d c8r = _mm512_set1_pd(0.624641);
                         const register __m512d c8i = _mm512_set1_pd(0.360637);
                         __m512d E1r,E1i,H1r,H1i,ix43r,ix43i,ear,eai,t3r,t3i;
                         __m512d t0r,t0i,t1r,t1i,t2r,t2i,cearr,ceari,tr4,t4i;
                         __m512d exp1r,exp1i,exp2r,exp2i,t5r,t5i;
                         ix43r = Inr;
                         ix43i = _mm512_mul_pd(Ini,px43);
                         ear = Ir;
                         eai = _mm512_add_pd(_mm512_mul_pd(Ii,_2pi),x16);
                         cexp_zmm8r8(ear,eai,&t3r,&t3i); 
                         cmul_zmm8r8(ix43r,ix43i,t3r,t3i,&cearr,&ceari);// ei2pi(x+1/6)
                         t0r = _mm512_fmadd_pd(pnx23,c1r,c0r);
                         t0i = _mm512_fmadd_pd(pnx23,c1i,c0i);
                         t1r = _mm512_sub_pd(t0r,_mm512_mul_pd(pnx43,c2r));
                         t1i = _mm512_sub_pd(t0i,_mm512_mul_pd(pnx43,c2i));
                         t2r = _mm512_mul_pd(ix43r,t1r); // first term re
                         t2i = _mm512_mul_pd(ix43i,t1i); // first term im
                         t0r = _mm512_fmadd_pd(xpn23,c3r,c0);
                         t0i = _mm512_fmadd_pd(xpn23,c3i,c0);
                         t4r = _mm512_fmadd_pd(npx13,c4r,
                                                    _mm512_mul_pd(pnx13,c5r));
                         t4i = _mm512_fmadd_pd(npx13,c4i,
                                                    _mm512_mul_pd(pnx13,c5i));
                         cexp_zmm8r8(t4r,t4i,&exp1r,&exp1i);
                         cmul_zmm8r8(t0r,t0i,exp1r,exp1i,&E1r,&E1i);//t3r,t3i hold second exp term
                         t1r = _mm512_fmadd_pd(pnx23,c6r,c1);
                         t1i = _mm512_fmadd_pd(pnx23,c6i,c1);
                         t5r = _mm512_fmsub_pd(npx13,c7r,
                                                     _mm512_mul_pd(pnx13,c8r));
                         t5i = _mm512_fmsub_pd(npx13,c7i,
                                                     _mm512_mul_pd(pnx13,c8i));
                         cexp_zmm8r8(t5r,t5i,&exp2r,&exp2i);
                         cmul_zmm8r8(t1r,t1i,exp2r,exp2i,&H1r,&H1i);
                         t3r = _mm512_add_pd(E1r,H1r);
                         t3i = _mm512_add_pd(E1i,H1i);
                         cmul_zmm8r8(cearr,ceari,t3r,t3i,&t1r,&t1i);
                         _mm512_storeu_pd(&FCr[0] ,_mm512_sub_pd(t2r,t1r));
                         _mm512_storeu_pd(&FCi[0] ,_mm512_sub_pd(t2i,t1i));
               } 


                 /*
                         Low frquency region (k0a < 0.4), Rayleigh approximation
                         for forward scattering function and cross-section.
                         Formulae: 3.2-31, 3.2-32
                         Forward scattering function.
                   */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d F_f3231_zmm8r8(const __m512d k0a) {

                          const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                          register __m512d k0a3,Fpi;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          Fpi  = _mm512_mul_pd(nhlf,k0a3);
                          return (Fpi); 
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void F_f3231_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                          double * __restrict __ATTR_ALIGN__(64) Fpi) {

                          const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                          register __m512d k0a3,Fpi;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          _mm512_store_pd(&Fpi[0] ,_mm512_mul_pd(nhlf,k0a3));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void F_f3231_zmm8r8_u(const double * __restrict  pk0a,
                                          double * __restrict  Fpi) {

                          const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const register __m512d nhlf = _mm512_set1_pd(-0.5f);
                          register __m512d k0a3,Fpi;
                          k0a3 = _mm512_mul_pd(k0a,_mm512_mul_pd(k0a,k0a));
                          _mm512_storeu_pd(&Fpi[0] ,_mm512_mul_pd(nhlf,k0a3));
                 }


                    /*
                         Low frquency region (k0a < 0.4), Rayleigh approximation
                         for forward scattering function and cross-section.
                         Formulae: 3.2-31, 3.2-32
                         RCS.
                   */


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3232_zmm8r8(const __m512d k0a,
                                            const __m512d a) {

                          register __m512d a2,ka04,rcs,t0;
                          a2   = _mm512_mul_pd(a,a);
                          t0   = _mm512_mul_pd(ka0,ka0);
                          ka04 = _mm512_mul_pd(t0,t0);
                          rcs  = _mm512_mul_pd(pi,_mm512_mul_pd(a2,ka04));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3232_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            double * __restrict __ATTR_ALIGN__(64) rcs) {

                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d a2,ka04,t0;
                          a2   = _mm512_mul_pd(a,a);
                          t0   = _mm512_mul_pd(ka0,ka0);
                          ka04 = _mm512_mul_pd(t0,t0);
                          _mm512_store_pd(&rcs[0] ,_mm512_mul_pd(
                                                     pi,_mm512_mul_pd(a2,ka04)));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3232_zmm8r8_u(  const double * __restrict  pk0a,
                                              const double * __restrict    pa,
                                              double * __restrict  rcs) {

                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d a2,ka04,t0;
                          a2   = _mm512_mul_pd(a,a);
                          t0   = _mm512_mul_pd(ka0,ka0);
                          ka04 = _mm512_mul_pd(t0,t0);
                          _mm512_storeu_pd(&rcs[0] ,_mm512_mul_pd(
                                                     pi,_mm512_mul_pd(a2,ka04)));
                 }


                   /*
                         High-frequency region -- forward scattering function and 
                         cross-section (k0a > 20)
                         Formula 3.2-33 (Forward scattering function).
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void F_f3233_zmm8r8(const __m512d k0a,
                                        __m512d * __restrict Fr,
                                        __m512d * __restrict Fi) {

                        const register __m512d hlf = _mm512_set1_pd(0.5f);
                        register __m512d k0a2;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        *Fr  = nIr;
                        *Fi  = _mm512_mul_pd(nIi,_mm512_mul_pd(hlf,k0a2));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void F_f3233_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                          double * __restrict __ATTR_ALIGN__(64) Fr,
                                          double * __restrict __ATTR_ALIGN__(64) Fi) {

                        register __m512d k0a       = _mm512_load_pd(&pk0a[0]);
                        const register __m512d hlf = _mm512_set1_pd(0.5f);
                        register __m512d k0a2;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        _mm512_store_pd(&Fr[0]  ,nIr);
                        _mm512_store_pd(&Fi[0]  ,_mm512_mul_pd(nIi,_mm512_mul_pd(hlf,k0a2)));
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void (const double * __restrict  pk0a,
                                          double * __restrict  Fr,
                                          double * __restrict  Fi) {

                        register __m512d k0a       = _mm512_loadu_pd(&pk0a[0]);
                        const register __m512d hlf = _mm512_set1_pd(0.5f);
                        register __m512d k0a2;
                        k0a2 = _mm512_mul_pd(k0a,k0a);
                        _mm512_storeu_pd(&Fr[0]  ,nIr);
                        _mm512_storeu_pd(&Fi[0]  ,_mm512_mul_pd(nIi,_mm512_mul_pd(hlf,k0a2)));
                 }

                   /*
                         High-frequency region -- forward scattering function and 
                         cross-section (k0a > 20)
                         Formula 3.2-34 (RCS).
                     */

 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3234_zmm8r8(const __m512d k0a,
                                            const __m512d a) {

                          register __m512d a2,ka02,rcs;
                          a2   = _mm512_mul_pd(a,a);
                          ka02   = _mm512_mul_pd(ka0,ka0);
                          rcs  = _mm512_mul_pd(pi,_mm512_mul_pd(a2,ka02));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3234_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            double * __restrict __ATTR_ALIGN__(64) rcs) {

                          register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                          register __m512d a   = _mm512_load_pd(&pa[0]);
                          register __m512d a2,ka02,rcs;
                          a2   = _mm512_mul_pd(a,a);
                          ka02   = _mm512_mul_pd(ka0,ka0);
                          _mm512_store_pd(&rcs[0], _mm512_mul_pd(pi,_mm512_mul_pd(a2,ka02)));
                          
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f3234_zmm8r8_u(const double * __restrict  pk0a,
                                            const double * __restrict  pa,
                                            double * __restrict  rcs) {

                          register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                          register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          register __m512d a2,ka02,rcs;
                          a2   = _mm512_mul_pd(a,a);
                          ka02   = _mm512_mul_pd(ka0,ka0);
                          _mm512_storeu_pd(&rcs[0], _mm512_mul_pd(pi,_mm512_mul_pd(a2,ka02)));
                          
                 }


                  /*
                          Low-frequency region (k1a < 0.8).
                          Expansion by two series terms i.e. A0,A1 and B0,B1.
                          Formula 3.3-5
                    */
                   
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A_coffs_f335_zmm8r8(const __m512d k0a5,
                                         const __m512d m1r,
                                         const __m512d m1i,
                                         __m512d * __restrict A1r,
                                         __m512d * __restrict A1i,
                                         __m512d * __restrict A2r,
                                         __m512d * __restrict A2i) {

                        const register __m512d c0 = _mm512_set1_pd(0.033333333333333333333333333333);
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d c0r,c0i;
                        c0r = _mm512_mul_pd(_mm512_sub_pd(m1r,_1),k0a5);
                        c0i = _mm512_mul_pd(_mm512_sub_pd(m1i,_1),k0a5);
                        *A1r = _mm512_mul_pd(c0r,c0);
                        *A1i = _mm512_mul_pd(c0i,c0);
                        *A2r = Ir;
                        *A2i = Ir
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A_coffs_f335_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0a5,
                                               const double * __restrict __ATTR_ALIGN__(64) pm1r,
                                               const double * __restrict __ATTR_ALIGN__(64) pm1i,
                                               double * __restrict __ATTR_ALIGN__(64) A1r,
                                               double * __restrict __ATTR_ALIGN__(64) A1i,
                                               double * __restrict __ATTR_ALIGN__(64) A2r,
                                               double * __restrict __ATTR_ALIGN__(64) A2i) {

                        register __m512d k0a5     = _mm512_load_pd(&pk0a5[0]);
                        register __m512d m1r      = _mm512_load_pd(&pm1r[0]); 
                        register __m512d m1i      = _mm512_load_pd(&pm1i[0]);
                        const register __m512d c0 = _mm512_set1_pd(0.033333333333333333333333333333);
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d c0r,c0i;
                        c0r = _mm512_mul_pd(_mm512_sub_pd(m1r,_1),k0a5);
                        c0i = _mm512_mul_pd(_mm512_sub_pd(m1i,_1),k0a5);
                        _mm512_store_pd(&A1r[0], _mm512_mul_pd(c0r,c0));
                        _mm512_store_pd(&A1i[0], _mm512_mul_pd(c0i,c0));
                        _mm512_store_pd(&A2r[0], Ir);
                        _mm512_store_pd(&A2i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void A_coffs_f335_zmm8r8_u(const double * __restrict  pk0a5,
                                               const double * __restrict  pm1r,
                                               const double * __restrict  pm1i,
                                               double * __restrict  A1r,
                                               double * __restrict  A1i,
                                               double * __restrict  A2r,
                                               double * __restrict  A2i) {

                        register __m512d k0a5     = _mm512_loadu_pd(&pk0a5[0]);
                        register __m512d m1r      = _mm512_loadu_pd(&pm1r[0]); 
                        register __m512d m1i      = _mm512_loadu_pd(&pm1i[0]);
                        const register __m512d c0 = _mm512_set1_pd(0.033333333333333333333333333333);
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        register __m512d c0r,c0i;
                        c0r = _mm512_mul_pd(_mm512_sub_pd(m1r,_1),k0a5);
                        c0i = _mm512_mul_pd(_mm512_sub_pd(m1i,_1),k0a5);
                        _mm512_storeu_pd(&A1r[0], _mm512_mul_pd(c0r,c0));
                        _mm512_storeu_pd(&A1i[0], _mm512_mul_pd(c0i,c0));
                        _mm512_storeu_pd(&A2r[0], Ir);
                        _mm512_storeu_pd(&A2i[0], Ir);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B_coffs_f335_zmm8r8(const __m512d k0a3,
                                             const __m512d k0a5,
                                             const __m512d m1r,
                                             const __m512d m1i,
                                             __m512d * __restrict B1r,
                                             __m512d * __restrict B1i,
                                             __m512d * __restrict B2r,
                                             __m512d * __restrict B2i) {
                 
                        register __m512d mm1r,mm1i,m1s1r,m1s1i,m1s2r,m1s2i;
                        register __m512d m1a2r,m1a2i,_2m1r,_2m1i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mm1r                     = _mm512_mul_pd(m1r,m1r);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mm1i                     = _mm512_mul_pd(m1i,m1i);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        m1s1r                    = _mm512_sub_pd(mm1r,_1);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        m1s1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d c1 = _mm512_set1_pd(0.055555555555555555555555555556);
                        m1s2r                    = _mm512_sub_pd(mm1r,_2);
                        register __m512d t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
                        m1s2i                    = _mm512_sub_pd(mm1i,_2);
                        m1a2r                    = _mm512_add_pd(mm1r,_2);
                        m1a2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(m1a2r,m1a2i,m1a2r,m1a2i,&t3r,&t3i);
                        _2m1r                    = _mm512_fmadd_pd(_2,mm1r,_3);
                        _2m1i                    = _mm512_fmadd_pd(_2,mm1i,_3);
                        cdiv_zmm8r8(m1s1r,m1s1i,m1a2r,m1a2i,&t0r,&t0i);
                        t0r = _mm512_mul_pd(t0r,k0a3);
                        t0i = _mm512_mul_pd(t0i,k0a3);
                        cmul_zmm8r8(m1s1r,m1s1i,m1s2r,m1s2i,&t1r,&t1i);
                        cdiv_zmm8r8(t1r,t1i,t3r,t3i,&t2r,&t2i);
                        t2r = _mm512_mul_pd(c0,_mm512_mul_pd(t2r,k0a5));
                        t2i = _mm512_mul_pd(c0,_mm512_mul_pd(t2i,k0a5));
                        *B1r = _mm512_add_pd(t0r,t2r);
                        *B1i = _mm512_add_pd(t0i,t2i);
                        cdiv_zmm8r8(m1s1r,m1s1i,_2m1r,_2m1i,&t0r,&t0i);
                        *B2r = _mm512_mul_pd(c1,_mm512_mul_pd(t0r,k0a5));
                        *B2i = _mm512_mul_pd(c1,_mm512_mul_pd(t0i,k0a5));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B_coffs_f335_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64)  pk0a3,
                                               const double * __restrict __ATTR_ALIGN__(64)  pk0a5,
                                               const double * __restrict __ATTR_ALIGN__(64)  pm1r,
                                               const double * __restrict __ATTR_ALIGN__(64)  pm1i,
                                               double * __restrict __ATTR_ALIGN__(64) B1r,
                                               double * __restrict __ATTR_ALIGN__(64) B1i,
                                               double * __restrict __ATTR_ALIGN__(64) B2r,
                                               double * __restrict __ATTR_ALIGN__(64) B2i) {
                 
                        register __m512d k0a3 = _mm512_load_pd(&pk0a3[0]);
                        register __m512d k0a5 = _mm512_load_pd(&pk0a5[0]);
                        register __m512d m1r  = _mm512_load_pd(&pm1r[0]);
                        register __m512d m1i  = _mm512_load_pd(&pm1i[0]);
                        register __m512d mm1r,mm1i,m1s1r,m1s1i,m1s2r,m1s2i;
                        register __m512d m1a2r,m1a2i,_2m1r,_2m1i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mm1r                     = _mm512_mul_pd(m1r,m1r);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mm1i                     = _mm512_mul_pd(m1i,m1i);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        m1s1r                    = _mm512_sub_pd(mm1r,_1);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        m1s1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d c1 = _mm512_set1_pd(0.055555555555555555555555555556);
                        m1s2r                    = _mm512_sub_pd(mm1r,_2);
                        register __m512d t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
                        m1s2i                    = _mm512_sub_pd(mm1i,_2);
                        m1a2r                    = _mm512_add_pd(mm1r,_2);
                        m1a2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(m1a2r,m1a2i,m1a2r,m1a2i,&t3r,&t3i);
                        _2m1r                    = _mm512_fmadd_pd(_2,mm1r,_3);
                        _2m1i                    = _mm512_fmadd_pd(_2,mm1i,_3);
                        cdiv_zmm8r8(m1s1r,m1s1i,m1a2r,m1a2i,&t0r,&t0i);
                        t0r = _mm512_mul_pd(t0r,k0a3);
                        t0i = _mm512_mul_pd(t0i,k0a3);
                        cmul_zmm8r8(m1s1r,m1s1i,m1s2r,m1s2i,&t1r,&t1i);
                        cdiv_zmm8r8(t1r,t1i,t3r,t3i,&t2r,&t2i);
                        t2r = _mm512_mul_pd(c0,_mm512_mul_pd(t2r,k0a5));
                        t2i = _mm512_mul_pd(c0,_mm512_mul_pd(t2i,k0a5));
                        _mm512_store_pd(&B1r[0], _mm512_add_pd(t0r,t2r));
                        _mm512_store_pd(&B1i[0], _mm512_add_pd(t0i,t2i));
                        cdiv_zmm8r8(m1s1r,m1s1i,_2m1r,_2m1i,&t0r,&t0i);
                        _mm512_store_pd(&B2r[0], _mm512_mul_pd(c1,_mm512_mul_pd(t0r,k0a5)));
                        _mm512_store_pd(&B2i[0], _mm512_mul_pd(c1,_mm512_mul_pd(t0i,k0a5)));
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void B_coffs_f335_zmm8r8_u(const double * __restrict   pk0a3,
                                               const double * __restrict   pk0a5,
                                               const double * __restrict   pm1r,
                                               const double * __restrict   pm1i,
                                               double * __restrict  B1r,
                                               double * __restrict  B1i,
                                               double * __restrict  B2r,
                                               double * __restrict  B2i) {
                 
                        register __m512d k0a3 = _mm512_loadu_pd(&pk0a3[0]);
                        register __m512d k0a5 = _mm512_loadu_pd(&pk0a5[0]);
                        register __m512d m1r  = _mm512_loadu_pd(&pm1r[0]);
                        register __m512d m1i  = _mm512_loadu_pd(&pm1i[0]);
                        register __m512d mm1r,mm1i,m1s1r,m1s1i,m1s2r,m1s2i;
                        register __m512d m1a2r,m1a2i,_2m1r,_2m1i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mm1r                     = _mm512_mul_pd(m1r,m1r);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mm1i                     = _mm512_mul_pd(m1i,m1i);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        m1s1r                    = _mm512_sub_pd(mm1r,_1);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        m1s1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d c1 = _mm512_set1_pd(0.055555555555555555555555555556);
                        m1s2r                    = _mm512_sub_pd(mm1r,_2);
                        register __m512d t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
                        m1s2i                    = _mm512_sub_pd(mm1i,_2);
                        m1a2r                    = _mm512_add_pd(mm1r,_2);
                        m1a2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(m1a2r,m1a2i,m1a2r,m1a2i,&t3r,&t3i);
                        _2m1r                    = _mm512_fmadd_pd(_2,mm1r,_3);
                        _2m1i                    = _mm512_fmadd_pd(_2,mm1i,_3);
                        cdiv_zmm8r8(m1s1r,m1s1i,m1a2r,m1a2i,&t0r,&t0i);
                        t0r = _mm512_mul_pd(t0r,k0a3);
                        t0i = _mm512_mul_pd(t0i,k0a3);
                        cmul_zmm8r8(m1s1r,m1s1i,m1s2r,m1s2i,&t1r,&t1i);
                        cdiv_zmm8r8(t1r,t1i,t3r,t3i,&t2r,&t2i);
                        t2r = _mm512_mul_pd(c0,_mm512_mul_pd(t2r,k0a5));
                        t2i = _mm512_mul_pd(c0,_mm512_mul_pd(t2i,k0a5));
                        _mm512_storeu_pd(&B1r[0], _mm512_add_pd(t0r,t2r));
                        _mm512_storeu_pd(&B1i[0], _mm512_add_pd(t0i,t2i));
                        cdiv_zmm8r8(m1s1r,m1s1i,_2m1r,_2m1i,&t0r,&t0i);
                        _mm512_storeu_pd(&B2r[0], _mm512_mul_pd(c1,_mm512_mul_pd(t0r,k0a5)));
                        _mm512_storeu_pd(&B2i[0], _mm512_mul_pd(c1,_mm512_mul_pd(t0i,k0a5)));
                }


                   /*
                         Rayleigh backscattering RCS for dielectric spheres at angle 0.
                         Formula 3.3-7
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f337_zmm8r8(const __m512d a,
                                           const __m512d k0a4,
                                           const __m512d m1r,
                                           const __m512d m1i) {

                          register __m512d aa       = _mm512_mul_pd(a,a);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d _2 = _mm512_set1_pd(2.0f);
                          const register __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d fac      = _mm512_mul_pd(_4,
                                                               _mm512_mul_pd(PI,aa));
                          register __m512d mm1r,mm1i,mma2r,mma2i,mms1r,mms1i,t0r,t0i;
                          register __m512d cabs,rcs;
                          cmul_zmm8r8(m1r,m1i,m1r,m1i,&mm1r,&mm1i);
                          mms1r = _mm512_sub_pd(mm1r,_1);
                          mms1i = _mm512_sub_pd(mm1i,_1);
                          mma2r = _mm512_add_pd(mm1r,_2);
                          mma2i = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mm1sr,mm1si,mma2r,mma2i,&t0r,&t0i);
                          cabs = cabs_zmm8r8(t0r,t0i);
                          rcs  = _mm512_mul_pd(fac,_mm512_mul_pd(cabs,k0a4));
                          return (rcs);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f337_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pa,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a4,
                                           const double * __restrict __ATTR_ALIGN__(64) pm1r,
                                           const double * __restrict __ATTR_ALIGN__(64) pm1i,
                                           double * __restrict __ATTR_ALIGN__(64) rcs ) {

                          register __m512d a        = _mm512_load_pd(&pa[0]);
                          register __m512d aa       = _mm512_mul_pd(a,a);
                          register __m512d k0a4     = _mm512_load_pd(&pk0a4[0]);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d m1r      = _mm512_load_pd(&pm1r[0]);
                          const register __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d m1i      = _mm512_load_pd(&pm1i[0]);
                          const register __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d fac      = _mm512_mul_pd(_4,
                                                               _mm512_mul_pd(PI,aa));
                          register __m512d mm1r,mm1i,mma2r,mma2i,mms1r,mms1i,t0r,t0i;
                          register __m512d cabs,rcs;
                          cmul_zmm8r8(m1r,m1i,m1r,m1i,&mm1r,&mm1i);
                          mms1r = _mm512_sub_pd(mm1r,_1);
                          mms1i = _mm512_sub_pd(mm1i,_1);
                          mma2r = _mm512_add_pd(mm1r,_2);
                          mma2i = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mm1sr,mm1si,mma2r,mma2i,&t0r,&t0i);
                          cabs = cabs_zmm8r8(t0r,t0i);
                          _mm512_store_pd(&rcs[0], _mm512_mul_pd(fac,_mm512_mul_pd(cabs,k0a4)));
                    
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void rcs_f337_zmm8r8_u(const double * __restrict  pa,
                                           const double * __restrict  pk0a4,
                                           const double * __restrict  pm1r,
                                           const double * __restrict  pm1i,
                                           double * __restrict  rcs ) {

                          register __m512d a        = _mm512_loadu_pd(&pa[0]);
                          register __m512d aa       = _mm512_mul_pd(a,a);
                          register __m512d k0a4     = _mm512_loadu_pd(&pk0a4[0]);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          register __m512d m1r      = _mm512_loadu_pd(&pm1r[0]);
                          const register __m512d _2 = _mm512_set1_pd(2.0f);
                          register __m512d m1i      = _mm512_loadu_pd(&pm1i[0]);
                          const register __m512d _4 = _mm512_set1_pd(4.0f);
                          register __m512d fac      = _mm512_mul_pd(_4,
                                                               _mm512_mul_pd(PI,aa));
                          register __m512d mm1r,mm1i,mma2r,mma2i,mms1r,mms1i,t0r,t0i;
                          register __m512d cabs,rcs;
                          cmul_zmm8r8(m1r,m1i,m1r,m1i,&mm1r,&mm1i);
                          mms1r = _mm512_sub_pd(mm1r,_1);
                          mms1i = _mm512_sub_pd(mm1i,_1);
                          mma2r = _mm512_add_pd(mm1r,_2);
                          mma2i = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mm1sr,mm1si,mma2r,mma2i,&t0r,&t0i);
                          cabs = cabs_zmm8r8(t0r,t0i);
                          _mm512_storeu_pd(&rcs[0], _mm512_mul_pd(fac,_mm512_mul_pd(cabs,k0a4)));
                    
               }


                 /*
                        Low-frequency bi-static scattering
                  */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f338_zmm8r8(const __m512d ka03,
                                        const __m512d ka05,
                                        const __m512d tht,
                                        const __m512d mm1r,
                                        const __m512d mm1i,
                                        __m512d * __restrict S1r,
                                        __m512d * __restrict S1i) {

                        register __m512d mms1r,mms1i,mms2r,mms2i;
                        register __m512d mma2r,mma2i,mma3r,mma3i;
                        register __m512d div1r,div1i,div2r,div2i;
                        register __m512d div3r,div3i,mulr,muli,t1r,t1i;
                        register __m512d cost,cos2t,msqr,msqi,t0r,t0i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mms1r                    = _mm512_sub_pd(mm1r,_1);
                        mms1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mms2r                    = _mm512_sub_pd(mm1r,_2);
                        mms2i                    = _mm512_sub_pd(mm1i,_2);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        mma2r                    = _mm512_add_pd(mm1r,_2);
                        mma2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        mma3r                    = _mm512_add_pd(mm1r,_3);
                        mma3i                    = _mm512_add_pd(mm1i,_3);
                        const register __m512d c1 = _mm512_set1_pd(0.166666666666666666666666666667);
                        cost                     = xcos(tht);
                        const register __m512d c2 = _mm512_set1_pd(0.033333333333333333333333333333);
                        cos2t                    = xcos(_mm512_add_pd(tht,tht));
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm512_mul_pd(div1r,_mm512_mul_pd(cost,ka03));
                        cmul_zmm8r8(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm512_mul_pd(div1i,_mm512_mul_pd(cost,ka03));
                        cdiv_zmm8r8(mulr,muli,msqr,msqi,&div2r,&div2i);
                        div2r = _mm512_mul_pd(div2r,cost);
                        div2i = _mm512_mul_pd(div2i,cost);
                        t0r   = _mm512_mul_pd(c2,mms1r);
                        t0i   = _mm512_mul_pd(c2,mms1i);
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm512_mul_pd(c1,_mm512_mul_pd(div3r,cos2t));
                        div3i = _mm512_mul_pd(c1,_mm512_mul_pd(div3i,cos2t));
                        t1r   = _mm512_sub_pd(div2r,_mm512_sub_pd(div3r,t0r));
                        t1i   = _mm512_sub_pd(div2i,_mm512_sub_pd(div3i,t0i));
                        *S1r  = _mm512_fmadd_pd(t1r,k0a5,div1r);
                        *S1i  = _mm512_fmadd_pd(t1i,k0a5,div1i);
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f338_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pka03,
                                        const double * __restrict __ATTR_ALIGN__(64) pka05,
                                        const double * __restrict __ATTR_ALIGN__(64) ptht,
                                        const double * __restrict __ATTR_ALIGN__(64) pmm1r,
                                        const double * __restrict __ATTR_ALIGN__(64) pmm1i,
                                        double * __restrict __ATTR_ALIGN__(64) S1r,
                                        double * __restrict __ATTR_ALIGN__(64) S1i) {

                        register __m512d ka03     = _mm512_load_pd(&pka03[0]);
                        register __m512d ka05     = _mm512_load_pd(&pka05[0]);
                        register __m512d tht      = _mm512_load_pd(&ptht[0]);
                        register __m512d mm1r     = _mm512_load_pd(&pmm1r[0]);
                        register __m512d mm1i     = _mm512_load_pd(&pmm1i[0]);
                        register __m512d mms1r,mms1i,mms2r,mms2i;
                        register __m512d mma2r,mma2i,mma3r,mma3i;
                        register __m512d div1r,div1i,div2r,div2i;
                        register __m512d div3r,div3i,mulr,muli,t1r,t1i;
                        register __m512d cost,cos2t,msqr,msqi,t0r,t0i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mms1r                    = _mm512_sub_pd(mm1r,_1);
                        mms1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mms2r                    = _mm512_sub_pd(mm1r,_2);
                        mms2i                    = _mm512_sub_pd(mm1i,_2);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        mma2r                    = _mm512_add_pd(mm1r,_2);
                        mma2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        mma3r                    = _mm512_add_pd(mm1r,_3);
                        mma3i                    = _mm512_add_pd(mm1i,_3);
                        const register __m512d c1 = _mm512_set1_pd(0.166666666666666666666666666667);
                        cost                     = xcos(tht);
                        const register __m512d c2 = _mm512_set1_pd(0.033333333333333333333333333333);
                        cos2t                    = xcos(_mm512_add_pd(tht,tht));
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm512_mul_pd(div1r,_mm512_mul_pd(cost,ka03));
                        cmul_zmm8r8(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm512_mul_pd(div1i,_mm512_mul_pd(cost,ka03));
                        cdiv_zmm8r8(mulr,muli,msqr,msqi,&div2r,&div2i);
                        div2r = _mm512_mul_pd(div2r,cost);
                        div2i = _mm512_mul_pd(div2i,cost);
                        t0r   = _mm512_mul_pd(c2,mms1r);
                        t0i   = _mm512_mul_pd(c2,mms1i);
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm512_mul_pd(c1,_mm512_mul_pd(div3r,cos2t));
                        div3i = _mm512_mul_pd(c1,_mm512_mul_pd(div3i,cos2t));
                        t1r   = _mm512_sub_pd(div2r,_mm512_sub_pd(div3r,t0r));
                        t1i   = _mm512_sub_pd(div2i,_mm512_sub_pd(div3i,t0i));
                        _mm512_store_pd(&S1r[0] ,_mm512_fmadd_pd(t1r,k0a5,div1r));
                        _mm512_store_pd(&S1i[0] ,_mm512_fmadd_pd(t1i,k0a5,div1i));
               }


                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f338_zmm8r8_u(const double * __restrict  pka03,
                                        const double * __restrict  pka05,
                                        const double * __restrict  ptht,
                                        const double * __restrict  pmm1r,
                                        const double * __restrict  pmm1i,
                                        double * __restrict  S1r,
                                        double * __restrict  S1i) {

                        register __m512d ka03     = _mm512_loadu_pd(&pka03[0]);
                        register __m512d ka05     = _mm512_loadu_pd(&pka05[0]);
                        register __m512d tht      = _mm512_loadu_pd(&ptht[0]);
                        register __m512d mm1r     = _mm512_loadu_pd(&pmm1r[0]);
                        register __m512d mm1i     = _mm512_loadu_pd(&pmm1i[0]);
                        register __m512d mms1r,mms1i,mms2r,mms2i;
                        register __m512d mma2r,mma2i,mma3r,mma3i;
                        register __m512d div1r,div1i,div2r,div2i;
                        register __m512d div3r,div3i,mulr,muli,t1r,t1i;
                        register __m512d cost,cos2t,msqr,msqi,t0r,t0i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mms1r                    = _mm512_sub_pd(mm1r,_1);
                        mms1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mms2r                    = _mm512_sub_pd(mm1r,_2);
                        mms2i                    = _mm512_sub_pd(mm1i,_2);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        mma2r                    = _mm512_add_pd(mm1r,_2);
                        mma2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        mma3r                    = _mm512_add_pd(mm1r,_3);
                        mma3i                    = _mm512_add_pd(mm1i,_3);
                        const register __m512d c1 = _mm512_set1_pd(0.166666666666666666666666666667);
                        cost                     = xcos(tht);
                        const register __m512d c2 = _mm512_set1_pd(0.033333333333333333333333333333);
                        cos2t                    = xcos(_mm512_add_pd(tht,tht));
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm512_mul_pd(div1r,_mm512_mul_pd(cost,ka03));
                        cmul_zmm8r8(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm512_mul_pd(div1i,_mm512_mul_pd(cost,ka03));
                        cdiv_zmm8r8(mulr,muli,msqr,msqi,&div2r,&div2i);
                        div2r = _mm512_mul_pd(div2r,cost);
                        div2i = _mm512_mul_pd(div2i,cost);
                        t0r   = _mm512_mul_pd(c2,mms1r);
                        t0i   = _mm512_mul_pd(c2,mms1i);
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm512_mul_pd(c1,_mm512_mul_pd(div3r,cos2t));
                        div3i = _mm512_mul_pd(c1,_mm512_mul_pd(div3i,cos2t));
                        t1r   = _mm512_sub_pd(div2r,_mm512_sub_pd(div3r,t0r));
                        t1i   = _mm512_sub_pd(div2i,_mm512_sub_pd(div3i,t0i));
                        _mm512_storeu_pd(&S1r[0] ,_mm512_fmadd_pd(t1r,k0a5,div1r));
                        _mm512_storeu_pd(&S1i[0] ,_mm512_fmadd_pd(t1i,k0a5,div1i));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f338_zmm8r8(const __m512d ka03,
                                        const __m512d ka05,
                                        const __m512d tht,
                                        const __m512d mm1r,
                                        const __m512d mm1i,
                                        __m512d * __restrict S2r,
                                        __m512d * __restrict S2i) {

                        register __m512d mms1r,mms1i,mms2r,mms2i;
                        register __m512d mma2r,mma2i,mma3r,mma3i;
                        register __m512d div1r,div1i,div2r,div2i;
                        register __m512d div3r,div3i,mulr,muli,t1r,t1i;
                        register __m512d cost,cos2t,msqr,msqi,t0r,t0i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mms1r                    = _mm512_sub_pd(mm1r,_1);
                        mms1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mms2r                    = _mm512_sub_pd(mm1r,_2);
                        mms2i                    = _mm512_sub_pd(mm1i,_2);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        mma2r                    = _mm512_add_pd(mm1r,_2);
                        mma2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        mma3r                    = _mm512_add_pd(mm1r,_3);
                        mma3i                    = _mm512_add_pd(mm1i,_3);
                        const register __m512d c1 = _mm512_set1_pd(0.166666666666666666666666666667);
                        cost                     = xcos(tht);
                        const register __m512d c2 = _mm512_set1_pd(0.033333333333333333333333333333);
                        cos2t                    = xcos(_mm512_add_pd(tht,tht));
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm512_mul_pd(div1r,ka03);
                        cmul_zmm8r8(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm512_mul_pd(div1i,ka03);
                        cdiv_zmm8r8(mulr,muli,msqr,msqi,&div2r,&div2i);
                        t0r   = _mm512_mul_pd(c2,mms1r);
                        t0i   = _mm512_mul_pd(c2,mms1i);
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm512_mul_pd(c1,_mm512_mul_pd(div3r,cost));
                        div3i = _mm512_mul_pd(c1,_mm512_mul_pd(div2i,cost));
                        t1r   = _mm512_sub_pd(div2r,_mm512_sub_pd(div3r,t0r));
                        t1i   = _mm512_sub_pd(div2i,_mm512_sub_pd(div3i,t0i));
                        *S2r  = _mm512_fmadd_pd(t1r,cost,div1r);
                        *S2i  = _mm512_fmadd_pd(t1i,cost,div1i);
               }


                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f338_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pka03,
                                          const double * __restrict __ATTR_ALIGN__(64) pka05,
                                          const double * __restrict __ATTR_ALIGN__(64) ptht,
                                          const double * __restrict __ATTR_ALIGN__(64) pmm1r,
                                          const double * __restrict __ATTR_ALIGN__(64) pmm1i,
                                          double * __restrict __ATTR_ALIGN__(64) S2r,
                                          double * __restrict __ATTR_ALIGN__(64) S2i) {

                        register __m512d ka03     = _mm512_load_pd(&pka03[0]);
                        register __m512d ka05     = _mm512_load_pd(&pka05[0]);
                        register __m512d tht      = _mm512_load_pd(&ptht[0]);
                        register __m512d mm1r     = _mm512_load_pd(&pmm1r[0]);
                        register __m512d mm1i     = _mm512_load_pd(&pmm1i[0]);
                        register __m512d mms1r,mms1i,mms2r,mms2i;
                        register __m512d mma2r,mma2i,mma3r,mma3i;
                        register __m512d div1r,div1i,div2r,div2i;
                        register __m512d div3r,div3i,mulr,muli,t1r,t1i;
                        register __m512d cost,cos2t,msqr,msqi,t0r,t0i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mms1r                    = _mm512_sub_pd(mm1r,_1);
                        mms1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mms2r                    = _mm512_sub_pd(mm1r,_2);
                        mms2i                    = _mm512_sub_pd(mm1i,_2);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        mma2r                    = _mm512_add_pd(mm1r,_2);
                        mma2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        mma3r                    = _mm512_add_pd(mm1r,_3);
                        mma3i                    = _mm512_add_pd(mm1i,_3);
                        const register __m512d c1 = _mm512_set1_pd(0.166666666666666666666666666667);
                        cost                     = xcos(tht);
                        const register __m512d c2 = _mm512_set1_pd(0.033333333333333333333333333333);
                        cos2t                    = xcos(_mm512_add_pd(tht,tht));
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm512_mul_pd(div1r,ka03);
                        cmul_zmm8r8(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm512_mul_pd(div1i,ka03);
                        cdiv_zmm8r8(mulr,muli,msqr,msqi,&div2r,&div2i);
                        t0r   = _mm512_mul_pd(c2,mms1r);
                        t0i   = _mm512_mul_pd(c2,mms1i);
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm512_mul_pd(c1,_mm512_mul_pd(div3r,cost));
                        div3i = _mm512_mul_pd(c1,_mm512_mul_pd(div2i,cost));
                        t1r   = _mm512_sub_pd(div2r,_mm512_sub_pd(div3r,t0r));
                        t1i   = _mm512_sub_pd(div2i,_mm512_sub_pd(div3i,t0i));
                        _mm512_store_pd(&S2r[0] ,_mm512_fmadd_pd(t1r,cost,div1r));
                        _mm512_store_pd(&S2i[0] ,_mm512_fmadd_pd(t1i,cost,div1i));
               }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f338_zmm8r8_u(const double * __restrict  pka03,
                                          const double * __restrict  pka05,
                                          const double * __restrict  ptht,
                                          const double * __restrict  pmm1r,
                                          const double * __restrict  pmm1i,
                                          double * __restrict  S2r,
                                          double * __restrict  S2i) {

                        register __m512d ka03     = _mm512_loadu_pd(&pka03[0]);
                        register __m512d ka05     = _mm512_loadu_pd(&pka05[0]);
                        register __m512d tht      = _mm512_loadu_pd(&ptht[0]);
                        register __m512d mm1r     = _mm512_loadu_pd(&pmm1r[0]);
                        register __m512d mm1i     = _mm512_loadu_pd(&pmm1i[0]);
                        register __m512d mms1r,mms1i,mms2r,mms2i;
                        register __m512d mma2r,mma2i,mma3r,mma3i;
                        register __m512d div1r,div1i,div2r,div2i;
                        register __m512d div3r,div3i,mulr,muli,t1r,t1i;
                        register __m512d cost,cos2t,msqr,msqi,t0r,t0i;
                        const register __m512d _1 = _mm512_set1_pd(1.0f);
                        mms1r                    = _mm512_sub_pd(mm1r,_1);
                        mms1i                    = _mm512_sub_pd(mm1i,_1);
                        const register __m512d _2 = _mm512_set1_pd(2.0f);
                        mms2r                    = _mm512_sub_pd(mm1r,_2);
                        mms2i                    = _mm512_sub_pd(mm1i,_2);
                        const register __m512d _3 = _mm512_set1_pd(3.0f);
                        mma2r                    = _mm512_add_pd(mm1r,_2);
                        mma2i                    = _mm512_add_pd(mm1i,_2);
                        cmul_zmm8r8(mma2r,mma2i,mma2r,mma2i,&msqr,&msqi);
                        const register __m512d c0 = _mm512_set1_pd(0.6f);
                        mma3r                    = _mm512_add_pd(mm1r,_3);
                        mma3i                    = _mm512_add_pd(mm1i,_3);
                        const register __m512d c1 = _mm512_set1_pd(0.166666666666666666666666666667);
                        cost                     = xcos(tht);
                        const register __m512d c2 = _mm512_set1_pd(0.033333333333333333333333333333);
                        cos2t                    = xcos(_mm512_add_pd(tht,tht));
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div1r,&div1i);
                        div1r = _mm512_mul_pd(div1r,ka03);
                        cmul_zmm8r8(mms1r,mms1i,mms2r,mms2i,&mulr,&muli);
                        div1i = _mm512_mul_pd(div1i,ka03);
                        cdiv_zmm8r8(mulr,muli,msqr,msqi,&div2r,&div2i);
                        t0r   = _mm512_mul_pd(c2,mms1r);
                        t0i   = _mm512_mul_pd(c2,mms1i);
                        cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&div3r,&div3i);
                        div3r = _mm512_mul_pd(c1,_mm512_mul_pd(div3r,cost));
                        div3i = _mm512_mul_pd(c1,_mm512_mul_pd(div2i,cost));
                        t1r   = _mm512_sub_pd(div2r,_mm512_sub_pd(div3r,t0r));
                        t1i   = _mm512_sub_pd(div2i,_mm512_sub_pd(div3i,t0i));
                        _mm512_storeu_pd(&S2r[0] ,_mm512_fmadd_pd(t1r,cost,div1r));
                        _mm512_storeu_pd(&S2i[0] ,_mm512_fmadd_pd(t1i,cost,div1i));
               }

                   /*
                         E-plane and H-plane RCS.
                         Formulae: 3.3-10,3.3-11
                     */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3310_zmm8r8(const __m512d tht,
                                            const __m512d a,
                                            const __m512d ka04,
                                            const __m512d mm1r,
                                            const __m512d mm1i) {

                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d aa = _mm512_mul_pd(a,a);
                          const register __m512d _2 = _mm512_set1_pd(2.0f); 
                          const register __m512d cost = xcos(tht);
                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d cos2t,mms1r,mms1i,mma2r,mma2i;
                          register __m512d divr,divi,rcs,frac,cabs;
                          cos2t   = _mm512_mul_pd(cost,cost);
                          mms1r   = _mm512_sub_pd(mm1r,_1);
                          mma2r   = _mm512_add_pd(mm1r,_2);
                          mms1i   = _mm512_sub_pd(mm1i,_1);
                          mma2i   = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm512_mul_pd(_4pi,aa);
                          cabs    = cabs_zmm8r8(divr,divi);
                          rcs     = _mm512_mul_pd(_mm512_mul_pd(frac,cabs),
                                                  _mm512_mul_pd(cos2t,ka04));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3310_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                            const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pka04,
                                            const double * __restrict __ATTR_ALIGN__(64) pmm1r,
                                            const double * __restrict __ATTR_ALIGN__(64) pmm1i) {

                          const register __m512d tht = _mm512_load_pd(&ptht[0]);
                          const register __m512d a   = _mm512_load_pd(&pa[0]);
                          const register __m512d ka04= _mm512_load_pd(&pka04[0]);
                          const register __m512d mm1r= _mm512_load_pd(&pmm1r[0]);
                          const register __m512d mm1i= _mm512_load_pd(&pmm1i[0]);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d aa = _mm512_mul_pd(a,a);
                          const register __m512d _2 = _mm512_set1_pd(2.0f); 
                          const register __m512d cost = xcos(tht);
                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d cos2t,mms1r,mms1i,mma2r,mma2i;
                          register __m512d divr,divi,rcs,frac,cabs;
                          cos2t   = _mm512_mul_pd(cost,cost);
                          mms1r   = _mm512_sub_pd(mm1r,_1);
                          mma2r   = _mm512_add_pd(mm1r,_2);
                          mms1i   = _mm512_sub_pd(mm1i,_1);
                          mma2i   = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm512_mul_pd(_4pi,aa);
                          cabs    = cabs_zmm8r8(divr,divi);
                          rcs     = _mm512_mul_pd(_mm512_mul_pd(frac,cabs),
                                                  _mm512_mul_pd(cos2t,ka04));
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3310_zmm8r8_u(const double * __restrict  ptht,
                                              const double * __restrict  pa,
                                              const double * __restrict  pka04,
                                              const double * __restrict  pmm1r,
                                              const double * __restrict  pmm1i) {

                          const register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          const register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          const register __m512d ka04= _mm512_loadu_pd(&pka04[0]);
                          const register __m512d mm1r= _mm512_loadu_pd(&pmm1r[0]);
                          const register __m512d mm1i= _mm512_loadu_pd(&pmm1i[0]);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d aa = _mm512_mul_pd(a,a);
                          const register __m512d _2 = _mm512_set1_pd(2.0f); 
                          const register __m512d cost = xcos(tht);
                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d cos2t,mms1r,mms1i,mma2r,mma2i;
                          register __m512d divr,divi,rcs,frac,cabs;
                          cos2t   = _mm512_mul_pd(cost,cost);
                          mms1r   = _mm512_sub_pd(mm1r,_1);
                          mma2r   = _mm512_add_pd(mm1r,_2);
                          mms1i   = _mm512_sub_pd(mm1i,_1);
                          mma2i   = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm512_mul_pd(_4pi,aa);
                          cabs    = cabs_zmm8r8(divr,divi);
                          rcs     = _mm512_mul_pd(_mm512_mul_pd(frac,cabs),
                                                  _mm512_mul_pd(cos2t,ka04));
                          return (rcs); 
                 }


                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3311_zmm8r8(const __m512d tht,
                                            const __m512d a,
                                            const __m512d ka04,
                                            const __m512d mm1r,
                                            const __m512d mm1i) {

                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d aa = _mm512_mul_pd(a,a);
                          const register __m512d _2 = _mm512_set1_pd(2.0f); 
                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d mms1r,mms1i,mma2r,mma2i;
                          register __m512d divr,divi,rcs,frac,cabs;
                          mms1r   = _mm512_sub_pd(mm1r,_1);
                          mma2r   = _mm512_add_pd(mm1r,_2);
                          mms1i   = _mm512_sub_pd(mm1i,_1);
                          mma2i   = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm512_mul_pd(_4pi,aa);
                          cabs    = cabs_zmm8r8(divr,divi);
                          rcs     = _mm512_mul_pd(ka04,_mm512_mul_pd(frac,cabs));
                                                  
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3311_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ptht,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pka04,
                                              const double * __restrict __ATTR_ALIGN__(64) pmm1r,
                                              const double * __restrict __ATTR_ALIGN__(64) pmm1i) {

                          const register __m512d tht = _mm512_load_pd(&ptht[0]);
                          const register __m512d a   = _mm512_load_pd(&pa[0]);
                          const register __m512d ka04= _mm512_load_pd(&pka04[0]);
                          const register __m512d mm1r= _mm512_load_pd(&pmm1r[0]);
                          const register __m512d mm1i= _mm512_load_pd(&pmm1i[0]);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d aa = _mm512_mul_pd(a,a);
                          const register __m512d _2 = _mm512_set1_pd(2.0f); 
                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d mms1r,mms1i,mma2r,mma2i;
                          register __m512d divr,divi,rcs,frac,cabs;
                          mms1r   = _mm512_sub_pd(mm1r,_1);
                          mma2r   = _mm512_add_pd(mm1r,_2);
                          mms1i   = _mm512_sub_pd(mm1i,_1);
                          mma2i   = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm512_mul_pd(_4pi,aa);
                          cabs    = cabs_zmm8r8(divr,divi);
                          rcs     = _mm512_mul_pd(ka04,_mm512_mul_pd(frac,cabs));
                                                  
                          return (rcs); 
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3311_zmm8r8_u(const double * __restrict  ptht,
                                              const double * __restrict  pa,
                                              const double * __restrict  pka04,
                                              const double * __restrict  pmm1r,
                                              const double * __restrict  pmm1i) {

                          const register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                          const register __m512d a   = _mm512_loadu_pd(&pa[0]);
                          const register __m512d ka04= _mm512_loadu_pd(&pka04[0]);
                          const register __m512d mm1r= _mm512_loadu_pd(&pmm1r[0]);
                          const register __m512d mm1i= _mm512_loadu_pd(&pmm1i[0]);
                          const register __m512d _1 = _mm512_set1_pd(1.0f);
                          const register __m512d aa = _mm512_mul_pd(a,a);
                          const register __m512d _2 = _mm512_set1_pd(2.0f); 
                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d mms1r,mms1i,mma2r,mma2i;
                          register __m512d divr,divi,rcs,frac,cabs;
                          mms1r   = _mm512_sub_pd(mm1r,_1);
                          mma2r   = _mm512_add_pd(mm1r,_2);
                          mms1i   = _mm512_sub_pd(mm1i,_1);
                          mma2i   = _mm512_add_pd(mm1i,_2);
                          cdiv_zmm8r8(mms1r,mms1i,mma2r,mma2i,&divr,&divi);
                          frac    = _mm512_mul_pd(_4pi,aa);
                          cabs    = cabs_zmm8r8(divr,divi);
                          rcs     = _mm512_mul_pd(ka04,_mm512_mul_pd(frac,cabs));
                                                  
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
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3314_zmm8r8(const __m512d * __restrict Snthr,
                                            const __m512d * __restrict Snthi,
                                            const __m512d * __restrict cjphr,
                                            const __m512d * __restrict cjphi,
                                            __m512d * __restrict wrkr,
                                            __m512d * __restrict wrki,
                                            const __m512d k02,
                                            const int32_t N) {

                          const register __m512d _4pi = _mm512_set1_pd(12.566370614359172953850573533118);
                          register __m512d tmpre,tmpim;
                          register __m512d rcs,cabs,frac;
                          frac = _mm512_div_pd(_4pi,k02);
                          accre = _mm512_setzero_pd();
                          accim = _mm512_setzero_pd();
                          for(int32_t i = 0; i != N; ++i) {
                              cmul_zmm8r8(Snthr[i],Snthi[i],cjphr[i],cjphi[i], 
                                           &wrkr[i],&wrki[i]);
                                                   
                          }
                          for(int32_t i = 0; i != N-1; ++i) {
                               tmpre = _mm512_add_pd(wrkr[i],wrkr[i+1]);
                               accre = _mm512_add_pd(accre,tmpre); 
                               tmpim = _mm512_add_pd(wrki[i],wrki[i+1]);
                               accim = _mm512_add_pd(accim,tmpim);
                          }
                          cabs = cabs_zmm8r8(accre,accim);
                          rcs  = _mm512_mul_pd(frac,cabs);
                          return (rcs);
                }


                  /*
                         Large sphere limit, k0a > 1.15/m1 (reflective region).
                         Backscattering RCS, formula 3.3-17
                    */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3317_zmm8r8(const __m512d m1r,
                                            const __m512d m1i,
                                            const __m512d a) {

                          const register __m512d _1   = _mm512_set1_pd(1.0f);
                          const register __m512d frac = _mm512_mul_pd(PI,
                                                                 _mm512_mul_pd(a,a));
                          register __m512d divr,divi,m1s1r,m1s1i;
                          register __m512d m1a1r,m1a1i,cabs,rcs;
                          m1s1r = _mm512_sub_pd(m1r,_1);
                          m1a1r = _mm512_add_pd(m1r,_1);
                          m1s1i = _mm512_sub_pd(m1i,_1);
                          m1a1i = _mm512_add_pd(m1i,_1);
                          cdiv_zmm8r8(m1s1r,m1s1i,m1a1r,m1a1i,&divr,&divi); 
                          cabs = _zmm8r8(divr,divi);
                          rcs  = _mm512_mul_pd(frac,cabs);
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3317_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pm1r,
                                              const double * __restrict __ATTR_ALIGN__(64) pm1i,
                                              const double * __restrict __ATTR_ALIGN__(64) pa) {

                          const register __m512d m1r  = _mm512_load_pd(&pm1r[0]);
                          const register __m512d m1i  = _mm512_load_pd(&pm1i[0]);
                          const register __m512d a    = _mm512_load_pd(&pa[0]);
                          const register __m512d _1   = _mm512_set1_pd(1.0f);
                          const register __m512d frac = _mm512_mul_pd(PI,
                                                                 _mm512_mul_pd(a,a));
                          register __m512d divr,divi,m1s1r,m1s1i;
                          register __m512d m1a1r,m1a1i,cabs,rcs;
                          m1s1r = _mm512_sub_pd(m1r,_1);
                          m1a1r = _mm512_add_pd(m1r,_1);
                          m1s1i = _mm512_sub_pd(m1i,_1);
                          m1a1i = _mm512_add_pd(m1i,_1);
                          cdiv_zmm8r8(m1s1r,m1s1i,m1a1r,m1a1i,&divr,&divi); 
                          cabs = _zmm8r8(divr,divi);
                          rcs  = _mm512_mul_pd(frac,cabs);
                          return (rcs);
                  }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3317_zmm8r8_u(const double * __restrict  pm1r,
                                              const double * __restrict  pm1i,
                                              const double * __restrict  pa) {

                          const register __m512d m1r  = _mm512_loadu_pd(&pm1r[0]);
                          const register __m512d m1i  = _mm512_loadu_pd(&pm1i[0]);
                          const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          const register __m512d _1   = _mm512_set1_pd(1.0f);
                          const register __m512d frac = _mm512_mul_pd(PI,
                                                                 _mm512_mul_pd(a,a));
                          register __m512d divr,divi,m1s1r,m1s1i;
                          register __m512d m1a1r,m1a1i,cabs,rcs;
                          m1s1r = _mm512_sub_pd(m1r,_1);
                          m1a1r = _mm512_add_pd(m1r,_1);
                          m1s1i = _mm512_sub_pd(m1i,_1);
                          m1a1i = _mm512_add_pd(m1i,_1);
                          cdiv_zmm8r8(m1s1r,m1s1i,m1a1r,m1a1i,&divr,&divi); 
                          cabs = _zmm8r8(divr,divi);
                          rcs  = _mm512_mul_pd(frac,cabs);
                          return (rcs);
                  }


                    /*
                       Forward scattering RCS.
                       Formula 3.3-19
                         */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3319_zmm8r8(const __m512d a,
                                            const __m512d k0a) {

                          const register __m512d aa   = _mm512_mul_pd(a,a);
                          const register __m512d c0   = _mm512_set1_pd(2.25);
                          const register __m512d k0a23= _mm512_pow_pd(k0a,
                                                                 _mm512_set1_pd(0.666666666666666666666666666667);
                          register __m512d rcs,fac;
                          fac = _mm512_mul_pd(PI,aa);
                          rcs = _mm512_mul_pd(fac,_mm512_mul_pd(c0,k0a23));
                          return (rcs);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3319_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                            const double * __restrict __ATTR_ALIGN__(64) pk0a) {

                          const register __m512d a    = _mm512_load_pd(&pa[0]);
                          const register __m512d aa   = _mm512_mul_pd(a,a);
                          const register __m512d k0a  = _mm512_load_pd(&pk0a[0]);
                          const register __m512d c0   = _mm512_set1_pd(2.25);
                          const register __m512d k0a23= _mm512_pow_pd(k0a,
                                                                 _mm512_set1_pd(0.666666666666666666666666666667);
                          register __m512d rcs,fac;
                          fac = _mm512_mul_pd(PI,aa);
                          rcs = _mm512_mul_pd(fac,_mm512_mul_pd(c0,k0a23));
                          return (rcs);
                 }


                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512d rcs_f3319_zmm8r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a) {

                          const register __m512d a    = _mm512_loadu_pd(&pa[0]);
                          const register __m512d aa   = _mm512_mul_pd(a,a);
                          const register __m512d k0a  = _mm512_loadu_pd(&pk0a[0]);
                          const register __m512d c0   = _mm512_set1_pd(2.25);
                          const register __m512d k0a23= _mm512_pow_pd(k0a,
                                                                 _mm512_set1_pd(0.666666666666666666666666666667);
                          register __m512d rcs,fac;
                          fac = _mm512_mul_pd(PI,aa);
                          rcs = _mm512_mul_pd(fac,_mm512_mul_pd(c0,k0a23));
                          return (rcs);
                 }


                    /*
                         Approximate solutions for far-field region (Rayleigh-Gans)
                         (abs(m1-1) << 1,2*k0a abs(m1-1) << 1)
                         Bistatic scattering formula 3.3-22
                     */

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3322_zmm8r8(const __m512d m1r,
                                         const __m512d m1i,
                                         const __m512d tht,
                                         const __m512d k0a,
                                         __m512d * __restrict S1r,
                                         __m512d * __restrict S1i) {
                       
                        register __m512d cosht,m1s1r,m1s1i,htht,cost,t0,t1;
                        register __m512d st,ct,carg,carg2,facr,faci,sinc;
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        cost  = xcos(tht);
                        const register __m512d k0a3 = _mm512_mul_pd(k0a,
                                                           _mm512_mul_pd(k0a,k0a));
                        m1s1r = _mm512_sub_pd(m1r,_1);
                        m1s1i = _mm512_sub_pd(m1i,_1);
                        const register __m512d _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d _2k0a3 = _mm512_add_pd(k0a3,k0a3)
                        htht = _mm512_mul_pd(_mm512_set1_pd(0.5f),tht);
                        cosht= xcos(htht);
                        facr = _mm512_mul_pd(_2k0a3,m1s1r);
                        faci = _mm512_mul_pd(_2k0a3,m1s1i);
                        carg = _mm512_mul_pd(_2k0a,cosht);
                        carg2= _mm512_mul_pd(carg,carg);
                        st   = xsin(carg);
                        sinc = _mm512_div_pd(st,carg);
                        ct   = xcos(carg);
                        t0   = _mm512_mul_pd(_mm512_sub_pd(st,ct),cost);
                        t1   = _mm512_div_pd(t0,carg2);
                        *S1r = _mm512_mul_pd(facr,t1);
                        *S1i = _mm512_mul_pd(faci,t1);
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3322_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pm1r,
                                           const double * __restrict __ATTR_ALIGN__(64) pm1i,
                                           const double * __restrict __ATTR_ALIGN__(64) ptht,
                                           const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                           double * __restrict __ATTR_ALIGN__(64) S1r,
                                           double * __restrict __ATTR_ALIGN__(64) S1i) {
                       
                        const register __m512d m1r = _mm512_load_pd(&pm1r[0]);
                        const register __m512d m1i = _mm512_load_pd(&pm1i[0]);
                        const register __m512d tht = _mm512_load_pd(&ptht[0]);
                        const register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d cosht,m1s1r,m1s1i,htht,cost,t0,t1;
                        register __m512d st,ct,carg,carg2,facr,faci,sinc;
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        cost  = xcos(tht);
                        const register __m512d k0a3 = _mm512_mul_pd(k0a,
                                                           _mm512_mul_pd(k0a,k0a));
                        m1s1r = _mm512_sub_pd(m1r,_1);
                        m1s1i = _mm512_sub_pd(m1i,_1);
                        const register __m512d _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d _2k0a3 = _mm512_add_pd(k0a3,k0a3)
                        htht = _mm512_mul_pd(_mm512_set1_pd(0.5f),tht);
                        cosht= xcos(htht);
                        facr = _mm512_mul_pd(_2k0a3,m1s1r);
                        faci = _mm512_mul_pd(_2k0a3,m1s1i);
                        carg = _mm512_mul_pd(_2k0a,cosht);
                        carg2= _mm512_mul_pd(carg,carg);
                        st   = xsin(carg);
                        sinc = _mm512_div_pd(st,carg);
                        ct   = xcos(carg);
                        t0   = _mm512_mul_pd(_mm512_sub_pd(st,ct),cost);
                        t1   = _mm512_div_pd(t0,carg2);
                        _mm512_store_pd(&S1r[0], _mm512_mul_pd(facr,t1));
                        _mm512_store_pd(&S1i[0], _mm512_mul_pd(faci,t1));
                 }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3322_zmm8r8_u(const double * __restrict  pm1r,
                                           const double * __restrict  pm1i,
                                           const double * __restrict  ptht,
                                           const double * __restrict  pk0a,
                                           double * __restrict S1r,
                                           double * __restrict  S1i) {
                       
                        const register __m512d m1r = _mm512_loadu_pd(&pm1r[0]);
                        const register __m512d m1i = _mm512_loadu_pd(&pm1i[0]);
                        const register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        const register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d cosht,m1s1r,m1s1i,htht,cost,t0,t1;
                        register __m512d st,ct,carg,carg2,facr,faci,sinc;
                        const register __m512d _1   = _mm512_set1_pd(1.0f);
                        cost  = xcos(tht);
                        const register __m512d k0a3 = _mm512_mul_pd(k0a,
                                                           _mm512_mul_pd(k0a,k0a));
                        m1s1r = _mm512_sub_pd(m1r,_1);
                        m1s1i = _mm512_sub_pd(m1i,_1);
                        const register __m512d _2k0a  = _mm512_add_pd(k0a,k0a);
                        const register __m512d _2k0a3 = _mm512_add_pd(k0a3,k0a3)
                        htht = _mm512_mul_pd(_mm512_set1_pd(0.5f),tht);
                        cosht= xcos(htht);
                        facr = _mm512_mul_pd(_2k0a3,m1s1r);
                        faci = _mm512_mul_pd(_2k0a3,m1s1i);
                        carg = _mm512_mul_pd(_2k0a,cosht);
                        carg2= _mm512_mul_pd(carg,carg);
                        st   = xsin(carg);
                        sinc = _mm512_div_pd(st,carg);
                        ct   = xcos(carg);
                        t0   = _mm512_mul_pd(_mm512_sub_pd(st,ct),cost);
                        t1   = _mm512_div_pd(t0,carg2);
                        _mm512_storeu_pd(&S1r[0], _mm512_mul_pd(facr,t1));
                        _mm512_storeu_pd(&S1i[0], _mm512_mul_pd(faci,t1));
                 }


                  

     } // radiolocation

} // gms











#endif /*__GMS_RCS_SPHERE_ZMM8R8_HPP__*/
