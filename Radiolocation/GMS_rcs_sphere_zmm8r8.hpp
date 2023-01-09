
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

                   const  static __m512d Ir  = _mm512_setzero_pd();
                   const  static __m512d Ii  = _mm512_set1_pd(1.0f);
                   const  static __m512d nIr = _mm512_set1_pd(-0.0f);
                   const  static __m512d nIi = _mm512_set1_pd(-1.0f);
                   const  static __m512d PI  = _mm512_set1_pd(3.14159265358979323846264338328f);

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


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm8r8_u(const double * __restrict  pk0a, // size of sphere expressed in wavenumber units
                                               __m512d * __restrict F0r, // the results
                                               __m512d * __restrict F0i) { // the results

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
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
                      t1r                       = _mm512_set1_pd(-0.0); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_pd(pn16,_mm512_set1_pd(-3.14159265358979323846264338328)); //-//-//
                      const register __m512d c5r = _mm512_set1_pd(0.444477);
                      const register __m512d c5i = _mm512_set1_pd(0.256619);
                      t2r                       = _mm512_fmadd_pd(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_pd(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512d c6r = _mm512_set1_pd(0.798821f;
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
                      t3r                       = _mm512_fmadd_pd(c3r,pn23,_mm512_set1_pd(0.695864f); //0.695864 + (0.964654 + /1.670829) at 2/3
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
                      const register __m512d c0r = _mm512_set1_pd(0.741196f);
                      const register __m512d c0i = _mm512_set1_pd(1.283788f);
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
                      const register __m512d c3i = _mm512_set1_pd(1.670829f;
                      pn13                      = _mm512_rcp14_pd(p13);
                      const register __m512d c4r = _mm512_set1_pd(7.014224);
                      const register __m512d c4i = _mm512_set1_pd(-4.049663);
                      t1r                       = _mm512_set1_pd(-0.0); // e^ipi(x-1/6)
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
                      t4i                       = _mm512_fmadd_pd(c4i,pn23,_mm512_set1_pd(0.807104f); //--//--//
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
                   void S1_f3213_zmm8r8( const __m512d k0a,
                                           const __m512d tht,
                                           __m512d * __restrict S1r,
                                           __m512d * __restrict S1i) {

                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(0.244444444444444444444444444444);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.083333333333333333333333333333f;
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/120.0f;
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(1807.0/70.0);
                        t0                         = _mm512_add_pd(half,cost);  // 0.5+cost
                        const register __m512d c5   = _mm512_set1_pd(2531.0/105.0);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmadd_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0/42.0);
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
                        *S1i                       = t6;
                        t4                         = _mm512_sub_pd(t0,_mm512_add_pd(t2,t3));
                        *S1r                       = _mm512_mul_pd(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3213_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht,
                                             double * __restrict __ATTR_ALIGN__(64) S1r,
                                             double * __restrict __ATTR_ALIGN__(64) S1i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333f;
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(0.244444444444444444444444444444);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.083333333333333333333333333333);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/120.0);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(1807.0/70.0);
                        t0                         = _mm512_add_pd(half,cost);  // 0.5+cost
                        const register __m512d c5   = _mm512_set1_pd(2531.0/105.0);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmadd_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0/42.0);
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
                   void S1_f3213_zmm8r8_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict S1r,
                                             double * __restrict S1i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d half = _mm512_set1_pd(0.5f;
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(0.244444444444444444444444444444);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.083333333333333333333333333333f;
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/120.0);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(1807.0/70.0);
                        t0                         = _mm512_add_pd(half,cost);  // 0.5+cost
                        const register __m512d c5   = _mm512_set1_pd(2531.0/105.0);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmadd_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0/42.0);
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
                   void S2_f3214_zmm8r8( const __m512d k0a,
                                           const __m512d tht,
                                           __m512d * __restrict S2r,
                                           __m512d * __restrict S2i) {

                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d c9   = _mm512_set1_pd(0.125);
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(23.0/60.0f;
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.055555555555555555555555555556);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/60.0);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(-1343.0/150.0f;
                        t0                         = _mm512_fmadd_pd(half,cost,_1);  // 0.5*cost+1
                        const register __m512d c5   = _mm512_set1_pd(3769.0/280.0);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmsub_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0/63.0);
                        t2                         = _mm512_add_pd(c4,mm512_fmadd_pd(c5,cost,
                                                                        _mm512_fmadd_pd(c6,cos2t,
                                                                                          _mm512_mul_pd(c9,cos3t))));
                        t3                         = _mm512_mul_pd(c3,_mm512_mul_pd(t2,k0a4));
                        const register __m512d c7   = _mm512_set1_pd(0.166666666666666666666666666667);
                        const register __m512d c8   = _mm512_set1_pd(0.2);
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
                   void S2_f3214_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(64) ptht,
                                             double * __restrict __ATTR_ALIGN__(64) S2r,
                                             double * __restrict __ATTR_ALIGN__(64) S2i) {

                        register __m512d k0a = _mm512_load_pd(&pk0a[0]);
                        register __m512d tht = _mm512_load_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d c9   = _mm512_set1_pd(0.125);
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333f;
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(23.0/60.0);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.055555555555555555555555555556);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/60.0);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(-1343.0/150.0);
                        t0                         = _mm512_fmadd_pd(half,cost,_1);  // 0.5*cost+1
                        const register __m512d c5   = _mm512_set1_pd(3769.0/280.0);
                        t1                         = _mm512_sub_pd(c0,_mm512_fmsub_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0/63.0);
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
                   void S2_f3214_zmm8r8_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict  S2r,
                                             double * __restrict  S2i) {

                        register __m512d k0a = _mm512_loadu_pd(&pk0a[0]);
                        register __m512d tht = _mm512_loadu_pd(&ptht[0]);
                        register __m512d cost,cos2t,cos3t;
                        register __m512d k0a3,k0a2,k0a4,k0a6;
                        register __m512d t0,t1,t2,t3,t4,t5,t6;
                        const register __m512d c9   = _mm512_set1_pd(0.125);
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
                        cost                       = xcos(tht);
                        const register __m512d _2   = _mm512_set1_pd(2.0);
                        cos2t                      = xcos(_mm512_add_pd(tht,tht));
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        cos3t                      = xcos(_mm512_add_pd(tht,_mm512_add_pd(tht,tht)));
                        const register __m512d c0   = _mm512_set1_pd(0.3333333333333333333333333333333);
                        k0a2                       = _mm512_mul_pd(k0a,k0a);
                        const register __m512d c1   = _mm512_set1_pd(23.0/60.0);
                        k0a3                       = _mm512_mul_pd(k0a2,k0a);
                        const register __m512d c2   = _mm512_set1_pd(0.055555555555555555555555555556);
                        k0a4                       = _mm512_mul_pd(k0a2,k0a2);
                        const register __m512d c3   = _mm512_set1_pd(1.0/60.0);
                        k0a6                       = _mm512_mul_pd(k0a3,k0a2);
                        const register __m512d c4   = _mm512_set1_pd(-1343.0/150.0);
                        t0                         = _mm512_fmadd_pd(half,cost,_1);  // 0.5*cost+1
                        const register __m512d c5   = _mm512_set1_pd(3769.0/280.0f;
                        t1                         = _mm512_sub_pd(c0,_mm512_fmsub_pd(c1,cost,
                                                                            _mm512_mul_pd(c2,cos2t)));
                        t1                         = _mm512_mul_pd(t1,k0a2);
                        const register __m512d c6   = _mm512_set1_pd(57.0/63.0);
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

                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        const register __m512d _7   = _mm512_set1_pd(7.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
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
                        carr  = _mm512_set1_pd(-0.0);
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
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        const register __m512d _7   = _mm512_set1_pd(7.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
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
                        carr  = _mm512_set1_pd(-0.0);
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
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d _4   = _mm512_set1_pd(4.0);
                        const register __m512d _7   = _mm512_set1_pd(7.0f;
                        const register __m512d half = _mm512_set1_pd(0.5f;
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

                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d quat = _mm512_set1_pd(0.25);
                        const register __m512d _6   = _mm512_set1_pd(6.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
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
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d quat = _mm512_set1_pd(0.25);
                        const register __m512d _6   = _mm512_set1_pd(6.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
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
                        carr  = _mm512_set1_pd(-0.0);
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
                        const register __m512d _1   = _mm512_set1_pd(1.0);
                        const register __m512d quat = _mm512_set1_pd(0.25);
                        const register __m512d _6   = _mm512_set1_pd(6.0);
                        const register __m512d half = _mm512_set1_pd(0.5);
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

                       Creeping wave contribution equations.
                       Eqs. (3.2-18) and (3.2-19) are not valid near either the backscattering or
                       forward scattering directions, but are valid in an intermediate region
                    */

                 


     } // radiolocation

} // gms











#endif /*__GMS_RCS_SPHERE_ZMM8R8_HPP__*/
