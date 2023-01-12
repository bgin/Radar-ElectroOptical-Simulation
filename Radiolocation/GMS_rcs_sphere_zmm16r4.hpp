
#ifndef __GMS_RCS_SPHERE_ZMM16R4_HPP__
#define __GMS_RCS_SPHERE_ZMM16R4_HPP__


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

    const unsigned int GMS_RCS_SPHERE_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_SPHERE_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_SPHERE_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_SPHERE_ZMM16R4_FULLVER =
      1000U*GMS_RCS_SPHERE_ZMM16R4_MAJOR+
      100U*GMS_RCS_SPHERE_ZMM16R4_MINOR+
      10U*GMS_RCS_SPHERE_ZMM16R4_MICRO;
    const char * const GMS_RCS_SPHERE_ZMM16R4_CREATION_DATE = "04-01-2023 12:45 AM +00200 (WED 04 JAN 2023 GMT+2)";
    const char * const GMS_RCS_SPHERE_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_SPHERE_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_SPHERE_ZMM16R4_DESCRIPTION   = "AVX512 optimized Sphere Radar Cross Section (analytic) functionality."

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

                   const  static __m512 Ir  = _mm512_setzero_ps();
                   const  static __m512 Ii  = _mm512_set1_ps(1.0f);
                   const  static __m512 nIr = _mm512_set1_ps(-0.0f);
                   const  static __m512 nIi = _mm512_set1_ps(-1.0f);
                   const  static __m512 PI  = _mm512_set1_ps(3.14159265358979323846264338328f);

                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm16r4(const __m512 k0a, // size of sphere expressed in wavenumber units
                                               __m512 * __restrict F0r, // the results
                                               __m512 * __restrict F0i) { // the results

                        const register __m512 _1 = _mm512_set1_ps(1.0f);
                        const register __m512 c0 = _mm512_set1_ps(1.5f);
                        const register __m512 c1 = _mm512_set1_ps(0.092592592592592592592592592593f);
                        const register __m512 c2 = _mm512_set1_ps(0.018888888888888888888888888889f);
                        const register __m512 c3 = _mm512_set1_ps(0.558656504577139497774418409339f);
                        const register __m512 c4 = _mm512_set1_ps(1.2f);
                        const register __m512 c5 = _mm512_set1_ps(0.5f);
                        __m512 k0a2,k0a3,k0a4,k0a6;
                        __m512 t0,t1,t2,t3;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        t1   = _mm512_sub_ps(_1,_mm512_mul_ps(c1,k0a2));
                        k0a3 = _mm512_mul_ps(k0a2,k0a);
                        t0   = _mm512_mul_ps(c0,k0a3);
                        k0a4 = _mm512_mul_ps(k0a3,k0a);
                        k0a6 = _mm512_mul_ps(k0a3,k0a2);
                        t2   = _mm512_fmsub_ps(c2,k0a4,_mm512_mul_ps(c3,k0a6));
                        t3   = _mm512_mul_ps(t0,_mm512_add_ps(t1,t1));
                        *F0r = t3;
                        t0   = _mm512_mul_ps(c5,k0a6);
                        t1   = _mm512_add_ps(_1,_mm512_mul_ps(c4,k0a2));
                        *F0i = _mm512_mul_ps(t0,t1);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a, // size of sphere expressed in wavenumber units
                                               __m512 * __restrict F0r, // the results
                                               __m512 * __restrict F0i) { // the results

                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        const register __m512 _1 = _mm512_set1_ps(1.0f);
                        const register __m512 c0 = _mm512_set1_ps(1.5f);
                        const register __m512 c1 = _mm512_set1_ps(0.092592592592592592592592592593f);
                        const register __m512 c2 = _mm512_set1_ps(0.018888888888888888888888888889f);
                        const register __m512 c3 = _mm512_set1_ps(0.558656504577139497774418409339f);
                        const register __m512 c4 = _mm512_set1_ps(1.2f);
                        const register __m512 c5 = _mm512_set1_ps(0.5f);
                        __m512 k0a2,k0a3,k0a4,k0a6;
                        __m512 t0,t1,t2,t3;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        t1   = _mm512_sub_ps(_1,_mm512_mul_ps(c1,k0a2));
                        k0a3 = _mm512_mul_ps(k0a2,k0a);
                        t0   = _mm512_mul_ps(c0,k0a3);
                        k0a4 = _mm512_mul_ps(k0a3,k0a);
                        k0a6 = _mm512_mul_ps(k0a3,k0a2);
                        t2   = _mm512_fmsub_ps(c2,k0a4,_mm512_mul_ps(c3,k0a6));
                        t3   = _mm512_mul_ps(t0,_mm512_add_ps(t1,t1));
                        *F0r = t3;
                        t0   = _mm512_mul_ps(c5,k0a6);
                        t1   = _mm512_add_ps(_1,_mm512_mul_ps(c4,k0a2));
                        *F0i = _mm512_mul_ps(t0,t1);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_func_324_zmm16r4_u(const float * __restrict  pk0a, // size of sphere expressed in wavenumber units
                                               __m512 * __restrict F0r, // the results
                                               __m512 * __restrict F0i) { // the results

                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 _1 = _mm512_set1_ps(1.0f);
                        const register __m512 c0 = _mm512_set1_ps(1.5f);
                        const register __m512 c1 = _mm512_set1_ps(0.092592592592592592592592592593f);
                        const register __m512 c2 = _mm512_set1_ps(0.018888888888888888888888888889f);
                        const register __m512 c3 = _mm512_set1_ps(0.558656504577139497774418409339f);
                        const register __m512 c4 = _mm512_set1_ps(1.2f);
                        const register __m512 c5 = _mm512_set1_ps(0.5f);
                        __m512 k0a2,k0a3,k0a4,k0a6;
                        __m512 t0,t1,t2,t3;
                        k0a2 = _mm512_mul_ps(k0a,k0a);
                        t1   = _mm512_sub_ps(_1,_mm512_mul_ps(c1,k0a2));
                        k0a3 = _mm512_mul_ps(k0a2,k0a);
                        t0   = _mm512_mul_ps(c0,k0a3);
                        k0a4 = _mm512_mul_ps(k0a3,k0a);
                        k0a6 = _mm512_mul_ps(k0a3,k0a2);
                        t2   = _mm512_fmsub_ps(c2,k0a4,_mm512_mul_ps(c3,k0a6));
                        t3   = _mm512_mul_ps(t0,_mm512_add_ps(t1,t1));
                        *F0r = t3;
                        t0   = _mm512_mul_ps(c5,k0a6);
                        t1   = _mm512_add_ps(_1,_mm512_mul_ps(c4,k0a2));
                        *F0i = _mm512_mul_ps(t0,t1);
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
                   __m512 rcs_f325_zmm16r4(const __m512 k0,
                                           const __m512 a ) {
                                        
                        register __m512 a2,k0a,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,sigma; 
                        a2  = _mm512_mul_ps(a,a);
                        const register __m512 pi9 = _mm512_set1_ps(28.274333882308139146163790449516);
                        k0a = _mm512_mul_ps(k0,a); 
                        const register __m512 _1 = _mm512_set1_ps(1.0f);
                        k0a2= _mm512_mul_ps(k0a,k0a);
                        const register __m512 c0 = _mm512_set1_ps(0.185185185185185185185185185185f);
                        k0a4= _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c1 = _mm512_set1_ps(0.04635116598079561042524005487f);
                        k0a6= _mm512_mul_ps(k0a4,_mm512_mul_ps(k0a,k0a));
                        const register __m512 c2 = _mm512_set1_ps(1.00969984042999916015789031662f);
                        t0  = _mm512_mul_ps(k0a4,_mm512_mul_ps(pi9,a2));
                        t1  = _mm512_sub_ps(_1,_mm512_mul_ps(c0,k0a2));
                        t2  = _mm512_fmsub_ps(c1,k0a4,_mm512_mul_ps(c2,k0a6));
                        sigma = _mm512_fmadd_ps(t1,t2,t0);
                        return (sigma);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f325_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pa ) {
                        
                        register __m512 k0= _mm512_load_ps(&pk0[0]); 
                        register __m512 a = _mm512_load_ps(&pa[0]);               
                        register __m512 a2,k0a,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,sigma; 
                        a2  = _mm512_mul_ps(a,a);
                        const register __m512 pi9 = _mm512_set1_ps(28.274333882308139146163790449516);
                        k0a = _mm512_mul_ps(k0,a); 
                        const register __m512 _1 = _mm512_set1_ps(1.0f);
                        k0a2= _mm512_mul_ps(k0a,k0a);
                        const register __m512 c0 = _mm512_set1_ps(0.185185185185185185185185185185f);
                        k0a4= _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c1 = _mm512_set1_ps(0.04635116598079561042524005487f);
                        k0a6= _mm512_mul_ps(k0a4,_mm512_mul_ps(k0a,k0a));
                        const register __m512 c2 = _mm512_set1_ps(1.00969984042999916015789031662f);
                        t0  = _mm512_mul_ps(k0a4,_mm512_mul_ps(pi9,a2));
                        t1  = _mm512_sub_ps(_1,_mm512_mul_ps(c0,k0a2));
                        t2  = _mm512_fmsub_ps(c1,k0a4,_mm512_mul_ps(c2,k0a6));
                        sigma = _mm512_fmadd_ps(t1,t2,t0);
                        return (sigma);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f325_zmm16r4_u(const float * __restrict pk0,
                                             const float * __restrict pa ) {
                        
                        register __m512 k0= _mm512_loadu_ps(&pk0[0]); 
                        register __m512 a = _mm512_loadu_ps(&pa[0]);               
                        register __m512 a2,k0a,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,sigma; 
                        a2  = _mm512_mul_ps(a,a);
                        const register __m512 pi9 = _mm512_set1_ps(28.274333882308139146163790449516);
                        k0a = _mm512_mul_ps(k0,a); 
                        const register __m512 _1 = _mm512_set1_ps(1.0f);
                        k0a2= _mm512_mul_ps(k0a,k0a);
                        const register __m512 c0 = _mm512_set1_ps(0.185185185185185185185185185185f);
                        k0a4= _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c1 = _mm512_set1_ps(0.04635116598079561042524005487f);
                        k0a6= _mm512_mul_ps(k0a4,_mm512_mul_ps(k0a,k0a));
                        const register __m512 c2 = _mm512_set1_ps(1.00969984042999916015789031662f);
                        t0  = _mm512_mul_ps(k0a4,_mm512_mul_ps(pi9,a2));
                        t1  = _mm512_sub_ps(_1,_mm512_mul_ps(c0,k0a2));
                        t2  = _mm512_fmsub_ps(c1,k0a4,_mm512_mul_ps(c2,k0a6));
                        sigma = _mm512_fmadd_ps(t1,t2,t0);
                        return (sigma);
                }



#include "GMS_complex_zmm16r4.hpp"


                  /*
                        Creeping wave term, F c(0) at the upper end of the resonance region.
                        Formula 3.2-8
                    */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_f328_zmm16r4(const __m512 x,//k0a
                                         __m512 * __restrict Fc0r,
                                         __m512 * __restrict Fc0i) {

                      register __m512 p43,pn16,pn13,pn23,p13;
                      register __m512 t0r,t0i,t1r,t1i,t2r,t2i;  
                      register __m512 t3r,t3i,t4r,t4i,t5r,t5i;
                      register __m512 t6r,t6i,e1r,e1i,e2r,e2i;
                      register __m512 e3r,e3i,t6r,t6i,t7r,t7i;
                      register __m512 argr,argi,resr,resi;              
                      const register __m512 c0  = _mm512_set1_ps(1.357588f);
                      p43                       = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(1.3333333333333333333333333333333f));
                      t0r                       = _mm512_mul_ps(nIr,p43); // - / x 4/3
                      t0i                       = _mm512_mul_ps(nIi,p43); // - / x 4/3
                      const register __m512 c0r = _mm512_set1_ps(0.741196f);
                      const register __m512 c0i = _mm512_set1_ps(1.283788f);
                      pn16                      = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(0.166666666666666666666666666667f));
                      pn16                      = _mm512_rcp14_ps(pn16);
                      const register __m512 c1r = _mm512_set1_ps(2.200002f);
                      const register __m512 c1i = _mm512_set1_ps(-1.270172f);
                      pn23                      = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(0.666666666666666666666666666667f));
                      pn23                      = _mm512_rcp14_ps(pn23);
                      const register __m512 c2r = _mm512_set1_ps(0.445396f);
                      const register __m512 c2i = _mm512_set1_ps(0.257150f);
                      p13                       = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(-0.3333333333333333333333333333333333f));
                      const register __m512 c3r = _mm512_set1_ps(0.964654f);
                      const register __m512 c3i = _mm512_set1_ps(1.670829f);
                      pn13                      = _mm512_rcp14_ps(p13);
                      const register __m512 c4r = _mm512_set1_ps(7.014224f);
                      const register __m512 c4i = _mm512_set1_ps(-4.049663f);
                      t1r                       = _mm512_set1_ps(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_ps(pn16,_mm512_set1_ps(-3.14159265358979323846264338328f)); //-//-//
                      const register __m512 c5r = _mm512_set1_ps(0.444477f);
                      const register __m512 c5i = _mm512_set1_ps(0.256619f);
                      t2r                       = _mm512_fmadd_ps(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_ps(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512 c6r = _mm512_set1_ps(0.798821f);
                      const register __m512 c6i = _mm512_set1_ps(1.383598f);
                      t3r                       = _mm512_mul_ps(pn13,c1r);
                      t3i                       = _mm512_mul_ps(pn13,c1i);
                      const register __m512 c7r = _mm512_set1_ps(5.048956f);
                      const register __m512 c7i = _mm512_set1_ps(-2.915016f);
                      t4r                       = _mm512_mul_ps(pn13,c2r);
                      t4i                       = _mm512_mul_ps(pn13,c2i);
                      const register __m512 c8r = _mm512_set1_ps(0.312321f);
                      const register __m512 c8i = _mm512_set1_ps(0.180319f);
                      t5r                       = _mm512_add_ps(t3r,t4r);
                      t5i                       = _mm512_add_ps(t3i,t4i);
                      cexp_zmm16r4(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm512_fmadd_ps(c3r,pn23,_mm512_set1_ps(0.695864f)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm512_fmadd_ps(c3i,pn23,_mm512_set1_ps(0.695864f)); //--//--//--//
                      t4r                       = _mm512_mul_ps(p13,c4r);
                      t4i                       = _mm512_mul_ps(p13,c4i);
                      t5r                       = _mm512_mul_ps(pn13,c5r);
                      t5i                       = _mm512_mul_ps(pn13,c5i);
                      t6r                       = _mm512_sub_ps(t4r,t5r);
                      t6i                       = _mm512_sub_ps(t4i,t5i);
                      cexp_zmm16r4(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm512_fmadd_ps(c4r,pn23,_mm512_set1_ps(0.807104f)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm512_fmadd_ps(c4i,pn23,_mm512_set1_ps(0.807104f)); //--//--//
                      t5r                       = _mm512_mul_ps(p13,c7r);
                      t5i                       = _mm512_mul_ps(p13,c7i);
                      t6r                       = _mm512_mul_ps(pn13,c8r);
                      t6i                       = _mm512_mul_ps(pn13,c8i);
                      t7r                       = _mm512_sub_ps(t5r,t6r);
                      t7i                       = _mm512_sub_ps(t5i,t6i);
                      cexp_zmm16r4(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm512_setzero_ps();
                      t5i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm512_setzero_ps();
                      t6i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm512_setzero_ps();
                      t7i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm512_add_ps(t5r,_mm512_sub_ps(t6r,t7r));
                      argi                      = _mm512_add_ps(t5i,_mm512_sub_ps(t6i,t7i));
                      t4r                       = _mm512_setzero_ps();
                      t4i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_zmm16r4(t4r,t4i,&resr,&resi);
                      cmul_zmm16r4(t04,t0i,resr,resi,&Fc0r,&Fc0i);                    
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_f328_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) px,//k0a
                                           float * __restrict __ATTR_ALIGN__(64) Fc0r,
                                           float * __restrict __ATTR_ALIGN__(64) Fc0i) {

                      register __m512 x = _mm512_load_ps(&px[0]);
                      register __m512 p43,pn16,pn13,pn23,p13;
                      register __m512 t0r,t0i,t1r,t1i,t2r,t2i;  
                      register __m512 t3r,t3i,t4r,t4i,t5r,t5i;
                      register __m512 t6r,t6i,e1r,e1i,e2r,e2i;
                      register __m512 e3r,e3i,t6r,t6i,t7r,t7i;
                      register __m512 argr,argi,resr,resi; 
                      register __m512 fcr,fci;             
                      const register __m512 c0  = _mm512_set1_ps(1.357588f);
                      p43                       = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(1.3333333333333333333333333333333f));
                      t0r                       = _mm512_mul_ps(nIr,p43); // - / x 4/3
                      t0i                       = _mm512_mul_ps(nIi,p43); // - / x 4/3
                      const register __m512 c0r = _mm512_set1_ps(0.741196f);
                      const register __m512 c0i = _mm512_set1_ps(1.283788f);
                      pn16                      = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(0.166666666666666666666666666667f));
                      pn16                      = _mm512_rcp14_ps(pn16);
                      const register __m512 c1r = _mm512_set1_ps(2.200002f);
                      const register __m512 c1i = _mm512_set1_ps(-1.270172f);
                      pn23                      = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(0.666666666666666666666666666667f));
                      pn23                      = _mm512_rcp14_ps(pn23);
                      const register __m512 c2r = _mm512_set1_ps(0.445396f);
                      const register __m512 c2i = _mm512_set1_ps(0.257150f);
                      p13                       = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(-0.3333333333333333333333333333333333f));
                      const register __m512 c3r = _mm512_set1_ps(0.964654f);
                      const register __m512 c3i = _mm512_set1_ps(1.670829f);
                      pn13                      = _mm512_rcp14_ps(p13);
                      const register __m512 c4r = _mm512_set1_ps(7.014224f);
                      const register __m512 c4i = _mm512_set1_ps(-4.049663f);
                      t1r                       = _mm512_set1_ps(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_ps(pn16,_mm512_set1_ps(-3.14159265358979323846264338328f)); //-//-//
                      const register __m512 c5r = _mm512_set1_ps(0.444477f);
                      const register __m512 c5i = _mm512_set1_ps(0.256619f);
                      t2r                       = _mm512_fmadd_ps(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_ps(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512 c6r = _mm512_set1_ps(0.798821f);
                      const register __m512 c6i = _mm512_set1_ps(1.383598f);
                      t3r                       = _mm512_mul_ps(pn13,c1r);
                      t3i                       = _mm512_mul_ps(pn13,c1i);
                      const register __m512 c7r = _mm512_set1_ps(5.048956f);
                      const register __m512 c7i = _mm512_set1_ps(-2.915016f);
                      t4r                       = _mm512_mul_ps(pn13,c2r);
                      t4i                       = _mm512_mul_ps(pn13,c2i);
                      const register __m512 c8r = _mm512_set1_ps(0.312321f);
                      const register __m512 c8i = _mm512_set1_ps(0.180319f);
                      t5r                       = _mm512_add_ps(t3r,t4r);
                      t5i                       = _mm512_add_ps(t3i,t4i);
                      cexp_zmm16r4(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm512_fmadd_ps(c3r,pn23,_mm512_set1_ps(0.695864f)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm512_fmadd_ps(c3i,pn23,_mm512_set1_ps(0.695864f)); //--//--//--//
                      t4r                       = _mm512_mul_ps(p13,c4r);
                      t4i                       = _mm512_mul_ps(p13,c4i);
                      t5r                       = _mm512_mul_ps(pn13,c5r);
                      t5i                       = _mm512_mul_ps(pn13,c5i);
                      t6r                       = _mm512_sub_ps(t4r,t5r);
                      t6i                       = _mm512_sub_ps(t4i,t5i);
                      cexp_zmm16r4(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm512_fmadd_ps(c4r,pn23,_mm512_set1_ps(0.807104f)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm512_fmadd_ps(c4i,pn23,_mm512_set1_ps(0.807104f)); //--//--//
                      t5r                       = _mm512_mul_ps(p13,c7r);
                      t5i                       = _mm512_mul_ps(p13,c7i);
                      t6r                       = _mm512_mul_ps(pn13,c8r);
                      t6i                       = _mm512_mul_ps(pn13,c8i);
                      t7r                       = _mm512_sub_ps(t5r,t6r);
                      t7i                       = _mm512_sub_ps(t5i,t6i);
                      cexp_zmm16r4(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm512_setzero_ps();
                      t5i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm512_setzero_ps();
                      t6i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm512_setzero_ps();
                      t7i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm512_add_ps(t5r,_mm512_sub_ps(t6r,t7r));
                      argi                      = _mm512_add_ps(t5i,_mm512_sub_ps(t6i,t7i));
                      t4r                       = _mm512_setzero_ps();
                      t4i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_zmm16r4(t4r,t4i,&resr,&resi);
                      cmul_zmm16r4(t04,t0i,resr,resi,&fcr,&fci);    
                      __m512_store_ps(&Fc0r[0], fcr);
                      __m512_store_ps(&Fc0i[0], fci);                
                }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void bsc_f328_zmm16r4_u(const float * __restrict px,//k0a
                                           float * __restrict  Fc0r,
                                           float * __restrict  Fc0i) {

                      register __m512 x = _mm512_loadu_ps(&px[0]);
                      register __m512 p43,pn16,pn13,pn23,p13;
                      register __m512 t0r,t0i,t1r,t1i,t2r,t2i;  
                      register __m512 t3r,t3i,t4r,t4i,t5r,t5i;
                      register __m512 t6r,t6i,e1r,e1i,e2r,e2i;
                      register __m512 e3r,e3i,t6r,t6i,t7r,t7i;
                      register __m512 argr,argi,resr,resi; 
                      register __m512 fcr,fci;             
                      const register __m512 c0  = _mm512_set1_ps(1.357588f);
                      p43                       = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(1.3333333333333333333333333333333f));
                      t0r                       = _mm512_mul_ps(nIr,p43); // - / x 4/3
                      t0i                       = _mm512_mul_ps(nIi,p43); // - / x 4/3
                      const register __m512 c0r = _mm512_set1_ps(0.741196f);
                      const register __m512 c0i = _mm512_set1_ps(1.283788f);
                      pn16                      = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(0.166666666666666666666666666667f));
                      pn16                      = _mm512_rcp14_ps(pn16);
                      const register __m512 c1r = _mm512_set1_ps(2.200002f);
                      const register __m512 c1i = _mm512_set1_ps(-1.270172f);
                      pn23                      = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(0.666666666666666666666666666667f));
                      pn23                      = _mm512_rcp14_ps(pn23);
                      const register __m512 c2r = _mm512_set1_ps(0.445396f);
                      const register __m512 c2i = _mm512_set1_ps(0.257150f);
                      p13                       = _mm512_pow_ps(x,
                                                        _mm512_set1_ps(-0.3333333333333333333333333333333333f));
                      const register __m512 c3r = _mm512_set1_ps(0.964654f);
                      const register __m512 c3i = _mm512_set1_ps(1.670829f);
                      pn13                      = _mm512_rcp14_ps(p13);
                      const register __m512 c4r = _mm512_set1_ps(7.014224f);
                      const register __m512 c4i = _mm512_set1_ps(-4.049663f);
                      t1r                       = _mm512_set1_ps(-0.0f); // e^ipi(x-1/6)
                      t1i                       = _mm512_mul_ps(pn16,_mm512_set1_ps(-3.14159265358979323846264338328f)); //-//-//
                      const register __m512 c5r = _mm512_set1_ps(0.444477f);
                      const register __m512 c5i = _mm512_set1_ps(0.256619f);
                      t2r                       = _mm512_fmadd_ps(c0r,pn23,c0);
                      t2i                       = _mm512_fmadd_ps(c0i,pn23,c0); // [1.357588 + (0.741196 + /1.283788) at 2/3]
                      const register __m512 c6r = _mm512_set1_ps(0.798821f);
                      const register __m512 c6i = _mm512_set1_ps(1.383598f);
                      t3r                       = _mm512_mul_ps(pn13,c1r);
                      t3i                       = _mm512_mul_ps(pn13,c1i);
                      const register __m512 c7r = _mm512_set1_ps(5.048956f);
                      const register __m512 c7i = _mm512_set1_ps(-2.915016f);
                      t4r                       = _mm512_mul_ps(pn13,c2r);
                      t4i                       = _mm512_mul_ps(pn13,c2i);
                      const register __m512 c8r = _mm512_set1_ps(0.312321f);
                      const register __m512 c8i = _mm512_set1_ps(0.180319f);
                      t5r                       = _mm512_add_ps(t3r,t4r);
                      t5i                       = _mm512_add_ps(t3i,t4i);
                      cexp_zmm16r4(t5r,t5i,&e1r,&e1i); // exp[—x1/*(2.200002 - /1.270172) + x“1/3(0.445396 + i 0.257150)]}E
                      t3r                       = _mm512_fmadd_ps(c3r,pn23,_mm512_set1_ps(0.695864f)); //0.695864 + (0.964654 + /1.670829) at 2/3
                      t3i                       = _mm512_fmadd_ps(c3i,pn23,_mm512_set1_ps(0.695864f)); //--//--//--//
                      t4r                       = _mm512_mul_ps(p13,c4r);
                      t4i                       = _mm512_mul_ps(p13,c4i);
                      t5r                       = _mm512_mul_ps(pn13,c5r);
                      t5i                       = _mm512_mul_ps(pn13,c5i);
                      t6r                       = _mm512_sub_ps(t4r,t5r);
                      t6i                       = _mm512_sub_ps(t4i,t5i);
                      cexp_zmm16r4(t6r,t6i,&e2r,&e2i); // exp[-x 1/*(7.014224 -/4.049663) - * -1/2(0.444477 + /0.256619)]
                      t4r                       = _mm512_fmadd_ps(c4r,pn23,_mm512_set1_ps(0.807104f)); // 0.807104 + (0.798821 + /1.383598) at 2/3
                      t4i                       = _mm512_fmadd_ps(c4i,pn23,_mm512_set1_ps(0.807104f)); //--//--//
                      t5r                       = _mm512_mul_ps(p13,c7r);
                      t5i                       = _mm512_mul_ps(p13,c7i);
                      t6r                       = _mm512_mul_ps(pn13,c8r);
                      t6i                       = _mm512_mul_ps(pn13,c8i);
                      t7r                       = _mm512_sub_ps(t5r,t6r);
                      t7i                       = _mm512_sub_ps(t5i,t6i);
                      cexp_zmm16r4(t7r,t7i,&e3r,&e3i); // exp[—x1/8(5.048956 - /2.915016) - Ar1/3(0.312321 + /0.180319)]
                      t5r                       = _mm512_setzero_ps();
                      t5i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t2r,t2i,e1r,e1i,&t5r,&t5i); 
                      t6r                       = _mm512_setzero_ps();
                      t6i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t3r,t3i,e2r,e2i,&t6r,&t6i);
                      t7r                       = _mm512_setzero_ps();
                      t7i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t4r,t4i,e3r,e3i,&t7r,&t7i);
                      argr                      = _mm512_add_ps(t5r,_mm512_sub_ps(t6r,t7r));
                      argi                      = _mm512_add_ps(t5i,_mm512_sub_ps(t6i,t7i));
                      t4r                       = _mm512_setzero_ps();
                      t4i                       = _mm512_setzero_ps();
                      cmul_zmm16r4(t1r,t1i,argr,argi,&t4r,&t4i);
                      cexp_zmm16r4(t4r,t4i,&resr,&resi);
                      cmul_zmm16r4(t04,t0i,resr,resi,&fcr,&fci);    
                      __m512_storeu_ps(&Fc0r[0], fcr);
                      __m512_storeu_ps(&Fc0i[0], fci);                
                }




                  

#include "GMS_sleefsimdsp.hpp"

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
                   void S1_f3213_zmm146r4( const __m512 k0a,
                                           const __m512 tht,
                                           __m512 * __restrict S1r,
                                           __m512 * __restrict S1i) {

                        register __m512 cost,cos2t,cos3t;
                        register __m512 k0a3,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,t4,t5,t6;
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cost                       = xcosf(tht);
                        const register __m512 _2   = _mm512_set1_ps(2.0f);
                        cos2t                      = xcosf(_mm512_add_ps(tht,tht));
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        cos3t                      = xcosf(_mm512_add_ps(tht,_mm512_add_ps(tht,tht)));
                        const register __m512 c0   = _mm512_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c1   = _mm512_set1_ps(0.244444444444444444444444444444f);
                        k0a3                       = _mm512_mul_ps(k0a2,k0a);
                        const register __m512 c2   = _mm512_set1_ps(0.083333333333333333333333333333f);
                        k0a4                       = _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c3   = _mm512_set1_ps(1.0f/120.0f);
                        k0a6                       = _mm512_mul_ps(k0a3,k0a2);
                        const register __m512 c4   = _mm512_set1_ps(1807.0f/70.0f);
                        t0                         = _mm512_add_ps(half,cost);  // 0.5+cost
                        const register __m512 c5   = _mm512_set1_ps(2531.0f/105.0f);
                        t1                         = _mm512_sub_ps(c0,_mm512_fmadd_ps(c1,cost,
                                                                            _mm512_mul_ps(c2,cos2t)));
                        t1                         = _mm512_mul_ps(t1,k0a2);
                        const register __m512 c6   = _mm512_set1_ps(57.0f/42.0f);
                        t2                         = _mm512_sub_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm512_fmadd_ps(c6,cos2t,
                                                                                          _mm512_mul_ps(c0,cos3t))));
                        t3                         = _mm512_mul_ps(c3,_mm512_mul_ps(t2,k0a4));
                        const register __m512 c7   = _mm512_set1_ps(0.166666666666666666666666666667f);
                        const register __m512 c8   = _mm512_set1_ps(0.2f);
                        t4                         = _mm512_mul_ps(c7,_mm512_fmsub_ps(_4,cost,_1));
                        t5                         = _mm512_mul_ps(c8,_mm512_fmadd_ps(_2,cost,_1));
                        t5                         = _mm512_mul_ps(t5,k0a2);
                        t6                         = _mm512_mul_ps(k0a6,_mm512_add_ps(t4,t5)); // imaginary part
                        *S1i                       = t6;
                        t4                         = _mm512_sub_ps(t0,_mm512_add_ps(t2,t3));
                        *S1r                       = _mm512_mul_ps(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3213_zmm146r4_a( const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht,
                                             float * __restrict __ATTR_ALIGN__(64) S1r,
                                             float * __restrict __ATTR_ALIGN__(64) S1i) {

                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 cost,cos2t,cos3t;
                        register __m512 k0a3,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,t4,t5,t6;
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cost                       = xcosf(tht);
                        const register __m512 _2   = _mm512_set1_ps(2.0f);
                        cos2t                      = xcosf(_mm512_add_ps(tht,tht));
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        cos3t                      = xcosf(_mm512_add_ps(tht,_mm512_add_ps(tht,tht)));
                        const register __m512 c0   = _mm512_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c1   = _mm512_set1_ps(0.244444444444444444444444444444f);
                        k0a3                       = _mm512_mul_ps(k0a2,k0a);
                        const register __m512 c2   = _mm512_set1_ps(0.083333333333333333333333333333f);
                        k0a4                       = _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c3   = _mm512_set1_ps(1.0f/120.0f);
                        k0a6                       = _mm512_mul_ps(k0a3,k0a2);
                        const register __m512 c4   = _mm512_set1_ps(1807.0f/70.0f);
                        t0                         = _mm512_add_ps(half,cost);  // 0.5+cost
                        const register __m512 c5   = _mm512_set1_ps(2531.0f/105.0f);
                        t1                         = _mm512_sub_ps(c0,_mm512_fmadd_ps(c1,cost,
                                                                            _mm512_mul_ps(c2,cos2t)));
                        t1                         = _mm512_mul_ps(t1,k0a2);
                        const register __m512 c6   = _mm512_set1_ps(57.0f/42.0f);
                        t2                         = _mm512_sub_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm512_fmadd_ps(c6,cos2t,
                                                                                          _mm512_mul_ps(c0,cos3t))));
                        t3                         = _mm512_mul_ps(c3,_mm512_mul_ps(t2,k0a4));
                        const register __m512 c7   = _mm512_set1_ps(0.166666666666666666666666666667f);
                        const register __m512 c8   = _mm512_set1_ps(0.2f);
                        t4                         = _mm512_mul_ps(c7,_mm512_fmsub_ps(_4,cost,_1));
                        t5                         = _mm512_mul_ps(c8,_mm512_fmadd_ps(_2,cost,_1));
                        t5                         = _mm512_mul_ps(t5,k0a2);
                        t6                         = _mm512_mul_ps(k0a6,_mm512_add_ps(t4,t5)); // imaginary part
                        _mm512_store_ps(&S1i[0],t6);
                        
                        t4                         = _mm512_sub_ps(t0,_mm512_add_ps(t2,t3));
                        _mm512_store_ps(&S1r[0],_mm512_mul_ps(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3213_zmm146r4_u( const float * __restrict  pk0a,
                                             const float * __restrict  ptht,
                                             float * __restrict S1r,
                                             float * __restrict S1i) {

                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 cost,cos2t,cos3t;
                        register __m512 k0a3,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,t4,t5,t6;
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cost                       = xcosf(tht);
                        const register __m512 _2   = _mm512_set1_ps(2.0f);
                        cos2t                      = xcosf(_mm512_add_ps(tht,tht));
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        cos3t                      = xcosf(_mm512_add_ps(tht,_mm512_add_ps(tht,tht)));
                        const register __m512 c0   = _mm512_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c1   = _mm512_set1_ps(0.244444444444444444444444444444f);
                        k0a3                       = _mm512_mul_ps(k0a2,k0a);
                        const register __m512 c2   = _mm512_set1_ps(0.083333333333333333333333333333f);
                        k0a4                       = _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c3   = _mm512_set1_ps(1.0f/120.0f);
                        k0a6                       = _mm512_mul_ps(k0a3,k0a2);
                        const register __m512 c4   = _mm512_set1_ps(1807.0f/70.0f);
                        t0                         = _mm512_add_ps(half,cost);  // 0.5+cost
                        const register __m512 c5   = _mm512_set1_ps(2531.0f/105.0f);
                        t1                         = _mm512_sub_ps(c0,_mm512_fmadd_ps(c1,cost,
                                                                            _mm512_mul_ps(c2,cos2t)));
                        t1                         = _mm512_mul_ps(t1,k0a2);
                        const register __m512 c6   = _mm512_set1_ps(57.0f/42.0f);
                        t2                         = _mm512_sub_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm512_fmadd_ps(c6,cos2t,
                                                                                          _mm512_mul_ps(c0,cos3t))));
                        t3                         = _mm512_mul_ps(c3,_mm512_mul_ps(t2,k0a4));
                        const register __m512 c7   = _mm512_set1_ps(0.166666666666666666666666666667f);
                        const register __m512 c8   = _mm512_set1_ps(0.2f);
                        t4                         = _mm512_mul_ps(c7,_mm512_fmsub_ps(_4,cost,_1));
                        t5                         = _mm512_mul_ps(c8,_mm512_fmadd_ps(_2,cost,_1));
                        t5                         = _mm512_mul_ps(t5,k0a2);
                        t6                         = _mm512_mul_ps(k0a6,_mm512_add_ps(t4,t5)); // imaginary part
                        _mm512_storeu_ps(&S1i[0],t6);
                        
                        t4                         = _mm512_sub_ps(t0,_mm512_add_ps(t2,t3));
                        _mm512_storeu_ps(&S1r[0],_mm512_mul_ps(k0a3,t4);
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
                   void S2_f3214_zmm146r4( const __m512 k0a,
                                           const __m512 tht,
                                           __m512 * __restrict S2r,
                                           __m512 * __restrict S2i) {

                        register __m512 cost,cos2t,cos3t;
                        register __m512 k0a3,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,t4,t5,t6;
                        const register __m512 c9   = _mm512_set1_ps(0.125f);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cost                       = xcosf(tht);
                        const register __m512 _2   = _mm512_set1_ps(2.0f);
                        cos2t                      = xcosf(_mm512_add_ps(tht,tht));
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        cos3t                      = xcosf(_mm512_add_ps(tht,_mm512_add_ps(tht,tht)));
                        const register __m512 c0   = _mm512_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c1   = _mm512_set1_ps(23.0f/60.0f);
                        k0a3                       = _mm512_mul_ps(k0a2,k0a);
                        const register __m512 c2   = _mm512_set1_ps(0.055555555555555555555555555556f);
                        k0a4                       = _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c3   = _mm512_set1_ps(1.0f/60.0f);
                        k0a6                       = _mm512_mul_ps(k0a3,k0a2);
                        const register __m512 c4   = _mm512_set1_ps(-1343.0f/150.0f);
                        t0                         = _mm512_fmadd_ps(half,cost,_1);  // 0.5*cost+1
                        const register __m512 c5   = _mm512_set1_ps(3769.0f/280.0f);
                        t1                         = _mm512_sub_ps(c0,_mm512_fmsub_ps(c1,cost,
                                                                            _mm512_mul_ps(c2,cos2t)));
                        t1                         = _mm512_mul_ps(t1,k0a2);
                        const register __m512 c6   = _mm512_set1_ps(57.0f/63.0f);
                        t2                         = _mm512_add_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm512_fmadd_ps(c6,cos2t,
                                                                                          _mm512_mul_ps(c9,cos3t))));
                        t3                         = _mm512_mul_ps(c3,_mm512_mul_ps(t2,k0a4));
                        const register __m512 c7   = _mm512_set1_ps(0.166666666666666666666666666667f);
                        const register __m512 c8   = _mm512_set1_ps(0.2f);
                        t4                         = _mm512_mul_ps(c7,_mm512_sub_ps(_4,cost));
                        t5                         = _mm512_mul_ps(c8,_mm512_add_ps(_2,cost));
                        t5                         = _mm512_mul_ps(t5,k0a2);
                        t6                         = _mm512_mul_ps(k0a6,_mm512_add_ps(t4,t5)); // imaginary part
                        *S2i                       = t6;
                        t4                         = _mm512_add_ps(t0,_mm512_add_ps(t2,t3));
                        *S2r                       = _mm512_mul_ps(k0a3,t4);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3214_zmm146r4_a( const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                             const float * __restrict __ATTR_ALIGN__(64) ptht,
                                             float * __restrict __ATTR_ALIGN__(64) S2r,
                                             float * __restrict __ATTR_ALIGN__(64) S2i) {

                        register __m512 k0a = _mm512_load_ps(&pk0a[0]);
                        register __m512 tht = _mm512_load_ps(&ptht[0]);
                        register __m512 cost,cos2t,cos3t;
                        register __m512 k0a3,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,t4,t5,t6;
                        const register __m512 c9   = _mm512_set1_ps(0.125f);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cost                       = xcosf(tht);
                        const register __m512 _2   = _mm512_set1_ps(2.0f);
                        cos2t                      = xcosf(_mm512_add_ps(tht,tht));
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        cos3t                      = xcosf(_mm512_add_ps(tht,_mm512_add_ps(tht,tht)));
                        const register __m512 c0   = _mm512_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c1   = _mm512_set1_ps(23.0f/60.0f);
                        k0a3                       = _mm512_mul_ps(k0a2,k0a);
                        const register __m512 c2   = _mm512_set1_ps(0.055555555555555555555555555556f);
                        k0a4                       = _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c3   = _mm512_set1_ps(1.0f/60.0f);
                        k0a6                       = _mm512_mul_ps(k0a3,k0a2);
                        const register __m512 c4   = _mm512_set1_ps(-1343.0f/150.0f);
                        t0                         = _mm512_fmadd_ps(half,cost,_1);  // 0.5*cost+1
                        const register __m512 c5   = _mm512_set1_ps(3769.0f/280.0f);
                        t1                         = _mm512_sub_ps(c0,_mm512_fmsub_ps(c1,cost,
                                                                            _mm512_mul_ps(c2,cos2t)));
                        t1                         = _mm512_mul_ps(t1,k0a2);
                        const register __m512 c6   = _mm512_set1_ps(57.0f/63.0f);
                        t2                         = _mm512_add_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm512_fmadd_ps(c6,cos2t,
                                                                                          _mm512_mul_ps(c9,cos3t))));
                        t3                         = _mm512_mul_ps(c3,_mm512_mul_ps(t2,k0a4));
                        const register __m512 c7   = _mm512_set1_ps(0.166666666666666666666666666667f);
                        const register __m512 c8   = _mm512_set1_ps(0.2f);
                        t4                         = _mm512_mul_ps(c7,_mm512_sub_ps(_4,cost));
                        t5                         = _mm512_mul_ps(c8,_mm512_add_ps(_2,cost));
                        t5                         = _mm512_mul_ps(t5,k0a2);
                        t6                         = _mm512_mul_ps(k0a6,_mm512_add_ps(t4,t5)); // imaginary part
                        _mm512_store_ps(&S2i[0],   t6);
                        t4                         = _mm512_add_ps(t0,_mm512_add_ps(t2,t3));
                        _mm512_store_ps(&S2r[0],    _mm512_mul_ps(k0a3,t4);
                }



                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3214_zmm146r4_u( const float * __restrict  pk0a,
                                             const float * __restrict  ptht,
                                             float * __restrict  S2r,
                                             float * __restrict  S2i) {

                        register __m512 k0a = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 tht = _mm512_loadu_ps(&ptht[0]);
                        register __m512 cost,cos2t,cos3t;
                        register __m512 k0a3,k0a2,k0a4,k0a6;
                        register __m512 t0,t1,t2,t3,t4,t5,t6;
                        const register __m512 c9   = _mm512_set1_ps(0.125f);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        cost                       = xcosf(tht);
                        const register __m512 _2   = _mm512_set1_ps(2.0f);
                        cos2t                      = xcosf(_mm512_add_ps(tht,tht));
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        cos3t                      = xcosf(_mm512_add_ps(tht,_mm512_add_ps(tht,tht)));
                        const register __m512 c0   = _mm512_set1_ps(0.3333333333333333333333333333333f);
                        k0a2                       = _mm512_mul_ps(k0a,k0a);
                        const register __m512 c1   = _mm512_set1_ps(23.0f/60.0f);
                        k0a3                       = _mm512_mul_ps(k0a2,k0a);
                        const register __m512 c2   = _mm512_set1_ps(0.055555555555555555555555555556f);
                        k0a4                       = _mm512_mul_ps(k0a2,k0a2);
                        const register __m512 c3   = _mm512_set1_ps(1.0f/60.0f);
                        k0a6                       = _mm512_mul_ps(k0a3,k0a2);
                        const register __m512 c4   = _mm512_set1_ps(-1343.0f/150.0f);
                        t0                         = _mm512_fmadd_ps(half,cost,_1);  // 0.5*cost+1
                        const register __m512 c5   = _mm512_set1_ps(3769.0f/280.0f);
                        t1                         = _mm512_sub_ps(c0,_mm512_fmsub_ps(c1,cost,
                                                                            _mm512_mul_ps(c2,cos2t)));
                        t1                         = _mm512_mul_ps(t1,k0a2);
                        const register __m512 c6   = _mm512_set1_ps(57.0f/63.0f);
                        t2                         = _mm512_add_ps(c4,mm512_fmadd_ps(c5,cost,
                                                                        _mm512_fmadd_ps(c6,cos2t,
                                                                                          _mm512_mul_ps(c9,cos3t))));
                        t3                         = _mm512_mul_ps(c3,_mm512_mul_ps(t2,k0a4));
                        const register __m512 c7   = _mm512_set1_ps(0.166666666666666666666666666667f);
                        const register __m512 c8   = _mm512_set1_ps(0.2f);
                        t4                         = _mm512_mul_ps(c7,_mm512_sub_ps(_4,cost));
                        t5                         = _mm512_mul_ps(c8,_mm512_add_ps(_2,cost));
                        t5                         = _mm512_mul_ps(t5,k0a2);
                        t6                         = _mm512_mul_ps(k0a6,_mm512_add_ps(t4,t5)); // imaginary part
                        _mm512_storeu_ps(&S2i[0],   t6);
                        t4                         = _mm512_add_ps(t0,_mm512_add_ps(t2,t3));
                        _mm512_storeu_ps(&S2r[0],    _mm512_mul_ps(k0a3,t4);
                }


                  
                  /*
                       Formula 3.2-16, optics contribution at upper end of resonance region
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3216_zmm16r4(const __m512 k0a,
                                         const __m512 tht,
                                         __m512 * __restrict S1r,
                                         __m512 * __restrict S1i) {

                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        const register __m512 _7   = _mm512_set1_ps(7.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        register __m512 k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512 t0r,t0i,t1r,t1i,t2,t3,cos6t;
                        k0ah  = _mm512_mul_ps(half,k0a);
                        htht  = _mm512_mul_ps(half,tht);
                        cost  = xcosf(tht);
                        k0aa  = _mm512_add_ps(k0a,k0a);
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        sint  = xsinf(htht);
                        sin2t = _mm512_mul_ps(sint,sint);
                        carr  = _mm512_set1_ps(-0.0f);
                        cari  = _mm512_mul_ps(k0aa,htht);
                        cexp_zmm16r4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_ps(cost,cost);
                        cos3t = _mm512_mul_ps(cost,cos2t);
                        cos6t = _mm512_mul_ps(cos3t,cos2t);
                        t3    = _mm512_div_ps(sin2t,cos6t);
                        t0r   = _mm512_setzero_ps();
                        t0i   = _mm512_div_ps(Ii,_mm512_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_ps(_1,t0i); // second term
                        t2    = _mm512_div_ps(_7,_mm512_mul_ps(_4,k0a2));
                        cerr  = _mm512_mul_ps(t2,t3);
                        t1r   = _mm512_sub_ps(t0r,cerr);
                        t1i   = _mm512_sub_ps(t0i,cerr);
                        cmul_zmm16r4(cexpr,cexpi,t1r,t1i,&S1r,&S1i);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3216_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           float * __restrict __ATTR_ALIGN__(64) S1r,
                                           float * __restrict __ATTR_ALIGN__(64) S1i) {
     
                        const register __m512 k0a  = _mm512_load_ps(&pk0a[0]);
                        const register __m512 tht  = _mm512_load_ps(&ptht[0]);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        const register __m512 _7   = _mm512_set1_ps(7.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        register __m512 k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_ps(half,k0a);
                        htht  = _mm512_mul_ps(half,tht);
                        cost  = xcosf(tht);
                        k0aa  = _mm512_add_ps(k0a,k0a);
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        sint  = xsinf(htht);
                        sin2t = _mm512_mul_ps(sint,sint);
                        carr  = _mm512_set1_ps(-0.0f);
                        cari  = _mm512_mul_ps(k0aa,htht);
                        cexp_zmm16r4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_ps(cost,cost);
                        cos3t = _mm512_mul_ps(cost,cos2t);
                        cos6t = _mm512_mul_ps(cos3t,cos2t);
                        t3    = _mm512_div_ps(sin2t,cos6t);
                        t0r   = _mm512_setzero_ps();
                        t0i   = _mm512_div_ps(Ii,_mm512_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_ps(_1,t0i); // second term
                        t2    = _mm512_div_ps(_7,_mm512_mul_ps(_4,k0a2));
                        cerr  = _mm512_mul_ps(t2,t3);
                        t1r   = _mm512_sub_ps(t0r,cerr);
                        t1i   = _mm512_sub_ps(t0i,cerr);
                        cmul_zmm16r4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_store_ps(&S1r[0], resr);
                        _mm512_store_ps(&S1i[0], resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S1_f3216_zmm16r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  ptht,
                                           float * __restrict  S1r,
                                           float * __restrict  S1i) {
     
                        const register __m512 k0a  = _mm512_loadu_ps(&pk0a[0]);
                        const register __m512 tht  = _mm512_loadu_ps(&ptht[0]);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 _4   = _mm512_set1_ps(4.0f);
                        const register __m512 _7   = _mm512_set1_ps(7.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        register __m512 k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_ps(half,k0a);
                        htht  = _mm512_mul_ps(half,tht);
                        cost  = xcosf(tht);
                        k0aa  = _mm512_add_ps(k0a,k0a);
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        sint  = xsinf(htht);
                        sin2t = _mm512_mul_ps(sint,sint);
                        carr  = _mm512_set1_ps(-0.0f);
                        cari  = _mm512_mul_ps(k0aa,htht);
                        cexp_zmm16r4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_ps(cost,cost);
                        cos3t = _mm512_mul_ps(cost,cos2t);
                        cos6t = _mm512_mul_ps(cos3t,cos2t);
                        t3    = _mm512_div_ps(sin2t,cos6t);
                        t0r   = _mm512_setzero_ps();
                        t0i   = _mm512_div_ps(Ii,_mm512_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_ps(_1,t0i); // second term
                        t2    = _mm512_div_ps(_7,_mm512_mul_ps(_4,k0a2));
                        cerr  = _mm512_mul_ps(t2,t3);
                        t1r   = _mm512_sub_ps(t0r,cerr);
                        t1i   = _mm512_sub_ps(t0i,cerr);
                        cmul_zmm16r4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_storeu_ps(&S1r[0], resr);
                        _mm512_storeu_ps(&S1i[0], resi);
                }


                   /*
                       Formula 3.2-17, optics contribution at upper end of resonance region
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3217_zmm16r4(const __m512 k0a,
                                         const __m512 tht,
                                         __m512 * __restrict S2r,
                                         __m512 * __restrict S2i) {

                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 quat = _mm512_set1_ps(0.25f);
                        const register __m512 _6   = _mm512_set1_ps(6.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        register __m512 k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512 t0r,t0i,t1r,t1i,t2,t3,cos6t;
                        k0ah  = _mm512_mul_ps(half,k0a);
                        htht  = _mm512_mul_ps(half,tht);
                        cost  = xcosf(tht);
                        k0aa  = _mm512_add_ps(k0a,k0a);
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        sint  = xsinf(htht);
                        sin2t = _mm512_mul_ps(sint,sint);
                        carr  = _mm512_set1_ps(-0.0f);
                        cari  = _mm512_mul_ps(k0aa,htht);
                        cexp_zmm16r4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_ps(cost,cost);
                        cos3t = _mm512_mul_ps(cost,cos2t);
                        cos6t = _mm512_mul_ps(cos3t,cos2t);
                        t3    = _mm512_div_ps(_mm512_mul_ps(_mm512_add_ps(_6,cost),sin2t),cos6t);
                        t0r   = _mm512_setzero_ps();
                        t0i   = _mm512_div_ps(_mm512_mul_ps(Ii,cost),_mm512_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_ps(_1,t0i); // second term
                        t2    = _mm512_mul_ps(quat,k0a2);
                        cerr  = _mm512_mul_ps(t2,t3);
                        t1r   = _mm512_sub_ps(t0r,cerr);
                        t1i   = _mm512_sub_ps(t0i,cerr);
                        cmul_zmm16r4(cexpr,cexpi,t1r,t1i,&S2r,&S2i);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3217_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0a,
                                           const float * __restrict __ATTR_ALIGN__(64) ptht,
                                           float * __restrict __ATTR_ALIGN__(64) S2r,
                                           float * __restrict __ATTR_ALIGN__(64) S2i) {
                        
                        register __m512 k0a        = _mm512_load_ps(&pk0a[0]);
                        register __m512 tht        = _mm512_load_ps(&tht[0]);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 quat = _mm512_set1_ps(0.25f);
                        const register __m512 _6   = _mm512_set1_ps(6.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        register __m512 k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_ps(half,k0a);
                        htht  = _mm512_mul_ps(half,tht);
                        cost  = xcosf(tht);
                        k0aa  = _mm512_add_ps(k0a,k0a);
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        sint  = xsinf(htht);
                        sin2t = _mm512_mul_ps(sint,sint);
                        carr  = _mm512_set1_ps(-0.0f);
                        cari  = _mm512_mul_ps(k0aa,htht);
                        cexp_zmm16r4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_ps(cost,cost);
                        cos3t = _mm512_mul_ps(cost,cos2t);
                        cos6t = _mm512_mul_ps(cos3t,cos2t);
                        t3    = _mm512_div_ps(_mm512_mul_ps(_mm512_add_ps(_6,cost),sin2t),cos6t);
                        t0r   = _mm512_setzero_ps();
                        t0i   = _mm512_div_ps(_mm512_mul_ps(Ii,cost),_mm512_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_ps(_1,t0i); // second term
                        t2    = _mm512_mul_ps(quat,k0a2);
                        cerr  = _mm512_mul_ps(t2,t3);
                        t1r   = _mm512_sub_ps(t0r,cerr);
                        t1i   = _mm512_sub_ps(t0i,cerr);
                        cmul_zmm16r4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_store_ps(&S2r[0],resr);
                        _mm512_store_ps(&S2i[0],resi);
                }


                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void S2_f3217_zmm16r4_u(const float * __restrict  pk0a,
                                           const float * __restrict  ptht,
                                           float * __restrict  S2r,
                                           float * __restrict  S2i) {
                        
                        register __m512 k0a        = _mm512_loadu_ps(&pk0a[0]);
                        register __m512 tht        = _mm512_loadu_ps(&tht[0]);
                        const register __m512 _1   = _mm512_set1_ps(1.0f);
                        const register __m512 quat = _mm512_set1_ps(0.25f);
                        const register __m512 _6   = _mm512_set1_ps(6.0f);
                        const register __m512 half = _mm512_set1_ps(0.5f);
                        register __m512 k0ah,k0a2,k0aa,cexpr,cexpi;
                        register __m512 sint,sin2t,cost,cos2t,cos3t,carr,cari,htht; 
                        register __m512 t0r,t0i,t1r,t1i,t2,t3,cos6t,resr,resi;
                        k0ah  = _mm512_mul_ps(half,k0a);
                        htht  = _mm512_mul_ps(half,tht);
                        cost  = xcosf(tht);
                        k0aa  = _mm512_add_ps(k0a,k0a);
                        k0a2  = _mm512_mul_ps(k0a,k0a);
                        sint  = xsinf(htht);
                        sin2t = _mm512_mul_ps(sint,sint);
                        carr  = _mm512_set1_ps(-0.0f);
                        cari  = _mm512_mul_ps(k0aa,htht);
                        cexp_zmm16r4(carr,cari,&cexpr,&cexpi);
                        cexpr = _mm512_mul_ps(k0ah,cexpr);// exp term
                        cexpi = _mm512_mul_ps(k0ah,cexpi);// exp term
                        cos2t = _mm512_mul_ps(cost,cost);
                        cos3t = _mm512_mul_ps(cost,cos2t);
                        cos6t = _mm512_mul_ps(cos3t,cos2t);
                        t3    = _mm512_div_ps(_mm512_mul_ps(_mm512_add_ps(_6,cost),sin2t),cos6t);
                        t0r   = _mm512_setzero_ps();
                        t0i   = _mm512_div_ps(_mm512_mul_ps(Ii,cost),_mm512_mul_ps(k0aa,cos3t));
                        t0r   = _1;                    // second term
                        t0i   = _mm512_sub_ps(_1,t0i); // second term
                        t2    = _mm512_mul_ps(quat,k0a2);
                        cerr  = _mm512_mul_ps(t2,t3);
                        t1r   = _mm512_sub_ps(t0r,cerr);
                        t1i   = _mm512_sub_ps(t0i,cerr);
                        cmul_zmm16r4(cexpr,cexpi,t1r,t1i,&resr,&resi);
                        _mm512_storeu_ps(&S2r[0],resr);
                        _mm512_storeu_ps(&S2i[0],resi);
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
                   __m512 S1_f3220_zmm16r4(const __m512 k0a,
                                           const __m512 tht) {

                          register const __m512 half = _mm512_set1_ps(0.5f);
                          register __m512 S1;
                          register __m512 k0a3,cost;
                          k0a3 = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cost = xcosf(tht);
                          S1   = _mm512_mul_ps(k0a3,_mm512_add_ps(half,cost));
                          return (S1);
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
                   __m512 S2_f3221_zmm16r4(const __m512 k0a,
                                           const __m512 tht) {

                          register const __m512 half = _mm512_set1_ps(0.5f);
                          register const __m512 _1   = _mm512_set1_ps(1.0f);
                          register __m512 S2;
                          register __m512 k0a3,cost,t0; 
                          k0a3  = _mm512_mul_ps(k0a,_mm512_mul_ps(k0a,k0a));
                          cost  = xcosf(tht);
                          t0    = _mm512_fmadd_ps(half,cost,_1);
                          S2    = _mm512_mul_ps(k0a3,t0);
                          return (S2);
               }



     } // radiolocation

} // gms











#endif /*__GMS_RCS_SPHERE_ZMM16R4_HPP__*/
