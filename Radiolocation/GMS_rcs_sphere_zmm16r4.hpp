
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



     } // radiolocation

} // gms











#endif /*__GMS_RCS_SPHERE_ZMM16R4_HPP__*/
