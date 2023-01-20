
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
                   __m512 rcs_f419_zmm16r4(const float * __restrict __ATTR_ALIGN__(64) pa,
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





      } // radiolocation


} // gms









#endif /*__GMS_RCS_CYLINDER_ZMM16R4_HPP__*/
