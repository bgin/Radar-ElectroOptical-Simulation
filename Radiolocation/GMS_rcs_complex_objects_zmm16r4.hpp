

#ifndef __GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__
#define __GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__

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

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_CREATION_DATE = "11-05-2023 10:53 PM +00200 (THR 11 MAY 2023 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"




namespace  gms {


         namespace radiolocation {
         
              
              /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Formula: 8.1-21
              */     
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8121_zmm16r4(const __m512 gam,
                                             const __m512 phi,
                                             const __m512 k0,
                                             __m512 * __restrict D1r,
                                             __m512 * __restrict D1i,
                                            __m512 * __restrict D2r,
                                            __m512 * __restrict D2i) {
                                            
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);
                        register __m512 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm512_setzero_ps();
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = xsinf(x1);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cos1 = xcosf(x1);
                        x0   = _mm512_mul_ps(invn,phi2);
                        cos2 = xcosf(x0);
                        cer  = _mm512_mul_ps(cer,sinp);
                        trm1 = _mm512_rcp14_ps(_mm512_sub_ps(cos1),C10);
                        cei  = _mm512_mul_ps(cei,sinp);
                        trm2 = _mm512_rcp14_ps(_mm512_sub_ps(cos1,cos2));
                        sqr  = _mm512_mul_ps(invn,sqr);
                        ear  = _mm512_mul_ps(cer,sqr);
                        eai  = _mm512_mul_ps(cei,sqr);
                        x0   = _mm512_sub_ps(trm1,trm2);
                        *D1r = _mm512_mul_ps(ear,x0);
                        *D1i = _mm512_mul_ps(eai,x0);
                        x1   = _mm512_add_ps(trm1,trm2);
                        *D2r = _mm512_mul_ps(ear,x1);
                        *D2i = _mm512_mul_ps(eai,x1);
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8121_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam,
                                               const float * __restrict __ATTR_ALIGN__(64) pphi,
                                               const float * __restrict __ATTR_ALIGN__(64) k0,
                                               float * __restrict __ATTR_ALIGN__(64) D1r,
                                               float * __restrict __ATTR_ALIGN__(64) D1i,
                                               float * __restrict __ATTR_ALIGN__(64) D2r,
                                               float * __restrict __ATTR_ALIGN__(64) D2i) {
                                   
                        register __m512 gam = _mm512_load_ps(&pgam[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);  
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);       
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);
                        register __m512 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm512_setzero_ps();
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = xsinf(x1);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cos1 = xcosf(x1);
                        x0   = _mm512_mul_ps(invn,phi2);
                        cos2 = xcosf(x0);
                        cer  = _mm512_mul_ps(cer,sinp);
                        trm1 = _mm512_rcp14_ps(_mm512_sub_ps(cos1),C10);
                        cei  = _mm512_mul_ps(cei,sinp);
                        trm2 = _mm512_rcp14_ps(_mm512_sub_ps(cos1,cos2));
                        sqr  = _mm512_mul_ps(invn,sqr);
                        ear  = _mm512_mul_ps(cer,sqr);
                        eai  = _mm512_mul_ps(cei,sqr);
                        x0   = _mm512_sub_ps(trm1,trm2);
                        _mm512_store_ps(&D1r[0] ,_mm512_mul_ps(ear,x0));
                        _mm512_store_ps(&D1i[0] ,_mm512_mul_ps(eai,x0));
                        x1   = _mm512_add_ps(trm1,trm2);
                        _mm512_store_ps(&D2r[0] ,_mm512_mul_ps(ear,x1));
                        _mm512_store_ps(&D2i[0] ,_mm512_mul_ps(eai,x1));
                }
                
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__*/
