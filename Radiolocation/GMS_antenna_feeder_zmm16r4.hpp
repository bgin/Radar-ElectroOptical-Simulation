

#ifndef __GMS_ANTENNA_FEEDER_ZMM16R4_HPP__
#define __GMS_ANTENNA_FEEDER_ZMM16R4_HPP__ 090620230852


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

    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_FULLVER =
      1000U*GMS_ANTENNA_FEEDER_ZMM16R4_MAJOR+
      100U*GMS_ANTENNA_FEEDER_ZMM16R4_MINOR+
      10U*GMS_ANTENNA_FEEDER_ZMM16R4_MICRO;
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_CREATION_DATE = "09-06-2023 08:52 PM +00200 (FRI 09 JUN 2023 GMT+2)";
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_DESCRIPTION   = "AVX512 (single) optimized antenna and feeder kernels.";

}


#include <immintrin.h>
#include <complex>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"
#include "GMS_cspint_quad.hpp"
#include "GMS_avint_quad.hpp"
#include "GMS_cephes.h"


namespace gms {

 
          namespace radiolocation {
          
           
                   /*
                       Spherical unit vectors.
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void spher_unitv_zmm16r4(const __m512 ex,
	                                   const __m512 ey,
	                                   const __m512 ez,
	                                   const __m512 tht,
	                                   const __m512 phi,
	                                   __m512 * __restrict er,
	                                   __m512 * __restrict eth,
	                                   __m512 * __restrict eph) {
	               
	               using namespace gms::math;                   
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 t0,t1;
	               stht = xsinf(tht);
	               sphi = xsinf(phi);
	               t0   = _mm512_mul_ps(sphi,ey);
	               cphi = xcosf(phi);
	               t1   = _mm512_mul_ps(cphi,ex);
	               ctht = xcosf(tht);
	               *er  = _mm512_fmadd_ps(stht,t1,
	                                 _mm512_fmadd_ps(stht,t0,
	                                             _mm512_mul_ps(ctht,ez)));    
	               *eth = _mm512_fmadd_ps(ctht,t1,
	                                 _mm512_fmsub_ps(ctht,t0,
	                                             _mm512_mul_ps(stht,ez)));  
	               *eph = _mm512_fmadd_ps(negate_zmm16r4(sphi),ex,
	                                             _mm512_mul_ps(cphi,ey));                           
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void spher_unitv_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pex,
	                                      const float * __restrict __ATTR_ALIGN__(64) pey,
	                                      const float * __restrict __ATTR_ALIGN__(64) pez,
	                                      const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                      const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                      float * __restrict __ATTR_ALIGN__(64)  er,
	                                      float * __restrict __ATTR_ALIGN__(64)  eth,
	                                      float * __restrict __ATTR_ALIGN__(64)  eph) {
	               
	               using namespace gms::math;     
	               register __m512 ex = _mm512_load_ps(&pex[0]);
	               register __m512 ey = _mm512_load_ps(&pey[0]);
	               register __m512 ez = _mm512_load_ps(&pez[0]);              
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 t0,t1;
	               stht = xsinf(tht);
	               sphi = xsinf(phi);
	               t0   = _mm512_mul_ps(sphi,ey);
	               cphi = xcosf(phi);
	               t1   = _mm512_mul_ps(cphi,ex);
	               ctht = xcosf(tht);
	               _mm512_store_ps(&er[0]  ,_mm512_fmadd_ps(stht,t1,
	                                                  _mm512_fmadd_ps(stht,t0,
	                                                       _mm512_mul_ps(ctht,ez))));    
	               _mm512_store_ps(&eth[0] ,_mm512_fmadd_ps(ctht,t1,
	                                                  _mm512_fmsub_ps(ctht,t0,
	                                                       _mm512_mul_ps(stht,ez))));  
	               _mm512_store_ps(&eph[0] ,_mm512_fmadd_ps(negate_zmm16r4(sphi),ex,
	                                                       _mm512_mul_ps(cphi,ey)));                           
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void spher_unitv_zmm16r4_u(const float * __restrict  pex,
	                                      const float * __restrict  pey,
	                                      const float * __restrict  pez,
	                                      const float * __restrict  ptht,
	                                      const float * __restrict  pphi,
	                                      float * __restrict  er,
	                                      float * __restrict  eth,
	                                      float * __restrict  eph) {
	               
	               using namespace gms::math;     
	               register __m512 ex = _mm512_load_ps(&pex[0]);
	               register __m512 ey = _mm512_load_ps(&pey[0]);
	               register __m512 ez = _mm512_load_ps(&pez[0]);              
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 t0,t1;
	               stht = xsinf(tht);
	               sphi = xsinf(phi);
	               t0   = _mm512_mul_ps(sphi,ey);
	               cphi = xcosf(phi);
	               t1   = _mm512_mul_ps(cphi,ex);
	               ctht = xcosf(tht);
	               _mm512_store_ps(&er[0]  ,_mm512_fmadd_ps(stht,t1,
	                                                  _mm512_fmadd_ps(stht,t0,
	                                                       _mm512_mul_ps(ctht,ez))));    
	               _mm512_store_ps(&eth[0] ,_mm512_fmadd_ps(ctht,t1,
	                                                  _mm512_fmsub_ps(ctht,t0,
	                                                       _mm512_mul_ps(stht,ez))));  
	               _mm512_store_ps(&eph[0] ,_mm512_fmadd_ps(negate_zmm16r4(sphi),ex,
	                                                       _mm512_mul_ps(cphi,ey)));                           
	       }
	       
	       /*
	             Function 'N' = Nth(theta,phi)*eth+Nphi(theta,phi)*ephi
	             Formula 1-3, p. 10
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Rad_pattern_f13_zmm16r4(const __m512 eth,
	                                        const __m512 eph,
	                                        const __m512 nthr,
	                                        const __m512 nthi,
	                                        const __m512 nphr,
	                                        const __m512 nphi,
	                                        __m512 * __restrict Nthr,
	                                        __m512 * __restrict Nthi,
	                                        __m512 * __restrict Nphr,
	                                        __m512 * __restrict Nphi) {
	                                      
	                *Nthr = _mm512_mul_ps(nthr,eth);
	                *Nphr = _mm512_mul_ps(nphr,eph);
	                *Nthi = _mm512_mul_ps(nthi,eth);
	                *Nphi = _mm512_mul_ps(nphi,eph);                             
	       }
	       
	       
	       
             
        } // radiolocation

} // gms












#endif /*__GMS_ANTENNA_FEEDER_ZMM16R4_HPP__*/
