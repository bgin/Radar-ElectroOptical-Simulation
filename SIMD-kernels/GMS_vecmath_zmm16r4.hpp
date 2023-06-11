

#ifndef __GMS_VECMATH_ZMM16R4_HPP__
#define __GMS_VECMATH_ZMM16R4_HPP__

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

    const unsigned int GMS_VECMATH_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_VECMATH_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_VECMATH_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_VECMATH_ZMM16R4_FULLVER =
      1000U*GMS_VECMATH_ZMM16R4_MAJOR+
      100U*GMS_VECMATH_ZMM16R4_MINOR+
      10U*GMS_VECMATH_ZMM16R4_MICRO;
    const char * const GMS_VECMATH_ZMM16R4_CREATION_DATE = "11-06-2023 09:36 AM +00200 (SUN 11 06 2023 GMT+2)";
    const char * const GMS_VECMATH_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_VECMATH_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_VECMATH_ZMM16R4_DESCRIPTION   = "Vector operations (dot,cross, ...etc) accelerated by avx512 single."

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_zmm16r4.hpp"

namespace gms {



          namespace math {
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 sdotv_zmm16r4(const __m512 v1x,
	                                const __m512 v1y,
	                                const __m512 v1z,
	                                const __m512 v2x,
	                                const __m512 v2y,
	                                const __m512 v2z) {
	                                
	                  register __m512 result;
	                  result = _mm512_fmadd_pd(v1x,v2x,
	                                      _mm512_fmadd_pd(v1y,v2y,
	                                                 _mm512_mul_pd(v1z,v2z)));
	                  return (result);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 sdotv_zmm16r4(const float * __restrict __ATTR_ALIGN__(64) pv1x,
	                                const float * __restrict __ATTR_ALIGN__(64) pv1y,
	                                const float * __restrict __ATTR_ALIGN__(64) pv1z,
	                                const float * __restrict __ATTR_ALIGN__(64) pv2x,
	                                const float * __restrict __ATTR_ALIGN__(64) pv2y,
	                                const float * __restrict __ATTR_ALIGN__(64) pv2z) {
	                          
	                  register __m512 v1x = _mm512_load_ps(&pv1x[0]);
	                  register __m512 v1y = _mm512_load_ps(&pv1y[0]);  
	                  register __m512 v1z = _mm512_load_ps(&pv1z[0]); 
	                  register __m512 v2x = _mm512_load_ps(&pv2x[0]);  
	                  register __m512 v2y = _mm512_load_ps(&pv2y[0]); 
	                  register __m512 v2z = _mm512_load_ps(&pv2z[0]);
	                  register __m512 result;
	                  result = _mm512_fmadd_pd(v1x,v2x,
	                                      _mm512_fmadd_pd(v1y,v2y,
	                                                 _mm512_mul_pd(v1z,v2z)));
	                  return (result);                       
	        }
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cdotv_zmm16r4(
                
                
        } // math

} // gms


#endif /*__GMS_VECMATH_ZMM16R4_HPP__*/
