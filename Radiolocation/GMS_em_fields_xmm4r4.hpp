

#ifndef __GMS_EM_FIELDS_XMM4R4_HPP__
#define __GMS_EM_FIELDS_XMM4R4_HPP__ 301020231146

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

    const unsigned int GMS_EM_FIELDS_XMM4R4_MAJOR = 1U;
    const unsigned int GMS_EM_FIELDS_XMM4R4_MINOR = 0U;
    const unsigned int GMS_EM_FIELDS_XMM4R4_MICRO = 0U;
    const unsigned int GMS_EM_FIELDS_XMM4R4_FULLVER =
      1000U*GMS_EM_FIELDS_XMM4R4_MAJOR+
      100U*GMS_EM_FIELDS_XMM4R4_MINOR+
      10U*GMS_EM_FIELDS_XMM4R4_MICRO;
    const char * const GMS_EM_FIELDS_XMM4R4_CREATION_DATE = "30-10-2023 11:46 AM +00200 (MON 30 10 2023 GMT+2)";
    const char * const GMS_EM_FIELDS_XMM4R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EM_FIELDS_XMM4R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EM_FIELDS_XMM4R4_DESCRIPTION   = " Computational ElectroMagnetics related helper routines."
                       

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_complex_xmm4r4.hpp"

#ifndef __EM_FIELDS_PF_CACHE_HINT__
#define __EM_FIELDS_PF_CACHE_HINT__ 1
#endif 

namespace gms {



          namespace radiolocation {
          
          
 
              
               
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 sdotv_xmm4c4(const __m128 v1x,
	                                const __m128 v1y,
	                                const __m128 v1z,
	                                const __m128 v2x,
	                                const __m128 v2y,
	                                const __m128 v2z) {
	                                
	                  register __m128 result;
	                  result = _mm_fmadd_ps(v1x,v2x,
	                                      _mm_fmadd_ps(v1y,v2y,
	                                                 _mm_mul_ps(v1z,v2z)));
	                  return (result);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void sdotv_xmm4c4_unroll16x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m128 v1x;
	                register __m128 v1y;
	                register __m128 v1z;
	                register __m128 v2x;
	                register __m128 v2y;
	                register __m128 v2z;
	                register __m128 dtv;
	                int32_t j,m,m1;
	                
	                m = n%16;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                   
	                        v1x = pv1x[j];
	                        v2x = pv2x[j];
	                        v1y = pv1y[j];
	                        v2y = pv2y[j];
	                        v1z = pv1z[j];
	                        v2z = pv2z[j];
	                        dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                            v2x,v2y,v2z);
	                        pdtv[j] = dtv;
	                   }
	                   if(n<16) { return;}
	                }  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 16) {
	                     v1x = pv1x[j+0];
	                     v2x = pv2x[j+0];
	                     v1y = pv1y[j+0];
	                     v2y = pv2y[j+0];
	                     v1z = pv1z[j+0];
	                     v2z = pv2z[j+0];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+0] = dtv;
	                     v1x = pv1x[j+1];
	                     v2x = pv2x[j+1];
	                     v1y = pv1y[j+1];
	                     v2y = pv2y[j+1];
	                     v1z = pv1z[j+1];
	                     v2z = pv2z[j+1];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+1] = dtv;
	                     v1x = pv1x[j+2];
	                     v2x = pv2x[j+2];
	                     v1y = pv1y[j+2];
	                     v2y = pv2y[j+2];
	                     v1z = pv1z[j+2];
	                     v2z = pv2z[j+2];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+2] = dtv;
	                     v1x = pv1x[j+3];
	                     v2x = pv2x[j+3];
	                     v1y = pv1y[j+3];
	                     v2y = pv2y[j+3];
	                     v1z = pv1z[j+3];
	                     v2z = pv2z[j+3];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+3] = dtv;
	                     v1x = pv1x[j+4];
	                     v2x = pv2x[j+4];
	                     v1y = pv1y[j+4];
	                     v2y = pv2y[j+4];
	                     v1z = pv1z[j+4];
	                     v2z = pv2z[j+4];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+4] = dtv;
	                     v1x = pv1x[j+5];
	                     v2x = pv2x[j+5];
	                     v1y = pv1y[j+5];
	                     v2y = pv2y[j+5];
	                     v1z = pv1z[j+5];
	                     v2z = pv2z[j+5];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+5] = dtv;
	                     v1x = pv1x[j+6];
	                     v2x = pv2x[j+6];
	                     v1y = pv1y[j+6];
	                     v2y = pv2y[j+6];
	                     v1z = pv1z[j+6];
	                     v2z = pv2z[j+6];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+6] = dtv;
	                     v1x = pv1x[j+7];
	                     v2x = pv2x[j+7];
	                     v1y = pv1y[j+7];
	                     v2y = pv2y[j+7];
	                     v1z = pv1z[j+7];
	                     v2z = pv2z[j+7];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+7] = dtv;
	                     v1x = pv1x[j+8];
	                     v2x = pv2x[j+8];
	                     v1y = pv1y[j+8];
	                     v2y = pv2y[j+8];
	                     v1z = pv1z[j+8];
	                     v2z = pv2z[j+8];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+8] = dtv;
	                     v1x = pv1x[j+9];
	                     v2x = pv2x[j+9];
	                     v1y = pv1y[j+9];
	                     v2y = pv2y[j+9];
	                     v1z = pv1z[j+9];
	                     v2z = pv2z[j+9];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+9] = dtv;
	                     v1x = pv1x[j+10];
	                     v2x = pv2x[j+10];
	                     v1y = pv1y[j+10];
	                     v2y = pv2y[j+10];
	                     v1z = pv1z[j+10];
	                     v2z = pv2z[j+10];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+10] = dtv;
	                     v1x = pv1x[j+11];
	                     v2x = pv2x[j+11];
	                     v1y = pv1y[j+11];
	                     v2y = pv2y[j+11];
	                     v1z = pv1z[j+11];
	                     v2z = pv2z[j+11];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+11] = dtv;
	                     v1x = pv1x[j+12];
	                     v2x = pv2x[j+12];
	                     v1y = pv1y[j+12];
	                     v2y = pv2y[j+12];
	                     v1z = pv1z[j+12];
	                     v2z = pv2z[j+12];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+12] = dtv;
	                     v1x = pv1x[j+13];
	                     v2x = pv2x[j+13];
	                     v1y = pv1y[j+13];
	                     v2y = pv2y[j+13];
	                     v1z = pv1z[j+13];
	                     v2z = pv2z[j+13];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+13] = dtv;
	                     v1x = pv1x[j+14];
	                     v2x = pv2x[j+14];
	                     v1y = pv1y[j+14];
	                     v2y = pv2y[j+14];
	                     v1z = pv1z[j+14];
	                     v2z = pv2z[j+14];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+14] = dtv;
	                     v1x = pv1x[j+15];
	                     v2x = pv2x[j+15];
	                     v1y = pv1y[j+15];
	                     v2y = pv2y[j+15];
	                     v1z = pv1z[j+15];
	                     v2z = pv2z[j+15];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+15] = dtv;
	                }
	                                        
	      }
	      
	      
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void sdotv_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m128 v1x;
	                register __m128 v1y;
	                register __m128 v1z;
	                register __m128 v2x;
	                register __m128 v2y;
	                register __m128 v2z;
	                register __m128 dtv;
	                int32_t j,m,m1;
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                   
	                        v1x = pv1x[j];
	                        v2x = pv2x[j];
	                        v1y = pv1y[j];
	                        v2y = pv2y[j];
	                        v1z = pv1z[j];
	                        v2z = pv2z[j];
	                        dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                            v2x,v2y,v2z);
	                        pdtv[j] = dtv;
	                   }
	                   if(n<10) { return;}
	                }  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 10) {
	                     v1x = pv1x[j+0];
	                     v2x = pv2x[j+0];
	                     v1y = pv1y[j+0];
	                     v2y = pv2y[j+0];
	                     v1z = pv1z[j+0];
	                     v2z = pv2z[j+0];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+0] = dtv;
	                     v1x = pv1x[j+1];
	                     v2x = pv2x[j+1];
	                     v1y = pv1y[j+1];
	                     v2y = pv2y[j+1];
	                     v1z = pv1z[j+1];
	                     v2z = pv2z[j+1];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+1] = dtv;
	                     v1x = pv1x[j+2];
	                     v2x = pv2x[j+2];
	                     v1y = pv1y[j+2];
	                     v2y = pv2y[j+2];
	                     v1z = pv1z[j+2];
	                     v2z = pv2z[j+2];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+2] = dtv;
	                     v1x = pv1x[j+3];
	                     v2x = pv2x[j+3];
	                     v1y = pv1y[j+3];
	                     v2y = pv2y[j+3];
	                     v1z = pv1z[j+3];
	                     v2z = pv2z[j+3];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+3] = dtv;
	                     v1x = pv1x[j+4];
	                     v2x = pv2x[j+4];
	                     v1y = pv1y[j+4];
	                     v2y = pv2y[j+4];
	                     v1z = pv1z[j+4];
	                     v2z = pv2z[j+4];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+4] = dtv;
	                     v1x = pv1x[j+5];
	                     v2x = pv2x[j+5];
	                     v1y = pv1y[j+5];
	                     v2y = pv2y[j+5];
	                     v1z = pv1z[j+5];
	                     v2z = pv2z[j+5];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+5] = dtv;
	                     v1x = pv1x[j+6];
	                     v2x = pv2x[j+6];
	                     v1y = pv1y[j+6];
	                     v2y = pv2y[j+6];
	                     v1z = pv1z[j+6];
	                     v2z = pv2z[j+6];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+6] = dtv;
	                     v1x = pv1x[j+7];
	                     v2x = pv2x[j+7];
	                     v1y = pv1y[j+7];
	                     v2y = pv2y[j+7];
	                     v1z = pv1z[j+7];
	                     v2z = pv2z[j+7];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+7] = dtv;
	                     v1x = pv1x[j+8];
	                     v2x = pv2x[j+8];
	                     v1y = pv1y[j+8];
	                     v2y = pv2y[j+8];
	                     v1z = pv1z[j+8];
	                     v2z = pv2z[j+8];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+8] = dtv;
	                     v1x = pv1x[j+9];
	                     v2x = pv2x[j+9];
	                     v1y = pv1y[j+9];
	                     v2y = pv2y[j+9];
	                     v1z = pv1z[j+9];
	                     v2z = pv2z[j+9];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+9] = dtv;
	                    
	                }
	                                        
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void sdotv_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m128 v1x;
	                register __m128 v1y;
	                register __m128 v1z;
	                register __m128 v2x;
	                register __m128 v2y;
	                register __m128 v2z;
	                register __m128 dtv;
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                   
	                        v1x = pv1x[j];
	                        v2x = pv2x[j];
	                        v1y = pv1y[j];
	                        v2y = pv2y[j];
	                        v1z = pv1z[j];
	                        v2z = pv2z[j];
	                        dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                            v2x,v2y,v2z);
	                        pdtv[j] = dtv;
	                   }
	                   if(n<6) { return;}
	                }  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 6) {
	                     v1x = pv1x[j+0];
	                     v2x = pv2x[j+0];
	                     v1y = pv1y[j+0];
	                     v2y = pv2y[j+0];
	                     v1z = pv1z[j+0];
	                     v2z = pv2z[j+0];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+0] = dtv;
	                     v1x = pv1x[j+1];
	                     v2x = pv2x[j+1];
	                     v1y = pv1y[j+1];
	                     v2y = pv2y[j+1];
	                     v1z = pv1z[j+1];
	                     v2z = pv2z[j+1];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+1] = dtv;
	                     v1x = pv1x[j+2];
	                     v2x = pv2x[j+2];
	                     v1y = pv1y[j+2];
	                     v2y = pv2y[j+2];
	                     v1z = pv1z[j+2];
	                     v2z = pv2z[j+2];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+2] = dtv;
	                     v1x = pv1x[j+3];
	                     v2x = pv2x[j+3];
	                     v1y = pv1y[j+3];
	                     v2y = pv2y[j+3];
	                     v1z = pv1z[j+3];
	                     v2z = pv2z[j+3];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+3] = dtv;
	                     v1x = pv1x[j+4];
	                     v2x = pv2x[j+4];
	                     v1y = pv1y[j+4];
	                     v2y = pv2y[j+4];
	                     v1z = pv1z[j+4];
	                     v2z = pv2z[j+4];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+4] = dtv;
	                     v1x = pv1x[j+5];
	                     v2x = pv2x[j+5];
	                     v1y = pv1y[j+5];
	                     v2y = pv2y[j+5];
	                     v1z = pv1z[j+5];
	                     v2z = pv2z[j+5];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+5] = dtv;
	                   	                    
	                }
	                                        
	      }
	      
	      
	         
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void sdotv_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m128 v1x;
	                register __m128 v1y;
	                register __m128 v1z;
	                register __m128 v2x;
	                register __m128 v2y;
	                register __m128 v2z;
	                register __m128 dtv;
	                int32_t j,m,m1;
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                   
	                        v1x = pv1x[j];
	                        v2x = pv2x[j];
	                        v1y = pv1y[j];
	                        v2y = pv2y[j];
	                        v1z = pv1z[j];
	                        v2z = pv2z[j];
	                        dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                            v2x,v2y,v2z);
	                        pdtv[j] = dtv;
	                   }
	                   if(n<2) { return;}
	                }  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 2) {
	                     v1x = pv1x[j+0];
	                     v2x = pv2x[j+0];
	                     v1y = pv1y[j+0];
	                     v2y = pv2y[j+0];
	                     v1z = pv1z[j+0];
	                     v2z = pv2z[j+0];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+0] = dtv;
	                     v1x = pv1x[j+1];
	                     v2x = pv2x[j+1];
	                     v1y = pv1y[j+1];
	                     v2y = pv2y[j+1];
	                     v1z = pv1z[j+1];
	                     v2z = pv2z[j+1];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+1] = dtv;
	                                      	                    
	                }
	                                        
	      }
	      
	      
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void sdotv_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m128 v1x;
	                register __m128 v1y;
	                register __m128 v1z;
	                register __m128 v2x;
	                register __m128 v2y;
	                register __m128 v2z;
	                register __m128 dtv;
	                int32_t j;
	              	for(j = 0; j != n; ++j) {
	                     v1x = pv1x[j];
	                     v2x = pv2x[j];
	                     v1y = pv1y[j];
	                     v2y = pv2y[j];
	                     v1z = pv1z[j];
	                     v2z = pv2z[j];
	                     dtv = sdotv_xmm4c4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j] = dtv;
	                               
	                }
	                                        
	      }
	      
	      
	      
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 sdotv_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) pv1x,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv1y,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv1z,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv2x,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv2y,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv2z) {
	                          
	                  register __m128 v1x = _mm_load_ps(&pv1x[0]);
	                  register __m128 v1y = _mm_load_ps(&pv1y[0]);  
	                  register __m128 v1z = _mm_load_ps(&pv1z[0]); 
	                  register __m128 v2x = _mm_load_ps(&pv2x[0]);  
	                  register __m128 v2y = _mm_load_ps(&pv2y[0]); 
	                  register __m128 v2z = _mm_load_ps(&pv2z[0]);
	                  register __m128 result;
	                  result = _mm_fmadd_ps(v1x,v2x,
	                                      _mm_fmadd_ps(v1y,v2y,
	                                                 _mm_mul_ps(v1z,v2z)));
	                  return (result);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 sdotv_xmm4c4_u(const float * __restrict  pv1x,
	                                const float * __restrict  pv1y,
	                                const float * __restrict  pv1z,
	                                const float * __restrict  pv2x,
	                                const float * __restrict  pv2y,
	                                const float * __restrict  pv2z) {
	                          
	                  register __m128 v1x = _mm_loadu_ps(&pv1x[0]);
	                  register __m128 v1y = _mm_loadu_ps(&pv1y[0]);  
	                  register __m128 v1z = _mm_loadu_ps(&pv1z[0]); 
	                  register __m128 v2x = _mm_loadu_ps(&pv2x[0]);  
	                  register __m128 v2y = _mm_loadu_ps(&pv2y[0]); 
	                  register __m128 v2z = _mm_loadu_ps(&pv2z[0]);
	                  register __m128 result;
	                  result = _mm_fmadd_ps(v1x,v2x,
	                                      _mm_fmadd_ps(v1y,v2y,
	                                                 _mm_mul_ps(v1z,v2z)));
	                  return (result);                       
	        }
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void cdotv_xmm4c4( const xmm4c4_t v1x,
	                              const xmm4c4_t v1y,
	                              const xmm4c4_t v1z,
	                              const xmm4c4_t v2x,
	                              const xmm4c4_t v2y,
	                              const xmm4c4_t v2z,
	                              xmm4c4_t & res) {
	                              
	                xmm4c4_t tx,ty,tz;
	                cmul_xmm4c4(v1x.re,v1x.im,v2x.re,
	                                  v2x.im,&tx.re,&tx.im); 
	                cmul_xmm4c4(v1y.re,v1y.im,v2y.re,
	                                  v2y.im,&ty.re,&ty.im);
	                cmul_xmm4c4(v1z.re,v1z.im,v2z.re,
	                                  v2z.im,&tz.re,&tz.im);
	                res.re = _mm_add_ps(tx.re,
	                                   _mm_add_ps(ty.re,tz.re));
	                res.im = _mm_add_ps(tx.im,
	                                   _mm_add_ps(ty.im,tz.im));                   
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cdotv_xmm4c4_unroll16x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 16;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t v2z;
	                __ATTR_ALIGN__(16) xmm4c4_t res;
	                int32_t j,m,m1;
	                
	                m = n%16;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       cdotv_xmm4c4(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<16) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 16) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1                
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif

	                    v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+0] = res;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+1] = res;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+2] = res;
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+3] = res;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+4] = res;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+5] = res;
	                    v1x = pv1x[j+6];
	                    v2x = pv2x[j+6];
	                    v1y = pv1y[j+6];
	                    v2y = pv2y[j+6];
	                    v1z = pv1z[j+6];
	                    v2z = pv2z[j+6];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+6] = res;
	                    v1x = pv1x[j+7];
	                    v2x = pv2x[j+7];
	                    v1y = pv1y[j+7];
	                    v2y = pv2y[j+7];
	                    v1z = pv1z[j+7];
	                    v2z = pv2z[j+7];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+7] = res;
	                    v1x = pv1x[j+8];
	                    v2x = pv2x[j+8];
	                    v1y = pv1y[j+8];
	                    v2y = pv2y[j+8];
	                    v1z = pv1z[j+8];
	                    v2z = pv2z[j+8];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+8] = res;
	                    v1x = pv1x[j+9];
	                    v2x = pv2x[j+9];
	                    v1y = pv1y[j+9];
	                    v2y = pv2y[j+9];
	                    v1z = pv1z[j+9];
	                    v2z = pv2z[j+9];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+9] = res;
	                    v1x = pv1x[j+10];
	                    v2x = pv2x[j+10];
	                    v1y = pv1y[j+10];
	                    v2y = pv2y[j+10];
	                    v1z = pv1z[j+10];
	                    v2z = pv2z[j+10];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+10] = res;
	                    v1x = pv1x[j+11];
	                    v2x = pv2x[j+11];
	                    v1y = pv1y[j+11];
	                    v2y = pv2y[j+11];
	                    v1z = pv1z[j+11];
	                    v2z = pv2z[j+11];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+11] = res;
	                    v1x = pv1x[j+12];
	                    v2x = pv2x[j+12];
	                    v1y = pv1y[j+12];
	                    v2y = pv2y[j+12];
	                    v1z = pv1z[j+12];
	                    v2z = pv2z[j+12];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+12] = res;
	                    v1x = pv1x[j+13];
	                    v2x = pv2x[j+13];
	                    v1y = pv1y[j+13];
	                    v2y = pv2y[j+13];
	                    v1z = pv1z[j+13];
	                    v2z = pv2z[j+13];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+13] = res; 
	                    v1x = pv1x[j+14];
	                    v2x = pv2x[j+14];
	                    v1y = pv1y[j+14];
	                    v2y = pv2y[j+14];
	                    v1z = pv1z[j+14];
	                    v2z = pv2z[j+14];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+14] = res;
	                    v1x = pv1x[j+15];
	                    v2x = pv2x[j+15];
	                    v1y = pv1y[j+15];
	                    v2y = pv2y[j+15];
	                    v1z = pv1z[j+15];
	                    v2z = pv2z[j+15];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+15] = res;
	                }          
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cdotv_xmm4c4_unroll10x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t v2z;
	                __ATTR_ALIGN__(16) xmm4c4_t res;
	                int32_t j,m,m1;
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       cdotv_xmm4c4(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<10) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1                
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                   
	                    v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+0] = res;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+1] = res;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+2] = res;
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+3] = res;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+4] = res;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+5] = res;
	                    v1x = pv1x[j+6];
	                    v2x = pv2x[j+6];
	                    v1y = pv1y[j+6];
	                    v2y = pv2y[j+6];
	                    v1z = pv1z[j+6];
	                    v2z = pv2z[j+6];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+6] = res;
	                    v1x = pv1x[j+7];
	                    v2x = pv2x[j+7];
	                    v1y = pv1y[j+7];
	                    v2y = pv2y[j+7];
	                    v1z = pv1z[j+7];
	                    v2z = pv2z[j+7];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+7] = res;
	                    v1x = pv1x[j+8];
	                    v2x = pv2x[j+8];
	                    v1y = pv1y[j+8];
	                    v2y = pv2y[j+8];
	                    v1z = pv1z[j+8];
	                    v2z = pv2z[j+8];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+8] = res;
	                    v1x = pv1x[j+9];
	                    v2x = pv2x[j+9];
	                    v1y = pv1y[j+9];
	                    v2y = pv2y[j+9];
	                    v1z = pv1z[j+9];
	                    v2z = pv2z[j+9];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+9] = res;
	                
	             }          
	       }
	       
	       
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cdotv_xmm4c4_unroll6x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t v2z;
	                __ATTR_ALIGN__(16) xmm4c4_t res;
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       cdotv_xmm4c4(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<6) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1                
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                   
	                    v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+0] = res;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+1] = res;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+2] = res;
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+3] = res;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+4] = res;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+5] = res;
	                  
	             }          
	       }
	       
	       
	          
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cdotv_xmm4c4_unroll2x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t v2z;
	                __ATTR_ALIGN__(16) xmm4c4_t res;
	                int32_t j,m,m1;
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       cdotv_xmm4c4(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<2) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1                
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                   
	                    v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+0] = res;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+1] = res;
	                                   
	             }          
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cdotv_xmm4c4_rolled(    const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t v2z;
	                __ATTR_ALIGN__(16) xmm4c4_t res;
	                int32_t j;
	                
	               for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1                
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                   
	                    v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    cdotv_xmm4c4(v1x,v1y,v1z,
	                                  v2x,v2y,v2z,
	                                  res);
	                    pres[j+0] = res;
	                  	                                   
	             }          
	       }
	       
	       
	       
	       
	        
	       
	       
	       
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           __m128 cnorm_xmm4c4(const xmm4c4_t vx,
	                                const xmm4c4_t vy,
	                                const xmm4c4_t vz) {
	                                
	                  xmm4c4_t t,cx,cy,cz;
	                  __m128 vs;
	                  cconj_xmm4c4_v2(vx.re,vx.im,&cx.re,&cx.im);
	                  cconj_xmm4c4_v2(vy.re,vy.im,&cy.re,&cy.im);
	                  cconj_xmm4c4_v2(vz.re,vz.im,&cz.re,&cz.im);
	                  cdotv_xmm4c4(vx,vy,vz,cx,cy,cz,t);
	                  vs = _mm_sqrt_ps(t.re);
	                  return (vs);                      
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cnorm_xmm4c4_unroll16x( const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 16;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t vs;
	                int32_t j,m,m1;
	                
	                m = n%16;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<16) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 16) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
#endif
                            v1x = pv1x[j+0];
	                    v1y = pv1y[j+0];
	                    v1z = pv1z[j+0];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+0] = vs;
	                    v1x = pv1x[j+1];
	                    v1y = pv1y[j+1];
	                    v1z = pv1z[j+1];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+1] = vs;
	                    v1x = pv1x[j+2];
	                    v1y = pv1y[j+2];
	                    v1z = pv1z[j+2];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+2] = vs;
	                    v1x = pv1x[j+3];
	                    v1y = pv1y[j+3];
	                    v1z = pv1z[j+3];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+3] = vs;
	                    v1x = pv1x[j+4];
	                    v1y = pv1y[j+4];
	                    v1z = pv1z[j+4];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+4] = vs;
	                    v1x = pv1x[j+5];
	                    v1y = pv1y[j+5];
	                    v1z = pv1z[j+5];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+5] = vs;
	                    v1x = pv1x[j+6];
	                    v1y = pv1y[j+6];
	                    v1z = pv1z[j+6];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+6] = vs;
	                    v1x = pv1x[j+7];
	                    v1y = pv1y[j+7];
	                    v1z = pv1z[j+7];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+7] = vs;
	                    v1x = pv1x[j+8];
	                    v1y = pv1y[j+8];
	                    v1z = pv1z[j+8];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+8] = vs;
	                    v1x = pv1x[j+9];
	                    v1y = pv1y[j+9];
	                    v1z = pv1z[j+9];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+9] = vs;
	                    v1x = pv1x[j+10];
	                    v1y = pv1y[j+10];
	                    v1z = pv1z[j+10];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+10] = vs;
	                    v1x = pv1x[j+11];
	                    v1y = pv1y[j+11];
	                    v1z = pv1z[j+11];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+11] = vs;
	                    v1x = pv1x[j+12];
	                    v1y = pv1y[j+12];
	                    v1z = pv1z[j+12];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+12] = vs;
	                    v1x = pv1x[j+13];
	                    v1y = pv1y[j+13];
	                    v1z = pv1z[j+13];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+13] = vs;
	                    v1x = pv1x[j+14];
	                    v1y = pv1y[j+14];
	                    v1z = pv1z[j+14];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+14] = vs;
	                    v1x = pv1x[j+15];
	                    v1y = pv1y[j+15];
	                    v1z = pv1z[j+15];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+15] = vs;
	                }             
	      }
	      
	      
	         
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cnorm_xmm4c4_unroll10x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 10;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t vs;
	                int32_t j,m,m1;
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<10) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
#endif
                            v1x = pv1x[j+0];
	                    v1y = pv1y[j+0];
	                    v1z = pv1z[j+0];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+0] = vs;
	                    v1x = pv1x[j+1];
	                    v1y = pv1y[j+1];
	                    v1z = pv1z[j+1];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+1] = vs;
	                    v1x = pv1x[j+2];
	                    v1y = pv1y[j+2];
	                    v1z = pv1z[j+2];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+2] = vs;
	                    v1x = pv1x[j+3];
	                    v1y = pv1y[j+3];
	                    v1z = pv1z[j+3];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+3] = vs;
	                    v1x = pv1x[j+4];
	                    v1y = pv1y[j+4];
	                    v1z = pv1z[j+4];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+4] = vs;
	                    v1x = pv1x[j+5];
	                    v1y = pv1y[j+5];
	                    v1z = pv1z[j+5];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+5] = vs;
	                    v1x = pv1x[j+6];
	                    v1y = pv1y[j+6];
	                    v1z = pv1z[j+6];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+6] = vs;
	                    v1x = pv1x[j+7];
	                    v1y = pv1y[j+7];
	                    v1z = pv1z[j+7];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+7] = vs;
	                    v1x = pv1x[j+8];
	                    v1y = pv1y[j+8];
	                    v1z = pv1z[j+8];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+8] = vs;
	                    v1x = pv1x[j+9];
	                    v1y = pv1y[j+9];
	                    v1z = pv1z[j+9];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+9] = vs;
	                  
	                }             
	      }
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cnorm_xmm4c4_unroll6x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 6;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t vs;
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<6) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
#endif
                            v1x = pv1x[j+0];
	                    v1y = pv1y[j+0];
	                    v1z = pv1z[j+0];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+0] = vs;
	                    v1x = pv1x[j+1];
	                    v1y = pv1y[j+1];
	                    v1z = pv1z[j+1];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+1] = vs;
	                    v1x = pv1x[j+2];
	                    v1y = pv1y[j+2];
	                    v1z = pv1z[j+2];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+2] = vs;
	                    v1x = pv1x[j+3];
	                    v1y = pv1y[j+3];
	                    v1z = pv1z[j+3];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+3] = vs;
	                    v1x = pv1x[j+4];
	                    v1y = pv1y[j+4];
	                    v1z = pv1z[j+4];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+4] = vs;
	                    v1x = pv1x[j+5];
	                    v1y = pv1y[j+5];
	                    v1z = pv1z[j+5];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+5] = vs;
	                   	                  
	                }             
	      }
	      
	         
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cnorm_xmm4c4_unroll2x( const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 2;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t vs;
	                int32_t j,m,m1;
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<2) { return;}
	                }                     
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
#endif
                            v1x = pv1x[j+0];
	                    v1y = pv1y[j+0];
	                    v1z = pv1z[j+0];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+0] = vs;
	                    v1x = pv1x[j+1];
	                    v1y = pv1y[j+1];
	                    v1z = pv1z[j+1];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+1] = vs;
	                   	                   	                  
	                }             
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void cnorm_xmm4c4_rolled(    const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 1;
	                __ATTR_ALIGN__(16) xmm4c4_t v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t vs;
	                int32_t j;  
	                
	              	for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
#endif
                            v1x = pv1x[j+0];
	                    v1y = pv1y[j+0];
	                    v1z = pv1z[j+0];
	                    vs  = cnorm_xmm4c4(v1x,v1y,v1z);
	                    pvs[j+0] = vs;
	                   	                   	                   	                  
	                }             
	      }
	      
	      
	      
	                                       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossc_xmm4c4(const xmm4c4_t v1x,
	                                const xmm4c4_t v1y,
	                                const xmm4c4_t v1z,
	                                const xmm4c4_t v2x,
	                                const xmm4c4_t v2y,
	                                const xmm4c4_t v2z,
	                                xmm4c4 & resx,
	                                xmm4c4 & resy,
	                                xmm4c4 & resz) {
	                                
	                 xmm4c4_t t0,t1,t2,t3,t4,t5,t6;
	                 cmul_xmm4c4(v1y.re,v1y.im,v2z.re,
	                              v2z.im,&t0.re,&t0.im); 
	                 cmul_xmm4c4(v1z.re,v1z.im,v2y.re,
	                              v2y.im,&t1.re,&t1.im);
	                 resx.re = _mm_sub_ps(t0.re,t1.re);
	                 resx.im = _mm_sub_ps(t0.im,t1.im);
	                 cmul_xmm4c4(v1z.re,v1z.im,v2x.re,
	                              v2x.im,&t2.re,&t2.im);
	                 cmul_xmm4c4(v1x.re,v1x.im,v2z.re,
	                              v2z.im,&t3.re,&t3.im);
	                 resy.re = _mm_sub_ps(t2.re,t3.re);
	                 resy.im = _mm_sub_ps(t2.im,t3.im);
	                 cmul_xmm4c4(v1x.re,v1x.im,v2y.re,
	                              v2y.im,&t4.re,&t4.im);
	                 cmul_xmm4c4(v1y.re,v1y.im,v2x.re,
	                              v2x.im,&t5.re,&t5.im);    
	                 resz.re = _mm_sub_ps(t4.re,t5.re);
	                 resz.im = _mm_sub_ps(t4.im,t5.im);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossc_xmm4c4_unroll16x( const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 16;
	                __ATTR_ALIGN__(16) xmm4c4_t resx;
	                __ATTR_ALIGN__(16) xmm4c4_t resy;
	                __ATTR_ALIGN__(16) xmm4c4_t resz;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2z;   
	                int32_t j,m,m1;
	                
	                m = n%16;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       scrossc_xmm4c4(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<16) return;
	                }                  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 16) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                    
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+0] = resx;
	                    presy[j+0] = resy;
	                    presz[j+0] = resz;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+1] = resx;
	                    presy[j+1] = resy;
	                    presz[j+1] = resz;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+2] = resx;
	                    presy[j+2] = resy;
	                    presz[j+2] = resz;
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+3] = resx;
	                    presy[j+3] = resy;
	                    presz[j+3] = resz;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+4] = resx;
	                    presy[j+4] = resy;
	                    presz[j+4] = resz;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+5] = resx;
	                    presy[j+5] = resy;
	                    presz[j+5] = resz;
	                    v1x = pv1x[j+6];
	                    v2x = pv2x[j+6];
	                    v1y = pv1y[j+6];
	                    v2y = pv2y[j+6];
	                    v1z = pv1z[j+6];
	                    v2z = pv2z[j+6];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+6] = resx;
	                    presy[j+6] = resy;
	                    presz[j+6] = resz;
	                    v1x = pv1x[j+7];
	                    v2x = pv2x[j+7];
	                    v1y = pv1y[j+7];
	                    v2y = pv2y[j+7];
	                    v1z = pv1z[j+7];
	                    v2z = pv2z[j+7];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+7] = resx;
	                    presy[j+7] = resy;
	                    presz[j+7] = resz;
	                    v1x = pv1x[j+8];
	                    v2x = pv2x[j+8];
	                    v1y = pv1y[j+8];
	                    v2y = pv2y[j+8];
	                    v1z = pv1z[j+8];
	                    v2z = pv2z[j+8];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+8] = resx;
	                    presy[j+8] = resy;
	                    presz[j+8] = resz;
	                    v1x = pv1x[j+9];
	                    v2x = pv2x[j+9];
	                    v1y = pv1y[j+9];
	                    v2y = pv2y[j+9];
	                    v1z = pv1z[j+9];
	                    v2z = pv2z[j+9];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+9] = resx;
	                    presy[j+9] = resy;
	                    presz[j+9] = resz;
	                    v1x = pv1x[j+10];
	                    v2x = pv2x[j+10];
	                    v1y = pv1y[j+10];
	                    v2y = pv2y[j+10];
	                    v1z = pv1z[j+10];
	                    v2z = pv2z[j+10];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+10] = resx;
	                    presy[j+10] = resy;
	                    presz[j+10] = resz;
	                    v1x = pv1x[j+11];
	                    v2x = pv2x[j+11];
	                    v1y = pv1y[j+11];
	                    v2y = pv2y[j+11];
	                    v1z = pv1z[j+11];
	                    v2z = pv2z[j+11];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+11] = resx;
	                    presy[j+11] = resy;
	                    presz[j+11] = resz;
	                    v1x = pv1x[j+12];
	                    v2x = pv2x[j+12];
	                    v1y = pv1y[j+12];
	                    v2y = pv2y[j+12];
	                    v1z = pv1z[j+12];
	                    v2z = pv2z[j+12];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+12] = resx;
	                    presy[j+12] = resy;
	                    presz[j+12] = resz;
	                    v1x = pv1x[j+13];
	                    v2x = pv2x[j+13];
	                    v1y = pv1y[j+13];
	                    v2y = pv2y[j+13];
	                    v1z = pv1z[j+13];
	                    v2z = pv2z[j+13];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+13] = resx;
	                    presy[j+13] = resy;
	                    presz[j+13] = resz;
	                    v1x = pv1x[j+14];
	                    v2x = pv2x[j+14];
	                    v1y = pv1y[j+14];
	                    v2y = pv2y[j+14];
	                    v1z = pv1z[j+14];
	                    v2z = pv2z[j+14];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+14] = resx;
	                    presy[j+14] = resy;
	                    presz[j+14] = resz;
	                    v1x = pv1x[j+15];
	                    v2x = pv2x[j+15];
	                    v1y = pv1y[j+15];
	                    v2y = pv2y[j+15];
	                    v1z = pv1z[j+15];
	                    v2z = pv2z[j+15];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+15] = resx;
	                    presy[j+15] = resy;
	                    presz[j+15] = resz;
	                }          
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossc_xmm4c4_unroll10x( const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                __ATTR_ALIGN__(16) xmm4c4_t resx;
	                __ATTR_ALIGN__(16) xmm4c4_t resy;
	                __ATTR_ALIGN__(16) xmm4c4_t resz;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2z;   
	                int32_t j,m,m1;
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       scrossc_xmm4c4(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<10) return;
	                }                  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                    
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+0] = resx;
	                    presy[j+0] = resy;
	                    presz[j+0] = resz;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+1] = resx;
	                    presy[j+1] = resy;
	                    presz[j+1] = resz;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+2] = resx;
	                    presy[j+2] = resy;
	                    presz[j+2] = resz;
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+3] = resx;
	                    presy[j+3] = resy;
	                    presz[j+3] = resz;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+4] = resx;
	                    presy[j+4] = resy;
	                    presz[j+4] = resz;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+5] = resx;
	                    presy[j+5] = resy;
	                    presz[j+5] = resz;
	                    v1x = pv1x[j+6];
	                    v2x = pv2x[j+6];
	                    v1y = pv1y[j+6];
	                    v2y = pv2y[j+6];
	                    v1z = pv1z[j+6];
	                    v2z = pv2z[j+6];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+6] = resx;
	                    presy[j+6] = resy;
	                    presz[j+6] = resz;
	                    v1x = pv1x[j+7];
	                    v2x = pv2x[j+7];
	                    v1y = pv1y[j+7];
	                    v2y = pv2y[j+7];
	                    v1z = pv1z[j+7];
	                    v2z = pv2z[j+7];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+7] = resx;
	                    presy[j+7] = resy;
	                    presz[j+7] = resz;
	                    v1x = pv1x[j+8];
	                    v2x = pv2x[j+8];
	                    v1y = pv1y[j+8];
	                    v2y = pv2y[j+8];
	                    v1z = pv1z[j+8];
	                    v2z = pv2z[j+8];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+8] = resx;
	                    presy[j+8] = resy;
	                    presz[j+8] = resz;
	                    v1x = pv1x[j+9];
	                    v2x = pv2x[j+9];
	                    v1y = pv1y[j+9];
	                    v2y = pv2y[j+9];
	                    v1z = pv1z[j+9];
	                    v2z = pv2z[j+9];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+9] = resx;
	                    presy[j+9] = resy;
	                    presz[j+9] = resz;
	                 
	                }          
	        }
	        
	        
	          
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossc_xmm4c4_unroll6x(  const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                __ATTR_ALIGN__(16) xmm4c4_t resx;
	                __ATTR_ALIGN__(16) xmm4c4_t resy;
	                __ATTR_ALIGN__(16) xmm4c4_t resz;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2z;   
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       scrossc_xmm4c4(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<6) return;
	                }                  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                    
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+0] = resx;
	                    presy[j+0] = resy;
	                    presz[j+0] = resz;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+1] = resx;
	                    presy[j+1] = resy;
	                    presz[j+1] = resz;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+2] = resx;
	                    presy[j+2] = resy;
	                    presz[j+2] = resz;
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+3] = resx;
	                    presy[j+3] = resy;
	                    presz[j+3] = resz;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+4] = resx;
	                    presy[j+4] = resy;
	                    presz[j+4] = resz;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+5] = resx;
	                    presy[j+5] = resy;
	                    presz[j+5] = resz;
	                  
	                }          
	        }
	        
	        
	         
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossc_xmm4c4_unroll2x(  const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __ATTR_ALIGN__(16) xmm4c4_t resx;
	                __ATTR_ALIGN__(16) xmm4c4_t resy;
	                __ATTR_ALIGN__(16) xmm4c4_t resz;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2z;   
	                int32_t j,m,m1;
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v2x = pv2x[j];
	                       v1y = pv1y[j];
	                       v2y = pv2y[j];
	                       v1z = pv1z[j];
	                       v2z = pv2z[j];
	                       scrossc_xmm4c4(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<2) return;
	                }                  
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                    
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+0] = resx;
	                    presy[j+0] = resy;
	                    presz[j+0] = resz;
	                    v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j+1] = resx;
	                    presy[j+1] = resy;
	                    presz[j+1] = resz;
	                  
	                }          
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossc_xmm4c4_rolled(    const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm4c4_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm4c4_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __ATTR_ALIGN__(16) xmm4c4_t resx;
	                __ATTR_ALIGN__(16) xmm4c4_t resy;
	                __ATTR_ALIGN__(16) xmm4c4_t resz;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v1z;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2x;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2y;
	                __ATTR_ALIGN__(16) xmm4c4_t  v2z;   
	                int32_t j;
	                
	              	for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2z[j+PF_DIST].im,_MM_HINT_NTA);
#endif	                    
                            v1x = pv1x[j];
	                    v2x = pv2x[j];
	                    v1y = pv1y[j];
	                    v2y = pv2y[j];
	                    v1z = pv1z[j];
	                    v2z = pv2z[j];
	                    scrossc_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    resx,resy,resz);
	                    presx[j] = resx;
	                    presy[j] = resy;
	                    presz[j] = resz;
	                 	                  
	                }          
	        }
	        
	        
	        
	        
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_xmm4c4(const __m128 v1x,
	                                const __m128 v1y,
	                                const __m128 v1z,
	                                const __m128 v2x,
	                                const __m128 v2y,
	                                const __m128 v2z,
	                                __m128 * __restrict vcx,
	                                __m128 * __restrict vcy,
	                                __m128 * __restrict vcz) {
	                                
	                *vcx = _mm_fmsub_ps(v1y,v2z,
	                                   _mm_mul_ps(v1x,v2y));
	                *vcy = _mm_fmsub_ps(v1z,v2x,
	                                   _mm_mul_ps(v1x,v2z));
	                *vcz = _mm_fmsub_ps(v1x,v2y,
	                                   _mm_mul_ps(v1y,v2x));
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossv_xmm4c4_unroll16x( const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 16;
	                  __ATTR_ALIGN__(16) __m128 v1x;
	                  __ATTR_ALIGN__(16) __m128 v1y; 
	                  __ATTR_ALIGN__(16) __m128 v1z;
	                  __ATTR_ALIGN__(16) __m128 v2x;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 vcx;
	                  __ATTR_ALIGN__(16) __m128 vcy;
	                  __ATTR_ALIGN__(16) __m128 vcz;
	                  int32_t j,m,m1;
	                  
	                  m = n%16;
	                  if(m!=0) {
	                     for(j = 0;j != m; ++j) {
	                          v1x = pv1x[j];
	                          v2x = pv2x[j];
	                          v1y = pv1y[j];
	                          v2y = pv2y[j];
	                          v1z = pv1z[j];
	                          v2z = pv2z[j];
	                          scrossv_xmm4c4(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<16) return;
	                  }                  
	                  
	                  m1 = m+1;
	                  for(j = m1; j != n; j += 16) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_NTA);
#endif	                    	
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+0] = vcx;
	                    pvcy[j+0] = vcy;
	                    pvcz[j+0] = vcz;
                            v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+1] = vcx;
	                    pvcy[j+1] = vcy;
	                    pvcz[j+1] = vcz;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+2] = vcx;
	                    pvcy[j+2] = vcy;
	                    pvcz[j+2] = vcz;   
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+3] = vcx;
	                    pvcy[j+3] = vcy;
	                    pvcz[j+3] = vcz;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+4] = vcx;
	                    pvcy[j+4] = vcy;
	                    pvcz[j+4] = vcz;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+5] = vcx;
	                    pvcy[j+5] = vcy;
	                    pvcz[j+5] = vcz;
	                    v1x = pv1x[j+6];
	                    v2x = pv2x[j+6];
	                    v1y = pv1y[j+6];
	                    v2y = pv2y[j+6];
	                    v1z = pv1z[j+6];
	                    v2z = pv2z[j+6];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+6] = vcx;
	                    pvcy[j+6] = vcy;
	                    pvcz[j+6] = vcz;
	                    v1x = pv1x[j+7];
	                    v2x = pv2x[j+7];
	                    v1y = pv1y[j+7];
	                    v2y = pv2y[j+7];
	                    v1z = pv1z[j+7];
	                    v2z = pv2z[j+7];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+7] = vcx;
	                    pvcy[j+7] = vcy;
	                    pvcz[j+7] = vcz;  
	                    v1x = pv1x[j+8];
	                    v2x = pv2x[j+8];
	                    v1y = pv1y[j+8];
	                    v2y = pv2y[j+8];
	                    v1z = pv1z[j+8];
	                    v2z = pv2z[j+8];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+8] = vcx;
	                    pvcy[j+8] = vcy;
	                    pvcz[j+8] = vcz;   
	                    v1x = pv1x[j+9];
	                    v2x = pv2x[j+9];
	                    v1y = pv1y[j+9];
	                    v2y = pv2y[j+9];
	                    v1z = pv1z[j+9];
	                    v2z = pv2z[j+9];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+9] = vcx;
	                    pvcy[j+9] = vcy;
	                    pvcz[j+9] = vcz;
	                    v1x = pv1x[j+10];
	                    v2x = pv2x[j+10];
	                    v1y = pv1y[j+10];
	                    v2y = pv2y[j+10];
	                    v1z = pv1z[j+10];
	                    v2z = pv2z[j+10];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+10] = vcx;
	                    pvcy[j+10] = vcy;
	                    pvcz[j+10] = vcz;   
	                    v1x = pv1x[j+11];
	                    v2x = pv2x[j+11];
	                    v1y = pv1y[j+11];
	                    v2y = pv2y[j+11];
	                    v1z = pv1z[j+11];
	                    v2z = pv2z[j+11];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+11] = vcx;
	                    pvcy[j+11] = vcy;
	                    pvcz[j+11] = vcz;
	                    v1x = pv1x[j+12];
	                    v2x = pv2x[j+12];
	                    v1y = pv1y[j+12];
	                    v2y = pv2y[j+12];
	                    v1z = pv1z[j+12];
	                    v2z = pv2z[j+12];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+12] = vcx;
	                    pvcy[j+12] = vcy;
	                    pvcz[j+12] = vcz;   
	                    v1x = pv1x[j+13];
	                    v2x = pv2x[j+13];
	                    v1y = pv1y[j+13];
	                    v2y = pv2y[j+13];
	                    v1z = pv1z[j+13];
	                    v2z = pv2z[j+13];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+13] = vcx;
	                    pvcy[j+13] = vcy;
	                    pvcz[j+13] = vcz;
	                    v1x = pv1x[j+14];
	                    v2x = pv2x[j+14];
	                    v1y = pv1y[j+14];
	                    v2y = pv2y[j+14];
	                    v1z = pv1z[j+14];
	                    v2z = pv2z[j+14];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+14] = vcx;
	                    pvcy[j+14] = vcy;
	                    pvcz[j+14] = vcz; 
	                    v1x = pv1x[j+15];
	                    v2x = pv2x[j+15];
	                    v1y = pv1y[j+15];
	                    v2y = pv2y[j+15];
	                    v1z = pv1z[j+15];
	                    v2z = pv2z[j+15];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+15] = vcx;
	                    pvcy[j+15] = vcy;
	                    pvcz[j+15] = vcz;       
	               }             
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossv_xmm4c4_unroll10x( const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                  __ATTR_ALIGN__(16) __m128 v1x;
	                  __ATTR_ALIGN__(16) __m128 v1y; 
	                  __ATTR_ALIGN__(16) __m128 v1z;
	                  __ATTR_ALIGN__(16) __m128 v2x;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 vcx;
	                  __ATTR_ALIGN__(16) __m128 vcy;
	                  __ATTR_ALIGN__(16) __m128 vcz;
	                  int32_t j,m,m1;
	                  
	                  m = n%10;
	                  if(m!=0) {
	                     for(j = 0;j != m; ++j) {
	                          v1x = pv1x[j];
	                          v2x = pv2x[j];
	                          v1y = pv1y[j];
	                          v2y = pv2y[j];
	                          v1z = pv1z[j];
	                          v2z = pv2z[j];
	                          scrossv_xmm4c4(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<10) return;
	                  }                  
	                  
	                  m1 = m+1;
	                  for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_NTA);
#endif	                    	
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+0] = vcx;
	                    pvcy[j+0] = vcy;
	                    pvcz[j+0] = vcz;
                            v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+1] = vcx;
	                    pvcy[j+1] = vcy;
	                    pvcz[j+1] = vcz;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+2] = vcx;
	                    pvcy[j+2] = vcy;
	                    pvcz[j+2] = vcz;   
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+3] = vcx;
	                    pvcy[j+3] = vcy;
	                    pvcz[j+3] = vcz;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+4] = vcx;
	                    pvcy[j+4] = vcy;
	                    pvcz[j+4] = vcz;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+5] = vcx;
	                    pvcy[j+5] = vcy;
	                    pvcz[j+5] = vcz;
	                    v1x = pv1x[j+6];
	                    v2x = pv2x[j+6];
	                    v1y = pv1y[j+6];
	                    v2y = pv2y[j+6];
	                    v1z = pv1z[j+6];
	                    v2z = pv2z[j+6];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+6] = vcx;
	                    pvcy[j+6] = vcy;
	                    pvcz[j+6] = vcz;
	                    v1x = pv1x[j+7];
	                    v2x = pv2x[j+7];
	                    v1y = pv1y[j+7];
	                    v2y = pv2y[j+7];
	                    v1z = pv1z[j+7];
	                    v2z = pv2z[j+7];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+7] = vcx;
	                    pvcy[j+7] = vcy;
	                    pvcz[j+7] = vcz;  
	                    v1x = pv1x[j+8];
	                    v2x = pv2x[j+8];
	                    v1y = pv1y[j+8];
	                    v2y = pv2y[j+8];
	                    v1z = pv1z[j+8];
	                    v2z = pv2z[j+8];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+8] = vcx;
	                    pvcy[j+8] = vcy;
	                    pvcz[j+8] = vcz;   
	                    v1x = pv1x[j+9];
	                    v2x = pv2x[j+9];
	                    v1y = pv1y[j+9];
	                    v2y = pv2y[j+9];
	                    v1z = pv1z[j+9];
	                    v2z = pv2z[j+9];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+9] = vcx;
	                    pvcy[j+9] = vcy;
	                    pvcz[j+9] = vcz;
	                
	               }             
	         }
	         
	         
	        
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossv_xmm4c4_unroll6x(  const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                  __ATTR_ALIGN__(16) __m128 v1x;
	                  __ATTR_ALIGN__(16) __m128 v1y; 
	                  __ATTR_ALIGN__(16) __m128 v1z;
	                  __ATTR_ALIGN__(16) __m128 v2x;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 vcx;
	                  __ATTR_ALIGN__(16) __m128 vcy;
	                  __ATTR_ALIGN__(16) __m128 vcz;
	                  int32_t j,m,m1;
	                  
	                  m = n%6;
	                  if(m!=0) {
	                     for(j = 0;j != m; ++j) {
	                          v1x = pv1x[j];
	                          v2x = pv2x[j];
	                          v1y = pv1y[j];
	                          v2y = pv2y[j];
	                          v1z = pv1z[j];
	                          v2z = pv2z[j];
	                          scrossv_xmm4c4(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<6) return;
	                  }                  
	                  
	                  m1 = m+1;
	                  for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_NTA);
#endif	                    	
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+0] = vcx;
	                    pvcy[j+0] = vcy;
	                    pvcz[j+0] = vcz;
                            v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+1] = vcx;
	                    pvcy[j+1] = vcy;
	                    pvcz[j+1] = vcz;
	                    v1x = pv1x[j+2];
	                    v2x = pv2x[j+2];
	                    v1y = pv1y[j+2];
	                    v2y = pv2y[j+2];
	                    v1z = pv1z[j+2];
	                    v2z = pv2z[j+2];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+2] = vcx;
	                    pvcy[j+2] = vcy;
	                    pvcz[j+2] = vcz;   
	                    v1x = pv1x[j+3];
	                    v2x = pv2x[j+3];
	                    v1y = pv1y[j+3];
	                    v2y = pv2y[j+3];
	                    v1z = pv1z[j+3];
	                    v2z = pv2z[j+3];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+3] = vcx;
	                    pvcy[j+3] = vcy;
	                    pvcz[j+3] = vcz;
	                    v1x = pv1x[j+4];
	                    v2x = pv2x[j+4];
	                    v1y = pv1y[j+4];
	                    v2y = pv2y[j+4];
	                    v1z = pv1z[j+4];
	                    v2z = pv2z[j+4];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+4] = vcx;
	                    pvcy[j+4] = vcy;
	                    pvcz[j+4] = vcz;
	                    v1x = pv1x[j+5];
	                    v2x = pv2x[j+5];
	                    v1y = pv1y[j+5];
	                    v2y = pv2y[j+5];
	                    v1z = pv1z[j+5];
	                    v2z = pv2z[j+5];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+5] = vcx;
	                    pvcy[j+5] = vcy;
	                    pvcz[j+5] = vcz;
	               
	               }             
	         }
	         
	         
	          
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossv_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                  __ATTR_ALIGN__(16) __m128 v1x;
	                  __ATTR_ALIGN__(16) __m128 v1y; 
	                  __ATTR_ALIGN__(16) __m128 v1z;
	                  __ATTR_ALIGN__(16) __m128 v2x;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 vcx;
	                  __ATTR_ALIGN__(16) __m128 vcy;
	                  __ATTR_ALIGN__(16) __m128 vcz;
	                  int32_t j,m,m1;
	                  
	                  m = n%2;
	                  if(m!=0) {
	                     for(j = 0;j != m; ++j) {
	                          v1x = pv1x[j];
	                          v2x = pv2x[j];
	                          v1y = pv1y[j];
	                          v2y = pv2y[j];
	                          v1z = pv1z[j];
	                          v2z = pv2z[j];
	                          scrossv_xmm4c4(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<2) return;
	                  }                  
	                  
	                  m1 = m+1;
	                  for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_NTA);
#endif	                    	
                            v1x = pv1x[j+0];
	                    v2x = pv2x[j+0];
	                    v1y = pv1y[j+0];
	                    v2y = pv2y[j+0];
	                    v1z = pv1z[j+0];
	                    v2z = pv2z[j+0];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+0] = vcx;
	                    pvcy[j+0] = vcy;
	                    pvcz[j+0] = vcz;
                            v1x = pv1x[j+1];
	                    v2x = pv2x[j+1];
	                    v1y = pv1y[j+1];
	                    v2y = pv2y[j+1];
	                    v1z = pv1z[j+1];
	                    v2z = pv2z[j+1];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j+1] = vcx;
	                    pvcy[j+1] = vcy;
	                    pvcz[j+1] = vcz;
	               	               
	               }             
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void scrossv_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 1;
	                  __ATTR_ALIGN__(16) __m128 v1x;
	                  __ATTR_ALIGN__(16) __m128 v1y; 
	                  __ATTR_ALIGN__(16) __m128 v1z;
	                  __ATTR_ALIGN__(16) __m128 v2x;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 v2y;
	                  __ATTR_ALIGN__(16) __m128 vcx;
	                  __ATTR_ALIGN__(16) __m128 vcy;
	                  __ATTR_ALIGN__(16) __m128 vcz;
	                  int32_t j;
	                  
	                  for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T0);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pv1x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1y[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv1z[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2x[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pv2y[j+PF_DIST],_MM_HINT_NTA);
#endif	                    	
                            v1x = pv1x[j];
	                    v2x = pv2x[j];
	                    v1y = pv1y[j];
	                    v2y = pv2y[j];
	                    v1z = pv1z[j];
	                    v2z = pv2z[j];
	                    scrossv_xmm4c4(v1x,v1y,v1z,
	                                    v2x,v2y,v2z,
	                                    &vcx,&vcy,&vcz);
	                    pvcx[j] = vcx;
	                    pvcy[j] = vcy;
	                    pvcz[j] = vcz;
                           	               	               
	               }             
	         }
	         
	         
	         
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) pv1x,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv1y,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv1z,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv2x,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv2y,
	                                  const float * __restrict __ATTR_ALIGN__(16) pv2z,
	                                  float * __restrict __ATTR_ALIGN__(16) vcx,
	                                  float * __restrict __ATTR_ALIGN__(16) vcy,
	                                  float * __restrict __ATTR_ALIGN__(16) vcz) {
	                      
	                 register __m128 v1x = _mm_load_ps(&pv1x[0]);
	                 register __m128 v1y = _mm_load_ps(&pv1y[0]);
	                 register __m128 v1z = _mm_load_ps(&pv1z[0]);
	                 register __m128 v2x = _mm_load_ps(&pv2x[0]);
	                 register __m128 v2y = _mm_load_ps(&pv2y[0]);
	                 register __m128 v2z = _mm_load_ps(&pv2z[0]);          
	                *vcx = _mm_fmsub_ps(v1y,v2z,
	                                   _mm_mul_ps(v1x,v2y));
	                *vcy = _mm_fmsub_ps(v1z,v2x,
	                                   _mm_mul_ps(v1x,v2z));
	                *vcz = _mm_fmsub_ps(v1x,v2y,
	                                   _mm_mul_ps(v1y,v2x));
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_xmm4c4_u(const float * __restrict pv1x,
	                                  const float * __restrict pv1y,
	                                  const float * __restrict pv1z,
	                                  const float * __restrict pv2x,
	                                  const float * __restrict pv2y,
	                                  const float * __restrict pv2z,
	                                  float * __restrict vcx,
	                                  float * __restrict vcy,
	                                  float * __restrict vcz) {
	                      
	                 register __m128 v1x = _mm_loadu_ps(&pv1x[0]);
	                 register __m128 v1y = _mm_loadu_ps(&pv1y[0]);
	                 register __m128 v1z = _mm_loadu_ps(&pv1z[0]);
	                 register __m128 v2x = _mm_loadu_ps(&pv2x[0]);
	                 register __m128 v2y = _mm_loadu_ps(&pv2y[0]);
	                 register __m128 v2z = _mm_loadu_ps(&pv2z[0]);          
	                *vcx = _mm_fmsub_ps(v1y,v2z,
	                                   _mm_mul_ps(v1x,v2y));
	                *vcy = _mm_fmsub_ps(v1z,v2x,
	                                   _mm_mul_ps(v1x,v2z));
	                *vcz = _mm_fmsub_ps(v1x,v2y,
	                                   _mm_mul_ps(v1y,v2x));
	         }
	         
	         
	         //! Direction Vector spherical [theta,phi] (SIMD data-types)
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_xmm4c4(  const __m128 tht,
	                                  const __m128 phi,
	                                  __m128 * __restrict dvx,
	                                  __m128 * __restrict dvy,
	                                  __m128 * __restrict dvz) {
	                  
	                        
	                register __m128 stht,cphi,sphi,ctht;
	                cphi = _mm_cos_ps(phi);
	                stht = _mm_sin_ps(tht);
	                *dvx = _mm_mul_ps(stht,cphi);
	                sphi = _mm_sin_ps(phi);
	                *dvy = _mm_mul_ps(stht,sphi);
	                ctht = _mm_cos_ps(tht);
	                *dvz = ctht;                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void dir_vec_xmm4c4_unroll16x( const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 16;
	               register __m128 tht;
	               register __m128 phi;
	               __m128 dvx;
	               __m128 dvy;
	               __m128 dvz;
	               int32_t j,m,m1;
	               
	               m = n%16;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_xmm4c4(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<16) {return;}
	               } 
	               
	               m1 = m+1;
	               for(j = m1; j != n; j += 16) 
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
#endif	              
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+0] = dvx;
	                    pdvy[j+0] = dvy;
	                    pdvz[j+0] = dvz;
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+1] = dvx;
	                    pdvy[j+1] = dvy;
	                    pdvz[j+1] = dvz;
	                    tht = ptht[j+2];
	                    phi = pphi[j+2];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+2] = dvx;
	                    pdvy[j+2] = dvy;
	                    pdvz[j+2] = dvz;
	                    tht = ptht[j+3];
	                    phi = pphi[j+3];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+3] = dvx;
	                    pdvy[j+3] = dvy;
	                    pdvz[j+3] = dvz;
	                    tht = ptht[j+4];
	                    phi = pphi[j+4];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+4] = dvx;
	                    pdvy[j+4] = dvy;
	                    pdvz[j+4] = dvz;
	                    tht = ptht[j+5];
	                    phi = pphi[j+5];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+5] = dvx;
	                    pdvy[j+5] = dvy;
	                    pdvz[j+5] = dvz;
	                    tht = ptht[j+6];
	                    phi = pphi[j+6];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+6] = dvx;
	                    pdvy[j+6] = dvy;
	                    pdvz[j+6] = dvz;
	                    tht = ptht[j+7];
	                    phi = pphi[j+7];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+7] = dvx;
	                    pdvy[j+7] = dvy;
	                    pdvz[j+7] = dvz;
	                    tht = ptht[j+8];
	                    phi = pphi[j+8];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+8] = dvx;
	                    pdvy[j+8] = dvy;
	                    pdvz[j+8] = dvz;
	                    tht = ptht[j+9];
	                    phi = pphi[j+9];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+9] = dvx;
	                    pdvy[j+9] = dvy;
	                    pdvz[j+9] = dvz;
	                    tht = ptht[j+10];
	                    phi = pphi[j+10];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+10] = dvx;
	                    pdvy[j+10] = dvy;
	                    pdvz[j+10] = dvz;
	                    tht = ptht[j+11];
	                    phi = pphi[j+11];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+11] = dvx;
	                    pdvy[j+11] = dvy;
	                    pdvz[j+11] = dvz;
	                    tht = ptht[j+12];
	                    phi = pphi[j+12];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+12] = dvx;
	                    pdvy[j+12] = dvy;
	                    pdvz[j+12] = dvz;
	                    tht = ptht[j+13];
	                    phi = pphi[j+13];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+13] = dvx;
	                    pdvy[j+13] = dvy;
	                    pdvz[j+13] = dvz;
	                    tht = ptht[j+14];
	                    phi = pphi[j+14];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+14] = dvx;
	                    pdvy[j+14] = dvy;
	                    pdvz[j+14] = dvz;
	                    tht = ptht[j+15];
	                    phi = pphi[j+15];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+15] = dvx;
	                    pdvy[j+15] = dvy;
	                    pdvz[j+15] = dvz;
	               }                                  
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void dir_vec_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 10;
	               register __m128 tht;
	               register __m128 phi;
	               __m128 dvx;
	               __m128 dvy;
	               __m128 dvz;
	               int32_t j,m,m1;
	               
	               m = n%10;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_xmm4c4(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<10) {return;}
	               } 
	               
	               m1 = m+1;
	               for(j = m1; j != n; j += 10) 
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
#endif	              
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+0] = dvx;
	                    pdvy[j+0] = dvy;
	                    pdvz[j+0] = dvz;
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+1] = dvx;
	                    pdvy[j+1] = dvy;
	                    pdvz[j+1] = dvz;
	                    tht = ptht[j+2];
	                    phi = pphi[j+2];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+2] = dvx;
	                    pdvy[j+2] = dvy;
	                    pdvz[j+2] = dvz;
	                    tht = ptht[j+3];
	                    phi = pphi[j+3];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+3] = dvx;
	                    pdvy[j+3] = dvy;
	                    pdvz[j+3] = dvz;
	                    tht = ptht[j+4];
	                    phi = pphi[j+4];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+4] = dvx;
	                    pdvy[j+4] = dvy;
	                    pdvz[j+4] = dvz;
	                    tht = ptht[j+5];
	                    phi = pphi[j+5];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+5] = dvx;
	                    pdvy[j+5] = dvy;
	                    pdvz[j+5] = dvz;
	                    tht = ptht[j+6];
	                    phi = pphi[j+6];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+6] = dvx;
	                    pdvy[j+6] = dvy;
	                    pdvz[j+6] = dvz;
	                    tht = ptht[j+7];
	                    phi = pphi[j+7];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+7] = dvx;
	                    pdvy[j+7] = dvy;
	                    pdvz[j+7] = dvz;
	                    tht = ptht[j+8];
	                    phi = pphi[j+8];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+8] = dvx;
	                    pdvy[j+8] = dvy;
	                    pdvz[j+8] = dvz;
	                    tht = ptht[j+9];
	                    phi = pphi[j+9];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+9] = dvx;
	                    pdvy[j+9] = dvy;
	                    pdvz[j+9] = dvz;
	                  
	               }                                  
	       }
	       
	           
	       
	        
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  static inline
	           void dir_vec_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 6;
	               register __m128 tht;
	               register __m128 phi;
	               __m128 dvx;
	               __m128 dvy;
	               __m128 dvz;
	               int32_t j,m,m1;
	               
	               m = n%6;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_xmm4c4(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<6) {return;}
	               } 
	               
	               m1 = m+1;
	               for(j = m1; j != n; j += 6) 
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
#endif	              
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+0] = dvx;
	                    pdvy[j+0] = dvy;
	                    pdvz[j+0] = dvz;
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+1] = dvx;
	                    pdvy[j+1] = dvy;
	                    pdvz[j+1] = dvz;
	                    tht = ptht[j+2];
	                    phi = pphi[j+2];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+2] = dvx;
	                    pdvy[j+2] = dvy;
	                    pdvz[j+2] = dvz;
	                    tht = ptht[j+3];
	                    phi = pphi[j+3];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+3] = dvx;
	                    pdvy[j+3] = dvy;
	                    pdvz[j+3] = dvz;
	                    tht = ptht[j+4];
	                    phi = pphi[j+4];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+4] = dvx;
	                    pdvy[j+4] = dvy;
	                    pdvz[j+4] = dvz;
	                    tht = ptht[j+5];
	                    phi = pphi[j+5];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+5] = dvx;
	                    pdvy[j+5] = dvy;
	                    pdvz[j+5] = dvz;
	                  	                  
	               }                                  
	       }
	       
	       
	          
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void dir_vec_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 2;
	               register __m128 tht;
	               register __m128 phi;
	                __m128 dvx;
	                __m128 dvy;
	                __m128 dvz;
	               int32_t j,m,m1;
	               
	               m = n%2;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_xmm4c4(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<2) {return;}
	               } 
	               
	               m1 = m+1;
	               for(j = m1; j != n; j += 2) 
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
#endif	              
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+0] = dvx;
	                    pdvy[j+0] = dvy;
	                    pdvz[j+0] = dvz;
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j+1] = dvx;
	                    pdvy[j+1] = dvy;
	                    pdvz[j+1] = dvz;
	                                    	                  
	               }                                  
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void dir_vec_xmm4c4_rolled(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 1;
	               register __m128 tht;
	               register __m128 phi;
	                __m128 dvx;
	                __m128 dvy;
	                __m128 dvz;
	               int32_t j;
	               
	              for(j = 0; j != n; ++j) 
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
#endif	              
                            tht = ptht[j];
	                    phi = pphi[j];
	                    dir_vec_xmm4c4(tht,phi,
	                                    &dvx,&dvy,&dvz);
	                    pdvx[j] = dvx;
	                    pdvy[j] = dvy;
	                    pdvz[j] = dvz;
	                   
	                                    	                  
	               }                                  
	       }
	       
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) ptht,
	                                  const float * __restrict __ATTR_ALIGN__(16) pphi,
	                                  float * __restrict __ATTR_ALIGN__(16) dvx,
	                                  float * __restrict __ATTR_ALIGN__(16) dvy,
	                                  float * __restrict __ATTR_ALIGN__(16) dvz) {
	                  
	                register __m128 tht = _mm_load_ps(&ptht[0]);
	                register __m128 phi = _mm_load_ps(&pphi[0]);              
	                register __m128 stht,cphi,sphi,ctht;
	                cphi = _mm_cos_ps(phi);
	                stht = _mm_sin_ps(tht);
	                _mm_store_ps(&dvx[0] , _mm_mul_ps(stht,cphi));
	                sphi = _mm_sin_ps(phi);
	                _mm_store_ps(&dvy[0] , _mm_mul_ps(stht,sphi));
	                ctht = _mm_cos_ps(tht);
	                _mm_store_ps(&dvz[0] , ctht);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_xmm4c4_u(const float * __restrict  ptht,
	                                  const float * __restrict  pphi,
	                                  float * __restrict  dvx,
	                                  float * __restrict  dvy,
	                                  float * __restrict  dvz) {
	                  
	                register __m128 tht = _mm_loadu_ps(&ptht[0]);
	                register __m128 phi = _mm_loadu_ps(&pphi[0]);              
	                register __m128 stht,cphi,sphi,ctht;
	                cphi = _mm_cos_ps(phi);
	                stht = _mm_sin_ps(tht);
	                _mm_storeu_ps(&dvx[0] , _mm_mul_ps(stht,cphi));
	                sphi = _mm_sin_ps(phi);
	                _mm_storeu_ps(&dvy[0] , _mm_mul_ps(stht,sphi));
	                ctht = _mm_cos_ps(tht);
	                _mm_storeu_ps(&dvz[0] , ctht);                       
	        }
	        
	        
	         //! Polarization Vector of plane-wave propagating into direction computed by
                 //! dir_vector_xmmxrx (SIMD data-types)
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_xmm4c4(const __m128 tht,
	                                const __m128 phi,
	                                const __m128 psi,
	                                __m128 * __restrict pvx,
	                                __m128 * __restrict pvy,
	                                __m128 * __restrict pvz) {
	                 
	                using namespace gms::math               
	                register __m128 cpsi,cphi,spsi,sphi,t0;
	                cpsi = _mm_cos_ps(psi);
	                cphi = _mm_cos_ps(phi);
	                spsi = _mm_sin_ps(psi);
	                sphi = _mm_sin_ps(phi);
	                t0   = _mm_mul_ps(spsi,_mm_cos_ps(tht));
	                *pvx = _mm_fmsub_ps(cpsi,sphi,
	                                   _mm_mul_ps(t0,cphi));
	                *pvy = _mm_fmsub_ps(negate_xmm4r4(cpsi),cphi,
	                                                    _mm_mul_ps(t0,sphi));
	                *pvz = _mm_mul_ps(spsi,_mm_sin_ps(tht));                         
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void pol_vec_xmm4c4_unroll16x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 16;
	                register __m128 tht;
	                register __m128 phi;
	                register __m128 psi;
	                 __m128 pvx;
	                 __m128 pvy;
	                 __m128 pvz;
	                int32_t j,m,m1;    
	                
	                m = n%16;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_xmm4c4(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<16) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 16) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);	    
#endif	        	    
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    psi = ppsi[j+0];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+0] = pvx;
	                    ppvy[j+0] = pvy;
	                    ppvz[j+0] = pvz;  
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    psi = ppsi[j+1];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+1] = pvx;
	                    ppvy[j+1] = pvy;
	                    ppvz[j+1] = pvz;   
	                    tht = ptht[j+2];
	                    phi = pphi[j+2];
	                    psi = ppsi[j+2];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+2] = pvx;
	                    ppvy[j+2] = pvy;
	                    ppvz[j+2] = pvz;
	                    tht = ptht[j+3];
	                    phi = pphi[j+3];
	                    psi = ppsi[j+3];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+3] = pvx;
	                    ppvy[j+3] = pvy;
	                    ppvz[j+3] = pvz;  
	                    tht = ptht[j+4];
	                    phi = pphi[j+4];
	                    psi = ppsi[j+4];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+4] = pvx;
	                    ppvy[j+4] = pvy;
	                    ppvz[j+4] = pvz;
	                    tht = ptht[j+5];
	                    phi = pphi[j+5];
	                    psi = ppsi[j+5];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+5] = pvx;
	                    ppvy[j+5] = pvy;
	                    ppvz[j+5] = pvz;  
	                    tht = ptht[j+6];
	                    phi = pphi[j+6];
	                    psi = ppsi[j+6];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+6] = pvx;
	                    ppvy[j+6] = pvy;
	                    ppvz[j+6] = pvz;  
	                    tht = ptht[j+7];
	                    phi = pphi[j+7];
	                    psi = ppsi[j+7];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+7] = pvx;
	                    ppvy[j+7] = pvy;
	                    ppvz[j+7] = pvz;  
	                    tht = ptht[j+8];
	                    phi = pphi[j+8];
	                    psi = ppsi[j+8];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+8] = pvx;
	                    ppvy[j+8] = pvy;
	                    ppvz[j+8] = pvz;  
	                    tht = ptht[j+9];
	                    phi = pphi[j+9];
	                    psi = ppsi[j+9];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+9] = pvx;
	                    ppvy[j+9] = pvy;
	                    ppvz[j+9] = pvz;   
	                    tht = ptht[j+10];
	                    phi = pphi[j+10];
	                    psi = ppsi[j+10];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+10] = pvx;
	                    ppvy[j+10] = pvy;
	                    ppvz[j+10] = pvz;  
	                    tht = ptht[j+11];
	                    phi = pphi[j+11];
	                    psi = ppsi[j+11];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+11] = pvx;
	                    ppvy[j+11] = pvy;
	                    ppvz[j+11] = pvz;   
	                    tht = ptht[j+12];
	                    phi = pphi[j+12];
	                    psi = ppsi[j+12];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+12] = pvx;
	                    ppvy[j+12] = pvy;
	                    ppvz[j+12] = pvz; 
	                    tht = ptht[j+13];
	                    phi = pphi[j+13];
	                    psi = ppsi[j+13];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+13] = pvx;
	                    ppvy[j+13] = pvy;
	                    ppvz[j+13] = pvz;
	                    tht = ptht[j+14];
	                    phi = pphi[j+14];
	                    psi = ppsi[j+14];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+14] = pvx;
	                    ppvy[j+14] = pvy;
	                    ppvz[j+14] = pvz;   
	                    tht = ptht[j+15];
	                    phi = pphi[j+15];
	                    psi = ppsi[j+15];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+15] = pvx;
	                    ppvy[j+15] = pvy;
	                    ppvz[j+15] = pvz;                     
	                }            
	      }
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void pol_vec_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                register __m128 tht;
	                register __m128 phi;
	                register __m128 psi;
	                __m128 pvx;
	                __m128 pvy;
	                __m128 pvz;
	                int32_t j,m,m1;    
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_xmm4c4(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<10) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);	    
#endif	        	    
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    psi = ppsi[j+0];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+0] = pvx;
	                    ppvy[j+0] = pvy;
	                    ppvz[j+0] = pvz;  
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    psi = ppsi[j+1];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+1] = pvx;
	                    ppvy[j+1] = pvy;
	                    ppvz[j+1] = pvz;   
	                    tht = ptht[j+2];
	                    phi = pphi[j+2];
	                    psi = ppsi[j+2];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+2] = pvx;
	                    ppvy[j+2] = pvy;
	                    ppvz[j+2] = pvz;
	                    tht = ptht[j+3];
	                    phi = pphi[j+3];
	                    psi = ppsi[j+3];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+3] = pvx;
	                    ppvy[j+3] = pvy;
	                    ppvz[j+3] = pvz;  
	                    tht = ptht[j+4];
	                    phi = pphi[j+4];
	                    psi = ppsi[j+4];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+4] = pvx;
	                    ppvy[j+4] = pvy;
	                    ppvz[j+4] = pvz;
	                    tht = ptht[j+5];
	                    phi = pphi[j+5];
	                    psi = ppsi[j+5];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+5] = pvx;
	                    ppvy[j+5] = pvy;
	                    ppvz[j+5] = pvz;  
	                    tht = ptht[j+6];
	                    phi = pphi[j+6];
	                    psi = ppsi[j+6];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+6] = pvx;
	                    ppvy[j+6] = pvy;
	                    ppvz[j+6] = pvz;  
	                    tht = ptht[j+7];
	                    phi = pphi[j+7];
	                    psi = ppsi[j+7];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+7] = pvx;
	                    ppvy[j+7] = pvy;
	                    ppvz[j+7] = pvz;  
	                    tht = ptht[j+8];
	                    phi = pphi[j+8];
	                    psi = ppsi[j+8];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+8] = pvx;
	                    ppvy[j+8] = pvy;
	                    ppvz[j+8] = pvz;  
	                    tht = ptht[j+9];
	                    phi = pphi[j+9];
	                    psi = ppsi[j+9];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+9] = pvx;
	                    ppvy[j+9] = pvy;
	                    ppvz[j+9] = pvz;   
	                  
	                }            
	      }
	      
	      
	          
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  static inline
	           void pol_vec_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                register __m128 tht;
	                register __m128 phi;
	                register __m128 psi;
	                 __m128 pvx;
	                 __m128 pvy;
	                 __m128 pvz;
	                int32_t j,m,m1;    
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_xmm4c4(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<6) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);	    
#endif	        	    
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    psi = ppsi[j+0];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+0] = pvx;
	                    ppvy[j+0] = pvy;
	                    ppvz[j+0] = pvz;  
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    psi = ppsi[j+1];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+1] = pvx;
	                    ppvy[j+1] = pvy;
	                    ppvz[j+1] = pvz;   
	                    tht = ptht[j+2];
	                    phi = pphi[j+2];
	                    psi = ppsi[j+2];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+2] = pvx;
	                    ppvy[j+2] = pvy;
	                    ppvz[j+2] = pvz;
	                    tht = ptht[j+3];
	                    phi = pphi[j+3];
	                    psi = ppsi[j+3];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+3] = pvx;
	                    ppvy[j+3] = pvy;
	                    ppvz[j+3] = pvz;  
	                    tht = ptht[j+4];
	                    phi = pphi[j+4];
	                    psi = ppsi[j+4];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+4] = pvx;
	                    ppvy[j+4] = pvy;
	                    ppvz[j+4] = pvz;
	                    tht = ptht[j+5];
	                    phi = pphi[j+5];
	                    psi = ppsi[j+5];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+5] = pvx;
	                    ppvy[j+5] = pvy;
	                    ppvz[j+5] = pvz;  
	                                   
	                }            
	                
	          }
	          
	          
	        
	      
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void pol_vec_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                register __m128 tht;
	                register __m128 phi;
	                register __m128 psi;
	                register __m128 pvx;
	                register __m128 pvy;
	                register __m128 pvz;
	                int32_t j,m,m1;    
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_xmm4c4(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<2) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);	    
#endif	        	    
                            tht = ptht[j+0];
	                    phi = pphi[j+0];
	                    psi = ppsi[j+0];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+0] = pvx;
	                    ppvy[j+0] = pvy;
	                    ppvz[j+0] = pvz;  
	                    tht = ptht[j+1];
	                    phi = pphi[j+1];
	                    psi = ppsi[j+1];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j+1] = pvx;
	                    ppvy[j+1] = pvy;
	                    ppvz[j+1] = pvz;   
	                  	                                   
	                }            
	                
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void pol_vec_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 1;
	                register __m128 tht;
	                register __m128 phi;
	                register __m128 psi;
	                __m128 pvx;
	                __m128 pvy;
	                __m128 pvz;
	                int32_t j;    
	                
	             	for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);	                   
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);	    
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);	    
#endif	        	    
                            tht = ptht[j];
	                    phi = pphi[j];
	                    psi = ppsi[j];
	                    pol_vec_xmm4c4(tht,phi,psi,
	                                    &pvx,&pvy,&pvz);
	                    ppvx[j] = pvx;
	                    ppvy[j] = pvy;
	                    ppvz[j] = pvz;  
	                  	                  	                                   
	                }            
	                
	          }
	          
	      
	      
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) ptht,
	                                  const float * __restrict __ATTR_ALIGN__(16) pphi,
	                                  const float * __restrict __ATTR_ALIGN__(16) psi,
	                                  float * __restrict __ATTR_ALIGN__(16) pvx,
	                                  float * __restrict __ATTR_ALIGN__(16) pvy,
	                                  float * __restrict __ATTR_ALIGN__(16) pvz) {
	                 
	                 using namespace gms::math     
	                register __m128 tht = _mm_load_ps(&ptht[0]);
	                register __m128 phi = _mm_load_ps(&pphi[0]);  
	                register __m128 psi = _mm_load_ps(&ppsi[0]);           
	                register __m128 cpsi,cphi,spsi,sphi,t0;
	                cpsi = _mm_cos_ps(psi);
	                cphi = _mm_cos_ps(phi);
	                spsi = _mm_sin_ps(psi);
	                sphi = _mm_sin_ps(phi);
	                t0   = _mm_mul_ps(spsi,_mm_cos_ps(tht));
	                _mm_store_ps(&pvx[0] ,_mm_fmsub_ps(cpsi,sphi,
	                                   _mm_mul_ps(t0,cphi)));
	                _mm_store_ps(&pvy[0] ,_mm_fmsub_ps(negate_xmm4r4(cpsi),cphi,
	                                                    _mm_mul_ps(t0,sphi)));
	                _mm_store_ps(&pvz[0] ,_mm_mul_ps(spsi,_mm_sin_ps(tht)));                         
	      } 
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_xmm4c4_u(const float * __restrict  ptht,
	                                  const float * __restrict  pphi,
	                                  const float * __restrict  psi,
	                                  float * __restrict  pvx,
	                                  float * __restrict  pvy,
	                                  float * __restrict  pvz) {
	                 
	                  using namespace gms::math    
	                register __m128 tht = _mm_loadu_ps(&ptht[0]);
	                register __m128 phi = _mm_loadu_ps(&pphi[0]);  
	                register __m128 psi = _mm_loadu_ps(&ppsi[0]);           
	                register __m128 cpsi,cphi,spsi,sphi,t0;
	                cpsi = _mm_cos_ps(psi);
	                cphi = _mm_cos_ps(phi);
	                spsi = _mm_sin_ps(psi);
	                sphi = _mm_sin_ps(phi);
	                t0   = _mm_mul_ps(spsi,_mm_cos_ps(tht));
	                _mm_storeu_ps(&pvx[0] ,_mm_fmsub_ps(cpsi,sphi,
	                                   _mm_mul_ps(t0,cphi)));
	                _mm_storeu_ps(&pvy[0] ,_mm_fmsub_ps(negate_xmm4r4(cpsi),cphi,
	                                                    _mm_mul_ps(t0,sphi)));
	                _mm_storeu_ps(&pvz[0] ,_mm_mul_ps(spsi,_mm_sin_ps(tht)));                         
	      } 
	      
	      
	      /*
	           
     ! Vectorized Electric-field at 16 points 'R'
     ! vpol -- vector of vertical polarization at point 'R'
     ! vdir -- direction vector
     ! vr   -- vector radius r
     ! Exyz -- resulting electrical field (3D) at sixteen points 'R', i.e. R(xyz), x0-x15,y0-y15,z0-z15
	      */
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_xmm4c4(const __m128 vpolx,
	                                 const __m128 vpoly,
	                                 const __m128 vpolz,
	                                 const __m128 vdirx,
	                                 const __m128 vdiry,
	                                 const __m128 vdirz,
	                                 const __m128 vrx,
	                                 const __m128 vry,
	                                 const __m128 vrz,
	                                 const xmm4c4_t k,
	                                 xmm4c4_t & H_x,
	                                 xmm4c4_t & H_y,
	                                 xmm4c4_t & H_z) {
	               
	               	register __m128 dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_xmm4c4(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm_set1_ps(1.0);
	                ir = _mm_setzero_ps();
	                cmul_xmm4c4(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm_mul_ps(dp,cer);
	                cei = _mm_mul_ps(dp,cei);
	                cexp_xmm4c4(cer,cei,&expr,&expi);
	                H_x.re = _mm_mul_ps(vpolx,expr);
	                H_x.im = _mm_mul_ps(vpolx,expi);
	                H_y.re = _mm_mul_ps(vpoly,expr);
	                H_y.im = _mm_mul_ps(vpoly,expi);
	                H_z.re = _mm_mul_ps(vpolz,expr);
	                H_z.im = _mm_mul_ps(vpolz,expi);
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  static inline
	           void H_XYZ_VP_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                __ATTR_ALIGN__(16) xmm4c4_t k;
	                __ATTR_ALIGN__(16) xmm4c4_t H_x;
	                __ATTR_ALIGN__(16) xmm4c4_t H_y;
	                __ATTR_ALIGN__(16) xmm4c4_t H_z;
	                register __m128 vpolx;
	                register __m128 vpoly;
	                register __m128 vpolz;
	                register __m128 vdirx;
	                register __m128 vdiry;
	                register __m128 vdirz;
	                register __m128 vrx;
	                register __m128 vry;
	                register __m128 vrz;
	                int32_t j,m,m1;
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       vpolx = pvpolx[j];
	                       vpoly = pvpoly[j];
	                       vpolz = pvpolz[j];
	                       vdirx = pvdirx[j];
	                       vdiry = pvdiry[j];
	                       vdirz = pvdirz[j];
	                       vrx   = pvrx[j];
	                       vry   = pvry[j];
	                       vrz   = pvrz[j];
	                       k     = pk[j];
	                       H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                        vdirx,vdiry,vdirz,
	                                        vrx,vry,vrz,
	                                        H_x,H_y,H_z);
	                       pH_x[j] = H_x;
	                       pH_y[j] = H_y;
	                       pH_z[j] = H_z;   
	                   }
	                   if(n<10) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           
                            vpolx = pvpolx[j+0];
	                    vpoly = pvpoly[j+0];
	                    vpolz = pvpolz[j+0];
	                    vdirx = pvdirx[j+0];
	                    vdiry = pvdiry[j+0];
	                    vdirz = pvdirz[j+0];
	                    vrx   = pvrx[j+0];
	                    vry   = pvry[j+0];
	                    vrz   = pvrz[j+0];
	                    k     = pk[j+0];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+0] = H_x;
	                    pH_y[j+0] = H_y;
	                    pH_z[j+0] = H_z;   
                            vpolx = pvpolx[j+1];
	                    vpoly = pvpoly[j+1];
	                    vpolz = pvpolz[j+1];
	                    vdirx = pvdirx[j+1];
	                    vdiry = pvdiry[j+1];
	                    vdirz = pvdirz[j+1];
	                    vrx   = pvrx[j+1];
	                    vry   = pvry[j+1];
	                    vrz   = pvrz[j+1];
	                    k     = pk[j+1];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+1] = H_x;
	                    pH_y[j+1] = H_y;
	                    pH_z[j+1] = H_z;
	                    vpolx = pvpolx[j+2];
	                    vpoly = pvpoly[j+2];
	                    vpolz = pvpolz[j+2];
	                    vdirx = pvdirx[j+2];
	                    vdiry = pvdiry[j+2];
	                    vdirz = pvdirz[j+2];
	                    vrx   = pvrx[j+2];
	                    vry   = pvry[j+2];
	                    vrz   = pvrz[j+2];
	                    k     = pk[j+2];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+2] = H_x;
	                    pH_y[j+2] = H_y;
	                    pH_z[j+2] = H_z;   
	                    vpolx = pvpolx[j+3];
	                    vpoly = pvpoly[j+3];
	                    vpolz = pvpolz[j+3];
	                    vdirx = pvdirx[j+3];
	                    vdiry = pvdiry[j+3];
	                    vdirz = pvdirz[j+3];
	                    vrx   = pvrx[j+3];
	                    vry   = pvry[j+3];
	                    vrz   = pvrz[j+3];
	                    k     = pk[j+3];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+3] = H_x;
	                    pH_y[j+3] = H_y;
	                    pH_z[j+3] = H_z;  
	                    vpolx = pvpolx[j+4];
	                    vpoly = pvpoly[j+4];
	                    vpolz = pvpolz[j+4];
	                    vdirx = pvdirx[j+4];
	                    vdiry = pvdiry[j+4];
	                    vdirz = pvdirz[j+4];
	                    vrx   = pvrx[j+4];
	                    vry   = pvry[j+4];
	                    vrz   = pvrz[j+4];
	                    k     = pk[j+4];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+4] = H_x;
	                    pH_y[j+4] = H_y;
	                    pH_z[j+4] = H_z;
	                    vpolx = pvpolx[j+5];
	                    vpoly = pvpoly[j+5];
	                    vpolz = pvpolz[j+5];
	                    vdirx = pvdirx[j+5];
	                    vdiry = pvdiry[j+5];
	                    vdirz = pvdirz[j+5];
	                    vrx   = pvrx[j+5];
	                    vry   = pvry[j+5];
	                    vrz   = pvrz[j+5];
	                    k     = pk[j+5];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+5] = H_x;
	                    pH_y[j+5] = H_y;
	                    pH_z[j+5] = H_z;  
	                    vpolx = pvpolx[j+6];
	                    vpoly = pvpoly[j+6];
	                    vpolz = pvpolz[j+6];
	                    vdirx = pvdirx[j+6];
	                    vdiry = pvdiry[j+6];
	                    vdirz = pvdirz[j+6];
	                    vrx   = pvrx[j+6];
	                    vry   = pvry[j+6];
	                    vrz   = pvrz[j+6];
	                    k     = pk[j+6];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+6] = H_x;
	                    pH_y[j+6] = H_y;
	                    pH_z[j+6] = H_z;
	                    vpolx = pvpolx[j+7];
	                    vpoly = pvpoly[j+7];
	                    vpolz = pvpolz[j+7];
	                    vdirx = pvdirx[j+7];
	                    vdiry = pvdiry[j+7];
	                    vdirz = pvdirz[j+7];
	                    vrx   = pvrx[j+7];
	                    vry   = pvry[j+7];
	                    vrz   = pvrz[j+7];
	                    k     = pk[j+7];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+7] = H_x;
	                    pH_y[j+7] = H_y;
	                    pH_z[j+7] = H_z;  
	                    vpolx = pvpolx[j+8];
	                    vpoly = pvpoly[j+8];
	                    vpolz = pvpolz[j+8];
	                    vdirx = pvdirx[j+8];
	                    vdiry = pvdiry[j+8];
	                    vdirz = pvdirz[j+8];
	                    vrx   = pvrx[j+8];
	                    vry   = pvry[j+8];
	                    vrz   = pvrz[j+8];
	                    k     = pk[j+8];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+8] = H_x;
	                    pH_y[j+8] = H_y;
	                    pH_z[j+8] = H_z;
	                    vpolx = pvpolx[j+9];
	                    vpoly = pvpoly[j+9];
	                    vpolz = pvpolz[j+9];
	                    vdirx = pvdirx[j+9];
	                    vdiry = pvdiry[j+9];
	                    vdirz = pvdirz[j+9];
	                    vrx   = pvrx[j+9];
	                    vry   = pvry[j+9];
	                    vrz   = pvrz[j+9];
	                    k     = pk[j+9];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+9] = H_x;
	                    pH_y[j+9] = H_y;
	                    pH_z[j+9] = H_z;         
	                }                
	      }
	      
	      
	         
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  static inline
	           void H_XYZ_VP_xmm4c4_unroll6x( const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                __ATTR_ALIGN__(16) xmm4c4_t k;
	                __ATTR_ALIGN__(16) xmm4c4_t H_x;
	                __ATTR_ALIGN__(16) xmm4c4_t H_y;
	                __ATTR_ALIGN__(16) xmm4c4_t H_z;
	                register __m128 vpolx;
	                register __m128 vpoly;
	                register __m128 vpolz;
	                register __m128 vdirx;
	                register __m128 vdiry;
	                register __m128 vdirz;
	                register __m128 vrx;
	                register __m128 vry;
	                register __m128 vrz;
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       vpolx = pvpolx[j];
	                       vpoly = pvpoly[j];
	                       vpolz = pvpolz[j];
	                       vdirx = pvdirx[j];
	                       vdiry = pvdiry[j];
	                       vdirz = pvdirz[j];
	                       vrx   = pvrx[j];
	                       vry   = pvry[j];
	                       vrz   = pvrz[j];
	                       k     = pk[j];
	                       H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                        vdirx,vdiry,vdirz,
	                                        vrx,vry,vrz,
	                                        H_x,H_y,H_z);
	                       pH_x[j] = H_x;
	                       pH_y[j] = H_y;
	                       pH_z[j] = H_z;   
	                   }
	                   if(n<6) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           
                            vpolx = pvpolx[j+0];
	                    vpoly = pvpoly[j+0];
	                    vpolz = pvpolz[j+0];
	                    vdirx = pvdirx[j+0];
	                    vdiry = pvdiry[j+0];
	                    vdirz = pvdirz[j+0];
	                    vrx   = pvrx[j+0];
	                    vry   = pvry[j+0];
	                    vrz   = pvrz[j+0];
	                    k     = pk[j+0];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+0] = H_x;
	                    pH_y[j+0] = H_y;
	                    pH_z[j+0] = H_z;   
                            vpolx = pvpolx[j+1];
	                    vpoly = pvpoly[j+1];
	                    vpolz = pvpolz[j+1];
	                    vdirx = pvdirx[j+1];
	                    vdiry = pvdiry[j+1];
	                    vdirz = pvdirz[j+1];
	                    vrx   = pvrx[j+1];
	                    vry   = pvry[j+1];
	                    vrz   = pvrz[j+1];
	                    k     = pk[j+1];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+1] = H_x;
	                    pH_y[j+1] = H_y;
	                    pH_z[j+1] = H_z;
	                    vpolx = pvpolx[j+2];
	                    vpoly = pvpoly[j+2];
	                    vpolz = pvpolz[j+2];
	                    vdirx = pvdirx[j+2];
	                    vdiry = pvdiry[j+2];
	                    vdirz = pvdirz[j+2];
	                    vrx   = pvrx[j+2];
	                    vry   = pvry[j+2];
	                    vrz   = pvrz[j+2];
	                    k     = pk[j+2];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+2] = H_x;
	                    pH_y[j+2] = H_y;
	                    pH_z[j+2] = H_z;   
	                    vpolx = pvpolx[j+3];
	                    vpoly = pvpoly[j+3];
	                    vpolz = pvpolz[j+3];
	                    vdirx = pvdirx[j+3];
	                    vdiry = pvdiry[j+3];
	                    vdirz = pvdirz[j+3];
	                    vrx   = pvrx[j+3];
	                    vry   = pvry[j+3];
	                    vrz   = pvrz[j+3];
	                    k     = pk[j+3];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+3] = H_x;
	                    pH_y[j+3] = H_y;
	                    pH_z[j+3] = H_z;  
	                    vpolx = pvpolx[j+4];
	                    vpoly = pvpoly[j+4];
	                    vpolz = pvpolz[j+4];
	                    vdirx = pvdirx[j+4];
	                    vdiry = pvdiry[j+4];
	                    vdirz = pvdirz[j+4];
	                    vrx   = pvrx[j+4];
	                    vry   = pvry[j+4];
	                    vrz   = pvrz[j+4];
	                    k     = pk[j+4];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+4] = H_x;
	                    pH_y[j+4] = H_y;
	                    pH_z[j+4] = H_z;
	                    vpolx = pvpolx[j+5];
	                    vpoly = pvpoly[j+5];
	                    vpolz = pvpolz[j+5];
	                    vdirx = pvdirx[j+5];
	                    vdiry = pvdiry[j+5];
	                    vdirz = pvdirz[j+5];
	                    vrx   = pvrx[j+5];
	                    vry   = pvry[j+5];
	                    vrz   = pvrz[j+5];
	                    k     = pk[j+5];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+5] = H_x;
	                    pH_y[j+5] = H_y;
	                    pH_z[j+5] = H_z;  
	                   
	                }                
	      }
	      
	      
	         
	      
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  static inline
	           void H_XYZ_VP_xmm4c4_unroll2x( const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __ATTR_ALIGN__(16) xmm4c4_t k;
	                __ATTR_ALIGN__(16) xmm4c4_t H_x;
	                __ATTR_ALIGN__(16) xmm4c4_t H_y;
	                __ATTR_ALIGN__(16) xmm4c4_t H_z;
	                register __m128 vpolx;
	                register __m128 vpoly;
	                register __m128 vpolz;
	                register __m128 vdirx;
	                register __m128 vdiry;
	                register __m128 vdirz;
	                register __m128 vrx;
	                register __m128 vry;
	                register __m128 vrz;
	                int32_t j,m,m1;
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       vpolx = pvpolx[j];
	                       vpoly = pvpoly[j];
	                       vpolz = pvpolz[j];
	                       vdirx = pvdirx[j];
	                       vdiry = pvdiry[j];
	                       vdirz = pvdirz[j];
	                       vrx   = pvrx[j];
	                       vry   = pvry[j];
	                       vrz   = pvrz[j];
	                       k     = pk[j];
	                       H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                        vdirx,vdiry,vdirz,
	                                        vrx,vry,vrz,
	                                        H_x,H_y,H_z);
	                       pH_x[j] = H_x;
	                       pH_y[j] = H_y;
	                       pH_z[j] = H_z;   
	                   }
	                   if(n<2) {return;}
	                }                    
	                
	                m1 = m+1;
	                for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           
                            vpolx = pvpolx[j+0];
	                    vpoly = pvpoly[j+0];
	                    vpolz = pvpolz[j+0];
	                    vdirx = pvdirx[j+0];
	                    vdiry = pvdiry[j+0];
	                    vdirz = pvdirz[j+0];
	                    vrx   = pvrx[j+0];
	                    vry   = pvry[j+0];
	                    vrz   = pvrz[j+0];
	                    k     = pk[j+0];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+0] = H_x;
	                    pH_y[j+0] = H_y;
	                    pH_z[j+0] = H_z;   
                            vpolx = pvpolx[j+1];
	                    vpoly = pvpoly[j+1];
	                    vpolz = pvpolz[j+1];
	                    vdirx = pvdirx[j+1];
	                    vdiry = pvdiry[j+1];
	                    vdirz = pvdirz[j+1];
	                    vrx   = pvrx[j+1];
	                    vry   = pvry[j+1];
	                    vrz   = pvrz[j+1];
	                    k     = pk[j+1];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j+1] = H_x;
	                    pH_y[j+1] = H_y;
	                    pH_z[j+1] = H_z;
	                   
	                }                
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void H_XYZ_VP_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN_(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 1;
	                __ATTR_ALIGN__(16) xmm4c4_t k;
	                __ATTR_ALIGN__(16) xmm4c4_t H_x;
	                __ATTR_ALIGN__(16) xmm4c4_t H_y;
	                __ATTR_ALIGN__(16) xmm4c4_t H_z;
	                register __m128 vpolx;
	                register __m128 vpoly;
	                register __m128 vpolz;
	                register __m128 vdirx;
	                register __m128 vdiry;
	                register __m128 vdirz;
	                register __m128 vrx;
	                register __m128 vry;
	                register __m128 vrz;
	                int32_t j;
	                
	              	for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           
                            vpolx = pvpolx[j];
	                    vpoly = pvpoly[j];
	                    vpolz = pvpolz[j];
	                    vdirx = pvdirx[j];
	                    vdiry = pvdiry[j];
	                    vdirz = pvdirz[j];
	                    vrx   = pvrx[j];
	                    vry   = pvry[j];
	                    vrz   = pvrz[j];
	                    k     = pk[j];
	                    H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                     vdirx,vdiry,vdirz,
	                                     vrx,vry,vrz,
	                                     H_x,H_y,H_z);
	                    pH_x[j] = H_x;
	                    pH_y[j] = H_y;
	                    pH_z[j] = H_z;   
                         	                   
	                }                
	      }
	      
	      
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) vpolx,
	                                 const float * __restrict __ATTR_ALIGN__(16) vpoly,
	                                 const float * __restrict __ATTR_ALIGN__(16) vpolz,
	                                 const float * __restrict __ATTR_ALIGN__(16) vdirx,
	                                 const float * __restrict __ATTR_ALIGN__(16) vdiry,
	                                 const float * __restrict __ATTR_ALIGN__(16) vdirz,
	                                 const float * __restrict __ATTR_ALIGN__(16) vrx,
	                                 const float * __restrict __ATTR_ALIGN__(16) vry,
	                                 const float * __restrict __ATTR_ALIGN__(16) vrz,
	                                 const xmm4c4_t k,
	                                 xmm4c4_t & H_x,
	                                 xmm4c4_t & H_y,
	                                 xmm4c4_t & H_z) {
	               
	                register __m128 vpolx = _mm_load_ps(&vpolx[0]);
	                register __m128 vpoly = _mm_load_ps(&vpoly[0]);
	                register __m128 vpolz = _mm_load_ps(&vpolz[0]);
	                register __m128 vdirx = _mm_load_ps(&vdirx[0]);
	                register __m128 vdiry = _mm_load_ps(&vdiry[0]);
	                register __m128 vdirz = _mm_load_ps(&vdirz[0]);
	                register __m128 vrx   = _mm_load_ps(&vrx[0]);
	                register __m128 vry   = _mm_load_ps(&vry[0]);
	                register __m128 vrz   = _mm_load_ps(&vrz[0]);
	               	__m128 dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_xmm4c4(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm_set1_ps(1.0f);
	                ir = _mm_setzero_ps();
	                cmul_xmm4c4(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm_mul_ps(dp,cer);
	                cei = _mm_mul_ps(dp,cei);
	                cexp_xmm4c4(cer,cei,&expr,&expi);
	                H_x.re = _mm_mul_ps(vpolx,expr);
	                H_x.im = _mm_mul_ps(vpolx,expi);
	                H_y.re = _mm_mul_ps(vpoly,expr);
	                H_y.im = _mm_mul_ps(vpoly,expi);
	                H_z.re = _mm_mul_ps(vpolz,expr);
	                H_z.im = _mm_mul_ps(vpolz,expi);
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_xmm4c4_u(const float * __restrict  vpolx,
	                                 const float * __restrict  vpoly,
	                                 const float * __restrict  vpolz,
	                                 const float * __restrict  vdirx,
	                                 const float * __restrict  vdiry,
	                                 const float * __restrict  vdirz,
	                                 const float * __restrict  vrx,
	                                 const float * __restrict  vry,
	                                 const float * __restrict  vrz,
	                                 const xmm4c4_t k,
	                                 xmm4c4_t & H_x,
	                                 xmm4c4_t & H_y,
	                                 xmm4c4_t & H_z) {
	               
	                register __m128 vpolx = _mm_loadu_ps(&vpolx[0]);
	                register __m128 vpoly = _mm_loadu_ps(&vpoly[0]);
	                register __m128 vpolz = _mm_loadu_ps(&vpolz[0]);
	                register __m128 vdirx = _mm_loadu_ps(&vdirx[0]);
	                register __m128 vdiry = _mm_loadu_ps(&vdiry[0]);
	                register __m128 vdirz = _mm_loadu_ps(&vdirz[0]);
	                register __m128 vrx   = _mm_loadu_ps(&vrx[0]);
	                register __m128 vry   = _mm_loadu_ps(&vry[0]);
	                register __m128 vrz   = _mm_loadu_ps(&vrz[0]);
	               	__m128 dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_xmm4c4(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm_set1_ps(1.0);
	                ir = _mm_setzero_ps();
	                cmul_xmm4c4(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm_mul_ps(dp,cer);
	                cei = _mm_mul_ps(dp,cei);
	                cexp_xmm4c4(cer,cei,&expr,&expi);
	                H_x.re = _mm_mul_ps(vpolx,expr);
	                H_x.im = _mm_mul_ps(vpolx,expi);
	                H_y.re = _mm_mul_ps(vpoly,expr);
	                H_y.im = _mm_mul_ps(vpoly,expi);
	                H_z.re = _mm_mul_ps(vpolz,expr);
	                H_z.im = _mm_mul_ps(vpolz,expi);
	        }
	        
	        
	        /*
	             
     ! Magnetic Field (SIMD data-types) [plane-wave], polarization 'vpol' of
     !  wave-vector argument:  vdir*k at sixteen points 'r'.
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_xmm4c4(const __m128 vpolx,
	                                 const __m128 vpoly,
	                                 const __m128 vpolz,
	                                 const __m128 vdirx,
	                                 const __m128 vdiry,
	                                 const __m128 vdirz,
	                                 const xmm4c4_t k,
	                                 const __m128 omega,
	                                 const __m128 vrx,
	                                 const __m128 vry,
	                                 const __m128 vrz,
	                                 xmm4c4_t & B_x,
	                                 xmm4c4_t & B_y,
	                                 xmm4c4_t & B_z) {
	                                 
	                const __m128 mu0 = _mm_set1_ps(0.0000012566370614359173);
	                xmm4c4_t cdirx;
	                xmm4c4_t cdiry;
	                xmm4c4_t cdirz;
	                xmm4c4_t H_x;
	                xmm4c4_t H_y;
	                xmm4c4_t H_z;
	                xmm4c4_t cpx;
	                xmm4c4_t cpy;
	                xmm4c4_t cpz;
	                xmm4c4_t t0;
	                __m128 zz0;
	                H_XYZ_VP_xmm4c4(vpolx,vpoy,vpolz,
	                               	 vdirx,vdiry,vdirz,
	                                 vrx,vry,vrz,
	                                 H_x,H_y,H_z);
	                                 	                                 
	                cdirx.re = vdirx;
	                cdirx.im = _mm_setzero_ps();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm_mul_ps(omega,mu0);
	                t0.re    = _mm_div_ps(k.re,zz0);
	                t0.im    = _mm_div_ps(k.im,zz0);
	                
	                scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                                H_x,H_y,H_z,
                                        cpx,cpy,cpz);                
	                     	                                
	                cmul_xmm4c4(t0.re,t0.im,
	                             cpx.re,cpx.im,
	                             &B_x.re,&B_x.im);
	                             
	                cmul_xmm4c4(t0.re,t0.im,
	                             cpy.re,cpy.im,
	                             &B_y.re,&B_y.im);
	                            
	                cmul_xmm4c4(t0.re,t0.im,
	                             cpz.re,cpz.im,
	                             &B_z.re,&B_z.im);
	                          
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  static inline
	           void B_XYZ_VP_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                  __ATTR_ALIGN__(16) xmm4c4_t k;
	                  __ATTR_ALIGN__(16) xmm4c4_t B_x;
	                  __ATTR_ALIGN__(16) xmm4c4_t B_y;
	                  __ATTR_ALIGN__(16) xmm4c4_t B_z;
	                 register __m128 vpolx;
	                 register __m128 vpoly;
	                 register __m128 vpolz;
	                 register __m128 vdirx;
	                 register __m128 vdiry;
	                 register __m128 vdirz;
	                 register __m128 vrx;
	                 register __m128 vry;
	                 register __m128 vrz;
	                 register __m128 omg;
	                 int32_t j,m,m1;   
	                 
	                 m = n%10;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                        vpolx = pvpolx[j];
	                        vpoly = pvpoly[j];
	                        vpolz = pvpolz[j];
	                        vdirx = pvdirx[j];
	                        vdiry = pvdiry[j];
	                        vdirz = pvdirz[j];
	                        vrx   = pvrx[j];
	                        vry   = pvry[j];
	                        vrz   = pvrz[j];
	                        omg   = pomega[j];
	                        k     = pk[j];
	                        B_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                         vdirx,vdiry,vdirz,
	                                         k,omg,vrx,vry,vrz,
	                                         B_x,B_y,B_z);
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                    }
	                    if(n<10) {return;}
	                 }
	                 
	                 m1 = m+1;
	                 for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T1); 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T2);  
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_NTA); 	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           	
                                vpolx = pvpolx[j+0];
	                        vpoly = pvpoly[j+0];
	                        vpolz = pvpolz[j+0];
	                        vdirx = pvdirx[j+0];
	                        vdiry = pvdiry[j+0];
	                        vdirz = pvdirz[j+0];
	                        vrx   = pvrx[j+0];
	                        vry   = pvry[j+0];
	                        vrz   = pvrz[j+0];
	                        omg   = pomega[j+0];
	                        k     = pk[j+0];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;    
	                        vpolx = pvpolx[j+1];
	                        vpoly = pvpoly[j+1];
	                        vpolz = pvpolz[j+1];
	                        vdirx = pvdirx[j+1];
	                        vdiry = pvdiry[j+1];
	                        vdirz = pvdirz[j+1];
	                        vrx   = pvrx[j+1];
	                        vry   = pvry[j+1];
	                        vrz   = pvrz[j+1];
	                        omg   = pomega[j+1];
	                        k     = pk[j+1];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+1] = B_x;
	                        pB_y[j+1] = B_y;
	                        pB_z[j+1] = B_z;
	                        vpolx = pvpolx[j+2];
	                        vpoly = pvpoly[j+2];
	                        vpolz = pvpolz[j+2];
	                        vdirx = pvdirx[j+2];
	                        vdiry = pvdiry[j+2];
	                        vdirz = pvdirz[j+2];
	                        vrx   = pvrx[j+2];
	                        vry   = pvry[j+2];
	                        vrz   = pvrz[j+2];
	                        omg   = pomega[j+2];
	                        k     = pk[j+2];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+2] = B_x;
	                        pB_y[j+2] = B_y;
	                        pB_z[j+2] = B_z;
	                        vpolx = pvpolx[j+3];
	                        vpoly = pvpoly[j+3];
	                        vpolz = pvpolz[j+3];
	                        vdirx = pvdirx[j+3];
	                        vdiry = pvdiry[j+3];
	                        vdirz = pvdirz[j+3];
	                        vrx   = pvrx[j+3];
	                        vry   = pvry[j+3];
	                        vrz   = pvrz[j+3];
	                        omg   = pomega[j+3];
	                        k     = pk[j+3];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+3] = B_x;
	                        pB_y[j+3] = B_y;
	                        pB_z[j+3] = B_z;
	                        vpolx = pvpolx[j+4];
	                        vpoly = pvpoly[j+4];
	                        vpolz = pvpolz[j+4];
	                        vdirx = pvdirx[j+4];
	                        vdiry = pvdiry[j+4];
	                        vdirz = pvdirz[j+4];
	                        vrx   = pvrx[j+4];
	                        vry   = pvry[j+4];
	                        vrz   = pvrz[j+4];
	                        omg   = pomega[j+4];
	                        k     = pk[j+4];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+4] = B_x;
	                        pB_y[j+4] = B_y;
	                        pB_z[j+4] = B_z;
	                        vpolx = pvpolx[j+5];
	                        vpoly = pvpoly[j+5];
	                        vpolz = pvpolz[j+5];
	                        vdirx = pvdirx[j+5];
	                        vdiry = pvdiry[j+5];
	                        vdirz = pvdirz[j+5];
	                        vrx   = pvrx[j+5];
	                        vry   = pvry[j+5];
	                        vrz   = pvrz[j+5];
	                        omg   = pomega[j+5];
	                        k     = pk[j+5];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+5] = B_x;
	                        pB_y[j+5] = B_y;
	                        pB_z[j+5] = B_z;
	                        vpolx = pvpolx[j+6];
	                        vpoly = pvpoly[j+6];
	                        vpolz = pvpolz[j+6];
	                        vdirx = pvdirx[j+6];
	                        vdiry = pvdiry[j+6];
	                        vdirz = pvdirz[j+6];
	                        vrx   = pvrx[j+6];
	                        vry   = pvry[j+6];
	                        vrz   = pvrz[j+6];
	                        omg   = pomega[j+6];
	                        k     = pk[j+6];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+6] = B_x;
	                        pB_y[j+6] = B_y;
	                        pB_z[j+6] = B_z;  
	                        vpolx = pvpolx[j+7];
	                        vpoly = pvpoly[j+7];
	                        vpolz = pvpolz[j+7];
	                        vdirx = pvdirx[j+7];
	                        vdiry = pvdiry[j+7];
	                        vdirz = pvdirz[j+7];
	                        vrx   = pvrx[j+7];
	                        vry   = pvry[j+7];
	                        vrz   = pvrz[j+7];
	                        omg   = pomega[j+7];
	                        k     = pk[j+7];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+7] = B_x;
	                        pB_y[j+7] = B_y;
	                        pB_z[j+7] = B_z;   
	                        vpolx = pvpolx[j+8];
	                        vpoly = pvpoly[j+8];
	                        vpolz = pvpolz[j+8];
	                        vdirx = pvdirx[j+8];
	                        vdiry = pvdiry[j+8];
	                        vdirz = pvdirz[j+8];
	                        vrx   = pvrx[j+8];
	                        vry   = pvry[j+8];
	                        vrz   = pvrz[j+8];
	                        omg   = pomega[j+8];
	                        k     = pk[j+8];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+8] = B_x;
	                        pB_y[j+8] = B_y;
	                        pB_z[j+8] = B_z; 
	                        vpolx = pvpolx[j+9];
	                        vpoly = pvpoly[j+9];
	                        vpolz = pvpolz[j+9];
	                        vdirx = pvdirx[j+9];
	                        vdiry = pvdiry[j+9];
	                        vdirz = pvdirz[j+9];
	                        vrx   = pvrx[j+9];
	                        vry   = pvry[j+9];
	                        vrz   = pvrz[j+9];
	                        omg   = pomega[j+9];
	                        k     = pk[j+9];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+9] = B_x;
	                        pB_y[j+9] = B_y;
	                        pB_z[j+9] = B_z;                                          
	                 }
	                                               
	         }
	         
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_VP_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                 xmm4c4_t k;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 vpolx;
	                 register __m128 vpoly;
	                 register __m128 vpolz;
	                 register __m128 vdirx;
	                 register __m128 vdiry;
	                 register __m128 vdirz;
	                 register __m128 vrx;
	                 register __m128 vry;
	                 register __m128 vrz;
	                 register __m128 omg;
	                 int32_t j,m,m1;   
	                 
	                 m = n%6;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                        vpolx = pvpolx[j];
	                        vpoly = pvpoly[j];
	                        vpolz = pvpolz[j];
	                        vdirx = pvdirx[j];
	                        vdiry = pvdiry[j];
	                        vdirz = pvdirz[j];
	                        vrx   = pvrx[j];
	                        vry   = pvry[j];
	                        vrz   = pvrz[j];
	                        omg   = pomega[j];
	                        k     = pk[j];
	                        B_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                         vdirx,vdiry,vdirz,
	                                         k,omg,vrx,vry,vrz,
	                                         B_x,B_y,B_z);
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                    }
	                    if(n<6) {return;}
	                 }
	                 
	                 m1 = m+1;
	                 for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T1); 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T2);  
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_NTA); 	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           	
                                vpolx = pvpolx[j+0];
	                        vpoly = pvpoly[j+0];
	                        vpolz = pvpolz[j+0];
	                        vdirx = pvdirx[j+0];
	                        vdiry = pvdiry[j+0];
	                        vdirz = pvdirz[j+0];
	                        vrx   = pvrx[j+0];
	                        vry   = pvry[j+0];
	                        vrz   = pvrz[j+0];
	                        omg   = pomega[j+0];
	                        k     = pk[j+0];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;    
	                        vpolx = pvpolx[j+1];
	                        vpoly = pvpoly[j+1];
	                        vpolz = pvpolz[j+1];
	                        vdirx = pvdirx[j+1];
	                        vdiry = pvdiry[j+1];
	                        vdirz = pvdirz[j+1];
	                        vrx   = pvrx[j+1];
	                        vry   = pvry[j+1];
	                        vrz   = pvrz[j+1];
	                        omg   = pomega[j+1];
	                        k     = pk[j+1];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+1] = B_x;
	                        pB_y[j+1] = B_y;
	                        pB_z[j+1] = B_z;
	                        vpolx = pvpolx[j+2];
	                        vpoly = pvpoly[j+2];
	                        vpolz = pvpolz[j+2];
	                        vdirx = pvdirx[j+2];
	                        vdiry = pvdiry[j+2];
	                        vdirz = pvdirz[j+2];
	                        vrx   = pvrx[j+2];
	                        vry   = pvry[j+2];
	                        vrz   = pvrz[j+2];
	                        omg   = pomega[j+2];
	                        k     = pk[j+2];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+2] = B_x;
	                        pB_y[j+2] = B_y;
	                        pB_z[j+2] = B_z;
	                        vpolx = pvpolx[j+3];
	                        vpoly = pvpoly[j+3];
	                        vpolz = pvpolz[j+3];
	                        vdirx = pvdirx[j+3];
	                        vdiry = pvdiry[j+3];
	                        vdirz = pvdirz[j+3];
	                        vrx   = pvrx[j+3];
	                        vry   = pvry[j+3];
	                        vrz   = pvrz[j+3];
	                        omg   = pomega[j+3];
	                        k     = pk[j+3];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+3] = B_x;
	                        pB_y[j+3] = B_y;
	                        pB_z[j+3] = B_z;
	                        vpolx = pvpolx[j+4];
	                        vpoly = pvpoly[j+4];
	                        vpolz = pvpolz[j+4];
	                        vdirx = pvdirx[j+4];
	                        vdiry = pvdiry[j+4];
	                        vdirz = pvdirz[j+4];
	                        vrx   = pvrx[j+4];
	                        vry   = pvry[j+4];
	                        vrz   = pvrz[j+4];
	                        omg   = pomega[j+4];
	                        k     = pk[j+4];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+4] = B_x;
	                        pB_y[j+4] = B_y;
	                        pB_z[j+4] = B_z;
	                        vpolx = pvpolx[j+5];
	                        vpoly = pvpoly[j+5];
	                        vpolz = pvpolz[j+5];
	                        vdirx = pvdirx[j+5];
	                        vdiry = pvdiry[j+5];
	                        vdirz = pvdirz[j+5];
	                        vrx   = pvrx[j+5];
	                        vry   = pvry[j+5];
	                        vrz   = pvrz[j+5];
	                        omg   = pomega[j+5];
	                        k     = pk[j+5];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+5] = B_x;
	                        pB_y[j+5] = B_y;
	                        pB_z[j+5] = B_z;
	                                                        
	                 }
	                                               
	         }
	         
	         
	        
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_VP_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                 xmm4c4_t k;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 vpolx;
	                 register __m128 vpoly;
	                 register __m128 vpolz;
	                 register __m128 vdirx;
	                 register __m128 vdiry;
	                 register __m128 vdirz;
	                 register __m128 vrx;
	                 register __m128 vry;
	                 register __m128 vrz;
	                 register __m128 omg;
	                 int32_t j,m,m1;   
	                 
	                 m = n%2;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                        vpolx = pvpolx[j];
	                        vpoly = pvpoly[j];
	                        vpolz = pvpolz[j];
	                        vdirx = pvdirx[j];
	                        vdiry = pvdiry[j];
	                        vdirz = pvdirz[j];
	                        vrx   = pvrx[j];
	                        vry   = pvry[j];
	                        vrz   = pvrz[j];
	                        omg   = pomega[j];
	                        k     = pk[j];
	                        B_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                         vdirx,vdiry,vdirz,
	                                         k,omg,vrx,vry,vrz,
	                                         B_x,B_y,B_z);
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                    }
	                    if(n<2) {return;}
	                 }
	                 
	                 m1 = m+1;
	                 for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T1); 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T2);  
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_NTA); 	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           	
                                vpolx = pvpolx[j+0];
	                        vpoly = pvpoly[j+0];
	                        vpolz = pvpolz[j+0];
	                        vdirx = pvdirx[j+0];
	                        vdiry = pvdiry[j+0];
	                        vdirz = pvdirz[j+0];
	                        vrx   = pvrx[j+0];
	                        vry   = pvry[j+0];
	                        vrz   = pvrz[j+0];
	                        omg   = pomega[j+0];
	                        k     = pk[j+0];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;    
	                        vpolx = pvpolx[j+1];
	                        vpoly = pvpoly[j+1];
	                        vpolz = pvpolz[j+1];
	                        vdirx = pvdirx[j+1];
	                        vdiry = pvdiry[j+1];
	                        vdirz = pvdirz[j+1];
	                        vrx   = pvrx[j+1];
	                        vry   = pvry[j+1];
	                        vrz   = pvrz[j+1];
	                        omg   = pomega[j+1];
	                        k     = pk[j+1];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+1] = B_x;
	                        pB_y[j+1] = B_y;
	                        pB_z[j+1] = B_z;
	                    	                                                        
	                 }
	                                               
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_VP_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128 * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 1;
	                 xmm4c4_t k;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 vpolx;
	                 register __m128 vpoly;
	                 register __m128 vpolz;
	                 register __m128 vdirx;
	                 register __m128 vdiry;
	                 register __m128 vdirz;
	                 register __m128 vrx;
	                 register __m128 vry;
	                 register __m128 vrz;
	                 register __m128 omg;
	                 int32_t j;   
	                 
	                 for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T0);	 
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T1);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T1); 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T1);          
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_T2);	
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_T2);  
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_T2);  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&pvpolx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpoly[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvpolz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdiry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvdirz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvry[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pvrz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomega[j+PF_DIST],_MM_HINT_NTA); 	 
	                    _mm_prefetch((char*)&pk[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pk[j+PF_DIST].im,_MM_HINT_NTA);    
#endif	     	           	
                                vpolx = pvpolx[j+0];
	                        vpoly = pvpoly[j+0];
	                        vpolz = pvpolz[j+0];
	                        vdirx = pvdirx[j+0];
	                        vdiry = pvdiry[j+0];
	                        vdirz = pvdirz[j+0];
	                        vrx   = pvrx[j+0];
	                        vry   = pvry[j+0];
	                        vrz   = pvrz[j+0];
	                        omg   = pomega[j+0];
	                        k     = pk[j+0];
	                        B_XYZ_xmm4c4(vpolx,vpoly,vpolz,
	                                      vdirx,vdiry,vdirz,
	                                      k,omg,vrx,vry,vrz,
	                                      B_x,B_y,B_z);
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;    
	                       	                    	                                                        
	                 }
	                                               
	         }
	         
	         
	                  
	         
	         
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_VP_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                   const float * __restrict __ATTR_ALIGN__(16) pomega,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvrx,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvry,
	                                   const float * __restrict __ATTR_ALIGN__(16) pvrz,
	                                   const xmm4c4_t k,
	                                   xmm4c4_t & B_x,
	                                   xmm4c4_t & B_y,
	                                   xmm4c4_t & B_z) {
	                         
	                register __m128 vpolx = _mm_load_ps(&pvpolx[0]);
	                register __m128 vpoly = _mm_load_ps(&pvpoly[0]);    
	                register __m128 vpolz = _mm_load_ps(&pvpolz[0]);  
	                register __m128 vdirx = _mm_load_ps(&pvdirx[0]);  
	                register __m128 vdiry = _mm_load_ps(&pvdiry[0]);
	                register __m128 vdirz = _mm_load_ps(&pvdirz[0]); 
	                register __m128 onega = _mm_load_ps(&pomega[0]);
	                register __m128 vrx   = _mm_load_ps(&pvrx[0]);
	                register __m128 vry   = _mm_load_ps(&pvry[0]);
	                register __m128 vrz   = _mm_load_ps(&pvrz[0]);        
	                const __m128 mu0 = _mm_set1_ps(0.0000012566370614359173f);
	                xmm4c4_t cdirx;
	                xmm4c4_t cdiry;
	                xmm4c4_t cdirz;
	                xmm4c4_t H_x;
	                xmm4c4_t H_y;
	                xmm4c4_t H_z;
	                xmm4c4_t cpx;
	                xmm4c4_t cpy;
	                xmm4c4_t cpz;
	                xmm4c4_t t0;
	                __m128 zz0;
	                H_XYZ_VP_xmm4c4(vpolx,
	                                 vpoly,
	                                 vpolz,
	                                 vdirx,
	                                 vdiry,
	                                 vdirz,
	                                 vrx,
	                                 vry,
	                                 vrz,
	                                 H_x,
	                                 H_y,
	                                 H_z);
	                                 
	                cdirx.re = vdirx;
	                cdirx.im = _mm_setzero_ps();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm_mul_ps(omega,mu0);
	                t0.re    = _mm_div_ps(k.re,zz0);
	                t0.im    = _mm_div_ps(k.im,zz0);
	                scrossc_xmm4c4(cdirx,
	                                cdiry,
	                                cdirz,
	                                H_x,
	                                H_y,
	                                H_z,
	                                cpx,
	                                cpy,
	                                cpz);
	                                
	                cmul_xmm4c4(t0.re,
	                             t0.im,
	                             cpx.re,
	                             cpx.im,
	                             &B_x.re,
	                             &B_x.im);
	                cmul_xmm4c4(t0.re,
	                             t0.im,
	                             cpy.re,
	                             cpy.im,
	                             &B_y.re,
	                             &B_y.im);
	                cmul_xmm4c4(t0.re,
	                             t0.im,
	                             cpz.re,
	                             cpz.im,
	                             &B_z.re,
	                             &B_z.im);
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_xmm4c4_u( const float * __restrict  pvpolx,
	                                   const float * __restrict  pvpoly,
	                                   const float * __restrict  pvpolz,
	                                   const float * __restrict  pvdirx,
	                                   const float * __restrict  pvdiry,
	                                   const float * __restrict  pvdirz,
	                                   const float * __restrict  pomega,
	                                   const float * __restrict  pvrx,
	                                   const float * __restrict  pvry,
	                                   const float * __restrict  pvrz,
	                                   const xmm4c4_t k,
	                                   xmm4c4_t & B_x,
	                                   xmm4c4_t & B_y,
	                                   xmm4c4_t & B_z) {
	                         
	                register __m128 vpolx = _mm_loadu_ps(&pvpolx[0]);
	                register __m128 vpoly = _mm_loadu_ps(&pvpoly[0]);    
	                register __m128 vpolz = _mm_loadu_ps(&pvpolz[0]);  
	                register __m128 vdirx = _mm_loadu_ps(&pvdirx[0]);  
	                register __m128 vdiry = _mm_loadu_ps(&pvdiry[0]);
	                register __m128 vdirz = _mm_loadu_ps(&pvdirz[0]); 
	                register __m128 onega = _mm_loadu_ps(&pomega[0]);
	                register __m128 vrx   = _mm_loadu_ps(&pvrx[0]);
	                register __m128 vry   = _mm_loadu_ps(&pvry[0]);
	                register __m128 vrz   = _mm_loadu_ps(&pvrz[0]);        
	                const __m128 mu0 = _mm_set1_ps(0.0000012566370614359173f);
	                xmm4c4_t cdirx;
	                xmm4c4_t cdiry;
	                xmm4c4_t cdirz;
	                xmm4c4_t H_x;
	                xmm4c4_t H_y;
	                xmm4c4_t H_z;
	                xmm4c4_t cpx;
	                xmm4c4_t cpy;
	                xmm4c4_t cpz;
	                xmm4c4_t t0;
	                __m128 zz0;
	                H_XYZ_VP_xmm4c4(vpolx,
	                                 vpoly,
	                                 vpolz,
	                                 vdirx,
	                                 vdiry,
	                                 vdirz,
	                                 vrx,
	                                 vry,
	                                 vrz,
	                                 H_x,
	                                 H_y,
	                                 H_z);
	                                 
	                cdirx.re = vdirx;
	                cdirx.im = _mm_setzero_ps();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm_mul_ps(omega,mu0);
	                t0.re    = _mm_div_ps(k.re,zz0);
	                t0.im    = _mm_div_ps(k.im,zz0);
	                scrossc_xmm4c4(cdirx,
	                                cdiry,
	                                cdirz,
	                                H_x,
	                                H_y,
	                                H_z,
	                                cpx,
	                                cpy,
	                                cpz);
	                                
	                cmul_xmm4c4(t0.re,
	                             t0.im,
	                             cpx.re,
	                             cpx.im,
	                             &B_x.re,
	                             &B_x.im);
	                cmul_xmm4c4(t0.re,
	                             t0.im,
	                             cpy.re,
	                             cpy.im,
	                             &B_y.re,
	                             &B_y.im);
	                cmul_xmm4c4(t0.re,
	                             t0.im,
	                             cpz.re,
	                             cpz.im,
	                             &B_z.re,
	                             &B_z.im);
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_xmm4c4(const __m128 tht,
	                                      const __m128 phi,
	                                      const __m128 psi,
	                                      const __m128 omg,
	                                      const __m128 px,
	                                      const __m128 py,
	                                      const __m128 pz,
	                                      const xmm4c4_t r,
	                                      xmm4c4_t & H_x,
	                                      xmm4c4_t & H_y,
	                                      xmm4c4_t & H_z,
	                                      xmm4c4_t & B_x,
	                                      xmm4c4_t & B_y,
	                                      xmm4c4_t & B_z) {
	                                      
	                
	                const __m128 c = _mm_set1_ps(299792458.0f);
	                xmm4c4_t k;
	                register __m128 vpolx,vpoly,vpolz;
	                register __m128 vdirx,vdiry,vdirz;
	                register __m128 t0;
	                
	                t0 = _mm_div_ps(omg,c);
	                dir_vec_xmm4c4(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm_mul_ps(r.re,c);
	                k.im = _mm_mul_ps(r.im,c);
	                pol_vec_xmm4c4(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_P_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;    
	                 xmm4c4_t r;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 psi;
	                 register __m128 omg;
	                 register __m128 px;
	                 register __m128 py;
	                 register __m128 pz;
	                 int32_t j,m,m1;
	                 
	                 m = n%10;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                         tht = ptht[j];
	                         phi = pphi[j];
	                         psi = ppsi[j];
	                         omg = pomg[j];
	                         px  = ppx[j];
	                         py  = ppy[j];
	                         pz  = ppz[j];
	                         r   = pr[j];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j] = H_x;
	                         pH_y[j] = H_y;
	                         pH_z[j] = H_z;
	                         pB_x[j] = B_x;
	                         pB_y[j] = B_y;
	                         pB_z[j] = B_z;
	                      }
	                      if(n<10) {return;}  
	                   }       
	                   
	                   m1 = m+1;
	                   for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T1);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T2);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_NTA);       
#endif	     	       	            
                                 tht = ptht[j+0];
	                         phi = pphi[j+0];
	                         psi = ppsi[j+0];
	                         omg = pomg[j+0];
	                         px  = ppx[j+0];
	                         py  = ppy[j+0];
	                         pz  = ppz[j+0];
	                         r   = pr[j+0];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+0] = H_x;
	                         pH_y[j+0] = H_y;
	                         pH_z[j+0] = H_z;
	                         pB_x[j+0] = B_x;
	                         pB_y[j+0] = B_y;
	                         pB_z[j+0] = B_z;  
	                         tht = ptht[j+1];
	                         phi = pphi[j+1];
	                         psi = ppsi[j+1];
	                         omg = pomg[j+1];
	                         px  = ppx[j+1];
	                         py  = ppy[j+1];
	                         pz  = ppz[j+1];
	                         r   = pr[j+1];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+1] = H_x;
	                         pH_y[j+1] = H_y;
	                         pH_z[j+1] = H_z;
	                         pB_x[j+1] = B_x;
	                         pB_y[j+1] = B_y;
	                         pB_z[j+1] = B_z;  
	                         tht = ptht[j+2];
	                         phi = pphi[j+2];
	                         psi = ppsi[j+2];
	                         omg = pomg[j+2];
	                         px  = ppx[j+2];
	                         py  = ppy[j+2];
	                         pz  = ppz[j+2];
	                         r   = pr[j+2];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+2] = H_x;
	                         pH_y[j+2] = H_y;
	                         pH_z[j+2] = H_z;
	                         pB_x[j+2] = B_x;
	                         pB_y[j+2] = B_y;
	                         pB_z[j+2] = B_z;  
	                         tht = ptht[j+3];
	                         phi = pphi[j+3];
	                         psi = ppsi[j+3];
	                         omg = pomg[j+3];
	                         px  = ppx[j+3];
	                         py  = ppy[j+3];
	                         pz  = ppz[j+3];
	                         r   = pr[j+3];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+3] = H_x;
	                         pH_y[j+3] = H_y;
	                         pH_z[j+3] = H_z;
	                         pB_x[j+3] = B_x;
	                         pB_y[j+3] = B_y;
	                         pB_z[j+3] = B_z;  
	                         tht = ptht[j+4];
	                         phi = pphi[j+4];
	                         psi = ppsi[j+4];
	                         omg = pomg[j+4];
	                         px  = ppx[j+4];
	                         py  = ppy[j+4];
	                         pz  = ppz[j+4];
	                         r   = pr[j+4];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+4] = H_x;
	                         pH_y[j+4] = H_y;
	                         pH_z[j+4] = H_z;
	                         pB_x[j+4] = B_x;
	                         pB_y[j+4] = B_y;
	                         pB_z[j+4] = B_z;  
	                         tht = ptht[j+5];
	                         phi = pphi[j+5];
	                         psi = ppsi[j+5];
	                         omg = pomg[j+5];
	                         px  = ppx[j+5];
	                         py  = ppy[j+5];
	                         pz  = ppz[j+5];
	                         r   = pr[j+5];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+5] = H_x;
	                         pH_y[j+5] = H_y;
	                         pH_z[j+5] = H_z;
	                         pB_x[j+5] = B_x;
	                         pB_y[j+5] = B_y;
	                         pB_z[j+5] = B_z;  
	                         tht = ptht[j+6];
	                         phi = pphi[j+6];
	                         psi = ppsi[j+6];
	                         omg = pomg[j+6];
	                         px  = ppx[j+6];
	                         py  = ppy[j+6];
	                         pz  = ppz[j+6];
	                         r   = pr[j+6];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+6] = H_x;
	                         pH_y[j+6] = H_y;
	                         pH_z[j+6] = H_z;
	                         pB_x[j+6] = B_x;
	                         pB_y[j+6] = B_y;
	                         pB_z[j+6] = B_z;  
	                         tht = ptht[j+7];
	                         phi = pphi[j+7];
	                         psi = ppsi[j+7];
	                         omg = pomg[j+7];
	                         px  = ppx[j+7];
	                         py  = ppy[j+7];
	                         pz  = ppz[j+7];
	                         r   = pr[j+7];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+7] = H_x;
	                         pH_y[j+7] = H_y;
	                         pH_z[j+7] = H_z;
	                         pB_x[j+7] = B_x;
	                         pB_y[j+7] = B_y;
	                         pB_z[j+7] = B_z;  
	                         tht = ptht[j+8];
	                         phi = pphi[j+8];
	                         psi = ppsi[j+8];
	                         omg = pomg[j+8];
	                         px  = ppx[j+8];
	                         py  = ppy[j+8];
	                         pz  = ppz[j+8];
	                         r   = pr[j+8];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+8] = H_x;
	                         pH_y[j+8] = H_y;
	                         pH_z[j+8] = H_z;
	                         pB_x[j+8] = B_x;
	                         pB_y[j+8] = B_y;
	                         pB_z[j+8] = B_z;  
	                         tht = ptht[j+9];
	                         phi = pphi[j+9];
	                         psi = ppsi[j+9];
	                         omg = pomg[j+9];
	                         px  = ppx[j+9];
	                         py  = ppy[j+9];
	                         pz  = ppz[j+9];
	                         r   = pr[j+9];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+9] = H_x;
	                         pH_y[j+9] = H_y;
	                         pH_z[j+9] = H_z;
	                         pB_x[j+9] = B_x;
	                         pB_y[j+9] = B_y;
	                         pB_z[j+9] = B_z;  
	                                                                               
	                   }
	                                             
	      }
	      
	      
	          
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_P_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;    
	                 xmm4c4_t r;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 psi;
	                 register __m128 omg;
	                 register __m128 px;
	                 register __m128 py;
	                 register __m128 pz;
	                 int32_t j,m,m1;
	                 
	                 m = n%6;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                         tht = ptht[j];
	                         phi = pphi[j];
	                         psi = ppsi[j];
	                         omg = pomg[j];
	                         px  = ppx[j];
	                         py  = ppy[j];
	                         pz  = ppz[j];
	                         r   = pr[j];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j] = H_x;
	                         pH_y[j] = H_y;
	                         pH_z[j] = H_z;
	                         pB_x[j] = B_x;
	                         pB_y[j] = B_y;
	                         pB_z[j] = B_z;
	                      }
	                      if(n<6) {return;}  
	                   }       
	                   
	                   m1 = m+1;
	                   for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T1);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T2);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_NTA);       
#endif	     	       	            
                                 tht = ptht[j+0];
	                         phi = pphi[j+0];
	                         psi = ppsi[j+0];
	                         omg = pomg[j+0];
	                         px  = ppx[j+0];
	                         py  = ppy[j+0];
	                         pz  = ppz[j+0];
	                         r   = pr[j+0];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+0] = H_x;
	                         pH_y[j+0] = H_y;
	                         pH_z[j+0] = H_z;
	                         pB_x[j+0] = B_x;
	                         pB_y[j+0] = B_y;
	                         pB_z[j+0] = B_z;  
	                         tht = ptht[j+1];
	                         phi = pphi[j+1];
	                         psi = ppsi[j+1];
	                         omg = pomg[j+1];
	                         px  = ppx[j+1];
	                         py  = ppy[j+1];
	                         pz  = ppz[j+1];
	                         r   = pr[j+1];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+1] = H_x;
	                         pH_y[j+1] = H_y;
	                         pH_z[j+1] = H_z;
	                         pB_x[j+1] = B_x;
	                         pB_y[j+1] = B_y;
	                         pB_z[j+1] = B_z;  
	                         tht = ptht[j+2];
	                         phi = pphi[j+2];
	                         psi = ppsi[j+2];
	                         omg = pomg[j+2];
	                         px  = ppx[j+2];
	                         py  = ppy[j+2];
	                         pz  = ppz[j+2];
	                         r   = pr[j+2];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+2] = H_x;
	                         pH_y[j+2] = H_y;
	                         pH_z[j+2] = H_z;
	                         pB_x[j+2] = B_x;
	                         pB_y[j+2] = B_y;
	                         pB_z[j+2] = B_z;  
	                         tht = ptht[j+3];
	                         phi = pphi[j+3];
	                         psi = ppsi[j+3];
	                         omg = pomg[j+3];
	                         px  = ppx[j+3];
	                         py  = ppy[j+3];
	                         pz  = ppz[j+3];
	                         r   = pr[j+3];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+3] = H_x;
	                         pH_y[j+3] = H_y;
	                         pH_z[j+3] = H_z;
	                         pB_x[j+3] = B_x;
	                         pB_y[j+3] = B_y;
	                         pB_z[j+3] = B_z;  
	                         tht = ptht[j+4];
	                         phi = pphi[j+4];
	                         psi = ppsi[j+4];
	                         omg = pomg[j+4];
	                         px  = ppx[j+4];
	                         py  = ppy[j+4];
	                         pz  = ppz[j+4];
	                         r   = pr[j+4];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+4] = H_x;
	                         pH_y[j+4] = H_y;
	                         pH_z[j+4] = H_z;
	                         pB_x[j+4] = B_x;
	                         pB_y[j+4] = B_y;
	                         pB_z[j+4] = B_z;  
	                         tht = ptht[j+5];
	                         phi = pphi[j+5];
	                         psi = ppsi[j+5];
	                         omg = pomg[j+5];
	                         px  = ppx[j+5];
	                         py  = ppy[j+5];
	                         pz  = ppz[j+5];
	                         r   = pr[j+5];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+5] = H_x;
	                         pH_y[j+5] = H_y;
	                         pH_z[j+5] = H_z;
	                         pB_x[j+5] = B_x;
	                         pB_y[j+5] = B_y;
	                         pB_z[j+5] = B_z;  
	                                            
	                   }
	                                            
	      }
	      
	      
	        
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_P_xmm4c4_unroll2x(  const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;    
	                 xmm4c4_t r;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 psi;
	                 register __m128 omg;
	                 register __m128 px;
	                 register __m128 py;
	                 register __m128 pz;
	                 int32_t j,m,m1;
	                 
	                 m = n%2;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                         tht = ptht[j];
	                         phi = pphi[j];
	                         psi = ppsi[j];
	                         omg = pomg[j];
	                         px  = ppx[j];
	                         py  = ppy[j];
	                         pz  = ppz[j];
	                         r   = pr[j];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j] = H_x;
	                         pH_y[j] = H_y;
	                         pH_z[j] = H_z;
	                         pB_x[j] = B_x;
	                         pB_y[j] = B_y;
	                         pB_z[j] = B_z;
	                      }
	                      if(n<2) {return;}  
	                   }       
	                   
	                   m1 = m+1;
	                   for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T1);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T2);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_NTA);       
#endif	     	       	            
                                 tht = ptht[j+0];
	                         phi = pphi[j+0];
	                         psi = ppsi[j+0];
	                         omg = pomg[j+0];
	                         px  = ppx[j+0];
	                         py  = ppy[j+0];
	                         pz  = ppz[j+0];
	                         r   = pr[j+0];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+0] = H_x;
	                         pH_y[j+0] = H_y;
	                         pH_z[j+0] = H_z;
	                         pB_x[j+0] = B_x;
	                         pB_y[j+0] = B_y;
	                         pB_z[j+0] = B_z;  
	                         tht = ptht[j+1];
	                         phi = pphi[j+1];
	                         psi = ppsi[j+1];
	                         omg = pomg[j+1];
	                         px  = ppx[j+1];
	                         py  = ppy[j+1];
	                         pz  = ppz[j+1];
	                         r   = pr[j+1];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j+1] = H_x;
	                         pH_y[j+1] = H_y;
	                         pH_z[j+1] = H_z;
	                         pB_x[j+1] = B_x;
	                         pB_y[j+1] = B_y;
	                         pB_z[j+1] = B_z;  
	               }
	                                            
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_P_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128 * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm4c4_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 1;    
	                 xmm4c4_t r;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 psi;
	                 register __m128 omg;
	                 register __m128 px;
	                 register __m128 py;
	                 register __m128 pz;
	                 int32_t j;
	                 
	               	 for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T0);                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T1);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_T2);       
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppsi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pr[j+PF_DIST].im,_MM_HINT_NTA);       
#endif	     	       	            
                                 tht = ptht[j];
	                         phi = pphi[j];
	                         psi = ppsi[j];
	                         omg = pomg[j];
	                         px  = ppx[j];
	                         py  = ppy[j];
	                         pz  = ppz[j];
	                         r   = pr[j];
	                         B_XYZ_H_XYZ_P_xmm4c4(tht,phi,psi,
	                                               omg,px,py,pz,r,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                         pH_x[j] = H_x;
	                         pH_y[j] = H_y;
	                         pH_z[j] = H_z;
	                         pB_x[j] = B_x;
	                         pB_y[j] = B_y;
	                         pB_z[j] = B_z;  
	                    }
	                                            
	      }
	      
	      
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) ptht,
	                                      const float * __restrict __ATTR_ALIGN__(16) pphi,
	                                      const float * __restrict __ATTR_ALIGN__(16) ppsi,
	                                      const float * __restrict __ATTR_ALIGN__(16) pomg,
	                                      const float * __restrict __ATTR_ALIGN__(16) ppx,
	                                      const float * __restrict __ATTR_ALIGN__(16) ppy,
	                                      const float * __restrict __ATTR_ALIGN__(16) ppz,
	                                      const xmm4c4_t r,
	                                      xmm4c4_t & H_x,
	                                      xmm4c4_t & H_y,
	                                      xmm4c4_t & H_z,
	                                      xmm4c4_t & B_x,
	                                      xmm4c4_t & B_y,
	                                      xmm4c4_t & B_z) {
	                                      
	                
	                register __m128 tht = _mm_load_ps(&ptht[0]);
	                register __m128 phi = _mm_load_ps(&pphi[0]);
	                register __m128 psi = _mm_load_ps(&ppsi[0]);
	                register __m128 omg = _mm_load_ps(&pomg[0]);
	                register __m128 px  = _mm_load_ps(&ppx[0]);
	                register __m128 py  = _mm_load_ps(&ppy[0]);
	                register __m128 pz  = _mm_load_ps(&ppz[0]);
	                const __m128 c = _mm_set1_ps(299792458.0);
	                xmm4c4_t k;
	                register __m128 vpolx,vpoly,vpolz;
	                register __m128 vdirx,vdiry,vdirz;
	                register __m128 t0;
	                
	                t0 = _mm_div_ps(omg,c);
	                dir_vec_xmm4c4(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm_mul_ps(r.re,c);
	                k.im = _mm_mul_ps(r.im,c);
	                pol_vec_xmm4c4(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_xmm4c4_u(const float * __restrict ptht,
	                                      const float * __restrict  pphi,
	                                      const float * __restrict ppsi,
	                                      const float * __restrict  pomg,
	                                      const float * __restrict  ppx,
	                                      const float * __restrict  ppy,
	                                      const float * __restrict ppz,
	                                      const xmm4c4_t r,
	                                      xmm4c4_t & H_x,
	                                      xmm4c4_t & H_y,
	                                      xmm4c4_t & H_z,
	                                      xmm4c4_t & B_x,
	                                      xmm4c4_t & B_y,
	                                      xmm4c4_t & B_z) {
	                                      
	                
	                register __m128 tht = _mm_loadu_ps(&ptht[0]);
	                register __m128 phi = _mm_loadu_ps(&pphi[0]);
	                register __m128 psi = _mm_loadu_ps(&ppsi[0]);
	                register __m128 omg = _mm_loadu_ps(&pomg[0]);
	                register __m128 px  = _mm_loadu_ps(&ppx[0]);
	                register __m128 py  = _mm_loadu_ps(&ppy[0]);
	                register __m128 pz  = _mm_loadu_ps(&ppz[0]);
	                const __m128 c = _mm_set1_ps(299792458.0);
	                xmm4c4_t k;
	                register __m128 vpolx,vpoly,vpolz;
	                register __m128 vdirx,vdiry,vdirz;
	                register __m128 t0;
	                
	                t0 = _mm_div_ps(omg,c);
	                dir_vec_xmm4c4(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm_mul_ps(r.re,c);
	                k.im = _mm_mul_ps(r.im,c);
	                pol_vec_xmm4c4(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	       /*
	              ! Electric and Magnetic Fields elliptically polarized
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4(const __m128 tht,
	                                       const __m128 phi,
	                                       const __m128 omg,
	                                       const xmm4c4_t phase,
	                                       const xmm4c4_t refi,
	                                       const xmm4c4_t px,
	                                       const xmm4c4_t py,
	                                       const xmm4c4_t pz,
	                                       xmm4c4_t & H_x,
	                                       xmm4c4_t & H_y,
	                                       xmm4c4_t & H_z,
	                                       xmm4c4_t & B_x,
	                                       xmm4c4_t & B_y,
	                                       xmm4c4_t & B_z) {
	                                   
	               const __m128 c   = _mm_set1_ps(299792458.0); 
	               const __m128 mu0 = _mm_set1_ps(0.0000012566370614359173f);   
	               const __m128 psi0 = _mm_setzero_ps();
	               const __m128 C00  = _mm_setzero_ps();
	               
	               xmm4c4_t H_x_1;
	               xmm4c4_t H_y_1;
	               xmm4c4_t H_z_1;
	               xmm4c4_t H_x_2;
	               xmm4c4_t H_y_2;
	               xmm4c4_t H_z_2;
	               xmm4c4_t k;
	               xmm4c4_t t0;
	               xmm4c4_t cdirx;
	               xmm4c4_t cdiry;
	               xmm4c4_t cdirz;
	               
	               register __m128 vpolx;
	               register __m128 vpoly;
	               register __m128 vpolz;
	               register __m128 vdirx;
	               register __m128 vdiry;
	               register __m128 vdirz;
	               register __m128 cn;
	               register __m128 x0;
	               register __m128 t0r,t0i;
	               register __m128 t1r,t1i;
	               register __m128 t2r,t2i;
	               
	               dir_vec_xmm4c4(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm_div_ps(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm_mul_ps(refi.re,x0);
	               k.im     = _mm_mul_ps(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_xmm4c4(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_xmm4c4( vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm_add_ps(H_x_1.re,t0r);
	               H_x.im = _mm_add_ps(H_x_1.im,t0i);
	               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm_add_ps(H_y_1.re,t1r);
	               H_y.im = _mm_add_ps(H_y_1.im,t1i);
	               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm_add_ps(H_z_2.re,t2r);
	              H_z.im = _mm_add_ps(H_z_2.im,t2i);
	              
	              cn = cnorm_xmm4c4(H_x,H_y,H_z);
	              
	              x0     = _mm_div_ps(omg,mu0);
	              H_x.re = _mm_div_ps(H_x.re,cn);
	              H_x.im = _mm_div_ps(H_x.im,cn);
	              H_y.re = _mm_div_ps(H_y.re,cn);
	              H_y.im = _mm_div_ps(H_y.im,cn);
	              H_z.re = _mm_div_ps(H_z.re,cn);
	              H_z.im = _mm_div_ps(H_z.im,cn); 
	              
	              t0.re  = _mm_div_ps(k.re,x0);
	              t0.im  = _mm_div_ps(k.im,x0);
	              
	              scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                 xmm4c4_t phase;
	                 xmm4c4_t refi;
	                 xmm4c4_t px;
	                 xmm4c4_t py;
	                 xmm4c4_t pz;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 omg;
	                 int32_t j,m,m1;
	                 
	                 m = n%10;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                        tht   = ptht[j];
	                        phi   = pphi[j];
	                        omg   = pomg[j];
	                        phase = pphase[j];
	                        refi  = prefi[j];
	                        px    = ppx[j];
	                        py    = ppy[j];
	                        pz    = ppz[j];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j] = H_x;
	                        pH_y[j] = H_y;
	                        pH_z[j] = H_z;
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                    }
	                    if(n<10) { return;}
	                 }                     
	                 
	                 m1 = m+1;
	                 for(j = m1; j != n; j += 10) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T0);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T0);
	                                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T1);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T2);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_NTA);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_NTA);
#endif	     	   	      
                                tht   = ptht[j+0];
	                        phi   = pphi[j+0];
	                        omg   = pomg[j+0];
	                        phase = pphase[j+0];
	                        refi  = prefi[j+0];
	                        px    = ppx[j+0];
	                        py    = ppy[j+0];
	                        pz    = ppz[j+0];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+0] = H_x;
	                        pH_y[j+0] = H_y;
	                        pH_z[j+0] = H_z;
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;
	                        tht   = ptht[j+1];
	                        phi   = pphi[j+1];
	                        omg   = pomg[j+1];
	                        phase = pphase[j+1];
	                        refi  = prefi[j+1];
	                        px    = ppx[j+1];
	                        py    = ppy[j+1];
	                        pz    = ppz[j+1];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+1] = H_x;
	                        pH_y[j+1] = H_y;
	                        pH_z[j+1] = H_z;
	                        pB_x[j+1] = B_x;
	                        pB_y[j+1] = B_y;
	                        pB_z[j+1] = B_z;  
	                        tht   = ptht[j+2];
	                        phi   = pphi[j+2];
	                        omg   = pomg[j+2];
	                        phase = pphase[j+2];
	                        refi  = prefi[j+2];
	                        px    = ppx[j+2];
	                        py    = ppy[j+2];
	                        pz    = ppz[j+2];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+2] = H_x;
	                        pH_y[j+2] = H_y;
	                        pH_z[j+2] = H_z;
	                        pB_x[j+2] = B_x;
	                        pB_y[j+2] = B_y;
	                        pB_z[j+2] = B_z;
	                        tht   = ptht[j+3];
	                        phi   = pphi[j+3];
	                        omg   = pomg[j+3];
	                        phase = pphase[j+3];
	                        refi  = prefi[j+3];
	                        px    = ppx[j+3];
	                        py    = ppy[j+3];
	                        pz    = ppz[j+3];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+3] = H_x;
	                        pH_y[j+3] = H_y;
	                        pH_z[j+3] = H_z;
	                        pB_x[j+3] = B_x;
	                        pB_y[j+3] = B_y;
	                        pB_z[j+3] = B_z;  
	                        tht   = ptht[j+4];
	                        phi   = pphi[j+4];
	                        omg   = pomg[j+4];
	                        phase = pphase[j+4];
	                        refi  = prefi[j+4];
	                        px    = ppx[j+4];
	                        py    = ppy[j+4];
	                        pz    = ppz[j+4];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+4] = H_x;
	                        pH_y[j+4] = H_y;
	                        pH_z[j+4] = H_z;
	                        pB_x[j+4] = B_x;
	                        pB_y[j+4] = B_y;
	                        pB_z[j+4] = B_z; 
	                        tht   = ptht[j+5];
	                        phi   = pphi[j+5];
	                        omg   = pomg[j+5];
	                        phase = pphase[j+5];
	                        refi  = prefi[j+5];
	                        px    = ppx[j+5];
	                        py    = ppy[j+5];
	                        pz    = ppz[j+5];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+5] = H_x;
	                        pH_y[j+5] = H_y;
	                        pH_z[j+5] = H_z;
	                        pB_x[j+5] = B_x;
	                        pB_y[j+5] = B_y;
	                        pB_z[j+5] = B_z;  
	                        tht   = ptht[j+6];
	                        phi   = pphi[j+6];
	                        omg   = pomg[j+6];
	                        phase = pphase[j+6];
	                        refi  = prefi[j+6];
	                        px    = ppx[j+6];
	                        py    = ppy[j+6];
	                        pz    = ppz[j+6];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+6] = H_x;
	                        pH_y[j+6] = H_y;
	                        pH_z[j+6] = H_z;
	                        pB_x[j+6] = B_x;
	                        pB_y[j+6] = B_y;
	                        pB_z[j+6] = B_z;  
	                        tht   = ptht[j+7];
	                        phi   = pphi[j+7];
	                        omg   = pomg[j+7];
	                        phase = pphase[j+7];
	                        refi  = prefi[j+7];
	                        px    = ppx[j+7];
	                        py    = ppy[j+7];
	                        pz    = ppz[j+7];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+7] = H_x;
	                        pH_y[j+7] = H_y;
	                        pH_z[j+7] = H_z;
	                        pB_x[j+7] = B_x;
	                        pB_y[j+7] = B_y;
	                        pB_z[j+7] = B_z;  
	                        tht   = ptht[j+8];
	                        phi   = pphi[j+8];
	                        omg   = pomg[j+8];
	                        phase = pphase[j+8];
	                        refi  = prefi[j+8];
	                        px    = ppx[j+8];
	                        py    = ppy[j+8];
	                        pz    = ppz[j+8];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+8] = H_x;
	                        pH_y[j+8] = H_y;
	                        pH_z[j+8] = H_z;
	                        pB_x[j+8] = B_x;
	                        pB_y[j+8] = B_y;
	                        pB_z[j+8] = B_z; 
	                        tht   = ptht[j+9];
	                        phi   = pphi[j+9];
	                        omg   = pomg[j+9];
	                        phase = pphase[j+9];
	                        refi  = prefi[j+9];
	                        px    = ppx[j+9];
	                        py    = ppy[j+9];
	                        pz    = ppz[j+9];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+9] = H_x;
	                        pH_y[j+9] = H_y;
	                        pH_z[j+9] = H_z;
	                        pB_x[j+9] = B_x;
	                        pB_y[j+9] = B_y;
	                        pB_z[j+9] = B_z;
	                 }               
	      }
	      
	      
	           
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4_unroll6x(  const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                 xmm4c4_t phase;
	                 xmm4c4_t refi;
	                 xmm4c4_t px;
	                 xmm4c4_t py;
	                 xmm4c4_t pz;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 omg;
	                 int32_t j,m,m1;
	                 
	                 m = n%6;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                        tht   = ptht[j];
	                        phi   = pphi[j];
	                        omg   = pomg[j];
	                        phase = pphase[j];
	                        refi  = prefi[j];
	                        px    = ppx[j];
	                        py    = ppy[j];
	                        pz    = ppz[j];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j] = H_x;
	                        pH_y[j] = H_y;
	                        pH_z[j] = H_z;
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                    }
	                    if(n<6) { return;}
	                 }                     
	                 
	                 m1 = m+1;
	                 for(j = m1; j != n; j += 6) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T0);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T0);
	                                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T1);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T2);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_NTA);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_NTA);
#endif	     	   	      
                                tht   = ptht[j+0];
	                        phi   = pphi[j+0];
	                        omg   = pomg[j+0];
	                        phase = pphase[j+0];
	                        refi  = prefi[j+0];
	                        px    = ppx[j+0];
	                        py    = ppy[j+0];
	                        pz    = ppz[j+0];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+0] = H_x;
	                        pH_y[j+0] = H_y;
	                        pH_z[j+0] = H_z;
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;
	                        tht   = ptht[j+1];
	                        phi   = pphi[j+1];
	                        omg   = pomg[j+1];
	                        phase = pphase[j+1];
	                        refi  = prefi[j+1];
	                        px    = ppx[j+1];
	                        py    = ppy[j+1];
	                        pz    = ppz[j+1];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+1] = H_x;
	                        pH_y[j+1] = H_y;
	                        pH_z[j+1] = H_z;
	                        pB_x[j+1] = B_x;
	                        pB_y[j+1] = B_y;
	                        pB_z[j+1] = B_z;  
	                        tht   = ptht[j+2];
	                        phi   = pphi[j+2];
	                        omg   = pomg[j+2];
	                        phase = pphase[j+2];
	                        refi  = prefi[j+2];
	                        px    = ppx[j+2];
	                        py    = ppy[j+2];
	                        pz    = ppz[j+2];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+2] = H_x;
	                        pH_y[j+2] = H_y;
	                        pH_z[j+2] = H_z;
	                        pB_x[j+2] = B_x;
	                        pB_y[j+2] = B_y;
	                        pB_z[j+2] = B_z;
	                        tht   = ptht[j+3];
	                        phi   = pphi[j+3];
	                        omg   = pomg[j+3];
	                        phase = pphase[j+3];
	                        refi  = prefi[j+3];
	                        px    = ppx[j+3];
	                        py    = ppy[j+3];
	                        pz    = ppz[j+3];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+3] = H_x;
	                        pH_y[j+3] = H_y;
	                        pH_z[j+3] = H_z;
	                        pB_x[j+3] = B_x;
	                        pB_y[j+3] = B_y;
	                        pB_z[j+3] = B_z;  
	                        tht   = ptht[j+4];
	                        phi   = pphi[j+4];
	                        omg   = pomg[j+4];
	                        phase = pphase[j+4];
	                        refi  = prefi[j+4];
	                        px    = ppx[j+4];
	                        py    = ppy[j+4];
	                        pz    = ppz[j+4];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+4] = H_x;
	                        pH_y[j+4] = H_y;
	                        pH_z[j+4] = H_z;
	                        pB_x[j+4] = B_x;
	                        pB_y[j+4] = B_y;
	                        pB_z[j+4] = B_z; 
	                        tht   = ptht[j+5];
	                        phi   = pphi[j+5];
	                        omg   = pomg[j+5];
	                        phase = pphase[j+5];
	                        refi  = prefi[j+5];
	                        px    = ppx[j+5];
	                        py    = ppy[j+5];
	                        pz    = ppz[j+5];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+5] = H_x;
	                        pH_y[j+5] = H_y;
	                        pH_z[j+5] = H_z;
	                        pB_x[j+5] = B_x;
	                        pB_y[j+5] = B_y;
	                        pB_z[j+5] = B_z;  
	                      
	                 }               
	      }
	      
	      
	          
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4_unroll2x(  const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                 xmm4c4_t phase;
	                 xmm4c4_t refi;
	                 xmm4c4_t px;
	                 xmm4c4_t py;
	                 xmm4c4_t pz;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 omg;
	                 int32_t j,m,m1;
	                 
	                 m = n%2;
	                 if(m!=0) {
	                    for(j = 0; j != m; ++j) {
	                        tht   = ptht[j];
	                        phi   = pphi[j];
	                        omg   = pomg[j];
	                        phase = pphase[j];
	                        refi  = prefi[j];
	                        px    = ppx[j];
	                        py    = ppy[j];
	                        pz    = ppz[j];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j] = H_x;
	                        pH_y[j] = H_y;
	                        pH_z[j] = H_z;
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                    }
	                    if(n<2) { return;}
	                 }                     
	                 
	                 m1 = m+1;
	                 for(j = m1; j != n; j += 2) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T0);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T0);
	                                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T1);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T2);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_NTA);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_NTA);
#endif	     	   	      
                                tht   = ptht[j+0];
	                        phi   = pphi[j+0];
	                        omg   = pomg[j+0];
	                        phase = pphase[j+0];
	                        refi  = prefi[j+0];
	                        px    = ppx[j+0];
	                        py    = ppy[j+0];
	                        pz    = ppz[j+0];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+0] = H_x;
	                        pH_y[j+0] = H_y;
	                        pH_z[j+0] = H_z;
	                        pB_x[j+0] = B_x;
	                        pB_y[j+0] = B_y;
	                        pB_z[j+0] = B_z;
	                        tht   = ptht[j+1];
	                        phi   = pphi[j+1];
	                        omg   = pomg[j+1];
	                        phase = pphase[j+1];
	                        refi  = prefi[j+1];
	                        px    = ppx[j+1];
	                        py    = ppy[j+1];
	                        pz    = ppz[j+1];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j+1] = H_x;
	                        pH_y[j+1] = H_y;
	                        pH_z[j+1] = H_z;
	                        pB_x[j+1] = B_x;
	                        pB_y[j+1] = B_y;
	                        pB_z[j+1] = B_z;  
	                      
	                      
	                 }               
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128 * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm4c4_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm4c4_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 1;
	                 xmm4c4_t phase;
	                 xmm4c4_t refi;
	                 xmm4c4_t px;
	                 xmm4c4_t py;
	                 xmm4c4_t pz;
	                 xmm4c4_t H_x;
	                 xmm4c4_t H_y;
	                 xmm4c4_t H_z;
	                 xmm4c4_t B_x;
	                 xmm4c4_t B_y;
	                 xmm4c4_t B_z;
	                 register __m128 tht;
	                 register __m128 phi;
	                 register __m128 omg;
	                 int32_t j;
	                 
	                 for(j = 0; j != n; ++j) {
#if (__EM_FIELDS_PF_CACHE_HINT__) == 1
	                    _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T0);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T0);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T0);
	                                  
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 2
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T1);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T1);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T1);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 3
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_T2);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_T2);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_T2);
#elif (__EM_FIELDS_PF_CACHE_HINT__) == 4
                            _mm_prefetch((char*)&ptht[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphi[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pomg[j+PF_DIST],_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&pphase[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&prefi[j+PF_DIST].im,_MM_HINT_NTA);    
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppx[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST]re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppy[j+PF_DIST].im,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].re,_MM_HINT_NTA);
	                    _mm_prefetch((char*)&ppz[j+PF_DIST].im,_MM_HINT_NTA);
#endif	     	   	      
                                tht   = ptht[j];
	                        phi   = pphi[j];
	                        omg   = pomg[j];
	                        phase = pphase[j];
	                        refi  = prefi[j];
	                        px    = ppx[j];
	                        py    = ppy[j];
	                        pz    = ppz[j];
	                        B_XYZ_H_XYZ_EP_xmm4c4(tht,phi,omg,
	                                               phase,refi,px,
	                                               py,pz,
	                                               H_x,H_y,H_z,
	                                               B_x,B_y,B_z);
	                        pH_x[j] = H_x;
	                        pH_y[j] = H_y;
	                        pH_z[j] = H_z;
	                        pB_x[j] = B_x;
	                        pB_y[j] = B_y;
	                        pB_z[j] = B_z;
	                                              
	               }               
	      }
	      
	      
	      
	      
	      
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4_a(const float * __restrict __ATTR_ALIGN__(16) ptht,
	                                         const float * __restrict __ATTR_ALIGN__(16) pphi,
	                                         const float * __restrict __ATTR_ALIGN__(16) pomg,
	                                         const xmm4c4_t phase,
	                                         const xmm4c4_t refi,
	                                         const xmm4c4_t px,
	                                         const xmm4c4_t py,
	                                         const xmm4c4_t pz,
	                                         xmm4c4_t & H_x,
	                                         xmm4c4_t & H_y,
	                                         xmm4c4_t & H_z,
	                                         xmm4c4_t & B_x,
	                                         xmm4c4_t & B_y,
	                                         xmm4c4_t & B_z) {
	                         
	               register __m128 tht = _mm_load_ps(&ptht[0]);
	               register __m128 phi = _mm_load_ps(&pphi[0]);
	               register __m128 omg = _mm_load_ps(&pomg[0]);            
	               const __m128 c   = _mm_set1_ps(299792458.0); 
	               const __m128 mu0 = _mm_set1_ps(0.0000012566370614359173f);   
	               const __m128 psi0 = _mm_setzero_ps();
	               const __m128 C00  = _mm_setzero_ps();
	               
	               xmm4c4_t H_x_1;
	               xmm4c4_t H_y_1;
	               xmm4c4_t H_z_1;
	               xmm4c4_t H_x_2;
	               xmm4c4_t H_y_2;
	               xmm4c4_t H_z_2;
	               xmm4c4_t k;
	               xmm4c4_t t0;
	               xmm4c4_t cdirx;
	               xmm4c4_t cdiry;
	               xmm4c4_t cdirz;
	               
	               register __m128 vpolx;
	               register __m128 vpoly;
	               register __m128 vpolz;
	               register __m128 vdirx;
	               register __m128 vdiry;
	               register __m128 vdirz;
	               register __m128 cn;
	               register __m128 x0;
	               register __m128 t0r,t0i;
	               register __m128 t1r,t1i;
	               register __m128 t2r,t2i;
	               
	               dir_vec_xmm4c4(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm_div_ps(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm_mul_ps(refi.re,x0);
	               k.im     = _mm_mul_ps(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_xmm4c4(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm_add_ps(H_x_1.re,t0r);
	               H_x.im = _mm_add_ps(H_x_1.im,t0i);
	               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm_add_ps(H_y_1.re,t1r);
	               H_y.im = _mm_add_ps(H_y_1.im,t1i);
	               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm_add_ps(H_z_2.re,t2r);
	              H_z.im = _mm_add_ps(H_z_2.im,t2i);
	              
	              cn = cnorm_xmm4c4(H_x,H_y,H_z);
	              
	              x0     = _mm_div_ps(omg,mu0);
	              H_x.re = _mm_div_ps(H_x.re,cn);
	              H_x.im = _mm_div_ps(H_x.im,cn);
	              H_y.re = _mm_div_ps(H_y.re,cn);
	              H_y.im = _mm_div_ps(H_y.im,cn);
	              H_z.re = _mm_div_ps(H_z.re,cn);
	              H_z.im = _mm_div_ps(H_z.im,cn); 
	              
	              t0.re  = _mm_div_ps(k.re,x0);
	              t0.im  = _mm_div_ps(k.im,x0);
	              
	              scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_xmm4c4_u(const float * __restrict  ptht,
	                                         const float * __restrict  pphi,
	                                         const float * __restrict  pomg,
	                                         const xmm4c4_t phase,
	                                         const xmm4c4_t refi,
	                                         const xmm4c4_t px,
	                                         const xmm4c4_t py,
	                                         const xmm4c4_t pz,
	                                         xmm4c4_t & H_x,
	                                         xmm4c4_t & H_y,
	                                         xmm4c4_t & H_z,
	                                         xmm4c4_t & B_x,
	                                         xmm4c4_t & B_y,
	                                         xmm4c4_t & B_z) {
	                         
	               register __m128 tht = _mm_loadu_ps(&ptht[0]);
	               register __m128 phi = _mm_loadu_ps(&pphi[0]);
	               register __m128 omg = _mm_loadu_ps(&pomg[0]);            
	               const __m128 c   = _mm_set1_ps(299792458.0); 
	               const __m128 mu0 = _mm_set1_ps(0.0000012566370614359173f);   
	               const __m128 psi0 = _mm_setzero_ps();
	               const __m128 C00  = _mm_setzero_ps();
	               
	               xmm4c4_t H_x_1;
	               xmm4c4_t H_y_1;
	               xmm4c4_t H_z_1;
	               xmm4c4_t H_x_2;
	               xmm4c4_t H_y_2;
	               xmm4c4_t H_z_2;
	               xmm4c4_t k;
	               xmm4c4_t t0;
	               xmm4c4_t cdirx;
	               xmm4c4_t cdiry;
	               xmm4c4_t cdirz;
	               
	               register __m128 vpolx;
	               register __m128 vpoly;
	               register __m128 vpolz;
	               register __m128 vdirx;
	               register __m128 vdiry;
	               register __m128 vdirz;
	               register __m128 cn;
	               register __m128 x0;
	               register __m128 t0r,t0i;
	               register __m128 t1r,t1i;
	               register __m128 t2r,t2i;
	               
	               dir_vec_xmm4c4(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm_div_ps(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm_mul_ps(refi.re,x0);
	               k.im     = _mm_mul_ps(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_xmm4c4(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_xmm4c4(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm_add_ps(H_x_1.re,t0r);
	               H_x.im = _mm_add_ps(H_x_1.im,t0i);
	               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm_add_ps(H_y_1.re,t1r);
	               H_y.im = _mm_add_ps(H_y_1.im,t1i);
	               
	               cmul_xmm4c4(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm_add_ps(H_z_2.re,t2r);
	              H_z.im = _mm_add_ps(H_z_2.im,t2i);
	              
	              cn = cnorm_xmm4c4(H_x,H_y,H_z);
	              
	              x0     = _mm_div_ps(omg,mu0);
	              H_x.re = _mm_div_ps(H_x.re,cn);
	              H_x.im = _mm_div_ps(H_x.im,cn);
	              H_y.re = _mm_div_ps(H_y.re,cn);
	              H_y.im = _mm_div_ps(H_y.im,cn);
	              H_z.re = _mm_div_ps(H_z.re,cn);
	              H_z.im = _mm_div_ps(H_z.im,cn); 
	              
	              t0.re  = _mm_div_ps(k.re,x0);
	              t0.im  = _mm_div_ps(k.im,x0);
	              
	              scrossc_xmm4c4(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_xmm4c4(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	                                      
	       
	       
	       
	     
                
                
        } // radiolocation

} // gms


#endif /*__GMS_EM_FIELDS_XMM4R4_HPP__*/
