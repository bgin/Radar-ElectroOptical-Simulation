

#ifndef __GMS_EM_FIELDS_XMM4R4_H__
#define __GMS_EM_FIELDS_XMM4R4_H__ 301020231146

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
	        
	        
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm4c4_unroll16x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	      
	      
	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void sdotv_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	      
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	      
	      
	         
	      
	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n);
	                                        
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128 * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128 * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	      
	      
	        
	        
	           __ATTR_ALWAYS_INLINE__
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
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void cdotv_xmm4c4_unroll16x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	       
	       
	       
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void cdotv_xmm4c4_unroll10x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	       
	        
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                void cdotv_xmm4c4_unroll6x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	       
	       
	          
	       
	       
	       
	     
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cdotv_xmm4c4_unroll2x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	       
	       
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
               	   void cdotv_xmm4c4_rolled(    const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	       
	       
	       
	        
	       
	       
	       
	        
	        
	           __ATTR_ALWAYS_INLINE__
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
	       
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void cnorm_xmm4c4_unroll16x( const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	      
	         
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                void cnorm_xmm4c4_unroll10x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	      
	      
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cnorm_xmm4c4_unroll6x(const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	         
	      
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cnorm_xmm4c4_unroll2x( const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void cnorm_xmm4c4_rolled(    const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm4c4_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm4c4_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	      
	      
	      
	                                       
	       
	       
	           __ATTR_ALWAYS_INLINE__
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
	          
	          
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	                                          
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST);
	        
	          
	        
	        
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	        
	        
	         
	        
	        
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	                                          
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	        
	        
  
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
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
	         
	         
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST);
	         
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	         
	        
	         
	         
	         
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	         
	          
	         
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	         
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                          int32_t & PF_DIST); 
	         
	         
	         
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
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
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void dir_vec_xmm4c4_unroll16x( const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	       
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void dir_vec_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	       
	        
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void dir_vec_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	          
	       
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void dir_vec_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	       
	       
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                    void dir_vec_xmm4c4_rolled(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	       
	         
	         
	           __ATTR_ALWAYS_INLINE__
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
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void pol_vec_xmm4c4_unroll16x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void pol_vec_xmm4c4_unroll10x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	      
	      
	          
	      
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                void pol_vec_xmm4c4_unroll6x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	          
	        
	      
	          
	          
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void pol_vec_xmm4c4_unroll2x(const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	          
	          
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void pol_vec_xmm4c4_rolled(    const __m128 * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128 * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128 * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	          
	      
	      
	          __ATTR_ALWAYS_INLINE__
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
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DISPATCH); 
	             
	      
	         
	      
	      
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DISPATCH); 
	         
	      
	        
	        
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DISPATCH); 
	                                           
	                                           
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DISPATCH);
	                                            
	                                           
	                                           
	              
	      
	      
	        
	           __ATTR_ALWAYS_INLINE__
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
	     
	     
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DIST);
	         
	         
	         
	         
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DIST);
	         
	        
	         
	         
	         
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DIST); 
	         
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                           int32_t & PF_DIST); 
	         
	                  
	         
	         
	     
	           __ATTR_ALWAYS_INLINE__
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
	       
	       
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                int32_t & PF_DIST); 
	      
	      
	          
	      
	      
	      
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                int32_t & PF_DIST);
	      
	        
	      
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                int32_t & PF_DIST); 
	                                                
	      
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                int32_t & PF_DIST); 
	      
	       
	       
	           __ATTR_ALWAYS_INLINE__
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
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                 int32_t & PF_DIST); 
	      
	           
	      
	      
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                 int32_t & PF_DIST); 
	          
	      
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                 int32_t & PF_DIST); 
	                                                 
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
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
	                                                 int32_t & PF_DIST); 
	      
	      
	      
	      
	        
	        
	           __ATTR_ALWAYS_INLINE__
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


#endif /*__GMS_EM_FIELDS_XMM4R4_H__*/
