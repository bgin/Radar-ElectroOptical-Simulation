

#ifndef __GMS_EM_FIELDS_XMM2R8_HPP__
#define __GMS_EM_FIELDS_XMM4R8_HPP__ 301020230900

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

    const unsigned int GMS_EM_FIELDS_XMM2R8_MAJOR = 1U;
    const unsigned int GMS_EM_FIELDS_XMM2R8_MINOR = 0U;
    const unsigned int GMS_EM_FIELDS_XMM2R8_MICRO = 0U;
    const unsigned int GMS_EM_FIELDS_XMM2R8_FULLVER =
      1000U*GMS_EM_FIELDS_XMM2R8_MAJOR+
      100U*GMS_EM_FIELDS_XMM2R8_MINOR+
      10U*GMS_EM_FIELDS_XMM2R8_MICRO;
    const char * const GMS_EM_FIELDS_XMM2R8_CREATION_DATE = "30-10-2023 09:00 AM +00200 (MON 30 10 2023 GMT+2)";
    const char * const GMS_EM_FIELDS_XMM2R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EM_FIELDS_XMM2R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EM_FIELDS_XMM2R8_DESCRIPTION   = " Computational ElectroMagnetics related helper routines."
                       

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_complex_xmm2r8.hpp"

#ifndef __EM_FIELDS_PF_CACHE_HINT__
#define __EM_FIELDS_PF_CACHE_HINT__ 1
#endif 

namespace gms {



          namespace radiolocation {
          
          
 
              
               
                
                   __ATTR_ALWAYS_INLINE__
	           static inline
	           __m128d sdotv_xmm2c8(const __m128d v1x,
	                                const __m128d v1y,
	                                const __m128d v1z,
	                                const __m128d v2x,
	                                const __m128d v2y,
	                                const __m128d v2z) {
	                                
	                  register __m128d result;
	                  result = _mm_fmadd_pd(v1x,v2x,
	                                      _mm_fmadd_pd(v1y,v2y,
	                                                 _mm_mul_pd(v1z,v2z)));
	                  return (result);                       
	        }
	        
	        
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                    void sdotv_xmm2c8_unroll16x(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128d * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	      
	      
	      
	      
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void sdotv_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128d * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm2c8_unroll6x(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128d * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	      
	      
	         
	      
	      
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm2c8_unroll2x(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128d * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void sdotv_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        __m128d * __restrict __ATTR_ALIGN__(16) pdtv,
	                                        const int32_t n);
	        
                   __ATTR_ALWAYS_INLINE__
	           static inline
	           __m128d sdotv_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) pv1x,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv1y,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv1z,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv2x,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv2y,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv2z) {
	                          
	                  register __m128d v1x = _mm_load_pd(&pv1x[0]);
	                  register __m128d v1y = _mm_load_pd(&pv1y[0]);  
	                  register __m128d v1z = _mm_load_pd(&pv1z[0]); 
	                  register __m128d v2x = _mm_load_pd(&pv2x[0]);  
	                  register __m128d v2y = _mm_load_pd(&pv2y[0]); 
	                  register __m128d v2z = _mm_load_pd(&pv2z[0]);
	                  register __m128d result;
	                  result = _mm_fmadd_pd(v1x,v2x,
	                                      _mm_fmadd_pd(v1y,v2y,
	                                                 _mm_mul_pd(v1z,v2z)));
	                  return (result);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           __m128d sdotv_xmm2c8_u(const double * __restrict  pv1x,
	                                const double * __restrict  pv1y,
	                                const double * __restrict  pv1z,
	                                const double * __restrict  pv2x,
	                                const double * __restrict  pv2y,
	                                const double * __restrict  pv2z) {
	                          
	                  register __m128d v1x = _mm_loadu_pd(&pv1x[0]);
	                  register __m128d v1y = _mm_loadu_pd(&pv1y[0]);  
	                  register __m128d v1z = _mm_loadu_pd(&pv1z[0]); 
	                  register __m128d v2x = _mm_loadu_pd(&pv2x[0]);  
	                  register __m128d v2y = _mm_loadu_pd(&pv2y[0]); 
	                  register __m128d v2z = _mm_loadu_pd(&pv2z[0]);
	                  register __m128d result;
	                  result = _mm_fmadd_pd(v1x,v2x,
	                                      _mm_fmadd_pd(v1y,v2y,
	                                                 _mm_mul_pd(v1z,v2z)));
	                  return (result);                       
	        }
	        
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void cdotv_xmm2c8( const xmm2c8_t v1x,
	                              const xmm2c8_t v1y,
	                              const xmm2c8_t v1z,
	                              const xmm2c8_t v2x,
	                              const xmm2c8_t v2y,
	                              const xmm2c8_t v2z,
	                              xmm2c8_t & res) {
	                              
	                xmm2c8_t tx,ty,tz;
	                cmul_xmm2c8(v1x.re,v1x.im,v2x.re,
	                                  v2x.im,&tx.re,&tx.im); 
	                cmul_xmm2c8(v1y.re,v1y.im,v2y.re,
	                                  v2y.im,&ty.re,&ty.im);
	                cmul_xmm2c8(v1z.re,v1z.im,v2z.re,
	                                  v2z.im,&tz.re,&tz.im);
	                res.re = _mm_add_pd(tx.re,
	                                   _mm_add_pd(ty.re,tz.re));
	                res.im = _mm_add_pd(tx.im,
	                                   _mm_add_pd(ty.im,tz.im));                   
	        }
	        
	        
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cdotv_xmm2c8_unroll16x(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void cdotv_xmm2c8_unroll10x(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cdotv_xmm2c8_unroll6x(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cdotv_xmm2c8_unroll2x(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cdotv_xmm2c8_rolled(    const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	              	       	       
	      	      	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_VECTORCALL__
	           static inline
	           __m128d cnorm_xmm2c8(const xmm2c8_t vx,
	                                const xmm2c8_t vy,
	                                const xmm2c8_t vz) {
	                                
	                  xmm2c8_t t,cx,cy,cz;
	                  __m128d vs;
	                  cconj_xmm2c8_v2(vx.re,vx.im,&cx.re,&cx.im);
	                  cconj_xmm2c8_v2(vy.re,vy.im,&cy.re,&cy.im);
	                  cconj_xmm2c8_v2(vz.re,vz.im,&cz.re,&cz.im);
	                  cdotv_xmm2c8(vx,vy,vz,cx,cy,cz,t);
	                  vs = _mm_sqrt_pd(t.re);
	                  return (vs);                      
	       }
	       
	       
	       
	       
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cnorm_xmm2c8_unroll16x( const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	      
	         
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cnorm_xmm2c8_unroll10x(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cnorm_xmm2c8_unroll6x(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void cnorm_xmm2c8_unroll2x( const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	                                        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void cnorm_xmm2c8_rolled(    const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                        const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                        xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                        const int32_t n,
	                                        int32_t & PF_DIST); 
	      
  	      	                                       	       
	       
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void scrossc_xmm2c8(const xmm2c8_t v1x,
	                                const xmm2c8_t v1y,
	                                const xmm2c8_t v1z,
	                                const xmm2c8_t v2x,
	                                const xmm2c8_t v2y,
	                                const xmm2c8_t v2z,
	                                xmm2c8 & resx,
	                                xmm2c8 & resy,
	                                xmm2c8 & resz) {
	                                
	                 xmm2c8_t t0,t1,t2,t3,t4,t5,t6;
	                 cmul_xmm2c8(v1y.re,v1y.im,v2z.re,
	                              v2z.im,&t0.re,&t0.im); 
	                 cmul_xmm2c8(v1z.re,v1z.im,v2y.re,
	                              v2y.im,&t1.re,&t1.im);
	                 resx.re = _mm_sub_pd(t0.re,t1.re);
	                 resx.im = _mm_sub_pd(t0.im,t1.im);
	                 cmul_xmm2c8(v1z.re,v1z.im,v2x.re,
	                              v2x.im,&t2.re,&t2.im);
	                 cmul_xmm2c8(v1x.re,v1x.im,v2z.re,
	                              v2z.im,&t3.re,&t3.im);
	                 resy.re = _mm_sub_pd(t2.re,t3.re);
	                 resy.im = _mm_sub_pd(t2.im,t3.im);
	                 cmul_xmm2c8(v1x.re,v1x.im,v2y.re,
	                              v2y.im,&t4.re,&t4.im);
	                 cmul_xmm2c8(v1y.re,v1y.im,v2x.re,
	                              v2x.im,&t5.re,&t5.im);    
	                 resz.re = _mm_sub_pd(t4.re,t5.re);
	                 resz.im = _mm_sub_pd(t4.im,t5.im);
	          }
	          
	          
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void scrossc_xmm2c8_unroll16x( const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void scrossc_xmm2c8_unroll10x( const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void scrossc_xmm2c8_unroll6x(const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void scrossc_xmm2c8_unroll2x(const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void scrossc_xmm2c8_rolled(    const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                          xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void scrossv_xmm2c8(const __m128d v1x,
	                                const __m128d v1y,
	                                const __m128d v1z,
	                                const __m128d v2x,
	                                const __m128d v2y,
	                                const __m128d v2z,
	                                __m128d * __restrict vcx,
	                                __m128d * __restrict vcy,
	                                __m128d * __restrict vcz) {
	                                
	                *vcx = _mm_fmsub_pd(v1y,v2z,
	                                   _mm_mul_pd(v1x,v2y));
	                *vcy = _mm_fmsub_pd(v1z,v2x,
	                                   _mm_mul_pd(v1x,v2z));
	                *vcz = _mm_fmsub_pd(v1x,v2y,
	                                   _mm_mul_pd(v1y,v2x));
	         }
	         
	         
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void scrossv_xmm2c8_unroll16x( const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void scrossv_xmm2c8_unroll10x( const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void scrossv_xmm2c8_unroll6x(  const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void scrossv_xmm2c8_unroll2x(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                    void scrossv_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	         	           	         	         
	         
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void scrossv_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) pv1x,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv1y,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv1z,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv2x,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv2y,
	                                  const double * __restrict __ATTR_ALIGN__(16) pv2z,
	                                  double * __restrict __ATTR_ALIGN__(16) vcx,
	                                  double * __restrict __ATTR_ALIGN__(16) vcy,
	                                  double * __restrict __ATTR_ALIGN__(16) vcz) {
	                      
	                 register __m128d v1x = _mm_load_pd(&pv1x[0]);
	                 register __m128d v1y = _mm_load_pd(&pv1y[0]);
	                 register __m128d v1z = _mm_load_pd(&pv1z[0]);
	                 register __m128d v2x = _mm_load_pd(&pv2x[0]);
	                 register __m128d v2y = _mm_load_pd(&pv2y[0]);
	                 register __m128d v2z = _mm_load_pd(&pv2z[0]);          
	                *vcx = _mm_fmsub_pd(v1y,v2z,
	                                   _mm_mul_pd(v1x,v2y));
	                *vcy = _mm_fmsub_pd(v1z,v2x,
	                                   _mm_mul_pd(v1x,v2z));
	                *vcz = _mm_fmsub_pd(v1x,v2y,
	                                   _mm_mul_pd(v1y,v2x));
	         }
	         
	         
	         
	          
	           static inline
	           void scrossv_xmm2c8_u(const double * __restrict pv1x,
	                                  const double * __restrict pv1y,
	                                  const double * __restrict pv1z,
	                                  const double * __restrict pv2x,
	                                  const double * __restrict pv2y,
	                                  const double * __restrict pv2z,
	                                  double * __restrict vcx,
	                                  double * __restrict vcy,
	                                  double * __restrict vcz) {
	                      
	                 register __m128d v1x = _mm_loadu_pd(&pv1x[0]);
	                 register __m128d v1y = _mm_loadu_pd(&pv1y[0]);
	                 register __m128d v1z = _mm_loadu_pd(&pv1z[0]);
	                 register __m128d v2x = _mm_loadu_pd(&pv2x[0]);
	                 register __m128d v2y = _mm_loadu_pd(&pv2y[0]);
	                 register __m128d v2z = _mm_loadu_pd(&pv2z[0]);          
	                *vcx = _mm_fmsub_pd(v1y,v2z,
	                                   _mm_mul_pd(v1x,v2y));
	                *vcy = _mm_fmsub_pd(v1z,v2x,
	                                   _mm_mul_pd(v1x,v2z));
	                *vcz = _mm_fmsub_pd(v1x,v2y,
	                                   _mm_mul_pd(v1y,v2x));
	         }
	         
	         
	         //! Direction Vector spherical [theta,phi] (SIMD data-types)
	         
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void dir_vec_xmm2c8(  const __m128d tht,
	                                  const __m128d phi,
	                                  __m128d * __restrict dvx,
	                                  __m128d * __restrict dvy,
	                                  __m128d * __restrict dvz) {
	                  
	                        
	                register __m128d stht,cphi,sphi,ctht;
	                cphi = _mm_cos_pd(phi);
	                stht = _mm_sin_pd(tht);
	                *dvx = _mm_mul_pd(stht,cphi);
	                sphi = _mm_sin_pd(phi);
	                *dvy = _mm_mul_pd(stht,sphi);
	                ctht = _mm_cos_pd(tht);
	                *dvz = ctht;                       
	        }
	        
	        
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void dir_vec_xmm2c8_unroll16x( const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void dir_vec_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void dir_vec_xmm2c8_unroll6x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void dir_vec_xmm2c8_unroll2x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void dir_vec_xmm2c8_rolled(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	          __ATTR_ALWAYS_INLINE__
	           static inline
	           void dir_vec_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) ptht,
	                                  const double * __restrict __ATTR_ALIGN__(16) pphi,
	                                  double * __restrict __ATTR_ALIGN__(16) dvx,
	                                  double * __restrict __ATTR_ALIGN__(16) dvy,
	                                  double * __restrict __ATTR_ALIGN__(16) dvz) {
	                  
	                register __m128d tht = _mm_load_pd(&ptht[0]);
	                register __m128d phi = _mm_load_pd(&pphi[0]);              
	                register __m128d stht,cphi,sphi,ctht;
	                cphi = _mm_cos_pd(phi);
	                stht = _mm_sin_pd(tht);
	                _mm_store_pd(&dvx[0] , _mm_mul_pd(stht,cphi));
	                sphi = _mm_sin_pd(phi);
	                _mm_store_pd(&dvy[0] , _mm_mul_pd(stht,sphi));
	                ctht = _mm_cos_pd(tht);
	                _mm_store_pd(&dvz[0] , ctht);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void dir_vec_xmm2c8_u(const double * __restrict  ptht,
	                                  const double * __restrict  pphi,
	                                  double * __restrict  dvx,
	                                  double * __restrict  dvy,
	                                  double * __restrict  dvz) {
	                  
	                register __m128d tht = _mm_loadu_pd(&ptht[0]);
	                register __m128d phi = _mm_loadu_pd(&pphi[0]);              
	                register __m128d stht,cphi,sphi,ctht;
	                cphi = _mm_cos_pd(phi);
	                stht = _mm_sin_pd(tht);
	                _mm_storeu_pd(&dvx[0] , _mm_mul_pd(stht,cphi));
	                sphi = _mm_sin_pd(phi);
	                _mm_storeu_pd(&dvy[0] , _mm_mul_pd(stht,sphi));
	                ctht = _mm_cos_pd(tht);
	                _mm_storeu_pd(&dvz[0] , ctht);                       
	        }
	        
	        
	         //! Polarization Vector of plane-wave propagating into direction computed by
                 //! dir_vector_xmmxrx (SIMD data-types)
                 
                   __ATTR_ALWAYS_INLINE__
	           static inline
	           void pol_vec_xmm2c8(const __m128d tht,
	                                const __m128d phi,
	                                const __m128d psi,
	                                __m128d * __restrict pvx,
	                                __m128d * __restrict pvy,
	                                __m128d * __restrict pvz) {
	                 
	                using namespace gms::math               
	                register __m128d cpsi,cphi,spsi,sphi,t0;
	                cpsi = _mm_cos_pd(psi);
	                cphi = _mm_cos_pd(phi);
	                spsi = _mm_sin_pd(psi);
	                sphi = _mm_sin_pd(phi);
	                t0   = _mm_mul_pd(spsi,_mm_cos_pd(tht));
	                *pvx = _mm_fmsub_pd(cpsi,sphi,
	                                   _mm_mul_pd(t0,cphi));
	                *pvy = _mm_fmsub_pd(negate_xmm2c8(cpsi),cphi,
	                                                    _mm_mul_pd(t0,sphi));
	                *pvz = _mm_mul_pd(spsi,_mm_sin_pd(tht));                         
	      }
	      
	      
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void pol_vec_xmm2c8_unroll16x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	       
	      	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void pol_vec_xmm2c8_unroll6x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	             	      	          
	          
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void pol_vec_xmm2c8_unroll2x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                    void pol_vec_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	          
	      
	      
	          __ATTR_ALWAYS_INLINE__
	          static inline
	           void pol_vec_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) ptht,
	                                  const double * __restrict __ATTR_ALIGN__(16) pphi,
	                                  const double * __restrict __ATTR_ALIGN__(16) psi,
	                                  double * __restrict __ATTR_ALIGN__(16) pvx,
	                                  double * __restrict __ATTR_ALIGN__(16) pvy,
	                                  double * __restrict __ATTR_ALIGN__(16) pvz) {
	                 
	                 using namespace gms::math     
	                register __m128d tht = _mm_load_pd(&ptht[0]);
	                register __m128d phi = _mm_load_pd(&pphi[0]);  
	                register __m128d psi = _mm_load_pd(&ppsi[0]);           
	                register __m128d cpsi,cphi,spsi,sphi,t0;
	                cpsi = _mm_cos_pd(psi);
	                cphi = _mm_cos_pd(phi);
	                spsi = _mm_sin_pd(psi);
	                sphi = _mm_sin_pd(phi);
	                t0   = _mm_mul_pd(spsi,_mm_cos_pd(tht));
	                _mm_store_pd(&pvx[0] ,_mm_fmsub_pd(cpsi,sphi,
	                                   _mm_mul_pd(t0,cphi)));
	                _mm_store_pd(&pvy[0] ,_mm_fmsub_pd(negate_xmm2c8(cpsi),cphi,
	                                                    _mm_mul_pd(t0,sphi)));
	                _mm_store_pd(&pvz[0] ,_mm_mul_pd(spsi,_mm_sin_pd(tht)));                         
	      } 
	        
	        
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void pol_vec_xmm2c8_u(const double * __restrict  ptht,
	                                  const double * __restrict  pphi,
	                                  const double * __restrict  psi,
	                                  double * __restrict  pvx,
	                                  double * __restrict  pvy,
	                                  double * __restrict  pvz) {
	                 
	                  using namespace gms::math    
	                register __m128d tht = _mm_loadu_pd(&ptht[0]);
	                register __m128d phi = _mm_loadu_pd(&pphi[0]);  
	                register __m128d psi = _mm_loadu_pd(&ppsi[0]);           
	                register __m128d cpsi,cphi,spsi,sphi,t0;
	                cpsi = _mm_cos_pd(psi);
	                cphi = _mm_cos_pd(phi);
	                spsi = _mm_sin_pd(psi);
	                sphi = _mm_sin_pd(phi);
	                t0   = _mm_mul_pd(spsi,_mm_cos_pd(tht));
	                _mm_storeu_pd(&pvx[0] ,_mm_fmsub_pd(cpsi,sphi,
	                                   _mm_mul_pd(t0,cphi)));
	                _mm_storeu_pd(&pvy[0] ,_mm_fmsub_pd(negate_xmm2c8(cpsi),cphi,
	                                                    _mm_mul_pd(t0,sphi)));
	                _mm_storeu_pd(&pvz[0] ,_mm_mul_pd(spsi,_mm_sin_pd(tht)));                         
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
	           void H_XYZ_VP_xmm2c8(const __m128d vpolx,
	                                 const __m128d vpoly,
	                                 const __m128d vpolz,
	                                 const __m128d vdirx,
	                                 const __m128d vdiry,
	                                 const __m128d vdirz,
	                                 const __m128d vrx,
	                                 const __m128d vry,
	                                 const __m128d vrz,
	                                 const xmm2c8_t k,
	                                 xmm2c8_t & H_x,
	                                 xmm2c8_t & H_y,
	                                 xmm2c8_t & H_z) {
	               
	               	register __m128d dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_xmm2c8(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm_set1_pd(1.0);
	                ir = _mm_setzero_pd();
	                cmul_xmm2c8(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm_mul_pd(dp,cer);
	                cei = _mm_mul_pd(dp,cei);
	                cexp_xmm2c8(cer,cei,&expr,&expi);
	                H_x.re = _mm_mul_pd(vpolx,expr);
	                H_x.im = _mm_mul_pd(vpolx,expi);
	                H_y.re = _mm_mul_pd(vpoly,expr);
	                H_y.im = _mm_mul_pd(vpoly,expi);
	                H_z.re = _mm_mul_pd(vpolz,expr);
	                H_z.im = _mm_mul_pd(vpolz,expi);
	        }
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                    void H_XYZ_VP_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH); 
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void H_XYZ_VP_xmm2c8_unroll6x( const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH); 
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void H_XYZ_VP_xmm2c8_unroll2x( const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH); 
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void H_XYZ_VP_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN_(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                           const int32_t n,
	                                           int32_t & PF_DISPATCH); 
	      
	        
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void H_XYZ_VP_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) vpolx,
	                                 const double * __restrict __ATTR_ALIGN__(16) vpoly,
	                                 const double * __restrict __ATTR_ALIGN__(16) vpolz,
	                                 const double * __restrict __ATTR_ALIGN__(16) vdirx,
	                                 const double * __restrict __ATTR_ALIGN__(16) vdiry,
	                                 const double * __restrict __ATTR_ALIGN__(16) vdirz,
	                                 const double * __restrict __ATTR_ALIGN__(16) vrx,
	                                 const double * __restrict __ATTR_ALIGN__(16) vry,
	                                 const double * __restrict __ATTR_ALIGN__(16) vrz,
	                                 const xmm2c8_t k,
	                                 xmm2c8_t & H_x,
	                                 xmm2c8_t & H_y,
	                                 xmm2c8_t & H_z) {
	               
	                register __m128d vpolx = _mm_load_pd(&vpolx[0]);
	                register __m128d vpoly = _mm_load_pd(&vpoly[0]);
	                register __m128d vpolz = _mm_load_pd(&vpolz[0]);
	                register __m128d vdirx = _mm_load_pd(&vdirx[0]);
	                register __m128d vdiry = _mm_load_pd(&vdiry[0]);
	                register __m128d vdirz = _mm_load_pd(&vdirz[0]);
	                register __m128d vrx   = _mm_load_pd(&vrx[0]);
	                register __m128d vry   = _mm_load_pd(&vry[0]);
	                register __m128d vrz   = _mm_load_pd(&vrz[0]);
	               	__m128d dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_xmm2c8(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm_set1_pd(1.0f);
	                ir = _mm_setzero_pd();
	                cmul_xmm2c8(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm_mul_pd(dp,cer);
	                cei = _mm_mul_pd(dp,cei);
	                cexp_xmm2c8(cer,cei,&expr,&expi);
	                H_x.re = _mm_mul_pd(vpolx,expr);
	                H_x.im = _mm_mul_pd(vpolx,expi);
	                H_y.re = _mm_mul_pd(vpoly,expr);
	                H_y.im = _mm_mul_pd(vpoly,expi);
	                H_z.re = _mm_mul_pd(vpolz,expr);
	                H_z.im = _mm_mul_pd(vpolz,expi);
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void H_XYZ_VP_xmm2c8_u(const double * __restrict  vpolx,
	                                 const double * __restrict  vpoly,
	                                 const double * __restrict  vpolz,
	                                 const double * __restrict  vdirx,
	                                 const double * __restrict  vdiry,
	                                 const double * __restrict  vdirz,
	                                 const double * __restrict  vrx,
	                                 const double * __restrict  vry,
	                                 const double * __restrict  vrz,
	                                 const xmm2c8_t k,
	                                 xmm2c8_t & H_x,
	                                 xmm2c8_t & H_y,
	                                 xmm2c8_t & H_z) {
	               
	                register __m128d vpolx = _mm_loadu_pd(&vpolx[0]);
	                register __m128d vpoly = _mm_loadu_pd(&vpoly[0]);
	                register __m128d vpolz = _mm_loadu_pd(&vpolz[0]);
	                register __m128d vdirx = _mm_loadu_pd(&vdirx[0]);
	                register __m128d vdiry = _mm_loadu_pd(&vdiry[0]);
	                register __m128d vdirz = _mm_loadu_pd(&vdirz[0]);
	                register __m128d vrx   = _mm_loadu_pd(&vrx[0]);
	                register __m128d vry   = _mm_loadu_pd(&vry[0]);
	                register __m128d vrz   = _mm_loadu_pd(&vrz[0]);
	               	__m128d dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_xmm2c8(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm_set1_pd(1.0);
	                ir = _mm_setzero_pd();
	                cmul_xmm2c8(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm_mul_pd(dp,cer);
	                cei = _mm_mul_pd(dp,cei);
	                cexp_xmm2c8(cer,cei,&expr,&expi);
	                H_x.re = _mm_mul_pd(vpolx,expr);
	                H_x.im = _mm_mul_pd(vpolx,expi);
	                H_y.re = _mm_mul_pd(vpoly,expr);
	                H_y.im = _mm_mul_pd(vpoly,expi);
	                H_z.re = _mm_mul_pd(vpolz,expr);
	                H_z.im = _mm_mul_pd(vpolz,expi);
	        }
	        
	        
	        /*
	             
     ! Magnetic Field (SIMD data-types) [plane-wave], polarization 'vpol' of
     !  wave-vector argument:  vdir*k at sixteen points 'r'.
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void B_XYZ_VP_xmm2c8(const __m128d vpolx,
	                                 const __m128d vpoly,
	                                 const __m128d vpolz,
	                                 const __m128d vdirx,
	                                 const __m128d vdiry,
	                                 const __m128d vdirz,
	                                 const xmm2c8_t k,
	                                 const __m128d omega,
	                                 const __m128d vrx,
	                                 const __m128d vry,
	                                 const __m128d vrz,
	                                 xmm2c8_t & B_x,
	                                 xmm2c8_t & B_y,
	                                 xmm2c8_t & B_z) {
	                                 
	                const __m128d mu0 = _mm_set1_pd(0.0000012566370614359173);
	                xmm2c8_t cdirx;
	                xmm2c8_t cdiry;
	                xmm2c8_t cdirz;
	                xmm2c8_t H_x;
	                xmm2c8_t H_y;
	                xmm2c8_t H_z;
	                xmm2c8_t cpx;
	                xmm2c8_t cpy;
	                xmm2c8_t cpz;
	                xmm2c8_t t0;
	                __m128d zz0;
	                H_XYZ_VP_xmm2c8(vpolx,vpoy,vpolz,
	                               	 vdirx,vdiry,vdirz,
	                                 vrx,vry,vrz,
	                                 H_x,H_y,H_z);
	                                 	                                 
	                cdirx.re = vdirx;
	                cdirx.im = _mm_setzero_pd();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm_mul_pd(omega,mu0);
	                t0.re    = _mm_div_pd(k.re,zz0);
	                t0.im    = _mm_div_pd(k.im,zz0);
	                
	                scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                                H_x,H_y,H_z,
                                        cpx,cpy,cpz);                
	                     	                                
	                cmul_xmm2c8(t0.re,t0.im,
	                             cpx.re,cpx.im,
	                             &B_x.re,&B_x.im);
	                             
	                cmul_xmm2c8(t0.re,t0.im,
	                             cpy.re,cpy.im,
	                             &B_y.re,&B_y.im);
	                            
	                cmul_xmm2c8(t0.re,t0.im,
	                             cpz.re,cpz.im,
	                             &B_z.re,&B_z.im);
	                          
	                                           
	     }
	     
	     
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void B_XYZ_VP_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST); 
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_VP_xmm2c8_unroll6x(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST); 
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_VP_xmm2c8_unroll2x(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST); 
	                                           
	         
	         
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void B_XYZ_VP_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                           const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                           const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                           xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                           const int32_t n,
	                                           int32_t & PF_DIST); 
	         
	         
	                  
	         
	         
	     
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void B_XYZ_VP_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                   const double * __restrict __ATTR_ALIGN__(16) pomega,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvrx,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvry,
	                                   const double * __restrict __ATTR_ALIGN__(16) pvrz,
	                                   const xmm2c8_t k,
	                                   xmm2c8_t & B_x,
	                                   xmm2c8_t & B_y,
	                                   xmm2c8_t & B_z) {
	                         
	                register __m128d vpolx = _mm_load_pd(&pvpolx[0]);
	                register __m128d vpoly = _mm_load_pd(&pvpoly[0]);    
	                register __m128d vpolz = _mm_load_pd(&pvpolz[0]);  
	                register __m128d vdirx = _mm_load_pd(&pvdirx[0]);  
	                register __m128d vdiry = _mm_load_pd(&pvdiry[0]);
	                register __m128d vdirz = _mm_load_pd(&pvdirz[0]); 
	                register __m128d onega = _mm_load_pd(&pomega[0]);
	                register __m128d vrx   = _mm_load_pd(&pvrx[0]);
	                register __m128d vry   = _mm_load_pd(&pvry[0]);
	                register __m128d vrz   = _mm_load_pd(&pvrz[0]);        
	                const __m128d mu0 = _mm_set1_pd(0.0000012566370614359173);
	                xmm2c8_t cdirx;
	                xmm2c8_t cdiry;
	                xmm2c8_t cdirz;
	                xmm2c8_t H_x;
	                xmm2c8_t H_y;
	                xmm2c8_t H_z;
	                xmm2c8_t cpx;
	                xmm2c8_t cpy;
	                xmm2c8_t cpz;
	                xmm2c8_t t0;
	                __m128d zz0;
	                H_XYZ_VP_xmm2c8(vpolx,
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
	                cdirx.im = _mm_setzero_pd();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm_mul_pd(omega,mu0);
	                t0.re    = _mm_div_pd(k.re,zz0);
	                t0.im    = _mm_div_pd(k.im,zz0);
	                scrossc_xmm2c8(cdirx,
	                                cdiry,
	                                cdirz,
	                                H_x,
	                                H_y,
	                                H_z,
	                                cpx,
	                                cpy,
	                                cpz);
	                                
	                cmul_xmm2c8(t0.re,
	                             t0.im,
	                             cpx.re,
	                             cpx.im,
	                             &B_x.re,
	                             &B_x.im);
	                cmul_xmm2c8(t0.re,
	                             t0.im,
	                             cpy.re,
	                             cpy.im,
	                             &B_y.re,
	                             &B_y.im);
	                cmul_xmm2c8(t0.re,
	                             t0.im,
	                             cpz.re,
	                             cpz.im,
	                             &B_z.re,
	                             &B_z.im);
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void B_XYZ_VP_xmm2c8_u(const double * __restrict  pvpolx,
	                                   const double * __restrict  pvpoly,
	                                   const double * __restrict  pvpolz,
	                                   const double * __restrict  pvdirx,
	                                   const double * __restrict  pvdiry,
	                                   const double * __restrict  pvdirz,
	                                   const double * __restrict  pomega,
	                                   const double * __restrict  pvrx,
	                                   const double * __restrict  pvry,
	                                   const double * __restrict  pvrz,
	                                   const xmm2c8_t k,
	                                   xmm2c8_t & B_x,
	                                   xmm2c8_t & B_y,
	                                   xmm2c8_t & B_z) {
	                         
	                register __m128d vpolx = _mm_loadu_pd(&pvpolx[0]);
	                register __m128d vpoly = _mm_loadu_pd(&pvpoly[0]);    
	                register __m128d vpolz = _mm_loadu_pd(&pvpolz[0]);  
	                register __m128d vdirx = _mm_loadu_pd(&pvdirx[0]);  
	                register __m128d vdiry = _mm_loadu_pd(&pvdiry[0]);
	                register __m128d vdirz = _mm_loadu_pd(&pvdirz[0]); 
	                register __m128d onega = _mm_loadu_pd(&pomega[0]);
	                register __m128d vrx   = _mm_loadu_pd(&pvrx[0]);
	                register __m128d vry   = _mm_loadu_pd(&pvry[0]);
	                register __m128d vrz   = _mm_loadu_pd(&pvrz[0]);        
	                const __m128d mu0 = _mm_set1_pd(0.0000012566370614359173);
	                xmm2c8_t cdirx;
	                xmm2c8_t cdiry;
	                xmm2c8_t cdirz;
	                xmm2c8_t H_x;
	                xmm2c8_t H_y;
	                xmm2c8_t H_z;
	                xmm2c8_t cpx;
	                xmm2c8_t cpy;
	                xmm2c8_t cpz;
	                xmm2c8_t t0;
	                __m128d zz0;
	                H_XYZ_VP_xmm2c8(vpolx,
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
	                cdirx.im = _mm_setzero_pd();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm_mul_pd(omega,mu0);
	                t0.re    = _mm_div_pd(k.re,zz0);
	                t0.im    = _mm_div_pd(k.im,zz0);
	                scrossc_xmm2c8(cdirx,
	                                cdiry,
	                                cdirz,
	                                H_x,
	                                H_y,
	                                H_z,
	                                cpx,
	                                cpy,
	                                cpz);
	                                
	                cmul_xmm2c8(t0.re,
	                             t0.im,
	                             cpx.re,
	                             cpx.im,
	                             &B_x.re,
	                             &B_x.im);
	                cmul_xmm2c8(t0.re,
	                             t0.im,
	                             cpy.re,
	                             cpy.im,
	                             &B_y.re,
	                             &B_y.im);
	                cmul_xmm2c8(t0.re,
	                             t0.im,
	                             cpz.re,
	                             cpz.im,
	                             &B_z.re,
	                             &B_z.im);
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void B_XYZ_H_XYZ_P_xmm2c8(const __m128d tht,
	                                      const __m128d phi,
	                                      const __m128d psi,
	                                      const __m128d omg,
	                                      const __m128d px,
	                                      const __m128d py,
	                                      const __m128d pz,
	                                      const xmm2c8_t r,
	                                      xmm2c8_t & H_x,
	                                      xmm2c8_t & H_y,
	                                      xmm2c8_t & H_z,
	                                      xmm2c8_t & B_x,
	                                      xmm2c8_t & B_y,
	                                      xmm2c8_t & B_z) {
	                                      
	                
	                const __m128d c = _mm_set1_pd(299792458.0f);
	                xmm2c8_t k;
	                register __m128d vpolx,vpoly,vpolz;
	                register __m128d vdirx,vdiry,vdirz;
	                register __m128d t0;
	                
	                t0 = _mm_div_pd(omg,c);
	                dir_vec_xmm2c8(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm_mul_pd(r.re,c);
	                k.im = _mm_mul_pd(r.im,c);
	                pol_vec_xmm2c8(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_H_XYZ_P_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	      
	          
	      
	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_H_XYZ_P_xmm2c8_unroll6x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	      	         	      	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                 void B_XYZ_H_XYZ_P_xmm2c8_unroll2x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void B_XYZ_H_XYZ_P_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	      
	       
	       
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void B_XYZ_H_XYZ_P_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) ptht,
	                                      const double * __restrict __ATTR_ALIGN__(16) pphi,
	                                      const double * __restrict __ATTR_ALIGN__(16) ppsi,
	                                      const double * __restrict __ATTR_ALIGN__(16) pomg,
	                                      const double * __restrict __ATTR_ALIGN__(16) ppx,
	                                      const double * __restrict __ATTR_ALIGN__(16) ppy,
	                                      const double * __restrict __ATTR_ALIGN__(16) ppz,
	                                      const xmm2c8_t r,
	                                      xmm2c8_t & H_x,
	                                      xmm2c8_t & H_y,
	                                      xmm2c8_t & H_z,
	                                      xmm2c8_t & B_x,
	                                      xmm2c8_t & B_y,
	                                      xmm2c8_t & B_z) {
	                                      
	                
	                register __m128d tht = _mm_load_pd(&ptht[0]);
	                register __m128d phi = _mm_load_pd(&pphi[0]);
	                register __m128d psi = _mm_load_pd(&ppsi[0]);
	                register __m128d omg = _mm_load_pd(&pomg[0]);
	                register __m128d px  = _mm_load_pd(&ppx[0]);
	                register __m128d py  = _mm_load_pd(&ppy[0]);
	                register __m128d pz  = _mm_load_pd(&ppz[0]);
	                const __m128d c = _mm_set1_pd(299792458.0);
	                xmm2c8_t k;
	                register __m128d vpolx,vpoly,vpolz;
	                register __m128d vdirx,vdiry,vdirz;
	                register __m128d t0;
	                
	                t0 = _mm_div_pd(omg,c);
	                dir_vec_xmm2c8(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm_mul_pd(r.re,c);
	                k.im = _mm_mul_pd(r.im,c);
	                pol_vec_xmm2c8(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void B_XYZ_H_XYZ_P_xmm2c8_u(const double * __restrict ptht,
	                                      const double * __restrict  pphi,
	                                      const double * __restrict ppsi,
	                                      const double * __restrict  pomg,
	                                      const double * __restrict  ppx,
	                                      const double * __restrict  ppy,
	                                      const double * __restrict ppz,
	                                      const xmm2c8_t r,
	                                      xmm2c8_t & H_x,
	                                      xmm2c8_t & H_y,
	                                      xmm2c8_t & H_z,
	                                      xmm2c8_t & B_x,
	                                      xmm2c8_t & B_y,
	                                      xmm2c8_t & B_z) {
	                                      
	                
	                register __m128d tht = _mm_loadu_pd(&ptht[0]);
	                register __m128d phi = _mm_loadu_pd(&pphi[0]);
	                register __m128d psi = _mm_loadu_pd(&ppsi[0]);
	                register __m128d omg = _mm_loadu_pd(&pomg[0]);
	                register __m128d px  = _mm_loadu_pd(&ppx[0]);
	                register __m128d py  = _mm_loadu_pd(&ppy[0]);
	                register __m128d pz  = _mm_loadu_pd(&ppz[0]);
	                const __m128d c = _mm_set1_pd(299792458.0);
	                xmm2c8_t k;
	                register __m128d vpolx,vpoly,vpolz;
	                register __m128d vdirx,vdiry,vdirz;
	                register __m128d t0;
	                
	                t0 = _mm_div_pd(omg,c);
	                dir_vec_xmm2c8(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm_mul_pd(r.re,c);
	                k.im = _mm_mul_pd(r.im,c);
	                pol_vec_xmm2c8(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	       /*
	              ! Electric and Magnetic Fields elliptically polarized
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	          static inline
	           void B_XYZ_H_XYZ_EP_xmm2c8(const __m128d tht,
	                                       const __m128d phi,
	                                       const __m128d omg,
	                                       const xmm2c8_t phase,
	                                       const xmm2c8_t refi,
	                                       const xmm2c8_t px,
	                                       const xmm2c8_t py,
	                                       const xmm2c8_t pz,
	                                       xmm2c8_t & H_x,
	                                       xmm2c8_t & H_y,
	                                       xmm2c8_t & H_z,
	                                       xmm2c8_t & B_x,
	                                       xmm2c8_t & B_y,
	                                       xmm2c8_t & B_z) {
	                                   
	               const __m128d c   = _mm_set1_pd(299792458.0); 
	               const __m128d mu0 = _mm_set1_pd(0.0000012566370614359173);   
	               const __m128d psi0 = _mm_setzero_pd();
	               const __m128d C00  = _mm_setzero_pd();
	               
	               xmm2c8_t H_x_1;
	               xmm2c8_t H_y_1;
	               xmm2c8_t H_z_1;
	               xmm2c8_t H_x_2;
	               xmm2c8_t H_y_2;
	               xmm2c8_t H_z_2;
	               xmm2c8_t k;
	               xmm2c8_t t0;
	               xmm2c8_t cdirx;
	               xmm2c8_t cdiry;
	               xmm2c8_t cdirz;
	               
	               register __m128d vpolx;
	               register __m128d vpoly;
	               register __m128d vpolz;
	               register __m128d vdirx;
	               register __m128d vdiry;
	               register __m128d vdirz;
	               register __m128d cn;
	               register __m128d x0;
	               register __m128d t0r,t0i;
	               register __m128d t1r,t1i;
	               register __m128d t2r,t2i;
	               
	               dir_vec_xmm2c8(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm_div_pd(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm_mul_pd(refi.re,x0);
	               k.im     = _mm_mul_pd(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_xmm2c8(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_xmm2c8( vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm_add_pd(H_x_1.re,t0r);
	               H_x.im = _mm_add_pd(H_x_1.im,t0i);
	               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm_add_pd(H_y_1.re,t1r);
	               H_y.im = _mm_add_pd(H_y_1.im,t1i);
	               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm_add_pd(H_z_2.re,t2r);
	              H_z.im = _mm_add_pd(H_z_2.im,t2i);
	              
	              cn = cnorm_xmm2c8(H_x,H_y,H_z);
	              
	              x0     = _mm_div_pd(omg,mu0);
	              H_x.re = _mm_div_pd(H_x.re,cn);
	              H_x.im = _mm_div_pd(H_x.im,cn);
	              H_y.re = _mm_div_pd(H_y.re,cn);
	              H_y.im = _mm_div_pd(H_y.im,cn);
	              H_z.re = _mm_div_pd(H_z.re,cn);
	              H_z.im = _mm_div_pd(H_z.im,cn); 
	              
	              t0.re  = _mm_div_pd(k.re,x0);
	              t0.im  = _mm_div_pd(k.im,x0);
	              
	              scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	        
	        
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_H_XYZ_EP_xmm2c8_unroll10x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST); 
	      
	           
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                   void B_XYZ_H_XYZ_EP_xmm2c8_unroll6x(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST); 
	                                                 
	                    	          	      	      
	      
	   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_H_XYZ_EP_xmm2c8_unroll2x(  const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST); 
	      
	      
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(16)
                  void B_XYZ_H_XYZ_EP_xmm2c8_rolled(    const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                 const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                 const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                 xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                 const int32_t n,
	                                                 int32_t & PF_DIST);
	      	      	      	      
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void B_XYZ_H_XYZ_EP_xmm2c8_a(const double * __restrict __ATTR_ALIGN__(16) ptht,
	                                         const double * __restrict __ATTR_ALIGN__(16) pphi,
	                                         const double * __restrict __ATTR_ALIGN__(16) pomg,
	                                         const xmm2c8_t phase,
	                                         const xmm2c8_t refi,
	                                         const xmm2c8_t px,
	                                         const xmm2c8_t py,
	                                         const xmm2c8_t pz,
	                                         xmm2c8_t & H_x,
	                                         xmm2c8_t & H_y,
	                                         xmm2c8_t & H_z,
	                                         xmm2c8_t & B_x,
	                                         xmm2c8_t & B_y,
	                                         xmm2c8_t & B_z) {
	                         
	               register __m128d tht = _mm_load_pd(&ptht[0]);
	               register __m128d phi = _mm_load_pd(&pphi[0]);
	               register __m128d omg = _mm_load_pd(&pomg[0]);            
	               const __m128d c   = _mm_set1_pd(299792458.0); 
	               const __m128d mu0 = _mm_set1_pd(0.0000012566370614359173);   
	               const __m128d psi0 = _mm_setzero_pd();
	               const __m128d C00  = _mm_setzero_pd();
	               
	               xmm2c8_t H_x_1;
	               xmm2c8_t H_y_1;
	               xmm2c8_t H_z_1;
	               xmm2c8_t H_x_2;
	               xmm2c8_t H_y_2;
	               xmm2c8_t H_z_2;
	               xmm2c8_t k;
	               xmm2c8_t t0;
	               xmm2c8_t cdirx;
	               xmm2c8_t cdiry;
	               xmm2c8_t cdirz;
	               
	               register __m128d vpolx;
	               register __m128d vpoly;
	               register __m128d vpolz;
	               register __m128d vdirx;
	               register __m128d vdiry;
	               register __m128d vdirz;
	               register __m128d cn;
	               register __m128d x0;
	               register __m128d t0r,t0i;
	               register __m128d t1r,t1i;
	               register __m128d t2r,t2i;
	               
	               dir_vec_xmm2c8(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm_div_pd(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm_mul_pd(refi.re,x0);
	               k.im     = _mm_mul_pd(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_xmm2c8(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm_add_pd(H_x_1.re,t0r);
	               H_x.im = _mm_add_pd(H_x_1.im,t0i);
	               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm_add_pd(H_y_1.re,t1r);
	               H_y.im = _mm_add_pd(H_y_1.im,t1i);
	               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm_add_pd(H_z_2.re,t2r);
	              H_z.im = _mm_add_pd(H_z_2.im,t2i);
	              
	              cn = cnorm_xmm2c8(H_x,H_y,H_z);
	              
	              x0     = _mm_div_pd(omg,mu0);
	              H_x.re = _mm_div_pd(H_x.re,cn);
	              H_x.im = _mm_div_pd(H_x.im,cn);
	              H_y.re = _mm_div_pd(H_y.re,cn);
	              H_y.im = _mm_div_pd(H_y.im,cn);
	              H_z.re = _mm_div_pd(H_z.re,cn);
	              H_z.im = _mm_div_pd(H_z.im,cn); 
	              
	              t0.re  = _mm_div_pd(k.re,x0);
	              t0.im  = _mm_div_pd(k.im,x0);
	              
	              scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           static inline
	           void B_XYZ_H_XYZ_EP_xmm2c8_u(const double * __restrict  ptht,
	                                         const double * __restrict  pphi,
	                                         const double * __restrict  pomg,
	                                         const xmm2c8_t phase,
	                                         const xmm2c8_t refi,
	                                         const xmm2c8_t px,
	                                         const xmm2c8_t py,
	                                         const xmm2c8_t pz,
	                                         xmm2c8_t & H_x,
	                                         xmm2c8_t & H_y,
	                                         xmm2c8_t & H_z,
	                                         xmm2c8_t & B_x,
	                                         xmm2c8_t & B_y,
	                                         xmm2c8_t & B_z) {
	                         
	               register __m128d tht = _mm_loadu_pd(&ptht[0]);
	               register __m128d phi = _mm_loadu_pd(&pphi[0]);
	               register __m128d omg = _mm_loadu_pd(&pomg[0]);            
	               const __m128d c   = _mm_set1_pd(299792458.0); 
	               const __m128d mu0 = _mm_set1_pd(0.0000012566370614359173);   
	               const __m128d psi0 = _mm_setzero_pd();
	               const __m128d C00  = _mm_setzero_pd();
	               
	               xmm2c8_t H_x_1;
	               xmm2c8_t H_y_1;
	               xmm2c8_t H_z_1;
	               xmm2c8_t H_x_2;
	               xmm2c8_t H_y_2;
	               xmm2c8_t H_z_2;
	               xmm2c8_t k;
	               xmm2c8_t t0;
	               xmm2c8_t cdirx;
	               xmm2c8_t cdiry;
	               xmm2c8_t cdirz;
	               
	               register __m128d vpolx;
	               register __m128d vpoly;
	               register __m128d vpolz;
	               register __m128d vdirx;
	               register __m128d vdiry;
	               register __m128d vdirz;
	               register __m128d cn;
	               register __m128d x0;
	               register __m128d t0r,t0i;
	               register __m128d t1r,t1i;
	               register __m128d t2r,t2i;
	               
	               dir_vec_xmm2c8(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm_div_pd(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm_mul_pd(refi.re,x0);
	               k.im     = _mm_mul_pd(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_xmm2c8(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_xmm2c8(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm_add_pd(H_x_1.re,t0r);
	               H_x.im = _mm_add_pd(H_x_1.im,t0i);
	               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm_add_pd(H_y_1.re,t1r);
	               H_y.im = _mm_add_pd(H_y_1.im,t1i);
	               
	               cmul_xmm2c8(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm_add_pd(H_z_2.re,t2r);
	              H_z.im = _mm_add_pd(H_z_2.im,t2i);
	              
	              cn = cnorm_xmm2c8(H_x,H_y,H_z);
	              
	              x0     = _mm_div_pd(omg,mu0);
	              H_x.re = _mm_div_pd(H_x.re,cn);
	              H_x.im = _mm_div_pd(H_x.im,cn);
	              H_y.re = _mm_div_pd(H_y.re,cn);
	              H_y.im = _mm_div_pd(H_y.im,cn);
	              H_z.re = _mm_div_pd(H_z.re,cn);
	              H_z.im = _mm_div_pd(H_z.im,cn); 
	              
	              t0.re  = _mm_div_pd(k.re,x0);
	              t0.im  = _mm_div_pd(k.im,x0);
	              
	              scrossc_xmm2c8(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_xmm2c8(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	                                      
	       
	       
	       
	     
                
                
        } // radiolocation

} // gms


#endif /*__GMS_EM_FIELDS_XMM2R8_H__*/
