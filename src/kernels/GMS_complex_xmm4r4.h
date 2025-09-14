
#ifndef __GMS_COMPLEX_XMM4R4_H__
#define __GMS_COMPLEX_XMM4R4_H__ 181020231557


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

namespace file_version 
{

    const unsigned int GMS_COMPLEX_XMM4R4_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_XMM4R4_MINOR = 0U;
    const unsigned int GMS_COMPLEX_XMM4R4_MICRO = 0U;
    const unsigned int GMS_COMPLEX_XMM4R4_FULLVER =
      1000U*GMS_COMPLEX_XMM4R4_MAJOR+
      100U*GMS_COMPLEX_XMM4R4_MINOR+
      10U*GMS_COMPLEX_XMM4R4_MICRO;
    const char * const GMS_COMPLEX_XMM4R4_CREATION_DATE = "22-10-2023 10:35 AM +00200 (SAT 22 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_XMM4R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_XMM4R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_XMM4R4_DESCRIPTION   = "SSE optimized complex number implementation.";

}

#include <cstdint>
#include <immintrin.h>
#include <utility>
#include "GMS_config.h"
#include "GMS_simd_utils.h"

namespace  gms 
{


namespace math 
{
       
       
                   struct alignas(16) xmm4c4_t 
                   {
                          __m128 re;
                          __m128 im;
                   };
                   
                   
      	          
             
                                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   xmm4c4_t cadd_xmm4r4(const xmm4c4_t x,
                                        const xmm4c4_t y) {
                                     
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0   = _mm_add_ps(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1   = _mm_add_ps(x.im,y.im);
                        cv.im  = xmm1;  
                        return (cv);            
               }


                  
		          
          	                   	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                 
	           static inline
	           xmm4c4_t cadd_xmm4r4(const xmm4c4_t x,
                                     const __m128 s) {
                      
                      xmm4c4_t cv;
                      cv.re =  _mm_add_ps(x.re,s);
                      cv.im =  x.im;
                      return (cv);                       
                }


                ////////////////////////////////////////////////////////////////////

               	          
	               
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   xmm4c4_t csub_xmm4r4(const xmm4c4_t x,
                                          const xmm4c4_t y) {
                                    
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0  = _mm_sub_ps(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1  = _mm_sub_ps(x.im,y.im);
                        cv.im  = xmm1;
                        return (cv);
                }
                
                
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   xmm4c4_t csub_xmm4r4(const xmm4c4_t x,
                                          const __m128 s) {
                                    
                        xmm4c4_t cv;
                        cv.re = _mm_sub_ps(x.re,s);
                        cv.im = x.im;
                        return (cv);
               }
               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                   xmm4c4_t cmul_xmm4r4(const xmm4c4_t x,
                                        const xmm4c4_t y) {
                                     
                         xmm4c4_t cv;
                         register __m128 xmm0,xmm1;
                         xmm0 = _mm_sub_ps(_mm_mul_ps(x.re,y.re),
                                              _mm_mul_ps(x.im,y.im));
                         cv.re  = xmm0;
                         xmm1 = _mm_mul_ps(_mm_mul_ps(x.im,y.re),
                                              _mm_mul_ps(x.re,y.im));
                         cv.im  = xmm1;
                         return (cv);
                }


                 
	        
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__ 	                            
	           static inline
                   xmm4c4_t cmul_xmm4r4(const xmm4c4_t x,
                                          const __m128 s) {
                                     
                        xmm4c4_t cv;
                        cv.re = _mm_mul_ps(x.re,s);
                        cv.im = _mm_mul_ps(x.im,s);
                        return (cv);
               }
               


               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                   xmm4c4_t cdiv_xmm4r4(const xmm4c4_t x,
                                          const xmm4c4_t y) {
                                     
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_ps(x.re,y.re,
                                           _mm_mul_ps(x.im,y.im));
                      xmm1 = _mm_fmsub_ps(x.im,y.re,
                                           _mm_mul_ps(x.re,y.im));
                      xmm2 = _mm_fmadd_ps(xmm0,xmm0,
                                           _mm_mul_ps(xmm1,xmm1));
                      cv.re  = _mm_div_ps(xmm0,xmm2);
                      cv.im  = _mm_div_ps(xmm1,xmm2);
                      return (cv);
                }


          	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                xmm4c4_t cdiv_xmm4r4(const xmm4c4_t x,
                                     const __m128 s) {
                                     
                         xmm4c4_t cv;
                         cv.re = _mm_div_ps(x.re,s);
                         cv.im = _mm_div_ps(x.im,s);
                         return (cv);
               }
               
               
               
          	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                   xmm4c4_t cdiv_xmm4r4_s(const __m128 s,
                                          const xmm4c4_t x) {
                                       
                        xmm4c4_t cv,tmp;
                        tmp.re = s;
                        tmp.im = _mm_setzero_ps();
                        cv = cdiv_xmm4r4(tmp,x);
                        return (cv);                     
                 }
                 
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   xmm4c4_t cdiv_smith_xmm4r4(const xmm4c4_t x,
                                                const xmm4c4_t y) {
                                           
                        xmm4c4_t cv;
                        register __m128 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_ps_mask(xmm4r4_abs(y.re),
                                                  xmm4r4_abs(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_ps(m,_mm_div_ps(y.re,y.im),
                                                      _mm_div_ps(y.im,y.re)); // r
                        den  = _mm_mask_blend_ps(m,_mm_fmadd_ps(r,y.re,y.im),
                                                      _mm_fmadd_ps(r,y.im,y.re));
                        cv.re  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmadd_ps(x.re,r,x.im),den),
                                                _mm_div_ps(_mm_fmadd_ps(x.im,r,x.re),den));
                        cv.im  =  _mm_mask_blend_ps(m,
                                                _mm_div_ps(_mm_fmsub_ps(x.im,r,x.re),den),
                                                _mm_div_ps(_mm_sub_ps(x.im,_mm_mul_ps(x.re,r)),den));
                        return (cv);
               }          
	        
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                xmm4c4_t cdiv_smith_xmm4r4_s(const __m128 s,
                                                const xmm4c4_t x) {
                                             
                        xmm4c4_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm_setzero_ps(); 
                        cv = cdiv_smith_xmm4r4(t0,x);
                        return (cv);                 
                 }
                 
                   
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m128 cabs_xmm4r4(const __m128 re,
                                       const __m128 im) {

                        register __m128 xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_ps(re,re);
                        xmm1 = _mm_mul_ps(im,im);
                        cabs = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        return (cabs);
                 }
                 
                 
                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m128 cabs_xmm4r4(const xmm4c4_t x) {

                        register __m128 xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_ps(x.re,x.re);
                        xmm1 = _mm_mul_ps(x.im,x.im);
                        cabs = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        return (cabs);
                 }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m128 carg_xmm4r4(const __m128 re,
                                       const __m128 im) {

                       register __m128 carg;
                       carg = _mm_atan2_ps(re,im);
                       return (carg);
                }
                
                
                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m128 carg_xmm4r4(xmm4c4_t x) {

                       register __m128 carg;
                       carg = _mm_atan2_ps(x.re,x.im);
                       return (carg);
                }
                
                
                 
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
	           xmm4c4_t clog_xmm4r4(const xmm4c4_t x){
	                                  
	                xmm4c4_t clog;                           
	                register __m128 t1,t2,ln;
	                t1  = cabs_xmm4r4(x.re,x.im);
	                t2  = carg_xmm4r4(x.re,x.im);
	                ln  = _mm_log_ps(t1);
	                clog.re = ln;
	                clog.im = t2;
	                return (clog);                    
	        }

                 

                  	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t cconj_xmm4r4_v2(const __m128 xre,
                                              const __m128 xim) {                                              
                         
                        //register __m128 c;              
                        //c = negate_xmm4r4(*im);
                        //*im = c;
                        xmm4c4_t cv;
                        cv.re = xre; 
                        cv.im = negate_xmm4r4(xim);
                        return (cv);
                   } 
                   
                   
                
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t cconj_xmm4r4_v2(const xmm4c4_t x) {                                              
                         
                        //register __m128 c;              
                        //c = negate_xmm4r4(*im);
                        //*im = c;
                        xmm4c4_t cv;
                        cv.re = x.re; 
                        cv.im = negate_xmm4r4(x.im);
                        return (cv);
                   } 
                   


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t ccos_xmm4r4(const __m128 xre,
                                          const __m128 xim) {
                                    
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cos_ps(xre),_mm_cosh_ps(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sin_ps(xre),_mm_sinh_ps(xim));
                      cv.im = xmm1;
                      return (cv); 
               }
               
               
                 	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t ccos_xmm4r4(const xmm4c4_t x) {
                                    
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cos_ps(x.re),_mm_cosh_ps(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sin_ps(x.re),_mm_sinh_ps(x.im));
                      cv.im = xmm1;
                      return (cv); 
               }


                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t ccosh_xmm4r4(const __m128 xre,
                                           const __m128 xim) {
                                          
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cosh_ps(xre),_mm_cos_ps(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sinh_ps(xre),_mm_sin_ps(xim));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   
	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t ccosh_xmm4r4(const xmm4c4_t x) {
                                          
                      xmm4c4_t cv;
                      register __m128 xmm0,xmm1;
                      xmm0  = _mm_mul_ps(_mm_cosh_ps(x.re),_mm_cos_ps(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_ps(_mm_sinh_ps(x.re),_mm_sin_ps(x.im));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
	           xmm4c4_t cpow_xmm4r4(const xmm4c4_t x,
	                                  const float n) {
	                   
	                xmm4c4_t cp;        
	                register __m128 xmm0,xmm1;
	                register __m128 r,tht;
	                register __m128 vn,pt;
	                register __m128 ta;
	                xmm0  = _mm_mul_ps(x.re,x.re);
	                vn    = _mm_set1_ps(n);
	                xmm1  = _mm_mul_ps(x.im,x.im);
	                r     = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
	                tht   = _mm_atan_ps(_mm_div_ps(x.im,x.re));
	                pt    = _mm_pow_ps(r,vn);
	                ta    = _mm_mul_ps(vn,tht);
	                cp.re = _mm_mul_ps(pt,_mm_cos_ps(ta));
	                cp.im = _mm_mul_ps(pt,_mm_sin_ps(ta));      
	                return (cp);              
	       }
	       


                
              
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_xmm4r4(      const xmm4c4_t x,
                                    const xmm4c4_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              


                 	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_xmm4r4(      const xmm4c4_t x,
                                    const xmm4c4_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_xmm4r4(     const xmm4c4_t x,
                                    const xmm4c4_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }



                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                       
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_xmm4r4(      const xmm4c4_t x,
                                     const xmm4c4_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm_cmp_ps_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm_cmp_ps_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm4c4_t cexp_xmm4r4(const __m128 xre,
                                          const __m128 xim) {
                                     
                        xmm4c4_t cv;
                        register const __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0;
                        xmm0   = _mm_exp_ps(xre);
                        cv.re = _mm_mul_ps(xmm0,_mm_cos_ps(xre));
                        cv.im = _mm_mul_ps(xmm0,_mm_mul_ps(_mm_sin_ps(xim),I));
                        return (cv);
              }
              
              
                  
	        
	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                       
	           static inline
                   xmm4c4_t cexp_xmm4r4(const xmm4c4_t x) {
                                     
                        xmm4c4_t cv;
                        register  __m128 I = _mm_set1_ps(1.0);
                        register __m128 xmm0;
                        xmm0   = _mm_exp_ps(x.re);
                        cv.re = _mm_mul_ps(xmm0,_mm_cos_ps(x.re));
                        cv.im = _mm_mul_ps(xmm0,_mm_mul_ps(_mm_sin_ps(x.im),I));
                        return (cv);
              }


                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm4c4_t cpolar_xmm4r4(const __m128 rho,
                                            const __m128 tht) {
                                      
                        xmm4c4_t cv;
                        register __m128 xmm0,xmm1;
                        xmm0 = _mm_mul_ps(rho,_mm_cos_ps(tht));
                        cv.re  = xmm0;
                        xmm1 = _mm_mul_ps(rho,_mm_sin_ps(tht));
                        cv.im  = xmm1;
                        return (cv);
              }
              


     	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                       
	           static inline
                   xmm4c4_t csqrt_xmm4r4(const __m128 xre,
                                           const __m128 xim,
                                          __m128 * __restrict wrkc) {
                                          
                       xmm4c4_t cv;
                       register __m128 xmm0,xmm1;
                       register __m128 half = _mm_set1_ps(0.5); 
                       cabs_xmm4r4(xre,xim);
                       xmm0  = _mm_mul_ps(half,_mm_add_ps(*wrkc,xre));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_ps(half,_mm_sub_ps(*wrkc,xre));
                       cv.im = xmm1; 
                       return (cv);
              }
              
              
                   
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm4c4_t csqrt_xmm4r4(const xmm4c4_t x,
                                          __m128 * __restrict wrkc) {
                                          
                       xmm4c4_t cv;
                       register __m128 xmm0,xmm1;
                       register __m128 half = _mm_set1_ps(0.5); 
                       cabs_xmm4r4(x.re,x.im);
                       xmm0  = _mm_mul_ps(half,_mm_add_ps(*wrkc,x.re));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_ps(half,_mm_sub_ps(*wrkc,x.re));
                       cv.im = xmm1; 
                       return (cv);
              }


                  
	        
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm4c4_t cnorm_prod_xmm4r4(  const __m128  xre,
                                                  const __m128  xim,
                                                  const __m128  yre,
                                                  const __m128  yim) {
                                             
                        xmm4c4_t cv;
                        register __m128 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        imp  = _mm_fmadd_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        cv.re = _mm_div_ps(rep,xmm2);
                        cv.im = _mm_div_ps(imp,xmm2);
                        return (cv);
             }
             
             
                 
	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm4c4_t cnorm_prod_xmm4r4(  const xmm4c4_t x,
                                                  const xmm4c4_t y) {
                                             
                        xmm4c4_t cv;
                        register __m128 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        imp  = _mm_fmadd_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        xmm0 = _mm_mul_ps(rep,rep);
                        xmm1 = _mm_mul_ps(imp,imp);
                        xmm2 = _mm_sqrt_ps(_mm_add_ps(xmm0,xmm1));
                        cv.re = _mm_div_ps(rep,xmm2);
                        cv.im = _mm_div_ps(imp,xmm2);
                        return (cv);
             }
             
#include "GMS_simd_utils.h"


                  	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   void cmean_prod_xmm4r4(const xmm4c4_t x,
                                           const xmm4c4_t y,
                                           float * __restrict mre,
                                           float * __restrict mim) {

                        register __m128 rep,imp;
                        __m128 t0,t1;
                        constexpr float inv16 = 0.25f;
                        float sre,sim;
                        sre = 0.0;
                        rep  = _mm_fmsub_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        t0   = _mm_hadd_ps(rep,rep);
                        sre  = *(float*)&t0;
                        *mre = sre*inv16;
                        sim  = 0.0;
                        imp  = _mm_fmadd_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        t1   = _mm_hadd_ps(imp,imp);
                        sim  = *(float*)&t1;
                        *mim = sim*inv16;
             }


               	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   void cmean_quot_xmm4r4(  const xmm4c4_t x,
                                             const xmm4c4_t y,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m128 rep,imp,den,rquot,iquot;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.25f;
                        float sre,sim;
                        sre  = 0.0;
                        rep  = _mm_fmsub_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        imp  = _mm_fmadd_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        sim  = 0.0;
                        den  = _mm_fmadd_ps(y.re,y.re,
                                               _mm_mul_ps(y.im,y.im));
                        rquot = _mm_div_ps(rep,den);
                        t0    = _mm_hadd_ps(rquot,rquot);
                        sre   = *(float*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_ps(imp,den);
                        t1    = _mm_hadd_ps(iquot,iquot);
                        sim   = *(float*)&t1;
                        *mim  = sre*inv2;
              }  


               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cnorm_cprod_xmm4r4(const xmm4c4_t x,
                                            const xmm4c4_t y,
                                            __m128 * __restrict mre,
                                            __m128 * __restrict mim) {

                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmadd_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        *mre = _mm_div_ps(rep,vcmag);
                        *mim = _mm_div_ps(imp,vcmag);
             }
             
             
               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t cnorm_cprod_xmm4r4(const __m128 xre,
                                                 const __m128 xim,
                                                 const __m128 yre,
                                                 const __m128 yim) {
                                               
                        xmm4c4_t cv;
                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmadd_ps(xre,yre,
                                               _mm_mul_ps(xim,yim));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(xim,yre,
                                               _mm_mul_ps(xre,yim));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        cv.re = _mm_div_ps(rep,vcmag);
                        cv.im = _mm_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t cnorm_cprod_xmm4r4(const xmm4c4_t x,
                                                 const xmm4c4_t y) {
                                               
                        xmm4c4_t cv;
                        register __m128 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmadd_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        magc1= _mm_mul_ps(rep,rep);
                        imp  = _mm_fmsub_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        magc2= _mm_mul_ps(imp,imp);
                        vcmag= _mm_sqrt_ps(_mm_add_ps(magc1,magc2));
                        cv.re = _mm_div_ps(rep,vcmag);
                        cv.im = _mm_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
             
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cmean_cprod_xmm4r4(const xmm4c4_t x,
                                            const xmm4c4_t y,
                                            float * __restrict mre,
                                            float * __restrict mim) {

                        register __m128 re,im;
                        __m128 t0,t1;
                        constexpr float inv2 = 0.25f;
                        float sre,sim;
                        re   = _mm_fmadd_ps(x.re,y.re,
                                               _mm_mul_ps(x.im,y.im));
                        t0   = _mm_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_ps(x.im,y.re,
                                               _mm_mul_ps(x.re,y.im));
                        t1   = _mm_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }


                 
	          
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void arith_cmean_xmm4r4(  const __m128 xre,
                                              const __m128 xim,
                                              float * __restrict mre,
                                              float * __restrict mim) {

                        constexpr float inv2 = 0.25f;
                        __m128 t0,t1;
                        float sre,sim;
                        t0   = _mm_hadd_ps(xre,xre);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_ps(xim,xim);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }
             
             
                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void arith_cmean_xmm4r4(  const xmm4c4_t x,
                                              float * __restrict mre,
                                              float * __restrict mim) {

                        constexpr float inv2 = 0.25f;
                        __m128 t0,t1;
                        float sre,sim;
                        t0   = _mm_hadd_ps(x.re,x.re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_ps(x.im,x.im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }


                       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cnormalize_xmm4r4( const xmm4c4_t x,
                                            const xmm4c4_t y,
                                            __m128 * __restrict mre,
                                            __m128 * __restrict mim ) {

                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(x.re,y.re,
                                                    _mm_mul_ps(x.im,y.im)));
                        *mre = _mm_div_ps(x.re,cvmag);
                        *mim =  _mm_div_ps(x.im,cvmag);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm4c4_t cnormalize_xmm4r4( const __m128 xre,
                                                 const __m128 xim,
                                                 const __m128 yre,
                                                 const __m128 yim) {
                                            
                        xmm4c4_t cv;
                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xre,yre,
                                                    _mm_mul_ps(xim,yim)));
                        cv.re = _mm_div_ps(xre,cvmag);
                        cv.im =  _mm_div_ps(xim,cvmag);
                        return (cv);
             }
             
             
                  	        
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif	          
                  
	           static inline
                   xmm4c4_t cnormalize_xmm4r4( const xmm4c4_t x,
                                                 const xmm4c4_t y) {
                                            
                        xmm4c4_t cv;
                        register __m128 re,im,cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(x.re,y.re,
                                                    _mm_mul_ps(x.im,y.im)));
                        cv.re = _mm_div_ps(x.re,cvmag);
                        cv.im =  _mm_div_ps(x.im,cvmag);
                        return (cv);
             }


               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cmagnitude_xmm4r4(   const __m128 xre,
                                              const __m128 xim,
                                              const __m128 yre,
                                              const __m128 yim,
                                              __m128 * __restrict  mre) {

                        register __m128 cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(xre,yre,
                                                          _mm_mul_ps(xim,yim)));
                        *mre = cvmag;
             }
             
             
               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cmagnitude_xmm4r4(   const xmm4c4_t x,
                                              const xmm4c4_t y,
                                              __m128 * __restrict  mre) {

                        register __m128 cvmag;
                        cvmag= _mm_sqrt_ps(_mm_fmadd_ps(x.re,y.re,
                                          _mm_mul_ps(x.im,y.im)));
                        *mre = cvmag;
             }







} // math


} // gms















#endif /*__GMS_COMPLEX_XMM2R8_H__*/
