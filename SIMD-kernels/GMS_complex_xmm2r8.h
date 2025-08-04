
#ifndef __GMS_COMPLEX_XMM2R8_HPP__
#define __GMS_COMPLEX_XMM2R8_HPP__ 181020231557


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

    const unsigned int GMS_COMPLEX_XMM2R8_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_XMM2R8_MINOR = 0U;
    const unsigned int GMS_COMPLEX_XMM2R8_MICRO = 0U;
    const unsigned int GMS_COMPLEX_XMM2R8_FULLVER =
      1000U*GMS_COMPLEX_XMM2R8_MAJOR+
      100U*GMS_COMPLEX_XMM2R8_MINOR+
      10U*GMS_COMPLEX_XMM2R8_MICRO;
    const char * const GMS_COMPLEX_XMM2R8_CREATION_DATE = "22-10-2023 10:35 AM +00200 (SAT 22 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_XMM2R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_XMM2R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_XMM2R8_DESCRIPTION   = "SSE optimized complex number implementation.";

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
       
       
                   struct alignas(16) xmm2c8_t 
                   {
                          __m128d re;
                          __m128d im;
                   };
                   
                   
      	          
             
                                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   xmm2c8_t cadd_xmm2c8(const xmm2c8_t x,
                                        const xmm2c8_t y) {
                                     
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0   = _mm_add_pd(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1   = _mm_add_pd(x.im,y.im);
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
	           xmm2c8_t cadd_xmm2c8(const xmm2c8_t x,
                                     const __m128d s) {
                      
                      xmm2c8_t cv;
                      cv.re =  _mm_add_pd(x.re,s);
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
                   xmm2c8_t csub_xmm2c8(const xmm2c8_t x,
                                          const xmm2c8_t y) {
                                    
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0  = _mm_sub_pd(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1  = _mm_sub_pd(x.im,y.im);
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
                   xmm2c8_t csub_xmm2c8(const xmm2c8_t x,
                                          const __m128d s) {
                                    
                        xmm2c8_t cv;
                        cv.re = _mm_sub_pd(x.re,s);
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
                   xmm2c8_t cmul_xmm2c8(const xmm2c8_t x,
                                        const xmm2c8_t y) {
                                     
                         xmm2c8_t cv;
                         register __m128d xmm0,xmm1;
                         xmm0 = _mm_sub_pd(_mm_mul_pd(x.re,y.re),
                                              _mm_mul_pd(x.im,y.im));
                         cv.re  = xmm0;
                         xmm1 = _mm_mul_pd(_mm_mul_pd(x.im,y.re),
                                              _mm_mul_pd(x.re,y.im));
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
                   xmm2c8_t cmul_xmm2c8(const xmm2c8_t x,
                                          const __m128d s) {
                                     
                        xmm2c8_t cv;
                        cv.re = _mm_mul_pd(x.re,s);
                        cv.im = _mm_mul_pd(x.im,s);
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
                   xmm2c8_t cdiv_xmm2c8(const xmm2c8_t x,
                                          const xmm2c8_t y) {
                                     
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1,xmm2;
                      xmm0 = _mm_fmadd_pd(x.re,y.re,
                                           _mm_mul_pd(x.im,y.im));
                      xmm1 = _mm_fmsub_pd(x.im,y.re,
                                           _mm_mul_pd(x.re,y.im));
                      xmm2 = _mm_fmadd_pd(xmm0,xmm0,
                                           _mm_mul_pd(xmm1,xmm1));
                      cv.re  = _mm_div_pd(xmm0,xmm2);
                      cv.im  = _mm_div_pd(xmm1,xmm2);
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
                xmm2c8_t cdiv_xmm2c8(const xmm2c8_t x,
                                     const __m128d s) {
                                     
                         xmm2c8_t cv;
                         cv.re = _mm_div_pd(x.re,s);
                         cv.im = _mm_div_pd(x.im,s);
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
                   xmm2c8_t cdiv_xmm2c8_s(const __m128d s,
                                          const xmm2c8_t x) {
                                       
                        xmm2c8_t cv,tmp;
                        tmp.re = s;
                        tmp.im = _mm_setzero_pd();
                        cv = cdiv_xmm2c8(tmp,x);
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
                   xmm2c8_t cdiv_smith_xmm2c8(const xmm2c8_t x,
                                                const xmm2c8_t y) {
                                           
                        xmm2c8_t cv;
                        register __m128d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm_cmp_pd_mask(xmm2r8_abs(y.re),
                                                  xmm2r8_abs(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm_mask_blend_pd(m,_mm_div_pd(y.re,y.im),
                                                      _mm_div_pd(y.im,y.re)); // r
                        den  = _mm_mask_blend_pd(m,_mm_fmadd_pd(r,y.re,y.im),
                                                      _mm_fmadd_pd(r,y.im,y.re));
                        cv.re  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmadd_pd(x.re,r,x.im),den),
                                                _mm_div_pd(_mm_fmadd_pd(x.im,r,x.re),den));
                        cv.im  =  _mm_mask_blend_pd(m,
                                                _mm_div_pd(_mm_fmsub_pd(x.im,r,x.re),den),
                                                _mm_div_pd(_mm_sub_pd(x.im,_mm_mul_pd(x.re,r)),den));
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
                xmm2c8_t cdiv_smith_xmm2c8_s(const __m128d s,
                                                const xmm2c8_t x) {
                                             
                        xmm2c8_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm_setzero_pd(); 
                        cv = cdiv_smith_xmm2c8(t0,x);
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
                   __m128d cabs_xmm2c8(const __m128d re,
                                       const __m128d im) {

                        register __m128d xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_pd(re,re);
                        xmm1 = _mm_mul_pd(im,im);
                        cabs = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
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
                   __m128d cabs_xmm2c8(const xmm2c8_t x) {

                        register __m128d xmm0,xmm1,cabs;
                        xmm0 = _mm_mul_pd(x.re,x.re);
                        xmm1 = _mm_mul_pd(x.im,x.im);
                        cabs = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
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
                   __m128d carg_xmm2c8(const __m128d re,
                                       const __m128d im) {

                       register __m128d carg;
                       carg = _mm_atan2_pd(re,im);
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
                   __m128d carg_xmm2c8(xmm2c8_t x) {

                       register __m128d carg;
                       carg = _mm_atan2_pd(x.re,x.im);
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
	           xmm2c8_t clog_xmm2c8(const xmm2c8_t x){
	                                  
	                xmm2c8_t clog;                           
	                register __m128d t1,t2,ln;
	                t1  = cabs_xmm2c8(x.re,x.im);
	                t2  = carg_xmm2c8(x.re,x.im);
	                ln  = _mm_log_pd(t1);
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
                   xmm2c8_t cconj_xmm2c8_v2(const __m128d xre,
                                              const __m128d xim) {                                              
                         
                        //register __m128d c;              
                        //c = negate_xmm2c8(*im);
                        //*im = c;
                        xmm2c8_t cv;
                        cv.re = xre; 
                        cv.im = negate_xmm2r8(xim);
                        return (cv);
                   } 
                   
                   
                
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm2c8_t cconj_xmm2c8_v2(const xmm2c8_t x) {                                              
                         
                        //register __m128d c;              
                        //c = negate_xmm2c8(*im);
                        //*im = c;
                        xmm2c8_t cv;
                        cv.re = x.re; 
                        cv.im = negate_xmm2r8(x.im);
                        return (cv);
                   } 
                   


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm2c8_t ccos_xmm2c8(const __m128d xre,
                                          const __m128d xim) {
                                    
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cos_pd(xre),_mm_cosh_pd(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sin_pd(xre),_mm_sinh_pd(xim));
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
                   xmm2c8_t ccos_xmm2c8(const xmm2c8_t x) {
                                    
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cos_pd(x.re),_mm_cosh_pd(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sin_pd(x.re),_mm_sinh_pd(x.im));
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
                   xmm2c8_t ccosh_xmm2c8(const __m128d xre,
                                           const __m128d xim) {
                                          
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cosh_pd(xre),_mm_cos_pd(xim));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sinh_pd(xre),_mm_sin_pd(xim));
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
                   xmm2c8_t ccosh_xmm2c8(const xmm2c8_t x) {
                                          
                      xmm2c8_t cv;
                      register __m128d xmm0,xmm1;
                      xmm0  = _mm_mul_pd(_mm_cosh_pd(x.re),_mm_cos_pd(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm_mul_pd(_mm_sinh_pd(x.re),_mm_sin_pd(x.im));
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
	           xmm2c8_t cpow_xmm2c8(const xmm2c8_t x,
	                                  const double n) {
	                   
	                xmm2c8_t cp;        
	                register __m128d xmm0,xmm1;
	                register __m128d r,tht;
	                register __m128d vn,pt;
	                register __m128d ta;
	                xmm0  = _mm_mul_pd(x.re,x.re);
	                vn    = _mm_set1_pd(n);
	                xmm1  = _mm_mul_pd(x.im,x.im);
	                r     = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
	                tht   = _mm_atan_pd(_mm_div_pd(x.im,x.re));
	                pt    = _mm_pow_pd(r,vn);
	                ta    = _mm_mul_pd(vn,tht);
	                cp.re = _mm_mul_pd(pt,_mm_cos_pd(ta));
	                cp.im = _mm_mul_pd(pt,_mm_sin_pd(ta));      
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
                   ceq_xmm2c8(      const xmm2c8_t x,
                                    const xmm2c8_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_EQ_OQ);
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
                   cgt_xmm2c8(      const xmm2c8_t x,
                                    const xmm2c8_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_GT_OQ);
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
                   clt_xmm2c8(     const xmm2c8_t x,
                                    const xmm2c8_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_LT_OQ);
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
                   cneq_xmm2c8(      const xmm2c8_t x,
                                     const xmm2c8_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm_cmp_pd_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm_cmp_pd_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm2c8_t cexp_xmm2c8(const __m128d xre,
                                          const __m128d xim) {
                                     
                        xmm2c8_t cv;
                        register const __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0;
                        xmm0   = _mm_exp_pd(xre);
                        cv.re = _mm_mul_pd(xmm0,_mm_cos_pd(xre));
                        cv.im = _mm_mul_pd(xmm0,_mm_mul_pd(_mm_sin_pd(xim),I));
                        return (cv);
              }
              
              
                  
	        
	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                       
	           static inline
                   xmm2c8_t cexp_xmm2c8(const xmm2c8_t x) {
                                     
                        xmm2c8_t cv;
                        register  __m128d I = _mm_set1_pd(1.0);
                        register __m128d xmm0;
                        xmm0   = _mm_exp_pd(x.re);
                        cv.re = _mm_mul_pd(xmm0,_mm_cos_pd(x.re));
                        cv.im = _mm_mul_pd(xmm0,_mm_mul_pd(_mm_sin_pd(x.im),I));
                        return (cv);
              }


                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm2c8_t cpolar_xmm2c8(const __m128d rho,
                                            const __m128d tht) {
                                      
                        xmm2c8_t cv;
                        register __m128d xmm0,xmm1;
                        xmm0 = _mm_mul_pd(rho,_mm_cos_pd(tht));
                        cv.re  = xmm0;
                        xmm1 = _mm_mul_pd(rho,_mm_sin_pd(tht));
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
                   xmm2c8_t csqrt_xmm2c8(const __m128d xre,
                                           const __m128d xim,
                                          __m128d * __restrict wrkc) {
                                          
                       xmm2c8_t cv;
                       register __m128d xmm0,xmm1;
                       register __m128d half = _mm_set1_pd(0.5); 
                       cabs_xmm2c8(xre,xim);
                       xmm0  = _mm_mul_pd(half,_mm_add_pd(*wrkc,xre));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_pd(half,_mm_sub_pd(*wrkc,xre));
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
                   xmm2c8_t csqrt_xmm2c8(const xmm2c8_t x,
                                          __m128d * __restrict wrkc) {
                                          
                       xmm2c8_t cv;
                       register __m128d xmm0,xmm1;
                       register __m128d half = _mm_set1_pd(0.5); 
                       cabs_xmm2c8(x.re,x.im);
                       xmm0  = _mm_mul_pd(half,_mm_add_pd(*wrkc,x.re));
                       cv.re = xmm0;
                       xmm1  = _mm_mul_pd(half,_mm_sub_pd(*wrkc,x.re));
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
                   xmm2c8_t cnorm_prod_xmm2c8(  const __m128d  xre,
                                                  const __m128d  xim,
                                                  const __m128d  yre,
                                                  const __m128d  yim) {
                                             
                        xmm2c8_t cv;
                        register __m128d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        imp  = _mm_fmadd_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        cv.re = _mm_div_pd(rep,xmm2);
                        cv.im = _mm_div_pd(imp,xmm2);
                        return (cv);
             }
             
             
                 
	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   xmm2c8_t cnorm_prod_xmm2c8(  const xmm2c8_t x,
                                                  const xmm2c8_t y) {
                                             
                        xmm2c8_t cv;
                        register __m128d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm_fmsub_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        imp  = _mm_fmadd_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        xmm0 = _mm_mul_pd(rep,rep);
                        xmm1 = _mm_mul_pd(imp,imp);
                        xmm2 = _mm_sqrt_pd(_mm_add_pd(xmm0,xmm1));
                        cv.re = _mm_div_pd(rep,xmm2);
                        cv.im = _mm_div_pd(imp,xmm2);
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
                   void cmean_prod_xmm2c8(const xmm2c8_t x,
                                           const xmm2c8_t y,
                                           double * __restrict mre,
                                           double * __restrict mim) {

                        register __m128d rep,imp;
                        __m128d t0,t1;
                        constexpr double inv16 = 0.5;
                        double sre,sim;
                        sre = 0.0;
                        rep  = _mm_fmsub_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        t0   = _mm_hadd_pd(rep,rep);
                        sre  = *(double*)&t0;
                        *mre = sre*inv16;
                        sim  = 0.0;
                        imp  = _mm_fmadd_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        t1   = _mm_hadd_pd(imp,imp);
                        sim  = *(double*)&t1;
                        *mim = sim*inv16;
             }


               	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                        
	           static inline
                   void cmean_quot_xmm2c8(  const xmm2c8_t x,
                                             const xmm2c8_t y,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m128d rep,imp,den,rquot,iquot;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        sre  = 0.0;
                        rep  = _mm_fmsub_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        imp  = _mm_fmadd_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        sim  = 0.0;
                        den  = _mm_fmadd_pd(y.re,y.re,
                                               _mm_mul_pd(y.im,y.im));
                        rquot = _mm_div_pd(rep,den);
                        t0    = _mm_hadd_pd(rquot,rquot);
                        sre   = *(double*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm_div_pd(imp,den);
                        t1    = _mm_hadd_pd(iquot,iquot);
                        sim   = *(double*)&t1;
                        *mim  = sre*inv2;
              }  


               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cnorm_cprod_xmm2c8(const xmm2c8_t x,
                                            const xmm2c8_t y,
                                            __m128d * __restrict mre,
                                            __m128d * __restrict mim) {

                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmadd_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        *mre = _mm_div_pd(rep,vcmag);
                        *mim = _mm_div_pd(imp,vcmag);
             }
             
             
               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm2c8_t cnorm_cprod_xmm2c8(const __m128d xre,
                                                 const __m128d xim,
                                                 const __m128d yre,
                                                 const __m128d yim) {
                                               
                        xmm2c8_t cv;
                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmadd_pd(xre,yre,
                                               _mm_mul_pd(xim,yim));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(xim,yre,
                                               _mm_mul_pd(xre,yim));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        cv.re = _mm_div_pd(rep,vcmag);
                        cv.im = _mm_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm2c8_t cnorm_cprod_xmm2c8(const xmm2c8_t x,
                                                 const xmm2c8_t y) {
                                               
                        xmm2c8_t cv;
                        register __m128d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm_fmadd_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        magc1= _mm_mul_pd(rep,rep);
                        imp  = _mm_fmsub_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        magc2= _mm_mul_pd(imp,imp);
                        vcmag= _mm_sqrt_pd(_mm_add_pd(magc1,magc2));
                        cv.re = _mm_div_pd(rep,vcmag);
                        cv.im = _mm_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
             
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cmean_cprod_xmm2c8(const xmm2c8_t x,
                                            const xmm2c8_t y,
                                            double * __restrict mre,
                                            double * __restrict mim) {

                        register __m128d re,im;
                        __m128d t0,t1;
                        constexpr double inv2 = 0.5;
                        double sre,sim;
                        re   = _mm_fmadd_pd(x.re,y.re,
                                               _mm_mul_pd(x.im,y.im));
                        t0   = _mm_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm_fmsub_pd(x.im,y.re,
                                               _mm_mul_pd(x.re,y.im));
                        t1   = _mm_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }


                 
	          
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void arith_cmean_xmm2c8(  const __m128d xre,
                                              const __m128d xim,
                                              double * __restrict mre,
                                              double * __restrict mim) {

                        constexpr double inv2 = 0.5;
                        __m128d t0,t1;
                        double sre,sim;
                        t0   = _mm_hadd_pd(xre,xre);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_pd(xim,xim);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }
             
             
                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void arith_cmean_xmm2c8(  const xmm2c8_t x,
                                              double * __restrict mre,
                                              double * __restrict mim) {

                        constexpr double inv2 = 0.5;
                        __m128d t0,t1;
                        double sre,sim;
                        t0   = _mm_hadd_pd(x.re,x.re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm_hadd_pd(x.im,x.im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }


                       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cnormalize_xmm2c8( const xmm2c8_t x,
                                            const xmm2c8_t y,
                                            __m128d * __restrict mre,
                                            __m128d * __restrict mim ) {

                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(x.re,y.re,
                                                    _mm_mul_pd(x.im,y.im)));
                        *mre = _mm_div_pd(x.re,cvmag);
                        *mim =  _mm_div_pd(x.im,cvmag);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   xmm2c8_t cnormalize_xmm2c8( const __m128d xre,
                                                 const __m128d xim,
                                                 const __m128d yre,
                                                 const __m128d yim) {
                                            
                        xmm2c8_t cv;
                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xre,yre,
                                                    _mm_mul_pd(xim,yim)));
                        cv.re = _mm_div_pd(xre,cvmag);
                        cv.im =  _mm_div_pd(xim,cvmag);
                        return (cv);
             }
             
             
                  	        
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif	          
                  
	           static inline
                   xmm2c8_t cnormalize_xmm2c8( const xmm2c8_t x,
                                                 const xmm2c8_t y) {
                                            
                        xmm2c8_t cv;
                        register __m128d re,im,cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(x.re,y.re,
                                                    _mm_mul_pd(x.im,y.im)));
                        cv.re = _mm_div_pd(x.re,cvmag);
                        cv.im =  _mm_div_pd(x.im,cvmag);
                        return (cv);
             }


               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cmagnitude_xmm2c8(   const __m128d xre,
                                              const __m128d xim,
                                              const __m128d yre,
                                              const __m128d yim,
                                              __m128d * __restrict  mre) {

                        register __m128d cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(xre,yre,
                                                          _mm_mul_pd(xim,yim)));
                        *mre = cvmag;
             }
             
             
               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=SSE
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("sse")
#endif                  
	           static inline
                   void cmagnitude_xmm2c8(   const xmm2c8_t x,
                                              const xmm2c8_t y,
                                              __m128d * __restrict  mre) {

                        register __m128d cvmag;
                        cvmag= _mm_sqrt_pd(_mm_fmadd_pd(x.re,y.re,
                                          _mm_mul_pd(x.im,y.im)));
                        *mre = cvmag;
             }







} // math


} // gms















#endif /*__GMS_COMPLEX_XMM2R8_H__*/
