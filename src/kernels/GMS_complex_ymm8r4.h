
#ifndef __GMS_COMPLEX_YMM8R4_H__
#define __GMS_COMPLEX_YMM8R4_H__ 201020230913


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

    const unsigned int GMS_COMPLEX_YMM8R4_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_YMM8R4_MINOR = 0U;
    const unsigned int GMS_COMPLEX_YMM8R4_MICRO = 0U;
    const unsigned int GMS_COMPLEX_YMM8R4_FULLVER =
      1000U*GMS_COMPLEX_YMM8R4_MAJOR+
      100U*GMS_COMPLEX_YMM8R4_MINOR+
      10U*GMS_COMPLEX_YMM8R4_MICRO;
    const char * const GMS_COMPLEX_YMM8R4_CREATION_DATE = "20-10-2023 09:13 AM +00200 (FRI 20 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_YMM8R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_YMM8R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_YMM8R4_DESCRIPTION   = "AVX/AVX2 optimized complex number implementation.";

}

#include <cstdint>
#include <immintrin.h>
#include <utility>
#include "GMS_config.h"
#include "GMS_simd_utils.h"

namespace  gms {


       namespace math {
       
       
                   struct  alignas(32) ymm8c4_t 
                   {
                          __m256 re;
                          __m256 im;
                   };
                   


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   ymm8c4_t cadd_ymm8c4(const ymm8c4_t x,
                                        const ymm8c4_t y) {
                                     
                        ymm8c4_t cv;
                        register __m256 xmm0,xmm1;
                        xmm0   = _mm256_add_ps(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1   = _mm256_add_ps(x.im,y.im);
                        cv.im  = xmm1;  
                        return (cv);            
               }


                  
		          
          	                   	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                 
	           static inline
	           ymm8c4_t cadd_ymm8c4(const ymm8c4_t x,
                                     const __m256 s) {
                      
                      ymm8c4_t cv;
                      cv.re =  _mm256_add_ps(x.re,s);
                      cv.im =  x.im;
                      return (cv);                       
                }


                ////////////////////////////////////////////////////////////////////

               	          
	               
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   ymm8c4_t csub_ymm8c4(const ymm8c4_t x,
                                          const ymm8c4_t y) {
                                    
                        ymm8c4_t cv;
                        register __m256 xmm0,xmm1;
                        xmm0  = _mm256_sub_ps(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1  = _mm256_sub_ps(x.im,y.im);
                        cv.im  = xmm1;
                        return (cv);
                }
                
                
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   ymm8c4_t csub_ymm8c4(const ymm8c4_t x,
                                          const __m256 s) {
                                    
                        ymm8c4_t cv;
                        cv.re = _mm256_sub_ps(x.re,s);
                        cv.im = x.im;
                        return (cv);
               }
               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                   ymm8c4_t cmul_ymm8c4(const ymm8c4_t x,
                                        const ymm8c4_t y) {
                                     
                         ymm8c4_t cv;
                         register __m256 xmm0,xmm1;
                         xmm0 = _mm256_sub_ps(_mm256_mul_ps(x.re,y.re),
                                              _mm256_mul_ps(x.im,y.im));
                         cv.re  = xmm0;
                         xmm1 = _mm256_mul_ps(_mm256_mul_ps(x.im,y.re),
                                              _mm256_mul_ps(x.re,y.im));
                         cv.im  = xmm1;
                         return (cv);
                }


                 
	        
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__ 	                            
	           static inline
                   ymm8c4_t cmul_ymm8c4(const ymm8c4_t x,
                                          const __m256 s) {
                                     
                        ymm8c4_t cv;
                        cv.re = _mm256_mul_ps(x.re,s);
                        cv.im = _mm256_mul_ps(x.im,s);
                        return (cv);
               }
               


               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                   ymm8c4_t cdiv_ymm8c4(const ymm8c4_t x,
                                          const ymm8c4_t y) {
                                     
                      ymm8c4_t cv;
                      register __m256 xmm0,xmm1,xmm2;
                      xmm0 = _mm256_fmadd_ps(x.re,y.re,
                                           _mm256_mul_ps(x.im,y.im));
                      xmm1 = _mm256_fmsub_ps(x.im,y.re,
                                           _mm256_mul_ps(x.re,y.im));
                      xmm2 = _mm256_fmadd_ps(xmm0,xmm0,
                                           _mm256_mul_ps(xmm1,xmm1));
                      cv.re  = _mm256_div_ps(xmm0,xmm2);
                      cv.im  = _mm256_div_ps(xmm1,xmm2);
                      return (cv);
                }


          	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                ymm8c4_t cdiv_ymm8c4(const ymm8c4_t x,
                                     const __m256 s) {
                                     
                         ymm8c4_t cv;
                         cv.re = _mm256_div_ps(x.re,s);
                         cv.im = _mm256_div_ps(x.im,s);
                         return (cv);
               }
               
               
               
          	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                   
	           static inline
                   ymm8c4_t cdiv_ymm8c4_s(const __m256 s,
                                          const ymm8c4_t x) {
                                       
                        ymm8c4_t cv,tmp;
                        tmp.re = s;
                        tmp.im = _mm256_setzero_ps();
                        cv = cdiv_ymm8c4(tmp,x);
                        return (cv);                     
                 }
                 
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   ymm8c4_t cdiv_smith_ymm8c4(const ymm8c4_t x,
                                                const ymm8c4_t y) {
                                           
                        ymm8c4_t cv;
                        register __m256 r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_ps_mask(ymm8r4_abs(y.re),
                                                  ymm8r4_abs(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_ps(m,_mm256_div_ps(y.re,y.im),
                                                      _mm256_div_ps(y.im,y.re)); // r
                        den  = _mm256_mask_blend_ps(m,_mm256_fmadd_ps(r,y.re,y.im),
                                                      _mm256_fmadd_ps(r,y.im,y.re));
                        cv.re  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmadd_ps(x.re,r,x.im),den),
                                                _mm256_div_ps(_mm256_fmadd_ps(x.im,r,x.re),den));
                        cv.im  =  _mm256_mask_blend_ps(m,
                                                _mm256_div_ps(_mm256_fmsub_ps(x.im,r,x.re),den),
                                                _mm256_div_ps(_mm256_sub_ps(x.im,_mm256_mul_ps(x.re,r)),den));
                        return (cv);
               }          
	        
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                ymm8c4_t cdiv_smith_ymm8c4_s(const __m256 s,
                                                const ymm8c4_t x) {
                                             
                        ymm8c4_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm256_setzero_ps(); 
                        cv = cdiv_smith_ymm8c4(t0,x);
                        return (cv);                 
                 }
                 
                   
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m256 cabs_ymm8c4(const __m256 re,
                                       const __m256 im) {

                        register __m256 xmm0,xmm1,cabs;
                        xmm0 = _mm256_mul_ps(re,re);
                        xmm1 = _mm256_mul_ps(im,im);
                        cabs = _mm256_sqrt_ps(_mm256_add_ps(xmm0,xmm1));
                        return (cabs);
                 }
                 
                 
                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m256 cabs_ymm8c4(const ymm8c4_t x) {

                        register __m256 xmm0,xmm1,cabs;
                        xmm0 = _mm256_mul_ps(x.re,x.re);
                        xmm1 = _mm256_mul_ps(x.im,x.im);
                        cabs = _mm256_sqrt_ps(_mm256_add_ps(xmm0,xmm1));
                        return (cabs);
                 }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m256 carg_ymm8c4(const __m256 re,
                                       const __m256 im) {

                       register __m256 carg;
                       carg = _mm256_atan2_ps(re,im);
                       return (carg);
                }
                
                
                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   __m256 carg_ymm8c4(ymm8c4_t x) {

                       register __m256 carg;
                       carg = _mm256_atan2_ps(x.re,x.im);
                       return (carg);
                }
                
                
                 
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
	           ymm8c4_t clog_ymm8c4(const ymm8c4_t x){
	                                  
	                ymm8c4_t clog;                           
	                register __m256 t1,t2,ln;
	                t1  = cabs_ymm8c4(x.re,x.im);
	                t2  = carg_ymm8c4(x.re,x.im);
	                ln  = _mm256_log_ps(t1);
	                clog.re = ln;
	                clog.im = t2;
	                return (clog);                    
	        }

                 

                  	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t cconj_ymm8c4_v2(const __m256 xre,
                                              const __m256 xim) {                                              
                         
                        //register __m256 c;              
                        //c = negate_ymm8c4(*im);
                        //*im = c;
                        ymm8c4_t cv;
                        cv.re = xre; 
                        cv.im = negate_ymm8r4(xim);
                        return (cv);
                   } 
                   
                   
                
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t cconj_ymm8c4_v2(const ymm8c4_t x) {                                              
                         
                        //register __m256 c;              
                        //c = negate_ymm8c4(*im);
                        //*im = c;
                        ymm8c4_t cv;
                        cv.re = x.re; 
                        cv.im = negate_ymm8r4(x.im);
                        return (cv);
                   } 
                   


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t ccos_ymm8c4(const __m256 xre,
                                          const __m256 xim) {
                                    
                      ymm8c4_t cv;
                      register __m256 xmm0,xmm1;
                      xmm0  = _mm256_mul_ps(_mm256_cos_ps(xre),_mm256_cosh_ps(xim));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_ps(_mm256_sin_ps(xre),_mm256_sinh_ps(xim));
                      cv.im = xmm1;
                      return (cv); 
               }
               
               
                 	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t ccos_ymm8c4(const ymm8c4_t x) {
                                    
                      ymm8c4_t cv;
                      register __m256 xmm0,xmm1;
                      xmm0  = _mm256_mul_ps(_mm256_cos_ps(x.re),_mm256_cosh_ps(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_ps(_mm256_sin_ps(x.re),_mm256_sinh_ps(x.im));
                      cv.im = xmm1;
                      return (cv); 
               }


                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t ccosh_ymm8c4(const __m256 xre,
                                           const __m256 xim) {
                                          
                      ymm8c4_t cv;
                      register __m256 xmm0,xmm1;
                      xmm0  = _mm256_mul_ps(_mm256_cosh_ps(xre),_mm256_cos_ps(xim));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_ps(_mm256_sinh_ps(xre),_mm256_sin_ps(xim));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   
	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t ccosh_ymm8c4(const ymm8c4_t x) {
                                          
                      ymm8c4_t cv;
                      register __m256 xmm0,xmm1;
                      xmm0  = _mm256_mul_ps(_mm256_cosh_ps(x.re),_mm256_cos_ps(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_ps(_mm256_sinh_ps(x.re),_mm256_sin_ps(x.im));
                      cv.im = xmm1; 
                      return (cv);
               }
               
               
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
	           ymm8c4_t cpow_ymm8c4(const ymm8c4_t x,
	                                  const float n) {
	                   
	                ymm8c4_t cp;        
	                register __m256 xmm0,xmm1;
	                register __m256 r,tht;
	                register __m256 vn,pt;
	                register __m256 ta;
	                xmm0  = _mm256_mul_ps(x.re,x.re);
	                vn    = _mm256_set1_ps(n);
	                xmm1  = _mm256_mul_ps(x.im,x.im);
	                r     = _mm256_sqrt_ps(_mm256_add_ps(xmm0,xmm1));
	                tht   = _mm256_atan_ps(_mm256_div_ps(x.im,x.re));
	                pt    = _mm256_pow_ps(r,vn);
	                ta    = _mm256_mul_ps(vn,tht);
	                cp.re = _mm256_mul_ps(pt,_mm256_cos_ps(ta));
	                cp.im = _mm256_mul_ps(pt,_mm256_sin_ps(ta));      
	                return (cp);              
	       }
	       


                
              
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   ceq_ymm8c4(      const ymm8c4_t x,
                                    const ymm8c4_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_EQ_OQ);
                         return std::make_pair(eqr,eqi);
              }
              


                 	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cgt_ymm8c4(      const ymm8c4_t x,
                                    const ymm8c4_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_GT_OQ);
                         return std::make_pair(eqr,eqi);
              }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   clt_ymm8c4(     const ymm8c4_t x,
                                    const ymm8c4_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_LT_OQ);
                         return std::make_pair(eqr,eqi);
              }



                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                       
	           static inline
	           std::pair<__mmask8,__mmask8> 
                   cneq_ymm8c4(      const ymm8c4_t x,
                                     const ymm8c4_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm256_cmp_ps_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm256_cmp_ps_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm8c4_t cexp_ymm8c4(const __m256 xre,
                                          const __m256 xim) {
                                     
                        ymm8c4_t cv;
                        register const __m256 I = _mm256_set1_ps(1.0);
                        register __m256 xmm0;
                        xmm0   = _mm256_exp_ps(xre);
                        cv.re = _mm256_mul_ps(xmm0,_mm256_cos_ps(xre));
                        cv.im = _mm256_mul_ps(xmm0,_mm256_mul_ps(_mm256_sin_ps(xim),I));
                        return (cv);
              }
              
              
                  
	        
	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                       
	           static inline
                   ymm8c4_t cexp_ymm8c4(const ymm8c4_t x) {
                                     
                        ymm8c4_t cv;
                        register  __m256 I = _mm256_set1_ps(1.0);
                        register __m256 xmm0;
                        xmm0   = _mm256_exp_ps(x.re);
                        cv.re = _mm256_mul_ps(xmm0,_mm256_cos_ps(x.re));
                        cv.im = _mm256_mul_ps(xmm0,_mm256_mul_ps(_mm256_sin_ps(x.im),I));
                        return (cv);
              }


                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm8c4_t cpolar_ymm8c4(const __m256 rho,
                                            const __m256 tht) {
                                      
                        ymm8c4_t cv;
                        register __m256 xmm0,xmm1;
                        xmm0 = _mm256_mul_ps(rho,_mm256_cos_ps(tht));
                        cv.re  = xmm0;
                        xmm1 = _mm256_mul_ps(rho,_mm256_sin_ps(tht));
                        cv.im  = xmm1;
                        return (cv);
              }
              


     	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                       
	           static inline
                   ymm8c4_t csqrt_ymm8c4(const __m256 xre,
                                           const __m256 xim,
                                          __m256 * __restrict wrkc) {
                                          
                       ymm8c4_t cv;
                       register __m256 xmm0,xmm1;
                       register __m256 half = _mm256_set1_ps(0.5); 
                       cabs_ymm8c4(xre,xim);
                       xmm0  = _mm256_mul_ps(half,_mm256_add_ps(*wrkc,xre));
                       cv.re = xmm0;
                       xmm1  = _mm256_mul_ps(half,_mm256_sub_ps(*wrkc,xre));
                       cv.im = xmm1; 
                       return (cv);
              }
              
              
                   
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm8c4_t csqrt_ymm8c4(const ymm8c4_t x,
                                          __m256 * __restrict wrkc) {
                                          
                       ymm8c4_t cv;
                       register __m256 xmm0,xmm1;
                       register __m256 half = _mm256_set1_ps(0.5); 
                       cabs_ymm8c4(x.re,x.im);
                       xmm0  = _mm256_mul_ps(half,_mm256_add_ps(*wrkc,x.re));
                       cv.re = xmm0;
                       xmm1  = _mm256_mul_ps(half,_mm256_sub_ps(*wrkc,x.re));
                       cv.im = xmm1; 
                       return (cv);
              }


                  
	        
	     	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm8c4_t cnorm_prod_ymm8c4(  const __m256  xre,
                                                  const __m256  xim,
                                                  const __m256  yre,
                                                  const __m256  yim) {
                                             
                        ymm8c4_t cv;
                        register __m256 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm256_fmsub_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        imp  = _mm256_fmadd_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        xmm0 = _mm256_mul_ps(rep,rep);
                        xmm1 = _mm256_mul_ps(imp,imp);
                        xmm2 = _mm256_sqrt_ps(_mm256_add_ps(xmm0,xmm1));
                        cv.re = _mm256_div_ps(rep,xmm2);
                        cv.im = _mm256_div_ps(imp,xmm2);
                        return (cv);
             }
             
             
                 
	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm8c4_t cnorm_prod_ymm8c4(  const ymm8c4_t x,
                                                  const ymm8c4_t y) {
                                             
                        ymm8c4_t cv;
                        register __m256 rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm256_fmsub_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        imp  = _mm256_fmadd_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        xmm0 = _mm256_mul_ps(rep,rep);
                        xmm1 = _mm256_mul_ps(imp,imp);
                        xmm2 = _mm256_sqrt_ps(_mm256_add_ps(xmm0,xmm1));
                        cv.re = _mm256_div_ps(rep,xmm2);
                        cv.im = _mm256_div_ps(imp,xmm2);
                        return (cv);
             }
             

                  	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   void cmean_prod_ymm8c4(const ymm8c4_t x,
                                           const ymm8c4_t y,
                                           float * __restrict mre,
                                           float * __restrict mim) {

                        register __m256 rep,imp;
                        __m256 t0,t1;
                        constexpr float inv16 = 0.125;
                        float sre,sim;
                        sre = 0.0;
                        rep  = _mm256_fmsub_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        t0   = _mm256_hadd_ps(rep,rep);
                        sre  = *(float*)&t0;
                        *mre = sre*inv16;
                        sim  = 0.0;
                        imp  = _mm256_fmadd_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        t1   = _mm256_hadd_ps(imp,imp);
                        sim  = *(float*)&t1;
                        *mim = sim*inv16;
             }


               	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   void cmean_quot_ymm8c4(  const ymm8c4_t x,
                                             const ymm8c4_t y,
                                             float * __restrict mre,
                                             float * __restrict mim) {

                        register __m256 rep,imp,den,rquot,iquot;
                        __m256 t0,t1;
                        constexpr float inv2 = 0.125;
                        float sre,sim;
                        sre  = 0.0;
                        rep  = _mm256_fmsub_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        imp  = _mm256_fmadd_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        sim  = 0.0;
                        den  = _mm256_fmadd_ps(y.re,y.re,
                                               _mm256_mul_ps(y.im,y.im));
                        rquot = _mm256_div_ps(rep,den);
                        t0    = _mm256_hadd_ps(rquot,rquot);
                        sre   = *(float*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm256_div_ps(imp,den);
                        t1    = _mm256_hadd_ps(iquot,iquot);
                        sim   = *(float*)&t1;
                        *mim  = sre*inv2;
              }  


               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cnorm_cprod_ymm8c4(const ymm8c4_t x,
                                            const ymm8c4_t y,
                                            __m256 * __restrict mre,
                                            __m256 * __restrict mim) {

                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmadd_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        *mre = _mm256_div_ps(rep,vcmag);
                        *mim = _mm256_div_ps(imp,vcmag);
             }
             
             
               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t cnorm_cprod_ymm8c4(const __m256 xre,
                                                 const __m256 xim,
                                                 const __m256 yre,
                                                 const __m256 yim) {
                                               
                        ymm8c4_t cv;
                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmadd_ps(xre,yre,
                                               _mm256_mul_ps(xim,yim));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(xim,yre,
                                               _mm256_mul_ps(xre,yim));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        cv.re = _mm256_div_ps(rep,vcmag);
                        cv.im = _mm256_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t cnorm_cprod_ymm8c4(const ymm8c4_t x,
                                                 const ymm8c4_t y) {
                                               
                        ymm8c4_t cv;
                        register __m256 rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmadd_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        magc1= _mm256_mul_ps(rep,rep);
                        imp  = _mm256_fmsub_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        magc2= _mm256_mul_ps(imp,imp);
                        vcmag= _mm256_sqrt_ps(_mm256_add_ps(magc1,magc2));
                        cv.re = _mm256_div_ps(rep,vcmag);
                        cv.im = _mm256_div_ps(imp,vcmag);
                        return (cv);
             }
             
             
             
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cmean_cprod_ymm8c4(const ymm8c4_t x,
                                            const ymm8c4_t y,
                                            float * __restrict mre,
                                            float * __restrict mim) {

                        register __m256 re,im;
                        __m256 t0,t1;
                        constexpr float inv2 = 0.125;
                        float sre,sim;
                        re   = _mm256_fmadd_ps(x.re,y.re,
                                               _mm256_mul_ps(x.im,y.im));
                        t0   = _mm256_hadd_ps(re,re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        im   = _mm256_fmsub_ps(x.im,y.re,
                                               _mm256_mul_ps(x.re,y.im));
                        t1   = _mm256_hadd_ps(im,im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2;
             }


                 
	          
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void arith_cmean_ymm8c4(  const __m256 xre,
                                              const __m256 xim,
                                              float * __restrict mre,
                                              float * __restrict mim) {

                        constexpr float inv2 = 0.125;
                        __m256 t0,t1;
                        float sre,sim;
                        t0   = _mm256_hadd_ps(xre,xre);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm256_hadd_ps(xim,xim);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }
             
             
                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void arith_cmean_ymm8c4(  const ymm8c4_t x,
                                              float * __restrict mre,
                                              float * __restrict mim) {

                        constexpr float inv2 = 0.125;
                        __m256 t0,t1;
                        float sre,sim;
                        t0   = _mm256_hadd_ps(x.re,x.re);
                        sre  = *(float*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm256_hadd_ps(x.im,x.im);
                        sim  = *(float*)&t1;
                        *mim = sim*inv2; 
             }


                       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cnormalize_ymm8c4( const ymm8c4_t x,
                                            const ymm8c4_t y,
                                            __m256 * __restrict mre,
                                            __m256 * __restrict mim ) {

                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(x.re,y.re,
                                                    _mm256_mul_ps(x.im,y.im)));
                        *mre = _mm256_div_ps(x.re,cvmag);
                        *mim =  _mm256_div_ps(x.im,cvmag);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm8c4_t cnormalize_ymm8c4( const __m256 xre,
                                                 const __m256 xim,
                                                 const __m256 yre,
                                                 const __m256 yim) {
                                            
                        ymm8c4_t cv;
                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(xre,yre,
                                                    _mm256_mul_ps(xim,yim)));
                        cv.re = _mm256_div_ps(xre,cvmag);
                        cv.im =  _mm256_div_ps(xim,cvmag);
                        return (cv);
             }
             
             
                  	        
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif	          
                  
	           static inline
                   ymm8c4_t cnormalize_ymm8c4( const ymm8c4_t x,
                                                 const ymm8c4_t y) {
                                            
                        ymm8c4_t cv;
                        register __m256 re,im,cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(x.re,y.re,
                                                    _mm256_mul_ps(x.im,y.im)));
                        cv.re = _mm256_div_ps(x.re,cvmag);
                        cv.im =  _mm256_div_ps(x.im,cvmag);
                        return (cv);
             }


               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cmagnitude_ymm8c4(   const __m256 xre,
                                              const __m256 xim,
                                              const __m256 yre,
                                              const __m256 yim,
                                              __m256 * __restrict  mre) {

                        register __m256 cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(xre,yre,
                                                          _mm256_mul_ps(xim,yim)));
                        *mre = cvmag;
             }
             
             
               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cmagnitude_ymm8c4(   const ymm8c4_t x,
                                              const ymm8c4_t y,
                                              __m256 * __restrict  mre) {

                        register __m256 cvmag;
                        cvmag= _mm256_sqrt_ps(_mm256_fmadd_ps(x.re,y.re,
                                          _mm256_mul_ps(x.im,y.im)));
                        *mre = cvmag;
             }







} // math


} // gms

#endif /*__GMS_COMPLEX_YMM4R8_H__*/