
#ifndef __GMS_COMPLEX_YMM4R8_H__
#define __GMS_COMPLEX_YMM4R8_H__ 181020231557


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

    const unsigned int GMS_COMPLEX_YMM4R8_MAJOR = 1U;
    const unsigned int GMS_COMPLEX_YMM4R8_MINOR = 0U;
    const unsigned int GMS_COMPLEX_YMM4R8_MICRO = 0U;
    const unsigned int GMS_COMPLEX_YMM4R8_FULLVER =
      1000U*GMS_COMPLEX_YMM4R8_MAJOR+
      100U*GMS_COMPLEX_YMM4R8_MINOR+
      10U*GMS_COMPLEX_YMM4R8_MICRO;
    const char * const GMS_COMPLEX_YMM4R8_CREATION_DATE = "20-10-2023 09:13 AM +00200 (FRI 20 OCT 2022 GMT+2)";
    const char * const GMS_COMPLEX_YMM4R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_COMPLEX_YMM4R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_COMPLEX_YMM4R8_DESCRIPTION   = "AVX/AVX2 optimized complex number implementation.";

}

#include <cstdint>
#include <immintrin.h>
#include <utility>
#include "GMS_config.h"
#include "GMS_simd_utils.h"

namespace  gms {


       namespace math {
       
       
                   struct  alignas(32) ymm4c8_t 
                   {
                          __m256d re;
                          __m256d im;
                   };
                   


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif
                __ATTR_ALWAYS_INLINE__                  
	           static inline
                   ymm4c8_t cadd_ymm4c8(const ymm4c8_t x,
                                        const ymm4c8_t y) {
                                     
                        ymm4c8_t cv;
                        register __m256d xmm0,xmm1;
                        xmm0   = _mm256_add_pd(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1   = _mm256_add_pd(x.im,y.im);
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
	           ymm4c8_t cadd_ymm4c8(const ymm4c8_t x,
                                     const __m256d s) {
                      
                      ymm4c8_t cv;
                      cv.re =  _mm256_add_pd(x.re,s);
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
                   ymm4c8_t csub_ymm4c8(const ymm4c8_t x,
                                          const ymm4c8_t y) {
                                    
                        ymm4c8_t cv;
                        register __m256d xmm0,xmm1;
                        xmm0  = _mm256_sub_pd(x.re,y.re);
                        cv.re  = xmm0;
                        xmm1  = _mm256_sub_pd(x.im,y.im);
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
                   ymm4c8_t csub_ymm4c8(const ymm4c8_t x,
                                          const __m256d s) {
                                    
                        ymm4c8_t cv;
                        cv.re = _mm256_sub_pd(x.re,s);
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
                   ymm4c8_t cmul_ymm4c8(const ymm4c8_t x,
                                        const ymm4c8_t y) {
                                     
                         ymm4c8_t cv;
                         register __m256d xmm0,xmm1;
                         xmm0 = _mm256_sub_pd(_mm256_mul_pd(x.re,y.re),
                                              _mm256_mul_pd(x.im,y.im));
                         cv.re  = xmm0;
                         xmm1 = _mm256_mul_pd(_mm256_mul_pd(x.im,y.re),
                                              _mm256_mul_pd(x.re,y.im));
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
                   ymm4c8_t cmul_ymm4c8(const ymm4c8_t x,
                                          const __m256d s) {
                                     
                        ymm4c8_t cv;
                        cv.re = _mm256_mul_pd(x.re,s);
                        cv.im = _mm256_mul_pd(x.im,s);
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
                   ymm4c8_t cdiv_ymm4c8(const ymm4c8_t x,
                                          const ymm4c8_t y) {
                                     
                      ymm4c8_t cv;
                      register __m256d xmm0,xmm1,xmm2;
                      xmm0 = _mm256_fmadd_pd(x.re,y.re,
                                           _mm256_mul_pd(x.im,y.im));
                      xmm1 = _mm256_fmsub_pd(x.im,y.re,
                                           _mm256_mul_pd(x.re,y.im));
                      xmm2 = _mm256_fmadd_pd(xmm0,xmm0,
                                           _mm256_mul_pd(xmm1,xmm1));
                      cv.re  = _mm256_div_pd(xmm0,xmm2);
                      cv.im  = _mm256_div_pd(xmm1,xmm2);
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
                ymm4c8_t cdiv_ymm4c8(const ymm4c8_t x,
                                     const __m256d s) {
                                     
                         ymm4c8_t cv;
                         cv.re = _mm256_div_pd(x.re,s);
                         cv.im = _mm256_div_pd(x.im,s);
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
                   ymm4c8_t cdiv_ymm4c8_s(const __m256d s,
                                          const ymm4c8_t x) {
                                       
                        ymm4c8_t cv,tmp;
                        tmp.re = s;
                        tmp.im = _mm256_setzero_pd();
                        cv = cdiv_ymm4c8(tmp,x);
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
                   ymm4c8_t cdiv_smith_ymm4c8(const ymm4c8_t x,
                                                const ymm4c8_t y) {
                                           
                        ymm4c8_t cv;
                        register __m256d r,den;
                        __mmask8 m = 0x0;
                        m    = _mm256_cmp_pd_mask(ymm4r8_abs(y.re),
                                                  ymm4r8_abs(y.im),
                                                  _CMP_GE_OQ);
                        r    = _mm256_mask_blend_pd(m,_mm256_div_pd(y.re,y.im),
                                                      _mm256_div_pd(y.im,y.re)); // r
                        den  = _mm256_mask_blend_pd(m,_mm256_fmadd_pd(r,y.re,y.im),
                                                      _mm256_fmadd_pd(r,y.im,y.re));
                        cv.re  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmadd_pd(x.re,r,x.im),den),
                                                _mm256_div_pd(_mm256_fmadd_pd(x.im,r,x.re),den));
                        cv.im  =  _mm256_mask_blend_pd(m,
                                                _mm256_div_pd(_mm256_fmsub_pd(x.im,r,x.re),den),
                                                _mm256_div_pd(_mm256_sub_pd(x.im,_mm256_mul_pd(x.re,r)),den));
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
                ymm4c8_t cdiv_smith_ymm4c8_s(const __m256d s,
                                                const ymm4c8_t x) {
                                             
                        ymm4c8_t cv,t0;                    
                        t0.re = s;
                        t0.im = _mm256_setzero_pd(); 
                        cv = cdiv_smith_ymm4c8(t0,x);
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
                   __m256d cabs_ymm4c8(const __m256d re,
                                       const __m256d im) {

                        register __m256d xmm0,xmm1,cabs;
                        xmm0 = _mm256_mul_pd(re,re);
                        xmm1 = _mm256_mul_pd(im,im);
                        cabs = _mm256_sqrt_pd(_mm256_add_pd(xmm0,xmm1));
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
                   __m256d cabs_ymm4c8(const ymm4c8_t x) {

                        register __m256d xmm0,xmm1,cabs;
                        xmm0 = _mm256_mul_pd(x.re,x.re);
                        xmm1 = _mm256_mul_pd(x.im,x.im);
                        cabs = _mm256_sqrt_pd(_mm256_add_pd(xmm0,xmm1));
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
                   __m256d carg_ymm4c8(const __m256d re,
                                       const __m256d im) {

                       register __m256d carg;
                       carg = _mm256_atan2_pd(re,im);
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
                   __m256d carg_ymm4c8(ymm4c8_t x) {

                       register __m256d carg;
                       carg = _mm256_atan2_pd(x.re,x.im);
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
	           ymm4c8_t clog_ymm4c8(const ymm4c8_t x){
	                                  
	                ymm4c8_t clog;                           
	                register __m256d t1,t2,ln;
	                t1  = cabs_ymm4c8(x.re,x.im);
	                t2  = carg_ymm4c8(x.re,x.im);
	                ln  = _mm256_log_pd(t1);
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
                   ymm4c8_t cconj_ymm4c8_v2(const __m256d xre,
                                              const __m256d xim) {                                              
                         
                        //register __m256d c;              
                        //c = negate_ymm4c8(*im);
                        //*im = c;
                        ymm4c8_t cv;
                        cv.re = xre; 
                        cv.im = negate_ymm4r8(xim);
                        return (cv);
                   } 
                   
                   
                
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm4c8_t cconj_ymm4c8_v2(const ymm4c8_t x) {                                              
                         
                        //register __m256d c;              
                        //c = negate_ymm4c8(*im);
                        //*im = c;
                        ymm4c8_t cv;
                        cv.re = x.re; 
                        cv.im = negate_ymm4r8(x.im);
                        return (cv);
                   } 
                   


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm4c8_t ccos_ymm4c8(const __m256d xre,
                                          const __m256d xim) {
                                    
                      ymm4c8_t cv;
                      register __m256d xmm0,xmm1;
                      xmm0  = _mm256_mul_pd(_mm256_cos_pd(xre),_mm256_cosh_pd(xim));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_pd(_mm256_sin_pd(xre),_mm256_sinh_pd(xim));
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
                   ymm4c8_t ccos_ymm4c8(const ymm4c8_t x) {
                                    
                      ymm4c8_t cv;
                      register __m256d xmm0,xmm1;
                      xmm0  = _mm256_mul_pd(_mm256_cos_pd(x.re),_mm256_cosh_pd(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_pd(_mm256_sin_pd(x.re),_mm256_sinh_pd(x.im));
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
                   ymm4c8_t ccosh_ymm4c8(const __m256d xre,
                                           const __m256d xim) {
                                          
                      ymm4c8_t cv;
                      register __m256d xmm0,xmm1;
                      xmm0  = _mm256_mul_pd(_mm256_cosh_pd(xre),_mm256_cos_pd(xim));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_pd(_mm256_sinh_pd(xre),_mm256_sin_pd(xim));
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
                   ymm4c8_t ccosh_ymm4c8(const ymm4c8_t x) {
                                          
                      ymm4c8_t cv;
                      register __m256d xmm0,xmm1;
                      xmm0  = _mm256_mul_pd(_mm256_cosh_pd(x.re),_mm256_cos_pd(x.im));
                      cv.re = xmm0;
                      xmm1  = _mm256_mul_pd(_mm256_sinh_pd(x.re),_mm256_sin_pd(x.im));
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
	           ymm4c8_t cpow_ymm4c8(const ymm4c8_t x,
	                                  const double n) {
	                   
	                ymm4c8_t cp;        
	                register __m256d xmm0,xmm1;
	                register __m256d r,tht;
	                register __m256d vn,pt;
	                register __m256d ta;
	                xmm0  = _mm256_mul_pd(x.re,x.re);
	                vn    = _mm256_set1_pd(n);
	                xmm1  = _mm256_mul_pd(x.im,x.im);
	                r     = _mm256_sqrt_pd(_mm256_add_pd(xmm0,xmm1));
	                tht   = _mm256_atan_pd(_mm256_div_pd(x.im,x.re));
	                pt    = _mm256_pow_pd(r,vn);
	                ta    = _mm256_mul_pd(vn,tht);
	                cp.re = _mm256_mul_pd(pt,_mm256_cos_pd(ta));
	                cp.im = _mm256_mul_pd(pt,_mm256_sin_pd(ta));      
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
                   ceq_ymm4c8(      const ymm4c8_t x,
                                    const ymm4c8_t y) {
                                   
                          __mmask8 eqr;
                          __mmask8 eqi;
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_EQ_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_EQ_OQ);
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
                   cgt_ymm4c8(      const ymm4c8_t x,
                                    const ymm4c8_t y) {
                           
                         __mmask8 eqr;
                         __mmask8 eqi;         
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_GT_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_GT_OQ);
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
                   clt_ymm4c8(     const ymm4c8_t x,
                                    const ymm4c8_t y) {
                                   
                         __mmask8 eqr;
                         __mmask8 eqi; 
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_LT_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_LT_OQ);
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
                   cneq_ymm4c8(      const ymm4c8_t x,
                                     const ymm4c8_t y) {
                                    
                         __mmask8 eqr;
                         __mmask8 eqi;
                         eqr = _mm256_cmp_pd_mask(x.re,y.re,_CMP_NEQ_OQ);
                         eqi = _mm256_cmp_pd_mask(x.im,y.im,_CMP_NEQ_OQ);
                         return std::make_pair(eqr,eqi);
              }


                  
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm4c8_t cexp_ymm4c8(const __m256d xre,
                                          const __m256d xim) {
                                     
                        ymm4c8_t cv;
                        register const __m256d I = _mm256_set1_pd(1.0);
                        register __m256d xmm0;
                        xmm0   = _mm256_exp_pd(xre);
                        cv.re = _mm256_mul_pd(xmm0,_mm256_cos_pd(xre));
                        cv.im = _mm256_mul_pd(xmm0,_mm256_mul_pd(_mm256_sin_pd(xim),I));
                        return (cv);
              }
              
              
                  
	        
	          
 #if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                       
	           static inline
                   ymm4c8_t cexp_ymm4c8(const ymm4c8_t x) {
                                     
                        ymm4c8_t cv;
                        register  __m256d I = _mm256_set1_pd(1.0);
                        register __m256d xmm0;
                        xmm0   = _mm256_exp_pd(x.re);
                        cv.re = _mm256_mul_pd(xmm0,_mm256_cos_pd(x.re));
                        cv.im = _mm256_mul_pd(xmm0,_mm256_mul_pd(_mm256_sin_pd(x.im),I));
                        return (cv);
              }


                 
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm4c8_t cpolar_ymm4c8(const __m256d rho,
                                            const __m256d tht) {
                                      
                        ymm4c8_t cv;
                        register __m256d xmm0,xmm1;
                        xmm0 = _mm256_mul_pd(rho,_mm256_cos_pd(tht));
                        cv.re  = xmm0;
                        xmm1 = _mm256_mul_pd(rho,_mm256_sin_pd(tht));
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
                   ymm4c8_t csqrt_ymm4c8(const __m256d xre,
                                           const __m256d xim,
                                          __m256d * __restrict wrkc) {
                                          
                       ymm4c8_t cv;
                       register __m256d xmm0,xmm1;
                       register __m256d half = _mm256_set1_pd(0.5); 
                       cabs_ymm4c8(xre,xim);
                       xmm0  = _mm256_mul_pd(half,_mm256_add_pd(*wrkc,xre));
                       cv.re = xmm0;
                       xmm1  = _mm256_mul_pd(half,_mm256_sub_pd(*wrkc,xre));
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
                   ymm4c8_t csqrt_ymm4c8(const ymm4c8_t x,
                                          __m256d * __restrict wrkc) {
                                          
                       ymm4c8_t cv;
                       register __m256d xmm0,xmm1;
                       register __m256d half = _mm256_set1_pd(0.5); 
                       cabs_ymm4c8(x.re,x.im);
                       xmm0  = _mm256_mul_pd(half,_mm256_add_pd(*wrkc,x.re));
                       cv.re = xmm0;
                       xmm1  = _mm256_mul_pd(half,_mm256_sub_pd(*wrkc,x.re));
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
                   ymm4c8_t cnorm_prod_ymm4c8(  const __m256d  xre,
                                                  const __m256d  xim,
                                                  const __m256d  yre,
                                                  const __m256d  yim) {
                                             
                        ymm4c8_t cv;
                        register __m256d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm256_fmsub_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        imp  = _mm256_fmadd_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        xmm0 = _mm256_mul_pd(rep,rep);
                        xmm1 = _mm256_mul_pd(imp,imp);
                        xmm2 = _mm256_sqrt_pd(_mm256_add_pd(xmm0,xmm1));
                        cv.re = _mm256_div_pd(rep,xmm2);
                        cv.im = _mm256_div_pd(imp,xmm2);
                        return (cv);
             }
             
             
                 
	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   ymm4c8_t cnorm_prod_ymm4c8(  const ymm4c8_t x,
                                                  const ymm4c8_t y) {
                                             
                        ymm4c8_t cv;
                        register __m256d rep,imp,xmm0,xmm1,xmm2;
                        rep  = _mm256_fmsub_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        imp  = _mm256_fmadd_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        xmm0 = _mm256_mul_pd(rep,rep);
                        xmm1 = _mm256_mul_pd(imp,imp);
                        xmm2 = _mm256_sqrt_pd(_mm256_add_pd(xmm0,xmm1));
                        cv.re = _mm256_div_pd(rep,xmm2);
                        cv.im = _mm256_div_pd(imp,xmm2);
                        return (cv);
             }
             

                  	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   void cmean_prod_ymm4c8(const ymm4c8_t x,
                                           const ymm4c8_t y,
                                           double * __restrict mre,
                                           double * __restrict mim) {

                        register __m256d rep,imp;
                        __m256d t0,t1;
                        constexpr double inv16 = 0.25;
                        double sre,sim;
                        sre = 0.0;
                        rep  = _mm256_fmsub_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        t0   = _mm256_hadd_pd(rep,rep);
                        sre  = *(double*)&t0;
                        *mre = sre*inv16;
                        sim  = 0.0;
                        imp  = _mm256_fmadd_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        t1   = _mm256_hadd_pd(imp,imp);
                        sim  = *(double*)&t1;
                        *mim = sim*inv16;
             }


               	        
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                        
	           static inline
                   void cmean_quot_ymm4c8(  const ymm4c8_t x,
                                             const ymm4c8_t y,
                                             double * __restrict mre,
                                             double * __restrict mim) {

                        register __m256d rep,imp,den,rquot,iquot;
                        __m256d t0,t1;
                        constexpr double inv2 = 0.25;
                        double sre,sim;
                        sre  = 0.0;
                        rep  = _mm256_fmsub_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        imp  = _mm256_fmadd_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        sim  = 0.0;
                        den  = _mm256_fmadd_pd(y.re,y.re,
                                               _mm256_mul_pd(y.im,y.im));
                        rquot = _mm256_div_pd(rep,den);
                        t0    = _mm256_hadd_pd(rquot,rquot);
                        sre   = *(double*)&t0;
                        *mre  = sre*inv2;
                        iquot = _mm256_div_pd(imp,den);
                        t1    = _mm256_hadd_pd(iquot,iquot);
                        sim   = *(double*)&t1;
                        *mim  = sre*inv2;
              }  


               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cnorm_cprod_ymm4c8(const ymm4c8_t x,
                                            const ymm4c8_t y,
                                            __m256d * __restrict mre,
                                            __m256d * __restrict mim) {

                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmadd_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        *mre = _mm256_div_pd(rep,vcmag);
                        *mim = _mm256_div_pd(imp,vcmag);
             }
             
             
               	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm4c8_t cnorm_cprod_ymm4c8(const __m256d xre,
                                                 const __m256d xim,
                                                 const __m256d yre,
                                                 const __m256d yim) {
                                               
                        ymm4c8_t cv;
                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmadd_pd(xre,yre,
                                               _mm256_mul_pd(xim,yim));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(xim,yre,
                                               _mm256_mul_pd(xre,yim));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        cv.re = _mm256_div_pd(rep,vcmag);
                        cv.im = _mm256_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm4c8_t cnorm_cprod_ymm4c8(const ymm4c8_t x,
                                                 const ymm4c8_t y) {
                                               
                        ymm4c8_t cv;
                        register __m256d rep,imp,magc1,magc2,vcmag;
                        rep  = _mm256_fmadd_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        magc1= _mm256_mul_pd(rep,rep);
                        imp  = _mm256_fmsub_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        magc2= _mm256_mul_pd(imp,imp);
                        vcmag= _mm256_sqrt_pd(_mm256_add_pd(magc1,magc2));
                        cv.re = _mm256_div_pd(rep,vcmag);
                        cv.im = _mm256_div_pd(imp,vcmag);
                        return (cv);
             }
             
             
             
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cmean_cprod_ymm4c8(const ymm4c8_t x,
                                            const ymm4c8_t y,
                                            double * __restrict mre,
                                            double * __restrict mim) {

                        register __m256d re,im;
                        __m256d t0,t1;
                        constexpr double inv2 = 0.25;
                        double sre,sim;
                        re   = _mm256_fmadd_pd(x.re,y.re,
                                               _mm256_mul_pd(x.im,y.im));
                        t0   = _mm256_hadd_pd(re,re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        im   = _mm256_fmsub_pd(x.im,y.re,
                                               _mm256_mul_pd(x.re,y.im));
                        t1   = _mm256_hadd_pd(im,im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2;
             }


                 
	          
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void arith_cmean_ymm4c8(  const __m256d xre,
                                              const __m256d xim,
                                              double * __restrict mre,
                                              double * __restrict mim) {

                        constexpr double inv2 = 0.25;
                        __m256d t0,t1;
                        double sre,sim;
                        t0   = _mm256_hadd_pd(xre,xre);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm256_hadd_pd(xim,xim);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }
             
             
                           
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void arith_cmean_ymm4c8(  const ymm4c8_t x,
                                              double * __restrict mre,
                                              double * __restrict mim) {

                        constexpr double inv2 = 0.25;
                        __m256d t0,t1;
                        double sre,sim;
                        t0   = _mm256_hadd_pd(x.re,x.re);
                        sre  = *(double*)&t0;
                        *mre = sre*inv2;
                        t1   = _mm256_hadd_pd(x.im,x.im);
                        sim  = *(double*)&t1;
                        *mim = sim*inv2; 
             }


                       
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cnormalize_ymm4c8( const ymm4c8_t x,
                                            const ymm4c8_t y,
                                            __m256d * __restrict mre,
                                            __m256d * __restrict mim ) {

                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(x.re,y.re,
                                                    _mm256_mul_pd(x.im,y.im)));
                        *mre = _mm256_div_pd(x.re,cvmag);
                        *mim =  _mm256_div_pd(x.im,cvmag);
             }
             
             
                   	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   ymm4c8_t cnormalize_ymm4c8( const __m256d xre,
                                                 const __m256d xim,
                                                 const __m256d yre,
                                                 const __m256d yim) {
                                            
                        ymm4c8_t cv;
                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(xre,yre,
                                                    _mm256_mul_pd(xim,yim)));
                        cv.re = _mm256_div_pd(xre,cvmag);
                        cv.im =  _mm256_div_pd(xim,cvmag);
                        return (cv);
             }
             
             
                  	        
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif	          
                  
	           static inline
                   ymm4c8_t cnormalize_ymm4c8( const ymm4c8_t x,
                                                 const ymm4c8_t y) {
                                            
                        ymm4c8_t cv;
                        register __m256d re,im,cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(x.re,y.re,
                                                    _mm256_mul_pd(x.im,y.im)));
                        cv.re = _mm256_div_pd(x.re,cvmag);
                        cv.im =  _mm256_div_pd(x.im,cvmag);
                        return (cv);
             }


               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cmagnitude_ymm4c8(   const __m256d xre,
                                              const __m256d xim,
                                              const __m256d yre,
                                              const __m256d yim,
                                              __m256d * __restrict  mre) {

                        register __m256d cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(xre,yre,
                                                          _mm256_mul_pd(xim,yim)));
                        *mre = cvmag;
             }
             
             
               
	          
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3 
#pragma intel optimization_parameter target_arch=AVX
#elif defined (__GNUC__) && (!defined (__INTEL_COMPILER)	|| !defined(__ICC))
#pragma optimize GCC target("avx")
#endif                  
	           static inline
                   void cmagnitude_ymm4c8(   const ymm4c8_t x,
                                              const ymm4c8_t y,
                                              __m256d * __restrict  mre) {

                        register __m256d cvmag;
                        cvmag= _mm256_sqrt_pd(_mm256_fmadd_pd(x.re,y.re,
                                          _mm256_mul_pd(x.im,y.im)));
                        *mre = cvmag;
             }







} // math


} // gms

#endif /*__GMS_COMPLEX_YMM4R8_H__*/