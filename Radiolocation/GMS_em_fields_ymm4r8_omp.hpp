

#ifndef __GMS_EM_FIELDS_YMM4R8_HPP__
#define __GMS_EM_FIELDS_YMM4R8_HPP__ 261020230804

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

    const unsigned int GMS_EM_FIELDS_YMM4R8_OMP_MAJOR = 1U;
    const unsigned int GMS_EM_FIELDS_YMM4R8_OMP_MINOR = 0U;
    const unsigned int GMS_EM_FIELDS_YMM4R8_OMP_MICRO = 0U;
    const unsigned int GMS_EM_FIELDS_YMM4R8_OMP_FULLVER =
      1000U*GMS_EM_FIELDS_YMM4R8_OMP_MAJOR+
      100U*GMS_EM_FIELDS_YMM4R8_OMP_MINOR+
      10U*GMS_EM_FIELDS_YMM4R8_OMP_MICRO;
    const char * const GMS_EM_FIELDS_YMM4R8_OMP_CREATION_DATE = "26-10-2023 08:04 AM +00200 (26 OCT 2023 GMT+2)";
    const char * const GMS_EM_FIELDS_YMM4R8_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EM_FIELDS_YMM4R8_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EM_FIELDS_YMM4R8_OMP_DESCRIPTION   = " Computational ElectroMagnetics related helper routines OpenMP-multithreaded."
                       

}


#include <immintrin.h>
#include <cstdint>
#include <omp.h>
#include "GMS_config.h"
#include "GMS_em_fields_ymm4r8.hpp"


namespace gms {



          namespace radiolocation {
          
          
 
              
               
	      
	      
	      
	          
	      
	      
	         
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void sdotv_ymm4r8_unroll6x_omp( const __m256d * __restrict __ATTR_ALIGN__(32) pv1x,
	                                           const __m256d * __restrict __ATTR_ALIGN__(32) pv1y,
	                                           const __m256d * __restrict __ATTR_ALIGN__(32) pv1z,
	                                           const __m256d * __restrict __ATTR_ALIGN__(32) pv2x,
	                                           const __m256d * __restrict __ATTR_ALIGN__(32) pv2y,
	                                           const __m256d * __restrict __ATTR_ALIGN__(32) pv2z,
	                                           __m256d * __restrict __ATTR_ALIGN__(32) pdtv,
	                                           const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m256d v1x1,v1x2,v1x3,v1x4,v1x5,v1x6;
	                register __m256d v1y1,v1y2,v1y3,v1y4,v1y5,v1y6;
	                register __m256d v1z1,v1z2,v1z3,v1z4,v1z5,v1z6;
	                register __m256d v2x1,v2x2,v2x3,v2x4,v2x5,v2x6;
	                register __m256d v2y1,v2y2,v2y3,v2y4,v2y5,v2y6;
	                register __m256d v2z1,v2z2,v2z3,v2z4,v2z5,v2z6;
	                register __m256d dtv1,dtv2,dtv3,dtv4,dtv5,dtv6;
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                   
	                        v1x1 = pv1x[j];
	                        v2x1 = pv2x[j];
	                        v1y1 = pv1y[j];
	                        v2y1 = pv2y[j];
	                        v1z1 = pv1z[j];
	                        v2z1 = pv2z[j];
	                        dtv1 = sdotv_ymm4r8(v1x,v1y,v1z,
	                                            v2x,v2y,v2z);
	                        pdtv[j] = dtv1;
	                   }
	                   if(n<6) { return;}
	                }  
	                
	                m1 = m+1;
#pragma omp parallel for schedule(dynamic) default(none)                  \
            firstprivate(m1) private(j,v1x1,v1x2,v1x3,v1x4,v1x5,v1x6)     \
                             private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6)       \
                             private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6)       \
                             private(v2x1,v2x2,v2x3,v2x4,v2x5,v2x6)       \
                             private(v2y1,v2y2,v2y3,v2y4,v2y5,v2y6)       \
                             private(v2z1,v2z2,v2z3,v2z4,v2z5,v2z6)       \
                             private(dtv1,dtv2,dtv3,dtv4,dtv5,dtv6)       \
                             shared(pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,n,pdtv)
	                for(j = m1; j != n; j += 6) {
	                     v1x1 = pv1x[j+0];
	                     v2x1 = pv2x[j+0];
	                     v1y1 = pv1y[j+0];
	                     v2y1 = pv2y[j+0];
	                     v1z1 = pv1z[j+0];
	                     v2z1 = pv2z[j+0];
	                     dtv1 = sdotv_ymm4r8(v1x1,v1y1,v1z1,
	                                         v2x1,v2y1,v2z1);
	                     pdtv[j+0] = dtv1;
	                     v1x2 = pv1x[j+1];
	                     v2x2 = pv2x[j+1];
	                     v1y2 = pv1y[j+1];
	                     v2y2 = pv2y[j+1];
	                     v1z2 = pv1z[j+1];
	                     v2z2 = pv2z[j+1];
	                     dtv2 = sdotv_ymm4r8(v1x2,v1y2,v1z2,
	                                         v2x2,v2y2,v2z2);
	                     pdtv[j+1] = dtv2;
	                     v1x3 = pv1x[j+2];
	                     v2x3 = pv2x[j+2];
	                     v1y3 = pv1y[j+2];
	                     v2y3 = pv2y[j+2];
	                     v1z3 = pv1z[j+2];
	                     v2z3 = pv2z[j+2];
	                     dtv3 = sdotv_ymm4r8(v1x3,v1y3,v1z3,
	                                         v2x3,v2y3,v2z3);
	                     pdtv[j+2] = dtv3;
	                     v1x4 = pv1x[j+3];
	                     v2x4 = pv2x[j+3];
	                     v1y4 = pv1y[j+3];
	                     v2y4 = pv2y[j+3];
	                     v1z4 = pv1z[j+3];
	                     v2z4 = pv2z[j+3];
	                     dtv4 = sdotv_ymm4r8(v1x4,v1y4,v1z4,
	                                         v2x4,v2y4,v2z4);
	                     pdtv[j+3] = dtv4;
	                     v1x5 = pv1x[j+4];
	                     v2x5 = pv2x[j+4];
	                     v1y5 = pv1y[j+4];
	                     v2y5 = pv2y[j+4];
	                     v1z5 = pv1z[j+4];
	                     v2z5 = pv2z[j+4];
	                     dtv5 = sdotv_ymm4r8(v1x5,v1y5,v1z5,
	                                         v2x5,v2y5,v2z5);
	                     pdtv[j+4] = dtv5;
	                     v1x6 = pv1x[j+5];
	                     v2x6 = pv2x[j+5];
	                     v1y6 = pv1y[j+5];
	                     v2y6 = pv2y[j+5];
	                     v1z6 = pv1z[j+5];
	                     v2z6 = pv2z[j+5];
	                     dtv6 = sdotv_ymm4r8(v1x6,v1y6,v1z6,
	                                         v2x6,v2y6,v2z6);
	                     pdtv[j+5] = dtv6;
	                   	                    
	                }
	                                        
	      }
	      
	      
	      
	      
	        
	        
	        
	      
	        
	       
	       
	        
	        
	       
	       
	       
	        
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cdotv_ymm4c8_unroll10x_omp(const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm4c8_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                __ATTR_ALIGN__(32) ymm4c8_t v1x1,v1x2,v1x3,v1x4,v1x5,v1x6,v1x7,v1x8,v1x9,v1x10;
	                __ATTR_ALIGN__(32) ymm4c8_t v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10;
	                __ATTR_ALIGN__(32) ymm4c8_t v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10;
	                __ATTR_ALIGN__(32) ymm4c8_t v2x1,v2x2,v2x3,v2x4,v2x5,v2x6,v2x7,v2x8,v2x9,v2x10;
	                __ATTR_ALIGN__(32) ymm4c8_t v2y1,v2y2,v2y3,v2y4,v2y5,v2y6,v2y7,v2y8,v2y9,v2y10;
	                __ATTR_ALIGN__(32) ymm4c8_t v2z1,v2z2,v2z3,v2z4,v2z5,v2z6,v2z7,v2z8,v2z9,v2z10;
	                __ATTR_ALIGN__(32) ymm4c8_t res1,res2,res3,res4,res5,res6,res7,res8,res9,res10;
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
	                       cdotv_ymm4c8(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<10) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                       \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2,v1x3,v1x4,v1x5,v1x6,v1x7,v1x8,v1x9,v1x10) \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10)   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10)   \
                                 private(v2x1,v2x2,v2x3,v2x4,v2x5,v2x6,v2x7,v2x8,v2x9,v2x10)   \
                                 private(v2y1,v2y2,v2y3,v2y4,v2y5,v2y6,v2y7,v2y8,v2y9,v2y10)   \
                                 private(v2z1,v2z2,v2z3,v2z4,v2z5,v2z6,v2z7,v2z8,v2z9,v2z10)   \
                                 private(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)   \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,pres)
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
	                    v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    cdotv_ymm4c8(v1x1,v1y1,v1z1,
	                                  v2x1,v2y1,v2z1,
	                                  res1);
	                    pres[j+0] = res1;
	                    v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    cdotv_ymm4c8(v1x2,v1y2,v1z2,
	                                  v2x2,v2y2,v2z2,
	                                  res2);
	                    pres[j+1] = res2;
	                    v1x3 = pv1x[j+2];
	                    v2x3 = pv2x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v2y3 = pv2y[j+2];
	                    v1z3 = pv1z[j+2];
	                    v2z3 = pv2z[j+2];
	                    cdotv_ymm4c8(v1x3,v1y3,v1z3,
	                                  v2x3,v2y3,v2z3,
	                                  res3);
	                    pres[j+2] = res3;
	                    v1x4 = pv1x[j+3];
	                    v2x4 = pv2x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v2y4 = pv2y[j+3];
	                    v1z4 = pv1z[j+3];
	                    v2z4 = pv2z[j+3];
	                    cdotv_ymm4c8(v1x4,v1y4,v1z4,
	                                  v2x4,v2y4,v2z4,
	                                  res4);
	                    pres[j+3] = res4;
	                    v1x5 = pv1x[j+4];
	                    v2x5 = pv2x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v2y5 = pv2y[j+4];
	                    v1z5 = pv1z[j+4];
	                    v2z5 = pv2z[j+4];
	                    cdotv_ymm4c8(v1x5,v1y5,v1z5,
	                                  v2x5,v2y5,v2z5,
	                                  res5);
	                    pres[j+4] = res5;
	                    v1x6 = pv1x[j+5];
	                    v2x6 = pv2x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v2y6 = pv2y[j+5];
	                    v1z6 = pv1z[j+5];
	                    v2z6 = pv2z[j+5];
	                    cdotv_ymm4c8(v1x6,v1y6,v1z6,
	                                  v2x6,v2y6,v2z6,
	                                  res6);
	                    pres[j+5] = res6;
	                    v1x7 = pv1x[j+6];
	                    v2x7 = pv2x[j+6];
	                    v1y7 = pv1y[j+6];
	                    v2y7 = pv2y[j+6];
	                    v1z7 = pv1z[j+6];
	                    v2z7 = pv2z[j+6];
	                    cdotv_ymm4c8(v1x7,v1y7,v1z7,
	                                  v2x7,v2y7,v2z7,
	                                  res7);
	                    pres[j+6] = res7;
	                    v1x8 = pv1x[j+7];
	                    v2x8 = pv2x[j+7];
	                    v1y8 = pv1y[j+7];
	                    v2y8 = pv2y[j+7];
	                    v1z8 = pv1z[j+7];
	                    v2z8 = pv2z[j+7];
	                    cdotv_ymm4c8(v1x8,v1y8,v1z8,
	                                  v2x8,v2y8,v2z8,
	                                  res8);
	                    pres[j+7] = res8;
	                    v1x9 = pv1x[j+8];
	                    v2x9 = pv2x[j+8];
	                    v1y9 = pv1y[j+8];
	                    v2y9 = pv2y[j+8];
	                    v1z9 = pv1z[j+8];
	                    v2z9 = pv2z[j+8];
	                    cdotv_ymm4c8(v1x9,v1y9,v1z9,
	                                  v2x9,v2y9,v2z9,
	                                  res9);
	                    pres[j+8] = res9;
	                    v1x10 = pv1x[j+9];
	                    v2x10 = pv2x[j+9];
	                    v1y10 = pv1y[j+9];
	                    v2y10 = pv2y[j+9];
	                    v1z10 = pv1z[j+9];
	                    v2z10 = pv2z[j+9];
	                    cdotv_ymm4c8(v1x10,v1y10,v1z10,
	                                  v2x10,v2y10,v2z10,
	                                  res10);
	                    pres[j+9] = res10;
	                
	             }          
	       }
	       
	        
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cdotv_ymm4c8_unroll6x_omp(const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm4c8_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1x1,v1x2,v1x3,v1x4,v1x5,v1x6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1y1,v1y2,v1y3,v1y4,v1y5,v1y6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1z1,v1z2,v1z3,v1z4,v1z5,v1z6;
	                __ATTR_ALIGN__(32) ymm4c8_t v2x1,v2x2,v2x3,v2x4,v2x5,v2x6;
	                __ATTR_ALIGN__(32) ymm4c8_t v2y1,v2y2,v2y3,v2y4,v2y5,v2y6;
	                __ATTR_ALIGN__(32) ymm4c8_t v2z1,v2z2,v2z3,v2z4,v2z5,v2z6;
	                __ATTR_ALIGN__(32) ymm4c8_t res1,res2,res3,res4,res5,res6;
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
	                       cdotv_ymm4c8(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<6) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                  \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2,v1x3,v1x4,v1x5,v1x6) \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6)   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6)   \
                                 private(v2x1,v2x2,v2x3,v2x4,v2x5,v2x6)   \
                                 private(v2y1,v2y2,v2y3,v2y4,v2y5,v2y6)   \
                                 private(v2z1,v2z2,v2z3,v2z4,v2z5,v2z6)   \
                                 private(res1,res2,res3,res4,res5,res6)   \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,pres)
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
	                    v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    cdotv_ymm4c8(v1x1,v1y1,v1z1,
	                                  v2x1,v2y1,v2z1,
	                                  res1);
	                    pres[j+0] = res1;
	                    v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    cdotv_ymm4c8(v1x2,v1y2,v1z2,
	                                  v2x2,v2y2,v2z2,
	                                  res2);
	                    pres[j+1] = res2;
	                    v1x3 = pv1x[j+2];
	                    v2x3 = pv2x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v2y3 = pv2y[j+2];
	                    v1z3 = pv1z[j+2];
	                    v2z3 = pv2z[j+2];
	                    cdotv_ymm4c8(v1x3,v1y3,v1z3,
	                                  v2x3,v2y3,v2z3,
	                                  res3);
	                    pres[j+2] = res3;
	                    v1x4 = pv1x[j+3];
	                    v2x4 = pv2x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v2y4 = pv2y[j+3];
	                    v1z4 = pv1z[j+3];
	                    v2z4 = pv2z[j+3];
	                    cdotv_ymm4c8(v1x4,v1y4,v1z4,
	                                  v2x4,v2y4,v2z4,
	                                  res4);
	                    pres[j+3] = res4;
	                    v1x5 = pv1x[j+4];
	                    v2x5 = pv2x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v2y5 = pv2y[j+4];
	                    v1z5 = pv1z[j+4];
	                    v2z5 = pv2z[j+4];
	                    cdotv_ymm4c8(v1x5,v1y5,v1z5,
	                                  v2x5,v2y5,v2z5,
	                                  res5);
	                    pres[j+4] = res5;
	                    v1x6 = pv1x[j+5];
	                    v2x6 = pv2x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v2y6 = pv2y[j+5];
	                    v1z6 = pv1z[j+5];
	                    v2z6 = pv2z[j+5];
	                    cdotv_ymm4c8(v1x6,v1y6,v1z6,
	                                  v2x6,v2y6,v2z6,
	                                  res6);
	                    pres[j+5] = res6;
	                                 
	             }          
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cdotv_ymm4c8_unroll2x_omp(const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm4c8_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(0>=n,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1x1,v1x2;
	                __ATTR_ALIGN__(32) ymm4c8_t v1y1,v1y2;
	                __ATTR_ALIGN__(32) ymm4c8_t v1z1,v1z2;
	                __ATTR_ALIGN__(32) ymm4c8_t v2x1,v2x2;
	                __ATTR_ALIGN__(32) ymm4c8_t v2y1,v2y2;
	                __ATTR_ALIGN__(32) ymm4c8_t v2z1,v2z2;
	                __ATTR_ALIGN__(32) ymm4c8_t res1,res2;
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
	                       cdotv_ymm4c8(v1x,v1y,v1z,
	                                     v2x,v2y,v2z,
	                                     res);
	                       pres[j] = res;
	                   }
	                   if(n<2) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                  \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2) \
                                 private(v1y1,v1y2)   \
                                 private(v1z1,v1z2)   \
                                 private(v2x1,v2x2)   \
                                 private(v2y1,v2y2)   \
                                 private(v2z1,v2z2)   \
                                 private(res1,res2)   \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,pres)
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
	                    v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    cdotv_ymm4c8(v1x1,v1y1,v1z1,
	                                  v2x1,v2y1,v2z1,
	                                  res1);
	                    pres[j+0] = res1;
	                    v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    cdotv_ymm4c8(v1x2,v1y2,v1z2,
	                                  v2x2,v2y2,v2z2,
	                                  res2);
	                    pres[j+1] = res2;
	                   	                                 
	             }          
	       }
	       
	       
	      
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cnorm_ymm4c8_unroll10x_omp(const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm4c8_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 10;
	                __ATTR_ALIGN__(32) ymm4c8_t v1x1,v1x2,v1x3,v1x4,v1x,v1x6,v1x7,v1x8,v1x9,v1x10;
	                __ATTR_ALIGN__(32) ymm4c8_t v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10;
	                __ATTR_ALIGN__(32) ymm4c8_t v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10;
	                __ATTR_ALIGN__(32) ymm4c8_t vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8,vs9,vs10;
	                int32_t j,m,m1;
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_ymm4c8(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<10) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                       \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2,v1x3,v1x4,v1x,v1x6,v1x7,v1x8,v1x9,v1x10)  \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10)   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10)   \
                                 private(vs1,vs2,vs3,vs4,vs5,vs6,vs7,vs8,vs9,vs10)             \
                                 shared(n,pv1x,pv1y,pv1z,pvs)	       
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
                            v1x1 = pv1x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v1z1 = pv1z[j+0];
	                    vs1  = cnorm_ymm4c8(v1x1,v1y1,v1z1);
	                    pvs[j+0] = vs1;
	                    v1x2 = pv1x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v1z2 = pv1z[j+1];
	                    vs2  = cnorm_ymm4c8(v1x2,v1y2,v1z2);
	                    pvs[j+1] = vs2;
	                    v1x3 = pv1x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v1z3 = pv1z[j+2];
	                    vs3  = cnorm_ymm4c8(v1x3,v1y3,v1z3);
	                    pvs[j+2] = vs3;
	                    v1x4 = pv1x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v1z4 = pv1z[j+3];
	                    vs4  = cnorm_ymm4c8(v1x4,v1y4,v1z4);
	                    pvs[j+3] = vs4;
	                    v1x5 = pv1x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v1z5 = pv1z[j+4];
	                    vs5  = cnorm_ymm4c8(v1x5,v1y5,v1z5);
	                    pvs[j+4] = vs5;
	                    v1x6 = pv1x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v1z6 = pv1z[j+5];
	                    vs6  = cnorm_ymm4c8(v1x6,v1y6,v1z6);
	                    pvs[j+5] = vs6;
	                    v1x7 = pv1x[j+6];
	                    v1y7 = pv1y[j+6];
	                    v1z7 = pv1z[j+6];
	                    vs7  = cnorm_ymm4c8(v1x7,v1y7,v1z7);
	                    pvs[j+6] = vs7;
	                    v1x8 = pv1x[j+7];
	                    v1y8 = pv1y[j+7];
	                    v1z8 = pv1z[j+7];
	                    vs8  = cnorm_ymm4c8(v1x8,v1y8,v1z8);
	                    pvs[j+7] = vs8;
	                    v1x9 = pv1x[j+8];
	                    v1y9 = pv1y[j+8];
	                    v1z9 = pv1z[j+8];
	                    vs9  = cnorm_ymm4c8(v1x9,v1y9,v1z9);
	                    pvs[j+8] = vs9;
	                    v1x10 = pv1x[j+9];
	                    v1y10 = pv1y[j+9];
	                    v1z10 = pv1z[j+9];
	                    vs10  = cnorm_ymm4c8(v1x10,v1y10,v1z10);
	                    pvs[j+9] = vs10;
	              }             
	      }
	      
	      
	         
	      
	      
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cnorm_ymm4c8_unroll6x_omp(const ymm4c8_t   * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm4c8_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1x1,v1x2,v1x3,v1x4,v1x,v1x6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1y1,v1y2,v1y3,v1y4,v1y5,v1y6;
	                __ATTR_ALIGN__(32) ymm4c8_t v1z1,v1z2,v1z3,v1z4,v1z5,v1z6;
	                __ATTR_ALIGN__(32) ymm4c8_t vs1,vs2,vs3,vs4,vs5,vs6;
	                int32_t j,m,m1;
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_ymm4c8(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<6) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                  \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2,v1x3,v1x4,v1x,v1x6)  \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6)   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6)   \
                                 private(vs1,vs2,vs3,vs4,vs5,vs6)        \
                                 shared(n,pv1x,pv1y,pv1z,pvs)	       
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
                            v1x1 = pv1x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v1z1 = pv1z[j+0];
	                    vs1  = cnorm_ymm4c8(v1x1,v1y1,v1z1);
	                    pvs[j+0] = vs1;
	                    v1x2 = pv1x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v1z2 = pv1z[j+1];
	                    vs2  = cnorm_ymm4c8(v1x2,v1y2,v1z2);
	                    pvs[j+1] = vs2;
	                    v1x3 = pv1x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v1z3 = pv1z[j+2];
	                    vs3  = cnorm_ymm4c8(v1x3,v1y3,v1z3);
	                    pvs[j+2] = vs3;
	                    v1x4 = pv1x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v1z4 = pv1z[j+3];
	                    vs4  = cnorm_ymm4c8(v1x4,v1y4,v1z4);
	                    pvs[j+3] = vs4;
	                    v1x5 = pv1x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v1z5 = pv1z[j+4];
	                    vs5  = cnorm_ymm4c8(v1x5,v1y5,v1z5);
	                    pvs[j+4] = vs5;
	                    v1x6 = pv1x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v1z6 = pv1z[j+5];
	                    vs6  = cnorm_ymm4c8(v1x6,v1y6,v1z6);
	                    pvs[j+5] = vs6;
	                 
	              }             
	      }
	      
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void cnorm_ymm4c8_unroll2x_omp(  const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm4c8_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 2;
	                __ATTR_ALIGN__(32) ymm4c8_t v1x1,v1x2;
	                __ATTR_ALIGN__(32) ymm4c8_t v1y1,v1y2;
	                __ATTR_ALIGN__(32) ymm4c8_t v1z1,v1z2;
	                __ATTR_ALIGN__(32) ymm4c8_t vs1,vs2;
	                int32_t j,m,m1;
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       v1x = pv1x[j];
	                       v1y = pv1y[j];
	                       v1z = pv1z[j];
	                       vs  = cnorm_ymm4c8(v1x,v1y,v1z);
	                       pvs[j] = vs;
	                   }
	                   if(n<2) { return;}
	                }                     
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)\
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2)  \
                                 private(v1y1,v1y2)   \
                                 private(v1z1,v1z2)   \
                                 private(vs1,vs2)    \
                                 shared(n,pv1x,pv1y,pv1z,pvs)	       
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
                            v1x1 = pv1x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v1z1 = pv1z[j+0];
	                    vs1  = cnorm_ymm4c8(v1x1,v1y1,v1z1);
	                    pvs[j+0] = vs1;
	                    v1x2 = pv1x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v1z2 = pv1z[j+1];
	                    vs2  = cnorm_ymm4c8(v1x2,v1y2,v1z2);
	                    pvs[j+1] = vs2;
	                   	                 
	              }             
	      }
	      
	      
	                                       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossc_ymm4r8_unroll10x_omp( const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                ymm4c8_t resx1,resx2,resx3,resx4,resx5,resx6,resx6,resx7,resx8,resx9,resx10;
	                ymm4c8_t resy1,resy2,resy3,resy4,resy5,resy6,resy6,resy7,resy8,resy9,resy10;
	                ymm4c8_t resz1,resz2,resz3,resz4,resz5,resz6,resz6,resz7,resz8,resz9,resz10;
	                ymm4c8_t  v1x1,v1x2,v1x3,v1x4,v1x5,v1x6,v1x7,v1x8,v1x9,v1x10;
	                ymm4c8_t  v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10;
	                ymm4c8_t  v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10;
	                ymm4c8_t  v2x1,v2x2,v2x3,v2x4,v2x4,v2x5,v2x6,v2x7,v2x8,v2x9,v2x10;
	                ymm4c8_t  v2y1,v2y2,v2y3,v2y4,v2y4,v2y5,v2y6,v2y7,v2y8,v2y9,v2y10;
	                ymm4c8_t  v2z1,v2z2,v2z3,v2z4,v2z4,v2z5,v2z6,v2z7,v2z8,v2z9,v2z10;   
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
	                       scrossc_ymm4r8(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<10) return;
	                }                  
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                                       \
        firstprivate(m1,PF_DIST) private(j,resx1,resx2,resx3,resx4,resx5,resx6,resx6,resx7,resx8,resx9,resx10) \
                                 private(resy1,resy2,resy3,resy4,resy5,resy6,resy6,resy7,resy8,resy9,resy10)   \
                                 private(resz1,resz2,resz3,resz4,resz5,resz6,resz6,resz7,resz8,resz9,resz10)   \
                                 private(v1x1,v1x2,v1x3,v1x4,v1x5,v1x6,v1x7,v1x8,v1x9,v1x10)                   \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10)                   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10)                   \
                                 private(v2x1,v2x2,v2x3,v2x4,v2x4,v2x5,v2x6,v2x7,v2x8,v2x9,v2x10)              \
                                 private(v2y1,v2y2,v2y3,v2y4,v2y4,v2y5,v2y6,v2y7,v2y8,v2y9,v2y10)              \
                                 private(v2z1,v2z2,v2z3,v2z4,v2z4,v2z5,v2z6,v2z7,v2z8,v2z9,v2z10)              \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,presx,presy,presz)
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
                            v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    scrossc_ymm4r8(v1x1,v1y1,v1z1,
	                                    v2x1,v2y1,v2z1,
	                                    resx1,resy1,resz1);
	                    presx[j+0] = resx1;
	                    presy[j+0] = resy1;
	                    presz[j+0] = resz1;
	                    v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    scrossc_ymm4r8(v1x2,v1y2,v1z2,
	                                    v2x2,v2y2,v2z2,
	                                    resx2,resy2,resz2);
	                    presx[j+1] = resx2;
	                    presy[j+1] = resy2;
	                    presz[j+1] = resz2;
	                    v1x3 = pv1x[j+2];
	                    v2x3 = pv2x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v2y3 = pv2y[j+2];
	                    v1z3 = pv1z[j+2];
	                    v2z3 = pv2z[j+2];
	                    scrossc_ymm4r8(v1x3,v1y3,v1z3,
	                                    v2x3,v2y3,v2z3,
	                                    resx3,resy3,resz3);
	                    presx[j+2] = resx3;
	                    presy[j+2] = resy3;
	                    presz[j+2] = resz3;
	                    v1x4 = pv1x[j+3];
	                    v2x4 = pv2x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v2y4 = pv2y[j+3];
	                    v1z4 = pv1z[j+3];
	                    v2z4 = pv2z[j+3];
	                    scrossc_ymm4r8(v1x4,v1y4,v1z4,
	                                    v2x5,v2y4,v2z4,
	                                    resx4,resy4,resz4);
	                    presx[j+3] = resx4;
	                    presy[j+3] = resy4;
	                    presz[j+3] = resz4;
	                    v1x5 = pv1x[j+4];
	                    v2x5 = pv2x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v2y5 = pv2y[j+4];
	                    v1z5 = pv1z[j+4];
	                    v2z5 = pv2z[j+4];
	                    scrossc_ymm4r8(v1x5,v1y5,v1z5,
	                                    v2x5,v2y5,v2z5,
	                                    resx5,resy5,resz5);
	                    presx[j+4] = resx5;
	                    presy[j+4] = resy5;
	                    presz[j+4] = resz5;
	                    v1x6 = pv1x[j+5];
	                    v2x6 = pv2x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v2y6 = pv2y[j+5];
	                    v1z6 = pv1z[j+5];
	                    v2z6 = pv2z[j+5];
	                    scrossc_ymm4r8(v1x6,v1y6,v1z6,
	                                    v2x6,v2y6,v2z6,
	                                    resx6,resy6,resz6);
	                    presx[j+5] = resx6;
	                    presy[j+5] = resy6;
	                    presz[j+5] = resz6;
	                    v1x7 = pv1x[j+6];
	                    v2x7 = pv2x[j+6];
	                    v1y7 = pv1y[j+6];
	                    v2y7 = pv2y[j+6];
	                    v1z7 = pv1z[j+6];
	                    v2z7 = pv2z[j+6];
	                    scrossc_ymm4r8(v1x7,v1y7,v1z7,
	                                    v2x7,v2y7,v2z7,
	                                    resx7,resy7,resz7);
	                    presx[j+6] = resx7;
	                    presy[j+6] = resy7;
	                    presz[j+6] = resz7;
	                    v1x8 = pv1x[j+7];
	                    v2x8 = pv2x[j+7];
	                    v1y8 = pv1y[j+7];
	                    v2y8 = pv2y[j+7];
	                    v1z8 = pv1z[j+7];
	                    v2z8 = pv2z[j+7];
	                    scrossc_ymm4r8(v1x8,v1y8,v1z8,
	                                    v2x8,v2y8,v2z8,
	                                    resx8,resy8,resz8);
	                    presx[j+7] = resx8;
	                    presy[j+7] = resy8;
	                    presz[j+7] = resz8;
	                    v1x9 = pv1x[j+8];
	                    v2x9 = pv2x[j+8];
	                    v1y9 = pv1y[j+8];
	                    v2y9 = pv2y[j+8];
	                    v1z9 = pv1z[j+8];
	                    v2z9 = pv2z[j+8];
	                    scrossc_ymm4r8(v1x9,v1y9,v1z9,
	                                    v2x9,v2y9,v2z9,
	                                    resx9,resy9,resz9);
	                    presx[j+8] = resx9;
	                    presy[j+8] = resy9;
	                    presz[j+8] = resz9;
	                    v1x10 = pv1x[j+9];
	                    v2x10 = pv2x[j+9];
	                    v1y10 = pv1y[j+9];
	                    v2y10 = pv2y[j+9];
	                    v1z10 = pv1z[j+9];
	                    v2z10 = pv2z[j+9];
	                    scrossc_ymm4r8(v1x10,v1y10,v1z10,
	                                    v2x10,v2y10,v2z10,
	                                    resx10,resy10,resz10);
	                    presx[j+9] = resx10;
	                    presy[j+9] = resy10;
	                    presz[j+9] = resz10;
	                 
	                }          
	        }
	        
	        
	        
	     
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossc_ymm4r8_unroll6x_omp(  const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                ymm4c8_t resx1,resx2,resx3,resx4,resx5,resx6,resx6;
	                ymm4c8_t resy1,resy2,resy3,resy4,resy5,resy6,resy6;
	                ymm4c8_t resz1,resz2,resz3,resz4,resz5,resz6,resz6;
	                ymm4c8_t  v1x1,v1x2,v1x3,v1x4,v1x5,v1x6;
	                ymm4c8_t  v1y1,v1y2,v1y3,v1y4,v1y5,v1y6;
	                ymm4c8_t  v1z1,v1z2,v1z3,v1z4,v1z5,v1z6;
	                ymm4c8_t  v2x1,v2x2,v2x3,v2x4,v2x4,v2x5,v2x6;
	                ymm4c8_t  v2y1,v2y2,v2y3,v2y4,v2y4,v2y5,v2y6;
	                ymm4c8_t  v2z1,v2z2,v2z3,v2z4,v2z4,v2z5,v2z6;   
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
	                       scrossc_ymm4r8(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<6) return;
	                }                  
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                              \
        firstprivate(m1,PF_DIST) private(j,resx1,resx2,resx3,resx4,resx5,resx6,resx6) \
                                 private(resy1,resy2,resy3,resy4,resy5,resy6,resy6)   \
                                 private(resz1,resz2,resz3,resz4,resz5,resz6,resz6)   \
                                 private(v1x1,v1x2,v1x3,v1x4,v1x5,v1x6)               \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6)               \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6)               \
                                 private(v2x1,v2x2,v2x3,v2x4,v2x4,v2x5,v2x6)          \
                                 private(v2y1,v2y2,v2y3,v2y4,v2y4,v2y5,v2y6)          \
                                 private(v2z1,v2z2,v2z3,v2z4,v2z4,v2z5,v2z6)          \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,presx,presy,presz)
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
                            v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    scrossc_ymm4r8(v1x1,v1y1,v1z1,
	                                    v2x1,v2y1,v2z1,
	                                    resx1,resy1,resz1);
	                    presx[j+0] = resx1;
	                    presy[j+0] = resy1;
	                    presz[j+0] = resz1;
	                    v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    scrossc_ymm4r8(v1x2,v1y2,v1z2,
	                                    v2x2,v2y2,v2z2,
	                                    resx2,resy2,resz2);
	                    presx[j+1] = resx2;
	                    presy[j+1] = resy2;
	                    presz[j+1] = resz2;
	                    v1x3 = pv1x[j+2];
	                    v2x3 = pv2x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v2y3 = pv2y[j+2];
	                    v1z3 = pv1z[j+2];
	                    v2z3 = pv2z[j+2];
	                    scrossc_ymm4r8(v1x3,v1y3,v1z3,
	                                    v2x3,v2y3,v2z3,
	                                    resx3,resy3,resz3);
	                    presx[j+2] = resx3;
	                    presy[j+2] = resy3;
	                    presz[j+2] = resz3;
	                    v1x4 = pv1x[j+3];
	                    v2x4 = pv2x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v2y4 = pv2y[j+3];
	                    v1z4 = pv1z[j+3];
	                    v2z4 = pv2z[j+3];
	                    scrossc_ymm4r8(v1x4,v1y4,v1z4,
	                                    v2x4,v2y4,v2z4,
	                                    resx4,resy4,resz4);
	                    presx[j+3] = resx4;
	                    presy[j+3] = resy4;
	                    presz[j+3] = resz4;
	                    v1x5 = pv1x[j+4];
	                    v2x5 = pv2x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v2y5 = pv2y[j+4];
	                    v1z5 = pv1z[j+4];
	                    v2z5 = pv2z[j+4];
	                    scrossc_ymm4r8(v1x5,v1y5,v1z5,
	                                    v2x5,v2y5,v2z5,
	                                    resx5,resy5,resz5);
	                    presx[j+4] = resx5;
	                    presy[j+4] = resy5;
	                    presz[j+4] = resz5;
	                    v1x6 = pv1x[j+5];
	                    v2x6 = pv2x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v2y6 = pv2y[j+5];
	                    v1z6 = pv1z[j+5];
	                    v2z6 = pv2z[j+5];
	                    scrossc_ymm4r8(v1x6,v1y6,v1z6,
	                                    v2x6,v2y6,v2z6,
	                                    resx6,resy6,resz6);
	                    presx[j+5] = resx6;
	                    presy[j+5] = resy6;
	                    presz[j+5] = resz6;
	                   	                 
	                }          
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossc_ymm4r8_unroll2x_omp(  const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm4c8_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm4c8_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                ymm4c8_t resx1,resx2;
	                ymm4c8_t resy1,resy2;
	                ymm4c8_t resz1,resz2;
	                ymm4c8_t  v1x1,v1x2;
	                ymm4c8_t  v1y1,v1y2;
	                ymm4c8_t  v1z1,v1z2;
	                ymm4c8_t  v2x1,v2x2;
	                ymm4c8_t  v2y1,v2y2;
	                ymm4c8_t  v2z1,v2z2;   
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
	                       scrossc_ymm4r8(v1x,v1y,v1z,
	                                       v2x,v2y,v2z,
	                                       resx,resy,resz);
	                       presx[j] = resx;
	                       presy[j] = resy;
	                       presz[j] = resz;
	                   }
	                   if(n<2) return;
	                }                  
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)\
        firstprivate(m1,PF_DIST) private(j,resx1,resx2) \
                                 private(resy1,resy2)   \
                                 private(resz1,resz2)   \
                                 private(v1x1,v1x2)     \
                                 private(v1y1,v1y2)     \
                                 private(v1z1,v1z2)     \
                                 private(v2x1,v2x2)     \
                                 private(v2y1,v2y2)     \
                                 private(v2z1,v2z2)     \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z,presx,presy,presz)
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
                            v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    scrossc_ymm4r8(v1x1,v1y1,v1z1,
	                                    v2x1,v2y1,v2z1,
	                                    resx1,resy1,resz1);
	                    presx[j+0] = resx1;
	                    presy[j+0] = resy1;
	                    presz[j+0] = resz1;
	                    v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    scrossc_ymm4r8(v1x2,v1y2,v1z2,
	                                    v2x2,v2y2,v2z2,
	                                    resx2,resy2,resz2);
	                    presx[j+1] = resx2;
	                    presy[j+1] = resy2;
	                    presz[j+1] = resz2;
	                  
	                }          
	        }
	  
	         
	       
	         
	         
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_ymm4r8_unroll10x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv1y,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              __m256d * __restrict __ATTR_ALIGN__(32) pvcx,
	                                              __m256d * __restrict __ATTR_ALIGN__(32) pvcy,
	                                              __m256d * __restrict __ATTR_ALIGN__(32) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                  register __m256d v1x1,v1x2,v1x3,v1x4,v1x5,v1x6,v1x7,v1x8,v1x9,v1x10;
	                  register __m256d v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10; 
	                  register __m256d v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10;
	                  register __m256d v2x1,v2x2,v2x3,v2x4,v2x5,v2x6,v2x7,v2x8,v2x9,v2x10;
	                  register __m256d v2y1,v2y2,v2y3,v2y4,v2y5,v2y6,v2y7,v2y8,v2y9,v2y10;
	                  register __m256d v2z1,v2z2,v2z3,v2z4,v2z5,v2z6,v2z7,v2z8,v2z9,v2z10;
	                  register __m256d vcx1,vcx2,vcx3,vcx4,vcx5,vcx6,vcx7,vcx8,vcx9,vcx10;
	                  register __m256d vcy1,vcy2,vcy3,vcy4,vcy5,vcy6,vcy7,vcy8,vcy9,vcy10;
	                  register __m256d vcz1,vcz2,vcz3,vcz4,vcz5,vcz6,vcz7,vcz8,vcz9,vcz10;
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
	                          scrossv_ymm4r8(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<10) return;
	                  }                  
	                  
	                  m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                       \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2,v1x3,v1x4,v1x5,v1x6,v1x7,v1x8,v1x9,v1x10) \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6,v1y7,v1y8,v1y9,v1y10)   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6,v1z7,v1z8,v1z9,v1z10)   \
                                 private(v2x1,v2x2,v2x3,v2x4,v2x5,v2x6,v2x7,v2x8,v2x9,v2x10)   \
                                 private(v2y1,v2y2,v2y3,v2y4,v2y5,v2y6,v2y7,v2y8,v2y9,v2y10)   \
                                 private(v2z1,v2z2,v2z3,v2z4,v2z5,v2z6,v2z7,v2z8,v2z9,v2z10)   \
                                 private(vcx1,vcx2,vcx3,vcx4,vcx5,vcx6,vcx7,vcx8,vcx9,vcx10)   \
                                 private(vcy1,vcy2,vcy3,vcy4,vcy5,vcy6,vcy7,vcy8,vcy9,vcy10)   \
                                 private(vcz1,vcz2,vcz3,vcz4,vcz5,vcz6,vcz7,vcz8,vcz9,vcz10)   \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z)
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
                            v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    scrossv_ymm4r8(v1x1,v1y1,v1z1,
	                                    v2x1,v2y1,v2z1,
	                                    &vcx1,&vcy1,&vcz1);
	                    pvcx[j+0] = vcx1;
	                    pvcy[j+0] = vcy1;
	                    pvcz[j+0] = vcz1;
                            v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    scrossv_ymm4r8(v1x2,v1y2,v1z2,
	                                    v2x2,v2y2,v2z2,
	                                    &vcx2,&vcy2,&vcz2);
	                    pvcx[j+1] = vcx2;
	                    pvcy[j+1] = vcy2;
	                    pvcz[j+1] = vcz2;
	                    v1x3 = pv1x[j+2];
	                    v2x3 = pv2x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v2y3 = pv2y[j+2];
	                    v1z3 = pv1z[j+2];
	                    v2z3 = pv2z[j+2];
	                    scrossv_ymm4r8(v1x3,v1y3,v1z3,
	                                    v2x3,v2y3,v2z3,
	                                    &vcx3,&vcy3,&vcz3);
	                    pvcx[j+2] = vcx3;
	                    pvcy[j+2] = vcy3;
	                    pvcz[j+2] = vcz3;   
	                    v1x4 = pv1x[j+3];
	                    v2x4 = pv2x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v2y4 = pv2y[j+3];
	                    v1z4 = pv1z[j+3];
	                    v2z4 = pv2z[j+3];
	                    scrossv_ymm4r8(v1x4,v1y4,v1z4,
	                                    v2x4,v2y4,v2z4,
	                                    &vcx4,&vcy4,&vcz4);
	                    pvcx[j+3] = vcx4;
	                    pvcy[j+3] = vcy4;
	                    pvcz[j+3] = vcz4;
	                    v1x5 = pv1x[j+4];
	                    v2x5 = pv2x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v2y5 = pv2y[j+4];
	                    v1z5 = pv1z[j+4];
	                    v2z5 = pv2z[j+4];
	                    scrossv_ymm4r8(v1x5,v1y5,v1z5,
	                                    v2x5,v2y5,v2z5,
	                                    &vcx5,&vcy5,&vcz5);
	                    pvcx[j+4] = vcx5;
	                    pvcy[j+4] = vcy5;
	                    pvcz[j+4] = vcz5;
	                    v1x6 = pv1x[j+5];
	                    v2x6 = pv2x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v2y6 = pv2y[j+5];
	                    v1z6 = pv1z[j+5];
	                    v2z6 = pv2z[j+5];
	                    scrossv_ymm4r8(v1x6,v1y6,v1z6,
	                                    v2x6,v2y6,v2z6,
	                                    &vcx6,&vcy6,&vcz6);
	                    pvcx[j+5] = vcx6;
	                    pvcy[j+5] = vcy6;
	                    pvcz[j+5] = vcz6;
	                    v1x7 = pv1x[j+6];
	                    v2x7 = pv2x[j+6];
	                    v1y7 = pv1y[j+6];
	                    v2y7 = pv2y[j+6];
	                    v1z7 = pv1z[j+6];
	                    v2z7 = pv2z[j+6];
	                    scrossv_ymm4r8(v1x7,v1y7,v1z7,
	                                    v2x7,v2y7,v2z7,
	                                    &vcx7,&vcy7,&vcz7);
	                    pvcx[j+6] = vcx7;
	                    pvcy[j+6] = vcy7;
	                    pvcz[j+6] = vcz7;
	                    v1x8 = pv1x[j+7];
	                    v2x8 = pv2x[j+7];
	                    v1y8 = pv1y[j+7];
	                    v2y8 = pv2y[j+7];
	                    v1z8 = pv1z[j+7];
	                    v2z8 = pv2z[j+7];
	                    scrossv_ymm4r8(v1x8,v1y8,v1z8,
	                                    v2x8,v2y8,v2z8,
	                                    &vcx8,&vcy8,&vcz8);
	                    pvcx[j+7] = vcx8;
	                    pvcy[j+7] = vcy8;
	                    pvcz[j+7] = vcz8;  
	                    v1x9 = pv1x[j+8];
	                    v2x9 = pv2x[j+8];
	                    v1y9 = pv1y[j+8];
	                    v2y9 = pv2y[j+8];
	                    v1z9 = pv1z[j+8];
	                    v2z9 = pv2z[j+8];
	                    scrossv_ymm4r8(v1x9,v1y9,v1z9,
	                                    v2x9,v2y9,v2z9,
	                                    &vcx9,&vcy9,&vcz9);
	                    pvcx[j+8] = vcx9;
	                    pvcy[j+8] = vcy9;
	                    pvcz[j+8] = vcz9;   
	                    v1x10 = pv1x[j+9];
	                    v2x10 = pv2x[j+9];
	                    v1y10 = pv1y[j+9];
	                    v2y10 = pv2y[j+9];
	                    v1z10 = pv1z[j+9];
	                    v2z10 = pv2z[j+9];
	                    scrossv_ymm4r8(v1x10,v1y10,v1z10,
	                                    v2x10,v2y10,v2z10,
	                                    &vcx10,&vcy10,&vcz10);
	                    pvcx[j+9] = vcx10;
	                    pvcy[j+9] = vcy10;
	                    pvcz[j+9] = vcz10;
	                
	               }             
	         }
	         
	         
	         
	         
	       
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_ymm4r8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m256d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m256d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m256d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m256d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m256d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m256d * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m256d * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m256d * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                  register __m256d v1x1,v1x2,v1x3,v1x4,v1x5,v1x6;
	                  register __m256d v1y1,v1y2,v1y3,v1y4,v1y5,v1y6; 
	                  register __m256d v1z1,v1z2,v1z3,v1z4,v1z5,v1z6;
	                  register __m256d v2x1,v2x2,v2x3,v2x4,v2x5,v2x6;
	                  register __m256d v2y1,v2y2,v2y3,v2y4,v2y5,v2y6;
	                  register __m256d v2z1,v2z2,v2z3,v2z4,v2z5,v2z6;
	                  register __m256d vcx1,vcx2,vcx3,vcx4,vcx5,vcx6;
	                  register __m256d vcy1,vcy2,vcy3,vcy4,vcy5,vcy6;
	                  register __m256d vcz1,vcz2,vcz3,vcz4,vcz5,vcz6;
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
	                          scrossv_ymm4r8(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<6) return;
	                  }                  
	                  
	                  m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                  \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2,v1x3,v1x4,v1x5,v1x6) \
                                 private(v1y1,v1y2,v1y3,v1y4,v1y5,v1y6)   \
                                 private(v1z1,v1z2,v1z3,v1z4,v1z5,v1z6)   \
                                 private(v2x1,v2x2,v2x3,v2x4,v2x5,v2x6)   \
                                 private(v2y1,v2y2,v2y3,v2y4,v2y5,v2y6)   \
                                 private(v2z1,v2z2,v2z3,v2z4,v2z5,v2z6)   \
                                 private(vcx1,vcx2,vcx3,vcx4,vcx5,vcx6)   \
                                 private(vcy1,vcy2,vcy3,vcy4,vcy5,vcy6)   \
                                 private(vcz1,vcz2,vcz3,vcz4,vcz5,vcz6)   \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z)
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
                            v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    scrossv_ymm4r8(v1x1,v1y1,v1z1,
	                                    v2x1,v2y1,v2z1,
	                                    &vcx1,&vcy1,&vcz1);
	                    pvcx[j+0] = vcx1;
	                    pvcy[j+0] = vcy1;
	                    pvcz[j+0] = vcz1;
                            v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    scrossv_ymm4r8(v1x2,v1y2,v1z2,
	                                    v2x2,v2y2,v2z2,
	                                    &vcx2,&vcy2,&vcz2);
	                    pvcx[j+1] = vcx2;
	                    pvcy[j+1] = vcy2;
	                    pvcz[j+1] = vcz2;
	                    v1x3 = pv1x[j+2];
	                    v2x3 = pv2x[j+2];
	                    v1y3 = pv1y[j+2];
	                    v2y3 = pv2y[j+2];
	                    v1z3 = pv1z[j+2];
	                    v2z3 = pv2z[j+2];
	                    scrossv_ymm4r8(v1x3,v1y3,v1z3,
	                                    v2x3,v2y3,v2z3,
	                                    &vcx3,&vcy3,&vcz3);
	                    pvcx[j+2] = vcx3;
	                    pvcy[j+2] = vcy3;
	                    pvcz[j+2] = vcz3;   
	                    v1x4 = pv1x[j+3];
	                    v2x4 = pv2x[j+3];
	                    v1y4 = pv1y[j+3];
	                    v2y4 = pv2y[j+3];
	                    v1z4 = pv1z[j+3];
	                    v2z4 = pv2z[j+3];
	                    scrossv_ymm4r8(v1x4,v1y4,v1z4,
	                                    v2x4,v2y4,v2z4,
	                                    &vcx4,&vcy4,&vcz4);
	                    pvcx[j+3] = vcx4;
	                    pvcy[j+3] = vcy4;
	                    pvcz[j+3] = vcz4;
	                    v1x5 = pv1x[j+4];
	                    v2x5 = pv2x[j+4];
	                    v1y5 = pv1y[j+4];
	                    v2y5 = pv2y[j+4];
	                    v1z5 = pv1z[j+4];
	                    v2z5 = pv2z[j+4];
	                    scrossv_ymm4r8(v1x5,v1y5,v1z5,
	                                    v2x5,v2y5,v2z5,
	                                    &vcx5,&vcy5,&vcz5);
	                    pvcx[j+4] = vcx5;
	                    pvcy[j+4] = vcy5;
	                    pvcz[j+4] = vcz5;
	                    v1x6 = pv1x[j+5];
	                    v2x6 = pv2x[j+5];
	                    v1y6 = pv1y[j+5];
	                    v2y6 = pv2y[j+5];
	                    v1z6 = pv1z[j+5];
	                    v2z6 = pv2z[j+5];
	                    scrossv_ymm4r8(v1x6,v1y6,v1z6,
	                                    v2x6,v2y6,v2z6,
	                                    &vcx6,&vcy6,&vcz6);
	                    pvcx[j+5] = vcx6;
	                    pvcy[j+5] = vcy6;
	                    pvcz[j+5] = vcz6;
	                                        
	               }             
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_ymm4r8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv1y,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const __m256d * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              __m256d * __restrict __ATTR_ALIGN__(32) pvcx,
	                                              __m256d * __restrict __ATTR_ALIGN__(32) pvcy,
	                                              __m256d * __restrict __ATTR_ALIGN__(32) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST) {
	                                          
	                  if(__builtin_expect(n<=0,0)) { return;}
	                  if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                  register __m256d v1x1,v1x2;
	                  register __m256d v1y1,v1y2; 
	                  register __m256d v1z1,v1z2;
	                  register __m256d v2x1,v2x2;
	                  register __m256d v2y1,v2y2;
	                  register __m256d v2z1,v2z2;
	                  register __m256d vcx1,vcx2;
	                  register __m256d vcy1,vcy2;
	                  register __m256d vcz1,vcz2;
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
	                          scrossv_ymm4r8(v1x,v1y,v1z,
	                                          v2x,v2y,v2z,
	                                          &vcx,&vcy,&vcz);
	                          pvcx[j] = vcx;
	                          pvcy[j] = vcy;
	                          pvcz[j] = vcz;
	                     }
	                     if(n<2) return;
	                  }                  
	                  
	                  m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                  \
        firstprivate(m1,PF_DIST) private(j,v1x1,v1x2) \
                                 private(v1y1,v1y2)   \
                                 private(v1z1,v1z2)   \
                                 private(v2x1,v2x2)   \
                                 private(v2y1,v2y2)   \
                                 private(v2z1,v2z2)   \
                                 private(vcx1,vcx2)   \
                                 private(vcy1,vcy2)   \
                                 private(vcz1,vcz2)   \
                                 shared(n,pv1x,pv1y,pv1z,pv2x,pv2y,pv2z)
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
                            v1x1 = pv1x[j+0];
	                    v2x1 = pv2x[j+0];
	                    v1y1 = pv1y[j+0];
	                    v2y1 = pv2y[j+0];
	                    v1z1 = pv1z[j+0];
	                    v2z1 = pv2z[j+0];
	                    scrossv_ymm4r8(v1x1,v1y1,v1z1,
	                                    v2x1,v2y1,v2z1,
	                                    &vcx1,&vcy1,&vcz1);
	                    pvcx[j+0] = vcx1;
	                    pvcy[j+0] = vcy1;
	                    pvcz[j+0] = vcz1;
                            v1x2 = pv1x[j+1];
	                    v2x2 = pv2x[j+1];
	                    v1y2 = pv1y[j+1];
	                    v2y2 = pv2y[j+1];
	                    v1z2 = pv1z[j+1];
	                    v2z2 = pv2z[j+1];
	                    scrossv_ymm4r8(v1x2,v1y2,v1z2,
	                                    v2x2,v2y2,v2z2,
	                                    &vcx2,&vcy2,&vcz2);
	                    pvcx[j+1] = vcx2;
	                    pvcy[j+1] = vcy2;
	                    pvcz[j+1] = vcz2;
	                  	                                        
	               }             
	         }
	         
	         
	         
	         
    
	         
	         
	         
	          //! Direction Vector spherical [theta,phi] (SIMD data-types)
	         
	              
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_ymm4r8_unroll10x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 10;
	                __m256d tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10;
	                __m256d phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10;
	                __m256d dvx1,dvx2,dvx3,dvx4,dvx5,dvx6,dvx7,dvx8,dvx9,dvx10;
	                __m256d dvy1,dvy2,dvy3,dvy4,dvy5,dvy6,dvy7,dvy8,dvy9,dvy10;
	                __m256d dvz1,dvz2,dvz3,dvz4,dvz5,dvz6,dvz7,dvz8,dvz9,dvz10;
	               int32_t j,m,m1;
	               
	               m = n%10;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_ymm4r8(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<10) {return;}
	               } 
	               
	               m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                          \
        firstprivate(m1,PF_DIST) private(j,tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10)    \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10)      \
                                 private(dvx1,dvx2,dvx3,dvx4,dvx5,dvx6,dvx7,dvx8,dvx9,dvx10)      \
                                 private(dvy1,dvy2,dvy3,dvy4,dvy5,dvy6,dvy7,dvy8,dvy9,dvy10)      \
                                 private(dvz1,dvz2,dvz3,dvz4,dvz5,dvz6,dvz7,dvz8,dvz9,dvz10)      \
                                 shared(n,ptht,pphi,pdvx,pdvy,pdvz)
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
                            tht1 = ptht[j+0];
	                    phi1 = pphi[j+0];
	                    dir_vec_ymm4r8(tht1,phi1,
	                                    &dvx1,&dvy1,&dvz1);
	                    pdvx[j+0] = dvx1;
	                    pdvy[j+0] = dvy1;
	                    pdvz[j+0] = dvz1;
	                    tht2 = ptht[j+1];
	                    phi2 = pphi[j+1];
	                    dir_vec_ymm4r8(tht2,phi2,
	                                    &dvx2,&dvy2,&dvz2);
	                    pdvx[j+1] = dvx2;
	                    pdvy[j+1] = dvy2;
	                    pdvz[j+1] = dvz2;
	                    tht3 = ptht[j+2];
	                    phi3 = pphi[j+2];
	                    dir_vec_ymm4r8(tht3,phi3,
	                                    &dvx3,&dvy3,&dvz3);
	                    pdvx[j+2] = dvx3;
	                    pdvy[j+2] = dvy3;
	                    pdvz[j+2] = dvz3;
	                    tht4 = ptht[j+3];
	                    phi4 = pphi[j+3];
	                    dir_vec_ymm4r8(tht4,phi4,
	                                    &dvx4,&dvy4,&dvz4);
	                    pdvx[j+3] = dvx4;
	                    pdvy[j+3] = dvy4;
	                    pdvz[j+3] = dvz4;
	                    tht5 = ptht[j+4];
	                    phi5 = pphi[j+4];
	                    dir_vec_ymm4r8(tht5,phi5,
	                                    &dvx5,&dvy5,&dvz5);
	                    pdvx[j+4] = dvx5;
	                    pdvy[j+4] = dvy5;
	                    pdvz[j+4] = dvz5;
	                    tht6 = ptht[j+5];
	                    phi6 = pphi[j+5];
	                    dir_vec_ymm4r8(tht6,phi6,
	                                    &dvx6,&dvy6,&dvz6);
	                    pdvx[j+5] = dvx6;
	                    pdvy[j+5] = dvy6;
	                    pdvz[j+5] = dvz6;
	                    tht7 = ptht[j+6];
	                    phi7 = pphi[j+6];
	                    dir_vec_ymm4r8(tht7,phi7,
	                                    &dvx7,&dvy7,&dvz7);
	                    pdvx[j+6] = dvx7;
	                    pdvy[j+6] = dvy7;
	                    pdvz[j+6] = dvz7;
	                    tht8 = ptht[j+7];
	                    phi8 = pphi[j+7];
	                    dir_vec_ymm4r8(tht8,phi8,
	                                    &dvx8,&dvy8,&dvz8);
	                    pdvx[j+7] = dvx8;
	                    pdvy[j+7] = dvy8;
	                    pdvz[j+7] = dvz8;
	                    tht9 = ptht[j+8];
	                    phi9 = pphi[j+8];
	                    dir_vec_ymm4r8(tht9,phi9,
	                                    &dvx9,&dvy9,&dvz9);
	                    pdvx[j+8] = dvx9;
	                    pdvy[j+8] = dvy9;
	                    pdvz[j+8] = dvz9;
	                    tht10 = ptht[j+9];
	                    phi10 = pphi[j+9];
	                    dir_vec_ymm4r8(tht10,phi10,
	                                    &dvx10,&dvy10,&dvz10);
	                    pdvx[j+9] = dvx10;
	                    pdvy[j+9] = dvy10;
	                    pdvz[j+9] = dvz10;
	                  
	               }                                  
	       }
	       
	       
	        
	       
	     
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_ymm4r8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 6;
	                __m256d tht1,tht2,tht3,tht4,tht5,tht6;
	                __m256d phi1,phi2,phi3,phi4,phi5,phi6;
	                __m256d dvx1,dvx2,dvx3,dvx4,dvx5,dvx6;
	                __m256d dvy1,dvy2,dvy3,dvy4,dvy5,dvy6;
	                __m256d dvz1,dvz2,dvz3,dvz4,dvz5,dvz6;
	               int32_t j,m,m1;
	               
	               m = n%6;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_ymm4r8(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<6) {return;}
	               } 
	               
	               m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                     \
        firstprivate(m1,PF_DIST) private(j,tht1,tht2,tht3,tht4,tht5,tht6)    \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6)      \
                                 private(dvx1,dvx2,dvx3,dvx4,dvx5,dvx6)      \
                                 private(dvy1,dvy2,dvy3,dvy4,dvy5,dvy6)      \
                                 private(dvz1,dvz2,dvz3,dvz4,dvz5,dvz6)     \
                                 shared(n,ptht,pphi,pdvx,pdvy,pdvz)
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
                            tht1 = ptht[j+0];
	                    phi1 = pphi[j+0];
	                    dir_vec_ymm4r8(tht1,phi1,
	                                    &dvx1,&dvy1,&dvz1);
	                    pdvx[j+0] = dvx1;
	                    pdvy[j+0] = dvy1;
	                    pdvz[j+0] = dvz1;
	                    tht2 = ptht[j+1];
	                    phi2 = pphi[j+1];
	                    dir_vec_ymm4r8(tht2,phi2,
	                                    &dvx2,&dvy2,&dvz2);
	                    pdvx[j+1] = dvx2;
	                    pdvy[j+1] = dvy2;
	                    pdvz[j+1] = dvz2;
	                    tht3 = ptht[j+2];
	                    phi3 = pphi[j+2];
	                    dir_vec_ymm4r8(tht3,phi3,
	                                    &dvx3,&dvy3,&dvz3);
	                    pdvx[j+2] = dvx3;
	                    pdvy[j+2] = dvy3;
	                    pdvz[j+2] = dvz3;
	                    tht4 = ptht[j+3];
	                    phi4 = pphi[j+3];
	                    dir_vec_ymm4r8(tht4,phi4,
	                                    &dvx4,&dvy4,&dvz4);
	                    pdvx[j+3] = dvx4;
	                    pdvy[j+3] = dvy4;
	                    pdvz[j+3] = dvz4;
	                    tht5 = ptht[j+4];
	                    phi5 = pphi[j+4];
	                    dir_vec_ymm4r8(tht5,phi5,
	                                    &dvx5,&dvy5,&dvz5);
	                    pdvx[j+4] = dvx5;
	                    pdvy[j+4] = dvy5;
	                    pdvz[j+4] = dvz5;
	                    tht6 = ptht[j+5];
	                    phi6 = pphi[j+5];
	                    dir_vec_ymm4r8(tht6,phi6,
	                                    &dvx6,&dvy6,&dvz6);
	                    pdvx[j+5] = dvx6;
	                    pdvy[j+5] = dvy6;
	                    pdvz[j+5] = dvz6;
	                  
	            }                                  
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_ymm4r8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	               if(__builtin_expect(n<=0,0)) {return;}
	               if(__builtin_expect(PF_DIST<=0,0) PF_DIST = 2;
	                __m256d tht1,tht2;
	                __m256d phi1,phi2;
	                __m256d dvx1,dvx2;
	                __m256d dvy1,dvy2;
	                __m256d dvz1,dvz2;
	               int32_t j,m,m1;
	               
	               m = n%2;
	               if(m!=0) {
	                  for(j = 0; j != m; ++j) {
	                      tht = ptht[j];
	                      phi = pphi[j];
	                      dir_vec_ymm4r8(tht,phi,
	                                      &dvx,&dvy,&dvz);
	                      pdvx[j] = dvx;
	                      pdvy[j] = dvy;
	                      pdvz[j] = dvz;
	                  }
	                  if(n<2) {return;}
	               } 
	               
	               m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                     \
        firstprivate(m1,PF_DIST) private(j,tht1,tht2)    \
                                 private(phi1,phi2)      \
                                 private(dvx1,dvx2)      \
                                 private(dvy1,dvy2)      \
                                 private(dvz1,dvz2)     \
                                 shared(n,ptht,pphi,pdvx,pdvy,pdvz)
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
                            tht1 = ptht[j+0];
	                    phi1 = pphi[j+0];
	                    dir_vec_ymm4r8(tht1,phi1,
	                                    &dvx1,&dvy1,&dvz1);
	                    pdvx[j+0] = dvx1;
	                    pdvy[j+0] = dvy1;
	                    pdvz[j+0] = dvz1;
	                    tht2 = ptht[j+1];
	                    phi2 = pphi[j+1];
	                    dir_vec_ymm4r8(tht2,phi2,
	                                    &dvx2,&dvy2,&dvz2);
	                    pdvx[j+1] = dvx2;
	                    pdvy[j+1] = dvy2;
	                    pdvz[j+1] = dvz2;
	                  	                  
	            }                                  
	       }
	       
	                     	       	     
	        	        
	         //! Polarization Vector of plane-wave propagating into direction computed by
                 //! dir_vector_xmmxrx (SIMD data-types)
                                                   	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_ymm4r8_unroll10x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                          const __m256d * __restrict __ATTR_ALIGN__(32) ppsi,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) ppvx,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) ppvy,
	                                          __m256d * __restrict __ATTR_ALIGN__(32) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                register __m256d tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10;
	                register __m256d phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10;
	                register __m256d psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8,psi9,psi10;
	                register __m256d pvx1,pvx2,pvx3,pvx4,pvx5,pvx6,pvx7,pvx8,pvx9,pvx10;
	                register __m256d pvy1,pvy2,pvy3,pvy4,pvy5,pvy6,pvy7,pvy8,pvy9,pvy10;
	                register __m256d pvz1,pvz2,pvz3,pvz4,pvz5,pvz6,pvz7,pvz8,pvz9,pvz10;
	                int32_t j,m,m1;    
	                
	                m = n%10;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_ymm4r8(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<10) {return;}
	                }                    
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                                        \
        firstprivate(m1,PF_DIST) private(j,tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10)  \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10)    \
                                 private(psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8,psi9,psi10)    \
                                 private(pvx1,pvx2,pvx3,pvx4,pvx5,pvx6,pvx7,pvx8,pvx9,pvx10)    \
                                 private(pvy1,pvy2,pvy3,pvy4,pvy5,pvy6,pvy7,pvy8,pvy9,pvy10)    \
                                 private(pvz1,pvz2,pvz3,pvz4,pvz5,pvz6,pvz7,pvz8,pvz9,pvz10)    \
                                 shared(n,ptht,pphi,ppsi,ppvx,ppvy,ppvz)
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
                            tht1 = ptht[j+0];
	                    phi1 = pphi[j+0];
	                    psi1 = ppsi[j+0];
	                    pol_vec_ymm4r8(tht1,phi1,psi1,
	                                    &pvx1,&pvy1,&pvz1);
	                    ppvx[j+0] = pvx1;
	                    ppvy[j+0] = pvy1;
	                    ppvz[j+0] = pvz1;  
	                    tht2 = ptht[j+1];
	                    phi2 = pphi[j+1];
	                    psi2 = ppsi[j+1];
	                    pol_vec_ymm4r8(tht2,phi2,psi2,
	                                    &pvx2,&pvy2,&pvz2);
	                    ppvx[j+1] = pvx2;
	                    ppvy[j+1] = pvy2;
	                    ppvz[j+1] = pvz2;   
	                    tht3 = ptht[j+2];
	                    phi3 = pphi[j+2];
	                    psi3 = ppsi[j+2];
	                    pol_vec_ymm4r8(tht3,phi3,psi3,
	                                    &pvx3,&pvy3,&pvz3);
	                    ppvx[j+2] = pvx3;
	                    ppvy[j+2] = pvy3;
	                    ppvz[j+2] = pvz3;
	                    tht4 = ptht[j+3];
	                    phi4 = pphi[j+3];
	                    psi4 = ppsi[j+3];
	                    pol_vec_ymm4r8(tht4,phi4,psi4,
	                                    &pvx4,&pvy4,&pvz4);
	                    ppvx[j+3] = pvx4;
	                    ppvy[j+3] = pvy4;
	                    ppvz[j+3] = pvz4;  
	                    tht5 = ptht[j+4];
	                    phi5 = pphi[j+4];
	                    psi5 = ppsi[j+4];
	                    pol_vec_ymm4r8(tht5,phi5,psi5,
	                                    &pvx5,&pvy5,&pvz5);
	                    ppvx[j+4] = pvx5;
	                    ppvy[j+4] = pvy5;
	                    ppvz[j+4] = pvz5;
	                    tht6 = ptht[j+5];
	                    phi6 = pphi[j+5];
	                    psi6 = ppsi[j+5];
	                    pol_vec_ymm4r8(tht6,phi6,psi6,
	                                    &pvx6,&pvy6,&pvz6);
	                    ppvx[j+5] = pvx6;
	                    ppvy[j+5] = pvy6;
	                    ppvz[j+5] = pvz6;  
	                    tht7 = ptht[j+6];
	                    phi7 = pphi[j+6];
	                    psi7 = ppsi[j+6];
	                    pol_vec_ymm4r8(tht7,phi7,psi7,
	                                    &pvx7,&pvy7,&pvz7);
	                    ppvx[j+6] = pvx7;
	                    ppvy[j+6] = pvy7;
	                    ppvz[j+6] = pvz7;  
	                    tht8 = ptht[j+7];
	                    phi8 = pphi[j+7];
	                    psi8 = ppsi[j+7];
	                    pol_vec_ymm4r8(tht8,phi8,psi8,
	                                    &pvx8,&pvy8,&pvz8);
	                    ppvx[j+7] = pvx8;
	                    ppvy[j+7] = pvy8;
	                    ppvz[j+7] = pvz8;  
	                    tht9 = ptht[j+8];
	                    phi9 = pphi[j+8];
	                    psi9 = ppsi[j+8];
	                    pol_vec_ymm4r8(tht9,phi9,psi9,
	                                    &pvx9,&pvy9,&pvz9);
	                    ppvx[j+8] = pvx9;
	                    ppvy[j+8] = pvy9;
	                    ppvz[j+8] = pvz9;  
	                    tht10 = ptht[j+9];
	                    phi10 = pphi[j+9];
	                    psi10 = ppsi[j+9];
	                    pol_vec_ymm4r8(tht10,phi10,psi10,
	                                    &pvx10,&pvy10,&pvz10);
	                    ppvx[j+9] = pvx10;
	                    ppvy[j+9] = pvy10;
	                    ppvz[j+9] = pvz10;   
	                  
	                }            
	      }
	      
	      
	      
	      
	    
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_ymm4r8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m256d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m256d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m256d * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m256d * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m256d * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                register __m256d tht1,tht2,tht3,tht4,tht5,tht6;
	                register __m256d phi1,phi2,phi3,phi4,phi5,phi6;
	                register __m256d psi1,psi2,psi3,psi4,psi5,psi6;
	                register __m256d pvx1,pvx2,pvx3,pvx4,pvx5,pvx6;
	                register __m256d pvy1,pvy2,pvy3,pvy4,pvy5,pvy6;
	                register __m256d pvz1,pvz2,pvz3,pvz4,pvz5,pvz6;
	                int32_t j,m,m1;    
	                
	                m = n%6;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_ymm4r8(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<6) {return;}
	                }                    
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                   \
        firstprivate(m1,PF_DIST) private(j,tht1,tht2,tht3,tht4,tht5,tht6)  \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6)    \
                                 private(psi1,psi2,psi3,psi4,psi5,psi6)    \
                                 private(pvx1,pvx2,pvx3,pvx4,pvx5,pvx6)    \
                                 private(pvy1,pvy2,pvy3,pvy4,pvy5,pvy6)    \
                                 private(pvz1,pvz2,pvz3,pvz4,pvz5,pvz6)    \
                                 shared(n,ptht,pphi,ppsi,ppvx,ppvy,ppvz)
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
                            tht1 = ptht[j+0];
	                    phi1 = pphi[j+0];
	                    psi1 = ppsi[j+0];
	                    pol_vec_ymm4r8(tht1,phi1,psi1,
	                                    &pvx1,&pvy1,&pvz1);
	                    ppvx[j+0] = pvx1;
	                    ppvy[j+0] = pvy1;
	                    ppvz[j+0] = pvz1;  
	                    tht2 = ptht[j+1];
	                    phi2 = pphi[j+1];
	                    psi2 = ppsi[j+1];
	                    pol_vec_ymm4r8(tht2,phi2,psi2,
	                                    &pvx2,&pvy2,&pvz2);
	                    ppvx[j+1] = pvx2;
	                    ppvy[j+1] = pvy2;
	                    ppvz[j+1] = pvz2;   
	                    tht3 = ptht[j+2];
	                    phi3 = pphi[j+2];
	                    psi3 = ppsi[j+2];
	                    pol_vec_ymm4r8(tht3,phi3,psi3,
	                                    &pvx3,&pvy3,&pvz3);
	                    ppvx[j+2] = pvx3;
	                    ppvy[j+2] = pvy3;
	                    ppvz[j+2] = pvz3;
	                    tht4 = ptht[j+3];
	                    phi4 = pphi[j+3];
	                    psi4 = ppsi[j+3];
	                    pol_vec_ymm4r8(tht4,phi4,psi4,
	                                    &pvx4,&pvy4,&pvz4);
	                    ppvx[j+3] = pvx4;
	                    ppvy[j+3] = pvy4;
	                    ppvz[j+3] = pvz4;  
	                    tht5 = ptht[j+4];
	                    phi5 = pphi[j+4];
	                    psi5 = ppsi[j+4];
	                    pol_vec_ymm4r8(tht5,phi5,psi5,
	                                    &pvx5,&pvy5,&pvz5);
	                    ppvx[j+4] = pvx5;
	                    ppvy[j+4] = pvy5;
	                    ppvz[j+4] = pvz5;
	                    tht6 = ptht[j+5];
	                    phi6 = pphi[j+5];
	                    psi6 = ppsi[j+5];
	                    pol_vec_ymm4r8(tht6,phi6,psi6,
	                                    &pvx6,&pvy6,&pvz6);
	                    ppvx[j+5] = pvx6;
	                    ppvy[j+5] = pvy6;
	                    ppvz[j+5] = pvz6;  
	            }            
	      }
	      
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_ymm4r8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                            const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                            const __m256d * __restrict __ATTR_ALIGN__(32) ppsi,
	                                            __m256d * __restrict __ATTR_ALIGN__(32) ppvx,
	                                            __m256d * __restrict __ATTR_ALIGN__(32) ppvy,
	                                            __m256d * __restrict __ATTR_ALIGN__(32) ppvz,
	                                            const int32_t n,
	                                            int32_t & PF_DIST) {
	                                          
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __m256d tht1,tht2;
	                __m256d phi1,phi2;
	                __m256d psi1,psi2;
	                __m256d pvx1,pvx2;
	                __m256d pvy1,pvy2;
	                __m256d pvz1,pvz2;
	                int32_t j,m,m1;    
	                
	                m = n%2;
	                if(m!=0) {
	                   for(j = 0; j != m; ++j) {
	                       tht = ptht[j];
	                       phi = pphi[j];
	                       psi = ppsi[j];
	                       pol_vec_ymm4r8(tht,phi,psi,
	                                       &pvx,&pvy,&pvz);
	                       ppvx[j] = pvx;
	                       ppvy[j] = pvy;
	                       ppvz[j] = pvz;
	                   }
	                   if(n<2) {return;}
	                }                    
	                
	                m1 = m+1;
#pragma omp parallel for schedule(runtime) default(none)                   \
        firstprivate(m1,PF_DIST) private(j,tht1,tht2)  \
                                 private(phi1,phi2)    \
                                 private(psi1,psi2)    \
                                 private(pvx1,pvx2)    \
                                 private(pvy1,pvy2)    \
                                 private(pvz1,pvz2)    \
                                 shared(n,ptht,pphi,ppsi,ppvx,ppvy,ppvz)
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
                            tht1 = ptht[j+0];
	                    phi1 = pphi[j+0];
	                    psi1 = ppsi[j+0];
	                    pol_vec_ymm4r8(tht1,phi1,psi1,
	                                    &pvx1,&pvy1,&pvz1);
	                    ppvx[j+0] = pvx1;
	                    ppvy[j+0] = pvy1;
	                    ppvz[j+0] = pvz1;  
	                    tht2 = ptht[j+1];
	                    phi2 = pphi[j+1];
	                    psi2 = ppsi[j+1];
	                    pol_vec_ymm4r8(tht2,phi2,psi2,
	                                    &pvx2,&pvy2,&pvz2);
	                    ppvx[j+1] = pvx2;
	                    ppvy[j+1] = pvy2;
	                    ppvz[j+1] = pvz2;   
	            }            
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_ymm4c8_unroll10x_omp( const __m256d * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                __ATTR_ALIGN__(32) ymm4c8_t k1,k2,k3,k4,k5,k6,k7,k8,k9,k10;
	                __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2,H_x3,H_x4,H_x5,H_x6,H_x7,H_x8,H_x9,H_x10;
	                __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2,H_y3,H_y4,H_y5,H_y6,H_y7,H_y8,H_y9,H_y10;
	                __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2,H_z3,H_z4,H_z5,H_z6,H_z7,H_z8,H_z9,H_z10;
	                 __m256d vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6,vpolx7,vpolx8,vpolx9,vpolx10;
	                 __m256d vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6,vpoly7,vpoly8,vpoly9,vpoly10;
	                 __m256d vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6,vpolz7,vpolz8,vpolz9,vpolz10;
	                 __m256d vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6,vdirx7,vdirx8,vdirx9,vdirx10;
	                 __m256d vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6,vdiry7,vdiry8,vdiry9,vdiry10;
	                 __m256d vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6,vdirz7,vdirz8,vdirz9,vdirz10;
	                 __m256d vrx1,vrx2,vrx3,vrx4,vrx5,vrx6,vrx7,vrx8,vrx9,vrx10;
	                 __m256d vry1,vry2,vry3,vry4,vry5,vry6,vry7,vry8,vry9,vry10;
	                 __m256d vrz1,vrz2,vrz3,vrz4,vrz5,vrz6,vrz7,vrz8,vrz9,vrz10;
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
	                       H_XYZ_VP_ymm4c8(vpolx,vpoly,vpolz,
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
#pragma omp parallel for schedule(runtime) default(none)                                                         \
        firstprivate(m1,PF_DIST) private(j,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)                                       \
                                 private(H_x1,H_x2,H_x3,H_x4,H_x5,H_x6,H_x7,H_x8,H_x9,H_x10)                     \
                                 private(H_y1,H_y2,H_y3,H_y4,H_y5,H_y6,H_y7,H_y8,H_y9,H_y10)                     \
                                 private(H_z1,H_z2,H_z3,H_z4,H_z5,H_z6,H_z7,H_z8,H_z9,H_z10)                     \
                                 private(vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6,vpolx7,vpolx8,vpolx9,vpolx10) \
                                 private(vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6,vpoly7,vpoly8,vpoly9,vpoly10) \
                                 private(vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6,vpolz7,vpolz8,vpolz9,vpolz10) \
                                 private(vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6,vdirx7,vdirx8,vdirx9,vdirx10) \
                                 private(vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6,vdiry7,vdiry8,vdiry9,vdiry10) \
                                 private(vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6,vdirz7,vdirz8,vdirz9,vdirz10) \
                                 private(vrx1,vrx2,vrx3,vrx4,vrx5,vrx6,vrx7,vrx8,vrx9,vrx10)                     \
                                 private(vry1,vry2,vry3,vry4,vry5,vry6,vry7,vry8,vry9,vry10)                     \
                                 private(vrz1,vrz2,vrz3,vrz4,vrz5,vrz6,vrz7,vrz8,vrz9,vrz10)                     \
                                 shared(n,pvpolx,pvpoly,pvpolz,pvdirx,pvdiry,pvdirz,pvrx,pvry,pvrz,pk)           \
                                 shared(pH_x,pH_y,pH_z)
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
                            vpolx1 = pvpolx[j+0];
	                    vpoly1 = pvpoly[j+0];
	                    vpolz1 = pvpolz[j+0];
	                    vdirx1 = pvdirx[j+0];
	                    vdiry1 = pvdiry[j+0];
	                    vdirz1 = pvdirz[j+0];
	                    vrx1   = pvrx[j+0];
	                    vry1   = pvry[j+0];
	                    vrz1   = pvrz[j+0];
	                    k1     = pk[j+0];
	                    H_XYZ_VP_ymm4c8(vpolx1,vpoly1,vpolz1,
	                                     vdirx1,vdiry1,vdirz1,
	                                     vrx1,vry1,vrz1,
	                                     H_x1,H_y1,H_z1);
	                    pH_x[j+0] = H_x1;
	                    pH_y[j+0] = H_y1;
	                    pH_z[j+0] = H_z1;   
                            vpolx2 = pvpolx[j+1];
	                    vpoly2 = pvpoly[j+1];
	                    vpolz2 = pvpolz[j+1];
	                    vdirx2 = pvdirx[j+1];
	                    vdiry2 = pvdiry[j+1];
	                    vdirz2 = pvdirz[j+1];
	                    vrx2   = pvrx[j+1];
	                    vry2   = pvry[j+1];
	                    vrz2   = pvrz[j+1];
	                    k2     = pk[j+1];
	                    H_XYZ_VP_ymm4c8(vpolx2,vpoly2,vpolz2,
	                                     vdirx2,vdiry2,vdirz2,
	                                     vrx2,vry2,vrz2,
	                                     H_x2,H_y2,H_z2);
	                    pH_x[j+1] = H_x2;
	                    pH_y[j+1] = H_y2;
	                    pH_z[j+1] = H_z2;
	                    vpolx3 = pvpolx[j+2];
	                    vpoly3 = pvpoly[j+2];
	                    vpolz3 = pvpolz[j+2];
	                    vdirx3 = pvdirx[j+2];
	                    vdiry3 = pvdiry[j+2];
	                    vdirz3 = pvdirz[j+2];
	                    vrx3   = pvrx[j+2];
	                    vry3   = pvry[j+2];
	                    vrz3   = pvrz[j+2];
	                    k3     = pk[j+2];
	                    H_XYZ_VP_ymm4c8(vpolx3,vpoly3,vpolz3,
	                                     vdirx3,vdiry3,vdirz3,
	                                     vrx3,vry3,vrz3,
	                                     H_x3,H_y3,H_z3);
	                    pH_x[j+2] = H_x3;
	                    pH_y[j+2] = H_y3;
	                    pH_z[j+2] = H_z3;   
	                    vpolx4 = pvpolx[j+3];
	                    vpoly4 = pvpoly[j+3];
	                    vpolz4 = pvpolz[j+3];
	                    vdirx4 = pvdirx[j+3];
	                    vdiry4 = pvdiry[j+3];
	                    vdirz4 = pvdirz[j+3];
	                    vrx4   = pvrx[j+3];
	                    vry4   = pvry[j+3];
	                    vrz4   = pvrz[j+3];
	                    k4     = pk[j+3];
	                    H_XYZ_VP_ymm4c8(vpolx4,vpoly4,vpolz4,
	                                     vdirx4,vdiry4,vdirz4,
	                                     vrx4,vry4,vrz4,
	                                     H_x4,H_y4,H_z4);
	                    pH_x[j+3] = H_x4;
	                    pH_y[j+3] = H_y4;
	                    pH_z[j+3] = H_z4;  
	                    vpolx5 = pvpolx[j+4];
	                    vpoly5 = pvpoly[j+4];
	                    vpolz5 = pvpolz[j+4];
	                    vdirx5 = pvdirx[j+4];
	                    vdiry5 = pvdiry[j+4];
	                    vdirz5 = pvdirz[j+4];
	                    vrx5   = pvrx[j+4];
	                    vry5   = pvry[j+4];
	                    vrz   = pvrz[j+4];
	                    k5     = pk[j+4];
	                    H_XYZ_VP_ymm4c8(vpolx5,vpoly5,vpolz5,
	                                     vdirx5,vdiry5,vdirz5,
	                                     vrx5,vry5,vrz5,
	                                     H_x5,H_y5,H_z5);
	                    pH_x[j+4] = H_x5;
	                    pH_y[j+4] = H_y5;
	                    pH_z[j+4] = H_z5;
	                    vpolx6 = pvpolx[j+5];
	                    vpoly6 = pvpoly[j+5];
	                    vpolz6 = pvpolz[j+5];
	                    vdirx6 = pvdirx[j+5];
	                    vdiry6 = pvdiry[j+5];
	                    vdirz6 = pvdirz[j+5];
	                    vrx6   = pvrx[j+5];
	                    vry6   = pvry[j+5];
	                    vrz6   = pvrz[j+5];
	                    k6     = pk[j+5];
	                    H_XYZ_VP_ymm4c8(vpolx6,vpoly6,vpolz6,
	                                     vdirx6,vdiry6,vdirz6,
	                                     vrx6,vry6,vrz6,
	                                     H_x6,H_y6,H_z6);
	                    pH_x[j+5] = H_x6;
	                    pH_y[j+5] = H_y6;
	                    pH_z[j+5] = H_z6;  
	                    vpolx7 = pvpolx[j+6];
	                    vpoly7 = pvpoly[j+6];
	                    vpolz7 = pvpolz[j+6];
	                    vdirx7 = pvdirx[j+6];
	                    vdiry7 = pvdiry[j+6];
	                    vdirz7 = pvdirz[j+6];
	                    vrx7   = pvrx[j+6];
	                    vry7   = pvry[j+6];
	                    vrz7   = pvrz[j+6];
	                    k7     = pk[j+6];
	                    H_XYZ_VP_ymm4c8(vpolx7,vpoly7,vpolz7,
	                                     vdirx7,vdiry7,vdirz7,
	                                     vrx7,vry7,vrz7,
	                                     H_x7,H_y7,H_z7);
	                    pH_x[j+6] = H_x7;
	                    pH_y[j+6] = H_y7;
	                    pH_z[j+6] = H_z7;
	                    vpolx8 = pvpolx[j+7];
	                    vpoly8 = pvpoly[j+7];
	                    vpolz8 = pvpolz[j+7];
	                    vdirx8 = pvdirx[j+7];
	                    vdiry8 = pvdiry[j+7];
	                    vdirz8 = pvdirz[j+7];
	                    vrx8   = pvrx[j+7];
	                    vry8   = pvry[j+7];
	                    vrz8   = pvrz[j+7];
	                    k8     = pk[j+7];
	                    H_XYZ_VP_ymm4c8(vpolx8,vpoly8,vpolz8,
	                                     vdirx8,vdiry8,vdirz8,
	                                     vrx8,vry8,vrz8,
	                                     H_x8,H_y8,H_z8);
	                    pH_x[j+7] = H_x8;
	                    pH_y[j+7] = H_y8;
	                    pH_z[j+7] = H_z8;  
	                    vpolx9 = pvpolx[j+8];
	                    vpoly9 = pvpoly[j+8];
	                    vpolz9 = pvpolz[j+8];
	                    vdirx9 = pvdirx[j+8];
	                    vdiry9 = pvdiry[j+8];
	                    vdirz9 = pvdirz[j+8];
	                    vrx9   = pvrx[j+8];
	                    vry9   = pvry[j+8];
	                    vrz9   = pvrz[j+8];
	                    k9     = pk[j+8];
	                    H_XYZ_VP_ymm4c8(vpolx9,vpoly9,vpolz9,
	                                     vdirx9,vdiry9,vdirz9,
	                                     vrx9,vry9,vrz9,
	                                     H_x9,H_y9,H_z9);
	                    pH_x[j+8] = H_x9;
	                    pH_y[j+8] = H_y9;
	                    pH_z[j+8] = H_z9;
	                    vpolx10 = pvpolx[j+9];
	                    vpoly10 = pvpoly[j+9];
	                    vpolz10 = pvpolz[j+9];
	                    vdirx10 = pvdirx[j+9];
	                    vdiry10 = pvdiry[j+9];
	                    vdirz10 = pvdirz[j+9];
	                    vrx10   = pvrx[j+9];
	                    vry10   = pvry[j+9];
	                    vrz10   = pvrz[j+9];
	                    k10     = pk[j+9];
	                    H_XYZ_VP_ymm4c8(vpolx10,vpoly10,vpolz10,
	                                     vdirx10,vdiry10,vdirz10,
	                                     vrx10,vry10,vrz10,
	                                     H_x10,H_y10,H_z10);
	                    pH_x[j+9] = H_x10;
	                    pH_y[j+9] = H_y10;
	                    pH_z[j+9] = H_z10;         
	                }                
	      }
	      
	      
	      
	    
	      
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_ymm4c8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                __ATTR_ALIGN__(32) ymm4c8_t k1,k2,k3,k4,k5,k6;
	                __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2,H_x3,H_x4,H_x5,H_x6;
	                __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2,H_y3,H_y4,H_y5,H_y6;
	                __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2,H_z3,H_z4,H_z5,H_z6;
	                 __m256d vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6;
	                 __m256d vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6;
	                 __m256d vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6;
	                 __m256d vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6;
	                 __m256d vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6;
	                 __m256d vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6;
	                 __m256d vrx1,vrx2,vrx3,vrx4,vrx5,vrx6;
	                 __m256d vry1,vry2,vry3,vry4,vry5,vry6;
	                 __m256d vrz1,vrz2,vrz3,vrz4,vrz5,vrz6;
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
	                       H_XYZ_VP_ymm4c8(vpolx,vpoly,vpolz,
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
#pragma omp parallel for schedule(runtime) default(none)                            \
        firstprivate(m1,PF_DIST) private(j,k1,k2,k3,k4,k5,k6)                       \
                                 private(H_x1,H_x2,H_x3,H_x4,H_x5,H_x6)             \
                                 private(H_y1,H_y2,H_y3,H_y4,H_y5,H_y6)             \
                                 private(H_z1,H_z2,H_z3,H_z4,H_z5,H_z6)             \
                                 private(vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6) \
                                 private(vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6) \
                                 private(vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6) \
                                 private(vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6) \
                                 private(vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6) \
                                 private(vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6) \
                                 private(vrx1,vrx2,vrx3,vrx4,vrx5,vrx6)             \
                                 private(vry1,vry2,vry3,vry4,vry5,vry6)             \
                                 private(vrz1,vrz2,vrz3,vrz4,vrz5,vrz6)             \
                                 shared(n,pvpolx,pvpoly,pvpolz,pvdirx,pvdiry,pvdirz,pvrx,pvry,pvrz,pk)\
                                 shared(pH_x,pH_y,pH_z)
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
                            vpolx1 = pvpolx[j+0];
	                    vpoly1 = pvpoly[j+0];
	                    vpolz1 = pvpolz[j+0];
	                    vdirx1 = pvdirx[j+0];
	                    vdiry1 = pvdiry[j+0];
	                    vdirz1 = pvdirz[j+0];
	                    vrx1   = pvrx[j+0];
	                    vry1   = pvry[j+0];
	                    vrz1   = pvrz[j+0];
	                    k1     = pk[j+0];
	                    H_XYZ_VP_ymm4c8(vpolx1,vpoly1,vpolz1,
	                                     vdirx1,vdiry1,vdirz1,
	                                     vrx1,vry1,vrz1,
	                                     H_x1,H_y1,H_z1);
	                    pH_x[j+0] = H_x1;
	                    pH_y[j+0] = H_y1;
	                    pH_z[j+0] = H_z1;   
                            vpolx2 = pvpolx[j+1];
	                    vpoly2 = pvpoly[j+1];
	                    vpolz2 = pvpolz[j+1];
	                    vdirx2 = pvdirx[j+1];
	                    vdiry2 = pvdiry[j+1];
	                    vdirz2 = pvdirz[j+1];
	                    vrx2   = pvrx[j+1];
	                    vry2   = pvry[j+1];
	                    vrz2   = pvrz[j+1];
	                    k2     = pk[j+1];
	                    H_XYZ_VP_ymm4c8(vpolx2,vpoly2,vpolz2,
	                                     vdirx2,vdiry2,vdirz2,
	                                     vrx2,vry2,vrz2,
	                                     H_x2,H_y2,H_z2);
	                    pH_x[j+1] = H_x2;
	                    pH_y[j+1] = H_y2;
	                    pH_z[j+1] = H_z2;
	                    vpolx3 = pvpolx[j+2];
	                    vpoly3 = pvpoly[j+2];
	                    vpolz3 = pvpolz[j+2];
	                    vdirx3 = pvdirx[j+2];
	                    vdiry3 = pvdiry[j+2];
	                    vdirz3 = pvdirz[j+2];
	                    vrx3   = pvrx[j+2];
	                    vry3   = pvry[j+2];
	                    vrz3   = pvrz[j+2];
	                    k3     = pk[j+2];
	                    H_XYZ_VP_ymm4c8(vpolx3,vpoly3,vpolz3,
	                                     vdirx3,vdiry3,vdirz3,
	                                     vrx3,vry3,vrz3,
	                                     H_x3,H_y3,H_z3);
	                    pH_x[j+2] = H_x3;
	                    pH_y[j+2] = H_y3;
	                    pH_z[j+2] = H_z3;   
	                    vpolx4 = pvpolx[j+3];
	                    vpoly4 = pvpoly[j+3];
	                    vpolz4 = pvpolz[j+3];
	                    vdirx4 = pvdirx[j+3];
	                    vdiry4 = pvdiry[j+3];
	                    vdirz4 = pvdirz[j+3];
	                    vrx4   = pvrx[j+3];
	                    vry4   = pvry[j+3];
	                    vrz4   = pvrz[j+3];
	                    k4     = pk[j+3];
	                    H_XYZ_VP_ymm4c8(vpolx4,vpoly4,vpolz4,
	                                     vdirx4,vdiry4,vdirz4,
	                                     vrx4,vry4,vrz4,
	                                     H_x4,H_y4,H_z4);
	                    pH_x[j+3] = H_x4;
	                    pH_y[j+3] = H_y4;
	                    pH_z[j+3] = H_z4;  
	                    vpolx5 = pvpolx[j+4];
	                    vpoly5 = pvpoly[j+4];
	                    vpolz5 = pvpolz[j+4];
	                    vdirx5 = pvdirx[j+4];
	                    vdiry5 = pvdiry[j+4];
	                    vdirz5 = pvdirz[j+4];
	                    vrx5   = pvrx[j+4];
	                    vry5   = pvry[j+4];
	                    vrz5   = pvrz[j+4];
	                    k5     = pk[j+4];
	                    H_XYZ_VP_ymm4c8(vpolx5,vpoly5,vpolz5,
	                                     vdirx5,vdiry5,vdirz5,
	                                     vrx5,vry5,vrz5,
	                                     H_x5,H_y5,H_z5);
	                    pH_x[j+4] = H_x5;
	                    pH_y[j+4] = H_y5;
	                    pH_z[j+4] = H_z5;
	                    vpolx6 = pvpolx[j+5];
	                    vpoly6 = pvpoly[j+5];
	                    vpolz6 = pvpolz[j+5];
	                    vdirx6 = pvdirx[j+5];
	                    vdiry6 = pvdiry[j+5];
	                    vdirz6 = pvdirz[j+5];
	                    vrx6   = pvrx[j+5];
	                    vry6   = pvry[j+5];
	                    vrz6   = pvrz[j+5];
	                    k6     = pk[j+5];
	                    H_XYZ_VP_ymm4c8(vpolx6,vpoly6,vpolz6,
	                                     vdirx6,vdiry6,vdirz6,
	                                     vrx6,vry6,vrz6,
	                                     H_x6,H_y6,H_z6);
	                    pH_x[j+5] = H_x6;
	                    pH_y[j+5] = H_y6;
	                    pH_z[j+5] = H_z6;  
	                  
	                }                
	      }
	      
	      
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_ymm4c8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH) {
	                                           
	                                           
	                if(__builtin_expect(n<=0,0)) {return;}
	                if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                __ATTR_ALIGN__(32) ymm4c8_t k1,k2;
	                __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2;
	                __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2;
	                __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2;
	                 __m256d vpolx1,vpolx2;
	                 __m256d vpoly1,vpoly2;
	                 __m256d vpolz1,vpolz2;
	                 __m256d vdirx1,vdirx2;
	                 __m256d vdiry1,vdiry2;
	                 __m256d vdirz1,vdirz2;
	                 __m256d vrx1,vrx2;
	                 __m256d vry1,vry2;
	                 __m256d vrz1,vrz2;
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
	                       H_XYZ_VP_ymm4c8(vpolx,vpoly,vpolz,
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
#pragma omp parallel for schedule(runtime) default(none)        \
        firstprivate(m1,PF_DIST) private(j,k1,k2)               \
                                 private(H_x1,H_x2)             \
                                 private(H_y1,H_y2)             \
                                 private(H_z1,H_z2)             \
                                 private(vpolx1,vpolx2) \
                                 private(vpoly1,vpoly2) \
                                 private(vpolz1,vpolz2) \
                                 private(vdirx1,vdirx2) \
                                 private(vdiry1,vdiry2) \
                                 private(vdirz1,vdirz2) \
                                 private(vrx1,vrx2)             \
                                 private(vry1,vry2)             \
                                 private(vrz1,vrz2)             \
                                 shared(n,pvpolx,pvpoly,pvpolz,pvdirx,pvdiry,pvdirz,pvrx,pvry,pvrz,pk)\
                                 shared(pH_x,pH_y,pH_z)
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
                            vpolx1 = pvpolx[j+0];
	                    vpoly1 = pvpoly[j+0];
	                    vpolz1 = pvpolz[j+0];
	                    vdirx1 = pvdirx[j+0];
	                    vdiry1 = pvdiry[j+0];
	                    vdirz1 = pvdirz[j+0];
	                    vrx1   = pvrx[j+0];
	                    vry1   = pvry[j+0];
	                    vrz1   = pvrz[j+0];
	                    k1     = pk[j+0];
	                    H_XYZ_VP_ymm4c8(vpolx1,vpoly1,vpolz1,
	                                     vdirx1,vdiry1,vdirz1,
	                                     vrx1,vry1,vrz1,
	                                     H_x1,H_y1,H_z1);
	                    pH_x[j+0] = H_x1;
	                    pH_y[j+0] = H_y1;
	                    pH_z[j+0] = H_z1;   
                            vpolx2 = pvpolx[j+1];
	                    vpoly2 = pvpoly[j+1];
	                    vpolz2 = pvpolz[j+1];
	                    vdirx2 = pvdirx[j+1];
	                    vdiry2 = pvdiry[j+1];
	                    vdirz2 = pvdirz[j+1];
	                    vrx2   = pvrx[j+1];
	                    vry2   = pvry[j+1];
	                    vrz2   = pvrz[j+1];
	                    k2     = pk[j+1];
	                    H_XYZ_VP_ymm4c8(vpolx2,vpoly2,vpolz2,
	                                     vdirx2,vdiry2,vdirz2,
	                                     vrx2,vry2,vrz2,
	                                     H_x2,H_y2,H_z2);
	                    pH_x[j+1] = H_x2;
	                    pH_y[j+1] = H_y2;
	                    pH_z[j+1] = H_z2;
	                           
	                }                
	      }
	      
	      
	        
	        /*
	             
     ! Magnetic Field (SIMD data-types) [plane-wave], polarization 'vpol' of
     !  wave-vector argument:  vdir*k at sixteen points 'r'.
	        */
	        
	        
	        
	     
	 
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_ymm4c8_unroll10x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                 __ATTR_ALIGN__(32) ymm4c8_t k1,k2,k3,k4,k5,k6,k7,k8,k9,k10;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2,B_x3,B_x4,B_x5,B_x6,B_x7,B_x8,B_x9,B_x10;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2,B_y3,B_y4,B_y5,B_y6,B_y7,B_y8,B_y9,B_y10;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2,B_z3,B_z4,B_z5,B_z6,B_z7,B_z8,B_z9,B_z10;
	                  __m256d vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6,vpolx7,vpolx8,vpolx9,vpolx10;
	                  __m256d vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6,vpoly8,vpoly9,vpoly9,vpoly10;
	                  __m256d vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6,vpolz7,vpolz8,vpolz9,vpolz10;
	                  __m256d vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6,vdirx7,vdirx8,vdirx9,vdirx10;
	                  __m256d vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6,vdiry7,vdiry8,vdiry9,vdiry10;
	                  __m256d vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6,vdirz7,vdirz8,vdirz9,vdirz10;
	                  __m256d vrx1,vrx2,vrx3,vrx4,vrx5,vrx6,vrx7,vrx8,vrx9,vrx10;
	                  __m256d vry1,vry2,vry3,vry4,vry5,vry6,vry7,vry8,vry9,vry10;
	                  __m256d vrz1,vrz2,vrz3,vrz4,vrz5,vrz6,vrz7,vrz8,vrz9,vrz10;
	                  __m256d omg1,omg2,omg3,omg4,omg5,omg6,omg7,omg8,omg9,omg10;
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
	                        B_XYZ_VP_ymm4c8(vpolx,vpoly,vpolz,
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
#pragma omp parallel for schedule(runtime) default(none)                                                         \
        firstprivate(m1,PF_DIST) private(j,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)                                       \
                                 private(B_x1,B_x2,B_x3,B_x4,B_x5,B_x6,B_x7,B_x8,B_x9,B_x10)                     \
                                 private(B_y1,B_y2,B_y3,B_y4,B_y5,B_y6,B_y7,B_y8,B_y9,B_y10)                     \
                                 private(B_z1,B_z2,B_z3,B_z4,B_z5,B_z6,B_z7,B_z8,B_z9,B_z10)                     \
                                 private(vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6,vpolx7,vpolx8,vpolx9,vpolx10) \
                                 private(vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6,vpoly8,vpoly9,vpoly9,vpoly10) \
                                 private(vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6,vpolz7,vpolz8,vpolz9,vpolz10) \
                                 private(vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6,vdirx7,vdirx8,vdirx9,vdirx10) \
                                 private(vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6,vdiry7,vdiry8,vdiry9,vdiry10) \
                                 private(vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6,vdirz7,vdirz8,vdirz9,vdirz10) \
                                 private(vrx1,vrx2,vrx3,vrx4,vrx5,vrx6,vrx7,vrx8,vrx9,vrx10)                     \
                                 private(vry1,vry2,vry3,vry4,vry5,vry6,vry7,vry8,vry9,vry10)                     \
                                 private(vrz1,vrz2,vrz3,vrz4,vrz5,vrz6,vrz7,vrz8,vrz9,vrz10)                     \
                                 private(omg1,omg2,omg3,omg4,omg5,omg6,omg7,omg8,omg9,omg10)                     \
                                 shared(n,pvpolx,pvpoly,pvpolz,pvdirx,pvdiry,pvdirz)                             \
                                 shared(pvrx,pvry,pvrz,pomega,pk,pB_x,pB_y,pB_z)
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
                                vpolx1 = pvpolx[j+0];
	                        vpoly1 = pvpoly[j+0];
	                        vpolz1 = pvpolz[j+0];
	                        vdirx1 = pvdirx[j+0];
	                        vdiry1 = pvdiry[j+0];
	                        vdirz1 = pvdirz[j+0];
	                        vrx1   = pvrx[j+0];
	                        vry1   = pvry[j+0];
	                        vrz1   = pvrz[j+0];
	                        omg1   = pomega[j+0];
	                        k1     = pk[j+0];
	                        B_XYZ_ymm4c8(vpolx1,vpoly1,vpolz1,
	                                      vdirx1,vdiry1,vdirz1,
	                                      k1,omg1,vrx1,vry1,vrz1,
	                                      B_x1,B_y1,B_z1);
	                        pB_x[j+0] = B_x1;
	                        pB_y[j+0] = B_y1;
	                        pB_z[j+0] = B_z1;    
	                        vpolx2 = pvpolx[j+1];
	                        vpoly2 = pvpoly[j+1];
	                        vpolz2 = pvpolz[j+1];
	                        vdirx2 = pvdirx[j+1];
	                        vdiry2 = pvdiry[j+1];
	                        vdirz2 = pvdirz[j+1];
	                        vrx2   = pvrx[j+1];
	                        vry2   = pvry[j+1];
	                        vrz2   = pvrz[j+1];
	                        omg2   = pomega[j+1];
	                        k2     = pk[j+1];
	                        B_XYZ_ymm4c8(vpolx2,vpoly2,vpolz2,
	                                      vdirx2,vdiry2,vdirz2,
	                                      k2,omg2,vrx2,vry2,vrz2,
	                                      B_x2,B_y2,B_z2);
	                        pB_x[j+1] = B_x2;
	                        pB_y[j+1] = B_y2;
	                        pB_z[j+1] = B_z2;
	                        vpolx3 = pvpolx[j+2];
	                        vpoly3 = pvpoly[j+2];
	                        vpolz3 = pvpolz[j+2];
	                        vdirx3 = pvdirx[j+2];
	                        vdiry3 = pvdiry[j+2];
	                        vdirz3 = pvdirz[j+2];
	                        vrx3   = pvrx[j+2];
	                        vry3   = pvry[j+2];
	                        vrz3   = pvrz[j+2];
	                        omg3   = pomega[j+2];
	                        k3     = pk[j+2];
	                        B_XYZ_ymm4c8(vpolx3,vpoly3,vpolz3,
	                                      vdirx3,vdiry3,vdirz3,
	                                      k3,omg3,vrx3,vry3,vrz3,
	                                      B_x3,B_y3,B_z3);
	                        pB_x[j+2] = B_x3;
	                        pB_y[j+2] = B_y3;
	                        pB_z[j+2] = B_z3;
	                        vpolx4 = pvpolx[j+3];
	                        vpoly4 = pvpoly[j+3];
	                        vpolz4 = pvpolz[j+3];
	                        vdirx4 = pvdirx[j+3];
	                        vdiry4 = pvdiry[j+3];
	                        vdirz4 = pvdirz[j+3];
	                        vrx4   = pvrx[j+3];
	                        vry4   = pvry[j+3];
	                        vrz4   = pvrz[j+3];
	                        omg4   = pomega[j+3];
	                        k4     = pk[j+3];
	                        B_XYZ_ymm4c8(vpolx4,vpoly4,vpolz4,
	                                      vdirx4,vdiry4,vdirz4,
	                                      k4,omg4,vrx4,vry4,vrz4,
	                                      B_x4,B_y4,B_z4);
	                        pB_x[j+3] = B_x4;
	                        pB_y[j+3] = B_y4;
	                        pB_z[j+3] = B_z4;
	                        vpolx5 = pvpolx[j+4];
	                        vpoly5 = pvpoly[j+4];
	                        vpolz5 = pvpolz[j+4];
	                        vdirx5 = pvdirx[j+4];
	                        vdiry5 = pvdiry[j+4];
	                        vdirz5 = pvdirz[j+4];
	                        vrx5   = pvrx[j+4];
	                        vry5   = pvry[j+4];
	                        vrz5   = pvrz[j+4];
	                        omg5   = pomega[j+4];
	                        k5     = pk[j+4];
	                        B_XYZ_ymm4c8(vpolx5,vpoly5,vpolz5,
	                                      vdirx5,vdiry5,vdirz5,
	                                      k5,omg5,vrx5,vry5,vrz5,
	                                      B_x5,B_y5,B_z5);
	                        pB_x[j+4] = B_x5;
	                        pB_y[j+4] = B_y5;
	                        pB_z[j+4] = B_z5;
	                        vpolx6 = pvpolx[j+5];
	                        vpoly6 = pvpoly[j+5];
	                        vpolz6 = pvpolz[j+5];
	                        vdirx6 = pvdirx[j+5];
	                        vdiry6 = pvdiry[j+5];
	                        vdirz6 = pvdirz[j+5];
	                        vrx6   = pvrx[j+5];
	                        vry6   = pvry[j+5];
	                        vrz6   = pvrz[j+5];
	                        omg6   = pomega[j+5];
	                        k6     = pk[j+5];
	                        B_XYZ_ymm4c8(vpolx6,vpoly6,vpolz6,
	                                      vdirx6,vdiry6,vdirz6,
	                                      k6,omg6,vrx6,vry6,vrz6,
	                                      B_x6,B_y6,B_z6);
	                        pB_x[j+5] = B_x6;
	                        pB_y[j+5] = B_y6;
	                        pB_z[j+5] = B_z6;
	                        vpolx7 = pvpolx[j+6];
	                        vpoly7 = pvpoly[j+6];
	                        vpolz7 = pvpolz[j+6];
	                        vdirx7 = pvdirx[j+6];
	                        vdiry7 = pvdiry[j+6];
	                        vdirz7 = pvdirz[j+6];
	                        vrx7   = pvrx[j+6];
	                        vry7   = pvry[j+6];
	                        vrz7   = pvrz[j+6];
	                        omg7   = pomega[j+6];
	                        k7    = pk[j+6];
	                        B_XYZ_ymm4c8(vpolx7,vpoly7,vpolz7,
	                                      vdirx7,vdiry7,vdirz7,
	                                      k7,omg7,vrx7,vry7,vrz7,
	                                      B_x7,B_y7,B_z7);
	                        pB_x[j+6] = B_x7;
	                        pB_y[j+6] = B_y7;
	                        pB_z[j+6] = B_z7;  
	                        vpolx8 = pvpolx[j+7];
	                        vpoly8 = pvpoly[j+7];
	                        vpolz8 = pvpolz[j+7];
	                        vdirx8 = pvdirx[j+7];
	                        vdiry8 = pvdiry[j+7];
	                        vdirz8 = pvdirz[j+7];
	                        vrx8   = pvrx[j+7];
	                        vry8   = pvry[j+7];
	                        vrz8   = pvrz[j+7];
	                        omg8   = pomega[j+7];
	                        k8     = pk[j+7];
	                        B_XYZ_ymm4c8(vpolx8,vpoly8,vpolz8,
	                                      vdirx8,vdiry8,vdirz8,
	                                      k8,omg8,vrx8,vry8,vrz8,
	                                      B_x8,B_y8,B_z8);
	                        pB_x[j+7] = B_x8;
	                        pB_y[j+7] = B_y8;
	                        pB_z[j+7] = B_z8;   
	                        vpolx9 = pvpolx[j+8];
	                        vpoly9 = pvpoly[j+8];
	                        vpolz9 = pvpolz[j+8];
	                        vdirx9 = pvdirx[j+8];
	                        vdiry9 = pvdiry[j+8];
	                        vdirz9 = pvdirz[j+8];
	                        vrx9   = pvrx[j+8];
	                        vry9   = pvry[j+8];
	                        vrz9   = pvrz[j+8];
	                        omg9   = pomega[j+8];
	                        k9     = pk[j+8];
	                        B_XYZ_ymm4c8(vpolx9,vpoly9,vpolz9,
	                                      vdirx9,vdiry9,vdirz9,
	                                      k9,omg9,vrx9,vry9,vrz9,
	                                      B_x9,B_y9,B_z9);
	                        pB_x[j+8] = B_x9;
	                        pB_y[j+8] = B_y9;
	                        pB_z[j+8] = B_z9; 
	                        vpolx10 = pvpolx[j+9];
	                        vpoly10 = pvpoly[j+9];
	                        vpolz10 = pvpolz[j+9];
	                        vdirx10 = pvdirx[j+9];
	                        vdiry10 = pvdiry[j+9];
	                        vdirz10 = pvdirz[j+9];
	                        vrx10   = pvrx[j+9];
	                        vry10   = pvry[j+9];
	                        vrz10   = pvrz[j+9];
	                        omg10   = pomega[j+9];
	                        k10     = pk[j+9];
	                        B_XYZ_ymm4c8(vpolx10,vpoly10,vpolz10,
	                                      vdirx10,vdiry10,vdirz10,
	                                      k10,omg10,vrx10,vry10,vrz10,
	                                      B_x10,B_y10,B_z10);
	                        pB_x[j+9] = B_x10;
	                        pB_y[j+9] = B_y10;
	                        pB_z[j+9] = B_z10;                                          
	                 }
	                                               
	         }
	         
	      
	         
	         
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_ymm4c8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                 __ATTR_ALIGN__(32) ymm4c8_t k1,k2,k3,k4,k5,k6,k7,k8,k9,k10;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2,B_x3,B_x4,B_x5,B_x6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2,B_y3,B_y4,B_y5,B_y6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2,B_z3,B_z4,B_z5,B_z6;
	                  __m256d vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6;
	                  __m256d vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6;
	                  __m256d vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6;
	                  __m256d vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6;
	                  __m256d vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6;
	                  __m256d vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6;
	                  __m256d vrx1,vrx2,vrx3,vrx4,vrx5,vrx6;
	                  __m256d vry1,vry2,vry3,vry4,vry5,vry6;
	                  __m256d vrz1,vrz2,vrz3,vrz4,vrz5,vrz6;
	                  __m256d omg1,omg2,omg3,omg4,omg5,omg6;
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
	                        B_XYZ_VP_ymm4c8(vpolx,vpoly,vpolz,
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
#pragma omp parallel for schedule(runtime) default(none)                            \
        firstprivate(m1,PF_DIST) private(j,k1,k2,k3,k4,k5,k6)                       \
                                 private(B_x1,B_x2,B_x3,B_x4,B_x5,B_x6)             \
                                 private(B_y1,B_y2,B_y3,B_y4,B_y5,B_y6)             \
                                 private(B_z1,B_z2,B_z3,B_z4,B_z5,B_z6)             \
                                 private(vpolx1,vpolx2,vpolx3,vpolx4,vpolx5,vpolx6) \
                                 private(vpoly1,vpoly2,vpoly3,vpoly4,vpoly5,vpoly6) \
                                 private(vpolz1,vpolz2,vpolz3,vpolz4,vpolz5,vpolz6) \
                                 private(vdirx1,vdirx2,vdirx3,vdirx4,vdirx5,vdirx6) \
                                 private(vdiry1,vdiry2,vdiry3,vdiry4,vdiry5,vdiry6) \
                                 private(vdirz1,vdirz2,vdirz3,vdirz4,vdirz5,vdirz6) \
                                 private(vrx1,vrx2,vrx3,vrx4,vrx5,vrx6)             \
                                 private(vry1,vry2,vry3,vry4,vry5,vry6)             \
                                 private(vrz1,vrz2,vrz3,vrz4,vrz5,vrz6)             \
                                 private(omg1,omg2,omg3,omg4,omg5,omg6)             \
                                 shared(n,pvpolx,pvpoly,pvpolz,pvdirx,pvdiry,pvdirz)\
                                 shared(pvrx,pvry,pvrz,pomega,pk,pB_x,pB_y,pB_z)
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
                                vpolx1 = pvpolx[j+0];
	                        vpoly1 = pvpoly[j+0];
	                        vpolz1 = pvpolz[j+0];
	                        vdirx1 = pvdirx[j+0];
	                        vdiry1 = pvdiry[j+0];
	                        vdirz1 = pvdirz[j+0];
	                        vrx1   = pvrx[j+0];
	                        vry1   = pvry[j+0];
	                        vrz1   = pvrz[j+0];
	                        omg1   = pomega[j+0];
	                        k1     = pk[j+0];
	                        B_XYZ_ymm4c8(vpolx1,vpoly1,vpolz1,
	                                      vdirx1,vdiry1,vdirz1,
	                                      k1,omg1,vrx1,vry1,vrz1,
	                                      B_x1,B_y1,B_z1);
	                        pB_x[j+0] = B_x1;
	                        pB_y[j+0] = B_y1;
	                        pB_z[j+0] = B_z1;    
	                        vpolx2 = pvpolx[j+1];
	                        vpoly2 = pvpoly[j+1];
	                        vpolz2 = pvpolz[j+1];
	                        vdirx2 = pvdirx[j+1];
	                        vdiry2 = pvdiry[j+1];
	                        vdirz2 = pvdirz[j+1];
	                        vrx2   = pvrx[j+1];
	                        vry2   = pvry[j+1];
	                        vrz2   = pvrz[j+1];
	                        omg2   = pomega[j+1];
	                        k2     = pk[j+1];
	                        B_XYZ_ymm4c8(vpolx2,vpoly2,vpolz2,
	                                      vdirx2,vdiry2,vdirz2,
	                                      k2,omg2,vrx2,vry2,vrz2,
	                                      B_x2,B_y2,B_z2);
	                        pB_x[j+1] = B_x2;
	                        pB_y[j+1] = B_y2;
	                        pB_z[j+1] = B_z2;
	                        vpolx3 = pvpolx[j+2];
	                        vpoly3 = pvpoly[j+2];
	                        vpolz3 = pvpolz[j+2];
	                        vdirx3 = pvdirx[j+2];
	                        vdiry3 = pvdiry[j+2];
	                        vdirz3 = pvdirz[j+2];
	                        vrx3   = pvrx[j+2];
	                        vry3   = pvry[j+2];
	                        vrz3   = pvrz[j+2];
	                        omg3   = pomega[j+2];
	                        k3     = pk[j+2];
	                        B_XYZ_ymm4c8(vpolx3,vpoly3,vpolz3,
	                                      vdirx3,vdiry3,vdirz3,
	                                      k3,omg3,vrx3,vry3,vrz3,
	                                      B_x3,B_y3,B_z3);
	                        pB_x[j+2] = B_x3;
	                        pB_y[j+2] = B_y3;
	                        pB_z[j+2] = B_z3;
	                        vpolx4 = pvpolx[j+3];
	                        vpoly4 = pvpoly[j+3];
	                        vpolz4 = pvpolz[j+3];
	                        vdirx4 = pvdirx[j+3];
	                        vdiry4 = pvdiry[j+3];
	                        vdirz4 = pvdirz[j+3];
	                        vrx4   = pvrx[j+3];
	                        vry4   = pvry[j+3];
	                        vrz4   = pvrz[j+3];
	                        omg4   = pomega[j+3];
	                        k4     = pk[j+3];
	                        B_XYZ_ymm4c8(vpolx4,vpoly4,vpolz4,
	                                      vdirx4,vdiry4,vdirz4,
	                                      k4,omg4,vrx4,vry4,vrz4,
	                                      B_x4,B_y4,B_z4);
	                        pB_x[j+3] = B_x4;
	                        pB_y[j+3] = B_y4;
	                        pB_z[j+3] = B_z4;
	                        vpolx5 = pvpolx[j+4];
	                        vpoly5 = pvpoly[j+4];
	                        vpolz5 = pvpolz[j+4];
	                        vdirx5 = pvdirx[j+4];
	                        vdiry5 = pvdiry[j+4];
	                        vdirz5 = pvdirz[j+4];
	                        vrx5   = pvrx[j+4];
	                        vry5   = pvry[j+4];
	                        vrz5   = pvrz[j+4];
	                        omg5   = pomega[j+4];
	                        k5     = pk[j+4];
	                        B_XYZ_ymm4c8(vpolx5,vpoly5,vpolz5,
	                                      vdirx5,vdiry5,vdirz5,
	                                      k5,omg5,vrx5,vry5,vrz5,
	                                      B_x5,B_y5,B_z5);
	                        pB_x[j+4] = B_x5;
	                        pB_y[j+4] = B_y5;
	                        pB_z[j+4] = B_z5;
	                        vpolx6 = pvpolx[j+5];
	                        vpoly6 = pvpoly[j+5];
	                        vpolz6 = pvpolz[j+5];
	                        vdirx6 = pvdirx[j+5];
	                        vdiry6 = pvdiry[j+5];
	                        vdirz6 = pvdirz[j+5];
	                        vrx6   = pvrx[j+5];
	                        vry6   = pvry[j+5];
	                        vrz6   = pvrz[j+5];
	                        omg6   = pomega[j+5];
	                        k6     = pk[j+5];
	                        B_XYZ_ymm4c8(vpolx6,vpoly6,vpolz6,
	                                      vdirx6,vdiry6,vdirz6,
	                                      k6,omg6,vrx6,vry6,vrz6,
	                                      B_x6,B_y6,B_z6);
	                        pB_x[j+5] = B_x6;
	                        pB_y[j+5] = B_y6;
	                        pB_z[j+5] = B_z6;
	                                  
	                 }
	                                               
	         }
	         
	         
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_ymm4c8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256d * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST) {
	                                           
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                 __ATTR_ALIGN__(32) ymm4c8_t k1,k2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2;
	                  __m256d vpolx1,vpolx2;
	                  __m256d vpoly1,vpoly2;
	                  __m256d vpolz1,vpolz2;
	                  __m256d vdirx1,vdirx2;
	                  __m256d vdiry1,vdiry2;
	                  __m256d vdirz1,vdirz2;
	                  __m256d vrx1,vrx2;
	                  __m256d vry1,vry2;
	                  __m256d vrz1,vrz2;
	                  __m256d omg1,omg2;
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
	                        B_XYZ_VP_ymm4c8(vpolx,vpoly,vpolz,
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
#pragma omp parallel for schedule(runtime) default(none)                            \
        firstprivate(m1,PF_DIST) private(j,k1,k2)                       \
                                 private(B_x1,B_x2)             \
                                 private(B_y1,B_y2)             \
                                 private(B_z1,B_z2)             \
                                 private(vpolx1,vpolx2) \
                                 private(vpoly1,vpoly2) \
                                 private(vpolz1,vpolz2) \
                                 private(vdirx1,vdirx2) \
                                 private(vdiry1,vdiry2) \
                                 private(vdirz1,vdirz2) \
                                 private(vrx1,vrx2)             \
                                 private(vry1,vry2)             \
                                 private(vrz1,vrz2)             \
                                 private(omg1,omg2)             \
                                 shared(n,pvpolx,pvpoly,pvpolz,pvdirx,pvdiry,pvdirz)\
                                 shared(pvrx,pvry,pvrz,pomega,pk,pB_x,pB_y,pB_z)
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
                                vpolx1 = pvpolx[j+0];
	                        vpoly1 = pvpoly[j+0];
	                        vpolz1 = pvpolz[j+0];
	                        vdirx1 = pvdirx[j+0];
	                        vdiry1 = pvdiry[j+0];
	                        vdirz1 = pvdirz[j+0];
	                        vrx1   = pvrx[j+0];
	                        vry1   = pvry[j+0];
	                        vrz1   = pvrz[j+0];
	                        omg1   = pomega[j+0];
	                        k1     = pk[j+0];
	                        B_XYZ_ymm4c8(vpolx1,vpoly1,vpolz1,
	                                      vdirx1,vdiry1,vdirz1,
	                                      k1,omg1,vrx1,vry1,vrz1,
	                                      B_x1,B_y1,B_z1);
	                        pB_x[j+0] = B_x1;
	                        pB_y[j+0] = B_y1;
	                        pB_z[j+0] = B_z1;    
	                        vpolx2 = pvpolx[j+1];
	                        vpoly2 = pvpoly[j+1];
	                        vpolz2 = pvpolz[j+1];
	                        vdirx2 = pvdirx[j+1];
	                        vdiry2 = pvdiry[j+1];
	                        vdirz2 = pvdirz[j+1];
	                        vrx2   = pvrx[j+1];
	                        vry2   = pvry[j+1];
	                        vrz2   = pvrz[j+1];
	                        omg2   = pomega[j+1];
	                        k2     = pk[j+1];
	                        B_XYZ_ymm4c8(vpolx2,vpoly2,vpolz2,
	                                      vdirx2,vdiry2,vdirz2,
	                                      k2,omg2,vrx2,vry2,vrz2,
	                                      B_x2,B_y2,B_z2);
	                        pB_x[j+1] = B_x2;
	                        pB_y[j+1] = B_y2;
	                        pB_z[j+1] = B_z2;
	                       	                                  
	                 }
	                                               
	         }
	 
	         
	 	       
	      
	      
	      
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_ymm4c8_unroll10x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;    
	                  __ATTR_ALIGN__(32) ymm4c8_t r1,r2,r3,r4,r5,r6,r7,r8,r9,r10;
	                  __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2,H_x3,H_x4,H_x5,H_x6,H_x7,H_x8,H_x9,H_x10;
	                  __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2,H_y3,H_y4,H_y5,H_y6,H_y7,H_y8,H_y9,H_y10;
	                  __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2,H_z3,H_z4,H_z5,H_z6,H_z7,H_z8,H_z9,H_z10;
	                  __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2,B_x3,B_x4,B_x5,B_x6,B_x7,B_x8,B_x9,B_x10;
	                  __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2,B_y3,B_y4,B_y5,B_y6,B_y7,B_y8,B_y9,B_y10;
	                  __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2,B_z3,B_z4,B_z5,B_z6,B_z7,B_z8,B_z9,B_z10;
	                  __m256d tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10;
	                  __m256d phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10;
	                  __m256d psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8,spi9,psi10;
	                  __m256d omg1,omg2,omg3,omg4,omg5,omg6,omg7,omg8,omg9,omg10;
	                  __m256d px1,px2,px3,px4,px5,px6,px7,px8,px9,px10;
	                  __m256d py1,py2,py3,py4,py5,py6,py7,py8,py9,py10;
	                  __m256d pz1,pz2,pz3,pz4,pz5,pz6,pz7,pz8,pz9,pz10;
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
	                         B_XYZ_H_XYZ_P_ymm4c8(tht,phi,psi,
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
#pragma omp parallel for schedule(runtime) defualt(none)                                          \
        firstprivate(m1,PF_DIST) private(j,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)                        \
                                 private(H_x1,H_x2,H_x3,H_x4,H_x5,H_x6,H_x7,H_x8,H_x9,H_x10)      \
                                 private(H_y1,H_y2,H_y3,H_y4,H_y5,H_y6,H_y7,H_y8,H_y9,H_y10)      \
                                 private(H_z1,H_z2,H_z3,H_z4,H_z5,H_z6,H_z7,H_z8,H_z9,H_z10)      \
                                 private(B_x1,B_x2,B_x3,B_x4,B_x5,B_x6,B_x7,B_x8,B_x9,B_x10)      \
                                 private(B_y1,B_y2,B_y3,B_y4,B_y5,B_y6,B_y7,B_y8,B_y9,B_y10)      \
                                 private(B_z1,B_z2,B_z3,B_z4,B_z5,B_z6,B_z7,B_z8,B_z9,B_z10)      \
                                 private(tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10)      \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10)      \
                                 private(psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8,spi9,psi10)      \
                                 private(omg1,omg2,omg3,omg4,omg5,omg6,omg7,omg8,omg9,omg10)      \
                                 private(px1,px2,px3,px4,px5,px6,px7,px8,px9,px10)                \
                                 private(py1,py2,py3,py4,py5,py6,py7,py8,py9,py10)                \
                                 private(pz1,pz2,pz3,pz4,pz5,pz6,pz7,pz8,pz9,pz10)                \
                                 shared(n,ptht,pphi,ppsi,pomg,ppx,ppy,ppz,pr)                     \
                                 shared(pH_x,pH_y,pH_z,pB_x,pB_y,pB_z)
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
                                 tht1 = ptht[j+0];
	                         phi1 = pphi[j+0];
	                         psi1 = ppsi[j+0];
	                         omg1 = pomg[j+0];
	                         px1  = ppx[j+0];
	                         py1  = ppy[j+0];
	                         pz1  = ppz[j+0];
	                         r1   = pr[j+0];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht1,phi1,psi1,
	                                               omg1,px1,py1,pz1,r1,
	                                               H_x1,H_y1,H_z1,
	                                               B_x1,B_y1,B_z1);
	                         pH_x[j+0] = H_x1;
	                         pH_y[j+0] = H_y1;
	                         pH_z[j+0] = H_z1;
	                         pB_x[j+0] = B_x1;
	                         pB_y[j+0] = B_y1;
	                         pB_z[j+0] = B_z1;  
	                         tht2 = ptht[j+1];
	                         phi2 = pphi[j+1];
	                         psi2 = ppsi[j+1];
	                         omg2 = pomg[j+1];
	                         px2  = ppx[j+1];
	                         py2  = ppy[j+1];
	                         pz2  = ppz[j+1];
	                         r2   = pr[j+1];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht2,phi2,psi2,
	                                               omg2,px2,py2,pz2,r2,
	                                               H_x2,H_y2,H_z2,
	                                               B_x2,B_y2,B_z2);
	                         pH_x[j+1] = H_x2;
	                         pH_y[j+1] = H_y2;
	                         pH_z[j+1] = H_z2;
	                         pB_x[j+1] = B_x2;
	                         pB_y[j+1] = B_y2;
	                         pB_z[j+1] = B_z2;  
	                         tht3 = ptht[j+2];
	                         phi3 = pphi[j+2];
	                         psi3 = ppsi[j+2];
	                         omg3 = pomg[j+2];
	                         px3  = ppx[j+2];
	                         py3  = ppy[j+2];
	                         pz3  = ppz[j+2];
	                         r3   = pr[j+2];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht3,phi3,psi3,
	                                               omg3,px3,py3,pz3,r3,
	                                               H_x3,H_y3,H_z3,
	                                               B_x3,B_y3,B_z3);
	                         pH_x[j+2] = H_x3;
	                         pH_y[j+2] = H_y3;
	                         pH_z[j+2] = H_z3;
	                         pB_x[j+2] = B_x3;
	                         pB_y[j+2] = B_y3;
	                         pB_z[j+2] = B_z3;  
	                         tht4 = ptht[j+3];
	                         phi4 = pphi[j+3];
	                         psi4 = ppsi[j+3];
	                         omg4 = pomg[j+3];
	                         px4  = ppx[j+3];
	                         py4  = ppy[j+3];
	                         pz4  = ppz[j+3];
	                         r4   = pr[j+3];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht4,phi4,psi4,
	                                               omg4,px4,py4,pz4,r4,
	                                               H_x4,H_y4,H_z4,
	                                               B_x4,B_y4,B_z4);
	                         pH_x[j+3] = H_x4;
	                         pH_y[j+3] = H_y4;
	                         pH_z[j+3] = H_z4;
	                         pB_x[j+3] = B_x4;
	                         pB_y[j+3] = B_y4;
	                         pB_z[j+3] = B_z4;  
	                         tht5 = ptht[j+4];
	                         phi5 = pphi[j+4];
	                         psi5 = ppsi[j+4];
	                         omg5 = pomg[j+4];
	                         px5  = ppx[j+4];
	                         py5  = ppy[j+4];
	                         pz5  = ppz[j+4];
	                         r5   = pr[j+4];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht5,phi5,psi5,
	                                               omg5,px5,py5,pz5,r5,
	                                               H_x5,H_y5,H_z5,
	                                               B_x5,B_y5,B_z5);
	                         pH_x[j+4] = H_x5;
	                         pH_y[j+4] = H_y5;
	                         pH_z[j+4] = H_z5;
	                         pB_x[j+4] = B_x5;
	                         pB_y[j+4] = B_y5;
	                         pB_z[j+4] = B_z5;  
	                         tht6 = ptht[j+5];
	                         phi6 = pphi[j+5];
	                         psi6 = ppsi[j+5];
	                         omg6 = pomg[j+5];
	                         px6  = ppx[j+5];
	                         py6  = ppy[j+5];
	                         pz6  = ppz[j+5];
	                         r6   = pr[j+5];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht6,phi6,psi6,
	                                               omg6,px6,py6,pz6,r6,
	                                               H_x6,H_y6,H_z6,
	                                               B_x6,B_y6,B_z6);
	                         pH_x[j+5] = H_x6;
	                         pH_y[j+5] = H_y6;
	                         pH_z[j+5] = H_z6;
	                         pB_x[j+5] = B_x6;
	                         pB_y[j+5] = B_y6;
	                         pB_z[j+5] = B_z6;  
	                         tht7 = ptht[j+6];
	                         phi7 = pphi[j+6];
	                         psi7 = ppsi[j+6];
	                         omg7 = pomg[j+6];
	                         px7  = ppx[j+6];
	                         py7  = ppy[j+6];
	                         pz7  = ppz[j+6];
	                         r7   = pr[j+6];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht7,phi7,psi7,
	                                               omg7,px7,py7,pz7,r7,
	                                               H_x7,H_y7,H_z7,
	                                               B_x7,B_y7,B_z7);
	                         pH_x[j+6] = H_x7;
	                         pH_y[j+6] = H_y7;
	                         pH_z[j+6] = H_z7;
	                         pB_x[j+6] = B_x7;
	                         pB_y[j+6] = B_y7;
	                         pB_z[j+6] = B_z7;  
	                         tht8 = ptht[j+7];
	                         phi8 = pphi[j+7];
	                         psi8 = ppsi[j+7];
	                         omg8 = pomg[j+7];
	                         px8  = ppx[j+7];
	                         py8  = ppy[j+7];
	                         pz8  = ppz[j+7];
	                         r8   = pr[j+7];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht8,phi8,psi8,
	                                               omg8,px8,py8,pz8,r8,
	                                               H_x8,H_y8,H_z8,
	                                               B_x8,B_y8,B_z8);
	                         pH_x[j+7] = H_x8;
	                         pH_y[j+7] = H_y8;
	                         pH_z[j+7] = H_z8;
	                         pB_x[j+7] = B_x8;
	                         pB_y[j+7] = B_y8;
	                         pB_z[j+7] = B_z8;  
	                         tht9 = ptht[j+8];
	                         phi9 = pphi[j+8];
	                         psi9 = ppsi[j+8];
	                         omg9 = pomg[j+8];
	                         px9  = ppx[j+8];
	                         py9  = ppy[j+8];
	                         pz9  = ppz[j+8];
	                         r9   = pr[j+8];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht9,phi9,psi9,
	                                               omg9,px9,py9,pz9,r9,
	                                               H_x9,H_y9,H_z9,
	                                               B_x9,B_y9,B_z9);
	                         pH_x[j+8] = H_x9;
	                         pH_y[j+8] = H_y9;
	                         pH_z[j+8] = H_z9;
	                         pB_x[j+8] = B_x9;
	                         pB_y[j+8] = B_y9;
	                         pB_z[j+8] = B_z9;  
	                         tht10 = ptht[j+9];
	                         phi10 = pphi[j+9];
	                         psi10 = ppsi[j+9];
	                         omg10 = pomg[j+9];
	                         px10  = ppx[j+9];
	                         py10  = ppy[j+9];
	                         pz10  = ppz[j+9];
	                         r10   = pr[j+9];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht10,phi10,psi10,
	                                               omg10,px10,py10,pz10,r10,
	                                               H_x10,H_y10,H_z10,
	                                               B_x10,B_y10,B_z10);
	                         pH_x[j+9] = H_x10;
	                         pH_y[j+9] = H_y10;
	                         pH_z[j+9] = H_z10;
	                         pB_x[j+9] = B_x10;
	                         pB_y[j+9] = B_y10;
	                         pB_z[j+9] = B_z10;  
	                                                                               
	                   }
	                                             
	      }
	      
	      
	      
	      

	      
	      
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_ymm4c8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;    
	                 __ATTR_ALIGN__(32) ymm4c8_t r1,r2,r3,r4,r5,r6;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2,H_x3,H_x4,H_x5,H_x6;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2,H_y3,H_y4,H_y5,H_y6;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2,H_z3,H_z4,H_z5,H_z6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2,B_x3,B_x4,B_x5,B_x6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2,B_y3,B_y4,B_y5,B_y6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2,B_z3,B_z4,B_z5,B_z6;
	                  __m256d tht1,tht2,tht3,tht4,tht5,tht6;
	                  __m256d phi1,phi2,phi3,phi4,phi5,phi6;
	                  __m256d psi1,psi2,psi3,psi4,psi5,psi6;
	                  __m256d omg1,omg2,omg3,omg4,omg5,omg6;
	                  __m256d px1,px2,px3,px4,px5,px6;
	                  __m256d py1,py2,py3,py4,py5,py6;
	                  __m256d pz1,pz2,pz3,pz4,pz5,pz6;
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
	                         B_XYZ_H_XYZ_P_ymm4c8(tht,phi,psi,
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
#pragma omp parallel for schedule(runtime) defualt(none)                     \
        firstprivate(m1,PF_DIST) private(j,r1,r2,r3,r4,r5,r6)                \
                                 private(H_x1,H_x2,H_x3,H_x4,H_x5,H_x6)      \
                                 private(H_y1,H_y2,H_y3,H_y4,H_y5,H_y6)      \
                                 private(H_z1,H_z2,H_z3,H_z4,H_z5,H_z6)      \
                                 private(B_x1,B_x2,B_x3,B_x4,B_x5,B_x6)      \
                                 private(B_y1,B_y2,B_y3,B_y4,B_y5,B_y6)      \
                                 private(B_z1,B_z2,B_z3,B_z4,B_z5,B_z6)      \
                                 private(tht1,tht2,tht3,tht4,tht5,tht6)      \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6)      \
                                 private(psi1,psi2,psi3,psi4,psi5,psi6)      \
                                 private(omg1,omg2,omg3,omg4,omg5,omg6)      \
                                 private(px1,px2,px3,px4,px5,px6)            \
                                 private(py1,py2,py3,py4,py5,py6)            \
                                 private(pz1,pz2,pz3,pz4,pz5,pz6)            \
                                 shared(n,ptht,pphi,ppsi,pomg,ppx,ppy,ppz,pr)\
                                 shared(pH_x,pH_y,pH_z,pB_x,pB_y,pB_z)
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
                                 tht1 = ptht[j+0];
	                         phi1 = pphi[j+0];
	                         psi1 = ppsi[j+0];
	                         omg1 = pomg[j+0];
	                         px1  = ppx[j+0];
	                         py1  = ppy[j+0];
	                         pz1  = ppz[j+0];
	                         r1   = pr[j+0];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht1,phi1,psi1,
	                                               omg1,px1,py1,pz1,r1,
	                                               H_x1,H_y1,H_z1,
	                                               B_x1,B_y1,B_z1);
	                         pH_x[j+0] = H_x1;
	                         pH_y[j+0] = H_y1;
	                         pH_z[j+0] = H_z1;
	                         pB_x[j+0] = B_x1;
	                         pB_y[j+0] = B_y1;
	                         pB_z[j+0] = B_z1;  
	                         tht2 = ptht[j+1];
	                         phi2 = pphi[j+1];
	                         psi2 = ppsi[j+1];
	                         omg2 = pomg[j+1];
	                         px2  = ppx[j+1];
	                         py2  = ppy[j+1];
	                         pz2  = ppz[j+1];
	                         r2   = pr[j+1];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht2,phi2,psi2,
	                                               omg2,px2,py2,pz2,r2,
	                                               H_x2,H_y2,H_z2,
	                                               B_x2,B_y2,B_z2);
	                         pH_x[j+1] = H_x2;
	                         pH_y[j+1] = H_y2;
	                         pH_z[j+1] = H_z2;
	                         pB_x[j+1] = B_x2;
	                         pB_y[j+1] = B_y2;
	                         pB_z[j+1] = B_z2;  
	                         tht3 = ptht[j+2];
	                         phi3 = pphi[j+2];
	                         psi3 = ppsi[j+2];
	                         omg3 = pomg[j+2];
	                         px3  = ppx[j+2];
	                         py3  = ppy[j+2];
	                         pz3  = ppz[j+2];
	                         r3   = pr[j+2];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht3,phi3,psi3,
	                                               omg3,px3,py3,pz3,r3,
	                                               H_x3,H_y3,H_z3,
	                                               B_x3,B_y3,B_z3);
	                         pH_x[j+2] = H_x3;
	                         pH_y[j+2] = H_y3;
	                         pH_z[j+2] = H_z3;
	                         pB_x[j+2] = B_x3;
	                         pB_y[j+2] = B_y3;
	                         pB_z[j+2] = B_z3;  
	                         tht4 = ptht[j+3];
	                         phi4 = pphi[j+3];
	                         psi4 = ppsi[j+3];
	                         omg4 = pomg[j+3];
	                         px4  = ppx[j+3];
	                         py4  = ppy[j+3];
	                         pz4  = ppz[j+3];
	                         r4   = pr[j+3];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht4,phi4,psi4,
	                                               omg4,px4,py4,pz4,r4,
	                                               H_x4,H_y4,H_z4,
	                                               B_x4,B_y4,B_z4);
	                         pH_x[j+3] = H_x4;
	                         pH_y[j+3] = H_y4;
	                         pH_z[j+3] = H_z4;
	                         pB_x[j+3] = B_x4;
	                         pB_y[j+3] = B_y4;
	                         pB_z[j+3] = B_z4;  
	                         tht5 = ptht[j+4];
	                         phi5 = pphi[j+4];
	                         psi5 = ppsi[j+4];
	                         omg5 = pomg[j+4];
	                         px5  = ppx[j+4];
	                         py5  = ppy[j+4];
	                         pz5  = ppz[j+4];
	                         r5   = pr[j+4];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht5,phi5,psi5,
	                                               omg5,px5,py5,pz5,r5,
	                                               H_x5,H_y5,H_z5,
	                                               B_x5,B_y5,B_z5);
	                         pH_x[j+4] = H_x5;
	                         pH_y[j+4] = H_y5;
	                         pH_z[j+4] = H_z5;
	                         pB_x[j+4] = B_x5;
	                         pB_y[j+4] = B_y5;
	                         pB_z[j+4] = B_z5;  
	                         tht6 = ptht[j+5];
	                         phi6 = pphi[j+5];
	                         psi6 = ppsi[j+5];
	                         omg6 = pomg[j+5];
	                         px6  = ppx[j+5];
	                         py6  = ppy[j+5];
	                         pz6  = ppz[j+5];
	                         r6   = pr[j+5];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht6,phi6,psi6,
	                                               omg6,px6,py6,pz6,r6,
	                                               H_x6,H_y6,H_z6,
	                                               B_x6,B_y6,B_z6);
	                         pH_x[j+5] = H_x6;
	                         pH_y[j+5] = H_y6;
	                         pH_z[j+5] = H_z6;
	                         pB_x[j+5] = B_x6;
	                         pB_y[j+5] = B_y6;
	                         pB_z[j+5] = B_z6;  
	                        	                                                                               
	                   }
	                                             
	      }
	      
	      
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_ymm4c8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256d * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm4c8_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST) {
	                                                
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;    
	                 __ATTR_ALIGN__(32) ymm4c8_t r1,r2;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2;
	                  __m256d tht1,tht2;
	                  __m256d phi1,phi2;
	                  __m256d psi1,psi2;
	                  __m256d omg1,omg2;
	                  __m256d px1,px2;
	                  __m256d py1,py2;
	                  __m256d pz1,pz2;
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
	                         B_XYZ_H_XYZ_P_ymm4c8(tht,phi,psi,
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
#pragma omp parallel for schedule(runtime) defualt(none) \
        firstprivate(m1,PF_DIST) private(j,r1,r2)        \
                                 private(H_x1,H_x2)      \
                                 private(H_y1,H_y2)      \
                                 private(H_z1,H_z2)      \
                                 private(B_x1,B_x2)      \
                                 private(B_y1,B_y2)      \
                                 private(B_z1,B_z2)      \
                                 private(tht1,tht2)      \
                                 private(phi1,phi2)      \
                                 private(psi1,psi2)      \
                                 private(omg1,omg2)      \
                                 private(px1,px2)        \
                                 private(py1,py2)        \
                                 private(pz1,pz2)        \
                                 shared(n,ptht,pphi,ppsi,pomg,ppx,ppy,ppz,pr)\
                                 shared(pH_x,pH_y,pH_z,pB_x,pB_y,pB_z)
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
                                 tht1 = ptht[j+0];
	                         phi1 = pphi[j+0];
	                         psi1 = ppsi[j+0];
	                         omg1 = pomg[j+0];
	                         px1  = ppx[j+0];
	                         py1  = ppy[j+0];
	                         pz1  = ppz[j+0];
	                         r1   = pr[j+0];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht1,phi1,psi1,
	                                               omg1,px1,py1,pz1,r1,
	                                               H_x1,H_y1,H_z1,
	                                               B_x1,B_y1,B_z1);
	                         pH_x[j+0] = H_x1;
	                         pH_y[j+0] = H_y1;
	                         pH_z[j+0] = H_z1;
	                         pB_x[j+0] = B_x1;
	                         pB_y[j+0] = B_y1;
	                         pB_z[j+0] = B_z1;  
	                         tht2 = ptht[j+1];
	                         phi2 = pphi[j+1];
	                         psi2 = ppsi[j+1];
	                         omg2 = pomg[j+1];
	                         px2  = ppx[j+1];
	                         py2  = ppy[j+1];
	                         pz2  = ppz[j+1];
	                         r2   = pr[j+1];
	                         B_XYZ_H_XYZ_P_ymm4c8(tht2,phi2,psi2,
	                                               omg2,px2,py2,pz2,r2,
	                                               H_x2,H_y2,H_z2,
	                                               B_x2,B_y2,B_z2);
	                         pH_x[j+1] = H_x2;
	                         pH_y[j+1] = H_y2;
	                         pH_z[j+1] = H_z2;
	                         pB_x[j+1] = B_x2;
	                         pB_y[j+1] = B_y2;
	                         pB_z[j+1] = B_z2;  
	                        
	                        	                                                                               
	                   }
	                                             
	      }
	  
	           
	       
	       
	       
	       /*
	              ! Electric and Magnetic Fields elliptically polarized
	       */
	       
	       
	       
	
	      
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_ymm4c8_unroll10x_omp( const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                                      const __m256d * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 10;
	                  __ATTR_ALIGN__(32) ymm4c8_t phase1,phase2,phase3,phase4,phase5,phase6,phase7,phase8,phase9,phase10;
	                  __ATTR_ALIGN__(32) ymm4c8_t refi1,refi2,refi3,refi4,refi5,refi6,refi7,refi8,refi9,refi10;
	                  __ATTR_ALIGN__(32) ymm4c8_t px1,px2,px3,px4,px5,px6,px7,px8,px9,px10;
	                  __ATTR_ALIGN__(32) ymm4c8_t py1,py2,py3,py4,py5,py6,py7,py8,py9,py10;
	                  __ATTR_ALIGN__(32) ymm4c8_t pz1,pz2,pz3,pz4,pz5,pz6,pz7,pz8,pz9,pz10;
	                  __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2,H_x3,H_x4,H_x5,H_x6,H_x7,h_x8,H_x9,H_x10;
	                  __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2,H_y3,H_y4,H_y5,H_y6,H_y7,h_y8,H_y9,H_y10;
	                  __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2,H_z3,H_z4,H_z5,H_z6,H_z7,h_z8,H_z9,H_z10;
	                  __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2,B_x3,B_x4,B_x5,B_x6,B_x7,B_x8,B_x9,B_x10;
	                  __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2,B_y3,B_y4,B_y5,B_y6,B_y7,B_y8,B_y9,B_y10;
	                  __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2,B_z3,B_z4,B_z5,B_z6,B_z7,B_z8,B_z9,B_z10;
	                 register __m256d tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10;
	                 register __m256d phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10;
	                 register __m256d omg1,omg2,omg3,omg4,omg5,omg6,omg7,omg8,omg9,omg10;
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
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht,phi,omg,
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
#pragma omp parallel for schedule(runtime) default(none)                                                            \
        firstprivate(m1,PF_DIST) private(j,phase1,phase2,phase3,phase4,phase5,phase6,phase7,phase8,phase9,phase10)  \
                                 private(refi1,refi2,refi3,refi4,refi5,refi6,refi7,refi8,refi9,refi10)              \
                                 private(px1,px2,px3,px4,px5,px6,px7,px8,px9,px10)                                  \
                                 private(py1,py2,py3,py4,py5,py6,py7,py8,py9,py10)                                  \
                                 private(pz1,pz2,pz3,pz4,pz5,pz6,pz7,pz8,pz9,pz10)                                  \
                                 private(H_x1,H_x2,H_x3,H_x4,H_x5,H_x6,H_x7,h_x8,H_x9,H_x10)                        \
                                 private(H_y1,H_y2,H_y3,H_y4,H_y5,H_y6,H_y7,h_y8,H_y9,H_y10)                        \
                                 private(H_z1,H_z2,H_z3,H_z4,H_z5,H_z6,H_z7,h_z8,H_z9,H_z10)                        \
                                 private(B_x1,B_x2,B_x3,B_x4,B_x5,B_x6,B_x7,B_x8,B_x9,B_x10)                        \
                                 private(B_y1,B_y2,B_y3,B_y4,B_y5,B_y6,B_y7,B_y8,B_y9,B_y10)                        \
                                 private(B_z1,B_z2,B_z3,B_z4,B_z5,B_z6,B_z7,B_z8,B_z9,B_z10)                        \
                                 private(tht1,tht2,tht3,tht4,tht5,tht6,tht7,tht8,tht9,tht10)                        \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10)                        \
                                 private(omg1,omg2,omg3,omg4,omg5,omg6,omg7,omg8,omg9,omg10)                        \
                                 shared(n,ptht,pphi,pomg,pphase,prefi)                                              \
                                 shared(ppx,ppy,ppz)                                                                \
                                 shared(pH_x,pH_y,pH_z,pB_x,pB_y,pB_z)
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
                                tht1   = ptht[j+0];
	                        phi1   = pphi[j+0];
	                        omg1   = pomg[j+0];
	                        phase1 = pphase[j+0];
	                        refi1  = prefi[j+0];
	                        px1    = ppx[j+0];
	                        py1    = ppy[j+0];
	                        pz1    = ppz[j+0];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht1,phi1,omg1,
	                                               phase1,refi1,px1,
	                                               py1,pz1,
	                                               H_x1,H_y1,H_z1,
	                                               B_x1,B_y1,B_z1);
	                        pH_x[j+0] = H_x1;
	                        pH_y[j+0] = H_y1;
	                        pH_z[j+0] = H_z1;
	                        pB_x[j+0] = B_x1;
	                        pB_y[j+0] = B_y1;
	                        pB_z[j+0] = B_z1;
	                        tht2   = ptht[j+1];
	                        phi2   = pphi[j+1];
	                        omg2   = pomg[j+1];
	                        phase2 = pphase[j+1];
	                        refi2  = prefi[j+1];
	                        px2    = ppx[j+1];
	                        py2    = ppy[j+1];
	                        pz2    = ppz[j+1];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht2,phi2,omg2,
	                                               phase2,refi2,px2,
	                                               py2,pz2,
	                                               H_x2,H_y2,H_z2,
	                                               B_x2,B_y2,B_z2);
	                        pH_x[j+1] = H_x2;
	                        pH_y[j+1] = H_y2;
	                        pH_z[j+1] = H_z2;
	                        pB_x[j+1] = B_x2;
	                        pB_y[j+1] = B_y2;
	                        pB_z[j+1] = B_z2;  
	                        tht3   = ptht[j+2];
	                        phi3   = pphi[j+2];
	                        omg3   = pomg[j+2];
	                        phase3 = pphase[j+2];
	                        refi3  = prefi[j+2];
	                        px3    = ppx[j+2];
	                        py3    = ppy[j+2];
	                        pz3    = ppz[j+2];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht3,phi3,omg3,
	                                               phase3,refi3,px3,
	                                               py3,pz3,
	                                               H_x3,H_y3,H_z3,
	                                               B_x3,B_y3,B_z3);
	                        pH_x[j+2] = H_x3;
	                        pH_y[j+2] = H_y3;
	                        pH_z[j+2] = H_z3;
	                        pB_x[j+2] = B_x3;
	                        pB_y[j+2] = B_y3;
	                        pB_z[j+2] = B_z3;
	                        tht4   = ptht[j+3];
	                        phi4   = pphi[j+3];
	                        omg4   = pomg[j+3];
	                        phase4 = pphase[j+3];
	                        refi4  = prefi[j+3];
	                        px4    = ppx[j+3];
	                        py4    = ppy[j+3];
	                        pz4    = ppz[j+3];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht4,phi4,omg4,
	                                               phase4,refi4,px4,
	                                               py4,pz4,
	                                               H_x4,H_y4,H_z4,
	                                               B_x4,B_y4,B_z4);
	                        pH_x[j+3] = H_x4;
	                        pH_y[j+3] = H_y4;
	                        pH_z[j+3] = H_z4;
	                        pB_x[j+3] = B_x4;
	                        pB_y[j+3] = B_y4;
	                        pB_z[j+3] = B_z4;  
	                        tht5   = ptht[j+4];
	                        phi5   = pphi[j+4];
	                        omg5   = pomg[j+4];
	                        phase5 = pphase[j+4];
	                        refi5  = prefi[j+4];
	                        px5    = ppx[j+4];
	                        py5    = ppy[j+4];
	                        pz5    = ppz[j+4];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht5,phi5,omg5,
	                                               phase5,refi5,px5,
	                                               py5,pz5,
	                                               H_x5,H_y5,H_z5,
	                                               B_x5,B_y5,B_z);
	                        pH_x[j+4] = H_x5;
	                        pH_y[j+4] = H_y;
	                        pH_z[j+4] = H_z5;
	                        pB_x[j+4] = B_x5;
	                        pB_y[j+4] = B_y5;
	                        pB_z[j+4] = B_z5; 
	                        tht6   = ptht[j+5];
	                        phi6   = pphi[j+5];
	                        omg6   = pomg[j+5];
	                        phase6 = pphase[j+5];
	                        refi6  = prefi[j+5];
	                        px6    = ppx[j+5];
	                        py6    = ppy[j+5];
	                        pz6    = ppz[j+5];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht6,phi6,omg6,
	                                               phase6,refi6,px6,
	                                               py6,pz6,
	                                               H_x6,H_y6,H_z6,
	                                               B_x6,B_y6,B_z6);
	                        pH_x[j+5] = H_x6;
	                        pH_y[j+5] = H_y6;
	                        pH_z[j+5] = H_z6;
	                        pB_x[j+5] = B_x6;
	                        pB_y[j+5] = B_y6;
	                        pB_z[j+5] = B_z6;  
	                        tht7   = ptht[j+6];
	                        phi7   = pphi[j+6];
	                        omg7   = pomg[j+6];
	                        phase7 = pphase[j+6];
	                        refi7  = prefi[j+6];
	                        px7    = ppx[j+6];
	                        py7    = ppy[j+6];
	                        pz7    = ppz[j+6];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht7,phi7,omg7,
	                                               phase7,refi7,px7,
	                                               py7,pz7,
	                                               H_x7,H_y7,H_z7,
	                                               B_x7,B_y7,B_z7);
	                        pH_x[j+6] = H_x7;
	                        pH_y[j+6] = H_y7;
	                        pH_z[j+6] = H_z7;
	                        pB_x[j+6] = B_x7;
	                        pB_y[j+6] = B_y7;
	                        pB_z[j+6] = B_z7;  
	                        tht8   = ptht[j+7];
	                        phi8   = pphi[j+7];
	                        omg8   = pomg[j+7];
	                        phase8 = pphase[j+7];
	                        refi8  = prefi[j+7];
	                        px8    = ppx[j+7];
	                        py8    = ppy[j+7];
	                        pz8    = ppz[j+7];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht8,phi8,omg8,
	                                               phase8,refi8,px8,
	                                               py8,pz8,
	                                               H_x8,H_y8,H_z8,
	                                               B_x8,B_y8,B_z8);
	                        pH_x[j+7] = H_x8;
	                        pH_y[j+7] = H_y8;
	                        pH_z[j+7] = H_z8;
	                        pB_x[j+7] = B_x8;
	                        pB_y[j+7] = B_y8;
	                        pB_z[j+7] = B_z8;  
	                        tht9   = ptht[j+8];
	                        phi9   = pphi[j+8];
	                        omg9   = pomg[j+8];
	                        phase9 = pphase[j+8];
	                        refi9  = prefi[j+8];
	                        px9    = ppx[j+8];
	                        py9    = ppy[j+8];
	                        pz9    = ppz[j+8];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht9,phi9,omg9,
	                                               phase9,refi9,px9,
	                                               py9,pz9,
	                                               H_x9,H_y9,H_z9,
	                                               B_x9,B_y9,B_z9);
	                        pH_x[j+8] = H_x9;
	                        pH_y[j+8] = H_y9;
	                        pH_z[j+8] = H_z9;
	                        pB_x[j+8] = B_x9;
	                        pB_y[j+8] = B_y9;
	                        pB_z[j+8] = B_z9; 
	                        tht10   = ptht[j+9];
	                        phi10   = pphi[j+9];
	                        omg10   = pomg[j+9];
	                        phase10 = pphase[j+9];
	                        refi10  = prefi[j+9];
	                        px10    = ppx[j+9];
	                        py10    = ppy[j+9];
	                        pz10    = ppz[j+9];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht10,phi10,omg10,
	                                               phase10,refi10,px10,
	                                               py10,pz10,
	                                               H_x10,H_y10,H_z10,
	                                               B_x10,B_y10,B_z10);
	                        pH_x[j+9] = H_x10;
	                        pH_y[j+9] = H_y10;
	                        pH_z[j+9] = H_z10;
	                        pB_x[j+9] = B_x10;
	                        pB_y[j+9] = B_y10;
	                        pB_z[j+9] = B_z10;
	                 }               
	      }
	        
	      
	      
	        
	           
	      
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_ymm4c8_unroll6x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                                      const __m256d * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 6;
	                 __ATTR_ALIGN__(32) ymm4c8_t phase1,phase2,phase3,phase4,phase5,phase6;
	                 __ATTR_ALIGN__(32) ymm4c8_t refi1,refi2,refi3,refi4,refi5,refi6;
	                 __ATTR_ALIGN__(32) ymm4c8_t px1,px2,px3,px4,px5,px6;
	                 __ATTR_ALIGN__(32) ymm4c8_t py1,py2,py3,py4,py5,py6;
	                 __ATTR_ALIGN__(32) ymm4c8_t pz1,pz2,pz3,pz4,pz5,pz6;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2,H_x3,H_x4,H_x5,H_x6;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2,H_y3,H_y4,H_y5,H_y6;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2,H_z3,H_z4,H_z5,H_z6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2,B_x3,B_x4,B_x5,B_x6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2,B_y3,B_y4,B_y5,B_y6;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2,B_z3,B_z4,B_z5,B_z6;
	                 register __m256d tht1,tht2,tht3,tht4,tht5,tht6;
	                 register __m256d phi1,phi2,phi3,phi4,phi5,phi6;
	                 register __m256d omg1,omg2,omg3,omg4,omg5,omg6;
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
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht,phi,omg,
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
#pragma omp parallel for schedule(runtime) default(none)                               \
        firstprivate(m1,PF_DIST) private(j,phase1,phase2,phase3,phase4,phase5,phase6)  \
                                 private(refi1,refi2,refi3,refi4,refi5,refi6)          \
                                 private(px1,px2,px3,px4,px5,px6)                      \
                                 private(py1,py2,py3,py4,py5,py6)                      \
                                 private(pz1,pz2,pz3,pz4,pz5,pz6)                      \
                                 private(H_x1,H_x2,H_x3,H_x4,H_x5,H_x6)                \
                                 private(H_y1,H_y2,H_y3,H_y4,H_y5,H_y6)                \
                                 private(H_z1,H_z2,H_z3,H_z4,H_z5,H_z6)                \
                                 private(B_x1,B_x2,B_x3,B_x4,B_x5,B_x6)                \
                                 private(B_y1,B_y2,B_y3,B_y4,B_y5,B_y6)                \
                                 private(B_z1,B_z2,B_z3,B_z4,B_z5,B_z6)                \
                                 private(tht1,tht2,tht3,tht4,tht5,tht6)                \
                                 private(phi1,phi2,phi3,phi4,phi5,phi6)                \
                                 private(omg1,omg2,omg3,omg4,omg5,omg6)                \
                                 shared(n,ptht,pphi,pomg,pphase,prefi)                 \
                                 shared(ppx,ppy,ppz)                                   \
                                 shared(pH_x,pH_y,pH_z,pB_x,pB_y,pB_z)
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
                                tht1   = ptht[j+0];
	                        phi1   = pphi[j+0];
	                        omg1   = pomg[j+0];
	                        phase1 = pphase[j+0];
	                        refi1  = prefi[j+0];
	                        px1    = ppx[j+0];
	                        py1    = ppy[j+0];
	                        pz1    = ppz[j+0];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht1,phi1,omg1,
	                                               phase1,refi1,px1,
	                                               py1,pz1,
	                                               H_x1,H_y1,H_z1,
	                                               B_x1,B_y1,B_z1);
	                        pH_x[j+0] = H_x1;
	                        pH_y[j+0] = H_y1;
	                        pH_z[j+0] = H_z1;
	                        pB_x[j+0] = B_x1;
	                        pB_y[j+0] = B_y1;
	                        pB_z[j+0] = B_z1;
	                        tht2   = ptht[j+1];
	                        phi2   = pphi[j+1];
	                        omg2   = pomg[j+1];
	                        phase2 = pphase[j+1];
	                        refi2  = prefi[j+1];
	                        px2    = ppx[j+1];
	                        py2    = ppy[j+1];
	                        pz2    = ppz[j+1];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht2,phi2,omg2,
	                                               phase2,refi2,px2,
	                                               py2,pz2,
	                                               H_x2,H_y2,H_z2,
	                                               B_x2,B_y2,B_z2);
	                        pH_x[j+1] = H_x2;
	                        pH_y[j+1] = H_y2;
	                        pH_z[j+1] = H_z2;
	                        pB_x[j+1] = B_x2;
	                        pB_y[j+1] = B_y2;
	                        pB_z[j+1] = B_z2;  
	                        tht3   = ptht[j+2];
	                        phi3   = pphi[j+2];
	                        omg3   = pomg[j+2];
	                        phase3 = pphase[j+2];
	                        refi3  = prefi[j+2];
	                        px3    = ppx[j+2];
	                        py3    = ppy[j+2];
	                        pz3    = ppz[j+2];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht3,phi3,omg3,
	                                               phase3,refi3,px3,
	                                               py3,pz3,
	                                               H_x3,H_y3,H_z3,
	                                               B_x3,B_y3,B_z3);
	                        pH_x[j+2] = H_x3;
	                        pH_y[j+2] = H_y3;
	                        pH_z[j+2] = H_z3;
	                        pB_x[j+2] = B_x3;
	                        pB_y[j+2] = B_y3;
	                        pB_z[j+2] = B_z3;
	                        tht4   = ptht[j+3];
	                        phi4   = pphi[j+3];
	                        omg4   = pomg[j+3];
	                        phase4 = pphase[j+3];
	                        refi4  = prefi[j+3];
	                        px4    = ppx[j+3];
	                        py4    = ppy[j+3];
	                        pz4    = ppz[j+3];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht4,phi4,omg4,
	                                               phase4,refi4,px4,
	                                               py4,pz4,
	                                               H_x4,H_y4,H_z4,
	                                               B_x4,B_y4,B_z4);
	                        pH_x[j+3] = H_x4;
	                        pH_y[j+3] = H_y4;
	                        pH_z[j+3] = H_z4;
	                        pB_x[j+3] = B_x4;
	                        pB_y[j+3] = B_y4;
	                        pB_z[j+3] = B_z4;  
	                        tht5   = ptht[j+4];
	                        phi5   = pphi[j+4];
	                        omg5   = pomg[j+4];
	                        phase5 = pphase[j+4];
	                        refi5  = prefi[j+4];
	                        px5    = ppx[j+4];
	                        py5    = ppy[j+4];
	                        pz5    = ppz[j+4];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht5,phi5,omg5,
	                                               phase5,refi5,px5,
	                                               py5,pz5,
	                                               H_x5,H_y5,H_z5,
	                                               B_x5,B_y5,B_z);
	                        pH_x[j+4] = H_x5;
	                        pH_y[j+4] = H_y;
	                        pH_z[j+4] = H_z5;
	                        pB_x[j+4] = B_x5;
	                        pB_y[j+4] = B_y5;
	                        pB_z[j+4] = B_z5; 
	                        tht6   = ptht[j+5];
	                        phi6   = pphi[j+5];
	                        omg6   = pomg[j+5];
	                        phase6 = pphase[j+5];
	                        refi6  = prefi[j+5];
	                        px6    = ppx[j+5];
	                        py6    = ppy[j+5];
	                        pz6    = ppz[j+5];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht6,phi6,omg6,
	                                               phase6,refi6,px6,
	                                               py6,pz6,
	                                               H_x6,H_y6,H_z6,
	                                               B_x6,B_y6,B_z6);
	                        pH_x[j+5] = H_x6;
	                        pH_y[j+5] = H_y6;
	                        pH_z[j+5] = H_z6;
	                        pB_x[j+5] = B_x6;
	                        pB_y[j+5] = B_y6;
	                        pB_z[j+5] = B_z6;  
	                        tht7   = ptht[j+6];
	                        phi7   = pphi[j+6];
	                        omg7   = pomg[j+6];
	                        phase7 = pphase[j+6];
	                        refi7  = prefi[j+6];
	                        px7    = ppx[j+6];
	                        py7    = ppy[j+6];
	                        pz7    = ppz[j+6];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht7,phi7,omg7,
	                                               phase7,refi7,px7,
	                                               py7,pz7,
	                                               H_x7,H_y7,H_z7,
	                                               B_x7,B_y7,B_z7);
	                        pH_x[j+6] = H_x7;
	                        pH_y[j+6] = H_y7;
	                        pH_z[j+6] = H_z7;
	                        pB_x[j+6] = B_x7;
	                        pB_y[j+6] = B_y7;
	                        pB_z[j+6] = B_z7;  
	                       
	                 }               
	      }
	        
	      
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_ymm4c8_unroll2x_omp(const __m256d * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256d * __restrict __ATTR_ALIGN__(32) pphi,
	                                                      const __m256d * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm4c8_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm4c8_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST) {
	                                                 
	                 if(__builtin_expect(n<=0,0)) {return;}
	                 if(__builtin_expect(PF_DIST<=0,0)) PF_DIST = 2;
	                 __ATTR_ALIGN__(32) ymm4c8_t phase1,phase2;
	                 __ATTR_ALIGN__(32) ymm4c8_t refi1,refi2;
	                 __ATTR_ALIGN__(32) ymm4c8_t px1,px2;
	                 __ATTR_ALIGN__(32) ymm4c8_t py1,py2;
	                 __ATTR_ALIGN__(32) ymm4c8_t pz1,pz2;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_x1,H_x2;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_y1,H_y2;
	                 __ATTR_ALIGN__(32) ymm4c8_t H_z1,H_z2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_x1,B_x2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_y1,B_y2;
	                 __ATTR_ALIGN__(32) ymm4c8_t B_z1,B_z2;
	                 register __m256d tht1,tht2;
	                 register __m256d phi1,phi2;
	                 register __m256d omg1,omg2;
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
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht,phi,omg,
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
#pragma omp parallel for schedule(runtime) default(none)           \
        firstprivate(m1,PF_DIST) private(j,phase1,phase2)          \
                                 private(refi1,refi2)              \
                                 private(px1,px2)                  \
                                 private(py1,py2)                  \
                                 private(pz1,pz2)                  \
                                 private(H_x1,H_x2)                \
                                 private(H_y1,H_y2)                \
                                 private(H_z1,H_z2)                \
                                 private(B_x1,B_x2)                \
                                 private(B_y1,B_y2)                \
                                 private(B_z1,B_z2)                \
                                 private(tht1,tht2)                \
                                 private(phi1,phi2)                \
                                 private(omg1,omg2)                \
                                 shared(n,ptht,pphi,pomg,pphase,prefi)                 \
                                 shared(ppx,ppy,ppz)                                   \
                                 shared(pH_x,pH_y,pH_z,pB_x,pB_y,pB_z)
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
                                tht1   = ptht[j+0];
	                        phi1   = pphi[j+0];
	                        omg1   = pomg[j+0];
	                        phase1 = pphase[j+0];
	                        refi1  = prefi[j+0];
	                        px1    = ppx[j+0];
	                        py1    = ppy[j+0];
	                        pz1    = ppz[j+0];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht1,phi1,omg1,
	                                               phase1,refi1,px1,
	                                               py1,pz1,
	                                               H_x1,H_y1,H_z1,
	                                               B_x1,B_y1,B_z1);
	                        pH_x[j+0] = H_x1;
	                        pH_y[j+0] = H_y1;
	                        pH_z[j+0] = H_z1;
	                        pB_x[j+0] = B_x1;
	                        pB_y[j+0] = B_y1;
	                        pB_z[j+0] = B_z1;
	                        tht2   = ptht[j+1];
	                        phi2   = pphi[j+1];
	                        omg2   = pomg[j+1];
	                        phase2 = pphase[j+1];
	                        refi2  = prefi[j+1];
	                        px2    = ppx[j+1];
	                        py2    = ppy[j+1];
	                        pz2    = ppz[j+1];
	                        B_XYZ_H_XYZ_EP_ymm4c8(tht2,phi2,omg2,
	                                               phase2,refi2,px2,
	                                               py2,pz2,
	                                               H_x2,H_y2,H_z2,
	                                               B_x2,B_y2,B_z2);
	                        pH_x[j+1] = H_x2;
	                        pH_y[j+1] = H_y2;
	                        pH_z[j+1] = H_z2;
	                        pB_x[j+1] = B_x2;
	                        pB_y[j+1] = B_y2;
	                        pB_z[j+1] = B_z2;  
	                       
	                 }               
	      }
	  
	      
	         
	        
	      
	      
	      
	        
	      
	        
	        
	        
	                                      
	       
	       
	       
	     
                
                
        } // radiolocation

} // gms


#endif /*__GMS_EM_FIELDS_YMM4R8_OMP_HPP__*/
