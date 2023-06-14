

#ifndef __GMS_EM_FIELDS_ZMM16R4_HPP__
#define __GMS_EM_FIELDS_ZMM16R4_HPP__

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
    const char * const GMS_VECMATH_ZMM16R4_DESCRIPTION   = " Computational ElectroMagnetics related helper routines."
                       

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_complex_zmm16r4.hpp"

namespace gms {



          namespace radiolocation {
          
          
                   
                
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
	           void sdotv_zmm16r4_unroll16x(const __m512 * __restrict __ATTR_ALIGN__(64) pv1x,
	                                        const __m512 * __restrict __ATTR_ALIGN__(64) pv1y,
	                                        const __m512 * __restrict __ATTR_ALIGN__(64) pv1z,
	                                        const __m512 * __restrict __ATTR_ALIGN__(64) pv2x,
	                                        const __m512 * __restrict __ATTR_ALIGN__(64) pv2y,
	                                        const __m512 * __restrict __ATTR_ALIGN__(64) pv2z,
	                                        __m512 * __restrict __ATTR_ALIGN__(64) pdtv,
	                                        const int32_t n) {
	                                        
	                if(__builtin_expect(n<=0,0)) { return;}
	                register __m512 v1x;
	                register __m512 v1y;
	                register __m512 v1z;
	                register __m512 v2x;
	                register __m512 v2y;
	                register __m512 v2z;
	                register __m512 dtv;
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
	                        dtv = sdotv_zmm16r4(v1x,v1y,v1z,
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
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+0] = dtv;
	                     v1x = pv1x[j+1];
	                     v2x = pv2x[j+1];
	                     v1y = pv1y[j+1];
	                     v2y = pv2y[j+1];
	                     v1z = pv1z[j+1];
	                     v2z = pv2z[j+1];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+1] = dtv;
	                     v1x = pv1x[j+2];
	                     v2x = pv2x[j+2];
	                     v1y = pv1y[j+2];
	                     v2y = pv2y[j+2];
	                     v1z = pv1z[j+2];
	                     v2z = pv2z[j+2];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+2] = dtv;
	                     v1x = pv1x[j+3];
	                     v2x = pv2x[j+3];
	                     v1y = pv1y[j+3];
	                     v2y = pv2y[j+3];
	                     v1z = pv1z[j+3];
	                     v2z = pv2z[j+3];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+3] = dtv;
	                     v1x = pv1x[j+4];
	                     v2x = pv2x[j+4];
	                     v1y = pv1y[j+4];
	                     v2y = pv2y[j+4];
	                     v1z = pv1z[j+4];
	                     v2z = pv2z[j+4];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+4] = dtv;
	                     v1x = pv1x[j+5];
	                     v2x = pv2x[j+5];
	                     v1y = pv1y[j+5];
	                     v2y = pv2y[j+5];
	                     v1z = pv1z[j+5];
	                     v2z = pv2z[j+5];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+5] = dtv;
	                     v1x = pv1x[j+6];
	                     v2x = pv2x[j+6];
	                     v1y = pv1y[j+6];
	                     v2y = pv2y[j+6];
	                     v1z = pv1z[j+6];
	                     v2z = pv2z[j+6];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+6] = dtv;
	                     v1x = pv1x[j+7];
	                     v2x = pv2x[j+7];
	                     v1y = pv1y[j+7];
	                     v2y = pv2y[j+7];
	                     v1z = pv1z[j+7];
	                     v2z = pv2z[j+7];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+7] = dtv;
	                     v1x = pv1x[j+8];
	                     v2x = pv2x[j+8];
	                     v1y = pv1y[j+8];
	                     v2y = pv2y[j+8];
	                     v1z = pv1z[j+8];
	                     v2z = pv2z[j+8];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+8] = dtv;
	                     v1x = pv1x[j+9];
	                     v2x = pv2x[j+9];
	                     v1y = pv1y[j+9];
	                     v2y = pv2y[j+9];
	                     v1z = pv1z[j+9];
	                     v2z = pv2z[j+9];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+9] = dtv;
	                     v1x = pv1x[j+10];
	                     v2x = pv2x[j+10];
	                     v1y = pv1y[j+10];
	                     v2y = pv2y[j+10];
	                     v1z = pv1z[j+10];
	                     v2z = pv2z[j+10];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+10] = dtv;
	                     v1x = pv1x[j+11];
	                     v2x = pv2x[j+11];
	                     v1y = pv1y[j+11];
	                     v2y = pv2y[j+11];
	                     v1z = pv1z[j+11];
	                     v2z = pv2z[j+11];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+11] = dtv;
	                     v1x = pv1x[j+12];
	                     v2x = pv2x[j+12];
	                     v1y = pv1y[j+12];
	                     v2y = pv2y[j+12];
	                     v1z = pv1z[j+12];
	                     v2z = pv2z[j+12];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+12] = dtv;
	                     v1x = pv1x[j+13];
	                     v2x = pv2x[j+13];
	                     v1y = pv1y[j+13];
	                     v2y = pv2y[j+13];
	                     v1z = pv1z[j+13];
	                     v2z = pv2z[j+13];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+13] = dtv;
	                     v1x = pv1x[j+14];
	                     v2x = pv2x[j+14];
	                     v1y = pv1y[j+14];
	                     v2y = pv2y[j+14];
	                     v1z = pv1z[j+14];
	                     v2z = pv2z[j+14];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+14] = dtv;
	                     v1x = pv1x[j+15];
	                     v2x = pv2x[j+15];
	                     v1y = pv1y[j+15];
	                     v2y = pv2y[j+15];
	                     v1z = pv1z[j+15];
	                     v2z = pv2z[j+15];
	                     dtv = sdotv_zmm16r4(v1x,v1y,v1z,
	                                         v2x,v2y,v2z);
	                     pdtv[j+15] = dtv;
	                }
	                                        
	      }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 sdotv_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pv1x,
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
	           __m512 sdotv_zmm16r4_u(const float * __restrict  pv1x,
	                                const float * __restrict  pv1y,
	                                const float * __restrict  pv1z,
	                                const float * __restrict  pv2x,
	                                const float * __restrict  pv2y,
	                                const float * __restrict  pv2z) {
	                          
	                  register __m512 v1x = _mm512_loadu_ps(&pv1x[0]);
	                  register __m512 v1y = _mm512_loadu_ps(&pv1y[0]);  
	                  register __m512 v1z = _mm512_loadu_ps(&pv1z[0]); 
	                  register __m512 v2x = _mm512_loadu_ps(&pv2x[0]);  
	                  register __m512 v2y = _mm512_loadu_ps(&pv2y[0]); 
	                  register __m512 v2z = _mm512_loadu_ps(&pv2z[0]);
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
	           void cdotv_zmm16c4(const zmm16c4_t v1x,
	                              const zmm16c4_t v1y,
	                              const zmm16c4_t v1z,
	                              const zmm16c4_t v2x,
	                              const zmm16c4_t v2y,
	                              const zmm16c4_t v2z,
	                              zmm16c4_t & res) {
	                              
	                zmm16c4_t tx,ty,tz;
	                tx = cmul_zmm16r4(v1x.re,v1x.im,v2x.re,
	                                  v2x.im,&tx.re,&tx.im); 
	                ty = cmul_zmm16r4(v1y.re,v1y.im,v2y.re,
	                                  v2y.im,&ty.re,&ty.im);
	                tz = cmul_zmm16r4(v1z.re,v1z.im,v2z.re,
	                                  v2z.im,&tz.re,&tz.im);
	                res.re = _mm512_add_ps(tx.re,
	                                   _mm512_add_ps(ty.re,tz.re));
	                res.im = _mm512_add_ps(tx.im,
	                                   _mm512_add_ps(ty.im,tz.im));                   
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 cnorm_zmm16c4(const zmm16c4_t vx,
	                                const zmm16c4_t vy,
	                                const zmm16c4_t vz) {
	                                
	                  zmm16c4_t t,cx,cy,cz;
	                  __m512 vs;
	                  cconj_zmm16r4_v2(vx.re,vx.im,&cx.re,&cx.im);
	                  cconj_zmm16r4_v2(vy.re,vy.im,&cy.re,&cy.im);
	                  cconj_zmm16r4_v2(vz.re,vz.im,&cz.re,&cz.im);
	                  cdotv_zmm16c4(vx,vy,vz,cx,cy,cz,t);
	                  vs = _mm512_sqrt_ps(t.re);
	                  return (vs);                      
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossc_zmm16c4(const zmm16c4_t v1x,
	                                const zmm16c4_t v1y,
	                                const zmm16c4_t v1z,
	                                const zmm16c4_t v2x,
	                                const zmm16c4_t v2y,
	                                const zmm16c4_t v2z,
	                                zmm16c4 & resx,
	                                zmm16c4 & resy,
	                                zmm16c4 & resz) {
	                                
	                 zmm16c4_t t0,t1,t2,t3,t4,t5,t6;
	                 cmul_zmm16r4(v1y.re,v1y.im,v2z.re,
	                              v2z.im,&t0.re,&t0.im); 
	                 cmul_zmm16r4(v1z.re,v1z.im,v2y.re,
	                              v2y.im,&t1.re,&t1.im);
	                 resx.re = _mm512_sub_ps(t0.re,t1.re);
	                 resx.im = _mm512_sub_ps(t0.im,t1.im);
	                 cmul_zmm16r4(v1z.re,v1z.im,v2x.re,
	                              v2x.im,&t2.re,&t2.im);
	                 cmul_zmm16r4(v1x.re,v1x.im,v2z.re,
	                              v2z.im,&t3.re,&t3.im);
	                 resy.re = _mm512_sub_ps(t2.re,t3.re);
	                 resy.im = _mm512_sub_ps(t2.im,t3.im);
	                 cmul_zmm16r4(v1x.re,v1x.im,v2y.re,
	                              v2y.im,&t4.re,&t4.im);
	                 cmul_zmm16r4(v1y.re,v1y.im,v2x.re,
	                              v2x.im,&t5.re,&t5.im);    
	                 resz.re = _mm512_sub_ps(t4.re,t5.re);
	                 resz.im = _mm512_sub_ps(t4.im,t5.im);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_zmm16r4(const __m512 v1x,
	                                const __m512 v1y,
	                                const __m512 v1z,
	                                const __m512 v2x,
	                                const __m512 v2y,
	                                const __m512 v2z,
	                                __m512 * __restrict vcx,
	                                __m512 * __restrict vcy,
	                                __m512 * __restrict vcz) {
	                                
	                *vcx = _mm512_fmsub_ps(v1y,v2z,
	                                   _mm512_mul_ps(v1x,v2y));
	                *vcy = _mm512_fmsub_ps(v1z,v2x,
	                                   _mm512_mul_ps(v1x,v2z));
	                *vcz = _mm512_fmsub_ps(v1x,v2y,
	                                   _mm512_mul_ps(v1y,v2x));
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pv1x,
	                                  const float * __restrict __ATTR_ALIGN__(64) pv1y,
	                                  const float * __restrict __ATTR_ALIGN__(64) pv1z,
	                                  const float * __restrict __ATTR_ALIGN__(64) pv2x,
	                                  const float * __restrict __ATTR_ALIGN__(64) pv2y,
	                                  const float * __restrict __ATTR_ALIGN__(64) pv2z,
	                                  float * __restrict __ATTR_ALIGN__(64) vcx,
	                                  float * __restrict __ATTR_ALIGN__(64) vcy,
	                                  float * __restrict __ATTR_ALIGN__(64) vcz) {
	                      
	                 register __m512 v1x = _mm512_load_ps(&pv1x[0]);
	                 register __m512 v1y = _mm512_load_ps(&pv1y[0]);
	                 register __m512 v1z = _mm512_load_ps(&pv1z[0]);
	                 register __m512 v2x = _mm512_load_ps(&pv2x[0]);
	                 register __m512 v2y = _mm512_load_ps(&pv2y[0]);
	                 register __m512 v2z = _mm512_load_ps(&pv2z[0]);          
	                *vcx = _mm512_fmsub_ps(v1y,v2z,
	                                   _mm512_mul_ps(v1x,v2y));
	                *vcy = _mm512_fmsub_ps(v1z,v2x,
	                                   _mm512_mul_ps(v1x,v2z));
	                *vcz = _mm512_fmsub_ps(v1x,v2y,
	                                   _mm512_mul_ps(v1y,v2x));
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void scrossv_zmm16r4_u(const float * __restrict pv1x,
	                                  const float * __restrict pv1y,
	                                  const float * __restrict pv1z,
	                                  const float * __restrict pv2x,
	                                  const float * __restrict pv2y,
	                                  const float * __restrict pv2z,
	                                  float * __restrict vcx,
	                                  float * __restrict vcy,
	                                  float * __restrict vcz) {
	                      
	                 register __m512 v1x = _mm512_loadu_ps(&pv1x[0]);
	                 register __m512 v1y = _mm512_loadu_ps(&pv1y[0]);
	                 register __m512 v1z = _mm512_loadu_ps(&pv1z[0]);
	                 register __m512 v2x = _mm512_loadu_ps(&pv2x[0]);
	                 register __m512 v2y = _mm512_loadu_ps(&pv2y[0]);
	                 register __m512 v2z = _mm512_loadu_ps(&pv2z[0]);          
	                *vcx = _mm512_fmsub_ps(v1y,v2z,
	                                   _mm512_mul_ps(v1x,v2y));
	                *vcy = _mm512_fmsub_ps(v1z,v2x,
	                                   _mm512_mul_ps(v1x,v2z));
	                *vcz = _mm512_fmsub_ps(v1x,v2y,
	                                   _mm512_mul_ps(v1y,v2x));
	         }
	         
	         
	         //! Direction Vector spherical [theta,phi] (SIMD data-types)
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_zmm16r4(  const __m512 tht,
	                                  const __m512 phi,
	                                  __m512 * __restrict dvx,
	                                  __m512 * __restrict dvy,
	                                  __m512 * __restrict dvz) {
	                  
	                        
	                register __m512 stht,cphi,sphi,ctht;
	                cphi = xcosf(phi);
	                stht = xsinf(tht);
	                *dvx = _mm512_mul_ps(stht,cphi);
	                sphi = xsinf(phi);
	                *dvy = _mm512_mul_ps(stht,sphi);
	                ctht = xcosf(tht);
	                *dvz = ctht;                       
	        }
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                  const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                  float * __restrict __ATTR_ALIGN__(64) dvx,
	                                  float * __restrict __ATTR_ALIGN__(64) dvy,
	                                  float * __restrict __ATTR_ALIGN__(64) dvz) {
	                  
	                register __m512 tht = _mm512_load_ps(&ptht[0]);
	                register __m512 phi = _mm512_load_ps(&pphi[0]);              
	                register __m512 stht,cphi,sphi,ctht;
	                cphi = xcosf(phi);
	                stht = xsinf(tht);
	                _mm512_store_ps(&dvx[0] , _mm512_mul_ps(stht,cphi));
	                sphi = xsinf(phi);
	                _mm512_store_ps(&dvy[0] , _mm512_mul_ps(stht,sphi));
	                ctht = xcosf(tht);
	                _mm512_store_ps(&dvz[0] , ctht);                       
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void dir_vec_zmm16r4_u(const float * __restrict  ptht,
	                                  const float * __restrict  pphi,
	                                  float * __restrict  dvx,
	                                  float * __restrict  dvy,
	                                  float * __restrict  dvz) {
	                  
	                register __m512 tht = _mm512_loadu_ps(&ptht[0]);
	                register __m512 phi = _mm512_loadu_ps(&pphi[0]);              
	                register __m512 stht,cphi,sphi,ctht;
	                cphi = xcosf(phi);
	                stht = xsinf(tht);
	                _mm512_storeu_ps(&dvx[0] , _mm512_mul_ps(stht,cphi));
	                sphi = xsinf(phi);
	                _mm512_storeu_ps(&dvy[0] , _mm512_mul_ps(stht,sphi));
	                ctht = xcosf(tht);
	                _mm512_storeu_ps(&dvz[0] , ctht);                       
	        }
	        
	        
	         //! Polarization Vector of plane-wave propagating into direction computed by
                 //! dir_vector_xmmxrx (SIMD data-types)
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_zmm16r4(const __m512 tht,
	                                const __m512 phi,
	                                const __m512 psi,
	                                __m512 * __restrict pvx,
	                                __m512 * __restrict pvy,
	                                __m512 * __restrict pvz) {
	                 
	                using namespace gms::math               
	                register __m512 cpsi,cphi,spsi,sphi,t0;
	                cpsi = xcosf(psi);
	                cphi = xcosf(phi);
	                spsi = xsinf(psi);
	                sphi = xsinf(phi);
	                t0   = _mm512_mul_ps(spsi,xcosf(tht));
	                *pvx = _mm512_fmsub_ps(cpsi,sphi,
	                                   _mm512_mul_ps(t0,cphi));
	                *pvy = _mm512_fmsub_ps(negate_zmm16r4(cpsi),cphi,
	                                                    _mm512_mul_ps(t0,sphi));
	                *pvz = _mm512_mul_ps(spsi,xsinf(tht));                         
	      }
	      
	      
	      
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                  const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                  const float * __restrict __ATTR_ALIGN__(64) psi,
	                                  float * __restrict __ATTR_ALIGN__(64) pvx,
	                                  float * __restrict __ATTR_ALIGN__(64) pvy,
	                                  float * __restrict __ATTR_ALIGN__(64) pvz) {
	                 
	                 using namespace gms::math     
	                register __m512 tht = _mm512_load_ps(&ptht[0]);
	                register __m512 phi = _mm512_load_ps(&pphi[0]);  
	                register __m512 psi = _mm512_load_ps(&ppsi[0]);           
	                register __m512 cpsi,cphi,spsi,sphi,t0;
	                cpsi = xcosf(psi);
	                cphi = xcosf(phi);
	                spsi = xsinf(psi);
	                sphi = xsinf(phi);
	                t0   = _mm512_mul_ps(spsi,xcosf(tht));
	                _mm512_store_ps(&pvx[0] ,_mm512_fmsub_ps(cpsi,sphi,
	                                   _mm512_mul_ps(t0,cphi)));
	                _mm512_store_ps(&pvy[0] ,_mm512_fmsub_ps(negate_zmm16r4(cpsi),cphi,
	                                                    _mm512_mul_ps(t0,sphi)));
	                _mm512_store_ps(&pvz[0] ,_mm512_mul_ps(spsi,xsinf(tht)));                         
	      } 
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void pol_vec_zmm16r4_u(const float * __restrict  ptht,
	                                  const float * __restrict  pphi,
	                                  const float * __restrict  psi,
	                                  float * __restrict  pvx,
	                                  float * __restrict  pvy,
	                                  float * __restrict  pvz) {
	                 
	                  using namespace gms::math    
	                register __m512 tht = _mm512_loadu_ps(&ptht[0]);
	                register __m512 phi = _mm512_loadu_ps(&pphi[0]);  
	                register __m512 psi = _mm512_loadu_ps(&ppsi[0]);           
	                register __m512 cpsi,cphi,spsi,sphi,t0;
	                cpsi = xcosf(psi);
	                cphi = xcosf(phi);
	                spsi = xsinf(psi);
	                sphi = xsinf(phi);
	                t0   = _mm512_mul_ps(spsi,xcosf(tht));
	                _mm512_storeu_ps(&pvx[0] ,_mm512_fmsub_ps(cpsi,sphi,
	                                   _mm512_mul_ps(t0,cphi)));
	                _mm512_storeu_ps(&pvy[0] ,_mm512_fmsub_ps(negate_zmm16r4(cpsi),cphi,
	                                                    _mm512_mul_ps(t0,sphi)));
	                _mm512_storeu_ps(&pvz[0] ,_mm512_mul_ps(spsi,xsinf(tht)));                         
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
	           void H_XYZ_VP_zmm16c4(const __m512 vpolx,
	                                 const __m512 vpoly,
	                                 const __m512 vpolz,
	                                 const __m512 vdirx,
	                                 const __m512 vdiry,
	                                 const __m512 vdirz,
	                                 const __m512 vrx,
	                                 const __m512 vry,
	                                 const __m512 vrz,
	                                 const zmm16c4_t k,
	                                 zmm16c4_t & H_x,
	                                 zmm16c4_t & H_y,
	                                 zmm16c4_t & H_z) {
	               
	               	__m512 dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_zmm16r4(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm512_set1_ps(1.0f);
	                ir = _mm512_setzero_ps();
	                cmul_zmm16r4(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm512_mul_ps(dp,cer);
	                cei = _mm512_mul_ps(dp,cei);
	                cexp_zmm16r4(cer,cei,&expr,&expi);
	                H_x.re = _mm512_mul_ps(vpolx,expr);
	                H_x.im = _mm512_mul_ps(vpolx,expi);
	                H_y.re = _mm512_mul_ps(vpoly,expr);
	                H_y.im = _mm512_mul_ps(vpoly,expi);
	                H_z.re = _mm512_mul_ps(vpolz,expr);
	                H_z.im = _mm512_mul_ps(vpolz,expi);
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_zmm16c4_a(const float * __restrict __ATTR_ALIGN__(64) vpolx,
	                                 const float * __restrict __ATTR_ALIGN__(64) vpoly,
	                                 const float * __restrict __ATTR_ALIGN__(64) vpolz,
	                                 const float * __restrict __ATTR_ALIGN__(64) vdirx,
	                                 const float * __restrict __ATTR_ALIGN__(64) vdiry,
	                                 const float * __restrict __ATTR_ALIGN__(64) vdirz,
	                                 const float * __restrict __ATTR_ALIGN__(64) vrx,
	                                 const float * __restrict __ATTR_ALIGN__(64) vry,
	                                 const float * __restrict __ATTR_ALIGN__(64) vrz,
	                                 const zmm16c4_t k,
	                                 zmm16c4_t & H_x,
	                                 zmm16c4_t & H_y,
	                                 zmm16c4_t & H_z) {
	               
	                register __m512 vpolx = _mm512_load_ps(&vpolx[0]);
	                register __m512 vpoly = _mm512_load_ps(&vpoly[0]);
	                register __m512 vpolz = _mm512_load_ps(&vpolz[0]);
	                register __m512 vdirx = _mm512_load_ps(&vdirx[0]);
	                register __m512 vdiry = _mm512_load_ps(&vdiry[0]);
	                register __m512 vdirz = _mm512_load_ps(&vdirz[0]);
	                register __m512 vrx   = _mm512_load_ps(&vrx[0]);
	                register __m512 vry   = _mm512_load_ps(&vry[0]);
	                register __m512 vrz   = _mm512_load_ps(&vrz[0]);
	               	__m512 dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_zmm16r4(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm512_set1_ps(1.0f);
	                ir = _mm512_setzero_ps();
	                cmul_zmm16r4(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm512_mul_ps(dp,cer);
	                cei = _mm512_mul_ps(dp,cei);
	                cexp_zmm16r4(cer,cei,&expr,&expi);
	                H_x.re = _mm512_mul_ps(vpolx,expr);
	                H_x.im = _mm512_mul_ps(vpolx,expi);
	                H_y.re = _mm512_mul_ps(vpoly,expr);
	                H_y.im = _mm512_mul_ps(vpoly,expi);
	                H_z.re = _mm512_mul_ps(vpolz,expr);
	                H_z.im = _mm512_mul_ps(vpolz,expi);
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void H_XYZ_VP_zmm16c4_u(const float * __restrict  vpolx,
	                                 const float * __restrict  vpoly,
	                                 const float * __restrict  vpolz,
	                                 const float * __restrict  vdirx,
	                                 const float * __restrict  vdiry,
	                                 const float * __restrict  vdirz,
	                                 const float * __restrict  vrx,
	                                 const float * __restrict  vry,
	                                 const float * __restrict  vrz,
	                                 const zmm16c4_t k,
	                                 zmm16c4_t & H_x,
	                                 zmm16c4_t & H_y,
	                                 zmm16c4_t & H_z) {
	               
	                register __m512 vpolx = _mm512_loadu_ps(&vpolx[0]);
	                register __m512 vpoly = _mm512_loadu_ps(&vpoly[0]);
	                register __m512 vpolz = _mm512_loadu_ps(&vpolz[0]);
	                register __m512 vdirx = _mm512_loadu_ps(&vdirx[0]);
	                register __m512 vdiry = _mm512_loadu_ps(&vdiry[0]);
	                register __m512 vdirz = _mm512_loadu_ps(&vdirz[0]);
	                register __m512 vrx   = _mm512_loadu_ps(&vrx[0]);
	                register __m512 vry   = _mm512_loadu_ps(&vry[0]);
	                register __m512 vrz   = _mm512_loadu_ps(&vrz[0]);
	               	__m512 dp,cer,cei,ii,ir,expr,expi;
	                dp = sdotv_zmm16r4(vdirx,vdiry,vdirz,
	                                   vrx,vry,vrz);
	                ii = _mm512_set1_ps(1.0f);
	                ir = _mm512_setzero_ps();
	                cmul_zmm16r4(ir,ii,k.re,k.im,&cer,&cei);
	                cer = _mm512_mul_ps(dp,cer);
	                cei = _mm512_mul_ps(dp,cei);
	                cexp_zmm16r4(cer,cei,&expr,&expi);
	                H_x.re = _mm512_mul_ps(vpolx,expr);
	                H_x.im = _mm512_mul_ps(vpolx,expi);
	                H_y.re = _mm512_mul_ps(vpoly,expr);
	                H_y.im = _mm512_mul_ps(vpoly,expi);
	                H_z.re = _mm512_mul_ps(vpolz,expr);
	                H_z.im = _mm512_mul_ps(vpolz,expi);
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
	           void B_XYZ_VP_zmm16c4(const __m512 vpolx,
	                                 const __m512 vpoly,
	                                 const __m512 vpolz,
	                                 const __m512 vdirx,
	                                 const __m512 vdiry,
	                                 const __m512 vdirz,
	                                 const zmm16c4_t k,
	                                 const __m512 omega,
	                                 const __m512 vrx,
	                                 const __m512 vry,
	                                 const __m512 vrz,
	                                 zmm16c4_t & B_x,
	                                 zmm16c4_t & B_y,
	                                 zmm16c4_t & B_z) {
	                                 
	                const __m512 mu0 = _mm512_set1_ps(0.0000012566370614359173f);
	                zmm16c4_t cdirx;
	                zmm16c4_t cdiry;
	                zmm16c4_t cdirz;
	                zmm16c4_t H_x;
	                zmm16c4_t H_y;
	                zmm16c4_t H_z;
	                zmm16c4_t cpx;
	                zmm16c4_t cpy;
	                zmm16c4_t cpz;
	                zmm16c4_t t0;
	                __m512 zz0;
	                H_XYZ_VP_zmm16c4(vpolx,vpoy,vpolz,
	                               	 vdirx,vdiry,vdirz,
	                                 vrx,vry,vrz,
	                                 H_x,H_y,H_z);
	                                 	                                 
	                cdirx.re = vdirx;
	                cdirx.im = _mm512_setzero_ps();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm512_mul_ps(omega,mu0);
	                t0.re    = _mm512_div_ps(k.re,zz0);
	                t0.im    = _mm512_div_ps(k.im,zz0);
	                
	                scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                                H_x,H_y,H_z,
                                        cpx,cpy,cpz);                
	                     	                                
	                cmul_zmm16r4(t0.re,t0.im,
	                             cpx.re,cpx.im,
	                             &B_x.re,&B_x.im);
	                             
	                cmul_zmm16r4(t0.re,t0.im,
	                             cpy.re,cpy.im,
	                             &B_y.re,&B_y.im);
	                            
	                cmul_zmm16r4(t0.re,t0.im,
	                             cpz.re,cpz.im,
	                             &B_z.re,&B_z.im);
	                          
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_zmm16c4_a(const float * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                   const float * __restrict __ATTR_ALIGN__(64) pomega,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvrx,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvry,
	                                   const float * __restrict __ATTR_ALIGN__(64) pvrz,
	                                   const zmm16c4_t k,
	                                   zmm16c4_t & B_x,
	                                   zmm16c4_t & B_y,
	                                   zmm16c4_t & B_z) {
	                         
	                register __m512 vpolx = _mm512_load_ps(&pvpolx[0]);
	                register __m512 vpoly = _mm512_load_ps(&pvpoly[0]);    
	                register __m512 vpolz = _mm512_load_ps(&pvpolz[0]);  
	                register __m512 vdirx = _mm512_load_ps(&pvdirx[0]);  
	                register __m512 vdiry = _mm512_load_ps(&pvdiry[0]);
	                register __m512 vdirz = _mm512_load_ps(&pvdirz[0]); 
	                register __m512 onega = _mm512_load_ps(&pomega[0]);
	                register __m512 vrx   = _mm512_load_ps(&pvrx[0]);
	                register __m512 vry   = _mm512_load_ps(&pvry[0]);
	                register __m512 vrz   = _mm512_load_ps(&pvrz[0]);        
	                const __m512 mu0 = _mm512_set1_ps(0.0000012566370614359173f);
	                zmm16c4_t cdirx;
	                zmm16c4_t cdiry;
	                zmm16c4_t cdirz;
	                zmm16c4_t H_x;
	                zmm16c4_t H_y;
	                zmm16c4_t H_z;
	                zmm16c4_t cpx;
	                zmm16c4_t cpy;
	                zmm16c4_t cpz;
	                zmm16c4_t t0;
	                __m512 zz0;
	                H_XYZ_VP_zmm16c4(vpolx,
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
	                cdirx.im = _mm512_setzero_ps();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm512_mul_ps(omega,mu0);
	                t0.re    = _mm512_div_ps(k.re,zz0);
	                t0.im    = _mm512_div_ps(k.im,zz0);
	                scrossc_zmm16c4(cdirx,
	                                cdiry,
	                                cdirz,
	                                H_x,
	                                H_y,
	                                H_z,
	                                cpx,
	                                cpy,
	                                cpz);
	                                
	                cmul_zmm16r4(t0.re,
	                             t0.im,
	                             cpx.re,
	                             cpx.im,
	                             &B_x.re,
	                             &B_x.im);
	                cmul_zmm16r4(t0.re,
	                             t0.im,
	                             cpy.re,
	                             cpy.im,
	                             &B_y.re,
	                             &B_y.im);
	                cmul_zmm16r4(t0.re,
	                             t0.im,
	                             cpz.re,
	                             cpz.im,
	                             &B_z.re,
	                             &B_z.im);
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_VP_zmm16c4_u(const float * __restrict  pvpolx,
	                                   const float * __restrict  pvpoly,
	                                   const float * __restrict  pvpolz,
	                                   const float * __restrict  pvdirx,
	                                   const float * __restrict  pvdiry,
	                                   const float * __restrict  pvdirz,
	                                   const float * __restrict  pomega,
	                                   const float * __restrict  pvrx,
	                                   const float * __restrict  pvry,
	                                   const float * __restrict  pvrz,
	                                   const zmm16c4_t k,
	                                   zmm16c4_t & B_x,
	                                   zmm16c4_t & B_y,
	                                   zmm16c4_t & B_z) {
	                         
	                register __m512 vpolx = _mm512_loadu_ps(&pvpolx[0]);
	                register __m512 vpoly = _mm512_loadu_ps(&pvpoly[0]);    
	                register __m512 vpolz = _mm512_loadu_ps(&pvpolz[0]);  
	                register __m512 vdirx = _mm512_loadu_ps(&pvdirx[0]);  
	                register __m512 vdiry = _mm512_loadu_ps(&pvdiry[0]);
	                register __m512 vdirz = _mm512_loadu_ps(&pvdirz[0]); 
	                register __m512 onega = _mm512_loadu_ps(&pomega[0]);
	                register __m512 vrx   = _mm512_loadu_ps(&pvrx[0]);
	                register __m512 vry   = _mm512_loadu_ps(&pvry[0]);
	                register __m512 vrz   = _mm512_loadu_ps(&pvrz[0]);        
	                const __m512 mu0 = _mm512_set1_ps(0.0000012566370614359173f);
	                zmm16c4_t cdirx;
	                zmm16c4_t cdiry;
	                zmm16c4_t cdirz;
	                zmm16c4_t H_x;
	                zmm16c4_t H_y;
	                zmm16c4_t H_z;
	                zmm16c4_t cpx;
	                zmm16c4_t cpy;
	                zmm16c4_t cpz;
	                zmm16c4_t t0;
	                __m512 zz0;
	                H_XYZ_VP_zmm16c4(vpolx,
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
	                cdirx.im = _mm512_setzero_ps();
	                cdiry.re = vdiry;
	                cdiry.im = cdirx.im;
	                cdirz.re = vdirz;
	                cdirz.im = cdirx.im;
	                
	                zz0      = _mm512_mul_ps(omega,mu0);
	                t0.re    = _mm512_div_ps(k.re,zz0);
	                t0.im    = _mm512_div_ps(k.im,zz0);
	                scrossc_zmm16c4(cdirx,
	                                cdiry,
	                                cdirz,
	                                H_x,
	                                H_y,
	                                H_z,
	                                cpx,
	                                cpy,
	                                cpz);
	                                
	                cmul_zmm16r4(t0.re,
	                             t0.im,
	                             cpx.re,
	                             cpx.im,
	                             &B_x.re,
	                             &B_x.im);
	                cmul_zmm16r4(t0.re,
	                             t0.im,
	                             cpy.re,
	                             cpy.im,
	                             &B_y.re,
	                             &B_y.im);
	                cmul_zmm16r4(t0.re,
	                             t0.im,
	                             cpz.re,
	                             cpz.im,
	                             &B_z.re,
	                             &B_z.im);
	                                           
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_zmm16c4(const __m512 tht,
	                                      const __m512 phi,
	                                      const __m512 psi,
	                                      const __m512 omg,
	                                      const __m512 px,
	                                      const __m512 py,
	                                      const __m512 pz,
	                                      const zmm16c4_t r,
	                                      zmm16c4_t & H_x,
	                                      zmm16c4_t & H_y,
	                                      zmm16c4_t & H_z,
	                                      zmm16c4_t & B_x,
	                                      zmm16c4_t & B_y,
	                                      zmm16c4_t & B_z) {
	                                      
	                
	                const __m512 c = _mm512_set1_ps(299792458.0f);
	                zmm16c4_t k;
	                register __m512 vpolx,vpoly,vpolz;
	                register __m512 vdirx,vdiry,vdirz;
	                register __m512 t0;
	                
	                t0 = _mm512_div_ps(omg,c);
	                dir_vec_zmm16r4(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm512_mul_ps(r.re,c);
	                k.im = _mm512_mul_ps(r.im,c);
	                pol_vec_zmm16r4(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_zmm16c4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                      const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                      const float * __restrict __ATTR_ALIGN__(64) ppsi,
	                                      const float * __restrict __ATTR_ALIGN__(64) pomg,
	                                      const float * __restrict __ATTR_ALIGN__(64) ppx,
	                                      const float * __restrict __ATTR_ALIGN__(64) ppy,
	                                      const float * __restrict __ATTR_ALIGN__(64) ppz,
	                                      const zmm16c4_t r,
	                                      zmm16c4_t & H_x,
	                                      zmm16c4_t & H_y,
	                                      zmm16c4_t & H_z,
	                                      zmm16c4_t & B_x,
	                                      zmm16c4_t & B_y,
	                                      zmm16c4_t & B_z) {
	                                      
	                
	                register __m512 tht = _mm512_load_ps(&ptht[0]);
	                register __m512 phi = _mm512_load_ps(&pphi[0]);
	                register __m512 psi = _mm512_load_ps(&ppsi[0]);
	                register __m512 omg = _mm512_load_ps(&pomg[0]);
	                register __m512 px  = _mm512_load_ps(&ppx[0]);
	                register __m512 py  = _mm512_load_ps(&ppy[0]);
	                register __m512 pz  = _mm512_load_ps(&ppz[0]);
	                const __m512 c = _mm512_set1_ps(299792458.0f);
	                zmm16c4_t k;
	                register __m512 vpolx,vpoly,vpolz;
	                register __m512 vdirx,vdiry,vdirz;
	                register __m512 t0;
	                
	                t0 = _mm512_div_ps(omg,c);
	                dir_vec_zmm16r4(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm512_mul_ps(r.re,c);
	                k.im = _mm512_mul_ps(r.im,c);
	                pol_vec_zmm16r4(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_P_zmm16c4_u(const float * __restrict ptht,
	                                      const float * __restrict  pphi,
	                                      const float * __restrict ppsi,
	                                      const float * __restrict  pomg,
	                                      const float * __restrict  ppx,
	                                      const float * __restrict  ppy,
	                                      const float * __restrict ppz,
	                                      const zmm16c4_t r,
	                                      zmm16c4_t & H_x,
	                                      zmm16c4_t & H_y,
	                                      zmm16c4_t & H_z,
	                                      zmm16c4_t & B_x,
	                                      zmm16c4_t & B_y,
	                                      zmm16c4_t & B_z) {
	                                      
	                
	                register __m512 tht = _mm512_loadu_ps(&ptht[0]);
	                register __m512 phi = _mm512_loadu_ps(&pphi[0]);
	                register __m512 psi = _mm512_loadu_ps(&ppsi[0]);
	                register __m512 omg = _mm512_loadu_ps(&pomg[0]);
	                register __m512 px  = _mm512_loadu_ps(&ppx[0]);
	                register __m512 py  = _mm512_loadu_ps(&ppy[0]);
	                register __m512 pz  = _mm512_loadu_ps(&ppz[0]);
	                const __m512 c = _mm512_set1_ps(299792458.0f);
	                zmm16c4_t k;
	                register __m512 vpolx,vpoly,vpolz;
	                register __m512 vdirx,vdiry,vdirz;
	                register __m512 t0;
	                
	                t0 = _mm512_div_ps(omg,c);
	                dir_vec_zmm16r4(tht,phi,&vdirx,
	                                &vdiry,&vdirz);
	                k.re = _mm512_mul_ps(r.re,c);
	                k.im = _mm512_mul_ps(r.im,c);
	                pol_vec_zmm16r4(tht,psi,psi,
	                                &vpolx,&vpoly,&vpolz);
	                H_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,px,py,pz,
	                                 H_x,H_y,H_z);
	                B_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                 vdirx,vdiry,vdirz,
	                                 k,omg,px,py,pz,
	                                 B_x,B_y,B_z);
	                                    
	       }
	       
	       
	       /*
	              ! Electric and Magnetic Fields elliptically polarized
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_zmm16c4(const __m512 tht,
	                                       const __m512 phi,
	                                       const __m512 omg,
	                                       const zmm16c4_t phase,
	                                       const zmm16c4_t refi,
	                                       const zmm16c4_t px,
	                                       const zmm16c4_t py,
	                                       const zmm16c4_t pz,
	                                       zmm16c4_t & H_x,
	                                       zmm16c4_t & H_y,
	                                       zmm16c4_t & H_z,
	                                       zmm16c4_t & B_x,
	                                       zmm16c4_t & B_y,
	                                       zmm16c4_t & B_z) {
	                                   
	               const __m512 c   = _mm512_set1_ps(299792458.0f); 
	               const __m512 mu0 = _mm512_set1_ps(0.0000012566370614359173f);   
	               const __m512 psi0 = _mm512_setzero_ps();
	               const __m512 C00  = _mm512_setzero_ps();
	               
	               zmm16c4_t H_x_1;
	               zmm16c4_t H_y_1;
	               zmm16c4_t H_z_1;
	               zmm16c4_t H_x_2;
	               zmm16c4_t H_y_2;
	               zmm16c4_t H_z_2;
	               zmm16c4_t k;
	               zmm16c4_t t0;
	               zmm16c4_t cdirx;
	               zmm16c4_t cdiry;
	               zmm16c4_t cdirz;
	               
	               register __m512 vpolx;
	               register __m512 vpoly;
	               register __m512 vpolz;
	               register __m512 vdirx;
	               register __m512 vdiry;
	               register __m512 vdirz;
	               register __m512 cn;
	               register __m512 x0;
	               register __m512 t0r,t0i;
	               register __m512 t1r,t1i;
	               register __m512 t2r,t2i;
	               
	               dir_vec_zmm16r4(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm512_div_ps(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm512_mul_ps(refi.re,x0);
	               k.im     = _mm512_mul_ps(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_zmm16r4(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm512_add_ps(H_x_1.re,t0r);
	               H_x.im = _mm512_add_ps(H_x_1.im,t0i);
	               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm512_add_ps(H_y_1.re,t1r);
	               H_y.im = _mm512_add_ps(H_y_1.im,t1i);
	               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm512_add_ps(H_z_2.re,t2r);
	              H_z.im = _mm512_add_ps(H_z_2.im,t2i);
	              
	              cn = cnorm_zmm16c4(H_x,H_y,H_z);
	              
	              x0     = _mm512_div_ps(omg,mu0);
	              H_x.re = _mm512_div_ps(H_x.re,cn);
	              H_x.im = _mm512_div_ps(H_x.im,cn);
	              H_y.re = _mm512_div_ps(H_y.re,cn);
	              H_y.im = _mm512_div_ps(H_y.im,cn);
	              H_z.re = _mm512_div_ps(H_z.re,cn);
	              H_z.im = _mm512_div_ps(H_z.im,cn); 
	              
	              t0.re  = _mm512_div_ps(k.re,x0);
	              t0.im  = _mm512_div_ps(k.im,x0);
	              
	              scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_zmm16c4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                         const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pomg,
	                                         const zmm16c4_t phase,
	                                         const zmm16c4_t refi,
	                                         const zmm16c4_t px,
	                                         const zmm16c4_t py,
	                                         const zmm16c4_t pz,
	                                         zmm16c4_t & H_x,
	                                         zmm16c4_t & H_y,
	                                         zmm16c4_t & H_z,
	                                         zmm16c4_t & B_x,
	                                         zmm16c4_t & B_y,
	                                         zmm16c4_t & B_z) {
	                         
	               register __m512 tht = _mm512_load_ps(&ptht[0]);
	               register __m512 phi = _mm512_load_ps(&pphi[0]);
	               register __m512 omg = _mm512_load_ps(&pomg[0]);            
	               const __m512 c   = _mm512_set1_ps(299792458.0f); 
	               const __m512 mu0 = _mm512_set1_ps(0.0000012566370614359173f);   
	               const __m512 psi0 = _mm512_setzero_ps();
	               const __m512 C00  = _mm512_setzero_ps();
	               
	               zmm16c4_t H_x_1;
	               zmm16c4_t H_y_1;
	               zmm16c4_t H_z_1;
	               zmm16c4_t H_x_2;
	               zmm16c4_t H_y_2;
	               zmm16c4_t H_z_2;
	               zmm16c4_t k;
	               zmm16c4_t t0;
	               zmm16c4_t cdirx;
	               zmm16c4_t cdiry;
	               zmm16c4_t cdirz;
	               
	               register __m512 vpolx;
	               register __m512 vpoly;
	               register __m512 vpolz;
	               register __m512 vdirx;
	               register __m512 vdiry;
	               register __m512 vdirz;
	               register __m512 cn;
	               register __m512 x0;
	               register __m512 t0r,t0i;
	               register __m512 t1r,t1i;
	               register __m512 t2r,t2i;
	               
	               dir_vec_zmm16r4(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm512_div_ps(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm512_mul_ps(refi.re,x0);
	               k.im     = _mm512_mul_ps(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_zmm16r4(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm512_add_ps(H_x_1.re,t0r);
	               H_x.im = _mm512_add_ps(H_x_1.im,t0i);
	               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm512_add_ps(H_y_1.re,t1r);
	               H_y.im = _mm512_add_ps(H_y_1.im,t1i);
	               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm512_add_ps(H_z_2.re,t2r);
	              H_z.im = _mm512_add_ps(H_z_2.im,t2i);
	              
	              cn = cnorm_zmm16c4(H_x,H_y,H_z);
	              
	              x0     = _mm512_div_ps(omg,mu0);
	              H_x.re = _mm512_div_ps(H_x.re,cn);
	              H_x.im = _mm512_div_ps(H_x.im,cn);
	              H_y.re = _mm512_div_ps(H_y.re,cn);
	              H_y.im = _mm512_div_ps(H_y.im,cn);
	              H_z.re = _mm512_div_ps(H_z.re,cn);
	              H_z.im = _mm512_div_ps(H_z.im,cn); 
	              
	              t0.re  = _mm512_div_ps(k.re,x0);
	              t0.im  = _mm512_div_ps(k.im,x0);
	              
	              scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void B_XYZ_H_XYZ_EP_zmm16c4_u(const float * __restrict  ptht,
	                                         const float * __restrict  pphi,
	                                         const float * __restrict  pomg,
	                                         const zmm16c4_t phase,
	                                         const zmm16c4_t refi,
	                                         const zmm16c4_t px,
	                                         const zmm16c4_t py,
	                                         const zmm16c4_t pz,
	                                         zmm16c4_t & H_x,
	                                         zmm16c4_t & H_y,
	                                         zmm16c4_t & H_z,
	                                         zmm16c4_t & B_x,
	                                         zmm16c4_t & B_y,
	                                         zmm16c4_t & B_z) {
	                         
	               register __m512 tht = _mm512_loadu_ps(&ptht[0]);
	               register __m512 phi = _mm512_loadu_ps(&pphi[0]);
	               register __m512 omg = _mm512_loadu_ps(&pomg[0]);            
	               const __m512 c   = _mm512_set1_ps(299792458.0f); 
	               const __m512 mu0 = _mm512_set1_ps(0.0000012566370614359173f);   
	               const __m512 psi0 = _mm512_setzero_ps();
	               const __m512 C00  = _mm512_setzero_ps();
	               
	               zmm16c4_t H_x_1;
	               zmm16c4_t H_y_1;
	               zmm16c4_t H_z_1;
	               zmm16c4_t H_x_2;
	               zmm16c4_t H_y_2;
	               zmm16c4_t H_z_2;
	               zmm16c4_t k;
	               zmm16c4_t t0;
	               zmm16c4_t cdirx;
	               zmm16c4_t cdiry;
	               zmm16c4_t cdirz;
	               
	               register __m512 vpolx;
	               register __m512 vpoly;
	               register __m512 vpolz;
	               register __m512 vdirx;
	               register __m512 vdiry;
	               register __m512 vdirz;
	               register __m512 cn;
	               register __m512 x0;
	               register __m512 t0r,t0i;
	               register __m512 t1r,t1i;
	               register __m512 t2r,t2i;
	               
	               dir_vec_zmm16r4(tht,phi,&vdirx,
	                               &vdiry,&vdirz);
	               cdirx.re = vdirx;
	               x0       = _mm512_div_ps(omg,c);
	               cdirx.im = C00;
	               k.re     = _mm512_mul_ps(refi.re,x0);
	               k.im     = _mm512_mul_ps(refi.im,x0);
	               cdiry.re = vdiry;
	               cdiry.im = C00;
	               
	               pol_vec_zmm16r4(tht,phi,psi_0,
	                               &vpolx,&vpoly,&vpolz);
	               cdirz.re = vdirz;
	               cdirz.im = C00;
	               
	               H_XYZ_VP_zmm16c4(vpolx,vpoly,vpolz,
	                                vdirx,vdiry,vdirz,
	                                px,py,pz,
	                                H_x_1,H_y_1,H_z_1);
	                                
	               scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                               H_x_1,H_y_1,H_z_1,
	                               H_x_2,H_y_2,H_z_2);
	                               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_x_2.re,H_x_2.im,
	                            &t0r,&t0i);
	               H_x.re = _mm512_add_ps(H_x_1.re,t0r);
	               H_x.im = _mm512_add_ps(H_x_1.im,t0i);
	               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_y_2.re,H_y_2.im,
	                            &t1r,&t1i);
	               H_y.re = _mm512_add_ps(H_y_1.re,t1r);
	               H_y.im = _mm512_add_ps(H_y_1.im,t1i);
	               
	               cmul_zmm16r4(phase.re,phase.im,
	                            H_z_2.re,H_z_2.im,
	                            &t2r,&t2i);
	              H_z.re = _mm512_add_ps(H_z_2.re,t2r);
	              H_z.im = _mm512_add_ps(H_z_2.im,t2i);
	              
	              cn = cnorm_zmm16c4(H_x,H_y,H_z);
	              
	              x0     = _mm512_div_ps(omg,mu0);
	              H_x.re = _mm512_div_ps(H_x.re,cn);
	              H_x.im = _mm512_div_ps(H_x.im,cn);
	              H_y.re = _mm512_div_ps(H_y.re,cn);
	              H_y.im = _mm512_div_ps(H_y.im,cn);
	              H_z.re = _mm512_div_ps(H_z.re,cn);
	              H_z.im = _mm512_div_ps(H_z.im,cn); 
	              
	              t0.re  = _mm512_div_ps(k.re,x0);
	              t0.im  = _mm512_div_ps(k.im,x0);
	              
	              scrossc_zmm16c4(cdirx,cdiry,cdirz,
	                              H_x,H_y,H_z,
	                              H_x_2,H_y_2,H_z_2);
	                              
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_x_2.re,H_x_2.im,
	                           &B_x.re,&B_x.im);
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_y_2.re,H_y_2.im,
	                           &B_y.re,&B_y.im);
	              cmul_zmm16r4(t0.re,t0.im,
	                           H_z_2.re,H_z_2.im,
	                           &B_z.re,&B_z.im);
	              
	        }
	                                      
	       
	       
	       
	     
                
                
        } // radilocation

} // gms


#endif /*__GMS_EM_FIELDS_ZMM16R4_HPP__*/
