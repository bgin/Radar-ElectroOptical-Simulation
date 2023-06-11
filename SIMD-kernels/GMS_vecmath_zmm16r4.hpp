

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
	        
                
                
        } // math

} // gms


#endif /*__GMS_VECMATH_ZMM16R4_HPP__*/
