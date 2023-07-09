

#ifndef __GMS_ANTENNA_FEEDER_ZMM16R4_HPP__
#define __GMS_ANTENNA_FEEDER_ZMM16R4_HPP__ 090620230852


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

    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_ANTENNA_FEEDER_ZMM16R4_FULLVER =
      1000U*GMS_ANTENNA_FEEDER_ZMM16R4_MAJOR+
      100U*GMS_ANTENNA_FEEDER_ZMM16R4_MINOR+
      10U*GMS_ANTENNA_FEEDER_ZMM16R4_MICRO;
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_CREATION_DATE = "09-06-2023 08:52 PM +00200 (FRI 09 JUN 2023 GMT+2)";
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_ANTENNA_FEEDER_ZMM16R4_DESCRIPTION   = "AVX512 (single) optimized antenna and feeder kernels.";

}


#include <immintrin.h>
#include <complex>
#include <cstdint>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"
#include "GMS_cspint_quad.hpp"
#include "GMS_avint_quad.hpp"
#include "GMS_cubint_quad.hpp"
#include "GMS_filon_cos_quad.hpp"
#include "GMS_filon_sin_quad.hpp"
#include "GMS_hiordq_quad.hpp"
#include "GMS_plint_quad.hpp"
#include "GMS_wedint_quad.hpp"
#include "GMS_cephes.h"


namespace gms {

 
          namespace radiolocation {
          
          
                   struct __ATTR_ALIGN__(64) fwork_t {
                   
                          float * __restrict pxr;
                          float * __restrict pxi;
                          float * __restrict pyr;
                          float * __restrict pyi;
                          float * __restrict pzr;
                          float * __restrict pzi;
                          PAD_TO(0,16)
                   }
                   
                   
                  // Spherical unit vector.
                  struct __ATTR_ALIGN__(64) SUV_zmm16r4_t {
                  
                         __m512 x;
                         __m512 y;
                         __m512 z;
                  }
          
          
           
                   /*
                       Spherical unit vectors.
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void sphuv_zmm16r4(     const __m512 tht,
	                                   const __m512 phi,
	                                   SUV_zmm16r4_t & er,
	                                   SUV_zmm16r4_t & eth,
	                                   SUV_zmm16r4_t & eph) {
	               
	               using namespace gms::math;                   
	               register __m512 stht,sphi,cphi,ctht;
	               stht  = xsinf(tht);
	               ctht  = xcosf(tht);
	               cphi  = xcosf(phi);
	               sphi  = xsinf(phi);
	               er.x  = _mm512_mul_ps(stht,cphi);
	               er.y  = _mm512_mul_ps(stht,sphi);
	               er.z  = ctht; 
	               eth.x = _mm512_mul_ps(ctht,cphi);
	               eth.y = _mm512_mul_ps(ctht,sphi);
	               eth.z = negate_zmm16r4(stht);
	               eph.x = negate_zmm16r4(sphi);
	               eph.y = cphi;
	               eph.z = _mm512_setzero_ps();     
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void sphuv_zmm16r4_a(   const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                   const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                   SUV_zmm16r4_t & er,
	                                   SUV_zmm16r4_t & eth,
	                                   SUV_zmm16r4_t & eph) {
	               
	               using namespace gms::math;                   
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 tht,phi;
	               tht   = _mm512_load_ps(&ptht[0]);
	               phi   = _mm512_load_ps(&pphi[0]);
	               stht  = xsinf(tht);
	               ctht  = xcosf(tht);
	               cphi  = xcosf(phi);
	               sphi  = xsinf(phi);
	               er.x  = _mm512_mul_ps(stht,cphi);
	               er.y  = _mm512_mul_ps(stht,sphi);
	               er.z  = ctht; 
	               eth.x = _mm512_mul_ps(ctht,cphi);
	               eth.y = _mm512_mul_ps(ctht,sphi);
	               eth.z = negate_zmm16r4(stht);
	               eph.x = negate_zmm16r4(sphi);
	               eph.y = cphi;
	               eph.z = _mm512_setzero_ps();     
	       } 
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void sphuv_zmm16r4_u(   const float * __restrict  ptht,
	                                   const float * __restrict  pphi,
	                                   SUV_zmm16r4_t & er,
	                                   SUV_zmm16r4_t & eth,
	                                   SUV_zmm16r4_t & eph) {
	               
	               using namespace gms::math;                   
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 tht,phi;
	               tht   = _mm512_loadu_ps(&ptht[0]);
	               phi   = _mm512_loadu_ps(&pphi[0]);
	               stht  = xsinf(tht);
	               ctht  = xcosf(tht);
	               cphi  = xcosf(phi);
	               sphi  = xsinf(phi);
	               er.x  = _mm512_mul_ps(stht,cphi);
	               er.y  = _mm512_mul_ps(stht,sphi);
	               er.z  = ctht; 
	               eth.x = _mm512_mul_ps(ctht,cphi);
	               eth.y = _mm512_mul_ps(ctht,sphi);
	               eth.z = negate_zmm16r4(stht);
	               eph.x = negate_zmm16r4(sphi);
	               eph.y = cphi;
	               eph.z = _mm512_setzero_ps();     
	       } 
	       
	       
	         
	       
	       /*
	             Function 'N' = Nth(theta,phi)*eth+Nphi(theta,phi)*ephi
	             Formula 1-3, p. 10
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void N_f13_zmm16r4(const __m512 ntr,
	                              const __m512 nti,
	                              const __m512 npr,
	                              const __m512 npi,
	                              const SUV_zmm16r4_t eth,
	                              const SUV_zmm16r4_t eph,
	                              __m512 * __restrict Nthr,
	                              __m512 * __restrict Nthi,
	                              __m512 * __restrict Nphr,
	                              __m512 * __restrict Nphi) {
	                              
	                register __m512 t0r,t0i,t1r,t1i;
	                t0r   = _mm512_fmadd_ps(ntr,eth.x,
	                                  _mm512_fmsub_ps(ntr,eth.y,
	                                              _mm512_mul_ps(ntr,eth.z)));
	                *Nthr = t0r;
	                t0i   = _mm512_fmadd_ps(nti,eth.x,
	                                  _mm512_fmsub_ps(nti,eth.y,
	                                              _mm512_mul_ps(nti,eth.z)));   
	                *Nthi = t0i;
	                t1r   = _mm512_fmadd_ps(npr,eph.x,
	                                    _mm512_mul_ps(npr,eph.y));
	                *Nphr = t1r;
	                t1i   = _mm512_fmadd_ps(npi,eph.x,
	                                    _mm512_mul_ps(npi,eph.y));   
	                *Nphi = t1i;              
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void N_f13_zmm16r4_a(const float *  __restrict __ATTR_ALIGN__(64) pntr,
	                                const float *  __restrict __ATTR_ALIGN__(64) pnti,
	                                const float *  __restrict __ATTR_ALIGN__(64) pnpr,
	                                const float *  __restrict __ATTR_ALIGN__(64) pnpi,
	                                const SUV_zmm16r4_t eth,
	                                const SUV_zmm16r4_t eph,
	                                float * __restrict __ATTR_ALIGN__(64) pNthr,
	                                float * __restrict __ATTR_ALIGN__(64) pNthi,
	                                float * __restrict __ATTR_ALIGN__(64) pNphr,
	                                float * __restrict __ATTR_ALIGN__(64) pNphi) {
	                              
	                register __m512 t0r,t0i,t1r,t1i;
	                register __m512 ntr,nti,npr,npi;
	                ntr   = _mm512_load_ps(&pntr[0]);
	                nti   = _mm512_load_ps(&pnti[0]);
	                npr   = _mm512_load_ps(&pnpr[0]);
	                npi   = _mm512_load_ps(&pnpi[0]);
	                t0r   = _mm512_fmadd_ps(ntr,eth.x,
	                                  _mm512_fmsub_ps(ntr,eth.y,
	                                              _mm512_mul_ps(ntr,eth.z)));
	                _mm512_store_ps(&Nthr[0] ,t0r);
	                t0i   = _mm512_fmadd_ps(nti,eth.x,
	                                  _mm512_fmsub_ps(nti,eth.y,
	                                              _mm512_mul_ps(nti,eth.z)));   
	                _mm512_store_ps(&Nthi[0] ,t0i);
	                t1r   = _mm512_fmadd_ps(npr,eph.x,
	                                    _mm512_mul_ps(npr,eph.y));
	                _mm512_store_ps(&Nphr[0] ,t1r);
	                t1i   = _mm512_fmadd_ps(npi,eph.x,
	                                    _mm512_mul_ps(npi,eph.y));   
	                _mm512_store_ps(&Nphi[0] ,t1i);              
	       }
	                              
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void N_f13_zmm16r4_u(const float *  __restrict  pntr,
	                                const float *  __restrict  pnti,
	                                const float *  __restrict  pnpr,
	                                const float *  __restrict  pnpi,
	                                const SUV_zmm16r4_t eth,
	                                const SUV_zmm16r4_t eph,
	                                float * __restrict  pNthr,
	                                float * __restrict  pNthi,
	                                float * __restrict  pNphr,
	                                float * __restrict  pNphi) {
	                              
	                register __m512 t0r,t0i,t1r,t1i;
	                register __m512 ntr,nti,npr,npi;
	                ntr   = _mm512_loadu_ps(&pntr[0]);
	                nti   = _mm512_loadu_ps(&pnti[0]);
	                npr   = _mm512_loadu_ps(&pnpr[0]);
	                npi   = _mm512_loadu_ps(&pnpi[0]);
	                t0r   = _mm512_fmadd_ps(ntr,eth.x,
	                                  _mm512_fmsub_ps(ntr,eth.y,
	                                              _mm512_mul_ps(ntr,eth.z)));
	                _mm512_storeu_ps(&Nthr[0] ,t0r);
	                t0i   = _mm512_fmadd_ps(nti,eth.x,
	                                  _mm512_fmsub_ps(nti,eth.y,
	                                              _mm512_mul_ps(nti,eth.z)));   
	                _mm512_storeu_ps(&Nthi[0] ,t0i);
	                t1r   = _mm512_fmadd_ps(npr,eph.x,
	                                    _mm512_mul_ps(npr,eph.y));
	                _mm512_storeu_ps(&Nphr[0] ,t1r);
	                t1i   = _mm512_fmadd_ps(npi,eph.x,
	                                    _mm512_mul_ps(npi,eph.y));   
	                _mm512_storeu_ps(&Nphi[0] ,t1i);              
	       }
	       
	       
	       /*
	           Formula 1-11, p. 14
	       */
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f111_zmm16r4(const __m512 Nthr,
	                               const __m512 Nthi,
	                               const __m512 Nphr,
	                               const __m512 Nphi,
	                               const float nmaxr,
	                               const float nmaxi) {
	                  
	                  register __m512 t0,t1,st0,st1;
	                  register __m512 x0,x1,P;
	                  x0 = _mm512_set1_ps(nmaxr);
	                  t0 = cabs_zmm16r4(Nthr,Nthi);
	                  x1 = _mm512_set1_ps(nmaxi);
	                  t1 = cabs_zmm16r4(Nphr,Nphi);
	                  st0= _mm512_fmadd_ps(t0,t0,
	                                   _mm512_mul_ps(t1,t1));
	                  st1= _mm512_add_ps(x0,x1);
	                  P  = _mm512_div_ps(st0,st1);
	                  return (P);
	               }  
	               
	               
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f111_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pNthr,
	                                 const float * __restrict __ATTR_ALIGN__(64) pNthi,
	                                 const float * __restrict __ATTR_ALIGN__(64) pNphr,
	                                 const float * __restrict __ATTR_ALIGN__(64) pNphi,
	                                 const float nmaxr,
	                                 const float nmaxi) {
	                  
	                  register __m512 Nthr = _mm512_load_ps(&pNthr[0]);
	                  register __m512 Nthi = _mm512_load_ps(&pNthi[0]);
	                  register __m512 Nphr = _mm512_load_ps(&pNphr[0]);
	                  register __m512 Nphi = _mm512_load_ps(&pNphi[0]);
	                  register __m512 t0,t1,st0,st1;
	                  register __m512 x0,x1,P;
	                  x0 = _mm512_set1_ps(nmaxr);
	                  t0 = cabs_zmm16r4(Nthr,Nthi);
	                  x1 = _mm512_set1_ps(nmaxi);
	                  t1 = cabs_zmm16r4(Nphr,Nphi);
	                  st0= _mm512_fmadd_ps(t0,t0,
	                                   _mm512_mul_ps(t1,t1));
	                  st1= _mm512_add_ps(x0,x1);
	                  P  = _mm512_div_ps(st0,st1);
	                  return (P);
	               }   
	               
	               
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f111_zmm16r4_u(const float * __restrict  pNthr,
	                                 const float * __restrict  pNthi,
	                                 const float * __restrict  pNphr,
	                                 const float * __restrict  pNphi,
	                                 const float nmaxr,
	                                 const float nmaxi) {
	                  
	                  register __m512 Nthr = _mm512_loadu_ps(&pNthr[0]);
	                  register __m512 Nthi = _mm512_loadu_ps(&pNthi[0]);
	                  register __m512 Nphr = _mm512_loadu_ps(&pNphr[0]);
	                  register __m512 Nphi = _mm512_loadu_ps(&pNphi[0]);
	                  register __m512 t0,t1,st0,st1;
	                  register __m512 x0,x1,P;
	                  x0 = _mm512_set1_ps(nmaxr);
	                  t0 = cabs_zmm16r4(Nthr,Nthi);
	                  x1 = _mm512_set1_ps(nmaxi);
	                  t1 = cabs_zmm16r4(Nphr,Nphi);
	                  st0= _mm512_fmadd_ps(t0,t0,
	                                   _mm512_mul_ps(t1,t1));
	                  st1= _mm512_add_ps(x0,x1);
	                  P  = _mm512_div_ps(st0,st1);
	                  return (P);
	               }    
	               
	               
	               /*
	                  Formula 1-14, p. 18
	                  
	               */    
	               
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline       
	           __m512 f113_zmm16r4(const __m512 tht,
	                               const __m512 phi,
	                               const float  a,
	                               const float  b,
	                               const float  g,
	                               const float  d) {
	                               
	                  register __m512 va,vb,vg,vd;
	                  register __m512 stht,cphi,sphi,ctht;
	                  register __m512 N,t0,t1;
	                  va   = _mm512_set1_ps(a);
	                  stht = xsinf(tht);
	                  t0   = _mm512_mul_ps(va,stht);
	                  vb   = _mm512_set1_ps(b);
	                  cphi = xcosf(phi);
	                  t1   = _mm512_mul_ps(vb,stht);
	                  vg   = _mm512_set1_ps(g);
	                  sphi = xsinf(phi);
	                  vd   = _mm512_set1_ps(d);
	                  ctht = xcosf(tht);
	                  N    = _mm512_fmadd_ps(t0,cphi,
	                                    _mm512_fmadd_ps(t1,sphi,
	                                                _mm512_fmadd_ps(vg,ctht,vd)));
	                  return (N);                                                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline       
	           __m512 f113_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                 const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                 const float * __restrict __ATTR_ALIGN__(64) pa,
	                               const float  b,
	                               const float  g,
	                               const float  d) {
	                  
	                  register __m512 tht = _mm512_load_ps(&ptht[0]);
	                  register __m512 phi = _mm512_load_ps(&pphi[0]);
	                  register __m512 a   = _mm512_load_ps(&pa[0]);             
	                  register __m512 va,vb,vg,vd;
	                  register __m512 stht,cphi,sphi,ctht;
	                  register __m512 N,t0,t1;
	                  va   = _mm512_set1_ps(a);
	                  stht = xsinf(tht);
	                  t0   = _mm512_mul_ps(va,stht);
	                  vb   = _mm512_set1_ps(b);
	                  cphi = xcosf(phi);
	                  t1   = _mm512_mul_ps(vb,stht);
	                  vg   = _mm512_set1_ps(g);
	                  sphi = xsinf(phi);
	                  vd   = _mm512_set1_ps(d);
	                  ctht = xcosf(tht);
	                  N    = _mm512_fmadd_ps(t0,cphi,
	                                    _mm512_fmadd_ps(t1,sphi,
	                                                _mm512_fmadd_ps(vg,ctht,vd)));
	                  return (N);                                                    
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline       
	           __m512 f113_zmm16r4_u(const float * __restrict  ptht,
	                                 const float * __restrict  pphi,
	                                 const float * __restrict  pa,
	                                 const float  b,
	                                 const float  g,
	                                 const float  d) {
	                  
	                  register __m512 tht = _mm512_loadu_ps(&ptht[0]);
	                  register __m512 phi = _mm512_loadu_ps(&pphi[0]);
	                  register __m512 a   = _mm512_loadu_ps(&pa[0]);             
	                  register __m512 va,vb,vg,vd;
	                  register __m512 stht,cphi,sphi,ctht;
	                  register __m512 N,t0,t1;
	                  va   = _mm512_set1_ps(a);
	                  stht = xsinf(tht);
	                  t0   = _mm512_mul_ps(va,stht);
	                  vb   = _mm512_set1_ps(b);
	                  cphi = xcosf(phi);
	                  t1   = _mm512_mul_ps(vb,stht);
	                  vg   = _mm512_set1_ps(g);
	                  sphi = xsinf(phi);
	                  vd   = _mm512_set1_ps(d);
	                  ctht = xcosf(tht);
	                  N    = _mm512_fmadd_ps(t0,cphi,
	                                    _mm512_fmadd_ps(t1,sphi,
	                                                _mm512_fmadd_ps(vg,ctht,vd)));
	                  return (N);                                                    
	        }
	        
	        
	        /*
	           Formula 1-14, p. 18
	        */
	        
	           __ATTR_ALIGN__(32)
                   __ATTR_HOT__
	           static inline  
	           void f114_r4(const float gam,
	                        const float a,
	                        const float b,
	                        const float g,
	                        float * __restrict x0,
	                        float * __restrict y0,
	                        float * __restrict z0) {
	           
	                constexpr float C6283185307179586476925286766559 =
	                                     6.283185307179586476925286766559f;
	                const float rat = gam/C6283185307179586476925286766559;
	                *x0 = rat*a;
	                *y0 = rat*b;
	                *z0 = rat*g;         
	         }
	         
	         
	         /*
	            Formula 1-16, p. 19
	         */
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f116_zmm16r4(const __m512 Nthr,
	                               const __m512 Nthi,
	                               const __m512 Nphr,
	                               const __m512 Nphi,
	                               const float nmaxr,
	                               const float nmaxi,
	                               const float maxp) {
	                               
	                  register __m512 P,D,t0;
	                  t0 = _mm512_set1_ps(maxp);
	                  P  = f114_zmm16r4(Nthr,Nthi,Nphr,Nphi,nmaxr,nmaxi);
	                  D  = _mm512_div_pd(P,t0);
	                  return (D);                      
	        }
	        
	        
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f116_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) pNthr,
	                                 const float * __restrict __ATTR_ALIGN__(64) pNthi,
	                                 const float * __restrict __ATTR_ALIGN__(64) pNphr,
	                                 const float * __restrict __ATTR_ALIGN__(64) pNphi,
	                                 const float nmaxr,
	                                 const float nmaxi
	                                 const float maxp) {
	                         
	                        
	                  register __m512 P,D,t0;
	                  t0 = _mm512_set1_ps(maxp);
	                  P  = f114_zmm16r4_a(Nthr,Nthi,Nphr,Nphi,nmaxr,nmaxi);
	                  D  = _mm512_div_pd(P,t0);
	                  return (D);                      
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f116_zmm16r4_u(  const float * __restrict  pNthr,
	                                 const float * __restrict  pNthi,
	                                 const float * __restrict  pNphr,
	                                 const float * __restrict  pNphi,
	                                 const float nmaxr,
	                                 const float nmaxi
	                                 const float maxp) {
	                         
	                        
	                  register __m512 P,D,t0;
	                  t0 = _mm512_set1_ps(maxp);
	                  P  = f114_zmm16r4_u(Nthr,Nthi,Nphr,Nphi,nmaxr,nmaxi);
	                  D  = _mm512_div_pd(P,t0);
	                  return (D);                      
	        }
	        
	        
	        /*
	              Formula 1-9, p. 13
	        */
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 fth_f19_zmm16r4(const __m512 nthr,
	                               const __m512 nthi,
	                               const float nmax) {
	                 
	                 register __m512 vmax,abs,fth;
	                 vmax = _mm512_set1_ps(nmax);
	                 abs  = cabs_zmm16r4(nthr,nthi);
	                 fth  = _mm512_div_ps(abs,vmax);
	                 return (fth);                       
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 fth_f19_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pnthr,
	                                    const float * __restrict __ATTR_ALIGN__(64) pnthi,
	                                    const float nmax) {
	                 
	                 register __m512 nthr = _mm512_load_ps(&pnthr[0]);
	                 register __m512 nthi = _mm512_load_ps(&pnthi[0]);
	                 register __m512 vmax,abs,fth;
	                 vmax = _mm512_set1_ps(nmax);
	                 abs  = cabs_zmm16r4(nthr,nthi);
	                 fth  = _mm512_div_ps(abs,vmax);
	                 return (fth);                       
	         }
	         
	         
	          
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 fth_f19_zmm16r4_u(const float * __restrict  pnthr,
	                                    const float * __restrict  pnthi,
	                                    const float nmax) {
	                 
	                 register __m512 nthr = _mm512_loadu_ps(&pnthr[0]);
	                 register __m512 nthi = _mm512_loadu_ps(&pnthi[0]);
	                 register __m512 vmax,abs,fth;
	                 vmax = _mm512_set1_ps(nmax);
	                 abs  = cabs_zmm16r4(nthr,nthi);
	                 fth  = _mm512_div_ps(abs,vmax);
	                 return (fth);                       
	         }
	         
	         
	        /*
	              Formula 1-9, p. 13
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 fph_f19_zmm16r4(const __m512 nphr,
	                                  const __m512 nphi,
	                                  const float nmax) {
	                 
	                 register __m512 vmax,abs,fph;
	                 vmax = _mm512_set1_ps(nmax);
	                 abs  = cabs_zmm16r4(nphr,nphi);
	                 fth  = _mm512_div_ps(abs,vmax);
	                 return (fph);                       
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 fph_f19_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pnphr,
	                                    const float * __restrict __ATTR_ALIGN__(64) pnphi,
	                                    const float nmax) {
	                 
	                 register __m512 nphr = _mm512_load_ps(&pnphr[0]);
	                 register __m512 nphi = _mm512_load_ps(&pnphi[0]);
	                 register __m512 vmax,abs,fph;
	                 vmax = _mm512_set1_ps(nmax);
	                 abs  = cabs_zmm16r4(nphr,nphi);
	                 fth  = _mm512_div_ps(abs,vmax);
	                 return (fph);                       
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 fph_f19_zmm16r4_u(const float * __restrict  pnphr,
	                                    const float * __restrict  pnphi,
	                                    const float nmax) {
	                 
	                 register __m512 nphr = _mm512_loadu_ps(&pnphr[0]);
	                 register __m512 nphi = _mm512_loadu_ps(&pnphi[0]);
	                 register __m512 vmax,abs,fph;
	                 vmax = _mm512_set1_ps(nmax);
	                 abs  = cabs_zmm16r4(nphr,nphi);
	                 fth  = _mm512_div_ps(abs,vmax);
	                 return (fph);                       
	         }
	         
	         
	       /*
	            Hertz vector (electrical,magnetic), avint integrator.
	            Formula 2-13, 2-15, p. 35
	       */
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvme_f2135_zmm16r4_avint(const __m512 xre,
	                                       const __m512 xim,
	                                       const __m512 yre,
	                                       const __m512 yim,
	                                       const __m512 zre,
	                                       const __m512_zim,
	                                       const __m512 xd,
	                                       const __m512 yd,
	                                       const __m512 zd,
	                                       const float arg[10],
	                                       std::complex<float> & hx,                        
                                               std::complex<float> & hy,
                                               std::complex<float> & hz,
                                               int32_t & ierr) {
                            
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi                
                      
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                       
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        _mm512_store_ps(&dx[0],xd);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        _mm512_store_ps(&dy[0],yd);
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        _mm512_store_ps(&dz[0],zd);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&dx[0],&pxr[0],xa,xb,ier1);
                        sxi = avint(&dx[0],&pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&dy[0],&pyr[0],ya,yb,ier3);
                        syi = avint(&dy[0],&pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&dz[0],&pzr[0],za,zb,ier5);
                        szi = avint(&dz[0],&pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {
                           ierr = 3;
                           return;
                        }         
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_avint_a(const float * __restrict __ATTR_ALIGN__(64) pxre,
	                                         const float * __restrict __ATTR_ALIGN__(64) pxim,
	                                         const float * __restrict __ATTR_ALIGN__(64) pyre,
	                                         const float * __restrict __ATTR_ALIGN__(64) pyim,
	                                         const float * __restrict __ATTR_ALIGN__(64) pzre,
	                                         const float * __restrict __ATTR_ALIGN__(64) pzim,
	                                         float * __restrict __ATTR_ALIGN__(64) pxd,
	                                         float * __restrict __ATTR_ALIGN__(64) pyd,
	                                         float * __restrict __ATTR_ALIGN__(64) pzd,
	                                         const float arg[10],
	                                         std::complex<float> & hx,                        
                                                 std::complex<float> & hy,
                                                 std::complex<float> & hz,
                                                 int32_t & ierr) {
                            
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi                
                     
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xre,xim,yre,yim,zre,zim;
                        register __m512 xd,yd,zd;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                        xre = _mm512_load_ps(&pxre[0]);
                        xim = _mm512_load_ps(&pxim[0]);
                        yre = _mm512_load_ps(&pyre[0]);
                        yim = _mm512_load_ps(&pyim[0]);
                        zre = _mm512_load_ps(&pzre[0]);
                        zim = _mm512_load_ps(&pzim[0]);
                        
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&pxd[0],&pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&pyd[0],&pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&pzd[0],&pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {
                           ierr = 3;
                           return;
                        }         
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
#if !defined(__ANTENNA_FEEDER_PF_CACHE_HINT__)
#define  __ANTENNA_FEEDER_PF_CACHE_HINT__ 1
#endif      

                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void  f2135_integrand_zmm16r4_u6x_a(const float * __restrict __ATTR_ALIGN__(64) pxre,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pxim,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pyre,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pyim,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pzre,
	                                              const float * __restrict  __ATTR_ALIGN__(64) pzim,
	                                              fwork_t fw, //work arrays (caller allocated)
	                                              const __m512 cer,
	                                              const __m512 cei,
	                                              const std::complex<float> tmp,
                                                      const int32_t n,
                                                      const int32_t  PF_DIST) {
                                                      
                           if(__builtin_expect(n<=0,0)) { return;}
                           register __m512 xre,xim,yre,yim,zre,zim;
                            __m512 t0r,t0i,t1r,t1i,t2r,t2i;
                           int32_t i;
                           
                           for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_NTA);
#endif
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_load_ps(&pxre[i+16]);
                             xim = _mm512_load_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_load_ps(&pyre[i+16]);
                             yim = _mm512_load_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_load_ps(&pzre[i+16]);
                             zim = _mm512_load_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+16],t2r);
                             _mm512_store_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_load_ps(&pxre[i+32]);
                             xim = _mm512_load_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+32],t0r);
                             _mm512_store_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_load_ps(&pyre[i+32]);
                             yim = _mm512_load_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+32],t1r);
                             _mm512_store_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_load_ps(&pzre[i+32]);
                             zim = _mm512_load_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+32],t2r);
                             _mm512_store_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_load_ps(&pxre[i+48]);
                             xim = _mm512_load_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+48],t0r);
                             _mm512_store_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_load_ps(&pyre[i+48]);
                             yim = _mm512_load_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+48],t1r);
                             _mm512_store_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_load_ps(&pzre[i+48]);
                             zim = _mm512_load_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+48],t2r);
                             _mm512_store_ps(&fw.zyi[i+48],t2i)
                             xre = _mm512_load_ps(&pxre[i+64]);
                             xim = _mm512_load_ps(&pxim[i+64]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+64],t0r);
                             _mm512_store_ps(&fw.pxi[i+64],t0i);
                             yre = _mm512_load_ps(&pyre[i+64]);
                             yim = _mm512_load_ps(&pyim[i+64]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+64],t1r);
                             _mm512_store_ps(&fw.pyi[i+64],t1i);
                             zre = _mm512_load_ps(&pzre[i+64]);
                             zim = _mm512_load_ps(&pzim[i+64]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+64],t2r);
                             _mm512_store_ps(&fw.zyi[i+64],t2i)
                             xre = _mm512_load_ps(&pxre[i+80]);
                             xim = _mm512_load_ps(&pxim[i+80]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+80],t0r);
                             _mm512_store_ps(&fw.pxi[i+80],t0i);
                             yre = _mm512_load_ps(&pyre[i+80]);
                             yim = _mm512_load_ps(&pyim[i+80]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+80],t1r);
                             _mm512_store_ps(&fw.pyi[i+80],t1i);
                             zre = _mm512_load_ps(&pzre[i+80]);
                             zim = _mm512_load_ps(&pzim[i+80]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+80],t2r);
                             _mm512_store_ps(&fw.zyi[i+80],t2i)
                        }  
                        
                        for(; (i+63) < n; i += 64) {
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_load_ps(&pxre[i+16]);
                             xim = _mm512_load_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_load_ps(&pyre[i+16]);
                             yim = _mm512_load_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_load_ps(&pzre[i+16]);
                             zim = _mm512_load_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+16],t2r);
                             _mm512_store_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_load_ps(&pxre[i+32]);
                             xim = _mm512_load_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+32],t0r);
                             _mm512_store_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_load_ps(&pyre[i+32]);
                             yim = _mm512_load_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+32],t1r);
                             _mm512_store_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_load_ps(&pzre[i+32]);
                             zim = _mm512_load_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+32],t2r);
                             _mm512_store_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_load_ps(&pxre[i+48]);
                             xim = _mm512_load_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+48],t0r);
                             _mm512_store_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_load_ps(&pyre[i+48]);
                             yim = _mm512_load_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+48],t1r);
                             _mm512_store_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_load_ps(&pzre[i+48]);
                             zim = _mm512_load_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+48],t2r);
                             _mm512_store_ps(&fw.zyi[i+48],t2i)
                      }   
                      
                      for(; (i+31) < n; i += 32) {
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_load_ps(&pxre[i+16]);
                             xim = _mm512_load_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+16],t0r);
                             _mm512_store_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_load_ps(&pyre[i+16]);
                             yim = _mm512_load_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+16],t1r);
                             _mm512_store_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_load_ps(&pzre[i+16]);
                             zim = _mm512_load_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+16],t2r);
                             _mm512_store_ps(&fw.zyi[i+16],t2i) 
                      }   
                      
                      for(; (i+15) < n; i += 16) {
                             xre = _mm512_load_ps(&pxre[i+0]);
                             xim = _mm512_load_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_store_ps(&fw.pxr[i+0],t0r);
                             _mm512_store_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_load_ps(&pyre[i+0]);
                             yim = _mm512_load_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_store_ps(&fw.pyr[i+0],t1r);
                             _mm512_store_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_load_ps(&pzre[i+0]);
                             zim = _mm512_load_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_store_ps(&fw.zyr[i+0],t2r);
                             _mm512_store_ps(&fw.zyi[i+0],t2i);
                      }  
                      
                     
                      for(; (i+0) < n; i += 1) {
                           register float xr               = pxre[i];
                           register float xi               = pxim[i];
                           register std::complex<float> cx = {xr,xi};
                           register std::complex<float> xx = cx*tmp2;
                           fw.pxr[i]                       = xx.real(); 
                           fw.pxi[i]                       = xx.imag();
                           register float yr               = pyre[i];
                           register float yi               = pyim[i];
                           register std::complex<float> cy = {yr,yi};
                           register std::complex<float> yy = cy*tmp2;
                           fw.pyr[i]                       = yy.real(); 
                           fw.pyi[i]                       = yy.imag();
                           register float zr               = pzre[i];
                           register float zi               = pzim[i];
                           register std::complex<float> cz = {zr,zi};
                           register std::complex<float> zz = cz*tmp2;
                           fw.pzr[i]                       = zz.real(); 
                           fw.pzi[i]                       = zz.imag();
                      }  
                      
                                                       
                  }
                  
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void  f2135_integrand_zmm16r4_u6x_u(const float * __restrict  pxre,
	                                              const float * __restrict  pxim,
	                                              const float * __restrict  pyre,
	                                              const float * __restrict  pyim,
	                                              const float * __restrict  pzre,
	                                              const float * __restrict  pzim,
	                                              fwork_t fw, //work arrays (caller allocated)
	                                              const __m512 cer,
	                                              const __m512 cei,
	                                              const std::complex<float> tmp,
                                                      const int32_t n,
                                                      const int32_t  PF_DIST) {
                                                      
                                                      
                           if(__builtin_expect(n<=0,0)) { return;}
                           register __m512 xre,xim,yre,yim,zre,zim;
                            __m512 t0r,t0i,t1r,t1i,t2r,t2i;
                           int32_t i;
                           
                           for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pxre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pxim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pyim[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzre[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pzim[i+PF_DIST],_MM_HINT_NTA);
#endif
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+16]);
                             xim = _mm512_loadu_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+16]);
                             yim = _mm512_loadu_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+16]);
                             zim = _mm512_loadu_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+16],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+32]);
                             xim = _mm512_loadu_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+32]);
                             yim = _mm512_loadu_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+32]);
                             zim = _mm512_loadu_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+32],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+48]);
                             xim = _mm512_loadu_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+48]);
                             yim = _mm512_loadu_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+48]);
                             zim = _mm512_loadu_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+48],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+48],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+64]);
                             xim = _mm512_loadu_ps(&pxim[i+64]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+64],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+64],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+64]);
                             yim = _mm512_loadu_ps(&pyim[i+64]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+64],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+64],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+64]);
                             zim = _mm512_loadu_ps(&pzim[i+64]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+64],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+64],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+80]);
                             xim = _mm512_loadu_ps(&pxim[i+80]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+80],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+80],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+80]);
                             yim = _mm512_loadu_ps(&pyim[i+80]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+80],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+80],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+80]);
                             zim = _mm512_loadu_ps(&pzim[i+80]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+80],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+80],t2i)
                        }  
                        
                        for(; (i+63) < n; i += 64) {
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+16]);
                             xim = _mm512_loadu_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+16]);
                             yim = _mm512_loadu_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+16]);
                             zim = _mm512_loadu_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+16],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+16],t2i)
                             xre = _mm512_loadu_ps(&pxre[i+32]);
                             xim = _mm512_loadu_ps(&pxim[i+32]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+32]);
                             yim = _mm512_loadu_ps(&pyim[i+32]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+32]);
                             zim = _mm512_loadu_ps(&pzim[i+32]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+32],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+32],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+48]);
                             xim = _mm512_loadu_ps(&pxim[i+48]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+48]);
                             yim = _mm512_loadu_ps(&pyim[i+48]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+48]);
                             zim = _mm512_loadu_ps(&pzim[i+48]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+48],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+48],t2i)
                      }   
                      
                      for(; (i+31) < n; i += 32) {
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);
                             xre = _mm512_loadu_ps(&pxre[i+16]);
                             xim = _mm512_loadu_ps(&pxim[i+16]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+16]);
                             yim = _mm512_loadu_ps(&pyim[i+16]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+16]);
                             zim = _mm512_loadu_ps(&pzim[i+16]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+16],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+16],t2i) 
                      }   
                      
                      for(; (i+15) < n; i += 16) {
                             xre = _mm512_loadu_ps(&pxre[i+0]);
                             xim = _mm512_loadu_ps(&pxim[i+0]);
                             cmul_zmm16r4(xre,xim,cer,cei,&t0r,&t0i);
                             _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                             _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                             yre = _mm512_loadu_ps(&pyre[i+0]);
                             yim = _mm512_loadu_ps(&pyim[i+0]);
                             cmul_zmm16r4(yre,yim,cer,cei,&t1r,&t1i);
                             _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                             _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                             zre = _mm512_loadu_ps(&pzre[i+0]);
                             zim = _mm512_loadu_ps(&pzim[i+0]);
                             cmul_zmm16r4(zre,zim,cer,cei,&t2r,&t2i);
                             _mm512_storeu_ps(&fw.zyr[i+0],t2r);
                             _mm512_storeu_ps(&fw.zyi[i+0],t2i);  
                      }  
                      
                     
                      for(; (i+0) < n; i += 1) {
                           register float xr               = pxre[i];
                           register float xi               = pxim[i];
                           register std::complex<float> cx = {xr,xi};
                           register std::complex<float> xx = cx*tmp2;
                           fw.pxr[i]                       = xx.real(); 
                           fw.pxi[i]                       = xx.imag();
                           register float yr               = pyre[i];
                           register float yi               = pyim[i];
                           register std::complex<float> cy = {yr,yi};
                           register std::complex<float> yy = cy*tmp2;
                           fw.pyr[i]                       = yy.real(); 
                           fw.pyi[i]                       = yy.imag();
                           register float zr               = pzre[i];
                           register float zi               = pzim[i];
                           register std::complex<float> cz = {zr,zi};
                           register std::complex<float> zz = cz*tmp2;
                           fw.pzr[i]                       = zz.real(); 
                           fw.pzi[i]                       = zz.imag();
                      }  
                      
                                                       
                  }
                  
                  
                  
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem5_f213_zmm16r4_avint_u6x(const float * __restrict pxre,
	                                             const float * __restrict pxim,
	                                             const float * __restrict pyre,
	                                             const float * __restrict pyim,
	                                             const float * __restrict pzre,
	                                             const float * __restrict pzim,
	                                             float * __restrict  pxd,
	                                             float * __restrict  pyd,
	                                             float * __restrict  pzd,
	                                             fwork_t fw, //work arrays (caller allocated)
	                                             const float arg[10],
	                                             std::complex<float> & hx,                        
                                                     std::complex<float> & hy,
                                                     std::complex<float> & hz,
                                                     const int32_t n,
                                                     const int32_t  PF_DIST,
                                                     int32_t & ierr,
                                                     const bool aligned) {
                                                 
                       
                        
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi  
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        //register __m512 xre,xim,yre,yim,zre,zim;
                        register __m512 xd,yd,zd;
                        std::complex<float> tmp1,tmp2;
                        register float k,r,xa,xb,ya,yb,za,zb,inv;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6; 
                        //int32_t i;
                        k = arg[0];
                        r = arg[1];
                        inv  = 1.0f/r;
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        tmp1 = {0.0f,-1.0f*k*r};
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        tmp2 = tmp1*inv;
                        eps  = arg[9];   
                        if(aligned) {
                           f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                         fw,cer,cei,tmp2,n,PF_DIST);
                        }
                        else {
                           f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                         fw,cer,cei,tmp2,n,PF_DIST); 
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&pxd[0],&fw.pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&fw.pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&pyd[0],&fw.pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&fw.pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&pzd[0],&fw.pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&fw.pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {
                           ierr = 3;
                           return;
                        }         
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                     
               }
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_avint_u(const float * __restrict  pxre,
	                                         const float * __restrict  pxim,
	                                         const float * __restrict  pyre,
	                                         const float * __restrict  pyim,
	                                         const float * __restrict  pzre,
	                                         const float * __restrict  pzim,
	                                         float * __restrict  pxd,
	                                         float * __restrict  pyd,
	                                         float * __restrict  pzd,
	                                         const float arg[10],
	                                         std::complex<float> & hx,                        
                                                 std::complex<float> & hy,
                                                 std::complex<float> & hz,
                                                 int32_t & ierr) {
                            
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi                
                     
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xre,xim,yre,yim,zre,zim;
                        register __m512 xd,yd,zd;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                        xre = _mm512_loadu_ps(&pxre[0]);
                        xim = _mm512_loadu_ps(&pxim[0]);
                        yre = _mm512_loadu_ps(&pyre[0]);
                        yim = _mm512_loadu_ps(&pyim[0]);
                        zre = _mm512_loadu_ps(&pzre[0]);
                        zim = _mm512_loadu_ps(&pzim[0]);
                        
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&pxd[0],&pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&pyd[0],&pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&pzd[0],&pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {
                           ierr = 3;
                           return;
                        }         
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
                  
	       /*
	            Hertz vector (electrical), cubint integrator.
	            Formula 2-13, p. 35
	       */
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_cubint(const __m512 xre,
	                                       const __m512 xim,
	                                       const __m512 yre,
	                                       const __m512 yim,
	                                       const __m512 zre,
	                                       const __m512 zim,
	                                       const __m512 xd,
	                                       const __m512 yd,
	                                       const __m512 zd,
	                                       const float arg[10],
	                                       std::complex<float> & hx,                        
                                               std::complex<float> & hy,
                                               std::complex<float> & hz,
                                               float err[6]) { // integration error.
                            
                        
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        register float er1,er2,er3,er4,er5,er6;
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        _mm512_store_ps(&dx[0],xd);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        _mm512_store_ps(&dy[0],yd);
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        _mm512_store_ps(&dz[0],zd);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(ntab,&dx[0],&pxr[0],xa,xb,sxr,er1);
                        cubint(ntab,&dx[0],&pxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(ntab,&dy[0],&pyr[0],ya,yb,syr,er3);
                        cubint(ntab,&dy[0],&pyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(ntab,&dz[0],&pzr[0],za,zb,szr,er5);
                        cubint(ntab,&dz[0],&pzi[0],za,zb,szi,er6);
                        err[4] = er5;
                        err[5] = er6;
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_cubint_a(const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          float * __restrict __ATTR_ALIGN__(64)  pxd,
	                                          float * __restrict __ATTR_ALIGN__(64)  pyd,
	                                          float * __restrict __ATTR_ALIGN__(64)  pzd,
	                                          const float arg[10],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  float err[6]) { // integration error.
                            
                        
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi   
                        constexpr int32_t ntab = 16; 
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        register float er1,er2,er3,er4,er5,er6;
                        xr = _mm512_load_ps(&pxre[0]);
                        xi = _mm512_load_ps(&pxim[0]);
                        yr = _mm512_load_ps(&pyre[0]);
                        yi = _mm512_load_ps(&pyim[0]);
                        zr = _mm512_load_ps(&pzre[0]);
                        zi = _mm512_load_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(ntab,&pxd[0],&pxr[0],xa,xb,sxr,er1);
                        cubint(ntab,&pxd[0],&pxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(ntab,&pyd[0],&pyr[0],ya,yb,syr,er3);
                        cubint(ntab,&pyd[0],&pyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(ntab,&pzd[0],&pzr[0],za,zb,szr,er5);
                        cubint(ntab,&pzd[0],&pzi[0],za,zb,szi,er6);
                        err[4] = er5;
                        err[5] = er6;
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvemm_f2135_zmm16r4_cubint_u6x(const float * __restrict  pxre,
	                                               const float * __restrict  pxim,
	                                               const float * __restrict  pyre,
	                                               const float * __restrict  pyim,
	                                               const float * __restrict  pzre,
	                                               const float * __restrict  pzim,
	                                               float * __restrict  pxd,
	                                               float * __restrict  pyd,
	                                               float * __restrict  pzd,
	                                               fwork_t fw,
	                                               const float arg[10],
	                                               std::complex<float> & hx,                        
                                                       std::complex<float> & hy,
                                                       std::complex<float> & hz,
                                                       float err[6],
                                                       const int32_t n,
                                                       int32_t & PF_DIST,
                                                       const bool aligned) { 
                                                  
                        
                         constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi  
                        //register __m512 xr,xi,yr,yi,zr,zi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        std::complex<float> tmp1,tmp2;
                        register float k,r,xa,xb,ya,yb,za,zb,inv;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        register float er1,er2,er3,er4,er5,er6; 
                        //int32_t i;   
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        inv  = 1.0f/r;
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        tmp1 = {0.0f,-1.0f*k*r};
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        tmp2 = tmp1*inv;
                        eps  = arg[9];  
                        if(aligned) {
                           f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                     fw,cer,cei,tmp2,n,PF_DIST);
                        } else {
                           f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                     fw,cer,cei,tmp2,n,PF_DIST);
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(n,&pxd[0],&fw.pxr[0],xa,xb,sxr,er1);
                        cubint(n,&pxd[0],&fw.pxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(n,&pyd[0],&fw.pyr[0],ya,yb,syr,er3);
                        cubint(n,&pyd[0],&fw.pyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(n,&pzd[0],&fw.pzr[0],za,zb,szr,er5);
                        cubint(n,&pzd[0],&fw.pzi[0],za,zb,szi,er6);
                        err[4] = er5;
                        err[5] = er6;
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                           
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_cubint_u(const float * __restrict  pxre,
	                                          const float * __restrict  pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          float * __restrict  pxd,
	                                          float * __restrict  pyd,
	                                          float * __restrict  pzd,
	                                          const float arg[10],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  float err[6]) { // integration error.
                            
                        
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi   
                        constexpr int32_t ntab = 16; 
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        register float er1,er2,er3,er4,er5,er6;
                        xr = _mm512_loadu_ps(&pxre[0]);
                        xi = _mm512_loadu_ps(&pxim[0]);
                        yr = _mm512_loadu_ps(&pyre[0]);
                        yi = _mm512_loadu_ps(&pyim[0]);
                        zr = _mm512_loadu_ps(&pzre[0]);
                        zi = _mm512_loadu_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(ntab,&pxd[0],&pxr[0],xa,xb,sxr,er1);
                        cubint(ntab,&pxd[0],&pxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(ntab,&pyd[0],&pyr[0],ya,yb,syr,er3);
                        cubint(ntab,&pyd[0],&pyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(ntab,&pzd[0],&pzr[0],za,zb,szr,er5);
                        cubint(ntab,&pzd[0],&pzi[0],za,zb,szi,er6);
                        err[4] = er5;
                        err[5] = er6;
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
                 
               
                /*
	            Hertz vector (electrical,magnetic), hiordq integrator.
	            Formula 2-13, p. 35
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_hiordq(const __m512 xre,
	                                       const __m512 xim,
	                                       const __m512 yre,
	                                       const __m512 yim,
	                                       const __m512 zre,
	                                       const __m512 zim,
	                                       const float arg[7],
	                                       std::complex<float> & hx,                        
                                               std::complex<float> & hy,
                                               std::complex<float> & hz) {
                                               
                                              
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        __ATTR_ALIGN__(64) float work[32];
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,deltx,delty,deltz;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        deltx = arg[2];
                        delty = arg[3];
                        eai   = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        deltz = arg[4];
                        omg   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        eps  = arg[6];
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&pxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&pxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&pyr[0],&work[0],syr);
                        hiordq(ntab,delty,&pyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&pzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&pzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_hiordq_a(const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          const float arg[7],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                               
                                              
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        __ATTR_ALIGN__(64) float work[32];
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,deltx,delty,deltz;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_load_ps(&pxre[0]);
                        xi = _mm512_load_ps(&pxim[0]);
                        yr = _mm512_load_ps(&pyre[0]);
                        yi = _mm512_load_ps(&pyim[0]);
                        zr = _mm512_load_ps(&pzre[0]);
                        zi = _mm512_load_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        deltx = arg[2];
                        delty = arg[3];
                        eai   = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        deltz = arg[4];
                        omg   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        eps  = arg[6];
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&pxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&pxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&pyr[0],&work[0],syr);
                        hiordq(ntab,delty,&pyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&pzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&pzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_hiordq_u6x(const float * __restrict  pxre,
	                                          const float * __restrict  pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          float * __restrict  work, // size of work is 2*(n-1)
	                                          fwork_t fw,
	                                          const float arg[7],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  const int32_t n,
                                                  const int32_t PF_DIST,
                                                  const bool aligned) { 
                                                  
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register float k,r,deltx,delty,deltz,inv;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;  
                        std::complex<float> tmp1,tmp2;
                        k     = arg[0];
                        r     = arg[1];
                        inv   =  1.0f/r;
                        vk    = _mm512_set1_ps(k);
                        vr    = _mm512_set1_ps(r);
                        ir    = _mm512_setzero_ps();
                        invr  = _mm512_rcp14_ps(vr);
                        ii    = _mm512_set1_ps(-1.0f);
                        deltx = arg[2];
                        tmp1  = {0.0f,-1.0f*k*r};
                        delty = arg[3];
                        eai   = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        deltz = arg[4];
                        omg   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        tmp2  = tmp1*inv;
                        eps   = arg[6];
                        cer   = _mm512_mul_ps(cer,invr);
                        cei   = _mm512_mul_ps(cei,invr); 
                        if(aligned) {
                           f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                        fw,cer,cei,tmp2,n,PF_DIST);
                        } else {
                           f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                        fw,cer,cei,tmp2,n,PF_DIST); 
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&fw.pxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&fw.pxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&fw.pyr[0],&work[0],syr);
                        hiordq(ntab,delty,&fw.pyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&fw.pzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&fw.pzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                              
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_hiordq_u(const float * __restrict  pxre,
	                                          const float * __restrict  pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          const float arg[7],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                               
                                              
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        float work[32];
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,deltx,delty,deltz;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_loadu_ps(&pxre[0]);
                        xi = _mm512_loadu_ps(&pxim[0]);
                        yr = _mm512_loadu_ps(&pyre[0]);
                        yi = _mm512_loadu_ps(&pyim[0]);
                        zr = _mm512_loadu_ps(&pzre[0]);
                        zi = _mm512_loadu_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        deltx = arg[2];
                        delty = arg[3];
                        eai   = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        deltz = arg[4];
                        omg   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        eps  = arg[6];
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&pxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&pxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&pyr[0],&work[0],syr);
                        hiordq(ntab,delty,&pyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&pzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&pzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
                 
               
                /*
	            Hertz vector (electrical,magnetic), plint integrator.
	            Formula 2-13, p. 35
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_plint(const __m512 xre,
	                                       const __m512 xim,
	                                       const __m512 yre,
	                                       const __m512 yim,
	                                       const __m512 zre,
	                                       const __m512 zim,
	                                       const __m512 xd,
	                                       const __m512 yd,
	                                       const __m512 zd,
	                                       const float arg[10],
	                                       std::complex<float> & hx,                        
                                               std::complex<float> & hy,
                                               std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        _mm512_store_ps(&dx[0],xd);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        _mm512_store_ps(&dy[0],yd);
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        _mm512_store_ps(&dz[0],zd);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        plint(ntab,&dx[0],&pxr[0],xa,xb,sxr);
                        plint(ntab,&dx[0],&pxi[0],xa,xb,sxi);
                        plint(ntab,&dy[0],&pyr[0],ya,yb,syr);
                        plint(ntab,&dy[0],&pyi[0],ya,yb,syi);
                        plint(ntab,&dz[0],&pzr[0],za,zb,szr);
                        plint(ntab,&dz[0],&pzi[0],za,zb,szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
	       
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_plint_a(const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          float * __restrict __ATTR_ALIGN__(64)  pxd,
	                                          float * __restrict __ATTR_ALIGN__(64)  pyd,
	                                          float * __restrict __ATTR_ALIGN__(64)  pzd,
	                                          const float arg[10],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_load_ps(&pxre[0]);
                        xi = _mm512_load_ps(&pxim[0]);
                        yr = _mm512_load_ps(&pyre[0]);
                        yi = _mm512_load_ps(&pyim[0]);
                        zr = _mm512_load_ps(&pzre[0]);
                        zi = _mm512_load_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        plint(ntab,&pxd[0],&pxr[0],xa,xb,sxr);
                        plint(ntab,&pxd[0],&pxi[0],xa,xb,sxi);
                        plint(ntab,&pyd[0],&pyr[0],ya,yb,syr);
                        plint(ntab,&pyd[0],&pyi[0],ya,yb,syi);
                        plint(ntab,&pzd[0],&pzr[0],za,zb,szr);
                        plint(ntab,&pzd[0],&pzi[0],za,zb,szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_plint_u6x(const float * __restrict  pxre,
	                                             const float * __restrict  pxim,
	                                             const float * __restrict  pyre,
	                                             const float * __restrict  pyim,
	                                             const float * __restrict  pzre,
	                                             const float * __restrict pzim,
	                                             float * __restrict   pxd,
	                                             float * __restrict   pyd,
	                                             float * __restrict   pzd,
	                                             fwork_t fw,
	                                             const float arg[10],
	                                             std::complex<float> & hx,                        
                                                     std::complex<float> & hy,
                                                     std::complex<float> & hz,
                                                     const int32_t n,
                                                     const int32_t PF_DIST,
                                                     const bool aligned) {
                                                     
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        std::complex<float> tmp1,tmp2;
                        register float k,r,xa,xb,ya,yb,za,zb,inv;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;    
                        k    = arg[0];
                        r    = arg[1];
                        inv  = 1.0f/r;
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        tmp1 = {0.0f,-1.0f*k*r};
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        tmp2 = tmp1*inv;
                        eps  = arg[9];  
                        if(aligned) {
                           f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                        fw,cer,cei,tmp2,n,PF_DIST);
                        } else {
                           f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                        fw,cer,cei,tmp2,n,PF_DIST);
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        plint(n,&pxd[0],&fw.pxr[0],xa,xb,sxr);
                        plint(n,&pxd[0],&fw.pxi[0],xa,xb,sxi);
                        plint(n,&pyd[0],&fw.pyr[0],ya,yb,syr);
                        plint(n,&pyd[0],&fw.pyi[0],ya,yb,syi);
                        plint(n,&pzd[0],&fw.pzr[0],za,zb,szr);
                        plint(n,&pzd[0],&fw.pzi[0],za,zb,szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                                         
              }
                                              
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_plint_u(const float * __restrict  pxre,
	                                          const float * __restrict   pxim,
	                                          const float * __restrict   pyre,
	                                          const float * __restrict   pyim,
	                                          const float * __restrict   pzre,
	                                          const float * __restrict   pzim,
	                                          float * __restrict   pxd,
	                                          float * __restrict   pyd,
	                                          float * __restrict   pzd,
	                                          const float arg[10],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_loadu_ps(&pxre[0]);
                        xi = _mm512_loadu_ps(&pxim[0]);
                        yr = _mm512_loadu_ps(&pyre[0]);
                        yi = _mm512_loadu_ps(&pyim[0]);
                        zr = _mm512_loadu_ps(&pzre[0]);
                        zi = _mm512_loadu_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        xa   = arg[2];
                        xb   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        ya   = arg[4];
                        yb   = arg[5];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        za   = arg[6];
                        zb   = arg[7];
                        cer  = _mm512_mul_ps(cer,invr);
                        omg  = arg[8];
                        cei  = _mm512_mul_ps(cei,invr);
                        eps  = arg[9];
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        plint(ntab,&pxd[0],&pxr[0],xa,xb,sxr);
                        plint(ntab,&pxd[0],&pxi[0],xa,xb,sxi);
                        plint(ntab,&pyd[0],&pyr[0],ya,yb,syr);
                        plint(ntab,&pyd[0],&pyi[0],ya,yb,syi);
                        plint(ntab,&pzd[0],&pzr[0],za,zb,szr);
                        plint(ntab,&pzd[0],&pzi[0],za,zb,szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
                  
               
                /*
	            Hertz vector (electrical,magnetic), simpne integrator.
	            Formula 2-13, p. 35
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpne(const __m512 xre,
	                                          const __m512 xim,
	                                          const __m512 yre,
	                                          const __m512 yim,
	                                          const __m512 zre,
	                                          const __m512 zim,
	                                          const __m512 xd,
	                                          const __m512 yd,
	                                          const __m512 zd,
	                                          const float arg[4],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        _mm512_store_ps(&dx[0],xd);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        _mm512_store_ps(&dy[0],yd);
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        _mm512_store_ps(&dz[0],zd);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpne(ntab,&dx[0],&pxr[0],sxr);
                        simpne(ntab,&dx[0],&pxi[0],sxi);
                        simpne(ntab,&dy[0],&pyr[0],syr);
                        simpne(ntab,&dy[0],&pyi[0],syi);
                        simpne(ntab,&dz[0],&pzr[0],szr);
                        simpne(ntab,&dz[0],&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpne_a(const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          float * __restrict __ATTR_ALIGN__(64)  pxd,
	                                          float * __restrict __ATTR_ALIGN__(64)  pyd,
	                                          float * __restrict __ATTR_ALIGN__(64)  pzd,
	                                          const float arg[4],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_load_ps(&pxre[0]);
                        xi = _mm512_load_ps(&pxim[0]);
                        yr = _mm512_load_ps(&pyre[0]);
                        yi = _mm512_load_ps(&pyim[0]);
                        zr = _mm512_load_ps(&pzre[0]);
                        zi = _mm512_load_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpne(ntab,&pxd[0],&pxr[0],sxr);
                        simpne(ntab,&pxd[0],&pxi[0],sxi);
                        simpne(ntab,&pyd[0],&pyr[0],syr);
                        simpne(ntab,&pyd[0],&pyi[0],syi);
                        simpne(ntab,&pzd[0],&pzr[0],szr);
                        simpne(ntab,&pzd[0],&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpne_u6x(const float * __restrict  pxre,
	                                          const float * __restrict  pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          float * __restrict pxd,
	                                          float * __restrict pyd,
	                                          float * __restrict pzd,
	                                          fwork_t fw,
	                                          const float arg[4],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  const int32_t n,
                                                  const int32_t PF_DIST,
                                                  const bool aligned) {
                                                  
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        std::complex<float> tmp1,tmp2;
                        register float k,r,omg,eps,inv;
                        register float sxr,sxi,syr,syi,szr,szi,frac;     
                        k    = arg[0];
                        r    = arg[1];
                        inv  = 1.0f/r;
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        tmp1 = {0.0f,-1.0f*k*r};
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        tmp2 = tmp1*inv;
                        cei  = _mm512_mul_ps(cei,invr);  
                        if(aligned) {
                           f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                     fw,cer,cei,tmp2,n,PF_DIST);
                        } else {
                           f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                     fw,cer,cei,tmp2,n,PF_DIST); 
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpne(n,&pxd[0],&fw.pxr[0],sxr);
                        simpne(n,&pxd[0],&fw.pxi[0],sxi);
                        simpne(n,&pyd[0],&fw.pyr[0],syr);
                        simpne(n,&pyd[0],&fw.pyi[0],syi);
                        simpne(n,&pzd[0],&fw.pzr[0],szr);
                        simpne(n,&pzd[0],&fw.pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                                
            }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpne_u(const float * __restrict  pxre,
	                                          const float * __restrict   pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          float * __restrict   pxd,
	                                          float * __restrict   pyd,
	                                          float * __restrict   pzd,
	                                          const float arg[4],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_loadu_ps(&pxre[0]);
                        xi = _mm512_loadu_ps(&pxim[0]);
                        yr = _mm512_loadu_ps(&pyre[0]);
                        yi = _mm512_loadu_ps(&pyim[0]);
                        zr = _mm512_loadu_ps(&pzre[0]);
                        zi = _mm512_loadu_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpne(ntab,&pxd[0],&pxr[0],sxr);
                        simpne(ntab,&pxd[0],&pxi[0],sxi);
                        simpne(ntab,&pyd[0],&pyr[0],syr);
                        simpne(ntab,&pyd[0],&pyi[0],syi);
                        simpne(ntab,&pzd[0],&pzr[0],szr);
                        simpne(ntab,&pzd[0],&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               /*
                    Hertz vector (electrical,magnetic), simpn integrator.
	            Formula 2-13, p. 35  
               */
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpn(const __m512 xre,
	                                       const __m512 xim,
	                                       const __m512 yre,
	                                       const __m512 yim,
	                                       const __m512 zre,
	                                       const __m512 zim,
	                                       const float arg[5],
	                                       std::complex<float> & hx,                        
                                               std::complex<float> & hy,
                                               std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps,h;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpn(ntab,h,&pxr[0],sxr);
                        simpn(ntab,h,&pxi[0],sxi);
                        simpn(ntab,h,&pyr[0],syr);
                        simpn(ntab,h,&pyi[0],syi);
                        simpn(ntab,h,&pzr[0],szr);
                        simpn(ntab,h,&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpn_a( const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          const float arg[5],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps,h;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_load_ps(&pxre[0]);
                        xi = _mm512_load_ps(&pxim[0]);
                        yr = _mm512_load_ps(&pyre[0]);
                        yi = _mm512_load_ps(&pyim[0]);
                        zr = _mm512_load_ps(&pzre[0]);
                        zi = _mm512_load_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpn(ntab,h,&pxr[0],sxr);
                        simpn(ntab,h,&pxi[0],sxi);
                        simpn(ntab,h,&pyr[0],syr);
                        simpn(ntab,h,&pyi[0],syi);
                        simpn(ntab,h,&pzr[0],szr);
                        simpn(ntab,h,&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpn_u6x( const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          fwork_t fw,
	                                          const float arg[5],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  const int32_t n,
                                                  const int32_t PF_DIST,
                                                  const bool aligned) {
                                                  
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        std::complex<float> tmp1,tmp2;
                        register float k,r,omg,eps,h,inv;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        k    = arg[0];
                        r    = arg[1];
                        inv  = 1.0f/r;
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg  = arg[2];
                        tmp2 = {0.0f,-1.0f*k*r};
                        eps  = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        tmp2 = tmp1*inv;
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);  
                        if(aligned) {
                           f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                     fw,cer,cei,tmp2,n,PF_DIST);
                        } else {
                           f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                     fw,cer,cei,tmp2,n,PF_DIST);
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpn(n,h,&fw.pxr[0],sxr);
                        simpn(n,h,&fw.pxi[0],sxi);
                        simpn(n,h,&fw.pyr[0],syr);
                        simpn(n,h,&fw.pyi[0],syi);
                        simpn(n,h,&fw.pzr[0],szr);
                        simpn(n,h,&fw.pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                                            
             }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_simpn_u( const float * __restrict pxre,
	                                          const float * __restrict   pxim,
	                                          const float * __restrict   pyre,
	                                          const float * __restrict   pyim,
	                                          const float * __restrict   pzre,
	                                          const float * __restrict   pzim,
	                                          const float arg[5],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps,h;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_loadu_ps(&pxre[0]);
                        xi = _mm512_loadu_ps(&pxim[0]);
                        yr = _mm512_loadu_ps(&pyre[0]);
                        yi = _mm512_loadu_ps(&pyim[0]);
                        zr = _mm512_loadu_ps(&pzre[0]);
                        zi = _mm512_loadu_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        simpn(ntab,h,&pxr[0],sxr);
                        simpn(ntab,h,&pxi[0],sxi);
                        simpn(ntab,h,&pyr[0],syr);
                        simpn(ntab,h,&pyi[0],syi);
                        simpn(ntab,h,&pzr[0],szr);
                        simpn(ntab,h,&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
                 
               
                /*
                    Hertz vector (electrical,magnetic), wedint integrator.
	            Formula 2-13, p. 35  
               */
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_wedint(const __m512 xre,
	                                       const __m512 xim,
	                                       const __m512 yre,
	                                       const __m512 yim,
	                                       const __m512 zre,
	                                       const __m512 zim,
	                                       const float arg[5],
	                                       std::complex<float> & hx,                        
                                               std::complex<float> & hy,
                                               std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps,h;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                     
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg  = arg[2];
                        eps  = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        wedint(ntab,h,&pxr[0],sxr);
                        wedint(ntab,h,&pxi[0],sxi);
                        wedint(ntab,h,&pyr[0],syr);
                        wedint(ntab,h,&pyi[0],syi);
                        wedint(ntab,h,&pzr[0],szr);
                        wedint(ntab,h,&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_wedint_a( const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          const float arg[5],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps,h;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_load_ps(&pxre[0]);
                        xi = _mm512_load_ps(&pxim[0]);
                        yr = _mm512_load_ps(&pyre[0]);
                        yi = _mm512_load_ps(&pyim[0]);
                        zr = _mm512_load_ps(&pzre[0]);
                        zi = _mm512_load_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        wedint(ntab,h,&pxr[0],sxr);
                        wedint(ntab,h,&pxi[0],sxi);
                        wedint(ntab,h,&pyr[0],syr);
                        wedint(ntab,h,&pyi[0],syi);
                        wedint(ntab,h,&pzr[0],szr);
                        wedint(ntab,h,&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_wedint_u6x( const float * __restrict __ATTR_ALIGN__(64)  pxre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pxim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pyim,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzre,
	                                          const float * __restrict __ATTR_ALIGN__(64)  pzim,
	                                          fwork_t fw,
	                                          const float arg[5],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  const int32_t n,
                                                  const int32_t PF_DIST,
                                                  const bool aligned) {
                                                  
                          constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                          register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                          std::complex<float> tmp1,tmp2;
                          register float k,r,omg,eps,h,inv;
                          register float sxr,sxi,syr,syi,szr,szi,frac;    
                          k   = arg[0];
                          r   = arg[1];
                          inv = 1.0f/r;
                          vk  = _mm512_set1_ps(k);
                          vr  = _mm512_set1_ps(r);
                          ir  = _mm512_setzero_ps();
                          invr= _mm512_rcp14_ps(vr);
                          ii  = _mm512_set1_ps(-1.0f);
                          omg = arg[2];
                          tmp1= {0.0f,-1.0f*k*r};
                          eps = arg[3];
                          eai = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                          h   = arg[4];
                          tmp2 = tmp1*inv;
                          cexp_zmm16r4(ir,eai,&cer,&cei);
                          cer = _mm512_mul_ps(cer,invr);
                          cei = _mm512_mul_ps(cei,invr);  
                          if(aligned) {
                             f2135_integrand_zmm16r4_u6x_a(pxre,pxim,pyre,pyim,pzre,pzim,
                                                         fw,cer,cei,tmp2,n,PF_DIST);
                          } else {
                             f2135_integrand_zmm16r4_u6x_u(pxre,pxim,pyre,pyim,pzre,pzim,
                                                         fw,cer,cei,tmp2,n,PF_DIST);
                          }
                          sxr = 0.0f;
                          sxi = sxr;
                          syi = sxr;
                          syr = sxr;
                          szr = sxr;
                          szi = sxr;
                          float tmp = C12566370614359172953850573533118*omg*eps;
                          frac = 1.0f/tmp;
                          wedint(n,h,&fw.pxr[0],sxr);
                          wedint(n,h,&fw.pxi[0],sxi);
                          wedint(n,h,&fw.pyr[0],syr);
                          wedint(n,h,&fw.pyi[0],syi);
                          wedint(n,h,&fw.pzr[0],szr);
                          wedint(n,h,&fw.pzi[0],szi);
                          hx = {sxr*frac,sxi*frac};
                          hy = {syr*frac,syi*frac};
                          hz = {szr*frac,szi*frac};                                        
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hvem_f2135_zmm16r4_wedint_u( const float * __restrict  pxre,
	                                          const float * __restrict  pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          const float arg[5],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                              
                                           
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;
                        register float k,r,omg,eps,h;
                        register float sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_loadu_ps(&pxre[0]);
                        xi = _mm512_loadu_ps(&pxim[0]);
                        yr = _mm512_loadu_ps(&pyre[0]);
                        yi = _mm512_loadu_ps(&pyim[0]);
                        zr = _mm512_loadu_ps(&pzre[0]);
                        zi = _mm512_loadu_ps(&pzim[0]);
                        k = arg[0];
                        r = arg[1];
                        vk   = _mm512_set1_ps(k);
                        vr   = _mm512_set1_ps(r);
                        ir   = _mm512_setzero_ps();
                        invr = _mm512_rcp14_ps(vr);
                        ii   = _mm512_set1_ps(-1.0f);
                        omg   = arg[2];
                        eps   = arg[3];
                        eai  = _mm512_mul_ps(ii,_mm512_mul_ps(vk,vr));
                        h    = arg[4];
                        cexp_zmm16r4(ir,eai,&cer,&cei);
                        cer  = _mm512_mul_ps(cer,invr);
                        cei  = _mm512_mul_ps(cei,invr);
                        cmul_zmm16r4(xre,xim,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0]
                        cmul_zmm16r4(yre,yim,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(zre,zim,cer,cei,&intzr,&intzi);
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        wedint(ntab,h,&pxr[0],sxr);
                        wedint(ntab,h,&pxi[0],sxi);
                        wedint(ntab,h,&pyr[0],syr);
                        wedint(ntab,h,&pyi[0],syi);
                        wedint(ntab,h,&pzr[0],szr);
                        wedint(ntab,h,&pzi[0],szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
               /*
                   Formula 2.21, p. 36
               */
               
                __ATTR_ALWAYS_INLINE__
	        __ATTR_HOT__
	        static inline
	        bool f221_r4(const float R,
	                     const float D,
	                     const float gam) {
	            
	            register float DD  = D*D;
	            register gloat rat = DD/gam;
	            bool bres          = R>=(2.0*rat);
	            return (bres);               
	      }
	      
	      
	       /*
                   Formula 2.21, p. 36
                   (AVX512) version.
               */
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __mmask16 f221_zmm16r4(const __m512 R,
                                          const float D,
                                          const float gam) {
                     
                        
                        register __m512 vD,vgam,DD,rat,t0;
                        __mmask16 res = 0x0;
                        vD  = _mm512_set1_ps(D);
                        DD  = _mm512_mul_ps(vD,vD);
                        vgam= _mm512_set1_ps(gam);
                        rat = _mm512_div_ps(DD,vgam);
                        t0  = _mm512_add_ps(rat,rat);
                        res = _mm512_cmp_ps_mask(R,t0,_CMP_GE_OQ);
                        return (res);                        
                 }
                 
                 
                 /*
                     Formula 2-23, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'Avint' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_avint(const __m512 jxr,
	                                        const __m512 jxi,
	                                        const __m512 jyr,
	                                        const __m512 jyi,
	                                        const __m512 jzr,
	                                        const __m512 jzi,
	                                       __m512 xd,
	                                       __m512 yd,
	                                       __m512 zd,
	                                       const __m512 rho,
	                                       const __m512 cost,
	                                       const float args[7],
	                                       std::complex<float> & Nx,
	                                       std::complex<float> & Ny,
	                                       std::complex<float> & Nz,
	                                       int32_t & ierr) {
	                                       
	                register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cost));
                        ya  = args[3];
                        yb  = args[4];
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        sxr = avint(&pxd[0],&pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {goto ERROR;}
                           goto CORRECT;
                        syr = avint(&pyd[0],&pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {goto ERROR;}
                           goto CORRECT;
                        szr = avint(&pzd[0],&pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {goto ERROR;}
                           goto CORRECT;
                        ERROR:
                           {
                               ierr = 3;
                               return;
                        }  
                        CORRECT: {
                           Nx = {sxr,sxi};
                           Ny = {syr,syi};
                           Nz = {szr,szi};  
                        }                
	        }
	        
	        
	         
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_avint_a(const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                         float * __restrict __ATTR_ALIGN__(64) pxd,
	                                         float * __restrict __ATTR_ALIGN__(64) pyd,
	                                         float * __restrict __ATTR_ALIGN__(64) pzd,
	                                         const float * __restrict __ATTR_ALIGN__(64) prho,
	                                         const float * __restrict __ATTR_ALIGN__(64) pcost,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz,
	                                         int32_t & ierr ) {
	                                         
	                register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,rho,cost;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                                              
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        cost= _mm512_load_ps(&pcost[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        rho = _mm512_load_ps(&prho[0]);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        jxr = _mm512_load_ps(&pjxr[0]);
                        jxi = _mm512_load_ps(&pjxi[0]);
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cost));
                        jyr = _mm512_load_ps(&pjyr[0]);
                        jyi = _mm512_load_ps(&pjyi[0]);
                        ya  = args[3];
                        yb  = args[4];
                        jzr = _mm512_load_ps(&pjzr[0]);
                        jzi = _mm512_load_ps(&pjzi[0]);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        sxr = avint(&pxd[0],&pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {goto ERROR;}
                           goto CORRECT;
                        syr = avint(&pyd[0],&pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {goto ERROR;}
                           goto CORRECT;
                        szr = avint(&pzd[0],&pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {goto ERROR;}
                           goto CORRECT;
                        ERROR:
                           {
                               ierr = 3;
                               return;
                        }  
                        CORRECT: {
                           Nx = {sxr,sxi};
                           Ny = {syr,syi};
                           Nz = {szr,szi};  
                        }                                                 
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void f2235_integrand_zmm16r4_u6x_a(const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                              const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                              const float * __restrict __ATTR_ALIGN__(64) prho,
	                                              const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                              fwork_t fw, ////work arrays (caller allocated)
	                                              const float k,
	                                              const int32_t n,
	                                              const int32_t PF_DIST) {
	                                              
	                 if(__builtin_expect(n<=0,0)) { return;}
	                 register __m512 jxr,jxi,jyr,jyi,jzr,jzi;
	                 register __m512 rho,cst,ir,ii;
	                 register __m512 ear,eai,cer,cei;
	                 register __m512 vk,t0;
	                 __m512   t0r,t0i,t1r,t1i,t2r,t2i;
	                 int32_t i;
	                 
	                 vk = _mm512_set1_ps(k);
	                 ii = _mm512_set1_ps(1.0f);
	                 ir = _mm512_setzero_ps();
	                 ear= ir;
	                 t0 = _mm512_mul_ps(ii,vk); // shall be removed.
	                 for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_NTA);
#endif	                 
                              cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_load_ps(&pcst[i+16]);
                              rho = _mm512_load_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+16]);
                              jxi = _mm512_load_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+16],t0r);
                              _mm512_store_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+16]);
                              jyi = _mm512_load_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+16],t1r);
                              _mm512_store_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+16]);
                              jzi = _mm512_load_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+16],t2r);
                              _mm512_store_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_load_ps(&pcst[i+32]);
                              rho = _mm512_load_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+32]);
                              jxi = _mm512_load_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+32],t0r);
                              _mm512_store_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+32]);
                              jyi = _mm512_load_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+32],t1r);
                              _mm512_store_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+32]);
                              jzi = _mm512_load_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+32],t2r);
                              _mm512_store_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_load_ps(&pcst[i+48]);
                              rho = _mm512_load_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+48]);
                              jxi = _mm512_load_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+48],t0r);
                              _mm512_store_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+48]);
                              jyi = _mm512_load_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+48],t1r);
                              _mm512_store_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+48]);
                              jzi = _mm512_load_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+48],t2r);
                              _mm512_store_ps(&fw.pzi[i+48],t2i);
                              cst = _mm512_load_ps(&pcst[i+64]);
                              rho = _mm512_load_ps(&prho[i+64]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+64]);
                              jxi = _mm512_load_ps(&pjxi[i+64]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+64],t0r);
                              _mm512_store_ps(&fw.pxi[i+64],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+64]);
                              jyi = _mm512_load_ps(&pjyi[i+64]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+64],t1r);
                              _mm512_store_ps(&fw.pyi[i+64],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+64]);
                              jzi = _mm512_load_ps(&pjzi[i+64]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+64],t2r);
                              _mm512_store_ps(&fw.pzi[i+64],t2i);
                              cst = _mm512_load_ps(&pcst[i+80]);
                              rho = _mm512_load_ps(&prho[i+80]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+80]);
                              jxi = _mm512_load_ps(&pjxi[i+80]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+80],t0r);
                              _mm512_store_ps(&fw.pxi[i+80],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+80]);
                              jyi = _mm512_load_ps(&pjyi[i+80]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+80],t1r);
                              _mm512_store_ps(&fw.pyi[i+80],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+80]);
                              jzi = _mm512_load_ps(&pjzi[i+80]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+80],t2r);
                              _mm512_store_ps(&fw.pzi[i+80],t2i);
	                 }  
	                 
	                 for(; (i+63) < n; i += 64) {
	                      cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_load_ps(&pcst[i+16]);
                              rho = _mm512_load_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+16]);
                              jxi = _mm512_load_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+16],t0r);
                              _mm512_store_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+16]);
                              jyi = _mm512_load_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+16],t1r);
                              _mm512_store_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+16]);
                              jzi = _mm512_load_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+16],t2r);
                              _mm512_store_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_load_ps(&pcst[i+32]);
                              rho = _mm512_load_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+32]);
                              jxi = _mm512_load_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+32],t0r);
                              _mm512_store_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+32]);
                              jyi = _mm512_load_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+32],t1r);
                              _mm512_store_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+32]);
                              jzi = _mm512_load_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+32],t2r);
                              _mm512_store_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_load_ps(&pcst[i+48]);
                              rho = _mm512_load_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+48]);
                              jxi = _mm512_load_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+48],t0r);
                              _mm512_store_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+48]);
                              jyi = _mm512_load_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+48],t1r);
                              _mm512_store_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+48]);
                              jzi = _mm512_load_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+48],t2r);
                              _mm512_store_ps(&fw.pzi[i+48],t2i);
	                 }   
	                 
	                 for(; (i+31) < n; i += 32) {
	                      cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_load_ps(&pcst[i+16]);
                              rho = _mm512_load_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+16]);
                              jxi = _mm512_load_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+16],t0r);
                              _mm512_store_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+16]);
                              jyi = _mm512_load_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+16],t1r);
                              _mm512_store_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+16]);
                              jzi = _mm512_load_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+16],t2r);
                              _mm512_store_ps(&fw.pzi[i+16],t2i);
	                 }  
	                 
	                 for(; (i+15) < n; i += 16) {
	                      cst = _mm512_load_ps(&pcst[i+0]);
                              rho = _mm512_load_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_load_ps(&pjxr[i+0]);
                              jxi = _mm512_load_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_store_ps(&fw.pxr[i+0],t0r);
                              _mm512_store_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_load_ps(&pjyr[i+0]);
                              jyi = _mm512_load_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_store_ps(&fw.pyr[i+0],t1r);
                              _mm512_store_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_load_ps(&pjzr[i+0]);
                              jzi = _mm512_load_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_store_ps(&fw.pzr[i+0],t2r);
                              _mm512_store_ps(&fw.pzi[i+0],t2i);
	                 }  
	                 
	                 for(; (i+0) < n; i += 1) {
	                     float cst               = pcst[i];
	                     float rho               = prho[i];
	                     float eai               = k*cst*rho;
	                     std::complex<float> ce = std::cexp({0.0f,eai});
	                     float jxr               = pjxr[i];
	                     float jxi               = pjxi[i];
	                     std::complex<float> cx  = {jxr,jxi};
	                     std::complex<float> jx  = cx*ce;
	                     fw.pxr[i]               = jx.real();
	                     fw.pxi[i]               = jx.imag();
	                     float jyr               = pjyr[i];
	                     float jyi               = pjyi[i];
	                     std::complex<float> cy  = {jyr,jyi};
	                     std::complex<float> jy  = cy*ce;
	                     fw.pyr[i]               = jy.real();
	                     fw.pyi[i]               = jy.imag();
	                     float jzr               = pjzr[i];
	                     float jzi               = pjzi[i];
	                     std::complex<float> cz  = {jzr,jzi};
	                     std::complex<float> jz  = ce*cz;
	                     fw.pzr[i]               = jz.real();
	                     fw.pzi[i]               = jz.imag();
	                 }                             
	       } 
	        
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void f2235_integrand_zmm16r4_u6x_u(const float * __restrict  pjxr,
	                                              const float * __restrict  pjxi,
	                                              const float * __restrict  pjyr,
	                                              const float * __restrict  pjyi,
	                                              const float * __restrict  pjzr,
	                                              const float * __restrict  pjzi,
	                                              const float * __restrict  prho,
	                                              const float * __restrict  pcst,
	                                              fwork_t fw, ////work arrays (caller allocated)
	                                              const float k,
	                                              const int32_t n,
	                                              const int32_t PF_DIST) {
	                                              
	                 if(__builtin_expect(n<=0,0)) { return;}
	                 register __m512 jxr,jxi,jyr,jyi,jzr,jzi;
	                 register __m512 rho,cst,ir,ii;
	                 register __m512 ear,eai,cer,cei;
	                 register __m512 vk,t0;
	                 __m512   t0r,t0i,t1r,t1i,t2r,t2i;
	                 int32_t i;
	                 
	                 vk = _mm512_set1_ps(k);
	                 ii = _mm512_set1_ps(1.0f);
	                 ir = _mm512_setzero_ps();
	                 ear= ir;
	                 t0 = _mm512_mul_ps(ii,vk); // shall be removed.
	                 for(i = 0; (i+95) < n; i += 96) {
#if (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 1
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T0);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T0);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 2                       
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T1);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T1);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 3
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_T2);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_T2);
#elif (__ANTENNA_FEEDER_PF_CACHE_HINT__) == 4
                             _mm_prefetch((char*)&pjxr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjxi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjyi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzr[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pjzi[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&prho[i+PF_DIST],_MM_HINT_NTA);
                             _mm_prefetch((char*)&pcst[i+PF_DIST],_MM_HINT_NTA);
#endif	                 
                              cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+16]);
                              rho = _mm512_loadu_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+16]);
                              jxi = _mm512_loadu_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+16]);
                              jyi = _mm512_loadu_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+16]);
                              jzi = _mm512_loadu_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+16],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+32]);
                              rho = _mm512_loadu_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+32]);
                              jxi = _mm512_loadu_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+32]);
                              jyi = _mm512_loadu_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+32]);
                              jzi = _mm512_loadu_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+32],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+48]);
                              rho = _mm512_loadu_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+48]);
                              jxi = _mm512_loadu_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+48]);
                              jyi = _mm512_loadu_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+48]);
                              jzi = _mm512_loadu_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+48],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+48],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+64]);
                              rho = _mm512_loadu_ps(&prho[i+64]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+64]);
                              jxi = _mm512_loadu_ps(&pjxi[i+64]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+64],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+64],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+64]);
                              jyi = _mm512_loadu_ps(&pjyi[i+64]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+64],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+64],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+64]);
                              jzi = _mm512_loadu_ps(&pjzi[i+64]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+64],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+64],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+80]);
                              rho = _mm512_loadu_ps(&prho[i+80]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+80]);
                              jxi = _mm512_loadu_ps(&pjxi[i+80]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+80],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+80],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+80]);
                              jyi = _mm512_loadu_ps(&pjyi[i+80]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+80],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+80],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+80]);
                              jzi = _mm512_loadu_ps(&pjzi[i+80]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+80],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+80],t2i);
	                 }  
	                 
	                 for(; (i+63) < n; i += 64) {
	                      cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+16]);
                              rho = _mm512_loadu_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+16]);
                              jxi = _mm512_loadu_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+16]);
                              jyi = _mm512_loadu_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+16]);
                              jzi = _mm512_loadu_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+16],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+16],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+32]);
                              rho = _mm512_loadu_ps(&prho[i+32]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+32]);
                              jxi = _mm512_loadu_ps(&pjxi[i+32]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+32],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+32],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+32]);
                              jyi = _mm512_loadu_ps(&pjyi[i+32]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+32],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+32],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+32]);
                              jzi = _mm512_loadu_ps(&pjzi[i+32]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+32],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+32],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+48]);
                              rho = _mm512_loadu_ps(&prho[i+48]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+48]);
                              jxi = _mm512_loadu_ps(&pjxi[i+48]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+48],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+48],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+48]);
                              jyi = _mm512_loadu_ps(&pjyi[i+48]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+48],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+48],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+48]);
                              jzi = _mm512_loadu_ps(&pjzi[i+48]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+48],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+48],t2i);
	                 }   
	                 
	                 for(; (i+31) < n; i += 32) {
	                       cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
                              cst = _mm512_loadu_ps(&pcst[i+16]);
                              rho = _mm512_loadu_ps(&prho[i+16]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+16]);
                              jxi = _mm512_loadu_ps(&pjxi[i+16]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+16],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+16],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+16]);
                              jyi = _mm512_loadu_ps(&pjyi[i+16]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+16],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+16],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+16]);
                              jzi = _mm512_loadu_ps(&pjzi[i+16]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+16],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+16],t2i);
	                 }  
	                 
	                 for(; (i+15) < n; i += 16) {
	                      cst = _mm512_loadu_ps(&pcst[i+0]);
                              rho = _mm512_loadu_ps(&prho[i+0]);
                              eai = _mm512_mul_ps(t0,
                                              _mm512_mul_ps(rho,cst));
                              cexp_zmm16r4(ear,eai,&cer,&cei);  
                              jxr = _mm512_loadu_ps(&pjxr[i+0]);
                              jxi = _mm512_loadu_ps(&pjxi[i+0]);
                              cmul_zmm16r4(jxr,jxi,cer,cei,&t0r,&t0i);
                              _mm512_storeu_ps(&fw.pxr[i+0],t0r);
                              _mm512_storeu_ps(&fw.pxi[i+0],t0i);
                              jyr = _mm512_loadu_ps(&pjyr[i+0]);
                              jyi = _mm512_loadu_ps(&pjyi[i+0]);
                              cmul_zmm16r4(jyr,jyi,cer,cei,&t1r,&t1i);
                              _mm512_storeu_ps(&fw.pyr[i+0],t1r);
                              _mm512_storeu_ps(&fw.pyi[i+0],t1i);
                              jzr = _mm512_loadu_ps(&pjzr[i+0]);
                              jzi = _mm512_loadu_ps(&pjzi[i+0]);
                              cmul_zmm16r4(jzr,jzi,cer,cei,&t2r,&t2i);
                              _mm512_storeu_ps(&fw.pzr[i+0],t2r);
                              _mm512_storeu_ps(&fw.pzi[i+0],t2i);
	                 }  
	                 
	                 for(; (i+0) < n; i += 1) {
	                     float cst               = pcst[i];
	                     float rho               = prho[i];
	                     float eai               = k*cst*rho;
	                     std::complex<float> ce = std::cexp({0.0f,eai});
	                     float jxr               = pjxr[i];
	                     float jxi               = pjxi[i];
	                     std::complex<float> cx  = {jxr,jxi};
	                     std::complex<float> jx  = cx*ce;
	                     fw.pxr[i]               = jx.real();
	                     fw.pxi[i]               = jx.imag();
	                     float jyr               = pjyr[i];
	                     float jyi               = pjyi[i];
	                     std::complex<float> cy  = {jyr,jyi};
	                     std::complex<float> jy  = cy*ce;
	                     fw.pyr[i]               = jy.real();
	                     fw.pyi[i]               = jy.imag();
	                     float jzr               = pjzr[i];
	                     float jzi               = pjzi[i];
	                     std::complex<float> cz  = {jzr,jzi};
	                     std::complex<float> jz  = ce*cz;
	                     fw.pzr[i]               = jz.real();
	                     fw.pzi[i]               = jz.imag();
	                 }                             
	       } 
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_avint_u6x(const float * __restrict  pjxr,
	                                            const float * __restrict  pjxi,
	                                            const float * __restrict  pjyr,
	                                            const float * __restrict  pjyi,
	                                            const float * __restrict  pjzr,
	                                            const float * __restrict  pjzi,
	                                            float * __restrict  pxd,
	                                            float * __restrict  pyd,
	                                            float * __restrict  pzd,
	                                            const float * __restrict  prho,
	                                            const float * __restrict  pcst,
	                                            fwork_t fw,
	                                            const float args[7],
	                                            std::complex<float> & Nx,
	                                            std::complex<float> & Ny,
	                                            std::complex<float> & Nz,
	                                            const int32_t n,
	                                            const int32_t PF_DIST
	                                            int32_t & ierr,
	                                            const bool aligned ) {
	                                         
	                float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                            
                        k   = args[0];
                        xa  = args[1];
                        xb  = args[2];
                        ya  = args[3];
                        yb  = args[4];
                        za  = args[5];
                        zb  = args[6];
                        if(aligned) {
                           f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                           f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        sxr = avint(&pxd[0],&fw.pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&fw.pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {goto ERROR;}
                           goto CORRECT;
                        syr = avint(&pyd[0],&fw.pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&fw.pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {goto ERROR;}
                           goto CORRECT;
                        szr = avint(&pzd[0],&fw.pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&fw.pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {goto ERROR;}
                           goto CORRECT;
                        ERROR:
                           {
                               ierr = 3;
                               return;
                        }  
                        CORRECT: {
                           Nx = {sxr,sxi};
                           Ny = {syr,syi};
                           Nz = {szr,szi};  
                        }                                                 
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_avint_u(const float * __restrict pjxr,
	                                         const float * __restrict  pjxi,
	                                         const float * __restrict  pjyr,
	                                         const float * __restrict  pjyi,
	                                         const float * __restrict  pjzr,
	                                         const float * __restrict  pjzi,
	                                         float * __restrict  pxd,
	                                         float * __restrict  pyd,
	                                         float * __restrict  pzd,
	                                         const float * __restrict  prho,
	                                         const float * __restrict  pcost,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz,
	                                         int32_t & ierr ) {
	                                         
	                register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,rho,cost;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;
                        int32_t ier1,ier2,ier3,ier4,ier5,ier6;
                                              
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        cost= _mm512_loadu_ps(&pcost[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        rho = _mm512_loadu_ps(&prho[0]);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        jxr = _mm512_loadu_ps(&pjxr[0]);
                        jxi = _mm512_loadu_ps(&pjxi[0]);
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cost));
                        jyr = _mm512_loadu_ps(&pjyr[0]);
                        jyi = _mm512_loadu_ps(&pjyi[0]);
                        ya  = args[3];
                        yb  = args[4];
                        jzr = _mm512_loadu_ps(&pjzr[0]);
                        jzi = _mm512_loadu_ps(&pjzi[0]);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        sxr = avint(&pxd[0],&pxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&pxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {goto ERROR;}
                           goto CORRECT;
                        syr = avint(&pyd[0],&pyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&pyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {goto ERROR;}
                           goto CORRECT;
                        szr = avint(&pzd[0],&pzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&pzi[0],za,zb,ier6);
                        if(ier5==3 || ier6==3) {goto ERROR;}
                           goto CORRECT;
                        ERROR:
                           {
                               ierr = 3;
                               return;
                        }  
                        CORRECT: {
                           Nx = {sxr,sxi};
                           Ny = {syr,syi};
                           Nz = {szr,szi};  
                        }                                                 
	      }
	      
	      
	       /*
                     Formula 2-23, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'Cubint' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_cubint(const __m512 jxr,
	                                         const __m512 jxi,
	                                         const __m512 jyr,
	                                         const __m512 jyi,
	                                         const __m512 jzr,
	                                         const __m512 jzi,
	                                        __m512 xd,
	                                        __m512 yd,
	                                        __m512 zd,
	                                        const __m512 rho,
	                                        const __m512 cost,
	                                        const float args[7],
	                                        std::complex<float> & Nx,
	                                        std::complex<float> & Ny,
	                                        std::complex<float> & Nz,
	                                        float err[6]) {
	                                        
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;   
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cost));
                        ya  = args[3];
                        yb  = args[4];
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        cubint(16,&pxd[0],&pxr[0],xa,xb,sxr,err[0]);
                        cubint(16,&pxd[0],&pxi[0],xa,xb,sxi,err[1]);
                        cubint(16,&pyd[0],&pyr[0],ya,yb,syr,err[2]);
                        cubint(16,&pyd[0],&pyi[0],ya,yb,syi,err[3]);
                        cubint(16,&pzd[0],&pzr[0],za,zb,szr,err[4]);
                        cubint(16,&pzd[0],&pzi[0],za,zb,szi,err[5]);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                            
	        }
                 
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_cubint_a(const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                         float * __restrict __ATTR_ALIGN__(64) pxd,
	                                         float * __restrict __ATTR_ALIGN__(64) pyd,
	                                         float * __restrict __ATTR_ALIGN__(64) pzd,
	                                         const float * __restrict __ATTR_ALIGN__(64) prho,
	                                         const float * __restrict __ATTR_ALIGN__(64) pcost,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz,
	                                         float err[6]) {
	                                        
	                register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        register __m512 jxr,jxi,jyr;
                        register __m512 jyi,jzr,jzi;
                        register __m512 cst,rho;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;   
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        cst = _mm512_load_ps(&pcost[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        rho = _mm512_load_ps(&prho[0]);
                        xa  = args[1];
                        xb  = args[2];
                        jxr = _mm512_load_ps(&pjxr[0]);
                        jxi = _mm512_load_ps(&pjxi[0]);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        jyr = _mm512_load_ps(&pjyr[0]);
                        jyi = _mm512_load_ps(&pjyi[0]);
                        ya  = args[3];
                        yb  = args[4];
                        jzr = _mm512_load_ps(&pjzr[0]);
                        jzi = _mm512_load_ps(&pjzi[0]);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        cubint(16,&pxd[0],&pxr[0],xa,xb,sxr,err[0]);
                        cubint(16,&pxd[0],&pxi[0],xa,xb,sxi,err[1]);
                        cubint(16,&pyd[0],&pyr[0],ya,yb,syr,err[2]);
                        cubint(16,&pyd[0],&pyi[0],ya,yb,syi,err[3]);
                        cubint(16,&pzd[0],&pzr[0],za,zb,szr,err[4]);
                        cubint(16,&pzd[0],&pzi[0],za,zb,szi,err[5]);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                            
	        }
	        
	        
	        
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_cubint_u6x(const float * __restrict  pjxr,
	                                         const float * __restrict  pjxi,
	                                         const float * __restrict  pjyr,
	                                         const float * __restrict  pjyi,
	                                         const float * __restrict pjzr,
	                                         const float * __restrict  pjzi,
	                                         float * __restrict  pxd,
	                                         float * __restrict  pyd,
	                                         float * __restrict  pzd,
	                                         const float * __restrict  prho,
	                                         const float * __restrict  pcost,
	                                         fwork_t fw,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz,
	                                         const int32_t n,
	                                         const int32_t PF_DIST, 
	                                         float err[6],
	                                         const bool aligned) {
	                                        
	                
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;   
                       
                        k   = args[0];
                        xa  = args[1];
                        xb  = args[2];
                        ya  = args[3];
                        yb  = args[4];
                        za  = args[5];
                        zb  = args[6];
                        if(aligned) {
                            f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                            f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        } 
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        cubint(16,&pxd[0],&fw.pxr[0],xa,xb,sxr,err[0]);
                        cubint(16,&pxd[0],&fw.pxi[0],xa,xb,sxi,err[1]);
                        cubint(16,&pyd[0],&fw.pyr[0],ya,yb,syr,err[2]);
                        cubint(16,&pyd[0],&fw.pyi[0],ya,yb,syi,err[3]);
                        cubint(16,&pzd[0],&fw.pzr[0],za,zb,szr,err[4]);
                        cubint(16,&pzd[0],&fw.pzi[0],za,zb,szi,err[5]);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                            
	        }
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_cubint_u(const float * __restrict pjxr,
	                                         const float * __restrict  pjxi,
	                                         const float * __restrict  pjyr,
	                                         const float * __restrict  pjyi,
	                                         const float * __restrict  pjzr,
	                                         const float * __restrict  pjzi,
	                                         float * __restrict  pxd,
	                                         float * __restrict  pyd,
	                                         float * __restrict  pzd,
	                                         const float * __restrict  prho,
	                                         const float * __restrict  pcost,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz,
	                                         float err[6]) {
	                                        
	                register __m512 intxr,intxi;
                        register __m512 intyr,intyi;
                        register __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        register __m512 jxr,jxi,jyr;
                        register __m512 jyi,jzr,jzi;
                        register __m512 cst,rho;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;   
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        cst = _mm512_loadu_ps(&pcost[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        rho = _mm512_loadu_ps(&prho[0]);
                        xa  = args[1];
                        xb  = args[2];
                        jxr = _mm512_loadu_ps(&pjxr[0]);
                        jxi = _mm512_loadu_ps(&pjxi[0]);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        jyr = _mm512_loadu_ps(&pjyr[0]);
                        jyi = _mm512_loadu_ps(&pjyi[0]);
                        ya  = args[3];
                        yb  = args[4];
                        jzr = _mm512_loadu_ps(&pjzr[0]);
                        jzi = _mm512_loadu_ps(&pjzi[0]);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        cubint(16,&pxd[0],&pxr[0],xa,xb,sxr,err[0]);
                        cubint(16,&pxd[0],&pxi[0],xa,xb,sxi,err[1]);
                        cubint(16,&pyd[0],&pyr[0],ya,yb,syr,err[2]);
                        cubint(16,&pyd[0],&pyi[0],ya,yb,syi,err[3]);
                        cubint(16,&pzd[0],&pzr[0],za,zb,szr,err[4]);
                        cubint(16,&pzd[0],&pzi[0],za,zb,szi,err[5]);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                            
	        }
	        
	        
	         /*
                     Formula 2-23,2-25, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'hiordq' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_hiordq(const __m512 jxr,
	                                         const __m512 jxi,
	                                         const __m512 jyr,
	                                         const __m512 jyi,
	                                         const __m512 jzr,
	                                         const __m512 jzi,
	                                         const __m512 rho,
	                                         const __m512 cst,
	                                         float args[4],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                 
	                __ATTR_ALIGN__(64) float work[32];                        
	                __m512 intxr,intxi;
                        __m512 intyr,intyi;
                        __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;   
                        float k,deltx,delty,deltz;
                        float sxr,sxi,syr,syi,szr,szi;
                        k     = args[0];
                        ir    = _mm512_setzero_ps();
                        ii    = _mm512_set1_ps(1.0f);
                        vk    = _mm512_set1_ps(k);
                        deltx = args[1];
                        ear   = ir;
                        delty = args[2];
                        eai   = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                              _mm512_mul_ps(rho,cst));    
                        deltz = args[3];  
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr   = (float*)&intxr[0];
                        pxi   = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr   = (float*)&intyr[0];
                        pyi   = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);
                        pzr   = (float*)&intzr[0];
                        pzi   = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;    
                        hiordq(16,deltx,&pxr[0],&work[0],sxr);
                        hiordq(16,deltx,&pxi[0],&work[0],sxi);
                        hiordq(16,delty,&pyr[0],&work[0],syr);
                        hiordq(16,delty,&pyi[0],&work[0],syi);
                        hiordq(16,deltz,&pzr[0],&work[0],szr);
                        hiordq(16,deltz,&pzi[0],&work[0],szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                   
	     } 
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_hiordq_a(const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                           const float * __restrict __ATTR_ALIGN__(64) prho,
	                                           const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                           float args[4],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz) {
	                 
	                __ATTR_ALIGN__(64) float work[32];  
	                register jxr,jxi,jyr,jyi;
	                register jzr,jzi,rho,cst;                     
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;   
                        float k,deltx,delty,deltz;
                        float sxr,sxi,syr,syi,szr,szi;
                        jxr   = _mm512_load_ps(&pjxr[0]);
                        jxi   = _mm512_load_ps(&pjxi[0]);
                        jyr   = _mm512_load_ps(&pjyr[0]);
                        jyi   = _mm512_load_ps(&pjyi[0]);
                        jzr   = _mm512_load_ps(&pjzr[0]);
                        jzi   = _mm512_load_ps(&pjzi[0]);
                        rho   = _mm512_load_ps(&prho[0]);
                        cst   = _mm512_load_ps(&pcst[0]);
                        k     = args[0];
                        ir    = _mm512_setzero_ps();
                        ii    = _mm512_set1_ps(1.0f);
                        vk    = _mm512_set1_ps(k);
                        deltx = args[1];
                        ear   = ir;
                        delty = args[2];
                        eai   = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                              _mm512_mul_ps(rho,cst));    
                        deltz = args[3];  
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr   = (float*)&intxr[0];
                        pxi   = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr   = (float*)&intyr[0];
                        pyi   = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);
                        pzr   = (float*)&intzr[0];
                        pzi   = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;    
                        hiordq(16,deltx,&pxr[0],&work[0],sxr);
                        hiordq(16,deltx,&pxi[0],&work[0],sxi);
                        hiordq(16,delty,&pyr[0],&work[0],syr);
                        hiordq(16,delty,&pyi[0],&work[0],syi);
                        hiordq(16,deltz,&pzr[0],&work[0],szr);
                        hiordq(16,deltz,&pzi[0],&work[0],szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                   
	     } 
	     
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_hiordq_u6x(const float * __restrict  pjxr,
	                                           const float * __restrict  pjxi,
	                                           const float * __restrict  pjyr,
	                                           const float * __restrict  pjyi,
	                                           const float * __restrict  pjzr,
	                                           const float * __restrict  pjzi,
	                                           const float * __restrict  prho,
	                                           const float * __restrict  pcst,
	                                           float * __restrict work, // size of work is 2*(n-1)
	                                           float args[4],
	                                           fwork_t fw,
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz,
	                                           const int32_t n,
	                                           const int32_t PF_DIST,
	                                           const bool aligned) {
	                 
	                
                        float k,deltx,delty,deltz;
                        float sxr,sxi,syr,syi,szr,szi;
                        
                        k     = args[0];
                        deltx = args[1];
                        delty = args[2];
                        deltz = args[3];  
                        if(aligned) {
                           f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                            f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;    
                        hiordq(n,deltx,&fw.pxr[0],&work[0],sxr);
                        hiordq(n,deltx,&fw.pxi[0],&work[0],sxi);
                        hiordq(n,delty,&fw.pyr[0],&work[0],syr);
                        hiordq(n,delty,&fw.pyi[0],&work[0],syi);
                        hiordq(n,deltz,&fw.pzr[0],&work[0],szr);
                        hiordq(n,deltz,&fw.pzi[0],&work[0],szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                   
	     } 
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_hiordq_u(const float * __restrict  pjxr,
	                                           const float * __restrict  pjxi,
	                                           const float * __restrict  pjyr,
	                                           const float * __restrict  pjyi,
	                                           const float * __restrict  pjzr,
	                                           const float * __restrict  pjzi,
	                                           const float * __restrict  prho,
	                                           const float * __restrict  pcst,
	                                           float args[4],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz) {
	                 
	                __ATTR_ALIGN__(64) float work[32];  
	                register jxr,jxi,jyr,jyi;
	                register jzr,jzi,rho,cst;                     
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr;   
                        float k,deltx,delty,deltz;
                        float sxr,sxi,syr,syi,szr,szi;
                        jxr   = _mm512_loadu_ps(&pjxr[0]);
                        jxi   = _mm512_loadu_ps(&pjxi[0]);
                        jyr   = _mm512_loadu_ps(&pjyr[0]);
                        jyi   = _mm512_loadu_ps(&pjyi[0]);
                        jzr   = _mm512_loadu_ps(&pjzr[0]);
                        jzi   = _mm512_loadu_ps(&pjzi[0]);
                        rho   = _mm512_loadu_ps(&prho[0]);
                        cst   = _mm512_loadu_ps(&pcst[0]);
                        k     = args[0];
                        ir    = _mm512_setzero_ps();
                        ii    = _mm512_set1_ps(1.0f);
                        vk    = _mm512_set1_ps(k);
                        deltx = args[1];
                        ear   = ir;
                        delty = args[2];
                        eai   = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                              _mm512_mul_ps(rho,cst));    
                        deltz = args[3];  
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr   = (float*)&intxr[0];
                        pxi   = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr   = (float*)&intyr[0];
                        pyi   = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);
                        pzr   = (float*)&intzr[0];
                        pzi   = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;    
                        hiordq(16,deltx,&pxr[0],&work[0],sxr);
                        hiordq(16,deltx,&pxi[0],&work[0],sxi);
                        hiordq(16,delty,&pyr[0],&work[0],syr);
                        hiordq(16,delty,&pyi[0],&work[0],syi);
                        hiordq(16,deltz,&pzr[0],&work[0],szr);
                        hiordq(16,deltz,&pzi[0],&work[0],szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                   
	     } 
	     
	     
                 /*
                     Formula 2-23,2-25, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'plint' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_plint( const __m512 jxr,
	                                         const __m512 jxi,
	                                         const __m512 jyr,
	                                         const __m512 jyi,
	                                         const __m512 jzr,
	                                         const __m512 jzi,
	                                         const __m512 rho,
	                                         const __m512 cst,
	                                         __m512 xd,
	                                         __m512 yd,
	                                         __m512 zd,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;   
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        ya  = args[3];
                        yb  = args[4];
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;   
                        plint(16,&pxd[0],&pxr[0],xa,xb,sxr);
                        plint(16,&pxd[0],&pxi[0],xa,xb,sxi);  
                        plint(16,&pyd[0],&pyr[0],ya,yb,syr);
                        plint(16,&pyd[0],&pyi[0],ya,yb,syi);
                        plint(16,&pzd[0],&pzr[0],za,zb,szr);
                        plint(16,&pzd[0],&pzi[0],za,zb,szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	     }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_plint_a( const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                         const float * __restrict __ATTR_ALIGN__(64) prho,
	                                         const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                         __m512 xd,
	                                         __m512 yd,
	                                         __m512 zd,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,rho,cst;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;  
                        jxr = _mm512_load_ps(&pjxr[0]);
                        jxi = _mm512_load_ps(&pjxi[0]);
                        jyr = _mm512_load_ps(&pjyr[0]);
                        jyi = _mm512_load_ps(&pjyi[0]);
                        jzr = _mm512_load_ps(&pjzr[0]);
                        jzi = _mm512_load_ps(&pjzi[0]);
                        rho = _mm512_load_ps(&prho[0]);
                        cst = _mm512_load_ps(&pcst[0]);
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        ya  = args[3];
                        yb  = args[4];
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;   
                        plint(16,&pxd[0],&pxr[0],xa,xb,sxr);
                        plint(16,&pxd[0],&pxi[0],xa,xb,sxi);  
                        plint(16,&pyd[0],&pyr[0],ya,yb,syr);
                        plint(16,&pyd[0],&pyi[0],ya,yb,syi);
                        plint(16,&pzd[0],&pzr[0],za,zb,szr);
                        plint(16,&pzd[0],&pzi[0],za,zb,szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	     }   
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_plint_u6x( const float * __restrict  pxre,
	                                             const float * __restrict  pxim,
	                                             const float * __restrict  pyre,
	                                             const float * __restrict  pyim,
	                                             const float * __restrict  pzre,
	                                             const float * __restrict  pzim,
	                                             const float * __restrict  prho,
	                                             const float * __restrict  pcst,
	                                             float * __restrict   pxd,
	                                             float * __restrict   pyd,
	                                             float * __restrict   pzd,
	                                             fwork_t fw,
	                                             const float args[7],
	                                             std::complex<float> & Nx,
	                                             std::complex<float> & Ny,
	                                             std::complex<float> & Nz
	                                             const int32_t n,
                                                     const int32_t PF_DIST,
                                                     const bool aligned) {
	                                            
	                                       
	                 
	                float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;   
                        k   = args[0];
                       
                        xa  = args[1];
                        xb  = args[2];
                        ya  = args[3];
                        yb  = args[4];
                        za  = args[5];
                        zb  = args[6];
                        if(aligned) {
                            f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                            f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;   
                        plint(n,&pxd[0],&fw.pxr[0],xa,xb,sxr);
                        plint(n,&pxd[0],&fw.pxi[0],xa,xb,sxi);  
                        plint(n,&pyd[0],&fw.pyr[0],ya,yb,syr);
                        plint(n,&pyd[0],&fw.pyi[0],ya,yb,syi);
                        plint(n,&pzd[0],&fw.pzr[0],za,zb,szr);
                        plint(n,&pzd[0],&fw.pzi[0],za,zb,szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	     }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_plint_u( const float * __restrict  pjxr,
	                                         const float * __restrict  pjxi,
	                                         const float * __restrict  pjyr,
	                                         const float * __restrict  pjyi,
	                                         const float * __restrict  pjzr,
	                                         const float * __restrict  pjzi,
	                                         const float * __restrict  prho,
	                                         const float * __restrict  pcst,
	                                         __m512 xd,
	                                         __m512 yd,
	                                         __m512 zd,
	                                         const float args[7],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,rho,cst;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float k,xa,xb,ya,yb,za,zb;
                        float sxr,sxi,syr,syi,szr,szi;  
                        jxr = _mm512_loadu_ps(&pjxr[0]);
                        jxi = _mm512_loadu_ps(&pjxi[0]);
                        jyr = _mm512_loadu_ps(&pjyr[0]);
                        jyi = _mm512_loadu_ps(&pjyi[0]);
                        jzr = _mm512_loadu_ps(&pjzr[0]);
                        jzi = _mm512_loadu_ps(&pjzi[0]);
                        rho = _mm512_loadu_ps(&prho[0]);
                        cst = _mm512_loadu_ps(&pcst[0]);
                        pxd = (float*)&xd[0];
                        k   = args[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        xa  = args[1];
                        xb  = args[2];
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        ya  = args[3];
                        yb  = args[4];
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        za  = args[5];
                        zb  = args[6];
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;   
                        plint(16,&pxd[0],&pxr[0],xa,xb,sxr);
                        plint(16,&pxd[0],&pxi[0],xa,xb,sxi);  
                        plint(16,&pyd[0],&pyr[0],ya,yb,syr);
                        plint(16,&pyd[0],&pyi[0],ya,yb,syi);
                        plint(16,&pzd[0],&pzr[0],za,zb,szr);
                        plint(16,&pzd[0],&pzi[0],za,zb,szi);
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	     }   
	     
	     
	       /*
                     Formula 2-23,2-25, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'simpne' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpne(const __m512 jxr,
	                                         const __m512 jxi,
	                                         const __m512 jyr,
	                                         const __m512 jyi,
	                                         const __m512 jzr,
	                                         const __m512 jzi,
	                                         const __m512 rho,
	                                         const __m512 cst,
	                                         __m512 xd,
	                                         __m512 yd,
	                                         __m512 zd,
	                                         const float k,
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float sxr,sxi,syr,syi,szr,szi;  
                        pxd = (float*)&xd[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        simpne(16,&pxd[0],&pxr[0],sxr);
                        simpne(16,&pxd[0],&pxi[0],sxi);
                        simpne(16,&pyd[0],&pyr[0],syr);
                        simpne(16,&pyd[0],&pyi[0],syi);
                        simpme(16,&pyd[0],&pzr[0],szr);
                        simpne(16,&pyd[0],&pzi[0],szi);  
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                          
	      }
	     
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpne_a(const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                         const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                         const float * __restrict __ATTR_ALIGN__(64) prho,
	                                         const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                         __m512 xd,
	                                         __m512 yd,
	                                         __m512 zd,
	                                         const float k,
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,cst,rho;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float sxr,sxi,syr,syi,szr,szi;  
                        jxr = _mm512_load_ps(&pjxr[0]);
                        jxi = _mm512_load_ps(&pjxi[0]);
                        jyr = _mm512_load_ps(&pjyr[0]); 
                        jyi = _mm512_load_ps(&pjyi[0]);
                        jzr = _mm512_load_ps(&pjzr[0]);
                        jzi = _mm512_load_ps(&pjzi[0]);
                        cst = _mm512_load_ps(&pcst[0]);
                        rho = _mm512_load_ps(&prho[0]);
                        pxd = (float*)&xd[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        simpne(16,&pxd[0],&pxr[0],sxr);
                        simpne(16,&pxd[0],&pxi[0],sxi);
                        simpne(16,&pyd[0],&pyr[0],syr);
                        simpne(16,&pyd[0],&pyi[0],syi);
                        simpme(16,&pyd[0],&pzr[0],szr);
                        simpne(16,&pyd[0],&pzi[0],szi);  
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                          
	      }
	     
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpne_u6x(const float * __restrict  pjxr,
	                                             const float * __restrict  pjxi,
	                                             const float * __restrict  pjyr,
	                                             const float * __restrict  pjyi,
	                                             const float * __restrict  pjzr,
	                                             const float * __restrict  pjzi,
	                                             const float * __restrict  prho,
	                                             const float * __restrict  pcst,
	                                             float * __restrict  pxd,
	                                             float * __restrict  pyd,
	                                             float * __restrict  pzd,
	                                             fwork_t fw,
	                                             const float k,
	                                             std::complex<float> & Nx,
	                                             std::complex<float> & Ny,
	                                             std::complex<float> & Nz,
	                                             const int32_t n,
	                                             const int32_t PF_DIST,
	                                             const bool aligned) {
	                                         
	             
                        float sxr,sxi,syr,syi,szr,szi;  
                        if(aligned) {
                            f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                            f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        simpne(n,&pxd[0],&fw.pxr[0],sxr);
                        simpne(n,&pxd[0],&fw.pxi[0],sxi);
                        simpne(n,&pyd[0],&fw.pyr[0],syr);
                        simpne(n,&pyd[0],&fw.pyi[0],syi);
                        simpme(n,&pzd[0],&fw.pzr[0],szr);
                        simpne(n,&pzd[0],&fw.pzi[0],szi);  
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                          
	      }
	      
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpne_u(const float * __restrict  pjxr,
	                                         const float * __restrict  pjxi,
	                                         const float * __restrict  pjyr,
	                                         const float * __restrict  pjyi,
	                                         const float * __restrict  pjzr,
	                                         const float * __restrict  pjzi,
	                                         const float * __restrict  prho,
	                                         const float * __restrict  pcst,
	                                         __m512 xd,
	                                         __m512 yd,
	                                         __m512 zd,
	                                         const float k,
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,cst,rho;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float * __restrict pxd = nullptr;
                        float * __restrict pyd = nullptr;
                        float * __restrict pzd = nullptr;
                        float sxr,sxi,syr,syi,szr,szi;  
                        jxr = _mm512_loadu_ps(&pjxr[0]);
                        jxi = _mm512_loadu_ps(&pjxi[0]);
                        jyr = _mm512_loadu_ps(&pjyr[0]); 
                        jyi = _mm512_loadu_ps(&pjyi[0]);
                        jzr = _mm512_loadu_ps(&pjzr[0]);
                        jzi = _mm512_loadu_ps(&pjzi[0]);
                        cst = _mm512_loadu_ps(&pcst[0]);
                        rho = _mm512_loadu_ps(&prho[0]);
                        pxd = (float*)&xd[0];
                        pyd = (float*)&yd[0];
                        vk  = _mm512_set1_ps(k);
                        pzd = (float*)&zd[0];
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        simpne(16,&pxd[0],&pxr[0],sxr);
                        simpne(16,&pxd[0],&pxi[0],sxi);
                        simpne(16,&pyd[0],&pyr[0],syr);
                        simpne(16,&pyd[0],&pyi[0],syi);
                        simpme(16,&pyd[0],&pzr[0],szr);
                        simpne(16,&pyd[0],&pzi[0],szi);  
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                          
	      }
	      
	      
	       /*
                     Formula 2-23,2-25, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'simpn' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpn( const __m512 jxr,
	                                         const __m512 jxi,
	                                         const __m512 jyr,
	                                         const __m512 jyi,
	                                         const __m512 jzr,
	                                         const __m512 jzi,
	                                         const __m512 rho,
	                                         const __m512 cst,
	                                         const float args[2],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        vk  = _mm512_set1_ps(k);
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        simpn(16,h,&pxr[0],sxr);
                        simpn(16,h,&pxi[0],sxi);
                        simpn(16,h,&pyr[0],syr);
                        simpn(16,h,&pyi[0],syi);
                        simpn(16,h,&pzr[0],szr);
                        simpn(16,h,&pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	     
	      
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpn_a( const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                           const float * __restrict __ATTR_ALIGN__(64) prho,
	                                           const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                           const float args[2],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,cst,rho;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        jxr = _mm512_load_ps(&pjxr[0]);
                        vk  = _mm512_set1_ps(k);
                        jxi = _mm512_load_ps(&pjxi[0]);
                        ir  = _mm512_setzero_ps();
                        jyr = _mm512_load_ps(&pjyr[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        jyi = _mm512_load_ps(&pjyi[0]);
                        ear = ir;
                        jzr = _mm512_load_ps(&pjzr[0]);
                        jzi = _mm512_load_ps(&pjzi[0]);
                        cst = _mm512_load_ps(&pcst[0]);
                        rho = _mm512_load_ps(&prho[0]);
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        simpn(16,h,&pxr[0],sxr);
                        simpn(16,h,&pxi[0],sxi);
                        simpn(16,h,&pyr[0],syr);
                        simpn(16,h,&pyi[0],syi);
                        simpn(16,h,&pzr[0],szr);
                        simpn(16,h,&pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpn_u6x( const float * __restrict  pjxr,
	                                           const float * __restrict  pjxi,
	                                           const float * __restrict  pjyr,
	                                           const float * __restrict  pjyi,
	                                           const float * __restrict  pjzr,
	                                           const float * __restrict  pjzi,
	                                           const float * __restrict  prho,
	                                           const float * __restrict  pcst,
	                                           fwork_t fw,
	                                           const float args[2],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz,
	                                           const int32_t n,
	                                           const int32_t PF_DIST,
	                                           const bool aligned) {
	                                         
	                
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        if(aligned) {
                            f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                            f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        simpn(n,h,&fw.pxr[0],sxr);
                        simpn(n,h,&fw.pxi[0],sxi);
                        simpn(n,h,&fw.pyr[0],syr);
                        simpn(n,h,&fw.pyi[0],syi);
                        simpn(n,h,&fw.pzr[0],szr);
                        simpn(n,h,&fw.pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_simpn_u( const float * __restrict  pjxr,
	                                           const float * __restrict  pjxi,
	                                           const float * __restrict  pjyr,
	                                           const float * __restrict  pjyi,
	                                           const float * __restrict  pjzr,
	                                           const float * __restrict  pjzi,
	                                           const float * __restrict  prho,
	                                           const float * __restrict  pcst,
	                                           const float args[2],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,cst,rho;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        jxr = _mm512_loadu_ps(&pjxr[0]);
                        vk  = _mm512_set1_ps(k);
                        jxi = _mm512_loadu_ps(&pjxi[0]);
                        ir  = _mm512_setzero_ps();
                        jyr = _mm512_loadu_ps(&pjyr[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        jyi = _mm512_loadu_ps(&pjyi[0]);
                        ear = ir;
                        jzr = _mm512_loadu_ps(&pjzr[0]);
                        jzi = _mm512_loadu_ps(&pjzi[0]);
                        cst = _mm512_loadu_ps(&pcst[0]);
                        rho = _mm512_loadu_ps(&prho[0]);
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        simpn(16,h,&pxr[0],sxr);
                        simpn(16,h,&pxi[0],sxi);
                        simpn(16,h,&pyr[0],syr);
                        simpn(16,h,&pyi[0],syi);
                        simpn(16,h,&pzr[0],szr);
                        simpn(16,h,&pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	        /*
                     Formula 2-23,2-25, p. 36
                     Both electric and magnetic quantities (field amplitudes) are computed
                     by the single kernel (different surface currents shall be passed only).
                     'wedint' integrator.
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_wedint( const __m512 jxr,
	                                         const __m512 jxi,
	                                         const __m512 jyr,
	                                         const __m512 jyi,
	                                         const __m512 jzr,
	                                         const __m512 jzi,
	                                         const __m512 rho,
	                                         const __m512 cst,
	                                         const float args[2],
	                                         std::complex<float> & Nx,
	                                         std::complex<float> & Ny,
	                                         std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        vk  = _mm512_set1_ps(k);
                        ir  = _mm512_setzero_ps();
                        ii  = _mm512_set1_ps(1.0f);
                        ear = ir;
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        wedint(16,h,&pxr[0],sxr);
                        wedint(16,h,&pxi[0],sxi);
                        wedint(16,h,&pyr[0],syr);
                        wedint(16,h,&pyi[0],syi);
                        wedint(16,h,&pzr[0],szr);
                        wedint(16,h,&pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_wedint_a( const float * __restrict __ATTR_ALIGN__(64) pjxr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjxi,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjyr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjyi,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjzr,
	                                           const float * __restrict __ATTR_ALIGN__(64) pjzi,
	                                           const float * __restrict __ATTR_ALIGN__(64) prho,
	                                           const float * __restrict __ATTR_ALIGN__(64) pcst,
	                                           const float args[2],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,cst,rho;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        jxr = _mm512_load_ps(&pjxr[0]);
                        vk  = _mm512_set1_ps(k);
                        jxi = _mm512_load_ps(&pjxi[0]);
                        ir  = _mm512_setzero_ps();
                        jyr = _mm512_load_ps(&pjyr[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        jyi = _mm512_load_ps(&pjyi[0]);
                        ear = ir;
                        jzr = _mm512_load_ps(&pjzr[0]);
                        jzi = _mm512_load_ps(&pjzi[0]);
                        cst = _mm512_load_ps(&pcst[0]);
                        rho = _mm512_load_ps(&prho[0]);
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        wedint(16,h,&pxr[0],sxr);
                        wedint(16,h,&pxi[0],sxi);
                        wedint(16,h,&pyr[0],syr);
                        wedint(16,h,&pyi[0],syi);
                        wedint(16,h,&pzr[0],szr);
                        wedint(16,h,&pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   static inline
	           void Nem_f2235_zmm16r4_wedint_u6x( const float * __restrict  pjxr,
	                                           const float * __restrict  pjxi,
	                                           const float * __restrict  pjyr,
	                                           const float * __restrict  pjyi,
	                                           const float * __restrict  pjzr,
	                                           const float * __restrict  pjzi,
	                                           const float * __restrict  prho,
	                                           const float * __restrict  pcst,
	                                           fwork_t fw,
	                                           const float args[2],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz,
	                                           const int32_t n,
	                                           const int32_t PF_DIST,
	                                           const bool aligned) {
	                                         
	                
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        if(aligned) {
                            f2235_integrand_zmm16r4_u6x_a(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }
                        else {
                            f2235_integrand_zmm16r4_u6x_u(pjxr,pjxi,pjyr,pjyi,
                                                         pjzr,pjzi,prho,pcst,
                                                         fw,k,n,PF_DIST);
                        }  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        wedint(n,h,&fw.pxr[0],sxr);
                        wedint(n,h,&fw.pxi[0],sxi);
                        wedint(n,h,&fw.pyr[0],syr);
                        wedint(n,h,&fw.pyi[0],syi);
                        wedint(n,h,&fw.pzr[0],szr);
                        wedint(n,h,&fw.pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void Nem_f2235_zmm16r4_wedint_u( const float * __restrict  pjxr,
	                                           const float * __restrict  pjxi,
	                                           const float * __restrict  pjyr,
	                                           const float * __restrict  pjyi,
	                                           const float * __restrict  pjzr,
	                                           const float * __restrict  pjzi,
	                                           const float * __restrict  prho,
	                                           const float * __restrict  pcst,
	                                           const float args[2],
	                                           std::complex<float> & Nx,
	                                           std::complex<float> & Ny,
	                                           std::complex<float> & Nz) {
	                                         
	                 __m512 intxr,intxi;
                         __m512 intyr,intyi;
                         __m512 intzr,intzi;
                        register __m512 jxr,jxi,jyr,jyi;
                        register __m512 jzr,jzi,cst,rho;
                        register __m512 vk,ii,ir,ear,eai;
                        register __m512 cer,cei,t0r,t0i;
                        float * __restrict pxr = nullptr;
                        float * __restrict pxi = nullptr;
                        float * __restrict pyr = nullptr;
                        float * __restrict pyi = nullptr;
                        float * __restrict pzr = nullptr;
                        float * __restrict pzi = nullptr; 
                        float sxr,sxi,syr,syi,szr,szi;
                        float k,h;
                        k = args[0];
                        h = args[1];
                        jxr = _mm512_loadu_ps(&pjxr[0]);
                        vk  = _mm512_set1_ps(k);
                        jxi = _mm512_loadu_ps(&pjxi[0]);
                        ir  = _mm512_setzero_ps();
                        jyr = _mm512_loadu_ps(&pjyr[0]);
                        ii  = _mm512_set1_ps(1.0f);
                        jyi = _mm512_loadu_ps(&pjyi[0]);
                        ear = ir;
                        jzr = _mm512_loadu_ps(&pjzr[0]);
                        jzi = _mm512_loadu_ps(&pjzi[0]);
                        cst = _mm512_loadu_ps(&pcst[0]);
                        rho = _mm512_loadu_ps(&prho[0]);
                        eai = _mm512_mul_ps(_mm512_mul_ps(ii,vk),
                                            _mm512_mul_ps(rho,cst));
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cmul_zmm16r4(jxr,jxi,cer,cei,&intxr,&intxi);
                        pxr = (float*)&intxr[0];
                        pxi = (float*)&intxi[0];
                        cmul_zmm16r4(jyr,jyi,cer,cei,&intyr,&intyi);
                        pyr = (float*)&intyr[0];
                        pyi = (float*)&intyi[0];
                        cmul_zmm16r4(jzr,jzi,cer,cei,&intzr,&intzi);  
                        pzr = (float*)&intzr[0];
                        pzi = (float*)&intzi[0];  
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;  
                        wedint(16,h,&pxr[0],sxr);
                        wedint(16,h,&pxi[0],sxi);
                        wedint(16,h,&pyr[0],syr);
                        wedint(16,h,&pyi[0],syi);
                        wedint(16,h,&pzr[0],szr);
                        wedint(16,h,&pzi[0],szi);
                        
                        Nx = {sxr,sxi};
                        Ny = {syr,syi};
                        Nz = {szr,szi};                                
	       }
	       
	       
	       /*
	           Formula 2-26, p. 37
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f226_zmm16r4(const __m512 tht1,
	                               const __m512 tht2,
	                               const __m512 phi1,
	                               const __m512 phi2) {
	                               
	                register __m512 ctht1,ctht2,stht1,stht2;
	                register __m512 cdif,cpsi;
	                dif   = _mm512_sub_ps(phi1,phi2);
	                stht1 = xsinf(tht1);
	                ctht1 = xcosf(tht1);
	                stht2 = xsinf(tht2);
	                cdif  = xcosf(dif);
	                ctht2 = xcosf(tht2);
	                cpsi  = _mm512_fmadd_ps(ctht1,ctht2,
	                                    _mm512_mul_ps(stht1,
	                                             _mm512_mul_ps(stht2,cdif)));
	                return (cpsi);   
	        }
	        
	        
	        
	          __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f226_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ptht1,
	                                 const float * __restrict __ATTR_ALIGN__(64) ptht2,
	                                 const float * __restrict __ATTR_ALIGN__(64) pphi1,
	                                 const float * __restrict __ATTR_ALIGN__(64) pphi2) {
	                 
	                register __m512 tht1,tht2,phi1,phi2;
	                register __m512 ctht1,ctht2,stht1,stht2;
	                register __m512 cdif,cpsi;
	                tht1  = _mm512_load_ps(&ptht1[0]);
	                tht2  = _mm512_load_ps(&ptht2[0]);
	                phi1  = _mm512_load_ps(&pphi1[0]);
	                phi2  = _mm512_load_ps(&pphi2[0]);
	                dif   = _mm512_sub_ps(phi1,phi2);
	                stht1 = xsinf(tht1);
	                ctht1 = xcosf(tht1);
	                stht2 = xsinf(tht2);
	                cdif  = xcosf(dif);
	                ctht2 = xcosf(tht2);
	                cpsi  = _mm512_fmadd_ps(ctht1,ctht2,
	                                    _mm512_mul_ps(stht1,
	                                             _mm512_mul_ps(stht2,cdif)));
	                return (cpsi);   
	        }
	        
	     
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 f226_zmm16r4_u(const float * __restrict  ptht1,
	                                 const float * __restrict  ptht2,
	                                 const float * __restrict  pphi1,
	                                 const float * __restrict  pphi2) {
	                 
	                register __m512 tht1,tht2,phi1,phi2;
	                register __m512 ctht1,ctht2,stht1,stht2;
	                register __m512 cdif,cpsi;
	                tht1  = _mm512_loadu_ps(&ptht1[0]);
	                tht2  = _mm512_loadu_ps(&ptht2[0]);
	                phi1  = _mm512_loadu_ps(&pphi1[0]);
	                phi2  = _mm512_loadu_ps(&pphi2[0]);
	                dif   = _mm512_sub_ps(phi1,phi2);
	                stht1 = xsinf(tht1);
	                ctht1 = xcosf(tht1);
	                stht2 = xsinf(tht2);
	                cdif  = xcosf(dif);
	                ctht2 = xcosf(tht2);
	                cpsi  = _mm512_fmadd_ps(ctht1,ctht2,
	                                    _mm512_mul_ps(stht1,
	                                             _mm512_mul_ps(stht2,cdif)));
	                return (cpsi);   
	        }
	        
	        
	        /*
	            Formula: 2-36, p. 39
	            
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void N_f236_zmm16r4(const __m512 nxr,
	                               const __m512 nxi,
	                               const __m512 nyr,
	                               const __m512 nyi,
	                               const __m512 nzr,
	                               const __m512 nzi,
	                               const __m512 phi,
	                               const __m512 tht,
	                               __m512 * __restrict Nthr,
	                               __m512 * __restrict Nthi,
	                               __m512 * __restrict Nphr,
	                               __m512 * __restrict Nphi) {
	                using namespace gms::math;              
	                register __m512 tr,ti,t0r,t0i,t1r,t1i;
	                register __m512 cphi,sphi,ctht,stht;
	                cphi   = xcosf(phi);
	                stht   = xsinf(tht);
	                ctht   = xcosf(tht);
	                sphi   = xsinf(tht);    
	                tr     = _mm512_fmadd_ps(nxr,cphi,
	                                   _mm512_mul_ps(nyr,sphi));
	                t0r    = _mm512_fmsub_ps(tr,ctht,
	                                   _mm512_mul_ps(nzr,stht));
	                *Nthr  = t0r;
	                ti     = _mm512_fmadd_ps(nxi,cphi,
	                                   _mm512_mul_ps(nyi,sphi)); 
	                t0i    = _mm512_fmsub_ps(ti,ctht,
	                                   _mm512_mul_ps(nzi,stht));   
	                *Nthi  = t0i;
	                tr     = negate_zmm16r4(nxr);
	                t1r    = _mm512_fmadd_ps(tr,sphi,
	                                     _mm512_mul_ps(nyr,cphi));
	                *Nphr  = t1r;
	                ti     = negate_zmm16r4(nxi);
	                t1i    = _mm512_fmadd_ps(ti,sphi,
	                                     _mm512_mul_ps(nyi,cphi));
	                Nphi   = t1i;
	       }
	     
	      
	     
                 
               
        } // radiolocation

} // gms












#endif /*__GMS_ANTENNA_FEEDER_ZMM16R4_HPP__*/
