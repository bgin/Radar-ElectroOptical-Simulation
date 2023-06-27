

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
          
           
                   /*
                       Spherical unit vectors.
                   */
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void spher_unitv_zmm16r4(const __m512 ex,
	                                   const __m512 ey,
	                                   const __m512 ez,
	                                   const __m512 tht,
	                                   const __m512 phi,
	                                   __m512 * __restrict er,
	                                   __m512 * __restrict eth,
	                                   __m512 * __restrict eph) {
	               
	               using namespace gms::math;                   
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 t0,t1;
	               stht = xsinf(tht);
	               sphi = xsinf(phi);
	               t0   = _mm512_mul_ps(sphi,ey);
	               cphi = xcosf(phi);
	               t1   = _mm512_mul_ps(cphi,ex);
	               ctht = xcosf(tht);
	               *er  = _mm512_fmadd_ps(stht,t1,
	                                 _mm512_fmadd_ps(stht,t0,
	                                             _mm512_mul_ps(ctht,ez)));    
	               *eth = _mm512_fmadd_ps(ctht,t1,
	                                 _mm512_fmsub_ps(ctht,t0,
	                                             _mm512_mul_ps(stht,ez)));  
	               *eph = _mm512_fmadd_ps(negate_zmm16r4(sphi),ex,
	                                             _mm512_mul_ps(cphi,ey));                           
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void spher_unitv_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pex,
	                                      const float * __restrict __ATTR_ALIGN__(64) pey,
	                                      const float * __restrict __ATTR_ALIGN__(64) pez,
	                                      const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                      const float * __restrict __ATTR_ALIGN__(64) pphi,
	                                      float * __restrict __ATTR_ALIGN__(64)  er,
	                                      float * __restrict __ATTR_ALIGN__(64)  eth,
	                                      float * __restrict __ATTR_ALIGN__(64)  eph) {
	               
	               using namespace gms::math;     
	               register __m512 ex = _mm512_load_ps(&pex[0]);
	               register __m512 ey = _mm512_load_ps(&pey[0]);
	               register __m512 ez = _mm512_load_ps(&pez[0]);              
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 t0,t1;
	               stht = xsinf(tht);
	               sphi = xsinf(phi);
	               t0   = _mm512_mul_ps(sphi,ey);
	               cphi = xcosf(phi);
	               t1   = _mm512_mul_ps(cphi,ex);
	               ctht = xcosf(tht);
	               _mm512_store_ps(&er[0]  ,_mm512_fmadd_ps(stht,t1,
	                                                  _mm512_fmadd_ps(stht,t0,
	                                                       _mm512_mul_ps(ctht,ez))));    
	               _mm512_store_ps(&eth[0] ,_mm512_fmadd_ps(ctht,t1,
	                                                  _mm512_fmsub_ps(ctht,t0,
	                                                       _mm512_mul_ps(stht,ez))));  
	               _mm512_store_ps(&eph[0] ,_mm512_fmadd_ps(negate_zmm16r4(sphi),ex,
	                                                       _mm512_mul_ps(cphi,ey)));                           
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void spher_unitv_zmm16r4_u(const float * __restrict  pex,
	                                      const float * __restrict  pey,
	                                      const float * __restrict  pez,
	                                      const float * __restrict  ptht,
	                                      const float * __restrict  pphi,
	                                      float * __restrict  er,
	                                      float * __restrict  eth,
	                                      float * __restrict  eph) {
	               
	               using namespace gms::math;     
	               register __m512 ex = _mm512_load_ps(&pex[0]);
	               register __m512 ey = _mm512_load_ps(&pey[0]);
	               register __m512 ez = _mm512_load_ps(&pez[0]);              
	               register __m512 stht,sphi,cphi,ctht;
	               register __m512 t0,t1;
	               stht = xsinf(tht);
	               sphi = xsinf(phi);
	               t0   = _mm512_mul_ps(sphi,ey);
	               cphi = xcosf(phi);
	               t1   = _mm512_mul_ps(cphi,ex);
	               ctht = xcosf(tht);
	               _mm512_store_ps(&er[0]  ,_mm512_fmadd_ps(stht,t1,
	                                                  _mm512_fmadd_ps(stht,t0,
	                                                       _mm512_mul_ps(ctht,ez))));    
	               _mm512_store_ps(&eth[0] ,_mm512_fmadd_ps(ctht,t1,
	                                                  _mm512_fmsub_ps(ctht,t0,
	                                                       _mm512_mul_ps(stht,ez))));  
	               _mm512_store_ps(&eph[0] ,_mm512_fmadd_ps(negate_zmm16r4(sphi),ex,
	                                                       _mm512_mul_ps(cphi,ey)));                           
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
	           void f13_zmm16r4(            const __m512 eth,
	                                        const __m512 eph,
	                                        const __m512 nthr,
	                                        const __m512 nthi,
	                                        const __m512 nphr,
	                                        const __m512 nphi,
	                                        __m512 * __restrict Nthr,
	                                        __m512 * __restrict Nthi,
	                                        __m512 * __restrict Nphr,
	                                        __m512 * __restrict Nphi) {
	                                      
	                *Nthr = _mm512_mul_ps(nthr,eth);
	                *Nphr = _mm512_mul_ps(nphr,eph);
	                *Nthi = _mm512_mul_ps(nthi,eth);
	                *Nphi = _mm512_mul_ps(nphi,eph);                             
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void f13_zmm16r4_a( const float * __restrict __ATTR_ALIGN__(64) peth,
	                               const float * __restrict __ATTR_ALIGN__(64) peph,
	                               const float * __restrict __ATTR_ALIGN__(64) pnthr,
	                               const float * __restrict __ATTR_ALIGN__(64) pnthi,
	                               const float * __restrict __ATTR_ALIGN__(64) pnphr,
	                               const float * __restrict __ATTR_ALIGN__(64) pnphi,
	                               float * __restrict __ATTR_ALIGN__(64) Nthr,
	                               float * __restrict __ATTR_ALIGN__(64) Nthi,
	                               float * __restrict __ATTR_ALIGN__(64) Nphr,
	                               float * __restrict __ATTR_ALIGN__(64) Nphi) {
	                              
	                 register __m512 eth = _mm512_load_ps(&peth[0]);
	                 register __m512 eph = _mm512_load_ps(&peph[0]); 
	                 register __m512 nthr= _mm512_load_ps(&pnthr[0]);  
	                 register __m512 nthi= _mm512_load_ps(&pnthi[0]);  
	                 register __m512 nphr= _mm512_load_ps(&pnphr[0]);  
	                 register __m512 nphi= _mm512_load_ps(&pnphi[0]); 
	                 _mm512_store_ps(&Nthr[0] ,_mm512_mul_ps(nthr,eth));
	                 _mm512_store_ps(&Nphr[0] ,_mm512_mul_ps(nphr,eph));
	                 _mm512_store_ps(&Nthi[0] ,_mm512_mul_ps(nthi,eth));
	                 _mm512_store_ps(&Nphi[0] ,_mm512_mul_ps(nphi,eph));                             
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void f13_zmm16r4_u( const float * __restrict  peth,
	                               const float * __restrict  peph,
	                               const float * __restrict  pnthr,
	                               const float * __restrict  pnthi,
	                               const float * __restrict  pnphr,
	                               const float * __restrict  pnphi,
	                               float * __restrict  Nthr,
	                               float * __restrict  Nthi,
	                               float * __restrict  Nphr,
	                               float * __restrict  Nphi) {
	                              
	                 register __m512 eth = _mm512_loadu_ps(&peth[0]);
	                 register __m512 eph = _mm512_loadu_ps(&peph[0]); 
	                 register __m512 nthr= _mm512_loadu_ps(&pnthr[0]);  
	                 register __m512 nthi= _mm512_loadu_ps(&pnthi[0]);  
	                 register __m512 nphr= _mm512_loadu_ps(&pnphr[0]);  
	                 register __m512 nphi= _mm512_loadu_ps(&pnphi[0]); 
	                 _mm512_storeu_ps(&Nthr[0] ,_mm512_mul_ps(nthr,eth));
	                 _mm512_storeu_ps(&Nphr[0] ,_mm512_mul_ps(nphr,eph));
	                 _mm512_storeu_ps(&Nthi[0] ,_mm512_mul_ps(nthi,eth));
	                 _mm512_storeu_ps(&Nphi[0] ,_mm512_mul_ps(nphi,eph));                             
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
	            Hertz vector (electrical), avint integrator.
	            Formula 2-13, p. 35
	       */
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hve_f213_zmm16r4_avint(const __m512 xre,
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
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xre,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yre,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zre,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xim,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yim,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zim,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&dx[0],&intxr[0],xa,xb,ier1);
                        sxi = avint(&dx[0],&intxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&dy[0],&intyr[0],ya,yb,ier3);
                        syi = avint(&dy[0],&intyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&dz[0],&intzr[0],za,zb,ier5);
                        szi = avint(&dz[0],&intzi[0],za,zb,ier6);
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
	           void hve_f213_zmm16r4_avint_a(const float * __restrict __ATTR_ALIGN__(64) pxre,
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
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                       
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xre,xim,yre,yim,zre,zim;
                        register __m512 xd,yd,zd;
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xre,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yre,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zre,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xim,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yim,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zim,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&pxd[0],&intxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&intxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&pyd[0],&intyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&intyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&pzd[0],&intzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&intzi[0],za,zb,ier6);
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
	           void hve_f213_zmm16r4_avint_u(const float * __restrict  pxre,
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
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                       
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xre,xim,yre,yim,zre,zim;
                        register __m512 xd,yd,zd;
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
                        _mm512_storeu_ps(&intxr[0],_mm512_mul_ps(xre,cer));
                        _mm512_storeu_ps(&intyr[0],_mm512_mul_ps(yre,cer));
                        _mm512_storeu_ps(&intzr[0],_mm512_mul_ps(zre,cer));
                        _mm512_storeu_ps(&intxi[0],_mm512_mul_ps(xim,cei));
                        _mm512_storeu_ps(&intyi[0],_mm512_mul_ps(yim,cei));
                        _mm512_storeu_ps(&intzi[0],_mm512_mul_ps(zim,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        sxr = avint(&pxd[0],&intxr[0],xa,xb,ier1);
                        sxi = avint(&pxd[0],&intxi[0],xa,xb,ier2);
                        if(ier1==3 || ier2==3) {
                           ierr = 3;
                           return;
                        }  
                        syr = avint(&pyd[0],&intyr[0],ya,yb,ier3);
                        syi = avint(&pyd[0],&intyi[0],ya,yb,ier4);
                        if(ier3==3 || ier4==3) {
                           ierr = 3;
                           return;
                        }  
                        szr = avint(&pzd[0],&intzr[0],za,zb,ier5);
                        szi = avint(&pzd[0],&intzi[0],za,zb,ier6);
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
	           void hve_f213_zmm16r4_cubint(const __m512 xre,
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
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xre,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yre,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zre,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xim,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yim,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zim,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(ntab,&dx[0],&intxr[0],xa,xb,sxr,er1);
                        cubint(ntab,&dx[0],&intxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(ntab,&dy[0],&intyr[0],ya,yb,syr,er3);
                        cubint(ntab,&dy[0],&intyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(ntab,&dz[0],&intzr[0],za,zb,szr,er5);
                        cubint(ntab,&dz[0],&intzi[0],za,zb,szi,er6);
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
	           void hve_f213_zmm16r4_cubint_a(const float * __restrict __ATTR_ALIGN__(64)  pxre,
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
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        register __m512 xr,xi,yr,yi,zr,zi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        register float er1,er2,er3,er4,er5,er6;
                        xr = _mm512_load_ps(&xre[0]);
                        xi = _mm512_load_ps(&xim[0]);
                        yr = _mm512_load_ps(&yre[0]);
                        yi = _mm512_load_ps(&yim[0]);
                        zr = _mm512_load_ps(&zre[0]);
                        zi = _mm512_load_ps(&zim[0]);
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xr,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yr,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zr,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xi,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yi,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zi,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(ntab,&pxd[0],&intxr[0],xa,xb,sxr,er1);
                        cubint(ntab,&pxd[0],&intxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(ntab,&pyd[0],&intyr[0],ya,yb,syr,er3);
                        cubint(ntab,&pyd[0],&intyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(ntab,&pzd[0],&intzr[0],za,zb,szr,er5);
                        cubint(ntab,&pzd[0],&intzi[0],za,zb,szi,er6);
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
	           void hve_f213_zmm16r4_cubint_u(const float * __restrict  pxre,
	                                          const float * __restrict  pxim,
	                                          const float * __restrict  pyre,
	                                          const float * __restrict  pyim,
	                                          const float * __restrict  pzre,
	                                          const float * __restrict  pzim,
	                                          float * __restrict   pxd,
	                                          float * __restrict   pyd,
	                                          float * __restrict   pzd,
	                                          const float arg[10],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz,
                                                  float err[6]) { // integration error.
                            
                        
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi   
                        constexpr int32_t ntab = 16; 
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        register __m512 xr,xi,yr,yi,zr,zi;
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register float k,r,xa,xb,ya,yb,za,zb;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        register float er1,er2,er3,er4,er5,er6;
                        xr = _mm512_loadu_ps(&xre[0]);
                        xi = _mm512_loadu_ps(&xim[0]);
                        yr = _mm512_loadu_ps(&yre[0]);
                        yi = _mm512_loadu_ps(&yim[0]);
                        zr = _mm512_loadu_ps(&zre[0]);
                        zi = _mm512_loadu_ps(&zim[0]);
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
                        _mm512_storeu_ps(&intxr[0],_mm512_mul_ps(xr,cer));
                        _mm512_storeu_ps(&intyr[0],_mm512_mul_ps(yr,cer));
                        _mm512_storeu_ps(&intzr[0],_mm512_mul_ps(zr,cer));
                        _mm512_storeu_ps(&intxi[0],_mm512_mul_ps(xi,cei));
                        _mm512_storeu_ps(&intyi[0],_mm512_mul_ps(yi,cei));
                        _mm512_storeu_ps(&intzi[0],_mm512_mul_ps(zi,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        cubint(ntab,&pxd[0],&intxr[0],xa,xb,sxr,er1);
                        cubint(ntab,&pxd[0],&intxi[0],xa,xb,sxi,er2);
                        err[0] = er1;
                        err[1] = er2;
                        cubint(ntab,&pyd[0],&intyr[0],ya,yb,syr,er3);
                        cubint(ntab,&pyd[0],&intyi[0],ya,yb,syi,er4);
                        err[2] = er3;
                        err[3] = er4;
                        cubint(ntab,&pzd[0],&intzr[0],za,zb,szr,er5);
                        cubint(ntab,&pzd[0],&intzi[0],za,zb,szi,er6);
                        err[4] = er5;
                        err[5] = er6;
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                /*
	            Hertz vector (electrical), hiordq integrator.
	            Formula 2-13, p. 35
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hve_f213_zmm16r4_hiordq(const __m512 xre,
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
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        __ATTR_ALIGN__(64) float work[32];
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xre,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yre,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zre,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xim,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yim,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zim,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&intxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&intxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&intyr[0],&work[0],syr);
                        hiordq(ntab,delty,&intyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&intzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&intzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hve_f213_zmm16r4_hiordq_a(const float * __restrict __ATTR_ALIGN__(64)  pxre,
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
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        __ATTR_ALIGN__(64) float work[32];
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        register float k,r,deltx,delty,deltz;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_load_ps(&xre[0]);
                        xi = _mm512_load_ps(&xim[0]);
                        yr = _mm512_load_ps(&yre[0]);
                        yi = _mm512_load_ps(&yim[0]);
                        zr = _mm512_load_ps(&zre[0]);
                        zi = _mm512_load_ps(&zim[0]);
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xr,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yr,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zr,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xi,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yi,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zi,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&intxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&intxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&intyr[0],&work[0],syr);
                        hiordq(ntab,delty,&intyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&intzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&intzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hve_f213_zmm16r4_hiordq_u(const float * __restrict   pxre,
	                                          const float * __restrict   pxim,
	                                          const float * __restrict   pyre,
	                                          const float * __restrict   pyim,
	                                          const float * __restrict   pzre,
	                                          const float * __restrict   pzim,
	                                          const float arg[7],
	                                          std::complex<float> & hx,                        
                                                  std::complex<float> & hy,
                                                  std::complex<float> & hz) {
                                               
                                              
                        constexpr float C12566370614359172953850573533118 = 
                                              12.566370614359172953850573533118f; //4*pi 
                        constexpr int32_t ntab = 16;   
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        __ATTR_ALIGN__(64) float work[32];
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
                        register __m512 xr,xi,yr,yi,zr,zi;
                        register float k,r,deltx,delty,deltz;
                        register float omg,eps,sxr,sxi,syr,syi,szr,szi,frac;
                        xr = _mm512_loadu_ps(&xre[0]);
                        xi = _mm512_loadu_ps(&xim[0]);
                        yr = _mm512_loadu_ps(&yre[0]);
                        yi = _mm512_loadu_ps(&yim[0]);
                        zr = _mm512_loadu_ps(&zre[0]);
                        zi = _mm512_loadu_ps(&zim[0]);
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
                        _mm512_storeu_ps(&intxr[0],_mm512_mul_ps(xr,cer));
                        _mm512_storeu_ps(&intyr[0],_mm512_mul_ps(yr,cer));
                        _mm512_storeu_ps(&intzr[0],_mm512_mul_ps(zr,cer));
                        _mm512_storeu_ps(&intxi[0],_mm512_mul_ps(xi,cei));
                        _mm512_storeu_ps(&intyi[0],_mm512_mul_ps(yi,cei));
                        _mm512_storeu_ps(&intzi[0],_mm512_mul_ps(zi,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        hiordq(ntab,deltx,&intxr[0],&work[0],sxr);
                        hiordq(ntab,deltx,&intxi[0],&work[0],sxi);
                        hiordq(ntab,delty,&intyr[0],&work[0],syr);
                        hiordq(ntab,delty,&intyi[0],&work[0],syi);
                        hiordq(ntab,deltz,&intzr[0],&work[0],szr);
                        hiordq(ntab,deltz,&intzi[0],&work[0],szi);
                      
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
               
                /*
	            Hertz vector (electrical), plint integrator.
	            Formula 2-13, p. 35
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           void hve_f213_zmm16r4_plint(const __m512 xre,
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
                        const float NAN =  std::numeric_limits<float>::quiet_NaN();           
                        __ATTR_ALIGN__(64) float intxr[16];
                        __ATTR_ALIGN__(64) float intxi[16];
                        __ATTR_ALIGN__(64) float intyr[16];
                        __ATTR_ALIGN__(64) float intyi[16];
                        __ATTR_ALIGN__(64) float intzr[16];
                        __ATTR_ALIGN__(64) float intzi[16];
                        __ATTR_ALIGN__(64) float dx[16];
                        __ATTR_ALIGN__(64) float dy[16];
                        __ATTR_ALIGN__(64) float dz[16];
                        register __m512 vk,vr,ii,ir,invr,cer,cei,eai;
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
                        _mm512_store_ps(&intxr[0],_mm512_mul_ps(xre,cer));
                        _mm512_store_ps(&intyr[0],_mm512_mul_ps(yre,cer));
                        _mm512_store_ps(&intzr[0],_mm512_mul_ps(zre,cer));
                        _mm512_store_ps(&intxi[0],_mm512_mul_ps(xim,cei));
                        _mm512_store_ps(&intyi[0],_mm512_mul_ps(yim,cei));
                        _mm512_store_ps(&intzi[0],_mm512_mul_ps(zim,cei));
                        sxr = 0.0f;
                        sxi = sxr;
                        syi = sxr;
                        syr = sxr;
                        szr = sxr;
                        szi = sxr;
                        float tmp = C12566370614359172953850573533118*omg*eps;
                        frac = 1.0f/tmp;
                        plint(ntab,&dx[0],&intxr[0],xa,xb,sxr);
                        plint(ntab,&dx[0],&intxi[0],xa,xb,sxi);
                        plint(ntab,&dy[0],&intyr[0],ya,yb,syr);
                        plint(ntab,&dy[0],&intyi[0],ya,yb,syi);
                        plint(ntab,&dz[0],&intzr[0],za,zb,szr);
                        plint(ntab,&dz[0],&intzi[0],za,zb,szi);
                        hx = {sxr*frac,sxi*frac};
                        hy = {syr*frac,syi*frac};
                        hz = {szr*frac,szi*frac};                 
               }
               
	       
               
               
               
               
               
        } // radiolocation

} // gms












#endif /*__GMS_ANTENNA_FEEDER_ZMM16R4_HPP__*/
