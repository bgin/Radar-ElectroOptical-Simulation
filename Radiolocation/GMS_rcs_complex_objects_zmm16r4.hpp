

#ifndef __GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__
#define __GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__

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

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_CREATION_DATE = "11-05-2023 10:53 PM +00200 (THR 11 MAY 2023 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm16r4.hpp"
#include "GMS_simd_utils.hpp"




namespace  gms {


         namespace radiolocation {
         
              
              /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Formula: 8.1-21
              */     
              
              
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8121_zmm16r4(const __m512 gam,
                                             const __m512 phi,
                                             const __m512 k0,
                                             __m512 * __restrict D1r,
                                             __m512 * __restrict D1i,
                                            __m512 * __restrict D2r,
                                            __m512 * __restrict D2i) {
                                            
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);
                        register __m512 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm512_setzero_ps();
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = xsinf(x1);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cos1 = xcosf(x1);
                        x0   = _mm512_mul_ps(invn,phi2);
                        cos2 = xcosf(x0);
                        cer  = _mm512_mul_ps(cer,sinp);
                        trm1 = _mm512_rcp14_ps(_mm512_sub_ps(cos1),C10);
                        cei  = _mm512_mul_ps(cei,sinp);
                        trm2 = _mm512_rcp14_ps(_mm512_sub_ps(cos1,cos2));
                        sqr  = _mm512_mul_ps(invn,sqr);
                        ear  = _mm512_mul_ps(cer,sqr);
                        eai  = _mm512_mul_ps(cei,sqr);
                        x0   = _mm512_sub_ps(trm1,trm2);
                        *D1r = _mm512_mul_ps(ear,x0);
                        *D1i = _mm512_mul_ps(eai,x0);
                        x1   = _mm512_add_ps(trm1,trm2);
                        *D2r = _mm512_mul_ps(ear,x1);
                        *D2i = _mm512_mul_ps(eai,x1);
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8121_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam,
                                               const float * __restrict __ATTR_ALIGN__(64) pphi,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0,
                                               float * __restrict __ATTR_ALIGN__(64) D1r,
                                               float * __restrict __ATTR_ALIGN__(64) D1i,
                                               float * __restrict __ATTR_ALIGN__(64) D2r,
                                               float * __restrict __ATTR_ALIGN__(64) D2i) {
                                   
                        register __m512 gam = _mm512_load_ps(&pgam[0]);
                        register __m512 phi = _mm512_load_ps(&pphi[0]);  
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);       
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);
                        register __m512 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm512_setzero_ps();
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = xsinf(x1);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cos1 = xcosf(x1);
                        x0   = _mm512_mul_ps(invn,phi2);
                        cos2 = xcosf(x0);
                        cer  = _mm512_mul_ps(cer,sinp);
                        trm1 = _mm512_rcp14_ps(_mm512_sub_ps(cos1),C10);
                        cei  = _mm512_mul_ps(cei,sinp);
                        trm2 = _mm512_rcp14_ps(_mm512_sub_ps(cos1,cos2));
                        sqr  = _mm512_mul_ps(invn,sqr);
                        ear  = _mm512_mul_ps(cer,sqr);
                        eai  = _mm512_mul_ps(cei,sqr);
                        x0   = _mm512_sub_ps(trm1,trm2);
                        _mm512_store_ps(&D1r[0] ,_mm512_mul_ps(ear,x0));
                        _mm512_store_ps(&D1i[0] ,_mm512_mul_ps(eai,x0));
                        x1   = _mm512_add_ps(trm1,trm2);
                        _mm512_store_ps(&D2r[0] ,_mm512_mul_ps(ear,x1));
                        _mm512_store_ps(&D2i[0] ,_mm512_mul_ps(eai,x1));
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8121_zmm16r4_u(const float * __restrict  pgam,
                                               const float * __restrict  pphi,
                                               const float * __restrict  pk0,
                                               float * __restrict  D1r,
                                               float * __restrict  D1i,
                                               float * __restrict  D2r,
                                               float * __restrict  D2i) {
                                   
                        register __m512 gam = _mm512_loadu_ps(&pgam[0]);
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]);  
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);       
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);
                        register __m512 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm512_setzero_ps();
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = xsinf(x1);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        cos1 = xcosf(x1);
                        x0   = _mm512_mul_ps(invn,phi2);
                        cos2 = xcosf(x0);
                        cer  = _mm512_mul_ps(cer,sinp);
                        trm1 = _mm512_rcp14_ps(_mm512_sub_ps(cos1),C10);
                        cei  = _mm512_mul_ps(cei,sinp);
                        trm2 = _mm512_rcp14_ps(_mm512_sub_ps(cos1,cos2));
                        sqr  = _mm512_mul_ps(invn,sqr);
                        ear  = _mm512_mul_ps(cer,sqr);
                        eai  = _mm512_mul_ps(cei,sqr);
                        x0   = _mm512_sub_ps(trm1,trm2);
                        _mm512_storeu_ps(&D1r[0] ,_mm512_mul_ps(ear,x0));
                        _mm512_storeu_ps(&D1i[0] ,_mm512_mul_ps(eai,x0));
                        x1   = _mm512_add_ps(trm1,trm2);
                        _mm512_storeu_ps(&D2r[0] ,_mm512_mul_ps(ear,x1));
                        _mm512_storeu_ps(&D2i[0] ,_mm512_mul_ps(eai,x1));
                }
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter singly diffracted far-zone fields (E,H).
                    Formula: 8.1-19, 8.1-20
                
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EsHs_f811920_zmm16r4(    const __m512 betai,
                                                 const __m512 betas,
                                                 const __m512 gam,
                                                 const __m512 phi,
                                                 const __m512 k0,
                                                 const __m512 r,
                                                 const __m512 rho,
                                                 const __m512 psi,
                                                 __m512 * __restrict Esr,
                                                 __m512 * __restrict Esi,
                                                 __m512 * __restrict Hsr,
                                                 __m512 * __restrict Hsi) {
                                                 
                       register __m512 ear,eai,cer,cei;
                       register __m512 D1r,D1i,D2r,D2i,x0,x1;
                       register __m512 rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm512_fmadd_ps(k0,r,psi);
                       coef_D12_f8121_zmm16r4(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm512_rcp14_ps(r);
                       ear   = _mm512_setzero_ps();
                       cosbi = xcosf(betai);
                       eai   = k0rp;
                       cosbs = xcosf(betas);
                       cexp_zmm16r4(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_ps(cer,invr);    
                       rhos  = _mm512_div_ps(rho,_mm512_add_ps(cosbi,cosbs));
                       cei   = _mm512_mul_ps(cei,invr);
                       sqrho = _mm512_sqrt_ps(rhos);
                       x0    = _mm512_mul_ps(sqrho,cer);
                       x1    = _mm512_mul_ps(sqrho,cei);
                       *Esr  = _mm512_mul_ps(D1r,x0);
                       *Hsr  = _mm512_mul_ps(D2r,x0);
                       *Esi  = _mm512_mul_ps(D1i,x1);
                       *Hsi  = _mm512_mul_ps(D2i,x1);                               
            }
            
            
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EsHs_f811920_zmm16r4_a(  const float * __restrict __ATTR_ALIGN__(64) pbetai,
                                                 const float * __restrict __ATTR_ALIGN__(64) pbetas,
                                                 const float * __restrict __ATTR_ALIGN__(64) pgam,
                                                 const float * __restrict __ATTR_ALIGN__(64) pphi,
                                                 const float * __restrict __ATTR_ALIGN__(64) pk0,
                                                 const float * __restrict __ATTR_ALIGN__(64) pr,
                                                 const float * __restrict __ATTR_ALIGN__(64) prho,
                                                 const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                                 float * __restrict __ATTR_ALIGN__(64) Esr,
                                                 float * __restrict __ATTR_ALIGN__(64) Esi,
                                                 float * __restrict __ATTR_ALIGN__(64) Hsr,
                                                 float * __restrict __ATTR_ALIGN__(64) Hsi) {
                              
                       register __m512 betai = _mm512_load_ps(&pbetai[0]);
                       register __m512 betas = _mm512_load_ps(&pbetas[0]); 
                       register __m512 gam   = _mm512_load_ps(&pgam[0]);   
                       register __m512 k0    = _mm512_load_ps(&pk0[0]); 
                       register __m512 r     = _mm512_load_ps(&pr[0]);
                       register __m512 rho   = _mm512_load_ps(&prho[0]); 
                       register __m512 psi   = _mm512_load_ps(&ppsi[0]);             
                       register __m512 ear,eai,cer,cei;
                       register __m512 D1r,D1i,D2r,D2i,x0,x1;
                       register __m512 rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm512_fmadd_ps(k0,r,psi);
                       coef_D12_f8121_zmm16r4(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm512_rcp14_ps(r);
                       ear   = _mm512_setzero_ps();
                       cosbi = xcosf(betai);
                       eai   = k0rp;
                       cosbs = xcosf(betas);
                       cexp_zmm16r4(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_ps(cer,invr);    
                       rhos  = _mm512_div_ps(rho,_mm512_add_ps(cosbi,cosbs));
                       cei   = _mm512_mul_ps(cei,invr);
                       sqrho = _mm512_sqrt_ps(rhos);
                       x0    = _mm512_mul_ps(sqrho,cer);
                       x1    = _mm512_mul_ps(sqrho,cei);
                       _mm512_store_ps(&Esr[0] ,_mm512_mul_ps(D1r,x0));
                       _mm512_store_ps(&Hsr[0] ,_mm512_mul_ps(D2r,x0));
                       _mm512_store_ps(&Esi[0] ,_mm512_mul_ps(D1i,x1));
                       _mm512_store_ps(&Hsi[0] ,_mm512_mul_ps(D2i,x1));                               
            }
            
            
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EsHs_f811920_zmm16r4_u(  const float * __restrict  pbetai,
                                                 const float * __restrict  pbetas,
                                                 const float * __restrict  pgam,
                                                 const float * __restrict  pphi,
                                                 const float * __restrict  pk0,
                                                 const float * __restrict  pr,
                                                 const float * __restrict  prho,
                                                 const float * __restrict ppsi,
                                                 float * __restrict  Esr,
                                                 float * __restrict  Esi,
                                                 float * __restrict  Hsr,
                                                 float * __restrict  Hsi) {
                              
                       register __m512 betai = _mm512_loadu_ps(&pbetai[0]);
                       register __m512 betas = _mm512_loadu_ps(&pbetas[0]); 
                       register __m512 gam   = _mm512_loadu_ps(&pgam[0]);   
                       register __m512 k0    = _mm512_loadu_ps(&pk0[0]); 
                       register __m512 r     = _mm512_loadu_ps(&pr[0]);
                       register __m512 rho   = _mm512_loadu_ps(&prho[0]); 
                       register __m512 psi   = _mm512_loadu_ps(&ppsi[0]);             
                       register __m512 ear,eai,cer,cei;
                       register __m512 D1r,D1i,D2r,D2i,x0,x1;
                       register __m512 rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm512_fmadd_ps(k0,r,psi);
                       coef_D12_f8121_zmm16r4(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm512_rcp14_ps(r);
                       ear   = _mm512_setzero_ps();
                       cosbi = xcosf(betai);
                       eai   = k0rp;
                       cosbs = xcosf(betas);
                       cexp_zmm16r4(ear,eai,&cer,&cei);
                       cer   = _mm512_mul_ps(cer,invr);    
                       rhos  = _mm512_div_ps(rho,_mm512_add_ps(cosbi,cosbs));
                       cei   = _mm512_mul_ps(cei,invr);
                       sqrho = _mm512_sqrt_ps(rhos);
                       x0    = _mm512_mul_ps(sqrho,cer);
                       x1    = _mm512_mul_ps(sqrho,cei);
                       _mm512_storeu_ps(&Esr[0] ,_mm512_mul_ps(D1r,x0));
                       _mm512_storeu_ps(&Hsr[0] ,_mm512_mul_ps(D2r,x0));
                       _mm512_storeu_ps(&Esi[0] ,_mm512_mul_ps(D1i,x1));
                       _mm512_storeu_ps(&Hsi[0] ,_mm512_mul_ps(D2i,x1));                               
            }
            
            
            /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Ray normal-incidence to one of edge faces.
                    Formula: 8.1-24
            */
            
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8124_zmm16r4(const __m512 k0,
                                               const __m512 gam,
                                               __m512 * __restrict D1r,
                                               __m512 * __restrict D1i,
                                               __m512 * __restrict D2r,
                                               __m512 * __restrict D2i) {
                                               
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);  
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm512_setzero_ps();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328);  
                        x1   = _mm512_mul_ps(k0,C6283185307179586476925286766559);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        invn = _mm512_rcp14_ps(x0);
                        cer  = gms::math::negate_zmm16r4(cer);
                        invn2= _mm512_add_ps(invn,invn);      
                        cei  = gms::math::negate_zmm16r4(cei);
                        arg1 = _mm512_mul_ps(C314159265358979323846264338328,invn2)
                        sin1 = xsinf(arg1);
                        arg2 = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        t0r  = _mm512_div_ps(cer,x1);
                        cos1 = xcosf(arg1);
                        x0   = _mm512_div_ps(cos1,sin1);
                        t0i  = _mm512_div_ps(cei,x1);   
                        sin2 = xsinf(arg2);
                        cos2 = xcosf(arg2);   
                        x1   = _mm512_div_ps(cos2,sin2);
                        ear  = _mm512_fmsub_ps(invn,x0,_mm512_mul_ps(invn2,x1));
                        eai  = _mm512_fmadd_ps(invn,x0,_mm512_mul_ps(invn2,x1));
                        cer  = _mm512_mul_ps(t0r,ear);
                        cei  = _mm512_mul_ps(t0i,eai);
                        *D1r = cer;
                        *D2r = cer;
                        *D1i = cei;
                        *D2i = cei;       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8124_zmm16r4_a(const  float * __restrict __ATTR_ALIGN__(64) pk0,
                                                 const  float * __restrict __ATTR_ALIGN__(64) pgam,
                                                 float * __restrict __ATTR_ALIGN__(64) D1r,
                                                 float * __restrict __ATTR_ALIGN__(64) D1i,
                                                 float * __restrict __ATTR_ALIGN__(64) D2r,
                                                 float * __restrict __ATTR_ALIGN__(64) D2i) {
                                    
                        register __m512 k0 = _mm512_load_ps(&pk0[0]);
                        register __m512 gam= _mm512_load_ps(&pgam[0]);           
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);  
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm512_setzero_ps();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328);  
                        x1   = _mm512_mul_ps(k0,C6283185307179586476925286766559);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        invn = _mm512_rcp14_ps(x0);
                        cer  = gms::math::negate_zmm16r4(cer);
                        invn2= _mm512_add_ps(invn,invn);      
                        cei  = gms::math::negate_zmm16r4(cei);
                        arg1 = _mm512_mul_ps(C314159265358979323846264338328,invn2)
                        sin1 = xsinf(arg1);
                        arg2 = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        t0r  = _mm512_div_ps(cer,x1);
                        cos1 = xcosf(arg1);
                        x0   = _mm512_div_ps(cos1,sin1);
                        t0i  = _mm512_div_ps(cei,x1);   
                        sin2 = xsinf(arg2);
                        cos2 = xcosf(arg2);   
                        x1   = _mm512_div_ps(cos2,sin2);
                        ear  = _mm512_fmsub_ps(invn,x0,_mm512_mul_ps(invn2,x1));
                        eai  = _mm512_fmadd_ps(invn,x0,_mm512_mul_ps(invn2,x1));
                        cer  = _mm512_mul_ps(t0r,ear);
                        cei  = _mm512_mul_ps(t0i,eai);
                        _mm512_store_ps(&D1r[0] ,cer);
                        _mm512_store_ps(&D2r[0] ,cer);
                        _mm512_store_ps(&D1i[0] ,cei);
                        _mm512_store_ps(&D2i[0] ,cei);       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8124_zmm16r4_u(const  float * __restrict  pk0,
                                                 const  float * __restrict pgam,
                                                 float * __restrict  D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict  D2r,
                                                 float * __restrict  D2i) {
                                    
                        register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                        register __m512 gam= _mm512_loadu_ps(&pgam[0]);           
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);
                        const __m512 C10                              = _mm512_set1_ps(1.0f);  
                        register __m512 ear,eai,cer,cei,t0r,t0i;
                        register __m512 x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm512_setzero_ps();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328);  
                        x1   = _mm512_mul_ps(k0,C6283185307179586476925286766559);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        invn = _mm512_rcp14_ps(x0);
                        cer  = gms::math::negate_zmm16r4(cer);
                        invn2= _mm512_add_ps(invn,invn);      
                        cei  = gms::math::negate_zmm16r4(cei);
                        arg1 = _mm512_mul_ps(C314159265358979323846264338328,invn2)
                        sin1 = xsinf(arg1);
                        arg2 = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        t0r  = _mm512_div_ps(cer,x1);
                        cos1 = xcosf(arg1);
                        x0   = _mm512_div_ps(cos1,sin1);
                        t0i  = _mm512_div_ps(cei,x1);   
                        sin2 = xsinf(arg2);
                        cos2 = xcosf(arg2);   
                        x1   = _mm512_div_ps(cos2,sin2);
                        ear  = _mm512_fmsub_ps(invn,x0,_mm512_mul_ps(invn2,x1));
                        eai  = _mm512_fmadd_ps(invn,x0,_mm512_mul_ps(invn2,x1));
                        cer  = _mm512_mul_ps(t0r,ear);
                        cei  = _mm512_mul_ps(t0i,eai);
                        _mm512_storeu_ps(&D1r[0] ,cer);
                        _mm512_storeu_ps(&D2r[0] ,cer);
                        _mm512_storeu_ps(&D1i[0] ,cei);
                        _mm512_storeu_ps(&D2i[0] ,cei);       
                }
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Backscatter direction axial caustic (for slightly diffracted rays).
                    Formula: 8.1-26
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_Ddiff_f8126_zmm16r4(const __m512 gam,
                                             const __m512 phi,
                                             const __m512 k0,
                                             __m512 * __restrict Dr,
                                             __m512 * __restrict Di) {
                                             
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);  
                        const __m512 C20                              = 
                                                     _mm512_set1_ps(2.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_mul_ps(C6283185307179586476925286766559,k0);
                        ear  = _mm512_setzero_ps();
                        n    = _mm512_div_ps(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm512_rcp14_ps(n); 
                        sqr  = _mm512_sqrt_ps(x0);
                        pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        spin = xsinf(pin);
                        cer  = _mm512_div_ps(cer,sqr);
                        cpin = xcosf(pin);
                        cei  = _mm512_div_ps(cei,sqr);
                        cphin= xcosf(_mm512_mul_ps(phi2,invn));
                        x0   = _mm512_mul_ps(_mm512_mul_ps(C20,invn),spin);
                        x1   = _mm512_sub_ps(cpin,cphin);
                        n    = _mm512_div_ps(x0,x1);
                        *Dr  = _mm512_mul_ps(cer,n);
                        *Di  = _mm512_mul_ps(cei,n);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_Ddiff_f8126_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pgam,
                                               const float * __restrict __ATTR_ALIGN__(64) pphi,
                                               const float * __restrict __ATTR_ALIGN__(64) pk0,
                                               float * __restrict __ATTR_ALIGN__(64)  Dr,
                                               float * __restrict __ATTR_ALIGN__(64)  Di) {
                                    
                        register __m512 gam  = _mm512_load_ps(&pgam[0]);
                        register __m512 phi  = _mm512_load_ps(&pphi[0]);
                        register __m512 k0   = _mm512_load_ps(&pk0[0]);         
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);  
                        const __m512 C20                              = 
                                                     _mm512_set1_ps(2.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_mul_ps(C6283185307179586476925286766559,k0);
                        ear  = _mm512_setzero_ps();
                        n    = _mm512_div_ps(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm512_rcp14_ps(n); 
                        sqr  = _mm512_sqrt_ps(x0);
                        pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        spin = xsinf(pin);
                        cer  = _mm512_div_ps(cer,sqr);
                        cpin = xcosf(pin);
                        cei  = _mm512_div_ps(cei,sqr);
                        cphin= xcosf(_mm512_mul_ps(phi2,invn));
                        x0   = _mm512_mul_ps(_mm512_mul_ps(C20,invn),spin);
                        x1   = _mm512_sub_ps(cpin,cphin);
                        n    = _mm512_div_ps(x0,x1);
                        _mm512_store_ps(&Dr[0] ,_mm512_mul_ps(cer,n));
                        _mm512_store_ps(&Di[0] ,_mm512_mul_ps(cei,n));
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_Ddiff_f8126_zmm16r4_u(const float * __restrict  pgam,
                                                   const float * __restrict pphi,
                                                   const float * __restrict  pk0,
                                                   float * __restrict  Dr,
                                                   float * __restrict  Di) {
                                    
                        register __m512 gam  = _mm512_loadu_ps(&pgam[0]);
                        register __m512 phi  = _mm512_loadu_ps(&pphi[0]);
                        register __m512 k0   = _mm512_loadu_ps(&pk0[0]);         
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f);  
                        const __m512 C20                              = 
                                                     _mm512_set1_ps(2.0f);
                        register __m512 ear,eai,cer,cei;
                        register __m512 invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm512_add_ps(phi,phi);
                        x0   = _mm512_mul_ps(C6283185307179586476925286766559,k0);
                        ear  = _mm512_setzero_ps();
                        n    = _mm512_div_ps(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm512_rcp14_ps(n); 
                        sqr  = _mm512_sqrt_ps(x0);
                        pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        spin = xsinf(pin);
                        cer  = _mm512_div_ps(cer,sqr);
                        cpin = xcosf(pin);
                        cei  = _mm512_div_ps(cei,sqr);
                        cphin= xcosf(_mm512_mul_ps(phi2,invn));
                        x0   = _mm512_mul_ps(_mm512_mul_ps(C20,invn),spin);
                        x1   = _mm512_sub_ps(cpin,cphin);
                        n    = _mm512_div_ps(x0,x1);
                        _mm512_storeu_ps(&Dr[0] ,_mm512_mul_ps(cer,n));
                        _mm512_storeu_ps(&Di[0] ,_mm512_mul_ps(cei,n));
                }
                
                
                   /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Backscatter direction axial caustic (for slightly diffracted rays).
                    Scattered Electric and Magnetic fields.
                    Formula: 8.1-25
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EsHs_f8125_zmm16r4(const __m512 a,
                                           const __m512 k0,
                                           const __m512 r,
                                           const __m512 gam,
                                           const __m512 phi,
                                           const __m512 psi,
                                           __m512 * __restrict Esr,
                                           __m512 * __restrict Esi,
                                           __m512 * __restrict Hsr,
                                           __m512 * __restrict Hsi) {
                                           
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f); 
                        const __m512 C05                              = 
                                                     _mm512_set1_ps(0.5f);
                        register __m512 Dr,Di,x0,x1;
                        register __m512 ear,eai,cer,cei;
                        register __m512 k0r,invr,t0r,t0i;
                        register __m512 sqr;
                        k0r  = _mm512_mul_ps(k0,r);
                        ear  = _mm512_setzero_ps();
                        eai  = _mm512_add_ps(k0r,_mm512_sub_ps(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_zmm16r4(gam,phi,k0,Dr,Di); 
                        invr = _mm512_rcp_ps(r);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(C6283185307179586476925286766559,k0));
                        t0r  = _mm512_mul_ps(Dr,C05);
                        cer  = _mm512_mul_ps(cer,invr);
                        sqr  = _mm512_mul_ps(a,sqr);
                        t0i  = _mm512_mul_ps(Di,C05);  
                        cei  = _mm512_mul_ps(cei,invr);
                        x0   = _mm512_mul_ps(t0r,_mm512_mul_ps(sqr,cer));
                        x1   = _mm512_mul_ps(t0i,_mm512_mul_ps(sqr,cei));
                        *Esr = gms::math::negate_zmm16r4(x0);
                        *Esi = gms::math::negate_zmm16r4(x1); 
                        *Hsr = x0;
                        *Hsi = x1;             
                                    
                }
                
            
                
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__*/
