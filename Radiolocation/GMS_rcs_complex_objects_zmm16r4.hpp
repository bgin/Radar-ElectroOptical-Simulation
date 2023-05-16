

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
                   Work (input) arrays for kernel rcs_f8162_zmm16r4_2t_u and
                   rcs_f8162_zmm16r4_2t_a.
               */
               __ATTR_ALIGN__(64) struct RCS_F8162_DATA {
               
                       float * __restrict  Ya1; 
                       float * __restrict  Ya2; 
                       float * __restrict  Ya3; 
                       float * __restrict  Ea;  
                       float * __restrict  WRKa; 
                       float * __restrict  Yb1;
                       float * __restrict  Yb2; 
                       float * __restrict  Yb3; 
                       float * __restrict  Eb; 
                       float * __restrict  WRKb;  
               };
         
              
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
                
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EsHs_f8125_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pa,
                                             const float * __restrict __ATTR_ALIGN__(64) pk0,
                                             const float * __restrict __ATTR_ALIGN__(64) pr,
                                             const float * __restrict __ATTR_ALIGN__(64) pgam,
                                             const float * __restrict __ATTR_ALIGN__(64) pphi,
                                             const float * __restrict __ATTR_ALIGN__(64) ppsi,
                                             float * __restrict __ATTR_ALIGN__(64) Esr,
                                             float * __restrict __ATTR_ALIGN__(64) Esi,
                                             float * __restrict __ATTR_ALIGN__(64) Hsr,
                                             float * __restrict __ATTR_ALIGN__(64) Hsi) {
                            
                        register __m512 a   = _mm512_load_ps(&pa[0]);
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);   
                        register __m512 r   = _mm512_load_ps(&pr[0]);   
                        register __m512 gam = _mm512_load_ps(&pgam[0]);   
                        register __m512 phi = _mm512_load_ps(&pphi[0]); 
                        register __m512 psi = _mm512_load_ps(&ppsi[0]);     
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
                        _mm512_store_ps(&Esr[0] ,gms::math::negate_zmm16r4(x0));
                        _mm512_store_ps(&Esi[0] ,gms::math::negate_zmm16r4(x1)); 
                        _mm512_store_ps(&Hsr[0] ,x0);
                        _mm512_store_ps(&Hsi[0] ,x1);             
                      
                }
                
            
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void EsHs_f8125_zmm16r4_u(const float * __restrict  pa,
                                             const float * __restrict  pk0,
                                             const float * __restrict  pr,
                                             const float * __restrict pgam,
                                             const float * __restrict  pphi,
                                             const float * __restrict  ppsi,
                                             float * __restrict  Esr,
                                             float * __restrict  Esi,
                                             float * __restrict  Hsr,
                                             float * __restrict  Hsi) {
                            
                        register __m512 a   = _mm512_loadu_ps(&pa[0]);
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);   
                        register __m512 r   = _mm512_loadu_ps(&pr[0]);   
                        register __m512 gam = _mm512_loadu_ps(&pgam[0]);   
                        register __m512 phi = _mm512_loadu_ps(&pphi[0]); 
                        register __m512 psi = _mm512_loadu_ps(&ppsi[0]);     
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
                        _mm512_storeu_ps(&Esr[0] ,gms::math::negate_zmm16r4(x0));
                        _mm512_storeu_ps(&Esi[0] ,gms::math::negate_zmm16r4(x1)); 
                        _mm512_storeu_ps(&Hsr[0] ,x0);
                        _mm512_storeu_ps(&Hsi[0] ,x1);             
                      
                }
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Doubly and high order diffracted rays --
                    bistatic diffraction coefficients.
                    Formula: 8.1-27  
                
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8127_zmm16r4(const __m512 k0,
                                               const __m512 gam,
                                               const __m512 phi1,
                                               const __m512 phi2,
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
                        register __m512 invn,x0,x1,ear,eai,cer,cei;
                        register __m512 sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328);
                        phid = _mm512_mul_ps(_mm512_sub_ps(phi1,phi2),invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm512_setzero_ps();
                        phis = _mm512_mul_ps(_mm512_add_ps(phi1,phi2),invn);
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        x0   = xsinf(pin);
                        cer  = _mm512_div_ps(cer,sqr);
                        spin = _mm512_mul_ps(x0,invn);
                        cei  = _mm512_div_ps(cei,sqr);
                        cpin = xcosf(pin);
                        x0   = xcosf(phid);
                        invc1= _mm512_rcp14_ps(_mm512_sub_ps(cpin,x0));
                        x1   = xcosf(phis);
                        invc2= _mm512_rcp14_ps(_mm512_sub_ps(cpin,x1));
                        ear  = _mm512_mul_ps(cer,spin);
                        phis = _mm512_sub_ps(invc1,invc2);
                        eai  = _mm512_mul_ps(cei,spin);
                        phid = _mm512_add_ps(invc1,invc2);
                        *D1r = _mm512_mul_ps(ear,phis);
                        *D1i = _mm512_mul_ps(eai,phis);
                        *D2r = _mm512_mul_ps(ear,phid);
                        *D2i = _mm512_mul_ps(eai,phid);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8127_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                                 const float * __restrict __ATTR_ALIGN__(64) pgam,
                                                 const float * __restrict __ATTR_ALIGN__(64) pphi1,
                                                 const float * __restrict __ATTR_ALIGN__(64) pphi2,
                                                 float * __restrict __ATTR_ALIGN__(64)  D1r,
                                                 float * __restrict __ATTR_ALIGN__(64)  D1i,
                                                 float * __restrict __ATTR_ALIGN__(64)  D2r,
                                                 float * __restrict __ATTR_ALIGN__(64)  D2i) {
                              
                        register __m512 k0  = _mm512_load_ps(&pk0[0]);
                        register __m512 gam = _mm512_load_ps(&pgam[0]);
                        register __m512 phi1= _mm512_load_ps(&pphi1[0]);
                        register __m512 phi2= _mm512_load_ps(&pphi2[0]);                 
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                        register __m512 invn,x0,x1,ear,eai,cer,cei;
                        register __m512 sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328);
                        phid = _mm512_mul_ps(_mm512_sub_ps(phi1,phi2),invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm512_setzero_ps();
                        phis = _mm512_mul_ps(_mm512_add_ps(phi1,phi2),invn);
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        x0   = xsinf(pin);
                        cer  = _mm512_div_ps(cer,sqr);
                        spin = _mm512_mul_ps(x0,invn);
                        cei  = _mm512_div_ps(cei,sqr);
                        cpin = xcosf(pin);
                        x0   = xcosf(phid);
                        invc1= _mm512_rcp14_ps(_mm512_sub_ps(cpin,x0));
                        x1   = xcosf(phis);
                        invc2= _mm512_rcp14_ps(_mm512_sub_ps(cpin,x1));
                        ear  = _mm512_mul_ps(cer,spin);
                        phis = _mm512_sub_ps(invc1,invc2);
                        eai  = _mm512_mul_ps(cei,spin);
                        phid = _mm512_add_ps(invc1,invc2);
                        _mm512_store_ps(&D1r[0] ,_mm512_mul_ps(ear,phis));
                        _mm512_store_ps(&D1i[0] ,_mm512_mul_ps(eai,phis));
                        _mm512_store_ps(&D2r[0] ,_mm512_mul_ps(ear,phid));
                        _mm512_store_ps(&D2i[0] ,_mm512_mul_ps(eai,phid));
                }
                
                
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   void coef_D12_f8127_zmm16r4_u(const float * __restrict  pk0,
                                                 const float * __restrict  pgam,
                                                 const float * __restrict  pphi1,
                                                 const float * __restrict  pphi2,
                                                 float * __restrict   D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict   D2r,
                                                 float * __restrict  D2i) {
                              
                        register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                        register __m512 gam = _mm512_loadu_ps(&pgam[0]);
                        register __m512 phi1= _mm512_loadu_ps(&pphi1[0]);
                        register __m512 phi2= _mm512_loadu_ps(&pphi2[0]);                 
                        const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                        const __m512 C6283185307179586476925286766559 = 
                                                     _mm512_set1_ps(6.283185307179586476925286766559f);
                        const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                        register __m512 invn,x0,x1,ear,eai,cer,cei;
                        register __m512 sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm512_div_ps(gam,C314159265358979323846264338328);
                        phid = _mm512_mul_ps(_mm512_sub_ps(phi1,phi2),invn);
                        sqr  = _mm512_sqrt_ps(_mm512_mul_ps(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm512_setzero_ps();
                        phis = _mm512_mul_ps(_mm512_add_ps(phi1,phi2),invn);
                        invn = _mm512_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                        cexp_zmm16r4(ear,eai,&cer,&cei);
                        x0   = xsinf(pin);
                        cer  = _mm512_div_ps(cer,sqr);
                        spin = _mm512_mul_ps(x0,invn);
                        cei  = _mm512_div_ps(cei,sqr);
                        cpin = xcosf(pin);
                        x0   = xcosf(phid);
                        invc1= _mm512_rcp14_ps(_mm512_sub_ps(cpin,x0));
                        x1   = xcosf(phis);
                        invc2= _mm512_rcp14_ps(_mm512_sub_ps(cpin,x1));
                        ear  = _mm512_mul_ps(cer,spin);
                        phis = _mm512_sub_ps(invc1,invc2);
                        eai  = _mm512_mul_ps(cei,spin);
                        phid = _mm512_add_ps(invc1,invc2);
                        _mm512_storeu_ps(&D1r[0] ,_mm512_mul_ps(ear,phis));
                        _mm512_storeu_ps(&D1i[0] ,_mm512_mul_ps(eai,phis));
                        _mm512_storeu_ps(&D2r[0] ,_mm512_mul_ps(ear,phid));
                        _mm512_storeu_ps(&D2i[0] ,_mm512_mul_ps(eai,phid));
                }
                
                
                /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Formula 8.1-62
                */
#include <cstdint>                
#include <complex>               
#include "GMS_cspint_quad.hpp"
#include "GMS_avint_quad.hpp"

                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_u(const float * __restrict pdAdl,
                                             const float *  __restrict pdl,
                                             const float   k0,
                                             const float   l) {
                          
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16] = {};
                         
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         register __m512 dAdl = _mm512_loadu_ps(&pdAdl[0]);
                         register __m512 dl   = _mm512_loadu_ps(&pdl[0]);                                 
                         register __m512 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                         register float rcs,k02,frac,sumr,sumi;
                         vk0  = _mm512_set1_ps(k0);
                         k0l  = _mm512_mul_ps(vk0,dl);
                         ear  = _mm512_setzero_ps();
                         eai  = _mm512_add_ps(k0l,k0l);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         _mm512_store_ps(&intr[0], _mm512_mul_ps(cer,dAdl);
                         _mm512_store_ps(&inti[0], _mm512_mul_ps(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         cspint(NTAB,pdl,&intr[0],0.0f,l,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,pdl,&inti[0],0.0f,l,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         c = {sumr,sumi};
                         k02   = k0*k0;   
                         frac  = k02/C314159265358979323846264338328;
                         rcs   = frac*std::abs(c);
                         return (rcs);                         
                  }
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pdAdl,
                                             const float * __restrict __ATTR_ALIGN__(64) pdl,
                                             const float   k0,
                                             const float   l) {
                          
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                         __ATTR_ALIGN__(64) float Y1[16]   = {};
                         __ATTR_ALIGN__(64) float Y2[16]   = {};
                         __ATTR_ALIGN__(64) float Y3[16]   = {}; 
                         __ATTR_ALIGN__(64) float E[16]    = {};
                         __ATTR_ALIGN__(64) float WRK[16] = {};
                         
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         register __m512 dAdl = _mm512_load_ps(&pdAdl[0]);
                         register __m512 dl   = _mm512_load_ps(&pdl[0]);                                 
                         register __m512 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                         register float rcs,k02,frac,sumr,sumi;
                         vk0  = _mm512_set1_ps(k0);
                         k0l  = _mm512_mul_ps(vk0,dl);
                         ear  = _mm512_setzero_ps();
                         eai  = _mm512_add_ps(k0l,k0l);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         _mm512_store_ps(&intr[0], _mm512_mul_ps(cer,dAdl);
                         _mm512_store_ps(&inti[0], _mm512_mul_ps(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         cspint(NTAB,pdl,&intr[0],0.0f,l,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumr);
                         cspint(NTAB,pdl,&inti[0],0.0f,l,&Y1[0],&Y2[0],&Y3[0],&E[0],&WRK[0],sumi);
                         c = {sumr,sumi};
                         k02   = k0*k0;   
                         frac  = k02/C314159265358979323846264338328;
                         rcs   = frac*std::abs(c);
                         return (rcs);                         
                  }
                  
                  
                   /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of small integrand. 
                       Integrator 'avint' i.e. irregular abscissas
                       Formula 8.1-62
                */
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_avint_u(const float * __restrict pdAdl,
                                                   const float *  __restrict pdl,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri) {
                          
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                                         
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         register __m512 dAdl = _mm512_loadu_ps(&pdAdl[0]);
                         register __m512 dl   = _mm512_loadu_ps(&pdl[0]);                                 
                         register __m512 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                         register float rcs,k02,frac,sumr,sumi;
                         int32_t err,eri;
                         vk0  = _mm512_set1_ps(k0);
                         k0l  = _mm512_mul_ps(vk0,dl);
                         ear  = _mm512_setzero_ps();
                         eai  = _mm512_add_ps(k0l,k0l);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         _mm512_store_ps(&intr[0], _mm512_mul_ps(cer,dAdl);
                         _mm512_store_ps(&inti[0], _mm512_mul_ps(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         sumr = avint(pdl,&intr[0],0.0f,l,err);
                         sumi = avint(pdl,&inti[0],0.0f,l,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr == 3 || ieri == 3) {
                            std::numeric_limits<float>::quiet_NaN();
                         }
                         c = {sumr,sumi};
                         k02   = k0*k0;   
                         frac  = k02/C314159265358979323846264338328;
                         rcs   = frac*std::abs(c);
                         return (rcs);                         
                  }
                  
                  
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_avint_a(const float * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                   const float * __restrict __ATTR_ALIGN__(64) pdl,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri) {
                          
                         __ATTR_ALIGN__(64) float intr[16] = {};
                         __ATTR_ALIGN__(64) float inti[16] = {}; 
                                         
                         constexpr int32_t NTAB = 16;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         register __m512 dAdl = _mm512_load_ps(&pdAdl[0]);
                         register __m512 dl   = _mm512_load_ps(&pdl[0]);                                 
                         register __m512 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                         register float rcs,k02,frac,sumr,sumi;
                         int32_t err,eri;
                         vk0  = _mm512_set1_ps(k0);
                         k0l  = _mm512_mul_ps(vk0,dl);
                         ear  = _mm512_setzero_ps();
                         eai  = _mm512_add_ps(k0l,k0l);
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         _mm512_stor_ps(&intr[0], _mm512_mul_ps(cer,dAdl);
                         _mm512_stor_ps(&inti[0], _mm512_mul_ps(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         sumr = avint(pdl,&intr[0],0.0f,l,err);
                         sumi = avint(pdl,&inti[0],0.0f,l,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr == 3 || ieri == 3) {
                            std::numeric_limits<float>::quiet_NaN();
                         }
                         c = {sumr,sumi};
                         k02   = k0*k0;   
                         frac  = k02/C314159265358979323846264338328;
                         rcs   = frac*std::abs(c);
                         return (rcs);                         
                  }
                  
                  
                  
                  
                  
                   /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of small integrand -- single-threaded execution.
                       Formula 8.1-62
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_u(const float * __restrict  pdAdl,
                                             const float * __restrict  pdl,
                                             float * __restrict  intr,
                                             float * __restrict  inti,
                                             float * __restrict  Y1,
                                             float * __restrict  Y2,
                                             float * __restrict  Y3,
                                             float * __restrict  E,
                                             float * __restrict  WRK
                                             const float   k0,
                                             const float   l,
                                             const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==16,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_zmm16r4_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_loadu_ps(&pdAdl[i]);
                             register __m512 y = _mm512_loadu_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_storeu_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_storeu_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
                      cspint(NTAB,pdl,intr,0.0f,l,Y1,Y2,Y3,E,WRK,sumr); 
                      cspint(NTAB,pdl,inti,0.0f,l,Y1,Y2,Y3,E,WRK,sumi);  
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pdAdl,
                                             const float * __restrict __ATTR_ALIGN__(64) pdl,
                                             float * __restrict __ATTR_ALIGN__(64) intr,
                                             float * __restrict __ATTR_ALIGN__(64) inti,
                                             float * __restrict __ATTR_ALIGN__(64) Y1,
                                             float * __restrict __ATTR_ALIGN__(64) Y2,
                                             float * __restrict __ATTR_ALIGN__(64) Y3,
                                             float * __restrict __ATTR_ALIGN__(64) E,
                                             float * __restrict __ATTR_ALIGN__(64) WRK
                                             const float   k0,
                                             const float   l,
                                             const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==16,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_zmm16r4_a(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_load_ps(&pdAdl[i]);
                             register __m512 y = _mm512_load_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_store_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_store_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
                      cspint(NTAB,pdl,intr,0.0f,l,Y1,Y2,Y3,E,WRK,sumr); 
                      cspint(NTAB,pdl,inti,0.0f,l,Y1,Y2,Y3,E,WRK,sumi);  
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                   /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of large integrand - single thread.
                       Integrator 'avint' i.e. irregular abscissas
                       Formula 8.1-62
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_avint_u(const float * __restrict  pdAdl,
                                                   const float * __restrict  pdl,
                                                   float * __restrict  intr,
                                                   float * __restrict  inti,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB) {
                                             
                                                
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_loadu_ps(&pdAdl[i]);
                             register __m512 y = _mm512_loadu_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_storeu_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_storeu_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
                      sumr = avint(pdl,intr,NTAB,0.0f,l,err); 
                      sumi = avint(pdl,inti,NTAB,0.0f,l,eri); 
                      ierr = err;
                      ieri = eri;
                      if(ierr == 3 || ieri == 3) {
                         return std::numerical_limits<float>::quiet_NaN();
                      } 
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_avint_a(const float * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                   const float * __restrict __ATTR_ALIGN__(64) pdl,
                                                   float * __restrict __ATTR_ALIGN__(64) intr,
                                                   float * __restrict __ATTR_ALIGN__(64) inti,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB) {
                                             
                                                
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_load_ps(&pdAdl[i]);
                             register __m512 y = _mm512_load_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_store_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_store_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
                      sumr = avint(pdl,intr,NTAB,0.0f,l,err); 
                      sumi = avint(pdl,inti,NTAB,0.0f,l,eri); 
                      ierr = err;
                      ieri = eri;
                      if(ierr == 3 || ieri == 3) {
                         return std::numerical_limits<float>::quiet_NaN();
                      } 
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
               
               
                 /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of large integrand -- two-threaded execution of integrator.
                       Integrator 'cspint'
                       Formula 8.1-62
                */
                
#include <omp.h>
    
                
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_cspint2t_u(const float * __restrict  pdAdl,
                                                     const float * __restrict  pdl,
                                                     float * __restrict  intr,
                                                     float * __restrict  inti,
                                                     struct RCS_F8162_DATA w,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==16,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_zmm16r4_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        float * __restrict px1 = w.Ya1;
                        float * __restrict py1 = w.Yb1;
                        float * __restrict px2 = w.Ya2;
                        float * __restrict py2 = w.Yb2;
                        float * __restrict px3 = w.Ya3;
                        float * __restrict py3 = w.Yb3;
                        float * __restrict px4 = w.Ea;
                        float * __restrict py4 = w.Eb;
                        float * __restrict px5 = w.WRKa;
                        float * __restrict py5 = w.WRKb;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_loadu_ps(&pdAdl[i]);
                             register __m512 y = _mm512_loadu_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_storeu_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_storeu_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
   #pragma omp parallel sctions
                {
                     #pragma omp section
                       {
                             cspint(NTAB,pdl,intr,0.0f,l,px1,px2,px3,px4,px5,sumr); 
                       }
                     #pragma omp section
                       {
                             cspint(NTAB,pdl,inti,0.0f,l,py1,py2,py3,py4,py5,sumi); 
                       }
                } 
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_cspint2t_a(const float * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                     const float * __restrict __ATTR_ALIGN__(64) pdl,
                                                     float * __restrict __ATTR_ALIGN__(64) intr,
                                                     float * __restrict __ATTR_ALIGN__(64) inti,
                                                     struct RCS_F8162_DATA w,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==16,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_zmm16r4_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        float * __restrict __ATTR_ALIGN__(64) px1 = w.Ya1;
                        float * __restrict __ATTR_ALIGN__(64) py1 = w.Yb1;
                        float * __restrict __ATTR_ALIGN__(64) px2 = w.Ya2;
                        float * __restrict __ATTR_ALIGN__(64) py2 = w.Yb2;
                        float * __restrict __ATTR_ALIGN__(64) px3 = w.Ya3;
                        float * __restrict __ATTR_ALIGN__(64) py3 = w.Yb3;
                        float * __restrict __ATTR_ALIGN__(64) px4 = w.Ea;
                        float * __restrict __ATTR_ALIGN__(64) py4 = w.Eb;
                        float * __restrict __ATTR_ALIGN__(64) px5 = w.WRKa;
                        float * __restrict __ATTR_ALIGN__(64) py5 = w.WRKb;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_load_ps(&pdAdl[i]);
                             register __m512 y = _mm512_load_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_store_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_store_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
   #pragma omp parallel sctions
                {
                     #pragma omp section
                       {
                             cspint(NTAB,pdl,intr,0.0f,l,px1,px2,px3,px4,px5,sumr); 
                       }
                     #pragma omp section
                       {
                             cspint(NTAB,pdl,inti,0.0f,l,py1,py2,py3,py4,py5,sumi); 
                       }
                } 
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of large integrand -- two-threaded execution of integrator.
                       Integrator 'avint' (irregular abscissas).
                       Formula 8.1-62
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_avint2t_u(const float * __restrict  pdAdl,
                                                     const float * __restrict  pdl,
                                                     float * __restrict  intr,
                                                     float * __restrict  inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                                              
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_loadu_ps(&pdAdl[i]);
                             register __m512 y = _mm512_loadu_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_storeu_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_storeu_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
   #pragma omp parallel sctions
                {
                     #pragma omp section
                       {
                            sumr = avint(pdl,intr,NTAB,0.0f,l,err); 
                       }
                     #pragma omp section
                       {
                            sumi = avint(pdl,inti,NTAB,0.0f,l,eri); 
                       }
                } 
                      ierr = err;
                      ieri = eri
                      if(ierr == 3 || ieri == 3) {
                          return std::numeric_limits<float>::quiet_NaN();
                      }
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   float rcs_f8162_zmm16r4_avint2t_a(const float * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                     const float * __restrict __ATTR_ALIGN__(64) pdl,
                                                     float * __restrict __ATTR_ALIGN__(64) intr,
                                                     float * __restrict __ATTR_ALIGN__(64) inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                                              
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                        register __m512 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                        register float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm512_set1_ps(k0);
                        ear  = _mm512_setzero_ps();
                        for(i = 0; i != ROUND_TO_SIXTEEN(NTAB,15); i += 16) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m512 x = _mm512_load_ps(&pdAdl[i]);
                             register __m512 y = _mm512_load_ps(&pdl[i]);
                             k0l               = _mm512_mul_ps(vk0,y);
                             eai               = _mm512_add_ps(k0l,k0l);
                             cexp_zmm16r4(ear,eai,&cer,&cei);
                             register __m512 t0 = cer;
                             register __m512 t1 = cei;
                             _mm512_store_ps(&intr[i], _mm512_mul_ps(t0,x));
                             _mm512_store_ps(&inti[i], _mm512_mul_ps(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const float x  = pdAdl[i];
                           const float y  = pdl[i];
                           const float k0l= k0*y;
                           const float eai= k0l+k0l;
                           const std::complex<float> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
   #pragma omp parallel sctions
                {
                     #pragma omp section
                       {
                            sumr = avint(pdl,intr,NTAB,0.0f,l,err); 
                       }
                     #pragma omp section
                       {
                            sumi = avint(pdl,inti,NTAB,0.0f,l,eri); 
                       }
                } 
                      ierr = err;
                      ieri = eri
                      if(ierr == 3 || ieri == 3) {
                          return std::numeric_limits<float>::quiet_NaN();
                      }
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
               
               
               
               
               /*
                     High frequency approximations.
                     Rounded-tip cone total nose-on
                     backscatter RCS.
                     Formula 8.1-93
               */
               
               
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f8193_zmm16r4(const __m512 b,
                                            const __m512 a,
                                            const __m512 k0,
                                            const __m512 alp,
                                            const __m512 l) {
                                            
                         const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                         const __m512 C1772453850905516027298167483341 = 
                                                     _mm512_set1_ps(1.772453850905516027298167483341f);
                         const __m512 C10                              = 
                                                     _mm512_set1_ps(1.0f);
                         const __m512 C15                              =
                                                     _mm512_set1_ps(1.5f);
                         register __m512 sina,pin,n,invn,x0,x1,k0l,sqr;
                         register __m512 ear,eai,cer,cei,cosa;
                         register __m512 rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm512_mul_ps(k0,b);
                         ear  = _mm512_setzero_ps();
                         n    = _mm512_add_ps(C15,_mm512_div_ps(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm512_mul_ps(k0,l);
                         eai  = _mm512_add_ps(k0l,k0l);
                         invn = _mm512_rcp14_ps(n);
                         sina = _mm512_sub_ps(C10,xsinf(alp));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosa = xcosf(alp);
                         x0   = xsinf(_mm512_mul_ps(
                                              _mm512_add_ps(k0b,k0b),sina));
                         x1   = _mm512_mul_ps(k0b,_mm512_mul_ps(cosa,cosa));
                         sqr  = _mm512_sub_ps(C10,_mm512_div_ps(x0,x1));
                         pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                         x0   = _mm512_sqrt_ps(sqr);
                         spin = _mm512_mul_ps(xsinf(pin),invn);
                         cpin = xcosf(pin);
                         x1   = _mm512_mul_ps(_mm512_add_ps(alp,alp),invn);
                         cos2a= xcosf(x1);
                         t0   = _mm512_mul_ps(C1772453850905516027298167483341,
                                                              _mm512_mul_ps(b,x0)); // keep
                         x1   = _mm512_sub_ps(cpin,cos2a);
                         t1   = _mm512_mul_ps(_mm512_add_ps(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm512_mul_ps(a,spin)); // keep
                         x0   = _mm512_rcp14_ps(x1); // keep
                         ear  = _mm512_mul_ps(_mm512_add_ps(t0,t1),x0);
                         t0   = _mm512_mul_ps(ear,cer);
                         t1   = _mm512_mul_ps(ear,cei);
                         rcs  = cabs_zmm16r4(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f8193_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pb,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) palp,
                                              const float * __restrict __ATTR_ALIGN__(64) pl) {
                                     
                         register __m512 b  = _mm512_load_ps(&pb[0]);
                         register __m512 a  = _mm512_load_ps(&pa[0]);  
                         register __m512 k0 = _mm512_load_ps(&pk0[0]);
                         register __m512 alp= _mm512_load_ps(&palp[0]);  
                         register __m512 l  = _mm512_load_ps(&pl[0]);   
                         const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                         const __m512 C1772453850905516027298167483341 = 
                                                     _mm512_set1_ps(1.772453850905516027298167483341f);
                         const __m512 C10                              = 
                                                     _mm512_set1_ps(1.0f);
                         const __m512 C15                              =
                                                     _mm512_set1_ps(1.5f);
                         register __m512 sina,pin,n,invn,x0,x1,k0l,sqr;
                         register __m512 ear,eai,cer,cei,cosa;
                         register __m512 rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm512_mul_ps(k0,b);
                         ear  = _mm512_setzero_ps();
                         n    = _mm512_add_ps(C15,_mm512_div_ps(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm512_mul_ps(k0,l);
                         eai  = _mm512_add_ps(k0l,k0l);
                         invn = _mm512_rcp14_ps(n);
                         sina = _mm512_sub_ps(C10,xsinf(alp));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosa = xcosf(alp);
                         x0   = xsinf(_mm512_mul_ps(
                                              _mm512_add_ps(k0b,k0b),sina));
                         x1   = _mm512_mul_ps(k0b,_mm512_mul_ps(cosa,cosa));
                         sqr  = _mm512_sub_ps(C10,_mm512_div_ps(x0,x1));
                         pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                         x0   = _mm512_sqrt_ps(sqr);
                         spin = _mm512_mul_ps(xsinf(pin),invn);
                         cpin = xcosf(pin);
                         x1   = _mm512_mul_ps(_mm512_add_ps(alp,alp),invn);
                         cos2a= xcosf(x1);
                         t0   = _mm512_mul_ps(C1772453850905516027298167483341,
                                                              _mm512_mul_ps(b,x0)); // keep
                         x1   = _mm512_sub_ps(cpin,cos2a);
                         t1   = _mm512_mul_ps(_mm512_add_ps(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm512_mul_ps(a,spin)); // keep
                         x0   = _mm512_rcp14_ps(x1); // keep
                         ear  = _mm512_mul_ps(_mm512_add_ps(t0,t1),x0);
                         t0   = _mm512_mul_ps(ear,cer);
                         t1   = _mm512_mul_ps(ear,cei);
                         rcs  = cabs_zmm16r4(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f8193_zmm16r4_u(const float * __restrict  pb,
                                              const float * __restrict  pa,
                                              const float * __restrict  pk0,
                                              const float * __restrict  palp,
                                              const float * __restrict _pl) {
                                     
                         register __m512 b  = _mm512_loadu_ps(&pb[0]);
                         register __m512 a  = _mm512_loadu_ps(&pa[0]);  
                         register __m512 k0 = _mm512_loadu_ps(&pk0[0]);
                         register __m512 alp= _mm512_loadu_ps(&palp[0]);  
                         register __m512 l  = _mm512_loadu_ps(&pl[0]);   
                         const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                         const __m512 C1772453850905516027298167483341 = 
                                                     _mm512_set1_ps(1.772453850905516027298167483341f);
                         const __m512 C10                              = 
                                                     _mm512_set1_ps(1.0f);
                         const __m512 C15                              =
                                                     _mm512_set1_ps(1.5f);
                         register __m512 sina,pin,n,invn,x0,x1,k0l,sqr;
                         register __m512 ear,eai,cer,cei,cosa;
                         register __m512 rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm512_mul_ps(k0,b);
                         ear  = _mm512_setzero_ps();
                         n    = _mm512_add_ps(C15,_mm512_div_ps(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm512_mul_ps(k0,l);
                         eai  = _mm512_add_ps(k0l,k0l);
                         invn = _mm512_rcp14_ps(n);
                         sina = _mm512_sub_ps(C10,xsinf(alp));
                         cexp_zmm16r4(ear,eai,&cer,&cei);
                         cosa = xcosf(alp);
                         x0   = xsinf(_mm512_mul_ps(
                                              _mm512_add_ps(k0b,k0b),sina));
                         x1   = _mm512_mul_ps(k0b,_mm512_mul_ps(cosa,cosa));
                         sqr  = _mm512_sub_ps(C10,_mm512_div_ps(x0,x1));
                         pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                         x0   = _mm512_sqrt_ps(sqr);
                         spin = _mm512_mul_ps(xsinf(pin),invn);
                         cpin = xcosf(pin);
                         x1   = _mm512_mul_ps(_mm512_add_ps(alp,alp),invn);
                         cos2a= xcosf(x1);
                         t0   = _mm512_mul_ps(C1772453850905516027298167483341,
                                                              _mm512_mul_ps(b,x0)); // keep
                         x1   = _mm512_sub_ps(cpin,cos2a);
                         t1   = _mm512_mul_ps(_mm512_add_ps(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm512_mul_ps(a,spin)); // keep
                         x0   = _mm512_rcp14_ps(x1); // keep
                         ear  = _mm512_mul_ps(_mm512_add_ps(t0,t1),x0);
                         t0   = _mm512_mul_ps(ear,cer);
                         t1   = _mm512_mul_ps(ear,cei);
                         rcs  = cabs_zmm16r4(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for |theta| = PI/2-alpha
                     Formula 8.1-96
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f8196_zmm16r4(const __m512 k0,
                                            const __m512 alp,
                                            const __m512 a,
                                            const __m512 b) {
                                            
                         
                          const __m512 C0444444444444444444444444444444 = 
                                                        _mm512_set1_ps(0.444444444444444444444444444444f);
                          register __m512 rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm512_mul_ps(C0444444444444444444444444444444,k0);
                          cosa  = xcosf(alp);
                          a32   = _mm512_mul_ps(a,_mm512_sqrt_ps(a));
                          x0    = _mm512_mul_ps(b,cosa);
                          sina  = xsinf(alp);
                          bca32 = _mm512_mul_ps(x0,_mm512_sqrt_ps(x0));
                          x1    = _mm512_div_ps(cosa,sina);
                          cota  = _mm512_mul_ps(x1,x1);
                          t1    = _mm512_mul_ps(t0,_mm512_mul_ps(cosa,cota));
                          x0    = _mm512_sub_ps(a32,bca32);
                          rcs   = _mm512_mul_ps(t1,_mm512_mul_ps(x0,x0));
                          return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f8196_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) pk0,
                                              const float * __restrict __ATTR_ALIGN__(64) palp,
                                              const float * __restrict __ATTR_ALIGN__(64) pa,
                                              const float * __restrict __ATTR_ALIGN__(64) pb) {
                                            
                         
                          register __m512 k0  = _mm512_load_ps(&pk0[0]);
                          register __m512 alp = _mm512_load_ps(&palp[0]);
                          register __m512 a   = _mm512_load_ps(&pa[0]);
                          register __m512 b   = _mm512_load_ps(&pb[0]);
                          const __m512 C0444444444444444444444444444444 = 
                                                        _mm512_set1_ps(0.444444444444444444444444444444f);
                          register __m512 rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm512_mul_ps(C0444444444444444444444444444444,k0);
                          cosa  = xcosf(alp);
                          a32   = _mm512_mul_ps(a,_mm512_sqrt_ps(a));
                          x0    = _mm512_mul_ps(b,cosa);
                          sina  = xsinf(alp);
                          bca32 = _mm512_mul_ps(x0,_mm512_sqrt_ps(x0));
                          x1    = _mm512_div_ps(cosa,sina);
                          cota  = _mm512_mul_ps(x1,x1);
                          t1    = _mm512_mul_ps(t0,_mm512_mul_ps(cosa,cota));
                          x0    = _mm512_sub_ps(a32,bca32);
                          rcs   = _mm512_mul_ps(t1,_mm512_mul_ps(x0,x0));
                          return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
                   __m512 rcs_f8196_zmm16r4_u(const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pa,
                                              const float * __restrict pb) {
                                            
                         
                          register __m512 k0  = _mm512_loadu_ps(&pk0[0]);
                          register __m512 alp = _mm512_loadu_ps(&palp[0]);
                          register __m512 a   = _mm512_loadu_ps(&pa[0]);
                          register __m512 b   = _mm512_loadu_ps(&pb[0]);
                          const __m512 C0444444444444444444444444444444 = 
                                                        _mm512_set1_ps(0.444444444444444444444444444444f);
                          register __m512 rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm512_mul_ps(C0444444444444444444444444444444,k0);
                          cosa  = xcosf(alp);
                          a32   = _mm512_mul_ps(a,_mm512_sqrt_ps(a));
                          x0    = _mm512_mul_ps(b,cosa);
                          sina  = xsinf(alp);
                          bca32 = _mm512_mul_ps(x0,_mm512_sqrt_ps(x0));
                          x1    = _mm512_div_ps(cosa,sina);
                          cota  = _mm512_mul_ps(x1,x1);
                          t1    = _mm512_mul_ps(t0,_mm512_mul_ps(cosa,cota));
                          x0    = _mm512_sub_ps(a32,bca32);
                          rcs   = _mm512_mul_ps(t1,_mm512_mul_ps(x0,x0));
                          return (rcs);
                 }
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for 0<|theta|<alpha
                     Perpendicular RCS.
                     Formula 8.1-94
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 rcs_perpendicular_f8194_zmm16r4(const __m512 h,
	                                        const __m512 l,
	                                        const __m512 b,
	                                        const __m512 a,
	                                        const __m512 k0,
	                                        const __m512 tht,
	                                        const __m512 alp) {
	                                 
	                                  
	                 const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                         const __m512 C1772453850905516027298167483341 = 
                                                     _mm512_set1_ps(1.772453850905516027298167483341f);
                         const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C10                              = 
                                                     _mm512_set1_ps(1.0f);  
                         const __m512 C15                              = 
                                                     _mm512_set1_ps(1.5f); 
                         const __m512 C05                              = 
                                                     _mm512_set1_ps(0.5f);
                         const __m512 C20                              =
                                                     _mm512_set1_ps(2.0f);
                         register __m512 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m512 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m512 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m512 cpin1,cpin2,trm1,trm2,rcs;
                         __m512 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm512_sub_ps(h,_mm512_add_ps(l,b));
                         sint = xsinf(tht);
                         k02  = _mm512_add_ps(k0,k0);
                         n    = _mm512_mul_ps(C15,_mm512_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm512_rcp14_ps(sint);
                         a2   = _mm512_mul_ps(a,C05);
                         ear  = _mm512_setzero_ps();
                         sk02 = _mm512_sqrt_ps(_mm512_mul_ps(k0,C05));
                         x0   = _mm512_mul_ps(hlb,_mm512_sub_ps(cost,b));
                         invn = _mm512_rcp14_ps(n);
                         //x2   = _mm512_mul_ps(a,C05);
                         eai  = _mm512_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm512_sqrt_ps(_mm512_mul_ps(a2,csct));
                         atant= _mm512_mul_ps(a,tant);
                         cost = xcosf(tht);
                         x0   = _mm512_mul_ps(b,C1772453850905516027298167483341);
                         cexp_zmm16r4(ear,eai,&cer1,&cei1);
                         cer1 = _mm512_mul_ps(x0,cer1);
                         spin = xsinf(pin);
                         cei1 = _mm512_mul_ps(x0,cei1);
                         cpin = xcosf(pin);
                         x1   = _mm512_mul_ps(_mm512_sub_ps(h,atant),cost);
                         eai2 = _mm512_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_zmm16r4(ear,eai2,&cer2,&cei);
                         x0   = _mm512_div_ps(spin,_mm512_mul_ps(n,sk02));
                         x1   = _mm512_mul_ps(x0,sacs);
                         cer2 = _mm512_mul_ps(cer2,x1);
                         cei2 = _mm512_mul_ps(cei2,x1);
                         cpin1= _mm512_rcp14_ps(_mm512_sub_ps(cpin,C10));
                         x2   = _mm512_mul_ps(C20,_mm512_add_ps(alp,tht));
                         cpin2= xcosf(_mm512_mul_ps(x2,invn));
                         x3   = _mm512_rcp14_ps(_mm512_sub_ps(cpin,cpin2));
                         trm1 = _mm512_sub_ps(cpin1,x3);
                         cmul_zmm16r4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm512_mul_ps(t0r,trm1);
                         t0i  = _mm512_mul_ps(t0i,trm1);
                         x0   = _mm512_mul_ps(C20,_mm512_sub_ps(alp,tht));
                         cpin2= xcosf(_mm512_mul_ps(x0,invn));
                         x1   = _mm512_rcp14_ps(cpin2);
                         trm2 = _mm512_sub_ps(cpin1,x1);
                         x2   = _mm512_fmadd_ps(cost,_mm512_mul_ps(k02,
                                                               _mm512_add_ps(h,atant)));
                         eai3 = _mm512_add_ps(C078539816339744830961566084582,x2);
                         cexp_zmm16r4(ear,ea3,&cer3,&cei3);
                         x0   = _mm512_div_ps(spin,_mm512_mul_ps(n,sk02));
                         x1   = _mm512_sqrt_ps(_mm512_mul_ps(gms::math::
                                                                  negate_zmm16r4(a2),csct));
                         x2   = _mm512_mul_ps(x0,x1);
                         cer3 = _mm512_mul_ps(_mm512_mul_ps(cer3,x2),trm2);
                         cei3 = _mm512_mul_ps(_mm512_mul_ps(cei3,x2),trm2);
                         cmul_zmm16r4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_zmm16r4(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	                                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512 rcs_perpendicular_f8194_zmm16r4_a(const float * __restrict __ATTR_ALIGN__(64) ph,
	                                         const float * __restrict __ATTR_ALIGN__(64) pl,
	                                         const float * __restrict __ATTR_ALIGN__(64) pb,
	                                         const float * __restrict __ATTR_ALIGN__(64) pa,
	                                         const float * __restrict __ATTR_ALIGN__(64) pk0,
	                                         const float * __restrict __ATTR_ALIGN__(64) ptht,
	                                         const float * __restrict __ATTR_ALIGN__(64) palp) {
	                                 
	                  
	                 register __m512 h  = _mm512_load_ps(&ph[0]);
	                 register __m512 l  = _mm512_load_ps(&pl[0]); 
	                 register __m512 b  = _mm512_load_ps(&pb[0]);   
	                 register __m512 a  = _mm512_load_ps(&pa[0]);  
	                 register __m512 k0 = _mm512_load_ps(&pk0[0]);
	                 register __m512 tht= _mm512_load_ps(&ptht[0]); 
	                 register __m512 alp= _mm512_load_ps(&palp[0]);        
	                 const __m512 C314159265358979323846264338328  = 
                                                     _mm512_set1_ps(3.14159265358979323846264338328f); 
                         const __m512 C1772453850905516027298167483341 = 
                                                     _mm512_set1_ps(1.772453850905516027298167483341f);
                         const __m512 C078539816339744830961566084582  = 
                                                     _mm512_set1_ps(0.78539816339744830961566084582f);
                         const __m512 C10                              = 
                                                     _mm512_set1_ps(1.0f);  
                         const __m512 C15                              = 
                                                     _mm512_set1_ps(1.5f); 
                         const __m512 C05                              = 
                                                     _mm512_set1_ps(0.5f);
                         const __m512 C20                              =
                                                     _mm512_set1_ps(2.0f);
                         __m512 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m512 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m512 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m512 cpin1,cpin2,trm1,trm2,rcs;
                         __m512 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm512_sub_ps(h,_mm512_add_ps(l,b));
                         sint = xsinf(tht);
                         k02  = _mm512_add_ps(k0,k0);
                         n    = _mm512_mul_ps(C15,_mm512_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm512_rcp14_ps(sint);
                         a2   = _mm512_mul_ps(a,C05);
                         ear  = _mm512_setzero_ps();
                         sk02 = _mm512_sqrt_ps(_mm512_mul_ps(k0,C05));
                         x0   = _mm512_mul_ps(hlb,_mm512_sub_ps(cost,b));
                         invn = _mm512_rcp14_ps(n);
                         //x2   = _mm512_mul_ps(a,C05);
                         eai  = _mm512_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm512_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm512_sqrt_ps(_mm512_mul_ps(a2,csct));
                         atant= _mm512_mul_ps(a,tant);
                         cost = xcosf(tht);
                         x0   = _mm512_mul_ps(b,C1772453850905516027298167483341);
                         cexp_zmm16r4(ear,eai,&cer1,&cei1);
                         cer1 = _mm512_mul_ps(x0,cer1);
                         spin = xsinf(pin);
                         cei1 = _mm512_mul_ps(x0,cei1);
                         cpin = xcosf(pin);
                         x1   = _mm512_mul_ps(_mm512_sub_ps(h,atant),cost);
                         eai2 = _mm512_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_zmm16r4(ear,eai2,&cer2,&cei);
                         x0   = _mm512_div_ps(spin,_mm512_mul_ps(n,sk02));
                         x1   = _mm512_mul_ps(x0,sacs);
                         cer2 = _mm512_mul_ps(cer2,x1);
                         cei2 = _mm512_mul_ps(cei2,x1);
                         cpin1= _mm512_rcp14_ps(_mm512_sub_ps(cpin,C10));
                         x2   = _mm512_mul_ps(C20,_mm512_add_ps(alp,tht));
                         cpin2= xcosf(_mm512_mul_ps(x2,invn));
                         x3   = _mm512_rcp14_ps(_mm512_sub_ps(cpin,cpin2));
                         trm1 = _mm512_sub_ps(cpin1,x3);
                         cmul_zmm16r4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm512_mul_ps(t0r,trm1);
                         t0i  = _mm512_mul_ps(t0i,trm1);
                         x0   = _mm512_mul_ps(C20,_mm512_sub_ps(alp,tht));
                         cpin2= xcosf(_mm512_mul_ps(x0,invn));
                         x1   = _mm512_rcp14_ps(cpin2);
                         trm2 = _mm512_sub_ps(cpin1,x1);
                         x2   = _mm512_fmadd_ps(cost,_mm512_mul_ps(k02,
                                                               _mm512_add_ps(h,atant)));
                         eai3 = _mm512_add_ps(C078539816339744830961566084582,x2);
                         cexp_zmm16r4(ear,ea3,&cer3,&cei3);
                         x0   = _mm512_div_ps(spin,_mm512_mul_ps(n,sk02));
                         x1   = _mm512_sqrt_ps(_mm512_mul_ps(gms::math::
                                                                  negate_zmm16r4(a2),csct));
                         x2   = _mm512_mul_ps(x0,x1);
                         cer3 = _mm512_mul_ps(_mm512_mul_ps(cer3,x2),trm2);
                         cei3 = _mm512_mul_ps(_mm512_mul_ps(cei3,x2),trm2);
                         cmul_zmm16r4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_zmm16r4(t1r,t1i);
                         return (rcs);
	        }
	        
                 
                 
                 
               
               
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_ZMM16R4_HPP__*/
