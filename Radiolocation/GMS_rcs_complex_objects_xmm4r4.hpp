

#ifndef __GMS_RCS_COMPLEX_OBJECTS_XMM4R4_HPP__
#define __GMS_RCS_COMPLEX_OBJECTS_XMM4R4_HPP__ 250920240807

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

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_XMM4R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_XMM4R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_XMM4R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_XMM4R4_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_XMM4R4_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_XMM4R4_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_XMM4R4_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_XMM4R4_CREATION_DATE = "25-09-2024 08:07 PM +00200 (WED 25 SEP 2024 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_XMM4R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_XMM4R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_XMM4R4_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_xmm4r4.hpp"
#include "GMS_simd_utils.hpp"

#if !defined(__AVX512F__) || !defined(__AVX512VL__)
#error "Support of ISA AVX512F or AVX512VL is required!!"
#endif



namespace  gms {


         namespace radiolocation {
         
         
               /*
                   Work (input) arrays for kernel rcs_f8162_xmm4r4_2t_u and
                   rcs_f8162_xmm4r4_2t_a.
               */
               struct RCS_F8162_DATA {
               
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
	           static inline
                   void coef_D12_f8121_xmm4r4(const __m128 gam,
                                             const __m128 phi,
                                             const __m128 k0,
                                             __m128 * __restrict D1r,
                                             __m128 * __restrict D1i,
                                            __m128 * __restrict D2r,
                                            __m128 * __restrict D2i) {
                                            
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 C10                              = _mm_set1_ps(1.0f);
                         __m128 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm_add_ps(phi,phi);
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm_setzero_ps();
                        invn = _mm_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = _mm_sin_ps(x1);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        cos1 = _mm_cos_ps(x1);
                        x0   = _mm_mul_ps(invn,phi2);
                        cos2 = _mm_cos_ps(x0);
                        cer  = _mm_mul_ps(cer,sinp);
                        trm1 = _mm_rcp14_ps(_mm_sub_ps(cos1),C10);
                        cei  = _mm_mul_ps(cei,sinp);
                        trm2 = _mm_rcp14_ps(_mm_sub_ps(cos1,cos2));
                        sqr  = _mm_mul_ps(invn,sqr);
                        ear  = _mm_mul_ps(cer,sqr);
                        eai  = _mm_mul_ps(cei,sqr);
                        x0   = _mm_sub_ps(trm1,trm2);
                        *D1r = _mm_mul_ps(ear,x0);
                        *D1i = _mm_mul_ps(eai,x0);
                        x1   = _mm_add_ps(trm1,trm2);
                        *D2r = _mm_mul_ps(ear,x1);
                        *D2i = _mm_mul_ps(eai,x1);
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8121_xmm4r4_a(const float * __restrict pgam,
                                               const float * __restrict pphi,
                                               const float * __restrict pk0,
                                               float * __restrict D1r,
                                               float * __restrict D1i,
                                               float * __restrict D2r,
                                               float * __restrict D2i) {
                                   
                         __m128 gam = _mm_load_ps(&pgam[0]);
                         __m128 phi = _mm_load_ps(&pphi[0]);  
                         __m128 k0  = _mm_load_ps(&pk0[0]);       
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 C10                              = _mm_set1_ps(1.0f);
                         __m128 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm_add_ps(phi,phi);
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm_setzero_ps();
                        invn = _mm_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = _mm_sin_ps(x1);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        cos1 = _mm_cos_ps(x1);
                        x0   = _mm_mul_ps(invn,phi2);
                        cos2 = _mm_cos_ps(x0);
                        cer  = _mm_mul_ps(cer,sinp);
                        trm1 = _mm_rcp14_ps(_mm_sub_ps(cos1),C10);
                        cei  = _mm_mul_ps(cei,sinp);
                        trm2 = _mm_rcp14_ps(_mm_sub_ps(cos1,cos2));
                        sqr  = _mm_mul_ps(invn,sqr);
                        ear  = _mm_mul_ps(cer,sqr);
                        eai  = _mm_mul_ps(cei,sqr);
                        x0   = _mm_sub_ps(trm1,trm2);
                        _mm_store_ps(&D1r[0] ,_mm_mul_ps(ear,x0));
                        _mm_store_ps(&D1i[0] ,_mm_mul_ps(eai,x0));
                        x1   = _mm_add_ps(trm1,trm2);
                        _mm_store_ps(&D2r[0] ,_mm_mul_ps(ear,x1));
                        _mm_store_ps(&D2i[0] ,_mm_mul_ps(eai,x1));
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8121_xmm4r4_u(const float * __restrict  pgam,
                                               const float * __restrict  pphi,
                                               const float * __restrict  pk0,
                                               float * __restrict  D1r,
                                               float * __restrict  D1i,
                                               float * __restrict  D2r,
                                               float * __restrict  D2i) {
                                   
                         __m128 gam = _mm_loadu_ps(&pgam[0]);
                         __m128 phi = _mm_loadu_ps(&pphi[0]);  
                         __m128 k0  = _mm_loadu_ps(&pk0[0]);       
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 C10                              = _mm_set1_ps(1.0f);
                         __m128 invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm_add_ps(phi,phi);
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328); 
                        ear  = _mm_setzero_ps();
                        invn = _mm_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm_mul_ps(C314159265358979323846264338328,invn);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(k0,
                                                C6283185307179586476925286766559));
                        sinp = _mm_sin_ps(x1);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        cos1 = _mm_cos_ps(x1);
                        x0   = _mm_mul_ps(invn,phi2);
                        cos2 = _mm_cos_ps(x0);
                        cer  = _mm_mul_ps(cer,sinp);
                        trm1 = _mm_rcp14_ps(_mm_sub_ps(cos1),C10);
                        cei  = _mm_mul_ps(cei,sinp);
                        trm2 = _mm_rcp14_ps(_mm_sub_ps(cos1,cos2));
                        sqr  = _mm_mul_ps(invn,sqr);
                        ear  = _mm_mul_ps(cer,sqr);
                        eai  = _mm_mul_ps(cei,sqr);
                        x0   = _mm_sub_ps(trm1,trm2);
                        _mm_storeu_ps(&D1r[0] ,_mm_mul_ps(ear,x0));
                        _mm_storeu_ps(&D1i[0] ,_mm_mul_ps(eai,x0));
                        x1   = _mm_add_ps(trm1,trm2);
                        _mm_storeu_ps(&D2r[0] ,_mm_mul_ps(ear,x1));
                        _mm_storeu_ps(&D2i[0] ,_mm_mul_ps(eai,x1));
                }
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter singly diffracted far-zone fields (E,H).
                    Formula: 8.1-19, 8.1-20
                
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f811920_xmm4r4(    const __m128 betai,
                                                 const __m128 betas,
                                                 const __m128 gam,
                                                 const __m128 phi,
                                                 const __m128 k0,
                                                 const __m128 r,
                                                 const __m128 rho,
                                                 const __m128 psi,
                                                 __m128 * __restrict Esr,
                                                 __m128 * __restrict Esi,
                                                 __m128 * __restrict Hsr,
                                                 __m128 * __restrict Hsi) {
                                                 
                        __m128 ear,eai,cer,cei;
                        __m128 D1r,D1i,D2r,D2i,x0,x1;
                        __m128 rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm_fmadd_ps(k0,r,psi);
                       coef_D12_f8121_xmm4r4(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm_rcp14_ps(r);
                       ear   = _mm_setzero_ps();
                       cosbi = _mm_cos_ps(betai);
                       eai   = k0rp;
                       cosbs = _mm_cos_ps(betas);
                       cexp_xmm4c4(ear,eai,&cer,&cei);
                       cer   = _mm_mul_ps(cer,invr);    
                       rhos  = _mm_div_ps(rho,_mm_add_ps(cosbi,cosbs));
                       cei   = _mm_mul_ps(cei,invr);
                       sqrho = _mm_sqrt_ps(rhos);
                       x0    = _mm_mul_ps(sqrho,cer);
                       x1    = _mm_mul_ps(sqrho,cei);
                       *Esr  = _mm_mul_ps(D1r,x0);
                       *Hsr  = _mm_mul_ps(D2r,x0);
                       *Esi  = _mm_mul_ps(D1i,x1);
                       *Hsi  = _mm_mul_ps(D2i,x1);                               
            }
            
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f811920_xmm4r4_a(  const float * __restrict pbetai,
                                                 const float * __restrict pbetas,
                                                 const float * __restrict pgam,
                                                 const float * __restrict pphi,
                                                 const float * __restrict pk0,
                                                 const float * __restrict pr,
                                                 const float * __restrict prho,
                                                 const float * __restrict ppsi,
                                                 float * __restrict Esr,
                                                 float * __restrict Esi,
                                                 float * __restrict Hsr,
                                                 float * __restrict Hsi) {
                              
                        __m128 betai = _mm_load_ps(&pbetai[0]);
                        __m128 betas = _mm_load_ps(&pbetas[0]); 
                        __m128 gam   = _mm_load_ps(&pgam[0]);   
                        __m128 k0    = _mm_load_ps(&pk0[0]); 
                        __m128 r     = _mm_load_ps(&pr[0]);
                        __m128 rho   = _mm_load_ps(&prho[0]); 
                        __m128 psi   = _mm_load_ps(&ppsi[0]);             
                        __m128 ear,eai,cer,cei;
                        __m128 D1r,D1i,D2r,D2i,x0,x1;
                        __m128 rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm_fmadd_ps(k0,r,psi);
                       coef_D12_f8121_xmm4r4(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm_rcp14_ps(r);
                       ear   = _mm_setzero_ps();
                       cosbi = _mm_cos_ps(betai);
                       eai   = k0rp;
                       cosbs = _mm_cos_ps(betas);
                       cexp_xmm4c4(ear,eai,&cer,&cei);
                       cer   = _mm_mul_ps(cer,invr);    
                       rhos  = _mm_div_ps(rho,_mm_add_ps(cosbi,cosbs));
                       cei   = _mm_mul_ps(cei,invr);
                       sqrho = _mm_sqrt_ps(rhos);
                       x0    = _mm_mul_ps(sqrho,cer);
                       x1    = _mm_mul_ps(sqrho,cei);
                       _mm_store_ps(&Esr[0] ,_mm_mul_ps(D1r,x0));
                       _mm_store_ps(&Hsr[0] ,_mm_mul_ps(D2r,x0));
                       _mm_store_ps(&Esi[0] ,_mm_mul_ps(D1i,x1));
                       _mm_store_ps(&Hsi[0] ,_mm_mul_ps(D2i,x1));                               
            }
            
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f811920_xmm4r4_u(  const float * __restrict  pbetai,
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
                              
                        __m128 betai = _mm_loadu_ps(&pbetai[0]);
                        __m128 betas = _mm_loadu_ps(&pbetas[0]); 
                        __m128 gam   = _mm_loadu_ps(&pgam[0]);   
                        __m128 k0    = _mm_loadu_ps(&pk0[0]); 
                        __m128 r     = _mm_loadu_ps(&pr[0]);
                        __m128 rho   = _mm_loadu_ps(&prho[0]); 
                        __m128 psi   = _mm_loadu_ps(&ppsi[0]);             
                        __m128 ear,eai,cer,cei;
                        __m128 D1r,D1i,D2r,D2i,x0,x1;
                        __m128 rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm_fmadd_ps(k0,r,psi);
                       coef_D12_f8121_xmm4r4(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm_rcp14_ps(r);
                       ear   = _mm_setzero_ps();
                       cosbi = _mm_cos_ps(betai);
                       eai   = k0rp;
                       cosbs = _mm_cos_ps(betas);
                       cexp_xmm4c4(ear,eai,&cer,&cei);
                       cer   = _mm_mul_ps(cer,invr);    
                       rhos  = _mm_div_ps(rho,_mm_add_ps(cosbi,cosbs));
                       cei   = _mm_mul_ps(cei,invr);
                       sqrho = _mm_sqrt_ps(rhos);
                       x0    = _mm_mul_ps(sqrho,cer);
                       x1    = _mm_mul_ps(sqrho,cei);
                       _mm_storeu_ps(&Esr[0] ,_mm_mul_ps(D1r,x0));
                       _mm_storeu_ps(&Hsr[0] ,_mm_mul_ps(D2r,x0));
                       _mm_storeu_ps(&Esi[0] ,_mm_mul_ps(D1i,x1));
                       _mm_storeu_ps(&Hsi[0] ,_mm_mul_ps(D2i,x1));                               
            }
            
            
            /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Ray normal-incidence to one of edge faces.
                    Formula: 8.1-24
            */
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8124_xmm4r4(const __m128 k0,
                                               const __m128 gam,
                                               __m128 * __restrict D1r,
                                               __m128 * __restrict D1i,
                                               __m128 * __restrict D2r,
                                               __m128 * __restrict D2i) {
                                               
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 C10                              = _mm_set1_ps(1.0f);  
                         __m128 ear,eai,cer,cei,t0r,t0i;
                         __m128 x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm_setzero_ps();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328);  
                        x1   = _mm_mul_ps(k0,C6283185307179586476925286766559);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        invn = _mm_rcp14_ps(x0);
                        cer  = gms::math::negate_xmm4r4(cer);
                        invn2= _mm_add_ps(invn,invn);      
                        cei  = gms::math::negate_xmm4r4(cei);
                        arg1 = _mm_mul_ps(C314159265358979323846264338328,invn2)
                        sin1 = _mm_sin_ps(arg1);
                        arg2 = _mm_mul_ps(C314159265358979323846264338328,invn);
                        t0r  = _mm_div_ps(cer,x1);
                        cos1 = _mm_cos_ps(arg1);
                        x0   = _mm_div_ps(cos1,sin1);
                        t0i  = _mm_div_ps(cei,x1);   
                        sin2 = _mm_sin_ps(arg2);
                        cos2 = _mm_cos_ps(arg2);   
                        x1   = _mm_div_ps(cos2,sin2);
                        ear  = _mm_fmsub_ps(invn,x0,_mm_mul_ps(invn2,x1));
                        eai  = _mm_fmadd_ps(invn,x0,_mm_mul_ps(invn2,x1));
                        cer  = _mm_mul_ps(t0r,ear);
                        cei  = _mm_mul_ps(t0i,eai);
                        *D1r = cer;
                        *D2r = cer;
                        *D1i = cei;
                        *D2i = cei;       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8124_xmm4r4_a(const  float * __restrict pk0,
                                                 const  float * __restrict pgam,
                                                 float * __restrict D1r,
                                                 float * __restrict D1i,
                                                 float * __restrict D2r,
                                                 float * __restrict D2i) {
                                    
                         __m128 k0 = _mm_load_ps(&pk0[0]);
                         __m128 gam= _mm_load_ps(&pgam[0]);           
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 C10                              = _mm_set1_ps(1.0f);  
                         __m128 ear,eai,cer,cei,t0r,t0i;
                         __m128 x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm_setzero_ps();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328);  
                        x1   = _mm_mul_ps(k0,C6283185307179586476925286766559);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        invn = _mm_rcp14_ps(x0);
                        cer  = gms::math::negate_xmm4r4(cer);
                        invn2= _mm_add_ps(invn,invn);      
                        cei  = gms::math::negate_xmm4r4(cei);
                        arg1 = _mm_mul_ps(C314159265358979323846264338328,invn2)
                        sin1 = _mm_sin_ps(arg1);
                        arg2 = _mm_mul_ps(C314159265358979323846264338328,invn);
                        t0r  = _mm_div_ps(cer,x1);
                        cos1 = _mm_cos_ps(arg1);
                        x0   = _mm_div_ps(cos1,sin1);
                        t0i  = _mm_div_ps(cei,x1);   
                        sin2 = _mm_sin_ps(arg2);
                        cos2 = _mm_cos_ps(arg2);   
                        x1   = _mm_div_ps(cos2,sin2);
                        ear  = _mm_fmsub_ps(invn,x0,_mm_mul_ps(invn2,x1));
                        eai  = _mm_fmadd_ps(invn,x0,_mm_mul_ps(invn2,x1));
                        cer  = _mm_mul_ps(t0r,ear);
                        cei  = _mm_mul_ps(t0i,eai);
                        _mm_store_ps(&D1r[0] ,cer);
                        _mm_store_ps(&D2r[0] ,cer);
                        _mm_store_ps(&D1i[0] ,cei);
                        _mm_store_ps(&D2i[0] ,cei);       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8124_xmm4r4_u(const  float * __restrict  pk0,
                                                 const  float * __restrict pgam,
                                                 float * __restrict  D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict  D2r,
                                                 float * __restrict  D2i) {
                                    
                         __m128 k0 = _mm_loadu_ps(&pk0[0]);
                         __m128 gam= _mm_loadu_ps(&pgam[0]);           
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 C10                              = _mm_set1_ps(1.0f);  
                         __m128 ear,eai,cer,cei,t0r,t0i;
                         __m128 x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm_setzero_ps();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328);  
                        x1   = _mm_mul_ps(k0,C6283185307179586476925286766559);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        invn = _mm_rcp14_ps(x0);
                        cer  = gms::math::negate_xmm4r4(cer);
                        invn2= _mm_add_ps(invn,invn);      
                        cei  = gms::math::negate_xmm4r4(cei);
                        arg1 = _mm_mul_ps(C314159265358979323846264338328,invn2)
                        sin1 = _mm_sin_ps(arg1);
                        arg2 = _mm_mul_ps(C314159265358979323846264338328,invn);
                        t0r  = _mm_div_ps(cer,x1);
                        cos1 = _mm_cos_ps(arg1);
                        x0   = _mm_div_ps(cos1,sin1);
                        t0i  = _mm_div_ps(cei,x1);   
                        sin2 = _mm_sin_ps(arg2);
                        cos2 = _mm_cos_ps(arg2);   
                        x1   = _mm_div_ps(cos2,sin2);
                        ear  = _mm_fmsub_ps(invn,x0,_mm_mul_ps(invn2,x1));
                        eai  = _mm_fmadd_ps(invn,x0,_mm_mul_ps(invn2,x1));
                        cer  = _mm_mul_ps(t0r,ear);
                        cei  = _mm_mul_ps(t0i,eai);
                        _mm_storeu_ps(&D1r[0] ,cer);
                        _mm_storeu_ps(&D2r[0] ,cer);
                        _mm_storeu_ps(&D1i[0] ,cei);
                        _mm_storeu_ps(&D2i[0] ,cei);       
                }
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Backscatter direction axial caustic (for slightly diffracted rays).
                    Formula: 8.1-26
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_Ddiff_f8126_xmm4r4(const __m128 gam,
                                             const __m128 phi,
                                             const __m128 k0,
                                             __m128 * __restrict Dr,
                                             __m128 * __restrict Di) {
                                             
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);  
                        const __m128 C20                              = 
                                                     _mm_set1_ps(2.0f);
                         __m128 ear,eai,cer,cei;
                         __m128 invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm_add_ps(phi,phi);
                        x0   = _mm_mul_ps(C6283185307179586476925286766559,k0);
                        ear  = _mm_setzero_ps();
                        n    = _mm_div_ps(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm_rcp14_ps(n); 
                        sqr  = _mm_sqrt_ps(x0);
                        pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        spin = _mm_sin_ps(pin);
                        cer  = _mm_div_ps(cer,sqr);
                        cpin = _mm_cos_ps(pin);
                        cei  = _mm_div_ps(cei,sqr);
                        cphin= _mm_cos_ps(_mm_mul_ps(phi2,invn));
                        x0   = _mm_mul_ps(_mm_mul_ps(C20,invn),spin);
                        x1   = _mm_sub_ps(cpin,cphin);
                        n    = _mm_div_ps(x0,x1);
                        *Dr  = _mm_mul_ps(cer,n);
                        *Di  = _mm_mul_ps(cei,n);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_Ddiff_f8126_xmm4r4_a(const float * __restrict pgam,
                                               const float * __restrict pphi,
                                               const float * __restrict pk0,
                                               float * __restrict  Dr,
                                               float * __restrict  Di) {
                                    
                         __m128 gam  = _mm_load_ps(&pgam[0]);
                         __m128 phi  = _mm_load_ps(&pphi[0]);
                         __m128 k0   = _mm_load_ps(&pk0[0]);         
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);  
                        const __m128 C20                              = 
                                                     _mm_set1_ps(2.0f);
                         __m128 ear,eai,cer,cei;
                         __m128 invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm_add_ps(phi,phi);
                        x0   = _mm_mul_ps(C6283185307179586476925286766559,k0);
                        ear  = _mm_setzero_ps();
                        n    = _mm_div_ps(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm_rcp14_ps(n); 
                        sqr  = _mm_sqrt_ps(x0);
                        pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        spin = _mm_sin_ps(pin);
                        cer  = _mm_div_ps(cer,sqr);
                        cpin = _mm_cos_ps(pin);
                        cei  = _mm_div_ps(cei,sqr);
                        cphin= _mm_cos_ps(_mm_mul_ps(phi2,invn));
                        x0   = _mm_mul_ps(_mm_mul_ps(C20,invn),spin);
                        x1   = _mm_sub_ps(cpin,cphin);
                        n    = _mm_div_ps(x0,x1);
                        _mm_store_ps(&Dr[0] ,_mm_mul_ps(cer,n));
                        _mm_store_ps(&Di[0] ,_mm_mul_ps(cei,n));
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_Ddiff_f8126_xmm4r4_u(const float * __restrict  pgam,
                                                   const float * __restrict pphi,
                                                   const float * __restrict  pk0,
                                                   float * __restrict  Dr,
                                                   float * __restrict  Di) {
                                    
                         __m128 gam  = _mm_loadu_ps(&pgam[0]);
                         __m128 phi  = _mm_loadu_ps(&pphi[0]);
                         __m128 k0   = _mm_loadu_ps(&pk0[0]);         
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f);  
                        const __m128 C20                              = 
                                                     _mm_set1_ps(2.0f);
                         __m128 ear,eai,cer,cei;
                         __m128 invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm_add_ps(phi,phi);
                        x0   = _mm_mul_ps(C6283185307179586476925286766559,k0);
                        ear  = _mm_setzero_ps();
                        n    = _mm_div_ps(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm_rcp14_ps(n); 
                        sqr  = _mm_sqrt_ps(x0);
                        pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        spin = _mm_sin_ps(pin);
                        cer  = _mm_div_ps(cer,sqr);
                        cpin = _mm_cos_ps(pin);
                        cei  = _mm_div_ps(cei,sqr);
                        cphin= _mm_cos_ps(_mm_mul_ps(phi2,invn));
                        x0   = _mm_mul_ps(_mm_mul_ps(C20,invn),spin);
                        x1   = _mm_sub_ps(cpin,cphin);
                        n    = _mm_div_ps(x0,x1);
                        _mm_storeu_ps(&Dr[0] ,_mm_mul_ps(cer,n));
                        _mm_storeu_ps(&Di[0] ,_mm_mul_ps(cei,n));
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
	           
	          
                   
	           static inline
                   void EsHs_f8125_xmm4r4(const __m128 a,
                                           const __m128 k0,
                                           const __m128 r,
                                           const __m128 gam,
                                           const __m128 phi,
                                           const __m128 psi,
                                           __m128 * __restrict Esr,
                                           __m128 * __restrict Esi,
                                           __m128 * __restrict Hsr,
                                           __m128 * __restrict Hsi) {
                                           
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f); 
                        const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         __m128 Dr,Di,x0,x1;
                         __m128 ear,eai,cer,cei;
                         __m128 k0r,invr,t0r,t0i;
                         __m128 sqr;
                        k0r  = _mm_mul_ps(k0,r);
                        ear  = _mm_setzero_ps();
                        eai  = _mm_add_ps(k0r,_mm_sub_ps(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_xmm4r4(gam,phi,k0,Dr,Di); 
                        invr = _mm_rcp_ps(r);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(C6283185307179586476925286766559,k0));
                        t0r  = _mm_mul_ps(Dr,C05);
                        cer  = _mm_mul_ps(cer,invr);
                        sqr  = _mm_mul_ps(a,sqr);
                        t0i  = _mm_mul_ps(Di,C05);  
                        cei  = _mm_mul_ps(cei,invr);
                        x0   = _mm_mul_ps(t0r,_mm_mul_ps(sqr,cer));
                        x1   = _mm_mul_ps(t0i,_mm_mul_ps(sqr,cei));
                        *Esr = gms::math::negate_xmm4r4(x0);
                        *Esi = gms::math::negate_xmm4r4(x1); 
                        *Hsr = x0;
                        *Hsi = x1;             
                      
                }
                
                
                  __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f8125_xmm4r4_a(const float * __restrict pa,
                                             const float * __restrict pk0,
                                             const float * __restrict pr,
                                             const float * __restrict pgam,
                                             const float * __restrict pphi,
                                             const float * __restrict ppsi,
                                             float * __restrict Esr,
                                             float * __restrict Esi,
                                             float * __restrict Hsr,
                                             float * __restrict Hsi) {
                            
                         __m128 a   = _mm_load_ps(&pa[0]);
                         __m128 k0  = _mm_load_ps(&pk0[0]);   
                         __m128 r   = _mm_load_ps(&pr[0]);   
                         __m128 gam = _mm_load_ps(&pgam[0]);   
                         __m128 phi = _mm_load_ps(&pphi[0]); 
                         __m128 psi = _mm_load_ps(&ppsi[0]);     
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f); 
                        const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         __m128 Dr,Di,x0,x1;
                         __m128 ear,eai,cer,cei;
                         __m128 k0r,invr,t0r,t0i;
                         __m128 sqr;
                        k0r  = _mm_mul_ps(k0,r);
                        ear  = _mm_setzero_ps();
                        eai  = _mm_add_ps(k0r,_mm_sub_ps(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_xmm4r4(gam,phi,k0,Dr,Di); 
                        invr = _mm_rcp_ps(r);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(C6283185307179586476925286766559,k0));
                        t0r  = _mm_mul_ps(Dr,C05);
                        cer  = _mm_mul_ps(cer,invr);
                        sqr  = _mm_mul_ps(a,sqr);
                        t0i  = _mm_mul_ps(Di,C05);  
                        cei  = _mm_mul_ps(cei,invr);
                        x0   = _mm_mul_ps(t0r,_mm_mul_ps(sqr,cer));
                        x1   = _mm_mul_ps(t0i,_mm_mul_ps(sqr,cei));
                        _mm_store_ps(&Esr[0] ,gms::math::negate_xmm4r4(x0));
                        _mm_store_ps(&Esi[0] ,gms::math::negate_xmm4r4(x1)); 
                        _mm_store_ps(&Hsr[0] ,x0);
                        _mm_store_ps(&Hsi[0] ,x1);             
                      
                }
                
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f8125_xmm4r4_u(const float * __restrict  pa,
                                             const float * __restrict  pk0,
                                             const float * __restrict  pr,
                                             const float * __restrict pgam,
                                             const float * __restrict  pphi,
                                             const float * __restrict  ppsi,
                                             float * __restrict  Esr,
                                             float * __restrict  Esi,
                                             float * __restrict  Hsr,
                                             float * __restrict  Hsi) {
                            
                         __m128 a   = _mm_loadu_ps(&pa[0]);
                         __m128 k0  = _mm_loadu_ps(&pk0[0]);   
                         __m128 r   = _mm_loadu_ps(&pr[0]);   
                         __m128 gam = _mm_loadu_ps(&pgam[0]);   
                         __m128 phi = _mm_loadu_ps(&pphi[0]); 
                         __m128 psi = _mm_loadu_ps(&ppsi[0]);     
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f); 
                        const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         __m128 Dr,Di,x0,x1;
                         __m128 ear,eai,cer,cei;
                         __m128 k0r,invr,t0r,t0i;
                         __m128 sqr;
                        k0r  = _mm_mul_ps(k0,r);
                        ear  = _mm_setzero_ps();
                        eai  = _mm_add_ps(k0r,_mm_sub_ps(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_xmm4r4(gam,phi,k0,Dr,Di); 
                        invr = _mm_rcp_ps(r);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(C6283185307179586476925286766559,k0));
                        t0r  = _mm_mul_ps(Dr,C05);
                        cer  = _mm_mul_ps(cer,invr);
                        sqr  = _mm_mul_ps(a,sqr);
                        t0i  = _mm_mul_ps(Di,C05);  
                        cei  = _mm_mul_ps(cei,invr);
                        x0   = _mm_mul_ps(t0r,_mm_mul_ps(sqr,cer));
                        x1   = _mm_mul_ps(t0i,_mm_mul_ps(sqr,cei));
                        _mm_storeu_ps(&Esr[0] ,gms::math::negate_xmm4r4(x0));
                        _mm_storeu_ps(&Esi[0] ,gms::math::negate_xmm4r4(x1)); 
                        _mm_storeu_ps(&Hsr[0] ,x0);
                        _mm_storeu_ps(&Hsi[0] ,x1);             
                      
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
	           
	          
                   
	           static inline
                   void coef_D12_f8127_xmm4r4(const __m128 k0,
                                               const __m128 gam,
                                               const __m128 phi1,
                                               const __m128 phi2,
                                               __m128 * __restrict D1r,
                                               __m128 * __restrict D1i,
                                               __m128 * __restrict D2r,
                                               __m128 * __restrict D2i) {
                                               
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         __m128 invn,x0,x1,ear,eai,cer,cei;
                         __m128 sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328);
                        phid = _mm_mul_ps(_mm_sub_ps(phi1,phi2),invn);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm_setzero_ps();
                        phis = _mm_mul_ps(_mm_add_ps(phi1,phi2),invn);
                        invn = _mm_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        x0   = _mm_sin_ps(pin);
                        cer  = _mm_div_ps(cer,sqr);
                        spin = _mm_mul_ps(x0,invn);
                        cei  = _mm_div_ps(cei,sqr);
                        cpin = _mm_cos_ps(pin);
                        x0   = _mm_cos_ps(phid);
                        invc1= _mm_rcp14_ps(_mm_sub_ps(cpin,x0));
                        x1   = _mm_cos_ps(phis);
                        invc2= _mm_rcp14_ps(_mm_sub_ps(cpin,x1));
                        ear  = _mm_mul_ps(cer,spin);
                        phis = _mm_sub_ps(invc1,invc2);
                        eai  = _mm_mul_ps(cei,spin);
                        phid = _mm_add_ps(invc1,invc2);
                        *D1r = _mm_mul_ps(ear,phis);
                        *D1i = _mm_mul_ps(eai,phis);
                        *D2r = _mm_mul_ps(ear,phid);
                        *D2i = _mm_mul_ps(eai,phid);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8127_xmm4r4_a(const float * __restrict pk0,
                                                 const float * __restrict pgam,
                                                 const float * __restrict pphi1,
                                                 const float * __restrict pphi2,
                                                 float * __restrict  D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict  D2r,
                                                 float * __restrict  D2i) {
                              
                         __m128 k0  = _mm_load_ps(&pk0[0]);
                         __m128 gam = _mm_load_ps(&pgam[0]);
                         __m128 phi1= _mm_load_ps(&pphi1[0]);
                         __m128 phi2= _mm_load_ps(&pphi2[0]);                 
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         __m128 invn,x0,x1,ear,eai,cer,cei;
                         __m128 sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328);
                        phid = _mm_mul_ps(_mm_sub_ps(phi1,phi2),invn);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm_setzero_ps();
                        phis = _mm_mul_ps(_mm_add_ps(phi1,phi2),invn);
                        invn = _mm_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        x0   = _mm_sin_ps(pin);
                        cer  = _mm_div_ps(cer,sqr);
                        spin = _mm_mul_ps(x0,invn);
                        cei  = _mm_div_ps(cei,sqr);
                        cpin = _mm_cos_ps(pin);
                        x0   = _mm_cos_ps(phid);
                        invc1= _mm_rcp14_ps(_mm_sub_ps(cpin,x0));
                        x1   = _mm_cos_ps(phis);
                        invc2= _mm_rcp14_ps(_mm_sub_ps(cpin,x1));
                        ear  = _mm_mul_ps(cer,spin);
                        phis = _mm_sub_ps(invc1,invc2);
                        eai  = _mm_mul_ps(cei,spin);
                        phid = _mm_add_ps(invc1,invc2);
                        _mm_store_ps(&D1r[0] ,_mm_mul_ps(ear,phis));
                        _mm_store_ps(&D1i[0] ,_mm_mul_ps(eai,phis));
                        _mm_store_ps(&D2r[0] ,_mm_mul_ps(ear,phid));
                        _mm_store_ps(&D2i[0] ,_mm_mul_ps(eai,phid));
                }
                
                
                
                  __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8127_xmm4r4_u(const float * __restrict  pk0,
                                                 const float * __restrict  pgam,
                                                 const float * __restrict  pphi1,
                                                 const float * __restrict  pphi2,
                                                 float * __restrict   D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict   D2r,
                                                 float * __restrict  D2i) {
                              
                         __m128 k0  = _mm_loadu_ps(&pk0[0]);
                         __m128 gam = _mm_loadu_ps(&pgam[0]);
                         __m128 phi1= _mm_loadu_ps(&pphi1[0]);
                         __m128 phi2= _mm_loadu_ps(&pphi2[0]);                 
                        const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                        const __m128 C6283185307179586476925286766559 = 
                                                     _mm_set1_ps(6.283185307179586476925286766559f);
                        const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         __m128 invn,x0,x1,ear,eai,cer,cei;
                         __m128 sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm_div_ps(gam,C314159265358979323846264338328);
                        phid = _mm_mul_ps(_mm_sub_ps(phi1,phi2),invn);
                        sqr  = _mm_sqrt_ps(_mm_mul_ps(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm_setzero_ps();
                        phis = _mm_mul_ps(_mm_add_ps(phi1,phi2),invn);
                        invn = _mm_rcp14_ps(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                        cexp_xmm4c4(ear,eai,&cer,&cei);
                        x0   = _mm_sin_ps(pin);
                        cer  = _mm_div_ps(cer,sqr);
                        spin = _mm_mul_ps(x0,invn);
                        cei  = _mm_div_ps(cei,sqr);
                        cpin = _mm_cos_ps(pin);
                        x0   = _mm_cos_ps(phid);
                        invc1= _mm_rcp14_ps(_mm_sub_ps(cpin,x0));
                        x1   = _mm_cos_ps(phis);
                        invc2= _mm_rcp14_ps(_mm_sub_ps(cpin,x1));
                        ear  = _mm_mul_ps(cer,spin);
                        phis = _mm_sub_ps(invc1,invc2);
                        eai  = _mm_mul_ps(cei,spin);
                        phid = _mm_add_ps(invc1,invc2);
                        _mm_storeu_ps(&D1r[0] ,_mm_mul_ps(ear,phis));
                        _mm_storeu_ps(&D1i[0] ,_mm_mul_ps(eai,phis));
                        _mm_storeu_ps(&D2r[0] ,_mm_mul_ps(ear,phid));
                        _mm_storeu_ps(&D2i[0] ,_mm_mul_ps(eai,phid));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_u(const float * __restrict pdAdl,
                                             const float *  __restrict pdl,
                                             const float   k0,
                                             const float   l) {
                          
                         float intr[4] = {};
                         float inti[4] = {}; 
                         float Y1[4]   = {};
                         float Y2[4]   = {};
                         float Y3[4]   = {}; 
                         float E[4]    = {};
                         float WRK[4] = {};
                         
                         constexpr int32_t NTAB = 4;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m128 dAdl = _mm_loadu_ps(&pdAdl[0]);
                          __m128 dl   = _mm_loadu_ps(&pdl[0]);                                 
                          __m128 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                          float rcs,k02,frac,sumr,sumi;
                         vk0  = _mm_set1_ps(k0);
                         k0l  = _mm_mul_ps(vk0,dl);
                         ear  = _mm_setzero_ps();
                         eai  = _mm_add_ps(k0l,k0l);
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         _mm_store_ps(&intr[0], _mm_mul_ps(cer,dAdl);
                         _mm_store_ps(&inti[0], _mm_mul_ps(cei,dAdl);
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_a(const float * __restrict pdAdl,
                                             const float * __restrict pdl,
                                             const float   k0,
                                             const float   l) {
                          
                         float intr[4] = {};
                         float inti[4] = {}; 
                         float Y1[4]   = {};
                         float Y2[4]   = {};
                         float Y3[4]   = {}; 
                         float E[4]    = {};
                         float WRK[4] = {};
                         
                         constexpr int32_t NTAB = 4;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m128 dAdl = _mm_load_ps(&pdAdl[0]);
                          __m128 dl   = _mm_load_ps(&pdl[0]);                                 
                          __m128 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                          float rcs,k02,frac,sumr,sumi;
                         vk0  = _mm_set1_ps(k0);
                         k0l  = _mm_mul_ps(vk0,dl);
                         ear  = _mm_setzero_ps();
                         eai  = _mm_add_ps(k0l,k0l);
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         _mm_store_ps(&intr[0], _mm_mul_ps(cer,dAdl);
                         _mm_store_ps(&inti[0], _mm_mul_ps(cei,dAdl);
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_avint_u(const float * __restrict pdAdl,
                                                   const float *  __restrict pdl,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri) {
                          
                         float intr[4] = {};
                         float inti[4] = {}; 
                                         
                         constexpr int32_t NTAB = 4;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m128 dAdl = _mm_loadu_ps(&pdAdl[0]);
                          __m128 dl   = _mm_loadu_ps(&pdl[0]);                                 
                          __m128 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                          float rcs,k02,frac,sumr,sumi;
                         int32_t err,eri;
                         vk0  = _mm_set1_ps(k0);
                         k0l  = _mm_mul_ps(vk0,dl);
                         ear  = _mm_setzero_ps();
                         eai  = _mm_add_ps(k0l,k0l);
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         _mm_store_ps(&intr[0], _mm_mul_ps(cer,dAdl);
                         _mm_store_ps(&inti[0], _mm_mul_ps(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         sumr = avint(pdl,&intr[0],0.0f,l,err);
                         sumi = avint(pdl,&inti[0],0.0f,l,eri);
                         ierr = err;
                         ieri = eri;
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_avint_a(const float * __restrict pdAdl,
                                                   const float * __restrict pdl,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri) {
                          
                         float intr[4] = {};
                         float inti[4] = {}; 
                                         
                         constexpr int32_t NTAB = 4;               
                         constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                          __m128 dAdl = _mm_load_ps(&pdAdl[0]);
                          __m128 dl   = _mm_load_ps(&pdl[0]);                                 
                          __m128 vk0,k0l,ear,eai,cer,cei;
                         std::complex<float> c;
                          float rcs,k02,frac,sumr,sumi;
                         int32_t err,eri;
                         vk0  = _mm_set1_ps(k0);
                         k0l  = _mm_mul_ps(vk0,dl);
                         ear  = _mm_setzero_ps();
                         eai  = _mm_add_ps(k0l,k0l);
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         _mm_stor_ps(&intr[0], _mm_mul_ps(cer,dAdl);
                         _mm_stor_ps(&inti[0], _mm_mul_ps(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         sumr = avint(pdl,&intr[0],0.0f,l,err);
                         sumi = avint(pdl,&inti[0],0.0f,l,eri);
                         ierr = err;
                         ieri = eri;
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
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of small integrand -- single-threaded execution.
                       Formula 8.1-62
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_u(const float * __restrict  pdAdl,
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
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_xmm4r4_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m128 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_loadu_ps(&pdAdl[i]);
                              __m128 y = _mm_loadu_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_storeu_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_storeu_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_a(const float * __restrict pdAdl,
                                             const float * __restrict pdl,
                                             float * __restrict intr,
                                             float * __restrict inti,
                                             float * __restrict Y1,
                                             float * __restrict Y2,
                                             float * __restrict Y3,
                                             float * __restrict E,
                                             float * __restrict WRK
                                             const float   k0,
                                             const float   l,
                                             const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_xmm4r4_a(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m128 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_load_ps(&pdAdl[i]);
                              __m128 y = _mm_load_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_store_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_store_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_avint_u(const float * __restrict  pdAdl,
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
                         __m128 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3);i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_loadu_ps(&pdAdl[i]);
                              __m128 y = _mm_loadu_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_storeu_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_storeu_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_avint_a(const float * __restrict pdAdl,
                                                   const float * __restrict pdl,
                                                   float * __restrict intr,
                                                   float * __restrict inti,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB) {
                                             
                                                
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m128 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3);i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_load_ps(&pdAdl[i]);
                              __m128 y = _mm_load_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_store_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_store_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_cspint2t_u(const float * __restrict  pdAdl,
                                                     const float * __restrict  pdl,
                                                     float * __restrict  intr,
                                                     float * __restrict  inti,
                                                     struct RCS_F8162_DATA w,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_xmm4r4_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m128 vk0,k0l,ear,eai,cer,cei;
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
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3);i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_loadu_ps(&pdAdl[i]);
                              __m128 y = _mm_loadu_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_storeu_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_storeu_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_cspint2t_a(const float * __restrict pdAdl,
                                                     const float * __restrict pdl,
                                                     float * __restrict intr,
                                                     float * __restrict inti,
                                                     struct RCS_F8162_DATA w,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            float rcs = 0.0f;
                            rcs = rcs_f8162_xmm4r4_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m128 vk0,k0l,ear,eai,cer,cei;
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
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3);i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_load_ps(&pdAdl[i]);
                              __m128 y = _mm_load_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_store_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_store_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_avint2t_u(const float * __restrict  pdAdl,
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
                         __m128 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3);i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_loadu_ps(&pdAdl[i]);
                              __m128 y = _mm_loadu_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_storeu_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_storeu_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   float rcs_f8162_xmm4r4_avint2t_a(const float * __restrict pdAdl,
                                                     const float * __restrict pdl,
                                                     float * __restrict intr,
                                                     float * __restrict inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB) {
                                             
                                              
                        constexpr float C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328f;
                         __m128 vk0,k0l,ear,eai,cer,cei;
                        std::complex<float> c;
                         float rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm_set1_ps(k0);
                        ear  = _mm_setzero_ps();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3);i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                              __m128 x = _mm_load_ps(&pdAdl[i]);
                              __m128 y = _mm_load_ps(&pdl[i]);
                             k0l               = _mm_mul_ps(vk0,y);
                             eai               = _mm_add_ps(k0l,k0l);
                             cexp_xmm4c4(ear,eai,&cer,&cei);
                              __m128 t0 = cer;
                              __m128 t1 = cei;
                             _mm_store_ps(&intr[i], _mm_mul_ps(t0,x));
                             _mm_store_ps(&inti[i], _mm_mul_ps(t1,x));
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
	           
	          
                   
	           static inline
                   __m128 rcs_f8193_xmm4r4(const __m128 b,
                                            const __m128 a,
                                            const __m128 k0,
                                            const __m128 alp,
                                            const __m128 l) {
                                            
                         const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);
                         const __m128 C15                              =
                                                     _mm_set1_ps(1.5f);
                          __m128 sina,pin,n,invn,x0,x1,k0l,sqr;
                          __m128 ear,eai,cer,cei,cosa;
                          __m128 rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm_mul_ps(k0,b);
                         ear  = _mm_setzero_ps();
                         n    = _mm_add_ps(C15,_mm_div_ps(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm_mul_ps(k0,l);
                         eai  = _mm_add_ps(k0l,k0l);
                         invn = _mm_rcp14_ps(n);
                         sina = _mm_sub_ps(C10,_mm_sin_ps(alp));
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         cosa = _mm_cos_ps(alp);
                         x0   = _mm_sin_ps(_mm_mul_ps(
                                              _mm_add_ps(k0b,k0b),sina));
                         x1   = _mm_mul_ps(k0b,_mm_mul_ps(cosa,cosa));
                         sqr  = _mm_sub_ps(C10,_mm_div_ps(x0,x1));
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         x0   = _mm_sqrt_ps(sqr);
                         spin = _mm_mul_ps(_mm_sin_ps(pin),invn);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_add_ps(alp,alp),invn);
                         cos2a= _mm_cos_ps(x1);
                         t0   = _mm_mul_ps(C1772453850905516027298167483341,
                                                              _mm_mul_ps(b,x0)); // keep
                         x1   = _mm_sub_ps(cpin,cos2a);
                         t1   = _mm_mul_ps(_mm_add_ps(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm_mul_ps(a,spin)); // keep
                         x0   = _mm_rcp14_ps(x1); // keep
                         ear  = _mm_mul_ps(_mm_add_ps(t0,t1),x0);
                         t0   = _mm_mul_ps(ear,cer);
                         t1   = _mm_mul_ps(ear,cei);
                         rcs  = cabs_xmm4c4(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m128 rcs_f8193_xmm4r4_a(const float * __restrict pb,
                                              const float * __restrict pa,
                                              const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pl) {
                                     
                          __m128 b  = _mm_load_ps(&pb[0]);
                          __m128 a  = _mm_load_ps(&pa[0]);  
                          __m128 k0 = _mm_load_ps(&pk0[0]);
                          __m128 alp= _mm_load_ps(&palp[0]);  
                          __m128 l  = _mm_load_ps(&pl[0]);   
                         const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);
                         const __m128 C15                              =
                                                     _mm_set1_ps(1.5f);
                          __m128 sina,pin,n,invn,x0,x1,k0l,sqr;
                          __m128 ear,eai,cer,cei,cosa;
                          __m128 rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm_mul_ps(k0,b);
                         ear  = _mm_setzero_ps();
                         n    = _mm_add_ps(C15,_mm_div_ps(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm_mul_ps(k0,l);
                         eai  = _mm_add_ps(k0l,k0l);
                         invn = _mm_rcp14_ps(n);
                         sina = _mm_sub_ps(C10,_mm_sin_ps(alp));
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         cosa = _mm_cos_ps(alp);
                         x0   = _mm_sin_ps(_mm_mul_ps(
                                              _mm_add_ps(k0b,k0b),sina));
                         x1   = _mm_mul_ps(k0b,_mm_mul_ps(cosa,cosa));
                         sqr  = _mm_sub_ps(C10,_mm_div_ps(x0,x1));
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         x0   = _mm_sqrt_ps(sqr);
                         spin = _mm_mul_ps(_mm_sin_ps(pin),invn);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_add_ps(alp,alp),invn);
                         cos2a= _mm_cos_ps(x1);
                         t0   = _mm_mul_ps(C1772453850905516027298167483341,
                                                              _mm_mul_ps(b,x0)); // keep
                         x1   = _mm_sub_ps(cpin,cos2a);
                         t1   = _mm_mul_ps(_mm_add_ps(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm_mul_ps(a,spin)); // keep
                         x0   = _mm_rcp14_ps(x1); // keep
                         ear  = _mm_mul_ps(_mm_add_ps(t0,t1),x0);
                         t0   = _mm_mul_ps(ear,cer);
                         t1   = _mm_mul_ps(ear,cei);
                         rcs  = cabs_xmm4c4(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m128 rcs_f8193_xmm4r4_u(const float * __restrict  pb,
                                              const float * __restrict  pa,
                                              const float * __restrict  pk0,
                                              const float * __restrict  palp,
                                              const float * __restrict _pl) {
                                     
                          __m128 b  = _mm_loadu_ps(&pb[0]);
                          __m128 a  = _mm_loadu_ps(&pa[0]);  
                          __m128 k0 = _mm_loadu_ps(&pk0[0]);
                          __m128 alp= _mm_loadu_ps(&palp[0]);  
                          __m128 l  = _mm_loadu_ps(&pl[0]);   
                         const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);
                         const __m128 C15                              =
                                                     _mm_set1_ps(1.5f);
                          __m128 sina,pin,n,invn,x0,x1,k0l,sqr;
                          __m128 ear,eai,cer,cei,cosa;
                          __m128 rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm_mul_ps(k0,b);
                         ear  = _mm_setzero_ps();
                         n    = _mm_add_ps(C15,_mm_div_ps(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm_mul_ps(k0,l);
                         eai  = _mm_add_ps(k0l,k0l);
                         invn = _mm_rcp14_ps(n);
                         sina = _mm_sub_ps(C10,_mm_sin_ps(alp));
                         cexp_xmm4c4(ear,eai,&cer,&cei);
                         cosa = _mm_cos_ps(alp);
                         x0   = _mm_sin_ps(_mm_mul_ps(
                                              _mm_add_ps(k0b,k0b),sina));
                         x1   = _mm_mul_ps(k0b,_mm_mul_ps(cosa,cosa));
                         sqr  = _mm_sub_ps(C10,_mm_div_ps(x0,x1));
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         x0   = _mm_sqrt_ps(sqr);
                         spin = _mm_mul_ps(_mm_sin_ps(pin),invn);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_add_ps(alp,alp),invn);
                         cos2a= _mm_cos_ps(x1);
                         t0   = _mm_mul_ps(C1772453850905516027298167483341,
                                                              _mm_mul_ps(b,x0)); // keep
                         x1   = _mm_sub_ps(cpin,cos2a);
                         t1   = _mm_mul_ps(_mm_add_ps(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm_mul_ps(a,spin)); // keep
                         x0   = _mm_rcp14_ps(x1); // keep
                         ear  = _mm_mul_ps(_mm_add_ps(t0,t1),x0);
                         t0   = _mm_mul_ps(ear,cer);
                         t1   = _mm_mul_ps(ear,cei);
                         rcs  = cabs_xmm4c4(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for |theta| = PI/2-alpha
                     Formula 8.1-96
                 */
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m128 rcs_f8196_xmm4r4(const __m128 k0,
                                            const __m128 alp,
                                            const __m128 a,
                                            const __m128 b) {
                                            
                         
                          const __m128 C0444444444444444444444444444444 = 
                                                        _mm_set1_ps(0.444444444444444444444444444444f);
                           __m128 rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm_mul_ps(C0444444444444444444444444444444,k0);
                          cosa  = _mm_cos_ps(alp);
                          a32   = _mm_mul_ps(a,_mm_sqrt_ps(a));
                          x0    = _mm_mul_ps(b,cosa);
                          sina  = _mm_sin_ps(alp);
                          bca32 = _mm_mul_ps(x0,_mm_sqrt_ps(x0));
                          x1    = _mm_div_ps(cosa,sina);
                          cota  = _mm_mul_ps(x1,x1);
                          t1    = _mm_mul_ps(t0,_mm_mul_ps(cosa,cota));
                          x0    = _mm_sub_ps(a32,bca32);
                          rcs   = _mm_mul_ps(t1,_mm_mul_ps(x0,x0));
                          return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m128 rcs_f8196_xmm4r4_a(const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pa,
                                              const float * __restrict pb) {
                                            
                         
                           __m128 k0  = _mm_load_ps(&pk0[0]);
                           __m128 alp = _mm_load_ps(&palp[0]);
                           __m128 a   = _mm_load_ps(&pa[0]);
                           __m128 b   = _mm_load_ps(&pb[0]);
                          const __m128 C0444444444444444444444444444444 = 
                                                        _mm_set1_ps(0.444444444444444444444444444444f);
                           __m128 rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm_mul_ps(C0444444444444444444444444444444,k0);
                          cosa  = _mm_cos_ps(alp);
                          a32   = _mm_mul_ps(a,_mm_sqrt_ps(a));
                          x0    = _mm_mul_ps(b,cosa);
                          sina  = _mm_sin_ps(alp);
                          bca32 = _mm_mul_ps(x0,_mm_sqrt_ps(x0));
                          x1    = _mm_div_ps(cosa,sina);
                          cota  = _mm_mul_ps(x1,x1);
                          t1    = _mm_mul_ps(t0,_mm_mul_ps(cosa,cota));
                          x0    = _mm_sub_ps(a32,bca32);
                          rcs   = _mm_mul_ps(t1,_mm_mul_ps(x0,x0));
                          return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m128 rcs_f8196_xmm4r4_u(const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pa,
                                              const float * __restrict pb) {
                                            
                         
                           __m128 k0  = _mm_loadu_ps(&pk0[0]);
                           __m128 alp = _mm_loadu_ps(&palp[0]);
                           __m128 a   = _mm_loadu_ps(&pa[0]);
                           __m128 b   = _mm_loadu_ps(&pb[0]);
                          const __m128 C0444444444444444444444444444444 = 
                                                        _mm_set1_ps(0.444444444444444444444444444444f);
                           __m128 rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm_mul_ps(C0444444444444444444444444444444,k0);
                          cosa  = _mm_cos_ps(alp);
                          a32   = _mm_mul_ps(a,_mm_sqrt_ps(a));
                          x0    = _mm_mul_ps(b,cosa);
                          sina  = _mm_sin_ps(alp);
                          bca32 = _mm_mul_ps(x0,_mm_sqrt_ps(x0));
                          x1    = _mm_div_ps(cosa,sina);
                          cota  = _mm_mul_ps(x1,x1);
                          t1    = _mm_mul_ps(t0,_mm_mul_ps(cosa,cota));
                          x0    = _mm_sub_ps(a32,bca32);
                          rcs   = _mm_mul_ps(t1,_mm_mul_ps(x0,x0));
                          return (rcs);
                 }
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for 0<|theta|<alpha
                     Perpendicular RCS.
                     Formula 8.1-94
                 */
                 
                 
             /*      __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_perpendic_f8194_xmm4r4(const __m128 h,
	                                        const __m128 l,
	                                        const __m128 b,
	                                        const __m128 a,
	                                        const __m128 k0,
	                                        const __m128 tht,
	                                        const __m128 alp) {
	                                 
	                                  
	                 const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);  
                         const __m128 C15                              = 
                                                     _mm_set1_ps(1.5f); 
                         const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         const __m128 C20                              =
                                                     _mm_set1_ps(2.0f);
                          __m128 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m128 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m128 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m128 cpin1,cpin2,trm1,trm2,rcs;
                         __m128 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm_sub_ps(h,_mm_add_ps(l,b));
                         sint = _mm_sin_ps(tht);
                         k02  = _mm_add_ps(k0,k0);
                         n    = _mm_mul_ps(C15,_mm_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm_rcp14_ps(sint);
                         a2   = _mm_mul_ps(a,C05);
                         ear  = _mm_setzero_ps();
                         sk02 = _mm_sqrt_ps(_mm_mul_ps(k0,C05));
                         x0   = _mm_mul_ps(hlb,_mm_sub_ps(cost,b));
                         invn = _mm_rcp14_ps(n);
                         //x2   = _mm_mul_ps(a,C05);
                         eai  = _mm_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm_sqrt_ps(_mm_mul_ps(a2,csct));
                         atant= _mm_mul_ps(a,tant);
                         cost = _mm_cos_ps(tht);
                         x0   = _mm_mul_ps(b,C1772453850905516027298167483341);
                         cexp_xmm4c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm_mul_ps(x0,cer1);
                         spin = _mm_sin_ps(pin);
                         cei1 = _mm_mul_ps(x0,cei1);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_sub_ps(h,atant),cost);
                         eai2 = _mm_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_xmm4c4(ear,eai2,&cer2,&cei);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_mul_ps(x0,sacs);
                         cer2 = _mm_mul_ps(cer2,x1);
                         cei2 = _mm_mul_ps(cei2,x1);
                         cpin1= _mm_rcp14_ps(_mm_sub_ps(cpin,C10));
                         x2   = _mm_mul_ps(C20,_mm_add_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x2,invn));
                         x3   = _mm_rcp14_ps(_mm_sub_ps(cpin,cpin2));
                         trm1 = _mm_sub_ps(cpin1,x3);
                         cmul_xmm4c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm_mul_ps(t0r,trm1);
                         t0i  = _mm_mul_ps(t0i,trm1);
                         x0   = _mm_mul_ps(C20,_mm_sub_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x0,invn));
                         x1   = _mm_rcp14_ps(cpin2);
                         trm2 = _mm_sub_ps(cpin1,x1);
                         x2   = _mm_fmadd_ps(cost,_mm_mul_ps(k02,
                                                               _mm_add_ps(h,atant)));
                         eai3 = _mm_add_ps(C078539816339744830961566084582,x2);
                         cexp_xmm4c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_sqrt_ps(_mm_mul_ps(gms::math::
                                                                  negate_xmm4r4(a2),csct));
                         x2   = _mm_mul_ps(x0,x1);
                         cer3 = _mm_mul_ps(_mm_mul_ps(cer3,x2),trm2);
                         cei3 = _mm_mul_ps(_mm_mul_ps(cei3,x2),trm2);
                         cmul_xmm4c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_xmm4c4(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	                                  
                    __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_perpendic_f8194_xmm4r4_a(const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pb,
	                                         const float * __restrict pa,
	                                         const float * __restrict pk0,
	                                         const float * __restrict ptht,
	                                         const float * __restrict palp) {
	                                 
	                  
	                  __m128 h  = _mm_load_ps(&ph[0]);
	                  __m128 l  = _mm_load_ps(&pl[0]); 
	                  __m128 b  = _mm_load_ps(&pb[0]);   
	                  __m128 a  = _mm_load_ps(&pa[0]);  
	                  __m128 k0 = _mm_load_ps(&pk0[0]);
	                  __m128 tht= _mm_load_ps(&ptht[0]); 
	                  __m128 alp= _mm_load_ps(&palp[0]);        
	                 const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);  
                         const __m128 C15                              = 
                                                     _mm_set1_ps(1.5f); 
                         const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         const __m128 C20                              =
                                                     _mm_set1_ps(2.0f);
                         __m128 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m128 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m128 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m128 cpin1,cpin2,trm1,trm2,rcs;
                         __m128 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm_sub_ps(h,_mm_add_ps(l,b));
                         sint = _mm_sin_ps(tht);
                         k02  = _mm_add_ps(k0,k0);
                         n    = _mm_mul_ps(C15,_mm_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm_rcp14_ps(sint);
                         a2   = _mm_mul_ps(a,C05);
                         ear  = _mm_setzero_ps();
                         sk02 = _mm_sqrt_ps(_mm_mul_ps(k0,C05));
                         x0   = _mm_mul_ps(hlb,_mm_sub_ps(cost,b));
                         invn = _mm_rcp14_ps(n);
                         //x2   = _mm_mul_ps(a,C05);
                         eai  = _mm_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm_sqrt_ps(_mm_mul_ps(a2,csct));
                         atant= _mm_mul_ps(a,tant);
                         cost = _mm_cos_ps(tht);
                         x0   = _mm_mul_ps(b,C1772453850905516027298167483341);
                         cexp_xmm4c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm_mul_ps(x0,cer1);
                         spin = _mm_sin_ps(pin);
                         cei1 = _mm_mul_ps(x0,cei1);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_sub_ps(h,atant),cost);
                         eai2 = _mm_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_xmm4c4(ear,eai2,&cer2,&cei);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_mul_ps(x0,sacs);
                         cer2 = _mm_mul_ps(cer2,x1);
                         cei2 = _mm_mul_ps(cei2,x1);
                         cpin1= _mm_rcp14_ps(_mm_sub_ps(cpin,C10));
                         x2   = _mm_mul_ps(C20,_mm_add_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x2,invn));
                         x3   = _mm_rcp14_ps(_mm_sub_ps(cpin,cpin2));
                         trm1 = _mm_sub_ps(cpin1,x3);
                         cmul_xmm4c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm_mul_ps(t0r,trm1);
                         t0i  = _mm_mul_ps(t0i,trm1);
                         x0   = _mm_mul_ps(C20,_mm_sub_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x0,invn));
                         x1   = _mm_rcp14_ps(cpin2);
                         trm2 = _mm_sub_ps(cpin1,x1);
                         x2   = _mm_fmadd_ps(cost,_mm_mul_ps(k02,
                                                               _mm_add_ps(h,atant)));
                         eai3 = _mm_add_ps(C078539816339744830961566084582,x2);
                         cexp_xmm4c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_sqrt_ps(_mm_mul_ps(gms::math::
                                                                  negate_xmm4r4(a2),csct));
                         x2   = _mm_mul_ps(x0,x1);
                         cer3 = _mm_mul_ps(_mm_mul_ps(cer3,x2),trm2);
                         cei3 = _mm_mul_ps(_mm_mul_ps(cei3,x2),trm2);
                         cmul_xmm4c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_xmm4c4(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_perpendic_f8194_xmm4r4_u(    const float * __restrict  ph,
	                                                    const float * __restrict  pl,
	                                                    const float * __restrict  pb,
	                                                    const float * __restrict  pa,
	                                                    const float * __restrict  pk0,
	                                                    const float * __restrict  ptht,
	                                                    const float * __restrict  palp) {
	                                 
	                  
	                  __m128 h  = _mm_loadu_ps(&ph[0]);
	                  __m128 l  = _mm_loadu_ps(&pl[0]); 
	                  __m128 b  = _mm_loadu_ps(&pb[0]);   
	                  __m128 a  = _mm_loadu_ps(&pa[0]);  
	                  __m128 k0 = _mm_loadu_ps(&pk0[0]);
	                  __m128 tht= _mm_loadu_ps(&ptht[0]); 
	                  __m128 alp= _mm_loadu_ps(&palp[0]);        
	                 const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);  
                         const __m128 C15                              = 
                                                     _mm_set1_ps(1.5f); 
                         const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         const __m128 C20                              =
                                                     _mm_set1_ps(2.0f);
                         __m128 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m128 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m128 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m128 cpin1,cpin2,trm1,trm2,rcs;
                         __m128 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm_sub_ps(h,_mm_add_ps(l,b));
                         sint = _mm_sin_ps(tht);
                         k02  = _mm_add_ps(k0,k0);
                         n    = _mm_mul_ps(C15,_mm_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm_rcp14_ps(sint);
                         a2   = _mm_mul_ps(a,C05);
                         ear  = _mm_setzero_ps();
                         sk02 = _mm_sqrt_ps(_mm_mul_ps(k0,C05));
                         x0   = _mm_mul_ps(hlb,_mm_sub_ps(cost,b));
                         invn = _mm_rcp14_ps(n);
                         //x2   = _mm_mul_ps(a,C05);
                         eai  = _mm_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm_sqrt_ps(_mm_mul_ps(a2,csct));
                         atant= _mm_mul_ps(a,tant);
                         cost = _mm_cos_ps(tht);
                         x0   = _mm_mul_ps(b,C1772453850905516027298167483341);
                         cexp_xmm4c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm_mul_ps(x0,cer1);
                         spin = _mm_sin_ps(pin);
                         cei1 = _mm_mul_ps(x0,cei1);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_sub_ps(h,atant),cost);
                         eai2 = _mm_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_xmm4c4(ear,eai2,&cer2,&cei);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_mul_ps(x0,sacs);
                         cer2 = _mm_mul_ps(cer2,x1);
                         cei2 = _mm_mul_ps(cei2,x1);
                         cpin1= _mm_rcp14_ps(_mm_sub_ps(cpin,C10));
                         x2   = _mm_mul_ps(C20,_mm_add_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x2,invn));
                         x3   = _mm_rcp14_ps(_mm_sub_ps(cpin,cpin2));
                         trm1 = _mm_sub_ps(cpin1,x3);
                         cmul_xmm4c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm_mul_ps(t0r,trm1);
                         t0i  = _mm_mul_ps(t0i,trm1);
                         x0   = _mm_mul_ps(C20,_mm_sub_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x0,invn));
                         x1   = _mm_rcp14_ps(cpin2);
                         trm2 = _mm_sub_ps(cpin1,x1);
                         x2   = _mm_fmadd_ps(cost,_mm_mul_ps(k02,
                                                               _mm_add_ps(h,atant)));
                         eai3 = _mm_add_ps(C078539816339744830961566084582,x2);
                         cexp_xmm4c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_sqrt_ps(_mm_mul_ps(gms::math::
                                                                  negate_xmm4r4(a2),csct));
                         x2   = _mm_mul_ps(x0,x1);
                         cer3 = _mm_mul_ps(_mm_mul_ps(cer3,x2),trm2);
                         cei3 = _mm_mul_ps(_mm_mul_ps(cei3,x2),trm2);
                         cmul_xmm4c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_xmm4c4(t1r,t1i);
                         return (rcs);
	        }*/
	        
	        
	        
	         /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for 0<|theta|<alpha
                     Parallel RCS.
                     Formula 8.1-94
                 */
                 
                 
             /*      __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_parallel_f8194_xmm4r4(const __m128 h,
	                                             const __m128 l,
	                                             const __m128 b,
	                                             const __m128 a,
	                                             const __m128 k0,
	                                             const __m128 tht,
	                                             const __m128 alp) {
	                                 
	                                  
	                 const __m128 C314159265358979323846264338328  = 
                                                     _mm_set1_ps(3.14159265358979323846264338328f); 
                         const __m128 C1772453850905516027298167483341 = 
                                                     _mm_set1_ps(1.772453850905516027298167483341f);
                         const __m128 C078539816339744830961566084582  = 
                                                     _mm_set1_ps(0.78539816339744830961566084582f);
                         const __m128 C10                              = 
                                                     _mm_set1_ps(1.0f);  
                         const __m128 C15                              = 
                                                     _mm_set1_ps(1.5f); 
                         const __m128 C05                              = 
                                                     _mm_set1_ps(0.5f);
                         const __m128 C20                              =
                                                     _mm_set1_ps(2.0f);
                          __m128 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m128 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m128 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m128 cpin1,cpin2,trm1,trm2,rcs;
                         __m128 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm_sub_ps(h,_mm_add_ps(l,b));
                         sint = _mm_sin_ps(tht);
                         k02  = _mm_add_ps(k0,k0);
                         n    = _mm_mul_ps(C15,_mm_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm_rcp14_ps(sint);
                         a2   = _mm_mul_ps(a,C05);
                         ear  = _mm_setzero_ps();
                         sk02 = _mm_sqrt_ps(_mm_mul_ps(k0,C05));
                         x0   = _mm_mul_ps(hlb,_mm_sub_ps(cost,b));
                         invn = _mm_rcp14_ps(n);
                         //x2   = _mm_mul_ps(a,C05);
                         eai  = _mm_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm_sqrt_ps(_mm_mul_ps(a2,csct));
                         atant= _mm_mul_ps(a,tant);
                         cost = _mm_cos_ps(tht);
                         x0   = _mm_mul_ps(b,C1772453850905516027298167483341);
                         cexp_xmm4c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm_mul_ps(x0,cer1);
                         spin = _mm_sin_ps(pin);
                         cei1 = _mm_mul_ps(x0,cei1);
                         cpin = _mm_cos_ps(pin);
                         x1   = _mm_mul_ps(_mm_sub_ps(h,atant),cost);
                         eai2 = _mm_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_xmm4c4(ear,eai2,&cer2,&cei);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_mul_ps(x0,sacs);
                         cer2 = _mm_mul_ps(cer2,x1);
                         cei2 = _mm_mul_ps(cei2,x1);
                         cpin1= _mm_rcp14_ps(_mm_sub_ps(cpin,C10));
                         x2   = _mm_mul_ps(C20,_mm_add_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x2,invn));
                         x3   = _mm_rcp14_ps(_mm_sub_ps(cpin,cpin2));
                         trm1 = _mm_add_ps(cpin1,x3);
                         cmul_xmm4c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm_mul_ps(t0r,trm1);
                         t0i  = _mm_mul_ps(t0i,trm1);
                         x0   = _mm_mul_ps(C20,_mm_sub_ps(alp,tht));
                         cpin2= _mm_cos_ps(_mm_mul_ps(x0,invn));
                         x1   = _mm_rcp14_ps(cpin2);
                         trm2 = _mm_add_ps(cpin1,x1);
                         x2   = _mm_fmadd_ps(cost,_mm_mul_ps(k02,
                                                               _mm_add_ps(h,atant)));
                         eai3 = _mm_add_ps(C078539816339744830961566084582,x2);
                         cexp_xmm4c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm_div_ps(spin,_mm_mul_ps(n,sk02));
                         x1   = _mm_sqrt_ps(_mm_mul_ps(gms::math::
                                                                  negate_xmm4r4(a2),csct));
                         x2   = _mm_mul_ps(x0,x1);
                         cer3 = _mm_mul_ps(_mm_mul_ps(cer3,x2),trm2);
                         cei3 = _mm_mul_ps(_mm_mul_ps(cei3,x2),trm2);
                         cmul_xmm4c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_xmm4c4(t1r,t1i);
                         return (rcs);
	        }
	        */
	        
	        
	        /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Helper formula coefficient B(gamma).
	             Formula 9.1-37
	        */
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 coef_Bg_f9137_xmm4r4(const __m128 A,
	                                        const __m128 N,
	                                        const __m128 k0,
	                                        const __m128 epsr,
	                                        const __m128 epsi,
	                                        const __m128 thti,
	                                        const __m128 thts,
	                                        const __m128 phis,
	                                        const int pol) {
	                                        
	                  const __m128 C06          = _mm_set1_ps(0.6f);
	                  const __m128 C01875       = _mm_set1_ps(0.1875f); // 3/16
	                  const __m128 C87964594300514210676954014731826 = 
	                                            = _mm_set1_ps(87.964594300514210676954014731826f);
	                  const __m128 C20          = _mm_set1_ps(2.0f);
	                  const __m128 C10          = _mm_set1_ps(1.0f);
	                  const __m128 C30          = _mm_set1_ps(3.0f);
	                   __m128 Bg,t,t2,cthti,cthts,sthti,sthts,cphis;
	                   __m128 AN,Ak0,num,den,alp,x0,x1,x2,x3,x4,secti;
	                  AN    = _mm_mul_ps(A,N);
	                  t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                  x1    = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  x2    = _mm_fmadd_ps(epsi,epsi,
	                                               _mm_mul_ps(x1,x1));
	                  t2    = _mm_mul_ps(t,t);
	                  secti = _mm_rcp14_ps(cthti);
	                  sthti = _mm_sin_ps(thti);
	                  x0    = _mm_fmadd_ps(t2,C30,C10);
	                  cthts = _mm_cos_ps(thts);
	                  Ak0   = _mm_mul_ps(A,k0);
	                  x1    = _mm_mul_ps(AN,_mm_mul_ps(Ak0,AK0));
	                  num   = _mm_mul_ps(x1,x2);
	                  if(pol == 0) {
	                     x1 = _mm_mul_ps(_mm_mul_ps(C01875,AN),
	                                                  _mm_mul_ps(epsi,secti));
	                     x2 = _mm_mul_ps(x1,x0);
	                     alp= _mm_mul_ps(k0,x2);
	                  }
	                  else if(pol == 1) {
	                     x1 = _mm_mul_ps(_mm_mul_ps(C01875,AN),
	                                                  _mm_mul_ps(epsi,secti));
	                     x2 = _mm_mul_ps(sthti,sthti);
	                     x3 = _mm_sub_ps(c10,t2);
	                     alp= _mm_fmadd_ps(x1,x0,_mm_mul_ps(x2,x3));
	                     alp= _mm_mul_ps(k0,alp);
	                  }
	                  x0 = _mm_div_ps(alp,k0);
	                  x2 = _mm_add_ps(cthti,cthts);
	                  x1 = _mm_fmadd_ps(C06,_mm_mul_ps(x0,x0),
	                                                     _mm_mul_ps(C30,
	                                                                   _mm_mul_ps(x2,x2)));
	                  x3 = _mm_fmadd_ps(sthti,sthti,_mm_mul_ps(sthts,sthts));
	                  x4 = _mm_mul_ps(C20,_mm_mul_ps(sthti,sthts));
	                  x0 = _mm_mul_ps(x4,cthts);
	                  den= _mm_mul_ps(C87964594300514210676954014731826,
	                                                             _mm_add_ps(x1,x0));
	                  
	                  Bg = _mm_div_ps(num,den);
	                  return (Bg);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 coef_Bg_f9137_xmm4r4_a(const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol) {
	                         
	                   __m128 A   = _mm_load_ps(&pA[0]);
	                   __m128 N   = _mm_load_ps(&pN[0]); 
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 epsr= _mm_load_ps(&pepsr[0]);
	                   __m128 epsi= _mm_load_ps(&pepsi[0]);
	                   __m128 thti= _mm_load_ps(&pthti[0]);
	                   __m128 thts= _mm_load_ps(&pthts[0]);
	                   __m128 phis= _mm_load_ps(&pphis[0]);              
	                  const __m128 C06          = _mm_set1_ps(0.6f);
	                  const __m128 C01875       = _mm_set1_ps(0.1875f); // 3/16
	                  const __m128 C87964594300514210676954014731826 = 
	                                            = _mm_set1_ps(87.964594300514210676954014731826f);
	                  const __m128 C20          = _mm_set1_ps(2.0f);
	                  const __m128 C10          = _mm_set1_ps(1.0f);
	                  const __m128 C30          = _mm_set1_ps(3.0f);
	                   __m128 Bg,t,t2,cthti,cthts,sthti,sthts,cphis;
	                   __m128 AN,Ak0,num,den,alp,x0,x1,x2,x3,x4,secti;
	                  AN    = _mm_mul_ps(A,N);
	                  t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                  x1    = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  x2    = _mm_fmadd_ps(epsi,epsi,
	                                               _mm_mul_ps(x1,x1));
	                  t2    = _mm_mul_ps(t,t);
	                  secti = _mm_rcp14_ps(cthti);
	                  sthti = _mm_sin_ps(thti);
	                  x0    = _mm_fmadd_ps(t2,C30,C10);
	                  cthts = _mm_cos_ps(thts);
	                  Ak0   = _mm_mul_ps(A,k0);
	                  x1    = _mm_mul_ps(AN,_mm_mul_ps(Ak0,AK0));
	                  num   = _mm_mul_ps(x1,x2);
	                  if(pol == 0) {
	                     x1 = _mm_mul_ps(_mm_mul_ps(C01875,AN),
	                                                  _mm_mul_ps(epsi,secti));
	                     x2 = _mm_mul_ps(x1,x0);
	                     alp= _mm_mul_ps(k0,x2);
	                  }
	                  else if(pol == 1) {
	                     x1 = _mm_mul_ps(_mm_mul_ps(C01875,AN),
	                                                  _mm_mul_ps(epsi,secti));
	                     x2 = _mm_mul_ps(sthti,sthti);
	                     x3 = _mm_sub_ps(c10,t2);
	                     alp= _mm_fmadd_ps(x1,x0,_mm_mul_ps(x2,x3));
	                     alp= _mm_mul_ps(k0,alp);
	                  }
	                  x0 = _mm_div_ps(alp,k0);
	                  x2 = _mm_add_ps(cthti,cthts);
	                  x1 = _mm_fmadd_ps(C06,_mm_mul_ps(x0,x0),
	                                                     _mm_mul_ps(C30,
	                                                                   _mm_mul_ps(x2,x2)));
	                  x3 = _mm_fmadd_ps(sthti,sthti,_mm_mul_ps(sthts,sthts));
	                  x4 = _mm_mul_ps(C20,_mm_mul_ps(sthti,sthts));
	                  x0 = _mm_mul_ps(x4,cthts);
	                  den= _mm_mul_ps(C87964594300514210676954014731826,
	                                                             _mm_add_ps(x1,x0));
	                  
	                  Bg = _mm_div_ps(num,den);
	                  return (Bg);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 coef_Bg_f9137_xmm4r4_u(const float * __restrict  pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol) {
	                         
	                   __m128 A   = _mm_loadu_ps(&pA[0]);
	                   __m128 N   = _mm_loadu_ps(&pN[0]); 
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                   __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                   __m128 thti= _mm_loadu_ps(&pthti[0]);
	                   __m128 thts= _mm_loadu_ps(&pthts[0]);
	                   __m128 phis= _mm_loadu_ps(&pphis[0]);              
	                  const __m128 C06          = _mm_set1_ps(0.6f);
	                  const __m128 C01875       = _mm_set1_ps(0.1875f); // 3/16
	                  const __m128 C87964594300514210676954014731826 = 
	                                            = _mm_set1_ps(87.964594300514210676954014731826f);
	                  const __m128 C20          = _mm_set1_ps(2.0f);
	                  const __m128 C10          = _mm_set1_ps(1.0f);
	                  const __m128 C30          = _mm_set1_ps(3.0f);
	                   __m128 Bg,t,t2,cthti,cthts,sthti,sthts,cphis;
	                   __m128 AN,Ak0,num,den,alp,x0,x1,x2,x3,x4,secti;
	                  AN    = _mm_mul_ps(A,N);
	                  t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                  x1    = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  x2    = _mm_fmadd_ps(epsi,epsi,
	                                               _mm_mul_ps(x1,x1));
	                  t2    = _mm_mul_ps(t,t);
	                  secti = _mm_rcp14_ps(cthti);
	                  sthti = _mm_sin_ps(thti);
	                  x0    = _mm_fmadd_ps(t2,C30,C10);
	                  cthts = _mm_cos_ps(thts);
	                  Ak0   = _mm_mul_ps(A,k0);
	                  x1    = _mm_mul_ps(AN,_mm_mul_ps(Ak0,AK0));
	                  num   = _mm_mul_ps(x1,x2);
	                  if(pol == 0) {
	                     x1 = _mm_mul_ps(_mm_mul_ps(C01875,AN),
	                                                  _mm_mul_ps(epsi,secti));
	                     x2 = _mm_mul_ps(x1,x0);
	                     alp= _mm_mul_ps(k0,x2);
	                  }
	                  else if(pol == 1) {
	                     x1 = _mm_mul_ps(_mm_mul_ps(C01875,AN),
	                                                  _mm_mul_ps(epsi,secti));
	                     x2 = _mm_mul_ps(sthti,sthti);
	                     x3 = _mm_sub_ps(c10,t2);
	                     alp= _mm_fmadd_ps(x1,x0,_mm_mul_ps(x2,x3));
	                     alp= _mm_mul_ps(k0,alp);
	                  }
	                  x0 = _mm_div_ps(alp,k0);
	                  x2 = _mm_add_ps(cthti,cthts);
	                  x1 = _mm_fmadd_ps(C06,_mm_mul_ps(x0,x0),
	                                                     _mm_mul_ps(C30,
	                                                                   _mm_mul_ps(x2,x2)));
	                  x3 = _mm_fmadd_ps(sthti,sthti,_mm_mul_ps(sthts,sthts));
	                  x4 = _mm_mul_ps(C20,_mm_mul_ps(sthti,sthts));
	                  x0 = _mm_mul_ps(x4,cthts);
	                  den= _mm_mul_ps(C87964594300514210676954014731826,
	                                                             _mm_add_ps(x1,x0));
	                  
	                  Bg = _mm_div_ps(num,den);
	                  return (Bg);
	       }
	       
	       
	       /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (hh) polarized per unit surface area.
	             Formula 9.1-33
	       
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9133_xmm4r4( const __m128 A,
	                                        const __m128 N,
	                                        const __m128 k0,
	                                        const __m128 epsr,
	                                        const __m128 epsi,
	                                        const __m128 thti,
	                                        const __m128 thts,
	                                        const __m128 phis,
	                                        const int pol) {
	                    
	                  const __m128 C10 = _mm_set1_ps(1.0f);                    
	                  const __m128 C30 = _mm_set1_ps(3.0f);
	                  const __m128 C20 = _mm_set1_ps(2.0f);
	                  const __m128 C80 = _mm_set1_ps(8.0f);
	                  const __m128 C100= _mm_set1_ps(10.0f);
	                  const __m128 C240= _mm_set1_ps(24.0f);
	                  const __m128 C230= _mm_set1_ps(23.0f);
	                   rcs,Bg,sphis,x0,t,t2,trm1,trm2,trm3;
	                  t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                  x0    = _mm_sin_ps(phis);
	                  Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                  t2    = _mm_mul_ps(t,t);
	                  sphis = _mm_mul_ps(x0,x0); 
	                  trm1  = _mm_sub_ps(C30,_mm_mul_ps(C20,sphis));
	                  trm2  = _mm_mul_ps(t,_mm_sub_ps(C80,
	                                                      _mm_mul_ps(C100,sphis)));
	                  trm3  = _mm_mul_ps(t2,_mm_sub_ps(C240,
	                                                      _mm_mul_ps(C230,sphis)));
	                  x0    = _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                  rcs   = _mm_mul_ps(Bg,x0);
	                  return (rcs);                             
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9133_xmm4r4_a( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol) {
	                         
	                   __m128 A   = _mm_load_ps(&pA[0]);
	                   __m128 N   = _mm_load_ps(&pN[0]); 
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 epsr= _mm_load_ps(&pepsr[0]);
	                   __m128 epsi= _mm_load_ps(&pepsi[0]);
	                   __m128 thti= _mm_load_ps(&pthti[0]);
	                   __m128 thts= _mm_load_ps(&pthts[0]);
	                   __m128 phis= _mm_load_ps(&pphis[0]); 
	                    
	                  const __m128 C10 = _mm_set1_ps(1.0f);                    
	                  const __m128 C30 = _mm_set1_ps(3.0f);
	                  const __m128 C20 = _mm_set1_ps(2.0f);
	                  const __m128 C80 = _mm_set1_ps(8.0f);
	                  const __m128 C100= _mm_set1_ps(10.0f);
	                  const __m128 C240= _mm_set1_ps(24.0f);
	                  const __m128 C230= _mm_set1_ps(23.0f);
	                   rcs,Bg,sphis,x0,t,t2,trm1,trm2,trm3;
	                  t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                  x0    = _mm_sin_ps(phis);
	                  Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                  t2    = _mm_mul_ps(t,t);
	                  sphis = _mm_mul_ps(x0,x0); 
	                  trm1  = _mm_sub_ps(C30,_mm_mul_ps(C20,sphis));
	                  trm2  = _mm_mul_ps(t,_mm_sub_ps(C80,
	                                                      _mm_mul_ps(C100,sphis)));
	                  trm3  = _mm_mul_ps(t2,_mm_sub_ps(C240,
	                                                      _mm_mul_ps(C230,sphis)));
	                  x0    = _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                  rcs   = _mm_mul_ps(Bg,x0);
	                  return (rcs);                             
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9133_xmm4r4_u( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol) {
	                         
	                   __m128 A   = _mm_loadu_ps(&pA[0]);
	                   __m128 N   = _mm_loadu_ps(&pN[0]); 
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                   __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                   __m128 thti= _mm_loadu_ps(&pthti[0]);
	                   __m128 thts= _mm_loadu_ps(&pthts[0]);
	                   __m128 phis= _mm_loadu_ps(&pphis[0]); 
	                    
	                  const __m128 C10 = _mm_set1_ps(1.0f);                    
	                  const __m128 C30 = _mm_set1_ps(3.0f);
	                  const __m128 C20 = _mm_set1_ps(2.0f);
	                  const __m128 C80 = _mm_set1_ps(8.0f);
	                  const __m128 C100= _mm_set1_ps(10.0f);
	                  const __m128 C240= _mm_set1_ps(24.0f);
	                  const __m128 C230= _mm_set1_ps(23.0f);
	                   rcs,Bg,sphis,x0,t,t2,trm1,trm2,trm3;
	                  t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                  x0    = _mm_sin_ps(phis);
	                  Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                  t2    = _mm_mul_ps(t,t);
	                  sphis = _mm_mul_ps(x0,x0); 
	                  trm1  = _mm_sub_ps(C30,_mm_mul_ps(C20,sphis));
	                  trm2  = _mm_mul_ps(t,_mm_sub_ps(C80,
	                                                      _mm_mul_ps(C100,sphis)));
	                  trm3  = _mm_mul_ps(t2,_mm_sub_ps(C240,
	                                                      _mm_mul_ps(C230,sphis)));
	                  x0    = _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                  rcs   = _mm_mul_ps(Bg,x0);
	                  return (rcs);                             
	       }
	       
	       
	        
	       /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (vh) polarized per unit surface area.
	             Formula 9.1-34
	       
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vh_f9134_xmm4r4( const __m128 A,
	                                        const __m128 N,
	                                        const __m128 k0,
	                                        const __m128 epsr,
	                                        const __m128 epsi,
	                                        const __m128 thti,
	                                        const __m128 thts,
	                                        const __m128 phis,
	                                        const int pol) {
	                     
	                    const __m128 C10 = _mm_set1_ps(1.0f);                             
	                    const __m128 C100= _mm_set1_ps(10.0f);
	                    const __m128 C20 = _mm_set1_ps(2.0f);   
	                    const __m128 C230= _mm_set1_ps(23.0f); 
	                     __m128 rcs,trm1,trm2,trm3,trm4;
	                     __m128 Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                    x0    = _mm_sin_ps(phis);
	                    Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm_mul_ps(t,t);  
	                    x1    = _mm_cos_ps(thts);
	                    sphis = _mm_mul_ps(x0,x0); 
	                    cthts = _mm_mul_ps(x1,x1);
	                    x0    = _mm_sin_ps(thts);
	                    sthts = _mm_mul_ps(x0,x0);
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);  
	                    trm1= _mm_fmadd_ps(sthts,_mm_mul_ps(C20,cthts),x1);
	                    x0  = _mm_sub_ps(gms::math::negate_xmm4r4(C20),
	                                                          _mm_mul_ps(C40,sthts));
	                    trm2= _mm_mul_ps(t,
	                                     _mm_fmadd_ps(sphis,
	                                                    _mm_mul_ps(C100,cthts),x0));
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);
	                    trm3= _mm_mul_ps(t2,
	                                      _mm_fmadd_ps(sphis,
	                                                     _mm_mul_ps(cthts,C230),x1));
	                    trm4= _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                    rcs = _mm_mul_ps(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vh_f9134_xmm4r4_a( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol) {
	                         
	                     __m128 A   = _mm_load_ps(&pA[0]);
	                     __m128 N   = _mm_load_ps(&pN[0]); 
	                     __m128 k0  = _mm_load_ps(&pk0[0]);
	                     __m128 epsr= _mm_load_ps(&pepsr[0]);
	                     __m128 epsi= _mm_load_ps(&pepsi[0]);
	                     __m128 thti= _mm_load_ps(&pthti[0]);
	                     __m128 thts= _mm_load_ps(&pthts[0]);
	                     __m128 phis= _mm_load_ps(&pphis[0]); 
	                     
	                    const __m128 C10 = _mm_set1_ps(1.0f);                             
	                    const __m128 C100= _mm_set1_ps(10.0f);
	                    const __m128 C20 = _mm_set1_ps(2.0f);   
	                    const __m128 C230= _mm_set1_ps(23.0f); 
	                     __m128 rcs,trm1,trm2,trm3,trm4;
	                     __m128 Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                    x0    = _mm_sin_ps(phis);
	                    Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm_mul_ps(t,t);  
	                    x1    = _mm_cos_ps(thts);
	                    sphis = _mm_mul_ps(x0,x0); 
	                    cthts = _mm_mul_ps(x1,x1);
	                    x0    = _mm_sin_ps(thts);
	                    sthts = _mm_mul_ps(x0,x0);
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);  
	                    trm1= _mm_fmadd_ps(sthts,_mm_mul_ps(C20,cthts),x1);
	                    x0  = _mm_sub_ps(gms::math::negate_xmm4r4(C20),
	                                                          _mm_mul_ps(C40,sthts));
	                    trm2= _mm_mul_ps(t,
	                                     _mm_fmadd_ps(sphis,
	                                                    _mm_mul_ps(C100,cthts),x0));
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);
	                    trm3= _mm_mul_ps(t2,
	                                      _mm_fmadd_ps(sphis,
	                                                     _mm_mul_ps(cthts,C230),x1));
	                    trm4= _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                    rcs = _mm_mul_ps(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vh_f9134_xmm4r4_u( const float * __restrict  pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol) {
	                         
	                     __m128 A   = _mm_loadu_ps(&pA[0]);
	                     __m128 N   = _mm_loadu_ps(&pN[0]); 
	                     __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                     __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                     __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                     __m128 thti= _mm_loadu_ps(&pthti[0]);
	                     __m128 thts= _mm_loadu_ps(&pthts[0]);
	                     __m128 phis= _mm_loadu_ps(&pphis[0]); 
	                     
	                    const __m128 C10 = _mm_set1_ps(1.0f);                             
	                    const __m128 C100= _mm_set1_ps(10.0f);
	                    const __m128 C20 = _mm_set1_ps(2.0f);   
	                    const __m128 C230= _mm_set1_ps(23.0f); 
	                     __m128 rcs,trm1,trm2,trm3,trm4;
	                     __m128 Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                    x0    = _mm_sin_ps(phis);
	                    Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm_mul_ps(t,t);  
	                    x1    = _mm_cos_ps(thts);
	                    sphis = _mm_mul_ps(x0,x0); 
	                    cthts = _mm_mul_ps(x1,x1);
	                    x0    = _mm_sin_ps(thts);
	                    sthts = _mm_mul_ps(x0,x0);
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);  
	                    trm1= _mm_fmadd_ps(sthts,_mm_mul_ps(C20,cthts),x1);
	                    x0  = _mm_sub_ps(gms::math::negate_xmm4r4(C20),
	                                                          _mm_mul_ps(C40,sthts));
	                    trm2= _mm_mul_ps(t,
	                                     _mm_fmadd_ps(sphis,
	                                                    _mm_mul_ps(C100,cthts),x0));
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);
	                    trm3= _mm_mul_ps(t2,
	                                      _mm_fmadd_ps(sphis,
	                                                     _mm_mul_ps(cthts,C230),x1));
	                    trm4= _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                    rcs = _mm_mul_ps(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	           /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (hv) polarized per unit surface area.
	             Formula 9.1-35
	       
	       */
	       
	       
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hv_f9135_xmm4r4( const __m128 A,
	                                        const __m128 N,
	                                        const __m128 k0,
	                                        const __m128 epsr,
	                                        const __m128 epsi,
	                                        const __m128 thti,
	                                        const __m128 thts,
	                                        const __m128 phis,
	                                        const int pol) {
	                     
	                    const __m128 C10 = _mm_set1_ps(1.0f);                             
	                    const __m128 C100= _mm_set1_ps(10.0f);
	                    const __m128 C20 = _mm_set1_ps(2.0f);   
	                    const __m128 C230= _mm_set1_ps(23.0f); 
	                     __m128 rcs,trm1,trm2,trm3,trm4;
	                     __m128 Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                    x0    = _mm_sin_ps(phis);
	                    Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm_mul_ps(t,t);  
	                    x1    = _mm_cos_ps(thti);
	                    sphis = _mm_mul_ps(x0,x0); 
	                    cthts = _mm_mul_ps(x1,x1);
	                    x0    = _mm_sin_ps(thti);
	                    sthts = _mm_mul_ps(x0,x0);
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);  
	                    trm1= _mm_fmadd_ps(sthts,_mm_mul_ps(C20,cthts),x1);
	                    x0  = _mm_sub_ps(gms::math::negate_xmm4r4(C20),
	                                                          _mm_mul_ps(C40,sthts));
	                    trm2= _mm_mul_ps(t,
	                                     _mm_fmadd_ps(sphis,
	                                                    _mm_mul_ps(C100,cthts),x0));
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);
	                    trm3= _mm_mul_ps(t2,
	                                      _mm_fmadd_ps(sphis,
	                                                     _mm_mul_ps(cthts,C230),x1));
	                    trm4= _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                    rcs = _mm_mul_ps(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hv_f9135_xmm4r4_a(  const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol) {
	                         
	                     __m128 A   = _mm_load_ps(&pA[0]);
	                     __m128 N   = _mm_load_ps(&pN[0]); 
	                     __m128 k0  = _mm_load_ps(&pk0[0]);
	                     __m128 epsr= _mm_load_ps(&pepsr[0]);
	                     __m128 epsi= _mm_load_ps(&pepsi[0]);
	                     __m128 thti= _mm_load_ps(&pthti[0]);
	                     __m128 thts= _mm_load_ps(&pthts[0]);
	                     __m128 phis= _mm_load_ps(&pphis[0]); 
	                     
	                    const __m128 C10 = _mm_set1_ps(1.0f);                             
	                    const __m128 C100= _mm_set1_ps(10.0f);
	                    const __m128 C20 = _mm_set1_ps(2.0f);   
	                    const __m128 C230= _mm_set1_ps(23.0f); 
	                     __m128 rcs,trm1,trm2,trm3,trm4;
	                     __m128 Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                    x0    = _mm_sin_ps(phis);
	                    Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm_mul_ps(t,t);  
	                    x1    = _mm_cos_ps(thti);
	                    sphis = _mm_mul_ps(x0,x0); 
	                    cthts = _mm_mul_ps(x1,x1);
	                    x0    = _mm_sin_ps(thti);
	                    sthts = _mm_mul_ps(x0,x0);
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);  
	                    trm1= _mm_fmadd_ps(sthts,_mm_mul_ps(C20,cthts),x1);
	                    x0  = _mm_sub_ps(gms::math::negate_xmm4r4(C20),
	                                                          _mm_mul_ps(C40,sthts));
	                    trm2= _mm_mul_ps(t,
	                                     _mm_fmadd_ps(sphis,
	                                                    _mm_mul_ps(C100,cthts),x0));
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);
	                    trm3= _mm_mul_ps(t2,
	                                      _mm_fmadd_ps(sphis,
	                                                     _mm_mul_ps(cthts,C230),x1));
	                    trm4= _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                    rcs = _mm_mul_ps(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128  rcs_hv_f9135_xmm4r4_u(const float * __restrict pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol) {
	                         
	                     __m128 A   = _mm_loadu_ps(&pA[0]);
	                     __m128 N   = _mm_loadu_ps(&pN[0]); 
	                     __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                     __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                     __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                     __m128 thti= _mm_loadu_ps(&pthti[0]);
	                     __m128 thts= _mm_loadu_ps(&pthts[0]);
	                     __m128 phis= _mm_loadu_ps(&pphis[0]); 
	                     
	                    const __m128 C10 = _mm_set1_ps(1.0f);                             
	                    const __m128 C100= _mm_set1_ps(10.0f);
	                    const __m128 C20 = _mm_set1_ps(2.0f);   
	                    const __m128 C230= _mm_set1_ps(23.0f); 
	                     __m128 rcs,trm1,trm2,trm3,trm4;
	                     __m128 Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));
	                    x0    = _mm_sin_ps(phis);
	                    Bg    = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm_mul_ps(t,t);  
	                    x1    = _mm_cos_ps(thti);
	                    sphis = _mm_mul_ps(x0,x0); 
	                    cthts = _mm_mul_ps(x1,x1);
	                    x0    = _mm_sin_ps(thti);
	                    sthts = _mm_mul_ps(x0,x0);
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);  
	                    trm1= _mm_fmadd_ps(sthts,_mm_mul_ps(C20,cthts),x1);
	                    x0  = _mm_sub_ps(gms::math::negate_xmm4r4(C20),
	                                                          _mm_mul_ps(C40,sthts));
	                    trm2= _mm_mul_ps(t,
	                                     _mm_fmadd_ps(sphis,
	                                                    _mm_mul_ps(C100,cthts),x0));
	                    x1  = _mm_fmadd_ps(sthts,C20,C10);
	                    trm3= _mm_mul_ps(t2,
	                                      _mm_fmadd_ps(sphis,
	                                                     _mm_mul_ps(cthts,C230),x1));
	                    trm4= _mm_add_ps(trm1,_mm_add_ps(trm2,trm3));
	                    rcs = _mm_mul_ps(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	            /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (vv) polarized per unit surface area.
	             Formula 9.1-36
	       
	       */
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9136_xmm4r4( const __m128 A,
	                                        const __m128 N,
	                                        const __m128 k0,
	                                        const __m128 epsr,
	                                        const __m128 epsi,
	                                        const __m128 thti,
	                                        const __m128 thts,
	                                        const __m128 phis,
	                                        const int pol) {
	                                        
	                const __m128 C10 = _mm_set1_ps(1.0f);     
	                const __m128 C20 = _mm_set1_ps(2.0f);   
	                const __m128 C30 = _mm_set1_ps(3.0f);
	                const __m128 C120= _mm_set1_ps(12.0f);
	                const __m128 C140= _mm_set1_ps(14.0f);
	                const __m128 C350= _mm_set1_ps(35.0f);
	                 __m128 t,t2,x0,x1,x2,x3,cthti,cthts;
	                 __m128 sphis,cphis,sthts,cthts;
	                 __m128 trm1,trm2,trm3,trm4;
	                 __m128 rcs,Bg,cterm,sterm,sctrm;
	                x0     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));   
	                t      = _mm_sub_ps(C10,x0);
	                t2     = _mm_sub_ps(C10,_mm_mul_ps(x0,x0));
	                Bg     = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol); 
	                sphis  = _mm_sin_ps(phis);
	                cphis  = _mm_cos_ps(phis);
	                cthti  = _mm_cos_ps(thti);
	                sthti  = _mm_sin_ps(thti);
	                cthts  = _mm_cos_ps(thts);
	                cterm  = _mm_mul_ps(cthti,_mm_mul_ps(cthts,cphis));
	                sthts  = _mm_sin_ps(thts);
	                sterm  = _mm_mul_ps(sthti,sthts);
	                sctrm  = _mm_sub_ps(sterm,cterm);
	                x0     = _mm_mul_ps(sphis,sphis);
	                x1     = _mm_mul_ps(cthts,cthts);
	                x2     = _mm_mul_ps(C30,_mm_mul_ps(cthti,cthti));
	                x3     = _mm_sub_ps(C30,_mm_mul_ps(x0,
	                                                    _mm_mul_ps(x1,x2)));
	                trm1   = _mm_fmadd_ps(t2,x3,sctrm);
	                trm2   = _mm_mul_ps(C120,_mm_mul_ps(t2,sterm));
	                x0     = _mm_fmsub_ps(C30,sterm,cterm);
	                trm3   = _mm_mul_ps(C140,_mm_mul_ps(t,x0));
	                trm4   = _mm_mul_ps(_mm_mul_ps(C350,t2),sctrm);
	                x1     = _mm_mul_ps(trm1,trm2);
	                x2     = _mm_add_ps(x1,_mm_add_ps(trm3,trm4));
	                rcs    = _mm_mul_ps(Bg,x2);
	                return (rcs);
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9136_xmm4r4_a( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol ) {
	                                      
	                 __m128 A   = _mm_load_ps(&pA[0]);
	                 __m128 N   = _mm_load_ps(&pN[0]); 
	                 __m128 k0  = _mm_load_ps(&pk0[0]);
	                 __m128 epsr= _mm_load_ps(&pepsr[0]);
	                 __m128 epsi= _mm_load_ps(&pepsi[0]);
	                 __m128 thti= _mm_load_ps(&pthti[0]);
	                 __m128 thts= _mm_load_ps(&pthts[0]);
	                 __m128 phis= _mm_load_ps(&pphis[0]);   
	                const __m128 C10 = _mm_set1_ps(1.0f);     
	                const __m128 C20 = _mm_set1_ps(2.0f);   
	                const __m128 C30 = _mm_set1_ps(3.0f);
	                const __m128 C120= _mm_set1_ps(12.0f);
	                const __m128 C140= _mm_set1_ps(14.0f);
	                const __m128 C350= _mm_set1_ps(35.0f);
	                 __m128 t,t2,x0,x1,x2,x3,cthti,cthts;
	                 __m128 sphis,cphis,sthts,cthts;
	                 __m128 trm1,trm2,trm3,trm4;
	                 __m128 rcs,Bg,cterm,sterm,sctrm;
	                x0     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));   
	                t      = _mm_sub_ps(C10,x0);
	                t2     = _mm_sub_ps(C10,_mm_mul_ps(x0,x0));
	                Bg     = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol); 
	                sphis  = _mm_sin_ps(phis);
	                cphis  = _mm_cos_ps(phis);
	                cthti  = _mm_cos_ps(thti);
	                sthti  = _mm_sin_ps(thti);
	                cthts  = _mm_cos_ps(thts);
	                cterm  = _mm_mul_ps(cthti,_mm_mul_ps(cthts,cphis));
	                sthts  = _mm_sin_ps(thts);
	                sterm  = _mm_mul_ps(sthti,sthts);
	                sctrm  = _mm_sub_ps(sterm,cterm);
	                x0     = _mm_mul_ps(sphis,sphis);
	                x1     = _mm_mul_ps(cthts,cthts);
	                x2     = _mm_mul_ps(C30,_mm_mul_ps(cthti,cthti));
	                x3     = _mm_sub_ps(C30,_mm_mul_ps(x0,
	                                                    _mm_mul_ps(x1,x2)));
	                trm1   = _mm_fmadd_ps(t2,x3,sctrm);
	                trm2   = _mm_mul_ps(C120,_mm_mul_ps(t2,sterm));
	                x0     = _mm_fmsub_ps(C30,sterm,cterm);
	                trm3   = _mm_mul_ps(C140,_mm_mul_ps(t,x0));
	                trm4   = _mm_mul_ps(_mm_mul_ps(C350,t2),sctrm);
	                x1     = _mm_mul_ps(trm1,trm2);
	                x2     = _mm_add_ps(x1,_mm_add_ps(trm3,trm4));
	                rcs    = _mm_mul_ps(Bg,x2);
	                return (rcs);
	         }
	         
	         
	         
	              __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9136_xmm4r4_u( const float * __restrict  pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol ) {
	                                      
	                 __m128 A   = _mm_loadu_ps(&pA[0]);
	                 __m128 N   = _mm_loadu_ps(&pN[0]); 
	                 __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                 __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                 __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                 __m128 thti= _mm_loadu_ps(&pthti[0]);
	                 __m128 thts= _mm_loadu_ps(&pthts[0]);
	                 __m128 phis= _mm_loadu_ps(&pphis[0]);   
	                const __m128 C10 = _mm_set1_ps(1.0f);     
	                const __m128 C20 = _mm_set1_ps(2.0f);   
	                const __m128 C30 = _mm_set1_ps(3.0f);
	                const __m128 C120= _mm_set1_ps(12.0f);
	                const __m128 C140= _mm_set1_ps(14.0f);
	                const __m128 C350= _mm_set1_ps(35.0f);
	                 __m128 t,t2,x0,x1,x2,x3,cthti,cthts;
	                 __m128 sphis,cphis,sthts,cthts;
	                 __m128 trm1,trm2,trm3,trm4;
	                 __m128 rcs,Bg,cterm,sterm,sctrm;
	                x0     = _mm_div_ps(C20,_mm_add_ps(C10,epsr));   
	                t      = _mm_sub_ps(C10,x0);
	                t2     = _mm_sub_ps(C10,_mm_mul_ps(x0,x0));
	                Bg     = coef_Bg_f9137_xmm4r4(A,N,k0,epsr,epsi,thti,thts,phis,pol); 
	                sphis  = _mm_sin_ps(phis);
	                cphis  = _mm_cos_ps(phis);
	                cthti  = _mm_cos_ps(thti);
	                sthti  = _mm_sin_ps(thti);
	                cthts  = _mm_cos_ps(thts);
	                cterm  = _mm_mul_ps(cthti,_mm_mul_ps(cthts,cphis));
	                sthts  = _mm_sin_ps(thts);
	                sterm  = _mm_mul_ps(sthti,sthts);
	                sctrm  = _mm_sub_ps(sterm,cterm);
	                x0     = _mm_mul_ps(sphis,sphis);
	                x1     = _mm_mul_ps(cthts,cthts);
	                x2     = _mm_mul_ps(C30,_mm_mul_ps(cthti,cthti));
	                x3     = _mm_sub_ps(C30,_mm_mul_ps(x0,
	                                                    _mm_mul_ps(x1,x2)));
	                trm1   = _mm_fmadd_ps(t2,x3,sctrm);
	                trm2   = _mm_mul_ps(C120,_mm_mul_ps(t2,sterm));
	                x0     = _mm_fmsub_ps(C30,sterm,cterm);
	                trm3   = _mm_mul_ps(C140,_mm_mul_ps(t,x0));
	                trm4   = _mm_mul_ps(_mm_mul_ps(C350,t2),sctrm);
	                x1     = _mm_mul_ps(trm1,trm2);
	                x2     = _mm_add_ps(x1,_mm_add_ps(trm3,trm4));
	                rcs    = _mm_mul_ps(Bg,x2);
	                return (rcs);
	         }
	         
	         
	         
	        /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (vv) polarization.
	            Formula 9.1-74
	        */
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9174_xmm4r4(const __m128 k0,
	                                       const __m128 h,
	                                       const __m128 l,
	                                       const __m128 thti,
	                                       const __m128 epsr,
	                                       const __m128 epsi,
	                                       const __m128 mur,
	                                       const __m128 mui) {
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 epsrm1,epsim1;
	                   __m128 epsr2,epsi2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsi,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm_mul_ps(k02,_mm_mul_ps(l2,sthti));   
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C40,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm_fmadd_ps(epsrm1,sthti,t1r);
	                  epsim1= _mm_fmadd_ps(epsim1,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(epsr,cthti);
	                  x1    = _mm_mul_ps(epsi,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm_sub_ps(epsrm1,x2);
	                  epsim1 = _mm_sub_ps(epsim1,x3);
	                  cdiv_xmm4c4(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9174_xmm4r4_a(const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui) {
	                        
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 h   = _mm_load_ps(&ph[0]);
	                   __m128 l   = _mm_load_ps(&pl[0]);
	                   __m128 epsr= _mm_load_ps(&pepsr[0]);
	                   __m128 epsi= _mm_load_ps(&pepsi[0]);
	                   __m128 mur = _mm_load_ps(&pmur[0]);
	                   __m128 mui = _mm_load_ps(&pmui[0]);               
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 epsrm1,epsim1;
	                   __m128 epsr2,epsi2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsi,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm_mul_ps(k02,_mm_mul_ps(l2,sthti));   
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C40,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm_fmadd_ps(epsrm1,sthti,t1r);
	                  epsim1= _mm_fmadd_ps(epsim1,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(epsr,cthti);
	                  x1    = _mm_mul_ps(epsi,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm_sub_ps(epsrm1,x2);
	                  epsim1 = _mm_sub_ps(epsim1,x3);
	                  cdiv_xmm4c4(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9174_xmm4r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui) {
	                        
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 h   = _mm_loadu_ps(&ph[0]);
	                   __m128 l   = _mm_loadu_ps(&pl[0]);
	                   __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                   __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                   __m128 mur = _mm_loadu_ps(&pmur[0]);
	                   __m128 mui = _mm_loadu_ps(&pmui[0]);               
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 epsrm1,epsim1;
	                   __m128 epsr2,epsi2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsi,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm_mul_ps(k02,_mm_mul_ps(l2,sthti));   
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C40,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm_fmadd_ps(epsrm1,sthti,t1r);
	                  epsim1= _mm_fmadd_ps(epsim1,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(epsr,cthti);
	                  x1    = _mm_mul_ps(epsi,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm_sub_ps(epsrm1,x2);
	                  epsim1 = _mm_sub_ps(epsim1,x3);
	                  cdiv_xmm4c4(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	        /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hh) polarization.
	            Formula 9.1-75
	        */
	        
	        
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9175_xmm4r4(const __m128 k0,
	                                       const __m128 h,
	                                       const __m128 l,
	                                       const __m128 thti,
	                                       const __m128 epsr,
	                                       const __m128 epsi,
	                                       const __m128 mur,
	                                       const __m128 mui) {
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 murm1,muim1,murm1s,muim1s;
	                   __m128 mur2,mui2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsi,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsr,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm_mul_ps(k02,_mm_mul_ps(l2,sthti));   
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C40,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm_fmadd_ps(murm1s,sthti,t1r);
	                  muim1s= _mm_fmadd_ps(muim1s,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(mur,cthti);
	                  x1    = _mm_mul_ps(mui,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm_sub_ps(murm1,x2);
	                  muim1 = _mm_sub_ps(muim1,x3);
	                  cdiv_xmm4c4(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9175_xmm4r4_a(  const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui) {
	                        
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 h   = _mm_load_ps(&ph[0]);
	                   __m128 l   = _mm_load_ps(&pl[0]);
	                   __m128 epsr= _mm_load_ps(&pepsr[0]);
	                   __m128 epsi= _mm_load_ps(&pepsi[0]);
	                   __m128 mur = _mm_load_ps(&pmur[0]);
	                   __m128 mui = _mm_load_ps(&pmui[0]);  
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 murm1,muim1,murm1s,muim1s;
	                   __m128 mur2,mui2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsi,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsr,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm_mul_ps(k02,_mm_mul_ps(l2,sthti));   
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C40,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm_fmadd_ps(murm1s,sthti,t1r);
	                  muim1s= _mm_fmadd_ps(muim1s,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(mur,cthti);
	                  x1    = _mm_mul_ps(mui,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm_sub_ps(murm1,x2);
	                  muim1 = _mm_sub_ps(muim1,x3);
	                  cdiv_xmm4c4(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9175_xmm4r4_u(  const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui) {
	                        
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 h   = _mm_loadu_ps(&ph[0]);
	                   __m128 l   = _mm_loadu_ps(&pl[0]);
	                   __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                   __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                   __m128 mur = _mm_loadu_ps(&pmur[0]);
	                   __m128 mui = _mm_loadu_ps(&pmui[0]);  
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 murm1,muim1,murm1s,muim1s;
	                   __m128 mur2,mui2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsi,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsr,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm_mul_ps(k02,_mm_mul_ps(l2,sthti));   
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C40,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm_fmadd_ps(murm1s,sthti,t1r);
	                  muim1s= _mm_fmadd_ps(muim1s,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(mur,cthti);
	                  x1    = _mm_mul_ps(mui,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm_sub_ps(murm1,x2);
	                  muim1 = _mm_sub_ps(muim1,x3);
	                  cdiv_xmm4c4(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	         /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hv) polarization.
	            Formula 9.1-76
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hv_f9176_xmm4r4() { 
	           
	                return _mm_setzero_ps();
	         } 
	         
	         
	         
	         /*
	            Exponential surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (vv) polarization.
	            Formula 9.1-77
	        */
	        
	        
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9177_xmm4r4(const __m128 k0,
	                                       const __m128 h,
	                                       const __m128 l,
	                                       const __m128 thti,
	                                       const __m128 epsr,
	                                       const __m128 epsi,
	                                       const __m128 mur,
	                                       const __m128 mui) {
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  const __m128 C80  = _mm_set1_ps(8.0f);
	                  const __m128 C150 = _mm_set1_ps(1.5f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 epsrm1,epsim1;
	                   __m128 epsr2,epsi2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsi,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm_fmadd_pd(_mm_mul_pd(sthti,l2),
	                                          _mm_mul_pd(k0,C40),C10);
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  x3    = _mm_mul_ps(arg,_mm_sqrt_ps(arg));
	                  inve  = _mm_rcp14_ps(x3);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C80,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm_fmadd_ps(epsrm1,sthti,t1r);
	                  epsim1= _mm_fmadd_ps(epsim1,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(epsr,cthti);
	                  x1    = _mm_mul_ps(epsi,cthti);
	                  cmul_xmm4c4(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm_sub_ps(epsrm1,x2);
	                  epsim1 = _mm_sub_ps(epsim1,x3);
	                  cdiv_xmm4c4(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9177_xmm4r4_a(const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui) {
	                        
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 h   = _mm_load_ps(&ph[0]);
	                   __m128 l   = _mm_load_ps(&pl[0]);
	                   __m128 epsr= _mm_load_ps(&pepsr[0]);
	                   __m128 epsi= _mm_load_ps(&pepsi[0]);
	                   __m128 mur = _mm_load_ps(&pmur[0]);
	                   __m128 mui = _mm_load_ps(&pmui[0]); 
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  const __m128 C80  = _mm_set1_ps(8.0f);
	                  const __m128 C150 = _mm_set1_ps(1.5f);
	                  __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 epsrm1,epsim1;
	                   __m128 epsr2,epsi2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsi,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm_fmadd_pd(_mm_mul_pd(sthti,l2),
	                                          _mm_mul_pd(k0,C40),C10);
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  x3    = _mm_mul_ps(arg,_mm_sqrt_ps(arg));
	                  inve  = _mm_rcp14_ps(x3);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C80,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm_fmadd_ps(epsrm1,sthti,t1r);
	                  epsim1= _mm_fmadd_ps(epsim1,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(epsr,cthti);
	                  x1    = _mm_mul_ps(epsi,cthti);
	                  cmul_xmm4c4(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm_sub_ps(epsrm1,x2);
	                  epsim1 = _mm_sub_ps(epsim1,x3);
	                  cdiv_xmm4c4(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9177_xmm4r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui) {
	                        
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 h   = _mm_loadu_ps(&ph[0]);
	                   __m128 l   = _mm_loadu_ps(&pl[0]);
	                   __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                   __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                   __m128 mur = _mm_loadu_ps(&pmur[0]);
	                   __m128 mui = _mm_loadu_ps(&pmui[0]); 
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  const __m128 C80  = _mm_set1_ps(8.0f);
	                  const __m128 C150 = _mm_set1_ps(1.5f);
	                  __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 epsrm1,epsim1;
	                   __m128 epsr2,epsi2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsr,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsi,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm_fmadd_pd(_mm_mul_pd(sthti,l2),
	                                          _mm_mul_pd(k0,C40),C10);
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  x3    = _mm_mul_ps(arg,_mm_sqrt_ps(arg));
	                  inve  = _mm_rcp14_ps(x3);
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C80,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm_fmadd_ps(epsrm1,sthti,t1r);
	                  epsim1= _mm_fmadd_ps(epsim1,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(epsr,cthti);
	                  x1    = _mm_mul_ps(epsi,cthti);
	                  cmul_xmm4c4(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm_sub_ps(epsrm1,x2);
	                  epsim1 = _mm_sub_ps(epsim1,x3);
	                  cdiv_xmm4c4(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	         /*
	            Exponential surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hh) polarization.
	            Formula 9.1-78
	        */
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9178_xmm4r4(const __m128 k0,
	                                       const __m128 h,
	                                       const __m128 l,
	                                       const __m128 thti,
	                                       const __m128 epsr,
	                                       const __m128 epsi,
	                                       const __m128 mur,
	                                       const __m128 mui) {
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  const __m128 C80  = _mm_set1_ps(8.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 murm1,muim1,murm1s,muim1s;
	                   __m128 mur2,mui2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsi,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsr,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  x3    = _mm_mul_ps(C40,x0);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm_fmadd_ps(_mm_mul_ps(sthti,l2),x3,C10);
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_mul_ps(arg,_mm_sqrt_ps(arg));
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C80,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm_fmadd_ps(murm1s,sthti,t1r);
	                  muim1s= _mm_fmadd_ps(muim1s,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(mur,cthti);
	                  x1    = _mm_mul_ps(mui,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm_sub_ps(murm1,x2);
	                  muim1 = _mm_sub_ps(muim1,x3);
	                  cdiv_xmm4c4(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9178_xmm4r4_a(const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui) {
	                        
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 h   = _mm_load_ps(&ph[0]);
	                   __m128 l   = _mm_load_ps(&pl[0]);
	                   __m128 epsr= _mm_load_ps(&pepsr[0]);
	                   __m128 epsi= _mm_load_ps(&pepsi[0]);
	                   __m128 mur = _mm_load_ps(&pmur[0]);
	                   __m128 mui = _mm_load_ps(&pmui[0]); 
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  const __m128 C80  = _mm_set1_ps(8.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 murm1,muim1,murm1s,muim1s;
	                   __m128 mur2,mui2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsi,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsr,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  x3    = _mm_mul_ps(C40,x0);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm_fmadd_ps(_mm_mul_ps(sthti,l2),x3,C10);
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_mul_ps(arg,_mm_sqrt_ps(arg));
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C80,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm_fmadd_ps(murm1s,sthti,t1r);
	                  muim1s= _mm_fmadd_ps(muim1s,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(mur,cthti);
	                  x1    = _mm_mul_ps(mui,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm_sub_ps(murm1,x2);
	                  muim1 = _mm_sub_ps(muim1,x3);
	                  cdiv_xmm4c4(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9178_xmm4r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui) {
	                        
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 h   = _mm_loadu_ps(&ph[0]);
	                   __m128 l   = _mm_loadu_ps(&pl[0]);
	                   __m128 epsr= _mm_loadu_ps(&pepsr[0]);
	                   __m128 epsi= _mm_loadu_ps(&pepsi[0]);
	                   __m128 mur = _mm_loadu_ps(&pmur[0]);
	                   __m128 mui = _mm_loadu_ps(&pmui[0]); 
	                                       
	                  const __m128 C10  = _mm_set1_ps(1.0f);
	                  const __m128 C40  = _mm_set1_ps(4.0f);
	                  const __m128 C80  = _mm_set1_ps(8.0f);
	                   __m128 k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                   __m128 murm1,muim1,murm1s,muim1s;
	                   __m128 mur2,mui2,t0r,t0i,t1r,t1i;
	                   __m128 t2r,t2i,inve,t3r,t3i;
	                   __m128 rcs,cabs,arg,earg,frac;
	                  h2    = _mm_mul_ps(h,h);
	                  t0r   = _mm_sub_ps(epsi,C10);
	                  cthti = _mm_cos_ps(thti);
	                  t0i   = _mm_sub_ps(epsr,C10);
	                  murm1 = _mm_sub_ps(mur,C10);
	                  x0    = _mm_mul_ps(k0,k0);
	                  muim1 = _mm_mul_ps(mui,C10);
	                  x3    = _mm_mul_ps(C40,x0);
	                  l2    = _mm_mul_ps(l,l);
	                  sthti = _mm_sin_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  sthti = _mm_mul_ps(sthti,sthti);
	                  cmul_xmm4c4(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm_fmadd_ps(_mm_mul_ps(sthti,l2),x3,C10);
	                  x0    = _mm_mul_ps(cthti,cthti);
	                  earg  = _mm_mul_ps(arg,_mm_sqrt_ps(arg));
	                  x1    = _mm_mul_ps(x0,x0);
	                  x2    = _mm_mul_ps(C80,k04);
	                  frac  = _mm_mul_ps(_mm_mul_ps(x2,h2),
	                                        _mm_mul_ps(l2,x1));
	                  cmul_xmm4c4(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm_fmadd_ps(murm1s,sthti,t1r);
	                  muim1s= _mm_fmadd_ps(muim1s,sthti,t1i);
	                  t1r   = _mm_sub_ps(t1r,sthti);
	                  t1i   = _mm_sub_ps(t1i,sthti);
	                  csqrt_xmm4c4(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm_mul_ps(mur,cthti);
	                  x1    = _mm_mul_ps(mui,cthti);
	                  inve  = _mm_rcp14_ps(earg);
	                  cmul_xmm4c4(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm_add_ps(x0,t2r);
	                  x1    = _mm_add_ps(x1,t2i);
	                  cmul_xmm4c4(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_xmm4c4(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm_sub_ps(murm1,x2);
	                  muim1 = _mm_sub_ps(muim1,x3);
	                  cdiv_xmm4c4(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_xmm4c4(t0r,t0i);
	                  rcs  = _mm_mul_ps(frac,_mm_mul_ps(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	         /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hv) polarization.
	            Formula 9.1-79
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hv_f9179_xmm4r4() { 
	           
	                return _mm_setzero_ps();
	         } 
	         
	         
	         /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-80
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vv_f9180_xmm4r4(const __m128 thti,
	                                         const __m128 thts,
	                                         const __m128 phis) {
	                                         
	                   __m128 sthti,sthts,cthti,cthts,cphis;
	                   __m128 num,den,avv;  
	                  sthti = _mm_sin_ps(thti);
	                  cthti = _mm_cos_ps(thti);
	                  cphis = _mm_cos_ps(phis);
	                  cthts = _mm_cos_ps(thts);
	                  sthts = _mm_sin_ps(thts);
	                  num   = _mm_fmsub_ps(sthti,sthts,cphis);
	                  den   = _mm_mul_ps(cthti,cthts);
	                  avv   = _mm_div_ps(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vv_f9180_xmm4r4_a(  const float * __restrict pthti,
	                                         const float * __restrict pthts,
	                                         const float * __restrict pphis) {
	                     
	                   __m128 thti = _mm_load_ps(&pthti[0]);
	                   __m128 thts = _mm_load_ps(&pthts[0]);
	                   __m128 phis = _mm_load_ps(&phis[0]);                    
	                   __m128 sthti,sthts,cthti,cthts,cphis;
	                   __m128 num,den,avv;  
	                  sthti = _mm_sin_ps(thti);
	                  cthti = _mm_cos_ps(thti);
	                  cphis = _mm_cos_ps(phis);
	                  cthts = _mm_cos_ps(thts);
	                  sthts = _mm_sin_ps(thts);
	                  num   = _mm_fmsub_ps(sthti,sthts,cphis);
	                  den   = _mm_mul_ps(cthti,cthts);
	                  avv   = _mm_div_ps(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vv_f9180_xmm4r4_u(  const float * __restrict  pthti,
	                                         const float * __restrict  pthts,
	                                         const float * __restrict  pphis) {
	                     
	                   __m128 thti = _mm_loadu_ps(&pthti[0]);
	                   __m128 thts = _mm_loadu_ps(&pthts[0]);
	                   __m128 phis = _mm_loadu_ps(&phis[0]);                    
	                   __m128 sthti,sthts,cthti,cthts,cphis;
	                   __m128 num,den,avv;  
	                  sthti = _mm_sin_ps(thti);
	                  cthti = _mm_cos_ps(thti);
	                  cphis = _mm_cos_ps(phis);
	                  cthts = _mm_cos_ps(thts);
	                  sthts = _mm_sin_ps(thts);
	                  num   = _mm_fmsub_ps(sthti,sthts,cphis);
	                  den   = _mm_mul_ps(cthti,cthts);
	                  avv   = _mm_div_ps(num,den);
	                  return (avv);        
	          }
	          
	          
	          /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-81
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hv_f9181_xmm4r4(const __m128 phis,
	                                     const __m128 thti) {
	                                     
	                   __m128 sphis,cthti;
	                   __m128 ahv;
	                  sphis = _mm_sin_ps(phis);
	                  cthti = _mm_cos_ps(thti);
	                  ahv   = _mm_div_ps(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hv_f9181_xmm4r4_a(const float * __restrict pphis,
	                                     const float * __restrict pthti) {
	                                     
	                   __m128 thti = _mm_load_ps(&pthti[0]);
	                   __m128 phis = _mm_load_ps(&phis[0]); 
	                   __m128 sphis,cthti;
	                   __m128 ahv;
	                  sphis = _mm_sin_ps(phis);
	                  cthti = _mm_cos_ps(thti);
	                  ahv   = _mm_div_ps(sphis,cthti);
	                  return (ahv);                           
	         }
	          
	       
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hv_f9181_xmm4r4_u(const float * __restrict pphis,
	                                     const float * __restrict  pthti) {
	                                     
	                   __m128 thti = _mm_loadu_ps(&pthti[0]);
	                   __m128 phis = _mm_loadu_ps(&phis[0]); 
	                   __m128 sphis,cthti;
	                   __m128 ahv;
	                  sphis = _mm_sin_ps(phis);
	                  cthti = _mm_cos_ps(thti);
	                  ahv   = _mm_div_ps(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	            /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-82
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hh_f9182_xmm4r4(const __m128 phis) {
	                 
	                   __m128 cphis,ahh;
	                  cphis = _mm_cos_ps(phis);
	                  ahh   = negate_xmm4r4(cphis);
	                  return (ahh);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hh_f9182_xmm4r4_a(const float * __restrict pphis) {
	                 
	                   __m128 phis = _mm_load_ps(&pphis[0]);
	                   __m128 cphis,ahh;
	                  cphis = _mm_cos_ps(phis);
	                  ahh   = negate_xmm4r4(cphis);
	                  return (ahh);
	          }
	        
	       
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hh_f9182_xmm4r4_u(const float * __restrict  pphis) {
	                 
	                   __m128 phis = _mm_loadu_ps(&pphis[0]);
	                   __m128 cphis,ahh;
	                  cphis = _mm_cos_ps(phis);
	                  ahh   = negate_xmm4r4(cphis);
	                  return (ahh);
	          }
	          
	          
	             /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-83
	         */
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vh_f9183_xmm4r4(const __m128 phis,
	                                     const __m128 thts) {
	                                     
	                   __m128 sphis,cthts;
	                   __m128 ahv;
	                  sphis = _mm_sin_ps(phis);
	                  cthts = _mm_cos_ps(thti);
	                  ahv   = gms::math::negate_xmm4r4(_mm_div_ps(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vh_f9183_xmm4r4_a(const float * __restrict pphis,
	                                       const float * __restrict pthts) {
	                         
	                   __m128 phis = _mm_load_ps(&pphis[0]);
	                   __m128 thts = _mm_load_ps(&pphis[0]);            
	                   __m128 sphis,cthts;
	                   __m128 ahv;
	                  sphis = _mm_sin_ps(phis);
	                  cthts = _mm_cos_ps(thti);
	                  ahv   = gms::math::negate_xmm4r4(_mm_div_ps(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vh_f9183_xmm4r4_u(const float * __restrict  pphis,
	                                       const float * __restrict pthts) {
	                         
	                   __m128 phis = _mm_loadu_ps(&pphis[0]);
	                   __m128 thts = _mm_loadu_ps(&pphis[0]);            
	                   __m128 sphis,cthts;
	                   __m128 ahv;
	                  sphis = _mm_sin_ps(phis);
	                  cthts = _mm_cos_ps(thti);
	                  ahv   = gms::math::negate_xmm4r4(_mm_div_ps(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-84 
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vv_f9184_xmm4r4(const __m128 thti) {
	           
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 sthti,cthti,num,den;
	                   __m128 avv;
	                  sthti = _mm_sin_ps(thti);
	                  cthti = _mm_cos_ps(thti);
	                  num   = _mm_fmadd_ps(sthti,sthti,C10);
	                  den   = _mm_mul_ps(cthti,cthti);
	                  avv   = _mm_div_ps(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vv_f9184_xmm4r4_a(const float * __restrict pthti) {
	           
	                   __m128 thti = _mm_load_ps(&pthti[0]);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 sthti,cthti,num,den;
	                   __m128 avv;
	                  sthti = _mm_sin_ps(thti);
	                  cthti = _mm_cos_ps(thti);
	                  num   = _mm_fmadd_ps(sthti,sthti,C10);
	                  den   = _mm_mul_ps(cthti,cthti);
	                  avv   = _mm_div_ps(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vv_f9184_xmm4r4_u(const float * __restrict  pthti) {
	           
	                   __m128 thti = _mm_loadu_ps(&pthti[0]);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 sthti,cthti,num,den;
	                   __m128 avv;
	                  sthti = _mm_sin_ps(thti);
	                  cthti = _mm_cos_ps(thti);
	                  num   = _mm_fmadd_ps(sthti,sthti,C10);
	                  den   = _mm_mul_ps(cthti,cthti);
	                  avv   = _mm_div_ps(num,den);
	                  return (avv);
	           }
	           
	           
	           
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-85
	         */                  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_hh_f9185_xmm4r4() {
	           
	                  return _mm_set1_ps(1.0f);
	           }  
	           
	           
	            /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-86
	         */   
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 a_vh_f9186_xmm4r4() {
	           
	                  return _mm_setzero_ps();
	           }  
	           
	           
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Average backscattering RCS per unit area.
	              Gaussian surface height correlation coefficient.
	              Formula: 9.1-87
	         */
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9187_xmm4r4(const __m128 k0,
	                                       const __m128 h,
	                                       const __m128 l,
	                                       const __m128 thti) {
	                                       
	                  const __m128 C40 = _mm_set1_ps(4.0f);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 k04,x0,x1,l2,h2,sthti;
	                   __m128 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm_mul_ps(k0,k0);
	                  l2    = _mm_mul_ps(l,l);
	                  k04   = _mm_mul_ps(x0,x0);
	                  h2    = _mm_mul_ps(h,h);
	                  x1    = _mm_sin_ps(thti);
	                  strm  = _mm_fmadd_pd(x1,x1,C10);
	                  x0    = _mm_mul_ps(C40,k04);
	                  sthti = _mm_mul_ps(x1,x1);
	                  arg   = _mm_mul_ps(_mm_mul_ps(k0,k0),
	                                        _mm_mul_ps(sthti,l2));
	                  trm   = _mm_mul_ps(_mm_mul_ps(l2,h2,x0));
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(strm,strm);
	                  inve  = _mm_rcp14_ps(earg);
	                  rcs   = _mm_mul_ps(trm,_mm_mul_ps(x1,inve));
	                  return (rcs);
	         }
	           
	           
	           
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9187_xmm4r4_a(const float * __restrict pk0,
	                                       const float * __restrict ph,
	                                       const float * __restrict pl,
	                                       const float * __restrict pthti) {
	                             
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 h   = _mm_load_ps(&ph[0]);
	                   __m128 l   = _mm_load_ps(&pl[0]);   
	                   __m128 thti= _mm_load_ps(&pthti[0]);       
	                  const __m128 C40 = _mm_set1_ps(4.0f);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 k04,x0,x1,l2,h2,sthti;
	                   __m128 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm_mul_ps(k0,k0);
	                  l2    = _mm_mul_ps(l,l);
	                  k04   = _mm_mul_ps(x0,x0);
	                  h2    = _mm_mul_ps(h,h);
	                  x1    = _mm_sin_ps(thti);
	                  strm  = _mm_fmadd_pd(x1,x1,C10);
	                  x0    = _mm_mul_ps(C40,k04);
	                  sthti = _mm_mul_ps(x1,x1);
	                  arg   = _mm_mul_ps(_mm_mul_ps(k0,k0),
	                                        _mm_mul_ps(sthti,l2));
	                  trm   = _mm_mul_ps(_mm_mul_ps(l2,h2,x0));
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(strm,strm);
	                  inve  = _mm_rcp14_ps(earg);
	                  rcs   = _mm_mul_ps(trm,_mm_mul_ps(x1,inve));
	                  return (rcs);
	         }
	           
	            
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vv_f9187_xmm4r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti) {
	                             
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 h   = _mm_loadu_ps(&ph[0]);
	                   __m128 l   = _mm_loadu_ps(&pl[0]);   
	                   __m128 thti= _mm_loadu_ps(&pthti[0]);       
	                  const __m128 C40 = _mm_set1_ps(4.0f);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 k04,x0,x1,l2,h2,sthti;
	                   __m128 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm_mul_ps(k0,k0);
	                  l2    = _mm_mul_ps(l,l);
	                  k04   = _mm_mul_ps(x0,x0);
	                  h2    = _mm_mul_ps(h,h);
	                  x1    = _mm_sin_ps(thti);
	                  strm  = _mm_fmadd_pd(x1,x1,C10);
	                  x0    = _mm_mul_ps(C40,k04);
	                  sthti = _mm_mul_ps(x1,x1);
	                  arg   = _mm_mul_ps(_mm_mul_ps(k0,k0),
	                                        _mm_mul_ps(sthti,l2));
	                  trm   = _mm_mul_ps(_mm_mul_ps(l2,h2,x0));
	                  earg  = _mm_exp_ps(arg);
	                  x1    = _mm_mul_ps(strm,strm);
	                  inve  = _mm_rcp14_ps(earg);
	                  rcs   = _mm_mul_ps(trm,_mm_mul_ps(x1,inve));
	                  return (rcs);
	         }
	         
	         
	          /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Average backscattering RCS per unit area.
	              Gaussian surface height correlation coefficient.
	              Formula: 9.1-88
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9188_xmm4r4(const __m128 k0,
	                                       const __m128 h,
	                                       const __m128 l,
	                                       const __m128 thti) {
	                                       
	                  const __m128 C40 = _mm_set1_ps(4.0f);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 k04,x0,x1,l2,h2,sthti,cthti;
	                   __m128 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm_mul_ps(k0,k0);
	                  l2    = _mm_mul_ps(l,l);
	                  strm  = _mm_cos_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  h2    = _mm_mul_ps(h,h);
	                  x1    = _mm_sin_ps(thti);
	                  x0    = _mm_mul_ps(C40,k04);
	                  sthti = _mm_mul_ps(x1,x1);
	                  arg   = _mm_mul_ps(_mm_mul_ps(k0,k0),
	                                        _mm_mul_ps(sthti,l2));
	                  trm   = _mm_mul_ps(_mm_mul_ps(l2,h2,x0));
	                  x1    = _mm_mul_ps(strm,strm);
	                  earg  = _mm_exp_ps(arg);
	                  cthti = _mm_mul_ps(x1,x1);
	                  inve  = _mm_rcp14_ps(earg);
	                  rcs   = _mm_mul_ps(trm,_mm_mul_ps(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9188_xmm4r4_a(const float * __restrict pk0,
	                                       const float * __restrict ph,
	                                       const float * __restrict pl,
	                                       const float * __restrict pthti) {
	                             
	                   __m128 k0  = _mm_load_ps(&pk0[0]);
	                   __m128 h   = _mm_load_ps(&ph[0]);
	                   __m128 l   = _mm_load_ps(&pl[0]);   
	                   __m128 thti= _mm_load_ps(&pthti[0]); 
	                                       
	                  const __m128 C40 = _mm_set1_ps(4.0f);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 k04,x0,x1,l2,h2,sthti,cthti;
	                   __m128 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm_mul_ps(k0,k0);
	                  l2    = _mm_mul_ps(l,l);
	                  strm  = _mm_cos_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  h2    = _mm_mul_ps(h,h);
	                  x1    = _mm_sin_ps(thti);
	                  x0    = _mm_mul_ps(C40,k04);
	                  sthti = _mm_mul_ps(x1,x1);
	                  arg   = _mm_mul_ps(_mm_mul_ps(k0,k0),
	                                        _mm_mul_ps(sthti,l2));
	                  trm   = _mm_mul_ps(_mm_mul_ps(l2,h2,x0));
	                  x1    = _mm_mul_ps(strm,strm);
	                  earg  = _mm_exp_ps(arg);
	                  cthti = _mm_mul_ps(x1,x1);
	                  inve  = _mm_rcp14_ps(earg);
	                  rcs   = _mm_mul_ps(trm,_mm_mul_ps(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_hh_f9188_xmm4r4_u(const float * __restrict  pk0,
	                                       const float * __restrict  ph,
	                                       const float * __restrict  pl,
	                                       const float * __restrict  pthti) {
	                             
	                   __m128 k0  = _mm_loadu_ps(&pk0[0]);
	                   __m128 h   = _mm_loadu_ps(&ph[0]);
	                   __m128 l   = _mm_loadu_ps(&pl[0]);   
	                   __m128 thti= _mm_loadu_ps(&pthti[0]); 
	                                       
	                  const __m128 C40 = _mm_set1_ps(4.0f);
	                  const __m128 C10 = _mm_set1_ps(1.0f);
	                   __m128 k04,x0,x1,l2,h2,sthti,cthti;
	                   __m128 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm_mul_ps(k0,k0);
	                  l2    = _mm_mul_ps(l,l);
	                  strm  = _mm_cos_ps(thti);
	                  k04   = _mm_mul_ps(x0,x0);
	                  h2    = _mm_mul_ps(h,h);
	                  x1    = _mm_sin_ps(thti);
	                  x0    = _mm_mul_ps(C40,k04);
	                  sthti = _mm_mul_ps(x1,x1);
	                  arg   = _mm_mul_ps(_mm_mul_ps(k0,k0),
	                                        _mm_mul_ps(sthti,l2));
	                  trm   = _mm_mul_ps(_mm_mul_ps(l2,h2,x0));
	                  x1    = _mm_mul_ps(strm,strm);
	                  earg  = _mm_exp_ps(arg);
	                  cthti = _mm_mul_ps(x1,x1);
	                  inve  = _mm_rcp14_ps(earg);
	                  rcs   = _mm_mul_ps(trm,_mm_mul_ps(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	         /*
	              
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Average backscattering RCS per unit area.
	              Gaussian surface height correlation coefficient.
	              Formula: 9.1-89
	           
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m128 rcs_vhhv_f9189_xmm4r4() {
	           
	                  return (_mm_setzero_ps()); 
	          }
	          
	          
	          
	           
	      
	                                    
	        
                 
                 
                 
               
               
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_XMM4R4_HPP__*/
