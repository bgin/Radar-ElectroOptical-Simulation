

#ifndef __GMS_RCS_COMPLEX_OBJECTS_YMM4R8_HPP__
#define __GMS_RCS_COMPLEX_OBJECTS_YMM4R8_HPP__ 260920240740

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

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM4R8_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM4R8_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM4R8_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM4R8_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_YMM4R8_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_YMM4R8_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_YMM4R8_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM4R8_CREATION_DATE = "26-09-2024 07:40 PM +00200 (FRI 26 SEP 2024 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM4R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM4R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM4R8_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_ymm4r8.hpp"
#include "GMS_simd_utils.hpp"




namespace  gms {


         namespace radiolocation {
         
         
               /*
                   Work (input) arrays for kernel rcs_f8162_ymm4r8_2t_u and
                   rcs_f8162_ymm4r8_2t_a.
               */
             __ATTR_ALIGN__(32)  struct RCS_F8162_DATA {
               
                       double * __restrict  Ya1; 
                       double * __restrict  Ya2; 
                       double * __restrict  Ya3; 
                       double * __restrict  Ea;  
                       double * __restrict  WRKa; 
                       double * __restrict  Yb1;
                       double * __restrict  Yb2; 
                       double * __restrict  Yb3; 
                       double * __restrict  Eb; 
                       double * __restrict  WRKb;  
               };
         
              
              /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Formula: 8.1-21
              */     
              
              
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8121_ymm4r8(const __m256d gam,
                                             const __m256d phi,
                                             const __m256d k0,
                                             __m256d * __restrict D1r,
                                             __m256d * __restrict D1i,
                                            __m256d * __restrict D2r,
                                            __m256d * __restrict D2i) {
                                            
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);
                        const __m256d C10                              = _mm256_set1_pd(1.0f);
                        register __m256d invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm256_add_pd(phi,phi);
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328); 
                        ear  = _mm256_setzero_pd();
                        invn = _mm256_rcp14_pd(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(k0,
                                                C6283185307179586476925286766559));
                        sinp = _mm256_sin_pd(x1);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        cos1 = _mm256_cos_pd(x1);
                        x0   = _mm256_mul_pd(invn,phi2);
                        cos2 = _mm256_cos_pd(x0);
                        cer  = _mm256_mul_pd(cer,sinp);
                        trm1 = _mm256_rcp14_pd(_mm256_sub_pd(cos1),C10);
                        cei  = _mm256_mul_pd(cei,sinp);
                        trm2 = _mm256_rcp14_pd(_mm256_sub_pd(cos1,cos2));
                        sqr  = _mm256_mul_pd(invn,sqr);
                        ear  = _mm256_mul_pd(cer,sqr);
                        eai  = _mm256_mul_pd(cei,sqr);
                        x0   = _mm256_sub_pd(trm1,trm2);
                        *D1r = _mm256_mul_pd(ear,x0);
                        *D1i = _mm256_mul_pd(eai,x0);
                        x1   = _mm256_add_pd(trm1,trm2);
                        *D2r = _mm256_mul_pd(ear,x1);
                        *D2i = _mm256_mul_pd(eai,x1);
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8121_ymm4r8_a(const double * __restrict pgam,
                                               const double * __restrict pphi,
                                               const double * __restrict pk0,
                                               double * __restrict D1r,
                                               double * __restrict D1i,
                                               double * __restrict D2r,
                                               double * __restrict D2i) {
                                   
                        register __m256d gam = _mm256_load_pd(&pgam[0]);
                        register __m256d phi = _mm256_load_pd(&pphi[0]);  
                        register __m256d k0  = _mm256_load_pd(&pk0[0]);       
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);
                        const __m256d C10                              = _mm256_set1_pd(1.0f);
                        register __m256d invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm256_add_pd(phi,phi);
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328); 
                        ear  = _mm256_setzero_pd();
                        invn = _mm256_rcp14_pd(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(k0,
                                                C6283185307179586476925286766559));
                        sinp = _mm256_sin_pd(x1);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        cos1 = _mm256_cos_pd(x1);
                        x0   = _mm256_mul_pd(invn,phi2);
                        cos2 = _mm256_cos_pd(x0);
                        cer  = _mm256_mul_pd(cer,sinp);
                        trm1 = _mm256_rcp14_pd(_mm256_sub_pd(cos1),C10);
                        cei  = _mm256_mul_pd(cei,sinp);
                        trm2 = _mm256_rcp14_pd(_mm256_sub_pd(cos1,cos2));
                        sqr  = _mm256_mul_pd(invn,sqr);
                        ear  = _mm256_mul_pd(cer,sqr);
                        eai  = _mm256_mul_pd(cei,sqr);
                        x0   = _mm256_sub_pd(trm1,trm2);
                        _mm256_store_pd(&D1r[0] ,_mm256_mul_pd(ear,x0));
                        _mm256_store_pd(&D1i[0] ,_mm256_mul_pd(eai,x0));
                        x1   = _mm256_add_pd(trm1,trm2);
                        _mm256_store_pd(&D2r[0] ,_mm256_mul_pd(ear,x1));
                        _mm256_store_pd(&D2i[0] ,_mm256_mul_pd(eai,x1));
                }
                
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8121_ymm4r8_u(const double * __restrict  pgam,
                                               const double * __restrict  pphi,
                                               const double * __restrict  pk0,
                                               double * __restrict  D1r,
                                               double * __restrict  D1i,
                                               double * __restrict  D2r,
                                               double * __restrict  D2i) {
                                   
                        register __m256d gam = _mm256_loadu_pd(&pgam[0]);
                        register __m256d phi = _mm256_loadu_pd(&pphi[0]);  
                        register __m256d k0  = _mm256_loadu_pd(&pk0[0]);       
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);
                        const __m256d C10                              = _mm256_set1_pd(1.0f);
                        register __m256d invn,x0,x1,ear,eai,cer,cei,phi2,sqr,sinp,cos1,cos2,trm1,trm2;
                        phi2 = _mm256_add_pd(phi,phi);
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328); 
                        ear  = _mm256_setzero_pd();
                        invn = _mm256_rcp14_pd(x0);
                        eai  = C078539816339744830961566084582;
                        x1   = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(k0,
                                                C6283185307179586476925286766559));
                        sinp = _mm256_sin_pd(x1);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        cos1 = _mm256_cos_pd(x1);
                        x0   = _mm256_mul_pd(invn,phi2);
                        cos2 = _mm256_cos_pd(x0);
                        cer  = _mm256_mul_pd(cer,sinp);
                        trm1 = _mm256_rcp14_pd(_mm256_sub_pd(cos1),C10);
                        cei  = _mm256_mul_pd(cei,sinp);
                        trm2 = _mm256_rcp14_pd(_mm256_sub_pd(cos1,cos2));
                        sqr  = _mm256_mul_pd(invn,sqr);
                        ear  = _mm256_mul_pd(cer,sqr);
                        eai  = _mm256_mul_pd(cei,sqr);
                        x0   = _mm256_sub_pd(trm1,trm2);
                        _mm256_storeu_pd(&D1r[0] ,_mm256_mul_pd(ear,x0));
                        _mm256_storeu_pd(&D1i[0] ,_mm256_mul_pd(eai,x0));
                        x1   = _mm256_add_pd(trm1,trm2);
                        _mm256_storeu_pd(&D2r[0] ,_mm256_mul_pd(ear,x1));
                        _mm256_storeu_pd(&D2i[0] ,_mm256_mul_pd(eai,x1));
                }
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter singly diffracted far-zone fields (E,H).
                    Formula: 8.1-19, 8.1-20
                
                */
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f811920_ymm4r8(    const __m256d betai,
                                                 const __m256d betas,
                                                 const __m256d gam,
                                                 const __m256d phi,
                                                 const __m256d k0,
                                                 const __m256d r,
                                                 const __m256d rho,
                                                 const __m256d psi,
                                                 __m256d * __restrict Esr,
                                                 __m256d * __restrict Esi,
                                                 __m256d * __restrict Hsr,
                                                 __m256d * __restrict Hsi) {
                                                 
                       register __m256d ear,eai,cer,cei;
                       register __m256d D1r,D1i,D2r,D2i,x0,x1;
                       register __m256d rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm256_fmadd_pd(k0,r,psi);
                       coef_D12_f8121_ymm4r8(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm256_rcp14_pd(r);
                       ear   = _mm256_setzero_pd();
                       cosbi = _mm256_cos_pd(betai);
                       eai   = k0rp;
                       cosbs = _mm256_cos_pd(betas);
                       cexp_ymm4c8(ear,eai,&cer,&cei);
                       cer   = _mm256_mul_pd(cer,invr);    
                       rhos  = _mm256_div_pd(rho,_mm256_add_pd(cosbi,cosbs));
                       cei   = _mm256_mul_pd(cei,invr);
                       sqrho = _mm256_sqrt_pd(rhos);
                       x0    = _mm256_mul_pd(sqrho,cer);
                       x1    = _mm256_mul_pd(sqrho,cei);
                       *Esr  = _mm256_mul_pd(D1r,x0);
                       *Hsr  = _mm256_mul_pd(D2r,x0);
                       *Esi  = _mm256_mul_pd(D1i,x1);
                       *Hsi  = _mm256_mul_pd(D2i,x1);                               
            }
            
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f811920_ymm4r8_a(  const double * __restrict pbetai,
                                                 const double * __restrict pbetas,
                                                 const double * __restrict pgam,
                                                 const double * __restrict pphi,
                                                 const double * __restrict pk0,
                                                 const double * __restrict pr,
                                                 const double * __restrict prho,
                                                 const double * __restrict ppsi,
                                                 double * __restrict Esr,
                                                 double * __restrict Esi,
                                                 double * __restrict Hsr,
                                                 double * __restrict Hsi) {
                              
                       register __m256d betai = _mm256_load_pd(&pbetai[0]);
                       register __m256d betas = _mm256_load_pd(&pbetas[0]); 
                       register __m256d gam   = _mm256_load_pd(&pgam[0]);   
                       register __m256d k0    = _mm256_load_pd(&pk0[0]); 
                       register __m256d r     = _mm256_load_pd(&pr[0]);
                       register __m256d rho   = _mm256_load_pd(&prho[0]); 
                       register __m256d psi   = _mm256_load_pd(&ppsi[0]);             
                       register __m256d ear,eai,cer,cei;
                       register __m256d D1r,D1i,D2r,D2i,x0,x1;
                       register __m256d rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm256_fmadd_pd(k0,r,psi);
                       coef_D12_f8121_ymm4r8(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm256_rcp14_pd(r);
                       ear   = _mm256_setzero_pd();
                       cosbi = _mm256_cos_pd(betai);
                       eai   = k0rp;
                       cosbs = _mm256_cos_pd(betas);
                       cexp_ymm4c8(ear,eai,&cer,&cei);
                       cer   = _mm256_mul_pd(cer,invr);    
                       rhos  = _mm256_div_pd(rho,_mm256_add_pd(cosbi,cosbs));
                       cei   = _mm256_mul_pd(cei,invr);
                       sqrho = _mm256_sqrt_pd(rhos);
                       x0    = _mm256_mul_pd(sqrho,cer);
                       x1    = _mm256_mul_pd(sqrho,cei);
                       _mm256_store_pd(&Esr[0] ,_mm256_mul_pd(D1r,x0));
                       _mm256_store_pd(&Hsr[0] ,_mm256_mul_pd(D2r,x0));
                       _mm256_store_pd(&Esi[0] ,_mm256_mul_pd(D1i,x1));
                       _mm256_store_pd(&Hsi[0] ,_mm256_mul_pd(D2i,x1));                               
            }
            
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f811920_ymm4r8_u(  const double * __restrict  pbetai,
                                                 const double * __restrict  pbetas,
                                                 const double * __restrict  pgam,
                                                 const double * __restrict  pphi,
                                                 const double * __restrict  pk0,
                                                 const double * __restrict  pr,
                                                 const double * __restrict  prho,
                                                 const double * __restrict ppsi,
                                                 double * __restrict  Esr,
                                                 double * __restrict  Esi,
                                                 double * __restrict  Hsr,
                                                 double * __restrict  Hsi) {
                              
                       register __m256d betai = _mm256_loadu_pd(&pbetai[0]);
                       register __m256d betas = _mm256_loadu_pd(&pbetas[0]); 
                       register __m256d gam   = _mm256_loadu_pd(&pgam[0]);   
                       register __m256d k0    = _mm256_loadu_pd(&pk0[0]); 
                       register __m256d r     = _mm256_loadu_pd(&pr[0]);
                       register __m256d rho   = _mm256_loadu_pd(&prho[0]); 
                       register __m256d psi   = _mm256_loadu_pd(&ppsi[0]);             
                       register __m256d ear,eai,cer,cei;
                       register __m256d D1r,D1i,D2r,D2i,x0,x1;
                       register __m256d rhos,cosb1,cosbs,sqrho,k0rp,invr;
                       k0rp  = _mm256_fmadd_pd(k0,r,psi);
                       coef_D12_f8121_ymm4r8(gam,phi,k0,&D1r,&D1i,&D2r,&D2i);
                       invr  = _mm256_rcp14_pd(r);
                       ear   = _mm256_setzero_pd();
                       cosbi = _mm256_cos_pd(betai);
                       eai   = k0rp;
                       cosbs = _mm256_cos_pd(betas);
                       cexp_ymm4c8(ear,eai,&cer,&cei);
                       cer   = _mm256_mul_pd(cer,invr);    
                       rhos  = _mm256_div_pd(rho,_mm256_add_pd(cosbi,cosbs));
                       cei   = _mm256_mul_pd(cei,invr);
                       sqrho = _mm256_sqrt_pd(rhos);
                       x0    = _mm256_mul_pd(sqrho,cer);
                       x1    = _mm256_mul_pd(sqrho,cei);
                       _mm256_storeu_pd(&Esr[0] ,_mm256_mul_pd(D1r,x0));
                       _mm256_storeu_pd(&Hsr[0] ,_mm256_mul_pd(D2r,x0));
                       _mm256_storeu_pd(&Esi[0] ,_mm256_mul_pd(D1i,x1));
                       _mm256_storeu_pd(&Hsi[0] ,_mm256_mul_pd(D2i,x1));                               
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
                   void coef_D12_f8124_ymm4r8(const __m256d k0,
                                               const __m256d gam,
                                               __m256d * __restrict D1r,
                                               __m256d * __restrict D1i,
                                               __m256d * __restrict D2r,
                                               __m256d * __restrict D2i) {
                                               
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);
                        const __m256d C10                              = _mm256_set1_pd(1.0f);  
                        register __m256d ear,eai,cer,cei,t0r,t0i;
                        register __m256d x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm256_setzero_pd();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328);  
                        x1   = _mm256_mul_pd(k0,C6283185307179586476925286766559);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        invn = _mm256_rcp14_pd(x0);
                        cer  = gms::math::negate_ymm4r8(cer);
                        invn2= _mm256_add_pd(invn,invn);      
                        cei  = gms::math::negate_ymm4r8(cei);
                        arg1 = _mm256_mul_pd(C314159265358979323846264338328,invn2)
                        sin1 = _mm256_sin_pd(arg1);
                        arg2 = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        t0r  = _mm256_div_pd(cer,x1);
                        cos1 = _mm256_cos_pd(arg1);
                        x0   = _mm256_div_pd(cos1,sin1);
                        t0i  = _mm256_div_pd(cei,x1);   
                        sin2 = _mm256_sin_pd(arg2);
                        cos2 = _mm256_cos_pd(arg2);   
                        x1   = _mm256_div_pd(cos2,sin2);
                        ear  = _mm256_fmsub_pd(invn,x0,_mm256_mul_pd(invn2,x1));
                        eai  = _mm256_fmadd_pd(invn,x0,_mm256_mul_pd(invn2,x1));
                        cer  = _mm256_mul_pd(t0r,ear);
                        cei  = _mm256_mul_pd(t0i,eai);
                        *D1r = cer;
                        *D2r = cer;
                        *D1i = cei;
                        *D2i = cei;       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8124_ymm4r8_a(const  double * __restrict pk0,
                                                 const  double * __restrict pgam,
                                                 double * __restrict D1r,
                                                 double * __restrict D1i,
                                                 double * __restrict D2r,
                                                 double * __restrict D2i) {
                                    
                        register __m256d k0 = _mm256_load_pd(&pk0[0]);
                        register __m256d gam= _mm256_load_pd(&pgam[0]);           
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);
                        const __m256d C10                              = _mm256_set1_pd(1.0f);  
                        register __m256d ear,eai,cer,cei,t0r,t0i;
                        register __m256d x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm256_setzero_pd();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328);  
                        x1   = _mm256_mul_pd(k0,C6283185307179586476925286766559);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        invn = _mm256_rcp14_pd(x0);
                        cer  = gms::math::negate_ymm4r8(cer);
                        invn2= _mm256_add_pd(invn,invn);      
                        cei  = gms::math::negate_ymm4r8(cei);
                        arg1 = _mm256_mul_pd(C314159265358979323846264338328,invn2)
                        sin1 = _mm256_sin_pd(arg1);
                        arg2 = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        t0r  = _mm256_div_pd(cer,x1);
                        cos1 = _mm256_cos_pd(arg1);
                        x0   = _mm256_div_pd(cos1,sin1);
                        t0i  = _mm256_div_pd(cei,x1);   
                        sin2 = _mm256_sin_pd(arg2);
                        cos2 = _mm256_cos_pd(arg2);   
                        x1   = _mm256_div_pd(cos2,sin2);
                        ear  = _mm256_fmsub_pd(invn,x0,_mm256_mul_pd(invn2,x1));
                        eai  = _mm256_fmadd_pd(invn,x0,_mm256_mul_pd(invn2,x1));
                        cer  = _mm256_mul_pd(t0r,ear);
                        cei  = _mm256_mul_pd(t0i,eai);
                        _mm256_store_pd(&D1r[0] ,cer);
                        _mm256_store_pd(&D2r[0] ,cer);
                        _mm256_store_pd(&D1i[0] ,cei);
                        _mm256_store_pd(&D2i[0] ,cei);       
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8124_ymm4r8_u(const  double * __restrict  pk0,
                                                 const  double * __restrict pgam,
                                                 double * __restrict  D1r,
                                                 double * __restrict  D1i,
                                                 double * __restrict  D2r,
                                                 double * __restrict  D2i) {
                                    
                        register __m256d k0 = _mm256_loadu_pd(&pk0[0]);
                        register __m256d gam= _mm256_loadu_pd(&pgam[0]);           
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);
                        const __m256d C10                              = _mm256_set1_pd(1.0f);  
                        register __m256d ear,eai,cer,cei,t0r,t0i;
                        register __m256d x0,invn,x1,sin1,cos1,sin2,cos2,x2,invn2,arg1,arg2;
                        ear  = _mm256_setzero_pd();
                        eai  = C078539816339744830961566084582; 
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328);  
                        x1   = _mm256_mul_pd(k0,C6283185307179586476925286766559);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        invn = _mm256_rcp14_pd(x0);
                        cer  = gms::math::negate_ymm4r8(cer);
                        invn2= _mm256_add_pd(invn,invn);      
                        cei  = gms::math::negate_ymm4r8(cei);
                        arg1 = _mm256_mul_pd(C314159265358979323846264338328,invn2)
                        sin1 = _mm256_sin_pd(arg1);
                        arg2 = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        t0r  = _mm256_div_pd(cer,x1);
                        cos1 = _mm256_cos_pd(arg1);
                        x0   = _mm256_div_pd(cos1,sin1);
                        t0i  = _mm256_div_pd(cei,x1);   
                        sin2 = _mm256_sin_pd(arg2);
                        cos2 = _mm256_cos_pd(arg2);   
                        x1   = _mm256_div_pd(cos2,sin2);
                        ear  = _mm256_fmsub_pd(invn,x0,_mm256_mul_pd(invn2,x1));
                        eai  = _mm256_fmadd_pd(invn,x0,_mm256_mul_pd(invn2,x1));
                        cer  = _mm256_mul_pd(t0r,ear);
                        cei  = _mm256_mul_pd(t0i,eai);
                        _mm256_storeu_pd(&D1r[0] ,cer);
                        _mm256_storeu_pd(&D2r[0] ,cer);
                        _mm256_storeu_pd(&D1i[0] ,cei);
                        _mm256_storeu_pd(&D2i[0] ,cei);       
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
                   void coef_Ddiff_f8126_ymm4r8(const __m256d gam,
                                             const __m256d phi,
                                             const __m256d k0,
                                             __m256d * __restrict Dr,
                                             __m256d * __restrict Di) {
                                             
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);  
                        const __m256d C20                              = 
                                                     _mm256_set1_pd(2.0f);
                        register __m256d ear,eai,cer,cei;
                        register __m256d invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm256_add_pd(phi,phi);
                        x0   = _mm256_mul_pd(C6283185307179586476925286766559,k0);
                        ear  = _mm256_setzero_pd();
                        n    = _mm256_div_pd(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm256_rcp14_pd(n); 
                        sqr  = _mm256_sqrt_pd(x0);
                        pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        spin = _mm256_sin_pd(pin);
                        cer  = _mm256_div_pd(cer,sqr);
                        cpin = _mm256_cos_pd(pin);
                        cei  = _mm256_div_pd(cei,sqr);
                        cphin= _mm256_cos_pd(_mm256_mul_pd(phi2,invn));
                        x0   = _mm256_mul_pd(_mm256_mul_pd(C20,invn),spin);
                        x1   = _mm256_sub_pd(cpin,cphin);
                        n    = _mm256_div_pd(x0,x1);
                        *Dr  = _mm256_mul_pd(cer,n);
                        *Di  = _mm256_mul_pd(cei,n);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_Ddiff_f8126_ymm4r8_a(const double * __restrict pgam,
                                               const double * __restrict pphi,
                                               const double * __restrict pk0,
                                               double * __restrict  Dr,
                                               double * __restrict  Di) {
                                    
                        register __m256d gam  = _mm256_load_pd(&pgam[0]);
                        register __m256d phi  = _mm256_load_pd(&pphi[0]);
                        register __m256d k0   = _mm256_load_pd(&pk0[0]);         
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);  
                        const __m256d C20                              = 
                                                     _mm256_set1_pd(2.0f);
                        register __m256d ear,eai,cer,cei;
                        register __m256d invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm256_add_pd(phi,phi);
                        x0   = _mm256_mul_pd(C6283185307179586476925286766559,k0);
                        ear  = _mm256_setzero_pd();
                        n    = _mm256_div_pd(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm256_rcp14_pd(n); 
                        sqr  = _mm256_sqrt_pd(x0);
                        pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        spin = _mm256_sin_pd(pin);
                        cer  = _mm256_div_pd(cer,sqr);
                        cpin = _mm256_cos_pd(pin);
                        cei  = _mm256_div_pd(cei,sqr);
                        cphin= _mm256_cos_pd(_mm256_mul_pd(phi2,invn));
                        x0   = _mm256_mul_pd(_mm256_mul_pd(C20,invn),spin);
                        x1   = _mm256_sub_pd(cpin,cphin);
                        n    = _mm256_div_pd(x0,x1);
                        _mm256_store_pd(&Dr[0] ,_mm256_mul_pd(cer,n));
                        _mm256_store_pd(&Di[0] ,_mm256_mul_pd(cei,n));
                }
                
                
                
                    __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_Ddiff_f8126_ymm4r8_u(const double * __restrict  pgam,
                                                   const double * __restrict pphi,
                                                   const double * __restrict  pk0,
                                                   double * __restrict  Dr,
                                                   double * __restrict  Di) {
                                    
                        register __m256d gam  = _mm256_loadu_pd(&pgam[0]);
                        register __m256d phi  = _mm256_loadu_pd(&pphi[0]);
                        register __m256d k0   = _mm256_loadu_pd(&pk0[0]);         
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328);  
                        const __m256d C20                              = 
                                                     _mm256_set1_pd(2.0f);
                        register __m256d ear,eai,cer,cei;
                        register __m256d invn,n,phi2,pin,spin,cpin,cphin,sqr,x0,x1;
                        phi2 = _mm256_add_pd(phi,phi);
                        x0   = _mm256_mul_pd(C6283185307179586476925286766559,k0);
                        ear  = _mm256_setzero_pd();
                        n    = _mm256_div_pd(gam,C314159265358979323846264338328);
                        eai  = C078539816339744830961566084582;
                        invn = _mm256_rcp14_pd(n); 
                        sqr  = _mm256_sqrt_pd(x0);
                        pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        spin = _mm256_sin_pd(pin);
                        cer  = _mm256_div_pd(cer,sqr);
                        cpin = _mm256_cos_pd(pin);
                        cei  = _mm256_div_pd(cei,sqr);
                        cphin= _mm256_cos_pd(_mm256_mul_pd(phi2,invn));
                        x0   = _mm256_mul_pd(_mm256_mul_pd(C20,invn),spin);
                        x1   = _mm256_sub_pd(cpin,cphin);
                        n    = _mm256_div_pd(x0,x1);
                        _mm256_storeu_pd(&Dr[0] ,_mm256_mul_pd(cer,n));
                        _mm256_storeu_pd(&Di[0] ,_mm256_mul_pd(cei,n));
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
                   void EsHs_f8125_ymm4r8(const __m256d a,
                                           const __m256d k0,
                                           const __m256d r,
                                           const __m256d gam,
                                           const __m256d phi,
                                           const __m256d psi,
                                           __m256d * __restrict Esr,
                                           __m256d * __restrict Esi,
                                           __m256d * __restrict Hsr,
                                           __m256d * __restrict Hsi) {
                                           
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559); 
                        const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                        register __m256d Dr,Di,x0,x1;
                        register __m256d ear,eai,cer,cei;
                        register __m256d k0r,invr,t0r,t0i;
                        register __m256d sqr;
                        k0r  = _mm256_mul_pd(k0,r);
                        ear  = _mm256_setzero_pd();
                        eai  = _mm256_add_pd(k0r,_mm256_sub_pd(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_ymm4r8(gam,phi,k0,Dr,Di); 
                        invr = _mm256_rcp_pd(r);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(C6283185307179586476925286766559,k0));
                        t0r  = _mm256_mul_pd(Dr,C05);
                        cer  = _mm256_mul_pd(cer,invr);
                        sqr  = _mm256_mul_pd(a,sqr);
                        t0i  = _mm256_mul_pd(Di,C05);  
                        cei  = _mm256_mul_pd(cei,invr);
                        x0   = _mm256_mul_pd(t0r,_mm256_mul_pd(sqr,cer));
                        x1   = _mm256_mul_pd(t0i,_mm256_mul_pd(sqr,cei));
                        *Esr = gms::math::negate_ymm4r8(x0);
                        *Esi = gms::math::negate_ymm4r8(x1); 
                        *Hsr = x0;
                        *Hsi = x1;             
                      
                }
                
                
                  __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f8125_ymm4r8_a(const double * __restrict pa,
                                             const double * __restrict pk0,
                                             const double * __restrict pr,
                                             const double * __restrict pgam,
                                             const double * __restrict pphi,
                                             const double * __restrict ppsi,
                                             double * __restrict Esr,
                                             double * __restrict Esi,
                                             double * __restrict Hsr,
                                             double * __restrict Hsi) {
                            
                        register __m256d a   = _mm256_load_pd(&pa[0]);
                        register __m256d k0  = _mm256_load_pd(&pk0[0]);   
                        register __m256d r   = _mm256_load_pd(&pr[0]);   
                        register __m256d gam = _mm256_load_pd(&pgam[0]);   
                        register __m256d phi = _mm256_load_pd(&pphi[0]); 
                        register __m256d psi = _mm256_load_pd(&ppsi[0]);     
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559); 
                        const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                        register __m256d Dr,Di,x0,x1;
                        register __m256d ear,eai,cer,cei;
                        register __m256d k0r,invr,t0r,t0i;
                        register __m256d sqr;
                        k0r  = _mm256_mul_pd(k0,r);
                        ear  = _mm256_setzero_pd();
                        eai  = _mm256_add_pd(k0r,_mm256_sub_pd(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_ymm4r8(gam,phi,k0,Dr,Di); 
                        invr = _mm256_rcp_pd(r);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(C6283185307179586476925286766559,k0));
                        t0r  = _mm256_mul_pd(Dr,C05);
                        cer  = _mm256_mul_pd(cer,invr);
                        sqr  = _mm256_mul_pd(a,sqr);
                        t0i  = _mm256_mul_pd(Di,C05);  
                        cei  = _mm256_mul_pd(cei,invr);
                        x0   = _mm256_mul_pd(t0r,_mm256_mul_pd(sqr,cer));
                        x1   = _mm256_mul_pd(t0i,_mm256_mul_pd(sqr,cei));
                        _mm256_store_pd(&Esr[0] ,gms::math::negate_ymm4r8(x0));
                        _mm256_store_pd(&Esi[0] ,gms::math::negate_ymm4r8(x1)); 
                        _mm256_store_pd(&Hsr[0] ,x0);
                        _mm256_store_pd(&Hsi[0] ,x1);             
                      
                }
                
            
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void EsHs_f8125_ymm4r8_u(const double * __restrict  pa,
                                             const double * __restrict  pk0,
                                             const double * __restrict  pr,
                                             const double * __restrict pgam,
                                             const double * __restrict  pphi,
                                             const double * __restrict  ppsi,
                                             double * __restrict  Esr,
                                             double * __restrict  Esi,
                                             double * __restrict  Hsr,
                                             double * __restrict  Hsi) {
                            
                        register __m256d a   = _mm256_loadu_pd(&pa[0]);
                        register __m256d k0  = _mm256_loadu_pd(&pk0[0]);   
                        register __m256d r   = _mm256_loadu_pd(&pr[0]);   
                        register __m256d gam = _mm256_loadu_pd(&pgam[0]);   
                        register __m256d phi = _mm256_loadu_pd(&pphi[0]); 
                        register __m256d psi = _mm256_loadu_pd(&ppsi[0]);     
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559); 
                        const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                        register __m256d Dr,Di,x0,x1;
                        register __m256d ear,eai,cer,cei;
                        register __m256d k0r,invr,t0r,t0i;
                        register __m256d sqr;
                        k0r  = _mm256_mul_pd(k0,r);
                        ear  = _mm256_setzero_pd();
                        eai  = _mm256_add_pd(k0r,_mm256_sub_pd(psi,
                                                   C078539816339744830961566084582));
                        coef_Ddiff_f8126_ymm4r8(gam,phi,k0,Dr,Di); 
                        invr = _mm256_rcp_pd(r);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(C6283185307179586476925286766559,k0));
                        t0r  = _mm256_mul_pd(Dr,C05);
                        cer  = _mm256_mul_pd(cer,invr);
                        sqr  = _mm256_mul_pd(a,sqr);
                        t0i  = _mm256_mul_pd(Di,C05);  
                        cei  = _mm256_mul_pd(cei,invr);
                        x0   = _mm256_mul_pd(t0r,_mm256_mul_pd(sqr,cer));
                        x1   = _mm256_mul_pd(t0i,_mm256_mul_pd(sqr,cei));
                        _mm256_storeu_pd(&Esr[0] ,gms::math::negate_ymm4r8(x0));
                        _mm256_storeu_pd(&Esi[0] ,gms::math::negate_ymm4r8(x1)); 
                        _mm256_storeu_pd(&Hsr[0] ,x0);
                        _mm256_storeu_pd(&Hsi[0] ,x1);             
                      
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
                   void coef_D12_f8127_ymm4r8(const __m256d k0,
                                               const __m256d gam,
                                               const __m256d phi1,
                                               const __m256d phi2,
                                               __m256d * __restrict D1r,
                                               __m256d * __restrict D1i,
                                               __m256d * __restrict D2r,
                                               __m256d * __restrict D2i) {
                                               
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328); 
                        register __m256d invn,x0,x1,ear,eai,cer,cei;
                        register __m256d sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328);
                        phid = _mm256_mul_pd(_mm256_sub_pd(phi1,phi2),invn);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm256_setzero_pd();
                        phis = _mm256_mul_pd(_mm256_add_pd(phi1,phi2),invn);
                        invn = _mm256_rcp14_pd(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        x0   = _mm256_sin_pd(pin);
                        cer  = _mm256_div_pd(cer,sqr);
                        spin = _mm256_mul_pd(x0,invn);
                        cei  = _mm256_div_pd(cei,sqr);
                        cpin = _mm256_cos_pd(pin);
                        x0   = _mm256_cos_pd(phid);
                        invc1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,x0));
                        x1   = _mm256_cos_pd(phis);
                        invc2= _mm256_rcp14_pd(_mm256_sub_pd(cpin,x1));
                        ear  = _mm256_mul_pd(cer,spin);
                        phis = _mm256_sub_pd(invc1,invc2);
                        eai  = _mm256_mul_pd(cei,spin);
                        phid = _mm256_add_pd(invc1,invc2);
                        *D1r = _mm256_mul_pd(ear,phis);
                        *D1i = _mm256_mul_pd(eai,phis);
                        *D2r = _mm256_mul_pd(ear,phid);
                        *D2i = _mm256_mul_pd(eai,phid);
                }
                
                
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8127_ymm4r8_a(const double * __restrict pk0,
                                                 const double * __restrict pgam,
                                                 const double * __restrict pphi1,
                                                 const double * __restrict pphi2,
                                                 double * __restrict  D1r,
                                                 double * __restrict  D1i,
                                                 double * __restrict  D2r,
                                                 double * __restrict  D2i) {
                              
                        register __m256d k0  = _mm256_load_pd(&pk0[0]);
                        register __m256d gam = _mm256_load_pd(&pgam[0]);
                        register __m256d phi1= _mm256_load_pd(&pphi1[0]);
                        register __m256d phi2= _mm256_load_pd(&pphi2[0]);                 
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328); 
                        register __m256d invn,x0,x1,ear,eai,cer,cei;
                        register __m256d sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328);
                        phid = _mm256_mul_pd(_mm256_sub_pd(phi1,phi2),invn);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm256_setzero_pd();
                        phis = _mm256_mul_pd(_mm256_add_pd(phi1,phi2),invn);
                        invn = _mm256_rcp14_pd(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        x0   = _mm256_sin_pd(pin);
                        cer  = _mm256_div_pd(cer,sqr);
                        spin = _mm256_mul_pd(x0,invn);
                        cei  = _mm256_div_pd(cei,sqr);
                        cpin = _mm256_cos_pd(pin);
                        x0   = _mm256_cos_pd(phid);
                        invc1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,x0));
                        x1   = _mm256_cos_pd(phis);
                        invc2= _mm256_rcp14_pd(_mm256_sub_pd(cpin,x1));
                        ear  = _mm256_mul_pd(cer,spin);
                        phis = _mm256_sub_pd(invc1,invc2);
                        eai  = _mm256_mul_pd(cei,spin);
                        phid = _mm256_add_pd(invc1,invc2);
                        _mm256_store_pd(&D1r[0] ,_mm256_mul_pd(ear,phis));
                        _mm256_store_pd(&D1i[0] ,_mm256_mul_pd(eai,phis));
                        _mm256_store_pd(&D2r[0] ,_mm256_mul_pd(ear,phid));
                        _mm256_store_pd(&D2i[0] ,_mm256_mul_pd(eai,phid));
                }
                
                
                
                  __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   void coef_D12_f8127_ymm4r8_u(const double * __restrict  pk0,
                                                 const double * __restrict  pgam,
                                                 const double * __restrict  pphi1,
                                                 const double * __restrict  pphi2,
                                                 double * __restrict   D1r,
                                                 double * __restrict  D1i,
                                                 double * __restrict   D2r,
                                                 double * __restrict  D2i) {
                              
                        register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
                        register __m256d gam = _mm256_loadu_pd(&pgam[0]);
                        register __m256d phi1= _mm256_loadu_pd(&pphi1[0]);
                        register __m256d phi2= _mm256_loadu_pd(&pphi2[0]);                 
                        const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582);
                        const __m256d C6283185307179586476925286766559 = 
                                                     _mm256_set1_pd(6.283185307179586476925286766559);
                        const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328); 
                        register __m256d invn,x0,x1,ear,eai,cer,cei;
                        register __m256d sqr,pin,spin,cpin,phis,phid,invc1,invc2;
                        x0   = _mm256_div_pd(gam,C314159265358979323846264338328);
                        phid = _mm256_mul_pd(_mm256_sub_pd(phi1,phi2),invn);
                        sqr  = _mm256_sqrt_pd(_mm256_mul_pd(
                                             C6283185307179586476925286766559,k0)); 
                        ear  = _mm256_setzero_pd();
                        phis = _mm256_mul_pd(_mm256_add_pd(phi1,phi2),invn);
                        invn = _mm256_rcp14_pd(x0);
                        eai  = C078539816339744830961566084582;
                        pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                        cexp_ymm4c8(ear,eai,&cer,&cei);
                        x0   = _mm256_sin_pd(pin);
                        cer  = _mm256_div_pd(cer,sqr);
                        spin = _mm256_mul_pd(x0,invn);
                        cei  = _mm256_div_pd(cei,sqr);
                        cpin = _mm256_cos_pd(pin);
                        x0   = _mm256_cos_pd(phid);
                        invc1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,x0));
                        x1   = _mm256_cos_pd(phis);
                        invc2= _mm256_rcp14_pd(_mm256_sub_pd(cpin,x1));
                        ear  = _mm256_mul_pd(cer,spin);
                        phis = _mm256_sub_pd(invc1,invc2);
                        eai  = _mm256_mul_pd(cei,spin);
                        phid = _mm256_add_pd(invc1,invc2);
                        _mm256_storeu_pd(&D1r[0] ,_mm256_mul_pd(ear,phis));
                        _mm256_storeu_pd(&D1i[0] ,_mm256_mul_pd(eai,phis));
                        _mm256_storeu_pd(&D2r[0] ,_mm256_mul_pd(ear,phid));
                        _mm256_storeu_pd(&D2i[0] ,_mm256_mul_pd(eai,phid));
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
                   double rcs_f8162_ymm4r8_u(const double * __restrict pdAdl,
                                             const double *  __restrict pdl,
                                             const double   k0,
                                             const double   l) {
                          
                         double intr[4] = {};
                         double inti[4] = {}; 
                         double Y1[4]   = {};
                         double Y2[4]   = {};
                         double Y3[4]   = {}; 
                         double E[4]    = {};
                         double WRK[4] = {};
                         
                         constexpr int32_t NTAB = 4;               
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                         register __m256d dAdl = _mm256_loadu_pd(&pdAdl[0]);
                         register __m256d dl   = _mm256_loadu_pd(&pdl[0]);                                 
                         register __m256d vk0,k0l,ear,eai,cer,cei;
                         std::complex<double> c;
                         register double rcs,k02,frac,sumr,sumi;
                         vk0  = _mm256_set1_pd(k0);
                         k0l  = _mm256_mul_pd(vk0,dl);
                         ear  = _mm256_setzero_pd();
                         eai  = _mm256_add_pd(k0l,k0l);
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         _mm256_store_pd(&intr[0], _mm256_mul_pd(cer,dAdl);
                         _mm256_store_pd(&inti[0], _mm256_mul_pd(cei,dAdl);
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
                   double rcs_f8162_ymm4r8_a(const double * __restrict pdAdl,
                                             const double * __restrict pdl,
                                             const double   k0,
                                             const double   l) {
                          
                         double intr[4] = {};
                         double inti[4] = {}; 
                         double Y1[4]   = {};
                         double Y2[4]   = {};
                         double Y3[4]   = {}; 
                         double E[4]    = {};
                         double WRK[4] = {};
                         
                         constexpr int32_t NTAB = 4;               
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                         register __m256d dAdl = _mm256_load_pd(&pdAdl[0]);
                         register __m256d dl   = _mm256_load_pd(&pdl[0]);                                 
                         register __m256d vk0,k0l,ear,eai,cer,cei;
                         std::complex<double> c;
                         register double rcs,k02,frac,sumr,sumi;
                         vk0  = _mm256_set1_pd(k0);
                         k0l  = _mm256_mul_pd(vk0,dl);
                         ear  = _mm256_setzero_pd();
                         eai  = _mm256_add_pd(k0l,k0l);
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         _mm256_store_pd(&intr[0], _mm256_mul_pd(cer,dAdl);
                         _mm256_store_pd(&inti[0], _mm256_mul_pd(cei,dAdl);
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
                   double rcs_f8162_ymm4r8_avint_u(const double * __restrict pdAdl,
                                                   const double *  __restrict pdl,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri) {
                          
                         double intr[4] = {};
                         double inti[4] = {}; 
                                         
                         constexpr int32_t NTAB = 4;               
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                         register __m256d dAdl = _mm256_loadu_pd(&pdAdl[0]);
                         register __m256d dl   = _mm256_loadu_pd(&pdl[0]);                                 
                         register __m256d vk0,k0l,ear,eai,cer,cei;
                         std::complex<double> c;
                         register double rcs,k02,frac,sumr,sumi;
                         int32_t err,eri;
                         vk0  = _mm256_set1_pd(k0);
                         k0l  = _mm256_mul_pd(vk0,dl);
                         ear  = _mm256_setzero_pd();
                         eai  = _mm256_add_pd(k0l,k0l);
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         _mm256_store_pd(&intr[0], _mm256_mul_pd(cer,dAdl);
                         _mm256_store_pd(&inti[0], _mm256_mul_pd(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         sumr = avint(pdl,&intr[0],0.0f,l,err);
                         sumi = avint(pdl,&inti[0],0.0f,l,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr == 3 || ieri == 3) {
                            return std::numeric_limits<double>::quiet_NaN();
                         }
                         c = {sumr,sumi};
                         k02   = k0*k0;   
                         frac  = k02/C314159265358979323846264338328;
                         rcs   = frac*std::abs(c);
                         return (rcs);                         
                  }
                  
                  
                   __ATTR_ALWAYS_INLINE__
	          
	          
                   
	           static inline
                   double rcs_f8162_ymm4r8_avint_a(const double * __restrict pdAdl,
                                                   const double * __restrict pdl,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri) {
                          
                         double intr[4] = {};
                         double inti[4] = {}; 
                                         
                         constexpr int32_t NTAB = 4;               
                         constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                         register __m256d dAdl = _mm256_load_pd(&pdAdl[0]);
                         register __m256d dl   = _mm256_load_pd(&pdl[0]);                                 
                         register __m256d vk0,k0l,ear,eai,cer,cei;
                         std::complex<double> c;
                         register double rcs,k02,frac,sumr,sumi;
                         int32_t err,eri;
                         vk0  = _mm256_set1_pd(k0);
                         k0l  = _mm256_mul_pd(vk0,dl);
                         ear  = _mm256_setzero_pd();
                         eai  = _mm256_add_pd(k0l,k0l);
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         _mm256_stor_pd(&intr[0], _mm256_mul_pd(cer,dAdl);
                         _mm256_stor_pd(&inti[0], _mm256_mul_pd(cei,dAdl);
                         sumr = 0.0f;
                         sumi = 0.0f;
                         sumr = avint(pdl,&intr[0],0.0f,l,err);
                         sumi = avint(pdl,&inti[0],0.0f,l,eri);
                         ierr = err;
                         ieri = eri;
                         if(ierr == 3 || ieri == 3) {
                            return std::numeric_limits<double>::quiet_NaN();
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
                   double rcs_f8162_ymm4r8_u(const double * __restrict  pdAdl,
                                             const double * __restrict  pdl,
                                             double * __restrict  intr,
                                             double * __restrict  inti,
                                             double * __restrict  Y1,
                                             double * __restrict  Y2,
                                             double * __restrict  Y3,
                                             double * __restrict  E,
                                             double * __restrict  WRK
                                             const double   k0,
                                             const double   l,
                                             const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            double rcs = 0.0f;
                            rcs = rcs_f8162_ymm4r8_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_loadu_pd(&pdAdl[i]);
                             register __m256d y = _mm256_loadu_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_storeu_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_storeu_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
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
                   double rcs_f8162_ymm4r8_a(const double * __restrict pdAdl,
                                             const double * __restrict pdl,
                                             double * __restrict intr,
                                             double * __restrict inti,
                                             double * __restrict Y1,
                                             double * __restrict Y2,
                                             double * __restrict Y3,
                                             double * __restrict E,
                                             double * __restrict WRK
                                             const double   k0,
                                             const double   l,
                                             const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            double rcs = 0.0f;
                            rcs = rcs_f8162_ymm4r8_a(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_load_pd(&pdAdl[i]);
                             register __m256d y = _mm256_load_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_store_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_store_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
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
                   double rcs_f8162_ymm4r8_avint_u(const double * __restrict  pdAdl,
                                                   const double * __restrict  pdl,
                                                   double * __restrict  intr,
                                                   double * __restrict  inti,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB) {
                                             
                                                
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_loadu_pd(&pdAdl[i]);
                             register __m256d y = _mm256_loadu_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_storeu_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_storeu_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
                      sumr = avint(pdl,intr,NTAB,0.0f,l,err); 
                      sumi = avint(pdl,inti,NTAB,0.0f,l,eri); 
                      ierr = err;
                      ieri = eri;
                      if(ierr == 3 || ieri == 3) {
                         return std::numerical_limits<double>::quiet_NaN();
                      } 
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   double rcs_f8162_ymm4r8_avint_a(const double * __restrict pdAdl,
                                                   const double * __restrict pdl,
                                                   double * __restrict intr,
                                                   double * __restrict inti,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB) {
                                             
                                                
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_load_pd(&pdAdl[i]);
                             register __m256d y = _mm256_load_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_store_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_store_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
                           intr[i]        = c.real()*x;
                           inti[i]        = c.imag()*x;
                      }   
                      sumr = avint(pdl,intr,NTAB,0.0f,l,err); 
                      sumi = avint(pdl,inti,NTAB,0.0f,l,eri); 
                      ierr = err;
                      ieri = eri;
                      if(ierr == 3 || ieri == 3) {
                         return std::numerical_limits<double>::quiet_NaN();
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
                   double rcs_f8162_ymm4r8_cspint2t_u(const double * __restrict  pdAdl,
                                                     const double * __restrict  pdl,
                                                     double * __restrict  intr,
                                                     double * __restrict  inti,
                                                     struct RCS_F8162_DATA w,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            double rcs = 0.0f;
                            rcs = rcs_f8162_ymm4r8_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        double * __restrict px1 = w.Ya1;
                        double * __restrict py1 = w.Yb1;
                        double * __restrict px2 = w.Ya2;
                        double * __restrict py2 = w.Yb2;
                        double * __restrict px3 = w.Ya3;
                        double * __restrict py3 = w.Yb3;
                        double * __restrict px4 = w.Ea;
                        double * __restrict py4 = w.Eb;
                        double * __restrict px5 = w.WRKa;
                        double * __restrict py5 = w.WRKb;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_loadu_pd(&pdAdl[i]);
                             register __m256d y = _mm256_loadu_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_storeu_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_storeu_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
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
                   double rcs_f8162_ymm4r8_cspint2t_a(const double * __restrict pdAdl,
                                                     const double * __restrict pdl,
                                                     double * __restrict intr,
                                                     double * __restrict inti,
                                                     struct RCS_F8162_DATA w,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB) {
                                             
                        if(__builtin_expect(NTAB==4,0)) {
                            double rcs = 0.0f;
                            rcs = rcs_f8162_ymm4r8_u(pdAdl,pdl,k0,l);
                            return (rcs);
                        }    
                        
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        double * __restrict px1 = w.Ya1;
                        double * __restrict py1 = w.Yb1;
                        double * __restrict px2 = w.Ya2;
                        double * __restrict py2 = w.Yb2;
                        double * __restrict px3 = w.Ya3;
                        double * __restrict py3 = w.Yb3;
                        double * __restrict px4 = w.Ea;
                        double * __restrict py4 = w.Eb;
                        double * __restrict px5 = w.WRKa;
                        double * __restrict py5 = w.WRKb;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_load_pd(&pdAdl[i]);
                             register __m256d y = _mm256_load_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_store_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_store_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
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
                   double rcs_f8162_ymm4r8_avint2t_u(const double * __restrict  pdAdl,
                                                     const double * __restrict  pdl,
                                                     double * __restrict  intr,
                                                     double * __restrict  inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB) {
                                             
                                              
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_loadu_pd(&pdAdl[i]);
                             register __m256d y = _mm256_loadu_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_storeu_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_storeu_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
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
                          return std::numeric_limits<double>::quiet_NaN();
                      }
                      c = {sumr,sumi};
                      k02   = k0*k0;   
                      frac  = k02/C314159265358979323846264338328;
                      rcs   = frac*std::abs(c);
                      return (rcs);               
               }
               
               
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   double rcs_f8162_ymm4r8_avint2t_a(const double * __restrict pdAdl,
                                                     const double * __restrict pdl,
                                                     double * __restrict intr,
                                                     double * __restrict inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB) {
                                             
                                              
                        constexpr double C314159265358979323846264338328 = 
                                                        3.14159265358979323846264338328;
                        register __m256d vk0,k0l,ear,eai,cer,cei;
                        std::complex<double> c;
                        register double rcs,k02,frac,sumr,sumi; 
                        int32_t i,err,eri; 
                        vk0  = _mm256_set1_pd(k0);
                        ear  = _mm256_setzero_pd();
                        for(i = 0; i != ROUND_TO_FOUR(NTAB,3); i += 4) {
                             _mm_prefetch((const char*)&pdAdl[i],_MM_NTA_T0);
                             _mm_prefetch((const char*)&pdl[i],  _MM_NTA_T0);
                             register __m256d x = _mm256_load_pd(&pdAdl[i]);
                             register __m256d y = _mm256_load_pd(&pdl[i]);
                             k0l               = _mm256_mul_pd(vk0,y);
                             eai               = _mm256_add_pd(k0l,k0l);
                             cexp_ymm4c8(ear,eai,&cer,&cei);
                             register __m256d t0 = cer;
                             register __m256d t1 = cei;
                             _mm256_store_pd(&intr[i], _mm256_mul_pd(t0,x));
                             _mm256_store_pd(&inti[i], _mm256_mul_pd(t1,x));
                       } 
                       sumr = 0.0f;
                       sumi = 0.0f;
                       for(; i != NTAB; ++i) {
                           const double x  = pdAdl[i];
                           const double y  = pdl[i];
                           const double k0l= k0*y;
                           const double eai= k0l+k0l;
                           const std::complex<double> c = std::exp({0.0f,eai});
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
                          return std::numeric_limits<double>::quiet_NaN();
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
                   __m256d rcs_f8193_ymm4r8(const __m256d b,
                                            const __m256d a,
                                            const __m256d k0,
                                            const __m256d alp,
                                            const __m256d l) {
                                            
                         const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);
                         const __m256d C15                              =
                                                     _mm256_set1_pd(1.5f);
                         register __m256d sina,pin,n,invn,x0,x1,k0l,sqr;
                         register __m256d ear,eai,cer,cei,cosa;
                         register __m256d rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm256_mul_pd(k0,b);
                         ear  = _mm256_setzero_pd();
                         n    = _mm256_add_pd(C15,_mm256_div_pd(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm256_mul_pd(k0,l);
                         eai  = _mm256_add_pd(k0l,k0l);
                         invn = _mm256_rcp14_pd(n);
                         sina = _mm256_sub_pd(C10,_mm256_sin_pd(alp));
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         cosa = _mm256_cos_pd(alp);
                         x0   = _mm256_sin_pd(_mm256_mul_pd(
                                              _mm256_add_pd(k0b,k0b),sina));
                         x1   = _mm256_mul_pd(k0b,_mm256_mul_pd(cosa,cosa));
                         sqr  = _mm256_sub_pd(C10,_mm256_div_pd(x0,x1));
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         x0   = _mm256_sqrt_pd(sqr);
                         spin = _mm256_mul_pd(_mm256_sin_pd(pin),invn);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_add_pd(alp,alp),invn);
                         cos2a= _mm256_cos_pd(x1);
                         t0   = _mm256_mul_pd(C1772453850905516027298167483341,
                                                              _mm256_mul_pd(b,x0)); // keep
                         x1   = _mm256_sub_pd(cpin,cos2a);
                         t1   = _mm256_mul_pd(_mm256_add_pd(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm256_mul_pd(a,spin)); // keep
                         x0   = _mm256_rcp14_pd(x1); // keep
                         ear  = _mm256_mul_pd(_mm256_add_pd(t0,t1),x0);
                         t0   = _mm256_mul_pd(ear,cer);
                         t1   = _mm256_mul_pd(ear,cei);
                         rcs  = cabs_ymm4c8(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m256d rcs_f8193_ymm4r8_a(const double * __restrict pb,
                                              const double * __restrict pa,
                                              const double * __restrict pk0,
                                              const double * __restrict palp,
                                              const double * __restrict pl) {
                                     
                         register __m256d b  = _mm256_load_pd(&pb[0]);
                         register __m256d a  = _mm256_load_pd(&pa[0]);  
                         register __m256d k0 = _mm256_load_pd(&pk0[0]);
                         register __m256d alp= _mm256_load_pd(&palp[0]);  
                         register __m256d l  = _mm256_load_pd(&pl[0]);   
                         const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);
                         const __m256d C15                              =
                                                     _mm256_set1_pd(1.5f);
                         register __m256d sina,pin,n,invn,x0,x1,k0l,sqr;
                         register __m256d ear,eai,cer,cei,cosa;
                         register __m256d rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm256_mul_pd(k0,b);
                         ear  = _mm256_setzero_pd();
                         n    = _mm256_add_pd(C15,_mm256_div_pd(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm256_mul_pd(k0,l);
                         eai  = _mm256_add_pd(k0l,k0l);
                         invn = _mm256_rcp14_pd(n);
                         sina = _mm256_sub_pd(C10,_mm256_sin_pd(alp));
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         cosa = _mm256_cos_pd(alp);
                         x0   = _mm256_sin_pd(_mm256_mul_pd(
                                              _mm256_add_pd(k0b,k0b),sina));
                         x1   = _mm256_mul_pd(k0b,_mm256_mul_pd(cosa,cosa));
                         sqr  = _mm256_sub_pd(C10,_mm256_div_pd(x0,x1));
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         x0   = _mm256_sqrt_pd(sqr);
                         spin = _mm256_mul_pd(_mm256_sin_pd(pin),invn);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_add_pd(alp,alp),invn);
                         cos2a= _mm256_cos_pd(x1);
                         t0   = _mm256_mul_pd(C1772453850905516027298167483341,
                                                              _mm256_mul_pd(b,x0)); // keep
                         x1   = _mm256_sub_pd(cpin,cos2a);
                         t1   = _mm256_mul_pd(_mm256_add_pd(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm256_mul_pd(a,spin)); // keep
                         x0   = _mm256_rcp14_pd(x1); // keep
                         ear  = _mm256_mul_pd(_mm256_add_pd(t0,t1),x0);
                         t0   = _mm256_mul_pd(ear,cer);
                         t1   = _mm256_mul_pd(ear,cei);
                         rcs  = cabs_ymm4c8(t0,t1);
                         return (rcs);                      
                 }
                 
                 
                 
                  __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m256d rcs_f8193_ymm4r8_u(const double * __restrict  pb,
                                              const double * __restrict  pa,
                                              const double * __restrict  pk0,
                                              const double * __restrict  palp,
                                              const double * __restrict _pl) {
                                     
                         register __m256d b  = _mm256_loadu_pd(&pb[0]);
                         register __m256d a  = _mm256_loadu_pd(&pa[0]);  
                         register __m256d k0 = _mm256_loadu_pd(&pk0[0]);
                         register __m256d alp= _mm256_loadu_pd(&palp[0]);  
                         register __m256d l  = _mm256_loadu_pd(&pl[0]);   
                         const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);
                         const __m256d C15                              =
                                                     _mm256_set1_pd(1.5f);
                         register __m256d sina,pin,n,invn,x0,x1,k0l,sqr;
                         register __m256d ear,eai,cer,cei,cosa;
                         register __m256d rcs,cpin,cos2a,k0b,t0,t1;
                         k0b  = _mm256_mul_pd(k0,b);
                         ear  = _mm256_setzero_pd();
                         n    = _mm256_add_pd(C15,_mm256_div_pd(alp,   
                                                          C314159265358979323846264338328));
                         k0l  = _mm256_mul_pd(k0,l);
                         eai  = _mm256_add_pd(k0l,k0l);
                         invn = _mm256_rcp14_pd(n);
                         sina = _mm256_sub_pd(C10,_mm256_sin_pd(alp));
                         cexp_ymm4c8(ear,eai,&cer,&cei);
                         cosa = _mm256_cos_pd(alp);
                         x0   = _mm256_sin_pd(_mm256_mul_pd(
                                              _mm256_add_pd(k0b,k0b),sina));
                         x1   = _mm256_mul_pd(k0b,_mm256_mul_pd(cosa,cosa));
                         sqr  = _mm256_sub_pd(C10,_mm256_div_pd(x0,x1));
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         x0   = _mm256_sqrt_pd(sqr);
                         spin = _mm256_mul_pd(_mm256_sin_pd(pin),invn);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_add_pd(alp,alp),invn);
                         cos2a= _mm256_cos_pd(x1);
                         t0   = _mm256_mul_pd(C1772453850905516027298167483341,
                                                              _mm256_mul_pd(b,x0)); // keep
                         x1   = _mm256_sub_pd(cpin,cos2a);
                         t1   = _mm256_mul_pd(_mm256_add_pd(C1772453850905516027298167483341,
                                                            C1772453850905516027298167483341),
                                                                           _mm256_mul_pd(a,spin)); // keep
                         x0   = _mm256_rcp14_pd(x1); // keep
                         ear  = _mm256_mul_pd(_mm256_add_pd(t0,t1),x0);
                         t0   = _mm256_mul_pd(ear,cer);
                         t1   = _mm256_mul_pd(ear,cei);
                         rcs  = cabs_ymm4c8(t0,t1);
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
                   __m256d rcs_f8196_ymm4r8(const __m256d k0,
                                            const __m256d alp,
                                            const __m256d a,
                                            const __m256d b) {
                                            
                         
                          const __m256d C0444444444444444444444444444444 = 
                                                        _mm256_set1_pd(0.444444444444444444444444444444);
                          register __m256d rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm256_mul_pd(C0444444444444444444444444444444,k0);
                          cosa  = _mm256_cos_pd(alp);
                          a32   = _mm256_mul_pd(a,_mm256_sqrt_pd(a));
                          x0    = _mm256_mul_pd(b,cosa);
                          sina  = _mm256_sin_pd(alp);
                          bca32 = _mm256_mul_pd(x0,_mm256_sqrt_pd(x0));
                          x1    = _mm256_div_pd(cosa,sina);
                          cota  = _mm256_mul_pd(x1,x1);
                          t1    = _mm256_mul_pd(t0,_mm256_mul_pd(cosa,cota));
                          x0    = _mm256_sub_pd(a32,bca32);
                          rcs   = _mm256_mul_pd(t1,_mm256_mul_pd(x0,x0));
                          return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m256d rcs_f8196_ymm4r8_a(const double * __restrict pk0,
                                              const double * __restrict palp,
                                              const double * __restrict pa,
                                              const double * __restrict pb) {
                                            
                         
                          register __m256d k0  = _mm256_load_pd(&pk0[0]);
                          register __m256d alp = _mm256_load_pd(&palp[0]);
                          register __m256d a   = _mm256_load_pd(&pa[0]);
                          register __m256d b   = _mm256_load_pd(&pb[0]);
                          const __m256d C0444444444444444444444444444444 = 
                                                        _mm256_set1_pd(0.444444444444444444444444444444);
                          register __m256d rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm256_mul_pd(C0444444444444444444444444444444,k0);
                          cosa  = _mm256_cos_pd(alp);
                          a32   = _mm256_mul_pd(a,_mm256_sqrt_pd(a));
                          x0    = _mm256_mul_pd(b,cosa);
                          sina  = _mm256_sin_pd(alp);
                          bca32 = _mm256_mul_pd(x0,_mm256_sqrt_pd(x0));
                          x1    = _mm256_div_pd(cosa,sina);
                          cota  = _mm256_mul_pd(x1,x1);
                          t1    = _mm256_mul_pd(t0,_mm256_mul_pd(cosa,cota));
                          x0    = _mm256_sub_pd(a32,bca32);
                          rcs   = _mm256_mul_pd(t1,_mm256_mul_pd(x0,x0));
                          return (rcs);
                 }
                 
                 
                   __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
                   __m256d rcs_f8196_ymm4r8_u(const double * __restrict pk0,
                                              const double * __restrict palp,
                                              const double * __restrict pa,
                                              const double * __restrict pb) {
                                            
                         
                          register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
                          register __m256d alp = _mm256_loadu_pd(&palp[0]);
                          register __m256d a   = _mm256_loadu_pd(&pa[0]);
                          register __m256d b   = _mm256_loadu_pd(&pb[0]);
                          const __m256d C0444444444444444444444444444444 = 
                                                        _mm256_set1_pd(0.444444444444444444444444444444);
                          register __m256d rcs,t0,cosa,cota,x0,x1,a32,bca32,sina,t1;
                          t0    = _mm256_mul_pd(C0444444444444444444444444444444,k0);
                          cosa  = _mm256_cos_pd(alp);
                          a32   = _mm256_mul_pd(a,_mm256_sqrt_pd(a));
                          x0    = _mm256_mul_pd(b,cosa);
                          sina  = _mm256_sin_pd(alp);
                          bca32 = _mm256_mul_pd(x0,_mm256_sqrt_pd(x0));
                          x1    = _mm256_div_pd(cosa,sina);
                          cota  = _mm256_mul_pd(x1,x1);
                          t1    = _mm256_mul_pd(t0,_mm256_mul_pd(cosa,cota));
                          x0    = _mm256_sub_pd(a32,bca32);
                          rcs   = _mm256_mul_pd(t1,_mm256_mul_pd(x0,x0));
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
	           _
	          
                   
	           static inline
	           __m256d rcs_perpendic_f8194_ymm4r8(const __m256d h,
	                                        const __m256d l,
	                                        const __m256d b,
	                                        const __m256d a,
	                                        const __m256d k0,
	                                        const __m256d tht,
	                                        const __m256d alp) {
	                                 
	                                  
	                 const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328f); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341f);
                         const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582f);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);  
                         const __m256d C15                              = 
                                                     _mm256_set1_pd(1.5f); 
                         const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                         const __m256d C20                              =
                                                     _mm256_set1_pd(2.0f);
                         register __m256d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m256d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m256d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m256d cpin1,cpin2,trm1,trm2,rcs;
                         __m256d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_pd(h,_mm256_add_pd(l,b));
                         sint = _mm256_sin_pd(tht);
                         k02  = _mm256_add_pd(k0,k0);
                         n    = _mm256_mul_pd(C15,_mm256_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_pd(sint);
                         a2   = _mm256_mul_pd(a,C05);
                         ear  = _mm256_setzero_pd();
                         sk02 = _mm256_sqrt_pd(_mm256_mul_pd(k0,C05));
                         x0   = _mm256_mul_pd(hlb,_mm256_sub_pd(cost,b));
                         invn = _mm256_rcp14_pd(n);
                         //x2   = _mm256_mul_pd(a,C05);
                         eai  = _mm256_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_pd(_mm256_mul_pd(a2,csct));
                         atant= _mm256_mul_pd(a,tant);
                         cost = _mm256_cos_pd(tht);
                         x0   = _mm256_mul_pd(b,C1772453850905516027298167483341);
                         cexp_ymm4c8(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_pd(x0,cer1);
                         spin = _mm256_sin_pd(pin);
                         cei1 = _mm256_mul_pd(x0,cei1);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_sub_pd(h,atant),cost);
                         eai2 = _mm256_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_ymm4c8(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_mul_pd(x0,sacs);
                         cer2 = _mm256_mul_pd(cer2,x1);
                         cei2 = _mm256_mul_pd(cei2,x1);
                         cpin1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,C10));
                         x2   = _mm256_mul_pd(C20,_mm256_add_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x2,invn));
                         x3   = _mm256_rcp14_pd(_mm256_sub_pd(cpin,cpin2));
                         trm1 = _mm256_sub_pd(cpin1,x3);
                         cmul_ymm4c8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_pd(t0r,trm1);
                         t0i  = _mm256_mul_pd(t0i,trm1);
                         x0   = _mm256_mul_pd(C20,_mm256_sub_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x0,invn));
                         x1   = _mm256_rcp14_pd(cpin2);
                         trm2 = _mm256_sub_pd(cpin1,x1);
                         x2   = _mm256_fmadd_pd(cost,_mm256_mul_pd(k02,
                                                               _mm256_add_pd(h,atant)));
                         eai3 = _mm256_add_pd(C078539816339744830961566084582,x2);
                         cexp_ymm4c8(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_sqrt_pd(_mm256_mul_pd(gms::math::
                                                                  negate_ymm4r8(a2),csct));
                         x2   = _mm256_mul_pd(x0,x1);
                         cer3 = _mm256_mul_pd(_mm256_mul_pd(cer3,x2),trm2);
                         cei3 = _mm256_mul_pd(_mm256_mul_pd(cei3,x2),trm2);
                         cmul_ymm4c8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm4c8(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	                                  
                    __ATTR_ALWAYS_INLINE__
	           _
	          
                   
	           static inline
	           __m256d rcs_perpendic_f8194_ymm4r8_a(const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pb,
	                                         const double * __restrict pa,
	                                         const double * __restrict pk0,
	                                         const double * __restrict ptht,
	                                         const double * __restrict palp) {
	                                 
	                  
	                 register __m256d h  = _mm256_load_pd(&ph[0]);
	                 register __m256d l  = _mm256_load_pd(&pl[0]); 
	                 register __m256d b  = _mm256_load_pd(&pb[0]);   
	                 register __m256d a  = _mm256_load_pd(&pa[0]);  
	                 register __m256d k0 = _mm256_load_pd(&pk0[0]);
	                 register __m256d tht= _mm256_load_pd(&ptht[0]); 
	                 register __m256d alp= _mm256_load_pd(&palp[0]);        
	                 const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328f); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341f);
                         const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582f);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);  
                         const __m256d C15                              = 
                                                     _mm256_set1_pd(1.5f); 
                         const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                         const __m256d C20                              =
                                                     _mm256_set1_pd(2.0f);
                         __m256d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m256d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m256d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m256d cpin1,cpin2,trm1,trm2,rcs;
                         __m256d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_pd(h,_mm256_add_pd(l,b));
                         sint = _mm256_sin_pd(tht);
                         k02  = _mm256_add_pd(k0,k0);
                         n    = _mm256_mul_pd(C15,_mm256_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_pd(sint);
                         a2   = _mm256_mul_pd(a,C05);
                         ear  = _mm256_setzero_pd();
                         sk02 = _mm256_sqrt_pd(_mm256_mul_pd(k0,C05));
                         x0   = _mm256_mul_pd(hlb,_mm256_sub_pd(cost,b));
                         invn = _mm256_rcp14_pd(n);
                         //x2   = _mm256_mul_pd(a,C05);
                         eai  = _mm256_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_pd(_mm256_mul_pd(a2,csct));
                         atant= _mm256_mul_pd(a,tant);
                         cost = _mm256_cos_pd(tht);
                         x0   = _mm256_mul_pd(b,C1772453850905516027298167483341);
                         cexp_ymm4c8(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_pd(x0,cer1);
                         spin = _mm256_sin_pd(pin);
                         cei1 = _mm256_mul_pd(x0,cei1);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_sub_pd(h,atant),cost);
                         eai2 = _mm256_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_ymm4c8(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_mul_pd(x0,sacs);
                         cer2 = _mm256_mul_pd(cer2,x1);
                         cei2 = _mm256_mul_pd(cei2,x1);
                         cpin1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,C10));
                         x2   = _mm256_mul_pd(C20,_mm256_add_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x2,invn));
                         x3   = _mm256_rcp14_pd(_mm256_sub_pd(cpin,cpin2));
                         trm1 = _mm256_sub_pd(cpin1,x3);
                         cmul_ymm4c8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_pd(t0r,trm1);
                         t0i  = _mm256_mul_pd(t0i,trm1);
                         x0   = _mm256_mul_pd(C20,_mm256_sub_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x0,invn));
                         x1   = _mm256_rcp14_pd(cpin2);
                         trm2 = _mm256_sub_pd(cpin1,x1);
                         x2   = _mm256_fmadd_pd(cost,_mm256_mul_pd(k02,
                                                               _mm256_add_pd(h,atant)));
                         eai3 = _mm256_add_pd(C078539816339744830961566084582,x2);
                         cexp_ymm4c8(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_sqrt_pd(_mm256_mul_pd(gms::math::
                                                                  negate_ymm4r8(a2),csct));
                         x2   = _mm256_mul_pd(x0,x1);
                         cer3 = _mm256_mul_pd(_mm256_mul_pd(cer3,x2),trm2);
                         cei3 = _mm256_mul_pd(_mm256_mul_pd(cei3,x2),trm2);
                         cmul_ymm4c8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm4c8(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           _
	          
                   
	           static inline
	           __m256d rcs_perpendic_f8194_ymm4r8_u(    const double * __restrict  ph,
	                                                    const double * __restrict  pl,
	                                                    const double * __restrict  pb,
	                                                    const double * __restrict  pa,
	                                                    const double * __restrict  pk0,
	                                                    const double * __restrict  ptht,
	                                                    const double * __restrict  palp) {
	                                 
	                  
	                 register __m256d h  = _mm256_loadu_pd(&ph[0]);
	                 register __m256d l  = _mm256_loadu_pd(&pl[0]); 
	                 register __m256d b  = _mm256_loadu_pd(&pb[0]);   
	                 register __m256d a  = _mm256_loadu_pd(&pa[0]);  
	                 register __m256d k0 = _mm256_loadu_pd(&pk0[0]);
	                 register __m256d tht= _mm256_loadu_pd(&ptht[0]); 
	                 register __m256d alp= _mm256_loadu_pd(&palp[0]);        
	                 const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328f); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341f);
                         const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582f);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);  
                         const __m256d C15                              = 
                                                     _mm256_set1_pd(1.5f); 
                         const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                         const __m256d C20                              =
                                                     _mm256_set1_pd(2.0f);
                         __m256d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m256d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m256d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m256d cpin1,cpin2,trm1,trm2,rcs;
                         __m256d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_pd(h,_mm256_add_pd(l,b));
                         sint = _mm256_sin_pd(tht);
                         k02  = _mm256_add_pd(k0,k0);
                         n    = _mm256_mul_pd(C15,_mm256_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_pd(sint);
                         a2   = _mm256_mul_pd(a,C05);
                         ear  = _mm256_setzero_pd();
                         sk02 = _mm256_sqrt_pd(_mm256_mul_pd(k0,C05));
                         x0   = _mm256_mul_pd(hlb,_mm256_sub_pd(cost,b));
                         invn = _mm256_rcp14_pd(n);
                         //x2   = _mm256_mul_pd(a,C05);
                         eai  = _mm256_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_pd(_mm256_mul_pd(a2,csct));
                         atant= _mm256_mul_pd(a,tant);
                         cost = _mm256_cos_pd(tht);
                         x0   = _mm256_mul_pd(b,C1772453850905516027298167483341);
                         cexp_ymm4c8(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_pd(x0,cer1);
                         spin = _mm256_sin_pd(pin);
                         cei1 = _mm256_mul_pd(x0,cei1);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_sub_pd(h,atant),cost);
                         eai2 = _mm256_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_ymm4c8(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_mul_pd(x0,sacs);
                         cer2 = _mm256_mul_pd(cer2,x1);
                         cei2 = _mm256_mul_pd(cei2,x1);
                         cpin1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,C10));
                         x2   = _mm256_mul_pd(C20,_mm256_add_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x2,invn));
                         x3   = _mm256_rcp14_pd(_mm256_sub_pd(cpin,cpin2));
                         trm1 = _mm256_sub_pd(cpin1,x3);
                         cmul_ymm4c8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_pd(t0r,trm1);
                         t0i  = _mm256_mul_pd(t0i,trm1);
                         x0   = _mm256_mul_pd(C20,_mm256_sub_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x0,invn));
                         x1   = _mm256_rcp14_pd(cpin2);
                         trm2 = _mm256_sub_pd(cpin1,x1);
                         x2   = _mm256_fmadd_pd(cost,_mm256_mul_pd(k02,
                                                               _mm256_add_pd(h,atant)));
                         eai3 = _mm256_add_pd(C078539816339744830961566084582,x2);
                         cexp_ymm4c8(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_sqrt_pd(_mm256_mul_pd(gms::math::
                                                                  negate_ymm4r8(a2),csct));
                         x2   = _mm256_mul_pd(x0,x1);
                         cer3 = _mm256_mul_pd(_mm256_mul_pd(cer3,x2),trm2);
                         cei3 = _mm256_mul_pd(_mm256_mul_pd(cei3,x2),trm2);
                         cmul_ymm4c8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm4c8(t1r,t1i);
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
	           _
	          
                   
	           static inline
	           __m256d rcs_parallel_f8194_ymm4r8(const __m256d h,
	                                             const __m256d l,
	                                             const __m256d b,
	                                             const __m256d a,
	                                             const __m256d k0,
	                                             const __m256d tht,
	                                             const __m256d alp) {
	                                 
	                                  
	                 const __m256d C314159265358979323846264338328  = 
                                                     _mm256_set1_pd(3.14159265358979323846264338328f); 
                         const __m256d C1772453850905516027298167483341 = 
                                                     _mm256_set1_pd(1.772453850905516027298167483341f);
                         const __m256d C078539816339744830961566084582  = 
                                                     _mm256_set1_pd(0.78539816339744830961566084582f);
                         const __m256d C10                              = 
                                                     _mm256_set1_pd(1.0f);  
                         const __m256d C15                              = 
                                                     _mm256_set1_pd(1.5f); 
                         const __m256d C05                              = 
                                                     _mm256_set1_pd(0.5f);
                         const __m256d C20                              =
                                                     _mm256_set1_pd(2.0f);
                         register __m256d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m256d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m256d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m256d cpin1,cpin2,trm1,trm2,rcs;
                         __m256d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_pd(h,_mm256_add_pd(l,b));
                         sint = _mm256_sin_pd(tht);
                         k02  = _mm256_add_pd(k0,k0);
                         n    = _mm256_mul_pd(C15,_mm256_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_pd(sint);
                         a2   = _mm256_mul_pd(a,C05);
                         ear  = _mm256_setzero_pd();
                         sk02 = _mm256_sqrt_pd(_mm256_mul_pd(k0,C05));
                         x0   = _mm256_mul_pd(hlb,_mm256_sub_pd(cost,b));
                         invn = _mm256_rcp14_pd(n);
                         //x2   = _mm256_mul_pd(a,C05);
                         eai  = _mm256_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_pd(_mm256_mul_pd(a2,csct));
                         atant= _mm256_mul_pd(a,tant);
                         cost = _mm256_cos_pd(tht);
                         x0   = _mm256_mul_pd(b,C1772453850905516027298167483341);
                         cexp_ymm4c8(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_pd(x0,cer1);
                         spin = _mm256_sin_pd(pin);
                         cei1 = _mm256_mul_pd(x0,cei1);
                         cpin = _mm256_cos_pd(pin);
                         x1   = _mm256_mul_pd(_mm256_sub_pd(h,atant),cost);
                         eai2 = _mm256_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_ymm4c8(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_mul_pd(x0,sacs);
                         cer2 = _mm256_mul_pd(cer2,x1);
                         cei2 = _mm256_mul_pd(cei2,x1);
                         cpin1= _mm256_rcp14_pd(_mm256_sub_pd(cpin,C10));
                         x2   = _mm256_mul_pd(C20,_mm256_add_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x2,invn));
                         x3   = _mm256_rcp14_pd(_mm256_sub_pd(cpin,cpin2));
                         trm1 = _mm256_add_pd(cpin1,x3);
                         cmul_ymm4c8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_pd(t0r,trm1);
                         t0i  = _mm256_mul_pd(t0i,trm1);
                         x0   = _mm256_mul_pd(C20,_mm256_sub_pd(alp,tht));
                         cpin2= _mm256_cos_pd(_mm256_mul_pd(x0,invn));
                         x1   = _mm256_rcp14_pd(cpin2);
                         trm2 = _mm256_add_pd(cpin1,x1);
                         x2   = _mm256_fmadd_pd(cost,_mm256_mul_pd(k02,
                                                               _mm256_add_pd(h,atant)));
                         eai3 = _mm256_add_pd(C078539816339744830961566084582,x2);
                         cexp_ymm4c8(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_pd(spin,_mm256_mul_pd(n,sk02));
                         x1   = _mm256_sqrt_pd(_mm256_mul_pd(gms::math::
                                                                  negate_ymm4r8(a2),csct));
                         x2   = _mm256_mul_pd(x0,x1);
                         cer3 = _mm256_mul_pd(_mm256_mul_pd(cer3,x2),trm2);
                         cei3 = _mm256_mul_pd(_mm256_mul_pd(cei3,x2),trm2);
                         cmul_ymm4c8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm4c8(t1r,t1i);
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
	           __m256d coef_Bg_f9137_ymm4r8(const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol) {
	                                        
	                  const __m256d C06          = _mm256_set1_pd(0.6f);
	                  const __m256d C01875       = _mm256_set1_pd(0.1875f); // 3/16
	                  const __m256d C87964594300514210676954014731826 = 
	                                            = _mm256_set1_pd(87.964594300514210676954014731826);
	                  const __m256d C20          = _mm256_set1_pd(2.0f);
	                  const __m256d C10          = _mm256_set1_pd(1.0f);
	                  const __m256d C30          = _mm256_set1_pd(3.0f);
	                  register __m256d Bg,t,t2,cthti,cthts,sthti,sthts,cphis;
	                  register __m256d AN,Ak0,num,den,alp,x0,x1,x2,x3,x4,secti;
	                  AN    = _mm256_mul_pd(A,N);
	                  t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                  x1    = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  x2    = _mm256_fmadd_pd(epsi,epsi,
	                                               _mm256_mul_pd(x1,x1));
	                  t2    = _mm256_mul_pd(t,t);
	                  secti = _mm256_rcp14_pd(cthti);
	                  sthti = _mm256_sin_pd(thti);
	                  x0    = _mm256_fmadd_pd(t2,C30,C10);
	                  cthts = _mm256_cos_pd(thts);
	                  Ak0   = _mm256_mul_pd(A,k0);
	                  x1    = _mm256_mul_pd(AN,_mm256_mul_pd(Ak0,AK0));
	                  num   = _mm256_mul_pd(x1,x2);
	                  if(pol == 0) {
	                     x1 = _mm256_mul_pd(_mm256_mul_pd(C01875,AN),
	                                                  _mm256_mul_pd(epsi,secti));
	                     x2 = _mm256_mul_pd(x1,x0);
	                     alp= _mm256_mul_pd(k0,x2);
	                  }
	                  else if(pol == 1) {
	                     x1 = mm512_mul_pd(_mm256_mul_pd(C01875,AN),
	                                                  _mm256_mul_pd(epsi,secti));
	                     x2 = _mm256_mul_pd(sthti,sthti);
	                     x3 = _mm256_sub_pd(c10,t2);
	                     alp= _mm256_fmadd_pd(x1,x0,_mm256_mul_pd(x2,x3));
	                     alp= _mm256_mul_pd(k0,alp);
	                  }
	                  x0 = _mm256_div_pd(alp,k0);
	                  x2 = _mm256_add_pd(cthti,cthts);
	                  x1 = _mm256_fmadd_pd(C06,_mm256_mul_pd(x0,x0),
	                                                     _mm256_mul_pd(C30,
	                                                                   _mm256_mul_pd(x2,x2)));
	                  x3 = _mm256_fmadd_pd(sthti,sthti,_mm256_mul_pd(sthts,sthts));
	                  x4 = _mm256_mul_pd(C20,_mm256_mul_pd(sthti,sthts));
	                  x0 = _mm256_mul_pd(x4,cthts);
	                  den= _mm256_mul_pd(C87964594300514210676954014731826,
	                                                             _mm256_add_pd(x1,x0));
	                  
	                  Bg = _mm256_div_pd(num,den);
	                  return (Bg);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d coef_Bg_f9137_ymm4r8_a(const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol) {
	                         
	                  register __m256d A   = _mm256_load_pd(&pA[0]);
	                  register __m256d N   = _mm256_load_pd(&pN[0]); 
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                  register __m256d thti= _mm256_load_pd(&pthti[0]);
	                  register __m256d thts= _mm256_load_pd(&pthts[0]);
	                  register __m256d phis= _mm256_load_pd(&pphis[0]);              
	                  const __m256d C06          = _mm256_set1_pd(0.6f);
	                  const __m256d C01875       = _mm256_set1_pd(0.1875f); // 3/16
	                  const __m256d C87964594300514210676954014731826 = 
	                                            = _mm256_set1_pd(87.964594300514210676954014731826);
	                  const __m256d C20          = _mm256_set1_pd(2.0f);
	                  const __m256d C10          = _mm256_set1_pd(1.0f);
	                  const __m256d C30          = _mm256_set1_pd(3.0f);
	                  register __m256d Bg,t,t2,cthti,cthts,sthti,sthts,cphis;
	                  register __m256d AN,Ak0,num,den,alp,x0,x1,x2,x3,x4,secti;
	                  AN    = _mm256_mul_pd(A,N);
	                  t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                  x1    = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  x2    = _mm256_fmadd_pd(epsi,epsi,
	                                               _mm256_mul_pd(x1,x1));
	                  t2    = _mm256_mul_pd(t,t);
	                  secti = _mm256_rcp14_pd(cthti);
	                  sthti = _mm256_sin_pd(thti);
	                  x0    = _mm256_fmadd_pd(t2,C30,C10);
	                  cthts = _mm256_cos_pd(thts);
	                  Ak0   = _mm256_mul_pd(A,k0);
	                  x1    = _mm256_mul_pd(AN,_mm256_mul_pd(Ak0,AK0));
	                  num   = _mm256_mul_pd(x1,x2);
	                  if(pol == 0) {
	                     x1 = _mm256_mul_pd(_mm256_mul_pd(C01875,AN),
	                                                  _mm256_mul_pd(epsi,secti));
	                     x2 = _mm256_mul_pd(x1,x0);
	                     alp= _mm256_mul_pd(k0,x2);
	                  }
	                  else if(pol == 1) {
	                     x1 = mm512_mul_pd(_mm256_mul_pd(C01875,AN),
	                                                  _mm256_mul_pd(epsi,secti));
	                     x2 = _mm256_mul_pd(sthti,sthti);
	                     x3 = _mm256_sub_pd(c10,t2);
	                     alp= _mm256_fmadd_pd(x1,x0,_mm256_mul_pd(x2,x3));
	                     alp= _mm256_mul_pd(k0,alp);
	                  }
	                  x0 = _mm256_div_pd(alp,k0);
	                  x2 = _mm256_add_pd(cthti,cthts);
	                  x1 = _mm256_fmadd_pd(C06,_mm256_mul_pd(x0,x0),
	                                                     _mm256_mul_pd(C30,
	                                                                   _mm256_mul_pd(x2,x2)));
	                  x3 = _mm256_fmadd_pd(sthti,sthti,_mm256_mul_pd(sthts,sthts));
	                  x4 = _mm256_mul_pd(C20,_mm256_mul_pd(sthti,sthts));
	                  x0 = _mm256_mul_pd(x4,cthts);
	                  den= _mm256_mul_pd(C87964594300514210676954014731826,
	                                                             _mm256_add_pd(x1,x0));
	                  
	                  Bg = _mm256_div_pd(num,den);
	                  return (Bg);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d coef_Bg_f9137_ymm4r8_u(const double * __restrict  pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
	                                          const int pol) {
	                         
	                  register __m256d A   = _mm256_loadu_pd(&pA[0]);
	                  register __m256d N   = _mm256_loadu_pd(&pN[0]); 
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                  register __m256d thti= _mm256_loadu_pd(&pthti[0]);
	                  register __m256d thts= _mm256_loadu_pd(&pthts[0]);
	                  register __m256d phis= _mm256_loadu_pd(&pphis[0]);              
	                  const __m256d C06          = _mm256_set1_pd(0.6f);
	                  const __m256d C01875       = _mm256_set1_pd(0.1875f); // 3/16
	                  const __m256d C87964594300514210676954014731826 = 
	                                            = _mm256_set1_pd(87.964594300514210676954014731826);
	                  const __m256d C20          = _mm256_set1_pd(2.0f);
	                  const __m256d C10          = _mm256_set1_pd(1.0f);
	                  const __m256d C30          = _mm256_set1_pd(3.0f);
	                  register __m256d Bg,t,t2,cthti,cthts,sthti,sthts,cphis;
	                  register __m256d AN,Ak0,num,den,alp,x0,x1,x2,x3,x4,secti;
	                  AN    = _mm256_mul_pd(A,N);
	                  t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                  x1    = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  x2    = _mm256_fmadd_pd(epsi,epsi,
	                                               _mm256_mul_pd(x1,x1));
	                  t2    = _mm256_mul_pd(t,t);
	                  secti = _mm256_rcp14_pd(cthti);
	                  sthti = _mm256_sin_pd(thti);
	                  x0    = _mm256_fmadd_pd(t2,C30,C10);
	                  cthts = _mm256_cos_pd(thts);
	                  Ak0   = _mm256_mul_pd(A,k0);
	                  x1    = _mm256_mul_pd(AN,_mm256_mul_pd(Ak0,AK0));
	                  num   = _mm256_mul_pd(x1,x2);
	                  if(pol == 0) {
	                     x1 = _mm256_mul_pd(_mm256_mul_pd(C01875,AN),
	                                                  _mm256_mul_pd(epsi,secti));
	                     x2 = _mm256_mul_pd(x1,x0);
	                     alp= _mm256_mul_pd(k0,x2);
	                  }
	                  else if(pol == 1) {
	                     x1 = mm512_mul_pd(_mm256_mul_pd(C01875,AN),
	                                                  _mm256_mul_pd(epsi,secti));
	                     x2 = _mm256_mul_pd(sthti,sthti);
	                     x3 = _mm256_sub_pd(c10,t2);
	                     alp= _mm256_fmadd_pd(x1,x0,_mm256_mul_pd(x2,x3));
	                     alp= _mm256_mul_pd(k0,alp);
	                  }
	                  x0 = _mm256_div_pd(alp,k0);
	                  x2 = _mm256_add_pd(cthti,cthts);
	                  x1 = _mm256_fmadd_pd(C06,_mm256_mul_pd(x0,x0),
	                                                     _mm256_mul_pd(C30,
	                                                                   _mm256_mul_pd(x2,x2)));
	                  x3 = _mm256_fmadd_pd(sthti,sthti,_mm256_mul_pd(sthts,sthts));
	                  x4 = _mm256_mul_pd(C20,_mm256_mul_pd(sthti,sthts));
	                  x0 = _mm256_mul_pd(x4,cthts);
	                  den= _mm256_mul_pd(C87964594300514210676954014731826,
	                                                             _mm256_add_pd(x1,x0));
	                  
	                  Bg = _mm256_div_pd(num,den);
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
	           __m256d rcs_hh_f9133_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol) {
	                    
	                  const __m256d C10 = _mm256_set1_pd(1.0f);                    
	                  const __m256d C30 = _mm256_set1_pd(3.0f);
	                  const __m256d C20 = _mm256_set1_pd(2.0f);
	                  const __m256d C80 = _mm256_set1_pd(8.0f);
	                  const __m256d C100= _mm256_set1_pd(10.0f);
	                  const __m256d C240= _mm256_set1_pd(24.0f);
	                  const __m256d C230= _mm256_set1_pd(23.0f);
	                  register rcs,Bg,sphis,x0,t,t2,trm1,trm2,trm3;
	                  t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                  x0    = _mm256_sin_pd(phis);
	                  Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                  t2    = _mm256_mul_pd(t,t);
	                  sphis = _mm256_mul_pd(x0,x0); 
	                  trm1  = _mm256_sub_pd(C30,_mm256_mul_pd(C20,sphis));
	                  trm2  = _mm256_mul_pd(t,_mm256_sub_pd(C80,
	                                                      _mm256_mul_pd(C100,sphis)));
	                  trm3  = _mm256_mul_pd(t2,_mm256_sub_pd(C240,
	                                                      _mm256_mul_pd(C230,sphis)));
	                  x0    = _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                  rcs   = _mm256_mul_pd(Bg,x0);
	                  return (rcs);                             
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9133_ymm4r8_a( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol) {
	                         
	                  register __m256d A   = _mm256_load_pd(&pA[0]);
	                  register __m256d N   = _mm256_load_pd(&pN[0]); 
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                  register __m256d thti= _mm256_load_pd(&pthti[0]);
	                  register __m256d thts= _mm256_load_pd(&pthts[0]);
	                  register __m256d phis= _mm256_load_pd(&pphis[0]); 
	                    
	                  const __m256d C10 = _mm256_set1_pd(1.0f);                    
	                  const __m256d C30 = _mm256_set1_pd(3.0f);
	                  const __m256d C20 = _mm256_set1_pd(2.0f);
	                  const __m256d C80 = _mm256_set1_pd(8.0f);
	                  const __m256d C100= _mm256_set1_pd(10.0f);
	                  const __m256d C240= _mm256_set1_pd(24.0f);
	                  const __m256d C230= _mm256_set1_pd(23.0f);
	                  register rcs,Bg,sphis,x0,t,t2,trm1,trm2,trm3;
	                  t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                  x0    = _mm256_sin_pd(phis);
	                  Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                  t2    = _mm256_mul_pd(t,t);
	                  sphis = _mm256_mul_pd(x0,x0); 
	                  trm1  = _mm256_sub_pd(C30,_mm256_mul_pd(C20,sphis));
	                  trm2  = _mm256_mul_pd(t,_mm256_sub_pd(C80,
	                                                      _mm256_mul_pd(C100,sphis)));
	                  trm3  = _mm256_mul_pd(t2,_mm256_sub_pd(C240,
	                                                      _mm256_mul_pd(C230,sphis)));
	                  x0    = _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                  rcs   = _mm256_mul_pd(Bg,x0);
	                  return (rcs);                             
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9133_ymm4r8_u( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol) {
	                         
	                  register __m256d A   = _mm256_loadu_pd(&pA[0]);
	                  register __m256d N   = _mm256_loadu_pd(&pN[0]); 
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                  register __m256d thti= _mm256_loadu_pd(&pthti[0]);
	                  register __m256d thts= _mm256_loadu_pd(&pthts[0]);
	                  register __m256d phis= _mm256_loadu_pd(&pphis[0]); 
	                    
	                  const __m256d C10 = _mm256_set1_pd(1.0f);                    
	                  const __m256d C30 = _mm256_set1_pd(3.0f);
	                  const __m256d C20 = _mm256_set1_pd(2.0f);
	                  const __m256d C80 = _mm256_set1_pd(8.0f);
	                  const __m256d C100= _mm256_set1_pd(10.0f);
	                  const __m256d C240= _mm256_set1_pd(24.0f);
	                  const __m256d C230= _mm256_set1_pd(23.0f);
	                  register rcs,Bg,sphis,x0,t,t2,trm1,trm2,trm3;
	                  t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                  x0    = _mm256_sin_pd(phis);
	                  Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                  t2    = _mm256_mul_pd(t,t);
	                  sphis = _mm256_mul_pd(x0,x0); 
	                  trm1  = _mm256_sub_pd(C30,_mm256_mul_pd(C20,sphis));
	                  trm2  = _mm256_mul_pd(t,_mm256_sub_pd(C80,
	                                                      _mm256_mul_pd(C100,sphis)));
	                  trm3  = _mm256_mul_pd(t2,_mm256_sub_pd(C240,
	                                                      _mm256_mul_pd(C230,sphis)));
	                  x0    = _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                  rcs   = _mm256_mul_pd(Bg,x0);
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
	           __m256d rcs_vh_f9134_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol) {
	                     
	                    const __m256d C10 = _mm256_set1_pd(1.0f);                             
	                    const __m256d C100= _mm256_set1_pd(10.0f);
	                    const __m256d C20 = _mm256_set1_pd(2.0f);   
	                    const __m256d C230= _mm256_set1_pd(23.0f); 
	                    register __m256d rcs,trm1,trm2,trm3,trm4;
	                    register __m256d Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                    x0    = _mm256_sin_pd(phis);
	                    Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm256_mul_pd(t,t);  
	                    x1    = _mm256_cos_pd(thts);
	                    sphis = _mm256_mul_pd(x0,x0); 
	                    cthts = _mm256_mul_pd(x1,x1);
	                    x0    = _mm256_sin_pd(thts);
	                    sthts = _mm256_mul_pd(x0,x0);
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);  
	                    trm1= _mm256_fmadd_pd(sthts,_mm256_mul_pd(C20,cthts),x1);
	                    x0  = _mm256_sub_pd(gms::math::negate_ymm4r8(C20),
	                                                          _mm256_mul_pd(C40,sthts));
	                    trm2= _mm256_mul_pd(t,
	                                     _mm256_fmadd_pd(sphis,
	                                                    _mm256_mul_pd(C100,cthts),x0));
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);
	                    trm3= _mm256_mul_pd(t2,
	                                      _mm256_fmadd_pd(sphis,
	                                                     _mm256_mul_pd(cthts,C230),x1));
	                    trm4= _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                    rcs = _mm256_mul_pd(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vh_f9134_ymm4r8_a( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol) {
	                         
	                    register __m256d A   = _mm256_load_pd(&pA[0]);
	                    register __m256d N   = _mm256_load_pd(&pN[0]); 
	                    register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                    register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                    register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                    register __m256d thti= _mm256_load_pd(&pthti[0]);
	                    register __m256d thts= _mm256_load_pd(&pthts[0]);
	                    register __m256d phis= _mm256_load_pd(&pphis[0]); 
	                     
	                    const __m256d C10 = _mm256_set1_pd(1.0f);                             
	                    const __m256d C100= _mm256_set1_pd(10.0f);
	                    const __m256d C20 = _mm256_set1_pd(2.0f);   
	                    const __m256d C230= _mm256_set1_pd(23.0f); 
	                    register __m256d rcs,trm1,trm2,trm3,trm4;
	                    register __m256d Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                    x0    = _mm256_sin_pd(phis);
	                    Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm256_mul_pd(t,t);  
	                    x1    = _mm256_cos_pd(thts);
	                    sphis = _mm256_mul_pd(x0,x0); 
	                    cthts = _mm256_mul_pd(x1,x1);
	                    x0    = _mm256_sin_pd(thts);
	                    sthts = _mm256_mul_pd(x0,x0);
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);  
	                    trm1= _mm256_fmadd_pd(sthts,_mm256_mul_pd(C20,cthts),x1);
	                    x0  = _mm256_sub_pd(gms::math::negate_ymm4r8(C20),
	                                                          _mm256_mul_pd(C40,sthts));
	                    trm2= _mm256_mul_pd(t,
	                                     _mm256_fmadd_pd(sphis,
	                                                    _mm256_mul_pd(C100,cthts),x0));
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);
	                    trm3= _mm256_mul_pd(t2,
	                                      _mm256_fmadd_pd(sphis,
	                                                     _mm256_mul_pd(cthts,C230),x1));
	                    trm4= _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                    rcs = _mm256_mul_pd(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vh_f9134_ymm4r8_u( const double * __restrict  pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
	                                          const int pol) {
	                         
	                    register __m256d A   = _mm256_loadu_pd(&pA[0]);
	                    register __m256d N   = _mm256_loadu_pd(&pN[0]); 
	                    register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                    register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                    register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                    register __m256d thti= _mm256_loadu_pd(&pthti[0]);
	                    register __m256d thts= _mm256_loadu_pd(&pthts[0]);
	                    register __m256d phis= _mm256_loadu_pd(&pphis[0]); 
	                     
	                    const __m256d C10 = _mm256_set1_pd(1.0f);                             
	                    const __m256d C100= _mm256_set1_pd(10.0f);
	                    const __m256d C20 = _mm256_set1_pd(2.0f);   
	                    const __m256d C230= _mm256_set1_pd(23.0f); 
	                    register __m256d rcs,trm1,trm2,trm3,trm4;
	                    register __m256d Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                    x0    = _mm256_sin_pd(phis);
	                    Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm256_mul_pd(t,t);  
	                    x1    = _mm256_cos_pd(thts);
	                    sphis = _mm256_mul_pd(x0,x0); 
	                    cthts = _mm256_mul_pd(x1,x1);
	                    x0    = _mm256_sin_pd(thts);
	                    sthts = _mm256_mul_pd(x0,x0);
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);  
	                    trm1= _mm256_fmadd_pd(sthts,_mm256_mul_pd(C20,cthts),x1);
	                    x0  = _mm256_sub_pd(gms::math::negate_ymm4r8(C20),
	                                                          _mm256_mul_pd(C40,sthts));
	                    trm2= _mm256_mul_pd(t,
	                                     _mm256_fmadd_pd(sphis,
	                                                    _mm256_mul_pd(C100,cthts),x0));
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);
	                    trm3= _mm256_mul_pd(t2,
	                                      _mm256_fmadd_pd(sphis,
	                                                     _mm256_mul_pd(cthts,C230),x1));
	                    trm4= _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                    rcs = _mm256_mul_pd(Bg,trm4);
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
	           __m256d rcs_hv_f9135_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol) {
	                     
	                    const __m256d C10 = _mm256_set1_pd(1.0f);                             
	                    const __m256d C100= _mm256_set1_pd(10.0f);
	                    const __m256d C20 = _mm256_set1_pd(2.0f);   
	                    const __m256d C230= _mm256_set1_pd(23.0f); 
	                    register __m256d rcs,trm1,trm2,trm3,trm4;
	                    register __m256d Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                    x0    = _mm256_sin_pd(phis);
	                    Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm256_mul_pd(t,t);  
	                    x1    = _mm256_cos_pd(thti);
	                    sphis = _mm256_mul_pd(x0,x0); 
	                    cthts = _mm256_mul_pd(x1,x1);
	                    x0    = _mm256_sin_pd(thti);
	                    sthts = _mm256_mul_pd(x0,x0);
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);  
	                    trm1= _mm256_fmadd_pd(sthts,_mm256_mul_pd(C20,cthts),x1);
	                    x0  = _mm256_sub_pd(gms::math::negate_ymm4r8(C20),
	                                                          _mm256_mul_pd(C40,sthts));
	                    trm2= _mm256_mul_pd(t,
	                                     _mm256_fmadd_pd(sphis,
	                                                    _mm256_mul_pd(C100,cthts),x0));
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);
	                    trm3= _mm256_mul_pd(t2,
	                                      _mm256_fmadd_pd(sphis,
	                                                     _mm256_mul_pd(cthts,C230),x1));
	                    trm4= _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                    rcs = _mm256_mul_pd(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hv_f9135_ymm4r8_a(  const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol) {
	                         
	                    register __m256d A   = _mm256_load_pd(&pA[0]);
	                    register __m256d N   = _mm256_load_pd(&pN[0]); 
	                    register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                    register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                    register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                    register __m256d thti= _mm256_load_pd(&pthti[0]);
	                    register __m256d thts= _mm256_load_pd(&pthts[0]);
	                    register __m256d phis= _mm256_load_pd(&pphis[0]); 
	                     
	                    const __m256d C10 = _mm256_set1_pd(1.0f);                             
	                    const __m256d C100= _mm256_set1_pd(10.0f);
	                    const __m256d C20 = _mm256_set1_pd(2.0f);   
	                    const __m256d C230= _mm256_set1_pd(23.0f); 
	                    register __m256d rcs,trm1,trm2,trm3,trm4;
	                    register __m256d Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                    x0    = _mm256_sin_pd(phis);
	                    Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm256_mul_pd(t,t);  
	                    x1    = _mm256_cos_pd(thti);
	                    sphis = _mm256_mul_pd(x0,x0); 
	                    cthts = _mm256_mul_pd(x1,x1);
	                    x0    = _mm256_sin_pd(thti);
	                    sthts = _mm256_mul_pd(x0,x0);
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);  
	                    trm1= _mm256_fmadd_pd(sthts,_mm256_mul_pd(C20,cthts),x1);
	                    x0  = _mm256_sub_pd(gms::math::negate_ymm4r8(C20),
	                                                          _mm256_mul_pd(C40,sthts));
	                    trm2= _mm256_mul_pd(t,
	                                     _mm256_fmadd_pd(sphis,
	                                                    _mm256_mul_pd(C100,cthts),x0));
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);
	                    trm3= _mm256_mul_pd(t2,
	                                      _mm256_fmadd_pd(sphis,
	                                                     _mm256_mul_pd(cthts,C230),x1));
	                    trm4= _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                    rcs = _mm256_mul_pd(Bg,trm4);
	                    return (rcs);                
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d  rcs_hv_f9135_ymm4r8_u(const double * __restrict pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
	                                          const int pol) {
	                         
	                    register __m256d A   = _mm256_loadu_pd(&pA[0]);
	                    register __m256d N   = _mm256_loadu_pd(&pN[0]); 
	                    register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                    register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                    register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                    register __m256d thti= _mm256_loadu_pd(&pthti[0]);
	                    register __m256d thts= _mm256_loadu_pd(&pthts[0]);
	                    register __m256d phis= _mm256_loadu_pd(&pphis[0]); 
	                     
	                    const __m256d C10 = _mm256_set1_pd(1.0f);                             
	                    const __m256d C100= _mm256_set1_pd(10.0f);
	                    const __m256d C20 = _mm256_set1_pd(2.0f);   
	                    const __m256d C230= _mm256_set1_pd(23.0f); 
	                    register __m256d rcs,trm1,trm2,trm3,trm4;
	                    register __m256d Bg,sthts,cthts,sphis,x0,x1,t,t2;
	                    t     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));
	                    x0    = _mm256_sin_pd(phis);
	                    Bg    = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol);
	                    t2    = _mm256_mul_pd(t,t);  
	                    x1    = _mm256_cos_pd(thti);
	                    sphis = _mm256_mul_pd(x0,x0); 
	                    cthts = _mm256_mul_pd(x1,x1);
	                    x0    = _mm256_sin_pd(thti);
	                    sthts = _mm256_mul_pd(x0,x0);
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);  
	                    trm1= _mm256_fmadd_pd(sthts,_mm256_mul_pd(C20,cthts),x1);
	                    x0  = _mm256_sub_pd(gms::math::negate_ymm4r8(C20),
	                                                          _mm256_mul_pd(C40,sthts));
	                    trm2= _mm256_mul_pd(t,
	                                     _mm256_fmadd_pd(sphis,
	                                                    _mm256_mul_pd(C100,cthts),x0));
	                    x1  = _mm256_fmadd_pd(sthts,C20,C10);
	                    trm3= _mm256_mul_pd(t2,
	                                      _mm256_fmadd_pd(sphis,
	                                                     _mm256_mul_pd(cthts,C230),x1));
	                    trm4= _mm256_add_pd(trm1,_mm256_add_pd(trm2,trm3));
	                    rcs = _mm256_mul_pd(Bg,trm4);
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
	           __m256d rcs_vv_f9136_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol) {
	                                        
	                const __m256d C10 = _mm256_set1_pd(1.0f);     
	                const __m256d C20 = _mm256_set1_pd(2.0f);   
	                const __m256d C30 = _mm256_set1_pd(3.0f);
	                const __m256d C120= _mm256_set1_pd(12.0f);
	                const __m256d C140= _mm256_set1_pd(14.0f);
	                const __m256d C350= _mm256_set1_pd(35.0f);
	                register __m256d t,t2,x0,x1,x2,x3,cthti,cthts;
	                register __m256d sphis,cphis,sthts,cthts;
	                register __m256d trm1,trm2,trm3,trm4;
	                register __m256d rcs,Bg,cterm,sterm,sctrm;
	                x0     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));   
	                t      = _mm256_sub_pd(C10,x0);
	                t2     = _mm256_sub_pd(C10,_mm256_mul_pd(x0,x0));
	                Bg     = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol); 
	                sphis  = _mm256_sin_pd(phis);
	                cphis  = _mm256_cos_pd(phis);
	                cthti  = _mm256_cos_pd(thti);
	                sthti  = _mm256_sin_pd(thti);
	                cthts  = _mm256_cos_pd(thts);
	                cterm  = _mm256_mul_pd(cthti,_mm256_mul_pd(cthts,cphis));
	                sthts  = _mm256_sin_pd(thts);
	                sterm  = _mm256_mul_pd(sthti,sthts);
	                sctrm  = _mm256_sub_pd(sterm,cterm);
	                x0     = _mm256_mul_pd(sphis,sphis);
	                x1     = _mm256_mul_pd(cthts,cthts);
	                x2     = _mm256_mul_pd(C30,_mm256_mul_pd(cthti,cthti));
	                x3     = _mm256_sub_pd(C30,_mm256_mul_pd(x0,
	                                                    _mm256_mul_pd(x1,x2)));
	                trm1   = _mm256_fmadd_pd(t2,x3,sctrm);
	                trm2   = _mm256_mul_pd(C120,_mm256_mul_pd(t2,sterm));
	                x0     = _mm256_fmsub_pd(C30,sterm,cterm);
	                trm3   = _mm256_mul_pd(C140,_mm256_mul_pd(t,x0));
	                trm4   = _mm256_mul_pd(_mm256_mul_pd(C350,t2),sctrm);
	                x1     = _mm256_mul_pd(trm1,trm2);
	                x2     = _mm256_add_pd(x1,_mm256_add_pd(trm3,trm4));
	                rcs    = _mm256_mul_pd(Bg,x2);
	                return (rcs);
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9136_ymm4r8_a( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol ) {
	                                      
	                register __m256d A   = _mm256_load_pd(&pA[0]);
	                register __m256d N   = _mm256_load_pd(&pN[0]); 
	                register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                register __m256d thti= _mm256_load_pd(&pthti[0]);
	                register __m256d thts= _mm256_load_pd(&pthts[0]);
	                register __m256d phis= _mm256_load_pd(&pphis[0]);   
	                const __m256d C10 = _mm256_set1_pd(1.0f);     
	                const __m256d C20 = _mm256_set1_pd(2.0f);   
	                const __m256d C30 = _mm256_set1_pd(3.0f);
	                const __m256d C120= _mm256_set1_pd(12.0f);
	                const __m256d C140= _mm256_set1_pd(14.0f);
	                const __m256d C350= _mm256_set1_pd(35.0f);
	                register __m256d t,t2,x0,x1,x2,x3,cthti,cthts;
	                register __m256d sphis,cphis,sthts,cthts;
	                register __m256d trm1,trm2,trm3,trm4;
	                register __m256d rcs,Bg,cterm,sterm,sctrm;
	                x0     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));   
	                t      = _mm256_sub_pd(C10,x0);
	                t2     = _mm256_sub_pd(C10,_mm256_mul_pd(x0,x0));
	                Bg     = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol); 
	                sphis  = _mm256_sin_pd(phis);
	                cphis  = _mm256_cos_pd(phis);
	                cthti  = _mm256_cos_pd(thti);
	                sthti  = _mm256_sin_pd(thti);
	                cthts  = _mm256_cos_pd(thts);
	                cterm  = _mm256_mul_pd(cthti,_mm256_mul_pd(cthts,cphis));
	                sthts  = _mm256_sin_pd(thts);
	                sterm  = _mm256_mul_pd(sthti,sthts);
	                sctrm  = _mm256_sub_pd(sterm,cterm);
	                x0     = _mm256_mul_pd(sphis,sphis);
	                x1     = _mm256_mul_pd(cthts,cthts);
	                x2     = _mm256_mul_pd(C30,_mm256_mul_pd(cthti,cthti));
	                x3     = _mm256_sub_pd(C30,_mm256_mul_pd(x0,
	                                                    _mm256_mul_pd(x1,x2)));
	                trm1   = _mm256_fmadd_pd(t2,x3,sctrm);
	                trm2   = _mm256_mul_pd(C120,_mm256_mul_pd(t2,sterm));
	                x0     = _mm256_fmsub_pd(C30,sterm,cterm);
	                trm3   = _mm256_mul_pd(C140,_mm256_mul_pd(t,x0));
	                trm4   = _mm256_mul_pd(_mm256_mul_pd(C350,t2),sctrm);
	                x1     = _mm256_mul_pd(trm1,trm2);
	                x2     = _mm256_add_pd(x1,_mm256_add_pd(trm3,trm4));
	                rcs    = _mm256_mul_pd(Bg,x2);
	                return (rcs);
	         }
	         
	         
	         
	              __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9136_ymm4r8_u( const double * __restrict  pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
	                                          const int pol ) {
	                                      
	                register __m256d A   = _mm256_loadu_pd(&pA[0]);
	                register __m256d N   = _mm256_loadu_pd(&pN[0]); 
	                register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                register __m256d thti= _mm256_loadu_pd(&pthti[0]);
	                register __m256d thts= _mm256_loadu_pd(&pthts[0]);
	                register __m256d phis= _mm256_loadu_pd(&pphis[0]);   
	                const __m256d C10 = _mm256_set1_pd(1.0f);     
	                const __m256d C20 = _mm256_set1_pd(2.0f);   
	                const __m256d C30 = _mm256_set1_pd(3.0f);
	                const __m256d C120= _mm256_set1_pd(12.0f);
	                const __m256d C140= _mm256_set1_pd(14.0f);
	                const __m256d C350= _mm256_set1_pd(35.0f);
	                register __m256d t,t2,x0,x1,x2,x3,cthti,cthts;
	                register __m256d sphis,cphis,sthts,cthts;
	                register __m256d trm1,trm2,trm3,trm4;
	                register __m256d rcs,Bg,cterm,sterm,sctrm;
	                x0     = _mm256_div_pd(C20,_mm256_add_pd(C10,epsr));   
	                t      = _mm256_sub_pd(C10,x0);
	                t2     = _mm256_sub_pd(C10,_mm256_mul_pd(x0,x0));
	                Bg     = coef_Bg_f9137_ymm4r8(A,N,k0,epsr,epsi,thti,thts,phis,pol); 
	                sphis  = _mm256_sin_pd(phis);
	                cphis  = _mm256_cos_pd(phis);
	                cthti  = _mm256_cos_pd(thti);
	                sthti  = _mm256_sin_pd(thti);
	                cthts  = _mm256_cos_pd(thts);
	                cterm  = _mm256_mul_pd(cthti,_mm256_mul_pd(cthts,cphis));
	                sthts  = _mm256_sin_pd(thts);
	                sterm  = _mm256_mul_pd(sthti,sthts);
	                sctrm  = _mm256_sub_pd(sterm,cterm);
	                x0     = _mm256_mul_pd(sphis,sphis);
	                x1     = _mm256_mul_pd(cthts,cthts);
	                x2     = _mm256_mul_pd(C30,_mm256_mul_pd(cthti,cthti));
	                x3     = _mm256_sub_pd(C30,_mm256_mul_pd(x0,
	                                                    _mm256_mul_pd(x1,x2)));
	                trm1   = _mm256_fmadd_pd(t2,x3,sctrm);
	                trm2   = _mm256_mul_pd(C120,_mm256_mul_pd(t2,sterm));
	                x0     = _mm256_fmsub_pd(C30,sterm,cterm);
	                trm3   = _mm256_mul_pd(C140,_mm256_mul_pd(t,x0));
	                trm4   = _mm256_mul_pd(_mm256_mul_pd(C350,t2),sctrm);
	                x1     = _mm256_mul_pd(trm1,trm2);
	                x2     = _mm256_add_pd(x1,_mm256_add_pd(trm3,trm4));
	                rcs    = _mm256_mul_pd(Bg,x2);
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
	           __m256d rcs_vv_f9174_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui) {
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  register __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d epsrm1,epsim1;
	                  register __m256d epsr2,epsi2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsi,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm256_mul_pd(k02,_mm256_mul_pd(l2,sthti));   
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C40,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm256_fmadd_pd(epsrm1,sthti,t1r);
	                  epsim1= _mm256_fmadd_pd(epsim1,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(epsr,cthti);
	                  x1    = _mm256_mul_pd(epsi,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm256_sub_pd(epsrm1,x2);
	                  epsim1 = _mm256_sub_pd(epsim1,x3);
	                  cdiv_ymm4c8(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9174_ymm4r8_a(const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui) {
	                        
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d h   = _mm256_load_pd(&ph[0]);
	                  register __m256d l   = _mm256_load_pd(&pl[0]);
	                  register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_load_pd(&pmur[0]);
	                  register __m256d mui = _mm256_load_pd(&pmui[0]);               
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                   __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d epsrm1,epsim1;
	                  register __m256d epsr2,epsi2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsi,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm256_mul_pd(k02,_mm256_mul_pd(l2,sthti));   
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C40,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm256_fmadd_pd(epsrm1,sthti,t1r);
	                  epsim1= _mm256_fmadd_pd(epsim1,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(epsr,cthti);
	                  x1    = _mm256_mul_pd(epsi,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm256_sub_pd(epsrm1,x2);
	                  epsim1 = _mm256_sub_pd(epsim1,x3);
	                  cdiv_ymm4c8(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9174_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui) {
	                        
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d h   = _mm256_loadu_pd(&ph[0]);
	                  register __m256d l   = _mm256_loadu_pd(&pl[0]);
	                  register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_loadu_pd(&pmur[0]);
	                  register __m256d mui = _mm256_loadu_pd(&pmui[0]);               
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                   __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d epsrm1,epsim1;
	                  register __m256d epsr2,epsi2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsi,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm256_mul_pd(k02,_mm256_mul_pd(l2,sthti));   
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C40,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm256_fmadd_pd(epsrm1,sthti,t1r);
	                  epsim1= _mm256_fmadd_pd(epsim1,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(epsr,cthti);
	                  x1    = _mm256_mul_pd(epsi,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm256_sub_pd(epsrm1,x2);
	                  epsim1 = _mm256_sub_pd(epsim1,x3);
	                  cdiv_ymm4c8(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
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
	           __m256d rcs_hh_f9175_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui) {
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  register __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d murm1,muim1,murm1s,muim1s;
	                  register __m256d mur2,mui2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsi,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsr,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm256_mul_pd(k02,_mm256_mul_pd(l2,sthti));   
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C40,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm256_fmadd_pd(murm1s,sthti,t1r);
	                  muim1s= _mm256_fmadd_pd(muim1s,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(mur,cthti);
	                  x1    = _mm256_mul_pd(mui,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm256_sub_pd(murm1,x2);
	                  muim1 = _mm256_sub_pd(muim1,x3);
	                  cdiv_ymm4c8(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9175_ymm4r8_a(  const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui) {
	                        
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d h   = _mm256_load_pd(&ph[0]);
	                  register __m256d l   = _mm256_load_pd(&pl[0]);
	                  register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_load_pd(&pmur[0]);
	                  register __m256d mui = _mm256_load_pd(&pmui[0]);  
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d murm1,muim1,murm1s,muim1s;
	                  register __m256d mur2,mui2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsi,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsr,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm256_mul_pd(k02,_mm256_mul_pd(l2,sthti));   
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C40,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm256_fmadd_pd(murm1s,sthti,t1r);
	                  muim1s= _mm256_fmadd_pd(muim1s,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(mur,cthti);
	                  x1    = _mm256_mul_pd(mui,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm256_sub_pd(murm1,x2);
	                  muim1 = _mm256_sub_pd(muim1,x3);
	                  cdiv_ymm4c8(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9175_ymm4r8_u(  const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui) {
	                        
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d h   = _mm256_loadu_pd(&ph[0]);
	                  register __m256d l   = _mm256_loadu_pd(&pl[0]);
	                  register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_loadu_pd(&pmur[0]);
	                  register __m256d mui = _mm256_loadu_pd(&pmui[0]);  
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d murm1,muim1,murm1s,muim1s;
	                  register __m256d mur2,mui2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsi,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsr,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm256_mul_pd(k02,_mm256_mul_pd(l2,sthti));   
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C40,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm256_fmadd_pd(murm1s,sthti,t1r);
	                  muim1s= _mm256_fmadd_pd(muim1s,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(mur,cthti);
	                  x1    = _mm256_mul_pd(mui,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm256_sub_pd(murm1,x2);
	                  muim1 = _mm256_sub_pd(muim1,x3);
	                  cdiv_ymm4c8(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
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
	           __m256d rcs_hv_f9176_ymm4r8() { 
	           
	                return _mm256_setzero_pd();
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
	           __m256d rcs_vv_f9177_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui) {
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  const __m256d C80  = _mm256_set1_pd(8.0f);
	                  const __m256d C150 = _mm256_set1_pd(1.5f);
	                  register __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d epsrm1,epsim1;
	                  register __m256d epsr2,epsi2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsi,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm256_fmadd_pd(_mm256_mul_pd(sthti,l2),
	                                          _mm256_mul_pd(k0,C40),C10);
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  x3    = _mm256_mul_pd(arg,_mm256_sqrt_pd(arg));
	                  inve  = _mm256_rcp14_pd(x3);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C80,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm256_fmadd_pd(epsrm1,sthti,t1r);
	                  epsim1= _mm256_fmadd_pd(epsim1,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(epsr,cthti);
	                  x1    = _mm256_mul_pd(epsi,cthti);
	                  cmul_ymm4c8(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm256_sub_pd(epsrm1,x2);
	                  epsim1 = _mm256_sub_pd(epsim1,x3);
	                  cdiv_ymm4c8(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9177_ymm4r8_a(const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui) {
	                        
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d h   = _mm256_load_pd(&ph[0]);
	                  register __m256d l   = _mm256_load_pd(&pl[0]);
	                  register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_load_pd(&pmur[0]);
	                  register __m256d mui = _mm256_load_pd(&pmui[0]); 
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  const __m256d C80  = _mm256_set1_pd(8.0f);
	                  const __m256d C150 = _mm256_set1_pd(1.5f);
	                  __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d epsrm1,epsim1;
	                  register __m256d epsr2,epsi2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsi,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm256_fmadd_pd(_mm256_mul_pd(sthti,l2),
	                                          _mm256_mul_pd(k0,C40),C10);
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  x3    = _mm256_mul_pd(arg,_mm256_sqrt_pd(arg));
	                  inve  = _mm256_rcp14_pd(x3);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C80,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm256_fmadd_pd(epsrm1,sthti,t1r);
	                  epsim1= _mm256_fmadd_pd(epsim1,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(epsr,cthti);
	                  x1    = _mm256_mul_pd(epsi,cthti);
	                  cmul_ymm4c8(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm256_sub_pd(epsrm1,x2);
	                  epsim1 = _mm256_sub_pd(epsim1,x3);
	                  cdiv_ymm4c8(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9177_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui) {
	                        
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d h   = _mm256_loadu_pd(&ph[0]);
	                  register __m256d l   = _mm256_loadu_pd(&pl[0]);
	                  register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_loadu_pd(&pmur[0]);
	                  register __m256d mui = _mm256_loadu_pd(&pmui[0]); 
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  const __m256d C80  = _mm256_set1_pd(8.0f);
	                  const __m256d C150 = _mm256_set1_pd(1.5f);
	                  __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d epsrm1,epsim1;
	                  register __m256d epsr2,epsi2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,murm1,muim1,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsr,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsi,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(t0r,t0i,t0r,t0i,epsrm1,epsim1);
	                  arg   = _mm256_fmadd_pd(_mm256_mul_pd(sthti,l2),
	                                          _mm256_mul_pd(k0,C40),C10);
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  x3    = _mm256_mul_pd(arg,_mm256_sqrt_pd(arg));
	                  inve  = _mm256_rcp14_pd(x3);
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C80,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  epsrm1= _mm256_fmadd_pd(epsrm1,sthti,t1r);
	                  epsim1= _mm256_fmadd_pd(epsim1,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(epsr,cthti);
	                  x1    = _mm256_mul_pd(epsi,cthti);
	                  cmul_ymm4c8(epsr,epsi,epsr,epsi,&epsr2,&epsi2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(epsr2,epsi2,murm1,muim1,&x2,&x3);
	                  epsrm1 = _mm256_sub_pd(epsrm1,x2);
	                  epsim1 = _mm256_sub_pd(epsim1,x3);
	                  cdiv_ymm4c8(epsrm1,epsim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
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
	           __m256d rcs_hh_f9178_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui) {
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  const __m256d C80  = _mm256_set1_pd(8.0f);
	                  register __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d murm1,muim1,murm1s,muim1s;
	                  register __m256d mur2,mui2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsi,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsr,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  x3    = _mm256_mul_pd(C40,x0);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm256_fmadd_pd(_mm256_mul_pd(sthti,l2),x3,C10);
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_mul_pd(arg,_mm256_sqrt_pd(arg));
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C80,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm256_fmadd_pd(murm1s,sthti,t1r);
	                  muim1s= _mm256_fmadd_pd(muim1s,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(mur,cthti);
	                  x1    = _mm256_mul_pd(mui,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm256_sub_pd(murm1,x2);
	                  muim1 = _mm256_sub_pd(muim1,x3);
	                  cdiv_ymm4c8(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9178_ymm4r8_a(const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui) {
	                        
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d h   = _mm256_load_pd(&ph[0]);
	                  register __m256d l   = _mm256_load_pd(&pl[0]);
	                  register __m256d epsr= _mm256_load_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_load_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_load_pd(&pmur[0]);
	                  register __m256d mui = _mm256_load_pd(&pmui[0]); 
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  const __m256d C80  = _mm256_set1_pd(8.0f);
	                  register __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d murm1,muim1,murm1s,muim1s;
	                  register __m256d mur2,mui2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsi,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsr,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  x3    = _mm256_mul_pd(C40,x0);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm256_fmadd_pd(_mm256_mul_pd(sthti,l2),x3,C10);
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_mul_pd(arg,_mm256_sqrt_pd(arg));
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C80,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm256_fmadd_pd(murm1s,sthti,t1r);
	                  muim1s= _mm256_fmadd_pd(muim1s,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(mur,cthti);
	                  x1    = _mm256_mul_pd(mui,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm256_sub_pd(murm1,x2);
	                  muim1 = _mm256_sub_pd(muim1,x3);
	                  cdiv_ymm4c8(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
	                  return (rcs);
	       }
	       
	       
	       
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9178_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui) {
	                        
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d h   = _mm256_loadu_pd(&ph[0]);
	                  register __m256d l   = _mm256_loadu_pd(&pl[0]);
	                  register __m256d epsr= _mm256_loadu_pd(&pepsr[0]);
	                  register __m256d epsi= _mm256_loadu_pd(&pepsi[0]);
	                  register __m256d mur = _mm256_loadu_pd(&pmur[0]);
	                  register __m256d mui = _mm256_loadu_pd(&pmui[0]); 
	                                       
	                  const __m256d C10  = _mm256_set1_pd(1.0f);
	                  const __m256d C40  = _mm256_set1_pd(4.0f);
	                  const __m256d C80  = _mm256_set1_pd(8.0f);
	                  register __m256d k04,h2,l2,x0,x1,x2,x3,cthti,sthti;
	                  register __m256d murm1,muim1,murm1s,muim1s;
	                  register __m256d mur2,mui2,t0r,t0i,t1r,t1i;
	                  register __m256d t2r,t2i,inve,t3r,t3i;
	                  register __m256d rcs,cabs,arg,earg,frac;
	                  h2    = _mm256_mul_pd(h,h);
	                  t0r   = _mm256_sub_pd(epsi,C10);
	                  cthti = _mm256_cos_pd(thti);
	                  t0i   = _mm256_sub_pd(epsr,C10);
	                  murm1 = _mm256_sub_pd(mur,C10);
	                  x0    = _mm256_mul_pd(k0,k0);
	                  muim1 = _mm256_mul_pd(mui,C10);
	                  x3    = _mm256_mul_pd(C40,x0);
	                  l2    = _mm256_mul_pd(l,l);
	                  sthti = _mm256_sin_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  sthti = _mm256_mul_pd(sthti,sthti);
	                  cmul_ymm4c8(murm1,muim1,murm1,muim1,&murm1s,&muim1s);
	                  arg   = _mm256_fmadd_pd(_mm256_mul_pd(sthti,l2),x3,C10);
	                  x0    = _mm256_mul_pd(cthti,cthti);
	                  earg  = _mm256_mul_pd(arg,_mm256_sqrt_pd(arg));
	                  x1    = _mm256_mul_pd(x0,x0);
	                  x2    = _mm256_mul_pd(C80,k04);
	                  frac  = _mm256_mul_pd(_mm256_mul_pd(x2,h2),
	                                        _mm256_mul_pd(l2,x1));
	                  cmul_ymm4c8(epsr,epsi,mur,mui,&t1r,&t1i);
	                  murm1s= _mm256_fmadd_pd(murm1s,sthti,t1r);
	                  muim1s= _mm256_fmadd_pd(muim1s,sthti,t1i);
	                  t1r   = _mm256_sub_pd(t1r,sthti);
	                  t1i   = _mm256_sub_pd(t1i,sthti);
	                  csqrt_ymm4c8(t1r,t1i,&t2r,&t2i);
	                  x0    = _mm256_mul_pd(mur,cthti);
	                  x1    = _mm256_mul_pd(mui,cthti);
	                  inve  = _mm256_rcp14_pd(earg);
	                  cmul_ymm4c8(mur,mui,mur,mui,&mur2,&mui2);
	                  x0    = _mm256_add_pd(x0,t2r);
	                  x1    = _mm256_add_pd(x1,t2i);
	                  cmul_ymm4c8(x0,x1,x0,x1,&t3r,&t3i);//denom
	                  cmul_ymm4c8(mur2,mui2,t0r,t0i,&x2,&x3);
	                  murm1 = _mm256_sub_pd(murm1,x2);
	                  muim1 = _mm256_sub_pd(muim1,x3);
	                  cdiv_ymm4c8(murm1,muim1,t3r,t3i,&t0r,&t0i); // ratio (complex).
	                  cabs = cabs_ymm4c8(t0r,t0i);
	                  rcs  = _mm256_mul_pd(frac,_mm256_mul_pd(cabs,inve));
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
	           __m256d rcs_hv_f9179_ymm4r8() { 
	           
	                return _mm256_setzero_pd();
	         } 
	         
	         
	         /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-80
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vv_f9180_ymm4r8(const __m256d thti,
	                                         const __m256d thts,
	                                         const __m256d phis) {
	                                         
	                  register __m256d sthti,sthts,cthti,cthts,cphis;
	                  register __m256d num,den,avv;  
	                  sthti = _mm256_sin_pd(thti);
	                  cthti = _mm256_cos_pd(thti);
	                  cphis = _mm256_cos_pd(phis);
	                  cthts = _mm256_cos_pd(thts);
	                  sthts = _mm256_sin_pd(thts);
	                  num   = _mm256_fmsub_pd(sthti,sthts,cphis);
	                  den   = _mm256_mul_pd(cthti,cthts);
	                  avv   = _mm256_div_pd(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vv_f9180_ymm4r8_a(  const double * __restrict pthti,
	                                         const double * __restrict pthts,
	                                         const double * __restrict pphis) {
	                     
	                  register __m256d thti = _mm256_load_pd(&pthti[0]);
	                  register __m256d thts = _mm256_load_pd(&pthts[0]);
	                  register __m256d phis = _mm256_load_pd(&phis[0]);                    
	                  register __m256d sthti,sthts,cthti,cthts,cphis;
	                  register __m256d num,den,avv;  
	                  sthti = _mm256_sin_pd(thti);
	                  cthti = _mm256_cos_pd(thti);
	                  cphis = _mm256_cos_pd(phis);
	                  cthts = _mm256_cos_pd(thts);
	                  sthts = _mm256_sin_pd(thts);
	                  num   = _mm256_fmsub_pd(sthti,sthts,cphis);
	                  den   = _mm256_mul_pd(cthti,cthts);
	                  avv   = _mm256_div_pd(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vv_f9180_ymm4r8_u(  const double * __restrict  pthti,
	                                         const double * __restrict  pthts,
	                                         const double * __restrict  pphis) {
	                     
	                  register __m256d thti = _mm256_loadu_pd(&pthti[0]);
	                  register __m256d thts = _mm256_loadu_pd(&pthts[0]);
	                  register __m256d phis = _mm256_loadu_pd(&phis[0]);                    
	                  register __m256d sthti,sthts,cthti,cthts,cphis;
	                  register __m256d num,den,avv;  
	                  sthti = _mm256_sin_pd(thti);
	                  cthti = _mm256_cos_pd(thti);
	                  cphis = _mm256_cos_pd(phis);
	                  cthts = _mm256_cos_pd(thts);
	                  sthts = _mm256_sin_pd(thts);
	                  num   = _mm256_fmsub_pd(sthti,sthts,cphis);
	                  den   = _mm256_mul_pd(cthti,cthts);
	                  avv   = _mm256_div_pd(num,den);
	                  return (avv);        
	          }
	          
	          
	          /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-81
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_hv_f9181_ymm4r8(const __m256d phis,
	                                     const __m256d thti) {
	                                     
	                  register __m256d sphis,cthti;
	                  register __m256d ahv;
	                  sphis = _mm256_sin_pd(phis);
	                  cthti = _mm256_cos_pd(thti);
	                  ahv   = _mm256_div_pd(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_hv_f9181_ymm4r8_a(const double * __restrict pphis,
	                                     const double * __restrict pthti) {
	                                     
	                  register __m256d thti = _mm256_load_pd(&pthti[0]);
	                  register __m256d phis = _mm256_load_pd(&phis[0]); 
	                  register __m256d sphis,cthti;
	                  register __m256d ahv;
	                  sphis = _mm256_sin_pd(phis);
	                  cthti = _mm256_cos_pd(thti);
	                  ahv   = _mm256_div_pd(sphis,cthti);
	                  return (ahv);                           
	         }
	          
	       
	        
	           __ATTR_ALWAYS_INLINE__
	           
	         
                   
	           static inline
	           __m256d a_hv_f9181_ymm4r8_u(const double * __restrict pphis,
	                                     const double * __restrict  pthti) {
	                                     
	                  register __m256d thti = _mm256_loadu_pd(&pthti[0]);
	                  register __m256d phis = _mm256_loadu_pd(&phis[0]); 
	                  register __m256d sphis,cthti;
	                  register __m256d ahv;
	                  sphis = _mm256_sin_pd(phis);
	                  cthti = _mm256_cos_pd(thti);
	                  ahv   = _mm256_div_pd(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	            /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-82
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_hh_f9182_ymm4r8(const __m256d phis) {
	                 
	                  register __m256d cphis,ahh;
	                  cphis = _mm256_cos_pd(phis);
	                  ahh   = gms::math::negate_ymm4r8(cphis);
	                  return (ahh);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_hh_f9182_ymm4r8_a(const double * __restrict pphis) {
	                 
	                  register __m256d phis = _mm256_load_pd(&pphis[0]);
	                  register __m256d cphis,ahh;
	                  cphis = _mm256_cos_pd(phis);
	                  ahh   = gms::math::negate_ymm4r8(cphis);
	                  return (ahh);
	          }
	        
	       
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_hh_f9182_ymm4r8_u(const double * __restrict  pphis) {
	                 
	                  register __m256d phis = _mm256_loadu_pd(&pphis[0]);
	                  register __m256d cphis,ahh;
	                  cphis = _mm256_cos_pd(phis);
	                  ahh   = gms::math::negate_ymm4r8(cphis);
	                  return (ahh);
	          }
	          
	          
	             /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-83
	         */
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vh_f9183_ymm4r8(const __m256d phis,
	                                     const __m256d thts) {
	                                     
	                  register __m256d sphis,cthts;
	                  register __m256d ahv;
	                  sphis = _mm256_sin_pd(phis);
	                  cthts = _mm256_cos_pd(thti);
	                  ahv   = gms::math::negate_ymm4r8(_mm256_div_pd(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vh_f9183_ymm4r8_a(const double * __restrict pphis,
	                                       const double * __restrict pthts) {
	                         
	                  register __m256d phis = _mm256_load_pd(&pphis[0]);
	                  register __m256d thts = _mm256_load_pd(&pphis[0]);            
	                  register __m256d sphis,cthts;
	                  register __m256d ahv;
	                  sphis = _mm256_sin_pd(phis);
	                  cthts = _mm256_cos_pd(thti);
	                  ahv   = gms::math::negate_ymm4r8(_mm256_div_pd(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vh_f9183_ymm4r8_u(const double * __restrict  pphis,
	                                       const double * __restrict pthts) {
	                         
	                  register __m256d phis = _mm256_loadu_pd(&pphis[0]);
	                  register __m256d thts = _mm256_loadu_pd(&pphis[0]);            
	                  register __m256d sphis,cthts;
	                  register __m256d ahv;
	                  sphis = _mm256_sin_pd(phis);
	                  cthts = _mm256_cos_pd(thti);
	                  ahv   = gms::math::negate_ymm4r8(_mm256_div_pd(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-84 
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vv_f9184_ymm4r8(const __m256d thti) {
	           
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d sthti,cthti,num,den;
	                  register __m256d avv;
	                  sthti = _mm256_sin_pd(thti);
	                  cthti = _mm256_cos_pd(thti);
	                  num   = _mm256_fmadd_pd(sthti,sthti,C10);
	                  den   = _mm256_mul_pd(cthti,cthti);
	                  avv   = _mm256_div_pd(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vv_f9184_ymm4r8_a(const double * __restrict pthti) {
	           
	                  register __m256d thti = _mm256_load_pd(&pthti[0]);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d sthti,cthti,num,den;
	                  register __m256d avv;
	                  sthti = _mm256_sin_pd(thti);
	                  cthti = _mm256_cos_pd(thti);
	                  num   = _mm256_fmadd_pd(sthti,sthti,C10);
	                  den   = _mm256_mul_pd(cthti,cthti);
	                  avv   = _mm256_div_pd(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vv_f9184_ymm4r8_u(const double * __restrict  pthti) {
	           
	                  register __m256d thti = _mm256_loadu_pd(&pthti[0]);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d sthti,cthti,num,den;
	                  register __m256d avv;
	                  sthti = _mm256_sin_pd(thti);
	                  cthti = _mm256_cos_pd(thti);
	                  num   = _mm256_fmadd_pd(sthti,sthti,C10);
	                  den   = _mm256_mul_pd(cthti,cthti);
	                  avv   = _mm256_div_pd(num,den);
	                  return (avv);
	           }
	           
	           
	           
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-85
	         */                  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_hh_f9185_ymm4r8() {
	           
	                  return _mm256_set1_pd(1.0f);
	           }  
	           
	           
	            /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-86
	         */   
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d a_vh_f9186_ymm4r8() {
	           
	                  return _mm256_setzero_pd();
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
	           __m256d rcs_vv_f9187_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti) {
	                                       
	                  const __m256d C40 = _mm256_set1_pd(4.0f);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d k04,x0,x1,l2,h2,sthti;
	                  register __m256d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_pd(k0,k0);
	                  l2    = _mm256_mul_pd(l,l);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  h2    = _mm256_mul_pd(h,h);
	                  x1    = _mm256_sin_pd(thti);
	                  strm  = _mm256_fmadd_pd(x1,x1,C10);
	                  x0    = _mm256_mul_pd(C40,k04);
	                  sthti = _mm256_mul_pd(x1,x1);
	                  arg   = _mm256_mul_pd(_mm256_mul_pd(k0,k0),
	                                        _mm256_mul_pd(sthti,l2));
	                  trm   = _mm256_mul_pd(_mm256_mul_pd(l2,h2,x0));
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(strm,strm);
	                  inve  = _mm256_rcp14_pd(earg);
	                  rcs   = _mm256_mul_pd(trm,_mm256_mul_pd(x1,inve));
	                  return (rcs);
	         }
	           
	           
	           
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9187_ymm4r8_a(const double * __restrict pk0,
	                                       const double * __restrict ph,
	                                       const double * __restrict pl,
	                                       const double * __restrict pthti) {
	                             
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d h   = _mm256_load_pd(&ph[0]);
	                  register __m256d l   = _mm256_load_pd(&pl[0]);   
	                  register __m256d thti= _mm256_load_pd(&pthti[0]);       
	                  const __m256d C40 = _mm256_set1_pd(4.0f);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d k04,x0,x1,l2,h2,sthti;
	                  register __m256d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_pd(k0,k0);
	                  l2    = _mm256_mul_pd(l,l);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  h2    = _mm256_mul_pd(h,h);
	                  x1    = _mm256_sin_pd(thti);
	                  strm  = _mm256_fmadd_pd(x1,x1,C10);
	                  x0    = _mm256_mul_pd(C40,k04);
	                  sthti = _mm256_mul_pd(x1,x1);
	                  arg   = _mm256_mul_pd(_mm256_mul_pd(k0,k0),
	                                        _mm256_mul_pd(sthti,l2));
	                  trm   = _mm256_mul_pd(_mm256_mul_pd(l2,h2,x0));
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(strm,strm);
	                  inve  = _mm256_rcp14_pd(earg);
	                  rcs   = _mm256_mul_pd(trm,_mm256_mul_pd(x1,inve));
	                  return (rcs);
	         }
	           
	            
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_vv_f9187_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti) {
	                             
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d h   = _mm256_loadu_pd(&ph[0]);
	                  register __m256d l   = _mm256_loadu_pd(&pl[0]);   
	                  register __m256d thti= _mm256_loadu_pd(&pthti[0]);       
	                  const __m256d C40 = _mm256_set1_pd(4.0f);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d k04,x0,x1,l2,h2,sthti;
	                  register __m256d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_pd(k0,k0);
	                  l2    = _mm256_mul_pd(l,l);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  h2    = _mm256_mul_pd(h,h);
	                  x1    = _mm256_sin_pd(thti);
	                  strm  = _mm256_fmadd_pd(x1,x1,C10);
	                  x0    = _mm256_mul_pd(C40,k04);
	                  sthti = _mm256_mul_pd(x1,x1);
	                  arg   = _mm256_mul_pd(_mm256_mul_pd(k0,k0),
	                                        _mm256_mul_pd(sthti,l2));
	                  trm   = _mm256_mul_pd(_mm256_mul_pd(l2,h2,x0));
	                  earg  = _mm256_exp_pd(arg);
	                  x1    = _mm256_mul_pd(strm,strm);
	                  inve  = _mm256_rcp14_pd(earg);
	                  rcs   = _mm256_mul_pd(trm,_mm256_mul_pd(x1,inve));
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
	           __m256d rcs_hh_f9188_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti) {
	                                       
	                  const __m256d C40 = _mm256_set1_pd(4.0f);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d k04,x0,x1,l2,h2,sthti,cthti;
	                  register __m256d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_pd(k0,k0);
	                  l2    = _mm256_mul_pd(l,l);
	                  strm  = _mm256_cos_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  h2    = _mm256_mul_pd(h,h);
	                  x1    = _mm256_sin_pd(thti);
	                  x0    = _mm256_mul_pd(C40,k04);
	                  sthti = _mm256_mul_pd(x1,x1);
	                  arg   = _mm256_mul_pd(_mm256_mul_pd(k0,k0),
	                                        _mm256_mul_pd(sthti,l2));
	                  trm   = _mm256_mul_pd(_mm256_mul_pd(l2,h2,x0));
	                  x1    = _mm256_mul_pd(strm,strm);
	                  earg  = _mm256_exp_pd(arg);
	                  cthti = _mm256_mul_pd(x1,x1);
	                  inve  = _mm256_rcp14_pd(earg);
	                  rcs   = _mm256_mul_pd(trm,_mm256_mul_pd(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9188_ymm4r8_a(const double * __restrict pk0,
	                                       const double * __restrict ph,
	                                       const double * __restrict pl,
	                                       const double * __restrict pthti) {
	                             
	                  register __m256d k0  = _mm256_load_pd(&pk0[0]);
	                  register __m256d h   = _mm256_load_pd(&ph[0]);
	                  register __m256d l   = _mm256_load_pd(&pl[0]);   
	                  register __m256d thti= _mm256_load_pd(&pthti[0]); 
	                                       
	                  const __m256d C40 = _mm256_set1_pd(4.0f);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d k04,x0,x1,l2,h2,sthti,cthti;
	                  register __m256d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_pd(k0,k0);
	                  l2    = _mm256_mul_pd(l,l);
	                  strm  = _mm256_cos_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  h2    = _mm256_mul_pd(h,h);
	                  x1    = _mm256_sin_pd(thti);
	                  x0    = _mm256_mul_pd(C40,k04);
	                  sthti = _mm256_mul_pd(x1,x1);
	                  arg   = _mm256_mul_pd(_mm256_mul_pd(k0,k0),
	                                        _mm256_mul_pd(sthti,l2));
	                  trm   = _mm256_mul_pd(_mm256_mul_pd(l2,h2,x0));
	                  x1    = _mm256_mul_pd(strm,strm);
	                  earg  = _mm256_exp_pd(arg);
	                  cthti = _mm256_mul_pd(x1,x1);
	                  inve  = _mm256_rcp14_pd(earg);
	                  rcs   = _mm256_mul_pd(trm,_mm256_mul_pd(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256d rcs_hh_f9188_ymm4r8_u(const double * __restrict  pk0,
	                                       const double * __restrict  ph,
	                                       const double * __restrict  pl,
	                                       const double * __restrict  pthti) {
	                             
	                  register __m256d k0  = _mm256_loadu_pd(&pk0[0]);
	                  register __m256d h   = _mm256_loadu_pd(&ph[0]);
	                  register __m256d l   = _mm256_loadu_pd(&pl[0]);   
	                  register __m256d thti= _mm256_loadu_pd(&pthti[0]); 
	                                       
	                  const __m256d C40 = _mm256_set1_pd(4.0f);
	                  const __m256d C10 = _mm256_set1_pd(1.0f);
	                  register __m256d k04,x0,x1,l2,h2,sthti,cthti;
	                  register __m256d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_pd(k0,k0);
	                  l2    = _mm256_mul_pd(l,l);
	                  strm  = _mm256_cos_pd(thti);
	                  k04   = _mm256_mul_pd(x0,x0);
	                  h2    = _mm256_mul_pd(h,h);
	                  x1    = _mm256_sin_pd(thti);
	                  x0    = _mm256_mul_pd(C40,k04);
	                  sthti = _mm256_mul_pd(x1,x1);
	                  arg   = _mm256_mul_pd(_mm256_mul_pd(k0,k0),
	                                        _mm256_mul_pd(sthti,l2));
	                  trm   = _mm256_mul_pd(_mm256_mul_pd(l2,h2,x0));
	                  x1    = _mm256_mul_pd(strm,strm);
	                  earg  = _mm256_exp_pd(arg);
	                  cthti = _mm256_mul_pd(x1,x1);
	                  inve  = _mm256_rcp14_pd(earg);
	                  rcs   = _mm256_mul_pd(trm,_mm256_mul_pd(cthti,inve));
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
	           __m256d rcs_vhhv_f9189_ymm4r8() {
	           
	                  return (_mm256_setzero_pd()); 
	          }
	          
	          
	          
	           
	      
	                                    
	        
                 
                 
                 
               
               
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_YMM4R8_HPP__*/
