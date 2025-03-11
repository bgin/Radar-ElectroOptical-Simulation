

#ifndef __GMS_RCS_COMPLEX_OBJECTS_YMM8R4_H__
#define __GMS_RCS_COMPLEX_OBJECTS_YMM8R4_H__ 240920240814

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

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM8R4_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM8R4_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM8R4_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_YMM8R4_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_YMM8R4_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_YMM8R4_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_YMM8R4_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM8R4_CREATION_DATE = "24-09-2024 08:14 PM +00200 (TUE 24 SEP 2024 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM8R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM8R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_YMM8R4_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_complex_ymm8r4.hpp"


#if !defined(__AVX512F__) || !defined(__AVX512VL__)
#error "Support of ISA AVX512F or AVX512VL is required!!"
#endif



namespace  gms {


         namespace radiolocation {
         
         
               /*
                   Work (input) arrays for kernel rcs_f8162_ymm8r4_2t_u and
                   rcs_f8162_ymm8r4_2t_a.
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
              
              
                   __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8121_ymm8r4(const __m256 gam,
                                             const __m256 phi,
                                             const __m256 k0,
                                             __m256 * __restrict D1r,
                                             __m256 * __restrict D1i,
                                            __m256 * __restrict D2r,
                                            __m256 * __restrict D2i); 
                
                
                
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8121_ymm8r4_a(const float * __restrict pgam,
                                               const float * __restrict pphi,
                                               const float * __restrict pk0,
                                               float * __restrict D1r,
                                               float * __restrict D1i,
                                               float * __restrict D2r,
                                               float * __restrict D2i); 
                
                
                
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8121_ymm8r4_u(const float * __restrict  pgam,
                                               const float * __restrict  pphi,
                                               const float * __restrict  pk0,
                                               float * __restrict  D1r,
                                               float * __restrict  D1i,
                                               float * __restrict  D2r,
                                               float * __restrict  D2i); 
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter singly diffracted far-zone fields (E,H).
                    Formula: 8.1-19, 8.1-20
                
                */
                
                
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void EsHs_f811920_ymm8r4(    const __m256 betai,
                                                 const __m256 betas,
                                                 const __m256 gam,
                                                 const __m256 phi,
                                                 const __m256 k0,
                                                 const __m256 r,
                                                 const __m256 rho,
                                                 const __m256 psi,
                                                 __m256 * __restrict Esr,
                                                 __m256 * __restrict Esi,
                                                 __m256 * __restrict Hsr,
                                                 __m256 * __restrict Hsi); 
            
                
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void EsHs_f811920_ymm8r4_a(  const float * __restrict pbetai,
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
                                                 float * __restrict Hsi); 
            
            
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void EsHs_f811920_ymm8r4_u(  const float * __restrict  pbetai,
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
                                                 float * __restrict  Hsi); 
            
            
            /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Ray normal-incidence to one of edge faces.
                    Formula: 8.1-24
            */
            
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8124_ymm8r4(const __m256 k0,
                                               const __m256 gam,
                                               __m256 * __restrict D1r,
                                               __m256 * __restrict D1i,
                                               __m256 * __restrict D2r,
                                               __m256 * __restrict D2i); 
                
                
                   
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8124_ymm8r4_a(const  float * __restrict pk0,
                                                 const  float * __restrict pgam,
                                                 float * __restrict D1r,
                                                 float * __restrict D1i,
                                                 float * __restrict D2r,
                                                 float * __restrict D2i); 
                
                
                 
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8124_ymm8r4_u(const  float * __restrict  pk0,
                                                 const  float * __restrict pgam,
                                                 float * __restrict  D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict  D2r,
                                                 float * __restrict  D2i); 
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Backscatter direction axial caustic (for slightly diffracted rays).
                    Formula: 8.1-26
                */
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_Ddiff_f8126_ymm8r4(const __m256 gam,
                                             const __m256 phi,
                                             const __m256 k0,
                                             __m256 * __restrict Dr,
                                             __m256 * __restrict Di); 
                
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_Ddiff_f8126_ymm8r4_a(const float * __restrict pgam,
                                               const float * __restrict pphi,
                                               const float * __restrict pk0,
                                               float * __restrict  Dr,
                                               float * __restrict  Di); 
                
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_Ddiff_f8126_ymm8r4_u(const float * __restrict  pgam,
                                                   const float * __restrict pphi,
                                                   const float * __restrict  pk0,
                                                   float * __restrict  Dr,
                                                   float * __restrict  Di); 
                
                
                   /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Backscatter direction axial caustic (for slightly diffracted rays).
                    Scattered Electric and Magnetic fields.
                    Formula: 8.1-25
                */
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void EsHs_f8125_ymm8r4(const __m256 a,
                                           const __m256 k0,
                                           const __m256 r,
                                           const __m256 gam,
                                           const __m256 phi,
                                           const __m256 psi,
                                           __m256 * __restrict Esr,
                                           __m256 * __restrict Esi,
                                           __m256 * __restrict Hsr,
                                           __m256 * __restrict Hsi); 
                
                
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void EsHs_f8125_ymm8r4_a(const float * __restrict pa,
                                             const float * __restrict pk0,
                                             const float * __restrict pr,
                                             const float * __restrict pgam,
                                             const float * __restrict pphi,
                                             const float * __restrict ppsi,
                                             float * __restrict Esr,
                                             float * __restrict Esi,
                                             float * __restrict Hsr,
                                             float * __restrict Hsi); 
                
            
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void EsHs_f8125_ymm8r4_u(const float * __restrict  pa,
                                             const float * __restrict  pk0,
                                             const float * __restrict  pr,
                                             const float * __restrict pgam,
                                             const float * __restrict  pphi,
                                             const float * __restrict  ppsi,
                                             float * __restrict  Esr,
                                             float * __restrict  Esi,
                                             float * __restrict  Hsr,
                                             float * __restrict  Hsi); 
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Doubly and high order diffracted rays --
                    bistatic diffraction coefficients.
                    Formula: 8.1-27  
                
                */
                
                
                   __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8127_ymm8r4(const __m256 k0,
                                               const __m256 gam,
                                               const __m256 phi1,
                                               const __m256 phi2,
                                               __m256 * __restrict D1r,
                                               __m256 * __restrict D1i,
                                               __m256 * __restrict D2r,
                                               __m256 * __restrict D2i); 
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8127_ymm8r4_a(const float * __restrict pk0,
                                                 const float * __restrict pgam,
                                                 const float * __restrict pphi1,
                                                 const float * __restrict pphi2,
                                                 float * __restrict  D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict  D2r,
                                                 float * __restrict  D2i); 
                
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   void coef_D12_f8127_ymm8r4_u(const float * __restrict  pk0,
                                                 const float * __restrict  pgam,
                                                 const float * __restrict  pphi1,
                                                 const float * __restrict  pphi2,
                                                 float * __restrict   D1r,
                                                 float * __restrict  D1i,
                                                 float * __restrict   D2r,
                                                 float * __restrict  D2i); 
                
                
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


                 
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_u(const float * __restrict pdAdl,
                                             const float *  __restrict pdl,
                                             const float   k0,
                                             const float   l); 
                  
                  
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_a(const float * __restrict pdAdl,
                                             const float * __restrict pdl,
                                             const float   k0,
                                             const float   l); 
                  
                  
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
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_avint_u(const float * __restrict pdAdl,
                                                   const float *  __restrict pdl,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri); 
                  
                  
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_avint_a(const float * __restrict pdAdl,
                                                   const float * __restrict pdl,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri); 
                  
                  
                  
                  
                  
                   /*
                       Adachi expression for axial-incidence
                       of backscatter RCS for entire scatterer length.
                       Shall be used in case of thin long axially symetric 
                       bodies e.g. 'ogives,double-cones, etc.,'
                       Vectorization of an integrand.
                       Case of small integrand -- single-threaded execution.
                       Formula 8.1-62
                */
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_u(const float * __restrict  pdAdl,
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
                                             const int32_t NTAB); 
               
               
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_a(const float * __restrict pdAdl,
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
                                             const int32_t NTAB); 
               
               
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
                
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_avint_u(const float * __restrict  pdAdl,
                                                   const float * __restrict  pdl,
                                                   float * __restrict  intr,
                                                   float * __restrict  inti,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB); 
               
               
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_avint_a(const float * __restrict pdAdl,
                                                   const float * __restrict pdl,
                                                   float * __restrict intr,
                                                   float * __restrict inti,
                                                   const float   k0,
                                                   const float   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB); 
               
               
               
               
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
                

    
                
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_cspint2t_u(const float * __restrict  pdAdl,
                                                     const float * __restrict  pdl,
                                                     float * __restrict  intr,
                                                     float * __restrict  inti,
                                                     struct RCS_F8162_DATA w,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB); 
               
               
               
                   
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_cspint2t_a(const float * __restrict pdAdl,
                                                     const float * __restrict pdl,
                                                     float * __restrict intr,
                                                     float * __restrict inti,
                                                     struct RCS_F8162_DATA w,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB); 
               
               
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
                
                
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_avint2t_u(const float * __restrict  pdAdl,
                                                     const float * __restrict  pdl,
                                                     float * __restrict  intr,
                                                     float * __restrict  inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB); 
               
                  
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   float rcs_f8162_ymm8r4_avint2t_a(const float * __restrict pdAdl,
                                                     const float * __restrict pdl,
                                                     float * __restrict intr,
                                                     float * __restrict inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const float   k0,
                                                     const float   l,
                                                     const int32_t NTAB); 
               
               
               
               
               
               /*
                     High frequency approximations.
                     Rounded-tip cone total nose-on
                     backscatter RCS.
                     Formula 8.1-93
               */
               
               
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m256 rcs_f8193_ymm8r4(const __m256 b,
                                            const __m256 a,
                                            const __m256 k0,
                                            const __m256 alp,
                                            const __m256 l); 
                                            
                 
                 
                 
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m256 rcs_f8193_ymm8r4_a(const float * __restrict pb,
                                              const float * __restrict pa,
                                              const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pl);
                 
                 
                 
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m256 rcs_f8193_ymm8r4_u(const float * __restrict  pb,
                                              const float * __restrict  pa,
                                              const float * __restrict  pk0,
                                              const float * __restrict  palp,
                                              const float * __restrict _pl); 
                                              
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for |theta| = PI/2-alpha
                     Formula 8.1-96
                 */
                 
                 
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m256 rcs_f8196_ymm8r4(const __m256 k0,
                                            const __m256 alp,
                                            const __m256 a,
                                            const __m256 b); 
                 
                 
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m256 rcs_f8196_ymm8r4_a(const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pa,
                                              const float * __restrict pb); 
                 
                 
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __m256 rcs_f8196_ymm8r4_u(const float * __restrict pk0,
                                              const float * __restrict palp,
                                              const float * __restrict pa,
                                              const float * __restrict pb); 
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for 0<|theta|<alpha
                     Perpendicular RCS.
                     Formula 8.1-94
                 */
                 
                 
             /*      __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_perpendic_f8194_ymm8r4(const __m256 h,
	                                        const __m256 l,
	                                        const __m256 b,
	                                        const __m256 a,
	                                        const __m256 k0,
	                                        const __m256 tht,
	                                        const __m256 alp) {
	                                 
	                                  
	                 const __m256 C314159265358979323846264338328  = 
                                                     _mm256_set1_ps(3.14159265358979323846264338328f); 
                         const __m256 C1772453850905516027298167483341 = 
                                                     _mm256_set1_ps(1.772453850905516027298167483341f);
                         const __m256 C078539816339744830961566084582  = 
                                                     _mm256_set1_ps(0.78539816339744830961566084582f);
                         const __m256 C10                              = 
                                                     _mm256_set1_ps(1.0f);  
                         const __m256 C15                              = 
                                                     _mm256_set1_ps(1.5f); 
                         const __m256 C05                              = 
                                                     _mm256_set1_ps(0.5f);
                         const __m256 C20                              =
                                                     _mm256_set1_ps(2.0f);
                          __m256 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m256 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m256 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m256 cpin1,cpin2,trm1,trm2,rcs;
                         __m256 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_ps(h,_mm256_add_ps(l,b));
                         sint = _mm256_sin_ps(tht);
                         k02  = _mm256_add_ps(k0,k0);
                         n    = _mm256_mul_ps(C15,_mm256_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_ps(sint);
                         a2   = _mm256_mul_ps(a,C05);
                         ear  = _mm256_setzero_ps();
                         sk02 = _mm256_sqrt_ps(_mm256_mul_ps(k0,C05));
                         x0   = _mm256_mul_ps(hlb,_mm256_sub_ps(cost,b));
                         invn = _mm256_rcp14_ps(n);
                         //x2   = _mm256_mul_ps(a,C05);
                         eai  = _mm256_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_ps(_mm256_mul_ps(a2,csct));
                         atant= _mm256_mul_ps(a,tant);
                         cost = _mm256_cos_ps(tht);
                         x0   = _mm256_mul_ps(b,C1772453850905516027298167483341);
                         cexp_ymm8c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_ps(x0,cer1);
                         spin = _mm256_sin_ps(pin);
                         cei1 = _mm256_mul_ps(x0,cei1);
                         cpin = _mm256_cos_ps(pin);
                         x1   = _mm256_mul_ps(_mm256_sub_ps(h,atant),cost);
                         eai2 = _mm256_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_ymm8c4(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_mul_ps(x0,sacs);
                         cer2 = _mm256_mul_ps(cer2,x1);
                         cei2 = _mm256_mul_ps(cei2,x1);
                         cpin1= _mm256_rcp14_ps(_mm256_sub_ps(cpin,C10));
                         x2   = _mm256_mul_ps(C20,_mm256_add_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x2,invn));
                         x3   = _mm256_rcp14_ps(_mm256_sub_ps(cpin,cpin2));
                         trm1 = _mm256_sub_ps(cpin1,x3);
                         cmul_ymm8c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_ps(t0r,trm1);
                         t0i  = _mm256_mul_ps(t0i,trm1);
                         x0   = _mm256_mul_ps(C20,_mm256_sub_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x0,invn));
                         x1   = _mm256_rcp14_ps(cpin2);
                         trm2 = _mm256_sub_ps(cpin1,x1);
                         x2   = _mm256_fmadd_ps(cost,_mm256_mul_ps(k02,
                                                               _mm256_add_ps(h,atant)));
                         eai3 = _mm256_add_ps(C078539816339744830961566084582,x2);
                         cexp_ymm8c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_sqrt_ps(_mm256_mul_ps(gms::math::
                                                                  negate_ymm8r4(a2),csct));
                         x2   = _mm256_mul_ps(x0,x1);
                         cer3 = _mm256_mul_ps(_mm256_mul_ps(cer3,x2),trm2);
                         cei3 = _mm256_mul_ps(_mm256_mul_ps(cei3,x2),trm2);
                         cmul_ymm8c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm8c4(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	                                  
                    __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_perpendic_f8194_ymm8r4_a(const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pb,
	                                         const float * __restrict pa,
	                                         const float * __restrict pk0,
	                                         const float * __restrict ptht,
	                                         const float * __restrict palp) {
	                                 
	                  
	                  __m256 h  = _mm256_load_ps(&ph[0]);
	                  __m256 l  = _mm256_load_ps(&pl[0]); 
	                  __m256 b  = _mm256_load_ps(&pb[0]);   
	                  __m256 a  = _mm256_load_ps(&pa[0]);  
	                  __m256 k0 = _mm256_load_ps(&pk0[0]);
	                  __m256 tht= _mm256_load_ps(&ptht[0]); 
	                  __m256 alp= _mm256_load_ps(&palp[0]);        
	                 const __m256 C314159265358979323846264338328  = 
                                                     _mm256_set1_ps(3.14159265358979323846264338328f); 
                         const __m256 C1772453850905516027298167483341 = 
                                                     _mm256_set1_ps(1.772453850905516027298167483341f);
                         const __m256 C078539816339744830961566084582  = 
                                                     _mm256_set1_ps(0.78539816339744830961566084582f);
                         const __m256 C10                              = 
                                                     _mm256_set1_ps(1.0f);  
                         const __m256 C15                              = 
                                                     _mm256_set1_ps(1.5f); 
                         const __m256 C05                              = 
                                                     _mm256_set1_ps(0.5f);
                         const __m256 C20                              =
                                                     _mm256_set1_ps(2.0f);
                         __m256 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m256 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m256 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m256 cpin1,cpin2,trm1,trm2,rcs;
                         __m256 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_ps(h,_mm256_add_ps(l,b));
                         sint = _mm256_sin_ps(tht);
                         k02  = _mm256_add_ps(k0,k0);
                         n    = _mm256_mul_ps(C15,_mm256_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_ps(sint);
                         a2   = _mm256_mul_ps(a,C05);
                         ear  = _mm256_setzero_ps();
                         sk02 = _mm256_sqrt_ps(_mm256_mul_ps(k0,C05));
                         x0   = _mm256_mul_ps(hlb,_mm256_sub_ps(cost,b));
                         invn = _mm256_rcp14_ps(n);
                         //x2   = _mm256_mul_ps(a,C05);
                         eai  = _mm256_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_ps(_mm256_mul_ps(a2,csct));
                         atant= _mm256_mul_ps(a,tant);
                         cost = _mm256_cos_ps(tht);
                         x0   = _mm256_mul_ps(b,C1772453850905516027298167483341);
                         cexp_ymm8c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_ps(x0,cer1);
                         spin = _mm256_sin_ps(pin);
                         cei1 = _mm256_mul_ps(x0,cei1);
                         cpin = _mm256_cos_ps(pin);
                         x1   = _mm256_mul_ps(_mm256_sub_ps(h,atant),cost);
                         eai2 = _mm256_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_ymm8c4(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_mul_ps(x0,sacs);
                         cer2 = _mm256_mul_ps(cer2,x1);
                         cei2 = _mm256_mul_ps(cei2,x1);
                         cpin1= _mm256_rcp14_ps(_mm256_sub_ps(cpin,C10));
                         x2   = _mm256_mul_ps(C20,_mm256_add_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x2,invn));
                         x3   = _mm256_rcp14_ps(_mm256_sub_ps(cpin,cpin2));
                         trm1 = _mm256_sub_ps(cpin1,x3);
                         cmul_ymm8c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_ps(t0r,trm1);
                         t0i  = _mm256_mul_ps(t0i,trm1);
                         x0   = _mm256_mul_ps(C20,_mm256_sub_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x0,invn));
                         x1   = _mm256_rcp14_ps(cpin2);
                         trm2 = _mm256_sub_ps(cpin1,x1);
                         x2   = _mm256_fmadd_ps(cost,_mm256_mul_ps(k02,
                                                               _mm256_add_ps(h,atant)));
                         eai3 = _mm256_add_ps(C078539816339744830961566084582,x2);
                         cexp_ymm8c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_sqrt_ps(_mm256_mul_ps(gms::math::
                                                                  negate_ymm8r4(a2),csct));
                         x2   = _mm256_mul_ps(x0,x1);
                         cer3 = _mm256_mul_ps(_mm256_mul_ps(cer3,x2),trm2);
                         cei3 = _mm256_mul_ps(_mm256_mul_ps(cei3,x2),trm2);
                         cmul_ymm8c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm8c4(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_perpendic_f8194_ymm8r4_u(    const float * __restrict  ph,
	                                                    const float * __restrict  pl,
	                                                    const float * __restrict  pb,
	                                                    const float * __restrict  pa,
	                                                    const float * __restrict  pk0,
	                                                    const float * __restrict  ptht,
	                                                    const float * __restrict  palp) {
	                                 
	                  
	                  __m256 h  = _mm256_loadu_ps(&ph[0]);
	                  __m256 l  = _mm256_loadu_ps(&pl[0]); 
	                  __m256 b  = _mm256_loadu_ps(&pb[0]);   
	                  __m256 a  = _mm256_loadu_ps(&pa[0]);  
	                  __m256 k0 = _mm256_loadu_ps(&pk0[0]);
	                  __m256 tht= _mm256_loadu_ps(&ptht[0]); 
	                  __m256 alp= _mm256_loadu_ps(&palp[0]);        
	                 const __m256 C314159265358979323846264338328  = 
                                                     _mm256_set1_ps(3.14159265358979323846264338328f); 
                         const __m256 C1772453850905516027298167483341 = 
                                                     _mm256_set1_ps(1.772453850905516027298167483341f);
                         const __m256 C078539816339744830961566084582  = 
                                                     _mm256_set1_ps(0.78539816339744830961566084582f);
                         const __m256 C10                              = 
                                                     _mm256_set1_ps(1.0f);  
                         const __m256 C15                              = 
                                                     _mm256_set1_ps(1.5f); 
                         const __m256 C05                              = 
                                                     _mm256_set1_ps(0.5f);
                         const __m256 C20                              =
                                                     _mm256_set1_ps(2.0f);
                         __m256 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m256 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m256 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m256 cpin1,cpin2,trm1,trm2,rcs;
                         __m256 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_ps(h,_mm256_add_ps(l,b));
                         sint = _mm256_sin_ps(tht);
                         k02  = _mm256_add_ps(k0,k0);
                         n    = _mm256_mul_ps(C15,_mm256_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_ps(sint);
                         a2   = _mm256_mul_ps(a,C05);
                         ear  = _mm256_setzero_ps();
                         sk02 = _mm256_sqrt_ps(_mm256_mul_ps(k0,C05));
                         x0   = _mm256_mul_ps(hlb,_mm256_sub_ps(cost,b));
                         invn = _mm256_rcp14_ps(n);
                         //x2   = _mm256_mul_ps(a,C05);
                         eai  = _mm256_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_ps(_mm256_mul_ps(a2,csct));
                         atant= _mm256_mul_ps(a,tant);
                         cost = _mm256_cos_ps(tht);
                         x0   = _mm256_mul_ps(b,C1772453850905516027298167483341);
                         cexp_ymm8c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_ps(x0,cer1);
                         spin = _mm256_sin_ps(pin);
                         cei1 = _mm256_mul_ps(x0,cei1);
                         cpin = _mm256_cos_ps(pin);
                         x1   = _mm256_mul_ps(_mm256_sub_ps(h,atant),cost);
                         eai2 = _mm256_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_ymm8c4(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_mul_ps(x0,sacs);
                         cer2 = _mm256_mul_ps(cer2,x1);
                         cei2 = _mm256_mul_ps(cei2,x1);
                         cpin1= _mm256_rcp14_ps(_mm256_sub_ps(cpin,C10));
                         x2   = _mm256_mul_ps(C20,_mm256_add_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x2,invn));
                         x3   = _mm256_rcp14_ps(_mm256_sub_ps(cpin,cpin2));
                         trm1 = _mm256_sub_ps(cpin1,x3);
                         cmul_ymm8c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_ps(t0r,trm1);
                         t0i  = _mm256_mul_ps(t0i,trm1);
                         x0   = _mm256_mul_ps(C20,_mm256_sub_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x0,invn));
                         x1   = _mm256_rcp14_ps(cpin2);
                         trm2 = _mm256_sub_ps(cpin1,x1);
                         x2   = _mm256_fmadd_ps(cost,_mm256_mul_ps(k02,
                                                               _mm256_add_ps(h,atant)));
                         eai3 = _mm256_add_ps(C078539816339744830961566084582,x2);
                         cexp_ymm8c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_sqrt_ps(_mm256_mul_ps(gms::math::
                                                                  negate_ymm8r4(a2),csct));
                         x2   = _mm256_mul_ps(x0,x1);
                         cer3 = _mm256_mul_ps(_mm256_mul_ps(cer3,x2),trm2);
                         cei3 = _mm256_mul_ps(_mm256_mul_ps(cei3,x2),trm2);
                         cmul_ymm8c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm8c4(t1r,t1i);
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
	           __m256 rcs_parallel_f8194_ymm8r4(const __m256 h,
	                                             const __m256 l,
	                                             const __m256 b,
	                                             const __m256 a,
	                                             const __m256 k0,
	                                             const __m256 tht,
	                                             const __m256 alp) {
	                                 
	                                  
	                 const __m256 C314159265358979323846264338328  = 
                                                     _mm256_set1_ps(3.14159265358979323846264338328f); 
                         const __m256 C1772453850905516027298167483341 = 
                                                     _mm256_set1_ps(1.772453850905516027298167483341f);
                         const __m256 C078539816339744830961566084582  = 
                                                     _mm256_set1_ps(0.78539816339744830961566084582f);
                         const __m256 C10                              = 
                                                     _mm256_set1_ps(1.0f);  
                         const __m256 C15                              = 
                                                     _mm256_set1_ps(1.5f); 
                         const __m256 C05                              = 
                                                     _mm256_set1_ps(0.5f);
                         const __m256 C20                              =
                                                     _mm256_set1_ps(2.0f);
                          __m256 pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                          __m256 ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                          __m256 cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                          __m256 cpin1,cpin2,trm1,trm2,rcs;
                         __m256 t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm256_sub_ps(h,_mm256_add_ps(l,b));
                         sint = _mm256_sin_ps(tht);
                         k02  = _mm256_add_ps(k0,k0);
                         n    = _mm256_mul_ps(C15,_mm256_div_ps(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm256_rcp14_ps(sint);
                         a2   = _mm256_mul_ps(a,C05);
                         ear  = _mm256_setzero_ps();
                         sk02 = _mm256_sqrt_ps(_mm256_mul_ps(k0,C05));
                         x0   = _mm256_mul_ps(hlb,_mm256_sub_ps(cost,b));
                         invn = _mm256_rcp14_ps(n);
                         //x2   = _mm256_mul_ps(a,C05);
                         eai  = _mm256_mul_ps(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm256_mul_ps(C314159265358979323846264338328,invn);
                         sacs = _mm256_sqrt_ps(_mm256_mul_ps(a2,csct));
                         atant= _mm256_mul_ps(a,tant);
                         cost = _mm256_cos_ps(tht);
                         x0   = _mm256_mul_ps(b,C1772453850905516027298167483341);
                         cexp_ymm8c4(ear,eai,&cer1,&cei1);
                         cer1 = _mm256_mul_ps(x0,cer1);
                         spin = _mm256_sin_ps(pin);
                         cei1 = _mm256_mul_ps(x0,cei1);
                         cpin = _mm256_cos_ps(pin);
                         x1   = _mm256_mul_ps(_mm256_sub_ps(h,atant),cost);
                         eai2 = _mm256_fmadd_ps(k02,x1,C078539816339744830961566084582);
                         cexp_ymm8c4(ear,eai2,&cer2,&cei);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_mul_ps(x0,sacs);
                         cer2 = _mm256_mul_ps(cer2,x1);
                         cei2 = _mm256_mul_ps(cei2,x1);
                         cpin1= _mm256_rcp14_ps(_mm256_sub_ps(cpin,C10));
                         x2   = _mm256_mul_ps(C20,_mm256_add_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x2,invn));
                         x3   = _mm256_rcp14_ps(_mm256_sub_ps(cpin,cpin2));
                         trm1 = _mm256_add_ps(cpin1,x3);
                         cmul_ymm8c4(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm256_mul_ps(t0r,trm1);
                         t0i  = _mm256_mul_ps(t0i,trm1);
                         x0   = _mm256_mul_ps(C20,_mm256_sub_ps(alp,tht));
                         cpin2= _mm256_cos_ps(_mm256_mul_ps(x0,invn));
                         x1   = _mm256_rcp14_ps(cpin2);
                         trm2 = _mm256_add_ps(cpin1,x1);
                         x2   = _mm256_fmadd_ps(cost,_mm256_mul_ps(k02,
                                                               _mm256_add_ps(h,atant)));
                         eai3 = _mm256_add_ps(C078539816339744830961566084582,x2);
                         cexp_ymm8c4(ear,ea3,&cer3,&cei3);
                         x0   = _mm256_div_ps(spin,_mm256_mul_ps(n,sk02));
                         x1   = _mm256_sqrt_ps(_mm256_mul_ps(gms::math::
                                                                  negate_ymm8r4(a2),csct));
                         x2   = _mm256_mul_ps(x0,x1);
                         cer3 = _mm256_mul_ps(_mm256_mul_ps(cer3,x2),trm2);
                         cei3 = _mm256_mul_ps(_mm256_mul_ps(cei3,x2),trm2);
                         cmul_ymm8c4(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_ymm8c4(t1r,t1i);
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
	        
	        
	        
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 coef_Bg_f9137_ymm8r4(const __m256 A,
	                                        const __m256 N,
	                                        const __m256 k0,
	                                        const __m256 epsr,
	                                        const __m256 epsi,
	                                        const __m256 thti,
	                                        const __m256 thts,
	                                        const __m256 phis,
	                                        const int pol); 
	       
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 coef_Bg_f9137_ymm8r4_a(const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol); 
	                                          
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 coef_Bg_f9137_ymm8r4_u(const float * __restrict  pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol); 
	       
	       
	       /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (hh) polarized per unit surface area.
	             Formula 9.1-33
	       
	       */
	       
	       
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9133_ymm8r4( const __m256 A,
	                                        const __m256 N,
	                                        const __m256 k0,
	                                        const __m256 epsr,
	                                        const __m256 epsi,
	                                        const __m256 thti,
	                                        const __m256 thts,
	                                        const __m256 phis,
	                                        const int pol); 
	       
	       
	       
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9133_ymm8r4_a( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol); 
	       
	       
	       
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9133_ymm8r4_u( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol); 
	       
	       
	        
	       /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (vh) polarized per unit surface area.
	             Formula 9.1-34
	       
	       */
	       
	       
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vh_f9134_ymm8r4( const __m256 A,
	                                        const __m256 N,
	                                        const __m256 k0,
	                                        const __m256 epsr,
	                                        const __m256 epsi,
	                                        const __m256 thti,
	                                        const __m256 thts,
	                                        const __m256 phis,
	                                        const int pol); 
	         
	         
	         
	          
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vh_f9134_ymm8r4_a( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol); 
	         
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vh_f9134_ymm8r4_u( const float * __restrict  pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol); 
	         
	         
	           /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (hv) polarized per unit surface area.
	             Formula 9.1-35
	       
	       */
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hv_f9135_ymm8r4( const __m256 A,
	                                        const __m256 N,
	                                        const __m256 k0,
	                                        const __m256 epsr,
	                                        const __m256 epsi,
	                                        const __m256 thti,
	                                        const __m256 thts,
	                                        const __m256 phis,
	                                        const int pol); 
	         
	         
	         
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hv_f9135_ymm8r4_a(  const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol); 
	         
	         
	            
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256  rcs_hv_f9135_ymm8r4_u(const float * __restrict pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol); 
	         
	         
	            /*
	             Model 9B4 (Peake's Model)
	             Model resembling many natural grass-like structures like
	             forests,grass,wheat fields, etc.
	             Bistatic RCS (vv) polarized per unit surface area.
	             Formula 9.1-36
	       
	       */
	       
	       
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9136_ymm8r4( const __m256 A,
	                                        const __m256 N,
	                                        const __m256 k0,
	                                        const __m256 epsr,
	                                        const __m256 epsi,
	                                        const __m256 thti,
	                                        const __m256 thts,
	                                        const __m256 phis,
	                                        const int pol); 
	         
	         
	         
	          
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9136_ymm8r4_a( const float * __restrict pA,
	                                          const float * __restrict pN,
	                                          const float * __restrict pk0,
	                                          const float * __restrict pepsr,
	                                          const float * __restrict pepsi,
	                                          const float * __restrict pthti,
	                                          const float * __restrict pthts,
	                                          const float * __restrict pphis,
	                                          const int pol ); 
	         
	         
	             
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9136_ymm8r4_u( const float * __restrict  pA,
	                                          const float * __restrict  pN,
	                                          const float * __restrict  pk0,
	                                          const float * __restrict  pepsr,
	                                          const float * __restrict  pepsi,
	                                          const float * __restrict  pthti,
	                                          const float * __restrict  pthts,
	                                          const float * __restrict  pphis,
	                                          const int pol ); 
	         
	         
	         
	        /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (vv) polarization.
	            Formula 9.1-74
	        */
	        
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9174_ymm8r4(const __m256 k0,
	                                       const __m256 h,
	                                       const __m256 l,
	                                       const __m256 thti,
	                                       const __m256 epsr,
	                                       const __m256 epsi,
	                                       const __m256 mur,
	                                       const __m256 mui); 
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9174_ymm8r4_a(const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui); 
	       
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9174_ymm8r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui); 
	                                         
	       
	       
	        /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hh) polarization.
	            Formula 9.1-75
	        */
	        
	        
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9175_ymm8r4(const __m256 k0,
	                                       const __m256 h,
	                                       const __m256 l,
	                                       const __m256 thti,
	                                       const __m256 epsr,
	                                       const __m256 epsi,
	                                       const __m256 mur,
	                                       const __m256 mui); 
	       
	            
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9175_ymm8r4_a(  const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui); 
	                                         
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9175_ymm8r4_u(  const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui);
	       
	         /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hv) polarization.
	            Formula 9.1-76
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_hv_f9176_ymm8r4() { 
	           
	                return _mm256_setzero_ps();
	         } 
	         
	         
	         
	         /*
	            Exponential surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (vv) polarization.
	            Formula 9.1-77
	        */
	        
	        
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9177_ymm8r4(const __m256 k0,
	                                       const __m256 h,
	                                       const __m256 l,
	                                       const __m256 thti,
	                                       const __m256 epsr,
	                                       const __m256 epsi,
	                                       const __m256 mur,
	                                       const __m256 mui); 
	       
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9177_ymm8r4_a(const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui); 
	       
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_vv_f9177_ymm8r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui); 
	                                         
	       
	       
	         /*
	            Exponential surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hh) polarization.
	            Formula 9.1-78
	        */
	        
	        
	        
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9178_ymm8r4(const __m256 k0,
	                                       const __m256 h,
	                                       const __m256 l,
	                                       const __m256 thti,
	                                       const __m256 epsr,
	                                       const __m256 epsi,
	                                       const __m256 mur,
	                                       const __m256 mui); 
	       
	       
	             __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9178_ymm8r4_a(const float * __restrict pk0,
	                                         const float * __restrict ph,
	                                         const float * __restrict pl,
	                                         const float * __restrict pthti,
	                                         const float * __restrict pepsr,
	                                         const float * __restrict pepsi,
	                                         const float * __restrict pmur,
	                                         const float * __restrict pmui); 
	       
	       
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
	           __m256 rcs_hh_f9178_ymm8r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti,
	                                         const float * __restrict  pepsr,
	                                         const float * __restrict  pepsi,
	                                         const float * __restrict  pmur,
	                                         const float * __restrict  pmui);
	       
	         /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hv) polarization.
	            Formula 9.1-79
	        */
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_hv_f9179_ymm8r4() { 
	           
	                return _mm256_setzero_ps();
	         } 
	         
	         
	         /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-80
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vv_f9180_ymm8r4(const __m256 thti,
	                                         const __m256 thts,
	                                         const __m256 phis) {
	                                         
	                   __m256 sthti,sthts,cthti,cthts,cphis;
	                   __m256 num,den,avv;  
	                  sthti = _mm256_sin_ps(thti);
	                  cthti = _mm256_cos_ps(thti);
	                  cphis = _mm256_cos_ps(phis);
	                  cthts = _mm256_cos_ps(thts);
	                  sthts = _mm256_sin_ps(thts);
	                  num   = _mm256_fmsub_ps(sthti,sthts,cphis);
	                  den   = _mm256_mul_ps(cthti,cthts);
	                  avv   = _mm256_div_ps(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vv_f9180_ymm8r4_a(  const float * __restrict pthti,
	                                         const float * __restrict pthts,
	                                         const float * __restrict pphis) {
	                     
	                   __m256 thti = _mm256_load_ps(&pthti[0]);
	                   __m256 thts = _mm256_load_ps(&pthts[0]);
	                   __m256 phis = _mm256_load_ps(&phis[0]);                    
	                   __m256 sthti,sthts,cthti,cthts,cphis;
	                   __m256 num,den,avv;  
	                  sthti = _mm256_sin_ps(thti);
	                  cthti = _mm256_cos_ps(thti);
	                  cphis = _mm256_cos_ps(phis);
	                  cthts = _mm256_cos_ps(thts);
	                  sthts = _mm256_sin_ps(thts);
	                  num   = _mm256_fmsub_ps(sthti,sthts,cphis);
	                  den   = _mm256_mul_ps(cthti,cthts);
	                  avv   = _mm256_div_ps(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vv_f9180_ymm8r4_u(  const float * __restrict  pthti,
	                                         const float * __restrict  pthts,
	                                         const float * __restrict  pphis) {
	                     
	                   __m256 thti = _mm256_loadu_ps(&pthti[0]);
	                   __m256 thts = _mm256_loadu_ps(&pthts[0]);
	                   __m256 phis = _mm256_loadu_ps(&phis[0]);                    
	                   __m256 sthti,sthts,cthti,cthts,cphis;
	                   __m256 num,den,avv;  
	                  sthti = _mm256_sin_ps(thti);
	                  cthti = _mm256_cos_ps(thti);
	                  cphis = _mm256_cos_ps(phis);
	                  cthts = _mm256_cos_ps(thts);
	                  sthts = _mm256_sin_ps(thts);
	                  num   = _mm256_fmsub_ps(sthti,sthts,cphis);
	                  den   = _mm256_mul_ps(cthti,cthts);
	                  avv   = _mm256_div_ps(num,den);
	                  return (avv);        
	          }
	          
	          
	          /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-81
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hv_f9181_ymm8r4(const __m256 phis,
	                                     const __m256 thti) {
	                                     
	                   __m256 sphis,cthti;
	                   __m256 ahv;
	                  sphis = _mm256_sin_ps(phis);
	                  cthti = _mm256_cos_ps(thti);
	                  ahv   = _mm256_div_ps(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hv_f9181_ymm8r4_a(const float * __restrict pphis,
	                                     const float * __restrict pthti) {
	                                     
	                   __m256 thti = _mm256_load_ps(&pthti[0]);
	                   __m256 phis = _mm256_load_ps(&phis[0]); 
	                   __m256 sphis,cthti;
	                   __m256 ahv;
	                  sphis = _mm256_sin_ps(phis);
	                  cthti = _mm256_cos_ps(thti);
	                  ahv   = _mm256_div_ps(sphis,cthti);
	                  return (ahv);                           
	         }
	          
	       
	        
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hv_f9181_ymm8r4_u(const float * __restrict pphis,
	                                     const float * __restrict  pthti) {
	                                     
	                   __m256 thti = _mm256_loadu_ps(&pthti[0]);
	                   __m256 phis = _mm256_loadu_ps(&phis[0]); 
	                   __m256 sphis,cthti;
	                   __m256 ahv;
	                  sphis = _mm256_sin_ps(phis);
	                  cthti = _mm256_cos_ps(thti);
	                  ahv   = _mm256_div_ps(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	            /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-82
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hh_f9182_ymm8r4(const __m256 phis) {
	                 
	                   __m256 cphis,ahh;
	                  cphis = _mm256_cos_ps(phis);
	                  ahh   = negate_ymm8r4(cphis);
	                  return (ahh);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hh_f9182_ymm8r4_a(const float * __restrict pphis) {
	                 
	                   __m256 phis = _mm256_load_ps(&pphis[0]);
	                   __m256 cphis,ahh;
	                  cphis = _mm256_cos_ps(phis);
	                  ahh   = negate_ymm8r4(cphis);
	                  return (ahh);
	          }
	        
	       
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hh_f9182_ymm8r4_u(const float * __restrict  pphis) {
	                 
	                   __m256 phis = _mm256_loadu_ps(&pphis[0]);
	                   __m256 cphis,ahh;
	                  cphis = _mm256_cos_ps(phis);
	                  ahh   = negate_ymm8r4(cphis);
	                  return (ahh);
	          }
	          
	          
	             /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-83
	         */
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vh_f9183_ymm8r4(const __m256 phis,
	                                     const __m256 thts) {
	                                     
	                   __m256 sphis,cthts;
	                   __m256 ahv;
	                  sphis = _mm256_sin_ps(phis);
	                  cthts = _mm256_cos_ps(thti);
	                  ahv   = gms::math::negate_ymm8r4(_mm256_div_ps(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vh_f9183_ymm8r4_a(const float * __restrict pphis,
	                                       const float * __restrict pthts) {
	                         
	                   __m256 phis = _mm256_load_ps(&pphis[0]);
	                   __m256 thts = _mm256_load_ps(&pphis[0]);            
	                   __m256 sphis,cthts;
	                   __m256 ahv;
	                  sphis = _mm256_sin_ps(phis);
	                  cthts = _mm256_cos_ps(thti);
	                  ahv   = gms::math::negate_ymm8r4(_mm256_div_ps(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vh_f9183_ymm8r4_u(const float * __restrict  pphis,
	                                       const float * __restrict pthts) {
	                         
	                   __m256 phis = _mm256_loadu_ps(&pphis[0]);
	                   __m256 thts = _mm256_loadu_ps(&pphis[0]);            
	                   __m256 sphis,cthts;
	                   __m256 ahv;
	                  sphis = _mm256_sin_ps(phis);
	                  cthts = _mm256_cos_ps(thti);
	                  ahv   = gms::math::negate_ymm8r4(_mm256_div_ps(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-84 
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vv_f9184_ymm8r4(const __m256 thti) {
	           
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 sthti,cthti,num,den;
	                   __m256 avv;
	                  sthti = _mm256_sin_ps(thti);
	                  cthti = _mm256_cos_ps(thti);
	                  num   = _mm256_fmadd_ps(sthti,sthti,C10);
	                  den   = _mm256_mul_ps(cthti,cthti);
	                  avv   = _mm256_div_ps(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vv_f9184_ymm8r4_a(const float * __restrict pthti) {
	           
	                   __m256 thti = _mm256_load_ps(&pthti[0]);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 sthti,cthti,num,den;
	                   __m256 avv;
	                  sthti = _mm256_sin_ps(thti);
	                  cthti = _mm256_cos_ps(thti);
	                  num   = _mm256_fmadd_ps(sthti,sthti,C10);
	                  den   = _mm256_mul_ps(cthti,cthti);
	                  avv   = _mm256_div_ps(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vv_f9184_ymm8r4_u(const float * __restrict  pthti) {
	           
	                   __m256 thti = _mm256_loadu_ps(&pthti[0]);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 sthti,cthti,num,den;
	                   __m256 avv;
	                  sthti = _mm256_sin_ps(thti);
	                  cthti = _mm256_cos_ps(thti);
	                  num   = _mm256_fmadd_ps(sthti,sthti,C10);
	                  den   = _mm256_mul_ps(cthti,cthti);
	                  avv   = _mm256_div_ps(num,den);
	                  return (avv);
	           }
	           
	           
	           
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-85
	         */                  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_hh_f9185_ymm8r4() {
	           
	                  return _mm256_set1_ps(1.0f);
	           }  
	           
	           
	            /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-86
	         */   
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 a_vh_f9186_ymm8r4() {
	           
	                  return _mm256_setzero_ps();
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
	           __m256 rcs_vv_f9187_ymm8r4(const __m256 k0,
	                                       const __m256 h,
	                                       const __m256 l,
	                                       const __m256 thti) {
	                                       
	                  const __m256 C40 = _mm256_set1_ps(4.0f);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 k04,x0,x1,l2,h2,sthti;
	                   __m256 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_ps(k0,k0);
	                  l2    = _mm256_mul_ps(l,l);
	                  k04   = _mm256_mul_ps(x0,x0);
	                  h2    = _mm256_mul_ps(h,h);
	                  x1    = _mm256_sin_ps(thti);
	                  strm  = _mm256_fmadd_pd(x1,x1,C10);
	                  x0    = _mm256_mul_ps(C40,k04);
	                  sthti = _mm256_mul_ps(x1,x1);
	                  arg   = _mm256_mul_ps(_mm256_mul_ps(k0,k0),
	                                        _mm256_mul_ps(sthti,l2));
	                  trm   = _mm256_mul_ps(_mm256_mul_ps(l2,h2,x0));
	                  earg  = _mm256_exp_ps(arg);
	                  x1    = _mm256_mul_ps(strm,strm);
	                  inve  = _mm256_rcp14_ps(earg);
	                  rcs   = _mm256_mul_ps(trm,_mm256_mul_ps(x1,inve));
	                  return (rcs);
	         }
	           
	           
	           
	          __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_vv_f9187_ymm8r4_a(const float * __restrict pk0,
	                                       const float * __restrict ph,
	                                       const float * __restrict pl,
	                                       const float * __restrict pthti) {
	                             
	                   __m256 k0  = _mm256_load_ps(&pk0[0]);
	                   __m256 h   = _mm256_load_ps(&ph[0]);
	                   __m256 l   = _mm256_load_ps(&pl[0]);   
	                   __m256 thti= _mm256_load_ps(&pthti[0]);       
	                  const __m256 C40 = _mm256_set1_ps(4.0f);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 k04,x0,x1,l2,h2,sthti;
	                   __m256 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_ps(k0,k0);
	                  l2    = _mm256_mul_ps(l,l);
	                  k04   = _mm256_mul_ps(x0,x0);
	                  h2    = _mm256_mul_ps(h,h);
	                  x1    = _mm256_sin_ps(thti);
	                  strm  = _mm256_fmadd_pd(x1,x1,C10);
	                  x0    = _mm256_mul_ps(C40,k04);
	                  sthti = _mm256_mul_ps(x1,x1);
	                  arg   = _mm256_mul_ps(_mm256_mul_ps(k0,k0),
	                                        _mm256_mul_ps(sthti,l2));
	                  trm   = _mm256_mul_ps(_mm256_mul_ps(l2,h2,x0));
	                  earg  = _mm256_exp_ps(arg);
	                  x1    = _mm256_mul_ps(strm,strm);
	                  inve  = _mm256_rcp14_ps(earg);
	                  rcs   = _mm256_mul_ps(trm,_mm256_mul_ps(x1,inve));
	                  return (rcs);
	         }
	           
	            
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_vv_f9187_ymm8r4_u(const float * __restrict  pk0,
	                                         const float * __restrict  ph,
	                                         const float * __restrict  pl,
	                                         const float * __restrict  pthti) {
	                             
	                   __m256 k0  = _mm256_loadu_ps(&pk0[0]);
	                   __m256 h   = _mm256_loadu_ps(&ph[0]);
	                   __m256 l   = _mm256_loadu_ps(&pl[0]);   
	                   __m256 thti= _mm256_loadu_ps(&pthti[0]);       
	                  const __m256 C40 = _mm256_set1_ps(4.0f);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 k04,x0,x1,l2,h2,sthti;
	                   __m256 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_ps(k0,k0);
	                  l2    = _mm256_mul_ps(l,l);
	                  k04   = _mm256_mul_ps(x0,x0);
	                  h2    = _mm256_mul_ps(h,h);
	                  x1    = _mm256_sin_ps(thti);
	                  strm  = _mm256_fmadd_pd(x1,x1,C10);
	                  x0    = _mm256_mul_ps(C40,k04);
	                  sthti = _mm256_mul_ps(x1,x1);
	                  arg   = _mm256_mul_ps(_mm256_mul_ps(k0,k0),
	                                        _mm256_mul_ps(sthti,l2));
	                  trm   = _mm256_mul_ps(_mm256_mul_ps(l2,h2,x0));
	                  earg  = _mm256_exp_ps(arg);
	                  x1    = _mm256_mul_ps(strm,strm);
	                  inve  = _mm256_rcp14_ps(earg);
	                  rcs   = _mm256_mul_ps(trm,_mm256_mul_ps(x1,inve));
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
	           __m256 rcs_hh_f9188_ymm8r4(const __m256 k0,
	                                       const __m256 h,
	                                       const __m256 l,
	                                       const __m256 thti) {
	                                       
	                  const __m256 C40 = _mm256_set1_ps(4.0f);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 k04,x0,x1,l2,h2,sthti,cthti;
	                   __m256 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_ps(k0,k0);
	                  l2    = _mm256_mul_ps(l,l);
	                  strm  = _mm256_cos_ps(thti);
	                  k04   = _mm256_mul_ps(x0,x0);
	                  h2    = _mm256_mul_ps(h,h);
	                  x1    = _mm256_sin_ps(thti);
	                  x0    = _mm256_mul_ps(C40,k04);
	                  sthti = _mm256_mul_ps(x1,x1);
	                  arg   = _mm256_mul_ps(_mm256_mul_ps(k0,k0),
	                                        _mm256_mul_ps(sthti,l2));
	                  trm   = _mm256_mul_ps(_mm256_mul_ps(l2,h2,x0));
	                  x1    = _mm256_mul_ps(strm,strm);
	                  earg  = _mm256_exp_ps(arg);
	                  cthti = _mm256_mul_ps(x1,x1);
	                  inve  = _mm256_rcp14_ps(earg);
	                  rcs   = _mm256_mul_ps(trm,_mm256_mul_ps(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_hh_f9188_ymm8r4_a(const float * __restrict pk0,
	                                       const float * __restrict ph,
	                                       const float * __restrict pl,
	                                       const float * __restrict pthti) {
	                             
	                   __m256 k0  = _mm256_load_ps(&pk0[0]);
	                   __m256 h   = _mm256_load_ps(&ph[0]);
	                   __m256 l   = _mm256_load_ps(&pl[0]);   
	                   __m256 thti= _mm256_load_ps(&pthti[0]); 
	                                       
	                  const __m256 C40 = _mm256_set1_ps(4.0f);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 k04,x0,x1,l2,h2,sthti,cthti;
	                   __m256 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_ps(k0,k0);
	                  l2    = _mm256_mul_ps(l,l);
	                  strm  = _mm256_cos_ps(thti);
	                  k04   = _mm256_mul_ps(x0,x0);
	                  h2    = _mm256_mul_ps(h,h);
	                  x1    = _mm256_sin_ps(thti);
	                  x0    = _mm256_mul_ps(C40,k04);
	                  sthti = _mm256_mul_ps(x1,x1);
	                  arg   = _mm256_mul_ps(_mm256_mul_ps(k0,k0),
	                                        _mm256_mul_ps(sthti,l2));
	                  trm   = _mm256_mul_ps(_mm256_mul_ps(l2,h2,x0));
	                  x1    = _mm256_mul_ps(strm,strm);
	                  earg  = _mm256_exp_ps(arg);
	                  cthti = _mm256_mul_ps(x1,x1);
	                  inve  = _mm256_rcp14_ps(earg);
	                  rcs   = _mm256_mul_ps(trm,_mm256_mul_ps(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	            __ATTR_ALWAYS_INLINE__
	           
	          
                   
	           static inline
	           __m256 rcs_hh_f9188_ymm8r4_u(const float * __restrict  pk0,
	                                       const float * __restrict  ph,
	                                       const float * __restrict  pl,
	                                       const float * __restrict  pthti) {
	                             
	                   __m256 k0  = _mm256_loadu_ps(&pk0[0]);
	                   __m256 h   = _mm256_loadu_ps(&ph[0]);
	                   __m256 l   = _mm256_loadu_ps(&pl[0]);   
	                   __m256 thti= _mm256_loadu_ps(&pthti[0]); 
	                                       
	                  const __m256 C40 = _mm256_set1_ps(4.0f);
	                  const __m256 C10 = _mm256_set1_ps(1.0f);
	                   __m256 k04,x0,x1,l2,h2,sthti,cthti;
	                   __m256 rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm256_mul_ps(k0,k0);
	                  l2    = _mm256_mul_ps(l,l);
	                  strm  = _mm256_cos_ps(thti);
	                  k04   = _mm256_mul_ps(x0,x0);
	                  h2    = _mm256_mul_ps(h,h);
	                  x1    = _mm256_sin_ps(thti);
	                  x0    = _mm256_mul_ps(C40,k04);
	                  sthti = _mm256_mul_ps(x1,x1);
	                  arg   = _mm256_mul_ps(_mm256_mul_ps(k0,k0),
	                                        _mm256_mul_ps(sthti,l2));
	                  trm   = _mm256_mul_ps(_mm256_mul_ps(l2,h2,x0));
	                  x1    = _mm256_mul_ps(strm,strm);
	                  earg  = _mm256_exp_ps(arg);
	                  cthti = _mm256_mul_ps(x1,x1);
	                  inve  = _mm256_rcp14_ps(earg);
	                  rcs   = _mm256_mul_ps(trm,_mm256_mul_ps(cthti,inve));
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
	           __m256 rcs_vhhv_f9189_ymm8r4() {
	           
	                  return (_mm256_setzero_ps()); 
	          }
	          
	          
	          
	           
	      
	                                    
	        
                 
                 
                 
               
               
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_YMM8R4_H__*/
