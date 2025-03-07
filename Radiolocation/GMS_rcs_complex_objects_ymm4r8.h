

#ifndef __GMS_RCS_COMPLEX_OBJECTS_YMM4R8_H__
#define __GMS_RCS_COMPLEX_OBJECTS_YMM4R8_H__ 260920240740

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
              
              
                  __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8121_ymm4r8(const __m256d gam,
                                             const __m256d phi,
                                             const __m256d k0,
                                             __m256d * __restrict D1r,
                                             __m256d * __restrict D1i,
                                            __m256d * __restrict D2r,
                                            __m256d * __restrict D2i); 
                
                
                
                    __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8121_ymm4r8_a(const double * __restrict pgam,
                                               const double * __restrict pphi,
                                               const double * __restrict pk0,
                                               double * __restrict D1r,
                                               double * __restrict D1i,
                                               double * __restrict D2r,
                                               double * __restrict D2i); 
                
                
                   __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8121_ymm4r8_u(const double * __restrict  pgam,
                                               const double * __restrict  pphi,
                                               const double * __restrict  pk0,
                                               double * __restrict  D1r,
                                               double * __restrict  D1i,
                                               double * __restrict  D2r,
                                               double * __restrict  D2i); 
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter singly diffracted far-zone fields (E,H).
                    Formula: 8.1-19, 8.1-20
                
                */
                
                
                    __ATTR_HOT__
                  __ATTR_VECTORCALL__
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
                                                 __m256d * __restrict Hsi); 
                                                 
                                                 
                                                 
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__                         
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
                                                 double * __restrict Hsi); 
            
            
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
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
                                                 double * __restrict  Hsi); 
                                                 
            
            
            /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Ray normal-incidence to one of edge faces.
                    Formula: 8.1-24
            */
            
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8124_ymm4r8(const __m256d k0,
                                               const __m256d gam,
                                               __m256d * __restrict D1r,
                                               __m256d * __restrict D1i,
                                               __m256d * __restrict D2r,
                                               __m256d * __restrict D2i); 
                
                
                    __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8124_ymm4r8_a(const  double * __restrict pk0,
                                                 const  double * __restrict pgam,
                                                 double * __restrict D1r,
                                                 double * __restrict D1i,
                                                 double * __restrict D2r,
                                                 double * __restrict D2i); 
                
                
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8124_ymm4r8_u(const  double * __restrict  pk0,
                                                 const  double * __restrict pgam,
                                                 double * __restrict  D1r,
                                                 double * __restrict  D1i,
                                                 double * __restrict  D2r,
                                                 double * __restrict  D2i); 
                
                
                /*
                    Surface discontinuities.
                    General perfectly conducting convex edge.
                    Backscatter diffraction coefficient 'D'.
                    Backscatter direction axial caustic (for slightly diffracted rays).
                    Formula: 8.1-26
                */
                
                
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_Ddiff_f8126_ymm4r8(const __m256d gam,
                                             const __m256d phi,
                                             const __m256d k0,
                                             __m256d * __restrict Dr,
                                             __m256d * __restrict Di); 
                                             
                
                    
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_Ddiff_f8126_ymm4r8_a(const double * __restrict pgam,
                                               const double * __restrict pphi,
                                               const double * __restrict pk0,
                                               double * __restrict  Dr,
                                               double * __restrict  Di); 
                
                
                
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_Ddiff_f8126_ymm4r8_u(const double * __restrict  pgam,
                                                   const double * __restrict pphi,
                                                   const double * __restrict  pk0,
                                                   double * __restrict  Dr,
                                                   double * __restrict  Di); 
                
                
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
                   void EsHs_f8125_ymm4r8(const __m256d a,
                                           const __m256d k0,
                                           const __m256d r,
                                           const __m256d gam,
                                           const __m256d phi,
                                           const __m256d psi,
                                           __m256d * __restrict Esr,
                                           __m256d * __restrict Esi,
                                           __m256d * __restrict Hsr,
                                           __m256d * __restrict Hsi); 
                
                
                    
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void EsHs_f8125_ymm4r8_a(const double * __restrict pa,
                                             const double * __restrict pk0,
                                             const double * __restrict pr,
                                             const double * __restrict pgam,
                                             const double * __restrict pphi,
                                             const double * __restrict ppsi,
                                             double * __restrict Esr,
                                             double * __restrict Esi,
                                             double * __restrict Hsr,
                                             double * __restrict Hsi); 
                
            
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void EsHs_f8125_ymm4r8_u(const double * __restrict  pa,
                                             const double * __restrict  pk0,
                                             const double * __restrict  pr,
                                             const double * __restrict pgam,
                                             const double * __restrict  pphi,
                                             const double * __restrict  ppsi,
                                             double * __restrict  Esr,
                                             double * __restrict  Esi,
                                             double * __restrict  Hsr,
                                             double * __restrict  Hsi); 
                
                
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
                   void coef_D12_f8127_ymm4r8(const __m256d k0,
                                               const __m256d gam,
                                               const __m256d phi1,
                                               const __m256d phi2,
                                               __m256d * __restrict D1r,
                                               __m256d * __restrict D1i,
                                               __m256d * __restrict D2r,
                                               __m256d * __restrict D2i); 
                
                
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8127_ymm4r8_a(const double * __restrict pk0,
                                                 const double * __restrict pgam,
                                                 const double * __restrict pphi1,
                                                 const double * __restrict pphi2,
                                                 double * __restrict  D1r,
                                                 double * __restrict  D1i,
                                                 double * __restrict  D2r,
                                                 double * __restrict  D2i); 
                
                
                
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   void coef_D12_f8127_ymm4r8_u(const double * __restrict  pk0,
                                                 const double * __restrict  pgam,
                                                 const double * __restrict  pphi1,
                                                 const double * __restrict  pphi2,
                                                 double * __restrict   D1r,
                                                 double * __restrict  D1i,
                                                 double * __restrict   D2r,
                                                 double * __restrict  D2i); 
                
                
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
                   double rcs_f8162_ymm4r8_u(const double * __restrict pdAdl,
                                             const double *  __restrict pdl,
                                             const double   k0,
                                             const double   l); 
                  
                  
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   double rcs_f8162_ymm4r8_a(const double * __restrict pdAdl,
                                             const double * __restrict pdl,
                                             const double   k0,
                                             const double   l); 
                  
                  
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
                   double rcs_f8162_ymm4r8_avint_u(const double * __restrict pdAdl,
                                                   const double *  __restrict pdl,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri); 
                  
                  
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   double rcs_f8162_ymm4r8_avint_a(const double * __restrict pdAdl,
                                                   const double * __restrict pdl,
                                                   const double   k0,
                                                   const double   l,
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
                                             const int32_t NTAB); 
               
               
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
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
                   double rcs_f8162_ymm4r8_avint_u(const double * __restrict  pdAdl,
                                                   const double * __restrict  pdl,
                                                   double * __restrict  intr,
                                                   double * __restrict  inti,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB); 
               
               
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   double rcs_f8162_ymm4r8_avint_a(const double * __restrict pdAdl,
                                                   const double * __restrict pdl,
                                                   double * __restrict intr,
                                                   double * __restrict inti,
                                                   const double   k0,
                                                   const double   l,
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
                   double rcs_f8162_ymm4r8_cspint2t_u(const double * __restrict  pdAdl,
                                                     const double * __restrict  pdl,
                                                     double * __restrict  intr,
                                                     double * __restrict  inti,
                                                     struct RCS_F8162_DATA w,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB); 
               
               
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   double rcs_f8162_ymm4r8_cspint2t_a(const double * __restrict pdAdl,
                                                     const double * __restrict pdl,
                                                     double * __restrict intr,
                                                     double * __restrict inti,
                                                     struct RCS_F8162_DATA w,
                                                     const double   k0,
                                                     const double   l,
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
                   double rcs_f8162_ymm4r8_avint2t_u(const double * __restrict  pdAdl,
                                                     const double * __restrict  pdl,
                                                     double * __restrict  intr,
                                                     double * __restrict  inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB); 
               
               
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   double rcs_f8162_ymm4r8_avint2t_a(const double * __restrict pdAdl,
                                                     const double * __restrict pdl,
                                                     double * __restrict intr,
                                                     double * __restrict inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB); 
               
               
               
               
               
               /*
                     High frequency approximations.
                     Rounded-tip cone total nose-on
                     backscatter RCS.
                     Formula 8.1-93
               */
               
               
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   __m256d rcs_f8193_ymm4r8(const __m256d b,
                                            const __m256d a,
                                            const __m256d k0,
                                            const __m256d alp,
                                            const __m256d l); 
                 
                 
                 
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   __m256d rcs_f8193_ymm4r8_a(const double * __restrict pb,
                                              const double * __restrict pa,
                                              const double * __restrict pk0,
                                              const double * __restrict palp,
                                              const double * __restrict pl); 
                 
                 
                 
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   __m256d rcs_f8193_ymm4r8_u(const double * __restrict  pb,
                                              const double * __restrict  pa,
                                              const double * __restrict  pk0,
                                              const double * __restrict  palp,
                                              const double * __restrict _pl);
                 
                 
                 /*
                     High frequency approximations.
                     Backscatter RCS of conical frustum
                     for |theta| = PI/2-alpha
                     Formula 8.1-96
                 */
                 
                 
                    __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   __m256d rcs_f8196_ymm4r8(const __m256d k0,
                                            const __m256d alp,
                                            const __m256d a,
                                            const __m256d b); 
                 
                 
                     __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   __m256d rcs_f8196_ymm4r8_a(const double * __restrict pk0,
                                              const double * __restrict palp,
                                              const double * __restrict pa,
                                              const double * __restrict pb); 
                 
                 
                    __ATTR_HOT__
                  __ATTR_VECTORCALL__
                   __m256d rcs_f8196_ymm4r8_u(const double * __restrict pk0,
                                              const double * __restrict palp,
                                              const double * __restrict pa,
                                              const double * __restrict pb); 
                 
                 
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
	        
	        
	        
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d coef_Bg_f9137_ymm4r8(const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol); 
	       
	       
	            
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d coef_Bg_f9137_ymm4r8_a(const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol); 
	       
	       
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d coef_Bg_f9137_ymm4r8_u(const double * __restrict  pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
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
	           __m256d rcs_hh_f9133_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol); 
	       
	       
	       
	            __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9133_ymm4r8_a( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol); 
	       
	       
	       
	            __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9133_ymm4r8_u( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
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
	           __m256d rcs_vh_f9134_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol); 
	         
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vh_f9134_ymm4r8_a( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol); 
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vh_f9134_ymm4r8_u( const double * __restrict  pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
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
	           __m256d rcs_hv_f9135_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol); 
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hv_f9135_ymm4r8_a(  const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol); 
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d  rcs_hv_f9135_ymm4r8_u(const double * __restrict pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
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
	           __m256d rcs_vv_f9136_ymm4r8( const __m256d A,
	                                        const __m256d N,
	                                        const __m256d k0,
	                                        const __m256d epsr,
	                                        const __m256d epsi,
	                                        const __m256d thti,
	                                        const __m256d thts,
	                                        const __m256d phis,
	                                        const int pol); 
	         
	         
	         
	              __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9136_ymm4r8_a( const double * __restrict pA,
	                                          const double * __restrict pN,
	                                          const double * __restrict pk0,
	                                          const double * __restrict pepsr,
	                                          const double * __restrict pepsi,
	                                          const double * __restrict pthti,
	                                          const double * __restrict pthts,
	                                          const double * __restrict pphis,
	                                          const int pol ); 
	         
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9136_ymm4r8_u( const double * __restrict  pA,
	                                          const double * __restrict  pN,
	                                          const double * __restrict  pk0,
	                                          const double * __restrict  pepsr,
	                                          const double * __restrict  pepsi,
	                                          const double * __restrict  pthti,
	                                          const double * __restrict  pthts,
	                                          const double * __restrict  pphis,
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
	           __m256d rcs_vv_f9174_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui); 
	       
	       
	            __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9174_ymm4r8_a(const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui);
	       
	       
	       
	              __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9174_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui);
	       
	       
	        /*
	            Gaussian surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hh) polarization.
	            Formula 9.1-75
	        */
	        
	        
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9175_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui); 
	       
	       
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9175_ymm4r8_a(  const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui);
	       
	       
	              __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9175_ymm4r8_u(  const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui); 
	       
	       
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
	        
	        
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9177_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui); 
	       
	       
	       
	           __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9177_ymm4r8_a(const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui); 
	       
	       
	            __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9177_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui); 
	       
	       
	         /*
	            Exponential surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (hh) polarization.
	            Formula 9.1-78
	        */
	        
	        
	        
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9178_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti,
	                                       const __m256d epsr,
	                                       const __m256d epsi,
	                                       const __m256d mur,
	                                       const __m256d mui); 
	                                       
	       
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9178_ymm4r8_a(const double * __restrict pk0,
	                                         const double * __restrict ph,
	                                         const double * __restrict pl,
	                                         const double * __restrict pthti,
	                                         const double * __restrict pepsr,
	                                         const double * __restrict pepsi,
	                                         const double * __restrict pmur,
	                                         const double * __restrict pmui); 
	       
	       
	       
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9178_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti,
	                                         const double * __restrict  pepsr,
	                                         const double * __restrict  pepsi,
	                                         const double * __restrict  pmur,
	                                         const double * __restrict  pmui); 
	       
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
	         
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9187_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti); 
	           
	           
	           
	            __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9187_ymm4r8_a(const double * __restrict pk0,
	                                       const double * __restrict ph,
	                                       const double * __restrict pl,
	                                       const double * __restrict pthti); 
	           
	            
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_vv_f9187_ymm4r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti); 
	         
	         
	          /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Average backscattering RCS per unit area.
	              Gaussian surface height correlation coefficient.
	              Formula: 9.1-88
	         */
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9188_ymm4r8(const __m256d k0,
	                                       const __m256d h,
	                                       const __m256d l,
	                                       const __m256d thti); 
	         
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9188_ymm4r8_a(const double * __restrict pk0,
	                                       const double * __restrict ph,
	                                       const double * __restrict pl,
	                                       const double * __restrict pthti); 
	         
	             __ATTR_HOT__
                  __ATTR_VECTORCALL__
	           __m256d rcs_hh_f9188_ymm4r8_u(const double * __restrict  pk0,
	                                       const double * __restrict  ph,
	                                       const double * __restrict  pl,
	                                       const double * __restrict  pthti); 
	         
	         
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




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_YMM4R8_H__*/
