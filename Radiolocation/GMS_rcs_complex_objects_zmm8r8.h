

#ifndef __GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_H__
#define __GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_H__

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

    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_FULLVER =
      1000U*GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_MAJOR+
      100U*GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_MINOR+
      10U*GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_MICRO;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_CREATION_DATE = "11-05-2023 10:53 PM +00200 (THR 11 MAY 2023 GMT+2)";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_DESCRIPTION   = "AVX512 optimized Complex Surfaces Radar Cross Section functionality.";

}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"
#include "GMS_complex_zmm8r8.hpp"
#include "GMS_simd_utils.hpp"




namespace  gms {


         namespace radiolocation {
         
         
               /*
                   Work (input) arrays for kernel rcs_f8162_zmm8r8_2t_u and
                   rcs_f8162_zmm8r8_2t_a.
               */
               __ATTR_ALIGN__(64) struct RCS_F8162_DATA {
               
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void coef_D12_f8121_zmm8r8(const __m512d gam,
                                             const __m512d phi,
                                             const __m512d k0,
                                             __m512d * __restrict D1r,
                                             __m512d * __restrict D1i,
                                            __m512d * __restrict D2r,
                                            __m512d * __restrict D2i); 
                
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	        
                   void coef_D12_f8121_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgam,
                                               const double * __restrict __ATTR_ALIGN__(64) pphi,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0,
                                               double * __restrict __ATTR_ALIGN__(64) D1r,
                                               double * __restrict __ATTR_ALIGN__(64) D1i,
                                               double * __restrict __ATTR_ALIGN__(64) D2r,
                                               double * __restrict __ATTR_ALIGN__(64) D2i); 
                
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void coef_D12_f8121_zmm8r8_u(const double * __restrict  pgam,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void EsHs_f811920_zmm8r8(    const __m512d betai,
                                                 const __m512d betas,
                                                 const __m512d gam,
                                                 const __m512d phi,
                                                 const __m512d k0,
                                                 const __m512d r,
                                                 const __m512d rho,
                                                 const __m512d psi,
                                                 __m512d * __restrict Esr,
                                                 __m512d * __restrict Esi,
                                                 __m512d * __restrict Hsr,
                                                 __m512d * __restrict Hsi); 
                                                 
            
            
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void EsHs_f811920_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pbetai,
                                                 const double * __restrict __ATTR_ALIGN__(64) pbetas,
                                                 const double * __restrict __ATTR_ALIGN__(64) pgam,
                                                 const double * __restrict __ATTR_ALIGN__(64) pphi,
                                                 const double * __restrict __ATTR_ALIGN__(64) pk0,
                                                 const double * __restrict __ATTR_ALIGN__(64) pr,
                                                 const double * __restrict __ATTR_ALIGN__(64) prho,
                                                 const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                                 double * __restrict __ATTR_ALIGN__(64) Esr,
                                                 double * __restrict __ATTR_ALIGN__(64) Esi,
                                                 double * __restrict __ATTR_ALIGN__(64) Hsr,
                                                 double * __restrict __ATTR_ALIGN__(64) Hsi); 
            
            
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void EsHs_f811920_zmm8r8_u(  const double * __restrict  pbetai,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   void coef_D12_f8124_zmm8r8(const __m512d k0,
                                               const __m512d gam,
                                               __m512d * __restrict D1r,
                                               __m512d * __restrict D1i,
                                               __m512d * __restrict D2r,
                                               __m512d * __restrict D2i); 
                                               
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void coef_D12_f8124_zmm8r8_a(const  double * __restrict __ATTR_ALIGN__(64) pk0,
                                                 const  double * __restrict __ATTR_ALIGN__(64) pgam,
                                                 double * __restrict __ATTR_ALIGN__(64) D1r,
                                                 double * __restrict __ATTR_ALIGN__(64) D1i,
                                                 double * __restrict __ATTR_ALIGN__(64) D2r,
                                                 double * __restrict __ATTR_ALIGN__(64) D2i); 
                                                 
                
                
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void coef_D12_f8124_zmm8r8_u(const  double * __restrict  pk0,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void coef_Ddiff_f8126_zmm8r8(const __m512d gam,
                                             const __m512d phi,
                                             const __m512d k0,
                                             __m512d * __restrict Dr,
                                             __m512d * __restrict Di); 
                                             
                
                
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void coef_Ddiff_f8126_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pgam,
                                               const double * __restrict __ATTR_ALIGN__(64) pphi,
                                               const double * __restrict __ATTR_ALIGN__(64) pk0,
                                               double * __restrict __ATTR_ALIGN__(64)  Dr,
                                               double * __restrict __ATTR_ALIGN__(64)  Di); 
                                               
                
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   void coef_Ddiff_f8126_zmm8r8_u(const double * __restrict  pgam,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void EsHs_f8125_zmm8r8(const __m512d a,
                                           const __m512d k0,
                                           const __m512d r,
                                           const __m512d gam,
                                           const __m512d phi,
                                           const __m512d psi,
                                           __m512d * __restrict Esr,
                                           __m512d * __restrict Esi,
                                           __m512d * __restrict Hsr,
                                           __m512d * __restrict Hsi); 
                                           
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void EsHs_f8125_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pa,
                                             const double * __restrict __ATTR_ALIGN__(64) pk0,
                                             const double * __restrict __ATTR_ALIGN__(64) pr,
                                             const double * __restrict __ATTR_ALIGN__(64) pgam,
                                             const double * __restrict __ATTR_ALIGN__(64) pphi,
                                             const double * __restrict __ATTR_ALIGN__(64) ppsi,
                                             double * __restrict __ATTR_ALIGN__(64) Esr,
                                             double * __restrict __ATTR_ALIGN__(64) Esi,
                                             double * __restrict __ATTR_ALIGN__(64) Hsr,
                                             double * __restrict __ATTR_ALIGN__(64) Hsi); 
                                             
                
            
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	        
                   void EsHs_f8125_zmm8r8_u(const double * __restrict  pa,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void coef_D12_f8127_zmm8r8(const __m512d k0,
                                               const __m512d gam,
                                               const __m512d phi1,
                                               const __m512d phi2,
                                               __m512d * __restrict D1r,
                                               __m512d * __restrict D1i,
                                               __m512d * __restrict D2r,
                                               __m512d * __restrict D2i); 
                                               
                
                
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   void coef_D12_f8127_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                                 const double * __restrict __ATTR_ALIGN__(64) pgam,
                                                 const double * __restrict __ATTR_ALIGN__(64) pphi1,
                                                 const double * __restrict __ATTR_ALIGN__(64) pphi2,
                                                 double * __restrict __ATTR_ALIGN__(64)  D1r,
                                                 double * __restrict __ATTR_ALIGN__(64)  D1i,
                                                 double * __restrict __ATTR_ALIGN__(64)  D2r,
                                                 double * __restrict __ATTR_ALIGN__(64)  D2i); 
                                                 
                
                
                
               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   void coef_D12_f8127_zmm8r8_u(const double * __restrict  pk0,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	        
                   double rcs_f8162_zmm8r8_u(const double * __restrict pdAdl,
                                             const double *  __restrict pdl,
                                             const double   k0,
                                             const double   l); 
                                             
                  
                  
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double rcs_f8162_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pdAdl,
                                             const double * __restrict __ATTR_ALIGN__(64) pdl,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double rcs_f8162_zmm8r8_avint_u(const double * __restrict pdAdl,
                                                   const double *  __restrict pdl,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri); 
                                                   
                  
                  
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double rcs_f8162_zmm8r8_avint_a(const double * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                   const double * __restrict __ATTR_ALIGN__(64) pdl,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double rcs_f8162_zmm8r8_u(const double * __restrict  pdAdl,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double rcs_f8162_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pdAdl,
                                             const double * __restrict __ATTR_ALIGN__(64) pdl,
                                             double * __restrict __ATTR_ALIGN__(64) intr,
                                             double * __restrict __ATTR_ALIGN__(64) inti,
                                             double * __restrict __ATTR_ALIGN__(64) Y1,
                                             double * __restrict __ATTR_ALIGN__(64) Y2,
                                             double * __restrict __ATTR_ALIGN__(64) Y3,
                                             double * __restrict __ATTR_ALIGN__(64) E,
                                             double * __restrict __ATTR_ALIGN__(64) WRK
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double rcs_f8162_zmm8r8_avint_u(const double * __restrict  pdAdl,
                                                   const double * __restrict  pdl,
                                                   double * __restrict  intr,
                                                   double * __restrict  inti,
                                                   const double   k0,
                                                   const double   l,
                                                   int32_t & ierr,
                                                   int32_t & ieri,
                                                   const int32_t NTAB); 
                                                   
               
               
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double rcs_f8162_zmm8r8_avint_a(const double * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                   const double * __restrict __ATTR_ALIGN__(64) pdl,
                                                   double * __restrict __ATTR_ALIGN__(64) intr,
                                                   double * __restrict __ATTR_ALIGN__(64) inti,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
                   double rcs_f8162_zmm8r8_cspint2t_u(const double * __restrict  pdAdl,
                                                     const double * __restrict  pdl,
                                                     double * __restrict  intr,
                                                     double * __restrict  inti,
                                                     struct RCS_F8162_DATA w,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB); 
               
               
               
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double rcs_f8162_zmm8r8_cspint2t_a(const double * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                     const double * __restrict __ATTR_ALIGN__(64) pdl,
                                                     double * __restrict __ATTR_ALIGN__(64) intr,
                                                     double * __restrict __ATTR_ALIGN__(64) inti,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   double rcs_f8162_zmm8r8_avint2t_u(const double * __restrict  pdAdl,
                                                     const double * __restrict  pdl,
                                                     double * __restrict  intr,
                                                     double * __restrict  inti,
                                                     int32_t & ierr,
                                                     int32_t & ieri,
                                                     const double   k0,
                                                     const double   l,
                                                     const int32_t NTAB); 
                                                     
               
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   double rcs_f8162_zmm8r8_avint2t_a(const double * __restrict __ATTR_ALIGN__(64) pdAdl,
                                                     const double * __restrict __ATTR_ALIGN__(64) pdl,
                                                     double * __restrict __ATTR_ALIGN__(64) intr,
                                                     double * __restrict __ATTR_ALIGN__(64) inti,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   __m512d rcs_f8193_zmm8r8(const __m512d b,
                                            const __m512d a,
                                            const __m512d k0,
                                            const __m512d alp,
                                            const __m512d l); 
                                            
                 
                 
                 
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
                   __m512d rcs_f8193_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pb,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) palp,
                                              const double * __restrict __ATTR_ALIGN__(64) pl); 
                                              
                 
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   __m512d rcs_f8193_zmm8r8_u(const double * __restrict  pb,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   __m512d rcs_f8196_zmm8r8(const __m512d k0,
                                            const __m512d alp,
                                            const __m512d a,
                                            const __m512d b); 
                                            
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   __m512d rcs_f8196_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
                                              const double * __restrict __ATTR_ALIGN__(64) palp,
                                              const double * __restrict __ATTR_ALIGN__(64) pa,
                                              const double * __restrict __ATTR_ALIGN__(64) pb); 
                                              
                 
                 
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
                   __m512d rcs_f8196_zmm8r8_u(const double * __restrict pk0,
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
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d rcs_perpendic_f8194_zmm8r8(const __m512d h,
	                                        const __m512d l,
	                                        const __m512d b,
	                                        const __m512d a,
	                                        const __m512d k0,
	                                        const __m512d tht,
	                                        const __m512d alp) {
	                                 
	                                  
	                 const __m512d C314159265358979323846264338328  = 
                                                     _mm512_set1_pd(3.14159265358979323846264338328f); 
                         const __m512d C1772453850905516027298167483341 = 
                                                     _mm512_set1_pd(1.772453850905516027298167483341f);
                         const __m512d C078539816339744830961566084582  = 
                                                     _mm512_set1_pd(0.78539816339744830961566084582f);
                         const __m512d C10                              = 
                                                     _mm512_set1_pd(1.0f);  
                         const __m512d C15                              = 
                                                     _mm512_set1_pd(1.5f); 
                         const __m512d C05                              = 
                                                     _mm512_set1_pd(0.5f);
                         const __m512d C20                              =
                                                     _mm512_set1_pd(2.0f);
                         register __m512d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m512d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m512d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m512d cpin1,cpin2,trm1,trm2,rcs;
                         __m512d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm512_sub_pd(h,_mm512_add_pd(l,b));
                         sint = xsin(tht);
                         k02  = _mm512_add_pd(k0,k0);
                         n    = _mm512_mul_pd(C15,_mm512_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm512_rcp14_pd(sint);
                         a2   = _mm512_mul_pd(a,C05);
                         ear  = _mm512_setzero_pd();
                         sk02 = _mm512_sqrt_pd(_mm512_mul_pd(k0,C05));
                         x0   = _mm512_mul_pd(hlb,_mm512_sub_pd(cost,b));
                         invn = _mm512_rcp14_pd(n);
                         //x2   = _mm512_mul_pd(a,C05);
                         eai  = _mm512_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm512_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm512_sqrt_pd(_mm512_mul_pd(a2,csct));
                         atant= _mm512_mul_pd(a,tant);
                         cost = xcos(tht);
                         x0   = _mm512_mul_pd(b,C1772453850905516027298167483341);
                         cexp_zmm8r8(ear,eai,&cer1,&cei1);
                         cer1 = _mm512_mul_pd(x0,cer1);
                         spin = xsin(pin);
                         cei1 = _mm512_mul_pd(x0,cei1);
                         cpin = xcos(pin);
                         x1   = _mm512_mul_pd(_mm512_sub_pd(h,atant),cost);
                         eai2 = _mm512_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear,eai2,&cer2,&cei);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_mul_pd(x0,sacs);
                         cer2 = _mm512_mul_pd(cer2,x1);
                         cei2 = _mm512_mul_pd(cei2,x1);
                         cpin1= _mm512_rcp14_pd(_mm512_sub_pd(cpin,C10));
                         x2   = _mm512_mul_pd(C20,_mm512_add_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x2,invn));
                         x3   = _mm512_rcp14_pd(_mm512_sub_pd(cpin,cpin2));
                         trm1 = _mm512_sub_pd(cpin1,x3);
                         cmul_zmm8r8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm512_mul_pd(t0r,trm1);
                         t0i  = _mm512_mul_pd(t0i,trm1);
                         x0   = _mm512_mul_pd(C20,_mm512_sub_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x0,invn));
                         x1   = _mm512_rcp14_pd(cpin2);
                         trm2 = _mm512_sub_pd(cpin1,x1);
                         x2   = _mm512_fmadd_pd(cost,_mm512_mul_pd(k02,
                                                               _mm512_add_pd(h,atant)));
                         eai3 = _mm512_add_pd(C078539816339744830961566084582,x2);
                         cexp_zmm8r8(ear,ea3,&cer3,&cei3);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_sqrt_pd(_mm512_mul_pd(gms::math::
                                                                  negate_zmm8r8(a2),csct));
                         x2   = _mm512_mul_pd(x0,x1);
                         cer3 = _mm512_mul_pd(_mm512_mul_pd(cer3,x2),trm2);
                         cei3 = _mm512_mul_pd(_mm512_mul_pd(cei3,x2),trm2);
                         cmul_zmm8r8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_zmm8r8(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	                                  
                    __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d rcs_perpendic_f8194_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) ph,
	                                         const double * __restrict __ATTR_ALIGN__(64) pl,
	                                         const double * __restrict __ATTR_ALIGN__(64) pb,
	                                         const double * __restrict __ATTR_ALIGN__(64) pa,
	                                         const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                         const double * __restrict __ATTR_ALIGN__(64) ptht,
	                                         const double * __restrict __ATTR_ALIGN__(64) palp) {
	                                 
	                  
	                 register __m512d h  = _mm512_load_pd(&ph[0]);
	                 register __m512d l  = _mm512_load_pd(&pl[0]); 
	                 register __m512d b  = _mm512_load_pd(&pb[0]);   
	                 register __m512d a  = _mm512_load_pd(&pa[0]);  
	                 register __m512d k0 = _mm512_load_pd(&pk0[0]);
	                 register __m512d tht= _mm512_load_pd(&ptht[0]); 
	                 register __m512d alp= _mm512_load_pd(&palp[0]);        
	                 const __m512d C314159265358979323846264338328  = 
                                                     _mm512_set1_pd(3.14159265358979323846264338328f); 
                         const __m512d C1772453850905516027298167483341 = 
                                                     _mm512_set1_pd(1.772453850905516027298167483341f);
                         const __m512d C078539816339744830961566084582  = 
                                                     _mm512_set1_pd(0.78539816339744830961566084582f);
                         const __m512d C10                              = 
                                                     _mm512_set1_pd(1.0f);  
                         const __m512d C15                              = 
                                                     _mm512_set1_pd(1.5f); 
                         const __m512d C05                              = 
                                                     _mm512_set1_pd(0.5f);
                         const __m512d C20                              =
                                                     _mm512_set1_pd(2.0f);
                         __m512d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m512d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m512d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m512d cpin1,cpin2,trm1,trm2,rcs;
                         __m512d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm512_sub_pd(h,_mm512_add_pd(l,b));
                         sint = xsin(tht);
                         k02  = _mm512_add_pd(k0,k0);
                         n    = _mm512_mul_pd(C15,_mm512_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm512_rcp14_pd(sint);
                         a2   = _mm512_mul_pd(a,C05);
                         ear  = _mm512_setzero_pd();
                         sk02 = _mm512_sqrt_pd(_mm512_mul_pd(k0,C05));
                         x0   = _mm512_mul_pd(hlb,_mm512_sub_pd(cost,b));
                         invn = _mm512_rcp14_pd(n);
                         //x2   = _mm512_mul_pd(a,C05);
                         eai  = _mm512_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm512_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm512_sqrt_pd(_mm512_mul_pd(a2,csct));
                         atant= _mm512_mul_pd(a,tant);
                         cost = xcos(tht);
                         x0   = _mm512_mul_pd(b,C1772453850905516027298167483341);
                         cexp_zmm8r8(ear,eai,&cer1,&cei1);
                         cer1 = _mm512_mul_pd(x0,cer1);
                         spin = xsin(pin);
                         cei1 = _mm512_mul_pd(x0,cei1);
                         cpin = xcos(pin);
                         x1   = _mm512_mul_pd(_mm512_sub_pd(h,atant),cost);
                         eai2 = _mm512_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear,eai2,&cer2,&cei);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_mul_pd(x0,sacs);
                         cer2 = _mm512_mul_pd(cer2,x1);
                         cei2 = _mm512_mul_pd(cei2,x1);
                         cpin1= _mm512_rcp14_pd(_mm512_sub_pd(cpin,C10));
                         x2   = _mm512_mul_pd(C20,_mm512_add_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x2,invn));
                         x3   = _mm512_rcp14_pd(_mm512_sub_pd(cpin,cpin2));
                         trm1 = _mm512_sub_pd(cpin1,x3);
                         cmul_zmm8r8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm512_mul_pd(t0r,trm1);
                         t0i  = _mm512_mul_pd(t0i,trm1);
                         x0   = _mm512_mul_pd(C20,_mm512_sub_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x0,invn));
                         x1   = _mm512_rcp14_pd(cpin2);
                         trm2 = _mm512_sub_pd(cpin1,x1);
                         x2   = _mm512_fmadd_pd(cost,_mm512_mul_pd(k02,
                                                               _mm512_add_pd(h,atant)));
                         eai3 = _mm512_add_pd(C078539816339744830961566084582,x2);
                         cexp_zmm8r8(ear,ea3,&cer3,&cei3);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_sqrt_pd(_mm512_mul_pd(gms::math::
                                                                  negate_zmm8r8(a2),csct));
                         x2   = _mm512_mul_pd(x0,x1);
                         cer3 = _mm512_mul_pd(_mm512_mul_pd(cer3,x2),trm2);
                         cei3 = _mm512_mul_pd(_mm512_mul_pd(cei3,x2),trm2);
                         cmul_zmm8r8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_zmm8r8(t1r,t1i);
                         return (rcs);
	        }
	        
	        
	        
	           __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d rcs_perpendic_f8194_zmm8r8_u(    const double * __restrict  ph,
	                                                    const double * __restrict  pl,
	                                                    const double * __restrict  pb,
	                                                    const double * __restrict  pa,
	                                                    const double * __restrict  pk0,
	                                                    const double * __restrict  ptht,
	                                                    const double * __restrict  palp) {
	                                 
	                  
	                 register __m512d h  = _mm512_loadu_pd(&ph[0]);
	                 register __m512d l  = _mm512_loadu_pd(&pl[0]); 
	                 register __m512d b  = _mm512_loadu_pd(&pb[0]);   
	                 register __m512d a  = _mm512_loadu_pd(&pa[0]);  
	                 register __m512d k0 = _mm512_loadu_pd(&pk0[0]);
	                 register __m512d tht= _mm512_loadu_pd(&ptht[0]); 
	                 register __m512d alp= _mm512_loadu_pd(&palp[0]);        
	                 const __m512d C314159265358979323846264338328  = 
                                                     _mm512_set1_pd(3.14159265358979323846264338328f); 
                         const __m512d C1772453850905516027298167483341 = 
                                                     _mm512_set1_pd(1.772453850905516027298167483341f);
                         const __m512d C078539816339744830961566084582  = 
                                                     _mm512_set1_pd(0.78539816339744830961566084582f);
                         const __m512d C10                              = 
                                                     _mm512_set1_pd(1.0f);  
                         const __m512d C15                              = 
                                                     _mm512_set1_pd(1.5f); 
                         const __m512d C05                              = 
                                                     _mm512_set1_pd(0.5f);
                         const __m512d C20                              =
                                                     _mm512_set1_pd(2.0f);
                         __m512d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m512d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m512d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m512d cpin1,cpin2,trm1,trm2,rcs;
                         __m512d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm512_sub_pd(h,_mm512_add_pd(l,b));
                         sint = xsin(tht);
                         k02  = _mm512_add_pd(k0,k0);
                         n    = _mm512_mul_pd(C15,_mm512_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm512_rcp14_pd(sint);
                         a2   = _mm512_mul_pd(a,C05);
                         ear  = _mm512_setzero_pd();
                         sk02 = _mm512_sqrt_pd(_mm512_mul_pd(k0,C05));
                         x0   = _mm512_mul_pd(hlb,_mm512_sub_pd(cost,b));
                         invn = _mm512_rcp14_pd(n);
                         //x2   = _mm512_mul_pd(a,C05);
                         eai  = _mm512_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm512_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm512_sqrt_pd(_mm512_mul_pd(a2,csct));
                         atant= _mm512_mul_pd(a,tant);
                         cost = xcos(tht);
                         x0   = _mm512_mul_pd(b,C1772453850905516027298167483341);
                         cexp_zmm8r8(ear,eai,&cer1,&cei1);
                         cer1 = _mm512_mul_pd(x0,cer1);
                         spin = xsin(pin);
                         cei1 = _mm512_mul_pd(x0,cei1);
                         cpin = xcos(pin);
                         x1   = _mm512_mul_pd(_mm512_sub_pd(h,atant),cost);
                         eai2 = _mm512_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear,eai2,&cer2,&cei);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_mul_pd(x0,sacs);
                         cer2 = _mm512_mul_pd(cer2,x1);
                         cei2 = _mm512_mul_pd(cei2,x1);
                         cpin1= _mm512_rcp14_pd(_mm512_sub_pd(cpin,C10));
                         x2   = _mm512_mul_pd(C20,_mm512_add_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x2,invn));
                         x3   = _mm512_rcp14_pd(_mm512_sub_pd(cpin,cpin2));
                         trm1 = _mm512_sub_pd(cpin1,x3);
                         cmul_zmm8r8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm512_mul_pd(t0r,trm1);
                         t0i  = _mm512_mul_pd(t0i,trm1);
                         x0   = _mm512_mul_pd(C20,_mm512_sub_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x0,invn));
                         x1   = _mm512_rcp14_pd(cpin2);
                         trm2 = _mm512_sub_pd(cpin1,x1);
                         x2   = _mm512_fmadd_pd(cost,_mm512_mul_pd(k02,
                                                               _mm512_add_pd(h,atant)));
                         eai3 = _mm512_add_pd(C078539816339744830961566084582,x2);
                         cexp_zmm8r8(ear,ea3,&cer3,&cei3);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_sqrt_pd(_mm512_mul_pd(gms::math::
                                                                  negate_zmm8r8(a2),csct));
                         x2   = _mm512_mul_pd(x0,x1);
                         cer3 = _mm512_mul_pd(_mm512_mul_pd(cer3,x2),trm2);
                         cei3 = _mm512_mul_pd(_mm512_mul_pd(cei3,x2),trm2);
                         cmul_zmm8r8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_zmm8r8(t1r,t1i);
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
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           __m512d rcs_parallel_f8194_zmm8r8(const __m512d h,
	                                             const __m512d l,
	                                             const __m512d b,
	                                             const __m512d a,
	                                             const __m512d k0,
	                                             const __m512d tht,
	                                             const __m512d alp) {
	                                 
	                                  
	                 const __m512d C314159265358979323846264338328  = 
                                                     _mm512_set1_pd(3.14159265358979323846264338328f); 
                         const __m512d C1772453850905516027298167483341 = 
                                                     _mm512_set1_pd(1.772453850905516027298167483341f);
                         const __m512d C078539816339744830961566084582  = 
                                                     _mm512_set1_pd(0.78539816339744830961566084582f);
                         const __m512d C10                              = 
                                                     _mm512_set1_pd(1.0f);  
                         const __m512d C15                              = 
                                                     _mm512_set1_pd(1.5f); 
                         const __m512d C05                              = 
                                                     _mm512_set1_pd(0.5f);
                         const __m512d C20                              =
                                                     _mm512_set1_pd(2.0f);
                         register __m512d pin,n,invn,spin,cos1,k02,cos2,sint,cost;
                         register __m512d ear,eai1,eai2,eai3,cer1,cei1,sk02,sacs;
                         register __m512d cer2,cei2,cer3,cei3,x0,x1,x2,x3,atant;
                         register __m512d cpin1,cpin2,trm1,trm2,rcs;
                         __m512d t0r,t0i,t1r,t1i,a2;
                         hlb  = _mm512_sub_pd(h,_mm512_add_pd(l,b));
                         sint = xsin(tht);
                         k02  = _mm512_add_pd(k0,k0);
                         n    = _mm512_mul_pd(C15,_mm512_div_pd(alp,  
                                                           C314159265358979323846264338328));   
                         csct = _mm512_rcp14_pd(sint);
                         a2   = _mm512_mul_pd(a,C05);
                         ear  = _mm512_setzero_pd();
                         sk02 = _mm512_sqrt_pd(_mm512_mul_pd(k0,C05));
                         x0   = _mm512_mul_pd(hlb,_mm512_sub_pd(cost,b));
                         invn = _mm512_rcp14_pd(n);
                         //x2   = _mm512_mul_pd(a,C05);
                         eai  = _mm512_mul_pd(k02,x0);
                         tant = xtanf(tht);
                         pin  = _mm512_mul_pd(C314159265358979323846264338328,invn);
                         sacs = _mm512_sqrt_pd(_mm512_mul_pd(a2,csct));
                         atant= _mm512_mul_pd(a,tant);
                         cost = xcos(tht);
                         x0   = _mm512_mul_pd(b,C1772453850905516027298167483341);
                         cexp_zmm8r8(ear,eai,&cer1,&cei1);
                         cer1 = _mm512_mul_pd(x0,cer1);
                         spin = xsin(pin);
                         cei1 = _mm512_mul_pd(x0,cei1);
                         cpin = xcos(pin);
                         x1   = _mm512_mul_pd(_mm512_sub_pd(h,atant),cost);
                         eai2 = _mm512_fmadd_pd(k02,x1,C078539816339744830961566084582);
                         cexp_zmm8r8(ear,eai2,&cer2,&cei);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_mul_pd(x0,sacs);
                         cer2 = _mm512_mul_pd(cer2,x1);
                         cei2 = _mm512_mul_pd(cei2,x1);
                         cpin1= _mm512_rcp14_pd(_mm512_sub_pd(cpin,C10));
                         x2   = _mm512_mul_pd(C20,_mm512_add_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x2,invn));
                         x3   = _mm512_rcp14_pd(_mm512_sub_pd(cpin,cpin2));
                         trm1 = _mm512_add_pd(cpin1,x3);
                         cmul_zmm8r8(cer1,cei1,cer2,cei2,&t0r,&t0i);
                         t0r  = _mm512_mul_pd(t0r,trm1);
                         t0i  = _mm512_mul_pd(t0i,trm1);
                         x0   = _mm512_mul_pd(C20,_mm512_sub_pd(alp,tht));
                         cpin2= xcos(_mm512_mul_pd(x0,invn));
                         x1   = _mm512_rcp14_pd(cpin2);
                         trm2 = _mm512_add_pd(cpin1,x1);
                         x2   = _mm512_fmadd_pd(cost,_mm512_mul_pd(k02,
                                                               _mm512_add_pd(h,atant)));
                         eai3 = _mm512_add_pd(C078539816339744830961566084582,x2);
                         cexp_zmm8r8(ear,ea3,&cer3,&cei3);
                         x0   = _mm512_div_pd(spin,_mm512_mul_pd(n,sk02));
                         x1   = _mm512_sqrt_pd(_mm512_mul_pd(gms::math::
                                                                  negate_zmm8r8(a2),csct));
                         x2   = _mm512_mul_pd(x0,x1);
                         cer3 = _mm512_mul_pd(_mm512_mul_pd(cer3,x2),trm2);
                         cei3 = _mm512_mul_pd(_mm512_mul_pd(cei3,x2),trm2);
                         cmul_zmm8r8(t0r,t0i,cer3,cei3,&t1r,&t1i);
                         rcs  = cabs_zmm8r8(t1r,t1i);
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d coef_Bg_f9137_zmm8r8(const __m512d A,
	                                        const __m512d N,
	                                        const __m512d k0,
	                                        const __m512d epsr,
	                                        const __m512d epsi,
	                                        const __m512d thti,
	                                        const __m512d thts,
	                                        const __m512d phis,
	                                        const int pol); 
	       
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d coef_Bg_f9137_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pA,
	                                          const double * __restrict __ATTR_ALIGN__(64) pN,
	                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthts,
	                                          const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                          const int pol); 
	                                          
	       
	       
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	       
	           __m512d coef_Bg_f9137_zmm8r8_u(const double * __restrict  pA,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_hh_f9133_zmm8r8( const __m512d A,
	                                        const __m512d N,
	                                        const __m512d k0,
	                                        const __m512d epsr,
	                                        const __m512d epsi,
	                                        const __m512d thti,
	                                        const __m512d thts,
	                                        const __m512d phis,
	                                        const int pol); 
	                                        
	       
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_hh_f9133_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pA,
	                                          const double * __restrict __ATTR_ALIGN__(64) pN,
	                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthts,
	                                          const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                          const int pol); 
	                                          
	       
	       
	       
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_hh_f9133_zmm8r8_u( const double * __restrict pA,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d rcs_vh_f9134_zmm8r8( const __m512d A,
	                                        const __m512d N,
	                                        const __m512d k0,
	                                        const __m512d epsr,
	                                        const __m512d epsi,
	                                        const __m512d thti,
	                                        const __m512d thts,
	                                        const __m512d phis,
	                                        const int pol); 
	                                        
	         
	         
	         
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d rcs_vh_f9134_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pA,
	                                          const double * __restrict __ATTR_ALIGN__(64) pN,
	                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthts,
	                                          const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                          const int pol); 
	                                          
	         
	         
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_vh_f9134_zmm8r8_u( const double * __restrict  pA,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d rcs_hv_f9135_zmm8r8( const __m512d A,
	                                        const __m512d N,
	                                        const __m512d k0,
	                                        const __m512d epsr,
	                                        const __m512d epsi,
	                                        const __m512d thti,
	                                        const __m512d thts,
	                                        const __m512d phis,
	                                        const int pol); 
	                                        
	         
	         
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d rcs_hv_f9135_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pA,
	                                          const double * __restrict __ATTR_ALIGN__(64) pN,
	                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthts,
	                                          const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                          const int pol); 
	                                          
	         
	         
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d  rcs_hv_f9135_zmm8r8_u(const double * __restrict pA,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_vv_f9136_zmm8r8( const __m512d A,
	                                        const __m512d N,
	                                        const __m512d k0,
	                                        const __m512d epsr,
	                                        const __m512d epsi,
	                                        const __m512d thti,
	                                        const __m512d thts,
	                                        const __m512d phis,
	                                        const int pol); 
	         
	         
	         
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_vv_f9136_zmm8r8_a( const double * __restrict __ATTR_ALIGN__(64) pA,
	                                          const double * __restrict __ATTR_ALIGN__(64) pN,
	                                          const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                          const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                          const double * __restrict __ATTR_ALIGN__(64) pthts,
	                                          const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                          const int pol ); 
	         
	         
	         
	            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_vv_f9136_zmm8r8_u( const double * __restrict  pA,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_vv_f9174_zmm8r8(const __m512d k0,
	                                       const __m512d h,
	                                       const __m512d l,
	                                       const __m512d thti,
	                                       const __m512d epsr,
	                                       const __m512d epsi,
	                                       const __m512d mur,
	                                       const __m512d mui); 
	                                       
	       
	       
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_vv_f9174_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                         const double * __restrict __ATTR_ALIGN__(64) ph,
	                                         const double * __restrict __ATTR_ALIGN__(64) pl,
	                                         const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmur,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmui); 
	       
	       
	       
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_vv_f9174_zmm8r8_u(const double * __restrict  pk0,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_hh_f9175_zmm8r8(const __m512d k0,
	                                       const __m512d h,
	                                       const __m512d l,
	                                       const __m512d thti,
	                                       const __m512d epsr,
	                                       const __m512d epsi,
	                                       const __m512d mur,
	                                       const __m512d mui); 
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_hh_f9175_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                         const double * __restrict __ATTR_ALIGN__(64) ph,
	                                         const double * __restrict __ATTR_ALIGN__(64) pl,
	                                         const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmur,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmui);
	                                         
	       
	       
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_hh_f9175_zmm8r8_u(  const double * __restrict  pk0,
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
	           __m512d rcs_hv_f9176_zmm8r8() { 
	           
	                return _mm512_setzero_pd();
	         } 
	         
	         
	         
	         /*
	            Exponential surface-height correlation
	            coefficient of average backscattering RCS 
	            per unit area.
	            RCS (vv) polarization.
	            Formula 9.1-77
	        */
	        
	        
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_vv_f9177_zmm8r8(const __m512d k0,
	                                       const __m512d h,
	                                       const __m512d l,
	                                       const __m512d thti,
	                                       const __m512d epsr,
	                                       const __m512d epsi,
	                                       const __m512d mur,
	                                       const __m512d mui); 
	       
	       
	       
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d rcs_vv_f9177_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                         const double * __restrict __ATTR_ALIGN__(64) ph,
	                                         const double * __restrict __ATTR_ALIGN__(64) pl,
	                                         const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmur,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmui); 
	       
	       
	       
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           
	           __m512d rcs_vv_f9177_zmm8r8_u(const double * __restrict  pk0,
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
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_hh_f9178_zmm8r8(const __m512d k0,
	                                       const __m512d h,
	                                       const __m512d l,
	                                       const __m512d thti,
	                                       const __m512d epsr,
	                                       const __m512d epsi,
	                                       const __m512d mur,
	                                       const __m512d mui); 
	       
	       
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	         
	           __m512d rcs_hh_f9178_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                         const double * __restrict __ATTR_ALIGN__(64) ph,
	                                         const double * __restrict __ATTR_ALIGN__(64) pl,
	                                         const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsr,
	                                         const double * __restrict __ATTR_ALIGN__(64) pepsi,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmur,
	                                         const double * __restrict __ATTR_ALIGN__(64) pmui); 
	                                         
	       
	       
	       
	           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	          
	           __m512d rcs_hh_f9178_zmm8r8_u(const double * __restrict  pk0,
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
	           __m512d rcs_hv_f9179_zmm8r8() { 
	           
	                return _mm512_setzero_pd();
	         } 
	         
	         
	         /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-80
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_vv_f9180_zmm8r8(const __m512d thti,
	                                         const __m512d thts,
	                                         const __m512d phis) {
	                                         
	                  register __m512d sthti,sthts,cthti,cthts,cphis;
	                  register __m512d num,den,avv;  
	                  sthti = xsin(thti);
	                  cthti = xcos(thti);
	                  cphis = xcos(phis);
	                  cthts = xcos(thts);
	                  sthts = xsin(thts);
	                  num   = _mm512_fmsub_pd(sthti,sthts,cphis);
	                  den   = _mm512_mul_pd(cthti,cthts);
	                  avv   = _mm512_div_pd(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_vv_f9180_zmm8r8_a(  const double * __restrict __ATTR_ALIGN__(64) pthti,
	                                         const double * __restrict __ATTR_ALIGN__(64) pthts,
	                                         const double * __restrict __ATTR_ALIGN__(64) pphis) {
	                     
	                  register __m512d thti = _mm512_load_pd(&pthti[0]);
	                  register __m512d thts = _mm512_load_pd(&pthts[0]);
	                  register __m512d phis = _mm512_load_pd(&phis[0]);                    
	                  register __m512d sthti,sthts,cthti,cthts,cphis;
	                  register __m512d num,den,avv;  
	                  sthti = xsin(thti);
	                  cthti = xcos(thti);
	                  cphis = xcos(phis);
	                  cthts = xcos(thts);
	                  sthts = xsin(thts);
	                  num   = _mm512_fmsub_pd(sthti,sthts,cphis);
	                  den   = _mm512_mul_pd(cthti,cthts);
	                  avv   = _mm512_div_pd(num,den);
	                  return (avv);        
	          }
	          
	          
	          
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_vv_f9180_zmm8r8_u(  const double * __restrict  pthti,
	                                         const double * __restrict  pthts,
	                                         const double * __restrict  pphis) {
	                     
	                  register __m512d thti = _mm512_loadu_pd(&pthti[0]);
	                  register __m512d thts = _mm512_loadu_pd(&pthts[0]);
	                  register __m512d phis = _mm512_loadu_pd(&phis[0]);                    
	                  register __m512d sthti,sthts,cthti,cthts,cphis;
	                  register __m512d num,den,avv;  
	                  sthti = xsin(thti);
	                  cthti = xcos(thti);
	                  cphis = xcos(phis);
	                  cthts = xcos(thts);
	                  sthts = xsin(thts);
	                  num   = _mm512_fmsub_pd(sthti,sthts,cphis);
	                  den   = _mm512_mul_pd(cthti,cthts);
	                  avv   = _mm512_div_pd(num,den);
	                  return (avv);        
	          }
	          
	          
	          /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-81
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_hv_f9181_zmm8r8(const __m512d phis,
	                                     const __m512d thti) {
	                                     
	                  register __m512d sphis,cthti;
	                  register __m512d ahv;
	                  sphis = xsin(phis);
	                  cthti = xcos(thti);
	                  ahv   = _mm512_div_pd(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_hv_f9181_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                     const double * __restrict __ATTR_ALIGN__(64) pthti) {
	                                     
	                  register __m512d thti = _mm512_load_pd(&pthti[0]);
	                  register __m512d phis = _mm512_load_pd(&phis[0]); 
	                  register __m512d sphis,cthti;
	                  register __m512d ahv;
	                  sphis = xsin(phis);
	                  cthti = xcos(thti);
	                  ahv   = _mm512_div_pd(sphis,cthti);
	                  return (ahv);                           
	         }
	          
	       
	        
	           __ATTR_ALWAYS_INLINE__
	           
	           static inline
	           __m512d a_hv_f9181_zmm8r8_u(const double * __restrict pphis,
	                                     const double * __restrict  pthti) {
	                                     
	                  register __m512d thti = _mm512_loadu_pd(&pthti[0]);
	                  register __m512d phis = _mm512_loadu_pd(&phis[0]); 
	                  register __m512d sphis,cthti;
	                  register __m512d ahv;
	                  sphis = xsin(phis);
	                  cthti = xcos(thti);
	                  ahv   = _mm512_div_pd(sphis,cthti);
	                  return (ahv);                           
	         }
	         
	         
	            /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-82
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_hh_f9182_zmm8r8(const __m512d phis) {
	                 
	                  register __m512d cphis,ahh;
	                  cphis = xcos(phis);
	                  ahh   = gms::math::negate_zmm8r8(cphis);
	                  return (ahh);
	          }
	          
	          
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_hh_f9182_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphis) {
	                 
	                  register __m512d phis = _mm512_load_pd(&pphis[0]);
	                  register __m512d cphis,ahh;
	                  cphis = xcos(phis);
	                  ahh   = gms::math::negate_zmm8r8(cphis);
	                  return (ahh);
	          }
	        
	       
	          __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_hh_f9182_zmm8r8_u(const double * __restrict  pphis) {
	                 
	                  register __m512d phis = _mm512_loadu_pd(&pphis[0]);
	                  register __m512d cphis,ahh;
	                  cphis = xcos(phis);
	                  ahh   = gms::math::negate_zmm8r8(cphis);
	                  return (ahh);
	          }
	          
	          
	             /*
	                Scattering matrix elements for
	                a perfectly conducting surface.
	                Formula: 9.1-83
	         */
	         
	         
	            __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_vh_f9183_zmm8r8(const __m512d phis,
	                                     const __m512d thts) {
	                                     
	                  register __m512d sphis,cthts;
	                  register __m512d ahv;
	                  sphis = xsin(phis);
	                  cthts = xcos(thti);
	                  ahv   = gms::math::negate_zmm8r8(_mm512_div_pd(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_vh_f9183_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pphis,
	                                       const double * __restrict __ATTR_ALIGN__(64) pthts) {
	                         
	                  register __m512d phis = _mm512_load_pd(&pphis[0]);
	                  register __m512d thts = _mm512_load_pd(&pphis[0]);            
	                  register __m512d sphis,cthts;
	                  register __m512d ahv;
	                  sphis = xsin(phis);
	                  cthts = xcos(thti);
	                  ahv   = gms::math::negate_zmm8r8(_mm512_div_pd(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	          __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_vh_f9183_zmm8r8_u(const double * __restrict  pphis,
	                                       const double * __restrict pthts) {
	                         
	                  register __m512d phis = _mm512_loadu_pd(&pphis[0]);
	                  register __m512d thts = _mm512_loadu_pd(&pphis[0]);            
	                  register __m512d sphis,cthts;
	                  register __m512d ahv;
	                  sphis = xsin(phis);
	                  cthts = xcos(thti);
	                  ahv   = gms::math::negate_zmm8r8(_mm512_div_pd(sphis,cthts));
	                  return (ahv);                           
	         }
	         
	         
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-84 
	         */
	         
	         
	           __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_vv_f9184_zmm8r8(const __m512d thti) {
	           
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d sthti,cthti,num,den;
	                  register __m512d avv;
	                  sthti = xsin(thti);
	                  cthti = xcos(thti);
	                  num   = _mm512_fmadd_pd(sthti,sthti,C10);
	                  den   = _mm512_mul_pd(cthti,cthti);
	                  avv   = _mm512_div_pd(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_vv_f9184_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pthti) {
	           
	                  register __m512d thti = _mm512_load_pd(&pthti[0]);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d sthti,cthti,num,den;
	                  register __m512d avv;
	                  sthti = xsin(thti);
	                  cthti = xcos(thti);
	                  num   = _mm512_fmadd_pd(sthti,sthti,C10);
	                  den   = _mm512_mul_pd(cthti,cthti);
	                  avv   = _mm512_div_pd(num,den);
	                  return (avv);
	           }
	           
	           
	           __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_vv_f9184_zmm8r8_u(const double * __restrict  pthti) {
	           
	                  register __m512d thti = _mm512_loadu_pd(&pthti[0]);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d sthti,cthti,num,den;
	                  register __m512d avv;
	                  sthti = xsin(thti);
	                  cthti = xcos(thti);
	                  num   = _mm512_fmadd_pd(sthti,sthti,C10);
	                  den   = _mm512_mul_pd(cthti,cthti);
	                  avv   = _mm512_div_pd(num,den);
	                  return (avv);
	           }
	           
	           
	           
	         /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-85
	         */                  
	         
	         
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d a_hh_f9185_zmm8r8() {
	           
	                  return _mm512_set1_pd(1.0f);
	           }  
	           
	           
	            /*
	              Backscattering from a perfectly conducting surface
	              Theta (inc) == theta (scat) , phi (scat) = 180 (grad).
	              Formula: 9.1-86
	         */   
	         
	         
	            __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d a_vh_f9186_zmm8r8() {
	           
	                  return _mm512_setzero_pd();
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
	           __m512d rcs_vv_f9187_zmm8r8(const __m512d k0,
	                                       const __m512d h,
	                                       const __m512d l,
	                                       const __m512d thti) {
	                                       
	                  const __m512d C40 = _mm512_set1_pd(4.0f);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d k04,x0,x1,l2,h2,sthti;
	                  register __m512d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm512_mul_pd(k0,k0);
	                  l2    = _mm512_mul_pd(l,l);
	                  k04   = _mm512_mul_pd(x0,x0);
	                  h2    = _mm512_mul_pd(h,h);
	                  x1    = xsin(thti);
	                  strm  = _mm512_fmadd_pd(x1,x1,C10);
	                  x0    = _mm512_mul_pd(C40,k04);
	                  sthti = _mm512_mul_pd(x1,x1);
	                  arg   = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
	                                        _mm512_mul_pd(sthti,l2));
	                  trm   = _mm512_mul_pd(_mm512_mul_pd(l2,h2,x0));
	                  earg  = xexp(arg);
	                  x1    = _mm512_mul_pd(strm,strm);
	                  inve  = _mm512_rcp14_pd(earg);
	                  rcs   = _mm512_mul_pd(trm,_mm512_mul_pd(x1,inve));
	                  return (rcs);
	         }
	           
	           
	           
	          __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d rcs_vv_f9187_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                       const double * __restrict __ATTR_ALIGN__(64) ph,
	                                       const double * __restrict __ATTR_ALIGN__(64) pl,
	                                       const double * __restrict __ATTR_ALIGN__(64) pthti) {
	                             
	                  register __m512d k0  = _mm512_load_pd(&pk0[0]);
	                  register __m512d h   = _mm512_load_pd(&ph[0]);
	                  register __m512d l   = _mm512_load_pd(&pl[0]);   
	                  register __m512d thti= _mm512_load_pd(&pthti[0]);       
	                  const __m512d C40 = _mm512_set1_pd(4.0f);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d k04,x0,x1,l2,h2,sthti;
	                  register __m512d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm512_mul_pd(k0,k0);
	                  l2    = _mm512_mul_pd(l,l);
	                  k04   = _mm512_mul_pd(x0,x0);
	                  h2    = _mm512_mul_pd(h,h);
	                  x1    = xsin(thti);
	                  strm  = _mm512_fmadd_pd(x1,x1,C10);
	                  x0    = _mm512_mul_pd(C40,k04);
	                  sthti = _mm512_mul_pd(x1,x1);
	                  arg   = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
	                                        _mm512_mul_pd(sthti,l2));
	                  trm   = _mm512_mul_pd(_mm512_mul_pd(l2,h2,x0));
	                  earg  = xexp(arg);
	                  x1    = _mm512_mul_pd(strm,strm);
	                  inve  = _mm512_rcp14_pd(earg);
	                  rcs   = _mm512_mul_pd(trm,_mm512_mul_pd(x1,inve));
	                  return (rcs);
	         }
	           
	            
	         
	           __ATTR_ALWAYS_INLINE__
	           
	           static inline
	           __m512d rcs_vv_f9187_zmm8r8_u(const double * __restrict  pk0,
	                                         const double * __restrict  ph,
	                                         const double * __restrict  pl,
	                                         const double * __restrict  pthti) {
	                             
	                  register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
	                  register __m512d h   = _mm512_loadu_pd(&ph[0]);
	                  register __m512d l   = _mm512_loadu_pd(&pl[0]);   
	                  register __m512d thti= _mm512_loadu_pd(&pthti[0]);       
	                  const __m512d C40 = _mm512_set1_pd(4.0f);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d k04,x0,x1,l2,h2,sthti;
	                  register __m512d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm512_mul_pd(k0,k0);
	                  l2    = _mm512_mul_pd(l,l);
	                  k04   = _mm512_mul_pd(x0,x0);
	                  h2    = _mm512_mul_pd(h,h);
	                  x1    = xsin(thti);
	                  strm  = _mm512_fmadd_pd(x1,x1,C10);
	                  x0    = _mm512_mul_pd(C40,k04);
	                  sthti = _mm512_mul_pd(x1,x1);
	                  arg   = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
	                                        _mm512_mul_pd(sthti,l2));
	                  trm   = _mm512_mul_pd(_mm512_mul_pd(l2,h2,x0));
	                  earg  = xexp(arg);
	                  x1    = _mm512_mul_pd(strm,strm);
	                  inve  = _mm512_rcp14_pd(earg);
	                  rcs   = _mm512_mul_pd(trm,_mm512_mul_pd(x1,inve));
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
	           __m512d rcs_hh_f9188_zmm8r8(const __m512d k0,
	                                       const __m512d h,
	                                       const __m512d l,
	                                       const __m512d thti) {
	                                       
	                  const __m512d C40 = _mm512_set1_pd(4.0f);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d k04,x0,x1,l2,h2,sthti,cthti;
	                  register __m512d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm512_mul_pd(k0,k0);
	                  l2    = _mm512_mul_pd(l,l);
	                  strm  = xcos(thti);
	                  k04   = _mm512_mul_pd(x0,x0);
	                  h2    = _mm512_mul_pd(h,h);
	                  x1    = xsin(thti);
	                  x0    = _mm512_mul_pd(C40,k04);
	                  sthti = _mm512_mul_pd(x1,x1);
	                  arg   = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
	                                        _mm512_mul_pd(sthti,l2));
	                  trm   = _mm512_mul_pd(_mm512_mul_pd(l2,h2,x0));
	                  x1    = _mm512_mul_pd(strm,strm);
	                  earg  = xexp(arg);
	                  cthti = _mm512_mul_pd(x1,x1);
	                  inve  = _mm512_rcp14_pd(earg);
	                  rcs   = _mm512_mul_pd(trm,_mm512_mul_pd(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	           __ATTR_ALWAYS_INLINE__
	          
	           static inline
	           __m512d rcs_hh_f9188_zmm8r8_a(const double * __restrict __ATTR_ALIGN__(64) pk0,
	                                       const double * __restrict __ATTR_ALIGN__(64) ph,
	                                       const double * __restrict __ATTR_ALIGN__(64) pl,
	                                       const double * __restrict __ATTR_ALIGN__(64) pthti) {
	                             
	                  register __m512d k0  = _mm512_load_pd(&pk0[0]);
	                  register __m512d h   = _mm512_load_pd(&ph[0]);
	                  register __m512d l   = _mm512_load_pd(&pl[0]);   
	                  register __m512d thti= _mm512_load_pd(&pthti[0]); 
	                                       
	                  const __m512d C40 = _mm512_set1_pd(4.0f);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d k04,x0,x1,l2,h2,sthti,cthti;
	                  register __m512d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm512_mul_pd(k0,k0);
	                  l2    = _mm512_mul_pd(l,l);
	                  strm  = xcos(thti);
	                  k04   = _mm512_mul_pd(x0,x0);
	                  h2    = _mm512_mul_pd(h,h);
	                  x1    = xsin(thti);
	                  x0    = _mm512_mul_pd(C40,k04);
	                  sthti = _mm512_mul_pd(x1,x1);
	                  arg   = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
	                                        _mm512_mul_pd(sthti,l2));
	                  trm   = _mm512_mul_pd(_mm512_mul_pd(l2,h2,x0));
	                  x1    = _mm512_mul_pd(strm,strm);
	                  earg  = xexp(arg);
	                  cthti = _mm512_mul_pd(x1,x1);
	                  inve  = _mm512_rcp14_pd(earg);
	                  rcs   = _mm512_mul_pd(trm,_mm512_mul_pd(cthti,inve));
	                  return (rcs);
	         }
	         
	         
	            __ATTR_ALWAYS_INLINE__
	         
	           static inline
	           __m512d rcs_hh_f9188_zmm8r8_u(const double * __restrict  pk0,
	                                       const double * __restrict  ph,
	                                       const double * __restrict  pl,
	                                       const double * __restrict  pthti) {
	                             
	                  register __m512d k0  = _mm512_loadu_pd(&pk0[0]);
	                  register __m512d h   = _mm512_loadu_pd(&ph[0]);
	                  register __m512d l   = _mm512_loadu_pd(&pl[0]);   
	                  register __m512d thti= _mm512_loadu_pd(&pthti[0]); 
	                                       
	                  const __m512d C40 = _mm512_set1_pd(4.0f);
	                  const __m512d C10 = _mm512_set1_pd(1.0f);
	                  register __m512d k04,x0,x1,l2,h2,sthti,cthti;
	                  register __m512d rcs,arg,earg,inve,strm,trm;
	                  x0    = _mm512_mul_pd(k0,k0);
	                  l2    = _mm512_mul_pd(l,l);
	                  strm  = xcos(thti);
	                  k04   = _mm512_mul_pd(x0,x0);
	                  h2    = _mm512_mul_pd(h,h);
	                  x1    = xsin(thti);
	                  x0    = _mm512_mul_pd(C40,k04);
	                  sthti = _mm512_mul_pd(x1,x1);
	                  arg   = _mm512_mul_pd(_mm512_mul_pd(k0,k0),
	                                        _mm512_mul_pd(sthti,l2));
	                  trm   = _mm512_mul_pd(_mm512_mul_pd(l2,h2,x0));
	                  x1    = _mm512_mul_pd(strm,strm);
	                  earg  = xexp(arg);
	                  cthti = _mm512_mul_pd(x1,x1);
	                  inve  = _mm512_rcp14_pd(earg);
	                  rcs   = _mm512_mul_pd(trm,_mm512_mul_pd(cthti,inve));
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
	           __m512d rcs_vhhv_f9189_zmm8r8() {
	           
	                  return (_mm512_setzero_pd()); 
	          }
	          
	          
	          
	           
	      
	                                    
	        
                 
                 
                 
               
               
         
     } // radiolocation


} // gms




















#endif /*__GMS_RCS_COMPLEX_OBJECTS_ZMM8R8_H__*/
