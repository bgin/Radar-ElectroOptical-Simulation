
#ifndef __GMS_RCS_SPHERE_XMM2R8_H__
#define __GMS_RCS_SPHERE_XMM2R8_H__ 110920240639


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

    const unsigned int GMS_RCS_SPHERE_XMM2R8_MAJOR = 1U;
    const unsigned int GMS_RCS_SPHERE_XMM2R8_MINOR = 0U;
    const unsigned int GMS_RCS_SPHERE_XMM2R8_MICRO = 0U;
    const unsigned int GMS_RCS_SPHERE_XMM2R8_FULLVER =
      1000U*GMS_RCS_SPHERE_XMM2R8_MAJOR+
      100U*GMS_RCS_SPHERE_XMM2R8_MINOR+
      10U*GMS_RCS_SPHERE_XMM2R8_MICRO;
    const char * const GMS_RCS_SPHERE_XMM2R8_CREATION_DATE = "11-09-2024 06:39 PM +00200 (WED 11 SEP 2024 GMT+2)";
    const char * const GMS_RCS_SPHERE_XMM2R8_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_SPHERE_XMM2R8_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_SPHERE_XMM2R8_DESCRIPTION   = "SSE optimized Sphere Radar Cross Section (analytic) functionality."

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"


namespace gms {

         namespace radiolocation {

                  /*
                       Radar Cross Section Handbook 1, page 147, formula 3.2-4
                       Backscattering function ,resonance region 0.4 .le. k0a .le. 20.0
                       Theta = 0, far-field
                       Valid for k0a < 1 only!!
                   */

               
                   __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void bsc_func_324_xmm2r8(const __m128d k0a, // size of sphere expressed in wavenumber units
                                               __m128d * __restrict F0r, // the results
                                               __m128d * __restrict F0i); 
                                               
	           
                 
	            __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void bsc_func_324_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a, // size of sphere expressed in wavenumber units
                                               __m128d * __restrict F0r, // the results
                                               __m128d * __restrict F0i); 
                                               
                                               
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)                        
                   void bsc_func_324_xmm2r8_u(const double * __restrict  pk0a, // size of sphere expressed in wavenumber units
                                               __m128d * __restrict F0r, // the results
                                               __m128d * __restrict F0i); 


                  /*
                        Radar Cross Section Handbook 1, page 147, formula 3.2-5
                        Backscattering cross section
                        
                    */
                  
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f325_xmm2r8(const __m128d k0,
                                           const __m128d a ); 
                                           

                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f325_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0,
                                             const double * __restrict __ATTR_ALIGN__(16) pa ); 
                                             


                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f325_xmm2r8_u(const double * __restrict pk0,
                                             const double * __restrict pa ); 






                  /*
                        Creeping wave term, F c(0) at the upper end of the resonance region.
                        Formula 3.2-8
                    */
                  
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void bsc_f328_xmm2r8(const __m128d x,//k0a
                                         __m128d * __restrict Fc0r,
                                         __m128d * __restrict Fc0i); 

                  
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void bsc_f328_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) px,//k0a
                                           double * __restrict __ATTR_ALIGN__(16) Fc0r,
                                           double * __restrict __ATTR_ALIGN__(16) Fc0i); 


                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void bsc_f328_xmm2r8_u(const double * __restrict px,//k0a
                                           double * __restrict  Fc0r,
                                           double * __restrict  Fc0i);



                  



                  /*
                       The complex scattering amplitudes near the lower end of the resonance
                       region (say, 0.4 < k^a < 1) 
                       E-plane approximation
                       These equations are not valid at all for k0a > 1. They are
                       valid for all theta angles.
                   */

                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3213_xmm2r8( const __m128d k0a,
                                           const __m128d tht,
                                           __m128d * __restrict S1r,
                                           __m128d * __restrict S1i); 


                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3213_xmm2r8_a( const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(16) ptht,
                                             double * __restrict __ATTR_ALIGN__(16) S1r,
                                             double * __restrict __ATTR_ALIGN__(16) S1i); 


                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3213_xmm2r8_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict S1r,
                                             double * __restrict S1i); 


                   /*
                       The complex scattering amplitudes near the lower end of the resonance
                       region (say, 0.4 < k^a < 1) 
                       H-plane approximation
                       These equations are not valid at all for k0a > 1. They are
                       valid for all theta angles.
                   */

                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3214_xmm2r8( const __m128d k0a,
                                           const __m128d tht,
                                           __m128d * __restrict S2r,
                                           __m128d * __restrict S2i);


                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3214_xmm2r8_a( const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(16) ptht,
                                             double * __restrict __ATTR_ALIGN__(16) S2r,
                                             double * __restrict __ATTR_ALIGN__(16) S2i); 


                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3214_xmm2r8_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict  S2r,
                                             double * __restrict  S2i); 


                  
                  /*
                       Formula 3.2-16, optics contribution at upper end of resonance region
                   */
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3216_xmm2r8(const __m128d k0a,
                                         const __m128d tht,
                                         __m128d * __restrict S1r,
                                         __m128d * __restrict S1i);
                                         

                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3216_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(16) ptht,
                                           double * __restrict __ATTR_ALIGN__(16) S1r,
                                           double * __restrict __ATTR_ALIGN__(16) S1i); 


                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3216_xmm2r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  ptht,
                                           double * __restrict  S1r,
                                           double * __restrict  S1i); 


                   /*
                       Formula 3.2-17, optics contribution at upper end of resonance region
                   */
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3217_xmm2r8(const __m128d k0a,
                                         const __m128d tht,
                                         __m128d * __restrict S2r,
                                         __m128d * __restrict S2i); 


                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3217_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(16) ptht,
                                           double * __restrict __ATTR_ALIGN__(16) S2r,
                                           double * __restrict __ATTR_ALIGN__(16) S2i); 


                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3217_xmm2r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  ptht,
                                           double * __restrict  S2r,
                                           double * __restrict  S2i); 


                 /*
                        Low frequency region (k0a < 0.4) i.e. Rayleigh-region complex scattering
                        amplitudes.
                        Formula 3.2-20
                    */
                 
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d S1_f3220_xmm2r8(const __m128d k0a,
                                           const __m128d tht); 


                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3220_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(16) ptht,
                                           double * __restrict __ATTR_ALIGN__(16) S1); 

                 
                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3220_xmm2r8_u(const double * __restrict pk0a,
                                           const double * __restrict ptht,
                                           double * __restrict S1);


                   /*
                        Low frequency region (k0a < 0.4) i.e. Rayleigh-region complex scattering
                        amplitudes.
                        Formula 3.2-21
                    */

                     __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d S2_f3221_xmm2r8(const __m128d k0a,
                                           const __m128d tht); 

                   
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3221_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                           const double * __restrict __ATTR_ALIGN__(16) ptht,
                                           double * __restrict S2); 

                  
                    __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f3221_xmm2r8_u(const double * __restrict  pk0a,
                                           const double * __restrict  ptht,
                                           double * __restrict S2); 

                 /*
                         The E-plane and H-plane scattering cross sections.
                         Formula 3.2-22
                    */

                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3222_xmm2r8(const __m128d k0a,
                                            const __m128d a,
                                            const __m128d theta); 
                                            
                                            
                       __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)                       
                   void rcs_f3222_xmm2r8_a(  const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(16) pa,
                                              const double * __restrict __ATTR_ALIGN__(16) ptheta,
                                              double * __restrict __ATTR_ALIGN__(16) rcs ); 


                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3222_xmm2r8_u(  const double * __restrict  pk0a,
                                              const double * __restrict  pa,
                                              const double * __restrict ptheta,
                                              double * __restrict  rcs ); 


                 /*
                         The E-plane and H-plane scattering cross sections.
                         Formula 3.2-23
                    */

                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3223_xmm2r8(const __m128d k0a,
                                            const __m128d a,
                                            const __m128d theta); 

                    
                      __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3223_xmm2r8_a(  const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                              const double * __restrict __ATTR_ALIGN__(16) pa,
                                              const double * __restrict __ATTR_ALIGN__(16) ptheta,
                                              double * __restrict __ATTR_ALIGN__(16) rcs  ); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3223_xmm2r8_u(  const double * __restrict  pk0a,
                                              const double * __restrict  pa,
                                              const double * __restrict  ptheta,
                                              double * __restrict  rcs  ); 


                 /*
                        High frequency region (k0a > 20).
                        Complex scattering amplitudes.
                        Formula 3.2-24
                   */


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S12_f3224_xmm2r8( const __m128d k0a,
                                           const __m128d tht,
                                           __m128d * __restrict S12r,
                                           __m128d * __restrict S12i); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S12_f3224_xmm2r8_a( const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                             const double * __restrict __ATTR_ALIGN__(16) ptht,
                                             double * __restrict __ATTR_ALIGN__(16) S12r,
                                             double * __restrict __ATTR_ALIGN__(16) S12i);


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S12_f3224_xmm2r8_u( const double * __restrict  pk0a,
                                             const double * __restrict  ptht,
                                             double * __restrict  S12r,
                                             double * __restrict  S12i); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3225_xmm2r8(const __m128d a);


                /*
                       Resonance region (0.4 < k0a < 20.0), equations are valid only for k0a < 1.0.
                       Complex scattering amplitude represented as a scattering function -- formula 3.2-26
                 */

                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S12_f3226_xmm2r8(const __m128d k0a,
                                          __m128d * __restrict S12r,
                                          __m128d * __restrict S12i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S12_f3226_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                            double * __restrict __ATTR_ALIGN__(16) S12r,
                                            double * __restrict __ATTR_ALIGN__(16) S12i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S12_f3226_xmm2r8_u(const double * __restrict  pk0a,
                                            double * __restrict  S12r,
                                            double * __restrict  S12i); 


                 /*
                           Resonance region (0.4 < k0a < 20.0), equations are valid only for k0a < 1.0.
                           Radar cross-section, formula 3.2-27
                     */

                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3227_xmm2r8(const __m128d k0a,
                                            const __m128d a); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3227_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(16) pa,
                                            double * __restrict __ATTR_ALIGN__(16) rcs); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3227_xmm2r8_u(const double * __restrict  pk0a,
                                            const double * __restrict  pa,
                                            double * __restrict  rcs); 


                 /*
                       Scattering functions equation at upper end of resonance region.
                       Optics wave term and creeping wave term.
                       Optics wave term, formula 3.2-28
                   */

                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void FO_f3228_xmm2r8(const __m128d k0a,
                                         __m128d * __restrict FOr,
                                         __m128d * __restrict FOi); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void FO_f3228_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                           double * __restrict __ATTR_ALIGN__(16) FOr,
                                           double * __restrict __ATTR_ALIGN__(16) FOi); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void FO_f3228_xmm2r8_u(const double * __restrict pk0a,
                                           double * __restrict FOr,
                                           double * __restrict FOi); 


                    /*
                       Scattering functions equation at upper end of resonance region.
                       Optics wave term and creeping wave term.
                       Creeping wave term, formula 3.2-29
                   */
                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void FC_f3229_xmm2r8(const __m128d x,
                                         __m128d * __restrict FCr,
                                         __m128d * __restrict FCi); 
                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void FC_f3229_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) px,
                                         double * __restrict __ATTR_ALIGN__(16) FCr,
                                         double * __restrict __ATTR_ALIGN__(16) FCi); 
                 
                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void FC_f3229_xmm2r8_u(const double * __restrict  px,
                                           double * __restrict  FCr,
                                           double * __restrict  FCi); 

                 /*
                         Low frquency region (k0a < 0.4), Rayleigh approximation
                         for forward scattering function and cross-section.
                         Formulae: 3.2-31, 3.2-32
                         Forward scattering function.
                   */

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d F_f3231_xmm2r8(const __m128d k0a); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void F_f3231_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                          double * __restrict __ATTR_ALIGN__(16) Fpi); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void F_f3231_xmm2r8_u(const double * __restrict  pk0a,
                                          double * __restrict  Fpi); 


                    /*
                         Low frquency region (k0a < 0.4), Rayleigh approximation
                         for forward scattering function and cross-section.
                         Formulae: 3.2-31, 3.2-32
                         RCS.
                   */


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3232_xmm2r8(const __m128d k0a,
                                            const __m128d a);


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3232_xmm2r8_a(  const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(16) pa,
                                            double * __restrict __ATTR_ALIGN__(16) rcs); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3232_xmm2r8_u(  const double * __restrict  pk0a,
                                              const double * __restrict    pa,
                                              double * __restrict  rcs); 


                   /*
                         High-frequency region -- forward scattering function and 
                         cross-section (k0a > 20)
                         Formula 3.2-33 (Forward scattering function).
                     */

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void F_f3233_xmm2r8(const __m128d k0a,
                                        __m128d * __restrict Fr,
                                        __m128d * __restrict Fi); 

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void F_f3233_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                          double * __restrict __ATTR_ALIGN__(16) Fr,
                                          double * __restrict __ATTR_ALIGN__(16) Fi); 


                    
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void F_f3233_xmm2r8_u(const double * __restrict  pk0a,
                                          double * __restrict  Fr,
                                          double * __restrict  Fi); 

                   /*
                         High-frequency region -- forward scattering function and 
                         cross-section (k0a > 20)
                         Formula 3.2-34 (RCS).
                     */

 
                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3234_xmm2r8(const __m128d k0a,
                                            const __m128d a); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3234_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                            const double * __restrict __ATTR_ALIGN__(16) pa,
                                            double * __restrict __ATTR_ALIGN__(16) rcs); 


                   
	           
	           
                  
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f3234_xmm2r8_u(const double * __restrict  pk0a,
                                            const double * __restrict  pa,
                                            double * __restrict  rcs); 


                  /*
                          Low-frequency region (k1a < 0.8).
                          Expansion by two series terms i.e. A0,A1 and B0,B1.
                          Formula 3.3-5
                    */
                   
                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void A_coffs_f335_xmm2r8(const __m128d k0a5,
                                         const __m128d m1r,
                                         const __m128d m1i,
                                         __m128d * __restrict A1r,
                                         __m128d * __restrict A1i,
                                         __m128d * __restrict A2r,
                                         __m128d * __restrict A2i); 

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void A_coffs_f335_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pk0a5,
                                               const double * __restrict __ATTR_ALIGN__(16) pm1r,
                                               const double * __restrict __ATTR_ALIGN__(16) pm1i,
                                               double * __restrict __ATTR_ALIGN__(16) A1r,
                                               double * __restrict __ATTR_ALIGN__(16) A1i,
                                               double * __restrict __ATTR_ALIGN__(16) A2r,
                                               double * __restrict __ATTR_ALIGN__(16) A2i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void A_coffs_f335_xmm2r8_u(const double * __restrict  pk0a5,
                                               const double * __restrict  pm1r,
                                               const double * __restrict  pm1i,
                                               double * __restrict  A1r,
                                               double * __restrict  A1i,
                                               double * __restrict  A2r,
                                               double * __restrict  A2i);


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void B_coffs_f335_xmm2r8(const __m128d k0a3,
                                             const __m128d k0a5,
                                             const __m128d m1r,
                                             const __m128d m1i,
                                             __m128d * __restrict B1r,
                                             __m128d * __restrict B1i,
                                             __m128d * __restrict B2r,
                                             __m128d * __restrict B2i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void B_coffs_f335_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16)  pk0a3,
                                               const double * __restrict __ATTR_ALIGN__(16)  pk0a5,
                                               const double * __restrict __ATTR_ALIGN__(16)  pm1r,
                                               const double * __restrict __ATTR_ALIGN__(16)  pm1i,
                                               double * __restrict __ATTR_ALIGN__(16) B1r,
                                               double * __restrict __ATTR_ALIGN__(16) B1i,
                                               double * __restrict __ATTR_ALIGN__(16) B2r,
                                               double * __restrict __ATTR_ALIGN__(16) B2i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void B_coffs_f335_xmm2r8_u(const double * __restrict   pk0a3,
                                               const double * __restrict   pk0a5,
                                               const double * __restrict   pm1r,
                                               const double * __restrict   pm1i,
                                               double * __restrict  B1r,
                                               double * __restrict  B1i,
                                               double * __restrict  B2r,
                                               double * __restrict  B2i); 


                   /*
                         Rayleigh backscattering RCS for dielectric spheres at angle 0.
                         Formula 3.3-7
                     */

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f337_xmm2r8(const __m128d a,
                                           const __m128d k0a4,
                                           const __m128d m1r,
                                           const __m128d m1i); 
                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f337_xmm2r8_a(  const double * __restrict __ATTR_ALIGN__(16) pa,
                                           const double * __restrict __ATTR_ALIGN__(16) pk0a4,
                                           const double * __restrict __ATTR_ALIGN__(16) pm1r,
                                           const double * __restrict __ATTR_ALIGN__(16) pm1i,
                                           double * __restrict __ATTR_ALIGN__(16) rcs ); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void rcs_f337_xmm2r8_u(const double * __restrict  pa,
                                           const double * __restrict  pk0a4,
                                           const double * __restrict  pm1r,
                                           const double * __restrict  pm1i,
                                           double * __restrict  rcs );


                 /*
                        Low-frequency bi-static scattering
                  */

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f338_xmm2r8(const __m128d ka03,
                                        const __m128d ka05,
                                        const __m128d tht,
                                        const __m128d mm1r,
                                        const __m128d mm1i,
                                        __m128d * __restrict S1r,
                                        __m128d * __restrict S1i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f338_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pka03,
                                        const double * __restrict __ATTR_ALIGN__(16) pka05,
                                        const double * __restrict __ATTR_ALIGN__(16) ptht,
                                        const double * __restrict __ATTR_ALIGN__(16) pmm1r,
                                        const double * __restrict __ATTR_ALIGN__(16) pmm1i,
                                        double * __restrict __ATTR_ALIGN__(16) S1r,
                                        double * __restrict __ATTR_ALIGN__(16) S1i);


                  
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f338_xmm2r8_u(const double * __restrict  pka03,
                                        const double * __restrict  pka05,
                                        const double * __restrict  ptht,
                                        const double * __restrict  pmm1r,
                                        const double * __restrict  pmm1i,
                                        double * __restrict  S1r,
                                        double * __restrict  S1i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f338_xmm2r8(const __m128d ka03,
                                        const __m128d ka05,
                                        const __m128d tht,
                                        const __m128d mm1r,
                                        const __m128d mm1i,
                                        __m128d * __restrict S2r,
                                        __m128d * __restrict S2i); 


                 
                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f338_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pka03,
                                          const double * __restrict __ATTR_ALIGN__(16) pka05,
                                          const double * __restrict __ATTR_ALIGN__(16) ptht,
                                          const double * __restrict __ATTR_ALIGN__(16) pmm1r,
                                          const double * __restrict __ATTR_ALIGN__(16) pmm1i,
                                          double * __restrict __ATTR_ALIGN__(16) S2r,
                                          double * __restrict __ATTR_ALIGN__(16) S2i); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S2_f338_xmm2r8_u(const double * __restrict  pka03,
                                          const double * __restrict  pka05,
                                          const double * __restrict  ptht,
                                          const double * __restrict  pmm1r,
                                          const double * __restrict  pmm1i,
                                          double * __restrict  S2r,
                                          double * __restrict  S2i); 

                   /*
                         E-plane and H-plane RCS.
                         Formulae: 3.3-10,3.3-11
                     */
                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3310_xmm2r8(const __m128d tht,
                                            const __m128d a,
                                            const __m128d ka04,
                                            const __m128d mm1r,
                                            const __m128d mm1i); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3310_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) ptht,
                                            const double * __restrict __ATTR_ALIGN__(16) pa,
                                            const double * __restrict __ATTR_ALIGN__(16) pka04,
                                            const double * __restrict __ATTR_ALIGN__(16) pmm1r,
                                            const double * __restrict __ATTR_ALIGN__(16) pmm1i); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3310_xmm2r8_u(const double * __restrict  ptht,
                                              const double * __restrict  pa,
                                              const double * __restrict  pka04,
                                              const double * __restrict  pmm1r,
                                              const double * __restrict  pmm1i); 

                    
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3311_xmm2r8(const __m128d tht,
                                            const __m128d a,
                                            const __m128d ka04,
                                            const __m128d mm1r,
                                            const __m128d mm1i); 


                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3311_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) ptht,
                                              const double * __restrict __ATTR_ALIGN__(16) pa,
                                              const double * __restrict __ATTR_ALIGN__(16) pka04,
                                              const double * __restrict __ATTR_ALIGN__(16) pmm1r,
                                              const double * __restrict __ATTR_ALIGN__(16) pmm1i); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3311_xmm2r8_u(const double * __restrict  ptht,
                                              const double * __restrict  pa,
                                              const double * __restrict  pka04,
                                              const double * __restrict  pmm1r,
                                              const double * __restrict  pmm1i); 


                    /*
                          Bistatic Geometric Optics Rays.
                          The RCS of sphere included N-rays.
                          E-plane or H-plane, for backscattering
                          and forward-scattering E and H RCS are 
                          identical
                     */

                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3314_xmm2r8(const __m128d * __restrict Snthr,
                                            const __m128d * __restrict Snthi,
                                            const __m128d * __restrict cjphr,
                                            const __m128d * __restrict cjphi,
                                            __m128d * __restrict wrkr,
                                            __m128d * __restrict wrki,
                                            const __m128d k02,
                                            const int32_t N); 


                  /*
                         Large sphere limit, k0a > 1.15/m1 (reflective region).
                         Backscattering RCS, formula 3.3-17
                    */

                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3317_xmm2r8(const __m128d m1r,
                                            const __m128d m1i,
                                            const __m128d a); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3317_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pm1r,
                                              const double * __restrict __ATTR_ALIGN__(16) pm1i,
                                              const double * __restrict __ATTR_ALIGN__(16) pa); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3317_xmm2r8_u(const double * __restrict  pm1r,
                                              const double * __restrict  pm1i,
                                              const double * __restrict  pa); 


                    /*
                       Forward scattering RCS.
                       Formula 3.3-19
                         */

                   
	           
	           
                 
	               __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3319_xmm2r8(const __m128d a,
                                            const __m128d k0a); 


                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3319_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pa,
                                            const double * __restrict __ATTR_ALIGN__(16) pk0a); 

                 
                   
	           
	           
                 
	                __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   __m128d rcs_f3319_xmm2r8_u(const double * __restrict  pa,
                                              const double * __restrict  pk0a); 


                    /*
                         Approximate solutions for far-field region (Rayleigh-Gans)
                         (abs(m1-1) << 1,2*k0a abs(m1-1) << 1)
                         Bistatic scattering formula 3.3-22
                     */

                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3322_xmm2r8(const __m128d m1r,
                                         const __m128d m1i,
                                         const __m128d tht,
                                         const __m128d k0a,
                                         __m128d * __restrict S1r,
                                         __m128d * __restrict S1i); 

                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3322_xmm2r8_a(const double * __restrict __ATTR_ALIGN__(16) pm1r,
                                           const double * __restrict __ATTR_ALIGN__(16) pm1i,
                                           const double * __restrict __ATTR_ALIGN__(16) ptht,
                                           const double * __restrict __ATTR_ALIGN__(16) pk0a,
                                           double * __restrict __ATTR_ALIGN__(16) S1r,
                                           double * __restrict __ATTR_ALIGN__(16) S1i); 


                   
	           
	           
                 
	              __ATTR_HOT__
                   __ATTR_VECTORCALL__
                   __ATTR_ALIGN__(32)
                   void S1_f3322_xmm2r8_u(const double * __restrict  pm1r,
                                           const double * __restrict  pm1i,
                                           const double * __restrict  ptht,
                                           const double * __restrict  pk0a,
                                           double * __restrict S1r,
                                           double * __restrict  S1i); 


                  

     } // radiolocation

} // gms











#endif /*__GMS_RCS_SPHERE_XMM2R8_H__*/
