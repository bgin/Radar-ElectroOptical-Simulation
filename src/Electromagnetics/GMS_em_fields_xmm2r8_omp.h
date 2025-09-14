
#ifndef __GMS_EM_FIELDS_XMM2R8_OMP_H__
#define __GMS_EM_FIELDS_XMM2R8_OMP_H__ 

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

    const unsigned int GMS_EM_FIELDS_XMM2R8_OMP_MAJOR = 1U;
    const unsigned int GMS_EM_FIELDS_XMM2R8_OMP_MINOR = 0U;
    const unsigned int GMS_EM_FIELDS_XMM2R8_OMP_MICRO = 0U;
    const unsigned int GMS_EM_FIELDS_XMM2R8_OMP_FULLVER =
      1000U*GMS_EM_FIELDS_XMM2R8_OMP_MAJOR+
      100U*GMS_EM_FIELDS_XMM2R8_OMP_MINOR+
      10U*GMS_EM_FIELDS_XMM2R8_OMP_MICRO;
    const char * const GMS_EM_FIELDS_XMM2R8_OMP_CREATION_DATE = "30-10-2023 11:16 AM +00200 (30 OCT 2023 GMT+2)";
    const char * const GMS_EM_FIELDS_XMM2R8_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EM_FIELDS_XMM2R8_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EM_FIELDS_XMM2R8_OMP_DESCRIPTION   = "Computational ElectroMagnetics related helper routines OpenMP-multithreaded."
                       

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {


            namespace radiolocation {
            
                    __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_NO_STACK_PROTECTOR__ 
	            void sdotv_xmm2r8_unroll6x_omp( const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) pdtv,
	                                            const int32_t n); 
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__   
	           void cdotv_xmm2c8_unroll10x_omp(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);   
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                 
	          void cdotv_xmm2c8_unroll6x_omp(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	          void cdotv_xmm2c8_unroll2x_omp(  const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	          void cdotv_xmm2c8_rolled_omp(    const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv2z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                
	          void cnorm_xmm2c8_unroll10x_omp(const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                
	          void cnorm_xmm2c8_unroll6x_omp(const xmm2c8_t   * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                             
	          void cnorm_xmm2c8_unroll2x_omp(   const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                  
	           void cnorm_xmm2c8_rolled_omp(    const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1x,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1y,
	                                            const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pv1z,
	                                            xmm2c8_t * __restrict __ATTR_ALIGN__(16) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                  
	           void scrossc_xmm2r8_unroll10x_omp( const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	           void scrossc_xmm2r8_unroll6x_omp(  const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                       
	          void scrossc_xmm2r8_unroll2x_omp(  const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                   
	          void scrossc_xmm2r8_rolled_omp(    const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1y, 
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const xmm2c8_t  * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presx,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presy,
	                                              xmm2c8_t * __restrict __ATTR_ALIGN__(16) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	          void scrossv_xmm2r8_unroll10x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                   
	          void scrossv_xmm2r8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m128d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m128d * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m128d * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m128d * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                              
	          void scrossv_xmm2r8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                  
	          void scrossv_xmm2r8_rolled_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pv1x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv1y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv1z,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2x,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2y,
	                                              const __m128d * __restrict __ATTR_ALIGN__(16) pv2z,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcx,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcy,
	                                              __m128d * __restrict __ATTR_ALIGN__(16) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                  
	          void dir_vec_xmm2r8_unroll10x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                
	           void dir_vec_xmm2r8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                
	           void dir_vec_xmm2r8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                 
	          void dir_vec_xmm2r8_rolled_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST); 
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                 
	          void pol_vec_xmm2r8_unroll10x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                
	           void pol_vec_xmm2r8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m128d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m128d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m128d * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m128d * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m128d * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                
	          void pol_vec_xmm2r8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                   
	           void pol_vec_xmm2r8_rolled_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                            const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) ppvx,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) ppvy,
	                                            __m128d * __restrict __ATTR_ALIGN__(16) ppvz,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                 
	           void H_XYZ_VP_xmm2c8_unroll10x_omp( const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                    
	           void H_XYZ_VP_xmm2c8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	           void H_XYZ_VP_xmm2c8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                    
	           void H_XYZ_VP_xmm2c8_rolled_omp(    const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	           void B_XYZ_VP_xmm2c8_unroll10x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                   
	           void B_XYZ_VP_xmm2c8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                   
	           void B_XYZ_VP_xmm2c8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                    
	           void B_XYZ_VP_xmm2c8_rolled_omp(const __m128d * __restrict __ATTR_ALIGN__(16) pvpolx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpoly,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvpolz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdiry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvdirz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrx,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvry,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pvrz,
	                                               const __m128d * __restrict __ATTR_ALIGN__(16) pomega,
	                                               const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pk,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                               xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_NO_STACK_PROTECTOR__                                     
	           void B_XYZ_H_XYZ_P_xmm2c8_unroll10x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_NO_STACK_PROTECTOR__                                    
	            void B_XYZ_H_XYZ_P_xmm2c8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_NO_STACK_PROTECTOR__                                      
	          void B_XYZ_H_XYZ_P_xmm2c8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)  
	        __ATTR_NO_STACK_PROTECTOR__                                       
	        void B_XYZ_H_XYZ_P_xmm2c8_rolled_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppsi,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppx,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppy,
	                                                const __m128d * __restrict __ATTR_ALIGN__(16) ppz,
	                                                const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pr,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pH_z,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_x,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_y,
	                                                xmm2c8_t * __restrict __ATTR_ALIGN__(16) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)  
	        __ATTR_NO_STACK_PROTECTOR__                                           
	       void B_XYZ_H_XYZ_EP_xmm2c8_unroll10x_omp( const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                     const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                      const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)  
	        __ATTR_NO_STACK_PROTECTOR__                                              
	       void B_XYZ_H_XYZ_EP_xmm2c8_unroll6x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                     const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                      const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)  
	        __ATTR_NO_STACK_PROTECTOR__                                              
	       void B_XYZ_H_XYZ_EP_xmm2c8_unroll2x_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                     const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                      const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)  
	        __ATTR_NO_STACK_PROTECTOR__                                             
	       void B_XYZ_H_XYZ_EP_xmm2c8_rolled_omp(const __m128d * __restrict __ATTR_ALIGN__(16) ptht,
	                                                     const __m128d * __restrict __ATTR_ALIGN__(16) pphi,
	                                                      const __m128d * __restrict __ATTR_ALIGN__(16) pomg,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) pphase,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) prefi,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppx,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppy,
	                                                     const xmm2c8_t * __restrict __ATTR_ALIGN__(16) ppz,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pH_z,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_x,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_y,
	                                                     xmm2c8_t * __restrict __ATTR_ALIGN__(16)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	                                                     
	         
	                                           
	                                            
	                                            
	                                           
	                                                     
         } // radiolocation

} // gms











#endif /*__GMS_EM_FIELDS_XMM2R8_OMP_H__*/
