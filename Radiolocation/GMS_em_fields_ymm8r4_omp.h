
#ifndef __GMS_EM_FIELDS_YMM8R4_OMP_H__
#define __GMS_EM_FIELDS_YMM8R4_OMP_H__ 281020231539

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

    const unsigned int GMS_EM_FIELDS_YMM8R4_OMP_MAJOR = 1U;
    const unsigned int GMS_EM_FIELDS_YMM8R4_OMP_MINOR = 0U;
    const unsigned int GMS_EM_FIELDS_YMM8R4_OMP_MICRO = 0U;
    const unsigned int GMS_EM_FIELDS_YMM8R4_OMP_FULLVER =
      1000U*GMS_EM_FIELDS_YMM8R4_OMP_MAJOR+
      100U*GMS_EM_FIELDS_YMM8R4_OMP_MINOR+
      10U*GMS_EM_FIELDS_YMM8R4_OMP_MICRO;
    const char * const GMS_EM_FIELDS_YMM8R4_OMP_CREATION_DATE = "28-10-2023 15:39  +00200 (28 OCT 2023 GMT+2)";
    const char * const GMS_EM_FIELDS_YMM8R4_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EM_FIELDS_YMM8R4_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EM_FIELDS_YMM8R4_OMP_DESCRIPTION   = " Computational ElectroMagnetics related helper routines OpenMP-multithreaded."
                       

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"



namespace gms {



          namespace radiolocation {
          
                   __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__  
                   void sdotv_ymm8r4_unroll6x_omp( const __m256 * __restrict __ATTR_ALIGN__(32) pv1x,
	                                           const __m256 * __restrict __ATTR_ALIGN__(32) pv1y,
	                                           const __m256 * __restrict __ATTR_ALIGN__(32) pv1z,
	                                           const __m256 * __restrict __ATTR_ALIGN__(32) pv2x,
	                                           const __m256 * __restrict __ATTR_ALIGN__(32) pv2y,
	                                           const __m256 * __restrict __ATTR_ALIGN__(32) pv2z,
	                                           __m256 * __restrict __ATTR_ALIGN__(32) pdtv,
	                                           const int32_t n);
	                                           
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                 
	           void cdotv_ymm8c4_unroll10x_omp(const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                
	           void cdotv_ymm8c4_unroll6x_omp(const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST); 
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                 
	           void cdotv_ymm8c4_unroll2x_omp(  const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                 
	           void cdotv_ymm8c4_rolled_omp(    const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv2z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	                                               
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                
	           void cnorm_ymm8c4_unroll10x_omp(const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                   
	           void cnorm_ymm8c4_unroll6x_omp(const ymm8c4_t   * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
	           __ATTR_OPTIMIZE_03__ 
	           __ATTR_NO_STACK_PROTECTOR__                                  
	           void cnorm_ymm8c4_unroll2x_omp(  const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
	           __ATTR_OPTIMIZE_03__ 
	           __ATTR_NO_STACK_PROTECTOR__                                   
	           void cnorm_ymm8c4_rolled_omp(    const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1x,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1y,
	                                            const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pv1z,
	                                            ymm8c4_t * __restrict __ATTR_ALIGN__(32) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                     
	           void scrossc_ymm8r4_unroll10x_omp( const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                      
	           void scrossc_ymm8r4_unroll6x_omp(  const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)  
	           __ATTR_OPTIMIZE_03__
	           __ATTR_NO_STACK_PROTECTOR__                                    
	           void scrossc_ymm8r4_unroll2x_omp(  const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32) 
	           __ATTR_OPTIMIZE_03__ 
	           __ATTR_NO_STACK_PROTECTOR__                                    
	           void scrossc_ymm8r4_rolled_omp(    const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1y, 
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const ymm8c4_t  * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presx,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presy,
	                                              ymm8c4_t * __restrict __ATTR_ALIGN__(32) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                      
	          void scrossv_ymm8r4_unroll10x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv1y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcx,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcy,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                      
	          void scrossv_ymm8r4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m256 * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m256 * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m256 * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m256 * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                      
	          void scrossv_ymm8r4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv1y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcx,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcy,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                        
	          void scrossv_ymm8r4_rolled_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pv1x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv1y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv1z,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2x,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2y,
	                                              const __m256 * __restrict __ATTR_ALIGN__(32) pv2z,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcx,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcy,
	                                              __m256 * __restrict __ATTR_ALIGN__(32) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                     
	          void dir_vec_ymm8r4_unroll10x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                 
	          void dir_vec_ymm8r4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                 
	          void dir_vec_ymm8r4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                   
	          void dir_vec_ymm8r4_rolled_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvx,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvy,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                
	          void pol_vec_ymm8r4_unroll10x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                          const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                          const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) ppvx,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) ppvy,
	                                          __m256 * __restrict __ATTR_ALIGN__(32) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                              
	          void pol_vec_ymm8r4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m256 * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m256 * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m256 * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m256 * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m256 * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                  
	          void pol_vec_ymm8r4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                            const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                            const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                            __m256 * __restrict __ATTR_ALIGN__(32) ppvx,
	                                            __m256 * __restrict __ATTR_ALIGN__(32) ppvy,
	                                            __m256 * __restrict __ATTR_ALIGN__(32) ppvz,
	                                            const int32_t n,
	                                            int32_t & PF_DIST); 
	                                            
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                   
	          void pol_vec_ymm8r4_rolled_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                            const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                            const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                            __m256 * __restrict __ATTR_ALIGN__(32) ppvx,
	                                            __m256 * __restrict __ATTR_ALIGN__(32) ppvy,
	                                            __m256 * __restrict __ATTR_ALIGN__(32) ppvz,
	                                            const int32_t n,
	                                            int32_t & PF_DIST); 
	                                            
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                   
	          void H_XYZ_VP_ymm8c4_unroll10x_omp( const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                      
	          void H_XYZ_VP_ymm8c4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                      
	          void H_XYZ_VP_ymm8c4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                       
	          void H_XYZ_VP_ymm8c4_rolled_omp(    const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                       
	          void B_XYZ_VP_ymm8c4_unroll10x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                       
	          void B_XYZ_VP_ymm8c4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	                                                 
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                        
	          void B_XYZ_VP_ymm8c4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	           __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                          
	          void B_XYZ_VP_ymm8c4_rolled_omp(const __m256 * __restrict __ATTR_ALIGN__(32) pvpolx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpoly,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvpolz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdiry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvdirz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrx,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvry,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pvrz,
	                                               const __m256 * __restrict __ATTR_ALIGN__(32) pomega,
	                                               const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pk,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                               ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                    
	          void B_XYZ_H_XYZ_P_ymm8c4_unroll10x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                         
	          void B_XYZ_H_XYZ_P_ymm8c4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);  
	                                                
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                       
	          void B_XYZ_H_XYZ_P_ymm8c4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                       
	          void B_XYZ_H_XYZ_P_ymm8c4_rolled_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppsi,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppx,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppy,
	                                                const __m256 * __restrict __ATTR_ALIGN__(32) ppz,
	                                                const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pr,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pH_z,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_x,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_y,
	                                                ymm8c4_t * __restrict __ATTR_ALIGN__(32) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	                                                
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                       
	          void B_XYZ_H_XYZ_EP_ymm8c4_unroll10x_omp( const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                              
	          void B_XYZ_H_XYZ_EP_ymm8c4_unroll6x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	          __ATTR_HOT__
	          __ATTR_ALIGN__(32)  
	          __ATTR_OPTIMIZE_03__
	          __ATTR_NO_STACK_PROTECTOR__                                             
	          void B_XYZ_H_XYZ_EP_ymm8c4_unroll2x_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                      const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	        __ATTR_HOT__
	        __ATTR_ALIGN__(32)  
	        __ATTR_OPTIMIZE_03__
	        __ATTR_NO_STACK_PROTECTOR__                                             
	        void B_XYZ_H_XYZ_EP_ymm8c4_rolled_omp(const __m256 * __restrict __ATTR_ALIGN__(32) ptht,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pphi,
	                                                     const __m256 * __restrict __ATTR_ALIGN__(32) pomg,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) pphase,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) prefi,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppx,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppy,
	                                                     const ymm8c4_t * __restrict __ATTR_ALIGN__(32) ppz,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pH_z,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_x,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_y,
	                                                     ymm8c4_t * __restrict __ATTR_ALIGN__(32)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);                                        
          
          
         } // radiolocation  
          
          
          
          
          
          
          
          
          
          
          
          
} // gms          
          
          
          
          
          
          
          
#endif /*__GMS_EM_FIELDS_YMM8R4_OMP_H__*/
