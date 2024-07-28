

#ifndef __GMS_EM_FIELDS_ZMM8R8_OMP_H__
#define __GMS_EM_FIELDS_ZMM8R8_OMP_H__ 110620230936

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

    const unsigned int GMS_EM_FIELDS_ZMM8R8_OMP_MAJOR = 1U;
    const unsigned int GMS_EM_FIELDS_ZMM8R8_OMP_MINOR = 0U;
    const unsigned int GMS_EM_FIELDS_ZMM8R8_OMP_MICRO = 0U;
    const unsigned int GMS_EM_FIELDS_ZMM8R8_OMP_FULLVER =
      1000U*GMS_EM_FIELDS_ZMM8R8_OMP_MAJOR+
      100U*GMS_EM_FIELDS_ZMM8R8_OMP_MINOR+
      10U*GMS_EM_FIELDS_ZMM8R8_OMP_MICRO;
    const char * const GMS_EM_FIELDS_ZMM8R8_OMP_CREATION_DATE = "11-06-2023 09:36 AM +00200 (SUN 11 06 2023 GMT+2)";
    const char * const GMS_EM_FIELDS_ZMM8R8_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_EM_FIELDS_ZMM8R8_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_EM_FIELDS_ZMM8R8_OMP_DESCRIPTION   = " Computational ElectroMagnetics related helper routines OpenMP-multithreaded."
                       

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
                    void sdotv_zmm8r8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                           const __m512d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                           const __m512d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                           const __m512d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                           const __m512d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                           const __m512d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                           __m512d * __restrict __ATTR_ALIGN__(64) pdtv,
	                                           const int32_t n);
	                                           
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                               
	            void cdotv_zmm8c8_unroll10x_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                  
	            void cdotv_zmm8c8_unroll6x_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                 
	            void cdotv_zmm8c8_unroll2x_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	             __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                
	            void cdotv_zmm8c8_rolled_omp(  const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv2z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pres,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                 
	            void cnorm_zmm8c8_unroll10x_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void cnorm_zmm8c8_unroll6x_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                               
	            void cnorm_zmm8c8_unroll2x_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                
	            void cnorm_zmm8c8_rolled_omp(const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1x,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1y,
	                                            const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pv1z,
	                                            zmm8c8_t * __restrict __ATTR_ALIGN__(64) pvs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST);
	                                            
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void scrossc_zmm8r8_unroll10x_omp(const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1y, 
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presx,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presy,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                    
	            void scrossc_zmm8r8_unroll6x_omp(const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1y, 
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presx,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presy,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void scrossc_zmm8r8_unroll2x_omp(const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1y, 
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presx,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presy,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                    
	            void scrossc_zmm8r8_rolled_omp(const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1y, 
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const zmm8c8_t  * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presx,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presy,
	                                              zmm8c8_t * __restrict __ATTR_ALIGN__(64) presz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	             __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void scrossv_zmm8r8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                     
	            void scrossv_zmm8r8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void scrossv_zmm8r8_unroll2x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void scrossv_zmm8r8_rolled_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pv1x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv1z,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2x,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2y,
	                                              const __m512d * __restrict __ATTR_ALIGN__(64) pv2z,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcx,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcy,
	                                              __m512d * __restrict __ATTR_ALIGN__(64) pvcz,
	                                              const int32_t n,
	                                              int32_t & PF_DIST);
	                                              
	             __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                     
	            void dir_vec_zmm8r8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	             __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                  
	            void dir_vec_zmm8r8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                               
	            void dir_vec_zmm8r8_unroll2x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                
	            void dir_vec_zmm8r8_rolled_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) pdvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	                                             
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                
	            void pol_vec_zmm8r8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                               
	            void pol_vec_zmm8r8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void pol_vec_zmm8r8_unroll2x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                              
	            void pol_vec_zmm8r8_rolled_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                          const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvx,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvy,
	                                          __m512d * __restrict __ATTR_ALIGN__(64) ppvz,
	                                          const int32_t n,
	                                          int32_t & PF_DIST);
	                                          
	                                           
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                
	            void H_XYZ_VP_zmm8c8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void H_XYZ_VP_zmm8c8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                      
	            void H_XYZ_VP_zmm8c8_rolled_omp(    const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                               const int32_t n,
	                                               int32_t & PF_DISPATCH);
	                                               
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                       
	            void B_XYZ_VP_zmm8c8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pomega,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void B_XYZ_VP_zmm8c8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pomega,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	                                                 
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                      
	            void B_XYZ_VP_zmm8c8_unroll2x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pomega,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                      
	            void B_XYZ_VP_zmm8c8_rolled_omp(    const __m512d * __restrict __ATTR_ALIGN__(64) pvpolx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpoly,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvpolz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdiry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvdirz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrx,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvry,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pvrz,
	                                               const __m512d * __restrict __ATTR_ALIGN__(64) pomega,
	                                               const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pk,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                               zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                               const int32_t n,
	                                               int32_t & PF_DIST);
	                                               
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                    
	            void B_XYZ_H_XYZ_P_zmm8c8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppx,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppy,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppz,
	                                                const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pr,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	               __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                        
	            void B_XYZ_H_XYZ_P_zmm8c8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppx,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppy,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppz,
	                                                const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pr,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	                                                
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                          
	            void B_XYZ_H_XYZ_P_zmm8c8_unroll2x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppx,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppy,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppz,
	                                                const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pr,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	             __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                        
	            void B_XYZ_H_XYZ_P_zmm8c8_rolled_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppsi,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppx,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppy,
	                                                const __m512d * __restrict __ATTR_ALIGN__(64) ppz,
	                                                const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pr,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pH_z,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_x,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_y,
	                                                zmm8c8_t * __restrict __ATTR_ALIGN__(64) pB_z,
	                                                const int32_t n,
	                                                int32_t & PF_DIST);
	                                                
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                   
	            void B_XYZ_H_XYZ_EP_zmm8c8_unroll10x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                     const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                      const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pphase,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) prefi,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppx,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppy,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppz,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_z,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	            __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                          
	            void B_XYZ_H_XYZ_EP_zmm8c8_unroll6x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                     const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                      const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pphase,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) prefi,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppx,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppy,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppz,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_z,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	             __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                          
	            void B_XYZ_H_XYZ_EP_zmm8c8_unroll2x_omp(const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                     const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                      const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pphase,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) prefi,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppx,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppy,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppz,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_z,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
	                                                     
	              __ATTR_HOT__
	            __ATTR_ALIGN__(32)  
	            __ATTR_OPTIMIZE_03__
	            __ATTR_NO_STACK_PROTECTOR__                                         
	            void B_XYZ_H_XYZ_EP_zmm8c8_rolled_omp(    const __m512d * __restrict __ATTR_ALIGN__(64) ptht,
	                                                     const __m512d * __restrict __ATTR_ALIGN__(64) pphi,
	                                                     const __m512d * __restrict __ATTR_ALIGN__(64) pomg,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) pphase,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) prefi,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppx,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppy,
	                                                     const zmm8c8_t * __restrict __ATTR_ALIGN__(64) ppz,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pH_z,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_x,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_y,
	                                                     zmm8c8_t * __restrict __ATTR_ALIGN__(64)  pB_z,
	                                                     const int32_t n,
	                                                     int32_t & PF_DIST);
          
          
          
         } // radiolocation  
        
        
        
        
        
} // gms
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
#endif /*__GMS_EM_FIELDS_ZMM8R8_OMP_H__*/
