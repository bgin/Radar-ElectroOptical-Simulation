

#ifndef __GMS_RCS_CYLINDER_ZMM16R4_OMP_H__
#define __GMS_RCS_CYLINDER_ZMM16R4_OMP_H__

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

    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_MAJOR = 1U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_MINOR = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_MICRO = 0U;
    const unsigned int GMS_RCS_CYLINDER_ZMM16R4_OMP_FULLVER =
      1000U*GMS_RCS_CYLINDER_ZMM16R4_OMP_MAJOR+
      100U*GMS_RCS_CYLINDER_ZMM16R4_OMP_MINOR+
      10U*GMS_RCS_CYLINDER_ZMM16R4_OMP_MICRO;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_CREATION_DATE = "20-01-2023 16:36 PM +00200 (FRI 20 JAN 2023 GMT+2)";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_CYLINDER_ZMM16R4_OMP_DESCRIPTION   = "AVX512 optimized Cylinder Radar Cross Section (analytic) functionality OpenMP accelerated.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"



namespace gms {


           namespace radiolocation {
           
                        
                  
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                  
	           void rcs_f419_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST);
                                                       
              
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                  	        
                   void rcs_f419_zmm16r4_unroll10x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST); 
              
              
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
               	        
                   void rcs_f419_zmm16r4_unroll6x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST); 
                                                       
              
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   void rcs_f4120_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST); 
              
              
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   
                   void rcs_f4120_zmm16r4_unroll10x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST); 
              
              
              
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                 
                   void rcs_f4120_zmm16r4_unroll6x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                   const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                   __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                   const int32_t n,
                                                   int32_t & PF_DIST); 
                                                   
              
              
                 
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                 
	           void rcs_f4122_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pphi,
	                                            const __m512 * __restrict __ATTR_ALIGN__(64) pa,
	                                            const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
	                                            __m512 * __restrict __ATTR_ALIGN__(64) prcs,
	                                            const int32_t n,
	                                            int32_t & PF_DIST); 
	       
	       
	        
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                  
	           void rcs_f4122_zmm16r4_unroll10x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pphi,
	                                                const __m512 * __restrict __ATTR_ALIGN__(64) pa,
	                                                const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
	                                                __m512 * __restrict __ATTR_ALIGN__(64) prcs,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	       
	       
	         
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                 
	           void rcs_f4122_zmm16r4_unroll6x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pphi,
	                                                const __m512 * __restrict __ATTR_ALIGN__(64) pa,
	                                                const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
	                                                __m512 * __restrict __ATTR_ALIGN__(64) prcs,
	                                                const int32_t n,
	                                                int32_t & PF_DIST); 
	                                                
	       
	       
	          
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   
                   void rcs_f4124_zmm16r4_unroll16x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST); 
              
              
                   
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                 
                   void rcs_f4124_zmm16r4_unroll10x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST); 
                                                       
              
                
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                 
                   void rcs_f4124_zmm16r4_unroll6x_omp(const __m512 * __restrict __ATTR_ALIGN__(64) pa,
                                                       const __m512 * __restrict __ATTR_ALIGN__(64) pk0a,
                                                       __m512 * __restrict __ATTR_ALIGN__(64) prcs,
                                                       const int32_t n,
                                                       int32_t & PF_DIST); 
              
              
              
                
              
              
              
              
              
	       
	       
	       
	       
	       
              
              
              
              
              
           
         } // radiolocation

} // gms










#endif /*__GMS_RCS_CYLINDER_ZMM16R4_OMP_H__*/
