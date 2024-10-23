

#ifndef __GMS_SIMD_MEMOPS_H__
#define __GMS_SIMD_MEMOPS_H__


namespace file_info {

   const unsigned int gGMS_SIMD_MEMOPS_MAJOR = 1U;
   const unsigned int gGMS_SIMD_MEMOPS_MINOR = 0U;
   const unsigned int gGMS_SIMD_MEMOPS_MICRO = 0U;
   const unsigned int gGMS_SIMD_MEMOPS_FULLVER =
         1000U*gGMS_SIMD_MEMOPS_MAJOR+100U*gGMS_SIMD_MEMOPS_MINOR+10U*gGMS_SIMD_MEMOPS_MICRO;
   const char * const pgGMS_SIMD_MEMOPS_CREATE_DATE = "26-09-2020 3:47PM +00200 (SAT 26 SEP 2020 GMT+2)";
   const char * const pgGMS_SIMD_MEMOPS_BUILD_DATE  = __DATE__ ":" __TIME__;
   const char * const pgGMS_SIMD_MEMOPS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"

namespace gms {

         namespace common {

	         __ATTR_HOT__
	         __ATTR_ALIGN__(32)
                 void init_unroll2x_cmplxr4(std::complex<float> * __restrict  vc4,
			                    const int32_t len,
			                    const std::complex<float> val); 

	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    void init_unroll4x_cmplxr4(std::complex<float> * __restrict  vc4,
			               const int32_t len,
			               const std::complex<float> val); 

	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 void init_unroll8x_cmplxr4(std::complex<float> * __restrict  vc4,
			            const int32_t len,
			            const std::complex<float> val);
	
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 void avxvec8_init_unroll2x(AVXVec8 * __restrict  vec8,
		                    const int32_t len,
		                    const AVXVec8 v); 
		                    
          __ATTR_HOT__
	 __ATTR_ALIGN__(32)                  
	 void avxvec8_init_unroll4x(AVXVec8 * __restrict  vec8,
                                    const int32_t len,
		                    const AVXVec8 v); 
		                    
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                    
	 void avxvec8_init_unroll8x(AVXVec8 * __restrict  vec8,
                                    const int32_t len,
		                    const AVXVec8 v); 
		                    
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                    
	 void avxvec8_copy_unroll2x(AVXVec8 * __restrict  dst,
		                    const AVXVec8 * __restrict  src,
		                    const int32_t len); 
		                    
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                    
	 void avxvec8_copy_unroll4x(AVXVec8 * __restrict  dst,
		                    const AVXVec8 * __restrict  src,
		                    const int32_t len);
		                    
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                    
	 void avxvec8_copy_unroll8x(AVXVec8 * __restrict  dst,
		                    const AVXVec8 * __restrict  src,
		                    const int32_t len); 
		                    

          __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 void avx512vec16_init_unroll2x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) vec16,
			                const int32_t len,
			                const AVX512Vec16 v);
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avx512vec16_init_unroll4x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) vec16,
			                const int32_t len,
			                const AVX512Vec16 v);
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avx512vec16_init_unroll8x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) vec16,
			                const int32_t len,
			                const AVX512Vec16 v); 
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avx512vec16_copy_unroll2x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const AVX512Vec16 * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len); 
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avx512vec16_copy_unroll4x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const AVX512Vec16 * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len); 
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avx512vec16_copy_unroll8x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const AVX512Vec16 * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len); 
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avxvec8_copy_from_r4(AVXVec8 * __restrict  dst,
		                   const float * __restrict  src,
		                   const int32_t len); 
		                   
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                   
	 void r4_copy_from_avxvec8(float * __restrict _ATTR_ALIGN__(32) dst,
                                   const AVXVec8 * __restrict __ATTR_ALIGN__(64) src,
		                   const int32_t len); 
		                   
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                   
	 void r4_copy_from_avx512vec16(float * __restrict __ATTR_ALIGN__(64) dst,
                                       const AVX512Vec16 * __restrict _ATTTR_ALIGN__(64) src,
			               const int32_t len); 
			               
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		               
	 void  avx512vec16_copy_from_r4(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const float * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len); 
			                
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                
	 void avxc8f32_copy_from_r4(AVXc8f32 * __restrict  dst,
                                    const float * __restrict  src_re,
		                    const float * __restrict _ATTR_ALIGN__(32) src_im,
		                    const int32_t len); 
		                    
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                    
	 void r4_copy_from_avxc8f32(float * __restrict  dst_re,
                                    float * __restrict  dst_im,
		                    const AVXc8f32 * __restrict _ATTR_ALIGN__(32) src,
		                    const int32_t len); 
		                    
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                    
	 void avx256_init_unroll2x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                               double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       double * __restrict __ATTR_ALIGN__(64),
#endif
		                       const int32_t vlen,
		                       const double val); 
		                       
	  __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                       
	 void avx256_init_unroll2x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                               float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       float * __restrict  v,
#endif
			               const int32_t vlen,
			               const float val); 
			               
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		               
	 void avx256_init_unroll4x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict  v,
#endif
		                      const int32_t vlen,
		                      const double val); 
		                      
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                      
	 void avx256_init_unroll4x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict  v,
#endif
			              const int32_t vlen,
			              const float val); 
			              
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		              
	 void avx256_init_unroll8x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict  v,
#endif
			              const int32_t vlen,
			              const double val); 
			              
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)			              
	 void avx256_init_unroll8x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                             float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                     float * __restrict  v,
#endif
			             const int32_t vlen,
			             const float val); 


	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)	
	 void avx512_init_unroll2x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(64) v,
#endif
			             const int32_t vlen,
			             const double val); 
			             
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			             
	 void avx512_init_unroll2x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val);
			              
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)			              
	 void avx512_init_unroll4x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const double val); 
			              
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			              
	 void avx512_init_unroll4x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val); 
			              
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			              
	 void avx512_init_unroll8x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(64) v,
#endif
		                      const int32_t vlen,
		                      const double val); 
		                      
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)		                      
	 void avx512_init_unroll8x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float  __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val); 
			              
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			              
	 void avx512_init_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val); 
			              
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			              
	 void avx256_init_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict  v,
#endif
			              const int32_t vlen,
			              const float val); 

#include <complex.h>


       


	 
         __ATTR_HOT__
	 __ATTR_ALIGN__(32)	
	 void  avx256_memcpy_ps(
#if defined __ICC || defined __INTEL_COMPILER
                                   float * __restrict dst,
				   const float * __restrict src,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                   float * __restrict  dst,
				   const float * __restrict  src,
#endif
                                   const int32_t len); 
                                   
                                   
          __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                         
	 void  avx512_memcpy_ps(
#if defined __ICC || defined __INTEL_COMPILER
                                   float * __restrict dst,
				   const float * __restrict src,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                   float * __restrict __ATTR_ALIGN__(64) dst,
				   const float * __restrict __ATTR_ALIGN__(64) src,
#endif
                                   const int32_t len);
                                   
         __ATTR_HOT__
	 __ATTR_ALIGN__(32)	                          
	 void avx256_cached_memmove(void * __restrict _Dst,
			            const void * __restrict _Src,
			            const int32_t _nelems); 
			            
			            
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			            
	 void avx256_uncached_memmove(void * __restrict _Dst,
			              const void * __restrict _Src,
			              const int32_t _nelems); 


        __ATTR_HOT__
	 __ATTR_ALIGN__(32)	
	 void avx512_cached_memmove(void * __restrict _Dst,
			            const void * __restrict _Src,
			            const int32_t _nelems); 
			            
			            
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			            
	 void avx512_cached_memmove_a(void * __restrict _Dst,
			            const void * __restrict _Src,
			            const int32_t _nelems); 
			            
			            
	__ATTR_HOT__
	 __ATTR_ALIGN__(32)			            
	 void avx512_uncached_memmove(void * __restrict _Dst,
			             const void * __restrict _Src,
			             const int32_t _nelems); 
			             
			             
        __ATTR_HOT__
	 __ATTR_ALIGN__(32)		             
         static inline
	 void avx512_uncached_memmove_a(void * __restrict _Dst,
			             const void * __restrict _Src,
			             const int32_t _nelems); 


            
			 
    } // common

} // gms



#endif /*__GMS_SIMD_MEMOPS_H__*/
