

#ifndef __GMS_MEMOPS_H__
#define __GMS_MEMOPS_H__


namespace file_info {

   const unsigned int gGMS_MEMOPS_MAJOR = 1U;
   const unsigned int gGMS_MEMOPS_MINOR = 0U;
   const unsigned int gGMS_MEMOPS_MICRO = 0U;
   const unsigned int gGMS_MEMOPS_FULLVER =
         1000U*gGMS_MEMOPS_MAJOR+100U*gGMS_MEMOPS_MINOR+10U*gGMS_MEMOPS_MICRO;
   const char * const pgGMS_MEMOPS_CREATE_DATE = "26-09-2020 3:47PM +00200 (SAT 26 SEP 2020 GMT+2)";
   const char * const pgGMS_MEMOPS_BUILD_DATE  = __DATE__ ":" __TIME__;
   const char * const pgGMS_MEMOPS_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
}


#include <immintrin.h>
#include <cstdint>
#if defined __GNUC__ && !defined __INTEL_COMPILER
   #include <omp.h>
#endif
#include "GMS_config.h"
#include "GMS_avxvecf32.h"
#include "GMS_avxc8f32.h"
#include "GMS_avx512vec16.h"


namespace gms {

         namespace common {

	         __ATTR_ALWAYS_INLINE__
                 __ATTR_HOT__
		 __ATTR_ALIGN__(32)
		 static inline
                 void init_unroll2x_cmplxr4(std::complex<float> * __restrict __ATTR_ALIGN__(64) vc4,
			                    const int32_t len,
			                    const std::complex<float> val) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  vc4 = (std::complex<float>*)__builtin_assume_aligned(vc4,64);
#pragma omp simd aligned(vc4:64)
#elif defined __INTEL_COMPILER
                  __assume_aligned(vc4,64);
#pragma vector always		  
#endif
                 for(int32_t i = 0; i != len-1; i += 2) {
                       vc4[i+0] = val;
	               vc4[i+1] = val;
                 }
            }

	    __ATTR_ALWAYS_INLINE__
	    __ATTR_HOT__
	    __ATTR_ALIGN__(32)
	    static inline
	    void init_unroll4x_cmplxr4(std::complex<float> * __restrict __ATTR_ALIGN__(64) vc4,
			               const int32_t len,
			               const std::complex<float> val) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                 vc4 = (std::complex<float>*)__builtin_assume_aligned(vc4,64);
#pragma omp simd aligned(vc4:64)
#elif defined __INTEL_COMPILER
                 __assume_aligned(vc4,64);
#pragma vector always
#endif
              for(int32_t i = 0; i != len-3; i += 4) {
                  vc4[i+0] = val;
	          vc4[i+1] = val;
	          vc4[i+2] = val;
	          vc4[i+3] = val;
              }
         }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void init_unroll8x_cmplxr4(std::complex<float> * __restrict __ATTR_ALIGN__(64) vc4,
			            const int32_t len,
			            const std::complex<float> val) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                vc4 = (std::complex<float>*)__builtin_assume_aligned(vc4,64);
#pragma omp simd aligned(vc4:64)
#elif defined __INTEL_COMPILER
              __assume_aligned(vc4,64);
#pragma vector always
#endif
            for(int32_t i = 0; i != len-7; i += 8) {
                vc4[i+0LL] = val;
	        vc4[i+1LL] = val;
	        vc4[i+2LL] = val;
	        vc4[i+3LL] = val;
	        vc4[i+4LL] = val;
	        vc4[i+5LL] = val;
	        vc4[i+6LL] = val;
	        vc4[i+7LL] = val;
             }
         }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_init_unroll2x(AVXVec8 * __restrict __ATTR_ALIGN__(64) vec8,
		                    const int32_t len,
		                    const AVXVec8 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,64);
#elif defined __INTEL_COMPILER
       assume_aligned(vec8,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec8:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
               for(int32_t i = 0; i != len-1; i += 2) {
                    vec8[i+0] = v;
	            vec8[i+1] = v;
               }
           }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_init_unroll4x(AVXVec8 * __restrict __ATTR_ALIGN__(64) vec8,
                                    const int32_t len,
		                    const AVXVec8 v) {
#if defined __GNUC__
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,64);
#elif defined __INTEL_COMPILER
       assume_aligned(vec8,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec8:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
               for(int32_t i = 0; i != len-3; i += 4) {
                   vec8[i+0] = v;
	           vec8[i+1] = v;
	           vec8[i+2] = v;
	           vec8[i+3] = v;
               }
          }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_init_unroll8x(AVXVec8 * __restrict __ATTR_ALIGN__(64) vec8,
                                    const int32_t len,
		                    const AVXVec8 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,64);
#elif defined __INTEL_COMPILER
       __assume_aligned(vec8,64)
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec8:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
                for(int32_t i = 0; i != len-7; i += 8) {
                     vec8[i+0LL] = v;
	             vec8[i+1LL] = v;
	             vec8[i+2LL] = v;
	             vec8[i+3LL] = v;
	             vec8[i+4LL] = v;
	             vec8[i+5LL] = v;
	             vec8[i+6LL] = v;
	             vec8[i+7LL] = v;
	       } 
          }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_copy_unroll2x(AVXVec8 * __restrict __ATTR_ALIGN__(64) dst,
		                    const AVXVec8 * __restrict __ATTR_ALIGN__(64) src,
		                    const int32_t len) {
   
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
                  src = (const AVXVec8*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
                  __assume_aligned(dst,64);
                  __assume_aligned(src,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(src,dst:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
                for(int32_t i = 0; i != len-1; i += 2) {
                    dst[i+0] = src[i+0];
                    dst[i+1] = src[i+1];
               }
          }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_copy_unroll4x(AVXVec8 * __restrict __ATTR_ALIGN__(64) dst,
		                    const AVXVec8 * __restrict __ATTR_ALIGN__(64) src,
		                    const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                   dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
                   src = (const AVXVec8*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(dst,src:64);
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
                for(int32_t i = 0; i != len-3; i += 4) {
                    dst[i+0] = src[i+0];
                    dst[i+1] = src[i+1];
                    dst[i+2] = src[i+2];
                    dst[i+3] = src[i+3];
                }
            }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_copy_unroll8x(AVXVec8 * __restrict __ATTR_ALIGN__(64) dst,
		                    const AVXVec8 * __restrict __ATTR_ALIGN__(64) src,
		                    const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                    dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
                    src = (const AVXVec8*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(dst,src:64);
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
               for(int32_t i = 0; i != len-7; i += 8) {
                   dst[i+0] = src[i+0];
                   dst[i+1] = src[i+1];
                   dst[i+2] = src[i+2];
                   dst[i+3] = src[i+3];
	           dst[i+4] = src[i+4];
	           dst[i+5] = src[i+5];
	           dst[i+6] = src[i+6];
	           dst[i+7] = src[i+7];
              }
        }
	  
    }

}



#endif /*__GMS_MEMOPS_H__*/
