

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
#if defined __GNUC__ && !defined __INTEL_COMPILER
   #include <omp.h>
#endif
#include "GMS_config.h"
#include "GMS_avxvecf32.h"
#include "GMS_avxc8f32.h"
#include "GMS_avx512vec16.h"

#if !defined (MEMMOVE_1ELEM)
#define MEMMOVE_1ELEM 1
#endif

#if !defined (MEMMOVE_16ELEMS)
#define MEMMOVE_16ELEMS 16
#endif

#if !defined (MEMMOVE_32ELEMS)
#define MEMMOVE_32ELEMS 32
#endif

#if !defined (MEMMOVE_64ELEMS)
#define MEMMOVE_64ELEMS 64
#endif

#if !defined (MEMMOVE_128ELEMS)
#define MEMMOVE_128ELEMS 128
#endif

#if !defined (MEMMOVE_256ELEMS)
#define MEMMOVE_256ELEMS 256
#endif
// float type (4-bytes)
#define YMM_LEN (8)
#define ZMM_LEN (16)

#if !defined (PAGE4KiB)
#define PAGE4KiB 4096
#endif

#if !defined (MAXFLOATSPERPAGE4KiB)
#define MAXFLOATSPERPAGE4KiB 1024
#endif



#if !defined (min_val)
#define min_val(A,B) ((A)<(B)?(A):(B))
#endif

namespace gms {

         namespace common {

	         __ATTR_ALWAYS_INLINE__
                 __ATTR_HOT__
		 __ATTR_ALIGN__(32)
		 static inline
                 void init_unroll2x_cmplxr4(std::complex<float> * __restrict __ATTR_ALIGN__(32) vc4,
			                    const int32_t len,
			                    const std::complex<float> val) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  vc4 = (std::complex<float>*)__builtin_assume_aligned(vc4,32);
#pragma omp simd aligned(vc4:32)
#elif defined __INTEL_COMPILER
                  __assume_aligned(vc4,32);
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
	    void init_unroll4x_cmplxr4(std::complex<float> * __restrict __ATTR_ALIGN__(32) vc4,
			               const int32_t len,
			               const std::complex<float> val) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                 vc4 = (std::complex<float>*)__builtin_assume_aligned(vc4,32);
#pragma omp simd aligned(vc4:32)
#elif defined __INTEL_COMPILER
                 __assume_aligned(vc4,32);
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
	 void init_unroll8x_cmplxr4(std::complex<float> * __restrict __ATTR_ALIGN__(32) vc4,
			            const int32_t len,
			            const std::complex<float> val) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                vc4 = (std::complex<float>*)__builtin_assume_aligned(vc4,32);
#pragma omp simd aligned(vc4:32)
#elif defined __INTEL_COMPILER
              __assume_aligned(vc4,32);
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
	 void avxvec8_init_unroll2x(AVXVec8 * __restrict __ATTR_ALIGN__(32) vec8,
		                    const int32_t len,
		                    const AVXVec8 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,32);
#elif defined __INTEL_COMPILER
       assume_aligned(vec8,32);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec8:32)
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
	 void avxvec8_init_unroll4x(AVXVec8 * __restrict __ATTR_ALIGN__(32) vec8,
                                    const int32_t len,
		                    const AVXVec8 v) {
#if defined __GNUC__
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,32);
#elif defined __INTEL_COMPILER
       assume_aligned(vec8,32);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec8:32)
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
	 void avxvec8_init_unroll8x(AVXVec8 * __restrict __ATTR_ALIGN__(32) vec8,
                                    const int32_t len,
		                    const AVXVec8 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,32);
#elif defined __INTEL_COMPILER
       __assume_aligned(vec8,32)
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec8:32)
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
	 void avxvec8_copy_unroll2x(AVXVec8 * __restrict __ATTR_ALIGN__(32) dst,
		                    const AVXVec8 * __restrict __ATTR_ALIGN__(32) src,
		                    const int32_t len) {
   
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  dst = (AVXVec8*)__builtin_assume_aligned(dst,32);
                  src = (const AVXVec8*)__builtin_assume_aligned(src,32);
#elif defined __INTEL_COMPILER
                  __assume_aligned(dst,32);
                  __assume_aligned(src,32);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(src,dst:32)
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
	 void avxvec8_copy_unroll4x(AVXVec8 * __restrict __ATTR_ALIGN__(32) dst,
		                    const AVXVec8 * __restrict __ATTR_ALIGN__(32) src,
		                    const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                   dst = (AVXVec8*)__builtin_assume_aligned(dst,32);
                   src = (const AVXVec8*)__builtin_assume_aligned(src,32);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,32);
     __assume_aligned(src,32);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(dst,src:32);
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
	 void avxvec8_copy_unroll8x(AVXVec8 * __restrict __ATTR_ALIGN__(32) dst,
		                    const AVXVec8 * __restrict __ATTR_ALIGN__(32) src,
		                    const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                    dst = (AVXVec8*)__builtin_assume_aligned(dst,32);
                    src = (const AVXVec8*)__builtin_assume_aligned(src,32);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,32);
     __assume_aligned(src,32);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(dst,src:32);
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


#if defined __AVX512F__

       	 __ATTR_ALWAYS_INLINE__
	 __ATTR_VECTORCALL__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512vec16_init_unroll2x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) vec16,
			                const int32_t len,
			                const AVX512Vec16 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     vec16 = (AVX512Vec16*)__builtin_assume_aligned(vec16,64);
#elif defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(vec16,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec16:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
              for(int32_t i = 0; i != len-1; i += 2) {
                  vec16[i+0] = v;
	          vec16[i+1] = v;
              }
         }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_VECTORCALL__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512vec16_init_unroll4x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) vec16,
			                const int32_t len,
			                const AVX512Vec16 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                vec16 = (AVX512Vec16*)__builtin_assume_aligned(vec16,64);
#elif defined __ICC || defined __INTEL_COMPILER
               __assume_aligned(vec16,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec16:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
              for(int64_t i = 0; i != len-3; i += 4) {
                   vec16[i+0] = v;
	           vec16[i+1] = v;
	           vec16[i+2] = v;
	           vec16[i+3] = v;
              }
         }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_VECTORCALL__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512vec16_init_unroll8x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) vec16,
			                const int32_t len,
			                const AVX512Vec16 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                vec16 = (AVX512Vec16*)__builtin_assume_aligned(vec16,64);
#elif defined __ICC || defined __INTEL_COMPILER
            __assume_aligned(vec16,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(vec16:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#pragma code_align(32)
#endif
              for(int64_t i = 0; i != len-7; i += 8) {
                   vec16[i+0] = v;
	           vec16[i+1] = v;
	           vec16[i+2] = v;
	           vec16[i+3] = v;
	           vec16[i+4] = v;
	           vec16[i+5] = v;
	           vec16[i+6] = v;
	           vec16[i+7] = v;
                }
         }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512vec16_copy_unroll2x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const AVX512Vec16 * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
            dst = (AVXV512ec16*)__builtin_assume_aligned(dst,64);
            src = (const AVX512Vec16*)__builtin_assume_aligned(src,64);
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
                for(int32_t i = 0; i != len-1; i += 2) {
                    dst[i+0] = src[i+0];
	            dst[i+1] = src[i+1];
                }

          }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512vec16_copy_unroll4x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const AVX512Vec16 * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  dst = (AVXV512ec16*)__builtin_assume_aligned(dst,64);
                  src = (const AVX512Vec16*)__builtin_assume_aligned(src,64);
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
	 void avx512vec16_copy_unroll8x(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const AVX512Vec16 * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                      dst = (AVXV512ec16*)__builtin_assume_aligned(dst,64);
                      src = (const AVX512Vec16*)__builtin_assume_aligned(src,64);
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
                  for(int64_t i = 0; i != len-7; i += 8) {
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

	    

#endif

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxvec8_copy_from_r4(AVXVec8 * __restrict __ATTR_ALIGN__(32) dst,
		                   const float * __restrict __ATTR_ALIGN__(32) src,
		                   const int32_t len) {
     
#if defined __GNUC__ && !defined __INTEL_COMPILER
              dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
              src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
              __assume_aligned(dst,64);
              __assume_aligned(src,64);
#endif
             int32_t j;
             j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
               for(int32_t i = 0; i != len; i += 8) {
		   const __m256 t = _mm256_load_ps(&src[i]);
	           dst[j].m_v8 = t;
                   ++j;
                }

         }
          
	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void r4_copy_from_avxvec8(float * __restrict _ATTR_ALIGN__(32) dst,
                                   const AVXVec8 * __restrict __ATTR_ALIGN__(64) src,
		                   const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                 dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
                 src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
                __assume_aligned(dst,64);
                __assume_aligned(src,64);
#endif
               int32_t j;
               j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif    
              for(int32_t i = 0; i != len; i += 8) {
                  const __m256 t = src[j].m_v8;
	          _mm256_store_ps(&dst[i],t);
	          ++j;
              }

         }

#if defined __AVX512F__

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void r4_copy_from_avx512vec16(float * __restrict __ATTR_ALIGN__(64) dst,
                                       const AVX512Vec16 * __restrict _ATTTR_ALIGN__(64) src,
			               const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                  dst = (AVX512Vec16*)__builtin_assume_aligned(dst,64);
                  src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
                 __assume_aligned(dst,64);
                 __assume_aligned(src,64);
#endif
                int32_t j;
                j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif    
                 for(int32_t i = 0; i != len; i += 16) {
                      const __m512 t = src[j].m_v16;
	              _mm512_store_ps(&dst[i],t);
	              ++j;
                    }
               }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void  avx512vec16_copy_from_r4(AVX512Vec16 * __restrict __ATTR_ALIGN__(64) dst,
                                        const float * __restrict __ATTR_ALIGN__(64) src,
			                const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                     dst = (AVX512Vec16*)__builtin_assume_aligned(dst,64);
                     src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
                     __assume_aligned(dst,64);
                     __assume_aligned(src,64);
#endif
                     int32_t j;
                     j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif    
                  for(int32_t i = 0; i != len; i += 16){
		       const __m512 t = _mm512_load_ps(&src[i]);
	               dst[j].m_v16 = t;
 	               ++j;
                 }
             }

#endif

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avxc8f32_copy_from_r4(AVXc8f32 * __restrict __ATTR_ALIGN__(32) dst,
                                    const float * __restrict __ATTR_ALIGN__(32) src_re,
		                    const float * __restrict _ATTR_ALIGN__(32) src_im,
		                    const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                   dst    = (AVXc8f32*)__builtin_assume_aligned(dst,64);
                   src_re = (const float*)__builtin_assume_aligned(src_re,64);
                   src_im = (const float*)__builtin_assume_aligned(src_im,64);
#elif defined __ICC || defined __INTEL_COMPILER
                   __assume_aligned(dst,64);
                   __assume_aligned(src_re,64);
                   __assume_aligned(src_im,64);
#endif
                  int32_t j;
                  j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                for(int64_t i = 0LL; i != len; i += 8LL) {
                    const __m256 tre = _mm256_load_ps(&src_re[i]);
	            dst[j].m_re = tre;
	            const __m256 tim = _mm256_load_ps(&src_im[i]);
	            dst[j].m_im = tim;
	            ++j;
                 }
           }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void r4_copy_from_avxc8f32(float * __restrict __ATTR_ALIGN__(32) dst_re,
                                    float * __restrict __ATTR_ALIGN__(32) dst_im,
		                    const AVXc8f32 * __restrict _ATTR_ALIGN__(32) src,
		                    const int32_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
                     dst_re = (float*)__builtin_assume_aligned(dst_re,64);
                     dst_im = (float*)__builtin_assume_aligned(dst_im,64);
                     src    = (const AVXc8f32*)__builtin_assume_aligned(src,64);
#elif defined __ICC || defined __INTEL_COMPILER
                     __assume_aligned(dst_re,64);
                     __assume_aligned(dst_im,64);
                     __assume_aligned(src,64);
#endif
                    int32_t j;
                    j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                    for(int32_t i = 0; i != len; i += 8) {
                        const __m256 tre = src[j].m_re;
	                _mm256_store_ps(&dst_re[i],tre;
	                const __m256 tim = src[j].m_im;
	                _mm256_store_ps(&dst_im[i],tim);
                   }
               }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_unroll2x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                               double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       double * __restrict __ATTR_ALIGN__(64),
#endif
		                       const int32_t vlen,
		                       const double val) {
	
 	     __m256d vec = _mm256_set1_pd(val);
	     int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
              for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 8) {
			_mm256_storeu_pd(&v[i + 0], vec);
			_mm256_storeu_pd(&v[i + 4], vec);
	        }
#pragma loop_count min(1),avg(2),max(3)
               for(; i != vlen; ++i) {
		   v[i] = val;
	        }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
       
	if ((reinterpret_cast<uintptr_t>(v)& 0x1F) != 0ULL) {
		for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 8) {
			_mm256_storeu_pd(&v[i + 0], vec);
			_mm256_storeu_pd(&v[i + 4], vec);
		}
	}
	else {
	     v = (double*)__builtin_assume_aligned(v,32);
     	      for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 8) {
		    _mm256_store_pd(&v[i+0], vec);
		    _mm256_store_pd(&v[i+4], vec);
	      }
	  }
	     for (; i != vlen; ++i) {
		v[i] = val;
	     }
#endif
          }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_unroll2x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                               float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                       float * __restrict __ATTR_ALIGN__(32) v,
#endif
			               const int32_t vlen,
			               const float val) {
                __m256 ymm0 = _mm256_set1_ps(val);
	        int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 16) {
	             _mm256_storeu_ps(&v[i+0], ymm0);
	             _mm256_storeu_ps(&v[i+8], ymm0);
	         }
#pragma loop_count min(1),avg(4),max(7)
                 for(; i != vlen; ++i) {
	            v[i] = val;
	         }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        
	      if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {
	          for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 16) {
	               _mm256_storeu_ps(&v[i+0], ymm0);
	               _mm256_storeu_ps(&v[i+8], ymm0);
	          }
	      } else {
	          v = (float*)__builtin_assume_aligned(v,32);
                  for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 16) {
                      _mm256_store_ps(&v[i+0], ymm0);
	              _mm256_store_ps(&v[i+8], ymm0);
	          }
	      }
	       for(; i != vlen; ++i) {
	           v[i] = val;
	       }
#endif
            }

	  __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_unroll4x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(32) v,
#endif
		                      const int32_t vlen,
		                      const double val) {
	
	          __m256d vec = _mm256_set1_pd(val);
	          int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                        for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 16) {
		 	     _mm256_storeu_pd(&v[i + 0], vec);
			     _mm256_storeu_pd(&v[i + 4], vec);
			     _mm256_storeu_pd(&v[i + 8], vec);
			     _mm256_storeu_pd(&v[i + 12], vec);
		}
#pragma loop_count min(1),avg(2),max(3)
               	        for (; i != vlen; ++i) {
		             v[i] = val;
 	                }    
#elif defined __GNUC__ && !defined __INTEL_COMPILER
       
	                if ((reinterpret_cast<uintptr_t>(v)& 0x1F) != 0ULL) {
		             for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 16) {
		 	          _mm256_storeu_pd(&v[i + 0], vec);
			          _mm256_storeu_pd(&v[i + 4], vec);
			          _mm256_storeu_pd(&v[i + 8], vec);
			          _mm256_storeu_pd(&v[i + 12], vec);
		            }
	              }
 	               else {
	                      v = (double*)__builtin_assume_aligned(v,32);
	                     for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 16) {
		                  _mm256_store_pd(&v[i+0],  vec);
		                  _mm256_store_pd(&v[i+4],  vec);
		                  _mm256_store_pd(&v[i+8],  vec);
		                  _mm256_store_pd(&v[i+12], vec);
	                      }
	              }
 	                      for (; i != vlen; ++i) {
		                   v[i] = val;
	                       }
#endif
                    }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_unroll4x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIGN__(32) v,
#endif
			              const int32_t vlen,
			              const float val) {
                __m256 ymm0 = _mm256_set1_ps(val);
                int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                     for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 32) {
	                  _mm256_storeu_ps(&v[i+0],  ymm0);
	                  _mm256_storeu_ps(&v[i+8],  ymm0);
	                  _mm256_storeu_ps(&v[i+16], ymm0);
	                  _mm256_storeu_ps(&v[i+24], ymm0);
	            }
#pragma loop_count min(1),avg(4),max(7)
                    for(; i != vlen; ++i) {
	                v[i] = val;
                    }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     
                   if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {
	               for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 32) {
	                    _mm256_storeu_ps(&v[i+0],  ymm0);
	                    _mm256_storeu_ps(&v[i+8],  ymm0);
	                    _mm256_storeu_ps(&v[i+16], ymm0);
	                    _mm256_storeu_ps(&v[i+24], ymm0);
	               }
                   } else {
                        v = (float*)__builtin_assume_aligned(v,32);
	                for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 32) {
                             _mm256_store_ps(&v[i+0],  ymm0);
	                     _mm256_store_ps(&v[i+8],  ymm0);
	                     _mm256_store_ps(&v[i+16], ymm0);
	                     _mm256_store_ps(&v[i+24], ymm0);
	                }
                   }
                      for(; i != vlen; ++i) {
	                  v[i] = val;
                      }
#endif
              }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_unroll8x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(32) v,
#endif
			              const int32_t vlen,
			              const double val) {
	
                   __m256d vec = _mm256_set1_pd(val);
	           int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                           for(i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 32) {
			        _mm256_storeu_pd(&v[i + 0], vec);
			        _mm256_storeu_pd(&v[i + 4], vec);
			        _mm256_storeu_pd(&v[i + 8], vec);
			        _mm256_storeu_pd(&v[i + 12], vec);
			        _mm256_storeu_pd(&v[i + 16], vec);
			        _mm256_storeu_pd(&v[i + 20], vec);
			        _mm256_storeu_pd(&v[i + 24], vec);
			        _mm256_storeu_pd(&v[i + 28], vec);
		            }
#pragma loop_count min(1),avg(2),max(3)
	                  for (; i != vlen; ++i){
		               v[i] = val;
	                  }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
        
	                  if ((reinterpret_cast<uintptr_t>(v)& 0x1F) != 0ULL) {
		                for (i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 32) {
			               _mm256_storeu_pd(&v[i + 0], vec);
			               _mm256_storeu_pd(&v[i + 4], vec);
			               _mm256_storeu_pd(&v[i + 8], vec);
			               _mm256_storeu_pd(&v[i + 12], vec);
			               _mm256_storeu_pd(&v[i + 16], vec);
			               _mm256_storeu_pd(&v[i + 20], vec);
			               _mm256_storeu_pd(&v[i + 24], vec);
			               _mm256_storeu_pd(&v[i + 28], vec);
		                  }
	                  }
	                   else {
	                        v = (double*)__builtin_assume_aligned(v,32);
	                        for ( i = 0; i != ROUND_TO_FOUR(vlen, 4); i += 32) {
		                      _mm256_store_pd(&v[i+0], vec);
		                      _mm256_store_pd(&v[i+4], vec);
		                      _mm256_store_pd(&v[i+8], vec);
		                      _mm256_store_pd(&v[i+12],vec);
		                      _mm256_store_pd(&v[i+16],vec);
		                      _mm256_store_pd(&v[i+20],vec);
		                      _mm256_store_pd(&v[i+24],vec);
		                      _mm256_store_pd(&v[i+28],vec);
	                       }
	                  }
	                    for (; i != vlen; ++i){
		                 v[i] = val;
	                    }
#endif
                      }


	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_unroll8x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                             float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                     float * __restrict __ATTR_ALIGN__(32) v,
#endif
			             const int32_t vlen,
			             const float val) {
                      __m256 ymm0 = _mm256_set1_ps(val);
                      int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                      for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 64) {
	                    _mm256_storeu_ps(&v[i+0],  ymm0);
	                    _mm256_storeu_ps(&v[i+8],  ymm0);
	                    _mm256_storeu_ps(&v[i+16], ymm0);
	                    _mm256_storeu_ps(&v[i+24], ymm0);
	                    _mm256_storeu_ps(&v[i+32], ymm0);
	                    _mm256_storeu_ps(&v[i+40], ymm0);
	                    _mm256_storeu_ps(&v[i+48], ymm0);
	                    _mm256_storeu_ps(&v[i+56], ymm0);
	              }
#pragma loop_count min(1),avg(4),max(7)
                        for(; i != vlen; ++i) {
	                    v[i] = val;
                           }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
     
                 if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {
	               for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 64) {
	                    _mm256_storeu_ps(&v[i+0],  ymm0);
	                    _mm256_storeu_ps(&v[i+8],  ymm0);
	                    _mm256_storeu_ps(&v[i+16], ymm0);
	                    _mm256_storeu_ps(&v[i+24], ymm0);
	                    _mm256_storeu_ps(&v[i+32], ymm0);
	                    _mm256_storeu_ps(&v[i+40], ymm0);
	                    _mm256_storeu_ps(&v[i+48], ymm0);
	                    _mm256_storeu_ps(&v[i+56], ymm0);
	          }
                }else {
                        v = (float*)__builtin_assume_aligned(v,32)l
	                for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 64) {
                             _mm256_store_ps(&v[i+0],  ymm0);
	                     _mm256_store_ps(&v[i+8],  ymm0);
	                     _mm256_store_ps(&v[i+16], ymm0);
	                     _mm256_store_ps(&v[i+24], ymm0);
	                     _mm256_store_ps(&v[i+32], ymm0);
	                     _mm256_store_ps(&v[i+40], ymm0);
	                     _mm256_store_ps(&v[i+48], ymm0);
	                     _mm256_store_ps(&v[i+56], ymm0);
	                 }
                 }
                       for(; i != vlen; ++i) {
	                   v[i] = val;
                       }
#endif
           }

#if defined __AVX512F__

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_unroll2x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(64) v,
#endif
			             const int32_t vlen,
			             const double val) {
                     __m512d zmm0 = _mm512_set1_pd(val);
                     int32_t i;
#if defined ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                    for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 16) {
	                 _mm512_storeu_pd(&v[i+0], zmm0);
	                 _mm512_storeu_pd(&v[i+8], zmm0);
	            }
#pragma loop_count min(1),avg(4),max(7)
                    for(; i != vlen; ++i) {
	                v[i] = val;
                    }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
       
               if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	           for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 16) {
	                _mm512_storeu_pd(&v[i+0], zmm0);
	                _mm512_storeu_pd(&v[i+8], zmm0);
	            }
                } else {
                       v = (double*)__builtin_assume_aligned(v,64);
	              for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 16) {
                          _mm512_store_pd(&v[i+0], zmm0);
	                  _mm512_store_pd(&v[i+8], zmm0);
	               }
                }
                  for(; i != vlen; ++i) {
	              v[i] = val;
                  }
#endif
             }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_unroll2x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val) {
                  __m512 zmm0 = _mm512_set1_ps(val);
                  int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
                  for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 32) {
	                _mm512_storeu_ps(&v[i+0],  zmm0);
	                _mm512_storeu_ps(&v[i+16], zmm0);
	           }
#pragma loop_count min(1),avg(4),max(7)
                  for(; i != vlen; ++i) {
	              v[i] = val;
                  }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
      
             if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	           for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 32) {
	                 _mm512_storeu_ps(&v[i+0],  zmm0);
	                 _mm512_storeu_ps(&v[i+16], zmm0);
	           }
              } else {
                      v = (float*)__builtin_assume_aligned(v,64);
	           for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 32) {
                       _mm512_store_ps(&v[i+0],  zmm0);
	               _mm512_store_ps(&v[i+16], zmm0);
	           }
               }
                  for(; i != vlen; ++i) {
	              v[i] = val;
                  }
#endif
          }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_unroll4x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const double val) {
                  __m512d zmm0 = _mm512_set1_pd(val);
                  int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                  for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 32) {
	              _mm512_storeu_pd(&v[i+0],  zmm0);
	              _mm512_storeu_pd(&v[i+8],  zmm0);
	              _mm512_storeu_pd(&v[i+16], zmm0);
	              _mm512_storeu_pd(&v[i+24], zmm0);
	          }
#pragma loop_count min(1),avg(4),max(7)
                  for(; i != vlen; ++i) {
	              v[i] = val;
                  }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
       
                 if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	             for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 32) {
	                _mm512_storeu_pd(&v[i+0],  zmm0);
	                _mm512_storeu_pd(&v[i+8],  zmm0);
	                _mm512_storeu_pd(&v[i+16], zmm0);
	                _mm512_storeu_pd(&v[i+24], zmm0);
	             }
                  } else {
                        v = (double*)__builtin_assume_aligned(v,64);
                     for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 32) {
	                 _mm512_store_pd(&v[i+0],  zmm0);
	                 _mm512_store_pd(&v[i+8],  zmm0);
	                 _mm512_store_pd(&v[i+16], zmm0);
	                 _mm512_store_pd(&v[i+24], zmm0);
	              }
                  }
                     for(; i != vlen; ++i) {
	                 v[i] = val;
                     }
#endif
              }

       	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_unroll4x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val) {
                    __m512 zmm0 = _mm512_set1_ps(val);
                    int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                    for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 64) {
	                 _mm512_storeu_ps(&v[i+0], zmm0);
	                 _mm512_storeu_ps(&v[i+16], zmm0);
	                 _mm512_storeu_ps(&v[i+32], zmm0);
	                 _mm512_storeu_ps(&v[i+48], zmm0);
	            }
#pragma loop_count min(1),avg(8),max(15)
                    for(; i != vlen; ++i) {
	                v[i] = val;
                    }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                    if((reinterpret_cast<uintptr_t>(v) & 0x3F) ! = 0ULL) {
	                 for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 64LL) {
	                      _mm512_storeu_ps(&v[i+0LL], zmm0);
	                      _mm512_storeu_ps(&v[i+16LL], zmm0);
	                      _mm512_storeu_ps(&v[i+32LL], zmm0);
	                      _mm512_storeu_ps(&v[i+48LL], zmm0);
	                 }
                     } else {
                          v = (float*)__builtin_assume_aligned(v,64);
                          for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 64LL) {
	                        _mm512_store_ps(&v[i+0LL], zmm0);
	                        _mm512_store_ps(&v[i+16LL], zmm0);
	                        _mm512_store_ps(&v[i+32LL], zmm0);
	                        _mm512_store_ps(&v[i+48LL], zmm0);
	                   }
                     }
                          for(; i != vlen; ++i) {
	                      v[i] = val;
                          }
#endif
                 }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_unroll8x_pd(
#if defined __ICC || defined __INTEL_COMPILER
	                              double * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      double * __restrict __ATTR_ALIGN__(64) v,
#endif
		                      const int32_t vlen,
		                      const double val) {
	
                  __m512d vec = _mm512_set1_pd(val);
	          int32_t i;
#if defined __ICC || defined __INTEL_COMPILER

                  for (i = 0; i != ROUND_TO_EIGHT(vlen, 8); i += 64) {
			_mm512_storeu_pd(&v[i + 0], vec);
			_mm512_storeu_pd(&v[i + 8], vec);
			_mm512_storeu_pd(&v[i + 16], vec);
			_mm512_storeu_pd(&v[i + 24], vec);
			_mm512_storeu_pd(&v[i + 32], vec);
			_mm512_storeu_pd(&v[i + 40], vec);
			_mm512_storeu_pd(&v[i + 48], vec);
			_mm512_storeu_pd(&v[i + 56], vec);
		   }
#pragma loop_count min(1),avg(4),max(7)
	           for (; i != vlen; ++i){
		        v[i] = val;
		   }
#elif defined __GNUC__ && !defined(__INTEL_COMPILER)

	          if ((reinterpret_cast<uintptr_t>(v)& 0x3F) != 0ULL) {
		       for (i = 0; i != ROUND_TO_EIGHT(vlen, 8); i += 64) {
		   	     _mm512_storeu_pd(&v[i + 0], vec);
			     _mm512_storeu_pd(&v[i + 8], vec);
			     _mm512_storeu_pd(&v[i + 16], vec);
			     _mm512_storeu_pd(&v[i + 24], vec);
			     _mm512_storeu_pd(&v[i + 32], vec);
			     _mm512_storeu_pd(&v[i + 40], vec);
			     _mm512_storeu_pd(&v[i + 48], vec);
			     _mm512_storeu_pd(&v[i + 56], vec);
		       }
	            }
 	             else {
	                 for (i = 0; i != ROUND_TO_EIGHT(vlen, 8); i += 64) {
		            _mm512_store_pd(&v[i+0],vec);
		            _mm512_store_pd(&v[i+8],vec);
		            _mm512_store_pd(&v[i+16],vec);
		            _mm512_store_pd(&v[i+24],vec);
		            _mm512_store_pd(&v[i+32],vec);
		            _mm512_store_pd(&v[i+40],vec);
		            _mm512_store_pd(&v[i+48],vec);
		            _mm512_store_pd(&v[i+56],vec);
	                  }
	             }
	                 for (; i != vlen; ++i){
		              v[i] = val;
		         }
#endif
                  }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_unroll8x_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float  __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val) {
                    __m512 zmm0 = _mm512_set1_ps(val);
                    int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
                    for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 128) {
	                  _mm512_storeu_ps(&v[i+0],  zmm0);
	                  _mm512_storeu_ps(&v[i+16], zmm0);
	                  _mm512_storeu_ps(&v[i+32], zmm0);
	                  _mm512_storeu_ps(&v[i+48], zmm0);
	                  _mm512_storeu_ps(&v[i+64], zmm0);
	                  _mm512_storeu_ps(&v[i+80], zmm0);
	                  _mm512_storeu_ps(&v[i+96], zmm0);
	                  _mm512_storeu_ps(&v[i+112], zmm0);
	            }
#pragma loop_count min(1),avg(8),max(15)
                    for(; i != vlen; ++i) {
	                v[i] = val;
                    }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                     if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	                 for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 128) {
	                       _mm512_storeu_ps(&v[i+0],  zmm0);
	                       _mm512_storeu_ps(&v[i+16], zmm0);
	                       _mm512_storeu_ps(&v[i+32], zmm0);
	                       _mm512_storeu_ps(&v[i+48], zmm0);
	                       _mm512_storeu_ps(&v[i+64], zmm0);
	                       _mm512_storeu_ps(&v[i+80], zmm0);
	                       _mm512_storeu_ps(&v[i+96], zmm0);
	                       _mm512_storeu_ps(&v[i+112], zmm0);
	                 }
                     } else {
                         v = (float*)__builtin_assume_aligned(v,64);
                         for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 128) {
	                      _mm512_store_ps(&v[i+0LL],  zmm0);
	                      _mm512_store_ps(&v[i+16LL], zmm0);
	                      _mm512_store_ps(&v[i+32LL], zmm0);
	                      _mm512_store_ps(&v[i+48LL], zmm0);
	                      _mm512_store_ps(&v[i+64LL], zmm0);
	                      _mm512_store_ps(&v[i+80LL], zmm0);
	                      _mm512_store_ps(&v[i+96LL], zmm0);
	                      _mm512_store_ps(&v[i+112LL], zmm0);
	                 }
                      }
                         for(; i != vlen; ++i) {
	                     v[i] = val;
                            }
#endif
               }

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_init_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIGN__(64) v,
#endif
			              const int32_t vlen,
			              const float val) {
                    __m512 zmm = _mm512_set1_ps(val);
#if defined __ICC || defined __INTEL_COMPILER
                    if(vlen <= MEMMOVE_1ELEM) {
                       return;
		    }
		    else if(vlen <= MEMMOVE_16ELEMS) {
                        _mm512_storeu_ps(&v[0],zmm);
			return;
		    }
		    else if(vlen <= MEMMOVE_32ELEMS) {
                        _mm512_storeu_ps(&v[0],zmm);
			_mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			return;
		    }
		    else if(vlen <= MEMMOVE_64ELEMS) {
                        _mm512_storeu_ps(&v[0],zmm);
			_mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[2*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[3*ZMM_LEN],zmm);
			return;
		    }
		    else if(vlen <= MEMMOVE_128ELEMS) {
                        _mm512_storeu_ps(&v[0],zmm);
			_mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[2*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[3*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[4*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[5*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[6*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[7*ZMM_LEN],zmm);
			return;
		    }
		    else if(vlen <= MEMMOVE_256ELEMS) {
                        _mm512_storeu_ps(&v[0],zmm);
			_mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[2*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[3*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[4*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[5*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[6*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[7*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[8*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[9*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[10*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[11*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[12*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[13*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[14*ZMM_LEN],zmm);
			_mm512_storeu_ps(&v[15*ZMM_LEN],zmm);
			return;
		    }
		    else if(vlen > MEMMOVE_256ELEMS) {
                         int32_t i;
#pragma code_align(32)
                         for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 128) {
	                      _mm512_storeu_ps(&v[i+0],  zmm);
	                      _mm512_storeu_ps(&v[i+16], zmm);
	                      _mm512_storeu_ps(&v[i+32], zmm);
	                      _mm512_storeu_ps(&v[i+48], zmm);
	                      _mm512_storeu_ps(&v[i+64], zmm);
	                      _mm512_storeu_ps(&v[i+80], zmm);
	                      _mm512_storeu_ps(&v[i+96], zmm);
	                      _mm512_storeu_ps(&v[i+112], zmm);
	                 }
#pragma loop_count min(1),avg(8),max(15)
                        for(; i != vlen; ++i) {
	                    v[i] = val;
                        }
			return;
		    }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                     if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {

		           if(vlen <= MEMMOVE_1ELEM) {
                              return;
		           }
		           else if(vlen <= MEMMOVE_16ELEMS) {
                                 _mm512_storeu_ps(&v[0],zmm);
			         return;
		           }
		           else if(vlen <= MEMMOVE_32ELEMS) {
                                _mm512_storeu_ps(&v[0],zmm);
		 	        _mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			         return;
		           }
		           else if(vlen <= MEMMOVE_64ELEMS) {
                                _mm512_storeu_ps(&v[0],zmm);
		 	        _mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[2*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[3*ZMM_LEN],zmm);
			        return;
		           }
		           else if(vlen <= MEMMOVE_128ELEMS) {
                                _mm512_storeu_ps(&v[0],zmm);
		 	        _mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[2*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[3*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[4*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[5*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[6*ZMM_LEN],zmm);
			        _mm512_storeu_ps(&v[7*ZMM_LEN],zmm);
			        return;
		           }
			   else if(vlen <= MEMMOVE_256ELEMS) {
                                 _mm512_storeu_ps(&v[0],zmm);
		         	 _mm512_storeu_ps(&v[1*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[2*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[3*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[4*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[5*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[6*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[7*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[8*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[9*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[10*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[11*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[12*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[13*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[14*ZMM_LEN],zmm);
			         _mm512_storeu_ps(&v[15*ZMM_LEN],zmm);
			         return;
		            }
			     else if(vlen > MEMMOVE_256ELEMS) {
                                  for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 128) {
	                              _mm512_storeu_ps(&v[i+0],  zmm);
	                              _mm512_storeu_ps(&v[i+16], zmm);
	                              _mm512_storeu_ps(&v[i+32], zmm);
	                              _mm512_storeu_ps(&v[i+48], zmm);
	                              _mm512_storeu_ps(&v[i+64], zmm);
	                              _mm512_storeu_ps(&v[i+80], zmm);
	                              _mm512_storeu_ps(&v[i+96], zmm);
	                              _mm512_storeu_ps(&v[i+112], zmm);
	                          }
                                  for(; i != vlen; ++i) {
	                              v[i] = val;
                                  }
			          return;  
		    }
		     else {

		          if(vlen <= MEMMOVE_1ELEM) {
                              return;
		           }
		           else if(vlen <= MEMMOVE_16ELEMS) {
                                 _mm512_store_ps(&v[0],zmm);
			         return;
		           }
		           else if(vlen <= MEMMOVE_32ELEMS) {
                                _mm512_store_ps(&v[0],zmm);
		 	        _mm512_store_ps(&v[1*ZMM_LEN],zmm);
			         return;
		           }
		           else if(vlen <= MEMMOVE_64ELEMS) {
                                _mm512_store_ps(&v[0],zmm);
		 	        _mm512_store_ps(&v[1*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[2*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[3*ZMM_LEN],zmm);
			        return;
		           }
		           else if(vlen <= MEMMOVE_128ELEMS) {
                                _mm512_store_ps(&v[0],zmm);
		 	        _mm512_store_ps(&v[1*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[2*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[3*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[4*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[5*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[6*ZMM_LEN],zmm);
			        _mm512_store_ps(&v[7*ZMM_LEN],zmm);
			        return;
		           }
			   else if(vlen <= MEMMOVE_256ELEMS) {
                                 _mm512_store_ps(&v[0],zmm);
		         	 _mm512_store_ps(&v[1*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[2*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[3*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[4*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[5*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[6*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[7*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[8*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[9*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[10*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[11*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[12*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[13*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[14*ZMM_LEN],zmm);
			         _mm512_store_ps(&v[15*ZMM_LEN],zmm);
			         return;
		            }
			     else if(vlen > MEMMOVE_256ELEMS) {
			          v = (float*)__builtin_assume_aligned(v,64);
                                  for(i = 0; i != ROUND_TO_SIXTEEN(vlen,16); i += 128) {
	                              _mm512_store_ps(&v[i+0],  zmm);
	                              _mm512_store_ps(&v[i+16], zmm);
	                              _mm512_store_ps(&v[i+32], zmm);
	                              _mm512_store_ps(&v[i+48], zmm);
	                              _mm512_store_ps(&v[i+64], zmm);
	                              _mm512_store_ps(&v[i+80], zmm);
	                              _mm512_store_ps(&v[i+96], zmm);
	                              _mm512_store_ps(&v[i+112], zmm);
	                          }
                                  for(; i != vlen; ++i) {
	                              v[i] = val;
                                  }
			          return;  
		         }
#endif
                    }

#endif // __AVX512F__

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_init_ps(
#if defined __ICC || defined __INTEL_COMPILER
	                              float * __restrict v,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                      float * __restrict __ATTR_ALIGN__(32) v,
#endif
			              const int32_t vlen,
			              const float val) {
                 __m256 ymm = _mm256_set1_ps(val);
#if defined __ICC || defined __INTEL_COMPILER
                 if(vlen <= MEMMOVE_1ELEM) {
                    return;
		 }
		 else if(vlen <= MEMMOVE_16ELEMS) {
                      _mm256_storeu_ps(&v[0],ymm);
		      _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		      return;
		 }
		 else if(vlen <= MEMMOVE_32ELEMS) {
                      _mm256_storeu_ps(&v[0],ymm);
		      _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[2*YMM_LEN],ymm)
		      _mm256_storeu_ps(&v[3*YMM_LEN],ymm);
		      return;
		 }
		 else if(vlen <= MEMMOVE_64ELEMS) {
                      _mm256_storeu_ps(&v[0],ymm);
		      _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[2*YMM_LEN],ymm)
		      _mm256_storeu_ps(&v[3*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[4*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[5*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[6*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[7*YMM_LEN],ymm);
		      return;
		 }
		 else if(vlen <= MEMMOVE_128ELEMS) {
                      _mm256_storeu_ps(&v[0],ymm);
		      _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[2*YMM_LEN],ymm)
		      _mm256_storeu_ps(&v[3*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[4*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[5*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[6*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[7*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[8*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[9*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[10*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[11*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[12*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[13*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[14*YMM_LEN],ymm);
		      _mm256_storeu_ps(&v[15*YMM_LEN],ymm);
		      return; 
		 }
		 else if(vlen > MEMMOVE_128ELEMS) {
                      int32_t i;
#pragma code_align(32)
                      for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 64) {
	                    _mm256_storeu_ps(&v[i+0],  ymm0);
	                    _mm256_storeu_ps(&v[i+8],  ymm0);
	                    _mm256_storeu_ps(&v[i+16], ymm0);
	                    _mm256_storeu_ps(&v[i+24], ymm0);
	                    _mm256_storeu_ps(&v[i+32], ymm0);
	                    _mm256_storeu_ps(&v[i+40], ymm0);
	                    _mm256_storeu_ps(&v[i+48], ymm0);
	                    _mm256_storeu_ps(&v[i+56], ymm0);
	              }
#pragma loop_count min(1),avg(4),max(7)
                        for(; i != vlen; ++i) {
	                    v[i] = val;
                           }
		 }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                  if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {

		      if(vlen <= MEMMOVE_1ELEM) {
                         return;
		      }
		      else if(vlen <= MEMMOVE_16ELEMS) {
                           _mm256_storeu_ps(&v[0],ymm);
		           _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		         return;
		      }
		      else if(vlen <= MEMMOVE_32ELEMS) {
                           _mm256_storeu_ps(&v[0],ymm);
		           _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		           _mm256_storeu_ps(&v[2*YMM_LEN],ymm)
		           _mm256_storeu_ps(&v[3*YMM_LEN],ymm);
		           return;
		      }
		      else if(vlen <= MEMMOVE_64ELEMS) {
                            _mm256_storeu_ps(&v[0],ymm);
		            _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[2*YMM_LEN],ymm)
		            _mm256_storeu_ps(&v[3*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[4*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[5*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[6*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[7*YMM_LEN],ymm);
		            return;
		      }
		      else if(vlen <= MEMMOVE_128ELEMS) {
                            _mm256_storeu_ps(&v[0],ymm);
		            _mm256_storeu_ps(&v[1*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[2*YMM_LEN],ymm)
		            _mm256_storeu_ps(&v[3*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[4*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[5*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[6*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[7*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[8*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[9*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[10*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[11*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[12*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[13*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[14*YMM_LEN],ymm);
		            _mm256_storeu_ps(&v[15*YMM_LEN],ymm);
		            return; 
		 }
		 else if(vlen > MEMMOVE_128ELEMS) {
                      int32_t i;
                        for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 64) {
	                    _mm256_storeu_ps(&v[i+0],  ymm0);
	                    _mm256_storeu_ps(&v[i+8],  ymm0);
	                    _mm256_storeu_ps(&v[i+16], ymm0);
	                    _mm256_storeu_ps(&v[i+24], ymm0);
	                    _mm256_storeu_ps(&v[i+32], ymm0);
	                    _mm256_storeu_ps(&v[i+40], ymm0);
	                    _mm256_storeu_ps(&v[i+48], ymm0);
	                    _mm256_storeu_ps(&v[i+56], ymm0);
	                }
                        for(; i != vlen; ++i) {
	                    v[i] = val;
                           }
		  }else {

		         if(vlen <= MEMMOVE_1ELEM) {
                               return;
		         }
		         else if(vlen <= MEMMOVE_16ELEMS) {
                                 _mm256_store_ps(&v[0],ymm);
		                 _mm256_store_ps(&v[1*YMM_LEN],ymm);
		                 return;
		        }
		        else if(vlen <= MEMMOVE_32ELEMS) {
                                _mm256_store_ps(&v[0],ymm);
		                _mm256_store_ps(&v[1*YMM_LEN],ymm);
		                _mm256_store_ps(&v[2*YMM_LEN],ymm)
		                _mm256_store_ps(&v[3*YMM_LEN],ymm);
		                return;
		         }
		         else if(vlen <= MEMMOVE_64ELEMS) {
                                _mm256_store_ps(&v[0],ymm);
		                _mm256_store_ps(&v[1*YMM_LEN],ymm);
		                _mm256_store_ps(&v[2*YMM_LEN],ymm)
		                _mm256_store_ps(&v[3*YMM_LEN],ymm);
		                _mm256_store_ps(&v[4*YMM_LEN],ymm);
		                _mm256_store_ps(&v[5*YMM_LEN],ymm);
		                _mm256_store_ps(&v[6*YMM_LEN],ymm);
		                _mm256_store_ps(&v[7*YMM_LEN],ymm);
		                return;
		        }
		        else if(vlen <= MEMMOVE_128ELEMS) {
                                _mm256_store_ps(&v[0],ymm);
		                _mm256_store_ps(&v[1*YMM_LEN],ymm);
		                _mm256_store_ps(&v[2*YMM_LEN],ymm)
		                _mm256_store_ps(&v[3*YMM_LEN],ymm);
		                _mm256_store_ps(&v[4*YMM_LEN],ymm);
		                _mm256_store_ps(&v[5*YMM_LEN],ymm);
		                _mm256_store_ps(&v[6*YMM_LEN],ymm);
		                _mm256_store_ps(&v[7*YMM_LEN],ymm);
		                _mm256_store_ps(&v[8*YMM_LEN],ymm);
		                _mm256_store_ps(&v[9*YMM_LEN],ymm);
		                _mm256_store_ps(&v[10*YMM_LEN],ymm);
		                _mm256_store_ps(&v[11*YMM_LEN],ymm);
		                _mm256_store_ps(&v[12*YMM_LEN],ymm);
		                _mm256_store_ps(&v[13*YMM_LEN],ymm);
		                _mm256_store_ps(&v[14*YMM_LEN],ymm);
		                _mm256_store_ps(&v[15*YMM_LEN],ymm);
		                return; 
		       }
		       else if(vlen > MEMMOVE_128ELEMS) {
		              v = (float*)__builtin_assume_aligned(v,32);
                              int32_t i;
                              for(i = 0; i != ROUND_TO_EIGHT(vlen,8); i += 64) {
	                          _mm256_store_ps(&v[i+0],  ymm0);
	                          _mm256_store_ps(&v[i+8],  ymm0);
	                          _mm256_store_ps(&v[i+16], ymm0);
	                          _mm256_store_ps(&v[i+24], ymm0);
	                          _mm256_store_ps(&v[i+32], ymm0);
	                          _mm256_store_ps(&v[i+40], ymm0);
	                          _mm256_store_ps(&v[i+48], ymm0);
	                          _mm256_store_ps(&v[i+56], ymm0);
	                      }
                              for(; i != vlen; ++i) {
	                          v[i] = val;
                              }
		       }
#endif
                  }

#include <complex.h>


         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_copy_2r4_c4_u(const float * __restrict re,
                                   const float * __restrict im,
                                   std::complex<float> * __restrict c,
                                   const int32_t len) {

               float * __restrict fptr = reinterpret_cast<float*>(&c[0]);
               const int32_t len2      = 2*len;
               if(len <= MEMMOVE_1ELEM) { 
                  return;
               }
               else
        }


	 
         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void  avx256_memcpy_ps(
#if defined __ICC || defined __INTEL_COMPILER
                                   float * __restrict dst,
				   const float * __restrict src,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                   float * __restrict __ATTR_ALIGN__(32) dst,
				   const float * __restrict __ATTR_ALIGN__(32) src,
#endif
                                   const int32_t len) {
#if defined __ICC || defined __INTEL_COMPILER
                    if(len <= MEMMOVE_1ELEM) {
                       return;
		    }
		    else if(len <= MEMMOVE_16ELEMS) {
                            const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			    _mm256_storeu_ps(&dst[0],ymm0);
			    const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			    _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			    return;
		    }
		    else if(len <= MEMMOVE_32ELEMS) {
                            const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			    _mm256_storeu_ps(&dst[0],ymm0);
			    const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			    _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			    const __m256 ymm2 = _mm256_loadu_ps(&src[2*YMM_LEN]);
			    _mm256_storeu_ps(&dst[2*YMM_LEN],ymm2);
			    const _mm256 ymm3 = _mm256_loadu_ps(&src[3*YMM_LEN]);
			    _mm256_storeu_ps(&dst[3*YMM_LEN],ymm3);
			    return;
		    }
		    else if(len <= MEMMOVE_64ELEMS) {
                            const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			    _mm256_storeu_ps(&dst[0],ymm0);
			    const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			    _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			    const __m256 ymm2 = _mm256_loadu_ps(&src[2*YMM_LEN]);
			    _mm256_storeu_ps(&dst[2*YMM_LEN],ymm2);
			    const _mm256 ymm3 = _mm256_loadu_ps(&src[3*YMM_LEN]);
			    _mm256_storeu_ps(&dst[3*YMM_LEN],ymm3);
			    const __m256 ymm4 = _mm256_loadu_ps(&src[4*YMM_LEN]);
			    _mm256_storeu_ps(&dst[4*YMM_LEN],ymm4);
			    const __m256 ymm5 = _mm256_loadu_ps(&src[5*YMM_LEN]);
			    _mm256_storeu_ps(&dst[5*YMM_LEN],ymm5);
			    const __m256 ymm6 = _mm256_loadu_ps(&src[6*YMM_LEN]);
			    _mm256_storeu_ps(&dst[6*YMM_LEN],ymm6);
			    const __m256 ymm7 = _mm256_loadu_ps(&src[7*YMM_LEN]);
			    _mm256_storeu_ps(&dst[7*YMM_LEN],ymm7);
			    return;
		    }
		    else if(len <= MEMMOVE_128ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
		            _mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		            _mm_prefetch((const char*)&src[2*YMM_LEN],_MM_HINT_T0);
		            _mm_prefetch((const char*)&src[4*YMM_LEN],_MM_HINT_T0);
		            _mm_prefetch((const char*)&src[6*YMM_LEN],_MM_HINT_T0);
		            _mm_prefetch((const char*)&src[8*YMM_LEN],_MM_HINT_T0);
#endif
                            const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			    _mm256_storeu_ps(&dst[0],ymm0);
			    const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			    _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			    const __m256 ymm2 = _mm256_loadu_ps(&src[2*YMM_LEN]);
			    _mm256_storeu_ps(&dst[2*YMM_LEN],ymm2);
			    const _mm256 ymm3 = _mm256_loadu_ps(&src[3*YMM_LEN]);
			    _mm256_storeu_ps(&dst[3*YMM_LEN],ymm3);
			    const __m256 ymm4 = _mm256_loadu_ps(&src[4*YMM_LEN]);
			    _mm256_storeu_ps(&dst[4*YMM_LEN],ymm4);
			    const __m256 ymm5 = _mm256_loadu_ps(&src[5*YMM_LEN]);
			    _mm256_storeu_ps(&dst[5*YMM_LEN],ymm5);
			    const __m256 ymm6 = _mm256_loadu_ps(&src[6*YMM_LEN]);
			    _mm256_storeu_ps(&dst[6*YMM_LEN],ymm6);
			    const __m256 ymm7 = _mm256_loadu_ps(&src[7*YMM_LEN]);
			    _mm256_storeu_ps(&dst[7*YMM_LEN],ymm7);
			    const __m256 ymm8 = _mm256_loadu_ps(&src[8*YMM_LEN]);
			    _mm256_storeu_ps(&dst[8*YMM_LEN],ymm8);
			    const __m256 ymm9 = _mm256_loadu_ps(&src[9*YMM_LEN]);
			    _mm256_storeu_ps(&dst[9*YMM_LEN],ymm9);
			    const __m256 ymm10 = _mm256_loadu_ps(&src[10*YMM_LEN]);
			    _mm256_storeu_ps(&dst[10*YMM_LEN],ymm10);
			    const __m256 ymm11 = _mm256_loadu_ps(&src[11*YMM_LEN]);
			    _mm256_storeu_ps(&dst[11*YMM_LEN],ymm11);
			    const __m256 ymm12 = _mm256_loadu_ps(&src[12*YMM_LEN]);
			    _mm256_storeu_ps(&dst[12*YMM_LEN],ymm12);
			    const __m256 ymm13 = _mm256_loadu_ps(&src[13*YMM_LEN]);
			    _mm256_storeu_ps(&dst[13*YMM_LEN],ymm13);
			    const __m256 ymm14 = _mm256_loadu_ps(&src[14*YMM_LEN]);
			    _mm256_storeu_ps(&dst[14*YMM_LEN],ymm14);
			    const __m256 ymm15 = _mm256_loadu_ps(&src[15*YMM_LEN]);
			    _mm256_storeu_ps(&dst[15*YMM_LEN],ymm155);
			    return;
		    }
		    else if(len > MEMMOVE_128ELEMS) {
		            int32_t i;
#pragma code_align(32)
                            for(i = 0; i != ROUND_TO_EIGHT(len,8); i += 64) {
#if (GMS_MAN_PREFETCH) == 1
                            _mm_prefetch((const char*)&src[i+0], _MM_HINT_T0);
		            _mm_prefetch((const char*)&src[i+2*YMM_LEN],_MM_HINT_T0);
		            _mm_prefetch((const char*)&src[i+4*YMM_LEN],_MM_HINT_T0);
#endif
                            const __m256 ymm0 = _mm256_loadu_ps(&src[i+0]);
			    _mm256_storeu_ps(&dst[i+0],ymm0);
			    const __m256 ymm1 = _mm256_loadu_ps(&src[i+1*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+1*YMM_LEN],ymm1);
			    const __m256 ymm2 = _mm256_loadu_ps(&src[i+2*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+2*YMM_LEN],ymm2);
			    const _mm256 ymm3 = _mm256_loadu_ps(&src[i+3*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+3*YMM_LEN],ymm3);
			    const __m256 ymm4 = _mm256_loadu_ps(&src[i+4*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+4*YMM_LEN],ymm4);
			    const __m256 ymm5 = _mm256_loadu_ps(&src[i+5*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+5*YMM_LEN],ymm5);
			    const __m256 ymm6 = _mm256_loadu_ps(&src[i+6*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+6*YMM_LEN],ymm6);
			    const __m256 ymm7 = _mm256_loadu_ps(&src[i+7*YMM_LEN]);
			    _mm256_storeu_ps(&dst[i+7*YMM_LEN],ymm7);  
			    }
#pragma loop_count min(1),avg(4),max(7)
                           for(; i != len; ++i) {
                               dst[i] = src[i];
			   }
			   return;
		    }
		    
 
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                   	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
	                    (reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {

			        if(len <= MEMMOVE_1ELEM) {
                                   return;
		                }
		                else if(len <= MEMMOVE_16ELEMS) {
                                   const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			           _mm256_storeu_ps(&dst[0],ymm0);
			           const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			           _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			           return;
		                }
		                else if(len <= MEMMOVE_32ELEMS) {
                                  const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			          _mm256_storeu_ps(&dst[0],ymm0);
			          const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			          _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			          const __m256 ymm2 = _mm256_loadu_ps(&src[2*YMM_LEN]);
			          _mm256_storeu_ps(&dst[2*YMM_LEN],ymm2);
			          const _mm256 ymm3 = _mm256_loadu_ps(&src[3*YMM_LEN]);
			          _mm256_storeu_ps(&dst[3*YMM_LEN],ymm3);
			          return;
		                }
		                else if(len <= MEMMOVE_64ELEMS) {
                                    const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			            _mm256_storeu_ps(&dst[0],ymm0);
			            const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			            _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			            const __m256 ymm2 = _mm256_loadu_ps(&src[2*YMM_LEN]);
			            _mm256_storeu_ps(&dst[2*YMM_LEN],ymm2);
			            const _mm256 ymm3 = _mm256_loadu_ps(&src[3*YMM_LEN]);
			            _mm256_storeu_ps(&dst[3*YMM_LEN],ymm3);
			            const __m256 ymm4 = _mm256_loadu_ps(&src[4*YMM_LEN]);
			            _mm256_storeu_ps(&dst[4*YMM_LEN],ymm4);
			            const __m256 ymm5 = _mm256_loadu_ps(&src[5*YMM_LEN]);
			            _mm256_storeu_ps(&dst[5*YMM_LEN],ymm5);
			            const __m256 ymm6 = _mm256_loadu_ps(&src[6*YMM_LEN]);
			            _mm256_storeu_ps(&dst[6*YMM_LEN],ymm6);
			            const __m256 ymm7 = _mm256_loadu_ps(&src[7*YMM_LEN]);
			            _mm256_storeu_ps(&dst[7*YMM_LEN],ymm7);
			            return;
		                }
		                else if(len <= MEMMOVE_128ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
		                    _mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[2*YMM_LEN],_MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[4*YMM_LEN],_MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[6*YMM_LEN],_MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[8*YMM_LEN],_MM_HINT_T0);
#endif
                                    const __m256 ymm0 = _mm256_loadu_ps(&src[0]);
			            _mm256_storeu_ps(&dst[0],ymm0);
			            const __m256 ymm1 = _mm256_loadu_ps(&src[1*YMM_LEN]);
			            _mm256_storeu_ps(&dst[1*YMM_LEN],ymm1);
			            const __m256 ymm2 = _mm256_loadu_ps(&src[2*YMM_LEN]);
			            _mm256_storeu_ps(&dst[2*YMM_LEN],ymm2);
			            const _mm256 ymm3 = _mm256_loadu_ps(&src[3*YMM_LEN]);
			            _mm256_storeu_ps(&dst[3*YMM_LEN],ymm3);
			            const __m256 ymm4 = _mm256_loadu_ps(&src[4*YMM_LEN]);
			            _mm256_storeu_ps(&dst[4*YMM_LEN],ymm4);
			            const __m256 ymm5 = _mm256_loadu_ps(&src[5*YMM_LEN]);
			            _mm256_storeu_ps(&dst[5*YMM_LEN],ymm5);
			            const __m256 ymm6 = _mm256_loadu_ps(&src[6*YMM_LEN]);
			            _mm256_storeu_ps(&dst[6*YMM_LEN],ymm6);
			            const __m256 ymm7 = _mm256_loadu_ps(&src[7*YMM_LEN]);
			            _mm256_storeu_ps(&dst[7*YMM_LEN],ymm7);
			            const __m256 ymm8 = _mm256_loadu_ps(&src[8*YMM_LEN]);
			            _mm256_storeu_ps(&dst[8*YMM_LEN],ymm8);
			            const __m256 ymm9 = _mm256_loadu_ps(&src[9*YMM_LEN]);
			            _mm256_storeu_ps(&dst[9*YMM_LEN],ymm9);
			            const __m256 ymm10 = _mm256_loadu_ps(&src[10*YMM_LEN]);
			            _mm256_storeu_ps(&dst[10*YMM_LEN],ymm10);
			            const __m256 ymm11 = _mm256_loadu_ps(&src[11*YMM_LEN]);
			            _mm256_storeu_ps(&dst[11*YMM_LEN],ymm11);
			            const __m256 ymm12 = _mm256_loadu_ps(&src[12*YMM_LEN]);
			            _mm256_storeu_ps(&dst[12*YMM_LEN],ymm12);
			            const __m256 ymm13 = _mm256_loadu_ps(&src[13*YMM_LEN]);
			            _mm256_storeu_ps(&dst[13*YMM_LEN],ymm13);
			             const __m256 ymm14 = _mm256_loadu_ps(&src[14*YMM_LEN]);
			            _mm256_storeu_ps(&dst[14*YMM_LEN],ymm14);
			            const __m256 ymm15 = _mm256_loadu_ps(&src[15*YMM_LEN]);
			            _mm256_storeu_ps(&dst[15*YMM_LEN],ymm155);
			             return;
		             }
		             else if(len > MEMMOVE_128ELEMS) {
		                     int32_t i;

                                      for(i = 0; i != ROUND_TO_EIGHT(len,8); i += 64) {
#if (GMS_MAN_PREFETCH) == 1
                                        _mm_prefetch((const char*)&src[i+0], _MM_HINT_T0);
		                        _mm_prefetch((const char*)&src[i+2*YMM_LEN],_MM_HINT_T0);
		                        _mm_prefetch((const char*)&src[i+4*YMM_LEN],_MM_HINT_T0);
#endif
                                        const __m256 ymm0 = _mm256_loadu_ps(&src[i+0]);
			                _mm256_storeu_ps(&dst[i+0],ymm0);
			                const __m256 ymm1 = _mm256_loadu_ps(&src[i+1*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+1*YMM_LEN],ymm1);
			                const __m256 ymm2 = _mm256_loadu_ps(&src[i+2*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+2*YMM_LEN],ymm2);
			                const _mm256 ymm3 = _mm256_loadu_ps(&src[i+3*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+3*YMM_LEN],ymm3);
			                const __m256 ymm4 = _mm256_loadu_ps(&src[i+4*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+4*YMM_LEN],ymm4);
			                const __m256 ymm5 = _mm256_loadu_ps(&src[i+5*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+5*YMM_LEN],ymm5);
			                const __m256 ymm6 = _mm256_loadu_ps(&src[i+6*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+6*YMM_LEN],ymm6);
			                const __m256 ymm7 = _mm256_loadu_ps(&src[i+7*YMM_LEN]);
			                _mm256_storeu_ps(&dst[i+7*YMM_LEN],ymm7);  
			               }
 
                                       for(; i != len; ++i) {
                                           dst[i] = src[i];
			                }
			                return;
		                    }
			      }
			       else {

			        if(len <= MEMMOVE_1ELEM) {
                                         return;
		                }
		                else if(len <= MEMMOVE_16ELEMS) {
                                   const __m256 ymm0 = _mm256_load_ps(&src[0]);
			           _mm256_store_ps(&dst[0],ymm0);
			           const __m256 ymm1 = _mm256_load_ps(&src[1*YMM_LEN]);
			           _mm256_store_ps(&dst[1*YMM_LEN],ymm1);
			           return;
		                }
		                else if(len <= MEMMOVE_32ELEMS) {
                                  const __m256 ymm0 = _mm256_load_ps(&src[0]);
			          _mm256_store_ps(&dst[0],ymm0);
			          const __m256 ymm1 = _mm256_load_ps(&src[1*YMM_LEN]);
			          _mm256_store_ps(&dst[1*YMM_LEN],ymm1);
			          const __m256 ymm2 = _mm256_load_ps(&src[2*YMM_LEN]);
			          _mm256_store_ps(&dst[2*YMM_LEN],ymm2);
			          const _mm256 ymm3 = _mm256_load_ps(&src[3*YMM_LEN]);
			          _mm256_store_ps(&dst[3*YMM_LEN],ymm3);
			          return;
		                }
		                else if(len <= MEMMOVE_64ELEMS) {
                                    const __m256 ymm0 = _mm256_load_ps(&src[0]);
			            _mm256_store_ps(&dst[0],ymm0);
			            const __m256 ymm1 = _mm256_load_ps(&src[1*YMM_LEN]);
			            _mm256_store_ps(&dst[1*YMM_LEN],ymm1);
			            const __m256 ymm2 = _mm256_load_ps(&src[2*YMM_LEN]);
			            _mm256_store_ps(&dst[2*YMM_LEN],ymm2);
			            const _mm256 ymm3 = _mm256_load_ps(&src[3*YMM_LEN]);
			            _mm256_store_ps(&dst[3*YMM_LEN],ymm3);
			            const __m256 ymm4 = _mm256_load_ps(&src[4*YMM_LEN]);
			            _mm256_store_ps(&dst[4*YMM_LEN],ymm4);
			            const __m256 ymm5 = _mm256_load_ps(&src[5*YMM_LEN]);
			            _mm256_store_ps(&dst[5*YMM_LEN],ymm5);
			            const __m256 ymm6 = _mm256_load_ps(&src[6*YMM_LEN]);
			            _mm256_store_ps(&dst[6*YMM_LEN],ymm6);
			            const __m256 ymm7 = _mm256_load_ps(&src[7*YMM_LEN]);
			            _mm256_store_ps(&dst[7*YMM_LEN],ymm7);
			            return;
		                }
		                else if(len <= MEMMOVE_128ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
		                    _mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[2*YMM_LEN],_MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[4*YMM_LEN],_MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[6*YMM_LEN],_MM_HINT_T0);
		                    _mm_prefetch((const char*)&src[8*YMM_LEN],_MM_HINT_T0);
#endif
                                    const __m256 ymm0 = _mm256_load_ps(&src[0]);
			            _mm256_store_ps(&dst[0],ymm0);
			            const __m256 ymm1 = _mm256_load_ps(&src[1*YMM_LEN]);
			            _mm256_store_ps(&dst[1*YMM_LEN],ymm1);
			            const __m256 ymm2 = _mm256_load_ps(&src[2*YMM_LEN]);
			            _mm256_store_ps(&dst[2*YMM_LEN],ymm2);
			            const _mm256 ymm3 = _mm256_load_ps(&src[3*YMM_LEN]);
			            _mm256_store_ps(&dst[3*YMM_LEN],ymm3);
			            const __m256 ymm4 = _mm256_load_ps(&src[4*YMM_LEN]);
			            _mm256_store_ps(&dst[4*YMM_LEN],ymm4);
			            const __m256 ymm5 = _mm256_load_ps(&src[5*YMM_LEN]);
			            _mm256_store_ps(&dst[5*YMM_LEN],ymm5);
			            const __m256 ymm6 = _mm256_load_ps(&src[6*YMM_LEN]);
			            _mm256_store_ps(&dst[6*YMM_LEN],ymm6);
			            const __m256 ymm7 = _mm256_load_ps(&src[7*YMM_LEN]);
			            _mm256_store_ps(&dst[7*YMM_LEN],ymm7);
			            const __m256 ymm8 = _mm256_load_ps(&src[8*YMM_LEN]);
			            _mm256_store_ps(&dst[8*YMM_LEN],ymm8);
			            const __m256 ymm9 = _mm256_load_ps(&src[9*YMM_LEN]);
			            _mm256_store_ps(&dst[9*YMM_LEN],ymm9);
			            const __m256 ymm10 = _mm256_load_ps(&src[10*YMM_LEN]);
			            _mm256_store_ps(&dst[10*YMM_LEN],ymm10);
			            const __m256 ymm11 = _mm256_load_ps(&src[11*YMM_LEN]);
			            _mm256_store_ps(&dst[11*YMM_LEN],ymm11);
			            const __m256 ymm12 = _mm256_load_ps(&src[12*YMM_LEN]);
			            _mm256_store_ps(&dst[12*YMM_LEN],ymm12);
			            const __m256 ymm13 = _mm256_load_ps(&src[13*YMM_LEN]);
			            _mm256_store_ps(&dst[13*YMM_LEN],ymm13);
			             const __m256 ymm14 = _mm256_load_ps(&src[14*YMM_LEN]);
			            _mm256_store_ps(&dst[14*YMM_LEN],ymm14);
			            const __m256 ymm15 = _mm256_load_ps(&src[15*YMM_LEN]);
			            _mm256_store_ps(&dst[15*YMM_LEN],ymm155);
			             return;
		             }
		             else if(len > MEMMOVE_128ELEMS) {
		                      int32_t i;
                                      src = (float*)__builtin_assume_aligned(src,32);
				      dst = (const float*)__builtin_assume_aligned(dst,32);
                                      for(i = 0; i != ROUND_TO_EIGHT(len,8); i += 64) {
#if (GMS_MAN_PREFETCH) == 1
                                        _mm_prefetch((const char*)&src[i+0], _MM_HINT_T0);
		                        _mm_prefetch((const char*)&src[i+2*YMM_LEN],_MM_HINT_T0);
		                        _mm_prefetch((const char*)&src[i+4*YMM_LEN],_MM_HINT_T0);
#endif
                                        const __m256 ymm0 = _mm256_load_ps(&src[i+0]);
			                _mm256_store_ps(&dst[i+0],ymm0);
			                const __m256 ymm1 = _mm256_load_ps(&src[i+1*YMM_LEN]);
			                _mm256_store_ps(&dst[i+1*YMM_LEN],ymm1);
			                const __m256 ymm2 = _mm256_load_ps(&src[i+2*YMM_LEN]);
			                _mm256_store_ps(&dst[i+2*YMM_LEN],ymm2);
			                const _mm256 ymm3 = _mm256_load_ps(&src[i+3*YMM_LEN]);
			                _mm256_store_ps(&dst[i+3*YMM_LEN],ymm3);
			                const __m256 ymm4 = _mm256_load_ps(&src[i+4*YMM_LEN]);
			                _mm256_store_ps(&dst[i+4*YMM_LEN],ymm4);
			                const __m256 ymm5 = _mm256_load_ps(&src[i+5*YMM_LEN]);
			                _mm256_store_ps(&dst[i+5*YMM_LEN],ymm5);
			                const __m256 ymm6 = _mm256_load_ps(&src[i+6*YMM_LEN]);
			                _mm256_store_ps(&dst[i+6*YMM_LEN],ymm6);
			                const __m256 ymm7 = _mm256_load_ps(&src[i+7*YMM_LEN]);
			                _mm256_storeups(&dst[i+7*YMM_LEN],ymm7);  
			               }
 
                                         for(; i != len; ++i) {
                                             dst[i] = src[i];
			                }
			                return;
			     }
#endif
		  
                         }

#if defined __AVX512F__

          __ATTR_ALWAYS_INLINE__
	  __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void  avx512_memcpy_ps(
#if defined __ICC || defined __INTEL_COMPILER
                                   float * __restrict dst,
				   const float * __restrict src,
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                   float * __restrict __ATTR_ALIGN__(64) dst,
				   const float * __restrict __ATTR_ALIGN__(64) src,
#endif
                                   const int32_t len) {
                   
#endif
#if defined __ICC || defined __INTEL_COMPILER
                     if(len <= MEMMOVE_1ELEM) {
                        return;
		     }
		     else if(len <= MEMMOVE_16ELEMS) {
                        const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
			_mm512_storeu_ps(&dst[0],zmm0);
			return;
		     }
		     else if(len <= MEMMOVE_32ELEMS) {
                        const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
			_mm512_storeu_ps(&dst[0],zmm0);
			const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			_mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			return;
		     }
		     else if(len <= MEMMOVE_64ELEMS) {
                        const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
			_mm512_storeu_ps(&dst[0],zmm0);
			const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			_mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			const __m512 zmm2 = _mm512_loadu_ps(&src[2*ZMM_LEN]);
			_mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			const __m512 zmm3 = _mm512_loadu_ps(&src[3*ZMM_LEN]);
			_mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			return;
		     }
		     else if(len <= MEMMOVE_128ELEMS) {
                         const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
		 	 _mm512_storeu_ps(&dst[0],zmm0);
			 const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			 const __m512 zmm2 = _mm512_loadu_ps(&src[2*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			 const __m512 zmm3 = _mm512_loadu_ps(&src[3*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			 const __m512 zmm4 = _mm512_loadu_ps(&src[4*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[4*ZMM_LEN],zmm4);
			 const __m512 zmm5 = _mm512_loadu_ps(&src[5*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[5*ZMM_LEN],zmm5);
			 const __m512 zmm6 = _mm512_loadu_ps(&src[6*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[6*ZMM_LEN],zmm6);
			 const __m512 zmm7 = _m512_loadu_ps(&src[7*ZMM_LEN]);
			 _mm512_storeu_ps(&dst[7*ZMM_LEN],zmm7);
			return;
		     }
		     else if(len <= MEMMOVE_256ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
                             	 _mm_prefetch((const char *)&src[0],          _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[1*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[2*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[3*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[4*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[5*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[6*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[7*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[8*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[9*ZMM_LEN],  _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[10*ZMM_LEN], _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[11*ZMM_LEN], _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[12*ZMM_LEN], _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[13*ZMM_LEN], _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[14*ZMM_LEN], _MM_HINT_T0);
		                 _mm_prefetch((const char *)&src[15*ZMM_LEN], _MM_HINT_T0);
#endif
                                 const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
		 	         _mm512_storeu_ps(&dst[0],zmm0);
			         const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			         const __m512 zmm2 = _mm512_loadu_ps(&src[2*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			         const __m512 zmm3 = _mm512_loadu_ps(&src[3*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			         const __m512 zmm4 = _mm512_loadu_ps(&src[4*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[4*ZMM_LEN],zmm4);
			         const __m512 zmm5 = _mm512_loadu_ps(&src[5*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[5*ZMM_LEN],zmm5);
			         const __m512 zmm6 = _mm512_loadu_ps(&src[6*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[6*ZMM_LEN],zmm6);
			         const __m512 zmm7 = _mm512_loadu_ps(&src[7*ZMM_LEN]);
			         _mm512_storeu_ps(&dst[7*ZMM_LEN],zmm7);
				 const __m512 zmm8 = _mm512_loadu_ps(&src[8*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[8*ZMM_LEN],zmm8);
				 const __m512 zmm9 = _mm512_loadu_ps(&src[9*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[9*ZMM_LEN],zmm9);
				 const __m512 zmm10 = _mm512_loadu_ps(&src[10*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[10*ZMM_LEN],zmm10);
				 const __m512 zmm11 = _mm512_loadu_ps(&src[11*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[11*ZMM_LEN],zmm11);
				 const __m512 zmm12 = _mm512_loadu_ps(&src[12*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[12*ZMM_LEN],zmm12);
				 const __m512 zmm13 = _mm512_loadu_ps(&src[13*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[13*ZMM_LEN],zmm13);
				 const __m512 zmm14 = _mm512_loadu_ps(&src[14*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[14*ZMM_LEN],zmm14);
				 const __m512 zmm15 = _mm512_loadu_ps(&src[15*ZMM_LEN]);
				 _mm512_storeu_ps(&dst[15*ZMM_LEN],zmm15);
			         return;
		     }
		     else if(len > MEMMOVE_256ELEMS) {
                                 int32_t i;
#pragma code_align(32)
				 for(i = 0; i != ROUND_TO_SIXTEEN(len,16); i += 128) {
#if (GMS_MAN_PREFETCH) == 1
                                           _mm_prefetch((const char *)&src[i+0],          _MM_HINT_T0);
		                           _mm_prefetch((const char *)&src[i+1*ZMM_LEN],  _MM_HINT_T0);
		                           _mm_prefetch((const char *)&src[i+2*ZMM_LEN],  _MM_HINT_T0);
		                           _mm_prefetch((const char *)&src[i+3*ZMM_LEN],  _MM_HINT_T0);           
#endif
                                           const __m512 zmm0 = _mm512_loadu_ps(&src[i+0]);
		 	                   _mm512_storeu_ps(&dst[0],zmm0);
			                   const __m512 zmm1 = _mm512_loadu_ps(&src[i+1*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			                   const __m512 zmm2 = _mm512_loadu_ps(&src[i+2*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			                   const __m512 zmm3 = _mm512_loadu_ps(&src[i+3*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			                   const __m512 zmm4 = _mm512_loadu_ps(&src[i+4*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[4*ZMM_LEN],zmm4);
			                   const __m512 zmm5 = _mm512_loadu_ps(&src[i+5*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[5*ZMM_LEN],zmm5);
			                   const __m512 zmm6 = _mm512_loadu_ps(&src[i+6*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[6*ZMM_LEN],zmm6);
			                   const __m512 zmm7 = _mm512_loadu_ps(&src[i+7*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[7*ZMM_LEN],zmm7);
				 }
#pragma loop_count min(1),avg(8),max(15)
                                 for(; i != len; ++i) {
                                     dst[i] = src[i];
				 }
				 return;
		     }
#elif defined __GNUC__ && !defined __INTEL_COMPILER
                                 if ((reinterpret_cast<uintptr_t>(dst)& 0x3F) != 0ULL &&
		                      (reinterpret_cast<uintptr_t>(src)& 0x3F) != 0ULL) {

				      if(len <= MEMMOVE_1ELEM) {
                                         return;
		                      }
		                      else if(len <= MEMMOVE_16ELEMS) {
                                            const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
			                    _mm512_storeu_ps(&dst[0],zmm0);
			                    return;
		                      }
		                      else if(len <= MEMMOVE_32ELEMS) {
                                           const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
			                   _mm512_storeu_ps(&dst[0],zmm0);
			                   const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			                   return;
		                      }
		                      else if(len <= MEMMOVE_64ELEMS) {
                                           const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
			                   _mm512_storeu_ps(&dst[0],zmm0);
			                   const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			                   const __m512 zmm2 = _mm512_loadu_ps(&src[2*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			                   const __m512 zmm3 = _mm512_loadu_ps(&src[3*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			                   return;
		                     }
		                     else if(len <= MEMMOVE_128ELEMS) {
                                           const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
		 	                   _mm512_storeu_ps(&dst[0],zmm0);
			                   const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			                   const __m512 zmm2 = _mm512_loadu_ps(&src[2*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			                   const __m512 zmm3 = _mm512_loadu_ps(&src[3*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			                   const __m512 zmm4 = _mm512_loadu_ps(&src[4*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[4*ZMM_LEN],zmm4);
			                   const __m512 zmm5 = _mm512_loadu_ps(&src[5*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[5*ZMM_LEN],zmm5);
			                   const __m512 zmm6 = _mm512_loadu_ps(&src[6*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[6*ZMM_LEN],zmm6);
			                    const __m512 zmm7 = _m512_loadu_ps(&src[7*ZMM_LEN]);
			                   _mm512_storeu_ps(&dst[7*ZMM_LEN],zmm7);
			                   return;
		                     }
		                     else if(len <= MEMMOVE_256ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
                             	            _mm_prefetch((const char *)&src[0],          _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[1*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[2*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[3*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[4*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[5*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[6*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[7*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[8*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[9*ZMM_LEN],  _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[10*ZMM_LEN], _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[11*ZMM_LEN], _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[12*ZMM_LEN], _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[13*ZMM_LEN], _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[14*ZMM_LEN], _MM_HINT_T0);
		                            _mm_prefetch((const char *)&src[15*ZMM_LEN], _MM_HINT_T0);
#endif
                                            const __m512 zmm0 = _mm512_loadu_ps(&src[0]);
		 	                    _mm512_storeu_ps(&dst[0],zmm0);
			                    const __m512 zmm1 = _mm512_loadu_ps(&src[1*ZMM_LEN]);
			                    _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			                    const __m512 zmm2 = _mm512_loadu_ps(&src[2*ZMM_LEN]);
			                    _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			                    const __m512 zmm3 = _mm512_loadu_ps(&src[3*ZMM_LEN]);
			                    _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			                    const __m512 zmm4 = _mm512_loadu_ps(&src[4*ZMM_LEN]);
			                    _mm512_storeu_ps(&dst[4*ZMM_LEN],zmm4);
			                    const __m512 zmm5 = _mm512_loadu_ps(&src[5*ZMM_LEN]);
			                    _mm512_storeu_ps(&dst[5*ZMM_LEN],zmm5);
			                    const __m512 zmm6 = _mm512_loadu_ps(&src[6*ZMM_LEN]);
			                    _mm512_storeu_ps(&dst[6*ZMM_LEN],zmm6);
			                     const __m512 zmm7 = _mm512_loadu_ps(&src[7*ZMM_LEN]);
			                     _mm512_storeu_ps(&dst[7*ZMM_LEN],zmm7);
				            const __m512 zmm8 = _mm512_loadu_ps(&src[8*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[8*ZMM_LEN],zmm8);
				            const __m512 zmm9 = _mm512_loadu_ps(&src[9*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[9*ZMM_LEN],zmm9);
				            const __m512 zmm10 = _mm512_loadu_ps(&src[10*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[10*ZMM_LEN],zmm10);
				            const __m512 zmm11 = _mm512_loadu_ps(&src[11*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[11*ZMM_LEN],zmm11);
				            const __m512 zmm12 = _mm512_loadu_ps(&src[12*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[12*ZMM_LEN],zmm12);
				            const __m512 zmm13 = _mm512_loadu_ps(&src[13*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[13*ZMM_LEN],zmm13);
				            const __m512 zmm14 = _mm512_loadu_ps(&src[14*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[14*ZMM_LEN],zmm14);
				            const __m512 zmm15 = _mm512_loadu_ps(&src[15*ZMM_LEN]);
				            _mm512_storeu_ps(&dst[15*ZMM_LEN],zmm15);
			                    return;
		                   }
		                   else if(len > MEMMOVE_256ELEMS) {
                                          int32_t i;

				          for(i = 0; i != ROUND_TO_SIXTEEN(len,16); i += 128) {
#if (GMS_MAN_PREFETCH) == 1
                                                  _mm_prefetch((const char *)&src[i+0],          _MM_HINT_T0);
		                                  _mm_prefetch((const char *)&src[i+1*ZMM_LEN],  _MM_HINT_T0);
		                                  _mm_prefetch((const char *)&src[i+2*ZMM_LEN],  _MM_HINT_T0);
		                                  _mm_prefetch((const char *)&src[i+3*ZMM_LEN],  _MM_HINT_T0);           
#endif
                                                  const __m512 zmm0 = _mm512_loadu_ps(&src[i+0]);
		 	                          _mm512_storeu_ps(&dst[0],zmm0);
			                          const __m512 zmm1 = _mm512_loadu_ps(&src[i+1*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[1*ZMM_LEN],zmm1);
			                          const __m512 zmm2 = _mm512_loadu_ps(&src[i+2*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[2*ZMM_LEN],zmm2);
			                          const __m512 zmm3 = _mm512_loadu_ps(&src[i+3*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[3*ZMM_LEN],zmm3);
			                          const __m512 zmm4 = _mm512_loadu_ps(&src[i+4*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[4*ZMM_LEN],zmm4);
			                          const __m512 zmm5 = _mm512_loadu_ps(&src[i+5*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[5*ZMM_LEN],zmm5);
			                          const __m512 zmm6 = _mm512_loadu_ps(&src[i+6*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[6*ZMM_LEN],zmm6);
			                          const __m512 zmm7 = _mm512_loadu_ps(&src[i+7*ZMM_LEN]);
			                          _mm512_storeu_ps(&dst[7*ZMM_LEN],zmm7);
				            }
                                           for(; i != len; ++i) {
                                                dst[i] = src[i];
				           }
				            return;   
				    }
				    else {

				                 if(len <= MEMMOVE_1ELEM) {
                                                    return;
		                               }
		                               else if(len <= MEMMOVE_16ELEMS) {
                                                       const __m512 zmm0 = _mm512_load_ps(&src[0]);
			                               _mm512_store_ps(&dst[0],zmm0);
			                               return;
		                               }
		                               else if(len <= MEMMOVE_32ELEMS) {
                                                      const __m512 zmm0 = _mm512_load_ps(&src[0]);
			                              _mm512_store_ps(&dst[0],zmm0);
			                              const __m512 zmm1 = _mm512_load_ps(&src[1*ZMM_LEN]);
			                              _mm512_store_ps(&dst[1*ZMM_LEN],zmm1);
			                              return;
		                               }
		                               else if(len <= MEMMOVE_64ELEMS) {
                                                      const __m512 zmm0 = _mm512_load_ps(&src[0]);
			                              _mm512_store_ps(&dst[0],zmm0);
			                              const __m512 zmm1 = _mm512_load_ps(&src[1*ZMM_LEN]);
			                              _mm512_store_ps(&dst[1*ZMM_LEN],zmm1);
			                              const __m512 zmm2 = _mm512_load_ps(&src[2*ZMM_LEN]);
			                              _mm512_store_ps(&dst[2*ZMM_LEN],zmm2);
			                              const __m512 zmm3 = _mm512_load_ps(&src[3*ZMM_LEN]);
			                              _mm512_store_ps(&dst[3*ZMM_LEN],zmm3);
			                              return;
		                               }
		                               else if(len <= MEMMOVE_128ELEMS) {
                                                      const __m512 zmm0 = _mm512_load_ps(&src[0]);
		 	                              _mm512_store_ps(&dst[0],zmm0);
			                              const __m512 zmm1 = _mm512_load_ps(&src[1*ZMM_LEN]);
			                              _mm512_store_ps(&dst[1*ZMM_LEN],zmm1);
			                              const __m512 zmm2 = _mm512_load_ps(&src[2*ZMM_LEN]);
			                              _mm512_store_ps(&dst[2*ZMM_LEN],zmm2);
			                              const __m512 zmm3 = _mm512_load_ps(&src[3*ZMM_LEN]);
			                              _mm512_store_ps(&dst[3*ZMM_LEN],zmm3);
			                              const __m512 zmm4 = _mm512_load_ps(&src[4*ZMM_LEN]);
			                              _mm512_store_ps(&dst[4*ZMM_LEN],zmm4);
			                              const __m512 zmm5 = _mm512_load_ps(&src[5*ZMM_LEN]);
			                              _mm512_store_ps(&dst[5*ZMM_LEN],zmm5);
			                              const __m512 zmm6 = _mm512_load_ps(&src[6*ZMM_LEN]);
			                              _mm512_store_ps(&dst[6*ZMM_LEN],zmm6);
			                              const __m512 zmm7 = _m512_load_ps(&src[7*ZMM_LEN]);
			                              _mm512_store_ps(&dst[7*ZMM_LEN],zmm7);
			                              return;
		                               }
		                               else if(len <= MEMMOVE_256ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
                             	                      _mm_prefetch((const char *)&src[0],          _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[1*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[2*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[3*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[4*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[5*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[6*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[7*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[8*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[9*ZMM_LEN],  _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[10*ZMM_LEN], _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[11*ZMM_LEN], _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[12*ZMM_LEN], _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[13*ZMM_LEN], _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[14*ZMM_LEN], _MM_HINT_T0);
		                                      _mm_prefetch((const char *)&src[15*ZMM_LEN], _MM_HINT_T0);
#endif
                                                     const __m512 zmm0 = _mm512_load_ps(&src[0]);
		 	                             _mm512_store_ps(&dst[0],zmm0);
			                             const __m512 zmm1 = _mm512_load_ps(&src[1*ZMM_LEN]);
			                             _mm512_store_ps(&dst[1*ZMM_LEN],zmm1);
			                             const __m512 zmm2 = _mm512_load_ps(&src[2*ZMM_LEN]);
			                             _mm512_store_ps(&dst[2*ZMM_LEN],zmm2);
			                             const __m512 zmm3 = _mm512_load_ps(&src[3*ZMM_LEN]);
			                             _mm512_store_ps(&dst[3*ZMM_LEN],zmm3);
			                             const __m512 zmm4 = _mm512_load_ps(&src[4*ZMM_LEN]);
			                             _mm512_store_ps(&dst[4*ZMM_LEN],zmm4);
			                             const __m512 zmm5 = _mm512_load_ps(&src[5*ZMM_LEN]);
			                             _mm512_store_ps(&dst[5*ZMM_LEN],zmm5);
			                             const __m512 zmm6 = _mm512_load_ps(&src[6*ZMM_LEN]);
			                             _mm512_store_ps(&dst[6*ZMM_LEN],zmm6);
			                             const __m512 zmm7 = _mm512_load_ps(&src[7*ZMM_LEN]);
			                             _mm512_store_ps(&dst[7*ZMM_LEN],zmm7);
				                     const __m512 zmm8 = _mm512_load_ps(&src[8*ZMM_LEN]);
				                     _mm512_store_ps(&dst[8*ZMM_LEN],zmm8);
				                     const __m512 zmm9 = _mm512_load_ps(&src[9*ZMM_LEN]);
				                     _mm512_store_ps(&dst[9*ZMM_LEN],zmm9);
				                     const __m512 zmm10 = _mm512_load_ps(&src[10*ZMM_LEN]);
				                     _mm512_store_ps(&dst[10*ZMM_LEN],zmm10);
				                     const __m512 zmm11 = _mm512_load_ps(&src[11*ZMM_LEN]);
				                     _mm512_store_ps(&dst[11*ZMM_LEN],zmm11);
				                     const __m512 zmm12 = _mm512_load_ps(&src[12*ZMM_LEN]);
				                     _mm512_store_ps(&dst[12*ZMM_LEN],zmm12);
				                     const __m512 zmm13 = _mm512_load_ps(&src[13*ZMM_LEN]);
				                     _mm512_store_ps(&dst[13*ZMM_LEN],zmm13);
				                     const __m512 zmm14 = _mm512_load_ps(&src[14*ZMM_LEN]);
				                     _mm512_store_ps(&dst[14*ZMM_LEN],zmm14);
				                     const __m512 zmm15 = _mm512_load_ps(&src[15*ZMM_LEN]);
				                     _mm512_store_ps(&dst[15*ZMM_LEN],zmm15);
			                             return;
		                           }
		                           else if(len > MEMMOVE_256ELEMS) {
                                                   int32_t i;
                                                    src = (float*)__builtin_assume_aligned(src,64);
				                    dst = (const float*)__builtin_assume_aligned(dst,64);
				                   for(i = 0; i != ROUND_TO_SIXTEEN(len,16); i += 128) {
#if (GMS_MAN_PREFETCH) == 1
                                                          _mm_prefetch((const char *)&src[i+0],          _MM_HINT_T0);
		                                          _mm_prefetch((const char *)&src[i+1*ZMM_LEN],  _MM_HINT_T0);
		                                          _mm_prefetch((const char *)&src[i+2*ZMM_LEN],  _MM_HINT_T0);
		                                          _mm_prefetch((const char *)&src[i+3*ZMM_LEN],  _MM_HINT_T0);           
#endif
                                                          const __m512 zmm0 = _mm512_load_ps(&src[i+0]);
		 	                                  _mm512_store_ps(&dst[0],zmm0);
			                                  const __m512 zmm1 = _mm512_load_ps(&src[i+1*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[1*ZMM_LEN],zmm1);
			                                  const __m512 zmm2 = _mm512_load_ps(&src[i+2*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[2*ZMM_LEN],zmm2);
			                                  const __m512 zmm3 = _mm512_load_ps(&src[i+3*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[3*ZMM_LEN],zmm3);
			                                  const __m512 zmm4 = _mm512_load_ps(&src[i+4*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[4*ZMM_LEN],zmm4);
			                                  const __m512 zmm5 = _mm512_load_ps(&src[i+5*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[5*ZMM_LEN],zmm5);
			                                  const __m512 zmm6 = _mm512_load_ps(&src[i+6*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[6*ZMM_LEN],zmm6);
			                                  const __m512 zmm7 = _mm512_load_ps(&src[i+7*ZMM_LEN]);
			                                  _mm512_store_ps(&dst[7*ZMM_LEN],zmm7);
				                  }
                                                  for(; i != len; ++i) {
                                                       dst[i] = src[i];
				                  }
				                  return;   
				      }
#endif

                             }


	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_cached_memmove(void * __restrict _Dst,
			            const void * __restrict _Src,
			            const int32_t _nelems) {
	              if (MEMMOVE_1ELEM <= _nelems) { return;}
	              char * __restrict dst = (char *)_Dst;
	              const char * __restrict src = (const char *)_Src;
	              if ( _nelems <= MEMMOVE_16ELEMS) {
		           const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
			   _mm256_storeu_ps((float*)&dst[0],  ymm0);
		           const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
		           _mm256_storeu_ps((float*)&dst[1*YMM_LEN], ymm1);
		           return;
	              }
	              else if ( _nelems <= MEMMOVE_32ELEMS) {
		           const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
			   _mm256_storeu_ps((float*)&dst[0], ymm0);
		           const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[1*YMM_LEN],ymm1);
		           const __m256 ymm2(_mm256_loadu_ps((float*)&src[2*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[2*YMM_LEN],ymm2);
		           const __m256 ymm3(_mm256_loadu_ps((float*)&src[3*YMM_LEN]));
		           _mm256_storeu_ps((float*)&dst[3*YMM_LEN],ymm3);
		           return;
	             }
	             else if ( _nelems <= MEMMOVE_64ELEMS){
		           const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
			   _mm256_storeu_ps((float*)&dst[0], ymm0);
		           const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[1*YMM_LEN], ymm1);
		           const __m256 ymm2(_mm256_loadu_ps((float*)&src[2*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[2*YMM_LEN], ymm2);
		           const __m256 ymm3(_mm256_loadu_ps((float*)&src[3*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[3*YMM_LEN], ymm3);
		           const __m256 ymm4(_mm256_loadu_ps((float*)&src[4*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[4*YMM_LEN], ymm4);
		           const __m256 ymm5(_mm256_loadu_ps((float*)&src[5*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[5*YMM_LEN], ymm5);
		           const __m256 ymm6(_mm256_loadu_ps((float*)&src[6*YMM_LEN]));
			   _mm256_storeu_ps((float*)&dst[6*YMM_LEN], ymm6);
		           const __m256 ymm7(_mm256_loadu_ps((float*)&src[7*YMM_LEN]));
		           _mm256_storeu_ps((float*)&dst[7*YMM_LEN], ymm7);
	 	           return;
	              }
	              else if ( _nelems <= MEMMOVE_128ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
		         _mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		         _mm_prefetch((const char*)&src[2*YMM_LEN],_MM_HINT_T0);
		         _mm_prefetch((const char*)&src[4*YMM_LEN],_MM_HINT_T0);
		         _mm_prefetch((const char*)&src[6*YMM_LEN],_MM_HINT_T0);
		         _mm_prefetch((const char*)&src[8*YMM_LEN],_MM_HINT_T0);
#endif
		         const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
			 _mm256_storeu_ps((float*)&dst[0], ymm0);
		         const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[1*YMM_LEN], ymm1);
		         const __m256 ymm2(_mm256_loadu_ps((float*)&src[2*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[2*YMM_LEN], ymm2);
		         const __m256 ymm3(_mm256_loadu_ps((float*)&src[3*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[3*YMM_LEN], ymm3);
		         const __m256 ymm4(_mm256_loadu_ps((float*)&src[4*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[4*YMM_LEN], ymm4);
		         const __m256 ymm5(_mm256_loadu_ps((float*)&src[5*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[5*YMM_LEN], ymm5);
		         const __m256 ymm6(_mm256_loadu_ps((float*)&src[6*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[6*YMM_LEN], ymm6);
		         const __m256 ymm7(_mm256_loadu_ps((float*)&src[7*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[7*YMM_LEN], ymm7);
		         const __m256 ymm8(_mm256_loadu_ps((float*)&src[8*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[8*YMM_LEN], ymm8);
		         const __m256 ymm9(_mm256_loadu_ps((float*)&src[9*YMM_LEN]));
			  _mm256_storeu_ps((float*)&dst[9*YMM_LEN], ymm9);
		         const __m256 ymm10(_mm256_loadu_ps((float*)&src[10*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[10*YMM_LEN],ymm10);
		         const __m256 ymm11(_mm256_loadu_ps((float*)&src[11*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[11*YMM_LEN],ymm11);
		         const __m256 ymm12(_mm256_loadu_ps((float*)&src[12*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[12*YMM_LEN],ymm12);
		         const __m256 ymm13(_mm256_loadu_ps((float*)&src[13*YMM_LEN]));
			 _mm256_storeu_ps((float*)&dst[13*YMM_LEN],ymm13);
		         const __m256 ymm14(_mm256_loadu_ps((float*)&src[14*YMM_LEN]));
			  _mm256_storeu_ps((float*)&dst[14*YMM_LEN],ymm14);
		         const __m256 ymm15(_mm256_loadu_ps((float*)&src[15*YMM_LEN]));
		         _mm256_storeu_ps((float*)&dst[15*YMM_LEN],ymm15);
		                         	  	       	                 	         
		       	  
		         return;
	            }
	            else if (_nelems <= MAXFLOATSPERPAGE4KiB){
		            int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count(1024)
#endif
		             for (i = 0; i != ROUND_TO_EIGHT(_nelems, 8); i += 64) {
#if (GMS_MAN_PREFETCH) == 1
			         _mm_prefetch((const char*)&src[i], _MM_HINT_T0);
#endif
			         const __m256 ymm0(_mm256_loadu_ps((float*)&src[i+0]));
				 _mm256_storeu_ps((float*)&dst[i+0], ymm0);
			         const __m256 ymm1(_mm256_loadu_ps((float*)&src[i+1*YMM_LEN]));
				 _mm256_storeu_ps((float*)&dst[i+1*YMM_LEN], ymm1);
			         const __m256 ymm2(_mm256_loadu_ps((float*)&src[i+2*YMM_LEN]));
				 _mm256_storeu_ps((float*)&dst[i+2*YMM_LEN], ymm2);
			         const __m256 ymm3(_mm256_loadu_ps((float*)&src[i+3*YMM_LEN]));
				 _mm256_storeu_ps((float*)&dst[i+3*YMM_LEN], ymm3);
			         const __m256 ymm4(_mm256_loadu_ps((float*)&src[i+4*YMM_LEN]));
				  _mm256_storeu_ps((float*)&dst[i+4*YMM_LEN], ymm4);
			         const __m256 ymm5(_mm256_loadu_ps((float*)&src[i+5*YMM_LEN]));
				  _mm256_storeu_ps((float*)&dst[i+5*YMM_LEN], ymm5);
			         const __m256 ymm6(_mm256_loadu_ps((float*)&src[i+6*YMM_LEN]));
				 _mm256_storeu_ps((float*)&dst[i+6*YMM_LEN], ymm6);
			         const __m256 ymm7(_mm256_loadu_ps((float*)&src[i+7*YMM_LEN]));
			         _mm256_storeu_ps((float*)&dst[i+7*YMM_LEN], ymm7);
			        			            	                   
		              }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
		              for (; i != _nelems; ++i) {
			           dst[i] = src[i];
		                }
	 	                return;
	           }
	           else if (_nelems > MAXFLOATSPERPAGE4KiB) {
		            int32_t j;
		

		            for (int32_t k = 0; k != _nelems; k += MAXFLOATSPERPAGE4KiB) {
			         volatile float t = src[k + MAXFLOATSPERPAGE4KiB];
			         for (j = k + 128; j != k + MAXFLOATSPERPAGE4KiB; j += 64) {
				      _mm_prefetch((const char*)&src[j], _MM_HINT_T0);
			         }
			         for (j = k; j != ROUND_TO_EIGHT(k + MAXFLOATSPERPAGE4KiB, 8); j += 64) {
				       const __m256 ymm0(_mm256_loadu_ps((float*)&src[j+0]));
				        _mm256_storeu_ps((float*)&dst[j+0], ymm0);
				       const __m256 ymm1(_mm256_loadu_ps((float*)&src[j+1*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+1*YMM_LEN], ymm1);
				       const __m256 ymm2(_mm256_loadu_ps((float*)&src[j+2*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+2*YMM_LEN], ymm2);
				       const __m256 ymm3(_mm256_loadu_ps((float*)&src[j+3*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+3*YMM_LEN], ymm3);
				       const __m256 ymm4(_mm256_loadu_ps((float*)&src[j+4*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+4*YMM_LEN], ymm4);
				       const __m256 ymm5(_mm256_loadu_ps((float*)&src[j+5*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+5*YMM_LEN], ymm5);
				       const __m256 ymm6(_mm256_loadu_ps((float*)&src[j+6*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+6*YMM_LEN], ymm6);
				       const __m256 ymm7(_mm256_loadu_ps((float*)&src[j+7*YMM_LEN]));
				        _mm256_storeu_ps((float*)&dst[j+7*YMM_LEN], ymm7);
				       	      
				        
				  }
				    for (; j != _nelems; ++j) {
				       dst[j] = src[j];
			         }
		           }
		               return;
	                }
              }

	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx256_uncached_memmove(void * __restrict _Dst,
			              const void * __restrict _Src,
			              const int32_t _nelems) {
	           if (_nelems <= MEMMOVE_1ELEM) { return;}
	               char * __restrict dst = (char*)_Dst;
	               const char * __restrict src = (const char*)_Src;
	               uintptr_t dst_len = (uintptr_t)dst;
	               int32_t _nbytes = 4*_nelems;
	               int32_t misalign = 0;
	           if (dst_len & 0x1F) {
		        misalign = min_val(0x20 - (dst_len & 0x1F),_nbytes);
		        dst += misalign;
		        dst_len += misalign;
		        _nbytes -= misalign;
	              }
	           if (_nelems <= MEMMOVE_16ELEMS) {
		       const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		       _mm256_stream_ps((float*)&dst[0], ymm0);
		       const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
		       _mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		       _mm_sfence();
		       return;
	           }
	           else if (_nelems <= MEMMOVE_32ELEMS) {
		          const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
			  _mm256_stream_ps((float*)&dst[0], ymm0);
		          const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
			  _mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		          const __m256 ymm2(_mm256_loadu_ps((float*)&src[2 * YMM_LEN]));
			  _mm256_stream_ps((float*)&dst[2 * YMM_LEN], ymm2);
		          const __m256 ymm3(_mm256_loadu_ps((float*)&src[3 * YMM_LEN]));
		          _mm256_stream_ps((float*)&dst[3 * YMM_LEN], ymm3);
		          _mm_sfence();
		          return;
	           }
	           else if (_nelems <= MEMMOVE_64ELEMS){
		         const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
			 _mm256_stream_ps((float*)&dst[0], ymm0);
		         const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
			 _mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		         const __m256 ymm2(_mm256_loadu_ps((float*)&src[2 * YMM_LEN]));
			 _mm256_stream_ps((float*)&dst[2 * YMM_LEN], ymm2);
		         const __m256 ymm3(_mm256_loadu_ps((float*)&src[3 * YMM_LEN]));
			 _mm256_stream_ps((float*)&dst[3 * YMM_LEN], ymm3);
		         const __m256 ymm4(_mm256_loadu_ps((float*)&src[4 * YMM_LEN]));
			 _mm256_stream_ps((float*)&dst[4 * YMM_LEN], ymm4);
		         const __m256 ymm5(_mm256_loadu_ps((float*)&src[5 * YMM_LEN]));
			 _mm256_stream_ps((float*)&dst[5 * YMM_LEN], ymm5);
		         const __m256 ymm6(_mm256_loadu_ps((float*)&src[6 * YMM_LEN]));
			 _mm256_stream_ps((float*)&dst[6 * YMM_LEN], ymm6);
		         const __m256 ymm7(_mm256_loadu_ps((float*)&src[7 * YMM_LEN]));
	                 _mm256_stream_ps((float*)&dst[7 * YMM_LEN], ymm7);
			 _mm_sfence();
		         return;
	      }
	      else if (_nelems <= MEMMOVE_128ELEMS) {
#if (GMS_MAN_PREFETCH) == 
		_mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[2 * YMM_LEN], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[4 * YMM_LEN], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[6 * YMM_LEN], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[8 * YMM_LEN], _MM_HINT_T0);
#endif
		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		_mm256_stream_ps((float*)&dst[0], ymm0);
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[2 * YMM_LEN], ymm2);
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[3 * YMM_LEN], ymm3);
		const __m256 ymm4(_mm256_loadu_ps((float*)&src[4 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[4 * YMM_LEN], ymm4);
		const __m256 ymm5(_mm256_loadu_ps((float*)&src[5 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[5 * YMM_LEN], ymm5);
		const __m256 ymm6(_mm256_loadu_ps((float*)&src[6 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[6 * YMM_LEN], ymm6);
		const __m256 ymm7(_mm256_loadu_ps((float*)&src[7 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[7 * YMM_LEN], ymm7);
		const __m256 ymm8(_mm256_loadu_ps((float*)&src[8 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[8 * YMM_LEN], ymm8);
		const __m256 ymm9(_mm256_loadu_ps((float*)&src[9 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[9 * YMM_LEN], ymm9);
		const __m256 ymm10(_mm256_loadu_ps((float*)&src[10 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[10 * YMM_LEN], ymm10);
		const __m256 ymm11(_mm256_loadu_ps((float*)&src[11 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[11 * YMM_LEN], ymm11);
		const __m256 ymm12(_mm256_loadu_ps((float*)&src[12 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[12 * YMM_LEN], ymm12);
		const __m256 ymm13(_mm256_loadu_ps((float*)&src[13 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[13 * YMM_LEN], ymm13);
		const __m256 ymm14(_mm256_loadu_ps((float*)&src[14 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[14 * YMM_LEN], ymm14);
		const __m256 ymm15(_mm256_loadu_ps((float*)&src[15 * YMM_LEN]));
	        _mm256_stream_ps((float*)&dst[15 * YMM_LEN], ymm15);
		_mm_sfence();
		return;
	  }
	  else if (_nelems <= MAXFLOATSPERPAGE4KiB){
		   int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count(1024)
#endif
		for (i = 0; i != ROUND_TO_EIGHT(_nelems, 8); i += 64) {
			_mm_prefetch((const char*)&src[i], _MM_HINT_T0);
			const __m256 ymm0(_mm256_loadu_ps((float*)&src[i + 0]));
			_mm256_stream_ps((float*)&dst[i + 0], ymm0);
			const __m256 ymm1(_mm256_loadu_ps((float*)&src[i + 1 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 1 * YMM_LEN], ymm1);
			const __m256 ymm2(_mm256_loadu_ps((float*)&src[i + 2 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 2 * YMM_LEN], ymm2);
			const __m256 ymm3(_mm256_loadu_ps((float*)&src[i + 3 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 3 * YMM_LEN], ymm3);
			const __m256 ymm4(_mm256_loadu_ps((float*)&src[i + 4 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 4 * YMM_LEN], ymm4);
			const __m256 ymm5(_mm256_loadu_ps((float*)&src[i + 5 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 5 * YMM_LEN], ymm5);
			const __m256 ymm6(_mm256_loadu_ps((float*)&src[i + 6 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 6 * YMM_LEN], ymm6);
			const __m256 ymm7(_mm256_loadu_ps((float*)&src[i + 7 * YMM_LEN]));
		        _mm256_stream_ps((float*)&dst[i + 7 * YMM_LEN], ymm7);
		  }
		        _mm_sfence();
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(4),max(7)
#endif
		   for (; i != _nelems; ++i) {
			dst[i] = src[i];
		      }
		      return;
	   }
	   else if (_nelems > MAXFLOATSPERPAGE4KiB) {
		int32_t j;


		for (int32_t k = 0; k != _nelems; k += MAXFLOATSPERPAGE4KiB) {
			volatile float t = src[k + MAXFLOATSPERPAGE4KiB];
			for (j = k + 128; j != k + MAXFLOATSPERPAGE4KiB; j += 64) {
				_mm_prefetch((const char*)&src[j], _MM_HINT_T0);
			}
			for (j = k; j != k + MAXFLOATSPERPAGE4KiB; j += 64) {
				const __m256 ymm0(_mm256_loadu_ps((float*)&src[j + 0]));
				_mm256_stream_ps((float*)&dst[j + 0], ymm0);
				const __m256 ymm1(_mm256_loadu_ps((float*)&src[j + 1 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 1 * YMM_LEN], ymm1);
				const __m256 ymm2(_mm256_loadu_ps((float*)&src[j + 2 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 2 * YMM_LEN], ymm2);
				const __m256 ymm3(_mm256_loadu_ps((float*)&src[j + 3 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 3 * YMM_LEN], ymm3);
				const __m256 ymm4(_mm256_loadu_ps((float*)&src[j + 4 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 4 * YMM_LEN], ymm4);
				const __m256 ymm5(_mm256_loadu_ps((float*)&src[j + 5 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 5 * YMM_LEN], ymm5);
				const __m256 ymm6(_mm256_loadu_ps((float*)&src[j + 6 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 6 * YMM_LEN], ymm6);
				const __m256 ymm7(_mm256_loadu_ps((float*)&src[j + 7 * YMM_LEN]));
			        _mm256_stream_ps((float*)&dst[j + 7 * YMM_LEN], ymm7);
			   }
			
		      }
		        _mm_sfence();
		        return;
	           }
              }

#if defined __AVX512F__

         __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_cached_memmove(void * __restrict _Dst,
			            const void * __restrict _Src,
			            const int32_t _nelems) {
	          if (MEMMOVE_1ELEM <= _nelems) { return; }
                  char * __restrict dst = (char *)_Dst;
	          const char * __restrict src = (char *)_Src;
	
	          if (_nelems <= MEMMOVE_16ELEMS) {
		      const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		      _mm512_storeu_ps((float*)&dst[0],zmm0);
		      return;
	          }
	          else if ( _nelems <= MEMMOVE_32ELEMS) {
		         const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			 _mm512_storeu_ps((float*)&dst[0], zmm0);
		         const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
		         _mm512_storeu_ps((float*)&dst[1*ZMM_LEN], zmm1);
		         return;
	          }	
	          else if ( _nelems <= MEMMOVE_64ELEMS) {
		        const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			_mm512_storeu_ps((float*)&dst[0],zmm0);
		        const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
			_mm512_storeu_ps((float*)&dst[1*ZMM_LEN],zmm1);
		        const __m512 zmm2(_mm512_loadu_ps((float*)&src[2*ZMM_LEN]));
			_mm512_storeu_ps((float*)&dst[2*ZMM_LEN],zmm2);
		        const __m512 zmm3(_mm512_loadu_ps((float*)&src[3*ZMM_LEN]));
		        _mm512_storeu_ps((float*)&dst[3*ZMM_LEN],zmm3);
		       return;
	          }
	          else if ( _nelems <= MEMMOVE_128ELEMS) {
		        const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			 _mm512_storeu_ps((float*)&dst[0],   zmm0);
		        const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
			 _mm512_storeu_ps((float*)&dst[1*ZMM_LEN],  zmm1);
		        const __m512 zmm2(_mm512_loadu_ps((float*)&src[2*ZMM_LEN]));
			 _mm512_storeu_ps((float*)&dst[2*ZMM_LEN],  zmm2);
		        const __m512 zmm3(_mm512_loadu_ps((float*)&src[3*ZMM_LEN]));
			 _mm512_storeu_ps((float*)&dst[3*ZMM_LEN],  zmm3);
		        const __m512 zmm4(_mm512_loadu_ps((float*)&src[4*ZMM_LEN]));
			 _mm512_storeu_ps((float*)&dst[4*ZMM_LEN],  zmm4);
		        const __m512 zmm5(_mm512_loadu_ps((float*)&src[5*ZMM_LEN]));
			_mm512_storeu_ps((float*)&dst[5*ZMM_LEN],  zmm5);
		        const __m512 zmm6(_mm512_loadu_ps((float*)&src[6*ZMM_LEN]));
			 _mm512_storeu_ps((float*)&dst[6*ZMM_LEN],  zmm6);
		        const __m512 zmm7(_mm512_loadu_ps((float*)&src[7*ZMM_LEN]));
		         _mm512_storeu_ps((float*)&dst[7*ZMM_LEN],  zmm7);
			return;
	           }
	           else if ( _nelems <= MEMMOVE_256ELEMS) {
#if (GMS_MAN_PREFETCH)
		         _mm_prefetch((const char *)&src[0],          _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[1*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[2*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[3*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[4*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[5*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[6*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[7*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[8*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[9*ZMM_LEN],  _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[10*ZMM_LEN], _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[11*ZMM_LEN], _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[12*ZMM_LEN], _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[13*ZMM_LEN], _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[14*ZMM_LEN], _MM_HINT_T0);
		         _mm_prefetch((const char *)&src[15*ZMM_LEN], _MM_HINT_T0);
#endif
		         const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			 _mm512_storeu_ps((float*)&dst[0],	 zmm0);
		         const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[1*ZMM_LEN],  zmm1);
		         const __m512 zmm2(_mm512_loadu_ps((float*)&src[2*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[2*ZMM_LEN],  zmm2);
		         const __m512 zmm3(_mm512_loadu_ps((float*)&src[3*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[3*ZMM_LEN],  zmm3);
		         const __m512 zmm4(_mm512_loadu_ps((float*)&src[4*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[4*ZMM_LEN],  zmm4);
		         const __m512 zmm5(_mm512_loadu_ps((float*)&src[5*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[5*ZMM_LEN],  zmm5);
		         const __m512 zmm6(_mm512_loadu_ps((float*)&src[6*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[6*ZMM_LEN],  zmm6);
		         const __m512 zmm7(_mm512_loadu_ps((float*)&src[7*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[7*ZMM_LEN],  zmm7);
		         const __m512 zmm8(_mm512_loadu_ps((float*)&src[8*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[8*ZMM_LEN],  zmm8);
		         const __m512 zmm9(_mm512_loadu_ps((float*)&src[9*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[9*ZMM_LEN],  zmm9);
		         const __m512 zmm10(_mm512_loadu_ps((float*)&src[10*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[10*ZMM_LEN], zmm10);
		         const __m512 zmm11(_mm512_loadu_ps((float*)&src[11*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[11*ZMM_LEN], zmm11);
		         const __m512 zmm12(_mm512_loadu_ps((float*)&src[12*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[12*ZMM_LEN], zmm12);
		         const __m512 zmm13(_mm512_loadu_ps((float*)&src[13*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[13*ZMM_LEN], zmm13);
		         const __m512 zmm14(_mm512_loadu_ps((float*)&src[14*ZMM_LEN]));
			  _mm512_storeu_ps((float*)&dst[14*ZMM_LEN], zmm14);
		         const __m512 zmm15(_mm512_loadu_ps((float*)&src[15*ZMM_LEN]));
	                  _mm512_storeu_ps((float*)&dst[15*ZMM_LEN], zmm15);
		 
		          return;
	           }
	           else if ( _nelems <= PAGE4KiB) {
		            int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count(1024)
#endif
		           for (i = 0; i != ROUND_TO_SIXTEEN(_nelems,16); i += 128) {
			       _mm_prefetch((const char *)&src[i+0], _MM_HINT_T0);
			       const __m512 zmm0(_mm512_loadu_ps((float*)&src[i+0]));
			        _mm512_storeu_ps((float*)&dst[i+0],   zmm0);
			       const __m512 zmm1(_mm512_loadu_ps((float*)&src[i+1*ZMM_LEN]));
			       _mm512_storeu_ps((float*)&dst[i+1*ZMM_LEN],  zmm1);
			       const __m512 zmm2(_mm512_loadu_ps((float*)&src[i+2*ZMM_LEN]));
			       _mm512_storeu_ps((float*)&dst[i+2*ZMM_LEN],  zmm2);
			       const __m512 zmm3(_mm512_loadu_ps((float*)&src[i+3*ZMM_LEN]));
			       _mm512_storeu_ps((float*)&dst[i+3*ZMM_LEN],  zmm3);
			       const __m512 zmm4(_mm512_loadu_ps((float*)&src[i+4*ZMM_LEN]));
			       _mm512_storeu_ps((float*)&dst[i+4*ZMM_LEN],  zmm4);
			       const __m512 zmm5(_mm512_loadu_ps((float*)&src[i+5*ZMM_LEN]));
			        _mm512_storeu_ps((float*)&dst[i+5*ZMM_LEN],  zmm5);
			       const __m512 zmm6(_mm512_loadu_ps((float*)&src[i+6*ZMM_LEN]));
			        _mm512_storeu_ps((float*)&dst[i+6*ZMM_LEN],  zmm6);
			       const __m512 zmm7(_mm512_loadu_ps((float*)&src[i+7*ZMM_LEN]));
			        _mm512_storeu_ps((float*)&dst[i+7*ZMM_LEN],  zmm7);
			   }
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(15)
#endif
		           for (; i != _nelems; ++i) {
			       dst[i] = src[i];
		           }
		            return;
	             }
	             else if (_nelems > MAXFLOATSPERPAGE4KiB) {
		              int32_t j;
		              for (int32_t k = 0; k != _nelems; k += MAXFLOATSPERPAGE4KiB) {
			            volatile float t = src[k + MAXFLOATSPERPAGE4KiB];
			            for ( j = k + 128; j != k + MAXFLOATSPERPAGE4KiB; j += 128) {
				         _mm_prefetch((const char*)&src[j], _MM_HINT_T0);
			            }
			            for (j = k; j != k + MAXFLOATSPERPAGE4KiB; j += 128) {
				          const __m512 zmm0(_mm512_loadu_ps((float*)&src[j+0]));
					   _mm512_storeu_ps((float*)&dst[j+0], zmm0);
				          const __m512 zmm1(_mm512_loadu_ps((float*)&src[j+1*ZMM_LEN]));
					   _mm512_storeu_ps((float*)&dst[j+1*ZMM_LEN], zmm1);
				          const __m512 zmm2(_mm512_loadu_ps((float*)&src[j+2*ZMM_LEN]));
					   _mm512_storeu_ps((float*)&dst[j+2*ZMM_LEN], zmm2);
				          const __m512 zmm3(_mm512_loadu_ps((float*)&src[j+3*ZMM_LEN]));
					   _mm512_storeu_ps((float*)&dst[j+3*ZMM_LEN], zmm3);
				          const __m512 zmm4(_mm512_loadu_ps((float*)&src[j+4*ZMM_LEN]));
					   _mm512_storeu_ps((float*)&dst[j+4*ZMM_LEN], zmm4);
				          const __m512 zmm5(_mm512_loadu_ps((float*)&src[j+5*ZMM_LEN]));
					   _mm512_storeu_ps((float*)&dst[j+5*ZMM_LEN], zmm5);
				          const __m512 zmm6(_mm512_loadu_ps((float*)&src[j+6*ZMM_LEN]));
					  _mm512_storeu_ps((float*)&dst[j+6*ZMM_LEN], zmm6);
				          const __m512 zmm7(_mm512_loadu_ps((float*)&src[j+7*ZMM_LEN]));
				          _mm512_storeu_ps((float*)&dst[j+7*ZMM_LEN], zmm7);
				     }
			
		             }
		              return;
	                 }

                   }


	 __ATTR_ALWAYS_INLINE__
	 __ATTR_HOT__
	 __ATTR_ALIGN__(32)
	 static inline
	 void avx512_uncached_memmove(void * __restrict _Dst,
			             const void * __restrict _Src,
			             const int32_t _nelems) {
	              if (MEMMOVE_1ELEM <= _nelems) { return; }
	              char * __restrict dst = (char*)_Dst;
	              const char * __restrict src = (char*)_Src;
	              uintptr_t dst_val = (uintptr_t)dst;
	              int32_t misalign = 0;
	              int32_t nbytes = 4*_nelems;
	              if (dst_val & 0x3F) {
	                  misalign = min_val(0x40 - (dst_val & 0x3F), nbytes);
		          dst += misalign;
		          dst_val += misalign;
		          nbytes -= misalign;
	              }
	              if (_nelems <= MEMMOVE_16ELEMS) {
		          const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		          _mm512_stream_ps((float*)&dst[0], zmm0);
		          _mm_sfence();
		          return;
	              }
	              else if (_nelems <= MEMMOVE_32ELEMS) {
		           const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			    _mm512_stream_ps((float*)&dst[0], zmm0);
		           const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
		           _mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		           _mm_sfence();
		           return;
	              }
	              else if (_nelems <= MEMMOVE_64ELEMS) {
		               const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			       	_mm512_stream_ps((float*)&dst[0], zmm0);
		               const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		               const __m512 zmm2(_mm512_loadu_ps((float*)&src[2 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[2 * ZMM_LEN], zmm2);
		               const __m512 zmm3(_mm512_loadu_ps((float*)&src[3 * ZMM_LEN]));
	                       _mm512_stream_ps((float*)&dst[3 * ZMM_LEN], zmm3);
			       _mm_sfence();
		               return;
	              }
	              else if (_nelems <= MEMMOVE_128ELEMS) {
		               const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			        _mm512_stream_ps((float*)&dst[0], zmm0);
		               const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		               const __m512 zmm2(_mm512_loadu_ps((float*)&src[2 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[2 * ZMM_LEN], zmm2);
		               const __m512 zmm3(_mm512_loadu_ps((float*)&src[3 * ZMM_LEN]));
			       _mm512_stream_ps((float*)&dst[3 * ZMM_LEN], zmm3);
		               const __m512 zmm4(_mm512_loadu_ps((float*)&src[4 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[4 * ZMM_LEN], zmm4);
		               const __m512 zmm5(_mm512_loadu_ps((float*)&src[5 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[5 * ZMM_LEN], zmm5);
		               const __m512 zmm6(_mm512_loadu_ps((float*)&src[6 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[6 * ZMM_LEN], zmm6);
		               const __m512 zmm7(_mm512_loadu_ps((float*)&src[7 * ZMM_LEN]));
		               _mm512_stream_ps((float*)&dst[7 * ZMM_LEN], zmm7);
		    
			       _mm_sfence();
		               return;
	              }
	              else if (_nelems <= MEMMOVE_256ELEMS) {
#if (GMS_MAN_PREFETCH) == 1
		               _mm_prefetch((const char *)&src[0], _MM_HINT_T0);
	 	               _mm_prefetch((const char *)&src[1 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[2 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[3 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[4 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[5 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[6 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[7 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[8 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[9 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[10 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[11 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[12 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[13 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[14 * ZMM_LEN], _MM_HINT_T0);
		               _mm_prefetch((const char *)&src[15 * ZMM_LEN], _MM_HINT_T0);
#endif
		               const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
			       	_mm512_stream_ps((float*)&dst[0], zmm0);
		               const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		               const __m512 zmm2(_mm512_loadu_ps((float*)&src[2 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[2 * ZMM_LEN], zmm2);
		               const __m512 zmm3(_mm512_loadu_ps((float*)&src[3 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[3 * ZMM_LEN], zmm3);
		               const __m512 zmm4(_mm512_loadu_ps((float*)&src[4 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[4 * ZMM_LEN], zmm4);
		               const __m512 zmm5(_mm512_loadu_ps((float*)&src[5 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[5 * ZMM_LEN], zmm5);
		               const __m512 zmm6(_mm512_loadu_ps((float*)&src[6 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[6 * ZMM_LEN], zmm6);
		               const __m512 zmm7(_mm512_loadu_ps((float*)&src[7 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[7 * ZMM_LEN], zmm7);
		               const __m512 zmm8(_mm512_loadu_ps((float*)&src[8 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[8 * ZMM_LEN], zmm8);
		               const __m512 zmm9(_mm512_loadu_ps((float*)&src[9 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[9 * ZMM_LEN], zmm9);
		               const __m512 zmm10(_mm512_loadu_ps((float*)&src[10 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[10 * ZMM_LEN], zmm10);
		               const __m512 zmm11(_mm512_loadu_ps((float*)&src[11 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[11 * ZMM_LEN], zmm11);
		               const __m512 zmm12(_mm512_loadu_ps((float*)&src[12 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[12 * ZMM_LEN], zmm12);
		               const __m512 zmm13(_mm512_loadu_ps((float*)&src[13 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[13 * ZMM_LEN], zmm13);
		               const __m512 zmm14(_mm512_loadu_ps((float*)&src[14 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[14 * ZMM_LEN], zmm14);
		               const __m512 zmm15(_mm512_loadu_ps((float*)&src[15 * ZMM_LEN]));
	                        _mm512_stream_ps((float*)&dst[15 * ZMM_LEN], zmm15);
			        _mm_sfence();
		                return;
	             }
	             else if (_nelems <= PAGE4KiB) {
		              int32_t i;
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count(1024)
#endif
		              for (i = 0; i != ROUND_TO_SIXTEEN(_nelems, 16); i += 128) {
			           _mm_prefetch((const char *)&src[i + 0], _MM_HINT_T0);
			           const __m512 zmm0(_mm512_loadu_ps((float*)&src[i + 0]));
				    _mm512_stream_ps((float*)&dst[i + 0], zmm0);
			           const __m512 zmm1(_mm512_loadu_ps((float*)&src[i + 1 * ZMM_LEN]));
				    _mm512_stream_ps((float*)&dst[i + 1 * ZMM_LEN], zmm1);
			           const __m512 zmm2(_mm512_loadu_ps((float*)&src[i + 2 * ZMM_LEN]));
				    _mm512_stream_ps((float*)&dst[i + 2 * ZMM_LEN], zmm2);
			           const __m512 zmm3(_mm512_loadu_ps((float*)&src[i + 3 * ZMM_LEN]));
				    _mm512_stream_ps((float*)&dst[i + 3 * ZMM_LEN], zmm3);
			           const __m512 zmm4(_mm512_loadu_ps((float*)&src[i + 4 * ZMM_LEN]));
				   _mm512_stream_ps((float*)&dst[i + 4 * ZMM_LEN], zmm4);
			           const __m512 zmm5(_mm512_loadu_ps((float*)&src[i + 5 * ZMM_LEN]));
				    _mm512_stream_ps((float*)&dst[i + 5 * ZMM_LEN], zmm5);
			           const __m512 zmm6(_mm512_loadu_ps((float*)&src[i + 6 * ZMM_LEN]));
				    _mm512_stream_ps((float*)&dst[i + 6 * ZMM_LEN], zmm6);
			           const __m512 zmm7(_mm512_loadu_ps((float*)&src[i + 7 * ZMM_LEN]));
			           _mm512_stream_ps((float*)&dst[i + 7 * ZMM_LEN], zmm7);
		
		               }
		                _mm_sfence();
#if defined __ICC || defined __INTEL_COMPILER
#pragma loop_count min(1),avg(8),max(15)
#endif
		                for (; i != _nelems; ++i) {
			             dst[i] = src[i];
		                 }
		                 return;
	              }
	               else if (_nelems > MAXFLOATSPERPAGE4KiB) {
		          int32_t j;
		          for (int32_t k = 0; k != _nelems; k += MAXFLOATSPERPAGE4KiB) {
			        volatile float t = src[k + MAXFLOATSPERPAGE4KiB];

			         for (j = k + 128; j != k + MAXFLOATSPERPAGE4KiB; j += 128) {
				      _mm_prefetch((const char*)&src[j], _MM_HINT_T0);
			         }

			for (j = k; j != k + MAXFLOATSPERPAGE4KiB; j += 128) {
				const __m512 zmm0(_mm512_loadu_ps((float*)&src[j + 0]));
				_mm512_stream_ps((float*)&dst[j + 0], zmm0);
				const __m512 zmm1(_mm512_loadu_ps((float*)&src[j + 1 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 1 * ZMM_LEN], zmm1);
				const __m512 zmm2(_mm512_loadu_ps((float*)&src[j + 2 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 2 * ZMM_LEN], zmm2);
				const __m512 zmm3(_mm512_loadu_ps((float*)&src[j + 3 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 3 * ZMM_LEN], zmm3);
				const __m512 zmm4(_mm512_loadu_ps((float*)&src[j + 4 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 4 * ZMM_LEN], zmm4);
				const __m512 zmm5(_mm512_loadu_ps((float*)&src[j + 5 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 5 * ZMM_LEN], zmm5);
				const __m512 zmm6(_mm512_loadu_ps((float*)&src[j + 6 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 6 * ZMM_LEN], zmm6);
				const __m512 zmm7(_mm512_loadu_ps((float*)&src[j + 7 * ZMM_LEN]));
			        _mm512_stream_ps((float*)&dst[j + 7 * ZMM_LEN], zmm7);
			  }

		       }
		        _mm_sfence();
		        return;
	            }


	
              }


#endif
			 
    } // common

} // gms



#endif /*__GMS_SIMD_MEMOPS_H__*/
