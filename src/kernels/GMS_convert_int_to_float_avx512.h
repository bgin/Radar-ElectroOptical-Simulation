


#ifndef __GMS_CONVERT_INT_TO_FLOAT_AVX512_H__
#define __GMS_CONVERT_INT_TO_FLOAT_AVX512_H__


#include <cstdint>





__attribute__((hot))
__attribute__((aligned(32)))
void
cvrt_int64_double_avx512_omp_ptr1(uint64_t * __restrict __attribute__((aligned(64))) a,
				  double  * __restrict __attribute__((aligned(64))) b,
				  const int32_t data_len); 
				  
__attribute__((hot))
__attribute__((aligned(32)))

void cvrt_double_float_avx512_omp_ptr1(double * __restrict __attribute__((aligned(64))) a,
				       float * __restrict __attribute__((aligned(64))) b,
				       const int32_t data_len); 



__attribute__((hot))
__attribute__((aligned(32)))
void
cvrt_int64_double_avx512_omp_ptr4(uint64_t * __restrict __attribute__((aligned(64))) a1,
                                  uint64_t * __restrict __attribute__((aligned(64))) a2,
				  uint64_t * __restrict __attribute__((aligned(64))) a3,
				  uint64_t * __restrict __attribute__((aligned(64))) a4,
				  double  * __restrict __attribute__((aligned(64))) b1,
				  double  * __restrict __attribute__((aligned(64))) b2,
				  double  * __restrict __attribute__((aligned(64))) b3,
				  double  * __restrict __attribute__((aligned(64))) b4,
				  const int32_t data_len); 











#endif /*__GMS_CONVERT_INT_TO_FLOAT_AVX512_H__*/
