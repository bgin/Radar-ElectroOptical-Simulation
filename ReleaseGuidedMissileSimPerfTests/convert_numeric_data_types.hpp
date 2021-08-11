


#ifndef __CONVERT_NUMERIC_DATA_TYPES_HPP__
#define __CONVERT_NUMERIC_DATA_TYPES_HPP__


#include <cstdint>
#include <omp.h>
#include <immintrin.h>


__attribute__((always_inline))
__attribute__((hot))
__attribute__((aligned(32)))
static inline
void
cvrt_int64_double_avx512_omp_ptr1(const uint64_t * __restrict __attribute__((aligned(64))) a,
				  double  * __restrict __attribute__((aligned(64))) b,
				  const int32_t data_len) {

     
     int32_t i, last_i;
     last_i = 0;
#if(defined __INTEL_COMPILER) || defined(__ICC)
    __assume_aligned(a,64);
    __assume_aligned(b,64);
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,64) default(none)   \
        lastprivate(last_i) private(i) shared(data_len,a,b)  \
	if(data_len>=100000)
     for(i = 0; i != ROUND_TO_EIGHT(data_len,8); i += 64) {
          last_i = i;
          _mm512_store_pd(&b[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+0])));
          _mm512_store_pd(&b[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+8])));
	  _mm512_store_pd(&b[i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+16])));
	  _mm512_store_pd(&b[i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+24])));
	  _mm512_store_pd(&b[i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+32])));
	  _mm512_store_pd(&b[i+40],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+40])));
	  _mm512_store_pd(&b[i+48],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+48])));
	  _mm512_store_pd(&b[i+56],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+56])));
     }
     for(; last_i != data_len; ++last_i) {
           b[last_i] = (double)a[last_i];
     }
}

__attribute__((always_inline))
__attribute__((hot))
__attribute__((aligned(32)))
static inline
void cvrt_double_float_avx512_omp_ptr1(const double * __restrict __attribute__((aligned(64))) a,
				       float * __restrict __attribute__((aligned(64))) b,
				       const int32_t data_len) {
        int32_t i,j, last_i;
	last_i = 0;
#if(defined __INTEL_COMPILER) || defined(__ICC)
        __assume_aligned(a,64);
	__assume_aligned(b,64);
#pragma code_align(32)
#endif
#pragma omp parallel for schedule(static,64) default(none) \
        lastprivate(last_i) private(i,j) shared(data_len,a,b) \
	if(data_len >= 100000)
     for(i = 0,j = 0; i != ROUND_TO_EIGHT(data_len,8); i +=64, j += 128) {
         last_i = i;
	   _mm512_store_ps(&b[j+0],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
          _mm512_store_ps(&b[j+16],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+8])));
	  _mm512_store_ps(&b[j+32],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+16])));
	  _mm512_store_ps(&b[j+48],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+24])));
	  _mm512_store_ps(&b[j+64],
	                  _mm512_castsipd_ps(_mm512_load_pd(&a[i+32])));
	  _mm512_store_ps(&b[j+80],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+40])));
	  _mm512_store_ps(&b[j+96],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+48])));
	  _mm512_store_ps(&b[j+112],
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+56])));
	 
     }
      for(; last_i != data_len; ++last_i) {
           b[last_i] = (double)a[last_i];
     }
}


__attribute__((always_inline))
__attribute__((hot))
__attribute__((aligned(32)))
static inline
void
cvrt_int64_double_avx512_omp_ptr4(const uint64_t * __restrict __attribute__((aligned(64))) a1,
                                  const uint64_t * __restrict __attribute__((aligned(64))) a2,
				  const uint64_t * __restrict __attribute__((aligned(64))) a3,
				  const uint64_t * __restrict __attribute__((aligned(64))) a4,
				  double  * __restrict __attribute__((aligned(64))) b1,
				  double  * __restrict __attribute__((aligned(64))) b2,
				  double  * __restrict __attribute__((aligned(64))) b3,
				  double  * __restrict __attribute__((aligned(64))) b4,
				  const int32_t data_len) {

        int32_t i, last_i;
	last_i = 0;
#if(defined __INTEL_COMPILER) || defined(__ICC)
        __assume_aligned(a1,64);
	__assume_aligned(a2,64);
	__assume_aligned(a3,64);
	__assume_aligned(a4,64);
	__assume_aligned(b1,64);
	__assume_aligned(b2,64);
	__assume_aligned(b3,64);
	__assume_aligned(b4,64)
#pragma code_align(32)
#endif
#pragma parallel for schedule(static,64) default(none) lastprivate(last_i) \
         private(i) shared(data_len,a1,a2,a3,a4,b1,b2,b3,b4) if(data_len>=10000)
         for(i = 0; i != ROUND_TO_EIGHT(data_len,8); i += 64) {
	      last_i = i;
          _mm512_store_pd(&b1[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+0])));
	  _mm512_store_pd(&b2[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+0])));
	  _mm512_store_pd(&b3[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+0])));
	  _mm512_store_pd(&b4[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+0])));
          _mm512_store_pd(&b1[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+8])));
	  _mm512_store_pd(&b2[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+8])));
	  _mm512_store_pd(&b3[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+8])));
	  _mm512_store_pd(&b4[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+8])));
	  _mm512_store_pd(&b1[i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+16])));
	  _mm512_store_pd(&b2[i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+16])));
	   _mm512_store_pd(&b3[i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+16])));
	   _mm512_store_pd(&b4[i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+16])));
	  _mm512_store_pd(&b1[i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+24])));
	   _mm512_store_pd(&b2[i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+24])));
	  _mm512_store_pd(&b3[i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+24])));
	   _mm512_store_pd(&b4[i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+24])));
	  _mm512_store_pd(&b1[i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+32])));
	   _mm512_store_pd(&b2[i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+32])));
	   _mm512_store_pd(&b3[i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+32])));
	   _mm512_store_pd(&b4[i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+32])));
	  _mm512_store_pd(&b1[i+40],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+40])));
	   _mm512_store_pd(&b2[i+40],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+40])));
	   _mm512_store_pd(&b3[i+40],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+40])));
	   _mm512_store_pd(&b4[i+40],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+40])));
	  _mm512_store_pd(&b1[i+48],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+48])));
	   _mm512_store_pd(&b2[i+48],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+48])));
	  _mm512_store_pd(&b2[i+48],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+48])));
	   _mm512_store_pd(&b3[i+48],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+48])));
	   _mm512_store_pd(&b4[i+48],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+48])));
	  _mm512_store_pd(&b1[i+56],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+56])));
	   _mm512_store_pd(&b2[i+56],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+56])));
	   _mm512_store_pd(&b3[i+56],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+56])));
	    _mm512_store_pd(&b4[i+56],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+56])));
           
	 }
	 for(; last_i != data_len; ++last_i) {
             b1[last_i] = (double)a1[last_i];
	     b2[last_i] = (double)a2[last_i];
	     b3[last_i] = (double)a3[last_i];
	     b4[last_i] = (double)a4[last_i];
	 }
  
}












#endif /*__CONVERT_NUMERIC_DATA_TYPES_HPP__*/
