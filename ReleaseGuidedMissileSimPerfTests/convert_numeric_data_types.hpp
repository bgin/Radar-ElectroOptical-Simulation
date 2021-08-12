


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
cvrt_int64_double_avx512_omp_ptr1(uint64_t * __restrict __attribute__((aligned(64))) a,
				  double  * __restrict __attribute__((aligned(64))) b,
				  const int32_t data_len) {

     
     int32_t i, last_i;
     last_i = 0;

    if(data_len>=100000) {
#if(defined __INTEL_COMPILER) || defined(__ICC)
    __assume_aligned(a,64);
    __assume_aligned(b,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
    a = (uint64_t*)__builtin_assume_aligned(a,64);
    b = (double*)__builtin_assume_aligned(b,64)
#endif
#pragma omp parallel for schedule(static,64) default(none)   \
        lastprivate(last_i) private(i) shared(data_len,a,b)  
    for(i = 0; (i+63) < data_len; i += 64) {
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
      // last_i must be available here!!
    for(; (last_i+39) < data_len; last_i += 40) {
          _mm512_store_pd(&b[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+0])));
          _mm512_store_pd(&b[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+8])));
	  _mm512_store_pd(&b[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+16])));
	  _mm512_store_pd(&b[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+24])));
	  _mm512_store_pd(&b[last_i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+32])));
    }
    for(; (last_i+31) < data_len; last_i += 32) {
          _mm512_store_pd(&b[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+0])));
          _mm512_store_pd(&b[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+8])));
	  _mm512_store_pd(&b[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+16])));
	  _mm512_store_pd(&b[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+24])));
    }
    for(; (last_i+15) < data_len; last_i += 16) {
           _mm512_store_pd(&b[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+0])));
          _mm512_store_pd(&b[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+8])));
    }
    for(; (last_i+7) < data_len; last_i += 8) {
           _mm512_store_pd(&b[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[last_i+0])));
    }
    
     for(; (last_i+0) < data_len; last_i += 1) {
           b[last_i] = (double)a[last_i];
     }
   } else // data_len < 100000
#if(defined __INTEL_COMPILER) || defined(__ICC)
    __assume_aligned(a,64);
    __assume_aligned(b,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__ICC) || !defined(__INTEL_COMPILER))
    a = (uint64_t*)__builtin_assume_aligned(a,64);
    b = (double*)__builtin_assume_aligned(b,64)
#endif   
       for(i = 0; (i+63) < data_len; i += 64) {
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
      // last_i must be available here!!
    for(; (i+39) < data_len; i += 40) {
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
    }
    for(; (i+31) < data_len; i += 32) {
          _mm512_store_pd(&b[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+0])));
          _mm512_store_pd(&b[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+8])));
	  _mm512_store_pd(&b[i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+16])));
	  _mm512_store_pd(&b[i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+24])));
    }
    for(; (i+15) < data_len; i += 16) {
           _mm512_store_pd(&b[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+0])));
          _mm512_store_pd(&b[i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+8])));
    }
    for(; (i+7) < data_len; i += 8) {
           _mm512_store_pd(&b[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a[i+0])));
    }
    
     for(; (i+0) < data_len; i += 1) {
           b[i] = (double)a[i];
     }  
   }
}

__attribute__((always_inline))
__attribute__((hot))
__attribute__((aligned(32)))
static inline
void cvrt_double_float_avx512_omp_ptr1(double * __restrict __attribute__((aligned(64))) a,
				       float * __restrict __attribute__((aligned(64))) b,
				       const int32_t data_len) {
        int32_t i,j, last_i;
	const int32_t vsf_len = 16;
	//const int32_t vdf_len = 8;
	last_i = 0;
     if(data_len > 100000) {
#if(defined __INTEL_COMPILER) || defined(__ICC)
        __assume_aligned(a,64);
	__assume_aligned(b,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
        a = (double*)__builtin_assume_aligned(a,64);
	b = (float*)__builtin_assume_aligned(b,64);
#endif
#pragma omp parallel for schedule(static,64) default(none) \
        lastprivate(last_i) private(i,j) shared(data_len,a,b,vsf_len) 
        for(i = 0; (i+63) < data_len; i +=64) {
             last_i = i;
	      _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+8])));
	      _mm512_store_ps(b+2*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+16])));
	      _mm512_store_ps(b+3*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+24])));
	      _mm512_store_ps(b+4*vsf_len,
	                  _mm512_castsipd_ps(_mm512_load_pd(&a[i+32])));
	      _mm512_store_ps(b+5*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+40])));
	      _mm512_store_ps(b+6*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+48])));
	      _mm512_store_ps(b+7*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+56])));
	      b += 8*vsf_len;
         }
         for(; (last_i+39) < data_len; last_i += 40) {
              _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+8])));
	      _mm512_store_ps(b+2*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+16])));
	      _mm512_store_ps(b+3*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+24])));
	      _mm512_store_ps(b+4*vsf_len,
	                  _mm512_castsipd_ps(_mm512_load_pd(&a[last_i+32])));
	      b += 5*vsf_len;
	 }
	 for(; (last_i+31) < data_len; last_i += 32) {
              _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+8])));
	      _mm512_store_ps(b+2*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+16])));
	      _mm512_store_ps(b+3*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+24])));
	      b += 4*vsf_len;
	 }
	 for(; (last_i+15) < data_len; last_i += 16) {
              _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+8])));
	      b += 2*vsf_len;
	 }
	 for(; (last_i+7) < data_len; last_i += 8) {
               _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[last_i+0])));
	       b += 1*vsf_len;
	 }
         for(; (last_i+0) < data_len; last_i += 1) {
             b[last_i] = (double)a[last_i];
         }
    }
    else { // data_len < 100000
#if(defined __INTEL_COMPILER) || defined(__ICC)
        __assume_aligned(a,64);
	__assume_aligned(b,64);
#pragma code_align(32)
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
        a = (double*)__builtin_assume_aligned(a,64);
	b = (float*)__builtin_assume_aligned(b,64);
#endif
         for(i = 0; (i+63) < data_len; i +=64) {
           
	      _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+8])));
	      _mm512_store_ps(b+2*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+16])));
	      _mm512_store_ps(b+3*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+24])));
	      _mm512_store_ps(b+4*vsf_len,
	                  _mm512_castsipd_ps(_mm512_load_pd(&a[i+32])));
	      _mm512_store_ps(b+5*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+40])));
	      _mm512_store_ps(b+6*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+48])));
	      _mm512_store_ps(b+7*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+56])));
	      b += 8*vsf_len;
         }
         for(; (i+39) < data_len; i += 40) {
              _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+8])));
	      _mm512_store_ps(b+2*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+16])));
	      _mm512_store_ps(b+3*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+24])));
	      _mm512_store_ps(b+4*vsf_len,
	                  _mm512_castsipd_ps(_mm512_load_pd(&a[i+32])));
	      b += 5*vsf_len;
	 }
	 for(; (i+31) < data_len; i += 32) {
              _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+8])));
	      _mm512_store_ps(b+2*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+16])));
	      _mm512_store_ps(b+3*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+24])));
	      b += 4*vsf_len;
	 }
	 for(; (i+15) < data_len; i += 16) {
              _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
              _mm512_store_ps(b+1*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+8])));
	      b += 2*vsf_len;
	 }
	 for(; (i+7) < data_len; last_i += 8) {
               _mm512_store_ps(b+0*vsf_len,
	                  _mm512_castspd_ps(_mm512_load_pd(&a[i+0])));
	       b += 1*vsf_len;
	 }
         for(; (i+0) < data_len; last_i += 1) {
             b[i] = (double)a[i];
         }
  
    }
   
}


__attribute__((always_inline))
__attribute__((hot))
__attribute__((aligned(32)))
static inline
void
cvrt_int64_double_avx512_omp_ptr4(uint64_t * __restrict __attribute__((aligned(64))) a1,
                                  uint64_t * __restrict __attribute__((aligned(64))) a2,
				  uint64_t * __restrict __attribute__((aligned(64))) a3,
				  uint64_t * __restrict __attribute__((aligned(64))) a4,
				  double  * __restrict __attribute__((aligned(64))) b1,
				  double  * __restrict __attribute__((aligned(64))) b2,
				  double  * __restrict __attribute__((aligned(64))) b3,
				  double  * __restrict __attribute__((aligned(64))) b4,
				  const int32_t data_len) {

        int32_t i, last_i;
	last_i = 0;
      if(data_len>=100000) {
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
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
        a1 = (uint64_t*)__builtin_assume_aligned(a1,64);
	a2 = (uint64_t*)__builtin_assume_aligned(a2,64);
	a3 = (uint64_t*)__builtin_assume_aligned(a3,64);
	a4 = (uint64_t*)__builtin_assume_aligned(a4,64);
	b1 = (double*)__builtin_assume_aligned(b1,64);
	b2 = (double*)__builtin_assume_aligned(b2,64);
	b3 = (double*)__builtin_assume_aligned(b3,64);
	b4 = (double*)__builtin_assume_aligned(b4,64);
#endif
#pragma parallel for schedule(static,64) default(none) lastprivate(last_i) \
         private(i) shared(data_len,a1,a2,a3,a4,b1,b2,b3,b4) 
         for(i = 0; (i+63) < data_len; i += 64) {
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
         for(; (last_i+39); last_i += 40) {
               _mm512_store_pd(&b1[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+0])));
	       _mm512_store_pd(&b2[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+0])));
	       _mm512_store_pd(&b3[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+0])));
	       _mm512_store_pd(&b4[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+0])));
               _mm512_store_pd(&b1[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+8])));
	       _mm512_store_pd(&b2[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+8])));
	       _mm512_store_pd(&b3[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+8])));
	       _mm512_store_pd(&b4[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+8])));
	       _mm512_store_pd(&b1[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+16])));
	       _mm512_store_pd(&b2[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+16])));
	       _mm512_store_pd(&b3[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+16])));
	       _mm512_store_pd(&b4[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+16])));
	       _mm512_store_pd(&b1[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+24])));
	       _mm512_store_pd(&b2[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+24])));
	       _mm512_store_pd(&b3[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+24])));
	       _mm512_store_pd(&b4[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+24])));
	       _mm512_store_pd(&b1[last_i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+32])));
	       _mm512_store_pd(&b2[last_i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+32])));
	       _mm512_store_pd(&b3[last_i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+32])));
	       _mm512_store_pd(&b4[last_i+32],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+32])));
	 }
	 for(; (last_i+31) < data_len; last_i += 32) {
               _mm512_store_pd(&b1[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+0])));
	       _mm512_store_pd(&b2[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+0])));
	       _mm512_store_pd(&b3[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+0])));
	       _mm512_store_pd(&b4[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+0])));
               _mm512_store_pd(&b1[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+8])));
	       _mm512_store_pd(&b2[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+8])));
	       _mm512_store_pd(&b3[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+8])));
	       _mm512_store_pd(&b4[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+8])));
	       _mm512_store_pd(&b1[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+16])));
	       _mm512_store_pd(&b2[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+16])));
	       _mm512_store_pd(&b3[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+16])));
	       _mm512_store_pd(&b4[last_i+16],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+16])));
	       _mm512_store_pd(&b1[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+24])));
	       _mm512_store_pd(&b2[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+24])));
	       _mm512_store_pd(&b3[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+24])));
	       _mm512_store_pd(&b4[last_i+24],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+24])));
	 }
	 for(; (last_i+15) < data_len; last_i += 16) {
                _mm512_store_pd(&b1[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+0])));
	       _mm512_store_pd(&b2[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+0])));
	       _mm512_store_pd(&b3[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+0])));
	       _mm512_store_pd(&b4[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+0])));
               _mm512_store_pd(&b1[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+8])));
	       _mm512_store_pd(&b2[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+8])));
	       _mm512_store_pd(&b3[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+8])));
	       _mm512_store_pd(&b4[last_i+8],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+8])));
	 }
	 for(; (last_i+7) < data_len; last_i += 8) {
                _mm512_store_pd(&b1[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[last_i+0])));
	       _mm512_store_pd(&b2[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[last_i+0])));
	       _mm512_store_pd(&b3[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[last_i+0])));
	       _mm512_store_pd(&b4[last_i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[last_i+0])));
	 }
	 for(; (last_i+0) < data_len; last_i += 1) {
             b1[last_i] = (double)a1[last_i];
	     b2[last_i] = (double)a2[last_i];
	     b3[last_i] = (double)a3[last_i];
	     b4[last_i] = (double)a4[last_i];
	 }

     }
     else {  // data_len < 100000
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
#elif defined(__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
        a1 = (uint64_t*)__builtin_assume_aligned(a1,64);
	a2 = (uint64_t*)__builtin_assume_aligned(a2,64);
	a3 = (uint64_t*)__builtin_assume_aligned(a3,64);
	a4 = (uint64_t*)__builtin_assume_aligned(a4,64);
	b1 = (double*)__builtin_assume_aligned(b1,64);
	b2 = (double*)__builtin_assume_aligned(b2,64);
	b3 = (double*)__builtin_assume_aligned(b3,64);
	b4 = (double*)__builtin_assume_aligned(b4,64);
#endif
         for(i = 0; (i+63) < data_len; i += 64) {
	     
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
         for(; (i+39); i += 40) {
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
	 }
	 for(; (i+31) < data_len; i += 32) {
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
	 }
	 for(; (i+15) < data_len; i += 16) {
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
	 }
	 for(; (i+7) < data_len; i += 8) {
                _mm512_store_pd(&b1[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a1[i+0])));
	       _mm512_store_pd(&b2[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a2[i+0])));
	       _mm512_store_pd(&b3[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a3[i+0])));
	       _mm512_store_pd(&b4[i+0],
	                  _mm512_castsi512_pd(_mm512_load_epi64(&a4[i+0])));
	 }
	 for(; (i+0) < data_len; i += 1) {
             b1[i] = (double)a1[i];
	     b2[i] = (double)a2[i];
	     b3[i] = (double)a3[i];
	     b4[i] = (double)a4[i];
	 }

 
     }
  
}












#endif /*__CONVERT_NUMERIC_DATA_TYPES_HPP__*/
