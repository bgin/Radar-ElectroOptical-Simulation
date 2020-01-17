
#include <math.h>
#include "GMS_common.h"
#include "GMS_simd_defs.h"
#include "GMS_config.h"
#include "GMS_avxvecf32.h"
#include "GMS_avxc8f32.h"
//
//	Implementation
//


bool
gms::common
::approximately_equalf64(const double a,
			 const double b,
			 const double eps) {
	return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::essentialy_equalf64(const double a,
		      const double b,
		      const double eps) {
	return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::definitely_greaterf64(const double a,
			const double b,
			const double eps) {
	return fabs(a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::definitely_lessf64(const double a,
		     const double b,
		     const double eps) {
	return fabs(b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::approximately_equalf32(const float a,
			 const float b,
			 const float eps) {
	return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::essentialy_equalf32(const float a,
		      const float b,
		      const float eps) {
	return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::definitely_greaterf32(const double a,
			const double b,
			const double eps) {
	return fabs(a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * eps);
}

bool
gms::common
::Is_ptr_aligned32(const double * __restrict x) {
	if ((reinterpret_cast<uintptr_t>(x)& 0x1F) == 0ULL) {
		return (true);
	}
	else { return (false); }
}

bool
gms::common
::Is_ptr_aligned32(const int64_t * __restrict x) {
	if ((reinterpret_cast<uintptr_t>(x) & 0x1F) == 0ULL) {
		return (true);
	}
	else { return (false); }
}

bool
gms::common::
Is_ptr_aligned64(const double * __restrict x) {
	if ((reinterpret_cast<uintptr_t>(x)& 0x3F) == 0ULL) {
		return (true);
	}
	else { return (false); }
}

bool
gms::common
::Is_ptr_aligned64(const int64_t * __restrict x) {
	if ((reinterpret_cast<uintptr_t>(x)& 0x3F) == 0ULL) {
		return (true);
	}
	else { return (false); }
}

void gms::common::
init_array(double* __restrict ar,
           const int64_t arlen,
           const double val) {
	__assume_aligned(ar,64);
	for (int64_t i = 0LL; i != arlen; ++i) {
		 ar[i] = val;
	}
}

void
gms::common
::init_array(float* __restrict ar,
	     const int64_t arlen,
	     const float val) {
  __assume_aligned(ar,64);
  for(int64_t i = 0LL; i != arlen; ++i) {
       ar[i] = val;
  }
}

void
gms::common::
init_array(  int64_t * __restrict in,
	     const int64_t len,
	     const int64_t val) {
	__assume_aligned(in,64);
	for (int64_t i = 0LL; i != len; ++i) {
		 in[i] = val;
	}
}

void
gms::common::
init_array(int32_t * __restrict in,
	   const int64_t len,
	   const int32_t val) {
	__assume_aligned(in,64);
	for (int64_t i = 0LL; i != len; ++i) {
		in[i] = val;
	}
}

void gms::common::
init_array(unsigned char* __restrict ar,
	   const int64_t arlen,
	   const unsigned char val) {
	__assume_aligned(ar,64);
	for (int64_t i = 0LL; i != arlen; ++i) {
		 ar[i] = val;
	}
}

void gms::common::
init_varray(double* __restrict ar,
	    const int64_t arlen,
	    const double val) {
	if ((reinterpret_cast<uintptr_t>(ar)& 0x1F) != 0ULL) {
		for (int64_t i = 0LL; i != arlen; i += 8LL) {
			_mm256_storeu_pd(&ar[i + 0LL], _mm256_set1_pd(val));
			_mm256_storeu_pd(&ar[i + 4LL], _mm256_set1_pd(val));
		}
	}
	else {
	     for (int64_t i = 0LL; i != arlen; i += 8LL) {
		     _mm256_store_pd(&ar[i + 0LL], _mm256_set1_pd(val));
		     _mm256_store_pd(&ar[i + 4LL], _mm256_set1_pd(val));
	  }
   }
}

void
gms::common
::init_varray(float* __restrict ar,
	      const int64_t arlen,
	      const float val) {
  if((reinterpret_cast<uintptr_t>(ar) & 0x1F) != 0ULL) {
    for(int64_t i = 0LL; i != arlen; i += 16) {
      _mm256_storeu_ps(&ar[i+0LL], _mm256_set1_ps(val));
      _mm256_storeu_ps(&ar[i+8LL], _mm256_set1_ps(val));
    }
  }
  else {
    for(int64_t i = 0LL; i != arlen; i += 16) {
      _mm256_store_ps(&ar[i+0LL], _mm256_set1_ps(val));
      _mm256_store_ps(&ar[i+8LL], _mm256_set1_ps(val));
    }
  }
}

void
gms::common::
avxvec8_init_unroll2x(AVXVec8 * __restrict vec8,
		      const int64_t len,
		      const AVXVec8 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,64);
#elif defined __INTEL_COMPILER
       assume_aligned(vec8,64)
#pragma vector always
#endif
       for(int64_t i = 0; i != len-1LL; i += 2LL) {
           vec8[i+0LL] = v;
	   vec8[i+1LL] = v;
       }
}

void
gms::common::
avxvec8_init_unroll4x(AVXVec8 * __restrict vec8,
                      const int64_t len,
		      const AVXVec8 v) {
#if defined __GNUC__
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,64);
#elif defined __INTEL_COMPILER
       assume_aligned(vec8,64)
#pragma vector always
#endif
       for(int64_t i = 0; i != len-3LL; i += 4LL) {
           vec8[i+0LL] = v;
	   vec8[i+1LL] = v;
	   vec8[i+2LL] = v;
	   vec8[i+3LL] = v;
       }
}

void
gms::common::
avxvec8_init_unroll8x(AVXVec8 * __restrict vec8,
                      const int64_t len,
		      const AVXVec8 v) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
       vec8 = (AVXVec8*)__builtin_assume_aligned(vec8,64);
#elif defined __INTEL_COMPILER
       __assume_aligned(vec8,64)
#pragma vector always
#endif
        for(int64_t i = 0LL; i != len-7LL; i += 8LL) {
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

void
gms::common::
avxvec8_copy_unroll2x(AVXVec8 * __restrict dst,
		      const AVXVec8 * __restrict src,
		      const int64_t len) {
   
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
     src = (AVXVec8*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd
#elif defined __INTEL_COMPILER
#pragma vector always
#endif
   for(int64_t i = 0LL; i != len-1LL; i += 2LL) {
       dst[i+0LL] = src[i+0LL];
       dst[i+1LL] = src[i+1LL];
   }
}

void
gms::common::
avxvec8_copy_unroll4x(AVXVec8 * __restrict dst,
		      const AVXVec8 * __restrict src,
		      const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
     src = (AVXVec8*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd
#elif defined __INTEL_COMPILER
#pragma vector always
#endif
   for(int64_t i = 0LL; i != len-3LL; i += 4LL) {
       dst[i+0LL] = src[i+0LL];
       dst[i+1LL] = src[i+1LL];
       dst[i+2LL] = src[i+2LL];
       dst[i+3LL] = src[i+3LL];
   }
}

void
gms::common::
avxvec8_copy_unroll8x(AVXVec8 * __restrict dst,
		      const AVXVec8 * __restrict src,
		      const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
     src = (AVXVec8*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
#if defined __GNUC__ && !defined __INTEL_COMPILER
#pragma omp simd aligned(dst,src:64)
#elif defined __INTEL_COMPILER
#pragma vector always
#endif
   for(int64_t i = 0LL; i != len-7LL; i += 8LL) {
        dst[i+0LL] = src[i+0LL];
        dst[i+1LL] = src[i+1LL];
        dst[i+2LL] = src[i+2LL];
        dst[i+3LL] = src[i+3LL];
	dst[i+4LL] = src[i+4LL];
	dst[i+5LL] = src[i+5LL];
	dst[i+6LL] = src[i+6LL];
	dst[i+7LL] = src[i+7LL];
   }
}

void
gms::common::
avxvec8_copy_from_r4(AVXVec8 * __restrict dst,
		     const float * __restrict src,
		     const int64_t len) {
     
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
     src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
    int64_t j;
    j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
    for(int64_t i = 0LL; i != len; i += 8LL) {
        const __m256 t = _mm256_load_ps(&src[i]);
	dst[j].m_v8 = t;
	++j;
    }

}

void
gms::common::
r4_copy_from_avxvec8(float * __restrict dst,
                     const AVXVec8 * __restrict src,
		     const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVXVec8*)__builtin_assume_aligned(dst,64);
     src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
    int64_t j;
    j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif    
    for(int64_t i = 0LL; i != len; i += 8LL) {
        const __m256 t = src[j].m_v8;
	_mm256_store_ps(&dst[i],t);
	++j;
    }

}

#if defined __AVX512F__
void
gms::common::
r4_copy_from_avx512vec16(float * __restrict dst,
                         const AVX512Vec16 * __restrict src,
			 const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVX512Vec16*)__builtin_assume_aligned(dst,64);
     src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
    int64_t j;
    j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif    
    for(int64_t i = 0LL; i != len; i += 16LL) {
         const __m512 t = src[j].m_v16;
	 _mm512_store_ps(&dst[i],t);
	 ++j;
    }
}

void
gms::common::
avx512vec16_copy_from_r4(AVX512Vec16 * __restrict dst,
                         const float * __restrict src,
			 const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst = (AVX512Vec16*)__builtin_assume_aligned(dst,64);
     src = (const float*)__builtin_assume_aligned(src,64);
#elif defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src,64);
#endif
    int64_t j;
    j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif    
    for(int64_t i = 0LL; i != len; i += 16LL) {
        const __m512 t = _mm512_load_ps(&src[i]);
	dst[j].m_v16 = t;
	++j;
    }
}

#endif

void
gms::common::
avxc8f32_copy_from_r4(AVXc8f32 * __restrict dst,
                      const float * __restrict src_re,
		      const float * __restrict src_im,
		      const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
     dst    = (AVXc8f32*)__builtin_assume_aligned(dst,64);
     src_re = (const float*)__builtin_assume_aligned(src_re,64);
     src_im = (const float*)__builtin_assume_aligned(src_im,64);
#elif defined __ICC || defined __INTEL_COMPILER
     __assume_aligned(dst,64);
     __assume_aligned(src_re,64);
     __assume_aligned(src_im,64);
#endif
     int64_t j;
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

void
gms::common::
r4_copy_from_avxc8f32(float * __restrict dst_re,
                      float * __restrict dst_im,
		      const AVXc8f32 * __restrict src,
		      const int64_t len) {
#if defined __GNUC__ && !defined __INTEL_COMPILER
    dst_re = (float*)__builtin_assume_aligned(dst_re,64);
    dst_im = (float*)__builtin_assume_aligned(dst_im,64);
    src    = (const AVXc8f32*)__builtin_assume_aligned(src,64);
#elif defined __ICC || defined __INTEL_COMPILER
    __assume_aligned(dst_re,64);
    __assume_aligned(dst_im,64);
    __assume_aligned(src,64);
#endif
    int64_t j;
    j = 0;
#if defined __ICC || defined __INTEL_COMPILER
#pragma code_align(32)
#endif
    for(int64_t i = 0LL; i != len; i += 8LL) {
        const __m256 tre = src[j].m_re;
	_mm256_store_ps(&dst_re[i],tre;
	const __m256 tim = src[j].m_im;
	_mm256_store_ps(&dst_im[i],tim);
    }
}

void
gms::common::
avx256_init_unroll2x_pd( double * __restrict v,
		      const int64_t vlen,
		      const double val) {
	
	__m256d vec = _mm256_set1_pd(val);
	int64_t i;
	
	if ((reinterpret_cast<uintptr_t>(v)& 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_FOUR(vlen, 4LL); i += 8LL) {
			_mm256_storeu_pd(&v[i + 0LL], vec);
			_mm256_storeu_pd(&v[i + 4LL], vec);
		}
	}
	else {
     	for (i = 0LL; i != ROUND_TO_FOUR(vlen, 4LL); i += 8LL) {
		    _mm256_store_pd(&v[i+0LL], vec);
		    _mm256_store_pd(&v[i+4LL], vec);
	   }
	}
	for (; i != vlen; ++i) {
		v[i] = val;
	}
}

void
gms::common
::avx256_init_unroll2x_ps( float * __restrict v,
			   const int64_t vlen,
			   const float val) {
          __m256 xmm0 = _mm256_set1_ps(val);
	 int64_t i;
	 if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {
	   for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 16LL) {
	     _mm256_storeu_ps(&v[i+0LL], xmm0);
	     _mm256_storeu_ps(&v[i+8LL], xmm0);
	   }
	 } else {
           for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 16LL) {
             _mm256_store_ps(&v[i+0LL], xmm0);
	     _mm256_store_ps(&v[i+8LL], xmm0);
	   }
	 }
	 for(; i != vlen; ++i) {
	   v[i] = val;
	 }
}

void
gms::common::
avx256_init_unroll4x_pd(double * __restrict v,
		     const int64_t vlen,
		     const double val) {
	
	 __m256d vec = _mm256_set1_pd(val);
	int64_t i;
	if ((reinterpret_cast<uintptr_t>(v)& 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_FOUR(vlen, 4LL); i += 16LL) {
			_mm256_storeu_pd(&v[i + 0LL], vec);
			_mm256_storeu_pd(&v[i + 4LL], vec);
			_mm256_storeu_pd(&v[i + 8LL], vec);
			_mm256_storeu_pd(&v[i + 12LL], vec);
		}
	}
	else {
	      for (i = 0LL; i != ROUND_TO_FOUR(vlen, 4LL); i += 16LL) {
		       _mm256_store_pd(&v[i+0LL],  vec);
		       _mm256_store_pd(&v[i+4LL],  vec);
		       _mm256_store_pd(&v[i+8LL],  vec);
		       _mm256_store_pd(&v[i+12LL], vec);
	      }
	   }
	for (; i != vlen; ++i) {
		v[i] = val;
	}
}

void
gms::common::
avx256_init_unroll4x_ps(float * __restrict v,
			const int64_t vlen,
			const float val) {
       __m256 xmm0 = _mm256_set1_ps(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {
	 for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 32LL) {
	     _mm256_storeu_ps(&v[i+0LL],  xmm0);
	     _mm256_storeu_ps(&v[i+8LL],  xmm0);
	     _mm256_storeu_ps(&v[i+16LL], xmm0);
	     _mm256_storeu_ps(&v[i+24LL], xmm0);
	 }
       } else {
	 for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 32LL) {
             _mm256_store_ps(&v[i+0LL],  xmm0);
	     _mm256_store_ps(&v[i+8LL],  xmm0);
	     _mm256_store_ps(&v[i+16LL], xmm0);
	     _mm256_store_ps(&v[i+24LL], xmm0);
	 }
       }
       for(; i != vlen; ++i) {
	 v[i] = val;
       }
}


void
gms::common::
avx256_init_unroll8x_pd(double * __restrict v,
			const int64_t vlen,
			const double val) {
	
        __m256d vec = _mm256_set1_pd(val);
	int64_t i; 
	if ((reinterpret_cast<uintptr_t>(v)& 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_FOUR(vlen, 4LL); i += 32LL) {
			    _mm256_storeu_pd(&v[i + 0LL], vec);
			    _mm256_storeu_pd(&v[i + 4LL], vec);
			    _mm256_storeu_pd(&v[i + 8LL], vec);
			    _mm256_storeu_pd(&v[i + 12LL], vec);
			    _mm256_storeu_pd(&v[i + 16LL], vec);
			    _mm256_storeu_pd(&v[i + 20LL], vec);
			    _mm256_storeu_pd(&v[i + 24LL], vec);
			    _mm256_storeu_pd(&v[i + 28LL], vec);
		}
	}
	else {
	      for ( i = 0LL; i != ROUND_TO_FOUR(vlen, 4LL); i += 32LL) {
		         _mm256_store_pd(&v[i+0LL], vec);
		         _mm256_store_pd(&v[i+4LL], vec);
		         _mm256_store_pd(&v[i+8LL], vec);
		         _mm256_store_pd(&v[i+12LL],vec);
		         _mm256_store_pd(&v[i+16LL],vec);
		         _mm256_store_pd(&v[i+20LL],vec);
		         _mm256_store_pd(&v[i+24LL],vec);
		         _mm256_store_pd(&v[i+28LL],vec);
	     }
	}
	for (; i != vlen; ++i)
		v[i] = val;
}

void
gms::common::
avx256_init_unroll8x_ps(float * __restrict v,
			const int64_t vlen,
			const float val) {
       __m256 xmm0 = _mm256_set1_ps(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x1F) != 0ULL) {
	 for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 64LL) {
	   _mm256_storeu_ps(&v[i+0LL],  xmm0);
	   _mm256_storeu_ps(&v[i+8LL],  xmm0);
	   _mm256_storeu_ps(&v[i+16LL], xmm0);
	   _mm256_storeu_ps(&v[i+24LL], xmm0);
	   _mm256_storeu_ps(&v[i+32LL], xmm0);
	   _mm256_storeu_ps(&v[i+40LL], xmm0);
	   _mm256_storeu_ps(&v[i+48LL], xmm0);
	   _mm256_storeu_ps(&v[i+56LL], xmm0);
	 }
       }else {
	 for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 64LL) {
           _mm256_store_ps(&v[i+0LL],  xmm0);
	   _mm256_store_ps(&v[i+8LL],  xmm0);
	   _mm256_store_ps(&v[i+16LL], xmm0);
	   _mm256_store_ps(&v[i+24LL], xmm0);
	   _mm256_store_ps(&v[i+32LL], xmm0);
	   _mm256_store_ps(&v[i+40LL], xmm0);
	   _mm256_store_ps(&v[i+48LL], xmm0);
	   _mm256_store_ps(&v[i+56LL], xmm0);
	 }
       }
       for(; i != vlen; ++i) {
	 v[i] = val;
       }
}

#if defined __AVX512F__

void
gms::common::
avx512_init_unroll2x_pd(double * __restrict v,
			const int64_t vlen,
			const double val) {
       __m512d xmm0 = _mm512_set1_pd(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	   for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 16LL) {
	        _mm512_storeu_pd(&v[i+0LL], xmm0);
	        _mm512_storeu_pd(&v[i+8LL], xmm0);
	 }
       } else {
	 for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 16LL) {
                _mm512_store_pd(&v[i+0LL], xmm0);
	        _mm512_store_pd(&v[i+8LL], xmm0);
	 }
       }
       for(; i != vlen; ++i) {
	 v[i] = val;
       }
}

void
gms::common::
avx512_init_unroll2x_ps(float * __restrict v,
			const int64_t vlen,
			const float val) {
       __m512 xmm0 = _mm512_set1_ps(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	   for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 32LL) {
	        _mm512_storeu_ps(&v[i+0LL],  xmm0);
	        _mm512_storeu_ps(&v[i+16LL], xmm0);
	 }
       } else {
	    for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 32LL) {
                _mm512_store_ps(&v[i+0LL],  xmm0);
	        _mm512_store_ps(&v[i+16LL], xmm0);
	 }
       }
       for(; i != vlen; ++i) {
	 v[i] = val;
       }
}

void
gms::common::
avx512_init_unroll4x_pd(double * __restrict v,
			const int64_t vlen,
			const double val) {
       __m512d zmm0 = _mm512_set1_pd(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	 for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 32LL) {
	     _mm512_storeu_pd(&v[i+0LL],  zmm0);
	     _mm512_storeu_pd(&v[i+8LL],  zmm0);
	     _mm512_storeu_pd(&v[i+16LL], zmm0);
	     _mm512_storeu_pd(&v[i+24LL], zmm0);
	 }
       } else {
           for(i = 0LL; i != ROUND_TO_EIGHT(vlen,8LL); i += 32LL) {
	     _mm512_store_pd(&v[i+0LL],  zmm0);
	     _mm512_store_pd(&v[i+8LL],  zmm0);
	     _mm512_store_pd(&v[i+16LL], zmm0);
	     _mm512_store_pd(&v[i+24LL], zmm0);
	 }
       }
       for(; i != vlen; ++i) {
	 v[i] = val;
       }
}

void
gms::common::
avx512_init_unroll4x_ps(float * __restrict v,
			const int64_t vlen,
			const float val) {
       __m512 zmm0 = _mm512_set1_ps(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x3F) ! = 0ULL) {
	 for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 64LL) {
	      _mm512_storeu_ps(&v[i+0LL], zmm0);
	      _mm512_storeu_ps(&v[i+16LL], zmm0);
	      _mm512_storeu_ps(&v[i+32LL], zmm0);
	      _mm512_storeu_ps(&v[i+48LL], zmm0);
	 }
       } else {
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
}

void
gms::common::
avx512_init_unroll8x_pd(double * __restrict v,
		        const int64_t vlen,
		        const double val) {
	
        __m512d vec = _mm512_set1_pd(val);
	int64_t i;
	if ((reinterpret_cast<uintptr_t>(v)& 0x3F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_EIGHT(vlen, 8LL); i += 64LL) {
			_mm512_storeu_pd(&v[i + 0LL], vec);
			_mm512_storeu_pd(&v[i + 8LL], vec);
			_mm512_storeu_pd(&v[i + 16LL], vec);
			_mm512_storeu_pd(&v[i + 24LL], vec);
			_mm512_storeu_pd(&v[i + 32LL], vec);
			_mm512_storeu_pd(&v[i + 40LL], vec);
			_mm512_storeu_pd(&v[i + 48LL], vec);
			_mm512_storeu_pd(&v[i + 56LL], vec);
		}
	}
	else {
	   for (i = 0LL; i != ROUND_TO_EIGHT(vlen, 8LL); i += 64LL) {
		   _mm512_store_pd(&v[i+0LL],vec);
		   _mm512_store_pd(&v[i+8LL],vec);
		   _mm512_store_pd(&v[i+16LL],vec);
		   _mm512_store_pd(&v[i+24LL],vec);
		   _mm512_store_pd(&v[i+32LL],vec);
		   _mm512_store_pd(&v[i+40LL],vec);
		   _mm512_store_pd(&v[i+48LL],vec);
		   _mm512_store_pd(&v[i+56LL],vec);
	   }
	}
	for (; i != vlen; ++i)
		v[i] = val;
}

void
gms::common::
avx512_init_unroll8x_ps(float * __restrict v,
			const int64_t vlen,
			const float val) {
       __m512 zmm0 = _mm512_set1_ps(val);
       int64_t i;
       if((reinterpret_cast<uintptr_t>(v) & 0x3F) != 0ULL) {
	 for(i = 0LL; i != ROUND_TO_SIXTEEN(vlen,16LL); i += 128) {
	      _mm512_storeu_ps(&v[i+0LL],  zmm0);
	      _mm512_storeu_ps(&v[i+16LL], zmm0);
	      _mm512_storeu_ps(&v[i+32LL], zmm0);
	      _mm512_storeu_ps(&v[i+48LL], zmm0);
	      _mm512_storeu_ps(&v[i+64LL], zmm0);
	      _mm512_storeu_ps(&v[i+80LL], zmm0);
	      _mm512_storeu_ps(&v[i+96LL], zmm0);
	      _mm512_storeu_ps(&v[i+112LL], zmm0);
	 }
       } else {
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
}

#endif

bool gms::common::
vzero_check(_In_ double* __restrict re,
			_In_ double* __restrict im,
			_In_ const int64_t len) {
	//using namespace lam::common;
	__m256d zmask1(_mm256_setzero_pd());
	__m256d zmask2(_mm256_setzero_pd());
	//__m256i vres(_mm256_setr_epi64x(0LL,0LL,0LL,0LL));
	const __m256d vzero = _mm256_set1_pd(0.0);
	
	bool bres = false;
	for (int64_t i = 1LL; i != len; i += 4LL) {
		zmask1 = _mm256_cmp_pd(vzero, _mm256_loadu_pd(&re[i]),_CMP_EQ_OQ);
		zmask2 = _mm256_cmp_pd(vzero, _mm256_loadu_pd(&im[i]),_CMP_EQ_OQ);
		if (_mm256_testc_pd(zmask1, _mm256_setzero_pd()) ||
			_mm256_testc_pd(zmask2, _mm256_setzero_pd())) {
			//vres = _mm256_setr_epi64x(1LL,1LL,1LL,1LL);
			goto found;
		}
	}
	found:
	bres = true;
	return (bres);
}

bool gms::common::
vzero_check1(_In_ const double* __restrict re,
			_In_ const int64_t len) {
	__m256d zmask1(_mm256_setzero_pd());
	//__m256i vres(_mm256_setr_epi64x(0LL,0LL,0LL,0LL));
	const __m256d vzero = _mm256_setzero_pd();
	bool bres = false;
	for (int64_t i = 0LL; i != len; i += 4LL) {
		zmask1 = _mm256_cmp_pd(vzero, _mm256_loadu_pd(&re[i]),_CMP_EQ_OQ);
		if (_mm256_testc_pd(zmask1, _mm256_setzero_pd())) {
			//vres = _mm256_setr_epi64x(1LL,1LL,1LL,1LL);
			goto found;
		}
	}
	
	found:
	bres = true;
	return (bres);
}

bool gms::common::
vzero_check3D(_In_ const double* __restrict data,
			  _In_ const int64_t nx,
			  _In_ const int64_t ny,
			  _In_ const int64_t nz) {
	__m256d zmask1(_mm256_setzero_pd());
	//__m256i vres(_mm256_setr_epi64x(0LL,0LL,0LL,0LL));
	const __m256d vzero = _mm256_setzero_pd();
	bool bres = false;
	for (int64_t i = 0; i != nx; ++i) {
		for (int64_t j = 0; j != ny; ++j) {
			for (int64_t k = 0; k != nz; k += 4LL) {
				zmask1 = _mm256_cmp_pd(vzero, _mm256_loadu_pd(&data[i+nx*j+ny*k]),_CMP_EQ_OQ);
				if (_mm256_testc_pd(zmask1, _mm256_setzero_pd())) {
					//vres = _mm256_setr_epi64x(1LL,1LL,1LL,1LL);
					goto found;
				}
			}
		}
	}
	found:
	bres = true;
	return (bres);
}

bool gms::common::
vzero_check3D(_In_ const double* __restrict data,
			  _In_ const int64_t its,
		      _In_ const int64_t ite,
              _In_ const int64_t kts,
              _In_ const int64_t kte,
              _In_ const int64_t jts,
              _In_ const int64_t jte) {
	__m256d zmask1(_mm256_setzero_pd());
	//__m256i vres(_mm256_setr_epi64x(0LL,0LL,0LL,0LL));
	const __m256d vzero = _mm256_setzero_pd();
	bool bres = false;
	int64_t j,k,i;
	for (j = jts; j != jte; ++j) {
		for (k = kts; k != kte; ++k) {
			for (i = its; i != ROUND_TO_FOUR(ite, 4ULL); i += 4ULL) {
				zmask1 = _mm256_cmp_pd(vzero, _mm256_loadu_pd(&data[j+jte*k+kte*i]),_CMP_EQ_OQ);
				if (_mm256_testc_pd(zmask1, _mm256_setzero_pd())) {
					goto found;
				}
			}
		}
	}
	// Scalar remainder
	for (; j != jte; ++j) {
		for (; k != kte; ++k) {
			for (; i != ite; ++i) {
				if (data[j + jte*k + kte*i] == 0.0) {
					goto found;
				}
			}
		}
	}
	found:
	bres = true;
	return (bres);
}



void
gms::common::
avx256_memcpy2x_nt_pd(double * __restrict dst,
		      const double * __restrict src,
		      const int64_t len) {

	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {
		for (i = 0ULL; i != ROUND_TO_FOUR(len, 4ULL); i += 8LL) {
			_mm256_stream_pd(&dst[i + 0LL], _mm256_loadu_pd(&src[i + 0LL]));
			_mm256_stream_pd(&dst[i + 0LL], _mm256_loadu_pd(&src[i + 4LL]));
		}
		_mm_sfence();
	}
         else   {
	
	   for (i = 0ULL; i != ROUND_TO_FOUR(len, 4ULL); i += 8LL) {
		    _mm256_stream_pd(&dst[i + 0LL], _mm256_load_pd(&src[i + 0LL]));
		    _mm256_stream_pd(&dst[i + 0LL], _mm256_load_pd(&src[i + 4LL]));
	   }
	  _mm_sfence();
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy2x_nt_ps(float * __restrict dst,
		      const float * __restrict src,
		      const int64_t len) {
       	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst) & 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src) & 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 16LL) {
			_mm256_stream_ps(&dst[i + 0LL], _mm256_loadu_ps(&src[i + 0LL]));
			_mm256_stream_ps(&dst[i + 8LL], _mm256_loadu_ps(&src[i + 8LL]));
		}
		_mm_sfence();
	}
         else   {
	
	   for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 16LL) {
		    _mm256_stream_ps(&dst[i + 0LL], _mm256_load_ps(&src[i + 0LL]));
		    _mm256_stream_ps(&dst[i + 8LL], _mm256_load_ps(&src[i + 8LL]));
	   }
	  _mm_sfence();
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy4x_nt_pd(double * __restrict dst,
		   const double * __restrict src,
		   const int64_t len) {

	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst) & 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src) & 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 16LL) {
			_mm256_stream_pd(&dst[i + 0LL],  _mm256_loadu_pd(&src[i + 0LL]));
			_mm256_stream_pd(&dst[i + 4LL],  _mm256_loadu_pd(&src[i + 4LL]));
			_mm256_stream_pd(&dst[i + 8LL],  _mm256_loadu_pd(&src[i + 8LL]));
			_mm256_stream_pd(&dst[i + 12LL], _mm256_loadu_pd(&src[i + 12LL]));
		}
		   _mm_sfence();
	} 
	else {
	    for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 16LL) {
		    _mm256_stream_pd(&dst[i + 0LL],  _mm256_load_pd(&src[i + 0LL]));
		    _mm256_stream_pd(&dst[i + 4LL],  _mm256_load_pd(&src[i + 4LL]));
		    _mm256_stream_pd(&dst[i + 8LL],  _mm256_load_pd(&src[i + 8LL]));
		    _mm256_stream_pd(&dst[i + 12LL], _mm256_load_pd(&src[i + 12LL]));
	   }
	      _mm_sfence();
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy4x_nt_ps(float * __restrict dst,
		      const float * __restrict src,
		      const int64_t len) {
      	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst) & 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src) & 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 32LL) {
			_mm256_stream_ps(&dst[i + 0LL],    _mm256_loadu_ps(&src[i + 0LL]));
			_mm256_stream_ps(&dst[i + 8LL],    _mm256_loadu_ps(&src[i + 8LL]));
			_mm256_stream_ps(&dst[i + 16LL],   _mm256_loadu_ps(&src[i + 16LL]));
			_mm256_stream_ps(&dst[i + 24LL],   _mm256_loadu_ps(&src[i + 24LL]));
		}
		   _mm_sfence();
	} 
	else {
	    for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 32LL) {
		    _mm256_stream_ps(&dst[i + 0LL],    _mm256_load_ps(&src[i + 0LL]));
		    _mm256_stream_ps(&dst[i + 8LL],    _mm256_load_ps(&src[i + 8LL]));
		    _mm256_stream_ps(&dst[i + 16LL],   _mm256_load_ps(&src[i + 16LL]));
		    _mm256_stream_ps(&dst[i + 24LL],   _mm256_load_ps(&src[i + 24LL]));
	   }
	      _mm_sfence();
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy8x_nt_pd(double * __restrict dst,
		      const double * __restrict src,
		      const int64_t len) {

	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
		(reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 32LL) {
			_mm256_stream_pd(&dst[i + 0LL],  _mm256_loadu_pd(&src[i + 0LL]));
			_mm256_stream_pd(&dst[i + 4LL],  _mm256_loadu_pd(&src[i + 4LL]));
			_mm256_stream_pd(&dst[i + 8LL],  _mm256_loadu_pd(&src[i + 8LL]));
			_mm256_stream_pd(&dst[i + 12LL], _mm256_loadu_pd(&src[i + 12LL]));
			_mm256_stream_pd(&dst[i + 16LL], _mm256_loadu_pd(&src[i + 16LL]));
			_mm256_stream_pd(&dst[i + 20LL], _mm256_loadu_pd(&src[i + 20LL]));
			_mm256_stream_pd(&dst[i + 24LL], _mm256_loadu_pd(&src[i + 24LL]));
			_mm256_stream_pd(&dst[i + 28LL], _mm256_loadu_pd(&src[i + 28LL]));
		}
		_mm_sfence();
	}
	else {
	     for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 32LL) {
		     _mm256_stream_pd(&dst[i + 0LL],  _mm256_load_pd(&src[i+0LL]));
		     _mm256_stream_pd(&dst[i + 4LL],  _mm256_load_pd(&src[i+4LL]));
		     _mm256_stream_pd(&dst[i + 8LL],  _mm256_load_pd(&src[i+8LL]));
		     _mm256_stream_pd(&dst[i + 12LL], _mm256_load_pd(&src[i+12LL]));
		     _mm256_stream_pd(&dst[i + 16LL], _mm256_load_pd(&src[i+16LL]));
		     _mm256_stream_pd(&dst[i + 20LL], _mm256_load_pd(&src[i+20LL]));
		     _mm256_stream_pd(&dst[i + 24LL], _mm256_load_pd(&src[i+24LL]));
		     _mm256_stream_pd(&dst[i + 28LL], _mm256_load_pd(&src[i+28LL]));
	}
	_mm_sfence();
  }
	for (; i != len ; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy8x_nt_ps(float * __restrict dst,
		      const float * __restrict src,
		      const int64_t len) {
     	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 64LL) {
			_mm256_stream_ps(&dst[i + 0LL],  _mm256_loadu_ps(&src[i + 0LL]));
			_mm256_stream_ps(&dst[i + 8LL],  _mm256_loadu_ps(&src[i + 8LL]));
			_mm256_stream_ps(&dst[i + 16LL], _mm256_loadu_ps(&src[i + 16LL]));
			_mm256_stream_ps(&dst[i + 24LL], _mm256_loadu_ps(&src[i + 24LL]));
			_mm256_stream_ps(&dst[i + 32LL], _mm256_loadu_ps(&src[i + 32LL]));
			_mm256_stream_ps(&dst[i + 40LL], _mm256_loadu_ps(&src[i + 40LL]));
			_mm256_stream_ps(&dst[i + 48LL], _mm256_loadu_ps(&src[i + 48LL]));
			_mm256_stream_ps(&dst[i + 56LL], _mm256_loadu_ps(&src[i + 56LL]));
		}
		_mm_sfence();
	}
	else {
	     for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 64LL) {
		     _mm256_stream_ps(&dst[i + 0LL],  _mm256_load_ps(&src[i+0LL]));
		     _mm256_stream_ps(&dst[i + 8LL],  _mm256_load_ps(&src[i+8LL]));
		     _mm256_stream_ps(&dst[i + 16LL],  _mm256_load_ps(&src[i+16LL]));
		     _mm256_stream_ps(&dst[i + 24LL], _mm256_load_ps(&src[i+24LL]));
		     _mm256_stream_ps(&dst[i + 32LL], _mm256_load_ps(&src[i+32LL]));
		     _mm256_stream_ps(&dst[i + 40LL], _mm256_load_ps(&src[i+40LL]));
		     _mm256_stream_ps(&dst[i + 48LL], _mm256_load_ps(&src[i+48LL]));
		     _mm256_stream_ps(&dst[i + 56LL], _mm256_load_ps(&src[i+56LL]));
	}
	_mm_sfence();
  }
	for (; i != len ; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy2x_pd(double * __restrict dst,
		   const double * __restrict src,
		   const int64_t len) {

	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst) & 0x1F) != 0ULL &&
		(reinterpret_cast<uintptr_t>(src) & 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 8LL) {
			_mm256_storeu_pd(&dst[i + 0LL], _mm256_loadu_pd(&src[i + 0LL]));
			_mm256_storeu_pd(&dst[i + 4LL], _mm256_loadu_pd(&src[i + 4LL]));
		}
	}
	else {
	    for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 8LL) {
		     _mm256_store_pd(&dst[i + 0LL], _mm256_load_pd(&src[i+0LL]));
		     _mm256_store_pd(&dst[i + 4LL], _mm256_load_pd(&src[i+4LL]));
	    }
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common
avx256_memcpy2x_ps(float * __restrict dst,
		   const float * __restrict src,
		   const int64_t len) {
       	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst) & 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src) & 0x1F) != 0ULL) {
		for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 16LL) {
			_mm256_storeu_ps(&dst[i + 0LL], _mm256_loadu_ps(&src[i + 0LL]));
			_mm256_storeu_ps(&dst[i + 8LL], _mm256_loadu_ps(&src[i + 8LL]));
		}
	}
	else {
	    for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 16LL) {
		     _mm256_store_ps(&dst[i + 0LL], _mm256_load_ps(&src[i+0LL]));
		     _mm256_store_ps(&dst[i + 8LL], _mm256_load_ps(&src[i+8LL]));
	    }
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy4x_pd(double * __restrict dst,
		   const double * __restrict src,
		   const int64_t len) {

	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
		(reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {
	
	    for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 16LL) {
		      _mm256_storeu_pd(&dst[i + 0LL],  _mm256_loadu_pd(&src[i+0LL]));
		      _mm256_storeu_pd(&dst[i + 4LL],  _mm256_loadu_pd(&src[i+4LL]));
		      _mm256_storeu_pd(&dst[i + 8LL],  _mm256_loadu_pd(&src[i+8LL]));
		      _mm256_storeu_pd(&dst[i + 12LL], _mm256_loadu_pd(&src[i+12LL]));
	     }
	 }
	else {
		for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 16LL) {
			_mm256_store_pd(&dst[i + 0LL],  _mm256_load_pd(&src[i + 0LL]));
			_mm256_store_pd(&dst[i + 4LL],  _mm256_load_pd(&src[i + 4LL]));
			_mm256_store_pd(&dst[i + 8LL],  _mm256_load_pd(&src[i + 8LL]));
			_mm256_store_pd(&dst[i + 12LL], _mm256_load_pd(&src[i + 12LL]));
		}
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy4x_ps(float * __restrict dst,
		   const float * __restrict src,
		   const int64_t len) {
      	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
	    (reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {
	
	    for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 32LL) {
		      _mm256_storeu_ps(&dst[i + 0LL],   _mm256_loadu_ps(&src[i+0LL]));
		      _mm256_storeu_ps(&dst[i + 8LL],   _mm256_loadu_ps(&src[i+8LL]));
		      _mm256_storeu_ps(&dst[i + 16LL],  _mm256_loadu_ps(&src[i+16LL]));
		      _mm256_storeu_ps(&dst[i + 24LL],  _mm256_loadu_ps(&src[i+24LL]));
	     }
	 }
	else {
		for (i = 0LL; i != ROUND_TO_EIGHT(len, 8LL); i += 32LL) {
			_mm256_store_ps(&dst[i + 0LL],  _mm256_load_ps(&src[i + 0LL]));
			_mm256_store_ps(&dst[i + 8LL],  _mm256_load_ps(&src[i + 8LL]));
			_mm256_store_ps(&dst[i + 16LL],  _mm256_load_ps(&src[i + 16LL]));
			_mm256_store_ps(&dst[i + 24LL], _mm256_load_ps(&src[i + 24LL]));
		}
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common::
avx256_memcpy8x(double    *    __restrict dst,
		const double *    __restrict src,
		const int64_t len) {

	int64_t i;
	if ((reinterpret_cast<uintptr_t>(dst)& 0x1F) != 0ULL &&
		(reinterpret_cast<uintptr_t>(src)& 0x1F) != 0ULL) {
	
	       for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 32LL) {
		       _mm256_storeu_pd(&dst[i + 0LL], _mm256_loadu_pd(&src[i + 0LL]));
		       _mm256_storeu_pd(&dst[i + 4LL], _mm256_loadu_pd(&src[i + 4LL]));
		       _mm256_storeu_pd(&dst[i + 8LL], _mm256_loadu_pd(&src[i + 8LL]));
		       _mm256_storeu_pd(&dst[i + 12LL], _mm256_loadu_pd(&src[i + 12LL]));
		       _mm256_storeu_pd(&dst[i + 16LL], _mm256_loadu_pd(&src[i+16LL]));
		       _mm256_storeu_pd(&dst[i + 20LL], _mm256_loadu_pd(&src[i+20LL]));
		       _mm256_storeu_pd(&dst[i + 24LL], _mm256_loadu_pd(&src[i+24LL]));
		       _mm256_storeu_pd(&dst[i + 28LL], _mm256_loadu_pd(&src[i+28LL]));
	      }
	}
	else {
		for (i = 0LL; i != ROUND_TO_FOUR(len, 4LL); i += 32LL) {
			_mm256_store_pd(&dst[i + 0LL], _mm256_load_pd(&src[i + 0LL]));
			_mm256_store_pd(&dst[i + 4LL], _mm256_load_pd(&src[i + 4LL]));
			_mm256_store_pd(&dst[i + 8LL], _mm256_load_pd(&src[i + 8LL]));
			_mm256_store_pd(&dst[i + 12LL], _mm256_load_pd(&src[i + 12LL]));
			_mm256_store_pd(&dst[i + 16LL], _mm256_load_pd(&src[i + 16LL]));
			_mm256_store_pd(&dst[i + 20LL], _mm256_load_pd(&src[i + 20LL]));
			_mm256_store_pd(&dst[i + 24LL], _mm256_load_pd(&src[i + 24LL]));
			_mm256_store_pd(&dst[i + 28LL], _mm256_load_pd(&src[i + 28LL]));
		}
	}
	for (; i != len; ++i)
		dst[i] = src[i];
}

void
gms::common
::avx256_cached_memmove(void * __restrict _Dst,
			const void * __restrict _Src,
			const int32_t _nelems) {
	if (MEMMOVE_1ELEM <= _nelems) { return;}
	char * __restrict dst = (char *)_Dst;
	const char * __restrict src = (const char *)_Src;
	if ( _nelems <= MEMMOVE_16ELEMS) {
		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
		_mm256_storeu_ps((float*)&dst[0],  ymm0);
		_mm256_storeu_ps((float*)&dst[1*YMM_LEN], ymm1);
		return;
	}
	else if ( _nelems <= MEMMOVE_32ELEMS) {
		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2*YMM_LEN]));
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3*YMM_LEN]));
		_mm256_storeu_ps((float*)&dst[0], ymm0);
		_mm256_storeu_ps((float*)&dst[1*YMM_LEN],ymm1);
		_mm256_storeu_ps((float*)&dst[2*YMM_LEN],ymm2);
		_mm256_storeu_ps((float*)&dst[3*YMM_LEN],ymm3);
		return;
	}
	else if ( _nelems <= MEMMOVE_64ELEMS){
		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2*YMM_LEN]));
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3*YMM_LEN]));
		const __m256 ymm4(_mm256_loadu_ps((float*)&src[4*YMM_LEN]));
		const __m256 ymm5(_mm256_loadu_ps((float*)&src[5*YMM_LEN]));
		const __m256 ymm6(_mm256_loadu_ps((float*)&src[6*YMM_LEN]));
		const __m256 ymm7(_mm256_loadu_ps((float*)&src[7*YMM_LEN]));
		_mm256_storeu_ps((float*)&dst[0], ymm0);
		_mm256_storeu_ps((float*)&dst[1*YMM_LEN], ymm1);
		_mm256_storeu_ps((float*)&dst[2*YMM_LEN], ymm2);
		_mm256_storeu_ps((float*)&dst[3*YMM_LEN], ymm3);
		_mm256_storeu_ps((float*)&dst[4*YMM_LEN], ymm4);
		_mm256_storeu_ps((float*)&dst[5*YMM_LEN], ymm5);
		_mm256_storeu_ps((float*)&dst[6*YMM_LEN], ymm6);
		_mm256_storeu_ps((float*)&dst[7*YMM_LEN], ymm7);
		return;
	}
	else if ( _nelems <= MEMMOVE_128ELEMS) {

		_mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[2*YMM_LEN],_MM_HINT_T0);
		_mm_prefetch((const char*)&src[4*YMM_LEN],_MM_HINT_T0);
		_mm_prefetch((const char*)&src[6*YMM_LEN],_MM_HINT_T0);
		_mm_prefetch((const char*)&src[8*YMM_LEN],_MM_HINT_T0);

		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1*YMM_LEN]));
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2*YMM_LEN]));
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3*YMM_LEN]));
		const __m256 ymm4(_mm256_loadu_ps((float*)&src[4*YMM_LEN]));
		const __m256 ymm5(_mm256_loadu_ps((float*)&src[5*YMM_LEN]));
		const __m256 ymm6(_mm256_loadu_ps((float*)&src[6*YMM_LEN]));
		const __m256 ymm7(_mm256_loadu_ps((float*)&src[7*YMM_LEN]));
		const __m256 ymm8(_mm256_loadu_ps((float*)&src[8*YMM_LEN]));
		const __m256 ymm9(_mm256_loadu_ps((float*)&src[9*YMM_LEN]));
		const __m256 ymm10(_mm256_loadu_ps((float*)&src[10*YMM_LEN]));
		const __m256 ymm11(_mm256_loadu_ps((float*)&src[11*YMM_LEN]));
		const __m256 ymm12(_mm256_loadu_ps((float*)&src[12*YMM_LEN]));
		const __m256 ymm13(_mm256_loadu_ps((float*)&src[13*YMM_LEN]));
		const __m256 ymm14(_mm256_loadu_ps((float*)&src[14*YMM_LEN]));
		const __m256 ymm15(_mm256_loadu_ps((float*)&src[15*YMM_LEN]));
		_mm256_storeu_ps((float*)&dst[0], ymm0);
		_mm256_storeu_ps((float*)&dst[1*YMM_LEN], ymm1);
		_mm256_storeu_ps((float*)&dst[2*YMM_LEN], ymm2);
		_mm256_storeu_ps((float*)&dst[3*YMM_LEN], ymm3);
		_mm256_storeu_ps((float*)&dst[4*YMM_LEN], ymm4);
		_mm256_storeu_ps((float*)&dst[5*YMM_LEN], ymm5);
		_mm256_storeu_ps((float*)&dst[6*YMM_LEN], ymm6);
		_mm256_storeu_ps((float*)&dst[7*YMM_LEN], ymm7);
		_mm256_storeu_ps((float*)&dst[8*YMM_LEN], ymm8);
		_mm256_storeu_ps((float*)&dst[9*YMM_LEN], ymm9);
		_mm256_storeu_ps((float*)&dst[10*YMM_LEN],ymm10);
		_mm256_storeu_ps((float*)&dst[11*YMM_LEN],ymm11);
		_mm256_storeu_ps((float*)&dst[12*YMM_LEN],ymm12);
		_mm256_storeu_ps((float*)&dst[13*YMM_LEN],ymm13);
		_mm256_storeu_ps((float*)&dst[14*YMM_LEN],ymm14);
		_mm256_storeu_ps((float*)&dst[15*YMM_LEN],ymm15);
		return;
	}
	else if (_nelems <= MAXFLOATSPERPAGE4KiB){
		int32_t i;
		for (i = 0; i != ROUND_TO_EIGHT(_nelems, 8); i += 64) {
			_mm_prefetch((const char*)&src[i], _MM_HINT_T0);
			const __m256 ymm0(_mm256_loadu_ps((float*)&src[i+0]));
			const __m256 ymm1(_mm256_loadu_ps((float*)&src[i+1*YMM_LEN]));
			const __m256 ymm2(_mm256_loadu_ps((float*)&src[i+2*YMM_LEN]));
			const __m256 ymm3(_mm256_loadu_ps((float*)&src[i+3*YMM_LEN]));
			const __m256 ymm4(_mm256_loadu_ps((float*)&src[i+4*YMM_LEN]));
			const __m256 ymm5(_mm256_loadu_ps((float*)&src[i+5*YMM_LEN]));
			const __m256 ymm6(_mm256_loadu_ps((float*)&src[i+6*YMM_LEN]));
			const __m256 ymm7(_mm256_loadu_ps((float*)&src[i+7*YMM_LEN]));
			_mm256_storeu_ps((float*)&dst[i+0], ymm0);
			_mm256_storeu_ps((float*)&dst[i+1*YMM_LEN], ymm1);
			_mm256_storeu_ps((float*)&dst[i+2*YMM_LEN], ymm2);
			_mm256_storeu_ps((float*)&dst[i+3*YMM_LEN], ymm3);
			_mm256_storeu_ps((float*)&dst[i+4*YMM_LEN], ymm4);
			_mm256_storeu_ps((float*)&dst[i+5*YMM_LEN], ymm5);
			_mm256_storeu_ps((float*)&dst[i+6*YMM_LEN], ymm6);
			_mm256_storeu_ps((float*)&dst[i+7*YMM_LEN], ymm7);
		}
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
				const __m256 ymm1(_mm256_loadu_ps((float*)&src[j+1*YMM_LEN]));
				const __m256 ymm2(_mm256_loadu_ps((float*)&src[j+2*YMM_LEN]));
				const __m256 ymm3(_mm256_loadu_ps((float*)&src[j+3*YMM_LEN]));
				const __m256 ymm4(_mm256_loadu_ps((float*)&src[j+4*YMM_LEN]));
				const __m256 ymm5(_mm256_loadu_ps((float*)&src[j+5*YMM_LEN]));
				const __m256 ymm6(_mm256_loadu_ps((float*)&src[j+6*YMM_LEN]));
				const __m256 ymm7(_mm256_loadu_ps((float*)&src[j+7*YMM_LEN]));
				_mm256_storeu_ps((float*)&dst[j+0], ymm0);
				_mm256_storeu_ps((float*)&dst[j+1*YMM_LEN], ymm1);
				_mm256_storeu_ps((float*)&dst[j+2*YMM_LEN], ymm2);
				_mm256_storeu_ps((float*)&dst[j+3*YMM_LEN], ymm3);
				_mm256_storeu_ps((float*)&dst[j+4*YMM_LEN], ymm4);
				_mm256_storeu_ps((float*)&dst[j+5*YMM_LEN], ymm5);
				_mm256_storeu_ps((float*)&dst[j+6*YMM_LEN], ymm6);
				_mm256_storeu_ps((float*)&dst[j+7*YMM_LEN], ymm7);
			}
			for (; j != _nelems; ++j) {
				dst[j] = src[j];
			}
		}
		return;
	}
}

void
gms::common
::avx256_uncached_memmove(void * __restrict _Dst,
			  const void * __restrict _Src,
			  const int32_t _nelems) {
	if (_nelems < MEMMOVE_1ELEM) { return;}
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
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[0], ymm0);
		_mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MEMMOVE_32ELEMS) {
		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2 * YMM_LEN]));
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[0], ymm0);
		_mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		_mm256_stream_ps((float*)&dst[2 * YMM_LEN], ymm2);
		_mm256_stream_ps((float*)&dst[3 * YMM_LEN], ymm3);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MEMMOVE_64ELEMS){
		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2 * YMM_LEN]));
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3 * YMM_LEN]));
		const __m256 ymm4(_mm256_loadu_ps((float*)&src[4 * YMM_LEN]));
		const __m256 ymm5(_mm256_loadu_ps((float*)&src[5 * YMM_LEN]));
		const __m256 ymm6(_mm256_loadu_ps((float*)&src[6 * YMM_LEN]));
		const __m256 ymm7(_mm256_loadu_ps((float*)&src[7 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[0], ymm0);
		_mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		_mm256_stream_ps((float*)&dst[2 * YMM_LEN], ymm2);
		_mm256_stream_ps((float*)&dst[3 * YMM_LEN], ymm3);
		_mm256_stream_ps((float*)&dst[4 * YMM_LEN], ymm4);
		_mm256_stream_ps((float*)&dst[5 * YMM_LEN], ymm5);
		_mm256_stream_ps((float*)&dst[6 * YMM_LEN], ymm6);
		_mm256_stream_ps((float*)&dst[7 * YMM_LEN], ymm7);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MEMMOVE_128ELEMS) {

		_mm_prefetch((const char*)&src[0], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[2 * YMM_LEN], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[4 * YMM_LEN], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[6 * YMM_LEN], _MM_HINT_T0);
		_mm_prefetch((const char*)&src[8 * YMM_LEN], _MM_HINT_T0);

		const __m256 ymm0(_mm256_loadu_ps((float*)&src[0]));
		const __m256 ymm1(_mm256_loadu_ps((float*)&src[1 * YMM_LEN]));
		const __m256 ymm2(_mm256_loadu_ps((float*)&src[2 * YMM_LEN]));
		const __m256 ymm3(_mm256_loadu_ps((float*)&src[3 * YMM_LEN]));
		const __m256 ymm4(_mm256_loadu_ps((float*)&src[4 * YMM_LEN]));
		const __m256 ymm5(_mm256_loadu_ps((float*)&src[5 * YMM_LEN]));
		const __m256 ymm6(_mm256_loadu_ps((float*)&src[6 * YMM_LEN]));
		const __m256 ymm7(_mm256_loadu_ps((float*)&src[7 * YMM_LEN]));
		const __m256 ymm8(_mm256_loadu_ps((float*)&src[8 * YMM_LEN]));
		const __m256 ymm9(_mm256_loadu_ps((float*)&src[9 * YMM_LEN]));
		const __m256 ymm10(_mm256_loadu_ps((float*)&src[10 * YMM_LEN]));
		const __m256 ymm11(_mm256_loadu_ps((float*)&src[11 * YMM_LEN]));
		const __m256 ymm12(_mm256_loadu_ps((float*)&src[12 * YMM_LEN]));
		const __m256 ymm13(_mm256_loadu_ps((float*)&src[13 * YMM_LEN]));
		const __m256 ymm14(_mm256_loadu_ps((float*)&src[14 * YMM_LEN]));
		const __m256 ymm15(_mm256_loadu_ps((float*)&src[15 * YMM_LEN]));
		_mm256_stream_ps((float*)&dst[0], ymm0);
		_mm256_stream_ps((float*)&dst[1 * YMM_LEN], ymm1);
		_mm256_stream_ps((float*)&dst[2 * YMM_LEN], ymm2);
		_mm256_stream_ps((float*)&dst[3 * YMM_LEN], ymm3);
		_mm256_stream_ps((float*)&dst[4 * YMM_LEN], ymm4);
		_mm256_stream_ps((float*)&dst[5 * YMM_LEN], ymm5);
		_mm256_stream_ps((float*)&dst[6 * YMM_LEN], ymm6);
		_mm256_stream_ps((float*)&dst[7 * YMM_LEN], ymm7);
		_mm256_stream_ps((float*)&dst[8 * YMM_LEN], ymm8);
		_mm256_stream_ps((float*)&dst[9 * YMM_LEN], ymm9);
		_mm256_stream_ps((float*)&dst[10 * YMM_LEN], ymm10);
		_mm256_stream_ps((float*)&dst[11 * YMM_LEN], ymm11);
		_mm256_stream_ps((float*)&dst[12 * YMM_LEN], ymm12);
		_mm256_stream_ps((float*)&dst[13 * YMM_LEN], ymm13);
		_mm256_stream_ps((float*)&dst[14 * YMM_LEN], ymm14);
		_mm256_stream_ps((float*)&dst[15 * YMM_LEN], ymm15);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MAXFLOATSPERPAGE4KiB){
		int32_t i;
		for (i = 0; i != ROUND_TO_EIGHT(_nelems, 8); i += 64) {
			_mm_prefetch((const char*)&src[i], _MM_HINT_T0);
			const __m256 ymm0(_mm256_loadu_ps((float*)&src[i + 0]));
			const __m256 ymm1(_mm256_loadu_ps((float*)&src[i + 1 * YMM_LEN]));
			const __m256 ymm2(_mm256_loadu_ps((float*)&src[i + 2 * YMM_LEN]));
			const __m256 ymm3(_mm256_loadu_ps((float*)&src[i + 3 * YMM_LEN]));
			const __m256 ymm4(_mm256_loadu_ps((float*)&src[i + 4 * YMM_LEN]));
			const __m256 ymm5(_mm256_loadu_ps((float*)&src[i + 5 * YMM_LEN]));
			const __m256 ymm6(_mm256_loadu_ps((float*)&src[i + 6 * YMM_LEN]));
			const __m256 ymm7(_mm256_loadu_ps((float*)&src[i + 7 * YMM_LEN]));
			_mm256_stream_ps((float*)&dst[i + 0], ymm0);
			_mm256_stream_ps((float*)&dst[i + 1 * YMM_LEN], ymm1);
			_mm256_stream_ps((float*)&dst[i + 2 * YMM_LEN], ymm2);
			_mm256_stream_ps((float*)&dst[i + 3 * YMM_LEN], ymm3);
			_mm256_stream_ps((float*)&dst[i + 4 * YMM_LEN], ymm4);
			_mm256_stream_ps((float*)&dst[i + 5 * YMM_LEN], ymm5);
			_mm256_stream_ps((float*)&dst[i + 6 * YMM_LEN], ymm6);
			_mm256_stream_ps((float*)&dst[i + 7 * YMM_LEN], ymm7);
		}
		_mm_sfence();
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
				const __m256 ymm1(_mm256_loadu_ps((float*)&src[j + 1 * YMM_LEN]));
				const __m256 ymm2(_mm256_loadu_ps((float*)&src[j + 2 * YMM_LEN]));
				const __m256 ymm3(_mm256_loadu_ps((float*)&src[j + 3 * YMM_LEN]));
				const __m256 ymm4(_mm256_loadu_ps((float*)&src[j + 4 * YMM_LEN]));
				const __m256 ymm5(_mm256_loadu_ps((float*)&src[j + 5 * YMM_LEN]));
				const __m256 ymm6(_mm256_loadu_ps((float*)&src[j + 6 * YMM_LEN]));
				const __m256 ymm7(_mm256_loadu_ps((float*)&src[j + 7 * YMM_LEN]));
				_mm256_stream_ps((float*)&dst[j + 0], ymm0);
				_mm256_stream_ps((float*)&dst[j + 1 * YMM_LEN], ymm1);
				_mm256_stream_ps((float*)&dst[j + 2 * YMM_LEN], ymm2);
				_mm256_stream_ps((float*)&dst[j + 3 * YMM_LEN], ymm3);
				_mm256_stream_ps((float*)&dst[j + 4 * YMM_LEN], ymm4);
				_mm256_stream_ps((float*)&dst[j + 5 * YMM_LEN], ymm5);
				_mm256_stream_ps((float*)&dst[j + 6 * YMM_LEN], ymm6);
				_mm256_stream_ps((float*)&dst[j + 7 * YMM_LEN], ymm7);
			}
			
		}
		_mm_sfence();
		return;
	}
}

#if defined __AVX512F__
void
gms::common
::avx512_cached_memmove(void * __restrict _Dst,
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
		const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
		_mm512_storeu_ps((float*)&dst[0], zmm0);
		_mm512_storeu_ps((float*)&dst[1*ZMM_LEN], zmm1);
		return;
	}	
	 else if ( _nelems <= MEMMOVE_64ELEMS) {
		 const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		 const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
		 const __m512 zmm2(_mm512_loadu_ps((float*)&src[2*ZMM_LEN]));
		 const __m512 zmm3(_mm512_loadu_ps((float*)&src[3*ZMM_LEN]));
		 _mm512_storeu_ps((float*)&dst[0],zmm0);
		 _mm512_storeu_ps((float*)&dst[1*ZMM_LEN],zmm1);
		 _mm512_storeu_ps((float*)&dst[2*ZMM_LEN],zmm2);
		 _mm512_storeu_ps((float*)&dst[3*ZMM_LEN],zmm3);
		 return;
	 }
	 else if ( _nelems <= MEMMOVE_128ELEMS) {
		 const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		 const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
		 const __m512 zmm2(_mm512_loadu_ps((float*)&src[2*ZMM_LEN]));
		 const __m512 zmm3(_mm512_loadu_ps((float*)&src[3*ZMM_LEN]));
		 const __m512 zmm4(_mm512_loadu_ps((float*)&src[4*ZMM_LEN]));
		 const __m512 zmm5(_mm512_loadu_ps((float*)&src[5*ZMM_LEN]));
		 const __m512 zmm6(_mm512_loadu_ps((float*)&src[6*ZMM_LEN]));
		 const __m512 zmm7(_mm512_loadu_ps((float*)&src[7*ZMM_LEN]));
		 _mm512_storeu_ps((float*)&dst[0],   zmm0);
		 _mm512_storeu_ps((float*)&dst[1*ZMM_LEN],  zmm1);
		 _mm512_storeu_ps((float*)&dst[2*ZMM_LEN],  zmm2);
		 _mm512_storeu_ps((float*)&dst[3*ZMM_LEN],  zmm3);
		 _mm512_storeu_ps((float*)&dst[4*ZMM_LEN],  zmm4);
		 _mm512_storeu_ps((float*)&dst[5*ZMM_LEN],  zmm5);
		 _mm512_storeu_ps((float*)&dst[6*ZMM_LEN],  zmm6);
		 _mm512_storeu_ps((float*)&dst[7*ZMM_LEN],  zmm7);
		 return;
	 }
	 else if ( _nelems <= MEMMOVE_256ELEMS) {
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
		 const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		 const __m512 zmm1(_mm512_loadu_ps((float*)&src[1*ZMM_LEN]));
		 const __m512 zmm2(_mm512_loadu_ps((float*)&src[2*ZMM_LEN]));
		 const __m512 zmm3(_mm512_loadu_ps((float*)&src[3*ZMM_LEN]));
		 const __m512 zmm4(_mm512_loadu_ps((float*)&src[4*ZMM_LEN]));
		 const __m512 zmm5(_mm512_loadu_ps((float*)&src[5*ZMM_LEN]));
		 const __m512 zmm6(_mm512_loadu_ps((float*)&src[6*ZMM_LEN]));
		 const __m512 zmm7(_mm512_loadu_ps((float*)&src[7*ZMM_LEN]));
		 const __m512 zmm8(_mm512_loadu_ps((float*)&src[8*ZMM_LEN]));
		 const __m512 zmm9(_mm512_loadu_ps((float*)&src[9*ZMM_LEN]));
		 const __m512 zmm10(_mm512_loadu_ps((float*)&src[10*ZMM_LEN]));
		 const __m512 zmm11(_mm512_loadu_ps((float*)&src[11*ZMM_LEN]));
		 const __m512 zmm12(_mm512_loadu_ps((float*)&src[12*ZMM_LEN]));
		 const __m512 zmm13(_mm512_loadu_ps((float*)&src[13*ZMM_LEN]));
		 const __m512 zmm14(_mm512_loadu_ps((float*)&src[14*ZMM_LEN]));
		 const __m512 zmm15(_mm512_loadu_ps((float*)&src[15*ZMM_LEN]));
		 _mm512_storeu_ps((float*)&dst[0],	 zmm0);
		 _mm512_storeu_ps((float*)&dst[1*ZMM_LEN],  zmm1);
		 _mm512_storeu_ps((float*)&dst[2*ZMM_LEN],  zmm2);
		 _mm512_storeu_ps((float*)&dst[3*ZMM_LEN],  zmm3);
		 _mm512_storeu_ps((float*)&dst[4*ZMM_LEN],  zmm4);
		 _mm512_storeu_ps((float*)&dst[5*ZMM_LEN],  zmm5);
		 _mm512_storeu_ps((float*)&dst[6*ZMM_LEN],  zmm6);
		 _mm512_storeu_ps((float*)&dst[7*ZMM_LEN],  zmm7);
		 _mm512_storeu_ps((float*)&dst[8*ZMM_LEN],  zmm8);
		 _mm512_storeu_ps((float*)&dst[9*ZMM_LEN],  zmm9);
		 _mm512_storeu_ps((float*)&dst[10*ZMM_LEN], zmm10);
		 _mm512_storeu_ps((float*)&dst[11*ZMM_LEN], zmm11);
		 _mm512_storeu_ps((float*)&dst[12*ZMM_LEN], zmm12);
		 _mm512_storeu_ps((float*)&dst[13*ZMM_LEN], zmm13);
		 _mm512_storeu_ps((float*)&dst[14*ZMM_LEN], zmm14);
		 _mm512_storeu_ps((float*)&dst[15*ZMM_LEN], zmm15);
		 return;
	 }
	 else if ( _nelems <= PAGE4KiB) {
		 int32_t i;
		 for (i = 0; i != ROUND_TO_SIXTEEN(_nelems,16); i += 128) {
			 _mm_prefetch((const char *)&src[i+0], _MM_HINT_T0);
			 const __m512 zmm0(_mm512_loadu_ps((float*)&src[i+0]));
			 const __m512 zmm1(_mm512_loadu_ps((float*)&src[i+1*ZMM_LEN]));
			 const __m512 zmm2(_mm512_loadu_ps((float*)&src[i+2*ZMM_LEN]));
			 const __m512 zmm3(_mm512_loadu_ps((float*)&src[i+3*ZMM_LEN]));
			 const __m512 zmm4(_mm512_loadu_ps((float*)&src[i+4*ZMM_LEN]));
			 const __m512 zmm5(_mm512_loadu_ps((float*)&src[i+5*ZMM_LEN]));
			 const __m512 zmm6(_mm512_loadu_ps((float*)&src[i+6*ZMM_LEN]));
			 const __m512 zmm7(_mm512_loadu_ps((float*)&src[i+7*ZMM_LEN]));
			 _mm512_storeu_ps((float*)&dst[i+0],   zmm0);
			 _mm512_storeu_ps((float*)&dst[i+1*ZMM_LEN],  zmm1);
			 _mm512_storeu_ps((float*)&dst[i+2*ZMM_LEN],  zmm2);
			 _mm512_storeu_ps((float*)&dst[i+3*ZMM_LEN],  zmm3);
			 _mm512_storeu_ps((float*)&dst[i+4*ZMM_LEN],  zmm4);
			 _mm512_storeu_ps((float*)&dst[i+5*ZMM_LEN],  zmm5);
			 _mm512_storeu_ps((float*)&dst[i+6*ZMM_LEN],  zmm6);
			 _mm512_storeu_ps((float*)&dst[i+7*ZMM_LEN],  zmm7);
		 }
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
				 const __m512 zmm1(_mm512_loadu_ps((float*)&src[j+1*ZMM_LEN]));
				 const __m512 zmm2(_mm512_loadu_ps((float*)&src[j+2*ZMM_LEN]));
				 const __m512 zmm3(_mm512_loadu_ps((float*)&src[j+3*ZMM_LEN]));
				 const __m512 zmm4(_mm512_loadu_ps((float*)&src[j+4*ZMM_LEN]));
				 const __m512 zmm5(_mm512_loadu_ps((float*)&src[j+5*ZMM_LEN]));
				 const __m512 zmm6(_mm512_loadu_ps((float*)&src[j+6*ZMM_LEN]));
				 const __m512 zmm7(_mm512_loadu_ps((float*)&src[j+7*ZMM_LEN]));
				 _mm512_storeu_ps((float*)&dst[j+0], zmm0);
				 _mm512_storeu_ps((float*)&dst[j+1*ZMM_LEN], zmm1);
				 _mm512_storeu_ps((float*)&dst[j+2*ZMM_LEN], zmm2);
				 _mm512_storeu_ps((float*)&dst[j+3*ZMM_LEN], zmm3);
				 _mm512_storeu_ps((float*)&dst[j+4*ZMM_LEN], zmm4);
				 _mm512_storeu_ps((float*)&dst[j+5*ZMM_LEN], zmm5);
				 _mm512_storeu_ps((float*)&dst[j+6*ZMM_LEN], zmm6);
				 _mm512_storeu_ps((float*)&dst[j+7*ZMM_LEN], zmm7);
			 }
			
		 }
		 return;
	 }

}

void
gms::common
::avx512_uncached_memmove(void * __restrict _Dst,
			  const void * __restrict _Src,
			  int32_t _nelems) {
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
		const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
		_mm512_stream_ps((float*)&dst[0], zmm0);
		_mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MEMMOVE_64ELEMS) {
		const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
		const __m512 zmm2(_mm512_loadu_ps((float*)&src[2 * ZMM_LEN]));
		const __m512 zmm3(_mm512_loadu_ps((float*)&src[3 * ZMM_LEN]));
		_mm512_stream_ps((float*)&dst[0], zmm0);
		_mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		_mm512_stream_ps((float*)&dst[2 * ZMM_LEN], zmm2);
		_mm512_stream_ps((float*)&dst[3 * ZMM_LEN], zmm3);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MEMMOVE_128ELEMS) {
		const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
		const __m512 zmm2(_mm512_loadu_ps((float*)&src[2 * ZMM_LEN]));
		const __m512 zmm3(_mm512_loadu_ps((float*)&src[3 * ZMM_LEN]));
		const __m512 zmm4(_mm512_loadu_ps((float*)&src[4 * ZMM_LEN]));
		const __m512 zmm5(_mm512_loadu_ps((float*)&src[5 * ZMM_LEN]));
		const __m512 zmm6(_mm512_loadu_ps((float*)&src[6 * ZMM_LEN]));
		const __m512 zmm7(_mm512_loadu_ps((float*)&src[7 * ZMM_LEN]));
		_mm512_stream_ps((float*)&dst[0], zmm0);
		_mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		_mm512_stream_ps((float*)&dst[2 * ZMM_LEN], zmm2);
		_mm512_stream_ps((float*)&dst[3 * ZMM_LEN], zmm3);
		_mm512_stream_ps((float*)&dst[4 * ZMM_LEN], zmm4);
		_mm512_stream_ps((float*)&dst[5 * ZMM_LEN], zmm5);
		_mm512_stream_ps((float*)&dst[6 * ZMM_LEN], zmm6);
		_mm512_stream_ps((float*)&dst[7 * ZMM_LEN], zmm7);
		_mm_sfence();
		return;
	}
	else if (_nelems <= MEMMOVE_256ELEMS) {
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
		const __m512 zmm0(_mm512_loadu_ps((float*)&src[0]));
		const __m512 zmm1(_mm512_loadu_ps((float*)&src[1 * ZMM_LEN]));
		const __m512 zmm2(_mm512_loadu_ps((float*)&src[2 * ZMM_LEN]));
		const __m512 zmm3(_mm512_loadu_ps((float*)&src[3 * ZMM_LEN]));
		const __m512 zmm4(_mm512_loadu_ps((float*)&src[4 * ZMM_LEN]));
		const __m512 zmm5(_mm512_loadu_ps((float*)&src[5 * ZMM_LEN]));
		const __m512 zmm6(_mm512_loadu_ps((float*)&src[6 * ZMM_LEN]));
		const __m512 zmm7(_mm512_loadu_ps((float*)&src[7 * ZMM_LEN]));
		const __m512 zmm8(_mm512_loadu_ps((float*)&src[8 * ZMM_LEN]));
		const __m512 zmm9(_mm512_loadu_ps((float*)&src[9 * ZMM_LEN]));
		const __m512 zmm10(_mm512_loadu_ps((float*)&src[10 * ZMM_LEN]));
		const __m512 zmm11(_mm512_loadu_ps((float*)&src[11 * ZMM_LEN]));
		const __m512 zmm12(_mm512_loadu_ps((float*)&src[12 * ZMM_LEN]));
		const __m512 zmm13(_mm512_loadu_ps((float*)&src[13 * ZMM_LEN]));
		const __m512 zmm14(_mm512_loadu_ps((float*)&src[14 * ZMM_LEN]));
		const __m512 zmm15(_mm512_loadu_ps((float*)&src[15 * ZMM_LEN]));
		_mm512_stream_ps((float*)&dst[0], zmm0);
		_mm512_stream_ps((float*)&dst[1 * ZMM_LEN], zmm1);
		_mm512_stream_ps((float*)&dst[2 * ZMM_LEN], zmm2);
		_mm512_stream_ps((float*)&dst[3 * ZMM_LEN], zmm3);
		_mm512_stream_ps((float*)&dst[4 * ZMM_LEN], zmm4);
		_mm512_stream_ps((float*)&dst[5 * ZMM_LEN], zmm5);
		_mm512_stream_ps((float*)&dst[6 * ZMM_LEN], zmm6);
		_mm512_stream_ps((float*)&dst[7 * ZMM_LEN], zmm7);
		_mm512_stream_ps((float*)&dst[8 * ZMM_LEN], zmm8);
		_mm512_stream_ps((float*)&dst[9 * ZMM_LEN], zmm9);
		_mm512_stream_ps((float*)&dst[10 * ZMM_LEN], zmm10);
		_mm512_stream_ps((float*)&dst[11 * ZMM_LEN], zmm11);
		_mm512_stream_ps((float*)&dst[12 * ZMM_LEN], zmm12);
		_mm512_stream_ps((float*)&dst[13 * ZMM_LEN], zmm13);
		_mm512_stream_ps((float*)&dst[14 * ZMM_LEN], zmm14);
		_mm512_stream_ps((float*)&dst[15 * ZMM_LEN], zmm15);
		_mm_sfence();
		return;
	}
	else if (_nelems <= PAGE4KiB) {
		int32_t i;
		for (i = 0; i != ROUND_TO_SIXTEEN(_nelems, 16); i += 128) {
			_mm_prefetch((const char *)&src[i + 0], _MM_HINT_T0);
			const __m512 zmm0(_mm512_loadu_ps((float*)&src[i + 0]));
			const __m512 zmm1(_mm512_loadu_ps((float*)&src[i + 1 * ZMM_LEN]));
			const __m512 zmm2(_mm512_loadu_ps((float*)&src[i + 2 * ZMM_LEN]));
			const __m512 zmm3(_mm512_loadu_ps((float*)&src[i + 3 * ZMM_LEN]));
			const __m512 zmm4(_mm512_loadu_ps((float*)&src[i + 4 * ZMM_LEN]));
			const __m512 zmm5(_mm512_loadu_ps((float*)&src[i + 5 * ZMM_LEN]));
			const __m512 zmm6(_mm512_loadu_ps((float*)&src[i + 6 * ZMM_LEN]));
			const __m512 zmm7(_mm512_loadu_ps((float*)&src[i + 7 * ZMM_LEN]));
			_mm512_stream_ps((float*)&dst[i + 0], zmm0);
			_mm512_stream_ps((float*)&dst[i + 1 * ZMM_LEN], zmm1);
			_mm512_stream_ps((float*)&dst[i + 2 * ZMM_LEN], zmm2);
			_mm512_stream_ps((float*)&dst[i + 3 * ZMM_LEN], zmm3);
			_mm512_stream_ps((float*)&dst[i + 4 * ZMM_LEN], zmm4);
			_mm512_stream_ps((float*)&dst[i + 5 * ZMM_LEN], zmm5);
			_mm512_stream_ps((float*)&dst[i + 6 * ZMM_LEN], zmm6);
			_mm512_stream_ps((float*)&dst[i + 7 * ZMM_LEN], zmm7);
		}
		_mm_sfence();
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
				const __m512 zmm1(_mm512_loadu_ps((float*)&src[j + 1 * ZMM_LEN]));
				const __m512 zmm2(_mm512_loadu_ps((float*)&src[j + 2 * ZMM_LEN]));
				const __m512 zmm3(_mm512_loadu_ps((float*)&src[j + 3 * ZMM_LEN]));
				const __m512 zmm4(_mm512_loadu_ps((float*)&src[j + 4 * ZMM_LEN]));
				const __m512 zmm5(_mm512_loadu_ps((float*)&src[j + 5 * ZMM_LEN]));
				const __m512 zmm6(_mm512_loadu_ps((float*)&src[j + 6 * ZMM_LEN]));
				const __m512 zmm7(_mm512_loadu_ps((float*)&src[j + 7 * ZMM_LEN]));
				_mm512_stream_ps((float*)&dst[j + 0], zmm0);
				_mm512_stream_ps((float*)&dst[j + 1 * ZMM_LEN], zmm1);
				_mm512_stream_ps((float*)&dst[j + 2 * ZMM_LEN], zmm2);
				_mm512_stream_ps((float*)&dst[j + 3 * ZMM_LEN], zmm3);
				_mm512_stream_ps((float*)&dst[j + 4 * ZMM_LEN], zmm4);
				_mm512_stream_ps((float*)&dst[j + 5 * ZMM_LEN], zmm5);
				_mm512_stream_ps((float*)&dst[j + 6 * ZMM_LEN], zmm6);
				_mm512_stream_ps((float*)&dst[j + 7 * ZMM_LEN], zmm7);
			}

		}
		_mm_sfence();
		return;
	}


	
}

#endif 


