
#include <iomanip>
#include "GMS_complex_vec1D_zmm8r8.h"   
#include "GMS_complex_common_zmm8r8.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"

//
// Implementation
//

//
//	Parametrized macros
//
/*
#if (GMS_MAN_PREFETCH) == 1

#if !defined (AVX512_COMPLEX_PREFETCH_FROM_OBJ)
#define AVX512_COMPLEX_PREFETCH_FROM_OBJ(obj,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Re[(idx)+(off)]),(hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Im[(idx)+(off)]),(hint));
#endif

#if !defined (AVX512_COMPLEX_PREFETCH_FROM_PTR)
#define AVX512_COMPLEX_PREFETCH_FROM_PTR(ptr,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(ptr)[(idx)+(off)]),(hint));
#endif
#endif


/*#if !defined (AVX512_COMPLEX_CHECK_FATAL_ERROR)
#define AVX512_COMPLEX_CHECK_FATAL_ERROR(ptr1,ptr2,nsize,msg) \
	 do {													  \
		if ((NULL  == (ptr1) && (nsize) != 0ULL) ||        \
		     (NULL == (ptr2) && (nsize) != 0ULL) ) {      \
			    StackWalker sw{};						  \
			    sw.ShowCallstack();						  \
			    ABORT_ON_ERROR(msg,MALLOC_FAILED)		  \
	 }											          \
  } while (0);
#endif*/

/*

#if !defined (AVX512_COMPLEX_ADDITION)
#define AVX512_COMPLEX_ADDITION(out,v1,v2,idx,off) \
	(out) = _mm512_add_pd(_mm512_mul_pd(_mm512_load_pd(&(v1).m_Re[(idx)+(off)]), \
	_mm512_load_pd(&(v2).m_Re[(idx)+(off)])), _mm512_mul_pd(_mm512_load_pd(&(v1).m_Im[(idx)+(off)]), \
	_mm512_load_pd(&(v2).m_Im[(idx)+(off)])));
#endif

#if !defined (AVX512_COMPLEX_SUBTRACTION)
#define AVX512_COMPLEX_SUBTRACTION(out,v1,v2,idx,off) \
	(out) = _mm512_sub_pd(_mm512_mul_pd(_mm512_load_pd(&(v1).m_Im[(idx)+(off)]), \
	_mm512_load_pd(&(v2).m_Re[(idx)+(off)])), _mm512_mul_pd(_mm512_load_pd(&(v1).m_Re[(idx)+(off)]), \
	_mm512_load_pd(&(v2).m_Im[(idx)+(off)])));
#endif

	// Warning macro parameter v2 must be an exact copy
	// of parameter v1. This should done by calling class (CV1DZMM8r8)
	// Move Constructor.
#if !defined (AVX512_COMPLEX_MAGNITUDE)
#define AVX512_COMPLEX_MAGNITUDE(out,v1,v2,idx,off) \
	(out) = _mm512_sqrt_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_load_pd(&(v1).m_Re[(idx)+(off)]), \
	_mm512_load_pd(&(v2).m_Re[(idx)+(off)])), _mm512_mul_pd(_mm512_load_pd(&(v1).m_Im[(idx)+(off)]), \
	_mm512_load_pd(&(v2).m_Im[(idx)+(off)]))));
#endif
*/

gms::math
::CV1DZMM8r8::
CV1DZMM8r8() {
   m_Re = NULL;
   m_Im = NULL;
   m_nsize = 0;
}


gms::math
::CV1DZMM8r8::
CV1DZMM8r8(const int32_t nsize) {
	using namespace gms::common;
        m_nsize = nsize;
	m_Re = (double*)gms_mm_malloc(static_cast<size_t>(m_nsize),64ULL);
	m_Im = (double*)gms_mm_malloc(static_cast<size_t>(m_nsize),64ULL);
}

gms::math::
CV1DZMM8r8
::CV1DZMM8r8(const double * __restrict Re,
	     const double * __restrict Im,
	     const int32_t nsize) {
	using namespace gms::common;
        m_nsize = nsize;
	m_Re = (double*)gms_mm_malloc(static_cast<size_t>(m_nsize),64ULL);
	m_Im = (double*)gms_mm_malloc(static_cast<size_t>(m_nsize),64ULL);
#if (USE_NT_STORES) == 1
	avx512_uncached_memmove(&m_Re[0], &Re[0], m_nsize);
	avx512_uncached_memmove(&m_Im[0], &Im[0], m_nsize);
#else
	avx512_cached_memmove(&m_Re[0], &Re[0], m_nsize);
	avx512_cached_memmove(&m_Im[0], &Im[0], m_nsize);
#endif
}



gms::math::CV1DZMM8r8::
CV1DZMM8r8(const CV1DZMM8r8 &x) {
	using namespace gms::common;
        m_nsize = x.m_nsize;
	m_Re = (double*)gms_mm_malloc(static_cast<size_t>(m_nsize),64ULL);
	m_Im = (double*)gms_mm_malloc(static_cast<size_t>(m_nsize),64ULL);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
	avx512_uncached_memmove(&m_Re[0], &x.m_m_Re[0], m_nsize);
	avx512_uncached_memmove(&m_Im[0], &x.m_m_Im[0], m_nsize);
#else
	avx512_cached_memmove(&m_Re[0], &x.m_Re[0], m_nsize);
	avx512_cached_memmove(&m_Im[0], &x.m_Im[0], m_nsize);
#endif
}

gms::math
::CV1DZMM8r8::
CV1DZMM8r8(CV1DZMM8r8 &&x) {
	m_Re    = &x.m_Re[0];
	m_Im    = &x.m_Im[0];
	m_nsize = x.m_nsize;
	x.m_Re = NULL;
	x.m_Im = NULL;
	x.m_nsize = 0;
}

	


gms::math
::CV1DZMM8r8::
~CV1DZMM8r8() {

     gms_mm_free(m_Re); m_Re = NULL;
     gms_mm_free(m_Im); m_Im = NULL;

}	
	

/*gms::math::CV1DZMM8r8 &
gms::math::CV1DZMM8r8
::operator=(const CV1DZMM8r8 &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	if (m_nsize != x.m_nsize) {
#if defined _WIN64
  #if (GMS_DEBUG_ON) == 1
		_aligned_free_dbg(m_Re);
		_aligned_free_dbg(m_Im);
  #else
		_mm_free(m_Re);
		_mm_free(m_Im);
  #endif
#elif defined __linux
		_mm_free(m_Re);
		_mm_free(m_Im);
#endif
		
		m_nsize = 0; // Preserve an invariant here!!
		m_Re = NULL;
		m_Im = NULL;
#if defined _WIN64	
   #if (GMS_DEBUG_ON) == 1
		m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.m_nsize),align64B,__FILE__,__LINE__);
		m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.m_nsize),align64B,__FILE__,__LINE__);
   #else
		m_Re = gms_edmalloca(static_cast<size_t>(x.m_nsize), align64B);
		m_Im = gms_edmalloca(static_cast<size_t>(x.m_nsize), align64B);
   #endif
#elif defined __linux
   #if (GMS_DEBUG_ON) == 1
		m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.m_nsize),align64B);
		m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.m_nsize),align64B);
   #else
		m_Re = gms_edmalloca(static_cast<size_t>(x.m_nsize),align64B);
		m_Im = gms_edmalloca(static_cast<size_t>(x.m_nsize),align64B);
   #endif
#endif
   }
   else {
       
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		avx512_uncached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
		avx512_uncached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#else
		avx512_cached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
		avx512_cached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#endif
   }
	
	return (*this);
}*/

gms::math::CV1DZMM8r8 &
gms::math::CV1DZMM8r8
::operator=(CV1DZMM8r8 &&x) {
	if (this == &x) return (*this);

	gms_mm_free(m_Re);
	gms_mm_free(m_Im);
	m_Re      = x.m_Re;
	m_Im      = x.m_Im;
	m_nsize   = x.m_nsize;
	x.m_Re    = NULL;
	x.m_Im    = NULL;
	x.m_nsize = 0;
	return (*this);
}	


gms::math
::CV1DZMM8r8
gms::math::operator+(const CV1DZMM8r8 &x,
		     const CV1DZMM8r8 &y) {
	
	if (x.m_nsize != y.m_nsize) {
		return (CV1DZMM8r8{});
	}
	CV1DZMM8r8 ret_vec{x.m_nsize};
	int32_t i;
       
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16)  {
		 // Linearly growing indices, no need for software prefetch.


		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
		const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
#endif
		const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
		const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
#endif
		const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
		const __m512d zmm5 = _mm512_load_pd(&y.m_Im[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(zmm4, zmm5));
#endif
		const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i + 8]);
		const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(zmm6, zmm7));
#endif
		}	
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
		ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
	}
	return (ret_vec);
}		


	
	


gms::math::CV1DZMM8r8
gms::math::operator+(const CV1DZMM8r8 &x,
		     const double * __restrict Re) {
	
	CV1DZMM8r8 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {

		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
		const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
#else
		_mm512_storeu_pd(&ret_vec.m_Re[i + 0], _mm512_add_pd(zmm0, zmm1));
#endif
		const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
		const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_add_pd(zmm2, zmm3));
#endif
	}	
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] + Re[i];
	}
	return (ret_vec);
}	
	


gms::math
::CV1DZMM8r8
gms::math::operator-(const CV1DZMM8r8 &x,
		     const CV1DZMM8r8 &y) {
	if (x.m_nsize != y.m_nsize) {
		return (CV1DZMM8r8{});
	}
	CV1DZMM8r8 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {


		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
		const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
#endif
		const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
		const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
#endif
		const __m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
		const __m512d zmm5 = _mm512_load_pd(&y.m_Im[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_sub_pd(zmm4, zmm5));
#endif
		const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i + 8]);
		const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_sub_pd(zmm6, zmm7));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
		ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
	}
	return (ret_vec);
}		

	


gms::math::CV1DZMM8r8
gms::math::operator-(const CV1DZMM8r8 &x,
		     const double * __restrict Re) {
	
	CV1DZMM8r8 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {

	    
		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
		const __m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(zmm0, zmm1));
#endif
		const __m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
		const __m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(zmm2, zmm3));
#endif		

	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] - Re[i];
	}
	return (ret_vec);
}

gms::math::CV1DZMM8r8
gms::math::operator*(const CV1DZMM8r8 &x,
		     const CV1DZMM8r8 &y) {
	if (x.m_nsize != y.m_nsize) {
		return (CV1DZMM8r8{});	
	}
	CV1DZMM8r8 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {

	
		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i+0]);
		const __m512d zmm1 = _mm512_load_pd(&y.m_Re[i+0]);
		const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i+0]);
		const __m512d zmm3 = _mm512_load_pd(&y.m_Im[i+0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(
			_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2,zmm3)));
		_mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(
			_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_sub_pd(
			_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm3)));
		_mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_add_pd(
			_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm3)));

#endif
		const __m512d zmm4 = _mm512_load_pd(&x.m_Re[i+8]);
		const __m512d zmm5 = _mm512_load_pd(&y.m_Re[i+8]);
		const __m512d zmm6 = _mm512_load_pd(&x.m_Im[i+8]);
		const __m512d zmm7 = _mm512_load_pd(&y.m_Im[i+8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(
			_mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6,zmm7)));
		_mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(
			_mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_sub_pd(
			_mm512_mul_pd(zmm4, zmm5), _mm512_mul_pd(zmm6, zmm7)));
		_mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_add_pd(
			_mm512_mul_pd(zmm6, zmm5), _mm512_mul_pd(zmm4, zmm7)));
#endif		

	}

	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
		ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
	}
	return (ret_vec);
}

gms::math::CV1DZMM8r8
gms::math::operator*(const CV1DZMM8r8 &x,
		     const double * __restrict Re) {
	
	CV1DZMM8r8 ret_vec{x.m_nsize}; 
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {

		__m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
		__m512d zmm1 = _mm512_load_pd(&Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_mul_pd(zmm0, zmm1));
#endif
		__m512d zmm2 = _mm512_load_pd(&x.m_Re[i + 8]);
		__m512d zmm3 = _mm512_load_pd(&Re[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 8], _mm512_mul_pd(zmm2, zmm3));
#endif
		__m512d zmm4 = _mm512_load_pd(&x.m_Im[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_mul_pd(zmm4, zmm1));
#endif
		__m512d zmm5 = _mm512_load_pd(&x.m_Im[i + 8]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 8], _mm512_mul_pd(zmm5, zmm3));
#endif
		

	}

	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] * Re[i];
		ret_vec.m_Im[i] = x.m_Im[i] * Re[i];
	}
	return (ret_vec);
}

gms::math::CV1DZMM8r8
gms::math::operator/(const CV1DZMM8r8 &x,
		     const CV1DZMM8r8 &y) {
	if (x.m_nsize != y.m_nsize) { return (CV1DZMM8r8{}); }
	CV1DZMM8r8 ret_vec{x.m_nsize}; 
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 8) {
			// Will unrolling 2x not saturate divider unit.
			// We have two parallel division so at least second
			// operation will be pipelined at divider level.

		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i + 0]);
		const __m512d zmm1 = _mm512_load_pd(&y.m_Im[i + 0]);
		const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i + 0]);
		const __m512d re_term1 = _mm512_add_pd(
			_mm512_mul_pd(zmm0, zmm1), _mm512_mul_pd(zmm2, zmm1));
		const __m512d re_term2 = _mm512_add_pd(
			_mm512_mul_pd(zmm2, zmm1), _mm512_mul_pd(zmm0, zmm1));
		const __m512d zmm3 = _mm512_load_pd(&y.m_Re[i + 0]);
		const __m512d den_term = _mm512_add_pd(
			_mm512_mul_pd(zmm3, zmm3), _mm512_mul_pd(zmm1, zmm1));
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
#else
		_mm512_store_pd(&ret_vec.m_Re[i + 0], _mm512_div_pd(re_term1, den_term));
#endif
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));
#else
		_mm512_store_pd(&ret_vec.m_Im[i + 0], _mm512_div_pd(re_term2, den_term));
#endif

	}
	
	for (; i != ret_vec.m_nsize; ++i) {
		const double tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
		const double tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
		const double den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
		ret_vec.m_Re[i] = tre / den;
		ret_vec.m_Im[i] = tim / den;
	}
	return (ret_vec);
}

gms::math::CV1DZMM8r8
gms::math::operator/(const  CV1DZMM8r8 &x,
		     const double * __restrict Re) {
	
        CV1DZMM8r8 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 8) {
		// Will unrolling 2x not saturate divider unit.
		// We have two parallel division so at least second
		// operation will be pipelined at divider level.

		const __m512d zmm0 = _mm512_load_pd(&x.m_Re[i]);
		const __m512d zmm1 = _mm512_load_pd(&Re[i]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Re[i], _mm512_div_pd(zmm0, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Re[i], _mm512_div_pd(zmm0, zmm1));
#endif
		const __m512d zmm2 = _mm512_load_pd(&x.m_Im[i]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_pd(&ret_vec.m_Im[i], _mm512_div_pd(zmm2, zmm1));
#else
		_mm512_store_pd(&ret_vec.m_Im[i], _mm512_div_pd(zmm2, zmm1));
#endif
	}

	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] / Re[i];
		ret_vec.m_Im[i] = x.m_Im[i] / Re[i];
	}
	return (ret_vec);
}






	
std::ostream &
gms::math::operator<<(std::ostream &os,
		      const CV1DZMM8r8 &x) {

	for (int64_t i = 0LL; i != x.m_nsize; ++i) {
		os << std::fixed << std::showpoint << std::setprecision(15) <<
			std::setw(4) << "Re: " << "{" << x.m_Re[i] << "}" <<
			std::setw(12) << "Im: " << "{" << x.m_Im[i] << "}" << std::endl;
	}
	return (os);
}



	
	




// vout size must be equal to v1 and v2 size.
void
gms::math::v512cnormalize_product(CV1DZMM8r8 &vout, 
				  const CV1DZMM8r8 &v1,
				  const CV1DZMM8r8 &v2,
				  const bool do_nt_stream) {
	
	avx512_cnormalize_prod<CV1DZMM8r8>(vout,v1,v2,do_nt_stream);
}	




void
gms::math::v512cmean_product(std::complex<double> &mean, 
			     const CV1DZMM8r8 &v1,
			     const CV1DZMM8r8 &v2) {
	
	avx512_cmean_prod<CV1DZMM8r8>(mean,v1,v2);
}
	


void
gms::math::v512cmean_quotient(std::complex<double> &mean,
			      const CV1DZMM8r8 &v1,
			      const CV1DZMM8r8 &v2) {
	
	avx512_cmean_quot<CV1DZMM8r8>(mean,v1,v2);
}

void
gms::math::v512cconj_product(CV1DZMM8r8    &vout,
			     const CV1DZMM8r8 &v1,
			     const CV1DZMM8r8 &v2,
			     const bool do_nt_store) {
	
	avx512_cconj_prod<CV1DZMM8r8>(vout,v1,v2,do_nt_store);
}      


void
gms::math::v512cnorm_conjprod(CV1DZMM8r8    &vout,
			      const CV1DZMM8r8 &v1,
			      const CV1DZMM8r8 &v2,
			      const bool do_nt_store) {
	
	avx512_cnorm_conjprod<CV1DZMM8r8>(vout,v1,v2,do_nt_store);
}


void
gms::math::v512cmean_conjprod(std::complex<double> &mean,
			      const CV1DZMM8r8 &v1,
			      const CV1DZMM8r8 &v2) {
	
	avx512_cmean_conjprod<CV1DZMM8r8>(mean,v1,v2);
}

void
gms::math::v512c_arithmean(std::complex<double> &mean,
			   const CV1DZMM8r8 &v1) {
	
	avx512_arith_mean<CV1DZMM8r8>(mean,v1);
	
}

void
gms::math::v512c_normalize(CV1DZMM8r8 &vnorm,
			   const CV1DZMM8r8 &v,
			   const CV1DZMM8r8 &cv,
			   const bool do_nt_store) {
	
	avx512_cnormalize<CV1DZMM8r8>(vnorm,v,cv,do_nt_store);
}

void
gms::math::v512c_magnitude(double * __restrict vmag,
			   const CV1DZMM8r8 &v,
			   const CV1DZMM8r8 &cv) {
	
	avx512_cmagnitude<CV1DZMM8r8>(vmag,v,cv);
}	







