

#include "GMS_complex_vec1D_zmm16r4.h"   
#include "GMS_complex_common_zmm16r4.h"
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
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).data.m_Re[(idx)+(off)]),(hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).data.m_Im[(idx)+(off)]),(hint));
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
	(out) = _mm512_add_ps(_mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Re[(idx)+(off)])), _mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Im[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Im[(idx)+(off)])));
#endif

#if !defined (AVX512_COMPLEX_SUBTRACTION)
#define AVX512_COMPLEX_SUBTRACTION(out,v1,v2,idx,off) \
	(out) = _mm512_sub_ps(_mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Im[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Re[(idx)+(off)])), _mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Im[(idx)+(off)])));
#endif

	// Warning macro parameter v2 must be an exact copy
	// of parameter v1. This should done by calling class (CV1DZMM16r4)
	// Move Constructor.
#if !defined (AVX512_COMPLEX_MAGNITUDE)
#define AVX512_COMPLEX_MAGNITUDE(out,v1,v2,idx,off) \
	(out) = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Re[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Re[(idx)+(off)])), _mm512_mul_ps(_mm512_load_ps(&(v1).data.m_Im[(idx)+(off)]), \
	_mm512_load_ps(&(v2).data.m_Im[(idx)+(off)]))));
#endif
*/

gms::math
::CV1DZMM16r4::
CV1DZMM16r4() {
   data.m_Re = NULL;
   data.m_Im = NULL;
   data.m_nsize = 0;
}


gms::math
::CV1DZMM16r4::
CV1DZMM16r4(const int32_t nsize) {
	using namespace gms::common;
        data.m_nsize = nsize;
	data.m_Re = (float*)gms_mm_malloc(static_cast<size_t>(data.m_nsize),64ULL);
	data.m_Im = (float*)gms_mm_malloc(static_cast<size_t>(data.m_nsize),64ULL);
}

gms::math::
CV1DZMM16r4
::CV1DZMM16r4(const float * __restrict Re,
	     const float * __restrict Im,
	     const int32_t nsize) {
	using namespace gms::common;
        data.m_nsize = nsize;
	data.m_Re = (float*)gms_mm_malloc(static_cast<size_t>(data.m_nsize),64ULL);
	data.m_Im = (float*)gms_mm_malloc(static_cast<size_t>(data.m_nsize),64ULL);
#if (USE_NT_STORES) == 1
	avx512_uncached_memmove(&data.m_Re[0], &Re[0], data.m_nsize);
	avx512_uncached_memmove(&data.m_Im[0], &Im[0], data.m_nsize);
#else
	avx512_cached_memmove(&data.m_Re[0], &Re[0], data.m_nsize);
	avx512_cached_memmove(&data.m_Im[0], &Im[0], data.m_nsize);
#endif
}



gms::math::CV1DZMM16r4::
CV1DZMM16r4(const CV1DZMM16r4 &x) {
	using namespace gms::common;
        data.m_nsize = x.data.m_nsize;
	data.m_Re = (float*)gms_mm_malloc(static_cast<size_t>(data.m_nsize),64ULL);
	data.m_Im = (float*)gms_mm_malloc(static_cast<size_t>(data.m_nsize),64ULL);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
	avx512_uncached_memmove(&data.m_Re[0], &x.m_data.m_Re[0], data.m_nsize);
	avx512_uncached_memmove(&data.m_Im[0], &x.m_data.m_Im[0], data.m_nsize);
#else
	avx512_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], data.m_nsize);
	avx512_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], data.m_nsize);
#endif
}

gms::math
::CV1DZMM16r4::
CV1DZMM16r4(CV1DZMM16r4 &&x) {
	data.m_Re    = &x.data.m_Re[0];
	data.m_Im    = &x.data.m_Im[0];
	data.m_nsize = x.data.m_nsize;
	x.data.m_Re = NULL;
	x.data.m_Im = NULL;
	x.data.m_nsize = 0;
}

	


gms::math
::CV1DZMM16r4::
~CV1DZMM16r4() {

     gms_mm_free(data.m_Re); data.m_Re = NULL;
     gms_mm_free(data.m_Im); data.m_Im = NULL;

}	
	

/*gms::math::CV1DZMM16r4 &
gms::math::CV1DZMM16r4
::operator=(const CV1DZMM16r4 &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	if (data.m_nsize != x.data.m_nsize) {
#if defined _WIN64
  #if (GMS_DEBUG_ON) == 1
		_aligned_free_dbg(data.m_Re);
		_aligned_free_dbg(data.m_Im);
  #else
		_mm_free(data.m_Re);
		_mm_free(data.m_Im);
  #endif
#elif defined __linux
		_mm_free(data.m_Re);
		_mm_free(data.m_Im);
#endif
		
		data.m_nsize = 0; // Preserve an invariant here!!
		data.m_Re = NULL;
		data.m_Im = NULL;
#if defined _WIN64	
   #if (GMS_DEBUG_ON) == 1
		data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize),align64B,__FILE__,__LINE__);
		data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize),align64B,__FILE__,__LINE__);
   #else
		data.m_Re = gms_edmalloca(static_cast<size_t>(x.data.m_nsize), align64B);
		data.m_Im = gms_edmalloca(static_cast<size_t>(x.data.m_nsize), align64B);
   #endif
#elif defined __linux
   #if (GMS_DEBUG_ON) == 1
		data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize),align64B);
		data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize),align64B);
   #else
		data.m_Re = gms_edmalloca(static_cast<size_t>(x.data.m_nsize),align64B);
		data.m_Im = gms_edmalloca(static_cast<size_t>(x.data.m_nsize),align64B);
   #endif
#endif
   }
   else {
       
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		avx512_uncached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
		avx512_uncached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#else
		avx512_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
		avx512_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#endif
   }
	
	return (*this);
}*/

gms::math::CV1DZMM16r4 &
gms::math::CV1DZMM16r4
::operator=(CV1DZMM16r4 &&x) {
	if (this == &x) return (*this);

	gms_mm_free(data.m_Re);
	gms_mm_free(data.m_Im);
	data.m_Re      = x.data.m_Re;
	data.m_Im      = x.data.m_Im;
	data.m_nsize   = x.data.m_nsize;
	x.data.m_Re    = NULL;
	x.data.m_Im    = NULL;
	x.data.m_nsize = 0;
	return (*this);
}	


gms::math
::CV1DZMM16r4
gms::math::operator+(const CV1DZMM16r4 &x,
		     const CV1DZMM16r4 &y) {
	
	if (x.data.m_nsize != y.data.m_nsize) {
		return (CV1DZMM16r4{});
	}
	CV1DZMM16r4 ret_vec{x.data.m_nsize};
	int32_t i;
       
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize, 16); i += 32)  {
		 // Linearly growing indices, no need for software prefetch.


		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i + 0]);
		const __m512 zmm1 = _mm512_load_ps(&y.data.m_Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_add_ps(zmm0, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 0], _mm512_add_ps(zmm0, zmm1));
#endif
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Re[i + 16]);
		const __m512 zmm3 = _mm512_load_ps(&y.data.m_Re[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 16], _mm512_add_ps(zmm2, zmm3));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 16], _mm512_add_ps(zmm2, zmm3));
#endif
		const __m512 zmm4 = _mm512_load_ps(&x.data.m_Im[i + 0]);
		const __m512 zmm5 = _mm512_load_ps(&y.data.m_Im[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 0], _mm512_add_ps(zmm4, zmm5));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 0], _mm512_add_ps(zmm4, zmm5));
#endif
		const __m512 zmm6 = _mm512_load_ps(&x.data.m_Im[i + 16]);
		const __m512 zmm7 = _mm512_load_ps(&y.data.m_Im[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 16], _mm512_add_ps(zmm6, zmm7));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 16], _mm512_add_ps(zmm6, zmm7));
#endif
		}	
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] + y.data.m_Re[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] + y.data.m_Re[i];
	}
	return (ret_vec);
}		


	
	


gms::math::CV1DZMM16r4
gms::math::operator+(const CV1DZMM16r4 &x,
		     const float * __restrict Re) {
	
	CV1DZMM16r4 ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize,16); i += 32) {

		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i + 0]);
		const __m512 zmm1 = _mm512_load_ps(&Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_add_ps(zmm0, zmm1));
#else
		_mm512_storeu_ps(&ret_vec.data.m_Re[i + 0], _mm512_add_ps(zmm0, zmm1));
#endif
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Re[i + 16]);
		const __m512 zmm3 = _mm512_load_ps(&Re[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 16], _mm512_add_ps(zmm2, zmm3));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 16], _mm512_add_ps(zmm2, zmm3));
#endif
	}	
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] + Re[i];
	}
	return (ret_vec);
}	
	


gms::math
::CV1DZMM16r4
gms::math::operator-(const CV1DZMM16r4 &x,
		     const CV1DZMM16r4 &y) {
	if (x.data.m_nsize != y.data.m_nsize) {
		return (CV1DZMM16r4{});
	}
	CV1DZMM16r4 ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize,16); i += 32) {


		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i + 0]);
		const __m512 zmm1 = _mm512_load_ps(&y.data.m_Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_sub_ps(zmm0, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 0], _mm512_sub_ps(zmm0, zmm1));
#endif
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Re[i + 16]);
		const __m512 zmm3 = _mm512_load_ps(&y.data.m_Re[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 16], _mm512_sub_ps(zmm2, zmm3));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 16], _mm512_sub_ps(zmm2, zmm3));
#endif
		const __m512 zmm4 = _mm512_load_ps(&x.data.m_Im[i + 0]);
		const __m512 zmm5 = _mm512_load_ps(&y.data.m_Im[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 0], _mm512_sub_ps(zmm4, zmm5));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 0], _mm512_sub_ps(zmm4, zmm5));
#endif
		const __m512 zmm6 = _mm512_load_ps(&x.data.m_Im[i + 16]);
		const __m512 zmm7 = _mm512_load_ps(&y.data.m_Im[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 16], _mm512_sub_ps(zmm6, zmm7));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 16], _mm512_sub_ps(zmm6, zmm7));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] - y.data.m_Re[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] - y.data.m_Im[i];
	}
	return (ret_vec);
}		

	


gms::math::CV1DZMM16r4
gms::math::operator-(const CV1DZMM16r4 &x,
		     const float * __restrict Re) {
	
	CV1DZMM16r4 ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize,16); i += 32) {

	    
		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i + 0]);
		const __m512 zmm1 = _mm512_load_ps(&Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_sub_ps(zmm0, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 0], _mm512_sub_ps(zmm0, zmm1));
#endif
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Re[i + 16]);
		const __m512 zmm3 = _mm512_load_ps(&Re[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 16], _mm512_sub_ps(zmm2, zmm3));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 16], _mm512_sub_ps(zmm2, zmm3));
#endif		

	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] - Re[i];
	}
	return (ret_vec);
}

gms::math::CV1DZMM16r4
gms::math::operator*(const CV1DZMM16r4 &x,
		     const CV1DZMM16r4 &y) {
	if (x.data.m_nsize != y.data.m_nsize) {
		return (CV1DZMM16r4{});	
	}
	CV1DZMM16r4 ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize,16); i += 32) {

	
		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i+0]);
		const __m512 zmm1 = _mm512_load_ps(&y.data.m_Re[i+0]);
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Im[i+0]);
		const __m512 zmm3 = _mm512_load_ps(&y.data.m_Im[i+0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_sub_ps(
			_mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2,zmm3)));
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 0], _mm512_add_ps(
			_mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0, zmm3)));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 0], _mm512_sub_ps(
			_mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2, zmm3)));
		_mm512_store_ps(&ret_vec.data.m_Im[i + 0], _mm512_add_ps(
			_mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0, zmm3)));

#endif
		const __m512 zmm4 = _mm512_load_ps(&x.data.m_Re[i+16]);
		const __m512 zmm5 = _mm512_load_ps(&y.data.m_Re[i+16]);
		const __m512 zmm6 = _mm512_load_ps(&x.data.m_Im[i+16]);
		const __m512 zmm7 = _mm512_load_ps(&y.data.m_Im[i+16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 16], _mm512_sub_ps(
			_mm512_mul_ps(zmm4, zmm5), _mm512_mul_ps(zmm6,zmm7)));
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 16], _mm512_add_ps(
			_mm512_mul_ps(zmm6, zmm5), _mm512_mul_ps(zmm4, zmm7)));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 16], _mm512_sub_ps(
			_mm512_mul_ps(zmm4, zmm5), _mm512_mul_ps(zmm6, zmm7)));
		_mm512_store_ps(&ret_vec.data.m_Im[i + 16], _mm512_add_ps(
			_mm512_mul_ps(zmm6, zmm5), _mm512_mul_ps(zmm4, zmm7)));
#endif		

	}

	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = (x.data.m_Re[i] * y.data.m_Re[i]) - (x.data.m_Im[i] * y.data.m_Im[i]);
		ret_vec.data.m_Im[i] = (x.data.m_Im[i] * y.data.m_Re[i]) + (x.data.m_Re[i] * y.data.m_Im[i]);
	}
	return (ret_vec);
}

gms::math::CV1DZMM16r4
gms::math::operator*(const CV1DZMM16r4 &x,
		     const float * __restrict Re) {
	
	CV1DZMM16r4 ret_vec{x.data.m_nsize}; 
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize,16); i += 32) {

		__m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i + 0]);
		__m512 zmm1 = _mm512_load_ps(&Re[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_mul_ps(zmm0, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 0], _mm512_mul_ps(zmm0, zmm1));
#endif
		__m512 zmm2 = _mm512_load_ps(&x.data.m_Re[i + 16]);
		__m512 zmm3 = _mm512_load_ps(&Re[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 16], _mm512_mul_ps(zmm2, zmm3));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 16], _mm512_mul_ps(zmm2, zmm3));
#endif
		__m512 zmm4 = _mm512_load_ps(&x.data.m_Im[i + 0]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 0], _mm512_mul_ps(zmm4, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 0], _mm512_mul_ps(zmm4, zmm1));
#endif
		__m512 zmm5 = _mm512_load_ps(&x.data.m_Im[i + 16]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 16], _mm512_mul_ps(zmm5, zmm3));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 16], _mm512_mul_ps(zmm5, zmm3));
#endif
		

	}

	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] * Re[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] * Re[i];
	}
	return (ret_vec);
}

gms::math::CV1DZMM16r4
gms::math::operator/(const CV1DZMM16r4 &x,
		     const CV1DZMM16r4 &y) {
	if (x.data.m_nsize != y.data.m_nsize) { return (CV1DZMM16r4{}); }
	CV1DZMM16r4 ret_vec{x.data.m_nsize}; 
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize,16); i += 16) {
			// Will unrolling 2x not saturate divider unit.
			// We have two parallel division so at least second
			// operation will be pipelined at divider level.

		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i + 0]);
		const __m512 zmm1 = _mm512_load_ps(&y.data.m_Im[i + 0]);
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Im[i + 0]);
		const __m512 re_term1 = _mm512_add_ps(
			_mm512_mul_ps(zmm0, zmm1), _mm512_mul_ps(zmm2, zmm1));
		const __m512 re_term2 = _mm512_add_ps(
			_mm512_mul_ps(zmm2, zmm1), _mm512_mul_ps(zmm0, zmm1));
		const __m512 zmm3 = _mm512_load_ps(&y.data.m_Re[i + 0]);
		const __m512 den_term = _mm512_add_ps(
			_mm512_mul_ps(zmm3, zmm3), _mm512_mul_ps(zmm1, zmm1));
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i + 0], _mm512_div_ps(re_term1, den_term));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i + 0], _mm512_div_ps(re_term1, den_term));
#endif
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i + 0], _mm512_div_ps(re_term2, den_term));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i + 0], _mm512_div_ps(re_term2, den_term));
#endif

	}
	
	for (; i != ret_vec.data.m_nsize; ++i) {
		const float tre = (x.data.m_Re[i] * y.data.m_Im[i]) + (x.data.m_Im[i] * y.data.m_Im[i]);
		const float tim = (x.data.m_Im[i] * y.data.m_Im[i]) - (x.data.m_Re[i] * y.data.m_Im[i]);
		const float den = (y.data.m_Re[i] * y.data.m_Re[i]) + (y.data.m_Im[i] * y.data.m_Im[i]);
		ret_vec.data.m_Re[i] = tre / den;
		ret_vec.data.m_Im[i] = tim / den;
	}
	return (ret_vec);
}

gms::math::CV1DZMM16r4
gms::math::operator/(const  CV1DZMM16r4 &x,
		     const float * __restrict Re) {
	
        CV1DZMM16r4 ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_SIXTEEN(ret_vec.data.m_nsize, 16); i += 16) {
		// Will unrolling 2x not saturate divider unit.
		// We have two parallel division so at least second
		// operation will be pipelined at divider level.

		const __m512 zmm0 = _mm512_load_ps(&x.data.m_Re[i]);
		const __m512 zmm1 = _mm512_load_ps(&Re[i]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Re[i], _mm512_div_ps(zmm0, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Re[i], _mm512_div_ps(zmm0, zmm1));
#endif
		const __m512 zmm2 = _mm512_load_ps(&x.data.m_Im[i]);
#if (USE_AVX512COMPLEX_NT_STORES) == 1
		_mm512_stream_ps(&ret_vec.data.m_Im[i], _mm512_div_ps(zmm2, zmm1));
#else
		_mm512_store_ps(&ret_vec.data.m_Im[i], _mm512_div_ps(zmm2, zmm1));
#endif
	}

	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] / Re[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] / Re[i];
	}
	return (ret_vec);
}






	
std::ostream &
gms::math::operator<<(std::ostream &os,
		      const CV1DZMM16r4 &x) {

	for (int64_t i = 0LL; i != x.data.m_nsize; ++i) {
		os << std::fixed << std::showpoint << std::setprecision(15) <<
			std::setw(4) << "Re: " << "{" << x.data.m_Re[i] << "}" <<
			std::setw(12) << "Im: " << "{" << x.data.m_Im[i] << "}" << std::endl;
	}
	return (os);
}



	
	




// vout size must be equal to v1 and v2 size.
void
gms::math::v512cnormalize_product(CV1DZMM16r4 &vout, 
				  const CV1DZMM16r4 &v1,
				  const CV1DZMM16r4 &v2,
				  const bool do_nt_stream) {
	
	avx512_cnormalize_prod<CV1DZMM16r4>(vout,v1,v2,do_nt_stream);
}	




void
gms::math::v512cmean_product(std::complex<float> &mean, 
			     const CV1DZMM16r4 &v1,
			     const CV1DZMM16r4 &v2) {
	
	avx512_cmean_prod<CV1DZMM16r4>(mean,v1,v2);
}
	


void
gms::math::v512cmean_quotient(std::complex<float> &mean,
			      const CV1DZMM16r4 &v1,
			      const CV1DZMM16r4 &v2) {
	
	avx512_cmean_quot<CV1DZMM16r4>(mean,v1,v2);
}

void
gms::math::v512cconj_product(CV1DZMM16r4    &vout,
			     const CV1DZMM16r4 &v1,
			     const CV1DZMM16r4 &v2,
			     const bool do_nt_store) {
	
	avx512_cconj_prod<CV1DZMM16r4>(vout,v1,v2,do_nt_store);
}      


void
gms::math::v512cnorm_conjprod(CV1DZMM16r4    &vout,
			      const CV1DZMM16r4 &v1,
			      const CV1DZMM16r4 &v2,
			      const bool do_nt_store) {
	
	avx512_cnorm_conjprod<CV1DZMM16r4>(vout,v1,v2,do_nt_store);
}


void
gms::math::v512cmean_conjprod(std::complex<float> &mean,
			      const CV1DZMM16r4 &v1,
			      const CV1DZMM16r4 &v2) {
	
	avx512_cmean_conjprod<CV1DZMM16r4>(mean,v1,v2);
}

void
gms::math::v512c_arithmean(std::complex<float> &mean,
			   const CV1DZMM16r4 &v1) {
	
	avx512_arith_mean<CV1DZMM16r4>(mean,v1);
	
}

void
gms::math::v512c_normalize(CV1DZMM16r4 &vnorm,
			   const CV1DZMM16r4 &v,
			   const CV1DZMM16r4 &cv,
			   const bool do_nt_store) {
	
	avx512_cnormalize<CV1DZMM16r4>(vnorm,v,cv,do_nt_store);
}

void
gms::math::v512c_magnitude(float * __restrict vmag,
			   const CV1DZMM16r4 &v,
			   const CV1DZMM16r4 &cv) {
	
	avx512_cmagnitude<CV1DZMM16r4>(vmag,v,cv);
}	







