
#include <iomanip>
#include <immintrin.h>
#include "GMS_complex_vec1D_ymm8r4.h"   
#include "GMS_complex_common_ymm8r4.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"
#include "GMS_common.h"
   

//
//	Implementation
//

/*#if (GMS_MAN_PREFETCH) == 1

#if !defined (AVXCOMPLEX_PREFETCH_FROM_OBJ)
#define AVXCOMPLEX_PREFETCH_FROM_OBJ(obj,idx,off,hint) \
    _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Re[(idx)+(off)],(hint)); \
    _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Im[(idx)+(off)],(hint));
#endif

#if !defined (AVXCOMPLEX_PREFETCH_FROM_PTR)
#define AVXCOMPLEX_PREFETCH_FROM_PTR(ptr,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(ptr)[(idx)+(off)],(hint));
#endif

#if !defined (AVCOMPLEX_PREFETCH_FROM_THIS)
#define AVXCOMPLEX_PREFETCH_FROM_THIS(obj,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj)->m_Re[(idx)+(off)],(hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj)->m_Im[(idx)+(off)],(hint));
#endif

#endif

/*#if !defined (AVXCOMPLEX_CHECK_FATAL_ERROR)
#define AVXCOMPLEX_CHECK_FATAL_ERROR(ptr1,ptr2,nsize,msg) \
  do {													  \
		if ((NULL  == (ptr1) && (nsize) != 0ULL) ||        \
		     (NULL == (ptr2) && (nsize) != 0ULL) ) {      \
			    StackWalker sw{};						  \
			    sw.ShowCallstack();						  \
			    ABORT_ON_ERROR(msg,MALLOC_FAILED)		  \
	 }											          \
  } while (0);
#endif*/



gms::math::CV1DYMM8r4::
CV1DYMM8r4(){
	m_Re = NULL;
	m_Im = NULL;
	m_nsize = 0;
}



gms::math::CV1DYMM8r4::
CV1DYMM8r4(const int32_t nsize){
	using namespace gms::common;
        m_nsize = nsize;
	m_Re = (float*)gms_mm_malloc(static_cast<size_t>(m_nsize),32ULL);
	m_Im = (float*)gms_mm_malloc(static_cast<size_t>(m_nsize),32ULL);
}	
	
	
gms::math::CV1DYMM8r4
::CV1DYMM8r4(   const float * __restrict Re,
		const float * __restrict Im,
		const int32_t nsize) {
	using namespace gms::common;
        m_nsize = nsize;
        m_Re = (float*)gms_mm_malloc(static_cast<size_t>(m_nsize),32ULL);
	m_Im = (float*)gms_mm_malloc(static_cast<size_t>(m_nsize),32ULL);
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
	avx256_uncached_memmove(&m_Re[0], &Re[0],m_nsize);
	avx256_uncached_memmove(&m_Im[0], &Im[0],m_nsize);
#else
	avx256_cached_memmove(&m_Re[0], &Re[0], m_nsize);
	avx256_cached_memmove(&m_Im[0], &Im[0], m_nsize);
#endif
}


	




gms::math::CV1DYMM8r4::
CV1DYMM8r4(const CV1DYMM8r4 &x){
       using namespace gms::common;
       m_nsize = x.m_nsize;
       m_Re = (float*)gms_mm_malloc(static_cast<size_t>(m_nsize),32ULL);
       m_Im = (float*)gms_mm_malloc(static_cast<size_t>(m_nsize),32ULL);
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
	avx256_uncached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
	avx256_uncached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#else
	avx256_cached_memmove(&m_Re[0], &x.m_Re[0], x.m_nsize);
	avx256_cached_memmove(&m_Im[0], &x.m_Im[0], x.m_nsize);
#endif
}

gms::math::CV1DYMM8r4::
CV1DYMM8r4(CV1DYMM8r4 &&x) {
	m_Re = &x.m_Re[0];
	m_Im = &x.m_Im[0];
	m_nsize = x.m_nsize;
	x.m_Re = NULL;
	x.m_Im = NULL;
	x.m_nsize = 0;
}


gms::math::CV1DYMM8r4::
		  ~CV1DYMM8r4() {
         using namespace gms::common;
	 gms_mm_free(m_Re); m_Re = NULL;
	 gms_mm_free(m_Im); m_Im = NULL;
}
		
	
	

/*gms::math::CV1DYMM8r4 &
gms::math::CV1DYMM8r4
::operator=(const CV1DYMM8r4 &x) {
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
	m_Re = NULL;
	m_Im = NULL;
	m_nsize = 0;
#if defined _WIN64
    #if (GMS_DEBUG_ON) == 1
	  m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.m_nsize), align64B, __FILE__, __LINE__);
	  m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.m_nsize), align64B, __FILE__, __LINE__);
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
	
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
	 avx256_uncached_memmove(&m_Re[0],&x.m_Re[0],x.m_nsize);
	 avx256_uncached_memmove(&m_Im[0],&x.m_Im[0],x.m_nsize);
#else
	 avx256_cached_memmove(&m_Re[0], &x.m_Re[0],x.m_nsize);
	 avx256_cached_memmove(&m_Im[0], &x.m_Im[0],x.m_nsize);
#endif		
	}
	return (*this);
}*/

gms::math::CV1DYMM8r4 &
gms::math::CV1DYMM8r4::operator=(CV1DYMM8r4 &&x) {
        using namespace gms::common;
	if (this == &x) return (*this);
	gms_mm_free(m_Re);
	gms_mm_free(m_Im);

	m_Re    = &x.m_Re[0];
	m_Im    = &x.m_Im[0];
	m_nsize = x.m_nsize;
	x.m_Re = NULL;
	x.m_Im = NULL;
	x.m_nsize = 0;
	return (*this);
}	

std::ostream &
gms::math::operator<<(std::ostream &os,
		      const CV1DYMM8r4 &x) {
	for (int64_t i = 0LL; i != x.m_nsize; ++i) {
		os << std::fixed << std::showpoint << std::setprecision(15) <<
			  std::setw(4) <<  "Re: " << "{" << x.m_Re[i] << "}" <<
			  std::setw(12) << "Im: " << "{" << x.m_Im[i] << "}" << std::endl;
	}
	return (os);
}
	
gms::math::CV1DYMM8r4
gms::math::operator+(const CV1DYMM8r4 &x,
		     const CV1DYMM8r4 &y) {
	if (x.m_nsize != y.m_nsize) {
		return (CV1DYMM8r4{});
	}
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {
		// Linearly growing indices, no need for software prefetch.
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i+0]));
		const __m256 ymm1(_mm256_load_ps(&y.m_Re[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+0], _mm256_add_ps(ymm0,ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Re[i + 0], _mm256_add_ps(ymm0, ymm1));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Re[i+8]));
		const __m256 ymm3(_mm256_load_ps(&y.m_Re[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+8], _mm256_add_ps(ymm2,ymm3));
#else
		_mm256_store_ps(&ret_vec.m_Re[i + 8], _mm256_add_ps(ymm2, ymm3));
#endif
		const __m256 ymm4(_mm256_load_ps(&x.m_Im[i+0]));
		const __m256 ymm5(_mm256_load_ps(&y.m_Im[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i+0], _mm256_add_ps(ymm4,ymm5));
#else
		_mm256_store_ps(&ret_vec.m_Im[i+0], _mm256_add_ps(ymm4, ymm5));
#endif
		const __m256 ymm6(_mm256_load_ps(&x.m_Im[i+8]));
		const __m256 ymm7(_mm256_load_ps(&y.m_Im[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i+8], _mm256_add_ps(ymm6,ymm7));
#else
		_mm256_store_ps(&ret_vec.m_Im[i+8], _mm256_add_ps(ymm6, ymm7));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] + y.m_Re[i];
		ret_vec.m_Im[i] = x.m_Im[i] + y.m_Re[i];
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator+(const CV1DYMM8r4 &x,
		     const float * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) {
		return (CV1DYMM8r4{});
	}
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i+0]));
		const __m256 ymm1(_mm256_load_ps(&y[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+0],_mm256_add_ps(ymm0,ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+0], _mm256_add_ps(ymm0, ymm1));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Re[i+8]));
		const __m256 ymm3(_mm256_load_ps(&y[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+8],_mm256_add_ps(ymm2,ymm3));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+8], _mm256_add_ps(ymm2, ymm3));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] + y[i];
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator-(const CV1DYMM8r4 &x,
		     const CV1DYMM8r4 &y) {
	if (x.m_nsize != y.m_nsize) { return (CV1DYMM8r4{}); }
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i+0]));
		const __m256 ymm1(_mm256_load_ps(&y.m_Re[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+0], _mm256_sub_ps(ymm0,ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+0], _mm256_sub_ps(ymm0, ymm1));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Re[i+8]));
		const __m256 ymm3(_mm256_load_ps(&y.m_Re[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+8], _mm256_sub_ps(ymm2,ymm3));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+8], _mm256_sub_ps(ymm2, ymm3));
#endif
		const __m256 ymm4(_mm256_load_ps(&x.m_Im[i+0]));
		const __m256 ymm5(_mm256_load_ps(&y.m_Im[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i+0], _mm256_sub_ps(ymm4,ymm5));
#else
		_mm256_store_ps(&ret_vec.m_Im[i+0], _mm256_sub_ps(ymm4, ymm5));
#endif
		const __m256 ymm6(_mm256_load_ps(&x.m_Im[i+8]));
		const __m256 ymm7(_mm256_load_ps(&y.m_Im[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i+8], _mm256_sub_ps(ymm6,ymm7));
#else
		_mm256_store_ps(&ret_vec.m_Im[i+8], _mm256_sub_ps(ymm6, ymm7));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] - y.m_Re[i];
		ret_vec.m_Im[i] = x.m_Im[i] - y.m_Im[i];
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator-(const CV1DYMM8r4 &x,
		     const float * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) {return (CV1DYMM8r4{});}
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize, 8); i += 16) {
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i+0]));
		const __m256 ymm1(_mm256_load_ps(&y[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+0], _mm256_sub_ps(ymm0,ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+0], _mm256_sub_ps(ymm0, ymm1));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Re[i+8]));
		const __m256 ymm3(_mm256_load_ps(&y[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+8], _mm256_sub_ps(ymm2,ymm3));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+8], _mm256_sub_ps(ymm2, ymm3));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] - y[i];
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator*(const CV1DYMM8r4 &x,
		     const CV1DYMM8r4 &y) {
	if (x.m_nsize != y.m_nsize) { return (CV1DYMM8r4{}); }
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT	(ret_vec.m_nsize, 8); i += 16) {
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i+0]));
		const __m256 ymm1(_mm256_load_ps(&y.m_Re[i+0]));
		const __m256 ymm2(_mm256_load_ps(&x.m_Im[i+0]));
		const __m256 ymm3(_mm256_load_ps(&y.m_Im[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+0], _mm256_sub_ps(
			_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2, ymm3)));
		_mm256_stream_ps(&ret_vec.m_Im[i+0], _mm256_add_ps(
			_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0, ymm3)));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+0], _mm256_sub_ps(
			_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2, ymm3)));
		_mm256_store_ps(&ret_vec.m_Im[i+0], _mm256_add_ps(
			_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0, ymm3)));
#endif
		const __m256 ymm4(_mm256_load_ps(&x.m_Re[i+8]));
		const __m256 ymm5(_mm256_load_ps(&y.m_Im[i+8]));
		const __m256 ymm6(_mm256_load_ps(&x.m_Im[i+8]));
		const __m256 ymm7(_mm256_load_ps(&y.m_Im[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+8], _mm256_sub_ps(
			_mm256_mul_ps(ymm4, ymm5), _mm256_mul_ps(ymm6, ymm7)));
		_mm256_stream_ps(&ret_vec.m_Im[i+8], _mm256_add_ps(
			_mm256_mul_ps(ymm6, ymm5), _mm256_mul_ps(ymm4, ymm7)));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+8], _mm256_sub_ps(
			_mm256_mul_ps(ymm4, ymm5), _mm256_mul_ps(ymm6, ymm7)));
		_mm256_store_ps(&ret_vec.m_Im[i+8], _mm256_add_ps(
			_mm256_mul_ps(ymm6, ymm5), _mm256_mul_ps(ymm4, ymm7)));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = (x.m_Re[i] * y.m_Re[i]) - (x.m_Im[i] * y.m_Im[i]);
		ret_vec.m_Im[i] = (x.m_Im[i] * y.m_Re[i]) + (x.m_Re[i] * y.m_Im[i]);
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator*(const CV1DYMM8r4 &x,
		     const float * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) {return (CV1DYMM8r4{});}
	CV1DYMM8r4 ret_vec{x.m_nsize}; 
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 16) {
		__m256 ymm0(_mm256_load_ps(&x.m_Re[i+0]));
		__m256 ymm1(_mm256_load_ps(&y[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+0], _mm256_mul_ps(ymm0, ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+0], _mm256_mul_ps(ymm0, ymm1));
#endif
		__m256 ymm2(_mm256_load_ps(&x.m_Re[i+8]));
		__m256 ymm3(_mm256_load_ps(&y[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i+8], _mm256_mul_ps(ymm2, ymm3));
#else
		_mm256_store_ps(&ret_vec.m_Re[i+8], _mm256_mul_ps(ymm2, ymm3));
#endif
		__m256 ymm4(_mm256_load_ps(&x.m_Im[i+0]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i+0], _mm256_mul_ps(ymm4, ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Im[i+0], _mm256_mul_ps(ymm4, ymm1));
#endif
		__m256 ymm5(_mm256_load_ps(&x.m_Im[i+8]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i+8], _mm256_mul_ps(ymm5, ymm3));
#else
		_mm256_store_ps(&ret_vec.m_Im[i+8], _mm256_mul_ps(ymm5, ymm3));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] * y[i];
		ret_vec.m_Im[i] = x.m_Im[i] * y[i];
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator/(const CV1DYMM8r4 &x,
		     const CV1DYMM8r4 &y) {
	if (x.m_nsize != y.m_nsize) { return (CV1DYMM8r4{}); }
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 8) {
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i]));
		const __m256 ymm1(_mm256_load_ps(&y.m_Im[i]));
		const __m256 ymm2(_mm256_load_ps(&x.m_Im[i]));
		const __m256 re_term1 = _mm256_add_ps(
			_mm256_mul_ps(ymm0, ymm1), _mm256_mul_ps(ymm2,ymm1));
		const __m256 re_term2 = _mm256_add_ps(
			_mm256_mul_ps(ymm2, ymm1), _mm256_mul_ps(ymm0,ymm1));
		const __m256 ymm3(_mm256_load_ps(&y.m_Re[i]));
		const __m256 den_term = _mm256_add_ps(
			_mm256_mul_ps(ymm3, ymm3), _mm256_mul_ps(ymm1,ymm1));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i], _mm256_div_ps(re_term1,den_term));
#else
		_mm256_store_ps(&ret_vec.m_Re[i], _mm256_div_ps(re_term1, den_term));
#endif
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i], _mm256_div_ps(re_term2,den_term));
#else
		_mm256_store_ps(&ret_vec.m_Im[i], _mm256_div_ps(re_term2, den_term));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		const float tre = (x.m_Re[i] * y.m_Im[i]) + (x.m_Im[i] * y.m_Im[i]);
		const float tim = (x.m_Im[i] * y.m_Im[i]) - (x.m_Re[i] * y.m_Im[i]);
		const float den = (y.m_Re[i] * y.m_Re[i]) + (y.m_Im[i] * y.m_Im[i]);
		ret_vec.m_Re[i] = tre / den;
		ret_vec.m_Im[i] = tim / den;
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator/(const CV1DYMM8r4 &x,
		     const float * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) { return (CV1DYMM8r4{}); }
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(x.m_nsize,8); i += 8) {
		const __m256 ymm0(_mm256_load_ps(&x.m_Re[i]));
		const __m256 ymm1(_mm256_load_ps(&y[i]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i], _mm256_div_ps(ymm0, ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Re[i],  _mm256_div_ps(ymm0, ymm1));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Im[i]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i], _mm256_div_ps(ymm2, ymm1));
#else
		_mm256_store_ps(&ret_vec.m_Im[i], _mm256_div_ps(ymm2, ymm1));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		ret_vec.m_Re[i] = x.m_Re[i] / y[i];
		ret_vec.m_Im[i] = x.m_Im[i] / y[i];
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator==(const CV1DYMM8r4 &x,
		      const CV1DYMM8r4 &y) {
	using namespace gms::common;
	using namespace gms::math::constants;
	if (x.m_nsize != y.m_nsize) { return (CV1DYMM8r4{}); }
	CV1DYMM8r4 ret_vec{x.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT(ret_vec.m_nsize,8); i += 8) {
		const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i]);
		const __m256 ymm1 = _mm256_load_ps(&y.m_Re[i]);
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i], _mm256_cmp_ps(ymm0,ymm1,_CMP_EQ_OQ));
#else
		_mm256_store_ps(&ret_vec.m_Re[i], _mm256_cmp_ps(ymm0, ymm1, _CMP_EQ_OQ));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Im[i]));
		const __m256 ymm3(_mm256_load_ps(&y.m_Im[i]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i], _mm256_cmp_ps(ymm2,ymm3, _CMP_EQ_OQ));
#else
		_mm256_store_ps(&ret_vec.m_Im[i], _mm256_cmp_ps(ymm2, ymm3, _CMP_EQ_OQ));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		if (approximately_equalf32(x.m_Re[i], 
					y.m_Re[i],std::numeric_limits<float>::epsilon())) {
			ret_vec.m_Re[i] = 1.0;
		 }
		 else {
			 ret_vec.m_Re[i] = 0.0;
		 }
		 if (approximately_equalf32(x.m_Im[i],
					y.m_Im[i],std::numeric_limits<float>::epsilon())) {
			 ret_vec.m_Im[i] = 1.0;
		 }
		 else {
			 ret_vec.m_Im[i] = 0.0;
		 }
	}
	return (ret_vec);
}

gms::math::CV1DYMM8r4
gms::math::operator!=(const CV1DYMM8r4 &x,
		      const CV1DYMM8r4 &y){
	using namespace gms::common;
	using namespace gms::math::constants;
	if (x.m_nsize != y.m_nsize) { return (CV1DYMM8r4{}); }
	CV1DYMM8r4 ret_vec{ x.m_nsize };
	int32_t i;
	for (i = 0; i != ROUND_TO_EIGHT	(ret_vec.m_nsize, 8); i += 8) {
		const __m256 ymm0 = _mm256_load_ps(&x.m_Re[i]);
		const __m256 ymm1 = _mm256_load_ps(&y.m_Re[i]);
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Re[i], _mm256_cmp_ps(ymm0, ymm1, _CMP_NEQ_OQ));
#else
		_mm256_store_ps(&ret_vec.m_Re[i],  _mm256_cmp_ps(ymm0, ymm1,  _CMP_NEQ_OQ));
#endif
		const __m256 ymm2(_mm256_load_ps(&x.m_Im[i]));
		const __m256 ymm3(_mm256_load_ps(&y.m_Im[i]));
#if (USE_COMPLEX_VEC1D_YMM4R8_NT_STORES) == 1
		_mm256_stream_ps(&ret_vec.m_Im[i], _mm256_cmp_ps(ymm2, ymm3, _CMP_NEQ_OQ));
#else
		_mm256_store_ps(&ret_vec.m_Im[i],  _mm256_cmp_ps(ymm2, ymm3, _CMP_NEQ_OQ));
#endif
	}
	for (; i != ret_vec.m_nsize; ++i) {
		if (!approximately_equalf32(x.m_Re[i],
			y.m_Re[i], std::numeric_limits<float>::epsilon())) {
			ret_vec.m_Re[i] = 1.0;
		}
		else {
			ret_vec.m_Re[i] = 0.0;
		}
		if (!approximately_equalf32(x.m_Im[i],
			y.m_Im[i], std::numeric_limits<float>::epsilon())) {
			ret_vec.m_Im[i] = 1.0;
		}
		else {
			ret_vec.m_Im[i] = 0.0;
		}
	}
}




void
gms::math
::v256cnormalize_product( CV1DYMM8r4 &out,
			  const CV1DYMM8r4 &v1,
			  const CV1DYMM8r4 &v2,
			  const bool do_nt_store) {
	avx256_cnormalize_prod<CV1DYMM8r4>(out,v1,v2,do_nt_store);
}

void
gms::math::v256cmean_product(std::complex<float> &mean,
			     const CV1DYMM8r4 &v1,
			     const CV1DYMM8r4 &v2) {
	avx256_cmean_prod<CV1DYMM8r4>(mean,v1,v2);
}

void
gms::math::v256cmean_quotient(std::complex<float> &mean,
			      const CV1DYMM8r4 &v1,
			      const CV1DYMM8r4 &v2) {
	avx256_cmean_quot<CV1DYMM8r4>(mean,v1,v2);
}

void
gms::math::v256cconj_product(CV1DYMM8r4 &out,
			     const CV1DYMM8r4 &v1,
			     const CV1DYMM8r4 &v2,
			     const bool do_nt_store) {
	avx256_cconj_prod<CV1DYMM8r4>(out,v1,v2,do_nt_store);
}

void
gms::math::v256cnorm_conjprod(CV1DYMM8r4 &out,
			      const CV1DYMM8r4 &v1,
			      const CV1DYMM8r4 &v2,
			      const bool do_nt_store) {
	avx256_cnorm_conjprod<CV1DYMM8r4>(out,v1,v2,do_nt_store);
}

void
gms::math::v256cmean_conjprod(std::complex<float> &mean,
			      const CV1DYMM8r4 &v1,
			      const CV1DYMM8r4 &v2) {
	avx256_cmean_conjprod<CV1DYMM8r4>(mean,v1,v2);
}

void
gms::math::v256c_arithmean(std::complex<float> &mean,
			   const CV1DYMM8r4 &v) {
	avx256_arith_mean<CV1DYMM8r4>(mean,v);
}

void
gms::math::v256c_normalize(CV1DYMM8r4 &out,
			   const CV1DYMM8r4 &v1,
			   const CV1DYMM8r4 &v2,
			   const bool do_nt_store) {
	avx256_cnormalize<CV1DYMM8r4>(out,v1,v2,do_nt_store);
}

void
gms::math::v256c_magnitude( CV1DYMM8r4 &out,
			    const CV1DYMM8r4 &v1,
			    const CV1DYMM8r4 &v2) {
	avx256_cmagnitude<CV1DYMM8r4>(out,v1,v2);
}























	








