
#include <iomanip>
#include <complex>
#include <immintrin.h>
#include "GMS_avxcomplex.h"
#if defined _WIN64
    #if (GMS_DEBUG_ON) == 1
        #include "../GMS_debug.h"
    #else
        #include "../GMS_malloc.h"
    #endif
#elif defined __linux
    #if (GMS_DEBUG_ON) == 1
        #include "GMS_debug.h"
    #else
        #include "GMS_malloc.h"
    #endif
#endif
#if defined _WIN64
    #include "../GMS_common.h"
    #include "../Math/GMS_constants.h"
#elif defined __linux
    #include "GMS_common.h"
    #include "GMS_constants.h"
#endif
//
//	Implementation
//

// Parametrized macros
#if (GMS_MAN_PREFETCH) == 1

#if !defined (AVXCOMPLEX_PREFETCH_FROM_OBJ)
#define AVXCOMPLEX_PREFETCH_FROM_OBJ(obj,idx,off,hint) \
    _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Re[(idx)+(off)]),(hint)); \
    _mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Im[(idx)+(off)]),(hint));
#endif

#if !defined (AVXCOMPLEX_PREFETCH_FROM_PTR)
#define AVXCOMPLEX_PREFETCH_FROM_PTR(ptr,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(ptr)[(idx)+(off)]),(hint));
#endif

#if !defined (AVCOMPLEX_PREFETCH_FROM_THIS)
#define AVXCOMPLEX_PREFETCH_FROM_THIS(obj,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj)->m_Re[(idx)+(off)]),(hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj)->m_Im[(idx)+(off)]),(hint));
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



gms::math::AVXVComplex1D::
AVXVComplex1D(){
	data.m_Re = NULL;
	data.m_Im = NULL;
	data.m_nsize = 0;
}



gms::math::AVXVComplex1D::
AVXVComplex1D(const int32_t nsize){
	using namespace gms::common;
#if defined _WIN64
    #if (GMS_DEBUG_ON) == 1
	data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B,__FILE__,__LINE__);
	data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B,__FILE__,__LINE__);
    #else
	data.m_Re = gms_edmalloca(static_cast<size_t>(nsize),align64B);
	data.m_Im = gms_edmalloca(static_cast<size_t>(nsize),align64B);
    #endif
#elif defined __linux
    #if (GMS_DEBUG_ON) == 1
	data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B);
	data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B);
    #else
	data.m_Re = gms_dmalloca(static_cast<size_t>(nsize),align64B);
	data.m_Im = gms_dmalloca(static_cast<size_t>(nsize),align64B);
    #endif
#endif
	data.m_nsize = nsize;
	// Memory first touch.
	avx256_init_unroll8x_pd(&data.m_Re[0],data.m_nsize,0.0); 
	avx256_init_unroll8x_pd(&data.m_Im[0],data.m_nsize,0.0);
}	
	
	
gms::math::AVXVComplex1D
::AVXVComplex1D(const double * __restrict Re,
		const double * __restrict Im,
		const int32_t nsize) {
	using namespace gms::common;
#if defined _WIN64
     #if (GMS_DEBUG_ON) == 1
	   data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B,__FILE__,__LINE__);
	   data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B,__FILE__,__LINE__);
     #else
	   data.m_Re = gms_edmalloca(static_cast<size_t>(nsize), align64B);
	   data.m_Im = gms_edmalloca(static_cast<size_t>(nsize), align64B);
     #endif
#elif defined __linux
     #if (GMS_DEBUG_ON) == 1
	   data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B);
	   data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(nsize),align64B);
     #else
	   data.m_Re = gms_edmalloca(static_cast<size_t>(nsize),align64B);
	   data.m_Im = gms_edmalloca(static_cast<size_t>(nsize),align64B);
     #endif
#endif
        data.m_nsize = nsize;
#if (USE_AVXCOMPLEX_NT_STORES) == 1
	avx256_uncached_memmove(&data.m_Re[0], &Re[0],data.m_nsize);
	avx256_uncached_memmove(&data.m_Im[0], &Im[0],data.m_nsize);
#else
	avx256_cached_memmove(&data.m_Re[0], &Re[0], data.m_nsize);
	avx256_cached_memmove(&data.m_Im[0], &Im[0], data.m_nsize);
#endif
}


	




gms::math::AVXVComplex1D::
AVXVComplex1D(const AVXVComplex1D &x){
	using namespace gms::common;
#if defined _WIN64	
    #if (GMS_DEBUG_ON) == 1
	data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize), align64B, __FILE__, __LINE__);
	data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize), align64B, __FILE__, __LINE__);
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
    data.m_nsize = x.data.m_nsize;
#if (USE_AVXCOMPLEX_NT_STORES) == 1
	avx256_uncached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
	avx256_uncached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#else
	avx256_cached_memmove(&data.m_Re[0], &x.data.m_Re[0], x.data.m_nsize);
	avx256_cached_memmove(&data.m_Im[0], &x.data.m_Im[0], x.data.m_nsize);
#endif
}

gms::math::AVXVComplex1D::
AVXVComplex1D(AVXVComplex1D &&x) {
	data.m_Re = &x.data.m_Re[0];
	data.m_Im = &x.data.m_Im[0];
	data.m_nsize = x.data.m_nsize;
	x.data.m_Re = NULL;
	x.data.m_Im = NULL;
	x.data.m_nsize = 0;
}


gms::math::AVXVComplex1D::
		  ~AVXVComplex1D() {
#if (GMS_DEBUG_ON) == 1
	if (NULL != data.m_Re) _aligned_free_dbg(data.m_Re); data.m_Re = NULL;
	if (NULL != data.m_Im) _aligned_free_dbg(data.m_Im); data.m_Im = NULL;
#else
	if (data.m_Re != NULL) _mm_free(data.m_Re); data.m_Re = NULL;
	if (data.m_Im != NULL) _mm_free(data.m_Im); data.m_Im = NULL;
#endif		  
}
		
	
	

gms::math::AVXVComplex1D &
gms::math::AVXVComplex1D
::operator=(const AVXVComplex1D &x) {
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
	data.m_Re = NULL;
	data.m_Im = NULL;
	data.m_nsize = 0;
#if defined _WIN64
    #if (GMS_DEBUG_ON) == 1
	  data.m_Re = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize), align64B, __FILE__, __LINE__);
	  data.m_Im = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nsize), align64B, __FILE__, __LINE__);
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
	
#if (USE_AVXCOMPLEX_NT_STORES) == 1
	 avx256_uncached_memmove(&data.m_Re[0],&x.data.m_Re[0],x.data.m_nsize);
	 avx256_uncached_memmove(&data.m_Im[0],&x.data.m_Im[0],x.data.m_nsize);
#else
	 avx256_cached_memmove(&data.m_Re[0], &x.data.m_Re[0],x.data.m_nsize);
	 avx256_cached_memmove(&data.m_Im[0], &x.data.m_Im[0],x.data.m_nsize);
#endif		
	}
	return (*this);
}

gms::math::AVXVComplex1D &
gms::math::AVXVComplex1D::operator=(AVXVComplex1D &&x) {
	if (this == &x) return (*this);
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
	data.m_Re    = &x.data.m_Re[0];
	data.m_Im    = &x.data.m_Im[0];
	data.m_nsize = x.data.m_nsize;
	x.data.m_Re = NULL;
	x.data.m_Im = NULL;
	x.data.m_nsize = 0;
	return (*this);
}	

std::ostream &
gms::math::operator<<(std::ostream &os,
		      const AVXVComplex1D &x) {
	for (int64_t i = 0LL; i != x.data.m_nsize; ++i) {
		os << std::fixed << std::showpoint << std::setprecision(15) <<
			  std::setw(4) <<  "Re: " << "{" << x.data.m_Re[i] << "}" <<
			  std::setw(12) << "Im: " << "{" << x.data.m_Im[i] << "}" << std::endl;
	}
	return (os);
}
	
gms::math::AVXVComplex1D
gms::math::operator+(const AVXVComplex1D &x,
		     const AVXVComplex1D &y) {
	if (x.data.m_nsize != y.data.m_nsize) {
		return (AVXVComplex1D{});
	}
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize,4); i += 8) {
		// Linearly growing indices, no need for software prefetch.
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y.data.m_Re[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+0], _mm256_add_pd(ymm0,ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i+4]));
		const __m256d ymm3(_mm256_load_pd(&y.data.m_Re[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+4], _mm256_add_pd(ymm2,ymm3));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));
#endif
		const __m256d ymm4(_mm256_load_pd(&x.data.m_Im[i+0]));
		const __m256d ymm5(_mm256_load_pd(&y.data.m_Im[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i+0], _mm256_add_pd(ymm4,ymm5));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i+0], _mm256_add_pd(ymm4, ymm5));
#endif
		const __m256d ymm6(_mm256_load_pd(&x.data.m_Im[i+4]));
		const __m256d ymm7(_mm256_load_pd(&y.data.m_Im[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i+4], _mm256_add_pd(ymm6,ymm7));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(ymm6, ymm7));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] + y.data.m_Re[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] + y.data.m_Re[i];
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator+(const AVXVComplex1D &x,
		     const double * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) {
		return (AVXVComplex1D{});
	}
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize,4); i += 8) {
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+0],_mm256_add_pd(ymm0,ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_add_pd(ymm0, ymm1));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i+4]));
		const __m256d ymm3(_mm256_load_pd(&y[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+4],_mm256_add_pd(ymm2,ymm3));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_add_pd(ymm2, ymm3));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] + y[i];
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator-(const AVXVComplex1D &x,
		     const AVXVComplex1D &y) {
	if (x.data.m_nsize != y.data.m_nsize) { return (AVXVComplex1D{}); }
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y.data.m_Re[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+0], _mm256_sub_pd(ymm0,ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i+4]));
		const __m256d ymm3(_mm256_load_pd(&y.data.m_Re[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+4], _mm256_sub_pd(ymm2,ymm3));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));
#endif
		const __m256d ymm4(_mm256_load_pd(&x.data.m_Im[i+0]));
		const __m256d ymm5(_mm256_load_pd(&y.data.m_Im[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i+0], _mm256_sub_pd(ymm4,ymm5));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_sub_pd(ymm4, ymm5));
#endif
		const __m256d ymm6(_mm256_load_pd(&x.data.m_Im[i+4]));
		const __m256d ymm7(_mm256_load_pd(&y.data.m_Im[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i+4], _mm256_sub_pd(ymm6,ymm7));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_sub_pd(ymm6, ymm7));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] - y.data.m_Re[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] - y.data.m_Im[i];
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator-(const AVXVComplex1D &x,
		     const double * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) {return (AVXVComplex1D{});}
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+0], _mm256_sub_pd(ymm0,ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(ymm0, ymm1));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Re[i+4]));
		const __m256d ymm3(_mm256_load_pd(&y[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+4], _mm256_sub_pd(ymm2,ymm3));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(ymm2, ymm3));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] - y[i];
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator*(const AVXVComplex1D &x,
		     const AVXVComplex1D &y) {
	if (x.data.m_nsize != y.data.m_nsize) { return (AVXVComplex1D{}); }
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 8) {
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y.data.m_Re[i+0]));
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Im[i+0]));
		const __m256d ymm3(_mm256_load_pd(&y.data.m_Im[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(
			_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm3)));
		_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_add_pd(
			_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm3)));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_sub_pd(
			_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2, ymm3)));
		_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_add_pd(
			_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0, ymm3)));
#endif
		const __m256d ymm4(_mm256_load_pd(&x.data.m_Re[i+4]));
		const __m256d ymm5(_mm256_load_pd(&y.data.m_Im[i+4]));
		const __m256d ymm6(_mm256_load_pd(&x.data.m_Im[i+4]));
		const __m256d ymm7(_mm256_load_pd(&y.data.m_Im[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(
			_mm256_mul_pd(ymm4, ymm5), _mm256_mul_pd(ymm6, ymm7)));
		_mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(
			_mm256_mul_pd(ymm6, ymm5), _mm256_mul_pd(ymm4, ymm7)));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_sub_pd(
			_mm256_mul_pd(ymm4, ymm5), _mm256_mul_pd(ymm6, ymm7)));
		_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_add_pd(
			_mm256_mul_pd(ymm6, ymm5), _mm256_mul_pd(ymm4, ymm7)));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = (x.data.m_Re[i] * y.data.m_Re[i]) - (x.data.m_Im[i] * y.data.m_Im[i]);
		ret_vec.data.m_Im[i] = (x.data.m_Im[i] * y.data.m_Re[i]) + (x.data.m_Re[i] * y.data.m_Im[i]);
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator*(const AVXVComplex1D &x,
		     const double * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) {return (AVXVComplex1D{});}
	AVXVComplex1D ret_vec{x.data.m_nsize}; // <<-- here occurs Memory first touch.
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize,4); i += 8) {
		__m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		__m256d ymm1(_mm256_load_pd(&y[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i + 0], _mm256_mul_pd(ymm0, ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_mul_pd(ymm0, ymm1));
#endif
		__m256d ymm2(_mm256_load_pd(&x.data.m_Re[i+4]));
		__m256d ymm3(_mm256_load_pd(&y[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm3));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 4], _mm256_mul_pd(ymm2, ymm3));
#endif
		__m256d ymm4(_mm256_load_pd(&x.data.m_Im[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i + 0], _mm256_mul_pd(ymm4, ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_mul_pd(ymm4, ymm1));
#endif
		__m256d ymm5(_mm256_load_pd(&x.data.m_Im[i+4]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i + 4], _mm256_mul_pd(ymm5, ymm3));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i + 4], _mm256_mul_pd(ymm5, ymm3));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] * y[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] * y[i];
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator/(const AVXVComplex1D &x,
		     const AVXVComplex1D &y) {
	if (x.data.m_nsize != y.data.m_nsize) { return (AVXVComplex1D{}); }
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 4) {
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y.data.m_Im[i+0]));
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Im[i+0]));
		const __m256d re_term1 = _mm256_add_pd(
			_mm256_mul_pd(ymm0, ymm1), _mm256_mul_pd(ymm2,ymm1));
		const __m256d re_term2 = _mm256_add_pd(
			_mm256_mul_pd(ymm2, ymm1), _mm256_mul_pd(ymm0,ymm1));
		const __m256d ymm3(_mm256_load_pd(&y.data.m_Re[i+0]));
		const __m256d den_term = _mm256_add_pd(
			_mm256_mul_pd(ymm3, ymm3), _mm256_mul_pd(ymm1,ymm1));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i+0], _mm256_div_pd(re_term1,den_term));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i + 0], _mm256_div_pd(re_term1, den_term));
#endif
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i+0], _mm256_div_pd(re_term2,den_term));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i + 0], _mm256_div_pd(re_term2, den_term));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		const double tre = (x.data.m_Re[i] * y.data.m_Im[i]) + (x.data.m_Im[i] * y.data.m_Im[i]);
		const double tim = (x.data.m_Im[i] * y.data.m_Im[i]) - (x.data.m_Re[i] * y.data.m_Im[i]);
		const double den = (y.data.m_Re[i] * y.data.m_Re[i]) + (y.data.m_Im[i] * y.data.m_Im[i]);
		ret_vec.data.m_Re[i] = tre / den;
		ret_vec.data.m_Im[i] = tim / den;
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator/(const AVXVComplex1D &x,
		     const double * __restrict y) {
	using namespace gms::common;
	if (!Is_ptr_aligned32(y)) { return (AVXVComplex1D{}); }
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(x.data.m_nsize,4); i += 4) {
		const __m256d ymm0(_mm256_load_pd(&x.data.m_Re[i+0]));
		const __m256d ymm1(_mm256_load_pd(&y[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i], _mm256_div_pd(ymm0, ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i],  _mm256_div_pd(ymm0, ymm1));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Im[i+0]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i], _mm256_div_pd(ymm2, ymm1));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i], _mm256_div_pd(ymm2, ymm1));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		ret_vec.data.m_Re[i] = x.data.m_Re[i] / y[i];
		ret_vec.data.m_Im[i] = x.data.m_Im[i] / y[i];
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator==(const AVXVComplex1D &x,
		      const AVXVComplex1D &y) {
	using namespace gms::common;
	using namespace gms::math::constants;
	if (x.data.m_nsize != y.data.m_nsize) { return (AVXVComplex1D{}); }
	AVXVComplex1D ret_vec{x.data.m_nsize};
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize,4); i += 4) {
		const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i]);
		const __m256d ymm1 = _mm256_load_pd(&y.data.m_Re[i]);
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i], _mm256_cmp_pd(ymm0,ymm1,_CMP_EQ_OQ));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i], _mm256_cmp_pd(ymm0, ymm1, _CMP_EQ_OQ));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Im[i]));
		const __m256d ymm3(_mm256_load_pd(&y.data.m_Im[i]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i], _mm256_cmp_pd(ymm2,ymm3, _CMP_EQ_OQ));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i], _mm256_cmp_pd(ymm2, ymm3, _CMP_EQ_OQ));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		if (approximately_equalf64(x.data.m_Re[i], 
								   y.data.m_Re[i], DEPS)) {
			ret_vec.data.m_Re[i] = 1.0;
		 }
		 else {
			 ret_vec.data.m_Re[i] = 0.0;
		 }
		 if (approximately_equalf64(x.data.m_Im[i],
									y.data.m_Im[i],DEPS)) {
			 ret_vec.data.m_Im[i] = 1.0;
		 }
		 else {
			 ret_vec.data.m_Im[i] = 0.0;
		 }
	}
	return (ret_vec);
}

gms::math::AVXVComplex1D
gms::math::operator!=(const AVXVComplex1D &x,
		      const AVXVComplex1D &y){
	using namespace gms::common;
	using namespace gms::math::constants;
	if (x.data.m_nsize != y.data.m_nsize) { return (AVXVComplex1D{}); }
	AVXVComplex1D ret_vec{ x.data.m_nsize };
	int32_t i;
	for (i = 0; i != ROUND_TO_FOUR(ret_vec.data.m_nsize, 4); i += 4) {
		const __m256d ymm0 = _mm256_load_pd(&x.data.m_Re[i]);
		const __m256d ymm1 = _mm256_load_pd(&y.data.m_Re[i]);
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Re[i], _mm256_cmp_pd(ymm0, ymm1, _CMP_NEQ_OQ));
#else
		_mm256_store_pd(&ret_vec.data.m_Re[i],  _mm256_cmp_pd(ymm0, ymm1,  _CMP_NEQ_OQ));
#endif
		const __m256d ymm2(_mm256_load_pd(&x.data.m_Im[i]));
		const __m256d ymm3(_mm256_load_pd(&y.data.m_Im[i]));
#if (USE_AVXCOMPLEX_NT_STORES) == 1
		_mm256_stream_pd(&ret_vec.data.m_Im[i], _mm256_cmp_pd(ymm2, ymm3, _CMP_NEQ_OQ));
#else
		_mm256_store_pd(&ret_vec.data.m_Im[i],  _mm256_cmp_pd(ymm2, ymm3, _CMP_NEQ_OQ));
#endif
	}
	for (; i != ret_vec.data.m_nsize; ++i) {
		if (!approximately_equalf64(x.data.m_Re[i],
			y.data.m_Re[i], DEPS)) {
			ret_vec.data.m_Re[i] = 1.0;
		}
		else {
			ret_vec.data.m_Re[i] = 0.0;
		}
		if (!approximately_equalf64(x.data.m_Im[i],
			y.data.m_Im[i], DEPS)) {
			ret_vec.data.m_Im[i] = 1.0;
		}
		else {
			ret_vec.data.m_Im[i] = 0.0;
		}
	}
}


#include "LAM_avxcomplex_common.h"

void
gms::math
::v256cnormalize_product( AVXVComplex1D &out,
			  const AVXVComplex1D &v1,
			  const AVXVComplex1D &v2,
			  const bool do_nt_store) {
	avx256_cnormalize_prod<AVXVComplex1D>(out,v1,v2,do_nt_store);
}

void
gms::math::v256cmean_product(std::complex<double> &mean,
			     const AVXVComplex1D &v1,
			     const AVXVComplex1D &v2) {
	avx256_cmean_prod<AVXVComplex1D>(mean,v1,v2);
}

void
gms::math::v256cmean_quotient(std::complex<double> &mean,
			      const AVXVComplex1D &v1,
			      const AVXVComplex1D &v2) {
	avx256_cmean_quot<AVXVComplex1D>(mean,v1,v2);
}

void
gms::math::v256cconj_product(AVXVComplex1D &out,
			     const AVXVComplex1D &v1,
			     const AVXVComplex1D &v2,
			     const bool do_nt_store) {
	avx256_cconj_prod<AVXVComplex1D>(out,v1,v2,do_nt_store);
}

void
gms::math::v256cnorm_conjprod(AVXVComplex1D &out,
			      const AVXVComplex1D &v1,
			      const AVXVComplex1D &v2,
			      const bool do_nt_store) {
	avx256_cnorm_conjprod<AVXVComplex1D>(out,v1,v2,do_nt_store);
}

void
gms::math::v256cmean_conjprod(std::complex<double> &mean,
			      const AVXVComplex1D &v1,
			      const AVXVComplex1D &v2) {
	avx256_cmean_conjprod<AVXVComplex1D>(mean,v1,v2);
}

void
gms::math::v256c_arithmean(std::complex<double> &mean,
			   const AVXVComplex1D &v) {
	avx256_arith_mean<AVXVComplex1D>(mean,v);
}

void
gms::math::v256c_normalize(AVXVComplex1D &out,
			   const AVXVComplex1D &v1,
			   const AVXVComplex1D &v2,
			   const bool do_nt_store) {
	avx256_cnormalize<AVXVComplex1D>(out,v1,v2,do_nt_store);
}

void
gms::math::v256c_magnitude( AVXVComplex1D &out,
			    const AVXVComplex1D &v1,
			    const AVXVComplex1D &v2) {
	avx256_cmagnitude<AVXVComplex1D>(out,v1,v2);
}























	








