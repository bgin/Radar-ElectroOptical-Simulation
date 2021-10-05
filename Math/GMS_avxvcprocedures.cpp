
//
//	Implementation
//

#include "GMS_avxvcprocedures.h"

 #include "GMS_config.h"


#if !defined (AVXVCPROCEDURES_PREFETCH_FROM_OBJS)
#define AVXVCPROCEDURES_PREFETCH_FROM_OBJS(obj1,obj2,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj1).m_Re[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj2).m_Re[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj1).m_Im[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj2).m_Im[(idx)+(off)]), (hint));
#endif

#if !defined (AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR)
#define AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(obj,ptr,idx,off,hint) \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Re[(idx)+(off)]), (hint)); \
	_mm_prefetch(reinterpret_cast<const char*>(&(ptr)[(idx)+(off)]), (hint));	   \
	_mm_prefetch(reinterpret_cast<const char*>(&(obj).m_Im[(idx)+(off)]), (hint));
#endif



namespace {
	
	constexpr uint64_t OFF0 = 0ULL;
	constexpr uint64_t OFF4 = 4ULL;
}
void
gms::math::avxvcomplex_add(AVXVComplex1D &c,
			   const AVXVComplex1D &a,
			   const AVXVComplex1D &b) {
	uint64_t i;
#if (GMS_MAN_PREFETCH) == 1
size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#endif
	
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1
#pragma noprefetch
#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a,b,i,4ULL,_MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a,b,i,4ULL,_MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a,b,i,8ULL,_MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
			// Using AVX for comparison.
			// Unrolled 2x in order to match perfectly
			// per core 2 load operations per cycle.
			_mm256_store_pd(&c.m_Re[i+OFF0], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
		    _mm256_load_pd(&b.m_Re[i+OFF0])));
		    _mm256_store_pd(&c.m_Im[i + OFF4], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i+OFF4]),
			_mm256_load_pd(&b.m_Re[i+OFF4])));
			_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_add_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i+OFF0])));
			_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_add_pd(_mm256_load_pd(&a.m_Im[i+OFF4]),
			_mm256_load_pd(&b.m_Im[i+OFF4])));
#else
			//	Do not pollute cache on memory stores
			//  You must be sure that stores will not be
			// used later, so do not keep them in cache.
			_mm256_stream_pd(&c.m_Re[i + OFF0], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b.m_Re[i + OFF0])));
		    _mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Re[i + OFF4])));
		    _mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_add_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0])));
		    _mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_add_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
#endif
	}
	// Scalar remainder
	for (; i != c.m_nsize; ++i) {
		c.m_Re[i] = a.m_Re[i] + b.m_Re[i];
		c.m_Im[i] = a.m_Im[i] + b.m_Im[i];
	}
}

void
gms::math::avxvcomplex_add(AVXVComplex1D &c,
			   const AVXVComplex1D &a,
			   const double * __restrict b) {
	std::uint64_t i;
#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif	
	
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1

#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a,b,i,4ULL,_MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a,b,i,4ULL,_MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a,b,i,8ULL,_MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
			// Using AVX for comparison.
			// Unrolled 2x in order to match perfectly
			// per core 2 load operations per cycle.
			_mm256_store_pd(&c.m_Re[i+OFF0], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b[i+OFF0])));
		    _mm256_store_pd(&c.m_Re[i+OFF4], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i+OFF4]),
			_mm256_load_pd(&b[i+OFF4])));
#else
			//	Do not pollute cache on memory stores
			//  You must be sure that stores will not be
			// used later, so do not keep them in cache.
			_mm256_stream_pd(&c.m_Re[i + OFF0], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b[i + OFF0])));
		    _mm256_stream_pd(&c.m_Re[i + OFF4], _mm256_add_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));
#endif
	}
	for (; i != c.m_nsize; ++i)
		c.m_Re[i] = a.m_Re[i] + b[i];
}

void
gms::math::avxvcomplex_sub(AVXVComplex1D &c,
			   const AVXVComplex1D &a,
			   const AVXVComplex1D &b) {
#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#endif		
	uint64_t i;
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1
#pragma noprefetch
#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a,b,i,4ULL,_MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a,b,i,4ULL,_MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a,b,i,4ULL,_MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
			// Using AVX for comparison.
			// Unrolled 2x in order to match perfectly
			// per core 2 load operations per cycle.
			_mm256_store_pd(&c.m_Re[i + OFF0], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b.m_Re[i+OFF0])));
			_mm256_store_pd(&c.m_Re[i + OFF4], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i+OFF4]),
			_mm256_load_pd(&b.m_Re[i+OFF4])));

			_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_sub_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i+OFF0])));
			_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_sub_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
#else
			//	Do not pollute cache on memory stores
			//  You must be sure that stores will not be
			// used later, so do not keep them in cache.
			_mm256_stream_pd(&c.m_Re[i + OFF0], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b.m_Re[i + OFF0])));
			_mm256_stream_pd(&c.m_Re[i + OFF4], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Re[i + OFF4])));

			_mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_sub_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0])));
			_mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_sub_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
#endif
	}
	for (; i != c.m_nsize; ++i) {
		c.m_Re[i] = a.m_Re[i] - b.m_Re[i];
		c.m_Im[i] = a.m_Im[i] - b.m_Im[i];
	}
}

void
gms::math::avxvcomplex_sub(AVXVComplex1D &c, 
			   const AVXVComplex1D &a,
			   const double * __restrict b)  {
	std::uint64_t i;						
#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif	
	
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1

#if (CACHE_SIZE) < (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a,b,i,4ULL,_MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a,b,i,4ULL,_MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a,b,i,4ULL,_MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
			// Using AVX for comparison.
			// Unrolled 2x in order to match perfectly
			// per core 2 load operations per cycle.
			_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b[i+OFF0])));
			_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));
#else
			//	Do not pollute cache on memory stores
			//  You must be sure that stores will not be
			// used later, so do not keep them in cache.
			_mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b[i + OFF0])));
			_mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_sub_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));
#endif
	}
	for (; i != c.m_nsize; ++i) 
		c.m_Re[i] = a.m_Re[i] - b[i];
}

void
gms::math::avxvcomplex_mul(AVXVComplex1D &c,
			   const AVXVComplex1D &a,
			   const AVXVComplex1D &b) {
	std::uint64_t i;

#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch 
#endif	
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1 

#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
		// Using AVX for comparison.
		// Unrolled 2x in order to match perfectly
		// per core 2 load operations per cycle.
		_mm256_store_pd(&c.m_Re[i + OFF0], 
				_mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
					_mm256_load_pd(&b.m_Re[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
									    _mm256_load_pd(&b.m_Im[i+OFF0]))));
		_mm256_store_pd(&c.m_Re[i + OFF4], 
				_mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
					_mm256_load_pd(&b.m_Re[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
									    _mm256_load_pd(&b.m_Im[i + OFF4]))));

		_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
		_mm256_load_pd(&b.m_Im[i+OFF0]))));
		_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4]))));
#else
		//	Do not pollute cache on memory stores
		//  You must be sure that stores will not be
		// used later, so do not keep them in cache.
		_mm256_stream_pd(&c.m_Re[i + OFF0], _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
		_mm256_load_pd(&b.m_Re[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0]))));
		_mm256_stream_pd(&c.m_Re[i + OFF4], _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
		_mm256_load_pd(&b.m_Re[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4]))));

		_mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0]))));
		_mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4]))));
#endif
	}
	for (; i != c.m_nsize; ++i) {
		c.m_Re[i] = (a.m_Re[i] * b.m_Re[i]) - (a.m_Im[i] * b.m_Im[i]);
		c.m_Im[i] = (a.m_Im[i] * b.m_Re[i]) + (a.m_Re[i] * b.m_Im[i]);
	}
}

void
gms::math::avxvcomplex_mul(AVXVComplex1D &c,
			   const AVXVComplex1D &a,
			   const double * __restrict b) {
	std::uint64_t i;

#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif	
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1

#if (CACHE_SIZE) < (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a, b, i, 4ULL, _MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a, b, i, 4ULL, _MM_HINT_T0)
			AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a, b, i, 4ULL, _MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
			// Using AVX for comparison.
			// Unrolled 2x in order to match perfectly
			// per core 2 load operations per cycle.
			_mm256_store_pd(&c.m_Im[i+OFF0], _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b[i+OFF0])));
			_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));

			_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
			_mm256_load_pd(&b[i+OFF0])));
			_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));
#else
			//	Do not pollute cache on memory stores
			//  You must be sure that stores will not be
			// used later, so do not keep them in cache.
			_mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b[i + OFF0])));
			_mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));

		    _mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
			_mm256_load_pd(&b[i + OFF0])));
		    _mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b[i + OFF4])));
#endif
	}
	for (; i != c.m_nsize; ++i) {
		c.m_Re[i] = a.m_Re[i] * b[i];
		c.m_Im[i] = a.m_Im[i] * b[i];
	}
}

void
::math::avxcomplex_div(AVXVComplex1D &c,
		       const AVXVComplex1D &a,
		       const AVXVComplex1D &b,
		       const AVXVComplex1D &tmp) {
	std::uint64_t i;
#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize,4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1 

#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
		// Using AVX for comparison.
		// Unrolled 2x in order to match perfectly
		// per core 2 load operations per cycle.
		const __m256d tre1 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i+OFF0])));
		const __m256d tre2 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
		const __m256d tim1 = _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i+OFF0])));
		const __m256d tim2 = _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
		const __m256d den1 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&b.m_Re[i+OFF0]),
			_mm256_load_pd(&tmp.m_Re[i+OFF0])), _mm256_mul_pd(_mm256_load_pd(&b.m_Im[i+OFF0]),
			_mm256_load_pd(&tmp.m_Im[i+OFF0])));
		const __m256d den2 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&b.m_Re[i + OFF4]),
			_mm256_load_pd(&tmp.m_Re[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&b.m_Im[i + OFF4]),
			_mm256_load_pd(&tmp.m_Im[i + OFF4])));

		_mm256_store_pd(&c.m_Re[i + OFF0], _mm256_div_pd(tre1,den1));
		_mm256_store_pd(&c.m_Re[i + OFF4], _mm256_div_pd(tre2,den2));
		_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_div_pd(tim1,den1));
		_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_div_pd(tim2,den2));
#else
		//	Do not pollute cache on memory stores
		//  You must be sure that stores will not be
		// used later, so do not keep them in cache.
		const __m256d tre1 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0])));
		const __m256d tre2 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
		const __m256d tim1 = _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0])));
		const __m256d tim2 = _mm256_sub_pd(_mm256_mul_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4])));
		const __m256d den1 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&b.m_Re[i + OFF0]),
			_mm256_load_pd(&tmp.m_Re[i + OFF0])), _mm256_mul_pd(_mm256_load_pd(&b.m_Im[i + OFF0]),
			_mm256_load_pd(&tmp.m_Im[i + OFF0])));
		const __m256d den2 = _mm256_add_pd(_mm256_mul_pd(_mm256_load_pd(&b.m_Re[i + OFF4]),
			_mm256_load_pd(&tmp.m_Re[i + OFF4])), _mm256_mul_pd(_mm256_load_pd(&b.m_Im[i + OFF4]),
			_mm256_load_pd(&tmp.m_Im[i + OFF4])));

		_mm256_stream_pd(&c.m_Re[i + OFF0], _mm256_div_pd(tre1, den1));
		_mm256_stream_pd(&c.m_Re[i + OFF4], _mm256_div_pd(tre2, den2));
		_mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_div_pd(tim1, den1));
		_mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_div_pd(tim2, den2));
#endif		
	}
	for (; i != c.m_nsize; ++i) {
		const double re =  (a.m_Re[i] * b.m_Im[i])   + (a.m_Im[i] * b.m_Im[i]);
		const double im =  (a.m_Im[i] * b.m_Im[i])   - (a.m_Re[i] * b.m_Im[i]);
		const double den = (b.m_Re[i] * tmp.m_Re[i]) + (b.m_Im[i] * tmp.m_Im[i]);
		c.m_Re[i] = re / den;
		c.m_Im[i] = im / den;
	}
}

void
gms::math::avxcomplex_div(AVXVComplex1D &c,
			  const AVXVComplex1D &a,
			  const double * __restrict b) {
	std::uint64_t i;
#if (GMS_MAN_PREFETCH) == 1
	size_t len = c.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif
	for (i = 0ULL; i != ROUND_TO_FOUR(c.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1
#if (GMS_MAN_PREFETCH) == 1

#if (CACHE_SIZE) < (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a, b, i, 4ULL, _MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a, b, i, 4ULL, _MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJ_PTR(a, b, i, 4ULL, _MM_HINT_T1)
#endif
#endif	 

#endif
#if (GMS_CACHE_MEM_STORES) == 1
		// Unrolled 2x in order to match perfectly
		// per core 2 load operations per cycle.
		// At expense of divider stall.
		_mm256_store_pd(&c.m_Re[i + OFF0], _mm256_div_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
		_mm256_load_pd(&b[i+OFF0])));
		_mm256_store_pd(&c.m_Re[i + OFF4], _mm256_div_pd(_mm256_load_pd(&a.m_Re[i+OFF4]),
		_mm256_load_pd(&b[i+OFF4])));

		_mm256_store_pd(&c.m_Im[i + OFF0], _mm256_div_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
		_mm256_load_pd(&b[i+OFF0])));
		_mm256_store_pd(&c.m_Im[i + OFF4], _mm256_div_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b[i + OFF4])));
#else
		//	Do not pollute cache on memory stores
		//  You must be sure that stores will not be
		// used later, so do not keep them in cache.
		_mm256_stream_pd(&c.m_Re[i + OFF0], _mm256_div_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
		_mm256_load_pd(&b[i + OFF0])));
		_mm256_stream_pd(&c.m_Re[i + OFF4], _mm256_div_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
		_mm256_load_pd(&b[i + OFF4])));

		_mm256_stream_pd(&c.m_Im[i + OFF0], _mm256_div_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
			_mm256_load_pd(&b[i + OFF0])));
		_mm256_stream_pd(&c.m_Im[i + OFF4], _mm256_div_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b[i + OFF4])));
#endif
	}
	for (; i != c.m_nsize; ++i) {
		c.m_Re[i] = a.m_Re[i] / b[i];
		c.m_Im[i] = a.m_Im[i] / b[i];
	}
}

void
gms::math::avxvcomplex_eq(double * __restrict req,
			  double * __restrict imeq,
			  const AVXVComplex1D &a,
			  const AVXVComplex1D &b) {
	uint64_t i;

#if (GMS_MAN_PREFETCH) == 1
	size_t len = a.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif	

	for (i = 0ULL; i != ROUND_TO_FOUR(a.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1 

#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
			AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
			// Unrolled 2x in order to match perfectly
			// per core 2 load operations per cycle.
			_mm256_store_pd(&req[i+OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
			_mm256_load_pd(&b.m_Re[i+OFF0]),_CMP_EQ_OQ));
			_mm256_store_pd(&req[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Re[i + OFF4]), _CMP_EQ_OQ));

			_mm256_store_pd(&imeq[i + OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
			_mm256_load_pd(&b.m_Im[i+OFF0]),_CMP_EQ_OQ));
			_mm256_store_pd(&imeq[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4]), _CMP_EQ_OQ));
#else
			//	Do not pollute cache on memory stores
			//  You must be sure that stores will not be
			// used later, so do not keep them in cache.
			_mm256_stream_pd(&req[i + OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
			_mm256_load_pd(&b.m_Re[i + OFF0]), _CMP_EQ_OQ));
			_mm256_stream_pd(&req[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
			_mm256_load_pd(&b.m_Re[i + OFF4]), _CMP_EQ_OQ));

			_mm256_stream_pd(&imeq[i + OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
			_mm256_load_pd(&b.m_Im[i + OFF0]), _CMP_EQ_OQ));
			_mm256_stream_pd(&imeq[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
			_mm256_load_pd(&b.m_Im[i + OFF4]), _CMP_EQ_OQ));
#endif
	}
	for (; i != a.m_nsize; ++i) {
		if (a.m_Re[i] == b.m_Re[i])
			req[i] = 1.0;
		req[i] = 0.0;
		if (a.m_Im[i] == b.m_Im[i])
			imeq[i] = 1.0;
		imeq[i] = 0.0;
		}
}		
			
	


void
gms::math::avxvcomplex_neq(double * __restrict rneq,
			   double * __restrict imneq,
			   const AVXVComplex1D &a,
			   const AVXVComplex1D &b) {
	std::uint64_t i;

#if (GMS_MAN_PREFETCH) == 1
	size_t len = a.m_nsize;
#define CACHE_SIZE (len)
#pragma noprefetch
#endif	
	for (i = 0ULL; i != ROUND_TO_FOUR(a.m_nsize, 4ULL); i += 8ULL) {
#if (GMS_MAN_PREFETCH) == 1 

#if (CACHE_SIZE) <= (L1_MAX_DOUBLES)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
#else
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T0)
		AVXVCPROCEDURES_PREFETCH_FROM_OBJS(a, b, i, 4ULL, _MM_HINT_T1)
#endif
#endif
#if (GMS_CACHE_MEM_STORES) == 1
		// Unrolled 2x in order to match perfectly
		// per core 2 load operations per cycle.
		_mm256_store_pd(&rneq[i+OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i+OFF0]),
		_mm256_load_pd(&b.m_Re[i+OFF0]), _CMP_NEQ_OQ));
		_mm256_store_pd(&rneq[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
		_mm256_load_pd(&b.m_Re[i + OFF4]), _CMP_NEQ_OQ));

		_mm256_store_pd(&imneq[i + OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i+OFF0]),
		_mm256_load_pd(&b.m_Im[i+OFF0]),_CMP_NEQ_OQ));
		_mm256_store_pd(&imneq[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4]), _CMP_NEQ_OQ));
#else
		//	Do not pollute cache on memory stores
		//  You must be sure that stores will not be
		// used later, so do not keep them in cache.
		_mm256_stream_pd(&rneq[i + OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i + OFF0]),
		_mm256_load_pd(&b.m_Re[i + OFF0]), _CMP_NEQ_OQ));
		_mm256_stream_pd(&rneq[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Re[i + OFF4]),
		_mm256_load_pd(&b.m_Re[i + OFF4]), _CMP_NEQ_OQ));

		_mm256_stream_pd(&imneq[i + OFF0], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i + OFF0]),
		_mm256_load_pd(&b.m_Im[i + OFF0]), _CMP_NEQ_OQ));
		_mm256_stream_pd(&imneq[i + OFF4], _mm256_cmp_pd(_mm256_load_pd(&a.m_Im[i + OFF4]),
		_mm256_load_pd(&b.m_Im[i + OFF4]), _CMP_NEQ_OQ));
#endif
	}
	for (; i != a.m_nsize; ++i) {
		if (a.m_Re[i] != b.m_Re[i])
			rneq[i] = 1.0;
		else
			rneq[i] = 0.0;

		if (a.m_Im[i] != b.m_Im[i])
			imneq[i] = 1.0;
		else
			imneq[i] = 0.0;
	}
}
