
#include <malloc.h>
#include "GMS_error_macros.h"
#include "GMS_malloc.h"
#include "GMS_config.h"

#if defined _WIN64

double *
gms::common::
gms_dmallocu(_In_ const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
	_ASSERT(len > 0ULL);
#endif
	return (reinterpret_cast<double*>(malloc(len*sizeof(double))));
}

float *
gms::common::
gms_fmallocu(_In_ const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
	_ASSERT(len > 0ULL);
#endif
	return (reinterpret_cast<float*>(malloc(len*sizeof(float))));
}

int32_t *
gms::common::
gms_imallocu(_In_ const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
	_ASSERT(len > 0ULL);
#endif
	return (reinterpret_cast<int32_t*>(malloc(len*sizeof(int32_t))));
}

double * 
gms::common::
gms_dmalloca(_In_ const std::size_t len,
			 _In_ const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
	_ASSERT(len > 0ULL);
#endif
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<double*>(_mm_malloc(len*sizeof(double) + padding64B,alignment)));
#else
	return (reinterpret_cast<double*>(_mm_malloc(len*sizeof(double),alignment)));
#endif
}

float * 
gms::common::
gms_fmalloca(_In_ const std::size_t len,
			 _In_ const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
	_ASSERT(len > 0ULL);
#endif
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<float*>(_mm_malloc(len*sizeof(float) + padding64B,alignment)));
#else
	return (reinterpret_cast<float*>(_mm_malloc(len*sizeof(float),alignment)));
#endif
}

int32_t * 
gms::common::
gms_imalloca(_In_ const std::size_t len,
	     _In_ const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
	_ASSERT(len > 0ULL);
#endif
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t) + padding64B,alignment)));
#else
	return (reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t),alignment)));
#endif
}

// ----------------- IMPORTANT ------------------------------
// ALLOCA can not be used that way!! It must be used to allocate
// dynamically on stack only inside computational functions.!!
// ---------------------------------------------------------------
/*double *
gms::common::
gms_dalloca_u(_In_ const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
      _ASSERTE(len > 0ULL);
#endif
      return (reinterpret_cast<double*>(_alloca(len*sizeof(double)));// aligned on 16-bytes
}

float *
gms::common::
gms_falloca_u(_In_ const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
       _ASSERTE(len > 0ULL);
#endif
      return (reinterpret_cast<float*>(_alloca(len*sizeof(float))));
}

int32_t *
gms::common::
gms_ialloca_u(_In_ const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
      _ASSERTE(len > 0ULL)
#endif
      return (reinterpret_cast<int32_t*>(_alloca(len*sizeof(int32_t))));
}

double *
gms::common::
gms_dalloca_a(_In_ const std::size_t len,
              _In_ const std::size_t align) { // align must have a value of  32 or 64
#if (GMS_DEBUG_ON) == 1
      _ASSERTE(len > 0ULL);
#endif
    const std::size_t nbytes = len * sizeof(double);
    void * p = _alloca(nbytes+align-1);
    double * dptr = (double*)(((UINT_PTR)p+(align-1)) & ~(align-1));
    return (dptr);
}

float *
gms::common::
gms_falloca_a(_In_ const std::size_t len,
              _In_ const std::size_t align) {
#if (GMS_DEBUG_ON) == 1
      _ASSERTE(len > 0ULL);
#endif
    const std::size_t nbytes = len * sizeof(float);
    void * p = _alloca(nbytes+align-1);
    float * fptr = (float*)(((UINT_PTR)p + (align-1)) & ~(align-1));
    return (fptr);
}*/

double * 
gms::common::
gms_edmalloca(_In_ const std::size_t len, 
			   _In_ const int32_t alignment) {
	typedef double * __restrict __declspec(align_value(64)) aligned64_r8ptr;
	aligned64_r8ptr real_ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	 real_ptr = reinterpret_cast<double*>(_mm_malloc(len*sizeof(double) + padding64B,alignment));
#else
	 real_ptr = reinterpret_cast<double*>(_mm_malloc(len*sizeof(double), alignment));
#endif
	
	if (NULL == real_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
	    ABORT_ON_ERROR("gms_edmalloca -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
	}
	return (real_ptr);
}

float * 
gms::common::
gms_efmalloca(_In_ const std::size_t len,
			  _In_ const int32_t alignment) {
	typedef float * __restrict __declspec(align_value(64)) aligned64_r4ptr;
	aligned64_r4ptr real_ptr = NULL; 
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	real_ptr = reinterpret_cast<float*>(_mm_malloc(len*sizeof(float) + padding64B,alignment));
#else
	real_ptr = reinterpret_cast<float*>(_mm_malloc(len*sizeof(float), alignment));
#endif
	if (NULL == real_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
		ABORT_ON_ERROR("gms_efmalloca -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
	return (real_ptr);
}

int32_t * 
gms::common::
gms_eimalloca4(_In_ const std::size_t len,
			   _In_ const int32_t alignment) {
	typedef int32_t * __restrict __declspec(align_value(64)) aligned64_i4ptr;
	aligned64_i4ptr i4_ptr = NULL; 
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	i4_ptr = 	reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t) + padding64B,alignment));
#else
	i4_ptr =    reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t), alignment));
#endif
	if (NULL == i4_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
			ABORT_ON_ERROR("gms_eimalloca4 -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
	return (i4_ptr);
}

int64_t * 
gms::common::
gms_eimalloca(_In_ const std::size_t len,
			  _In_ const int32_t alignment) {
	typedef int64_t * __restrict __declspec(align_value(64)) aligned64_i8ptr;
	aligned64_i8ptr i8_ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	 i8_ptr = reinterpret_cast<int64_t*>(_mm_malloc(len*sizeof(int64_t) + padding64B,alignment)); // padding for vectorization of loop peeling of remainder
#else
	 i8_ptr = reinterpret_cast<int64_t*>(_mm_malloc(len*sizeof(int64_t), alignment));
#endif
	if (NULL == i8_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
		ABORT_ON_ERROR("gms_eimalloca -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
	return (i8_ptr);
}

float *
gms::common
::gms_efmalloca_padded2D(_In_ const int32_t dim1,
					    _In_ const int32_t dim2,
						_In_ const int32_t alignment,
					    _Inout_ int32_t & pad) { // pad will contain either 0 or padding, use it later in flat array indexing (i+pad*j)
	
	_ASSERT(0 == pad && 64 == alignment);
	typedef float * __restrict __declspec(align_value(64)) align64_r4ptr;
	align64_r4ptr ptr = NULL;
	const int32_t misalign = dim2 % 16;
	if (misalign) {
		pad = dim2 + (misalign == 0 ? 0 : 16 - misalign);
		ptr = reinterpret_cast<float*>(_mm_malloc(sizeof(float)*dim1*dim2*pad,alignment));
		if (NULL == ptr &&  (dim1*dim2*pad) != 0) { goto failed;}
		return (ptr);
	}	
	else {
		ptr = reinterpret_cast<float*>(_mm_malloc(sizeof(float)*dim1*dim2, alignment));
		if (NULL == ptr && ((dim1*dim2) != 0)) { goto failed;}
		return (ptr);
	}
failed: {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	DUMP_CALLSTACK_ON_ERROR
#endif
		ABORT_ON_ERROR("gms_efmalloca_padded2D -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
     }
}	

double * 
gms::common::
gms_edmalloca_padded2D(_In_ const int32_t dim1,
					  _In_ const int32_t dim2,
					  _In_ const int32_t alignment,
					  _Inout_ int32_t & pad) { // pad will contain either 0 or padding, use it later in flat array indexing (i+pad*j)
	_ASSERT(0 == pad && 64 == alignment);
	typedef double * __restrict __declspec(aligned_value(64)) align64_r8ptr;
	align64_r8ptr ptr = NULL;
	const int32_t misalign = dim2 % 8;
	if (misalign) {
		pad = dim2 + (misalign == 0 ? 0 : 8 - misalign);
		ptr = reinterpret_cast<double*>(_mm_malloc(sizeof(double)*dim1*dim2*pad,alignment));
		if (NULL == ptr && (dim1*dim2*pad) != 0) {goto failed;}
		return (ptr);
	}
	else {
		ptr = reinterpret_cast<double*>(_mm_malloc(sizeof(double)*dim1*dim2, alignment));
		if (NULL == ptr && (dim1*dim2) != 0) { goto failed;}
		return (ptr);
	}

failed:   {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	DUMP_CALLSTACK_ON_ERROR
#endif
		ABORT_ON_ERROR("gms_edmalloca_padded2D -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
}

int32_t *
gms::common
::gms_eimalloca4_padded2D(_In_ const int32_t dim1,
						  _In_ const int32_t dim2,
					      _In_ const int32_t alignment,
						  _Inout_ int32_t & pad) {
	_ASSERT(0 == pad && 64 == alignment);
	typedef int32_t * __restrict __declspec(aligned_value(64)) align64_i4ptr;
	align64_i4ptr ptr = NULL;
	const int32_t misalign = dim2 % 16;
	if (misalign) {
		pad = dim2 + (misalign == 0 ? 0 : 16 - misalign);
		ptr = reinterpret_cast<int32_t*>(_mm_malloc(sizeof(int32_t)*dim1*dim2*pad,alignment));
		if (NULL == ptr && (dim1*dim2*pad) != 0) { goto failed;}
		return (ptr);
	}
	else {
		ptr = reinterpret_cast<int32_t*>(_mm_malloc(sizeof(int32_t)*dim1*dim2,alignment));
		if (NULL == ptr && (dim1*dim2) != 0) { goto failed; }
		return (ptr);
	}
failed:  {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	DUMP_CALLSTACK_ON_ERROR
#endif
		ABORT_ON_ERROR("gms_eimalloca4_padded2D -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
}

#elif defined __linux
#include <assert>
#include "GMS_avxvecf32.h"
#include "GMS_avx512vecf32.h"

double *
gms::common::
gms_dmallocu(const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
	assert(len > 0ULL);
#endif
	return (reinterpret_cast<double*>(malloc(len*sizeof(double))));
}

float *
gms::common::
gms_fmallocu(const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
	assert(len > 0ULL);
#endif
	return (reinterpret_cast<float*>(malloc(len*sizeof(float))));
}     

int32_t *
gms::common::
gms_imallocu(const std::size_t len) {
#if (GMS_DEBUG_ON) == 1
	assert(len > 0ULL);
#endif
	return (reinterpret_cast<int32_t*>(malloc(len*sizeof(int32_t))));
}

double * 
gms::common::
gms_dmalloca(const std::size_t len,
	     const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
	assert(len > 0ULL);
#endif
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<double*>(_mm_malloc(len*sizeof(double) + padding64B,alignment)));
#else
	return (reinterpret_cast<double*>(_mm_malloc(len*sizeof(double),alignment)));
#endif
}

float * 
gms::common::
gms_fmalloca(const std::size_t len,
	     const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
	assert(len > 0ULL);
#endif
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<float*>(_mm_malloc(len*sizeof(float) + padding64B,alignment)));
#else
	return (reinterpret_cast<float*>(_mm_malloc(len*sizeof(float),alignment)));
#endif
}

int32_t * 
gms::common::
gms_imalloca(const std::size_t len,
	     const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
	assert(len > 0ULL);
#endif
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t) + padding64B,alignment)));
#else
	return (reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t),alignment)));
#endif
}

AVXVec8 * 
gms::common::
gms_avxvec8_malloca(const std::size_t len,
		    const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
      assert(len > 0ULL);
#endif
     return (reinterpret_cast<AVXVec8*>(_mm_malloc(len*sizeof(AVXVec8),alignment)));
}

AVX512Vec16 *
gms::common::
gms_avx512vec16_malloca(const std::size_t len,
			const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
      assert(len > 0ULL);
#endif
     return (reinterpret_cast<AVX512Vec16*>(_mm_malloc(len*sizeof(AVX512Vec16),alignment)));
}

double * 
gms::common::
gms_edmalloca(const std::size_t len, 
	      const int32_t alignment) {
  typedef double * __restrict __attribute__((align(64))) aligned64_r8ptr;
	aligned64_r8ptr real_ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	 real_ptr = reinterpret_cast<double*>(_mm_malloc(len*sizeof(double) + padding64B,alignment));
#else
	 real_ptr = reinterpret_cast<double*>(_mm_malloc(len*sizeof(double), alignment));
#endif
	
	if (NULL == real_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!";
#endif
	    ABORT_ON_ERROR("gms_edmalloca -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
	}
	return (real_ptr);
}

float * 
gms::common::
gms_efmalloca(const std::size_t len,
	      const int32_t alignment) {
  typedef float * __restrict __attribute__((align(64)))  aligned64_r4ptr;
	aligned64_r4ptr real_ptr = NULL; 
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	real_ptr = reinterpret_cast<float*>(_mm_malloc(len*sizeof(float) + padding64B,alignment));
#else
	real_ptr = reinterpret_cast<float*>(_mm_malloc(len*sizeof(float), alignment));
#endif
	if (NULL == real_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!"
#endif
		ABORT_ON_ERROR("gms_efmalloca -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
	return (real_ptr);
}

int32_t * 
gms::common::
gms_eimalloca4(const std::size_t len,
	       const int32_t alignment) {
  typedef int32_t * __restrict __attribute__((align(64))) aligned64_i4ptr;
	aligned64_i4ptr i4_ptr = NULL; 
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	i4_ptr = 	reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t) + padding64B,alignment));
#else
	i4_ptr =    reinterpret_cast<int32_t*>(_mm_malloc(len*sizeof(int32_t), alignment));
#endif
	if (NULL == i4_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!";
#endif
			ABORT_ON_ERROR("gms_eimalloca4 -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
	return (i4_ptr);
}

int64_t * 
gms::common::
gms_eimalloca(const std::size_t len,
	      const int32_t alignment) {
  typedef int64_t * __restrict __attribute__((align(64))) aligned64_i8ptr;
	aligned64_i8ptr i8_ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	 i8_ptr = reinterpret_cast<int64_t*>(_mm_malloc(len*sizeof(int64_t) + padding64B,alignment)); // padding for vectorization of loop peeling of remainder
#else
	 i8_ptr = reinterpret_cast<int64_t*>(_mm_malloc(len*sizeof(int64_t), alignment));
#endif
	if (NULL == i8_ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!";
#endif
		ABORT_ON_ERROR("gms_eimalloca -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
	return (i8_ptr);
}

AVXVec8 *
gms::common::
gms_avxvec8_emalloca(const std::size_t len,
                     const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
     assert(len > 0ULL);
#endif
     AVXVec8 * ptr = NULL;
     ptr = reinterpret_cast<AVXVec8*>(_mm_malloc(len*sizeof(AVXVec8),alignment));
     if(NULL == ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
         std::cerr << "Not implemented yet!!" << std::endl;
#endif
         ABORT_ON_ERROR("gms_avxvec8_emalloca -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
     }
      return (ptr);
}

AVX512Vec16 *
gms::common::
gms_avx512vec16_emalloca(const std::size_t len,
                         const int32_t alignment) {
#if (GMS_DEBUG_ON) == 1
     assert(len > 0ULL);
#endif
     AVX512Vec16 * ptr = NULL;
     ptr = reinterpret_cast<AVX512Vec16*>(_mm_malloc(len*sizeof(AVX512Vec16),alignment));
     if(NULL == ptr && len != 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
         std::cerr << "Not implemented yet!!" << std::endl;
#endif
         ABORT_ON_ERROR("gms_avx512vec16_emalloca -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
     }
      return (ptr);
}

float *
gms::common
::gms_efmalloca_padded2D(const int32_t dim1,
			 const int32_t dim2,
			 const int32_t alignment,
			 int32_t & pad) { // pad will contain either 0 or padding, use it later in flat array indexing (i+pad*j)
	
	asssert(0 == pad && 64 == alignment);
	typedef float * __restrict __attribute__((align(64))) align64_r4ptr;
	align64_r4ptr ptr = NULL;
	const int32_t misalign = dim2 % 16;
	if (misalign) {
		pad = dim2 + (misalign == 0 ? 0 : 16 - misalign);
		ptr = reinterpret_cast<float*>(_mm_malloc(sizeof(float)*dim1*dim2*pad,alignment));
		if (NULL == ptr &&  (dim1*dim2*pad) != 0) { goto failed;}
		return (ptr);
	}	
	else {
		ptr = reinterpret_cast<float*>(_mm_malloc(sizeof(float)*dim1*dim2, alignment));
		if (NULL == ptr && ((dim1*dim2) != 0)) { goto failed;}
		return (ptr);
	}
failed: {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!"
#endif
		ABORT_ON_ERROR("gms_efmalloca_padded2D -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
     }
}


double * 
gms::common::
gms_edmalloca_padded2D(const int32_t dim1,
		       const int32_t dim2,
		       const int32_t alignment,
		       int32_t & pad) { // pad will contain either 0 or padding, use it later in flat array indexing (i+pad*j)
	assert(0 == pad && 64 == alignment);
	typedef double * __restrict __attribute__((align(64))) align64_r8ptr;
	align64_r8ptr ptr = NULL;
	const int32_t misalign = dim2 % 8;
	if (misalign) {
		pad = dim2 + (misalign == 0 ? 0 : 8 - misalign);
		ptr = reinterpret_cast<double*>(_mm_malloc(sizeof(double)*dim1*dim2*pad,alignment));
		if (NULL == ptr && (dim1*dim2*pad) != 0) {goto failed;}
		return (ptr);
	}
	else {
		ptr = reinterpret_cast<double*>(_mm_malloc(sizeof(double)*dim1*dim2, alignment));
		if (NULL == ptr && (dim1*dim2) != 0) { goto failed;}
		return (ptr);
	}

failed:   {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!"
#endif
		ABORT_ON_ERROR("gms_edmalloca_padded2D -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
}

int32_t *
gms::common
::gms_eimalloca4_padded2D(const int32_t dim1,
			  const int32_t dim2,
			  const int32_t alignment,
			  int32_t & pad) {
	assert(0 == pad && 64 == alignment);
	typedef int32_t * __restrict __attribute__((align(64))) align64_i4ptr;
	align64_i4ptr ptr = NULL;
	const int32_t misalign = dim2 % 16;
	if (misalign) {
		pad = dim2 + (misalign == 0 ? 0 : 16 - misalign);
		ptr = reinterpret_cast<int32_t*>(_mm_malloc(sizeof(int32_t)*dim1*dim2*pad,alignment));
		if (NULL == ptr && (dim1*dim2*pad) != 0) { goto failed;}
		return (ptr);
	}
	else {
		ptr = reinterpret_cast<int32_t*>(_mm_malloc(sizeof(int32_t)*dim1*dim2,alignment));
		if (NULL == ptr && (dim1*dim2) != 0) { goto failed; }
		return (ptr);
	}
failed:  {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	  std::cerr << " Not implemented yet!!"
#endif
		ABORT_ON_ERROR("gms_eimalloca4_padded2D -- !!! Memory Allocation Failure !!! ", MALLOC_FAILED)
	}
}

#endif



