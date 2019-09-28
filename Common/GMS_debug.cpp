
#include "GMS_debug.h"
#include "GMS_error_macros.h"


#if (GMS_DEBUG_ON) == 1
     #if defined (_WIN64)

double *
gms::common
::gms_dmallocu_dbg(_In_ const size_t Size,
		   _In_ int32_t BlockType,
		   _In_ const char * filename,
		   _In_ int32_t linenum) {
	_ASSERTE(Size > 0ULL);
	return (reinterpret_cast<double*>(_malloc_dbg(Size * sizeof(double),BlockType,filename,linenum)));
}

float * 
gms::common
::gms_fmallocu_dbg(_In_ const size_t Size,
				   _In_ int32_t BlockType,
				   _In_ const char * filename,
				   _In_ int32_t linenum) {
	_ASSERTE(Size > 0ULL);
	return (reinterpret_cast<float*>(_malloc_dbg(Size * sizeof(float),BlockType,filename,linenum)));
}

int32_t *
gms::common
::gms_imallocu_dbg(_In_ const size_t Size,
				   _In_ int32_t BlockType,
				   _In_ const char * filename,
				   _In_ int32_t linenum) {
	_ASSERTE(Size > 0ULL);
	return (reinterpret_cast<int32_t*>(_malloc_dbg(Size * sizeof(int32_t),BlockType,filename,linenum)));
}

double * 
gms::common
::gms_dmalloca_dbg(_In_ const size_t Size,
				   _In_ size_t Alignment,
				   _In_ const char * filename,
				   _In_ int32_t linenum) {
	_ASSERTE(Size > 0ULL);
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<double*>(
				_aligned_malloc_dbg(Size * sizeof(double) + padding64B,Alignment,filename,linenum)));
#else
	return (reinterpret_cast<double*>(
			   _aligned_malloc_dbg(Size * sizeof(double),Alignment,filename,linenum)));
#endif
}
					 
float * 
gms::common
::gms_fmalloca_dbg(_In_ const size_t Size,
				   _In_ size_t Alignment,
				   _In_ const char * filename,
				   _In_ int32_t linenum) {
	_ASSERTE(Size > 0ULL);
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<float*>(
				_aligned_malloc_dbg(Size * sizeof(float) + padding64B, Alignment,filename,linenum)));
#else
	return (reinterpret_cast<float*>(
				_aligned_malloc_dbg(Size * sizeof(float), Alignment, filename,linenum)));
#endif
}

int32_t *
gms::common
::gms_imalloca_dbg(_In_ const size_t Size,
				   _In_ size_t Alignment,
				   _In_ const char * filename,
				   _In_ int32_t linenum) {
	_ASSERTE(Size > 0ULL);
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	return (reinterpret_cast<int32_t*>(
				 _aligned_malloc_dbg(Size * sizeof(int32_t) + padding64B, Alignment,filename,linenum)));
#else
	return (reinterpret_cast<int32_t*>(
		        _aligned_malloc_dbg(Size * sizeof(int32_t), Alignment, filename, linenum)));
#endif
}

double * 
gms::common
::gms_edmalloca_dbg(_In_ const size_t Size,
					_In_ size_t Alignment,
					_In_ const char * filename,
					_In_ int32_t linenum) {
	typedef double * __restrict __declspec(align_value(64)) r8ptra;
	r8ptra ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	ptr = reinterpret_cast<double*>(
				_aligned_malloc_dbg(Size * sizeof(double) + padding64B,Alignment,filename,linenum));
#else
	ptr = reinterpret_cast<double*>(
		       _aligned_malloc_dbg(Size * sizeof(double), Alignment, filename, linenum));
#endif
	if (NULL == ptr && 0ULL != Size) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
			ABORT_ON_ERROR("gms_edmalloca_dbg -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
	}
	return (ptr);
}

float * 
gms::common::
gms_efmalloca_dbg(_In_ const size_t Size,
				  _In_ size_t Alignment,
				  _In_ const char * filename,
				  _In_ int32_t linenum) {
	typedef float * __restrict __declspec(align_value(64)) r4ptr;
	r4ptr ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
	ptr = reinterpret_cast<float*>(
					_aligned_malloc_dbg(Size * sizeof(float) + padding64B, Alignment,filename,linenum));
#else
	ptr = reinterpret_cast<float*>(
				   _aligned_malloc_dbg(Size * sizeof(float), Alignment, filename, linenum));
#endif
	if (NULL == ptr && 0ULL != Size) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
			ABORT_ON_ERROR("gms_efmalloca_dbg -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
	}
	return (ptr);
}

int64_t * 
gms::common
::gms_eimalloca_dbg(_In_ const size_t Size,
				    _In_ size_t Alignment,
					_In_ const char * filename,
				    _In_ int32_t linenum) {
	typedef int32_t * __restrict __declspec(align_value(64)) i4ptr;
	i4ptr ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
		ptr = reinterpret_cast<int32_t*>(
					_aligned_malloc_dbg(Size * sizeof(int32_t) + padding64B,Alignment,filename,linenum));
#else
	ptr = reinterpret_cast<int32_t*>(
		_aligned_malloc_dbg(Size * sizeof(int32_t), Alignment, filename, linenum));
#endif
	if (NULL == ptr && 0ULL != Size) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
		DUMP_CALLSTACK_ON_ERROR
#endif
			ABORT_ON_ERROR("gms_eimalloca_dbg -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
	}
	return (ptr);
}

     #elif defined (__linux)

     // use mallopt with value 3

double *
gms_common
::gms_dmallocu_dbg(const size_t Size) {
      if(0ULL >= Size) return (NULL);
      return(reinterpret_cast<double*>(malloc(Size*sizeof(double))));
}

float *
gms_common
::gms_fmallocu_dbg(const size_t Size) {
      if(0ULL >= Size) return (NULL);
      return (reinterpret_cast<float*>(malloc(Size*sizeof(float))));
}

int32_t *
gms_common
::gms_imallocu_dbg(const size_t Size) {
      if(0ULL >= Size) return (NULL);
      return (reinterpret_cast<int32_t*>(malloc(Size*sizeof(int32_t))));
}

double *
gms_common
::gms_dmalloca_dbg(const size_t Size, const size_t Alignment) {
      if(0ULL >= Size) return (NULL);
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
      return (reinterpret_cast<double*>(
				 _mm_malloc(Size*sizeof(double)+padding64B,Alignment)));
#else
      return (reinterpret_cast<double*>(
				 _mm_malloc(Size*sizeof(double),Alignment)));
#endif
}

float *
gms_common
::gms_fmalloca_dbg(const size_t Size, const size_t Alignment) {
      if(0ULL >= Size) return (NULL);
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
      return (reinterpret_cast<float*>(
				       _mm_malloc(Size*sizeof(float)+padding64B,Alignment)));
#else
      return (reinterpret_cast<float*>(
				       _mm_malloc(Size*sizeof(float),Alignment)));
#endif
}

int32_t *
gms_common
::gms_imalloca_dbg(const size_t Size, const size_t Alignment) {
      if(0ULL >= Size) return (NULL);
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
      return (reinterpret_cast<int32_t*>(
					 _mm_malloc(Size*sizeof(int32_t)+padding64B,Alignment)));
#else
      return (reinterpret_cast<int32_t*>(
					 _mm_malloc(Size*sizeof(int32_t),Alignment)));
#endif
}

double *
gms_common
::gms_edmalloca_dbg(const size_t Size, const size_t Alignment) {
      if(0ULL >= Size) return (NULL);
      typedef double * __restrict __attribute__((align(64))) r8ptr;
      r8ptr ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
      ptr = reinterpret_cast<double*>(
				      _mm_malloc(Size*sizeof(double)+padding64B,Alignment));
#else
      ptr = reinterpret_cast<double*>(
				      _mm_malloc(Size*sizeof(double),Alignment));
#endif
      if(NULL == ptr && Size > 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	// will be added later
	std::cerr << "Not implemented yet!!";
#endif
	  ABORT_ON_ERROR("gms_edmalloca_dbg -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
      }
      return (ptr);
}

float *
gms_common
::gms_efmalloca_dbg(const size_t Size, const size_t Alignment) {
      if(0ULL >= Size) return (NULL);
      typedef float * __restrict __attribute__((align(64))) r4ptr;
      r4ptr ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
       ptr = reinterpret_cast<float*>(
				      _mm_malloc(Size*sizeof(float)+padding64B,Alignment));
#else
      ptr = reinterpret_cast<float*>(
				      _mm_malloc(Size*sizeof(float),Alignment));
#endif
      if(NULL == ptr && Size > 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	// will be added later
	std::cerr << "Not implemented yet!!";
#endif
	ABORT_ON_ERROR("gms_efmalloca_dbg -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
      }
      return (ptr);
}

int32_t *
gms_common
::gms_eimalloca_dbg(const size_t Size, const size_t Alignment) {
      if(0ULL >= Size) return (NULL);
      typedef int32_t * __restrict __attribute__((align(64))) i4ptr;
      i4ptr ptr = NULL;
#if (ADD_PADDING_64B_LOOP_PEEL) == 1
      ptr = reinterpret_cast<int32_t*>(
				      _mm_malloc(Size*sizeof(int32_t)+padding64B,Alignment));
#else
      ptr = reinterpret_cast<int32_t*>(
				      _mm_malloc(Size*sizeof(int32_t),Alignment));
#endif
      if(NULL == ptr && Size > 0ULL) {
#if (PRINT_CALLSTACK_ON_ERROR) == 1
	// will be added later
	std::cerr << "Not implemented yet!!";
#endif
	ABORT_ON_ERROR("gms_eimalloca_dbg -- !!! Memory Allocation Failure !!!", MALLOC_FAILED)
      }
      return (ptr);      
}

     #endif


#else
#error "*****ERROR***** -->  This file compiles only in debug build!!"
#endif
