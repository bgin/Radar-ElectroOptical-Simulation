

#ifndef __GMS_SSE_MEMCPY_H__
#define __GMS_SSE_MEMCPY_H__ 

namespace  file_info
{

    const unsigned int GMS_SSE_MEMCPY_MAJOR = 1U;
    const unsigned int GMS_SSE_MEMCPY_MINOR = 0U;
    const unsigned int GMS_SSE_MEMCPY_MICRO = 0U;
    const unsigned int GMS_SSE_MEMCPY_FULLVER =
      1000U*GMS_SSE_MEMCPY_MAJOR+100U*GMS_SSE_MEMCPY_MINOR+10U*GMS_SSE_MEMCPY_MICRO;
    const char * const GMS_SSE_MEMCPY_CREATE_DATE = "26-05-2025 08:33 +00200 (MON 26 MAY 2025 GMT+2)";
    const char * const GMS_SSE_MEMCPY_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const GMS_SSE_MEMCPY_AUTHOR      =  "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";
    const char * const GMS_SSE_MEMCPY_DESCRIPT    =  "SSE memcpy kernels.";
}

#include <cstdint>
#include "GMS_config.h"

#if !defined(SSE_MEMCPY_INTERLEAVE_SIMD_OPS)
#define SSE_MEMCPY_INTERLEAVE_SIMD_OPS 1
#endif 

namespace gms
{
namespace common
{

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__ 
__ATTR_ALIGN__(32)
void sse_memcpy_unroll8x_ps(float * __restrict__,float * __restrict__,std::size_t);

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__ 
__ATTR_ALIGN__(32)
void sse_memcpy_unroll16x_ps(float * __restrict__,float * __restrict__,std::size_t);

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void sse_memcpy_unroll8x_pd(double * __restrict__,double * __restrict__,std::size_t);

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void sse_memcpy_unroll16x_pd(double * __restrict__,double * __restrict__,std::size_t);



} // common
}// gms


















#endif /*__GMS_SSE_MEMCPY_H__*/