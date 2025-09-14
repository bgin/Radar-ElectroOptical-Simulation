

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!*/

#ifndef __GMS_AVX_MEMCPY_H__
#define __GMS_AVX_MEMCPY_H__ 030620250652

namespace  file_info
{

    const unsigned int GMS_AVX_MEMCPY_MAJOR = 1U;
    const unsigned int GMS_AVX_MEMCPY_MINOR = 0U;
    const unsigned int GMS_AVX_MEMCPY_MICRO = 0U;
    const unsigned int GMS_AVX_MEMCPY_FULLVER =
      1000U*GMS_AVX_MEMCPY_MAJOR+100U*GMS_AVX_MEMCPY_MINOR+10U*GMS_AVX_MEMCPY_MICRO;
    const char * const GMS_AVX_MEMCPY_CREATE_DATE = "03-06-2025 06:52 +00200 (TUE 03 JUN 2025 GMT+2)";
    const char * const GMS_AVX_MEMCPY_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const GMS_AVX_MEMCPY_AUTHOR      =  "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";
    const char * const GMS_AVX_MEMCPY_DESCRIPT    =  "AVX memcpy kernels.";
}

#include <cstdint>
#include "GMS_config.h"

#if !defined(AVX_MEMCPY_INTERLEAVE_SIMD_OPS)
#define AVX_MEMCPY_INTERLEAVE_SIMD_OPS 1
#endif 

#if !defined(AVX_MEMCPY_SOFT_PREFETCHING)
#define AVX_MEMCPY_SOFT_PREFETCHING 1
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
void avx_memcpy_unroll8x_ps(float * __restrict__,float * __restrict__,std::size_t);

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__ 
__ATTR_ALIGN__(32)
void avx_memcpy_unroll16x_ps(float * __restrict__,float * __restrict__,std::size_t);

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void avx_memcpy_unroll8x_pd(double * __restrict__,double * __restrict__,std::size_t);

#if defined (__GNUC__) && (!defined(__INTEL_COMPILER) || !defined(__ICC))
__attribute__((no_sanitize("coverage")))
#endif
__ATTR_HOT__
__ATTR_ALIGN__(32)
void avx_memcpy_unroll16x_pd(double * __restrict__,double * __restrict__,std::size_t);


} // common
} // gms















#endif /*__GMS_AVX_MEMCPY_H__*/