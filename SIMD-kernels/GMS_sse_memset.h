
#ifndef __GMS_SSE_MEMSET_H__
#define __GMS_SSE_MEMSET_H__

namespace  file_info
{

    const unsigned int GMS_SSE_MEMSET_MAJOR = 1U;
    const unsigned int GMS_SSE_MEMSET_MINOR = 0U;
    const unsigned int GMS_SSE_MEMSET_MICRO = 0U;
    const unsigned int GMS_SSE_MEMSET_FULLVER =
      1000U*GMS_SSE_MEMSET_MAJOR+100U*GMS_SSE_MEMSET_MINOR+10U*GMS_SSE_MEMSET_MICRO;
    const char * const GMS_SSE_MEMSET_CREATE_DATE = "24-05-2025 06:08 +00200 (SAT 24 MAY 2025 GMT+2)";
    const char * const GMS_SSE_MEMSET_BUILD_DATE  = __DATE__ ":" __TIME__;
    const char * const GMS_SSE_MEMSET_AUTHOR      =  "Programmer: Bernard Gingold, e-mail: beniekg@gmail.com";
    const char * const GMS_SSE_MEMSET_DESCRIPT    =  "SSE memset kernels.";
}

#include <cstdint>
#include "GMS_config.h"

namespace gms
{
namespace common
{


__ATTR_HOT__ 
__ATTR_ALIGN__(32)
void sse_memset_unroll8x_ps(float *,const float,std::size_t);

__ATTR_HOT__
__ATTR_ALIGN__(32)
void sse_memset_unroll8x_pd(double *,const double,std::size_t);




} // common
} // gms



































#endif /*__GMS_SSE_MEMSET_H__*/