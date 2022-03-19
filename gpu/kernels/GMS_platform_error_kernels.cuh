

#ifndef __GMS_PLATFORM_ERROR_KERNELS_CUH__
#define __GMS_PLATFORM_ERROR_KERNELS_CUH__

#include <cstdint>

extern "C"
void platform_orient_err_cuda(const float * __restrict__ ,
                              const float * __restrict__ ,
                              const float * __restrict__ , 
                              const float * __restrict__ ,
                              const float * __restrict__ ,
                              float * __restrict__ ,//ang,min, azimuth measurement error
                              float * __restrict__ , //ang,min, elevation measurement error
                              const uint32_t,
                              const uint32_t,
                              const uint32_t);

















#endif /*__GMS_PLATFORM_ERROR_KERNELS_CUH__*/
