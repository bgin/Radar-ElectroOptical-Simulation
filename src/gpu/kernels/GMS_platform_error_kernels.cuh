

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


extern "C"
void platform_pos_err_cuda(const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           const float * __restrict__,
                           float * __restrict__ , //min, azimuth error
                           float * __restrict__ ,  //min, elevation error
                           float * __restrict__,    //m,   range error
                           const uint32_t,
                           const uint32_t,
                           const uint32_t);


extern "C"
void elev_mpath_err_cuda(   const float * __restrict__,
                            const float * __restrict__,
                            const float * __restrict__,
                            const float,
                            const float,
                            const float,
                            const float * __restrict__,
                            const float * __restrict__,
                            const float * __restrict__,
                            float * __restrict__, // deg, elevation multipath error
                            const uint32_t,
                            const uint32_t,
                            const uint32_t);


extern "C"
void elev_refract_err_cuda(     const float * __restrict__,
                                const float * __restrict__,
                                const float * __restrict__,
                                const float,
                                const float,
                                const float,
                                const float * __restrict__,
                                const float * __restrict__,
                                const float * __restrict__, 
                                const float,
                                const int32_t,
                                float * __restrict__, //deg, rms error of elevation measurement
                                const uint32_t,
                                const uint32_t,
                                const uint32_t);














#endif /*__GMS_PLATFORM_ERROR_KERNELS_CUH__*/
