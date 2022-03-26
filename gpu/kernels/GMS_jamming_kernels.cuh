

#ifndef __GMS_JAMMING_KERNELS_CUH__
#define __GMS_JAMMING_KERNELS_CUH__

#include <cstdint>



extern "C"
void therm_noise_range_cuda(       const float,
                                   const float,
                                   const float,
                                   const float,
                                   const float,
                                   const float,
				   const float,
				   const float * __restrict,
				   const float,
				   const float * __restrict,
				   const float * __restrict,
				   const float,
				   const float * __restrict,
				   const float,
				   const float,
				   const float,
				   const float,
				   const float,
				   const float * __restrict,
				   const float * __restrict,
				   float * __restrict,
				   const uint32_t,
                                   const uint32_t,
                                   const uint32_t );


extern "C"
void tropo_range_loss_cuda(        const float,
                                   const float,
                                   const float,
                                   const float,
                                   const float,
				   const float,
				   const float * __restrict,
				   const float * __restrict,
				   const float,
				   const float * __restrict,
				   const float * __restrict,
				   const float,
				   const float * __restrict,
				   const float,
				   const float,
				   const float,
				   const float,
				   const float,
				   const float * __restrict,
				   const float * __restrict,
				   float * __restrict,
				   float * __restrict,
				   const uint32_t,
                                   const uint32_t,
                                   const uint32_t);
				  











#endif /*__GMS_JAMMING_KERNELS_CUH__*/
