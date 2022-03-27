

#ifndef __GMS_RADAR_JAMMING_KERNELS_CUH__
#define __GMS_RADAR_JAMMING_KERNELS_CUH__

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


extern "C"
void jammer_erp_cuda(const float * __restrict__,
                     const float,
                     const float,
                     float * __restrict__,
                     const uint32_t,
                     const uint32_t,
                     const uint32_t);


extern "C"
void jammer_ernp_cuda(const float * __restrict__,
                      const float,
                      const float,
                      const float,
                      const float,
                      float * __restrict__,
                      const uint32_t,
                      const uint32_t,
                      const uint32_t);
				  

extern "C"
void jammer_spectr_dens_cuda(  const float * __restrict__,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float * __restrict__,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                float * __restrict__,
                                const uint32_t,
                                const uint32_t,
                                const uint32_t);


void single_jammer_temp_cuda(   const float * __restrict__,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float * __restrict_,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                const float,
                                float * __restrict__,
                                const uint32_t,
                                const uint32_t,
                                const uint32_t);




























#endif /*__GMS_RADAR_JAMMING_KERNELS_CUH__*/
