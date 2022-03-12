

#ifndef __GMS_RCS_CUH__
#define __GMS_RCS_CUH__

#include <cstdint>

extern "C"
void empirical_K_cuda(const float * __restrict__,
                      const float * __restrict__,
                      const float,
                      float * __restrict__,
                      const uint32_t,
                      const uint32_t,
                      const uint32_t);

extern "C"
void effective_rcs_cuda(const float 
                        const float * __restrict__, 
                        const float, 
                        const float, 
                        float * __restrict__,
                        float * __restrict__,
                        const uint32_t,
                        const uint32_t,
                        const uint32_t);


extern "C"
void bistatic_target_rcs_cuda(const float, 
		              const float, 
	                      const float * __restrict__,
		              const float * __restrict__,
			      const float * __restrict__,
			      const float,   
			      const int32_t, 
			      float * __restrict,               
                              float * __restrict,
                              const uint32_t,
                              const uint32_t,
                              const uint32_t);


extern "C"
void antenna_rcs_cuda(   const float * __restrict__,
                         const float,            
                         const float,                
                         float * __restrict,   
                         float * __restrict__ ,   
                         const uint32_t,
                         const uint32_t,
                         const uint32_t);


extern "C"
void bird_insect_rcs_cuda(const float * __restrict__,
                          float * __restrict__,
                          const uint32_t);


extern "C"
void cone_ogive_rcs_cuda(const float, 
                         const float * __restrict__,
                         float * __restrict__,
                         const uint32_t,
                         const uint32_t,
                         const uint32_t);

extern "C"
void cylinder_rcs_cuda(const float * __restrict__,
                       const float * __restrict__,
                       const float,              
                       const float,            
                       float * __restrict__,
                       const uint32_t,
                       const uint32_t,
                       const uint32_t);





#endif /*__GMS_RCS_CUH__*/
