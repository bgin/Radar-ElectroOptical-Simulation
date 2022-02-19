

#ifndef __GMS_RCS_CUH__
#define __GMS_RCS_CUH__

#include <cstdint>

extern "C"
void empirical_K_cuda(const float * __restrict,
                      const float * __restrict,
                      const float,
                      float * __restrict,
                      const uint32_t);





#endif /*__GMS_RCS_CUH__*/
