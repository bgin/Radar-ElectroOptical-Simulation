

#ifndef __GMS_CUBLAS_L1_WRAPPERS_CUH__
#define __GMS_CUBLAS_L1_WRAPPERS_CUH__


/*
     cuBLAS Level-1 wrappers 
*/

#include <cstdint>

#include "GMS_config.h"


void
cuBLAS_Isamax_iface(const float * __restrict,
                    const int32_t,
                    const int32_t,
                    int32_t *     __restrict,
                    cublasStatus_t * __restrict,
                    int32_t * __restrict,
                    uint64_t * __restrict) 
                                       __ATTR_HOT__
                                       __ATTR_ALIGN__(32);

void
cuBLAS_Icamax_iface(const cuComplex * __restrict,
                    const int32_t,
                    const int32_t,
                    int32_t *  __restrict,
                    cublasStatus_t * __restrict,
                    int32_t * __restrict,
                    uint64_t * __restrict)  
                                       __ATTR_HOT__
                                       __ATTR_ALIGN__(32);








#endif /*__GMS_CUBLAS_L1_WRAPPERS_CUH__*/
