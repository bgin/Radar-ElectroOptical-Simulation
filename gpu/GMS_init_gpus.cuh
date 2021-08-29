
#ifndef __GMS_INIT_GPUS_CUH__
#define __GMS_INIT_GPUS_CUH__


#include <stdint.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


//
//	Declarations
//

//
//	Print device capabilities to screen.
//
void devcaps_to_screen(const int32_t, 
                       int32_t * , const bool);


//
//	Print device capabilities to file.
//
void devcaps_to_file(const int32_t, const char *, 
		     int32_t * , const bool );


//
// Run integer and real arithmentic test.
// Host code.
//
cudaError_t  gpu_vec_add_tests(const32_t int);

// 
// Device code.
// simple vector addition.
//
__global__ void kvec_add_int32(int32_t * __restrict,
			       const int32_t * __restrict,
			       const int32_t * __restrict);

__global__ void kvec_add_float(float * __restrict,
			       const float * __restrict ,
			       const float * __restrict);

__global__ void kvec_add_double(double * __restrict,
				const double * __restrict,
				const double * __restrict);




#endif /*__GMS_INIT_GPUS_CUH__*/
