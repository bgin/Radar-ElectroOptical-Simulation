
#ifndef __GMS_INIT_GPUS_CUH__
#define __GMS_INIT_GPUS_CUH__



#include "cuda_runtime.h"
#include "device_launch_parameters.h"


//
//	Declarations
//

//
//	Print device capabilities to screen.
//
void devcaps_to_screen(const int, int * , const bool);


//
//	Print device capabilities to file.
//
void devcaps_to_file(const int, const char *, 
				   int * , const bool );


//
// Run integer and real arithmentic test.
// Host code.
//
cudaError_t  gpu_vec_add_tests(const int);

// 
// Device code.
// simple vector addition int32_t,REAL4 and REAL8.
//
__global__ void kvec_add_int32(int32_t * __restrict,
							const int32_t * __restrict,
							const int32_t * __restrict);

__global__ void kvec_add_real4(float * __restrict,
						   const REAL4 * __restrict ,
						   const REAL4 * __restrict);

__global__ void kvec_add_real8(double * __restrict,
						    const REAL8 * __restrict,
							const REAL8 * __restrict);




#endif /*__GMS_INIT_GPUS_CUH__*/
