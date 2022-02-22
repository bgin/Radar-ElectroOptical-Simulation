
#ifndef __GMS_CUDA_UTILS_CU__
#define __GMS_CUDA_UTILS_CU__

/*

     A collection of various short cuda utility
     functions both device and kernels.
     Part of them adapted from the SO answers.
*/


#include <cub/cub.cuh>
#include <cuda.h>
#include <cstdint>


//
// Block reduction kernel
// 

__global__
void
block_reduction_r4(const float * __restrict__ in,
                   float       * __restrict__ out,
                   const int32_t n) {
     
       int tid = blockIdx.x*blockDIm.x+threadIdx.x;
       typedef cub::BlockReduce<float,32> br; // block-size == 32
       __shared__ br::TempStorage ts;
       float sum = 0.0f;
       if(tid < n) sum = br(ts).Sum(in[tid]);
       if(threadIdx.x == 0) out[blockIdx.x] = sum;
}

//Adapted from CUDA Handbook, A Comprehensive Guide to GPU Programming by Nicholas Wilt

__global__
void atomic_reduction(float  * __restrict__  out,
                      const float * __restrict__ in,
                      const int32_t n) {
     const int32_t tid = threadIdx.x;
     float partsum = 0.0f;
     int32_t i = blockIdx.x*blockDIm.x+tid;
     const int32_t stride = blockDIm.x*gridDIm.x;
     for(; i < n; i += stride) {
         partsum += in[i];
    }
    atomicAdd(&out[i],partsum);
}



#define R4_BUFF_LEN  256
#define C4_BUFF_LEN  128


__global__ void
device_copy_r4(float * __restrict__ to
               const float * __restrict__ from,
               const int32_t n) {
     float buf[R4_BUFF_LEN];
     int32_t tid = R4_BUFF_LEN*blockIdx.x*blockDim.x+threadIdx.x;
     const int32_t stride = R4_BUFF_LEN*blockDIm.x*gridDIm.x;
     int32_t i;
     for(i = tid; i < n-stride i += stride) {
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             buf[j] = from[idx];
         }
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             to[idx] = buf[j];
         }
     }
     
     for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDIm.x;
             if(idx<n) buf[j] = from[idx];
         }
         for(int32_t j = 0; j < R4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             if(idx<n) to[idx] = buf[j];
         }
     }
}


__global__ void
device_copy_c4(cuComplex * __restrict__ to,
               const cuComplex * __restrict__ from,
               const int32_t n) {

     cuComplex buf[C4_BUFF_LEN];
     int32_t tid = C4_BUFF_LEN*blockIdx.x*blockDim.x+threadIdx.x;
     const int32_t stride = C4_BUFF_LEN*blockDIm.x*gridDIm.x;
     int32_t i;
     for(i = tid; i < n-stride i += stride) {
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             buf[j] = from[idx];
         }
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             to[idx] = buf[j];
         }
     }
     
     for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDIm.x;
             if(idx<n) buf[j] = from[idx];
         }
         for(int32_t j = 0; j < C4_BUFF_LEN; ++j) {
             const int32_t idx = i+j*blockDim.x;
             if(idx<n) to[idx] = buf[j];
         }
     }
}

#endif /*__GMS_CUDA_UTILS_CU__*/
