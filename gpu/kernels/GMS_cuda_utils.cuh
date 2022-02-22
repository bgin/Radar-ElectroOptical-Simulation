
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
//===========================================================================//
// Adapted from the Cuda-samples and modified by removal of templated code.
//===========================================================================//
/* Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

// Utility class used to avoid the linker errors.
// Templating removed created two three versions for the following data types
// 1) float
// 2) int
// 3) cuComplex

struct SharedMemoryR4 {
 
   __device__ inline operator float *() {
         extern __shared__ float __smem[];
         return (float*)__smem;
  }

  __device__ inline operator const float *() const {
         extern __shared__ float __smem[];
         return (float*)__smem;
  }

};


struct SharedMemoryI4 {

   __device__ inline operator int *() {
        extern __shared__ int __smem[];
        return (int*)__smem;
  }

  __device__ inline operator const int *() const {
        extern __shared__ int __smem[];
        return (int*)__smem;
  }
};


struct SharedMemoryC4 {

    __device__ inline operator cuComplex *() {
            extern __shared__ cuComplex __smem[];
            return (cuComplex*)__smem;
    }

    __device__ inline operator const cuComplex *() const {
            extern __shared__ cuComplex __smem[];
            return (cuComplex*)__smem;
    }
};


__device__ __forceinline__
float warpReduceSum(uint32_t mask, 
                    float sum) {
  for(int32_t off = warpSize/2; off > 0; off /= 2) {
      sum += __shfl_down_sync(mask,sum,off);
  }
   return (sum);
}

  /*
    This version uses sequential addressing -- no divergence or bank conflicts.
*/

__global__
void reduce_r4_v1(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n) {
 
   // Create thread block group
   namespace cg = cooperative_groups;
   cg::thread_block tb = cg::this_thread_block();
   float * __restrict__ shmem  = SharedMemoryR4();
   
   // Shmem loading
   uint32_t tid = threadIdx.x;
   uint32_t i   = blockIdx.x*blockDIm.x+threadIdx.x;
   shmem[tid]   = (i<n) ? in[i] : 0;
   cg::sync(tb);
   
   // Reduction in shared memory
   for(uint32_t s = blockDIm.x/2; s > 0; s >>= 1) {
       if(tid<s){
          shmem[tid] += shmem[tid+s];
       }
    cg::sync(tb);
   }
   // SToring result in global memory
   if(tid == 0) out[blockIdx.x] = shmem[0];
}

__global__
void reduce_i4_v1(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n) {
 
   // Create thread block group
   namespace cg = cooperative_groups;
   cg::thread_block tb = cg::this_thread_block();
   int32_t * __restrict__ shmem  = SharedMemoryI4();
   
   // Shmem loading
   uint32_t tid = threadIdx.x;
   uint32_t i   = blockIdx.x*blockDIm.x+threadIdx.x;
   shmem[tid]   = (i<n) ? in[i] : 0;
   cg::sync(tb);
   
   // Reduction in shared memory
   for(uint32_t s = blockDIm.x/2; s > 0; s >>= 1) {
       if(tid<s){
          shmem[tid] += shmem[tid+s];
       }
    cg::sync(tb);
   }
   // SToring result in global memory
   if(tid == 0) out[blockIdx.x] = shmem[0];
}

__global__
void reduce_c4_v1(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n) {
 
   // Create thread block group
   namespace cg = cooperative_groups;
   cg::thread_block tb = cg::this_thread_block();
   cuComplex * __restrict__ shmem  = SharedMemoryC4();
   
   // Shmem loading
   uint32_t tid = threadIdx.x;
   uint32_t i   = blockIdx.x*blockDIm.x+threadIdx.x;
   shmem[tid]   = (i<n) ? in[i] : 0;
   cg::sync(tb);
   
   // Reduction in shared memory
   for(uint32_t s = blockDIm.x/2; s > 0; s >>= 1) {
       if(tid<s){
          shmem[tid] += shmem[tid+s];
       }
    cg::sync(tb);
   }
   // SToring result in global memory
   if(tid == 0) out[blockIdx.x] = shmem[0];
}


/*
   This version uses n/2 threads --
    it performs the first level of reduction when reading from global memory. 
*/

__global__
void reduce_r4_v2(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n) {
     
     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     float * __restrict__ shmem = SharedMemoryR4();
     
     // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDimx.x/2)+threadIdx.x;
     float sum    = (i<n) ? in[i] : 0.0f;
     if(i+blockDim.x<n) sum += in[i+blockDim.x];
     shmem[tid] = sum;
     cg::sync(tb); // barier before an entry to shared memory reduce ops
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const float tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::synx(tb);
     // Result is written back to global memory
     if(tid==0) out[blockIdx.x] = sum;
}


__global__
void reduce_i4_v2(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n) {
     
     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     int32_t * __restrict__ shmem = SharedMemoryI4();
     
     // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDimx.x/2)+threadIdx.x;
     int32_t sum    = (i<n) ? in[i] : 0;
     if(i+blockDim.x<n) sum += in[i+blockDim.x];
     shmem[tid] = sum;
     cg::sync(tb); // barier before an entry to shared memory reduce ops
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const int32_t tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::synx(tb);
     // Result is written back to global memory
     if(tid==0) out[blockIdx.x] = sum;
}


__global__
void reduce_c4_v2(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n) {
     
     // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     cuComplex * __restrict__ shmem = SharedMemoryC4();
     
     // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDimx.x/2)+threadIdx.x;
     cuComplex sum    = (i<n) ? in[i] : 0;
     if(i+blockDim.x<n) sum += in[i+blockDim.x];
     shmem[tid] = sum;
     cg::sync(tb); // barier before an entry to shared memory reduce ops
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const cuComplex tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::synx(tb);
     // Result is written back to global memory
     if(tid==0) out[blockIdx.x] = sum;
}


/*
    This version uses the warp shuffle operation if available to reduce
    warp synchronization. When shuffle is not available the final warp's
    worth of work is unrolled to reduce looping overhead.

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/

__global__
void reduce_r4_v3(float * __restrict__ in,
                  float * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {
  
      // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     float * __restrict__ shmem = SharedMemoryR4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDim.x/2)+threadIdx.x;
     float sum    = (i<n) ? in[i] : 0.0f;
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum;
     cg::sync(tb);
     // Reduction in shared memory
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const float tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::sync(tb);
     
     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if(blocksize >= 64) sum += shmem[tid+32];
        // Final warp reduction (shuffle)
        for(int32_t off = tile32.size()/2; off > 0; off /= 2) {
             sum += tile32.shfl_down(sum,off);
        }
     }
     // Result is written to global memory
     if(tb.thread_rank() == 0) out[blockIdx.x] = sum;
}


__global__
void reduce_i4_v3(int32_t * __restrict__ in,
                  int32_t * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {
  
      // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     int32_t * __restrict__ shmem = SharedMemoryI4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDim.x/2)+threadIdx.x;
     int32_t sum    = (i<n) ? in[i] : 0;
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum;
     cg::sync(tb);
     // Reduction in shared memory
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const int32_t tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::sync(tb);
     
     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if(blocksize >= 64) sum += shmem[tid+32];
        // Final warp reduction (shuffle)
        for(int32_t off = tile32.size()/2; off > 0; off /= 2) {
             sum += tile32.shfl_down(sum,off);
        }
     }
     // Result is written to global memory
     if(tb.thread_rank() == 0) out[blockIdx.x] = sum;
}


__global__
void reduce_c4_v3(cuComplex * __restrict__ in,
                  cuComplex * __restrict__ out,
                  const uint32_t n,
                  const uint32_t blocksize) {
  
      // Thread block group instantiation
     namespace cg = cooperative_groups;
     cg::thread_block tb = cg::this_thread_block();
     cuComplex * __restrict__ shmem = SharedMemoryC4();
      // first level of reduction
     // global memory read, write to shared memory
     uint32_t tid = threadIdx.x;
     uint32_t i   = blockIdx.x*(blockDim.x/2)+threadIdx.x;
     cuComplex sum    = (i<n) ? in[i] : {0.0f,0.0f};
     if(i+blocksize<n) sum += in[i+blocksize];
     shmem[tid] = sum;
     cg::sync(tb);
     // Reduction in shared memory
     for(uint32_t s = blockDim.x/2; s > 0; s >>= 1) {
         if(tid<s) {
             const cuComplex tmp = sum
             shmem[tid] = tmp + shmem[tid+s];
         }
     }
     cg::sync(tb);
     
     cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(tb);
     if(tb.thread_rank() < 32) {
        // Fetch final intermediate sum from 2nd warp
        if(blocksize >= 64) sum += shmem[tid+32];
        // Final warp reduction (shuffle)
        for(int32_t off = tile32.size()/2; off > 0; off /= 2) {
             sum += tile32.shfl_down(sum,off);
        }
     }
     // Result is written to global memory
     if(tb.thread_rank() == 0) out[blockIdx.x] = sum;
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
