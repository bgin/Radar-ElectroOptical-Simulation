
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "GMS_cuBLAS_L1_wrappers.cuh"
#include "GMS_gpu_config.cuh"
#if (PROFILE_HOST_TO_DEVICE) == 1
#include <immintrin.h> //rdtscp
#endif 
#include "GMS_cuda_memops.cuh"

const uint64_t rdtscp_cost = 42; // Skylake uarch
void
cuBLAS_Isamax_iface(const float * __restrict h_ptr,
                    const int32_t n,
                    const int32_t incx,
                    int32_t * result,
                    cublasStatus * __restrict err,
                    int32_t * __restrict ierr
                    uint64_t * __restrict tsc_delta) { // for profiling usage 
                                  
    if(__builtin_expect(n<=1,0)) {
       *ierr = -1;
       return;
      }
    *err    = 0;
    *ierr   = 0;
    *result = 0;
#if (PROFILE_HOST_TO_DEVICE) == 1
     uint64_t dummy1;
     uint32_t dummy2;
     uint64_t tsc_start,tsc_end;
     uint32_t coreid;
#endif
    float * __restrict d_ptr = NULL;
    cublasHandle_t handle    = NULL;
    cudaStream_t   stream    = NULL;
    cublasStatus_t status;
    GMS_CUBLAS_STAT_CHECK(cublasCreate(&handle));
    GMS_CUDA_DEBUG_CHECK(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
    GMS_CUBLAS_STAT_CHECK(cublasSetStream(handle,stream));
    alloc_float_gpu(&d_ptr[0],(size_t)n,ierr);
    GMS_CUDA_DEBUG_CHECK(cudaMemcpyAsync(d_ptr,h_ptr,sizeof(float)*((size_t)n),
                         cudaMemcpyHostToDevice,stream));
#if (PROFILE_HOST_TO_DEVICE) == 1
      dummy1    = __rdtscp(&dummy2);
      tsc_start = __rdtscp(&coreid);
#endif
     GMS_CUBLAS_STAT_CHECK(cublasIsamax(cublasH,n,d_ptr,incx,&result));
#if (PROFILE_HOST_TO_DEVICE) == 1
      tsc_end     = __rdtscp(&coreid);
      *tsc_delta  = tsc_end-tsc_start-rdtscp_cost;
      printf("cublasIsamax executed in: %llu reference cycles\n",*tsc_delta);
#endif
      GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
      GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
      GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));
      *ierr = 0;
      *err  = status;
      Error: 
            {
                GMS_CUDA_DEBUG_CHECK(cudaFree(d_ptr));
                GMS_CUBLAS_STAT_CHECK(cublasDestroy(handle));
                GMS_CUDA_DEBUG_CHECK(cudaStreamDestroy(stream));  
                ierr = -2;
                *err = status;    
      }
}

