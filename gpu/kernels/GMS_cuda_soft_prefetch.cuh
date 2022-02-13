
#ifndef __GMS_CUDA_SOFT_PREFETCH_CUH__
#define __GMS_CUDA_SOFT_PREFETCH_CUH__



#define PXL_GLOBAL_PTR   "r"

// Software prefetching.
inline __device__ 
void __prefetch_global_l1(const void* const ptr)
{
  asm("prefetch.global.L1 [%0];" : : PXL_GLOBAL_PTR(ptr));
}

inline __device__
void __prefetch_global_uniform(const void* const ptr)
{
  asm("prefetchu.L1 [%0];" : : PXL_GLOBAL_PTR(ptr));
}

inline __device__ 
void __prefetch_global_l2(const void* const ptr)
{
  asm("prefetch.global.L2 [%0];" : : PXL_GLOBAL_PTR(ptr));
}












#endif /*__GMS_CUDA_SOFT_PREFETCH_CUH__*/
