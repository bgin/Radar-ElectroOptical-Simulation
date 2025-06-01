
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
#include <sys/mman.h>
#include <cerrno>
#include "GMS_sse_memset.h"
#include "GMS_sse_memcpy.h"
#include "GMS_malloc.h"

/*
   icpc -o unit_test_sse_memcpy_pd_par -fp-model fast=2 -fno-exceptions -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp unit_test_sse_memcpy_pd_par.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 GMS_config.h  GMS_sse_memcpy.h GMS_sse_memcpy.cpp unit_test_sse_memcpy_pd_par.cpp

*/

void unit_test_sse_memcpy_unroll8x_pd_par();

void unit_test_sse_memcpy_unroll8x_pd_par()
{
     using namespace gms::common;
     constexpr std::size_t buf_len{100000ULL};
     constexpr std::size_t n_threads{12ULL};
     constexpr double fill{3.14159265358979323846264338328};
     std::vector<std::thread> threads;
     double * __restrict__ p_src{NULL};
     double * __restrict__ p_dst{NULL};
     std::size_t chunk_size{0ULL};
     std::size_t rem_chunk{0ULL};
     int32_t mlock_ret{-2};
     bool b_fail{false};
     printf("[UNIT-TEST --: Parallel-Memcpy]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST --: Parallel-Memcpy]: -- fill buffer of size: %llu with value=%.16f\n",buf_len,fill);
     p_src = reinterpret_cast<double*>(gms_mm_malloc((std::size_t)(buf_len*sizeof(double)),64ULL));
     p_dst = reinterpret_cast<double*>(gms_mm_malloc((std::size_t)(buf_len*sizeof(double)),64ULL));
     mlock_ret = mlock(p_src,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- mlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     mlock_ret = -2;
     mlock_ret = mlock(p_dst,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- mlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     sse_memset_unroll8x_pd(p_src,fill,buf_len);
     chunk_size = buf_len/n_threads;
     for(std::size_t __i{0ULL}; __i != n_threads; ++__i)
     {
          double * __restrict chunk_src = p_src+__i*chunk_size;
          double * __restrict chunk_dst = p_dst+__i*chunk_size;
          threads.push_back(std::thread(sse_memcpy_unroll8x_pd,chunk_dst,chunk_src,chunk_size));
          rem_chunk += chunk_size;
     }
     std::size_t remaining_chunks{buf_len-rem_chunk};
     double * __restrict rem_src{p_src+rem_chunk};
     double * __restrict rem_dst{p_dst+rem_chunk};
     //threads.push_back(std::thread(sse_memcpy_unroll16x_ps,rem_dst,rem_src,remaining_chunks));
     sse_memcpy_unroll8x_pd(rem_dst,rem_src,remaining_chunks);
     for(std::size_t __i{0ULL}; __i != n_threads; ++__i)
     {
          threads[__i].join();
     }

     for(std::size_t __i{0ULL}; __i != buf_len; ++__i)
     {
          const double s{p_src[__i]};
          const double d{p_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST --:Parallel-Memcpy]: FAILED, found value of %.16f at pos: %llu, but expected: %.16f" ANSI_RESET_ALL "\n",d,__i,s);
               b_fail = true;
               break;
          }
     }
     if(b_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST --:Parallel-Memcpy]: PASSED!!" ANSI_RESET_ALL "\n");}
    
     printf("[UNIT-TEST --:Parallel-Memcpy]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
     mlock_ret = -2;
     mlock_ret = munlock(p_src,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- munlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     mlock_ret = -2;
     mlock_ret = munlock(p_dst,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- munlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     if(p_src) gms_mm_free(p_src);
     if(p_dst) gms_mm_free(p_dst);

}

void unit_test_sse_memcpy_unroll16x_pd_par();

void unit_test_sse_memcpy_unroll16x_pd_par()
{
     using namespace gms::common;
     constexpr std::size_t buf_len{100000ULL};
     constexpr std::size_t n_threads{12ULL};
     constexpr double fill{3.14159265358979323846264338328};
     std::vector<std::thread> threads;
     double * __restrict__ p_src{NULL};
     double * __restrict__ p_dst{NULL};
     std::size_t chunk_size{0ULL};
     std::size_t rem_chunk{0ULL};
     int32_t mlock_ret{-2};
     bool b_fail{false};
     printf("[UNIT-TEST --: Parallel-Memcpy]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST --: Parallel-Memcpy]: -- fill buffer of size: %llu with value=%.16f\n",buf_len,fill);
     p_src = reinterpret_cast<double*>(gms_mm_malloc((std::size_t)(buf_len*sizeof(double)),64ULL));
     p_dst = reinterpret_cast<double*>(gms_mm_malloc((std::size_t)(buf_len*sizeof(double)),64ULL));
     mlock_ret = mlock(p_src,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- mlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     mlock_ret = -2;
     mlock_ret = mlock(p_dst,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- mlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     sse_memset_unroll16x_pd(p_src,fill,buf_len);
     chunk_size = buf_len/n_threads;
     for(std::size_t __i{0ULL}; __i != n_threads; ++__i)
     {
          double * __restrict chunk_src = p_src+__i*chunk_size;
          double * __restrict chunk_dst = p_dst+__i*chunk_size;
          threads.push_back(std::thread(sse_memcpy_unroll16x_pd,chunk_dst,chunk_src,chunk_size));
          rem_chunk += chunk_size;
     }
     std::size_t remaining_chunks{buf_len-rem_chunk};
     double * __restrict rem_src{p_src+rem_chunk};
     double * __restrict rem_dst{p_dst+rem_chunk};
     //threads.push_back(std::thread(sse_memcpy_unroll16x_ps,rem_dst,rem_src,remaining_chunks));
     sse_memcpy_unroll16x_pd(rem_dst,rem_src,remaining_chunks);
     for(std::size_t __i{0ULL}; __i != n_threads; ++__i)
     {
          threads[__i].join();
     }

     for(std::size_t __i{0ULL}; __i != buf_len; ++__i)
     {
          const double s{p_src[__i]};
          const double d{p_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST --:Parallel-Memcpy]: FAILED, found value of %.16f at pos: %llu, but expected: %.16f" ANSI_RESET_ALL "\n",d,__i,s);
               b_fail = true;
               break;
          }
     }
     if(b_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST --:Parallel-Memcpy]: PASSED!!" ANSI_RESET_ALL "\n");}
    
     printf("[UNIT-TEST --:Parallel-Memcpy]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
     mlock_ret = -2;
     mlock_ret = munlock(p_src,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- munlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     mlock_ret = -2;
     mlock_ret = munlock(p_dst,buf_len);
     if(mlock_ret!=0) 
     {
         std::printf("!!WARNING!!: -- munlock() failed with: %d\n",mlock_ret);
         std::perror("errno:");
     }
     if(p_src) gms_mm_free(p_src);
     if(p_dst) gms_mm_free(p_dst);

}


int main()
{
     unit_test_sse_memcpy_unroll8x_pd_par();
     unit_test_sse_memcpy_unroll16x_pd_par();
     return 0;
}