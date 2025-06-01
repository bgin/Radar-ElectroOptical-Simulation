

#include <stdio.h>
#include <cstdlib>
#include <random>
#include <functional>
#include <ctime>
#include "GMS_malloc.h"
#include "GMS_sse_memset.h"

/*
   icpc -o unit_test_sse_memset_pd_rand_alloc -fp-model fast=2 -fno-exceptions -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp unit_test_sse_memset_pd_rand_alloc.cpp

   g++ -o unit_test_sse_memset_pd_rand_alloc -fp-model fast=2 -fno-exceptions -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp unit_test_sse_memset_pd_rand_alloc.cpp

   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 GMS_config.h  GMS_sse_memset.h GMS_sse_memset.cpp unit_test_sse_memset_pd_rand_alloc.cpp

*/

void unit_test_sse_memset_unroll8x_pd_rand_alloc();

void unit_test_sse_memset_unroll8x_pd_rand_alloc()
{
     using namespace gms::common;

     constexpr int32_t n_fillers{100};
     double rand_fillers[n_fillers] = {};
     double * fptrs[n_fillers] = {};
     bool fails_mask[n_fillers] = {};
     std::size_t buf_len_r{};
     const std::size_t lo_limit{501ull};
     const std::size_t hi_limit{25145ull};
     std::clock_t seed_buf{};
     std::clock_t seed_fill{};
     const int32_t n_iters{100};
     std::size_t fails_count{};
     
     seed_buf = std::clock();
     auto rbuff_len{std::bind(std::uniform_int_distribution<std::size_t>(lo_limit,hi_limit),
                          std::mt19937(seed_buf))};
     seed_fill = std::clock();
     auto rbuff_fill{std::bind(std::uniform_real_distribution<double>(0.1f,1.0f),
                            std::mt19937(seed_fill))};
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     for(int32_t __i{0}; __i != n_iters; ++__i)
     {
         double r{rbuff_fill()};
         rand_fillers[__i] = r;
         std::size_t rlen{rbuff_len()};
         if(rlen<lo_limit) rlen = lo_limit;
         else if(rlen>hi_limit) rlen = hi_limit;
         //printf("rlen=%llu\n",rlen);
         fptrs[__i] = reinterpret_cast<double*>(gms_mm_malloc((std::size_t)(sizeof(double)*rlen),64ull));
         __builtin_memset(&fptrs[__i][0],0,(unsigned long)rlen);
         const double filler{rand_fillers[__i]};
         sse_memset_unroll8x_pd(fptrs[__i],filler,rlen);
       
         for(std::size_t __j{0}; __j != rlen; ++__j)
         {
            const double curr_val{fptrs[__i][__j]};
            if(curr_val==0.0)
            {
                 printf(ANSI_COLOR_RED "[UNIT-TEST: random length of buffer]: FAILED,value of %.16f at pos: %llu,iter=%d,buf_len=%llu" ANSI_RESET_ALL "\n",curr_val,__j,__i,rlen);
                 fails_mask[__i] = true;
                 fails_count += 1ull;
            }
         }
        
     }
     if(fails_count == 0ull)
     {
         printf(ANSI_COLOR_GREEN "[UNIT-TEST: random length of buffer]: PASSED!!" ANSI_RESET_ALL "\n");
     }
     else
     {
         for(int32_t __i{0}; __i != n_iters; ++__i)
         {
            printf(ANSI_COLOR_RED "[UNIT-TEST: random length of buffer]: fails mask, iter=%d,mask=%d" ANSI_RESET_ALL "\n",__i,fails_mask[__i]);
         }
         printf(ANSI_COLOR_RED "[UNIT-TEST: random length of buffer]: fails count=%llu" ANSI_RESET_ALL "\n",fails_count);
     }
     for(int32_t __i{0}; __i != n_iters; ++__i)
     {
         if(fptrs[__i] != NULL)
         {
            gms_mm_free(fptrs[__i]);
         }
     }
     printf("[UNIT-TEST]: %s() ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}


void unit_test_sse_memset_unroll16x_pd_rand_alloc();

void unit_test_sse_memset_unroll16x_pd_rand_alloc()
{
     using namespace gms::common;

     constexpr int32_t n_fillers{100};
     double rand_fillers[n_fillers] = {};
     double * fptrs[n_fillers] = {};
     bool fails_mask[n_fillers] = {};
     std::size_t buf_len_r{};
     const std::size_t lo_limit{501ull};
     const std::size_t hi_limit{25145ull};
     std::clock_t seed_buf{};
     std::clock_t seed_fill{};
     const int32_t n_iters{100};
     std::size_t fails_count{};
     
     seed_buf = std::clock();
     auto rbuff_len{std::bind(std::uniform_int_distribution<std::size_t>(lo_limit,hi_limit),
                          std::mt19937(seed_buf))};
     seed_fill = std::clock();
     auto rbuff_fill{std::bind(std::uniform_real_distribution<double>(0.1,1.0),
                            std::mt19937(seed_fill))};
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     for(int32_t __i{0}; __i != n_iters; ++__i)
     {
         double r{rbuff_fill()};
         rand_fillers[__i] = r;
         std::size_t rlen{rbuff_len()};
         if(rlen<lo_limit) rlen = lo_limit;
         else if(rlen>hi_limit) rlen = hi_limit;
         //printf("rlen=%llu\n",rlen);
         fptrs[__i] = reinterpret_cast<double*>(gms_mm_malloc((std::size_t)(sizeof(double)*rlen),64ull));
         __builtin_memset(&fptrs[__i][0],0,(unsigned long)rlen);
         const double filler{rand_fillers[__i]};
         sse_memset_unroll16x_pd(fptrs[__i],filler,rlen);
       
         for(std::size_t __j{0}; __j != rlen; ++__j)
         {
            const double curr_val{fptrs[__i][__j]};
            if(curr_val==0.0)
            {
                 printf(ANSI_COLOR_RED "[UNIT-TEST: random length of buffer]: FAILED,value of %.16f at pos: %llu,iter=%d,buf_len=%llu" ANSI_RESET_ALL "\n",curr_val,__j,__i,rlen);
                 fails_mask[__i] = true;
                 fails_count += 1ull;
            }
         }
        
     }
     if(fails_count == 0ull)
     {
         printf(ANSI_COLOR_GREEN "[UNIT-TEST: random length of buffer]: PASSED!!" ANSI_RESET_ALL "\n");
     }
     else
     {
         for(int32_t __i{0}; __i != n_iters; ++__i)
         {
            printf(ANSI_COLOR_RED "[UNIT-TEST: random length of buffer]: fails mask, iter=%d,mask=%d" ANSI_RESET_ALL "\n",__i,fails_mask[__i]);
         }
         printf(ANSI_COLOR_RED "[UNIT-TEST: random length of buffer]: fails count=%llu" ANSI_RESET_ALL "\n",fails_count);
     }
     for(int32_t __i{0}; __i != n_iters; ++__i)
     {
         if(fptrs[__i] != NULL)
         {
            gms_mm_free(fptrs[__i]);
         }
     }
     printf("[UNIT-TEST]: %s() ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}




int main()
{
    unit_test_sse_memset_unroll8x_pd_rand_alloc();
    unit_test_sse_memset_unroll16x_pd_rand_alloc();
    return 0;
}