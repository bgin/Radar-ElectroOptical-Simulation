


#include <cstdio>
#include <cstdlib>
#include "GMS_sse_memset.h"
#include "GMS_avx512_uncached_memcpy.h"
#include "GMS_malloc.h"

/*
   icpc -o unit_test_avx512_uc_memcpy_ps -fp-model fast=2 -fno-exceptions -qopt-zmm-usage=high -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_avx512_uncached_memcpy.h GMS_avx512_uncached_memcpy.cpp unit_test_avx512_uc_memcpy_ps.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 -qopt-zmm-usage=high \ 
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_avx512_uncached_memcpy.h GMS_avx512_uncached_memcpy.cpp unit_test_avx512_uc_memcpy_ps.cpp

*/

void unit_test_avx512_uc_memcpy_unroll8x_ps();

void unit_test_avx512_uc_memcpy_unroll8x_ps()
{
    using namespace gms::common;

     constexpr std::size_t size_1{1ULL};
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     float data4000_src[size_4000];
     float data4000_dst[size_4000] = {};
     float data1024_src[size_1024];
     float data1024_dst[size_1024] = {};
     float data489_src[size_489];
     float data489_dst[size_489]   = {};
     float data145_src[size_145];
     float data145_dst[size_145]   = {};
     float data16_src[size_16];
     float data16_dst[size_16]     = {};
     float data5_src[size_5];
     float data5_dst[size_5]       = {};
     float data1_src[size_1];
     float data1_dst[size_1]       = {};    
     struct alignas(64) data_ptrs
     {
            float * __restrict__ pdata4000_src{nullptr};
            float * __restrict__ pdata4000_dst{nullptr};
            float * __restrict__ pdata1024_src{nullptr};
            float * __restrict__ pdata1024_dst{nullptr};
            float * __restrict__ pdata489_src{ nullptr};
            float * __restrict__ pdata489_dst{ nullptr};
            float * __restrict__ pdata145_src{ nullptr};
            float * __restrict__ pdata145_dst{ nullptr};
            float * __restrict__ pdata16_src{  nullptr};
            float * __restrict__ pdata16_dst{  nullptr};
            float * __restrict__ pdata5_src{   nullptr};
            float * __restrict__ pdata5_dst{   nullptr};
            float * __restrict__ pdata1_src{   nullptr};
            float * __restrict__ pdata1_dst{   nullptr};
     };

     const float fill{3.14159265358979323846264338328F};

     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     data_ptrs dptrs;
     dptrs.pdata4000_src = &data4000_src[0];
     dptrs.pdata4000_dst = &data4000_dst[0];
     dptrs.pdata1024_src = &data1024_src[0];
     dptrs.pdata1024_dst = &data1024_dst[0];
     dptrs.pdata489_src  = &data489_src[0];
     dptrs.pdata489_dst  = &data489_dst[0];
     dptrs.pdata145_src  = &data145_src[0];
     dptrs.pdata145_dst  = &data145_dst[0];
     dptrs.pdata16_src   = &data16_src[0];
     dptrs.pdata16_dst   = &data16_dst[0];
     dptrs.pdata5_src    = &data5_src[0];
     dptrs.pdata5_dst    = &data5_dst[0];
     dptrs.pdata1_src    = &data1_src[0];
     dptrs.pdata1_dst    = &data1_dst[0];
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.7f\n",size_4000,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata4000_src[0],fill,size_4000);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata4000_dst[0],&dptrs.pdata4000_src[0],size_4000);
     for(std::size_t __i{0ULL}; __i != size_4000; ++__i)
     {
          const float s{dptrs.pdata4000_src[__i]};
          const float d{dptrs.pdata4000_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b4000_fail = true;
               break;
          }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.7f\n",size_1024,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata1024_src[0],fill,size_1024);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata1024_dst[0],&dptrs.pdata1024_src[0],size_1024);
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i)
     {
          const float s{dptrs.pdata1024_src[__i]};
          const float d{dptrs.pdata1024_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b1024_fail = true;
               break;
          }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.7f\n",size_489,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata489_src[0],fill,size_489);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata489_dst[0],&dptrs.pdata489_src[0],size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i)
     {
          const float s{dptrs.pdata489_src[__i]};
          const float d{dptrs.pdata489_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b489_fail = true;
               break;
          }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.7f\n",size_145,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata145_src[0],fill,size_145);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata145_dst[0],&dptrs.pdata145_src[0],size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i)
     {
          const float s{dptrs.pdata145_src[__i]};
          const float d{dptrs.pdata145_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b145_fail = true;
               break;
          }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.7f\n",size_16,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata16_src[0],fill,size_16);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata16_dst[0],&dptrs.pdata16_src[0],size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i)
     {
          const float s{dptrs.pdata16_src[__i]};
          const float d{dptrs.pdata16_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b16_fail = true;
               break;
          }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.7f\n",size_5,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata5_src[0],fill,size_5);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata5_dst[0],&dptrs.pdata5_src[0],size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i)
     {
          const float s{dptrs.pdata5_src[__i]};
          const float d{dptrs.pdata5_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b5_fail = true;
               break;
          }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#7]: -- fill buffer of size: %llu with value=%.7f\n",size_1,fill);
     sse_memset_unroll8x_ps(&dptrs.pdata1_src[0],fill,size_1);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll8x_ps(&dptrs.pdata1_dst[0],&dptrs.pdata1_src[0],size_1);
     for(std::size_t __i{0ULL}; __i != size_1; ++__i)
     {
          const float s{dptrs.pdata1_src[__i]};
          const float d{dptrs.pdata1_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b1_fail = true;
               break;
          }
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}

void unit_test_avx512_uc_memcpy_unroll16x_ps();

void unit_test_avx512_uc_memcpy_unroll16x_ps()
{
       using namespace gms::common;

     constexpr std::size_t size_1{1ULL};
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     float data4000_src[size_4000];
     float data4000_dst[size_4000] = {};
     float data1024_src[size_1024];
     float data1024_dst[size_1024] = {};
     float data489_src[size_489];
     float data489_dst[size_489]   = {};
     float data145_src[size_145];
     float data145_dst[size_145]   = {};
     float data16_src[size_16];
     float data16_dst[size_16]     = {};
     float data5_src[size_5];
     float data5_dst[size_5]       = {};
     float data1_src[size_1];
     float data1_dst[size_1]       = {};    
     struct alignas(64) data_ptrs
     {
            float * __restrict__ pdata4000_src{nullptr};
            float * __restrict__ pdata4000_dst{nullptr};
            float * __restrict__ pdata1024_src{nullptr};
            float * __restrict__ pdata1024_dst{nullptr};
            float * __restrict__ pdata489_src{ nullptr};
            float * __restrict__ pdata489_dst{ nullptr};
            float * __restrict__ pdata145_src{ nullptr};
            float * __restrict__ pdata145_dst{ nullptr};
            float * __restrict__ pdata16_src{  nullptr};
            float * __restrict__ pdata16_dst{  nullptr};
            float * __restrict__ pdata5_src{   nullptr};
            float * __restrict__ pdata5_dst{   nullptr};
            float * __restrict__ pdata1_src{   nullptr};
            float * __restrict__ pdata1_dst{   nullptr};
     };

     const float fill{3.14159265358979323846264338328F};

     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     bool  b1_fail{false};
     data_ptrs dptrs;
     dptrs.pdata4000_src = &data4000_src[0];
     dptrs.pdata4000_dst = &data4000_dst[0];
     dptrs.pdata1024_src = &data1024_src[0];
     dptrs.pdata1024_dst = &data1024_dst[0];
     dptrs.pdata489_src  = &data489_src[0];
     dptrs.pdata489_dst  = &data489_dst[0];
     dptrs.pdata145_src  = &data145_src[0];
     dptrs.pdata145_dst  = &data145_dst[0];
     dptrs.pdata16_src   = &data16_src[0];
     dptrs.pdata16_dst   = &data16_dst[0];
     dptrs.pdata5_src    = &data5_src[0];
     dptrs.pdata5_dst    = &data5_dst[0];
     dptrs.pdata1_src    = &data1_src[0];
     dptrs.pdata1_dst    = &data1_dst[0];
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.7f\n",size_4000,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata4000_src[0],fill,size_4000);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata4000_dst[0],&dptrs.pdata4000_src[0],size_4000);
     for(std::size_t __i{0ULL}; __i != size_4000; ++__i)
     {
          const float s{dptrs.pdata4000_src[__i]};
          const float d{dptrs.pdata4000_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b4000_fail = true;
               break;
          }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.7f\n",size_1024,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata1024_src[0],fill,size_1024);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata1024_dst[0],&dptrs.pdata1024_src[0],size_1024);
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i)
     {
          const float s{dptrs.pdata1024_src[__i]};
          const float d{dptrs.pdata1024_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b1024_fail = true;
               break;
          }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.7f\n",size_489,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata489_src[0],fill,size_489);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata489_dst[0],&dptrs.pdata489_src[0],size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i)
     {
          const float s{dptrs.pdata489_src[__i]};
          const float d{dptrs.pdata489_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b489_fail = true;
               break;
          }
     }
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.7f\n",size_145,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata145_src[0],fill,size_145);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata145_dst[0],&dptrs.pdata145_src[0],size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i)
     {
          const float s{dptrs.pdata145_src[__i]};
          const float d{dptrs.pdata145_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b145_fail = true;
               break;
          }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.7f\n",size_16,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata16_src[0],fill,size_16);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata16_dst[0],&dptrs.pdata16_src[0],size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i)
     {
          const float s{dptrs.pdata16_src[__i]};
          const float d{dptrs.pdata16_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b16_fail = true;
               break;
          }
     }
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.7f\n",size_5,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata5_src[0],fill,size_5);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata5_dst[0],&dptrs.pdata5_src[0],size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i)
     {
          const float s{dptrs.pdata5_src[__i]};
          const float d{dptrs.pdata5_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b5_fail = true;
               break;
          }
     }
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#7]: -- fill buffer of size: %llu with value=%.7f\n",size_1,fill);
     sse_memset_unroll16x_ps(&dptrs.pdata1_src[0],fill,size_1);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     avx512_uncached_memcpy_unroll16x_ps(&dptrs.pdata1_dst[0],&dptrs.pdata1_src[0],size_1);
     for(std::size_t __i{0ULL}; __i != size_1; ++__i)
     {
          const float s{dptrs.pdata1_src[__i]};
          const float d{dptrs.pdata1_dst[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#7]: FAILED, found value of %.7f at pos: %llu, but expected: %.7f" ANSI_RESET_ALL "\n",d,__i,s);
               b1_fail = true;
               break;
          }
     }
     if(b1_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#7]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}




int main()
{
    
    unit_test_avx512_uc_memcpy_unroll8x_ps();
    unit_test_avx512_uc_memcpy_unroll16x_ps();
    return 0;
}