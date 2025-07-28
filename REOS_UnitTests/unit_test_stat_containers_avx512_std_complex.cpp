


#include <cstdio>
#include <cstdlib>
#include "GMS_stat_containers.h"
#include "GMS_sse_memset.h"
#include "GMS_avx512_uncached_memcpy.h"
#include "GMS_avx512_memcpy.h"
#include "GMS_malloc.h"

/*
   icpc -o unit_test_stat_containers_avx512_std_complex -std=c++17 -fp-model fast=2 -fno-exceptions -qopt-zmm-usage=high -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_avx512_memcpy.h GMS_avx512_memcpy.cpp GMS_avx512_uncached_memcpy.h GMS_avx512_uncached_memcpy.cpp GMS_stat_containers.hpp unit_test_stat_containers_avx512_std_complex
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 -qopt-zmm-usage=high \ 
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_avx512_memcpy.h GMS_avx512_memcpy.cpp GMS_avx512_uncached_memcpy.h GMS_avx512_uncached_memcpy.cpp GMS_stat_containers.hpp unit_test_stat_containers_avx512_std_complex

*/

void unit_test_stat_containers_StatC1D_c4_t();

void unit_test_stat_containers_StatC1D_c4_t()
{
    using namespace gms;
    using namespace gms::common;

     
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     
     std::complex<float> data4000_src[size_4000];
     std::complex<float> data1024_src[size_1024];
     std::complex<float> data489_src[size_489];
     std::complex<float> data145_src[size_145];
     std::complex<float> data16_src[size_16];
     std::complex<float> data5_src[size_5];
     
     StatC1D_c4_t<size_4000,0> IQ_size4000;
     StatC1D_c4_t<size_1024,0> IQ_size1024;
     StatC1D_c4_t<size_489,0>  IQ_size489;
     StatC1D_c4_t<size_145,0>  IQ_size145;
     StatC1D_c4_t<size_16,0>   IQ_size16;
     StatC1D_c4_t<size_5,0>    IQ_size5;
     struct alignas(64) data_ptrs
     {
            std::complex<float> * __restrict__ pdata4000_src{nullptr};
            std::complex<float> * __restrict__ pdata1024_src{nullptr};
            std::complex<float> * __restrict__ pdata489_src{ nullptr};
            std::complex<float> * __restrict__ pdata145_src{ nullptr};
            std::complex<float> * __restrict__ pdata16_src{  nullptr};
            std::complex<float> * __restrict__ pdata5_src{   nullptr};
                       
     };

     constexpr float fill{3.14159265358979323846264338328};

     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     data_ptrs dptrs;
     dptrs.pdata4000_src = &data4000_src[0];
     dptrs.pdata1024_src = &data1024_src[0];
     dptrs.pdata489_src  = &data489_src[0];
     dptrs.pdata145_src  = &data145_src[0];
     dptrs.pdata16_src   = &data16_src[0];
     dptrs.pdata5_src    = &data5_src[0];
     
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.17f\n",size_4000,fill);
     sse_memset_unroll8x_ps(reinterpret_cast<float * __restrict__>(&dptrs.pdata4000_src[0]),fill,2ull*size_4000);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size4000.copy_fill(&dptrs.pdata4000_src[0],size_4000);

     for(std::size_t __i{0ULL}; __i != size_4000; ++__i)
     {
          std::complex<float> s{dptrs.pdata4000_src[__i]};
          std::complex<float> d{IQ_size4000.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b4000_fail = true;
               break;
          }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.17f\n",size_1024,fill);
     sse_memset_unroll8x_ps(reinterpret_cast<float * __restrict>(&dptrs.pdata1024_src[0]),fill,2ull*size_1024);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size1024.copy_fill(&dptrs.pdata1024_src[0],size_1024);
     
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i)
     {
          std::complex<float> s{dptrs.pdata1024_src[__i]};
          std::complex<float> d{IQ_size1024.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b1024_fail = true;
               break;
          }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.17f\n",size_489,fill);
     sse_memset_unroll8x_ps(reinterpret_cast<float * __restrict>(&dptrs.pdata489_src[0]),fill,2ull*size_489);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size489.copy_fill(&dptrs.pdata489_src[0],size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i)
     {
          std::complex<float> s{dptrs.pdata489_src[__i]};
          std::complex<float> d{IQ_size489.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b489_fail = true;
               break;
          }
     }
    
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.17f\n",size_145,fill);
     sse_memset_unroll8x_ps(reinterpret_cast<float * __restrict>(&dptrs.pdata145_src[0]),fill,2ull*size_145);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size145.copy_fill(&dptrs.pdata145_src[0],size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i)
     {
          std::complex<float> s{dptrs.pdata145_src[__i]};
          std::complex<float> d{IQ_size145.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b145_fail = true;
               break;
          }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.17f\n",size_16,fill);
     sse_memset_unroll8x_ps(reinterpret_cast<float * __restrict>(&dptrs.pdata16_src[0]),fill,2ull*size_16);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size16.copy_fill(&dptrs.pdata16_src[0],size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i)
     {
          std::complex<float> s{dptrs.pdata16_src[__i]};
          std::complex<float> d{IQ_size16.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b16_fail = true;
               break;
          }
     }
     
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.17f\n",size_5,fill);
     sse_memset_unroll8x_ps(reinterpret_cast<float * __restrict>(&dptrs.pdata5_src[0]),fill,2ull*size_5);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size5.copy_fill(&dptrs.pdata5_src[0],size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i)
     {
          std::complex<float> s{dptrs.pdata5_src[__i]};
          std::complex<float> d{IQ_size5.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b5_fail = true;
               break;
          }
     }
     
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}

void unit_test_stat_containers_StatC1D_c8_t();

void unit_test_stat_containers_StatC1D_c8_t()
{
    using namespace gms;
    using namespace gms::common;

     
     constexpr std::size_t size_5{5ULL};
     constexpr std::size_t size_16{16ULL};
     constexpr std::size_t size_145{145ULL};
     constexpr std::size_t size_489{489ULL};
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     
     std::complex<double> data4000_src[size_4000];
     std::complex<double> data1024_src[size_1024];
     std::complex<double> data489_src[size_489];
     std::complex<double> data145_src[size_145];
     std::complex<double> data16_src[size_16];
     std::complex<double> data5_src[size_5];
     
     StatC1D_c8_t<size_4000,0> IQ_size4000;
     StatC1D_c8_t<size_1024,0> IQ_size1024;
     StatC1D_c8_t<size_489,0>  IQ_size489;
     StatC1D_c8_t<size_145,0>  IQ_size145;
     StatC1D_c8_t<size_16,0>   IQ_size16;
     StatC1D_c8_t<size_5,0>    IQ_size5;
     struct alignas(64) data_ptrs
     {
            std::complex<double> * __restrict__ pdata4000_src{nullptr};
            std::complex<double> * __restrict__ pdata1024_src{nullptr};
            std::complex<double> * __restrict__ pdata489_src{ nullptr};
            std::complex<double> * __restrict__ pdata145_src{ nullptr};
            std::complex<double> * __restrict__ pdata16_src{  nullptr};
            std::complex<double> * __restrict__ pdata5_src{   nullptr};
                       
     };

     constexpr float fill{3.14159265358979323846264338328};

     bool  b4000_fail{false};
     bool  b1024_fail{false};
     bool  b489_fail{false};
     bool  b145_fail{false};
     bool  b16_fail{false};
     bool  b5_fail{false};
     data_ptrs dptrs;
     dptrs.pdata4000_src = &data4000_src[0];
     dptrs.pdata1024_src = &data1024_src[0];
     dptrs.pdata489_src  = &data489_src[0];
     dptrs.pdata145_src  = &data145_src[0];
     dptrs.pdata16_src   = &data16_src[0];
     dptrs.pdata5_src    = &data5_src[0];
     
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.17f\n",size_4000,fill);
     sse_memset_unroll8x_pd(reinterpret_cast<double * __restrict__>(&dptrs.pdata4000_src[0]),fill,2ull*size_4000);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size4000.copy_fill(&dptrs.pdata4000_src[0],size_4000);

     for(std::size_t __i{0ULL}; __i != size_4000; ++__i)
     {
          std::complex<double> s{dptrs.pdata4000_src[__i]};
          std::complex<double> d{IQ_size4000.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#1]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b4000_fail = true;
               break;
          }
     }
     if(b4000_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.17f\n",size_1024,fill);
     sse_memset_unroll8x_pd(reinterpret_cast<double * __restrict>(&dptrs.pdata1024_src[0]),fill,2ull*size_1024);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size1024.copy_fill(&dptrs.pdata1024_src[0],size_1024);
     
     for(std::size_t __i{0ULL}; __i != size_1024; ++__i)
     {
          std::complex<double> s{dptrs.pdata1024_src[__i]};
          std::complex<double> d{IQ_size1024.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#2]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b1024_fail = true;
               break;
          }
     }
     if(b1024_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#2]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.17f\n",size_489,fill);
     sse_memset_unroll8x_pd(reinterpret_cast<double * __restrict>(&dptrs.pdata489_src[0]),fill,2ull*size_489);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size489.copy_fill(&dptrs.pdata489_src[0],size_489);
     for(std::size_t __i{0ULL}; __i != size_489; ++__i)
     {
          std::complex<double> s{dptrs.pdata489_src[__i]};
          std::complex<double> d{IQ_size489.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#3]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b489_fail = true;
               break;
          }
     }
    
     if(b489_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#1]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#4]: -- fill buffer of size: %llu with value=%.17f\n",size_145,fill);
     sse_memset_unroll8x_pd(reinterpret_cast<double * __restrict>(&dptrs.pdata145_src[0]),fill,2ull*size_145);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size145.copy_fill(&dptrs.pdata145_src[0],size_145);
     for(std::size_t __i{0ULL}; __i != size_145; ++__i)
     {
          std::complex<double> s{dptrs.pdata145_src[__i]};
          std::complex<double> d{IQ_size145.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#4]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b145_fail = true;
               break;
          }
     }
     if(b145_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#4]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#5]: -- fill buffer of size: %llu with value=%.17f\n",size_16,fill);
     sse_memset_unroll8x_pd(reinterpret_cast<double * __restrict>(&dptrs.pdata16_src[0]),fill,2ull*size_16);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size16.copy_fill(&dptrs.pdata16_src[0],size_16);
     for(std::size_t __i{0ULL}; __i != size_16; ++__i)
     {
          std::complex<double> s{dptrs.pdata16_src[__i]};
          std::complex<double> d{IQ_size16.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#5]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b16_fail = true;
               break;
          }
     }
     
     if(b16_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#5]: PASSED!!" ANSI_RESET_ALL "\n");}
     printf("[UNIT-TEST--#6]: -- fill buffer of size: %llu with value=%.17f\n",size_5,fill);
     sse_memset_unroll8x_pd(reinterpret_cast<double * __restrict>(&dptrs.pdata5_src[0]),fill,2ull*size_5);
     __asm__ __volatile__ ( "vzeroupper" : : : );
     IQ_size5.copy_fill(&dptrs.pdata5_src[0],size_5);
     for(std::size_t __i{0ULL}; __i != size_5; ++__i)
     {
          std::complex<double> s{dptrs.pdata5_src[__i]};
          std::complex<double> d{IQ_size5.mx[__i]};
          if(d != s)
          {
               printf(ANSI_COLOR_RED "[UNIT-TEST--#6]: FAILED, found value of (%.17f,%.17f) at pos: %llu, but expected: (%.17f,%.17f)" 
                      ANSI_RESET_ALL "\n",d.real(),d.imag(),__i,s.real(),s.imag());
               b5_fail = true;
               break;
          }
     }
     
     if(b5_fail==false) {printf(ANSI_COLOR_GREEN "[UNIT-TEST--#6]: PASSED!!" ANSI_RESET_ALL "\n");}
     
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}


int main()
{
    
    unit_test_stat_containers_StatC1D_c4_t();
    unit_test_stat_containers_StatC1D_c8_t();
    return 0;
}