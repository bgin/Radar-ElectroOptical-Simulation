


#include <cstdio>
#include <cstdlib>
#include <vector>
#include <valarray>
#include "GMS_stat_containers_xmm2r8.h"



/*
   icpc -o unit_test_stat_containers_xmm2r8_2 -std=c++17 -fp-model fast=2 -fno-exceptions -qopt-zmm-usage=high -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp GMS_sse_uncached_memcpy.h GMS_sse_uncached_memcpy.cpp GMS_simd_utils.h GMS_complex_xmm2r8.h GMS_stat_containers_xmm2r8.h unit_test_stat_containers_xmm2r8_2.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 -qopt-zmm-usage=high \ 
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp GMS_sse_uncached_memcpy.h GMS_sse_uncached_memcpy.cpp GMS_simd_utils.h GMS_complex_xmm2r8.h GMS_stat_containers_xmm2r8.h unit_test_stat_containers_xmm2r8_2.cpp

*/

void unit_test_stat_containers_SC1D_xmm2c8_t_nt_2();

void unit_test_stat_containers_SC1D_xmm2c8_t_nt_2()
{
    using namespace gms;
    using namespace gms::common;

     
     
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     constexpr std::size_t size_65536{65536ull};
     constexpr double fill{3.14159265358979323846264338328};
     std::vector<gms::math::xmm2c8_t> data65536_src;
     std::vector<gms::math::xmm2c8_t> data4000_src;  
     std::vector<gms::math::xmm2c8_t> data1024_src;
     
     SC1D_xmm2c8_t<size_65536,1> IQ_size65536;
     SC1D_xmm2c8_t<size_4000,1> IQ_size4000;
     SC1D_xmm2c8_t<size_1024,1> IQ_size1024;
     
             
     
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.17f\n",size_65536,fill);
     data65536_src.reserve(size_65536);
     for(std::size_t __i{0ull}; __i != data65536_src.size(); ++__i) 
     {
        data65536_src[__i].re = _mm_set1_pd(fill);
        data65536_src[__i].im = _mm_set1_pd(fill);
     }
     IQ_size65536.copy_fill(&data65536_src[0],size_65536);

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.17f\n",size_4000,fill);
     data4000_src.reserve(size_4000);
     for(std::size_t __i{0ull}; __i != data4000_src.size(); ++__i) 
     {
        data4000_src[__i].re = _mm_set1_pd(fill);
        data4000_src[__i].im = _mm_set1_pd(fill);
     }
     IQ_size4000.copy_fill(&data4000_src[0],size_4000);

     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.17f\n",size_1024,fill);
     data1024_src.reserve(size_1024);
     for(std::size_t __i{0ull}; __i != data1024_src.size(); ++__i) 
     {
        data1024_src[__i].re = _mm_set1_pd(fill);
        data1024_src[__i].im = _mm_set1_pd(fill);
     }
     IQ_size1024.copy_fill(&data1024_src[0],size_1024);
       
     
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}


void unit_test_stat_containers_SC1D_xmm2c8_t_2();

void unit_test_stat_containers_SC1D_xmm2c8_t_2()
{
    using namespace gms;
    using namespace gms::common;

     
     
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     constexpr std::size_t size_65536{65536ull};
     constexpr double fill{3.14159265358979323846264338328};
     std::valarray<gms::math::xmm2c8_t> data65536_src(size_65536);
     std::valarray<gms::math::xmm2c8_t> data4000_src(size_4000);
     std::valarray<gms::math::xmm2c8_t> data1024_src(size_1024);
     
     SC1D_xmm2c8_t<size_65536,0> IQ_size65536;
     SC1D_xmm2c8_t<size_4000,0> IQ_size4000;
     SC1D_xmm2c8_t<size_1024,0> IQ_size1024;
     
       
     
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.17f\n",size_65536,fill);
     for(std::size_t __i{0ull}; __i != data65536_src.size(); ++__i) 
     {
        data65536_src[__i].re = _mm_set1_pd(fill);
        data65536_src[__i].im = _mm_set1_pd(fill);
     }
     IQ_size65536.copy_fill(&data65536_src[0],size_65536);

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.17f\n",size_4000,fill);
     for(std::size_t __i{0ull}; __i != data4000_src.size(); ++__i) 
     {
        data4000_src[__i].re = _mm_set1_pd(fill);
        data4000_src[__i].im = _mm_set1_pd(fill);
     }
     IQ_size4000.copy_fill(&data4000_src[0],size_4000);

     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.17f\n",size_1024,fill);
     for(std::size_t __i{0ull}; __i != data1024_src.size(); ++__i) 
     {
        data1024_src[__i].re = _mm_set1_pd(fill);
        data1024_src[__i].im = _mm_set1_pd(fill);
     }
     IQ_size1024.copy_fill(&data1024_src[0],size_1024);
     
       
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}




int main()
{
    
    unit_test_stat_containers_SC1D_xmm2c8_t_nt_2();
    unit_test_stat_containers_SC1D_xmm2c8_t_2();
   
    return 0;
}