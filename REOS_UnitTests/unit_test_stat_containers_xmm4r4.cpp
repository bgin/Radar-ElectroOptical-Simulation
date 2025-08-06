


#include <cstdio>
#include <cstdlib>
#include "GMS_stat_containers_xmm4r4.h"

/*
   icpc -o unit_test_stat_containers_xmm4r4 -std=c++17 -fp-model fast=2 -fno-exceptions -qopt-zmm-usage=high -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp GMS_sse_uncached_memcpy.h GMS_sse_uncached_memcpy.cpp GMS_simd_utils.h GMS_complex_xmm4r4.h GMS_stat_containers_xmm4r4.h unit_test_stat_containers_xmm4r4.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -march=skylake-avx512 -mavx512f -fno-exceptions -falign-functions=32 -qopt-zmm-usage=high \ 
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_sse_memcpy.h GMS_sse_memcpy.cpp GMS_sse_uncached_memcpy.h GMS_sse_uncached_memcpy.cpp GMS_simd_utils.h GMS_complex_xmm4r4.h GMS_stat_containers_xmm4r4.h unit_test_stat_containers_xmm4r4.cpp

*/

void unit_test_stat_containers_SC1D_xmm4c4_t_nt();

void unit_test_stat_containers_SC1D_xmm4c4_t_nt()
{
    using namespace gms;
    using namespace gms::common;

     
     
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     constexpr std::size_t size_65536{65536ull};
     constexpr float fill{3.14159265358979323846264338328f};
     gms::math::xmm4c4_t data65536_src[size_65536] = {_mm_set1_ps(fill),_mm_set1_ps(fill)};
     gms::math::xmm4c4_t data4000_src[size_4000]   = {_mm_set1_ps(fill),_mm_set1_ps(fill)};
     gms::math::xmm4c4_t data1024_src[size_1024]   = {_mm_set1_ps(fill),_mm_set1_ps(fill)};
     
     SC1D_xmm4c4_t<size_65536,1> IQ_size65536;
     SC1D_xmm4c4_t<size_4000,1> IQ_size4000;
     SC1D_xmm4c4_t<size_1024,1> IQ_size1024;
     
             
     
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.7f\n",size_65536,fill);
     IQ_size65536.copy_fill(&data65536_src[0],size_65536);

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.7f\n",size_4000,fill);
     IQ_size4000.copy_fill(&data4000_src[0],size_4000);

     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.7f\n",size_1024,fill);
     IQ_size1024.copy_fill(&data1024_src[0],size_1024);
       
     
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}


void unit_test_stat_containers_SC1D_xmm4c4_t();

void unit_test_stat_containers_SC1D_xmm4c4_t()
{
    using namespace gms;
    using namespace gms::common;

     
     
     constexpr std::size_t size_1024{1024ULL};
     constexpr std::size_t size_4000{4000ULL};
     constexpr std::size_t size_65536{65536ull};
     constexpr float fill{3.14159265358979323846264338328f};
     gms::math::xmm4c4_t data65536_src[size_65536] = {_mm_set1_ps(fill), _mm_set1_ps(fill)};
     gms::math::xmm4c4_t data4000_src[size_4000]   = {_mm_set1_ps(fill), _mm_set1_ps(fill)};
     gms::math::xmm4c4_t data1024_src[size_1024]   = {_mm_set1_ps(fill), _mm_set1_ps(fill)};
     
     SC1D_xmm4c4_t<size_65536,0> IQ_size65536;
     SC1D_xmm4c4_t<size_4000,0> IQ_size4000;
     SC1D_xmm4c4_t<size_1024,0> IQ_size1024;
     
       
     
     printf("[UNIT-TEST]: %s ---> STARTED\n", __PRETTY_FUNCTION__);
     printf("[UNIT-TEST--#1]: -- fill buffer of size: %llu with value=%.7f\n",size_65536,fill);
     IQ_size65536.copy_fill(&data65536_src[0],size_65536);

     printf("[UNIT-TEST--#2]: -- fill buffer of size: %llu with value=%.7f\n",size_4000,fill);
    
     IQ_size4000.copy_fill(&data4000_src[0],size_4000);

     printf("[UNIT-TEST--#3]: -- fill buffer of size: %llu with value=%.7f\n",size_1024,fill);
     IQ_size1024.copy_fill(&data1024_src[0],size_1024);
     
       
     printf("[UNIT-TEST]: %s ---> ENDED CORRECTLY\n", __PRETTY_FUNCTION__);
}




int main()
{
    
    unit_test_stat_containers_SC1D_xmm4c4_t_nt();
    unit_test_stat_containers_SC1D_xmm4c4_t();
   
    return 0;
}