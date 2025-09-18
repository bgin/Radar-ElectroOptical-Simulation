#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_dyn_array.h"
#include "GMS_trapezoid_waveform.h"

/*
   icpc -o unit_test_trapezw_single_v2 -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_trapezoid_waveform.h GMS_trapezoid_waveform.cpp unit_test_single_trapezoid_wave_v2.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_trapezoid_waveform.h GMS_trapezoid_waveform.cpp unit_test_single_trapezoid_wave_v2.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_single_v2_1024();

void unit_test_create_trapezw_single_v2_1024()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr std::uint32_t n_a_param{32ull};
     constexpr std::uint32_t n_l_param{32ull};
     constexpr std::uint32_t n_c_param{32ull};
     constexpr std::uint32_t n_m_param{32ull};
     constexpr float a{10.0f};
     constexpr float m{5.0f};
     constexpr float l{2.0f};
     constexpr float c{2.0f};
          
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     trapezoid_waveform_t __trapezw_1__ = trapezoid_waveform_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating series of trapezoid waves coded.\n");
     __trapezw_1__.single_trapezoid_wave_v2(a,m,l,c,2);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "single_trapezoid_wave_v2_1024_test_1","Series_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_create_trapezw_single_v2_128();

void unit_test_create_trapezw_single_v2_128()
{
     using namespace gms::radiolocation;
     using namespace gms;
     constexpr std::size_t   n_samples{128ull};
     constexpr std::uint32_t n_waves{16ull};
     constexpr std::uint32_t n_a_param{32ull}; // not used here
     constexpr std::uint32_t n_l_param{32ull}; // same
     constexpr std::uint32_t n_c_param{32ull}; // same 
     constexpr std::uint32_t n_m_param{32ull}; // same 
     constexpr float a{1.0f};
     constexpr float m{0.5f};
     constexpr float l{0.1f};
     constexpr float c{0.5f};
          
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     trapezoid_waveform_t __trapezw_1__ = trapezoid_waveform_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);
     char * ctor_name{gms::common::demangle(typeid(__trapezw_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__trapezw_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }   
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_1__.__trapezw_samples__.info_size_alignment();
     std::printf("[UNIT-TEST:] -- Creating series of trapezoid waves coded.\n");
     __trapezw_1__.single_trapezoid_wave_v2(a,m,l,c,2);
     std::printf("[UNIT-TEST:] -- Creating gnuplot plotting command file.\n");
     trapezoid_waveform_t::create_signal_plot(__trapezw_1__.__n_samples__,__trapezw_1__.__trapezw_samples__.m_data,nullptr,
                                              "single_trapezoid_wave_v2_128_test_1","Single_Trapezoid_Waveform",false);
     printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}


int main()
{
     std::clock_t seed = std::clock();
     int32_t which{-1};
     auto rand_d{std::bind(std::uniform_int_distribution<int32_t>(0,1),std::mt19937(seed))};
     which = rand_d();
     switch (which)
     {
         case 0 : 
            unit_test_create_trapezw_single_v2_128();
         break;
         case 1 : 
            unit_test_create_trapezw_single_v2_1024();
         break;
         default : 
            std::printf("[UNIT-TEST:] -- Invalid switch variable = %d\n",which);
            std::terminate();
     }
     
     
     return 0;
}