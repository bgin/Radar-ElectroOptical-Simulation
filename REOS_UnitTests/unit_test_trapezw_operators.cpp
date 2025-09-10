#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_trapezoid_waveform.h"

/*
   icpc -o unit_test_trapezw_operators -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_trapezoid_waveform.h GMS_trapezoid_waveform.cpp unit_test_trapezw_operators.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_trapezoid_waveform.h GMS_trapezoid_waveform.cpp unit_test_trapezw_operators.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_trapezw_move_op();

void unit_test_trapezw_move_op()
{
     using namespace gms::radiolocation;
     constexpr std::size_t   n_samples{1024ull};
     constexpr std::uint32_t n_waves{32ull};
     constexpr std::uint32_t n_a_param{32ull};
     constexpr std::uint32_t n_l_param{32ull};
     constexpr std::uint32_t n_c_param{32ull};
     constexpr std::uint32_t n_m_param{32ull};

     constexpr std::size_t   n_samples_2{512ull};
     constexpr std::uint32_t n_waves_2{16ull};
     constexpr std::uint32_t n_a_param_2{16ull};
     constexpr std::uint32_t n_l_param_2{16ull};
     constexpr std::uint32_t n_c_param_2{16ull};
     constexpr std::uint32_t n_m_param_2{16ull};
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     trapezoid_waveform_t __trapezw_1__ = trapezoid_waveform_t(n_samples,n_waves,
                                                               n_a_param,n_l_param,
                                                               n_c_param,n_m_param);

     trapezoid_waveform_t __trapezw_2__ = trapezoid_waveform_t(n_samples_2,n_waves_2,
                                                               n_a_param_2,n_l_param_2,
                                                               n_c_param_2,n_m_param_2);

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
     __trapezw_2__.operator=(std::move(__trapezw_1__));
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __trapezw_2__.__trapezw_samples__.info_size_alignment();
     
    std::printf("[UNIT-TEST]: Finished a dump of array members.\n");
    printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);

}

int main()
{
    unit_test_trapezw_move_op();
    return 0;
}
