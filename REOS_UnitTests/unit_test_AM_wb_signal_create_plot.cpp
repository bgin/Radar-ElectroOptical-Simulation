
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <limits>
#include "GMS_AM_wideband_signal.h"

/*
   icpc -o unit_test_AM_wb_signal_create_plot -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_AM_wideband_signal.h GMS_AM_wideband_signal.cpp unit_test_AM_wb_signal_create_plot.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 \ 
   GMS_config.h GMS_malloc.h GMS_sse_memset.h GMS_sse_memset.cpp GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_AM_wideband_signal.h GMS_AM_wideband_signal.cpp unit_test_AM_wb_signal_create_plot.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_AM_wideband_signal_create_plot();

void unit_test_AM_wideband_signal_create_plot()
{
     using namespace gms::radiolocation;
     constexpr std::complex<float> cfill{std::numeric_limits<float>::quiet_NaN(),
                                         std::numeric_limits<float>::quiet_NaN()};
     constexpr float rfill{std::numeric_limits<float>::quiet_NaN()};
     constexpr std::size_t T{128ull};
     constexpr std::size_t N{128ull};
     constexpr std::size_t baude_rate{N};
     constexpr std::size_t nfreqs{512ull};
     constexpr std::size_t nomegs{128ull};
     constexpr std::size_t nthets{128ull};
     std::complex<float> Ac{1.0f,0.0f};
     std::complex<float> A0{1.25f,1.25f};
    
     const int32_t id{1};
     float fc{3.1256e+9f};
     const bool sym_dep{false};
     const bool layout{true};
     const bool set_order{true};
         
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     AM_wb_signal_t __symbol_id_1__ = AM_wb_signal_t(Ac,A0,baude_rate,N,T,
                                                     nfreqs,nomegs,nthets,
                                                     set_order,id,fc,sym_dep);

     char * ctor_name{gms::common::demangle(typeid(__symbol_id_1__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__symbol_id_1__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __symbol_id_1__.m_ambig.info_size_alignment();
     __symbol_id_1__.m_carrier.info_size_alignment();
     __symbol_id_1__.m_carrier_i.info_size_alignment();
     __symbol_id_1__.m_carrier_q.info_size_alignment();
     __symbol_id_1__.m_cenv.info_size_alignment();
     __symbol_id_1__.m_cenv_corr.info_size_alignment();
     __symbol_id_1__.m_cenv_corr_i.info_size_alignment();
     __symbol_id_1__.m_cenv_corr_q.info_size_alignment();
     __symbol_id_1__.m_cenv_i.info_size_alignment();
     __symbol_id_1__.m_cenv_q.info_size_alignment();
     __symbol_id_1__.m_cenv_spec.info_size_alignment();
     __symbol_id_1__.m_code_seq.info_size_alignment();
     __symbol_id_1__.m_mod_ambig.info_size_alignment();
     __symbol_id_1__.m_sig_samp.info_size_alignment();
     __symbol_id_1__.m_signal.info_size_alignment();
     __symbol_id_1__.m_signal_i.info_size_alignment();
     __symbol_id_1__.m_signal_q.info_size_alignment();
    std::printf("[UNIT-TEST]: Finished a dump of array members.\n");
    std::printf("[UNIT_TEST]: Calling an init_storage member function.\n");
    __symbol_id_1__.init_storage(cfill,rfill);
    std::printf("[UNIT-TEST]:Calling a plotting function.\n");
    AM_wb_signal_t::ceate_signal_plot(__symbol_id_1__.m_T,__symbol_id_1__.m_signal_i.m_data,
                                      __symbol_id_1__.m_signal_q.m_data,"create_signal_plot_test","Test_NaN");
    printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);

}






int main()
{
    unit_test_AM_wideband_signal_create_plot();
    return 0;
}