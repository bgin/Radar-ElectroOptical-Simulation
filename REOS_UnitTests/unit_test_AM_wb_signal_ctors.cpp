
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "GMS_AM_wideband_signal.h"

/*
   icpc -o unit_test_AM_wb_signal_ctors -fp-model -std=c++17 fast=2 -ftz -ggdb -ipo -march=skylake-avx512 -mavx512f -falign-functions=32 -w1 -qopt-report=5  \
   GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_AM_wideband_signal.h GMS_AM_wideband_signal.cpp unit_test_AM_wb_signal_ctors.cpp
   ASM: 
   icpc -S -fverbose-asm -masm=intel  -std=c++17 -march=skylake-avx512 -mavx512f -falign-functions=32 GMS_config.h GMS_malloc.h GMS_fast_pmc_access.h  GMS_dyn_array.h GMS_AM_wideband_signal.h GMS_AM_wideband_signal.cpp unit_test_AM_wb_signal_ctors.cpp

*/

__attribute__((hot))
__attribute__((noinline))
void unit_test_AM_wideband_signal_ctor_2();

void unit_test_AM_wideband_signal_ctor_2()
{
     using namespace gms::radiolocation;
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
    printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);

}


__attribute__((hot))
__attribute__((noinline))
void unit_test_AM_wideband_signal_ctor_3();

void unit_test_AM_wideband_signal_ctor_3()
{
     using namespace gms::radiolocation;
     constexpr std::size_t T{128ull};
     constexpr std::size_t N{128ull};
     constexpr std::size_t baude_rate{N};
     constexpr std::size_t nfreqs{512ull};
     constexpr std::size_t nomegs{128ull};
     constexpr std::size_t nthets{128ull};
     AM_wb_signal_darray_params_t dp;
     std::complex<float> Ac{1.0f,0.0f};
     std::complex<float> A0{1.25f,1.25f};
     const int32_t id{1};
     float fc{3.1256e+9f};
     const bool sym_dep{false};
     const bool layout{true};
     const bool set_order{true};
     dp.prot   = PROT_READ | PROT_WRITE;
     dp.flags  = MAP_ANONYMOUS | MAP_PRIVATE;
     dp.fd     = -1;
     dp.offset = 0;
     dp.fsize  = 0;
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     AM_wb_signal_t __symbol_id_1__ = AM_wb_signal_t(Ac,A0,dp,baude_rate,N,T,
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
    printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}

__attribute__((hot))
__attribute__((noinline))
void unit_test_AM_wideband_signal_move_ctor();

void unit_test_AM_wideband_signal_move_ctor()
{
     using namespace gms::radiolocation;
     constexpr std::size_t T{128ull};
     constexpr std::size_t N{128ull};
     constexpr std::size_t baude_rate{N};
     constexpr std::size_t nfreqs{512ull};
     constexpr std::size_t nomegs{128ull};
     constexpr std::size_t nthets{128ull};
     AM_wb_signal_darray_params_t dp;
     std::complex<float> Ac{1.0f,0.0f};
     std::complex<float> A0{1.25f,1.25f};
     const int32_t id{1};
     float fc{3.1256e+9f};
     const bool sym_dep{false};
     const bool layout{true};
     const bool set_order{true};
     dp.prot   = PROT_READ | PROT_WRITE;
     dp.flags  = MAP_ANONYMOUS | MAP_PRIVATE;
     dp.fd     = -1;
     dp.offset = 0;
     dp.fsize  = 0;
     int32_t status{};
     //__asm__ ("int3");
     printf("[UNIT_TEST]: function=%s -- **START**\n", __PRETTY_FUNCTION__);
     AM_wb_signal_t __symbol_id_1__ = AM_wb_signal_t(Ac,A0,dp,baude_rate,N,T,
                                                     nfreqs,nomegs,nthets,
                                                     set_order,id,fc,sym_dep);

     AM_wb_signal_t __symbol_id_1_cpy__{std::move(__symbol_id_1__)};

     char * ctor_name{gms::common::demangle(typeid(__symbol_id_1_cpy__).name(),status)};
     if(status==0 && ctor_name != NULL)
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", ctor_name);
          gms::common::gms_mm_free(ctor_name);
     }
     else
     {
          printf("[UNIT-TEST]: Instantiation of object Constructor of type: %s\n\n", typeid(__symbol_id_1_cpy__).name());
          if(ctor_name != NULL) gms::common::gms_mm_free(ctor_name);
     }
     std::printf("[UNIT-TEST:] -- Dumping an array members info and alignment.\n");
     __symbol_id_1_cpy__.m_ambig.info_size_alignment();
     __symbol_id_1_cpy__.m_carrier.info_size_alignment();
     __symbol_id_1_cpy__.m_carrier_i.info_size_alignment();
     __symbol_id_1_cpy__.m_carrier_q.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv_corr.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv_corr_i.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv_corr_q.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv_i.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv_q.info_size_alignment();
     __symbol_id_1_cpy__.m_cenv_spec.info_size_alignment();
     __symbol_id_1_cpy__.m_code_seq.info_size_alignment();
     __symbol_id_1_cpy__.m_mod_ambig.info_size_alignment();
     __symbol_id_1_cpy__.m_sig_samp.info_size_alignment();
     __symbol_id_1_cpy__.m_signal.info_size_alignment();
     __symbol_id_1_cpy__.m_signal_i.info_size_alignment();
     __symbol_id_1_cpy__.m_signal_q.info_size_alignment();
    std::printf("[UNIT-TEST]: Finished a dump of array members.\n");
    printf("[UNIT_TEST]: function=%s -- **END**\n", __PRETTY_FUNCTION__);
}



int main()
{
    //unit_test_AM_wideband_signal_ctor_2();
    //unit_test_AM_wideband_signal_ctor_3();
    unit_test_AM_wideband_signal_move_ctor();
    return 0;
}