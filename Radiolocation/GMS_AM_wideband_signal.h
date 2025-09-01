

/*MIT License
!Copyright (c) 2020 Bernard Gingold
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
*/

#ifndef __GMS_AM_WIDEBAND_SIGNAL_H__
#define __GMS_AM_WIDEBAND_SIGNAL_H__ 270820251227


namespace file_info 
{

     static const unsigned int GMS_AM_WIDEBAND_SIGNAL_MAJOR = 1;
     static const unsigned int GMS_AM_WIDEBAND_SIGNAL_MINOR = 1;
     static const unsigned int GMS_AM_WIDEBAND_SIGNAL_MICRO = 0;
     static const unsigned int GMS_AM_WIDEBAND_SIGNAL_FULLVER =
       1000U*GMS_AM_WIDEBAND_SIGNAL_MAJOR+100U*GMS_AM_WIDEBAND_SIGNAL_MINOR+
       10U*GMS_AM_WIDEBAND_SIGNAL_MICRO;
     static const char GMS_AM_WIDEBAND_SIGNAL_CREATION_DATE[] = "27-08-2025 12:27 +00200 (WED  27 AUG 2025 GMT+2)";
     static const char GMS_AM_WIDEBAND_SIGNAL_BUILD_DATE[]    = __DATE__; 
     static const char GMS_AM_WIDEBAND_SIGNAL_BUILD_TIME[]  = __TIME__;
     static const char GMS_AM_WIDEBAND_SIGNAL_SYNOPSIS[]    = "Amplitude-Modulated wideband signal.";

}





#include <cstdint>
#include <string>
#include <complex>
#include "GMS_config.h"
#include "GMS_dyn_array.h"

// Enable non-temporal stores for this class only( used with free-standing operators)
// defaulted to 0.
// To be added.
#if !defined (USE_GMS_AM_WIDEBAND_SIGNAL_NT_STORES)
#define USE_GMS_AM_WIDEBAND_SIGNAL_NT_STORES 0
#endif

// Enable for the basic PMC tracing (wall-clock) readout (not statistically rigorous)!!
// *** Warning *** -- An access for the PM hardware counters must be enabled for the user-mode space!!
// 
#if !defined (AM_WIDEBAND_SIGNAL_USE_PMC_INSTRUMENTATION)
#define AM_WIDEBAND_SIGNAL_USE_PMC_INSTRUMENTATION 1
#endif 

#if (AM_WIDEBAND_SIGNAL_USE_PMC_INSTRUMENTATION) == 1
#include "GMS_hw_perf_macros.h"

#define PMC_VARS                      \
uint64_t prog_counters_start[4] = {}; \
uint64_t prog_counters_stop[4]  = {}; \
uint64_t tsc_start,tsc_stop;          \
uint64_t act_cyc_start,act_cyc_stop;  \
uint64_t ref_cyc_start,ref_cyc_stop;  \
[[maybe_unused]] uint64_t dummy1;     \
[[maybe_unused]] uint64_t dummy2;     \
[[maybe_unused]] uint64_t dummy3;     \
int32_t core_counter_width;           \
double utilization,nom_ghz,avg_ghz;
#endif 

namespace gms 
{

namespace radiolocation 
{
        

          enum class AM_wb_signal_env_type : int32_t 
          {
                     trapezoidal,
                     sine_squared,
                     sine
          };

          enum class AM_wb_signal_rand_gens : int32_t 
          {
                     rg_minstd_rand0,
                     rg_minstd_rand,
                     rg_mt19937,
                     rg_mt19937_64,
                     rg_ranlux24_base,
                     rg_ranlux48_base,
                     rg_ranlux24,
                     rg_ranlux48,
                     rg_knuth_b
          };

          enum class AM_wb_signal_omega_pdf : int32_t 
          {
                     uniform,
                     normal,
                     cauchy,
                     log_normal,
                     exponential,
                     poisson,
                     weibull,
                     gamma
          };

          enum class AM_wb_signal_theta_pdf : int32_t 
          {
                     uniform,
                     normal,
                     cauchy,
                     log_normal,
                     exponential,
                     poisson,
                     weibull,
                     gamma
          };

          enum class AM_wb_signal_code_seq : int32_t 
          {
                     alternate_zero_one,
                     alternate_pos_one_neg_one
          };

          enum class AM_wb_signal_FT_processing : int32_t 
          {
                     quadpack_integration,
                     mkl_fft
          };

         
          struct alignas(64) AM_wb_signal_rdistro_params_t 
          {
                 float uni_real_a_r;
                 float uni_real_a_i;
                 float uni_real_b_r;
                 float uni_real_b_i;
                 float norm_mu_r;
                 float norm_mu_i;
                 float norm_sigma_r;
                 float norm_sigma_i;
                 float cauchy_a_r;
                 float cauchy_a_i;
                 float cauchy_b_r;
                 float cauchy_b_i;
                 float log_norm_m_r;
                 float log_norm_m_i;
                 float log_norm_s_r;
                 float log_norm_s_i;
                 float expo_gamma_r;
                 float expo_gamma_i;
                 float weibull_a_r;
                 float weibull_a_i;
                 float weibull_b_r;
                 float weibull_b_i;
                 float gamma_alph_r;
                 float gamma_alph_i;
                 float gamma_bet_r;
                 float gamma_bet_i;
                 int32_t poisson_mu;

          };

          // Holder for darray ctor parameters
          struct alignas(16) AM_wb_signal_darray_params_t 
          {
                 int32_t prot;
                 int32_t flags;
                 int32_t fd;
                 long offset;
                 int32_t fsize;
          };

           
          struct alignas(64) AM_wb_signal_t final 
          {
                
                
                 darray_r4_t         m_code_seq;    // [0,1 or -1,1] 
                 darray_c4_t         m_carrier;     // carrier signal part
                 darray_c4_t         m_cenv;        // complex envelope 
                 darray_c4_t         m_signal;      // wideband signal (whole)
                 darray_c4_t         m_cenv_spec;   // complex envelope spectrum
                 darray_c4_t         m_sig_samp;    // signal samples
                 darray_c4_t         m_cenv_corr;   // complex envelope correlation function
                 darray_c4_t         m_ambig;      // ambiguity function for: (omega,theta)
                 darray_r4_t         m_mod_ambig;  // module of ambiuity function 
                 darray_r4_t         m_carrier_i;   // I/Q parts 
                 darray_r4_t         m_carrier_q;
                 darray_r4_t         m_cenv_i;
                 darray_r4_t         m_cenv_q;
                 darray_r4_t         m_signal_i;
                 darray_r4_t         m_signal_q;
                 darray_r4_t         m_cenv_corr_i;
                 darray_r4_t         m_cenv_corr_q;

                 std::complex<float> m_A0;         // complex amplitude value 
                 std::complex<float> m_Ac;         // complex carrier amplitude 
                 std::complex<float> m_SNR;        // signal-to-noise ratio
                 std::size_t         m_baude_rate; //number of baude bit changes i.e. -1,1,or 0,1 usually the same as m_Ne value
                 std::size_t         m_N;         //number of narrowband signals (start value)
                 std::size_t         m_T;         // symbol length (period) (end value)
                 std::size_t         m_nfreqs;     // number of frequencies (Spectral analysis) -- start value
                                                 //                                                    N-1
                 std::size_t         m_nsamples;    // number of accumulated samples per symbol i.e. y(t)=Sum a(t-kT)*dr
                                       //                                                    !k=0
                 std::size_t         m_nomegs;     // number of doppler frequency shifts (start value)
                 std::size_t         m_nthets;     // number of time delays (for return signal) (start value)
                 int32_t             m_id;         // this signal ID 
                 int32_t             m_order;      // 1 for: (omega,theta), 2 for: (theta,omega)
                 float               m_carrier_ph;    // initial (randomly) carrier phase
                 float               m_invT;       // inverse of symbol period 'T'
                 float               m_fc;         // carrier frequency
                 float               m_fs;         // sampling frequency
                 float               m_sig_width;  // signal width (integral)
                 float               m_sig_energy; // signal energy (integral)
                 float               m_Ps;         // signal(symbol) error Probability function 
                 bool                m_sym_dep;    // previous-next symbol dependency (page: 50, formula: 2.7)
                 
                 static std::string  m_signal_name;    // signal/symbol  plain name

                 AM_wb_signal_t() = delete;

                 AM_wb_signal_t(const std::complex<float>,
                                const std::complex<float>,
                                const std::size_t,
                                const std::size_t,
                                const std::size_t,
                                const std::size_t,
                                const std::size_t,
                                const std::size_t,
                                const int32_t,
                                const int32_t,
                                const float,
                                const bool) noexcept(false);
                             

               ~AM_wb_signal_t();

          }; 


} //radiolocation

}// gms





















#endif /*__GMS_AM_WIDEBAND_SIGNAL_H__*/