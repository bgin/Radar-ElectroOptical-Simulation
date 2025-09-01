
#include "GMS_AM_wideband_signal.h"


std::string gms::radiolocation::AM_wb_signal_t::m_signal_name = 
                    std::string("Amplitude_Modulated_wideband_signal");


gms::radiolocation
::AM_wb_signal_t
::AM_wb_signal_t(const std::complex<float> A0,
                 const std::complex<float> Ac,
                 const std::size_t         baude_rate,
                 const std::size_t         N,
                 const std::size_t         T,
                 const std::size_t         nfreqs,
                 const std::size_t         nomegs,
                 const std::size_t         nthets,
                 const int32_t             id,
                 const int32_t             order,
                 const float               fc,
                 const bool                sym_dep)    
:
m_A0{A0},
m_Ac{Ac},
m_SNR{std::complex<float>(0.0f,0.0f)},
m_baude_rate{baude_rate},
m_N{N},
m_T{T},
m_nfreqs{nfreqs},
m_nsamples{m_T*m_N},
m_nomegs{nomegs},
m_nthets{nthets},
m_id{id},
m_order{order},
m_carrier_ph{0.0f},
m_invT{0.0f},
m_fc{fc},
m_fs{0.0f},
m_sig_width{0.0f},
m_sig_energy{0.0f},
m_Ps{-0.0f},
m_sym_dep{sym_dep}
{
    
    this->m_code_seq   = darray_r4_t(m_baude_rate);
    this->m_carrier    = darray_c4_t(m_T);
    this->m_cenv       = darray_c4_t(m_T);
    this->m_signal     = darray_c4_t(m_T);
    this->m_cenv_spec  = darray_c4_t(m_T);
    this->m_sig_samp   = darray_c4_t(m_T*m_N);
    this->m_cenv_corr  = darray_c4_t(m_nomegs*m_nthets);
    this->m_ambig      = darray_c4_t(m_nomegs*m_nthets);
    this->m_mod_ambig  = darray_r4_t(m_nomegs*m_nthets);    
    this->m_carrier_i  = darray_r4_t(m_T);
    this->m_carrier_q  = darray_r4_t(m_T);
    this->m_cenv_i     = darray_r4_t(m_T);
    this->m_cenv_q     = darray_r4_t(m_T);
    this->m_signal_i   = darray_r4_t(m_T);     
    this->m_signal_q   = darray_r4_t(m_T);
    this->m_cenv_corr_i=darray_r4_t(m_nomegs*m_nthets);
    this->m_cenv_corr_q=darray_r4_t(m_nomegs*m_nthets);
   
}

gms::radiolocation
::AM_wb_signal_t
::~AM_wb_signal_t()
{}





