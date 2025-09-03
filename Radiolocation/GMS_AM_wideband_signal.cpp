
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
                 const bool                set_order,  /*  true for: (omega,theta), false for: (theta,omega)*/
                 const int32_t             id,
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
m_order{set_order?1:2},
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
    this->m_cenv_corr  = darray_c4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs);
    this->m_ambig      = darray_c4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs);
    this->m_mod_ambig  = darray_r4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs);    
    this->m_carrier_i  = darray_r4_t(m_T);
    this->m_carrier_q  = darray_r4_t(m_T);
    this->m_cenv_i     = darray_r4_t(m_T);
    this->m_cenv_q     = darray_r4_t(m_T);
    this->m_signal_i   = darray_r4_t(m_T);     
    this->m_signal_q   = darray_r4_t(m_T);
    this->m_cenv_corr_i=darray_r4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs);
    this->m_cenv_corr_q=darray_r4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs);
   
}

gms::radiolocation
::AM_wb_signal_t
::AM_wb_signal_t(const std::complex<float> A0,
                 const std::complex<float> Ac,
                 const AM_wb_signal_darray_params_t & dp,
                 const std::size_t         baude_rate,
                 const std::size_t         N,
                 const std::size_t         T,
                 const std::size_t         nfreqs,
                 const std::size_t         nomegs,
                 const std::size_t         nthets,
                 const bool                set_order, /*  true for: (omega,theta), false for: (theta,omega)*/
                 const int32_t             id,
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
m_order{set_order?1:2}, //1 for: (omega,theta), 2 for: (theta,omega)
m_carrier_ph{0.0f},
m_invT{0.0f},
m_fc{fc},
m_fs{0.0f},
m_sig_width{0.0f},
m_sig_energy{0.0f},
m_Ps{-0.0f},
m_sym_dep{sym_dep}
{
    int32_t prot_   = dp.prot;
    int32_t flags_  = dp.flags;
    int32_t fd_     = dp.fd;
    long    offset_ = dp.offset;
    int32_t fsize_  = dp.offset;

    this->m_code_seq   = darray_r4_t(m_baude_rate,prot_,flags_,fd_,offset_,fsize_);
    this->m_carrier    = darray_c4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv       = darray_c4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_signal     = darray_c4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv_spec  = darray_c4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_sig_samp   = darray_c4_t(m_T*m_N,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv_corr  = darray_c4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs,prot_,flags_,fd_,offset_,fsize_);
    this->m_ambig      = darray_c4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs,prot_,flags_,fd_,offset_,fsize_);
    this->m_mod_ambig  = darray_r4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs,prot_,flags_,fd_,offset_,fsize_);    
    this->m_carrier_i  = darray_r4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_carrier_q  = darray_r4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv_i     = darray_r4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv_q     = darray_r4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_signal_i   = darray_r4_t(m_T,prot_,flags_,fd_,offset_,fsize_);     
    this->m_signal_q   = darray_r4_t(m_T,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv_corr_i=darray_r4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs,prot_,flags_,fd_,offset_,fsize_);
    this->m_cenv_corr_q=darray_r4_t(set_order?m_nomegs*m_nthets:m_nthets*m_nomegs,prot_,flags_,fd_,offset_,fsize_);
}

gms::radiolocation
::AM_wb_signal_t
::AM_wb_signal_t(AM_wb_signal_t &&other)
:
m_A0{std::move(other.m_A0)},
m_Ac{std::move(other.m_A0)},
m_SNR{std::move(other.m_SNR)},
m_baude_rate{std::move(other.m_baude_rate)},
m_N{std::move(other.m_N)},
m_T{std::move(other.m_T)},
m_nfreqs{std::move(other.m_nfreqs)},
m_nsamples{std::move(other.m_nsamples)},
m_nomegs{std::move(other.m_nomegs)},
m_nthets{std::move(other.m_nthets)},
m_id{std::move(other.m_id)},
m_order{std::move(other.m_id)}, //1 for: (omega,theta), 2 for: (theta,omega)
m_carrier_ph{std::move(other.m_carrier_ph)},
m_invT{std::move(other.m_invT)},
m_fc{std::move(other.m_fc)},
m_fs{std::move(other.m_fs)},
m_sig_width{std::move(other.m_sig_width)},
m_sig_energy{std::move(other.m_sig_energy)},
m_Ps{std::move(other.m_Ps)},
m_sym_dep{std::move(other.m_sym_dep)},
m_code_seq{std::move(other.m_code_seq)},
m_carrier{std::move(other.m_carrier)},
m_cenv{std::move(other.m_cenv)},
m_signal{std::move(other.m_signal)},
m_cenv_spec{std::move(other.m_cenv_spec)},
m_sig_samp{std::move(other.m_sig_samp)},
m_cenv_corr{std::move(other.m_cenv_corr)},
m_ambig{std::move(other.m_ambig)},
m_mod_ambig{std::move(other.m_mod_ambig)},
m_carrier_i{std::move(other.m_carrier_i)},
m_carrier_q{std::move(other.m_carrier_q)},
m_cenv_i{std::move(other.m_cenv_i)},
m_cenv_q{std::move(other.m_cenv_q)},
m_signal_i{std::move(other.m_signal_i)},
m_signal_q{std::move(other.m_signal_q)},
m_cenv_corr_i{std::move(other.m_cenv_corr_i)},
m_cenv_corr_q{std::move(other.m_cenv_corr_q)}
{ 

}


gms::radiolocation
::AM_wb_signal_t &
gms::radiolocation
::AM_wb_signal_t::operator=(AM_wb_signal_t &&other)
{
    if(this==&other) { return (*this);}
    
    this->m_A0         = other.m_A0;
    this->m_Ac         = other.m_A0;
    this->m_SNR        = other.m_SNR;
    this->m_baude_rate = other.m_baude_rate;
    this->m_N          = other.m_N;
    this->m_T          = other.m_T;
    this->m_nfreqs     = other.m_nfreqs;   
    this->m_nsamples   = other.m_nsamples;
    this->m_nomegs     = other.m_nomegs;
    this->m_nthets     = other.m_nthets;
    this->m_id         = other.m_id;
    this->m_order      = other.m_id; 
    this->m_carrier_ph = other.m_carrier_ph;   
    this->m_invT       = other.m_invT;
    this->m_fc         = other.m_fc;
    this->m_fs         = other.m_fs; 
    this->m_sig_width  = other.m_sig_width;
    this->m_sig_energy = other.m_sig_energy;
    this->m_Ps         = other.m_Ps;
    this->m_sym_dep    = other.m_sym_dep;
    this->m_code_seq.operator=(std::move(other.m_code_seq));
    this->m_carrier.operator=(std::move(other.m_carrier)); 
    this->m_cenv.operator=(std::move(other.m_cenv));
    this->m_signal.operator=(std::move(other.m_signal));
    this->m_cenv_spec.operator=(std::move(other.m_cenv_spec));
    this->m_sig_samp.operator=(std::move(other.m_sig_samp));
    this->m_cenv_corr.operator=(std::move(other.m_cenv_corr));
    this->m_ambig.operator=(std::move(other.m_ambig));
    this->m_mod_ambig.operator=(std::move(other.m_mod_ambig));
    this->m_carrier_i.operator=(std::move(other.m_carrier_i));
    this->m_carrier_q.operator=(std::move(other.m_carrier_q));
    this->m_cenv_i.operator=(std::move(other.m_cenv_i));  
    this->m_cenv_q.operator=(std::move(other.m_cenv_q));
    this->m_signal_i.operator=(std::move(other.m_signal_i));  
    this->m_signal_q.operator=(std::move(other.m_signal_q));
    this->m_cenv_corr_i.operator=(std::move(other.m_cenv_corr_i));
    this->m_cenv_corr_q.operator=(std::move(other.m_cenv_corr_q));
    return (*this);
}

gms::radiolocation
::AM_wb_signal_t
::~AM_wb_signal_t()
{

}





