
#include <iostream>
#include <fstream>
#include "GMS_AM_wideband_signal.h"
#include "GMS_sse_memset.h"

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

void 
gms::radiolocation
::AM_wb_signal_t
::init_storage(const std::complex<float> cval,
               const float               rval)
{
    using namespace gms::common;
#if (INIT_BY_STD_FILL) == 0
    sse_memset_unroll8x_ps(&this->m_code_seq.m_data[0],rval,this->m_baude_rate);
    sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_carrier.m_data[0]),rval,2ull*this->m_T);
    sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_cenv.m_data[0]),rval,2ull*this->m_T);
    sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_signal.m_data[0]),rval,2ull*this->m_T);
    sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_cenv_spec.m_data[0]),rval,2ull*this->m_T);
    sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_sig_samp.m_data[0]),rval,2ull*(this->m_T*this->m_N));
    if(this->m_order == 1)
    {
       sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_cenv_corr.m_data[0]),rval,2ull*(this->m_nomegs*this->m_nthets));
       sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_ambig.m_data[0]),rval,2ull*(this->m_nomegs*this->m_nthets));
       sse_memset_unroll8x_ps(&this->m_mod_ambig.m_data[0],rval,this->m_nomegs*this->m_nthets);
       sse_memset_unroll8x_ps(&this->m_cenv_corr_i.m_data[0],rval,this->m_nomegs*this->m_nthets);
       sse_memset_unroll8x_ps(&this->m_cenv_corr_q.m_data[0],rval,this->m_nomegs*this->m_nthets);
    }
    else if(this->m_order == 2)
    {
       sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_cenv_corr.m_data[0]),rval,2ull*(this->m_nthets*this->m_nomegs));
       sse_memset_unroll8x_ps(reinterpret_cast<float *>(&this->m_ambig.m_data[0]),rval,2ull*(this->m_nthets*this->m_nomegs));
       sse_memset_unroll8x_ps(&this->m_mod_ambig.m_data[0],rval,this->m_nthets*this->m_nomegs);
       sse_memset_unroll8x_ps(&this->m_cenv_corr_i.m_data[0],rval,this->m_nthets*this->m_nomegs);
       sse_memset_unroll8x_ps(&this->m_cenv_corr_q.m_data[0],rval,this->m_nthets*this->m_nomegs);
    }
    sse_memset_unroll8x_ps(&this->m_carrier_i.m_data[0],rval,this->m_T);
    sse_memset_unroll8x_ps(&this->m_carrier_q.m_data[0],rval,this->m_T);
    sse_memset_unroll8x_ps(&this->m_cenv_i.m_data[0],rval,this->m_T);
    sse_memset_unroll8x_ps(&this->m_cenv_q.m_data[0],rval,this->m_T);
    sse_memset_unroll8x_ps(&this->m_signal_i.m_data[0],rval,this->m_T);
    sse_memset_unroll8x_ps(&this->m_signal_q.m_data[0],rval,this->m_T);
#else 
    std::fill(this->m_code_seq.m_data,this->m_code_seq.m_data+this->m_T,rval);
    std::fill(this->m_carrier.m_data,this->m_carrier.m_data+this->m_T,cval);
    std::fill(this->m_cenv.m_data,this->m_cenv.m_data+this->m_T,cval);
    std::fill(this->m_signal.m_data,this->m_signal.m_data+this->m_T,cval);
    std::fill(this->m_cenv_spec.m_data,this->m_cenv_spec.m_data+this->m_T,cval);
    std::fill(this->m_sig_samp.m_data,this->m_sig_samp.m_data+(this->m_T*this->m_N),cval);
    if(this->m_order == 1)
    {
        std::fill(this->m_cenv_corr.m_data,this->m_cenv_corr.m_data+(this->m_nomegs*this->m_nthets),cval);
        std::fill(this->m_ambig.m_data,this->m_ambig.m_data+(this->m_nomegs*this->m_nthets),cval);
        std::fill(this->m_mod_ambig.m_data,this->m_mod_ambig.m_data+(this->m_nomegs*this->m_nthets),rval);
        std::fill(this->m_cenv_corr_i.m_data,this->m_cenv_corr_i.m_data+(this->m_nomegs*this->m_nthets),rval);
        std::fill(this->m_cenv_corr_q.m_data,this->m_cenv_corr_q.m_data+(this->m_nomegs*this->m_nthets),rval);
    }
    else if(this->m_order == 2)
    {
        std::fill(this->m_cenv_corr.m_data,this->m_cenv_corr.m_data+(this->m_nthets*this->m_nomegs),cval);
        std::fill(this->m_ambig.m_data,this->m_ambig.m_data+(this->m_nthets*this->m_nomegs),cval);
        std::fill(this->m_mod_ambig.m_data,this->m_mod_ambig.m_data+(this->m_nthets*this->m_nomegs),rval);
        std::fill(this->m_cenv_corr_i.m_data,this->m_cenv_corr_i.m_data+(this->m_nthets*this->m_nomegs),rval);
        std::fill(this->m_cenv_corr_q.m_data,this->m_cenv_corr_q.m_data+(this->m_nthets*this->m_nomegs),rval);
    }
    std::fill(this->m_carrier_i.m_data,this->m_carrier_i.m_data+this->m_T,rval);
    std::fill(this->m_carrier_q.m_data,this->m_carrier_i.m_data+this->m_T,rval);
    std::fill(this->m_cenv_i.m_data,this->m_cenv_i.m_data+this->m_T,rval);
    std::fill(this->m_cenv_q.m_data,this->m_cenv_q.m_data+this->m_T,rval);
    std::fill(this->m_signal_i.m_data,this->m_signal_i.m_data+this->m_T,rval);
    std::fill(this->m_signal_q.m_data,this->m_signal_q.m_data+this->m_T,rval);
#endif 
    
}

void
gms::radiolocation
::AM_wb_signal_t
::ceate_signal_plot(const std::size_t n_samp,
                    const float * __restrict sig_arg,
                    const float * __restrict sig_val,
                    const std::string &header,
                    const std::string &title)
{
    std::string plot_fname;
    std::string sig_fname;
    std::ofstream plot_unit;
    std::ofstream sig_unit;
    sig_fname = header+"_plot.txt";
    sig_unit.open(sig_fname.c_str());
    for(std::size_t __i{0ull}; __i != n_samp; ++__i)
    {
        sig_unit.operator<<(" ").operator<<(sig_arg[__i]).operator<<(" ").operator<<(sig_val[__i]).operator<<("\n");
    }
    sig_unit.close();
    std::cout.operator<<("Created signal data file \"").operator<<(sig_fname.c_str()).operator<<("\".\n");
    plot_fname.operator=(header+"plot_commands.txt");
    plot_unit.open(plot_fname.c_str());
    plot_unit.operator<<("#").operator<<(plot_fname.c_str()).operator<<("\n");
    plot_unit.operator<<("#\n");
    plot_unit.operator<<("# Usage:\n");
    plot_unit.operator<<("# gnuplot < ").operator<<(plot_fname.c_str()).operator<<("\n");
    plot_unit.operator<<("#\n");
    plot_unit.operator<<("set term png\n");
    plot_unit.operator<<("set output \"").operator<<(header.c_str()).operator<<(".png\"\n");
    plot_unit.operator<<("set xlabel 't'\n");
    plot_unit.operator<<("set ylabel 'y(t)'\n");
    plot_unit.operator<<("set title '").operator<<(title.c_str()).operator<<("'\n");
    plot_unit.operator<<("set grid\n");
    plot_unit.operator<<("set style data lines\n");
    plot_unit.operator<<("plot \"").operator<<(plot_fname.c_str()).operator<<("\" using 1:2 lw 1 linecolor rgb \"blue\"\n");
    plot_unit.operator<<("quit\n");
    plot_unit.close();
    std::cout << " Created signal data file \"" << plot_fname <<("\"\n");
}





