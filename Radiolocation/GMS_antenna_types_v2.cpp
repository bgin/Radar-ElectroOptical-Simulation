

gms::
radiolocation::
E_c4_t::E_c4_t(const int32_t npts) {
    
      m_npts = npts;
      m_ex   = vc4(m_npts);
      m_ey   = vc4(m_npts);
      m_ez   = vc4(m_npts);
}


gms::
radiolocation::
E_c4_t::E_c4_t(const E_c4_t &x) {

      m_npts = x.m_npts;
      m_ex   = x.m_ex;
      m_ey   = x.m_ey;
      m_ez   = x.m_ez;
}


gms::
radiolocation::
E_c4_t::E_c4_t(E_c4_t &&x) {

      m_npts  = std::move(x.m_npts);
      m_ex    = std::move(x.m_ex);
      m_ey    = std::move(x.m_ey);
      m_ez    = std::move(x.m_ez);
}


gms::
radiolocation::
E_c4_t &
gms::
radiolocation::
E_c4_t::operator=(const E_c4_t &x) {

      if(this == &x) return (*this);
      E_c4_t tmp(x);
      std::swap(*this,tmp);
      return (*this);
}


gms::
radiolocation::
E_c4_t &
gms::
radiolocation::
E_c4_t::operator=(E_c4_t &&x) {

      if(this == &x) return (*this);
      *this = std::move(x);
      return (*this);
}



