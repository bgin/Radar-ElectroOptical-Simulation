/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
file Polarization.h
class JonesVector
@aulthor: Bernard Gingold
@version:  1.0  19/10/2015

*/


// JonesVector inline function implementation
inline  std::complex<double>   gms::radiolocation::JonesVector::h()
{
	return this->m_h;
}

inline  const  std::complex<double>  gms::radiolocation::JonesVector::h() const
{
	return this->m_h;
}

inline   std::complex<double>   gms::radiolocation::JonesVector::v()
{
	return this->m_v;
}

inline  const  std::complex<double>  gms::radiolocation::JonesVector::v() const
{
	return this->m_v;
}


//  JonesMatrix inline functions implementation.
inline  std::complex<double>         gms::radiolocation::JonesMatrix::s0()
{
	return this->m_matrix[0];
}

inline  const   std::complex<double>   gms::radiolocation::JonesMatrix::s0() const
{
	return this->m_matrix[0];
}

inline  std::complex<double>         gms::radiolocation::JonesMatrix::s1()
{
	return this->m_matrix[1];
}

inline  const  std::complex<double>  gms::radiolocation::JonesMatrix::s1() const
{
	return this->m_matrix[1];
}

inline  std::complex<double>         gms::radiolocation::JonesMatrix::s2()
{
	return this->m_matrix[2];
}

inline  const  std::complex<double>  gms::radiolocation::JonesMatrix::s2() const
{
	return this->m_matrix[2];
}

inline  std::complex<double>          gms::radiolocation::JonesMatrix::s3()
{
	return this->m_matrix[3];
}

inline  const  std::complex<double>   gms::radiolocation::JonesMatrix::s3() const
{
	return this->m_matrix[3];
}

inline  std::complex<double>          *gms::radiolocation::JonesMatrix::matrix()
{
	return this->m_matrix;
}

inline  const  std::complex<double>    *gms::radiolocation::JonesMatrix::matrix() const
{
	return this->m_matrix;
}

inline   double    gms::radiolocation::StokesVector::s0() const
{
	return this->m_stokes[0];
}

inline   double    gms::radiolocation::StokesVector::s1() const
{
	return this->m_stokes[1];
}

inline   double    gms::radiolocation::StokesVector::s2() const
{
	return this->m_stokes[2];
}

inline   double    gms::radiolocation::StokesVector::s3() const
{
	return this->m_stokes[3];
}

inline const  double   *gms::radiolocation::StokesVector::stokes() const
{
	return this->m_stokes;
}

inline   double       gms::radiolocation::StokesVector::normalized_s1() const
{
	return s1() / s0();
}

inline   double       gms::radiolocation::StokesVector::normalized_s2() const
{
	return s2() / s0();
}

inline   double       gms::radiolocation::StokesVector::normalized_s3() const
{
	return s3() / s0();
}

inline   double      gms::radiolocation::StokesVector::dot(_In_ StokesVector const& lhs) const
{
	return (this->m_stokes[0] * lhs.m_stokes[0]) + (this->m_stokes[1] * lhs.m_stokes[1]) + (this->m_stokes[2] * lhs.m_stokes[2]) +
		(this->m_stokes[3] * lhs.m_stokes[3]);
}
