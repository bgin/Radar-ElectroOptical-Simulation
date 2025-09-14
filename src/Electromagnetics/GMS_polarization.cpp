
/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
file Polarization.h
class JonesVector
@aulthor: Bernard Gingold
@version:  1.0  19/10/2015

*/


#include "GMS_polarization.h"

radiolocation::JonesVector::JonesVector( std::complex<double> const h,  std::complex<double> const v) :
m_h(h), m_v(v)
{

}

radiolocation::JonesVector::JonesVector( struct JonesVectorParams const& params) : m_h(params.ps0), m_v(params.ps1)
{
	
}

radiolocation::JonesVector::JonesVector( JonesVector const& rhs) :
m_h(rhs.m_h), m_v(rhs.m_v)
{

}

radiolocation::JonesVector::JonesVector( JonesVector &&rhs) :
m_h(std::move(rhs.m_h)), m_v(std::move(rhs.m_v))
{

}

radiolocation::JonesVector&        radiolocation::JonesVector::operator=( JonesVector const& rhs)
{
	this->m_h = rhs.m_h;
	this->m_v = rhs.m_v;

	return *this;
}

radiolocation::JonesVector&         radiolocation::JonesVector::operator=( JonesVector &&rhs)
{
	this->m_h = std::move(rhs.m_h);
	this->m_v = std::move(rhs.m_v);

	return *this;
}

std::ostream&                      radiolocation::operator<<( std::ostream &os,  JonesVector const& rhs)
{
	os.scientific;
	os << "Ex=" << rhs.m_h.imag() << "e^iph=" << rhs.m_h.real() << "Ey=" << rhs.m_v.imag() << "e^iph=" << rhs.m_v.real() << std::endl;
	return os;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator*=( JonesVector const& rhs)
{
	this->m_h.operator*=(rhs.m_h);
	this->m_v.operator*=(rhs.m_v);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator*=( std::complex<double> const& c)
{
	this->m_h.operator*=(c);
	this->m_v.operator*=(c);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator/=( JonesVector const& rhs)
{
	this->m_h.operator/=(rhs.m_h);
	this->m_v.operator/=(rhs.m_v);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator/=( std::complex<double> const& c)
{
	this->m_h.operator/=(c);
	this->m_v.operator/=(c);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator-=( JonesVector const& rhs)
{
	this->m_h.operator-=(rhs.m_h);
	this->m_v.operator-=(rhs.m_v);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator+=( JonesVector const& rhs)
{
	this->m_h.operator+=(rhs.m_h);
	this->m_v.operator+=(rhs.m_v);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator+=( std::complex<double> const& c)
{
	this->m_h.operator+=(c);
	this->m_v.operator+=(c);

	return *this;
}

radiolocation::JonesVector&        radiolocation::JonesVector::operator-=( std::complex<double> const& c)
{
	this->m_h.operator-=(c);
	this->m_v.operator-=(c);

	return *this;
}

radiolocation::JonesVector         radiolocation::JonesVector::operator-() const
{
	return JonesVector(-m_h, -m_v);
}

                            radiolocation::JonesVector::operator==( JonesVector const& lhs) const


radiolocation::JonesVector         radiolocation::operator*( JonesVector const& lhs,  JonesVector const& rhs)
{
	return JonesVector(lhs.m_h * rhs.m_h, lhs.m_v * rhs.m_v);
}

radiolocation::JonesVector         radiolocation::operator*( JonesVector const& lhs,  std::complex<double> const& rhs)
{
	return JonesVector(lhs.m_h * rhs, lhs.m_v * rhs);
}

std::complex<double>               radiolocation::JonesVector::operator*( JonesVector const& lhs)
{
	
	return ((this->m_h * lhs.m_h) + (this->m_v * lhs.m_v));
}

radiolocation::JonesVector         radiolocation::operator/( JonesVector const& lhs,  JonesVector const& rhs)
{
	return JonesVector(lhs.m_h / rhs.m_h, lhs.m_v / rhs.m_v);
}

radiolocation::JonesVector         radiolocation::operator/( JonesVector const& lhs,  std::complex<double> const& rhs)
{
	return JonesVector(lhs.m_h / rhs, lhs.m_v / rhs);
}

radiolocation::JonesVector          radiolocation::operator-( JonesVector const& lhs,  JonesVector const& rhs)
{
	return  JonesVector(lhs.m_h - rhs.m_h, lhs.m_v - rhs.m_v);

}

radiolocation::JonesVector           radiolocation::operator-( JonesVector const& lhs,  std::complex<double> const& rhs)
{
	return JonesVector(lhs.m_h - rhs, lhs.m_v - rhs);
}



radiolocation::JonesVector          radiolocation::operator+( JonesVector const& lhs,  JonesVector const& rhs)
{
	return JonesVector(lhs.m_h + rhs.m_h, lhs.m_v + rhs.m_v);
}

radiolocation::JonesVector          radiolocation::operator+( JonesVector const& lhs,  std::complex<double> const& rhs)
{
	return JonesVector(lhs.m_h + rhs, lhs.m_v + rhs);
}

double            radiolocation::JonesVector::field_intensity() const
{
	return std::norm(this->m_h + std::norm(this->m_v));
}

double            radiolocation::JonesVector::field_atan() const
{
	return std::atan(std::abs(m_v) / std::abs(m_h));
	
}

double            radiolocation::JonesVector::field_phase_diff() const
{
	return std::arg(this->m_v) - std::arg(this->m_h);
}

double            radiolocation::JonesVector::degree_polarization() const
{
	return 1.0;
}


//  JonesMatrix Implementation.

radiolocation::JonesMatrix::JonesMatrix( std::complex<double> const& _s0,  std::complex<double> const& _s1,  std::complex<double> const& _s2,
	 std::complex<double> const& _s3)
{
	this->m_matrix[0].operator=( _s0);
	this->m_matrix[1].operator=( _s1);
	this->m_matrix[2].operator=(_s2);
	this->m_matrix[3].operator=(_s3);
}


radiolocation::JonesMatrix::JonesMatrix( struct JonesMatrixParams const& params)
{
	this->m_matrix[0].operator=(params.ps0);
	this->m_matrix[1].operator=(params.ps1);
	this->m_matrix[2].operator=(params.ps2);
	this->m_matrix[3].operator=(params.ps3);
}

radiolocation::JonesMatrix::JonesMatrix( JonesMatrix const& rhs)
{
	this->m_matrix[0].operator=(rhs.m_matrix[0]);
	this->m_matrix[1].operator=(rhs.m_matrix[1]);
	this->m_matrix[2].operator=(rhs.m_matrix[2]);
	this->m_matrix[3].operator=(rhs.m_matrix[3]);
}

radiolocation::JonesMatrix::JonesMatrix( JonesMatrix &&rhs)
{
	this->m_matrix[0].operator=(std::move(rhs.m_matrix[0]));
	this->m_matrix[1].operator=(std::move(rhs.m_matrix[1]));
	this->m_matrix[2].operator=(std::move(rhs.m_matrix[2]));
	this->m_matrix[3].operator=(std::move(rhs.m_matrix[3]));
}

radiolocation::JonesMatrix&        radiolocation::JonesMatrix::operator=( JonesMatrix const& rhs)
{
	this->m_matrix[0].operator=(rhs.m_matrix[0]);
	this->m_matrix[1].operator=(rhs.m_matrix[1]);
	this->m_matrix[2].operator=(rhs.m_matrix[2]);
	this->m_matrix[3].operator=(rhs.m_matrix[3]);

	return *this;
}

radiolocation::JonesMatrix&        radiolocation::JonesMatrix::operator=( JonesMatrix &&rhs)
{
	this->m_matrix[0].operator=(std::move(rhs.m_matrix[0]));
	this->m_matrix[1].operator=(std::move(rhs.m_matrix[1]));
	this->m_matrix[2].operator=(std::move(rhs.m_matrix[2]));
	this->m_matrix[3].operator=(std::move(rhs.m_matrix[3]));

	return *this;
}

radiolocation::JonesMatrix&        radiolocation::JonesMatrix::operator+=( JonesMatrix const& rhs)
{
	this->m_matrix[0].operator+=(rhs.m_matrix[0]);
	this->m_matrix[1].operator+=(rhs.m_matrix[1]);
	this->m_matrix[2].operator+=(rhs.m_matrix[2]);
	this->m_matrix[3].operator+=(rhs.m_matrix[3]);

	return *this;
}

radiolocation::JonesMatrix&        radiolocation::JonesMatrix::operator+=( std::complex<double> const& c)
{
	this->m_matrix[0].operator+=(c);
	this->m_matrix[1].operator+=(c);
	this->m_matrix[2].operator+=(c);
	this->m_matrix[3].operator+=(c);

	return *this;
}

radiolocation::JonesMatrix&        radiolocation::JonesMatrix::operator-=( JonesMatrix const& rhs)
{
	this->m_matrix[0].operator-=(rhs.m_matrix[0]);
	this->m_matrix[1].operator-=(rhs.m_matrix[1]);
	this->m_matrix[2].operator-=(rhs.m_matrix[2]);
	this->m_matrix[3].operator-=(rhs.m_matrix[3]);

	return *this;
}

radiolocation::JonesMatrix&        radiolocation::JonesMatrix::operator-=( std::complex<double> const& c)
{
	this->m_matrix[0].operator-=(c);
	this->m_matrix[1].operator-=(c);
	this->m_matrix[2].operator-=(c);
	this->m_matrix[3].operator-=(c);

	return *this;
}

radiolocation::JonesMatrix&         radiolocation::JonesMatrix::operator*=( JonesMatrix const& rhs)
{
	this->m_matrix[0].operator*=(rhs.m_matrix[0]);
	this->m_matrix[1].operator*=(rhs.m_matrix[1]);
	this->m_matrix[2].operator*=(rhs.m_matrix[2]);
	this->m_matrix[3].operator*=(rhs.m_matrix[3]);

	return *this;
}

radiolocation::JonesMatrix&         radiolocation::JonesMatrix::operator*=(  std::complex<double> const& c)
{
	this->m_matrix[0].operator*=(c);
	this->m_matrix[1].operator*=(c);
	this->m_matrix[2].operator*=(c);
	this->m_matrix[3].operator*=(c);

	return *this;
}

radiolocation::JonesMatrix&          radiolocation::JonesMatrix::operator/=( JonesMatrix const& rhs)
{
	this->m_matrix[0].operator/=(rhs.m_matrix[0]);
	this->m_matrix[1].operator/=(rhs.m_matrix[1]);
	this->m_matrix[2].operator/=(rhs.m_matrix[2]);
	this->m_matrix[3].operator/=(rhs.m_matrix[3]);

	return *this;
}

radiolocation::JonesMatrix&          radiolocation::JonesMatrix::operator/=( std::complex<double> const& c)
{
	this->m_matrix[0].operator/=(c);
	this->m_matrix[1].operator/=(c);
	this->m_matrix[2].operator/=(c);
	this->m_matrix[3].operator/=(c);

	return *this;
}

radiolocation::JonesMatrix&          radiolocation::JonesMatrix::operator+=( JonesMatrix &&rhs)
{
	this->m_matrix[0].operator+=(std::move(rhs.m_matrix[0]));
	this->m_matrix[1].operator+=(std::move(rhs.m_matrix[1]));
	this->m_matrix[2].operator+=(std::move(rhs.m_matrix[2]));
	this->m_matrix[3].operator+=(std::move(rhs.m_matrix[3]));

	return *this;
}

radiolocation::JonesMatrix&          radiolocation::JonesMatrix::operator-=(  JonesMatrix &&rhs)
{
	this->m_matrix[0].operator-=(std::move(rhs.m_matrix[0]));
	this->m_matrix[1].operator-=(std::move(rhs.m_matrix[1]));
	this->m_matrix[2].operator-=(std::move(rhs.m_matrix[2]));
	this->m_matrix[3].operator-=(std::move(rhs.m_matrix[3]));

	return *this;
}

radiolocation::JonesMatrix&          radiolocation::JonesMatrix::operator*=(  JonesMatrix &&rhs)
{
	this->m_matrix[0].operator*=(std::move(rhs.m_matrix[0]));
	this->m_matrix[1].operator*=(std::move(rhs.m_matrix[1]));
	this->m_matrix[2].operator*=(std::move(rhs.m_matrix[2]));
	this->m_matrix[3].operator*=(std::move(rhs.m_matrix[3]));

	return *this;
}

radiolocation::JonesMatrix&          radiolocation::JonesMatrix::operator/=(  JonesMatrix &&rhs)
{
	this->m_matrix[0].operator/=(std::move(rhs.m_matrix[0]));
	this->m_matrix[1].operator/=(std::move(rhs.m_matrix[1]));
	this->m_matrix[2].operator/=(std::move(rhs.m_matrix[2]));
	this->m_matrix[3].operator/=(std::move(rhs.m_matrix[3]));

	return *this;
}

std::complex<double>                 radiolocation::JonesMatrix::operator[]( const int index)
{
	_ASSERT((0 <= index) && (3 <= index));
	return this->m_matrix[index];
}

const  std::complex<double>          radiolocation::JonesMatrix::operator[]( const int index) const
{
	_ASSERT((0 <= index) && (3 <= index));
	return this->m_matrix[index];
}





radiolocation::JonesMatrix           radiolocation::JonesMatrix::operator-() const
{
	return JonesMatrix(-m_matrix[0], -m_matrix[1], -m_matrix[2], -m_matrix[3]);
}

radiolocation::JonesMatrix           radiolocation::operator+( JonesMatrix const& lhs,  JonesMatrix const& rhs)
{
	return JonesMatrix(lhs.m_matrix[0] + rhs.m_matrix[0], lhs.m_matrix[1] + rhs.m_matrix[1],
		lhs.m_matrix[2] + rhs.m_matrix[2], lhs.m_matrix[3] + rhs.m_matrix[3]);
}

radiolocation::JonesMatrix           radiolocation::operator+( JonesMatrix const& lhs, std::complex<double> const& c)
{
	return JonesMatrix(lhs.m_matrix[0] + c, lhs.m_matrix[1] + c, lhs.m_matrix[2] + c, lhs.m_matrix[3] + c);
}

radiolocation::JonesMatrix           radiolocation::operator-( JonesMatrix const& lhs,  JonesMatrix const& rhs)
{
	return JonesMatrix(lhs.m_matrix[0] - rhs.m_matrix[0], lhs.m_matrix[1] - rhs.m_matrix[1],
		lhs.m_matrix[2] - rhs.m_matrix[2], lhs.m_matrix[3] - rhs.m_matrix[3]);
}

radiolocation::JonesMatrix           radiolocation::operator-( JonesMatrix const& lhs,  std::complex<double> const& c)
{
	return JonesMatrix(lhs.m_matrix[0] - c, lhs.m_matrix[1] - c, lhs.m_matrix[2] - c, lhs.m_matrix[3] - c);
}

radiolocation::JonesMatrix           radiolocation::operator*( JonesMatrix const& lhs,  JonesMatrix const& rhs)
{
	return JonesMatrix(lhs.m_matrix[0] * rhs.m_matrix[0], lhs.m_matrix[1] * rhs.m_matrix[1],
		lhs.m_matrix[2] * rhs.m_matrix[2], lhs.m_matrix[3] * rhs.m_matrix[3]);
}

radiolocation::JonesMatrix           radiolocation::operator*( JonesMatrix const& lhs,  std::complex<double> const& c)
{
	return JonesMatrix(lhs.m_matrix[0] * c, lhs.m_matrix[1] * c, lhs.m_matrix[2] * c, lhs.m_matrix[3] * c);
}

radiolocation::JonesMatrix           radiolocation::operator/( JonesMatrix const& lhs,  JonesMatrix const& rhs)
{
	return JonesMatrix(lhs.m_matrix[0] / rhs.m_matrix[0], lhs.m_matrix[1] / rhs.m_matrix[1],
		lhs.m_matrix[2] / rhs.m_matrix[2], lhs.m_matrix[3] / rhs.m_matrix[3]);
}

radiolocation::JonesMatrix           radiolocation::operator/( JonesMatrix const& lhs,  std::complex<double> const& c)
{
	return JonesMatrix(lhs.m_matrix[0] / c, lhs.m_matrix[1] / c, lhs.m_matrix[2] / c, lhs.m_matrix[3] / c);
}

radiolocation::JonesMatrix         radiolocation::JonesMatrix::matrix_transpose() const
{
	return JonesMatrix(this->m_matrix[0], this->m_matrix[1], this->m_matrix[3], this->m_matrix[2]);
}

radiolocation::JonesMatrix         radiolocation::JonesMatrix::matrix_hermitian() const
{
	return JonesMatrix(this->m_matrix[0], this->m_matrix[1], std::conj(this->m_matrix[3]), std::conj(this->m_matrix[2]));
}

std::ostream&                     radiolocation::operator<<( std::ostream& os,  JonesMatrix const& rhs)
{
	os.scientific;
	os << "s0=" << rhs.m_matrix[0] << "s1=" << rhs.m_matrix[1] << "s2=" << rhs.m_matrix[2] << "s3=" << rhs.m_matrix[3] << std::endl;
	return os;
}




/*
@Brief: Stokes Vector implementation
*/

radiolocation::StokesVector::StokesVector( const double s0,  const double s1,  const double s2,  const double s3)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_set_pd(s0, s1, s2, s3));
}

radiolocation::StokesVector::StokesVector( struct StokesVectorParams const& params)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_set_pd(params.ps0, params.ps1, params.ps2, params.ps3));
}

radiolocation::StokesVector::StokesVector( __m256d const& v)
{
	_mm256_storeu_pd(&this->m_stokes[0], v);
}

radiolocation::StokesVector::StokesVector( StokesVector const& rhs) 
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_loadu_pd(&rhs.m_stokes[0]));
}

radiolocation::StokesVector::StokesVector( StokesVector &&rhs)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_loadu_pd(std::move(rhs.m_stokes)));
}

radiolocation::StokesVector&        radiolocation::StokesVector::operator=( StokesVector const& lhs)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_loadu_pd(&lhs.m_stokes[0]));
	return *this;
}

radiolocation::StokesVector&        radiolocation::StokesVector::operator=( StokesVector &&lhs)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_loadu_pd(std::move(lhs.m_stokes)));
	return *this;
}

std::ostream&                       radiolocation::operator<<( std::ostream& os,  StokesVector const& rhs)
{
	os.scientific;
	os << "I:" << rhs.m_stokes[0] << "Q:" << rhs.m_stokes[1] << "U:" << rhs.m_stokes[2] << "V:" << rhs.m_stokes[3] << std::endl;
	return os;
}

double                               radiolocation::StokesVector::operator[]( const int index)
{
	_ASSERT((index >= 0) && (index <= 3));
	return this->m_stokes[index];
}

const  double                        radiolocation::StokesVector::operator[]( const int index) const
{
	_ASSERT((index >= 0) && (index <= 3));
	return this->m_stokes[index];
}

double                                radiolocation::StokesVector::poa() const
{
	return std::atan2(this->m_stokes[3], this->m_stokes[0]) * 0.5;
}

double                                radiolocation::StokesVector::dcp() const
{
	return this->m_stokes[3] / this->m_stokes[0];
}

double                                radiolocation::StokesVector::polarization_ellipticity() const
{
	return m_stokes[3] / (m_stokes[0] + std::sqrt((m_stokes[1] * m_stokes[1]) + (m_stokes[2] * m_stokes[2])));
}

double                                radiolocation::StokesVector::polarization_eccentrity() const
{
	return std::sqrt(1.0 - (polarization_ellipticity() * polarization_ellipticity()));
}

radiolocation::StokesVector           radiolocation::StokesVector::zero_stokes_vector()
{
	return StokesVector(0.0, 0.0, 0.0, 0.0);
}

radiolocation::StokesVector           radiolocation::StokesVector::nonpolarized_stokes_vector()
{
	return StokesVector(1.0, 0.0, 0.0, 0.0);
}


radiolocation::StokesVector&         radiolocation::StokesVector::operator+=( StokesVector const& lhs)
{
	
	
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_add_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_loadu_pd(&lhs.m_stokes[0])));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator+=( const double s)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_add_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_set1_pd(s)));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator-=( StokesVector const& lhs)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_sub_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_loadu_pd(&lhs.m_stokes[0])));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator-=( const double s)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_sub_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_set1_pd(s)));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator*=( StokesVector const& lhs)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_mul_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_loadu_pd(&lhs.m_stokes[0])));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator*=( const double s)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_mul_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_set1_pd(s)));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator/=( StokesVector const& lhs)
{
	
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_div_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_loadu_pd(&lhs.m_stokes[0])));
	return *this;
}

radiolocation::StokesVector&         radiolocation::StokesVector::operator/=( const double s)
{
	_mm256_storeu_pd(&this->m_stokes[0], _mm256_div_pd(_mm256_loadu_pd(&this->m_stokes[0]), _mm256_set1_pd(s)));
	return *this;
}




radiolocation::StokesVector           radiolocation::StokesVector::operator-() const
{
	return StokesVector(-this->m_stokes[0], -this->m_stokes[1], -this->m_stokes[2], -this->m_stokes[3]);
}

radiolocation::StokesVector           radiolocation::operator+( StokesVector const& lhs,  StokesVector const& rhs)
{
	return  StokesVector(_mm256_add_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_loadu_pd(&rhs.m_stokes[0])));
}



radiolocation::StokesVector           radiolocation::operator+( StokesVector const& lhs,  const double rhs)
{
	return StokesVector(_mm256_add_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_set1_pd(rhs)));
}

radiolocation::StokesVector           radiolocation::operator-( StokesVector const& lhs,  StokesVector const& rhs)
{
	return  StokesVector(_mm256_sub_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_loadu_pd(&rhs.m_stokes[0])));
}

radiolocation::StokesVector           radiolocation::operator-( StokesVector const& lhs,  const double rhs)
{
	return StokesVector(_mm256_sub_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_set1_pd(rhs)));
}

radiolocation::StokesVector           radiolocation::operator*( StokesVector const& lhs,  StokesVector const& rhs)
{
	return  StokesVector(_mm256_mul_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_loadu_pd(&rhs.m_stokes[0])));
}

radiolocation::StokesVector           radiolocation::operator*( StokesVector const& lhs,  const double rhs)
{
	return StokesVector(_mm256_mul_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_set1_pd(rhs)));
}

radiolocation::StokesVector           radiolocation::operator/( StokesVector const& lhs,  StokesVector const& rhs)
{
	return StokesVector(_mm256_div_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_loadu_pd(&rhs.m_stokes[0])));
}

radiolocation::StokesVector           radiolocation::operator/( StokesVector const& lhs,  const double rhs)
{
	return StokesVector(_mm256_div_pd(_mm256_loadu_pd(&lhs.m_stokes[0]), _mm256_set1_pd(rhs)));
}

/*radiolocation::StokesVector         radiolocation::operator+( mathlib::VecF64AVX const& lhs,  mathlib::VecF64AVX const& rhs)
{
	return StokesVector(lhs + rhs);
}

radiolocation::StokesVector         radiolocation::operator-( mathlib::VecF64AVX const& lhs,   mathlib::VecF64AVX const& rhs)
{
	return StokesVector(lhs - rhs);
}

radiolocation::StokesVector         radiolocation::operator*( mathlib::VecF64AVX const& lhs,   mathlib::VecF64AVX const& rhs)
{
	return StokesVector(lhs * rhs);
}

radiolocation::StokesVector         radiolocation::operator/( mathlib::VecF64AVX const& lhs,   mathlib::VecF64AVX const& rhs)
{
	return StokesVector(lhs / rhs);
}*/



