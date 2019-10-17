
//
//	Implementation
//

gms::math::
Quaternion::Quaternion()
:
m_x{0.0E+00},
m_y{0.0E+00},
m_z{0.0e+00},
m_w{0.0E+00} {}

gms::math::
Quaternion::Quaternion(const double a)
:
m_x{ a },
m_y{ 0.0E+00 },
m_z{ 0.0E+00 },
m_w{ 0.0E+00 } {}

gms::math::
Quaternion::Quaternion(const double b,
		       const double c,
	               const double d)
:
m_x{ 0.0E+00 },
m_y{ b },
m_z{ c },
m_w{ d } {}

gms::math::
Quaternion::Quaternion(  const double x,
			 const double y,
			 const double z,
			 const double w)
:
m_x{ x },
m_y{ y },
m_z{ z },
m_w{ w }  {}

gms::math::
Quaternion::Quaternion(const std::complex<double> &c1,
			const std::complex<double> &c2)
:
m_x{ c1.real() },
m_y{ c1.imag() },
m_z{ c2.real() },
m_w{ c2.imag() } {}

gms::math::
Quaternion::Quaternion(const Quaternion &q)
:
m_x{ q.m_x },
m_y{ q.m_y },
m_z{ q.m_z },
m_w{ q.m_w }  {}

gms::math::Quaternion &
gms::math::Quaternion::operator=(const Quaternion &q) {
	if (this == &q) return (*this);
	m_x = q.m_x;
	m_y = q.m_y;
	m_z = q.m_z;
	m_w = q.m_w;
	return (*this);
}

gms::math::Quaternion &
gms::math::Quaternion::operator=(const double s) {
	m_x = s;
	m_y = 0.0E+00;
	m_z = 0.0E+00;
	m_w = 0.0E+00;
	return (*this);
}

gms::math::Quaternion &
gms::math::Quaternion::operator=(const C64 &c) {
	m_x = c.real();
	m_y = c.imag();
	m_z = 0.0E+00;
	m_w = 0.0E+00;
	return (*this);
}

std::ostream &
gms::math::operator<<(std::ostream &os,
		     const Quaternion &q) {
	os << std::fixed << std::setprecision(15) <<
		"x = " << "{" << q.get_x() << "}" <<
		"y = " << "{" << q.get_y() << "}" <<
		"z = " << "{" << q.get_z() << "}" <<
		"w = " << "{" << q.get_w() << "}" << "\n";
	return (os);
}

const double
gms::math::Quaternion::operator[](const uint32_t idx) const {
   
	return (reinterpret_cast<const double*>(this)[idx]);
}

gms::math::Quaternion
gms::math::operator+(const Quaternion &q1,
		     const Quaternion &q2) {
	return (Quaternion{ Q_x(q1) + Q_x(q2),
						Q_y(q1) + Q_y(q2),
						Q_z(q1) + Q_z(q2),
						Q_w(q1) + Q_w(q2) });
}

gms::math::Quaternion
gms::math::operator+(const Quaternion &q,
                     const std::complex<double> &c) {
	return (Quaternion{ Q_x(q) + c.real(),
						Q_y(q) + c.imag(),
						Q_z(q) + c.real(),
						Q_w(q) + c.imag() });
}

gms::math::Quaternion
gms::math::operator+(const Quaternion &q,
		     const double s) {
	return (Quaternion{ Q_x(q) + s,
						Q_y(q) + s,
						Q_z(q) + s,
						Q_w(q) + s });
}

gms::math::Quaternion
gms::math::operator+(const std::complex<double> &c,
		     const Quaternion &q) {
	return (Quaternion{ c.real() + Q_x(q),
						c.imag() + Q_y(q),
						c.real() + Q_z(q),
						c.imag() + Q_w(q) });
}

gms::math::Quaternion
gms::math::operator+(const double s,
		     const Quaternion &q) {
	return (Quaternion{ s + Q_x(q),
					    s + Q_y(q),
						s + Q_z(q),
						s + Q_w(q) });
}

gms::math::Quaternion
gms::math::operator+=(Quaternion &q1,
			const Quaternion &q2) {
	q1 = q1 + q2;
	return (q1);
}

gms::math::Quaternion
gms::math::operator+=(Quaternion &q,
					 const std::complex<double> &c) {
	q = q + c;
	return (q);
}

gms::math::Quaternion
gms::math::operator+=(const C64 &c,
			Quaternion &q) {
	q = c / q;
	return (q);
}

gms::math::Quaternion
gms::math::operator+=(Quaternion &q,
				 const double s) {
	q = q + s;
	return (q);
}

gms::math::Quaternion
gms::math::operator+=(const double s,
		      Quaternion &q) {
	q = s / q;
	return (q);
}

gms::math::Quaternion
gms::math::operator-(const Quaternion &q1,
		     const Quaternion &q2) {
	return (Quaternion{ Q_x(q1) - Q_x(q2),
					    Q_y(q1) - Q_y(q2),
						Q_z(q1) - Q_z(q2),
						Q_w(q1) - Q_w(q2) });
}

gms::math::Quaternion
gms::math::operator-(const Quaternion &q,
		     const C64 &c) {
	return (Quaternion{ Q_x(q) - c.real(),
						Q_y(q) - c.imag(),
						Q_z(q) - c.real(),
						Q_w(q) - c.imag() });
}

gms::math::Quaternion
gms::math::operator-(const Quaternion &q,
			 const double s) {
	return (Quaternion{ Q_x(q)-s,
						Q_y(q) - s,
						Q_z(q) - s,
						Q_w(q) - s });
}

gms::math::Quaternion
gms::math::operator-(const C64 &c,
		     const Quaternion &q) {
	return (Quaternion{ c.real() - Q_x(q),
						c.imag() - Q_y(q),
						c.real() - Q_z(q),
						c.imag() - Q_w(q) });
}

gms::math::Quaternion
gms::math::operator-(const double s,
		     const Quaternion &q) {
	return (Quaternion{ s - Q_x(q),
						s - Q_y(q),
						s - Q_z(q),
						s - Q_w(q) });
}

gms::math::Quaternion
gms::math::operator-=(Quaternion &q1,
		       const Quaternion &q2) {
	q1 = q1 - q2;
	return (q1);
}

gms::math::Quaternion
gms::math::operator-=(Quaternion &q,
		      const C64 &c) {
	q = q - c;
	return (q);
}

gms::math::Quaternion
gms::math::operator-=(const C64 &c,
		      Quaternion &q) {
	q = c - q;
	return (q);
}

gms::math::Quaternion
gms::math::operator-=(Quaternion &q,
		      const double s) {
	q = q - s;
	return (q);
}

gms::math::Quaternion
gms::math::operator-=(const double s,
		     Quaternion &q) {
	q = s - q;
	return (q);
}

gms::math::Quaternion
gms::math::operator*(const Quaternion &q1,
			const Quaternion &q2) {
	const double x = Q_x(q1)*Q_x(q2) - Q_y(q1)*Q_y(q2) - 
					 Q_z(q1)*Q_z(q2) - Q_w(q1)*Q_w(q2);
	const double y = Q_x(q1)*Q_y(q2) + Q_y(q1)*Q_x(q2) +
					 Q_z(q1)*Q_w(q2) - Q_w(q1)*Q_w(q2);
	const double z = Q_x(q1)*Q_z(q2) - Q_y(q1)*Q_w(q2) +
					 Q_z(q1)*Q_x(q2) + Q_w(q1)*Q_y(q2);
	const double w = Q_x(q1)*Q_z(q2) + Q_y(q1)*Q_z(q2) - 
					 Q_z(q1)*Q_y(q2) + Q_z(q1)*Q_x(q2);
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::operator*(const Quaternion &q,
		     const C64 &c) {
	const double x = Q_x(q)*c.real() - Q_y(q)*c.imag();
	const double y = Q_x(q)*c.imag() + Q_y(q)*c.real();
	const double z = Q_z(q)*c.real() + Q_w(q)*c.imag();
	const double w = -Q_z(q)*c.imag() + Q_w(q)*c.real();
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::operator*(const Quaternion &q,
			const double s) {
	return (Quaternion{ Q_x(q)*s,
						Q_y(q)*s,
						Q_z(q)*s,
						Q_w(q)*s });
}

gms::math::Quaternion
gms::math::operator*(const C64 &c,
			const Quaternion &q) {
	const double x = c.real()*Q_x(q) - c.imag()*Q_y(q);
	const double y = c.imag()*Q_x(q) + c.real()*Q_y(q);
	const double z = c.real()*Q_z(q) + c.imag()*Q_w(q);
	const double w = c.imag()*(-Q_z(q)) + c.real()*Q_w(q);
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::operator*(const double s,
		      const Quaternion &q) {
	return (Quaternion{ s*Q_x(q),
						s*Q_y(q),
						s*Q_z(q),
						s*Q_w(q) });
}

gms::math::Quaternion
gms::math::operator*=(Quaternion &q1,
		      const Quaternion &q2) {
	q1 = q1 * q2;
	return (q1);
}

gms::math::Quaternion
gms::math::operator*=(Quaternion &q,
			 const C64 &c) {
	q = q * c;
	return (q);
}

gms::math::Quaternion
gms::math::operator*=(const C64 &c,
			Quaternion &q) {
	q = c * q;
	return (q);
}

gms::math::Quaternion
gms::math::operator*=( Quaternion &q,
			 const double s) {
	q = q * s;
	return (q);
}

gms::math::Quaternion
gms::math::operator*=(const double s,
			 Quaternion &q) {
	q = s * q;
	return (q);
}

gms::math::Quaternion
gms::math::operator/(const Quaternion &q1,
			 const Quaternion &q2) {
	const double denom = Q_x(q2)*Q_x(q2) + Q_y(q2)*Q_y(q2) + 
						 Q_z(q2)*Q_z(q2) + Q_w(q2)*Q_w(q2);
	const double x = (Q_x(q1)*Q_x(q1) + Q_y(q1)*Q_y(q2)+
				      Q_z(q1)*Q_z(q2) + Q_w(q1)*Q_w(q2)) / denom;
	const double y = (-Q_x(q1)*Q_y(q2) + Q_y(q1)*Q_x(q2)-
					  Q_z(q1)*Q_w(q2) + Q_w(q1)*Q_z(q2)) / denom;
	const double z = (-Q_x(q1)*Q_z(q2) + Q_y(q1)*Q_w(q2)+
					  Q_z(q1)*Q_x(q2) - Q_w(q1)*Q_y(q2)) / denom;
	const double w = (-Q_x(q1)*Q_w(q2) - Q_y(q1)*Q_z(q2)+
					  Q_z(q1)*Q_y(q2) + Q_w(q1)*Q_x(q2)) / denom;
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::operator/(const Quaternion &q,
			 const C64 &c) {
	const double denom = c.real()*c.real() + 
						 c.imag()*c.imag();
	const double x = (Q_x(q)*c.real() + Q_y(q)*c.imag())/denom;
	const double y = (-Q_x(q)*c.imag() + Q_y(q)*c.real())/denom;
	const double z = (Q_z(q)*c.real() - Q_w(q)*c.imag())/denom;
	const double w = (Q_z(q)*c.imag() + Q_w(q)*c.real()) / denom;
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::operator/(const Quaternion &q,
			 const double s) {
	return (Quaternion{ Q_x(q) / s,
						Q_y(q) / s,
						Q_z(q) / s,
						Q_w(q) / s });
}

gms::math::Quaternion
gms::math::operator/(const C64 &c,
			 const Quaternion &q) {
	const double denom = c.real()*c.real() +
						 c.imag()*c.imag();
	const double x = (c.real()*Q_x(q) + Q_y(q)*c.imag())/denom;
	const double y = (c.imag()*(-Q_x(q)) + Q_y(q)*c.real())/denom;
	const double z = (c.real()*Q_z(q) - c.imag()*Q_w(q))/denom;
	const double w = (c.imag()*Q_z(q) + c.real()*Q_w(q)) / denom;
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::operator/(const double s,
		     const Quaternion &q) {
	return (Quaternion{ s / Q_x(q),
						s / Q_y(q),
						s / Q_z(q),
						s / Q_w(q) });
}

gms::math::Quaternion
gms::math::operator/=(Quaternion &q1,
		      const Quaternion &q2) {
	q1 = q1 / q2;
	return (q1);
}

gms::math::Quaternion
gms::math::operator/=(Quaternion &q,
		     const C64 &c) {
	q = q / c;
	return (q);
}

gms::math::Quaternion
gms::math::operator/=(const C64 &c,
		      Quaternion &q) {
	q = c / q;
	return (q);
}

gms::math::Quaternion
gms::math::operator/=(Quaternion &q,
			 const double s) {
	q = q / s;
	return (q);
}

gms::math::Quaternion
gms::math::operator/=(const double s,
		      Quaternion &q) {
	q = s / q;
	return (q);
}

bool
gms::math::operator==(const Quaternion &q1,
			 const Quaternion &q2) {
	return (Q_x(q1) == Q_x(q2)     &&
		    Q_y(q1) == Q_y(q2)     &&
		    Q_z(q1) == Q_z(q2)     &&
		    Q_w(q1) == Q_w(q2));
}

bool
gms::math::operator==(const Quaternion &q,
			 const C64 &c) {
	return (Q_x(q) == c.real()     &&
			Q_y(q) == c.imag()     &&
		    Q_z(q) == 0.0          &&
		    Q_w(q) == 0.0);
}

bool
gms::math::operator==(const Quaternion &q,
			 const double s) {
	return (Q_x(q) == s      &&
		    Q_y(q) == 0.0    &&
		    Q_z(q) == 0.0    &&
		    Q_w(q) == 0.0);
}

bool
gms::math::operator==(const C64 &c,
			const Quaternion &q) {
	return (c.real() == Q_x(q) &&
		    c.imag() == Q_y(q) &&
		    0.0 == Q_z(q)      &&
		    0.0 == Q_w(q));
}

bool
gms::math::operator==(const double s,
			 const Quaternion &q) {
	return (s == Q_x(q) &&
		    0.0 == Q_y(q) &&
		    0.0 == Q_z(q) &&
		    0.0 == Q_w(q));
}

bool
gms::math::operator!=(const Quaternion &q1,
			 const Quaternion &q2) {
	return (!(q1 == q2));
}

bool
gms::math::operator!=(const Quaternion &q1,
			 const C64 &c) {
	return (!(q1 == c));
}

bool
gms::math::operator!=(const Quaternion &q1,
		      const double s) {
	return(!(q1 == s));
}

bool
gms::math::operator!=(const C64 &c,
		      const Quaternion &q) {
	return (!(c == q));
}

bool
gms::math::operator!=(const double s,
		      const Quaternion &q) {
	return (!(s == q));
}

bool
gms::math::operator>(const Quaternion &q1,
		     const Quaternion &q2) {
	return (Q_x(q1) > Q_x(q2) &&
			Q_y(q1) > Q_y(q2) &&
			Q_z(q1) > Q_z(q2) &&
			Q_w(q1) > Q_w(q2));
}

bool
gms::math::operator<(const Quaternion &q1,
		     const Quaternion &q2) {
	return (Q_x(q1) < Q_x(q2) &&
			Q_y(q1) < Q_y(q2) &&
			Q_z(q1) < Q_z(q2) &&
			Q_w(q1) < Q_w(q2));
}

bool
gms::math::operator>=(const Quaternion &q1,
		      const Quaternion &q2) {
	return (Q_x(q1) >= Q_y(q2) &&
			Q_y(q1) >= Q_y(q2)    &&
			Q_z(q1) >= Q_z(q2)   &&
			Q_w(q1) >= Q_w(q2) );
}

bool
gms::math::operator<=(const Quaternion &q1,
		      const Quaternion &q2) {
	return (Q_x(q1) <= Q_y(q2) &&
		    Q_y(q1) <= Q_y(q2) &&
		    Q_z(q1) <= Q_z(q2) &&
		    Q_w(q1) <= Q_w(q2));
}

gms::math::Quaternion
gms::math::conjugate(const Quaternion &q) {
	return (Quaternion{ Q_x(q), -Q_y(q),
					   -Q_z(q), -Q_w(q) });
}

double
gms::math::norm(const Quaternion &q) {
	return (std::sqrt(Q_SQR(Q_x(q))+
					  Q_SQR(Q_y(q))+
					  Q_SQR(Q_z(q))+
					  Q_SQR(Q_w(q))));
}

double
gms::math::vnorm(const Quaternion &q) {
	return (std::sqrt(Q_SQR(Q_y(q))+
					  Q_SQR(Q_z(q))+
					  Q_SQR(Q_w(q))));
}

double
gms::math::distance(const Quaternion &q1,
		     const Quaternion &q2) {
	double dx = (Q_x(q1) - Q_x(q2)) * (Q_x(q1) - Q_x(q2));
	double dy = (Q_y(q1) - Q_y(q2)) * (Q_y(q1) - Q_y(q2));
	double dz = (Q_z(q1) - Q_z(q2)) * (Q_z(q1) - Q_z(q2));
	double dw = (Q_w(q1) - Q_w(q2)) * (Q_w(q1) - Q_w(q2));
	return (std::sqrt(dx+dy+dz+dw));
}

gms::math::Quaternion
gms::math::unit(const Quaternion &q) {
	double t = norm(q);
	if (t < std::numeric_limits<double>::epsilon())
		t = std::numeric_limits<double>::epsilon();
	return (q / t);
}

gms::math::Quaternion
gms::math::polar_decomp(const Quaternion &q) {
	Quaternion tq = unit(q);
	double t = norm(q);
	return (tq * t);
}

gms::math::Quaternion
gms::math::reciprocal(const Quaternion &q) {
	Quaternion tq = conjugate(q);
	double t = norm(q);
	return (tq / (t*t));
}

void
gms::math::mat4x4(const Quaternion &q,
		  double(&m4x4)[4][4],
		  const int32_t mtype) {
	
	if (mtype == 0) {
		m4x4[0][0] =  Q_x(q);
		m4x4[0][1] = -Q_y(q);
		m4x4[0][2] = -Q_z(q);
		m4x4[0][3] = -Q_w(q);
		m4x4[1][0] =  Q_y(q);
		m4x4[1][1] =  Q_x(q);
		m4x4[1][2] = -Q_w(q);
		m4x4[1][3] =  Q_z(q);
		m4x4[2][0] =  Q_z(q);
		m4x4[2][1] =  Q_w(q);
		m4x4[2][2] =  Q_x(q);
		m4x4[2][3] = -Q_y(q);
		m4x4[3][0] =  Q_w(q);
		m4x4[3][1] = -Q_z(q);
		m4x4[3][2] =  Q_y(q);
		m4x4[3][3] =  Q_x(q);
	}
	else if (mtype == 1) {
		m4x4[0][0] =  Q_x(q);
		m4x4[0][1] =  Q_y(q);
		m4x4[0][2] =  Q_z(q);
		m4x4[0][3] =  Q_w(q);
		m4x4[1][0] = -Q_y(q);
		m4x4[1][1] =  Q_x(q);
		m4x4[1][2] = -Q_w(q);
		m4x4[1][3] =  Q_z(q);
		m4x4[2][0] = -Q_z(q);
		m4x4[2][1] =  Q_w(q);
		m4x4[2][2] =  Q_x(q);
		m4x4[2][3] = -Q_y(q);
		m4x4[3][0] = -Q_w(q);
		m4x4[3][1] = -Q_z(q);
		m4x4[3][2] =  Q_y(q);
		m4x4[3][3] =  Q_x(q);
	}
}

gms::math::Quaternion
gms::math::exp(const Quaternion &q) {
	double t1 = std::exp(q.real());
	double t2 = vnorm(q);
	CHECK_QUATERNION_LEPS(t2)
	double s1 = q.get_y();
	double s2 = q.get_z();
	double s3 = q.get_w();
	s1 = s1 / t2 * std::sin(t2);
	s2 = s2 / t2 * std::sin(t2);
	s3 = s3 / t2 * std::sin(t2);
	return (Quaternion{ t1*std::cos(t2),
					    t1*s1,
						t1*s2,
						t1*s3 });
}

gms::math::Quaternion
gms::math::ln(const Quaternion &q) {
	double t1 = std::log(norm(q));
	CHECK_QUATERNION_LEPS(t1)
	double t2 = vnorm(q);
	double s1 = Q_y(q);
	double s2 = Q_z(q);
	double s3 = Q_w(q);
	s1 = s1 / t2 * std::acos(Q_x(q)/t2);
	s2 = s2 / t2 * std::acos(Q_x(q)/t2);
	s3 = s3 / t2 * std::acos(Q_x(q) / t2);
	return (Quaternion{t1, s1,s2,s3});
}

gms::math::Quaternion
gms::math::spherical(const double rho,
		     const double theta,
		     const double phi1,
		      const double phi2) {
	
	double ONE = 1.0;
	const double w = std::sin(phi2);
	ONE *= std::cos(phi2);
	double z = std::sin(phi1) * ONE;
	ONE *= std::cos(phi1);
	return (Quaternion{ std::cos(theta)*ONE,
						std::sin(theta)*ONE,
						z,w });
}

gms::math::Quaternion
gms::math::semipolar(const double rho,
		     const double alpha,
		     const double theta1,
		     const double theta2) {
	
	const double x = std::cos(alpha)*std::cos(theta1);
	const double y = std::cos(alpha)*std::sin(theta1);
	const double z = std::sin(alpha)*std::cos(theta2);
	const double w = std::sin(alpha)*std::sin(theta2);
	return (rho*Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::multipolar(const double rho1,
		      const double theta1,
		      const double rho2,
		      const double theta2) {
	
	const double x = rho1*std::cos(theta1);
	const double y = rho1*std::sin(theta1);
	const double z = rho2*std::cos(theta2);
	const double w = rho2*std::sin(theta2);
	return (Quaternion{x,y,z,w});
}

gms::math::Quaternion
gms::math::cylindrospherical(    const double t,
				 const double rad,
				 const double lat,
				 const double lon) {
	
	const double y = rad*std::cos(lon)*std::cos(lat);
	const double z = rad*std::sin(lon)*std::cos(lat);
	const double w = rad*std::sin(lat);
	return (Quaternion{t,y,z,w});
}

gms::math::Quaternion
gms::math::cylindrical( const double r,
			const double angle,
		        const double h1,
		        const double h2) {
	
	const double x = r*std::cos(angle);
	const double y = r*std::sin(angle);
	return (Quaternion{x,y,h1,h2});
}


