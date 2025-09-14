
//
//	Implementation
//

gms::math::
Vector3D::Vector3D()
:
m_x{},
m_y{},
m_z{} {}

gms::math::
Vector3D::Vector3D(const double x,
		 const double y,
		 const double z)
:
m_x{ x },
m_y{ y },
m_z{ z }   {}

gms::math::
Vector3D::Vector3D(const double x)
:
m_x{ x },
m_y{},
m_z{} {}

gms::math::
Vector3D::Vector3D(const Vector3D &v)
:
m_x{ v.m_x },
m_y{ v.m_y },
m_z{ v.m_z } {}

gms::math::Vector3D &
gms::math::Vector3D::operator=(const Vector3D &v) {
	if (this == &v) return (*this);
	m_x = v.m_x;
	m_y = v.m_y;
	m_z = v.m_z;
	return (*this);
}

const double
gms::math::Vector3D::operator[](const uint32_t idx) const {

	return (reinterpret_cast<const double*>(&this->m_x)[idx]);
}

std::ostream &
gms::math::operator<<(std::ostream &os,
	              const Vector3D &v) {
	os << std::fixed << std::scientific <<
		  std::setprecision(15) <<
		  "Component x= " << V3D_X(v) <<
		  "Component y= " << V3D_Y(v) <<
		  "Component z= " << V3D_Z(v) << "\n";
	return (os);
}

gms::math::Vector3D
gms::math::operator+(const Vector3D &v1,
		     const Vector3D &v2) {
	return (Vector3D{ V3D_X(v1) + V3D_X(v2),
					  V3D_Y(v1) + V3D_Y(v2),
					  V3D_Z(v1) + V3D_Z(v2)});
}

gms::math::Vector3D
gms::math::operator+(const Vector3D &v,
		     const double s) {
	return (Vector3D{ V3D_X(v) + s,
					  V3D_Y(v) + s,
					  V3D_Z(v) + s });
}

gms::math::Vector3D
gms::math::operator+(const double s,
		     const Vector3D &v) {
	return (Vector3D{ s + V3D_X(v),
		              s + V3D_Y(v),
		              s + V3D_Z(v) });
}

gms::math::Vector3D
gms::math::operator+=(Vector3D &v1,
		      const Vector3D &v2) {
	v1 = v1 + v2;
	return (v1);
}

gms::math::Vector3D
gms::math::operator+=(Vector3D &v,
		      const double s) {
	v = v + s;
	return (v);
}

gms::math::Vector3D
gms::math::operator+=(const double s,
		Vector3D &v) {
	v = s + v;
	return (v);
}

gms::math::Vector3D
gms::math::operator-(const Vector3D &v1,
		 const Vector3D &v2) {
	return (Vector3D{ V3D_X(v1) - V3D_X(v2),
					  V3D_Y(v1) - V3D_Y(v2),
					  V3D_Z(v1) - V3D_Z(v2) });
}

gms::math::Vector3D
gms::math::operator-(const Vector3D &v,
		 const double s) {
	return (Vector3D{ V3D_X(v) - s,
		              V3D_Y(v) - s ,
		              V3D_Z(v) - s});
}

gms::math::Vector3D
gms::math::operator-(const double s,
		 const Vector3D &v) {
	return (Vector3D{ s - V3D_X(v),
					  s - V3D_Y(v),
					  s - V3D_Z(v) });
}

gms::math::Vector3D
gms::math::operator-=(Vector3D &v1,
			 const Vector3D &v2) {
	v1 = v1 - v2;
	return (v1);
}

gms::math::Vector3D
gms::math::operator-=(Vector3D &v1,
		const double s) {
	v1 = v1 - s;
	return (v1);
}

gms::math::Vector3D
gms::math::operator-=(const double s,
		Vector3D &v) {
	v = s - v;
	return (v);
}

gms::math::Vector3D
gms::math::operator*(const Vector3D &v1,
		const Vector3D &v2) {
	return (Vector3D{ V3D_X(v1) * V3D_X(v2),
					  V3D_Y(v1) * V3D_Y(v2),
					  V3D_Z(v1) * V3D_Z(v2) });
}

gms::math::Vector3D
gms::math::operator*(const Vector3D &v,
		     const double s) {
	return (Vector3D{ V3D_X(v) * s,
					  V3D_Y(v) * s,
					  V3D_Z(v) * s });
}

gms::math::Vector3D
gms::math::operator*(const double s,
			 const Vector3D &v) {
	return (Vector3D{ s * V3D_X(v),
					  s * V3D_Y(v),
					  s * V3D_Z(v) });
}

gms::math::Vector3D
gms::math::operator*=(Vector3D &v1,
		const Vector3D &v2) {
	v1 = v1 * v2;
	return (v1);
}

gms::math::Vector3D
gms::math::operator*=(Vector3D &v,
			const double s) {
	v = v * s;
	return (v);
}

gms::math::Vector3D
gms::math::operator*=(const double s,
		Vector3D &v) {
	v = s * v;
	return (v);
}

gms::math::Vector3D
gms::math::operator/(const Vector3D &v1,
		const Vector3D &v2) {
	return (Vector3D{V3D_X(v1) / V3D_X(v2),
					 V3D_Y(v1) / V3D_Y(v2),
					 V3D_Z(v1) / V3D_Z(v2)});
}

gms::math::Vector3D
gms::math::operator/(const Vector3D &v,
		const double s) {
	return (Vector3D{V3D_X(v) / s,
					 V3D_Y(v) / s,
					 V3D_Z(v) / s});
}

gms::math::Vector3D
gms::math::operator/(const double s,
		const Vector3D &v) {
	return (Vector3D{ s / V3D_X(v),
					  s / V3D_Y(v),
					  s / V3D_Z(v) });
}

gms::math::Vector3D
gms::math::operator/=(Vector3D &v1,
			 const Vector3D &v2) {
	v1 = v1 / v2;
	return (v1);
}

gms::math::Vector3D
gms::math::operator/=(Vector3D &v,
		const double s) {
	v = v / s;
	return (v);
}

gms::math::Vector3D
gms::math::operator/=(const double s,
			Vector3D &v) {
	v = s / v;
	return (v);
}

bool
gms::math::operator==(const Vector3D &v1,
		      const Vector3D &v2) {
	return (V3D_X(v1) == V3D_X(v2) &&
			V3D_Y(v1) == V3D_Y(v2) &&
			V3D_Z(v1) == V3D_Z(v2));
}

bool
gms::math::operator==(const Vector3D &v,
		      const double s) {
	return (V3D_X(v) == s &&
			V3D_Y(v) == 0.0E+00 &&
			V3D_Z(v) == 0.0E+00);
}

bool
gms::math::operator==(const double s,
		     const Vector3D &v) {
	return (s == V3D_X(v) &&
		0.0E+00 == V3D_Y(v) &&
		0.0E+00 == V3D_Z(v));
}

bool
gms::math::operator!=(const Vector3D &v1,
		      const Vector3D &v2) {
	
	return (!(v1 == v2));
}

bool
gms::math::operator!=(const Vector3D &v,
		      const double s) {
	return (!(v == s));
}

bool
gms::math::operator!=(const double s,
		      const Vector3D &v) {
	return (!(s == v));
}

bool
gms::math::operator>(const Vector3D &v1,
		 const Vector3D &v2) {
	return (V3D_X(v1) > V3D_X(v2) &&
			V3D_Y(v1) > V3D_Y(v2) &&
			V3D_Z(v1) > V3D_Z(v2));
}

bool
gms::math::operator>(const Vector3D &v1,
		     const double s) {
	return (V3D_X(v1) > s &&
			V3D_Y(v1) > 0.0E+00 &&
			V3D_Z(v1) > 0.0E+00);
}

bool
gms::math::operator>(const double s,
		     const Vector3D &v) {
	return (s > V3D_X(v) &&
			0.0E+00 > V3D_Y(v) &&
			0.0E+00 > V3D_Z(v));
}

bool
gms::math::operator<(const Vector3D &v1,
		     const Vector3D &v2) {
	return (V3D_X(v1) < V3D_X(v2) &&
			V3D_Y(v1) < V3D_Y(v2) &&
			V3D_Z(v1) < V3D_Z(v2));
}

bool
gms::math::operator<(const Vector3D &v,
		     const double s) {
	return (V3D_X(v) < s &&
			V3D_Y(v) < 0.0E+00 &&
			V3D_Z(v) < 0.0E+00);
}

bool
gms::math::operator<(const double s,
			const Vector3D &v) {
	return (s < V3D_X(v)		&&
			0.0E+00 < V3D_Y(v)  &&
			0.0E+00 < V3D_Z(v));
}

bool
gms::math::operator>=(const Vector3D &v1,
		      const Vector3D &v2) {
	return (V3D_X(v1) >= V3D_X(v2) &&
			V3D_Y(v1) >= V3D_Y(v2) &&
			V3D_Z(v1) >= V3D_Z(v2));
}

bool
gms::math::operator>=(const Vector3D &v,
		 const double s) {
	return (V3D_X(v) >= s &&
		    V3D_Y(v) >= 0.0E+00 &&
		    V3D_Z(v) >= 0.0E+00);
}

bool
gms::math::operator>=(const double s,
		      const Vector3D &v) {
	return (s >= V3D_X(v) &&
			0.0E+00 >= V3D_Y(v) &&
			0.0E+00 >= V3D_Z(v));
}

bool
gms::math::operator<=(const Vector3D &v1,
		      const Vector3D &v2) {
	return (V3D_X(v1) <= V3D_X(v2) &&
			V3D_Y(v1) <= V3D_Y(v2) &&
			V3D_Z(v1) <= V3D_Z(v2));
}

bool
gms::math::operator<=(const Vector3D &v,
		      const double s) {
	return (V3D_X(v) <= s &&
			V3D_Y(v) <= 0.0E+00 &&
			V3D_Z(v) <= 0.0E+00);
}

bool
gms::math::operator<=(const double s,
	              const Vector3D &v) {
	return (s <= V3D_X(v) &&
		0.0E+00 <= V3D_Y(v) &&
		0.0E+00 <= V3D_Z(v));
}

double
gms::math::dot(_In_ const Vector3D &v1,
			   _In_ const Vector3D &v2) {
	return (V3D_X(v1)*V3D_X(v2) + 
		    V3D_Y(v1)*V3D_Y(v2) + 
		    V3D_Z(v1)*V3D_Z(v2));
}

double
gms::math::abs_dot(const Vector3D &v1,
		   const Vector3D &v2) {
	return (std::fabs(dot(v1,v2)));
}

gms::math::Vector3D
gms::math::cross(const Vector3D &v1,
		const Vector3D &v2) {
	const double s0 = V3D_Y(v1)*V3D_Z(v2) - 
			          V3D_Z(v1)*V3D_Y(v2);
	const double s1 = V3D_Z(v1)*V3D_X(v2) - 
					  V3D_X(v1)*V3D_Z(v2);
	const double s2 = V3D_X(v1)*V3D_Y(v2) - 
					  V3D_Y(v1)*V3D_X(v2);
	return (Vector3D{s0,s1,s2});
}

double
gms::math::tri_prod(const Vector3D &v1,
			const Vector3D &v2,
			 const Vector3D &v3) {
	const double s0 = V3D_X(v1)*(V3D_Y(v2)*V3D_Z(v3) - 
					             V3D_Z(v2)*V3D_Y(v3));
	const double s1 = V3D_Y(v1)*(V3D_Z(v2)*V3D_X(v3) -
								 V3D_X(v1)*V3D_Z(v3));
	const double s2 = V3D_Z(v2)*(V3D_X(v2)*V3D_Y(v3) - 
								 V3D_Y(v2)*V3D_X(v3));
	return (s0+s1+s2);
}

gms::math::Vector3D
gms::math::dir_cos(const Vector3D &v) {
	double n = norm(v);
	return (Vector3D{ V3D_X(v)/n,
					  V3D_Y(v) / n,
					  V3D_Z(v)/n });
}

double
gms::math::norm(const Vector3D &v) {
	return (std::sqrt(V3D_X(v)*V3D_X(v) + 
					  V3D_Y(v)*V3D_Y(v) + 
					  V3D_Z(v)*V3D_Z(v)));
}

gms::math::Vector3D
gms::math::normalize(const Vector3D &v) {
	double n = norm(v);
	return (Vector3D{ V3D_X(v) / n,
		V3D_Y(v) / n,
		V3D_Z(v) / n });
}
