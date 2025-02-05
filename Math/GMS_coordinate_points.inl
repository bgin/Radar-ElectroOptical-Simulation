
gms::math::CoordPoints<__m128>::CoordPoints() {
	/*Initialize coordinates to zero*/
	this->m_F32Coord3D = _mm_setzero_ps();
}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ const float x, _In_ const float y, _In_ const float z) :
m_F32Coord3D(_mm_set_ps(NAN_FLT, z, y, x)) {

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ const float s) :
m_F32Coord3D(_mm_set_ps(NAN_FLT, s, s, s)) {

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ const double x, _In_ const double y, _In_ const double z) :
/*  Warning: Precion is lost from ~16 digits to ~8*/
m_F32Coord3D(_mm256_cvtpd_ps(_mm256_set_pd(NAN_DBL, z, y, x))) {

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ float(&ar)[4]) :
m_F32Coord3D(_mm_loadu_ps(&ar[0])) {

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ const std::initializer_list<float> &list) :
m_F32Coord3D(_mm_loadu_ps(&list.begin()[0])){

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ const __m128 &v) :
m_F32Coord3D(v) {

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ const CoordPoints &rhs) :
m_F32Coord3D(rhs.m_F32Coord3D){

}

gms::math::CoordPoints<__m128>::CoordPoints(_In_ CoordPoints &&rhs) :
m_F32Coord3D(std::move(rhs.m_F32Coord3D)) {

}


auto    gms::math::CoordPoints<__m128>::operator=(_In_ const CoordPoints &rhs)->gms::math::CoordPoints<__m128> &{

	if (this == &rhs) return (*this);
	this->m_F32Coord3D = rhs.m_F32Coord3D;
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator=(_In_ const CoordPoints &&rhs)->gms::math::CoordPoints<__m128> & {

	if (this == &rhs) return (*this);
	this->m_F32Coord3D = std::move(rhs.m_F32Coord3D);
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator+=(_In_ const CoordPoints &rhs)->gms::math::CoordPoints<__m128> & {

	this->m_F32Coord3D = _mm_add_ps(this->m_F32Coord3D, rhs.m_F32Coord3D);
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator+=(_In_ const float s)->gms::math::CoordPoints<__m128> & {

	this->m_F32Coord3D = _mm_add_ps(this->m_F32Coord3D,_mm_set_ps(NAN_FLT, s, s, s));
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator-=(_In_ const CoordPoints &rhs)->gms::math::CoordPoints<__m128> & {

	this->m_F32Coord3D = _mm_sub_ps(this->m_F32Coord3D, rhs.m_F32Coord3D);
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator-=(_In_ const float s)->gms::math::CoordPoints<__m128> & {

	this->m_F32Coord3D = _mm_sub_ps(this->m_F32Coord3D, _mm_set_ps(NAN_FLT, s, s, s));
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator*=(_In_ const CoordPoints &rhs)->gms::math::CoordPoints<__m128> & {

	this->m_F32Coord3D = _mm_mul_ps(this->m_F32Coord3D, rhs.m_F32Coord3D);
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator*=(_In_ const float s)->gms::math::CoordPoints<__m128> & {

	this->m_F32Coord3D = _mm_mul_ps(this->m_F32Coord3D, _mm_set_ps(NAN_FLT, s, s, s));
	return (*this);
}

auto    gms::math::CoordPoints<__m128>::operator/=(_In_ const CoordPoints &rhs)->gms::math::CoordPoints<__m128> & {

	if (!(_mm_testz_ps(rhs.m_F32Coord3D, _mm_set_ps(NAN_FLT, 0.f, 0.f, 0.f)))){

		this->m_F32Coord3D = _mm_mul_ps(this->m_F32Coord3D, _mm_rcp_ps(rhs.m_F32Coord3D));
		return (*this);
	}
}

auto    gms::math::CoordPoints<__m128>::operator/=(_In_ const float s)->gms::math::CoordPoints<__m128> & {

	if (s != 0.f){
		float inv_s{ 1.f / s };
		this->m_F32Coord3D = _mm_mul_ps(this->m_F32Coord3D, _mm_set_ps(NAN_FLT, inv_s, inv_s, inv_s));
		return (*this);
	}
}

auto    gms::math::CoordPoints<__m128>::operator==(_In_ const CoordPoints &rhs)const-> __m128 {

	return (_mm_cmpeq_ps(this->m_F32Coord3D, rhs.m_F32Coord3D));
}

auto    gms::math::CoordPoints<__m128>::operator==(_In_ const float s)const-> __m128 {

	// 3rd scalar counting from 0 = Don't care
	return (_mm_cmpeq_ps(this->m_F32Coord3D, _mm_set_ps(NAN_FLT, s, s, s)));
}

auto    gms::math::CoordPoints<__m128>::operator!=(_In_ const CoordPoints &rhs)const-> __m128 {

	return (_mm_cmpneq_ps(this->m_F32Coord3D, rhs.m_F32Coord3D));
}

auto    gms::math::CoordPoints<__m128>::operator!=(_In_ const float s)const-> __m128 {
	// 3rd scalar counting from 0 = Don't care
	return (_mm_cmpneq_ps(this->m_F32Coord3D, _mm_set_ps(NAN_FLT, s, s, s)));
}

gms::math::CoordPoints<__m128>::operator __m128 () {

	return this->m_F32Coord3D;
}

gms::math::CoordPoints<__m128>::operator  __m128 const () const {

	return this->m_F32Coord3D;
}

auto    gms::math::operator<<(_In_ std::ostream &os, _In_ const gms::math::CoordPoints<__m128> &rhs)->std::ostream & {

	os << "X: " << reinterpret_cast<const double*>(&rhs.m_F32Coord3D)[0] << std::endl;
	os << "Y: " << reinterpret_cast<const double*>(&rhs.m_F32Coord3D)[1] << std::endl;
	os << "Z: " << reinterpret_cast<const double*>(&rhs.m_F32Coord3D)[2] << std::endl;
	return os;
}

inline  auto  gms::math::CoordPoints<__m128>::getF32Coord3D()const-> __m128 {
	return this->m_F32Coord3D;
}

inline  auto  gms::math::CoordPoints<__m128>::X()const->float {
	return (reinterpret_cast<const double*>(&this->m_F32Coord3D)[0]);
}

inline  auto  gms::math::CoordPoints<__m128>::Y()const->float {
	return (reinterpret_cast<const double*>(&this->m_F32Coord3D)[1]);
}

inline  auto  gms::math::CoordPoints<__m128>::Z()const->float {
	return (reinterpret_cast<const double*>(&this->m_F32Coord3D)[2]);
}

inline  auto  gms::math::CoordPoints<__m128>::CArray()->double * {

	double* pF32Coord3D = nullptr;
	pF32Coord3D = reinterpret_cast<double*>(&this->m_F32Coord3D);
	return (pF32Coord3D);
}

inline  auto  gms::math::CoordPoints<__m128>::CArray()const->const double * {

	const double* pF32Coord3D = nullptr;
	pF32Coord3D = reinterpret_cast<const double*>(&this->m_F32Coord3D);
	return (pF32Coord3D);
}

inline  auto  gms::math::CoordPoints<__m128>::setX(_In_ const float s)->void {

	this->m_F32Coord3D.m128_f32[0] = s;
}

inline  auto  gms::math::CoordPoints<__m128>::setY(_In_ const float s)->void {

	this->m_F32Coord3D.m128_f32[1] = s;
}

inline  auto  gms::math::CoordPoints<__m128>::setZ(_In_ const float s)->void {

	this->m_F32Coord3D.m128_f32[2] = s;
}

inline  auto  gms::math::CoordPoints<__m128>::magnitude()const->float {
	auto temp(_mm_mul_ps(this->m_F32Coord3D, this->m_F32Coord3D));
	return (reinterpret_cast<const double*>(&temp)[0] +
		reinterpret_cast<const double*>(&temp)[1] +
		reinterpret_cast<const double*>(&temp)[2]);
}

inline  auto  gms::math::CoordPoints<__m128>::perpendicular()const->float {
	auto temp(_mm_mul_ps(this->m_F32Coord3D, this->m_F32Coord3D));
	return (reinterpret_cast<const double*>(&temp)[0] + reinterpret_cast<const double*>(&temp)[1]);
}

inline  auto  gms::math::CoordPoints<__m128>::rho()const->float {

	return (std::sqrt(this->perpendicular()));
}

inline  auto  gms::math::CoordPoints<__m128>::phi()const->float {

	return((this->X() == 0.f && this->Y() == 0.f) ? 0.f : std::atan2(this->X(), this->Y())); 
}

inline  auto  gms::math::CoordPoints<__m128>::theta()const->float {

	return ((this->X() == 0.f && this->Y() == 0.f && this->Z()) ? 0.f : std::atan2(this->rho(), this->Z()));
}
