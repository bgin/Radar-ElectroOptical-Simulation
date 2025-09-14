
gms::math::
AVXVecF64::AVXVecF64()
:
m_vf64(_mm256_setzero_pd()) {}

gms::math::
AVXVecF64::AVXVecF64(const double v[4])
:
m_vf64(_mm256_loadu_pd(v)) {}

gms::math::
AVXVecF64::AVXVecF64(const double s1,
		    const double s2,
		    const double s3,
		    const double s4)
:
m_vf64(_mm256_setr_pd(s1, s2, s3, s4)){}

gms::math::
AVXVecF64::AVXVecF64(const double s)
:
m_vf64(_mm256_set1_pd(s)) {}

gms::math::
AVXVecF64::AVXVecF64(const __m256d v)
:
m_vf64(v) {}

gms::math::
AVXVecF64::AVXVecF64(const __m256i v)
:
m_vf64(_mm256_castsi256_pd(v)) {}

gms::math::
AVXVecF64::AVXVecF64(const AVXVecF64 &x)
:
m_vf64(x.m_vf64) {}

gms::math::
AVXVecF64::AVXVecF64(const __m128d v1,
			const __m128d v2){
	m_vf64 = _mm256_insertf128_pd(m_vf64,v1,0);
	m_vf64 = _mm256_insertf128_pd(m_vf64,v2,1);
}

const __m256d 
gms::math::AVXVecF64::get_vf64() const {
	return (m_vf64);
}

__m256d
gms::math::AVXVecF64::get_vf64() {
	return (m_vf64);
}

__m128d
gms::math::AVXVecF64::lo_part() const {
	return (_mm256_extractf128_pd(m_vf64,0));
}

__m128d
gms::math::AVXVecF64::hi_part() const {
	return (_mm256_extractf128_pd(m_vf64,1));
}

gms::math::AVXVecF64 &
gms::math::AVXVecF64::load_a(const double* __restrict adress) {
	m_vf64 = _mm256_load_pd(adress);
	return (*this);
}

gms::math::AVXVecF64 &
gms::math::AVXVecF64::load_u(const double* __restrict address) {
	m_vf64 = _mm256_loadu_pd(address);
	return (*this);
}

void
gms::math::AVXVecF64::store_a(double* __restrict address) const {
	_mm256_store_pd(address,m_vf64);
}

void
gms::math::AVXVecF64::store_u(double* __restrict address) const {
	_mm256_storeu_pd(address,m_vf64);
}

void
gms::math::AVXVecF64::stream_store(double* __restrict address) const {
	_mm256_stream_pd(address,m_vf64);
}

double
gms::math::AVXVecF64::extract_scalar(const uint32_t idx) const {
	__declspec(align(32))double t[4] = {0.0};
	store_a(t);
	return (t[idx & 3]);
}

gms::math::AVXVecF64 const &
gms::math::AVXVecF64::insert(const uint32_t idx,
			const double value) {
	__m256d ymm = _mm256_broadcast_sd(&value);
	switch (idx) {
		case 0: 
			m_vf64 = _mm256_blend_pd(m_vf64,ymm,1);
			break;
	    case 1:
			m_vf64 = _mm256_blend_pd(m_vf64,ymm,2);
			break;
	    case 2:
			m_vf64 = _mm256_blend_pd(m_vf64,ymm,4);
			break;
	    default:
			m_vf64 = _mm256_blend_pd(m_vf64,ymm,8);
		    break;
	}
	return (*this);
}

gms::math::AVXVecF64 &
gms::math::AVXVecF64::operator=(const AVXVecF64 &x) {
	if (this == &x)
		return (*this);
	m_vf64 = x.m_vf64;
	return (*this);
}

gms::math::AVXVecF64::operator __m256d () const {
	return (m_vf64);
}

double
gms::math::AVXVecF64::operator[](const uint32_t idx) const {
	return (reinterpret_cast<const double*>(&m_vf64)[idx]);
}

std::ostream &
gms::math::operator<<(std::ostream &os, 
			 const AVXVecF64 &x) {
	os << std::fixed << std::showpoint << "x=" << x.m_vf64.m256d_f64[0] <<
										  "y=" << x.m_vf64.m256d_f64[1] <<
										  "z=" << x.m_vf64.m256d_f64[2] <<
										  "w=" << x.m_vf64.m256d_f64[3] << "\n";
	return (os);
}

__m128d
gms::math::extract(AVXVecF64 &x, const int32_t idx) {
	return (_mm256_extractf128_pd(x,idx));
}

gms::math::AVXVecF64
gms::math::select_vec(const AVXVecF64 &a,
		      const AVXVecF64 &b,
		      const __m256d pred) {

	return (_mm256_blendv_pd(a,b,pred));
}

gms::math::AVXVecF64
gms::math::max(const AVXVecF64 &x,
	       const AVXVecF64 &y) {
	return (_mm256_max_pd(x,y));
}

gms::math::AVXVecF64
gms::math::min(const AVXVecF64 &x,
		const AVXVecF64 &y) {
	return (_mm256_min_pd(x,y));
}

gms::math::AVXVecF64
gms::math::abs(const AVXVecF64 &x) {
	const __m256d mask = _mm256_set1_pd(0x7FFFFFFFFFFFFFF);
	return (_mm256_and_pd(x,mask));
}

gms::math::AVXVecF64
gms::math::sqrt(const AVXVecF64 &x) {
	return (_mm256_sqrt_pd(x));
}

gms::math::AVXVecF64
gms::math::sqr(const AVXVecF64 &x) {
	return (x*x);
}

gms::math::AVXVecF64
gms::math::ceil(const AVXVecF64 &x) {
	return (_mm256_ceil_pd(x));
}

gms::math::AVXVecF64
gms::math::floor(const AVXVecF64 &x) {
	return (_mm256_floor_pd(x));
}

gms::math::AVXVecF64
gms::math::round(const AVXVecF64 &x,
		const int32_t by) {
	return (_mm256_round_pd(x,by));
}

gms::math::AVXVecF64
gms::math::sin(const AVXVecF64 &x) {
	return (_mm256_sin_pd(x));
}

gms::math::AVXVecF64
gms::math::cos(const AVXVecF64 &x) {
	return (_mm256_cos_pd(x));
}

gms::math::AVXVecF64
gms::math::sinh(const AVXVecF64 &x) {
	return (_mm256_sinh_pd(x));
}

gms::math::AVXVecF64
gms::math::cosh(const AVXVecF64 &x) {
	return (_mm256_cosh_pd(x));
}

gms::math::AVXVecF64
gms::math::tan(const AVXVecF64 &x) {
	return (_mm256_tan_pd(x));
}

gms::math::AVXVecF64
gms::math::tanh(const AVXVecF64 &x) {
	return (_mm256_tanh_pd(x));
}

gms::math::AVXVecF64
gms::math::asin(const AVXVecF64 &x) {
	return (_mm256_asin_pd(x));
}

gms::math::AVXVecF64
gms::math::asinh(const AVXVecF64 &x) {
	return (_mm256_asinh_pd(x));
}

lam::math::AVXVecF64
lam::math::acos(_In_ const AVXVecF64 &x) {
	return (_mm256_acos_pd(x));
}

gms::math::AVXVecF64
gms::math::acosh(const AVXVecF64 &x) {
	return (_mm256_acosh_pd(x));
}

gms::math::AVXVecF64
gms::math::atan(const AVXVecF64 &x) {
	return (_mm256_atan_pd(x));
}

gms::math::AVXVecF64
gms::math::atanh(const AVXVecF64 &x) {
	return (_mm256_atanh_pd(x));
}

gms::math::AVXVecF64
gms::math::atan2(const AVXVecF64 &x,
			   const AVXVecF64 &y) {
	return (_mm256_atan2_pd(x,y));
}

gms::math::AVXVecF64
gms::math::unary_minus(const AVXVecF64 &x) {
	return (_mm256_sub_pd(_mm256_setzero_pd(),x));
}

gms::math::AVXVecF64
gms::math::exp(AVXVecF64 &x) {
	return (_mm256_exp_pd(x));
}

gms::math::AVXVecF64
gms::math::log10(const AVXVecF64 &x) {
	return (_mm256_log10_pd(x));
}

gms::math::AVXVecF64
gms::math::log(const AVXVecF64 &x) {
	return (_mm256_log_pd(x));
}

gms::math::AVXVecF64
gms::math::pow(const AVXVecF64 &x,
		const AVXVecF64 &y) {
	return (_mm256_pow_pd(x,y));
}

gms::math::AVXVecF64
gms::math::fmadd(const AVXVecF64 &x,
		const AVXVecF64 &y,
		const AVXVecF64 &z) {
	return (_mm256_fmadd_pd(x,y,z));
}

gms::math::AVXVecF64
gms::math::fmadsubb(const AVXVecF64 &x,
	        const AVXVecF64 &y,
		const AVXVecF64 &z) {
	return (_mm256_fmaddsub_pd(x,y,z));
}

gms::math::AVXVecF64
gms::math::fmsub(const AVXVecF64 &x,
		 const AVXVecF64 &y,
		 const AVXVecF64 &z) {
	return (_mm256_fmsub_pd(x,y,z));
}

gms::math::AVXVecF64
gms::math::fmsubadd(const AVXVecF64 &x,
			const AVXVecF64 &y,
			const AVXVecF64 &z) {
	return (_mm256_fmsubadd_pd(x,y,z));
}

gms::math::AVXVecF64
gms::math::fnmadd(const AVXVecF64 &x,
			const AVXVecF64 &y,
			const AVXVecF64 &z) {
	return (_mm256_fnmadd_pd(x,y,z));
}

gms::math::AVXVecF64
gms::math::fnmsub(const AVXVecF64 &x,
		const AVXVecF64 &y,
		const AVXVecF64 &z) {
	return (_mm256_fnmsub_pd(x,y,z));
}

gms::math::AVXVecF64
gms::math::operator+(const AVXVecF64 &x,
			const AVXVecF64 &y) {
	return (_mm256_add_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator+(const AVXVecF64 &x,
			const double s) {
	return (x + AVXVecF64(s));
}

gms::math::AVXVecF64
gms::math::operator+(const double s,
			const AVXVecF64 &x) {
	return (AVXVecF64(s) + x);
}

gms::math::AVXVecF64
gms::math::operator+=(AVXVecF64 &x,
		 const AVXVecF64 &y) {
	x = x + y;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator++(AVXVecF64 &x) {
	x = x + 1.0L;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator-(const AVXVecF64 &x,
		    const AVXVecF64 &y) {
	return (_mm256_sub_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator-(const AVXVecF64 &x,
		const double s) {
	return (x - AVXVecF64(s));
}

gms::math::AVXVecF64
gms::math::operator-(const double s,
			const AVXVecF64 &x) {
	return (AVXVecF64(s) - x);
}

gms::math::AVXVecF64
gms::math::operator-=(AVXVecF64 &x,
			const AVXVecF64 &y) {
	x = x - y;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator--(AVXVecF64 &x) {
	x = x - 1.0L;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator*(const AVXVecF64 &x,
			const AVXVecF64 &y) {
	return (_mm256_mul_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator*(const AVXVecF64 &x,
			const double s) {
	return (x * AVXVecF64(s));
}

gms::math::AVXVecF64
gms::math::operator*(const double s,
			 const AVXVecF64 &x) {
	return (AVXVecF64(s) * x);
}

gms::math::AVXVecF64
gms::math::operator*=(AVXVecF64 &x,
			const AVXVecF64 &y) {
	x = x * y;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator/(const AVXVecF64 &x,
		 const AVXVecF64 &y) {
	return (_mm256_div_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator/(const AVXVecF64 &x,
		     const double s) {
	return (x / AVXVecF64(s));
}

gms::math::AVXVecF64
gms::math::operator/(const double s,
			 const AVXVecF64 &x) {
	return (AVXVecF64(s) / x);
}

gms::math::AVXVecF64
gms::math::operator/=(AVXVecF64 &x,
		      const AVXVecF64 &y) {
	x = x / y;
	return (x);
}

__m256d
gms::math::operator==(const AVXVecF64 &x,
		      const AVXVecF64 &y) {
	return (_mm256_cmp_pd(x,y,_CMP_EQ_OQ));
}

__m256d
gms::math::operator==(const AVXVecF64 &x,
		     const double y) {
	return (_mm256_cmp_pd(x, AVXVecF64{y},_CMP_EQ_OQ));
}

__m256d
gms::math::operator==(const double x,
                      const AVXVecF64 &y) {
	return (_mm256_cmp_pd(AVXVecF64{x},y,_CMP_EQ_OQ));
}

__m256d
gms::math::operator!=(const AVXVecF64 &x,
		      const AVXVecF64 &y) {
	return (_mm256_cmp_pd(x,y,_CMP_NEQ_OQ));
}

__m256d
gms::math::operator>(const AVXVecF64 &x,
			const AVXVecF64 &y) {
	return (_mm256_cmp_pd(x,y,_CMP_GT_OQ));
}

__m256d
gms::math::operator<(const AVXVecF64 &x,
			const AVXVecF64 &y) {
	return (_mm256_cmp_pd(x,y,_CMP_LT_OQ));
}

__m256d
lam::math::operator>=(_In_ const AVXVecF64 &x,
					  _In_ const AVXVecF64 &y) {
	return (_mm256_cmp_pd(x, y, _CMP_GE_OQ));
}

__m256d
gms::math::operator<=(const AVXVecF64 &x,
			const AVXVecF64 &y) {
	return (_mm256_cmp_pd(x, y, _CMP_LE_OQ));
}

gms::math::AVXVecF64
gms::math::operator&(const AVXVecF64 &x,
		      const AVXVecF64 &y) {
	return (_mm256_and_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator&=(AVXVecF64 &x,
				 const AVXVecF64 &y) {
	x = x & y;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator|(const AVXVecF64 &x,
			const AVXVecF64 &y) {
	return (_mm256_or_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator|=(AVXVecF64 &x,
					 const AVXVecF64 &y) {
	x = x | y;
	return (x);
}

gms::math::AVXVecF64
gms::math::operator^(const AVXVecF64 &x,
		const AVXVecF64 &y) {
	return (_mm256_xor_pd(x,y));
}

gms::math::AVXVecF64
gms::math::operator^=(AVXVecF64 &x,
					  const AVXVecF64 &y) {
	x = x ^ y;
	return (x);
}
