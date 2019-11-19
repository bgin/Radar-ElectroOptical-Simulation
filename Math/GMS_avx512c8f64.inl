
const gms::math::AVX512c8f64
gms::math::AVX512c8f64::CZERO = AVX512c8f64{};

gms::math::AVX512c8f64
::AVX512c8f64() {
	m_re = _mm512_setzero_pd();
	m_im = _mm512_setzero_pd();
}

gms::math::AVX512c8f64
::AVX512c8f64(const AVX512c8Payload x) {
	m_re = _mm512_set_pd(x.re0, x.re1, x.re2, x.re3, x.re4, x.re5, x.re6, x.re7);
	m_im = _mm512_set_pd(x.im0, x.im1, x.im2, x.im3, x.im4, x.im5, x.im6, x.im7);
}

gms::math::AVX512c8f64
::AVX512c8f64(const double * __restrict Re,
              const double * __restrict Im) {
	m_re = _mm512_load_pd(&Re[0]);
	m_im = _mm512_load_pd(&Im[0]);
}

gms::math::AVX512c8f64
::AVX512c8f64(const double re,
	      const double im) {
	m_re = _mm512_set1_pd(re);
	m_im = _mm512_set1_pd(im);
}

gms::math::AVX512c8f64
::AVX512c8f64(const double c) {
	m_re = _mm512_set1_pd(c);
	m_im = _mm512_setzero_pd();
}

gms::math::AVX512c8f64
::AVX512c8f64(const double re0,
	     const double re1,
	     const double re2,
             const double re3,
	     const double re4,
	     const double re5,
	     const double re6,
	     const double re7) {
	m_re = _mm512_set_pd(re0, re1, re2, re3, re4, re5, re6, re7);
	m_im = _mm512_setzero_pd();
}

gms::math::AVX512c8f64
::AVX512c8f64(const __m512d re,
	      const __m512d im) {
	m_re = re;
	m_im = im;
}

gms::math::AVX512c8f64
::AVX512c8f64(const AVX512c8f64 &x) {
	m_re = x.m_re;
	m_im = x.m_im;
}

gms::math::AVX512c8f64 &
gms::math::AVX512c8f64::load_a(const double * __restrict Re,
			       const double * __restrict Im) {
	m_re = _mm512_load_pd(&Re[0]);
	m_im = _mm512_load_pd(&Im[0]);
	return (*this);
}

gms::math::AVX512c8f64 &
gms::math::AVX512c8f64::load_u(const double * __restrict Re,
			       const double * __restrict Im) {
	m_re = _mm512_loadu_pd(&Re[0]);
	m_im = _mm512_loadu_pd(&Im[0]);
	return (*this);
}

void 
gms::math::AVX512c8f64::store_a(double * __restrict Re,
				double * __restrict Im) const {
	_mm512_store_pd(&Re[0], m_re);
	_mm512_store_pd(&Im[0], m_im);
}

void 
gms::math::AVX512c8f64::store_u(double * __restrict Re,
			      double * __restrict Im) const {
	_mm512_storeu_pd(&Re[0], m_re);
	_mm512_storeu_pd(&Im[0], m_im);
}

void 
gms::math::AVX512c8f64::stream(double * __restrict Re,
			       double * __restrict Im) const {
	_mm512_stream_pd(&Re[0], m_re);
	_mm512_stream_pd(&Im[0], m_im);
	_mm_sfence();
}

double 
gms::math::AVX512c8f64::extract_f64(const int32_t idx) const {
#if defined _WIN64
	__declspec(align(64)) double mem[16] = {};
#elif defined __linux
        __attribute__((align(64))) double mem[16] = {};
#endif
	store_a(&mem[0], &mem[8]);
	return (mem[idx & 0xF]);
}

std::pair<double, double>
gms::math::AVX512c8f64::extract_2f64(const int32_t re_idx,
				     const int32_t im_idx) {
#if defined _WIN64
	__declspec(align(64)) double re_mem[8] = {};
	__declspec(align(64)) double im_mem[8] = {};
#elif defined __linux
        __attribute__((align(64))) double re_mem[8] = {};
	__attribute__((align(64))) double im_mem[8] = {};
#endif
	store_a(&re_mem[0], &im_mem[0]);
	return (std::make_pair(re_mem[re_idx & 0x7], im_mem[im_idx & 0x7]));
}

gms::math::AVX512c8f64 & 
gms::math::AVX512c8f64::insert_1f64(const int32_t idx,
				    const double val) {
#if defined _WIN64
	__declspec(align(64)) double mem[16] = {};
#elif defined __linux
        __attribute__((align(64))) double mem[16] = {};
#endif
	store_a(&mem[0], &mem[8]);
	mem[idx & 0xF] = val;
	load_a(&mem[0], &mem[8]);
	return (*this);
}

gms::math::AVX512c8f64 & 
gms::math::AVX512c8f64::insert_2f64(const int32_t re_idx,
				    const int32_t im_idx,
				    const double re,
				    const double im) {
#if defined _WIN64				    
	__declspec(align(64)) double mem_re[8] = {};
	__declspec(align(64)) double mem_im[8] = {};
#elif defined __linux
        __attribute__((align(64))) double mem_re[8] = {};
	__attribute__((align(64))) double mem_im[8] = {};
#endif
	store_a(&mem_re[0], &mem_im[0]);
	mem_re[re_idx & 0x7] = re;
	mem_im[im_idx & 0x7] = im;
	load_a(&mem_re[0], &mem_im[0]);
	return (*this);
}

// Length of 16 doubles
void 
gms::math::AVX512c8f64
::concatenate_a(double * __restrict out) const {
	store_a(&out[0], &out[8]);
}

// Length of 16 doubles
void 
gms::math::AVX512C8f64
::concatenate_u(double * __restrict out) const {
	store_u(&out[0], &out[8]);
}

gms::math::AVX512c8f64 &
gms::math::AVX512c8f64
::partial_loadu(const double * __restrict Re,
		const int32_t n_re,
		const double * __restrict Im,
		const int32_t n_im) {
	m_re = _mm512_maskz_loadu_pd(__mmask8((1 << n_re)-1),Re);
	m_im = _mm512_maskz_loadu_pd(__mmask8((1 << n_im)-1),Im);
	return (*this);
}

gms::math::AVX512c8f64 &
gms::math::AVX512C8f64
::partial_loada(const double * __restrict Re,
		const int32_t n_re,
		const double * __restrict Im,
		const int32_t n_im) {
	m_re = _mm512_maskz_load_pd(__mmask8((1 << n_re) - 1), Re);
	m_im = _mm512_maskz_load_pd(__mmask8((1 << n_im) - 1), Im);
	return (*this);
}

void
gms::math::AVX512c8f64
::partial_storeu(double * __restrict Re,
		 const int32_t n_re,
		 double * __restrict Im,
		 const int32_t n_im) {
	_mm512_mask_storeu_pd(&Re[0],__mmask8((1<<n_re)-1),m_re);
	_mm512_mask_storeu_pd(&Im[0],__mmask8((1<<n_im)-1),m_im);
}

void
gms::math::AVX512c8f64
::partial_storea(double * __restrict Re,
		 const int32_t n_re,
		 double * __restrict Im,
		 const int32_t n_im) {
	_mm512_mask_store_pd(&Re[0], __mmask8((1 << n_re) - 1), m_re);
	_mm512_mask_store_pd(&Im[0], __mmask8((1 << n_im) - 1), m_im);
}

gms::math::AVX512c8f64 & 
gms::math::AVX512c8f64
::expand(const AVX512c8f64 x,
         const __mmask8 mask) {
	m_re = _mm512_maskz_expand_pd(mask, x.m_re);
	m_im = _mm512_maskz_expand_pd(mask, x.m_im);
	return (*this);
}

gms::math::AVX512c8f64 & 
gms::math::AVX512c8f64
::expand_load(const AVX512c8f64 x,
	      const __mmask8 mask,
	      const double * __restrict re,
	      const double * __restrict im) {
	m_re = _mm512_mask_expandloadu_pd(x.m_re, mask, &re[0]);
	m_im = _mm512_mask_expandloadu_pd(x.m_im, mask, &im[0]);
	return (*this);
}

gms::math::AVX512c8f64 & 
gms::math::AVX512c8f64
::permute(const __mmask8 mask,
	  const int32_t imm) {
	m_re = _mm512_mask_permute_pd(m_re, mask, m_im, imm);
	m_im = _mm512_mask_permute_pd(m_im, mask, m_re, imm);
	return (*this);
}

__m256d 
gms::math::AVX512c8f64::re_low2() const {
	return (_mm512_extractf64x4_pd(m_re, 0));
}



__m256d 
gms::math::AVX512c8f64::re_hi2() const {
	return (_mm512_extractf64x4_pd(m_re, 1));
}

__m256d 
gms::math::AVX512c8f64::im_low2() const {
	return (_mm512_extractf64x4_pd(m_im, 0));
}

__m256d 
gms::math::AVX512c8f64::im_hi2() const {
	return (_mm512_extractf64x4_pd(m_im, 1));
}

gms::math::AVX512c8f64 & 
gms::math::AVX512c8f64
::operator=(const AVX512c8f64 x) {
	if (this == &x) { return (*this); }
	m_re = x.m_re;
	m_im = x.m_im;
	return (*this);
}

static inline gms::math::AVX512c8f64
gms::math::conj(const AVX512c8f64 x) {
     auto tmp = ~x
     const __m512d im_part = tmp.m_im;
     return (AVX512c8f64{x.m_re,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::polar(const __m512d rho,
                 const __m512d theta) {
       const __m512d re_part =
           _mm512_mul_pd(rho,_mm512_cos_pd(theta));
       const __m512d im_part =
           _mm512_mul_pd(rho,_mm512_sin_pd(theta));
       return (AVX512c8f64{re_part,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::carg(const AVX512c8f64 x) {
       const __m512d re_part =
             _mm512_atan2_pd(x.m_im,x.m_re);
       return (AVX512c8f64{re_part,_mm512_setzero_pd()});
}

static inline gms::math::AVX512c8f64
gms::math::carg(const double re,
                const double im) {
       const __m512d real = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,re);
       const __m512d imag = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,im);
       const __m512d re_part =
               _mm512_atan2_pd(imag,real);
       return (AVX512c8f64{re_part,_mm512_setzero_pd()});
}

static inline gms::math::AVX512c8f64
gms::math::csin(const AVX512c8f64 x) {
	const __m512d re_part =
		_mm512_mul_pd(_mm512_sin_pd(x.m_re),
		_mm512_cosh_pd(x.m_im));
	const __m512d im_part =
		_mm512_mul_pd(_mm512_cos_pd(x.m_re),
		_mm512_sinh_pd(x.m_im));
	return (AVX512c8f64{ re_part, im_part });
}

static inline gms::math::AVX512c8f64
gms::math::csin(const double re,
                const double im) {
       const __m512d real = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,re);
       const __m512d imag = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,im);
       const __m512d re_part =
               _mm512_mul_pd(_mm512_sin_pd(real),
	                     _mm512_cosh_pd(imag));
       const __m512d im_part =
               _mm512_mul_pd(_mm512_cos_pd(real),
	                     _mm512_sinh_pd(imag));
       return (AVX512c8f64{re_part,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::csinh(const AVX512c8f64 x) {
        const __m512d re_part =
	       _mm512_mul_pd(_mm512_sinh_pd(x.m_re),
	       _mm512_cos_pd(x.m_im));
	const __m512d im_part =
	       _mm512_mul_pd(_mm512_cosh_pd(x.m_re),
	       _mm512_sin_pd(x.m_im));
	return (AVX512c8f64{re_part,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::csinh(const double re,
                 const double im) {
       const __m512d real = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,re);
       const __m512d imag = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,im);
       const __m512d re_part =
	       _mm512_mul_pd(_mm512_sinh_pd(real),
	       _mm512_cos_pd(imag));
       const __m512d im_part =
	       _mm512_mul_pd(_mm512_cosh_pd(real),
	       _mm512_sin_pd(imag));
	return (AVX512c8f64{re_part,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::ccos(const AVX512c8f64 x) {
	const __m512d re_part =
		_mm512_mul_pd(_mm512_cos_pd(x.m_re),
		_mm512_cosh_pd(x.m_im));
	const __m512d im_part =
		_mm512_mul_pd(_mm512_sin_pd(x.m_re),
		_mm512_sinh_pd(x.m_im));
	
	return (AVX512c8f64{ re_part, im_part });
}

static inline gms::math::AVX512c8f64
gms::math::ccos(const double re,
                const double im) {
       const __m512d real = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,re);
       const __m512d imag = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,im);
       const __m512d re_part =
		_mm512_mul_pd(_mm512_cos_pd(real),
		_mm512_cosh_pd(imag));
       const __m512d im_part =
		_mm512_mul_pd(_mm512_sin_pd(real),
		_mm512_sinh_pd(imag));
	
       return (AVX512c8f64{ re_part, im_part });
}

static inline gms::math::AVX512c8f64
gms::math::ccosh(const AVX512c8f64 x) {
        const __m512d re_part =
	       _mm512_mul_pd(_mm512_cosh_pd(x.m_re),
	       _mm512_cos_pd(x.m_im));
	const __m512d im_part =
	       _mm512_mul_pd(_mm512_sinh_pd(x.m_re),
	       _mm512_sin_pd(x.m_im));

        return (AVX512c8f64{re_part,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::ccosh(const double re,
                 const double im) {
       const __m512d real = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,re);
       const __m512d imag = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,im);
       const __m512d re_part =
	       _mm512_mul_pd(_mm512_cosh_pd(real),
	       _mm512_cos_pd(imag));
       const __m512d im_part =
	       _mm512_mul_pd(_mm512_sinh_pd(real),
	       _mm512_sin_pd(imag));

       return (AVX512c8f64{re_part,im_part});
}

static inline gms::math::AVX512c8f64
gms::math::cexp(const AVX512c8f64 x) {
	const __m512d re_part =
		_mm512_mul_pd(_mm512_exp_pd(x.m_re),
		_mm512_cos_pd(x.m_im));
	const __m512d im_part =
		_mm512_mul_pd(_mm512_exp_pd(x.m_re),
		_mm512_sin_pd(x.m_im));

	return (AVX512c8f64{ re_part, im_part });
}

static inline __m512d
gms::math::cabs(const AVX512c8f64 x) {
	const __m512d re_part =
		_mm512_mul_pd(x.m_re, x.m_re);
	const __m512d im_part =
		_mm512_mul_pd(x.m_im, x.m_im);
	return ( _mm512_sqrt_pd(_mm512_add_pd(re_part, im_part)));
	
}

static inline __m512d
gms::math::cabs(const double re,
                const double im) {
       const __m512d real = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,re);
       const __m512d imag = _mm512_set_pd(1.0,1.0,1.0,1.0,1.0,1.0,1.0,im);
       const __m512d re_part =
		_mm512_mul_pd(real,real);
       const __m512d im_part =
		_mm512_mul_pd(imag,imag);
       return ( _mm512_sqrt_pd(_mm512_add_pd(re_part, im_part)));
		
}

static inline gms::math::AVX512c8f64
gms::math::cpowi(const AVX512c8f64 x,
                 const double n) {
	const __m512d re_part =
		_mm512_mul_pd(x.m_re, x.m_re);
	const __m512d im_part =
		_mm512_mul_pd(x.m_im, x.m_im);
	const __m512d r =
		_mm512_sqrt_pd(_mm512_add_pd(re_part, im_part));
	const __m512d theta =
		_mm512_atan_pd(_mm512_div_pd(x.m_im, x.m_re));
	const __m512d vn = _mm512_set1_pd(n);
	const __m512d pow_term = _mm512_pow_pd(r, vn);
	const __m512d trig_arg = _mm512_mul_pd(vn, theta);
	return (AVX512c8f64{ _mm512_mul_pd(pow_term, _mm512_cos_pd(trig_arg)),
		_mm512_mul_pd(pow_term, _mm512_sin_pd(trig_arg)) });
}

static inline gms::math::AVX512c8f64
gms::math::clog(const AVX512c8f64 x) {
        auto tmp1 = cabs(x);
	auto tmp2 = carg(x);
        const __m512d re_part =
             _mm512_log_pd(tmp1.m_re);

        return (AVX512c8f64{re_part,tmp2});    
}

static inline gms::math::AVX512c8f64
gms::math::csqrt(const AVX512c8f64 x) {
       auto tmp = cabs(x);
       const __m512d re_part =
             _mm512_mul_pd(_mm512_set1_pd(0.5),
	                     _mm512_add_pd(tmp.m_re,x.m_re));
       const __m512d im_part =
             _mm512_mul_pd(_mm512_set1_pd(0.5),
	                     _mm512_sub_pd(tmp.m_re,x.m_re));
       return (AVX512c8f64{_mm512_sqrt_pd(re_part),
                           _mm512_sqrt_pd(im_part)});
}

static inline gms::math::AVX512c8f64
gms::math::ctan(const AVX512c8f64 x) {
    return (csin(x)/ccos(x));
}

static inline gms::math::AVX512c8f64
gms::math::ctan(const double re,
                const double im) {
    return (csin(re,im)/ccos(re,im));
}

static inline gms::math::AVX512c8f64
gms::math::ctanh(const AVX512c8f64 x) {
    return (csinh(x)/ccosh(x));
}

static inline gms::math::AVX512c8f64
gms::math::ctanh(const double re,
                 const double im) {
    return (csinh(re,im)/ccosh(re,im));
}



static inline gms::math::AVX512c8f64
gms::math::select(const AVX512c8f64 x,
		  const AVX512c8f64 y,
		  const __mmask8 mask) {

	return (AVX512c8f64{ _mm512_mask_blend_pd(mask, x.m_re, y.m_re),
		_mm512_mask_blend_pd(mask, x.m_im, y.m_im) });
}

static inline gms::math::AVX512c8f64
gms::math::cdiv_smith(const AVX512c8f64 x,
                      const AVX512c8f64 y) {

    __m512d ratio,denom,re_part,im_part;
    constexpr __mmask8 all_ones = 0xFF;
    re_part = _mm512_setzero_pd();
    im_part = _mm512_setzero_pd();
    ratio = _mm512_setzero_pd();
    denom = _mm512_setzero_pd();
    __mask8 is_gte = _mm512_cmp_pd_mask(
                           _mm512_abs_pd(y.m_re),
			     _mm512_abs_pd(y.m_im),_CMP_GE_OQ);
   if(is_gte == all_ones) {
      ratio = _mm512_div_pd(y.m_im,y.m_re);
      denom = _mm512_add_pd(y.m_re,_mm512_mul_pd(ratio,y.m_im));
      re_part = _mm512_div_pd(_mm512_add_pd(x.m_re,
                            _mm512_mul_pd(x.m_im,ratio)),denom);
      im_part = _mm512_div_pd(_mm512_sub_pd(x.m_im,
                            _mm512_mul_pd(x.m_re,ratio)),denom);
      return (AVX512c8f64{re_part,im_part});
   }
   else {
      ratio = _mm512_div_pd(y.m_re,y.m_im);
      denom = _mm512_add_pd(y.m_im,_mm512_mul_pd(ratio,y.m_re));
      re_part = _mm512_div_pd(_mm512_add_pd(
                                _mm512_mul_pd(x.m_re,ratio)),denom);
      im_part = _mm512_div_pd(_mm512_sub_pd(
                                _mm512_mul_pd(x.m_im,ratio)),denom);
      return (AVX512c8f64{re_part,im_part});
   }
}

static inline gms::math::AVX512c8f64
gms::math::operator+(const AVX512c8f64 x,
		     const AVX512c8f64 y) {
	return (AVX512c8f64{ _mm512_add_pd(x.m_re, y.m_re),
		_mm512_add_pd(x.m_im, y.m_im) });
}

static inline gms::math::AVX512c8f64
gms::math::operator+(const AVX512c8f64 x,
                     const __m512d v) {
        return (AVX512c8f64{_mm512_add_pd(x.m_re,v),
	                    x.m_im)});
}

static inline gms::math::AVX512c8f64
gms::math::operator+(const __m512d v,
                     const AVX512c8f64 x) {
        return (AVX512c8f64{_mm512_add_pd(v,x.m_re),
	                   x.m_im});
}

static inline gms::math::AVX512c8f64
gms::math::operator+(const AVX512c8f64 x,
		    const double s) {
	return (x + AVX512c8f64{ s });
}

static inline gms::math::AVX512C8f64
gms::math::operator+(const double s,
		    const AVX512c8f64 x) {
	return (AVX512c8f64{ s } +x);
}

static inline gms::math::AVX512c8f64
gms::math::operator+=(AVX512c8f64 x,
		     const AVX512c8f64 y) {
	x = x + y;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator+=(AVX512c8f64 x,
                      const __m512d v) {
        x = x + v;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator+=(const __m512d v,
                      AVX512c8f64 x) {
        x = v + x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator+=(AVX512c8f64 x,
		     const double s) {
	x = x + AVX512c8f64{s};
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator+=(const double s,
		      AVX512c8f64 x) {
	x = AVX512c8f64{s} + x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-(const AVX512c8f64 x,
		     const AVX512c8f64 y) {
	return (AVX512c8f64{ _mm512_sub_pd(x.m_re, y.m_re),
		_mm512_sub_pd(x.m_im, y.m_im) });
}

static inline gms::math::AVX512c8f64
gms::math::operator-(const AVX512c8f64 x,
                     const __m512d v ) {
        return (AVX512c8f64{_mm512_sub_pd(x.m_re,v),
	                   x.m_im});
}

static inline gms::math::AVX512c8f64
gms::math::operator-(const __m512d v,
                     const AVX512c8f64 x) {
        return (AVX512c8f64{_mm512_sub_pd(v,x.m_re),
	                   x.m_im});
}

static inline gms::math::AVX512c8f64
gms::math::operator-=(AVX512c8f64 x,
                      const __m512d v) {
        x = x - v;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-=(const __m512d v,
                      AVX512c8f64 x) {
        x = v - x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-(const AVX512c8f64 x,
		    const double s) {
	return (x - AVX512c8f64{ s });
}

static inline gms::math::AVX512c8f64
gms::math::operator-(const double s,
		     const AVX512c8f64 x) {
	return(AVX512c8f64{ s } - x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-(AVX512c8f64 x) {
	x =  gms::math::AVX512c8f64::CZERO - x; 
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-=(AVX512c8f64 x,
		      const AVX512c8f64 y) {
	x = x - y;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-=(AVX512c8f64 x,
		      const double s) {
	x = x - AVX512c8f64{s};
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator-=(const double s,
			 AVX512c8f64 x) {
	x = AVX512c8f64{s} - x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator*(const AVX512c8f64 x,
		     const AVX512c8f64 y) {
	const __m512d ymm0(_mm512_mul_pd(x.m_re,y.m_re));
	const __m512d ymm1(_mm512_mul_pd(x.m_im,y.m_im));
	const __m512d ymm2(_mm512_mul_pd(x.m_im,y.m_re));
	const __m512d ymm3(_mm512_mul_pd(x.m_re, y.m_im));
	return (AVX512c8f64{ _mm512_sub_pd(ymm0,ymm1),
			_mm512_add_pd(ymm2,ymm3) });
}

static inline gms::math::AVX512c8f64
gms::math::operator*(const AVX512c8f64 x,
                     const __m512d v) {
        return (AVX512c8f64{_mm512_mul_pd(x.m_re,v),
	                    _mm512_mul_pd(x.m_im,v)});
}

static inline gms::math::AVX512c8f64
gms::math::operator*(const __m512d v,
                     const AVX512c8f64 x) {
        return (AVX512c8f64{_mm512_mul_pd(v,x.m_re),
	                    _mm512_mul_pd(v,x.m_im)});
}

static inline gms::math::AVX512c8f64
gms::math::operator*(const AVX512c8f64 x,
		    const double s) {
	const __m512d zmm0(_mm512_set1_pd(s));
	return (AVX512c8f64{_mm512_mul_pd(x.m_re,zmm0),
				      _mm512_mul_pd(x.m_im,zmm0)});
}

static inline gms::math::AVX512c8f64
gms::math::operator*(const double s,
		     const AVX512c8f64 x) {
	const __m512d zmm0(_mm512_set1_pd(s));
	return (AVX512c8f64{ _mm512_mul_pd(zmm0,x.m_re),
					   _mm512_mul_pd(zmm0,x.m_im) });
}

static inline gms::math::AVX512c8f64
gms::math::operator*=(AVX512c8f64 x,
		     const AVX512c8f64 y) {
	x = x * y;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator*=(AVX512c8f64 x,
                      const __m512d v) {
        x = x * v;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator*=(const __m512d v,
                      AVX512c8f64 x) {
       x = v * x;
       return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator*=(AVX512c8f64 x,
			const double s) {
         
	x = x * s;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator*=(const double s,
		      AVX512c8f64 &x) {
	x = s * x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator/(const AVX512c8f64 &,
		     const AVX512c8f64 y) {
	const __m512d ymm0(_mm512_mul_pd(x.m_re, y.m_re));
	const __m512d ymm1(_mm512_mul_pd(x.m_im, y.m_im));
	const __m512d ymm2(_mm512_mul_pd(x.m_im, y.m_re));
	const __m512d ymm3(_mm512_mul_pd(x.m_re, y.m_im));
	const __m512d den(_mm512_add_pd(_mm512_mul_pd(y.m_re,y.m_re),
		_mm512_mul_pd(y.m_im, y.m_im)));
	const __m512d re_part(_mm512_add_pd(ymm0,ymm1));
	const __m512d im_part(_mm512_sub_pd(ymm2, ymm3));
	return (AVX512c8f64{ _mm512_div_pd(re_part,den),
		_mm512_div_pd(im_part,den) });
}

static inline gms::math::AVX512c8f64
gms::math::operator/(const AVX512c8f64 x,
                     const __m512d v) {
        return (AVX512c8f64{_mm512_div_pd(x.m_re,v),
	                    _mm512_div_pd(x.m_im,v)});
}

static inline gms::math::AVX512c8f64
gms::math::operator/(const __m512d v,
                     const AVX512c8f64 x) {
        return (AVX512c8f64{_mm512_div_pd(v,x.m_re),
	                    _mm512_div_pd(v,x.m_im)});
}

static inline gms::math::AVX512c8f64
gms::math::operator/(const AVX512c8f64 x,
			 const double s) {
	const __m512d zmm0(_mm512_set1_pd(s));
	return (AVX512c8f64{_mm512_div_pd(x.m_re,zmm0),
			_mm512_div_pd(x.m_im,zmm0)});
}

static inline gms::math::AVX512c8f64
gms::math::operator/(const double s,
		   const AVX512c8f64 x) {
	const __m512d zmm0(_mm512_set1_pd(s));
	return (AVX512c8f64{ _mm512_div_pd(zmm0,x.m_re),
		_mm512_div_pd(zmm0,x.m_im) });
}

static inline gms::math::AVX512c8f64
gms::math::operator/=(AVX512c8f64 x,
                      const __m512d v) {
        x = x / v;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator/=(const __m512d v,
                      AVX512c8f64 x) {
        x = v / x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator/=(AVX512c8f64 x,
			const AVX512c8f64 y) {
	x = x / y;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator/=(AVX512c8f64 x,
			 const double s) {
	x = x / s;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator/=(const double s,
			AVX512c8f64 x) {
	x = s / x;
	return (x);
}

static inline gms::math::AVX512c8f64
gms::math::operator~(AVX512c8f64 x) {
	x.m_re = _mm512_sub_pd(_mm512_setzero_pd(),x.m_re);
	return (x);
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator==(const AVX512c8f64 x,
			const AVX512c8f64 y) {
	__mmask8 m1(_mm512_cmp_pd_mask(x.m_re,y.m_re,_CMP_EQ_OQ));
	__mmask8 m2(_mm512_cmp_pd_mask(x.m_im,y.m_im,_CMP_EQ_OQ));
	return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator==(const AVX512c8f64 x,
                      std::complex<double> c) {
       __mmask8 m1(_mm512_cmp_pd_mask(x.m_re,_mm512_set1_pd(c.real()),_CMP_EQ_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(x.m_im,_mm512_set1_pd(c.imag()),_CMP_EQ_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator==(const std::complex<double> c,
                      const AVX512c8f64 x) {
       __mmask8 m1(_mm512_cmp_pd_mask(_mm512_set1_pd(c.real(),x.m_re),_CMP_EQ_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(_mm512_set1_pd(c.imag(),x.m_im),_CMP_EQ_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator!=(const AVX512c8f64 x,
		      const AVX512c8f64 y) {
	__mmask8 m1(_mm512_cmp_pd_mask(x.m_re, y.m_re, _CMP_NEQ_OQ));
	__mmask8 m2(_mm512_cmp_pd_mask(x.m_im, y.m_im, _CMP_NEQ_OQ));
	return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator!=(const AVX512c8f64 x,
                      const std::complex<double> c) {
       __mmask8 m1(_mm512_cmp_pd_mask(x.m_re,_mm512_set1_pd(c.real()),_CMP_NEQ_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(x.m_im,_mm512_set1_pd(c.imag()),_CMP_NEQ_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator!=(const std::complex<double> c,
                      const AVX512c8f64 x) {
       __mmask8 m1(_mm512_cmp_pd_mask(_mm512_set1_pd(c.real(),x.m_re),_CMP_NEQ_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(_mm512_set1_pd(c.imag(),x.m_im),_CMP_NEQ_OQ));
       return (std::make_pair(m1,m2));
}



static inline std::pair<__mmask8,__mmask8>
gms::math::operator>(const AVX512c8f64 x,
		     const AVX512c8f64 y) {
	__mmask8 m1(_mm512_cmp_pd_mask(x.m_re, y.m_re, _CMP_GT_OQ));
	__mmask8 m2(_mm512_cmp_pd_mask(x.m_im, y.m_im, _CMP_GT_OQ));
	return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator>(const AVX512c8f64 x,
                     const std::complex<double> c) {
       __mmask8 m1(_mm512_cmp_pd_mask(x.m_re,_mm512_set1_pd(c.real()),_CMP_GT_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(x.m_im,_mm512_set1_pd(c.imag()),_CMP_GT_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator>(const std::complex<double> c,
                     const AVX512c8f64 x) {
       __mmask8 m1(_mm512_cmp_pd_mask(_mm512_set1_pd(c.real(),x.m_re),_CMP_GT_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(_mm512_set1_pd(c.imag(),x.m_im),_CMP_GT_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator<(const AVX512c8f64 x,
		     const AVX512c8f64 y) {
	__mmask8 m1(_mm512_cmp_pd_mask(x.m_re, y.m_re, _CMP_LT_OQ));
	__mmask8 m2(_mm512_cmp_pd_mask(x.m_im, y.m_im, _CMP_LT_OQ));
	return (std::make_pair(m1, m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator<(const AVX512c8f64 x,
                     const std::complex<double> c) {
       __mmask8 m1(_mm512_cmp_pd_mask(x.m_re,_mm512_set1_pd(c.real()),_CMP_LT_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(x.m_im,_mm512_set1_pd(c.imag()),_CMP_LT_OQ));
       return (std::make_pair(m1,m2)); 
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator<(const std::complex<double> c,
                     const AVX512c8f64 x) {
       __mmask8 m1(_mm512_cmp_pd_mask(_mm512_set1_pd(c.real(),x.m_re),_CMP_LT_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(_mm512_set1_pd(c.imag(),x.m_im),_CMP_LT_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8, __mmask8>
gms::math::operator>=(const AVX512c8f64 x,
		      const AVX512c8f64 y) {
	__mmask8 m1(_mm512_cmp_pd_mask(x.m_re, y.m_re, _CMP_GE_OQ));
	__mmask8 m2(_mm512_cmp_pd_mask(x.m_im, y.m_im, _CMP_GE_OQ));
	return (std::make_pair(m1, m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator>=(const AVX512c8f64 x,
                      const std::complex<double> c) {
       __mmask8 m1(_mm512_cmp_pd_mask(x.m_re,_mm512_set1_pd(c.real()),_CMP_GE_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(x.m_im,_mm512_set1_pd(c.imag()),_CMP_GE_OQ));
       return (std::make_pair(m1,m2)); 
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator>=(const std::complex<double> c,
                      const AVX512c8f64 x) {
       __mmask8 m1(_mm512_cmp_pd_mask(_mm512_set1_pd(c.real(),x.m_re),_CMP_GE_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(_mm512_set1_pd(c.imag(),x.m_im),_CMP_GE_OQ));
       return (std::make_pair(m1,m2));
}

static inline std::pair<__mmask8, __mmask8>
gms::math::operator<=(const AVX512c8f64 x,
                      const AVX512c8f64 y) {
	__mmask8 m1(_mm512_cmp_pd_mask(x.m_re, y.m_re, _CMP_LE_OQ));
	__mmask8 m2(_mm512_cmp_pd_mask(x.m_im, y.m_im, _CMP_LE_OQ));
	return (std::make_pair(m1, m2));
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator<=(const AVX512c8f64 x,
                      const std::complex<double> c) {
       __mmask8 m1(_mm512_cmp_pd_mask(x.m_re,_mm512_set1_pd(c.real()),_CMP_LE_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(x.m_im,_mm512_set1_pd(c.imag()),_CMP_LE_OQ));
       return (std::make_pair(m1,m2)); 
}

static inline std::pair<__mmask8,__mmask8>
gms::math::operator<=(const std::complex<double> c,
                      const AVX512c8f64 x) {
       __mmask8 m1(_mm512_cmp_pd_mask(_mm512_set1_pd(c.real(),x.m_re),_CMP_LE_OQ));
       __mmask8 m2(_mm512_cmp_pd_mask(_mm512_set1_pd(c.imag(),x.m_im),_CMP_LE_OQ));
       return (std::make_pair(m1,m2));
}



