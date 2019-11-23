
const gms::math::AVX512c4f32
gms::math::AVX512c4f32::C4F32_ZERO = AVX512c4f32{};

gms::math::AVX512c4f32
::AVX512c4f32() {
   m_re = _mm512_setzero_ps();
   m_im = _mm512_setzero_ps();
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const AVX512c4Payload x) {
   m_re = _mm512_set_ps(x.re0,x.re1,x.re2,x.re3,x.re4,
                        x.re5,x.re6,x.re7,x.re8,x.re9,
			x.re10,x.re11,x.re12,x.re13,
			x.re14,x.re15);
   m_im = _mm512_set_ps(x.im0,x.im1,x.im2,x.im3,x.im4,
                        x.im5,x.im6,x.im7,x.im8,x.im9,
			x.im10,x.im11,x.im12,x.im13,
			x.im14,x.im15);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const float * __restrict re,
              const float * __restrict im) {
   m_re = _mm512_load_ps(&re[0]);
   m_im = _mm512_load_ps(&im[0]);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const float re,
              const float im) {
   m_re = _mm512_set1_ps(re);
   m_im = _mm512_set1_ps(im);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const std::complex<float> c) {
   m_re = _mm512_set1_ps(c.real());
   m_im = _mm512_set1_ps(c.imag());
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const float re) {
   m_re = _mm512_set1_ps(re);
   m_im = _mm512_setzero_ps();
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const float re0,
              const float re1,
	      const float re2,
	      const float re3,
	      const float re4,
	      const float re5,
	      const float re6,
	      const float re7,
	      const float re8
	      const float re9,
	      const float re10,
	      const float re11,
	      const float re12,
	      const float re13,
	      const float re14,
	      const float re15) {
     m_re = _mm512_set_ps(re0,re1,re2,re3,re4,
                          re5,re6,re7,re8,re9,
			  re10,re11,re12,re13,
			  re14,re15);
     m_im = _mm512_setzero_ps();
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const __m512 re,
              const __m512 im) {
     m_re = re;
     m_im = im;
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32
::AVX512c4f32(const AVX512c4f32 &x) {
     m_re = x.m_re;
     m_im = x.m_im;
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::load_a(const float * __restrict re,
                               const float * __restrict im) {
     m_re = _mm512_load_ps(&re[0]);
     m_im = _mm512_load_ps(&im[0]);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::load_u(const float * __restrict re,
                               const float * __restrict im) {
     m_re = _mm512_loadu_ps(&re[0]);
     m_im = _mm512_loadu_ps(&im[0]);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

void
gms::math::AVX512c4f32::store_a(float * __restrict re,
                                float * __restrict im) {
     _mm512_store_ps(&re[0],m_re);
     _mm512_store_ps(&im[0],m_im);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

void
gms::math::AVX512c4f32::store_u(float * __restrict re,
                                float * __restrict im) {
     _mm512_storeu_ps(&re[0],m_re);
     _mm512_storeu_ps(&im[0],m_im);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

void
gms::math::AVX512c4f32::stream_nt(float * __restrict re,
                                  float * __restrict im) {
     _mm512_stream_ps(&re[0],m_re);
     _mm512_stream_ps(&im[0],m_im);
     _mm_sfence();
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

float
gms::math::AVX512c4f32::extract_1f32(const int32_t idx) {
#if defined _WIN64
    __declspec(align(64))        float mem[32] = {};
#elif defined __linux
    __attribute__((aligned(64))) float mem[32] = {};
#endif
    store_a(&mem[0],&mem[16]);
    return (mem[idx & 0x1F]);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

std::pair<float,float>
gms::math::AVX512c4f32::extract_2f32(const int32_t re_idx,
                                     const int32_t im_idx) {
#if defined _WIN64
    __declspec(align(64)) float re_mem[16] = {};
    __declspec(align(64)) float im_mem[16] = {};
#elif defined __linux
    __attribute__((aligned(64))) float re_mem[16] = {};
    __attribute__((aligned(64))) float im_mem[16] = {};
#endif
    store_a(&re_mem[0],&im_mem[0]);
    return (std::make_pair(re_mem[re_idx & 0x1F],im_mem[im_idx & 0x1F]));
} 
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::insert_1f32(const int32_t idx,
                                    const float value) {
#if defined _WIN64
    __declspec(align(64))        float mem[32] = {};
#elif defined __linux
    __attribute__((aligned(64))) float mem[32] = {};
#endif
    store_a(&mem[0],&mem[16]);
    mem[idx & 0x1F] = value;
    load_a(&mem[0],&mem[16]);
    return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::insert_2f32(const int32_t re_idx,
                                    const int32_t im_idx,
				    const float re,
				    const float im) {
#if defined _WIN64
    __declspec(align(64)) float mem_re[16] = {};
    __declspec(align(64)) float mem_im[16] = {};
#elif defined __linux
    __attribute__((aligned(64))) float mem_re[16] = {};
    __attribute__((aligned(64))) float mem_im[16] = {};
#endif
    store_a(&mem_re[0],&mem_im[0]);
    mem_re[re_idx & 0x1F] = re;
    mem_im[im_idx & 0x1F] = im;
    load_a(&mem_re[0],&mem_im[0]);
    return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

void
gms::math::AVX512c4f32::concatenate_a(float * __restrict out) {
     store_a(&out[0],&out[16]);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

void
gms::math::AVX512c4f32::concatenate_u(float * __restrict out) {
     store_u(&out[0],&out[16]);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::partial_loadu(const float const * __restrict re,
                                      const int32_t n_re,
				      const float const * __restrict im,
				      const int32_t n_im) {
     m_re = _mm512_maskz_loadu_ps(__mmask16((1 << n_re)-1),re);
     m_im = _mm512_maskz_loadu_ps(__mmask16((1 << n_im)-1),im);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::partial_loada(const float const * __restrict re,
                                      const int32_t n_re,
				      const float const * __restrict im,
				      const int32_t n_im) {
     m_re = _mm512_maskz_load_ps(__mmask16((1 << n_re)-1),re);
     m_im = _mm512_maskz_load_ps(__mmask16((1 << n_im)-1),im);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

void
gms::math::AVX512c4f32::partial_storeu(float * __restrict re,
                                       const int32_t n_re,
				       float * __restrict im,
				       const int32_t n_im) {
     _mm512_mask_storeu_ps(&re[0],__mmask16((1 << n_re)-1),m_re);
     _mm512_mask_storeu_ps(&im[0],__mmask16((1 << m_im)-1),m_im);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

void
gms::math::AVX512c4f32::partial_storea(float * __restrict re,
                                       const int32_t n_re,
				       float * __restrict im,
				       const int32_t n_im) {
     _mm512_mask_store_ps(&re[0],__mmask16((1 << n_re)-1),m_re);
     _mm512_mask_store_ps(&im[0],__mmask16((1 << m_im)-1),m_im);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::expand(const AVX512c4f32 x,
                               const __mmask16 mask) {
     m_re = _mm512_maskz_expand_ps(mask,x.m_re);
     m_im = _mm512_maskz_expand_ps(mask.x.m_im);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::expand_loadu(const AVX512c4f32 x,
                                    const __mmask16 mask,
				    const float * __restrict re,
				    const float * __restrict im) {
     m_re = _mm512_mask_expandloadu_ps(x.m_re,mask,&re[0]);
     m_im = _mm512_mask_expandloadu_ps(x.m_im,mask,&im[0]);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::permute(const __mmask16 mask,
                                const int32_t imm) {
     m_re = _mm512_mask_permute_ps(m_re,mask,m_im,imm);
     m_im = _mm512_mask_permute_ps(m_im,mask,m_re,imm);
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

__m256
gms::math::AVX512c4f32::re_low2() const {
     return (_mm512_extractf32x8_ps(m_re,0));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

__m256
gms::math::AVX512c4f32::re_hi2() const {
     return (_mm512_extract32x8_ps(m_re,1));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

__m256
gms::math::AVX512c4f32::im_low2() const {
     return (_mm512_extract32x8_ps(m_im,0));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

__m256
gms::math::AVX512c4f32::im_hi2() const {
     return (_mm512_extract32x8_ps(m_im,1));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

gms::math::AVX512c4f32 &
gms::math::AVX512c4f32::operator=(const AVX512c4f32 x) {
     if(this == &x) return (*this);
     m_re = x.m_re;
     m_im = x.m_im;
     return (*this);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::conj(const AVX512c4f32 x) {
     auto tmp = ~x;
     const __m512 im_part = tmp.m_im;
     return (AVX512c4f32{x.m_re,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::polar(const __m512 ro,
                 const __m512 theta) {
     const __m512 re_part =
           _mm512_mul_ps(rho,_mm512_cos_ps(theta));
     const __m512 im_part =
           _mm512_mul_ps(rho,_mm512_sin_ps(theta));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline __m512
gms::math::carg(const AVX512c4f32 x) {
     return (_mm512_atan2_ps(x.m_re,x.m_im));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline __m512
gms::math::carg(const float re,
                const float im) {
     const __m512 real =  AVX512C4F32_SETPS(re)
     const __m512 imag =  AVX512C4F32_SETPS(im)
     return (_mm512_atan2_ps(real,imag));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csin(const AVX512c4f32 x) {
     const __m512 re_part =
           _mm512_mul_ps(_mm512_sin_ps(x.m_re),
	                 _mm512_cosh_ps(x.m_im));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_cos_ps(x.m_re),
	                 _mm512_sinh_ps(x.m_im));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csin(const std::complex<float> c) {
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(_mm512_sin_ps(real),
	                 _mm512_cosh_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_cos_ps(real),
	                 _mm512_sinh_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csin(const float re,
                const float im) {
     const __m512 real = AVX512C4F32_SETPS(re)
     const __m512 imag = AVX512C4F32_SETPS(im)
     const __m512 re_part =
           _mm512_mul_ps(_mm512_sin_ps(real),
	                 _mm512_cosh_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_cos_ps(real),
	                 _mm512_sinh_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csinh(const AVX512c4f32 x) {
     const __m512 re_part =
           _mm512_mul_ps(_mm512_sinh_ps(x.m_re),
	                 _mm512_cos_ps(x.m_im));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_cosh_ps(x.m_re),
	                 _mm512_sin_ps(x.m_im));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csinh(const std::complex<float> c) {
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(_mm512_sinh_ps(real),
	                 _mm512_cos_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_cosh_ps(real),
	                 _mm512_sin_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csinh(const float re,
                 const float im) {
     const __m512 real = AVX512C4F32_SETPS(re)
     const __m512 imag = AVX512C4F32_SETPS(im)
     const __m512 re_part =
           _mm512_mul_ps(_mm512_sinh_ps(real),
	                 _mm512_cos_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_cosh_ps(real),
	                 _mm512_sin_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ccos(const AVX512c4f32 x) {
     const __m512 re_part =
           _mm512_mul_ps(_mm512_cos_ps(x.m_re),
	                 _mm512_cosh_ps(x.m_im));
     const __m512 im_part =
           _mm512_mul_part(_mm512_sin_ps(x.m_re),
	                   _mm512_sinh_ps(x.m_im));
      return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ccos(const std::complex<float> c) {
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(_mm512_cos_ps(real),
	                 _mm512_cosh_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_sin_ps(real),
	                 _mm512_sinh_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ccos(const float re,
                const float im) {
     const __m512 real = AVX512C4F32_SETPS(re)
     const __m512 imag = AVX512C4F32_SETPS(im)
     const __m512 re_part =
           _mm512_mul_ps(_mm512_cos_ps(real),
	                 _mm512_cosh_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_sin_ps(real),
	                 _mm512_sinh_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ccosh(const AVX512c4f32 x) {
     const __m512 re_part =
           _mm512_mul_ps(_mm512_cosh_ps(x.m_re),
	                 _mm512_cos_ps(x.m_im));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_sinh_ps(x.m_re),
	                 _mm512_sin_ps(x.m_im));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ccosh(const std::complex<float> c) {
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(_mm512_cosh_ps(real),
	                 _mm512_cos_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_sinh_ps(real),
	                 _mm512_sin_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ccosh(const float re,
                 const float im) {
     const __m512 real = AVX512C4F32_SETPS(re)
     const __m512 imag = AVX512C4F32_SETPS(im)
     const __m512 re_part =
           _mm512_mul_ps(_mm512_cosh_ps(real),
	                 _mm512_cos_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_sinh_ps(real),
	                 _mm512_sin_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::cexp(const AVX512c4f32 x) {
     const __m512 re_part =
           _mm512_mul_ps(_mm512_exp_ps(x.m_re),
	                 _mm512_cos_ps(x.m_im));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_exp_ps(x.m_re),
	                 _mm512_sin_ps(x.m_im));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::cexp(const std::complex<float> c) {
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(_mm512_exp_ps(real),
	                 _mm512_cos_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_exp_ps(real),
	                 _mm512_sin_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::cexp(const float re,
                const float im) {
     const __m512 real = AVX512C4F32_SETPS(re)
     const __m512 imag = AVX512C4F32_SETPS(im)
     const __m512 re_part =
           _mm512_mul_ps(_mm512_exp_ps(real),
	                 _mm512_cos_ps(imag));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_exp_ps(real),
	                 _mm512_sin_ps(imag));
     return (AVX512c4f32{re_part,im_part});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline __m512
gms::math::cabs(const AVX512c4f32 x) {
     const __m512 re_part =
           _mm512_mul_ps(x.m_re,x.m_re);
     const __m512 im_part =
           _mm512_mul_ps(x.m_im,x.m_im);
     return (_mm512_sqrt_ps(_mm512_add_ps(re_part,im_part)));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline __m512
gms::math::cabs(const std::complex<float> c) {
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(real,real);
     const __m512 im_part =
           _mm512_mul_ps(imag,imag);
     return (_mm512_sqrt_ps(_mm512_add_ps(re_part,im_part)));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline __m512
gms::math::cabs(const float re,
                const float im) {
     const __m512 real = AVX512C4F32_SETPS(real)
     const __m512 imag = AVX512C4F32_SETPS(imag)
     const __m512 re_part =
           _mm512_mul_ps(real,real);
     const __m512 im_part =
           _mm512_mul_ps(imag,imag);
     return (_mm512_sqrt_ps(_mm512_add_ps(re_part,im_part)));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::cpow(const AVX512c4f32 x,
                const float n) {
     const __m512 re_part =
           _mm512_mul_ps(x.m_re,x.m_re);
     const __m512 im_part =
           _mm512_mul_ps(x.m_im,x.m_im);
     const __m512 r =
           _mm512_sqrt_ps(_mm512_add_ps(re_part,im_part));
     const __m512 theta =
           _mm512_atan_ps(_mm512_div_ps(x.m_im,x.m_re));
     const __m512 vn = _mm512_set1_ps(n);
     const __m512 pow_term = _mm512_pow_ps(r,vn);
     const __m512 trig_arg = _mm512_mul_ps(vn,theta);
     return (AVX512c4f32{_mm512_mul_ps(pow_term,_mm512_cos_ps(trig_arg),
                         _mm512_mul_ps(pow_term,_mm512_sin_ps(trig_arg))}));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::clog(const AVX512c4f32 x) {
    const __m512 t1  = cabs(x);
    const __m512 t2  = carg(x);
    const __m512 re_part = _mm512_log_ps(t1);
    return (AVX512c4f32{re_part,t1});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::clog(const std::complex<float> c) {
     const __m512 t1 = cabs(c);
     const __m512 t2 = carg(c.real(),c.imag());
     const __m512 re_part = _mm512_log_ps(t1);
     return (AVX512c4f32{re_part,t1});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::clog(const float re,
                const float im) {
     const __m512 t1 = cabs(re,im);
     const __m512 t2 = carg(re,im);
     const __m512 re_part = _mm512_log_ps(t1);
     return (AVX512c4f32{re_part,t1});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csqrt(const AVX512c4f32 x) {
     const __m512 t = cabs(x);
     const __m512 re_part =
           _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                 _mm512_add_ps(t,x.m_re));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                 _mm512_sub_ps(t,x.m_re));
     return (AVX512c4f32{_mm512_sqrt_ps(re_part),
                         _mm512_sqrt_ps(im_part)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csqrt(const std::complex<float> c) {
     const __m512 t = cabs(c);
     const __m512 real = AVX512C4F32_SETPS(c.real())
     const __m512 imag = AVX512C4F32_SETPS(c.imag())
     const __m512 re_part =
           _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                 _mm512_add_ps(t,real));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                 _mm512_sub_ps(t,real));
     return (AVX512c4f32{_mm512_sqrt_ps(re_part),
                         _mm512_sqrt_ps(im_part)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::csqrt(const float re,
                 const float im) {
     const __m512 t = cabs(c);
     const __m512 real = AVX512C4F32_SETPS(re)
     const __m512 imag = AVX512C4F32_SETPS(im)
     const __m512 re_part =
           _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                 _mm512_add_ps(t,real));
     const __m512 im_part =
           _mm512_mul_ps(_mm512_set1_ps(0.5f),
	                 _mm512_sub_ps(t,real));
     return (AVX512c4f32{_mm512_sqrt_ps(re_part),
                         _mm512_sqrt_ps(im_part)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ctan(const AVX512c4f32 x) {
     return (csin(x)/ccos(x));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ctan(const std::complex<float> c) {
     return (csin(c)/ccos(c));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ctan(const float re,
                const float im) {
     return (csin(re,im)/ccos(re,im));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ctanh(const AVX512c4f32 x) {
     return (csinh(x)/ccosh(x));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ctanh(const std::complex<float> c) {
     return (csinh(c)/ccosh(c));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::ctanh(const float re,
                 const float im) {
     return (csinh(re,im)/ccosh(re,im));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::select(const AVX512c4f32 x,
                  const AVX512c4f32 y,
		  __mmask16 mask) {
     return (AVX512c4f32{_mm512_mask_blend_ps(mask,x.m_re,x.m_im),
                         _mm512_mask_blend_ps(mask,x.m_im,x.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::cdiv_smith(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     __m512 ratio,denom,re_part,im_part;
     constexpr __mmask16 all_ones = 0xFFFF;
     re_part = _mm512_setzero_ps();
     im_part = _mm512_setzero_ps();
     __mmask16 is_gte = _mm512_cmp_ps_mask(
                              _mm512_abs_ps(y.m_re),
			      _mm512_abs_ps(y.m_im),_CMP_GE_OQ);
     ratio = _mm512_setzero_ps();
     denom = _mm512_setzero_ps();
     if(is_gte == all_ones) {
          ratio = _mm512_div_ps(y.m_im,y.m_re);
	  denom = _mm512_add_ps(y.m_re,
	                _mm512_mul_ps(ratio,y.m_im));
	  re_part = _mm512_div_ps(_mm512_add_ps(x.m_re,
	                _mm512_mul_ps(x.m_im,ratio)),denom);
	  im_part = _mm512_div_ps(_mm512_sub_ps(x.m_im,
	                _mm512_mul_ps(x.m_re,ratio)),denom);
	  return (AVX512c4f32{re_part,im_part});
      }
      else {
          ratio = _mm512_div_ps(y.m_re,y.m_im);
	  denom = _mm512_add_ps(y.m_im,_mm512_mul_ps(ratio,y.m_re));
	  re_part = _mm512_div_ps(_mm512_add_ps(
	                          _mm512_mul_ps(x.m_re,ratio)),denom);
	  im_part = _mm512_div_ps(_mm512_sub_ps(
	                          _mm512_mul_ps(x.m_im,ratio)),denom);
	  return (AVX512c4f32{re_part,im_part});
      }
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+(const AVX512c4f32 x,
                     const AVX512c4f32 y) {
     return (AVX512c4f32{_mm512_add_ps(x.m_re,y.m_re),
                         _mm512_add_ps(x.m_im,x.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+(const AVX512c4f32 x,
                     const __m512 v) {
     return (AVX512c4f32{_mm512_add_ps(x.m_re,v),
                         x.m_im});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+(const __m512 v,
                     const AVX512c4f32 x) {
     return (AVX512c4f32{_mm512_add_ps(v,x.m_re),x.m_im});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+(const AVX512c4f32 x,
                     const float s) {
     return (x + AVX512c4f32{s});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+(const float s,
                     const AVX512c4f32 x) {
     return (AVX512c4f32{s} + x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+=(AVX512c4f32 x,
                      const AVX512c4f32 y) {
     x = x + y;
     return (x)
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+=(AVX512c4f32 x,
                      const __m512 v) {
     x = x + v;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+=(const __m512 v,
                      AVX512c4f32 x) {
     x = v + x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+=(AVX512c4f32 x,
                      const float s) {
     x = x + AVX512c4f32{s};
     return (x)
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator+=(const float s,
                      AVX512c4f32 x) {
     x = AVX512c4f32{s} + x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-(const AVX512c4f32 x,
                     const AVX512c4f32 y) {
     return (AVX512c4f32{_mm512_sub_ps(x.m_re,y.m_re),
                         _mm512_sub_ps(x.m_im,y.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-(const AVX512c4f32 x,
                     const __m512 v) {
     return (AVX512c4f32{_mm512_sub_ps(x.m_re,v),x.m_im});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-(const __m512 v,
                     const AVX512c4f32 x) {
     return (AVX512c4f32{_mm512_sub_ps(v,x.m_re),x.m_im});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-(const AVX512c4f32 x,
                     const float s) {
     return (x + AVX512c4f32{s});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-(const float s,
                     const AVX512c4f32 x) {
     return (AVX512c4f32{s} - x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-(AVX512c4f32 x) {
     x = gms::math::AVX512c4f32::C4F32_ZERO - x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-=(AVX512c4f32 x,
                      const AVX512c4f32 y) {
     x = x - y;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-=(AVX512c4f32 x,
                      const __m512 v) {
     x = x - v;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-=(const __m512 v,
                      AVX512c4f32 x) {
     x = v - x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-=(AVX512c4f32 x,
                      const float s) {
     x = x - AVX512c4f32{s};
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator-=(const float s,
                      AVX512c4f32 x) {
     x = AVX512C4F32{s} - x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*(const AVX512c4f32 x,
                     const AVX512c4f32 y) {
     const __m512 zmm0 = _mm512_mul_ps(x.m_re,y.m_re);
     const __m512 zmm1 = _mm512_mul_ps(x.m_im,y.m_im);
     const __m512 zmm2 = _mm512_mul_ps(x.m_im,y.m_re);
     const __m512 zmm3 = _mm512_mul_ps(x.m_re,y.m_im);
     return (AVX512c4f32{_mm512_sub_ps(zmm0,zmm1),
                         _mm512_sub_ps(zmm2,zmm3)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*(const AVX512c4f32 x,
                     const __m512 v) {
     return (AVX512c4f32{_mm512_mul_ps(x.m_re,v),
                         _mm512_mul_ps(x.m_im,v)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*(const __m512 v,
                     const AVX512c4f32 x) {
     return (AVX512c4f32{_mm512_mul_ps(v,x.m_re),
                         _mm512_mul_ps(v,x.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*(const AVX512c4f32 x,
                     const float s) {
     const __m512 zmm0 = AVX512C4F32_SETPS(s)
     return (AVX512c4f32{_mm512_mul_ps(x.m_re,zmm0),
                         _mm512_mul_ps(x.m_im,zmm0)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*(const float s,
                     const AVX512c4f32 x) {
     const __m512 zmm0 = AVX512C4F32_SETPS(s)
     return (AVX512c4f32{_mm512_mul_ps(zmm0,x.m_re),
                         _mm512_mul_ps(zmm0,x.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*=(AVX512c4f32 x,
                      const AVX512c4f32 y) {
     x = x * y;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*=(AVX512c4f32 x,
                      const __m512 v) {
     x = x * v;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*=(const __m512 v,
                      AVX512c4f32 x) {
     x = v * x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*=(AVX512c4f32 x,
                      const float s) {
     x = x * s;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator*=(const float s,
                      AVX512c4f32 x) {
     x = s * x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/(const AVX512c4f32 x,
                     const AVX512c4f32 y) {
#if defined USE_SAFE_COMPLEX_DIVISION && (USE_SAFE_COMPLEX_DIVISION) == 1
    return (cdiv_smith(x,y));
#else
    const __m512 zmm0 = _mm512_mul_ps(x.m_re,y.m_re);
    const __m512 zmm1 = _mm512_mul_ps(x.m_im,y.m_im);
    const __m512 zmm2 = _mm512_mul_ps(x.m_im,y.m_re);
    const __m512 zmm3 = _mm512_mul_ps(x.m_re,y.m_im);
    const __m512 den  = _mm512_add_ps(_mm512_mul_ps(y.m_re,y.m_re),
                                      _mm512_mul_ps(y.m_im,y.m_im));
    const __m512 re_part = _mm512_add_ps(zmm0,zmm1);
    const __m512 im_part = _mm512_sub_ps(zmm2,zmm3);
    return (AVX512c4f32{_mm512_div_ps(re_part,den),
                        _mm512_div_ps(im_part,den)});
#endif
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/(const AVX512c4f32 x,
                     const __m512 v) {
     return (AVX512c4f32{_mm512_div_ps(x.m_re,v),
                         _mm512_div_ps(x.m_im,v)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/(const __m512 v,
                     const AVX512c4f32 x) {
     return (AVX512c4f32{_mm512_div_ps(v,x.m_re),
                         _mm512_div_ps(v,x.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/(const AVX512c4f32 x,
                     const float s) {
     const __m512 zmm0 = AVX512C4F32_SETPS(s)
     return (AVX512c4f32{_mm512_div_ps(x.m_re,zmm0),
                         _mm512_div_ps(x.m_im,zmm0)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/(const float s,
                     const AVX512c4f32 x) {
     const __m512 zmm0 = AVX512C4F32_SETPS(s)
     return (AVX512c4f32{_mm512_div_ps(zmm0,x.m_re),
                         _mm512_div_ps(zmm0,x.m_im)});
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/=(AVX512c4f32 x,
                      const AVX512c4f32 y) {
     x = x / y;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/=(AVX512c4f32 x,
                      const __m512 v) {
     x = x / v;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/=(const __m512 v,
                      AVX512c4f32 x) {
     x = v / x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/=(AVX512c4f32 x,
                      const float s) {
     x = x / s;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator/=(const float s,
                      AVX512c4f32 x) {
     x = s / x;
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(32)
    #elif defined __linux
        CODE_ALIGN_LINUX(32)
    #endif
#endif

static inline gms::math::AVX512c4f32
gms::math::operator~(AVX512c4f32 x) {
     x.m_re = _mm512_sub_ps(_mm512_setzero_ps(),x.m_re);
     return (x);
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator==(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_EQ_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_EQ_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator==(const AVX512c4f32 x,
                      const std::complex<float> c) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_EQ_OQ);
     const __mmask16 m2 =
           _mm512-cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_EQ_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator==(const std::complex<float> c,
                      const AVX512c4f32 x) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_EQ_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_EQ_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator!=(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_NEQ_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_NEQ_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator!=(const AVX512c4f32 x,
                      const std::complex<float> c) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_NEQ_OQ);
     const __mmask16 m2 =
           _mm512-cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_NEQ_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator!=(const std::complex<float> c,
                      const AVX512c4f32 x) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_NEQ_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_NEQ_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator>(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_GT_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_GT_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator>(const AVX512c4f32 x,
                      const std::complex<float> c) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_GT_OQ);
     const __mmask16 m2 =
           _mm512-cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_GT_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator>(const std::complex<float> c,
                      const AVX512c4f32 x) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_GT_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_GT_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator<(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_LT_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_LT_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator<(const AVX512c4f32 x,
                      const std::complex<float> c) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_LT_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_LT_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator<(const std::complex<float> c,
                      const AVX512c4f32 x) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_LT_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_LT_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator>=(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_GE_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_GE_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator>=(const AVX512c4f32 x,
                      const std::complex<float> c) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_GE_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_GE_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator>=(const std::complex<float> c,
                      const AVX512c4f32 x) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_GE_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_GE_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator<=(const AVX512c4f32 x,
                      const AVX512c4f32 y) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,y.m_re,_CMP_LE_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(x.m_im,y.m_im,_CMP_LE_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator<=(const AVX512c4f32 x,
                      const std::complex<float> c) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(x.m_re,_mm512_set1_ps(c.real()),_CMP_LE_OQ);
     const __mmask16 m2 =
           _mm512-cmp_ps_mask(x.m_im,_mm512_set1_ps(c.imag()),_CMP_LE_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif

static inline MMASK16_RETTYPE
gms::math::operator<=(const std::complex<float> c,
                      const AVX512c4f32 x) {
     const __mmask16 m1 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.real()),x.m_re,_CMP_LE_OQ);
     const __mmask16 m2 =
           _mm512_cmp_ps_mask(_mm512_set1_ps(c.imag()),x.m_im,_CMP_LE_OQ);
     return (std::make_pair(m1,m2));
}
#if defined USE_CODE_ALIGNMENT && (USE_CODE_ALIGNMENT) == 1
    #if defined _WIN64
        CODE_ALIGN_WIN(16)
    #elif defined __linux
        CODE_ALIGN_LINUX(16)
    #endif
#endif











