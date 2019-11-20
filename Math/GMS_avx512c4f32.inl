
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
                       