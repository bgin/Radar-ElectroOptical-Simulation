
#include "GMS_mkl_uniformrng.h"
#if defined _WIN64
    #include "../GMS_common.h"
    #include "../GMS_malloc.h"
    #include "../GMS_error_macros.h"
    #include "../Math/GMS_constants.h"
#elif defined __linux
    #include "GMS_common.h"
    #include "GMS_malloc.h"
    #include "GMS_error_macros.h"
    #include "GMS_constants.h"
#endif
//
//	Implementation
//

gms::math::stat::
MKLUniformRNG::MKLUniformRNG() {
	using namespace gms::math::constants;
	m_rvec    = NULL;
	m_a       = dinf;
	m_b       = dinf;
	m_nvalues = 0;
	m_brng    = 0;
	m_seed    = 0;
	m_error   = 1;
}



gms::math::stat::
MKLUniformRNG::
MKLUniformRNG(const MKL_INT nvalues,
	      const MKL_UINT brng,
	      const MKL_INT seed,
	      const double a,
              const double b) {
	using namespace gms::common;
#if defined _WIN64 
	    m_rvec = gms_edmalloca(static_cast<size_t>(nvalues), align64B);
#elif defined __linux
	    m_rvec = gms_edmalloca(static_cast<size_t>(nvalues), align64B);
#endif
	m_a = a;
	m_b = b;
	m_nvalues = nvalues;
	m_brng    = brng;
	m_seed    = seed;
	m_error   = 1;
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&m_rvec[0], static_cast<size_t>(m_nvalues),0.0);
#elif defined __AVX__
	avx256_init_unroll8x_pd(&m_rvec[0], static_cast<size_t>(m_nvalues),0.0);
#else
#error Unsupported SIMD ISA
#endif
}


gms::math::stat::
MKLUniformRNG::
MKLUniformRNG(const MKLUniformRNG &x) {
	using namespace gms::common;
#if defined _WIN64
	m_rvec    = gms_edmalloca(static_cast<size_t>(x.m_nvalues),align64B);
#elif defined __linux
	m_rvec    = gms_edmalloca(static_cast<size_t>(x.m_nvalues),align64B);
#endif
	m_a       = x.m_a;
	m_b       = x.m_b;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;
#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	    avx512_memcpy8x_nt_pd(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx512_memcpy8x_pd(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#elif defined __AVX__
     #if (USE_NT_STORES) == 1
	    avx256_memcpy8x_nt_pd(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx256_memcpy8x_pd(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
}


gms::math::stat::
MKLUniformRNG::
MKLUniformRNG(MKLUniformRNG &&x) {
	m_rvec      = &x.m_rvec[0];
	m_a         = x.m_a;
	m_b         = x.m_b;
	m_nvalues   = x.m_nvalues;
	m_brng      = x.m_brng;
	m_seed      = x.m_seed;
	m_error     = x.m_error;
	x.m_rvec    = NULL;
	x.m_nvalues = 0;
}



gms::math::stat::
MKLUniformRNG::
~MKLUniformRNG() {
	if (NULL != m_rvec) _mm_free(m_rvec); m_rvec = NULL;
}		
	
	


gms::math::stat::MKLUniformRNG &
gms::math::stat::MKLUniformRNG
::operator=(const MKLUniformRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	_mm_free(m_rvec);
	m_a       = x.m_a;
	m_b       = x.m_b;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;
#if defined _WIN64
	_Field_size_(m_nvalues)double * __restrict rvec{
		lam_edmalloca(static_cast<size_t>(m_nvalues), align64B) };
#elif defined __linux
	double * __restrict rvec{
		lam_edmalloca(static_cast<size_t>(m_nvalues), align64B) };
#endif
#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	    avx512_memcpy8x_nt_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx512_memcpy8x_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#elif defined __AVX__
     #if (USE_NT_STORES) == 1
	    avx256_memcpy8x_nt_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx256_memcpy8x_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
	m_rvec = &rvec[0];
	return (*this);
}	
	


gms::math::stat::MKLUniformRNG &
gms::math::stat::MKLUniformRNG::operator=(MKLUniformRNG &&x) {
	if (this == &x) return (*this);
	_mm_free(m_rvec);
	m_rvec = &x.m_rvec[0];
	m_a = x.m_a;
	m_b = x.m_b;
	m_nvalues = x.m_nvalues;
	m_brng = x.m_brng;
	m_seed = x.m_seed;
	m_error = x.m_error;
	
	x.m_nvalues = 0;
	x.m_rvec = NULL;
	return (*this);
}	
	


void
gms::math::stat::
MKLUniformRNG::compute_rand_distribution(const MKL_INT method) {
	VSLStreamStatePtr stream;
	m_error = vslNewStream(&stream,m_brng,m_seed);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vdRngUniform(method, stream, m_nvalues, &m_rvec[0],m_a,m_b);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vdRngUniform -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vslDeleteStream(&stream);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslDeleteStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

void
gms::math::stat::
MKLUniformRNG::
compute_rand_distribution(VSLStreamStatePtr stream, 
			  const MKL_INT method) {
	// caller deinitializes VSLStreamStatePtr
	m_error = vslNewStream(&stream,m_brng,m_seed);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vdRngUniform(method, stream, m_nvalues, &m_rvec[0],m_a,m_b);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vdRngUniform -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

std::ostream &
gms::math::stat::
operator<<(std::ostream &os,
           const MKLUniformRNG &x) {
	
	for (MKL_INT i = 0; i != x.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Randomly Uniform distributed: m_rvec = " << x.m_rvec[i] << std::endl;
	}
	return (os);
}
