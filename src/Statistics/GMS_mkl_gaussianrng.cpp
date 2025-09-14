
#include "GMS_mkl_gaussianrng.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"
//
//	Implementation
//

gms::math::
MKLGaussianRNG::MKLGaussianRNG() {
	
	m_rvec  = NULL;
	m_a     = 0.0;
	m_sigma = 0.0;
	m_nvalues = 0;
	m_brng    = 0;
	m_seed    = 0;
	m_error   = 1;
}


gms::math::
MKLGaussianRNG::
MKLGaussianRNG(	const MKL_INT nvalues,
		const MKL_UINT brng,
		const MKL_INT seed,
		const double a,
		const double sigma) {
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
	m_rvec  = gms_mm_malloc(static_cast<size_t>(nvalues), align64B);

	m_a     = a;
	m_sigma = sigma;
	m_nvalues = nvalues;
	m_brng    = brng;
	m_seed    = seed;
	m_error   = 1;
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&m_rvec[0], static_cast<int64_t>(m_nvalues), 0.0);
#else
	avx256_init_unroll8x_pd(&m_rvec[0], static_cast<int64_t>(m_nvalues), 0.0);
#endif
}


gms::math::
MKLGaussianRNG::
MKLGaussianRNG(const MKLGaussianRNG &x) {
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
	m_rvec    = gms_mm_malloc(static_cast<size_t>(x.m_nvalues), align64B);

	m_a       = x.m_a;
	m_sigma   = x.m_sigma;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;
#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	 avx512_uncached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	 avx512_cached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#else
     #if (USE_NT_STORES) == 1
	avx256_uncached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	avx256_cached_memmove(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
}


gms::math::
MKLGaussianRNG::
MKLGaussianRNG(MKLGaussianRNG &&x) {
	m_rvec = &x.m_rvec[0];
	m_a = x.m_a;
	m_sigma = x.m_sigma;
	m_nvalues = x.m_nvalues;
	m_brng = x.m_brng;
	m_seed = x.m_seed;
	m_error = x.m_error;
	x.m_rvec = NULL;
	x.m_nvalues = 0;
}


gms::math::
MKLGaussianRNG::
~MKLGaussianRNG() {
        using namespace gms::common;
	if (NULL != m_rvec) gms_mm_free(m_rvec); m_rvec = NULL;
}		
		
	


gms::math::MKLGaussianRNG &
gms::math::MKLGaussianRNG
::operator=(const MKLGaussianRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	gms_mm_free(m_rvec);
	m_a       = x.m_a; // Should be present in L1D for small tests.
	m_sigma   = x.m_sigma;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;
        constexpr size_t align64B = 64ULL;
	double * __restrict rvec{gms_edmalloca(static_cast<size_t>(m_nvalues), align64B) };

#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	    avx512_uncached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx512_cached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#else
     #if (USE_NT_STORES) == 1
	    avx256_uncached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx256_cached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
	m_rvec = &rvec[0];
	return (*this);
}	
	


gms::math::MKLGaussianRNG &
gms::math::MKLGaussianRNG
::operator=(MKLGaussianRNG &&x) {
        using namespace gms::common;
	if (this == &x) return (*this);
        gms_mm_free(m_rvec);
	m_rvec = &x.m_rvec[0];
	m_a = x.m_a;
	m_sigma = x.m_sigma;
	m_nvalues = x.m_nvalues;
	m_brng = x.m_brng;
	m_seed = x.m_seed;
	m_error = x.m_error;
	x.m_nvalues = 0;
	x.m_rvec = NULL;
	return (*this);
}	
	


void
gms::math::MKLGaussianRNG::
compute_rand_distribution(const MKL_INT method){
	VSLStreamStatePtr streamptr;
	m_error = vslNewStream(&streamptr, m_brng, m_seed);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vdRngGaussian(method, streamptr, m_nvalues, &m_rvec[0],m_a,m_sigma);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vdRngGaussian -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vslDeleteStream(&streamptr);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslDeleteStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
	    return;
	}
}

void
gms::math::MKLGaussianRNG::
compute_rand_distribution(VSLStreamStatePtr stream,
			   const MKL_INT method) {
	// VSLStreamStatePtr must deallocated by the caller of compute_rand_distribution.
	m_error = vslNewStream(&stream,m_brng,m_seed);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	m_error = vdRngGaussian(method, stream, m_nvalues, &m_rvec[0],m_a,m_sigma);
	if (VSL_ERROR_OK != m_error) {
		PRINT_ERROR_VALUE("vdRngGaussian -- failed with an error value: ", m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	// stream is deleted by the caller.
}

std::ostream &
gms::math::
operator<<(std::ostream &os,
	   const MKLGaussianRNG &x) {
	
	for (MKL_INT i = 0; i != x.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Randomly gaussian distributed: m_rvec = " << x.m_rvec[i] << std::endl;
	}
	return (os);
}		



