
#include "GMS_mkl_gaussianrng.h"

#include "GMS_common.h"
#include "GMS_malloc.h"
#include "GMS_error_macros.h"
#include "GMS_constants.h"

//
//	Implementation
//

gms::math::stat::
MKLGaussianRNG::MKLGaussianRNG() {
	using namespace gms::math::constants;
	m_rvec  = NULL;
	m_a     = dinf;
	m_sigma = dinf;
	m_nvalues = 0;
	m_brng    = 0;
	m_seed    = 0;
	m_error   = 1;
}


gms::math::stat::
MKLGaussianRNG::
MKLGaussianRNG(	const MKL_INT nvalues,
		const MKL_UINT brng,
		const MKL_INT seed,
		const double a,
		const double sigma) {
	using namespace gms::common;

	m_rvec  = (double*)gms_mm_malloc(static_cast<size_t>(nvalues), align64B);

	m_a     = a;
	m_sigma = sigma;
	m_nvalues = nvalues;
	m_brng    = brng;
	m_seed    = seed;
	m_error   = 1;
#if (GMS_INIT_ARRAYS) == 1
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&m_rvec[0], static_cast<int64_t>(m_nvalues), 0.0);
#else
	avx256_init_unroll8x_pd(&m_rvec[0], static_cast<int64_t>(m_nvalues), 0.0);
#endif
#endif
}


gms::math::stat::
MKLGaussianRNG::
MKLGaussianRNG(const MKLGaussianRNG &x) {
	using namespace gms::common;

	m_rvec    = (double*)gms_mm_malloc(static_cast<size_t>(x.m_nvalues), align64B);

	m_a       = x.m_a;
	m_sigma   = x.m_sigma;
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
#else
     #if (USE_NT_STORES) == 1
	avx256_memcpy8x_nt_pd(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	avx256_memcpy8x_pd(&m_rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
}


gms::math::stat::
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


gms::math::stat::
MKLGaussianRNG::
~MKLGaussianRNG() {
	if (NULL != m_rvec) gms_mm_free(m_rvec); m_rvec = NULL;
}		
		
	


gms::math::stat::MKLGaussianRNG &
gms::math::stat::MKLGaussianRNG
::operator=(const MKLGaussianRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	_mm_free(m_rvec);
	m_a       = x.m_a; // Should be present in L1D for small tests.
	m_sigma   = x.m_sigma;
	m_nvalues = x.m_nvalues;
	m_brng    = x.m_brng;
	m_seed    = x.m_seed;
	m_error   = x.m_error;

	double * __restrict rvec = (double*)gms_mm_malloc(static_cast<size_t>(m_nvalues), align64B) };

#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	    avx512_memcpy8x_nt_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx512_memcpy8x_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#else
     #if (USE_NT_STORES) == 1
	    avx256_memcpy8x_nt_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx256_memcpy8x_pd(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
	m_rvec = &rvec[0];
	return (*this);
}	
	


gms::math::stat::MKLGaussianRNG &
gms::math::stat::MKLGaussianRNG
::operator=(MKLGaussianRNG &&x) {
	if (this == &x) return (*this);
	_mm_free(m_rvec);
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
gms::math::stat::MKLGaussianRNG::
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
gms::math::stat::MKLGaussianRNG::
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
gms::math::stat::
operator<<(std::ostream &os,
	   const MKLGaussianRNG &x) {
	
	for (MKL_INT i = 0; i != x.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Randomly gaussian distributed: m_rvec = " << x.m_rvec[i] << std::endl;
	}
	return (os);
}		



