
#include "GMS_mkl_rayleighrng.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"

//
//	Implementation
//

gms::math::
MKLRayleighRNG::
MKLRayleighRNG() {
	
	datum.m_rvec = NULL;
	datum.m_a    = dinf;
	datum.m_beta = dinf;
	datum.m_nvalues = 0;
	datum.m_brng    = 0;
	datum.m_seed    = 0;
	datum.m_error   = 1;
}


gms::math::
MKLRayleighRNG::
MKLRayleighRNG( _In_ const MKL_INT nvalues,
			    _In_ const MKL_UINT brng,
				_In_ const MKL_INT seed,
				_In_ const double a,
				_In_ const double beta) {
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
	datum.m_rvec    = gms_mm_malloc(static_cast<size_t>(nvalues), align64B);
	datum.m_a       = a;
	datum.m_beta    = beta;
	datum.m_nvalues = nvalues;
	datum.m_brng    = brng;
	datum.m_seed    = seed;
	datum.m_error   = 1;
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&datum.m_rvec[0], static_cast<size_t>(datum.nvalues),0.0);
#elif defined __AVX__
	avx256_init_unroll8x_pd(&datum.m_rvec[0], static_cast<size_t>(datum.nvalues),0.0);
#else
#error Unsupported SIMD ISA
#endif
}


gms::math::
MKLRayleighRNG::
MKLRayleighRNG(_In_ const MKLRayleighRNG &x) {
	
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
	datum.m_rvec = gms_mm_malloc(static_cast<size_t>(x.datum.m_nvalues), align64B);

	datum.m_a       = x.datum.m_a;
	datum.m_beta    = x.datum.m_beta;
	datum.m_nvalues = x.datum.m_nvalues;
	datum.m_brng    = x.datum.m_brng;
	datum.m_seed    = x.datum.m_seed;
	datum.m_error   = x.datum.m_error;
#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	    avx512_uncached_memmove(&datum.m_rvec[0], &x.datum.m_rvec[0], static_cast<size_t>(datum.nvalues));
     #else
	    avx512_cached_memmove(&datum.m_rvec[0], &x.datum.m_rvec[0], static_cast<size_t>(datum.nvalues));
     #endif
#elif defined __AVX__
     #if (USE_NT_STORES) == 1
	    avx256_uncached_memmove(&datum.m_rvec[0], &x.datum.m_rvec[0], static_cast<size_t>(datum.nvalues));
     #else
	    avx256_cached_memmove(&datum.m_rvec[0], &x.datum.m_rvec[0], static_cast<size_t>(datum.nvalues));
     #endif
#endif
}


gms::math::
MKLRayleighRNG::
MKLRayleighRNG(_In_ MKLRayleighRNG &&x) {
	datum.m_rvec      = &x.datum.m_rvec[0];
	datum.m_a         = x.datum.m_a;
	datum.m_beta      = x.datum.m_beta;
	datum.m_nvalues   = x.datum.m_nvalues;
	datum.m_brng      = x.datum.m_brng;
	datum.m_seed      = x.datum.m_seed;
	datum.m_error     = x.datum.m_error;
	x.datum.m_rvec    =  NULL;
	x.datum.m_nvalues = 0;
}


gms::math::
MKLRayleighRNG::
~MKLRayleighRNG() {

        using namespace gms::common;
	if (NULL != datum.m_rvec) gms_mm_free(datum.m_rvec); datum.m_rvec = NULL;

}		
		
	


gms::math::MKLRayleighRNG &
gms::math::MKLRayleighRNG::
operator=(_In_ const MKLRayleighRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
        constexpr size_t align64B = 64ULL;
	gms_mm_free(datum.m_rvec);
	datum.m_a       = x.datum.m_a;
	datum.m_beta    = x.datum.m_beta;
	datum.m_nvalues = x.datum.m_nvalues;
	datum.m_brng    = x.datum.m_brng;
	datum.m_seed    = x.datum.m_seed;
	datum.m_error   = x.datum.m_error;
	double * __restrict rvec = NULL;
	rvec =  gms_mm_malloc(static_cast<size_t>(datum.m_nvalues), align64B) };
#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	    avx512_uncached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx512_cached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#elif defined __AVX__
     #if (USE_NT_STORES) == 1
	    avx256_uncached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #else
	    avx256_cached_memmove(&rvec[0], &x.m_rvec[0], static_cast<size_t>(m_nvalues));
     #endif
#endif
	datum.m_rvec = &rvec[0];
	return (*this);
}	
	


gms::math::MKLRayleighRNG &
gms::math::MKLRayleighRNG::
operator=(_In_ MKLRayleighRNG &&x) {
	if (this == &x) return (*this);
        using namespace gms::common;
	gms_mm_free(datum.m_rvec);

	datum.m_rvec = &x.datum.m_rvec[0];
	datum.m_a = x.datum.m_a;
	datum.m_beta = x.datum.m_beta;
	datum.m_nvalues = x.datum.m_nvalues;
	datum.m_brng = x.datum.m_brng;
	datum.m_seed = x.datum.m_seed;
	datum.m_error = x.datum.m_error;
	
	x.datum.m_nvalues = 0;
	x.datum.m_rvec = NULL;
	return (*this);
}	
	


void
gms::math::
MKLRayleighRNG::
compute_rand_distribution(_In_ const MKL_INT method) {
	VSLStreamStatePtr stream;
	datum.m_error = vslNewStream(&stream,datum.m_brng,datum.m_seed);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	datum.m_error = vdRngRayleigh(method, stream, datum.m_nvalues, &datum.m_rvec[0],datum.m_a,datum.m_beta);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vdRngRayleigh -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	datum.m_error = vslDeleteStream(&stream);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vslDeleteStream -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

void
gms::math::
MKLRayleighRNG::
compute_rand_distribution(_Inout_ VSLStreamStatePtr stream,
						  _In_ const MKL_INT method) {
	// Stream must be deallocated (deinitialized) by the caller of this procedure.
	datum.m_error = vslNewStream(&stream, datum.m_brng, datum.m_seed);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	datum.m_error = vdRngRayleigh(method, stream, datum.m_nvalues, &datum.m_rvec[0],datum.m_a,datum.m_beta);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vdRngRayleigh -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

std::ostream &
gms::math::
operator<<(_Inout_ std::ostream &os,
		   _In_ const MKLRayleighRNG &x) {
	for (MKL_INT i = 0; i != x.datum.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Random Rayleigh distribution: m_rvec = " << x.datum.m_rvec[i] << std::endl;
	}
	return (os);
}
