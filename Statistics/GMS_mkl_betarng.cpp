
#include "GMS_mkl_betarng.h"
#include "GMS_common.h"
#include "GMS_malloc.h"
#include "GMS_error_macros.h"
#include "GMS_constants"

//
//	Implementation
//






gms::math::stat::
MKLBetaRNG::
MKLBetaRNG() {
    using namespace gms::math::constants;
	datum.m_rvec    = NULL;
	datum.m_p       = dinf;
	datum.m_q       = dinf;
	datum.m_a       = dinf;
	datum.m_beta    = dinf;
	datum.m_nvalues = 0;
	datum.m_brng    = 0;
	datum.m_seed    = 0;
	datum.m_error   = 1;
}

gms::math::stat::
MKLBetaRNG::
MKLBetaRNG( const MKL_INT nvalues,
	    const MKL_UINT brng,
            const MKL_INT seed,
	    const double p,
	    const double q,
	    const double a,
	    const double beta) {
	using namespace gms::common;
 	datum.m_rvec = (double*)gms_mm_malloc(static_cast<size_t>(nvalues),align64B);
 
	datum.m_p    = p;
	datum.m_q    = q;
	datum.m_a    = a;
	datum.m_beta = beta;
	datum.m_nvalues = nvalues;
	datum.m_brng    = brng;
	datum.m_seed    = seed;
	datum.m_error   = 1;
#if (GMS_INIT_ARRAYS) == 1	
	avx256_init_unroll8x_pd(&datum.m_rvec[0], datum.m_nvalues, 0.0);
#endif
}

gms::math::stat::
MKLBetaRNG
::MKLBetaRNG(const MKLBetaRNG &x) {
	using namespace gms::common;

        datum.m_rvec    = (double*)gms_mm_malloc(static_cast<size_t>(x.datum.m_nvalues), align64B);

	datum.m_p       = x.datum.m_p;
	datum.m_q       = x.datum.m_q;
	datum.m_a       = x.datum.m_a;
	datum.m_beta    = x.datum.m_beta;
	datum.m_nvalues = x.datum.m_nvalues;
	datum.m_brng    = x.datum.m_brng;
	datum.m_seed    = x.datum.m_seed;
	datum.m_error   = x.datum.m_error;
#if (USE_NT_STORES) == 1
	avx256_memcpy8x_nt_pd(&datum.m_rvec[0], &datum.x.m_rvec[0], datum.m_nvalues);
#else
	avx256_memcpy8x_pd(&datum.m_rvec[0], &x.datum.m_rvec[0], datum.m_nvalues);
#endif
}


gms::math::stat::
MKLBetaRNG
::MKLBetaRNG(MKLBetaRNG &&x) {
	datum.m_rvec      = &x.datum.m_rvec[0];
	datum.m_p         = x.datum.m_p;
	datum.m_q         = x.datum.m_q;
	datum.m_a         = x.datum.m_a;
	datum.m_beta      = x.datum.m_beta;
	datum.m_nvalues   = x.datum.m_nvalues;
	datum.m_brng      = x.datum.m_brng;
	datum.m_seed      = x.datum.m_seed;
	datum.m_error     = x.datum.m_error;
	x.datum.m_rvec    = NULL;
	x.datum.m_nvalues = 0;
}


gms::math::stat::
MKLBetaRNG
::~MKLBetaRNG() {

   if (NULL != datum.m_rvec) gms_mm_free(datum.m_rvec); datum.m_rvec = NULL;

}		
		
	


gms::math::stat::MKLBetaRNG &
gms::math::stat::MKLBetaRNG::
operator=(const MKLBetaRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	if (datum.m_nvalues != x.datum.m_nvalues){

	   gms_mm_free(datum.m_rvec);

	datum.m_p = 0.0;
	datum.m_q = 0.0;
	datum.m_a = 0.0;
	datum.m_beta = 0.0;
	datum.m_nvalues = 0;
	datum.m_brng = 0;
	datum.m_seed = 0;
	datum.m_error = 1;
  	 datum.m_rvec = (double*)gms_mm_malloc(static_cast<size_t>(x.datum.m_nvalues),align64B);
   
   }
	else {
	datum.m_p = x.datum.m_a;
	datum.m_q = x.datum.m_q;
	datum.m_a = x.datum.m_a;
	datum.m_beta = x.datum.m_beta;
	datum.m_nvalues = x.datum.m_nvalues;
	datum.m_brng = x.datum.m_brng;
	datum.m_seed = x.datum.m_seed;
	datum.m_error = x.datum.m_error;
#if (USE_NT_STORES) == 1
	avx256_uncached_memmove(&datum.m_rvec[0], &x.datum.m_rvec[0], x.datum.m_nvalues);
#else
	avx256_cached_memmove(&datum.m_rvec[0], &x.datum.m_rvec[0],x.datum.m_nvalues);
#endif
	}
	
	return (*this);
}	
	


gms::math::stat::MKLBetaRNG &
gms::math::stat::MKLBetaRNG::
operator=(MKLBetaRNG &&x) {
	if (this == &x) return (*this);

        gms_mm_free(datum.m_rvec);
	datum.m_rvec      = &x.datum.m_rvec[0];
	datum.m_p         =  x.datum.m_p;
	datum.m_q         =  x.datum.m_q;
	datum.m_a         =  x.datum.m_a;
	datum.m_beta      =  x.datum.m_beta;
	datum.m_nvalues   = x.datum.m_nvalues;
	datum.m_brng      = x.datum.m_brng;
	datum.m_seed      = x.datum.m_seed;
	datum.m_error     = x.datum.m_error;
	x.datum.m_nvalues = 0;
	x.datum.m_rvec    = NULL;
	return (*this);
}	
	


void
gms::math::stat::MKLBetaRNG::
compute_rand_distribution(const MKL_INT method) {
	VSLStreamStatePtr stream;
	datum.m_error = vslNewStream(&stream,datum.m_brng,datum.m_seed);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	datum.m_error = vdRngBeta(method, stream, datum.m_nvalues, &datum.m_rvec[0],datum.m_p,datum.m_q,datum.m_a,datum.m_beta);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vdRngBeta -- failed with an error value: ", datum.m_error)
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
gms::math::stat::MKLBetaRNG::
compute_rand_distribution(VSLStreamStatePtr stream, 
			  const MKL_INT method) {
	// Stream must be deallocated (deinitialized) by the caller of this procedure.
	datum.m_error = vslNewStream(&stream,datum.m_brng,datum.m_seed);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	datum.m_error = vdRngBeta(method, stream, datum.m_nvalues, &datum.m_rvec[0], datum.m_p, datum.m_q, datum.m_a, datum.m_beta);
	if (VSL_ERROR_OK != datum.m_error) {
		PRINT_ERROR_VALUE("vdRngBeta -- failed with an error value: ", datum.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

std::ostream &
gms::math::stat::
operator<<(std::ostream &os,
	   const MKLBetaRNG &x) {
	for (MKL_INT i = 0; i != x.datum.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Random Beta distribution: m_rvec = " << x.datum.m_rvec[i] << std::endl;
	}
	return (os);
}
