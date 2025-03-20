
#include "GMS_mkl_lognormalrng.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h"  

      

//
//	Implementation
//

gms::math::
MKLLognormalRNG::
MKLLognormalRNG() {
	
	data.m_rvec  = NULL;
	data.m_a     = 0.0;
	data.m_sigma = 0.0;
	data.m_b     = 0.0;
	data.m_beta  = 0.0;
	data.m_nvalues = 0;
	data.m_brng    = 0;
	data.m_seed    = 0;
	data.m_error   = 1; 
}


gms::math::
MKLLognormalRNG::
MKLLognormalRNG(const MKL_INT nvalues,
	        const MKL_UINT brng,
	        const MKL_INT seed,
	        const double a,
	        const double sigma,
	        const double b,
	        const double beta) {
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
	data.m_rvec    = gms_mm_malloc(static_cast<size_t>(nvalues), align64B);
  

	data.m_a       = a;
	data.m_sigma   = sigma;
	data.m_b       = b;
	data.m_beta    = beta;
	data.m_nvalues = nvalues;
	data.m_brng    = brng;
	data.m_seed    = seed;
	data.m_error   = 1;
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&data.m_rvec[0], data.m_nvalues,0.0);
#elif defined __AVX__ 
	avx256_init_unroll8x_pd(&data.m_rvec[0], data.m_nvalues,0.0);
#else
#error Unsupported ISA SIMD
#endif

}


gms::math::
MKLLognormalRNG::
MKLLognormalRNG(const MKLLognormalRNG &x) {
	using namespace gms::common;
        constexpr size_t align64B = 64ULL;
        data.m_rvec    = gms_mm_malloc(static_cast<size_t>(x.data.m_nvalues), align64B);
  
	data.m_a       = x.data.m_a;
	data.m_sigma   = x.data.m_sigma;
	data.m_b       = x.data.m_b;
	data.m_beta    = x.data.m_beta;
	data.m_nvalues = x.data.m_nvalues;
	data.m_brng    = x.data.m_brng;
	data.m_seed    = x.data.m_seed;
	data.m_error   = x.data.m_error;
#if defined __AVX512F__
       #if (USE_NT_STORES) == 1
	     avx512_uncached_memove(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #else
	     avx512_cached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #endif
#else
       #if (USE_NT_STORES) == 1
	     avx256_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #else
	     avx256_cached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #endif
#endif
}


gms::math::
MKLLognormalRNG::
MKLLognormalRNG(MKLLognormalRNG &&x) {
	data.m_rvec = &x.data.m_rvec[0];
	data.m_a = x.data.m_a;
	data.m_sigma = x.data.m_sigma;
	data.m_b = x.data.m_b;
	data.m_beta = x.data.m_beta;
	data.m_nvalues = x.data.m_nvalues;
	data.m_brng = x.data.m_brng;
	data.m_seed = x.data.m_seed;
	data.m_error = x.data.m_error;
	x.data.m_rvec = NULL;
	x.data.m_nvalues = 0;
}


gms::math::
MKLLognormalRNG::
~MKLLognormalRNG() {
        using namespace gms::common;
        if (NULL != data.m_rvec) gms_mm_free(data.m_rvec); 	data.m_rvec = NULL;

}		
	
	


gms::math::MKLLognormalRNG &
gms::math::MKLLognormalRNG::
operator=(const MKLLognormalRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	if (data.m_nvalues != x.data.m_nvalues) { // Handle size mismatch
                constexpr size_t align64B = 64ULL;  
                gms_mm_free(data.m_rvec);

           
    // Preserve an invariant
		data.m_a = 0.0; data.m_sigma = 0.0;
		data.m_b = 0.0; data.m_beta = 0.0;
		data.m_nvalues = 0; data.m_brng = 0;
		data.m_seed = 0; data.m_error = 1;
                data.m_rvec = gms_mm_malloc(static_cast<size_t>(x.data.m_nvalues), align64B);
   
	}
	else {
		// Copy state
		data.m_a = x.data.m_a; data.m_sigma = x.data.m_sigma;
		data.m_b = x.data.m_b; data.m_beta  = x.data.m_beta;
		data.m_nvalues = x.data.m_nvalues; data.m_brng  = x.data.m_brng;
		data.m_seed    = x.data.m_seed;    data.m_error = x.data.m_error;
#if defined __AVX512F__
        #if (USE_NT_STORES) == 1
		avx512_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues);
        #else
		avx512_cached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues); 
        #endif
#else
        #if (USE_NT_STORES) == 1
		avx256_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues);
        #else
		avx256_cached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues); 
        #endif
#endif
	}
	return (*this);
}	
	
	

	
	
	
	


gms::math::MKLLognormalRNG &
gms::math::MKLLognormalRNG::
operator=(MKLLognormalRNG &&x) {
	if (this == &x) return (*this);
         using namespace gms::common;
        gms_mm_free(data.m_rvec);

	data.m_rvec = &x.data.m_rvec[0];
	data.m_a = x.data.m_a;
	data.m_sigma = x.data.m_sigma;
	data.m_b = x.data.m_b;
	data.m_beta = x.data.m_beta;
	data.m_nvalues = x.data.m_nvalues;
	data.m_brng = x.data.m_brng;
	data.m_seed = x.data.m_seed;
	data.m_error = x.data.m_error;
	
	x.data.m_nvalues = 0;
	x.data.m_rvec = NULL;
	return (*this);
}
	


void
gms::math::
MKLLognormalRNG::
compute_rand_distribution(const MKL_INT method) {
	VSLStreamStatePtr stream;
	data.m_error = vslNewStream(&stream,data.m_brng,data.m_seed);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	data.m_error = vdRngLognormal(method, stream, data.m_nvalues, &data.m_rvec[0],data.m_a,data.m_sigma,data.m_b,data.m_beta);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vdRngLognormal -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	data.m_error = vslDeleteStream(&stream);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vslDeleteStream -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

void
gms::math::
MKLLognormalRNG::
compute_rand_distribution(VSLStreamStatePtr stream, 
		          const MKL_INT method) {
	// Stream must be deallocated (deinitialized) by the caller of this procedure.
	data.m_error = vslNewStream(&stream,data.m_brng,data.m_seed); 
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	data.m_error = vdRngLognormal(method, stream, data.m_nvalues, &data.m_rvec[0], data.m_a, data.m_sigma, data.m_b, data.m_beta);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vdRngLognormal -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

std::ostream &
gms::math::operator<<(std::ostream &os,
			    const MKLLognormalRNG &x) {
	for (MKL_INT i = 0; i != x.data.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Random Lognormal distribution: m_rvec = " << x.data.m_rvec[i] << std::endl;
	}
	return (os);
}
