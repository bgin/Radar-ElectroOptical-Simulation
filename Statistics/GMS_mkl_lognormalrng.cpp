
#include "GMS_mkl_lognormalrng.h"
#if defined _WIN64
    #include "../GMS_common.h"
    #if (GMS_DEBUG_ON) == 1
        #include "../GMS_debug.h"
    #else
        #include "../GMS_malloc.h"
    #endif
#include "../GMS_error_macros.h"
#include "../Math/GMS_constants.h"
#elif defined __linux
    #include "GMS_common.h"
    #if (GMS_DEBUG_ON) == 1
        #include "GMS_debug.h"
    #else
        #include "GMS_malloc.h"
    #endif
#include "GMS_error_macros.h"
#include "GMS_constants.h"

//
//	Implementation
//

gms::math::stat::
MKLLognormalRNG::
MKLLognormalRNG() {
	using namespace gms::math::constants;
	data.m_rvec  = NULL;
	data.m_a     = dinf;
	data.m_sigma = dinf;
	data.m_b     = dinf;
	data.m_beta  = dinf;
	data.m_nvalues = 0;
	data.m_brng    = 0;
	data.m_seed    = 0;
	data.m_error   = 1; 
}


gms::math::stat::
MKLLognormalRNG::
MKLLognormalRNG(const MKL_INT nvalues,
	        const MKL_UINT brng,
	        const MKL_INT seed,
	        const double a,
	        const double sigma,
	        const double b,
	        const double beta) {
	using namespace gms::common;
#if defined _WIN64
   #if (GMS_DEBUG_ON) == 1
	    data.m_rvec    = gms_edmalloca_dbg(static_cast<size_t>(nvalues),align64B,__FILE__,__LINE__);
   #else
	    data.m_rvec    = gms_edmalloca(static_cast<size_t>(nvalues), align64B);
   #endif
#elif defined __linux
   #if (GMS_DEBUG_ON) == 1
            data.m_rvec    = gms_edmalloca(static_cast<size_t>(nvalues), align64B);
   #else
            data.m_rvec    = gms_edmalloca(static_cast<size_t>(nvalues), align64B);
   #endif
#endif
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


gms::math::stat::
MKLLognormalRNG::
MKLLognormalRNG(const MKLLognormalRNG &x) {
	using namespace gms::common;
#if defined _WIN64
   #if (GMS_DEBUG_ON) == 1
	    data.m_rvec    = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nvalues),align64B,__FILE__,__LINE__);
   #else
	    data.m_rvec    = gms_edmalloca(static_cast<size_t>(x.data.m_nvalues), align64B);
   #endif
#elif defined __linux
   #if (GMS_DEBUG_ON) == 1
            data.m_rvec    = gms_edmalloca(static_cast<size_t>(x.data.m_nvalues), align64B);
   #else
            data.m_rvec    = gms_edmalloca(static_cast<size_t>(x.data.m_nvalues), align64B);
   #endif
#endif
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
	     avx512_memcpy8x_nt_pd(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #else
	     avx512_memcpy8x_pd(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #endif
#else
       #if (USE_NT_STORES) == 1
	     avx256_memcpy8x_nt_pd(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #else
	     avx256_memcpy8x_pd(&data.m_rvec[0], &x.data.m_rvec[0], data.m_nvalues);
       #endif
#endif
}


gms::math::stat::
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


gms::math::stat::
MKLLognormalRNG::
~MKLLognormalRNG() {
#if defined _WIN64
   #if (GMS_DEBUG_ON) == 1
	if (NULL != data.m_rvec) _aligned_free_dbg(data.m_rvec); data.m_rvec = NULL;
   #else
	if (NULL != data.m_rvec) _mm_free(data.m_rvec); 	data.m_rvec = NULL;
   #endif
#elif defined __linux
        if (NULL != data.m_rvec) _mm_free(data.m_rvec); 	data.m_rvec = NULL;
#endif
}		
	
	


gms::math::stat::MKLLognormalRNG &
gms::math::stat::MKLLognormalRNG::
operator=(const MKLLognormalRNG &x) {
	using namespace gms::common;
	if (this == &x) return (*this);
	if (data.m_nvalues != x.data.m_nvalues) { // Handle size mismatch
#if defined _WIN64   
    #if (GMS_DEBUG_ON) == 1
		_aligned_free_dbg(data.m_rvec);
    #else
		_mm_free(data.m_rvec);
    #endif
#elif defined __linux
                _mm_free(data.m_rvec);
#endif
           
    // Preserve an invariant
		data.m_a = 0.0; data.m_sigma = 0.0;
		data.m_b = 0.0; data.m_beta = 0.0;
		data.m_nvalues = 0; data.m_brng = 0;
		data.m_seed = 0; data.m_error = 1;
#if defined _WIN64
    #if (GMS_DEBUG_ON) == 1
		  data.m_rvec = gms_edmalloca_dbg(static_cast<size_t>(x.data.m_nvalues), align64B, __FILE__, __LINE__);
    #else
		  data.m_rvec = gms_edmalloca(static_cast<size_t>(x.data.m_nvalues), align64B);
    #endif
#elif defined __linux
    #if (GMS_DEBUG_ON) == 1
                  data.m_rvec = gms_edmalloca(static_cast<size_t>(x.data.m_nvalues), align64B);
    #else
                  data.m_rvec = gms_edmalloca(static_cast<size_t>(x.data.m_nvalues), align64B);
    #endif
#endif
	}
	else {
		// Copy state
		data.m_a = x.data.m_a; data.m_sigma = x.data.m_sigma;
		data.m_b = x.data.m_b; data.m_beta  = x.data.m_beta;
		data.m_nvalues = x.data.m_nvalues; data.m_brng  = x.data.m_brng;
		data.m_seed    = x.data.m_seed;    data.m_error = x.data.m_error;
#if defined __AVX512F__
        #if (USE_NT_STORES) == 1
		avx512_memcpy8x_nt_pd(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues);
        #else
		avx512_memcpy8x_pd(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues); 
        #endif
#else
        #if (USE_NT_STORES) == 1
		avx256_memcpy8x_nt_pd(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues);
        #else
		avx256_memcpy8x_pd(&data.m_rvec[0], &x.data.m_rvec[0], x.data.m_nvalues); 
        #endif
#endif
	}
	return (*this);
}	
	
	

	
	
	
	


gms::math::stat::MKLLognormalRNG &
gms::math::stat::MKLLognormalRNG::
operator=(MKLLognormalRNG &&x) {
	if (this == &x) return (*this);
#if defined _WIN64
    #if (GMS_DEBUG_ON) == 1
	  _aligned_free_dbg(data.m_rvec);
    #else
	_mm_free(data.m_rvec);
    #endif
#elif defined __linux
        _mm_free(data.m_rvec);
#endif
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
gms::math::stat::
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
gms::math::stat::
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
gms::math::stat::operator<<(std::ostream &os,
			    const MKLLognormalRNG &x) {
	for (MKL_INT i = 0; i != x.data.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Random Lognormal distribution: m_rvec = " << x.data.m_rvec[i] << std::endl;
	}
	return (os);
}
