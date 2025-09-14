
#include "GMS_mkl_weibullrng.h"
#include "GMS_common.h"
#include "GMS_malloc.h"
#include "GMS_simd_memops.h" 


//
//	Implementation
//

gms::math::
MKLWeibullRNG::
MKLWeibullRNG() {
	using namespace gms::math::constants;
	data.m_rvec    = NULL;
	data.m_alpha   = dinf;
	data.m_a       = dinf;
	data.m_beta    = dinf;
	data.m_nvalues = 0;
	data.m_brng    = 0;
	data.m_seed    = 0;
	data.m_error   = 1;
}



gms::math::
MKLWeibullRNG::
MKLWeibullRNG(const MKL_INT nvalues,
	      const MKL_UINT brng,
	      const MKL_INT seed,
	      const double alpha,
	      const double a,
	      const double beta) {
	using namespace gms::common;

        constexpr size_t align64b = 64ULL;
        data.m_rvec = gms_mm_malloc(static_cast<size_t>(nvalues),align64b ); // <-- Upon touching m_rvec load 64byte(whole object) to cache
    

	data.m_alpha   = alpha;
	data.m_a       = a;
	data.m_beta    = beta;
	data.m_nvalues = nvalues;
	data.m_brng    = brng;
	data.m_seed    = seed;
	data.m_error   = 0;
#if defined __AVX512F__
        avx512_init_unroll8x_pd(&data.m_rvec[0], static_cast<int64_t>(data.m_nvalues), 0.0);
#elif defined __AVX__
	avx256_init_unroll8x_pd(&data.m_rvec[0], static_cast<int64_t>(data.m_nvalues), 0.0); // m_rvec should be in cache
#else
#error Unsupported SIMD ISA
#endif
}
	





gms::math::
MKLWeibullRNG::
MKLWeibullRNG(const MKLWeibullRNG &x) {
	using namespace gms::common;
        
        constexpr size_t align64b = 64ULL;
        data.m_rvec = gms_mm_malloc(static_cast<size_t>(x.m_nvalues), align64b );
   
	data.m_alpha   = x.data.m_alpha;
	data.m_a       = x.data.m_a;
	data.m_beta    = x.data.m_beta;
	data.m_nvalues = x.data.m_nvalues;
	data.m_brng    = x.data.m_brng;
	data.m_seed    = x.data.m_seed;
	data.m_error   = x.data.m_error;
#if defined __AVX512F__
     #if (USE_NT_STORES) == 1
	 avx512_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(data.m_nvalues));
     #else
	 avx512_cached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(data.m_nvalues));
     #endif
#elif defined __AVX__
     #if (USE_NT_STORES) == 1
	   avx256_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(data.m_nvalues));
     #else
	   avx256_cached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(data.m_nvalues));
     #endif
#else
#error Unsupported SIMD ISA
#endif
}	
	




gms::math::
MKLWeibullRNG::
MKLWeibullRNG(MKLWeibullRNG &&x){
	
	data.m_rvec      = &x.data.m_rvec[0];
	data.m_alpha     = x.data.m_alpha;
	data.m_a         = x.data.m_a;
	data.m_beta      = x.data.m_beta;
	data.m_nvalues   = x.data.m_nvalues;
	data.m_brng      = x.data.m_brng;
	data.m_seed      = x.data.m_seed;
	data.m_error     = x.data.m_error;
	x.data.m_nvalues = 0;
	x.data.m_rvec    = NULL;
}	
	




gms::math::
MKLWeibullRNG::
~MKLWeibullRNG() {
       if (NULL != data.m_rvec) gms_mm_free(data.m_rvec); data.m_rvec = NULL;
}		
		
	


gms::math::MKLWeibullRNG &
gms::math::MKLWeibullRNG::
operator=(const MKLWeibullRNG &x){
    using namespace gms::common;
	if (this == &x) return (*this);
	if (data.m_nvalues != x.data.m_nvalues) { // Handle size mismatch
	     constexpr size_t align64b = 64ULL;
             gms_mm_free(data.m_rvec);
	     data.m_alpha = 0.0; data.m_a = 0.0;
	     data.m_beta = 0.0;  data.m_nvalues = 0;
	     data.m_brng = 0;    data.m_seed = 0;
	     data.m_error = 1;
             data.m_rvec = gms_mm_malloc(static_cast<size_t>(x.data.m_nvalues), align64b);
 
	}
	else {
		// Copy state
		data.m_alpha = x.data.m_alpha;data.m_a = x.data.m_a;
		data.m_beta  = x.data.m_beta; data.m_nvalues = x.data.m_nvalues;
		data.m_brng  = x.data.m_brng; data.m_seed   = x.data.m_seed;
		data.m_error = x.data.m_error;
#if defined __AVX512F__
      #if (USE_NT_STORES) == 1
		avx512_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(x.data.m_nvalues));
      #else
		avx512_cached_memove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(x.data.m_nvalues));
      #endif
#elif defined __AVX__
      #if (USE_NT_STORES) == 1
		avx256_uncached_memmove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(x.data.m_nvalues));
      #else
		avx256_cached_memove(&data.m_rvec[0], &x.data.m_rvec[0], static_cast<size_t>(x.data.m_nvalues));
      #endif
#else
#error Unsupported SIMD ISA
#endif
	}
	return (*this);
}	
	

	

	
	
	
	



gms::math::MKLWeibullRNG &
gms::math::MKLWeibullRNG::
operator=(MKLWeibullRNG &&x) {
	if (this == &x) return (*this);

          gms_mm_free(data.m_rvec);
	  data.m_rvec      = &x.data.m_rvec[0];
	  data.m_alpha     = x.data.m_alpha;
	  data.m_a         = x.data.m_a;
	  data.m_beta      = x.data.m_beta;
	  data.m_nvalues   = x.data.m_nvalues;
	  data.m_brng      = x.data.m_brng;
	  data.m_seed      = x.data.m_seed;
	  data.m_error     = x.data.m_error;
	  x.data.m_nvalues = 0;
	  x.data.m_rvec    = NULL;
	  return (*this);
}	
	
	


void
gms::math::MKLWeibullRNG::
compute_rand_distribution(const MKL_INT method) {
	VSLStreamStatePtr stream;
	data.m_error = vslNewStream(&stream,data.m_brng,data.m_seed);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
	data.m_error = vdRngWeibull(method, stream, data.m_nvalues, &data.m_rvec[0],data.m_alpha,data.m_a,data.m_beta);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vdRngWeibull -- failed with an error value: ", data.m_error)
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
gms::math::MKLWeibullRNG::
compute_rand_distribution(VSLStreamStatePtr stream, 
		          const MKL_INT method) {
	// Stream must be deallocated (deinitialized) by the caller of this procedure.
	data.m_error = vslNewStream(&stream, data.m_brng, data.m_seed);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vslNewStream -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
			return;
	}
	data.m_error = vdRngWeibull(method, stream, data.m_nvalues, &data.m_rvec[0], data.m_alpha, data.m_a, data.m_beta);
	if (VSL_ERROR_OK != data.m_error) {
		PRINT_ERROR_VALUE("vdRngWeibull -- failed with an error value: ", data.m_error)
		DUMP_CALLSTACK_ON_ERROR
		return;
	}
}

std::ostream &
gms::math::stat::
operator<<(std::ostream &os,
	   const MKLWeibullRNG &x) {
	for (MKL_INT i = 0; i != x.data.m_nvalues; ++i) {
		os << std::setprecision(15) << std::showpoint <<
			"Random Weibull distribution: m_rvec = " << x.data.m_rvec[i] << std::endl;
	}
	return (os);
}
