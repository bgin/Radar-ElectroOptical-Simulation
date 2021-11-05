

#ifndef __GMS_STOCHASTIC_RK_AVX512_HPP__
#define __GMS_STOCHASTIC_RK_AVX512_HPP__ 051120211416


namespace file_info {

      const unsigned int gGMS_STOCHASTIC_RK_AVX512_MAJOR = 1;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_MINOR = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_MICRO = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_FULLVER =
        1000*gGMS_STOCHASTIC_RK_AVX512_MAJOR+100*gGMS_STOCHASTIC_RK_AVX512_MINOR+
	10*gGMS_STOCHASTIC_RK_AVX512_MICRO;
      const char * const pgGMS_STOCHASTIC_RK_AVX512_BUILD_DATE = __DATE__":"__TIME__;
      const char * const pgGMS_STOCHASTIC_RK_AVX512_CREATION_DATE = "05-11-2021 14:16  +00200 (FRI 05 NOV 2021 GMT+2)";
      const char * const pgGMS_STOCHASTIC_RK_AVX512_DESCRIPTION   = "Stochastic Runge-Kutte AVX512 vectorized."

}

/*
     Modified:

    07 July 2010

  Author:

    John Burkardt

  Modified:  
    
     Bernard Gingold on 05-11-2021 14:16  +00200 (FRI 05 NOV 2021 GMT+2)
     Original implementation manually vectorized by means of AVX512 intrinsics.
     Removed J. Burkardt pseudorandom (scalar) generators.

  Reference:

    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.

    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
*/


#include <immintrin.h>
#if defined(__INTEL_COMPILER) || defined(__ICC)
#include <svrng.h>
#else
    #error "Compilation of this translation unit requires an Intel Compiler and implementation of SVRNG library."
#endif
#include <cstdint>
#include "GMS_config.h"

namespace  gms {

           namespace math {

                                    /*
                                           Parameters:

    Input, __m512d X, the values at the current time.

    Input, __m512d T, the current time (8 values).

    Input, __m512d H, the time step (8 values).

    Input, __m512d Q, the spectral density of the input white noise (8 values).

    Input, double FI ( double X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, double GI ( double X ), the name of the stochastic
    right hand side function vector SIMD.

    Input/output, int *SEED, a seed for the random 
    number generator.

    Output, double RK1_TI_STEP, the value at time T+H.
                                      */
	                            __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512d rk1_ti_step_zmm8r8(const __m512d x,
				                               const __m512d t,
							       const __m512d h,
							       const __m512d q,
							       __m512d (*fi) (const __m512d),
							       __m512d (*gi) (const __m512d),
							       int32_t &seed) {
							       
                                           register const  __m512d a21 = _mm512_set1_pd(1.0);
					   register const  __m512d q1  = _mm512_set1_pd(1.0);
					   register        __m512d w1;
					   register        __m512d val;
					   register        __m512d k1;
					   register        __m512d tmp0;
					   register        __m512d tmp1;
					   register        __m512d tmp2;
					   register        __m512d tmp3;
					   register        __m512d vrand1;
					   svrng_engine_t engine;
					   svrng_distribution_t normal;   
		                           //int32_t err = -9999;                          
                                     	   int32_t result = -9999;                       
		                           result = _rdrand32_step(&seed);                
		                           if(__builtin_expect(!result,0)) seed = (int32_t)(__rdtsc());
					   tmp2 = gi(x);
		                           engine = svrng_new_mt19937_engine(seed);
					   // Currently 'error handling' code deactivated!!
		                           //err = svrng_get_status();                      
		                           //if(err != SVRNG_STATUS_OK) {                   
                                           //   status = err;                              
		                           //    return;                                    
		                           //}
					   tmp1  = fi(x);                   
                                  	   normal  = svrng_new_normal_distribution(0.0,1.0);
					   const double * __restrict ptr = (const double*)&svrng_generate8_double(engine,normal);
					   vrand1 = _mm512_loadu_pd(&ptr[0]);
					   tmp0   = _mm512_mul_pd(q1,_mm512_div_pd(q,h));
					   w1     = _mm512_mul_pd(vrand1,tmp0);
					   tmp3   = _mm512_mul_pd(h,_mm512_mul_pd(tmp2,w1);
					   k1     = _mm512_fmadd_pd(h,tmp1,tmp3);
					   val    = _mm512_fmadd_pd(a21,k1,x);
					   svrng_delete_distribution(normal);
		                           svrng_delete_engine(engine);
					   return (val);
				 }
    }

}




#endif /*__GMS_STOCHASTIC_RK_AVX512_HPP__*/
