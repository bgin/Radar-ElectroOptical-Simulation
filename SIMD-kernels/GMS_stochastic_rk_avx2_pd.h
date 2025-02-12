



#ifndef __GMS_STOCHASTIC_RK_AVX2_PD_H__
#define __GMS_STOCHASTIC_RK_AVX2_PD_H__ 071120211530


namespace file_info {

      const unsigned int gGMS_STOCHASTIC_RK_AVX2_PD_MAJOR = 1;
      const unsigned int gGMS_STOCHASTIC_RK_AVX2_PD_MINOR = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX2_PD_MICRO = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX2_PD_FULLVER =
        1000*gGMS_STOCHASTIC_RK_AVX2_PD_MAJOR+100*gGMS_STOCHASTIC_RK_AVX2_PD_MINOR+
	10*gGMS_STOCHASTIC_RK_AVX2_PD_MICRO;
      const char * const pgGMS_STOCHASTIC_RK_AVX2_PD_BUILD_DATE = __DATE__":"__TIME__;
      const char * const pgGMS_STOCHASTIC_RK_AVX2_PD_CREATION_DATE = "07-11-2021 15:30  +00200 (SUN 07 NOV 2021 GMT+2)";
      const char * const pgGMS_STOCHASTIC_RK_AVX2_PD_DESCRIPTION   = "Stochastic Runge-Kutte AVX512 vectorized."

}

/*
     Modified:

    07 July 2010

  Author:

    John Burkardt

  Modified:  
    
     Bernard Gingold on 05-11-2021 14:16  +00200 (FRI 05 NOV 2021 GMT+2)
     Original implementation manually vectorized by means of AVX2 intrinsics.
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
#include <cstdint>
#include "GMS_config.h"

namespace  gms {

           namespace math {

                                    /*
                                           Calling svrng_generate4_double function!!
                                           normal  = svrng_new_normal_distribution(0.0,1.0);
					   const double * __restrict ptr = (const double*)&svrng_generate4_double(engine,normal);
					   vrand1 = _mm256_loadu_pd(&ptr[0]);

    The Runge-Kutta scheme is first-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)

                                           Parameters:

    Input, __m256d X, the values at the current time.

    Input, __m256d T, the current time (8 values).

    Input, __m256d H, the time step (8 values).

    Input, __m256d Q, the spectral density of the input white noise (8 values).

    Input, __m256d *FI ( __m256d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m256d  *GI ( __m256d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m256d STEP, the 4 values at time T+H.
                                      */
	                           
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                           
				    __m256d rk1_ti_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran, // result of call to svrng_generate4_double(engine,normal)
							       __m256d (*fi) (const __m256d),
							       __m256d (*gi) (const __m256d));

		   /*

    The Runge-Kutta scheme is second-order, and suitable for time-invariant
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                           Parameters:

    Input, __m256d X, the values at the current time.

    Input, __m256d T, the current time (4 values).

    Input, __m256d H, the time step (4 values).

    Input, __m256d Q, the spectral density of the input white noise (8 values).

    Input, __m256d *FI ( __m256d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m256d  *GI ( __m256d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m256d STEP, 4 values at time T+H.
                        */

			         
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                           
				    __m256d rk2_ti_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran1, // result of call to svrng_generate4_double(engine,normal)
							       const __m256d vran2, // result of call to svrng_generate4_double(engine,normal)
							       __m256d (*fi) (const __m256d),
							       __m256d (*gi) (const __m256d));

/*
    The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                            Parameters:

    Input, __m256d X, the values at the current time.

    Input, __m256d T, the current time (4 values).

    Input, __m256d H, the time step (4 values).

    Input, __m256d Q, the spectral density of the input white noise (4 values).

    Input, __m256d *FI ( __m256d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m256d  *GI ( __m256d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m256d STEP, the 4 values at time T+H.
*/

				  
                                   __ATTR_HOT__
                                   __ATTR_ALIGN__(32)
			           __ATTR_REGCALL__
	                         
				   __m256d rk3_ti_step_ymm4r8(const __m256d x,
				                              const __m256d t,
							      const __m256d h,
							      const __m256d q,
							      const __m256d vran1,
							      const __m256d vran2,
							      const __m256d vran3,
							      __m256d (*fi) (const __m256d),
							      __m256d (*gi) (const __m256d)); 

/*
       The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                            Parameters:

    Input, __m256d X, the values at the current time.

    Input, __m256d T, the current time (4 values).

    Input, __m256d H, the time step (4 values).

    Input, __m256d Q, the spectral density of the input white noise (4 values).

    Input, __m256d *FI ( __m256d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m256d  *GI ( __m256d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m256d STEP, the 4 values at time T+H.                       
*/

				 
                                   __ATTR_HOT__
                                   __ATTR_ALIGN__(32)
			           __ATTR_REGCALL__
	                         
				   __m256d rk4_ti_step_ymm4r8(const __m256d x,
				                              const __m256d t,
							      const __m256d h,
							      const __m256d q,
							      const __m256d vran1,
							      const __m256d vran2,
							      const __m256d vran3,
							      const __m256d vran4,
							      __m256d (*fi) (const __m256d),
							      __m256d (*gi) (const __m256d)); 
/*
            The Runge-Kutta scheme is fourth-order, and suitable for time-varying
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
*/


                                 
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                         
				    __m256d rk1_tv_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran,
							       __m256d (*fi) (const __m256d, const __m256d),
							       __m256d (*gi) (const __m256d, const __m256d) ); 


				 
                                 
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                           
				    __m256d rk2_tv_step_ymm4r8(const __m256d x,
				                               const __m2556d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran1,
							       const __m256d vran2,
							       __m256d (*fi) (const __m256d, const __m256d),
							       __m256d (*gi) (const __m256d, const __m256d));

	                       
                                  
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                           
				    __m256d rk4_tv_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran1,
							       const __m256d vran2,
							       const __m256d vran3,
							       const __m256d vran4,
							       __m256d (*fv) (const __m256d, const __m256d),
							       __m256d (*gv) (const __m256d, const __m256d));

				    
    }

}




#endif /*__GMS_STOCHASTIC_RK_AVX2_PD_H__*/
