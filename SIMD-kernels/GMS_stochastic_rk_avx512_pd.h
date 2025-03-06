

#ifndef __GMS_STOCHASTIC_RK_AVX512_PD_H__
#define __GMS_STOCHASTIC_RK_AVX512_PD_H__ 051120211416


namespace file_info {

      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PD_MAJOR = 1;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PD_MINOR = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PD_MICRO = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PD_FULLVER =
        1000*gGMS_STOCHASTIC_RK_AVX512_PD_MAJOR+100*gGMS_STOCHASTIC_RK_AVX512_PD_MINOR+
	10*gGMS_STOCHASTIC_RK_AVX512_PD_MICRO;
      const char * const pgGMS_STOCHASTIC_RK_AVX512_PD_BUILD_DATE = __DATE__":"__TIME__;
      const char * const pgGMS_STOCHASTIC_RK_AVX512_PD_CREATION_DATE = "05-11-2021 14:16  +00200 (FRI 05 NOV 2021 GMT+2)";
      const char * const pgGMS_STOCHASTIC_RK_AVX512_PD_DESCRIPTION   = "Stochastic Runge-Kutte AVX512 vectorized."

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
#include <cstdint>
#include "GMS_config.h"

namespace  gms {

           namespace math {

                                    /*
                                           Calling svrng_generate8_double function!!
                                           normal  = svrng_new_normal_distribution(0.0,1.0);
					   const double * __restrict ptr = (const double*)&svrng_generate8_double(engine,normal);
					   vrand1 = _mm512_loadu_pd(&ptr[0]);

    The Runge-Kutta scheme is first-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)

                                           Parameters:

    Input, __m512d X, the values at the current time.

    Input, __m512d T, the current time (8 values).

    Input, __m512d H, the time step (8 values).

    Input, __m512d Q, the spectral density of the input white noise (8 values).

    Input, __m512d *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512d  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512d RK1_TI_STEP, the value at time T+H.
                                      */
	                          
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                          
				    __m512d rk1_ti_step_zmm8r8(const __m512d x,
				                               const __m512d t,
							       const __m512d h,
							       const __m512d q,
							       const __m512d vran, // result of call to svrng_generate8_double(engine,normal)
							       __m512d (*fi) (const __m512d),
							       __m512d (*gi) (const __m512d));

		   /*

    The Runge-Kutta scheme is second-order, and suitable for time-invariant
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                           Parameters:

    Input, __m512d X, the values at the current time.

    Input, __m512d T, the current time (8 values).

    Input, __m512d H, the time step (8 values).

    Input, __m512d Q, the spectral density of the input white noise (8 values).

    Input, __m512d *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512d  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512d STEP, the value at time T+H.
                        */

			       
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                          
				    __m512d rk2_ti_step_zmm8r8(const __m512d x,
				                               const __m512d t,
							       const __m512d h,
							       const __m512d q,
							       const __m512d vran1,
							       const __m512d vran2,
							       __m512d (*fi) (const __m512d),
							       __m512d (*gi) (const __m512d));

/*
    The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                            Parameters:

    Input, __m512d X, the values at the current time.

    Input, __m512d T, the current time (8 values).

    Input, __m512d H, the time step (8 values).

    Input, __m512d Q, the spectral density of the input white noise (8 values).

    Input, __m512d *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512d  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512d STEP, the value at time T+H.
*/

				  
                                   __ATTR_HOT__
                                   __ATTR_ALIGN__(32)
			           __ATTR_REGCALL__
	                          
				   __m512d rk3_ti_step_zmm8r8(const __m512d x,
				                              const __m512d t,
							      const __m512d h,
							      const __m512d q,
							      const __m512d vran1,
							      const __m512d vran2,
							      const __m512d vran3,
							      __m512d (*fi) (const __m512d),
							      __m512d (*gi) (const __m512d));
/*
       The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                            Parameters:

    Input, __m512d X, the values at the current time.

    Input, __m512d T, the current time (8 values).

    Input, __m512d H, the time step (8 values).

    Input, __m512d Q, the spectral density of the input white noise (8 values).

    Input, __m512d *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512d  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512d STEP, the values at time T+H.                       
*/

				  
                                   __ATTR_HOT__
                                   __ATTR_ALIGN__(32)
			           __ATTR_REGCALL__
	                          
				   __m512d rk4_ti_step_zmm8r8(const __m512d x,
				                              const __m512d t,
							      const __m512d h,
							      const __m512d q,
							      const __m512d vran1,
							      const __m512d vran2,
							      const __m512d vran3,
							      const __m512d vran4,
							      __m512d (*fi) (const __m512d),
							      __m512d (*gi) (const __m512d));

/*
            The Runge-Kutta scheme is fourth-order, and suitable for time-varying
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
*/


                                   
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                           
				    __m512d rk1_tv_step_zmm8r8(const __m512d x,
				                               const __m512d t,
							       const __m512d h,
							       const __m512d q,
							       const __m512d vran,
							       __m512d (*fi) (const __m512d, const __m512d),
							       __m512d (*gi) (const __m512d, const __m512d) ); 


				 
                                  
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                           
				    __m512d rk2_tv_step_zmm8r8(const __m512d x,
				                               const __m512d t,
							       const __m512d h,
							       const __m512d q,
							       const __m512d vran1,
							       const __m512d vran2,
							       __m512d (*fi) (const __m512d, const __m512d),
							       __m512d (*gi) (const __m512d, const __m512d));
							       

	                       
                                  
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            
				    __m512d rk4_tv_step_zmm8r8(const __m512d x,
				                               const __m512d t,
							       const __m512d h,
							       const __m512d q,
							       const __m512d vran1,
							       const __m512d vran2,
							       const __m512d vran3,
							       const __m512d vran4,
							       __m512d (*fv) (const __m512d, const __m512d),
							       __m512d (*gv) (const __m512d, const __m512d));

				    
    }

}




#endif /*__GMS_STOCHASTIC_RK_AVX512_PD_H__*/
