

#ifndef __GMS_STOCHASTIC_RK_AVX512_PS_HPP__
#define __GMS_STOCHASTIC_RK_AVX512_PS_HPP__ 071120211318
/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_info {

      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PS_MAJOR = 1;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PS_MINOR = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PS_MICRO = 0;
      const unsigned int gGMS_STOCHASTIC_RK_AVX512_PS_FULLVER =
        1000*gGMS_STOCHASTIC_RK_AVX512_PS_MAJOR+100*gGMS_STOCHASTIC_RK_AVX512_PS_MINOR+
	10*gGMS_STOCHASTIC_RK_AVX512_PS_MICRO;
      const char * const pgGMS_STOCHASTIC_RK_AVX512_PS_BUILD_DATE = __DATE__":"__TIME__;
      const char * const pgGMS_STOCHASTIC_RK_AVX512_PS_CREATION_DATE = "07-11-2021 13:18  +00200 (SUN 07 NOV 2021 GMT+2)";
      const char * const pgGMS_STOCHASTIC_RK_AVX512_PS_DESCRIPTION   = "Stochastic Runge-Kutte AVX512 vectorized."

}

/*
     Modified:

    07 July 2010

  Author:

    John Burkardt

  Modified:  
    
     Bernard Gingold on 07-11-2021 13:16  +00200 (SUN 07 NOV 2021 GMT+2)
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
					   const float * __restrict ptr = (const float*)&svrng_generate16_float(engine,normal);
					   vrand1 = _mm512_loadu_ps(&ptr[0]);

    The Runge-Kutta scheme is first-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)

                                           Parameters:

    Input, __m512 X, the values at the current time.

    Input, __m512 T, the current time (8 values).

    Input, __m512 H, the time step (8 values).

    Input, __m512 Q, the spectral density of the input white noise (8 values).

    Input, __m512 *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512 RK1_TI_STEP, the value at time T+H.
                                      */
	                            __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512 rk1_ti_step_zmm16r4(const __m512 x,
				                                const __m512 t,
							        const __m512 h,
							        const __m512 q,
							        const __m512 vran, // result of call to svrng_generate8_double(engine,normal)
							       __m512 (*fi) (const __m512),
							       __m512 (*gi) (const __m512)){
							       
							       
                                           register const  __m512 a21 = _mm512_set1_ps(1.0);
					   register const  __m512 q1  = _mm512_set1_ps(1.0);
					   register        __m512 w1;
					   register        __m512 step;
					   register        __m512 k1;
					   register        __m512 tmp0;
					   register        __m512 tmp1;
					   register        __m512 tmp2;
					   register        __m512 tmp3;
					   tmp2   = gi(x);
		                           tmp0   = _mm512_mul_ps(q1,_mm512_div_ps(q,h));
					   w1     = _mm512_mul_ps(vran,tmp0);
					   tmp1   = fi(x);    
					   tmp3   = _mm512_mul_ps(h,_mm512_mul_pd(tmp2,w1));
					   k1     = _mm512_fmadd_ps(h,tmp1,tmp3);
					   step   = _mm512_fmadd_ps(a21,k1,x);
					   return (step);
				 }

		   /*

    The Runge-Kutta scheme is second-order, and suitable for time-invariant
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                           Parameters:

    Input, __m512 X, the values at the current time.

    Input, __m512 T, the current time (8 values).

    Input, __m512 H, the time step (8 values).

    Input, __m512 Q, the spectral density of the input white noise (8 values).

    Input, __m512 *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512 STEP, the value at time T+H.
                        */

			            __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512 rk2_ti_step_zmm16r4(const __m512 x,
				                                const __m512 t,
							        const __m512 h,
							        const __m512 q,
							        const __m512 vran1,
							        const __m512 vran2,
							        __m512 (*fi) (const __m512),
							        __m512 (*gi) (const __m512)){
							      
                                        const __m512 _0  = _mm512_setzero_ps();
                                        const __m512 a21 = _mm512_set1_ps(1.0F);
					const __m512 a31 = _mm512_set1_ps(0.5F);
					const __m512 a32 = a31;
					const __m512 q1  = _mm512_set1_ps(2.0F);
					const __m512 q2  = q1;
					const __m512 qh  = _mm512_div_ps(q,h);
					register __m512 tfi       = _0;
					register __m512 tgi       = _0;
					register __m512 w1        = _0;
					register __m512 w2        = _0;
					register __m512 x2        = _0;
					register __m512 k1        = _0;
					register __m512 k2        = _0;
					register __m512 t0        = _0;
					register __m512 t1        = _0;
					register __m512 t2        = _0;
					register __m512 step      = _0;
				        
					tfi  = fi(x);
					t0   = _mm512_sqrt_ps(_mm512_mul_ps(q1,qh));
					tgi  = gi(x);
					w1   = _mm512_mul_ps(vran1,t0);
					t1   = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w1));
					k1   = _mm512_fmadd_ps(h,tfi,t1);
					x2   = _mm512_fmadd_ps(a21,k1,x);
					tfi  = fi(x2);
					t0   = _mm512_sqrt_ps(_mm512_mul_ps(q2,qh));
					w2   = _mm512_mul_ps(vran2,t0);
					tgi  = gi(x2);
					t2   = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w2));
					k2   = _mm512_fmadd_ps(h,tfi,t2);
					step = _mm512_fmadd_ps(a32,k2,_mm512_fmadd_ps(a31,k1,x));
					return (step);
				}

/*
    The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                            Parameters:

    Input, __m512 X, the values at the current time.

    Input, __m512 T, the current time (8 values).

    Input, __m512 H, the time step (8 values).

    Input, __m512 Q, the spectral density of the input white noise (8 values).

    Input, __m512 *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512 STEP, the value at time T+H.
*/

				   __ATTR_ALWAYS_INLINE__
                                   __ATTR_HOT__
                                   __ATTR_ALIGN__(32)
			           __ATTR_REGCALL__
	                           static inline
				   __m512 rk3_ti_step_zmm16r4(const __m512 x,
				                               const __m512 t,
							       const __m512 h,
							       const __m512 q,
							       const __m512 vran1,
							       const __m512 vran2,
							       const __m512 vran3,
							      __m512 (*fi) (const __m512),
							      __m512 (*gi) (const __m512)) {
							     
 
                                         const __m512 _0      = _mm512_setzero_ps();
					 const __m512 a21     = _mm512_set1_ps(1.52880952525675);
					 const __m512 a31     = _0;
					 const __m512 a32     = _mm512_set1_ps(0.51578733443615);
					 const __m512 a41     = _mm512_set1_ps(0.53289582961739);
					 const __m512 a42     = _mm512_set1_ps(0.25574324768195);
					 const __m512 a43     = _mm512_set1_ps(0.21136092270067);
					 const __m512 q1      = _mm512_set1_ps(1.87653936176981);
					 const __m512 q2      = _mm512_set1_ps(3.91017166264989);
					 const __m512 q3      = _mm512_set1_ps(4.73124353935667);
					 const __m512 qh      = _mm512_div_ps(q,h);
					 register __m512 k1   = _0;
					 register __m512 k2   = _0;
					 register __m512 k3   = _0;
					 register __m512 t0   = _0;
					 register __m512 t1   = _0;
					 register __m512 t2   = _0;
					 register __m512 w1   = _0;
					 register __m512 w2   = _0;
					 register __m512 w3   = _0;
					 register __m512 x2   = _0;
					 register __m512 x3   = _0;
					 register __m512 step = _0;
					 __m512          tfi  = _0;
					 __m512          tgi  = _0;
					 __m512          tmp  = _0;
					 __m512          tmp2 = _0;
					
		                      	 tfi  = fi(x);
					 t0   = _mm512_sqrt_ps(_mm512_mul_ps(q1,qh));
					 t1   = _mm512_sqrt_ps(_mm512_mul_ps(q2,qh));
					 t2   = _mm512_sqrt_ps(_mm512_mul_ps(q3,qh));
					 tgi  = gi(x);
					 w1   = _mm512_mul_ps(vran1,t0);
					 tmp  = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w1));
					 k1   = _mm512_fmadd_ps(h,tfi,tmp);
					 tmp2 = _mm512_fmadd_ps(a41,k1,x);
					 x2   = _mm512_fmadd_ps(a21,k1,x);
					 tfi  = fi(x2);
					 w2   = _mm512_mul_ps(vran2,t1);
					 tgi  = gi(x2);
					 tmp  = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w2));
					 k2   = _mm512_fmadd_ps(h,tfi,tmp);
					 x3   = _mm512_fmadd_ps(a32,k2,_mm512_fmadd_ps(a31,k1,x));
					 tfi  = fi(x3);
					 w3   = _mm512_mul_ps(vran3,t2);
					 tgi  = gi(x3);
					 tmp  = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w3));
					 k3   = _mm512_fmadd_ps(h,tfi,tmp);
					 step = _mm512_fmadd_ps(a43,k3,_mm512_fmadd_ps(a42,k2,tmp2));
					 return (step);
				   }

/*
       The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.

    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
                                            Parameters:

    Input, __m512 X, the values at the current time.

    Input, __m512 T, the current time (8 values).

    Input, __m512 H, the time step (8 values).

    Input, __m512 Q, the spectral density of the input white noise (8 values).

    Input, __m512 *FI ( __m512d X ), the name of the deterministic
    right hand side function vector SIMD.

    Input, __m512  *GI ( __m512d X ), the name of the stochastic
    right hand side function vector SIMD.

   

    Output, __m512 STEP, the values at time T+H.                       
*/

				   __ATTR_ALWAYS_INLINE__
                                   __ATTR_HOT__
                                   __ATTR_ALIGN__(32)
			           __ATTR_REGCALL__
	                           static inline
				   __m512 rk4_ti_step_zmm16r4(const __m512 x,
				                               const __m512 t,
							       const __m512 h,
							       const __m512 q,
							       const __m512 vran1,
							       const __m512 vran2,
							       const __m512 vran3,
							       const __m512 vran4,
							      __m512d (*fi) (const __m512),
							      __m512d (*gi) (const __m512)) {

                                           const __m512 _0    = _mm512_setzero_ps();
					   const __m512 a21   = _mm512_set1_ps(2.71644396264860);
					   const __m512 a31   = _mm512_set1_ps(-6.95653259006152);
					   const __m512 a32   = _mm512_set1_ps(0.78313689457981);
					   const __m512 a41   = _0;
					   const __m512 a42   = _mm512_set1_ps(0.48257353309214);
					   const __m512 a43   = _mm512_set1_ps(0.26171080165848);
					   const __m512 a51   = _mm512_set1_ps(0.47012396888046);
					   const __m512 a52   = _mm512_set1_ps(0.36597075368373);
					   const __m512 a53   = _mm512_set1_ps(0.08906615686702);
					   const __m512 a54   = _mm512_set1_ps(0.07483912056879);
					   const __m512 q1    = _mm512_set1_ps(2.12709852335625);
					   const __m512 q2    = _mm512_set1_ps(2.73245878238737);
					   const __m512 q3    = _mm512_set1_ps(11.22760917474960);
					   const __m512 q4    = _mm512_set1_ps(13.36199560336697);
					   const __m512 qh    = _mm512_div_ps(q,h);
					   register __m512 k1 = _0;
					   register __m512 k2 = _0;
					   register __m512 k3 = _0;
					   register __m512 k4 = _0;
					   register __m512 t1 = _0;
					   register __m512 t2 = _0;
					   register __m512 t3 = _0;
					   register __m512 t4 = _0;
					   register __m512 w1 = _0;
					   register __m512 w2 = _0;
					   register __m512 w3 = _0;
					   register __m512 w4 = _0;
					   register __m512 x2 = _0;
					   register __m512 x3 = _0;
					   register __m512 x4 = _0;
					   register __m512 step = _0;
					   __m512          tfi  = _0;
					   __m512          tgi  = _0;
					   __m512          tmp1 = _0;
					   __m512         tmp2 = _0;
					   tfi  = fi(x);
					   w1   = _mm512_mul_ps(vran1,_mm512_sqrt_ps(_mm512_mul_ps(q1,qh)));
					   tgi  = gi(x);
					   tmp1 = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w1));
					   k1   = _mm512_fmadd_ps(h,tfi,tmp1);
					   x2   = _mm512_fmadd_ps(a21,k1,x);
					   tfi  = fi(x2);
					   w2   = _mm512_mul_ps(vran2,_mm512_sqrt_ps(_mm512_mul_ps(q2,qh)));
					   tgi  = gi(x2);
					   tmp1 = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w2));
					   k2   = _mm512_fmadd_ps(h,tfi,tmp1);
					   x3   = _mm512_fmadd_ps(a32,k2,_mm512_fmadd_ps(a31,k1,x1));
					   tfi  = fi(x3);
					   w3   = _mm512_mul_ps(vran3,_mm512_sqrt_ps(_mm512_mul_ps(q3,qh)));
					   tgi  = gi(x3);
					   tmp1 = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w3));
					   k3   = _mm512_fmadd_ps(h,tfi,tmp1);
					   x4   = _mm512_fmadd_ps(a42,k2,_mm512_fmadd_ps(a41,k1,x1));
					   tfi  = fi(x4);
					   w4   = _mm512_mul_ps(vran4,_mm512_sqrt_ps(_mm512_mul_ps(q4,qh)));
					   tgi  = gi(x4);
					   tmp1 = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w4));
					   k4   = _mm512_fmadd_ps(h,tfi,tmp1);
					   step = _mm512_fmadd_ps(a51,k1,x1);
					   step = _mm512_fmadd_ps(a52,k2,step);
					   step = _mm512_fmadd_ps(a53,k3,step);
					   step = _mm512_fmadd_ps(a54,k4,step);
					   return (step);
				  }

/*
            The Runge-Kutta scheme is fourth-order, and suitable for time-varying
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
*/


                                    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512 rk1_tv_step_zmm16r4(const __m512 x,
				                                const __m512 t,
							        const __m512 h,
							        const __m512 q,
							        const __m512 vran,
							        __m512 (*fi) (const __m512, const __m512),
							        __m512 (*gi) (const __m512, const __m512) ) {
							        
							       
                                           register const  __m512 a21 = _mm512_set1_ps(1.0F);
					   register const  __m512 q1  = _mm512_set1_ps(1.0F);
					   register        __m512 w1;
					   register        __m512 step;
					   register        __m512 k1;
					   register        __m512 tmp0;
					   register        __m512 tmp1;
					   register        __m512 tmp2;
					   register        __m512 tmp3;
					 
		                           tmp1   = fi(t,x);                   
                                  	   tmp0   = _mm512_mul_ps(q1,_mm512_div_ps(q,h));
					   w1     = _mm512_mul_ps(vran,tmp0);
					   tmp2   = gi(t,x);
					   tmp3   = _mm512_mul_ps(h,_mm512_mul_ps(tmp2,w1));
					   k1     = _mm512_fmadd_ps(h,tmp1,tmp3);
					   step   = _mm512_fmadd_ps(a21,k1,x);
					   return (step);
				 }


				 
                                    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512 rk2_tv_step_zmm16r4(const __m512 x,
				                                const __m512 t,
							        const __m512 h,
							        const __m512 q,
							        const __m512 vran1,
							        const __m512 vran2,
							        __m512 (*fv) (const __m512, const __m512),
							        __m512 (*gv) (const __m512, const __m512)){
							      
                                        register const __m512 _0  = _mm512_setzero_ps();
                                        register const __m512 a21 = _mm512_set1_ps(1.0F);
					register const __m512 a31 = _mm512_set1_ps(0.5F);
					register const __m512 a32 = a31;
					register const __m512 q1  = _mm512_set1_ps(2.0F);
					register const __m512 q2  = q1;
					const __m512d         qh  = _mm512_div_ps(q,h);
					register __m512 tfv       = _0;
					register __m512 tgv       = _0;
					register __m512 w1        = _0;
					register __m512 w2        = _0;
					register __m512 x2        = _0;
					register __m512 k1        = _0;
					register __m512 k2        = _0;
					register __m512 t0        = _0;
					register __m512 t1        = _0;
					register __m512 t2        = _0;
					register __m512 tt        = _0;
					register __m512 step      = _0;
				        
					tfv  = fv(t,x);
					t0   = _mm512_sqrt_ps(_mm512_mul_ps(q1,qh));
					tgv  = gv(t,x);
					w1   = _mm512_mul_ps(vran1,t0);
					t1   = _mm512_mul_ps(h,_mm512_mul_pd(tgi,w1));
					k1   = _mm512_fmadd_ps(h,tfi,t1);
					x2   = _mm512_fmadd_ps(a21,k1,x);
					tt   = _mm512_fmadd_ps(a21,h,t1);
					tfv  = fv(tt,x2);
					t0   = _mm512_sqrt_ps(_mm512_mul_ps(q2,qh));
					w2   = _mm512_mul_ps(vran2,t0);
					tgv  = gv(tt,x2);
					t2   = _mm512_mul_ps(h,_mm512_mul_ps(tgi,w2));
					k2   = _mm512_fmadd_ps(h,tfi,t2);
					step = _mm512_fmadd_ps(a32,k2,_mm512_fmadd_ps(a31,k1,x));
					return (step);
				}


	                       
                                    __ATTR_ALWAYS_INLINE__
                                    __ATTR_HOT__
                                    __ATTR_ALIGN__(32)
			            __ATTR_REGCALL__
	                            static inline
				    __m512 rk4_tv_step_zmm16r4(const __m512 x,
				                                const __m512 t,
							        const __m512 h,
							        const __m512 q,
							        const __m512 vran1,
							        const __m512 vran2,
							        const __m512 vran3,
							        const __m512 vran4,
							       __m512 (*fv) (const __m512, const __m512),
							       __m512 (*gv) (const __m512, const __m512)){

					    const __m512 _0   = _mm512_setzero_ps();
					    const __m512 a21  = _mm512_set1_ps(0.66667754298442);
					    const __m512 a31  = _mm512_set1_ps(0.63493935027993);
					    const __m512 a32  = _mm512_set1_ps(0.00342761715422);
					    const __m512 a41  = _mm512_set1_ps(-2.32428921184321);
					    const __m512 a42  = _mm512_set1_ps(2.69723745129487);
					    const __m512 a43  = _mm512_set1_ps(0.29093673271592);
					    const __m512 a51  = _mm512_set1_ps(0.25001351164789);
					    const __m512 a52  = _mm512_set1_ps(0.67428574806272);
					    const __m512 a53  = _mm512_set1_ps(-0.00831795169360);
					    const __m512 a54  = _mm512_set1_ps(0.08401868181222);
					    const __m512 q1   = _mm512_set1_ps(3.99956364361748);
					    const __m512 q2   = _mm512_set1_ps(1.64524970733585);
					    const __m512 q3   = _mm512_set1_ps(1.59330355118722);
					    const __m512 q4   = _mm512_set1_ps(0.26330006501868);
					    const __m512 qh   = _mm512_div_ps(q,h);
					    register __m512 k1 = _0;
					    register __m512 k2 = _0;
					    register __m512 k3 = _0;
					    register __m512 k4 = _0;
					    register __m512 t2 = _0;
					    register __m512 t3 = _0;
					    register __m512 t4 = _0;
					    register __m512 x2 = _0;
					    register __m512 x3 = _0;
					    register __m512 x4 = _0;
					    register __m512 tt2=_0;
					    register __m512 tt3=_0;
					    register __m512 tt4=_0;
					    __m512          tgv=_0;
					    __m512          tfv=_0;
					    __m512          step=_0;
					    __m512          tmp =_0;
					    tfv = fv(t,x);
					    w1  = _mm512_mul_ps(vran1,_mm512_sqrt_ps(_mm512_mul_ps(q1,qh)));
					    tgv = gv(t,x);
					    tmp = _mm512_mul_ps(h,_mm512_mul_ps(tgv,w1));
					    k1  = _mm512_fmadd_ps(h,tfv,tmp);
					    tt2 = _mm512_fmadd_ps(a21,h,t);
					    x2  = _mm512_fmadd_ps(a21,k1,x);
					    tfv = fv(tt2,x2);
					    w2  = _mm512_mul_ps(vran2,_mm512_sqrt_ps(_mm512_mul_ps(q2,qh)));
					    tgv = gv(tt2,x2);
					    tmp = _mm512_mul_ps(h,_mm512_mul_ps(tgv,w2));
					    k2  = _mm512_fmadd_ps(h,tfv,tmp);
					    tt3 = _mm512_fmadd_ps(a32,h,_mm512_fmadd_ps(a31,h,t));
					    x3  = _mm512_fmadd_ps(a32,k2,_mm512_fmadd_ps(a31,k1,x));
					    tfv = fv(tt3,x3);
					    w3  = _mm512_mul_ps(vran3,_mm512_sqrt_ps(_mm512_mul_ps(q3,qh)));
					    tgv = gv(tt3,x3);
					    tmp = _mm512_mul_ps(h,_mm512_mul_ps(tgv,w3));
					    k3  = _mm512_fmadd_ps(h,tfv,tmp);
					    tt4 = _mm512_fmadd_ps(a43,h,_mm512_fmadd_ps(a42,h,_mm512_fmadd_ps(a41,h,t)));
					    x4  = _mm512_fmadd_ps(a43,k3,_mm512_fmadd_ps(a42,k2,_mm512_fmadd_ps(a41,k1,x)));
					    tgv = gv(tt4,x4);
					    w4  = _mm512_mul_ps(vran4,_mm512_sqrt_ps(_mm512_mul_ps(q4,qh)));
					    tfv = fv(tt4,x4);
					    tmp = _mm512_mul_ps(h,_mm512_mul_ps(tgv,w4));
					    k4  = _mm512_fmadd_ps(h,tfv,tmp);
					    step = _mm512_fmadd_ps(a54,k4,_mm512_fmadd_ps(a53,k3,
					                       _mm512_fmadd_ps(a52,k2,_mm512_fmadd_ps(a51,k1,x))));
					    return (step);
				   }


				    
    }

}




#endif /*__GMS_STOCHASTIC_RK_AVX512_PS_HPP__*/
