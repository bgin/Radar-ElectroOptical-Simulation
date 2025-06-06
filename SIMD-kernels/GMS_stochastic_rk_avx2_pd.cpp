





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



#include "GMS_stochastic_rk_avx2.h"



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
	                           
                                   
				    __m256d gms::math::rk1_ti_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran, // result of call to svrng_generate4_double(engine,normal)
							       __m256d (*fi) (const __m256d),
							       __m256d (*gi) (const __m256d)){
							       
							       
                                           register const  __m256d a21 = _mm256_set1_pd(1.0);
					   register const  __m256d q1  = _mm256_set1_pd(1.0);
					   register        __m256d w1;
					   register        __m256d step;
					   register        __m256d k1;
					   register        __m256d tmp0;
					   register        __m256d tmp1;
					   register        __m256d tmp2;
					   register        __m256d tmp3;
					   tmp2 = gi(x);
		                           tmp1  = fi(x);                   
                                  	   tmp0   = _mm256_mul_pd(q1,_mm256_div_pd(q,h));
					   w1     = _mm256_mul_pd(vran,tmp0);
					   tmp3   = _mm256_mul_pd(h,_mm256_mul_pd(tmp2,w1));
					   k1     = _mm256_fmadd_pd(h,tmp1,tmp3);
					   step   = _mm256_fmadd_pd(a21,k1,x);
					   return (step);
				 }

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

			         
				    __m256d gms::math::rk2_ti_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran1, // result of call to svrng_generate4_double(engine,normal)
							       const __m256d vran2, // result of call to svrng_generate4_double(engine,normal)
							       __m256d (*fi) (const __m256d),
							       __m256d (*gi) (const __m256d)){
							      
                                        register const __m256d _0  = _mm256_setzero_pd();
                                        register const __m256d a21 = _mm256_set1_pd(1.0);
					register const __m256d a31 = _mm256_set1_pd(0.5);
					register const __m256d a32 = a31;
					register const __m256d q1  = _mm256_set1_pd(2.0);
					register const __m256d q2  = q1;
					register __m256d tfi       = _0;
					register __m256d tgi       = _0;
					register __m256d w1        = _0;
					register __m256d w2        = _0;
					register __m256d x2        = _0;
					register __m256d k1        = _0;
					register __m256d k2        = _0;
					register __m256d t0        = _0;
					register __m256d t1        = _0;
					register __m256d t2        = _0;
					register __m256d step      = _0;
				        
					tfi = fi(x);
					t0  = _mm256_sqrt_pd(_mm256_mul_pd(q1,_mm256_div_pd(q,h)));
					tgi  = gi(x);
					w1   = _mm256_mul_pd(vran1,t0);
					t1   = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w1));
					k1   = _mm256_fmadd_pd(h,tfi,t1);
					x2   = _mm256_fmadd_pd(a21,k1,x);
					tfi  = fi(x2);
					t0   = _mm256_sqrt_pd(_mm256_mul_pd(q2,_mm256_div_pd(q,h)));
					w2   = _mm256_mul_pd(vran2,t0);
					tgi  = gi(x2);
					t2   = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w2));
					k2   = _mm256_fmadd_pd(h,tfi,t2);
					step = _mm256_fmadd_pd(a32,k2,_mm256_fmadd_pd(a31,k1,x));
					return (step);
				}

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

				  
				   __m256d gms::math::rk3_ti_step_ymm4r8(const __m256d x,
				                              const __m256d t,
							      const __m256d h,
							      const __m256d q,
							      const __m256d vran1,
							      const __m256d vran2,
							      const __m256d vran3,
							      __m256d (*fi) (const __m256d),
							      __m256d (*gi) (const __m256d)) {
							     
 
                                         const __m256d _0      = _mm256_setzero_pd();
					 const __m256d a21     = _mm256_set1_pd(1.52880952525675);
					 const __m256d a31     = _0;
					 const __m256d a32     = _mm256_set1_pd(0.51578733443615);
					 const __m256d a41     = _mm256_set1_pd(0.53289582961739);
					 const __m256d a42     = _mm256_set1_pd(0.25574324768195);
					 const __m256d a43     = _mm256_set1_pd(0.21136092270067);
					 const __m256d q1      = _mm256_set1_pd(1.87653936176981);
					 const __m256d q2      = _mm256_set1_pd(3.91017166264989);
					 const __m256d q3      = _mm256_set1_pd(4.73124353935667);
					 const __m256d qh      = _mm256_div_pd(q,h);
					 register __m256d k1   = _0;
					 register __m256d k2   = _0;
					 register __m256d k3   = _0;
					 register __m256d t0   = _0;
					 register __m256d t1   = _0;
					 register __m256d t2   = _0;
					 register __m256d w1   = _0;
					 register __m256d w2   = _0;
					 register __m256d w3   = _0;
					 register __m256d x2   = _0;
					 register __m256d x3   = _0;
					 register __m256d step = _0;
					 __m256d          tfi  = _0;
					 __m256d          tgi  = _0;
					 __m256d          tmp  = _0;
					 __m256d          tmp2 = _0;
					
		                      	 tfi = fi(x);
					 t0  = _mm256_sqrt_pd(_mm256_mul_pd(q1,qh));
					 t1   = _mm256_sqrt_pd(_mm256_mul_pd(q2,qh));
					 t2   = _mm256_sqrt_pd(_mm256_mul_pd(q3,qh));
					 tgi  = gi(x);
					 w1   = _mm256_mul_pd(vran1,t0);
					 tmp  = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w1));
					 k1   = _mm256_fmadd_pd(h,tfi,tmp);
					 tmp2 = _mm256_fmadd_pd(a41,k1,x);
					 x2   = _mm256_fmadd_pd(a21,k1,x);
					 tfi  = fi(x2);
					 w2   = _mm256_mul_pd(vran2,t1);
					 tgi  = gi(x2);
					 tmp  = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w2));
					 k2   = _mm256_fmadd_pd(h,tfi,tmp);
					 x3   = _mm256_fmadd_pd(a32,k2,_mm256_fmadd_pd(a31,k1,x));
					 tfi  = fi(x3);
					 w3   = _mm256_mul_pd(vran3,t2);
					 tgi  = gi(x3);
					 tmp  = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w3));
					 k3   = _mm256_fmadd_pd(h,tfi,tmp);
					 step = _mm256_fmadd_pd(a43,k3,_mm256_fmadd_pd(a42,k2,tmp2));
					 return (step);
				   }

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

				
				   __m256d gms::math::rk4_ti_step_ymm4r8(const __m256d x,
				                              const __m256d t,
							      const __m256d h,
							      const __m256d q,
							      const __m256d vran1,
							      const __m256d vran2,
							      const __m256d vran3,
							      const __m256d vran4,
							      __m256d (*fi) (const __m256d),
							      __m256d (*gi) (const __m256d)) {

                                           const __m256d _0    = _mm256_setzero_pd();
					   const __m256d a21   = _mm256_set1_pd(2.71644396264860);
					   const __m256d a31   = _mm256_set1_pd(-6.95653259006152);
					   const __m256d a32   = _mm256_set1_pd(0.78313689457981);
					   const __m256d a41   = _0;
					   const __m256d a42   = _mm256_set1_pd(0.48257353309214);
					   const __m256d a43   = _mm256_set1_pd(0.26171080165848);
					   const __m256d a51   = _mm256_set1_pd(0.47012396888046);
					   const __m256d a52   = _mm256_set1_pd(0.36597075368373);
					   const __m256d a53   = _mm256_set1_pd(0.08906615686702);
					   const __m256d a54   = _mm256_set1_pd(0.07483912056879);
					   const __m256d q1    = _mm256_set1_pd(2.12709852335625);
					   const __m256d q2    = _mm256_set1_pd(2.73245878238737);
					   const __m256d q3    = _mm256_set1_pd(11.22760917474960);
					   const __m256d q4    = _mm256_set1_pd(13.36199560336697);
					   const __m256d qh    = _mm256_div_pd(q,h);
					   register __m256d k1 = _0;
					   register __m256d k2 = _0;
					   register __m256d k3 = _0;
					   register __m256d k4 = _0;
					   register __m256d t1 = _0;
					   register __m256d t2 = _0;
					   register __m256d t3 = _0;
					   register __m256d t4 = _0;
					   register __m256d w1 = _0;
					   register __m256d w2 = _0;
					   register __m256d w3 = _0;
					   register __m256d w4 = _0;
					   register __m256d x2 = _0;
					   register __m256d x3 = _0;
					   register __m256d x4 = _0;
					   register __m256d step = _0;
					   __m256d          tfi  = _0;
					   __m256d          tgi  = _0;
					   __m256d          tmp1 = _0;
					   __m256d          tmp2 = _0;
					   tfi  = fi(x);
					   w1   = _mm256_mul_pd(vran1,_mm256_sqrt_pd(_mm256_mul_pd(q1,qh)));
					   tgi  = gi(x);
					   tmp1 = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w1));
					   k1   = _mm256_fmadd_pd(h,tfi,tmp1);
					   x2   = _mm256_fmadd_pd(a21,k1,x);
					   tfi  = fi(x2);
					   w2   = _mm256_mul_pd(vran2,_mm256_sqrt_pd(_mm256_mul_pd(q2,qh)));
					   tgi  = gi(x2);
					   tmp1 = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w2));
					   k2   = _mm256_fmadd_pd(h,tfi,tmp1);
					   x3   = _mm256_fmadd_pd(a32,k2,_mm256_fmadd_pd(a31,k1,x1));
					   tfi  = fi(x3);
					   w3   = _mm256_mul_pd(vran3,_mm256_sqrt_pd(_mm256_mul_pd(q3,qh)));
					   tgi  = gi(x3);
					   tmp1 = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w3));
					   k3   = _mm256_fmadd_pd(h,tfi,tmp1);
					   x4   = _mm256_fmadd_pd(a42,k2,_mm256_fmadd_pd(a41,k1,x1));
					   tfi  = fi(x4);
					   w4   = _mm256_mul_pd(vran4,_mm256_sqrt_pd(_mm256_mul_pd(q4,qh)));
					   tgi  = gi(x4);
					   tmp1 = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w4));
					   k4   = _mm256_fmadd_pd(h,tfi,tmp1);
					   step = _mm256_fmadd_pd(a51,k1,x1);
					   step = _mm256_fmadd_pd(a52,k2,step);
					   step = _mm256_fmadd_pd(a53,k3,step);
					   step = _mm256_fmadd_pd(a54,k4,step);
					   return (step);
				  }

/*
            The Runge-Kutta scheme is fourth-order, and suitable for time-varying
    systems.

    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
*/


				    __m256d gms::math::rk1_tv_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran,
							       __m256d (*fi) (const __m256d, const __m256d),
							       __m256d (*gi) (const __m256d, const __m256d) ) {
							       
							       
                                           register const  __m256d a21 = _mm256_set1_pd(1.0);
					   register const  __m256d q1  = _mm256_set1_pd(1.0);
					   register        __m256d w1;
					   register        __m256d step;
					   register        __m256d k1;
					   register        __m256d tmp0;
					   register        __m256d tmp1;
					   register        __m256d tmp2;
					   register        __m256d tmp3;
					 
		                           tmp1  = fi(t,x);                   
                                  	   tmp0   = _mm256_mul_pd(q1,_mm256_div_pd(q,h));
					   w1     = _mm256_mul_pd(vran,tmp0);
					   tmp2  = gi(t,x);
					   tmp3   = _mm256_mul_pd(h,_mm256_mul_pd(tmp2,w1));
					   k1     = _mm256_fmadd_pd(h,tmp1,tmp3);
					   step   = _mm256_fmadd_pd(a21,k1,x);
					   return (step);
				 }


				 
                                  
				    __m256d gms::math::rk2_tv_step_ymm4r8(const __m256d x,
				                               const __m2556d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran1,
							       const __m256d vran2,
							       __m256d (*fi) (const __m256d, const __m256d),
							       __m256d (*gi) (const __m256d, const __m256d)){
							      
                                        register const __m256d _0  = _mm256_setzero_pd();
                                        register const __m256d a21 = _mm256_set1_pd(1.0);
					register const __m256d a31 = _mm2556_set1_pd(0.5);
					register const __m256d a32 = a31;
					register const __m256d q1  = _mm256_set1_pd(2.0);
					register const __m256d q2  = q1;
					const __m256d          qh  = _mm256_div_pd(q,h);
					register __m256d tfi       = _0;
					register __m256d tgi       = _0;
					register __m256d w1        = _0;
					register __m256d w2        = _0;
					register __m256d x2        = _0;
					register __m256d k1        = _0;
					register __m256d k2        = _0;
					register __m256d t0        = _0;
					register __m256d t1        = _0;
					register __m256d t2        = _0;
					register __m256d tt        = _0;
					register __m256d step      = _0;
				        
					tfi  = fi(t,x);
					t0   = _mm256_sqrt_pd(_mm256_mul_pd(q1,qh));
					tgi  = gi(t,x);
					w1   = _mm256_mul_pd(vran1,t0);
					t1   = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w1));
					k1   = _mm256_fmadd_pd(h,tfi,t1);
					x2   = _mm256_fmadd_pd(a21,k1,x);
					tt   = _mm256_fmadd_pd(a21,h,t1);
					tfi  = fi(tt,x2);
					t0   = _mm256_sqrt_pd(_mm256_mul_pd(q2,qh));
					w2   = _mm256_mul_pd(vran2,t0);
					tgi  = gi(tt,x2);
					t2   = _mm256_mul_pd(h,_mm256_mul_pd(tgi,w2));
					k2   = _mm256_fmadd_pd(h,tfi,t2);
					step = _mm256_fmadd_pd(a32,k2,_mm256_fmadd_pd(a31,k1,x));
					return (step);
				}


	                       
                                  
				    __m256d gms::math::rk4_tv_step_ymm4r8(const __m256d x,
				                               const __m256d t,
							       const __m256d h,
							       const __m256d q,
							       const __m256d vran1,
							       const __m256d vran2,
							       const __m256d vran3,
							       const __m256d vran4,
							       __m256d (*fv) (const __m256d, const __m256d),
							       __m256d (*gv) (const __m256d, const __m256d)){

					    const __m256d _0   = _mm256_setzero_pd();
					    const __m256d a21  = _mm256_set1_pd(0.66667754298442);
					    const __m256d a31  = _mm256_set1_pd(0.63493935027993);
					    const __m256d a32  = _mm256_set1_pd(0.00342761715422);
					    const __m2556d a41  = _mm256_set1_pd(-2.32428921184321);
					    const __m256d a42  = _mm256_set1_pd(2.69723745129487);
					    const __m256d a43  = _mm256_set1_pd(0.29093673271592);
					    const __m256d a51  = _mm256_set1_pd(0.25001351164789);
					    const __m256d a52  = _mm256_set1_pd(0.67428574806272);
					    const __m256d a53  = _mm256_set1_pd(-0.00831795169360);
					    const __m256d a54  = _mm256_set1_pd(0.08401868181222);
					    const __m256d q1   = _mm256_set1_pd(3.99956364361748);
					    const __m256d q2   = _mm256_set1_pd(1.64524970733585);
					    const __m256d q3   = _mm256_set1_pd(1.59330355118722);
					    const __m256d q4   = _mm256_set1_pd(0.26330006501868);
					    const __m256d qh   = _mm256_div_pd(q,h);
					    register __m256d k1 = _0;
					    register __m256d k2 = _0;
					    register __m256d k3 = _0;
					    register __m256d k4 = _0;
					    register __m256d t2 = _0;
					    register __m256d t3 = _0;
					    register __m256d t4 = _0;
					    register __m256d x2 = _0;
					    register __m256d x3 = _0;
					    register __m256d x4 = _0;
					    register __m256d tt2=_0;
					    register __m256d tt3=_0;
					    register __m256d tt4=_0;
					    __m256d          tgv=_0;
					    __m256d          tfv=_0;
					    __m256d         step=_0;
					    __m256d         tmp =_0;
					    tfv = fv(t,x);
					    w1  = _mm256_mul_pd(vran1,_mm256_sqrt_pd(_mm256_mul_pd(q1,qh)));
					    tgv = gv(t,x);
					    tmp = _mm256_mul_pd(h,_mm256_mul_pd(tgv,w1));
					    k1  = _mm256_fmadd_pd(h,tfv,tmp);
					    tt2 = _mm256_fmadd_pd(a21,h,t);
					    x2  = _mm256_fmadd_pd(a21,k1,x);
					    tfv = fv(tt2,x2);
					    w2  = _mm256_mul_pd(vran2,_mm256_sqrt_pd(_mm256_mul_pd(q2,qh)));
					    tgv = gv(tt2,x2);
					    tmp = _mm256_mul_pd(h,_mm256_mul_pd(tgv,w2));
					    k2  = _mm256_fmadd_pd(h,tfv,tmp);
					    tt3 = _mm256_fmadd_pd(a32,h,_mm256_fmadd_pd(a31,h,t));
					    x3  = _mm256_fmadd_pd(a32,k2,_mm256_fmadd_pd(a31,k1,x));
					    tfv = fv(tt3,x3);
					    w3  = _mm256_mul_pd(vran3,_mm256_sqrt_pd(_mm256_mul_pd(q3,qh)));
					    tgv = gv(tt3,x3);
					    tmp = _mm256_mul_pd(h,_mm256_mul_pd(tgv,w3));
					    k3  = _mm256_fmadd_pd(h,tfv,tmp);
					    tt4 = _mm256_fmadd_pd(a43,h,_mm256_fmadd_pd(a42,h,_mm256_fmadd_pd(a41,h,t)));
					    x4  = _mm256_fmadd_pd(a43,k3,_mm256_fmadd_pd(a42,k2,_mm256_fmadd_pd(a41,k1,x)));
					    tgv = gv(tt4,x4);
					    w4  = _mm256_mul_pd(vran4,_mm256_sqrt_pd(_mm256_mul_pd(q4,qh)));
					    tfv = fv(tt4,x4);
					    tmp = _mm256_mul_pd(h,_mm256_mul_pd(tgv,w4));
					    k4  = _mm256_fmadd_pd(h,tfv,tmp);
					    step = _mm256_fmadd_pd(a54,k4,_mm256_fmadd_pd(a53,k3,
					                       _mm256_fmadd_pd(a52,k2,_mm256_fmadd_pd(a51,k1,x))));
					    return (step);
				   }


				    
   
