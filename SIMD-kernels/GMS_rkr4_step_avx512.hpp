

#ifndef __GMS_RKR4_STEP_AVX512_HPP__
#define __GMS_RKR4_STEP_AVX512_HPP__ 020620221447

/*
     Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_ralston_4.html
     Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_RKR2_STEP_AVX512_MAJOR = 1U;
    const unsigned int GMS_RKR2_STEP_AVX512_MINOR = 0U;
    const unsigned int GMS_RKR2_STEP_AVX512_MICRO = 0U;
    const unsigned int GMS_RKR2_STEP_AVX512_FULLVER =
      1000U*GMS_RKR2_STEP_AVX512_MAJOR+
      100U*GMS_RKR2_STEP_AVX512_MINOR+
      10U*GMS_RKR2_STEP_AVX512_MICRO;
    const char * const GMS_RKR2_STEP_AVX512_CREATION_DATE = "02-06-2022 14:47 PM +00200 (THR 02 JUN 2022 GMT+2)";
    const char * const GMS_RKR2_STEPAVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RKR2_STEP_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RKR2_STEP_AVX512_DESCRIPTION   = "Vectorized (AVX512) Runge-Kutta-Ralston order 4 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {


           namespace math {

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The fourth order Ralston's method is a four stage Runge-Kutta method   //
//     for approximating the solution of the differential equation            //
//     y'(x) = f(x,y) with initial condition y(x0) = c.  The fourth order     //
//     Ralston's method evaluates f(x,y) four times per step. For step i+1,   //
//     y[i+1] = y[i] + (263 + 24s5)/1812  k1 + (125-1000s5)/3828 k2           //
//                   + 1024*(3346 + 1623s5)/5924787 k3 + (30 - 4s5)/123 k4,   //
//     where                                                                  //
//     k1 = h * f(x[i],y[i]),                                                 //
//     k2 = h * f(x[i]+2/5 h, y[i]+2/5 k1),                                   //
//     k3 = h * f(x[i] + (7/8 - 3s5/16)h, y[i] + (-2889+1428s5)/1024 k1       //
//                                             + (3785-1620s5)/1024 k2)       //
//     k4 = h * f(x[i] +h, y[i] + (-3365+2094s5)/6040 k1                      //
//                                (-975-3046s5)/2552 k2                       //
//                                (467040+203968s5)/2400845 k3) )             //
//     and s5 = sqrt(5), x[i] = x0 + i * h.                                   //
//                                                                            //
//     This version of the Runge-Kutta method is a fourth order method.       //
//     Richardson extrapolation can be used to increase the order and         //
//     accuracy.                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d
			rkr4_step_zmm8r8(__m512d(*f)(__m512d,
			                             __m512d),
					 __m512d y0,
					 __m512d x0,
					 const __m512d h,
					 const int32_t n) {

			     const __m512d a2    = _mm512_set1_pd(0.40);
			     const __m512d h2    = _mm512_mul_pd(a2,h);
			     const __m512d a3    = _mm512_set1_pd(0.4557372542187894323578);
			     const __m512d h3    = _mm512_mul_pd(a3,h);
			     const __m512d b31   = _mm512_set1_pd(0.2969776092477535968389);
			     const __m512d b32   = _mm512_set1_pd(0.1587596449710358355189);
			     const __m512d b41   = _mm512_set1_pd(0.2181003882259204667927);
			     const __m512d b42   = _mm512_set1_pd(−3.0509651486929308025876);
			     const __m512d b43   = _mm512_set1_pd(3.8328647604670103357948);
			     const __m512d g1    = _mm512_set1_pd(0.1747602822626903712242);
			     const __m512d g2    = _mm512_set1_pd(−0.5514806628787329399404);
			     const __m512d g3    = _mm512_set1_pd(1.2055355993965235343777);
			     const __m512d g4    = _mm512_set1_pd(0.1711847812195190343385);
			     __m512d k1,k2,k3,k4;
			     __m512d t0,t1,t2;

			     while(--n >= 0) {
                                   k1 = f(x0,y0);
				   k2 = f(_mm512_add_pd(x0,h2),
				          _mm512_fmadd_pd(h2,k1,y0));
				   t0 = _mm512_fmadd_pd(h,_mm512_fmadd_pd(b31,k2,_mm512_mul_pd(b32,k2)),y0);
				   k3 = f(_mm512_add_pd(x0,h3),t0);
				   x0 = _mm512_add_pd(x0,h);
				   t1 = _mm512_fmadd_pd(h,_mm512_fmadd_pd(b41,k1,
				                                      _mm512_fmadd_pd(b42,k2,
								                  _mm512_mul_pd(b43,k3)),y0));
				   k4 = f(x0,t0);
				   t2 = _mm512_fmadd_pd(g1,k1,
				                    _mm512_fmadd_pd(g2,k2,
						                _mm512_fmadd_pd(g3,k3,
								            _mm512_mul_pd(g4,k4))));
				   y0 = _mm512_fmadd_pd(h0,t2,y0);
			     }
			     return (y0);
		      }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			rkr4_step_zmm16r4(__m512(*f)(__m512,
			                             __m512),
					 __m512 y0,
					 __m512 x0,
					 const __m512 h,
					 const int32_t n) {

			     const __m512 a2    = _mm512_set1_ps(0.40F);
			     const __m512 h2    = _mm512_mul_ps(a2,h);
			     const __m512 a3    = _mm512_set1_ps(0.4557372542187894323578F);
			     const __m512 h3    = _mm512_mul_ps(a3,h);
			     const __m512 b31   = _mm512_set1_ps(0.2969776092477535968389F);
			     const __m512 b32   = _mm512_set1_ps(0.1587596449710358355189F);
			     const __m512 b41   = _mm512_set1_ps(0.2181003882259204667927F);
			     const __m512 b42   = _mm512_set1_ps(−3.0509651486929308025876F);
			     const __m512 b43   = _mm512_set1_ps(3.8328647604670103357948F);
			     const __m512 g1    = _mm512_set1_ps(0.1747602822626903712242F);
			     const __m512 g2    = _mm512_set1_ps(−0.5514806628787329399404F);
			     const __m512 g3    = _mm512_set1_ps(1.2055355993965235343777F);
			     const __m512 g4    = _mm512_set1_ps(0.1711847812195190343385F);
			     __m512 k1,k2,k3,k4;
			     __m512 t0,t1,t2;

			     while(--n >= 0) {
                                   k1 = f(x0,y0);
				   k2 = f(_mm512_add_ps(x0,h2),
				          _mm512_fmadd_ps(h2,k1,y0));
				   t0 = _mm512_fmadd_ps(h,_mm512_fmadd_ps(b31,k2,_mm512_mul_ps(b32,k2)),y0);
				   k3 = f(_mm512_add_ps(x0,h3),t0);
				   x0 = _mm512_add_ps(x0,h);
				   t1 = _mm512_fmadd_ps(h,_mm512_fmadd_ps(b41,k1,
				                                      _mm512_fmadd_ps(b42,k2,
								                  _mm512_mul_ps(b43,k3)),y0));
				   k4 = f(x0,t0);
				   t2 = _mm512_fmadd_ps(g1,k1,
				                    _mm512_fmadd_ps(g2,k2,
						                _mm512_fmadd_ps(g3,k3,
								            _mm512_mul_ps(g4,k4))));
				   y0 = _mm512_fmadd_ps(h0,t2,y0);
			     }
			     return (y0);
		      }

	   
     }


}






#endif /*__GMS_RKR4_STEP_AVX512_HPP__*/
