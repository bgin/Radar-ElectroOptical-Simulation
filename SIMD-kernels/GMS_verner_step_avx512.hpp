

#ifndef __GMS_VERNER_STEP_AVX512_HPP__
#define __GMS_VERNER_STEP_AVX512_HPP__ 050620221002


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



/*
   Adapted from: http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_verner.html
   Manually vectorized by Bernard Gingold, beniekg@gmail.com
*/

namespace file_version {

    const unsigned int GMS_VERNER_STEP_AVX512_MAJOR = 1U;
    const unsigned int GMS_VERNER_STEP_AVX512_MINOR = 0U;
    const unsigned int GMS_VERNER_STEP_AVX512_MICRO = 0U;
    const unsigned int GMS_VERNER_STEP_AVX512_FULLVER =
      1000U*GMS_VERNER_STEP_AVX512_MAJOR+
      100U*GMS_VERNER_STEP_AVX512_MINOR+
      10U*GMS_VERNER_STEP_AVX512_MICRO;
    const char * const GMS_VERNER_STEP_AVX512_CREATION_DATE = "05-06-2022 10:002 AM +00200 (SUN 05 JUN 2022 GMT+2)";
    const char * const GMS_VERNER_STEPAVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_VERNER_STEP_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_VERNER_STEP_AVX512_DESCRIPTION   = "Vectorized (AVX512) Runge-Kutta-Verner order 7 step."

}


#include <immintrin.h>
#include <cstdint>
#include "GMS_config.h"


namespace gms {

       namespace math {
////////////////////////////////////////////////////////////////////////////////
//  const register __m512d Runge_Kutta_Verner( const register __m512d (*f)(const register __m512d, const register __m512d), const register __m512d y0,        //
//                               const register __m512d x0, const register __m512d h, int number_of_steps );  //
//                                                                            //
//  Description:                                                              //
//     This routine uses the Runge-Kutta_Verner method described above to     //
//     approximate the solution at x = x0 + h * number_of_steps of the initial//
//     value problem y'=f(x,y), y(x0) = y0.                                   //
//                                                                            //
//  Arguments:                                                                //
//     const register __m512d *f                                                              //
//            Pointer to the function which returns the slope at (x,y) of the //
//            integral curve of the differential equation y' = f(x,y) which   //
//            passes through the point (x0,y0).                               //
//     const register __m512d y0                                                              //
//            The initial value of y at x = x0.                               //
//     const register __m512d x0                                                              //
//            The initial value of x.                                         //
//     const register __m512d h                                                               //
//            The step size.                                                  //
//     int    number_of_steps                                                 //
//            The number of steps. Must be a nonnegative integer.             //
//                                                                            //
//  Return Values:                                                            //
//     The solution of the initial value problem y' = f(x,y), y(x0) = y0 at   //
//     x = x0 + number_of_steps * h.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


                        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512d
			verner_step_zmm8r8(__m512d(*f)(__m512d,
			                             __m512d),
					    __m512d y0,
					    __m512d x0,
					    const __m512d h,
					    const int32_t n) {
#define sqrt21 4.58257569495584000680
                               const register __m512d c1 = _mm512_set1_pd(1.0 / 2.0);    
                               const register __m512d c2 = _mm512_set1_pd((7.0 + sqrt21 ) / 14.0);
                               const register __m512d c3 = _mm512_set1_pd((7.0 - sqrt21 ) / 14.0);

                               const register __m512d a21 =  _mm512_set1_pd(1.0 / 2.0);
                               const register __m512d a31 =  _mm512_set1_pd(1.0 / 4.0);
                               const register __m512d a32 =  _mm512_set1_pd(1.0 / 4.0);
                               const register __m512d a41 =  _mm512_set1_pd(1.0 / 7.0);
                               const register __m512d a42 =  _mm512_set1_pd(-(7.0 + 3.0 * sqrt21) / 98.0);
                               const register __m512d a43 =  _mm512_set1_pd((21.0 + 5.0 * sqrt21) / 49.0);
                               const register __m512d a51 =  _mm512_set1_pd((11.0 + sqrt21) / 84.0);
                               const register __m512d a53 =  _mm512_set1_pd((18.0 + 4.0 * sqrt21) / 63.0);
                               const register __m512d a54 =  _mm512_set1_pd((21.0 - sqrt21) / 252.0);
                               const register __m512d a61 =  _mm512_set1_pd((5.0 + sqrt21) / 48.0);
                               const register __m512d a63 =  _mm512_set1_pd((9.0 + sqrt21) / 36.0);
                               const register __m512d a64 =  _mm512_set1_pd((-231.0 + 14.0 * sqrt21) / 360.0);
                               const register __m512d a65 =  _mm512_set1_pd((63.0 - 7.0 * sqrt21) / 80.0);
                               const register __m512d a71 =  _mm512_set1_pd((10.0 - sqrt21) / 42.0);
                               const register __m512d a73 =  _mm512_set1_pd((-432.0 + 92.0 * sqrt21) / 315.0);
                               const register __m512d a74 =  _mm512_set1_pd((633.0 - 145.0 * sqrt21) / 90.0);
                               const register __m512d a75 =  _mm512_set1_pd((-504.0 + 115.0 * sqrt21) / 70.0);
                               const register __m512d a76 =  _mm512_set1_pd((63.0 - 13.0 * sqrt21) / 35.0);
                               const register __m512d a81 =  _mm512_set1_pd(1.0 / 14.0);
                               const register __m512d a85 =  _mm512_set1_pd((14.0 - 3.0 * sqrt21) / 126.0);
                               const register __m512d a86 =  _mm512_set1_pd((13.0 - 3.0 * sqrt21) / 63.0);
                               const register __m512d a87 =  _mm512_set1_pd(1.0 / 9.0);
                               const register __m512d a91 =  _mm512_set1_pd(1.0 / 32.0);
                               const register __m512d a95 =  _mm512_set1_pd((91.0 - 21.0 * sqrt21) / 576.0);
                               const register __m512d a96 =  _mm512_set1_pd(11.0 / 72.0);
                               const register __m512d a97 =  _mm512_set1_pd(-(385.0 + 75.0 * sqrt21) / 1152.0);
                               const register __m512d a98 =  _mm512_set1_pd((63.0 + 13.0 * sqrt21) / 128.0);
                               const register __m512d a10_1 =  _mm512_set1_pd(1.0 / 14.0);
                               const register __m512d a10_5 =  _mm512_set1_pd(1.0 / 9.0);
                               const register __m512d a10_6 =  _mm512_set1_pd(-(733.0 + 147.0 * sqrt21) / 2205.0);
                               const register __m512d a10_7 =  _mm512_set1_pd((515.0 + 111.0 * sqrt21) / 504.0);
                               const register __m512d a10_8 =  _mm512_set1_pd(-(51.0 + 11.0 * sqrt21) / 56.0);
                               const register __m512d a10_9 =  _mm512_set1_pd((132.0 + 28.0 * sqrt21) / 245.0);
                               const register __m512d a11_5 =  _mm512_set1_pd((-42.0 + 7.0 * sqrt21) / 18.0);
                               const register __m512d a11_6 =  _mm512_set1_pd((-18.0 + 28.0 * sqrt21) / 45.0);
                               const register __m512d a11_7 =  _mm512_set1_pd(-(273.0 + 53.0 * sqrt21) / 72.0);
                               const register __m512d a11_8 =  _mm512_set1_pd((301.0 + 53.0 * sqrt21) / 72.0);
                               const register __m512d a11_9 =  _mm512_set1_pd((28.0 - 28.0 * sqrt21) / 45.0);
                               const register __m512d a11_10 = _mm512_set1_pd((49.0 - 7.0 * sqrt21) / 18.0);

                               const register __m512d  b1  = _mm512_set1_pd(9.0 / 180.0);
                               const register __m512d  b8  = _mm512_set1_pd(49.0 / 180.0);
                               const register __m512d  b9  = _mm512_set1_pd(64.0 / 180.0);

                               const register __m512d c1h = _mm512_mul_pd(c1,h);
			       const register __m512d c2h = _mm512_mul_pd(c2,h);
			       const register __m512d c3h = _mm512_mul_pd(c3,h);
			       __m512d k1,k2,k3,k4,k,k6,k7,k8,k9,k10;
			      
			       while(--n >= 0) {
                                     k1 = _mm512_mul_pd(h,f(x0,y0));
				     k2 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c1h),
				            _mm512_fmadd_pd(a21,k1,y0)));
				     k3 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c1h),
				            _mm512_add_pd(y0,
					              _mm512_fmadd_pd(a31,k1,
						                  _mm512_mul_pd(a32,k2)))));
				     k4 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c2h),
				            		  _mm512_add_pd(y0,
						                 _mm512_fmadd_pd(a41,k1,
								             _mm512_fmadd_pd(a42,k2,
									                 _mm512_mul_pd(a43,k3))))));
				     k5 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c2h),
				             		   _mm512_add_pd(y0,
						                 _mm512_fmadd_pd(a51,k1,
								             _mm512_fmadd_pd(a53,k3,
									                 _mm512_mul_pd(a54,k4))))));
				     k6 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c1h),
				                          _mm512_add_pd(y0,
							         _mm512_fmadd_pd(a61,k1,
								           _mm512_fmadd_pd(a63,k3,
									           _mm512_fmadd_pd(a64,k4,
										              _mm512_mul_pd(a65,k5)))))));
				     k7 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c3h),
				                          _mm512_add_pd(y0,
							         _mm512_fmadd_pd(a71,k1,
								       _mm512_fmadd_pd(a73,k3,
								              _mm512_fmadd_pd(a74,k4,
									             _mm512_fmadd_pd(a75,k5,
										                _mm512_mul_pd(a76,k6))))))));
				     k8 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c3h),
				                          _mm512_add_pd(y0,
							          _mm512_fmadd_pd(a81,k1,
								           _mm512_fmadd_pd(a85,k5,
									           _mm512_fmadd_pd(a86,k6,
										              _mm512_mul_pd(a87,k7)))))));
				     k9 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c1h),
				                          _mm512_add_pd(y0,
							         _mm512_fmadd_pd(a91,k1,
								       _mm512_fmadd_pd(a95,k5,
								              _mm512_fmadd_pd(a96,k6,
									             _mm512_fmadd_pd(a97,k7,
										                _mm512_mul_pd(a98,k8))))))));
				     k10 = _mm512_mul_pd(h,f(_mm512_add_pd(x0,c2h),
				                           _mm512_add_pd(y0,
							           _mm512_fmadd_pd(a10_1,k1,
								         _mm512_fmadd_pd(a10_5,k5,
									       _mm512_fmadd_pd(a10_6,k6,
									             _mm512_fmadd_pd(a10_7,k7,
										          _mm512_fmadd_pd(a10_8,k8,
											        _mm512_mul_pd(a10_9,k9)))))))));
				     x0 = _mm512_add_pd(x0,h);
				     k10 = _mm512_mul_pd(h,f(x0,
				                           _mm512_add_pd(y0,
							           _mm512_fmadd_pd(a11_5,k5,
								         _mm512_fmadd_pd(a11_6,k6,
									       _mm512_fmadd_pd(a11_7,k7,
									             _mm512_fmadd_pd(a11_8,k8,
										          _mm512_fmadd_pd(a11_9,k9,
											        _mm512_mul_pd(a11_10,k10)))))))));
				     y0 = _mm512_add_pd(y0,
				                    _mm512_fmadd_pd(b1,k1,
						          _mm512_fmadd_pd(b8,k8,
							        _mm512_fmadd_pd(b9,k9,
								      _mm512_fmadd_pd(b8,k10,
								            _mm512_mul_pd(b1,k11))))));
			       }
			       return (y0);
		      }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_VECTORCALL__
	                static inline
			__m512
			verner_step_zmm16r4(__m512(*f)(__m512,
			                             __m512),
					    __m512 y0,
					    __m512 x0,
					    const __m512 h,
					    const int32_t n) {
#define sqrt21 4.58257569495584000680f
                               const register __m512 c1 = _mm512_set1_ps(1.0f / 2.0f);    
                               const register __m512 c2 = _mm512_set1_ps((7.0f + sqrt21 ) / 14.0f);
                               const register __m512 c3 = _mm512_set1_ps((7.0f - sqrt21 ) / 14.0f);

                               const register __m512 a21 =  _mm512_set1_ps(1.0f / 2.0f);
                               const register __m512 a31 =  _mm512_set1_ps(1.0f / 4.0f);
                               const register __m512 a32 =  _mm512_set1_ps(1.0f / 4.0f);
                               const register __m512 a41 =  _mm512_set1_ps(1.0f / 7.0f);
                               const register __m512 a42 =  _mm512_set1_ps(-(7.0f + 3.0f * sqrt21) / 98.0f);
                               const register __m512 a43 =  _mm512_set1_ps((21.0f + 5.0f * sqrt21) / 49.0f);
                               const register __m512 a51 =  _mm512_set1_ps((11.0f + sqrt21) / 84.0f);
                               const register __m512 a53 =  _mm512_set1_ps((18.0f + 4.0f * sqrt21) / 63.0f);
                               const register __m512 a54 =  _mm512_set1_ps((21.0f - sqrt21) / 252.0f);
                               const register __m512 a61 =  _mm512_set1_ps((5.0f + sqrt21) / 48.0f);
                               const register __m512 a63 =  _mm512_set1_ps((9.0f + sqrt21) / 36.0f);
                               const register __m512 a64 =  _mm512_set1_ps((-231.0f + 14.0f * sqrt21) / 360.0f);
                               const register __m512 a65 =  _mm512_set1_ps((63.0f - 7.0f * sqrt21) / 80.0f);
                               const register __m512 a71 =  _mm512_set1_ps((10.0f - sqrt21) / 42.0f);
                               const register __m512 a73 =  _mm512_set1_ps((-432.0f + 92.0f * sqrt21) / 315.0f);
                               const register __m512 a74 =  _mm512_set1_ps((633.0f - 145.0f * sqrt21) / 90.0f);
                               const register __m512 a75 =  _mm512_set1_ps((-504.0f + 115.0f * sqrt21) / 70.0f);
                               const register __m512 a76 =  _mm512_set1_ps((63.0f - 13.0f * sqrt21) / 35.0f);
                               const register __m512 a81 =  _mm512_set1_ps(1.0f / 14.0f);
                               const register __m512 a85 =  _mm512_set1_ps((14.0f - 3.0f * sqrt21) / 126.0f);
                               const register __m512 a86 =  _mm512_set1_ps((13.0f - 3.0f * sqrt21) / 63.0f);
                               const register __m512 a87 =  _mm512_set1_ps(1.0f / 9.0f);
                               const register __m512 a91 =  _mm512_set1_ps(1.0f / 32.0f);
                               const register __m512 a95 =  _mm512_set1_ps((91.0f - 21.0f * sqrt21) / 576.0f);
                               const register __m512 a96 =  _mm512_set1_ps(11.0f / 72.0f);
                               const register __m512 a97 =  _mm512_set1_ps(-(385.0f + 75.0f * sqrt21) / 1152.0f);
                               const register __m512 a98 =  _mm512_set1_ps((63.0f + 13.0f * sqrt21) / 128.0f);
                               const register __m512 a10_1 =  _mm512_set1_ps(1.0f / 14.0f);
                               const register __m512 a10_5 =  _mm512_set1_ps(1.0f / 9.0f);
                               const register __m512 a10_6 =  _mm512_set1_ps(-(733.0f + 147.0f * sqrt21) / 2205.0f);
                               const register __m512 a10_7 =  _mm512_set1_ps((515.0f + 111.0f * sqrt21) / 504.0f);
                               const register __m512 a10_8 =  _mm512_set1_ps(-(51.0f + 11.0f * sqrt21) / 56.0f);
                               const register __m512 a10_9 =  _mm512_set1_ps((132.0f + 28.0f * sqrt21) / 245.0f);
                               const register __m512 a11_5 =  _mm512_set1_ps((-42.0f + 7.0f * sqrt21) / 18.0f);
                               const register __m512 a11_6 =  _mm512_set1_ps((-18.0f + 28.0f * sqrt21) / 45.0f);
                               const register __m512 a11_7 =  _mm512_set1_ps(-(273.0f + 53.0f * sqrt21) / 72.0f);
                               const register __m512 a11_8 =  _mm512_set1_ps((301.0f + 53.0f * sqrt21) / 72.0f);
                               const register __m512 a11_9 =  _mm512_set1_ps((28.0f - 28.0f * sqrt21) / 45.0f);
                               const register __m512 a11_10 = _mm512_set1_ps((49.0f - 7.0f * sqrt21) / 18.0f);

                               const register __m512  b1  = _mm512_set1_ps(9.0f / 180.0f);
                               const register __m512  b8  = _mm512_set1_ps(49.0f / 180.0f);
                               const register __m512  b9  = _mm512_set1_ps(64.0f / 180.0f);

                               const register __m512 c1h = _mm512_mul_ps(c1,h);
			       const register __m512 c2h = _mm512_mul_ps(c2,h);
			       const register __m512 c3h = _mm512_mul_ps(c3,h);
			       __m512 k1,k2,k3,k4,k,k6,k7,k8,k9,k10;
			      
			       while(--n >= 0) {
                                     k1 = _mm512_mul_ps(h,f(x0,y0));
				     k2 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c1h),
				            _mm512_fmadd_ps(a21,k1,y0)));
				     k3 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c1h),
				            _mm512_add_ps(y0,
					              _mm512_fmadd_ps(a31,k1,
						                  _mm512_mul_ps(a32,k2)))));
				     k4 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c2h),
				            		  _mm512_add_ps(y0,
						                 _mm512_fmadd_ps(a41,k1,
								             _mm512_fmadd_ps(a42,k2,
									                 _mm512_mul_ps(a43,k3))))));
				     k5 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c2h),
				             		   _mm512_add_ps(y0,
						                 _mm512_fmadd_ps(a51,k1,
								             _mm512_fmadd_ps(a53,k3,
									                 _mm512_mul_ps(a54,k4))))));
				     k6 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c1h),
				                          _mm512_add_ps(y0,
							         _mm512_fmadd_ps(a61,k1,
								           _mm512_fmadd_ps(a63,k3,
									           _mm512_fmadd_ps(a64,k4,
										              _mm512_mul_ps(a65,k5)))))));
				     k7 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c3h),
				                          _mm512_add_ps(y0,
							         _mm512_fmadd_ps(a71,k1,
								       _mm512_fmadd_ps(a73,k3,
								              _mm512_fmadd_ps(a74,k4,
									             _mm512_fmadd_ps(a75,k5,
										                _mm512_mul_ps(a76,k6))))))));
				     k8 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c3h),
				                          _mm512_add_ps(y0,
							          _mm512_fmadd_ps(a81,k1,
								           _mm512_fmadd_ps(a85,k5,
									           _mm512_fmadd_ps(a86,k6,
										              _mm512_mul_ps(a87,k7)))))));
				     k9 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c1h),
				                          _mm512_add_ps(y0,
							         _mm512_fmadd_ps(a91,k1,
								       _mm512_fmadd_ps(a95,k5,
								              _mm512_fmadd_ps(a96,k6,
									             _mm512_fmadd_ps(a97,k7,
										                _mm512_mul_ps(a98,k8))))))));
				     k10 = _mm512_mul_ps(h,f(_mm512_add_ps(x0,c2h),
				                           _mm512_add_ps(y0,
							           _mm512_fmadd_ps(a10_1,k1,
								         _mm512_fmadd_ps(a10_5,k5,
									       _mm512_fmadd_ps(a10_6,k6,
									             _mm512_fmadd_ps(a10_7,k7,
										          _mm512_fmadd_ps(a10_8,k8,
											        _mm512_mul_ps(a10_9,k9)))))))));
				     x0 = _mm512_add_ps(x0,h);
				     k10 = _mm512_mul_ps(h,f(x0,
				                           _mm512_add_ps(y0,
							           _mm512_fmadd_ps(a11_5,k5,
								         _mm512_fmadd_ps(a11_6,k6,
									       _mm512_fmadd_ps(a11_7,k7,
									             _mm512_fmadd_ps(a11_8,k8,
										          _mm512_fmadd_ps(a11_9,k9,
											        _mm512_mul_ps(a11_10,k10)))))))));
				     y0 = _mm512_add_ps(y0,
				                    _mm512_fmadd_ps(b1,k1,
						          _mm512_fmadd_ps(b8,k8,
							        _mm512_fmadd_ps(b9,k9,
								      _mm512_fmadd_ps(b8,k10,
								            _mm512_mul_ps(b1,k11))))));
			       }
			       return (y0);
		      }


    }


}






#endif /*__GMS_VERNER_STEP_AVX512_HPP__*/
