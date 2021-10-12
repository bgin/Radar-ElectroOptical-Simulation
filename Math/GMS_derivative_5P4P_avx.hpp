
#ifndef __GMS_DERIVATIVE_5P4P_AVX_HPP__
#define __GMS_DERIVATIVE_5P4P_AVX_HPP__ 101220211522



namespace file_version {

    const unsigned int gGMS_DERIVATIVE_5P4P_AVX_MAJOR = 1U;
    const unsigned int gGMS_DERIVATIVE_5P4P_AVX_MINOR = 0U;
    const unsigned int gGMS_DERIVATIVE_5P4P_AVX_MICRO = 0U;
    const unsigned int gGMS_DERIVATIVE_5P4P_AVX_FULLVER =
      1000U*gGMS_DERIVATIVE_5P4P_AVX_MAJOR+
      100U*gGMS_DERIVATIVE_5P4P_AVX_MINOR+
      10U*gGMS_DERIVATIVE_5P4P_AVX_MICRO;
    const char * const pgGMS_DERIVATIVE_5P4P_AVX_CREATION_DATE = "10-12-2021 15:22  +00200 (TUE 12 OCT 2021 GMT+2)";
    const char * const pgGMS_DERIVATIVE_5P4P_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DERIVATIVE_5P4P_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DERIVATIVE_5P4P_AVX_DESCRIPTION   = "Vectorized (AVX/AVX2) derivative implementation (kernel)."

}


#include <immintrin.h>
#include <limits>
#include "GMS_config.h"


namespace  gms {

          namespace math {

                            namespace {
			    
                                 __ATTR_ALWAYS_INLINE__
                                 __ATTR_HOT__
				 __ATTR_ALIGN__(32)
				 __attribute__((regcall))
		                 static inline
				 __m256d  abs_ymm4r8(const __m256d x)  {
                                     const __m256d mask = _mm256_set1_pd(0x7FFFFFFFFFFFFFFF);
			             return (_mm256_and_pd(x,mask));
		                  }
			    }

	              /*
                            Central 5-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
	                __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__attribute__((regcall))
	                static inline
			__m256d stencil_5P_central_ymm4r8(__m256d (*f) (__m256d),
			                                  __m256d vx,
						          __m256d vh,
						          __m256d &verr_ro,
						          __m256d &verr_tr) {
						  
                               const __m256d vn0   = _mm256_setzero_pd();
                               const __m256d vn0_5 = _m256_set1_pd(0.5);
			       const __m256d vn1_3 = _mm256_set1_pd(0.33333333333333333333333333);
			       const __m256d vn4_3 = _mm256_set1_pd(1.33333333333333333333333);
			       const __m256d vn2   = _mm256_set1_pd(2.0);
			       const __m256d veps  = _mm256_set1_pd(std::numeric_limits<double>::epsilon());
			       __m256d vp1         = vn0;
			       __m256d vp2         = vn0;
			       __m256d vp1h        = vn0;
			       __m256d vp2h        = vn0;
			       __m256d vt0         = vn0;
			       __m256d vt1         = vn0;
			       __m256d vt2         = vn0;
			       __m256d vt3         = vn0;
			       __m256d vt4         = vn0;
			       __m256 vtmp0        = vn0;
			       __m256d vtmp1       = vn0;
			       __m256d vtmp2       = vn0;
			       __m256d vtmp3       = vn0;
			       vp1   = f(_mm256_sub_pd(vx,vh));
			       vp2   = f(_mm256_add_pd(vx,vh));
			       vt0   = _mm256_mul_pd(vn0_,_mm256_sub_pd(vp1,vp2));
			       vp1h  = f(_mm256_sub_pd(vx,_mm256_mul_pd(vh,vn0_5)));
			       vp2h  = f(_mm256_add_pd(vx,_mm256_mul_pd(vh,vn0_5)));
			       vtmp3 = _mm256_div_pd(_mm256_abs_pd(vx),vh);
			       vt1   = _mm256_sub_pd(
			                       _mm256_mul_pd(vn4_3,
					              _mm256_sub_pd(vp2h,vp1h)),
						                    _mm256_mul_pd(vn3_0,vt0));
			       
			       vt2   = _mm256_mul_pd(_mm256_add_pd(
			                                   _mm256_abs_pd(vp2),
							          _mm256_abs_pd(vp1)),veps);
			       vtmp0 = _mm256_div_pd(vt1,vh);
			       vtmp2 = _mm256_mul_pd(_mm256_abs_pd(vp2h),
			                                       _mm256_abs_pd(vp1h));
			       vt3   = _mm256_fmadd_pd(vtmp2,veps,vt2);
			       vtmp1 = _mm256_div_pd(vt2,h);
			       vt4   = _mm256_mul_pd(_mm256_max_pd(vtmp0,vtmp1),
			                                  _mm256_mul_pd(vtmp3,veps));
								          
			       verr_ro = _mm256_mul_pd(
			                        abs_ymm4r8(
						         _mm256_div_pd(vt3,vh)),vt4);
			       verr_tr = abs_ymm4r8(
			                         _mm256_div_pd(
						            _mm256_sub_pd(vt1,vt0),vh));
			       return (vtmp0);
		       }


                       /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__attribute__((regcall))
	                static inline
			__m256d stencil_4P_forward_ymm4r8(__m256d (*f)(__m256d),
			                                  __m256d vx,
						          __m256d vh,
						          __m256d &verr_ro,
						          __m256d &verr_tr) {
                                 const __m256d vn0    = _mm256_setzero_pd();
				 const __m256d vn1_4  = _mm256_set1_pd(0.25);
				 const __m256d vn1_2  = _mm256_set1_pd(0.5);
				 const __m256d vn3_4  = _mm256_set1_pd(0.75);
				 const __m256d vn22_3 = _mm256_set1_pd(7.3333333333333333333333);
				 const __m256d vn62_3 = _mm256_set1_pd(20.6666666666666666666667);
				 const __m256d vn52_3 = _mm256_set1_pd(17.3333333333333333333333);
				 const __m256d vn2    = _mm256_set1_pd(2.0);
				 const __m256d vn4134 = _mm256_set1_pd(41.34);
				 const __m256d veps   = _mm256_set1_pd(std::numeric_limits<double>::epsilon());
				 __m256d dydx = vn0;
				 __m256d vp1  = vn0;
				 __m256d vp2  = vn0;
				 __m256d vp3  = vn0;
				 __m256d vp4  = vn0;
				 __m256d vt0  = vn0;
				 __m256d vt1  = vn0;
				 __m256d vt2  = vn0;
				 __m256d vt3  = vn0;
				 __m256d vtmp0 = vn0;
				 __m256d vtmp1 = vn0;
				 __m256d vtmp2 = vn0;
				 __m256d vtmp3 = vn0;
				 __m256d vtmp4 = vn0;
				 __m256d vtmp5 = vn0;
				 __m256d vtmp6 = vn0;
				 __m256d vtmp7 = vn0;
				 vp1 = f(_mm256_fmadd_pd(vh,vn1_4,vx));
				 vtmp7 = _mm256_div_pd(vx,vh);
				 vp2 = f(_mm256_fmadd_pd(vh,vn1_2,vx));
				 vtmp0 = _mm256_mul_pd(vv52_3,_mm256_sub_pd(vp2,vp1));
				 vp3 = f(_mm256_fmadd_pd(vh,vn3_4,vx));
				 vtmp1 = _mm256_mul_pd(vn62_3,_mm256_sub_pd(vp3,vp2));
				 vp4 = f(_mm256_add_pd(vx,vh));
				 vtmp2 = _mm256_mul_pd(vn22_3,_mm256_sub_pd(vp4,vp3));
				 vt0 = _mm256_mul_pd(vn2,_mm256_sub_pd(vp4,vp2));
				 vtmp5 = _mm256_div_pd(vt0,vh);
				 vt1 = _mm256_add_pd(_mm256_sub_pd(vtmp2,vtmp1),vtmp0);
				 vtmp6 = _mm256_div_pd(vt1,vh);
				 vtmp0 = abs_ymm4r8(vp4);
				 vtmp1 = abs_ymm4r8(vp3);
				 vtmp2 = abs_ymm4r8(vp2);
				 vtmp3 = abs_ymm4r8(vp1);
				 vtmp4 = _mm256_add_pd(_mm256_add_pd(vtmp0,vtmp1),
				                             _mm256_add_pd(vtmp2,vtmp3));
				 vt2   = _mm256_mul_pd(vn4134,
				                           _mm256_mul_pd(vtmp4,veps));
				 vtmp0 = _mm256_max_pd(_mm256_abs_pd(vtmp7),
				                              _mm256_abs_pd(vtmp6));
				 vt3   = _mm256_mul_pd(vtmp0,
				                     _mm256_mul_pd(vtmp7,veps));
				 dydx  = vtmp6;
				 verr_tr = abs_ymm4r8(_mm256_div_pd(
				                                _mm256_sub_pd(vt2,vt0),vh));
				 verr_ro = _mm256_add_pd(_mm256_abs_pd(
				                                _mm256_div_pd(vt2,vh),vt3));
				 return (dydx);
			}


			 /*
                            Backward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__attribute__((regcall))
	                static inline
			__m256d
			stencil_4P_backward_ymm4r8(__m256d (*f) (__m256d),
			                           __m256d vx,
					           __m256d vh,
						   __m256d &verr_ro,
						   __m256d &verr_tr ) {
                              const __m256d vxhn = _mm256_sub_pd(_mm256_setzero_pd(),vh);
			      return (stencil_4P_forward_ymm4r8(f,vx,vxhn,verr_ro,verr_tr));
			}
		       
      }

}














#endif /* __GMS_DERIVATIVE_5P4P_AVX_HPP__*/
