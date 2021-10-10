
#ifndef __GMS_DERIVATIVE_AVX512_HPP__
#define __GMS_DERIVATIVE_AVX512_HPP__ 101020210944

namespace file_version {

    const unsigned int gGMS_DERIVATIVE_AVX512_MAJOR = 1U;
    const unsigned int gGMS_DERIVATIVE_AVX512_MINOR = 0U;
    const unsigned int gGMS_DERIVATIVE_AVX512_MICRO = 0U;
    const unsigned int gGMS_DERIVATIVE_AVX512_FULLVER =
      1000U*gGMS_DERIVATIVE_AVX512_MAJOR+
      100U*gGMS_DERIVATIVE_AVX512_MINOR+
      10U*gGMS_DERIVATIVE_AVX512_MICRO;
    const char * const pgGMS_DERIVATIVE_AVX512_CREATION_DATE = "10-10-2021 09:44 AM +00200 (SUN 10 OCT 2021 GMT+2)";
    const char * const pgGMS_DERIVATIVE_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DERIVATIVE_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DERIVATIVE_AVX512_DESCRIPTION   = "Vectorized (AVX512) derivative implementation."

}


#include <immintrin.h>
#include <limits>
#include "GMS_config.h"



namespace  gms {

            namespace math {

	               /*
                            Central 5-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
	                __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
			__m512d stencil_5P_zmm8r8(__m512d (*f) (__m512d),
			                          __m512d vx,
						  __m512d vh,
						  __m512d &verr_ro,
						  __m512d &verr_tr) {
						  
                               const __m512d vn0   = _mm512_setzero_pd();
                               const __m512d vn0_5 = _mm512_set1_pd(0.5);
			       const __m512d vn1_3 = _mm512_set1_pd(0.33333333333333333333333333);
			       const __m512d vn4_3 = _mm512_set1_pd(1.33333333333333333333333);
			       const __m512d vn2   = _mm512_set1_pd(2.0);
			       const __m512d veps  = _mm512_set1_pd(std::numeric_limits<double>::eps());
			       __m512d vp1         = vn0;
			       __m512d vp2         = vn0;
			       __m512d vp1h        = vn0;
			       __m512d vp2h        = vn0;
			       __m512d vt0         = vn0;
			       __m512d vt1         = vn0;
			       __m512d vt2         = vn0;
			       __m512d vt3         = vn0;
			       __m512d vt4         = vn0;
			       __m512d vtmp0       = vn0;
			       __m512d vtmp1       = vn0;
			       __m512d vtmp2       = vn0;
			       __m512d vtmp3       = vn0;
			       vp1   = f(_mm512_sub_pd(vx,vh));
			       vp2   = f(_mm512_add_pd(vx,vh));
			       vt0   = _mm512_mul_pd(vn0_,_mm512_sub_pd(vp1,vp2));
			       vp1h  = f(_mm512_sub_pd(vx,_mm512_mul_pd(vh,vn0_5)));
			       vp2h  = f(_mm512_add_pd(vx,_mm512_mul_pd(vh,vn0_5)));
			       vtmp3 = _mm512_div_pd(_mm512_abs_pd(vx),vh);
			       vt1   = _mm512_sub_pd(
			                       _mm512_mul_pd(vn4_3,
					              _mm512_sub_pd(vp2h,vp1h)),
						                    _mm512_mul_pd(vn3_0,vt0));
			       
			       vt2   = _mm512_mul_pd(_mm512_add_pd(
			                                   _mm512_abs_pd(vp2),
							          _mm512_abs_pd(vp1)),veps);
			       vtmp0 = _mm512_div_pd(vt1,vh);
			       vtmp2 = _mm512_mul_pd(_mm512_abs_pd(vp2h),
			                                       _mm512_abs_pd(vp1h));
			       vt3   = _mm512_fmadd_pd(vtmp2,veps,vt2);
			       vtmp1 = _mm512_div_pd(vt2,h);
			       vt4   = _mm512_mul_pd(_mm512_max_pd(vtmp0,vtmp1),
			                                  _mm512_mul_pd(vtmp3,veps));
								          
			       verr_ro = _mm512_mul_pd(
			                        _mm512_fabs_pd(
						         _mm512_div_pd(vt3,vh)),vt4);
			       verr_tr = _mm512_fabs_pd(
			                         _mm512_div_pd(
						            _mm512_sub_pd(vt1,vt0),vh));
			       return (vtmp0);
		       }


		        __ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
			__m512d
			stencil_5P_central_zmm8r8_optim(__m512d (*f) (__m512d),
			                        __m512d vx,
						__m512d vh,
						__m512d &vabserr) {
                            const __m512d vn0   = _mm512_setzero_pd();
                            const __m512d vn2   = _mm512_set1_pd(2.0);
			    const __m512d vn1_3 = _mm512_set1_pd(0.33333333333333333333333333333333);
			    const __m512d vn4   = _mm512_set1_pd(4.0);
			    __m512d vx0     = vn0;
			    __m512d verr_ro = vn0;
			    __m512d verr_tr = vn0;
			    __m512d verr    = vn0;
			    __mmask8 vb1    = 0;
			    __mmask8 vb2    = 0;
			    __mmask8 vb3    = 0;
			    __mmask8 vb4    = 0;
			    vx0 = stencil_5P_zmm8r8(f,vx,vh,verr_ro,verr_tr);
			    verr = _mm512_add_pd(verr_ro,verr_tr);
			    vb1  = _mm512_cmp_pd_mask(verr_ro,vn0,_CMP_GT_OQ) && 
			           _mm512_cmp_pd_mask(verr_tr,vn0,_CMP_GT_OQ);
			    vb2  = _mm512_cmp_pd_mask(verr_ro,verr_tr,_CMP_LT_OQ);
			    if(vb1 && vb2) {
                                __m512d vx02     = vn0;
			        __m512d verr_ro2 = vn0;
			        __m512d verr_tr2 = vn0;
			        __m512d verr2    = vn0;
				__m512d vh2      = vn0;
				__m512d vt0;
				__m512d tmp0;
				vt0   = _mm512_div_pd(verr_r0,
				                       _mm512_add_pd(verr_ro,verr_tr));

				vh2   = _mm512_mul_pd(vh,_mm512_pow_pd(vt0,vn1_3));
				vx02  = stencil_5P_zmm8r8(f,vx,vh,verr_ro2,verr_tr2);
				verr2 = _mm512_add_pd(verr_ro2,verr_tr2);
				vb3   = _mm512_cmp_pd_mask(verr2,verr,_CMP_LT_OQ);
				tmp0  = _mm512_abs_pd(_mm512_sub_pd(vx02,vx0));
				vb4   = _mm512_cmp_pd_mask(tmp0,
				                     _mm512_mul_pd(vn4,verr),_CMP_LT_OQ);
				if(vb3 && vb4) {
                                   vx0 = vx02;
				   verr = verr2;
				}
			    }
			    vabserr = err;
			    return (vx0);
			}


			  /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
			__m512d stencil_4P_zmm8r8(__m512d (*f)(__m512d),
			                          __m512d vx,
						  __m512d vh,
						  __m512d &verr_ro,
						  __m512d &verr_tr) {
                                 const __m512d vn0   = _mm512_setzero_pd();
				 const __m512d vn1_4  = _mm512_set1_pd(0.25);
				 const __m512d vn1_2  = _mm512_set1_pd(0.5);
				 const __m512d vn3_4  = _mm512_set1_pd(0.75);
				 const __m512d vn22_3 = _mm512_set1_pd(7.3333333333333333333333);
				 const __m512d vn62_3 = _mm512_set1_pd(20.6666666666666666666667);
				 const __m512d vn52_3 = _mm512_set1_pd(17.3333333333333333333333);
				 const __m512d vn2    = _mm512_set1_pd(2.0);
				 const __m512d vn4134 = _mm512_set1_pd(41.34);
				 const __m512d veps  = _mm512_set1_pd(std::numeric_limits<double>::eps());
				 __m512d dydx = vn0;
				 __m512d vp1  = vn0;
				 __m512d vp2  = vn0;
				 __m512d vp3  = vn0;
				 __m512d vp4  = vn0;
				 __m512d vt0  = vn0;
				 __m512d vt1  = vn0;
				 __m512d vt2  = vn0;
				 __m512d vt3  = vn0;
				 __m512d vtmp0 = vn0;
				 __m512d vtmp1 = vn0;
				 __m512d vtmp2 = vn0;
				 __m512d vtmp3 = vn0;
				 __m512d vtmp4 = vn0;
				 __m512d vtmp5 = vn0;
				 __m512d vtmp6 = vn0;
				 __m512d vtmp7 = vn0;
				 vp1 = f(_mm512_fmadd_pd(vh,vn1_4,vx));
				 vtmp7 = _mm512_div_pd(vx,vh);
				 vp2 = f(_mm512_fmadd_pd(vh,vn1_2,vx));
				 vtmp0 = _mm512_mul_pd(vv52_3,_mm512_sub_pd(vp2,vp1));
				 vp3 = f(_mm512_fmadd_pd(vh,vn3_4,vx));
				 vtmp1 = _mm512_mul_pd(vn62_3,_mm512_sub_pd(vp3,vp2));
				 vp4 = f(_mm512_add_pd(vx,vh));
				 vtmp2 = _mm512_mul_pd(vn22_3,_mm512_sub_pd(vp4,vp3));
				 vt0 = _mm512_mul_pd(vn2,_mm512_sub_pd(vp4,vp2));
				 vtmp5 = _mm512_div_pd(vt0,vh);
				 vt1 = _mm512_add_pd(_mm512_sub_pd(vtmp2,vtmp1),vtmp0);
				 vtmp6 = _mm512_div_pd(vt1,vh);
				 vtmp0 = _mm512_abs_pd(vp4);
				 vtmp1 = _mm512_abs_pd(vp3);
				 vtmp2 = _mm512_abs_pd(vp2);
				 vtmp3 = _mm512_abs_pd(vp1);
				 vtmp4 = _mm512_add_pd(_mm512_add_pd(vtmp0,vtmp1),
				                             _mm512_add_pd(vtmp2,vtmp3));
				 vt2   = _mm512_mul_pd(vn4134,
				                           _mm512_mul_pd(vtmp4,veps));
				 vtmp0 = _mm512_max_pd(_mm512_abs_pd(vtmp7),
				                              _mm512_abs_pd(vtmp6));
				 vt3   = _mm512_mul_pd(vtmp0,
				                     _mm512_mul_pd(vtmp7,veps));
				 dydx  = vtmp6;
				 verr_tr = _mm512_abs_pd(_mm512_div_pd(
				                                _mm512_sub_pd(vt2,vt0),vh));
				 verr_ro = _mm512_add_pd(_mm512_abs_pd(
				                                _mm512_div_pd(vt2,vh),vt3));
				 return (dydx);
			}


			
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
			__m512d
			stencil_4P_forward_zmm8r8_optim(__m512d (*f) (__m512d),
			                                __m512d vx,
							__m512d vh,
							__m512d &vabserr) {

                                 const __m512d vn0   = _mm512_setzero_pd();
				 const __m512d v1_2  = _mm512_set1_pd(0.5);
				 const __m512d vn4   = _mm512_set1_pd(4.0);
				 __m512d vx0     = vn0;
			         __m512d verr_ro = vn0;
			         __m512d verr_tr = vn0;
			         __m512d verr    = vn0;
			         __mmask8 vb1    = 0;
			         __mmask8 vb2    = 0;
			         __mmask8 vb3    = 0;
			         __mmask8 vb4    = 0;
			         vx0 = stencil_4P_zmm8r8(f,vx,vh,verr_ro,verr_tr);
			         verr = _mm512_add_pd(verr_ro,verr_tr);
			         vb1  = _mm512_cmp_pd_mask(verr_ro,vn0,_CMP_GT_OQ) && 
			                _mm512_cmp_pd_mask(verr_tr,vn0,_CMP_GT_OQ);
			         vb2  = _mm512_cmp_pd_mask(verr_ro,verr_tr,_CMP_LT_OQ);
			         if(vb1 && vb2) {
                                    __m512d vx02     = vn0;
			            __m512d verr_ro2 = vn0;
			            __m512d verr_tr2 = vn0;
			            __m512d verr2    = vn0;
				    __m512d vh2      = vn0;
				    __m512d vt0;
				    __m512d vtmp0;
				    __m512d vtmp1;
				    vtmp0 = _mm512_div_pd(verr_ro,verr_tr);
				    vh2   = _mm512_mul_pd(vh,tmp0);
				    vx02  = stencil_4P_zmm8r8(f,vx,vh,verr_r02,verr_tr2);
                                    verr2 = _mm512_add_pd(verr_ro2,verr_tr2);
				    vb3   = _mm512_cmp_pd_mask(verr2,verr,_CMP_LT_OQ);
				    vtmp1 = _mm512_abs_pd(_mm512_sub_pd(vx02,vx0));
				    vb4   = _mm512_cmp_pd_mask(vtmp1,_mm512_mul_pd(vn4,verr),_CMP_LT_OQ);
				    if(vb3 && vb4) {
                                       vx0 = vx02;
				       verr = verr2;
				    }
				}
				vabserr = verr;
				return (vx0);
		        }


			  /*
                            Backward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			__ATTR_ALWAYS_INLINE__
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                static inline
			__m512d
			stencil_4P_backward_zmm8r8_optim(__m512d (*f) (__m512d),
			                                 __m512d vx,
							 __m512d vh,
							 __m512d &vabserr) {
                              const __m512d vxhn = _mm512_sub_pd(_mm512_setzero_pd(),vh);
			      return (stencil_4P_forward_zmm8r8_optim(f,vx,vxhn,vabserr));
			}

     }

}


#endif /*__GMS_DERIVATIVE_AVX512_HPP__*/
