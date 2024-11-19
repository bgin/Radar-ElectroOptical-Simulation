
#ifndef __GMS_DERIVATIVE_AVX512_H__
#define __GMS_DERIVATIVE_AVX512_H__ 101020210944

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
	              
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
			__ATTR_REGCALL__
	             	__m512d stencil_5P_zmm8r8(__m512d (*f) (__m512d),
			                          __m512d vx,
						  __m512d vh,
						  __m512d &verr_ro,
						  __m512d &verr_tr); 


		      
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	             	__m512d
			stencil_5P_central_zmm8r8_optim(__m512d (*f) (__m512d),
			                        __m512d vx,
						__m512d vh,
						__m512d &vabserr); 


			  /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	              	__m512d stencil_4P_zmm8r8(__m512d (*f)(__m512d),
			                          __m512d vx,
						  __m512d vh,
						  __m512d &verr_ro,
						  __m512d &verr_tr); 


			
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	             	__m512d
			stencil_4P_forward_zmm8r8_optim(__m512d (*f) (__m512d),
			                                __m512d vx,
							__m512d vh,
							__m512d &vabserr);

			  /*
                            Backward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	              	__m512d
			stencil_4P_backward_zmm8r8_optim(__m512d (*f) (__m512d),
			                                 __m512d vx,
							 __m512d vh,
							 __m512d &vabserr);

     }

}


#endif /*__GMS_DERIVATIVE_AVX512_H__*/
