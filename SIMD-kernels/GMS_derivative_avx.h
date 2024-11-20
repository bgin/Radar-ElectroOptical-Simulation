
#ifndef __GMS_DERIVATIVE_AVX_HPP__
#define __GMS_DERIVATIVE_AVX_HPP__ 011120231351

namespace file_version {

    const unsigned int gGMS_DERIVATIVE_AVX_MAJOR = 1U;
    const unsigned int gGMS_DERIVATIVE_AVX_MINOR = 0U;
    const unsigned int gGMS_DERIVATIVE_AVX_MICRO = 0U;
    const unsigned int gGMS_DERIVATIVE_AVX_FULLVER =
      1000U*gGMS_DERIVATIVE_AVX_MAJOR+
      100U*gGMS_DERIVATIVE_AVX_MINOR+
      10U*gGMS_DERIVATIVE_AVX_MICRO;
    const char * const pgGMS_DERIVATIVE_AVX_CREATION_DATE = "01-11-2023 13:51 AM +00200 (WED 01 nov 2023 GMT+2)";
    const char * const pgGMS_DERIVATIVE_AVX_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DERIVATIVE_AVX_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DERIVATIVE_AVX_DESCRIPTION   = "Vectorized (AVX) derivative implementation."

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
	               	__m256d stencil_5P_ymm4r8(__m256d (*f) (__m256d),
			                          __m256d vx,
						  __m256d vh,
						  __m256d &verr_ro,
						  __m256d &verr_tr); 


		     
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	             	__m256d
			stencil_5P_central_ymm4r8_optim(__m256d (*f) (__m256d),
			                        __m256d vx,
						__m256d vh,
						__m256d &vabserr); 
						

			  /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	                __m256d stencil_4P_ymm4r8(__m256d (*f)(__m256d),
			                          __m256d vx,
						  __m256d vh,
						  __m256d &verr_ro,
						  __m256d &verr_tr); 


			
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	              	__m256d
			stencil_4P_forward_ymm4r8_optim(__m256d (*f) (__m256d),
			                                __m256d vx,
							__m256d vh,
							__m256d &vabserr); 

			  /*
                            Backward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	             	__m256d
			stencil_4P_backward_ymm4r8_optim(__m256d (*f) (__m256d),
			                                 __m256d vx,
							 __m256d vh,
							 __m256d &vabserr); 

     }

}


#endif /*__GMS_DERIVATIVE_AVX_H__*/
