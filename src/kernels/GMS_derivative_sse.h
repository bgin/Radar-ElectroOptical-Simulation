
#ifndef __GMS_DERIVATIVE_SSE_H__
#define __GMS_DERIVATIVE_SSE_H__ 011120231451

namespace file_version {

    const unsigned int gGMS_DERIVATIVE_SSE_MAJOR = 1U;
    const unsigned int gGMS_DERIVATIVE_SSE_MINOR = 0U;
    const unsigned int gGMS_DERIVATIVE_SSE_MICRO = 0U;
    const unsigned int gGMS_DERIVATIVE_SSE_FULLVER =
      1000U*gGMS_DERIVATIVE_SSE_MAJOR+
      100U*gGMS_DERIVATIVE_SSE_MINOR+
      10U*gGMS_DERIVATIVE_SSE_MICRO;
    const char * const pgGMS_DERIVATIVE_SSE_CREATION_DATE = "01-11-2023 14:51 AM +00200 (WED 01 nov 2023 GMT+2)";
    const char * const pgGMS_DERIVATIVE_SSE_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const pgGMS_DERIVATIVE_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const pgGMS_DERIVATIVE_SSE_DESCRIPTION   = "Vectorized (SSE) derivative implementation."

}


#include <immintrin.h>
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
	            	__m128d stencil_5P_xmm2r8(__m128d (*f) (__m128d),
			                          __m128d vx,
						  __m128d vh,
						  __m128d &verr_ro,
						  __m128d &verr_tr); 


		       
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	              	__m128d
			stencil_5P_central_xmm2r8_optim(__m128d (*f) (__m128d),
			                        __m128d vx,
						__m128d vh,
						__m128d &vabserr); 
						

			  /*
                            Forward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	               	__m128d stencil_4P_xmm2r8(__m128d (*f)(__m128d),
			                          __m128d vx,
						  __m128d vh,
						  __m128d &verr_ro,
						  __m128d &verr_tr); 

			
			
                        __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	             	__m128d
			stencil_4P_forward_xmm2r8_optim(__m128d (*f) (__m128d),
			                                __m128d vx,
							__m128d vh,
							__m128d &vabserr); 


			  /*
                            Backward derivative i.e. 4-point rule function derivative computation.
                            https://en.wikipedia.org/wiki/Five-point_stencil
                        */
			 __ATTR_HOT__
                        __ATTR_ALIGN__(32)
	               	__m128d
			stencil_4P_backward_xmm2r8_optim(__m128d (*f) (__m128d),
			                                 __m128d vx,
							 __m128d vh,
							 __m128d &vabserr); 

     }

}


#endif /*__GMS_DERIVATIVE_SSE_H__*/
