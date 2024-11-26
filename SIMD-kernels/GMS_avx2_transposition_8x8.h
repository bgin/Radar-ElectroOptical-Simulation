


#ifndef __GMS_AVX2_TRANSPOSITION_8X8_HPP__
#define __GMS_AVX2_TRANSPOSITION_8X8_HPP__


#include <immintrin.h>
#include "GMS_config.h"



namespace gms {


          namespace math {

	             // In-place
	              __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     void transpose_ymm8r4_8x8_ip(__m256 &x0,
		                                   __m256 &x1,
						   __m256 &x2,
						   __m256 &x3,
						   __m256 &x4,
						   __m256 &x5,
						   __m256 &x6,
						   __m256 &x7); 
						   

	               // In-place
	              __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		       void transpose_u_ymm8r4_8x8_ip(float * __restrict x0,
		                                     float * __restrict x1,
						     float * __restrict x2,
						     float * __restrict x3,
						     float * __restrict x4,
						     float * __restrict x5,
						     float * __restrict x6,
						     float * __restrict x7); 


		    // In-place
	              __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		       void transpose_a_ymm8r4_8x8_ip(float * __restrict __ATTR_ALIGN__(32) x0,
		                                     float * __restrict __ATTR_ALIGN__(32) x1,
						     float * __restrict __ATTR_ALIGN__(32) x2,
						     float * __restrict __ATTR_ALIGN__(32) x3,
						     float * __restrict __ATTR_ALIGN__(32) x4,
						     float * __restrict __ATTR_ALIGN__(32) x5,
						     float * __restrict __ATTR_ALIGN__(32) x6,
						     float * __restrict __ATTR_ALIGN__(32) x7); 



		        // In-place
	              __ATTR_REGCALL__
                       __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_u_ymm8r4_8x8_ip(float * __restrict x); 


		          // In-place
	              __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_u_ymm8r4_8x8_ip(float * __restrict __ATTR_ALIGN__(32) x); 



		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		    void transpose_u_ymm8r4_8x8(float * __restrict x,
		                                  float * __restrict y); 


		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      void transpose_a_ymm8r4_8x8(float * __restrict __ATTR_ALIGN__(32) x,
		                                  float * __restrict __ATTR_ALIGN__(32) y); 
		  }

#include <cstdint>


                      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     void transpose_u_ymm8r4_8x8_v2(float * __restrict x,
		                                     float * __restrict y,
						     const int32_t n); 




		      __ATTR_REGCALL__
                      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     void transpose_a_ymm8r4_8x8_v2(float * __restrict __ATTR_ALIGN__(32) x,
		                                     float * __restrict __ATTR_ALIGN__(32) y,
						     const int32_t n); 



		 

     }


}









#endif /*__GMS_AVX2_TRANSPOSITION_8X8_H__*/
