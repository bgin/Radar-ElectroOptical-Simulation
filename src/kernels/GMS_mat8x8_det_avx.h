

#ifndef __GMS_MAT8X8_DET_AVX_H__
#define __GMS_MAT8X8_DET_AVX_H__

#include <immintrin.h>
#include "GMS_config.h"

namespace gms {


         namespace  math {


#define NUM_MATS 8
#define NUM_MATS_HALF 4
#define MAT_WIDTH 4
#define MAT_HEIGHT 4
#define MAT_SIZE (MAT_WIDTH * MAT_HEIGHT)
	              /*
                           Based on Intel example.
                       */
                     
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     
		      void
		      mat8x8_det_a_ymm8r4(const float * __restrict input,
		                          float * __restrict output); 
		                          

	   // Unaligned version
                   
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		    
		      void
		      mat8x8_det_u_ymm8r4(const float * __restrict input,
		                          float * __restrict output); 
	   
	            
     }

}

















#endif /*__GMS_MAT8X8_DET_AVX_H__*/
