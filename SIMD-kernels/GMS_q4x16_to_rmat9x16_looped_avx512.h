

#ifndef __GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_H__
#define __GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_H__ 141120211425


namespace file_info {

   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MAJOR = 1;
   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MINOR = 0;
   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MICRO = 0;
   const unsigned int GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_FULLVER =
    1000*GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MAJOR+
    100*GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MINOR+
    10*GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_MICRO;
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_CREATE_DATE = "14-11-2021 14:25 +00200 (SUN 14 NOV 2021 GMT+2)";
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_BUILD_DATE  = __DATE__":"__TIME__;
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_DESCRIPTION = "Quaternion 4x16 to Rotation Matrix 9x16 kernel AVX512.";
}



#include <cstdint>
#include "GMS_config.h"



namespace gms {

          namespace  math {


                  
                    
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		   
                      void
                      q4x16_rm9x16_looped_u_nonunroll_zmm16r4(float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      const float * __restrict __ATTR_ALIGN__(64) q_x,
							      const float * __restrict __ATTR_ALIGN__(64) q_y,
							      const float * __restrict __ATTR_ALIGN__(64) q_z,
							      const float * __restrict __ATTR_ALIGN__(64) q_w,
							      const int32_t n); 



                    
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		    
                      void
                      q4x16_rm9x16_looped_a_nonunroll_zmm16r4(float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      float * __restrict __ATTR_ALIGN__(64) q_x,
							      float * __restrict __ATTR_ALIGN__(64) q_y,
							      float * __restrict __ATTR_ALIGN__(64) q_z,
							      float * __restrict __ATTR_ALIGN__(64) q_w,
							      const int32_t n); 


                     
		     
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		   
                      void
                      q4x16_rm9x16_looped_u_unroll4x_zmm16r4(float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      const float * __restrict __ATTR_ALIGN__(64) q_x,
							      const float * __restrict __ATTR_ALIGN__(64) q_y,
							      const float * __restrict __ATTR_ALIGN__(64) q_z,
							      const float * __restrict __ATTR_ALIGN__(64) q_w,
							      const int32_t n); 


		     
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     
                      void
                      q4x16_rm9x16_looped_a_unroll4x_zmm16r4( float * __restrict __ATTR_ALIGN__(64) row1,
		                                              float * __restrict __ATTR_ALIGN__(64) row2,
							      float * __restrict __ATTR_ALIGN__(64) row3,
							      float * __restrict __ATTR_ALIGN__(64) row4,
							      float * __restrict __ATTR_ALIGN__(64) row5,
							      float * __restrict __ATTR_ALIGN__(64) row6,
							      float * __restrict __ATTR_ALIGN__(64) row7,
							      float * __restrict __ATTR_ALIGN__(64) row8,
							      float * __restrict __ATTR_ALIGN__(64) row9,
							      float * __restrict __ATTR_ALIGN__(64) q_x,
							      float * __restrict __ATTR_ALIGN__(64) q_y,
							      float * __restrict __ATTR_ALIGN__(64) q_z,
							      float * __restrict __ATTR_ALIGN__(64) q_w,
							      const int32_t n);
	  
     }

}










#endif /*__GMS_Q4X16_TO_RMAT9X16_LOOPED_AVX512_H__*/
