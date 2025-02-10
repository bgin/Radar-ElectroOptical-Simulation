

#ifndef __GMS_RAND_RMAT9X8_LOOPED_AVX512_H__
#define __GMS_RAND_RMAT9X8_LOOPED_AVX512_H__


namespace file_info {

   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_MAJOR = 1;
   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_MINOR = 0;
   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_MICRO = 0;
   const unsigned int GMS_RAND_RMAT9X8_LOOPED_AVX512_FULLVER =
    1000*GMS_RAND_RMAT9X8_LOOPED_AVX512_MAJOR+
    100*GMS_RAND_RMAT9X8_LOOPED_AVX512_MINOR+
    10*GMS_RAND_RMAT9X8_LOOPED_AVX512_MICRO;
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_CREATE_DATE = "28-11-2021 13357 +00200 (SUN 28 NOV 2021 GMT+2)";
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_BUILD_DATE  = __DATE__":"__TIME__;
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_AUTHOR      = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
   const char * const GMS_RAND_RMAT9X8_LOOPED_AVX512_DESCRIPTION = "Random Rotation Matrix 9x8 (of sphere) kernel AVX512.";
}



#include <cstdint>
#include "GMS_config.h"
.




namespace  gms {


        namespace  math {

	

                   
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     
                      void
                      rand_rm9x8_looped_u_nonunroll_zmm8r8(   double * __restrict row1,
		                                              double * __restrict row2,
							      double * __restrict row3,
							      double * __restrict row4,
							      double * __restrict row5,
							      double * __restrict  row6,
							      double * __restrict  row7,
							      double * __restrict  row8,
							      double * __restrict  row9,
							      const double * __restrict  rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  rz, // ramdomly normally distributed vector [0,1]
							      const int32_t n); 


		     
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		     
                      void
                      rand_rm9x8_looped_u_nonunroll_zmm8r8(   double * __restrict __ATTR_ALIGN__(64) row1,
		                                              double * __restrict __ATTR_ALIGN__(64) row2,
							      double * __restrict __ATTR_ALIGN__(64) row3,
							      double * __restrict __ATTR_ALIGN__(64) row4,
							      double * __restrict __ATTR_ALIGN__(64) row5,
							      double * __restrict __ATTR_ALIGN__(64) row6,
							      double * __restrict __ATTR_ALIGN__(64) row7,
							      double * __restrict __ATTR_ALIGN__(64) row8,
							      double * __restrict __ATTR_ALIGN__(64) row9,
							      const double * __restrict __ATTR_ALIGN__(64) rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict __ATTR_ALIGN__(64) ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict __ATTR_ALIGN__(64) rz, // ramdomly normally distributed vector [0,1]
							      const int32_t n); 
							      

		     
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		    
                      void
                      rand_rm9x8_looped_u_unroll4x_zmm8r8(    double * __restrict row1,
		                                              double * __restrict row2,
							      double * __restrict row3,
							      double * __restrict row4,
							      double * __restrict row5,
							      double * __restrict  row6,
							      double * __restrict  row7,
							      double * __restrict  row8,
							      double * __restrict  row9,
							      const double * __restrict  rx, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  ry, // ramdomly normally distributed vector [0,1]
							      const double * __restrict  rz, // ramdomly normally distributed vector [0,1]
							      const int32_t n); 


		     

         }






}








#endif /*__GMS_RAND_RMAT9X8_LOOPED_AVX512_H__*/
