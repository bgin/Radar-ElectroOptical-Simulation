
#ifndef __GMS_ROTATION_KERNELS_AVX512_H__
#define __GMS_ROTATION_KERNELS_AVX512_H__ 090820210225


namespace file_info {

const unsigned int gGMS_ROTATION_KERNELS_AVX512_MAJOR = 1U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_MINOR = 0U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_MICRO = 1U;
const unsigned int gGMS_ROTATION_KERNELS_AVX512_FULLVER =
       1000U*gGMS_ROTATION_KERNELS_AVX512_MAJOR+
       100U*gGMS_ROTATION_KERNELS_AVX512_MINOR +
       10U*gGMS_ROTATION_KERNELS_AVX512_MICRO;
const char * const pgGMS_ROTATION_KERNELS_AVX512_CREATION_DATE = "09-08-2021 02:25 PM +00200 (SUN 09 AUG 2021 GMT+2)";
const char * const pgGMS_ROTATION_KERNELS_AVX512_BUILD_DATE    = __DATE__ ":" __TIME__;
const char * const pgGMS_ROTATION_KERNELS_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
const char * const pgGMS_ROTATION_KERNELS_AVX512_DESCRIPTION   = "AVX512 vectorized basic rotation operations.";
}


#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_dcm_avx512.hpp"


namespace gms {

         namespace math {

                     


		  

			       

                               

				 


	             


	           


		   


		   


		      

	             
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      DCM9x16
		      q4x16_to_rmat9x16_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
						const __m512 q_z,
						const __m512 q_w); 


	             
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
                      DCM9x8
		      q4x8_to_rmat9x8_zmm8r8(const __m512d q_x,
		                             const __m512d q_y,
					     const __m512d q_z,
					     const __m512d q_w); 


	          /*
                       Random rotation matrix 3x3 of sphere.
                       Based on Jim Arvo, 1991 implementation.
                       Original Matrix 3x3 is represented as
                       SIMD Matrix 9x16
                   */

		    
                    __ATTR_VECTORCALL__
                    __ATTR_HOT__
		    DCM9x16
		    random_sphere_rm9x16_zmm16r4(const __m512 vr1,
		                                 const __m512 vr2,
                                                 const __m512 vr3); 


		      /*
                       Random rotation matrix 3x3 of sphere.
                       Based on Jim Arvo, 1991 implementation.
                       Original Matrix 3x3 is represented as
                       SIMD Matrix 9x8
                   */

		    
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
		    DCM9x8
		    random_sphere_rm9x8_zmm8r8(  const __m512d vr1,
		                                 const __m512d vr2,
                                                 const __m512d vr3); 


		    
                            /*  This algorithm generates a gaussian deviate for each coordinate, so
                                *  the total effect is to generate a symmetric 4-D gaussian distribution,
                                *  by separability. Projecting onto the surface of the hypersphere gives
                                *  a uniform distribution.
                                Based on  Ken Shoemake, September 1991 implementation.
                                Manually vectorized.
                            */

		       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      urand_q4x16_a_zmm16r4(const __m512 vrx,   // random gaussian vector uniformly distributed [0,1]
		                            const __m512 vry,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrz,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrw,   // random gaussian vector uniformly distributed [0,1]
		                            float * __restrict __ATTR_ALIGN__(64) q_x,
					    float * __restrict __ATTR_ALIGN__(64) q_y,
					    float * __restrict __ATTR_ALIGN__(64) q_z,
					    float * __restrict __ATTR_ALIGN__(64) q_w); 


		    
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      urand_q4x16_u_zmm16r4(const __m512 vrx,   // random gaussian vector uniformly distributed [0,1]
		                            const __m512 vry,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrz,   // random gaussian vector uniformly distributed [0,1]
					    const __m512 vrw,   // random gaussian vector uniformly distributed [0,1]
		                            float * __restrict q_x,
					    float * __restrict q_y,
					    float * __restrict q_z,
					    float * __restrict q_w); 



		    
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      urand_q4x8_a_zmm8r8(  const __m512d vrx,   // random gaussian vector uniformly distributed [0,1]
		                          const __m512d vry,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrz,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrw,   // random gaussian vector uniformly distributed [0,1]
		                          double * __restrict  __ATTR_ALIGN__(64) q_x,
					  double * __restrict  __ATTR_ALIGN__(64) q_y,
					  double * __restrict  __ATTR_ALIGN__(64) q_z,
					  double * __restrict  __ATTR_ALIGN__(64) q_w); 

		    
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      urand_q4x8_u_zmm8r8(  const __m512d vrx,   // random gaussian vector uniformly distributed [0,1]
		                          const __m512d vry,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrz,   // random gaussian vector uniformly distributed [0,1]
					  const __m512d vrw,   // random gaussian vector uniformly distributed [0,1]
		                          double * __restrict   q_x,
					  double * __restrict   q_y,
					  double * __restrict   q_z,
					  double * __restrict   q_w); 

	   

		    /*
                            Convert unit quaternion to Euler angles
                      */

		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      q4x16_to_ea3x16_a_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      float * __restrict __ATTR_ALIGN__(64) alpha,
					      float * __restrict __ATTR_ALIGN__(64) beta,
					      float * __restrict __ATTR_ALIGN__(64) gamma); 

		  
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      q4x16_to_ea3x16_u_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
					        const __m512 q_z,
					        const __m512 q_w,
					        float * __restrict  alpha,
					        float * __restrict  beta,
					        float * __restrict  gamma); 




		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      q4x8_to_ea3x8_a_zmm8r8(   const __m512d q_x,
		                                const __m512d q_y,
					        const __m512d q_z,
					        const __m512d q_w,
					        double * __restrict __ATTR_ALIGN__(64) alpha,
					        double * __restrict __ATTR_ALIGN__(64) beta,
					        double * __restrict __ATTR_ALIGN__(64) gamma); 
					        

		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
	              void
		      q4x8_to_ea3x8_u_zmm8r8(   const __m512d q_x,
		                                const __m512d q_y,
					        const __m512d q_z,
					        const __m512d q_w,
					        double * __restrict  alpha,
					        double * __restrict  beta,
					        double * __restrict  gamma); 




		   /*
                       Convert unit quaternion to axis angle pair
                    */


		    
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x16_to_ax4x16_a_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
					        const __m512 q_z,
					        const __m512 q_w,
					        float * __restrict __ATTR_ALIGN__(64) ax_1,
					        float * __restrict __ATTR_ALIGN__(64) ax_2,
					        float * __restrict __ATTR_ALIGN__(64) ax_3,
					        float * __restrict __ATTR_ALIGN__(64) ax_4); 


		      
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x16_to_ax4x16_u_zmm16r4(const __m512 q_x,
		                                const __m512 q_y,
					        const __m512 q_z,
					        const __m512 q_w,
					        float * __restrict  ax_1,
					        float * __restrict  ax_2,
					        float * __restrict  ax_3,
					        float * __restrict  ax_4); 





		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x8_to_ax4x8_a_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict __ATTR_ALIGN__(64) ax_1,
					      double * __restrict __ATTR_ALIGN__(64) ax_2,
					      double * __restrict __ATTR_ALIGN__(64) ax_3,
					      double * __restrict __ATTR_ALIGN__(64) ax_4); 


		      
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x8_to_ax4x8_u_zmm8r8( const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict  ax_1,
					      double * __restrict  ax_2,
					      double * __restrict  ax_3,
					      double * __restrict  ax_4); 



		    /*
                        Convert unit quaternion to Rodrigues vector
                     */
                      
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x16_to_rv4x16_a_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      float * __restrict __ATTR_ALIGN__(64) r_x,
					      float * __restrict __ATTR_ALIGN__(64) r_y,
					      float * __restrict __ATTR_ALIGN__(64) r_z,
					      float * __restrict __ATTR_ALIGN__(64) r_w); 
					      


		      
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x16_to_rv4x16_u_zmm16r4(const __m512 q_x,
		                              const __m512 q_y,
					      const __m512 q_z,
					      const __m512 q_w,
					      float * __restrict r_x,
					      float * __restrict r_y,
					      float * __restrict r_z,
					      float * __restrict r_w);

		     
	              
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x8_to_rv4x8_a_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict __ATTR_ALIGN__(64) r_x,
					      double * __restrict __ATTR_ALIGN__(64) r_y,
					      double * __restrict __ATTR_ALIGN__(64) r_z,
					      double * __restrict __ATTR_ALIGN__(64) r_w); 
					      

		     
                      __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      q4x8_to_rv4x8_u_zmm8r8(   const __m512d q_x,
		                              const __m512d q_y,
					      const __m512d q_z,
					      const __m512d q_w,
					      double * __restrict  r_x,
					      double * __restrict  r_y,
					      double * __restrict  r_z,
					      double * __restrict  r_w); 


		     /*
                           Orientation i.e. (Direct Cosine Matrix)  matrix to Euler angles.
                       */


		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      rmat9x16_to_ea3x16_a_zmm16r4(const DCM9x16 rm,
		                                 float * __restrict __ATTR_ALIGN__(64) alpha,
						 float * __restrict __ATTR_ALIGN__(64) beta,
						 float * __restrict __ATTR_ALIGN__(64) gamma);


		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      rmat9x16_to_ea3x16_u_zmm16r4(const DCM9x16 rm,
		                                 float * __restrict  alpha,
						 float * __restrict  beta,
						 float * __restrict  gamma);



		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      rmat9x8_to_ea3x8_a_zmm8r8(const DCM9x8 rm,
		                              double * __restrict __ATTR_ALIGN__(64) alpha,
					      double * __restrict __ATTR_ALIGN__(64) beta,
					      double * __restrict __ATTR_ALIGN__(64) gamma); 


		     
                       __ATTR_VECTORCALL__
                      __ATTR_HOT__
		      void
		      rmat9x8_to_ea3x8_u_zmm8r8(const DCM9x8 rm,
		                              double * __restrict  alpha,
					      double * __restrict  beta,
					      double * __restrict  gamma);
					
		      
     } // math

} // gms













#endif /* __GMS_ROTATION_KERNELS_AVX512_H__*/
