

#ifndef __GMS_POS_TO_STATE_SSE_HPP__
#define __GMS_POS_TO_STATE_SSE_HPP__ 280820240727


/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
@@Modified by Bernard Gingold, on 22-05-2022 15:52 +00200 (SUN 22 MAY 2022 15:25 GMT+2)
  contact: beniekg@gmail.com
*/

namespace file_info {

 const unsigned int GMS_POS_TO_STATE_SSE_MAJOR = 1U;
 const unsigned int GMS_POS_TO_STATE_SSE_MINOR = 0U;
 const unsigned int GMS_POS_TO_STATE_SSE_MICRO = 0U;
 const unsigned int GMS_POS_TO_STATE_SSE_FULLVER =
  1000U*GMS_POS_TO_STATE_SSE_MAJOR+100U*GMS_POS_TO_STATE_SSE_MINOR+10U*GMS_POS_TO_STATE_SSE_MICRO;
 const char * const GMS_POS_TO_STATE_SSE_CREATION_DATE = "28-08-2024 07:27PM +00200 (WED 28 AUG 2024 07:27 GMT+2)";
 const char * const GMS_POS_TO_STATE_SSE_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const GMS_POS_TO_STATE_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const GMS_POS_TO_STATE_SSE_SYNOPSIS      = "SSE based position [2D] to state convertion functions (vectorized)."


}

#if !defined(__AVX512F__) || !defined(__AVX512VL__)
#error "Support of AVX512F or AVX512VL required!!"
#endif

#include <immintrin.h>
#include "GMS_config.h"

#include "GMS_simd_utils.hpp"

namespace gms {

        namespace math {
	

/*
https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary/tree/27317126c57c19f57f6d14eb6aa35700f57b869a/Coordinate_Systems/State_Conversion
 %%CART2DSSTATE2POLARSTATE Transform a 2D Cartesian state into a state
%                         consisting of position, heading and speed as well
%                         as possibly a turn rate and a linear
%                         acceleration, depending on the choice of
%                         systemType.
%
%INPUTS: xCart A Cartesian state vector consisting of position velocity and
%              possibly acceleration into a state where heading and speed
%              have been separated. xCart has the form
%              [x;y;xdot;ydot;xddot;yddot], where the acceleration terms
%              xddot;yddot can be omitted if the system type is 'ConstVel'.
%   systemType A string constant specifying the desired type of output. In
%              all instances, the heading is measured in terms of radians
%              counterclockwise from the x-axis. Possible values are:
%              'ConstVel'     The target state is [position;heading;speed]
%                             and xCart is [position;velocity]
%              'ConstAccel'   The target state is [position;heading;speed;
%                             speed derivative] and xCart is
%                             [position;velocity;acceleration]
%              'ConstTurn'    The target state is [position;heading;speed;
%                             turn rate] and xCart is
%                             [position;velocity;acceleration]
%              'TurnAndAccel' The target state is [position;heading;speed;
%                             turnrate; speed derivative] and xCart is
%                             [position;velocity;acceleration]
%
%OUTPUTS: xPol The state converted from 2D Cartesian coordinates into the
%              selected 2D coordinate system.
%
%When the system type is 'ConstVel' or 'TurnAndAccel', only a single
%solution is mathematically observable. When the system type is
%'ConstAccel' or 'ConstTurn', the system is overdetermined, but only a
%simple solution is used, not a least squares solution.
%
%The use of 2D states where the heading and speed have been separated is
%discussed in [1] and [2].
%
%The opposite of this function is polar2DState2CartState.
%
%REFERENCES:
%[1] M. Busch and S. Blackman, "Evaluation of IMM filtering for an air
%    defense system application," in Proceedings of SPIE: Signal and Data
%    Processing of Small Targets, vol. 2561, 9 Jul. 1995, pp. 435-447.
%[1] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
**@@Modified Bernard Gingold May 2022 ,beniekg@gmail.com
%%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
*/
	             
                      __ATTR_ALWAYS_INLINE__
		      static inline void
                      const_velocity_xmm2r8(const __m128d xDot,
					    const __m128d yDot,
					    __m128d * __restrict s_a,
					    __m128d * __restrict s_b) {

                         *s_a = _mm_atan2_pd(yDot,xDot);

                         *s_b = _mm_fmadd_pd(yDot,yDot,_mm_mul_pd(xDot,xDot));
		     }


		   
                      __ATTR_ALWAYS_INLINE__
		       static inline void
                      const_velocity_xmm2r8_a(const __m128d xDot,
					      const __m128d yDot,
					      double * __restrict __ATTR_ALIGN__(16)  s_a,
					      double * __restrict __ATTR_ALIGN__(16)  s_b) {

                         _mm_store_pd(&s_a[0],_mm_atan2_pd(yDot,xDot));

                         _mm_store_pd(&s_b[0],_mm_fmadd_pd(yDot,yDot,_mm_mul_pd(xDot,xDot)));
		     }


		   
                      __ATTR_ALWAYS_INLINE__
		      static inline void
                      const_velocity_xmm2r8_u(const __m128d xDot,
					      const __m128d yDot,
					      double * __restrict s_a,
					      double * __restrict s_b) {

                         _mm_storeu_pd(&s_a[0],_mm_atan2_pd(yDot,xDot));

                         _mm_storeu_pd(&s_b[0],_mm_fmadd_pd(yDot,yDot,_mm_mul_pd(xDot,xDot)));
		     }


		     
                      __ATTR_ALWAYS_INLINE__
		      static inline void
                      const_velocity_xmm4r4(const __m128 xDot,
					     const __m128 yDot,
					     __m128 * __restrict s_a,
					     __m128 * __restrict s_b) {

                         *s_a = _mm_atan2_ps(yDot,xDot);

                         *s_b = _mm_fmadd_ps(yDot,yDot,_mm_mul_ps(xDot,xDot));
		     }


		    
                      __ATTR_ALWAYS_INLINE__
		      static inline void
                      const_velocity_xmm4r4_a(const __m128 xDot,
					       const __m128 yDot,
					       float * __restrict __ATTR_ALIGN__(16)  s_a,
					       float * __restrict __ATTR_ALIGN__(16)  s_b) {

                         _mm_store_pd(&s_a[0],_mm_atan2_ps(yDot,xDot));

                         _mm_store_ps(&s_b[0],_mm_fmadd_ps(yDot,yDot,_mm_mul_ps(xDot,xDot)));
		     }


		    
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_velocity_xmm4r4_u(const __m128 xDot,
					       const __m128 yDot,
					       float * __restrict  s_a,
					       float * __restrict  s_b) {

                         _mm_storeu_pd(&s_a[0],_mm_atan2_ps(yDot,xDot));

                         _mm_storeu_ps(&s_b[0],_mm_fmadd_ps(yDot,yDot,_mm_mul_ps(xDot,xDot)));
		     }



		   
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
                      const_acceleration_xmm2r8(const __m128d xDot,
					        const __m128d yDot,
		                                const __m128d xDdot,
		                                const __m128d yDdot,
						__m128d * __restrict s_a,
						__m128d * __restrict s_b,
						__m128d * __restrict s_c) {
       
                         __m128d theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;

                         theta = _mm_atan2_pd(yDot,xDot);
			 *s_a  = theta;
			 costh = _mm_cos_pd(theta);
			 sinth = _mm_sin_pd(theta);

                         *s_b  = _mm_fmadd_pd(yDot,yDot,_mm_mul_pd(xDot,xDot));
                         vDot  = _mm_sqrt_pd(_mm_fmadd_pd(yDdot,yDdot,
			                                    _mm_mul_pd(xDdot,xDdot)));
			 t0    = _mm_sub_pd(_mm_mul_pd(vDot,costh),xDdot);
			 t1    = _mm_sub_pd(_mm_mul_pd(vDot,sinth),yDdot);
			 t2    = _mm_sub_pd(_mm_mul_pd(xmm2r8_negate(vDot),costh),xDdot);
			 t3    = _mm_sub_pd(_mm_mul_pd(xmm2r8_negate(vDot),sinth),yDdot);
			 diff1 = _mm_fmadd_pd(t0,t0,_mm_mul_pd(t1,t1));
			 diff2 = _mm_fmadd_pd(t2,t2,_mm_mul_pd(t3,t3));
			 m     = _mm_cmp_pd_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm_mask_blend_pd(m,vDot,xmm2r8_negate(vDot));
			 *s_c  = vdot;
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm2r8_a(const __m128d xDot,
					          const __m128d yDot,
		                                  const __m128d xDdot,
		                                  const __m128d yDdot,
						  double * __restrict __ATTR_ALIGN__(16)  s_a,
						  double * __restrict __ATTR_ALIGN__(16)  s_b,
						  double * __restrict __ATTR_ALIGN__(16)  s_c) {
       
                         __m128d theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;

                         theta = _mm_atan2_pd(yDot,xDot);
			 _mm_store_pd(&s_a[0],theta);
			 costh = _mm_cos_pd(theta);
			 sinth = _mm_sin_pd(theta);
                         _mm_store_pd(&s_b[0],_mm_fmadd_pd(yDot,yDot,_mm_mul_pd(xDot,xDot)));
                         vDot  = _mm_sqrt_pd(_mm_fmadd_pd(yDdot,yDdot,
			                                    _mm_mul_pd(xDdot,xDdot)));
			 t0    = _mm_sub_pd(_mm_mul_pd(vDot,costh),xDdot);
			 t1    = _mm_sub_pd(_mm_mul_pd(vDot,sinth),yDdot);
			 t2    = _mm_sub_pd(_mm_mul_pd(xmm2r8_negate(vDot),costh),xDdot);
			 t3    = _mm_sub_pd(_mm_mul_pd(xmm2r8_negate(vDot),sinth),yDdot);
			 diff1 = _mm_fmadd_pd(t0,t0,_mm_mul_pd(t1,t1));
			 diff2 = _mm_fmadd_pd(t2,t2,_mm_mul_pd(t3,t3));
			 m     = _mm_cmp_pd_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm_mask_blend_pd(m,vDot,xmm2r8_negate(vDot));
			 _mm_store_pd(&s_c[0],vdot);
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
                      const_acceleration_xmm2r8_u(const __m128d xDot,
					          const __m128d yDot,
		                                  const __m128d xDdot,
		                                  const __m128d yDdot,
						  double * __restrict s_a,
						  double * __restrict s_b,
						  double * __restrict s_c) {
       
                         __m128d theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;

                         theta = _mm_atan2_pd(yDot,xDot);
			 _mm_storeu_pd(&s_a[0],theta);
			 costh = _mm_cos_pd(theta);
			 sinth = _mm_sin_pd(theta);

                         _mm_storeu_pd(&s_b[0],_mm_fmadd_pd(yDot,yDot,_mm_mul_pd(xDot,xDot)));
                         vDot  = _mm_sqrt_pd(_mm_fmadd_pd(yDdot,yDdot,
			                                    _mm_mul_pd(xDdot,xDdot)));
			 t0    = _mm_sub_pd(_mm_mul_pd(vDot,costh),xDdot);
			 t1    = _mm_sub_pd(_mm_mul_pd(vDot,sinth),yDdot);
			 t2    = _mm_sub_pd(_mm_mul_pd(xmm2r8_negate(vDot),costh),xDdot);
			 t3    = _mm_sub_pd(_mm_mul_pd(xmm2r8_negate(vDot),sinth),yDdot);
			 diff1 = _mm_fmadd_pd(t0,t0,_mm_mul_pd(t1,t1));
			 diff2 = _mm_fmadd_pd(t2,t2,_mm_mul_pd(t3,t3));
			 m     = _mm_cmp_pd_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm_mask_blend_pd(m,vDot,xmm2r8_negate(vDot));
			 _mm_storeu_pd(&s_c[0],vdot);
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm4r4(const __m128 xDot,
					         const __m128 yDot,
		                                 const __m128 xDdot,
		                                 const __m128 yDdot,
						 __m128 * __restrict s_a,
						 __m128 * __restrict s_b,
						 __m128 * __restrict s_c) {
       
                         __m128 theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;

                         theta = _mm_atan2_ps(yDot,xDot);
			 *s_a  = theta;
			 costh = _mm_cos_ps(theta);
			 sinth = _mm_sin_ps(theta);

                         *s_b  = _mm_fmadd_ps(yDot,yDot,_mm_mul_ps(xDot,xDot));
                         vDot  = _mm_sqrt_ps(_mm_fmadd_ps(yDdot,yDdot,
			                                    _mm_mul_ps(xDdot,xDdot)));
			 t0    = _mm_sub_ps(_mm_mul_ps(vDot,costh),xDdot);
			 t1    = _mm_sub_ps(_mm_mul_ps(vDot,sinth),yDdot);
			 t2    = _mm_sub_ps(_mm_mul_ps(xmm4r4_negate(vDot),costh),xDdot);
			 t3    = _mm_sub_ps(_mm_mul_ps(xmm4r4_negate(vDot),sinth),yDdot);
			 diff1 = _mm_fmadd_ps(t0,t0,_mm_mul_ps(t1,t1));
			 diff2 = _mm_fmadd_ps(t2,t2,_mm_mul_ps(t3,t3));
			 m     = _mm_cmp_ps_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm_mask_blend_ps(m,vDot,xmm4r4_negate(vDot));
			 *s_c  = vdot;
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
                      const_acceleration_xmm4r4_a(const __m128 xDot,
					           const __m128 yDot,
		                                   const __m128 xDdot,
		                                   const __m128 yDdot,
						   float * __restrict __ATTR_ALIGN__(16)  s_a,
						   float * __restrict __ATTR_ALIGN__(16)  s_b,
						   float * __restrict __ATTR_ALIGN__(16)  s_c) {
       
                         __m128 theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;

                         theta = _mm_atan2_ps(yDot,xDot);
			 _mm_store_ps(&s_a[0],theta);
			 costh = _mm_cos_ps(theta);
			 sinth = _mm_sin_ps(theta);

                         _mm_store_ps(&s_b[0],_mm_fmadd_ps(yDot,yDot,_mm_mul_ps(xDot,xDot)));
                         vDot  = _mm_sqrt_ps(_mm_fmadd_ps(yDdot,yDdot,
			                                    _mm_mul_ps(xDdot,xDdot)));
			 t0    = _mm_sub_ps(_mm_mul_ps(vDot,costh),xDdot);
			 t1    = _mm_sub_ps(_mm_mul_ps(vDot,sinth),yDdot);
			 t2    = _mm_sub_ps(_mm_mul_ps(xmm4r4_negate(vDot),costh),xDdot);
			 t3    = _mm_sub_ps(_mm_mul_ps(xmm4r4_negate(vDot),sinth),yDdot);
			 diff1 = _mm_fmadd_ps(t0,t0,_mm_mul_ps(t1,t1));
			 diff2 = _mm_fmadd_ps(t2,t2,_mm_mul_ps(t3,t3));
			 m     = _mm_cmp_ps_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm_mask_blend_ps(m,vDot,xmm4r4_negate(vDot));
			 _mm_store_ps(&s_c[0],vdot);
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm4r4_u(const __m128 xDot,
					           const __m128 yDot,
		                                   const __m128 xDdot,
		                                   const __m128 yDdot,
						   float * __restrict __ATTR_ALIGN__(16)  s_a,
						   float * __restrict __ATTR_ALIGN__(16)  s_b,
						   float * __restrict __ATTR_ALIGN__(16)  s_c) {
       
                         __m128 theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;

                         theta = _mm_atan2_ps(yDot,xDot);
			 _mm_storeu_ps(&s_a[0],theta);
			 costh = _mm_cos_ps(theta);
			 sinth = _mm_sin_ps(theta);

                         _mm_storeu_ps(&s_b[0],_mm_fmadd_ps(yDot,yDot,_mm_mul_ps(xDot,xDot)));
                         vDot  = _mm_sqrt_ps(_mm_fmadd_ps(yDdot,yDdot,
			                                    _mm_mul_ps(xDdot,xDdot)));
			 t0    = _mm_sub_ps(_mm_mul_ps(vDot,costh),xDdot);
			 t1    = _mm_sub_ps(_mm_mul_ps(vDot,sinth),yDdot);
			 t2    = _mm_sub_ps(_mm_mul_ps(xmm4r4_negate(vDot),costh),xDdot);
			 t3    = _mm_sub_ps(_mm_mul_ps(xmm4r4_negate(vDot),sinth),yDdot);
			 diff1 = _mm_fmadd_ps(t0,t0,_mm_mul_ps(t1,t1));
			 diff2 = _mm_fmadd_ps(t2,t2,_mm_mul_ps(t3,t3));
			 m     = _mm_cmp_ps_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm_mask_blend_ps(m,vDot,xmm4r4_negate(vDot));
			 _mm_storeu_ps(&s_c[0],vdot);
		    }




		      
                      __ATTR_ALWAYS_INLINE__
		    
		      static inline void
		      const_turn_xmm2r8(const __m128d xDot,
					const __m128d yDot,
		                        const __m128d xDdot,
		                        const __m128d yDdot,
					__m128d * __restrict __ATTR_ALIGN__(16)  a,
					__m128d * __restrict __ATTR_ALIGN__(16)  s,
					__m128d * __restrict __ATTR_ALIGN__(16)  omega) {

                          const __m128d t0 = _mm_fmsub_pd(xDot,yDot,
			                                 _mm_mul_pd(yDot,xDot));
			  const __m128d t1 = _mm_fmadd_pd(xDot,xDot,
			                                 _mm_mul_pd(yDot,yDot));
			  *omega           = _mm_div_pd(t0,t1);

                          *a               = _mm_atan2_pd(yDot,xDot);

                          *s               = _mm_sqrt_pd(_mm_fmadd_pd(yDot,yDot,
			                                                _mm_mul_pd(xDot,xDot)));
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      const_turn_xmm2r8_a(const __m128d xDot,
					  const __m128d yDot,
		                          const __m128d xDdot,
		                          const __m128d yDdot,
					  double * __restrict __ATTR_ALIGN__(16)  a,
					  double * __restrict __ATTR_ALIGN__(16)  s,
					  double * __restrict __ATTR_ALIGN__(16)  omega) {

                          const __m128d t0 = _mm_fmsub_pd(xDot,yDot,
			                                 _mm_mul_pd(yDot,xDot));
			  const __m128d t1 = _mm_fmadd_pd(xDot,xDot,
			                                 _mm_mul_pd(yDot,yDot));
			  _mm_store_pd(&omega[0],_mm_div_pd(t0,t1));

                          _mm_store_pd(&a[0],_mm_atan2_pd(yDot,xDot));

                          _mm_store_pd(&s[0],_mm_sqrt_pd(_mm_fmadd_pd(yDot,yDot,
			                                                _mm_mul_pd(xDot,xDot))));
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      const_turn_xmm2r8_u(const __m128d xDot,
					  const __m128d yDot,
		                          const __m128d xDdot,
		                          const __m128d yDdot,
					  double * __restrict  a,
					  double * __restrict  s,
					  double * __restrict  omega) {

                          const __m128d t0 = _mm_fmsub_pd(xDot,yDot,
			                                 _mm_mul_pd(yDot,xDot));
			  const __m128d t1 = _mm_fmadd_pd(xDot,xDot,
			                                 _mm_mul_pd(yDot,yDot));
			  _mm_storeu_pd(&omega[0],_mm_div_pd(t0,t1));

                          _mm_storeu_pd(&a[0],_mm_atan2_pd(yDot,xDot));

                          _mm_storeu_pd(&s[0],_mm_sqrt_pd(_mm_fmadd_pd(yDot,yDot,
			                                                _mm_mul_pd(xDot,xDot))));
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      const_turn_xmm4r4(const __m128 xDot,
					 const __m128 yDot,
		                         const __m128 xDdot,
		                         const __m128 yDdot,
					__m128 * __restrict __ATTR_ALIGN__(16)  a,
					__m128 * __restrict __ATTR_ALIGN__(16)  s,
					__m128 * __restrict __ATTR_ALIGN__(16)  omega) {

                          const __m128 t0 = _mm_fmsub_ps(xDot,yDdot,
			                                 _mm_mul_ps(yDot,xDdot));
			  const __m128 t1 = _mm_fmadd_ps(xDot,xDot,
			                                 _mm_mul_ps(yDot,yDot));
			  *omega           = _mm_div_ps(t0,t1);

                          *a               = _mm_atan2_ps(yDot,xDot);

                          *s               = _mm_sqrt_ps(_mm_fmadd_ps(yDot,yDot,
			                                                _mm_mul_ps(xDot,xDot)));
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      const_turn_xmm4r4_a(const __m128 xDot,
					   const __m128 yDot,
		                           const __m128 xDdot,
		                           const __m128 yDdot,
					   float * __restrict __ATTR_ALIGN__(16)  a,
					   float * __restrict __ATTR_ALIGN__(16)  s,
					   float * __restrict __ATTR_ALIGN__(16)  omega) {

                          const __m128 t0 = _mm_fmsub_ps(xDot,yDdot,
			                                 _mm_mul_ps(yDot,xDdot));
			  const __m128 t1 = _mm_fmadd_ps(xDot,xDot,
			                                 _mm_mul_ps(yDot,yDot));
			  _mm_store_ps(&omega[0],_mm_div_ps(t0,t1));

                          _mm_store_ps(&a[0],_mm_atan2_ps(yDot,xDot));

                          _mm_store_ps(&s[0],_mm_sqrt_ps(_mm_fmadd_ps(yDot,yDot,
			                                                _mm_mul_ps(xDot,xDot))));
		    }


		    
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      const_turn_xmm4r4_a(const __m128 xDot,
					   const __m128 yDot,
		                           const __m128 xDdot,
		                           const __m128 yDdot,
					   float * __restrict  a,
					   float * __restrict  s,
					   float * __restrict  omega) {

                          const __m128 t0 = _mm_fmsub_ps(xDot,yDdot,
			                                 _mm_mul_ps(yDot,xDdot));
			  const __m128 t1 = _mm_fmadd_ps(xDot,xDot,
			                                 _mm_mul_ps(yDot,yDot));
			  _mm_storeu_ps(&omega[0],_mm_div_ps(t0,t1));

                          _mm_storeu_ps(&a[0],_mm_atan2_ps(yDot,xDot));

                          _mm_storeu_ps(&s[0],_mm_sqrt_ps(_mm_fmadd_ps(yDot,yDot,
			                                                _mm_mul_ps(xDot,xDot))));
		    }




	            
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      turn_accelerate_xmm2r8(const __m128d xDot,
					     const __m128d yDot,
		                             const __m128d xDdot,
		                             const __m128d yDdot,
					     __m128d * __restrict __ATTR_ALIGN__(16)  theta,
					     __m128d * __restrict __ATTR_ALIGN__(16)  v,
					     __m128d * __restrict __ATTR_ALIGN__(16)  omega,
					     __m128d * __restrict __ATTR_ALIGN__(16)  vDot) {

			

                         const __m128d th    = _mm_atan2_pd(yDot,xDot);
                         const __m128d costh = _mm_cos_pd(th);
			 const __m128d sinth = _mm_sin_pd(th);

                         *theta              = th;
			 *v                  = _mm_sqrt_pd(_mm_fmadd_pd(yDot,yDot,
			                                               _mm_mul_pd(xDot,xDot)));
			 *omega              = _mm_div_pd(_mm_fmsub_pd(yDot,costh,
			                                               _mm_mul_pd(xDot,sinth)),*v);
			 *vDot               = _mm_fmadd_pd(xDot,costh,
			                                               _mm_mul_pd(yDot,sinth));
                         
		    }


		      
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      turn_accelerate_xmm2r8_a(const __m128d xDot,
					       const __m128d yDot,
		                               const __m128d xDdot,
		                               const __m128d yDdot,
					       double * __restrict __ATTR_ALIGN__(16)  theta,
					       double * __restrict __ATTR_ALIGN__(16)  v,
					       double * __restrict __ATTR_ALIGN__(16)  omega,
					       double * __restrict __ATTR_ALIGN__(16)  vDot) {

			

                         const __m128d th    = _mm_atan2_pd(yDot,xDot);
                         const __m128d costh = _mm_cos_pd(th);
			 const __m128d sinth = _mm_sin_pd(th);

                         _mm_store_pd(&theta[0],th);
			 _mm_store_pd(&v[0],_mm_sqrt_pd(_mm_fmadd_pd(yDot,yDot,
			                                               _mm_mul_pd(xDot,xDot))));
			 _mm_store_pd(&omega[0],_mm_div_pd(_mm_fmsub_pd(yDot,costh,
			                                               _mm_mul_pd(xDot,sinth)),
								              _mm_load_pd(&v[0])));
			 _mm_store_pd(&vDot[0],_mm_fmadd_pd(xDot,costh,
			                                               _mm_mul_pd(yDot,sinth)));
                         
		    }


		    
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      turn_accelerate_xmm2r8_u(const __m128d xDot,
					       const __m128d yDot,
		                               const __m128d xDdot,
		                               const __m128d yDdot,
					       double * __restrict theta,
					       double * __restrict v,
					       double * __restrict omega,
					       double * __restrict vDot) {

			

                         const __m128d th    = _mm_atan2_pd(yDot,xDot);
                         const __m128d costh = _mm_cos_pd(th);
			 const __m128d sinth = _mm_sin_pd(th);

                         _mm_storeu_pd(&theta[0],th);
			 _mm_storeu_pd(&v[0],_mm_sqrt_pd(_mm_fmadd_pd(yDot,yDot,
			                                               _mm_mul_pd(xDot,xDot))));
			 _mm_storeu_pd(&omega[0],_mm_div_pd(_mm_fmsub_pd(yDot,costh,
			                                               _mm_mul_pd(xDot,sinth)),
								              _mm_loadu_pd(&v[0])));
			 _mm_storeu_pd(&vDot[0],_mm_fmadd_pd(xDot,costh,
			                                               _mm_mul_pd(yDot,sinth)));
                         
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      turn_accelerate_xmm4r4(const __m128 xDot,
					      const __m128 yDot,
		                              const __m128 xDdot,
		                              const __m128 yDdot,
					      __m128 * __restrict __ATTR_ALIGN__(16)  theta,
					      __m128 * __restrict __ATTR_ALIGN__(16)  v,
					      __m128 * __restrict __ATTR_ALIGN__(16)  omega,
					      __m128 * __restrict __ATTR_ALIGN__(16)  vDot) {

			

                         const __m128 th    = _mm_atan2_ps(yDot,xDot);
                         const __m128 costh = _mm_cos_ps(th);
			 const __m128 sinth = _mm_sin_ps(th);

                         *theta              = th;
			 *v                  = _mm_sqrt_ps(_mm_fmadd_ps(yDot,yDot,
			                                               _mm_mul_ps(xDot,xDot)));
			 *omega              = _mm_div_ps(_mm_fmsub_ps(yDot,costh,
			                                               _mm_mul_ps(xDot,sinth)),*v);
			 *vDot               = _mm_fmadd_ps(xDot,costh,
			                                               _mm_mul_ps(yDot,sinth));
                         
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      turn_accelerate_xmm4r4_a(const __m128 xDot,
					        const __m128 yDot,
		                                const __m128 xDdot,
		                                const __m128 yDdot,
					        float * __restrict __ATTR_ALIGN__(16)  theta,
					        float * __restrict __ATTR_ALIGN__(16)  v,
					        float * __restrict __ATTR_ALIGN__(16)  omega,
					        float * __restrict __ATTR_ALIGN__(16)  vDot) {

			

                         const __m128 th    = _mm_atan2_ps(yDot,xDot);
                         const __m128 costh = _mm_cos_ps(th);
			 const __m128 sinth = _mm_sin_ps(th);

                         _mm_store_ps(&theta[0],th);
			 _mm_store_ps(&v[0],_mm_sqrt_ps(_mm_fmadd_ps(yDot,yDot,
			                                               _mm_mul_ps(xDot,xDot))));
			 _mm_store_ps(&omega[0],_mm_div_ps(_mm_fmsub_ps(yDot,costh,
			                                               _mm_mul_ps(xDot,sinth)),
								              _mm_load_ps(&v[0])));
			 _mm_store_ps(&vDot[0],_mm_fmadd_ps(xDot,costh,
			                                               _mm_mul_ps(yDot,sinth)));
                         
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      turn_accelerate_xmm4r4_u(const __m128 xDot,
					        const __m128 yDot,
		                                const __m128 xDdot,
		                                const __m128 yDdot,
					        float * __restrict  theta,
					        float * __restrict  v,
					        float * __restrict  omega,
					        float * __restrict  vDot) {

			

                         const __m128 th    = _mm_atan2_ps(yDot,xDot);
                         const __m128 costh = _mm_cos_ps(th);
			 const __m128 sinth = _mm_sin_ps(th);

                         _mm_storeu_ps(&theta[0],th);
			 _mm_storeu_ps(&v[0],_mm_sqrt_ps(_mm_fmadd_ps(yDot,yDot,
			                                               _mm_mul_ps(xDot,xDot))));
			 _mm_storeu_ps(&omega[0],_mm_div_ps(_mm_fmsub_ps(yDot,costh,
			                                               _mm_mul_ps(xDot,sinth)),
								              _mm_loadu_ps(&v[0])));
			 _mm_storeu_ps(&vDot[0],_mm_fmadd_ps(xDot,costh,
			                                               _mm_mul_ps(yDot,sinth)));
                         
		    }


/*
  %%POLAR2DSTATE2CARTSTATE Convert a 2D target state where the velocity had
%               been decomposed into a direction angle (heading) and speed
%               components into Cartesian components. Depending on the
%               system type chosen, the state can have components for a
%               linear acceleration and/ or a turn rate.
%
%INPUTS: xPol  A 4X1, 5X1 or 6X1 target state where the first four
%              components are [position;heading;speed] in 2D. The other
%              components depend on the value of systemType.
%   systemType A string constant specifying the desired type of input and
%              output. In all instances, the heading is measured in terms
%              of radians counterclockwise from the x-axis. Possible values
%              are:
%              'ConstVel'     The target state is [position;heading;speed]
%                             and xCart is [position;velocity]
%              'ConstAccel'   The target state is [position;heading;speed;
%                             speed derivative] and xCart is
%                             [position;velocity;acceleration]
%              'ConstTurn'    The target state is [position;heading;speed;
%                             turn rate] and xCart is
%                             [position;velocity;acceleration]
%              'TurnAndAccel' The target state is [position;heading;speed;
%                             turnrate; speed derivative] and xCart is
%                             [position;velocity;acceleration]
%
%%OUTPUTS: xCart The state converted into 2D Cartesian coordinates
%                consisting of position and velocity and, depending on
%                systemType, possibly acceleration components.
%
%The use of 2D states where the heading and speed have been separated is
%discussed in [1] and [2].
%
%The opposite of this function is Cart2DState2PolarState.
%
%REFERENCES:
%[1] M. Busch and S. Blackman, "Evaluation of IMM filtering for an air
%    defense system application," in Proceedings of SPIE: Signal and Data
%    Processing of Small Targets, vol. 2561, 9 Jul. 1995, pp. 435-447.
%[2] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
@@Modified by Bernard Gingold, on 29-05-2022 09:12 +00200 (SUN 29 MAY 2022 09:12 GMT+2)
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
*/	           


	             
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void			     
		      const_velocity_xmm2r8(const __m128d theta,// heading
		                            const __m128d v,    // speed
                                            __m128d &vcth,
					    __m128d &vsth) {


                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         vcth = _mm_mul_pd(v,cth);
			 vsth = _mm_mul_pd(v,sth);
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void			     
		      const_velocity_xmm2r8_a(const __m128d theta,// heading
		                              const __m128d v,    // speed
                                              double * __restrict __ATTR_ALIGN__(16)  vcth,
					      double * __restrict __ATTR_ALIGN__(16)  vsth) {


                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         _mm_store_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_store_pd(&vsth[0],_mm_mul_pd(v,sth));
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void			     
		      const_velocity_xmm2r8_u(const __m128d theta,// heading
		                              const __m128d v,    // speed
                                              double * __restrict  vcth,
					      double * __restrict  vsth) {


                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         _mm_storeu_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_storeu_pd(&vsth[0],_mm_mul_pd(v,sth));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void			     
		      const_velocity_xmm4r4(const __m128 theta,// heading
		                             const __m128 v,    // speed
                                             __m128 &vcth,
					     __m128 &vsth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         vcth = _mm_mul_ps(v,cth);
			 vsth = _mm_mul_ps(v,sth);
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void			     
		      const_velocity_xmm4r4_a(const __m128 theta,// heading
		                               const __m128 v,    // speed
                                               float * __restrict __ATTR_ALIGN__(16)  vcth,
					       float * __restrict __ATTR_ALIGN__(16)  vsth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         _mm_store_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_store_ps(&vsth[0],_mm_mul_ps(v,sth));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void			     
		      const_velocity_xmm4r4_u(const __m128 theta,// heading
		                               const __m128 v,    // speed
                                               float * __restrict vcth,
					       float * __restrict vsth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         _mm_storeu_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_storeu_ps(&vsth[0],_mm_mul_ps(v,sth));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm2r8(const __m128d theta,
		                                const __m128d v,
						const __m128d vDot, //linear acceleration
						__m128d &vcth,
						__m128d &vsth,
						__m128d &vdcth,
						__m128d &vdsth) {


                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         vcth              = _mm_mul_pd(v,cth);
			 vsth              = _mm_mul_pd(v,sth);
			 vdcth             = _mm_mul_pd(vDot,cth);
			 vdsth             = _mm_mul_pd(vDot,sth);
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm2r8_a(const __m128d theta,
		                                  const __m128d v,
						  const __m128d vDot, //linear acceleration
						  double * __restrict __ATTR_ALIGN__(16)  vcth,
						  double * __restrict __ATTR_ALIGN__(16)  vsth,
						  double * __restrict __ATTR_ALIGN__(16)  vdcth,
						  double * __restrict __ATTR_ALIGN__(16)  vdsth) {


                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         _mm_store_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_store_pd(&vsth[0],_mm_mul_pd(v,sth));
			 _mm_store_pd(&vdcth[0],_mm_mul_pd(vDot,cth));
			 _mm_store_pd(&vdsth[0],_mm_mul_pd(vDot,sth));
		   }


		    
                      __ATTR_ALWAYS_INLINE__
		    
		      static inline void
                      const_acceleration_xmm2r8_u(const __m128d theta,
		                                  const __m128d v,
						  const __m128d vDot, //linear acceleration
						  double * __restrict  vcth,
						  double * __restrict  vsth,
						  double * __restrict  vdcth,
						  double * __restrict  vdsth) {


                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         _mm_storeu_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_storeu_pd(&vsth[0],_mm_mul_pd(v,sth));
			 _mm_storeu_pd(&vdcth[0],_mm_mul_pd(vDot,cth));
			 _mm_storeu_pd(&vdsth[0],_mm_mul_pd(vDot,sth));
		   }


		      
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm4r4(const __m128 theta,
		                                 const __m128 v,
						 const __m128 vDot, //linear acceleration
						__m128 &vcth,
						__m128 &vsth,
						__m128 &vdcth,
						__m128 &vdsth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         vcth              = _mm_mul_ps(v,cth);
			 vsth              = _mm_mul_ps(v,sth);
			 vdcth             = _mm_mul_ps(vDot,cth);
			 vdsth             = _mm_mul_ps(vDot,sth);
		   }


		      
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      const_acceleration_xmm4r4_a(const __m128 theta,
		                                   const __m128 v,
						   const __m128 vDot, //linear acceleration
						   float * __restrict __ATTR_ALIGN__(16)  vcth,
						   float * __restrict __ATTR_ALIGN__(16)  vsth,
						   float * __restrict __ATTR_ALIGN__(16)  vdcth,
						   float * __restrict __ATTR_ALIGN__(16)  vdsth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         _mm_store_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_store_ps(&vsth[0],_mm_mul_ps(v,sth));
			 _mm_store_ps(&vdcth[0],_mm_mul_ps(vDot,cth));
			 _mm_store_ps(&vdsth[0],_mm_mul_ps(vDot,sth));
		   }



		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
                      const_acceleration_xmm4r4_u(const __m128 theta,
		                                   const __m128 v,
						   const __m128 vDot, //linear acceleration
						   float * __restrict vcth,
						   float * __restric  vsth,
						   float * __restrict vdcth,
						   float * __restrict vdsth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         _mm_storeu_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_storeu_ps(&vsth[0],_mm_mul_ps(v,sth));
			 _mm_storeu_ps(&vdcth[0],_mm_mul_ps(vDot,cth));
			 _mm_storeu_ps(&vdsth[0],_mm_mul_ps(vDot,sth));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      const_turn_xmm2r8(const __m128d theta,
		                        const __m128d v,
					const __m128d omega, // turn rate
                                        __m128d &vcth,
					__m128d &vsth,
					__m128d &vomsth,
					__m128d &vomcth) {

                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         vcth              = _mm_mul_pd(v,cth);
			 vsth              = _mm_mul_pd(v,sth);
                         vomsth            = _mm_mul_pd(xmm2r8_negate(v),
			                               _mm_mul_pd(omega,sth));
			 vomcth            = _mm_mul_pd(v,_mm_mul_pd(omega,cth));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		    
		      static inline void
		      const_turn_xmm2r8_a(const __m128d theta,
		                          const __m128d v,
					  const __m128d omega, // turn rate
                                          double * __restrict __ATTR_ALIGN__(16)  vcth,
					  double * __restrict __ATTR_ALIGN__(16)  vsth,
					  double * __restrict __ATTR_ALIGN__(16)  vomsth,
					  double * __restrict __ATTR_ALIGN__(16)  vomcth) {

                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         _mm_store_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_store_pd(&vsth[0],_mm_mul_pd(v,sth));
                         _mm_store_pd(&vomsth[0],_mm_mul_pd(xmm2r8_negate(v),
			                               _mm_mul_pd(omega,sth)));
			 _mm_store_pd(&vomcth[0],_mm_mul_pd(v,_mm_mul_pd(omega,cth)));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      const_turn_xmm2r8_u(const __m128d theta,
		                          const __m128d v,
					  const __m128d omega, // turn rate
                                          double * __restrict  vcth,
					  double * __restrict  vsth,
					  double * __restrict  vomsth,
					  double * __restrict  vomcth) {

                         const __m128d cth = _mm_cos_pd(theta);
			 const __m128d sth = _mm_sin_pd(theta);

                         _mm_storeu_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_storeu_pd(&vsth[0],_mm_mul_pd(v,sth));
                         _mm_storeu_pd(&vomsth[0],_mm_mul_pd(xmm2r8_negate(v),
			                               _mm_mul_pd(omega,sth)));
			 _mm_storeu_pd(&vomcth[0],_mm_mul_pd(v,_mm_mul_pd(omega,cth)));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      const_turn_xmm4r4(const __m128 theta,
		                         const __m128 v,
					 const __m128 omega, // turn rate
                                         __m128 &vcth,
					 __m128 &vsth,
					 __m128 &vomsth,
					 __m128 &vomcth) {

                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         vcth              = _mm_mul_ps(v,cth);
			 vsth              = _mm_mul_ps(v,sth);
                         vomsth            = _mm_mul_ps(xmm4r4_negate(v),
			                               _mm_mul_ps(omega,sth));
			 vomcth            = _mm_mul_ps(v,_mm_mul_ps(omega,cth));
		   }



		    
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
		      const_turn_xmm4r4_a(const __m128 theta,
		                         const __m128 v,
					 const __m128 omega, // turn rate
                                         float * __restrict __ATTR_ALIGN__(16)  vcth,
					 float * __restrict __ATTR_ALIGN__(16)  vsth,
					 float * __restrict __ATTR_ALIGN__(16)  vomsth,
					 float * __restrict __ATTR_ALIGN__(16)  vomcth) {

                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         _mm_store_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_store_ps(&vsth[0],_mm_mul_ps(v,sth));
                         _mm_store_ps(&vomsth[0],_mm_mul_ps(xmm4r4_negate(v),
			                               _mm_mul_ps(omega,sth)));
			 _mm_store_ps(&vomcth[0],_mm_mul_ps(v,_mm_mul_ps(omega,cth)));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
		      const_turn_xmm4r4_u(const __m128 theta,
		                         const __m128 v,
					 const __m128 omega, // turn rate
                                         float * __restrict  vcth,
					 float * __restrict  vsth,
					 float * __restrict  vomsth,
					 float * __restrict  vomcth) {

                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         _mm_storeu_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_storeu_ps(&vsth[0],_mm_mul_ps(v,sth));
                         _mm_storeu_ps(&vomsth[0],_mm_mul_ps(xmm4r4_negate(v),
			                               _mm_mul_ps(omega,sth)));
			 _mm_storeu_ps(&vomcth[0],_mm_mul_ps(v,_mm_mul_ps(omega,cth)));
		   }


		     
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
                      turn_accelerate_xmm2r8(const __m128d theta,
		                             const __m128d omega,
		                             const __m128d vDot,
					     const __m128d v,
					     __m128d &vcth,
					     __m128d &vsth,
					     __m128d &vomsth,
					     __m128d &vomcth) {


                         const __m128 cth = _mm_cos_pd(theta);
			 const __m128 sth = _mm_sin_pd(theta);

                         const __m128d vom= _mm_mul_pd(v,omega);
                         vcth             = _mm_mul_pd(v,cth);
			 vsth             = _mm_mul_pd(v,sth);
			 vomsth           = _mm_fmsub_pd(vDot,cth,
			                                 _mm_mul_pd(vom,sth));
			 vomcth           = _mm_fmadd_pd(vDot,sth,
			                                 _mm_mul_pd(vom,cth));
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      turn_accelerate_xmm2r8_a(const __m128d theta,
		                               const __m128d omega,
		                               const __m128d vDot,
					       const __m128d v,
					       double * __restrict __ATTR_ALIGN__(16)  vcth,
					       double * __restrict __ATTR_ALIGN__(16)  vsth,
					       double * __restrict __ATTR_ALIGN__(16)  vomsth,
					       double * __restrict __ATTR_ALIGN__(16)  vomcth) {


                         const __m128 cth = _mm_cos_pd(theta);
			 const __m128 sth = _mm_sin_pd(theta);

                         const __m128d vom= _mm_mul_pd(v,omega);
                         _mm_store_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_store_pd(&vsth[0],_mm_mul_pd(v,sth));
			 _mm_store_pd(&vomsth[0],_mm_fmsub_pd(vDot,cth,
			                                 _mm_mul_pd(vom,sth)));
			 _mm_store_pd(&vomcth[0],_mm_fmadd_pd(vDot,sth,
			                                 _mm_mul_pd(vom,cth)));
		    }

 
		    
                      __ATTR_ALWAYS_INLINE__
		      
		      static inline void
                      turn_accelerate_xmm2r8_u(const __m128d theta,
		                               const __m128d omega,
		                               const __m128d vDot,
					       const __m128d v,
					       double * __restrict  vcth,
					       double * __restrict  vsth,
					       double * __restrict  vomsth,
					       double * __restrict  vomcth) {


                         const __m128 cth = _mm_cos_pd(theta);
			 const __m128 sth = _mm_sin_pd(theta);

                         const __m128d vom= _mm_mul_pd(v,omega);
                         _mm_storeu_pd(&vcth[0],_mm_mul_pd(v,cth));
			 _mm_storeu_pd(&vsth[0],_mm_mul_pd(v,sth));
			 _mm_storeu_pd(&vomsth[0],_mm_fmsub_pd(vDot,cth,
			                                 _mm_mul_pd(vom,sth)));
			 _mm_storeu_pd(&vomcth[0],_mm_fmadd_pd(vDot,sth,
			                                 _mm_mul_pd(vom,cth)));
		    }


		   
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      turn_accelerate_xmm4r4(const __m128 theta,
		                             const __m128 omega,
		                             const __m128 vDot,
					     const __m128d v,
					     __m128 &vcth,
					     __m128 &vsth,
					     __m128 &vomsth,
					     __m128 &vomcth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         const __m128d vom= _mm_mul_ps(v,omega);
                         vcth             = _mm_mul_ps(v,cth);
			 vsth             = _mm_mul_ps(v,sth);
			 vomsth           = _mm_fmsub_ps(vDot,cth,
			                                 _mm_mul_ps(vom,sth));
			 vomcth           = _mm_fmadd_pd(vDot,sth,
			                                 _mm_mul_ps(vom,cth));
		    }


		     
                      __ATTR_ALWAYS_INLINE__
		     
		      static inline void
                      turn_accelerate_xmm4r4_a(const __m128 theta,
		                             const __m128 omega,
		                             const __m128 vDot,
					     const __m128d v,
					     float * __restrict __ATTR_ALIGN__(16)  vcth,
					     float * __restrict __ATTR_ALIGN__(16)  vsth,
					     float * __restrict __ATTR_ALIGN__(16)  vomsth,
					     float * __restrict __ATTR_ALIGN__(16)  vomcth) {


                         const __m128 cth = _mm_cos_ps(theta);
			 const __m128 sth = _mm_sin_ps(theta);

                         const __m128d vom= _mm_mul_ps(v,omega);
                         _mm_store_ps(&vcth[0],_mm_mul_ps(v,cth));
			 _mm_store_ps(&vsth[0],_mm_mul_ps(v,sth));
			 _mm_store_ps(&vomsth[0],_mm_fmsub_ps(vDot,cth,
			                                 _mm_mul_ps(vom,sth)));
			 _mm_store_ps(&vomcth[0],_mm_fmadd_pd(vDot,sth,
			                                 _mm_mul_ps(vom,cth)));
		    }


		   





		   
     } // math

} //gms








#endif /*__GMS_POS_TO_STATE_SSE_HPP__*/
