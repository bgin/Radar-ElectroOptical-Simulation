

#ifndef __GMS_POS_TO_STATE_AVX512_HPP__
#define __GMS_POS_TO_STATE_AVX512_HPP__ 220520221552


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

 const unsigned int GMS_POS_TO_STATE_AVX512_MAJOR = 1U;
 const unsigned int GMS_POS_TO_STATE_AVX512_MINOR = 0U;
 const unsigned int GMS_POS_TO_STATE_AVX512_MICRO = 0U;
 const unsigned int GMS_POS_TO_STATE_AVX512_FULLVER =
  1000U*GMS_POS_TO_STATE_AVX512_MAJOR+100U*GMS_POS_TO_STATE_AVX512_MINOR+10U*GMS_POS_TO_STATE_AVX512_MICRO;
 const char * const GMS_POS_TO_STATE_AVX512_CREATION_DATE = "22-05-2022 15:52 +00200 (SUN 22 MAY 2022 15:52 GMT+2)";
 const char * const GMS_POS_TO_STATE_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const GMS_POS_TO_STATE_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const GMS_POS_TO_STATE_AVX512_SYNOPSIS      = "AVX512 based position [2D] to state convertion functions (vectorized)."


}


#include <immintrin.h>
#include "GMS_config.h"
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif
#include "GMS_simd_utils.hpp"

namespace gms {

        namespace math {

/*
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
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline 
                      const_velocity_zmm8r8(const __m512d xDot,
					    const __m512d yDot,
					    __m512d * __restrict s_a,
					    __m512d * __restrict s_b) {
#if (USE_SLEEF_LIB) == 1					    
                         *s_a = atan2k(yDot,xDot);
#else
                         *s_a = _mm512_atan2_pd(yDot,xDot);
#endif
                         *s_b = _mm512_fmadd_pd(yDot,yDot,_mm512_mul_pd(xDot,xDot));
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline 
                      const_velocity_zmm16r4(const __m512 xDot,
					     const __m512 yDot,
					     __m512 * __restrict s_a,
					     __m512 * __restrict s_b) {
#if (USE_SLEEF_LIB) == 1					    
                         *s_a = atan2kf(yDot,xDot);
#else
                         *s_a = _mm512_atan2_ps(yDot,xDot);
#endif
                         *s_b = _mm512_fmadd_ps(yDot,yDot,_mm512_mul_ps(xDot,xDot));
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline 
                      const_acceleration_zmm8r8(const __m512d xDot,
					        const __m512d yDot,
		                                const __m512d xDdot,
		                                const __m512d yDdot,
						__m512d * __restrict s_a,
						__m512d * __restrict s_b,
						__m512d * __restrict s_c) {
       
                         __m512d theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask8 m = 0x0;
#if (USE_SLEEF_LIB) == 1
                         theta = atan2k(yDot,xDot);
			 *s_a  = theta;
			 costh = xcos(theta);
			 sinth = xsin(theta);
#else
                         theta = _mm512_atan2_pd(yDot,xDot);
			 *s_a  = theta;
			 costh = _mm512_cos_pd(theta);
			 sinth = _mm512_sin_pd(theta);
#endif
                         *s_b  = _mm512_fmadd_pd(yDot,yDot,_mm512_mul_pd(xDot,xDot));
                         vDot  = _mm512_sqrt_pd(_mm512_fmadd_pd(yDdot,yDdot,
			                                    _mm512_mul_pd(xDdot,xDdot)));
			 t0    = _mm512_sub_pd(_mm512_mul_pd(vDot,costh),xDdot);
			 t1    = _mm512_sub_pd(_mm512_mul_pd(vDot,sinth),yDdot);
			 t2    = _mm512_sub_pd(_mm512_mul_pd(zmm8r8_negate(vDot),costh),xDdot);
			 t3    = _mm512_sub_pd(_mm512_mul_pd(zmm8r8_negate(vDot),sinth),yDdot);
			 diff1 = _mm512_fmadd_pd(t0,t0,_mm512_mul_pd(t1,t1));
			 diff2 = _mm512_fmadd_pd(t2,t2,_mm512_mul_pd(t3,t3));
			 m     = _mm512_cmp_pd_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm512_mask_blend_pd(m,vDot,zmm8r8_negate(vDot));
			 *s_c  = vdot;
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline 
                      const_acceleration_zmm16r4(const __m512 xDot,
					         const __m512 yDot,
		                                 const __m512 xDdot,
		                                 const __m512 yDdot,
						 __m512 * __restrict s_a,
						 __m512 * __restrict s_b,
						 __m512 * __restrict s_c) {
       
                         __m512 theta,costh,sinth,diff1,diff2,vDot,t0,t1,t2,t3;
			 __mmask16 m = 0x0;
#if (USE_SLEEF_LIB) == 1
                         theta = atan2kf(yDot,xDot);
			 *s_a  = theta;
			 costh = xcosf(theta);
			 sinth = xsinf(theta);
#else
                         theta = _mm512_atan2_ps(yDot,xDot);
			 *s_a  = theta;
			 costh = _mm512_cos_ps(theta);
			 sinth = _mm512_sin_ps(theta);
#endif
                         *s_b  = _mm512_fmadd_ps(yDot,yDot,_mm512_mul_ps(xDot,xDot));
                         vDot  = _mm512_sqrt_ps(_mm512_fmadd_ps(yDdot,yDdot,
			                                    _mm512_mul_ps(xDdot,xDdot)));
			 t0    = _mm512_sub_ps(_mm512_mul_ps(vDot,costh),xDdot);
			 t1    = _mm512_sub_ps(_mm512_mul_ps(vDot,sinth),yDdot);
			 t2    = _mm512_sub_ps(_mm512_mul_ps(zmm16r4_negate(vDot),costh),xDdot);
			 t3    = _mm512_sub_ps(_mm512_mul_ps(zmm16r4_negate(vDot),sinth),yDdot);
			 diff1 = _mm512_fmadd_ps(t0,t0,_mm512_mul_ps(t1,t1));
			 diff2 = _mm512_fmadd_ps(t2,t2,_mm512_mul_ps(t3,t3));
			 m     = _mm512_cmp_ps_mask(diff1,diff2,_CMP_LT_OQ);
			 vdot  = _mm512_mask_blend_ps(m,vDot,zmm16r4_negate(vDot));
			 *s_c  = vdot;
		    }
		     

     } // math

} //gms








#endif /*__GMS_POS_TO_STATE_AVX512_HPP__*/
