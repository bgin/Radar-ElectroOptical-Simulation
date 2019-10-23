
#ifndef __GMS_VECSPECFUNCS_H__
#define __GMS_VECSPECFUNCS_H__

namespace file_info {
#if defined _WIN64
    #include "../GMS_version.h"
#elif defined __linux
    #include "GMS_version.h"
#endif
	const unsigned int gGMS_VECSPECFUNCS_MAJOR = gms::common::gVersionInfo.m_VersionMajor;

	const unsigned int gGMS_VECSPECFUNCS_MINOR = gms::common::gVersionInfo.m_VersionMinor;

	const unsigned int gGMS_VECSPECFUNCS_MICRO = gms::common::gVersionInfo.m_VersionMicro;

	const unsigned int gGMS_VECSPECFUNCS_FULLVER = 
			1000U*gGMS_VECSPECFUNCS_MAJOR+100U*gGMS_VECSPECFUNCS_MINOR+10U*gGMS_VECSPECFUNCS_MICRO;

	const char * const  pgGMS_VECSPECFUNCS_CREATE_DATE = "01-10-2018 10:41 +00200 (MON 01 OCT 2018 GMT+2)";

	const char * const  pgGMS_VECSPECFUNCS_BUILD_DATE  = "00-00-0000 00:00";

	const char * const  pgGMS_VECSPECFUNCS_AUTHOR = "Adapted by Bernard Gingold e-mail: beniekg@gmail.com";

	const char * const  pgGMS_VECSPECFUNCS_SYNOPSIS = "Vector SIMD version of special functions library written by: Shanjie Zhang and Jianming Jin.";
}

/*
	Original authors license:
	The original FORTRAN77 version of this routine is copyrighted by
	Shanjie Zhang and Jianming Jin.  However, they give permission to
	incorporate this routine into a user program that the copyright
	is acknowledged.

	Modified:

	Porting to C++ by Bernard Gingold 01/10/2018
	(This is second version which now is based solely on
	 AVX and AVX3 intrinsic C interface. Previous version
	 was based on C++ class semantics with construction
	 default destruction and extensive usage of overloaded
	 operators -- this resulting in suboptimal code generation
	 [operators prolog and epilog was not eliminated by ICC])

	Author:

	Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
	FORTRAN90 version by John Burkardt.

	Reference:

	Shanjie Zhang, Jianming Jin,
	Computation of Special Functions,
	Wiley, 1996,
	ISBN: 0-471-11963-6,
	LC: QA351.C45.
*/

/*
	Bernard Gingold copyright notice:
	This file is a part of LibAtmosModelCPP library
	Copyright(C) 2017 Bernard Gingold
	License : GNU General Public Locense version 3 or later,
	for additional information check file LICENSE.txt in
	project directory.
*/


#include <limits>
#include <cstdint>
//#include <immintrin.h>
//#include <zmmintrin.h>
#include "GMS_simd_macros.h"
#if defined _WIN64
    #include "../GMS_config.h"
    #include "../GMS_common.h"
#elif defined __linux
    #include "GMS_config.h"
    #include "GMS_common.h"
#endif
#include "GMS_avxc2f64.h"
namespace gms {
	namespace math {
#if !defined (GMS_VECSPECFUNCS_UNROLL2X)
#define GMS_VECSPECFUNCS_UNROLL2X 0
#endif

typedef __m256d Vec4;

//#if defined __AVX512F__
typedef __m512d Vec8;

//#endif

		namespace  {
			
			const Vec4 v4PIR(Vec4_SET1(0.318309886183891));
			const Vec4 v4PI(Vec4_SET1(3.141592653589793));
			const Vec4 v4_n0(Vec4_SETZERO);
			const Vec4 v4_n1(Vec4_SET1(1.0));
			const Vec4 v4_nneg1(Vec4_SET1(-1.0));
			const Vec4 v4AIRC1(Vec4_SET1(0.355028053887817));
			const Vec4 v4AIRC2(Vec4_SET1(0.258819403792807));
			const Vec4 v4EPS1(Vec4_SET1(0.000000000000001));
			const Vec4 v4EPS2(Vec4_SET1(0.000000001));
			const Vec4 v4EPS3(Vec4_SET1(1.0E-60));
			const Vec4 v4_1over2(Vec4_SET1(0.5));
			const Vec4 v4_1over3(Vec4_SET1(0.3333333333333333));
			const Vec4 v4_1over4(Vec4_SET1(0.25));
			const Vec4 v4_n1over4(Vec4_SET1(-0.25));
			const Vec4 v4_1over8(Vec4_SET1(0.125));
			const Vec4 v4_n2(Vec4_SET1(2.0));
			const Vec4 v4_nneg2(Vec4_SET1(-2.0));
			const Vec4 v4_n4(Vec4_SET1(4.0));
			const Vec4 v4_n3(Vec4_SET1(3.0));
			const Vec4 v4_n6(Vec4_SET1(6.0));
			const Vec4 v4_n9(Vec4_SET1(9.0));
			const Vec4 v4_n15(Vec4_SET1(15.0));
			const Vec4 v4_ninv15(Vec4_SET1(0.66666666666666666666666666666667));
			const Vec4 v4SR3(Vec4_SET1(1.732050807568877));
			const Vec4 v4HUGE(Vec4_SET1(std::numeric_limits<double>::max()));
			const Vec4 v4NHUGE(Vec4_SET1(std::numeric_limits<double>::lowest()));
			const Vec4 v4TINY(Vec4_SET1(1.0E-35));
			const Vec4 v4_n35(Vec4_SET1(35.0));
			const Vec4 v4_n50(Vec4_SET1(50.0));
			const Vec4 v4_n12(Vec4_SET1(12.0));
			const Vec4 v4_absmask(Vec4_SET1(0x7FFFFFFFFFFFFFF));
			const Vec4 v4_1over180(Vec4_SET1(0.00555555555555555555555555555556));

			const Vec8 v8PIR(Vec8_SET1(0.318309886183891));
			const Vec8 v8PI(Vec8_SET1(3.141592653589793));
			const Vec8 v8_n0(Vec8_SET1(0.0));
			const Vec8 v8_n1(Vec8_SET1(1.0));
			const Vec8 v8_nneg1(Vec8_SET1(-1.0));
			const Vec8 v8AIRC1(Vec8_SET1(0.355028053887817));
			const Vec8 v8AIRC2(Vec8_SET1(0.258819403792807));
			const Vec8 v8EPS1(Vec8_SET1(0.000000000000001));
			const Vec8 v8EPS2(Vec8_SET1(0.000000001));
			const Vec8 v8EPS4(Vec8_SET1(1.0E-60));
			const Vec8 v8EPS1ton12(Vec8_SET1(0.000000000001));
			const Vec8 v8_1over2(Vec8_SET1(0.5));
			const Vec8 v8_1over3(Vec8_SET1(0.33333333333333333333333333));
			const Vec8 v8_1over4(Vec8_SET1(0.25));
			const Vec8 v8_n1over4(Vec8_SET1(-0.25));
			const Vec8 v8_1over8(Vec8_SET1(0.125));
			const Vec8 v8_n1over8(Vec8_SET1(-0.125));
			const Vec8 v8_n2(Vec8_SET1(2.0));
			const Vec8 v8_nneg2(Vec8_SET1(-2.0));
			const Vec8 v8_n3(Vec8_SET1(3.0));
			const Vec8 v8_n4(Vec8_SET1(4.0));
			const Vec8 v8_n6(Vec8_SET1(6.0));
			const Vec8 v8_n9(Vec8_SET1(9.0));
			const Vec8 v8_n15(Vec8_SET1(15.0));
			const Vec8 v8_n18(Vec8_SET1(18.0));
			const Vec8 v8_ninv15(Vec8_SET1(0.66666666666666666666666666666667));
			const Vec8 v8SR3(Vec8_SET1(1.732050807568877));
			const Vec8 v8HUGEP(Vec8_SET1(std::numeric_limits<double>::max()));
			const Vec8 v8HUGEN(Vec8_SET1(std::numeric_limits<double>::lowest()));
			const Vec8 v8TINY(Vec8_SET1(1.0E-35));
			const Vec8 v8_n35(Vec8_SET1(35.0));
			const Vec8 v8_n50(Vec8_SET1(50.0));
			const Vec8 v8_n12(Vec8_SET1(12.0));
			const Vec8 v8_n1idx(Vec8_CVTI4(Vec8_SETI4(1)));
			const Vec8 v8_1over180(Vec8_SET1(0.00555555555555555555555555555556));
			constexpr double s_pi = 3.141592653589793;
		}

		// helpers

		inline Vec4 flip_sign(const Vec4 v) {
			return (Vec4_SUB(v4_n0,v));
		}

		inline Vec8 flip_sign(const Vec8 v) {
			return (Vec8_SUB(v8_n0,v));
		}

		inline Vec4 abs(const Vec4 x) {
			return (Vec4_AND(x,v4_absmask));
		}
		/*
			 AIRYA computes Airy functions and their derivatives.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    30 June 2012
!    01 October 2018 (Bernard Gingold)
!	 
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!    SIMD AVX,AVX3 version by Bernard Gingold
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!       
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Airy function.
!
!    Output, real ( kind = 8 ) AI, BI, AD, BD, the values of Ai(x), Bi(x),
!    Ai'(x), Bi'(x).
		*/
		inline void v4_airya_pd( const Vec4 x,
				         Vec4 &ai,
				         Vec4 &bi,
				         Vec4 &ad,
				         Vec4 &bd) {
		    Vec4 z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2;
			const Vec4 invsr3 = Vec4_SET1(0.577350269189625862351);
			const Vec4 xa = Vec4_SQRT(x);
			const Vec4 z = Vec4_MUL(Vec4_POW(xa,v4_ninv15),v4_ninv15);
			const Vec4 xq = Vec4_SQRT(xa);
			
			// call ajyik
			v4_ajyik_pd(z,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2);
			const Vec4 pirvk2 = Vec4_MUL(v4PIR,vk2);
			const __m256d mask1 = Vec4_CMP(x,v4_n0,_CMP_EQ_OQ);
			const __m256d mask2 = Vec4_CMP(x,v4_n0,_CMP_GT_OQ);
			if (!Vec4_TESTZ(mask1, mask1)) {
				ai = v4AIRC1;
				bi = Vec4_MUL(v4SR3,v4AIRC1);
				ad = flip_sign(v4AIRC2);
				bd = Vec4_MUL(v4SR3,v4AIRC2);
			}
			else if (!Vec4_TESTZ(mask2, mask2)) {
				ai = Vec4_MUL(Vec4_MUL(v4PIR, xq), Vec4_MUL(invsr3,vk1));
				const Vec4 t1 = Vec4_MUL(v4_n2, Vec4_MUL(invsr3,vi1));
				const Vec4 t2 = Vec4_MUL(v4PIR,vk1);
				bi = Vec4_MUL(xq, Vec4_ADD(t1,t2));
				ad = Vec4_MUL(flip_sign(xa), Vec4_MUL(invsr3,pirvk2));
				const Vec4 t3 = Vec4_MUL(v4_n2, Vec4_MUL(invsr3,vi2));
				bd = Vec4_MUL(xa, Vec4_ADD(pirvk2,t3));
			}
			else {
				const Vec4 t1 = Vec4_SUB(vj1, Vec4_MUL(vy1,invsr3));
				ai = Vec4_MUL(v4_1over2, Vec4_MUL(xq,t1));
				const Vec4 t2 = Vec4_ADD(Vec4_MUL(vj1,invsr3),vy1);
				bi = Vec4_MUL(flip_sign(v4_1over2),Vec4_MUL(xq,t2));
				const Vec4 t3 = Vec4_ADD(vj2, Vec4_MUL(vy2,invsr3));
				ad = Vec4_MUL(v4_1over2, Vec4_MUL(xa,t3));
				const Vec4 t4 = Vec4_SUB(Vec4_MUL(vj2,invsr3),vy2);
				bd = Vec4_MUL(v4_1over2, Vec4_MUL(xa,t4));

			}
		}

		void v4_airya_pd_over_rfield1D(const double * __restrict vx,
					       double * __restrict vai,
					       double * __restrict vbi,
					       double * __restrict vad,
					       double * __restrict vbd,
					       const int32_t vlen) {
			 using namespace gms::common;
			if ((vlen % 4) != 0) { return;}
			const bool is_aligned = Is_ptr_aligned32(vx)   &&
				                Is_ptr_aligned32(vai)  &&
						Is_ptr_aligned32(vbi)  &&
						Is_ptr_aligned32(vad)  &&
					        Is_ptr_aligned32(vbd);
			if (is_aligned) {
#if (GMS_VECSPECFUNCS_UNROLL2X) == 1
				for (int32_t i = 0; i != vlen; i += 8) {
					v4_airya_pd(Vec4_LOAD(&vx[i+0]),
						    Vec4_LOAD(&vai[i+0]),
						    Vec4_LOAD(&vbi[i+0]),
						    Vec4_LOAD(&vad[i+0]),
						    Vec4_LOAD(&vbd[i+0]));

					v4_airya_pd(Vec4_LOAD(&vx[i+4]),
						    Vec4_LOAD(&vai[i+4]),
						    Vec4_LOAD(&vbi[i+4]),
						    Vec4_LOAD(&vad[i+4]),
						    Vec4_LOAD(&vbd[i+4]));

				}
#else
				for(int32_t i = 0; i != vlen; i += 4) {
                                      	v4_airya_pd(Vec4_LOAD(&vx[i+0]),
						    Vec4_LOAD(&vai[i+0]),
						    Vec4_LOAD(&vbi[i+0]),
						    Vec4_LOAD(&vad[i+0]),
						    Vec4_LOAD(&vbd[i+0]));
				}
#endif
			}
			else {
#if (GMS_VECSPECFUNCS_UNROLL2X) == 1
				for (int32_t i = 0; i != vlen; i += 8) {
					v4_airya_pd(Vec4_LOADU(&vx[i+0]),
						    Vec4_LOADU(&vai[i+0]),
						    Vec4_LOADU(&vbi[i+0]),
						    Vec4_LOADU(&vad[i+0]),
						    Vec4_LOADU(&vbd[i+0]));
					v4_airya_pd(Vec4_LOADU(&vx[i+4]),
						    Vec4_LOADU(&vai[i+4]),
						    Vec4_LOADU(&vbi[i+4]),
						    Vec4_LOADU(&vad[i+4]),
						    Vec4_LOADU(&vbd[i+4]));
					
				}
#else
				for(int32_t i = 0; i != vlen; i += 4) {
                                         v4_airya_pd(Vec4_LOADU(&vx[i+0]),
						    Vec4_LOADU(&vai[i+0]),
						    Vec4_LOADU(&vbi[i+0]),
						    Vec4_LOADU(&vad[i+0]),
						    Vec4_LOADU(&vbd[i+0]));
				}
#endif
			}

		}

		inline void v8_airya_pd(const Vec8 x,
				        Vec8 &ai,
				        Vec8 &bi,
				        Vec8 &ad,
				        Vec8 &bd) {
			Vec8 z, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2;
			const Vec8 invsr3 = Vec8_SET1(0.577350269189625862351);
			const Vec8 xa = Vec8_SQRT(x);
			const Vec8 z = Vec8_MUL(Vec8_POW(xa,v8_ninv15),v8_ninv15);
			const Vec8 xq = Vec8_SQRT(xa);
			// Call ajyik
			v8_ajyik_pd(x,vj1,vj2,vy1,vy2,vi1,vi2,vk1,vk2);
			const Vec8 pirvk2 = Vec8_MUL(v8PIR,vk2);
			if (Vec8_CMP(x, v8_n0, _CMP_EQ_OQ)) {
				ai = v8AIRC1;
				bi = Vec8_MUL(v8SR3,v8AIRC1);
				ad = flip_sign(v8AIRC2);
				bd = Vec8_MUL(v8SR3,v8AIRC2);
			}
			else if (Vec8_CMP(x, v8_n0, _CMP_GT_OQ)) {
				ai = Vec8_MUL(Vec8_MUL(v8PIR, xq), Vec8_MUL(invsr3, vk1));
				const Vec8 t1 = Vec8_MUL(v8_n2, Vec8_MUL(invsr3, vi1));
				const Vec8 t2 = Vec8_MUL(v8PIR, vk1);
				bi = Vec8_MUL(xq, Vec8_ADD(t1, t2));
				ad = Vec8_MUL(flip_sign(xa), Vec8_MUL(invsr3, pirvk2));
				const Vec8 t3 = Vec8_MUL(v8_n2, Vec8_MUL(invsr3, vi2));
				bd = Vec8_MUL(xa, Vec8_ADD(pirvk2, t3));
			}
			else {
				const Vec8 t1 = Vec8_SUB(vj1, Vec8_MUL(vy1, invsr3));
				ai = Vec8_MUL(v8_1over2, Vec8_MUL(xq, t1));
				const Vec8 t2 = Vec8_ADD(Vec8_MUL(vj1, invsr3), vy1);
				bi = Vec8_MUL(flip_sign(v8_1over2), Vec8_MUL(xq, t2));
				const Vec8 t3 = Vec8_ADD(vj2, Vec8_MUL(vy2, invsr3));
				ad = Vec8_MUL(v8_1over2, Vec8_MUL(xa, t3));
				const Vec8 t4 = Vec8_SUB(Vec8_MUL(vj2, invsr3), vy2);
				bd = Vec8_MUL(v8_1over2, Vec8_MUL(xa, t4));
			}
		}

		void v8_airya_pd_over_rfield1D(const double * __restrict vx,
					       double * __restrict vai,
					       double * __restrict vbi,
					       double * __restrict vad,
					       double * __restrict vbd,
					       const int32_t vlen) {
			using namespace gms::common;
			if ((vlen % 8) != 0) { return;}
			const bool is_aligned = Is_ptr_aligned64(vx)   &&
						Is_ptr_aligned64(vai)  &&
						Is_ptr_aligned64(vbi)  &&
						Is_ptr_aligned64(vad)  &&
						Is_ptr_aligned64(vbd);
			if (is_aligned) {
#if (GMS_VECSPECFUNCS_UNROLL2X) == 1			  
				for (int32_t i = 0; i != vlen; i += 16) {
					v8_airya_pd(Vec8_LOAD(&vx[i+0]),
						    Vec8_LOAD(&vai[i+0]),
						    Vec8_LOAD(&vbi[i+0]),
						    Vec8_LOAD(&vad[i+0]),
						    Vec8_LOAD(&vbd[i+0]));
					 v8_airya_pd(Vec8_LOAD(&vx[i+8]),
						    Vec8_LOAD(&vai[i+8]),
						    Vec8_LOAD(&vbi[i+8]),
						    Vec8_LOAD(&vad[i+8]),
						    Vec8_LOAD(&vbd[i+8]));
				}
#else
                                for (int32_t i = 0; i != vlen; i += 8) {
                                         v8_airya_pd(Vec8_LOAD(&vx[i+0]),
						    Vec8_LOAD(&vai[i+0]),
						    Vec8_LOAD(&vbi[i+0]),
						    Vec8_LOAD(&vad[i+0]),
						    Vec8_LOAD(&vbd[i+0]));
				}
#endif
			}
			else {
#if (GMS_VSPECFUNCS_UNROLL2X) == 1			 
				for (int32_t i = 0; i != vlen; i += 16) {
					v8_airya_pd(Vec8_LOADU(&vx[i+0]),
						        Vec8_LOADU(&vai[i+0]),
						        Vec8_LOADU(&vbi[i+0]),
						        Vec8_LOADU(&vad[i+0]),
						        Vec8_LOADU(&vbd[i+0]));
					v8_airya_pd(Vec8_LOADU(&vx[i+8]),
						        Vec8_LOADU(&vai[i+8]),
						        Vec8_LOADU(&vbi[i+8]),
						        Vec8_LOADU(&vad[i+8]),
						        Vec8_LOADU(&vbd[i+8]));
				}
#else
                                for (int32_t i = 0; i != vlen; i += 8) {
                                         v8_airya_pd(Vec8_LOADU(&vx[i+0]),
						        Vec8_LOADU(&vai[i+0]),
						        Vec8_LOADU(&vbi[i+0]),
						        Vec8_LOADU(&vad[i+0]),
						        Vec8_LOADU(&vbd[i+0]));
				}
#endif
			}
		}

		/*
				AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
!
!  Discussion: 
!
!    Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions 
!    Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!    01 October 2018 (Bernard Gingold)

!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    SIMD AVX,AVX3 versions by Bernard Gingold
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.  X should not be zero.
!
!    Output, real ( kind = 8 ) VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2,
!    the values of J1/3(x), J2/3(x), Y1/3(x), Y2/3(x), I1/3(x), I2/3(x),
!    K1/3(x), K2/3(x).
		*/
		inline void v4_ajyik_pd(const Vec4 x,
				        Vec4 &vj1,
				        Vec4 &vj2,
				        Vec4 &vy1,
				        Vec4 &vy2,
				        Vec4 &vi1,
				        Vec4 &vi2,
				        Vec4 &vk1,
				        Vec4 &vk2) {
			const __m256d mask_xeq0 = Vec4_CMP(x,v4_n0,_CMP_EQ_OQ);
			if (!Vec4_TESTZ(mask_xeq0, mask_xeq0)) {
				vj1 = v4_n0;
				vj2 = v4_n0;
				vy1 = v4HUGE;
				vy2 = v4NHUGE;
				vi1 = v4_n0;
				vi2 = v4_n0;
				vk1 = v4NHUGE;
				vk2 = v4NHUGE;
				return;
			}
		    Vec4 r,vl,vsl,vjl,a0,b0,vv,px,rp,qx,rq,xk,ck,sk,
			     ifl1_t1,ifl1_t2,uj1,uj2,pv1,pv2,vil,c1,gn,
				 sum;
			r       = v4_n0; 
			vl      = v4_n0; 
			vsl     = v4_n0; 
			vjl     = v4_n0;
			a0      = v4_n0; 
			b0      = v4_n0; 
			vv      = v4_n0; 
			px      = v4_n0;
			rp      = v4_n0; 
			qx      = v4_n0; 
			rq      = v4_n0; 
			xk      = v4_n0;
			ck      = v4_n0; 
			sk      = v4_n0; 
			ifl1_t1 = v4_n0; 
			ifl1_t2 = v4_n0; 
			uj1     = v4_n0; 
			uj2     = v4_n0;
			pv1     = v4_n0; 
			pv2     = v4_n0; 
			vil     = v4_n0;
			c1      = v4_n0; 
			gn      = v4_n0; 
			sum     = v4_n0; // Memory first touch (if stored on stack)
			const Vec4 rp2 = Vec4_SET1(0.63661977236758);
			const Vec4 gp1 = Vec4_SET1(0.892979511569249);
			const Vec4 gp2 = Vec4_SET1(0.902745292950934);
			const Vec4 gn1 = Vec4_SET1(1.3541179394264);
			const Vec4 gn2 = Vec4_SET1(2.678938534707747);
			const Vec4 vv0 = Vec4_SET1(0.444444444444444);
			const Vec4 uu0 = Vec4_SET1(1.1547005383793);
			const Vec4 c0 = flip_sign(Vec4_SET1(0.78125E-02));
			const Vec4 x2 = Vec4_MUL(x,x);
			const Vec4 invx = Vec4_DIV(v4_n1,x);
			int32_t k,k0,l;
			const Vec4 mask_xlt35 = Vec4_CMP(x,v4_n35,_CMP_LT_OQ);
			const Vec4 mask_xlt50 = Vec4_CMP(x,v4_n50,_CMP_LT_OQ);
			if (!Vec4_TESTZ(mask_xlt35, mask_xlt35)) {
				k0 = 12;
			}
			else if (!Vec4_TESTZ(mask_xlt50, mask_xlt50)) {
				k0 = 10;
			}
			else {
				k0 = 8;
			}
			const Vec4 mask_xlt12 = Vec4_CMP(x,v4_n12,_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_xlt12, mask_xlt12)) {
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vl = Vec4_MUL(tl,v4_1over3);
					vjl = v4_n1;
					r = v4_n1;
					for (k = 1; k != 40; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 rhs = Vec4_MUL(tk, Vec4_ADD(tk,vl));
						const Vec4 lhs = Vec4_MUL(v4_n1over4, Vec4_MUL(r,x2));
						r = Vec4_DIV(lhs,rhs);
						vjl = Vec4_ADD(vjl,r);
						const Vec4 mask_rleps = Vec4_CMP(abs(r),v4EPS1,_CMP_LT_OQ);
						if (!Vec4_TESTZ(mask_rleps, mask_rleps)) {
							break;
						}
					}
					a0 = Vec4_POW(Vec4_MUL(v4_1over2,x),vl);
					if (l == 1) {
						vj1 = Vec4_DIV(a0, Vec4_MUL(gp1,vjl));
					}
					else {
						vj2 = Vec4_DIV(a0, Vec4_MUL(gp2,vjl));
					}
				}
			}
			else {
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vv = Vec4_MUL(vv0, Vec4_MUL(tl,tl));
					px = v4_n1;
					rp = v4_n1;
					for (k = 1; k != k0; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 st1 = Vec4_SUB(Vec4_MUL(v4_n4,tk),v4_n3);
						const Vec4 st1p = Vec4_MUL(st1,st1);
						const Vec4 t1 = Vec4_SUB(vv,st1p);
						const Vec4 st2 = Vec4_SUB(Vec4_MUL(v4_n4,tk),v4_n1);
						const Vec4 st2p = Vec4_MUL(st2,st2);
						const Vec4 t2 = Vec4_SUB(vv,st2p);
						const Vec4 st3 = Vec4_SUB(Vec4_MUL(v4_n2,tk),v4_n1);
						const Vec4 t3 = Vec4_MUL(tk, Vec4_MUL(st3,x2));
						const Vec4 t4 = Vec4_DIV(Vec4_MUL(t1,t2),t3);
						rp = Vec4_MUL(Vec4_MUL(c0,rp),t4);
						px = Vec4_ADD(px,rp);
					}
					qx = v4_n1;
					rq = v4_n1;
					for (k = 1; k != k0; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 st1 = Vec4_SUB(Vec4_MUL(v4_n4, tk), v4_n1);
						const Vec4 st1p = Vec4_MUL(st1, st1);
						const Vec4 t1 = Vec4_SUB(vv, st1p);
						const Vec4 st2 = Vec4_ADD(Vec4_MUL(v4_n4, tk), v4_n1);
						const Vec4 st2p = Vec4_MUL(st2, st2);
						const Vec4 t2 = Vec4_SUB(vv, st2p);
						const Vec4 st3 = Vec4_SUB(Vec4_MUL(v4_n2, tk), v4_n1);
						const Vec4 t3 = Vec4_MUL(tk, Vec4_MUL(st3, x2));
						const Vec4 t4 = Vec4_DIV(Vec4_MUL(t1,t2),t3);
						rq = Vec4_MUL(Vec4_MUL(c0,rq),t4);
						qx = Vec4_ADD(qx,rq);
					}
					const Vec4 t1 = Vec4_MUL(v4_1over8,Vec4_SUB(vv,v4_n1));
					qx = Vec4_MUL(t1, Vec4_MUL(qx,x));
					const Vec4 t2 = Vec4_ADD(Vec4_MUL(v4_1over2, Vec4_MUL(tl,v4_1over3)),v4_1over4);
					xk = Vec4_SUB(x, Vec4_MUL(t2,v4PI));
					a0 = Vec4_SQRT(Vec4_MUL(rp2,invx));
					ck = Vec4_COS(xk);
					sk = Vec4_SIN(xk);
					ifl1_t1 = Vec4_SUB(Vec4_MUL(px, ck), Vec4_MUL(qx, sk));
					ifl1_t2 = Vec4_ADD(Vec4_MUL(px, sk), Vec4_MUL(qx, ck));
					if (l == 1) {
						vj1 = Vec4_MUL(a0, ifl1_t1);
						vy1 = Vec4_MUL(a0, ifl1_t2);
					}	
					else {
						vj2 = Vec4_MUL(a0,ifl1_t1);
						vy2 = Vec4_MUL(a0,ifl1_t2);
					}
				}	
			}		
			
			if (!Vec4_TESTZ(mask_xlt12, mask_xlt12)) {
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vl = Vec4_MUL(tl, v4_1over3);
					vjl = v4_n1;
					r = v4_n1;
					for (k = 1; k != 40; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 rhs = Vec4_MUL(tk, Vec4_SUB(tk,vl));
						const Vec4 lhs = Vec4_MUL(flip_sign(v4_1over4), Vec4_MUL(r,x2));
						r = Vec4_DIV(lhs,rhs);
						vjl = Vec4_ADD(vjl,r);
						const Vec4 mask_rleps = Vec4_CMP(abs(r), v4EPS1, _CMP_LT_OQ);
						if (!Vec4_TESTZ(mask_rleps, mask_rleps)) {
							break;
						}
					}
					b0 = Vec4_POW(Vec4_MUL(v4_n2,invx),vl);
					if (l == 1) {
						uj1 = Vec4_MUL(b0, Vec4_DIV(vjl,gn1));
					}
					else {
						uj2 = Vec4_MUL(b0, Vec4_DIV(vjl,gn2));
					}
				}
				pv1 = Vec4_MUL(v4PI,v4_1over2);
				pv2 = Vec4_MUL(v4PI, Vec4_SET1(0.66666666666666666667));
				vy1 = Vec4_MUL(uu0, Vec4_MUL(vj1, Vec4_SUB(Vec4_COS(pv1),uj1)));
				vy2 = Vec4_MUL(uu0, Vec4_MUL(vj2, Vec4_SUB(Vec4_COS(pv2),uj2)));
			}
			
			const Vec4 mask_xlt18 = Vec4_CMP(x, Vec4_SET1(18.0),_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_xlt18, mask_xlt18)) {
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vl = Vec4_MUL(tl, v4_1over3);
					vil = v4_n1;
					r = v4_n1;
					for (k = 1; k != 40; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 rhs = Vec4_MUL(tk, Vec4_ADD(tk,vl));
						const Vec4 lhs = Vec4_MUL(v4_1over4, Vec4_MUL(r,x2));
						r = Vec4_DIV(lhs,rhs);
						vil = Vec4_ADD(vil,r);
						const Vec4 mask_rleps = Vec4_CMP(abs(r), v4EPS1, _CMP_LT_OQ);
						if (!Vec4_TESTZ(mask_rleps, mask_rleps)) {
							break;
						}
					}
					a0 = Vec4_POW(Vec4_MUL(v4_1over2,x),vl);
					if (l == 1) {
						vi1 = Vec4_MUL(Vec4_DIV(a0,gp1),vil);
					}
					else {
						vi2 = Vec4_MUL(Vec4_DIV(a0,gp2),vil);
					}
				}
			}
			else {
				c1 = Vec4_DIV(Vec4_EXP(x), Vec4_SQRT(Vec4_MUL(Vec4_SET1(5.1415926535897),x)));
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vv = Vec4_MUL(vv0, Vec4_MUL(tl,tl));
					vsl = v4_n1;
					r = v4_n1;
					for (k = 1; k != k0; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 st1 = Vec4_SUB(vv,Vec4_SUB(Vec4_MUL(v4_n2,tk),v4_n1));
						const Vec4 st1p = Vec4_MUL(st1,st1);
						const Vec4 rhs = Vec4_DIV(st1p, Vec4_MUL(tk,x));
						r = Vec4_MUL(Vec4_MUL(v4_1over8,r),rhs);
						vsl = Vec4_ADD(vsl,r);
					}
					if (l == 1) {
						vi1 = Vec4_MUL(c1,vsl);
					}
					else {
						vi2 = Vec4_MUL(c1,vsl);
					}
				}
			}
			const Vec4 mask_xlte9 = Vec4_CMP(x,v4_n9,_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_xlte9, mask_xlte9)) {
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vl = Vec4_MUL(tl,v4_1over3);
					if (l == 1) {
						gn = gn1;
					}
					else {
						gn = gn2;
					}
					a0 = Vec4_DIV(Vec4_POW(Vec4_MUL(v4_n2,invx),vl),gn);
					sum = v4_n0;
					r = v4_n0;
					for (k = 1; k != 60; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 rhs = Vec4_MUL(tk, Vec4_SUB(tk,vl));
						r = Vec4_DIV(Vec4_MUL(v4_1over4, Vec4_MUL(r,x2)),rhs);
						sum = Vec4_ADD(sum,r);
						const Vec4 mask_rleps = Vec4_CMP(abs(r), v4EPS1, _CMP_LT_OQ);
						if (!Vec4_TESTZ(mask_rleps, mask_rleps)) {
							break;
						}
					}
					if (l == 1) {
						const Vec4 rhs = Vec4_SUB(Vec4_MUL(sum,a0),vi1);
						vk1 = Vec4_MUL(Vec4_MUL(v4_1over2, uu0), Vec4_MUL(v4PI,rhs));
					}
					else {
						const Vec4 rhs = Vec4_SUB(Vec4_MUL(sum, a0), vi1);
						vk2 = Vec4_MUL(Vec4_MUL(v4_1over2, uu0), Vec4_MUL(v4PI, rhs));
					}
				}
			}
			else {
				c1 = Vec4_MUL(Vec4_EXP(flip_sign(x)), Vec4_SQRT(Vec4_MUL(v4_1over2, Vec4_MUL(v4PI,invx))));
				for (l = 1; l <= 2; ++l) {
					const Vec4 tl = Vec4_CVTI4(Vec4_SETI4(l));
					vv = Vec4_MUL(vv0, Vec4_MUL(tl,tl));
					sum = v4_n1;
					r = v4_n1;
					for (k = 1; k != k0; ++k) {
						const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
						const Vec4 st1 = Vec4_SUB(Vec4_MUL(v4_n2,tk),v4_n1);
						const Vec4 st1p = Vec4_MUL(st1,st1);
						const Vec4 rhs = Vec4_DIV(Vec4_SUB(vv, st1p), Vec4_MUL(tk,x));
						r = Vec4_MUL(v4_1over8, Vec4_MUL(r,rhs));
						sum = Vec4_ADD(sum,r);
					}
					if (l == 1) {
						vk1 = Vec4_MUL(c1,sum);
					}
					else {
						vk2 = Vec4_MUL(c1,sum);
					}
				}
			}
		}

		void v4_ajyik_pd_over_rfield1D(const double * __restrict x,
					       double * __restrict  j1,
					       double * __restrict  j2,
					       double * __restrict  y1,
					       double * __restrict  y2,
					       double * __restrict  i1,
					       double * __restrict  i2,
					       double * __restrict  k1,
					       double * __restrict  k2
					       const int32_t vlen   ) {
		             using namespace gms::common;
			     if((vlen % 4) != 0) {return;}
			     const bool is_aligned =  Is_ptr_aligned32(x)   &&
			                              Is_ptr_aligned32(j1)  &&
			                              Is_ptr_aligned32(j2)  &&
			                              Is_ptr_aligned32(y1)  &&
			                              Is_ptr_aligned32(y2)  &&
			                              Is_ptr_aligned32(i1)  &&
			                              Is_ptr_aligned32(i2)  &&
			                              Is_ptr_aligned32(k1)  &&
			                              Is_ptr_aligned32(k2);
                             if(is_aligned) {
#if (GMS_VECSPECFUNCS_UNROLL2X) == 1
			       for (int32_t i = 0; i != vlen; i += 8) {
				 v4_ajyik_pd(Vec4_LOAD(&x[i+0]),
					     Vec4_LOAD(&j1[i+0]),
					     Vec4_LOAD(&j2[i+0]),
					     Vec4_LOAD(&y1[i+0]),
					     Vec4_LOAD(&y2[i+0]),
					     Vec4_LOAD(&i1[i+0]),
					     Vec4_LOAD(&i2[i+0]),
					     Vec4_LOAD(&k1[i+0]),
					     Vec4_LOAD(&k2[i+0]));
				  v4_ajyik_pd(Vec4_LOAD(&x[i+4]),
					     Vec4_LOAD(&j1[i+4]),
					     Vec4_LOAD(&j2[i+4]),
					     Vec4_LOAD(&y1[i+4]),
					     Vec4_LOAD(&y2[i+4]),
					     Vec4_LOAD(&i1[i+4]),
					     Vec4_LOAD(&i2[i+4]),
					     Vec4_LOAD(&k1[i+4]),
					     Vec4_LOAD(&k2[i+4]));
			       }
#else 
                               for (int32_t i = 0; i != vlen; i += 4) {
                                    v4_ajyik_pd(Vec4_LOAD(&x[i+0]),
					     Vec4_LOAD(&j1[i+0]),
					     Vec4_LOAD(&j2[i+0]),
					     Vec4_LOAD(&y1[i+0]),
					     Vec4_LOAD(&y2[i+0]),
					     Vec4_LOAD(&i1[i+0]),
					     Vec4_LOAD(&i2[i+0]),
					     Vec4_LOAD(&k1[i+0]),
					     Vec4_LOAD(&k2[i+0]));
			       }
#endif
			     }
			     else {
#if (GMS_VECSPECFUNCS_UNROLL2X) == 1
                             for (int32_t i = 0; i != vlen; i += 8) {
				 v4_ajyik_pd(Vec4_LOADU(&x[i+0]),
					     Vec4_LOADU(&j1[i+0]),
					     Vec4_LOADU(&j2[i+0]),
					     Vec4_LOADU(&y1[i+0]),
					     Vec4_LOADU(&y2[i+0]),
					     Vec4_LOADU(&i1[i+0]),
					     Vec4_LOADU(&i2[i+0]),
					     Vec4_LOADU(&k1[i+0]),
					     Vec4_LOADU(&k2[i+0]));
				  v4_ajyik_pd(Vec4_LOADU(&x[i+4]),
					     Vec4_LOADU(&j1[i+4]),
					     Vec4_LOADU(&j2[i+4]),
					     Vec4_LOADU(&y1[i+4]),
					     Vec4_LOADU(&y2[i+4]),
					     Vec4_LOADU(&i1[i+4]),
					     Vec4_LOADU(&i2[i+4]),
					     Vec4_LOADU(&k1[i+4]),
					     Vec4_LOADU(&k2[i+4]));
			       }
#else
                              for (int32_t i = 0; i != vlen; i += 4) {
                                    v4_ajyik_pd(Vec4_LOADU(&x[i+0]),
					     Vec4_LOADU(&j1[i+0]),
					     Vec4_LOADU(&j2[i+0]),
					     Vec4_LOADU(&y1[i+0]),
					     Vec4_LOADU(&y2[i+0]),
					     Vec4_LOADU(&i1[i+0]),
					     Vec4_LOADU(&i2[i+0]),
					     Vec4_LOADU(&k1[i+0]),
					     Vec4_LOADU(&k2[i+0]));
			       }
#endif
			     }

		}

		void v8_ajyik_pd(const Vec8 x,
			         Vec8 &vj1,
			         Vec8 &vj2,
			         Vec8 &vy1,
			         Vec8 &vy2,
			         Vec8 &vi1,
			         Vec8 &vi2,
			         Vec8 &vk1,
			         Vec8 &vk2) {
			if (Vec8_CMP(x, v8_n0, _CMP_EQ_OQ)) {
				vj1 = v8_n0;
				vj2 = v8_n0;
				vy1 = v8HUGEN;
				vy2 = v8HUGEP;
				vi1 = v8_n0;
				vi2 = v8_n0;
				vk1 = v8HUGEN;
				vk2 = v8HUGEN;
			}
			const Vec8 rp2 = Vec8_SET1(0.63661977236758);
			const Vec8 gp1 = Vec8_SET1(0.892979511569249);
			const Vec8 gp2 = Vec8_SET1(0.902745292950934);
			const Vec8 gn1 = Vec8_SET1(1.3541179394264);
			const Vec8 gn2 = Vec8_SET1(2.678938534707747);
			const Vec8 vv0 = Vec8_SET1(0.444444444444444);
			const Vec8 uu0 = Vec8_SET1(1.1547005383793);
			const Vec8 x2 = Vec8_MUL(x,x);
			const Vec8 invx = Vec8_DIV(v8_n1,x);
			const Vec8 coeff = Vec8_SET1(-0.78125E-02);
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 a0,b0,c0,ck,pv1,pv2,px,qx,r,rp,rp2,
				     rq,sk,sum,uj1,uj2,uu0,vl,vsl,vv,xk,
					 iflt_t1,iflt_t2,vjl,vil,gn;
			}ca;
#elif defined __linux
			__attribute__((align(64))) struct {
                                 Vec8 a0,b0,c0,ck,pv1,pv2,px,qx,r,rp,rp2,
				     rq,sk,sum,uj1,uj2,uu0,vl,vsl,vv,xk,
					 iflt_t1,iflt_t2,vjl,vil,gn;
			}ca;  
			int32_t k0,l,k;
			ca.a0      = v8_n0;
			ca.b0      = v8_n0;
			ca.c0      = v8_n0;
			ca.ck      = v8_n0;
			ca.pv1     = v8_n0;
			ca.pv2     = v8_n0;
			ca.px      = v8_n0;
			ca.qx      = v8_n0;
			ca.r       = v8_n0;
			ca.rp      = v8_n0;
			ca.rp2     = v8_n0;
			ca.rq      = v8_n0;
			ca.sk      = v8_n0;
			ca.sum     = v8_n0;
			ca.uj1     = v8_n0;
			ca.uj2     = v8_n0;
			ca.uu0     = v8_n0;
			ca.vl      = v8_n0;
			ca.vsl     = v8_n0;
			ca.vv      = v8_n0;
			ca.xk      = v8_n0;
			ca.iflt_t1 = v8_n0;
			ca.iflt_t2 = v8_n0;
			ca.vjl     = v8_n0;
			ca.vil     = v8_n0;
			ca.gn      = v8_n0;
			if (Vec8_CMP(x, v8_n35, _CMP_LT_OQ)) {
			   k0 = 12;
			}
			else if (Vec8_CMP(x, v8_n50, _CMP_LT_OQ)) {
				k0 = 10;
			}
			else {
				k0 = 8;
			}

			if (Vec8_CMP(x, v8_n12, _CMP_LE_OQ)) {
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vl = Vec8_MUL(tl,v8_1over3);
					ca.vjl = v8_n1;
					ca.r = v8_n1;
					for (k = 1; k <= 40; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 rhs = Vec8_MUL(tk, Vec8_ADD(tk,ca.vl));
						ca.r = Vec8_DIV(Vec8_MUL(v8_n1over4, Vec8_MUL(ca.r,x2)),rhs);
						ca.vjl = Vec8_ADD(ca.vjl,ca.r);
						if (Vec8_CMP(Vec8_ABS(ca.r), v8EPS1, _CMP_LT_OQ)) {
							break;
						}
					}
					ca.a0 = Vec8_POW(Vec8_MUL(v8_1over2,x),ca.vl);
					if (l == 1) {
						vj1 = Vec8_MUL(Vec8_DIV(ca.a0,gp1),ca.vjl);
					}
					else {
						vj2 = Vec8_MUL(Vec8_DIV(ca.a0,gp2),ca.vjl);
					}

				}
			}
			else {
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vv = Vec8_MUL(vv0, Vec8_MUL(tl,tl));
					ca.px = v8_n1;
					ca.rp = v8_n1;
					for (k = 1; k != k0; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 st1 = Vec8_SUB(Vec8_MUL(v8_n4,tk),v8_n3);
						const Vec8 st1p = Vec8_SUB(ca.vv,Vec8_MUL(st1,st1));
						const Vec8 st2 = Vec8_SUB(Vec8_MUL(v8_n4,tk),v8_n1);
						const Vec8 st2p = Vec8_SUB(ca.vv, Vec8_MUL(st2,st2));
						const Vec8 st3 = Vec8_SUB(Vec8_MUL(v8_n2,tk),v8_n1);
						const Vec8 st3t = Vec8_MUL(tk, Vec8_MUL(st3,x2));
						const Vec8 t1 = Vec8_DIV(Vec8_MUL(st1p,st2p),st3t);
						ca.rp = Vec8_MUL(Vec8_MUL(coeff,ca.rp),t1);
						ca.px = Vec8_ADD(ca.px,ca.rp);
					}
					ca.qx = v8_n1;
					ca.rq = v8_n1;
					for (k = 1; k != k0; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 st1 = Vec8_SUB(Vec8_MUL(v8_n4, tk), v8_n1);
						const Vec8 st1p = Vec8_SUB(ca.vv, Vec8_MUL(st1, st1));
						const Vec8 st2 = Vec8_ADD(Vec8_MUL(v8_n4, tk), v8_n1);
						const Vec8 st2p = Vec8_SUB(ca.vv, Vec8_MUL(st2, st2));
						const Vec8 st3 = Vec8_ADD(Vec8_MUL(v8_n2, tk), v8_n1);
						const Vec8 st3t = Vec8_MUL(tk, Vec8_MUL(st3, x2));
						const Vec8 t1 = Vec8_DIV(Vec8_MUL(st1p, st2p), st3t);
						ca.rq = Vec8_MUL(Vec8_MUL(coeff,ca.rq),t1);
						ca.qx = Vec8_ADD(ca.qx,ca.rq);
					}
					 ca.qx =  Vec8_MUL(v8_1over8,Vec8_MUL(Vec8_SUB(ca.vv, v8_n1), Vec8_MUL(ca.qx,invx)));
					 const Vec8 rhs = Vec8_ADD(Vec8_MUL(v8_1over2, Vec8_MUL(tl,v8_1over3)),v8_1over4);
					 ca.xk = Vec8_SUB(x, Vec8_MUL(rhs,v8PI));
					 ca.a0 = Vec8_SQRT(Vec8_MUL(rp2,invx));
					 ca.ck = Vec8_SIN(ca.xk);
					 ca.sk = Vec8_COS(ca.xk);
					 ca.iflt_t1 = Vec8_MUL(ca.a0, Vec8_SUB(Vec8_MUL(ca.px, ca.ck), Vec8_MUL(ca.qx,ca.sk)));
					 ca.iflt_t2 = Vec8_MUL(ca.a0, Vec8_ADD(Vec8_MUL(ca.px, ca.sk), Vec8_MUL(ca.qx,ca.ck)));
					 if (l == 1) {
					    vj1 = ca.iflt_t1;
						vy1 = ca.iflt_t2;
					 }
					 else {
						 vj2 = ca.iflt_t1;
						 vy2 = ca.iflt_t2;
					 }
				}
			}

			if (Vec8_CMP(x, v8_n12, _CMP_LE_OQ)) {
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vl = Vec8_MUL(tl,v8_1over3);
					ca.vjl = v8_n1;
					ca.r = v8_n1;
					for (k = 1; k != 40; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 rhs = Vec8_MUL(tk, Vec8_SUB(tk,ca.vl));
						ca.r = Vec8_DIV(Vec8_MUL(v8_n1over4, Vec8_MUL(ca.r,x2)),rhs);
						ca.vjl = Vec8_ADD(ca.vjl,ca.r);
						if (Vec8_CMP(Vec8_ABS(ca.r), v8EPS1, _CMP_LT_OQ)) {
							break;
						}
					}
					ca.b0 = Vec8_POW(Vec8_MUL(v8_n2,invx),ca.vl);
					if (l == 1) {
						ca.uj1 = Vec8_MUL(ca.b0, Vec8_DIV(ca.vjl,gn1));
					}
					else {
						ca.uj2 = Vec8_MUL(ca.b0, Vec8_DIV(ca.vjl,gn2));
					}
					
				}
				ca.pv1 = Vec8_MUL(v8PI,v8_1over3);
				ca.pv2 = Vec8_MUL(v8PI, Vec8_SET1(0.66666666666666666667));
				vy1 = Vec8_MUL(ca.uu0, Vec8_SUB(Vec8_MUL(vj1, Vec8_COS(ca.pv1)),ca.uj1));
				vy2 = Vec8_MUL(ca.uu0, Vec8_SUB(Vec8_MUL(vj2, Vec8_COS(ca.pv2)),ca.uj2));
			}

			if (Vec8_CMP(x, v8_n18, _CMP_LE_OQ)) {
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vl = Vec8_MUL(tl, v8_1over3);
					ca.vil = v8_n1;
					ca.r = v8_n1;
					for (k = 1; k != 41; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 rhs = Vec8_MUL(tk, Vec8_ADD(tk,ca.vl));
						ca.r = Vec8_DIV(Vec8_MUL(v8_1over4, Vec8_MUL(ca.r,x2)),rhs);
						ca.vil = Vec8_ADD(ca.vil,ca.r);
						if (Vec8_CMP(Vec8_ABS(ca.r), v8EPS1, _CMP_LT_OQ)) {
							break;
						}
					}
					ca.a0 = Vec8_POW(Vec8_MUL(v8_1over2,x),ca.vl);
					if (l == 1) {
						vi1 = Vec8_MUL(Vec8_DIV(ca.a0,gp1),ca.vil);
					}
					else {
						vi2 = Vec8_MUL(Vec8_DIV(ca.a0,gp2),ca.vil);
					}
				}
			}
			else {
				ca.c0 = Vec8_DIV(Vec8_EXP(x), Vec8_SQRT(Vec8_MUL(v8_n2, Vec8_MUL(v8PI,x))));
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vv = Vec8_MUL(vv0, Vec8_MUL(tl,tl));
					ca.vsl = v8_n1;
					ca.r = v8_n1;
					for (k = 1; k != k0; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 st1 = Vec8_SUB(Vec8_MUL(v8_n2,tk),v8_n1);
						const Vec8 st1p = Vec8_DIV(Vec8_SUB(ca.vv,Vec8_MUL(st1,st1)),Vec8_MUL(tk,x));
						ca.r = Vec8_MUL(Vec8_MUL(v8_n1over8,ca.r),st1p);
						ca.vsl = Vec8_ADD(ca.vsl,ca.r);
					}
					const Vec8 tmp = Vec8_MUL(ca.c0, ca.vsl);
					if (l == 1) {
						vi1 = tmp;
					}
					else {
						vi2 = tmp;
					}
				}
			}

			if (Vec8_CMP(x, v8_n9, _CMP_LE_OQ)) {
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vl = Vec8_MUL(tl, v8_1over3);
					ca.gn = Vec8_BLEND(Vec8_CMP(tl, v8_n1idx,_CMP_EQ_OQ), gn1, gn2);
					ca.a0 = Vec8_DIV(Vec8_POW(Vec8_MUL(v8_n2,invx),ca.vl),ca.gn);
					ca.sum = v8_n1;
					ca.r = v8_n1;
					for (k = 1; k != 60; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 rhs = Vec8_MUL(tk, Vec8_SUB(tk,ca.vl));
						ca.r = Vec8_DIV(Vec8_MUL(Vec8_MUL(ca.r,x2),v8_1over4),rhs);
						ca.sum = Vec8_ADD(ca.sum,ca.r);
						if (Vec8_CMP(Vec8_ABS(ca.r), v8EPS1, _CMP_LT_OQ)) {
							break;
						}
					}
					const Vec8 tmp = Vec8_MUL(v8_1over2,Vec8_MUL(ca.uu0,v8PI));
					if (l == 1) {
						vk1 = Vec8_MUL(tmp, Vec8_SUB(Vec8_MUL(ca.sum,ca.a0),vi1));
					}
					else {
						vk2 = Vec8_MUL(tmp, Vec8_SUB(Vec8_MUL(ca.sum, ca.a0), vi2));
					}
				}
			}
			else {
				ca.c0 = Vec8_MUL(Vec8_EXP(flip_sign(x)), Vec8_SQRT(Vec8_MUL(v8_1over2, Vec8_MUL(v8PI, invx))));
				for (l = 1; l != 3; ++l) {
					const Vec8 tl = Vec8_CVTI4(Vec8_SETI4(l));
					ca.vv = Vec8_MUL(vv0, Vec8_MUL(tl, tl));
					ca.sum = v8_n1;
					ca.r = v8_n1;
					for (k = 1; k != k0; ++k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
						const Vec8 st1 = Vec8_SUB(Vec8_MUL(v8_n2, tk), v8_n1);
						const Vec8 st1p = Vec8_MUL(st1, st1);
						const Vec8 rhs = Vec8_DIV(Vec8_SUB(ca.vv, st1p), Vec8_MUL(tk, x));
						ca.r = Vec8_MUL(v8_1over8, Vec8_MUL(ca.r, rhs));
						ca.sum = Vec8_ADD(ca.sum, ca.r);
					}
					const Vec8 tmp = Vec8_MUL(ca.c0, ca.sum);
					if (l == 1) {
						vk1 = tmp;
					}
					else {
						vk2 = tmp;
					}
				}
			}
		}

		/*
			
!
!! SPHY computes spherical Bessel functions yn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!	 04 October 2018 Bernard Gingold

!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
! 
		*/
		template<int32_t N>
		void v4_sphy_pd(int32_t nx,
				Vec4 x,
				int32_t &nm,
				Vec4 __restrict sy[N],
				Vec4 __restrict dy[N]) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec4 f,f0,f1,invx,tsin,tcos;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#elif defined __linux
			__attribute__((align(64))) struct {
                                Vec4 f,f0,f1,invx,tsin,tcos;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#endif
			nm = nx;
			ca.invx = Vec4_DIV(v4_n1,x);
			Vec4 mask_lteps3 = Vec4_CMP(x,v4EPS3,_CMP_LT_OQ);
			if (!Vec4_TESTZ(mask_lteps3, mask_lteps3)) {
				for (ca.k = 0; ca.k != nx, ++ca.k) {
					sy[ca.k] = v4NHUGE;
					dy[ca.k] = v4HUGE;
				}
				return;
			}
			ca.tsin = Vec4_SIN(x);
			ca.tcos = Vec4_COS(x);
			sy[0] = flip_sign(Vec4_MUL(ca.tcos,invx));
			sy[1] = Vec4_MUL(Vec4_SUB(sy[0],ca.tsin,invx);
			ca.f0 = sy[0];
			ca.f1 = sy[1];
			for (ca.k = 2; ca.k != nx; ++ca.k) {
				const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
				const Vec4 lhs = Vec4_MUL(Vec4_SUB(Vec4_MUL(v4_n2,k),v4_n1),ca.f1);
				ca.f = Vec4_SUB(Vec4_MUL(lhs,invx),ca.f0);
				sy[ca.k] = f;
				const Vec4 mask_geeps = Vec4_CMP(Vec4_abs(x),v4HUGE,_CMP_GE_OQ);
				if (!Vec4_TESTZ(mask_geeps, mask_geeps)) {
					break;
				}
				ca.f0 = ca.f1;
				ca.f1 = ca.f;
			}
			nm = ca.k - 1;
			dy[0] = Vec4_MUL(Vec4_MUL(Vec4_ADD(ca.tsin,ca.tcos),invx),invx);
			for (ca.k = 1; ca.k != nm; ++v = ca.k) {
				const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
				const Vec4 lhs = Vec4_SUB(sy[ca.k - 1], Vec4_ADD(tk,v4_n1));
				dy[ca.k] = Vec4_MUL(lhs, Vec4_MUL(sy[ca.k],invx));
			}
		}

		template<int32_t N>
		void v8_sphy_pd(int32_t nx,
				Vec8 x,
			       _int32_t &nm,
			       Vec8 __restrict sy[N],
			       Vec8 __restrict dy[N]) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 f,f0,f1,invx,tsin,tcos;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
	                        Vec8 f,f0,f1,invx,tsin,tcos;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#endif
			nm = nx;
			ca.invx = Vec8_DIV(v8_n1,x);
			if (Vec8_CMP(x, v8EPS3, _CMP_LT_OQ)) {
				for (ca.k = 0; ca.k != nx; ++ca.k) {
					sy[ca.k] = v8HUGEP;
					dy[ca.k] = v8HUGEN;
				}
				return;
			}
			ca.tsin = Vec8_SIN(x);
			ca.tcos = Vec8_COS(x);
			sy[0] = Vec8_MUL(ca.tcos,invx);
			sy[1] = Vec8_MUL(Vec8_SUB(sy[0],ca.tsin),invx);
			ca.f0 = sy[0];
			ca.f1 = sy[1];
			for (ca.k = 2; ca.k != nx; ++ca.k) {
				const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
				const Vec8 rhs = Vec8_SUB(Vec8_MUL(ca.f1,invx),ca.f0);
				ca.f = Vec8_MUL(Vec8_SUB(Vec8_MUL(v8_n2,tk),v8_n1),rhs);
				sy[ca.k] = ca.f;
				if (Vec8_CMP(Vec8_ABS(ca.f), v8HUGEP, _CMP_GE_OQ)) {
					break;
				}
				ca.f0 = ca.f1;
				ca.f1 = ca.f;
			}
			nm = ca.k - 1;
			dy[0] = Vec8_MUL(Vec8_MUL(Vec8_ADD(ca.tsin,ca.tcos),invx),invx);
			for (ca.k = 1; ca.k != nm; ++ca.k) {
				const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
				const Vec8 lhs = Vec8_SUB(sy[ca.k - 1], Vec8_ADD(tk, v8_n1));
				dy[ca.k] = Vec8_MUL(lhs, Vec8_MUL(sy[ca.k], invx));
			}
		}

		/*
			SPHK computes modified spherical Bessel functions kn(x) and derivatives.
!
!  Discussion:
!
!    This procedure computes modified spherical Bessel functions
!    of the second kind, kn(x) and kn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!	 04 October 2018 Bernard Gingold

!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) SK(0:N), DK(0:N), the values of kn(x) and kn'(x).
!
		*/

		template<int32_t N>
		void v4_sphk_pd(const int32_t nx,
			        Vec4 x,
			        int32_t &nm,
			        Vec4 __restrict sk[N],
			        Vec4 __restrict dk[N]) {
#if defined _WIN64		  
			__declspec(align(64))  struct {
				Vec4 f,f0,f1,invx;
				int32_t k
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec4 f,f0,f1,invx;
				int32_t k
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#endif
			nm = nx;
			ca.invx = Vec4_DIV(v4_n1,x);
			const Vec4 mask_lteps = Vec4_CMP(x,v4EPS3,_CMP_LT_OQ);
			if (!Vec4_TESTZ(mask_lteps, mask_lteps)) {
				for (ca.k = 0; ca.k != nx; ++ca.k) {
					sk[ca.k] = v4HUGE;
					dk[ca.k] = v4NHUGE;
				}
				return;
			}
			sk[0] = Vec4_MUL(Vec4_MUL(v4_1over2, v4PI), Vec4_MUL(invx, Vec4_EXP(flip_sign(x))));
			sk[1] = Vec4_MUL(sk[0], Vec4_ADD(v4_n1, Vec4_MUL(v4_n1,invx)));
			ca.f0 = sk[0];
			ca.f1 = sk[1];
			for (ca.k = 2; ca.k != nx; ++ca.k) {
				const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
				const Vec4 lhs = Vec4_MUL(Vec4_SUB(Vec4_MUL(v4_n2,tk),v4_n1),ca.f1);
				ca.f = Vec4_MUL(lhs, Vec4_ADD(invx,ca.f0));
				sk[ca.k] = ca.f;
				const Vec4 mask_geeps = Vec4_CMP(abs(ca.f),v4HUGE,_CMP_GT_OQ);
				if (!Vec4_TESTZ(mask_geeps, mask_geeps)) {
					break;
				}
				ca.f0 = ca.f1;
				ca.f1 = ca.f;
			}
			nm = ca.k - 1;
			dk[0] = flip_sign(sk[1]);
			for (ca.k = 1; ca.k != nm; ++ca.k) {
				const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
				const Vec4 rhs = Vec4_MUL(Vec4_MUL(Vec4_ADD(tk, v4_n1), invx), sk[ca.k]);
				dk[ca.k] = Vec4_SUB(flip_sign(sk[ca.k-1]),rhs);
			}
		}

		template<int32_t N>
		void v8_sphk_pd(const int32_t nx,
			        Vec8 x,
			        int32_t &nm,
			        Vec8 __restrict sk[N],
			        Vec8 __restrict dk[N]) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 f,f0,f1,invx;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec8 f,f0,f1,invx;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#endif
			ca.invx = Vec8_DIV_(v8_n1,x);
			nm = nx;
			if (Vec8_CMP(x, v8EPS4, _CMP_LT_OQ)) {
				for (ca.k = 0; ca.k != nx; ++ca.k) {
					sk[k] = v8HUGEP;
					dk[k] = v8HUGEN;
				}
				return;
			}
			sk[0] = Vec8_MUL(Vec8_MUL(v8_1over2, v8PI), Vec8_MUL(invx, Vec8_EXP(flip_sign(x))));
			sk[1] = Vec8_MUL(sk[0], Vec8_ADD(v8_n1,Vec8_MUL(v8_n1,invx)));
			ca.f0 = sk[0];
			ca.f1 = sk[1];
			for (ca.k = 2; ca.k != nx; ++ca.k) {
				const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
				const Vec8 lhs = Vec8_SUB(Vec8_MUL(v8_n2,tk),v8_n1);
				ca.f = Vec8_MUL(lhs, Vec8_ADD(ca.f1,invx),ca.f0);
				sk[k] = ca.f;
				if (Vec8_CMP(Vec8_ABS(ca.f), v8HUGEP, _CMP_GT_OQ)) {
					break;
				}
				ca.f0 = ca.f1;
				ca.f1 = ca.f;
			}
			nm = ca.k - 1;
			dk[0] = flip_sign(sk[1]);
			for (ca.k = 1; ca.k != nm; ++ca.k) {
				const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
				dk[k] = Vec8_SUB(flip_sign(sk[ca.k - 1]), Vec8_MUL(Vec8_ADD(tk,v8_n1),Vec8_MUL(invx,sk[k])));
			}
		}

		/*
			RCTY computes Riccati-Bessel function of the second kind, and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!    05 October 2018
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of yn(x).
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) RY(0:N), the values of x yn(x).
!
!    Output, real ( kind = 8 ) DY(0:N), the values of [x yn(x)]'.
		*/

		template<int32_t N>
		void v4_rcty_pd(const Vec4 x,
			        int32_t & nm,
			        Vec4 __restrict ry[N],
			        Vec4 __restrict dy[N]) {
#if defined _WIN64		  
			__declspec(align(64)) struct {
				Vec4 rf0,rf1,rf2,invx,tsin;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec4 rf0,rf1,rf2,invx,tsin;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#endif
			ca.invx = Vec4_DIV(v4_n1,x);
			nm = N;
			const Vec4 mask_lt1e60 = Vec4_CMP(x,v4EPS3,_CMP_LT_OQ);
			if (!Vec4_TESTZ(mask_lt1e60, mask_lt1e60)) {
				for (ca.k = 0; ca.k != N; ++ca.k) {
					ry[ca.k] = v4HUGE;
					dy[ca.k]  = v4NHUGE;
				}
				ry[0] = v4_nneg1;
				dy[0] = v4_n0;
			}
			const Vec4 cosx = Vec4_COS(x);
			ca.tsin = Vec4_SIN(x);
			ry[0]  = flip_sign(cosx);
			ry[1]  = Vec4_SUB(Vec4_MUL(ry[0],invx),ca.tsin);
			ca.rf0 = ry[0];
			ca.rf1 = ry[1];
			for (ca.k = 2; ca.k != N; ++ca.k) {
				 const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
				 const Vec4 lhs = Vec4_SUB(Vec4_MUL(v4_n2,tk),v4_n1);
				 ca.rf2 = Vec4_MUL(lhs, Vec4_SUB(Vec4_MUL(ca.rf1,invx),ca.rf0));
				 Vec4 const mask_gt1e300 = Vec4_CMP(abs(ca.rf2), v4HUGE, _CMP_GT_OQ);
				 if (!Vec4_TESTZ(mask_gt1e300, mask_gt1e300)) {
					 break;
				 }
				 ry[ca.k] = ca.rf2;
				 ca.rf0 = ca.rf1;
				 ca.rf1 = ca.rf2;
			}
			nm = ca.k - 1;
			dy[0] = ca.tsin;
			for (ca.k = 1; ca.k != mn; ++ca.k) {
				const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
				dy[ca.k] = Vec4_ADD(Vec4_MUL(flip_sign(tk), Vec4_MUL(ry[ca.k], invx)), ry[ca.k-1]);
			}
		}

		


		template<int32_t N>
		void v8_rcty_pd(const Vec8 x,
			        int32_t & nm,
			        Vec8 __restrict ry[N],
			        Vec8 __restrict dy[N]) {
#if defined _WIN64		  
			__declspec(align(64)) struct {
				Vec8 rf0, rf1, rf2, invx, tsin;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0, 60)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec8 rf0, rf1, rf2, invx, tsin;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0, 60)
#endif
			}ca;
#endif
			ca.invx = Vec8_DIV(v8_n1, x);
			nm = N;
			
			if (Vec8_CMP(x,v8EPS4,_CMP_LT_OQ)) {
				for (ca.k = 0; ca.k != N; ++ca.k) {
					ry[ca.k] = v8HUGEP;
					dy[ca.k] = v4HUGEN;
				}
				ry[0] = v8_nneg1;
				dy[0] = v8_n0;
			}
			const Vec8 cosx = Vec8_COS(x);
			ca.tsin = Vec8_SIN(x);
			ry[0] = flip_sign(cosx);
			ry[1] = Vec8_SUB(Vec8_MUL(ry[0], invx), ca.tsin);
			ca.rf0 = ry[0];
			ca.rf1 = ry[1];
			for (ca.k = 2; ca.k != N; ++ca.k) {
				const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
				const Vec8 lhs = Vec8_SUB(Vec8_MUL(v8_n2, tk), v8_n1);
				ca.rf2 = Vec8_MUL(lhs, Vec8_SUB(Vec8_MUL(ca.rf1, invx), ca.rf0));
				
				if (Vec8_CMP(abs(ca.rf2), v8HUGEP, _CMP_GT_OQ)) {
					break;
				}
				ry[ca.k] = ca.rf2;
				ca.rf0 = ca.rf1;
				ca.rf1 = ca.rf2;
			}
			nm = ca.k - 1;
			dy[0] = ca.tsin;
			for (ca.k = 1; ca.k != mn; ++ca.k) {
				const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
				dy[ca.k] = Vec8_ADD(Vec8_MUL(flip_sign(tk), Vec8_MUL(ry[ca.k], invx)), ry[ca.k - 1]);
			}
		}

		/*
			 LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!	 05 October 2018
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of Qn(x).
!
!    Input, real ( kind = 8 ) X, the argument of Qn(x).
!
!    Output, real ( kind = 8 ) QN(0:N), QD(0:N), the values of
!    Qn(x) and Qn'(x).
		*/

		template<int32_t N>
		void v4_lqnb_pd(Vec4 x,
			        Vec4 __restrict qn[N],
			        Vec4 __restrict qd[N]) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec4 q0,q1,qc1,qc2,qf,qf0,qf1,qf2,qr,x2;
				int32_t j;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,4)
#endif
			    int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,4)
#endif
				int32_t l;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2,4)
#endif
			    int32_t nl;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3,4)
#endif
#if (USE_STRUCT_PADDING) == 1
			    PAD_TO(4,32)
#endif
			}ca;
#elif defined __linux
			  __attribute__((align(64))) struct {
                                Vec4 q0,q1,qc1,qc2,qf,qf0,qf1,qf2,qr,x2;
				int32_t j;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,4)
#endif
			    int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,4)
#endif
				int32_t l;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2,4)
#endif
			    int32_t nl;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3,4)
#endif
#if (USE_STRUCT_PADDING) == 1
			    PAD_TO(4,32)
#endif
			}ca;  
#endif
			const Vec4 maskeq1 = Vec_CMP(x,v4_n1,_CMP_EQ_OQ);
			if (Vec4_TESTZ(maskeq1, maskeq1)) {
				for (ca.k = 0; ca.k != N; ++ca.k) {
					qn[ca.k] = v4HUGE;
					qd[ca.k] = v4HUGE;
				}
				return;
			}
			ca.q0   = v4_n0;
			ca.q1   = v4_n0;
			ca.qc1  = v4_n0;
			ca.qc2  = v4_n0;
			ca.qf   = v4_n0;
			ca.qf0  = v4_n0;
			ca.qf1  = v4_n0;
			ca.qf2  = v4_n0;
			ca.qr   = v4_n0;
			ca.x2   = v4_n0;
			const Vec4 eps = Vec4_SET1(1.0E-14);
			const Vec4 v4_c1021 = Vec4_SET1(1.021E+00);
			const Vec4 mask_xlc1021 = Vec4_CMP(x,v4_c1021,_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_xlc1021, mask_xlc1021)) {
				ca.x2 = abs(Vec4_DIV(Vec4_ADD(v4_n1, x), Vec4_SUB(v4_n1,x)));
				ca.q0 = Vec4_MUL(v4_1over2, Vec4_LOG(x2));
				ca.q1 = Vec4_SUB(Vec4_MUL(x,ca.q0),v4_n1);
				qn[0] = ca.q0;
				qn[1] = ca.q1;
				qd[0] = Vec4_DIV(v4_n1, Vec4_SUB(v4_n1, Vec4_MUL(x,x)));
				qd[1] = Vec4_ADD(qn[0], Vec4_MUL(x, qd[0]));
				for (ca.k = 2; ca.k != N; ++ca.k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
					const Vec4 lhs = Vec4_MUL(Vec4_SUB(Vec4_MUL(v4_n2, tk), v4_n1), Vec4_MUL(x,ca.q1));
					const Vec4 rhs = Vec4_MUL(Vec4_SUB(tk,v4_n1),ca.q0);
					ca.qf = Vec4_DIV(Vec4_SUB(lhs,rhs),tk);
					qn[ca.k] = ca.qf;
					const Vec4 rhs2 = Vec4_SUB(v4_n1,Vec4_MUL(x,x));
					const Vec4 lhs2 = Vec4_MUL(Vec4_SUB(qn[ca.k - 1], Vec4_MUL(x,ca.qf)),tk);
					qd[ca.k] = Vec4_DIV(lhs2,rhs2);
					ca.q0 = ca.q1;
					ca.q1 = ca.qf;
				}
			}
			else {
				ca.qc2 = Vec4_DIV(v4_n1,x);
				for (ca.j = 0; ca.j != N; ++ca.j) {
					const Vec4 tj = Vec4_CVTI4(Vec4_SETI4(ca.j));
					const Vec4 rhs = Vec4_MUL(Vec4_ADD(Vec4_MUL(v4_n2,tj),v4_n1),x);
					ca.qc2 = Vec4_DIV(Vec4_MUL(ca.qc2,ca.j),rhs);
					if (j == N - 1) {
						ca.qc1 = ca.qc2;
					}
				}
				for (ca.l = 0; ca.l != 2; ++ca.l) {
				     ca.nl = N + ca.l;
					 ca.qf = v4_n1;
					 ca.qr = v4_n1;
					 const Vec4 tnl = Vec4_CVTI4(Vec4_SETI4(ca.nl));
					 for (ca.k = 1; ca.k != 500; ++ca.k) {
						 const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
						 const Vec4 st1 = Vec4_ADD(Vec4_MUL(v4_1over2, tnl), Vec4_SUB(tk,v4_n1));
						 const Vec4 st2 = Vec4_MUL(v4_1over2,Vec4_ADD(Vec4_SUB(tnl,v4_n1),tk));
						 const Vec4 st3 = Vec4_SUB(Vec4_ADD(tnl,tk),v4_1over2);
						 const Vec4 den = Vec4_MUL(st3, Vec4_MUL(tk,Vec4_MUL(x,x)));
						 ca.qr = Vec4_DIV(Vec4_MUL(ca.qr, Vec4_MUL(st1,st2)),den);
						 ca.qf = Vec4_ADD(ca.qf,ca.qr);
						 const Vec4 mask_lteps = Vec4_CMP(abs(Vec4_DIV(ca.qr,ca.)),eps,_CMP_LT_OQ);
						 if (!Vec4_TESTZ(mask_lteps, mask_lteps)) {
							 break;
						 }
					 }
					 if (ca.l == 0) {
						 qn[N - 1] = Vec4_MUL(ca.qf,ca.qc1);
					 }
					 else {
						 qn[N] = Vec4_MUL(ca.qf,ca.qc2);
					 }
				}
				ca.qf2 = qn[N];
				ca.qf1 = qn[N-1];
				for (ca.k = N; ca.k != 2; --ca.k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
					const Vec4 st1 = Vec4_SUB(Vec4_MUL(v4_n2,tk),v4_n1);
					const Vec4 st2 = Vec4_SUB(Vec4_MUL(st1, Vec4_MUL(x, ca.qf1)), Vec4_MUL(tk,ca.qf2));
					ca.qf0 = Vec4_DIV(st2, Vec4_SUB(tk,v4_n1));
					qn[ca.k-2] = ca.qf0;
					ca.qf2 = ca.qf1;
					ca.qf1 = ca.qf0;
				}
				qd[0] = Vec4_DIV(v4_n1, Vec4_SUB(v4_n1, Vec4_MUL(x,x)));
				for (ca.k = 1; ca.k != N; ++ca.k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
					const Vec4 lhs = Vec4_SUB(qn[ca.k - 1], Vec4_MUL(x, qn[ca.k]));
					const Vec4 rhs = Vec4_SUB(v4_n1, Vec4_MUL(x,x));
					qd[ca.k] = Vec4_MUL(tk, Vec4_DIV(lhs,rhs));
				}
			}
		}

		template<int32_t N>
		void v8_lqnb_pd(Vec8 x,
			        Vec8 __restrict qn[N],
			        Vec8 __restrict qd[N]) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 q0, q1, qc1, qc2, qf, qf0, qf1, qf2, qr, x2;
				int32_t j;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0, 4)
#endif
					int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1, 4)
#endif
					int32_t l;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2, 4)
#endif
					int32_t nl;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3, 4)
#endif
#if (USE_STRUCT_PADDING) == 1
			    PAD_TO(4,32)
#endif
			}ca;
#elif defined __linux
			__attribute__((align(64))) struct {
                                Vec8 q0, q1, qc1, qc2, qf, qf0, qf1, qf2, qr, x2;
				int32_t j;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0, 4)
#endif
					int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1, 4)
#endif
					int32_t l;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(2, 4)
#endif
					int32_t nl;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(3, 4)
#endif
#if (USE_STRUCT_PADDING) == 1
			    PAD_TO(4,32)
#endif 
			}ca;
#endif
			if (Vec8_CMP(x, v4_n1, _CMP_EQ_OQ)) {
				for (ca.k = 0; ca.k != N; ++ca.k) {
					qn[ca.k] = v8HUGEP;
					qd[ca.k] = v8HUGEP;
				}
				return;
			}
			ca.q0 = v8_n0;
			ca.q1 = v8_n0;
			ca.qc1 = v8_n0;
			ca.qc2 = v8_n0;
			ca.qf = v8_n0;
			ca.qf0 = v8_n0;
			ca.qf1 = v8_n0;
			ca.qf2 = v8_n0;
			ca.qr = v8_n0;
			ca.x2 = v8_n0;
			const Vec8 eps = Vec8_SET1(1.0E-14);
			const Vec8 v8_c1021 = Vec8_SET1(1.021E+00);
			
			if (Vec8_CMP(x, v8_c1021, _CMP_LE_OQ)) {
				ca.x2 = Vec8_ABS(Vec8_DIV(Vec8_ADD(v8_n1, x), Vec8_SUB(v8_n1, x)));
				ca.q0 = Vec8_MUL(v8_1over2, Vec8_LOG(x2));
				ca.q1 = Vec8_SUB(Vec8_MUL(x, ca.q0), v8_n1);
				qn[0] = ca.q0;
				qn[1] = ca.q1;
				qd[0] = Vec8_DIV(v8_n1, Vec8_SUB(v8_n1, Vec8_MUL(x, x)));
				qd[1] = Vec8_ADD(qn[0], Vec8_MUL(x, qd[0]));
				for (ca.k = 2; ca.k != N; ++ca.k) {
					const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
					const Vec8 lhs = Vec8_MUL(Vec8_SUB(Vec8_MUL(v8_n2, tk), v8_n1), Vec8_MUL(x, ca.q1));
					const Vec8 rhs = Vec8_MUL(Vec8_SUB(tk, v8_n1), ca.q0);
					ca.qf = Vec8_DIV(Vec8_SUB(lhs, rhs), tk);
					qn[ca.k] = ca.qf;
					const Vec8 rhs2 = Vec8_SUB(v8_n1, Vec8_MUL(x, x));
					const Vec8 lhs2 = Vec8_MUL(Vec8_SUB(qn[ca.k - 1], Vec8_MUL(x, ca.qf)), tk);
					qd[ca.k] = Vec8_DIV(lhs2, rhs2);
					ca.q0 = ca.q1;
					ca.q1 = ca.qf;
				}
			}
			else {
				ca.qc2 = Vec8_DIV(v8_n1, x);
				for (ca.j = 0; ca.j != N; ++ca.j) {
					const Vec8 tj = Vec8_CVTI4(Vec8_SETI4(ca.j));
					const Vec8 rhs = Vec8_MUL(Vec8_ADD(Vec8_MUL(v8_n2, tj), v8_n1), x);
					ca.qc2 = Vec8_DIV(Vec8_MUL(ca.qc2, ca.j), rhs);
					if (j == N - 1) {
						ca.qc1 = ca.qc2;
					}
				}
				for (ca.l = 0; ca.l != 2; ++ca.l) {
					ca.nl = N + ca.l;
					ca.qf = v8_n1;
					ca.qr = v8_n1;
					const Vec8 tnl = Vec8_CVTI4(Vec8_SETI4(ca.nl));
					for (ca.k = 1; ca.k != 500; ++ca.k) {
						const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
						const Vec8 st1 = Vec8_ADD(Vec8_MUL(v8_1over2, tnl), Vec8_SUB(tk, v8_n1));
						const Vec8 st2 = Vec8_MUL(v8_1over2, Vec8_ADD(Vec8_SUB(tnl, v8_n1), tk));
						const Vec8 st3 = Vec8_SUB(Vec8_ADD(tnl, tk), v8_1over2);
						const Vec8 den = Vec8_MUL(st3, Vec8_MUL(tk, Vec8_MUL(x, x)));
						ca.qr = Vec8_DIV(Vec8_MUL(ca.qr, Vec8_MUL(st1, st2)), den);
						ca.qf = Vec8_ADD(ca.qf, ca.qr);
						
						if (Vec8_CMP(Vec8_ABS(Vec4_DIV(ca.qr, ca.)), eps, _CMP_LT_OQ)) {
							break;
						}
					}
					if (ca.l == 0) {
						qn[N - 1] = Vec8_MUL(ca.qf, ca.qc1);
					}
					else {
						qn[N] = Vec8_MUL(ca.qf, ca.qc2);
					}
				}
				ca.qf2 = qn[N];
				ca.qf1 = qn[N - 1];
				for (ca.k = N; ca.k != 2; --ca.k) {
					const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
					const Vec8 st1 = Vec8_SUB(Vec8_MUL(v8_n2, tk), v8_n1);
					const Vec8 st2 = Vec8_SUB(Vec8_MUL(st1, Vec8_MUL(x, ca.qf1)), Vec8_MUL(tk, ca.qf2));
					ca.qf0 = Vec8_DIV(st2, Vec8_SUB(tk, v8_n1));
					qn[ca.k - 2] = ca.qf0;
					ca.qf2 = ca.qf1;
					ca.qf1 = ca.qf0;
				}
				qd[0] = Vec8_DIV(v8_n1, Vec8_SUB(v8_n1, Vec8_MUL(x, x)));
				for (ca.k = 1; ca.k != N; ++ca.k) {
					const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(ca.k));
					const Vec8 lhs = Vec8_SUB(qn[ca.k - 1], Vec8_MUL(x, qn[ca.k]));
					const Vec8 rhs = Vec8_SUB(v8_n1, Vec8_MUL(x, x));
					qd[ca.k] = Vec8_MUL(tk, Vec8_DIV(lhs, rhs));
				}
			}
		}

		/*
			 E1XA computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!    06 October 2018 by Bernard Gingold
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
		*/

		inline void v4_e1xa_pd(Vec4 x,
				      Vec4 &e1) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec4 es1,es2;

			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                 Vec4 es1,es2;
			}ca;
#endif
			const Vec4 mask_eq0 = Vec4_CMP(x,v4_n0,_CMP_EQ_OQ);
			const Vec4 mask_lt1 = Vec4_CMP(x,v4_n1,_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_eq0, mask_eq0)) {
			     x = v4HUGE;
				return;
			}
			else if (!Vec4_TESTZ(mask_lt1, mask_lt1)) {
				const Vec4 coef1 = Vec4_SET1(1.07857E-03);
				const Vec4 coef2 = Vec4_SET1(9.76004E-03);
				const Vec4 coef3 = Vec4_SET1(5.519968E-02);
				const Vec4 coef4 = Vec4_SET1(0.24991055E+00);
				const Vec4 coef5 = Vec4_SET1(0.99999193E+00);
				const Vec4 coef6 = Vec4_SET1(0.57721566E+00);
				const Vec4 t0 = flip_sign(Vec4_LOG(x));
				x = Vec4_FMSUB(coef1,x,coef2);
				x = Vec4_FMAD(x,x,coef3);
				x = Vec4_FMSUB(x,x,coef4);
				x = Vec4_FMAD(x,x,coef5);
				x = Vec4_FMSUB(x,xcoef6);
				e1 = Vec8_ADD(t0,x);
			}
			else {
				const Vec4 coef1 = Vec4_SET1(8.5733287401E+00);
				const Vec4 coef2 = Vec4_SET1(18.059016973E+00);
				const Vec4 coef3 = Vec4_SET1(8.6347608925E+00);
				const Vec4 coef4 = Vec4_SET1(0.2677737343E+00);
				x = Vec4_MUL(Vec4_ADD(x,coef1),x);
				x = Vec4_FMAD(x,coef2,x);
				x = Vec4_FMAD(x,coef3,x);
				ca.es1 = Vec4_ADD(x,coef4);
				const Vec4 coefa1 = Vec4_SET1(9.5733223454E+00);
				const Vec4 coefa2 = Vec4_SET1(25.6329561486E+00);
				const Vec4 coefa3 = Vec4_SET1(21.0996530827E+00);
				const Vec4 coefa4 = Vec4_SET1(3.9584969228E+00);
				x = Vec4_MUL(Vec4_ADD(x,coefa1),x);
				x = Vec4_FMAD(x,coefa2,x);
				x = Vec4_FMAD(x,coefa3,x);
				ca.es2 = Vec4_ADD(x,coefa4);
				e1 = Vec4_MUL(Vec4_DIV(Vec4_EXP(flip_sign(x)), x), Vec4_DIV(es1,es2));
			}
		}

		inline void v8_e1xa_pd(const Vec8 x,
				       Vec8 &el) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 es1,es2;
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec8 es1,es2;
			}ca;
#endif
			if (Vec8_CMP(x, v8_n0, _CMP_EQ_OQ)) {
				e1 = v8HUGEP;
				return;
			}
			else if (Vec8_CMP(x, v8_n1, _CMP_LE_OQ)) {
				const Vec8 coef1 = Vec8_SET1(1.07857E-03);
				const Vec8 coef2 = Vec8_SET1(9.76004E-03);
				const Vec8 coef3 = Vec8_SET1(5.519968E-02);
				const Vec8 coef4 = Vec8_SET1(0.24991055E+00);
				const Vec8 coef5 = Vec8_SET1(0.99999193E+00);
				const Vec8 coef6 = Vec8_SET1(0.57721566E+00);
				const Vec8 t0 = flip_sign(Vec8_LOG(x));
				x = Vec4_FMSUB(coef1, x, coef2);
				x = Vec4_FMAD(x, x, coef3);
				x = Vec4_FMSUB(x, x, coef4);
				x = Vec4_FMAD(x, x, coef5);
				x = Vec4_FMSUB(x, xcoef6);
				e1 = Vec8_ADD(t0,x);
			}
			else {
				const Vec8 coef1 = Vec8_SET1(8.5733287401E+00);
				const Vec8 coef2 = Vec8_SET1(18.059016973E+00);
				const Vec8 coef3 = Vec8_SET1(8.6347608925E+00);
				const Vec8 coef4 = Vec8_SET1(0.2677737343E+00);
				x = Vec8_MUL(Vec8_ADD(x, coef1), x);
				x = Vec8_FMAD(x, coef2, x);
				x = Vec8_FMAD(x, coef3, x);
				ca.es1 = Vec8_ADD(x, coef4);
				const Vec8 coefa1 = Vec8_SET1(9.5733223454E+00);
				const Vec8 coefa2 = Vec8_SET1(25.6329561486E+00);
				const Vec8 coefa3 = Vec8_SET1(21.0996530827E+00);
				const Vec8 coefa4 = Vec8_SET1(3.9584969228E+00);
				x = Vec8_MUL(Vec8_ADD(x, coefa1), x);
				x = Vec8_FMAD(x, coefa2, x);
				x = Vec8_FMAD(x, coefa3, x);
				ca.es2 = Vec8_ADD(x, coefa4);
				e1 = Vec8_MUL(Vec8_DIV(Vec8_EXP(flip_sign(x)), x), Vec8_DIV(es1, es2));
			}
		}

		/*
				EIX computes the exponential integral Ei(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!	 06 October 2018 Bernard Gingold
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) EI, the function value.
		*/

		inline void v4_eix_pd(const Vec4 x,
				      Vec4 &ei) {
			Vec4 r,invx;
			int32_t k;
			const Vec4 eps = Vec4_SET1(1.0E-15);
			const Vec4 mask_eq0 = Vec4_CMP(x,v4_n0,_CMP_EQ_OQ);
			const Vec4 mask_lt40 = Vec4_CMP(x, Vec4_SET1(40.0),_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_eq0, mask_eq0)) {
				ei = v4NHUGE;
			}
			else if (!Vec4_TESTZ(mask_eq0, mask_eq0)) {
				ei = v4_n1;
				r = v4_n1;
				for (k = 1; k != 100; ++k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
					const Vec4 rhs = Vec4_MUL(Vec4_ADD(tk, v4_n1), Vec4_ADD(tk,v4_n1));
					r = Vec4_DIV(Vec4_MUL(r, Vec4_MUL(tk,x)),rhs);
					ei = Vec4_ADD(ei,r);
					const Vec4 mask_lteps = Vec4_CMP(abs(Vec4_DIV(r,ei)),eps,_CMP_LE_OQ);
					if (!Vec4_TESTZ(mask_lteps, mask_lteps)) {
						break;
					}
				}
				Vec4 ga = Vec4_SET1(0.5772156649015328E+00);
				ei = Vec4_ADD(Vec4_ADD(ga, Vec4_LOG(x)), Vec4_MUL(x,ei));
			}
			else {
			    
				ei = v4_n1;
				r = v4_n1;
				invx = Vec4_DIV(v4_n1,x);
				for (k = 1; k != 20; ++k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(k));
					r = Vec4_MUL(r, Vec4_MUL(tk,invx));
					ei = Vec4_ADD(ei,r);
				}
				ei = Vec4_MUL(Vec4_MUL(Vec4_EXP(x),invx),ei);
			}
		}

		inline void v8_eix_pd(const Vec8 x,
				      Vec8 &ei) {
			Vec8 r,invx;
			int32_t k;
			if (Vec8_CMP(x, v8_n0, _CMP_EQ_OQ)) {
				ei = v8HUGEN;
				return;
			}
			else if (Vec8_CMP(x, Vec8_SET1(40.0), _CMP_LE_OQ)) {
				const Vec8 eps = Vec8_SET1(1.0E-15);
				ei = v8_n1;
				r = v8_n1;
				for (k = 1; k != 100; ++k) {
					const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
					const Vec8 rhs = Vec8_MUL(Vec8_ADD(tk, v8_n1), Vec8_ADD(tk, v8_n1));
					r = Vec8_DIV(Vec8_MUL(r, Vec8_MUL(tk, x)), rhs);
					ei = Vec8_ADD(ei, r);
					if (Vec8_CMP(Vec8_ABS(r, ei), eps, _CMP_LE_OQ)) {
						break;
					}
				}
				const Vec8 ga = Vec8_SET1(0.5772156649015328E+00);
				ei = Vec8_ADD(Vec8_ADD(ga, Vec8_LOG(x)), Vec8_MUL(x, ei));
			}
			else {
			    
				ei = v8_n1;
				r = v8_n1;
				invx = Vec8_DIV(v8_n1,x);
				for (k = 1; k != 20; ++k) {
					const Vec8 tk = Vec8_CVTI4(Vec8_SETI4(k));
					r = Vec8_MUL(r, Vec8_MUL(tk, invx));
					ei = Vec8_ADD(ei, r);
				}
				ei = Vec8_MUL(Vec8_MUL(Vec8_EXP(x), invx), ei);
			}
		}

		/*
			ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    12 July 2012
!	 06 October 2018 Bernard Gingold
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HK, the modulus, between 0 and 1.
!
!    Input, real ( kind = 8 ) PHI, the argument in degrees.
!
!    Output, real ( kind = 8 ) FE, EE, the values of F(k,phi) and E(k,phi).
		*/

		void v4_elit_pd(const Vec4 hk,
			        const Vec4 phi,
			        Vec4 &fe,
			        Vec4 &ee) {
#if defined _WIN64		  
			__declspec(align(64)) struct {
			  Vec4 a,a0,b,b0,c,c0,ck,d,d0,fac,g,r;
			  int32_t n;
#if (USE_STRUCT_PADDING) == 1
			  PAD_TO(0,28)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                          Vec4 a,a0,b,b0,c,c0,ck,d,d0,fac,g,r;
			  int32_t n;
#if (USE_STRUCT_PADDING) == 1
			  PAD_TO(0,28)
#endif
			}ca;
#endif
			ca.a    = v4_n0;
			ca.a0   = v4_n0;
			ca.b    = v4_n0;
			ca.b0   = v4_n0;
			ca.c    = v4_n0;
			ca.c0   = v4_n0;
			ca.ck   = v4_n0;
			ca.d    = v4_n0;
			ca.d0   = v4_n0;
			ca.fac  = v4_n0;
			ca.g    = v4_n0;
			ca.r    = v4_n0;
			
			ca.a0 = v4_n1;
			ca.b0 = Vec4_SQRT(Vec4_SUB(v4_n1, Vec4_MUL(hk,hk)));
			ca.d0 = Vec4_MUL(Vec4_MUL(v4PI,v4_1over180),phi);
			ca.r = Vec4_MUL(hk,hk);
			const Vec4 mask_eq1  = Vec4_CMP(hk,v4_n1,_CMP_EQ_OQ);
			const Vec4 mask_eq90 = Vec4_CMP(phi, Vec4_SET1(90.0),_CMP_EQ_OQ);
			
			if (!Vec4_TESTZ(mask_eq1, mask_eq1) &&
				!Vec4_TESTZ(mask_eq90, mask_eq90)) {
				fe = v4HUGE;
				ee = v4_n1;
				return;
			
		   }
		   else if (!Vec4_TESTZ(mask_eq1, mask_eq1)) {
			    const Vec4 cosd0 = Vec4_COS(ca.d0);
			    fe = Vec4_LOG(Vec4_DIV(Vec4_ADD(v4_n1, Vec4_SIN(ca.d0)),cosd0));
			    ee = Vec4_SIN(ca.d0);
		  }
		  else {
			   ca.fac = v4_n1;
			   const Vec4 eps = Vec4SET1(1.0E-7);
			   for (ca.n = 1; ca.n != 40; ++ca.n) {
				     ca.a = Vec4_MUL(Vec4_ADD(ca.a0,ca.b0),v4_1over2);
				     ca.b = Vec4_SQRT(Vec4_MUL(ca.a0,ca.b0));
				     ca.c = Vec4_MUL(Vec4_SUB(ca.a0,ca.b0),v4_1over2);
				     ca.fac = Vec4_ADD(v4_n2,ca.fac);
				     ca.r = Vec4_MUL(Vec4_ADD(ca.r, ca.fac), Vec4_MUL(ca.c,ca.c));
				     const Vec4 mask_neq90 = Vec4_CMP(phi, Vec4_SET1(90.0),_CMP_NEQ_OQ);
				     if (!Vec4_TESTZ(mask_neq90, mask_neq90)) {
					     ca.d = Vec4_ADD(ca.d0, Vec4_ATAN(Vec4_MUL(Vec4_DIV(ca.b0, ca.a0), Vec4_TAN(ca.d0))));
					     ca.g = Vec4_FMAD(ca.c, Vec4_SIN(ca.d),ca.g);
					     const  Vec4 rhs = Vec4_ADD(Vec4_DIV(ca.d,v4PI),v4_1over2); // int ( d / pi + 0.5D+00 )
					     ca.d0 = Vec4_MUL(Vec4_ADD(ca.d,v4PI),rhs);
				     }
				     ca.a0 = ca.a;
				     ca.b0 = ca.b;
			    }
			       const Vec4 mask_lteps = Vec4_CMP(ca.c,eps,_CMP_LT_OQ);
			       if (!Vec4_TESTZ(mask_lteps, mask_lteps)) {
				        break;
			      }
			      ca.ck = Vec4_DIV(v4PI, Vec4_MUL(v4_n2,ca.a));
			      ca.ce = Vec4_MUL(v4PI, Vec4_DIV(Vec4_SUB(v4_n2, ca.r), Vec4_MUL(v4_n4,ca.a)));
				  if (!Vec4_TESTZ(mask_eq90, mask_eq90)) {
					  fe = ca.ck;
					  ee = ca.ce;
				  }
				  else {
					  fe = Vec4_DIV(ca.d, Vec4_MUL(ca.fac,ca.a));
					  ee = Vec4_DIV(Vec4_MUL(fe, ca.ce), Vec4_ADD(ca.ck,ca.g));
				  }
		      }

	      }

		void v8_elit_pd(const Vec8 hk,
			        const Vec8 phi,
			        Vec8 &fe,
			        Vec8 &ee) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 a,a0,b,b0,c,ce,ck,d,d0,fac,g,r,t0;
				int32_t n;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec8 a,a0,b,b0,c,ce,ck,d,d0,fac,g,r,t0;
				int32_t n;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#endif
			ca.a = v8_n0;
			ca.a0 = v8_n0;
			ca.b = v8_n0;
			ca.b0 = v8_n0;
			ca.c  = v8_n0;
			ca.ce = v8_n0;
			ca.ck = v8_n0;
			ca.d  = v8_n0;
			ca.d0 = v8_n0;
			ca.fac = v8_n0;
			ca.g = v8_n0;
			ca.r = v8_n0;
			
			ca.a0 = v8_n1;
			ca.t0 = Vec8_MUL(hk,hk);
			ca.b0 = Vec8_SQRT(Vec8_SUB(v8_n1,ca.t0));
			ca.d0 = Vec8_MUL(Vec8_MUL(v8PI,v8_1over180),phi);
			ca.r = ca.t0;
			if (Vec8_CMP(hk, v8_n1, _CMP_EQ_OQ) &&
				Vec8_CMP(phi, Vec8_SET1(90.0), _CMP_EQ_OQ)) {
				fe = v8HUGEP;
				ee = v8_n1;
				return;
			}
			else if (Vec8_CMP(hk, v8_n1, _CMP_EQ_OQ)) {
				const Vec8 sind0 = Vec8_SIN(ca.d0);
				fe = Vec8_LOG(Vec8_DIV(Vec8_ADD(v8_n1, sind0), Vec8_COS(ca.d0)));
				ee = sind0;
			}
			else {
				ca.fac = v8_n1;
				const Vec8 eps = Vec8SET1(1.0E-7);
				for (ca.n = 1; ca.n != 40; ++ca.n) {
					ca.a = Vec8_MUL(Vec8_ADD(ca.a0, ca.b0), v8_1over2);
					ca.b = Vec8_SQRT(Vec8_MUL(ca.a0, ca.b0));
					ca.c = Vec8_MUL(Vec8_SUB(ca.a0, ca.b0), v8_1over2);
					ca.fac = Vec8_ADD(v8_n2, ca.fac);
					ca.r = Vec8_MUL(Vec8_ADD(ca.r, ca.fac), Vec8_MUL(ca.c, ca.c));
					
					if (Vec8_CMP(phi, Vec8_SET1(90.0), _CMP_NEQ_OQ)) {
						ca.d = Vec8_ADD(ca.d0, Vec8_ATAN(Vec8_MUL(Vec4_DIV(ca.b0, ca.a0), Vec8_TAN(ca.d0))));
						ca.g = Vec8_FMAD(ca.c, Vec8_SIN(ca.d), ca.g);
						const  Vec8 rhs = Vec8_ADD(Vec8_DIV(ca.d, v8PI), v8_1over2); // int ( d / pi + 0.5D+00 )
						ca.d0 = Vec8_MUL(Vec8_ADD(ca.d, v8PI), rhs);
					}
					ca.a0 = ca.a;
					ca.b0 = ca.b;
				}
				
				if (Vec8_CMP(ca.c, eps, _CMP_LT_OQ)) {
					break;
				}
				ca.ck = Vec8_DIV(v8PI, Vec8_MUL(v8_n2, ca.a));
				ca.ce = Vec8_MUL(v8PI, Vec8_DIV(Vec8_SUB(v8_n2, ca.r), Vec8_MUL(v8_n4, ca.a)));
				if (Vec8_CMP(phi, Vec4_SET1(90.0), _CMP_NEQ_OQ)) {
					fe = ca.ck;
					ee = ca.ce;
				}
				else {
					fe = Vec8_DIV(ca.d, Vec8_MUL(ca.fac, ca.a));
					ee = Vec8_DIV(Vec8_MUL(fe, ca.ce), Vec8_ADD(ca.ck, ca.g));
				}
			}
		}

		/*
			 ELIT3 computes the elliptic integral of the third kind.
!
!  Discussion:
!
!    Gauss-Legendre quadrature is used.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 July 2012
!    07 October 2018 Bernard Gingold

!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PHI, the argument in degrees.
!
!    Input, real ( kind = 8 ) HK, the modulus, between 0 and 1.
!
!    Input, real ( kind = 8 ) C, the parameter, between 0 and 1.
!
!    Output, real ( kind = 8 ) EL3, the value of the elliptic integral
!    of the third kind.
		*/
		 
		v4_elit3_pd(const Vec4 phi,
			    const Vec4 hk,
			    const Vec4 c,
			    Vec4 &el3) {
#if defined _WIN64		  
			__declspec(align(64)) struct {
				Vec4 c0,c1,c2,f1,f2,t1,t2;
				int32_t i;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align))) struct {
                                Vec4 c0,c1,c2,f1,f2,t1,t2;
				int32_t i;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,28)
#endif
			}ca;
#endif
#if defined _WIN64
			const __declspec(align(64)) double t[10] = { 
			    0.9931285991850949E+00, 0.9639719272779138E+00,
				0.9122344282513259E+00, 0.8391169718222188E+00,
				0.7463319064601508E+00, 0.6360536807265150E+00, 
				0.5108670019508271E+00, 0.3737060887154195E+00,
				0.2277858511416451E+00, 0.7652652113349734E-01 };
			const __declspec(align(64)) double w[10] = {
				0.1761400713915212E-01, 0.4060142980038694E-01,
				0.6267204833410907E-01, 0.8327674157670475E-01, 
				0.1019301198172404E+00, 0.1181945319615184E+00, 
				0.1316886384491766E+00, 0.1420961093183820E+00, 
				0.1491729864726037E+00, 0.1527533871307258E+00};
#elif defined
			const __attribute__((align(64))) double t[10] = { 
			        0.9931285991850949E+00, 0.9639719272779138E+00,
				0.9122344282513259E+00, 0.8391169718222188E+00,
				0.7463319064601508E+00, 0.6360536807265150E+00, 
				0.5108670019508271E+00, 0.3737060887154195E+00,
				0.2277858511416451E+00, 0.7652652113349734E-01 };
			const __attribute__((align(64))) double w[10] = {
				0.1761400713915212E-01, 0.4060142980038694E-01,
				0.6267204833410907E-01, 0.8327674157670475E-01, 
				0.1019301198172404E+00, 0.1181945319615184E+00, 
				0.1316886384491766E+00, 0.1420961093183820E+00, 
				0.1491729864726037E+00, 0.1527533871307258E+00};
#endif
			ca.c0 = v4_n0;
			ca.c1 = v4_n0;
			ca.c2 = v4_n0;
			ca.f1 = v4_n0;
			ca.f2 = v4_n0;
			ca.t1 = v4_n0;
			ca.t2 = v4_n0;
			const Vec4 eps = Vec4_SET1(1.0E-8);
			const Vec4 mask_hkeq1 = Vec4_CMP(hk,v4_n1,_CMP_EQ_OQ);
			const Vec4 mask_philteps = Vec4_CMP(abs(Vec4_SUB(phi, Vec4_SET1(90.0))),eps,_CMP_LE_OQ);
			const Vec4 mask_ceq1 = Vec4_CMP(c,v4_n1,_CMP_EQ_OQ);
			if ((!Vec4_TESTZ(mask_hkeq1, mask_hkeq1) &&
				!Vec4_TESTZ(mask_philteps, mask_philteps)) ||
				(!Vec4_TESTZ(mask_ceq1, mask_ceq1) &&
				!Vec4_TESTZ(mask_philteps, mask_philteps))) {
				el3 = v4HUGE;
				return;
			}
			ca.c1 = Vec4_MUL(Vec4_SET1(0.87266462599716E-02),phi);
			ca.c2 = ca.c1;
			el3 = v4_n0;
			const Vec4 hkhk = Vec4_MUL(hk,hk);
			for (ca.i = 0; ca.i != 10; ++ca.i) {
				ca.c0 = Vec4_MUL(ca.c2, Vec4_SET1(t[i]));
				ca.t1 = Vec4_ADD(ca.c1,ca.c0);
				ca.t2 = Vec4_SUB(ca.c1,ca.c0);
				const Vec4 sint1 = Vec4_SIN(ca.t1);
				const Vec4 lhs = Vec4_SUB(v4_n1, Vec4_MUL(ca.c, Vec4_MUL(sint1,sint1)));
				const Vec4 rhs = Vec4_SUB(v4_n1, Vec4_MUL(hkhk, Vec4_MUL(sint1,sint1)));
				ca.f1 = Vec4_DIV(v4_n1, Vec4_MUL(lhs, Vec4_SQRT(rhs)));
				const Vec4 sint2 = Vec4_SIN(ca.t2);
				const Vec4 lhs2 = Vec4_SUB(v4_n1, Vec4_MUL(ca.c, Vec4_MUL(sint2, sint2)));
				const Vec4 rhs2 = Vec4_SUB(v4_n1, Vec4_MUL(hkhk, Vec4_MUL(sint2, sint2)));
				ca.f2 = Vec4_DIV(v4_n1, Vec4_MUL(lhs2, Vec4_SQRT(rhs2)));
				el3 = Vec4_ADD(el3, Vec4_MUL(Vec4_SET1(w[i]), Vec4_ADD(ca.f1,ca.f2)));
			}
			el3 = Vec4_MUL(ca.c1,el3);
		}

		v8_elit3_pd(const Vec8 phi,
			    const Vec8 hk,
		            const Vec8 c,
			    Vec8 &el3) {
#if defined _WIN64
			__declspec(align(64)) struct {
				Vec8 c,c0,c1,c2,f1,f2,t1,t2;
				int32_t i;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec8 c,c0,c1,c2,f1,f2,t1,t2;
				int32_t i;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,60)
#endif
			}ca;
#endif
#if defined _WIN64			  
			const __declspec(align(64)) double t[10] = {
				0.9931285991850949E+00, 0.9639719272779138E+00,
				0.9122344282513259E+00, 0.8391169718222188E+00,
				0.7463319064601508E+00, 0.6360536807265150E+00,
				0.5108670019508271E+00, 0.3737060887154195E+00,
				0.2277858511416451E+00, 0.7652652113349734E-01 };
			const __declspec(align(64)) double w[10] = {
				0.1761400713915212E-01, 0.4060142980038694E-01,
				0.6267204833410907E-01, 0.8327674157670475E-01,
				0.1019301198172404E+00, 0.1181945319615184E+00,
				0.1316886384491766E+00, 0.1420961093183820E+00,
				0.1491729864726037E+00, 0.1527533871307258E+00 };
#elif defined __linux
			  const __attribute__((align(64))) double t[10] = {
				0.9931285991850949E+00, 0.9639719272779138E+00,
				0.9122344282513259E+00, 0.8391169718222188E+00,
				0.7463319064601508E+00, 0.6360536807265150E+00,
				0.5108670019508271E+00, 0.3737060887154195E+00,
				0.2277858511416451E+00, 0.7652652113349734E-01 };
			  const __attribute__((align(64))) double w[10] = {
				0.1761400713915212E-01, 0.4060142980038694E-01,
				0.6267204833410907E-01, 0.8327674157670475E-01,
				0.1019301198172404E+00, 0.1181945319615184E+00,
				0.1316886384491766E+00, 0.1420961093183820E+00,
				0.1491729864726037E+00, 0.1527533871307258E+00 };
#endif
				ca.c  = v8_n0;
				ca.c0 = v8_n0;
				ca.c1 = v8_n0;
				ca.c2 = v8_n0;
				ca.f1 = v8_n0;
				ca.f2 = v8_n0;
				ca.t1 = v8_n0;
				ca.t2 = v8_n0;
				const Vec8 eps = Vec8_SET1(1.0E-8);
				if ((Vec8_CMP(hk, v4_n1, _CMP_EQ_OQ) &&
					 Vec8_CMP(Vec8_ABS(Vec8_SUB(phi, Vec8_SET1(90.0))), eps, _CMP_LE_OQ)) ||
					(Vec8_CMP(c, v8_n1, _CMP_EQ_OQ) &&
					 Vec8_CMP(Vec8_ABS(Vec8_SUB(phi, Vec8_SET1(90.0))), eps, _CMP_LE_OQ))) {
					 el3 = v8HUGEP;
				}
				ca.c1 = Vec8_MUL(Vec8_SET1(0.87266462599716E-02), phi);
				ca.c2 = ca.c1;
				el3 = v8_n0;
				const Vec8 hkhk = Vec8_MUL(hk, hk);
				for (ca.i = 0; ca.i != 10; ++ca.i) {
					ca.c0 = Vec8_MUL(ca.c2, Vec8_SET1(t[i]));
					ca.t1 = Vec8_ADD(ca.c1, ca.c0);
					ca.t2 = Vec8_SUB(ca.c1, ca.c0);
					const Vec8 sint1 = Vec8_SIN(ca.t1);
					const Vec8 lhs = Vec8_SUB(v8_n1, Vec8_MUL(ca.c, Vec8_MUL(sint1, sint1)));
					const Vec8 rhs = Vec8_SUB(v8_n1, Vec8_MUL(hkhk, Vec8_MUL(sint1, sint1)));
					ca.f1 = Vec8_DIV(v8_n1, Vec8_MUL(lhs, Vec8_SQRT(rhs)));
					const Vec8 sint2 = Vec8_SIN(ca.t2);
					const Vec8 lhs2 = Vec8_SUB(v8_n1, Vec8_MUL(ca.c, Vec8_MUL(sint2, sint2)));
					const Vec8 rhs2 = Vec8_SUB(v8_n1, Vec8_MUL(hkhk, Vec8_MUL(sint2, sint2)));
					ca.f2 = Vec8_DIV(v8_n1, Vec8_MUL(lhs2, Vec8_SQRT(rhs2)));
					el3 = Vec8_ADD(el3, Vec8_MUL(Vec8_SET1(w[i]), Vec8_ADD(ca.f1, ca.f2)));
				}
				el3 = Vec8_MUL(ca.c1, el3);
		}

		/*
			JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!	 07 October 2018 Bernard Gingold
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
		*/

		void v4_jy01a_pd(const Vec4 x,
			         Vec4 &bj0,
			         Vec4 &dj0,
			         Vec4 &bj1,
			         Vec4 &dj1,
			         Vec4 &by0,
			         Vec4 &dy0,
			         Vec4 &by1,
			         Vec4 &dy1) {
#if defined _WIN64		  
			__declspec(align(64)) struct {
				Vec4 cs0,cs1,cu,ec,p0,p1,q0,q1,r,r0,
					rp2,t1,t2,w0,w1,x2;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,4)
#endif
			    int32_t k0;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,4)
#endif
#if (USE_STRUCT_PADDING) == 1
			    PAD_TO(2,16)
#endif
			}ca;
#elif defined __linux
                        __attribute__((align(64))) struct {
                                Vec4 cs0,cs1,cu,ec,p0,p1,q0,q1,r,r0,
					rp2,t1,t2,w0,w1,x2;
				int32_t k;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(0,4)
#endif
			    int32_t k0;
#if (USE_STRUCT_PADDING) == 1
				PAD_TO(1,4)
#endif
#if (USE_STRUCT_PADDING) == 1
			    PAD_TO(2,16)
#endif
			}ca;
#endif
#if defined _WIN64
			__declspec(align(64)) const double a[12] = {
				-0.7031250000000000E-01, 0.1121520996093750E+00,
				-0.5725014209747314E+00, 0.6074042001273483E+01,
				-0.1100171402692467E+03, 0.3038090510922384E+04,
				-0.1188384262567832E+06, 0.6252951493434797E+07,
				-0.4259392165047669E+09, 0.3646840080706556E+11,
				-0.3833534661393944E+13, 0.4854014686852901E+15};
			__declspec(align(64)) const double a1[12] = {
				0.1171875000000000E+00, -0.1441955566406250E+00,
				0.6765925884246826E+00, -0.6883914268109947E+01, 
				0.1215978918765359E+03, -0.3302272294480852E+04, 
				0.1276412726461746E+06, -0.6656367718817688E+07, 
				0.4502786003050393E+09, -0.3833857520742790E+11, 
				0.4011838599133198E+13, -0.5060568503314727E+15};
			__declspec(align(64)) const double b[12] = {
				0.7324218750000000E-01, -0.2271080017089844E+00, 
				0.1727727502584457E+01, -0.2438052969955606E+02, 
				0.5513358961220206E+03,  -0.1825775547429318E+05, 
				0.8328593040162893E+06, -0.5006958953198893E+08, 
				0.3836255180230433E+10, -0.3649010818849833E+12, 
				0.4218971570284096E+14, -0.5827244631566907E+16 };
			__declspec(align(64)) const double b1[12] = {
				-0.1025390625000000E+00, 0.2775764465332031E+00, 
				-0.1993531733751297E+01, 0.2724882731126854E+02, 
				-0.6038440767050702E+03, 0.1971837591223663E+05, 
				-0.8902978767070678E+06, 0.5310411010968522E+08, 
				-0.4043620325107754E+10, 0.3827011346598605E+12, 
				-0.4406481417852278E+14, 0.6065091351222699E+16 };
#elif defined __linux
			 __attribute__((align(64)))   const double a[12] = {
				-0.7031250000000000E-01, 0.1121520996093750E+00,
				-0.5725014209747314E+00, 0.6074042001273483E+01,
				-0.1100171402692467E+03, 0.3038090510922384E+04,
				-0.1188384262567832E+06, 0.6252951493434797E+07,
				-0.4259392165047669E+09, 0.3646840080706556E+11,
				-0.3833534661393944E+13, 0.4854014686852901E+15};
			 __attribute__((align(64)))  const double a1[12] = {
				0.1171875000000000E+00, -0.1441955566406250E+00,
				0.6765925884246826E+00, -0.6883914268109947E+01, 
				0.1215978918765359E+03, -0.3302272294480852E+04, 
				0.1276412726461746E+06, -0.6656367718817688E+07, 
				0.4502786003050393E+09, -0.3833857520742790E+11, 
				0.4011838599133198E+13, -0.5060568503314727E+15};
			  __attribute__((align(64))) const double b[12] = {
				0.7324218750000000E-01, -0.2271080017089844E+00, 
				0.1727727502584457E+01, -0.2438052969955606E+02, 
				0.5513358961220206E+03,  -0.1825775547429318E+05, 
				0.8328593040162893E+06, -0.5006958953198893E+08, 
				0.3836255180230433E+10, -0.3649010818849833E+12, 
				0.4218971570284096E+14, -0.5827244631566907E+16 };
		          __attribute__((align(64))) const double b1[12] = {
				-0.1025390625000000E+00, 0.2775764465332031E+00, 
				-0.1993531733751297E+01, 0.2724882731126854E+02, 
				-0.6038440767050702E+03, 0.1971837591223663E+05, 
				-0.8902978767070678E+06, 0.5310411010968522E+08, 
				-0.4043620325107754E+10, 0.3827011346598605E+12, 
				-0.4406481417852278E+14, 0.6065091351222699E+16 };
#endif
			ca.cs0 = v4_n0;
			ca.cs1 = v4_n0;
			ca.cu  = v4_n0;
			ca.ec  = v4_n0;
			ca.p0  = v4_n0;
			ca.p1  = v4_n0;
			ca.q0  = v4_n0;
			ca.q1  = v4_n0;
			ca.r   = v4_n0;
			ca.r0  = v4_n0;
			ca.rp2 = Vec4_SET1(0.63661977236758);
			ca.t1  = v4_n0;
			ca.t2  = v4_n0;
			ca.w0  = v4_n0;
			ca.w1  = v4_n0;
			ca.x2 = Vec4_MUL(x,x);
			const Vec4 mask_xeq0 = Vec4_CMP(x,v4_n0,_CMP_EQ_OQ);
			if (!Vec4_TESTZ(mask_xeq0, mask_xeq0)) {
				bj0 = v4_n1;
				bj1 = v4_n0;
				dj0 = v4_n0;
				dj1 = v4_1over2;
				by0 = v4HUGE;
				by1 = v4NHUGE;
				dy0 = v4HUGE;
				dy1 = v4HUGE;
				return;
			}
			const Vec4 mask_xle12 = Vec4_CMP(x,v4_n12,_CMP_LE_OQ);
			if (!Vec4_TESTZ(mask_xle12, mask_xle12)) {
				bj0 = v4_n1;
				ca.r = v4_n1;
				for (ca.k = 1; ca.k != 31; ++ca.k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));
					const Vec4 kk = Vec4_MUL(tk,tk);
					ca.r = Vec4_DIV(Vec4_MUL(Vec4_MUL(flip_sign(v4_1over4),ca.r),x2),kk);
					bj0 = Vec4_ADD(bj0,ca.r);
					const Vec4 mask_rltbj0 = Vec4_CMP(abs(ca.r), Vec4_MUL(abs(bj0),v4EPS1),_CMP_LT_OQ);
					if (!Vec4_TESTZ(mask_rltbj0, mask_rltbj0)) {
						break;
					}
				}
				bj1 = v4_n1;
				ca.r = v4_n1;
				for (ca.k = 1; ca.k != 31; ++ca.k) {
					const Vec4 tk = Vec4_CVTI4(Vec4_SETI4(ca.k));

				}
			}
		}

#include "GMS_avx512c8f64.h"


		void v8_cerf_pd(const AVX512c8f64 z,
				AVX512c8f64 &cer,
				AVX512c8f64 &cder ) {
#if defined _WIN64
		  __declspec(align(64)) struct {
		    AVX512c8f64 c0,cs; // 4*64 L1D lines
		  }dc; // data cold  
#elif defined __linux
		  __attribute__((align(64))) struct {
                    AVX512c8f64 c0,cs; // 4*64 L1D lines
		  }dc;  
#endif
#if defined _WIN64
		  __declspec(align(64)) struct {
                    Vec8 ei1,ei2; // 2*64 L1D
		    Vec8 er,er0;  // 2*64 L1D
		    Vec8 er1,er2; // 2*64 L1D
		    Vec8 eri,err; // 2*64 L1D
		    Vec8 r,ss;    // 2*64 L1D
		    Vec8 w,w1;    // 2*64 L1D
		    Vec8 w2,x;    // 2*64 L1D
		    Vec8 x2,y;    // 2*64 L1D
		    Vec8 expx2;
		    int32_t k,n;
#if (USE_STRUCT_PADDING) == 1
		    PAD_TO(0,58)
#endif
		  }dh; // data hot
#elif defined __linux
		  __attribute__((align(64))) struct {
                    Vec8 ei1,ei2; // 2*64 L1D
		    Vec8 er,er0;  // 2*64 L1D
		    Vec8 er1,er2; // 2*64 L1D
		    Vec8 eri,err; // 2*64 L1D
		    Vec8 r,ss;    // 2*64 L1D
		    Vec8 w,w1;    // 2*64 L1D
		    Vec8 w2,x;    // 2*64 L1D
		    Vec8 x2,y;    // 2*64 L1D
		    Vec8 expx2;
		    int32_t k,n;
#if (USE_STRUCT_PADDING) == 1
		    PAD_TO(0,58)
#endif
		  }dh; // data hot
#endif
		  // first touch
		  dh.ei1   = v8_n0;
		  dh.ei2   = v8_n0;
		  dh.er    = v8_n0;
		  dh.er0   = v8_n0;
		  dh.er1   = v8_n0;
		  dh.er2   = v8_n0;
		  dh.eri   = v8_n0;
		  dh.err   = v8_n0;
		  dh.r     = v8_n0;
		  dh.ss    = v8_n0;
		  dh.w     = v8_n0;
		  dh.w1    = v8_n0;
		  dh.w2    = v8_n0;
		  dh.x     = v8_n0;
		  dh.y     = v8_n0;
		  dh.x2    = v8_n0;
		  dh.expx2 = v8_n0;
		  dh.x = z.m_re;
		  dh.y = z.m_im;
		  dh.x2 = Vec8_MUL(dh.x,dh.x);
		  dh.expx2 = Vec8_EXP(flip_sign(dh.x2));
		  if (Vec8_CMP(dh.x,Vec8_SET1(3.5),_CMP_LE_OQ)) {
		      dh.er = v8_n1;
		      dh.r  = v8_n1;
		      dh.k = 0;
		      for(dh.k = 1; dh.k != 100; ++dh.k) {
                          const Vec8 vk = Vec8_CVTI4(Vec8_SETI4(dh.k));
			  dh.r = Vec8_DIV(Vec8_MUL(dh.r,dh.x2),Vec8_ADD(vk,v8_1over2));
			  dh.er = Vec8_ADD(dh.er,dh.r);
			  const vcmp1 = Vec8_ABS(Vec8_SUB(dh.er,dh.w));
			  const vcmp2 = Vec8_MUL(v8EPS1ton12,Vec8_ABS(dh.er));
			  if (Vec8_CMP(vcmp1,vcmp2,_CMP_LE_OQ)) {
                              return;
			  }
			  dh.w = dh.er;
		      }
		      dh.c0 = AVX512c8f64{};
		      dc.c0.m_re = Vec8_MUL(Vec8_SET1(1.128379167),Vec8_MUL(dh.x,dh.expx2)));
		      dh.er0     = Vec8_MUL(dc.c0.m_re,dh.er);
		}
		 else {
                      dh.er = v8_n1;
		      dh.r  = v8_n1;
		      dc.c0 = AVX512c8f64{};
		      for(dh.k = 1; dh.k != 12; ++dh.k) {
                          const Vec8 vk = Vec8_CVTI4(Vec8_SETI4(dh.k));
			  dh.r = Vec8_DIV(Vec8_MUL(flip_sign(dh.r),Vec8_SUB(vk,v8_1over2)),dh.x2);
			  dh.er = Vec8_ADD(dh.er,dh.r);
		      }
		      dc.c0.m_re = Vec8_DIV(dh.expx2,Vec8_MUL(dh.x,Vec8_SET1(1.128379167)));
		      dh.er0     = Vec8_SUB(v8_n1,Vec8_MUL(dc.c0.m_re,dh.er));
		 }

		 if (Vec8_CMP(dh.y,v8_n0,CMP_EQ_OQ)) {
                      dh.err = dh.er0;
		      dh.eri = v8_n0;
		 }
		  else {
                      dc.cs = AVX512c8f64{};
		      dc.cs.m_re = Vec8_COS(Vec8_MUL(v8_n2,Vec8_MUL(dh.x,dh.y)));
		      dh.ss = Vec8_SIN(Vec8_MUL(v8_n2,Vec8_MUL(dh.x,dh.y)));
		      dh.er1 = Vec8_MUL(dh.expx2,Vec8_DIV(Vec8_SUB(v8_n1,dc.cs.m_re),
		                       Vec8_MUL(Vec8_SET1(6.28318530718),dh.x)));
		      dh.ei1 = Vec8_MUL(dh.expx2,Vec8_DIV(dh.ss,
		                             Vec8_MUL(Vec8_SET1(6.28318530718),dh.x)));
		      dh.er2 = v8_n0;
		      dh.n = 0;
		      for(dh.n = 1; dh.n != 100; ++dh.n) {
                         const Vec8 vn     = Vec8_CVTI4(Vec8_SETI4(dh.n));
			 const Vec8 arg    = Vec8_MUL(vn,dh.y);
			 const Vec8 exp    = Vec8_EXP(Vec8_MUL(flip_sign(v8_1over4),Vec8_MUL(vn,vn)));
			 const Vec8 t0     = Vec8_ADD(Vec8_MUL(vn,vn),Vec8_MUL(v8_n4,dh.x2));
		         const Vec8 t1     = Vec8_MUL(v8_n2,dh.x);
			 const Vec8 t2     = Vec8_MUL(t1,Vec8_MUL(Vec8_COSH(arg),dc.cs.m_re));
			 const Vec8 t3     = Vec8_MUL(vn,Vec8_MUL(Vec8_SINH(arg),dh.ss));
			 const Vec8 t4     = Vec8_DIV(exp,t0);
			 const Vec8 t5     = Vec8_SUB(t1,t2);
			 const Vec8 t6     = Vec8_ADD(t5,t3);
			 dh.er2 = Vec8_ADD(dh.er2,Vec8_MUL(t4,t6));
			 const Vec8 vcmp1 = Vec8_DIV(Vec8_SUB(dh.er2,dh.w1),dh.er2);
			 if(Vec8_CMP(Vec8_ABS(vcmp1),v8EPS1ton12,_CMP_LT_OQ)) {
                             return;
			 }
			 dh.w1 = dh.er2;
		      }

		      dc.c0 = AVX512c8f64{};
		      dc.c0.m_re = Vec8_MUL(v8_n2,Vec8_MUL(dh.expx2,v8PIR));
		      dh.err = Vec8_ADD(dh.er0,Vec8_FMAD(dc.c0.m_re,dh.er2,dh.er1));
		      dh.ei2 = v8_n0;
		      for(dh.n = 1; dh.n != 100; ++dh.n) {
                           const Vec8 vn     = Vec8_CVTI4(Vec8_SETI4(dh.n));
			   const Vec8 arg    = Vec8_MUL(vn,dh.y);
			   const Vec8 exp    = Vec8_EXP(Vec8_MUL(flip_sign(v8_1over4),Vec8_MUL(vn,vn)));
			   const Vec8 t0     = Vec8_ADD(Vec8_MUL(vn,vn),Vec8_MUL(v8_n4,dh.x2));
			   const Vec8 t1     = Vec8_MUL(Vec8_MUL(v8_n2,dh.x),Vec8_MUL(Vec8_COSH(arg),dh.ss));
			   const Vec8 t2     = Vec8_MUL(vn,Vec8_MUL(Vec8_SINH(arg),dc.cs.m_re));
			   const Vec8 t3     = Vec8_DIV(exp,t0);
			   dh.ei2 = Vec8_ADD(dh.ei2,Vec8_FMAD(t3,t1,t2));
			   const Vec8 vcmp1 = Vec8_DIV(Vec8_SUB(dh.ei2,dh.w2),dh.ei2);
			   if(Vec8_CMP(Vec8_ABS(vcmp1,v8EPS1ton12,_CMP_LT_OQ)) {
                                 return;
			   }
			   dh.w2 = dh.ei2;
		      }
		  }
		  dh.eri = Vec8_FMAD(cs.c0.m_re,dh.ei2,dh.ei1);
	     }
	     cer  = AVX512c8f64{dh.err,dh.eri};
	     const Vec8 tmp = Vec8_SET1(1.128379167095512573896);
	     cder = tmp * cexp(~z*z);
      }

      void v8_cerror_pd(const AVX512c8f64 z,
                        AVX512c8f64 &cer ) {
#if defined _WIN64
              __declspec(align(64)) struct {
                   AVX512c8f64 c0,cl; // 4*64 L1D 
		   AVX512c8f64 cr,cs; // 4*64 l1D
		   AVX512c8f64 z1;    // 4*64 l1D
		   Vec8 a0;           // 1*64 l1D
		   int32_t k;         // 1/16 l1d
#if defined (USE_STRUCT_PADDING) == 1
                   PAD_TO(0,60)
#endif
	      }dh;
#elif defined __linux
              __attribute__((align(64))) struct {
                   AVX512c8f64 c0,cl; // 4*64 L1D 
		   AVX512c8f64 cr,cs; // 4*64 l1D
		   AVX512c8f64 z1;    // 4*64 l1D
		   Vec8 a0;           // 1*64 l1D
		   int32_t k;         // 1/16 l1d
#if defined (USE_STRUCT_PADDING) == 1
                   PAD_TO(0,60)
#endif
	      }dh;
#endif
              // First touch
	      dh.c0 = AVX512c8f64{};
	      dh.cl = AVX512c8f64{};
	      dh.cr = AVX512c8f64{};
	      dh.cs = AVX512c8f64{};
	      dh.z1 = AVX512c8f64{};
	      dh.a0 = v8_n0;
              dh.a0 = Vec8_ABS(z.m_re);
	      dh.c0 = cexp(~z*z);
	      dh.z1 = z;
	      if (Vec8_CMP(z.m_re,v8_n0,_CMP_LT_OQ)) {
                   dh.z1 = ~z;
	      }
	      if (Vec8_CMP(dh.a0,Vec8_SET1(5.8),_CMP_LE_OQ)) {
                   dh.cs = dh.z1;
		   dh.cr = dh.z1;
		   dh.k = 0;
		   for (dh.k = 1; dh.k != 120; ++dh.k) {
                        const Vec8 vk = Vec8_CVTI4(Vec8_SETI4(dh.k));
			dh.cr = dh.cr * dh.z1 * dh.z1 / Vec8_ADD(vk,v8_1over2);
			dh.cs = dh.cs + dh.cr;
			const AVX512c8f64 tmp = cabs(dh.cr/dh.cs);
			if (Vec8_CMP(tmp.m_re,Vec8_SET1(0.000000000000001),_CMP_LT_OQ)) {
                             return;
			}
		   }
		   cer = v8_n2 * dh.c0 * dh.cs * Vec8_SET1(0.564189583547756286948);
	      }
	       else {
                   dh.cl = v8_n1 / dh.z1;
		   dh.cr = dh.cl;
		   dh.k = 0;
		   for (dh.k = 1; dh.k != 13; ++dh.k) {
                        const Vec8 vk = Vec8_CVTI4(Vec8_SETI4(dh.k));
			dh.cr = ~dh.cr * Vec8_SUB(vk,v8_1over2) / (dh.z1*dh.z1);
			dh.cl = dh.cl + dh.cr;
			const AVX512c8f64 tmp = cabs(dh.cr/dh.cl);
			if (Vec8_CMP(tmp.m_re,Vec8_SET1(0.000000000000001),_CMP_LT_OQ)) {
                             return;
			}
		   }
		   cer = v8_n1 - dh.c0 * dh.cl / Vec8_SET1(0.564189583547756286948);
	       }

	       if (Vec8_CMP(z.m_re,v8_n0,_CMP_LT_OQ)) {
                   cer = ~cer;
	       }
      }

      template<int32_t nt>
           void v8_cerzo_pd(AVX512c8f64 __restrict zo[nt]) {
#if defined _WIN64
              __declspec(align(64)) struct {
                 Vec8 pu; // 1 * 64 L1D line
		 Vec8 pv; // 1 * 64 L1D line
		 Vec8 px; // 1 * 64 L1D line
		 Vec8 py; // 1 * 64 L1D line
		 Vec8 w;  // 1 * 64 L1D line
		 Vec8 w0; // 1 * 64 L1D line
	      }dc; // data cold
#elif defined __linux
              __attribute__((align(64))) struct {
                 Vec8 pu; // 1 * L1D (assumed single L1D cache line (64 bytes)
		 Vec8 pv; // 1 * L1D
		 Vec8 px; // 1 * L1D
		 Vec8 py; // 1 * L1D
		 Vec8 w;  // 1 * L1D
		 Vec8 w0; // 1 * L1D
	      }dc;
#endif
#if defined _WIN64
             __declspec(align(64)) struct {
                 AVX512c8f64 z;   // 2 * L1D lines
		 AVX512c8f64 zd;  // 2 * L1D lines
		 AVX512c8f64 zf;  // same as above
		 AVX512c8f64 zfd;
		 AVX512c8f64 zgd;
		 AVX512c8f64 zp;
		 AVX512c8f64 zq;
		 AVX512c8f64 zw;
		 int32_t i,it,j,nr;
#if (USE_STRUCT_PADDING) == 1
                 PAD_TO(0,48)
#endif
	     }dh; // data hot
#elif defined __linux
             __attribute__((align(64))) struct {
                 AVX512c8f64 z;   // 2 * L1D lines
		 AVX512c8f64 zd;  // 2 * L1D lines
		 AVX512c8f64 zf;  // same as above
		 AVX512c8f64 zfd;
		 AVX512c8f64 zgd;
		 AVX512c8f64 zp;
		 AVX512c8f64 zq;
		 AVX512c8f64 zw;
		 int32_t i,it,j,nr;
#if (USE_STRUCT_PADDING) == 1
                 PAD_TO(0,48)
#endif
	     }dh;
#endif
             dh.z   = AVX512c8f64{};
	     dh.zd  = AVX512c8f64{};
	     dh.zf  = AVX512c8f64{};
	     dh.zfd = AVX512c8f64{};
	     dh.zgd = AVX512c8f64{};
	     dh.zp  = AVX512c8f64{};
	     dh.zq  = AVX512c8f64{};
	     dh.zw  = AVX512c8f64{};
	     
             dc.pu  = v8_n0;
	     dc.pv  = v8_n0;
	     dc.px  = v8_n0;
	     dc.py  = v8_n0;
	     dc.w   = v8_n0;
	     dc.w0  = v8_n0;
	     dh.nr  = 0;
	     for (dh.nr = 1; dh.nr != nt; ++dh.nr) {
                  const Vec8 vnr = Vec8_CVTI4(Vec8_SETI4(dh.nr));
		  dc.pu = Vec8_SQRT(Vec8_MUL(v8PI,Vec8_SUB(Vec8_MUL(v8_n4,vnr),v8_1over2)));
		  dc.pv = Vec8_MUL(v8PI,Vec8_SQRT(Vec8_SUB(Vec8_MUL(v8_n2,vnr),v8_1over4)));
		  const Vec8 t0 = Vec8_MUL(v8_1over2,dc.pu);
		  const Vec8 t1 = Vec8_DIV(Vec8_LOG(dc.pv),dc.pu);
		  dc.px = Vec8_SUB(t0,Vec8_MUL(v8_1over2,t1));
		  dc.py = Vec8_ADD(t0,Vec8_MUL(v8_1over2,t1));
		  dh.z = AVX512c8f64{dc.px,dc.py};
		  dh.it = 0;

		  do {

		        dh.it += 1;
			v8_cerf_pd(dh.z,dh.zf,dh.zd);
			dh.zp = AVX512c8f64{v8_n1,v8_n0};
			for (dh.i = 1; dh.i != dh.nr-1; ++dh.i) {
                              dh.zp = dh.zp * (dh.z - zo[i]);
			}
                        dh.zfd = dh.zf / sh.zp;
			dh.zq = AVX512c8f64{};
			for (dh.i = 1; dh.i != dh.nr-1; ++dh.i) {
                               dh.zw = AVX512c8f64{v8_n1,v8_n0};
			       for (dh.j = 1; dh.j != dh.nr-1; ++dh.j) {
                                    if( dh.j != dh.i) {
                                        dh.zw = dh.zw * (dh.z - zo[j]);
				    }
			       }
			       dh.zq = dh.zq + dh.zw;
			}
			dh.zgd = (dh.zd - dh.zq * dh.zfd) / dh.zp;
			dh.z = dh.z - dh.zfd / dh.zgd;
			dc.w0 = dc.w;
			dc.w = Vec8_ABS(dh.z.m_re);
			const Vec8 vcmp1 = Vec8_ABS(Vec8_DIV(Vec8_SUB(dc.w,dc.w0),dc.w));
			
		  }while(50 < dh.it || Vec8_CMP(vcmp1,Vec8_SET1(0.00000000001),_CMP_LE_OQ) );

		   zo[dh.nr] = dh.z;
	    }  
       }

       void v8_cfc_pd(const AVX512c8f64 z,
                      AVX512c8f64 &zf,
		      AVX512c8f64 &zd)      {
#if defined _WIN64
	      __declspec(align(64)) struct {
                  AVX512c8f64 c;
		  AVX512c8f64 cf;
		  AVX512c8f64 cf0;
		  AVX512c8f64 cf1;
		  AVX512c8f64 cg;
		  AVX512c8f64 cr;
		  AVX512c8f64 z0;
		  AVX512c8f64 zp;
		  AVX512c8f64 zp2;
		  int32_t k,m;
#if (USE_STRUCT_PADDING) == 1
                  PAD_TO(0,58)
#endif
	      }dh;
#elif defined __linux
             __attribute__((align(64))) struct {
                  AVX512c8f64 c;
		  AVX512c8f64 cf;
		  AVX512c8f64 cf0;
		  AVX512c8f64 cf1;
		  AVX512c8f64 cg;
		  AVX512c8f64 cr;
		  AVX512c8f64 z0;
		  AVX512c8f64 zp;
		  AVX512c8f64 zp2;
		  int32_t k,m;
#if (USE_STRUCT_PADDING) == 1
                  PAD_TO(0,58)
#endif
	     }dh;
#endif
#if defined _WIN64
             __declspec(align(64)) struct {
                  Vec8 w0;
		  Vec8 wa;
		  Vec8 wa0;
	     }dc;
#elif defined __linux
             __attribute__((align(64))) struct {
                  Vec8 w0;
		  Vec8 wa;
		  Vec8 wa0;
	     }dc;
#endif
                  // First touch
		  dh.c   = AVX512c8f64{};
		  dh.cf  = AVX512c8f64{};
		  dh.cf0 = AVX512c8f64{};
		  dh.cf1 = AVX512c8f64{};
		  dh.cg  = AVX512c8f64{};
		  dh.cr  = AVX512c8f64{};
		  dh.z0  = AVX512c8f64{};
		  dh.zp  = AVX512c8f64{};
		  dh.zp2 = AVX512c8f64{};
		  dc.w0  = v8_n0;
		  dc.wa  = v8_n0;
		  dc.wa0 = v8_n0;

		  dc.w0  = cabs(z);
		  dh.zp  = Vec8_MUL(v8_1over2,v8PI) * dh.z * dh.z;
		  dh.zp2 = dh.z * dh.z;
		  const std::pair<__mmask8> cmp = z == dh.z0;
		  if(cmp.first && cmp.second) {
                     dh.c = dh.z0;
		}
		  else if (Vec8_CMP(dc.w0,Vec8_SET1(2.5),_CMP_LE_OQ)) {
                     dh.k = 0;
		     dh.cr = z;
		     dh.c  = dh.cr;
		     for (dh.k = 1; dh.k != 80; ++dh.k) {
                           const Vec8 vk =  Vec8_CVTI4(Vec8_SETI4(dh.k));
			   const Vec8 t0 =  Vec8_FMSUB(v8_n4,vk,v8_n3);
			   const Vec8 t1 =  Vec8_FMSUB(v8_n2,vk,v8_n1);
			   const Vec8 t2 =  Vec8_FMAD(v8_n4,vk,v8_n1);
			   dh.cr = flip_sign(v8_1over2) * dh.cr * t0 /
			           vk / t1 / t2 * dh.zp2;
			   dh.c += dh.cr;
			   dc.wa = cabs(dh.c);
			   const Vec8 vcmp = Vec8_ABS((Vec8_DIV(Vec8_SUB(dc.wa,dc.wa0),dc.wa)));
			   if (Vec8_CMP(vcmp,Vec8_SET1(0.00000000000001),_CMP_LT_OQ)) {
                                  return;
			   }
			   dc.wa0 = dc.wa;
		     }
		}
		 else if (Vec8_CMP(Vec8_SET1(2.5),dc.w0,_CMP_LT_OQ) &&
		          Vec8_CMP(dc.w0,Vec8_SET1(4.5),_CMP_LT_OQ))  {

                          dh.m = 85;
			  dh.c = dh.z0;
			  dh.cf1 = dh.z0;
			  dh.cf0  = AVX512c8f64{Vec8_SET1(1.0E-30),v8_n0};
			  for (dh.k = dh.m; dh.k != 0; --dh.k) {
                                const Vec8 vk =  Vec8_CVTI4(Vec8_SETI4(dh.k));
				const Vec8 t0 = Vec8_FMAD(v8_n2,vk.v8_n3);
				dh.cf = t0 * dh.cf0 / dh.zp - dh.cf1;
				if(dh.k == (dh.k/2)*2) {
                                    dh.c += dh.cf;             
				}
				dh.cf1 = dh.cf0;
				dh.cf0 = dh.cf;
			  }
			  const AVX512c8f64 tmp1 = csqrt(v8_n2/(v8PI*dh.zp));
			  const AVX512c8f64 tmp2 = csin(dh.zp)/dh.cf*dh.c;
			  dh.c = tmp1 * tmp2;
	       }
	        else {
                          dh.cr = AVX512c8f64{v8_n1,v8_n0}
			  dh.cf = AVX512c8f64{v8_n1,v8_n0};
			  for (dh.k = 1; dh.k != 20; ++dh.k) {
                                const Vec8 vk =  Vec8_CVTI4(Vec8_SETI4(dh.k));
				const Vec8 t0 = Vec8_FMSUB(v8_n4,vk,v8_n1);
				const Vec8 t1 = Vec8_FMSUB(v8_n4,vk,v8_n3);
				dh.cr = flip_sign(v8_1over2) * dh.cr * t0 * t1 / dh.zp2;
				dh.cf += dh.cr;
			  }
			  dh.cr = v8_n1 / (v8PI * z * z);
			  dh.cg = dh.cr;
			  for (dh.k = 1; dh.k != 12; ++dh.k) {
			         const Vec8 vk =  Vec8_CVTI4(Vec8_SETI4(dh.k));
                                 const Vec8 t0 = Vec8_FMAD(v8_n4,vk.v8_n1);
				 const Vec8 t1 = Vec8_FMSUB(v4_n8,vk,v8_n1);
				 dh.cr = flip_sign(v8_1over2) * dh.cr * t0 * t1 / dh.zp2;
				 dh.cg += dh.cr;
			  }
			  const AVX512c8f64 tmp3 = dh.cf * csin(dh.zp) - dh.cg * ccos(dh.zp);
			  dh.c = v8_1over2 + tmp3 / (v8PI * z);
		}
		zf = dh.c;
		zd = ccos(v8_1over2*v8PI*z*z);
       }



	
		
	 }  // math
} // gms


	
#endif /*__GMS_VECSPECFUNCS_H__*/
