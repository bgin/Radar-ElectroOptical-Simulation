
#ifndef __GMS_PDF_CDF_SSE_HPP__
#define __GMS_PDF_CDF_SSE_HPP__ 060120240832

namespace file_info {

 const unsigned int gGMS_PDF_CDF_SSE_MAJOR = 1U;
 const unsigned int gGMS_PDF_CDF_SSE_MINOR = 0U;
 const unsigned int gGMS_PDF_CDF_SSE_MICRO = 0U;
 const unsigned int gGMS_PDF_CDF_SSE_FULLVER =
  1000U*gGMS_PDF_CDF_SSE_MAJOR+100U*gGMS_PDF_CDF_SSE_MINOR+10U*gGMS_PDF_CDF_SSE_MICRO;
 const char * const pgGMS_PDF_CDF_SSE_CREATION_DATE = "06-01-2024 08:34 +00200 (SAT 06 JAN 2024 08:34 GMT+2)";
 const char * const pgGMS_PDF_CDF_SSE_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_PDF_CDF_SSE_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_PDF_CDF_SSE_SYNOPSIS      = "Manually vectorized [SSE] PDF,CDF functions"


}

#include <immintrin.h>
#include <limits>
#include "GMS_config.h"
#include "GMS_simd_utils.hpp"


namespace gms {

        namespace math {


 /*
      !*****************************************************************************80
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ).
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from Hart et al, while approximations
!    for X < 12.0D+00 are similar to those in Cody and Hillstrom,
!    but are unpublished.
!
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine dependent
!    constants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma
!    function of X.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) BETA, the radix for the floating-point
!    representation.
!
!    Local, integer MAXEXP, the smallest positive power of BETA that overflows.
!
!    Local, real ( kind = 8 ) XBIG, the largest argument for which
!    LN(GAMMA(X)) is representable in the machine, the solution to the equation
!      LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    Local, real ( kind = 8 ) FRTBIG, a rough estimate of the fourth root
!    of XBIG.
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG     FRTBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62D+2461  3.13D+615
!  Cyber 180/855 (S.P.)        2        1070       1.72D+319   6.44D+79
!  IEEE (IBM/XT) (S.P.)        2         128       4.08D+36    1.42D+9
!  IEEE (IBM/XT) (D.P.)        2        1024       2.55D+305   2.25D+76
!  IBM 3033      (D.P.)       16          63       4.29D+73    2.56D+18
!  VAX D-Format  (D.P.)        2         127       2.05D+36    1.20D+9
!  VAX G-Format  (D.P.)        2        1023       1.28D+305   1.89D+76
!    
*/

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d gamma_log_xmm2r8(const __m128d x) {
		      
                           // if(__builtin_expect(_mm_cmp_pd_mask(x,_0,_CMP_LE_OQ),0) ||
			   //    __builtin_expect(_mm_cmp_pd_mask(x,xbig,_CMP_GT_OQ),0)) {
                           //    return (huge);
			   // }
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(16) static __m128d c[7] = { _mm_set1_pd(-1.910444077728E-03),
			                                     _mm_set1_pd(8.4171387781295E-04),
                                                             _mm_set1_pd(-5.952379913043012E-04), 
                                                             _mm_set1_pd(7.93650793500350248E-04), 
                                                             _mm_set1_pd(-2.777777777777681622553E-03), 
                                                             _mm_set1_pd(8.333333333333333331554247E-02), 
                                                             _mm_set1_pd(5.7083835261E-03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128d p1[8] = {_mm_set1_pd(4.945235359296727046734888E+00), 
                                                             _mm_set1_pd(2.018112620856775083915565E+02), 
                                                             _mm_set1_pd(2.290838373831346393026739E+03), 
                                                             _mm_set1_pd(1.131967205903380828685045E+04),
                                                             _mm_set1_pd(2.855724635671635335736389E+04), 
                                                             _mm_set1_pd(3.848496228443793359990269E+04), 
                                                             _mm_set1_pd(2.637748787624195437963534E+04), 
                                                             _mm_set1_pd(7.225813979700288197698961E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128d p2[8] = {_mm_set1_pd(4.974607845568932035012064E+00), 
                                                             _mm_set1_pd(5.424138599891070494101986E+02), 
                                                             _mm_set1_pd(1.550693864978364947665077E+04), 
                                                             _mm_set1_pd(1.847932904445632425417223E+05), 
                                                             _mm_set1_pd(1.088204769468828767498470E+06), 
                                                             _mm_set1_pd(3.338152967987029735917223E+06), 
                                                             _mm_set1_pd(5.106661678927352456275255E+06), 
                                                             _mm_set1_pd(3.074109054850539556250927E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128d p4[8] = {_mm_set1_pd(1.474502166059939948905062E+04), 
                                                             _mm_set1_pd(2.426813369486704502836312E+06), 
                                                             _mm_set1_pd(1.214755574045093227939592E+08), 
                                                             _mm_set1_pd(2.663432449630976949898078E+09), 
                                                             _mm_set1_pd(2.940378956634553899906876E+10), 
                                                             _mm_set1_pd(1.702665737765398868392998E+11), 
                                                             _mm_set1_pd(4.926125793377430887588120E+11), 
                                                             _mm_set1_pd(5.606251856223951465078242E+11)};
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(16) static __m128d q1[8] = {_mm_set1_pd(6.748212550303777196073036E+01), 
                                                             _mm_set1_pd(1.113332393857199323513008E+03), 
                                                             _mm_set1_pd(7.738757056935398733233834E+03), 
                                                             _mm_set1_pd(2.763987074403340708898585E+04), 
                                                             _mm_set1_pd(5.499310206226157329794414E+04), 
                                                             _mm_set1_pd(6.161122180066002127833352E+04), 
                                                             _mm_set1_pd(3.635127591501940507276287E+04), 
                                                             _mm_set1_pd(8.785536302431013170870835E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128d q2[8] = {_mm_set1_pd(1.830328399370592604055942E+02),
                                                             _mm_set1_pd(7.765049321445005871323047E+03), 
                                                             _mm_set1_pd(1.331903827966074194402448E+05),
                                                             _mm_set1_pd(1.136705821321969608938755E+06), 
                                                             _mm_set1_pd(5.267964117437946917577538E+06), 
                                                             _mm_set1_pd(1.346701454311101692290052E+07), 
                                                             _mm_set1_pd(1.782736530353274213975932E+07), 
                                                             _mm_set1_pd(9.533095591844353613395747E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128d q4[8] = {_mm_set1_pd(2.690530175870899333379843E+03), 
                                                             _mm_set1_pd(6.393885654300092398984238E+05), 
                                                             _mm_set1_pd(4.135599930241388052042842E+07), 
                                                             _mm_set1_pd(1.120872109616147941376570E+09), 
                                                             _mm_set1_pd(1.488613728678813811542398E+10), 
                                                             _mm_set1_pd(1.016803586272438228077304E+11), 
                                                             _mm_set1_pd(3.417476345507377132798597E+11), 
                                                             _mm_set1_pd(4.463158187419713286462081E+11)};
			    const __m128d d1     = _mm_set1_pd(-5.772156649015328605195174E-01);
			    const __m128d d2     = _mm_set1_pd(4.227843350984671393993777E-01);
                            const __m128d d4     = _mm_set1_pd(1.791759469228055000094023E+00);
                            const __m128d frtbig = _mm_set1_pd(1.42E+09);
                            const __m128d pnt68  = _mm_set1_pd(0.6796875E+00);
			    const __m128d sqrtpi = _mm_set1_pd(0.9189385332046727417803297E+00);
			    const __m128d xbig   = _mm_set1_pd(4.08E+36);
			    const __m128d _0     = _mm_setzero_pd();
			    const __m128d _1_2   = _mm_set1_pd(0.5);
			    const __m128d _1_5   = _mm_set1_pd(1.5);
			    const __m128d _1     = _mm_set1_pd(1.0);
			    const __m128d _4     = _mm_set1_pd(4.0);
			    const __m128d _2     = _mm_set1_pd(2.0);
			    const __m128d _12    = _mm_set1_pd(12.0);
			    const __m128d huge   = _mm_set1_pd(std::numeric_limits<double>::max());
			    const __m128d eps    = _mm_set1_pd(std::numeric_limits<double>::epsilon());
			    __m128d gamlog,res,xden;
			    __m128d xm1,xm2,xm4;
			    __m128d xnum,xsq,corr;
			    gamlog = _mm_setzero_pd();
			   
			    if(_mm_cmp_pd_mask(x,eps,_CMP_LE_OQ)) {
                               res = negate_xmm2r8(_mm_log_pd(x));
			    }
			    else if(_mm_cmp_pd_mask(x,_1_5,_CMP_LE_OQ)) {
                               const __mmask8 m0 = _mm_cmp_pd_mask(x,pnt68,_CMP_LT_OQ);
			       corr = _mm_mask_blend_pd(m0,_0,xmm2r8_negate(_mm_log_pd(x)));
			       xm1  = _mm_mask_blend_pd(m0,_mm_sub_pd(
			                                                _mm_sub_pd(x,_1_2),_1_2));

			       if(_mm_cmp_pd_mask(x,_1_2,_CMP_LE_OQ) ||
			          _mm_cmp_pd_mask(pnt68,x,_CMP_LE_OQ)) {
                                   xden = _1;
				   xnum = _0;
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[0]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[0]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[1]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[1]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[2]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[2]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[3]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[3]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[4]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[4]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[5]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[5]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[6]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[6]);
				   xnum = _mm_fmadd_pd(xnum,xm1,p1[7]);
				   xden = _mm_fmadd_pd(xden,xm1,q1[7]);
				   const __m128d t0 = _mm_fmadd_pd(xm1,
				                                  _mm_div_pd(xnum,xden),d1);
				   res  = _mm_add_pd(corr,
				                    _mm_mul_pd(xm1,t0));
				}
				else {

                                   xm2  = _mm_sub_pd(_mm_sub_pd(x,_1_2),_1_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[0]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[0]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[1]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[1]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[2]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[2]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[3]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[3]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[4]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[4]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[5]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[5]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[6]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[6]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[7]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[7]);
				   const __m128d t0 = _mm_fmadd_pd(xm2,
				                                  _mm_div_pd(xnum,xden),d2);
				   res  = _mm_add_pd(corr,
				                    _mm_mul_pd(xm2,t0));
				}
			    }
			    else if(_mm_cmp_pd_mask(x,_4,_CMP_LE_OQ)) {
                                   xm2  = _mm_sub_pd(x,_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[0]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[0]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[1]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[1]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[2]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[2]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[3]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[3]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[4]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[4]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[5]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[5]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[6]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[6]);
				   xnum = _mm_fmadd_pd(xnum,xm2,p2[7]);
				   xden = _mm_fmadd_pd(xden,xm2,q2[7]);
				   res  = _mm_mul_pd(xm2,
				                    _mm_fmadd_pd(xm2,
						                _mm_div_pd(xnum,xden),d2));
			    }
			    else if(_mm_cmp_pd_mask(x,_12,_CMP_LE_OQ)) {
                                   xm4  = _mm_sub_pd(x,_4);
				   xden = negate_xmm2r8(_1);
				   xnum = _0;
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[0]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[0]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[1]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[1]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[2]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[2]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[3]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[3]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[4]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[4]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[5]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[5]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[6]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[6]);
				   xnum = _mm_fmadd_pd(xnum,xm4,p4[7]);
				   xden = _mm_fmadd_pd(xden,xm4,q4[7]);
				   res  = _mm_fmadd_pd(xm4,_mm_div_pd(xnum,xden),d4);
			    }
			    else {
                                   res  = _0;
				   if(_mm_cmp_pd_mask(x,frtbig,_CMP_LE_OQ)) {
                                      res = c[6];
				      xsq = _mm_mul_pd(x,x);
				      res = _mm_add_pd(_mm_div_pd(res,xsq),c[0]);
				      res = _mm_add_pd(_mm_div_pd(res,xsq),c[1]);
				      res = _mm_add_pd(_mm_div_pd(res,xsq),c[2]);
				      res = _mm_add_pd(_mm_div_pd(res,xsq),c[3]);
				      res = _mm_add_pd(_mm_div_pd(res,xsq),c[4]);
				      res = _mm_add_pd(_mm_div_pd(res,xsq),c[5]);
				   }
                                   res  = _mm_div_pd(res,x);
				   corr = _mm_log_pd(x);
				   res  = _mm_sub_pd(_mm_add_pd(res,sqrtpi),
				                        _mm_mul_pd(_1_2,corr));
				   res  = _mm_fmadd_pd(x,_mm_sub_pd(corr,_1),res);
				   
			    }

			    gamlog = res;
			    return (gamlog);
			    
		  }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 gamma_log_xmm4r4(const __m128 x) {
		      
                          //  if(__builtin_expect(_mm_cmp_ps_mask(x,_0,_CMP_LE_OQ),0) ||
			 //      __builtin_expect(_mm_cmp_ps_mask(x,xbig,_CMP_GT_OQ),0)) {
                         //      return (huge);
			 //   }
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(16) static __m128 c[7] = { _mm_set1_ps(-1.910444077728E-03),
			                                     _mm_set1_ps(8.4171387781295E-04),
                                                             _mm_set1_ps(-5.952379913043012E-04), 
                                                             _mm_set1_ps(7.93650793500350248E-04), 
                                                             _mm_set1_ps(-2.777777777777681622553E-03), 
                                                             _mm_set1_ps(8.333333333333333331554247E-02), 
                                                             _mm_set1_ps(5.7083835261E-03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128 p1[8] = {_mm_set1_ps(4.945235359296727046734888E+00), 
                                                             _mm_set1_ps(2.018112620856775083915565E+02), 
                                                             _mm_set1_ps(2.290838373831346393026739E+03), 
                                                             _mm_set1_ps(1.131967205903380828685045E+04),
                                                             _mm_set1_ps(2.855724635671635335736389E+04), 
                                                             _mm_set1_ps(3.848496228443793359990269E+04), 
                                                             _mm_set1_ps(2.637748787624195437963534E+04), 
                                                             _mm_set1_ps(7.225813979700288197698961E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128 p2[8] = {_mm_set1_ps(4.974607845568932035012064E+00), 
                                                             _mm_set1_ps(5.424138599891070494101986E+02), 
                                                             _mm_set1_ps(1.550693864978364947665077E+04), 
                                                             _mm_set1_ps(1.847932904445632425417223E+05), 
                                                             _mm_set1_ps(1.088204769468828767498470E+06), 
                                                             _mm_set1_ps(3.338152967987029735917223E+06), 
                                                             _mm_set1_ps(5.106661678927352456275255E+06), 
                                                             _mm_set1_ps(3.074109054850539556250927E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128 p4[8] = {_mm_set1_ps(1.474502166059939948905062E+04), 
                                                             _mm_set1_ps(2.426813369486704502836312E+06), 
                                                             _mm_set1_ps(1.214755574045093227939592E+08), 
                                                             _mm_set1_ps(2.663432449630976949898078E+09), 
                                                             _mm_set1_ps(2.940378956634553899906876E+10), 
                                                             _mm_set1_ps(1.702665737765398868392998E+11), 
                                                             _mm_set1_ps(4.926125793377430887588120E+11), 
                                                             _mm_set1_ps(5.606251856223951465078242E+11)};
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(16) static __m128 q1[8] = {_mm_set1_ps(6.748212550303777196073036E+01), 
                                                             _mm_set1_ps(1.113332393857199323513008E+03), 
                                                             _mm_set1_ps(7.738757056935398733233834E+03), 
                                                             _mm_set1_ps(2.763987074403340708898585E+04), 
                                                             _mm_set1_ps(5.499310206226157329794414E+04), 
                                                             _mm_set1_ps(6.161122180066002127833352E+04), 
                                                             _mm_set1_ps(3.635127591501940507276287E+04), 
                                                             _mm_set1_ps(8.785536302431013170870835E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128 q2[8] = {_mm_set1_ps(1.830328399370592604055942E+02),
                                                             _mm_set1_ps(7.765049321445005871323047E+03), 
                                                             _mm_set1_ps(1.331903827966074194402448E+05),
                                                             _mm_set1_ps(1.136705821321969608938755E+06), 
                                                             _mm_set1_ps(5.267964117437946917577538E+06), 
                                                             _mm_set1_ps(1.346701454311101692290052E+07), 
                                                             _mm_set1_ps(1.782736530353274213975932E+07), 
                                                             _mm_set1_ps(9.533095591844353613395747E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(16) static __m128 q4[8] = {_mm_set1_ps(2.690530175870899333379843E+03), 
                                                             _mm_set1_ps(6.393885654300092398984238E+05), 
                                                             _mm_set1_ps(4.135599930241388052042842E+07), 
                                                             _mm_set1_ps(1.120872109616147941376570E+09), 
                                                             _mm_set1_ps(1.488613728678813811542398E+10), 
                                                             _mm_set1_ps(1.016803586272438228077304E+11), 
                                                             _mm_set1_ps(3.417476345507377132798597E+11), 
                                                             _mm_set1_ps(4.463158187419713286462081E+11)};
			    const __m128 d1     = _mm_set1_ps(-5.772156649015328605195174E-01);
			    const __m128 d2     = _mm_set1_ps(4.227843350984671393993777E-01);
                            const __m128 d4     = _mm_set1_ps(1.791759469228055000094023E+00);
                            const __m128 frtbig = _mm_set1_ps(1.42E+09);
                            const __m128 pnt68  = _mm_set1_ps(0.6796875E+00);
			    const __m128 sqrtpi = _mm_set1_ps(0.9189385332046727417803297E+00);
			    const __m128 xbig   = _mm_set1_ps(4.08E+36);
			    const __m128 _0     = _mm_setzero_ps();
			    const __m128 _1_2   = _mm_set1_ps(0.5);
			    const __m128 _1_5   = _mm_set1_ps(1.5);
			    const __m128 _1     = _mm_set1_ps(1.0);
			    const __m128 _4     = _mm_set1_ps(4.0);
			    const __m128 _2     = _mm_set1_ps(2.0);
			    const __m128 _12    = _mm_set1_ps(12.0);
			    const __m128 huge   = _mm_set1_ps(std::numeric_limits<float>::max());
			    const __m128 eps    = _mm_set1_ps(std::numeric_limits<float>::epsilon());
			    __m128 gamlog,res,xden;
			    __m128 xm1,xm2,xm4;
			    __m128 xnum,xsq,corr;
			    gamlog = _mm_setzero_ps();
			   
			    if(_mm_cmp_ps_mask(x,eps,_CMP_LE_OQ)) {
                               res = negate_xmm4r4(_mm_log_ps(x));
			    }
			    else if(_mm_cmp_ps_mask(x,_1_5,_CMP_LE_OQ)) {
                               const __mmask8 m0 = _mm_cmp_ps_mask(x,pnt68,_CMP_LT_OQ);
			       corr = _mm_mask_blend_ps(m0,_0,negate_xmm4r4(_mm_log_ps(x)));
			       xm1  = _mm_mask_blend_ps(m0,_mm_sub_ps(
			                                                _mm_sub_ps(x,_1_2),_1_2));

			       if(_mm_cmp_ps_mask(x,_1_2,_CMP_LE_OQ) ||
			          _mm_cmp_ps_mask(pnt68,x,_CMP_LE_OQ)) {
                                   xden = _1;
				   xnum = _0;
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[0]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[0]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[1]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[1]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[2]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[2]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[3]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[3]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[4]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[4]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[5]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[5]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[6]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[6]);
				   xnum = _mm_fmadd_ps(xnum,xm1,p1[7]);
				   xden = _mm_fmadd_ps(xden,xm1,q1[7]);
				   const __m128 t0 = _mm_fmadd_ps(xm1,
				                                  _mm_div_ps(xnum,xden),d1);
				   res  = _mm_add_ps(corr,
				                    _mm_mul_ps(xm1,t0));
				}
				else {

                                   xm2  = _mm_sub_ps(_mm_sub_ps(x,_1_2),_1_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[0]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[0]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[1]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[1]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[2]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[2]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[3]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[3]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[4]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[4]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[5]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[5]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[6]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[6]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[7]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[7]);
				   const __m128 t0 = _mm_fmadd_ps(xm2,
				                                  _mm_div_ps(xnum,xden),d2);
				   res  = _mm_add_ps(corr,
				                    _mm_mul_ps(xm2,t0));
				}
			    }
			    else if(_mm_cmp_ps_mask(x,_4,_CMP_LE_OQ)) {
                                   xm2  = _mm_sub_ps(x,_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[0]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[0]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[1]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[1]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[2]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[2]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[3]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[3]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[4]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[4]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[5]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[5]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[6]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[6]);
				   xnum = _mm_fmadd_ps(xnum,xm2,p2[7]);
				   xden = _mm_fmadd_ps(xden,xm2,q2[7]);
				   res  = _mm_mul_ps(xm2,
				                    _mm_fmadd_ps(xm2,
						                _mm_div_ps(xnum,xden),d2));
			    }
			    else if(_mm_cmp_ps_mask(x,_12,_CMP_LE_OQ)) {
                                   xm4  = _mm_sub_ps(x,_4);
				   xden = xmm4r4_negate(_1);
				   xnum = _0;
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[0]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[0]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[1]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[1]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[2]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[2]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[3]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[3]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[4]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[4]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[5]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[5]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[6]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[6]);
				   xnum = _mm_fmadd_ps(xnum,xm4,p4[7]);
				   xden = _mm_fmadd_ps(xden,xm4,q4[7]);
				   res  = _mm_fmadd_ps(xm4,_mm_div_ps(xnum,xden),d4);
			    }
			    else {
                                   res  = _0;
				   if(_mm_cmp_ps_mask(x,frtbig,_CMP_LE_OQ)) {
                                      res = c[6];
				      xsq = _mm_mul_ps(x,x);
				      res = _mm_add_ps(_mm_div_ps(res,xsq),c[0]);
				      res = _mm_add_ps(_mm_div_ps(res,xsq),c[1]);
				      res = _mm_add_ps(_mm_div_ps(res,xsq),c[2]);
				      res = _mm_add_ps(_mm_div_ps(res,xsq),c[3]);
				      res = _mm_add_ps(_mm_div_ps(res,xsq),c[4]);
				      res = _mm_add_ps(_mm_div_ps(res,xsq),c[5]);
				   }
                                   res  = _mm_div_ps(res,x);
				   corr = _mm_log_ps(x);
				   res  = _mm_sub_ps(_mm_add_ps(res,sqrtpi),
				                        _mm_mul_ps(_1_2,corr));
				   res  = _mm_fmadd_ps(x,_mm_sub_ps(corr,_1),res);
				   
			    }

			    gamlog = res;
			    return (gamlog);
			    
		  }


/*
!*****************************************************************************80
!
!! GAMMA_INC computes the incomplete Gamma function.
!
!  Discussion:
!
!    GAMMA_INC(P,       0) = 0,
!    GAMMA_INC(P,Infinity) = 1.
!
!    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    Original FORTRAN77 version by B L Shea.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BL Shea,
!    Chi-squared and Incomplete Gamma Integral,
!    Algorithm AS239,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the exponent parameter.
!    0.0D+00 < P.
!
!    Input, real ( kind = 8 ) X, the integral limit parameter.
!    If X is less than or equal to 0, GAMMA_INC is returned as 0.
!
!    Output, real ( kind = 8 ) GAMMA_INC, the value of the function.
!
*/	

#include <limits>

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d  
                      gamma_incomplete_xmm2r8(const __m128d p,
                                              const __m128d x) {
                                              
                            const __m128d exp_arg_min = _mm_set1_pd(-88.0e+00);
                            const __m128d overflow    = _mm_set1_pd(1.0e+37);
                            const __m128d plimit      = _mm_set1_pd(1000.0e+00);
                            const __m128d tol         = _mm_set1_pd(1.0e-7);
                            const __m128d xbig        = _mm_set1_pd(1.0e+8);
                            const __m128d C0          = _mm_setzero_pd();
                            const __m128d C0333333333 = _mm_set1_pd(0.3333333333333333333333);
                            const __m128d C1          = _mm_set1_pd(1.0);
                            const __m128d C2          = _mm_set1_pd(2.0);
                            const __m128d C3          = _mm_set1_pd(3.0);
                            const __m128d C9          = _mm_set1_pd(9.0);
                            __m128d cdf,arg,b,c;
                            __m128d pn1,pn2,pn3,pn4;
                            __m128d pn5,pn6,rn,t0,t1;
                            __m128d gaminc;
                            __mmask8 m0,m1;
                            m0 = _mm_cmp_pd_mask(plimit,p,_CMP_LT_OQ);
                            if(m0) {
                               __m128d sqrp,xp,_9p1,t0,t1;
                               xp     = _mm_div_pd(x,p);
                               _9p1   = _mm_fmsub_pd(C9,p,C1);
                               sqrp   = _mm_mul_pd(C3,_mm_sqrt_pd(p));
                               t0     = _mm_pow_pd(xp,C0333333333);
                               t1     = _mm_add_pd(t0,
                                            _mm_div_pd(C1,_9p1));
                               pn1    = _mm_mul_pd(sqrp,t1);
                               gaminc = normal_01_cdf_xmm2r8(pn1);
                               return (gaminc);
                            }   
                            m0 = _mm_cmp_pd_mask(x,C1,_CMP_LE_OQ);
                            m1 = _mm_cmp_pd_mask(x,p,_CMP_LT_OQ);
                            if(m0 || m1) {

                               t0  = _mm_log_pd(x);
                          
                               t1  = gamma_log_xmm2r8(_mm_add_pd(p,C1));
                               arg = _mm_fmsub_pd(p,t0,_mm_sub_pd(x,t1));
                               c   = C1;
                               gaminc = C1;
                               a   = p; 
                               while(true) {
                                    a      = _mm_add_pd(a,C1);
                                    c      = _mm_mul_pd(c,_mm_div_pd(x,a));
                                    gaminc = _mm_add_pd(gaminc,c);
                                    m0     = _mm_cmp_pd_mask(c,tol,_CMP_LE_OQ);
                                    if(m0) break;
                               }

                               t0  = _mm_log_pd(x);
                               arg = _mm_add_pd(arg,t0);  
                               m1  = _mm_cmp_pd_mask(exp_arg_min,arg,_CMP_LE_OQ);
                               gaminc = _mm_mask_blend_pd(m1,C0,_mm_exp_pd(arg));  
                                
                           } 
                           else {

                               t0  = _mm_log_pd(x);
                             
                               t1  = gamma_log_xmm2r8(p);
                               arg = _mm_fmsub_pd(p,t0,_mm_sub_pd(x,t1));                               
                               a   = _mm_sub_pd(C1,p);
                               b   = _mm_add_pd(a,_mm_add_pd(x,C1));
                               c   = C0;
                               pn1 = C1;
                               pn2 = x;
                               pn3 = _mm_add_pd(x,C1);
                               pn4 = _mm_mul_pd(x,b);
                               gaminc = _mm_div_pd(pn3,pn4);
                               while(true) {
                                   a = _mm_add_pd(a,C1);
                                   b = _mm_add_pd(b,C2);
                                   c = _mm_add_pd(c,C1);
                                   pn5 = _mm_fmsub_pd(b,pn3,
                                                     _mm_mul_pd(a,
                                                           _mm_mul_pd(c,pn1)));
                                   pn6 = _mm_fmsub_pd(b,pn4,
                                                     _mm_mul_pd(a,
                                                           _mm_mul_pd(c,pn2)));
                                   if(_mm_cmp_pd_mask(C0,_mm_abs_pd(pn6),
                                                                       _CMP_LT_OQ)) {
                                        rn = _mm_div_pd(pn5,pn6);
                                        t0 = _mm_abs_pd(_mm_sub_pd(gaminc,rn));
                                        t1 = _mm_min_pd(tol,_mm_mul_pd(tol,rn));
                                        if(_mm_cmp_pd_mask(t0,t1,_CMP_LE_OQ)) {
                                           arg  = _mm_add_pd(_mm_log_pd(gaminc));       
                                           m1   = _mm_cmp_pd_mask(exp_arg_min,arg,_CMP_LE_OQ);
                                           gaminc = _mm_mask_blend_pd(m1,C1,_mm_sub_pd(C1,
                                                                                    _mm_exp_pd(arg)));
      
                                           return (gaminc);                               
                                        }    
                                        gaminc = rn;                               
                                   }
                                   pn1 = pn3;
                                   pn2 = pn4;
                                   pn3 = pn5;
                                   pn4 = pn6;
                                   if(_mm_cmp_pd_mask(overflow,
                                                   _mm_abs_pd(pn5),_CMP_LE_OQ)) {
                                      t0 = _mm_div_pd(C1,overflow);
                                      pn1= _mm_mul_pd(pn1,t0);
                                      pn2= _mm_mul_pd(pn2,t0);
                                      pn3= _mm_mul_pd(pn3,t0);
                                      pn4= _mm_mul_pd(pn4,t0);               
                                   }
                               }
                           } 
                           
                           return (gaminc);                   
                   }	
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128  
                      gamma_incomplete_xmm4r4(const __m128 p,
                                              const __m128 x) {
                                              
                            const __m128 exp_arg_min = _mm_set1_ps(-88.0e+00);
                            const __m128 overflow    = _mm_set1_ps(std::numeric_limits<float>::max());
                            const __m128 plimit      = _mm_set1_ps(1000.0e+00);
                            const __m128 tol         = _mm_set1_ps(1.0e-7f);
                            const __m128 xbig        = _mm_set1_ps(1.0e+8f);
                            const __m128 C0          = _mm_setzero_ps();
                            const __m128 C0333333333 = _mm_set1_ps(0.3333333333333333333333f);
                            const __m128 C1          = _mm_set1_ps(1.0);
                            const __m128 C2          = _mm_set1_ps(2.0);
                            const __m128 C3          = _mm_set1_ps(3.0);
                            const __m128 C9          = _mm_set1_ps(9.0);
                            __m128 cdf,arg,b,c;
                            __m128 pn1,pn2,pn3,pn4;
                            __m128 pn5,pn6,rn,t0,t1;
                            __m128 gaminc;
                            __mmask8 m0,m1;
                            m0 = _mm_cmp_ps_mask(plimit,p,_CMP_LT_OQ);
                            if(m0) {
                               __m128 sqrp,xp,_9p1,t0,t1;
                               xp     = _mm_div_ps(x,p);
                               _9p1   = _mm_fmsub_ps(C9,p,C1);
                               sqrp   = _mm_mul_ps(C3,_mm_sqrt_pd(p));
                               t0     = _mm_pow_ps(xp,C0333333333);
                               t1     = _mm_add_ps(t0,
                                            _mm_div_ps(C1,_9p1));
                               pn1    = _mm_mul_ps(sqrp,t1);
                               gaminc = normal_01_cdf_xmm4r4(pn1);
                               return (gaminc);
                            }   
                            m0 = _mm_cmp_ps_mask(x,C1,_CMP_LE_OQ);
                            m1 = _mm_cmp_ps_mask(x,p,_CMP_LT_OQ);
                            if(m0 || m1) {

                               t0  = _mm_log_ps(x);
                          
                               t1  = gamma_log_xmm4r4(_mm_add_ps(p,C1));
                               arg = _mm_fmsub_ps(p,t0,_mm_sub_ps(x,t1));
                               c   = C1;
                               gaminc = C1;
                               a   = p; 
                               while(true) {
                                    a      = _mm_add_ps(a,C1);
                                    c      = _mm_mul_ps(c,_mm_div_ps(x,a));
                                    gaminc = _mm_add_ps(gaminc,c);
                                    m0     = _mm_cmp_ps_mask(c,tol,_CMP_LE_OQ);
                                    if(m0) break;
                               }

                               t0  = _mm_log_ps(x);
                               arg = _mm_add_ps(arg,t0);  
                               m1  = _mm_cmp_ps_mask(exp_arg_min,arg,_CMP_LE_OQ);
                               gaminc = _mm_mask_blend_ps(m1,C0,_mm_exp_pd(arg));  
                                
                           } 
                           else {

                               t0  = _mm_log_ps(x);
                             
                               t1  = gamma_log_xmm4r4(p);
                               arg = _mm_fmsub_ps(p,t0,_mm_sub_pd(x,t1));                               
                               a   = _mm_sub_ps(C1,p);
                               b   = _mm_add_ps(a,_mm_add_pd(x,C1));
                               c   = C0;
                               pn1 = C1;
                               pn2 = x;
                               pn3 = _mm_add_ps(x,C1);
                               pn4 = _mm_mul_ps(x,b);
                               gaminc = _mm_div_ps(pn3,pn4);
                               while(true) {
                                   a = _mm_add_ps(a,C1);
                                   b = _mm_add_ps(b,C2);
                                   c = _mm_add_ps(c,C1);
                                   pn5 = _mm_fmsub_ps(b,pn3,
                                                     _mm_mul_ps(a,
                                                           _mm_mul_ps(c,pn1)));
                                   pn6 = _mm_fmsub_ps(b,pn4,
                                                     _mm_mul_ps(a,
                                                           _mm_mul_ps(c,pn2)));
                                   if(_mm_cmp_ps_mask(C0,_mm_abs_ps(pn6),
                                                                       _CMP_LT_OQ)) {
                                        rn = _mm_div_ps(pn5,pn6);
                                        t0 = _mm_abs_ps(_mm_sub_ps(gaminc,rn));
                                        t1 = _mm_min_ps(tol,_mm_mul_ps(tol,rn));
                                        if(_mm_cmp_ps_mask(t0,t1,_CMP_LE_OQ)) {
                                           arg  = _mm_add_ps(_mm_log_pd(gaminc));       
                                           m1   = _mm_cmp_ps_mask(exp_arg_min,arg,_CMP_LE_OQ);
                                           gaminc = _mm_mask_blend_ps(m1,C1,_mm_sub_ps(C1,
                                                                                    _mm_exp_ps(arg)));
      
                                           return (gaminc);                               
                                        }    
                                        gaminc = rn;                               
                                   }
                                   pn1 = pn3;
                                   pn2 = pn4;
                                   pn3 = pn5;
                                   pn4 = pn6;
                                   if(_mm_cmp_ps_mask(overflow,
                                                   _mm_abs_ps(pn5),_CMP_LE_OQ)) {
                                      t0 = _mm_div_ps(C1,overflow);
                                      pn1= _mm_mul_ps(pn1,t0);
                                      pn2= _mm_mul_ps(pn2,t0);
                                      pn3= _mm_mul_ps(pn3,t0);
                                      pn4= _mm_mul_ps(pn4,t0);               
                                   }
                               }
                           } 
                           
                           return (gaminc);                   
                   }	
                   
                   

       /*
    !*****************************************************************************80
!
!! R8POLY_VALUE evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE, the value of the polynomial at X.
!     
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d  
		      vpoly_eval_xmm2r8(const int32_t n,
		                        const __m128d * __restrict __ATTR_ALIGN__(16) a,
		                        const __m128d x) {
		         
		         register __m128d vpoly;
		         vpoly = _mm_load_pd(&a[n]);
		         for(int32_t i=n; i != 0; --i) {
		             register __m128d t0 = a[i];
		             vpoly = _mm_fmadd_pd(vpoly,x,t0);   
		         }  
		         return (vpoly);              
		    }
		    
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128  
		      vpoly_eval_xmm4r4(const int32_t n,
		                        const __m128 * __restrict __ATTR_ALIGN__(16) a,
		                        const __m128 x) {
		         
		         register __m128 vpoly;
		         vpoly = _mm_load_ps(&a[n]);
		         for(int32_t i=n; i != 0; --i) {
		             register __m128 t0 = a[i];
		             vpoly = _mm_fmadd_ps(vpoly,x,t0);   
		         }  
		         return (vpoly);              
		    }	 
		    
		   
		   
/*

       !*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.    

*/ 


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d    
		      normal_01_cdf_inv_xmm2r8(const __m128d p) {
		            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(16) static __m128d  a[8] = {
		                     _mm_set1_pd(3.3871328727963666080e+00),
                                     _mm_set1_pd(1.3314166789178437745e+02),
                                     _mm_set1_pd(1.9715909503065514427e+03),
                                     _mm_set1_pd(1.3731693765509461125e+04),
                                     _mm_set1_pd(4.5921953931549871457e+04),
                                     _mm_set1_pd(6.7265770927008700853e+04),
                                     _mm_set1_pd(3.3430575583588128105e+04),
                                     _mm_set1_pd(2.5090809287301226727e+03)};   
                            __attribute__((section(".rodata")))  
		            __ATTR_ALIGN__(16) static __m128d   b[8] = {
		                      _mm_set1_pd(1.0e+00),
                                      _mm_set1_pd(4.2313330701600911252e+01),
                                      _mm_set1_pd(6.8718700749205790830e+02),
                                      _mm_set1_pd(5.3941960214247511077e+03),
                                      _mm_set1_pd(2.1213794301586595867e+04),
                                      _mm_set1_pd(3.9307895800092710610e+04),
                                      _mm_set1_pd(2.8729085735721942674e+04),
                                      _mm_set1_pd(5.2264952788528545610e+03)}; 
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(16) static __m128d   c[8] = {
		                      _mm_set1_pd(1.42343711074968357734e+00),
                                      _mm_set1_pd(4.63033784615654529590e+00),
                                      _mm_set1_pd(5.76949722146069140550e+00),
                                      _mm_set1_pd(3.64784832476320460504e+00),
                                      _mm_set1_pd(1.27045825245236838258e+00),
                                      _mm_set1_pd(2.41780725177450611770e-01),
                                      _mm_set1_pd(2.27238449892691845833e-02),
                                      _mm_set1_pd(7.74545014278341407640e-04)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128d   d[8] = {
                                      _mm_set1_pd(1.0e+00),
                                      _mm_set1_pd(2.05319162663775882187e+00),
                                      _mm_set1_pd(1.67638483018380384940e+00),
                                      _mm_set1_pd(6.89767334985100004550e-01),
                                      _mm_set1_pd(1.48103976427480074590e-01),
                                      _mm_set1_pd(1.51986665636164571966e-02),
                                      _mm_set1_pd(5.47593808499534494600e-04),
                                      _mm_set1_pd(1.05075007164441684324e-09)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128d   e[8] = {
                                      _mm_set1_pd(6.65790464350110377720e+00),
                                      _mm_set1_pd(5.46378491116411436990e+00),
                                      _mm_set1_pd(1.78482653991729133580e+00),
                                      _mm_set1_pd(2.96560571828504891230e-01),
                                      _mm_set1_pd(2.65321895265761230930e-02),
                                      _mm_set1_pd(1.24266094738807843860e-03),
                                      _mm_set1_pd(2.71155556874348757815e-05),
                                      _mm_set1_pd(2.01033439929228813265e-07)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128d   f[8] = {
                                      _mm_set1_pd(1.0e+00),
                                      _mm_set1_pd(5.99832206555887937690e-01),
                                      _mm_set1_pd(1.36929880922735805310e-01),
                                      _mm_set1_pd(1.48753612908506148525e-02),
                                      _mm_set1_pd(7.86869131145613259100e-04), 
                                      _mm_set1_pd(1.84631831751005468180e-05),
                                      _mm_set1_pd(1.42151175831644588870e-07),
                                      _mm_set1_pd(2.04426310338993978564e-15)};
                          __m128d const1 = _mm_set1_pd(0.180625e+00);
                          __m128d const2 = _mm_set1_pd(1.6e+00);
                          __m128d split1 = _mm_set1_pd(0.425e+00);
                          __m128d split2 = _mm_set1_pd(5.0e+00);
                          __m128d C0     = _mm_setzero_pd();
                          __m128d C1     = _mm_set1_pd(1.0);
                          __m128d C05    = _mm_set1_pd(0.5);
                          register __m128d q,r,t0,t1;
                          register __m128d x;
                          q = _mm_sub_pd(p,C05);
                          if(_mm_cmp_pd_mask(q,split1,_CMP_LE_OQ)) {
                             r = _mm_sub_pd(const1,_mm_mul_pd(q,q));
                             t0= vpoly_eval_xmm2r8(8,a,r);
                             t1= vpoly_eval_xmm2r8(8,b,r);
                             x = _mm_div_pd(_mm_mul_pd(q,t0),t1);
                          } 
                          else {
                             const __mmask8 m = _mm_cmp_pd_mask(q,C0,_CMP_LT_OQ);
                             r                = _mm_mask_blend_pd(m,_mm_sub_pd(C1,p),p);
                             if(_mm_cmp_pd_mask(r,C0,_CMP_LE_OQ)) {
                                x = _mm_set1_pd(std::numeric_limits<double>::max());
                             }
                             else {

                                r = _mm_sqrt_pd(negate_xmm2r8(_mm_log_pd(r)));
                        
                                const __mmask8 m = _mm_cmp_pd_mask(r,split2,_CMP_LE_OQ);
                                r                = _mm_mask_blend_pd(m,_mm_sub_pd(r,split2),
                                                                          _mm_sub_pd(r,const2));
                                t0               = _mm_div_pd(vpoly_eval_xmm2r8(8,c,r),
                                                                 vpoly_eval_xmm2r8(8,d,r));
                                t1               = _mm_div_pd(vpoly_eval_xmm2r8(8,e,r),
                                                                 vpoly_eval_xmm2r8(8,f,r));
                                x                = _mm_mask_blend_pd(m,t1,t0);      
                             }
                             if(_mm_cmp_pd_mask(q,C0,_CMP_LT_OQ)) x = negate_xmm2r8(x);
                          }
                          return (x);
                          
		    }
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128    
		      normal_01_cdf_inv_xmm4r4(const __m128 p) {
		            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(16) static __m128  a[8] = {
		                     _mm_set1_ps(3.3871328727963666080e+00f),
                                     _mm_set1_ps(1.3314166789178437745e+02f),
                                     _mm_set1_ps(1.9715909503065514427e+03f),
                                     _mm_set1_ps(1.3731693765509461125e+04f),
                                     _mm_set1_ps(4.5921953931549871457e+04f),
                                     _mm_set1_ps(6.7265770927008700853e+04f),
                                     _mm_set1_ps(3.3430575583588128105e+04f),
                                     _mm_set1_ps(2.5090809287301226727e+03f)};   
                            __attribute__((section(".rodata")))  
		            __ATTR_ALIGN__(16) static __m128   b[8] = {
		                      _mm_set1_ps(1.0e+00),
                                      _mm_set1_ps(4.2313330701600911252e+01f),
                                      _mm_set1_ps(6.8718700749205790830e+02f),
                                      _mm_set1_ps(5.3941960214247511077e+03f),
                                      _mm_set1_ps(2.1213794301586595867e+04f),
                                      _mm_set1_ps(3.9307895800092710610e+04f),
                                      _mm_set1_ps(2.8729085735721942674e+04f),
                                      _mm_set1_ps(5.2264952788528545610e+03f)}; 
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(16) static __m128   c[8] = {
		                      _mm_set1_ps(1.42343711074968357734e+00f),
                                      _mm_set1_ps(4.63033784615654529590e+00f),
                                      _mm_set1_ps(5.76949722146069140550e+00f),
                                      _mm_set1_ps(3.64784832476320460504e+00f),
                                      _mm_set1_ps(1.27045825245236838258e+00f),
                                      _mm_set1_ps(2.41780725177450611770e-01f),
                                      _mm_set1_ps(2.27238449892691845833e-02f),
                                      _mm_set1_ps(7.74545014278341407640e-04f)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128   d[8] = {
                                      _mm_set1_ps(1.0e+00),
                                      _mm_set1_ps(2.05319162663775882187e+00f),
                                      _mm_set1_ps(1.67638483018380384940e+00f),
                                      _mm_set1_ps(6.89767334985100004550e-01f),
                                      _mm_set1_ps(1.48103976427480074590e-01f),
                                      _mm_set1_ps(1.51986665636164571966e-02f),
                                      _mm_set1_ps(5.47593808499534494600e-04f),
                                      _mm_set1_ps(1.05075007164441684324e-09f)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128   e[8] = {
                                      _mm_set1_ps(6.65790464350110377720e+00f),
                                      _mm_set1_ps(5.46378491116411436990e+00f),
                                      _mm_set1_ps(1.78482653991729133580e+00f),
                                      _mm_set1_ps(2.96560571828504891230e-01f),
                                      _mm_set1_ps(2.65321895265761230930e-02f),
                                      _mm_set1_ps(1.24266094738807843860e-03f),
                                      _mm_set1_ps(2.71155556874348757815e-05f),
                                      _mm_set1_ps(2.01033439929228813265e-07f)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128   f[8] = {
                                      _mm_set1_ps(1.0e+00),
                                      _mm_set1_ps(5.99832206555887937690e-01f),
                                      _mm_set1_ps(1.36929880922735805310e-01f),
                                      _mm_set1_ps(1.48753612908506148525e-02f),
                                      _mm_set1_ps(7.86869131145613259100e-04f), 
                                      _mm_set1_ps(1.84631831751005468180e-05f),
                                      _mm_set1_ps(1.42151175831644588870e-07f),
                                      _mm_set1_ps(2.04426310338993978564e-15f)};
                          __m128 const1 = _mm_set1_ps(0.180625e+00f);
                          __m128 const2 = _mm_set1_ps(1.6e+00f);
                          __m128 split1 = _mm_set1_ps(0.425e+00f);
                          __m128 split2 = _mm_set1_ps(5.0e+00f);
                          __m128 C0     = _mm_setzero_ps();
                          __m128 C1     = _mm_set1_ps(1.0f);
                          __m128 C05    = _mm_set1_ps(0.5f);
                          register __m128 q,r,t0,t1;
                          register __m128 x;
                          q = _mm_sub_ps(p,C05);
                          if(_mm_cmp_ps_mask(q,split1,_CMP_LE_OQ)) {
                             r = _mm_sub_ps(const1,_mm_mul_ps(q,q));
                             t0= vpoly_eval_xmm4r4(8,a,r);
                             t1= vpoly_eval_xmm4r4(8,b,r);
                             x = _mm_div_ps(_mm_mul_ps(q,t0),t1);
                          } 
                          else {
                             const __mmask8 m = _mm_cmp_ps_mask(q,C0,_CMP_LT_OQ);
                             r                = _mm_mask_blend_ps(m,_mm_sub_ps(C1,p),p);
                             if(_mm_cmp_ps_mask(r,C0,_CMP_LE_OQ)) {
                                x = _mm_set1_pd(std::numeric_limits<float>::max());
                             }
                             else {

                                r = _mm_sqrt_ps(negate_xmm4r4(_mm_log_ps(r)));
                       
                                const __mmask8 m = _mm_cmp_ps_mask(r,split2,_CMP_LE_OQ);
                                r                = _mm_mask_blend_ps(m,_mm_sub_ps(r,split2),
                                                                          _mm_sub_ps(r,const2));
                                t0               = _mm_div_ps(vpoly_eval_xmm4r4(8,c,r),
                                                                 vpoly_eval_xmm4r4(8,d,r));
                                t1               = _mm_div_ps(vpoly_eval_xmm4r4(8,e,r),
                                                                 vpoly_eval_xmm4r4(8,f,r));
                                x                = _mm_mask_blend_ps(m,t1,t0);      
                             }
                             if(_mm_cmp_ps_mask(q,C0,_CMP_LT_OQ)) x = negate_xmm4r4(x);
                          }
                          return (x);
                          
		    }
		    
		    
/*
       !*****************************************************************************80
!
!! RECIPROCAL_CDF evaluates the Reciprocal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/


               
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d 
		      reciprocal_cdf_xmm2r8(const __m128d x,
		                            const __m128d a,
		                            const __m128d b) {
		          
		          register __m128d ax,ab,l1,l2;
		          register __m128d cdf;       
		          ax = _mm_div_pd(a,x);
                          l1 = _mm_log_pd(ax);
	                  ab = _mm_div_pd(a,b);
                          l2 = _mm_log_pd(ab);
                          cdf= _mm_div_pd(l1,l2);
                          return (cdf);
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 
		      reciprocal_cdf_xmm4r4(const __m128 x,
		                             const __m128 a,
		                             const __m128 b) {
		          
		          register __m128 ax,ab,l1,l2;
		          register __m128 cdf;       
		          ax = _mm_div_ps(a,x);
                          l1 = _mm_log_ps(ax);
	                  ab = _mm_div_ps(a,b);
                          l2 = _mm_log_ps(ab);
                          cdf= _mm_div_ps(l1,l2);
                          return (cdf);
		     }
		     
		     
/*
          !*****************************************************************************80
!
!! RECIPROCAL_CDF_INV inverts the Reciprocal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!                      
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d 		     
		      reciprocal_cdf_inv_xmm2r8(const __m128d cdf,
		                                const __m128d a,
		                                const __m128d b) {
		         
		           register __m128d C1 = _mm_set1_pd(1.0);
		           register __m128d pow1,pow2,cdf1;
		           register __m128d inv;
		           cdf1 = _mm_sub_pd(cdf,C1);
		           pow2 = _mm_pow_pd(b,cdf);
		           pow1 = _mm_pow_pd(a,cdf1);
		           inv  = _mm_div_pd(pow2,pow1);
		           return (inv);                          
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128		     
		      reciprocal_cdf_inv_xmm4r4(const __m128 cdf,
		                                const __m128 a,
		                                const __m128 b) {
		         
		           register __m128 C1 = _mm_set1_ps(1.0f);
		           register __m128 pow1,pow2,cdf1;
		           register __m128 inv;
		           cdf1 = _mm_sub_ps(cdf,C1);
		           pow2 = _mm_pow_ps(b,cdf);
		           pow1 = _mm_pow_ps(a,cdf1);
		           inv  = _mm_div_ps(pow2,pow1);
		           return (inv);                          
		     }
		     
		     
		     
/*
     !*****************************************************************************80
!
!! RECIPROCAL_MEAN returns the mean of the Reciprocal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.
!    
*/
		         
	  
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d 
		      reciprocal_mean_xmm2r8(const __m128d a,
		                             const __m128d b) {
		           
		           register __m128d ab,amb,l1;
		           register __m128d mean;
		           amb = _mm_sub_pd(a,b);
		           ab  = _mm_div_pd(a,b);
                           l1  = _mm_log_pd(ab);
	                   mean= _mm_div_pd(amb,l1);
                           return (mean);
		     }	 
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 
		      reciprocal_mean_xmm4r4(const __m128 a,
		                             const __m128 b) {
		           
		           register __m128 ab,amb,l1;
		           register __m128 mean;
		           amb = _mm_sub_ps(a,b);
		           ab  = _mm_div_ps(a,b);
                           l1  = _mm_log_ps(ab);
		           mean= _mm_div_ps(amb,l1);
                           return (mean);
		     }	 
		     
		     
/*
           !*****************************************************************************80
!
!! RECIPROCAL_PDF evaluates the Reciprocal PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = 1.0D+00 / ( X * LOG ( B / A ) )
!    for 0.0D+00 <= X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A <= B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!                
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d
		      reciprocal_pdf_xmm2r8(    const __m128d x,
		                                const __m128d a,
		                                const __m128d b) {
		          
		          register __m128d C1 = _mm_set1_pd(1.0);
		          register __m128d ba,l1;
		          register __m128d pdf;
		          ba = _mm_div_pd(b,a);
                          l1 = _mm_mul_pd(x,_mm_log_pd(ba));
	                  pdf= _mm_div_pd(C1,l1);
                          return (pdf);                            
		    }
		    
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
		      reciprocal_pdf_xmm4r4(    const __m128 x,
		                                const __m128 a,
		                                const __m128 b) {
		          
		          register __m128 C1 = _mm_set1_ps(1.0f);
		          register __m128 ba,l1;
		          register __m128 pdf;
		          ba = _mm_div_ps(b,a);
                          l1 = _mm_mul_ps(x,_mm_log_ps(ba));
	                  pdf= _mm_div_ps(C1,l1);
                          return (pdf);                            
		    }
		    
		     


/*
 !*****************************************************************************80
!
!! BESSEL_I0 evaluates the modified Bessel function I0(X).
!
!  Discussion:
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.  This transportable program is patterned after
!    the machine dependent FUNPACK packet NATSI0, but cannot match
!    that version for efficiency or accuracy.  This version uses
!    rational functions that theoretically approximate I-SUB-0(X)
!    to at least 18 significant decimal digits.
!
!  Machine dependent constants:
!
!    beta   = Radix for the floating-point system
!    maxexp = Smallest power of beta that overflows
!    XMAX =   Largest argument acceptable to BESI0;  Solution to
!             equation:
!               W(X) * (1+1/(8*X)+9/(128*X^2) = beta^maxexp
!             where  W(X) = EXP(X)/sqrt(2*PI*X)
!
!    Approximate values for some important machines are:
!
!                             beta       maxexp       XMAX
!
!    CRAY-1        (S.P.)       2         8191       5682.810
!    Cyber 180/855
!      under NOS   (S.P.)       2         1070        745.893
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       2          128         91.900
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       2         1024        713.986
!    IBM 3033      (D.P.)      16           63        178.182
!    VAX           (S.P.)       2          127         91.203
!    VAX D-Format  (D.P.)       2          127         91.203
!    VAX G-Format  (D.P.)       2         1023        713.293
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.
!
!    Output, real ( kind = 8 ) BESSEL_I0, the value of the modified
!    Bessel function of the first kind.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d  bessesl_i0_xmm2r8(const __m128d arg) {
                            __attribute__((section(".rodata")))
                            __ATTR_ALIGN__(16) static __m128d p[15] = {_mm_set1_pd(-5.2487866627945699800E-18),
                                                                      _mm_set1_pd(-1.5982226675653184646E-14), 
                                                                      _mm_set1_pd(-2.6843448573468483278E-11), 
                                                                      _mm_set1_pd(-3.0517226450451067446E-08), 
                                                                      _mm_set1_pd(-2.5172644670688975051E-05), 
                                                                      _mm_set1_pd(-1.5453977791786851041E-02), 
                                                                      _mm_set1_pd(-7.0935347449210549190E+00), 
                                                                      _mm_set1_pd(-2.4125195876041896775E+03), 
                                                                      _mm_set1_pd(-5.9545626019847898221E+05), 
                                                                      _mm_set1_pd(-1.0313066708737980747E+08), 
                                                                      _mm_set1_pd(-1.1912746104985237192E+10), 
                                                                      _mm_set1_pd(-8.4925101247114157499E+11), 
                                                                      _mm_set1_pd(-3.2940087627407749166E+13), 
                                                                      _mm_set1_pd(-5.5050369673018427753E+14), 
                                                                      _mm_set1_pd(-2.2335582639474375249E+15)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(16) static __m128d pp[8] = {_mm_set1_pd(-3.9843750000000000000E-01), 
                                                                      _mm_set1_pd(2.9205384596336793945E+00), 
                                                                      _mm_set1_pd(-2.4708469169133954315E+00), 
                                                                      _mm_set1_pd(4.7914889422856814203E-01), 
                                                                      _mm_set1_pd(-3.7384991926068969150E-03), 
                                                                      _mm_set1_pd(-2.6801520353328635310E-03), 
                                                                      _mm_set1_pd(9.9168777670983678974E-05), 
                                                                      _mm_set1_pd(-2.1877128189032726730E-06)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(16) static __m128d q[5]  = {_mm_set1_pd(-3.7277560179962773046E+03), 
                                                                      _mm_set1_pd(6.5158506418655165707E+06), 
                                                                      _mm_set1_pd(-6.5626560740833869295E+09), 
                                                                      _mm_set1_pd(3.7604188704092954661E+12), 
                                                                      _mm_set1_pd(-9.7087946179594019126E+14)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(16) static __m128d qq[7] = {_mm_set1_pd(-3.1446690275135491500E+01), 
                                                                      _mm_set1_pd(8.5539563258012929600E+01), 
                                                                      _mm_set1_pd(-6.0228002066743340583E+01), 
                                                                      _mm_set1_pd(1.3982595353892851542E+01), 
                                                                      _mm_set1_pd(-1.1151759188741312645E+00), 
                                                                      _mm_set1_pd(3.2547697594819615062E-02), 
                                                                      _mm_set1_pd(-5.5194330231005480228E-04)};
			    const __m128d rec15                    =  _mm_set1_pd(6.6666666666666666666E-02);
			    const __m128d xmax                     =  _mm_set1_pd(91.9E+00);
			    const __m128d exp40                    =  _mm_set1_pd(2.353852668370199854E+17);
			    const __m128d _1                       =  _mm_set1_pd(1.0);
			    const __m128d _15                      =  _mm_set1_pd(15.0);
			    const __m128d _225                     =  _mm_set1_pd(225.0);
			    const __m128d _40                      =  _mm_set1_pd(40.0);
			    const __m128d eps                      =  _mm_set1_pd(std::numeric_limits<double>::epsilon());
			    const __m128d huge                     =  _mm_set1_pd(std::mumeric_limits<double>::max());
			    __m128d value,a,b,bessel_i0;
			    __m128d sump,sumq,x,xx;
                            x = _mm_abs_pd(arg);
			    if(_mm_cmp_pd_mask(x,eps,_CMP_LT_OQ)) {
                               value = _1;
			    }
			    else if(_mm_cmp_pd_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm_mul_pd(x,x);
			       sump = p[0];
			       sump = _mm_fmadd_pd(sump,xx,p[1]);
			       sump = _mm_fmadd_pd(sump,xx,p[2]);
			       sump = _mm_fmadd_pd(sump,xx,p[3]);
			       sump = _mm_fmadd_pd(sump,xx,p[4]);
			       sump = _mm_fmadd_pd(sump,xx,p[5]);
			       sump = _mm_fmadd_pd(sump,xx,p[6]);
			       sump = _mm_fmadd_pd(sump,xx,p[7]);
			       sump = _mm_fmadd_pd(sump,xx,p[8]);
			       sump = _mm_fmadd_pd(sump,xx,p[9]);
			       sump = _mm_fmadd_pd(sump,xx,p[10]);
			       sump = _mm_fmadd_pd(sump,xx,p[11]);
			       sump = _mm_fmadd_pd(sump,xx,p[12]);
			       sump = _mm_fmadd_pd(sump,xx,p[13]);
			       sump = _mm_fmadd_pd(sump,xx,p[14]);
			       xx   = _mm_sub_pd(xx,_225);
			       const __m128d xxq0 = _mm_add_pd(xx,q[0]);
			       const __m128d xxq1 = _mm_add_pd(xx,q[1]);
			       const __m128d xxq2 = _mm_add_pd(xx,q[2]);
			       const __m128d xxq3 = _mm_add_pd(xx,q[3]);
			       const __m128d xxq4 = _mm_add_pd(xx,q[4]);
			       sumq = _mm_mul_pd(xxq0,
			                        _mm_mul_pd(xxq1,
						          _mm_mul_pd(xxq2,
							            _mm_mul_pd(xxq3,xxq4))));
			       value = _mm_div_pd(sump,sumq);
			                                         
			    }
			    else if(_mm_cmp_pd_mask(_15,x,_CMP_LE_OQ)) {
                                    if(_mm_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                       value = huge;
				    }
				    else {
                                           xx = _mm_sub_pd(_mm_div_pd(_1,x),rec15);
					   const __m128d t0 = _mm_fmadd_pd(pp[0],xx,pp[1]);
					   const __m128d c0 = _mm_fmadd_pd(_mm_add_pd(xx,qq[0]),xx,qq[1]);
					   const __m128d t1 = _mm_fmadd_pd(t0,xx,pp[2]);
					   const __m128d c1 = _mm_fmadd_pd(c0,xx,qq[2]);
					   const __m128d t2 = _mm_fmadd_pd(t1,xx,pp[3]);
					   const __m128d c2 = _mm_fmadd_pd(c1,xx,qq[3]);
					   const __m128d t3 = _mm_fmadd_pd(t2,xx,pp[4]);
					   const __m128d c3 = _mm_fmadd_pd(c2,xx,qq[4]);
					   const __m128d t4 = _mm_fmadd_pd(t3,xx,pp[5]);
					   const __m128d c4 = _mm_fmadd_pd(c3,xx,qq[5]);
					   const __m128d t5 = _mm_fmadd_pd(t4,xx,pp[6]);
					   const __m128d c5 = _mm_fmadd_pd(c4,xx,qq[6]);
					   const __m128d t6 = _mm_fmadd_pd(t5,xx,pp[7]);
					   sump             = t6;
					   sumq             = c5;
					   value            = _mm_div_pd(sump,sumq);
					   const __mmask8 m = _mm_cmp_pd_mask(x,_mm_sub_pd(xmax,_15),_CMP_LE_OQ);
					   a                = _mm_mask_blend_pd(m,_mm_exp_pd(_mm_sub_pd(x,_40)),
					                                             _mm_exp_pd(x));
					   b                = _mm_mask_blend_pd(m,exp40,_1);
					   const __m128d tmp = _mm_sub_pd(_mm_mul_pd(value,a),
					                                    _mm_mul_pd(pp[0],a));
					   value            = _mm_mul_pd(_mm_div_pd(tmp,_mm_sqrt_pd(x)),b);
				    }
			    }
			   
			    bessel_i0 = value;
			    return (bessel_i0);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128  bessesl_i0_xmm4r4(const __m128 arg) {
                            __attribute__((section(".rodata")))
                            __ATTR_ALIGN__(16) static __m128 p[15] = {_mm_set1_ps(-5.2487866627945699800E-18f),
                                                                      _mm_set1_ps(-1.5982226675653184646E-14f), 
                                                                      _mm_set1_ps(-2.6843448573468483278E-11f), 
                                                                      _mm_set1_ps(-3.0517226450451067446E-08f), 
                                                                      _mm_set1_ps(-2.5172644670688975051E-05f), 
                                                                      _mm_set1_ps(-1.5453977791786851041E-02f), 
                                                                      _mm_set1_ps(-7.0935347449210549190E+00f), 
                                                                      _mm_set1_ps(-2.4125195876041896775E+03f), 
                                                                      _mm_set1_ps(-5.9545626019847898221E+05f), 
                                                                      _mm_set1_ps(-1.0313066708737980747E+08f), 
                                                                      _mm_set1_ps(-1.1912746104985237192E+10f), 
                                                                      _mm_set1_ps(-8.4925101247114157499E+11f), 
                                                                      _mm_set1_ps(-3.2940087627407749166E+13f), 
                                                                      _mm_set1_ps(-5.5050369673018427753E+14f), 
                                                                      _mm_set1_ps(-2.2335582639474375249E+15f)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(16) static __m128 pp[8] = {_mm_set1_ps(-3.9843750000000000000E-01f), 
                                                                      _mm_set1_ps(2.9205384596336793945E+00f), 
                                                                      _mm_set1_ps(-2.4708469169133954315E+00f), 
                                                                      _mm_set1_ps(4.7914889422856814203E-01f), 
                                                                      _mm_set1_ps(-3.7384991926068969150E-03f), 
                                                                      _mm_set1_ps(-2.6801520353328635310E-03f), 
                                                                      _mm_set1_ps(9.9168777670983678974E-05f), 
                                                                      _mm_set1_ps(-2.1877128189032726730E-06f)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(16) static __m128 q[5]  = {_mm_set1_ps(-3.7277560179962773046E+03f), 
                                                                      _mm_set1_ps(6.5158506418655165707E+06f), 
                                                                      _mm_set1_ps(-6.5626560740833869295E+09f), 
                                                                      _mm_set1_ps(3.7604188704092954661E+12f), 
                                                                      _mm_set1_ps(-9.7087946179594019126E+14f)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(16) static __m128 qq[7] = {_mm_set1_ps(-3.1446690275135491500E+01f), 
                                                                      _mm_set1_ps(8.5539563258012929600E+01f), 
                                                                      _mm_set1_ps(-6.0228002066743340583E+01f), 
                                                                      _mm_set1_ps(1.3982595353892851542E+01f), 
                                                                      _mm_set1_ps(-1.1151759188741312645E+00f), 
                                                                      _mm_set1_ps(3.2547697594819615062E-02f), 
                                                                      _mm_set1_ps(-5.5194330231005480228E-04f)};
			    const __m128 rec15                    =  _mm_set1_ps(6.6666666666666666666E-02f);
			    const __m128 xmax                     =  _mm_set1_ps(91.9E+00f);
			    const __m128 exp40                    =  _mm_set1_ps(2.353852668370199854E+17f);
			    const __m128 _1                       =  _mm_set1_ps(1.0f);
			    const __m128 _15                      =  _mm_set1_ps(15.0f);
			    const __m128 _225                     =  _mm_set1_ps(225.0f);
			    const __m128 _40                      =  _mm_set1_ps(40.0f);
			    const __m128 eps                      =  _mm_set1_pd(std::numeric_limits<float>::epsilon());
			    const __m128 huge                     =  _mm_set1_pd(std::mumeric_limits<float>::max());
			    __m128 value,a,b,bessel_i0;
			    __m128 sump,sumq,x,xx;
                            x = _mm_abs_ps(arg);
			    if(_mm_cmp_ps_mask(x,eps,_CMP_LT_OQ)) {
                               value = _1;
			    }
			    else if(_mm_cmp_ps_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm_mul_ps(x,x);
			       sump = p[0];
			       sump = _mm_fmadd_ps(sump,xx,p[1]);
			       sump = _mm_fmadd_ps(sump,xx,p[2]);
			       sump = _mm_fmadd_ps(sump,xx,p[3]);
			       sump = _mm_fmadd_ps(sump,xx,p[4]);
			       sump = _mm_fmadd_ps(sump,xx,p[5]);
			       sump = _mm_fmadd_ps(sump,xx,p[6]);
			       sump = _mm_fmadd_ps(sump,xx,p[7]);
			       sump = _mm_fmadd_ps(sump,xx,p[8]);
			       sump = _mm_fmadd_ps(sump,xx,p[9]);
			       sump = _mm_fmadd_ps(sump,xx,p[10]);
			       sump = _mm_fmadd_ps(sump,xx,p[11]);
			       sump = _mm_fmadd_ps(sump,xx,p[12]);
			       sump = _mm_fmadd_ps(sump,xx,p[13]);
			       sump = _mm_fmadd_ps(sump,xx,p[14]);
			       xx   = _mm_sub_ps(xx,_225);
			       const __m128 xxq0 = _mm_add_ps(xx,q[0]);
			       const __m128 xxq1 = _mm_add_ps(xx,q[1]);
			       const __m128 xxq2 = _mm_add_ps(xx,q[2]);
			       const __m128 xxq3 = _mm_add_ps(xx,q[3]);
			       const __m128 xxq4 = _mm_add_ps(xx,q[4]);
			       sumq = _mm_mul_ps(xxq0,
			                        _mm_mul_ps(xxq1,
						          _mm_mul_ps(xxq2,
							            _mm_mul_ps(xxq3,xxq4))));
			       value = _mm_div_ps(sump,sumq);
			                                         
			    }
			    else if(_mm_cmp_ps_mask(_15,x,_CMP_LE_OQ)) {
                                    if(_mm_cmp_ps_mask(xmax,x,_CMP_LT_OQ)) {
                                       value = huge;
				    }
				    else {
                                           xx = _mm_sub_ps(_mm_div_ps(_1,x),rec15);
					   const __m128 t0 = _mm_fmadd_ps(pp[0],xx,pp[1]);
					   const __m128 c0 = _mm_fmadd_ps(_mm_add_ps(xx,qq[0]),xx,qq[1]);
					   const __m128 t1 = _mm_fmadd_ps(t0,xx,pp[2]);
					   const __m128 c1 = _mm_fmadd_ps(c0,xx,qq[2]);
					   const __m128 t2 = _mm_fmadd_ps(t1,xx,pp[3]);
					   const __m128 c2 = _mm_fmadd_ps(c1,xx,qq[3]);
					   const __m128 t3 = _mm_fmadd_ps(t2,xx,pp[4]);
					   const __m128 c3 = _mm_fmadd_ps(c2,xx,qq[4]);
					   const __m128 t4 = _mm_fmadd_ps(t3,xx,pp[5]);
					   const __m128 c4 = _mm_fmadd_ps(c3,xx,qq[5]);
					   const __m128 t5 = _mm_fmadd_ps(t4,xx,pp[6]);
					   const __m128 c5 = _mm_fmadd_ps(c4,xx,qq[6]);
					   const __m128 t6 = _mm_fmadd_ps(t5,xx,pp[7]);
					   sump             = t6;
					   sumq             = c5;
					   value            = _mm_div_ps(sump,sumq);
					   const __mmask8 m = _mm_cmp_ps_mask(x,_mm_sub_ps(xmax,_15),_CMP_LE_OQ);

					   a                = _mm_mask_blend_ps(m,_mm_exp_ps(_mm_sub_ps(x,_40)),
					                                             _mm_exp_ps(x));
   					   b                = _mm_mask_blend_ps(m,exp40,_1);
					   const __m128 tmp = _mm_sub_ps(_mm_mul_ps(value,a),
					                                    _mm_mul_ps(pp[0],a));
					   value            = _mm_mul_ps(_mm_div_ps(tmp,_mm_sqrt_ps(x)),b);
				    }
			    }
			   
			    bessel_i0 = value;
			    return (bessel_i0);
		    }



/*
 !*****************************************************************************80
!
!! BESSEL_I1 evaluates the Bessel I function of order I.
!
!  Discussion:
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards.
!    This transportable program is patterned after the machine-dependent
!    FUNPACK packet NATSI1, but cannot match that version for efficiency
!    or accuracy.  This version uses rational functions that theoretically
!    approximate I-SUB-1(X) to at least 18 significant decimal digits.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    the intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Machine-dependent constants:
!
!    beta   = Radix for the floating-point system.
!    maxexp = Smallest power of beta that overflows.
!    XMAX =   Largest argument acceptable to BESI1;  Solution to
!             equation:
!               EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp
!
!
!    Approximate values for some important machines are:
!
!                            beta       maxexp    XMAX
!
!    CRAY-1        (S.P.)       2         8191    5682.810
!    Cyber 180/855
!      under NOS   (S.P.)       2         1070     745.894
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       2          128      91.906
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       2         1024     713.987
!    IBM 3033      (D.P.)      16           63     178.185
!    VAX           (S.P.)       2          127      91.209
!    VAX D-Format  (D.P.)       2          127      91.209
!    VAX G-Format  (D.P.)       2         1023     713.293
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Blair, Edwards,
!    Chalk River Report AECL-4928,
!    Atomic Energy of Canada, Limited,
!    October, 1974.
!
!  Parameters:
!
!    Input, real (kind = 8 ) ARG, the argument.
!
!    Output, real ( kind = 8 ) BESSEL_I1, the value of the Bessel
!    I1 function.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d bessel_i1_xmm2r8(const __m128d arg) {
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128d  p[15] = {_mm_set1_pd(-1.9705291802535139930E-19), 
                                                                      _mm_set1_pd(-6.5245515583151902910E-16), 
                                                                      _mm_set1_pd(-1.1928788903603238754E-12), 
                                                                      _mm_set1_pd(-1.4831904935994647675E-09), 
                                                                      _mm_set1_pd(-1.3466829827635152875E-06), 
                                                                      _mm_set1_pd(-9.1746443287817501309E-04), 
                                                                      _mm_set1_pd(-4.7207090827310162436E-01), 
                                                                      _mm_set1_pd(-1.8225946631657315931E+02), 
                                                                      _mm_set1_pd(-5.1894091982308017540E+04), 
                                                                      _mm_set1_pd(-1.0588550724769347106E+07), 
                                                                      _mm_set1_pd(-1.4828267606612366099E+09), 
                                                                      _mm_set1_pd(-1.3357437682275493024E+11), 
                                                                      _mm_set1_pd(-6.9876779648010090070E+12), 
                                                                      _mm_set1_pd(-1.7732037840791591320E+14), 
                                                                      _mm_set1_pd(-1.4577180278143463643E+15)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(16) static __m128d pp[8]  = {_mm_set1_pd(-6.0437159056137600000E-02), 
                                                                      _mm_set1_pd(4.5748122901933459000E-01), 
                                                                      _mm_set1_pd(-4.2843766903304806403E-01), 
                                                                      _mm_set1_pd(9.7356000150886612134E-02), 
                                                                      _mm_set1_pd(-3.2457723974465568321E-03), 
                                                                      _mm_set1_pd(-3.6395264712121795296E-04), 
                                                                      _mm_set1_pd(1.6258661867440836395E-05), 
                                                                      _mm_set1_pd(-3.6347578404608223492E-07)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(16) static __m128d q[5]   = {_mm_set1_pd(-4.0076864679904189921E+03), 
                                                                      _mm_set1_pd(7.4810580356655069138E+06), 
                                                                      _mm_set1_pd(-8.0059518998619764991E+09), 
                                                                      _mm_set1_pd(4.8544714258273622913E+12), 
                                                                      _mm_set1_pd(-1.3218168307321442305E+15)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(16) static __m128d qq[6]  = {_mm_set1_pd(-3.8806586721556593450E+00), 
                                                                      _mm_set1_pd(3.2593714889036996297E+00), 
                                                                      _mm_set1_pd(-8.5017476463217924408E-01), 
                                                                      _mm_set1_pd(7.4212010813186530069E-02), 
                                                                      _mm_set1_pd(-2.2835624489492512649E-03), 
                                                                      _mm_set1_pd(3.7510433111922824643E-05)};
			   const __m128d exp40                     =  _mm_set1_pd(2.353852668370199854E+17);
			   const __m128d _40                       =  _mm_set1_pd(40.0);
			   const __m128d _1_2                      =  _mm_set1_pd(0.5);
			   const __m128d _1                        =  _mm_set1_pd(1.0);
			   const __m128d _15                       =  _mm_set1_pd(15.0);
			   const __m128d pbar                      =  _mm_set1_pd(3.98437500E-01);
			   const __m128d rec15                     =  _mm_set1_pd(6.6666666666666666666E-02);
			   const __m128d _225                      =  _mm_set1_pd(225.0);
			   const __m128d xmax                      =  _mm_set1_pd(713.987E+00);
			   const __m128d _0                        =  _mm_setzero_pd();
			   const __m128d eps                       =  _mm_set1_pd(std::numeric_limits<double>::epsilon());
			   const __m128d huge                      =  _mm_set1_pd(std::mumeric_limits<double>::max());
			   __m128d a,b,bessel_i1,value;
			   __m128d sump,sumq,x,xx;

			   x  = _mm_abs_pd(arg);
			   if(_mm_cmp_pd_mask(x,eps,_CMP_LT_OQ)) {
                               value = _mm_mul_pd(_1_2,x);
			   }
			   else if(_mm_cmp_pd_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm_mul_pd(x,x);
			       sump = p[0];
			       sump = _mm_fmadd_pd(sump,xx,p[1]);
			       sump = _mm_fmadd_pd(sump,xx,p[2]);
			       sump = _mm_fmadd_pd(sump,xx,p[3]);
			       sump = _mm_fmadd_pd(sump,xx,p[4]);
			       sump = _mm_fmadd_pd(sump,xx,p[5]);
			       sump = _mm_fmadd_pd(sump,xx,p[6]);
			       sump = _mm_fmadd_pd(sump,xx,p[7]);
			       sump = _mm_fmadd_pd(sump,xx,p[8]);
			       sump = _mm_fmadd_pd(sump,xx,p[9]);
			       sump = _mm_fmadd_pd(sump,xx,p[10]);
			       sump = _mm_fmadd_pd(sump,xx,p[11]);
			       sump = _mm_fmadd_pd(sump,xx,p[12]);
			       sump = _mm_fmadd_pd(sump,xx,p[13]);
			       sump = _mm_fmadd_pd(sump,xx,p[14]);
			       xx   = _mm_sub_pd(xx,_225);
			       const __m128d t0 = _mm_fmadd_pd(_mm_add_pd(xx,q[0]),xx,q[1]);
			       const __m128d t1 = _mm_fmadd_pd(t0,xx,q[2]);
			       const __m128d t2 = _mm_fmadd_pd(t1,xx,q[3]);
			       const __m128d t3 = _mm_fmadd_pd(t2,xx,q[4]);
			       sumq             = t3;
			       value            = _mm_mul_pd(_mm_div_pd(sump,sumq),x);
			   }
			   else if(_mm_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                               value            = huge;
			   }
			   else {
                               xx               = _mm_sub_pd(_mm_div_pd(_1,x),rec15);
			       const __m128d t0 = _mm_fmadd_pd(pp[0],xx,pp[1]);
			       const __m128d c0 = _mm_fmadd_pd(_mm_add_pd(xx,qq[0]),xx,qq[1]);
			       const __m128d t1 = _mm_fmadd_pd(t0,xx,pp[2]);
			       const __m128d c1 = _mm_fmadd_pd(c0,xx,qq[2]);
			       const __m128d t2 = _mm_fmadd_pd(t1,xx,pp[3]);
			       const __m128d c2 = _mm_fmadd_pd(c1,xx,qq[3]);
			       const __m128d t3 = _mm_fmadd_pd(t2,xx,pp[4]);
			       const __m128d c3 = _mm_fmadd_pd(c2,xx,qq[4]);
			       const __m128d t4 = _mm_fmadd_pd(t3,xx,pp[5]);
			       const __m128d c4 = _mm_fmadd_pd(c3,xx,qq[5]);
			       const __m128d t5 = _mm_fmadd_pd(t4,xx,pp[6]);
			       const __m128d c5 = _mm_fmadd_pd(c4,xx,qq[6]);
			       const __m128d t6 = _mm_fmadd_pd(t5,xx,pp[7]);
			       sump             = t6;
			       sumq             = c5;
			       value            = _mm_div_pd(sump,sumq);
			       const __mmask8 m = _mm_cmp_pd_mask(_mm_sub_pd(xmax,_15),_CMP_LT_OQ);

			       a                = _mm_mask_blend_pd(m,_mm_exp_pd(x),
			                                                           _mm_exp_pd(_mm_sub_pd(x,_40)));
			       b                = _mm_mask_blend_pd(m,_1,_40);
			       const __m128d tmp= _mm_add_pd(_mm_mul_pd(value,a),
			                                        _mm_mul_pd(pbar,a));
			       value            = _mm_mul_pd(_mm_div_pd(tmp,_mm_sqrt_pd(x)),b);
			   }
			   if(_mm_cmp_pd_mask(arg,_0,_CMP_LT_OQ)) {
                              value             = negate_xmm2r8(value);
			   }
			   bessel_i1            = value
			   return (bessel_i1);
		    }




		     __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 bessel_i1_xmm4r4(const __m128 arg) {
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(16) static __m128  p[15] = {_mm_set1_ps(-1.9705291802535139930E-19f), 
                                                                      _mm_set1_ps(-6.5245515583151902910E-16f), 
                                                                      _mm_set1_ps(-1.1928788903603238754E-12f), 
                                                                      _mm_set1_ps(-1.4831904935994647675E-09f), 
                                                                      _mm_set1_ps(-1.3466829827635152875E-06f), 
                                                                      _mm_set1_ps(-9.1746443287817501309E-04f), 
                                                                      _mm_set1_ps(-4.7207090827310162436E-01f), 
                                                                      _mm_set1_ps(-1.8225946631657315931E+02f), 
                                                                      _mm_set1_ps(-5.1894091982308017540E+04f), 
                                                                      _mm_set1_ps(-1.0588550724769347106E+07f), 
                                                                      _mm_set1_ps(-1.4828267606612366099E+09f), 
                                                                      _mm_set1_ps(-1.3357437682275493024E+11f), 
                                                                      _mm_set1_ps(-6.9876779648010090070E+12f), 
                                                                      _mm_set1_ps(-1.7732037840791591320E+14f), 
                                                                      _mm_set1_ps(-1.4577180278143463643E+15f)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(16) static __m128 pp[8]  = {_mm_set1_ps(-6.0437159056137600000E-02f), 
                                                                      _mm_set1_ps(4.5748122901933459000E-01f), 
                                                                      _mm_set1_ps(-4.2843766903304806403E-01f), 
                                                                      _mm_set1_ps(9.7356000150886612134E-02f), 
                                                                      _mm_set1_ps(-3.2457723974465568321E-03f), 
                                                                      _mm_set1_ps(-3.6395264712121795296E-04f), 
                                                                      _mm_set1_ps(1.6258661867440836395E-05f), 
                                                                      _mm_set1_ps(-3.6347578404608223492E-07f)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(16) static __m128 q[5]   = {_mm_set1_ps(-4.0076864679904189921E+03f), 
                                                                      _mm_set1_ps(7.4810580356655069138E+06f), 
                                                                      _mm_set1_ps(-8.0059518998619764991E+09f), 
                                                                      _mm_set1_ps(4.8544714258273622913E+12f), 
                                                                      _mm_set1_ps(-1.3218168307321442305E+15f)};
                          __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(16) static __m128 qq[6]  = {_mm_set1_ps(-3.8806586721556593450E+00f), 
                                                                      _mm_set1_ps(3.2593714889036996297E+00f), 
                                                                      _mm_set1_ps(-8.5017476463217924408E-01f), 
                                                                      _mm_set1_ps(7.4212010813186530069E-02f), 
                                                                      _mm_set1_ps(-2.2835624489492512649E-03f), 
                                                                      _mm_set1_ps(3.7510433111922824643E-05f)};
			   const __m128 exp40                     =  _mm_set1_ps(2.353852668370199854E+17f);
			   const __m128 _40                       =  _mm_set1_ps(40.0f);
			   const __m128 _1_2                      =  _mm_set1_ps(0.5f);
			   const __m128 _1                        =  _mm_set1_ps(1.0f);
			   const __m128 _15                       =  _mm_set1_ps(15.0f);
			   const __m128 pbar                      =  _mm_set1_ps(3.98437500E-01f);
			   const __m128 rec15                     =  _mm_set1_ps(6.6666666666666666666E-02f);
			   const __m128 _225                      =  _mm_set1_ps(225.0f);
			   const __m128 xmax                      =  _mm_set1_ps(713.987E+00f);
			   const __m128 _0                        =  _mm_setzero_ps();
			   const __m128 eps                       =  _mm_set1_ps(std::numeric_limits<float>::epsilon());
			   const __m128 huge                      =  _mm_set1_ps(std::mumeric_limits<float>::max());
			   __m128 a,b,bessel_i1,value;
			   __m128 sump,sumq,x,xx;

			   x  = _mm_abs_ps(arg);
			   if(_mm_cmp_ps_mask(x,eps,_CMP_LT_OQ)) {
                               value = _mm_mul_ps(_1_2,x);
			   }
			   else if(_mm_cmp_ps_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm_mul_ps(x,x);
			       sump = p[0];
			       sump = _mm_fmadd_ps(sump,xx,p[1]);
			       sump = _mm_fmadd_ps(sump,xx,p[2]);
			       sump = _mm_fmadd_ps(sump,xx,p[3]);
			       sump = _mm_fmadd_ps(sump,xx,p[4]);
			       sump = _mm_fmadd_ps(sump,xx,p[5]);
			       sump = _mm_fmadd_ps(sump,xx,p[6]);
			       sump = _mm_fmadd_ps(sump,xx,p[7]);
			       sump = _mm_fmadd_ps(sump,xx,p[8]);
			       sump = _mm_fmadd_ps(sump,xx,p[9]);
			       sump = _mm_fmadd_ps(sump,xx,p[10]);
			       sump = _mm_fmadd_ps(sump,xx,p[11]);
			       sump = _mm_fmadd_ps(sump,xx,p[12]);
			       sump = _mm_fmadd_ps(sump,xx,p[13]);
			       sump = _mm_fmadd_ps(sump,xx,p[14]);
			       xx   = _mm_sub_ps(xx,_225);
			       const __m128 t0 = _mm_fmadd_ps(_mm_add_ps(xx,q[0]),xx,q[1]);
			       const __m128 t1 = _mm_fmadd_ps(t0,xx,q[2]);
			       const __m128 t2 = _mm_fmadd_ps(t1,xx,q[3]);
			       const __m128 t3 = _mm_fmadd_ps(t2,xx,q[4]);
			       sumq             = t3;
			       value            = _mm_mul_ps(_mm_div_ps(sump,sumq),x);
			   }
			   else if(_mm_cmp_ps_mask(xmax,x,_CMP_LT_OQ)) {
                               value            = huge;
			   }
			   else {
                               xx               = _mm_sub_ps(_mm_div_ps(_1,x),rec15);
			       const __m128 t0 = _mm_fmadd_ps(pp[0],xx,pp[1]);
			       const __m128 c0 = _mm_fmadd_ps(_mm_add_ps(xx,qq[0]),xx,qq[1]);
			       const __m128 t1 = _mm_fmadd_ps(t0,xx,pp[2]);
			       const __m128 c1 = _mm_fmadd_ps(c0,xx,qq[2]);
			       const __m128 t2 = _mm_fmadd_ps(t1,xx,pp[3]);
			       const __m128 c2 = _mm_fmadd_ps(c1,xx,qq[3]);
			       const __m128 t3 = _mm_fmadd_ps(t2,xx,pp[4]);
			       const __m128 c3 = _mm_fmadd_ps(c2,xx,qq[4]);
			       const __m128 t4 = _mm_fmadd_ps(t3,xx,pp[5]);
			       const __m128 c4 = _mm_fmadd_ps(c3,xx,qq[5]);
			       const __m128 t5 = _mm_fmadd_ps(t4,xx,pp[6]);
			       const __m128 c5 = _mm_fmadd_ps(c4,xx,qq[6]);
			       const __m128 t6 = _mm_fmadd_ps(t5,xx,pp[7]);
			       sump             = t6;
			       sumq             = c5;
			       value            = _mm_div_ps(sump,sumq);
			       const __mmask8 m = _mm_cmp_ps_mask(_mm_sub_ps(xmax,_15),_CMP_LT_OQ);

			       a                = _mm_mask_blend_ps(m,_mm_exp_ps(x),
			                                                           _mm_exp_ps(_mm_sub_ps(x,_40)));
			       b                = _mm_mask_blend_ps(m,_1,_40);
			       const __m128 tmp= _mm_add_ps(_mm_mul_ps(value,a),
			                                        _mm_mul_ps(pbar,a));
			       value            = _mm_mul_ps(_mm_div_ps(tmp,_mm_sqrt_ps(x)),b);
			   }
			   if(_mm_cmp_ps_mask(arg,_0,_CMP_LT_OQ)) {
                              value             = negate_xmm4r4(value);
			   }
			   bessel_i1            = value
			   return (bessel_i1);
		    }


		    

/*
!*****************************************************************************80
!
!! BETA returns the value of the Beta function.
!
!  Discussion:
!
!    The Beta function is defined as
!
!      BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
!                = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the function.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) BETA, the value of the function.
!

*/

		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d beta_xmm2r8(const __m128d a,
		                          const __m128d b) {

                       // const __m128d nan = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			const __m128d _0  = _mm_setzero_pd();
			//if(__builtin_expect(_mm_cmp_pd_mask(a,_0,_CMP_LE_OQ),0) ||
			//   __builtin_expect(_mm_cmp_pd_mask(b,_0,_CMP_LE_OQ),0)) {
                       //    return (nan);
			//}
			const __m128d ab  = _mm_add_pd(a,b);
			__m128d beta      = _mm_setzero_pd();
			beta              = _mm_exp_pd(
			                              _mm_sub_pd(
						                 _mm_add_pd(gamma_log_xmm2r8(a),
								               gamma_log_xmm2r8(b)),
			return (beta);						                     gamma_log_xmm2r8(ab)));
		    }
                 







/*!*****************************************************************************80
!
!! ANGLIT_SAMPLE samples the Anglit PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
*/

#if defined(__ICC) || defined(__INTEL_COMPILER)
#include <svrng.h>
#else
#error 'Required Intel Compiler distribution'
#endif
/*
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d anglit_sample_xmm2r8() {

                         __m128d cdf;
			 svrng_engine_t engine;
			 svrng_distribution_t uniform;
			 uint32_t seed    = 0U;
			 int32_t result   = -9999;
			 int32_t err      = -9999;
			 result           = _rdrand32_step(&seed);
			 if(!result) seed = 1563548129U;
			 engine           = svrng_new_mt19937_engine(seed);
			 err              = svrng_get_status();
			 if(err!=SVRNG_STATUS_OK) {
                            const __m128d nan = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_uniform_distribution_double(0.0,1.0);
			 const double * __restrict ptr = (const double*)(&svrng_generate8_double(engine,uniform));
			 cdf              = anglit_cdf_inv_xmm2r8(_mm_loadu_pd(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d anglit_sample_xmm2r8(const __m128 cdf) {

                            return (anglit_cdf_inv_xmm2r8(cdf));
		    }
*/

/*
      !*****************************************************************************80
!
!! ARCSIN_CDF evaluates the Arcsin CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d arcsin_cdf_xmm2r8(const __m128d x,
		                                const __m128d a) {

                         const __m128d invpi = _mm_set1_pd(0.318309886183790671537767526745);
			 const __m128d _0    = _mm_setzero_pd();
			 const __m128d _1_2  = _mm_set1_pd(0.5);
			 const __m128d _1    = _mm_set1_pd(1.0);
			 __m128d t0,cdf;
			 __mmask8 m0,m1;
			 m0  = _mm_cmp_pd_mask(x,negate_xmm2r8(a),_CMP_LE_OQ);
                         t0  = _mm_mul_pd(_mm_asin_pd(_mm_div_pd(x,a),invpi));
			 m1  = _mm_cmp_pd_mask(x,a,_CMP_LT_OQ);
                         cdf = _mm_mask_blend_pd(m0,_mm_add_pd(_1_2,t0),_0);
			 cdf = _mm_mask_blend_pd(m1,cdf,_1); 
                         return (cdf);
		   }


		   
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 arcsin_cdf_xmm4r4(const __m128 x,
		                                const __m128 a) {

                         const __m128 invpi = _mm_set1_ps(0.318309886183790671537767526745f);
			 const __m128 _0    = _mm_setzero_ps();
			 const __m128 _1_2  = _mm_set1_ps(0.5f);
			 const __m128 _1    = _mm_set1_ps(1.0f);
			 __m128 t0,cdf;
			 __mmask8 m0,m1;
			 m0  = _mm_cmp_ps_mask(x,negate_xmm4r4(a),_CMP_LE_OQ);
                         t0  = _mm_mul_ps(_mm_asin_ps(_mm_div_ps(x,a),invpi));
			 m1  = _mm_cmp_ps_mask(x,a,_CMP_LT_OQ);
                         cdf = _mm_mask_blend_ps(m0,_mm_add_pd(_1_2,t0),_0);
			 cdf = _mm_mask_blend_ps(m1,cdf,_1); 
                         return (cdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d arcsin_cdf_inv_xmm2r8(const __m128d cdf,
		                                    const __m128d a) {

                           const __m128d pi    = _mm_set1_pd(3.14159265358979323846264338328);
			   const __m128d _0    = _mm_setzero_pd();
			   const __m128d _1    = _mm_set1_pd(1.0);
			   //const __m128d nan   = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			   const __m128d _1_2  = _mm_set1_pd(0.5);
			   __m128d x;
			  // if(__builtin_expect(_mm_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			  //    __builtin_expect(_mm_cmp_pd_mask(cdf,_1,_CMP_GT_OQ),0)) {
                          //    return (nan);
			 //  }
                             x = _mm_mul_pd(_mm_sin_pd(_mm_mul_pd(pi,_mm_sub_pd(cdf,_1_2))));
                             return (x);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 arcsin_cdf_inv_xmm4r4(const __m128 cdf,
		                                    const __m128 a) {

                           const __m128 pi    = _mm_set1_ps(3.14159265358979323846264338328f);
			   const __m128 _0    = _mm_setzero_ps();
			   const __m128 _1    = _mm_set1_ps(1.0f);
			   //const __m128 nan   = _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());
			   const __m128 _1_2  = _mm_set1_ps(0.5f);
			   __m128 x;
			  // if(__builtin_expect(_mm_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			  //    __builtin_expect(_mm_cmp_ps_mask(cdf,_1,_CMP_GT_OQ),0)) {
                          //    return (nan);
			  // }
                             x = _mm_mul_ps(_mm_sin_ps(_mm_mul_ps(pi,_mm_sub_ps(cdf,_1_2))));
                             return (x);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m128d arcsin_mean_xmm2r8() {

		            return (_mm_setzero_pd());
		      }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m128 arcsin_mean_xmm4r4() {

		            return (_mm_setzero_ps());
		      }


/*
!*****************************************************************************80
!
!! ARCSIN_PDF evaluates the Arcsin PDF.
!
!  Discussion:
!
!    The LOGISTIC EQUATION has the form:
!
!      X(N+1) = 4.0D+00 * LAMBDA * ( 1.0D+00 - X(N) ).
!
!    where 0 < LAMBDA <= 1.  This nonlinear difference equation maps
!    the unit interval into itself, and is a simple example of a system
!    exhibiting chaotic behavior.  Ulam and von Neumann studied the
!    logistic equation with LAMBDA = 1, and showed that iterates of the
!    function generated a sequence of pseudorandom numbers with
!    the Arcsin probability density function.
!
!    The derived sequence
!
!      Y(N) = ( 2 / PI ) * Arcsin ( SQRT ( X(N) ) )
!
!    is a pseudorandom sequence with the uniform probability density
!    function on [0,1].  For certain starting values, such as X(0) = 0, 0.75,
!    or 1.0D+00, the sequence degenerates into a constant sequence, and for
!    values very near these, the sequence takes a while before becoming
!    chaotic.
!
!    The formula is:
!
!      PDF(X) = 1 / ( pi * sqrt ( A^2 - X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, Stephen Kokoska,
!    CRC Standard Probability and Statistics Tables and Formulae,
!    Chapman and Hall/CRC, 2000, pages 114-115.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    -A < X < A.
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d arcsin_pdf_xmm2r8(const __m128d x,
		                                const __m128d a) {

                           const __m128d pi    = _mm_set1_pd(3.14159265358979323846264338328);
			   const __m128d _0    = _mm_setzero_pd();
			   const __m128d _1    = _mm_set1_pd(1.0);
			  // const __m128d nan   = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			   __m128d pdf,t0;
			   __mmask8 m,m1;
			   //if(__builtin_expect(_mm_cmp_pd_mask(a,_0,_CMP_LE_OQ))) {
                          //     return (nan);
			  // }
			   m  =  _mm_cmp_pd_mask(x,negate_xmm2r8(a),_CMP_LE_OQ);
			   t0 =  _mm_sqrt_pd(_mm_sub_pd(_mm_mul_pd(a,a),
			                                      _mm_mul_pd(x,x)));
			   m1 = _mm_cmp_pd_mask(x,a,_CMP_GE_OQ);
			   __mmask8 m2 = m || m1;
			   pdf = _mm_mask_blend_pd(m2,_mm_div_pd(_1,
			                                           _mm_mul_pd(pi,t0)),_0);
			   return (pdf);
			   
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128 arcsin_pdf_xmm4r4(const __m128 x,
		                                const __m128 a) {

                           const __m128 pi    = _mm_set1_ps(3.14159265358979323846264338328f);
			   const __m128 _0    = _mm_setzero_ps();
			   const __m128 _1    = _mm_set1_ps(1.0f);
			  // const __m128 nan   = _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());
			   __m128 pdf,t0;
			   __mmask8 m,m1;
			   //if(__builtin_expect(_mm_cmp_ps_mask(a,_0,_CMP_LE_OQ))) {
                          //     return (nan);
			  // }
			   m  =  _mm_cmp_ps_mask(x,negate_xmm4r4(a),_CMP_LE_OQ);
			   t0 =  _mm_sqrt_ps(_mm_sub_ps(_mm_mul_ps(a,a),
			                                      _mm_mul_ps(x,x)));
			   m1 = _mm_cmp_ps_mask(x,a,_CMP_GE_OQ);
			   const __mmask8 m2 = m || m1;
			   pdf = _mm_mask_blend_ps(m2,_mm_div_ps(_1,
			                                           _mm_mul_ps(pi,t0)),_0);
			   return (pdf);
			   
		    }

/*
!*****************************************************************************80
!
!! ARCSIN_VARIANCE returns the variance of the Arcsin PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!		    
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                       __m128d arcsin_variance_xmm2r8(const __m128d a) {

                         const __m128d _1_2 = _mm_set1_pd(0.5);
			 __m128d variance;
			 variance = _mm_mul_pd(a,_mm_mul_pd(a,_1_2));
			 return (variance);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                       __m128 arcsin_variance_xmm4r4(const __m128 a) {

                         const __m128 _1_2 = _mm_set1_ps(0.5f);
			 __m128 variance;
			 variance = _mm_mul_ps(a,_mm_mul_ps(a,_1_2));
			 return (variance);
		     }

/*
!*****************************************************************************80
!
!! ARCSIN_SAMPLE samples the Arcsin PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the CDF.
!    A must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
*/
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d arcsin_sample_xmm2r8() {

                         __m128d cdf;
			 svrng_engine_t engine;
			 svrng_distribution_t uniform;
			 uint32_t seed    = 0U;
			 int32_t result   = -9999;
			 int32_t err      = -9999;
			 result           = _rdrand32_step(&seed);
			 if(!result) seed = 1043915199U;
			 engine           = svrng_new_mt19937_engine(seed);
			 err              = svrng_get_status();
			 if(err!=SVRNG_STATUS_OK) {
                            const __m128d nan = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_uniform_distribution_double(0.0,1.0);
			 const double * __restrict ptr = (const double*)(&svrng_generate2_double(engine,uniform));
			 cdf              = arcsin_cdf_inv_xmm2r8(_mm_loadu_pd(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d arcsin_sample_xmm2r8(const __m128 cdf) {

                            return (arcsin_cdf_inv_xmm2r8(cdf));
		    }
		    
		    
/*

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d  
		      normal_01_cdf_xmm2r8(const __m128d x) {
		          
		          __m128d a1 = _mm_set1_pd(0.398942280444e+00);
		          __m128d a2 = _mm_set1_pd(0.399903438504e+00);
		          __m128d a3 = _mm_set1_pd(5.75885480458e+00);
                          __m128d a4 = _mm_set1_pd(29.8213557808e+00);
                          __m128d a5 = _mm_set1_pd(2.62433121679e+00);
                          __m128d a6 = _mm_set1_pd(48.6959930692e+00);
                          __m128d a7 = _mm_set1_pd(5.92885724438e+00);
                          __m128d b0 = _mm_set1_pd(0.398942280385e+00);
                          __m128d b1 = _mm_set1_pd(3.8052e-08);
                          __m128d b2 = _mm_set1_pd(1.00000615302e+00);
                          __m128d b3 = _mm_set1_pd(3.98064794e-04);
                          __m128d b4 = _mm_set1_pd(1.98615381364e+00);
                          __m128d b5 = _mm_set1_pd(0.151679116635e+00);
                          __m128d b6 = _mm_set1_pd(5.29330324926e+00);
                          __m128d b7 = _mm_set1_pd(4.8385912808e+00);
                          __m128d b8 = _mm_set1_pd(15.1508972451e+00);
                          __m128d b9 = _mm_set1_pd(0.742380924027e+00);
                          __m128d b10= _mm_set1_pd(30.789933034e+00);
                          __m128d b11= _mm_set1_pd(3.99019417011e+00);
                          __m128d C1 = _mm_set1_pd(1.0);
                          __m128d C128 = _mm_set1_pd(1.28);
                          __m128d C05  = _mm_set1_pd(0.5);
                          __m128d C127 = _mm_set1_pd(12.7);
                          __m128d absx,y,q,cdf,t0,t1;
                          __mmask8 m0,m1,m2;
                          m2   = _mm_cmp_pd_mask(x,_mm_setzero_pd(),_CMP_LT_OQ);
                          absx = _mm_abs_pd(x);
                          m0   = _mm_cmp_pd_mask(x,C128,_CMP_LE_OQ);
                          y    = _mm_mul_pd(C05,
                                        _mm_mul_pd(x,x));
                          m1   = _mm_cmp_pd_mask(x,C127,_CMP_LE_OQ);
                          if(m0) {
                             register __m128d ya3;
                             register __m128d ya5a6
                             register __m128d ya7;
                             register __m128d a2y;
                             ya7   = _mm_add_pd(y,a7);
                             ya5a6 = _mm_add_pd(y,_mm_add_pd(a5,a6));
                             a2y   = _mm_mul_pd(a2,y);
                             ya3a4 = _mm_sub_pd(_mm_add_pd(y,a3),a4);
                             q     = _mm_sub_pd(a1,
                                           _mm_div_pd(a2y,
                                                  _mm_div_pd(ya3a4,
                                                        _mm_div_pd(ya5a6,ya7))));
                          }
                          else if(m1) {
                             register __m128d expmy;
                             register __m128d absb1;
                             register __m128d absb3;
                             register __m128d absb5;
                             register __m128d absb7;
                             register __m128d absb9;
                             register __m128d absb11;

                             expmy = _mm_mul_pd(_mm_exp_pd(negate_xmm2r8(y)),b0); 
                             absb1 = _mm_sub_pd(absx,b1);
                             absb3 = _mm_add_pd(absx,b3);
                             absb5 = _mm_sub_pd(absx,b5);
                             absb7 = _mm_add_pd(absx,b7);
                             absb9 = _mm_add_pd(absx,b9);
                             absb11= _mm_add_pd(absx,b11);
                             t0    = (absb1+b2/(absb3+b4/(absb5+b6/(absb7-b8/(absb9+b10/(absb11))))));
                             q     = _mm_div_pd(expmy,t0);
                          }
                          else {
                             q = _mm_setzero_pd();
                          }
                          
                          cdf = _mm_mask_blend_pd(m2,_mm_sub_pd(C1,q),q);
                          return (cdf);
		    }
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128  
		      normal_01_cdf_xmm4r4(const __m128 x) {
		          
		          __m128 a1 = _mm_set1_ps(0.398942280444f);
		          __m128 a2 = _mm_set1_ps(0.399903438504f);
		          __m128 a3 = _mm_set1_ps(5.75885480458f);
                          __m128 a4 = _mm_set1_ps(29.8213557808f);
                          __m128 a5 = _mm_set1_ps(2.62433121679f);
                          __m128 a6 = _mm_set1_ps(48.6959930692f);
                          __m128 a7 = _mm_set1_ps(5.92885724438f);
                          __m128 b0 = _mm_set1_ps(0.398942280385f);
                          __m128 b1 = _mm_set1_ps(3.8052e-08f);
                          __m128 b2 = _mm_set1_ps(1.00000615302f);
                          __m128 b3 = _mm_set1_ps(3.98064794e-04f);
                          __m128 b4 = _mm_set1_ps(1.98615381364f);
                          __m128 b5 = _mm_set1_ps(0.151679116635f);
                          __m128 b6 = _mm_set1_ps(5.29330324926f);
                          __m128 b7 = _mm_set1_ps(4.8385912808f);
                          __m128 b8 = _mm_set1_ps(15.1508972451f);
                          __m128 b9 = _mm_set1_ps(0.742380924027f);
                          __m128 b10= _mm_set1_ps(30.789933034f);
                          __m128 b11= _mm_set1_ps(3.99019417011f);
                          __m128 C1 = _mm_set1_ps(1.0);
                          __m128 C128 = _mm_set1_ps(1.28f);
                          __m128 C05  = _mm_set1_ps(0.5f);
                          __m128 C127 = _mm_set1_ps(12.7f);
                          __m128 absx,y,q,cdf,t0;
                          __mmask8 m0,m1,m2;
                          m2   = _mm_cmp_ps_mask(x,_mm_setzero_pd(),_CMP_LT_OQ);
                          absx = _mm_abs_ps(x);
                          m0   = _mm_cmp_ps_mask(x,C128,_CMP_LE_OQ);
                          y    = _mm_mul_ps(C05,
                                        _mm_mul_ps(x,x));
                          m1   = _mm_cmp_ps_mask(x,C127,_CMP_LE_OQ);
                          if(m0) {
                             register __m128 ya3;
                             register __m128 ya5a6
                             register __m128 ya7;
                             register __m128 a2y;
                             ya7   = _mm_add_ps(y,a7);
                             ya5a6 = _mm_add_ps(y,_mm_add_ps(a5,a6));
                             a2y   = _mm_mul_ps(a2,y);
                             ya3a4 = _mm_sub_ps(_mm_add_ps(y,a3),a4);
                             q     = _mm_sub_ps(a1,
                                           _mm_div_ps(a2y,
                                                  _mm_div_ps(ya3a4,
                                                        _mm_div_ps(ya5a6,ya7))));
                          }
                          else if(m1) {
                             register __m128 expmy;
                             register __m128 absb1;
                             register __m128 absb3;
                             register __m128 absb5;
                             register __m128 absb7;
                             register __m128 absb9;
                             register __m128 absb11;

                             expmy = _mm_mul_ps(_mm_exp_ps(negate_xmm4r4(y)),b0); 
                             absb1 = _mm_sub_ps(absx,b1);
                             absb3 = _mm_add_ps(absx,b3);
                             absb5 = _mm_sub_ps(absx,b5);
                             absb7 = _mm_add_ps(absx,b7);
                             absb9 = _mm_add_ps(absx,b9);
                             absb11= _mm_add_ps(absx,b11);
                             t0    = (absb1+b2/(absb3+b4/(absb5+b6/(absb7-b8/(absb9+b10/(absb11))))));
                             q     = _mm_div_ps(expmy,t0);
                          }
                          else {
                             q = _mm_setzero_ps();
                          }
                          
                          cdf = _mm_mask_blend_ps(m2,_mm_sub_ps(C1,q),q);
                          return (cdf);
		    }
		    
		    


/*
!*****************************************************************************80
!
!! BETA_BINOMIAL_CDF evaluates the Beta Binomial CDF.
!
!  Discussion:
!
!    A simple summing approach is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Input, integer ( kind = 4 ) C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d beta_binomial_cdf_xmm2r8(const int32_t x,
		                                       const int32_t c,
						       const __m128d a,
						       const __m128d b) {

			      const __m128d _0  = _mm_setzero_pd();
                              const __m128d _1  = _mm_set1_pd(1.0);
			      __m128d vx,vy,vcy,vc1,vy1,vcy1;
			      __m128d cdf,pdf;

			      if(x<0) {
                                 cdf = _0;
			      }
			      else if(x<c) {
                                 cdf = _0;
				 for(int32_t y = 0; y < x; ++y) {
                                     vy  = _mm_set1_pd((double)y);
				     vx  = _mm_set1_pd((double)x);
				     vcy = _mm_set1_pd((double)(c-y));
				     vc1 = _mm_set1_pd((double)(c+1));
				     vy1 = _mm_set1_pd((double)(y+1));
				     vcy1= _mm_set1_pd((double)(c-y+1));
				     const __m128d t0 = beta_xmm2r8(_mm_add_pd(a,vy),
				                                    _mm_add_pd(b,vcy));
				     const __m128d t1 = _mm_mul_pd(vc1,beta_xmm2r8(vy1,vcy1));
				     const __m128d t2 = beta_xmm2r8(a,b);
				     pdf              = _mm_div_pd(t0,_mm_mul_pd(t1,t2));
				     cdf              = _mm_add_pd(cdf,pdf);
				 }
			      }
			      else if(c<=x) {
                                  cdf = _1;
			      }
			      return (cdf);
		    }


/*!*****************************************************************************80
!
!! BETA_PDF evaluates the Beta PDF.
!
!  Discussion:
!
!    The formula for the PDF is:
!
!      PDF(A,B;X) = X**(A-1) * (1-X)**(B-1) / BETA(A,B).
!
!    A = B = 1 yields the Uniform distribution on [0,1].
!    A = B = 1/2 yields the Arcsin distribution.
!        B = 1 yields the power function distribution.
!    A = B -> Infinity tends to the Normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!*/

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d beta_pdf_xmm2r8(const __m128d x,
		                              const __m128d a,
					      const __m128d b) {

                         const __m128d _0 = _mm_setzero_pd();
			 const __m128d _1 = _mm_set1_pd(1.0);
			 const __m128d t0 = _mm_sub_pd(a,_1);
			 const __m128d t1 = _mm_sub_pd(_1,x);
			 const __m128d t2 = _mm_sub_pd(b,_1);
			 __m128d pdf,term1,term2,term3;
			 __mmask8 m0,m1,m2;
			 term1            = _mm_pow_pd(x,t0);
			 m0               = _mm_cmp_pd_mask(x,_0,_CMP_LT_OQ);
			 term2            = _mm_mul_pd(term1,_mm_pow_pd(t1,t2));
			 m1               = _mm_cmp_pd_mask(x,_1,_CMP_LT_OQ);
			 term3            = _mm_div_pd(term2,beta_xmm2r8(a,b));
			 m                = m1||m2;
			 pdf              = _mm_mask_blend_pd(m,term3,_0);
			 return (pdf);
		    }


		     __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
		      beta_pdf_xmm4r4(const __m128 x,
		                       const __m128 a,
				       const __m128 b) {

                         const __m128 _0 = _mm_setzero_ps();
			 const __m128 _1 = _mm_set1_ps(1.0);
			 const __m128 t0 = _mm_sub_ps(a,_1);
			 const __m128 t1 = _mm_sub_ps(_1,x);
			 const __m128 t2 = _mm_sub_ps(b,_1);
			 __m128 pdf,term1,term2,term3;
			 __mmask8 m0,m1,m2;
			 term1            = _mm_pow_ps(x,t0);
			 m0               = _mm_cmp_ps_mask(x,_0,_CMP_LT_OQ);
			 term2            = _mm_mul_ps(term1,_mm_pow_pd(t1,t2));
			 m1               = _mm_cmp_ps_mask(x,_1,_CMP_LT_OQ);
			 term3            = _mm_div_ps(term2,beta_xmm4r4(a,b));
			 m                = m1||m2;
			 pdf              = _mm_mask_blend_ps(m,term3,_0);
			 return (pdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d
		      beta_variance_xmm2r8(const __m128d a,
		                           const __m128d b) {

			  __m128d variance;
                          const __m128d _1  = _mm_set1_pd(1.0);
			  const __m128d ab  = _mm_add_pd(a,b);
			  const __m128d t0  = _mm_mul_pd(_mm_mul_pd(ab,ab),
			                                    _mm_add_pd(_1,ab));
			  variance          = _mm_div_pd(_mm_mul_pd(a,b),t0);				   
			  
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
		      beta_variance_xmm4r4(const __m128 a,
		                            const __m128 b) {

			  __m128 variance;
                          const __m128 _1  = _mm_set1_ps(1.0f);
			  const __m128 ab  = _mm_add_ps(a,b);
			  const __m128 t0  = _mm_mul_ps(_mm_mul_ps(ab,ab),
			                                    _mm_add_ps(_1,ab));
			  variance          = _mm_div_ps(_mm_mul_ps(a,b),t0);				   
			  
		    }

/*
!*****************************************************************************80
!
!! WEIBULL_CDF evaluates the Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    A <= X.
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!		    
*/

                      
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m128d
		      weibull_cdf_xmm2r8(const __m128d x,
		                         const __m128d a,
					 const __m128d b,
					 const __m128d c) {

                          const __m128d  _0 = _mm_setzero_pd();
			  const __m128d  _1 = _mm_set1_pd(1.0);
			  const __m128d  y  = _mm_div_pd(_mm_sub_pd(x,a),b);
			  const __m128d  exc= _mm_exp_pd(_mm_pow_pd(y,c));
			  __m128d cdf;
			  const __mmask8 m  = _mm_cmp_pd_mask(a,x,_CMP_LT_OQ);
			  cdf               = _mm_mask_blend_pd(m,_mm_sub_pd(_1,
			                                                       _mm_div_pd(_1,exc)),_0);
			  return (cdf);
		   }
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m128
		      weibull_cdf_xmm4r4(const __m128 x,
		                          const __m128 a,
					  const __m128 b,
					  const __m128 c) {

                          const __m128  _0 = _mm_setzero_ps();
			  const __m128  _1 = _mm_set1_ps(1.0f);
			  const __m128  y  = _mm_div_ps(_mm_sub_ps(x,a),b);
			  const __m128  exc= _mm_exp_ps(_mm_pow_ps(y,c));
			  __m128 cdf;
			  const __mmask8 m  = _mm_cmp_ps_mask(a,x,_CMP_LT_OQ);
			  cdf               = _mm_mask_blend_ps(m,_mm_sub_ps(_1,
			                                                       _mm_div_ps(_1,exc)),_0);
			  return (cdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m128d
		      weibull_cdf_inv_xmm2r8(const __m128d a,
		                             const __m128d b,
					     const __m128d c,
					     const __m128d cdf) {

                        const __m128d  _0  = _mm_setzero_pd();
			const __m128d  _1  = _mm_set1_pd(1.0);
			//const __m128d  nan = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			__m128d t0,t1,x;
			//if(__builtin_expect(_mm_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//   __builtin_expect(_mm_cmp_pd_mask(cdf,_1,_CMP_GT_OQ),0)) {
                       //    return (nan);
			//}
			t0                 = negate_xmm2r8(_mm_log_pd(_mm_sub_pd(_1,cdf)));
			t1                 = _mm_pow_pd(t0,_mm_div_pd(_1,c));
			x                  = _mm_fmadd_pd(a,b,t1);
			return (x);
			
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m128
		      weibull_cdf_inv_xmm4r4(const __m128 a,
		                             const __m128 b,
					     const __m128 c,
					     const __m128 cdf) {

                        const __m128  _0  = _mm_setzero_ps();
			const __m128  _1  = _mm_set1_ps(1.0f);
			//const __m128  nan = _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());
			__m128 t0,t1,x;
			//if(__builtin_expect(_mm_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//   __builtin_expect(_mm_cmp_ps_mask(cdf,_1,_CMP_GT_OQ),0)) {
                        //   return (nan);
			//}
			t0                 = negate_xmm4r4(_mm_log_pd(_mm_sub_ps(_1,cdf)));
			t1                 = _mm_pow_ps(t0,_mm_div_ps(_1,c));
			x                  = _mm_fmadd_ps(a,b,t1);
			return (x);
			
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d
		      weibull_sample_xmm2r8(const __m128d vrand,
		                            const __m128d a,
					    const __m128d b,
					    const __m128d c) {

                         return (weibull_cdf_xmm2r8(a,b,c,vrand));
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
		      weibull_sample_xmm4r4(const __m128 vrand,
		                            const __m128 a,
					    const __m128 b,
					    const __m128 c) {

                         return (weibull_cdf_xmm4r4(a,b,c,vrand));
		   }

/*
!*****************************************************************************80
!
!! WEIBULL_VARIANCE returns the variance of the Weibull PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d
                      weibull_discrete_cdf_xmm2r8(const __m128d x,
		                                  const __m128d a,
					          const __m128d b) {

			    __m128d cdf;
                            const __m128d  _0 = _mm_setzero_pd();
			    const __m128d  _1 = _mm_set1_pd(1.0);
			    const __m128d  t0 = _mm_pow_pd(_mm_add_pd(x,_1),b);
			    const __m128d  t1 = _mm_pow_pd(_mm_sub_pd(_1,a),t0);
			    const __mmask8 m  = _mm_cmp_pd_mask(x,_0,_CMP_LT_OQ);
			    cdf               = _mm_mask_blend_pd(m,_mm_sub_pd(_1,t1),_0);
			    return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
                      weibull_discrete_cdf_xmm4r4(const __m128 x,
		                                  const __m128 a,
					          const __m128 b) {

			    __m128 cdf;
                            const __m128  _0 = _mm_setzero_ps();
			    const __m128  _1 = _mm_set1_ps(1.0f);
			    const __m128  t0 = _mm_pow_ps(_mm_add_ps(x,_1),b);
			    const __m128  t1 = _mm_pow_ps(_mm_sub_ps(_1,a),t0);
			    const __mmask8 m  = _mm_cmp_ps_mask(x,_0,_CMP_LT_OQ);
			    cdf               = _mm_mask_blend_ps(m,_mm_sub_pd(_1,t1),_0);
			    return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d
		      weibull_discrete_pdf_xmm2r8(const __m128d x,
		                                  const __m128d a,
					          const __m128d b) {

                            __m128d pdf;
                            const __m128d  _0 = _mm_setzero_pd();
			    const __m128d  _1 = _mm_set1_pd(1.0);
			    const __m128d  t0 = _mm_pow_pd(_mm_add_pd(x,_1),b);
			    const __m128d  _1a= _mm_sub_pd(_1,a);
			    const __m128d  t1 = _mm_pow_pd(_1a,t0);
                            const __m128d  t2 = _mm_pow_pd(x,b);
			    const __m128d  t3 = _mm_pow_pd(_1a,t2);
			    pdf               = _mm_sub_pd(t3,t1);
			    return (pdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
		      weibull_discrete_pdf_xmm4r4(const __m128 x,
		                                  const __m128 a,
					          const __m128 b) {

                            __m128 pdf;
                            const __m128  _0 = _mm_setzero_ps();
			    const __m128  _1 = _mm_set1_ps(1.0);
			    const __m128  t0 = _mm_pow_ps(_mm_add_ps(x,_1),b);
			    const __m128  _1a= _mm_sub_ps(_1,a);
			    const __m128  t1 = _mm_pow_ps(_1a,t0);
                            const __m128  t2 = _mm_pow_ps(x,b);
			    const __m128  t3 = _mm_pow_ps(_1a,t2);
			    pdf               = _mm_sub_ps(t3,t1);
			    return (pdf);
		   }


/*
!*****************************************************************************80
!
!! WEIBULL_DISCRETE_CDF_INV inverts the Discrete Weibull CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A <= 1.0D+00,
!    0.0D+00 < B.
!
!    Output, integer ( kind = 4 ) X, the corresponding argument.
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128d
		      weibull_discr_icdf_xmm2r8(const __m128d cdf,
		                                const __m128d a,
						const __m128d b) {

                        //  const __m128d  nan = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			  const __m128d  _0  = _mm_setzero_pd();
			  const __m128d  _1  = _mm_set1_pd(1.0);
			 // if(__builtin_expect(_mm_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			 //    __builtin_expect(_mm_cmp_pd_mask(cdf,_1,_CMP_GT_OQ),0)) {
                         //  return (nan);
			 // }
			  const __m128d t0   =  _mm_log_pd(_mm_sub_pd(_1,cdf));
			  const __m128d t1   =  _mm_log_pd(_mm_sub_pd(_1,a));
			  const __m128d t2   =  _mm_div_pd(t1,t2)
			  const __m128d t3   =  _mm_pow_pd(t2,_mm_div_pd(_1,b));
			  __m128d x;
			  x                  =  _mm_ceil_pd(_mm_sub_pd(t3,_1));
			  return (x);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m128
		      weibull_discr_icdf_xmm4r4(const __m128 cdf,
		                                const __m128 a,
						const __m128 b) {

                         // const __m128  nan = _mm_set1_pd(std::numeric_limits<float>::quiet_NaN());
			  const __m128  _0  = _mm_setzero_ps();
			  const __m128  _1  = _mm_set1_ps(1.0f);
			 // if(__builtin_expect(_mm_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			 //    __builtin_expect(_mm_cmp_ps_mask(cdf,_1,_CMP_GT_OQ),0)) {
                         //  return (nan);
			 // }
			  const __m128 t0   =  _mm_log_ps(_mm_sub_ps(_1,cdf));
			  const __m128 t1   =  _mm_log_ps(_mm_sub_ps(_1,a));
			  const __m128 t2   =  _mm_div_ps(t1,t2)
			  const __m128 t3   =  _mm_pow_ps(t2,_mm_div_ps(_1,b));
			  __m128 x;
			  x                  =  _mm_ceil_ps(_mm_sub_ps(t3,_1));
			  return (x);
		    }


		    

/*
!*****************************************************************************80
!
!! WEIBULL_DISCRETE_SAMPLE samples the discrete Weibull PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 <= A <= 1.0D+00,
!    0.0D+00 < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) X, a sample of the PDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      weibull_discr_samp_xmm2r8(   const __m128d vrand,
		                                   const __m128d a,
						   const __m128d b) {

                         return (weibull_discr_icdf_xmm2r8(vrand,a,b));
		    }


/*
!*****************************************************************************80
!
!! VON_MISES_CDF evaluates the von Mises CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    Original FORTRAN77 version by Geoffrey Hill.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Geoffrey Hill,
!    Algorithm 518,
!    Incomplete Bessel Function I0: The von Mises Distribution,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 3, September 1977, pages 279-284.
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    QA276.M335
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/

#include "GMS_rotations_sse_helpers.hpp" // fmod


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      von_misses_cdf_xmm2r8(const __m128d x,
		                            const __m128d a,
					    const __m128d b) {

                        //Early exit.
			const __m128d   _0  = _mm_setzero_pd();
			const __m128d   _1  = _mm_set1_pd(1.0);
			const __m128d   pi  = _mm_set1_pd(3.14159265358979323846264338328);
			const __m128d   npi = _mm_set1_pd(-3.14159265358979323846264338328);
			const __m128d   xsa = _mm_sub_pd(x,a);
			if(__builtin_expect(_mm_cmp_pd_mask(xsa,npi,_CMP_LE_OQ),0)) {
		             return (_0);
			}
			if(__builtin_expect(_mm_cmp_pd_mask(npi,xsa,_CMP_LE_OQ),0)) {
                             return (_1); 
			}
			const __m128d  _2pi = _mm_set1_pd(6.283185307179586476925286766559);
			const __m128d  a1  = _mm_set1_pd(12.0);
			const __m128d  a2  = _mm_set1_pd(0.8);
			const __m128d  a3  = _mm_set1_pd(8.0);
			const __m128d  a4  = _mm_set1_pd(1.0);
			const __m128d  c1  = _mm_set1_pd(56.0);
			const __m128d  ck  = _mm_set1_pd(10.5);
			const __m128d  _2  = _mm_set1_pd(2.0);
			const __m128d  _1_2= _mm_set1_pd(0.5);
			__m128d arg,cdf,cn,p,r,s,sn,u,v,y,z,uprv,erfx;
			//  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
			z    = b;
			u    = _mm_castps_pd(fmod_xmm2r8(_mm_castpd_ps(_mm_add_pd(xsa,pi)),
			                                           _mm_castpd_ps(_2pi)));
			uprv = u;
			const __mmask8 m = _mm_cmp_pd_mask(u,_0,_CMP_LT_OQ);
			u    = _mm_add_pd(u,_2pi);
			u    = _mm_mask_blend_pd(m,uprv,u);
			y    = _mm_sub_pd(u,pi);
			
			//For small B, sum IP terms by backwards recursion.
			// Can not be vectorized manually, hence 0 is returned.
			// Only large B is computed.
			/*
                              This scalar code can not be vectorized.
                              ip = int ( z * a2 - a3 / ( z + a4 ) + a1 )
                              Used as loop control variable
                              do n = 2, ip
                         */
                        if(_mm_cmp_pd_mask(z,ck,_CMP_LE_OQ)) {
                           return (_0);
			}
			else {
                           const __m128d t0 = _mm_set1_pd(24.0);
			   const __m128d t1 = _mm_set1_pd(54.0);
			   const __m128d t2 = _mm_set1_pd(347.0);
			   const __m128d t3 = _mm_set1_pd(26.0);
			   const __m128d t4 = _mm_set1_pd(6.0);
			   const __m128d t5 = _mm_set1_pd(12.0);
			   const __m128d t6 = _mm_set1_pd(3.0);
			   const __m128d t7 = _mm_set1_pd(16.0);
			   const __m128d t8 = _mm_set1_pd(1.75);
			   const __m128d t9 = _mm_set1_pd(83.5);
			   c                = _mm_mul_pd(t0,z);
			   v                = _mm_sub_pd(c,c1);
			   const __m128d tmp1 = _mm_sub_pd(_mm_add_pd(v,t3),c);
			   const __m128d tmp2 = _mm_div_pd(t1,_mm_div_pd(t2,tmp1));
			   const __m128d tmp3 = _mm_add_pd(_mm_sub_pd(tmp2,t4),c);
			   r                  = _mm_sqrt_pd(_mm_div_pd(tmp3,t5));

			   z                  = _mm_mul_pd(_mm_sin_pd(
			                                                _mm_mul_pd(_1_2,y)),r);
                           s                  = _mm_mul_pd(_2,_mm_mul_pd(z,z));
			   v                  = _mm_sub_pd(v,_mm_add_pd(s,t6));
			   y                  = _mm_div_pd(_mm_sub_pd(_mm_sub_pd(c,s),
			                                                    _mm_sub_pd(s,t7)),t6);
			   tmp1               = _mm_sub_pd(v,y);
			   y                  = _mm_div_pd(_mm_fmadd_pd(_mm_add_pd(s,t8),s,t9),tmp1);
			   tmp2               = _mm_mul_pd(y,y);
			   arg                = _mm_mul_pd(z,_mm_sub_pd(_1,
			                                                  _mm_div_pd(s,tmp2)));
			   erfx               = _mm_erf_pd(arg);
			   cdf                = _mm_fmadd_pd(_1_2,erfx,_1_2);
			}
			cdf                   = _mm_max_pd(cdf,_0);
			cdf                   = _mm_min_pd(cdf,_1);
			return (cdf);
			
		   }


		      
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      von_misses_cdf_xmm4r4(const __m128 x,
		                            const __m128 a,
					    const __m128 b) {

                        //Early exit.
			const __m128   _0  = _mm_setzero_ps();
			const __m128   _1  = _mm_set1_ps(1.0f);
			const __m128   pi  = _mm_set1_ps(3.14159265358979323846264338328f);
			const __m128   npi = _mm_set1_ps(-3.14159265358979323846264338328f);
			const __m128   xsa = _mm_sub_ps(x,a);
			if(__builtin_expect(_mm_cmp_ps_mask(xsa,npi,_CMP_LE_OQ),0)) {
		             return (_0);
			}
			if(__builtin_expect(_mm_cmp_ps_mask(npi,xsa,_CMP_LE_OQ),0)) {
                             return (_1); 
			}
			const __m128  _2pi = _mm_set1_ps(6.283185307179586476925286766559f);
			const __m128  a1  = _mm_set1_ps(12.0f);
			const __m128  a2  = _mm_set1_ps(0.8f);
			const __m128  a3  = _mm_set1_ps(8.0f);
			const __m128  a4  = _mm_set1_ps(1.0f);
			const __m128  c1  = _mm_set1_ps(56.0f);
			const __m128  ck  = _mm_set1_ps(10.5f);
			const __m128  _2  = _mm_set1_ps(2.0f);
			const __m128  _1_2= _mm_set1_ps(0.5f);
			__m128 arg,cdf,cn,p,r,s,sn,u,v,y,z,uprv,erfx;
			//  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
			z    = b;
			u    = fmod_xmm4r4(_mm_add_ps(xsa,pi),_2pi);
			uprv = u;
			const __mmask8 m = _mm_cmp_ps_mask(u,_0,_CMP_LT_OQ);
			u    = _mm_add_ps(u,_2pi);
			u    = _mm_mask_blend_ps(m,uprv,u);
			y    = _mm_sub_ps(u,pi);
			
			//For small B, sum IP terms by backwards recursion.
			// Can not be vectorized manually, hence 0 is returned.
			// Only large B is computed.
			/*
                              This scalar code can not be vectorized.
                              ip = int ( z * a2 - a3 / ( z + a4 ) + a1 )
                              Used as loop control variable
                              do n = 2, ip
                         */
                        if(_mm_cmp_ps_mask(z,ck,_CMP_LE_OQ)) {
                           return (_0);
			}
			else {
                           const __m128 t0 = _mm_set1_ps(24.0f);
			   const __m128 t1 = _mm_set1_ps(54.0f);
			   const __m128 t2 = _mm_set1_ps(347.0f);
			   const __m128 t3 = _mm_set1_ps(26.0f);
			   const __m128 t4 = _mm_set1_ps(6.0f);
			   const __m128 t5 = _mm_set1_ps(12.0f);
			   const __m128 t6 = _mm_set1_ps(3.0f);
			   const __m128 t7 = _mm_set1_ps(16.0f);
			   const __m128 t8 = _mm_set1_ps(1.75f);
			   const __m128 t9 = _mm_set1_ps(83.5f);
			   c                = _mm_mul_ps(t0,z);
			   v                = _mm_sub_ps(c,c1);
			   const __m128 tmp1 = _mm_sub_ps(_mm_add_ps(v,t3),c);
			   const __m128 tmp2 = _mm_div_ps(t1,_mm_div_ps(t2,tmp1));
			   const __m128 tmp3 = _mm_add_ps(_mm_sub_ps(tmp2,t4),c);
			   r                  = _mm_sqrt_ps(_mm_div_ps(tmp3,t5));

			   z                  = _mm_mul_ps(_mm_sin_ps(
			                                                _mm_mul_ps(_1_2,y)),r);
                           s                  = _mm_mul_ps(_2,_mm_mul_ps(z,z));
			   v                  = _mm_sub_ps(v,_mm_add_ps(s,t6));
			   y                  = _mm_div_ps(_mm_sub_ps(_mm_sub_ps(c,s),
			                                                    _mm_sub_ps(s,t7)),t6);
			   tmp1               = _mm_sub_ps(v,y);
			   y                  = _mm_div_ps(_mm_fmadd_ps(_mm_add_ps(s,t8),s,t9),tmp1);
			   tmp2               = _mm_mul_ps(y,y);
			   arg                = _mm_mul_ps(z,_mm_sub_ps(_1,
			                                                  _mm_div_ps(s,tmp2)));
			   erfx               = _mm_erf_ps(arg);
			   cdf                = _mm_fmadd_ps(_1_2,erfx,_1_2);
			}
			cdf                   = _mm_max_ps(cdf,_0);
			cdf                   = _mm_min_ps(cdf,_1);
			return (cdf);
			
		   }


/*
   !*****************************************************************************80
!
!! VON_MISES_PDF evaluates the von Mises PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
!
!    where:
!
!      I0(*) is the modified Bessel function of the first
!      kind of order 0.
!
!    The von Mises distribution for points on the unit circle is
!    analogous to the normal distribution of points on a line.
!    The variable X is interpreted as a deviation from the angle A,
!    with B controlling the amount of dispersion.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 160.
!
!    Donald Best, Nicholas Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
!
!    Merran Evans, Nicholas Hastings, Brian Peacock,
!    Statistical Distributions,
!    Wiley, 2000,
!    LC: QA273.6.E92, pages 189-191.
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    LC: QA276.M335
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!              
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      von_misses_pdf_xmm2r8(const __m128d x,
		                            const __m128d a,
					    const __m128d b) {
 
                           const __m128d   pi  = _mm_set1_pd(3.14159265358979323846264338328);
			   const __m128d   _2pi= _mm_set1_pd(6.283185307179586476925286766559);
			   const __m128d   _0  = _mm_setzero_pd();
			   const __m128d   _2  = _mm_set1_pd(2.0);
			   const __m128d   t0  = _mm_sub_pd(a,pi);
			   const __m128d   t1  = _mm_add_pd(a,pi);
			   __m128d pdf;
			   __mmask8 m1,m2;
			   m1                  = _mm_cmp_pd_mask(x,t0,_CMP_LT_OQ);
			   pdf                 = _mm_mask_blend_pd(m1,_0,_0);
			   m2                  = _mm_cmp_pd_mask(x,t1,_CMP_LE_OQ);

                           const __m128d tmp1  = _mm_exp_pd(_mm_mul_pd(b,
			                                              _mm_cos_pd(
								                _mm_sub_pd(x,a))));
                           
			   pdf                 = _mm_mask_blend_pd(m2,_0,_mm_div_pd(tmp1,
			                                              _mm_mul_pd(_2pi,bessesl_i0_xmm2r8(b))));
			   return (pdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      von_misses_pdf_xmm4r4(const __m128 x,
		                            const __m128 a,
					    const __m128 b) {
 
                           const __m128   pi  = _mm_set1_pd(3.14159265358979323846264338328f);
			   const __m128   _2pi= _mm_set1_pd(6.283185307179586476925286766559f);
			   const __m128   _0  = _mm_setzero_pd();
			   const __m128   _2  = _mm_set1_pd(2.0);
			   const __m128   t0  = _mm_sub_pd(a,pi);
			   const __m128   t1  = _mm_add_pd(a,pi);
			   __m128 pdf;
			   __mmask8 m1,m2;
			   m1                  = _mm_cmp_pd_mask(x,t0,_CMP_LT_OQ);
			   pdf                 = _mm_mask_blend_pd(m1,_0,_0);
			   m2                  = _mm_cmp_pd_mask(x,t1,_CMP_LE_OQ);
                           const __m128 tmp1  = _mm_exp_ps(_mm_mul_pd(b,
			                                              _mm_cos_pd(
								                _mm_sub_pd(x,a))));
                           
			   pdf                 = _mm_mask_blend_pd(m2,_0,_mm_div_pd(tmp1,
			                                              _mm_mul_pd(_2pi,bessesl_i0_xmm4r4(b))));
			   return (pdf);
		   }

/*
!*****************************************************************************80
!
!! VON_MISES_SAMPLE samples the von Mises PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Best, Nicholas Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0D+00 <= B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
	              
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
                      von_misses_sample_xmm2r8(const __m128d a,
		                               const __m128d b) {

                          const __m128d  pi   = _mm_set1_pd(3.14159265358979323846264338328);
			  const __m128d  _1   = _mm_set1_pd(1.0);
			  const __m128d  _2   = _mm_set1_pd(2.0);
			  const __m128d  _4   = _mm_set1_pd(4.0);
			  const __m128d  _1_2 = _mm_set1_pd(0.5);
			  __m128d c,f,rho,tau,u1,r;
			  __m128d u2,u3,x,z;
			  __m128d t0,t1,t2;
			  svrng_engine_t engine;
			  svrng_distribution_t uniform;
			  uint32_t seed    = 0U;
			  int32_t result   = -9999;
			  int32_t err      = -9999;
			  result           = _rdrand32_step(&seed);
			  if(!result) seed = 1563548129U;
			  engine           = svrng_new_mt19937_engine(seed);
			  err              = svrng_get_status();
			  if(err!=SVRNG_STATUS_OK) {
                             const __m128d nan = _mm_set1_pd(std::numeric_limits<double>::quiet_NaN());
			     return (nan);
			  }
			  uniform             = svrng_new_normal_distribution_double(0.0,1.0);
			  t0                  = _mm_fmadd_pd(_4,_mm_mul_pd(b,b),_1);
			  tau                 = _mm_add_pd(_1,_mm_sqrt_pd(t0));
			  t1                  = _mm_add_pd(b,b);
			  rho                 = _mm_div_pd(_mm_sub_pd(tau,
			                                                _mm_sqrt_pd(_mm_add_pd(tau,tau))),t1);
			  t2                  = _mm_fmadd_pd(rho,rho,_1);
			  r                   = _mm_div_pd(t2,_mm_add_pd(rho,rho));
            
 			 while(true) {
                               
                              const double * __restrict ptr = (const double*)(&svrng_generate2_double(engine,uniform));
                              u1                            = _mm_loadu_pd(&ptr[0]);

                              z                             = _mm_cos_pd(_mm_mul_pd(pi,u1));
                              f                             = _mm_div_pd(_mm_fmadd_pd(r,z,_1),
			                                                    _mm_add_pd(r,z));
			      c                             = _mm_mul_pd(b,_mm_sub_pd(r,f));
			      t0                            = _mm_mul_pd(c,_mm_sub_pd(_2,c));
			                       
			      if(_mm_cmp_mask_pd(u2,t0,_CMP_LT_OQ)) break;
			      t1                            = _mm_add_pd(_mm_log_pd(
			                                                  _mm_div_pd(c,u2)),_1);
			      if(_mm_cmp_mask_pd(c,t1,_CMP_LE_OQ)) break;
			 }
			 const double * __restrict ptr2 =
			                    (const double*)(&svrng_generate2_double(engine,uniform));
			 u3                             = _mm_loadu_pd(&ptr2[0]);
		         t2                             = xmm2r8_sign_xmm2r8(_1,_mm_sub_pd(u3,_1_2));
			 x                              = _mm_fmadd_pd(t2,_mm_acos_pd(f),a);
			 svrng_delete_engine(engine);
			 return (x)
		   }


		   
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
                      von_misses_sample_xmm4r4(const __m128 a,
		                                const __m128 b) {

                          const __m128   pi   = _mm_set1_ps(3.14159265358979323846264338328f);
			  const __m128   _1   = _mm_set1_ps(1.0f);
			  const __m128  _2    = _mm_set1_ps(2.0f);
			  const __m128  _4    = _mm_set1_ps(4.0f);
			  const __m128  _1_2  = _mm_set1_ps(0.5f);
			  __m128 c,f,rho,tau,u1,r;
			  __m128 u2,u3,x,z;
			  __m128 t0,t1,t2;
			  svrng_engine_t engine;
			  svrng_distribution_t uniform;
			  uint32_t seed    = 0U;
			  int32_t result   = -9999;
			  int32_t err      = -9999;
			  result           = _rdrand32_step(&seed);
			  if(!result) seed = 1563548129U;
			  engine           = svrng_new_mt19937_engine(seed);
			  err              = svrng_get_status();
			  if(err!=SVRNG_STATUS_OK) {
                             const __m128 nan = _mm_set1_ps(std::numeric_limits<float>::quiet_NaN());
			     return (nan);
			  }
			  uniform             = svrng_new_normal_distribution_float(0.0f,1.0f);
			  t0                  = _mm_fmadd_ps(_4,_mm_mul_ps(b,b),_1);
			  tau                 = _mm_add_ps(_1,_mm_sqrt_ps(t0));
			  t1                  = _mm_add_ps(b,b);
			  rho                 = _mm_div_ps(_mm_sub_ps(tau,
			                                                _mm_sqrt_ps(_mm_add_ps(tau,tau))),t1);
			  t2                  = _mm_fmadd_ps(rho,rho,_1);
			  r                   = _mm_div_ps(t2,_mm_add_ps(rho,rho));
            
 			 while(true) {
                               
                              const float * __restrict ptr = (const float*)(&svrng_generate4_float(engine,uniform));
                              u1                            = _mm_loadu_ps(&ptr[0]);

                              z                             = _mm_cos_ps(_mm_mul_ps(pi,u1));
                              f                             = _mm_div_ps(_mm_fmadd_ps(r,z,_1),
			                                                    _mm_add_ps(r,z));
			      c                             = _mm_mul_ps(b,_mm_sub_ps(r,f));
			      t0                            = _mm_mul_ps(c,_mm_sub_ps(_2,c));
			                       
			      if(_mm_cmp_mask_ps(u2,t0,_CMP_LT_OQ)) break;
			      t1                            = _mm_add_ps(_mm_log_ps(
			                                                  _mm_div_ps(c,u2)),_1);
			      if(_mm_cmp_mask_ps(c,t1,_CMP_LE_OQ)) break;
			 }
			 const float * __restrict ptr2 =
			                    (const float*)(&svrng_generate4_float(engine,uniform));
			 u3                             = _mm_loadu_ps(&ptr2[0]);
		         t2                             = xmm4r4_sign_xmm4r4(_1,_mm_sub_ps(u3,_1_2));
			 x                              = _mm_fmadd_ps(t2,_mm_acos_ps(f),a);
			 svrng_delete_engine(engine);
			 return (x)
		   }

/*
!*****************************************************************************80
!
!! RAYLEIGH_PDF evaluates the Rayleigh PDF.
!
!  Discussion:
!
!    PDF(A;X) = ( X / A^2 ) * EXP ( - X^2 / ( 2 * A^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
                      
*/


                       __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      rayleigh_pdf_xmm2r8(const __m128d x,
		                          const __m128d a) {

                           const __m128d  _0 = _mm_setzero_pd();
			   __m128d t0,t1,t2,t3,pdf;
			   const __mmask8 m  = _mm_cmp_pd_mask(x,_0,_CMP_LT_OQ);
			   t0                = _mm_mul_pd(a,a);
			   t1                = negate_xmm2r8(_mm_div_pd(_mm_mul_pd(x,x),
			                                                   _mm_add_pd(t0,t0)));
			   t2                = _mm_div_pd(x,t0);
			   t3               = _mm_mul_pd(t2,_mm_exp_pd(t1));
                           pdf              = _mm_mask_blend_pd(m,t3,_0);
                           return (pdf);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      rayleigh_pdf_xmm4r4(const __m128 x,
		                           const __m128 a) {

                           const __m128  _0 = _mm_setzero_ps();
			   __m128 t0,t1,t2t3,pdf;
			   const __mmask8 m  = _mm_cmp_ps_mask(x,_0,_CMP_LT_OQ);
			   t0                = _mm_mul_ps(a,a);
			   t1                = negate_xmm4r4(_mm_div_ps(_mm_mul_ps(x,x),
			                                                   _mm_add_ps(t0,t0)));
			   t2                = _mm_div_ps(x,t0);
			   t3                = _mm_mul_ps(t2,_mm_exp_ps(t1));
                           pdf               = _mm_mask_blend_ps(m,_t3,_0);
                           return (pdf);
		     }

/*
!*****************************************************************************80
!
!! RAYLEIGH_MEAN returns the mean of the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the PDF.		     
*/


                      
      		       __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      rayleigh_mean_xmm2r8(const __m128d a) {

                          const __m128d hpi =  _mm_set1_pd(0.5*3.14159265358979323846264338328);
			  __m128d mean;
			  mean              =  _mm_mul_pd(a,_mm_sqrt_pd(hpi));
			  return (mean);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      rayleigh_mean_xmm4r4(const __m128d a) {

                          const __m128 hpi =  _mm_set1_ps(0.5f*3.14159265358979323846264338328f);
			  __m128 mean;
			  mean              =  _mm_mul_ps(a,_mm_sqrt_ps(hpi));
			  return (mean);
		   }


/*
!*****************************************************************************80
!
!! RAYLEIGH_CDF_INV inverts the Rayleigh CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      rayleigh_invcdf_xmm2r8(const __m128d cdf,
		                             const __m128d a) {

			 const __m128d _0 = _mm_setzero_pd();
			 const __m128d _1 = _mm_setzero_pd(1.0);
			 const __m128d n2 = _mm_setzero_pd(-2.0);
			 __m128d inv,t0,t1,;
                        // if(__builtin_expect(_mm_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//    __builtin_expect(_mm_cmp_pd_mask(_1,cdf,_CMP_LT_OQ),0)) {return;}
			 t0  = _mm_log_pd(_mm_sub_pd(_1,cdf));
			 t1  = _mm_mul_pd(_2,_mm_mul_pd(a,a));
                         inv = _mm_sqrt_pd(_mm_mul_pd(t0,t1));
			 return (inv);
			   
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      rayleigh_invcdf_xmm4r4(const __m128 cdf,
		                             const __m128 a) {

			 const __m128 _0 = _mm_setzero_ps();
			 const __m128 _1 = _mm_setzero_ps(1.0f);
			 const __m128 n2 = _mm_setzero_ps(-2.0f);
			 __m128 inv,t0,t1,;
                        // if(__builtin_expect(_mm_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//    __builtin_expect(_mm_cmp_ps_mask(_1,cdf,_CMP_LT_OQ),0)) {return;}
			 t0  = _mm_log_ps(_mm_sub_ps(_1,cdf));
			 t1  = _mm_mul_ps(_2,_mm_mul_ps(a,a));
                         inv = _mm_sqrt_ps(_mm_mul_ps(t0,t1));
			 return (inv);
			   
		     }


/*
!*****************************************************************************80
!
!! RAYLEIGH_CDF evaluates the Rayleigh CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!    0.0D+00 <= X.
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0.0D+00 < A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
*/


                       __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      rayleigh_cdf_xmm2r8(const __m128d x,
		                          const __m128d a) {

                         const __m128d _0 = _mm_setzero_pd();
			 const __m128d _1 = _mm_setzero_pd(1.0);
			 __m128d cdf,t0,t1;
			 t0              = _mm_mul_pd(_2,_mm_mul_pd(a,a));
			 t1              = negate_xmm2r8(_mm_mul_pd(x,x));
			 cdf             = _mm_sub_pd(_1,
			                             _mm_exp_pd(_mm_div_pd(t1,t0)));
                         return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      rayleigh_cdf_xmm4r4(const __m128 x,
		                           const __m128 a) {

                         const __m128 _0 = _mm_setzero_ps();
			 const __m128 _1 = _mm_setzero_ps(1.0f);
			 __m128 cdf,t0,t1;
			 t0              = _mm_mul_pd(_2,_mm_mul_ps(a,a));
			 t1              = negate_xmm4r4(_mm_mul_ps(x,x));
			 cdf             = _mm_sub_ps(_1,
			                             _mm_exp_ps(_mm_div_ps(t1,t0)));
                         return (cdf);
		    }


		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
		      rayleigh_sample_xmm2r8(const __m128d rand,
		                             const __m128d a) {

                          return (rayleigh_invcdf_xmm2r8(rand,a));
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
		      rayleigh_sample_xmm4r4(const __m128 rand,
		                             const __m128 a) {

                          return (rayleigh_invcdf_xmm4r4(rand,a));
		     }
		     
		          /*
!*****************************************************************************80
!
!! CAUCHY_CDF evaluates the Cauchy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
! 
      */
      
        	      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d    
                      cauchy_cdf_xmm2r8(const __m128d x,
                                        const __m128d a,
                                        const __m128d b) {
                        
                         const __m128d C314159265358979323846264 = 
                                               _mm_set1_pd(3.14159265358979323846264);
                         const __m128d C05 = _mm_set1_pd(0.5);
                         register __m128d cdf,y,t0,t1;
                         t0 = _mm_sub_pd(x,a);
                         t1 = _mm_div_pd(_mm_atan2_pd(t0,b),
                                     C314159265358979323846264);
                         cdf = _mm_add_pd(C05,t1);
                         return (cdf);  
                    }
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128    
                      cauchy_cdf_xmm4r4( const __m128 x,
                                         const __m128 a,
                                         const __m128 b) {
                        
                         const __m128 C314159265358979323846264 = 
                                               _mm_set1_ps(3.14159265358979323846264f);
                         const __m128 C05 = _mm_set1_ps(0.5f);
                         register __m128 cdf,y,t0,t1;
                         t0 = _mm_sub_ps(x,a);
                         t1 = _mm_div_ps(_mm_atan2_ps(t0,b),
                                     C314159265358979323846264);
                         cdf = _mm_add_ps(C05,t1);
                         return (cdf);  
                    }
                    
                    
/*
 !*****************************************************************************80
!
!! CAUCHY_CDF_INV inverts the Cauchy CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CDF, the value of the CDF.
!    0.0D+00 <= CDF <= 1.0.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument.
!
*/     


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d 
                      cauchy_cdf_inv_xmm2r8(const __m128d a,
                                            const __m128d b,
                                            const __m128d x) {
                           
                         const __m128d C314159265358979323846264 = 
                                               __m128_set1_pd(3.14159265358979323846264);
                         const __m128d C05 = _mm_set1_pd(0.5);    
                         register __m128d cdf,t0,t1;
                         t0 = _mm_mul_pd(C314159265358979323846264,
                                            _mm_sub_pd(cdf,C05));
                         t1 = _mm_tan_pd(t0);
                         cdf = _mm_fmadd_pd(a,b,t1);
                         return (cdf);
                   }    
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
                      cauchy_cdf_inv_xmm4r4(const __m128 a,
                                            const __m128 b,
                                            const __m128 x) {
                           
                         const __m128 C314159265358979323846264 = 
                                               __m128_set1_pd(3.14159265358979323846264f);
                         const __m128 C05 = _mm_set1_ps(0.5);    
                         register __m128 cdf,t0,t1;
                         t0 = _mm_mul_ps(C314159265358979323846264,
                                            _mm_sub_ps(cdf,C05));
                         t1 = _mm_tan_ps(t0);
                         cdf = _mm_fmadd_ps(a,b,t1);
                         return (cdf);
                   }  
                   
                   
/*
  !*****************************************************************************80
!
!! CAUCHY_PDF evaluates the Cauchy PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = 1 / ( PI * B * ( 1 + ( ( X - A ) / B )^2 ) )
!
!    The Cauchy PDF is also known as the Breit-Wigner PDF.  It
!    has some unusual properties.  In particular, the integrals for the
!    expected value and higher order moments are "singular", in the
!    sense that the limiting values do not exist.  A result can be
!    obtained if the upper and lower limits of integration are set
!    equal to +T and -T, and the limit as T=>INFINITY is taken, but
!    this is a very weak and unreliable sort of limit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
*/       


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d  
                      cauchy_pdf_xmm2r8( const __m128d x,
                                         const __m128d a,
                                         const __m128d b) {
                           
                         const __m128d C314159265358979323846264 = 
                                               __m128_set1_pd(3.14159265358979323846264);
                         const __m128d C1 = _mm_set1_pd(1.0);
                         register __m128d pdf,t0,t1,y,pib;
                         y   = _mm_div_pd(_mm_sub_pd(x,a),b);
                         pib = _mm_mul_pd(C314159265358979323846264,b);
                         t0  = _mm_fmadd_pd(y,y,C1);
                         t1  = _mm_mul_pd(pib,t0);
                         pdf = _mm_div_pd(C1,t1);
                         return (pdf);                     
                   }
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128 
                      cauchy_pdf_xmm4r4(const __m128 x,
                                         const __m128 a,
                                         const __m128 b) {
                           
                         const __m128 C314159265358979323846264 = 
                                               __m128_set1_ps(3.14159265358979323846264f);
                         const __m128 C1 = _mm_set1_ps(1.0);
                         register __m128 pdf,t0,t1,y,pib;
                         y   = _mm_div_ps(_mm_sub_ps(x,a),b);
                         pib = _mm_mul_ps(C314159265358979323846264,b);
                         t0  = _mm_fmadd_ps(y,y,C1);
                         t1  = _mm_mul_ps(pib,t0);
                         pdf = _mm_div_ps(C1,t1);
                         return (pdf);                     
                   }
                   
                  
/*
 !*****************************************************************************80
!
!! MAXWELL_CDF evaluates the Maxwell CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    0.0D+00 <= X
!
!    Input, real ( kind = 8 ) A, the parameter of the PDF.
!    0 < A.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!                    
*/


                
                       __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d 
                      maxwell_cdf_xmm2r8(const __m128d x,
                                         const __m128d a) {
                         
                         const __m128d C15 = _mm_set1_pd(1.5);
                         register __m128d x2,cdf;
                         x2 = _mm_div_pd(x,a);
                         cdf = gamma_incomplete_xmm2r8(C15,x2);
                         return (cdf);                      
                    }      
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128 
                      maxwell_cdf_xmm4r4(const __m128 x,
                                         const __m128 a) {
                         
                         const __m128 C15 = _mm_set1_ps(1.5f);
                         register __m128 x2,cdf;
                         x2 = _mm_div_ps(x,a);
                         cdf = gamma_incomplete_xmm4r4(C15,x2);
                         return (cdf);                      
                    }                           
                   
/*
!*****************************************************************************80
!
!! TFN calculates the T function of Owen.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        Integral ( 0 <= X <= A ) e^( -H^2 * (1+X^2) / 2 ) / (1+X^2) dX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by J C Young, C E Minder.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Donald Owen,
!    Tables for computing the bivariate normal distribution,
!    Annals of Mathematical Statistics,
!    Volume 27, pages 1075-1090, 1956.
!
!    JC Young, CE Minder,
!    Algorithm AS 76,
!    An Algorithm Useful in Calculating Non-Central T and
!    Bivariate Normal Distributions,
!    Applied Statistics,
!    Volume 23, Number 3, 1974, pages 455-457.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H, A, the arguments of the T function.
!
!    Output, real ( kind = 8 ) TFN, the value of the T function.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline  
		      __m128d owen_tfunc_xmm2r8(const __m128d h,
		                                const __m128d a) {
		             
		             __attribute__((section(".rodata")))
		             __ATTR_ALIGN__(16) static __m128d  weight[10] = {
		                                _mm_set1_pd(0.666713443086881375935688098933e-01),
                                                _mm_set1_pd(0.149451349150580593145776339658e+00),
                                                _mm_set1_pd(0.219086362515982043995534934228e+00),
                                                _mm_set1_pd(0.269266719309996355091226921569e+00),
                                                _mm_set1_pd(0.295524224714752870173892994651e+00),
                                                _mm_set1_pd(0.295524224714752870173892994651e+00),
                                                _mm_set1_pd(0.269266719309996355091226921569e+00),
                                                _mm_set1_pd(0.219086362515982043995534934228e+00),
                                                _mm_set1_pd(0.149451349150580593145776339658e+00), 
                                                _mm_set1_pd(0.666713443086881375935688098933e-01)};
                           __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(16) static __m128d  xtab[10] = {
		                                _mm_set1_pd(-0.973906528517171720077964012084e+00),
                                                _mm_set1_pd(-0.865063366688984510732096688423e+00),
                                                _mm_set1_pd(-0.679409568299024406234327365115e+00), 
                                                _mm_set1_pd(-0.433395394129247190799265943166e+00), 
                                                _mm_set1_pd(-0.148874338981631210884826001130e+00), 
                                                _mm_set1_pd(0.148874338981631210884826001130e+00), 
                                                _mm_set1_pd(0.433395394129247190799265943166e+00), 
                                                _mm_set1_pd(0.679409568299024406234327365115e+00), 
                                                _mm_set1_pd(0.865063366688984510732096688423e+00), 
                                                _mm_set1_pd(0.973906528517171720077964012084e+00)};
		           
		             const __m128d twopinv = _mm_set1_pd(0.15915494309189533576888e+00);   
		             const __m128d tv1     = _mm_set1_pd(1.0e-35);
		             const __m128d tv2     = _mm_set1_pd(15.0);
		             const __m128d tv3     = tv2;
		             const __m128d tv4     = _mm_set1_pd(1.0e-5);
		             const __m128d C05     = _mm_set1_pd(0.5);
		             const __m128d C1      = _mm_set1_pd(1.0);
		             const __m128d C2      = _mm_set1_pd(2.0);
		             const __m128d C025    = _mm_set1_pd(0.25);
		             const __m128d C0      = _mm_setzero_pd();
		             __m128d x,rt,as,h1,h2,hs,t0,t1,t2;
		             __m128d tfn;
		             if(_mm_cmp_pd_mask(_mm_abs_pd(h),tv1,_CMP_LT_OQ)) {

                                tfn = _mm_mul_pd(_mm_atan_pd(a),twopinv);
	               
		             }
		             else if(_mm_cmp_pd_mask(tv2,_mm_abs_pd(h),_CMP_LT_OQ)) {
		                tfn = C0;
		             }
		             else if(_mm_cmp_pd_mask(_mm_abs_pd(a),tv1,_CMP_LT_OQ)) {
		                 tfn = C0;
		             }
		             else {
		                 hs = _mm_mul_pd(negate_xmm2r8(C05),
		                            _mm_mul_pd(h,h));
		                 h2 = a;
		                 as = _mm_mul_pd(a,a);
                                 t0 = _mm_log_pd(_mm_add_pd(C1,as));
       
                                 __mmask8 m = _mm_cmp_pd_mask(tv3,_mm_sub_pd(t0,
                                                                     _mm_mul_pd(hs,as)),_CMP_LE_OQ);
                                 if(m) {
                                    h1 = _mm_mul_pd(C05,a);
                                    as = _mm_mul_pd(C025,as);
                                    while(true) {
                                          rt = _mm_add_pd(as,C1);
                                          t0 = _mm_add_pd(h1,_mm_fmadd_pd(hs,as,_mm_log_pd(rt)));                                 
                                          t1 = _mm_sub_pd(_mm_div_pd(C1,rt),hs);
                                          t2 = _mm_mul_pd(C2,_mm_mul_pd(h1,t1));
                                          h2 = _mm_div_pd(t0,t2);
                                          as = _mm_mul_pd(h2,h2);
                                          if(_mm_cmp_pd_mask(_mm_abs_pd(
                                                          _mm_mul_pd(h2,h1),tv4,_CMP_LT_OQ))) break;
                                          h1 = h2;                                               
                                    }
                                 }
                                 rt = C0;
                                 
                                 //for(int32_t i=0; i<10; ++i) 

                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[0],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[0],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[1],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[1],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[2],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[2],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[3],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[3],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[4],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[4],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[5],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[5],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[6],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[6],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[7],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[7],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[8],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[8],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0)));
                                 x  = _mm_fmadd_pd(C05,h2,_mm_add_pd(xtab[9],C1));
                                 t0 = _mm_add_pd(C1,_mm_mul_pd(x,x));
                                 rt = _mm_add_pd(rt,_mm_mul_pd(weight[9],
                                               _mm_div_pd(_mm_exp_pd(_mm_mul_pd(hs,t0)),t0))); 


                                 t1 = _mm_mul_pd(C05,h2);
                                 tfn= _mm_mul_pd(rt,_mm_mul_pd(t1,twopinv)); 
		             }
		             return (tfn);
		    }   
		                                           
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline  
		      __m128 owen_tfunc_xmm4r4(const __m128 h,
		                                const __m128 a) {
		             
		             __attribute__((section(".rodata")))
		             __ATTR_ALIGN__(16) static __m128  weight[10] = {
		                                _mm_set1_ps(0.666713443086881375935688098933e-01f),
                                                _mm_set1_ps(0.149451349150580593145776339658e+00f),
                                                _mm_set1_ps(0.219086362515982043995534934228e+00f),
                                                _mm_set1_ps(0.269266719309996355091226921569e+00f),
                                                _mm_set1_ps(0.295524224714752870173892994651e+00f),
                                                _mm_set1_ps(0.295524224714752870173892994651e+00f),
                                                _mm_set1_ps(0.269266719309996355091226921569e+00f),
                                                _mm_set1_ps(0.219086362515982043995534934228e+00f),
                                                _mm_set1_ps(0.149451349150580593145776339658e+00f), 
                                                _mm_set1_ps(0.666713443086881375935688098933e-01f)};
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(16) static __m128  xtab[10] = {
		                                _mm_set1_ps(-0.973906528517171720077964012084e+00f),
                                                _mm_set1_ps(-0.865063366688984510732096688423e+00f),
                                                _mm_set1_ps(-0.679409568299024406234327365115e+00f), 
                                                _mm_set1_ps(-0.433395394129247190799265943166e+00f), 
                                                _mm_set1_ps(-0.148874338981631210884826001130e+00f), 
                                                _mm_set1_ps(0.148874338981631210884826001130e+00f), 
                                                _mm_set1_ps(0.433395394129247190799265943166e+00f), 
                                                _mm_set1_ps(0.679409568299024406234327365115e+00f), 
                                                _mm_set1_ps(0.865063366688984510732096688423e+00f), 
                                                _mm_set1_ps(0.973906528517171720077964012084e+00f)};
		           
		             const __m128 twopinv = _mm_set1_ps(0.15915494309189533576888e+00f);   
		             const __m128 tv1     = _mm_set1_ps(std::numeric_limits<float>::min());
		             const __m128 tv2     = _mm_set1_ps(15.0f);
		             const __m128 tv3     = tv2;
		             const __m128 tv4     = _mm_set1_ps(1.0e-5f);
		             const __m128 C05     = _mm_set1_ps(0.5f);
		             const __m128 C1      = _mm_set1_ps(1.0f);
		             const __m128 C2      = _mm_set1_ps(2.0f);
		             const __m128 C025    = _mm_set1_ps(0.25f);
		             const __m128 C0      = _mm_setzero_ps();
		             __m128 x,rt,as,h1,h2,hs,t0,t1,t2;
		             __m128 tfn;
		             if(_mm_cmp_ps_mask(_mm_abs_ps(h),tv1,_CMP_LT_OQ)) {
                                tfn = _mm_mul_ps(_mm_atan_ps(a),twopinv);
              		     }
		             else if(_mm_cmp_ps_mask(tv2,_mm_abs_ps(h),_CMP_LT_OQ)) {
		                tfn = C0;
		             }
		             else if(_mm_cmp_pd_mask(_mm_abs_ps(a),tv1,_CMP_LT_OQ)) {
		                 tfn = C0;
		             }
		             else {
		                 hs = _mm_mul_ps(negate_zmm16r4(C05),
		                            _mm_mul_ps(h,h));
		                 h2 = a;
		                 as = _mm_mul_ps(a,a);
                                 t0 = _mm_log_ps(_mm_add_ps(C1,as));
           
                                 __mmask8 m = _mm_cmp_ps_mask(tv3,_mm_sub_ps(t0,
                                                                     _mm_mul_ps(hs,as)),_CMP_LE_OQ);
                                 if(m) {
                                    h1 = _mm_mul_ps(C05,a);
                                    as = _mm_mul_ps(C025,as);
                                    while(true) {
                                          rt = _mm_add_ps(as,C1);
                                          t0 = _mm_add_ps(h1,_mm_fmadd_ps(hs,as,_mm_log_ps(rt)));
                                          t1 = _mm_sub_ps(_mm_div_ps(C1,rt),hs);
                                          t2 = _mm_mul_ps(C2,_mm_mul_ps(h1,t1));
                                          h2 = _mm_div_ps(t0,t2);
                                          as = _mm_mul_ps(h2,h2);
                                          if(_mm_cmp_ps_mask(_mm_abs_ps(
                                                          _mm_mul_ps(h2,h1),tv4,_CMP_LT_OQ))) break;
                                          h1 = h2;                                               
                                    }
                                 }
                                 rt = C0;
                                 
                                 //for(int32_t i=0; i<10; ++i) 

                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[0],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[0],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[1],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[1],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[2],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[2],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[3],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[3],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[4],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[4],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[5],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[5],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[6],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[6],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[7],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[7],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[8],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[8],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0)));
                                 x  = _mm_fmadd_ps(C05,h2,_mm_add_ps(xtab[9],C1));
                                 t0 = _mm_add_ps(C1,_mm_mul_ps(x,x));
                                 rt = _mm_add_ps(rt,_mm_mul_ps(weight[9],
                                               _mm_div_ps(_mm_exp_ps(_mm_mul_ps(hs,t0)),t0))); 


                                 t1 = _mm_mul_ps(C05,h2);
                                 tfn= _mm_mul_ps(rt,_mm_mul_ps(t1,twopinv)); 
		             }
		             return (tfn);
		    }  
		    
		    
		    /*
       !*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!    Important notice: ******For argument interval: 2.0<=x<12.0 there is a scalarization in form
!    of 8 scalar for-loops*******!!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d 
                      gamma_xmm2r8(const __m128d x) {
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(16) static __m128d  c[7] = {
                                        _mm_set1_pd(-1.910444077728e-03), 
                                        _mm_set1_pd(8.4171387781295e-04), 
                                        _mm_set1_pd(-5.952379913043012e-04), 
                                        _mm_set1_pd(7.93650793500350248e-04),
                                        _mm_set1_pd(-2.777777777777681622553e-03),
                                        _mm_set1_pd(8.333333333333333331554247e-02),
                                        _mm_set1_pd(5.7083835261e-03)};
                        __attribute__((section(".rodata")))
                        __ATTR_ALIGN__(16) static __m128d  p[8]  = {
                                        _mm_set1_pd(-1.71618513886549492533811e+00),
                                        _mm_set1_pd(2.47656508055759199108314e+01),
                                        _mm_set1_pd(-3.79804256470945635097577e+02),
                                        _mm_set1_pd(6.29331155312818442661052e+02),
                                        _mm_set1_pd(8.66966202790413211295064e+02), 
                                        _mm_set1_pd(-3.14512729688483675254357e+04), 
                                        _mm_set1_pd(-3.61444134186911729807069e+04),
                                        _mm_set1_pd(6.64561438202405440627855e+04)};   
                       __attribute__((section(".rodata")))       
                       __ATTR_ALIGN__(16) static __m128d  q[8]   = {
                                        _mm_set1_pd(-3.08402300119738975254353e+01), 
                                        _mm_set1_pd(3.15350626979604161529144e+02),
                                        _mm_set1_pd(-1.01515636749021914166146e+03),
                                        _mm_set1_pd(-3.10777167157231109440444e+03),
                                        _mm_set1_pd(2.25381184209801510330112e+04),
                                        _mm_set1_pd(4.75584627752788110767815e+03),
                                        _mm_set1_pd(-1.34659959864969306392456e+05),
                                        _mm_set1_pd(1.15132259675553483497211e+05)};
                       __m128d pi                          = 
                                        _mm_set1_pd(3.1415926535897932384626434e+00);
                       __m128d eps                         =
                                        _mm_set1_pd(2.22e-16);
                       __m128d one                         =
                                        _mm_set1_pd(1.0);
                       __m128d half                        =
                                        _mm_set1_pd(0.5);
                       __m128d sqrtpi                      =
                                        _mm_set1_pd(0.9189385332046727417803297e+00);
                       __m128d twelve                      = 
                                        _mm_set1_pd(12.0);
                       __m128d two                         =
                                        _mm_set1_pd(2.0);
                       __m128d xbig                        =
                                        _mm_set1_pd(171.624);
                       __m128d xinf                        =
                                        _mm_set1_pd(1.0e+30);
                       __m128d xminin                      =
                                        _mm_set1_pd(2.23e-308);
                       
                       __m128d zero                        =
                                        _mm_setzero_pd();
                       register __m128d res,sum,xden,xnum;
                       register __m128d y,y1,ysq,z,fact;
                       register __m128i n;
                       
                       bool     parity;
                       parity = false;
                       y      = x;
                       // Negative argument
                       if(_mm_cmp_pd_mask(y,zero,_CMP_LE_OQ)) {
                          register __m128d t0,t1;
                          y  = negate_xmm2r8(x);
                          y1 = _mm_castsi_pd(_mm_cvttpd_epu64(y));
                          res= _mm_sub_pd(y,y1);
                          if(_mm_cmp_pd_mask(res,zero,_CMP_NEQ_OQ)) {
                            
                             t0 = _mm_mul_pd(_mm_mul_pd(y1,half),two);
                             t1 = _mm_castsi_pd(_mm_cvttpd_epu64(t0));
                             if(_mm_cmp_pd_mask(y1,t1,_CMP_NEQ_OQ)) parity = true;
                             t0 = _mm_sin_pd(_mm_mul_pd(pi,res));
                             fact = _mm_div_pd(negate_xmm2r8(pi),t0);
                             y    = _mm_add_pd(y,one);
                          }
                          else {
                             res = xinf;
                             return (res);
                          }
                       }
                       // Positive argument
                       if(_mm_cmp_pd_mask(y,eps,_CMP_LT_OQ)) {
                          __mmask8 m;
                          m = _mm_cmp_pd_mask(xminin,y,_CMP_LE_OQ);
                          res = _mm_mask_blend_pd(m,xinf,_mm_div_pd(one,y));
                          return (res);
                       }
                  }
                  else if(_mm_cmp_pd_mask(y,twelve,_CMP_LT_OQ)) {
                          y1 = y;
                          // 0.0 < argument < 1.0.
                          if(_mm_cmp_pd_mask(y,one,_CMP_LT_OQ)) {
                             z = y;
                             y = _mm_add_pd(y,one);
                          }
                          else {
                             //!  1.0 < argument < 12.0.
                             //!  Reduce argument if necessary.
                             n = _mm_sub_epi64(mm_castpd_si128(y),
                                                  _mm_set1_epi64(1LL));
                             y = _mm_sub_pd(y,_mm_castsi128_pd(n));
                             z = _mm_sub_pd(y,one);
                          }
                          //  Evaluate approximation for 1.0 < argument < 2.0.
                          xnum = zero;
                          xden = one;
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[0]),z);
                          xden = _mm_fmadd_pd(xden,z,q[0]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[1]),z);
                          xden = _mm_fmadd_pd(xden,z,q[1]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[2]),z);
                          xden = _mm_fmadd_pd(xden,z,q[2]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[3]),z);
                          xden = _mm_fmadd_pd(xden,z,q[3]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[4]),z);
                          xden = _mm_fmadd_pd(xden,z,q[4]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[5]),z);
                          xden = _mm_fmadd_pd(xden,z,q[5]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[6]),z);
                          xden = _mm_fmadd_pd(xden,z,q[7]);
                          xnum = _mm_mul_pd(_mm_add_pd(xnum,p[7]),z);
                          xden = _mm_fmadd_pd(xden,z,q[7]);
                          res  = _mm_add_pd(_mm_div_pd(xnum,xden),one);
                          // Adjust result for case  0.0 < argument < 1.0.
                          if(_mm_cmp_pd_mask(y1,y,_CMP_LT_OQ)) 
                             res = _mm_div_pd(res,y1);
                          else if(_mm_cmp_pd_mask(y,y1,_CMP_LT_OQ)) {
                          //  Important notice: ******For argument interval: 2.0<=x<12.0 there is a scalarization in form
                          //  of 8 scalar for-loops*******!!
                             __ATTR_ALIGN__(16) int64_t sn[2];
                             __ATTR_ALIGN__(16) double  sres[2];
                             __ATTR_ALIGN__(16) double  sy[2];
                             __ATTR_ALIGN__(16) double  sone[2];
                             int64_t i;
                             _mm_store_si128(&sn[0],n);
                             _mm_store_pd(&sres[0],res);
                             _mm_store_pd(&sy[0],y);
                             _mm_store_pd(&sone[0],one);
                             for(i=0; i != sn[0]; ++i) {
                                 sres[0] *= sy[0];
                                 sy[0]   += sone[0];
                             }
                             for(i=0; i != sn[1]; ++i) {
                                 sres[1] *= sy[1];
                                 sy[1]   += sone[1];
                             }
                           
                            
                             res = _mm_load_pd(&sres[0]);
                             y   = _mm_load_pd(&sy[0]);
                          }
                          
                       }
                       else {
                            //  Evaluate for 12.0 <= argument.
                            if(_mm_cmp_pd_mask(y,xbig,_CMP_LE_OQ)) {
                               ysq = _mm_mul_pd(y,y);
                               sum = c[6];
                               sum = _mm_add_pd(_mm_div_pd(sum,ysq),c[0]);
                               sum = _mm_add_pd(_mm_div_pd(sum,ysq),c[1]);
                               sum = _mm_add_pd(_mm_div_pd(sum,ysq),c[2]);
                               sum = _mm_add_pd(_mm_div_pd(sum,ysq),c[3]);
                               sum = _mm_add_pd(_mm_div_pd(sum,ysq),c[4]);
                               sum = _mm_add_pd(_mm_div_pd(sum,ysq),c[5]);
                               sum = _mm_sub_pd(_mm_div_pd(sum,y),
                                                   _mm_add_pd(y,sqrtpi));
                               sum = _mm_mul_pd(_mm_add_pd(sum,
                                                         _mm_sub_pd(y,half),_mm_log_pd(y)));
                               res = _mm_exp_pd(sum);
                                                    
                            }
                            else {
                               res = xinf;
                               return (res);
                            }
                       }
                       // !  Final adjustments and return.
                       if(parity) res = negate_xmm2r8(res);
                       if(_mm_cmp_pd_mask(fact,one,_CMP_NEQ_OQ)) res = _mm_div_pd(fact,res);
                       return (res);
                  }
                    
		     
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128 
                      gamma_xmm4r4(const __m128 x) {
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(16) static __m128  c[7] = {
                                        _mm_set1_ps(-1.910444077728e-03f), 
                                        _mm_set1_ps(8.4171387781295e-04f), 
                                        _mm_set1_ps(-5.952379913043012e-04f), 
                                        _mm_set1_ps(7.93650793500350248e-04f),
                                        _mm_set1_ps(-2.777777777777681622553e-03f),
                                        _mm_set1_ps(8.333333333333333331554247e-02f),
                                        _mm_set1_ps(5.7083835261e-03f)};
                        __attribute__((section(".rodata")))
                        __ATTR_ALIGN__(16) static __m128  p[8]  = {
                                        _mm_set1_ps(-1.71618513886549492533811e+00f),
                                        _mm_set1_ps(2.47656508055759199108314e+01f),
                                        _mm_set1_ps(-3.79804256470945635097577e+02f),
                                        _mm_set1_ps(6.29331155312818442661052e+02f),
                                        _mm_set1_ps(8.66966202790413211295064e+02f), 
                                        _mm_set1_ps(-3.14512729688483675254357e+04f), 
                                        _mm_set1_ps(-3.61444134186911729807069e+04f),
                                        _mm_set1_ps(6.64561438202405440627855e+04f)}; 
                       __attribute__((section(".rodata")))         
                       __ATTR_ALIGN__(16) static __m128  q[8]   = {
                                        _mm_set1_ps(-3.08402300119738975254353e+01f), 
                                        _mm_set1_ps(3.15350626979604161529144e+02f),
                                        _mm_set1_ps(-1.01515636749021914166146e+03f),
                                        _mm_set1_ps(-3.10777167157231109440444e+03f),
                                        _mm_set1_ps(2.25381184209801510330112e+04f),
                                        _mm_set1_ps(4.75584627752788110767815e+03f),
                                        _mm_set1_ps(-1.34659959864969306392456e+05f),
                                        _mm_set1_ps(1.15132259675553483497211e+05f)};
                       __m128 pi                          = 
                                        _mm_set1_ps(3.1415926535897932384626434e+00f);
                       __m128 eps                         =
                                        _mm_set1_ps(2.22e-16f);
                       __m128 one                         =
                                        _mm_set1_ps(1.0f);
                       __m128 half                        =
                                        _mm_set1_ps(0.5f);
                       __m128 sqrtpi                      =
                                        _mm_set1_ps(0.9189385332046727417803297e+00f);
                       __m128 twelve                      = 
                                        _mm_set1_ps(12.0f);
                       __m128 two                         =
                                        _mm_set1_ps(2.0f);
                       __m128 xbig                        =
                                        _mm_set1_ps(171.624f);
                       __m128 xinf                        =
                                        _mm_set1_ps(1.0e+30f);
                       __m128 xminin                      =
                                        _mm_set1_ps(std::numeric_limits<float>::min());
                       
                       __m128 zero                        =
                                        _mm_setzero_ps();
                       register __m128 res,sum,xden,xnum;
                       register __m128 y,y1,ysq,z,fact;
                       register __m128i n;
                       
                       bool     parity;
                       parity = false;
                       y      = x;
                       // Negative argument
                       if(_mm_cmp_ps_mask(y,zero,_CMP_LE_OQ)) {
                          register __m128 t0,t1;
                          y  = negate_xmm4r4(x);
                          y1 = _mm_castsi128_ps(_mm_cvttps_epu32(y));
                          res= _mm_sub_ps(y,y1);
                          if(_mm_cmp_ps_mask(res,zero,_CMP_NEQ_OQ)) {
                            
                             t0 = _mm_mul_ps(_mm_mul_ps(y1,half),two);
                             t1 = _mm_castsi128_ps(_mm_cvttps_epu32(t0));
                             if(_mm_cmp_ps_mask(y1,t1,_CMP_NEQ_OQ)) parity = true;
                             t0 = _mm_sin_ps(_mm_mul_ps(pi,res));
                             fact = _mm_div_ps(negate_xmm4r4(pi),t0);
                             y    = _mm_add_ps(y,one);
                          }
                          else {
                             res = xinf;
                             return (res);
                          }
                       }
                       // Positive argument
                       if(_mm_cmp_ps_mask(y,eps,_CMP_LT_OQ)) {
                          __mmask8 m;
                          m = _mm_cmp_ps_mask(xminin,y,_CMP_LE_OQ);
                          res = _mm_mask_blend_ps(m,xinf,_mm_div_ps(one,y));
                          return (res);
                       }
                  }
                  else if(_mm_cmp_ps_mask(y,twelve,_CMP_LT_OQ)) {
                          y1 = y;
                          // 0.0 < argument < 1.0.
                          if(_mm_cmp_ps_mask(y,one,_CMP_LT_OQ)) {
                             z = y;
                             y = _mm_add_ps(y,one);
                          }
                          else {
                             //!  1.0 < argument < 12.0.
                             //!  Reduce argument if necessary.
                             n = _mm_sub_epi32(_mm_castps_si128(y),
                                                  _mm_set1_epi32(1));
                             y = _mm_sub_ps(y,_mm_castsi128_ps(n));
                             z = _mm_sub_ps(y,one);
                          }
                          //  Evaluate approximation for 1.0 < argument < 2.0.
                          xnum = zero;
                          xden = one;
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[0]),z);
                          xden = _mm_fmadd_ps(xden,z,q[0]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[1]),z);
                          xden = _mm_fmadd_ps(xden,z,q[1]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[2]),z);
                          xden = _mm_fmadd_ps(xden,z,q[2]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[3]),z);
                          xden = _mm_fmadd_ps(xden,z,q[3]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[4]),z);
                          xden = _mm_fmadd_ps(xden,z,q[4]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[5]),z);
                          xden = _mm_fmadd_ps(xden,z,q[5]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[6]),z);
                          xden = _mm_fmadd_ps(xden,z,q[7]);
                          xnum = _mm_mul_ps(_mm_add_ps(xnum,p[7]),z);
                          xden = _mm_fmadd_ps(xden,z,q[7]);
                          res  = _mm_add_ps(_mm_div_ps(xnum,xden),one);
                          // Adjust result for case  0.0 < argument < 1.0.
                          if(_mm_cmp_ps_mask(y1,y,_CMP_LT_OQ)) 
                             res = _mm_div_ps(res,y1);
                          else if(_mm_cmp_ps_mask(y,y1,_CMP_LT_OQ)) {
                          //  Important notice: ******For argument interval: 2.0<=x<12.0 there is a scalarization in form
                          //  of 8 scalar for-loops*******!!
                             __ATTR_ALIGN__(16) int32_t sn[4];
                             __ATTR_ALIGN__(16) float  sres[4];
                             __ATTR_ALIGN__(16) float  sy[4];
                             __ATTR_ALIGN__(16) float  sone[4];
                             int32_t i;
                             _mm_store_si128(&sn[0],n);
                             _mm_store_ps(&sres[0],res);
                             _mm_store_ps(&sy[0],y);
                             _mm_store_ps(&sone[0],one);
                             for(i=0; i != sn[0]; ++i) {
                                 sres[0] *= sy[0];
                                 sy[0]   += sone[0];
                             }
                             for(i=0; i != sn[1]; ++i) {
                                 sres[1] *= sy[1];
                                 sy[1]   += sone[1];
                             }
                             for(i=0; i != sn[2]; ++i) {
                                 sres[2] *= sy[2];
                                 sy[2]   += sone[2];
                             }
                             for(i=0; i != sn[3]; ++i) {
                                 sres[3] *= sy[3];
                                 sy[3]   += sone[3];
                             }
                            
                              
                             res = _mm_load_ps(&sres[0]);
                             y   = _mm_load_ps(&sy[0]);
                          }
                          
                       }
                       else {
                            //  Evaluate for 12.0 <= argument.
                            if(_mm_cmp_ps_mask(y,xbig,_CMP_LE_OQ)) {
                               ysq = _mm_mul_ps(y,y);
                               sum = c[6];
                               sum = _mm_add_ps(_mm_div_ps(sum,ysq),c[0]);
                               sum = _mm_add_ps(_mm_div_ps(sum,ysq),c[1]);
                               sum = _mm_add_ps(_mm_div_ps(sum,ysq),c[2]);
                               sum = _mm_add_ps(_mm_div_ps(sum,ysq),c[3]);
                               sum = _mm_add_ps(_mm_div_ps(sum,ysq),c[4]);
                               sum = _mm_add_ps(_mm_div_ps(sum,ysq),c[5]);
                               sum = _mm_sub_ps(_mm_div_ps(sum,y),
                                                   _mm_add_ps(y,sqrtpi));

                               sum = _mm_mul_ps(_mm_add_ps(sum,
                                                         _mm_sub_ps(y,half),_mm_log_ps(y)));
                               res = _mm_exp_ps(sum);


                                                    
                            }
                            else {
                               res = xinf;
                               return (res);
                            }
                       }
                       // !  Final adjustments and return.
                       if(parity) res = negate_xmm4r4(res);
                       if(_mm_cmp_ps_mask(fact,one,_CMP_NEQ_OQ)) res = _mm_div_ps(fact,res);
                       return (res);
                  }
                  
                  
 /*                 
!*****************************************************************************80
!
!! STUDENT_VARIANCE returns the variance of the central Student T PDF.
!
!  Discussion:
!
!    The variance is not defined unless 2 < C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
*/  


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d 
                      student_variance_xmm2r8(const __m128d a,
                                              const __m128d b,
                                              const __m128d c) {
                                              
                          __m128d C2 = _mm_set1_pd(2.0);     
                          register __m128d bb,t1;
                          register __m128d var;
                          bb = _mm_mul_pd(b,b);
                          t1 = _mm_sub_pd(c,C2);
                          var= _mm_mul_pd(bb,_mm_div_pd(c,t1));
                          return (var);                    
                    }   
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
                      student_variance_xmm4r4(const __m128 a,
                                              const __m128 b,
                                              const __m128 c) {
                                              
                          __m128 C2 = _mm_set1_ps(2.0f);     
                          register __m128 bb,t1;
                          register __m128 var;
                          bb = _mm_mul_ps(b,b);
                          t1 = _mm_sub_ps(c,C2);
                          var= _mm_mul_ps(bb,_mm_div_ps(c,t1));
                          return (var);                    
                    }   
                    
                    
/*
  !*****************************************************************************80
!
!! WEIBULL_PDF evaluates the Weibull PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = ( C / B ) * ( ( X - A ) / B )**( C - 1 )
!     * EXP ( - ( ( X - A ) / B )**C ).
!
!    The Weibull PDF is also known as the Frechet PDF.
!
!    WEIBULL_PDF(A,B,1;X) is the Exponential PDF.
!
!    WEIBULL_PDF(0,1,2;X) is the Rayleigh PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A <= X
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B,
!    0.0D+00 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
! 
*/


          
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d	 
                      weibull_pdf_xmm2r8(const __m128d x,
                                         const __m128d a,
                                         const __m128d b,
                                         const __m128d c) {
                        
                         register __m128d C1 = _mm_set1_pd(1.0);
                         register __m128d y,t0,pow1,t1,exp;
                         register __m128d pdf;
                         t0 = _mm_div_pd(_mm_sub_pd(x,a),b);
                         pow1 = _mm_pow_pd(t0,_mm_sub_pd(c,C1));
                         exp  = _mm_exp_pd(_mm_pow_pd(y,c));
                         t1   = _mm_div_pd(c,b);
                         pdf  = _mm_div_pd(_mm_mul_pd(t1,pow1),exp); 
                         return (pdf);     
                   }
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128	 
                      weibull_pdf_xmm4r4(const __m128 x,
                                         const __m128 a,
                                         const __m128 b,
                                         const __m128 c) {
                        
                         register __m128 C1 = _mm_set1_ps(1.0f);
                         register __m128 y,t0,pow1,t1,exp;
                         register __m128 pdf;
                         t0 = _mm_div_ps(_mm_sub_ps(x,a),b);
                         pow1 = _mm_pow_ps(t0,_mm_sub_ps(c,C1));
                         exp  = _mm_exp_ps(_mm_pow_ps(y,c));
                         t1   = _mm_div_ps(c,b);
                         pdf  = _mm_div_ps(_mm_mul_ps(t1,pow1),exp); 
                         return (pdf);     
                   }
                    
                    
 /*
     !*****************************************************************************80
!
!! TRIGAMMA calculates the TriGamma function.
!
!  Discussion:
!
!    TriGamma(x) = d^2 log ( Gamma ( x ) ) / dx^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    FORTRAN77 original version by B Schneider
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BE Schneider,
!    Algorithm AS 121:
!    Trigamma Function,
!    Applied Statistics,
!    Volume 27, Number 1, page 97-99, 1978.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the trigamma function.
!    0 < X.
!
!    Output, real ( kind = 8 ) TRIGAMMA, the value of the
!    trigamma function at X.
!     
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128d
                      trigamma_xmm2r8(const __m128d x) {
                         
                         __m128d a  = _mm_setzero_pd();
                         __m128d C1 = _mm_set1_pd(1.0);
                         __m128d C05=_mm_set1_pd(0.5);
                         __m128d b  = _mm_set1_pd(5.0);
                         __m128d b2 = _mm_set1_pd(1.0/6.0);
                         __m128d b4 = _mm_set1_pd(-1.0/30.0);
                         __m128d b6 = _mm_set1_pd(1.0/42.0);
                         __m128d b8 = _mm_set1_pd(-1.0/30.0);
                         register __m128d y,z,t0,t1;
                         register __m128d trig;
                         
                         if(_mm_cmp_pd_mask(x,a,_CMP_LE_OQ)) {
                            trig = _mm_div_pd(C1,_mm_mul_pd(x,x));
                         }
                         else {
                            z = x;
                            trig = a;
                            while(_mm_cmp_pd_mask(z,b,_CMP_LT_OQ)) {
                                  trig = _mm_add_pd(_mm_div_pd(C1,
                                                         _mm_mul_pd(z,z)))
                                  z    = _mm_add_pd(z,C1);
                            }
                            y    = _mm_div_pd(C1,_mm_mul_pd(z,z));
                            trig = trig+C05*y+(C1+y*(b2+y*(b4+y*(b6+y*b8))))/z; 
                         }
                         return (trig);
                    } 
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m128
                      trigamma_xmm4r4(const __m128 x) {
                         
                         __m128 a  = _mm_setzero_ps();
                         __m128 C1 = _mm_set1_ps(1.0f);
                         __m128 C05=_mm_set1_ps(0.5f);
                         __m128 b  = _mm_set1_ps(5.0f);
                         __m128 b2 = _mm_set1_ps(1.0f/6.0f);
                         __m128 b4 = _mm_set1_ps(-1.0f/30.0f);
                         __m128 b6 = _mm_set1_ps(1.0f/42.0f);
                         __m128 b8 = _mm_set1_ps(-1.0f/30.0f);
                         register __m128 y,z,t0,t1;
                         register __m128 trig;
                         
                         if(_mm_cmp_ps_mask(x,a,_CMP_LE_OQ)) {
                            trig = _mm_div_ps(C1,_mm_mul_ps(x,x));
                         }
                         else {
                            z = x;
                            trig = a;
                            while(_mm_cmp_ps_mask(z,b,_CMP_LT_OQ)) {
                                  trig = _mm_add_ps(_mm_div_ps(C1,
                                                         _mm_mul_ps(z,z)))
                                  z    = _mm_add_ps(z,C1);
                            }
                            y    = _mm_div_ps(C1,_mm_mul_ps(z,z));
                            trig = trig+C05*y+(C1+y*(b2+y*(b4+y*(b6+y*b8))))/z; 
                         }
                         return (trig);
                    } 
                    
                    
                    
                   
                                     
                  

      } //math

} // gms














#endif /*__GMS_PDF_CDF_SSE_HPP__*/
