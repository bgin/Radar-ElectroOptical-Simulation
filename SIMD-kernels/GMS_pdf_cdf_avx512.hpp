
#ifndef __GMS_PDF_CDF_AVX512_HPP__
#define __GMS_PDF_CDF_AVX512_HPP__ 290520221332

namespace file_info {

 const unsigned int gGMS_PDF_CDF_AVX512_MAJOR = 1U;
 const unsigned int gGMS_PDF_CDF_AVX512_MINOR = 0U;
 const unsigned int gGMS_PDF_CDF_AVX512_MICRO = 0U;
 const unsigned int gGMS_PDF_CDF_AVX512_FULLVER =
  1000U*gGMS_PDF_CDF_AVX512_MAJOR+100U*gGMS_PDF_CDF_AVX512_MINOR+10U*gGMS_PDF_CDF_AVX512_MICRO;
 const char * const pgGMS_PDF_CDF_AVX512_CREATION_DATE = "29-05-2022 13:32 +00200 (SUN 29 MAY 2022 13:32 GMT+2)";
 const char * const pgGMS_PDF_CDF_AVX512_BUILD_DATE    = __DATE__ " " __TIME__ ;
 const char * const pgGMS_PDF_CDF_AVX512_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
 const char * const pgGMS_PDF_CDF_AVX512_SYNOPSIS      = "Manually vectorized [AVX512] PDF,CDF functions"


}

#include <immintrin.h>
#include <limits>
#include "GMS_config.h"
#if (USE_SLEEF_LIB) == 1
#include "GMS_sleefsimddp.hpp"
#include "GMS_sleefsimdsp.hpp"
#endif
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
		      __m512d gamma_log_zmm8r8(const __m512d x) {
		      
                            //if(__builtin_expect(_mm512_cmp_pd_mask(x,_0,_CMP_LE_OQ),0) ||
			   //    __builtin_expect(_mm512_cmp_pd_mask(x,xbig,_CMP_GT_OQ),0)) {
                           //    return (huge);
			   // }
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(64) static __m512d c[7] = { _mm512_set1_pd(-1.910444077728E-03),
			                                     _mm512_set1_pd(8.4171387781295E-04),
                                                             _mm512_set1_pd(-5.952379913043012E-04), 
                                                             _mm512_set1_pd(7.93650793500350248E-04), 
                                                             _mm512_set1_pd(-2.777777777777681622553E-03), 
                                                             _mm512_set1_pd(8.333333333333333331554247E-02), 
                                                             _mm512_set1_pd(5.7083835261E-03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512d p1[8] = {_mm512_set1_pd(4.945235359296727046734888E+00), 
                                                             _mm512_set1_pd(2.018112620856775083915565E+02), 
                                                             _mm512_set1_pd(2.290838373831346393026739E+03), 
                                                             _mm512_set1_pd(1.131967205903380828685045E+04),
                                                             _mm512_set1_pd(2.855724635671635335736389E+04), 
                                                             _mm512_set1_pd(3.848496228443793359990269E+04), 
                                                             _mm512_set1_pd(2.637748787624195437963534E+04), 
                                                             _mm512_set1_pd(7.225813979700288197698961E+03)};
                         __attribute__((section(".rodata")))                                    
			 __ATTR_ALIGN__(64) static __m512d p2[8] = {_mm512_set1_pd(4.974607845568932035012064E+00), 
                                                             _mm512_set1_pd(5.424138599891070494101986E+02), 
                                                             _mm512_set1_pd(1.550693864978364947665077E+04), 
                                                             _mm512_set1_pd(1.847932904445632425417223E+05), 
                                                             _mm512_set1_pd(1.088204769468828767498470E+06), 
                                                             _mm512_set1_pd(3.338152967987029735917223E+06), 
                                                             _mm512_set1_pd(5.106661678927352456275255E+06), 
                                                             _mm512_set1_pd(3.074109054850539556250927E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512d p4[8] = {_mm512_set1_pd(1.474502166059939948905062E+04), 
                                                             _mm512_set1_pd(2.426813369486704502836312E+06), 
                                                             _mm512_set1_pd(1.214755574045093227939592E+08), 
                                                             _mm512_set1_pd(2.663432449630976949898078E+09), 
                                                             _mm512_set1_pd(2.940378956634553899906876E+10), 
                                                             _mm512_set1_pd(1.702665737765398868392998E+11), 
                                                             _mm512_set1_pd(4.926125793377430887588120E+11), 
                                                             _mm512_set1_pd(5.606251856223951465078242E+11)};
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(64) static __m512d q1[8] = {_mm512_set1_pd(6.748212550303777196073036E+01), 
                                                             _mm512_set1_pd(1.113332393857199323513008E+03), 
                                                             _mm512_set1_pd(7.738757056935398733233834E+03), 
                                                             _mm512_set1_pd(2.763987074403340708898585E+04), 
                                                             _mm512_set1_pd(5.499310206226157329794414E+04), 
                                                             _mm512_set1_pd(6.161122180066002127833352E+04), 
                                                             _mm512_set1_pd(3.635127591501940507276287E+04), 
                                                             _mm512_set1_pd(8.785536302431013170870835E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512d q2[8] = {_mm512_set1_pd(1.830328399370592604055942E+02),
                                                             _mm512_set1_pd(7.765049321445005871323047E+03), 
                                                             _mm512_set1_pd(1.331903827966074194402448E+05),
                                                             _mm512_set1_pd(1.136705821321969608938755E+06), 
                                                             _mm512_set1_pd(5.267964117437946917577538E+06), 
                                                             _mm512_set1_pd(1.346701454311101692290052E+07), 
                                                             _mm512_set1_pd(1.782736530353274213975932E+07), 
                                                             _mm512_set1_pd(9.533095591844353613395747E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512d q4[8] = {_mm512_set1_pd(2.690530175870899333379843E+03), 
                                                             _mm512_set1_pd(6.393885654300092398984238E+05), 
                                                             _mm512_set1_pd(4.135599930241388052042842E+07), 
                                                             _mm512_set1_pd(1.120872109616147941376570E+09), 
                                                             _mm512_set1_pd(1.488613728678813811542398E+10), 
                                                             _mm512_set1_pd(1.016803586272438228077304E+11), 
                                                             _mm512_set1_pd(3.417476345507377132798597E+11), 
                                                             _mm512_set1_pd(4.463158187419713286462081E+11)};
			    const __m512d d1     = _mm512_set1_pd(-5.772156649015328605195174E-01);
			    const __m512d d2     = _mm512_set1_pd(4.227843350984671393993777E-01);
                            const __m512d d4     = _mm512_set1_pd(1.791759469228055000094023E+00);
                            const __m512d frtbig = _mm512_set1_pd(1.42E+09);
                            const __m512d pnt68  = _mm512_set1_pd(0.6796875E+00);
			    const __m512d sqrtpi = _mm512_set1_pd(0.9189385332046727417803297E+00);
			    const __m512d xbig   = _mm512_set1_pd(4.08E+36);
			    const __m512d _0     = _mm512_setzero_pd();
			    const __m512d _1_2   = _mm512_set1_pd(0.5);
			    const __m512d _1_5   = _mm512_set1_pd(1.5);
			    const __m512d _1     = _mm512_set1_pd(1.0);
			    const __m512d _4     = _mm512_set1_pd(4.0);
			    const __m512d _2     = _mm512_set1_pd(2.0);
			    const __m512d _12    = _mm512_set1_pd(12.0);
			    const __m512d huge   = _mm512_set1_pd(std::numeric_limits<double>::max());
			    const __m512d eps    = _mm512_set1_pd(std::numeric_limits<double>::epsilon());
			    __m512d gamlog,res,xden;
			    __m512d xm1,xm2,xm4;
			    __m512d xnum,xsq,corr;
			    gamlog = _mm512_setzero_pd();
			   
			    if(_mm512_cmp_pd_mask(x,eps,_CMP_LE_OQ)) {
                               res = zmm8r8_negate(_mm512_log_pd(x));
			    }
			    else if(_mm512_cmp_pd_mask(x,_1_5,_CMP_LE_OQ)) {
                               const __mmask8 m0 = _mm512_cmp_pd_mask(x,pnt68,_CMP_LT_OQ);
			       corr = _mm512_mask_blend_pd(m0,_0,zmm8r8_negate(_mm512_log_pd(x)));
			       xm1  = _mm512_mask_blend_pd(m0,_mm512_sub_pd(
			                                                _mm512_sub_pd(x,_1_2),_1_2));

			       if(_mm512_cmp_pd_mask(x,_1_2,_CMP_LE_OQ) ||
			          _mm512_cmp_pd_mask(pnt68,x,_CMP_LE_OQ)) {
                                   xden = _1;
				   xnum = _0;
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[0]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[0]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[1]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[1]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[2]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[2]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[3]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[3]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[4]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[4]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[5]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[5]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[6]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[6]);
				   xnum = _mm512_fmadd_pd(xnum,xm1,p1[7]);
				   xden = _mm512_fmadd_pd(xden,xm1,q1[7]);
				   const __m512d t0 = _mm512_fmadd_pd(xm1,
				                                  _mm512_div_pd(xnum,xden),d1);
				   res  = _mm512_add_pd(corr,
				                    _mm512_mul_pd(xm1,t0));
				}
				else {

                                   xm2  = _mm512_sub_pd(_mm512_sub_pd(x,_1_2),_1_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[0]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[0]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[1]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[1]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[2]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[2]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[3]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[3]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[4]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[4]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[5]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[5]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[6]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[6]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[7]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[7]);
				   const __m512d t0 = _mm512_fmadd_pd(xm2,
				                                  _mm512_div_pd(xnum,xden),d2);
				   res  = _mm512_add_pd(corr,
				                    _mm512_mul_pd(xm2,t0));
				}
			    }
			    else if(_mm512_cmp_pd_mask(x,_4,_CMP_LE_OQ)) {
                                   xm2  = _mm512_sub_pd(x,_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[0]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[0]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[1]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[1]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[2]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[2]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[3]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[3]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[4]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[4]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[5]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[5]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[6]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[6]);
				   xnum = _mm512_fmadd_pd(xnum,xm2,p2[7]);
				   xden = _mm512_fmadd_pd(xden,xm2,q2[7]);
				   res  = _mm512_mul_pd(xm2,
				                    _mm512_fmadd_pd(xm2,
						                _mm512_div_pd(xnum,xden),d2));
			    }
			    else if(_mm512_cmp_pd_mask(x,_12,_CMP_LE_OQ)) {
                                   xm4  = _mm512_sub_pd(x,_4);
				   xden = zmm8r8_negate(_1);
				   xnum = _0;
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[0]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[0]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[1]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[1]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[2]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[2]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[3]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[3]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[4]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[4]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[5]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[5]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[6]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[6]);
				   xnum = _mm512_fmadd_pd(xnum,xm4,p4[7]);
				   xden = _mm512_fmadd_pd(xden,xm4,q4[7]);
				   res  = _mm512_fmadd_pd(xm4,_mm512_div_pd(xnum,xden),d4);
			    }
			    else {
                                   res  = _0;
				   if(_mm512_cmp_pd_mask(x,frtbig,_CMP_LE_OQ)) {
                                      res = c[6];
				      xsq = _mm512_mul_pd(x,x);
				      res = _mm512_add_pd(_mm512_div_pd(res,xsq),c[0]);
				      res = _mm512_add_pd(_mm512_div_pd(res,xsq),c[1]);
				      res = _mm512_add_pd(_mm512_div_pd(res,xsq),c[2]);
				      res = _mm512_add_pd(_mm512_div_pd(res,xsq),c[3]);
				      res = _mm512_add_pd(_mm512_div_pd(res,xsq),c[4]);
				      res = _mm512_add_pd(_mm512_div_pd(res,xsq),c[5]);
				   }
                                   res  = _mm512_div_pd(res,x);
				   corr = _mm512_log_pd(x);
				   res  = _mm512_sub_pd(_mm512_add_pd(res,sqrtpi),
				                        _mm512_mul_pd(_1_2,corr));
				   res  = _mm512_fmadd_pd(x,_mm512_sub_pd(corr,_1),res);
				   
			    }

			    gamlog = res;
			    return (gamlog);
			    
		  }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 gamma_log_zmm16r4(const __m512 x) {
		      
                          //  if(__builtin_expect(_mm512_cmp_ps_mask(x,_0,_CMP_LE_OQ),0) ||
			  //     __builtin_expect(_mm512_cmp_ps_mask(x,xbig,_CMP_GT_OQ),0)) {
                          //     return (huge);
			  //  }
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(64) static __m512 c[7] = { _mm512_set1_ps(-1.910444077728E-03),
			                                     _mm512_set1_ps(8.4171387781295E-04),
                                                             _mm512_set1_ps(-5.952379913043012E-04), 
                                                             _mm512_set1_ps(7.93650793500350248E-04), 
                                                             _mm512_set1_ps(-2.777777777777681622553E-03), 
                                                             _mm512_set1_ps(8.333333333333333331554247E-02), 
                                                             _mm512_set1_ps(5.7083835261E-03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512 p1[8] = {_mm512_set1_ps(4.945235359296727046734888E+00), 
                                                             _mm512_set1_ps(2.018112620856775083915565E+02), 
                                                             _mm512_set1_ps(2.290838373831346393026739E+03), 
                                                             _mm512_set1_ps(1.131967205903380828685045E+04),
                                                             _mm512_set1_ps(2.855724635671635335736389E+04), 
                                                             _mm512_set1_ps(3.848496228443793359990269E+04), 
                                                             _mm512_set1_ps(2.637748787624195437963534E+04), 
                                                             _mm512_set1_ps(7.225813979700288197698961E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512 p2[8] = {_mm512_set1_ps(4.974607845568932035012064E+00), 
                                                             _mm512_set1_ps(5.424138599891070494101986E+02), 
                                                             _mm512_set1_ps(1.550693864978364947665077E+04), 
                                                             _mm512_set1_ps(1.847932904445632425417223E+05), 
                                                             _mm512_set1_ps(1.088204769468828767498470E+06), 
                                                             _mm512_set1_ps(3.338152967987029735917223E+06), 
                                                             _mm512_set1_ps(5.106661678927352456275255E+06), 
                                                             _mm512_set1_ps(3.074109054850539556250927E+06)};
                         __attribute__((section(".rodata")))                                    
			 __ATTR_ALIGN__(64) static __m512 p4[8] = {_mm512_set1_ps(1.474502166059939948905062E+04), 
                                                             _mm512_set1_ps(2.426813369486704502836312E+06), 
                                                             _mm512_set1_ps(1.214755574045093227939592E+08), 
                                                             _mm512_set1_ps(2.663432449630976949898078E+09), 
                                                             _mm512_set1_ps(2.940378956634553899906876E+10), 
                                                             _mm512_set1_ps(1.702665737765398868392998E+11), 
                                                             _mm512_set1_ps(4.926125793377430887588120E+11), 
                                                             _mm512_set1_ps(5.606251856223951465078242E+11)};
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(64) static __m512 q1[8] = {_mm512_set1_ps(6.748212550303777196073036E+01), 
                                                             _mm512_set1_ps(1.113332393857199323513008E+03), 
                                                             _mm512_set1_ps(7.738757056935398733233834E+03), 
                                                             _mm512_set1_ps(2.763987074403340708898585E+04), 
                                                             _mm512_set1_ps(5.499310206226157329794414E+04), 
                                                             _mm512_set1_ps(6.161122180066002127833352E+04), 
                                                             _mm512_set1_ps(3.635127591501940507276287E+04), 
                                                             _mm512_set1_ps(8.785536302431013170870835E+03)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512 q2[8] = {_mm512_set1_ps(1.830328399370592604055942E+02),
                                                             _mm512_set1_ps(7.765049321445005871323047E+03), 
                                                             _mm512_set1_ps(1.331903827966074194402448E+05),
                                                             _mm512_set1_ps(1.136705821321969608938755E+06), 
                                                             _mm512_set1_ps(5.267964117437946917577538E+06), 
                                                             _mm512_set1_ps(1.346701454311101692290052E+07), 
                                                             _mm512_set1_ps(1.782736530353274213975932E+07), 
                                                             _mm512_set1_ps(9.533095591844353613395747E+06)};
                         __attribute__((section(".rodata")))
			 __ATTR_ALIGN__(64) static __m512 q4[8] = {_mm512_set1_ps(2.690530175870899333379843E+03), 
                                                             _mm512_set1_ps(6.393885654300092398984238E+05), 
                                                             _mm512_set1_ps(4.135599930241388052042842E+07), 
                                                             _mm512_set1_ps(1.120872109616147941376570E+09), 
                                                             _mm512_set1_ps(1.488613728678813811542398E+10), 
                                                             _mm512_set1_ps(1.016803586272438228077304E+11), 
                                                             _mm512_set1_ps(3.417476345507377132798597E+11), 
                                                             _mm512_set1_ps(4.463158187419713286462081E+11)};
			    const __m512 d1     = _mm512_set1_ps(-5.772156649015328605195174E-01);
			    const __m512 d2     = _mm512_set1_ps(4.227843350984671393993777E-01);
                            const __m512 d4     = _mm512_set1_ps(1.791759469228055000094023E+00);
                            const __m512 frtbig = _mm512_set1_ps(1.42E+09);
                            const __m512 pnt68  = _mm512_set1_ps(0.6796875E+00);
			    const __m512 sqrtpi = _mm512_set1_ps(0.9189385332046727417803297E+00);
			    const __m512 xbig   = _mm512_set1_ps(4.08E+36);
			    const __m512 _0     = _mm512_setzero_ps();
			    const __m512 _1_2   = _mm512_set1_ps(0.5);
			    const __m512 _1_5   = _mm512_set1_ps(1.5);
			    const __m512 _1     = _mm512_set1_ps(1.0);
			    const __m512 _4     = _mm512_set1_ps(4.0);
			    const __m512 _2     = _mm512_set1_ps(2.0);
			    const __m512 _12    = _mm512_set1_ps(12.0);
			    const __m512 huge   = _mm512_set1_ps(std::numeric_limits<float>::max());
			    const __m512 eps    = _mm512_set1_ps(std::numeric_limits<float>::epsilon());
			    __m512 gamlog,res,xden;
			    __m512 xm1,xm2,xm4;
			    __m512 xnum,xsq,corr;
			    gamlog = _mm512_setzero_ps();
			   
			    if(_mm512_cmp_ps_mask(x,eps,_CMP_LE_OQ)) {
                               res = zmm16r4_negate(_mm512_log_ps(x));
			    }
			    else if(_mm512_cmp_ps_mask(x,_1_5,_CMP_LE_OQ)) {
                               const __mmask16 m0 = _mm512_cmp_ps_mask(x,pnt68,_CMP_LT_OQ);
			       corr = _mm512_mask_blend_ps(m0,_0,zmm16r4_negate(_mm512_log_ps(x)));
			       xm1  = _mm512_mask_blend_ps(m0,_mm512_sub_ps(
			                                                _mm512_sub_ps(x,_1_2),_1_2));

			       if(_mm512_cmp_ps_mask(x,_1_2,_CMP_LE_OQ) ||
			          _mm512_cmp_ps_mask(pnt68,x,_CMP_LE_OQ)) {
                                   xden = _1;
				   xnum = _0;
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[0]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[0]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[1]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[1]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[2]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[2]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[3]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[3]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[4]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[4]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[5]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[5]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[6]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[6]);
				   xnum = _mm512_fmadd_ps(xnum,xm1,p1[7]);
				   xden = _mm512_fmadd_ps(xden,xm1,q1[7]);
				   const __m512 t0 = _mm512_fmadd_ps(xm1,
				                                  _mm512_div_ps(xnum,xden),d1);
				   res  = _mm512_add_ps(corr,
				                    _mm512_mul_ps(xm1,t0));
				}
				else {

                                   xm2  = _mm512_sub_ps(_mm512_sub_ps(x,_1_2),_1_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[0]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[0]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[1]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[1]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[2]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[2]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[3]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[3]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[4]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[4]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[5]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[5]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[6]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[6]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[7]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[7]);
				   const __m512 t0 = _mm512_fmadd_ps(xm2,
				                                  _mm512_div_ps(xnum,xden),d2);
				   res  = _mm512_add_ps(corr,
				                    _mm512_mul_ps(xm2,t0));
				}
			    }
			    else if(_mm512_cmp_ps_mask(x,_4,_CMP_LE_OQ)) {
                                   xm2  = _mm512_sub_ps(x,_2);
				   xden = _1;
				   xnum = _0;
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[0]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[0]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[1]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[1]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[2]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[2]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[3]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[3]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[4]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[4]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[5]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[5]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[6]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[6]);
				   xnum = _mm512_fmadd_ps(xnum,xm2,p2[7]);
				   xden = _mm512_fmadd_ps(xden,xm2,q2[7]);
				   res  = _mm512_mul_ps(xm2,
				                    _mm512_fmadd_ps(xm2,
						                _mm512_div_ps(xnum,xden),d2));
			    }
			    else if(_mm512_cmp_ps_mask(x,_12,_CMP_LE_OQ)) {
                                   xm4  = _mm512_sub_ps(x,_4);
				   xden = zmm16r4_negate(_1);
				   xnum = _0;
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[0]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[0]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[1]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[1]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[2]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[2]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[3]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[3]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[4]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[4]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[5]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[5]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[6]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[6]);
				   xnum = _mm512_fmadd_ps(xnum,xm4,p4[7]);
				   xden = _mm512_fmadd_ps(xden,xm4,q4[7]);
				   res  = _mm512_fmadd_ps(xm4,_mm512_div_ps(xnum,xden),d4);
			    }
			    else {
                                   res  = _0;
				   if(_mm512_cmp_ps_mask(x,frtbig,_CMP_LE_OQ)) {
                                      res = c[6];
				      xsq = _mm512_mul_ps(x,x);
				      res = _mm512_add_ps(_mm512_div_ps(res,xsq),c[0]);
				      res = _mm512_add_ps(_mm512_div_ps(res,xsq),c[1]);
				      res = _mm512_add_ps(_mm512_div_ps(res,xsq),c[2]);
				      res = _mm512_add_ps(_mm512_div_ps(res,xsq),c[3]);
				      res = _mm512_add_ps(_mm512_div_ps(res,xsq),c[4]);
				      res = _mm512_add_ps(_mm512_div_ps(res,xsq),c[5]);
				   }
                                   res  = _mm512_div_ps(res,x);
				   corr = _mm512_log_ps(x);
				   res  = _mm512_sub_ps(_mm512_add_ps(res,sqrtpi),
				                        _mm512_mul_ps(_1_2,corr));
				   res  = _mm512_fmadd_ps(x,_mm512_sub_ps(corr,_1),res);
				   
			    }

			    gamlog = res;
			    return (gamlog);
			    
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
		      __m512d  bessesl_i0_zmm8r8(const __m512d arg) {

                            __attribute__((section(".rodata")))
                            __ATTR_ALIGN__(64) static __m512d p[15] = {_mm512_set1_pd(-5.2487866627945699800E-18),
                                                                      _mm512_set1_pd(-1.5982226675653184646E-14), 
                                                                      _mm512_set1_pd(-2.6843448573468483278E-11), 
                                                                      _mm512_set1_pd(-3.0517226450451067446E-08), 
                                                                      _mm512_set1_pd(-2.5172644670688975051E-05), 
                                                                      _mm512_set1_pd(-1.5453977791786851041E-02), 
                                                                      _mm512_set1_pd(-7.0935347449210549190E+00), 
                                                                      _mm512_set1_pd(-2.4125195876041896775E+03), 
                                                                      _mm512_set1_pd(-5.9545626019847898221E+05), 
                                                                      _mm512_set1_pd(-1.0313066708737980747E+08), 
                                                                      _mm512_set1_pd(-1.1912746104985237192E+10), 
                                                                      _mm512_set1_pd(-8.4925101247114157499E+11), 
                                                                      _mm512_set1_pd(-3.2940087627407749166E+13), 
                                                                      _mm512_set1_pd(-5.5050369673018427753E+14), 
                                                                      _mm512_set1_pd(-2.2335582639474375249E+15)};
                            __attribute__((section(".rodata")))                                          
			    __ATTR_ALIGN__(64) static __m512d pp[8] = {_mm512_set1_pd(-3.9843750000000000000E-01), 
                                                                      _mm512_set1_pd(2.9205384596336793945E+00), 
                                                                      _mm512_set1_pd(-2.4708469169133954315E+00), 
                                                                      _mm512_set1_pd(4.7914889422856814203E-01), 
                                                                      _mm512_set1_pd(-3.7384991926068969150E-03), 
                                                                      _mm512_set1_pd(-2.6801520353328635310E-03), 
                                                                      _mm512_set1_pd(9.9168777670983678974E-05), 
                                                                      _mm512_set1_pd(-2.1877128189032726730E-06)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(64) static __m512d q[5]  = {_mm512_set1_pd(-3.7277560179962773046E+03), 
                                                                      _mm512_set1_pd(6.5158506418655165707E+06), 
                                                                      _mm512_set1_pd(-6.5626560740833869295E+09), 
                                                                      _mm512_set1_pd(3.7604188704092954661E+12), 
                                                                      _mm512_set1_pd(-9.7087946179594019126E+14)};
                           __attribute__((section(".rodata")))                                           
			    __ATTR_ALIGN__(64) static __m512d qq[7] = {_mm512_set1_pd(-3.1446690275135491500E+01), 
                                                                      _mm512_set1_pd(8.5539563258012929600E+01), 
                                                                      _mm512_set1_pd(-6.0228002066743340583E+01), 
                                                                      _mm512_set1_pd(1.3982595353892851542E+01), 
                                                                      _mm512_set1_pd(-1.1151759188741312645E+00), 
                                                                      _mm512_set1_pd(3.2547697594819615062E-02), 
                                                                      _mm512_set1_pd(-5.5194330231005480228E-04)};
			    const __m512d rec15                    =  _mm512_set1_pd(6.6666666666666666666E-02);
			    const __m512d xmax                     =  _mm512_set1_pd(91.9E+00);
			    const __m512d exp40                    =  _mm512_set1_pd(2.353852668370199854E+17);
			    const __m512d _1                       =  _mm512_set1_pd(1.0);
			    const __m512d _15                      =  _mm512_set1_pd(15.0);
			    const __m512d _225                     =  _mm512_set1_pd(225.0);
			    const __m512d _40                      =  _mm512_set1_pd(40.0);
			    const __m512d eps                      =  _mm512_set1_pd(std::numeric_limits<double>::epsilon());
			    const __m512d huge                     =  _mm512_set1_pd(std::mumeric_limits<double>::max());
			    __m512d value,a,b,bessel_i0;
			    __m512d sump,sumq,x,xx;
                            x = _mm512_abs_pd(arg);
			    if(_mm512_cmp_pd_mask(x,eps,_CMP_LT_OQ)) {
                               value = _1;
			    }
			    else if(_mm512_cmp_pd_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm512_mul_pd(x,x);
			       sump = p[0];
			       sump = _mm512_fmadd_pd(sump,xx,p[1]);
			       sump = _mm512_fmadd_pd(sump,xx,p[2]);
			       sump = _mm512_fmadd_pd(sump,xx,p[3]);
			       sump = _mm512_fmadd_pd(sump,xx,p[4]);
			       sump = _mm512_fmadd_pd(sump,xx,p[5]);
			       sump = _mm512_fmadd_pd(sump,xx,p[6]);
			       sump = _mm512_fmadd_pd(sump,xx,p[7]);
			       sump = _mm512_fmadd_pd(sump,xx,p[8]);
			       sump = _mm512_fmadd_pd(sump,xx,p[9]);
			       sump = _mm512_fmadd_pd(sump,xx,p[10]);
			       sump = _mm512_fmadd_pd(sump,xx,p[11]);
			       sump = _mm512_fmadd_pd(sump,xx,p[12]);
			       sump = _mm512_fmadd_pd(sump,xx,p[13]);
			       sump = _mm512_fmadd_pd(sump,xx,p[14]);
			       xx   = _mm512_sub_pd(xx,_225);
			       const __m512d xxq0 = _mm512_add_pd(xx,q[0]);
			       const __m512d xxq1 = _mm512_add_pd(xx,q[1]);
			       const __m512d xxq2 = _mm512_add_pd(xx,q[2]);
			       const __m512d xxq3 = _mm512_add_pd(xx,q[3]);
			       const __m512d xxq4 = _mm512_add_pd(xx,q[4]);
			       sumq = _mm512_mul_pd(xxq0,
			                        _mm512_mul_pd(xxq1,
						          _mm512_mul_pd(xxq2,
							            _mm512_mul_pd(xxq3,xxq4))));
			       value = _mm512_div_pd(sump,sumq);
			                                         
			    }
			    else if(_mm512_cmp_pd_mask(_15,x,_CMP_LE_OQ)) {
                                    if(_mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                                       value = huge;
				    }
				    else {
                                           xx = _mm512_sub_pd(_mm512_div_pd(_1,x),rec15);
					   const __m512d t0 = _mm512_fmadd_pd(pp[0],xx,pp[1]);
					   const __m512d c0 = _mm512_fmadd_pd(_mm512_add_pd(xx,qq[0]),xx,qq[1]);
					   const __m512d t1 = _mm512_fmadd_pd(t0,xx,pp[2]);
					   const __m512d c1 = _mm512_fmadd_pd(c0,xx,qq[2]);
					   const __m512d t2 = _mm512_fmadd_pd(t1,xx,pp[3]);
					   const __m512d c2 = _mm512_fmadd_pd(c1,xx,qq[3]);
					   const __m512d t3 = _mm512_fmadd_pd(t2,xx,pp[4]);
					   const __m512d c3 = _mm512_fmadd_pd(c2,xx,qq[4]);
					   const __m512d t4 = _mm512_fmadd_pd(t3,xx,pp[5]);
					   const __m512d c4 = _mm512_fmadd_pd(c3,xx,qq[5]);
					   const __m512d t5 = _mm512_fmadd_pd(t4,xx,pp[6]);
					   const __m512d c5 = _mm512_fmadd_pd(c4,xx,qq[6]);
					   const __m512d t6 = _mm512_fmadd_pd(t5,xx,pp[7]);
					   sump             = t6;
					   sumq             = c5;
					   value            = _mm512_div_pd(sump,sumq);
					   const __mmask8 m = _mm512_cmp_pd_mask(x,_mm512_sub_pd(xmax,_15),_CMP_LE_OQ);
					   a                = _mm512_mask_blend_pd(m,_mm512_exp_pd(_mm512_sub_pd(x,_40)),
					                                             _mm512_exp_pd(x));
					   b                = _mm512_mask_blend_pd(m,exp40,_1);
					   const __m512 tmp = _mm512_sub_pd(_mm512_mul_pd(value,a),
					                                    _mm512_mul_pd(pp[0],a));
					   value            = _mm512_mul_pd(_mm512_div_pd(tmp,_mm512_sqrt_pd(x)),b);
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
		      __m512  bessesl_i0_zmm16r4(const __m512 arg) {

                            __attribute__((section(".rodata")))
                            __ATTR_ALIGN__(64) static __m512 p[15] = {_mm512_set1_ps(-5.2487866627945699800E-18f),
                                                                      _mm512_set1_ps(-1.5982226675653184646E-14f), 
                                                                      _mm512_set1_ps(-2.6843448573468483278E-11f), 
                                                                      _mm512_set1_ps(-3.0517226450451067446E-08f), 
                                                                      _mm512_set1_ps(-2.5172644670688975051E-05f), 
                                                                      _mm512_set1_ps(-1.5453977791786851041E-02f), 
                                                                      _mm512_set1_ps(-7.0935347449210549190E+00f), 
                                                                      _mm512_set1_ps(-2.4125195876041896775E+03f), 
                                                                      _mm512_set1_ps(-5.9545626019847898221E+05f), 
                                                                      _mm512_set1_ps(-1.0313066708737980747E+08f), 
                                                                      _mm512_set1_ps(-1.1912746104985237192E+10f), 
                                                                      _mm512_set1_ps(-8.4925101247114157499E+11f), 
                                                                      _mm512_set1_ps(-3.2940087627407749166E+13f), 
                                                                      _mm512_set1_ps(-5.5050369673018427753E+14f), 
                                                                      _mm512_set1_ps(-2.2335582639474375249E+15f)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(64) static __m512 pp[8] = {_mm512_set1_ps(-3.9843750000000000000E-01f), 
                                                                      _mm512_set1_ps(2.9205384596336793945E+00f), 
                                                                      _mm512_set1_ps(-2.4708469169133954315E+00f), 
                                                                      _mm512_set1_ps(4.7914889422856814203E-01f), 
                                                                      _mm512_set1_ps(-3.7384991926068969150E-03f), 
                                                                      _mm512_set1_ps(-2.6801520353328635310E-03f), 
                                                                      _mm512_set1_ps(9.9168777670983678974E-05f), 
                                                                      _mm512_set1_ps(-2.1877128189032726730E-06f)};
                            __attribute__((section(".rodata")))
			    __ATTR_ALIGN__(64) static __m512 q[5]  = {_mm512_set1_ps(-3.7277560179962773046E+03f), 
                                                                      _mm512_set1_ps(6.5158506418655165707E+06f), 
                                                                      _mm512_set1_ps(-6.5626560740833869295E+09f), 
                                                                      _mm512_set1_ps(3.7604188704092954661E+12f), 
                                                                      _mm512_set1_ps(-9.7087946179594019126E+14f)};
                            __attribute__((section(".rodata")))                                          
			    __ATTR_ALIGN__(64) static __m512 qq[7] = {_mm512_set1_ps(-3.1446690275135491500E+01f), 
                                                                      _mm512_set1_ps(8.5539563258012929600E+01f), 
                                                                      _mm512_set1_ps(-6.0228002066743340583E+01f), 
                                                                      _mm512_set1_ps(1.3982595353892851542E+01f), 
                                                                      _mm512_set1_ps(-1.1151759188741312645E+00f), 
                                                                      _mm512_set1_ps(3.2547697594819615062E-02f), 
                                                                      _mm512_set1_ps(-5.5194330231005480228E-04f)};
			    const __m512 rec15                    =  _mm512_set1_ps(6.6666666666666666666E-02f);
			    const __m512 xmax                     =  _mm512_set1_ps(91.9E+00f);
			    const __m512 exp40                    =  _mm512_set1_ps(2.353852668370199854E+17f);
			    const __m512 _1                       =  _mm512_set1_ps(1.0f);
			    const __m512 _15                      =  _mm512_set1_ps(15.0f);
			    const __m512 _225                     =  _mm512_set1_ps(225.0f);
			    const __m512 _40                      =  _mm512_set1_ps(40.0f);
			    const __m512 eps                      =  _mm512_set1_pd(std::numeric_limits<float>::epsilon());
			    const __m512 huge                     =  _mm512_set1_pd(std::mumeric_limits<float>::max());
			    __m512 value,a,b,bessel_i0;
			    __m512 sump,sumq,x,xx;
                            x = _mm512_abs_ps(arg);
			    if(_mm512_cmp_ps_mask(x,eps,_CMP_LT_OQ)) {
                               value = _1;
			    }
			    else if(_mm512_cmp_ps_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm512_mul_ps(x,x);
			       sump = p[0];
			       sump = _mm512_fmadd_ps(sump,xx,p[1]);
			       sump = _mm512_fmadd_ps(sump,xx,p[2]);
			       sump = _mm512_fmadd_ps(sump,xx,p[3]);
			       sump = _mm512_fmadd_ps(sump,xx,p[4]);
			       sump = _mm512_fmadd_ps(sump,xx,p[5]);
			       sump = _mm512_fmadd_ps(sump,xx,p[6]);
			       sump = _mm512_fmadd_ps(sump,xx,p[7]);
			       sump = _mm512_fmadd_ps(sump,xx,p[8]);
			       sump = _mm512_fmadd_ps(sump,xx,p[9]);
			       sump = _mm512_fmadd_ps(sump,xx,p[10]);
			       sump = _mm512_fmadd_ps(sump,xx,p[11]);
			       sump = _mm512_fmadd_ps(sump,xx,p[12]);
			       sump = _mm512_fmadd_ps(sump,xx,p[13]);
			       sump = _mm512_fmadd_ps(sump,xx,p[14]);
			       xx   = _mm512_sub_ps(xx,_225);
			       const __m512 xxq0 = _mm512_add_ps(xx,q[0]);
			       const __m512 xxq1 = _mm512_add_ps(xx,q[1]);
			       const __m512 xxq2 = _mm512_add_ps(xx,q[2]);
			       const __m512 xxq3 = _mm512_add_ps(xx,q[3]);
			       const __m512 xxq4 = _mm512_add_ps(xx,q[4]);
			       sumq = _mm512_mul_ps(xxq0,
			                        _mm512_mul_ps(xxq1,
						          _mm512_mul_ps(xxq2,
							            _mm512_mul_ps(xxq3,xxq4))));
			       value = _mm512_div_ps(sump,sumq);
			                                         
			    }
			    else if(_mm512_cmp_ps_mask(_15,x,_CMP_LE_OQ)) {
                                    if(_mm512_cmp_ps_mask(xmax,x,_CMP_LT_OQ)) {
                                       value = huge;
				    }
				    else {
                                           xx = _mm512_sub_ps(_mm512_div_ps(_1,x),rec15);
					   const __m512 t0 = _mm512_fmadd_ps(pp[0],xx,pp[1]);
					   const __m512 c0 = _mm512_fmadd_ps(_mm512_add_ps(xx,qq[0]),xx,qq[1]);
					   const __m512 t1 = _mm512_fmadd_ps(t0,xx,pp[2]);
					   const __m512 c1 = _mm512_fmadd_ps(c0,xx,qq[2]);
					   const __m512 t2 = _mm512_fmadd_ps(t1,xx,pp[3]);
					   const __m512 c2 = _mm512_fmadd_ps(c1,xx,qq[3]);
					   const __m512 t3 = _mm512_fmadd_ps(t2,xx,pp[4]);
					   const __m512 c3 = _mm512_fmadd_ps(c2,xx,qq[4]);
					   const __m512 t4 = _mm512_fmadd_ps(t3,xx,pp[5]);
					   const __m512 c4 = _mm512_fmadd_ps(c3,xx,qq[5]);
					   const __m512 t5 = _mm512_fmadd_ps(t4,xx,pp[6]);
					   const __m512 c5 = _mm512_fmadd_ps(c4,xx,qq[6]);
					   const __m512 t6 = _mm512_fmadd_ps(t5,xx,pp[7]);
					   sump             = t6;
					   sumq             = c5;
					   value            = _mm512_div_ps(sump,sumq);
					   const __mmask8 m = _mm512_cmp_ps_mask(x,_mm512_sub_ps(xmax,_15),_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1
                                           a                = _mm512_mask_blend_ps(m,xexpf(_mm512_sub_ps(x,_40)),
					                                             xexpf(x));
#else
					   a                = _mm512_mask_blend_ps(m,_mm512_exp_ps(_mm512_sub_ps(x,_40)),
					                                             _mm512_exp_ps(x));
#endif     
					   b                = _mm512_mask_blend_ps(m,exp40,_1);
					   const __m512 tmp = _mm512_sub_ps(_mm512_mul_ps(value,a),
					                                    _mm512_mul_ps(pp[0],a));
					   value            = _mm512_mul_ps(_mm512_div_ps(tmp,_mm512_sqrt_ps(x)),b);
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
		      __m512d bessel_i1_zmm8r8(const __m512d arg) {

                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512d  p[15] = {_mm512_set1_pd(-1.9705291802535139930E-19), 
                                                                      _mm512_set1_pd(-6.5245515583151902910E-16), 
                                                                      _mm512_set1_pd(-1.1928788903603238754E-12), 
                                                                      _mm512_set1_pd(-1.4831904935994647675E-09), 
                                                                      _mm512_set1_pd(-1.3466829827635152875E-06), 
                                                                      _mm512_set1_pd(-9.1746443287817501309E-04), 
                                                                      _mm512_set1_pd(-4.7207090827310162436E-01), 
                                                                      _mm512_set1_pd(-1.8225946631657315931E+02), 
                                                                      _mm512_set1_pd(-5.1894091982308017540E+04), 
                                                                      _mm512_set1_pd(-1.0588550724769347106E+07), 
                                                                      _mm512_set1_pd(-1.4828267606612366099E+09), 
                                                                      _mm512_set1_pd(-1.3357437682275493024E+11), 
                                                                      _mm512_set1_pd(-6.9876779648010090070E+12), 
                                                                      _mm512_set1_pd(-1.7732037840791591320E+14), 
                                                                      _mm512_set1_pd(-1.4577180278143463643E+15)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(64) static __m512d pp[8]  = {_mm512_set1_pd(-6.0437159056137600000E-02), 
                                                                      _mm512_set1_pd(4.5748122901933459000E-01), 
                                                                      _mm512_set1_pd(-4.2843766903304806403E-01), 
                                                                      _mm512_set1_pd(9.7356000150886612134E-02), 
                                                                      _mm512_set1_pd(-3.2457723974465568321E-03), 
                                                                      _mm512_set1_pd(-3.6395264712121795296E-04), 
                                                                      _mm512_set1_pd(1.6258661867440836395E-05), 
                                                                      _mm512_set1_pd(-3.6347578404608223492E-07)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(64) static __m512d q[5]   = {_mm512_set1_pd(-4.0076864679904189921E+03), 
                                                                      _mm512_set1_pd(7.4810580356655069138E+06), 
                                                                      _mm512_set1_pd(-8.0059518998619764991E+09), 
                                                                      _mm512_set1_pd(4.8544714258273622913E+12), 
                                                                      _mm512_set1_pd(-1.3218168307321442305E+15)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(64) static __m512d qq[6]  = {_mm512_set1_pd(-3.8806586721556593450E+00), 
                                                                      _mm512_set1_pd(3.2593714889036996297E+00), 
                                                                      _mm512_set1_pd(-8.5017476463217924408E-01), 
                                                                      _mm512_set1_pd(7.4212010813186530069E-02), 
                                                                      _mm512_set1_pd(-2.2835624489492512649E-03), 
                                                                      _mm512_set1_pd(3.7510433111922824643E-05)};
			   const __m512d exp40                     =  _mm512_set1_pd(2.353852668370199854E+17);
			   const __m512d _40                       =  _mm512_set1_pd(40.0);
			   const __m512d _1_2                      =  _mm512_set1_pd(0.5);
			   const __m512d _1                        =  _mm512_set1_pd(1.0);
			   const __m512d _15                       =  _mm512_set1_pd(15.0);
			   const __m512d pbar                      =  _mm512_set1_pd(3.98437500E-01);
			   const __m512d rec15                     =  _mm512_set1_pd(6.6666666666666666666E-02);
			   const __m512d _225                      =  _mm512_set1_pd(225.0);
			   const __m512d xmax                      =  _mm512_set1_pd(713.987E+00);
			   const __m512d _0                        =  _mm512_setzero_pd();
			   const __m512d eps                       =  _mm512_set1_pd(std::numeric_limits<double>::epsilon());
			   const __m512d huge                      =  _mm512_set1_pd(std::mumeric_limits<double>::max());
			   __m512d a,b,bessel_i1,value;
			   __m512d sump,sumq,x,xx;

			   x  = _mm512_abs_pd(arg);
			   if(_mm512_cmp_pd_mask(x,eps,_CMP_LT_OQ)) {
                               value = _mm512_mul_pd(_1_2,x);
			   }
			   else if(_mm512_cmp_pd_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm512_mul_pd(x,x);
			       sump = p[0];
			       sump = _mm512_fmadd_pd(sump,xx,p[1]);
			       sump = _mm512_fmadd_pd(sump,xx,p[2]);
			       sump = _mm512_fmadd_pd(sump,xx,p[3]);
			       sump = _mm512_fmadd_pd(sump,xx,p[4]);
			       sump = _mm512_fmadd_pd(sump,xx,p[5]);
			       sump = _mm512_fmadd_pd(sump,xx,p[6]);
			       sump = _mm512_fmadd_pd(sump,xx,p[7]);
			       sump = _mm512_fmadd_pd(sump,xx,p[8]);
			       sump = _mm512_fmadd_pd(sump,xx,p[9]);
			       sump = _mm512_fmadd_pd(sump,xx,p[10]);
			       sump = _mm512_fmadd_pd(sump,xx,p[11]);
			       sump = _mm512_fmadd_pd(sump,xx,p[12]);
			       sump = _mm512_fmadd_pd(sump,xx,p[13]);
			       sump = _mm512_fmadd_pd(sump,xx,p[14]);
			       xx   = _mm512_sub_pd(xx,_225);
			       const __m512d t0 = _mm512_fmadd_pd(_mm512_add_pd(xx,q[0]),xx,q[1]);
			       const __m512d t1 = _mm512_fmadd_pd(t0,xx,q[2]);
			       const __m512d t2 = _mm512_fmadd_pd(t1,xx,q[3]);
			       const __m512d t3 = _mm512_fmadd_pd(t2,xx,q[4]);
			       sumq             = t3;
			       value            = _mm512_mul_pd(_mm512_div_pd(sump,sumq),x);
			   }
			   else if(_mm512_cmp_pd_mask(xmax,x,_CMP_LT_OQ)) {
                               value            = huge;
			   }
			   else {
                               xx               = _mm512_sub_pd(_mm512_div_pd(_1,x),rec15);
			       const __m512d t0 = _mm512_fmadd_pd(pp[0],xx,pp[1]);
			       const __m512d c0 = _mm512_fmadd_pd(_mm512_add_pd(xx,qq[0]),xx,qq[1]);
			       const __m512d t1 = _mm512_fmadd_pd(t0,xx,pp[2]);
			       const __m512d c1 = _mm512_fmadd_pd(c0,xx,qq[2]);
			       const __m512d t2 = _mm512_fmadd_pd(t1,xx,pp[3]);
			       const __m512d c2 = _mm512_fmadd_pd(c1,xx,qq[3]);
			       const __m512d t3 = _mm512_fmadd_pd(t2,xx,pp[4]);
			       const __m512d c3 = _mm512_fmadd_pd(c2,xx,qq[4]);
			       const __m512d t4 = _mm512_fmadd_pd(t3,xx,pp[5]);
			       const __m512d c4 = _mm512_fmadd_pd(c3,xx,qq[5]);
			       const __m512d t5 = _mm512_fmadd_pd(t4,xx,pp[6]);
			       const __m512d c5 = _mm512_fmadd_pd(c4,xx,qq[6]);
			       const __m512d t6 = _mm512_fmadd_pd(t5,xx,pp[7]);
			       sump             = t6;
			       sumq             = c5;
			       value            = _mm512_div_pd(sump,sumq);
			       const __mmask8 m = _mm512_cmp_pd_mask(_mm512_sub_pd(xmax,_15),_CMP_LT_OQ);
#if (USE_SLEEF_LIB) == 1
                               a                = _mm512_mask_blend_pd(m,xexp(x),
			                                                           xexp(_mm512_sub_pd(x,_40)));
#else
			       a                = _mm512_mask_blend_pd(m,_mm512_exp_pd(x),
			                                                           _mm512_exp_pd(_mm512_sub_pd(x,_40)));
#endif
			       b                = _mm512_mask_blend_pd(m,_1,_40);
			       const __m512d tmp= _mm512_add_pd(_mm512_mul_pd(value,a),
			                                        _mm512_mul_pd(pbar,a));
			       value            = _mm512_mul_pd(_mm512_div_pd(tmp,_mm512_sqrt_pd(x)),b);
			   }
			   if(_mm512_cmp_pd_mask(arg,_0,_CMP_LT_OQ)) {
                              value             = zmm8r8_negate(value);
			   }
			   bessel_i1            = value
			   return (bessel_i1);
		    }




		     __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 bessel_i1_zmm16r4(const __m512 arg) {

                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512  p[15] = {_mm512_set1_ps(-1.9705291802535139930E-19f), 
                                                                      _mm512_set1_ps(-6.5245515583151902910E-16f), 
                                                                      _mm512_set1_ps(-1.1928788903603238754E-12f), 
                                                                      _mm512_set1_ps(-1.4831904935994647675E-09f), 
                                                                      _mm512_set1_ps(-1.3466829827635152875E-06f), 
                                                                      _mm512_set1_ps(-9.1746443287817501309E-04f), 
                                                                      _mm512_set1_ps(-4.7207090827310162436E-01f), 
                                                                      _mm512_set1_ps(-1.8225946631657315931E+02f), 
                                                                      _mm512_set1_ps(-5.1894091982308017540E+04f), 
                                                                      _mm512_set1_ps(-1.0588550724769347106E+07f), 
                                                                      _mm512_set1_ps(-1.4828267606612366099E+09f), 
                                                                      _mm512_set1_ps(-1.3357437682275493024E+11f), 
                                                                      _mm512_set1_ps(-6.9876779648010090070E+12f), 
                                                                      _mm512_set1_ps(-1.7732037840791591320E+14f), 
                                                                      _mm512_set1_ps(-1.4577180278143463643E+15f)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(64) static __m512 pp[8]  = {_mm512_set1_ps(-6.0437159056137600000E-02f), 
                                                                      _mm512_set1_ps(4.5748122901933459000E-01f), 
                                                                      _mm512_set1_ps(-4.2843766903304806403E-01f), 
                                                                      _mm512_set1_ps(9.7356000150886612134E-02f), 
                                                                      _mm512_set1_ps(-3.2457723974465568321E-03f), 
                                                                      _mm512_set1_ps(-3.6395264712121795296E-04f), 
                                                                      _mm512_set1_ps(1.6258661867440836395E-05f), 
                                                                      _mm512_set1_ps(-3.6347578404608223492E-07f)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(64) static __m512 q[5]   = {_mm512_set1_ps(-4.0076864679904189921E+03f), 
                                                                      _mm512_set1_ps(7.4810580356655069138E+06f), 
                                                                      _mm512_set1_ps(-8.0059518998619764991E+09f), 
                                                                      _mm512_set1_ps(4.8544714258273622913E+12f), 
                                                                      _mm512_set1_ps(-1.3218168307321442305E+15f)};
                           __attribute__((section(".rodata")))
			   __ATTR_ALIGN__(64) static __m512 qq[6]  = {_mm512_set1_ps(-3.8806586721556593450E+00f), 
                                                                      _mm512_set1_ps(3.2593714889036996297E+00f), 
                                                                      _mm512_set1_ps(-8.5017476463217924408E-01f), 
                                                                      _mm512_set1_ps(7.4212010813186530069E-02f), 
                                                                      _mm512_set1_ps(-2.2835624489492512649E-03f), 
                                                                      _mm512_set1_ps(3.7510433111922824643E-05f)};
			   const __m512 exp40                     =  _mm512_set1_ps(2.353852668370199854E+17f);
			   const __m512 _40                       =  _mm512_set1_ps(40.0f);
			   const __m512 _1_2                      =  _mm512_set1_ps(0.5f);
			   const __m512 _1                        =  _mm512_set1_ps(1.0f);
			   const __m512 _15                       =  _mm512_set1_ps(15.0f);
			   const __m512 pbar                      =  _mm512_set1_ps(3.98437500E-01f);
			   const __m512 rec15                     =  _mm512_set1_ps(6.6666666666666666666E-02f);
			   const __m512 _225                      =  _mm512_set1_ps(225.0f);
			   const __m512 xmax                      =  _mm512_set1_ps(713.987E+00f);
			   const __m512 _0                        =  _mm512_setzero_ps();
			   const __m512 eps                       =  _mm512_set1_ps(std::numeric_limits<float>::epsilon());
			   const __m512 huge                      =  _mm512_set1_ps(std::mumeric_limits<float>::max());
			   __m512 a,b,bessel_i1,value;
			   __m512 sump,sumq,x,xx;

			   x  = _mm512_abs_ps(arg);
			   if(_mm512_cmp_ps_mask(x,eps,_CMP_LT_OQ)) {
                               value = _mm512_mul_ps(_1_2,x);
			   }
			   else if(_mm512_cmp_ps_mask(x,_15,_CMP_LT_OQ)) {
                               xx   = _mm512_mul_ps(x,x);
			       sump = p[0];
			       sump = _mm512_fmadd_ps(sump,xx,p[1]);
			       sump = _mm512_fmadd_ps(sump,xx,p[2]);
			       sump = _mm512_fmadd_ps(sump,xx,p[3]);
			       sump = _mm512_fmadd_ps(sump,xx,p[4]);
			       sump = _mm512_fmadd_ps(sump,xx,p[5]);
			       sump = _mm512_fmadd_ps(sump,xx,p[6]);
			       sump = _mm512_fmadd_ps(sump,xx,p[7]);
			       sump = _mm512_fmadd_ps(sump,xx,p[8]);
			       sump = _mm512_fmadd_ps(sump,xx,p[9]);
			       sump = _mm512_fmadd_ps(sump,xx,p[10]);
			       sump = _mm512_fmadd_ps(sump,xx,p[11]);
			       sump = _mm512_fmadd_ps(sump,xx,p[12]);
			       sump = _mm512_fmadd_ps(sump,xx,p[13]);
			       sump = _mm512_fmadd_ps(sump,xx,p[14]);
			       xx   = _mm512_sub_ps(xx,_225);
			       const __m512 t0 = _mm512_fmadd_ps(_mm512_add_ps(xx,q[0]),xx,q[1]);
			       const __m512 t1 = _mm512_fmadd_ps(t0,xx,q[2]);
			       const __m512 t2 = _mm512_fmadd_ps(t1,xx,q[3]);
			       const __m512 t3 = _mm512_fmadd_ps(t2,xx,q[4]);
			       sumq             = t3;
			       value            = _mm512_mul_ps(_mm512_div_ps(sump,sumq),x);
			   }
			   else if(_mm512_cmp_ps_mask(xmax,x,_CMP_LT_OQ)) {
                               value            = huge;
			   }
			   else {
                               xx               = _mm512_sub_ps(_mm512_div_ps(_1,x),rec15);
			       const __m512 t0 = _mm512_fmadd_ps(pp[0],xx,pp[1]);
			       const __m512 c0 = _mm512_fmadd_ps(_mm512_add_ps(xx,qq[0]),xx,qq[1]);
			       const __m512 t1 = _mm512_fmadd_ps(t0,xx,pp[2]);
			       const __m512 c1 = _mm512_fmadd_ps(c0,xx,qq[2]);
			       const __m512 t2 = _mm512_fmadd_ps(t1,xx,pp[3]);
			       const __m512 c2 = _mm512_fmadd_ps(c1,xx,qq[3]);
			       const __m512 t3 = _mm512_fmadd_ps(t2,xx,pp[4]);
			       const __m512 c3 = _mm512_fmadd_ps(c2,xx,qq[4]);
			       const __m512 t4 = _mm512_fmadd_ps(t3,xx,pp[5]);
			       const __m512 c4 = _mm512_fmadd_ps(c3,xx,qq[5]);
			       const __m512 t5 = _mm512_fmadd_ps(t4,xx,pp[6]);
			       const __m512 c5 = _mm512_fmadd_ps(c4,xx,qq[6]);
			       const __m512 t6 = _mm512_fmadd_ps(t5,xx,pp[7]);
			       sump             = t6;
			       sumq             = c5;
			       value            = _mm512_div_ps(sump,sumq);
			       const __mmask16 m = _mm512_cmp_ps_mask(_mm512_sub_ps(xmax,_15),_CMP_LT_OQ);
#if (USE_SLEEF_LIB) == 1
                               a                = _mm512_mask_blend_ps(m,xexpf(x),
			                                                           xexpf(_mm512_sub_ps(x,_40)));
#else
			       a                = _mm512_mask_blend_ps(m,_mm512_exp_ps(x),
			                                                           _mm512_exp_ps(_mm512_sub_ps(x,_40)));
#endif
			       b                = _mm512_mask_blend_ps(m,_1,_40);
			       const __m512 tmp= _mm512_add_ps(_mm512_mul_ps(value,a),
			                                        _mm512_mul_ps(pbar,a));
			       value            = _mm512_mul_ps(_mm512_div_ps(tmp,_mm512_sqrt_ps(x)),b);
			   }
			   if(_mm512_cmp_ps_mask(arg,_0,_CMP_LT_OQ)) {
                              value             = zmm16r4_negate(value);
			   }
			   bessel_i1            = value
			   return (bessel_i1);
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
		      __m512d  
		      normal_01_cdf_zmm8r8(const __m512d x) {
		          
		          __m512d a1 = _mm512_set1_pd(0.398942280444e+00);
		          __m512d a2 = _mm512_set1_pd(0.399903438504e+00);
		          __m512d a3 = _mm512_set1_pd(5.75885480458e+00);
                          __m512d a4 = _mm512_set1_pd(29.8213557808e+00);
                          __m512d a5 = _mm512_set1_pd(2.62433121679e+00);
                          __m512d a6 = _mm512_set1_pd(48.6959930692e+00);
                          __m512d a7 = _mm512_set1_pd(5.92885724438e+00);
                          __m512d b0 = _mm512_set1_pd(0.398942280385e+00);
                          __m512d b1 = _mm512_set1_pd(3.8052e-08);
                          __m512d b2 = _mm512_set1_pd(1.00000615302e+00);
                          __m512d b3 = _mm512_set1_pd(3.98064794e-04);
                          __m512d b4 = _mm512_set1_pd(1.98615381364e+00);
                          __m512d b5 = _mm512_set1_pd(0.151679116635e+00);
                          __m512d b6 = _mm512_set1_pd(5.29330324926e+00);
                          __m512d b7 = _mm512_set1_pd(4.8385912808e+00);
                          __m512d b8 = _mm512_set1_pd(15.1508972451e+00);
                          __m512d b9 = _mm512_set1_pd(0.742380924027e+00);
                          __m512d b10= _mm512_set1_pd(30.789933034e+00);
                          __m512d b11= _mm512_set1_pd(3.99019417011e+00);
                          __m512d C1 = _mm512_set1_pd(1.0);
                          __m512d C128 = _mm512_set1_pd(1.28);
                          __m512d C05  = _mm512_set1_pd(0.5);
                          __m512d C127 = _mm512_set1_pd(12.7);
                          __m512d absx,y,q,cdf,t0,t1;
                          __mmask8 m0,m1,m2;
                          m2   = _mm512_cmp_pd_mask(x,_mm512_setzero_pd(),_CMP_LT_OQ);
                          absx = _mm512_abs_pd(x);
                          m0   = _mm512_cmp_pd_mask(x,C128,_CMP_LE_OQ);
                          y    = _mm512_mul_pd(C05,
                                        _mm512_mul_pd(x,x));
                          m1   = _mm512_cmp_pd_mask(x,C127,_CMP_LE_OQ);
                          if(m0) {
                             register __m512d ya3;
                             register __m512d ya5a6
                             register __m512d ya7;
                             register __m512d a2y;
                             ya7   = _mm512_add_pd(y,a7);
                             ya5a6 = _mm512_add_pd(y,_mm512_add_pd(a5,a6));
                             a2y   = _mm512_mul_pd(a2,y);
                             ya3a4 = _mm512_sub_pd(_mm512_add_pd(y,a3),a4);
                             q     = _mm512_sub_pd(a1,
                                           _mm512_div_pd(a2y,
                                                  _mm512_div_pd(ya3a4,
                                                        _mm512_div_pd(ya5a6,ya7))));
                          }
                          else if(m1) {
                             register __m512d expmy;
                             register __m512d absb1;
                             register __m512d absb3;
                             register __m512d absb5;
                             register __m512d absb7;
                             register __m512d absb9;
                             register __m512d absb11;
#if (USE_SLEEF_LIB) == 1                          
                             expmy = _mm512_mul_pd(xexp(negate_zmm8r8(y)),b0);
#else
                             expmy = _mm512_mul_pd(_mm512_exp_pd(negate_zmm8r8(y)),b0); 
#endif 
                             absb1 = _mm512_sub_pd(absx,b1);
                             absb3 = _mm512_add_pd(absx,b3);
                             absb5 = _mm512_sub_pd(absx,b5);
                             absb7 = _mm512_add_pd(absx,b7);
                             absb9 = _mm512_add_pd(absx,b9);
                             absb11= _mm512_add_pd(absx,b11);
                             t0    = (absb1+b2/(absb3+b4/(absb5+b6/(absb7-b8/(absb9+b10/(absb11))))));
                             q     = _mm512_div_pd(expmy,t0);
                          }
                          else {
                             q = _mm512_setzero_pd();
                          }
                          
                          cdf = _mm512_mask_blend_pd(m2,_mm512_sub_pd(C1,q),q);
                          return (cdf);
		    }
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512  
		      normal_01_cdf_zmm16r4(const __m512 x) {
		          
		          __m512 a1 = _mm512_set1_ps(0.398942280444f);
		          __m512 a2 = _mm512_set1_ps(0.399903438504f);
		          __m512 a3 = _mm512_set1_ps(5.75885480458f);
                          __m512 a4 = _mm512_set1_ps(29.8213557808f);
                          __m512 a5 = _mm512_set1_ps(2.62433121679f);
                          __m512 a6 = _mm512_set1_ps(48.6959930692f);
                          __m512 a7 = _mm512_set1_ps(5.92885724438f);
                          __m512 b0 = _mm512_set1_ps(0.398942280385f);
                          __m512 b1 = _mm512_set1_ps(3.8052e-08f);
                          __m512 b2 = _mm512_set1_ps(1.00000615302f);
                          __m512 b3 = _mm512_set1_ps(3.98064794e-04f);
                          __m512 b4 = _mm512_set1_ps(1.98615381364f);
                          __m512 b5 = _mm512_set1_ps(0.151679116635f);
                          __m512 b6 = _mm512_set1_ps(5.29330324926f);
                          __m512 b7 = _mm512_set1_ps(4.8385912808f);
                          __m512 b8 = _mm512_set1_ps(15.1508972451f);
                          __m512 b9 = _mm512_set1_ps(0.742380924027f);
                          __m512 b10= _mm512_set1_ps(30.789933034f);
                          __m512 b11= _mm512_set1_ps(3.99019417011f);
                          __m512 C1 = _mm512_set1_ps(1.0);
                          __m512 C128 = _mm512_set1_ps(1.28f);
                          __m512 C05  = _mm512_set1_ps(0.5f);
                          __m512 C127 = _mm512_set1_ps(12.7f);
                          __m512 absx,y,q,cdf,t0;
                          __mmask16 m0,m1,m2;
                          m2   = _mm512_cmp_ps_mask(x,_mm512_setzero_pd(),_CMP_LT_OQ);
                          absx = _mm512_abs_ps(x);
                          m0   = _mm512_cmp_ps_mask(x,C128,_CMP_LE_OQ);
                          y    = _mm512_mul_ps(C05,
                                        _mm512_mul_ps(x,x));
                          m1   = _mm512_cmp_ps_mask(x,C127,_CMP_LE_OQ);
                          if(m0) {
                             register __m512 ya3;
                             register __m512 ya5a6
                             register __m512 ya7;
                             register __m512 a2y;
                             ya7   = _mm512_add_ps(y,a7);
                             ya5a6 = _mm512_add_ps(y,_mm512_add_ps(a5,a6));
                             a2y   = _mm512_mul_ps(a2,y);
                             ya3a4 = _mm512_sub_ps(_mm512_add_ps(y,a3),a4);
                             q     = _mm512_sub_ps(a1,
                                           _mm512_div_ps(a2y,
                                                  _mm512_div_ps(ya3a4,
                                                        _mm512_div_ps(ya5a6,ya7))));
                          }
                          else if(m1) {
                             register __m512 expmy;
                             register __m512 absb1;
                             register __m512 absb3;
                             register __m512 absb5;
                             register __m512 absb7;
                             register __m512 absb9;
                             register __m512 absb11;
#if (USE_SLEEF_LIB) == 1                          
                             expmy = _mm512_mul_ps(xexpf(negate_zmm16r4(y)),b0);
#else
                             expmy = _mm512_mul_ps(_mm512_exp_ps(negate_zmm16r4(y)),b0); 
#endif 
                             absb1 = _mm512_sub_ps(absx,b1);
                             absb3 = _mm512_add_ps(absx,b3);
                             absb5 = _mm512_sub_ps(absx,b5);
                             absb7 = _mm512_add_ps(absx,b7);
                             absb9 = _mm512_add_ps(absx,b9);
                             absb11= _mm512_add_ps(absx,b11);
                             t0    = (absb1+b2/(absb3+b4/(absb5+b6/(absb7-b8/(absb9+b10/(absb11))))));
                             q     = _mm512_div_ps(expmy,t0);
                          }
                          else {
                             q = _mm512_setzero_ps();
                          }
                          
                          cdf = _mm512_mask_blend_ps(m2,_mm512_sub_ps(C1,q),q);
                          return (cdf);
		    }
		    
		    
		    
/*
          !*****************************************************************************80
!
!! NORMAL_01_SAMPLE samples the standard normal probability distribution.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but
!    generates two values at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the standard normal PDF.
!
*/
	
	
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d   	    
		      normal_01_sample_zmm8r8(__m512i & seed) {
		           
		           __m512d CN20 = _mm512_set1_pd(-2.0);
		           __m512d C628318530717958647692529  = 
		                                    _mm512_set1_pd(6.28318530717958647692529);
		           __m512d C00  = _mm512_setzero_pd();
		           register __m512d r1,r2,s,c,t0,t1,arg;
		           register __m512d x;
		           static __m512d y = _mm512_setzero_pd();
		           static int32_t used = -1;
		           if(used == -1) used = 0;
		           
		           if((used%2)==0) {
		           
		                 while(true) {
		                    r1 = uniform_01_zmm8r8(seed);
		                    if(_mm512_cmp_pd_mask(r1,C00,_CMP_NE_OQ)) break;
		               }
		               r2 = uniform_01_zmm8r8(seed);
		               arg= _mm512_mul_pd(r2,C628318530717958647692529)
#if (USE_SLEEF_LIB) == 1       
                               t0 = _mm512_mul_pd(CN20,xlog(r1));
                               c  = xcos(arg);
                               s  = xsin(arg);	
#else
                               t0 = _mm512_mul_pd(CN20,_mm512_log_pd(r1));
                               c  = _mm512_cos_pd(arg);
                               s  = _mm512_sin_pd(arg);
#endif                           	               
                               x  = _mm512_mul_pd(_mm512_sqrt_pd(t0),c);
                               y  = _mm512_mul_pd(_mm512_sqrt_pd(t0),s);
		           }
		           else {
		               x = y;
		           }
		           used += 1;
		           return (x);
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512  	    
		      normal_01_sample_zmm16r4(__m512i & seed) {
		          
		          return (_mm512_castpd_ps(normal_01_sample_zmm8r8(seed)));
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
		      __m512d 
		      reciprocal_cdf_zmm8r8(const __m512d x,
		                            const __m512d a,
		                            const __m512d b) {
		          
		          register __m512d ax,ab,l1,l2;
		          register __m512d cdf;       
		          ax = _mm512_div_pd(a,x);
#if (USE_SLEEF_LIB) == 1   
                          l1 = xlog(ax);
#else
                          l1 = _mm512_log_pd(ax);
#endif		                         
                          ab = _mm512_div_pd(a,b);
#if (USE_SLEEF_LIB) == 1 
                          l2 = xlog(ab);
#else
                          l2 = _mm512_log_pd(ab);
#endif                          
                          cdf= _mm512_div_pd(l1,l2);
                          return (cdf);
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 
		      reciprocal_cdf_zmm16r4(const __m512 x,
		                             const __m512 a,
		                             const __m512 b) {
		          
		          register __m512 ax,ab,l1,l2;
		          register __m512 cdf;       
		          ax = _mm512_div_ps(a,x);
#if (USE_SLEEF_LIB) == 1   
                          l1 = xlogf(ax);
#else
                          l1 = _mm512_log_ps(ax);
#endif		                         
                          ab = _mm512_div_ps(a,b);
#if (USE_SLEEF_LIB) == 1 
                          l2 = xlogf(ab);
#else
                          l2 = _mm512_log_ps(ab);
#endif                          
                          cdf= _mm512_div_ps(l1,l2);
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
		      __m512d 		     
		      reciprocal_cdf_inv_zmm8r8(const __m512d cdf,
		                                const __m512d a,
		                                const __m512d b) {
		         
		           register __m512d C1 = _mm512_set1_pd(1.0);
		           register __m512d pow1,pow2,cdf1;
		           register __m512d inv;
		           cdf1 = _mm512_sub_pd(cdf,C1);
		           pow2 = _mm512_pow_pd(b,cdf);
		           pow1 = _mm512_pow_pd(a,cdf1);
		           inv  = _mm512_div_pd(pow2,pow1);
		           return (inv);                          
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512		     
		      reciprocal_cdf_inv_zmm16r4(const __m512 cdf,
		                                const __m512 a,
		                                const __m512 b) {
		         
		           register __m512 C1 = _mm512_set1_ps(1.0f);
		           register __m512 pow1,pow2,cdf1;
		           register __m512 inv;
		           cdf1 = _mm512_sub_ps(cdf,C1);
		           pow2 = _mm512_pow_ps(b,cdf);
		           pow1 = _mm512_pow_ps(a,cdf1);
		           inv  = _mm512_div_ps(pow2,pow1);
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
		      __m512d 
		      reciprocal_mean_zmm8r8(const __m512d a,
		                             const __m512d b) {
		           
		           register __m512d ab,amb,l1;
		           register __m512d mean;
		           amb = _mm512_sub_pd(a,b);
		           ab  = _mm512_div_pd(a,b);
#if (USE_SLEEF_LIB) == 1  
                           l1  = xlog(ab);
#else
                           l1  = _mm512_log_pd(ab);
#endif		                                  
                           mean= _mm512_div_pd(amb,l1);
                           return (mean);
		     }	 
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 
		      reciprocal_mean_zmm16r4(const __m512 a,
		                             const __m512 b) {
		           
		           register __m512 ab,amb,l1;
		           register __m512 mean;
		           amb = _mm512_sub_ps(a,b);
		           ab  = _mm512_div_ps(a,b);
#if (USE_SLEEF_LIB) == 1  
                           l1  = xlogf(ab);
#else
                           l1  = _mm512_log_ps(ab);
#endif		                                  
                           mean= _mm512_div_ps(amb,l1);
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
		      __m512d
		      reciprocal_pdf_zmm8r8(    const __m512d x,
		                                const __m512d a,
		                                const __m512d b) {
		          
		          register __m512d C1 = _mm512_set1_pd(1.0);
		          register __m512d ba,l1;
		          register __m512d pdf;
		          ba = _mm512_div_pd(b,a);
#if (USE_SLEEF_LIB) == 1  
                          l1 = _mm512_mul_pd(x,xlog(ba));
#else
                          l1 = _mm512_mul_pd(x,_mm512_log_pd(ba));
#endif		         
                          pdf= _mm512_div_pd(C1,l1);
                          return (pdf);                            
		    }
		    
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      reciprocal_pdf_zmm16r4(    const __m512 x,
		                                const __m512 a,
		                                const __m512 b) {
		          
		          register __m512 C1 = _mm512_set1_ps(1.0f);
		          register __m512 ba,l1;
		          register __m512 pdf;
		          ba = _mm512_div_ps(b,a);
#if (USE_SLEEF_LIB) == 1  
                          l1 = _mm512_mul_ps(x,xlogf(ba));
#else
                          l1 = _mm512_mul_ps(x,_mm512_log_ps(ba));
#endif		         
                          pdf= _mm512_div_ps(C1,l1);
                          return (pdf);                            
		    }
		    
		    
/*
         !*****************************************************************************80
!
!! RECIPROCAL_SAMPLE samples the Reciprocal PDF.
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
		      __m512d
                      reciprocal_sample_zmm8r8( __m512i & seed,
                                               const __m512d a,
                                               const __m512d b) {
                           
                           register __m512d C1 = _mm512_set1_pd(1.0);
                           register __m512d pow1,pow2,arg,cdf;
                           register __m512d sample;
                           cdf = uniform_01_zmm8r8(seed);
                           arg = _mm512_sub_pd(cdf,C1);
                           pow1= _mm512_pow_pd(b,cdf);
                           pow2= _mm512_pow_pd(a,arg);
                           sample = _mm512_div_pd(pow1,pow2);
                           return (sample);                          
                    }
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
                      reciprocal_sample_zmm16r4( __m512i & seed,
                                               const __m512 a,
                                               const __m512 b) {
                         
                         return (_mm512_castpd_ps(reciprocal_sample_zmm8r8(seed,a,b)));                          
                    }
                    
                    
/*
     !*****************************************************************************80
!
!! RECIPROCAL_VARIANCE returns the variance of the Reciprocal PDF.
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
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!       
*/

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d  
                      reciprocal_var_zmm8r8(const __m512d a,
                                            const __m512d b) {
                        
                           register __m512d C2 = _mm512_set1_pd(2.0);
                           register __m512d ab,amb,dd,dm2,dp2,t0;
                           register var;
                           ab  = _mm512_div_pd(a,b);
                           amb = _mm512_sub_pd(a,b);
#if (USE_SLEEF_LIB) == 1  
                           d   = xlog(ab);
#else
                           d   = _mm512_log_pd(ab);
#endif                                             
                           dd  = _mm512_mul_pd(C2,_mm512_mul_pd(d,d));
                           dm2 = _mm512_mul_pd(a,_mm512_sub_pd(d,C2));
                           dp2 = _mm512_mul_pd(b,_mm512_add_pd(d,C2));
                           t0  = _mm512_fmadd_pd(amb,dm2,dp2);
                           var = _mm512_div_pd(t0,dd);
                           return (var);    
                     }
                     
                     
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512  
                      reciprocal_var_zmm16r4(const __m512 a,
                                            const __m512 b) {
                        
                           register __m512 C2 = _mm512_set1_ps(2.0f);
                           register __m512 ab,amb,dd,dm2,dp2,t0;
                           register var;
                           ab  = _mm512_div_ps(a,b);
                           amb = _mm512_sub_ps(a,b);
#if (USE_SLEEF_LIB) == 1  
                           d   = xlogf(ab);
#else
                           d   = _mm512_log_ps(ab);
#endif                                             
                           dd  = _mm512_mul_ps(C2,_mm512_mul_ps(d,d));
                           dm2 = _mm512_mul_ps(a,_mm512_sub_ps(d,C2));
                           dp2 = _mm512_mul_ps(b,_mm512_add_ps(d,C2));
                           t0  = _mm512_fmadd_ps(amb,dm2,dp2);
                           var = _mm512_div_ps(t0,dd);
                           return (var);    
                     }
                     
                     
/*
      !*****************************************************************************80
!
!! SECH returns the hyperbolic secant.
!
!  Discussion:
!
!    SECH ( X ) = 1.0D+00 / COSH ( X ) = 2.0D+00 / ( EXP ( X ) + EXP ( - X ) )
!
!    SECH is not a built-in function in FORTRAN, and occasionally it
!    is handier, or more concise, to be able to refer to it directly
!    rather than through its definition in terms of the sine function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) SECH, the hyperbolic secant of X.
!   
*/

          
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
                      sech_zmm8r8(const __m512d x) {
                          
                          register __m512d C1 = _mm512_set1_pd(1.0);
                          register __m512d csh;
                          register __m512d sch;
#if (USE_SLEEF_LIB) == 1  
                          csh = xcosh(x);
#else
                          csh = _mm512_cosh_pd(x);
#endif                          
                          sch = _mm512_div_pd(C1,csh);
                          return (sch);
                      }
                      
                      
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
                      sech_zmm16r4(const __m512 x) {
                          
                          register __m512 C1 = _mm512_set1_ps(1.0f);
                          register __m512 csh;
                          register __m512 sch;
#if (USE_SLEEF_LIB) == 1  
                          csh = xcoshf(x);
#else
                          csh = _mm512_cosh_ps(x);
#endif                          
                          sch = _mm512_div_ps(C1,csh);
                          return (sch);
                      }
                      
                      
/*
        !*****************************************************************************80
!
!! SECH_CDF evaluates the Hyperbolic Secant CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, the parameter of the PDF.
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
		      __m512d                     
                      sech_cdf_zmm8r8(const __m512d x,
                                      const __m512d a,
                                      const __m512d b) {
                         
                          register __m512d C2 = _mm512_set1_pd(2.0);
                          register __m512d C031830988618379067153777 = 
                                                _mm512_set1_pd(0.31830988618379067153777);
                          register __m512d y,expy,atn;
                          register __m512d cdf;
                          y = _mm512_div_pd(_mm512_sub_pd(x,a),b);
#if (USE_SLEEF_LIB) == 1 
                          atn = xatan(xexp(y));
#else
                          atn = _mm512_atan_pd(_mm512_exp_pd(y));
#endif                                           
                          cdf = _mm512_mul_pd(C2,_mm512_mul_pd(atn,
                                                 C031830988618379067153777));
                          return (cdf);
                     }
                     
                     
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512                     
                      sech_cdf_zmm16r4(const __m512 x,
                                      const __m512 a,
                                      const __m512 b) {
                         
                          register __m512 C2 = _mm512_set1_ps(2.0f);
                          register __m512 C031830988618379067153777 = 
                                                _mm512_set1_pd(0.31830988618379067153777f);
                          register __m512 y,expy,atn;
                          register __m512 cdf;
                          y = _mm512_div_ps(_mm512_sub_ps(x,a),b);
#if (USE_SLEEF_LIB) == 1 
                          atn = xatanf(xexpf(y));
#else
                          atn = _mm512_atan_ps(_mm512_exp_ps(y));
#endif                                           
                          cdf = _mm512_mul_ps(C2,_mm512_mul_ps(atn,
                                                 C031830988618379067153777));
                          return (cdf);
                     }
                     
                     
/*
         !*****************************************************************************80
!
!! SECH_CDF_INV inverts the Hyperbolic Secant CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
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
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) X, the corresponding argument of the CDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d  
		      sech_cdf_inv_zmm8r8(const __m512d cdf,
		                          const __m512d a,
		                          const __m512d b) {
		          
		           register __m512d C157079632679489661923132 = 
		                                 _mm512_set1_pd(1.57079632679489661923132);
		           register __m512d targ,ab,tan,log;
		           register __m512d x;
		           ab   = _mm512_add_pd(a,b);
		           targ = _mm512_mul_pd(cdf,C157079632679489661923132);
#if (USE_SLEEF_LIB) == 1 
                           tan  = xtan(targ);
                           log  = xlog(tan);
#else
                           tan  = _mm512_tan_pd(targ);
                           log  = _mm512_log_pd(tan);
#endif		                           
                           x    = _mm512_mul_pd(log,ab);
                           return (x);
		     }     
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512  
		      sech_cdf_inv_zmm16r4(const __m512 cdf,
		                          const __m512 a,
		                          const __m512 b) {
		          
		           register __m512 C157079632679489661923132 = 
		                                 _mm512_set1_pd(1.57079632679489661923132f);
		           register __m512 targ,ab,tan,log;
		           register __m512 x;
		           ab   = _mm512_add_ps(a,b);
		           targ = _mm512_mul_ps(cdf,C157079632679489661923132);
#if (USE_SLEEF_LIB) == 1 
                           tan  = xtanf(targ);
                           log  = xlogf(tan);
#else
                           tan  = _mm512_tan_ps(targ);
                           log  = _mm512_log_ps(tan);
#endif		                           
                           x    = _mm512_mul_ps(log,ab);
                           return (x);
		     }    
		     
		     
/*
       !*****************************************************************************80
!
!! SECH_PDF evaluates the Hypebolic Secant PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = sech ( ( X - A ) / B ) / ( PI * B )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
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
		      __m512d  
                      sech_pdf_zmm8r8(const __m512d x
                                      const __m512d a,
                                      const __m512d b) {
                          
                          
                          register __m512d C314159265358979323846264 = 
                                                           _mm512_set1_pd(3.14159265358979323846264);
                          register __m512d y,pib,sech;
                          pib = _mm512_mul_pd(b,C314159265358979323846264);
                          y   = _mm512_div_pd(_mm512_sub_pd(x,a),b);
                          sech= sech_zmm8r8(y);
                          pdf = _mm512_div_pd(sech,pib);
                          return (pdf);                
                    }
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512  
                      sech_pdf_zmm16r4(const __m512 x
                                      const __m512 a,
                                      const __m512 b) {
                          
                          
                          register __m512 C314159265358979323846264 = 
                                                           _mm512_set1_pd(3.14159265358979323846264f);
                          register __m512 y,pib,sech;
                          pib = _mm512_mul_ps(b,C314159265358979323846264);
                          y   = _mm512_div_ps(_mm512_sub_ps(x,a),b);
                          sech= sech_zmm16r4(y);
                          pdf = _mm512_div_ps(sech,pib);
                          return (pdf);                
                    }
                    
                    
/*
   !*****************************************************************************80
!
!! SECH_SAMPLE samples the Hyperbolic Secant PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
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
		      __m512d 
                      sech_sample_zmm8r8(const __m512d a,
                                         const __m512d b,
                                         const __m512i & seed) {
                        
                          register __m512d C157079632679489661923132 = 
		                                 _mm512_set1_pd(1.57079632679489661923132);
                          register __m512d cdf,ab,targ,tan,log;
                          register __m512d sample;
                          ab   = _mm512_add_pd(a,b);
                          cdf  = uniform_01_zmm8r8(seed);
                          targ = _mm512_mul_pd(cdf,C157079632679489661923132);
                                   
#if (USE_SLEEF_LIB) == 1 
                          tan = xtan(targ);
                          log = xlog(tan);
#else
                          tan = _mm512_tan_pd(targ);
                          log = _mm512_log_pd(tan);
#endif
                          sample = _mm512_mul_pd(ab,log);
                          return (sample);                                       
                     } 
                     
                     
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
                      sech_sample_zmm16r4(const __m512d a,
                                         const __m512d b,
                                         const __m512i & seed) {
                         
                          register __m512 sample;
                          sample = _mm512_castpd_ps(sech_sample_zmm8r8(a,b,seed));
                          return (sample);                   
                     }
                     
                     
/*
              !*****************************************************************************80
!
!! SECH_VARIANCE returns the variance of the Hyperbolic Secant PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the PDF.
!
*/

             
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d    
                      sech_variance_zmm8r8(const __m512d a,
                                           const __m512d b) {
                          
                            register __m512d C314159265358979323846264 = 
                                                           _mm512_set1_pd(3.14159265358979323846264);
                            register __m512d C025 = _mm512_set1_pd(0.25);
                            register __m512d pib,pow;
                            register __m512d var;
                            pib = _mm512_mul_pd(b,C314159265358979323846264);
                            pow = _mm512_mul_pd(pib,pib);
                            var = _mm512_mul_pd(C025,pow);
                            return (var);                   
                    }
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512   
                      sech_variance_zmm16r4(const __m512 a,
                                           const __m512 b) {
                          
                            register __m512 C314159265358979323846264 = 
                                                           _mm512_set1_pd(3.14159265358979323846264f);
                            register __m512 C025 = _mm512_set1_pd(0.25f);
                            register __m512 pib,pow;
                            register __m512 var;
                            pib = _mm512_mul_ps(b,C314159265358979323846264);
                            pow = _mm512_mul_ps(pib,pib);
                            var = _mm512_mul_ps(C025,pow);
                            return (var);                   
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
		      __m512d  
		      vpoly_eval_zmm8r8(const int32_t n,
		                        const __m512d * __restrict __ATTR_ALIGN__(64) a,
		                        const __m512d x) {
		         
		         register __m512d vpoly;
		         vpoly = _mm512_load_pd(&a[n]);
		         for(int32_t i=n; i != 0; --i) {
		             register __m512d t0 = a[i];
		             vpoly = _mm512_fmadd_pd(vpoly,x,t0);   
		         }  
		         return (vpoly);              
		    }
		    
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512  
		      vpoly_eval_zmm16r4(const int32_t n,
		                        const __m512 * __restrict __ATTR_ALIGN__(64) a,
		                        const __m512 x) {
		         
		         register __m512 vpoly;
		         vpoly = _mm512_load_ps(&a[n]);
		         for(int32_t i=n; i != 0; --i) {
		             register __m512 t0 = a[i];
		             vpoly = _mm512_fmadd_ps(vpoly,x,t0);   
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
		      __m512d    
		      normal_01_cdf_inv_zmm8r8(const __m512d p) {
		            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(64) static __m512d  a[8] = {
		                     _mm512_set1_pd(3.3871328727963666080e+00),
                                     _mm512_set1_pd(1.3314166789178437745e+02),
                                     _mm512_set1_pd(1.9715909503065514427e+03),
                                     _mm512_set1_pd(1.3731693765509461125e+04),
                                     _mm512_set1_pd(4.5921953931549871457e+04),
                                     _mm512_set1_pd(6.7265770927008700853e+04),
                                     _mm512_set1_pd(3.3430575583588128105e+04),
                                     _mm512_set1_pd(2.5090809287301226727e+03)};   
                            __attribute__((section(".rodata")))  
		            __ATTR_ALIGN__(64) static __m512d   b[8] = {
		                      _mm512_set1_pd(1.0e+00),
                                      _mm512_set1_pd(4.2313330701600911252e+01),
                                      _mm512_set1_pd(6.8718700749205790830e+02),
                                      _mm512_set1_pd(5.3941960214247511077e+03),
                                      _mm512_set1_pd(2.1213794301586595867e+04),
                                      _mm512_set1_pd(3.9307895800092710610e+04),
                                      _mm512_set1_pd(2.8729085735721942674e+04),
                                      _mm512_set1_pd(5.2264952788528545610e+03)}; 
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(64) static __m512d   c[8] = {
		                      _mm512_set1_pd(1.42343711074968357734e+00),
                                      _mm512_set1_pd(4.63033784615654529590e+00),
                                      _mm512_set1_pd(5.76949722146069140550e+00),
                                      _mm512_set1_pd(3.64784832476320460504e+00),
                                      _mm512_set1_pd(1.27045825245236838258e+00),
                                      _mm512_set1_pd(2.41780725177450611770e-01),
                                      _mm512_set1_pd(2.27238449892691845833e-02),
                                      _mm512_set1_pd(7.74545014278341407640e-04)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512d   d[8] = {
                                      _mm512_set1_pd(1.0e+00),
                                      _mm512_set1_pd(2.05319162663775882187e+00),
                                      _mm512_set1_pd(1.67638483018380384940e+00),
                                      _mm512_set1_pd(6.89767334985100004550e-01),
                                      _mm512_set1_pd(1.48103976427480074590e-01),
                                      _mm512_set1_pd(1.51986665636164571966e-02),
                                      _mm512_set1_pd(5.47593808499534494600e-04),
                                      _mm512_set1_pd(1.05075007164441684324e-09)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512d   e[8] = {
                                      _mm512_set1_pd(6.65790464350110377720e+00),
                                      _mm512_set1_pd(5.46378491116411436990e+00),
                                      _mm512_set1_pd(1.78482653991729133580e+00),
                                      _mm512_set1_pd(2.96560571828504891230e-01),
                                      _mm512_set1_pd(2.65321895265761230930e-02),
                                      _mm512_set1_pd(1.24266094738807843860e-03),
                                      _mm512_set1_pd(2.71155556874348757815e-05),
                                      _mm512_set1_pd(2.01033439929228813265e-07)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512d   f[8] = {
                                      _mm512_set1_pd(1.0e+00),
                                      _mm512_set1_pd(5.99832206555887937690e-01),
                                      _mm512_set1_pd(1.36929880922735805310e-01),
                                      _mm512_set1_pd(1.48753612908506148525e-02),
                                      _mm512_set1_pd(7.86869131145613259100e-04), 
                                      _mm512_set1_pd(1.84631831751005468180e-05),
                                      _mm512_set1_pd(1.42151175831644588870e-07),
                                      _mm512_set1_pd(2.04426310338993978564e-15)};
                          __m512d const1 = _mm512_set1_pd(0.180625e+00);
                          __m512d const2 = _mm512_set1_pd(1.6e+00);
                          __m512d split1 = _mm512_set1_pd(0.425e+00);
                          __m512d split2 = _mm512_set1_pd(5.0e+00);
                          __m512d C0     = _mm512_setzero_pd();
                          __m512d C1     = _mm512_set1_pd(1.0);
                          __m512d C05    = _mm512_set1_pd(0.5);
                          register __m512d q,r,t0,t1;
                          register __m512d x;
                          q = _mm512_sub_pd(p,C05);
                          if(_mm512_cmp_pd_mask(q,split1,_CMP_LE_OQ)) {
                             r = _mm512_sub_pd(const1,_mm512_mul_pd(q,q));
                             t0= vpoly_eval_zmm8r8(8,a,r);
                             t1= vpoly_eval_zmm8r8(8,b,r);
                             x = _mm512_div_pd(_mm512_mul_pd(q,t0),t1);
                          } 
                          else {
                             const __mmask8 m = _mm512_cmp_pd_mask(q,C0,_CMP_LT_OQ);
                             r                = _mm512_mask_blend_pd(m,_mm512_sub_pd(C1,p),p);
                             if(_mm512_cmp_pd_mask(r,C0,_CMP_LE_OQ)) {
                                x = _mm512_set1_pd(std::numeric_limits<double>::max());
                             }
                             else {
#if (USE_SLEEF_LIB) == 1     
                                r = _mm512_sqrt_pd(negate_zmm8r8(xlog(r)));
#else
                                r = _mm512_sqrt_pd(negate_zmm8r8(_mm512_log_pd(r)));
#endif                         
                                const __mmask8 m = _mm512_cmp_pd_mask(r,split2,_CMP_LE_OQ);
                                r                = _mm512_mask_blend_pd(m,_mm512_sub_pd(r,split2),
                                                                          _mm512_sub_pd(r,const2));
                                t0               = _mm512_div_pd(vpoly_eval_zmm8r8(8,c,r),
                                                                 vpoly_eval_zmm8r8(8,d,r));
                                t1               = _mm512_div_pd(vpoly_eval_zmm8r8(8,e,r),
                                                                 vpoly_eval_zmm8r8(8,f,r));
                                x                = _mm512_mask_blend_pd(m,t1,t0);      
                             }
                             if(_mm512_cmp_pd_mask(q,C0,_CMP_LT_OQ)) x = negate_zmm8r8(x);
                          }
                          return (x);
                          
		    }
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512    
		      normal_01_cdf_inv_zmm16r4(const __m512 p) {
		            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(64) static __m512  a[8] = {
		                     _mm512_set1_ps(3.3871328727963666080e+00f),
                                     _mm512_set1_ps(1.3314166789178437745e+02f),
                                     _mm512_set1_ps(1.9715909503065514427e+03f),
                                     _mm512_set1_ps(1.3731693765509461125e+04f),
                                     _mm512_set1_ps(4.5921953931549871457e+04f),
                                     _mm512_set1_ps(6.7265770927008700853e+04f),
                                     _mm512_set1_ps(3.3430575583588128105e+04f),
                                     _mm512_set1_ps(2.5090809287301226727e+03f)};   
                            __attribute__((section(".rodata")))  
		            __ATTR_ALIGN__(64) static __m512   b[8] = {
		                      _mm512_set1_ps(1.0e+00),
                                      _mm512_set1_ps(4.2313330701600911252e+01f),
                                      _mm512_set1_ps(6.8718700749205790830e+02f),
                                      _mm512_set1_ps(5.3941960214247511077e+03f),
                                      _mm512_set1_ps(2.1213794301586595867e+04f),
                                      _mm512_set1_ps(3.9307895800092710610e+04f),
                                      _mm512_set1_ps(2.8729085735721942674e+04f),
                                      _mm512_set1_ps(5.2264952788528545610e+03f)}; 
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(64) static __m512   c[8] = {
		                      _mm512_set1_ps(1.42343711074968357734e+00f),
                                      _mm512_set1_ps(4.63033784615654529590e+00f),
                                      _mm512_set1_ps(5.76949722146069140550e+00f),
                                      _mm512_set1_ps(3.64784832476320460504e+00f),
                                      _mm512_set1_ps(1.27045825245236838258e+00f),
                                      _mm512_set1_ps(2.41780725177450611770e-01f),
                                      _mm512_set1_ps(2.27238449892691845833e-02f),
                                      _mm512_set1_ps(7.74545014278341407640e-04f)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512   d[8] = {
                                      _mm512_set1_ps(1.0e+00),
                                      _mm512_set1_ps(2.05319162663775882187e+00f),
                                      _mm512_set1_ps(1.67638483018380384940e+00f),
                                      _mm512_set1_ps(6.89767334985100004550e-01f),
                                      _mm512_set1_ps(1.48103976427480074590e-01f),
                                      _mm512_set1_ps(1.51986665636164571966e-02f),
                                      _mm512_set1_ps(5.47593808499534494600e-04f),
                                      _mm512_set1_ps(1.05075007164441684324e-09f)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512   e[8] = {
                                      _mm512_set1_ps(6.65790464350110377720e+00f),
                                      _mm512_set1_ps(5.46378491116411436990e+00f),
                                      _mm512_set1_ps(1.78482653991729133580e+00f),
                                      _mm512_set1_ps(2.96560571828504891230e-01f),
                                      _mm512_set1_ps(2.65321895265761230930e-02f),
                                      _mm512_set1_ps(1.24266094738807843860e-03f),
                                      _mm512_set1_ps(2.71155556874348757815e-05f),
                                      _mm512_set1_ps(2.01033439929228813265e-07f)};
                           __attribute__((section(".rodata")))
                           __ATTR_ALIGN__(64) static __m512   f[8] = {
                                      _mm512_set1_ps(1.0e+00),
                                      _mm512_set1_ps(5.99832206555887937690e-01f),
                                      _mm512_set1_ps(1.36929880922735805310e-01f),
                                      _mm512_set1_ps(1.48753612908506148525e-02f),
                                      _mm512_set1_ps(7.86869131145613259100e-04f), 
                                      _mm512_set1_ps(1.84631831751005468180e-05f),
                                      _mm512_set1_ps(1.42151175831644588870e-07f),
                                      _mm512_set1_ps(2.04426310338993978564e-15f)};
                          __m512 const1 = _mm512_set1_ps(0.180625e+00f);
                          __m512 const2 = _mm512_set1_ps(1.6e+00f);
                          __m512 split1 = _mm512_set1_ps(0.425e+00f);
                          __m512 split2 = _mm512_set1_ps(5.0e+00f);
                          __m512 C0     = _mm512_setzero_ps();
                          __m512 C1     = _mm512_set1_ps(1.0f);
                          __m512 C05    = _mm512_set1_ps(0.5f);
                          register __m512 q,r,t0,t1;
                          register __m512 x;
                          q = _mm512_sub_ps(p,C05);
                          if(_mm512_cmp_ps_mask(q,split1,_CMP_LE_OQ)) {
                             r = _mm512_sub_ps(const1,_mm512_mul_ps(q,q));
                             t0= vpoly_eval_zmm16r4(8,a,r);
                             t1= vpoly_eval_zmm16r4(8,b,r);
                             x = _mm512_div_ps(_mm512_mul_ps(q,t0),t1);
                          } 
                          else {
                             const __mmask16 m = _mm512_cmp_ps_mask(q,C0,_CMP_LT_OQ);
                             r                = _mm512_mask_blend_ps(m,_mm512_sub_ps(C1,p),p);
                             if(_mm512_cmp_ps_mask(r,C0,_CMP_LE_OQ)) {
                                x = _mm512_set1_pd(std::numeric_limits<float>::max());
                             }
                             else {
#if (USE_SLEEF_LIB) == 1     
                                r = _mm512_sqrt_ps(negate_zmm16r4(xlogf(r)));
#else
                                r = _mm512_sqrt_ps(negate_zmm16r4(_mm512_log_ps(r)));
#endif                         
                                const __mmask16 m = _mm512_cmp_ps_mask(r,split2,_CMP_LE_OQ);
                                r                = _mm512_mask_blend_ps(m,_mm512_sub_ps(r,split2),
                                                                          _mm512_sub_ps(r,const2));
                                t0               = _mm512_div_ps(vpoly_eval_zmm16r4(8,c,r),
                                                                 vpoly_eval_zmm16r4(8,d,r));
                                t1               = _mm512_div_ps(vpoly_eval_zmm16r4(8,e,r),
                                                                 vpoly_eval_zmm16r4(8,f,r));
                                x                = _mm512_mask_blend_ps(m,t1,t0);      
                             }
                             if(_mm512_cmp_ps_mask(q,C0,_CMP_LT_OQ)) x = negate_zmm16r4(x);
                          }
                          return (x);
                          
		    }
		    
		    
/*
           !*****************************************************************************80
!
!! NORMAL_CDF_INV inverts the Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
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
		      __m512d 		   
		      normal_cdf_inv_zmm8r8(const __m512d cdf,
		                            const __m512d a,
		                            const __m512d b) {
		          
		          register __m512d x2;
		          register __m512d x;
		          x2 = normal_01_cdf_inv_zmm8r8(cdf);
		          x  = _mm512_add_pd(a,_mm512_mul_pd(b,x2));
		          return (x);      
		   }
		   
		   
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 		   
		      normal_cdf_inv_zmm16r4(const __m512 cdf,
		                            const __m512 a,
		                            const __m512 b) {
		          
		          register __m512 x2;
		          register __m512 x;
		          x2 = normal_01_cdf_inv_zmm16r4(cdf);
		          x  = _mm512_add_ps(a,_mm512_mul_ps(b,x2));
		          return (x);      
		   }
		   
		
/*
       !*****************************************************************************80
!
!! NORMAL_01_PDF evaluates the Normal 01 PDF.
!
!  Discussion:
!
!    The Normal 01 PDF is also called the "Standard Normal" PDF, or
!    the Normal PDF with 0 mean and variance 1.
!
!    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF. 
*/

  
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d  
		      normal_01_pdf_zmm8r8(const __m512d x) {
		         
		          register __m512d C039894228040143267793995 = 
		                                   _mm512_set1_pd(0.39894228040143267793995);
		          register __m512d C05 = _mm512_set1_pd(-0.5);
		          register __m512d earg,pdf;
		          earg = _mm512_mul_pd(C05,_mm512_mul_pd(x,x));
#if (USE_SLEEF_LIB) == 1 
                          pdf  = _mm512_mul_pd(xexp(earg),C039894228040143267793995);
#else
                          pdf  = _mm512_mul_pd(_mm512_exp_pd(earg),C039894228040143267793995);
#endif                         		   
                          return (pdf);
		     }  
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 
		      normal_01_pdf_zmm16r4(const __m512 x) {
		         
		          register __m512 C039894228040143267793995 = 
		                                   _mm512_set1_ps(0.39894228040143267793995f);
		          register __m512 C05 = _mm512_set1_ps(-0.5f);
		          register __m512 earg,pdf;
		          earg = _mm512_mul_ps(C05,_mm512_mul_ps(x,x));
#if (USE_SLEEF_LIB) == 1 
                          pdf  = _mm512_mul_ps(xexpf(earg),C039894228040143267793995);
#else
                          pdf  = _mm512_mul_ps(_mm512_exp_ps(earg),C039894228040143267793995);
#endif                         		   
                          return (pdf);
		     }  
		     
		     
/*
      !*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1   
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d 
		      uniform_01_zmm8r8( __m512i & seed) {
		         
		         register __m512i C127773 = _mm512_set1_epi32(127773);
		         register __m512i C16807  = _mm512_set1_epi32(16807);
		         register __m512i C2836   = _mm512_set1_epi32(2836);
		         register __m512d C4656612875 = _mm512_set1_pd(4.656612875e-10);
		         register __m512i k,t0,t1;
		         register __m512d uni01;
		         k  = _mm512_div_epi32(seed,C127773);
		         t0 = _mm512_mul_epi32(k,C2836);
		         t1 = _mm512_sub_epi32(seed,_mm512_mul_epi32(k,C127773));
		         seed = _mm512_mul_epi32(C16807,_mm512sub_epi32(t1,t0));
		         if(_mm512_cmp_epi32_mask(seed,_mm512_setzero_epi32(),_CMP_LT_OQ)) 
		            seed = _mm512_add_epi32(seed,_mm512_set1_epi32(2147483647));
		         uni01   = _mm512_mul_pd(_mm512_castsi512_pd(seed),C4656612875);
		         return (uni01);
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 
		      uniform_01_zmm16r4( __m512i & seed) {
		            
		            return (_mm512_castpd_ps(uniform_01_zmm8r8(seed)));
		      }
		     
/*
 !*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
*/	

 
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d 
	              uniform_zmm8r8(const __m512d a,
	                             const __m512d b,
	                             __m512i & seed) {
	                             
	                 register __m512i C127773 = _mm512_set1_epi32(127773);
		         register __m512i C16807  = _mm512_set1_epi32(16807);
		         register __m512i C2836   = _mm512_set1_epi32(2836);
		         register __m512d C4656612875 = _mm512_set1_pd(4.656612875e-10);
		         register __m512i k,t0,t1;
		         register __m512d uni,diff;
		         k  = _mm512_div_epi32(seed,C127773);
		         t0 = _mm512_mul_epi32(k,C2836);
		         t1 = _mm512_sub_epi32(seed,_mm512_mul_epi32(k,C127773));
		         seed = _mm512_mul_epi32(C16807,_mm512sub_epi32(t1,t0));
		         if(_mm512_cmp_epi32_mask(seed,_mm512_setzero_epi32(),_CMP_LT_OQ)) 
		            seed = _mm512_add_epi32(seed,_mm512_set1_epi32(2147483647));
		         diff    = _mm512_add_pd(_mm512_sub_pd(b,a));
		         uni     = _mm512_mul_pd(diff,_mm512_mul_pd(_mm512_castsi512_pd(seed),C4656612875);
		         return (uni);             
	            }
	            
	            
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 
	              uniform_zmm16r4(const __m512d a,
	                              const __m512d b,
	                              __m512i & seed) {
	                   
	                   return (_mm512_castpd_ps(uniform_zmm8r8(a,b,seed)));                 
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


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d  
                      gamma_incomplete_zmm8r8(const __m512d p,
                                              const __m512d x) {
                                              
                            const __m512d exp_arg_min = _mm512_set1_pd(-88.0e+00);
                            const __m512d overflow    = _mm512_set1_pd(1.0e+37);
                            const __m512d plimit      = _mm512_set1_pd(1000.0e+00);
                            const __m512d tol         = _mm512_set1_pd(1.0e-7);
                            const __m512d xbig        = _mm512_set1_pd(1.0e+8);
                            const __m512d C0          = _mm512_setzero_pd();
                            const __m512d C0333333333 = _mm512_set1_pd(0.3333333333333333333333);
                            const __m512d C1          = _mm512_set1_pd(1.0);
                            const __m512d C2          = _mm512_set1_pd(2.0);
                            const __m512d C3          = _mm512_set1_pd(3.0);
                            const __m512d C9          = _mm512_set1_pd(9.0);
                            __m512d cdf,arg,b,c;
                            __m512d pn1,pn2,pn3,pn4;
                            __m512d pn5,pn6,rn,t0,t1;
                            __m512d gaminc;
                            __mmask8 m0,m1;
                            m0 = _mm512_cmp_pd_mask(plimit,p,_CMP_LT_OQ);
                            if(m0) {
                               __m512d sqrp,xp,_9p1,t0,t1;
                               xp     = _mm512_div_pd(x,p);
                               _9p1   = _mm512_fmsub_pd(C9,p,C1);
                               sqrp   = _mm512_mul_pd(C3,_mm512_sqrt_pd(p));
                               t0     = _mm512_pow_pd(xp,C0333333333);
                               t1     = _mm512_add_pd(t0,
                                            _mm512_div_pd(C1,_9p1));
                               pn1    = _mm512_mul_pd(sqrp,t1);
                               gaminc = normal_01_cdf_zmm8r8(pn1);
                               return (gaminc);
                            }   
                            m0 = _mm512_cmp_pd_mask(x,C1,_CMP_LE_OQ);
                            m1 = _mm512_cmp_pd_mask(x,p,_CMP_LT_OQ);
                            if(m0 || m1) {
#if (USE_SLEEF_LIB) == 1 
                               t0  = xlog(x);
#else
                               t0  = _mm512_log_pd(x);
#endif                               
                               t1  = gamma_log_zmm8r8(_mm512_add_pd(p,C1));
                               arg = _mm512_fmsub_pd(p,t0,_mm512_sub_pd(x,t1));
                               c   = C1;
                               gaminc = C1;
                               a   = p; 
                               while(true) {
                                    a      = _mm512_add_pd(a,C1);
                                    c      = _mm512_mul_pd(c,_mm512_div_pd(x,a));
                                    gaminc = _mm512_add_pd(gaminc,c);
                                    m0     = _mm512_cmp_pd_mask(c,tol,_CMP_LE_OQ);
                                    if(m0) break;
                               }
#if (USE_SLEEF_LIB) == 1 
                               t0  = xlog(x);
#else
                               t0  = _mm512_log_pd(x);
#endif                         
                               arg = _mm512_add_pd(arg,t0);  
                               m1  = _mm512_cmp_pd_mask(exp_arg_min,arg,_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1 
                               gaminc = _mm512_mask_blend_pd(m1,C0,xexp(arg));  
#else
                               gaminc = _mm512_mask_blend_pd(m1,C0,_mm512_exp_pd(arg));  
#endif                                  
                           } 
                           else {
#if (USE_SLEEF_LIB) == 1 
                               t0  = xlog(x);
#else
                               t0  = _mm512_log_pd(x);
#endif                               
                               t1  = gamma_log_zmm8r8(p);
                               arg = _mm512_fmsub_pd(p,t0,_mm512_sub_pd(x,t1));                               
                               a   = _mm512_sub_pd(C1,p);
                               b   = _mm512_add_pd(a,_mm512_add_pd(x,C1));
                               c   = C0;
                               pn1 = C1;
                               pn2 = x;
                               pn3 = _mm512_add_pd(x,C1);
                               pn4 = _mm512_mul_pd(x,b);
                               gaminc = _mm512_div_pd(pn3,pn4);
                               while(true) {
                                   a = _mm512_add_pd(a,C1);
                                   b = _mm512_add_pd(b,C2);
                                   c = _mm512_add_pd(c,C1);
                                   pn5 = _mm512_fmsub_pd(b,pn3,
                                                     _mm512_mul_pd(a,
                                                           _mm512_mul_pd(c,pn1)));
                                   pn6 = _mm512_fmsub_pd(b,pn4,
                                                     _mm512_mul_pd(a,
                                                           _mm512_mul_pd(c,pn2)));
                                   if(_mm512_cmp_pd_mask(C0,_mm512_abs_pd(pn6),
                                                                       _CMP_LT_OQ)) {
                                        rn = _mm512_div_pd(pn5,pn6);
                                        t0 = _mm512_abs_pd(_mm512_sub_pd(gaminc,rn));
                                        t1 = _mm512_min_pd(tol,_mm512_mul_pd(tol,rn));
                                        if(_mm512_cmp_pd_mask(t0,t1,_CMP_LE_OQ)) {
#if (USE_SLEEF_LIB) == 1 
                                           arg  = _mm512_add_pd(arg,xlog(gaminc));
#else
                                           arg  = _mm512_add_pd(_mm512_log_pd(gaminc));
#endif       
                                           m1   = _mm512_cmp_pd_mask(exp_arg_min,arg,_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1 
                                           gaminc = _mm512_mask_blend_pd(m1,C1,_mm512_sub_pd(C1,xexp(arg)));
#else
                                           gaminc = _mm512_mask_blend_pd(m1,C1,_mm512_sub_pd(C1,
                                                                                        _mm512_exp_pd(arg)));
#endif       
                                           return (gaminc);                               
                                        }    
                                        gaminc = rn;                               
                                   }
                                   pn1 = pn3;
                                   pn2 = pn4;
                                   pn3 = pn5;
                                   pn4 = pn6;
                                   if(_mm512_cmp_pd_mask(overflow,
                                                   _mm512_abs_pd(pn5),_CMP_LE_OQ)) {
                                      t0 = _mm512_div_pd(C1,overflow);
                                      pn1= _mm512_mul_pd(pn1,t0);
                                      pn2= _mm512_mul_pd(pn2,t0);
                                      pn3= _mm512_mul_pd(pn3,t0);
                                      pn4= _mm512_mul_pd(pn4,t0);               
                                   }
                               }
                           } 
                           
                           return (gaminc);                   
                   }
                   
#include <limits>
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
                      gamma_incomplete_zmm16r4(const __m512 p,
                                              const __m512 x) {
                                              
                            const __m512 exp_arg_min = _mm512_set1_ps(-88.0e+00f);
                            const __m512 overflow    = _mm512_set1_ps(std::numeric_limits<float>::max());
                            const __m512 plimit      = _mm512_set1_ps(1000.0e+00f);
                            const __m512 tol         = _mm512_set1_ps(1.0e-7f);
                            const __m512 xbig        = _mm512_set1_ps(1.0e+8f);
                            const __m512 C0          = _mm512_setzero_ps();
                            const __m512 C0333333333 = _mm512_set1_ps(0.3333333333333333333333f);
                            const __m512 C1          = _mm512_set1_ps(1.0f);
                            const __m512 C2          = _mm512_set1_ps(2.0f);
                            const __m512 C3          = _mm512_set1_ps(3.0f);
                            const __m512 C9          = _mm512_set1_ps(9.0f);
                            __m512 cdf,arg,b,c;
                            __m512 pn1,pn2,pn3,pn4;
                            __m512 pn5,pn6,rn,t0,t1;
                            __m512 gaminc;
                            __mmask16 m0,m1;
                            m0 = _mm512_cmp_ps_mask(plimit,p,_CMP_LT_OQ);
                            if(m0) {
                               __m512 sqrp,xp,_9p1,t0,t1;
                               xp     = _mm512_div_ps(x,p);
                               _9p1   = _mm512_fmsub_ps(C9,p,C1);
                               sqrp   = _mm512_mul_ps(C3,_mm512_sqrt_ps(p));
                               t0     = _mm512_pow_ps(xp,C0333333333);
                               t1     = _mm512_add_ps(t0,
                                            _mm512_div_ps(C1,_9p1));
                               pn1    = _mm512_mul_ps(sqrp,t1);
                               gaminc = normal_01_cdf_zmm16r4(pn1);
                               return (gaminc);
                            }   
                            m0 = _mm512_cmp_ps_mask(x,C1,_CMP_LE_OQ);
                            m1 = _mm512_cmp_ps_mask(x,p,_CMP_LT_OQ);
                            if(m0 || m1) {
#if (USE_SLEEF_LIB) == 1 
                               t0  = xlogf(x);
#else
                               t0  = _mm512_log_ps(x);
#endif                               
                               t1  = gamma_log_zmm16r4(_mm512_add_ps(p,C1));
                               arg = _mm512_fmsub_ps(p,t0,_mm512_sub_ps(x,t1));
                               c   = C1;
                               gaminc = C1;
                               a   = p; 
                               while(true) {
                                    a      = _mm512_add_ps(a,C1);
                                    c      = _mm512_mul_ps(c,_mm512_div_ps(x,a));
                                    gaminc = _mm512_add_ps(gaminc,c);
                                    m0     = _mm512_cmp_ps_mask(c,tol,_CMP_LE_OQ);
                                    if(m0) break;
                               }
#if (USE_SLEEF_LIB) == 1 
                               t0  = xlogf(x);
#else
                               t0  = _mm512_log_ps(x);
#endif                         
                               arg = _mm512_add_ps(arg,t0);  
                               m1  = _mm512_cmp_ps_mask(exp_arg_min,arg,_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1 
                               gaminc = _mm512_mask_blend_ps(m1,C0,xexpf(arg));  
#else
                               gaminc = _mm512_mask_blend_ps(m1,C0,_mm512_exp_ps(arg));  
#endif                                  
                           } 
                           else {
#if (USE_SLEEF_LIB) == 1 
                               t0  = xlogf(x);
#else
                               t0  = _mm512_log_ps(x);
#endif                               
                               t1  = gamma_log_zmm16r4(p);
                               arg = _mm512_fmsub_ps(p,t0,_mm512_sub_ps(x,t1));                               
                               a   = _mm512_sub_ps(C1,p);
                               b   = _mm512_add_ps(a,_mm512_add_ps(x,C1));
                               c   = C0;
                               pn1 = C1;
                               pn2 = x;
                               pn3 = _mm512_add_ps(x,C1);
                               pn4 = _mm512_mul_ps(x,b);
                               gaminc = _mm512_div_ps(pn3,pn4);
                               while(true) {
                                   a = _mm512_add_ps(a,C1);
                                   b = _mm512_add_ps(b,C2);
                                   c = _mm512_add_ps(c,C1);
                                   pn5 = _mm512_fmsub_ps(b,pn3,
                                                     _mm512_mul_ps(a,
                                                           _mm512_mul_ps(c,pn1)));
                                   pn6 = _mm512_fmsub_ps(b,pn4,
                                                     _mm512_mul_ps(a,
                                                           _mm512_mul_ps(c,pn2)));
                                   if(_mm512_cmp_ps_mask(C0,_mm512_abs_ps(pn6),
                                                                       _CMP_LT_OQ)) {
                                        rn = _mm512_div_ps(pn5,pn6);
                                        t0 = _mm512_abs_ps(_mm512_sub_ps(gaminc,rn));
                                        t1 = _mm512_min_ps(tol,_mm512_mul_ps(tol,rn));
                                        if(_mm512_cmp_ps_mask(t0,t1,_CMP_LE_OQ)) {
#if (USE_SLEEF_LIB) == 1 
                                           arg  = _mm512_add_ps(arg,xlogf(gaminc));
#else
                                           arg  = _mm512_add_ps(_mm512_log_ps(gaminc));
#endif       
                                           m1   = _mm512_cmp_ps_mask(exp_arg_min,arg,_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1 
                                           gaminc = _mm512_mask_blend_ps(m1,C1,_mm512_sub_ps(C1,xexp(arg)));
#else
                                           gaminc = _mm512_mask_blend_ps(m1,C1,_mm512_sub_ps(C1,
                                                                                        _mm512_exp_ps(arg)));
#endif       
                                           return (gaminc);                               
                                        }    
                                        gaminc = rn;                               
                                   }
                                   pn1 = pn3;
                                   pn2 = pn4;
                                   pn3 = pn5;
                                   pn4 = pn6;
                                   if(_mm512_cmp_ps_mask(overflow,
                                                   _mm512_abs_ps(pn5),_CMP_LE_OQ)) {
                                      t0 = _mm512_div_ps(C1,overflow);
                                      pn1= _mm512_mul_ps(pn1,t0);
                                      pn2= _mm512_mul_ps(pn2,t0);
                                      pn3= _mm512_mul_ps(pn3,t0);
                                      pn4= _mm512_mul_ps(pn4,t0);               
                                   }
                               }
                           } 
                           
                           return (gaminc);                   
                   }
                   
                   
/*
 !*****************************************************************************80
!
!! GAMMA_CDF evaluates the Gamma CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2000
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
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      gamma_cdf_zmm8r8(const __m512d x,
		                       const __m512d a,
		                       const __m512d b,
		                       const __m512d c) {
		      
		          register __m512d x2,p2;
		          register __m512d cdf;
		          p2 = c;
		          x2 = _mm512_div_pd(_mm512_sub_pd(x,a),b);
		          cdf= gamma_incomplete_zmm8r8(p2,x2);
		          return (cdf);                  
		    }
		    
		    
		    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      gamma_cdf_zmm16r4(const __m512 x,
		                       const __m512 a,
		                       const __m512 b,
		                       const __m512 c) {
		      
		          register __m512 x2,p2;
		          register __m512 cdf;
		          p2 = c;
		          x2 = _mm512_div_ps(_mm512_sub_ps(x,a),b);
		          cdf= gamma_incomplete_zmm16r4(p2,x2);
		          return (cdf);                  
		    }
		    
         
                   
/*
   !*****************************************************************************80
!
!! CHI_CDF evaluates the Chi CDF.
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
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!                
*/


            
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      chi_cdf_zmm8r8(const __m512d x,
		                     const __m512d a,
		                     const __m512d b,
		                     const __m512d c) {
		       
		         const __m512d C05 = _mm512_set1_pd(0.5);
		         __m512d p2,x2,y;
		         __m512d cdf;
		         y  = _mm512_div_pd(_mm512_sub_pd(x,a),b);
		         x2 = _mm512_mul_pd(C05,_mm512_mul_pd(y,y));
		         p2 = _mm512_mul_pd(C05,c);
		         cdf= gamma_incomplete_zmm8r8(p2,x2);
		         return (cdf);               
		   }
		   
		   
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      chi_cdf_zmm16r4(const __m512 x,
		                      const __m512 a,
		                      const __m512 b,
		                      const __m512 c) {
		       
		         const __m512 C05 = _mm512_set1_ps(0.5f);
		         __m512 p2,x2,y;
		         __m512 cdf;
		         y  = _mm512_div_ps(_mm512_sub_ps(x,a),b);
		         x2 = _mm512_mul_ps(C05,_mm512_mul_ps(y,y));
		         p2 = _mm512_mul_ps(C05,c);
		         cdf= gamma_incomplete_zmm16r4(p2,x2);
		         return (cdf);               
		   }
		   
	
/*
      !*****************************************************************************80
!
!! CHI_SQUARE_CDF evaluates the Chi squared CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the random deviate.
!
!    Input, real ( kind = 8 ) A, the parameter of the distribution, usually
!    the number of degrees of freedom.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      chi_square_cdf_zmm8r8(const __m512d x,
		                            const __m512d a) {
		          
		          const __m512d C05 = _mm512_set1_pd(0.5);
		          const __m512d C1  = _mm512_set1_pd(1.0);
		          register __m512d x2,a2,b2,c2;
		          register __m512d cdf;
		          a2 =  _mm512_setzero_pd();
		          x2 =  _mm512_mul_pd(C05,C1);
		          b2 =  C1;
		          c2 =  _mm512_mul_pd(C05,a);
		          cdf=  gamma_cdf_zmm8r8(x2,a2,b2,c2);
		          return (cdf);                   
		     }
		     
		     
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      chi_square_cdf_zmm16r4(const __m512 x,
		                             const __m512 a) {
		          
		          const __m512 C05 = _mm512_set1_ps(0.5f);
		          const __m512 C1  = _mm512_set1_ps(1.0f);
		          register __m512 x2,a2,b2,c2;
		          register __m512 cdf;
		          a2 =  _mm512_setzero_ps();
		          x2 =  _mm512_mul_ps(C05,C1);
		          b2 =  C1;
		          c2 =  _mm512_mul_ps(C05,a);
		          cdf=  gamma_cdf_zmm16r4(x2,a2,b2,c2);
		          return (cdf);                   
		     }
		     
		   
/*
    !*****************************************************************************80
!
!! CHI_MEAN returns the mean of the Chi PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) MEAN, the mean value.
!   
*/

                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
                      chi_mean_zmm8r8(const __m512d a,
                                      const __m512d b,
                                      const __m512d c) {
                          
                          const __m512d C05 = _mm512_set1_pd(0.5);
                          const __m512d C1  = _mm512_set1_pd(1.0);
                          const __m512d C2  = _mm512_set1_pd(2.0);
                          const __m512d C141421356237309504880169 = 
                                              _mm512_set1_pd(1.41421356237309504880169);
                          register __m512d g0,g1,arg,t0;
                          register __m512d mean;
                          arg = _mm512_mul_pd(C05,_mm512_add_pd(c,C1));
                          g0  = gamma_zmm8r8(arg);
                          t0  = _mm512_mul_pd(C141421356237309504880169,b);
                          g1  = gamma_zmm8r8(_mm512_mul_pd(C05,c));
                          mean= _mm512_add_pd(a,_mm512_div_pd(_mm512_mul_pd(t0,g0),g1));
                          return (mean);
                    }	
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
                      chi_mean_zmm16r4(const __m512 a,
                                      const __m512 b,
                                      const __m512 c) {
                          
                          const __m512 C05 = _mm512_set1_ps(0.5f);
                          const __m512 C1  = _mm512_set1_ps(1.0f);
                          const __m512 C2  = _mm512_set1_ps(2.0f);
                          const __m512 C141421356237309504880169 = 
                                              _mm512_set1_ps(1.41421356237309504880169f);
                          register __m512 g0,g1,arg,t0;
                          register __m512 mean;
                          arg = _mm512_mul_ps(C05,_mm512_add_pd(c,C1));
                          g0  = gamma_zmm16r4(arg);
                          t0  = _mm512_mul_ps(C141421356237309504880169,b);
                          g1  = gamma_zmm16r4(_mm512_mul_ps(C05,c));
                          mean= _mm512_add_ps(a,_mm512_div_ps(_mm512_mul_ps(t0,g0),g1));
                          return (mean);
                    }	   
                    
                    
                    
/*
    !*****************************************************************************80
!
!! CHI_PDF evaluates the Chi PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = EXP ( - 0.5D+00 * ( ( X - A ) / B )^2 )
!      * ( ( X - A ) / B )^( C - 1 ) /
!      ( 2^( 0.5D+00 * C - 1 ) * B * GAMMA ( 0.5D+00 * C ) )
!
!    CHI(A,B,1) is the Half Normal PDF;
!    CHI(0,B,2) is the Rayleigh PDF;
!    CHI(0,B,3) is the Maxwell PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 December 1999
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
!    0 < B,
!    0 < C.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
*/


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d   
                      chi_pdf_zmm8r8(const __m512d x,
                                     const __m512d a,
                                     const __m512d b,
                                     const __m512d c) {
                          
                          const __m512d C05 = _mm512_set1_pd(0.5);
                          const __m512d C1  = _mm512_set1_pd(1.0);
                          const __m512d C2  = _mm512_set1_pd(2.0);
                          register __m512d y,t0,t1,t2,g0,exp;
                          register __m512d pdf;
                          t0 = _mm512_mul_pd(_mm512_mul_pd(negate_zmm8r8(C05,y),y),y);
                          t1 = _mm512_fmsub_pd(C05,c,C1);
                          y  = _mm512_div_pd(_mm512_sub_pd(x,a),b);
                          g0 = gamma_zmm8r8(_mm512_mul_pd(C05,c));  
#if (USE_SLEEF_LIB) == 1                           
                          exp= _mm512_mul_pd(xexp(t0),_mm512_pow_pd(_mm512_sub_pd(c,C1)));
#else
                          exp= _mm512_mul_pd(_mm512_exp_pd(t0),_mm512_pow_pd(_mm512_sub_pd(c,C1)));
#endif             
                          t2 = _mm512_mul_pd(_mm512_pow_pd(C2,t1),
                                             _mm512_mul_pd(b,g0));
                          pdf = _mm512_div_pd(exp,t2);
                          return (pdf);
                    }
                    
                    
                       __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512   
                      chi_pdf_zmm16r4(const __m512 x,
                                     const __m512 a,
                                     const __m512 b,
                                     const __m512 c) {
                          
                          const __m512 C05 = _mm512_set1_ps(0.5);
                          const __m512 C1  = _mm512_set1_ps(1.0);
                          const __m512 C2  = _mm512_set1_ps(2.0);
                          register __m512 y,t0,t1,t2,g0,exp;
                          register __m512 pdf;
                          t0 = _mm512_mul_ps(_mm512_mul_ps(negate_zmm16r4(C05,y),y),y);
                          t1 = _mm512_fmsub_ps(C05,c,C1);
                          y  = _mm512_div_ps(_mm512_sub_ps(x,a),b);
                          g0 = gamma_zmm16r4(_mm512_mul_ps(C05,c));  
#if (USE_SLEEF_LIB) == 1                           
                          exp= _mm512_mul_ps(xexpf(t0),_mm512_pow_ps(_mm512_sub_ps(c,C1)));
#else
                          exp= _mm512_mul_ps(_mm512_exp_ps(t0),_mm512_pow_ps(_mm512_sub_ps(c,C1)));
#endif             
                          t2 = _mm512_mul_ps(_mm512_pow_ps(C2,t1),
                                             _mm512_mul_ps(b,g0));
                          pdf = _mm512_div_ps(exp,t2);
                          return (pdf);
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
		      __m512d beta_zmm8r8(const __m512d a,
		                          const __m512d b) {

                       // const __m512d nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			const __m512d _0  = _mm512_setzero_pd();
			//if(__builtin_expect(_mm512_cmp_pd_mask(a,_0,_CMP_LE_OQ),0) ||
			//   __builtin_expect(_mm512_cmp_pd_mask(b,_0,_CMP_LE_OQ),0)) {
                       //    return (nan);
			//}
			const __m512d ab  = _mm512_add_pd(a,b);
			__m512d beta      = _mm512_setzero_pd();
			beta              = _mm512_exp_pd(
			                              _mm512_sub_pd(
						                 _mm512_add_pd(gamma_log_zmm8r8(a),
								               gamma_log_zmm8r8(b)),
			return (beta);						                     gamma_log_zmm8r8(ab)));
		    }
		    

/*
     !*****************************************************************************80
!
!! BETA_SAMPLE samples the Beta PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Algorithm BN,
!    Statistical Computing,
!    Dekker, 1980.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < A,
!    0.0D+00 < B.
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
		      __m512d 
                      beta_sample_zmm8r8(const __m512d a,
                                         const __m512d b,
                                         __m512i & seed) {
                         
                          __m512d C1 = _mm512_set1_pd(1.0);
                          __m512d C2 = _mm512_set1_pd(2.0);
                          __m512d C05= _mm512_set1_pd(0.5);
                          register __m512d mu,stdev,test,u,y;
                          register __m512d a1,b1,ab2,l0,l1,l2,t0;
                          register __m512d t1,t2,t3;
                          register __m512d x;
                          ab2 = _mm512_sub_pd(_mm512_add_pd(a,b),C2);
                          a1  = _mm512_sub_pd(a,C1);
                          b1  = _mm512_sub_pd(b,C1);
                          mu  = _mm512_div_pd(a1,ab2);
                          stdev = _mm512_div_pd(C05,_mm512_sqrt_pd(ab2));
                          while(true) {
                              
                              y = normal_01_sample_zmm8r8(seed);
                              x = _mm512_add_pd(mu,_mm512_mul_pd(stdev,y));
                              if(_mm512_cmp_pd_mask(x,_mm512_setzero_pd(),_CMP_LT_OQ) || 
                                 _mm512_cmp_pd_mask(C1,x,_CMP_LT_OQ)) continue;
                              t0= _mm512_mul_pd(C05,_mm512_mul_pd(y,y));
                              u = uniform_01_zmm8r8(seed);
#if (USE_SLEEF_LIB) == 1    
                              l0 = xlog(_mm512_div_pd(x,a1));
                              t1 = _mm512_mul_pd(a1,l0);
                              l1 = xlog(_mm512_div_pd(_mm512_sub_pd(C1,x),b1));
                              t2 = _mm512_mul_pd(b1,l1);
                              l2 = _mm512_add_pd(xlog(ab2),t0); 
                              t3 = _mm512_mul_pd(ab2,l2);
                              test = _mm512_add_pd(t1,_mm512_add_pd(t2,t3));
#else
                              l0 = _mm512_log_pd(_mm512_div_pd(x,a1));
                              t1 = _mm512_mul_pd(a1,l0);
                              l1 = _mm512_log_pd(_mm512_div_pd(_mm512_sub_pd(C1,x),b1));
                              t2 = _mm512_mul_pd(b1,l1);
                              l2 = _mm512_add_pd(_mm512_log_pd(ab2),t0); 
                              t3 = _mm512_mul_pd(ab2,l2);
                              test = _mm512_add_pd(t1,_mm512_add_pd(t2,t3));
#endif                      
#if (USE_SLEEF_LIB) == 1  
                              if(_mm512_cmp_pd_mask(xlog(u),test,_CMP_LE_OQ)) break;
#else
                              if(_mm512_cmp_pd_mask(_mm512_log_pd(u),test,_CMP_LE_OQ)) break;
#endif                                     
                          }  
                          return (x);                 
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
		      __m512d anglit_sample_zmm8r8() {

                         __m512d cdf;
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
                            const __m512d nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_uniform_distribution_double(0.0,1.0);
			 const double * __restrict ptr = (const double*)(&svrng_generate8_double(engine,uniform));
			 cdf              = anglit_cdf_inv_zmm8r8(_mm512_loadu_pd(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d anglit_sample_zmm8r8(const __m512 cdf) {

                            return (anglit_cdf_inv_zmm8r8(cdf));
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
		      __m512d arcsin_cdf_zmm8r8(const __m512d x,
		                                const __m512d a) {

                         const __m512d invpi = _mm512_set1_pd(0.318309886183790671537767526745);
			 const __m512d _0    = _mm512_setzero_pd();
			 const __m512d _1_2  = _mm512_set1_pd(0.5);
			 const __m512d _1    = _mm512_set1_pd(1.0);
			 __m512d t0,cdf;
			 __mmask8 m0,m1;
			 m0  = _mm512_cmp_pd_mask(x,zmm8r8_negate(a),_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1
                         t0  = _mm512_mul_pd(xasin(_mm512_div_pd(x,a),invpi));
#else
                         t0  = _mm512_mul_pd(_mm512_asin_pd(_mm512_div_pd(x,a),invpi));
#endif
			 m1  = _mm512_cmp_pd_mask(x,a,_CMP_LT_OQ);
                         cdf = _mm512_mask_blend_pd(m0,_mm512_add_pd(_1_2,t0),_0);
			 cdf = _mm512_mask_blend_pd(m1,cdf,_1); 
                         return (cdf);
		   }


		   
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 arcsin_cdf_zmm16r4(const __m512 x,
		                                const __m512 a) {

                         const __m512 invpi = _mm512_set1_ps(0.318309886183790671537767526745f);
			 const __m512 _0    = _mm512_setzero_ps();
			 const __m512 _1_2  = _mm512_set1_ps(0.5f);
			 const __m512 _1    = _mm512_set1_ps(1.0f);
			 __m512 t0,cdf;
			 __mmask16 m0,m1;
			 m0  = _mm512_cmp_ps_mask(x,zmm16r4_negate(a),_CMP_LE_OQ);
                         t0  = _mm512_mul_ps(_mm512_asin_ps(_mm512_div_ps(x,a),invpi));
			 m1  = _mm512_cmp_ps_mask(x,a,_CMP_LT_OQ);
                         cdf = _mm512_mask_blend_ps(m0,_mm512_add_pd(_1_2,t0),_0);
			 cdf = _mm512_mask_blend_ps(m1,cdf,_1); 
                         return (cdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d arcsin_cdf_inv_zmm8r8(const __m512d cdf,
		                                    const __m512d a) {

                           const __m512d pi    = _mm512_set1_pd(3.14159265358979323846264338328);
			   const __m512d _0    = _mm512_setzero_pd();
			   const __m512d _1    = _mm512_set1_pd(1.0);
			  // const __m512d nan   = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			   const __m512d _1_2  = _mm512_set1_pd(0.5);
			   __m512d x;
			  // if(__builtin_expect(_mm512_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			  //    __builtin_expect(_mm512_cmp_pd_mask(cdf,_1,_CMP_GT_OQ),0)) {
                          //    return (nan);
			  // }
#if (USE_SLEEF_LIB) == 1
                             x = _mm512_mul_pd(xsin(_mm512_mul_pd(pi,_mm512_sub_pd(cdf,_1_2))));
			     
#else
                             x = _mm512_mul_pd(_mm512_sin_pd(_mm512_mul_pd(pi,_mm512_sub_pd(cdf,_1_2))));
#endif
                             return (x);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 arcsin_cdf_inv_zmm16r4(const __m512 cdf,
		                                    const __m512 a) {

                           const __m512 pi    = _mm512_set1_ps(3.14159265358979323846264338328f);
			   const __m512 _0    = _mm512_setzero_ps();
			   const __m512 _1    = _mm512_set1_ps(1.0f);
			  // const __m512 nan   = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
			   const __m512 _1_2  = _mm512_set1_ps(0.5f);
			   __m512 x;
			 //  if(__builtin_expect(_mm512_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			 //     __builtin_expect(_mm512_cmp_ps_mask(cdf,_1,_CMP_GT_OQ),0)) {
                         //     return (nan);
			 //  }
#if (USE_SLEEF_LIB) == 1
                             x = _mm512_mul_ps(xsinf(_mm512_mul_ps(pi,_mm512_sub_ps(cdf,_1_2))));
			     
#else
                             x = _mm512_mul_ps(_mm512_sin_ps(_mm512_mul_ps(pi,_mm512_sub_ps(cdf,_1_2))));
#endif
                             return (x);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m512d arcsin_mean_zmm8r8() {

		            return (_mm512_setzero_pd());
		      }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m512 arcsin_mean_zmm16r4() {

		            return (_mm512_setzero_ps());
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
		      __m512d arcsin_pdf_zmm8r8(const __m512d x,
		                                const __m512d a) {

                           const __m512d pi    = _mm512_set1_pd(3.14159265358979323846264338328);
			   const __m512d _0    = _mm512_setzero_pd();
			   const __m512d _1    = _mm512_set1_pd(1.0);
			  // const __m512d nan   = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			   __m512d pdf,t0;
			   __mmask8 m,m1;
			 //  if(__builtin_expect(_mm512_cmp_pd_mask(a,_0,_CMP_LE_OQ))) {
                         //      return (nan);
			 //  }
			   m  =  _mm512_cmp_pd_mask(x,zmm8r8_negate(a),_CMP_LE_OQ);
			   t0 =  _mm512_sqrt_pd(_mm512_sub_pd(_mm512_mul_pd(a,a),
			                                      _mm512_mul_pd(x,x)));
			   m1 = _mm512_cmp_pd_mask(x,a,_CMP_GE_OQ);
			   __mmask8 m2 = m || m1;
			   pdf = _mm512_mask_blend_pd(m2,_mm512_div_pd(_1,
			                                           _mm512_mul_pd(pi,t0)),_0);
			   return (pdf);
			   
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512 arcsin_pdf_zmm16r4(const __m512 x,
		                                const __m512 a) {

                           const __m512 pi    = _mm512_set1_ps(3.14159265358979323846264338328f);
			   const __m512 _0    = _mm512_setzero_ps();
			   const __m512 _1    = _mm512_set1_ps(1.0f);
			  // const __m512 nan   = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
			   __m512 pdf,t0;
			   __mmask 16m,m1;
			  // if(__builtin_expect(_mm512_cmp_ps_mask(a,_0,_CMP_LE_OQ))) {
                         //      return (nan);
			 //  }
			   m  =  _mm512_cmp_ps_mask(x,zmm16r4_negate(a),_CMP_LE_OQ);
			   t0 =  _mm512_sqrt_ps(_mm512_sub_ps(_mm512_mul_ps(a,a),
			                                      _mm512_mul_ps(x,x)));
			   m1 = _mm512_cmp_ps_mask(x,a,_CMP_GE_OQ);
			   const __mmask16 m2 = m || m1;
			   pdf = _mm512_mask_blend_ps(m2,_mm512_div_ps(_1,
			                                           _mm512_mul_ps(pi,t0)),_0);
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
                       __m512d arcsin_variance_zmm8r8(const __m512d a) {

                         const __m512d _1_2 = _mm512_set1_pd(0.5);
			 __m512d variance;
			 variance = _mm512_mul_pd(a,_mm512_mul_pd(a,_1_2));
			 return (variance);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                       __m512 arcsin_variance_zmm16r4(const __m512 a) {

                         const __m512 _1_2 = _mm512_set1_ps(0.5f);
			 __m512 variance;
			 variance = _mm512_mul_ps(a,_mm512_mul_ps(a,_1_2));
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
		      __m512d arcsin_sample_zmm8r8() {

                         __m512d cdf;
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
                            const __m512d nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_normal_distribution_double(0.0,1.0);
			 const double * __restrict ptr = (const double*)(&svrng_generate8_double(engine,uniform));
			 cdf              = arcsin_cdf_inv_zmm8r8(_mm512_loadu_pd(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d arcsin_sample_zmm8r8(const __m512 cdf) {

                            return (arcsin_cdf_inv_zmm8r8(cdf));
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
		      __m512d beta_binomial_cdf_zmm8r8(const int32_t x,
		                                       const int32_t c,
						       const __m512d a,
						       const __m512d b) {

			      const __m512d _0  = _mm512_setzero_pd();
                              const __m512d _1  = _mm512_set1_pd(1.0);
			      __m512d vx,vy,vcy,vc1,vy1,vcy1;
			      __m512d cdf,pdf;

			      if(x<0) {
                                 cdf = _0;
			      }
			      else if(x<c) {
                                 cdf = _0;
				 for(int32_t y = 0; y < x; ++y) {
                                     vy  = _mm512_set1_pd((double)y);
				     vx  = _mm512_set1_pd((double)x);
				     vcy = _mm512_set1_pd((double)(c-y));
				     vc1 = _mm512_set1_pd((double)(c+1));
				     vy1 = _mm512_set1_pd((double)(y+1));
				     vcy1= _mm512_set1_pd((double)(c-y+1));
				     const __m512d t0 = beta_zmm8r8(_mm512_add_pd(a,vy),
				                                    _mm512_add_pd(b,vcy));
				     const __m512d t1 = _mm512_mul_pd(vc1,beta_zmm8r8(vy1,vcy1));
				     const __m512d t2 = beta_zmm8r8(a,b);
				     pdf              = _mm512_div_pd(t0,_mm512_mul_pd(t1,t2));
				     cdf              = _mm512_add_pd(cdf,pdf);
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
		      __m512d beta_pdf_zmm8r8(const __m512d x,
		                              const __m512d a,
					      const __m512d b) {

                         const __m512d _0 = _mm512_setzero_pd();
			 const __m512d _1 = _mm512_set1_pd(1.0);
			 const __m512d t0 = _mm512_sub_pd(a,_1);
			 const __m512d t1 = _mm512_sub_pd(_1,x);
			 const __m512d t2 = _mm512_sub_pd(b,_1);
			 __m512d pdf,term1,term2,term3;
			 __mmask8 m0,m1,m2;
			 term1            = _mm512_pow_pd(x,t0);
			 m0               = _mm512_cmp_pd_mask(x,_0,_CMP_LT_OQ);
			 term2            = _mm512_mul_pd(term1,_mm512_pow_pd(t1,t2));
			 m1               = _mm512_cmp_pd_mask(x,_1,_CMP_LT_OQ);
			 term3            = _mm512_div_pd(term2,beta_zmm8r8(a,b));
			 m                = m1||m2;
			 pdf              = _mm512_mask_blend_pd(m,term3,_0);
			 return (pdf);
		    }


		     __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      beta_pdf_zmm16r4(const __m512 x,
		                       const __m512 a,
				       const __m512 b) {

                         const __m512 _0 = _mm512_setzero_ps();
			 const __m512 _1 = _mm512_set1_ps(1.0);
			 const __m512 t0 = _mm512_sub_ps(a,_1);
			 const __m512 t1 = _mm512_sub_ps(_1,x);
			 const __m512 t2 = _mm512_sub_ps(b,_1);
			 __m512 pdf,term1,term2,term3;
			 __mmask16 m0,m1,m2;
			 term1            = _mm512_pow_ps(x,t0);
			 m0               = _mm512_cmp_ps_mask(x,_0,_CMP_LT_OQ);
			 term2            = _mm512_mul_ps(term1,_mm512_pow_pd(t1,t2));
			 m1               = _mm512_cmp_ps_mask(x,_1,_CMP_LT_OQ);
			 term3            = _mm512_div_ps(term2,beta_zmm16r4(a,b));
			 m                = m1||m2;
			 pdf              = _mm512_mask_blend_ps(m,term3,_0);
			 return (pdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      beta_variance_zmm8r8(const __m512d a,
		                           const __m512d b) {

			  __m512d variance;
                          const __m512d _1  = _mm512_set1_pd(1.0);
			  const __m512d ab  = _mm512_add_pd(a,b);
			  const __m512d t0  = _mm512_mul_pd(_mm512_mul_pd(ab,ab),
			                                    _mm512_add_pd(_1,ab));
			  variance          = _mm512_div_pd(_mm512_mul_pd(a,b),t0);				   
			  
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      beta_variance_zmm16r4(const __m512 a,
		                            const __m512 b) {

			  __m512 variance;
                          const __m512 _1  = _mm512_set1_ps(1.0f);
			  const __m512 ab  = _mm512_add_ps(a,b);
			  const __m512 t0  = _mm512_mul_ps(_mm512_mul_ps(ab,ab),
			                                    _mm512_add_ps(_1,ab));
			  variance          = _mm512_div_ps(_mm512_mul_ps(a,b),t0);				   
			  
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
                      __m512d
		      weibull_cdf_zmm8r8(const __m512d x,
		                         const __m512d a,
					 const __m512d b,
					 const __m512d c) {

                          const __m512d  _0 = _mm512_setzero_pd();
			  const __m512d  _1 = _mm512_set1_pd(1.0);
			  const __m512d  y  = _mm512_div_pd(_mm512_sub_pd(x,a),b);
			  const __m512d  exc= _mm512_exp_pd(_mm512_pow_pd(y,c));
			  __m512d cdf;
			  const __mmask8 m  = _mm512_cmp_pd_mask(a,x,_CMP_LT_OQ);
			  cdf               = _mm512_mask_blend_pd(m,_mm512_sub_pd(_1,
			                                                       _mm512_div_pd(_1,exc)),_0);
			  return (cdf);
		   }
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m512
		      weibull_cdf_zmm16r4(const __m512 x,
		                          const __m512 a,
					  const __m512 b,
					  const __m512 c) {

                          const __m512  _0 = _mm512_setzero_ps();
			  const __m512  _1 = _mm512_set1_ps(1.0f);
			  const __m512  y  = _mm512_div_ps(_mm512_sub_ps(x,a),b);
			  const __m512  exc= _mm512_exp_ps(_mm512_pow_ps(y,c));
			  __m512 cdf;
			  const __mmask16 m  = _mm512_cmp_ps_mask(a,x,_CMP_LT_OQ);
			  cdf               = _mm512_mask_blend_ps(m,_mm512_sub_ps(_1,
			                                                       _mm512_div_ps(_1,exc)),_0);
			  return (cdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m512d
		      weibull_cdf_inv_zmm8r8(const __m512d a,
		                             const __m512d b,
					     const __m512d c,
					     const __m512d cdf) {

                        const __m512d  _0  = _mm512_setzero_pd();
			const __m512d  _1  = _mm512_set1_pd(1.0);
			//const __m512d  nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			__m512d t0,t1,x;
			//if(__builtin_expect(_mm512_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//   __builtin_expect(_mm512_cmp_pd_mask(cdf,_1,_CMP_GT_OQ),0)) {
                       //    return (nan);
			//}
			t0                 = zmm8r8_negate(_mm512_log_pd(_mm512_sub_pd(_1,cdf)));
			t1                 = _mm512_pow_pd(t0,_mm512_div_pd(_1,c));
			x                  = _mm512_fmadd_pd(a,b,t1);
			return (x);
			
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
                      __m512
		      weibull_cdf_inv_zmm16r4(const __m512 a,
		                             const __m512 b,
					     const __m512 c,
					     const __m512 cdf) {

                        const __m512  _0  = _mm512_setzero_ps();
			const __m512  _1  = _mm512_set1_ps(1.0f);
			//const __m512  nan = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
			__m512 t0,t1,x;
			//if(__builtin_expect(_mm512_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//   __builtin_expect(_mm512_cmp_ps_mask(cdf,_1,_CMP_GT_OQ),0)) {
                       //    return (nan);
			//}
			t0                 = zmm16r4_negate(_mm512_log_pd(_mm512_sub_ps(_1,cdf)));
			t1                 = _mm512_pow_ps(t0,_mm512_div_ps(_1,c));
			x                  = _mm512_fmadd_ps(a,b,t1);
			return (x);
			
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      weibull_sample_zmm8r8(const __m512d vrand,
		                            const __m512d a,
					    const __m512d b,
					    const __m512d c) {

                         return (weibull_cdf_zmm8r8(a,b,c,vrand));
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      weibull_sample_zmm16r4(const __m512 vrand,
		                            const __m512 a,
					    const __m512 b,
					    const __m512 c) {

                         return (weibull_cdf_zmm16r4(a,b,c,vrand));
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
		      __m512d
                      weibull_discrete_cdf_zmm8r8(const __m512d x,
		                                  const __m512d a,
					          const __m512d b) {

			    __m512d cdf;
                            const __m512d  _0 = _mm512_setzero_pd();
			    const __m512d  _1 = _mm512_set1_pd(1.0);
			    const __m512d  t0 = _mm512_pow_pd(_mm512_add_pd(x,_1),b);
			    const __m512d  t1 = _mm512_pow_pd(_mm512_sub_pd(_1,a),t0);
			    const __mmask8 m  = _mm512_cmp_pd_mask(x,_0,_CMP_LT_OQ);
			    cdf               = _mm512_mask_blend_pd(m,_mm512_sub_pd(_1,t1),_0);
			    return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
                      weibull_discrete_cdf_zmm16r4(const __m512 x,
		                                  const __m512 a,
					          const __m512 b) {

			    __m512 cdf;
                            const __m512  _0 = _mm512_setzero_ps();
			    const __m512  _1 = _mm512_set1_ps(1.0f);
			    const __m512  t0 = _mm512_pow_ps(_mm512_add_ps(x,_1),b);
			    const __m512  t1 = _mm512_pow_ps(_mm512_sub_ps(_1,a),t0);
			    const __mmask16 m  = _mm512_cmp_ps_mask(x,_0,_CMP_LT_OQ);
			    cdf               = _mm512_mask_blend_ps(m,_mm512_sub_pd(_1,t1),_0);
			    return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512d
		      weibull_discrete_pdf_zmm8r8(const __m512d x,
		                                  const __m512d a,
					          const __m512d b) {

                            __m512d pdf;
                            const __m512d  _0 = _mm512_setzero_pd();
			    const __m512d  _1 = _mm512_set1_pd(1.0);
			    const __m512d  t0 = _mm512_pow_pd(_mm512_add_pd(x,_1),b);
			    const __m512d  _1a= _mm512_sub_pd(_1,a);
			    const __m512d  t1 = _mm512_pow_pd(_1a,t0);
                            const __m512d  t2 = _mm512_pow_pd(x,b);
			    const __m512d  t3 = _mm512_pow_pd(_1a,t2);
			    pdf               = _mm512_sub_pd(t3,t1);
			    return (pdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      weibull_discrete_pdf_zmm16r4(const __m512 x,
		                                  const __m512 a,
					          const __m512 b) {

                            __m512 pdf;
                            const __m512  _0 = _mm512_setzero_ps();
			    const __m512  _1 = _mm512_set1_ps(1.0);
			    const __m512  t0 = _mm512_pow_ps(_mm512_add_ps(x,_1),b);
			    const __m512  _1a= _mm512_sub_ps(_1,a);
			    const __m512  t1 = _mm512_pow_ps(_1a,t0);
                            const __m512  t2 = _mm512_pow_ps(x,b);
			    const __m512  t3 = _mm512_pow_ps(_1a,t2);
			    pdf               = _mm512_sub_ps(t3,t1);
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
		      __m512d
		      weibull_discr_icdf_zmm8r8(const __m512d cdf,
		                                const __m512d a,
						const __m512d b) {

                        //  const __m512d  nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			  const __m512d  _0  = _mm512_setzero_pd();
			  const __m512d  _1  = _mm512_set1_pd(1.0);
			 // if(__builtin_expect(_mm512_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			 //    __builtin_expect(_mm512_cmp_pd_mask(cdf,_1,_CMP_GT_OQ),0)) {
                         //  return (nan);
			 // }
			  const __m512d t0   =  _mm512_log_pd(_mm512_sub_pd(_1,cdf));
			  const __m512d t1   =  _mm512_log_pd(_mm512_sub_pd(_1,a));
			  const __m512d t2   =  _mm512_div_pd(t1,t2)
			  const __m512d t3   =  _mm512_pow_pd(t2,_mm512_div_pd(_1,b));
			  __m512d x;
			  x                  =  _mm512_ceil_pd(_mm512_sub_pd(t3,_1));
			  return (x);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline
		      __m512
		      weibull_discr_icdf_zmm16r4(const __m512 cdf,
		                                const __m512 a,
						const __m512 b) {

                       //   const __m512  nan = _mm512_set1_pd(std::numeric_limits<float>::quiet_NaN());
			  const __m512  _0  = _mm512_setzero_ps();
			  const __m512  _1  = _mm512_set1_ps(1.0f);
			//  if(__builtin_expect(_mm512_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//     __builtin_expect(_mm512_cmp_ps_mask(cdf,_1,_CMP_GT_OQ),0)) {
                        //   return (nan);
			//  }
			  const __m512 t0   =  _mm512_log_ps(_mm512_sub_ps(_1,cdf));
			  const __m512 t1   =  _mm512_log_ps(_mm512_sub_ps(_1,a));
			  const __m512 t2   =  _mm512_div_ps(t1,t2)
			  const __m512 t3   =  _mm512_pow_ps(t2,_mm512_div_ps(_1,b));
			  __m512 x;
			  x                  =  _mm512_ceil_ps(_mm512_sub_ps(t3,_1));
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
                      __m512d
		      weibull_discr_samp_zmm8r8(   const __m512d vrand,
		                                   const __m512d a,
						   const __m512d b) {

                         return (weibull_discr_icdf_zmm8r8(vrand,a,b));
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
                      __m512d 
                      gamma_zmm8r8(const __m512d x) {
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(64) static __m512d  c[7] = {
                                        _mm512_set1_pd(-1.910444077728e-03), 
                                        _mm512_set1_pd(8.4171387781295e-04), 
                                        _mm512_set1_pd(-5.952379913043012e-04), 
                                        _mm512_set1_pd(7.93650793500350248e-04),
                                        _mm512_set1_pd(-2.777777777777681622553e-03),
                                        _mm512_set1_pd(8.333333333333333331554247e-02),
                                        _mm512_set1_pd(5.7083835261e-03)};
                        __attribute__((section(".rodata")))
                        __ATTR_ALIGN__(64) static __m512d  p[8]  = {
                                        _mm512_set1_pd(-1.71618513886549492533811e+00),
                                        _mm512_set1_pd(2.47656508055759199108314e+01),
                                        _mm512_set1_pd(-3.79804256470945635097577e+02),
                                        _mm512_set1_pd(6.29331155312818442661052e+02),
                                        _mm512_set1_pd(8.66966202790413211295064e+02), 
                                        _mm512_set1_pd(-3.14512729688483675254357e+04), 
                                        _mm512_set1_pd(-3.61444134186911729807069e+04),
                                        _mm512_set1_pd(6.64561438202405440627855e+04)};
                       __attribute__((section(".rodata")))          
                       __ATTR_ALIGN__(64) static __m512d  q[8]   = {
                                        _mm512_set1_pd(-3.08402300119738975254353e+01), 
                                        _mm512_set1_pd(3.15350626979604161529144e+02),
                                        _mm512_set1_pd(-1.01515636749021914166146e+03),
                                        _mm512_set1_pd(-3.10777167157231109440444e+03),
                                        _mm512_set1_pd(2.25381184209801510330112e+04),
                                        _mm512_set1_pd(4.75584627752788110767815e+03),
                                        _mm512_set1_pd(-1.34659959864969306392456e+05),
                                        _mm512_set1_pd(1.15132259675553483497211e+05)};
                       __m512d pi                          = 
                                        _mm512_set1_pd(3.1415926535897932384626434e+00);
                       __m512d eps                         =
                                        _mm512_set1_pd(2.22e-16);
                       __m512d one                         =
                                        _mm512_set1_pd(1.0);
                       __m512d half                        =
                                        _mm512_set1_pd(0.5);
                       __m512d sqrtpi                      =
                                        _mm512_set1_pd(0.9189385332046727417803297e+00);
                       __m512d twelve                      = 
                                        _mm512_set1_pd(12.0);
                       __m512d two                         =
                                        _mm512_set1_pd(2.0);
                       __m512d xbig                        =
                                        _mm512_set1_pd(171.624);
                       __m512d xinf                        =
                                        _mm512_set1_pd(1.0e+30);
                       __m512d xminin                      =
                                        _mm512_set1_pd(2.23e-308);
                       
                       __m512d zero                        =
                                        _mm512_setzero_pd();
                       register __m512d res,sum,xden,xnum;
                       register __m512d y,y1,ysq,z,fact;
                       register __m512i n;
                       
                       bool     parity;
                       parity = false;
                       y      = x;
                       // Negative argument
                       if(_mm512_cmp_pd_mask(y,zero,_CMP_LE_OQ)) {
                          register __m512d t0,t1;
                          y  = negate_zmm8r8(x);
                          y1 = _mm512_castsi512_pd(_mm512_cvtt_roundpd_epu64(y,_MM_FROUND_NO_EXC));
                          res= _mm512_sub_pd(y,y1);
                          if(_mm512_cmp_pd_mask(res,zero,_CMP_NEQ_OQ)) {
                            
                             t0 = _mm512_mul_pd(_mm512_mul_pd(y1,half),two);
                             t1 = _mm512_castsi512_pd(_mm512_cvtt_roundpd_epu64(t0,_MM_FROUND_NO_EXC));
                             if(_mm512_cmp_pd_mask(y1,t1,_CMP_NEQ_OQ)) parity = true;
#if (USE_SLEEF_LIB) == 1
                             t0 = xsin(_mm512_mul_pd(pi,res));
#else
                             t0 = _mm512_sin_pd(_mm512_mul_pd(pi,res));
#endif                             
                             fact = _mm512_div_pd(negate_zmm8r8(pi),t0);
                             y    = _mm512_add_pd(y,one);
                          }
                          else {
                             res = xinf;
                             return (res);
                          }
                       }
                       // Positive argument
                       if(_mm512_cmp_pd_mask(y,eps,_CMP_LT_OQ)) {
                          __mmask8 m;
                          m = _mm512_cmp_pd_mask(xminin,y,_CMP_LE_OQ);
                          res = _mm512_mask_blend_pd(m,xinf,_mm512_div_pd(one,y));
                          return (res);
                       }
                  }
                  else if(_mm512_cmp_pd_mask(y,twelve,_CMP_LT_OQ)) {
                          y1 = y;
                          // 0.0 < argument < 1.0.
                          if(_mm512_cmp_pd_mask(y,one,_CMP_LT_OQ)) {
                             z = y;
                             y = _mm512_add_pd(y,one);
                          }
                          else {
                             //!  1.0 < argument < 12.0.
                             //!  Reduce argument if necessary.
                             n = _mm512_sub_epi64(mm512_castpd_si512(y),
                                                  _mm512_set1_epi64(1LL));
                             y = _mm512_sub_pd(y,_mm512_castsi512_pd(n));
                             z = _mm512_sub_pd(y,one);
                          }
                          //  Evaluate approximation for 1.0 < argument < 2.0.
                          xnum = zero;
                          xden = one;
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[0]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[0]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[1]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[1]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[2]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[2]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[3]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[3]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[4]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[4]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[5]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[5]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[6]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[7]);
                          xnum = _mm512_mul_pd(_mm512_add_pd(xnum,p[7]),z);
                          xden = _mm512_fmadd_pd(xden,z,q[7]);
                          res  = _mm512_add_pd(_mm512_div_pd(xnum,xden),one);
                          // Adjust result for case  0.0 < argument < 1.0.
                          if(_mm512_cmp_pd_mask(y1,y,_CMP_LT_OQ)) 
                             res = _mm512_div_pd(res,y1);
                          else if(_mm512_cmp_pd_mask(y,y1,_CMP_LT_OQ)) {
                          //  Important notice: ******For argument interval: 2.0<=x<12.0 there is a scalarization in form
                          //  of 8 scalar for-loops*******!!
                             __ATTR_ALIGN__(64) int64_t sn[8];
                             __ATTR_ALIGN__(64) double  sres[8];
                             __ATTR_ALIGN__(64) double  sy[8];
                             __ATTR_ALIGN__(64) double  sone[8];
                             int64_t i;
                             _mm512_store_si512(&sn[0],n);
                             _mm512_store_pd(&sres[0],res);
                             _mm512_store_pd(&sy[0],y);
                             _mm512_store_pd(&sone[0],one);
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
                             for(i=0; i != sn[4]; ++i) {
                                 sres[4] *= sy[4];
                                 sy[4]   += sone[4];
                             }
                             for(i=0; i != sn[5]; ++i) {
                                 sres[5] *= sy[5];
                                 sy[5]   += sone[5];
                             }
                             for(i=0; i != sn[6]; ++i) {
                                 sres[6] *= sy[6];
                                 sy[6]   += sone[6];
                             }
                             for(i=0; i != sn[7]; ++i) {
                                 sres[7] *= sy[7];
                                 sy[7]   += sone[7];
                             }
                             res = _mm512_load_pd(&sres[0]);
                             y   = _mm512_load_pd(&sy[0]);
                          }
                          
                       }
                       else {
                            //  Evaluate for 12.0 <= argument.
                            if(_mm512_cmp_pd_mask(y,xbig,_CMP_LE_OQ)) {
                               ysq = _mm512_mul_pd(y,y);
                               sum = c[6];
                               sum = _mm512_add_pd(_mm512_div_pd(sum,ysq),c[0]);
                               sum = _mm512_add_pd(_mm512_div_pd(sum,ysq),c[1]);
                               sum = _mm512_add_pd(_mm512_div_pd(sum,ysq),c[2]);
                               sum = _mm512_add_pd(_mm512_div_pd(sum,ysq),c[3]);
                               sum = _mm512_add_pd(_mm512_div_pd(sum,ysq),c[4]);
                               sum = _mm512_add_pd(_mm512_div_pd(sum,ysq),c[5]);
                               sum = _mm512_sub_pd(_mm512_div_pd(sum,y),
                                                   _mm512_add_pd(y,sqrtpi));
#if (USE_SLEEF_LIB) == 1
                               sum = _mm512_mul_pd(_mm512_add_pd(sum,
                                                         _mm512_sub_pd(y,half),xlog(y)));
                               res = xexp(sum);
#else
                               sum = _mm512_mul_pd(_mm512_add_pd(sum,
                                                         _mm512_sub_pd(y,half),_mm512_log_pd(y)));
                               res = _mm512_exp_pd(sum);
#endif

                                                    
                            }
                            else {
                               res = xinf;
                               return (res);
                            }
                       }
                       // !  Final adjustments and return.
                       if(parity) res = negate_zmm8r8(res);
                       if(_mm512_cmp_pd_mask(fact,one,_CMP_NEQ_OQ)) res = _mm512_div_pd(fact,res);
                       return (res);
                  }
                  
                  
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512 
                      gamma_zmm16r4(const __m512 x) {
                         __attribute__((section(".rodata")))
                         __ATTR_ALIGN__(64) static __m512  c[7] = {
                                        _mm512_set1_ps(-1.910444077728e-03f), 
                                        _mm512_set1_ps(8.4171387781295e-04f), 
                                        _mm512_set1_ps(-5.952379913043012e-04f), 
                                        _mm512_set1_ps(7.93650793500350248e-04f),
                                        _mm512_set1_ps(-2.777777777777681622553e-03f),
                                        _mm512_set1_ps(8.333333333333333331554247e-02f),
                                        _mm512_set1_ps(5.7083835261e-03f)};
                        __attribute__((section(".rodata")))
                        __ATTR_ALIGN__(64) static __m512  p[8]  = {
                                        _mm512_set1_ps(-1.71618513886549492533811e+00f),
                                        _mm512_set1_ps(2.47656508055759199108314e+01f),
                                        _mm512_set1_ps(-3.79804256470945635097577e+02f),
                                        _mm512_set1_ps(6.29331155312818442661052e+02f),
                                        _mm512_set1_ps(8.66966202790413211295064e+02f), 
                                        _mm512_set1_ps(-3.14512729688483675254357e+04f), 
                                        _mm512_set1_ps(-3.61444134186911729807069e+04f),
                                        _mm512_set1_ps(6.64561438202405440627855e+04f)}; 
                       __attribute__((section(".rodata")))         
                       __ATTR_ALIGN__(64) static __m512  q[8]   = {
                                        _mm512_set1_ps(-3.08402300119738975254353e+01f), 
                                        _mm512_set1_ps(3.15350626979604161529144e+02f),
                                        _mm512_set1_ps(-1.01515636749021914166146e+03f),
                                        _mm512_set1_ps(-3.10777167157231109440444e+03f),
                                        _mm512_set1_ps(2.25381184209801510330112e+04f),
                                        _mm512_set1_ps(4.75584627752788110767815e+03f),
                                        _mm512_set1_ps(-1.34659959864969306392456e+05f),
                                        _mm512_set1_ps(1.15132259675553483497211e+05f)};
                       __m512 pi                          = 
                                        _mm512_set1_ps(3.1415926535897932384626434e+00f);
                       __m512 eps                         =
                                        _mm512_set1_ps(2.22e-16f);
                       __m512 one                         =
                                        _mm512_set1_ps(1.0f);
                       __m512 half                        =
                                        _mm512_set1_ps(0.5f);
                       __m512 sqrtpi                      =
                                        _mm512_set1_ps(0.9189385332046727417803297e+00f);
                       __m512 twelve                      = 
                                        _mm512_set1_ps(12.0f);
                       __m512 two                         =
                                        _mm512_set1_ps(2.0f);
                       __m512 xbig                        =
                                        _mm512_set1_ps(171.624f);
                       __m512 xinf                        =
                                        _mm512_set1_ps(1.0e+30f);
                       __m512 xminin                      =
                                        _mm512_set1_ps(std::numeric_limits<float>::min());
                       
                       __m512 zero                        =
                                        _mm512_setzero_ps();
                       register __m512 res,sum,xden,xnum;
                       register __m512 y,y1,ysq,z,fact;
                       register __m512i n;
                       
                       bool     parity;
                       parity = false;
                       y      = x;
                       // Negative argument
                       if(_mm512_cmp_ps_mask(y,zero,_CMP_LE_OQ)) {
                          register __m512 t0,t1;
                          y  = negate_zmm16r4(x);
                          y1 = _mm512_castsi512_ps(_mm512_cvtt_roundps_epu32(y,_MM_FROUND_NO_EXC));
                          res= _mm512_sub_ps(y,y1);
                          if(_mm512_cmp_ps_mask(res,zero,_CMP_NEQ_OQ)) {
                            
                             t0 = _mm512_mul_ps(_mm512_mul_ps(y1,half),two);
                             t1 = _mm512_castsi512_ps(_mm512_cvtt_roundps_epu32(t0,_MM_FROUND_NO_EXC));
                             if(_mm512_cmp_ps_mask(y1,t1,_CMP_NEQ_OQ)) parity = true;
#if (USE_SLEEF_LIB) == 1
                             t0 = xsinf(_mm512_mul_ps(pi,res));
#else
                             t0 = _mm512_sin_ps(_mm512_mul_ps(pi,res));
#endif                             
                             fact = _mm512_div_ps(negate_zmm16r4(pi),t0);
                             y    = _mm512_add_ps(y,one);
                          }
                          else {
                             res = xinf;
                             return (res);
                          }
                       }
                       // Positive argument
                       if(_mm512_cmp_ps_mask(y,eps,_CMP_LT_OQ)) {
                          __mmask16 m;
                          m = _mm512_cmp_ps_mask(xminin,y,_CMP_LE_OQ);
                          res = _mm512_mask_blend_ps(m,xinf,_mm512_div_ps(one,y));
                          return (res);
                       }
                  }
                  else if(_mm512_cmp_ps_mask(y,twelve,_CMP_LT_OQ)) {
                          y1 = y;
                          // 0.0 < argument < 1.0.
                          if(_mm512_cmp_ps_mask(y,one,_CMP_LT_OQ)) {
                             z = y;
                             y = _mm512_add_ps(y,one);
                          }
                          else {
                             //!  1.0 < argument < 12.0.
                             //!  Reduce argument if necessary.
                             n = _mm512_sub_epi32(_mm512_castps_si512(y),
                                                  _mm512_set1_epi32(1));
                             y = _mm512_sub_ps(y,_mm512_castsi512_ps(n));
                             z = _mm512_sub_ps(y,one);
                          }
                          //  Evaluate approximation for 1.0 < argument < 2.0.
                          xnum = zero;
                          xden = one;
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[0]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[0]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[1]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[1]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[2]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[2]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[3]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[3]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[4]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[4]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[5]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[5]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[6]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[7]);
                          xnum = _mm512_mul_ps(_mm512_add_ps(xnum,p[7]),z);
                          xden = _mm512_fmadd_ps(xden,z,q[7]);
                          res  = _mm512_add_ps(_mm512_div_ps(xnum,xden),one);
                          // Adjust result for case  0.0 < argument < 1.0.
                          if(_mm512_cmp_ps_mask(y1,y,_CMP_LT_OQ)) 
                             res = _mm512_div_ps(res,y1);
                          else if(_mm512_cmp_ps_mask(y,y1,_CMP_LT_OQ)) {
                          //  Important notice: ******For argument interval: 2.0<=x<12.0 there is a scalarization in form
                          //  of 8 scalar for-loops*******!!
                             __ATTR_ALIGN__(64) int32_t sn[16];
                             __ATTR_ALIGN__(64) float  sres[16];
                             __ATTR_ALIGN__(64) float  sy[16];
                             __ATTR_ALIGN__(64) float  sone[16];
                             int32_t i;
                             _mm512_store_si512(&sn[0],n);
                             _mm512_store_ps(&sres[0],res);
                             _mm512_store_ps(&sy[0],y);
                             _mm512_store_ps(&sone[0],one);
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
                             for(i=0; i != sn[4]; ++i) {
                                 sres[4] *= sy[4];
                                 sy[4]   += sone[4];
                             }
                             for(i=0; i != sn[5]; ++i) {
                                 sres[5] *= sy[5];
                                 sy[5]   += sone[5];
                             }
                             for(i=0; i != sn[6]; ++i) {
                                 sres[6] *= sy[6];
                                 sy[6]   += sone[6];
                             }
                             for(i=0; i != sn[7]; ++i) {
                                 sres[7] *= sy[7];
                                 sy[7]   += sone[7];
                             }
                              for(i=0; i != sn[8]; ++i) {
                                 sres[8] *= sy[8];
                                 sy[8]   += sone[8];
                             }
                             for(i=0; i != sn[9]; ++i) {
                                 sres[9] *= sy[9];
                                 sy[9]   += sone[9];
                             }
                             for(i=0; i != sn[10]; ++i) {
                                 sres[10] *= sy[10];
                                 sy[10]   += sone[10];
                             }
                             for(i=0; i != sn[11]; ++i) {
                                 sres[11] *= sy[11];
                                 sy[11]   += sone[11];
                             }
                             for(i=0; i != sn[12]; ++i) {
                                 sres[12] *= sy[12];
                                 sy[12]   += sone[12];
                             }
                             for(i=0; i != sn[13]; ++i) {
                                 sres[13] *= sy[13];
                                 sy[13]   += sone[13];
                             }
                             for(i=0; i != sn[14]; ++i) {
                                 sres[14] *= sy[14];
                                 sy[14]   += sone[614];
                             }
                             for(i=0; i != sn[15]; ++i) {
                                 sres[15] *= sy[15];
                                 sy[15]   += sone[15];
                             }
                             res = _mm512_load_ps(&sres[0]);
                             y   = _mm512_load_ps(&sy[0]);
                          }
                          
                       }
                       else {
                            //  Evaluate for 12.0 <= argument.
                            if(_mm512_cmp_ps_mask(y,xbig,_CMP_LE_OQ)) {
                               ysq = _mm512_mul_ps(y,y);
                               sum = c[6];
                               sum = _mm512_add_ps(_mm512_div_ps(sum,ysq),c[0]);
                               sum = _mm512_add_ps(_mm512_div_ps(sum,ysq),c[1]);
                               sum = _mm512_add_ps(_mm512_div_ps(sum,ysq),c[2]);
                               sum = _mm512_add_ps(_mm512_div_ps(sum,ysq),c[3]);
                               sum = _mm512_add_ps(_mm512_div_ps(sum,ysq),c[4]);
                               sum = _mm512_add_ps(_mm512_div_ps(sum,ysq),c[5]);
                               sum = _mm512_sub_ps(_mm512_div_ps(sum,y),
                                                   _mm512_add_ps(y,sqrtpi));
#if (USE_SLEEF_LIB) == 1
                               sum = _mm512_mul_ps(_mm512_add_ps(sum,
                                                         _mm512_sub_ps(y,half),xlogf(y)));
                               res = xexpf(sum);
#else
                               sum = _mm512_mul_ps(_mm512_add_ps(sum,
                                                         _mm512_sub_ps(y,half),_mm512_log_ps(y)));
                               res = _mm512_exp_ps(sum);
#endif

                                                    
                            }
                            else {
                               res = xinf;
                               return (res);
                            }
                       }
                       // !  Final adjustments and return.
                       if(parity) res = negate_zmm16r4(res);
                       if(_mm512_cmp_ps_mask(fact,one,_CMP_NEQ_OQ)) res = _mm512_div_ps(fact,res);
                       return (res);
                  }
                  
                  
                  
/*
!*****************************************************************************80
!
!! STUDENT_PDF evaluates the central Student T PDF.
!
!  Discussion:
!
!    PDF(A,B,C;X) = Gamma ( (C+1)/2 ) /
!      ( Gamma ( C / 2 ) * Sqrt ( PI * C )
!      * ( 1 + ((X-A)/B)^2/C )^(C + 1/2 ) )
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
!    Input, real ( kind = 8 ) X, the argument of the PDF.
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
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
*/


          
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512d 
                      student_pdf_zmm8r8(const __m512d x,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d c) {
                         
                         __m512d C314159265358979323846264 = 
                                           _mm512_set1_pd(3.14159265358979323846264);   
                         __m512d C1 = _mm512_set1_pd(1.0);
                         __m512d C05= _mm512_set1_pd(0.5);  
                         __m512d C2 = _mm512_set1_pd(2.0); 
                         register __m512d y,t0,t1,t2,t3;
                         register __m512d r8g1,r8g2,pdf;
                         t0   = _mm512_mul_pd(C05,_mm512_add_pd(c,C1));
                         y    = _mm512_div_pd(_mm512_sub_pd(x,a),b);
                         r8g1 = gamma_zmm8r8(t0);
                         t1   = _mm512_fmadd_pd(C2,c,C1);
                         t2   = _mm512_add_pd(_mm512_div_pd(_mm512_mul_pd(y,y),c),C1);
                         t0   = _mm512_pow_pd(t2,t1); //used
                         r8g2 = gamma_zmm8r8(_mm512_mul_pd(C05,c));
                         y    = _mm512_sqrt_pd(_mm512_mul_pd(C314159265358979323846264,c));
                         t3   = _mm512_mul_pd(y,_mm512_mul_pd(r8g2,t0));
                         pdf  = _mm512_div_pd(r8g1,t3);
                         return (pdf);
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
                      __m512d 
                      student_variance_zmm8r8(const __m512d a,
                                              const __m512d b,
                                              const __m512d c) {
                                              
                          __m512d C2 = _mm512_set1_pd(2.0);     
                          register __m512d bb,t1;
                          register __m512d var;
                          bb = _mm512_mul_pd(b,b);
                          t1 = _mm512_sub_pd(c,C2);
                          var= _mm512_mul_pd(bb,_mm512_div_pd(c,t1));
                          return (var);                    
                    }   
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
                      student_variance_zmm16r4(const __m512 a,
                                              const __m512 b,
                                              const __m512 c) {
                                              
                          __m512 C2 = _mm512_set1_ps(2.0f);     
                          register __m512 bb,t1;
                          register __m512 var;
                          bb = _mm512_mul_ps(b,b);
                          t1 = _mm512_sub_ps(c,C2);
                          var= _mm512_mul_ps(bb,_mm512_div_ps(c,t1));
                          return (var);                    
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
                      __m512d
                      trigamma_zmm8r8(const __m512d x) {
                         
                         __m512d a  = _mm512_setzero_pd();
                         __m512d C1 = _mm512_set1_pd(1.0);
                         __m512d C05=_mm512_set1_pd(0.5);
                         __m512d b  = _mm512_set1_pd(5.0);
                         __m512d b2 = _mm512_set1_pd(1.0/6.0);
                         __m512d b4 = _mm512_set1_pd(-1.0/30.0);
                         __m512d b6 = _mm512_set1_pd(1.0/42.0);
                         __m512d b8 = _mm512_set1_pd(-1.0/30.0);
                         register __m512d y,z,t0,t1;
                         register __m512d trig;
                         
                         if(_mm512_cmp_pd_mask(x,a,_CMP_LE_OQ)) {
                            trig = _mm512_div_pd(C1,_mm512_mul_pd(x,x));
                         }
                         else {
                            z = x;
                            trig = a;
                            while(_mm512_cmp_pd_mask(z,b,_CMP_LT_OQ)) {
                                  trig = _mm512_add_pd(_mm512_div_pd(C1,
                                                         _mm512_mul_pd(z,z)))
                                  z    = _mm512_add_pd(z,C1);
                            }
                            y    = _mm512_div_pd(C1,_mm512_mul_pd(z,z));
                            trig = trig+C05*y+(C1+y*(b2+y*(b4+y*(b6+y*b8))))/z; 
                         }
                         return (trig);
                    } 
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
                      trigamma_zmm16r4(const __m512 x) {
                         
                         __m512 a  = _mm512_setzero_ps();
                         __m512 C1 = _mm512_set1_ps(1.0f);
                         __m512 C05=_mm512_set1_ps(0.5f);
                         __m512 b  = _mm512_set1_ps(5.0f);
                         __m512 b2 = _mm512_set1_ps(1.0f/6.0f);
                         __m512 b4 = _mm512_set1_ps(-1.0f/30.0f);
                         __m512 b6 = _mm512_set1_ps(1.0f/42.0f);
                         __m512 b8 = _mm512_set1_ps(-1.0f/30.0f);
                         register __m512 y,z,t0,t1;
                         register __m512 trig;
                         
                         if(_mm512_cmp_ps_mask(x,a,_CMP_LE_OQ)) {
                            trig = _mm512_div_ps(C1,_mm512_mul_ps(x,x));
                         }
                         else {
                            z = x;
                            trig = a;
                            while(_mm512_cmp_ps_mask(z,b,_CMP_LT_OQ)) {
                                  trig = _mm512_add_ps(_mm512_div_ps(C1,
                                                         _mm512_mul_ps(z,z)))
                                  z    = _mm512_add_ps(z,C1);
                            }
                            y    = _mm512_div_ps(C1,_mm512_mul_ps(z,z));
                            trig = trig+C05*y+(C1+y*(b2+y*(b4+y*(b6+y*b8))))/z; 
                         }
                         return (trig);
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
                      __m512d	 
                      weibull_pdf_zmm8r8(const __m512d x,
                                         const __m512d a,
                                         const __m512d b,
                                         const __m512d c) {
                        
                         register __m512d C1 = _mm512_set1_pd(1.0);
                         register __m512d y,t0,pow1,t1,exp;
                         register __m512d pdf;
                         t0 = _mm512_div_pd(_mm512_sub_pd(x,a),b);
                         pow1 = _mm512_pow_pd(t0,_mm512_sub_pd(c,C1));
#if (USE_SLEEF_LIB) == 1
                         exp  = xexp(_mm512_pow_pd(y,c));
#else
                         exp  = _mm512_exp_pd(_mm512_pow_pd(y,c));
#endif                   
                         t1   = _mm512_div_pd(c,b);
                         pdf  = _mm512_div_pd(_mm512_mul_pd(t1,pow1),exp); 
                         return (pdf);     
                   }
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512	 
                      weibull_pdf_zmm16r4(const __m512 x,
                                         const __m512 a,
                                         const __m512 b,
                                         const __m512 c) {
                        
                         register __m512 C1 = _mm512_set1_ps(1.0);
                         register __m512 y,t0,pow1,t1,exp;
                         register __m512 pdf;
                         t0 = _mm512_div_ps(_mm512_sub_ps(x,a),b);
                         pow1 = _mm512_pow_ps(t0,_mm512_sub_ps(c,C1));
#if (USE_SLEEF_LIB) == 1
                         exp  = xexp(_mm512_pow_ps(y,c));
#else
                         exp  = _mm512_exp_ps(_mm512_pow_ps(y,c));
#endif                   
                         t1   = _mm512_div_ps(c,b);
                         pdf  = _mm512_div_ps(_mm512_mul_ps(t1,pow1),exp); 
                         return (pdf);     
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

#include "GMS_rotations_avx512_helpers.hpp" // fmod


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512d
		      von_misses_cdf_zmm8r8(const __m512d x,
		                            const __m512d a,
					    const __m512d b) {

                        //Early exit.
			const __m512d   _0  = _mm512_setzero_pd();
			const __m512d   _1  = _mm512_set1_pd(1.0);
			const __m512d   pi  = _mm512_set1_pd(3.14159265358979323846264338328);
			const __m512d   npi = _mm512_set1_pd(-3.14159265358979323846264338328);
			const __m512d   xsa = _mm512_sub_pd(x,a);
			if(__builtin_expect(_mm512_cmp_pd_mask(xsa,npi,_CMP_LE_OQ),0)) {
		             return (_0);
			}
			if(__builtin_expect(_mm512_cmp_pd_mask(npi,xsa,_CMP_LE_OQ),0)) {
                             return (_1); 
			}
			const __m512d  _2pi = _mm512_set1_pd(6.283185307179586476925286766559);
			const __m512d  a1  = _mm512_set1_pd(12.0);
			const __m512d  a2  = _mm512_set1_pd(0.8);
			const __m512d  a3  = _mm512_set1_pd(8.0);
			const __m512d  a4  = _mm512_set1_pd(1.0);
			const __m512d  c1  = _mm512_set1_pd(56.0);
			const __m512d  ck  = _mm512_set1_pd(10.5);
			const __m512d  _2  = _mm512_set1_pd(2.0);
			const __m512d  _1_2= _mm512_set1_pd(0.5);
			__m512d arg,cdf,cn,p,r,s,sn,u,v,y,z,uprv,erfx;
			//  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
			z    = b;
			u    = fmod_zmm8r8(_mm512_add_pd(xsa,pi),_2pi);
			uprv = u;
			const __mmask8 m = _mm512_cmp_pd_mask(u,_0,_CMP_LT_OQ);
			u    = _mm512_add_pd(u,_2pi);
			u    = _mm512_mask_blend_pd(m,uprv,u);
			y    = _mm512_sub_pd(u,pi);
			
			//For small B, sum IP terms by backwards recursion.
			// Can not be vectorized manually, hence 0 is returned.
			// Only large B is computed.
			/*
                              This scalar code can not be vectorized.
                              ip = int ( z * a2 - a3 / ( z + a4 ) + a1 )
                              Used as loop control variable
                              do n = 2, ip
                         */
                        if(_mm512_cmp_pd_mask(z,ck,_CMP_LE_OQ)) {
                           return (_0);
			}
			else {
                           const __m512d t0 = _mm512_set1_pd(24.0);
			   const __m512d t1 = _mm512_set1_pd(54.0);
			   const __m512d t2 = _mm512_set1_pd(347.0);
			   const __m512d t3 = _mm512_set1_pd(26.0);
			   const __m512d t4 = _mm512_set1_pd(6.0);
			   const __m512d t5 = _mm512_set1_pd(12.0);
			   const __m512d t6 = _mm512_set1_pd(3.0);
			   const __m512d t7 = _mm512_set1_pd(16.0);
			   const __m512d t8 = _mm512_set1_pd(1.75);
			   const __m512d t9 = _mm512_set1_pd(83.5);
			   c                = _mm512_mul_pd(t0,z);
			   v                = _mm512_sub_pd(c,c1);
			   const __m512d tmp1 = _mm512_sub_pd(_mm512_add_pd(v,t3),c);
			   const __m512d tmp2 = _mm512_div_pd(t1,_mm512_div_pd(t2,tmp1));
			   const __m512d tmp3 = _mm512_add_pd(_mm512_sub_pd(tmp2,t4),c);
			   r                  = _mm512_sqrt_pd(_mm512_div_pd(tmp3,t5));
#if (USE_SLEEF_LIB) == 1
                           z                  = _mm512_mul_pd(xsin(_mm512_mul_pd(_1_2,y)),r);
#else
			   z                  = _mm512_mul_pd(_mm512_sin_pd(
			                                                _mm512_mul_pd(_1_2,y)),r);
#endif
                           s                  = _mm512_mul_pd(_2,_mm512_mul_pd(z,z));
			   v                  = _mm512_sub_pd(v,_mm512_add_pd(s,t6));
			   y                  = _mm512_div_pd(_mm512_sub_pd(_mm512_sub_pd(c,s),
			                                                    _mm512_sub_pd(s,t7)),t6);
			   tmp1               = _mm512_sub_pd(v,y);
			   y                  = _mm512_div_pd(_mm512_fmadd_pd(_mm512_add_pd(s,t8),s,t9),tmp1);
			   tmp2               = _mm512_mul_pd(y,y);
			   arg                = _mm512_mul_pd(z,_mm512_sub_pd(_1,
			                                                  _mm512_div_pd(s,tmp2)));
			   erfx               = _mm512_erf_pd(arg);
			   cdf                = _mm512_fmadd_pd(_1_2,erfx,_1_2);
			}
			cdf                   = _mm512_max_pd(cdf,_0);
			cdf                   = _mm512_min_pd(cdf,_1);
			return (cdf);
			
		   }


		      
	              __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      von_misses_cdf_zmm16r4(const __m512 x,
		                            const __m512 a,
					    const __m512 b) {

                        //Early exit.
			const __m512   _0  = _mm512_setzero_ps();
			const __m512   _1  = _mm512_set1_ps(1.0f);
			const __m512   pi  = _mm512_set1_ps(3.14159265358979323846264338328f);
			const __m512   npi = _mm512_set1_ps(-3.14159265358979323846264338328f);
			const __m512   xsa = _mm512_sub_ps(x,a);
			if(__builtin_expect(_mm512_cmp_ps_mask(xsa,npi,_CMP_LE_OQ),0)) {
		             return (_0);
			}
			if(__builtin_expect(_mm512_cmp_ps_mask(npi,xsa,_CMP_LE_OQ),0)) {
                             return (_1); 
			}
			const __m512  _2pi = _mm512_set1_ps(6.283185307179586476925286766559f);
			const __m512  a1  = _mm512_set1_ps(12.0f);
			const __m512  a2  = _mm512_set1_ps(0.8f);
			const __m512  a3  = _mm512_set1_ps(8.0f);
			const __m512  a4  = _mm512_set1_ps(1.0f);
			const __m512  c1  = _mm512_set1_ps(56.0f);
			const __m512  ck  = _mm512_set1_ps(10.5f);
			const __m512  _2  = _mm512_set1_ps(2.0f);
			const __m512  _1_2= _mm512_set1_ps(0.5f);
			__m512 arg,cdf,cn,p,r,s,sn,u,v,y,z,uprv,erfx;
			//  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
			z    = b;
			u    = fmod_zmm16r4(_mm512_add_ps(xsa,pi),_2pi);
			uprv = u;
			const __mmask16 m = _mm512_cmp_ps_mask(u,_0,_CMP_LT_OQ);
			u    = _mm512_add_ps(u,_2pi);
			u    = _mm512_mask_blend_ps(m,uprv,u);
			y    = _mm512_sub_ps(u,pi);
			
			//For small B, sum IP terms by backwards recursion.
			// Can not be vectorized manually, hence 0 is returned.
			// Only large B is computed.
			/*
                              This scalar code can not be vectorized.
                              ip = int ( z * a2 - a3 / ( z + a4 ) + a1 )
                              Used as loop control variable
                              do n = 2, ip
                         */
                        if(_mm512_cmp_ps_mask(z,ck,_CMP_LE_OQ)) {
                           return (_0);
			}
			else {
                           const __m512 t0 = _mm512_set1_ps(24.0f);
			   const __m512 t1 = _mm512_set1_ps(54.0f);
			   const __m512 t2 = _mm512_set1_ps(347.0f);
			   const __m512 t3 = _mm512_set1_ps(26.0f);
			   const __m512 t4 = _mm512_set1_ps(6.0f);
			   const __m512 t5 = _mm512_set1_ps(12.0f);
			   const __m512 t6 = _mm512_set1_ps(3.0f);
			   const __m512 t7 = _mm512_set1_ps(16.0f);
			   const __m512 t8 = _mm512_set1_ps(1.75f);
			   const __m512 t9 = _mm512_set1_ps(83.5f);
			   c                = _mm512_mul_ps(t0,z);
			   v                = _mm512_sub_ps(c,c1);
			   const __m512d tmp1 = _mm512_sub_ps(_mm512_add_ps(v,t3),c);
			   const __m512d tmp2 = _mm512_div_ps(t1,_mm512_div_ps(t2,tmp1));
			   const __m512d tmp3 = _mm512_add_ps(_mm512_sub_ps(tmp2,t4),c);
			   r                  = _mm512_sqrt_ps(_mm512_div_ps(tmp3,t5));
#if (USE_SLEEF_LIB) == 1
                           z                  = _mm512_mul_ps(xsinf(_mm512_mul_ps(_1_2,y)),r);
#else
			   z                  = _mm512_mul_ps(_mm512_sin_ps(
			                                                _mm512_mul_ps(_1_2,y)),r);
#endif
                           s                  = _mm512_mul_ps(_2,_mm512_mul_ps(z,z));
			   v                  = _mm512_sub_ps(v,_mm512_add_ps(s,t6));
			   y                  = _mm512_div_ps(_mm512_sub_ps(_mm512_sub_ps(c,s),
			                                                    _mm512_sub_ps(s,t7)),t6);
			   tmp1               = _mm512_sub_ps(v,y);
			   y                  = _mm512_div_ps(_mm512_fmadd_ps(_mm512_add_ps(s,t8),s,t9),tmp1);
			   tmp2               = _mm512_mul_ps(y,y);
			   arg                = _mm512_mul_ps(z,_mm512_sub_ps(_1,
			                                                  _mm512_div_ps(s,tmp2)));
			   erfx               = _mm512_erf_ps(arg);
			   cdf                = _mm512_fmadd_ps(_1_2,erfx,_1_2);
			}
			cdf                   = _mm512_max_ps(cdf,_0);
			cdf                   = _mm512_min_ps(cdf,_1);
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
                      __m512d
		      von_misses_pdf_zmm8r8(const __m512d x,
		                            const __m512d a,
					    const __m512d b) {
 
                           const __m512d   pi  = _mm512_set1_pd(3.14159265358979323846264338328);
			   const __m512d   _2pi= _mm512_set1_pd(6.283185307179586476925286766559);
			   const __m512d   _0  = _mm512_setzero_pd();
			   const __m512d   _2  = _mm512_set1_pd(2.0);
			   const __m512d   t0  = _mm512_sub_pd(a,pi);
			   const __m512d   t1  = _mm512_add_pd(a,pi);
			   __m512d pdf;
			   __mmask8 m1,m2;
			   m1                  = _mm512_cmp_pd_mask(x,t0,_CMP_LT_OQ);
			   pdf                 = _mm512_mask_blend_pd(m1,_0,_0);
			   m2                  = _mm512_cmp_pd_mask(x,t1,_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1
			   const __m512d tmp1  = xexp(_mm512_mul_pd(b,xcos(_mm512_sub_pd(x,a))));
#else
                           const __m512d tmp1  = _mm512_exp(_mm512_mul_pd(b,
			                                              _mm512_cos_pd(
								                _mm512_sub_pd(x,a))));
#endif
                           
			   pdf                 = _mm512_mask_blend_pd(m2,_0,_mm512_div_pd(tmp1,
			                                              _mm512_mul_pd(_2pi,bessesl_i0_zmm8r8(b))));
			   return (pdf);
		   }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      von_misses_pdf_zmm16r4(const __m512 x,
		                            const __m512 a,
					    const __m512 b) {
 
                           const __m512   pi  = _mm512_set1_pd(3.14159265358979323846264338328f);
			   const __m512   _2pi= _mm512_set1_pd(6.283185307179586476925286766559f);
			   const __m512   _0  = _mm512_setzero_pd();
			   const __m512   _2  = _mm512_set1_pd(2.0);
			   const __m512   t0  = _mm512_sub_pd(a,pi);
			   const __m512   t1  = _mm512_add_pd(a,pi);
			   __m512 pdf;
			   __mmask16 m1,m2;
			   m1                  = _mm512_cmp_pd_mask(x,t0,_CMP_LT_OQ);
			   pdf                 = _mm512_mask_blend_pd(m1,_0,_0);
			   m2                  = _mm512_cmp_pd_mask(x,t1,_CMP_LE_OQ);
#if (USE_SLEEF_LIB) == 1
			   const __m512 tmp1  = xexpf(_mm512_mul_pd(b,xcosf(_mm512_sub_pd(x,a))));
#else
                           const __m512 tmp1  = _mm512_exp(_mm512_mul_pd(b,
			                                              _mm512_cos_pd(
								                _mm512_sub_pd(x,a))));
#endif
                           
			   pdf                 = _mm512_mask_blend_pd(m2,_0,_mm512_div_pd(tmp1,
			                                              _mm512_mul_pd(_2pi,bessesl_i0_zmm8r8(b))));
			   return (pdf);
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
		      __m512d owen_tfunc_zmm8r8(const __m512d h,
		                                const __m512d a) {
		             
		             __attribute__((section(".rodata")))
		             __ATTR_ALIGN__(64) static __m512d  weight[10] = {
		                                _mm512_set1_pd(0.666713443086881375935688098933e-01),
                                                _mm512_set1_pd(0.149451349150580593145776339658e+00),
                                                _mm512_set1_pd(0.219086362515982043995534934228e+00),
                                                _mm512_set1_pd(0.269266719309996355091226921569e+00),
                                                _mm512_set1_pd(0.295524224714752870173892994651e+00),
                                                _mm512_set1_pd(0.295524224714752870173892994651e+00),
                                                _mm512_set1_pd(0.269266719309996355091226921569e+00),
                                                _mm512_set1_pd(0.219086362515982043995534934228e+00),
                                                _mm512_set1_pd(0.149451349150580593145776339658e+00), 
                                                _mm512_set1_pd(0.666713443086881375935688098933e-01)};
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(64) static __m512d  xtab[10] = {
		                                _mm512_set1_pd(-0.973906528517171720077964012084e+00),
                                                _mm512_set1_pd(-0.865063366688984510732096688423e+00),
                                                _mm512_set1_pd(-0.679409568299024406234327365115e+00), 
                                                _mm512_set1_pd(-0.433395394129247190799265943166e+00), 
                                                _mm512_set1_pd(-0.148874338981631210884826001130e+00), 
                                                _mm512_set1_pd(0.148874338981631210884826001130e+00), 
                                                _mm512_set1_pd(0.433395394129247190799265943166e+00), 
                                                _mm512_set1_pd(0.679409568299024406234327365115e+00), 
                                                _mm512_set1_pd(0.865063366688984510732096688423e+00), 
                                                _mm512_set1_pd(0.973906528517171720077964012084e+00)};
		           
		             const __m512d twopinv = _mm512_set1_pd(0.15915494309189533576888e+00);   
		             const __m512d tv1     = _mm512_set1_pd(1.0e-35);
		             const __m512d tv2     = _mm512_set1_pd(15.0);
		             const __m512d tv3     = tv2;
		             const __m512d tv4     = _mm512_set1_pd(1.0e-5);
		             const __m512d C05     = _mm512_set1_pd(0.5);
		             const __m512d C1      = _mm512_set1_pd(1.0);
		             const __m512d C2      = _mm512_set1_pd(2.0);
		             const __m512d C025    = _mm512_set1_pd(0.25);
		             const __m512d C0      = _mm512_setzero_pd();
		             __m512d x,rt,as,h1,h2,hs,t0,t1,t2;
		             __m512d tfn;
		             if(_mm512_cmp_pd_mask(_mm512_abs_pd(h),tv1,_CMP_LT_OQ)) {
#if (USE_SLEEF_LIB) == 1
                                tfn = _mm512_mul_pd(xatan(a),twopinv);	
#else
                                tfn = _mm512_mul_pd(_mm512_atan_pd(a),twopinv);
#endif	               
		             }
		             else if(_mm512_cmp_pd_mask(tv2,_mm512_abs_pd(h),_CMP_LT_OQ)) {
		                tfn = C0;
		             }
		             else if(_mm512_cmp_pd_mask(_mm512_abs_pd(a),tv1,_CMP_LT_OQ)) {
		                 tfn = C0;
		             }
		             else {
		                 hs = _mm512_mul_pd(negate_zmm8r8(C05),
		                            _mm512_mul_pd(h,h));
		                 h2 = a;
		                 as = _mm512_mul_pd(a,a);
#if (USE_SLEEF_LIB) == 1		                 
		                 t0 = xlog(_mm512_add_pd(C1,as));
#else
                                 t0 = _mm512_log_pd(_mm512_add_pd(C1,as));
#endif		           
                                 __mmask8 m = _mm512_cmp_pd_mask(tv3,_mm512_sub_pd(t0,
                                                                     _mm512_mul_pd(hs,as)),_CMP_LE_OQ);
                                 if(m) {
                                    h1 = _mm512_mul_pd(C05,a);
                                    as = _mm512_mul_pd(C025,as);
                                    while(true) {
                                          rt = _mm512_add_pd(as,C1);
#if (USE_SLEEF_LIB) == 1
                                          t0 = _mm512_add_pd(h1,_mm512_fmadd_pd(hs,as,xlog(rt)));
#else
                                          t0 = _mm512_add_pd(h1,_mm512_fmadd_pd(hs,as,_mm512_log_pd(rt)));
#endif                                   
                                          t1 = _mm512_sub_pd(_mm512_div_pd(C1,rt),hs);
                                          t2 = _mm512_mul_pd(C2,_mm512_mul_pd(h1,t1));
                                          h2 = _mm512_div_pd(t0,t2);
                                          as = _mm512_mul_pd(h2,h2);
                                          if(_mm512_cmp_pd_mask(_mm512_abs_pd(
                                                          _mm512_mul_pd(h2,h1),tv4,_CMP_LT_OQ))) break;
                                          h1 = h2;                                               
                                    }
                                 }
                                 rt = C0;
                                 
                                 //for(int32_t i=0; i<10; ++i) 
#if (USE_SLEEF_LIB) == 1
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[0],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[0],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[1],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[1],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[2],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[2],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[3],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[3],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[4],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[4],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[5],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[5],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[6],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[6],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[7],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[7],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[8],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[8],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[9],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[9],
                                               _mm512_div_pd(xexp(_mm512_mul_pd(hs,t0)),t0)));
                                 
#else
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[0],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[0],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[1],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[1],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[2],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[2],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[3],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[3],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[4],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[4],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[5],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[5],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[6],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[6],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[7],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[7],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[8],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[8],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0)));
                                 x  = _mm512_fmadd_pd(C05,h2,_mm512_add_pd(xtab[9],C1));
                                 t0 = _mm512_add_pd(C1,_mm512_mul_pd(x,x));
                                 rt = _mm512_add_pd(rt,_mm512_mul_pd(weight[9],
                                               _mm512_div_pd(_mm512_exp_pd(_mm512_mul_pd(hs,t0)),t0))); 

#endif
                                 t1 = _mm512_mul_pd(C05,h2);
                                 tfn= _mm512_mul_pd(rt,_mm512_mul_pd(t1,twopinv)); 
		             }
		             return (tfn);
		    }   
		    
		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline  
		      __m512 owen_tfunc_zmm16r4(const __m512 h,
		                                const __m512 a) {
		             
		             __attribute__((section(".rodata")))
		             __ATTR_ALIGN__(64) static __m512  weight[10] = {
		                                _mm512_set1_ps(0.666713443086881375935688098933e-01f),
                                                _mm512_set1_ps(0.149451349150580593145776339658e+00f),
                                                _mm512_set1_ps(0.219086362515982043995534934228e+00f),
                                                _mm512_set1_ps(0.269266719309996355091226921569e+00f),
                                                _mm512_set1_ps(0.295524224714752870173892994651e+00f),
                                                _mm512_set1_ps(0.295524224714752870173892994651e+00f),
                                                _mm512_set1_ps(0.269266719309996355091226921569e+00f),
                                                _mm512_set1_ps(0.219086362515982043995534934228e+00f),
                                                _mm512_set1_ps(0.149451349150580593145776339658e+00f), 
                                                _mm512_set1_ps(0.666713443086881375935688098933e-01f)};
                            __attribute__((section(".rodata")))
		            __ATTR_ALIGN__(64) static __m512  xtab[10] = {
		                                _mm512_set1_ps(-0.973906528517171720077964012084e+00f),
                                                _mm512_set1_ps(-0.865063366688984510732096688423e+00f),
                                                _mm512_set1_ps(-0.679409568299024406234327365115e+00f), 
                                                _mm512_set1_ps(-0.433395394129247190799265943166e+00f), 
                                                _mm512_set1_ps(-0.148874338981631210884826001130e+00f), 
                                                _mm512_set1_ps(0.148874338981631210884826001130e+00f), 
                                                _mm512_set1_ps(0.433395394129247190799265943166e+00f), 
                                                _mm512_set1_ps(0.679409568299024406234327365115e+00f), 
                                                _mm512_set1_ps(0.865063366688984510732096688423e+00f), 
                                                _mm512_set1_ps(0.973906528517171720077964012084e+00f)};
		           
		             const __m512 twopinv = _mm512_set1_ps(0.15915494309189533576888e+00f);   
		             const __m512 tv1     = _mm512_set1_ps(std::numeric_limits<float>::min());
		             const __m512 tv2     = _mm512_set1_ps(15.0f);
		             const __m512 tv3     = tv2;
		             const __m512 tv4     = _mm512_set1_ps(1.0e-5f);
		             const __m512 C05     = _mm512_set1_ps(0.5f);
		             const __m512 C1      = _mm512_set1_ps(1.0f);
		             const __m512 C2      = _mm512_set1_ps(2.0f);
		             const __m512 C025    = _mm512_set1_ps(0.25f);
		             const __m512 C0      = _mm512_setzero_ps();
		             __m512 x,rt,as,h1,h2,hs,t0,t1,t2;
		             __m512 tfn;
		             if(_mm512_cmp_ps_mask(_mm512_abs_ps(h),tv1,_CMP_LT_OQ)) {
#if (USE_SLEEF_LIB) == 1
                                tfn = _mm512_mul_ps(xatanf(a),twopinv);	
#else
                                tfn = _mm512_mul_ps(_mm512_atan_ps(a),twopinv);
#endif	               
		             }
		             else if(_mm512_cmp_ps_mask(tv2,_mm512_abs_ps(h),_CMP_LT_OQ)) {
		                tfn = C0;
		             }
		             else if(_mm512_cmp_pd_mask(_mm512_abs_ps(a),tv1,_CMP_LT_OQ)) {
		                 tfn = C0;
		             }
		             else {
		                 hs = _mm512_mul_ps(negate_zmm16r4(C05),
		                            _mm512_mul_ps(h,h));
		                 h2 = a;
		                 as = _mm512_mul_ps(a,a);
#if (USE_SLEEF_LIB) == 1		                 
		                 t0 = xlogf(_mm512_add_ps(C1,as));
#else
                                 t0 = _mm512_log_ps(_mm512_add_ps(C1,as));
#endif		           
                                 __mmask16 m = _mm512_cmp_ps_mask(tv3,_mm512_sub_ps(t0,
                                                                     _mm512_mul_ps(hs,as)),_CMP_LE_OQ);
                                 if(m) {
                                    h1 = _mm512_mul_ps(C05,a);
                                    as = _mm512_mul_ps(C025,as);
                                    while(true) {
                                          rt = _mm512_add_ps(as,C1);
#if (USE_SLEEF_LIB) == 1
                                          t0 = _mm512_add_ps(h1,_mm512_fmadd_ps(hs,as,xlogf(rt)));
#else
                                          t0 = _mm512_add_ps(h1,_mm512_fmadd_ps(hs,as,_mm512_log_ps(rt)));
#endif                                   
                                          t1 = _mm512_sub_ps(_mm512_div_ps(C1,rt),hs);
                                          t2 = _mm512_mul_ps(C2,_mm512_mul_ps(h1,t1));
                                          h2 = _mm512_div_ps(t0,t2);
                                          as = _mm512_mul_ps(h2,h2);
                                          if(_mm512_cmp_ps_mask(_mm512_abs_ps(
                                                          _mm512_mul_ps(h2,h1),tv4,_CMP_LT_OQ))) break;
                                          h1 = h2;                                               
                                    }
                                 }
                                 rt = C0;
                                 
                                 //for(int32_t i=0; i<10; ++i) 
#if (USE_SLEEF_LIB) == 1
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[0],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[0],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[1],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[1],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[2],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[2],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[3],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[3],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[4],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[4],
                                               _mm512_div_pd(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[5],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[5],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[6],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[6],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[7],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[7],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[8],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[8],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[9],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[9],
                                               _mm512_div_ps(xexpf(_mm512_mul_ps(hs,t0)),t0)));
                                 
#else
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[0],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[0],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[1],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[1],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[2],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[2],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[3],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[3],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[4],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[4],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[5],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[5],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[6],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[6],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[7],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[7],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[8],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[8],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0)));
                                 x  = _mm512_fmadd_ps(C05,h2,_mm512_add_ps(xtab[9],C1));
                                 t0 = _mm512_add_ps(C1,_mm512_mul_ps(x,x));
                                 rt = _mm512_add_ps(rt,_mm512_mul_ps(weight[9],
                                               _mm512_div_ps(_mm512_exp_ps(_mm512_mul_ps(hs,t0)),t0))); 

#endif
                                 t1 = _mm512_mul_ps(C05,h2);
                                 tfn= _mm512_mul_ps(rt,_mm512_mul_ps(t1,twopinv)); 
		             }
		             return (tfn);
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
                      __m512d
                      von_misses_sample_zmm8r8(const __m512d a,
		                               const __m512d b) {

                          const __m512d  pi   = _mm512_set1_pd(3.14159265358979323846264338328);
			  const __m512d  _1   = _mm512_set1_pd(1.0);
			  const __m512d  _2   = _mm512_set1_pd(2.0);
			  const __m512d  _4   = _mm512_set1_pd(4.0);
			  const __m512d  _1_2 = _mm512_set1_pd(0.5);
			  __m512d c,f,rho,tau,u1,r;
			  __m512d u2,u3,x,z;
			  __m512d t0,t1,t2;
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
                             const __m512d nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			     return (nan);
			  }
			  uniform             = svrng_new_normal_distribution_double(0.0,1.0);
			  t0                  = _mm512_fmadd_pd(_4,_mm512_mul_pd(b,b),_1);
			  tau                 = _mm512_add_pd(_1,_mm512_sqrt_pd(t0));
			  t1                  = _mm512_add_pd(b,b);
			  rho                 = _mm512_div_pd(_mm512_sub_pd(tau,
			                                                _mm512_sqrt_pd(_mm512_add_pd(tau,tau))),t1);
			  t2                  = _mm512_fmadd_pd(rho,rho,_1);
			  r                   = _mm512_div_pd(t2,_mm512_add_pd(rho,rho));
            
 			 while(true) {
                               
                              const double * __restrict ptr = (const double*)(&svrng_generate8_double(engine,uniform));
                              u1                            = _mm512_loadu_pd(&ptr[0]);
#if (USE_SLEEF_LIB) == 1
			      z                             = xcos(_mm512_mul_pd(pi,u1));
#else
                              z                             = _mm512_cos_pd(_mm512_mul_pd(pi,u1));
#endif
                              f                             = _mm512_div_pd(_mm512_fmadd_pd(r,z,_1),
			                                                    _mm512_add_pd(r,z));
			      c                             = _mm512_mul_pd(b,_mm512_sub_pd(r,f));
			      t0                            = _mm512_mul_pd(c,_mm512_sub_pd(_2,c));
			                       
			      if(_mm512_cmp_mask_pd(u2,t0,_CMP_LT_OQ)) break;
			      t1                            = _mm512_add_pd(_mm512_log_pd(
			                                                  _mm512_div_pd(c,u2)),_1);
			      if(_mm512_cmp_mask_pd(c,t1,_CMP_LE_OQ)) break;
			 }
			 const double * __restrict ptr2 =
			                    (const double*)(&svrng_generate8_double(engine,uniform));
			 u3                             = _mm512_loadu_pd(&ptr2[0]);
		         t2                             = zmm8r8_sign_zmm8r8(_1,_mm512_sub_pd(u3,_1_2));
			 x                              = _mm512_fmadd_pd(t2,_mm512_acos_pd(f),a);
			 svrng_delete_engine(engine);
			 return (x)
		   }


		   
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
                      von_misses_sample_zmm16r4(const __m512 a,
		                                const __m512 b) {

                          const __m512   pi   = _mm512_set1_ps(3.14159265358979323846264338328f);
			  const __m512   _1   = _mm512_set1_ps(1.0f);
			  const __m512  _2    = _mm512_set1_ps(2.0f);
			  const __m512  _4    = _mm512_set1_ps(4.0f);
			  const __m512  _1_2  = _mm512_set1_ps(0.5f);
			  __m512 c,f,rho,tau,u1,r;
			  __m512 u2,u3,x,z;
			  __m512 t0,t1,t2;
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
                             const __m512 nan = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
			     return (nan);
			  }
			  uniform             = svrng_new_normal_distribution_float(0.0f,1.0f);
			  t0                  = _mm512_fmadd_ps(_4,_mm512_mul_ps(b,b),_1);
			  tau                 = _mm512_add_ps(_1,_mm512_sqrt_ps(t0));
			  t1                  = _mm512_add_ps(b,b);
			  rho                 = _mm512_div_ps(_mm512_sub_ps(tau,
			                                                _mm512_sqrt_ps(_mm512_add_ps(tau,tau))),t1);
			  t2                  = _mm512_fmadd_ps(rho,rho,_1);
			  r                   = _mm512_div_ps(t2,_mm512_add_ps(rho,rho));
            
 			 while(true) {
                               
                              const float * __restrict ptr = (const float*)(&svrng_generate16_float(engine,uniform));
                              u1                            = _mm512_loadu_ps(&ptr[0]);
#if (USE_SLEEF_LIB) == 1
			      z                             = xcosf(_mm512_mul_ps(pi,u1));
#else
                              z                             = _mm512_cos_ps(_mm512_mul_ps(pi,u1));
#endif
                              f                             = _mm512_div_ps(_mm512_fmadd_ps(r,z,_1),
			                                                    _mm512_add_ps(r,z));
			      c                             = _mm512_mul_ps(b,_mm512_sub_ps(r,f));
			      t0                            = _mm512_mul_ps(c,_mm512_sub_ps(_2,c));
			                       
			      if(_mm512_cmp_mask_ps(u2,t0,_CMP_LT_OQ)) break;
			      t1                            = _mm512_add_ps(_mm512_log_ps(
			                                                  _mm512_div_ps(c,u2)),_1);
			      if(_mm512_cmp_mask_ps(c,t1,_CMP_LE_OQ)) break;
			 }
			 const float * __restrict ptr2 =
			                    (const float*)(&svrng_generate16_float(engine,uniform));
			 u3                             = _mm512_loadu_ps(&ptr2[0]);
		         t2                             = zmm16r4_sign_zmm16r4(_1,_mm512_sub_ps(u3,_1_2));
			 x                              = _mm512_fmadd_ps(t2,_mm512_acos_ps(f),a);
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
                      __m512d
		      rayleigh_pdf_zmm8r8(const __m512d x,
		                          const __m512d a) {

                           const __m512d  _0 = _mm512_setzero_pd();
			   __m512d t0,t1,t2,t3,pdf;
			   const __mmask8 m  = _mm512_cmp_pd_mask(x,_0,_CMP_LT_OQ);
			   t0                = _mm512_mul_pd(a,a);
			   t1                = zmm8r8_negate(_mm512_div_pd(_mm512_mul_pd(x,x),
			                                                   _mm512_add_pd(t0,t0)));
			   t2                = _mm512_div_pd(x,t0);
#if (USE_SLEEF_LIB) == 1
                           t3               = _mm512_mul_pd(t2,xexp(t1));
#else
			   t3               = _mm512_mul_pd(t2,_mm512_exp_pd(t1));
#endif
                           pdf              = _mm512_mask_blend_pd(m,t3,_0);
                           return (pdf);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      rayleigh_pdf_zmm16r4(const __m512 x,
		                           const __m512 a) {

                           const __m512  _0 = _mm512_setzero_ps();
			   __m512 t0,t1,t2t3,pdf;
			   const __mmask16 m  = _mm512_cmp_ps_mask(x,_0,_CMP_LT_OQ);
			   t0                = _mm512_mul_ps(a,a);
			   t1                = zmm16r4_negate(_mm512_div_ps(_mm512_mul_ps(x,x),
			                                                   _mm512_add_ps(t0,t0)));
			   t2                = _mm512_div_ps(x,t0);
#if (USE_SLEEF_LIB) == 1
                           t3                = _mm512_mul_ps(t2,xexpf(t1));
#else
			   t3                = _mm512_mul_ps(t2,_mm512_exp_ps(t1));
#endif
                           pdf               = _mm512_mask_blend_ps(m,_t3,_0);
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
                      __m512d
		      rayleigh_mean_zmm8r8(const __m512d a) {

                          const __m512d hpi =  _mm512_set1_pd(0.5*3.14159265358979323846264338328);
			  __m512d mean;
			  mean              =  _mm512_mul_pd(a,_mm512_sqrt_pd(hpi));
			  return (mean);
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      rayleigh_mean_zmmr16r4(const __m512d a) {

                          const __m512 hpi =  _mm512_set1_ps(0.5f*3.14159265358979323846264338328f);
			  __m512 mean;
			  mean              =  _mm512_mul_ps(a,_mm512_sqrt_ps(hpi));
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
                      __m512d
		      rayleigh_invcdf_zmm8r8(const __m512d cdf,
		                             const __m512d a) {

			 const __m512d _0 = _mm512_setzero_pd();
			 const __m512d _1 = _mm512_setzero_pd(1.0);
			 const __m512d n2 = _mm512_setzero_pd(-2.0);
			 __m512d inv,t0,t1,;
                       //  if(__builtin_expect(_mm512_cmp_pd_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//    __builtin_expect(_mm512_cmp_pd_mask(_1,cdf,_CMP_LT_OQ),0)) {return;}
			 t0  = _mm512_log_pd(_mm512_sub_pd(_1,cdf));
			 t1  = _mm512_mul_pd(_2,_mm512_mul_pd(a,a));
                         inv = _mm512_sqrt_pd(_mm512_mul_pd(t0,t1));
			 return (inv);
			   
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      rayleigh_invcdf_zmm16r4(const __m512 cdf,
		                             const __m512 a) {

			 const __m512 _0 = _mm512_setzero_ps();
			 const __m512 _1 = _mm512_setzero_ps(1.0f);
			 const __m512 n2 = _mm512_setzero_ps(-2.0f);
			 __m512 inv,t0,t1,;
                       //  if(__builtin_expect(_mm512_cmp_ps_mask(cdf,_0,_CMP_LT_OQ),0) ||
			//    __builtin_expect(_mm512_cmp_ps_mask(_1,cdf,_CMP_LT_OQ),0)) {return;}
			 t0  = _mm512_log_ps(_mm512_sub_ps(_1,cdf));
			 t1  = _mm512_mul_ps(_2,_mm512_mul_ps(a,a));
                         inv = _mm512_sqrt_ps(_mm512_mul_ps(t0,t1));
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
                      __m512d
		      rayleigh_cdf_zmm8r8(const __m512d x,
		                          const __m512d a) {

                         const __m512d _0 = _mm512_setzero_pd();
			 const __m512d _1 = _mm512_setzero_pd(1.0);
			 __m512d cdf,t0,t1;
			 t0              = _mm512_mul_pd(_2,_mm512_mul_pd(a,a));
			 t1              = zmm8r8_negate(_mm512_mul_pd(x,x));
#if (USE_SLEEF_LIB) == 1
                         cdf             = _mm512_sub_pd(_1,xexp(_mm512_div_pd(t1,t0)));
			                            
#else
			 cdf             = _mm512_sub_pd(_1,
			                             _mm512_exp_pd(_mm512_div_pd(t1,t0)));
#endif
                         return (cdf);
		    }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      rayleigh_cdf_zmm16r4(const __m512 x,
		                           const __m512 a) {

                         const __m512 _0 = _mm512_setzero_ps();
			 const __m512 _1 = _mm512_setzero_ps(1.0f);
			 __m512 cdf,t0,t1;
			 t0              = _mm512_mul_pd(_2,_mm512_mul_ps(a,a));
			 t1              = zmm16r4_negate(_mm512_mul_ps(x,x));
#if (USE_SLEEF_LIB) == 1
                         cdf             = _mm512_sub_ps(_1,xexpf(_mm512_div_ps(t1,t0)));
			                            
#else
			 cdf             = _mm512_sub_ps(_1,
			                             _mm512_exp_ps(_mm512_div_ps(t1,t0)));
#endif
                         return (cdf);
		    }


		    
		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512d
		      rayleigh_sample_zmm8r8(const __m512d rand,
		                             const __m512d a) {

                          return (rayleigh_invcdf_zmm8r8(rand,a));
		     }


		      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
		      rayleigh_sample_zmm16r4(const __m512 rand,
		                             const __m512 a) {

                          return (rayleigh_invcdf_zmm16r4(rand,a));
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
                      __m512d    
                      cauchy_cdf_zmm8r8(const __m512d x,
                                        const __m512d a,
                                        const __m512d b) {
                        
                         const __m512d C314159265358979323846264 = 
                                               __m512_set1_pd(3.14159265358979323846264);
                         const __m512d C05 = _mm512_set1_pd(0.5);
                         register __m512d cdf,y,t0,t1;
                         t0 = _mm512_sub_pd(x,a);
#if (USE_SLEEF_LIB) == 1
                         t1 = _mm512_div_pd(xatan2(t0,b),
                                     C314159265358979323846264);
#else
                         t1 = _mm512_div_pd(_mm512_atan2_pd(t0,b),
                                     C314159265358979323846264);
#endif                  
                         cdf = _mm512_add_pd(C05,t1);
                         return (cdf);  
                    }
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512    
                      cauchy_cdf_zmm16r4(const __m512 x,
                                         const __m512 a,
                                         const __m512 b) {
                        
                         const __m512 C314159265358979323846264 = 
                                               __m512_set1_ps(3.14159265358979323846264f);
                         const __m512 C05 = _mm512_set1_ps(0.5f);
                         register __m512 cdf,y,t0,t1;
                         t0 = _mm512_sub_ps(x,a);
#if (USE_SLEEF_LIB) == 1
                         t1 = _mm512_div_ps(xatan2f(t0,b),
                                     C314159265358979323846264);
#else
                         t1 = _mm512_div_ps(_mm512_atan2_ps(t0,b),
                                     C314159265358979323846264);
#endif                  
                         cdf = _mm512_add_ps(C05,t1);
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
                      __m512d 
                      cauchy_cdf_inv_zmm8r8(const __m512d a,
                                            const __m512d b,
                                            const __m512d x) {
                           
                         const __m512d C314159265358979323846264 = 
                                               __m512_set1_pd(3.14159265358979323846264);
                         const __m512d C05 = _mm512_set1_pd(0.5);    
                         register __m512d cdf,t0,t1;
                         t0 = _mm512_mul_pd(C314159265358979323846264,
                                            _mm512_sub_pd(cdf,C05));
#if (USE_SLEEF_LIB) == 1
                         t1 = xtan(t0);
#else
                         t1 = _mm512_tan_pd(t0);
#endif                                      
                         cdf = _mm512_fmadd_pd(a,b,t1);
                         return (cdf);
                   }    
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512
                      cauchy_cdf_inv_zmm16r4(const __m512 a,
                                            const __m512 b,
                                            const __m512 x) {
                           
                         const __m512 C314159265358979323846264 = 
                                               __m512_set1_pd(3.14159265358979323846264f);
                         const __m512 C05 = _mm512_set1_ps(0.5);    
                         register __m512 cdf,t0,t1;
                         t0 = _mm512_mul_ps(C314159265358979323846264,
                                            _mm512_sub_ps(cdf,C05));
#if (USE_SLEEF_LIB) == 1
                         t1 = xtanf(t0);
#else
                         t1 = _mm512_tan_ps(t0);
#endif                                      
                         cdf = _mm512_fmadd_ps(a,b,t1);
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
                      __m512d  
                      cauchy_pdf_zmm8r8(const __m512d x,
                                         const __m512d a,
                                         const __m512d b) {
                           
                         const __m512d C314159265358979323846264 = 
                                               __m512_set1_pd(3.14159265358979323846264);
                         const __m512d C1 = _mm512_set1_pd(1.0);
                         register __m512d pdf,t0,t1,y,pib;
                         y   = _mm512_div_pd(_mm512_sub_pd(x,a),b);
                         pib = _mm512_mul_pd(C314159265358979323846264,b);
                         t0  = _mm512_fmadd_pd(y,y,C1);
                         t1  = _mm512_mul_pd(pib,t0);
                         pdf = _mm512_div_pd(C1,t1);
                         return (pdf);                     
                   }
                   
                   
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512 
                      cauchy_pdf_zmm16r4(const __m512 x,
                                         const __m512 a,
                                         const __m512 b) {
                           
                         const __m512 C314159265358979323846264 = 
                                               __m512_set1_ps(3.14159265358979323846264f);
                         const __m512 C1 = _mm512_set1_ps(1.0);
                         register __m512 pdf,t0,t1,y,pib;
                         y   = _mm512_div_ps(_mm512_sub_ps(x,a),b);
                         pib = _mm512_mul_ps(C314159265358979323846264,b);
                         t0  = _mm512_fmadd_ps(y,y,C1);
                         t1  = _mm512_mul_ps(pib,t0);
                         pdf = _mm512_div_ps(C1,t1);
                         return (pdf);                     
                   }
                   
/*
   !*****************************************************************************80
!
!! CAUCHY_SAMPLE samples the Cauchy PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
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
                      __m512d 
                      cauchy_sample_zmm8r8(const __m512d a,
                                           const __m512d b) {
                         __m512d cdf;
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
                            const __m512d nan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_normal_distribution_double(0.0,1.0);
			 const double * __restrict ptr = (const double*)(&svrng_generate8_double(engine,uniform));
			 cdf              = cauchy_cdf_inv_zmm8r8(_mm512_loadu_pd(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);                 
                    }
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512 
                      cauchy_sample_zmm16r4(const __m512 a,
                                            const __m512 b) {
                         __m512 cdf;
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
                            const __m512 nan = _mm512_set1_ps(std::numeric_limits<float>::quiet_NaN());
			    return (nan);
			 }
			 uniform          = svrng_new_normal_distribution_float(0.0f,1.0f);
			 const float * __restrict ptr = (const float*)(&svrng_generate16_float(engine,uniform));
			 cdf              = cauchy_cdf_inv_zmm16r4(_mm512_loadu_ps(&ptr[0]));
			 svrng_delete_engine(engine);
			 return (cdf);                 
                    }
#if 0                    
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
#endif


                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512d 
                      maxwell_cdf_zmm8r8(const __m512d x,
                                         const __m512d a) {
                         
                         const __m512d C15 = _mm512_set1_pd(1.5);
                         register __m512d x2,cdf;
                         x2 = _mm512_div_pd(x,a);
                         cdf = gamma_incomplete_zmm8r8(C15,x2);
                         return (cdf);                      
                    }      
                    
                    
                      __ATTR_REGCALL__
                      __ATTR_ALWAYS_INLINE__
		      __ATTR_HOT__
		      __ATTR_ALIGN__(32)
		      static inline           
                      __m512 
                      maxwell_cdf_zmm16r4(const __m512 x,
                                         const __m512 a) {
                         
                         const __m512 C15 = _mm512_set1_ps(1.5f);
                         register __m512 x2,cdf;
                         x2 = _mm512_div_ps(x,a);
                         cdf = gamma_incomplete_zmm16r4(C15,x2);
                         return (cdf);                      
                    }   
                    



      } //math

} // gms














#endif /*__GMS_PDF_CDF_AVX512_HPP__*/
