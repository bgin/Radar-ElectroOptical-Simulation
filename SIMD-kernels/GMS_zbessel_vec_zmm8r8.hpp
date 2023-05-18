

#ifndef __GMS_ZBESSEL_VEC_ZMM8R8_HPP__
#define __GMS_ZBESSEL_VEC_ZMM8R8_HPP__


/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_version {

    const unsigned int GMS_ZBESSEL_VEC_ZMM8R8_MAJOR = 1U;
    const unsigned int GMS_ZBESSEL_VEC_ZMM8R8_MINOR = 0U;
    const unsigned int GMS_ZBESSEL_VEC_ZMM8R8_MICRO = 0U;
    const unsigned int GMS_ZBESSEL_VEC_ZMM8R8_FULLVER =
      1000U*GMS_ZBESSEL_VEC_ZMM8R8_MAJOR+
      100U*GMS_ZBESSEL_VEC_ZMM8R8_MINOR+
      10U*GMS_ZBESSEL_VEC_ZMM8R8_MICRO;
    const char * const GMS_ZBESSEL_VEC_ZMM8R8_CREATION_DATE = "18-05-2023 08:38 AM +00200 (THR 18 05 2023 GMT+2)";
    const char * const GMS_RCS_COMMON_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_COMMON_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_COMMON_ZMM16R4_DESCRIPTION   = "Amos BESSEL function package manually vectorized (avx512 double)."

}

#include <cstdint>
#include <immintrin.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "GMS_config.h"
#include "GMS_complex_zmm8r8.hpp"
#include "GMS_simd_utils.hpp"
#include "GMS_sleefsimddp.hpp"




namespace gms {


        namespace math {
        
        
               /*
               
                     /* ***BEGIN PROLOGUE  ZAIRY */
  /* ***PURPOSE  Compute the Airy function Ai(z) or its derivative dAi/dz */
  /*            for complex argument z.  A scaling option is available */
  /*            to help avoid underflow and overflow. */
  /* ***LIBRARY   SLATEC */
  /* ***CATEGORY  C10D */
  /* ***TYPE      COMPLEX (CAIRY-C, ZAIRY-C) */
  /* ***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD, */
  /*             BESSEL FUNCTION OF ORDER TWO THIRDS */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*                      ***A DOUBLE PRECISION ROUTINE*** */
  /*         On KODE=1, ZAIRY computes the complex Airy function Ai(z) */
  /*         or its derivative dAi/dz on ID=0 or ID=1 respectively. On */
  /*         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz */
  /*         is provided to remove the exponential decay in -pi/3<arg(z) */
  /*         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where */
  /*         zeta=(2/3)*z**(3/2). */

  /*         While the Airy functions Ai(z) and dAi/dz are analytic in */
  /*         the whole z-plane, the corresponding scaled functions defined */
  /*         for KODE=2 have a cut along the negative real axis. */

  /*         Input */
  /*           ZR     - DOUBLE PRECISION real part of argument Z */
  /*           ZI     - DOUBLE PRECISION imag part of argument Z */
  /*           ID     - Order of derivative, ID=0 or ID=1 */
  /*           KODE   - A parameter to indicate the scaling option */
  /*                    KODE=1  returns */
  /*                            AI=Ai(z)  on ID=0 */
  /*                            AI=dAi/dz on ID=1 */
  /*                            at z=Z */
  /*                        =2  returns */
  /*                            AI=exp(zeta)*Ai(z)  on ID=0 */
  /*                            AI=exp(zeta)*dAi/dz on ID=1 */
  /*                            at z=Z where zeta=(2/3)*z**(3/2) */
   /*         Output */
  /*           AIR    - DOUBLE PRECISION real part of result */
  /*           AII    - DOUBLE PRECISION imag part of result */
  /*           NZ     - Underflow indicator */
  /*                    NZ=0    Normal return */
  /*                    NZ=1    AI=0 due to underflow in */
  /*                            -pi/3<arg(Z)<pi/3 on KODE=1 */
  /*           IERR   - Error flag */
  /*                    IERR=0  Normal return     - COMPUTATION COMPLETED */
  /*                    IERR=1  Input error       - NO COMPUTATION */
  /*                    IERR=2  Overflow          - NO COMPUTATION */
  /*                            (Re(Z) too large with KODE=1) */
  /*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
  /*                            (Result has less than half precision) */
  /*                    IERR=4  Precision error   - NO COMPUTATION */
  /*                            (Result has no precision) */
  /*                    IERR=5  Algorithmic error - NO COMPUTATION */
  /*                            (Termination condition not met) */

  /* *Long Description: */

  /*         Ai(z) and dAi/dz are computed from K Bessel functions by */

  /*                Ai(z) =  c*sqrt(z)*K(1/3,zeta) */
  /*               dAi/dz = -c*   z   *K(2/3,zeta) */
  /*                    c =  1/(pi*sqrt(3)) */
  /*                 zeta =  (2/3)*z**(3/2) */

  /*         when abs(z)>1 and from power series when abs(z)<=1. */

  /*         In most complex variable computation, one must evaluate ele- */
  /*         mentary functions.  When the magnitude of Z is large, losses */
  /*         of significance by argument reduction occur.  Consequently, if */
  /*         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR), */
  /*         then losses exceeding half precision are likely and an error */
  /*         flag IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is */
  /*         double precision unit roundoff limited to 18 digits precision. */
  /*         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then */
  /*         all significance is lost and IERR=4.  In order to use the INT */
  /*         function, ZETA must be further restricted not to exceed */
  /*         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA */
  /*         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, */
  /*         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single */
  /*         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision. */
  /*         This makes U2 limiting is single precision and U3 limiting */
  /*         in double precision.  This means that the magnitude of Z */
  /*         cannot exceed approximately 3.4E+4 in single precision and */
  /*         2.1E+6 in double precision.  This also means that one can */
  /*         expect to retain, in the worst cases on 32-bit machines, */
  /*         no digits in single precision and only 6 digits in double */
  /*         precision. */

  /*         The approximate relative error in the magnitude of a complex */
  /*         Bessel function can be expressed as P*10**S where P=MAX(UNIT */
  /*         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre- */
  /*         sents the increase in error due to argument reduction in the */
  /*         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))), */
  /*         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF */
  /*         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may */
  /*         have only absolute accuracy.  This is most likely to occur */
  /*         when one component (in magnitude) is larger than the other by */
  /*         several orders of magnitude.  If one component is 10**K larger */
  /*         than the other, then one can expect only MAX(ABS(LOG10(P))-K, */
  /*         0) significant digits; or, stated another way, when K exceeds */
  /*         the exponent of P, no significant digits remain in the smaller */
  /*         component.  However, the phase angle retains absolute accuracy */
  /*         because, in complex arithmetic with precision P, the smaller */
  /*         component will not (as a rule) decrease below P times the */
  /*         magnitude of the larger component. In these extreme cases, */
  /*         the principal phase angle is on the order of +P, -P, PI/2-P, */
  /*         or -PI/2+P. */
   /* ***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe- */
  /*                 matical Functions, National Bureau of Standards */
  /*                 Applied Mathematics Series 55, U. S. Department */
  /*                 of Commerce, Tenth Printing (1972) or later. */
  /*               2. D. E. Amos, Computation of Bessel Functions of */
  /*                 Complex Argument and Large Order, Report SAND83-0643, */
  /*                 Sandia National Laboratories, Albuquerque, NM, May */
  /*                 1983. */
  /*               3. D. E. Amos, A Subroutine Package for Bessel Functions */
  /*                 of a Complex Argument and Nonnegative Order, Report */
  /*                 SAND85-1018, Sandia National Laboratory, Albuquerque, */
  /*                 NM, May 1985. */
  /*               4. D. E. Amos, A portable package for Bessel functions */
  /*                 of a complex argument and nonnegative order, ACM */
  /*                 Transactions on Mathematical Software, 12 (September */
  /*                 1986), pp. 265-273. */

  /* ***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACAI, ZBKNU, ZEXP, ZSQRT */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   890801  REVISION DATE from Version 3.2 */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   920128  Category corrected.  (WRB) */
  /*   920811  Prologue revised.  (DWL) */
  /*   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZAIRY */
  /*     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3 */
  /* ***FIRST EXECUTABLE STATEMENT  ZAIRY */
               */    
               
               
                   __ATTR_ALWAYS_INLINE__
	           __ATTR_HOT__
	           __ATTR_ALIGN__(32)
                   __ATTR_VECTORCALL__
	           static inline
	           int32_t zairy_zmm8r8(const __m512 zr,
	                             const __m512 zi,
	                             const int32_t id,
	                             const int32_t kode,
	                             __m512  * __restrict air,
	                             __m512  * __restrict aii,
	                             int32_t * __restrict nz) {
	                             
	               const __m512d R1M5 = _mm512_set1_ps(std::log10(
	                                                        std::numerical_limits<double>::radix));
	               const __m512d C0666666666666666666666666666667 = 
	                                   _mm512_set1_ps(0.666666666666666666666666666667); //2/3
	               const __m512d C035502805388781724              =            
	                                   _mm512_set1_ps(0.35502805388781724); //  const double c1 = .35502805388781724;
	               const __m512d C02588194037   	             =
	                                   _mm512_set1_ps(0.2588194037); //  const double c2 = .258819403792806799;
	               const __m512d C0183776298473930683             = 
	                                   _mm512_set1_ps(0.183776298473930683); // const double coef = .183776298473930683;
	               const __m512d C00                              =
	                                   _mm512_setzero_ps();
	               const __m512d C10                              =                   
	                                   _mm512_set1_ps(1.0);
	               const __m512d C20                              =
	                                   _mm512_set1_ps(2.0);
	               const __m512d C30                              =
	                                   _mm512_set1_ps(3.0);
	               const __m512d C40                              =
	                                   _mm512_set1_ps(4.0);
	               const __m512d C90                              =
	                                   _mm512_set1_ps(9.0);
	               const __m512d C300                             = 
	                                   _mm512_set1_ps(30.0);
	               const __m512d C240                             =
	                                   _mm512_set1_ps(24.0);
	               const __m512d C2303                            =
	                                   _mm512_set1_ps(2.303);
	               const __m512d C4145                            =
	                                   _mm512_set1_ps(-41.45);
	               const __m512d C12                              =
	                                   _mm512_set1_ps(1.0);
	               const __m512d C05                              = 
	                                   _mm512_set1_ps(0.5);
	               const __m512d C180                             =
	                                   _mm512_set1_ps(18.0);
	               const __m512d C033333333333333333333333        =
	                                   _mm512_set1_ps(0.333333333333333333333333333333);
	               const __m512d C025                             =
	                                   _mm512_set1_ps(0.25);
	               
	               __m512d cyi[1];
	               __m512d cyr[1];  
	               __m512d d1,d2;
	               __m512d aa,bb,ad,cc,ak,bk,ck,dk,az;
	               __m512d r1;
	               __m512d s1i,az3,s2i,s1r,s2r,z3i,dig,fid;
	               __m512d fnu,tol,sti,ptr,str,sfac,alim;
	               __m512d elim,alaz,csqi;
	               __m512d atrm,ztai,csqr,ztar;
	               __m512d trm1i,trm2i,trm1r,trm2r;
	               int32_t ierr,k,k1,k2,nn,mr,iflag;    
	               
	               ierr = 0;
	               *nz  = 0;
	               if(__builtin_expect(id<0,0) || 
	                  __builtin_expect(id>1,0)) {
	                  ierr = 1;
	               }  
	               if(__builtin_expect(kode<1,0) || 
	                  __builtin_expect(kode>2,0)) {
	                  ierr = 1;
	               }  
	               if(ierr!=0) {
	                   return (ierr);
	               }
	               az   = cabs_zmm8r8(zr,zi);
	               tol  = _mm512_set1_pd(std::numeric_limits<double>::epsilon(),1e-18);
	               fid  = _mm512_set1_pd((double)id);
	               if(_mm512_cmp_pd_mask(az,C10,_CMP_GT_OQ)) {
	                  goto L70;
	               }  
	               /* ----------------------------------------------------------------------- */
                       /*     POWER SERIES FOR ABS(Z).LE.1. */
                       /* ----------------------------------------------------------------------- */  
                       s1r = C10;
                       s1i = C00;
                       s2r = C10;
                       s2i = C00;
                       if(_mm512_cmp_pd_mask(az,tol,_CMP_LT_OQ)) {
                          goto L170;
                       } 
                       aa = _mm512_mul_pd(az,az);
                       if(_mm512_cmp_pd_mask(aa,
                                _mm512_div_pd(tol,az),_CMP_LT_OQ)) {
                          goto L40;        
                       }           
                       trm1r = C10;
                       trm1i = C00;
                       trm2r = C20;
                       trm2i = C00;
                       atrm  = C10;
                       str   = _mm512_fmsub_pd(zr,zr,_mm512_mul_ps(zi,zi));
                       sti   = _mm512_fmadd_pd(zr,zi,_mm512_mul_ps(zi,zr));
                       z3r   = _mm512_fmsub_pd(str,zr,_mm512_mul_ps(sti,zi));
                       z3i   = _mm512_fmadd_pd(str,zi,_mm512_mul_ps(sti,zr));
                       az3   = _mm512_mul_pd(az,aa);
                       ak    = _mm512_add_pd(fid,C20);
                       bk    = _mm512_sub_pd(C30,_mm512_sub_ps(fid,fid));
                       ck    = _mm512_sub_pd(C40,fid);
                       dk    = _mm512_add_pd(_mm512_add_ps(fid,C30),fid);
                       d1    = _mm512_mul_pd(ak,dk);
                       d2    = _mm512_mul_pd(bk,ck);
                       ad    = _mm512_min_pd(d1,d2);
                       ak    = _mm512_fmadd_pd(fid,C90,C240);
                       bk    = _mm512_sub_pd(C300,_mm512_mul_pd(fid,C90));
                       for(k = 1; k <= 25; ++k) {
                           str    = _mm512_div_pd(_mm512_fmsub_pd(trm1r,z3r,
                                                     _mm512_mul_pd(trm1i,z3i)),d1);
                           trm1i  = _mm512_div_pd(_mm512_fmadd_pd(trm1r,z3i,
                                                     _mm512_mul_pd(trm1i,z3r)),d1);
                           s1r    = _mm512_add_pd(s1r,trm1r);
                           s1i    = _mm512_add_pd(s1i,trm1i);
                           str    = _mm512_div_pd(_mm512_fmsub_pd(trm2r,z3r,
                                                     _mm512_mul_pd(trm2i,z3i)),d2);
                           trm2i  = _mm512_div_pd(_mm512_fmadd_pd(trm2r,z3i,
                                                     _mm512_mul_pd(trm2i,z3r)),d2);
                           trm2r  = str;
                           s2r    = _mm512_add_pd(s2r,trm2r);
                           s2i    = _mm512_add_pd(s2i,trm2i);
                           atrm   = _mm512_div_pd(_mm512_mul_pd(atrm,az),ad);
                           d1     = _mm512_add_pd(d1,ak);
                           d2     = _mm512_add_pd(d2,bk);
                           ad     = _mm512_min_pd(d1,d2);
                           if(_mm512_cmp_mask_pd(atrm,
                                             _mm512_mul_pd(tol,ad),_CMP_LT_OQ)) {
                               goto L40;                 
                           }
                           ak = _mm512_add_pd(ak,C180);
                           bk = _mm512_add_pd(bk,C180);
                       }
                     L40:
                           if(id==1) {
                               goto L50;
                           }
                           *air = _mm512_fmsub_pd(s1r,c1,
                                                    _mm512_mul_pd(c2,
                                                               _mm512_fmsub_pd(zr,s2r,
                                                                           _mm512_mul_pd(zi,s2i))));
                           aii  = _mm512_fmsub_pd(s1i,c1,
                                                    _mm512_mul_pd(c2,
                                                               _mm512_fmadd_pd(zr,s2i,
                                                                           _mm512_mul_pd(zi,s2r))));
                           if(kode==1) {
                              return ierr;
                           }
                           csqrt_zmm8r8(zr,zi,&str,&sti);
                           ztar = _mm512_mul_pd(C0666666666666666666666666666667,
                                                          _mm512_fmsub_pd(zr,str,
                                                                      _mm512_mul_pd(zi,sti)));
                           ztai = _mm512_mul_pd(C0666666666666666666666666666667,
                                                          _mm512_fmadd_pd(zr,sti,
                                                                      _mm512_mul_pd(zi,str)));
                           cexp_zmm8r8(ztar,ztai,&str,&sti);
                           ptr  = _mm512_fmsub_pd(*air,str,_mm512_mul_pd(*aii,sti));
                           *aii = _mm512_fmadd_pd(*air,sti,_mm512_mul_pd(*aii,str));
                           *air = ptrl
                           return ierr;
                       L50:
                           *air = _mm512_mul_pd(negate_zmm8r8(s2r),c2);
                           *aii = _mm512_mul_pd(negate_zmm8r8(s2i),c2);
                           if(_mm512_cmp_pd_mask(az,tol,_CMP_LE_OQ)) {
                              goto L60;
                           }
                           str = _mm512_fmsub_pd(zr,s1r,_mm512_mul_pd(zi,s1i));
                           sti = _mm512_fmadd_pd(zr,s1i,_mm512_mul_pd(zi,s1r));
                           cc  = _mm512_div_pd(c1,_mm512_add_pd(fid,C10));
                           *air= _mm512_fmadd_pd(_mm512_fmsub_pd(str,zr,
                                                             _mm512_mul_pd(sti,zi)),cc,*air);
                           *aii= _mm512_fmadd_pd(_mm512_fmadd_pd(str,zi,
                                                             _mm512_mul_pd(sti,zr)),cc,*aii);
                       L60:
                           if(kode==1) {
                              return ierr;
                           }
                           csqrt_zmm8r8(zr,zi,&str,&sti);
                           ztar = _mm512_mul_pd(C0666666666666666666666666666667,
                                                          _mm512_fmsub_pd(zr,str,
                                                                      _mm512_mul_pd(zi,sti)));
                           ztai = _mm512_mul_pd(C0666666666666666666666666666667,
                                                          _mm512_fmadd_pd(zr,sti,
                                                                      _mm512_mul_pd(zi,str)));
                           cexp_zmm8r8(ztar,ztai,&str,&sti);
                           ptr  = _mm512_fmsub_pd(str,*air,_mm512_mul_pd(sti,*aii));
                           *aii = _mm512_fmadd_pd(str,*aii,_mm512_mul_pd(sti,*air));
                           *air = ptr;
                           return ierr;
                      /* ----------------------------------------------------------------------- */
                      /*     CASE FOR ABS(Z).GT.1.0 */
                      /* ----------------------------------------------------------------------- */    
                      L70:
                           fnu = _mm512_mul_pd(_mm512_add_pd(fid,C10), 
                                                         C033333333333333333333333);
                       /* ----------------------------------------------------------------------- */
                       /*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
                       /*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18. */
                       /*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
                       /*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
                       /*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
                       /*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
                       /*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
                       /*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
                       /* ----------------------------------------------------------------------- */
                           k1   = std::numeric_limits<double>::min_exponent;
                           k2   = std::numeric_limits<double>::max_exponent;
                            /* Computing MIN */
                           k    = std::min(std::abs(k1), std::abs(k2));
                           elim = _mm512_mul_pd(_mm512_fmsub_pd(_mm512_set1_pd((double)k),
                                                                              R1M5,C30),C2303);
                           aa   = _mm512_mul_pd(R1M5,_mm512_set1_pd((double)k1));
                           dig  = _mm512_min_pd(aa,C180);
                           aa   = _mm512_mul_pd(aa,C2303);
                             /* Computing MAX */
                           alim = _mm512_add_pd(elim,_mm512_max_pd(negate_zmm8r8(aa),C4145));
                           r1   = _mm512_fmadd_pd(dig,C120,C30);
                           alaz = _mm512_log_pd(az);
                             /* ----------------------------------------------------------------------- */
                             /*     TEST FOR PROPER RANGE */
                             /* ----------------------------------------------------------------------- */
                           aa = _mm512_div_pd(C05,tol);
                           bb = _mm512_mul_pd(_mm512_set1_pd(
                                                  (std::numeric_limits<int>::max()*0.5)));
                           aa = _mm512_min_pd(aa,bb);
                           aa = _mm512_pow_pd(aa,C0666666666666666666666666666667);
                           if(_mm512_cmp_pd_mask(az,aa,_CMP_GT_OQ)) {
                              goto L260;
                           }
                           aa = _mm512_sqrt_pd(aa);
                           if(_mm512_cmp_pd_mask(az,aa,_CMP_GT_OQ)) {
                              ierr = 3;
                           }
                           csqrt_zmm8r8(zr,zi,&csqr,&csqi);
                           ztar = _mm512_mul_pd(C0666666666666666666666666666667,
                                                              _mm512_fmsub_pd(zr,csqr,
                                                                          _mm512_mul_pd(zi,csqi)));
                           ztai = _mm512_mul_pd(C0666666666666666666666666666667,
                                                              _mm512_fmadd_pd(zr,csqi,
                                                                          _mm512_mul_pd(zi,csqr)));
                           /* ----------------------------------------------------------------------- */
                           /*     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL */
                           /* ----------------------------------------------------------------------- */
                           iflag = 0; 
                           sfac  = C10;
                           ak    = ztai;
                           if(_mm512_cmp_pd_mask(zr,C00,_CMP_GE_OQ)) {
                              goto L80;
                           }
                           bk    = ztar;
                           ck    = negate_zmm8r8(_mm512_abs_pd(bk));
                           ztar  = ck;
                           ztai  = ci;
                       L80:
                           if(_mm512_cmp_pd_mask(zi,C00,_CMP_NE_OQ)) {
                              goto L90;
                           }
                           if(_mm512_cmp_pd_mask(zr,C00,_CMP_GT_OQ)) {
                              goto L90;
                           }
                           ztar  = C00;
                           ztai  = ak;
                       L90:
                           aa    = ztar;
                           if(_mm512_cmp_pd_mask(aa,C00,_CMP_GE_OQ) &&
                              _mm512_cmp_pd_mask(zr,C00,_CMP_GT_OQ)) {
                              goto L110;  
                           }
                           if(kode==2) {
                              goto L100;
                           }
                           /* ----------------------------------------------------------------------- */
                           /*     OVERFLOW TEST */
                           /* ----------------------------------------------------------------------- */ 
	       }
        
        
     } // math




} // gms































#endif /*__GMS_ZBESSEL_VEC_ZMM8R8_HPP__*/
