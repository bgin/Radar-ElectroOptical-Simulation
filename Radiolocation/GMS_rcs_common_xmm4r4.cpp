


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


#include "GMS_rcs_common_xmm4r4.h"






                       


                 /*

c*********************************************************************72
c
cc RD computes an incomplete elliptic integral of the second kind, RD(X,Y,Z).
c
c  Discussion:
c
c    This function computes an incomplete elliptic integral of the second kind.
c
c    RD(X,Y,Z) = Integral ( 0 <= T < oo )
c
c                                -1/2     -1/2     -3/2
c                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,
c
c    where X and Y are nonnegative, X + Y is positive, and Z is positive.
c
c    If X or Y is zero, the integral is complete.
c
c    The duplication theorem is iterated until the variables are
c    nearly equal, and the function is then expanded in Taylor
c    series to fifth order.  
c
c    Check: 
c
c      RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) = 3 / sqrt ( X * Y * Z ), 
c      where X, Y, and Z are positive.
c
c  Modified:
c
c    27 May 2018
c
c  Author:
c
c    Bille Carlson, Elaine Notis
c
c  Reference:
c
c    Bille Carlson,
c    Computing Elliptic Integrals by Duplication,
c    Numerische Mathematik,
c    Volume 33, 1979, pages 1-16.
c
c    Bille Carlson, Elaine Notis,
c    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 3, pages 398-403, September 1981.
c
c  Parameters:
c
c    Input, double precision X, Y, Z, the arguments in the integral.
c
c    Input, double precision ERRTOL, the error tolerance.
c    The relative error due to truncation is less than
c      3 * ERRTOL ^ 6 / (1-ERRTOL) ^ 3/2.
c    Sample choices:
c      ERRTOL   Relative truncation error less than
c      1.D-3    4.D-18
c      3.D-3    3.D-15
c      1.D-2    4.D-12
c      3.D-2    3.D-9
c      1.D-1    4.D-6
c
c    Output, integer IERR, the error flag.
c    0, no error occurred.
c    1, abnormal termination.
c
                   */  


                   
	          
	          
                  
	         
                   __m128 gms::radiolocation::rd_xmm4r4(const __m128 x,
                                     const __m128 y,
                                     const __m128 z,
                                     const __m128 errtot) {

                          const register __m128 _3 = _mm_set1_ps(3.0f);
                          const register __m128 _1 = _mm_set1_ps(1.0f);
                          const register __m128 c1 = _mm_set1_ps(-0.214285714285714285714285714286f);
                          const register __m128 c2 = _mm_set1_ps(0.166666666666666666666666666667f);
                          const register __m128 c3 = _mm_set1_ps(-0.409090909090909090909090909091f);
                          const register __m128 c4 = _mm_set1_ps(0.115384615384615384615384615385f);
                          const register __m128 c5 = _mm_set1_ps(6.0f);
                          const register __m128 c6 = _mm_set1_ps(1.5f);
                          const register __m128 c7 = _mm_set1_ps(0.2f);
                          const register __m128 c8 = _mm_set1_ps(0.25f);
                          register __m128 rd,xn,yn,zn,epslon,sigma,pow4,mu;
                          register __m128 xndev,yndev,zndev,ea,eb,ec,ed,ef;
                          register __m128 s1,s2,xnroot,ynroot,znroot,lamda;
                          register __m128 x0,x1,x2,x3,x4,x5;

                          xn    = x;
                          yn    = y;
                          zn    = z;
                          sigma = _mm_setzero_ps();
                          pow4  = _1; 
                          while(true) {
                                mu    = _mm_mul_ps(c7,_mm_fmadd_ps(zn,_3,
                                                          _mm_add_ps(xn,yn));
                                xndev = _mm_div_ps(_mm_sub_ps(mu,xn),mu);
                                yndev = _mm_div_ps(_mm_sub_ps(mu,yn),mu);
                                zndev = _mm_div_ps(_mm_sub_ps(mu,zn),mu);
                                epslon= _mm_abs_ps(xndev);
                                epslon= _mm_max_ps(epslon,_mm_abs_ps(yndev));
                                epslon= _mm_max_ps(epslon,_mm_abs_ps(zndev));

                                if(_mm_cmp_mask_ps(epslon,errtot,_CMP_LT_OQ)) {
                                   ea = _mm_mul_ps(xndev,yndev);
                                   eb = _mm_mul_ps(zndev,zndev);
                                   ec = _mm_sub_ps(ea,eb);
                                   ed = _mm_sub_ps(ea,_mm_mul_ps(c5,eb));
                                   ef = _mm_add_ps(ed,_mm_add_ps(ec,ec));
                                   x0 = _mm_fmadd_ps(c3,c8,c1);
                                   x1 = _mm_sub_ps(ed,_mm_sub_ps(c6,c4));
                                   x2 = _mm_mul_ps(zndev,ef);
                                   s1 = _mm_mul_ps(ed,_mm_mul_ps(x0,
                                                               _mm_mul_ps(x1,x2)));
                                   x3 = _mm_fmadd_ps(c3,ec,_mm_mul_ps(zndev,
                                                               _mm_mul_ps(c4,ea)));
                                   x4 = _mm_fmadd_ps(x3,zndev,_mm_mul_ps(ef,c2));
                                   s2 = _mm_mul_ps(zndev,x4);
                                   x0 = _mm_fmadd_ps(_3,sigma,pow4);
                                   x1 = _mm_add_ps(_1,_mm_add_ps(s1,s2));
                                   x2 = _mm_mul_ps(mu,_mm_sqrt_ps(mu));
                                   rd = _mm_div_ps(_mm_mul_ps(x0,x1),x2);
                                   return (rd);
                                } 

                                xnroot = _mm_sqrt_ps(xn);
                                ynroot = _mm_sqrt_ps(yn);
                                znroot = _mm_sqrt_ps(zn);
                                x0     = _mm_fmadd_ps(ynroot,znroot,_mm_add_ps(ynroot,znroot));
                                lamda  = _mm_mul_ps(xnroot,x0);
                                sigma  = _mm_div_ps(_mm_add_ps(sigma,pow4),
                                                       _mm_mul_ps(znroot,_mm_add_ps(zn,lamda)));
                                pow4   = _mm_mul_ps(pow4,c8);
                                xn     = _mm_mul_ps(_mm_add_ps(xn,lamda),c8);
                                yn     = _mm_mul_ps(_mm_add_ps(yn,lamda),c8);
                                zn     = _mm_mul_ps(_mm_add_ps(zn,lamda),c8);
                         }
                 }


                    
	          
	          
                  
	          
                   __m128 gms::radiolocation::rd_xmm4r4_a(const float * __restrict __ATTR_ALIGN__(16) px,
                                       const float * __restrict __ATTR_ALIGN__(16) py,
                                       const float * __restrict __ATTR_ALIGN__(16) pz,
                                       const float * __restrict __ATTR_ALIGN__(16) perrtot) {

                          register __m128 x       = _mm_load_ps(&px[0]);
                          register __m128 y       = _mm_load_ps(&py[0]);
                          register __m128 z       = _mm_load_ps(&pz[0]);
                          register __m128 errtot  = _mm_load_ps(&perrtot[0]);
                          const register __m128 _3 = _mm_set1_ps(3.0f);
                          const register __m128 _1 = _mm_set1_ps(1.0f);
                          const register __m128 c1 = _mm_set1_ps(-0.214285714285714285714285714286f);
                          const register __m128 c2 = _mm_set1_ps(0.166666666666666666666666666667f);
                          const register __m128 c3 = _mm_set1_ps(-0.409090909090909090909090909091f);
                          const register __m128 c4 = _mm_set1_ps(0.115384615384615384615384615385f);
                          const register __m128 c5 = _mm_set1_ps(6.0f);
                          const register __m128 c6 = _mm_set1_ps(1.5f);
                          const register __m128 c7 = _mm_set1_ps(0.2f);
                          const register __m128 c8 = _mm_set1_ps(0.25f);
                          register __m128 rd,xn,yn,zn,epslon,sigma,pow4,mu;
                          register __m128 xndev,yndev,zndev,ea,eb,ec,ed,ef;
                          register __m128 s1,s2,xnroot,ynroot,znroot,lamda;
                          register __m128 x0,x1,x2,x3,x4,x5;

                          xn    = x;
                          yn    = y;
                          zn    = z;
                          sigma = _mm_setzero_ps();
                          pow4  = _1; 
                          while(true) {
                                mu    = _mm_mul_ps(c7,_mm_fmadd_ps(zn,_3,
                                                          _mm_add_ps(xn,yn));
                                xndev = _mm_div_ps(_mm_sub_ps(mu,xn),mu);
                                yndev = _mm_div_ps(_mm_sub_ps(mu,yn),mu);
                                zndev = _mm_div_ps(_mm_sub_ps(mu,zn),mu);
                                epslon= _mm_abs_ps(xndev);
                                epslon= _mm_max_ps(epslon,_mm_abs_ps(yndev));
                                epslon= _mm_max_ps(epslon,_mm_abs_ps(zndev));

                                if(_mm_cmp_mask_ps(epslon,errtot,_CMP_LT_OQ)) {
                                   ea = _mm_mul_ps(xndev,yndev);
                                   eb = _mm_mul_ps(zndev,zndev);
                                   ec = _mm_sub_ps(ea,eb);
                                   ed = _mm_sub_ps(ea,_mm_mul_ps(c5,eb));
                                   ef = _mm_add_ps(ed,_mm_add_ps(ec,ec));
                                   x0 = _mm_fmadd_ps(c3,c8,c1);
                                   x1 = _mm_sub_ps(ed,_mm_sub_ps(c6,c4));
                                   x2 = _mm_mul_ps(zndev,ef);
                                   s1 = _mm_mul_ps(ed,_mm_mul_ps(x0,
                                                               _mm_mul_ps(x1,x2)));
                                   x3 = _mm_fmadd_ps(c3,ec,_mm_mul_ps(zndev,
                                                               _mm_mul_ps(c4,ea)));
                                   x4 = _mm_fmadd_ps(x3,zndev,_mm_mul_ps(ef,c2));
                                   s2 = _mm_mul_ps(zndev,x4);
                                   x0 = _mm_fmadd_ps(_3,sigma,pow4);
                                   x1 = _mm_add_ps(_1,_mm_add_ps(s1,s2));
                                   x2 = _mm_mul_ps(mu,_mm_sqrt_ps(mu));
                                   rd = _mm_div_ps(_mm_mul_ps(x0,x1),x2);
                                   return (rd);
                                } 

                                xnroot = _mm_sqrt_ps(xn);
                                ynroot = _mm_sqrt_ps(yn);
                                znroot = _mm_sqrt_ps(zn);
                                x0     = _mm_fmadd_ps(ynroot,znroot,_mm_add_ps(ynroot,znroot));
                                lamda  = _mm_mul_ps(xnroot,x0);
                                sigma  = _mm_div_ps(_mm_add_ps(sigma,pow4),
                                                       _mm_mul_ps(znroot,_mm_add_ps(zn,lamda)));
                                pow4   = _mm_mul_ps(pow4,c8);
                                xn     = _mm_mul_ps(_mm_add_ps(xn,lamda),c8);
                                yn     = _mm_mul_ps(_mm_add_ps(yn,lamda),c8);
                                zn     = _mm_mul_ps(_mm_add_ps(zn,lamda),c8);
                         }
                 }


                      
	          
	          
                  
	         
                   __m128 gms::radiolocation::rd_xmm4r4_u(const float * __restrict  px,
                                       const float * __restrict  py,
                                       const float * __restrict  pz,
                                       const float * __restrict  perrtot) {

                          register __m128 x       = _mm_loadu_ps(&px[0]);
                          register __m128 y       = _mm_loadu_ps(&py[0]);
                          register __m128 z       = _mm_loadu_ps(&pz[0]);
                          register __m128 errtot  = _mm_loadu_ps(&perrtot[0]);
                          const register __m128 _3 = _mm_set1_ps(3.0f);
                          const register __m128 _1 = _mm_set1_ps(1.0f);
                          const register __m128 c1 = _mm_set1_ps(-0.214285714285714285714285714286f);
                          const register __m128 c2 = _mm_set1_ps(0.166666666666666666666666666667f);
                          const register __m128 c3 = _mm_set1_ps(-0.409090909090909090909090909091f);
                          const register __m128 c4 = _mm_set1_ps(0.115384615384615384615384615385f);
                          const register __m128 c5 = _mm_set1_ps(6.0f);
                          const register __m128 c6 = _mm_set1_ps(1.5f);
                          const register __m128 c7 = _mm_set1_ps(0.2f);
                          const register __m128 c8 = _mm_set1_ps(0.25f);
                          register __m128 rd,xn,yn,zn,epslon,sigma,pow4,mu;
                          register __m128 xndev,yndev,zndev,ea,eb,ec,ed,ef;
                          register __m128 s1,s2,xnroot,ynroot,znroot,lamda;
                          register __m128 x0,x1,x2,x3,x4,x5;

                          xn    = x;
                          yn    = y;
                          zn    = z;
                          sigma = _mm_setzero_ps();
                          pow4  = _1; 
                          while(true) {
                                mu    = _mm_mul_ps(c7,_mm_fmadd_ps(zn,_3,
                                                          _mm_add_ps(xn,yn));
                                xndev = _mm_div_ps(_mm_sub_ps(mu,xn),mu);
                                yndev = _mm_div_ps(_mm_sub_ps(mu,yn),mu);
                                zndev = _mm_div_ps(_mm_sub_ps(mu,zn),mu);
                                epslon= _mm_abs_ps(xndev);
                                epslon= _mm_max_ps(epslon,_mm_abs_ps(yndev));
                                epslon= _mm_max_ps(epslon,_mm_abs_ps(zndev));

                                if(_mm_cmp_mask_ps(epslon,errtot,_CMP_LT_OQ)) {
                                   ea = _mm_mul_ps(xndev,yndev);
                                   eb = _mm_mul_ps(zndev,zndev);
                                   ec = _mm_sub_ps(ea,eb);
                                   ed = _mm_sub_ps(ea,_mm_mul_ps(c5,eb));
                                   ef = _mm_add_ps(ed,_mm_add_ps(ec,ec));
                                   x0 = _mm_fmadd_ps(c3,c8,c1);
                                   x1 = _mm_sub_ps(ed,_mm_sub_ps(c6,c4));
                                   x2 = _mm_mul_ps(zndev,ef);
                                   s1 = _mm_mul_ps(ed,_mm_mul_ps(x0,
                                                               _mm_mul_ps(x1,x2)));
                                   x3 = _mm_fmadd_ps(c3,ec,_mm_mul_ps(zndev,
                                                               _mm_mul_ps(c4,ea)));
                                   x4 = _mm_fmadd_ps(x3,zndev,_mm_mul_ps(ef,c2));
                                   s2 = _mm_mul_ps(zndev,x4);
                                   x0 = _mm_fmadd_ps(_3,sigma,pow4);
                                   x1 = _mm_add_ps(_1,_mm_add_ps(s1,s2));
                                   x2 = _mm_mul_ps(mu,_mm_sqrt_ps(mu));
                                   rd = _mm_div_ps(_mm_mul_ps(x0,x1),x2);
                                   return (rd);
                                } 

                                xnroot = _mm_sqrt_ps(xn);
                                ynroot = _mm_sqrt_ps(yn);
                                znroot = _mm_sqrt_ps(zn);
                                x0     = _mm_fmadd_ps(ynroot,znroot,_mm_add_ps(ynroot,znroot));
                                lamda  = _mm_mul_ps(xnroot,x0);
                                sigma  = _mm_div_ps(_mm_add_ps(sigma,pow4),
                                                       _mm_mul_ps(znroot,_mm_add_ps(zn,lamda)));
                                pow4   = _mm_mul_ps(pow4,c8);
                                xn     = _mm_mul_ps(_mm_add_ps(xn,lamda),c8);
                                yn     = _mm_mul_ps(_mm_add_ps(yn,lamda),c8);
                                zn     = _mm_mul_ps(_mm_add_ps(zn,lamda),c8);
                         }
                 }


/*							fresnl.c
 *
 *	Fresnel integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, S, C;
 * void fresnl();
 *
 * fresnl( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
 *
 * The integrals are evaluated by a power series for x < 1.
 * For x >= 1 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
 * S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
 *
 *
 *
 * ACCURACY:
 *
 *  Relative error.
 *
 * Arithmetic  function   domain     # trials      peak         rms
 *   IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
 *   IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
 *   DEC        S(x)      0, 10        6000       2.2e-16     3.9e-17
 *   DEC        C(x)      0, 10        5000       2.3e-16     3.9e-17
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/




                   
	          
	          
                  
	          
                   void gms::radiolocation::fresnel_xmm4r4(const __m128 xxa,
                                        __m128 * __restrict ssa,
                                        __m128 * __restrict cca) {

                        using namespace gms::math;
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,cc,ss,c,s,t,u,t0,t1;
                        register __m128 x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                           volatile __m128 prefcn = _mm_prefetch((const char*)&cn[0],_MM_HINT_T0);
                           volatile __m128 prefcd = _mm_prefetch((const char*)&cd[0],_MM_HINT_T0);
                           t = _mm_mul_ps(x2,x2);
                           acc1 = sn[0]; 
                           acc2 = _mm_add_ps(t,sd[0]);
                           acc3 = cn[0];
                           acc4 = cd[0];
                           acc1 = _mm_fmadd_ps(acc1,t,sn[1]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[1]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[2]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[2]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[3]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[3]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[5]);
                           t0   = _mm_div_ps(acc1,acc2);
                           ss   = _mm_mul_ps(_mm_mul_ps(x,x2),t0);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[1]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[1]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[2]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[2]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[3]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[3]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[5]);
                           t1   = _mm_div_ps(acc3,acc4);
                           cc   = _mm_mul_ps(x,t1);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          cc = hlf;
                          ss = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        c    = _mm_cos_ps(t);
                        s    = _mm_sin_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t0   = _mm_fmsub_ps(f,s,_mm_mul_ps(g,c));
                        cc   = _mm_add_ps(hlf,_mm_div_ps(t0,t));
                        t1   = _mm_fmadd_ps(f,c,_mm_mul_ps(g,s));
                        ss   = _mm_sub_ps(hlf,_mm_div_ps(t1,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         cc = negate_xmm4r4(cc);
                         ss = negate_xmm4r4(ss);
                     }
                     
                     *cca = cc;
                     *ssa = ss;
              }


                   
	          
	          
                  
	        
                   void gms::radiolocation::fresnel_xmm4r4_a(const float * __restrict __ATTR_ALIGN__(16) pxxa,
                                          float * __restrict __ATTR_ALIGN__(16) ssa,
                                          float * __restrict __ATTR_ALIGN__(16) cca) {

                        using namespace gms::math;
                        register __m128 xxa = _mm_load_ps(&pxxa[0]);
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,cc,ss,c,s,t,u,t0,t1;
                        register __m128 x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                           volatile __m128 prefcn = _mm_prefetch((const char*)&cn[0],_MM_HINT_T0);
                           volatile __m128 prefcd = _mm_prefetch((const char*)&cd[0],_MM_HINT_T0);
                           t = _mm_mul_ps(x2,x2);
                           acc1 = sn[0]; 
                           acc2 = _mm_add_ps(t,sd[0]);
                           acc3 = cn[0];
                           acc4 = cd[0];
                           acc1 = _mm_fmadd_ps(acc1,t,sn[1]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[1]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[2]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[2]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[3]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[3]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[5]);
                           t0   = _mm_div_ps(acc1,acc2);
                           ss   = _mm_mul_ps(_mm_mul_ps(x,x2),t0);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[1]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[1]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[2]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[2]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[3]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[3]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[5]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[6]);
                           t1   = _mm_div_ps(acc3,acc4);
                           cc   = _mm_mul_ps(x,t1);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          cc = hlf;
                          ss = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        c    = _mm_cos_ps(t);
                        s    = _mm_sin_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t0   = _mm_fmsub_ps(f,s,_mm_mul_ps(g,c));
                        cc   = _mm_add_ps(hlf,_mm_div_ps(t0,t));
                        t1   = _mm_fmadd_ps(f,c,_mm_mul_ps(g,s));
                        ss   = _mm_sub_ps(hlf,_mm_div_ps(t1,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         cc = negate_xmm4r4(cc);
                         ss = negate_xmm4r4(ss);
                     }
                     
                     _mm_store_ps(&cca[0] ,cc);
                     _mm_store_ps(&ssa[0] ,ss);
              }


                   
	          
	          
                  
	           
                   void gms::radiolocation::fresnel_xmm4r4_u(const float * __restrict  pxxa,
                                          float * __restrict  ssa,
                                          float * __restrict  cca) {

                        using namespace gms::math;
                        register __m128 xxa = _mm_loadu_ps(&pxxa[0]);
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,cc,ss,c,s,t,u,t0,t1;
                        register __m128 x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                           volatile __m128 prefcn = _mm_prefetch((const char*)&cn[0],_MM_HINT_T0);
                           volatile __m128 prefcd = _mm_prefetch((const char*)&cd[0],_MM_HINT_T0);
                           t = _mm_mul_ps(x2,x2);
                           acc1 = sn[0]; 
                           acc2 = _mm_add_ps(t,sd[0]);
                           acc3 = cn[0];
                           acc4 = cd[0];
                           acc1 = _mm_fmadd_ps(acc1,t,sn[1]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[1]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[2]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[2]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[3]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[3]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[5]);
                           t0   = _mm_div_ps(acc1,acc2);
                           ss   = _mm_mul_ps(_mm_mul_ps(x,x2),t0);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[1]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[1]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[2]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[2]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[3]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[3]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[5]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[6]);
                           t1   = _mm_div_ps(acc3,acc4);
                           cc   = _mm_mul_ps(x,t1);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          cc = hlf;
                          ss = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        c    = _mm_cos_ps(t);
                        s    = _mm_sin_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t0   = _mm_fmsub_ps(f,s,_mm_mul_ps(g,c));
                        cc   = _mm_add_ps(hlf,_mm_div_ps(t0,t));
                        t1   = _mm_fmadd_ps(f,c,_mm_mul_ps(g,s));
                        ss   = _mm_sub_ps(hlf,_mm_div_ps(t1,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         cc = negate_xmm4r4(cc);
                         ss = negate_xmm4r4(ss);
                     }
                     
                     _mm_storeu_ps(&cca[0] ,cc);
                     _mm_storeu_ps(&ssa[0] ,ss);
              }




                  /*
                           Same as above -- divided into Fresnel 'C' integral
                           and Fresnel 'S' integral.
                     */


                   
	          
	          
                  
	          
                   __m128 gms::radiolocation::fresnel_C_xmm4r4(const __m128 xxa) {
                                        
                        using namespace gms::math;
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,cc,c,t,u,t0,t1;
                        register __m128 cca,x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefcn = _mm_prefetch((const char*)&cn[0],_MM_HINT_T0);
                           volatile __m128 prefcd = _mm_prefetch((const char*)&cd[0],_MM_HINT_T0);
                          
                           t = _mm_mul_ps(x2,x2);
                           acc3 = cn[0];
                           acc4 = cd[0];
                           acc3 = _mm_fmadd_ps(acc3,t,cn[1]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[1]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[2]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[2]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[3]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[3]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[5]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[6]);
                           t1   = _mm_div_ps(acc3,acc4);
                           cc   = _mm_mul_ps(x,t1);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          cc = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        c    = _mm_cos_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t0   = _mm_fmsub_ps(f,s,_mm_mul_ps(g,c));
                        cc   = _mm_add_ps(hlf,_mm_div_ps(t0,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         cc = negate_xmm4r4(cc);
                     }
                     
                     cca = cc;
                     return (cca);
              }


                   
	          
	          
                  
	          
                   __m128 gms::radiolocation::fresnel_C_xmm4r4_a(const float * __restrict __ATTR_ALIGN__(16) pxxa) {
                                        
                        using namespace gms::math;
                        register __m128 xxa = _mm_load_ps(&pxxa[0]);
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,cc,c,t,u,t0,t1;
                        register __m128 cca,x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefcn = _mm_prefetch((const char*)&cn[0],_MM_HINT_T0);
                           volatile __m128 prefcd = _mm_prefetch((const char*)&cd[0],_MM_HINT_T0);
                          
                           t = _mm_mul_ps(x2,x2);
                           
                           acc3 = cn[0];
                           acc4 = cd[0];
                           acc3 = _mm_fmadd_ps(acc3,t,cn[1]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[1]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[2]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[2]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[3]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[3]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[5]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[6]);
                           t1   = _mm_div_ps(acc3,acc4);
                           cc   = _mm_mul_ps(x,t1);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          cc = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        c    = _mm_cos_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t0   = _mm_fmsub_ps(f,s,_mm_mul_ps(g,c));
                        cc   = _mm_add_ps(hlf,_mm_div_ps(t0,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         cc = negate_xmm4r4(cc);
                     }
                     
                     cca = cc;
                     return (cca);
              }


                  
	          
	          
                  
	          
                   __m128 gms::radiolocation::fresnel_C_xmm4r4_u(const float * __restrict  pxxa) {
                                        
                        using namespace gms::math;
                        register __m128 xxa = _mm_loadu_ps(&pxxa[0]);
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,cc,c,t,u,t0,t1;
                        register __m128 cca,x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                          
                           t = _mm_mul_ps(x2,x2);
                           acc3 = cn[0];
                           acc4 = cd[0];
                           acc3 = _mm_fmadd_ps(acc3,t,cn[1]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[1]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[2]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[2]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[3]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[3]);
                           acc3 = _mm_fmadd_ps(acc3,t,cn[4]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[5]);
                           acc4 = _mm_fmadd_ps(acc4,t,cd[6]);
                           t1   = _mm_div_ps(acc3,acc4);
                           cc   = _mm_mul_ps(x,t1);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          cc = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        c    = _mm_cos_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t0   = _mm_fmsub_ps(f,s,_mm_mul_ps(g,c));
                        cc   = _mm_add_ps(hlf,_mm_div_ps(t0,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         cc = negate_xmm4r4(cc);
                     }
                     
                     cca = cc;
                     return (cca);
              }


                   
	          
	          
                  
	          
                   __m128 gms::radiolocation::fresnel_S_xmm4r4(const __m128 xxa) {
                                        
                        using namespace gms::math;
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,ss,s,t,u,t0,t1;
                        register __m128 ssa,x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                           
                           t = _mm_mul_ps(x2,x2);
                           acc1 = sn[0]; 
                           acc2 = _mm_add_ps(t,sd[0]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[1]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[1]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[2]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[2]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[3]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[3]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[5]);
                           t0   = _mm_div_ps(acc1,acc2);
                           ss   = _mm_mul_ps(_mm_mul_ps(x,x2),t0);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          ss = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        s    = _mm_sin_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t1   = _mm_fmadd_ps(f,c,_mm_mul_ps(g,s));
                        ss   = _mm_sub_ps(hlf,_mm_div_ps(t1,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         ss = negate_xmm4r4(ss);
                     }
                     
                     ssa = ss;
                     return (ssa);
              }


                    
	          
	          
                  
	           
                   __m128 gms::radiolocation::fresnel_S_xmm4r4_a(const float * __restrict __ATTR_ALIGN__(16) pxxa) {
                                        
                        using namespace gms::math;
                        register __m128 xxa = _mm_load_ps(&pxxa[0]);
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,ss,s,t,u,t0,t1;
                        register __m128 ssa,x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                           
                           t = _mm_mul_ps(x2,x2);
                           acc1 = sn[0]; 
                           acc2 = _mm_add_ps(t,sd[0]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[1]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[1]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[2]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[2]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[3]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[3]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[5]);
                           t0   = _mm_div_ps(acc1,acc2);
                           ss   = _mm_mul_ps(_mm_mul_ps(x,x2),t0);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          ss = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        s    = _mm_sin_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t1   = _mm_fmadd_ps(f,c,_mm_mul_ps(g,s));
                        ss   = _mm_sub_ps(hlf,_mm_div_ps(t1,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         ss = negate_xmm4r4(ss);
                     }
                     
                     ssa = ss;
                     return (ssa);
              }
  


                   
	          
	          
                  
	          
                   __m128 gms::radiolocation::fresnel_S_xmm4r4_u(const float * __restrict  pxxa) {
                                        
                        using namespace gms::math;
                        register __m128 xxa = _mm_loadu_ps(&pxxa[0]);
                        const __m128 c0   = _mm_set1_ps(2.5625f);
                        const __m128 c1   = _mm_set1_ps(36974.0f);
                        const __m128 hlf  = _mm_set1_ps(0.5f);
                        const __m128 _0   = _mm_setzero_ps();
                        const __m128 _1   = _mm_set1_ps(1.0f); 
                        const __m128 pi   = _mm_set1_ps(3.14159265358979323846264338328f);
                        const __m128 pio2 = _mm_set1_ps(1.57079632679489661923132169164f);
                        register __m128 f,g,ss,s,t,u,t0,t1;
                        register __m128 ssa,x,x2,acc1,acc2,acc3,acc4;
                       
                        x   = _mm_abs_ps(xxa);
                        x2  = _mm_mul_ps(x,x);
                        if(_mm_cmp_ps_mask(x,c0,_CMP_LT_OQ)) {
			   volatile __m128 prefsn = _mm_prefetch((const char*)&sn[0],_MM_HINT_T0);
                           volatile __m128 prefsd = _mm_prefetch((const char*)&sd[0],_MM_HINT_T0);
                           
                           t = _mm_mul_ps(x2,x2);
                           acc1 = sn[0]; 
                           acc2 = _mm_add_ps(t,sd[0]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[1]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[1]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[2]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[2]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[3]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[3]);
                           acc1 = _mm_fmadd_ps(acc1,t,sn[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[4]);
                           acc2 = _mm_fmadd_ps(acc2,t,sd[5]);
                           t0   = _mm_div_ps(acc1,acc2);
                           ss   = _mm_mul_ps(_mm_mul_ps(x,x2),t0);
                           goto done;
                        }

                       if(_mm_cmp_ps_mask(x,c1,_CMP_GT_OQ)) {
                          ss = hlf;
                          goto done;
                      }

                      /*		Asymptotic power series auxiliary functions
                       *		for large argument
                       */

                        volatile __m128 prefsn = _mm_prefetch((const char*)&fn[0],_MM_HINT_T0);
                        volatile __m128 prefsd = _mm_prefetch((const char*)&fd[0],_MM_HINT_T0);
                        volatile __m128 prefcn = _mm_prefetch((const char*)&gn[0],_MM_HINT_T0);
                        volatile __m128 prefcd = _mm_prefetch((const char*)&gd[0],_MM_HINT_T0);
                        t = _mm_mul_ps(pi,x2);
                        u = _mm_div_ps(_1,_mm_mul_ps(t,t));
                        acc1 = fn[0];
                        acc2 = _mm_add_ps(u,fd[0]);
                        acc3 = gn[0];
                        acc4 = _mm_add_ps(u,gd[0]);
                        t = _mm_div_ps(_1,t);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[1]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[1]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[2]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[2]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[3]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[3]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[4]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[4]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[5]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[5]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[6]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[6]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[7]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[7]);
                        acc1 = _mm_fmadd_ps(acc1,u,fn[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[8]);
                        acc2 = _mm_fmadd_ps(acc2,u,fd[9]);
                        t0   = _mm_div_ps(acc1,acc2);
                        f    = _mm_sub_ps(_1,_mm_mul_ps(u,t0));
                        acc3 = _mm_fmadd_ps(acc3,u,gn[1]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[1]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[2]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[2]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[3]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[3]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[4]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[4]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[5]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[5]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[6]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[6]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[7]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[7]);
                        acc3 = _mm_fmadd_ps(acc3,u,gn[8]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[8]);
                        acc3 = _mm_fmadd_ps(acc3,u,gd[9]);
                        acc4 = _mm_fmadd_ps(acc4,u,gd[10]);
                        t1   = _mm_div_ps(acc3,acc4);
                        g    = _mm_mul_ps(t,t1);
                        
                        t    = _mm_mul_ps(pio2,x2);
                        s    = _mm_sin_ps(t);
                        t    = _mm_mul_ps(pi,x);
                        t1   = _mm_fmadd_ps(f,c,_mm_mul_ps(g,s));
                        ss   = _mm_sub_ps(hlf,_mm_div_ps(t1,t));
done:
                     if(_mm_cmp_ps_mask(xxa,
                                     _mm_setzero_ps(),_CMP_LT_OQ)) {
                         ss = negate_xmm4r4(ss);
                     }
                     
                     ssa = ss;
                     return (ssa);
              }
  

                   





   
