
#ifndef __GMS_ROOT_FINDING_HPP__
#define __GMS_ROOT_FINDING_HPP__ 170420221521

/*
    Converted by Bernard Gingold, beniekg@gmail.com from Fortran root_module.f90 implementation
   
    Root solver methods for:
!
!  * Bracked interval
!  * Without derivatives
!
!### Author
!  * Jacob Williams
*/

namespace file_info {

     const unsigned int GMS_ROOT_FINDING_MAJOR = 1;
     const unsigned int GMS_ROOT_FINDING_MINOR = 1;
     const unsigned int GMS_ROOT_FINDING_MICRO = 0;
     const unsigned int GMS_ROOT_FINDING_FULLVER =
       1000U*GMS_ROOT_FINDING_MAJOR+100U*GMS_ROOT_FINDING_MINOR+
       10U*GMS_ROOT_FINDING_MICRO;
     const char * const GMS_ROOT_FINDING_CREATION_DATE = "17-04-2022 15:21 +00200 (SUN 17 APR 2022 GMT+2)";
     const char * const GMS_ROOT_FINDING_BUILD_DATE    = __DATE__ " " __TIME__;
     const char * const GMS_ROOT_FINDING_SYNOPSIS      = "Function root finding algorithms."

}


#include <cstdint>
#include <limits>
#include <cmath>
#include "GMS_config.h"


#define FTOL8 0.0
#define FTOL4 0.0f
#define RTOL8 0.000001
#define RTOL4 0.000001f
#define ATOL8 0.000000000001
#define ATOL4 0.000000000001f
#define MAXITER 2000


namespace  gms {


           namespace math {


/*
                          Find a zero of the function \( f(x) \) in the given interval
!  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
!  where \( \epsilon \) is the relative machine precision defined as
!  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
!
!### References
!  * R. P. Brent, "[An algorithm with guaranteed convergence for
!    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
!    The Computer Journal, Vol 14, No. 4., 1971.
!  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
!    Prentice-Hall, Inc., 1973.
!
!### See also
!  * [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib
*/


                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void brent(float(*f)(float x),
		                const float ax,
				const float bx,
				const float fax,
				const float fbx,
				float & xzero,
				float & fzero,
				int32_t & iflag) {

                      float a,b,c,d,e,f,fa,fb,fc,tol1,xm,
		            p,q,r,s;
		      iflag = 0;
		      tol1  = std::numeric_limits<float>::epsilon()+1.0f;
		      a     = ax;
		      b     = bx;
		      fa    = fax;
		      fb    = fbx;
		      c     = a;
		      fc    = fa;
		      d     = b-a;
		      e     = d;

		      for(int32_t i = i; i != MAXITER; ++i) {
                          if(std::abs(fc)<std::abs(fb)) {
                             a = b;
			     b = c;
			     c = a;
			     fa=fb;
			     fb=fc;
			     fc=fa;
			  }
			  tol1 = 2.0f*std::numeric_limits<float>::epsilon()+
			         0.5f*RTOL4;
			  xm   = 0.5f*(c-b);
			  if(std::abs(xm<=tol1) break;
                          if(std::abs(e)>=tol1 && (std::abs(fa)>std::abs(fb))) {
                             s = fb/fa;
			     if(a!=c) {
                                //! inverse quadratic interpolation
                                q=fa/fc;
                                r=fb/fc;
                                p=s*(2.0f*xm*q*(q-r)-(b-a)*(r-1.0f));
                                q=(q-1.0f)*(r-1.0f)*(s-1.0f);
			     }
			     else {
                                //! linear interpolation
                                p=2.0f*xm*s;
                                q=1.0f-s;
			     }
			     if(p<=0.0f) {
                                p = -p;
			     }
			     else {
                                p = -q;
			     }
			     s = e;
			     e = d;
			     if((2.0f*p)>=(3.0f*xm*q-std::abs(tol1*q)) ||
			        (p>=std::abs(0.5f*s*q))) {
                                 d = xm;
				 e = d;
			     }
			     else {
                                 d = p/q;
			     }
			  }
			  else {
                             d = xm;
			     e = d;
			  }
			  a = b;
			  fa=fb;
			  if(std::abs(d)<=tol1) {
                             if(xm<=0.0f) {
                                b -= tol1;
			     }
			     else {
                                b += tol1;
			     }
			  }
			  else {
                            b += d;
			  }
			  fb = f(b);
			  if(std::abs(fb)<=FTOL4) break; // absolute convergence in f
			  if((fb*(fc/std::abs(fc)))>0.0f) {
                             c = a;
			     fc= fa;
			     d = b-a;
			     e = d;
			  }
			  if(i==MAXITER) iflag = -2; // max iterations reached!
		      }
		      xzero = b;
		      fzero = fb;
		 }


		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void brent(double(*f)(double x),
		                const double ax,
				const double bx,
				const double fax,
				const double fbx,
				double & xzero,
				double & fzero,
				int32_t & iflag) {

                      double a,b,c,d,e,f,fa,fb,fc,tol1,xm,
		            p,q,r,s;
		      iflag = 0;
		      tol1  = std::numeric_limits<double>::epsilon()+1.0;
		      a     = ax;
		      b     = bx;
		      fa    = fax;
		      fb    = fbx;
		      c     = a;
		      fc    = fa;
		      d     = b-a;
		      e     = d;

		      for(int32_t i = i; i != MAXITER; ++i) {
                          if(std::abs(fc)<std::abs(fb)) {
                             a = b;
			     b = c;
			     c = a;
			     fa=fb;
			     fb=fc;
			     fc=fa;
			  }
			  tol1 = 2.0*std::numeric_limits<double>::epsilon()+
			         0.5*RTOL8;
			  xm   = 0.5*(c-b);
			  if(std::abs(xm<=tol1) break;
                          if(std::abs(e)>=tol1 && (std::abs(fa)>std::abs(fb))) {
                             s = fb/fa;
			     if(a!=c) {
                                //! inverse quadratic interpolation
                                q=fa/fc;
                                r=fb/fc;
                                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                                q=(q-1.0)*(r-1.0)*(s-1.0);
			     }
			     else {
                                //! linear interpolation
                                p=2.0*xm*s;
                                q=1.0-s;
			     }
			     if(p<=0.0) {
                                p = -p;
			     }
			     else {
                                p = -q;
			     }
			     s = e;
			     e = d;
			     if((2.0*p)>=(3.0f*xm*q-std::abs(tol1*q)) ||
			        (p>=std::abs(0.5*s*q))) {
                                 d = xm;
				 e = d;
			     }
			     else {
                                 d = p/q;
			     }
			  }
			  else {
                             d = xm;
			     e = d;
			  }
			  a = b;
			  fa=fb;
			  if(std::abs(d)<=tol1) {
                             if(xm<=0.0) {
                                b -= tol1;
			     }
			     else {
                                b += tol1;
			     }
			  }
			  else {
                            b += d;
			  }
			  fb = f(b);
			  if(std::abs(fb)<=FTOL8) break; // absolute convergence in f
			  if((fb*(fc/std::abs(fc)))>0.0) {
                             c = a;
			     fc= fa;
			     d = b-a;
			     e = d;
			  }
			  if(i==MAXITER) iflag = -2; // max iterations reached!
		      }
		      xzero = b;
		      fzero = fb;
		 }



     } //math

} //gms







#endif /*__GMS_ROOT_FINDING_HPP__*/
