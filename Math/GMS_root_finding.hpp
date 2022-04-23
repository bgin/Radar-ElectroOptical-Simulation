
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


	              


	            // Local helper functions
		    //Determines convergence in x based on if the reltol or abstol is satisfied.
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     __ATTR_PURE__
		     static
		     inline
		     bool converged(const float a,
		                    const float b) {
                          
                          float d = 0.0f;
			  bool  converged = false;
			  d = std::abs(b-a);
			  if(d<=ATOL4) {
                             converged = true;
			  }
			  else {
                                 if(a!=0.0f) {
                                    converged = (d/std::abs(a))<=RTOL4;
				 }
				 else {
                                    converged = false;
				 }
			  }
			  return (converged);
		     }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
		     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     __ATTR_PURE__
		     static
		     inline
		     bool converged(const double a,
		                    const double b) {
                          
                          double d = 0.0;
			  bool  converged = false;
			  d = std::abs(b-a);
			  if(d<=ATOL8) {
                             converged = true;
			  }
			  else {
                                 if(a!=0.0) {
                                    converged = (d/std::abs(a))<=RTOL8;
				 }
				 else {
                                    converged = false;
				 }
			  }
			  return (converged);
		     }


		     //Given two points with two function evaluations, choose the best one
                     //!  (the one closest to the root).
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
                     void choose_best(const float x1,
		                      const float x2,
				      const float f1,
				      const float f2,
				      float & xbest,
				      float & fbest) {

			   if(std::abs(f1)<std::abs(f2)) {
                              xbest = x1;
			      fbest = f1;
			   }
			   else {
                              xbest = x2;
			      fbest = f2;
			   }
		    }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
                     void choose_best(const double x1,
		                      const double x2,
				      const double f1,
				      const double f2,
				      double & xbest,
				      double & fbest) {

			   if(std::abs(f1)<std::abs(f2)) {
                              xbest = x1;
			      fbest = f1;
			   }
			   else {
                              xbest = x2;
			      fbest = f2;
			   }
		    }


		    // Bisection step.
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     __ATTR_PURE__
		     static
		     inline
		     float bisect(const float x1,
		                  const float x2) {

			 float x3 = 0.0f;
			 x3 = (x1+x2)*0.5f;
			 return (x3);
		    }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     __ATTR_PURE__
		     static
		     inline
		     double bisect(const double x1,
		                   const double x2) {

			 double x3 = 0.0;
			 x3 = (x1+x2)*0.5;
			 return (x3);
		    }


		   //Regula Falsi step.
                   //!  With a protection to fall back to bisection if:
                   //!
                   //!   * the computed point is outside the original interval ([ax,bx]).
                   //!   * f2 == f1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                    
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    __ATTR_PURE__
		    static
		    inline
		    float regula_falsi_step(const float x1,
		                            const float x2,
					    const float f1,
					    const float f2,
					    const float ax,
					    const float bx) {

                         float delta = 0.0f;
			 float x3    = 0.0f;
			 delta       = f2-f1;
			 if(delta!=0.0f) {
                            //   ! intersection with x-axis of line connecting the two points:
			    x3 = x1-(f1/delta)*(x2-x1);
			    if(x3>ax && x3<bx) return (std::numeric_limits<float>::quiet_Nan());
			 }
			 x3 = bisect(x1,x2);
			 return (x3);
		   }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                    
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    __ATTR_PURE__
		    static
		    inline
		    double regula_falsi_step(const double x1,
		                             const double x2,
					     const double f1,
					     const double f2,
					     const double ax,
					     const double bx) {

                         double delta = 0.0;
			 double x3    = 0.0;
			 delta       = f2-f1;
			 if(delta!=0.0) {
                            //   ! intersection with x-axis of line connecting the two points:
			    x3 = x1-(f1/delta)*(x2-x1);
			    if(x3>ax && x3<bx) return (std::numeric_limits<double>::quiet_Nan());
			 }
			 x3 = bisect(x1,x2);
			 return (x3);
		   }


		    //Secent step.
                    //!  With a protection to fall back to bisection if:
                    //!
                    //!   * the computed point is outside the original interval ([ax,bx]).
                    //!   * f2 == f1
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                    
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    __ATTR_PURE__
		    static
		    inline
		    float secant(const float x1,
		                 const float x2,
				 const float f1,
				 const float f2,
				 const float ax,
				 const float bx) {

                        float x3 = 0.0f;
			if(f2==f1) {
                           x3 = bisect(x1,x2);
			}
			else {
                             x3 = x2-f2/(f2-f1)/(x2-x1);
			     if(x3<ax || x3>bx) x3 = bisect(x1,x2);
			}
			return (x3);
		   }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                    
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    __ATTR_PURE__
		    static
		    inline
		    double secant(const double x1,
		                  const double x2,
				  const double f1,
				  const double f2,
				  const double ax,
				  const double bx) {

                        double x3 = 0.0;
			if(f2==f1) {
                           x3 = bisect(x1,x2);
			}
			else {
                             x3 = x2-f2/((f2-f1)/(x2-x1));
			     if(x3<ax || x3>bx) x3 = bisect(x1,x2);
			}
			return (x3);
		   }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                    
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    __ATTR_PURE__
		    static
		    inline
		    bool solution(const float x,
		                  const float f,
				  const float ftol,
				  float & xzero,
				  float & fzero) {

                       bool result = false;
		       if(std::abs(f)<=ftol) {
                          xzero = x;
			  fzero = f;
			  result = true;
		       }
		       else {
                          result = false;
		       }
		       return (result);
		   }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endifcc optimization_level 3
#endif
                    
		    __ATTR_ALWAYS_INLINE
		    __ATTR_HOT__
		    __ATTR_ALIGN__(32)
		    __ATTR_PURE__
		    static
		    inline
		    bool solution(const double x,
		                  const double f,
				  const double ftol,
				  double & xzero,
				  double & fzero) {

                       bool result = false;
		       if(std::abs(f)<=ftol) {
                          xzero = x;
			  fzero = f;
			  result = true;
		       }
		       else {
                          result = false;
		       }
		       return (result);
		   }		   


		   
		   
		   

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
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
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

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
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

		 
/*
 BlendTF blended method of trisection and false position methods.
!
!### Reference
!  * E Badr, S Almotairi, A El Ghamry,
!    "A Comparative Study among New Hybrid Root Finding
!    Algorithms and Traditional Methods", Mathematics 2021, 9, 1306.
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void blendtf(float (*f)(float x),
		                  const float ax,
		                  const float bx,
				  const float fax,
				  const float fbx,
				  float & xzero,
				  float & fzero,
				  int32_t & iflag) {

			 constexpr float third = 0.33333333333333333333333333f;
                         float a1,a2,b1,b2,fa,fb,a,b,xt1,xt2,xf,x,fx,fxt2,
                               fxf,xprev,fxt1,fxprev,fa1,fa2,fb1,fb2;
			 iflag = 0;
			 a     = ax;
			 b     = bx;
			 fa    = fax;
			 fb    = fbx;
			 a1    = a;
			 a2    = a;
			 b1    = b;
			 b2    = b;
			 fa1   = fa;
			 fa2   = fa;
			 fb1   = fb;
			 fb2   = fb;
			 xprev = std::numeric_limits<float>::max();
			 fxprev= std::numeric_limits<float>::max();

			 for(int32_t i = 1; i < MAXITER; ++i) {

			      if(fa==fb) {
                                 iflag = -3;
				 return;
			      }
			     xt1  = (b + 2.0f * a) * third;
                             xt2  = (2.0f * b + a) * third;
                             xf   = a - (fa*(b-a))/(fb-fa);
                             x    = xt1;
			     fxt1 = f(xt1);
			     if(solution(xt1,fxt1,FTOL4,xzero,fzero)) return;
			     fxt2 = f(xt2);
			     if(solution(xt2,fxt2,FTOL4,xzero,fzero)) return;
			     fxf  = f(xf);
			     if(solution(xf,fxf,FTOL4,xzero,fzero)) return;
			     fx   = fxt1;

			     if(std::abs(fxf)<std::abs(fxt1)) {
                                x = xf;
				fx= fxf;
			     }
			     else if(std::abs(fxt2)<std::abs(fxt1)) {
                                x  = xt2;
				fx = fxt2;
			     }

			     if(converged(a,b) || i==MAXITER) {
                                choose_best(x,xprev,fx,fxprev,xzero,fzero);
				if(i==MAXITER) { iflag = -2; break;}
			     }
			     xprev = x;
			     fxprev= fx;

			     if((fx * fxt1) < 0.0f) {
                                 b1 = xt1;
				 fb1= fxt1;
			     }
			     else if((fxt1 * fxt2) < 0.0f) {
                                 a1 = xt1;
				 fa1= fxt1;
				 b1 = xt2;
				 fb1= fxt2;
			     }
			     else {
                                 a1 = xt2;
				 fa1= fxt2;
			     }

			     if(fa * fxf < 0.0f) {
                                b2 = xf;
				fb2= fxf;
			     }
			     else {
                                a2 = xf;
				fa2= fxf;
			     }

			     if(a1>a2) {
                                a = a1;
				fa= fa1;
			     }
			     else {
                                a = a2;
				fa= fa2;
			     }

			     if(b1<b2) {
                                b = b1;
				fb= fb1;
			     }
			     else {
                                b = b2;
				fb= fb2;
			     }
			 }
		   }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void blendtf(double (*f)(double x),
		                  const double ax,
		                  const double bx,
				  const double fax,
				  const double fbx,
				  double & xzero,
				  double & fzero,
				  int32_t & iflag) {

			 constexpr double third = 0.33333333333333333333333333;
                         double a1,a2,b1,b2,fa,fb,a,b,xt1,xt2,xf,x,fx,fxt2,
                                fxf,xprev,fxt1,fxprev,fa1,fa2,fb1,fb2;
			 iflag = 0;
			 a     = ax;
			 b     = bx;
			 fa    = fax;
			 fb    = fbx;
			 a1    = a;
			 a2    = a;
			 b1    = b;
			 b2    = b;
			 fa1   = fa;
			 fa2   = fa;
			 fb1   = fb;
			 fb2   = fb;
			 xprev = std::numeric_limits<double>::max();
			 fxprev= std::numeric_limits<double>::max();

			 for(int32_t i = 1; i < MAXITER; ++i) {

			      if(fa==fb) {
                                 iflag = -3;
				 return;
			      }
			     xt1  = (b + 2.0 * a) * third;
                             xt2  = (2.0 * b + a) * third;
                             xf   = a - (fa*(b-a))/(fb-fa);
                             x    = xt1;
			     fxt1 = f(xt1);
			     if(solution(xt1,fxt1,FTOL8,xzero,fzero)) return;
			     fxt2 = f(xt2);
			     if(solution(xt2,fxt2,FTOL8,xzero,fzero)) return;
			     fxf  = f(xf);
			     if(solution(xf,fxf,FTOL8,xzero,fzero)) return;
			     fx   = fxt1;

			     if(std::abs(fxf)<std::abs(fxt1)) {
                                x = xf;
				fx= fxf;
			     }
			     else if(std::abs(fxt2)<std::abs(fxt1)) {
                                x  = xt2;
				fx = fxt2;
			     }

			     if(converged(a,b) || i==MAXITER) {
                                choose_best(x,xprev,fx,fxprev,xzero,fzero);
				if(i==MAXITER) { iflag = -2; break;}
			     }
			     xprev = x;
			     fxprev= fx;

			     if((fx * fxt1) < 0.0) {
                                 b1 = xt1;
				 fb1= fxt1;
			     }
			     else if((fxt1 * fxt2) < 0.0) {
                                 a1 = xt1;
				 fa1= fxt1;
				 b1 = xt2;
				 fb1= fxt2;
			     }
			     else {
                                 a1 = xt2;
				 fa1= fxt2;
			     }

			     if(fa * fxf < 0.0) {
                                b2 = xf;
				fb2= fxf;
			     }
			     else {
                                a2 = xf;
				fa2= fxf;
			     }

			     if(a1>a2) {
                                a = a1;
				fa= fa1;
			     }
			     else {
                                a = a2;
				fa= fa2;
			     }

			     if(b1<b2) {
                                b = b1;
				fb= fb1;
			     }
			     else {
                                b = b2;
				fb= fb2;
			     }
			 }
		   }

/*
  Modified anderson-bjorck-king method. Same as [[anderson_bjorck]], but with
!  an extra initial bisection step.
!
!### See also
!  * Kroger & Torsten, "On-Line Trajectory Generation in Robotic Systems", 2010.
!    https://link.springer.com/content/pdf/bbm%3A978-3-642-05175-3%2F1.pdf
     real(wp),intent(in)    :: ax      !! left endpoint of initial interval
     real(wp),intent(in)    :: bx      !! right endpoint of initial interval
     real(wp),intent(in)    :: fax     !! `f(ax)`
     real(wp),intent(in)    :: fbx     !! `f(ax)`
     real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
     real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
     integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void anderson_bjorck_king(float(*f)(float x),
		                               const float ax,
		                               const float bx,
					       const float fax,
					       const float fbx,
					       float & xzero,
					       float & fzero,
					       int32_t & iflag) {

			    float x1,x2,x3,f1,f2,f3,g,f1tmp;
                            bool root_found;

			    iflag = 0;
			    x1    = ax;
			    x2    = bx;
			    f1    = fax;
			    f2    = fbx;

			    for(int32_t i = 1; i < MAXITER; ++i) {

			         x3 = bisect(x1,x2);
				 f3 = f(x3);
				 if(solution(x3,f3,FTOL4,xzero,fzero)) return;
				 if(f2*f3 < 0.0f) {
                                     x1 = x2;
                                     x2 = x3;
                                     f1 = f2;
                                     f2 = f3;
				 }
				 else {
                                     x2 = x3;
				     f2 = f3;
				 }
				 x3 = secant(x1,x2,f1,f2,ax,bx);
				 f3 = f(x3);
				 if(solution(x3,f3,FTOL4,xzero,fzero)) return;
				 if(f2*f3 < 0.0f) {
                                     x1 = x2;
                                     x2 = x3;
                                     f1 = f2;
                                     f2 = f3;
                                     f1tmp = f1;
				 }
				 else {
                                      //! zero lies between x1 and x3
                                     g = 1.0f - f3/f2;
                                     if (g<=0.0f) g = 0.5f;
                                     x2 = x3;
                                     f1tmp = f1;
                                     f1 = g*f1;
                                     f2 = f3;
				 }
				 root_found = converged(x1,x2);
				 if(root_found || i==MAXITER) {
                                    choose_best(x1,x2,f1tmp,f2,xzero,fzero);
				    if(!root_found) iflag = -2;
				    break;
				 }
			    }
		     }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void anderson_bjorck_king(double(*f)(double x),
		                               const double ax,
		                               const double bx,
					       const double fax,
					       const double fbx,
					       double & xzero,
					       double & fzero,
					       int32_t & iflag) {

			    double x1,x2,x3,f1,f2,f3,g,f1tmp;
                            bool root_found;

			    iflag = 0;
			    x1    = ax;
			    x2    = bx;
			    f1    = fax;
			    f2    = fbx;

			    for(int32_t i = 1; i < MAXITER; ++i) {

			         x3 = bisect(x1,x2);
				 f3 = f(x3);
				 if(solution(x3,f3,FTOL8,xzero,fzero)) return;
				 if(f2*f3 < 0.0) {
                                     x1 = x2;
                                     x2 = x3;
                                     f1 = f2;
                                     f2 = f3;
				 }
				 else {
                                     x2 = x3;
				     f2 = f3;
				 }
				 x3 = secant(x1,x2,f1,f2,ax,bx);
				 f3 = f(x3);
				 if(solution(x3,f3,FTOL8,xzero,fzero)) return;
				 if(f2*f3 < 0.0) {
                                     x1 = x2;
                                     x2 = x3;
                                     f1 = f2;
                                     f2 = f3;
                                     f1tmp = f1;
				 }
				 else {
                                      //! zero lies between x1 and x3
                                     g = 1.0 - f3/f2;
                                     if (g<=0.0f) g = 0.5;
                                     x2 = x3;
                                     f1tmp = f1;
                                     f1 = g*f1;
                                     f2 = f3;
				 }
				 root_found = converged(x1,x2);
				 if(root_found || i==MAXITER) {
                                    choose_best(x1,x2,f1tmp,f2,xzero,fzero);
				    if(!root_found) iflag = -2;
				    break;
				 }
			    }
		     }

#include <algorithm>		     
/*
  Zhang's method (with corrections from Stage).
!
!### Reference
!  * A. Zhang, "An Improvement to the Brent's Method",
!    International Journal of Experimental Algorithms (IJEA), Volume (2) : Issue (1) : 2011.
!    https://www.cscjournals.org/download/issuearchive/IJEA/Volume2/IJEA_V2_I1.pdf
!  * S. A. Stage, "Comments on An Improvement to the Brent's Method",
!    International Journal of Experimental Algorithms (IJEA), Volume (4) : Issue (1) : 2013.
!    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.740.923&rep=rep1&type=pdf
     real(wp),intent(in)  :: ax      !! left endpoint of initial interval
     real(wp),intent(in)  :: bx      !! right endpoint of initial interval
     real(wp),intent(in)  :: fax     !! `f(ax)`
     real(wp),intent(in)  :: fbx     !! `f(ax)`
     real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
     real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
     integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

*/	      
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void zhang(float(*f)(float x),
		                const float ax,
		                const float bx,
				const float fax,
				const float fbx,
				float & xzero,
				float & fzero,
				int32_t & iflag) {

                       float a,b,c,fa,fb,fc,s,fs;
		       iflag = 0;
                       a     = ax;
                       b     = bx;
                       fa    = fax;
                       fb    = fbx;

		       for(int32_t i = 1; i < MAXITER; ++i) {

		             c  = bisect(a,b);
                             fc = f(c);
                             if(solution(c,fc,FTOL4,xzero,fzero)) return;
                             if(fa!=fc && fb!=fc) {
                                  //! inverse quadratic interpolation
                                 s = a*fb*fc/((fa-fb)*(fa-fc)) + 
                                     b*fa*fc/((fb-fa)*(fb-fc)) + 
                                     c*fa*fb/((fc-fa)*(fc-fb));
                                 if(a<s && s<b)  {
                                    fs = f(s);
                                    if(abs(fs)<=FTOL4) {
                                       xzero = s;
                                       fzero = fs;
                                       return;
                                     }
		                 }		    
                                 else {
                                         //! s is not in (a,b)
                                         s = c; //! just use this (there are 3 options in the reference)
                                         fs = fc;
                                  }
			     }
                             else {
                                    //! secant
                                    if(fa*fc<0.0f) {   //! root in [a,c]
                                        s = secant(a,c,fa,fc,ax,bx);
				    }
                                    else {         //! root in [c,b]
                                        s = secant(c,b,fc,fb,ax,bx);
                                    }
                                    fs = f(s);
                                    if(solution(s,fs,FTOL4,xzero,fzero)) return;
                             }

                             if (c>s) {
                               //! ensures a <= c <= s <= b
                                 std::swap(s,c);
                                 std::swap(fs,fc);
                             }
                             if (fc*fs<0.0f){      //! root on [c,s]
                                  a = c;
                                  b = s;
                                  fa = fc;
                                  fb = fs;
			     }
                             else if(fa*fc<0.0f){  //! root on [a,c]
                                  b = c;
                                  fb = fc;
			     }
                             else {                        //! root on [s,b]
                                  a = s;
                                  fa = fs;
                            }

                            if(converged(a,b)) break;
                            if(i == MAXITER) iflag = -2; //! max iterations reached

		       }
		       choose_best(a,b,fa,fb,xzero,fzero);
		  }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void zhang(double(*f)(double x),
		                const double ax,
		                const double bx,
				const double fax,
				const double fbx,
				double & xzero,
				double & fzero,
				int32_t & iflag) {

                       double a,b,c,fa,fb,fc,s,fs;
		       iflag = 0;
                       a     = ax;
                       b     = bx;
                       fa    = fax;
                       fb    = fbx;

		       for(int32_t i = 1; i < MAXITER; ++i) {

		             c  = bisect(a,b);
                             fc = f(c);
                             if(solution(c,fc,FTOL8,xzero,fzero)) return;
                             if(fa!=fc && fb!=fc) {
                                  //! inverse quadratic interpolation
                                 s = a*fb*fc/((fa-fb)*(fa-fc)) + 
                                     b*fa*fc/((fb-fa)*(fb-fc)) + 
                                     c*fa*fb/((fc-fa)*(fc-fb));
                                 if(a<s && s<b)  {
                                    fs = f(s);
                                    if(abs(fs)<=FTOL8) {
                                       xzero = s;
                                       fzero = fs;
                                       return;
                                     }
		                 }		    
                                 else {
                                         //! s is not in (a,b)
                                         s = c; //! just use this (there are 3 options in the reference)
                                         fs = fc;
                                  }
			     }
                             else {
                                    //! secant
                                    if(fa*fc<0.0) {   //! root in [a,c]
                                        s = secant(a,c,fa,fc,ax,bx);
				    }
                                    else {         //! root in [c,b]
                                        s = secant(c,b,fc,fb,ax,bx);
                                    }
                                    fs = f(s);
                                    if(solution(s,fs,FTOL8,xzero,fzero)) return;
                             }

                             if (c>s) {
                               //! ensures a <= c <= s <= b
                                 std::swap(s,c);
                                 std::swap(fs,fc);
                             }
                             if (fc*fs<0.0){      //! root on [c,s]
                                  a = c;
                                  b = s;
                                  fa = fc;
                                  fb = fs;
			     }
                             else if(fa*fc<0.0){  //! root on [a,c]
                                  b = c;
                                  fb = fc;
			     }
                             else {                        //! root on [s,b]
                                  a = s;
                                  fa = fs;
                            }

                            if(converged(a,b)) break;
                            if(i == MAXITER) iflag = -2; //! max iterations reached

		       }
		       choose_best(a,b,fa,fb,xzero,fzero);
		  }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void bisection(float(*f)(float x),
		                    const float ax,
		                    const float bx,
				    const float fax,
				    const float fbx,
				    float & xzero,
				    float & fzero,
				    int32_t & iflag) {

			  float x1,x2,x3,f1,f2,f3;
			  bool  root_found = false;
			  //! initialize:
                          iflag = 0;
                          x1    = ax;
                          x2    = bx;
                          f1    = fax;
                          f2    = fbx;
			  for(int32_t i = 1; i != MAXITER; ++i) {
                                 //! bisection of the inclusion interval:
                                 //!  x1------x3------x2
                              x3 = bisect(x1,x2);

                               //! calculate the new function value:
                              f3 = f(x3);
                              if (solution(x3,f3,FTOL4,xzero,fzero)) return;

                              //! determine new inclusion interval:
                              if (f2*f3<0.0f){
                                  //! root lies between x2 and x3
                                  x1 = x3;
                                  x2 = x2;
                                  f1 = f3;
                                  f2 = f2;
			      }
                              else {
                                 //! root lies between x1 and x3
                                  x2 = x3;
                                  f2 = f3;
                             }
                             //! check for convergence:
                             root_found = converged(x1,x2);
                             if(root_found || i==MAXITER) {
                                  choose_best(x1,x2,f1,f2,xzero,fzero);
                                  if (!root_found) iflag = -2;  //! max iterations reached
                                  break;
                              }
			  }
		    }
/*
!  Compute the zero of the function f(x) in the interval ax,bx using the bisection method.
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.1, p 32-34.
*/


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void bisection(double(*f)(double x),
		                    const double ax,
		                    const double bx,
				    const double fax,
				    const double fbx,
				    double & xzero,
				    double & fzero,
				    int32_t & iflag) {

			  double x1,x2,x3,f1,f2,f3;
			  bool  root_found = false;
			  //! initialize:
                          iflag = 0;
                          x1    = ax;
                          x2    = bx;
                          f1    = fax;
                          f2    = fbx;
			  for(int32_t i = 1; i != MAXITER; ++i) {
                                 //! bisection of the inclusion interval:
                                 //!  x1------x3------x2
                              x3 = bisect(x1,x2);

                               //! calculate the new function value:
                              f3 = f(x3);
                              if (solution(x3,f3,FTOL8,xzero,fzero)) return;

                              //! determine new inclusion interval:
                              if (f2*f3<0.0){
                                  //! root lies between x2 and x3
                                  x1 = x3;
                                  x2 = x2;
                                  f1 = f3;
                                  f2 = f2;
			      }
                              else {
                                 //! root lies between x1 and x3
                                  x2 = x3;
                                  f2 = f3;
                             }
                             //! check for convergence:
                             root_found = converged(x1,x2);
                             if(root_found || i==MAXITER) {
                                  choose_best(x1,x2,f1,f2,xzero,fzero);
                                  if (!root_found) iflag = -2;  //! max iterations reached
                                  break;
                              }
			  }
		    }

/*
   !*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the regula falsi method.

*/
		    
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void regula_falsi(float(*f)(float x),
		                       const float ax,
		                       const float bx,
				       const float fax,
				       const float fbx,
				       float & xzero,
				       float & fzero,
				       int32_t & iflag) {

                         float  x1,x2,x3,f1,f2,f3;
			 bool   root_found;
			 //! initialize:
                         iflag = 0;
                         x1    = ax;
                         x2    = bx;
                         f1    = fax;
                         f2    = fbx;
			 for(int32_t i = 1; i != MAXITER; ++i) {
                              //! calculate the new function value:
                              f3 = f(x3);
                              if(solution(x3,f3,FTOL4,xzero,fzero)) return;
                              //! determine new inclusion interval:
				if (f2*f3<0.f) {
                                    //! root lies between x2 and x3
                                    x1 = x3;
                                    x2 = x2;
                                    f1 = f3;
                                    f2 = f2;
				}
                                 else {
                                   //! root lies between x1 and x3
                                   x2 = x3;
                                   f2 = f3;
                                }
                                //! check for convergence:
                               root_found = converged(x1,x2);
                               if (root_found || i==MAXITER) {
                                   choose_best(x1,x2,f1,f2,xzero,fzero);
                                   if(!root_found) iflag = -2;  //! max iterations reached
                                   break;
                               }
			 }
		    }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void regula_falsi(double(*f)(double x),
		                       const double ax,
		                       const double bx,
				       const double fax,
				       const double fbx,
				       double & xzero,
				       double & fzero,
				       int32_t & iflag) {

                         double  x1,x2,x3,f1,f2,f3;
			 bool   root_found;
			 //! initialize:
                         iflag = 0;
                         x1    = ax;
                         x2    = bx;
                         f1    = fax;
                         f2    = fbx;
			 for(int32_t i = 1; i != MAXITER; ++i) {
                              //! calculate the new function value:
                              f3 = f(x3);
                              if(solution(x3,f3,FTOL8,xzero,fzero)) return;
                              //! determine new inclusion interval:
				if (f2*f3<0.0) {
                                    //! root lies between x2 and x3
                                    x1 = x3;
                                    x2 = x2;
                                    f1 = f3;
                                    f2 = f2;
				}
                                 else {
                                   //! root lies between x1 and x3
                                   x2 = x3;
                                   f2 = f3;
                                }
                                //! check for convergence:
                               root_found = converged(x1,x2);
                               if (root_found || i==MAXITER) {
                                   choose_best(x1,x2,f1,f2,xzero,fzero);
                                   if(!root_found) iflag = -2;  //! max iterations reached
                                   break;
                               }
			 }
		    }
		    

		    
   
/*
     !  Illinois method.
!
!### Reference
!  * M. Dowell, P. Jarratt, "A modified regula falsi method for computing the root
!    of an equation', BIT 11 (1971), 168-174.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void illinois(float(*f)(float x),
		                       const float ax,
		                       const float bx,
				       const float fax,
				       const float fbx,
				       float & xzero,
				       float & fzero,
				       int32_t & iflag) {

                         float  x1,x2,x3,f1,f2,f3,delta,f1tmp;
			 bool   root_found;
			 //! initialize:
                         iflag = 0;
                         x1    = ax;
                         x2    = bx;
                         f1    = fax;
                         f2    = fbx;
			 for(int32 i = 1; i != MAXITER; ++i) {
                             x3 = regula_falsi_step(x1,x2,f1,f2,ax,bx);
                             //! calculate the new function value:
                             f3 = f(x3);
                             if(solution(x3,f3,FTOL4,xzero,fzero)) return;
                             //! determine new inclusion interval:
                             if(f2*f3<0.0f) {
                                 //! root lies between x2 and x3
                                 x1 = x2;
                                 x2 = x3;
                                 f1 = f2;
                                 f1tmp = f1;
                                 f2 = f3;
			     }
                             else {
                                 //! root lies between x1 and x3
                                 x2 = x3;
                                 f2 = f3;
                                 f1tmp = f1; //! actual function eval
                                 f1 = 0.5f * f1;
                             }
                             //! check for convergence:
                             root_found = converged(x1,x2);
                             if(root_found || i==MAXITER) {
                                  choose_best(x1,x2,f1tmp,f2,xzero,fzero);
                                  if (!root_found) iflag = -2;  //! max iterations reached
                                  break;
                             }
			 }
		     }

#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void illinois(double(*f)(double x),
		                       const double ax,
		                       const double bx,
				       const double fax,
				       const double fbx,
				       double & xzero,
				       double & fzero,
				       int32_t & iflag) {

                         double  x1,x2,x3,f1,f2,f3,delta,f1tmp;
			 bool   root_found;
			 //! initialize:
                         iflag = 0;
                         x1    = ax;
                         x2    = bx;
                         f1    = fax;
                         f2    = fbx;
			 for(int32 i = 1; i != MAXITER; ++i) {
                             x3 = regula_falsi_step(x1,x2,f1,f2,ax,bx);
                             //! calculate the new function value:
                             f3 = f(x3);
                             if(solution(x3,f3,FTOL8,xzero,fzero)) return;
                             //! determine new inclusion interval:
                             if(f2*f3<0.0) {
                                 //! root lies between x2 and x3
                                 x1 = x2;
                                 x2 = x3;
                                 f1 = f2;
                                 f1tmp = f1;
                                 f2 = f3;
			     }
                             else {
                                 //! root lies between x1 and x3
                                 x2 = x3;
                                 f2 = f3;
                                 f1tmp = f1; //! actual function eval
                                 f1 = 0.5 * f1;
                             }
                             //! check for convergence:
                             root_found = converged(x1,x2);
                             if(root_found || i==MAXITER) {
                                  choose_best(x1,x2,f1tmp,f2,xzero,fzero);
                                  if (!root_found) iflag = -2;  //! max iterations reached
                                  break;
                             }
			 }
		     }

/*
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the Anderson-Bjorck method.
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.2, p 36.
*/		     
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void anderson_bjorck(float(*f)(float x),
		                          const float ax,
		                          const float bx,
				          const float fax,
				          const float fbx,
				          float & xzero,
				          float & fzero,
				          int32_t & iflag) {

                          float  x1,x2,x3,f1,f2,f3,g,f1tmp;
			  bool   root_found = false;
			  //! initialize:
                          iflag = 0;
                          x1    = ax;
                          x2    = bx;
                          f1    = fax;
                          f2    = fbx;
			  for(int32_t i = 1; i != MAXITER; ++i) {
                               x3 = secant(x1,x2,f1,f2,ax,bx);
                               f3 = f(x3);
                               if (solution(x3,f3,FTOL4,xzero,fzero)) return;
                               //! determine a new inclusion interval:
                               if(f2*f3<0.0f) {
                                  //! zero lies between x2 and x3
                                  x1 = x2;
                                  x2 = x3;
                                  f1 = f2;
                                  f2 = f3;
                                  f1tmp = f1;
			       }
                               else {
	    
                                  // zero lies between x1 and x3
                                  g = 1.0f - f3/f2;
                                  if(g<=0.0f) g = 0.5f;
                                  x2 = x3;
                                  f1tmp = f1;
                                  f1 = g*f1;
                                  f2 = f3;
                               }
                               //! check for convergence:
                               root_found = converged(x1,x2);
                               if(root_found || i == MAXITER) {
                                   choose_best(x1,x2,f1tmp,f2,xzero,fzero);
                                   if (!root_found) iflag = -2;  //! max iterations reached
                                   break;
                               }
			  }
                   }

		   
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void anderson_bjorck(double(*f)(double x),
		                          const double ax,
		                          const double bx,
				          const double fax,
				          const double fbx,
				          double & xzero,
				          double & fzero,
				          int32_t & iflag) {

                          double  x1,x2,x3,f1,f2,f3,g,f1tmp;
			  bool   root_found = false;
			  //! initialize:
                          iflag = 0;
                          x1    = ax;
                          x2    = bx;
                          f1    = fax;
                          f2    = fbx;
			  for(int32_t i = 1; i != MAXITER; ++i) {
                               x3 = secant(x1,x2,f1,f2,ax,bx);
                               f3 = f(x3);
                               if (solution(x3,f3,FTOL8,xzero,fzero)) return;
                               //! determine a new inclusion interval:
                               if(f2*f3<0.0) {
                                  //! zero lies between x2 and x3
                                  x1 = x2;
                                  x2 = x3;
                                  f1 = f2;
                                  f2 = f3;
                                  f1tmp = f1;
			       }
                               else {
	    
                                  // zero lies between x1 and x3
                                  g = 1.0 - f3/f2;
                                  if(g<=0.0) g = 0.5;
                                  x2 = x3;
                                  f1tmp = f1;
                                  f1 = g*f1;
                                  f2 = f3;
                               }
                               //! check for convergence:
                               root_found = converged(x1,x2);
                               if(root_found || i == MAXITER) {
                                   choose_best(x1,x2,f1tmp,f2,xzero,fzero);
                                   if (!root_found) iflag = -2;  //! max iterations reached
                                   break;
                               }
			  }
                   }

/*
!  Ridders method to find a root of f(x).
!
!### See also
!  * Ridders, C., "A new algorithm for computing a single root of a real continuous function",
!    IEEE Trans. on Circuits and Systems, Vol 26, Issue 11, Nov 1979.
*/
#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void ridders(float(*f)(float x),
		                  const float ax,
		                  const float bx,
				  const float fax,
				  const float fbx,
				  float & xzero,
				  float & fzero,
				  int32_t & iflag) {

                        float  fh,fl,fm,fnew,denom,xh,xl,xm,xnew;
			//! initialize:
                        iflag = 0;
                        fl    = fax;
                        fh    = fbx;
                        xl    = ax;
                        xh    = bx;
                        xzero = std::numeric_limits<float>::max();
			for(int32_t i = 1; i != MAXITER; ++i) {
                            xm = bisect(xl,xh);
                            fm = f(xm);
                            if(solution(xm,fm,FTOL4,xzero,fzero)) return;
                            denom = std::sqrt(fm*fm-fl*fh);
                            if (denom == 0.0f) {
                                xzero = xm;
                                fzero = fm;
                                iflag = -3;       //! can't proceed: denominator is zero [TODO: add a bisection if this happens]
                                break;
                            }
                            xnew = xm + (xm-xl)*(std::copysign(1.0f,fl-fh)*fm/denom);
                            if(converged(xzero,xnew)) {  //! relative convergence in x
                               //! additional check to prevent false convergence
                                if (converged(xl,xm) || converged(xm,xh)) break;
                            }
                            xzero = xnew;
                            fnew  = f(xzero);
                            fzero = fnew;
                            if(std::abs(fnew) <= FTOL4) break;    //! abs convergence in f
                             //! to keep the root bracketed:
                            if(std::copysign(fm,fnew) != fm) {
                               xl = xm;
                               fl = fm;
                               xh = xzero;
                               fh = fnew;
			    }
                            else if(std::copysign(fl,fnew) != fl) {
                               xh = xzero;
                               fh = fnew;
			    }
                            else if(std::copysign(fh,fnew) /= fh) {
                               xl = xzero
                               fl = fnew
                            }
                            if(converged(xl,xh)) break;    //! relative convergence in x
                            if(i == MAXITER) iflag = -2;  //! max iterations exceeded

			}
		   }


#if defined(__INTEL_COMPILER) || defined(__ICC)
#pragma intel optimization_level 3
#pragma intel optimization_parameter target_arch=AVX
#endif
                     __ATTR_ALWAYS_INLINE
		     __ATTR_HOT__
		     __ATTR_ALIGN__(32)
		     static
		     inline
		     void ridders(double(*f)(double x),
		                  const double ax,
		                  const double bx,
				  const double fax,
				  const double fbx,
				  double & xzero,
				  double & fzero,
				  int32_t & iflag) {

                        double  fh,fl,fm,fnew,denom,xh,xl,xm,xnew;
			//! initialize:
                        iflag = 0;
                        fl    = fax;
                        fh    = fbx;
                        xl    = ax;
                        xh    = bx;
                        xzero = std::numeric_limits<double>::max();
			for(int32_t i = 1; i != MAXITER; ++i) {
                            xm = bisect(xl,xh);
                            fm = f(xm);
                            if(solution(xm,fm,FTOL8,xzero,fzero)) return;
                            denom = std::sqrt(fm*fm-fl*fh);
                            if (denom == 0.0) {
                                xzero = xm;
                                fzero = fm;
                                iflag = -3;       //! can't proceed: denominator is zero [TODO: add a bisection if this happens]
                                break;
                            }
                            xnew = xm + (xm-xl)*(std::copysign(1.0,fl-fh)*fm/denom);
                            if(converged(xzero,xnew)) {  //! relative convergence in x
                               //! additional check to prevent false convergence
                                if (converged(xl,xm) || converged(xm,xh)) break;
                            }
                            xzero = xnew;
                            fnew  = f(xzero);
                            fzero = fnew;
                            if(std::abs(fnew) <= FTOL8) break;    //! abs convergence in f
                             //! to keep the root bracketed:
                            if(std::copysign(fm,fnew) != fm) {
                               xl = xm;
                               fl = fm;
                               xh = xzero;
                               fh = fnew;
			    }
                            else if(std::copysign(fl,fnew) != fl) {
                               xh = xzero;
                               fh = fnew;
			    }
                            else if(std::copysign(fh,fnew) /= fh) {
                               xl = xzero
                               fl = fnew
                            }
                            if(converged(xl,xh)) break;    //! relative convergence in x
                            if(i == MAXITER) iflag = -2;  //! max iterations exceeded

			}
		   }
				       

     } //math

} //gms





#endif /*__GMS_ROOT_FINDING_HPP__*/
