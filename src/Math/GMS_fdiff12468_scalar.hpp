//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Modified by Bernard Gingold
//  License:
//  MIT License

/*Copyright (c) 2020 Bernard Gingold

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
#ifndef __GMS_FDIFF12468_SCALAR_HPP__
#define __GMS_FDIFF12468_SCALAR_HPP__ 110220220938

/*
 * Performs numerical differentiation by finite-differences.
 *
 * All numerical differentiation using finite-differences are ill-conditioned, and these routines are no exception.
 * A simple argument demonstrates that the error is unbounded as h->0.
 * Take the one sides finite difference formula f'(x) = (f(x+h)-f(x))/h.
 * The evaluation of f induces an error as well as the error from the finite-difference approximation, giving
 * |f'(x) - (f(x+h) -f(x))/h| < h|f''(x)|/2 + (|f(x)|+|f(x+h)|)eps/h =: g(h), where eps is the unit roundoff for the type.
 * It is reasonable to choose h in a way that minimizes the maximum error bound g(h).
 * The value of h that minimizes g is h = sqrt(2eps(|f(x)| + |f(x+h)|)/|f''(x)|), and for this value of h the error bound is
 * sqrt(2eps(|f(x+h) +f(x)||f''(x)|)).
 * In fact it is not necessary to compute the ratio (|f(x+h)| + |f(x)|)/|f''(x)|; the error bound of ~\sqrt{\epsilon} still holds if we set it to one.
 *
 *
 * For more details on this method of analysis, see
 *
 * http://www.uio.no/studier/emner/matnat/math/MAT-INF1100/h08/kompendiet/diffint.pdf
 * http://web.archive.org/web/20150420195907/http://www.uio.no/studier/emner/matnat/math/MAT-INF1100/h08/kompendiet/diffint.pdf
 *
 *
 * It can be shown on general grounds that when choosing the optimal h, the maximum error in f'(x) is ~(|f(x)|eps)^k/k+1|f^(k-1)(x)|^1/k+1.
 * From this we can see that full precision can be recovered in the limit k->infinity.
 *
 * References:
 *
 * 1) Fornberg, Bengt. "Generation of finite difference formulas on arbitrarily spaced grids." Mathematics of computation 51.184 (1988): 699-706.
 *
 *
 * The second algorithm, the complex step derivative, is not ill-conditioned.
 * However, it requires that your function can be evaluated at complex arguments.
 * The idea is that f(x+ih) = f(x) +ihf'(x) - h^2f''(x) + ... so f'(x) \approx Im[f(x+ih)]/h.
 * No subtractive cancellation occurs. The error is ~ eps|f'(x)| + eps^2|f'''(x)|/6; hard to beat that.
 *
 * References:
 *
 * 1) Squire, William, and George Trapp. "Using complex variables to estimate derivatives of real functions." Siam Review 40.1 (1998): 110-112.
 */

#include <complex>
#include <cmath>
#include "GMS_config.h"


namespace gms {

   namespace math {

    __ATTR_ALIGN__(32)
    template<class Real>
    Real make_xph_representable(Real x, Real h)
    {
        using std::numeric_limits;
        // Redefine h so that x + h is representable. Not using this trick leads to large error.
        // The compiler flag -ffast-math evaporates these operations . . .
        Real temp = x + h;
        h = temp - x;
        // Handle the case x + h == x:
        if (h == 0)
        {
	   h = std::nextafter(x, (numeric_limits<Real>::max)()) - x;
        }
        return h;
    }

__ATTR_ALIGN__(32)
template<class F, class Real>
Real complex_step_derivative(const F f, Real x)
{
    // Is it really this easy? Yes.
    // Note that some authors recommend taking the stepsize h to be smaller than epsilon(), some recommending use of the min().
    // This idea was tested over a few billion test cases and found the make the error *much* worse.
    // Even 2eps and eps/2 made the error worse, which was surprising.
    using std::complex;
    using std::numeric_limits;
    constexpr const Real step = (numeric_limits<Real>::epsilon)();
    constexpr const Real inv_step = 1/(numeric_limits<Real>::epsilon)();
    return f(complex<Real>(x, step)).imag()*inv_step;
}



  
   __ATTR_ALIGN__(32)
   template<class F, class Real>
   Real finite_diff_1(const F f, Real x, Real* error)
   {
      using std::sqrt;
      using std::pow;
      using std::abs;
      using std::numeric_limits;

      const Real eps = (numeric_limits<Real>::epsilon)();
      // Error bound ~eps^1/2
      // Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
      // Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
      // This approximation will get better as we move to higher orders of accuracy.
      Real h = (Real)2.0 * sqrt(eps);
      h = make_xph_representable(x, h);

      Real yh = f(x + h);
      Real y0 = f(x);
      Real diff = yh - y0;
      if (error)
      {
         Real ym = f(x - h);
         Real ypph = abs(yh - (Real)2.0 * y0 + ym) / h;
         // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
         *error = ypph / (Real)2.0 + (abs(yh) + abs(y0))*eps / h;
      }
      return diff / h;
   }

    __ATTR_ALIGN__(32)
   template<class F, class Real>
   Real finite_diff_2(const F f, Real x, Real* error)
   {
      using std::sqrt;
      using std::pow;
      using std::abs;
      using std::numeric_limits;

      const Real eps = (numeric_limits<Real>::epsilon)();
      // Error bound ~eps^2/3
      // See the previous discussion to understand determination of h and the error bound.
      // Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]
      Real h = pow((Real)3.0 * eps, static_cast<Real>(1.0) / static_cast<Real>(3.0));
      h = make_xph_representable(x, h);

      Real yh = f(x + h);
      Real ymh = f(x - h);
      Real diff = yh - ymh;
      if (error)
      {
         Real yth = f(x + (Real)2.0 * h);
         Real ymth = f(x - (Real)2.0 * h);
         *error = eps * (abs(yh) + abs(ymh)) / ((Real)2.0 * h) + abs((yth - ymth) / (Real)2.0 - diff) / ((Real)6.0 * h);
      }

      return diff / ((Real)2.0 * h);
   }

    __ATTR_ALIGN__(32)
   template<class F, class Real>
   Real finite_diff_4(const F f, Real x, Real* error)
   {
      using std::sqrt;
      using std::pow;
      using std::abs;
      using std::numeric_limits;

      const Real eps = (numeric_limits<Real>::epsilon)();
      // Error bound ~eps^4/5
      Real h = pow(11.25*eps, (Real)1 / (Real)5);
      h = make_xph_representable(x, h);
      Real ymth = f(x - 2 * h);
      Real yth = f(x + 2 * h);
      Real yh = f(x + h);
      Real ymh = f(x - h);
      Real y2 = ymth - yth;
      Real y1 = yh - ymh;
      if (error)
      {
         // Mathematica code to extract the remainder:
         // Series[(f[x-2*h]+ 8*f[x+h] - 8*f[x-h] - f[x+2*h])/(12*h), {h, 0, 7}]
         Real y_three_h = f(x + 3 * h);
         Real y_m_three_h = f(x - 3 * h);
         // Error from fifth derivative:
         *error = abs((y_three_h - y_m_three_h) / 2 + 2 * (ymth - yth) + 5 * (yh - ymh) / 2) / (30 * h);
         // Error from function evaluation:
         *error += eps * (abs(yth) + abs(ymth) + 8 * (abs(ymh) + abs(yh))) / (12 * h);
      }
      return (y2 + 8 * y1) / (12 * h);
   }

    __ATTR_ALIGN__(32)
   template<class F, class Real>
   Real finite_diff_6(const F f, Real x, Real* error)
   {
      using std::sqrt;
      using std::pow;
      using std::abs;
      using std::numeric_limits;

      const Real eps = (numeric_limits<Real>::epsilon)();
      // Error bound ~eps^6/7
      // Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
      Real h = pow(eps / 168, (Real)1 / (Real)7);
      h = make_xph_representable(x, h);

      Real yh = f(x + h);
      Real ymh = f(x - h);
      Real y1 = yh - ymh;
      Real y2 = f(x - 2 * h) - f(x + 2 * h);
      Real y3 = f(x + 3 * h) - f(x - 3 * h);

      if (error)
      {
         // Mathematica code to generate fd scheme for 7th derivative:
         // Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
         // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
         // Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
         Real y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
         *error = abs(y7) / (140 * h) + 5 * (abs(yh) + abs(ymh))*eps / h;
      }
      return (y3 + 9 * y2 + 45 * y1) / (60 * h);
   }

    __ATTR_ALIGN__(32)
   template<class F, class Real>
   Real finite_diff_8(const F f, Real x, Real* error)
   {
      using std::sqrt;
      using std::pow;
      using std::abs;
      using std::numeric_limits;

      const Real eps = (numeric_limits<Real>::epsilon)();
      // Error bound ~eps^8/9.
      // In double precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
      // Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
      // Mathematica code to get the error:
      // Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
      // If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.
      Real h = pow(551.25*eps, (Real)1 / (Real)9);
      h = detail::make_xph_representable(x, h);

      Real yh = f(x + h);
      Real ymh = f(x - h);
      Real y1 = yh - ymh;
      Real y2 = f(x - 2 * h) - f(x + 2 * h);
      Real y3 = f(x + 3 * h) - f(x - 3 * h);
      Real y4 = f(x - 4 * h) - f(x + 4 * h);

      Real tmp1 = 3 * y4 / 8 + 4 * y3;
      Real tmp2 = 21 * y2 + 84 * y1;

      if (error)
      {
         // Mathematica code to generate fd scheme for 7th derivative:
         // Sum[(-1)^i*Binomial[9, i]*(f[x+(4-i)*h] + f[x+(5-i)*h])/2, {i, 0, 9}]
         // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
         // Series[(f[x+5*h]-f[x- 5*h])/2 + 4*(f[x-4*h] - f[x+4*h]) + 27*(f[x+3*h] - f[x-3*h])/2 + 24*(f[x-2*h]  - f[x+2*h]) + 21*(f[x+h] - f[x-h]), {h, 0, 15}]
         Real f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
         *error = abs(f9) / (630 * h) + 7 * (abs(yh) + abs(ymh))*eps / h;
      }
      return (tmp1 + tmp2) / (105 * h);
   }


  } // math

} //gms  


#endif
